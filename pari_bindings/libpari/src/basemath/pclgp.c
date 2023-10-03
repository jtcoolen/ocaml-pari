/* Copyright (C) 2000, 2012  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_subcyclo

/* written by Takashi Fukuda */

#define onevec(x) const_vec(x,gen_1)
#define nullvec() cgetg(1, t_VEC)
#define order_f_x(f, x) Fl_order(x%f, eulerphiu(f), f)

#define USE_MLL         (1L<<0)
#define NO_PLUS_PART    (1L<<1)
#define NO_MINUS_PART   (1L<<2)
#define SKIP_PROPER     (1L<<3)
#define SAVE_MEMORY     (1L<<4)
#define USE_FULL_EL     (1L<<5)
#define USE_BASIS       (1L<<6)
#define USE_FACTOR      (1L<<7)
#define USE_GALOIS_POL  (1L<<8)
#define USE_F           (1L<<9)

static ulong
_get_d(GEN H) { return umael(H, 2, 1);}
static ulong
_get_f(GEN H) { return umael(H, 2, 2);}
static ulong
_get_h(GEN H) { return umael(H, 2, 3);}
static long
_get_s(GEN H) { return umael(H, 2, 4);}
static long
_get_g(GEN H) { return umael(H, 2, 5);}
static GEN
_get_H(GEN H) { return gel(H, 3);}
static ulong
K_get_d(GEN K) { return _get_d(gel(K,1)); }
static ulong
K_get_f(GEN K) { return _get_f(gel(K,1)); }
static ulong
K_get_h(GEN K) { return _get_h(gel(K,1)); }
static long
K_get_s(GEN K) { return _get_s(gel(K,1)); }
static ulong
K_get_g(GEN K) { return _get_g(gel(K,1)); }
static GEN
K_get_H(GEN K) { return _get_H(gel(K,1)); }
static ulong
K_get_dchi(GEN K) { return gel(K,6)[1]; }
static ulong
K_get_nconj(GEN K) { return gel(K,6)[2]; }

/* G=<s> is a cyclic group of order n and t=s^(-1).
 *  convert sum_i a_i*s^i to sum_i b_i*t^i */
static GEN
Flx_recip1_inplace(GEN x, long pn)
{
  long i, lx = lg(x);
  if(lx-2 != pn) /* This case scarcely occurs */
  {
    long ly = pn+2;
    GEN y = const_vecsmall(ly, 0);
    y[1] = x[1];y[2] = x[2];
    for(i=3;i<lx;i++) y[ly+2-i] = x[i];
    return Flx_renormalize(y, ly);
  }
  else /* almost all cases */
  {
    long t, mid = (lx+1)>>1;
    for(i=3;i<=mid;i++)
    {
      t = x[i];x[i] = x[lx+2-i];x[lx+2-i] = t;
    }
    return Flx_renormalize(x, lx);
  }
}

/* Return h^degpol(P) P(x / h) */
static GEN
Flx_rescale_inplace(GEN P, ulong h, ulong p)
{
  long i, l = lg(P);
  ulong hi = h;
  for (i=l-2; i>=2; i--)
  {
    P[i] = Fl_mul(P[i], hi, p);
    if (i == 2) break;
    hi = Fl_mul(hi,h, p);
  }
  return P;
}

static GEN
zx_to_Flx_inplace(GEN x, ulong p)
{
  long i, lx = lg(x);
  for (i=2; i<lx; i++) uel(x,i) = umodsu(x[i], p);
  return Flx_renormalize(x, lx);
}

/* zero pol of n components (i.e. deg=n-1). need to pass to ZX_renormalize */
INLINE GEN
pol_zero(long n)
{
  long i;
  GEN p = cgetg(n+2, t_POL);
  p[1] = evalsigne(1) | evalvarn(0);
  for (i = 2; i < n+2; i++) gel(p, i) = gen_0;
  return p;
}

/* e[i+1] = L*i + K for i >= n; determine K,L and reduce n if possible */
static GEN
vecsmall2vec2(GEN e, long n)
{
  long L = e[n+1] - e[n], K = e[n+1] - L*n;
  n--; while (n >= 0 && e[n+1] - L*n == K) n--;
  if (n < 0) e = nullvec(); else { setlg(e, n+2); e = zv_to_ZV(e); }
  return mkvec3(utoi(L), stoi(K), e); /* L >= 0 */
}

/* z=zeta_{p^n}; return k s.t. (z-1)^k || f(z) assuming deg(f)<phi(p^n) */
static long
zx_p_val(GEN f, ulong p, ulong n)
{
  pari_sp av = avma;
  ulong x = zx_lval(f, p);
  if (x) { f = zx_z_divexact(f, upowuu(p, x)); x *= (p-1)*upowuu(p, n-1); }
  x += Flx_val(Flx_translate1(zx_to_Flx(f, p), p));
  return gc_long(av, x);
}

static long
ZX_p_val(GEN f, ulong p, ulong n)
{
  pari_sp av = avma;
  ulong x = ZX_lval(f, p);
  if (x) { f = ZX_Z_divexact(f, powuu(p, x)); x *= (p-1)*upowuu(p, n-1); }
  x += Flx_val(Flx_translate1(ZX_to_Flx(f, p), p));
  return gc_long(av, x);
}

static GEN
set_A(GEN B, int *chi)
{
  long a, i, j, B1 = B[1], l = lg(B);
  GEN A = cgetg(l, t_VECSMALL);
  for (a = 0, j = 1; j < B1; j++) a += chi[j];
  A[1] = a;
  for (i = 2; i < l; i++)
  {
    long Bi = B[i];
    for (a = A[i-1], j = B[i-1]; j < Bi; j++) a += chi[j];
    A[i] = a;
  }
  return A;
}

/* g_n(a)=g_n(b) <==> a^2=b^2 mod 2^(n+2) <==> a=b,-b mod 2^(n+2)
 * g_n(a)=g_n(1+q0)^k <==> a=x(1+q0)^k x=1,-1
 * gam[1+a]=k, k<0  ==> g_n(a)=0
 *             k>=0 ==> g_n(a)^(-1)=gamma^k, gamma=g_n(1+q0) */
static GEN
set_gam2(long q01, long n)
{
  long i, x, x1, pn, pn2;
  GEN gam;
  pn = (1L<<n);
  pn2 = (pn<<2);
  gam = const_vecsmall(pn2, -1);
  x=Fl_inv(q01, pn2); x1=1;
  for (i=0; i<pn; i++)
  {
    gam[1+x1] = gam[1+Fl_neg(x1, pn2)] = i;
    x1 = Fl_mul(x1, x, pn2);
  }
  return gam;
}

/* g_n(a)=g_n(b) <==> a^(p-1)=b^(p-1) mod p^(n+1) <==> a=xb x=<g^(p^n)>
 * g_n(a)=g_n(1+q0)^k <==> a=x(1+q0)^k x=<g^(p^n)>
 * gam[1+a]=k, k<0  ==> g_n(a)=0
 *             k>=0 ==> g_n(a)^(-1)=gamma^k, gamma=g_n(1+q0) */
static GEN
set_gam(long q01, long p, long n)
{
  long i, j, g, g1, x, x1, p1, pn, pn1;
  GEN A, gam;
  p1 = p-1; pn = upowuu(p, n); pn1 = p*pn;
  gam = const_vecsmall(pn1, -1);
  g = pgener_Zl(p); g1 = Fl_powu(g, pn, pn1);
  A = Fl_powers(g1, p1-1, pn1);  /* A[1+i]=g^(i*p^n) mod p^(n+1), 0<=i<=p-2 */
  x = Fl_inv(q01, pn1); x1 = 1;
  for (i=0; i<pn; i++)
  {
    for (j=1; j<=p1; j++) gam[1+Fl_mul(x1, A[j], pn1)] = i;
    x1 = Fl_mul(x1, x, pn1);
  }
  return gam;
}

/* k=Q(sqrt(m)), A_n=p-class gr. of k_n, |A_n|=p^(e_n)
 * return e_n-e_(n-1)
 * essential assumption : m is not divisible by p
 * Gold, Acta Arith. XXVI (1974), p.22 formula (3) */
static long
ediff(ulong p, long m, ulong n, int *chi)
{
  pari_sp av = avma;
  long j, lx, *px;
  ulong i, d, s, y, g, p1, pn, pn1, pn_1, phipn, phipn1;
  GEN A, B, x, gs, cs;

  d=((m-1)%4==0)?labs(m):4*labs(m);
  p1=p-1; pn_1=upowuu(p, n-1); pn=p*pn_1; pn1=p*pn; phipn=p1*pn_1; phipn1=p1*pn;
  lx=2*p1*phipn;
  y=Fl_inv(pn1%d, d); g=pgener_Zl(p);  /* pn1 may > d */
  cs = cgetg(2+phipn, t_VECSMALL); cs[1] = evalvarn(0);
  x = cgetg(1+lx, t_VECSMALL);
  gs = Fl_powers(g, phipn1-1, pn1); /* gs[1+i]=g^i(mod p^(n+1)), 0<=i<p^(n+1) */

  for (px=x,i=0; i<p1; i++)
  {
    long ipn=i*pn+1,ipnpn=ipn+phipn;
    for (s=0; s<phipn; s++)
    {
      *++px = (y*gs[s+ipn])%d;  /* gs[s+ipn] may > d */
      *++px = (y*gs[(s%pn_1)+ipnpn])%d;
    }
  }
  B = vecsmall_uniq(x);
  A = set_A(B, chi);
  for (s=0; s<phipn; s++)
  {
    long a=0, ipn=1, spn1=s%pn_1;
    for (i=0; i<p1; i++)
    {
      if ((j=zv_search(B, (y*gs[s+ipn])%d))<=0)
        pari_err_BUG("zv_search failed\n");
      a+=A[j];
      if ((j=zv_search(B, (y*gs[spn1+ipn+phipn])%d))<=0)
        pari_err_BUG("zv_search failed\n");
      a-=A[j];
      ipn+=pn;
    }
    cs[2+s] = a;
  }
  cs = zx_renormalize(cs, lg(cs));
  y = (lg(cs)==3) ? phipn*z_lval(cs[2], p) : zx_p_val(cs, p, n);
  return gc_long(av, y);
}

static GEN
quadteichstk(GEN Chi, int *chi, GEN Gam, long p, long m, long n)
{
  GEN Gam1 = Gam+1, xi;
  long i, j, j0, d, f0, pn, pn1, deg, pn1d;

  d = ((m&3)==1)?m:m<<2;
  f0 = ulcm(p, d)/p;
  pn = upowuu(p, n); pn1 = p*pn; pn1d = pn1%d;
  xi = cgetg(pn+2, t_POL); xi[1] = evalsigne(1) | evalvarn(0);
  for (i=0; i<pn; i++) gel(xi, 2+i) = const_vecsmall(p, 0);
  for (j=1; j<pn1; j++)
  {
    long jp, ipn1d, *xij0;
    if ((j0 = Gam1[j])<0) continue;
    jp = j%p; ipn1d = j%d; xij0 = gel(xi, 2+j0)+2;
    for (i=1; i<f0; i++)
    {
      int sgn;
      if ((ipn1d += pn1d) >= d) ipn1d -= d;
      if ((sgn = chi[ipn1d])==0) continue;
      deg = Chi[jp];  /* jp!=0 because j0>=0 */
      if (sgn>0) xij0[deg] += i;
      else xij0[deg] -= i;
    }
  }
  for (i=0; i<pn; i++) gel(xi, 2+i) = zx_renormalize(gel(xi, 2+i), p+1);
  return FlxX_renormalize(xi, pn+2);  /* zxX_renormalize does not exist */
}

#ifdef DEBUG_QUADSTK
/* return f0*xi_n */
static GEN
quadstkp_by_def(int *chi, GEN gam, long n, long p, long f, long f0)
{
  long i, a, a1, pn, pn1, qn;
  GEN x, x2, gam1 = gam+1;
  pn = upowuu(p, n); pn1 = p*pn; qn = f0*pn1;
  x = const_vecsmall(pn+1, 0); x2 = x+2;
  for (a=1; a<qn; a++)
  {
    int sgn;
    if ((a1=gam1[a%pn1])<0 || (sgn=chi[a%f])==0) continue;
    if (sgn>0) x2[a1]+=a;
    else x2[a1]-=a;
  }
  for (i=0; i<pn; i++)
  {
    if (x2[i]%pn1) pari_err_BUG("stickel. ele. is not integral.\n");
    else x2[i]/=pn1;
  }
  return zx_renormalize(x, pn+2);
}
#endif

/* f!=p
 * xi_n = f0^(-1)*
 *   sum_{0<=j<pn1,(j,p)=1}(Q_n/Q,j)^(-1)*(sum_{0<=i<f0}i*chi^(-1)(pn1*i+j)) */
static GEN
quadstkp1(int *chi, GEN gam, long n, long p, long f, long f0)
{
  long i, j, j0, pn, pn1, pn1f, den;
  GEN x, x2;
  pn = upowuu(p, n); pn1 = p*pn; pn1f = pn1%f;
  x = const_vecsmall(pn+1, 0); x2 = x+2;
  if (f==3) den = (chi[p%f]>0)?f0<<1:2;
  else if (f==4) den = (chi[p%f]>0)?f0<<1:f0;
  else den = f0<<1;
  for (j=1; j<pn1; j++)
  {
    long ipn1;
    if (j%p==0) continue;
    j0 = gam[1+j]; ipn1 = j%f;
    for (i=1; i<f0; i++)
    {
      int sgn;
      if ((ipn1+=pn1f)>=f) ipn1-=f;
      if ((sgn = chi[ipn1])>0) x2[j0]+=i;
      else if (sgn<0) x2[j0]-=i;
    }
  }
  for (i=0; i<pn; i++)
  {
    if (x2[i]%den) pari_err_BUG("stickel. ele. is not integral.\n");
    else x2[i]/=den;
  }
  return zx_renormalize(x, pn+2);
}

/* f==p */
static GEN
quadstkp2(int *chi, GEN gam, long n, long p)
{
  long a, a1, i, pn, pn1, amodp;
  GEN x, x2, gam1 = gam+1;
  pn = upowuu(p, n); pn1 = p*pn;
  x = const_vecsmall(pn+1, 0); x2 = x+2;
  for (a=1,amodp=0; a<pn1; a++)
  {
    int sgn;
    if (++amodp==p) {amodp = 0; continue; }
    if ((sgn = chi[amodp])==0) continue;
    a1=gam1[a];
    if (sgn>0) x2[a1]+=a;
    else x2[a1]-=a;
  }
  for (i=0; i<pn; i++)
  {
    if (x2[i]%pn1) pari_err_BUG("stickel. ele. is not integral.\n");
    else x2[i]/=pn1;
  }
  return zx_renormalize(x, pn+2);
}

/*  p>=3
 *  f = conductor of Q(sqrt(m))
 *  q0 = lcm(f,p) = f0*p
 *  qn = q0*p^n = f0*p^(n+1)
 *  xi_n = qn^(-1)*sum_{1<=a<=qn,(a,qn)=1} a*chi(a)^(-1)*(Q_n/Q,a)^(-1) */
static GEN
quadstkp(long p, long m, long n, int *chi)
{
  long f, f0, pn, pn1, q0;
  GEN gam;
  f = ((m-1)%4==0)?labs(m):4*labs(m);
  pn = upowuu(p, n); pn1 = p*pn;
  if (f % p) { q0 = f * p; f0 = f; } else { q0 = f; f0 = f / p; }
  gam = set_gam((1+q0)%pn1, p, n);
#ifdef DEBUG_QUADSTK
  return quadstkp_by_def(chi, gam, n, p, f, f0);
#else
  return (f0!=1)?quadstkp1(chi, gam, n, p, f, f0):quadstkp2(chi, gam, n, p);
#endif
}

/* p=2 */
static GEN
quadstk2(long m, long n, int *chi)
{
  long i, j, j0, f, f0, pn, pn1, pn2, pn2f, q0;
  GEN x, x2, gam;
  f = ((m-1)%4==0)?labs(m):4*labs(m);
  pn = 1L<<n; pn1 = pn<<1; pn2 = pn1<<1; pn2f = pn2%f;
  q0 = (f&1)?f*4:f;
  f0 = (f&1)?f:f/4;
  x = const_vecsmall(pn+1, 0); x2 = x+2;
  gam = set_gam2((1+q0)%pn2, n);
  for (j=1; j<pn2; j++)
  {
    long ipn2;
    if (!(j&1)) continue;
    j0 = gam[1+j];
    ipn2 = j%f;
    /* for (i=1; i<f0; i++) x2[j0]+=i*chi[(i*pn2+j)%f]; */
    for (i=1; i<f0; i++)
    {
      int sgn;
      if ((ipn2+=pn2f)>=f) ipn2-=f;
      if ((sgn=chi[ipn2])>0) x2[j0]+=i;
      else if (sgn<0) x2[j0]-=i;
    }
  }
  for (f0<<=1, i=0; i<pn; i++)
  {
    if (x2[i]%f0) pari_err_BUG("stickel. ele. is not integral.\n");
    else x2[i]/=f0;
  }
  return zx_renormalize(x, pn+2);
}

/* Chin is a generator of the group of the characters of G(Q_n/Q).
 * chin[1+a]=k, k<0  ==> Chin(a)=0
 *              k>=0 ==> Chin(a)=zeta_{p^n}^k */
static GEN
set_chin(long p, long n)
{
  long i, j, x = 1, g, gpn, pn, pn1;
  GEN chin, chin1;
  pn = upowuu(p, n); pn1 = p*pn;
  chin = const_vecsmall(pn1, -1); chin1 = chin+1;
  g = pgener_Zl(p); gpn = Fl_powu(g, pn, pn1);
  for (i=0; i<pn; i++)
  {
    long y = x;
    for (j=1; j<p; j++)
    {
      chin1[y] = i;
      y = Fl_mul(y, gpn, pn1);
    }
    x = Fl_mul(x, g, pn1);
  }
  return chin;
}

/* k=Q(sqrt(m)), A_n=p-class gr. of k_n, |A_n|=p^(e_n), p|m
 * return e_n-e_(n-1).
 * There is an another method using the Stickelberger element based on
 * Coates-Lichtenbaum, Ann. Math. vol.98 No.3 (1973), 498-550, Lemma 2.15.
 * If kro(m,p)!=1, then orders of two groups coincide.
 * ediff_ber is faster than the Stickelberger element. */
static long
ediff_ber(ulong p, long m, ulong n, int *chi)
{
  pari_sp av = avma;
  long a, d, e, x, y, pn, pn1, qn1;
  GEN B, B2, chin = set_chin(p, n)+1;

  d = ((m-1)%4==0)?labs(m):4*labs(m);
  pn = upowuu(p, n); pn1 = p*pn; qn1 = (d*pn)>>1;
  B = const_vecsmall(pn+1, 0); B2 = B+2;
  for (a=x=y=1; a <= qn1; a++) /* x=a%d, y=a%pn1 */
  {
    int sgn = chi[x];
    if (sgn)
    {
      long k = chin[y];
      if (k >= 0) { if (sgn > 0) B2[k]++; else B2[k]--; }
    }
    if (++x == d) x = 0;
    if (++y == pn1) y = 0;
  }
  B = zx_renormalize(B, pn+2);
  e = (n==1)? zx_p_val(B, p, n)
            : ZX_p_val(ZX_rem(zx_to_ZX(B), polcyclo(pn, 0)), p, n);
  if (p==3 && chi[2] < 0) e--;  /* 2 is a primitive root of 3^n (n>=1) */
  return gc_long(av, e);
}

#ifdef DEBUG
/* slow */
static int*
set_quad_chi_1(long m)
{
  long a, d, f;
  int *chi;
  d=((m-1)%4==0)?m:4*m; f=labs(d);
  chi= (int*)stack_calloc(sizeof(int)*f);
  for (a=1; a<f; a++) chi[a]=kross(d, a);
  return chi;
}
#endif

/* chi[a]=kross(d, a)   0<=a<=f-1
 * d=discriminant of Q(sqrt(m)), f=abs(d)
 *
 * Algorithm: m=-p1*p2*...*pr ==> kross(d,gi)=-1 (1<=i<=r), gi=proot(pi)
 * set_quad_chi_1(m)=set_quad_chi_2(m) for all square-free m s.t. |m|<10^5. */
static int*
set_quad_chi_2(long m)
{
  long d = (m-1) % 4? 4*m: m, f = labs(d);
  GEN fa = factoru(f), P = gel(fa, 1), E = gel(fa,2), u, v;
  long i, j, np, nm, l = lg(P);
  int *chi = (int*)stack_calloc(sizeof(int)*f);
  pari_sp av = avma;
  int *plus = (int*)stack_calloc(sizeof(int)*f), *p0 = plus;
  int *minus = (int*)stack_calloc(sizeof(int)*f), *p1 = minus;

  u = cgetg(32, t_VECSMALL);
  v = cgetg(32, t_VECSMALL);
  for (i = 1; i < l; i++)
  {
    ulong p = upowuu(P[i], E[i]);
    u[i] = p * Fl_inv(p, f / p);
    v[i] = Fl_sub(1, u[i], f);
  }
  if (E[1]==2)       /* f=4*(-m) */
  {
    *p0++ = Fl_add(v[1], u[1], f);
    *p1++ = Fl_add(Fl_mul(3, v[1], f), u[1], f);
    i = 2;
  }
  else if (E[1]==3)  /* f=8*(-m) */
  {
    ulong a;
    *p0++ = Fl_add(v[1], u[1], f);
    a = Fl_add(Fl_mul(3, v[1], f), u[1], f);
    if (kross(d, a) > 0) *p0++ = a; else *p1++ = a;
    a = Fl_add(Fl_mul(5, v[1], f), u[1], f);
    if (kross(d, a) > 0) *p0++ = a; else *p1++ = a;
    a = Fl_add(Fl_mul(7, v[1], f), u[1], f);
    if (kross(d, a) > 0) *p0++ = a; else *p1++ = a;
    i = 2;
  }
  else              /* f=-m */
  {*p0++ = 1; i = 1; }
  for (; i < l; i++)
  {
    ulong gn, g = pgener_Fl(P[i]);
    gn = g = Fl_add(Fl_mul(g, v[i], f), u[i], f);
    np = p0-plus;
    nm = p1-minus;
    for (;;)
    {
      for (j = 0; j < np; j++) *p1++ = Fl_mul(plus[j], gn, f);
      for (j = 0; j < nm; j++) *p0++ = Fl_mul(minus[j], gn, f);
      gn = Fl_mul(gn, g, f); if (gn == 1) break;
      for (j= 0; j < np; j++) *p0++ = Fl_mul(plus[j], gn, f);
      for (j = 0; j < nm; j++) *p1++ = Fl_mul(minus[j], gn, f);
      gn = Fl_mul(gn, g, f); if (gn == 1) break;
    }
  }
  np = p0-plus;
  nm = p1-minus;
  for (i = 0; i < np; i++) chi[plus[i]] = 1;
  for (i = 0; i < nm; i++) chi[minus[i]] = -1;
  set_avma(av); return chi;
}

static long
srh_x(GEN T, long n, long x)
{
  for (; x<n; x++) if (!T[x]) return x;
  return -1;
}

/* G is a cyclic group of order d. hat(G)=<chi>.
 * chi, chi^p, ... , chi^(p^(d_chi-1)) are conjugate.
 * {chi^j | j in C} are repre. of Q_p-congacy classes of inj. chars.
 *
 * C is a set of representatives of H/<p>, where H=(Z/dZ)^* */
static GEN
set_C(long p, long d, long d_chi, long n_conj)
{
  long i, j, x, y, pmodd = p%d;
  GEN T = const_vecsmall(d, 0)+1;
  GEN C = cgetg(1+n_conj, t_VECSMALL);
  if (n_conj==1) { C[1] = 1; return C; }
  for (i=0, x=1; x >= 0; x = srh_x(T, d, x))
  {
    if (cgcd(x, d)==1) C[++i] = x;
    for (j=0, y=x; j<d_chi; j++) T[y = Fl_mul(y, pmodd, d)] = 1;
  }
  return C;
}

static GEN
FpX_one_cyclo(long n, GEN p)
{
  if (lgefint(p)==3)
    return Flx_to_ZX(Flx_factcyclo(n, p[2], 1));
  else
    return FpX_factcyclo(n, p, 1);
}

static void
Flx_red_inplace(GEN x, ulong p)
{
  long i, l = lg(x);
  for (i=2; i<l; i++) x[i] = uel(x, i)%p;
  Flx_renormalize(x, l);
}

/* x[i], T[i] < pn */
static GEN
Flxq_xi_conj(GEN x, GEN T, long j, long d, long pn)
{
  long i, deg = degpol(x);
  GEN z = const_vecsmall(d+1, 0);
  for (i=0; i<=deg; i++) z[2+Fl_mul(i, j, d)] = x[2+i];
  return Flx_rem(Flx_renormalize(z, d+2), T, pn);
}

static GEN
FlxqX_xi_conj(GEN x, GEN T, long j, long d, long pn)
{
  long i, l = lg(x);
  GEN z;
  z = cgetg(l, t_POL); z[1] = evalsigne(1) | evalvarn(0);
  for (i=2; i<l; i++) gel(z, i) = Flxq_xi_conj(gel(x, i), T, j, d, pn);
  return z;
}

static GEN
FlxqX_xi_norm(GEN x, GEN T, long p, long d, long pn)
{
  long i, d_chi = degpol(T);
  GEN z = x, z1 = x;
  for (i=1; i<d_chi; i++)
  {
    z1 = FlxqX_xi_conj(z1, T, p, d, pn);
    z = FlxqX_mul(z, z1, T, pn);
  }
  return z;
}

/* assume 0 <= x[i], y[j] <= m-1 */
static GEN
FpV_shift_add(GEN x, GEN y, GEN m, long start, long end)
{
  long i, j;
  for (i=start, j=1; i<=end; i++, j++)
  {
    pari_sp av = avma;
    GEN z = addii(gel(x, i), gel(y, j));
    gel(x, i) = (cmpii(z, m) >= 0)? gerepileuptoint(av, subii(z, m)): z;
  }
  return x;
}

/* assume 0 <= x[i], y[j] <= m-1 */
static GEN
FpV_shift_sub(GEN x, GEN y, GEN m, long start, long end)
{
  long i, j;
  for (i=start, j=1; i<=end; i++, j++)
  {
    pari_sp av = avma;
    GEN z = subii(gel(x, i), gel(y, j));
    gel(x, i) = (signe(z) < 0)? gerepileuptoint(av, addii(z, m)): z;
  }
  return x;
}

/* assume 0 <= x[i], y[j] <= m-1 */
static GEN
Flv_shift_add(GEN x, GEN y, ulong m, long start, long end)
{
  long i, j;
  for (i=start, j=1; i<=end; i++, j++)
  {
    ulong xi = x[i], yj = y[j];
    x[i] = Fl_add(xi, yj, m);
  }
  return x;
}

/* assume 0 <= x[i], y[j] <= m-1 */
static GEN
Flv_shift_sub(GEN x, GEN y, ulong m, long start, long end)
{
  long i, j;
  for (i=start, j=1; i<=end; i++, j++)
  {
    ulong xi = x[i], yj = y[j];
    x[i] = Fl_sub(xi, yj, m);
  }
  return x;
}

/* return 0 if p|x. else return 1 */
INLINE long
Flx_divcheck(GEN x, ulong p)
{
  long i, l = lg(x);
  for (i=2; i<l; i++) if (uel(x, i)%p) return 1;
  return 0;
}

static long
FlxX_weier_deg(GEN pol, long p)
{
  long i, l = lg(pol);
  for (i=2; i<l && Flx_divcheck(gel(pol, i), p)==0; i++);
  return (i<l)?i-2:-1;
}

static long
Flx_weier_deg(GEN pol, long p)
{
  long i, l = lg(pol);
  for (i=2; i<l && pol[i]%p==0; i++);
  return (i<l)?i-2:-1;
}

static GEN
Flxn_shift_mul(GEN g, long n, GEN p, long d, long m)
{
  return Flx_shift(Flxn_mul(g, p, d, m), n);
}

INLINE long
deg_trunc(long lam, long p, long n, long pn)
{
  long r, x, d;
  for (r=1,x=p; x<lam; r++) x *= p;  /* r is min int s.t. lam<=p^r */
  if ((d = (n-r+2)*lam+1)>=pn) d = pn;
  return d;
}

/*  Flx_translate1_basecase(g, pn) becomes slow when degpol(g)>1000.
 *  So I wrote Flxn_translate1().
 *  I need lambda to truncate pol.
 *  But I need to translate T --> 1+T to know lambda.
 *  Though the code has a little overhead, it is still fast. */
static GEN
Flxn_translate1(GEN g, long p, long n)
{
  long i, j, d, lam, pn, start;
  if (n==1) start = 3;
  else if (n==2) start = 9;
  else start = 10;
  pn = upowuu(p, n);
  for (lam=start; lam; lam<<=1)  /* least upper bound is 3 */
  {
    GEN z;
    d = deg_trunc(lam, p, n, pn);
    z = const_vecsmall(d+1, 0);  /* z[2],...,z[d+1] <--> a_0,...,a_{d-1} */
    for (i=degpol(g); i>=0; i--)
    {
      for (j=d+1; j>2; j--) z[j] = Fl_add(z[j], z[j-1], pn);  /* z = z*(1+T) */
      z[2] = Fl_add(z[2], g[2+i], pn);
    }
    z = Flx_renormalize(z, d+2);
    if (Flx_weier_deg(z, p) <= lam) return z;
  }
  return NULL; /*LCOV_EXCL_LINE*/
}

static GEN
FlxXn_translate1(GEN g, long p, long n)
{
  long i, j, d, lam, pn, start;
  GEN z;
  if (n==1) start = 3;
  else if (n==2) start = 9;
  else start = 10;
  pn = upowuu(p, n);
  for (lam=start; lam; lam<<=1)  /* least upper bound is 3 */
  {
    d = deg_trunc(lam, p, n, pn);
    z = const_vec(d+1, pol0_Flx(0));  /* z[2],...,z[d+1] <--> a_0,...,a_{d-1} */
    settyp(z, t_POL); z[1] = evalsigne(1) | evalvarn(0);
    for (i=degpol(g); i>=0; i--)
    {
      for (j=d+1; j>2; j--) gel(z, j) = Flx_add(gel(z, j), gel(z, j-1), pn);
      gel(z, 2) = Flx_add(gel(z, 2), gel(g, 2+i), pn);
    }
    z = FlxX_renormalize(z, d+2);
    if (FlxX_weier_deg(z, p) <= lam) return z;
  }
  return NULL; /*LCOV_EXCL_LINE*/
}

/* lam < 0 => error (lambda can't be determined)
 * lam = 0 => return 1
 * lam > 0 => return dist. poly. of degree lam. */
static GEN
Flxn_Weierstrass_prep(GEN g, long p, long n, long d_chi)
{
  long i, r0, d, dg = degpol(g), lam, pn, t;
  ulong lam0;
  GEN U, UINV, P, PU, g0, g1, gp, gU;
  if ((lam = Flx_weier_deg(g, p))==0) return(pol1_Flx(0));
  else if (lam<0)
    pari_err(e_MISC,"Flxn_Weierstrass_prep: precision too low. Increase n!");
  lam0 = lam/d_chi;
  pn = upowuu(p, n);
  d = deg_trunc(lam, p, n, pn);
  if (d>dg) d = dg;
  if (d<=lam) d=1+lam;
  for (r0=1; upowuu(p, r0)<lam0; r0++);
  g = Flxn_red(g, d);
  t = Fl_inv(g[2+lam], pn);
  g = Flx_Fl_mul(g, t, pn);  /* normalized so as g[2+lam]=1 */
  U = Flx_shift(g, -lam);
  UINV = Flxn_inv(U, d, pn);
  P = zx_z_divexact(Flxn_red(g, lam), p);  /* assume g[i] <= LONG_MAX */
  PU = Flxn_mul(P, UINV, d, pn);
  gU = Flxn_mul(g, UINV, d, pn);
  g0 = pol1_Flx(0);
  g1 = pol1_Flx(0);
  for (i=1; i<n; i++)
  {
    g1 = Flxn_shift_mul(g1, -lam, PU, d, pn);
    gp = Flx_Fl_mul(g1, upowuu(p, i), pn);
    g0 = (i&1)?Flx_sub(g0, gp, pn):Flx_add(g0, gp, pn);
  }
  g0 = Flxn_mul(g0, gU, lam+1, pn);
  g0 = Flx_red(g0, upowuu(p, (p==2)?n-r0:n+1-r0));
  return g0;
}

/* xi_n and Iwasawa pol. for Q(sqrt(m)) and p
 *
 * (flag&1)!=0 ==> output xi_n
 * (flag&2)!=0 ==> output power series
 * (flag&4)!=0 ==> output Iwasawa polynomial */
static GEN
imagquadstkpol(long p, long m, long n)
{
  long pn = upowuu(p, n);
  GEN pol, stk, stk2;
  int *chi;
  if (p==2 && (m==-1 || m==-2 || m==-3 || m==-6)) return nullvec();
  if (p==3 && m==-3) return nullvec();
  if (p==2 && m%2==0) m /= 2;
  chi = set_quad_chi_2(m);
  stk = (p==2)? quadstk2(m, n, chi): quadstkp(p, m, n, chi);
  stk2 = zx_to_Flx(stk, pn);
  pol = Flxn_Weierstrass_prep(zlx_translate1(stk2, p, n), p, n, 1);
  return degpol(pol)? mkvec(Flx_to_ZX(pol)): nullvec();
}

/* a mod p == g^i mod p ==> omega(a)=zeta_(p-1)^(-i)
 *  Chi[g^i mod p]=i (0 <= i <= p-2) */
static GEN
get_teich(long p, long g)
{
  long i, gi = 1, p1 = p-1;
  GEN Chi = cgetg(p, t_VECSMALL);
  for (i=0; i<p1; i++) { Chi[gi] = i; gi = Fl_mul(gi, g, p); }
  return Chi;
}

/* Ichimura-Sumida criterion for Greenberg conjecture for real quadratic field.
 * chi: character of Q(sqrt(m)), omega: Teichmuller character mod p or 4.
 * Get Stickelberger element from chi^* = omega*chi^(-1) and convert it to
 * power series by the correspondence (Q_n/Q,1+q0)^(-1) <-> (1+T)(1+q0)^(-1) */
static GEN
realquadstkpol(long p, long m, long n)
{
  int *chi;
  long pnm1 = upowuu(p, n-1),pn = p*pnm1, pn1 = p*pn, d, q0;
  GEN stk, ser, pol;
  if (m==1) pari_err_DOMAIN("quadstkpol", "m", "=", gen_1, gen_1);
  if (p==2 && (m&1)==0) m>>=1;
  d = ((m&3)==1)?m:m<<2;
  q0 = ulcm((p==2)?4:p, d);
  if (p==2)
  {
    chi = set_quad_chi_2(-m);
    stk = quadstk2(-m, n, chi);
    stk = zx_to_Flx_inplace(stk, pn);
  }
  else if (p==3 && m%3==0 && kross(-m/3,3)==1)
  {
    long m3 = m/3;
    chi = set_quad_chi_2(-m3);
    stk = quadstkp(3, -m3, n, chi);
    stk = zx_to_Flx_inplace(stk, pn);
  }
  else
  {
    long g = pgener_Zl(p);
    long x = Fl_powu(Fl_inv(g, p), pnm1, pn);
    GEN Chi = get_teich(p, g);
    GEN Gam = set_gam((1+q0)%pn1, p, n);
    chi = set_quad_chi_2(m);
    stk = quadteichstk(Chi, chi, Gam, p, m, n);  /* exact */
    stk = zxX_to_FlxX(stk, pn);  /* approx. */
    stk = FlxY_evalx(stk, x, pn);
  }
  stk = Flx_rescale_inplace(Flx_recip1_inplace(stk, pn), (1+q0)%pn, pn);
  ser = Flxn_translate1(stk, p, n);
  pol = Flxn_Weierstrass_prep(ser, p, n, 1);
  return degpol(pol)? mkvec(Flx_to_ZX(pol)): nullvec();
}

/* m > 0 square-free. lambda_2(Q(sqrt(-m)))
 * Kida, Tohoku Math. J. vol.31 (1979), 91-96, Theorem 1. */
static GEN
quadlambda2(long m)
{
  long i, l, L;
  GEN P;
  if ((m&1)==0) m >>= 1;  /* lam_2(Q(sqrt(-m)))=lam_2(Q(sqrt(-2*m))) */
  if (m <= 3) return mkvecs(0);
  P = gel(factoru(m), 1); l = lg(P);
  for (L = -1,i = 1; i < l; i++) L += 1L << (-3 + vals(P[i]-1) + vals(P[i]+1));
  return mkvecs(L);
}

/* Iwasawa lambda invariant of Q(sqrt(m)) (m<0) for p
 * |A_n|=p^(e[n])
 * kross(m,p)!=1 : e[n]-e[n-1]<eulerphi(p^n)  ==> lambda=e[n]-e[n-1]
 * kross(m,p)==1 : e[n]-e[n-1]<=eulerphi(p^n) ==> lambda=e[n]-e[n-1]
 * Gold, Acta Arith. XXVI (1974), p.25, Cor. 3
 * Gold, Acta Arith. XXVI (1975), p.237, Cor. */
static GEN
quadlambda(long p, long m)
{
  long flag, n, phipn;
  GEN e = cgetg(31, t_VECSMALL);
  int *chi;
  if (m>0) pari_err_IMPL("plus part of lambda invariant in quadlambda()");
  if (p==2) return quadlambda2(-m);
  if (p==3 && m==-3) return mkvec3(gen_0, gen_0, nullvec());
  flag = kross(m, p);
  e[1] = Z_lval(quadclassno(quaddisc(stoi(m))), p);
  if (flag!=1 && e[1]==0) return mkvec3(gen_0, gen_0, nullvec());
  chi = set_quad_chi_2(m);
  phipn = p-1;  /* phipn=phi(p^n) */
  for (n=1; n; n++, phipn *= p)
  {
    long L = flag? ediff(p, m, n, chi): ediff_ber(p, m, n, chi);
    e[n+1] = e[n] + L;
    if ((flag!=1 && (L < phipn))|| (flag==1 && (L <= phipn))) break;
  }
  return vecsmall2vec2(e, n);
}

/* factor n-th cyclotomic polynomial mod p^r and return a minimal
 *  polynomial of zeta_n over Q_p.
 *  phi(n)=deg*n_conj, n_conj == 1 <=> polcyclo(n) is irred mod p. */
static GEN
set_minpol(ulong n, GEN p, ulong r, long n_conj)
{
  GEN z, v, pol, pr;
  pari_timer ti;
  if (umodiu(p, n)==1) /* zeta_n in Z_p, faster than polcyclo() */
  {
    GEN prm1 = powiu(p, r-1), pr = mulii(prm1, p); /* pr=p^r */
    GEN prn = diviuexact(subii(pr, prm1), n);      /* prn=phi(p^r)/n */
    z = Fp_pow(pgener_Fp(p), prn, pr);
    return deg1pol_shallow(gen_1, Fp_neg(z, pr), 0);
  }
  pr = powiu(p, r);
  pol = polcyclo(n, 0);
  if (n_conj==1) return FpX_red(pol, pr);
  if (DEBUGLEVEL>3) timer_start(&ti);
  z = FpX_one_cyclo(n, p);
  if (DEBUGLEVEL>3) timer_printf(&ti, "FpX_one_cyclo:n=%ld  ", n);
  v = ZpX_liftfact(pol, mkvec2(z, FpX_div(pol, z, p)), pr, p, r);
  return gel(v, 1);
}

static GEN
set_minpol_teich(ulong g_K, GEN p, ulong r)
{
  GEN prm1 = powiu(p, r-1), pr = mulii(prm1, p), z;
  z = Fp_pow(Fp_inv(utoi(g_K), p), prm1, pr);
  return deg1pol_shallow(gen_1, Fp_neg(z, pr), 0);
}

static long
srh_1(GEN H)
{
  GEN bits = gel(H, 3);
  ulong f = bits[1];
  return F2v_coeff(bits, f-1);
}

/* (1/f)sum_{1<=a<=f}a*chi^{-1}(a) = -(1/(2-chi(a)))sum_{1<=a<=f/2} chi^{-1}(a)
 *  does not overflow */
static GEN
zx_ber_num(GEN Chi, long f, long d)
{
  long i, f2 = f>>1;
  GEN x = const_vecsmall(d+1, 0), x2 = x+2;
  for (i = 1; i <= f2; i++)
    if (Chi[i] >= 0) x2[Chi[i]] ++;
  return zx_renormalize(x, d+2);
}

/* x a zx
 * zx_ber_num is O(f). ZX[FpX,Flx]_ber_conj is O(d). Sometimes d<<f. */
static GEN
ZX_ber_conj(GEN x, long j, long d)
{
  long i, deg = degpol(x);
  GEN y = pol_zero(d), x2 = x+2, y2 = y+2;
  for (i=0; i<=deg; i++) gel(y2, Fl_mul(i, j, d)) = stoi(x2[i]);
  return ZX_renormalize(y, d+2);
}

/* x a zx */
static GEN
FpX_ber_conj(GEN x, long j, long d, GEN p)
{
  long i, deg = degpol(x);
  GEN y = pol_zero(d), x2 = x+2, y2 = y+2;
  for (i=0; i<=deg; i++) gel(y2, Fl_mul(i, j, d)) = modsi(x2[i], p);
  return FpX_renormalize(y, d+2);
}

/* x a zx */
static GEN
Flx_ber_conj(GEN x, long j, long d, ulong p)
{
  long i, deg = degpol(x);
  GEN y = const_vecsmall(d+1, 0), x2 = x+2, y2 = y+2;
  for (i=0; i<=deg; i++) y2[Fl_mul(i, j, d)] = umodsu(x2[i], p);
  return Flx_renormalize(y, d+2);
}

static GEN
ZX_ber_den(GEN Chi, long j, long d)
{
  GEN x = pol_zero(d), x2 = x+2;
  if (Chi[2]>=0) gel(x2, Fl_neg(Fl_mul(Chi[2], j, d), d)) = gen_1;
  gel(x2, 0) = subiu(gel(x2, 0), 2);
  return ZX_renormalize(x, d+2);
}

static GEN
Flx_ber_den(GEN Chi, long j, long d, ulong p)
{
  GEN x = const_vecsmall(d+1, 0), x2 = x+2;
  if (Chi[2]>=0) x2[Fl_neg(Fl_mul(Chi[2], j, d), d)] = 1;
  x2[0] = Fl_sub(x2[0], 2, p);
  return Flx_renormalize(x, d+2);
}

/* x is ZX of deg <= d-1 */
static GEN
ber_conj(GEN x, long k, long d)
{
  long i, deg = degpol(x);
  GEN z = pol_zero(d);
  if (k==1)
    for (i=0; i<=deg; i++) gel(z, 2+i) = gel(x, 2+i);
  else
    for (i=0; i<=deg; i++) gel(z, 2+Fl_mul(i, k, d)) = gel(x, 2+i);
  return ZX_renormalize(z, d+2);
}

/* The computation is fast when p^n and el=1+k*f*p^n are less than 2^64
 *  for m <= n <= M
 *  We believe M>=3 is enough when f%p=0 and M>=2 is enough for other case
 *  because we expect that p^2 does not divide |A_{K,psi}| for a large p.
 *  FIXME: M should be set according to p and f. */
static void
set_p_f(GEN pp, ulong f, long *pm, long *pM)
{
  ulong p = itou_or_0(pp);
  if (!p || p >= 2000000) { *pm=2; *pM = dvdui(f, pp)? 3: 2; }
  else if (p == 3)      { *pm=5; *pM=20; }
  else if (p == 5)      { *pm=5; *pM=13; }
  else if (p == 7)      { *pm=5; *pM=11; }
  else if (p == 11)     { *pm=5; *pM=9; }
  else if (p == 13)     { *pm=5; *pM=8; }
  else if (p < 400)     { *pm=5; *pM=7; }
  else if (p < 5000)    { *pm=3; *pM=5; }
  else if (p < 50000)   { *pm=2; *pM=4; }
  else                  { *pm=2; *pM=3; }
}

static GEN
subgp2ary(GEN H, long n)
{
  GEN v = gel(H, 3), w = cgetg(n+1, t_VECSMALL);
  long i, j, f = v[1];
  for (i = 1, j = 0; i <= f; i++)
    if (F2v_coeff(v,i)) w[++j] = i;
  return w;
}

static GEN
Flv_FlvV_factorback(GEN g, GEN x, ulong q)
{ pari_APPLY_ulong(Flv_factorback(g, gel(x,i), q)) }

/* lift chi character on G/H to character on G */
static GEN
zncharlift(GEN chi, GEN ncycGH, GEN U, GEN cycG)
{
  GEN nchi = char_normalize(chi, ncycGH);
  GEN c = ZV_ZM_mul(gel(nchi, 2), U), d = gel(nchi, 1);
  return char_denormalize(cycG, d, c);
}

/* 0 <= c[i] < d, i=1..r; (c[1],...,c[r], d) = 1; find e[i] such that
 * sum e[i]*c[i] = 1 mod d */
static GEN
Flv_extgcd(GEN c, ulong d)
{
  long i, j, u, f, l = lg(c);
  GEN e = zero_zv(l-1);
  if (l == 1) return e;
  for (f = d, i = 1; f != 1 && i < l; i++)
  {
    f = cbezout(f, itou(gel(c,i)), &u, &e[i]);
    if (!e[i]) continue;
    e[i] = umodsu(e[i], d);
    u = umodsu(u, d);
    if (u != 1) for (j = 1; j < i; j++) e[j] = Fl_mul(e[j], u, d);
  }
  return e;
}

/* f!=p; return exact xi. */
static GEN
get_xi_1(GEN Chi, GEN Gam, long p, long f, long n, long d, ulong pm)
{
  GEN Gam1 = Gam+1, xi;
  long i, j, j0, f0, pn, pn1, deg, pn1f;

  f0 = (f%p)?f:f/p;
  pn = upowuu(p, n); pn1 = p*pn; pn1f = pn1%f;
  xi = cgetg(pn+2, t_POL); xi[1] = evalsigne(1) | evalvarn(0);
  for (i=0; i<pn; i++) gel(xi, 2+i) = const_vecsmall(d+1, 0);
  for (j=1; j<pn1; j++)
  {
    long ipn1,*xij0;
    if ((j0 = Gam1[j])<0) continue;
    ipn1 = j%f; xij0 = gel(xi, 2+j0)+2;
    for (i=1; i<f0; i++)
    {
      if ((ipn1 += pn1f) >= f) ipn1 -= f;
      if (ipn1==0 || (deg = Chi[ipn1])<0) continue;
      xij0[deg] += i;
    }
  }
  for (i=0; i<pn; i++) Flx_red_inplace(gel(xi, 2+i), pm);
  return FlxX_renormalize(xi, pn+2);
}

/* f=p; return p^(n+1)*xi mod pm. */
static GEN
get_xi_2(GEN Chi, GEN Gam, long p, long f, long n, long d, ulong pm)
{
  long a, amodf, i, j0, pn, pn1, deg;
  GEN Gam1 = Gam+1, xi;

  pn = upowuu(p, n); pn1 = p*pn;
  xi = cgetg(pn+2, t_POL); xi[1] = evalsigne(1) | evalvarn(0);
  for (i=0; i<pn; i++) gel(xi, 2+i) = const_vecsmall(d+1, 0);
  for (a=1,amodf=0; a<pn1; a++)  /* xi is exact */
  {
    if (++amodf==f) amodf = 0;
    if ((j0=Gam1[a])<0 || amodf==0 || (deg=Chi[amodf])<0) continue;
    mael(xi, 2+j0, 2+deg) += a;
  }
  for (i=0; i<pn; i++) Flx_red_inplace(gel(xi, 2+i), pm);
  return FlxX_renormalize(xi, pn+2);
}

static GEN
pol_chi_xi(GEN K, long p, long j, long n)
{
  pari_sp av = avma;
  GEN MinPol2 = gel(K, 7), xi = gel(K, 8);
  long d = K_get_d(K), f = K_get_f(K), d_chi = K_get_dchi(K);
  long wd, minpolpow = (f==p)?2*n+1:n, pm = upowuu(p, minpolpow);
  GEN ser, pol, xi_conj;
  pari_timer ti;

  /* xi is FlxX mod p^m, MinPol2 is Flx mod p^m, xi_conj is FlxqX. */
  xi_conj = FlxqX_xi_conj(xi, MinPol2, j, d, pm);
  if (d_chi==1)  /* d_chi==1 if f==p */
  {
    xi_conj = FlxX_to_Flx(xi_conj);
    if (f==p) xi_conj = zx_z_divexact(xi_conj, upowuu(p, n+1));
  }
  /* Now xi_conj is mod p^n */
  if (DEBUGLEVEL>1) timer_start(&ti);
  ser = (d_chi==1) ? Flxn_translate1(xi_conj, p, n)
    : FlxXn_translate1(xi_conj, p, n);
  if (DEBUGLEVEL>1) timer_printf(&ti, "Flx%sn_translate1",(d_chi==1)?"":"X");
  wd = (d_chi==1)?Flx_weier_deg(ser, p):FlxX_weier_deg(ser, p);
  if (wd<0) pari_err(e_MISC,"pol_chi_xi: precision too low. Increase n!\n");
  else if (wd==0) return pol_1(0);
  /* wd>0, convert to dist. poly. */
  if (d_chi>1)  /* f!=p. minpolpow==n */
  {
    ser = FlxqX_xi_norm(ser, MinPol2, p, d, upowuu(p, n));
    ser = FlxX_to_Flx(ser);
  }
  pol = Flx_to_ZX(Flxn_Weierstrass_prep(ser, p, n, d_chi));
  setvarn(pol, fetch_user_var("T"));
#ifdef DEBUG
  if (wd>0 && d_chi>1)
    err_printf("(wd,d_chi,p,f,d,j,H)=(%ld,%ld,%ld,%ld,%ld,%ld,%Ps)\n",
        wd,d_chi,p,f,d,j,gmael3(K, 1, 1, 1));
#endif
  return gerepilecopy(av, pol);
}

/* return 0 if lam_psi (psi=chi^j) is determined to be zero.
 * else return -1.
 * If psi(p)!=1, then N_{Q(zeta_d)/Q}(1-psi(p))!=0 (mod p) */
static long
lam_chi_ber(GEN K, long p, long j)
{
  pari_sp av = avma;
  GEN B1, B2, Chi = gel(K, 2), MinPol2 = gel(K, 7), B_num = gel(K, 8);
  long x, p2 = p*p, d = K_get_d(K), f = K_get_f(K);

  if (f == d+1 && p == f && j == 1) return 0;  /* Teichmuller */

  B1 = Flx_rem(Flx_ber_conj(B_num, j, d, p2), MinPol2, p2);
  B2 = Flx_rem(Flx_ber_den(Chi, j, d, p2), MinPol2, p2);
  if (degpol(B1)<0 || degpol(B2)<0)
    return gc_long(av, -1); /* 0 mod p^2 */
  x = zx_lval(B1, p) - zx_lval(B2, p);
  if (x<0) pari_err_BUG("subcycloiwasawa [Bernoulli number]");
  return gc_long(av, x==0 ? 0: -1);
}

static long
lam_chi_xi(GEN K, long p, long j, long n)
{
  pari_sp av = avma;
  GEN xi_conj, z, MinPol2 = gel(K, 7), xi = gel(K, 8);
  long d = K_get_d(K), f = K_get_f(K), d_chi = K_get_dchi(K);
  long wd, minpolpow = (f==p)?n+2:1, pm = upowuu(p, minpolpow);

  /* xi is FlxX mod p^m, MinPol2 is Flx mod p^m, xi_conj is FlxqX. */
  xi_conj = FlxqX_xi_conj(xi, MinPol2, j, d, pm);
  if (d_chi==1)  /* d_chi==1 if f==p */
  {
    xi_conj = FlxX_to_Flx(xi_conj);
    if (f==p) xi_conj = zx_z_divexact(xi_conj, upowuu(p, n+1));
  }
  /* Now xi_conj is mod p^n */
  z = (d_chi==1) ? Flxn_translate1(xi_conj, p, n)
    : FlxXn_translate1(xi_conj, p, n);
  wd = (d_chi==1)?Flx_weier_deg(z, p):FlxX_weier_deg(z, p);
#ifdef DEBUG
  if (wd>0 && d_chi>1)
    err_printf("(wd,d_chi,p,f,d,j,H)=(%ld,%ld,%ld,%ld,%ld,%ld,%Ps)\n",
        wd,d_chi,p,f,d,j,gmael3(K, 1, 1, 1));
#endif
  return gc_long(av, wd<0 ? -1 : wd*d_chi);
}

/* K = [H1, Chi, Minpol, C, [d_chi, n_conj]] */
static GEN
imag_cyc_pol(GEN K, long p, long n)
{
  pari_sp av = avma;
  GEN Chi = gel(K, 2), MinPol = gel(K, 3), C = gel(K, 4), MinPol2;
  long d_K = K_get_d(K), f = K_get_f(K), n_conj = K_get_nconj(K);
  long i, q0, pn1, pM, pmodf = p%f, n_done = 0;
  GEN z = nullvec(), Gam, xi, Lam, K2;

  Lam = const_vecsmall(n_conj, -1);
  if (pmodf==0 || Chi[pmodf]) /* mark trivial chi-part using Bernoulli number */
  {
    MinPol2 = ZX_to_Flx(MinPol, p*p);  /* p^2 for B_{1,chi} */
    K2 = shallowconcat(K, mkvec2(MinPol2, zx_ber_num(Chi, f, d_K)));
    for (i=1; i<=n_conj; i++)
      if ((Lam[i] = lam_chi_ber(K2, p, C[i])) == 0) n_done++;
    if (n_conj==n_done) return gerepilecopy(av, z); /* all chi-parts trivial */
  }
  q0 = (f%p)? f*p: f;
  pn1 = upowuu(p, n+1);
  Gam = set_gam((1+q0)%pn1, p, n);
  pM = upowuu(p, (f==p)? 2*n+1: n);
  MinPol2 = ZX_to_Flx(MinPol, pM);
  xi = (f==p)? get_xi_2(Chi, Gam, p, f, n, d_K, pM)
             : get_xi_1(Chi, Gam, p, f, n, d_K, pM);
  K2 = shallowconcat(K, mkvec2(MinPol2, xi));
  for (i=1; i<=n_conj; i++)
  {
    GEN z1;
    if (Lam[i]>=0) continue;
    z1 = pol_chi_xi(K2, p, C[i], n);
    if (degpol(z1)) z = vec_append(z, z1);  /* degpol(z1) may be zero */
  }
  return gerepilecopy(av, z);
}

/* K is an imaginary cyclic extension of degree d contained in Q(zeta_f)
 * H is the subgr of G=(Z/fZ)^* corresponding to K
 * h=|H|, d*h=phi(f)
 * G/H=<g> i.e. g^d \in H
 * d_chi=[Q_p(zeta_d):Q_p], i.e. p^d_chi=1 (mod d)
 * An inj. char. of G(K/Q) is automatically imaginary.
 *
 * G(K/Q)=G/H=<g>, chi:G(K/Q) -> overline{Q_p} s.t. chi(g)=zeta_d^(-1)
 * Chi[a]=k, k<0  => chi(a)=0
 *           k>=0 => chi(a)=zeta_d^(-k)
 * psi=chi^j, j in C : repre. of inj. odd char.
 * psi(p)==1 <=> chi(p)^j==0 <=> j*Chi[p]=0 (mod d) <=> Chi[p]==0
 *
 * K = [H1, Chi, Minpol, C, [d_chi, n_conj]] */
static long
imag_cyc_lam(GEN K, long p)
{
  pari_sp av = avma;
  GEN Chi = gel(K, 2), MinPol = gel(K, 3), C = gel(K, 4), MinPol2;
  long d_K = K_get_d(K), f = K_get_f(K), n_conj = K_get_nconj(K);
  long i, q0, n, pmodf = p%f, n_done = 0;
  ulong pn1, pM;
  GEN p0 = utoi(p), Gam, Lam, xi, K2;

  q0 = (f%p)? f*p: f;
  Lam = const_vecsmall(n_conj, -1);
  if (pmodf==0 || Chi[pmodf])  /* 1st trial is Bernoulli number */
  {
    MinPol2 = ZX_to_Flx(MinPol, p*p);  /* p^2 for B_{1,chi} */
    K2 = shallowconcat(K, mkvec2(MinPol2, zx_ber_num(Chi, f, d_K)));
    for (i=1; i<=n_conj; i++)
      if ((Lam[i] = lam_chi_ber(K2, p, C[i])) == 0) n_done++;
    if (n_conj==n_done) return gc_long(av, 0);  /* all chi-parts trivial */
  }
  pM = pn1 = p;
  for (n=1; n>=0; n++)  /* 2nd trial is Stickelberger element */
  {
    pn1 *= p; /* p^(n+1) */
    if (f == p)
    { /* do not use set_minpol: it returns a new pol for each call */
      GEN fac, cofac, v, pol = polcyclo(d_K, 0);
      pM = pn1 * p; /* p^(n+2) */
      fac = FpX_red(MinPol, p0); cofac = FpX_div(pol, fac, p0);
      v = ZpX_liftfact(pol, mkvec2(fac, cofac), utoipos(pM), p0, n+2);
      MinPol2 = gel(v, 1);
    }
    Gam = set_gam((1+q0)%pn1, p, n);
    MinPol2 = ZX_to_Flx(MinPol, pM);
    xi = (f==p)? get_xi_2(Chi, Gam, p, f, n, d_K, pM)
               : get_xi_1(Chi, Gam, p, f, n, d_K, pM);
    K2 = shallowconcat(K, mkvec2(MinPol2, xi));
    for (i=1; i<=n_conj; i++)
      if (Lam[i]<0 && (Lam[i] = lam_chi_xi(K2, p, C[i], n)) >= 0) n_done++;
    if (n_conj==n_done) break;
  }
  return gc_long(av, zv_sum(Lam));
}
static GEN
GHinit(long f, GEN HH, GEN *pcycGH)
{
  GEN G = znstar0(utoipos(f), 1);
  GEN U, Ui, cycG, cycGH, ncycGH, gG, gGH, vChar, vH1, P, gH = gel(HH, 1);
  long i, expG, n_f, lgH = lg(gH); /* gens. of H */
  P = cgetg(lgH, t_MAT);
  for (i = 1; i < lgH; i++) gel(P,i) = Zideallog(G, utoi(gH[i]));

  /* group structure of G/H */
  cycG = znstar_get_cyc(G);
  expG = itou(gel(cycG, 1));
  /* gG generators of G, gGH generators of G/H: gGH = g.Ui, g = gGH.U */
  cycGH = ZM_snf_group(hnfmodid(P, cycG), &U, &Ui);
  ncycGH = cyc_normalize(cycGH);
  gG = ZV_to_Flv(znstar_get_gen(G), f); /* gens. of G */
  /* generators of G/H */
  gGH = Flv_FlvV_factorback(gG, ZM_to_Flm(Ui, expG), f);
  vChar = chargalois(cycGH, NULL); n_f = lg(vChar)-2;
  vH1 = cgetg(n_f+1, t_VEC);
  for (i = 1; i <= n_f; i++)
  { /* skip trivial character */
    GEN chi = gel(vChar,i+1), nchi = char_normalize(chi, ncycGH);
    GEN chiG, E, H1, C = gel(nchi, 2);
    long e, he, gen, d = itou(gel(nchi, 1));
    /* chi(prod g[i]^e[i]) = e(sum e[i]*C[i] / d), chi has order d = #(G/H1)*/
    E = Flv_extgcd(C, d); /* \sum C[i]*E[i] = 1 in Z/dZ */

    chiG = zncharlift(chi, ncycGH, U, cycG);
    H1 =  charker(cycG, chiG); /* H1 < G with G/H1 cyclic */
    e = itou( zncharconductor(G, chiG) ); /* cond H1 = cond chi */
    H1 = Flv_FlvV_factorback(zv_to_Flv(gG, e), ZM_to_Flm(H1, expG), e);
    gen = Flv_factorback(zv_to_Flv(gGH, e), E, e);
    H1 = znstar_generate(e, H1); /* G/H1 = <gen>, chi(gen) = e(1/d) */
    he = eulerphiu(e) / d;
    /* G/H1 = <gen> cyclic of index d, e = cond(H1) */
    gel(vH1,i) = mkvec3(H1, mkvecsmall5(d,e,he,srh_1(H1), gen),
                        subgp2ary(H1, he));
  }
  if (pcycGH) *pcycGH = cycGH;
  return vH1;
}

/* aH=g^iH ==> chi(a)=zeta_n^(-i); Chi[g^iH]=i; Chi[0] is never accessed */
static GEN
get_chi(GEN H1)
{
  GEN H = _get_H(H1);
  long i, j, gi, d = _get_d(H1), f = _get_f(H1), h = _get_h(H1), g = _get_g(H1);
  GEN Chi = const_vecsmall(f-1, -1);

  for (j=1; j<=h; j++) Chi[H[j]] = 0; /* i = 0 */
  for (i = 1, gi = g; i < d; i++)
  {
    for (j=1; j<=h; j++) Chi[Fl_mul(gi, H[j], f)] = i;
    gi = Fl_mul(gi, g, f);
  }
  return Chi;
}

static void
errpdiv(const char *f, GEN p, long d)
{
  pari_err_DOMAIN(f, "p", "divides",
                  strtoGENstr(stack_sprintf("[F:Q] = %ld", d)), p);
}
/* p odd doesn't divide degF; return lambda invariant if n==0 and
 * iwasawa polynomials if n>=1 */
static GEN
abeliwasawa(long p, long f, GEN HH, long degF, long n)
{
  long lam = 0, i, n_f;
  GEN vH1, vData, z = nullvec(), p0 = utoi(p) ;

  vH1 = GHinit(f, HH, NULL); n_f = lg(vH1)-1;
  vData = const_vec(degF, NULL);
  for (i=1; i<=n_f; i++) /* prescan. set Teichmuller */
  {
    GEN H1 = gel(vH1,i);
    long d_K = _get_d(H1), f_K = _get_f(H1), g_K = _get_g(H1);

    if (f_K == d_K+1 && p == f_K)  /* found K=Q(zeta_p) */
    {
      long d_chi = 1, n_conj = eulerphiu(d_K);
      GEN C = set_C(p, d_K, d_chi, n_conj);
      long minpow = n? 2*n+1: 2;
      GEN MinPol = set_minpol_teich(g_K, p0, minpow);
      gel(vData, d_K) = mkvec4(MinPol, C, NULL, mkvecsmall2(d_chi, n_conj));
      break;
    }
  }

  for (i=1; i<=n_f; i++)
  {
    GEN H1 = gel(vH1,i), z1, Chi, K;
    long d_K = _get_d(H1), s = _get_s(H1);

    if (s) continue;  /* F is real */
#ifdef DEBUG
    err_printf("  handling %s cyclic subfield K, deg(K)=%ld, cond(K)=%ld\n",
        s? "a real": "an imaginary", d_K, _get_f(H1));
#endif
    if (!gel(vData, d_K))
    {
      long d_chi = order_f_x(d_K, p), n_conj = eulerphiu(d_K)/d_chi;
      GEN C = set_C(p, d_K, d_chi, n_conj);
      long minpow = n? n+1: 2;
      GEN MinPol = set_minpol(d_K, p0, minpow, n_conj);
      gel(vData, d_K) = mkvec4(MinPol, C, NULL, mkvecsmall2(d_chi, n_conj));
    }
    Chi = get_chi(H1);
    K = shallowconcat(mkvec2(H1, Chi), gel(vData, d_K));
    if (n==0) lam += imag_cyc_lam(K, p);
    else if (lg(z1 = imag_cyc_pol(K, p, n)) > 1) z = shallowconcat(z, z1);
  }
  return n? z: mkvecs(lam);
}

static GEN
ary2mat(GEN x, long n)
{
  long i, j;
  GEN z = cgetg(n+1,t_MAT);
  for (i=1; i<=n; i++)
  {
    gel(z,i) = cgetg(n+1,t_COL);
    for (j=1; j<=n; j++) gmael(z, i, j) = utoi(x[(i-1)*n+j-1]);
  }
  return z;
}

static long
is_cyclic(GEN x)
{
  GEN y = gel(x, 2);
  long i, l = lg(y), n = 0;
  for (i = 1; i < l; i++) if (signe(gel(y,i))) n++;
  return n <= 1;
}

static GEN
make_p_part(GEN y, ulong p, long d_pow)
{
  long i, l = lg(y);
  GEN z = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) z[i] = signe(gel(y,i))? Z_lval(gel(y,i), p): d_pow;
  return z;
}

static GEN
structure_MLL(GEN y, long d_pow)
{
  long y0, i, l = lg(y);
  GEN x = gen_0, E = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    if ((y0 = d_pow-y[i]) < 0) y0 = 0;
    x = addiu(x, y0);
    gel(E, l-i) = utoi(y0);
  }
  return mkvec2(x, E);
}

static long
find_del_el(GEN *oldgr, GEN newgr, long n, long n_el, long d_chi)
{
  if (n_el==1) return 1;
  if (equalis(gmael(newgr, 2, 1), n_el)) return n;
  if (cmpii(gel(*oldgr, 1), gel(newgr, 1)) >= 0) return n;
  if (n > 1 && is_cyclic(newgr)) { *oldgr = newgr; return n-1; }
  if (n == n_el) return n;
  if (cmpis(gel(newgr, 1), n*d_chi) < 0) return n;
  return 0;
}

static GEN
subgr2vecsmall(GEN H, long h, long f)
{
  long i;
  GEN z = const_vecsmall(f-1, 0); /* f=lg(z) */
  for (i=1; i<=h; i++) z[H[i]] = 1; /* H[i]!=0 */
  return z;
}

/* K is the subfield of Q(zeta_f) with degree d corresponding to the subgroup
 * H in (Z/fZ)^*; for a divisor e of f, zeta_e \in K <=> H \subset He. */
static long
root_of_1(long f, GEN H)
{
  GEN g = gel(H, 1); /* generators */
  long e = f, i, l = lg(g);
  for (i = 1; i < l; i++)
  {
    e = cgcd(e, g[i] - 1);
    if (e <= 2) return 2;
  }
  return odd(e)? (e<<1): e;
}

static long
find_ele(GEN H)
{
  long i, f=lg(H);
  for (i=1; i<f; i++) if (H[i]) return i;
  return 0;
}

static void
delete_ele(GEN H, long j, long el)
{
  long f = lg(H), x = 1;
  do H[Fl_mul(j,x,f)] = 0;
  while ((x=Fl_mul(x,el,f))!=1);
}

static GEN
get_coset(GEN H, long h, long f, long el)
{
  long i, j, k = h/order_f_x(f, el);
  GEN H2, coset = const_vecsmall(k, 0);
  H2 = subgr2vecsmall(H, h, f);
  for (i=0; (j=find_ele(H2))>0; i++)
  {
    coset[1+i] = j;
    delete_ele(H2, j, el);
  }
  if (i != k) pari_err_BUG("failed to find coset\n");
  return coset;
}

static long
srh_pol(GEN xpows, GEN vn, GEN pols, long el, long k, long f)
{
  pari_sp av = avma;
  long i, j, l = lg(pols), d = degpol(gel(pols, 1));
  GEN pol = gel(pols, 1);

  for (i=1; i<l; i++)
  {
    GEN x, y, z;
    if (vn[i]==0) continue;
    y = gel(pols, vn[i]);
    z = pol0_Flx(0);
    for (j=0; j<=d; j++)
      z = Flx_add(z, Flx_Fl_mul(gel(xpows, 1+Fl_mul(j, k, f)), y[2+j], el), el);
    x = Flx_rem(z, pol, el);
    if (lg(x)==2)
    {vn[i] = 0; return gc_long(av, i); }  /* pols[i] is min pol of zeta_f^k */
  }
  pari_err_BUG("subcyclopclgp [minimal polinomial]");
  return 0;  /* to suppress warning */
}

/* e_chi[i mod dK] = chi(i*j), i = 0..dK-1; beware: e_chi is translated ! */
static GEN
get_e_chi(GEN K, ulong j, ulong d, ulong *pdK)
{
  ulong i, dK = K_get_d(K);
  GEN TR = gel(K,4) + 2, e_chi = cgetg(dK+1, t_VECSMALL) + 1;
  if (j == 1)
    for (i = 0; i < dK; i++) e_chi[i] = umodiu(gel(TR, i), d);
  else
    for (i = 0; i < dK; i++) e_chi[i] = umodiu(gel(TR, Fl_mul(i, j, dK)), d);
  *pdK = dK; return e_chi;
}
static GEN
get_e_chi_ll(GEN K, ulong j, GEN d)
{
  ulong i, dK = umael3(K, 1, 2, 1);
  GEN TR = gel(K,4) + 2, e_chi = cgetg(dK+1, t_VEC) + 1;
  for (i = 0; i < dK; i++) gel(e_chi,i) = modii(gel(TR, Fl_mul(i, j, dK)), d);
  return e_chi;
}

/* el=1 (mod f) */
static long
chk_el_real_f(GEN K, ulong p, ulong d_pow, ulong el, ulong j0)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong d_K, f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  ulong  i, j, gi, d = upowuu(p, d_pow), dp = d*p;
  ulong g_el, z_f, flag = 0, el1f = (el-1)/f, el1dp = (el-1)/dp;
  GEN e_chi = get_e_chi(K, j0, dp, &d_K);
  GEN vz_f, xi_el = cgetg(d_K+1, t_VECSMALL)+1;

  g_el = pgener_Fl(el);
  z_f = Fl_powu(g_el, el1f, el);
  vz_f = Fl_powers(z_f, f-1, el)+1;

  for (gi = 1, i = 0; i < d_K; i++)
  {
    ulong x = 1;
    for (j = 1; j <= h; j++)
    {
      ulong y = Fl_mul(H[j], gi, f);
      x = Fl_mul(x, vz_f[y]-1, el); /* vz_f[y] = z_f^y  */
    }
    gi = Fl_mul(gi, g_K, f);
    xi_el[i] = x;  /* xi_el[i] = xi^{g_K^i} mod el */
  }
  for (i=0; i<d_K; i++)
  {
    ulong x = 1;
    for (j=0; j<d_K; j++)
      x = Fl_mul(x, Fl_powu(xi_el[j], e_chi[(i+j)%d_K], el), el);
    if ((x = Fl_powu(x, el1dp, el))!=1) flag = 1;
    if (Fl_powu(x, p, el)!=1) return gc_long(av,0);
  }
  return gc_long(av, flag?1:0);
}

/* For a cyclic field K contained in Q(zeta_f),
 * computes minimal polynomial T of theta=Tr_{Q(zeta_f)/K}(zeta_f) over Q
 * and conjugates of theta */
static GEN
minpol_theta(GEN K)
{
  GEN HH = gmael3(K,1,1,1);
  return galoissubcyclo(utoi(K_get_f(K)), HH, 0, 0);
}

/*  xi[1+i] = i-th conj of xi = Tr_{Q(zeta_f)/K}(1-zeta_f).
 * |1-(cos(x)+i*sin(x))|^2 = 2(1-cos(x)) */
static GEN
xi_approx(GEN K, long prec)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  long d_K = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  GEN xi = cgetg(d_K+1, t_COL), vz_f = grootsof1(f, prec);
  long i, j, g = 1, h2 = h>>1;
  for (i=1; i<=d_K; i++)
  {
    GEN y = real_1(prec);
    for (j=1; j<=h2; j++)
    {
      GEN z = gmael(vz_f, 1+Fl_mul(H[j], g, f), 1);
      y = mulrr(y, shiftr(subsr(1, z), 1));
    }
    gel(xi, i) = y;
    g = Fl_mul(g, g_K, f);
  }
  return gerepilecopy(av, xi);
}

static GEN
theta_xi_el(GEN K, ulong el)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong d_K = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  GEN theta = cgetg(d_K+1, t_VECSMALL), xi = cgetg(d_K+1, t_VECSMALL), vz_f;
  ulong i, j, g = 1, x, y, g_el, z_f;

  g_el = pgener_Fl(el);
  z_f = Fl_powu(g_el, (el-1)/f, el);
  vz_f = Fl_powers(z_f, f-1, el);
  for (i=1; i<=d_K; i++)
  {
    x = 0; y = 1;
    for (j=1; j<=h; j++)
    {
      ulong z = vz_f[1+Fl_mul(H[j], g, f)];
      x = Fl_add(x, z, el);
      y = Fl_mul(y, z-1, el);
    }
    theta[i] = x;
    xi[i] = y;
    g = Fl_mul(g, g_K, f);
  }
  return gerepilecopy(av, mkvec2(theta, xi));
}

static GEN
make_Xi(GEN xi, long d)
{
  long i, j;
  GEN Xi = cgetg(d+1, t_MAT);
  for (j=0; j<d; j++)
  {
    GEN x = cgetg(d+1, t_VECSMALL);
    gel(Xi, 1+j) = x;
    for (i=0; i<d; i++) x[1+i] = xi[1+(i+j)%d];
  }
  return Xi;
}

static GEN
make_Theta(GEN theta, ulong d, ulong el)
{
  ulong i;
  GEN Theta = cgetg(d+1, t_MAT);
  for (i=1; i<=d; i++) gel(Theta, i) = Fl_powers(theta[i], d-1, el);
  return Flm_inv(Theta, el);
}

static GEN
Xi_el(GEN K, GEN tInvA, ulong el)
{
  pari_sp av = avma;
  ulong d_K = K_get_d(K);
  GEN tx = theta_xi_el(K, el), Theta, Xi, X;

  if ((Theta = make_Theta(gel(tx, 1), d_K, el))==NULL) return NULL;
  Xi = make_Xi(gel(tx, 2), d_K);
  X = Flm_mul(Flm_mul(Xi, Theta, el), ZM_to_Flm(tInvA, el), el);
  return gerepilecopy(av, X);
}

static GEN
pol_xi_el(GEN K, ulong el)
{
  pari_sp av = avma;
  ulong d_K = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  GEN H = K_get_H(K), xi = cgetg(d_K+1, t_VECSMALL), vz_f;
  ulong i, j, g = 1, y, g_el, z_f;

  g_el = pgener_Fl(el);
  z_f = Fl_powu(g_el, (el-1)/f, el);
  vz_f = Fl_powers(z_f, f-1, el);
  for (i=1; i<=d_K; i++)
  {
    y = 1;
    for (j=1; j<=h; j++)
    {
      ulong z = vz_f[1+Fl_mul(H[j], g, f)];
      y = Fl_mul(y, z-1, el);
    }
    xi[i] = y;
    g = Fl_mul(g, g_K, f);
  }
  return gerepilecopy(av, Flv_roots_to_pol(xi, el, 0));
}

/* theta[1+i] = i-th conj of theta; xi[1+i] = i-th conj of xi. */
static GEN
theta_xi_approx(GEN K, long prec)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  long d_K = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  GEN theta = cgetg(d_K+1, t_VEC), xi = cgetg(d_K+1, t_VEC);
  GEN vz_f = grootsof1(f, prec);
  long i, j, g = 1, h2 = h>>1;

  for (i=1; i<=d_K; i++)
  {
    GEN x = real_0(prec), y = real_1(prec);
    for (j=1; j<=h2; j++)
    {
      GEN z = gmael(vz_f, 1+Fl_mul(H[j], g, f), 1);
      x = addrr(x, z);
      y = mulrr(y, shiftr(subsr(1, z), 1));
    }
    gel(theta, i) = shiftr(x, 1);
    gel(xi, i) = y;
    g = Fl_mul(g, g_K, f);
  }
  return gerepilecopy(av, mkvec2(theta, xi));
}

static GEN
bound_coeff_xi(GEN K, GEN tInvA)
{
  pari_sp av = avma;
  long d_K = K_get_d(K), prec = MEDDEFAULTPREC, i;
  GEN tInvV, R = cgetg(d_K+1, t_MAT), theta_xi = theta_xi_approx(K, prec);
  GEN theta = gel(theta_xi, 1), xi = gel(theta_xi, 2), x1, x2, bound;

  for (i=1; i<=d_K; i++)
  {
    GEN z = gpowers(gel(theta, i), d_K-1);
    settyp(z, t_COL);
    gel(R, i) = z;
  }
  tInvV = RgM_mul(RgM_inv(R), tInvA);
  x1 = gsupnorm(tInvV, prec); x2 = gsupnorm(xi, prec);
  bound = mulrs(mulrr(x1, x2), 3*d_K);
  return gerepilecopy(av, bound);
}

static GEN
get_Xi(GEN K, GEN tInvA)
{
  pari_sp av = avma;
  GEN M0, XI, EL, Xi;
  ulong f = K_get_f(K), el, e, n, i;
  forprime_t T0;

  M0 = bound_coeff_xi(K, tInvA);
  e = expo(M0)+1; n = e/(BITS_IN_LONG-1); n++;
  EL = cgetg(1+n, t_VECSMALL);
  XI = cgetg(1+n, t_VEC);
  u_forprime_arith_init(&T0, LONG_MAX, ULONG_MAX, 1, f);

  for (i=1; i<=n; i++)
  {
    el = u_forprime_next(&T0);
    while ((Xi=Xi_el(K, tInvA, el))==NULL) el = u_forprime_next(&T0);
    gel(XI, i) = Xi;
    EL[i] = el;
  }
  return gerepileupto(av, nmV_chinese_center(XI, EL, NULL));
}

/* K is a cyclic field of conductor f with degree d=d_K
 * xi = Norm_{Q(zeta_f)/K}(1-zeta_f)
 * 1: T, min poly of a=Tr_{Q(zeta_f)/K}(zeta_f) over Q
 * 2: B, power basis of K with respect to a
 * 3: A, rational matrix s.t. t(v_1,...v_d) = A*t(1,a,...,a^{d-1})
 * 4: Xi, integer matrix s.t. t(xi^(1),...,xi^(d)) = Xi*t(v_1,...,v_d) */
static GEN
xi_data_basis(GEN K)
{
  pari_sp av = avma;
  GEN T = minpol_theta(K);
  GEN InvA, A, M, Xi, A_den;
  GEN D, B = nfbasis(T, &D);
  pari_timer ti;
  if (DEBUGLEVEL>1) timer_start(&ti);
  A = RgXV_to_RgM(B, lg(B)-1);
  M = gmael(A, 1, 1);
  if (!equali1(M)) A = RgM_Rg_div(A, M);
  InvA = QM_inv(A);
  A = Q_remove_denom(A, &A_den);
  if (A_den==NULL) A_den = gen_1;
  Xi = get_Xi(K, shallowtrans(InvA));
  if (DEBUGLEVEL>1) timer_printf(&ti, "xi_data_basis");
  return gerepilecopy(av, mkvec5(T, B, shallowtrans(A), Xi, A_den));
}

/* When factorization of polcyclo mod el is difficult, one can try to
 * check the condition of el using an integral basis of K.
 * This is useful when d_K is small. */
static long
chk_el_real_basis(GEN K, long p, long d_pow, long el, long j0)
{
  pari_sp av = avma;
  GEN xi = gel(K, 7), T = gel(xi, 1), A = gel(xi, 3), Xi = gel(xi, 4);
  GEN A_den = gel(xi, 5);
  ulong i, j, x, found = 0;
  GEN v_el, xi_el;
  GEN e_chi, xi_e_chi;
  ulong d_K, d, dp, el1dp;

  if (dvdiu(A_den, el)) return 0;

  d = upowuu(p, d_pow); dp = d*p; el1dp = (el-1)/dp;
  e_chi = get_e_chi(K, j0, dp, &d_K);
  xi_e_chi = cgetg(d_K+1, t_VECSMALL)+1;

  if (DEBUGLEVEL>1) err_printf("chk_el_real_basis: d_K=%ld el=%ld\n",d_K,el);
  A = ZM_to_Flm(A, el);
  A = Flm_Fl_mul(A, Fl_inv(umodiu(A_den, el), el), el);
  x = Flx_oneroot_split(ZX_to_Flx(T, el), el);
  v_el = Flm_Flc_mul(A, Fl_powers(x, d_K-1, el), el);
  xi_el = Flm_Flc_mul(ZM_to_Flm(Xi, el), v_el, el);
  if (DEBUGLEVEL>2) err_printf("el=%ld xi_el=%Ps\n", el, xi_el);
  for (i=0; i<d_K; i++)
  {
    long z = 1;
    for (j=0; j<d_K; j++)
      z = Fl_mul(z, Fl_powu(xi_el[1+j], e_chi[(i+j)%d_K], el), el);
    xi_e_chi[i] = z;
  }
  if (DEBUGLEVEL>2) err_printf("xi_e_chi=%Ps\n", xi_e_chi-1);
  for (i=0; i<d_K; i++)
  {
    long x = Fl_powu(xi_e_chi[i], el1dp, el);
    if (x!=1) found = 1;
    if (Fl_powu(x, p, el)!=1) return gc_long(av, 0);
  }
  return gc_long(av, found);
}

static GEN
bound_pol_xi(GEN K)
{
  pari_sp av = avma;
  GEN xi = xi_approx(K, MEDDEFAULTPREC);
  GEN M = real_1(MEDDEFAULTPREC), one = rtor(dbltor(1.0001), MEDDEFAULTPREC);
  long i, n = lg(xi);

  for (i=1; i<n; i++) M = mulrr(M, addrr(one, gel(xi, i)));
  M = mulrs(M, 3);
  return gerepilecopy(av, M);
}

static GEN
minpol_xi(GEN K)
{
  pari_sp av = avma;
  GEN M0, POL, EL;
  ulong f = K_get_f(K), el, e, n, i;
  forprime_t T0;

  M0 = bound_pol_xi(K);
  e = expo(M0)+1; n = e/(BITS_IN_LONG-1); n++;
  EL = cgetg(1+n, t_VECSMALL);
  POL = cgetg(1+n, t_VEC);
  u_forprime_arith_init(&T0, LONG_MAX, ULONG_MAX, 1, f);
  for (i=1; i<=n; i++)
  {
    el = u_forprime_next(&T0);
    gel(POL, i) = pol_xi_el(K, el);
    EL[i] = el;
  }
  return gerepileupto(av, nxV_chinese_center(POL, EL, NULL));
}

static long
find_conj_el(GEN K, GEN pol, GEN Den)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong d_K = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g = K_get_g(K);
  ulong j, k, el, g_el, z_f, xi = 1, xi_g = 1;
  GEN T = NULL, vz_f;

  for (el=f+1; el; el+=f)
    if (uisprime(el) && dvdiu(Den, el)==0)
    {
      T = ZX_to_Flx(pol, el);
      T = Flx_Fl_mul(T, Fl_inv(umodiu(Den, el), el), el);
      break;
    }
  g_el = pgener_Fl(el);
  z_f = Fl_powu(g_el, (el-1)/f, el);
  vz_f = Fl_powers(z_f, f-1, el);
  for (j=1; j<=h; j++)
    xi = Fl_mul(xi, vz_f[1+H[j]]-1, el);
  for (j=1; j<=h; j++)
    xi_g = Fl_mul(xi_g, vz_f[1+Fl_mul(H[j], g, f)]-1, el);
  for (k=1; k<=d_K; k++)
  {
    xi = Flx_eval(T, xi, el);
    if (xi == xi_g) break;
  }
  if (xi != xi_g) pari_err_BUG("find_conj_el");
  return gc_long(av, k);
}

/* G = H_1*H_2*...*H_m is cyclic of order n, H_i=<perms[i]>
 * G is not necessarily a direct product.
 * If p^e || n, then p^e || |H_i| for some i.
 * return a generator of G. */
static GEN
find_gen(GEN perms, long n)
{
  GEN fa = factoru(n), P = gel(fa, 1), E = gel(fa, 2);
  long i, j, l = lg(perms), r = lg(P);
  GEN gen = cgetg(r, t_VEC), orders = cgetg(l, t_VECSMALL), perm;

  for (i=1; i<l; i++) orders[i] = perm_orderu(gel(perms, i));
  for (i=1; i<r; i++)
  {
    long pe = upowuu(P[i], E[i]);
    for (j=1; j<l; j++) if (orders[j]%pe==0) break;
    gel(gen, i) = perm_powu(gel(perms, j), orders[j]/pe);
  }
  perm = gel(gen, 1);
  for (i=2; i<l; i++) perm = perm_mul(perm, gel(gen, i));
  return perm;
}

/* R is the roots of T. R[1+i] = R[1]^(g_K^i), 0 <= i <= d_K-1
 * 1: min poly T of xi over Q
 * 2: F(x)\in Q[x] s.t. xi^(g_K)=F(xi) */
static GEN
xi_data_galois(GEN K)
{
  pari_sp av = avma;
  GEN T, G, perms, perm, pol, pol2, Den;
  ulong k, d_K = K_get_d(K);
  pari_timer ti;

  if (DEBUGLEVEL>1) timer_start(&ti);
  T = minpol_xi(K);
  if (DEBUGLEVEL>1) timer_printf(&ti, "minpol_xi");
  G = galoisinit(T, NULL);
  if (DEBUGLEVEL>1) timer_printf(&ti, "galoisinit");
  perms = gal_get_gen(G);
  perm = (lg(perms)==2)?gel(perms, 1):find_gen(perms, d_K);
  if (DEBUGLEVEL>1) timer_start(&ti);
  pol = galoispermtopol(G, perm);
  if (DEBUGLEVEL>1) timer_printf(&ti, "galoispermtopol");
  pol = Q_remove_denom(pol, &Den);
  if (Den==NULL) Den = gen_1;
  k = find_conj_el(K, pol, Den);
  if (DEBUGLEVEL>1) timer_printf(&ti,"find_conj");
  pol2 = galoispermtopol(G, perm_powu(perm, k));
  pol2 = Q_remove_denom(pol2, &Den);
  if (Den==NULL) Den = gen_1;
  return gerepilecopy(av, mkvec3(T, pol2, Den));
}

/* If g(X)\in Q[X] s.t. g(xi)=xi^{g_K} was found,
 * then we fix an integer x_0 s.t. xi=x_0 (mod el) and construct x_i
 * s.t. xi^{g_K^i}=x_i (mod el) using g(X). */
static long
chk_el_real_galois(GEN K, long p, long d_pow, long el, long j0)
{
  pari_sp av = avma;
  GEN xi = gel(K, 7), T = gel(xi, 1), F = gel(xi, 2), Den = gel(xi, 3);
  GEN Fel, xi_el, xi_e_chi, e_chi;
  ulong i, j, x, found = 0;
  ulong d_K, d, dp, el1dp;

  if (dvdiu(Den, el)) return 0;

  d = upowuu(p, d_pow); dp = d*p; el1dp = (el-1)/dp;
  e_chi = get_e_chi(K, j0, dp, &d_K);
  xi_el = cgetg(d_K+1, t_VECSMALL)+1;
  xi_e_chi = cgetg(d_K+1, t_VECSMALL)+1;

  if (DEBUGLEVEL>1) err_printf("chk_el_real_galois: d_K=%ld el=%ld\n",d_K,el);
  Fel = ZX_to_Flx(F, el);
  Fel = Flx_Fl_mul(Fel, Fl_inv(umodiu(Den, el), el), el);
  x = Flx_oneroot_split(ZX_to_Flx(T, el), el);
  for (i=0; i<d_K; i++) { xi_el[i] = x; x = Flx_eval(Fel, x, el); }
  if (DEBUGLEVEL>2) err_printf("el=%ld xi_el=%Ps\n", el, xi_el-1);
  for (i=0; i<d_K; i++)
  {
    long z = 1;
    for (j=0; j<d_K; j++)
      z = Fl_mul(z, Fl_powu(xi_el[j], e_chi[(i+j)%d_K], el), el);
    xi_e_chi[i] = z;
  }
  if (DEBUGLEVEL>2) err_printf("xi_e_chi=%Ps\n", xi_e_chi-1);
  for (i=0; i<d_K; i++)
  {
    long x = Fl_powu(xi_e_chi[i], el1dp, el);
    if (x!=1) found = 1;
    if (Fl_powu(x, p, el)!=1) return gc_long(av, 0);
  }
  return gc_long(av, found);
}

/* checks the condition of el using the irreducible polynomial G_K(X) of zeta_f
 * over K. G_K(X) mod el is enough for our purpose and it is obtained by
 * factoring polcyclo(f) mod el */
static long
chk_el_real_factor(GEN K, long p, long d_pow, long el, long j0)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  ulong i, j, k, d_K, d = upowuu(p, d_pow), dp = d*p, found = 0;
  GEN pols, coset, vn_g, polnum, xpows, G_K;
  ulong el1dp = (el-1)/dp, n_coset, n_g, gi;
  GEN e_chi = get_e_chi(K, j0, dp, &d_K);
  pari_timer ti;

  if (DEBUGLEVEL>1) err_printf("chk_el_real_factor: f=%ld el=%ld\n",f,el);
  coset = get_coset(H, h, f, el);
  if (DEBUGLEVEL>1)
  {
    timer_start(&ti);
    err_printf("factoring polyclo(d) (mod %ld)\n",f, el);
  }
  pols = Flx_factcyclo(f, el, 0);
  if (DEBUGLEVEL>1) timer_printf(&ti,"Flx_factcyclo(%lu,%lu)",f,el);
  n_coset = lg(coset)-1;
  n_g = lg(pols)-1;
  vn_g = identity_perm(n_g);

  polnum = const_vec(d_K, NULL);
  for (i=1; i<=d_K; i++) gel(polnum, i) = const_vecsmall(n_coset, 0);
  xpows = Flxq_powers(polx_Flx(0), f-1, gel(pols, 1), el);
  for (gi=1,i=1; i<=d_K; i++)
  {
    for (j=1; j<=n_coset; j++)
    {
      long x, conj = Fl_mul(gi, coset[j], f);
      x = srh_pol(xpows, vn_g, pols, el, conj, f);
      gel(polnum, i)[j] = x;
    }
    gi = Fl_mul(gi, g_K, f);
  }
  G_K = const_vec(d_K, NULL);
  for (i=1; i<=d_K; i++)
  {
    GEN z = pol1_Flx(0);
    for (j=1; j<=n_coset; j++) z = Flx_mul(z, gel(pols, gel(polnum,i)[j]), el);
    gel(G_K, i) = z;
  }
  if (DEBUGLEVEL>2) err_printf("G_K(x)=%Ps\n",Flx_to_ZX(gel(G_K, 1)));
  for (k=0; k<d_K; k++)
  {
    long x = 1;
    for (i = 0; i < d_K; i++)
    {
      long x0, t;
      x0 = Flx_eval(gel(G_K, 1+i), 1, el);
      t = Fl_powu(x0, e_chi[(i+k)%d_K], el);
      x = Fl_mul(x, t, el);
    }
    x = Fl_powu(x, el1dp, el);
    if (x!=1) found = 1;
    if (Fl_powu(x, p, el)!=1) return gc_long(av, 0);
  }
  return gc_long(av, found);
}

static long
chk_el_real_chi(GEN K, ulong p, ulong d_pow, ulong el, ulong j0, long flag)
{
  ulong f = K_get_f(K);

  if (el%f == 1)
    return chk_el_real_f(K, p, d_pow, el, j0); /* must be faster */
  if (flag&USE_BASIS)
    return chk_el_real_basis(K, p, d_pow, el, j0);
  if (flag&USE_GALOIS_POL)
    return chk_el_real_galois(K, p, d_pow, el, j0);
  return chk_el_real_factor(K, p, d_pow, el, j0);
}

static long
chk_ell_real(GEN K, long d2, GEN ell, long j0)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  ulong d_K, i, j, gi;
  GEN e_chi = get_e_chi(K, j0, d2, &d_K);
  GEN g_ell, z_f, vz_f, xi_el = cgetg(d_K+1, t_VEC)+1;
  GEN ell_1 = subiu(ell,1), ell1d2 = diviuexact(ell_1, d2);

  g_ell = pgener_Fp(ell);
  z_f = Fp_pow(g_ell, diviuexact(subiu(ell, 1), f), ell);
  vz_f = Fp_powers(z_f, f-1, ell)+1;
  for (gi=1, i=0; i<d_K; i++)
  {
    GEN x = gen_1;
    for (j = 1; j <= h; j++)
    {
      ulong y = Fl_mul(H[j], gi, f);
      x = Fp_mul(x, subiu(gel(vz_f, y), 1), ell); /* vz_f[y] = z_f^y  */
    }
    gi = Fl_mul(gi, g_K, f);
    gel(xi_el, i) = x;  /* xi_el[i]=xi^{g_K^i} */
  }
  for (i=0; i<d_K; i++)
  {
    GEN x = gen_1;
    for (j=0; j<d_K; j++)
      x = Fp_mul(x, Fp_powu(gel(xi_el, j), e_chi[(i+j)%d_K], ell), ell);
    x = Fp_pow(x, ell1d2, ell);
    if (!equali1(x)) return gc_long(av, 0);
  }
  return gc_long(av, 1);
}

static GEN
next_el_real(GEN K, long p, long d_pow, GEN elg, long j0, long flag)
{
  GEN Chi = gel(K, 2);
  ulong f = K_get_f(K), h = K_get_h(K), d = upowuu(p, d_pow), d2 = d*d;
  ulong D = (flag & USE_F)? d2*f: d2<<1, el = elg[1] + D;

  /* O(el*h) -> O(el*log(el)) by FFT */
  if (1000 < h && el < h) { el = (h/D)*D+1; if (el < h) el += D; }
  if (flag&USE_F)  /* el=1 (mod f) */
  {
    for (;; el += D)
      if (uisprime(el) && chk_el_real_f(K, p, d_pow, el, j0)) break;
  }
  else
  {
    for (;; el += D)
      if (Chi[el%f]==0 && uisprime(el) &&
          chk_el_real_chi(K, p, d_pow, el, j0, flag)) break;
  }
  return mkvecsmall2(el, pgener_Fl(el));
}

static GEN
next_ell_real(GEN K, GEN ellg, long d2, GEN df0l, long j0)
{
  GEN ell = gel(ellg, 1);
  for (ell = addii(ell, df0l);; ell = addii(ell, df0l))
    if (BPSW_psp(ell) && chk_ell_real(K, d2, ell, j0))
      return mkvec2(ell, pgener_Fp(ell));
}

/* #velg >= n */
static long
delete_el(GEN velg, long n)
{
  long i, l;
  if (DEBUGLEVEL>1) err_printf("deleting %ld ...\n", gmael(velg, n, 1));
  for (l = lg(velg)-1; l >= 1; l--) if (gel(velg, l)) break;
  for (i = n; i < l; i++) gel(velg, i) = gel(velg, i+1);
  return l;
}

/* velg has n components */
static GEN
set_ell_real(GEN K, GEN velg, long n, long d_chi, long d2, long f0, long j0)
{
  long i, n_ell = n*d_chi;
  GEN z = cgetg(n_ell + 1, t_VEC);
  GEN df0l = muluu(d2, f0), ellg = mkvec2(gen_1, gen_1);
  for (i=1; i<=n; i++) df0l = muliu(df0l, gel(velg, i)[1]);
  for (i=1; i<=n_ell; i++) ellg = gel(z, i)= next_ell_real(K, ellg, d2, df0l, j0);
  return z;
}

static GEN
G_K_vell(GEN K, GEN vellg, ulong gk)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong f = K_get_f(K), h = K_get_h(K);
  GEN z_f, vz_f, A, P, M, z =  cgetg(h+1, t_VEC);
  ulong i, lv = lg(vellg);

  A=cgetg(lv, t_VEC);
  P=cgetg(lv, t_VEC);
  for (i=1; i<lv; i++)
  {
    GEN ell = gmael(vellg, i, 1), g_ell = gmael(vellg, i, 2);
    gel(A, i) = Fp_pow(g_ell, diviuexact(subiu(ell, 1), f), ell);
    gel(P, i) = ell;
  }
  z_f = ZV_chinese(A, P, &M);
  vz_f = Fp_powers(z_f, f-1, M)+1;
  for (i=1; i<=h; i++) gel(z, i) = gel(vz_f, Fl_mul(H[i], gk, f));
  return gerepilecopy(av, FpV_roots_to_pol(z, M, 0));
}

/* f=cond(K), M=product of ell in vell, G(K/Q)=<g_K>
 * G_K[1+i]=minimal polynomial of zeta_f^{g_k^i} over K mod M, 0 <= i < d_K */
static GEN
make_G_K(GEN K, GEN vellg)
{
  ulong d_K = K_get_d(K), f = K_get_f(K), g_K = K_get_g(K);
  GEN G_K = cgetg(d_K+1, t_VEC);
  ulong i, g = 1;

  for (i=0; i<d_K; i++)
  {
    gel(G_K, 1+i) = G_K_vell(K, vellg, g);
    g = Fl_mul(g, g_K, f);
  }
  return G_K;
}

static GEN
G_K_p(GEN K, GEN ellg, ulong gk)
{
  pari_sp av = avma;
  ulong i, f = K_get_f(K), h = K_get_h(K);
  GEN ell = gel(ellg, 1), g_ell = gel(ellg, 2);
  GEN H = K_get_H(K), z_f, vz_f, z = cgetg(h+1, t_VEC);

  z_f = Fp_pow(g_ell, diviuexact(subiu(ell, 1), f), ell);
  vz_f = Fp_powers(z_f, f-1, ell)+1;
  for (i=1; i<=h; i++) gel(z, i) = gel(vz_f, Fl_mul(H[i], gk, f));
  return gerepilecopy(av, FpV_roots_to_pol(z, ell, 0));
}

static GEN
G_K_l(GEN K, GEN ellg, ulong gk)
{
  pari_sp av = avma;
  ulong ell = itou(gel(ellg, 1)), g_ell = itou(gel(ellg, 2));
  ulong f = K_get_f(K), h = K_get_h(K), i, z_f;
  GEN H = K_get_H(K), vz_f, z = cgetg(h+1, t_VEC);

  z_f = Fl_powu(g_ell, (ell-1) / f, ell);
  vz_f = Fl_powers(z_f, f-1, ell)+1;
  for (i=1; i<=h; i++) z[i] = vz_f[Fl_mul(H[i], gk, f)];
  return gerepilecopy(av, Flv_roots_to_pol(z, ell, 0));
}

static GEN
vz_vell(long d, GEN vellg, GEN *pM)
{
  long i, l = lg(vellg);
  GEN A = cgetg(l, t_VEC), P = cgetg(l, t_VEC), z;

  for (i = 1; i < l; i++)
  {
    GEN ell = gmael(vellg, i, 1), g_ell = gmael(vellg, i, 2);
    gel(A, i) = Fp_pow(g_ell, diviuexact(subiu(ell, 1), d), ell);
    gel(P, i) = ell;
  }
  z = ZV_chinese(A, P, pM); return Fp_powers(z, d-1, *pM);
}

static GEN
D_xi_el_vell_FFT(GEN K, GEN elg, GEN vellg, ulong d, ulong j0, GEN vG_K)
{
  pari_sp av = avma;
  ulong d_K, h = K_get_h(K), d_chi = K_get_dchi(K);
  ulong el = elg[1], g_el = elg[2], el_1 = el-1;
  ulong i, j, i2, k, dwel;
  GEN u = cgetg(el+2, t_POL) , v = cgetg(h+3, t_POL);
  GEN w = cgetg(el+1, t_VEC), ww;
  GEN M, vz_el, G_K, z = const_vec(d_chi, gen_1);
  GEN e_chi = get_e_chi(K, j0, d, &d_K);

  vz_el = vz_vell(el, vellg, &M);
  u[1] = evalsigne(1) | evalvarn(0);
  v[1] = evalsigne(1) | evalvarn(0);

  for (i=i2=0; i<el; i++)
  {
    ulong j2 = i2?el-i2:i2; /* i2=(i*i)%el */
    gel(u, 2+i) = gel(vz_el, 1+j2);
    if ((i2+=i+i+1)>=el) i2%=el;
  }
  for (k=0; k<d_K; k++)
  {
    pari_sp av = avma;
    pari_timer ti;
    long gd, gi;
    GEN x1 = gen_1;
    G_K = gel(vG_K, 1+k);
    for (i=i2=0; i<=h; i++)
    {
      gel(v, 2+i) = Fp_mul(gel(G_K, 2+i), gel(vz_el, 1+i2), M);
      if ((i2+=i+i+1)>=el) i2%=el;
    }
    if (DEBUGLEVEL>2) timer_start(&ti);
    ww = ZX_mul(u, v);
    if (DEBUGLEVEL>2)
      timer_printf(&ti, "ZX_mul:%ld/%ld h*el=%ld*%ld", k, d_K, h, el);
    dwel = degpol(ww)-el;
    for (i=0; i<=dwel; i++) gel(w, 1+i) = addii(gel(ww, 2+i), gel(ww, 2+i+el));
    for (; i<el; i++) gel(w, 1+i) = gel(ww, 2+i);
    for (i=i2=1; i<el; i++)  /* w[i]=G_K(z_el^(2*i)) */
    {
      gel(w, i) = Fp_mul(gel(w, 1+i), gel(vz_el, 1+i2), M);
      if ((i2+=i+i+1)>=el) i2%=el;
    }
    gd = Fl_powu(g_el, d, el);  /* a bit faster */
    gi = g_el;
    for (i=1; i<d; i++)
    {
      GEN xi = gen_1;
      long gdi = gi;
      for (j=0; i+j<el_1; j+=d)
      {
        xi = Fp_mul(xi, gel(w, (gdi+gdi)%el), M);
        gdi = Fl_mul(gdi, gd, el);
      }
      gi = Fl_mul(gi, g_el, el);
      xi = Fp_powu(xi, i, M);
      x1 = Fp_mul(x1, xi, M);
    }
    for (i=1; i<=d_chi; i++)
    {
      GEN x2 = Fp_powu(x1, e_chi[(k+i-1)%d_K], M);
      gel(z, i) = Fp_mul(gel(z, i), x2, M);
    }
    z = gerepilecopy(av, z);
  }
  return gerepilecopy(av, z);
}

static GEN
D_xi_el_vell(GEN K, GEN elg, GEN vellg, ulong d, ulong j0)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  GEN z_f, z_el, vz_f, vz_el;
  ulong el = elg[1], g_el = elg[2], el_1 = el-1;
  ulong i, j, k, d_K, lv = lg(vellg), d_chi = K_get_dchi(K);
  GEN A, B, P, M, z = const_vec(d_chi, gen_1);
  GEN e_chi = get_e_chi(K, j0, d, &d_K);

  A=cgetg(lv, t_VEC);
  B=cgetg(lv, t_VEC);
  P=cgetg(lv, t_VEC);
  for (i = 1; i < lv; i++)
  {
    GEN ell = gmael(vellg, i, 1), g_ell = gmael(vellg, i, 2);
    GEN ell_1 = subiu(ell, 1);
    gel(A, i) = Fp_pow(g_ell, diviuexact(ell_1, f), ell);
    gel(B, i) = Fp_pow(g_ell, diviuexact(ell_1, el), ell);
    gel(P, i) = ell;
  }
  z_f = ZV_chinese(A, P, &M);
  z_el = ZV_chinese(B, P, NULL);
  vz_f = Fp_powers(z_f, f-1, M);
  vz_el = Fp_powers(z_el, el-1, M);
  for (k = 0; k < d_K; k++)
  {
    pari_sp av = avma;
    GEN x0 = gen_1;
    long gk = Fl_powu(g_K, k, f);
    for (i=1; i<el_1; i++)
    {
      long gi = Fl_powu(g_el, i, el);
      GEN x1 = gen_1;
      GEN x2 = gel(vz_el, 1+gi);
      for (j=1; j<=h; j++)
      {
        long y = Fl_mul(H[j], gk, f);
        x1 = Fp_mul(x1, Fp_sub(x2, gel(vz_f, 1+y), M), M);
      }
      x1 = Fp_powu(x1, i%d, M);
      x0 = Fp_mul(x0, x1, M);
    }
    for (i=1; i<=d_chi; i++)
    {
      GEN x2 = Fp_powu(x0, e_chi[(k+i-1)%d_K], M);
      gel(z, i) = Fp_mul(gel(z, i), x2, M);
    }
    z = gerepilecopy(av, z);
  }
  return gerepilecopy(av, z);
}

static GEN
D_xi_el_Flx_mul(GEN K, GEN elg, GEN ellg, GEN vG_K, ulong d, ulong j0)
{
  pari_sp av = avma;
  ulong d_K, f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  ulong el = elg[1], g_el = elg[2], el_1 = el-1, d_chi = K_get_dchi(K);
  ulong ell = itou(gel(ellg, 1)), g_ell = itou(gel(ellg, 2)), z_el;
  GEN u = cgetg(el+2, t_VECSMALL), v = cgetg(h+3, t_VECSMALL);
  GEN w = cgetg(el+1, t_VECSMALL), ww;
  GEN vz_el, G_K, z = const_vecsmall(d_chi, 1);
  GEN e_chi = get_e_chi(K, j0, d, &d_K);
  ulong i, j, i2, k, dwel;

  u[1] = evalvarn(0);
  v[1] = evalvarn(0);
  z_el = Fl_powu(g_ell, (ell - 1) / el, ell);
  vz_el = Fl_powers(z_el, el_1, ell)+1;

  for (i=i2=0; i<el; i++)
  {
    ulong j2 = i2?el-i2:i2;
    u[2+i] = vz_el[j2];
    if ((i2+=i+i+1)>=el) i2%=el;  /* i2=(i*i)%el */
  }
  for (k=0; k<d_K; k++)
  {
    pari_sp av = avma;
    pari_timer ti;
    ulong gk = Fl_powu(g_K, k, f);
    long gd, gi, x1 = 1;
    if (DEBUGLEVEL>2) timer_start(&ti);
    G_K = (vG_K==NULL)?G_K_l(K, ellg, gk):ZX_to_Flx(gel(vG_K, 1+k), ell);
    if (DEBUGLEVEL>2) timer_printf(&ti, "G_K_l");
    for (i=i2=0; i<=h; i++)
    {
      v[2+i] = Fl_mul(G_K[2+i], vz_el[i2], ell);
      if ((i2+=i+i+1)>=el) i2%=el;  /* i2=(i*i)%el */
    }
    if (DEBUGLEVEL>2) timer_start(&ti);
    ww = Flx_mul(u, v, ell);
    if (DEBUGLEVEL>2)
      timer_printf(&ti, "Flx_mul:%ld/%ld h*el=%ld*%ld", k, d_K, h, el);
    dwel=degpol(ww)-el; /* dwel=h-1 */
    for (i=0; i<=dwel; i++) w[1+i] = Fl_add(ww[2+i], ww[2+i+el], ell);
    for (; i<el; i++) w[1+i] = ww[2+i];
    for (i=i2=1; i<el; i++)  /* w[i]=G_K(z_el^(2*i)) */
    {
      w[i] = Fl_mul(w[1+i], vz_el[i2], ell);
      if ((i2+=i+i+1)>=el) i2%=el;  /* i2=(i*i)%el */
    }
    gd = Fl_powu(g_el, d, el);  /* a bit faster */
    gi = g_el;
    for (i=1; i<d; i++)
    {
      long xi = 1, gdi = gi;
      for (j=0; i+j<el_1; j+=d)
      {
        xi = Fl_mul(xi, w[(gdi+gdi)%el], ell);
        gdi = Fl_mul(gdi, gd, el);
      }
      gi = Fl_mul(gi, g_el, el);
      xi = Fl_powu(xi, i, ell);
      x1 = Fl_mul(x1, xi, ell);
    }
    for (i=1; i<=d_chi; i++)
    {
      long x2 = Fl_powu(x1, e_chi[(k+i-1)%d_K], ell);
      z[i] = Fl_mul(z[i], x2, ell);
    }
    set_avma(av);
  }
  return gerepilecopy(av, Flv_to_ZV(z));
}

static GEN
D_xi_el_ZX_mul(GEN K, GEN elg, GEN ellg, GEN vG_K, ulong d, ulong j0)
{
  pari_sp av = avma;
  GEN ell = gel(ellg,1), g_ell, u, v, w, ww, z_el, vz_el, G_K, z, e_chi;
  ulong d_K, f, h, g_K, el, g_el, el_1, d_chi, i, j, i2, k, dwel;

  if (lgefint(ell) == 3) return D_xi_el_Flx_mul(K, elg, ellg, vG_K, d, j0);
  f = K_get_f(K); h = K_get_h(K); g_K = K_get_g(K);
  el = elg[1]; g_el = elg[2]; el_1 = el-1; d_chi = K_get_dchi(K);
  g_ell = gel(ellg, 2);
  z = const_vec(d_chi, gen_1);
  e_chi = get_e_chi(K, j0, d, &d_K);

  u = cgetg(el+2,t_POL); u[1] = evalsigne(1) | evalvarn(0);
  v = cgetg(h+3, t_POL); v[1] = evalsigne(1) | evalvarn(0);
  w = cgetg(el+1, t_VEC);
  z_el = Fp_pow(g_ell, diviuexact(subiu(ell, 1), el), ell);
  vz_el = Fp_powers(z_el, el_1, ell)+1;

  for (i=i2=0; i<el; i++)
  {
    ulong j2 = i2?el-i2:i2; /* i2=(i*i)%el */
    gel(u, 2+i) = gel(vz_el, j2);
    if ((i2+=i+i+1)>=el) i2%=el;
  }
  for (k=0; k<d_K; k++)
  {
    pari_sp av = avma;
    pari_timer ti;
    long gd, gi, gk = Fl_powu(g_K, k, f);
    GEN x1 = gen_1;
    if (DEBUGLEVEL>2) timer_start(&ti);
    G_K = (vG_K==NULL) ? G_K_p(K, ellg, gk):RgX_to_FpX(gel(vG_K, 1+k), ell);
    if (DEBUGLEVEL>2) timer_printf(&ti, "G_K_p");
    for (i=i2=0; i<=h; i++)
    {
      gel(v, 2+i) = Fp_mul(gel(G_K, 2+i), gel(vz_el, i2), ell);
      if ((i2+=i+i+1)>=el) i2%=el;
    }
    if (DEBUGLEVEL>2) timer_start(&ti);
    ww = ZX_mul(u, v);
    if (DEBUGLEVEL>2)
      timer_printf(&ti, "ZX_mul:%ld/%ld h*el=%ld*%ld", k, d_K, h, el);
    dwel = degpol(ww)-el;
    for (i=0; i<=dwel; i++) gel(w, 1+i) = addii(gel(ww, 2+i), gel(ww, 2+i+el));
    for (; i<el; i++) gel(w, 1+i) = gel(ww, 2+i);
    for (i=i2=1; i<el; i++)  /* w[i]=G_K(z_el^(2*i)) */
    {
      gel(w, i) = Fp_mul(gel(w, 1+i), gel(vz_el, i2), ell);
      if ((i2+=i+i+1)>=el) i2%=el;
    }
    gd = Fl_powu(g_el, d, el);  /* a bit faster */
    gi = g_el;
    for (i=1; i<d; i++)
    {
      GEN xi = gen_1;
      long gdi = gi;
      for (j=0; i+j<el_1; j+=d)
      {
        xi = Fp_mul(xi, gel(w, (gdi+gdi)%el), ell);
        gdi = Fl_mul(gdi, gd, el);
      }
      gi = Fl_mul(gi, g_el, el);
      xi = Fp_powu(xi, i, ell);
      x1 = Fp_mul(x1, xi, ell);
    }
    for (i=1; i<=d_chi; i++)
    {
      GEN x2 = Fp_powu(x1, e_chi[(k+i-1)%d_K], ell);
      gel(z, i) = Fp_mul(gel(z, i), x2, ell);
    }
    z = gerepilecopy(av, z);
  }
  return gerepilecopy(av, z);
}

static GEN
D_xi_el_ss(GEN K, GEN elg, GEN ellg, ulong d, ulong j0)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong d_K, f = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  ulong el = elg[1], g_el = elg[2], el_1 = el-1;
  ulong ell = itou(gel(ellg, 1)), g_ell = itou(gel(ellg, 2));
  ulong i, j, k, gk, z_f, z_el, d_chi = K_get_dchi(K);
  GEN vz_f, vz_el, z = const_vecsmall(d_chi, 1);
  GEN e_chi = get_e_chi(K, j0, d, &d_K);

  z_f = Fl_powu(g_ell, (ell - 1) / f, ell);
  z_el = Fl_powu(g_ell, (ell - 1) / el, ell);
  vz_f = Fl_powers(z_f, f-1, ell)+1;
  vz_el = Fl_powers(z_el, el-1, ell)+1;
  gk = 1; /* g_K^k */
  for (k = 0; k < d_K; k++)
  {
    ulong x0 = 1, gi = g_el; /* g_el^i */
    for (i = 1; i < el_1; i++)
    {
      ulong x1 = 1, x2 = vz_el[gi];
      for (j=1; j<=h; j++)
      {
        ulong y = Fl_mul(H[j], gk, f);
        x1 = Fl_mul(x1, Fl_sub(x2, vz_f[y], ell), ell);
      }
      x1 = Fl_powu(x1, i%d, ell);
      x0 = Fl_mul(x0, x1, ell);
      gi = Fl_mul(gi, g_el, el);
    }
    for (i = 1; i <= d_chi; i++)
    {
      ulong x2 = Fl_powu(x0, e_chi[(k+i-1)%d_K], ell);
      z[i] = Fl_mul(z[i], x2, ell);
    }
    gk = Fl_mul(gk, g_K, f);
  }
  return gerepileupto(av, Flv_to_ZV(z));
}

static GEN
D_xi_el_sl(GEN K, GEN elg, GEN ellg, ulong d, ulong j0)
{
  pari_sp av = avma;
  GEN ell = gel(ellg, 1), H;
  GEN g_ell, ell_1, z_f, z_el, vz_f, vz_el, z, e_chi;
  ulong d_K, f, h, g_K, el, g_el, el_1, d_chi, i, j, k, gk;

  if (lgefint(ell) == 3) return D_xi_el_ss(K, elg, ellg, d, j0);
  H = K_get_H(K);
  f = K_get_f(K); h = K_get_h(K); g_K = K_get_g(K);
  el = elg[1]; g_el = elg[2]; el_1 = el-1; d_chi = K_get_dchi(K);
  g_ell = gel(ellg, 2); ell_1 = subiu(ell, 1);
  z = const_vec(d_chi, gen_1);
  e_chi = get_e_chi(K, j0, d, &d_K);

  z_f = Fp_pow(g_ell, diviuexact(ell_1, f), ell);
  z_el = Fp_pow(g_ell, diviuexact(ell_1, el), ell);
  vz_f = Fp_powers(z_f, f-1, ell) + 1;
  vz_el = Fp_powers(z_el, el-1, ell) + 1;
  gk = 1; /* g_K^k */
  for (k = 0; k < d_K; k++)
  {
    pari_sp av2 = avma;
    GEN x0 = gen_1;
    ulong gi = g_el; /* g_el^i */
    for (i = 1; i < el_1; i++)
    {
      pari_sp av3 = avma;
      GEN x1 = gen_1, x2 = gel(vz_el, gi);
      for (j = 1; j <= h; j++)
      {
        ulong y = Fl_neg(Fl_mul(H[j], gk, f), f);
        x1 = Fp_mul(x1, Fp_sub(x2, gel(vz_f, y), ell), ell);
      }
      x1 = Fp_powu(x1, i%d, ell);
      x0 = gerepileuptoint(av3, Fp_mul(x0, x1, ell));
      gi = Fl_mul(gi, g_el, el);
    }
    for (i=1; i<=d_chi; i++)
    {
      GEN x2 = Fp_powu(x0, e_chi[(k+i-1)%d_K], ell);
      gel(z, i) = Fp_mul(gel(z, i), x2, ell);
    }
    if (k == d_K-1) break;
    z = gerepilecopy(av2, z);
    gk = Fl_mul(gk, g_K, f);
  }
  return gerepilecopy(av, z);
}

static long
get_y(GEN z, GEN ellg, long d)
{
  GEN ell = gel(ellg, 1), g_ell = gel(ellg, 2);
  GEN elld = diviuexact(subiu(ell, 1), d);
  GEN g_elld = Fp_pow(g_ell, elld, ell);
  GEN x = Fp_pow(modii(z, ell), elld, ell);
  long k;
  for (k=0; k<d; k++)
  {
    if (equali1(x)) break;
    x = Fp_mul(x, g_elld, ell);
  }
  if (k==0) k=d;
  else if (d<=k) pari_err_BUG("subcyclopclgp [MLL]");
  return k;
}

static void
real_MLLn(long *y, GEN K, ulong p, ulong d_pow, ulong n,
    GEN velg, GEN vellg, GEN vG_K, ulong j0)
{
  pari_sp av = avma;
  ulong i, j, k, d = upowuu(p, d_pow), h = gmael(K, 1, 2)[3];
  ulong row = lg(vellg)-1;
  for (i=1; i<=n; i++)
  {
    GEN elg = gel(velg, i), z;
    ulong el = elg[1], nz;
    pari_timer ti;
    if (DEBUGLEVEL>1) timer_start(&ti);
    z = (h<el) ? D_xi_el_vell_FFT(K, elg, vellg, d, j0, vG_K)
               : D_xi_el_vell(K, elg, vellg, d, j0);
    if (DEBUGLEVEL>1) timer_printf(&ti, "subcyclopclgp:[D_xi_el]");
    if (DEBUGLEVEL>2) err_printf("z=%Ps\n", z);
    nz = lg(z)-1;
    for (k = 1; k <= nz; k++)
      for (j=1; j<=row; j++)
        y[(j-1)*row+(i-1)*nz+k-1] = get_y(gel(z, k), gel(vellg, j), d);
    set_avma(av);
  }
}

static void
real_MLL1(long *y, GEN K, ulong p, ulong d_pow, GEN velg, GEN vellg, ulong j0)
{
  ulong h = gmael(K, 1, 2)[3], d = upowuu(p, d_pow);
  GEN elg = gel(velg, 1), ellg = gel(vellg, 1), z;
  ulong el = elg[1];
  pari_timer ti;

  if (DEBUGLEVEL>2) timer_start(&ti);
  z = h < el? D_xi_el_ZX_mul(K, elg, ellg, NULL, d, j0)
            : D_xi_el_sl(K, elg, ellg, d, j0);
  if (DEBUGLEVEL>2) timer_printf(&ti, "subcyclopclgp:[D_xi_el]");
  if (DEBUGLEVEL>2) err_printf("z=%Ps\n", z);
  y[0] = get_y(gel(z, 1), ellg, d);
}

static void
real_MLL(long *y, GEN K, ulong p, ulong d_pow, ulong n,
    GEN velg, GEN vellg, GEN vG_K, ulong j0)
{
  ulong i, j, k, d = upowuu(p, d_pow), h = gmael(K, 1, 2)[3];
  ulong row = lg(vellg)-1;
  for (j=1; j<=row; j++)
  {
    GEN ellg = gel(vellg, j);
    for (i=1; i<=n; i++)
    {
      pari_sp av2 = avma;
      GEN elg = gel(velg, i), z;
      ulong el = elg[1], nz;
      pari_timer ti;
      if (DEBUGLEVEL>2) timer_start(&ti);
      z = h < el? D_xi_el_ZX_mul(K, elg, ellg, vG_K, d, j0)
                : D_xi_el_sl(K, elg, ellg, d, j0);
      if (DEBUGLEVEL>2) timer_printf(&ti, "subcyclopclgp:[D_xi_el]");
      if (DEBUGLEVEL>3) err_printf("z=%Ps\n", z);
      nz = lg(z)-1;
      for (k = 1; k <= nz; k++)
        y[(j-1)*row+(i-1)*nz+k-1] = get_y(gel(z, k), ellg, d);
      set_avma(av2);
    }
  }
}

static long
use_basis(long d_K, long f) { return (d_K<=10 || (d_K<=30 && f<=5000)); }

static long
use_factor(ulong f)
{ GEN fa = factoru(f), P = gel(fa, 1); return (P[lg(P)-1]<500); }

/* group structure, destroy gr */
static GEN
get_str(GEN gr)
{
  GEN z = gel(gr,2);
  long i, j, l = lg(z);
  for (i = j = 1; i < l; i++)
    if (lgefint(gel(z, i)) > 2) gel(z,j++) = gel(z,i);
  setlg(z, j); return z;
}

static GEN
cyc_real_MLL(GEN K, ulong p, long d_pow, long j0, long flag)
{
  ulong d_K = K_get_d(K), f = K_get_f(K), d_chi = K_get_dchi(K);
  ulong n, n0 = 1, f0, n_el = d_pow, d = upowuu(p, d_pow), rank = n_el*d_chi;
  GEN velg = const_vec(n_el, NULL), vellg = NULL;
  GEN oldgr = mkvec2(gen_0, NULL), newgr = mkvec2(gen_0, NULL);
  long *y0 = (long*)stack_calloc(sizeof(long)*rank*rank);

  if (DEBUGLEVEL>1)
    err_printf("cyc_real_MLL:p=%ld d_pow=%ld deg(K)=%ld cond(K)=%ld g_K=%ld\n",
        p, d_pow, d_K, f, K_get_g(K));
  gel(K, 2) = get_chi(gel(K,1));
  if (f-1 <= (d_K<<1)) flag |= USE_F;
  else if (use_basis(d_K, f)) flag |= USE_BASIS;
  else if (use_factor(f)) flag |= USE_FACTOR;
  else flag |= USE_GALOIS_POL;
  if (flag&USE_BASIS) K = vec_append(K, xi_data_basis(K));
  else if (flag&USE_GALOIS_POL) K = vec_append(K, xi_data_galois(K));
  f0 = f%p?f:f/p;
  gel(velg, 1) = next_el_real(K, p, d_pow, mkvecsmall2(1, 1), j0, flag);
  if (flag&USE_FULL_EL)
  {
    for (n=2; n<=n_el; n++)
      gel(velg, n) = next_el_real(K, p, d_pow, gel(velg, n+1), j0, flag);
    n0 = n_el;
  }

  for (n=n0; n<=n_el; n++) /* loop while structure is unknown */
  {
    pari_sp av2 = avma;
    long n_ell, m, M;
    GEN y;
    pari_timer ti;
    if (DEBUGLEVEL>2) timer_start(&ti);
    vellg = set_ell_real(K, velg, n, d_chi, d*d, f0, j0);
    n_ell = lg(vellg) -1; /* equal to n*d_chi */
    if (DEBUGLEVEL>2) timer_printf(&ti, "set_ell_real");
    if (DEBUGLEVEL>3) err_printf("vel=%Ps\nvell=%Ps\n", velg, vellg);
    if (n_ell==1)
      real_MLL1(y0, K, p, d_pow, velg, vellg, j0);
    else
    {
      GEN vG_K;
      if (DEBUGLEVEL>2) timer_start(&ti);
      vG_K = make_G_K(K, vellg);
      if (DEBUGLEVEL>2) timer_printf(&ti, "make_G_K");
      if (lgefint(gmael(vellg, n_ell, 1))<=3 || (flag&SAVE_MEMORY))
        real_MLL(y0, K, p, d_pow, n, velg, vellg, vG_K, j0);
      else
        real_MLLn(y0, K, p, d_pow, n, velg, vellg, vG_K, j0);
    }
    set_avma(av2);
    y = ary2mat(y0, n_ell);
    if (DEBUGLEVEL>3) err_printf("y=%Ps\n", y);
    y = ZM_snf(y);
    if (DEBUGLEVEL>3) err_printf("y=%Ps\n", y);
    y = make_p_part(y, p, d_pow);
    if (DEBUGLEVEL>3) err_printf("y=%Ps\n", y);
    newgr = structure_MLL(y, d_pow);
    if (DEBUGLEVEL>3)
      err_printf("d_pow=%ld d_chi=%ld old=%Ps new=%Ps\n",d_pow,d_chi,oldgr,newgr);
    if (equalsi(d_pow*d_chi, gel(newgr, 1))) break;
    if ((m = find_del_el(&oldgr, newgr, n, n_el, d_chi)))
    { M = m = delete_el(velg, m); n--; }
    else
    { M = n+1; m = n; }
    gel(velg, M) = next_el_real(K, p, d_pow, gel(velg, m), j0, flag);
  }
  return get_str(newgr);
}

static GEN
cyc_buch(long dK, GEN p, long d_pow)
{
  GEN z = Buchquad(stoi(dK), 0.0, 0.0, 0), cyc = gel(z,2);
  long i, l = lg(cyc);
  if (Z_pval(gel(z,1), p) != d_pow) pari_err_BUG("subcyclopclgp [Buchquad]");
  for (i = 1; i < l; i++)
  {
    long x = Z_pval(gel(cyc, i), p); if (!x) break;
    gel(cyc, i) = utoipos(x);
  }
  setlg(cyc, i); return cyc;
}

static void
verbose_output(GEN K, GEN p, long pow, long j)
{
  long d = K_get_d(K), f = K_get_f(K), s = K_get_s(K), d_chi = K_get_dchi(K);
  err_printf("|A_K_psi|=%Ps^%ld, psi=chi^%ld, d_psi=%ld, %s,\n\
    [K:Q]=%ld, [f,H]=[%ld, %Ps]\n",
    p,pow*d_chi,j,d_chi,s?"real":"imaginary",d,f,zv_to_ZV(gmael3(K,1,1,1)));
}

static int
cyc_real_pre(GEN K, GEN xi, ulong p, ulong j, long el)
{
  pari_sp av = avma;
  ulong i, d_K, x = 1;
  GEN e_chi = get_e_chi(K, j, p, &d_K);

  xi++;
  for (i = 0; i < d_K; i++) x = Fl_mul(x, Fl_powu(xi[i], e_chi[i], el), el);
  return gc_ulong(av, Fl_powu(x, (el-1)/p, el));
}

/* return vec[-1,[],0], vec[0,[],0], vec[1,[1],0], vec[2,[1,1],0] etc */
static GEN
cyc_real_ss(GEN K, GEN xi, ulong p, long j, long pow, long el, ulong pn, long flag)
{
  ulong d_chi = K_get_dchi(K);
  if (cyc_real_pre(K, xi, pn, j, el) == 1) return NULL; /* not determined */
  if (--pow==0) return mkvec3(gen_0, nullvec(), gen_0); /* trivial */
  if (DEBUGLEVEL) verbose_output(K, utoi(p), pow, j);
  if (flag&USE_MLL)
  {
    pari_sp av = avma;
    GEN gr = (K_get_d(K) == 2)? cyc_buch(K_get_f(K), utoi(p), pow)
                               : cyc_real_MLL(K, p, pow, j, flag);
    return gerepilecopy(av, mkvec3(utoipos(d_chi * pow), gr, gen_0));
  }
  if (pow==1) return mkvec3(utoi(d_chi), onevec(d_chi), gen_0);
  return mkvec3(utoi(pow*d_chi), nullvec(), gen_0);
}

static GEN
cyc_real_ll(GEN K, GEN xi, GEN p, long j, long pow, GEN el, GEN pn, long flag)
{
  pari_sp av = avma;
  ulong i, d_K = K_get_d(K), d_chi = K_get_dchi(K);
  GEN e_chi = get_e_chi_ll(K, j, pn), x = gen_1;

  xi++;
  for (i = 0; i < d_K; i++)
    x = Fp_mul(x, Fp_pow(gel(xi, i), gel(e_chi, i), el), el);
  x = Fp_pow(x, diviiexact(subiu(el, 1), pn), el); /* x = x^(el-1)/pn mod el */
  set_avma(av); if (equali1(x)) return NULL; /* not determined */
  if (--pow==0) return mkvec3(gen_0, nullvec(), gen_0); /* trivial */
  if (DEBUGLEVEL) verbose_output(K, p, pow, j);
  if (flag&USE_MLL)
    pari_err_IMPL(stack_sprintf("flag=%ld for large prime", USE_MLL));
  if (pow==1) return mkvec3(utoi(d_chi), onevec(d_chi), gen_0);
  return mkvec3(utoi(pow*d_chi), nullvec(), gen_0);
}

/* xi[1+i] = xi^(g^i), 0 <= i <= d-1 */
static GEN
xi_conj_s(GEN K, ulong el)
{
  pari_sp av = avma;
  GEN  H = K_get_H(K);
  long d = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g = K_get_g(K);
  long i, gi = 1, z = Fl_powu(pgener_Fl(el), (el-1)/f, el);
  GEN vz = Fl_powers(z, f, el)+1, xi = cgetg(d+1, t_VECSMALL);

  for (i=1; i<=d; i++)
  {
    long j, x = 1;
    for (j=1; j<=h; j++)
      x = Fl_mul(x, vz[Fl_mul(H[j], gi, f)]-1, el);
    xi[i] = x;
    gi = Fl_mul(gi, g, f);
  }
  return gerepilecopy(av, xi);
}

static GEN
xi_conj_l(GEN K, GEN el)
{
  pari_sp av = avma;
  GEN  H = K_get_H(K);
  long d = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g = K_get_g(K);
  long i, gi = 1;
  GEN z = Fp_pow(pgener_Fp(el), diviuexact(subiu(el, 1), f), el);
  GEN vz = Fp_powers(z, f, el)+1, xi = cgetg(d+1, t_VEC);

  for (i=1; i<=d; i++)
  {
    long j;
    GEN x = gen_1;
    for (j=1; j<=h; j++)
      x = Fp_mul(x, subiu(gel(vz, Fl_mul(H[j], gi, f)), 1), el);
    gel(xi, i) = x;
    gi = Fl_mul(gi, g, f);
  }
  return gerepilecopy(av, xi);
}

static GEN
pclgp_cyc_real(GEN K, GEN p, long max_pow, long flag)
{
  const long NUM_EL = 20;
  GEN C = gel(K, 5);
  long f_K = K_get_f(K), n_conj = K_get_nconj(K);
  long i, pow, n_el, n_done = 0;
  GEN gr = nullvec(), Done = const_vecsmall(n_conj, 0), xi;
  long first = 1;

  for (pow=1; pow<=max_pow; pow++)
  {
    GEN pn = powiu(p, pow), fpn = muliu(pn, f_K), el = addiu(fpn, 1);
    for (n_el = 0; n_el < NUM_EL; el = addii(el, fpn))
    {
      ulong uel;
      if (!BPSW_psp(el)) continue;
      n_el++; uel = itou_or_0(el);
      if (uel)
      {
        xi = xi_conj_s(K, uel);
        if (first && n_conj > 10) /* mark trivial chi-part */
        {
          for (i = 1; i <= n_conj; i++)
          {
            if (cyc_real_pre(K, xi, p[2], C[i], uel) == 1) continue;
            Done[i] = 1;
            if (++n_done == n_conj) return gr;
          }
          first = 0; continue;
        }
      }
      else
        xi = xi_conj_l(K, el);
      for (i = 1; i <= n_conj; i++)
      {
        GEN z;
        if (Done[i]) continue;
        if (uel)
          z = cyc_real_ss(K, xi, p[2], C[i], pow, uel, itou(pn), flag);
        else
          z = cyc_real_ll(K, xi, p, C[i], pow, el, pn, flag);
        if (!z) continue;
        Done[i] = 1;
        if (!isintzero(gel(z, 1))) gr = vec_append(gr, z);
        if (++n_done == n_conj) return gr;
      }
    }
  }
  pari_err_BUG("pclgp_cyc_real: max_pow is not enough");
  return NULL; /*LCOV_EXCL_LINE*/
}

/* return (el, g_el) */
static GEN
next_el_imag(GEN elg, long f)
{
  long el = elg[1];
  if (odd(f)) f<<=1;
  while (!uisprime(el+=f));
  return mkvecsmall2(el, pgener_Fl(el));
}

/* return (ell, g_ell) */
static GEN
next_ell_imag(GEN ellg, GEN df0l)
{
  GEN ell = gel(ellg, 1);
  while (!BPSW_psp(ell = addii(ell, df0l)));
  return mkvec2(ell, pgener_Fp(ell));
}

static GEN
set_ell_imag(GEN velg, long n, long d_chi, GEN df0)
{
  long i, n_ell = n*d_chi;
  GEN z = cgetg(n_ell + 1, t_VEC);
  GEN df0l = shifti(df0, 1), ellg = mkvec2(gen_1, gen_1);
  for (i=1; i<=n; i++) df0l = muliu(df0l, gel(velg, i)[1]);
  for (i=1; i<=n_ell; i++) ellg = gel(z, i)= next_ell_imag(ellg, df0l);
  return z;
}

/* U(X)=u(x)+u(X)*X^f+...+f(X)*X^((m-1)f) or u(x)-u(X)*X^f+...
 * U(X)V(X)=u(X)V(X)(1+X^f+...+X^((m-1)f))
 *         =w_0+w_1*X+...+w_{f+el-3}*X^(f+el-3)
 * w_i (1 <= i <= f+el-2) are needed.
 * w_{f+el-2}=0 if el-1 == f.
 * W_i = w_i + w_{i+el-1} (1 <= i <= f-1). */
static GEN
gauss_Flx_mul(ulong f, GEN elg, GEN ellg)
{
  pari_sp av = avma;
  ulong el = elg[1], g_el= elg[2];
  ulong el_1 = el-1, f2 = f<<1, lv = el_1, lu = f, m = el_1/f;
  ulong ell = itou(gel(ellg, 1)), g_ell = itou(gel(ellg, 2));
  ulong z_2f = Fl_powu(g_ell, (ell - 1) / f2, ell);
  ulong z_el = Fl_powu(g_ell, (ell - 1) / el, ell);
  ulong i, i2, gi;
  GEN W = cgetg(f+1, t_VECSMALL), vz_2f, vz_el;
  GEN u = cgetg(lu+2, t_VECSMALL), v = cgetg(lv+2, t_VECSMALL), w0;

  u[1] = evalsigne(1);
  v[1] = evalsigne(1);
  vz_2f = Fl_powers(z_2f, f2-1, ell);
  vz_el = Fl_powers(z_el, el_1, ell);
  for (i=i2=0; i<lu; i++)
  {
    long j2; /* i2=(i*i)%f2, gi=g_el^i */
    j2 = i2?f2-i2:i2;
    u[2+i] = vz_2f[1+j2];
    if ((i2+=i+i+1)>=f2) i2-=f2; /* same as i2%=f2 */
  }
  for (gi=1,i=i2=0; i<lv; i++)
  {
    v[2+i] = Fl_mul(vz_2f[1+i2], vz_el[1+gi], ell);
    gi = Fl_mul(gi, g_el, el);
    if ((i2+=i+i+1)>=f2) i2%=f2; /* i2-=f2 does not work */
  }
  w0 = Flx_mul(u, v, ell) + 1;
  if (m==1)
  { /* el_1=f */
    for (i=1; i<f; i++) W[i] = Fl_add(w0[i], w0[i+lv], ell);
    W[f] = w0[f];
  }
  else
  {
    ulong start = 1+f, end = f+el-1;
    GEN w = cgetg(end+1, t_VECSMALL);
    for (i=1; i<end; i++) w[i] = w0[i];
    w[end] = 0;
    for (i=1; i<m; i++, start+=f)
      w = both_odd(f,i)? Flv_shift_sub(w, w0, ell, start, end)
                       : Flv_shift_add(w, w0, ell, start, end);
    for (i=0; i<f; i++) W[1+i] = Fl_add(w[1+i], w[1+i+lv], ell);
  }
  for (i=i2=1; i<f; i++)
  {
    W[i]=Fl_mul(W[1+i], vz_2f[1+i2], ell);
    if ((i2+=i+i+1)>=f2) i2%=f2;
  }
  /* W[r]=tau_{LL}^{sigma_r}, 1<= r <= f-1 */
  return gerepilecopy(av, Flv_to_ZV(W));
}

static GEN
gauss_ZX_mul(ulong f, GEN elg, GEN ellg)
{
  pari_sp av = avma, av2;
  ulong el, g_el, el_1, f2, lv, lu, m, i, i2, gi;
  GEN  ell = gel(ellg, 1), g_ell, ell_1, z_2f, z_el, W, vz_2f, vz_el, u, v, w0;

  if (lgefint(ell) == 3) return gauss_Flx_mul(f, elg, ellg);
  g_ell = gel(ellg, 2); ell_1 = subiu(ell, 1);
  el = elg[1]; g_el = elg[2]; el_1 = el-1;
  f2 = f<<1; lv=el_1; lu=f; m=el_1/f;
  z_2f = Fp_pow(g_ell, diviuexact(ell_1, f2), ell);
  vz_2f = Fp_powers(z_2f, f2-1, ell);
  W = cgetg(f+1, t_VEC);
  av2 = avma;
  z_el = Fp_pow(g_ell, diviuexact(ell_1, el), ell);
  vz_el = Fp_powers(z_el, el_1, ell);
  u = cgetg(lu+2, t_POL); u[1] = evalsigne(1) | evalvarn(0);
  v = cgetg(lv+2, t_POL); v[1] = evalsigne(1) | evalvarn(0);
  for (gi=1,i=i2=0; i<lu; i++)
  {
    long j2; /* i2=(i*i)%f2, gi=g_el^i */
    j2 = i2?f2-i2:i2;
    gel(u, 2+i) = gel(vz_2f, 1+j2);
    if ((i2+=i+i+1)>=f2) i2-=f2;
  }
  for (gi=1,i=i2=0; i<lv; i++)
  {
    gel(v, 2+i) = Fp_mul(gel(vz_2f, 1+i2), gel(vz_el, 1+gi), ell);
    gi = Fl_mul(gi, g_el, el);
    if ((i2+=i+i+1)>=f2) i2%=f2;
  }
  w0 = gerepileupto(av2, FpX_mul(u, v, ell)) + 1; av2 = avma;
  if (m==1)
  {
    for (i=1; i < f; i++) gel(W,i) = Fp_add(gel(w0, i), gel(w0, i+lv), ell);
    gel(W, f) = gel(w0, f);
  }
  else
  {
    ulong start = 1+f, end = f+el-1;
    GEN w = cgetg(end+1, t_VEC);
    for (i=1; i<end; i++) gel(w, i) = gel(w0, i);
    gel(w, end) = gen_0;
    for (i=1; i<m; i++, start+=f)
    {
      w = both_odd(f,i)? FpV_shift_sub(w, w0, ell, start, end)
                       : FpV_shift_add(w, w0, ell, start, end);
      if ((i & 7) == 0) w = gerepilecopy(av2, w);
    }
    for (i = 1; i <= f; i++) gel(W, i) = addii(gel(w, i), gel(w, i+lv));
  }
  for (i = i2 = 1; i < f; i++)
  {
    gel(W, i) = Fp_mul(gel(W, 1+i), gel(vz_2f, 1+i2), ell);
    if ((i2+=i+i+1) >= f2) i2 %= f2;
  }
  return gerepilecopy(av, W);  /* W[r]=tau_{LL}^{sigma_r}, 1<= r <= f-1 */
}

/* fast but consumes memory */
static GEN
gauss_el_vell(ulong f, GEN elg, GEN vellg, GEN vz_2f)
{
  pari_sp av = avma, av2;
  ulong el = elg[1], g_el = elg[2], el_1 = el-1;
  ulong lv=el_1, f2=f<<1, lu=f, m=el_1/f;
  GEN W = cgetg(f+1, t_VEC), vz_el, u, v, w0, M;
  ulong i, i2, gi;

  av2 = avma;
  vz_el = vz_vell(el, vellg, &M);
  u = cgetg(lu+2, t_POL); u[1] = evalsigne(1) | evalvarn(0);
  v = cgetg(lv+2, t_POL); v[1] = evalsigne(1) | evalvarn(0);
  for (i=i2=0; i<lu; i++)
  {
    long j2; /* i2=(i*i)%f2, gi=g_el^i */
    j2 = i2?f2-i2:i2;
    gel(u, 2+i) = gel(vz_2f, 1+j2);
    if ((i2+=i+i+1)>=f2) i2%=f2;
  }
  for (gi=1,i=i2=0; i<lv; i++)
  {
    gel(v, 2+i) = Fp_mul(gel(vz_2f, 1+i2), gel(vz_el, 1+gi), M);
    gi = Fl_mul(gi, g_el, el);
    if ((i2+=i+i+1)>=f2) i2%=f2;
  }
  M = gclone(M);
  w0 = gerepileupto(av2, FpX_mul(u, v, M)) + 1;
  u = M; M = icopy(M); gunclone(u);
  av2 = avma;
  if (m==1)
  { /* el_1=f */
    for (i=1; i < f; i++) gel(W,i) = Fp_add(gel(w0, i), gel(w0, i+lv), M);
    gel(W, f) = gel(w0, f);
  }
  else
  {
    ulong start = 1+f, end = f+el-1;
    GEN w = cgetg(end+1, t_VEC);
    for (i=1; i<end; i++) gel(w, i) = gel(w0, i);
    gel(w, end) = gen_0;
    for (i=1; i<m; i++, start+=f)
    {
      w = both_odd(f,i)? FpV_shift_sub(w, w0, M, start, end)
                       : FpV_shift_add(w, w0, M, start, end);
      if ((i & 7) == 0) w = gerepilecopy(av2, w);
    }
    for (i = 1; i <= f; i++) gel(W, i) = Fp_add(gel(w, i), gel(w, i+lv), M);
  }
  for (i = i2 = 1; i < f; i++)
  {
    gel(W, i) = Fp_mul(gel(W, 1+i), gel(vz_2f, 1+i2), M);
    if ((i2+=i+i+1) >= f2) i2 %= f2;
  }
  return gerepilecopy(av, W);  /* W[r]=tau_{LL}^{sigma_r}, 1<= r <= f-1 */
}

static GEN
norm_chi(GEN K, GEN TAU, ulong p, long d_pow, GEN ell, long j0)
{
  pari_sp av = avma;
  GEN H = K_get_H(K);
  ulong d_K, f_K = K_get_f(K), h = K_get_h(K), g_K = K_get_g(K);
  ulong i, j, gi, pd = upowuu(p, d_pow), d_chi = K_get_dchi(K);
  GEN z = const_vec(d_chi, gen_1);
  GEN e_chi = get_e_chi(K, j0, pd, &d_K);

  for (gi=1, i=0; i<d_K; i++)
  {
    GEN y = gen_1;
    for (j=1; j<=h; j++)
      y = Fp_mul(y, gel(TAU, Fl_mul(gi, H[j], f_K)), ell);
    gi = Fl_mul(gi, g_K, f_K);
    for (j=1; j<=d_chi; j++)
    {
      GEN y2 = Fp_powu(y, e_chi[(i+j-1)%d_K], ell);
      gel(z, j) = Fp_mul(gel(z, j), y2, ell);
    }
  }
  return gerepilecopy(av, z);
}

static void
imag_MLLn(long *y, GEN K, ulong p, long d_pow, long n,
    GEN velg, GEN vellg, long j0)
{
  long f = K_get_f(K), d = upowuu(p, d_pow), row = lg(vellg)-1, i, j, k, nz;
  GEN g, z, M, vz_2f = vz_vell(f << 1, vellg, &M);
  for (i=1; i<=n; i++)
  {
    pari_sp av = avma;
    GEN elg = gel(velg, i);
    if (DEBUGLEVEL>1) err_printf("(f,el-1)=(%ld,%ld*%ld)\n", f,(elg[1]-1)/f,f);
    g = gauss_el_vell(f, elg, vellg, vz_2f);
    z = norm_chi(K, g, p, d_pow, M, j0);
    nz = lg(z)-1;
    for (k = 1; k <= nz; k++)
      for (j = 1; j <= row; j++)
        y[(j-1)*row+(i-1)*nz+k-1] = get_y(gel(z, k), gel(vellg, j), d);
    set_avma(av);
  }
}

static void
imag_MLL1(long *y, GEN K, ulong p, long d_pow, GEN velg, GEN vellg, long j0)
{
  long f = K_get_f(K), d = upowuu(p, d_pow);
  GEN elg = gel(velg, 1), ellg = gel(vellg, 1), ell = gel(ellg, 1), g, z;

  if (DEBUGLEVEL>1) err_printf("(f,el-1)=(%ld,%ld*%ld)\n", f, (elg[1]-1)/f, f);
  g = gauss_ZX_mul(f, elg, ellg);
  z = norm_chi(K, g, p, d_pow, ell, j0);
  y[0] = get_y(gel(z, 1), ellg, d);
}

static void
imag_MLL(long *y, GEN K, ulong p, long d_pow, long n, GEN velg, GEN vellg,
    long j0)
{
  pari_sp av = avma;
  long i, j, f = K_get_f(K), d = upowuu(p, d_pow), row = lg(vellg)-1;

  for (j=1; j<=row; j++)
  {
    GEN ellg = gel(vellg, j), ell = gel(ellg, 1);
    for (i=1; i<=n; i++)
    {
      GEN elg = gel(velg, i), g, z;
      ulong k, nz;
      if (DEBUGLEVEL>1) err_printf("(f,el-1)=(%ld,%ld*%ld)\n",f,(elg[1]-1)/f,f);
      g = gauss_ZX_mul(f, elg, ellg);
      z = norm_chi(K, g, p, d_pow, ell, j0);
      nz = lg(z)-1;
      for (k = 1; k <= nz; k++)
        y[(j-1)*row+(i-1)*nz+k-1] = get_y(gel(z, k), ellg, d);
      set_avma(av);
    }
  }
}

/* return an upper bound >= 0 if one was found, otherwise return -1.
 * set chi-part to be (1) if chi is Teichmuller character.
 * B_{1,omega^(-1)} is not p-adic integer. */
static GEN
cyc_imag_MLL(GEN K, ulong p, long d_pow, long j, long flag)
{
  long f = K_get_f(K), d_chi = K_get_dchi(K);
  long n, n0 = 1, n_el = d_pow, d = upowuu(p, d_pow), rank = n_el*d_chi;
  GEN df0, velg = const_vec(n_el, NULL), vellg = NULL;
  GEN oldgr = mkvec2(gen_0, NULL), newgr = mkvec2(gen_0, NULL);
  long *y0 = (long*)stack_calloc(sizeof(long)*rank*rank);

  if (DEBUGLEVEL>1)
    err_printf("cyc_imag_MLL:p=%ld d_pow=%ld deg(K)=%ld cond(K)=%ld avma=%ld\n",
        p, d_pow, K_get_d(K), f, avma);
  df0 = muluu(d, f%p?f:f/p);
  gel(velg, 1) = next_el_imag(mkvecsmall2(1, 1), f);
  if (flag&USE_FULL_EL)
  {
    for (n=2; n<=n_el; n++) gel(velg, n) = next_el_imag(gel(velg, n-1), f);
    n0 = n_el;
  }
  for (n=n0; n<=n_el; n++) /* loop while structure is unknown */
  {
    pari_sp av2 = avma;
    pari_timer ti;
    long n_ell, m, M;
    GEN y;
    vellg = set_ell_imag(velg, n, d_chi, df0);
    n_ell = lg(vellg)-1;  /* equal to n*d_chi */
    if (DEBUGLEVEL>2) err_printf("velg=%Ps\nvellg=%Ps\n", velg, vellg);
    if (DEBUGLEVEL>2) timer_start(&ti);
    if (n_ell==1)
      imag_MLL1(y0, K, p, d_pow, velg, vellg, j);
    else if (lgefint(gmael(vellg, n, 1))<=3 || (flag&SAVE_MEMORY))
      imag_MLL(y0, K, p, d_pow, n, velg, vellg, j);
    else
      imag_MLLn(y0, K, p, d_pow, n, velg, vellg, j);
    set_avma(av2);
    if (DEBUGLEVEL>2) timer_printf(&ti, "gauss sum");
    y = ary2mat(y0, n_ell);
    if (DEBUGLEVEL>3) err_printf("y=%Ps\n", y);
    y = ZM_snf(y);
    if (DEBUGLEVEL>3) err_printf("y=%Ps\n", y);
    y = make_p_part(y, p, d_pow);
    if (DEBUGLEVEL>3) err_printf("y=%Ps\n", y);
    newgr = structure_MLL(y, d_pow);
    if (DEBUGLEVEL>3)
      err_printf("d_pow=%ld d_chi=%ld old=%Ps new=%Ps\n",d_pow,d_chi,oldgr,newgr);
    if (equalsi(d_pow*d_chi, gel(newgr, 1))) break;
    if ((m = find_del_el(&oldgr, newgr, n, n_el, d_chi)))
    { M = m = delete_el(velg, m); n--; }
    else
    { M = n+1; m = n; }
    gel(velg, M) = next_el_imag(gel(velg, m), f);
  }
  return get_str(newgr);
}

/* When |A_psi|=p^e, A_psi=(p^e1,...,p^er) (psi=chi^j),
 *  return vec[e, [e1, ... ,er], 1].
 * If gr str is not determined, return vec[e, [], 1].
 * If |A_chi|=1, return vec[0, [], 1].
 * If |A_chi|=p, return vec[1, [1], 1].
 * If e is not determined, return vec[-1, [], 1].
 * If psi is Teichmuller, return vec[0, [], 1].
 * B_{1,omega^(-1)} is not p-adic integer. */
static GEN
cyc_imag(GEN K, GEN B, GEN p, long j, GEN powp, long flag)
{
  pari_sp av = avma;
  GEN MinPol = gel(K, 3), Chi = gel(K, 2), B1, B2, gr;
  long x, d_K = K_get_d(K), f_K = K_get_f(K), d_chi = K_get_dchi(K);

  if (f_K == d_K+1 && equaliu(p, f_K) && j == 1) /* Teichmuller */
    return mkvec3(gen_0, nullvec(), gen_1);
  B1 = FpX_rem(ZX_ber_conj(B, j, d_K), MinPol, powp);
  B2 = FpX_rem(ZX_ber_den(Chi, j, d_K), MinPol, powp);
  if (degpol(B1)<0 || degpol(B2)<0)
  {
    set_avma(av);
    return mkvec3(gen_m1, nullvec(), gen_1); /* B=0(mod p^pow) */
  }
  x = ZX_pval(B1, p) - ZX_pval(B2, p);
  set_avma(av);
  if (x<0) pari_err_BUG("subcyclopclgp [Bernoulli number]");
  if (DEBUGLEVEL && x) verbose_output(K, p, x, j);
  if (x==0) return mkvec3(gen_0, nullvec(), gen_1); /* trivial */
  if (x==1) return mkvec3(utoi(d_chi), onevec(d_chi), gen_1);
  if ((flag&USE_MLL)==0) return mkvec3(utoi(x*d_chi), nullvec(), gen_1);
  gr = d_K == 2? cyc_buch(-f_K, p, x): cyc_imag_MLL(K, itou(p), x, j, flag);
  return gerepilecopy(av, mkvec3(utoipos(d_chi * x), gr, gen_1));
}

/* handle representatives of all injective characters, d_chi=[Q_p(zeta_d):Q_p],
 * d=d_K */
static GEN
pclgp_cyc_imag(GEN K, GEN p, long start_pow, long max_pow, long flag)
{
  GEN C = gel(K, 5), Chi = gel(K, 2);
  long n_conj = K_get_nconj(K), d_K = K_get_d(K), f_K = K_get_f(K);
  long i, pow, n_done = 0;
  GEN gr = nullvec(), Done = const_vecsmall(n_conj, 0);
  GEN B = zx_ber_num(Chi, f_K, d_K), B_num;

  if (lgefint(p)==3 && n_conj>10) /* mark trivial chi-part by pre-calculation */
  {
    ulong up = itou(p);
    GEN minpol = ZX_to_Flx(gel(K, 3), up);
    for (i=1; i<=n_conj; i++)
    {
      pari_sp av = avma;
      long degB;
      B_num = Flx_rem(Flx_ber_conj(B, C[i], d_K, up), minpol, up);
      degB = degpol(B_num);
      set_avma(av);
      if (degB<0) continue;
      Done[i] = 1;
      if (++n_done == n_conj) return gr;
    }
  }
  for (pow = start_pow; pow<=max_pow; pow++)
  {
    GEN powp = powiu(p, pow);
    for (i = 1; i <= n_conj; i++)
    {
      GEN z;
      if (Done[i]) continue;
      z = cyc_imag(K, B, p, C[i], powp, flag);
      if (equalim1(gel(z, 1))) continue;
      Done[i] = 1;
      if (!isintzero(gel(z, 1))) gr = vec_append(gr, z);
      if (++n_done == n_conj) return gr;
    }
  }
  pari_err_BUG("pclgp_cyc_imag: max_pow is not enough");
  return NULL; /*LCOV_EXCL_LINE*/
}

static GEN
gather_part(GEN g, long sgn)
{
  long i, j, l = lg(g), ord = 0, flag = 1;
  GEN z2 = cgetg(l, t_VEC);
  for (i = j = 1; i < l; i++)
  {
    GEN t = gel(g,i);
    if (equaliu(gel(t, 3), sgn))
    {
      ord += itou(gel(t, 1));
      if (lg(gel(t, 2)) == 1) flag = 0;
      gel(z2, j++) = gel(t, 2);
    }
  }
  if (flag==0 || ord==0) z2 = nullvec();
  else
  {
    setlg(z2, j); z2 = shallowconcat1(z2);
    ZV_sort_inplace(z2); vecreverse_inplace(z2);
  }
  return mkvec2(utoi(ord), z2);
}

#ifdef DEBUG
static void
handling(GEN K)
{
  long d_K = K_get_d(K), f_K = K_get_f(K), s_K = K_get_s(K), g_K = K_get_g(K);
  long d_chi = K_get_dchi(K);
  err_printf("  handling %s cyclic subfield K,\
      deg(K)=%ld, cond(K)=%ld g_K=%ld d_chi=%ld H=%Ps\n",
      s_K? "a real": "an imaginary",d_K,f_K,g_K,d_chi,zv_to_ZV(gmael3(K,1,1,1)));
}
#endif

/* HH a t_VECSMALL listing group generators
 * Aoki and Fukuda, LNCS vol.4076 (2006), 56-74. */
static GEN
pclgp(GEN p0, long f, GEN HH, long degF, long flag)
{
  long start_pow, max_pow, ip, lp, i, n_f;
  GEN vH1, z, vData, cycGH, vp = typ(p0) == t_INT? mkvec(p0): p0;

  vH1 = GHinit(f, HH, &cycGH); n_f = lg(vH1)-1;
#ifdef DEBUG
  err_printf("F is %s, deg(F)=%ld, ", srh_1(HH)? "real": "imaginary", degF);
  err_printf("cond(F)=%ld, G(F/Q)=%Ps\n",f, cycGH);
  err_printf("F has %ld cyclic subfield%s except for Q.\n", n_f,n_f>1?"s":"");
#endif

  lp = lg(vp); z = cgetg(lp, t_MAT);
  for (ip = 1; ip < lp; ip++)
  {
    pari_sp av = avma;
    long n_sub=0, n_chi=0;
    GEN gr=nullvec(), p = gel(vp, ip), zi;
    /* find conductor e of cyclic subfield K and set the subgroup HE of (Z/eZ)^*
     * corresponding to K */
    set_p_f(p, f, &start_pow, &max_pow);
    vData = const_vec(degF, NULL);

    for (i=1; i<=n_f; i++) /* prescan. set Teichmuller */
    {
      GEN H1 = gel(vH1, i);
      long d_K = _get_d(H1), f_K = _get_f(H1), g_K = _get_g(H1);

      if (f_K == d_K+1 && equaliu(p, f_K)) /* found K=Q(zeta_p) */
      {
        pari_timer ti;
        GEN pnmax = powiu(p, max_pow), vNewton, C, MinPol;
        long d_chi = 1, n_conj = eulerphiu(d_K);
        ulong pmodd = umodiu(p, d_K);

        C = set_C(pmodd, d_K, d_chi, n_conj);
        MinPol = set_minpol_teich(g_K, p, max_pow);
        if (DEBUGLEVEL>3) timer_start(&ti);
        vNewton = FpX_Newton(MinPol, d_K+1, pnmax);
        if (DEBUGLEVEL>3)
          timer_printf(&ti, "FpX_Newton: teich: %ld %ld", degpol(MinPol), d_K);
        gel(vData, d_K) = mkvec4(MinPol, vNewton, C,
                                 mkvecsmall2(d_chi, n_conj));
        break;
      }
    }

    for (i=1; i<=n_f; i++) /* handle all cyclic K */
    {
      GEN H1 = gel(vH1, i), K, z1, Chi;
      long d_K = _get_d(H1), s_K = _get_s(H1);
      pari_sp av2;

      if ((flag&SKIP_PROPER) && degF != d_K) continue;
      if (!gel(vData, d_K))
      {
        pari_timer ti;
        GEN pnmax = powiu(p, max_pow), vNewton, C, MinPol;
        ulong pmodd = umodiu(p, d_K);
        long d_chi = order_f_x(d_K, pmodd), n_conj = eulerphiu(d_K)/d_chi;

        C = set_C(pmodd, d_K, d_chi, n_conj);
        MinPol = set_minpol(d_K, p, max_pow, n_conj);
        if (DEBUGLEVEL>3) timer_start(&ti);
        /* vNewton[2+i] = vNewton[2+i+d_K]. We need vNewton[2+i] for
         * 0 <= i < d_K. But vNewton[2+d_K-1] may be 0 and will be deleted.
         * So we need vNewton[2+d_K] not to delete vNewton[2+d_K-1]. */
        vNewton = FpX_Newton(MinPol, d_K+1, pnmax);
        if (DEBUGLEVEL>3)
          timer_printf(&ti, "FpX_Newton: %ld %ld", degpol(MinPol), d_K);
        gel(vData, d_K) = mkvec4(MinPol, vNewton, C,
                                 mkvecsmall2(d_chi, n_conj));
      }
      av2 = avma;
      Chi = s_K? NULL: get_chi(H1);
      K = shallowconcat(mkvec2(H1, Chi), gel(vData, d_K));
#ifdef DEBUG
      handling(K);
#endif
      if (s_K && !(flag&NO_PLUS_PART))
        z1 = pclgp_cyc_real(K, p, max_pow, flag);
      else if (!s_K && !(flag&NO_MINUS_PART))
        z1 = pclgp_cyc_imag(K, p, start_pow, max_pow, flag);
      else { set_avma(av2); continue; }
      n_sub++; n_chi += gmael(vData, d_K, 4)[2]; /* += n_conj */
      if (lg(z1) == 1) set_avma(av2);
      else gr = gerepilecopy(av2, shallowconcat(gr, z1));
    }
    zi = mkcol(p);
    zi = vec_append(zi, (flag&NO_PLUS_PART)?nullvec():gather_part(gr, 0));
    zi = vec_append(zi, (flag&NO_MINUS_PART)?nullvec():gather_part(gr, 1));
    zi = shallowconcat(zi, mkcol3(cycGH, utoi(n_sub), utoi(n_chi)));
    gel(z, ip) = gerepilecopy(av, zi);
  }
  return typ(p0) == t_INT? shallowtrans(gel(z,1)): shallowtrans(z);
}

static GEN
reduce_gcd(GEN x1, GEN x2)
{
  GEN d = gcdii(x1, x2);
  x1 = diviiexact(x1, d);
  x2 = diviiexact(x2, d);
  return mkvec2(x1, x2);
}

/* norm of x0 (= pol of zeta_d with deg <= d-1) by g of order n
 * x0^{1+g+g^2+...+g^(n-1)} */
static GEN
ber_norm_cyc(GEN x0, long g, long n, long d)
{
  pari_sp av = avma;
  long i, ei, di, fi = 0, l = ulogint(n, 2);
  GEN xi = x0;
  ei = 1L << l; di = n / ei;
  for (i = 1; i <= l; i++)
  {
    if (odd(di)) fi += ei;
    ei = 1L << (l-i); di = n / ei;
    xi = ZX_mod_Xnm1(ZX_mul(xi, ber_conj(xi, Fl_powu(g, ei, d), d)), d);
    if (odd(di))
      xi = ZX_mod_Xnm1(ZX_mul(xi, ber_conj(x0, Fl_powu(g, fi, d), d)), d);
  }
  return gerepilecopy(av, xi);
}

/* x0 a ZX of deg < d */
static GEN
ber_norm_by_cyc(GEN x0, long d, GEN MinPol)
{
  pari_sp av=avma;
  GEN x = x0, z = znstar(utoi(d)), cyc = gel(z, 2), gen = gel(z, 3);
  long i, l = lg(cyc);
  pari_timer ti;

  if (DEBUGLEVEL>1) timer_start(&ti);
  for (i = 1; i < l; i++)
    x = ber_norm_cyc(x, itou(gmael(gen, i, 2)), itou(gel(cyc, i)), d);
  if (DEBUGLEVEL>1) timer_printf(&ti, "ber_norm_by_cyc [ber_norm_cyc]");
  x = ZX_rem(x, MinPol);  /* slow */
  if (DEBUGLEVEL>1) timer_printf(&ti, "ber_norm_by_cyc [ZX_rem]");
  if (lg(x) != 3) pari_err_BUG("subcyclohminus [norm of Bernoulli number]");
  return gerepilecopy(av, gel(x, 2));
}

/* MinPol = polcyclo(d_K, 0).
 * MinPol = fac*cofac (mod p).
 * B is zv.
 * K : H1, MinPol, [fac, cofac], C, [d_chi, n_conj] */
static long
ber_norm_by_val(GEN K, GEN B, GEN p)
{
  pari_sp av = avma;
  GEN MinPol = gel(K, 2), C = gel(K, 4);
  GEN vfac = gel(K, 3), fac = gel(vfac, 1), cofac = gel(vfac, 2);
  long d_chi = K_get_dchi(K), n_conj = K_get_nconj(K), d_K = K_get_d(K);
  long i, r, n_done = 0, x = 0, dcofac = degpol(cofac);
  GEN pr, Done;

  Done = const_vecsmall(n_conj, 0);
  if (lgefint(p)==3)
  { /* mark trivial chi-part by pre-calculation */
    ulong up = itou(p);
    GEN facs = ZX_to_Flx(fac, up);
    for (i = 1; i <= n_conj; i++)
    {
      pari_sp av2 = avma;
      GEN B_conj = Flx_rem(Flx_ber_conj(B, C[i], d_K, up), facs, up);
      long degB = degpol(B_conj);
      set_avma(av2); if (degB < 0) continue;
      Done[i] = 1; if (++n_done == n_conj) return gc_long(av, x);
    }
  }
  else
  {
    for (i = 1; i <= n_conj; i++)
    {
      pari_sp av2 = avma;
      GEN B_conj = FpX_rem(FpX_ber_conj(B, C[i], d_K, p), fac, p);
      long degB = degpol(B_conj);
      set_avma(av2); if (degB < 0) continue;
      Done[i] = 1; if (++n_done == n_conj) return gc_long(av, x);
    }
  }
  for (pr = p, r = 2; r; r <<= 1)
  {
    GEN polr;
    pr = sqri(pr); /* p^r */
    polr = (dcofac==0)? FpX_red(MinPol, pr)
                      : gel(ZpX_liftfact(MinPol, vfac, pr, p, r), 1);
    for (i = 1; i <= n_conj; i++)
    {
      pari_sp av2 = avma;
      GEN B_conj;
      long degB;
      if (Done[i]) continue;
      B_conj = FpX_rem(FpX_ber_conj(B, C[i], d_K, pr), polr, pr);
      degB = degpol(B_conj);
      set_avma(av2); if (degB < 0) continue;
      x += d_chi * ZX_pval(B_conj, p);
      Done[i] = 1; if (++n_done == n_conj) return gc_long(av, x);
    }
  }
  pari_err_BUG("ber_norm_by_val"); return 0;/*LCOV_EXCL_LINE*/
}

/* n > 2, p = odd prime not dividing n, e > 0, pe = p^e; d = n*p^e
 * return generators of the subgroup H of (Z/dZ)^* corresponding to
 * Q(zeta_{p^e}): H = {1<=a<=d | gcd(a,n)=1, a=1(mod p^e)} */
static GEN
znstar_subgr(ulong n, ulong pe, ulong d)
{
  GEN z = znstar(utoi(n)), g = gel(z, 3), G;
  long i, l = lg(g);
  G = cgetg(l, t_VECSMALL);
  for (i=1; i<l; i++) G[i] = u_chinese_coprime(itou(gmael(g,i,2)), 1, n, pe, d);
  return mkvec2(gel(z,2), G);
}

/* K is a cyclic extension of degree n*p^e (n>=4 is even).
 * x a ZX of deg < n*p^e. */
static long
ber_norm_with_val(GEN x, long n, ulong p, ulong e)
{
  pari_sp av = avma;
  long i, j, r, degx, pe = upowuu(p, e), d = n*pe;
  GEN z, gr, gen, y = cgetg(pe+2, t_POL), MinPol = polcyclo(n, 0);
  y[1] = evalsigne(1) | evalvarn(0);
  z = znstar_subgr(n, pe, d);
  gr = gel(z, 1); gen = gel(z, 2); r = lg(gr)-1;
  for (i=1; i<=r; i++)
    x = ber_norm_cyc(x, itou(gel(gen, i)), itou(gel(gr, i)), d);
  degx = degpol(x);
  for (j=0; j<pe; j++)
  {
    GEN t = pol_zero(n), z;
    long a = j; /* a=i*pe+j */
    for (i=0; i<n; i++)
    {
      if (a>degx) break;
      gel(t, 2+a%n) = gel(x, 2+a);
      a += pe;
    }
    z = ZX_rem(ZX_renormalize(t, 2+n), MinPol);
    if (degpol(z)<0) gel(y, 2+j) = gen_0;
    else if (degpol(z)==0) gel(y, 2+j) = gel(z, 2);
    else pari_err_BUG("ber_norm_subgr");
  }
  y = ZX_renormalize(y, pe+2);
  if (e>1) y = ZX_rem(y, polcyclo(pe, 0));
  return gc_long(av,  ZX_p_val(y, p, e));
}

/* K is a cyclic extension of degree 2*p^e. x a ZX of deg < 2*p^e. In most
 * cases, deg(x)=2*p^e-1. But deg(x) can be any value in [0, 2*p^e-1]. */
static long
ber_norm_with_val2(GEN x, ulong p, ulong e)
{
  pari_sp av = avma;
  long i, d = degpol(x), pe = upowuu(p, e);
  GEN y = pol_zero(pe);
  if (d == 2*pe-1)
  {
    for (i = 0; i < pe; i++)
      gel(y, 2+i) = odd(i)? subii(gel(x, 2+i+pe), gel(x, 2+i))
                          : subii(gel(x, 2+i), gel(x, 2+i+pe));
  }
  else
  {
    for (i = 0; i < pe && i <= d; i++)
      gel(y, 2+i) = odd(i)? negi(gel(x, 2+i)): gel(x, 2+i);
    for (i = pe; i <= d; i++)
      gel(y, 2+i-pe) = odd(i)? subii(gel(y, 2+i-pe), gel(x, 2+i))
                             : addii(gel(y, 2+i-pe), gel(x, 2+i));
  }
  y = ZX_renormalize(y, 2+pe);
  if (e > 1) y = ZX_rem(y, polcyclo(pe, 0));
  return gc_long(av, ZX_p_val(y, p, e));
}

/* K : H1, MinPol, [fac, cofac], C, [d_chi, n_conj] */
static GEN
ber_cyc5(GEN K, GEN p)
{
  pari_sp av = avma;
  GEN MinPol = gel(K, 2), H = K_get_H(K);
  long d = K_get_d(K), f = K_get_f(K), h = K_get_h(K), g = K_get_g(K);
  GEN x, x1, x2, y, B = const_vecsmall(d+1, 0);
  long i, j, gi, e, f2 = f>>1, dMinPol = degpol(MinPol), chi2 = -1, *B2 = B+2;

  /* get_chi inlined here to save memory */
  for (j=1; j<=h; j++) /* i = 0 */
  {
    if (H[j] == 2) chi2 = 0;
    if (H[j] <= f2) B2[0]++; /* Chi[H[j]] = 0 */
  }
  for (i = 1, gi = g; i < d; i++)
  {
    for (j=1; j<=h; j++)
    {
      long t = Fl_mul(gi, H[j], f); /* Chi[t] = i */
      if (t == 2) chi2 = i;
      if (t <= f2) B2[i]++;
    }
    gi = Fl_mul(gi, g, f);
  }
  y = zx_to_ZX(zx_renormalize(B, d+2));

  if (p)
  {
    ulong n;
    e = u_pvalrem(d, p, &n);
    if (e == 0)
      x1 = utoi(ber_norm_by_val(K, B, p));
    else if (n > 2)
      x1 = utoi(ber_norm_with_val(y, n, itou(p), e));
    else
      x1 = utoi(ber_norm_with_val2(y, itou(p), e));
  }
  else
  {
    if (dMinPol > 100)
      x1 = ber_norm_by_cyc(y, d, MinPol);
    else
      x1 = ZX_resultant(MinPol, ZX_rem(y, MinPol));
  }

  if (chi2 < 0) /* chi2 = Chi[2] */
    x2 = shifti(gen_1, 2*dMinPol);
  else if (chi2 == 0)
    x2 = shifti(gen_1, dMinPol);
  else
  {
    long e = d/ugcd(chi2, d);
    x2 = powiu(polcyclo_eval(e, gen_2), eulerphiu(d)/eulerphiu(e));
    x2 = shifti(x2, dMinPol);
  }
  if (p) x = stoi(itou(x1)-Z_pval(x2, p)); else x = reduce_gcd(x1, x2);
  return gerepilecopy(av, x);
}

/*  Hirabayashi-Yoshino, Manuscripta Math. vol.60, 423-436 (1988), Theorem 1
 *
 *  F is a subfield of Q(zeta_f)
 *  f=p^a => Q=1
 *  If F=Q(zeta_f), Q=1 <=> f=p^a
 *  If f=4*p^a, p^a*q^b (p,q are odd primes), Q=2 <=> [Q(zeta_f):F] is odd */
static long
unit_index(ulong d, ulong f)
{
  ulong r, d_f;
  GEN fa = factoru(f), P = gel(fa, 1), E = gel(fa, 2); r = lg(P)-1;
  if (r==1) return 1;  /* f=P^a */
  d_f = eulerphiu_fact(fa);
  if (d==d_f) return 2;  /* F=Q(zeta_f) */
  if (r==2 && ((P[1]==2 && E[1]==2) || P[1]>2)) return odd(d_f/d)?2:1;
  return 0;
}

/* Compute relative class number h of the subfield K of Q(zeta_f)
 * corresponding to the subgroup HH of (Z/fZ)^*.
 * If p!=NULL, then return valuation(h,p). */
static GEN
rel_class_num(long f, GEN HH, long degF, GEN p)
{
  long i, n_f, W, Q;
  GEN vH1, vData, x, z = gen_1, z1 = gen_0, z2 = mkvec2(gen_1, gen_1);

  vH1 = GHinit(f, HH, NULL); n_f = lg(vH1)-1;
  vData = const_vec(degF, NULL);
  for (i=1; i<=n_f; i++)
  {
    GEN H1 = gel(vH1, i), K;
    long d_K = _get_d(H1), s = _get_s(H1);

    if (s) continue;  /* F is real */
#ifdef DEBUG
    err_printf("  handling %s cyclic subfield K, deg(K)=%ld, cond(K)=%ld\n",
        s? "a real": "an imaginary", d_K, _get_f(H1));
#endif
    if (!gel(vData, d_K))
    {
      GEN C, MinPol, fac, cofac;
      ulong d_chi, n_conj;
      MinPol = polcyclo(d_K,0);
      if (p && umodui(d_K, p))
      {
        ulong pmodd = umodiu(p, d_K);
        GEN MinPol_p = FpX_red(MinPol, p);
        d_chi = order_f_x(d_K, pmodd);
        n_conj = eulerphiu(d_K)/d_chi;
        if (n_conj==1) fac = MinPol_p;  /* polcyclo(d_K) is irred mod p */
        else fac = FpX_one_cyclo(d_K, p);
        cofac = FpX_div(MinPol_p, fac, p);
        C = set_C(pmodd, d_K, d_chi, n_conj);
      }
      else
      {
        fac = cofac = C = NULL;
        d_chi = n_conj = 0;
      }
      gel(vData, d_K) = mkvec5(MinPol, mkvec2(fac, cofac), C,
                               NULL, mkvecsmall2(d_chi, n_conj));
    }
    K = vec_prepend(gel(vData, d_K), H1);
    z = ber_cyc5(K, p);
    if (p) z1 = addii(z1, z);
    else
    {
      gel(z2, 1) = mulii(gel(z2, 1), gel(z, 1));
      gel(z2, 2) = mulii(gel(z2, 2), gel(z, 2));
    }
  }
  W = root_of_1(f, HH);
  if (p) return addiu(z1, z_pval(W, p));
  Q = unit_index(degF, f);
  x = dvmdii(muliu(gel(z2,1), 2 * W), gel(z2,2), &z1);
  if (signe(z1)) pari_err_BUG("subcyclohminus [norm of Bernoulli number]");
  if (!Q && mpodd(x)) Q = 2; /* FIXME: can this happen ? */
  if (Q == 1) x = shifti(x, -1);
  return mkvec2(x, utoi(Q));
}

static void
checkp(const char *fun, long degF, GEN p)
{
  if (!BPSW_psp(p)) pari_err_PRIME(fun, p);
  if (equaliu(p, 2)) pari_err_DOMAIN(fun,"p","=", gen_2, p);
  if (degF && dvdsi(degF, p)) errpdiv(fun, p, degF);
}

/* if flag is set, handle quadratic fields specially (don't set H) */
static long
subcyclo_init(const char *fun, GEN FH, long *pdegF, GEN *pH, long flag)
{
  long f = 0, degF = 2;
  GEN F = NULL, H = NULL;
  if (typ(FH) == t_POL)
  {
    degF = degpol(FH);
    if (degF < 1 || !RgX_is_ZX(FH)) pari_err_TYPE(fun, FH);
    if (flag && degF == 2)
    {
      F = coredisc(ZX_disc(FH));
      if (is_bigint(F))
        pari_err_IMPL(stack_sprintf("conductor f > %lu in %s", ULONG_MAX, fun));
      f = itos(F); if (f == 1) degF = 1;
    }
    else
    {
      GEN z, bnf = Buchall(pol_x(fetch_var()), 0, DEFAULTPREC);
      z = rnfconductor(bnf, FH); H = gel(z,3);
      f = subcyclo_nH(fun, gel(z,2), &H);
      delete_var();
      H = znstar_generate(f, H); /* group elements */
    }
  }
  else
  {
    long l = lg(FH), fH;
    if (typ(FH) == t_INT) F = FH;
    else if (typ(FH) == t_VEC && (l == 2 || l == 3))
    {
      F = gel(FH, 1);
      if (l == 3) H = gel(FH, 2);
    }
    else pari_err_TYPE(fun, FH);
    f = subcyclo_nH(fun, F, &H);
    H = znstar_generate(f, H); /* group elements */
    fH = znstar_conductor(H);
    if (fH == 1) degF = 1;
    else
    {
      if (fH != f) { H = znstar_reduce_modulus(H, fH); f = fH; }
      degF = eulerphiu(f) / zv_prod(gel(H, 2));
    }
  }
  *pH = H; *pdegF = degF; return f;
}

GEN
subcyclopclgp(GEN FH, GEN p, long flag)
{
  pari_sp av = avma;
  GEN H;
  long degF, f = subcyclo_init("subcyclopclgp", FH, &degF, &H, 0);
  if (typ(p) == t_VEC)
  {
    long i, l = lg(p);
    for (i = 1; i < l; i++) checkp("subcyclopclgp", degF, gel(p, i));
    if (f == 1) { set_avma(av); return const_vec(l-1, nullvec()); }
  }
  else
  {
    checkp("subcyclopclgp", degF, p);
    if (f == 1) { set_avma(av); return nullvec(); }
  }
  if (flag >= USE_BASIS) pari_err_FLAG("subcyclopclgp");
  return gerepilecopy(av, pclgp(p, f, H, degF, flag));
}

static GEN
subcycloiwasawa_i(GEN FH, GEN P, long n)
{
  long B, p, f, degF;
  GEN H;
  const char *fun = "subcycloiwasawa";

  if (typ(P) != t_INT) pari_err_TYPE(fun, P);
  if (n < 0) pari_err_DOMAIN(fun, "n", "<", gen_0, stoi(n));
  B = 1L << (BITS_IN_LONG/4);
  if (is_bigint(P) || cmpiu(P, B) > 0)
    pari_err_IMPL(stack_sprintf("prime p > %ld in %s", B, fun));
  p = itos(P);
  if (p <= 1 || !uisprime(p)) pari_err_PRIME(fun, P);
  f = subcyclo_init(fun, FH, &degF, &H, 1);
  if (degF == 1) return NULL;
  if (degF == 2)
  {
    long m = ((f & 3) == 0)? f / 4: f;
    if (H && !srh_1(H)) m = -m;
    if (!n) return quadlambda(p, m);
    return m < 0? imagquadstkpol(p, m, n): realquadstkpol(p, m, n);
  }
  if (p == 2) pari_err_DOMAIN(fun, "p", "=", gen_2, gen_2);
  if (srh_1(H)) return NULL;
  if (degF % p == 0) errpdiv("abeliwasawa", P, degF);
  return abeliwasawa(p, f, H, degF, n);
}
GEN
subcycloiwasawa(GEN FH, GEN P, long n)
{
  pari_sp av = avma;
  GEN z = subcycloiwasawa_i(FH, P, n);
  if (!z) { set_avma(av); return n? nullvec(): mkvec(gen_0); }
  return gerepilecopy(av, z);
}

GEN
subcyclohminus(GEN FH, GEN P)
{
  const char *fun = "subcyclohminus";
  pari_sp av = avma;
  GEN H;
  long degF, f = subcyclo_init(fun, FH, &degF, &H, 0);
  if (P)
  {
    if (typ(P) != t_INT) pari_err_TYPE(fun, P);
    if (isintzero(P)) P = NULL; else checkp(fun, 0, P);
  }
  if (degF == 1 ||  srh_1(H) == 1) return gen_1;
  return gerepilecopy(av, rel_class_num(f, H, degF, P));
}
