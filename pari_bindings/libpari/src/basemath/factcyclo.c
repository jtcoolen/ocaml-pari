/* Copyright (C) 2000, 2012  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/* written by Takashi Fukuda
 *  2019.10.27 : Flx_factcyclo_gen, FpX_factcyclo_gen
 *  2019.10.28 : Flx_factcyclo_lift, FpX_factcyclo_lift
 *  2019.11.3  : Flx_factcyclo_newton, FpX_factcyclo_newton
 *  2019.11.12 : gausspol for prime
 *  2019.11.13 : gausspol for prime power
 *  2019.11.14 : ZpX_roots_nonsep with ZX_Zp_root
 *  2019.11.15 : test ZpX_roots_nonsep with polrootspadic
 *  2019.11.16 : accept variable number
 *  2019.11.17 : gen_ascent
 *  2019.11.20 : ZpX_roots_nonsep with FpX_roots
 *  2021.7.19  : Flx_factcyclo_newton_general
 *  2021.7.22  : Flx_conductor_lift */

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_factcyclo

#define GENERAL            1
#define NEWTON_POWER       2
#define NEWTON_GENERAL     4
#define NEWTON_GENERAL_NEW 8
#define ASCENT            16

#define Flx_polcyclo(n, p) ZX_to_Flx(polcyclo(n, 0), p)
#define FpX_polcyclo(n, p) FpX_red(polcyclo(n, 0), p)

/* 0 <= z[i] <= ULONG_MAX */
static GEN
ZX_to_nx(GEN z)
{
  long i, r = lg(z);
  GEN x = cgetg(r, t_VECSMALL);
  for (i = 2; i < r; i++) x[i] = itou(gel(z, i));
  return x;
}

static long
QX_den_pval(GEN x, GEN p)
{
  long i, vmax = 0, l = lg(x);
  for (i = 2; i < l; i++)
  {
    GEN z = gel(x, i);
    long v;
    if (typ(z)==t_FRAC && (v = Z_pval(gel(z, 2), p)) > vmax) vmax = v;
  }
  return vmax;
}

static long
QXV_den_pval(GEN vT, GEN kT, GEN p)
{
  long k, vmax = 0, l = lg(kT);
  for (k = 1; k < l; k++)
  {
    long v = QX_den_pval(gel(vT, kT[k]), p);
    if (v > vmax) vmax = v;
  }
  return vmax;
}

/* n=el^e, p^d=1 (mod n), d is in [1,el-1], phi(n)=d*f.
 *  E[i] : 1 <= i <= n-1
 *  E[g^i*p^j mod n] = i+1  (0 <= i <= f-1, 0 <= j <= d-1)
 *  We use E[i] (1 <= i <= d). Namely i < el. */
static GEN
set_E(long pmodn, long n, long d, long f, long g)
{
  long i, j;
  GEN E = const_vecsmall(n-1, 0);
  pari_sp av = avma;
  GEN C = Fl_powers(g, f-1, n);
  for (i = 1; i <= f; i++)
  {
    ulong x = C[i];
    for (j = 1; j <= d; j++) { x = Fl_mul(x, pmodn, n); E[x] = i; }
  }
  return gc_const(av, E);
}

/* x1, x2 of the same degree; gcd(p1,p2) = 1, m = p1*p2, m2 = m >> 1*/
static GEN
ZX_chinese_center(GEN x1, GEN p1, GEN x2, GEN p2, GEN m, GEN m2)
{
  long i, l = lg(x1);
  GEN x = cgetg(l, t_POL);
  GEN y1, y2, q1, q2;
  (void)bezout(p1, p2, &q1, &q2);
  y1 = Fp_mul(p2, q2, m);
  y2 = Fp_mul(p1, q1, m);
  for (i = 2; i < l; i++)
  {
    GEN y = Fp_add(mulii(gel(x1, i), y1), mulii(gel(x2, i), y2), m);
    if (cmpii(y, m2) > 0) y = subii(y, m);
    gel(x, i) = y;
  }
  x[1] = x1[1]; return x;
}

/* find n_el primes el such that el=1 (mod n) and el does not divide d(T) */
static GEN
list_el_n(ulong el0, ulong n, GEN d1, long n_el)
{
  GEN v = cgetg(n_el + 1, t_VECSMALL);
  forprime_t T;
  long i;
  u_forprime_arith_init(&T, el0+n, ULONG_MAX, 1, n);
  for (i = 1; i <= n_el; i++)
  {
    ulong el;
    do el = u_forprime_next(&T); while (dvdiu(d1, el));
    v[i] = el;
  }
  return v;
}

/* return min el s.t. 2^63<el and el=1 (mod n). */
static ulong
start_el_n(ulong n)
{
  ulong MAXHLONG = 1UL<<(BITS_IN_LONG-1), el = (MAXHLONG/n)*n + 1;
  if ((el&1)==0) el += n; /* if el is even, then n is odd */
  return el + (n << 1);
}

/* start probably catches d0*T_k(x). So small second is enough. */
static ulong
get_n_el(GEN d0, ulong *psec)
{
  ulong start = ((lgefint(d0)-2)*BITS_IN_LONG)/(BITS_IN_LONG-1)+1, second = 1;
  if (start>10) second++;
  if (start>100)  { start++; second++; }
  if (start>500)  { start++; second++; }
  if (start>1000) { start++; second++; }
  *psec = second; return start;
}

static long
FpX_degsub(GEN P, GEN Q, GEN p)
{
  pari_sp av = avma;
  long d = degpol(FpX_sub(P, Q, p));
  return gc_long(av, d);
}

/* n=odd prime power, F=Q(zeta_n), G=G(F/Q)=(Z/nZ)^*, H=<p>, K <--> H,
 * t_k=Tr_{F/K}(zeta_n^k), T=minimal pol. of t_1 over Q.
 * g is a given generator of G(K/Q).
 * There exists a unique G(x) in Q[x] s.t. t_g=G(t_1).
 * return G(x) mod el assuming that el does not divide d(T), in which case
 * T(x) is separable over F_el and so Vandermonde mod el is regular. */
static GEN
gausspol_el(GEN H, ulong n, ulong d, ulong f, ulong g, ulong el)
{
  ulong j, k, z_n = rootsof1_Fl(n, el);
  GEN vz_n, L = cgetg(1+f, t_VECSMALL), x = cgetg(1+f, t_VECSMALL), X;

  vz_n = Fl_powers(z_n, n-1, el)+1;
  for (k = 0; k < f; k++)
  {
    ulong gk = Fl_powu(g, k, n), t = 0;
    for (j = 1; j <= d; j++)
      t = Fl_add(t, vz_n[Fl_mul(H[j], gk, n)], el);
    L[1+k] = t;
    x[1+(k+f-1)%f] = t;
  }
  X = Flv_invVandermonde(L, 1, el);
  return Flv_to_Flx(Flm_Flc_mul(X, x, el), 0);
}

static GEN
get_G(GEN H, GEN d0, GEN d1, GEN N, long k, ulong *pel, GEN *pM)
{
  long n = N[1], d = N[2], f = N[3], g = N[4], i;
  GEN POL = cgetg(1+k, t_VEC), EL, G, M, x;
  pari_timer ti;

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  EL = list_el_n(*pel, n, d1, k);
  for (i = 1; i <= k; i++)
  {
    ulong el = uel(EL,i);
    x = gausspol_el(H, n, d, f, g, el);
    gel(POL, i) = Flx_Fl_mul(x, umodiu(d0, el), el);
  }
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "get_G : make data k=%ld",k);
  G = nxV_chinese_center(POL, EL, &M);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "get_G : nxV_chinese_center k=%ld",k);
  *pel = EL[k]; *pM = M; return G;
}

static long
Q_size(GEN z)
{
  if (typ(z)==t_INT) return lgefint(z) - 2;
  return maxss(lgefint(gel(z,1)), lgefint(gel(z,2))) - 2;
}
/* return max log_a(|x[i]|), a=2^(BITS_IN_LONG-1) */
static long
ZX_size(GEN x)
{
  long i, l = lg(x), max = 0;
  for (i = 2; i < l; i++)
  {
    long y = lgefint(gel(x,i)) - 2;
    if (y > max) max = y;
  }
  return max;
}

/* d0 is a multiple of (O_K:Z[t_1]). i.e. d0*T_k(x) is in Z[x].
 * d1 has the same prime factors as d(T); d0 d1 = d(T)^2 */
static GEN
get_d0_d1(GEN T, GEN P)
{
  long i, l = lg(P);
  GEN x, y, dT = ZX_disc(T);
  x = y = dT; setsigne(dT, 1);
  for (i = 1; i < l; i++)
    if (odd(Z_lvalrem(dT, P[i], &dT)))
    {
      x = diviuexact(x, P[i]);
      y = muliu(y, P[i]);
    }
  return mkvec2(sqrti(x), sqrti(y));  /* x and y are square */
}

static GEN
gausspol(GEN T, GEN H, GEN N, ulong d, ulong f, ulong g)
{
  long n = N[1], el0 = N[2];
  GEN F, G1, M1, d0, d1, Data, d0d1 = get_d0_d1(T, mkvecsmall(el0));
  ulong el, n_el, start, second;
  pari_timer ti;

  d0 = gel(d0d1,1); /* d0*F is in Z[X] */
  d1 = gel(d0d1,2); /* d1 has same prime factors as disc(T) */
  Data = mkvecsmall4(n, d, f, g);
  start = get_n_el(d0, &second);
  el = start_el_n(n);

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  if (DEBUGLEVEL == 2) err_printf("gausspol:start=(%ld,%ld)\n",start,second);
  G1 = get_G(H, d0, d1, Data, start, &el, &M1);
  for (n_el=second; n_el; n_el++)
  {
    GEN m, G2, M2;
    G2 = get_G(H, d0, d1, Data, n_el, &el, &M2);
    if (FpX_degsub(G1, G2, M2) < 0) break;  /* G1 = G2 (mod M2) */
    if (DEBUGLEVEL == 2)
      err_printf("G1:%ld, G2:%ld\n",ZX_size(G1), ZX_size(G2));
    if (DEBUGLEVEL >= 6) timer_start(&ti);
    m = mulii(M1, M2);
    G2 = ZX_chinese_center(G1, M1, G2, M2, m, shifti(m,-1));
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "ZX_chinese_center");
    G1 = G2; M1 = m;
  }
  F = RgX_Rg_div(G1, d0);
  if (DEBUGLEVEL == 2)
    err_printf("G1:%ld, d0:%ld, M1:%ld, F:%ld\n",
               ZX_size(G1), Q_size(d0), Q_size(M1), ZX_size(F));
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "gausspol");
  return F;
}

/* Data = [H, GH, i_t, d0d1, kT, [n, d, f, n_T, mitk]]
 * v_t_el[k]=t_k mod el, 1<= k <= n-1 */
static GEN
mk_v_t_el(GEN vT, GEN Data, ulong el)
{
  pari_sp av = avma;
  GEN H = gel(Data, 1), GH = gel(Data,2), i_t = gel(Data, 3), N=gel(Data, 6);
  ulong n = N[1],  d = N[2], mitk = N[5], f = N[3], i, k;
  ulong z_n = rootsof1_Fl(n, el);
  GEN vz_n = Fl_powers(z_n, n-1, el)+1;
  GEN v_t_el = const_vecsmall(n-1, 0);

  for (k = 1; k <= mitk; k++)
  {
    if (k > 1 && !isintzero(gel(vT, k))) continue; /* k=1 is always handled */
    for (i=1; i<=f; i++)
    {
      ulong j, t = 0, x = Fl_mul(k, GH[i], n);
      long y = i_t[x]; /* x!=0, y!=0 */
      if (v_t_el[y]) continue;
      for (j = 1; j <= d; j++) t = Fl_add(t, vz_n[Fl_mul(x, H[j], n)], el);
      v_t_el[y] = t;
    }
  }
  return gerepileuptoleaf(av, v_t_el);
}

/* G=[[G_1,...,G_d],M,el]
 * Data = [H, GH, i_t, d0d1, kT, [n, d, f, n_T, mitk]] */
static GEN
get_vG(GEN vT, GEN Data, long n_el, ulong *pel, GEN *pM)
{
  GEN GH = gel(Data, 2), i_t = gel(Data, 3);
  GEN d0 = gmael(Data, 4, 1), d1 = gmael(Data, 4, 2);
  GEN kT = gel(Data, 5), N = gel(Data, 6);
  long n = N[1], f = N[3], n_T = N[4], mitk = N[5];
  GEN G = const_vec(mitk, gen_0), vPOL = cgetg(1+mitk, t_VEC);
  GEN EL, M, X, v_t_el;
  GEN L = cgetg(1+f, t_VECSMALL), x = cgetg(1+f, t_VECSMALL);
  long i, j, k;

  for (k=1; k<=mitk; k++) gel(vPOL, k) = cgetg(1+n_el, t_VEC);
  EL = list_el_n(*pel, n, d1, n_el);
  for (i=1; i<=n_el; i++)
  {
    ulong el = uel(EL,i), d0model = umodiu(d0, el);
    v_t_el = mk_v_t_el(vT, Data, el);
    for (j = 1; j <= f; j++) L[j] = v_t_el[i_t[GH[j]]];
    X = Flv_invVandermonde(L, 1, el);
    for (k = 1; k <= n_T; k++)
    {
      GEN y;
      long itk = kT[k];
      if (!isintzero(gel(vT, itk))) continue;
      for (j = 1; j <= f; j++) x[j] = v_t_el[i_t[Fl_mul(itk, GH[j], n)]];
      y = Flv_to_Flx(Flm_Flc_mul(X, x, el), 0);
      gmael(vPOL, itk, i) = Flx_Fl_mul(y, d0model, el);
    }
  }
  for (k = 1; k <= n_T; k++)
  {
    long itk = kT[k];
    if (!isintzero(gel(vT, itk))) continue;
    gel(G, itk) = nxV_chinese_center(gel(vPOL, itk), EL, &M);
  }
  *pel = EL[n_el]; *pM = M; return G;
}

/* F=Q(zeta_n), H=<p> in (Z/nZ)^*, K<-->H, t_k=Tr_{F/K}(zeta_n^k).
 * i_t[i]=k ==> iH=kH, i.e. t_i=t_k. We use t_k instead of t_i:
 * the number of k << the number of i. */
static GEN
get_i_t(long n, long p)
{
  long a;
  GEN v_t = const_vecsmall(n-1, 0);
  GEN i_t = cgetg(n, t_VECSMALL);  /* access i_t[a] with 1 <= a <= n-1 */
  for (a = 1; a < n; a++)
  {
    long b;
    while (a < n && v_t[a]) a++;
    if (a==n) break;
    b = a;
    do { v_t[b] = 1; i_t[b] = a; b = Fl_mul(b, p, n); } while (b != a);
  }
  return i_t;
}

/* get T_k(x) 1<= k <= d. d0*T_k(x) is in Z[x].
 * T_0(x)=T_n(x)=f.
 * Data = [H, GH, i_t, d0d1, kT, [n, d, f, n_T, mitk]] */
static GEN
get_vT(GEN Data, int NEW)
{
  pari_sp av = avma;
  GEN d0 = gmael(Data, 4, 1), kT = gel(Data, 5), N = gel(Data, 6);
  ulong k, n = N[1], n_T = N[4], mitk = N[5];
  GEN G1, M1, vT = const_vec(mitk, gen_0); /* vT[k]!=NULL ==> vT[k]=T_k */
  ulong n_k = 0, el, n_el, start, second;
  pari_timer ti;

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  if (!NEW) { gel(vT, 1) = pol_x(0); n_k++; }
  start = get_n_el(d0, &second);
  el = start_el_n(n);
  if (DEBUGLEVEL == 2) err_printf("get_vT: start=(%ld,%ld)\n",start,second);
  G1 = get_vG(vT, Data, start, &el, &M1);
  for (n_el = second;; n_el++)
  {
    GEN G2, M2, m, m2;
    G2 = get_vG(vT, Data, n_el, &el, &M2);
    m = mulii(M1, M2); m2 = shifti(m,-1);
    for (k = 1; k <= n_T; k++)
    {
      long j = kT[k];
      if (!isintzero(gel(vT,j))) continue;
      if (FpX_degsub(gel(G1,j), gel(G2,j), M2) < 0) /* G1=G2 (mod M2) */
      {
        gel(vT,j) = RgX_Rg_div(gel(G1,j), d0);
        n_k++;
        if (DEBUGLEVEL == 2)
          err_printf("G1:%ld, d0:%ld, M1:%ld, vT[%ld]:%ld words\n",
            ZX_size(gel(G1,j)), Q_size(d0), Q_size(M1), j, ZX_size(gel(vT,j)));
      }
      else
      {
        if (DEBUGLEVEL == 2)
          err_printf("G1:%ld, G2:%ld\n", ZX_size(gel(G1,j)),ZX_size(gel(G2,j)));
        gel(G1, j) = ZX_chinese_center(gel(G1,j),M1, gel(G2,j),M2, m,m2);
      }
    }
    if (n_k == n_T) break;
    M1 = m;
  }
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "get_vT");
  return gerepilecopy(av, vT);
}

/* return sorted kT={i_t[k] | 1<=k<=d}
 * {T_k(x) | k in kT} are all the different T_k(x) (1<=k<=d) */
static GEN
get_kT(GEN i_t, long d)
{ return vecsmall_uniq(vecsmall_shorten(i_t, d)); }

static GEN
get_kT_all(GEN GH, GEN i_t, long n, long d, long m)
{
  long i, j, k = 0;
  GEN x = const_vecsmall(d*m, 0);
  for (i = 1; i <= m; i++)
    for (j = 1; j <= d; j++) x[++k] = i_t[Fl_mul(GH[i], j, n)];
  return vecsmall_uniq(x);
}

static GEN
kT_from_kt_new(GEN gGH, GEN kt, GEN i_t, long n)
{
  GEN EL = gel(gGH, 1);
  long i, k = 0, lEL = lg(EL), lkt = lg(kt);
  GEN x = cgetg(lEL+lkt, t_VECSMALL);

  for (i = 1; i < lEL; i++) x[++k] = i_t[EL[i]];
  for (i = 2; i < lkt; i++) if (n%kt[i]==0) x[++k] = kt[i];
  setlg(x, k+1); return vecsmall_uniq(x);
}

static GEN
get_kTdiv(GEN kT, long n)
{
  long i, k = 0, l = lg(kT);
  GEN div = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) if (n%kT[i]==0) div[++k] = kT[i];
  setlg(div, k+1);
  return div;
}

/* T is separable over Zp but not separable over Fp.
 * receive all roots mod p^s and return all roots mod p^(s+1) */
static GEN
ZpX_roots_nonsep(GEN T, GEN R0, GEN p, GEN ps, GEN ps1)
{
  long i, j, n = 0, lr = lg(R0);
  GEN v = cgetg(lr, t_VEC), R;
  for (i = 1; i < lr; i++)
  {
    GEN z, f = ZX_unscale_div(ZX_translate(T, gel(R0, i)), ps);
    (void)ZX_pvalrem(f, p, &f);
    gel(v, i) = z = FpX_roots(f, p);
    n += lg(z)-1;
  }
  R = cgetg(n+1, t_VEC); n = 0;
  for (i = 1; i < lr; i++)
  {
    GEN z = gel(v, i);
    long lz = lg(z);
    for (j=1; j<lz; j++)
      gel(R, ++n) = Fp_add(gel(R0, i), mulii(gel(z, j), ps), ps1);
  }
  return ZV_sort_uniq_shallow(R);
}
static GEN
ZpX_roots_all(GEN T, GEN p, long f, long *ptrs)
{
  pari_sp av = avma;
  pari_timer ti;
  GEN v, ps, ps1;
  long s;

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  v = FpX_roots(T, p); /* FpX_roots returns sorted roots */
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpX_roots, deg=%ld", degpol(T));
  ps = NULL; ps1 = p;
  for (s = 1; lg(v) != f+1; s++)
  {
    ps = ps1; ps1 = mulii(ps1, p); /* p^s, p^(s+1) */
    v = ZpX_roots_nonsep(T, v, p, ps, ps1);
    if (gc_needed(av, 1)) gerepileall(av, 3, &v, &ps, &ps1);
  }
  *ptrs = s; return v;
}
/* x : roots of T in Zp. r < n.
 * receive vec of x mod p^r and return vec of x mod p^n.
 * under the assumtion lg(v)-1==degpol(T), x is uniquely determined by
 * x mod p^r. Namely, x mod p^n is unique. */
static GEN
ZX_Zp_liftroots(GEN T, GEN v, GEN p, long r, long n)
{
  long i, l = lg(v);
  GEN R = cgetg(l, t_VEC);
  GEN pr = powiu(p, r), pn = powiu(p, n);
  pari_timer ti;

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  for (i=1; i<l; i++)
  {
    GEN x = gel(v, i), y, z;
    GEN f = ZX_unscale_div(ZX_translate(T, x), pr);
    (void)ZX_pvalrem(f, p, &f);
    y = FpX_roots(f, p);
    z = (n==r+1) ? y : ZX_Zp_root(f, gel(y, 1), p, n-r);
    if (lg(y)!=2 || lg(z)!=2)
      pari_err_BUG("ZX_Zp_liftroots, roots are not separable");
    gel(R, i) = Fp_add(x, mulii(gel(z, 1), pr), pn);
  }
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "ZX_Zp_liftroots");
  return R;
}

static GEN
set_R(GEN T, GEN F, GEN Rs, GEN p, long nf, long r, long s, long u)
{
  long i;
  GEN x, pr = powiu(p, r), prs = powiu(p, r+s), R = cgetg(1+nf, t_VEC), Rrs;
  Rrs = r ? ZX_Zp_liftroots(T, Rs, p, s, r+s) : Rs;
  x = gel(Rrs, 1);
  for (i = 1; i <= nf; i++)
  {
    x = FpX_eval(F, x, prs);
    if (r) x = gel(Rrs, ZV_search(Rs, diviiexact(x, pr)));
    gel(R, i) = x;  /* R[i]=t_1^(g^i), t_1=Rrs[1], mod p^(r+s) */
  }
  if (r+s < u) R = ZX_Zp_liftroots(T, R, p, r+s, u);
  else if (r+s > u) R = FpV_red(R, powiu(p, u));
  return R;  /* mod p^u, accuracy for pol_newton */
}

/* Preparation for factoring polcyclo(el^e) mod p
 * f_0(x) : irred factor of polcyclo(el^e0) mod p
 * f_1(x)=f_0(x^(el^e1)) : irred factor of polcyclo(el^e) mod p
 *
 * p^d0 = 1 (mod el^e0), p^d = 1 (mod el^e)
 *
 * case el=2, 2^s || (p-1), s>=2
 * d=1 (1 <= e <= s), d=2^(e-s) (s < e)
 * e0=e, e1=0 if e <= s
 * e0=s, e1=e-s if s < e
 * d0=1
 *
 * case el=2, 2^s || (p+1), s>=2
 * d=1 (e == 1), d=2 (2 <= e <= s), d=2^(e-s) (s < e)
 * e0=e, e1=0 if e <= s+1
 * e0=s+1, e1=e-s-1 if s+1 < e
 * d0=1 if e==1,  d0=2 if e>1
 *
 * case el>2, el^s || (p^d0-1)
 * d=d0 (1 <= e <= s), d=d0*el^(e-s) (s < e)
 * e0=e, e1=0 if e <= s
 * e0=s, e1=e-s if s < e
 * d0 is min d s.t. p^d=1 (mod el)
 *
 * We do not need d. So d is not returned. */
static GEN
set_e0_e1(ulong el, ulong e, GEN p)
{
  ulong s, d0, e0, e1, f0, n, phin, g, up = itou_or_0(p);

  if (el == 2)
  {
    ulong pmod4 = up ? up&3 : mod4(p);
    if (pmod4 == 3)  /* p+1 = 0 (mod 4) */
    {
      s = up ? vals(up+1) : vali(addiu(p, 1));
      if (e <= s+1) { e0 = e; e1 = 0;}
      else { e0 = s+1; e1= e-s-1;}
      d0 = e == 1? 1: 2;
    }
    else  /* p-1 = 0 (mod 4) */
    {
      s = up ? vals(up-1) : vali(subiu(p, 1));
      if (e <= s) { e0 = e; e1 = 0; }
      else { e0 = s; e1 = e-s; }
      d0 = 1;
    }
    phin = 1L<<(e0-1);
  }
  else  /* el is odd */
  {
    ulong pmodel = up ? up%el : umodiu(p, el);
    d0 = Fl_order(pmodel, el-1, el);
    s = Z_lval(subiu(powiu(p, d0), 1), el);
    if (e <= s) { e0 = e; e1 = 0; } else { e0 = s; e1 = e-s; }
    phin = (el-1)*upowuu(el, e0-1);
  }
  n = upowuu(el, e0); f0 = phin/d0;
  g = (el==2) ? 1 : pgener_Zl(el);
  return mkvecsmalln(7, n, e0, e1, phin, g, d0, f0);
}

/* return 1 if newton is fast, return 0 if gen is fast */
static int
use_newton(long d, long f)
{
  if (2*d <= f) return 0;
  else if (f <= d) return 1;
  else if (d <= 50) return 0;
  else if (f <= 60) return 1;
  else if (d <= 90) return 0;
  else if (f <= 150) return 1;
  else if (d <= 150) return 0;
  else if (f <= 200) return 1;
  else if (200*d <= f*f) return 0;
  else return 1;
}

/* return 1 if newton_general is fast, return 0 otherwise. Assume f > 40 */
static int
use_newton_general(long d, long f, long maxdeg)
{
  if (maxdeg < 20) return 0;
  else if (f <= 50) return 1;
  else if (maxdeg < 30) return 0;
  else if (f <= 60) return 1;
  else if (maxdeg < 40) return 0;
  else if (f <= 70) return 1;
  else if (maxdeg < 50) return 0;
  else if (f <= 80) return 1;
  else if (d < 200) return 0;
  else if (f <= 100) return 1;
  else if (d < 300) return 0;
  else if (f <= 120) return 1;
  else if (6*maxdeg < f*f) return 0;
  else return 1;
}

static int
use_general(long d, long maxdeg)
{
  if (d <= 50) return 1;
  else if (maxdeg <= 3*d) return 0;
  else if (d <= 200) return 1;
  else if (maxdeg <= 10*d) return 0;
  else if (d <= 500) return 1;
  else if (maxdeg <= 20*d) return 0;
  else if (d <= 1000) return 1;
  else return 0;
}

static void
update_dfm(long *pd, long *pf, long *pm, long di, long fi)
{
  long c = ugcd(*pd,di), d1 = *pd * di, f1 = *pf * fi;
  *pd = d1 / c; *pf = c * f1; *pm += d1 * d1 * f1;
  if (DEBUGLEVEL == 4) err_printf("(%ld,%ld), ",d1,f1);
}
/* assume ord(p mod f) > 1 */
static ulong
set_action(GEN fn, GEN p, long d, long f)
{
  GEN EL = gel(fn, 1), E = gel(fn, 2);
  long i, d0, f0, m0, m1, maxdeg, max, l = lg(EL);
  ulong action = 0, up = itou_or_0(p);
  GEN D = cgetg(l, t_VECSMALL), F = cgetg(l, t_VECSMALL);

  d += 10*(lgefint(p)-3);
  if (l == 2)
  { /* n is a prime power */
    action |= (EL[1]==2 || !use_newton(d, f))? GENERAL: NEWTON_POWER;
    return action;
  }
  if (f <= 2) action |= NEWTON_GENERAL;
  else if (d <= 10) action |= GENERAL;
  else if (f <= 10) action |= NEWTON_GENERAL;
  else if (d <= 20) action |= GENERAL;
  else if (f <= 20) action |= NEWTON_GENERAL_NEW;
  else if (d <= 30) action |= GENERAL;
  else if (f <= 30) action |= NEWTON_GENERAL_NEW;
  else if (d <= 40) action |= GENERAL;
  else if (f <= 40) action |= NEWTON_GENERAL_NEW;
  if (action) return action;
  /* can assume that d > 40 and f > 40 */

  maxdeg = max = 1;
  for (i = 1; i < l; i++)
  {
    long x, el = EL[i], e = E[i];
    long q = upowuu(el, e-1), ni = q * el, phini = ni - q;
    long di = Fl_order(umodiu(p, ni), phini, ni);
    D[i] = di; F[i] = phini / di;
    x = ugcd(max, di); max = max * (di / x); /* = lcm(d1,..di) */
    if (x > 1) maxdeg = max*x;
    if (DEBUGLEVEL == 4) err_printf("(%ld,%ld), ", D[i], F[i]);
  }
  if (maxdeg == 1) return action;
  if (up != 2)
  {
    if (use_newton_general(d, f, maxdeg))
    { /* does not decompose n */
      action |= (20 < d)? NEWTON_GENERAL_NEW: NEWTON_GENERAL;
      return action;
    }
    if (use_general(d, maxdeg)) action |= GENERAL;
  }
  if (l < 4) return action; /* n has two factors */

  d0 = f0 = 1; m0 = 0;
  for (i = 1; i < l; i++) update_dfm(&d0, &f0, &m0, D[i], F[i]);
  if (DEBUGLEVEL == 4) err_printf("(d0,f0)=(%ld,%ld)\n",d0,f0);
  d0 = f0 = 1; m1 = 0;
  for (i = l-1; i >= 1; i--) update_dfm(&d0, &f0, &m1, D[i], D[i]);
  if (DEBUGLEVEL == 4) err_printf("(d0,f0)=(%ld,%ld)\n",d0,f0);
  if (DEBUGLEVEL == 4) err_printf("(m0,m1)=(%lu,%lu) %ld\n",m0,m1,m0<=m1);
  if (m0 <= m1) action |= ASCENT;
  return action;
}

static GEN
FpX_Newton_perm(long d, GEN R, GEN v, GEN pu, GEN p)
{
  GEN S = cgetg(d+2, t_VEC);
  long k;
  gel(S,1) = utoi(d); for (k = 1; k <= d; k++) gel(S, k+1) = gel(R, v[k]);
  return FpX_red(FpX_fromNewton(RgV_to_RgX(S, 0), pu), p);
}
static GEN
Flx_Newton_perm(long d, GEN R, GEN v, ulong pu, ulong p)
{
  GEN S = cgetg(d+2, t_VEC);
  long k;
  S[1] = d; for (k = 1; k <= d; k++) gel(S, k+1) = gel(R, v[k]);
  return Flx_red(Flx_fromNewton(Flv_to_Flx(S, 0), pu), p);
}

/*  Data = [H, GH, i_t, d0, kT, [n, d, f, n_T, mitk]]
    N2 = [p, pr, pu, pru] */
static GEN
FpX_pol_newton_general(GEN Data, GEN N2, GEN vT, GEN x)
{
  GEN i_t = gel(Data, 3), kT = gel(Data, 5), N = gel(Data, 6);
  long k, d = N[2], n_T = N[4], mitk = N[5];
  GEN p = gel(N2,1), pr = gel(N2,2), pu = gel(N2,3), pru = gel(N2,4);
  GEN R = cgetg(1+mitk, t_VEC);

  for (k = 1; k <= n_T; k++)
    gel(R, kT[k]) = diviiexact(FpX_eval(gel(vT, kT[k]), x, pru), pr);
  return FpX_Newton_perm(d, R, i_t, pu, p);
}

/* n is any integer prime to p, but must be equal to the conductor
 * of the splitting field of p in Q(zeta_n).
 * GH=G/H={g_1, ... ,g_f}
 * Data = [GHgen, GH, fn, p, [n, d, f, m]] */
static GEN
FpX_factcyclo_newton_general(GEN Data)
{
  GEN GH = gel(Data, 2), fn = gel(Data, 3), p = gel(Data, 4);
  long n = mael(Data, 5, 1), d = mael(Data, 5, 2), f = mael(Data, 5, 3);
  long m = mael(Data, 5, 4), pmodn = umodiu(p, n);
  long i, k, n_T, mitk, r, s = 0, u = 1;
  GEN vT, kT, H, i_t, T, d0d1, Data2, Data3, R, Rs, pr, pu, pru;
  pari_timer ti;

  if (m != 1) m = f;
  for (pu = p; cmpiu(pu,d) <= 0; u++) pu = mulii(pu, p);  /* d<pu, pu=p^n */

  H = Fl_powers(pmodn, d-1, n); /* H=<p> */
  i_t = get_i_t(n, pmodn); /* i_t[1+i]=k ==> t_i=t_k */
  kT = get_kT(i_t, d); n_T = lg(kT)-1; mitk = kT[n_T];
  if (DEBUGLEVEL == 2) err_printf("kT=%Ps %ld elements\n",kT,n_T);
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  T = galoissubcyclo(utoi(n), utoi(pmodn), 0, 0);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "galoissubcyclo");
  d0d1 = get_d0_d1(T, gel(fn,1)); /* d0*T_k(x) is in Z[x] */
  Data2 = mkvecn(6, H, GH, i_t, d0d1, kT, mkvecsmalln(5, n, d, f, n_T, mitk));
  vT = get_vT(Data2, 0);
  if (DEBUGLEVEL == 4) err_printf("vT=%Ps\n",vT);
  r = QXV_den_pval(vT, kT, p);
  Rs = ZpX_roots_all(T, p, f, &s);
  if (DEBUGLEVEL >= 2) err_printf("(u,s,r)=(%ld,%ld,%ld)\n",u,s,r);
  if (r+u < s) pari_err_BUG("FpX_factcyclo_newton_general (T is not separable mod p^(r+u))");
  /* R and vT are mod p^(r+u) */
  R = (r+u==s) ? Rs : ZX_Zp_liftroots(T, Rs, p, s, r+u);
  pr = powiu(p, r); pru = powiu(p, r+u); /* Usually, r=0, s=1, pr=1, pru=p */
  for (k = 1; k<=n_T; k++)
  {
    long itk = kT[k];
    GEN z = r? RgX_Rg_mul(gel(vT, itk), pr): gel(vT, itk);
    gel(vT, itk) = RgX_to_FpX(z, pru);
  }
  Data3 = mkvec4(p, pr, pu, pru);
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  for (i=1; i<=m; i++)
    gel(R,i) = FpX_pol_newton_general(Data2, Data3, vT, gel(R,i));
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpX_pol_newton_general");
  return R;
}

/* Data = [vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
          [n, r, s, n_t,mitk], div] */
static void
Fp_next_gen3(long x, long i, GEN v_t_p, GEN t, GEN Data)
{
  GEN vT = gel(Data, 1), gGH = gel(Data, 2), Rs = gel(Data, 3);
  GEN Rrs = gel(Data, 4), i_t = gel(Data, 5);
  GEN pu = gel(Data, 8), pr = gel(Data, 9), prs = gel(Data, 10);
  GEN EL = gel(gGH, 1), E = gel(gGH, 2), div = gel(Data, 12);
  long n = mael(Data, 11, 1), r = mael(Data, 11, 2), mitk = mael(Data, 11, 5);
  long j, k, l = lg(EL), ld = lg(div);
  if (i >= l) return;
  for (j = 0; j < E[i]; j++)
  {
    if (j > 0)
    {
      long itx = i_t[x];
      t = FpX_eval(gel(vT, i_t[EL[i]]), t, prs); /* mod p^(r+s) */
      if (r) t = gel(Rrs, ZV_search(Rs, diviiexact(t, pr))); /* mod p^(r+s) */
      if (itx <= mitk) gel(v_t_p, itx) = Fp_red(t, pu); /* mod p^u */
      for (k = 1; k<ld; k++)
      {
        ulong y = Fl_mul(x, div[k], n);
        long ity = i_t[y];
        GEN v;
        if (ity > mitk || !isintzero(gel(v_t_p, ity))) continue;
        v = FpX_eval(gel(vT, i_t[div[k]]), t, prs); /* mod p^(r+s) */
        if (r) v = diviiexact(v, pr); /* mod p^s */
        gel(v_t_p, ity) = Fp_red(v, pu);
      }
    }
    Fp_next_gen3(x, i+1, v_t_p, t, Data);
    x = Fl_mul(x, EL[i], n);
  }
}

/* Data = [vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
          [n, r, s, n_t, mitk], div] */
static GEN
Fp_mk_v_t_p3(GEN Data, long i)
{ /* v_t_p[k] != gen_0 => v_t_p[k] = t_k mod p^u */
  GEN Rs = gel(Data, 3), Rrs = gel(Data, 4);
  GEN pu = gel(Data, 8), pr = gel(Data, 9), prs = gel(Data, 10);
  GEN vT = gel(Data, 1), i_t = gel(Data, 5), div = gel(Data, 12);
  long k, r = mael(Data, 11, 2), mitk = mael(Data, 11, 5), ld = lg(div);
  GEN v_t_p = const_vec(mitk, gen_0);

  gel(v_t_p, 1) = Fp_red(gel(Rs, i), pu); /* mod p^u, guaranteed u<=s */
  Fp_next_gen3(1, 1, v_t_p, gel(Rrs, i), Data);
  for (k = 1; k<ld; k++)
  {
    ulong itk = i_t[div[k]];
    GEN x = FpX_eval(gel(vT, itk), gel(Rrs, i), prs);
    if (r) x = diviiexact(x, pr); /* mod p^s */
    gel(v_t_p, itk) = Fp_red(x, pu);
  }
  return v_t_p;
}

/* Data = [vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
          [n, r, s, n_t,mitk], div] */
static void
Fl_next_gen3(long x, long i, GEN v_t_p, ulong t, GEN Data)
{
  GEN vT = gel(Data, 1), gGH = gel(Data, 2), Rs = gel(Data, 3);
  GEN Rrs = gel(Data, 4), i_t = gel(Data, 5);
  long pu = mael(Data, 8, 2), pr = mael(Data, 9, 2), prs = mael(Data, 10, 2);
  GEN EL = gel(gGH, 1), E = gel(gGH, 2), div = gel(Data, 12);
  long n = mael(Data, 11, 1), r = mael(Data, 11, 2), mitk = mael(Data, 11, 5);
  long j, k, l = lg(EL), ld = lg(div);
  if (i >= l) return;
  for (j = 0; j < E[i]; j++)
  {
    if (j > 0)
    {
      long itx = i_t[x];
      t = Flx_eval(gel(vT, i_t[EL[i]]), t, prs); /* mod p^(r+s) */
      if (r) t = Rrs[zv_search(Rs, t/pr)]; /* mod p^(r+s) */
      if (itx <= mitk) v_t_p[itx] = t%pu; /* mod p^u */
      for (k = 1; k < ld; k++)
      {
        ulong y = Fl_mul(x, div[k], n), v;
        long ity = i_t[y];
        if (ity > mitk || v_t_p[ity]) continue;
        v = Flx_eval(gel(vT, i_t[div[k]]), t, prs); /* mod p^(r+s) */
        if (r) v /= pr; /* mod p^s */
        v_t_p[ity] = v%pu;
      }
    }
    Fl_next_gen3(x, i+1, v_t_p, t, Data);
    x = Fl_mul(x, EL[i], n);
  }
}

/* Data = [vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
          [n, r, s, n_t, mitk], div] */
static GEN
Fl_mk_v_t_p3(GEN Data, long i)
{ /* v_t_p[k] != 0 => v_t_p[k] = t_k mod p^u */
  GEN Rs = gel(Data, 3), Rrs = gel(Data, 4);
  ulong pu = mael(Data, 8, 2), pr = mael(Data, 9, 2), prs = mael(Data, 10, 2);
  GEN vT = gel(Data, 1), i_t = gel(Data, 5), div = gel(Data, 12);
  long k, r = mael(Data, 11, 2), mitk = mael(Data, 11, 5), ld = lg(div);
  GEN v_t_p = const_vecsmall(mitk, 0);

  v_t_p[1] = Rs[i] % pu; /* mod p^u, guaranteed u<=s */
  Fl_next_gen3(1, 1, v_t_p, Rrs[i], Data);
  for (k = 1; k < ld; k++)
  {
    ulong itk = i_t[div[k]], x = Flx_eval(gel(vT, itk), Rrs[i], prs);
    if (r) x /= pr; /* mod p^s */
    v_t_p[itk] = x % pu;
  }
  return v_t_p;
}

/* Data = [GHgen, GH, fn, p, [n, d, f, m]] */
static GEN
newton_general_new_pre3(GEN Data)
{
  GEN gGH = gel(Data, 1), GH = gel(Data, 2), fn = gel(Data, 3);
  GEN p = gel(Data, 4);
  long n = mael(Data, 5, 1), d = mael(Data, 5, 2), f = mael(Data, 5, 3);
  long pmodn = umodiu(p, n);
  long k, n_t, n_T, mitk, miTk, r, s = 0, u = 1;
  GEN vT, kt, kT, H, i_t, T, d0d1, Data2, Rs, Rrs, kTdiv;
  GEN pr, pu, prs;
  pari_timer ti;

  for (pu = p; cmpiu(pu,d)<=0; u++) pu = mulii(pu, p);  /* d<pu, pu=p^u */

  H = Fl_powers(pmodn, d-1, n); /* H=<p> */
  i_t = get_i_t(n, pmodn); /* i_t[1+i]=k ==> t_i=t_k */
  kt = get_kT_all(GH, i_t, n, d, 1); n_t = lg(kt)-1; mitk = kt[n_t];
  kT = kT_from_kt_new(gGH, kt, i_t, n); n_T = lg(kT)-1; miTk = kT[n_T];
  kTdiv = get_kTdiv(kT, n);
  if (DEBUGLEVEL == 2)
    err_printf("kt=%Ps %ld elements\nkT=%Ps %ld elements\n",kt,n_t,kT,n_T);
  if (DEBUGLEVEL == 2)
    err_printf("kTdiv=%Ps\n",kTdiv);

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  T = galoissubcyclo(utoi(n), utoi(pmodn), 0, 0);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "galoissubcyclo");
  d0d1 = get_d0_d1(T, gel(fn,1)); /* d0*T_k(x) is in Z[x] */
  Data2 = mkvecn(6, H, GH, i_t, d0d1, kT, mkvecsmalln(5, n, d, f, n_T, miTk));
  vT = get_vT(Data2, 1);
  if (DEBUGLEVEL == 4) err_printf("vT=%Ps\n",vT);
  r = QXV_den_pval(vT, kT, p);
  Rs = ZpX_roots_all(T, p, f, &s);
  if (DEBUGLEVEL >= 2) err_printf("(u,s,r)=(%ld,%ld,%ld)\n",u,s,r);
  if (s < u)
  {
    Rs = ZV_sort_shallow(ZX_Zp_liftroots(T, Rs, p, s, u));
    s = u;
  }
  /* Rs : mod p^s, Rrs : mod p^(r+s), vT : mod p^(r+s) */
  Rrs = r ? ZX_Zp_liftroots(T, Rs, p, s, r+s) : Rs;
  pr = powiu(p, r); prs = powiu(p, r+s); /* Usually, r=0, s=1, pr=1, pru=p */
  if (lgefint(prs)>3)  /* ULONG_MAX < p^(r+s), usually occur */
  {
    for (k = 1; k <= n_T; k++)
    {
      long itk = kT[k];
      GEN z = r? RgX_Rg_mul(gel(vT, itk), pr): gel(vT, itk);
      gel(vT, itk) = RgX_to_FpX(z, prs);
    }
  }
  else  /* p^(r+s) < ULONG_MAX, frequently occur */
  {
    ulong upr = itou(pr), uprs = itou(prs);
    for (k = 1; k <= n_T; k++)
    {
      long itk = kT[k];
      GEN z = r? RgX_muls(gel(vT, itk), upr): gel(vT, itk);
      gel(vT, itk) = RgX_to_Flx(z, uprs);
    }
    Rs = ZV_to_nv(Rs); Rrs = ZV_to_nv(Rrs);
  }
  return mkvecn(12, vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
                mkvecsmalln(6, n, r, s, n_t, mitk, d), kTdiv);
}

/* Data=[vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
 *      [n, r, s, n_t, mitk, d], div] */
static GEN
FpX_pol_newton_general_new3(GEN Data, long k)
{
  GEN i_t = gel(Data, 5), p = gel(Data, 7), pu = gel(Data, 8);
  long d = mael(Data, 11, 6);
  GEN v_t_p = Fp_mk_v_t_p3(Data, k);
  if (DEBUGLEVEL == 3) err_printf("v_t_p=%Ps\n",v_t_p);
  return FpX_Newton_perm(d, v_t_p, i_t, pu, p);
}

/* Data = [GHgen, GH, fn, p, [n, d, f, m]] */
static GEN
FpX_factcyclo_newton_general_new3(GEN Data)
{
  long i, f = mael(Data, 5, 3), m = mael(Data, 5, 4);
  GEN Data2, pols;
  pari_timer ti;

  if (m != 1) m = f;
  pols = cgetg(1+m, t_VEC);
  Data2 = newton_general_new_pre3(Data);
  /* Data2 = [vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
              [n, r, s, n_t, mitk, d], div] */
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  for (i = 1; i <= m; i++)
    gel(pols, i) = FpX_pol_newton_general_new3(Data2, i);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpX_pol_newton_general_new3");
  return pols;
}

/* return normalized z(-x) */
static GEN
FpX_1st_lift_2(GEN z, GEN p)
{
  long i, r = lg(z);
  GEN x = cgetg(r, t_POL);
  x[1] = evalsigne(1) | evalvarn(0);
  if (odd(r))
    for (i = 2; i < r; i++) gel(x,i) = odd(i)? Fp_neg(gel(z,i), p): gel(z,i);
  else
    for (i = 2; i < r; i++) gel(x,i) = odd(i)? gel(z,i): Fp_neg(gel(z,i), p);
  return x;
}

static GEN
FpX_1st_lift(GEN z, GEN p, ulong e, ulong el, GEN vP)
{
   GEN z2, z3, P = gel(vP, e);
   if (!gel(vP, e)) P = gel(vP, e) = FpX_polcyclo(e, p);
   z2 = RgX_inflate(z, el);
   z3 = FpX_normalize(FpX_gcd(P, z2, p), p);
   return FpX_div(z2, z3, p);
}

static GEN
FpX_lift(GEN z, GEN p, ulong e, ulong el, ulong r, ulong s, GEN vP)
{
  if (s == 0)
  {
    z = (el==2) ? FpX_1st_lift_2(z, p) : FpX_1st_lift(z, p, e, el, vP);
    if (r >= 2) z = RgX_inflate(z, upowuu(el, r-1));
  }
  else
    z = RgX_inflate(z, upowuu(el, r-s));
  return z;
}

/* e is the conductor of the splitting field of p in Q(zeta_n) */
static GEN
FpX_conductor_lift(GEN z, GEN p, GEN fn, ulong e, GEN vP)
{
  GEN EL = gel(fn, 1), En = gel(fn, 2);
  long i, r = lg(EL), new_e = e;

  for (i = 1; i < r; i++)
  {
    long el = EL[i], en = En[i], ee = u_lval(e, el);
    if (ee < en)
    {
      z = FpX_lift(z, p, new_e, el, en, ee, vP);
      new_e *= upowuu(el, en-ee);
    }
  }
  return z;
}

/* R0 is mod p^u, d < p^u */
static GEN
FpX_pol_newton(long j, GEN R0, GEN E, GEN D3, long d, long f, GEN p)
{
  long i, u = D3[3];
  GEN R = cgetg(1+f, t_VEC);
  for (i = 1; i <= f; i++) gel(R, i) = gel(R0, 1+(i+j)%f);
  return FpX_Newton_perm(d, R, E, powiu(p, u), p);
}

/* Data = [T, F, Rs, [d, nf, g, r, s, u]], nf>1 */
static GEN
FpX_factcyclo_newton_pre(GEN Data, GEN N, GEN p, ulong m)
{
  GEN T = gel(Data, 1), F = gel(Data, 2), Rs = gel(Data, 3), D4 = gel(Data, 4);
  long d = D4[1], nf = D4[2], g = D4[3], r = D4[4], s = D4[5], u = D4[6];
  GEN pols, E, R, Data3;
  ulong i, n = N[1], pmodn = umodiu(p, n);
  pari_timer ti;

  if (m!=1) m = nf;
  pols = cgetg(1+m, t_VEC);
  E = set_E(pmodn, n, d, nf, g);
  R = set_R(T, F, Rs, p, nf, r, s, u);
  Data3 = mkvecsmall3(r, s, u);
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  for (i = 1; i <= m; i++) gel(pols, i) = FpX_pol_newton(i, R,E,Data3, d,nf,p);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpX_pol_newton");
  return pols;
}

/* gcd(n1,n2)=1, e = gcd(d1,d2).
 * phi(n1) = d1*nf1, phi(n2) = d2*nf2.
 * d = lcm(d1,d2) = d1*d2/e, nf = phi(n1)*phi(n2)/d = nf1*nf2*e. */
static GEN
FpX_factcyclo_lift(long n1, GEN v1, long n2, GEN v2, GEN p, long m)
{
  pari_sp av = avma;
  long d1 = degpol(gel(v1, 1)), lv1 = lg(v1), nf1 = eulerphiu(n1)/d1;
  long d2 = degpol(gel(v2, 1)), lv2 = lg(v2), nf2 = eulerphiu(n2)/d2;
  long e = ugcd(d1, d2), nf = nf1*nf2*e, need_factor = (e!=1);
  long i, j, k = 0;
  GEN z = NULL, v;
  pari_timer ti;

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  if (m != 1) m = nf;
  v = cgetg(1+m, t_VEC);
  if (m == 1)
  {
    GEN x1 = gel(v1, 1), x2 = gel(v2, 1);
    z = FpX_gcd(RgX_inflate(x1, n2), RgX_inflate(x2, n1), p);
    z = FpX_normalize(z, p);  /* FpX_gcd sometimes returns non-monic */
    gel(v, 1) = (!need_factor)? z : gmael(FpX_factor(z, p), 1, 1);
  }
  else
  {
    for (i = 1; i < lv1; i++)
      for (j = 1; j < lv2; j++)
      {
        GEN x1 = gel(v1, i), x2 = gel(v2, j);
        z = FpX_gcd(RgX_inflate(x1, n2), RgX_inflate(x2, n1), p);
        z = FpX_normalize(z, p);  /* needed */
        if (!need_factor) gel(v, ++k) = z;
        else
        {
          GEN z1 = gel(FpX_factor(z, p), 1);
          long i, l = lg(z1);
          for (i = 1; i < l; i++) gel(v, ++k) = gel(z1, i);
        }
      }
  }
  if (DEBUGLEVEL >= 6)
    timer_printf(&ti, "FpX_factcyclo_lift (%ld,%ld)*(%ld,%ld)-->(%ld,%ld)-->(%ld,%ld)",
        d1, nf1, d2, nf2, degpol(z), nf1*nf2, d1*d2/e, nf);
  return gerepilecopy(av, v);
}

/* n is any integer prime to p; d>1 and f>1 */
static GEN
FpX_factcyclo_gen(GEN GH, ulong n, GEN p, long m)
{
  GEN fn = factoru(n), fa = zm_to_ZM(fn);
  GEN A, T, X, pd_n, v, pols;
  long i, j, pmodn = umodiu(p, n), phin = eulerphiu_fact(fn);
  long d = Fl_order(pmodn, phin, n), f = phin/d;
  pari_timer ti;

  if (m != 1) m = f;
  if (m != 1 && GH==NULL) /* FpX_factcyclo_fact is used */
  {
    GEN H = znstar_generate(n, mkvecsmall(pmodn));
    GH = znstar_cosets(n, phin, H); /* representatives of G/H */
  }

  pols = cgetg(1+m, t_VEC);
  v = cgetg(1+d, t_VEC);
  pd_n = diviuexact(subis(powiu(p, d), 1), n); /* (p^d-1)/n */
  T = init_Fq(p, d, 1);
  A = pol_x(1);  /* A is a generator of F_q, q=p^d */
  if (d == 1) A = FpX_rem(A, T, p);
  random_FpX(degpol(T), varn(T), p); /* skip 1st one */
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  do X = FpXQ_pow(random_FpX(degpol(T), varn(T), p), pd_n, T, p);
  while(lg(X)<3 || equaliu(FpXQ_order(X, fa, T, p), n)==0); /* find zeta_n */
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "find X");

  if (m == 1)
  {
    for (j = 1; j <= d; j++)
    {
      gel(v, j) = pol_x(0);
      gmael(v, j, 2) = FpX_neg(X, p);
      if (j < d) X = FpXQ_pow(X, p, T, p);
    }
    gel(pols, 1) = ZXX_evalx0(FpXQXV_prod(v, T, p));
  }
  else
  {
    GEN vX, vp;
    if (DEBUGLEVEL >= 6) timer_start(&ti);
    vX = FpXQ_powers(X, n, T, p)+1;
    vp = Fl_powers(pmodn, d, n);
    for (i = 1; i <= m; i++)
    {
      for (j = 1; j <= d; j++)
      {
        gel(v, j) = pol_x(0);
        gmael(v, j, 2) = FpX_neg(gel(vX, Fl_mul(GH[i], vp[j], n)), p);
      }
      gel(pols, i) = ZXX_evalx0(FpXQXV_prod(v, T, p));
    }
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpXQXV_prod");
  }
  return pols;
}

static GEN Flx_factcyclo_newton_pre(GEN Data, GEN N, ulong p, ulong m);
/* n is an odd prime power prime to p and equal to the conductor
 * of the splitting field of p in Q(zeta_n). d>1 and nf>1 */
static GEN
FpX_factcyclo_newton_power(GEN N, GEN p, ulong m, int small)
{
  GEN Rs, H, T, F, pr, prs, pu, Data;
  long n = N[1], el = N[2], phin = N[4], g = N[5];
  long pmodn = umodiu(p, n), pmodel = umodiu(p, el);
  long d = Fl_order(pmodel, el-1, el), nf = phin/d;
  long r, s = 0, u = 1;

  if (m != 1) m = nf;
  for (pu = p; cmpiu(pu,d) <= 0; u++) pu = mulii(pu,p);  /* d<p^u, pu=p^u */
  H = Fl_powers(pmodn, d, n);
  T = galoissubcyclo(utoi(n), utoi(pmodn), 0, 0);
  F = gausspol(T, H, N, d, nf, g);
  r = QX_den_pval(F, p);
  Rs = ZpX_roots_all(T, p, nf, &s);
  if (DEBUGLEVEL >= 2) err_printf("(u,s,r)=(%ld,%ld,%ld)\n",u,s,r);
  pr = powiu(p, r); prs = powiu(p, r+s); /* Usually, r=0, s=1, pr=1, prs=p */
  F = r ? RgX_to_FpX(RgX_Rg_mul(F, pr), prs) : RgX_to_FpX(F, prs);
  Data = mkvec4(T, F, Rs, mkvecsmalln(6, d, nf, g, r, s, u));
  if (small && lgefint(pu) == 3)
    return Flx_factcyclo_newton_pre(Data, N, uel(p,2), m);
  else
    return FpX_factcyclo_newton_pre(Data, N, p, m);
}

static GEN
FpX_split(ulong n, GEN p, ulong m)
{
  ulong i, j;
  GEN v, C, vx, z = rootsof1u_Fp(n, p);
  if (m == 1) return mkvec(deg1pol_shallow(gen_1, Fp_neg(z,p), 0));
  v = cgetg(m+1, t_VEC); C = coprimes_zv(n); vx = Fp_powers(z, n-1, p);
  for (i = j = 1; i <= n; i++)
    if (C[i]) gel(v, j++) = deg1pol_shallow(gen_1, Fp_neg(gel(vx,i+1), p), 0);
  return gen_sort(v, (void*)cmpii, &gen_cmp_RgX);
}

/* n is a prime power prime to p. n is not necessarily equal to the
 * conductor of the splitting field of p in Q(zeta_n). */
static GEN
FpX_factcyclo_prime_power_i(long el, long e, GEN p, long m)
{
  GEN z = set_e0_e1(el, e, p), v;
  long n = z[1], e0 = z[2], e1 = z[3], phin = z[4], g = z[5];
  long d = z[6], f = z[7]; /* d and f for n=el^e0 */

  if (f == 1) v = mkvec(FpX_polcyclo(n, p));
  else if (d == 1) v = FpX_split(n, p, (m==1)? 1: f);
  else if (el == 2) v = FpX_factcyclo_gen(NULL, n, p, m); /* d==2 in this case */
  else if (!use_newton(d,f)) v = FpX_factcyclo_gen(NULL, n, p, m);
  else
  {
    GEN N = mkvecsmall5(n, el, e0, phin, g);
    v = FpX_factcyclo_newton_power(N, p, m, 0);
  }
  if (e1)
  {
    long i, l = lg(v), ele1 = upowuu(el, e1);
    for (i = 1; i < l; i++) gel(v, i) = RgX_inflate(gel(v, i), ele1);
  }
  return v;
}
static GEN
FpX_factcyclo_prime_power(long el, long e, GEN p, long m)
{
  pari_sp av = avma;
  return gerepilecopy(av, FpX_factcyclo_prime_power_i(el, e, p, m));
}

static GEN
FpX_factcyclo_fact(GEN fn, GEN p, ulong m, long ascent)
{
  GEN EL = gel(fn, 1), E = gel(fn, 2), v1 = NULL;
  long i, l = lg(EL), n1 = 1;

  for (i = 1; i < l; i++)
  {
    long j = ascent? i: l-i, n2 = upowuu(EL[j], E[j]);
    GEN v2 = FpX_factcyclo_prime_power(EL[j], E[j], p, m);
    v1 = v1? FpX_factcyclo_lift(n1, v1, n2, v2, p, m): v2;
    n1 *= n2;
  }
  return v1;
}

static GEN
Flv_FlvV_factorback(GEN g, GEN x, ulong q)
{ pari_APPLY_ulong(Flv_factorback(g, gel(x,i), q)) }

/* return the structure and the generators of G/H. G=(Z/nZ)^, H=<p>.
 * For efficiency assume p mod n != 1 (trivial otherwise) */
static GEN
get_GH_gen(long n, long pmodn)
{
  GEN G = znstar0(utoipos(n), 1), cycG = znstar_get_cyc(G);
  GEN cycGH, gG, gGH, Ui, P;
  ulong expG;
  P = hnfmodid(mkmat(Zideallog(G, utoi(pmodn))), cycG);
  cycGH = ZV_to_nv(ZM_snf_group(P, NULL, &Ui));
  expG = itou(cyc_get_expo(cycG));
  gG = ZV_to_Flv(znstar_get_gen(G), n);
  gGH = Flv_FlvV_factorback(gG, ZM_to_Flm(Ui, expG), n);
  return mkvec2(gGH, cycGH);
}

/* 1st output */
static void
header(GEN fn, long n, long d, long f, GEN p)
{
  GEN EL = gel(fn, 1), E = gel(fn, 2);
  long i, l = lg(EL)-1;
  err_printf("n=%lu=", n);
  for (i = 1; i <= l; i++)
  {
    long el = EL[i], e = E[i];
    err_printf("%ld", el);
    if (e > 1) err_printf("^%ld", e);
    if (i < l) err_printf("*");
  }
  err_printf(", p=%Ps, phi(%lu)=%lu*%lu\n", p, n, d, f);
  err_printf("(n,d,f) : (%ld,%ld,%ld) --> ",n,d,f);
}

static ulong
FpX_factcyclo_just_conductor_init(GEN *pData, ulong n, GEN p, ulong m)
{
  GEN fn = factoru(n), GH = NULL, GHgen = NULL;
  long phin = eulerphiu_fact(fn), pmodn = umodiu(p, n);
  long d = Fl_order(pmodn, phin, n), f = phin/d; /* d > 1 */
  ulong action = set_action(fn, p, d, f);
  if (action & ~NEWTON_POWER)
  { /* needed for GENERAL* */
    GEN H = znstar_generate(n, mkvecsmall(pmodn));
    GH = znstar_cosets(n, phin, H); /* representatives of G/H */
    if (action & (NEWTON_GENERAL_NEW | NEWTON_GENERAL))
      GHgen = get_GH_gen(n, pmodn);  /* gen and order of G/H */
  }
  *pData = mkvec5(GHgen, GH, fn, p, mkvecsmall4(n, d, f, m));
  if (DEBUGLEVEL >= 1)
  {
    err_printf("(%ld,%ld,%ld)  action=%ld\n", n, d, f, action);
    if (GHgen)
    {
      GEN cycGH = gel(GHgen,2), gGH = gel(GHgen,1);
      err_printf("G(K/Q)=%Ps gen=%Ps\n", zv_to_ZV(cycGH), zv_to_ZV(gGH));
    }
  }
  return action;
}

static GEN
FpX_factcyclo_just_conductor(ulong n, GEN p, ulong m)
{
  GEN Data, fn;
  ulong action = FpX_factcyclo_just_conductor_init(&Data, n, p, m);
  fn = gel(Data,3);
  if (action & GENERAL)
    return FpX_factcyclo_gen(gel(Data,2), n, p, m);
  else if (action & NEWTON_POWER)
    return FpX_factcyclo_prime_power_i(ucoeff(fn,1,1), ucoeff(fn,1,2), p, m);
  else if (action & NEWTON_GENERAL)
    return FpX_factcyclo_newton_general(Data);
  else if (action & NEWTON_GENERAL_NEW)
    return FpX_factcyclo_newton_general_new3(Data);
  else
    return FpX_factcyclo_fact(fn, p, m, action & ASCENT);
}

static GEN
FpX_factcyclo_i(ulong n, GEN p, long fl)
{
  GEN fn = factoru(n), z;
  long phin = eulerphiu_fact(fn), pmodn = umodiu(p, n);
  ulong d = Fl_order(pmodn, phin, n), f = phin/d, fK;

  if (DEBUGLEVEL >= 1) header(fn, n, d, f, p);
  if (f == 1) { z = FpX_polcyclo(n, p); return fl? z: mkvec(z); }
  else if (d == 1) /* p=1 (mod n), zeta_n in Z_p */
  { z = FpX_split(n, p, fl? 1: f); return fl? gel(z,1): z; }
  fK = znstar_conductor(znstar_generate(n, mkvecsmall(pmodn)));
  if (fK != n && umodiu(p, fK) == 1)
    z = FpX_split(fK, p, fl? 1: f);
  else
    z = FpX_factcyclo_just_conductor(fK, p, fl? 1: f);
  if (n > fK)
  {
    GEN vP = const_vec(n, NULL);
    long i, l = fl? 2: lg(z);
    for (i = 1; i < l; i++)
      gel(z, i) = FpX_conductor_lift(gel(z, i), p, fn, fK, vP);
  }
  return fl? gel(z,1): gen_sort(z,(void*)cmpii, &gen_cmp_RgX);
}

GEN
FpX_factcyclo(ulong n, GEN p, ulong m)
{ pari_sp av = avma; return gerepilecopy(av, FpX_factcyclo_i(n, p, m)); }

/*  Data = [H, GH, i_t, d0, kT, [n, d, f, n_T, mitk]]
 *  N2 = [p, pr, pu, pru] */
static GEN
Flx_pol_newton_general(GEN Data, GEN N2, GEN vT, ulong x)
{
  GEN i_t = gel(Data, 3), kT = gel(Data, 5), N = gel(Data, 6);
  long k, d = N[2], n_T = N[4], mitk = N[5];
  long p = N2[1], pr = N2[2], pu = N2[3], pru = N2[4];
  GEN R = cgetg(1+mitk, t_VECSMALL);

  for (k = 1; k <= n_T; k++) uel(R,kT[k]) = Flx_eval(gel(vT, kT[k]), x, pru) / pr;
  return Flx_Newton_perm(d, R, i_t, pu, p);
}

/* n is any integer prime to p, but must be equal to the conductor
 * of the splitting field of p in Q(zeta_n).
 * GH=G/H={g_1, ... ,g_f}
 * Data = [GHgen, GH, fn, p, [n, d, f, m]] */
static GEN
Flx_factcyclo_newton_general(GEN Data)
{
  GEN GH = gel(Data, 2), fn = gel(Data, 3), p = gel(Data, 4);
  ulong up = p[2], n = mael(Data, 5, 1), pmodn = up%n;
  long d = mael(Data, 5, 2), f = mael(Data, 5, 3), m = mael(Data, 5, 4);
  long i, k, n_T, mitk, r, s = 0, u = 1;
  GEN vT, kT, H, i_t, T, d0d1, Data2, Data3, R, Rs, pr, pu, pru;
  pari_timer ti;

  if (m != 1) m = f;
  for (pu = p; cmpiu(pu,d) <= 0; u++) pu = muliu(pu, up);  /* d<pu, pu=p^u */

  H = Fl_powers(pmodn, d-1, n); /* H=<p> */
  i_t = get_i_t(n, pmodn); /* i_t[1+i]=k ==> t_i=t_k */
  kT = get_kT(i_t, d); n_T = lg(kT)-1; mitk = kT[n_T];
  if (DEBUGLEVEL == 2) err_printf("kT=%Ps %ld elements\n",kT,n_T);
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  T = galoissubcyclo(utoi(n), utoi(pmodn), 0, 0);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "galoissubcyclo");
  d0d1 = get_d0_d1(T, gel(fn,1)); /* d0*T_k(x) is in Z[x] */
  Data2 = mkvecn(6, H, GH, i_t, d0d1, kT, mkvecsmalln(5, n, d, f, n_T, mitk));
  vT = get_vT(Data2, 0);
  if (DEBUGLEVEL == 4) err_printf("vT=%Ps\n",vT);
  r = QXV_den_pval(vT, kT, p);
  Rs = ZpX_roots_all(T, p, f, &s);
  if (DEBUGLEVEL >= 2) err_printf("(u,s,r)=(%ld,%ld,%ld)\n",u,s,r);
  if (r+u < s) pari_err_BUG("Flx_factcyclo_newton_general, T is not separable mod p^(r+u)");
  /* R and vT are mod p^(r+u) */
  R = (r+u==s) ? Rs : ZV_sort_shallow(ZX_Zp_liftroots(T, Rs, p, s, r+u));
  pr = powiu(p, r); pru = powiu(p, r+u); /* Usually, r=0, s=1, pr=1, pru=p */
  if (lgefint(pru) > 3)  /* ULONG_MAX < p^(r+u), probably won't occur */
  {
    for (k = 1; k <= n_T; k++)
    {
      long itk = kT[k];
      GEN z = r? RgX_Rg_mul(gel(vT, itk), pr): gel(vT, itk);
      gel(vT, itk) = RgX_to_FpX(z, pru);
    }
    Data3 = mkvec4(p, pr, pu, pru);
    for (i = 1; i <= m; i++)
      gel(R,i) = ZX_to_nx(FpX_pol_newton_general(Data2, Data3, vT, gel(R,i)));
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpX_pol_newton_general");
  }
  else
  {
    ulong upr = itou(pr), upru = itou(pru), upu = itou(pu);
    for (k = 1; k <= n_T; k++)
    {
      long itk = kT[k];
      GEN z = r? RgX_muls(gel(vT, itk), upr): gel(vT, itk);
      gel(vT, itk) = RgX_to_Flx(z, upru);
    }
    Data3 = mkvecsmall4(up, upr, upu, upru);
    if (DEBUGLEVEL >= 6) timer_start(&ti);
    for (i = 1; i <= m; i++)
      gel(R,i) = Flx_pol_newton_general(Data2, Data3, vT, itou(gel(R,i)));
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "Flx_pol_newton_general");
  }
  return R;
}

/*  Data=[vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
 *       [n, r, s, n_t, mitk, d], div] */
static GEN
Flx_pol_newton_general_new3(GEN Data, long k)
{
  GEN i_t = gel(Data,5), p = gel(Data,7), pu = gel(Data,8), prs = gel(Data,10);
  long d = mael(Data, 11, 6);
  GEN v_t_p = (lgefint(prs)>3)? ZV_to_nv(Fp_mk_v_t_p3(Data, k))
                              : Fl_mk_v_t_p3(Data, k);
  if (DEBUGLEVEL == 3) err_printf("v_t_p=%Ps\n",v_t_p);
  return Flx_Newton_perm(d, v_t_p, i_t, pu[2], p[2]);
}

/* Data = [GHgen, GH, fn, p, [n, d, f, m]] */
static GEN
Flx_factcyclo_newton_general_new3(GEN Data)
{
  long i, f = mael(Data, 5, 3), m = mael(Data, 5, 4);
  GEN Data2, pu, pols;
  pari_timer ti;

  if (m != 1) m = f;
  pols = cgetg(1+m, t_VEC);
  Data2 = newton_general_new_pre3(Data); pu = gel(Data2, 8);
  /* Data2 = [vT, gGH, Rs, Rrs, i_t, kt, p, pu, pr, prs,
              [n, r, s, n_t, mitk, d], div] */
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  if (lgefint(pu) > 3)  /* ULONG_MAX < p^u, probably won't occur */
  { for (i = 1; i <= m; i++)
      gel(pols, i) = ZX_to_nx(FpX_pol_newton_general_new3(Data2, i));
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "FpX_pol_newton_general_new3"); }
  else
  { for (i = 1; i <= m; i++)
      gel(pols, i) = Flx_pol_newton_general_new3(Data2, i);
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "Flx_pol_newton_general_new3"); }
  return pols;
}

/* return normalized z(-x) */
static GEN
Flx_1st_lift_2(GEN z, ulong p)
{
  long i, r = lg(z);
  GEN x = cgetg(r, t_VECSMALL);
  if (odd(r))
    for (i = 2; i<r; i++) uel(x,i) = odd(i)? Fl_neg(uel(z,i), p) : uel(z,i);
  else
    for (i = 2; i<r; i++) uel(x,i) = odd(i)? uel(z,i): Fl_neg(uel(z,i), p);
  return x;
}

/* el does not divides e.
 * construct Phi_{e*el}(x) from Phi_e(x). */
static GEN
Flx_1st_lift(GEN z, ulong p, ulong e, ulong el, GEN vP)
{
   GEN z2, z3, P = gel(vP, e);
   if (!gel(vP, e)) P = gel(vP, e) = Flx_polcyclo(e, p);
   z2 = Flx_inflate(z, el);
   z3 = Flx_normalize(Flx_gcd(P, z2, p), p);
   return Flx_div(z2, z3, p);
}

static GEN
Flx_lift(GEN z, ulong p, ulong e, ulong el, ulong r, ulong s, GEN vP)
{
  if (s == 0)
  {
    z = (el==2) ? Flx_1st_lift_2(z, p) : Flx_1st_lift(z, p, e, el, vP);
    if (r >= 2) z = Flx_inflate(z, upowuu(el, r-1));
  }
  else
    z = Flx_inflate(z, upowuu(el, r-s));
  return z;
}

/* e is the conductor of the splitting field of p in Q(zeta_n) */
static GEN
Flx_conductor_lift(GEN z, ulong p, GEN fn, ulong e, GEN vP)
{
  GEN EL = gel(fn, 1), En = gel(fn, 2);
  long i, r = lg(EL), new_e = e;

  for (i = 1; i < r; i++)
  {
    long el = EL[i], en = En[i], ee = u_lval(e, el);
    if (ee < en)
    {
      z = Flx_lift(z, p, new_e, el, en, ee, vP);
      new_e *= upowuu(el, en-ee);
    }
  }
  return z;
}

/* R0 is mod p^u, d < p^u */
static GEN
Flx_pol_newton(long j, GEN R0, GEN E, GEN D3, long d, long f, ulong p)
{
  ulong u = D3[3];
  GEN R = cgetg(f+1, t_VECSMALL);
  long i;
  for (i = 1; i <= f; i++) R[i] = R0[1+(i+j)%f];
  return Flx_Newton_perm(d, R, E, upowuu(p,u), p);
}

/* Data = [T, F, Rs, [d, nf, g, r, s, u]], nf>1 */
static GEN
Flx_factcyclo_newton_pre(GEN Data, GEN N, ulong p, ulong m)
{
  GEN T = gel(Data, 1), F = gel(Data, 2), Rs = gel(Data, 3), D4 = gel(Data, 4);
  long d = D4[1], nf = D4[2], g = D4[3], r = D4[4], s = D4[5], u = D4[6];
  GEN pols, E, R, p0 = utoi(p), Data3;
  ulong i, n = N[1], pmodn = p%n;
  pari_timer ti;

  if (m != 1) m = nf;
  pols = cgetg(1+m, t_VEC);
  E = set_E(pmodn, n, d, nf, g);
  R = set_R(T, F, Rs, p0, nf, r, s, u);
  R = ZV_to_nv(R);
  Data3 = mkvecsmall3(r, s, u);
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  for (i = 1; i <= m; i++) gel(pols, i) = Flx_pol_newton(i, R,E,Data3, d,nf,p);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "Flx_pol_newton");
  return pols;
}

/* gcd(n1,n2)=1, e = gcd(d1,d2).
 * phi(n1) = d1*nf1, phi(n2) = d2*nf2.
 * d = lcm(d1,d2) = d1*d2/e, nf = phi(n1)*phi(n2)/d = nf1*nf2*e. */
static GEN
Flx_factcyclo_lift(long n1, GEN v1, long n2, GEN v2, long p, long m)
{
  pari_sp av = avma;
  long d1 = degpol(gel(v1, 1)), lv1 = lg(v1), nf1 = eulerphiu(n1)/d1;
  long d2 = degpol(gel(v2, 1)), lv2 = lg(v2), nf2 = eulerphiu(n2)/d2;
  long e = ugcd(d1, d2), nf = nf1*nf2*e, need_factor = (e!=1);
  long i, j, k = 0;
  GEN v, z = NULL;
  pari_timer ti;

  if (DEBUGLEVEL >= 6) timer_start(&ti);
  if (m != 1) m = nf;
  v = cgetg(1+m, t_VEC);
  if (m == 1)
  {
    GEN x1 = gel(v1, 1), x2 = gel(v2, 1);
    z = Flx_gcd(Flx_inflate(x1, n2), Flx_inflate(x2, n1), p);
    z = Flx_normalize(z, p);  /* Flx_gcd sometimes returns non-monic */
    gel(v, 1) = (!need_factor)? z : gmael(Flx_factor(z, p), 1, 1);
  }
  else
    for (i = 1; i < lv1; i++)
      for (j = 1; j < lv2; j++)
      {
        GEN x1 = gel(v1, i), x2 = gel(v2, j);
        z = Flx_gcd(Flx_inflate(x1, n2), Flx_inflate(x2, n1), p);
        z = Flx_normalize(z, p);  /* needed */
        if (!need_factor) gel(v, ++k) = z;
        else
        {
          GEN z1 = gel(Flx_factor(z, p), 1);
          long i, l = lg(z1);
          for (i = 1; i < l; i++) gel(v, ++k) = gel(z1, i);
        }
      }
  if (DEBUGLEVEL >= 6)
    timer_printf(&ti, "Flx_factcyclo_lift (%ld,%ld)*(%ld,%ld)-->(%ld,%ld)-->(%ld,%ld)",
        d1, nf1, d2, nf2, degpol(z), nf1*nf2, d1*d2/e, nf);
  return gerepilecopy(av, v);
}

/* factor polcyclo(n) mod p based on an idea of Bill Allombert; d>1 and nf>1 */
static GEN
Flx_factcyclo_gen(GEN GH, ulong n, ulong p, ulong m)
{
  GEN fu = factoru(n), fa = zm_to_ZM(fu), p0 = utoi(p);
  GEN T, X, pd_n, v, pols;
  ulong i, j, pmodn = p%n, phin = eulerphiu_fact(fu);
  ulong d = Fl_order(pmodn, phin, n), nf = phin/d;
  pari_timer ti;

  if (m != 1) m = nf;
  if (m != 1 && !GH) /* Flx_factcyclo_fact is used */
  {
    GEN H = znstar_generate(n, mkvecsmall(pmodn));
    GH = znstar_cosets(n, phin, H); /* representatives of G/H */
  }

  pols = cgetg(1+m, t_VEC);
  v = cgetg(1+d, t_VEC);
  pd_n = diviuexact(subis(powiu(p0, d), 1), n); /* (p^d-1)/n */
  T = ZX_to_Flx(init_Fq(p0, d, evalvarn(1)), p);
  random_Flx(degpol(T), T[1], p);  /* skip 1st one */
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  do X = Flxq_pow(random_Flx(degpol(T), T[1], p), pd_n, T, p);
  while (lg(X)<3 || equaliu(Flxq_order(X, fa, T, p), n)==0); /* find zeta_n */
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "find X");

  if (m == 1)
  {
    for (j = 1; j <= d; j++)
    {
      gel(v, j) = polx_FlxX(0, 1);
      gmael(v, j, 2) = Flx_neg(X, p);
      if (j < d) X = Flxq_powu(X, p, T, p);
    }
    gel(pols, 1) = FlxX_to_Flx(FlxqXV_prod(v, T, p));
  }
  else
  {
    GEN vX, vp;
    if (DEBUGLEVEL >= 6) timer_start(&ti);
    vX = Flxq_powers(X, n, T, p)+1;
    vp = Fl_powers(pmodn, d, n);
    for (i = 1; i <= m; i++)
    {
      for (j = 1; j <= d; j++)
      {
        gel(v, j) = polx_FlxX(0, 1);
        gmael(v, j, 2) = Flx_neg(gel(vX, Fl_mul(GH[i], vp[j], n)), p);
      }
      gel(pols, i) = FlxX_to_Flx(FlxqXV_prod(v, T, p));
    }
    if (DEBUGLEVEL >= 6) timer_printf(&ti, "FlxqXV_prod");
  }
  return pols;
}

static int
cmpGuGu(GEN a, GEN b) { return (ulong)a < (ulong)b? -1: (a == b? 0: 1); }

/* p=1 (mod n). If m!=1, then m=phi(n) */
static GEN
Flx_split(ulong n, ulong p, ulong m)
{
  ulong i, j, z = rootsof1_Fl(n, p);
  GEN v, C, vx;
  if (m == 1) return mkvec(mkvecsmall3(0, Fl_neg(z,p), 1));
  v = cgetg(m+1, t_VEC); C = coprimes_zv(n); vx = Fl_powers(z, n-1, p);
  for (i = j = 1; i <= n; i++)
    if (C[i]) gel(v, j++) = mkvecsmall3(0, Fl_neg(vx[i+1], p), 1);
  return gen_sort(v, (void*)cmpGuGu, &gen_cmp_RgX);
}

/* d==1 or f==1 occurs */
static GEN
Flx_factcyclo_prime_power_i(long el, long e, long p, long m)
{
  GEN p0 = utoipos(p), z = set_e0_e1(el, e, p0), v;
  long n = z[1], e0 = z[2], e1 = z[3], phin = z[4], g = z[5];
  long d = z[6], f = z[7]; /* d and f for n=el^e0 */

  if (f == 1) v = mkvec(Flx_polcyclo(n, p));
  else if (d == 1) v = Flx_split(n, p, (m==1)?1:f);
  else if (el == 2) v = Flx_factcyclo_gen(NULL, n, p, m);/* d==2 in this case */
  else if (!use_newton(d, f)) v = Flx_factcyclo_gen(NULL, n, p, m);
  else
  {
    GEN N = mkvecsmall5(n, el, e0, phin, g);
    v = FpX_factcyclo_newton_power(N, p0, m, 1);
    if (typ(gel(v,1)) == t_POL)
    { /* ZX: convert to Flx */
      long i, l = lg(v);
      for (i = 1; i < l; i++) gel(v,i) = ZX_to_nx(gel(v,i));
    }
  }
  if (e1)
  {
    long i, l = lg(v), ele1 = upowuu(el, e1);
    for (i = 1; i < l; i++) gel(v, i) = Flx_inflate(gel(v, i), ele1);
  }
  return v;
}
static GEN
Flx_factcyclo_prime_power(long el, long e, long p, long m)
{
  pari_sp av = avma;
  return gerepilecopy(av, Flx_factcyclo_prime_power_i(el, e, p, m));
}

static GEN
Flx_factcyclo_fact(GEN fn, ulong p, ulong m, long ascent)
{
  GEN EL = gel(fn, 1), E = gel(fn, 2), v1, v2;
  long l = lg(EL), i, j, n1, n2;

  j = ascent? 1: l-1;
  n1 = upowuu(EL[j], E[j]);
  v1 = Flx_factcyclo_prime_power(EL[j], E[j], p, m);
  for (i = 2; i < l; i++)
  {
    j = ascent? i: l-i;
    n2 = upowuu(EL[j], E[j]);
    v2 = Flx_factcyclo_prime_power(EL[j], E[j], p, m);
    v1 = Flx_factcyclo_lift(n1, v1, n2, v2, p, m);
    n1 *= n2;
  }
  return v1;
}

static GEN
Flx_factcyclo_just_conductor(ulong n, ulong p, ulong m)
{
  GEN Data, fn;
  ulong action = FpX_factcyclo_just_conductor_init(&Data, n, utoipos(p), m);
  fn = gel(Data,3);
  if (action & GENERAL)
    return Flx_factcyclo_gen(gel(Data,2), n, p, m);
  else if (action & NEWTON_POWER)
    return Flx_factcyclo_prime_power_i(ucoeff(fn,1,1), ucoeff(fn,1,2), p, m);
  else if (action & NEWTON_GENERAL)
    return Flx_factcyclo_newton_general(Data);
  else if (action & NEWTON_GENERAL_NEW)
    return Flx_factcyclo_newton_general_new3(Data);
  else
    return Flx_factcyclo_fact(fn, p, m, action & ASCENT);
}

static GEN
Flx_factcyclo_i(ulong n, ulong p, ulong fl)
{
  GEN fn = factoru(n), z;
  ulong phin = eulerphiu_fact(fn), pmodn = p%n;
  ulong d = Fl_order(pmodn, phin, n), f = phin/d, fK;

  if (DEBUGLEVEL >= 1) header(fn, n, d, f, utoi(p));
  if (f == 1) { z = Flx_polcyclo(n, p); return fl? z: mkvec(z); }
  if (d == 1) /* p=1 (mod n), zeta_n in Z_p */
  { z = Flx_split(n, p, fl? 1: f); return fl? gel(z,1): z; }
  fK = znstar_conductor(znstar_generate(n, mkvecsmall(pmodn)));
  if (fK != n && p % fK == 1)
    z = Flx_split(fK, p, fl? 1: f);
  else
    z = Flx_factcyclo_just_conductor(fK, p, fl? 1: f);
  if (n > fK)
  {
    GEN vP = const_vec(n, NULL);
    long i, l = fl? 2: lg(z);
    for (i = 1; i < l; i++)
      gel(z, i) = Flx_conductor_lift(gel(z, i), p, fn, fK, vP);
  }
  return fl? gel(z,1): gen_sort(z,(void*)cmpGuGu, &gen_cmp_RgX);
}

GEN
Flx_factcyclo(ulong n, ulong p, ulong m)
{ pari_sp av = avma; return gerepilecopy(av, Flx_factcyclo_i(n, p, m)); }

GEN
factormodcyclo(long n, GEN p, long fl, long v)
{
  const char *fun = "factormodcyclo";
  pari_sp av = avma;
  long i, l;
  GEN z;
  if (v < 0) v = 0;
  if (fl < 0 || fl > 1) pari_err_FLAG(fun);
  if (n <= 0) pari_err_DOMAIN(fun, "n", "<=", gen_0, stoi(n));
  if (typ(p) != t_INT) pari_err_TYPE(fun, p);
  if (dvdui(n, p)) pari_err_COPRIME(fun, stoi(n), p);
  if (fl)
  {
    if (lgefint(p) == 3)
      z = Flx_to_ZX(Flx_factcyclo_i(n, p[2], 1));
    else
      z = FpX_factcyclo_i(n, p, 1);
    setvarn(z, v);
    return gerepileupto(av, FpX_to_mod(z, p));
  }
  else
  {
    if (lgefint(p) == 3)
      z = FlxC_to_ZXC(Flx_factcyclo_i(n, p[2], 0));
    else
      z = FpX_factcyclo_i(n, p, 0);
    l = lg(z); for (i = 1; i < l; i++) setvarn(gel(z, i), v);
    return gerepileupto(av, FpXC_to_mod(z, p));
  }
}
