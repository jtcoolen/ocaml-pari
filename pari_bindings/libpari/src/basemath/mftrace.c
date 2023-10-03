/* Copyright (C) 2016  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/*************************************************************************/
/*                                                                       */
/*              Modular forms package based on trace formulas            */
/*                                                                       */
/*************************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_mf

enum {
  MF_SPLIT = 1,
  MF_EISENSPACE,
  MF_FRICKE,
  MF_MF2INIT,
  MF_SPLITN
};

typedef struct {
  GEN vnew, vfull, DATA, VCHIP;
  long n, newHIT, newTOTAL, cuspHIT, cuspTOTAL;
} cachenew_t;

static void init_cachenew(cachenew_t *c, long n, long N, GEN f);
static long mf1cuspdim_i(long N, GEN CHI, GEN TMP, GEN vSP, long *dih);
static GEN mfinit_i(GEN NK, long space);
static GEN mfinit_Nkchi(long N, long k, GEN CHI, long space, long flraw);
static GEN mf2init_Nkchi(long N, long k, GEN CHI, long space, long flraw);
static GEN mf2basis(long N, long r, GEN CHI, GEN *pCHI1, long space);
static GEN mfeisensteinbasis(long N, long k, GEN CHI);
static GEN mfeisensteindec(GEN mf, GEN F);
static GEN initwt1newtrace(GEN mf);
static GEN initwt1trace(GEN mf);
static GEN myfactoru(long N);
static GEN mydivisorsu(long N);
static GEN Qab_Czeta(long k, long ord, GEN C, long vt);
static GEN mfcoefs_i(GEN F, long n, long d);
static GEN bhnmat_extend(GEN M, long m,long l, GEN S, cachenew_t *cache);
static GEN initnewtrace(long N, GEN CHI);
static void dbg_cachenew(cachenew_t *C);
static GEN hecke_i(long m, long l, GEN V, GEN F, GEN DATA);
static GEN c_Ek(long n, long d, GEN F);
static GEN RgV_heckef2(long n, long d, GEN V, GEN F, GEN DATA);
static GEN mfcusptrace_i(long N, long k, long n, GEN Dn, GEN TDATA);
static GEN mfnewtracecache(long N, long k, long n, cachenew_t *cache);
static GEN colnewtrace(long n0, long n, long d, long N, long k, cachenew_t *c);
static GEN dihan(GEN bnr, GEN w, GEN k0j, long m, ulong n);
static GEN sigchi(long k, GEN CHI, long n);
static GEN sigchi2(long k, GEN CHI1, GEN CHI2, long n, long ord);
static GEN mflineardivtomat(long N, GEN vF, long n);
static GEN mfdihedralcusp(long N, GEN CHI, GEN vSP);
static long mfdihedralcuspdim(long N, GEN CHI, GEN vSP);
static GEN mfdihedralnew(long N, GEN CHI, GEN SP);
static GEN mfdihedral(long N);
static GEN mfdihedralall(long N);
static long mf1cuspdim(long N, GEN CHI, GEN vSP);
static long mf2dim_Nkchi(long N, long k, GEN CHI, ulong space);
static long mfdim_Nkchi(long N, long k, GEN CHI, long space);
static GEN charLFwtk(long N, long k, GEN CHI, long ord, long t);
static GEN mfeisensteingacx(GEN E,long w,GEN ga,long n,long prec);
static GEN mfgaexpansion(GEN mf, GEN F, GEN gamma, long n, long prec);
static GEN mfEHmat(long n, long r);
static GEN mfEHcoef(long r, long N);
static GEN mftobasis_i(GEN mf, GEN F);

static GEN
mkgNK(GEN N, GEN k, GEN CHI, GEN P) { return mkvec4(N, k, CHI, P); }
static GEN
mkNK(long N, long k, GEN CHI) { return mkgNK(stoi(N), stoi(k), CHI, pol_x(1)); }
GEN
MF_get_CHI(GEN mf) { return gmael(mf,1,3); }
GEN
MF_get_gN(GEN mf) { return gmael(mf,1,1); }
long
MF_get_N(GEN mf) { return itou(MF_get_gN(mf)); }
GEN
MF_get_gk(GEN mf) { return gmael(mf,1,2); }
long
MF_get_k(GEN mf)
{
  GEN gk = MF_get_gk(mf);
  if (typ(gk)!=t_INT) pari_err_IMPL("half-integral weight");
  return itou(gk);
}
long
MF_get_r(GEN mf)
{
  GEN gk = MF_get_gk(mf);
  if (typ(gk) == t_INT) pari_err_IMPL("integral weight");
  return itou(gel(gk, 1)) >> 1;
}
long
MF_get_space(GEN mf) { return itos(gmael(mf,1,4)); }
GEN
MF_get_E(GEN mf) { return gel(mf,2); }
GEN
MF_get_S(GEN mf) { return gel(mf,3); }
GEN
MF_get_basis(GEN mf) { return shallowconcat(gel(mf,2), gel(mf,3)); }
long
MF_get_dim(GEN mf)
{
  switch(MF_get_space(mf))
  {
    case mf_FULL:
      return lg(MF_get_S(mf)) - 1 + lg(MF_get_E(mf))-1;
    case mf_EISEN:
      return lg(MF_get_E(mf))-1;
    default: /* mf_NEW, mf_CUSP, mf_OLD */
      return lg(MF_get_S(mf)) - 1;
  }
}
GEN
MFnew_get_vj(GEN mf) { return gel(mf,4); }
GEN
MFcusp_get_vMjd(GEN mf) { return gel(mf,4); }
GEN
MF_get_M(GEN mf) { return gmael(mf,5,3); }
GEN
MF_get_Minv(GEN mf) { return gmael(mf,5,2); }
GEN
MF_get_Mindex(GEN mf) { return gmael(mf,5,1); }

/* ordinary gtocol forgets about initial 0s */
GEN
sertocol(GEN S) { return gtocol0(S, -(lg(S) - 2 + valser(S))); }
/*******************************************************************/
/*     Linear algebra in cyclotomic fields (TODO: export this)     */
/*******************************************************************/
/* return r and split prime p giving projection Q(zeta_n) -> Fp, zeta -> r */
static ulong
QabM_init(long n, ulong *p)
{
  ulong pinit = 1000000007;
  forprime_t T;
  if (n <= 1) { *p = pinit; return 0; }
  u_forprime_arith_init(&T, pinit, ULONG_MAX, 1, n);
  *p = u_forprime_next(&T);
  return Flx_oneroot(ZX_to_Flx(polcyclo(n, 0), *p), *p);
}
static ulong
Qab_to_Fl(GEN P, ulong r, ulong p)
{
  ulong t;
  GEN den;
  P = Q_remove_denom(liftpol_shallow(P), &den);
  if (typ(P) == t_POL) { GEN Pp = ZX_to_Flx(P, p); t = Flx_eval(Pp, r, p); }
  else t = umodiu(P, p);
  if (den) t = Fl_div(t, umodiu(den, p), p);
  return t;
}
static GEN
QabC_to_Flc(GEN C, ulong r, ulong p)
{
  long i, l = lg(C);
  GEN A = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) uel(A,i) = Qab_to_Fl(gel(C,i), r, p);
  return A;
}
static GEN
QabM_to_Flm(GEN M, ulong r, ulong p)
{
  long i, l;
  GEN A = cgetg_copy(M, &l);
  for (i = 1; i < l; i++)
    gel(A, i) = QabC_to_Flc(gel(M, i), r, p);
  return A;
}
/* A a t_POL */
static GEN
QabX_to_Flx(GEN A, ulong r, ulong p)
{
  long i, l = lg(A);
  GEN a = cgetg(l, t_VECSMALL);
  a[1] = ((ulong)A[1])&VARNBITS;
  for (i = 2; i < l; i++) uel(a,i) = Qab_to_Fl(gel(A,i), r, p);
  return Flx_renormalize(a, l);
}

/* FIXME: remove */
static GEN
ZabM_pseudoinv_i(GEN M, GEN P, long n, GEN *pv, GEN *den, int ratlift)
{
  GEN v = ZabM_indexrank(M, P, n);
  if (pv) *pv = v;
  M = shallowmatextract(M,gel(v,1),gel(v,2));
  return ratlift? ZabM_inv_ratlift(M, P, n, den): ZabM_inv(M, P, n, den);
}

/* M matrix with coeff in Q(\chi)), where Q(\chi) = Q(X)/(P) for
 * P = cyclotomic Phi_n. Assume M rational if n <= 2 */
static GEN
QabM_ker(GEN M, GEN P, long n)
{
  if (n <= 2) return QM_ker(M);
  return ZabM_ker(row_Q_primpart(liftpol_shallow(M)), P, n);
}
/* pseudo-inverse of M. FIXME: should replace QabM_pseudoinv */
static GEN
QabM_pseudoinv_i(GEN M, GEN P, long n, GEN *pv, GEN *pden)
{
  GEN cM, Mi;
  if (n <= 2)
  {
    M = Q_primitive_part(M, &cM);
    Mi = ZM_pseudoinv(M, pv, pden); /* M^(-1) = Mi / (cM * den) */
  }
  else
  {
    M = Q_primitive_part(liftpol_shallow(M), &cM);
    Mi = ZabM_pseudoinv(M, P, n, pv, pden);
  }
  *pden = mul_content(*pden, cM);
  return Mi;
}
/* FIXME: delete */
static GEN
QabM_pseudoinv(GEN M, GEN P, long n, GEN *pv, GEN *pden)
{
  GEN Mi = QabM_pseudoinv_i(M, P, n, pv, pden);
  return P? gmodulo(Mi, P): Mi;
}

static GEN
QabM_indexrank(GEN M, GEN P, long n)
{
  GEN z;
  if (n <= 2)
  {
    M = vec_Q_primpart(M);
    z = ZM_indexrank(M); /* M^(-1) = Mi / (cM * den) */
  }
  else
  {
    M = vec_Q_primpart(liftpol_shallow(M));
    z = ZabM_indexrank(M, P, n);
  }
  return z;
}

/*********************************************************************/
/*                    Simple arithmetic functions                    */
/*********************************************************************/
/* TODO: most of these should be exported and used in ifactor1.c */
/* phi(n) */
static ulong
myeulerphiu(ulong n)
{
  pari_sp av;
  if (n == 1) return 1;
  av = avma; return gc_ulong(av, eulerphiu_fact(myfactoru(n)));
}
static long
mymoebiusu(ulong n)
{
  pari_sp av;
  if (n == 1) return 1;
  av = avma; return gc_long(av, moebiusu_fact(myfactoru(n)));
}

static long
mynumdivu(long N)
{
  pari_sp av;
  if (N == 1) return 1;
  av = avma; return gc_long(av, numdivu_fact(myfactoru(N)));
}

/* N\prod_{p|N} (1+1/p) */
static long
mypsiu(ulong N)
{
  pari_sp av;
  GEN P;
  long j, l, a;
  if (N == 1) return 1;
  av = avma; P = gel(myfactoru(N), 1); l = lg(P);
  for (a = N, j = 1; j < l; j++) a += a / P[j];
  return gc_long(av, a);
}
/* write n = mf^2. Return m, set f. */
static ulong
mycore(ulong n, long *pf)
{
  pari_sp av = avma;
  GEN fa = myfactoru(n), P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), m = 1, f = 1;
  for (i = 1; i < l; i++)
  {
    long j, p = P[i], e = E[i];
    if (e & 1) m *= p;
    for (j = 2; j <= e; j+=2) f *= p;
  }
  *pf = f; return gc_long(av,m);
}

/* fa = factorization of -D > 0, return -D0 > 0 (where D0 is fundamental) */
static long
corediscs_fact(GEN fa)
{
  GEN P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), m = 1;
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (e & 1) m *= p;
  }
  if ((m&3L) != 3) m <<= 2;
  return m;
}
static long
mubeta(long n)
{
  pari_sp av = avma;
  GEN E = gel(myfactoru(n), 2);
  long i, s = 1, l = lg(E);
  for (i = 1; i < l; i++)
  {
    long e = E[i];
    if (e >= 3) return gc_long(av,0);
    if (e == 1) s *= -2;
  }
  return gc_long(av,s);
}

/* n = n1*n2, n1 = ppo(n, m); return mubeta(n1)*moebiusu(n2).
 * N.B. If n from newt_params we, in fact, never return 0 */
static long
mubeta2(long n, long m)
{
  pari_sp av = avma;
  GEN fa = myfactoru(n), P = gel(fa,1), E = gel(fa,2);
  long i, s = 1, l = lg(P);
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (m % p)
    { /* p^e in n1 */
      if (e >= 3) return gc_long(av,0);
      if (e == 1) s *= -2;
    }
    else
    { /* in n2 */
      if (e >= 2) return gc_long(av,0);
      s = -s;
    }
  }
  return gc_long(av,s);
}

/* write N = prod p^{ep} and n = df^2, d squarefree.
 * set g  = ppo(gcd(sqfpart(N), f), FC)
 *     N2 = prod p^if(e==1 || p|n, ep-1, ep-2) */
static void
newt_params(long N, long n, long FC, long *pg, long *pN2)
{
  GEN fa = myfactoru(N), P = gel(fa,1), E = gel(fa,2);
  long i, g = 1, N2 = 1, l = lg(P);
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (e == 1)
    { if (FC % p && n % (p*p) == 0) g *= p; }
    else
      N2 *= upowuu(p,(n % p)? e-2: e-1);
  }
  *pg = g; *pN2 = N2;
}
/* simplified version of newt_params for n = 1 (newdim) */
static void
newd_params(long N, long *pN2)
{
  GEN fa = myfactoru(N), P = gel(fa,1), E = gel(fa,2);
  long i, N2 = 1, l = lg(P);
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (e > 2) N2 *= upowuu(p, e-2);
  }
  *pN2 = N2;
}

static long
newd_params2(long N)
{
  GEN fa = myfactoru(N), P = gel(fa,1), E = gel(fa,2);
  long i, N2 = 1, l = lg(P);
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (e >= 2) N2 *= upowuu(p, e);
  }
  return N2;
}

/*******************************************************************/
/*   Relative trace between cyclotomic fields (TODO: export this)  */
/*******************************************************************/
/* g>=1; return g * prod_{p | g, (p,q) = 1} (1-1/p) */
static long
phipart(long g, long q)
{
  if (g > 1)
  {
    GEN P = gel(myfactoru(g), 1);
    long i, l = lg(P);
    for (i = 1; i < l; i++) { long p = P[i]; if (q % p) g -= g / p; }
  }
  return g;
}
/* Set s,v s.t. Trace(zeta_N^k) from Q(zeta_N) to Q(\zeta_N) = s * zeta_M^v
 * With k > 0, N = M*d and N, M != 2 mod 4 */
static long
tracerelz(long *pv, long d, long M, long k)
{
  long s, g, q, muq;
  if (d == 1) { *pv = k; return 1; }
  *pv = 0; g = ugcd(k, d); q = d / g;
  muq = mymoebiusu(q); if (!muq) return 0;
  if (M != 1)
  {
    long v = Fl_invsafe(q % M, M);
    if (!v) return 0;
    *pv = (v * (k/g)) % M;
  }
  s = phipart(g, M*q); if (muq < 0) s = -s;
  return s;
}
/* Pi = polcyclo(i), i = m or n. Let Ki = Q(zeta_i), initialize Tr_{Kn/Km} */
GEN
Qab_trace_init(long n, long m, GEN Pn, GEN Pm)
{
  long a, i, j, N, M, vt, d, D;
  GEN T, G;

  if (m == n || n <= 2) return mkvec(Pm);
  vt = varn(Pn);
  d = degpol(Pn);
  /* if (N != n) zeta_N = zeta_n^2 and zeta_n = - zeta_N^{(N+1)/2} */
  N = ((n & 3) == 2)? n >> 1: n;
  M = ((m & 3) == 2)? m >> 1: m; /* M | N | n */
  a = N / M;
  T = const_vec(d, NULL);
  D = d / degpol(Pm); /* relative degree */
  if (D == 1) G = NULL;
  else
  { /* zeta_M = zeta_n^A; s_j(zeta_M) = zeta_M <=> j = 1 (mod J) */
    long lG, A = (N == n)? a: (a << 1), J = n / ugcd(n, A);
    G = coprimes_zv(n);
    for (j = lG = 1; j < n; j += J)
      if (G[j]) G[lG++] = j;
    setlg(G, lG); /* Gal(Q(zeta_n) / Q(zeta_m)) */
  }
  T = const_vec(d, NULL);
  gel(T,1) = utoipos(D); /* Tr 1 */
  for (i = 1; i < d; i++)
  { /* if n = 2N, zeta_n^i = (-1)^i zeta_N^k */
    long s, v, k;
    GEN t;

    if (gel(T, i+1)) continue;
    k = (N == n)? i: ((odd(i)? i + N: i) >> 1);
    if ((s = tracerelz(&v, a, M, k)))
    {
      if (m != M) v *= 2;/* Tr = s * zeta_m^v */
      if (n != N && odd(i)) s = -s;
      t = Qab_Czeta(v, m, stoi(s), vt);
    }
    else
      t = gen_0;
    /* t = Tr_{Kn/Km} zeta_n^i; fill using Galois action */
    if (!G)
      gel(T, i + 1) = t;
    else
      for (j = 1; j <= D; j++)
      {
        long z = Fl_mul(i,G[j], n);
        if (z < d) gel(T, z + 1) = t;
      }
  }
  return mkvec3(Pm, Pn, T);
}
/* x a t_POL modulo Phi_n */
static GEN
tracerel_i(GEN T, GEN x)
{
  long k, l = lg(x);
  GEN S;
  if (l == 2) return gen_0;
  S = gmul(gel(T,1), gel(x,2));
  for (k = 3; k < l; k++) S = gadd(S, gmul(gel(T,k-1), gel(x,k)));
  return S;
}
static GEN
tracerel(GEN a, GEN v, GEN z)
{
  a = liftpol_shallow(a);
  a = simplify_shallow(z? gmul(z,a): a);
  if (typ(a) == t_POL)
  {
    GEN T = gel(v,3);
    long degrel = itou(gel(T,1));
    a = tracerel_i(T, RgX_rem(a, gel(v,2)));
    if (degrel != 1) a = gdivgu(a, degrel);
    if (typ(a) == t_POL) a = RgX_rem(a, gel(v,1));
  }
  return a;
}
static GEN
tracerel_z(GEN v, long t)
{
  GEN Pn = gel(v,2);
  return t? pol_xn(t, varn(Pn)): NULL;
}
/* v = Qab_trace_init(n,m); x is a t_VEC of polmodulo Phi_n; Kn = Q(zeta_n)
 * [Kn:Km]^(-1) Tr_{Kn/Km} (zeta_n^t * x); 0 <= t < [Kn:Km] */
GEN
Qab_tracerel(GEN v, long t, GEN a)
{
  if (lg(v) != 4) return a; /* => t = 0 */
  return tracerel(a, v, tracerel_z(v, t));
}
GEN
QabV_tracerel(GEN v, long t, GEN x)
{
  GEN z;
  if (lg(v) != 4) return x; /* => t = 0 */
  z = tracerel_z(v, t);
  pari_APPLY_same(tracerel(gel(x,i), v, z));
}
GEN
QabM_tracerel(GEN v, long t, GEN x)
{
  if (lg(v) != 4) return x;
  pari_APPLY_same(QabV_tracerel(v, t, gel(x,i)));
}

/* C*zeta_o^k mod X^o - 1 */
static GEN
Qab_Czeta(long k, long o, GEN C, long vt)
{
  if (!k) return C;
  if (!odd(o))
  { /* optimization: reduce max degree by a factor 2 for free */
    o >>= 1;
    if (k >= o) { k -= o; C = gneg(C); if (!k) return C; }
  }
  return monomial(C, k, vt);
}
/* zeta_o^k */
static GEN
Qab_zeta(long k, long o, long vt) { return Qab_Czeta(k, o, gen_1, vt); }

/*              Operations on Dirichlet characters                       */

/* A Dirichlet character can be given in GP in different formats, but in this
 * package, it will be a vector CHI=[G,chi,ord], where G is the (Z/MZ)^* to
 * which the character belongs, chi is the character in Conrey format, ord is
 * the order */

static GEN
gmfcharorder(GEN CHI) { return gel(CHI, 3); }
long
mfcharorder(GEN CHI) { return itou(gmfcharorder(CHI)); }
static long
mfcharistrivial(GEN CHI) { return !CHI || mfcharorder(CHI) == 1; }
static GEN
gmfcharmodulus(GEN CHI) { return gmael3(CHI, 1, 1, 1); }
long
mfcharmodulus(GEN CHI) { return itou(gmfcharmodulus(CHI)); }
GEN
mfcharpol(GEN CHI) { return gel(CHI,4); }

/* vz[i+1] = image of (zeta_o)^i in Fp */
static ulong
Qab_Czeta_Fl(long k, GEN vz, ulong C, ulong p)
{
  long o;
  if (!k) return C;
  o = lg(vz)-2;
  if ((k << 1) == o) return Fl_neg(C,p);
  return Fl_mul(C, vz[k+1], p);
}

static long
znchareval_i(GEN CHI, long n, GEN ord)
{ return itos(znchareval(gel(CHI,1), gel(CHI,2), stoi(n), ord)); }

/* n coprime with the modulus of CHI */
static GEN
mfchareval(GEN CHI, long n)
{
  GEN Pn, C, go = gmfcharorder(CHI);
  long k, o = go[2];
  if (o == 1) return gen_1;
  k = znchareval_i(CHI, n, go);
  Pn = mfcharpol(CHI);
  C = Qab_zeta(k, o, varn(Pn));
  if (typ(C) != t_POL) return C;
  return gmodulo(C, Pn);
}
/* d a multiple of ord(CHI); n coprime with char modulus;
 * return x s.t. CHI(n) = \zeta_d^x] */
static long
mfcharevalord(GEN CHI, long n, long d)
{
  if (mfcharorder(CHI) == 1) return 0;
  return znchareval_i(CHI, n, utoi(d));
}

/* G a znstar, L a Conrey log: return a 'mfchar' */
static GEN
mfcharGL(GEN G, GEN L)
{
  GEN o = zncharorder(G,L);
  long ord = itou(o), vt = fetch_user_var("t");
  return mkvec4(G, L, o, polcyclo(ord,vt));
}
static GEN
mfchartrivial()
{ return mfcharGL(znstar0(gen_1,1), cgetg(1,t_COL)); }
/* convert a generic character into an 'mfchar' */
static GEN
get_mfchar(GEN CHI)
{
  GEN G, L;
  if (typ(CHI) != t_VEC) CHI = znchar(CHI);
  else
  {
    long l = lg(CHI);
    if ((l != 3 && l != 5) || !checkznstar_i(gel(CHI,1)))
      pari_err_TYPE("checkNF [chi]", CHI);
    if (l == 5) return CHI;
  }
  G = gel(CHI,1);
  L = gel(CHI,2); if (typ(L) != t_COL) L = znconreylog(G,L);
  return mfcharGL(G, L);
}

/* parse [N], [N,k], [N,k,CHI]. If 'joker' is set, allow wildcard for CHI */
static GEN
checkCHI(GEN NK, long N, int joker)
{
  GEN CHI;
  if (lg(NK) == 3)
    CHI = mfchartrivial();
  else
  {
    long i, l;
    CHI = gel(NK,3); l = lg(CHI);
    if (isintzero(CHI) && joker)
      CHI = NULL; /* all character orbits */
    else if (isintm1(CHI) && joker > 1)
      CHI = gen_m1; /* sum over all character orbits */
    else if ((typ(CHI) == t_VEC &&
             (l == 1 || l != 3 || !checkznstar_i(gel(CHI,1)))) && joker)
    {
      CHI = shallowtrans(CHI); /* list of characters */
      for (i = 1; i < l; i++) gel(CHI,i) = get_mfchar(gel(CHI,i));
    }
    else
    {
      CHI = get_mfchar(CHI); /* single char */
      if (N % mfcharmodulus(CHI)) pari_err_TYPE("checkNF [chi]", NK);
    }
  }
  return CHI;
}
/* support half-integral weight */
static void
checkNK2(GEN NK, long *N, long *nk, long *dk, GEN *CHI, int joker)
{
  long l = lg(NK);
  GEN T;
  if (typ(NK) != t_VEC || l < 3 || l > 4) pari_err_TYPE("checkNK", NK);
  T = gel(NK,1); if (typ(T) != t_INT) pari_err_TYPE("checkNF [N]", NK);
  *N = itos(T); if (*N <= 0) pari_err_TYPE("checkNF [N <= 0]", NK);
  T = gel(NK,2);
  switch(typ(T))
  {
    case t_INT:  *nk = itos(T); *dk = 1; break;
    case t_FRAC:
      *nk = itos(gel(T,1));
      *dk = itou(gel(T,2)); if (*dk == 2) break;
    default: pari_err_TYPE("checkNF [k]", NK);
  }
  *CHI = checkCHI(NK, *N, joker);
}
/* don't support half-integral weight */
static void
checkNK(GEN NK, long *N, long *k, GEN *CHI, int joker)
{
  long d;
  checkNK2(NK, N, k, &d, CHI, joker);
  if (d != 1) pari_err_TYPE("checkNF [k]", NK);
}

static GEN
mfchargalois(long N, int odd, GEN flagorder)
{
  GEN G = znstar0(utoi(N), 1), L = chargalois(G, flagorder);
  long l = lg(L), i, j;
  for (i = j = 1; i < l; i++)
  {
    GEN chi = znconreyfromchar(G, gel(L,i));
    if (zncharisodd(G,chi) == odd) gel(L,j++) = mfcharGL(G,chi);
  }
  setlg(L, j); return L;
}
/* possible characters for nontrivial S_1(N, chi) */
static GEN
mf1chars(long N, GEN vCHI)
{
  if (vCHI) return vCHI; /*do not filter, user knows best*/
  /* Tate's theorem */
  return mfchargalois(N, 1, uisprime(N)? mkvecsmall2(2,4): NULL);
}
static GEN
mfchars(long N, long k, long dk, GEN vCHI)
{ return vCHI? vCHI: mfchargalois(N, (dk == 2)? 0: (k & 1), NULL); }

/* wrappers from mfchar to znchar */
static long
mfcharparity(GEN CHI)
{
  if (!CHI) return 1;
  return zncharisodd(gel(CHI,1), gel(CHI,2)) ? -1 : 1;
}
/* if CHI is primitive, return CHI itself, not a copy */
static GEN
mfchartoprimitive(GEN CHI, long *pF)
{
  pari_sp av;
  GEN chi, F;
  if (!CHI) { if (pF) *pF = 1; return mfchartrivial(); }
  av = avma; F = znconreyconductor(gel(CHI,1), gel(CHI,2), &chi);
  if (typ(F) == t_INT) set_avma(av);
  else
  {
    CHI = leafcopy(CHI);
    gel(CHI,1) = znstar0(F, 1);
    gel(CHI,2) = chi;
  }
  if (pF) *pF = mfcharmodulus(CHI);
  return CHI;
}
static long
mfcharconductor(GEN CHI)
{
  pari_sp av = avma;
  GEN res = znconreyconductor(gel(CHI,1), gel(CHI,2), NULL);
  if (typ(res) == t_VEC) res = gel(res, 1);
  return gc_long(av, itos(res));
}

/*                      Operations on mf closures                    */
static GEN
tagparams(long t, GEN NK) { return mkvec2(mkvecsmall(t), NK); }
static GEN
lfuntag(long t, GEN x) { return mkvec2(mkvecsmall(t), x); }
static GEN
tag0(long t, GEN NK) { retmkvec(tagparams(t,NK)); }
static GEN
tag(long t, GEN NK, GEN x) { retmkvec2(tagparams(t,NK), x); }
static GEN
tag2(long t, GEN NK, GEN x, GEN y) { retmkvec3(tagparams(t,NK), x,y); }
static GEN
tag3(long t, GEN NK, GEN x,GEN y,GEN z) { retmkvec4(tagparams(t,NK), x,y,z); }
static GEN
tag4(long t, GEN NK, GEN x,GEN y,GEN z,GEN a)
{ retmkvec5(tagparams(t,NK), x,y,z,a); }
/* is F a "modular form" ? */
int
checkmf_i(GEN F)
{ return typ(F) == t_VEC
    && lg(F) > 1 && typ(gel(F,1)) == t_VEC
    && lg(gel(F,1)) == 3
    && typ(gmael(F,1,1)) == t_VECSMALL
    && typ(gmael(F,1,2)) == t_VEC; }
long mf_get_type(GEN F) { return gmael(F,1,1)[1]; }
GEN mf_get_gN(GEN F) { return gmael3(F,1,2,1); }
GEN mf_get_gk(GEN F) { return gmael3(F,1,2,2); }
/* k - 1/2, assume k in 1/2 + Z */
long mf_get_r(GEN F) { return itou(gel(mf_get_gk(F),1)) >> 1; }
long mf_get_N(GEN F) { return itou(mf_get_gN(F)); }
long mf_get_k(GEN F)
{
  GEN gk = mf_get_gk(F);
  if (typ(gk)!=t_INT) pari_err_IMPL("half-integral weight");
  return itou(gk);
}
GEN mf_get_CHI(GEN F) { return gmael3(F,1,2,3); }
GEN mf_get_field(GEN F) { return gmael3(F,1,2,4); }
GEN mf_get_NK(GEN F) { return gmael(F,1,2); }
static void
mf_setfield(GEN f, GEN P)
{
  gel(f,1) = leafcopy(gel(f,1));
  gmael(f,1,2) = leafcopy(gmael(f,1,2));
  gmael3(f,1,2,4) = P;
}

/* UTILITY FUNCTIONS */
GEN
mftocol(GEN F, long lim, long d)
{ GEN c = mfcoefs_i(F, lim, d); settyp(c,t_COL); return c; }
GEN
mfvectomat(GEN vF, long lim, long d)
{
  long j, l = lg(vF);
  GEN M = cgetg(l, t_MAT);
  for (j = 1; j < l; j++) gel(M,j) = mftocol(gel(vF,j), lim, d);
  return M;
}

static GEN
RgV_to_ser_full(GEN x) { return RgV_to_ser(x, 0, lg(x)+1); }
/* TODO: delete */
static GEN
mfcoefsser(GEN F, long n) { return RgV_to_ser_full(mfcoefs_i(F,n,1)); }
static GEN
sertovecslice(GEN S, long n)
{
  GEN v = gtovec0(S, -(lg(S) - 2 + valser(S)));
  long l = lg(v), n2 = n + 2;
  if (l < n2) pari_err_BUG("sertovecslice [n too large]");
  return (l == n2)? v: vecslice(v, 1, n2-1);
}

/* a, b two RgV of the same length, multiply as truncated power series */
static GEN
RgV_mul_RgXn(GEN a, GEN b)
{
  long n = lg(a)-1;
  GEN c;
  a = RgV_to_RgX(a,0);
  b = RgV_to_RgX(b,0); c = RgXn_mul(a, b, n);
  c = RgX_to_RgC(c,n); settyp(c,t_VEC); return c;
}
/* divide as truncated power series */
static GEN
RgV_div_RgXn(GEN a, GEN b)
{
  long n = lg(a)-1;
  GEN c;
  a = RgV_to_RgX(a,0);
  b = RgV_to_RgX(b,0); c = RgXn_div_i(a, b, n);
  c = RgX_to_RgC(c,n); settyp(c,t_VEC); return c;
}
/* a^b */
static GEN
RgV_pows_RgXn(GEN a, long b)
{
  long n = lg(a)-1;
  GEN c;
  a = RgV_to_RgX(a,0);
  if (b < 0) { a = RgXn_inv(a, n); b = -b; }
  c = RgXn_powu_i(a,b,n);
  c = RgX_to_RgC(c,n); settyp(c,t_VEC); return c;
}

/* assume lg(V) >= n*d + 2 */
static GEN
c_deflate(long n, long d, GEN v)
{
  long i, id, l = n+2;
  GEN w;
  if (d == 1) return lg(v) == l ? v: vecslice(v, 1, l-1);
  w = cgetg(l, typ(v));
  for (i = id = 1; i < l; i++, id += d) gel(w, i) = gel(v, id);
  return w;
}

static void
err_cyclo(void)
{ pari_err_IMPL("changing cyclotomic fields in mf"); }
/* Q(zeta_a) = Q(zeta_b) ? */
static int
same_cyc(long a, long b)
{ return (a == b) || (odd(a) && b == (a<<1)) || (odd(b) && a == (b<<1)); }
/* need to combine elements in Q(CHI1) and Q(CHI2) with result in Q(CHI),
 * CHI = CHI1 * CHI2 or CHI / CHI2 times some character of order 2 */
static GEN
chicompat(GEN CHI, GEN CHI1, GEN CHI2)
{
  long o1 = mfcharorder(CHI1);
  long o2 = mfcharorder(CHI2), O, o;
  GEN T1, T2, P, Po;
  if (o1 <= 2 && o2 <= 2) return NULL;
  o = mfcharorder(CHI);
  Po = mfcharpol(CHI);
  P = mfcharpol(CHI1);
  if (o1 == o2)
  {
    if (o1 == o) return NULL;
    if (!same_cyc(o1,o)) err_cyclo();
    return mkvec4(P, gen_1,gen_1, Qab_trace_init(o1, o, P, Po));
  }
  O = ulcm(o1, o2);
  if (!same_cyc(O,o)) err_cyclo();
  if (O != o1) P = (O == o2)? mfcharpol(CHI2): polcyclo(O, varn(P));
  T1 = o1 <= 2? gen_1: utoipos(O / o1);
  T2 = o2 <= 2? gen_1: utoipos(O / o2);
  return mkvec4(P, T1, T2, O == o? gen_1: Qab_trace_init(O, o, P, Po));
}
/* *F a vector of cyclotomic numbers */
static void
compatlift(GEN *F, long o, GEN P)
{
  long i, l;
  GEN f = *F, g = cgetg_copy(f,&l);
  for (i = 1; i < l; i++)
  {
    GEN fi = lift_shallow(gel(f,i));
    gel(g,i) = gmodulo(typ(fi)==t_POL? RgX_inflate(fi,o): fi, P);
  }
  *F = g;
}
static void
chicompatlift(GEN T, GEN *F, GEN *G)
{
  long o1 = itou(gel(T,2)), o2 = itou(gel(T,3));
  GEN P = gel(T,1);
  if (o1 != 1) compatlift(F, o1, P);
  if (o2 != 1 && G) compatlift(G, o2, P);
}
static GEN
chicompatfix(GEN T, GEN F)
{
  GEN V = gel(T,4);
  if (typ(V) == t_VEC) F = gmodulo(QabV_tracerel(V, 0, F), gel(V,1));
  return F;
}

static GEN
c_mul(long n, long d, GEN S)
{
  pari_sp av = avma;
  long nd = n*d;
  GEN F = gel(S,2), G = gel(S,3);
  F = mfcoefs_i(F, nd, 1);
  G = mfcoefs_i(G, nd, 1);
  if (lg(S) == 5) chicompatlift(gel(S,4),&F,&G);
  F = c_deflate(n, d, RgV_mul_RgXn(F,G));
  if (lg(S) == 5) F = chicompatfix(gel(S,4), F);
  return gerepilecopy(av, F);
}
static GEN
c_pow(long n, long d, GEN S)
{
  pari_sp av = avma;
  long nd = n*d;
  GEN F = gel(S,2), a = gel(S,3), f = mfcoefs_i(F,nd,1);
  if (lg(S) == 5) chicompatlift(gel(S,4),&F, NULL);
  f = RgV_pows_RgXn(f, itos(a));
  f = c_deflate(n, d, f);
  if (lg(S) == 5) f = chicompatfix(gel(S,4), f);
  return gerepilecopy(av, f);
}

/* F * Theta */
static GEN
mfmultheta(GEN F)
{
  if (typ(mf_get_gk(F)) == t_FRAC && mf_get_type(F) == t_MF_DIV)
  {
    GEN T = gel(F,3); /* hopefully mfTheta() */
    if (mf_get_type(T) == t_MF_THETA && mf_get_N(T) == 4) return gel(F,2);
  }
  return mfmul(F, mfTheta(NULL));
}

static GEN
c_bracket(long n, long d, GEN S)
{
  pari_sp av = avma;
  long i, nd = n*d;
  GEN F = gel(S,2), G = gel(S,3), tF, tG, C, mpow, res, gk, gl;
  GEN VF = mfcoefs_i(F, nd, 1);
  GEN VG = mfcoefs_i(G, nd, 1);
  ulong j, m = itou(gel(S,4));

  if (!n)
  {
    if (m > 0) { set_avma(av); return mkvec(gen_0); }
    return gerepilecopy(av, mkvec(gmul(gel(VF, 1), gel(VG, 1))));
  }
  tF = cgetg(nd+2, t_VEC);
  tG = cgetg(nd+2, t_VEC);
  res = NULL; gk = mf_get_gk(F); gl = mf_get_gk(G);
  /* pow[i,j+1] = i^j */
  if (lg(S) == 6) chicompatlift(gel(S,5),&VF,&VG);
  mpow = cgetg(m+2, t_MAT);
  gel(mpow,1) = const_col(nd, gen_1);
  for (j = 1; j <= m; j++)
  {
    GEN c = cgetg(nd+1, t_COL);
    gel(mpow,j+1) = c;
    for (i = 1; i <= nd; i++) gel(c,i) = muliu(gcoeff(mpow,i,j), i);
  }
  C = binomial(gaddgs(gk, m-1), m);
  if (odd(m)) C = gneg(C);
  for (j = 0; j <= m; j++)
  { /* C = (-1)^(m-j) binom(m+l-1, j) binom(m+k-1,m-j) */
    GEN c;
    gel(tF,1) = j == 0? gel(VF,1): gen_0;
    gel(tG,1) = j == m? gel(VG,1): gen_0;
    gel(tF,2) = gel(VF,2); /* assume nd >= 1 */
    gel(tG,2) = gel(VG,2);
    for (i = 2; i <= nd; i++)
    {
      gel(tF, i+1) = gmul(gcoeff(mpow,i,j+1),   gel(VF, i+1));
      gel(tG, i+1) = gmul(gcoeff(mpow,i,m-j+1), gel(VG, i+1));
    }
    c = gmul(C, c_deflate(n, d, RgV_mul_RgXn(tF, tG)));
    res = res? gadd(res, c): c;
    if (j < m)
      C = gdiv(gmul(C, gmulsg(m-j, gaddgs(gl,m-j-1))),
               gmulsg(-(j+1), gaddgs(gk,j)));
  }
  if (lg(S) == 6) res = chicompatfix(gel(S,5), res);
  return gerepileupto(av, res);
}
/* linear combination \sum L[j] vecF[j] */
static GEN
c_linear(long n, long d, GEN F, GEN L, GEN dL)
{
  pari_sp av = avma;
  long j, l = lg(L);
  GEN S = NULL;
  for (j = 1; j < l; j++)
  {
    GEN c = gel(L,j);
    if (gequal0(c)) continue;
    c = gmul(c, mfcoefs_i(gel(F,j), n, d));
    S = S? gadd(S,c): c;
  }
  if (!S) return zerovec(n+1);
  if (!is_pm1(dL)) S = gdiv(S, dL);
  return gerepileupto(av, S);
}

/* B_d(T_j Trace^new) as t_MF_BD(t_MF_HECKE(t_MF_NEWTRACE)) or
 * t_MF_HECKE(t_MF_NEWTRACE)
 * or t_MF_NEWTRACE in level N. Set d and j, return t_MF_NEWTRACE component*/
static GEN
bhn_parse(GEN f, long *d, long *j)
{
  long t = mf_get_type(f);
  *d = *j = 1;
  if (t == t_MF_BD) { *d = itos(gel(f,3)); f = gel(f,2); t = mf_get_type(f); }
  if (t == t_MF_HECKE) { *j = gel(f,2)[1]; f = gel(f,3); }
  return f;
}
/* f as above, return the t_MF_NEWTRACE component */
static GEN
bhn_newtrace(GEN f)
{
  long t = mf_get_type(f);
  if (t == t_MF_BD) { f = gel(f,2); t = mf_get_type(f); }
  if (t == t_MF_HECKE) f = gel(f,3);
  return f;
}
static int
ok_bhn_linear(GEN vf)
{
  long i, N0 = 0, l = lg(vf);
  GEN CHI, gk;
  if (l == 1) return 1;
  gk = mf_get_gk(gel(vf,1));
  CHI = mf_get_CHI(gel(vf,1));
  for (i = 1; i < l; i++)
  {
    GEN f = bhn_newtrace(gel(vf,i));
    long N = mf_get_N(f);
    if (mf_get_type(f) != t_MF_NEWTRACE) return 0;
    if (N < N0) return 0; /* largest level must come last */
    N0 = N;
    if (!gequal(gk,mf_get_gk(f))) return 0; /* same k */
    if (!gequal(gel(mf_get_CHI(f),2), gel(CHI,2))) return 0; /* same CHI */
  }
  return 1;
}

/* vF not empty, same hypotheses as bhnmat_extend */
static GEN
bhnmat_extend_nocache(GEN M, long N, long n, long d, GEN vF)
{
  cachenew_t cache;
  long l = lg(vF);
  GEN f;
  if (l == 1) return M? M: cgetg(1, t_MAT);
  f = bhn_newtrace(gel(vF,1)); /* N.B. mf_get_N(f) divides N */
  init_cachenew(&cache, n*d, N, f);
  M = bhnmat_extend(M, n, d, vF, &cache);
  dbg_cachenew(&cache); return M;
}
/* c_linear of "bhn" mf closures, same hypotheses as bhnmat_extend */
static GEN
c_linear_bhn(long n, long d, GEN F)
{
  pari_sp av;
  GEN M, v, vF = gel(F,2), L = gel(F,3), dL = gel(F,4);
  if (lg(L) == 1) return zerovec(n+1);
  av = avma;
  M = bhnmat_extend_nocache(NULL, mf_get_N(F), n, d, vF);
  v = RgM_RgC_mul(M,L); settyp(v, t_VEC);
  if (!is_pm1(dL)) v = gdiv(v, dL);
  return gerepileupto(av, v);
}

/* c in K, K := Q[X]/(T) vz = vector of consecutive powers of root z of T
 * attached to an embedding s: K -> C. Return s(c) in C */
static GEN
Rg_embed1(GEN c, GEN vz)
{
  long t = typ(c);
  if (t == t_POLMOD) { c = gel(c,2); t = typ(c); }
  if (t == t_POL) c = RgX_RgV_eval(c, vz);
  return c;
}
/* return s(P) in C[X] */
static GEN
RgX_embed1(GEN P, GEN vz)
{
  long i, l;
  GEN Q = cgetg_copy(P, &l);
  Q[1] = P[1];
  for (i = 2; i < l; i++) gel(Q,i) = Rg_embed1(gel(P,i), vz);
  return normalizepol_lg(Q,l); /* normally a no-op */
}
/* return s(P) in C^n */
static GEN
vecembed1(GEN P, GEN vz)
{
  long i, l;
  GEN Q = cgetg_copy(P, &l);
  for (i = 1; i < l; i++) gel(Q,i) = Rg_embed1(gel(P,i), vz);
  return Q;
}
/* P in L = K[X]/(U), K = Q[t]/T; s an embedding of K -> C attached
 * to a root of T, extended to an embedding of L -> C attached to a root
 * of s(U); vT powers of the root of T, vU powers of the root of s(U).
 * Return s(P) in C^n */
static GEN
Rg_embed2(GEN P, long vt, GEN vT, GEN vU)
{
  long i, l;
  GEN Q;
  P = liftpol_shallow(P);
  if (typ(P) != t_POL) return P;
  if (varn(P) == vt) return Rg_embed1(P, vT);
  /* varn(P) == vx */
  Q = cgetg_copy(P, &l); Q[1] = P[1];
  for (i = 2; i < l; i++) gel(Q,i) = Rg_embed1(gel(P,i), vT);
  return Rg_embed1(Q, vU);
}
static GEN
vecembed2(GEN P, long vt, GEN vT, GEN vU)
{
  long i, l;
  GEN Q = cgetg_copy(P, &l);
  for (i = 1; i < l; i++) gel(Q,i) = Rg_embed2(gel(P,i), vt, vT, vU);
  return Q;
}
static GEN
RgX_embed2(GEN P, long vt, GEN vT, GEN vU)
{
  long i, l;
  GEN Q = cgetg_copy(P, &l);
  for (i = 2; i < l; i++) gel(Q,i) = Rg_embed2(gel(P,i), vt, vT, vU);
  Q[1] = P[1]; return normalizepol_lg(Q,l);
}
/* embed polynomial f in variable 0 [ may be a scalar ], E from getembed */
static GEN
RgX_embed(GEN f, GEN E)
{
  GEN vT;
  if (typ(f) != t_POL || varn(f) != 0) return mfembed(E, f);
  if (lg(E) == 1) return f;
  vT = gel(E,2);
  if (lg(E) == 3)
    f = RgX_embed1(f, vT);
  else
    f = RgX_embed2(f, varn(gel(E,1)), vT, gel(E,3));
  return f;
}
/* embed vector, E from getembed */
GEN
mfvecembed(GEN E, GEN v)
{
  GEN vT;
  if (lg(E) == 1) return v;
  vT = gel(E,2);
  if (lg(E) == 3)
    v = vecembed1(v, vT);
  else
    v = vecembed2(v, varn(gel(E,1)), vT, gel(E,3));
  return v;
}
GEN
mfmatembed(GEN E, GEN f)
{
  long i, l;
  GEN g;
  if (lg(E) == 1) return f;
  g = cgetg_copy(f, &l);
  for (i = 1; i < l; i++) gel(g,i) = mfvecembed(E, gel(f,i));
  return g;
}
/* embed vector of polynomials in var 0 */
static GEN
RgXV_embed(GEN f, GEN E)
{
  long i, l;
  GEN v;
  if (lg(E) == 1) return f;
  v = cgetg_copy(f, &l);
  for (i = 1; i < l; i++) gel(v,i) = RgX_embed(gel(f,i), E);
  return v;
}

/* embed scalar */
GEN
mfembed(GEN E, GEN f)
{
  GEN vT;
  if (lg(E) == 1) return f;
  vT = gel(E,2);
  if (lg(E) == 3)
    f = Rg_embed1(f, vT);
  else
    f = Rg_embed2(f, varn(gel(E,1)), vT, gel(E,3));
  return f;
}
/* vector of the sigma(f), sigma in vE */
static GEN
RgX_embedall(GEN f, GEN vE)
{
  long i, l = lg(vE);
  GEN v;
  if (l == 2) return RgX_embed(f, gel(vE,1));
  v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(v,i) = RgX_embed(f, gel(vE,i));
  return v;
}
/* matrix whose colums are the sigma(v), sigma in vE */
static GEN
RgC_embedall(GEN v, GEN vE)
{
  long j, l = lg(vE);
  GEN M = cgetg(l, t_MAT);
  for (j = 1; j < l; j++) gel(M,j) = mfvecembed(gel(vE,j), v);
  return M;
}
/* vector of the sigma(v), sigma in vE */
static GEN
Rg_embedall_i(GEN v, GEN vE)
{
  long j, l = lg(vE);
  GEN M = cgetg(l, t_VEC);
  for (j = 1; j < l; j++) gel(M,j) = mfembed(gel(vE,j), v);
  return M;
}
/* vector of the sigma(v), sigma in vE; if #vE == 1, return v */
static GEN
Rg_embedall(GEN v, GEN vE)
{ return (lg(vE) == 2)? mfembed(gel(vE,1), v): Rg_embedall_i(v, vE); }

static GEN
c_div_i(long n, GEN S)
{
  GEN F = gel(S,2), G = gel(S,3);
  GEN a0, a0i, H;
  F = mfcoefs_i(F, n, 1);
  G = mfcoefs_i(G, n, 1);
  if (lg(S) == 5) chicompatlift(gel(S,4),&F,&G);
  F = RgV_to_ser_full(F);
  G = RgV_to_ser_full(G);
  a0 = polcoef_i(G, 0, -1); /* != 0 */
  if (gequal1(a0)) a0 = a0i = NULL;
  else
  {
    a0i = ginv(a0);
    G = gmul(ser_unscale(G,a0), a0i);
    F = gmul(ser_unscale(F,a0), a0i);
  }
  H = gdiv(F, G);
  if (a0) H = ser_unscale(H,a0i);
  H = sertovecslice(H, n);
  if (lg(S) == 5) H = chicompatfix(gel(S,4), H);
  return H;
}
static GEN
c_div(long n, long d, GEN S)
{
  pari_sp av = avma;
  GEN D = (d==1)? c_div_i(n, S): c_deflate(n, d, c_div_i(n*d, S));
  return gerepilecopy(av, D);
}

static GEN
c_shift(long n, long d, GEN F, GEN gsh)
{
  pari_sp av = avma;
  GEN vF;
  long sh = itos(gsh), n1 = n*d + sh;
  if (n1 < 0) return zerovec(n+1);
  vF = mfcoefs_i(F, n1, 1);
  if (sh < 0) vF = shallowconcat(zerovec(-sh), vF);
  else vF = vecslice(vF, sh+1, n1+1);
  return gerepilecopy(av, c_deflate(n, d, vF));
}

static GEN
c_deriv(long n, long d, GEN F, GEN gm)
{
  pari_sp av = avma;
  GEN V = mfcoefs_i(F, n, d), res;
  long i, m = itos(gm);
  if (!m) return V;
  res = cgetg(n+2, t_VEC); gel(res,1) = gen_0;
  if (m < 0)
  { for (i=1; i <= n; i++) gel(res, i+1) = gdiv(gel(V, i+1), powuu(i,-m)); }
  else
  { for (i=1; i <= n; i++) gel(res, i+1) = gmul(gel(V,i+1), powuu(i,m)); }
  return gerepileupto(av, res);
}

static GEN
c_derivE2(long n, long d, GEN F, GEN gm)
{
  pari_sp av = avma;
  GEN VF, VE, res, tmp, gk;
  long i, m = itos(gm), nd;
  if (m == 0) return mfcoefs_i(F, n, d);
  nd = n*d;
  VF = mfcoefs_i(F, nd, 1); VE = mfcoefs_i(mfEk(2), nd, 1);
  gk = mf_get_gk(F);
  if (m == 1)
  {
    res = cgetg(n+2, t_VEC);
    for (i = 0; i <= n; i++) gel(res, i+1) = gmulsg(i, gel(VF, i*d+1));
    tmp = c_deflate(n, d, RgV_mul_RgXn(VF, VE));
    return gerepileupto(av, gsub(res, gmul(gdivgu(gk, 12), tmp)));
  }
  else
  {
    long j;
    for (j = 1; j <= m; j++)
    {
      tmp = RgV_mul_RgXn(VF, VE);
      for (i = 0; i <= nd; i++) gel(VF, i+1) = gmulsg(i, gel(VF, i+1));
      VF = gsub(VF, gmul(gdivgu(gaddgs(gk, 2*(j-1)), 12), tmp));
    }
    return gerepilecopy(av, c_deflate(n, d, VF));
  }
}

/* Twist by the character (D/.) */
static GEN
c_twist(long n, long d, GEN F, GEN D)
{
  pari_sp av = avma;
  GEN v = mfcoefs_i(F, n, d), z = cgetg(n+2, t_VEC);
  long i;
  for (i = 0; i <= n; i++)
  {
    long s;
    GEN a = gel(v, i+1);
    if (d == 1) s = krois(D, i);
    else
    {
      pari_sp av2 = avma;
      s = kronecker(D, muluu(i, d)); set_avma(av2);
    }
    switch(s)
    {
      case 1: a = gcopy(a); break;
      case -1: a = gneg(a); break;
      default: a = gen_0; break;
    }
    gel(z, i+1) = a;
  }
  return gerepileupto(av, z);
}

/* form F given by closure, compute T(n)(F) as closure */
static GEN
c_hecke(long m, long l, GEN DATA, GEN F)
{
  pari_sp av = avma;
  return gerepilecopy(av, hecke_i(m, l, NULL, F, DATA));
}
static GEN
c_const(long n, long d, GEN C)
{
  GEN V = zerovec(n+1);
  long i, j, l = lg(C);
  if (l > d*n+2) l = d*n+2;
  for (i = j = 1; i < l; i+=d, j++) gel(V, j) = gcopy(gel(C,i));
  return V;
}

/* m > 0 */
static GEN
eta3_ZXn(long m)
{
  long l = m+2, n, k;
  GEN P = cgetg(l,t_POL);
  P[1] = evalsigne(1)|evalvarn(0);
  for (n = 2; n < l; n++) gel(P,n) = gen_0;
  for (n = k = 0;; n++)
  {
    if (k + n >= m) { setlg(P, k+3); return P; }
    k += n;
    /* now k = n(n+1) / 2 */
    gel(P, k+2) = odd(n)? utoineg(2*n+1): utoipos(2*n+1);
  }
}

static GEN
c_delta(long n, long d)
{
  pari_sp ltop = avma;
  long N = n*d;
  GEN e;
  if (!N) return mkvec(gen_0);
  e = eta3_ZXn(N);
  e = ZXn_sqr(e,N);
  e = ZXn_sqr(e,N);
  e = ZXn_sqr(e,N); /* eta(x)^24 */
  settyp(e, t_VEC);
  gel(e,1) = gen_0; /* Delta(x) = x*eta(x)^24 as a t_VEC */
  return gerepilecopy(ltop, c_deflate(n, d, e));
}

/* return s(d) such that s|f <=> d | f^2 */
static long
mysqrtu(ulong d)
{
  GEN fa = myfactoru(d), P = gel(fa,1), E = gel(fa,2);
  long l = lg(P), i, s = 1;
  for (i = 1; i < l; i++) s *= upowuu(P[i], (E[i]+1)>>1);
  return s;
}
static GEN
c_theta(long n, long d, GEN psi)
{
  long lim = usqrt(n*d), F = mfcharmodulus(psi), par = mfcharparity(psi);
  long f, d2 = d == 1? 1: mysqrtu(d);
  GEN V = zerovec(n + 1);
  for (f = d2; f <= lim; f += d2)
    if (ugcd(F, f) == 1)
    {
      pari_sp av = avma;
      GEN c = mfchareval(psi, f);
      gel(V, f*f/d + 1) = gerepileupto(av, par < 0? gmulgu(c,2*f): gmul2n(c,1));
    }
  if (F == 1) gel(V, 1) = gen_1;
  return V; /* no gerepile needed */
}

static GEN
c_etaquo(long n, long d, GEN eta, GEN gs)
{
  pari_sp av = avma;
  long s = itos(gs), nd = n*d, nds = nd - s + 1;
  GEN c;
  if (nds <= 0) return zerovec(n+1);
  c = RgX_to_RgC(eta_product_ZXn(eta, nds), nds); settyp(c, t_VEC);
  if (s > 0) c = shallowconcat(zerovec(s), c);
  return gerepilecopy(av, c_deflate(n, d, c));
}

static GEN
c_ell(long n, long d, GEN E)
{
  pari_sp av = avma;
  GEN v;
  if (d == 1) return gconcat(gen_0, ellan(E, n));
  v = vec_prepend(ellan(E, n*d), gen_0);
  return gerepilecopy(av, c_deflate(n, d, v));
}

static GEN
c_cusptrace(long n, long d, GEN F)
{
  pari_sp av = avma;
  GEN D = gel(F,2), res = cgetg(n+2, t_VEC);
  long i, N = mf_get_N(F), k = mf_get_k(F);
  gel(res, 1) = gen_0;
  for (i = 1; i <= n; i++)
    gel(res, i+1) = mfcusptrace_i(N, k, i*d, mydivisorsu(i*d), D);
  return gerepilecopy(av, res);
}

static GEN
c_newtrace(long n, long d, GEN F)
{
  pari_sp av = avma;
  cachenew_t cache;
  long N = mf_get_N(F);
  GEN v;
  init_cachenew(&cache, n == 1? 1: n*d, N, F);
  v = colnewtrace(0, n, d, N, mf_get_k(F), &cache);
  settyp(v, t_VEC); return gerepilecopy(av, v);
}

static GEN
c_Bd(long n, long d, GEN F, GEN A)
{
  pari_sp av = avma;
  long a = itou(A), ad = ugcd(a,d), aad = a/ad, i, j;
  GEN w, v = mfcoefs_i(F, n/aad, d/ad);
  if (a == 1) return v;
  n++; w = zerovec(n);
  for (i = j = 1; j <= n; i++, j += aad) gel(w,j) = gcopy(gel(v,i));
  return gerepileupto(av, w);
}

static GEN
c_dihedral(long n, long d, GEN F)
{
  pari_sp av = avma;
  GEN CHI = mf_get_CHI(F);
  GEN w = gel(F,3), V = dihan(gel(F,2), w, gel(F,4), mfcharorder(CHI), n*d);
  GEN Tinit = gel(w,3), Pm = gel(Tinit,1);
  GEN A = c_deflate(n, d, V);
  if (degpol(Pm) == 1 || RgV_is_ZV(A)) return gerepilecopy(av, A);
  return gerepileupto(av, gmodulo(A, Pm));
}

static GEN
c_mfEH(long n, long d, GEN F)
{
  pari_sp av = avma;
  GEN v, M, A;
  long i, r = mf_get_r(F);
  if (n == 1)
    return gerepilecopy(av, mkvec2(mfEHcoef(r,0),mfEHcoef(r,d)));
  /* speedup mfcoef */
  if (r == 1)
  {
    v = cgetg(n+2, t_VEC);
    gel(v,1) = sstoQ(-1,12);
    for (i = 1; i <= n; i++)
    {
      long id = i*d, a = id & 3;
      gel(v,i+1) = (a==1 || a==2)? gen_0: uutoQ(hclassno6u(id), 6);
    }
    return v; /* no gerepile needed */
  }
  M = mfEHmat(n*d+1,r);
  if (d > 1)
  {
    long l = lg(M);
    for (i = 1; i < l; i++) gel(M,i) = c_deflate(n, d, gel(M,i));
  }
  A = gel(F,2); /* [num(B), den(B)] */
  v = RgC_Rg_div(RgM_RgC_mul(M, gel(A,1)), gel(A,2));
  settyp(v,t_VEC); return gerepileupto(av, v);
}

static GEN
c_mfeisen(long n, long d, GEN F)
{
  pari_sp av = avma;
  GEN v, vchi, E0, P, T, CHI, gk = mf_get_gk(F);
  long i, k;
  if (typ(gk) != t_INT) return c_mfEH(n, d, F);
  k = itou(gk);
  vchi = gel(F,2);
  E0 = gel(vchi,1);
  T = gel(vchi,2);
  P = gel(T,1);
  CHI = gel(vchi,3);
  v = cgetg(n+2, t_VEC);
  gel(v, 1) = gcopy(E0); /* E(0) */
  if (lg(vchi) == 5)
  { /* E_k(chi1,chi2) */
    GEN CHI2 = gel(vchi,4), F3 = gel(F,3);
    long ord = F3[1], j = F3[2];
    for (i = 1; i <= n; i++) gel(v, i+1) = sigchi2(k, CHI, CHI2, i*d, ord);
    v = QabV_tracerel(T, j, v);
  }
  else
  { /* E_k(chi) */
    for (i = 1; i <= n; i++) gel(v, i+1) = sigchi(k, CHI, i*d);
  }
  if (degpol(P) != 1 && !RgV_is_QV(v)) return gerepileupto(av, gmodulo(v, P));
  return gerepilecopy(av, v);
}

/* N^k * (D * B_k)(x/N), set D = denom(B_k) */
static GEN
bern_init(long N, long k, GEN *pD)
{ return ZX_rescale(Q_remove_denom(bernpol(k, 0), pD), utoi(N)); }

/* L(chi_D, 1-k) */
static GEN
lfunquadneg_naive(long D, long k)
{
  GEN B, dS, S;
  long r, N = labs(D);
  pari_sp av;
  if (k == 1 && N == 1) return gneg(ghalf);
  B = bern_init(N, k, &dS);
  dS = mul_denom(dS, stoi(-N*k));
  av = avma;
  for (r = 0, S = gen_0; r < N; r++)
  {
    long c = kross(D, r);
    if (c)
    {
      GEN t = ZX_Z_eval(B, utoi(r));
      S = c > 0 ? addii(S, t) : subii(S, t);
      S = gerepileuptoint(av, S);
    }
  }
  return gdiv(S, dS);
}

/* Returns vector of coeffs from F[0], F[d], ..., F[d*n] */
static GEN
mfcoefs_i(GEN F, long n, long d)
{
  if (n < 0) return gen_0;
  switch(mf_get_type(F))
  {
    case t_MF_CONST: return c_const(n, d, gel(F,2));
    case t_MF_EISEN: return c_mfeisen(n, d, F);
    case t_MF_Ek: return c_Ek(n, d, F);
    case t_MF_DELTA: return c_delta(n, d);
    case t_MF_THETA: return c_theta(n, d, gel(F,2));
    case t_MF_ETAQUO: return c_etaquo(n, d, gel(F,2), gel(F,3));
    case t_MF_ELL: return c_ell(n, d, gel(F,2));
    case t_MF_MUL: return c_mul(n, d, F);
    case t_MF_POW: return c_pow(n, d, F);
    case t_MF_BRACKET: return c_bracket(n, d, F);
    case t_MF_LINEAR: return c_linear(n, d, gel(F,2), gel(F,3), gel(F,4));
    case t_MF_LINEAR_BHN: return c_linear_bhn(n, d, F);
    case t_MF_DIV: return c_div(n, d, F);
    case t_MF_SHIFT: return c_shift(n, d, gel(F,2), gel(F,3));
    case t_MF_DERIV: return c_deriv(n, d, gel(F,2), gel(F,3));
    case t_MF_DERIVE2: return c_derivE2(n, d, gel(F,2), gel(F,3));
    case t_MF_TWIST: return c_twist(n, d, gel(F,2), gel(F,3));
    case t_MF_HECKE: return c_hecke(n, d, gel(F,2), gel(F,3));
    case t_MF_BD: return c_Bd(n, d, gel(F,2), gel(F,3));
    case t_MF_TRACE: return c_cusptrace(n, d, F);
    case t_MF_NEWTRACE: return c_newtrace(n, d, F);
    case t_MF_DIHEDRAL: return c_dihedral(n, d, F);
    default: pari_err_TYPE("mfcoefs",F); return NULL;/*LCOV_EXCL_LINE*/
  }
}

static GEN
matdeflate(long n, long d, GEN M)
{
  long i, l;
  GEN A;
  /*  if (d == 1) return M; */
  A = cgetg_copy(M,&l);
  for (i = 1; i < l; i++) gel(A,i) = c_deflate(n,d,gel(M,i));
  return A;
}
static int
space_is_cusp(long space) { return space != mf_FULL && space != mf_EISEN; }
/* safe with flraw mf */
static GEN
mfcoefs_mf(GEN mf, long n, long d)
{
  GEN MS, ME, E = MF_get_E(mf), S = MF_get_S(mf), M = MF_get_M(mf);
  long lE = lg(E), lS = lg(S), l = lE+lS-1;

  if (l == 1) return cgetg(1, t_MAT);
  if (typ(M) == t_MAT && lg(M) != 1 && (n+1)*d < nbrows(M))
    return matdeflate(n, d, M); /*cached; lg = 1 is possible from mfinit */
  ME = (lE == 1)? cgetg(1, t_MAT): mfvectomat(E, n, d);
  if (lS == 1)
    MS = cgetg(1, t_MAT);
  else if (mf_get_type(gel(S,1)) == t_MF_DIV) /*k 1/2-integer or k=1 (exotic)*/
    MS = matdeflate(n,d, mflineardivtomat(MF_get_N(mf), S, n*d));
  else if (MF_get_k(mf) == 1) /* k = 1 (dihedral) */
  {
    GEN M = mfvectomat(gmael(S,1,2), n, d);
    long i;
    MS = cgetg(lS, t_MAT);
    for (i = 1; i < lS; i++)
    {
      GEN f = gel(S,i), dc = gel(f,4), c = RgM_RgC_mul(M, gel(f,3));
      if (!equali1(dc)) c = RgC_Rg_div(c,dc);
      gel(MS,i) = c;
    }
  }
  else /* k >= 2 integer */
    MS = bhnmat_extend_nocache(NULL, MF_get_N(mf), n, d, S);
  return shallowconcat(ME,MS);
}
GEN
mfcoefs(GEN F, long n, long d)
{
  if (!checkmf_i(F))
  {
    pari_sp av = avma;
    GEN mf = checkMF_i(F); if (!mf) pari_err_TYPE("mfcoefs", F);
    return gerepilecopy(av, mfcoefs_mf(mf,n,d));
  }
  if (d <= 0) pari_err_DOMAIN("mfcoefs", "d", "<=", gen_0, stoi(d));
  if (n < 0) return cgetg(1, t_VEC);
  return mfcoefs_i(F, n, d);
}

/* assume k >= 0 */
static GEN
mfak_i(GEN F, long k)
{
  if (!k) return gel(mfcoefs_i(F,0,1), 1);
  return gel(mfcoefs_i(F,1,k), 2);
}
GEN
mfcoef(GEN F, long n)
{
  pari_sp av = avma;
  if (!checkmf_i(F)) pari_err_TYPE("mfcoef",F);
  return n < 0? gen_0: gerepilecopy(av, mfak_i(F, n));
}

static GEN
paramconst() { return tagparams(t_MF_CONST, mkNK(1,0,mfchartrivial())); }
static GEN
mftrivial(void) { retmkvec2(paramconst(), cgetg(1,t_VEC)); }
static GEN
mf1(void) { retmkvec2(paramconst(), mkvec(gen_1)); }

/* induce mfchar CHI to G */
static GEN
induce(GEN G, GEN CHI)
{
  GEN o, chi;
  if (typ(CHI) == t_INT) /* Kronecker */
  {
    chi = znchar_quad(G, CHI);
    o = ZV_equal0(chi)? gen_1: gen_2;
    CHI = mkvec4(G,chi,o,cgetg(1,t_VEC));
  }
  else
  {
    if (mfcharmodulus(CHI) == itos(znstar_get_N(G))) return CHI;
    CHI = leafcopy(CHI);
    chi = zncharinduce(gel(CHI,1), gel(CHI,2), G);
    gel(CHI,1) = G;
    gel(CHI,2) = chi;
  }
  return CHI;
}
/* induce mfchar CHI to znstar(N) */
static GEN
induceN(long N, GEN CHI)
{
  if (mfcharmodulus(CHI) != N) CHI = induce(znstar0(utoipos(N),1), CHI);
  return CHI;
}
/* *pCHI1 and *pCHI2 are mfchar, induce to common modulus */
static void
char2(GEN *pCHI1, GEN *pCHI2)
{
  GEN CHI1 = *pCHI1, G1 = gel(CHI1,1), N1 = znstar_get_N(G1);
  GEN CHI2 = *pCHI2, G2 = gel(CHI2,1), N2 = znstar_get_N(G2);
  if (!equalii(N1,N2))
  {
    GEN G, d = gcdii(N1,N2);
    if      (equalii(N2,d)) *pCHI2 = induce(G1, CHI2);
    else if (equalii(N1,d)) *pCHI1 = induce(G2, CHI1);
    else
    {
      if (!equali1(d)) N2 = diviiexact(N2,d);
      G = znstar0(mulii(N1,N2), 1);
      *pCHI1 = induce(G, CHI1);
      *pCHI2 = induce(G, CHI2);
    }
  }
}
/* mfchar or charinit wrt same modulus; outputs a mfchar */
static GEN
mfcharmul_i(GEN CHI1, GEN CHI2)
{
  GEN G = gel(CHI1,1), chi3 = zncharmul(G, gel(CHI1,2), gel(CHI2,2));
  return mfcharGL(G, chi3);
}
/* mfchar or charinit; outputs a mfchar */
static GEN
mfcharmul(GEN CHI1, GEN CHI2)
{
  char2(&CHI1, &CHI2); return mfcharmul_i(CHI1,CHI2);
}
/* mfchar or charinit; outputs a mfchar */
static GEN
mfcharpow(GEN CHI, GEN n)
{
  GEN G, chi;
  G = gel(CHI,1); chi = zncharpow(G, gel(CHI,2), n);
  return mfchartoprimitive(mfcharGL(G, chi), NULL);
}
/* mfchar or charinit wrt same modulus; outputs a mfchar */
static GEN
mfchardiv_i(GEN CHI1, GEN CHI2)
{
  GEN G = gel(CHI1,1), chi3 = znchardiv(G, gel(CHI1,2), gel(CHI2,2));
  return mfcharGL(G, chi3);
}
/* mfchar or charinit; outputs a mfchar */
static GEN
mfchardiv(GEN CHI1, GEN CHI2)
{
  char2(&CHI1, &CHI2); return mfchardiv_i(CHI1,CHI2);
}
static GEN
mfcharconj(GEN CHI)
{
  CHI = leafcopy(CHI);
  gel(CHI,2) = zncharconj(gel(CHI,1), gel(CHI,2));
  return CHI;
}

/* CHI mfchar, assume 4 | N. Multiply CHI by \chi_{-4} */
static GEN
mfchilift(GEN CHI, long N)
{
  CHI = induceN(N, CHI);
  return mfcharmul_i(CHI, induce(gel(CHI,1), stoi(-4)));
}
/* CHI defined mod N, N4 = N/4;
 * if CHI is defined mod N4 return CHI;
 * else if CHI' = CHI*(-4,.) is defined mod N4, return CHI' (primitive)
 * else error */
static GEN
mfcharchiliftprim(GEN CHI, long N4)
{
  long FC = mfcharconductor(CHI);
  GEN CHIP;
  if (N4 % FC == 0) return CHI;
  CHIP = mfchartoprimitive(mfchilift(CHI, N4 << 2), &FC);
  if (N4 % FC) pari_err_TYPE("mfkohnenbasis [incorrect CHI]", CHI);
  return CHIP;
}
/* ensure CHI(-1) = (-1)^k [k integer] or 1 [half-integer], by multiplying
 * by (-4/.) if needed */
static GEN
mfchiadjust(GEN CHI, GEN gk, long N)
{
  long par = mfcharparity(CHI);
  if (typ(gk) == t_INT &&  mpodd(gk)) par = -par;
  return par == 1 ? CHI : mfchilift(CHI, N);
}

static GEN
mfsamefield(GEN T, GEN P, GEN Q)
{
  if (degpol(P) == 1) return Q;
  if (degpol(Q) == 1) return P;
  if (!gequal(P,Q)) pari_err_TYPE("mfsamefield [different fields]",mkvec2(P,Q));
  if (T) err_cyclo();
  return P;
}

GEN
mfmul(GEN f, GEN g)
{
  pari_sp av = avma;
  GEN T, N, K, NK, CHI, CHIf, CHIg;
  if (!checkmf_i(f)) pari_err_TYPE("mfmul",f);
  if (!checkmf_i(g)) pari_err_TYPE("mfmul",g);
  N = lcmii(mf_get_gN(f), mf_get_gN(g));
  K = gadd(mf_get_gk(f), mf_get_gk(g));
  CHIf = mf_get_CHI(f);
  CHIg = mf_get_CHI(g);
  CHI = mfchiadjust(mfcharmul(CHIf,CHIg), K, itos(N));
  T = chicompat(CHI, CHIf, CHIg);
  NK = mkgNK(N, K, CHI, mfsamefield(T, mf_get_field(f), mf_get_field(g)));
  return gerepilecopy(av, T? tag3(t_MF_MUL,NK,f,g,T): tag2(t_MF_MUL,NK,f,g));
}
GEN
mfpow(GEN f, long n)
{
  pari_sp av = avma;
  GEN T, KK, NK, gn, CHI, CHIf;
  if (!checkmf_i(f)) pari_err_TYPE("mfpow",f);
  if (!n) return mf1();
  if (n == 1) return gcopy(f);
  KK = gmulsg(n,mf_get_gk(f));
  gn = stoi(n);
  CHIf = mf_get_CHI(f);
  CHI = mfchiadjust(mfcharpow(CHIf,gn), KK, mf_get_N(f));
  T = chicompat(CHI, CHIf, CHIf);
  NK = mkgNK(mf_get_gN(f), KK, CHI, mf_get_field(f));
  return gerepilecopy(av, T? tag3(t_MF_POW,NK,f,gn,T): tag2(t_MF_POW,NK,f,gn));
}
GEN
mfbracket(GEN f, GEN g, long m)
{
  pari_sp av = avma;
  GEN T, N, K, NK, CHI, CHIf, CHIg;
  if (!checkmf_i(f)) pari_err_TYPE("mfbracket",f);
  if (!checkmf_i(g)) pari_err_TYPE("mfbracket",g);
  if (m < 0) pari_err_TYPE("mfbracket [m<0]",stoi(m));
  K = gaddgs(gadd(mf_get_gk(f), mf_get_gk(g)), 2*m);
  if (gsigne(K) < 0) pari_err_IMPL("mfbracket for this form");
  N = lcmii(mf_get_gN(f), mf_get_gN(g));
  CHIf = mf_get_CHI(f);
  CHIg = mf_get_CHI(g);
  CHI = mfcharmul(CHIf, CHIg);
  CHI = mfchiadjust(CHI, K, itou(N));
  T = chicompat(CHI, CHIf, CHIg);
  NK = mkgNK(N, K, CHI, mfsamefield(T, mf_get_field(f), mf_get_field(g)));
  return gerepilecopy(av, T? tag4(t_MF_BRACKET, NK, f, g, utoi(m), T)
                           : tag3(t_MF_BRACKET, NK, f, g, utoi(m)));
}

/* remove 0 entries in L */
static int
mflinear_strip(GEN *pF, GEN *pL)
{
  pari_sp av = avma;
  GEN F = *pF, L = *pL;
  long i, j, l = lg(L);
  GEN F2 = cgetg(l, t_VEC), L2 = cgetg(l, t_VEC);
  for (i = j = 1; i < l; i++)
  {
    if (gequal0(gel(L,i))) continue;
    gel(F2,j) = gel(F,i);
    gel(L2,j) = gel(L,i); j++;
  }
  if (j == l) set_avma(av);
  else
  {
    setlg(F2,j); *pF = F2;
    setlg(L2,j); *pL = L2;
  }
  return (j > 1);
}
static GEN
taglinear_i(long t, GEN NK, GEN F, GEN L)
{
  GEN dL;
  L = Q_remove_denom(L, &dL); if (!dL) dL = gen_1;
  return tag3(t, NK, F, L, dL);
}
static GEN
taglinear(GEN NK, GEN F, GEN L)
{
  long t = ok_bhn_linear(F)? t_MF_LINEAR_BHN: t_MF_LINEAR;
   return taglinear_i(t, NK, F, L);
}
/* assume F has parameters NK = [N,K,CHI] */
static GEN
mflinear_i(GEN NK, GEN F, GEN L)
{
  if (!mflinear_strip(&F,&L)) return mftrivial();
  return taglinear(NK, F,L);
}
static GEN
mflinear_bhn(GEN mf, GEN L)
{
  long i, l;
  GEN P, NK, F = MF_get_S(mf);
  if (!mflinear_strip(&F,&L)) return mftrivial();
  l = lg(L); P = pol_x(1);
  for (i = 1; i < l; i++)
  {
    GEN c = gel(L,i);
    if (typ(c) == t_POLMOD && varn(gel(c,1)) == 1)
      P = mfsamefield(NULL, P, gel(c,1));
  }
  NK = mkgNK(MF_get_gN(mf), MF_get_gk(mf), MF_get_CHI(mf), P);
  return taglinear_i(t_MF_LINEAR_BHN,  NK, F,L);
}

/* F vector of forms with same weight and character but varying level, return
 * global [N,k,chi,P] */
static GEN
vecmfNK(GEN F)
{
  long i, l = lg(F);
  GEN N, f;
  if (l == 1) return mkNK(1, 0, mfchartrivial());
  f = gel(F,1); N = mf_get_gN(f);
  for (i = 2; i < l; i++) N = lcmii(N, mf_get_gN(gel(F,i)));
  return mkgNK(N, mf_get_gk(f), mf_get_CHI(f), mf_get_field(f));
}
/* do not use mflinear: mflineardivtomat rely on F being constant across the
 * basis where mflinear strips the ones matched by 0 coeffs. Assume k and CHI
 * constant, N is allowed to vary. */
static GEN
vecmflinear(GEN F, GEN C)
{
  long i, t, l = lg(C);
  GEN NK, v = cgetg(l, t_VEC);
  if (l == 1) return v;
  t = ok_bhn_linear(F)? t_MF_LINEAR_BHN: t_MF_LINEAR;
  NK = vecmfNK(F);
  for (i = 1; i < l; i++) gel(v,i) = taglinear_i(t, NK, F, gel(C,i));
  return v;
}
/* vecmflinear(F,C), then divide everything by E, which has valuation 0 */
static GEN
vecmflineardiv0(GEN F, GEN C, GEN E)
{
  GEN v = vecmflinear(F, C);
  long i, l = lg(v);
  if (l == 1) return v;
  gel(v,1) = mfdiv_val(gel(v,1), E, 0);
  for (i = 2; i < l; i++)
  { /* v[i] /= E */
    GEN f = shallowcopy(gel(v,1));
    gel(f,2) = gel(v,i);
    gel(v,i) = f;
  }
  return v;
}

/* Non empty linear combination of linear combinations of same
 * F_j=\sum_i \mu_{i,j}G_i so R = \sum_i (\sum_j(\la_j\mu_{i,j})) G_i */
static GEN
mflinear_linear(GEN F, GEN L, int strip)
{
  long l = lg(F), j;
  GEN vF, M = cgetg(l, t_MAT);
  L = shallowcopy(L);
  for (j = 1; j < l; j++)
  {
    GEN f = gel(F,j), c = gel(f,3), d = gel(f,4);
    if (typ(c) == t_VEC) c = shallowtrans(c);
    if (!isint1(d)) gel(L,j) = gdiv(gel(L,j),d);
    gel(M,j) = c;
  }
  vF = gmael(F,1,2); L = RgM_RgC_mul(M,L);
  if (strip && !mflinear_strip(&vF,&L)) return mftrivial();
  return taglinear(vecmfNK(vF), vF, L);
}
/* F nonempty vector of forms of the form mfdiv(mflinear(B,v), E) where E
 * does not vanish at oo, or mflinear(B,v). Apply mflinear(F, L) */
static GEN
mflineardiv_linear(GEN F, GEN L, int strip)
{
  long l = lg(F), j;
  GEN v, E, f;
  if (lg(L) != l) pari_err_DIM("mflineardiv_linear");
  f = gel(F,1); /* l > 1 */
  if (mf_get_type(f) != t_MF_DIV) return mflinear_linear(F,L,strip);
  E = gel(f,3);
  v = cgetg(l, t_VEC);
  for (j = 1; j < l; j++) { GEN f = gel(F,j); gel(v,j) = gel(f,2); }
  return mfdiv_val(mflinear_linear(v,L,strip), E, 0);
}
static GEN
vecmflineardiv_linear(GEN F, GEN M)
{
  long i, l = lg(M);
  GEN v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(v,i) = mflineardiv_linear(F, gel(M,i), 0);
  return v;
}

static GEN
tobasis(GEN mf, GEN F, GEN L)
{
  if (checkmf_i(L) && mf) return mftobasis(mf, L, 0);
  if (typ(F) != t_VEC) pari_err_TYPE("mflinear",F);
  if (!is_vec_t(typ(L))) pari_err_TYPE("mflinear",L);
  if (lg(L) != lg(F)) pari_err_DIM("mflinear");
  return L;
}
GEN
mflinear(GEN F, GEN L)
{
  pari_sp av = avma;
  GEN G, NK, P, mf = checkMF_i(F), N = NULL, K = NULL, CHI = NULL;
  long i, l;
  if (mf)
  {
    GEN gk = MF_get_gk(mf);
    F = MF_get_basis(F);
    if (typ(gk) != t_INT)
      return gerepilecopy(av, mflineardiv_linear(F, L, 1));
    if (itou(gk) > 1 && space_is_cusp(MF_get_space(mf)))
    {
      L = tobasis(mf, F, L);
      return gerepilecopy(av, mflinear_bhn(mf, L));
    }
  }
  L = tobasis(mf, F, L);
  if (!mflinear_strip(&F,&L)) return mftrivial();

  l = lg(F);
  if (l == 2 && gequal1(gel(L,1))) return gerepilecopy(av, gel(F,1));
  P = pol_x(1);
  for (i = 1; i < l; i++)
  {
    GEN f = gel(F,i), c = gel(L,i), Ni, Ki;
    if (!checkmf_i(f)) pari_err_TYPE("mflinear", f);
    Ni = mf_get_gN(f); N = N? lcmii(N, Ni): Ni;
    Ki = mf_get_gk(f);
    if (!K) K = Ki;
    else if (!gequal(K, Ki))
      pari_err_TYPE("mflinear [different weights]", mkvec2(K,Ki));
    P = mfsamefield(NULL, P, mf_get_field(f));
    if (typ(c) == t_POLMOD && varn(gel(c,1)) == 1)
      P = mfsamefield(NULL, P, gel(c,1));
  }
  G = znstar0(N,1);
  for (i = 1; i < l; i++)
  {
    GEN CHI2 = mf_get_CHI(gel(F,i));
    CHI2 = induce(G, CHI2);
    if (!CHI) CHI = CHI2;
    else if (!gequal(CHI, CHI2))
      pari_err_TYPE("mflinear [different characters]", mkvec2(CHI,CHI2));
  }
  NK = mkgNK(N, K, CHI, P);
  return gerepilecopy(av, taglinear(NK,F,L));
}

GEN
mfshift(GEN F, long sh)
{
  pari_sp av = avma;
  if (!checkmf_i(F)) pari_err_TYPE("mfshift",F);
  return gerepilecopy(av, tag2(t_MF_SHIFT, mf_get_NK(F), F, stoi(sh)));
}
static long
mfval(GEN F)
{
  pari_sp av = avma;
  long i = 0, n, sb;
  GEN gk, gN;
  if (!checkmf_i(F)) pari_err_TYPE("mfval", F);
  gN = mf_get_gN(F);
  gk = mf_get_gk(F);
  sb = mfsturmNgk(itou(gN), gk);
  for (n = 1; n <= sb;)
  {
    GEN v;
    if (n > 0.5*sb) n = sb+1;
    v = mfcoefs_i(F, n, 1);
    for (; i <= n; i++)
      if (!gequal0(gel(v, i+1))) return gc_long(av,i);
    n <<= 1;
  }
  return gc_long(av,-1);
}

GEN
mfdiv_val(GEN f, GEN g, long vg)
{
  GEN T, N, K, NK, CHI, CHIf, CHIg;
  if (vg) { f = mfshift(f,vg); g = mfshift(g,vg); }
  N = lcmii(mf_get_gN(f), mf_get_gN(g));
  K = gsub(mf_get_gk(f), mf_get_gk(g));
  CHIf = mf_get_CHI(f);
  CHIg = mf_get_CHI(g);
  CHI = mfchiadjust(mfchardiv(CHIf, CHIg), K, itos(N));
  T = chicompat(CHI, CHIf, CHIg);
  NK = mkgNK(N, K, CHI, mfsamefield(T, mf_get_field(f), mf_get_field(g)));
  return T? tag3(t_MF_DIV, NK, f, g, T): tag2(t_MF_DIV, NK, f, g);
}
GEN
mfdiv(GEN F, GEN G)
{
  pari_sp av = avma;
  long v = mfval(G);
  if (!checkmf_i(F)) pari_err_TYPE("mfdiv", F);
  if (v < 0 || (v && !gequal0(mfcoefs(F, v-1, 1))))
    pari_err_DOMAIN("mfdiv", "ord(G)", ">", strtoGENstr("ord(F)"),
                    mkvec2(F, G));
  return gerepilecopy(av, mfdiv_val(F, G, v));
}
GEN
mfderiv(GEN F, long m)
{
  pari_sp av = avma;
  GEN NK, gk;
  if (!checkmf_i(F)) pari_err_TYPE("mfderiv",F);
  gk = gaddgs(mf_get_gk(F), 2*m);
  NK = mkgNK(mf_get_gN(F), gk, mf_get_CHI(F), mf_get_field(F));
  return gerepilecopy(av, tag2(t_MF_DERIV, NK, F, stoi(m)));
}
GEN
mfderivE2(GEN F, long m)
{
  pari_sp av = avma;
  GEN NK, gk;
  if (!checkmf_i(F)) pari_err_TYPE("mfderivE2",F);
  if (m < 0) pari_err_DOMAIN("mfderivE2","m","<",gen_0,stoi(m));
  gk = gaddgs(mf_get_gk(F), 2*m);
  NK = mkgNK(mf_get_gN(F), gk, mf_get_CHI(F), mf_get_field(F));
  return gerepilecopy(av, tag2(t_MF_DERIVE2, NK, F, stoi(m)));
}

GEN
mftwist(GEN F, GEN D)
{
  pari_sp av = avma;
  GEN NK, CHI, NT, Da;
  long q;
  if (!checkmf_i(F)) pari_err_TYPE("mftwist", F);
  if (typ(D) != t_INT) pari_err_TYPE("mftwist", D);
  Da = mpabs_shallow(D);
  CHI = mf_get_CHI(F); q = mfcharconductor(CHI);
  NT = glcm(glcm(mf_get_gN(F), mulsi(q, Da)), sqri(Da));
  NK = mkgNK(NT, mf_get_gk(F), CHI, mf_get_field(F));
  return gerepilecopy(av, tag2(t_MF_TWIST, NK, F, D));
}

/***************************************************************/
/*                 Generic cache handling                      */
/***************************************************************/
enum { cache_FACT, cache_DIV, cache_H, cache_D, cache_DIH };
typedef struct {
  const char *name;
  GEN cache;
  ulong minself, maxself;
  void (*init)(long);
  ulong miss, maxmiss;
  long compressed;
} cache;

static void constfact(long lim);
static void constdiv(long lim);
static void consttabh(long lim);
static void consttabdihedral(long lim);
static void constcoredisc(long lim);
static THREAD cache caches[] = {
{ "Factors",  NULL,  50000,    50000, &constfact, 0, 0, 0 },
{ "Divisors", NULL,  50000,    50000, &constdiv, 0, 0, 0 },
{ "H",        NULL, 100000, 10000000, &consttabh, 0, 0, 1 },
{ "CorediscF",NULL, 100000, 10000000, &constcoredisc, 0, 0, 0 },
{ "Dihedral", NULL,   1000,     3000, &consttabdihedral, 0, 0, 0 },
};

static void
cache_reset(long id) { caches[id].miss = caches[id].maxmiss = 0; }
static void
cache_delete(long id) { guncloneNULL(caches[id].cache); }
static void
cache_set(long id, GEN S)
{
  GEN old = caches[id].cache;
  caches[id].cache = gclone(S);
  guncloneNULL(old);
}

/* handle a cache miss: store stats, possibly reset table; return value
 * if (now) cached; return NULL on failure. HACK: some caches contain an
 * ulong where the 0 value is impossible, and return it (typecast to GEN) */
static GEN
cache_get(long id, ulong D)
{
  cache *S = &caches[id];
  const ulong d = S->compressed? D>>1: D;
  ulong max, l;

  if (!S->cache)
  {
    max = maxuu(minuu(D, S->maxself), S->minself);
    S->init(max);
    l = lg(S->cache);
  }
  else
  {
    l = lg(S->cache);
    if (l <= d)
    {
      if (D > S->maxmiss) S->maxmiss = D;
      if (DEBUGLEVEL >= 3)
        err_printf("miss in cache %s: %lu, max = %lu\n",
                   S->name, D, S->maxmiss);
      if (S->miss++ >= 5 && D < S->maxself)
      {
        max = minuu(S->maxself, (long)(S->maxmiss * 1.2));
        if (max <= S->maxself)
        {
          if (DEBUGLEVEL >= 3)
            err_printf("resetting cache %s to %lu\n", S->name, max);
          S->init(max); l = lg(S->cache);
        }
      }
    }
  }
  return (l <= d)? NULL: gel(S->cache, d);
}
static GEN
cache_report(long id)
{
  cache *S = &caches[id];
  GEN v = zerocol(5);
  gel(v,1) = strtoGENstr(S->name);
  if (S->cache)
  {
    gel(v,2) = utoi(lg(S->cache)-1);
    gel(v,3) = utoi(S->miss);
    gel(v,4) = utoi(S->maxmiss);
    gel(v,5) = utoi(gsizebyte(S->cache));
  }
  return v;
}
GEN
getcache(void)
{
  pari_sp av = avma;
  GEN M = cgetg(6, t_MAT);
  gel(M,1) = cache_report(cache_FACT);
  gel(M,2) = cache_report(cache_DIV);
  gel(M,3) = cache_report(cache_H);
  gel(M,4) = cache_report(cache_D);
  gel(M,5) = cache_report(cache_DIH);
  return gerepilecopy(av, shallowtrans(M));
}

void
pari_close_mf(void)
{
  cache_delete(cache_FACT);
  cache_delete(cache_DIV);
  cache_delete(cache_H);
  cache_delete(cache_D);
  cache_delete(cache_DIH);
}

/*************************************************************************/
/* a odd, update local cache (recycle memory) */
static GEN
update_factor_cache(long a, long lim, long *pb)
{
  const long step = 16000; /* even; don't increase this: RAM cache thrashing */
  if (a + 2*step > lim)
    *pb = lim; /* fuse last 2 chunks */
  else
    *pb = a + step;
  return vecfactoroddu_i(a, *pb);
}
/* assume lim < MAX_LONG/8 */
static void
constcoredisc(long lim)
{
  pari_sp av2, av = avma;
  GEN D = caches[cache_D].cache, CACHE = NULL;
  long cachea, cacheb, N, LIM = !D ? 4 : lg(D)-1;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_D);
  D = zero_zv(lim);
  av2 = avma;
  cachea = cacheb = 0;
  for (N = 1; N <= lim; N+=2)
  { /* N odd */
    long i, d, d2;
    GEN F;
    if (N > cacheb)
    {
      set_avma(av2); cachea = N;
      CACHE = update_factor_cache(N, lim, &cacheb);
    }
    F = gel(CACHE, ((N-cachea)>>1)+1); /* factoru(N) */
    D[N] = d = corediscs_fact(F); /* = 3 mod 4 or 4 mod 16 */
    d2 = odd(d)? d<<3: d<<1;
    for (i = 1;;)
    {
      if ((N << i) > lim) break;
      D[N<<i] = d2; i++;
      if ((N << i) > lim) break;
      D[N<<i] = d; i++;
    }
  }
  cache_set(cache_D, D);
  set_avma(av);
}

static void
constfact(long lim)
{
  pari_sp av;
  GEN VFACT = caches[cache_FACT].cache;
  long LIM = VFACT? lg(VFACT)-1: 4;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_FACT); av = avma;
  cache_set(cache_FACT, vecfactoru_i(1,lim)); set_avma(av);
}
static void
constdiv(long lim)
{
  pari_sp av;
  GEN VFACT, VDIV = caches[cache_DIV].cache;
  long N, LIM = VDIV? lg(VDIV)-1: 4;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  constfact(lim);
  VFACT = caches[cache_FACT].cache;
  cache_reset(cache_DIV); av = avma;
  VDIV  = cgetg(lim+1, t_VEC);
  for (N = 1; N <= lim; N++) gel(VDIV,N) = divisorsu_fact(gel(VFACT,N));
  cache_set(cache_DIV, VDIV); set_avma(av);
}

/* n > 1, D = divisors(n); sets L = 2*lambda(n), S = sigma(n) */
static void
lamsig(GEN D, long *pL, long *pS)
{
  pari_sp av = avma;
  long i, l = lg(D), L = 1, S = D[l-1]+1;
  for (i = 2; i < l; i++) /* skip d = 1 */
  {
    long d = D[i], nd = D[l-i]; /* nd = n/d */
    if (d < nd) { L += d; S += d + nd; }
    else
    {
      L <<= 1; if (d == nd) { L += d; S += d; }
      break;
    }
  }
  set_avma(av); *pL = L; *pS = S;
}
/* table of 6 * Hurwitz class numbers D <= lim */
static void
consttabh(long lim)
{
  pari_sp av = avma, av2;
  GEN VHDH0, VDIV, CACHE = NULL;
  GEN VHDH = caches[cache_H].cache;
  long r, N, cachea, cacheb, lim0 = VHDH? lg(VHDH)-1: 2, LIM = lim0 << 1;

  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_H);
  r = lim&3L; if (r) lim += 4-r;
  cache_get(cache_DIV, lim);
  VDIV = caches[cache_DIV].cache;
  VHDH0 = cgetg(lim/2 + 1, t_VECSMALL);
  VHDH0[1] = 2;
  VHDH0[2] = 3;
  for (N = 3; N <= lim0; N++) VHDH0[N] = VHDH[N];
  av2 = avma;
  cachea = cacheb = 0;
  for (N = LIM + 3; N <= lim; N += 4)
  {
    long s = 0, limt = usqrt(N>>2), flsq = 0, ind, t, L, S;
    GEN DN, DN2;
    if (N + 2 >= lg(VDIV))
    { /* use local cache */
      GEN F;
      if (N + 2 > cacheb)
      {
        set_avma(av2); cachea = N;
        CACHE = update_factor_cache(N, lim+2, &cacheb);
      }
      F = gel(CACHE, ((N-cachea)>>1)+1); /* factoru(N) */
      DN = divisorsu_fact(F);
      F = gel(CACHE, ((N-cachea)>>1)+2); /* factoru(N+2) */
      DN2 = divisorsu_fact(F);
    }
    else
    { /* use global cache */
      DN = gel(VDIV,N);
      DN2 = gel(VDIV,N+2);
    }
    ind = N >> 1;
    for (t = 1; t <= limt; t++)
    {
      ind -= (t<<2)-2; /* N/2 - 2t^2 */
      if (ind) s += VHDH0[ind]; else flsq = 1;
    }
    lamsig(DN, &L,&S);
    VHDH0[N >> 1] = 2*S - 3*L - 2*s + flsq;
    s = 0; flsq = 0; limt = (usqrt(N+2) - 1) >> 1;
    ind = (N+1) >> 1;
    for (t = 1; t <= limt; t++)
    {
      ind -= t<<2; /* (N+1)/2 - 2t(t+1) */
      if (ind) s += VHDH0[ind]; else flsq = 1;
    }
    lamsig(DN2, &L,&S);
    VHDH0[(N+1) >> 1] = S - 3*(L >> 1) - s - flsq;
  }
  cache_set(cache_H, VHDH0); set_avma(av);
}

/*************************************************************************/
/* Core functions using factorizations, divisors of class numbers caches */
/* TODO: myfactoru and factorization cache should be exported */
static GEN
myfactoru(long N)
{
  GEN z = cache_get(cache_FACT, N);
  return z? gcopy(z): factoru(N);
}
static GEN
mydivisorsu(long N)
{
  GEN z = cache_get(cache_DIV, N);
  return z? leafcopy(z): divisorsu(N);
}
/* write -n = Df^2, D < 0 fundamental discriminant. Return D, set f. */
static long
mycoredisc2neg(ulong n, long *pf)
{
  ulong m, D = (ulong)cache_get(cache_D, n);
  if (D) { *pf = usqrt(n/D); return -(long)D; }
  m = mycore(n, pf);
  if ((m&3) != 3) { m <<= 2; *pf >>= 1; }
  return (long)-m;
}
/* write n = Df^2, D > 0 fundamental discriminant. Return D, set f. */
static long
mycoredisc2pos(ulong n, long *pf)
{
  ulong m = mycore(n, pf);
  if ((m&3) != 1) { m <<= 2; *pf >>= 1; }
  return (long)m;
}

/* D < 0 fundamental. Return 6*hclassno(-D); faster than quadclassunit up
 * to 5*10^5 or so */
static ulong
hclassno6_count(long D)
{
  ulong a, b, b2, h = 0, d = -D;
  int f = 0;

  if (d > 500000) return 6 * quadclassnos(D);
  /* this part would work with -d non fundamental */
  b = d&1; b2 = (1+d)>>2;
  if (!b)
  {
    for (a=1; a*a<b2; a++)
      if (b2%a == 0) h++;
    f = (a*a==b2); b=2; b2=(4+d)>>2;
  }
  while (b2*3 < d)
  {
    if (b2%b == 0) h++;
    for (a=b+1; a*a < b2; a++)
      if (b2%a == 0) h += 2;
    if (a*a == b2) h++;
    b += 2; b2 = (b*b+d)>>2;
  }
  if (b2*3 == d) return 6*h+2;
  if (f) return 6*h+3;
  return 6*h;
}
/* D0 < 0; 6 * hclassno(-D), using D = D0*F^2 */
static long
hclassno6u_2(long D0, long F)
{
  long h;
  if (F == 1) h = hclassno6_count(D0);
  else
  { /* second chance */
    h = (ulong)cache_get(cache_H, -D0);
    if (!h) h = hclassno6_count(D0);
    h *= uhclassnoF_fact(myfactoru(F), D0);
  }
  return h;
}
/* D > 0; 6 * hclassno(D) (6*Hurwitz). Beware, cached value for D (=0,3 mod 4)
 * is stored at D>>1 */
ulong
hclassno6u(ulong D)
{
  ulong z = (ulong)cache_get(cache_H, D);
  long D0, F;
  if (z) return z;
  D0 = mycoredisc2neg(D, &F);
  return hclassno6u_2(D0,F);
}
/* same as hclassno6u without creating caches */
ulong
hclassno6u_no_cache(ulong D)
{
  cache *S = &caches[cache_H];
  long D0, F;
  if (S->cache)
  {
    const ulong d = D>>1; /* compressed */
    if ((ulong)lg(S->cache) > d) return S->cache[d];
  }
  S = &caches[cache_D];
  if (!S->cache || (ulong)lg(S->cache) <= D) return 0;
  D0 = mycoredisc2neg(D, &F);
  return hclassno6u_2(D0,F);
}
/* same, where the decomposition D = D0*F^2 is already known */
static ulong
hclassno6u_i(ulong D, long D0, long F)
{
  ulong z = (ulong)cache_get(cache_H, D);
  if (z) return z;
  return hclassno6u_2(D0,F);
}

/* D < -4 fundamental, h(D), ordinary class number */
static long
myh(long D)
{
  ulong z = (ulong)cache_get(cache_H, -D);
  return z? z / 6: quadclassnos(D);
}

/*************************************************************************/
/*                          TRACE FORMULAS                               */
/* CHIP primitive, initialize for t_POLMOD output */
static GEN
mfcharinit(GEN CHIP)
{
  long n, o, l, vt, N = mfcharmodulus(CHIP);
  GEN c, v, V, G, Pn;
  if (N == 1) return mkvec2(mkvec(gen_1), pol_x(0));
  G = gel(CHIP,1);
  v = ncharvecexpo(G, znconrey_normalized(G, gel(CHIP,2)));
  l = lg(v); V = cgetg(l, t_VEC);
  o = mfcharorder(CHIP);
  Pn = mfcharpol(CHIP); vt = varn(Pn);
  if (o <= 2)
  {
    for (n = 1; n < l; n++)
    {
      if (v[n] < 0) c = gen_0; else c = v[n]? gen_m1: gen_1;
      gel(V,n) = c;
    }
  }
  else
  {
    for (n = 1; n < l; n++)
    {
      if (v[n] < 0) c = gen_0;
      else
      {
        c = Qab_zeta(v[n], o, vt);
        if (typ(c) == t_POL && lg(c) >= lg(Pn)) c = RgX_rem(c, Pn);
      }
      gel(V,n) = c;
    }
  }
  return mkvec2(V, Pn);
}
static GEN
vchip_lift(GEN VCHI, long x, GEN C)
{
  GEN V = gel(VCHI,1);
  long F = lg(V)-1;
  if (F == 1) return C;
  x %= F;
  if (!x) return C;
  if (x <= 0) x += F;
  return gmul(C, gel(V, x));
}
static long
vchip_FC(GEN VCHI) { return lg(gel(VCHI,1))-1; }
static GEN
vchip_mod(GEN VCHI, GEN S)
{ return (typ(S) == t_POL)? RgX_rem(S, gel(VCHI,2)): S; }
static GEN
vchip_polmod(GEN VCHI, GEN S)
{ return (typ(S) == t_POL)? mkpolmod(S, gel(VCHI,2)): S; }

/* contribution of scalar matrices in dimension formula */
static GEN
A1(long N, long k) { return uutoQ(mypsiu(N)*(k-1), 12); }
static long
ceilA1(long N, long k) { return ceildivuu(mypsiu(N) * (k-1), 12); }

/* sturm bound, slightly larger than dimension */
long
mfsturmNk(long N, long k) { return (mypsiu(N) * k) / 12; }
long
mfsturmNgk(long N, GEN k)
{
  long n,d; Qtoss(k,&n,&d);
  return 1 + (mypsiu(N)*n)/(d == 1? 12: 24);
}
static long
mfsturmmf(GEN F) { return mfsturmNgk(mf_get_N(F), mf_get_gk(F)); }

/* List of all solutions of x^2 + x + 1 = 0 modulo N, x modulo N */
static GEN
sqrtm3modN(long N)
{
  pari_sp av;
  GEN fa, P, E, B, mB, A, Q, T, R, v, gen_m3;
  long l, i, n, ct, fl3 = 0, Ninit;
  if (!odd(N) || (N%9) == 0) return cgetg(1,t_VECSMALL);
  Ninit = N;
  if ((N%3) == 0) { N /= 3; fl3 = 1; }
  fa = myfactoru(N); P = gel(fa, 1); E = gel(fa, 2);
  l = lg(P);
  for (i = 1; i < l; i++)
    if ((P[i]%3) == 2) return cgetg(1,t_VECSMALL);
  A = cgetg(l, t_VECSMALL);
  B = cgetg(l, t_VECSMALL);
  mB= cgetg(l, t_VECSMALL);
  Q = cgetg(l, t_VECSMALL); gen_m3 = utoineg(3);
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    Q[i] = upowuu(p,e);
    B[i] = itou( Zp_sqrt(gen_m3, utoipos(p), e) );
    mB[i]= Q[i] - B[i];
  }
  ct = 1 << (l-1);
  T = ZV_producttree(Q);
  R = ZV_chinesetree(Q,T);
  v = cgetg(ct+1, t_VECSMALL);
  av = avma;
  for (n = 1; n <= ct; n++)
  {
    long m = n-1, r;
    for (i = 1; i < l; i++)
    {
      A[i] = (m&1L)? mB[i]: B[i];
      m >>= 1;
    }
    r = itou( ZV_chinese_tree(A, Q, T, R) );
    if (fl3) while (r%3) r += N;
    set_avma(av); v[n] = odd(r) ? (r-1) >> 1 : (r+Ninit-1) >> 1;
  }
  return v;
}

/* number of elliptic points of order 3 in X0(N) */
static long
nu3(long N)
{
  long i, l;
  GEN P;
  if (!odd(N) || (N%9) == 0) return 0;
  if ((N%3) == 0) N /= 3;
  P = gel(myfactoru(N), 1); l = lg(P);
  for (i = 1; i < l; i++) if ((P[i]%3) == 2) return 0;
  return 1L<<(l-1);
}
/* number of elliptic points of order 2 in X0(N) */
static long
nu2(long N)
{
  long i, l;
  GEN P;
  if ((N&3L) == 0) return 0;
  if (!odd(N)) N >>= 1;
  P = gel(myfactoru(N), 1); l = lg(P);
  for (i = 1; i < l; i++) if ((P[i]&3L) == 3) return 0;
  return 1L<<(l-1);
}

/* contribution of elliptic matrices of order 3 in dimension formula
 * Only depends on CHIP the primitive char attached to CHI */
static GEN
A21(long N, long k, GEN CHI)
{
  GEN res, G, chi, o;
  long a21, i, limx, S;
  if ((N&1L) == 0) return gen_0;
  a21 = k%3 - 1;
  if (!a21) return gen_0;
  if (N <= 3) return sstoQ(a21, 3);
  if (!CHI) return sstoQ(nu3(N) * a21, 3);
  res = sqrtm3modN(N); limx = (N - 1) >> 1;
  G = gel(CHI,1); chi = gel(CHI,2);
  o = gmfcharorder(CHI);
  for (S = 0, i = 1; i < lg(res); i++)
  { /* (x,N) = 1; S += chi(x) + chi(x^2) */
    long x = res[i];
    if (x <= limx)
    { /* CHI(x)=e(c/o), 3rd-root of 1 */
      GEN c = znchareval(G, chi, utoi(x), o);
      if (!signe(c)) S += 2; else S--;
    }
  }
  return sstoQ(a21 * S, 3);
}

/* List of all square roots of -1 modulo N */
static GEN
sqrtm1modN(long N)
{
  pari_sp av;
  GEN fa, P, E, B, mB, A, Q, T, R, v;
  long l, i, n, ct, fleven = 0;
  if ((N&3L) == 0) return cgetg(1,t_VECSMALL);
  if ((N&1L) == 0) { N >>= 1; fleven = 1; }
  fa = myfactoru(N); P = gel(fa,1); E = gel(fa,2);
  l = lg(P);
  for (i = 1; i < l; i++)
    if ((P[i]&3L) == 3) return cgetg(1,t_VECSMALL);
  A = cgetg(l, t_VECSMALL);
  B = cgetg(l, t_VECSMALL);
  mB= cgetg(l, t_VECSMALL);
  Q = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    Q[i] = upowuu(p,e);
    B[i] = itou( Zp_sqrt(gen_m1, utoipos(p), e) );
    mB[i]= Q[i] - B[i];
  }
  ct = 1 << (l-1);
  T = ZV_producttree(Q);
  R = ZV_chinesetree(Q,T);
  v = cgetg(ct+1, t_VECSMALL);
  av = avma;
  for (n = 1; n <= ct; n++)
  {
    long m = n-1, r;
    for (i = 1; i < l; i++)
    {
      A[i] = (m&1L)? mB[i]: B[i];
      m >>= 1;
    }
    r = itou( ZV_chinese_tree(A, Q, T, R) );
    if (fleven && !odd(r)) r += N;
    set_avma(av); v[n] = r;
  }
  return v;
}

/* contribution of elliptic matrices of order 4 in dimension formula.
 * Only depends on CHIP the primitive char attached to CHI */
static GEN
A22(long N, long k, GEN CHI)
{
  GEN G, chi, o, res;
  long S, a22, i, limx, o2;
  if ((N&3L) == 0) return gen_0;
  a22 = (k & 3L) - 1; /* (k % 4) - 1 */
  if (!a22) return gen_0;
  if (N <= 2) return sstoQ(a22, 4);
  if (!CHI) return sstoQ(nu2(N)*a22, 4);
  if (mfcharparity(CHI) == -1) return gen_0;
  res = sqrtm1modN(N); limx = (N - 1) >> 1;
  G = gel(CHI,1); chi = gel(CHI,2);
  o = gmfcharorder(CHI);
  o2 = itou(o)>>1;
  for (S = 0, i = 1; i < lg(res); i++)
  { /* (x,N) = 1, S += real(chi(x)) */
    long x = res[i];
    if (x <= limx)
    { /* CHI(x)=e(c/o), 4th-root of 1 */
      long c = itou( znchareval(G, chi, utoi(x), o) );
      if (!c) S++; else if (c == o2) S--;
    }
  }
  return sstoQ(a22 * S, 2);
}

/* sumdiv(N,d,eulerphi(gcd(d,N/d))) */
static long
nuinf(long N)
{
  GEN fa = myfactoru(N), P = gel(fa,1), E = gel(fa,2);
  long i, t = 1, l = lg(P);
  for (i=1; i<l; i++)
  {
    long p = P[i], e = E[i];
    if (odd(e))
      t *= upowuu(p,e>>1) << 1;
    else
      t *= upowuu(p,(e>>1)-1) * (p+1);
  }
  return t;
}

/* contribution of hyperbolic matrices in dimension formula */
static GEN
A3(long N, long FC)
{
  long i, S, NF, l;
  GEN D;
  if (FC == 1) return uutoQ(nuinf(N),2);
  D = mydivisorsu(N); l = lg(D);
  S = 0; NF = N/FC;
  for (i = 1; i < l; i++)
  {
    long g = ugcd(D[i], D[l-i]);
    if (NF%g == 0) S += myeulerphiu(g);
  }
  return uutoQ(S, 2);
}

/* special contribution in weight 2 in dimension formula */
static long
A4(long k, long FC)
{ return (k==2 && FC==1)? 1: 0; }
/* gcd(x,N) */
static long
myugcd(GEN GCD, ulong x)
{
  ulong N = lg(GCD)-1;
  if (x >= N) x %= N;
  return GCD[x+1];
}
/* 1_{gcd(x,N) = 1} * chi(x), return NULL if 0 */
static GEN
mychicgcd(GEN GCD, GEN VCHI, long x)
{
  long N = lg(GCD)-1;
  if (N == 1) return gen_1;
  x = umodsu(x, N);
  if (GCD[x+1] != 1) return NULL;
  x %= vchip_FC(VCHI); if (!x) return gen_1;
  return gel(gel(VCHI,1), x);
}

/* contribution of scalar matrices to trace formula */
static GEN
TA1(long N, long k, GEN VCHI, GEN GCD, long n)
{
  GEN S;
  ulong m;
  if (!uissquareall(n, &m)) return gen_0;
  if (m == 1) return A1(N,k); /* common */
  S = mychicgcd(GCD, VCHI, m);
  return S? gmul(gmul(powuu(m, k-2), A1(N,k)), S): gen_0;
}

/* All square roots modulo 4N, x modulo 2N, precomputed to accelerate TA2 */
static GEN
mksqr(long N)
{
  pari_sp av = avma;
  long x, N2 = N << 1, N4 = N << 2;
  GEN v = const_vec(N2, cgetg(1, t_VECSMALL));
  gel(v, N2) = mkvecsmall(0); /* x = 0 */
  for (x = 1; x <= N; x++)
  {
    long r = (((x*x - 1)%N4) >> 1) + 1;
    gel(v,r) = vecsmall_append(gel(v,r), x);
  }
  return gerepilecopy(av, v);
}

static GEN
mkgcd(long N)
{
  GEN GCD, d;
  long i, N2;
  if (N == 1) return mkvecsmall(N);
  GCD = cgetg(N + 1, t_VECSMALL);
  d = GCD+1; /* GCD[i+1] = d[i] = gcd(i,N) = gcd(N-i,N), i = 0..N-1 */
  d[0] = N; d[1] = d[N-1] = 1; N2 = N>>1;
  for (i = 2; i <= N2; i++) d[i] = d[N-i] = ugcd(N, i);
  return GCD;
}

/* Table of \sum_{x^2-tx+n=0 mod Ng}chi(x) for all g dividing gcd(N,F),
 * F^2 largest such that (t^2-4n)/F^2=0 or 1 mod 4; t >= 0 */
static GEN
mutglistall(long t, long N, long NF, GEN VCHI, long n, GEN MUP, GEN li, GEN GCD)
{
  long i, lx = lg(li);
  GEN DNF = mydivisorsu(NF), v = zerovec(NF);
  long j, g, lDNF = lg(DNF);
  for (i = 1; i < lx; i++)
  {
    long x = (li[i] + t) >> 1, y, lD;
    GEN D, c = mychicgcd(GCD, VCHI, x);
    if (li[i] && li[i] != N)
    {
      GEN c2 = mychicgcd(GCD, VCHI, t - x);
      if (c2) c = c? gadd(c, c2): c2;
    }
    if (!c) continue;
    y = (x*(x - t) + n) / N; /* exact division */
    D = mydivisorsu(ugcd(labs(y), NF)); lD = lg(D);
    for (j=1; j < lD; j++) { g = D[j]; gel(v,g) = gadd(gel(v,g), c); }
  }
  /* j = 1 corresponds to g = 1, and MUP[1] = 1 */
  for (j=2; j < lDNF; j++) { g = DNF[j]; gel(v,g) = gmulsg(MUP[g], gel(v,g)); }
  return v;
}

/* special case (N,F) = 1: easier */
static GEN
mutg1(long t, long N, GEN VCHI, GEN li, GEN GCD)
{ /* (N,F) = 1 */
  GEN S = NULL;
  long i, lx = lg(li);
  for (i = 1; i < lx; i++)
  {
    long x = (li[i] + t) >> 1;
    GEN c = mychicgcd(GCD, VCHI, x);
    if (c) S = S? gadd(S, c): c;
    if (li[i] && li[i] != N)
    {
      c = mychicgcd(GCD, VCHI, t - x);
      if (c) S = S? gadd(S, c): c;
    }
    if (S && !signe(S)) S = NULL; /* strive hard to add gen_0 */
  }
  return S; /* single value */
}

/* Gegenbauer pol; n > 2, P = \sum_{0<=j<=n/2} (-1)^j (n-j)!/j!(n-2*j)! X^j */
GEN
mfrhopol(long n)
{
#ifdef LONG_IS_64BIT
  const long M = 2642249;
#else
  const long M = 1629;
#endif
  long j, d = n >> 1; /* >= 1 */
  GEN P = cgetg(d + 3, t_POL);

  if (n > M) pari_err_IMPL("mfrhopol for large weight"); /* avoid overflow */
  P[1] = evalvarn(0)|evalsigne(1);
  gel(P,2) = gen_1;
  gel(P,3) = utoineg(n-1); /* j = 1 */
  if (d > 1) gel(P,4) = utoipos(((n-3)*(n-2)) >> 1); /* j = 2 */
  if (d > 2) gel(P,5) = utoineg(((n-5)*(n-4)*(n-3)) / 6); /* j = 3 */
  for (j = 4; j <= d; j++)
    gel(P,j+2) = divis(mulis(gel(P,j+1), (n-2*j+1)*(n-2*j+2)), (n-j+1)*(-j));
  return P;
}

/* polrecip(Q)(t2), assume Q(0) = 1 */
GEN
mfrhopol_u_eval(GEN Q, ulong t2)
{
  GEN T = addiu(gel(Q,3), t2);
  long l = lg(Q), j;
  for (j = 4; j < l; j++) T = addii(gel(Q,j), mului(t2, T));
  return T;
}
GEN
mfrhopol_eval(GEN Q, GEN t2)
{
  long l, j;
  GEN T;
  if (lgefint(t2) == 3) return mfrhopol_u_eval(Q, t2[2]);
  l = lg(Q); T = addii(gel(Q,3), t2);
  for (j = 4; j < l; j++) T = addii(gel(Q,j), mulii(t2, T));
  return T;
}
/* return sh * sqrt(n)^nu * G_nu(t/(2*sqrt(n))) for t != 0
 * else (sh/2) * sqrt(n)^nu * G_nu(0) [ implies nu is even ]
 * G_nu(z) = \sum_{0<=j<=nu/2} (-1)^j (nu-j)!/j!(nu-2*j)! * (2z)^(nu-2*j)) */
static GEN
mfrhopowsimp(GEN Q, GEN sh, long nu, long t, long t2, long n)
{
  GEN T;
  switch (nu)
  {
    case 0: return t? sh: gmul2n(sh,-1);
    case 1: return gmulsg(t, sh);
    case 2: return t? gmulsg(t2 - n, sh): gmul(gmul2n(stoi(-n), -1), sh);
    case 3: return gmul(mulss(t, t2 - 2*n), sh);
    default:
      if (!t) return gmul(gmul2n(gel(Q, lg(Q) - 1), -1), sh);
      T = mfrhopol_u_eval(Q, t2); if (odd(nu)) T = mulsi(t, T);
      return gmul(T, sh);
  }
}

/* contribution of elliptic matrices to trace formula */
static GEN
TA2(long N, long k, GEN VCHI, long n, GEN SQRTS, GEN MUP, GEN GCD)
{
  const long n4 = n << 2, N4 = N << 2, nu = k - 2;
  const long st = (!odd(N) && odd(n)) ? 2 : 1;
  long limt, t;
  GEN S, Q;

  limt = usqrt(n4);
  if (limt*limt == n4) limt--;
  Q = nu > 3 ? ZX_z_unscale(mfrhopol(nu), n) : NULL;
  S = gen_0;
  for (t = odd(k)? st: 0; t <= limt; t += st) /* t^2 < 4n */
  {
    pari_sp av = avma;
    long t2 = t*t, D = n4 - t2, F, D0, NF;
    GEN sh, li;

    li = gel(SQRTS, (umodsu(-D - 1, N4) >> 1) + 1);
    if (lg(li) == 1) continue;
    D0 = mycoredisc2neg(D, &F);
    NF = myugcd(GCD, F);
    if (NF == 1)
    { /* (N,F) = 1 => single value in mutglistall */
      GEN mut = mutg1(t, N, VCHI, li, GCD);
      if (!mut) { set_avma(av); continue; }
      sh = gmul(uutoQ(hclassno6u_i(D,D0,F),6), mut);
    }
    else
    {
      GEN v = mutglistall(t, N, NF, VCHI, n, MUP, li, GCD);
      GEN DF = mydivisorsu(F);
      long i, lDF = lg(DF);
      sh = gen_0;
      for (i = 1; i < lDF; i++)
      {
        long Ff, f = DF[i], g = myugcd(GCD, f);
        GEN mut = gel(v, g);
        if (gequal0(mut)) continue;
        Ff = DF[lDF-i]; /* F/f */
        if (Ff == 1) sh = gadd(sh, mut);
        else
        {
          GEN P = gel(myfactoru(Ff), 1);
          long j, lP = lg(P);
          for (j = 1; j < lP; j++) { long p = P[j]; Ff -= kross(D0, p)*Ff/p; }
          sh = gadd(sh, gmulsg(Ff, mut));
        }
      }
      if (gequal0(sh)) { set_avma(av); continue; }
      if (D0 == -3) sh = gdivgu(sh, 3);
      else if (D0 == -4) sh = gdivgu(sh, 2);
      else sh = gmulgu(sh, myh(D0));
    }
    S = gerepileupto(av, gadd(S, mfrhopowsimp(Q,sh,nu,t,t2,n)));
  }
  return S;
}

/* compute global auxiliary data for TA3 */
static GEN
mkbez(long N, long FC)
{
  long ct, i, NF = N/FC;
  GEN w, D = mydivisorsu(N);
  long l = lg(D);

  w = cgetg(l, t_VEC);
  for (i = ct = 1; i < l; i++)
  {
    long u, v, h, c = D[i], Nc = D[l-i];
    if (c > Nc) break;
    h = cbezout(c, Nc, &u, &v);
    if (h == 1) /* shortcut */
      gel(w, ct++) = mkvecsmall4(1,u*c,1,i);
    else if (!(NF%h))
      gel(w, ct++) = mkvecsmall4(h,u*(c/h),myeulerphiu(h),i);
  }
  setlg(w,ct); stackdummy((pari_sp)(w+ct),(pari_sp)(w+l));
  return w;
}

/* contribution of hyperbolic matrices to trace formula, d * nd = n,
 * DN = divisorsu(N) */
static GEN
auxsum(GEN VCHI, GEN GCD, long d, long nd, GEN DN, GEN BEZ)
{
  GEN S = gen_0;
  long ct, g = nd - d, lDN = lg(DN), lBEZ = lg(BEZ);
  for (ct = 1; ct < lBEZ; ct++)
  {
    GEN y, B = gel(BEZ, ct);
    long ic, c, Nc, uch, h = B[1];
    if (g%h) continue;
    uch = B[2];
    ic  = B[4];
    c = DN[ic];
    Nc= DN[lDN - ic]; /* Nc = N/c */
    if (ugcd(Nc, nd) == 1)
      y = mychicgcd(GCD, VCHI, d + uch*g); /* 0 if (c,d) > 1 */
    else
      y = NULL;
    if (c != Nc && ugcd(Nc, d) == 1)
    {
      GEN y2 = mychicgcd(GCD, VCHI, nd - uch*g); /* 0 if (c,nd) > 1 */
      if (y2) y = y? gadd(y, y2): y2;
    }
    if (y) S = gadd(S, gmulsg(B[3], y));
  }
  return S;
}

static GEN
TA3(long N, long k, GEN VCHI, GEN GCD, GEN Dn, GEN BEZ)
{
  GEN S = gen_0, DN = mydivisorsu(N);
  long i, l = lg(Dn);
  for (i = 1; i < l; i++)
  {
    long d = Dn[i], nd = Dn[l-i]; /* = n/d */
    GEN t, u;
    if (d > nd) break;
    t = auxsum(VCHI, GCD, d, nd, DN, BEZ);
    if (isintzero(t)) continue;
    u = powuu(d,k-1); if (d == nd) u = gmul2n(u,-1);
    S = gadd(S, gmul(u,t));
  }
  return S;
}

/* special contribution in weight 2 in trace formula */
static long
TA4(long k, GEN VCHIP, GEN Dn, GEN GCD)
{
  long i, l, S;
  if (k != 2 || vchip_FC(VCHIP) != 1) return 0;
  l = lg(Dn); S = 0;
  for (i = 1; i < l; i++)
  {
    long d = Dn[i]; /* gcd(N,n/d) == 1? */
    if (myugcd(GCD, Dn[l-i]) == 1) S += d;
  }
  return S;
}

/* precomputation of products occurring im mutg, again to accelerate TA2 */
static GEN
mkmup(long N)
{
  GEN fa = myfactoru(N), P = gel(fa,1), D = divisorsu_fact(fa);
  long i, lP = lg(P), lD = lg(D);
  GEN MUP = zero_zv(N);
  MUP[1] = 1;
  for (i = 2; i < lD; i++)
  {
    long j, g = D[i], Ng = D[lD-i]; /*  N/g */
    for (j = 1; j < lP; j++) { long p = P[j]; if (Ng%p) g += g/p; }
    MUP[D[i]] = g;
  }
  return MUP;
}

/* quadratic nonresidues mod p; p odd prime, p^2 fits in a long */
static GEN
non_residues(long p)
{
  long i, j, p2 = p >> 1;
  GEN v = cgetg(p2+1, t_VECSMALL), w = const_vecsmall(p-1, 1);
  for (i = 2; i <= p2; i++) w[(i*i) % p] = 0; /* no need to check 1 */
  for (i = 2, j = 1; i < p; i++) if (w[i]) v[j++] = i;
  return v;
}

/* CHIP primitive. Return t_VECSMALL v of length q such that
 * Tr^new_{N,CHIP}(n) = 0 whenever v[(n%q) + 1] is nonzero */
static GEN
mfnewzerodata(long N, GEN CHIP)
{
  GEN V, M, L, faN = myfactoru(N), PN = gel(faN,1), EN = gel(faN,2);
  GEN G = gel(CHIP,1), chi = gel(CHIP,2);
  GEN fa = znstar_get_faN(G), P = ZV_to_zv(gel(fa,1)), E = gel(fa,2);
  long i, mod, j = 1, l = lg(PN);

  M = cgetg(l, t_VECSMALL); M[1] = 0;
  V = cgetg(l, t_VEC);
  /* Tr^new(n) = 0 if (n mod M[i]) in V[i]  */
  if ((N & 3) == 0)
  {
    long e = EN[1];
    long c = (lg(P) > 1 && P[1] == 2)? E[1]: 0; /* c = v_2(FC) */
    /* e >= 2 */
    if (c == e-1) return NULL; /* Tr^new = 0 */
    if (c == e)
    {
      if (e == 2)
      { /* sc: -4 */
        gel(V,1) = mkvecsmall(3);
        M[1] = 4;
      }
      else if (e == 3)
      { /* sc: -8 (CHI_2(-1)=-1<=>chi[1]=1) and 8 (CHI_2(-1)=1 <=> chi[1]=0) */
        long t = signe(gel(chi,1))? 7: 3;
        gel(V,1) = mkvecsmall2(5, t);
        M[1] = 8;
      }
    }
    else if (e == 5 && c == 3)
    { /* sc: -8 (CHI_2(-1)=-1<=>chi[1]=1) and 8 (CHI_2(-1)=1 <=> chi[1]=0) */
      long t = signe(gel(chi,1))? 7: 3;
      gel(V,1) = mkvecsmalln(6, 2L,4L,5L,6L,8L,t);
      M[1] = 8;
    }
    else if ((e == 4 && c == 2) || (e == 5 && c <= 2) || (e == 6 && c <= 2)
         || (e >= 7 && c == e - 3))
    { /* sc: 4 */
      gel(V,1) = mkvecsmall3(0,2,3);
      M[1] = 4;
    }
    else if ((e <= 4 && c == 0) || (e >= 5 && c == e - 2))
    { /* sc: 2 */
      gel(V,1) = mkvecsmall(0);
      M[1] = 2;
    }
    else if ((e == 6 && c == 3) || (e >= 7 && c <= e - 4))
    { /* sc: -2 */
      gel(V,1) = mkvecsmalln(7, 0L,2L,3L,4L,5L,6L,7L);
      M[1] = 8;
    }
  }
  j = M[1]? 2: 1;
  for (i = odd(N)? 1: 2; i < l; i++) /* skip p=2, done above */
  {
    long p = PN[i], e = EN[i];
    long z = zv_search(P, p), c = z? E[z]: 0; /* c = v_p(FC) */
    if ((e <= 2 && c == 1 && itos(gel(chi,z)) == (p>>1)) /* ord(CHI_p)=2 */
        || (e >= 3 && c <= e - 2))
    { /* sc: -p */
      GEN v = non_residues(p);
      if (e != 1) v = vecsmall_prepend(v, 0);
      gel(V,j) = v;
      M[j] = p; j++;
    }
    else if (e >= 2 && c < e)
    { /* sc: p */
      gel(V,j) = mkvecsmall(0);
      M[j] = p; j++;
    }
  }
  if (j == 1) return cgetg(1, t_VECSMALL);
  setlg(V,j); setlg(M,j); mod = zv_prod(M);
  L = zero_zv(mod);
  for (i = 1; i < j; i++)
  {
    GEN v = gel(V,i);
    long s, m = M[i], lv = lg(v);
    for (s = 1; s < lv; s++)
    {
      long a = v[s] + 1;
      do { L[a] = 1; a += m; } while (a <= mod);
    }
  }
  return L;
}
/* v=mfnewzerodata(N,CHI); returns TRUE if newtrace(n) must be zero,
 * (but newtrace(n) may still be zero if we return FALSE) */
static long
mfnewchkzero(GEN v, long n) { long q = lg(v)-1; return q && v[(n%q) + 1]; }

/* if (!VCHIP): from mftraceform_cusp;
 * else from initnewtrace and CHI is known to be primitive */
static GEN
inittrace(long N, GEN CHI, GEN VCHIP)
{
  long FC;
  if (VCHIP)
    FC = mfcharmodulus(CHI);
  else
    VCHIP = mfcharinit(mfchartoprimitive(CHI, &FC));
  return mkvecn(5, mksqr(N), mkmup(N), mkgcd(N), VCHIP, mkbez(N, FC));
}

/* p > 2 prime; return a sorted t_VECSMALL of primes s.t Tr^new(p) = 0 for all
 * weights > 2 */
static GEN
inittrconj(long N, long FC)
{
  GEN fa, P, E, v;
  long i, k, l;

  if (FC != 1) return cgetg(1,t_VECSMALL);

  fa = myfactoru(N >> vals(N));
  P = gel(fa,1); l = lg(P);
  E = gel(fa,2);
  v = cgetg(l, t_VECSMALL);
  for (i = k = 1; i < l; i++)
  {
    long j, p = P[i]; /* > 2 */
    for (j = 1; j < l; j++)
      if (j != i && E[j] == 1 && kross(-p, P[j]) == 1) v[k++] = p;
  }
  setlg(v,k); return v;
}

/* assume CHIP primitive, f(CHIP) | N; NZ = mfnewzerodata(N,CHIP) */
static GEN
initnewtrace_i(long N, GEN CHIP, GEN NZ)
{
  GEN T = const_vec(N, cgetg(1,t_VEC)), D, VCHIP;
  long FC = mfcharmodulus(CHIP), N1, N2, i, l;

  if (!NZ) NZ = mkvecsmall(1); /*Tr^new = 0; initialize data nevertheless*/
  VCHIP = mfcharinit(CHIP);
  N1 = N/FC; newd_params(N1, &N2);
  D = mydivisorsu(N1/N2); l = lg(D);
  N2 *= FC;
  for (i = 1; i < l; i++)
  {
    long M = D[i]*N2;
    gel(T,M) = inittrace(M, CHIP, VCHIP);
  }
  gel(T,N) = shallowconcat(gel(T,N), mkvec2(NZ, inittrconj(N,FC)));
  return T;
}
/* don't initialize if Tr^new = 0, return NULL */
static GEN
initnewtrace(long N, GEN CHI)
{
  GEN CHIP = mfchartoprimitive(CHI, NULL), NZ = mfnewzerodata(N,CHIP);
  return NZ? initnewtrace_i(N, CHIP, NZ): NULL;
}

/* (-1)^k */
static long
m1pk(long k) { return odd(k)? -1 : 1; }
static long
badchar(long N, long k, GEN CHI)
{ return mfcharparity(CHI) != m1pk(k) || (CHI && N % mfcharconductor(CHI)); }


static long
mfcuspdim_i(long N, long k, GEN CHI, GEN vSP)
{
  pari_sp av = avma;
  long FC;
  GEN s;
  if (k <= 0) return 0;
  if (k == 1) return CHI? mf1cuspdim(N, CHI, vSP): 0;
  FC = CHI? mfcharconductor(CHI): 1;
  if (FC == 1) CHI = NULL;
  s = gsub(A1(N, k), gadd(A21(N, k, CHI), A22(N, k, CHI)));
  s = gadd(s, gsubsg(A4(k, FC), A3(N, FC)));
  return gc_long(av, itos(s));
}
/* dimension of space of cusp forms S_k(\G_0(N),CHI)
 * Only depends on CHIP the primitive char attached to CHI */
long
mfcuspdim(long N, long k, GEN CHI) { return mfcuspdim_i(N, k, CHI, NULL); }

/* dimension of whole space M_k(\G_0(N),CHI)
 * Only depends on CHIP the primitive char attached to CHI; assumes !badchar */
long
mffulldim(long N, long k, GEN CHI)
{
  pari_sp av = avma;
  long FC = CHI? mfcharconductor(CHI): 1;
  GEN s;
  if (k <= 0) return (k == 0 && FC == 1)? 1: 0;
  if (k == 1) return gc_long(av, itos(A3(N, FC)) + mf1cuspdim(N, CHI, NULL));
  if (FC == 1) CHI = NULL;
  s = gsub(A1(N, k), gadd(A21(N, k, CHI), A22(N, k, CHI)));
  s = gadd(s, A3(N, FC));
  return gc_long(av, itos(s));
}

/* Dimension of the space of Eisenstein series */
long
mfeisensteindim(long N, long k, GEN CHI)
{
  pari_sp av = avma;
  long s, FC = CHI? mfcharconductor(CHI): 1;
  if (k <= 0) return (k == 0 && FC == 1)? 1: 0;
  s = itos(gmul2n(A3(N, FC), 1));
  if (k > 1) s -= A4(k, FC); else s >>= 1;
  return gc_long(av,s);
}

enum { _SQRTS = 1, _MUP, _GCD, _VCHIP, _BEZ, _NEWLZ, _TRCONJ };
/* Trace of T(n) on space of cuspforms; only depends on CHIP the primitive char
 * attached to CHI */
static GEN
mfcusptrace_i(long N, long k, long n, GEN Dn, GEN S)
{
  pari_sp av = avma;
  GEN a, b, VCHIP, GCD;
  long t;
  if (!n) return gen_0;
  VCHIP = gel(S,_VCHIP);
  GCD = gel(S,_GCD);
  t = TA4(k, VCHIP, Dn, GCD);
  a = TA1(N, k, VCHIP, GCD, n); if (t) a = gaddgs(a,t);
  b = TA2(N, k, VCHIP, n, gel(S,_SQRTS), gel(S,_MUP), GCD);
  b = gadd(b, TA3(N, k, VCHIP, GCD, Dn, gel(S,_BEZ)));
  b = gsub(a,b);
  if (typ(b) != t_POL) return gerepileupto(av, b);
  return gerepilecopy(av, vchip_polmod(VCHIP, b));
}

static GEN
mfcusptracecache(long N, long k, long n, GEN Dn, GEN S, cachenew_t *cache)
{
  GEN C = NULL, T = gel(cache->vfull,N);
  long lcache = lg(T);
  if (n < lcache) C = gel(T, n);
  if (C) cache->cuspHIT++; else C = mfcusptrace_i(N, k, n, Dn, S);
  cache->cuspTOTAL++;
  if (n < lcache) gel(T,n) = C;
  return C;
}

/* return the divisors of n, known to be among the elements of D */
static GEN
div_restrict(GEN D, ulong n)
{
  long i, j, l;
  GEN v, VDIV = caches[cache_DIV].cache;
  if (lg(VDIV) > n) return gel(VDIV,n);
  l = lg(D);
  v = cgetg(l, t_VECSMALL);
  for (i = j = 1; i < l; i++)
  {
    ulong d = D[i];
    if (n % d == 0) v[j++] = d;
  }
  setlg(v,j); return v;
}

/* for some prime divisors of N, Tr^new(p) = 0 */
static int
trconj(GEN T, long N, long n)
{ return (lg(T) > 1 && N % n == 0 && zv_search(T, n)); }

/* n > 0; trace formula on new space */
static GEN
mfnewtrace_i(long N, long k, long n, cachenew_t *cache)
{
  GEN VCHIP, s, Dn, DN1, SN, S = cache->DATA;
  long FC, N1, N2, N1N2, g, i, j, lDN1;

  if (!S) return gen_0;
  SN = gel(S,N);
  if (mfnewchkzero(gel(SN,_NEWLZ), n)) return gen_0;
  if (k > 2 && trconj(gel(SN,_TRCONJ), N, n)) return gen_0;
  VCHIP = gel(SN, _VCHIP); FC = vchip_FC(VCHIP);
  N1 = N/FC; newt_params(N1, n, FC, &g, &N2);
  N1N2 = N1/N2;
  DN1 = mydivisorsu(N1N2); lDN1 = lg(DN1);
  N2 *= FC;
  Dn = mydivisorsu(n); /* this one is probably out of cache */
  s = gmulsg(mubeta2(N1N2,n), mfcusptracecache(N2, k, n, Dn, gel(S,N2), cache));
  for (i = 2; i < lDN1; i++)
  { /* skip M1 = 1, done above */
    long M1 = DN1[i], N1M1 = DN1[lDN1-i];
    GEN Dg = mydivisorsu(ugcd(M1, g));
    M1 *= N2;
    s = gadd(s, gmulsg(mubeta2(N1M1,n),
                       mfcusptracecache(M1, k, n, Dn, gel(S,M1), cache)));
    for (j = 2; j < lg(Dg); j++) /* skip d = 1, done above */
    {
      long d = Dg[j], ndd = n/(d*d), M = M1/d;
      GEN z = mulsi(mubeta2(N1M1,ndd), powuu(d,k-1)), C = vchip_lift(VCHIP,d,z);
      GEN Dndd = div_restrict(Dn, ndd);
      s = gadd(s, gmul(C, mfcusptracecache(M, k, ndd, Dndd, gel(S,M), cache)));
    }
    s = vchip_mod(VCHIP, s);
  }
  return vchip_polmod(VCHIP, s);
}

static GEN
get_DIH(long N)
{
  GEN x = cache_get(cache_DIH, N);
  return x? gcopy(x): mfdihedral(N);
}
static GEN
get_vDIH(long N, GEN D)
{
  GEN x = const_vec(N, NULL);
  long i, l;
  if (!D) D = mydivisorsu(N);
  l = lg(D);
  for (i = 1; i < l; i++) { long d = D[i]; gel(x, d) = get_DIH(d); }
  return x;
}

/* divisors of N which are multiple of F */
static GEN
divisorsNF(long N, long F)
{
  GEN D = mydivisorsu(N / F);
  long l = lg(D), i;
  for (i = 1; i < l; i++) D[i] = N / D[i];
  return D;
}
/* mfcuspdim(N,k,CHI) - mfnewdim(N,k,CHI); CHIP primitive (for efficiency) */
static long
mfolddim_i(long N, long k, GEN CHIP, GEN vSP)
{
  long S, i, l, F = mfcharmodulus(CHIP), N1 = N / F, N2;
  GEN D;
  newd_params(N1, &N2); /* will ensure mubeta != 0 */
  D = mydivisorsu(N1/N2); l = lg(D); S = 0;
  if (k == 1 && !vSP) vSP = get_vDIH(N, divisorsNF(N, F));
  for (i = 2; i < l; i++)
  {
    long d = mfcuspdim_i(N / D[i], k, CHIP, vSP);
    if (d) S -= mubeta(D[i]) * d;
  }
  return S;
}
long
mfolddim(long N, long k, GEN CHI)
{
  pari_sp av = avma;
  GEN CHIP = mfchartoprimitive(CHI, NULL);
  return gc_long(av, mfolddim_i(N, k, CHIP, NULL));
}
/* Only depends on CHIP the primitive char attached to CHI; assumes !badchar */
long
mfnewdim(long N, long k, GEN CHI)
{
  pari_sp av;
  long S, F;
  GEN vSP, CHIP = mfchartoprimitive(CHI, &F);
  vSP = (k == 1)? get_vDIH(N, divisorsNF(N, F)): NULL;
  S = mfcuspdim_i(N, k, CHIP, vSP); if (!S) return 0;
  av = avma; return gc_long(av, S - mfolddim_i(N, k, CHIP, vSP));
}

/* trace form, given as closure */
static GEN
mftraceform_new(long N, long k, GEN CHI)
{
  GEN T;
  if (k == 1) return initwt1newtrace(mfinit_Nkchi(N, 1, CHI, mf_CUSP, 0));
  T = initnewtrace(N,CHI); if (!T) return mftrivial();
  return tag(t_MF_NEWTRACE, mkNK(N,k,CHI), T);
}
static GEN
mftraceform_cusp(long N, long k, GEN CHI)
{
  if (k == 1) return initwt1trace(mfinit_Nkchi(N, 1, CHI, mf_CUSP, 0));
  return tag(t_MF_TRACE, mkNK(N,k,CHI), inittrace(N,CHI,NULL));
}
static GEN
mftraceform_i(GEN NK, long space)
{
  GEN CHI;
  long N, k;
  checkNK(NK, &N, &k, &CHI, 0);
  if (!mfdim_Nkchi(N, k, CHI, space)) return mftrivial();
  switch(space)
  {
    case mf_NEW: return mftraceform_new(N, k, CHI);
    case mf_CUSP:return mftraceform_cusp(N, k, CHI);
  }
  pari_err_DOMAIN("mftraceform", "space", "=", utoi(space), NK);
  return NULL;/*LCOV_EXCL_LINE*/
}
GEN
mftraceform(GEN NK, long space)
{ pari_sp av = avma; return gerepilecopy(av, mftraceform_i(NK,space)); }

static GEN
hecke_data(long N, long n)
{ return mkvecsmall3(n, u_ppo(n, N), N); }
/* 1/2-integral weight */
static GEN
heckef2_data(long N, long n)
{
  ulong f, fN, fN2;
  if (!uissquareall(n, &f)) return NULL;
  fN = u_ppo(f, N); fN2 = fN*fN;
  return mkvec2(myfactoru(fN), mkvecsmall4(n, N, fN2, n/fN2));
}
/* N = mf_get_N(F) or a multiple */
static GEN
mfhecke_i(long n, long N, GEN F)
{
  if (n == 1) return F;
  return tag2(t_MF_HECKE, mf_get_NK(F), hecke_data(N,n), F);
}

GEN
mfhecke(GEN mf, GEN F, long n)
{
  pari_sp av = avma;
  GEN NK, CHI, gk, DATA;
  long N, nk, dk;
  mf = checkMF(mf);
  if (!checkmf_i(F)) pari_err_TYPE("mfhecke",F);
  if (n <= 0) pari_err_TYPE("mfhecke [n <= 0]", stoi(n));
  if (n == 1) return gcopy(F);
  gk = mf_get_gk(F);
  Qtoss(gk,&nk,&dk);
  CHI = mf_get_CHI(F);
  N = MF_get_N(mf);
  if (dk == 2)
  {
    DATA = heckef2_data(N,n);
    if (!DATA) return mftrivial();
  }
  else
    DATA = hecke_data(N,n);
  NK = mkgNK(lcmii(stoi(N), mf_get_gN(F)), gk, CHI, mf_get_field(F));
  return gerepilecopy(av, tag2(t_MF_HECKE, NK, DATA, F));
}

/* form F given by closure, compute B(d)(F) as closure (q -> q^d) */
static GEN
mfbd_i(GEN F, long d)
{
  GEN D, NK, gk, CHI;
  if (d == 1) return F;
  if (d <= 0) pari_err_TYPE("mfbd [d <= 0]", stoi(d));
  if (mf_get_type(F) != t_MF_BD) D = utoi(d);
  else { D = mului(d, gel(F,3)); F = gel(F,2); }
  gk = mf_get_gk(F); CHI = mf_get_CHI(F);
  if (typ(gk) != t_INT) CHI = mfcharmul(CHI, get_mfchar(utoi(d << 2)));
  NK = mkgNK(muliu(mf_get_gN(F), d), gk, CHI, mf_get_field(F));
  return tag2(t_MF_BD, NK, F, D);
}
GEN
mfbd(GEN F, long d)
{
  pari_sp av = avma;
  if (!checkmf_i(F)) pari_err_TYPE("mfbd",F);
  return gerepilecopy(av, mfbd_i(F, d));
}

/* A[i+1] = a(t*i^2) */
static GEN
RgV_shimura(GEN A, long n, long t, long N, long r, GEN CHI)
{
  GEN R, a0, Pn = mfcharpol(CHI);
  long m, st, ord = mfcharorder(CHI), vt = varn(Pn), Nt = t == 1? N: ulcm(N,t);

  R = cgetg(n + 2, t_VEC);
  st = odd(r)? -t: t;
  a0 = gel(A, 1);
  if (!gequal0(a0))
  {
    long o = mfcharorder(CHI);
    if (st != 1 && odd(o)) o <<= 1;
    a0 = gmul(a0, charLFwtk(Nt, r, CHI, o, st));
  }
  gel(R, 1) = a0;
  for (m = 1; m <= n; m++)
  {
    GEN Dm = mydivisorsu(u_ppo(m, Nt)), S = gel(A, m*m + 1);
    long i, l = lg(Dm);
    for (i = 2; i < l; i++)
    { /* (e,Nt) = 1; skip i = 1: e = 1, done above */
      long e = Dm[i], me = m / e, a = mfcharevalord(CHI, e, ord);
      GEN c, C = powuu(e, r - 1);
      if (kross(st, e) == -1) C = negi(C);
      c = Qab_Czeta(a, ord, C, vt);
      S = gadd(S, gmul(c, gel(A, me*me + 1)));
    }
    gel(R, m+1) = S;
  }
  return degpol(Pn) > 1? gmodulo(R, Pn): R;
}

static long
mfisinkohnen(GEN mf, GEN F)
{
  GEN v, gk = MF_get_gk(mf), CHI = MF_get_CHI(mf);
  long i, eps, N4 = MF_get_N(mf) >> 2, sb = mfsturmNgk(N4 << 4, gk) + 1;
  eps = N4 % mfcharconductor(CHI)? -1 : 1;
  if (odd(MF_get_r(mf))) eps = -eps;
  v = mfcoefs(F, sb, 1);
  for (i = 2;     i <= sb; i+=4) if (!gequal0(gel(v,i+1))) return 0;
  for (i = 2+eps; i <= sb; i+=4) if (!gequal0(gel(v,i+1))) return 0;
  return 1;
}

static long
mfshimura_space_cusp(GEN mf)
{
  long N4;
  if (MF_get_r(mf) == 1 && (N4 = MF_get_N(mf) >> 2) >= 4)
  {
    GEN E = gel(myfactoru(N4), 2);
    long ma = vecsmall_max(E);
    if (ma > 2 || (ma == 2 && !mfcharistrivial(MF_get_CHI(mf)))) return 0;
  }
  return 1;
}

/* D is either a discriminant (not necessarily fundamental) with
   sign(D)=(-1)^{k-1/2}*eps, or a positive squarefree integer t, which is then
   transformed into a fundamental discriminant of the correct sign. */
GEN
mfshimura(GEN mf, GEN F, long t)
{
  pari_sp av = avma;
  GEN G, res, mf2, CHI;
  long sb, M, r, N, space = mf_FULL;

  if (!checkmf_i(F)) pari_err_TYPE("mfshimura",F);
  mf = checkMF(mf);
  r = MF_get_r(mf);
  if (r <= 0) pari_err_DOMAIN("mfshimura", "weight", "<=", ghalf, mf_get_gk(F));
  if (t <= 0 || !uissquarefree(t)) pari_err_TYPE("mfshimura [t]", stoi(t));
  N = MF_get_N(mf); M = N >> 1;
  if (mfiscuspidal(mf,F))
  {
    if (mfshimura_space_cusp(mf)) space = mf_CUSP;
    if (mfisinkohnen(mf,F)) M = N >> 2;
  }
  CHI = MF_get_CHI(mf);
  mf2 = mfinit_Nkchi(M, r << 1, mfcharpow(CHI, gen_2), space, 0);
  sb = mfsturm(mf2);
  G = RgV_shimura(mfcoefs_i(F, sb*sb, t), sb, t, N, r, CHI);
  res = mftobasis_i(mf2, G);
  /* not mflinear(mf2,): we want lowest possible level */
  G = mflinear(MF_get_basis(mf2), res);
  return gerepilecopy(av, mkvec3(mf2, G, res));
}

/* W ZabM (ZM if n = 1), a t_INT or NULL, b t_INT, ZXQ mod P or NULL.
 * Write a/b = A/d with d t_INT and A Zab return [W,d,A,P] */
static GEN
mkMinv(GEN W, GEN a, GEN b, GEN P)
{
  GEN A = (b && typ(b) == t_POL)? Q_remove_denom(QXQ_inv(b,P), &b): NULL;
  if (a && b)
  {
    a = Qdivii(a,b);
    if (typ(a) == t_INT) b = gen_1; else { b = gel(a,2); a = gel(a,1); }
    if (is_pm1(a)) a = NULL;
  }
  if (a) A = A? ZX_Z_mul(A,a): a; else if (!A) A = gen_1;
  if (!b) b = gen_1;
  if (!P) P = gen_0;
  return mkvec4(W,b,A,P);
}
/* M square invertible QabM, return [M',d], M*M' = d*Id */
static GEN
QabM_Minv(GEN M, GEN P, long n)
{
  GEN dW, W, dM;
  M = Q_remove_denom(M, &dM);
  W = P? ZabM_inv(liftpol_shallow(M), P, n, &dW): ZM_inv(M, &dW);
  return mkMinv(W, dM, dW, P);
}
/* Simplified form of mfclean, after a QabM_indexrank: M a ZabM with full
 * column rank and z = indexrank(M) is known */
static GEN
mfclean2(GEN M, GEN z, GEN P, long n)
{
  GEN d, Minv, y = gel(z,1), W = rowpermute(M, y);
  W = P? ZabM_inv(liftpol_shallow(W), P, n, &d): ZM_inv(W, &d);
  M = rowslice(M, 1, y[lg(y)-1]);
  Minv = mkMinv(W, NULL, d, P);
  return mkvec3(y, Minv, M);
}
/* M QabM, lg(M)>1 and [y,z] its rank profile. Let Minv be the inverse of the
 * invertible square matrix in mkMinv format. Return [y,Minv, M[..y[#y],]]
 * P cyclotomic polynomial of order n > 2 or NULL */
static GEN
mfclean(GEN M, GEN P, long n, int ratlift)
{
  GEN W, v, y, z, d, Minv, dM, MdM = Q_remove_denom(M, &dM);
  if (n <= 2)
    W = ZM_pseudoinv(MdM, &v, &d);
  else
    W = ZabM_pseudoinv_i(liftpol_shallow(MdM), P, n, &v, &d, ratlift);
  y = gel(v,1);
  z = gel(v,2);
  if (lg(z) != lg(MdM)) M = vecpermute(M,z);
  M = rowslice(M, 1, y[lg(y)-1]);
  Minv = mkMinv(W, dM, d, P);
  return mkvec3(y, Minv, M);
}
/* call mfclean using only CHI */
static GEN
mfcleanCHI(GEN M, GEN CHI, int ratlift)
{
  long n = mfcharorder(CHI);
  GEN P = (n <= 2)? NULL: mfcharpol(CHI);
  return mfclean(M, P, n, ratlift);
}

/* DATA component of a t_MF_NEWTRACE. Was it stripped to save memory ? */
static int
newtrace_stripped(GEN DATA)
{ return DATA && (lg(DATA) == 5 && typ(gel(DATA,3)) == t_INT); }
/* f a t_MF_NEWTRACE */
static GEN
newtrace_DATA(long N, GEN f)
{
  GEN DATA = gel(f,2);
  return newtrace_stripped(DATA)? initnewtrace(N, DATA): DATA;
}
/* reset cachenew for new level incorporating new DATA, tf a t_MF_NEWTRACE
 * (+ possibly initialize 'full' for new allowed levels) */
static void
reset_cachenew(cachenew_t *cache, long N, GEN tf)
{
  long i, n, l;
  GEN v, DATA = newtrace_DATA(N,tf);
  cache->DATA = DATA;
  if (!DATA) return;
  n = cache->n;
  v = cache->vfull; l = N+1; /* = lg(DATA) */
  for (i = 1; i < l; i++)
    if (typ(gel(v,i)) == t_INT && lg(gel(DATA,i)) != 1)
      gel(v,i) = const_vec(n, NULL);
  cache->VCHIP = gel(gel(DATA,N),_VCHIP);
}
/* initialize a cache of newtrace / cusptrace up to index n and level | N;
 * DATA may be NULL (<=> Tr^new = 0). tf a t_MF_NEWTRACE */
static void
init_cachenew(cachenew_t *cache, long n, long N, GEN tf)
{
  long i, l = N+1; /* = lg(tf.DATA) when DATA != NULL */
  GEN v;
  cache->n = n;
  cache->vnew = v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(v,i) = (N % i)? gen_0: const_vec(n, NULL);
  cache->newHIT = cache->newTOTAL = cache->cuspHIT = cache->cuspTOTAL = 0;
  cache->vfull = v = zerovec(N);
  reset_cachenew(cache, N, tf);
}
static void
dbg_cachenew(cachenew_t *C)
{
  if (DEBUGLEVEL >= 2 && C)
    err_printf("newtrace cache hits: new = %ld/%ld, cusp = %ld/%ld\n",
                    C->newHIT, C->newTOTAL, C->cuspHIT, C->cuspTOTAL);
}

/* newtrace_{N,k}(d*i), i = n0, ..., n */
static GEN
colnewtrace(long n0, long n, long d, long N, long k, cachenew_t *cache)
{
  GEN v = cgetg(n-n0+2, t_COL);
  long i;
  for (i = n0; i <= n; i++) gel(v, i-n0+1) = mfnewtracecache(N, k, i*d, cache);
  return v;
}
/* T_n(l*m0, l*(m0+1), ..., l*m) F, F = t_MF_NEWTRACE [N,k],DATA, cache
 * contains DATA != NULL as well as cached values of F */
static GEN
heckenewtrace(long m0, long m, long l, long N, long NBIG, long k, long n, cachenew_t *cache)
{
  long lD, a, k1, nl = n*l;
  GEN D, V, v = colnewtrace(m0, m, nl, N, k, cache); /* d=1 */
  GEN VCHIP;
  if (n == 1) return v;
  VCHIP = cache->VCHIP;
  D = mydivisorsu(u_ppo(n, NBIG)); lD = lg(D);
  k1 = k - 1;
  for (a = 2; a < lD; a++)
  { /* d > 1, (d,NBIG) = 1 */
    long i, j, d = D[a], c = ugcd(l, d), dl = d/c, m0d = ceildivuu(m0, dl);
    GEN C = vchip_lift(VCHIP, d, powuu(d, k1));
    /* m0=0: i = 1 => skip F(0) = 0 */
    if (!m0) { i = 1; j = dl; } else { i = 0; j = m0d*dl; }
    V = colnewtrace(m0d, m/dl, nl/(d*c), N, k, cache);
    /* C = chi(d) d^(k-1) */
    for (; j <= m; i++, j += dl)
      gel(v,j-m0+1) = gadd(gel(v,j-m0+1), vchip_mod(VCHIP, gmul(C,gel(V,i+1))));
  }
  return v;
}

/* Given v = an[i], return an[d*i], i=0..n */
static GEN
anextract(GEN v, long n, long d)
{
  long i, id, l = n + 2;
  GEN w = cgetg(l, t_VEC);
  if (d == 1)
    for (i = 1; i < l; i++) gel(w, i) = gel(v, i);
  else
    for (i = id = 1; i < l; i++, id += d) gel(w, i) = gel(v, id);
  return w;
}
/* T_n(F)(0, l, ..., l*m) */
static GEN
hecke_i(long m, long l, GEN V, GEN F, GEN DATA)
{
  long k, n, nNBIG, NBIG, lD, M, a, t, nl;
  GEN D, v, CHI;
  if (typ(DATA) == t_VEC)
  { /* 1/2-integral k */
    if (!V) { GEN S = gel(DATA,2); V = mfcoefs_i(F, m*l*S[3], S[4]); }
    return RgV_heckef2(m, l, V, F, DATA);
  }
  k = mf_get_k(F);
  n = DATA[1]; nl = n*l;
  nNBIG = DATA[2];
  NBIG = DATA[3];
  if (nNBIG == 1) return V? V: mfcoefs_i(F,m,nl);
  if (!V && mf_get_type(F) == t_MF_NEWTRACE)
  { /* inline F to allow cache, T_n at level NBIG acting on Tr^new(N,k,CHI) */
    cachenew_t cache;
    long N = mf_get_N(F);
    init_cachenew(&cache, m*nl, N, F);
    v = heckenewtrace(0, m, l, N, NBIG, k, n, &cache);
    dbg_cachenew(&cache);
    settyp(v, t_VEC); return v;
  }
  CHI = mf_get_CHI(F);
  D = mydivisorsu(nNBIG); lD = lg(D);
  M = m + 1;
  t = nNBIG * ugcd(nNBIG, l);
  if (!V) V = mfcoefs_i(F, m * t, nl / t); /* usually nl = t */
  v = anextract(V, m, t); /* mfcoefs(F, m, nl); d = 1 */
  for (a = 2; a < lD; a++)
  { /* d > 1, (d, NBIG) = 1 */
    long d = D[a], c = ugcd(l, d), dl = d/c, i, idl;
    GEN C = gmul(mfchareval(CHI, d), powuu(d, k-1));
    GEN w = anextract(V, m/dl, t/(d*c)); /* mfcoefs(F, m/dl, nl/(d*c)) */
    for (i = idl = 1; idl <= M; i++, idl += dl)
      gel(v,idl) = gadd(gel(v,idl), gmul(C, gel(w,i)));
  }
  return v;
}

static GEN
mkmf(GEN x1, GEN x2, GEN x3, GEN x4, GEN x5)
{
  GEN MF = obj_init(5, MF_SPLITN);
  gel(MF,1) = x1;
  gel(MF,2) = x2;
  gel(MF,3) = x3;
  gel(MF,4) = x4;
  gel(MF,5) = x5; return MF;
}

/* return an integer b such that p | b => T_p^k Tr^new = 0, for all k > 0 */
static long
get_badj(long N, long FC)
{
  GEN fa = myfactoru(N), P = gel(fa,1), E = gel(fa,2);
  long i, b = 1, l = lg(P);
  for (i = 1; i < l; i++)
    if (E[i] > 1 && u_lval(FC, P[i]) < E[i]) b *= P[i];
  return b;
}
/* in place, assume perm strictly increasing */
static void
vecpermute_inplace(GEN v, GEN perm)
{
  long i, l = lg(perm);
  for (i = 1; i < l; i++) gel(v,i) = gel(v,perm[i]);
}

/* Find basis of newspace using closures; assume k >= 2 and !badchar.
 * Return NULL if space is empty, else
 * [mf1, list of closures T(j)traceform, list of corresponding j, matrix] */
static GEN
mfnewinit(long N, long k, GEN CHI, cachenew_t *cache, long init)
{
  GEN S, vj, M, CHIP, mf1, listj, P, tf;
  long j, ct, ctlj, dim, jin, SB, sb, two, ord, FC, badj;

  dim = mfnewdim(N, k, CHI);
  if (!dim && !init) return NULL;
  sb = mfsturmNk(N, k);
  CHIP = mfchartoprimitive(CHI, &FC);
  /* remove newtrace data from S to save space in output: negligible slowdown */
  tf = tag(t_MF_NEWTRACE, mkNK(N,k,CHIP), CHIP);
  badj = get_badj(N, FC);
  /* try sbsmall first: Sturm bound not sharp for new space */
  SB = ceilA1(N, k);
  listj = cgetg(2*sb + 3, t_VECSMALL);
  for (j = ctlj = 1; ctlj < 2*sb + 3; j++)
    if (ugcd(j, badj) == 1) listj[ctlj++] = j;
  if (init)
  {
    init_cachenew(cache, (SB+1)*listj[dim+1], N, tf);
    if (init == -1 || !dim) return NULL; /* old space or dim = 0 */
  }
  else
    reset_cachenew(cache, N, tf);
  /* cache.DATA is not NULL */
  ord = mfcharorder(CHIP);
  P = ord <= 2? NULL: mfcharpol(CHIP);
  vj = cgetg(dim+1, t_VECSMALL);
  M = cgetg(dim+1, t_MAT);
  for (two = 1, ct = 0, jin = 1; two <= 2; two++)
  {
    long a, jlim = jin + sb;
    for (a = jin; a <= jlim; a++)
    {
      GEN z, vecz;
      ct++; vj[ct] = listj[a];
      gel(M, ct) = heckenewtrace(0, SB, 1, N, N, k, vj[ct], cache);
      if (ct < dim) continue;

      z = QabM_indexrank(M, P, ord);
      vecz = gel(z, 2); ct = lg(vecz) - 1;
      if (ct == dim) { M = mkvec3(z, gen_0, M); break; } /*maximal rank, done*/
      vecpermute_inplace(M, vecz);
      vecpermute_inplace(vj, vecz);
    }
    if (a <= jlim) break;
    /* sbsmall was not sufficient, use Sturm bound: must extend M */
    for (j = 1; j <= ct; j++)
    {
      GEN t = heckenewtrace(SB + 1, sb, 1, N, N, k, vj[j], cache);
      gel(M,j) = shallowconcat(gel(M, j), t);
    }
    jin = jlim + 1; SB = sb;
  }
  S = cgetg(dim + 1, t_VEC);
  for (j = 1; j <= dim; j++) gel(S, j) = mfhecke_i(vj[j], N, tf);
  dbg_cachenew(cache);
  mf1 = mkvec4(utoipos(N), utoipos(k), CHI, utoi(mf_NEW));
  return mkmf(mf1, cgetg(1,t_VEC), S, vj, M);
}
/* k > 1 integral, mf space is mf_CUSP or mf_FULL */
static GEN
mfinittonew(GEN mf)
{
  GEN CHI = MF_get_CHI(mf), S = MF_get_S(mf), vMjd = MFcusp_get_vMjd(mf);
  GEN M = MF_get_M(mf), vj, mf1;
  long i, j, l, l0 = lg(S), N0 = MF_get_N(mf);
  for (i = l0-1; i > 0; i--)
  {
    long N = gel(vMjd,i)[1];
    if (N != N0) break;
  }
  if (i == l0-1) return NULL;
  S = vecslice(S, i+1, l0-1); /* forms of conductor N0 */
  l = lg(S); vj = cgetg(l, t_VECSMALL);
  for (j = 1; j < l; j++) vj[j] = gel(vMjd,j+i)[2];
  M = vecslice(M, lg(M)-lg(S)+1, lg(M)-1); /* their coefficients */
  M = mfcleanCHI(M, CHI, 0);
  mf1 = mkvec4(utoipos(N0), MF_get_gk(mf), CHI, utoi(mf_NEW));
  return mkmf(mf1, cgetg(1,t_VEC), S, vj, M);
}

/* Bd(f)[m0..m], v = f[ceil(m0/d)..floor(m/d)], m0d = ceil(m0/d) */
static GEN
RgC_Bd_expand(long m0, long m, GEN v, long d, long m0d)
{
  long i, j;
  GEN w;
  if (d == 1) return v;
  w = zerocol(m-m0+1);
  if (!m0) { i = 1; j = d; } else { i = 0; j = m0d*d; }
  for (; j <= m; i++, j += d) gel(w,j-m0+1) = gel(v,i+1);
  return w;
}
/* S a nonempty vector of t_MF_BD(t_MF_HECKE(t_MF_NEWTRACE)); M the matrix
 * of their coefficients r*0, r*1, ..., r*m0 (~ mfvectomat) or NULL (empty),
 * extend it to coeffs up to m > m0. The forms B_d(T_j(tf_N))in S should be
 * sorted by level N, then j, then increasing d. No reordering here. */
static GEN
bhnmat_extend(GEN M, long m, long r, GEN S, cachenew_t *cache)
{
  long i, mr, m0, m0r, Nold = 0, jold = 0, l = lg(S);
  GEN MAT = cgetg(l, t_MAT), v = NULL;
  if (M) { m0 = nbrows(M); m0r = m0 * r; } else m0 = m0r = 0;
  mr = m*r;
  for (i = 1; i < l; i++)
  {
    long d, j, md, N;
    GEN c, f = bhn_parse(gel(S,i), &d,&j); /* t_MF_NEWTRACE */
    N = mf_get_N(f);
    md = ceildivuu(m0r,d);
    if (N != Nold) { reset_cachenew(cache, N, f); Nold = N; jold = 0; }
    if (!cache->DATA) { gel(MAT,i) = zerocol(m+1); continue; }
    if (j != jold || md)
    { v = heckenewtrace(md, mr/d, 1, N, N, mf_get_k(f), j,cache); jold=j; }
    c = RgC_Bd_expand(m0r, mr, v, d, md);
    if (r > 1) c = c_deflate(m-m0, r, c);
    if (M) c = shallowconcat(gel(M,i), c);
    gel(MAT,i) = c;
  }
  return MAT;
}

/* k > 1 */
static GEN
mfinitcusp(long N, long k, GEN CHI, cachenew_t *cache, long space)
{
  long L, l, lDN1, FC, N1, d1, i, init;
  GEN vS, vMjd, DN1, vmf, CHIP = mfchartoprimitive(CHI, &FC);

  d1 = (space == mf_OLD)? mfolddim_i(N, k, CHIP, NULL): mfcuspdim(N, k, CHIP);
  if (!d1) return NULL;
  N1 = N/FC; DN1 = mydivisorsu(N1); lDN1 = lg(DN1);
  init = (space == mf_OLD)? -1: 1;
  vmf = cgetg(lDN1, t_VEC);
  for (i = lDN1 - 1, l = 1; i; i--)
  { /* by decreasing level to allow cache */
    GEN mf = mfnewinit(FC*DN1[i], k, CHIP, cache, init);
    if (mf) gel(vmf, l++) = mf;
    init = 0;
  }
  setlg(vmf,l); vmf = vecreverse(vmf); /* reorder by increasing level */

  L = mfsturmNk(N, k)+1;
  vS = vectrunc_init(L);
  vMjd = vectrunc_init(L);
  for (i = 1; i < l; i++)
  {
    GEN DNM, mf = gel(vmf,i), S = MF_get_S(mf), vj = MFnew_get_vj(mf);
    long a, lDNM, lS = lg(S), M = MF_get_N(mf);
    DNM = mydivisorsu(N / M); lDNM = lg(DNM);
    for (a = 1; a < lS; a++)
    {
      GEN tf = gel(S,a);
      long b, j = vj[a];
      for (b = 1; b < lDNM; b++)
      {
        long d = DNM[b];
        vectrunc_append(vS, mfbd_i(tf, d));
        vectrunc_append(vMjd, mkvecsmall3(M, j, d));
      }
    }
  }
  return mkmf(NULL, cgetg(1, t_VEC), vS, vMjd, NULL);
}

long
mfsturm_mf(GEN mf)
{
  GEN Mindex = MF_get_Mindex(mf);
  long n = lg(Mindex)-1;
  return n? Mindex[n]-1: 0;
}

long
mfsturm(GEN T)
{
  long N, nk, dk;
  GEN CHI, mf = checkMF_i(T);
  if (mf) return mfsturm_mf(mf);
  checkNK2(T, &N, &nk, &dk, &CHI, 0);
  return dk == 1 ? mfsturmNk(N, nk) : mfsturmNk(N, (nk + 1) >> 1);
}
long
mfisequal(GEN F, GEN G, long lim)
{
  pari_sp av = avma;
  long b;
  if (!checkmf_i(F)) pari_err_TYPE("mfisequal",F);
  if (!checkmf_i(G)) pari_err_TYPE("mfisequal",G);
  b = lim? lim: maxss(mfsturmmf(F), mfsturmmf(G));
  return gc_long(av, gequal(mfcoefs_i(F, b, 1), mfcoefs_i(G, b, 1)));
}

GEN
mffields(GEN mf)
{
  if (checkmf_i(mf)) return gcopy(mf_get_field(mf));
  mf = checkMF(mf); return gcopy(MF_get_fields(mf));
}

GEN
mfeigenbasis(GEN mf)
{
  pari_sp ltop = avma;
  GEN F, S, v, vP;
  long i, l, k, dS;

  mf = checkMF(mf);
  k = MF_get_k(mf);
  S = MF_get_S(mf); dS = lg(S)-1;
  if (!dS) return cgetg(1, t_VEC);
  F = MF_get_newforms(mf);
  vP = MF_get_fields(mf);
  if (k == 1)
  {
    if (MF_get_space(mf) == mf_FULL)
    {
      long dE = lg(MF_get_E(mf)) - 1;
      if (dE) F = rowslice(F, dE+1, dE+dS);
    }
    v = vecmflineardiv_linear(S, F);
    l = lg(v);
  }
  else
  {
    GEN (*L)(GEN, GEN) = (MF_get_space(mf) == mf_FULL)? mflinear: mflinear_bhn;
    l = lg(F); v = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(v,i) = L(mf, gel(F,i));
  }
  for (i = 1; i < l; i++) mf_setfield(gel(v,i), gel(vP,i));
  return gerepilecopy(ltop, v);
}

/* Minv = [M, d, A], v a t_COL; A a Zab, d a t_INT; return (A/d) * M*v */
static GEN
Minv_RgC_mul(GEN Minv, GEN v)
{
  GEN M = gel(Minv,1), d = gel(Minv,2), A = gel(Minv,3);
  v = RgM_RgC_mul(M, v);
  if (!equali1(A))
  {
    if (typ(A) == t_POL && degpol(A) > 0) A = mkpolmod(A, gel(Minv,4));
    v = RgC_Rg_mul(v, A);
  }
  if (!equali1(d)) v = RgC_Rg_div(v, d);
  return v;
}
static GEN
Minv_RgM_mul(GEN Minv, GEN B)
{
  long j, l = lg(B);
  GEN M = cgetg(l, t_MAT);
  for (j = 1; j < l; j++) gel(M,j) = Minv_RgC_mul(Minv, gel(B,j));
  return M;
}
/* B * Minv; allow B = NULL for Id */
static GEN
RgM_Minv_mul(GEN B, GEN Minv)
{
  GEN M = gel(Minv,1), d = gel(Minv,2), A = gel(Minv,3);
  if (B) M = RgM_mul(B, M);
  if (!equali1(A))
  {
    if (typ(A) == t_POL) A = mkpolmod(A, gel(Minv,4));
    M = RgM_Rg_mul(M, A);
  }
  if (!equali1(d)) M = RgM_Rg_div(M,d);
  return M;
}

/* perm vector of strictly increasing indices, v a vector or arbitrary length;
 * the last r entries of perm fall beyond v.
 * Return v o perm[1..(-r)], discarding the last r entries of v */
static GEN
vecpermute_partial(GEN v, GEN perm, long *r)
{
  long i, n = lg(v)-1, l = lg(perm);
  GEN w;
  if (perm[l-1] <= n) { *r = 0; return vecpermute(v,perm); }
  for (i = 1; i < l; i++)
    if (perm[i] > n) break;
  *r = l - i; l = i;
  w = cgetg(l, typ(v));
  for (i = 1; i < l; i++) gel(w,i) = gel(v,perm[i]);
  return w;
}

/* given form F, find coeffs of F on mfbasis(mf). If power series, not
 * guaranteed correct if precision less than Sturm bound */
static GEN
mftobasis_i(GEN mf, GEN F)
{
  GEN v, Mindex, Minv;
  if (!MF_get_dim(mf)) return cgetg(1, t_COL);
  Mindex = MF_get_Mindex(mf);
  Minv = MF_get_Minv(mf);
  if (checkmf_i(F))
  {
    long n = Mindex[lg(Mindex)-1];
    v = vecpermute(mfcoefs_i(F, n, 1), Mindex);
    return Minv_RgC_mul(Minv, v);
  }
  else
  {
    GEN A = gel(Minv,1), d = gel(Minv,2);
    long r;
    v = F;
    switch(typ(F))
    {
      case t_SER: v = sertocol(v);
      case t_VEC: case t_COL: break;
      default: pari_err_TYPE("mftobasis", F);
    }
    if (lg(v) == 1) pari_err_TYPE("mftobasis",v);
    v = vecpermute_partial(v, Mindex, &r);
    if (!r) return Minv_RgC_mul(Minv, v); /* single solution */
    /* affine space of dimension r */
    v = RgM_RgC_mul(vecslice(A, 1, lg(v)-1), v);
    if (!equali1(d)) v = RgC_Rg_div(v,d);
    return mkvec2(v, vecslice(A, lg(A)-r, lg(A)-1));
  }
}

static GEN
const_mat(long n, GEN x)
{
  long j, l = n+1;
  GEN A = cgetg(l,t_MAT);
  for (j = 1; j < l; j++) gel(A,j) = const_col(n, x);
  return A;
}

/* L is the mftobasis of a form on CUSP space. We allow mf_FULL or mf_CUSP */
static GEN
mftonew_i(GEN mf, GEN L, long *plevel)
{
  GEN S, listMjd, CHI, res, Aclos, Acoef, D, perm;
  long N1, LC, lD, i, l, t, level, N = MF_get_N(mf);

  if (MF_get_k(mf) == 1) pari_err_IMPL("mftonew in weight 1");
  listMjd = MFcusp_get_vMjd(mf);
  CHI = MF_get_CHI(mf); LC = mfcharconductor(CHI);
  S = MF_get_S(mf);

  N1 = N/LC;
  D = mydivisorsu(N1); lD = lg(D);
  perm = cgetg(N1+1, t_VECSMALL);
  for (i = 1; i < lD; i++) perm[D[i]] = i;
  Aclos = const_mat(lD-1, cgetg(1,t_VEC));
  Acoef = const_mat(lD-1, cgetg(1,t_VEC));
  l = lg(listMjd);
  for (i = 1; i < l; i++)
  {
    long M, d;
    GEN v;
    if (gequal0(gel(L,i))) continue;
    v = gel(listMjd, i);
    M = perm[ v[1]/LC ];
    d = perm[ v[3] ];
    gcoeff(Aclos,M,d) = vec_append(gcoeff(Aclos,M,d), gel(S,i));
    gcoeff(Acoef,M,d) = shallowconcat(gcoeff(Acoef,M,d), gel(L,i));
  }
  res = cgetg(l, t_VEC); level = 1;
  for (i = t = 1; i < lD; i++)
  {
    long j, M = D[i]*LC;
    GEN gM = utoipos(M);
    for (j = 1; j < lD; j++)
    {
      GEN f = gcoeff(Aclos,i,j), C, NK;
      long d;
      if (lg(f) == 1) continue;
      NK = mf_get_NK(gel(f,1));
      d = D[j];
      C = gcoeff(Acoef,i,j);
      level = ulcm(level, M*d);
      gel(res,t++) = mkvec3(gM, utoipos(d), mflinear_i(NK,f,C));
    }
  }
  if (plevel) *plevel = level;
  setlg(res, t); return res;
}
GEN
mftonew(GEN mf, GEN F)
{
  pari_sp av = avma;
  GEN ES;
  long s;
  mf = checkMF(mf);
  s = MF_get_space(mf);
  if (s != mf_FULL && s != mf_CUSP)
    pari_err_TYPE("mftonew [not a full or cuspidal space]", mf);
  ES = mftobasisES(mf,F);
  if (!gequal0(gel(ES,1)))
    pari_err_TYPE("mftonew [not a cuspidal form]", F);
  F = gel(ES,2);
  return gerepilecopy(av, mftonew_i(mf,F, NULL));
}

static GEN mfeisenstein_i(long k, GEN CHI1, GEN CHI2);

/* mfinit(F * Theta) */
static GEN
mf2init(GEN mf)
{
  GEN CHI = MF_get_CHI(mf), gk = gadd(MF_get_gk(mf), ghalf);
  long N = MF_get_N(mf);
  return mfinit_Nkchi(N, itou(gk), mfchiadjust(CHI, gk, N), mf_FULL, 0);
}

static long
mfvec_first_cusp(GEN v)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN F = gel(v,i);
    long t = mf_get_type(F);
    if (t == t_MF_BD) { F = gel(F,2); t = mf_get_type(F); }
    if (t == t_MF_HECKE) { F = gel(F,3); t = mf_get_type(F); }
    if (t == t_MF_NEWTRACE) break;
  }
  return i;
}
/* vF a vector of mf F of type DIV(LINEAR(BAS,L), f) in (lcm) level N,
 * F[2]=LINEAR(BAS,L), F[2][2]=BAS=fixed basis (Eisenstein or bhn type),
 * F[2][3]=L, F[3]=f; mfvectomat(vF, n) */
static GEN
mflineardivtomat(long N, GEN vF, long n)
{
  GEN F, M, f, fc, ME, dB, B, a0, V = NULL;
  long lM, lF = lg(vF), j;

  if (lF == 1) return cgetg(1,t_MAT);
  F = gel(vF,1);
  if (lg(F) == 5)
  { /* chicompat */
    V = gmael(F,4,4);
    if (typ(V) == t_INT) V = NULL;
  }
  M = gmael(F,2,2); /* BAS */
  lM = lg(M);
  j = mfvec_first_cusp(M);
  if (j == 1) ME = NULL;
  else
  { /* BAS starts by Eisenstein */
    ME = mfvectomat(vecslice(M,1,j-1), n, 1);
    M = vecslice(M, j,lM-1);
  }
  M = bhnmat_extend_nocache(NULL, N, n, 1, M);
  if (ME) M = shallowconcat(ME,M);
  /* M = mfcoefs of BAS */
  B = cgetg(lF, t_MAT);
  dB= cgetg(lF, t_VEC);
  for (j = 1; j < lF; j++)
  {
    GEN g = gel(vF, j); /* t_MF_DIV */
    gel(B,j) = RgM_RgC_mul(M, gmael(g,2,3));
    gel(dB,j)= gmael(g,2,4);
  }
  f = mfcoefsser(gel(F,3),n);
  a0 = polcoef_i(f, 0, -1);
  if (gequal0(a0) || gequal1(a0))
    a0 = NULL;
  else
    f = gdiv(ser_unscale(f, a0), a0);
  fc = ginv(f);
  for (j = 1; j < lF; j++)
  {
    pari_sp av = avma;
    GEN LISer = RgV_to_ser_full(gel(B,j)), f;
    if (a0) LISer = gdiv(ser_unscale(LISer, a0), a0);
    f = gmul(LISer, fc);
    if (a0) f = ser_unscale(f, ginv(a0));
    f = sertocol(f); setlg(f, n+2);
    if (!gequal1(gel(dB,j))) f = RgC_Rg_div(f, gel(dB,j));
    gel(B,j) = gerepileupto(av,f);
  }
  if (V) B = gmodulo(QabM_tracerel(V, 0, B), gel(V,1));
  return B;
}

static GEN
mfheckemat_mfcoefs(GEN mf, GEN B, GEN DATA)
{
  GEN Mindex = MF_get_Mindex(mf), Minv = MF_get_Minv(mf);
  long j, l = lg(B), sb = mfsturm_mf(mf);
  GEN b = MF_get_basis(mf), Q = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
  {
    GEN v = hecke_i(sb, 1, gel(B,j), gel(b,j), DATA); /* Tn b[j] */
    settyp(v,t_COL); gel(Q,j) = vecpermute(v, Mindex);
  }
  return Minv_RgM_mul(Minv,Q);
}
/* T_p^2, p prime, 1/2-integral weight; B = mfcoefs(mf,sb*p^2,1) or (mf,sb,p^2)
 * if p|N */
static GEN
mfheckemat_mfcoefs_p2(GEN mf, long p, GEN B)
{
  pari_sp av = avma;
  GEN DATA = heckef2_data(MF_get_N(mf), p*p);
  return gerepileupto(av, mfheckemat_mfcoefs(mf, B, DATA));
}
/* convert Mindex from row-index to mfcoef indexation: a(n) is stored in
 * mfcoefs()[n+1], so subtract 1 from all indices */
static GEN
Mindex_as_coef(GEN mf)
{
  GEN v, Mindex = MF_get_Mindex(mf);
  long i, l = lg(Mindex);
  v = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) v[i] = Mindex[i]-1;
  return v;
}
/* T_p, p prime; B = mfcoefs(mf,sb*p,1) or (mf,sb,p) if p|N; integral weight */
static GEN
mfheckemat_mfcoefs_p(GEN mf, long p, GEN B)
{
  pari_sp av = avma;
  GEN vm, Q, C, Minv = MF_get_Minv(mf);
  long lm, k, i, j, l = lg(B), N = MF_get_N(mf);

  if (N % p == 0) return Minv_RgM_mul(Minv, rowpermute(B, MF_get_Mindex(mf)));
  k = MF_get_k(mf);
  C = gmul(mfchareval(MF_get_CHI(mf), p), powuu(p, k-1));
  vm = Mindex_as_coef(mf); lm = lg(vm);
  Q = cgetg(l, t_MAT);
  for (j = 1; j < l; j++) gel(Q,j) = cgetg(lm, t_COL);
  for (i = 1; i < lm; i++)
  {
    long m = vm[i], mp = m*p;
    GEN Cm = (m % p) == 0? C : NULL;
    for (j = 1; j < l; j++)
    {
      GEN S = gel(B,j), s = gel(S, mp + 1);
      if (Cm) s = gadd(s, gmul(C, gel(S, m/p + 1)));
      gcoeff(Q, i, j) = s;
    }
  }
  return gerepileupto(av, Minv_RgM_mul(Minv,Q));
}
/* Matrix of T(p), p prime, dim(mf) > 0 and integral weight */
static GEN
mfheckemat_p(GEN mf, long p)
{
  pari_sp av = avma;
  long N = MF_get_N(mf), sb = mfsturm_mf(mf);
  GEN B = (N % p)? mfcoefs_mf(mf, sb * p, 1): mfcoefs_mf(mf, sb, p);
  return gerepileupto(av, mfheckemat_mfcoefs(mf, B, hecke_data(N,p)));
}

/* mf_NEW != (0), weight > 1, p prime. Use
 * T(p) T(j) = T(j*p) + p^{k-1} \chi(p) 1_{p | j, p \nmid N} T(j/p) */
static GEN
mfnewmathecke_p(GEN mf, long p)
{
  pari_sp av = avma;
  GEN tf, vj = MFnew_get_vj(mf), CHI = MF_get_CHI(mf);
  GEN Mindex = MF_get_Mindex(mf), Minv = MF_get_Minv(mf);
  long N = MF_get_N(mf), k = MF_get_k(mf);
  long i, j, lvj = lg(vj), lim = vj[lvj-1] * p;
  GEN M, perm, V, need = zero_zv(lim);
  GEN C = (N % p)? gmul(mfchareval(CHI,p), powuu(p,k-1)): NULL;
  tf = mftraceform_new(N, k, CHI);
  for (i = 1; i < lvj; i++)
  {
    j = vj[i]; need[j*p] = 1;
    if (N % p && j % p == 0) need[j/p] = 1;
  }
  perm = zero_zv(lim);
  V = cgetg(lim+1, t_VEC);
  for (i = j = 1; i <= lim; i++)
    if (need[i]) { gel(V,j) = mfhecke_i(i, N, tf); perm[i] = j; j++; }
  setlg(V, j);
  V = bhnmat_extend_nocache(NULL, N, mfsturm_mf(mf), 1, V);
  V = rowpermute(V, Mindex); /* V[perm[i]] = coeffs(T_i newtrace) */
  M = cgetg(lvj, t_MAT);
  for (i = 1; i < lvj; i++)
  {
    GEN t;
    j = vj[i]; t = gel(V, perm[j*p]);
    if (C && j % p == 0) t = RgC_add(t, RgC_Rg_mul(gel(V, perm[j/p]),C));
    gel(M,i) = t;
  }
  return gerepileupto(av, Minv_RgM_mul(Minv, M));
}

GEN
mfheckemat(GEN mf, GEN vn)
{
  pari_sp av = avma;
  long lv, lvP, i, N, dim, nk, dk, p, sb, flint = (typ(vn)==t_INT);
  GEN CHI, res, vT, FA, B, vP;

  mf = checkMF(mf);
  if (typ(vn) != t_VECSMALL) vn = gtovecsmall(vn);
  N = MF_get_N(mf); CHI = MF_get_CHI(mf); Qtoss(MF_get_gk(mf), &nk, &dk);
  dim = MF_get_dim(mf);
  lv = lg(vn);
  res = cgetg(lv, t_VEC);
  FA = cgetg(lv, t_VEC);
  vP = cgetg(lv, t_VEC);
  vT = const_vec(vecsmall_max(vn), NULL);
  for (i = 1; i < lv; i++)
  {
    ulong n = (ulong)labs(vn[i]);
    GEN fa;
    if (!n) pari_err_TYPE("mfheckemat", vn);
    if (dk == 1 || uissquareall(n, &n)) fa = myfactoru(n);
    else { n = 0; fa = myfactoru(1); } /* dummy: T_{vn[i]} = 0 */
    vn[i] = n;
    gel(FA,i) = fa;
    gel(vP,i) = gel(fa,1);
  }
  vP = shallowconcat1(vP); vecsmall_sort(vP);
  vP = vecsmall_uniq_sorted(vP); /* all primes occurring in vn */
  lvP = lg(vP); if (lvP == 1) goto END;
  p = vP[lvP-1];
  sb = mfsturm_mf(mf);
  if (dk == 1 && nk != 1 && MF_get_space(mf) == mf_NEW)
    B = NULL; /* special purpose mfnewmathecke_p is faster */
  else if (lvP == 2 && N % p == 0)
    B = mfcoefs_mf(mf, sb, dk==2? p*p: p); /* single prime | N, can optimize */
  else
    B = mfcoefs_mf(mf, sb * (dk==2? p*p: p), 1); /* general initialization */
  for (i = 1; i < lvP; i++)
  {
    long j, l, q, e = 1;
    GEN C, Tp, u1, u0;
    p = vP[i];
    for (j = 1; j < lv; j++) e = maxss(e, z_lval(vn[j], p));
    if (!B)
      Tp = mfnewmathecke_p(mf, p);
    else if (dk == 2)
      Tp = mfheckemat_mfcoefs_p2(mf,p, (lvP==2||N%p)? B: matdeflate(sb,p*p,B));
    else
      Tp = mfheckemat_mfcoefs_p(mf, p, (lvP==2||N%p)? B: matdeflate(sb,p,B));
    gel(vT, p) = Tp;
    if (e == 1) continue;
    u0 = gen_1;
    if (dk == 2)
    {
      C = N % p? gmul(mfchareval(CHI,p*p), powuu(p, nk-2)): NULL;
      if (e == 2) u0 = uutoQ(p+1,p); /* special case T_{p^4} */
    }
    else
      C = N % p? gmul(mfchareval(CHI,p),   powuu(p, nk-1)): NULL;
    for (u1=Tp, q=p, l=2; l <= e; l++)
    { /* u0 = T_{p^{l-2}}, u1 = T_{p^{l-1}} for l > 2 */
      GEN v = gmul(Tp, u1);
      if (C) v = gsub(v, gmul(C, u0));
      /* q = p^l, vT[q] = T_q for k integer else T_{q^2} */
      q *= p; u0 = u1; gel(vT, q) = u1 = v;
    }
  }
END:
  /* vT[p^e] = T_{p^e} for all p^e occurring below */
  for (i = 1; i < lv; i++)
  {
    long n = vn[i], j, lP;
    GEN fa, P, E, M;
    if (n == 0) { gel(res,i) = zeromat(dim,dim); continue; }
    if (n == 1) { gel(res,i) = matid(dim); continue; }
    fa = gel(FA,i);
    P = gel(fa,1); lP = lg(P);
    E = gel(fa,2); M = gel(vT, upowuu(P[1], E[1]));
    for (j = 2; j < lP; j++) M = RgM_mul(M, gel(vT, upowuu(P[j], E[j])));
    gel(res,i) = M;
  }
  if (flint) res = gel(res,1);
  return gerepilecopy(av, res);
}

/* f = \sum_i v[i] T_listj[i] (Trace Form) attached to v; replace by f/a_1(f) */
static GEN
mf_normalize(GEN mf, GEN v)
{
  GEN c, dc = NULL, M = MF_get_M(mf), Mindex = MF_get_Mindex(mf);
  v = Q_primpart(v);
  c = RgMrow_RgC_mul(M, v, 2); /* a_1(f) */
  if (gequal1(c)) return v;
  if (typ(c) == t_POL) c = gmodulo(c, mfcharpol(MF_get_CHI(mf)));
  if (typ(c) == t_POLMOD && varn(gel(c,1)) == 1 && degpol(gel(c,1)) >= 40
                         && Mindex[1] == 2
                         && mfcharorder(MF_get_CHI(mf)) <= 2)
  { /* normalize using expansion at infinity (small coefficients) */
    GEN w, P = gel(c,1), a1 = gel(c,2);
    long i, l = lg(Mindex);
    w = cgetg(l, t_COL);
    gel(w,1) = gen_1;
    for (i = 2; i < l; i++)
    {
      c = liftpol_shallow(RgMrow_RgC_mul(M, v, Mindex[i]));
      gel(w,i) = QXQ_div(c, a1, P);
    }
    /* w = expansion at oo of normalized form */
    v = Minv_RgC_mul(MF_get_Minv(mf), Q_remove_denom(w, &dc));
    v = gmodulo(v, P); /* back to mfbasis coefficients */
  }
  else
  {
    c = ginv(c);
    if (typ(c) == t_POLMOD) c = Q_remove_denom(c, &dc);
    v = RgC_Rg_mul(v, c);
  }
  if (dc) v = RgC_Rg_div(v, dc);
  return v;
}
static void
pol_red(GEN NF, GEN *pP, GEN *pa, long flag)
{
  GEN dP, a, P = *pP;
  long d = degpol(P);

  *pa = a = pol_x(varn(P));
  if (d * (NF ? nf_get_degree(NF): 1) > 30) return;

  dP = RgX_disc(P);
  if (typ(dP) != t_INT)
  { dP = gnorm(dP); if (typ(dP) != t_INT) pari_err_BUG("mfnewsplit"); }
  if (d == 2 || expi(dP) < 62)
  {
    if (expi(dP) < 31)
      P = NF? rnfpolredabs(NF, P,flag): polredabs0(P,flag);
    else
      P = NF? rnfpolredbest(NF,P,flag): polredbest(P,flag);
    if (flag)
    {
      a = gel(P,2); if (typ(a) == t_POLMOD) a = gel(a,2);
      P = gel(P,1);
    }
  }
  *pP = P;
  *pa = a;
}

/* Diagonalize and normalize. See mfsplit for meaning of flag. */
static GEN
mfspclean(GEN mf, GEN mf0, GEN NF, long ord, GEN simplesp, long flag)
{
  const long vz = 1;
  long i, l = lg(simplesp), dim = MF_get_dim(mf);
  GEN res = cgetg(l, t_MAT), pols = cgetg(l, t_VEC);
  GEN zeros = (mf == mf0)? NULL: zerocol(dim - MF_get_dim(mf0));
  for (i = 1; i < l; i++)
  {
    GEN ATP = gel(simplesp, i), A = gel(ATP,1), P = gel(ATP,3);
    long d = degpol(P);
    GEN a, v = (flag && d > flag)? NULL: gel(A,1);
    if (d == 1) P = pol_x(vz);
    else
    {
      pol_red(NF, &P, &a, !!v);
      if (v)
      { /* Mod(a,P) root of charpoly(T), K*gpowers(a) = eigenvector of T */
        GEN K, den, M = cgetg(d+1, t_MAT), T = gel(ATP,2);
        long j;
        T = shallowtrans(T);
        gel(M,1) = vec_ei(d,1); /* basis of cyclic vectors */
        for (j = 2; j <= d; j++) gel(M,j) = RgM_RgC_mul(T, gel(M,j-1));
        M = Q_primpart(M);
        K = NF? ZabM_inv(liftpol_shallow(M), nf_get_pol(NF), ord, &den)
              : ZM_inv(M,&den);
        K = shallowtrans(K);
        v = gequalX(a)? pol_x_powers(d, vz): RgXQ_powers(a, d-1, P);
        v = gmodulo(RgM_RgC_mul(A, RgM_RgC_mul(K,v)), P);
      }
    }
    if (v)
    {
      v = mf_normalize(mf0, v); if (zeros) v = shallowconcat(zeros,v);
      gel(res,i) = v; if (flag) setlg(res,i+1);
    }
    else
      gel(res,i) = zerocol(dim);
    gel(pols,i) = P;
  }
  return mkvec2(res, pols);
}

/* return v = v_{X-r}(P), and set Z = P / (X-r)^v */
static long
RgX_valrem_root(GEN P, GEN r, GEN *Z)
{
  long v;
  for (v = 0; degpol(P); v++)
  {
    GEN t, Q = RgX_div_by_X_x(P, r, &t);
    if (!gequal0(t)) break;
    P = Q;
  }
  *Z = P; return v;
}
static GEN
mynffactor(GEN NF, GEN P, long dimlim)
{
  long i, l, v;
  GEN R, E;
  if (dimlim != 1)
  {
    R = NF? nffactor(NF, P): QX_factor(P);
    if (!dimlim) return R;
    E = gel(R,2);
    R = gel(R,1); l = lg(R);
    for (i = 1; i < l; i++)
      if (degpol(gel(R,i)) > dimlim) break;
    if (i == 1) return NULL;
    setlg(E,i);
    setlg(R,i); return mkmat2(R, E);
  }
  /* dimlim = 1 */
  R = nfroots(NF, P); l = lg(R);
  if (l == 1) return NULL;
  v = varn(P);
  settyp(R, t_COL);
  if (degpol(P) == l-1)
    E = const_col(l-1, gen_1);
  else
  {
    E = cgetg(l, t_COL);
    for (i = 1; i < l; i++) gel(E,i) = utoi(RgX_valrem_root(P, gel(R,i), &P));
  }
  R = deg1_from_roots(R, v);
  return mkmat2(R, E);
}

/* Let K be a number field attached to NF (Q if NF = NULL). A K-vector
 * space of dimension d > 0 is given by a t_MAT A (n x d, full column rank)
 * giving a K-basis, X a section (d x n: left pseudo-inverse of A). Return a
 * pair (T, fa), where T is an element of the Hecke algebra (a sum of Tp taken
 * from vector vTp) acting on A (a d x d t_MAT) and fa is the factorization of
 * its characteristic polynomial, limited to factors of degree <= dimlim if
 * dimlim != 0 (return NULL if there are no factors of degree <= dimlim) */
static GEN
findbestsplit(GEN NF, GEN vTp, GEN A, GEN X, long dimlim, long vz)
{
  GEN T = NULL, Tkeep = NULL, fakeep = NULL;
  long lmax = 0, i, lT = lg(vTp);
  for (i = 1; i < lT; i++)
  {
    GEN D, P, E, fa, TpA = gel(vTp,i);
    long l;
    if (typ(TpA) == t_INT) break;
    if (lg(TpA) > lg(A)) TpA = RgM_mul(X, RgM_mul(TpA, A)); /* Tp | A */
    T = T ? RgM_add(T, TpA) : TpA;
    if (!NF) { P = QM_charpoly_ZX(T); setvarn(P, vz); }
    else
    {
      P = charpoly(Q_remove_denom(T, &D), vz);
      if (D) P = gdiv(RgX_unscale(P, D), powiu(D, degpol(P)));
    }
    fa = mynffactor(NF, P, dimlim);
    if (!fa) return NULL;
    E = gel(fa, 2);
    /* characteristic polynomial is separable ? */
    if (isint1(vecmax(E))) { Tkeep = T; fakeep = fa; break; }
    l = lg(E);
    /* characteristic polynomial has more factors than before ? */
    if (l > lmax) { lmax = l; Tkeep = T; fakeep = fa; }
  }
  return mkvec2(Tkeep, fakeep);
}

static GEN
nfcontent(GEN nf, GEN v)
{
  long i, l = lg(v);
  GEN c = gel(v,1);
  for (i = 2; i < l; i++) c = idealadd(nf, c, gel(v,i));
  if (typ(c) == t_MAT && gequal1(gcoeff(c,1,1))) c = gen_1;
  return c;
}
static GEN
nf_primpart(GEN nf, GEN B)
{
  switch(typ(B))
  {
    case t_COL:
    {
      GEN A = matalgtobasis(nf, B), c = nfcontent(nf, A);
      if (typ(c) == t_INT) return B;
      c = idealred_elt(nf,c);
      A = Q_primpart( nfC_nf_mul(nf, A, Q_primpart(nfinv(nf,c))) );
      A = liftpol_shallow( matbasistoalg(nf, A) );
      if (gexpo(A) > gexpo(B)) A = B;
      return A;
    }
    case t_MAT:
    {
      long i, l;
      GEN A = cgetg_copy(B, &l);
      for (i = 1; i < l; i++) gel(A,i) = nf_primpart(nf, gel(B,i));
      return A;
    }
    default:
      pari_err_TYPE("nf_primpart", B);
      return NULL; /*LCOV_EXCL_LINE*/
  }
}

/* rotate entries of v to accomodate new entry 'x' (push out oldest entry) */
static void
vecpush(GEN v, GEN x)
{
  long i;
  for (i = lg(v)-1; i > 1; i--) gel(v,i) = gel(v,i-1);
  gel(v,1) = x;
}

/* sort t_VEC of vector spaces by increasing dimension */
static GEN
sort_by_dim(GEN v)
{
  long i, l = lg(v);
  GEN D = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) D[i] = lg(gmael(v,i,2));
  return vecpermute(v , vecsmall_indexsort(D));
}
static GEN
split_starting_space(GEN mf)
{
  long d = MF_get_dim(mf), d2;
  GEN id = matid(d);
  switch(MF_get_space(mf))
  {
    case mf_NEW:
    case mf_CUSP: return mkvec2(id, id);
  }
  d2 = lg(MF_get_S(mf))-1;
  return mkvec2(vecslice(id, d-d2+1,d),
                shallowconcat(zeromat(d2,d-d2),matid(d2)));
}
/* If dimlim > 0, keep only the dimension <= dimlim eigenspaces.
 * See mfsplit for the meaning of flag. */
static GEN
split_ii(GEN mf, long dimlim, long flag, GEN vSP, long *pnewd)
{
  forprime_t iter;
  GEN CHI = MF_get_CHI(mf), empty = cgetg(1, t_VEC), mf0 = mf;
  GEN NF, POLCYC, todosp, Tpbigvec, simplesp;
  long N = MF_get_N(mf), k = MF_get_k(mf);
  long ord, FC, NEWT, dimsimple = 0, newd = -1;
  const long NBH = 5, vz = 1;
  ulong p;

  switch(MF_get_space(mf))
  {
    case mf_NEW: break;
    case mf_CUSP:
    case mf_FULL:
    {
      GEN CHIP;
      if (k > 1) { mf0 = mfinittonew(mf); break; }
      CHIP = mfchartoprimitive(CHI, NULL);
      newd = lg(MF_get_S(mf))-1 - mfolddim_i(N, k, CHIP, vSP);
      break;
    }
    default: pari_err_TYPE("mfsplit [space does not contain newspace]", mf);
      return NULL; /*LCOV_EXCL_LINE*/
  }
  if (newd < 0) newd = mf0? MF_get_dim(mf0): 0;
  *pnewd = newd;
  if (!newd) return mkvec2(cgetg(1, t_MAT), empty);

  NEWT = (k > 1 && MF_get_space(mf0) == mf_NEW);
  todosp = mkvec( split_starting_space(mf0) );
  simplesp = empty;
  FC = mfcharconductor(CHI);
  ord = mfcharorder(CHI);
  if (ord <= 2) NF = POLCYC = NULL;
  else
  {
    POLCYC = mfcharpol(CHI);
    NF = nfinit(POLCYC,DEFAULTPREC);
  }
  Tpbigvec = zerovec(NBH);
  u_forprime_init(&iter, 2, ULONG_MAX);
  while (dimsimple < newd && (p = u_forprime_next(&iter)))
  {
    GEN nextsp;
    long ind;
    if (N % (p*p) == 0 && N/p % FC == 0) continue; /* T_p = 0 in this case */
    vecpush(Tpbigvec, NEWT? mfnewmathecke_p(mf0,p): mfheckemat_p(mf0,p));
    nextsp = empty;
    for (ind = 1; ind < lg(todosp); ind++)
    {
      GEN tmp = gel(todosp, ind), fa, P, E, D, Tp, DTp;
      GEN A = gel(tmp, 1);
      GEN X = gel(tmp, 2);
      long lP, i;
      tmp = findbestsplit(NF, Tpbigvec, A, X, dimlim, vz);
      if (!tmp) continue; /* nothing there */
      Tp = gel(tmp, 1);
      fa = gel(tmp, 2);
      P = gel(fa, 1);
      E = gel(fa, 2); lP = lg(P);
      /* lP > 1 */
      if (DEBUGLEVEL) err_printf("Exponents = %Ps\n", E);
      if (lP == 2)
      {
        GEN P1 = gel(P,1);
        long e1 = itos(gel(E,1)), d1 = degpol(P1);
        if (e1 * d1 == lg(Tp)-1)
        {
          if (e1 > 1) nextsp = vec_append(nextsp, mkvec2(A,X));
          else
          { /* simple module */
            simplesp = vec_append(simplesp, mkvec3(A,Tp,P1));
            if ((dimsimple += d1) == newd) goto END;
          }
          continue;
        }
      }
      /* Found splitting */
      DTp = Q_remove_denom(Tp, &D);
      for (i = 1; i < lP; i++)
      {
        GEN Ai, Xi, dXi, AAi, v, y, Pi = gel(P,i);
        Ai = RgX_RgM_eval(D? RgX_rescale(Pi,D): Pi, DTp);
        Ai = QabM_ker(Ai, POLCYC, ord);
        if (NF) Ai = nf_primpart(NF, Ai);

        AAi = RgM_mul(A, Ai);
        /* gives section, works on nonsquare matrices */
        Xi = QabM_pseudoinv(Ai, POLCYC, ord, &v, &dXi);
        Xi = RgM_Rg_div(Xi, dXi);
        y = gel(v,1);
        if (isint1(gel(E,i)))
        {
          GEN Tpi = RgM_mul(Xi, RgM_mul(rowpermute(Tp,y), Ai));
          simplesp = vec_append(simplesp, mkvec3(AAi, Tpi, Pi));
          if ((dimsimple += degpol(Pi)) == newd) goto END;
        }
        else
        {
          Xi = RgM_mul(Xi, rowpermute(X,y));
          nextsp = vec_append(nextsp, mkvec2(AAi, Xi));
        }
      }
    }
    todosp = nextsp; if (lg(todosp) == 1) break;
  }
END:
  if (DEBUGLEVEL) err_printf("end split, need to clean\n");
  return mfspclean(mf, mf0, NF, ord, sort_by_dim(simplesp), flag);
}
static GEN
dim_filter(GEN v, long dim)
{
  GEN P = gel(v,2);
  long j, l = lg(P);
  for (j = 1; j < l; j++)
    if (degpol(gel(P,j)) > dim)
    {
      v = mkvec2(vecslice(gel(v,1),1,j-1), vecslice(P,1,j-1));
      break;
    }
  return v;
}
static long
dim_sum(GEN v)
{
  GEN P = gel(v,2);
  long j, l = lg(P), d = 0;
  for (j = 1; j < l; j++) d += degpol(gel(P,j));
  return d;
}
static GEN
split_i(GEN mf, long dimlim, long flag)
{ long junk; return split_ii(mf, dimlim, flag, NULL, &junk); }
/* mf is either already split or output by mfinit. Splitting is done only for
 * newspace except in weight 1. If flag = 0 (default) split completely.
 * If flag = d > 0, only give the Galois polynomials in degree > d
 * Flag is ignored if dimlim = 1. */
GEN
mfsplit(GEN mf0, long dimlim, long flag)
{
  pari_sp av = avma;
  GEN v, mf = checkMF_i(mf0);
  if (!mf) pari_err_TYPE("mfsplit", mf0);
  if ((v = obj_check(mf, MF_SPLIT)))
  { if (dimlim) v = dim_filter(v, dimlim); }
  else if (dimlim && (v = obj_check(mf, MF_SPLITN)))
  { v = (itos(gel(v,1)) >= dimlim)? dim_filter(gel(v,2), dimlim): NULL; }
  if (!v)
  {
    long newd;
    v = split_ii(mf, dimlim, flag, NULL, &newd);
    if (lg(v) == 1) obj_insert(mf, MF_SPLITN, mkvec2(utoi(dimlim), v));
    else if (!flag)
    {
      if (dim_sum(v) == newd) obj_insert(mf, MF_SPLIT,v);
      else obj_insert(mf, MF_SPLITN, mkvec2(utoi(dimlim), v));
    }
  }
  return gerepilecopy(av, v);
}
static GEN
split(GEN mf) { return split_i(mf,0,0); }
GEN
MF_get_newforms(GEN mf) { return gel(obj_checkbuild(mf,MF_SPLIT,&split),1); }
GEN
MF_get_fields(GEN mf) { return gel(obj_checkbuild(mf,MF_SPLIT,&split),2); }

/*************************************************************************/
/*                     Modular forms of Weight 1                         */
/*************************************************************************/
/* S_1(G_0(N)), small N. Return 1 if definitely empty; return 0 if maybe
 * nonempty  */
static int
wt1empty(long N)
{
  if (N <= 100) switch (N)
  { /* nonempty [32/100] */
    case 23: case 31: case 39: case 44: case 46:
    case 47: case 52: case 55: case 56: case 57:
    case 59: case 62: case 63: case 68: case 69:
    case 71: case 72: case 76: case 77: case 78:
    case 79: case 80: case 83: case 84: case 87:
    case 88: case 92: case 93: case 94: case 95:
    case 99: case 100: return 0;
    default: return 1;
  }
  if (N <= 600) switch(N)
  { /* empty [111/500] */
    case 101: case 102: case 105: case 106: case 109:
    case 113: case 121: case 122: case 123: case 125:
    case 130: case 134: case 137: case 146: case 149:
    case 150: case 153: case 157: case 162: case 163:
    case 169: case 170: case 173: case 178: case 181:
    case 182: case 185: case 187: case 193: case 194:
    case 197: case 202: case 205: case 210: case 218:
    case 221: case 226: case 233: case 241: case 242:
    case 245: case 246: case 250: case 257: case 265:
    case 267: case 269: case 274: case 277: case 281:
    case 289: case 293: case 298: case 305: case 306:
    case 313: case 314: case 317: case 326: case 337:
    case 338: case 346: case 349: case 353: case 361:
    case 362: case 365: case 369: case 370: case 373:
    case 374: case 377: case 386: case 389: case 394:
    case 397: case 401: case 409: case 410: case 421:
    case 425: case 427: case 433: case 442: case 449:
    case 457: case 461: case 466: case 481: case 482:
    case 485: case 490: case 493: case 509: case 514:
    case 521: case 530: case 533: case 534: case 538:
    case 541: case 545: case 554: case 557: case 562:
    case 565: case 569: case 577: case 578: case 586:
    case 593: return 1;
    default: return 0;
  }
  return 0;
}

static GEN
initwt1trace(GEN mf)
{
  GEN S = MF_get_S(mf), v, H;
  long l, i;
  if (lg(S) == 1) return mftrivial();
  H = mfheckemat(mf, Mindex_as_coef(mf));
  l = lg(H); v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(v,i) = gtrace(gel(H,i));
  v = Minv_RgC_mul(MF_get_Minv(mf), v);
  return mflineardiv_linear(S, v, 1);
}
static GEN
initwt1newtrace(GEN mf)
{
  GEN v, D, S, Mindex, CHI = MF_get_CHI(mf);
  long FC, lD, i, sb, N1, N2, lM, N = MF_get_N(mf);
  CHI = mfchartoprimitive(CHI, &FC);
  if (N % FC || mfcharparity(CHI) == 1) return mftrivial();
  D = mydivisorsu(N/FC); lD = lg(D);
  S = MF_get_S(mf);
  if (lg(S) == 1) return mftrivial();
  N2 = newd_params2(N);
  N1 = N / N2;
  Mindex = MF_get_Mindex(mf);
  lM = lg(Mindex);
  sb = Mindex[lM-1];
  v = zerovec(sb+1);
  for (i = 1; i < lD; i++)
  {
    long M = FC*D[i], j;
    GEN tf = initwt1trace(M == N? mf: mfinit_Nkchi(M, 1, CHI, mf_CUSP, 0));
    GEN listd, w;
    if (mf_get_type(tf) == t_MF_CONST) continue;
    w = mfcoefs_i(tf, sb, 1);
    if (M == N) { v = gadd(v, w); continue; }
    listd = mydivisorsu(u_ppo(ugcd(N/M, N1), FC));
    for (j = 1; j < lg(listd); j++)
    {
      long d = listd[j], d2 = d*d; /* coprime to FC */
      GEN dk = mfchareval(CHI, d);
      long NMd = N/(M*d), m;
      for (m = 1; m <= sb/d2; m++)
      {
        long be = mubeta2(NMd, m);
        if (be)
        {
          GEN c = gmul(dk, gmulsg(be, gel(w, m+1)));
          long n = m*d2;
          gel(v, n+1) = gadd(gel(v, n+1), c);
        }
      }
    }
  }
  if (gequal0(gel(v,2))) return mftrivial();
  v = vecpermute(v,Mindex);
  v = Minv_RgC_mul(MF_get_Minv(mf), v);
  return mflineardiv_linear(S, v, 1);
}

/* i*p + 1, i*p < lim corresponding to a_p(f_j), a_{2p}(f_j)...  */
static GEN
pindices(long p, long lim)
{
  GEN v = cgetg(lim, t_VECSMALL);
  long i, ip;
  for (i = 1, ip = p + 1; ip < lim; i++, ip += p) v[i] = ip;
  setlg(v, i); return v;
}

/* assume !wt1empty(N), in particular N>25 */
/* Returns [[lim,p], mf (weight 2), p*lim x dim matrix] */
static GEN
mf1_pre(long N)
{
  pari_timer tt;
  GEN mf, v, L, I, M, Minv, den;
  long B, lim, LIM, p;

  if (DEBUGLEVEL) timer_start(&tt);
  mf = mfinit_Nkchi(N, 2, mfchartrivial(), mf_CUSP, 0);
  if (DEBUGLEVEL)
    timer_printf(&tt, "mf1basis [pre]: S_2(%ld), dim = %ld",
                 N, MF_get_dim(mf));
  M = MF_get_M(mf); Minv = MF_get_Minv(mf); den = gel(Minv,2);
  B = mfsturm_mf(mf);
  if (uisprime(N))
  {
    lim = 2 * MF_get_dim(mf); /* ensure mfstabiter's first kernel ~ square */
    p = 2;
  }
  else
  {
    forprime_t S;
    u_forprime_init(&S, 2, N);
    while ((p = u_forprime_next(&S)))
      if (N % p) break;
    lim = B + 1;
  }
  LIM = (N & (N - 1))? 2 * lim: 3 * lim; /* N power of 2 ? */
  L = mkvecsmall4(lim, LIM, mfsturmNk(N,1), p);
  M = bhnmat_extend_nocache(M, N, LIM-1, 1, MF_get_S(mf));
  if (DEBUGLEVEL) timer_printf(&tt, "mf1basis [pre]: bnfmat_extend");
  v = pindices(p, LIM);
  if (!LIM) return mkvec4(L, mf, M, v);
  I = RgM_Rg_div(ZM_mul(rowslice(M, B+2, LIM), gel(Minv,1)), den);
  I = Q_remove_denom(I, &den);
  if (DEBUGLEVEL) timer_printf(&tt, "mf1basis [prec]: Iden");
  return mkvec5(L, mf, M, v, mkvec2(I, den));
}

/* lg(A) > 1, E a t_POL */
static GEN
mfmatsermul(GEN A, GEN E)
{
  long j, l = lg(A), r = nbrows(A);
  GEN M = cgetg(l, t_MAT);
  E = RgXn_red_shallow(E, r+1);
  for (j = 1; j < l; j++)
  {
    GEN c = RgV_to_RgX(gel(A,j), 0);
    gel(M, j) = RgX_to_RgC(RgXn_mul(c, E, r+1), r);
  }
  return M;
}
/* lg(Ap) > 1, Ep an Flxn */
static GEN
mfmatsermul_Fl(GEN Ap, GEN Ep, ulong p)
{
  long j, l = lg(Ap), r = nbrows(Ap);
  GEN M = cgetg(l, t_MAT);
  for (j = 1; j < l; j++)
  {
    GEN c = Flv_to_Flx(gel(Ap,j), 0);
    gel(M,j) = Flx_to_Flv(Flxn_mul(c, Ep, r+1, p), r);
  }
  return M;
}

/* CHI mod F | N, return mfchar of modulus N.
 * FIXME: wasteful, G should be precomputed  */
static GEN
mfcharinduce(GEN CHI, long N)
{
  GEN G, chi;
  if (mfcharmodulus(CHI) == N) return CHI;
  G = znstar0(utoipos(N), 1);
  chi = zncharinduce(gel(CHI,1), gel(CHI,2), G);
  CHI = leafcopy(CHI);
  gel(CHI,1) = G;
  gel(CHI,2) = chi; return CHI;
}

static GEN
gmfcharno(GEN CHI)
{
  GEN G = gel(CHI,1), chi = gel(CHI,2);
  return mkintmod(znconreyexp(G, chi), znstar_get_N(G));
}
static long
mfcharno(GEN CHI)
{
  GEN n = znconreyexp(gel(CHI,1), gel(CHI,2));
  return itou(n);
}

/* return k such that minimal mfcharacter in Galois orbit of CHI is CHI^k */
static long
mfconreyminimize(GEN CHI)
{
  GEN G = gel(CHI,1), cyc, chi;
  cyc = ZV_to_zv(znstar_get_cyc(G));
  chi = ZV_to_zv(znconreychar(G, gel(CHI,2)));
  return zv_cyc_minimize(cyc, chi, coprimes_zv(mfcharorder(CHI)));
}

/* find scalar c such that first nonzero entry of c*v is 1; return c*v */
static GEN
RgV_normalize(GEN v, GEN *pc)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN c = gel(v,i);
    if (!gequal0(c))
    {
      if (gequal1(c)) break;
      *pc = ginv(c); return RgV_Rg_mul(v, *pc);
    }
  }
  *pc = gen_1; return v;
}
/* pS != NULL; dim > 0 */
static GEN
mftreatdihedral(long N, GEN DIH, GEN POLCYC, long ordchi, GEN *pS)
{
  long l = lg(DIH), lim = mfsturmNk(N, 1), i;
  GEN Minv, C = cgetg(l, t_VEC), M = cgetg(l, t_MAT);
  for (i = 1; i < l; i++)
  {
    GEN c, v = mfcoefs_i(gel(DIH,i), lim, 1);
    gel(M,i) = RgV_normalize(v, &c);
    gel(C,i) = Rg_col_ei(c, l-1, i);
  }
  Minv = gel(mfclean(M,POLCYC,ordchi,0),2);
  M = RgM_Minv_mul(M, Minv);
  C = RgM_Minv_mul(C, Minv);
  *pS = vecmflinear(DIH, C); return M;
}

/* same mode a maximal ideal above q */
static GEN
Tpmod(GEN Ap, GEN A, ulong chip, long p, ulong q)
{
  GEN B = leafcopy(Ap);
  long i, ip, l = lg(B);
  for (i = 1, ip = p; ip < l; i++, ip += p)
    B[ip] = Fl_add(B[ip], Fl_mul(A[i], chip, q), q);
  return B;
}
/* Tp(f_1), ..., Tp(f_d) mod q */
static GEN
matTpmod(GEN Ap, GEN A, ulong chip, long p, ulong q)
{
  long i, l;
  GEN B = cgetg_copy(A, &l);
  for (i = 1; i < l; i++) gel(B,i) = Tpmod(gel(Ap,i), gel(A,i), chip, p, q);
  return B;
}

/* Ap[i] = a_{p*i}(F), A[i] = a_i(F), i = 1..lim
 * Tp(f)[n] = a_{p*n}(f) + chi(p) a_{n/p}(f) * 1_{p | n} */
static GEN
Tp(GEN Ap, GEN A, GEN chip, long p)
{
  GEN B = leafcopy(Ap);
  long i, ip, l = lg(B);
  for (i = 1, ip = p; ip < l; i++, ip += p)
    gel(B,ip) = gadd(gel(B,ip), gmul(gel(A,i), chip));
  return B;
}
/* Tp(f_1), ..., Tp(f_d) mod q */
static GEN
matTp(GEN Ap, GEN A, GEN chip, long p)
{
  long i, l;
  GEN B = cgetg_copy(A, &l);
  for (i = 1; i < l; i++) gel(B,i) = Tp(gel(Ap,i), gel(A,i), chip, p);
  return B;
}

static GEN
_RgXQM_mul(GEN x, GEN y, GEN T)
{ return T? RgXQM_mul(x, y, T): RgM_mul(x, y); }
/* largest T-stable Q(CHI)-subspace of Q(CHI)-vector space spanned by columns
 * of A */
static GEN
mfstabiter(GEN *pC, GEN A0, GEN chip, GEN TMP, GEN P, long ordchi)
{
  GEN A, Ap, vp = gel(TMP,4), C = NULL;
  long i, lA, lim1 = gel(TMP,1)[3], p = gel(TMP,1)[4];
  pari_timer tt;

  Ap = rowpermute(A0, vp);
  A = rowslice(A0, 2, nbrows(Ap)+1); /* remove a0 */
  for(;;)
  {
    GEN R = shallowconcat(matTp(Ap, A, chip, p), A);
    GEN B = QabM_ker(R, P, ordchi);
    long lB = lg(B);
    if (DEBUGLEVEL)
      timer_printf(&tt, "mf1basis: Hecke intersection (dim %ld)", lB-1);
    if (lB == 1) return NULL;
    lA = lg(A); if (lB == lA) break;
    B = rowslice(B, 1, lA-1);
    Ap = _RgXQM_mul(Ap, B, P);
    A = _RgXQM_mul(A, B, P);
    C = C? _RgXQM_mul(C, B, P): B;
  }
  if (nbrows(A) < lim1)
  {
    A0 = rowslice(A0, 2, lim1);
    A = C? _RgXQM_mul(A0, C, P): A0;
  }
  else /* all needed coefs computed */
    A = rowslice(A, 1, lim1-1);
  if (*pC) C = C? _RgXQM_mul(*pC, C, P): *pC;
  /* put back a0 */
  for (i = 1; i < lA; i++) gel(A,i) = vec_prepend(gel(A,i), gen_0);
  *pC = C; return A;
}

static long
mfstabitermod(GEN A, GEN vp, ulong chip, long p, ulong q)
{
  GEN Ap, C = NULL;
  Ap = rowpermute(A, vp);
  A = rowslice(A, 2, nbrows(Ap)+1);
  while (1)
  {
    GEN Rp = shallowconcat(matTpmod(Ap, A, chip, p, q), A);
    GEN B = Flm_ker(Rp, q);
    long lA = lg(A), lB = lg(B);
    if (lB == 1) return 0;
    if (lB == lA) return lA-1;
    B = rowslice(B, 1, lA-1);
    Ap = Flm_mul(Ap, B, q);
    A = Flm_mul(A, B, q);
    C = C? Flm_mul(C, B, q): B;
  }
}

static GEN
mfcharinv_i(GEN CHI)
{
  GEN G = gel(CHI,1);
  CHI = leafcopy(CHI); gel(CHI,2) =  zncharconj(G, gel(CHI,2)); return CHI;
}

/* upper bound dim S_1(Gamma_0(N),chi) performing the linear algebra mod p */
static long
mf1dimmod(GEN E1, GEN E, GEN chip, long ordchi, long dih, GEN TMP)
{
  GEN E1i, A, vp, mf, C = NULL;
  ulong q, r = QabM_init(ordchi, &q);
  long lim, LIM, p;

  LIM = gel(TMP,1)[2]; lim = gel(TMP,1)[1];
  mf= gel(TMP,2);
  A = gel(TMP,3);
  A = QabM_to_Flm(A, r, q);
  E1 = QabX_to_Flx(E1, r, q);
  E1i = Flxn_inv(E1, nbrows(A), q);
  if (E)
  {
    GEN Iden = gel(TMP,5), I = gel(Iden,1), den = gel(Iden,2);
    GEN Mindex = MF_get_Mindex(mf), F = rowslice(A, 1, LIM);
    GEN E1ip = Flxn_red(E1i, LIM);
    ulong d = den? umodiu(den, q): 1;
    long i, nE = lg(E) - 1;
    pari_sp av;

    I = ZM_to_Flm(I, q);
    if (d != 1) I = Flm_Fl_mul(I, Fl_inv(d, q), q);
    av = avma;
    for (i = 1; i <= nE; i++)
    {
      GEN e = Flxn_mul(E1ip, QabX_to_Flx(gel(E,i), r, q), LIM, q);
      GEN B = mfmatsermul_Fl(F, e, q), z;
      GEN B2 = Flm_mul(I, rowpermute(B, Mindex), q);
      B = rowslice(B, lim+1,LIM);
      z = Flm_ker(Flm_sub(B2, B, q), q);
      if (lg(z)-1 == dih) return dih;
      C = C? Flm_mul(C, z, q): z;
      F = Flm_mul(F, z, q);
      gerepileall(av, 2, &F,&C);
    }
    A = F;
  }
  /* use Schaeffer */
  p = gel(TMP,1)[4]; vp = gel(TMP,4);
  A = mfmatsermul_Fl(A, E1i, q);
  return mfstabitermod(A, vp, Qab_to_Fl(chip, r, q), p, q);
}

static GEN
mf1intermat(GEN A, GEN Mindex, GEN e, GEN Iden, long lim, GEN POLCYC)
{
  long j, l = lg(A), LIM = nbrows(A);
  GEN I = gel(Iden,1), den = gel(Iden,2), B = cgetg(l, t_MAT);

  for (j = 1; j < l; j++)
  {
    pari_sp av = avma;
    GEN c = RgV_to_RgX(gel(A,j), 0), c1, c2;
    c = RgX_to_RgC(RgXn_mul(c, e, LIM), LIM);
    if (POLCYC) c = liftpol_shallow(c);
    c1 = vecslice(c, lim+1, LIM);
    if (den) c1 = RgC_Rg_mul(c1, den);
    c2 = RgM_RgC_mul(I, vecpermute(c, Mindex));
    gel(B, j) = gerepileupto(av, RgC_sub(c2, c1));
  }
  return B;
}
/* Compute the full S_1(\G_0(N),\chi); return NULL if space is empty; else
 * if pS is NULL, return stoi(dim), where dim is the dimension; else *pS is
 * set to a vector of forms making up a basis, and return the matrix of their
 * Fourier expansions. pdih gives the dimension of the subspace generated by
 * dihedral forms; TMP is from mf1_pre or NULL. */
static GEN
mf1basis(long N, GEN CHI, GEN TMP, GEN vSP, GEN *pS, long *pdih)
{
  GEN E = NULL, EB, E1, E1i, dE1i, mf, A, C, POLCYC, DIH, Minv, chip;
  long nE = 0, p, LIM, lim, lim1, i, lA, dimp, ordchi, dih;
  pari_timer tt;
  pari_sp av;

  if (pdih) *pdih = 0;
  if (pS) *pS = NULL;
  if (wt1empty(N) || mfcharparity(CHI) != -1) return NULL;
  ordchi = mfcharorder(CHI);
  if (uisprime(N) && ordchi > 4) return NULL;
  if (pS)
  {
    DIH = mfdihedralcusp(N, CHI, vSP);
    dih = lg(DIH) - 1;
  }
  else
  {
    DIH = NULL;
    dih = mfdihedralcuspdim(N, CHI, vSP);
  }
  POLCYC = (ordchi <= 2)? NULL: mfcharpol(CHI);
  if (pdih) *pdih = dih;
  if (N <= 600) switch(N)
  {
    long m;
    case 219: case 273: case 283: case 331: case 333: case 344: case 416:
    case 438: case 468: case 491: case 504: case 546: case 553: case 563:
    case 566: case 581: case 592:
      break; /* one chi with both exotic and dihedral forms */
    default: /* only dihedral forms */
      if (!dih) return NULL;
      /* fall through */
    case 124: case 133: case 148: case 171: case 201: case 209: case 224:
    case 229: case 248: case 261: case 266: case 288: case 296: case 301:
    case 309: case 325: case 342: case 371: case 372: case 380: case 399:
    case 402: case 403: case 404: case 408: case 418: case 432: case 444:
    case 448: case 451: case 453: case 458: case 496: case 497: case 513:
    case 522: case 527: case 532: case 576: case 579:
      /* no chi with both exotic and dihedral; one chi with exotic forms */
      if (dih)
      {
        if (!pS) return utoipos(dih);
        return mftreatdihedral(N, DIH, POLCYC, ordchi, pS) ;
      }
      m = mfcharno(mfcharinduce(CHI,N));
      if (N == 124 && (m != 67 && m != 87)) return NULL;
      if (N == 133 && (m != 83 && m !=125)) return NULL;
      if (N == 148 && (m !=105 && m !=117)) return NULL;
      if (N == 171 && (m != 94 && m !=151)) return NULL;
      if (N == 201 && (m != 29 && m !=104)) return NULL;
      if (N == 209 && (m != 87 && m !=197)) return NULL;
      if (N == 224 && (m != 95 && m !=191)) return NULL;
      if (N == 229 && (m !=107 && m !=122)) return NULL;
      if (N == 248 && (m != 87 && m !=191)) return NULL;
      if (N == 261 && (m != 46 && m !=244)) return NULL;
      if (N == 266 && (m != 83 && m !=125)) return NULL;
      if (N == 288 && (m != 31 && m !=223)) return NULL;
      if (N == 296 && (m !=105 && m !=265)) return NULL;
  }
  if (DEBUGLEVEL)
    err_printf("mf1basis: start character %Ps, conductor = %ld, order = %ld\n",
               gmfcharno(CHI), mfcharconductor(CHI), ordchi);
  if (!TMP) TMP = mf1_pre(N);
  lim = gel(TMP,1)[1]; LIM = gel(TMP,1)[2]; lim1 = gel(TMP,1)[3];
  p = gel(TMP,1)[4];
  mf  = gel(TMP,2);
  A   = gel(TMP,3);
  EB = mfeisensteinbasis(N, 1, mfcharinv_i(CHI));
  nE = lg(EB) - 1;
  E1 = RgV_to_RgX(mftocol(gel(EB,1), LIM-1, 1), 0); /* + O(x^LIM) */
  if (--nE)
    E = RgM_to_RgXV(mfvectomat(vecslice(EB, 2, nE+1), LIM-1, 1), 0);
  chip = mfchareval(CHI, p); /* != 0 */
  if (DEBUGLEVEL) timer_start(&tt);
  av = avma; dimp = mf1dimmod(E1, E, chip, ordchi, dih, TMP);
  set_avma(av);
  if (DEBUGLEVEL) timer_printf(&tt, "mf1basis: dim mod p is %ld", dimp);
  if (!dimp) return NULL;
  if (!pS) return utoi(dimp);
  if (dimp == dih) return mftreatdihedral(N, DIH, POLCYC, ordchi, pS);
  E1i = RgXn_inv(E1, LIM); /* E[1] does not vanish at oo */
  if (POLCYC) E1i = liftpol_shallow(E1i);
  E1i = Q_remove_denom(E1i, &dE1i);
  if (DEBUGLEVEL)
  {
    GEN a0 = gel(E1,2);
    if (typ(a0) == t_POLMOD) a0 = gnorm(a0);
    a0 = Q_abs_shallow(a0);
    timer_printf(&tt, "mf1basis: invert E; norm(a0(E)) = %Ps", a0);
  }
  C = NULL;
  if (nE)
  { /* mf attached to S2(N), fi = mfbasis(mf)
     * M = coefs(f1,...,fd) up to LIM
     * F = coefs(F1,...,FD) = M * C, for some matrix C over Q(chi),
     * initially 1, eventually giving \cap_E S2 / E; D <= d.
     * B = coefs(E/E1 F1, .., E/E1 FD); we want X in Q(CHI)^d and
     * Y in Q(CHI)^D such that
     *   B * X = M * Y, i.e. Minv * rowpermute(B, Mindex * X) = Y
     *(B  - I * rowpermute(B, Mindex)) * X = 0.
     * where I = M * Minv. Rows of (B - I * ...) are 0 up to lim so
     * are not included */
    GEN Mindex = MF_get_Mindex(mf), Iden  = gel(TMP,5);
    pari_timer TT;
    pari_sp av = avma;
    if (DEBUGLEVEL) timer_start(&TT);
    for (i = 1; i <= nE; i++)
    {
      pari_sp av2 = avma;
      GEN e, z, B;

      e = Q_primpart(RgXn_mul(E1i, gel(E,i), LIM));
      if (DEBUGLEVEL) timer_printf(&TT, "mf1basis: E[%ld] / E[1]", i+1);
      /* the first time A is over Z and it is more efficient to lift than
         * to let RgXn_mul use Kronecker's trick */
      if (POLCYC && i == 1) e = liftpol_shallow(e);
      B = mf1intermat(A, Mindex, e, Iden, lim, i == 1? NULL: POLCYC);
      if (DEBUGLEVEL) timer_printf(&TT, "mf1basis: ... intermat");
      z = gerepileupto(av2, QabM_ker(B, POLCYC, ordchi));
      if (DEBUGLEVEL)
        timer_printf(&TT, "mf1basis: ... kernel (dim %ld)",lg(z)-1);
      if (lg(z) == 1) return NULL;
      if (lg(z) == lg(A)) { set_avma(av2); continue; } /* no progress */
      C = C? _RgXQM_mul(C, z, POLCYC): z;
      A = _RgXQM_mul(A, z, POLCYC);
      if (DEBUGLEVEL) timer_printf(&TT, "mf1basis: ... updates");
      if (lg(z)-1 == dimp) break;
      if (gc_needed(av, 1))
      {
        if (DEBUGMEM > 1) pari_warn(warnmem,"mf1basis i = %ld", i);
        gerepileall(av, 2, &A, &C);
      }
    }
    if (DEBUGLEVEL) timer_printf(&tt, "mf1basis: intersection [total]");
  }
  lA = lg(A);
  if (lA-1 == dimp)
  {
    A = mfmatsermul(rowslice(A, 1, lim1), E1i);
    if (POLCYC) A = RgXQM_red(A, POLCYC);
    if (DEBUGLEVEL) timer_printf(&tt, "mf1basis: matsermul [1]");
  }
  else
  {
    A = mfmatsermul(A, E1i);
    if (POLCYC) A = RgXQM_red(A, POLCYC);
    if (DEBUGLEVEL) timer_printf(&tt, "mf1basis: matsermul [2]");
    A = mfstabiter(&C, A, chip, TMP, POLCYC, ordchi);
    if (DEBUGLEVEL) timer_printf(&tt, "mf1basis: Hecke stability");
    if (!A) return NULL;
  }
  if (dE1i) C = RgM_Rg_mul(C, dE1i);
  if (POLCYC)
  {
    A = QXQM_to_mod_shallow(A, POLCYC);
    C = QXQM_to_mod_shallow(C, POLCYC);
  }
  lA = lg(A);
  for (i = 1; i < lA; i++)
  {
    GEN c, v = gel(A,i);
    gel(A,i) = RgV_normalize(v, &c);
    gel(C,i) = RgC_Rg_mul(gel(C,i), c);
  }
  Minv = gel(mfclean(A, POLCYC, ordchi, 0), 2);
  A = RgM_Minv_mul(A, Minv);
  C = RgM_Minv_mul(C, Minv);
  *pS = vecmflineardiv0(MF_get_S(mf), C, gel(EB,1));
  return A;
}

static void
MF_set_space(GEN mf, long x) { gmael(mf,1,4) = utoi(x); }
static GEN
mf1_cusptonew(GEN mf, GEN vSP)
{
  const long vy = 1;
  long i, lP, dSnew, ct;
  GEN vP, F, S, Snew, vF, v = split_ii(mf, 0, 0, vSP, &i);

  F = gel(v,1);
  vP= gel(v,2); lP = lg(vP);
  if (lP == 1) { obj_insert(mf, MF_SPLIT, v); return NULL; }
  MF_set_space(mf, mf_NEW);
  S = MF_get_S(mf);
  dSnew = dim_sum(v);
  Snew = cgetg(dSnew + 1, t_VEC); ct = 0;
  vF = cgetg(lP, t_MAT);
  for (i = 1; i < lP; i++)
  {
    GEN V, P = gel(vP,i), f = liftpol_shallow(gel(F,i));
    long j, d = degpol(P);
    gel(vF,i) = V = zerocol(dSnew);
    if (d == 1)
    {
      gel(Snew, ct+1) = mflineardiv_linear(S, f, 0);
      gel(V, ct+1) = gen_1;
    }
    else
    {
      f = RgXV_to_RgM(f,d);
      for (j = 1; j <= d; j++)
      {
        gel(Snew, ct+j) = mflineardiv_linear(S, row(f,j), 0);
        gel(V, ct+j) = mkpolmod(pol_xn(j-1,vy), P);
      }
    }
    ct += d;
  }
  obj_insert(mf, MF_SPLIT, mkvec2(vF, vP));
  gel(mf,3) = Snew; return mf;
}
static GEN
mf1init(long N, GEN CHI, GEN TMP, GEN vSP, long space, long flraw)
{
  GEN mf, mf1, S, M = mf1basis(N, CHI, TMP, vSP, &S, NULL);
  if (!M) return NULL;
  mf1 = mkvec4(stoi(N), gen_1, CHI, utoi(mf_CUSP));
  mf = mkmf(mf1, cgetg(1,t_VEC), S, gen_0, NULL);
  if (space == mf_NEW)
  {
    gel(mf,5) = mfcleanCHI(M,CHI, 0);
    mf = mf1_cusptonew(mf, vSP); if (!mf) return NULL;
    if (!flraw) M = mfcoefs_mf(mf, mfsturmNk(N,1)+1, 1);
  }
  gel(mf,5) = flraw? zerovec(3): mfcleanCHI(M, CHI, 0);
  return mf;
}

static GEN
mfEMPTY(GEN mf1)
{
  GEN Minv = mkMinv(cgetg(1,t_MAT), NULL,NULL,NULL);
  GEN M = mkvec3(cgetg(1,t_VECSMALL), Minv, cgetg(1,t_MAT));
  return mkmf(mf1, cgetg(1,t_VEC), cgetg(1,t_VEC), cgetg(1,t_VEC), M);
}
static GEN
mfEMPTYall(long N, GEN gk, GEN vCHI, long space)
{
  long i, l;
  GEN v, gN, gs;
  if (!vCHI) return cgetg(1, t_VEC);
  gN = utoipos(N); gs = utoi(space);
  l = lg(vCHI); v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(v,i) = mfEMPTY(mkvec4(gN,gk,gel(vCHI,i),gs));
  return v;
}

static GEN
fmt_dim(GEN CHI, long d, long dih)
{ return mkvec4(gmfcharorder(CHI), gmfcharno(CHI), utoi(d), stoi(dih)); }
/* merge two vector of fmt_dim's for the same vector of characters. If CHI
 * is not NULL, remove dim-0 spaces and add character from CHI */
static GEN
merge_dims(GEN V, GEN W, GEN CHI)
{
  long i, j, id, l = lg(V);
  GEN A = cgetg(l, t_VEC);
  if (l == 1) return A;
  id = CHI? 1: 3;
  for (i = j = 1; i < l; i++)
  {
    GEN v = gel(V,i), w = gel(W,i);
    long dv = itou(gel(v,id)), dvh = itou(gel(v,id+1)), d;
    long dw = itou(gel(w,id)), dwh = itou(gel(w,id+1));
    d = dv + dw;
    if (d || CHI)
      gel(A,j++) = CHI? fmt_dim(gel(CHI,i),d, dvh+dwh)
                      : mkvec2s(d,dvh+dwh);
  }
  setlg(A, j); return A;
}
static GEN
mfdim0all(GEN w)
{
  if (w) retconst_vec(lg(w)-1, zerovec(2));
  return cgetg(1,t_VEC);
}
static long
mf1cuspdim_i(long N, GEN CHI, GEN TMP, GEN vSP, long *dih)
{
  pari_sp av = avma;
  GEN b = mf1basis(N, CHI, TMP, vSP, NULL, dih);
  return gc_long(av, b? itou(b): 0);
}

static long
mf1cuspdim(long N, GEN CHI, GEN vSP)
{
  if (!vSP) vSP = get_vDIH(N, divisorsNF(N, mfcharconductor(CHI)));
  return mf1cuspdim_i(N, CHI, NULL, vSP, NULL);
}
static GEN
mf1cuspdimall(long N, GEN vCHI)
{
  GEN z, TMP, w, vSP;
  long i, j, l;
  if (wt1empty(N)) return mfdim0all(vCHI);
  w = mf1chars(N,vCHI);
  l = lg(w); if (l == 1) return cgetg(1,t_VEC);
  z = cgetg(l, t_VEC);
  TMP = mf1_pre(N); vSP = get_vDIH(N, NULL);
  for (i = j = 1; i < l; i++)
  {
    GEN CHI = gel(w,i);
    long dih, d = mf1cuspdim_i(N, CHI, TMP, vSP, &dih);
    if (vCHI)
      gel(z,j++) = mkvec2s(d, dih);
    else if (d)
      gel(z,j++) = fmt_dim(CHI, d, dih);
  }
  setlg(z,j); return z;
}

/* dimension of S_1(Gamma_1(N)) */
static long
mf1cuspdimsum(long N)
{
  pari_sp av = avma;
  GEN v = mf1cuspdimall(N, NULL);
  long i, ct = 0, l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN w = gel(v,i); /* [ord(CHI),*,dim,*] */
    ct += itou(gel(w,3))*myeulerphiu(itou(gel(w,1)));
  }
  return gc_long(av,ct);
}

static GEN
mf1newdimall(long N, GEN vCHI)
{
  GEN z, w, vTMP, vSP, fa, P, E;
  long i, c, l, lw, P1;
  if (wt1empty(N)) return mfdim0all(vCHI);
  w = mf1chars(N,vCHI);
  lw = lg(w); if (lw == 1) return cgetg(1,t_VEC);
  vTMP = const_vec(N, NULL);
  vSP = get_vDIH(N, NULL);
  gel(vTMP,N) = mf1_pre(N);
  /* if p || N and p \nmid F(CHI), S_1^new(G0(N),chi) = 0 */
  fa = znstar_get_faN(gmael(w,1,1));
  P = gel(fa,1); l = lg(P);
  E = gel(fa,2);
  for (i = P1 = 1; i < l; i++)
    if (E[i] == 1) P1 *= itou(gel(P,i));
  /* P1 = \prod_{v_p(N) = 1} p */
  z = cgetg(lw, t_VEC);
  for (i = c = 1; i < lw; i++)
  {
    long S, j, l, F, dihnew;
    GEN D, CHI = gel(w,i), CHIP = mfchartoprimitive(CHI,&F);

    S = F % P1? 0: mf1cuspdim_i(N, CHI, gel(vTMP,N), vSP, &dihnew);
    if (!S)
    {
      if (vCHI) gel(z, c++) = zerovec(2);
      continue;
    }
    D = mydivisorsu(N/F); l = lg(D);
    for (j = l-2; j > 0; j--) /* skip last M = N */
    {
      long M = D[j]*F, m, s, dih;
      GEN TMP = gel(vTMP,M);
      if (wt1empty(M) || !(m = mubeta(D[l-j]))) continue; /*m = mubeta(N/M)*/
      if (!TMP) gel(vTMP,M) = TMP = mf1_pre(M);
      s = mf1cuspdim_i(M, CHIP, TMP, vSP, &dih);
      if (s) { S += m * s; dihnew += m * dih; }
    }
    if (vCHI)
      gel(z,c++) = mkvec2s(S, dihnew);
    else if (S)
      gel(z, c++) = fmt_dim(CHI, S, dihnew);
  }
  setlg(z,c); return z;
}

static GEN
mf1olddimall(long N, GEN vCHI)
{
  long i, j, l;
  GEN z, w;
  if (wt1empty(N)) return mfdim0all(vCHI);
  w = mf1chars(N,vCHI);
  l = lg(w); z = cgetg(l, t_VEC);
  for (i = j = 1; i < l; i++)
  {
    GEN CHI = gel(w,i);
    long d = mfolddim(N, 1, CHI);
    if (vCHI)
      gel(z,j++) = mkvec2s(d,d?-1:0);
    else if (d)
      gel(z, j++) = fmt_dim(CHI, d, -1);
  }
  setlg(z,j); return z;
}

static long
mf1olddimsum(long N)
{
  GEN D;
  long N2, i, l, S = 0;
  newd_params(N, &N2); /* will ensure mubeta != 0 */
  D = mydivisorsu(N/N2); l = lg(D);
  for (i = 2; i < l; i++)
  {
    long M = D[l-i]*N2, d = mf1cuspdimsum(M);
    if (d) S -= mubeta(D[i]) * d;
  }
  return S;
}
static long
mf1newdimsum(long N)
{
  long S = mf1cuspdimsum(N);
  return S? S - mf1olddimsum(N): 0;
}

/* return the automorphism of a degree-2 nf */
static GEN
nf2_get_conj(GEN nf)
{
  GEN pol = nf_get_pol(nf);
  return deg1pol_shallow(gen_m1, negi(gel(pol,3)), varn(pol));
}
static int
foo_stable(GEN foo)
{ return lg(foo) != 3 || equalii(gel(foo,1), gel(foo,2)); }

static long
mfisdihedral(GEN vF, GEN DIH)
{
  GEN vG = gel(DIH,1), M = gel(DIH,2), v, G, bnr, w, gen, D, f, nf, tau;
  GEN bnr0 = NULL, f0, f0b, xin, foo;
  long i, l, e, j, L, n;
  if (lg(M) == 1) return 0;
  v = RgM_RgC_invimage(M, vF);
  if (!v) return 0;
  l = lg(v);
  for (i = 1; i < l; i++)
    if (!gequal0(gel(v,i))) break;
  if (i == l) return 0;
  G = gel(vG,i);
  bnr = gel(G,2); D = cyc_get_expo(bnr_get_cyc(bnr));
  w = gel(G,3);
  f = bnr_get_mod(bnr);
  nf = bnr_get_nf(bnr);
  tau = nf2_get_conj(nf);
  f0 = gel(f,1); foo = gel(f,2);
  f0b = galoisapply(nf, tau, f0);
  xin = zv_to_ZV(gel(w,2)); /* xi(bnr.gen[i]) = e(xin[i] / D) */
  if (!foo_stable(foo)) { foo = mkvec2(gen_1, gen_1); bnr0 = bnr; }
  if (!gequal(f0, f0b))
  {
    f0 = idealmul(nf, f0, idealdivexact(nf, f0b, idealadd(nf, f0, f0b)));
    bnr0 = bnr;
  }
  if (bnr0)
  { /* conductor not ambiguous */
    GEN S;
    bnr = Buchray(bnr_get_bnf(bnr), mkvec2(f0, foo), nf_INIT | nf_GEN);
    S = bnrsurjection(bnr, bnr0);
    xin = FpV_red(RgV_RgM_mul(xin, gel(S,1)), D);
    /* still xi(gen[i]) = e(xin[i] / D), for the new generators; D stays
     * the same, not exponent(bnr.cyc) ! */
  }
  gen = bnr_get_gen(bnr); L = lg(gen);
  for (j = 1, e = itou(D); j < L; j++)
  {
    GEN Ng = idealnorm(nf, gel(gen,j));
    GEN a = shifti(gel(xin,j), 1); /* xi(g_j^2) = e(a/D) */
    GEN b = FpV_dotproduct(xin, isprincipalray(bnr,Ng), D);
    GEN m = Fp_sub(a, b, D); /* xi(g_j/g_j^\tau) = e(m/D) */
    e = ugcd(e, itou(m)); if (e == 1) break;
  }
  n = itou(D) / e;
  return n == 1? 4: 2*n;
}

static ulong
myradicalu(ulong n) { return zv_prod(gel(myfactoru(n),1)); }

/* list of fundamental discriminants unramified outside N, with sign s
 * [s = 0 => no sign condition] */
static GEN
mfunram(long N, long s)
{
  long cN = myradicalu(N >> vals(N)), p = 1, m = 1, l, c, i;
  GEN D = mydivisorsu(cN), res;
  l = lg(D);
  if (s == 1) m = 0; else if (s == -1) p = 0;
  res = cgetg(6*l - 5, t_VECSMALL);
  c = 1;
  if (!odd(N))
  { /* d = 1 */
    if (p) res[c++] = 8;
    if (m) { res[c++] =-8; res[c++] =-4; }
  }
  for (i = 2; i < l; i++)
  { /* skip d = 1, done above */
    long d = D[i], d4 = d & 3L; /* d odd, squarefree, d4 = 1 or 3 */
    if (d4 == 1) { if (p) res[c++] = d; }
    else         { if (m) res[c++] =-d; }
    if (!odd(N))
    {
      if (p) { res[c++] = 8*d; if (d4 == 3) res[c++] = 4*d; }
      if (m) { res[c++] =-8*d; if (d4 == 1) res[c++] =-4*d; }
    }
  }
  setlg(res, c); return res;
}

/* Return 1 if F is definitely not S4 type; return 0 on failure. */
static long
mfisnotS4(long N, GEN w)
{
  GEN D = mfunram(N, 0);
  long i, lD = lg(D), lw = lg(w);
  for (i = 1; i < lD; i++)
  {
    long p, d = D[i], ok = 0;
    for (p = 2; p < lw; p++)
      if (w[p] && kross(d,p) == -1) { ok = 1; break; }
    if (!ok) return 0;
  }
  return 1;
}

/* Return 1 if Q(sqrt(5)) \not\subset Q(F), i.e. F is definitely not A5 type;
 * return 0 on failure. */
static long
mfisnotA5(GEN F)
{
  GEN CHI = mf_get_CHI(F), P = mfcharpol(CHI), T, Q;

  if (mfcharorder(CHI) % 5 == 0) return 0;
  T = mf_get_field(F); if (degpol(T) == 1) return 1;
  if (degpol(P) > 1) T = rnfequation(P,T);
  Q = gsubgs(pol_xn(2,varn(T)), 5);
  return (typ(nfisincl(Q, T)) == t_INT);
}

/* v[p+1]^2 / chi(p) - 2 = z + 1/z with z primitive root of unity of order n,
 * return n */
static long
mffindrootof1(GEN v, long p, GEN CHI)
{
  GEN ap = gel(v,p+1), u0, u1, u1k, u2;
  long c = 1;
  if (gequal0(ap)) return 2;
  u0 = gen_2; u1k = u1 = gsubgs(gdiv(gsqr(ap), mfchareval(CHI, p)), 2);
  while (!gequalsg(2, liftpol_shallow(u1))) /* u1 = z^c + z^-c */
  {
    u2 = gsub(gmul(u1k, u1), u0);
    u0 = u1; u1 = u2; c++;
  }
  return c;
}

/* we known that F is not dihedral */
static long
mfgaloistype_i(long N, GEN CHI, GEN F, GEN v)
{
  forprime_t iter;
  long lim = lg(v)-2;
  GEN w = zero_zv(lim);
  pari_sp av;
  ulong p;
  u_forprime_init(&iter, 2, lim);
  av = avma;
  while((p = u_forprime_next(&iter))) if (N%p) switch(mffindrootof1(v, p, CHI))
  {
    case 1: case 2: continue;
    case 3: w[p] = 1; break;
    case 4: return -24; /* S4 */
    case 5: return -60; /* A5 */
    default: pari_err_DOMAIN("mfgaloistype", "form", "not a",
                             strtoGENstr("cuspidal eigenform"), F);
    set_avma(av);
  }
  if (mfisnotS4(N,w) && mfisnotA5(F)) return -12; /* A4 */
  return 0; /* FAILURE */
}

static GEN
mfgaloistype0(long N, GEN CHI, GEN F, GEN DIH, long lim)
{
  pari_sp av = avma;
  GEN vF = mftocol(F, lim, 1);
  long t = mfisdihedral(vF, DIH), bound;
  if (t) return gc_stoi(av,t);
  bound = maxss(200, 5*expu(N)*expu(N));
  for(;;)
  {
    t = mfgaloistype_i(N, CHI, F, vF);
    set_avma(av); if (t) return stoi(t);
    if (lim > bound) return gen_0;
    lim += lim >> 1;
    vF = mfcoefs_i(F,lim,1);
  }
}

/* If f is NULL, give all the galoistypes, otherwise just for f */
/* Return 0 to indicate failure; in this case the type is either -12 or -60,
 * most likely -12. FIXME using the Galois representation. */
GEN
mfgaloistype(GEN NK, GEN f)
{
  pari_sp av = avma;
  GEN CHI, T, F, DIH, SP, mf = checkMF_i(NK);
  long N, k, lL, i, lim, SB;

  if (f && !checkmf_i(f)) pari_err_TYPE("mfgaloistype", f);
  if (mf)
  {
    N = MF_get_N(mf);
    k = MF_get_k(mf);
    CHI = MF_get_CHI(mf);
  }
  else
  {
    checkNK(NK, &N, &k, &CHI, 0);
    mf = f? NULL: mfinit_i(NK, mf_NEW);
  }
  if (k != 1) pari_err_DOMAIN("mfgaloistype", "k", "!=", gen_1, stoi(k));
  SB = mf? mfsturm_mf(mf): mfsturmNk(N,1);
  SP = get_DIH(N);
  DIH = mfdihedralnew(N, CHI, SP);
  lim = lg(DIH) == 1? 200: SB;
  DIH = mkvec2(DIH, mfvectomat(DIH,SB,1));
  if (f) return gerepileuptoint(av, mfgaloistype0(N,CHI, f, DIH, lim));
  F = mfeigenbasis(mf); lL = lg(F);
  T = cgetg(lL, t_VEC);
  for (i=1; i < lL; i++) gel(T,i) = mfgaloistype0(N, CHI, gel(F,i), DIH, lim);
  return gerepileupto(av, T);
}

/******************************************************************/
/*                   Find all dihedral forms.                     */
/******************************************************************/
/* lim >= 2 */
static void
consttabdihedral(long lim) { cache_set(cache_DIH, mfdihedralall(lim)); }

/* a ideal coprime to bnr modulus */
static long
mfdiheval(GEN bnr, GEN w, GEN a)
{
  GEN L, cycn = gel(w,1), chin = gel(w,2);
  long ordmax = cycn[1];
  L = ZV_to_Flv(isprincipalray(bnr,a), ordmax);
  return Flv_dotproduct(chin, L, ordmax);
}

/* A(x^k) mod T = polcyclo(m), 0 <= k < m */
static GEN
Galois(GEN A, long k, GEN T, long m)
{
  GEN B;
  long i, ik, d;
  if (typ(A) != t_POL) return A;
  if (varn(A) != varn(T))
  {
    B = cgetg_copy(A, &d); B[1] = A[1];
    for (i = 2; i < d; i++) gel(B, i) = Galois(gel(A, i), k, T, m);
    return B;
  }
  if ((d = degpol(A)) <= 0) return A;
  B = cgetg(m + 2, t_POL); B[1] = A[1]; gel(B,2) = gel(A,2);
  for (i = 1; i < m; i++) gel(B, i+2) = gen_0;
  for (i = 1, ik = k; i <= d; i++, ik = Fl_add(ik, k, m))
    gel(B, ik + 2) = gel(A, i+2);
  return QX_ZX_rem(normalizepol(B), T);
}
static GEN
vecGalois(GEN v, long k, GEN T, long m)
{
  long i, l;
  GEN w = cgetg_copy(v,&l);
  for (i = 1; i < l; i++) gel(w,i) = Galois(gel(v,i), k, T, m);
  return w;
}

static GEN
fix_pol(GEN S, GEN Pn, int *trace)
{
  if (typ(S) != t_POL) return S;
  S = RgX_rem(S, Pn);
  if (typ(S) == t_POL)
  {
    switch(lg(S))
    {
      case 2: return gen_0;
      case 3: return gel(S,2);
    }
    *trace = 1;
  }
  return S;
}

static GEN
dihan(GEN bnr, GEN w, GEN k0j, long m, ulong lim)
{
  GEN nf = bnr_get_nf(bnr), f = bid_get_ideal(bnr_get_bid(bnr));
  GEN v = zerovec(lim+1), cycn = gel(w,1), Tinit = gel(w,3);
  GEN Pn = gel(Tinit,lg(Tinit)==4? 2: 1);
  long j, ordmax = cycn[1];
  long D = itos(nf_get_disc(nf)), vt = varn(Pn);
  int trace = 0;
  ulong p, n;
  forprime_t T;

  if (!lim) return v;
  gel(v,2) = gen_1;
  u_forprime_init(&T, 2, lim);
  /* fill in prime powers first */
  while ((p = u_forprime_next(&T)))
  {
    GEN vP, vchiP, S;
    long k, lP;
    ulong q, qk;
    if (kross(D,p) >= 0) q = p;
    else if (!(q = umuluu_le(p,p,lim))) continue;
    /* q = Norm P */
    vP = idealprimedec(nf, utoipos(p));
    lP = lg(vP);
    vchiP = cgetg(lP, t_VECSMALL);
    for (j = k = 1; j < lP; j++)
    {
      GEN P = gel(vP,j);
      if (!idealval(nf, f, P)) vchiP[k++] = mfdiheval(bnr,w,P);
    }
    if (k == 1) continue;
    setlg(vchiP, k); lP = k;
    if (lP == 2)
    { /* one prime above p not dividing f */
      long s, s0 = vchiP[1];
      for (qk=q, s = s0;; s = Fl_add(s,s0,ordmax))
      {
        S = Qab_zeta(s, ordmax, vt);
        gel(v, qk+1) = fix_pol(S, Pn, &trace);
        if (!(qk = umuluu_le(qk,q,lim))) break;
      }
    }
    else /* two primes above p not dividing f */
    {
      long s, s0 = vchiP[1], s1 = vchiP[2];
      for (qk=q, k = 1;; k++)
      { /* sum over a,b s.t. Norm( P1^a P2^b ) = q^k, i.e. a+b = k */
        long a;
        GEN S = gen_0;
        for (a = 0; a <= k; a++)
        {
          s = Fl_add(Fl_mul(a, s0, ordmax), Fl_mul(k-a, s1, ordmax), ordmax);
          S = gadd(S, Qab_zeta(s, ordmax, vt));
        }
        gel(v, qk+1) = fix_pol(S, Pn, &trace);
        if (!(qk = umuluu_le(qk,q,lim))) break;
      }
    }
  }
  /* complete with nonprime powers */
  for (n = 2; n <= lim; n++)
  {
    GEN S, fa = myfactoru(n), P = gel(fa, 1), E = gel(fa, 2);
    long q;
    if (lg(P) == 2) continue;
    /* not a prime power */
    q = upowuu(P[1],E[1]);
    S = gmul(gel(v, q + 1), gel(v, n/q + 1));
    gel(v, n+1) = fix_pol(S, Pn, &trace);
  }
  if (trace)
  {
    long k0 = k0j[1], jdeg = k0j[2];
    v = QabV_tracerel(Tinit, jdeg, v); /* Apply Galois Mod(k0, ordw) */
    if (k0 > 1) v = vecGalois(v, k0, gel(Tinit,1), m);
  }
  return v;
}

/* as cyc_normalize for t_VECSMALL cyc */
static GEN
cyc_normalize_zv(GEN cyc)
{
  long i, o = cyc[1], l = lg(cyc); /* > 1 */
  GEN D = cgetg(l, t_VECSMALL);
  D[1] = o; for (i = 2; i < l; i++) D[i] = o / cyc[i];
  return D;
}
/* as char_normalize for t_VECSMALLs */
static GEN
char_normalize_zv(GEN chi, GEN ncyc)
{
  long i, l = lg(chi);
  GEN c = cgetg(l, t_VECSMALL);
  if (l > 1) {
    c[1] = chi[1];
    for (i = 2; i < l; i++) c[i] = chi[i] * ncyc[i];
  }
  return c;
}

static GEN
dihan_bnf(long D)
{
  GEN c = getrand(), bnf;
  setrand(gen_1);
  bnf = Buchall(quadpoly_i(stoi(D)), nf_FORCE, LOWDEFAULTPREC);
  setrand(c);
  return bnf;
}
static GEN
dihan_bnr(GEN bnf, GEN A)
{
  GEN c = getrand(), bnr;
  setrand(gen_1);
  bnr = Buchray(bnf, A, nf_INIT|nf_GEN);
  setrand(c);
  return bnr;
}
/* Hecke xi * (D/.) = Dirichlet chi, return v in Q^r st chi(g_i) = e(v[i]).
 * cycn = cyc_normalize_zv(bnr.cyc), chin = char_normalize_zv(chi,cyc) */
static GEN
bnrchartwist2conrey(GEN chin, GEN cycn, GEN bnrconreyN, GEN kroconreyN)
{
  long l = lg(bnrconreyN), c1 = cycn[1], i;
  GEN v = cgetg(l, t_COL);
  for (i = 1; i < l; i++)
  {
    GEN d = sstoQ(zv_dotproduct(chin, gel(bnrconreyN,i)), c1);
    if (kroconreyN[i] < 0) d = gadd(d, ghalf);
    gel(v,i) = d;
  }
  return v;
}

/* chi(g_i) = e(v[i]) denormalize wrt Conrey generators orders */
static GEN
conreydenormalize(GEN znN, GEN v)
{
  GEN gcyc = znstar_get_conreycyc(znN), w;
  long l = lg(v), i;
  w = cgetg(l, t_COL);
  for (i = 1; i < l; i++)
    gel(w,i) = modii(gmul(gel(v,i), gel(gcyc,i)), gel(gcyc,i));
  return w;
}

static long
Miyake(GEN vchi, GEN gb, GEN cycn)
{
  long i, e = cycn[1], lb = lg(gb);
  GEN v = char_normalize_zv(vchi, cycn);
  for (i = 1; i < lb; i++)
    if ((zv_dotproduct(v, gel(gb,i)) -  v[i]) % e) return 1;
  return 0;
}

/* list of Hecke characters not induced by a Dirichlet character up to Galois
 * conjugation, whose conductor is bnr.cond; cycn = cyc_normalize(bnr.cyc)*/
static GEN
mklvchi(GEN bnr, GEN cycn, GEN gb)
{
  GEN cyc = bnr_get_cyc(bnr), cycsmall = ZV_to_zv(cyc);
  GEN vchi = cyc2elts(cycsmall);
  long ordmax = cycsmall[1], c, i, l;
  l = lg(vchi);
  for (i = c = 1; i < l; i++)
  {
    GEN chi = gel(vchi,i);
    if (!gb || Miyake(chi, gb, cycn)) gel(vchi, c++) = Flv_to_ZV(chi);
  }
  setlg(vchi, c); l = c;
  for (i = 1; i < l; i++)
  {
    GEN chi = gel(vchi,i);
    long n;
    if (!chi) continue;
    for (n = 2; n < ordmax; n++)
      if (ugcd(n, ordmax) == 1)
      {
        GEN tmp = ZV_ZV_mod(gmulsg(n, chi), cyc);
        long j;
        for (j = i+1; j < l; j++)
          if (gel(vchi,j) && gequal(gel(vchi,j), tmp)) gel(vchi,j) = NULL;
      }
  }
  for (i = c = 1; i < l; i++)
  {
    GEN chi = gel(vchi,i);
    if (chi && bnrisconductor(bnr, chi)) gel(vchi, c++) = chi;
  }
  setlg(vchi, c); return vchi;
}

static GEN
get_gb(GEN bnr, GEN con)
{
  GEN gb, g = bnr_get_gen(bnr), nf = bnr_get_nf(bnr);
  long i, l = lg(g);
  gb = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(gb,i) = ZV_to_zv(isprincipalray(bnr, galoisapply(nf, con, gel(g,i))));
  return gb;
}
static GEN
get_bnrconreyN(GEN bnr, GEN znN)
{
  GEN z, g = znstar_get_conreygen(znN);
  long i, l = lg(g);
  z = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(z,i) = ZV_to_zv(isprincipalray(bnr,gel(g,i)));
  return z;
}
/* con = NULL if D > 0 or if D < 0 and id != idcon. */
static GEN
mfdihedralcommon(GEN bnf, GEN id, GEN znN, GEN kroconreyN, long N, long D, GEN con)
{
  GEN bnr = dihan_bnr(bnf, id), cyc = ZV_to_zv( bnr_get_cyc(bnr) );
  GEN bnrconreyN, cycn, cycN, Lvchi, res, P, vT;
  long j, ordmax, l, lc, deghecke, vt;

  lc = lg(cyc); if (lc == 1) return NULL;
  cycn = cyc_normalize_zv(cyc);
  Lvchi = mklvchi(bnr, cycn, con? get_gb(bnr, con): NULL);
  l = lg(Lvchi);
  if (l == 1) return NULL;

  bnrconreyN = get_bnrconreyN(bnr, znN);
  cycN = ZV_to_zv(znstar_get_cyc(znN));
  ordmax = cyc[1];
  vT = const_vec(odd(ordmax)? ordmax << 1: ordmax, NULL);
  vt = fetch_user_var("t");
  P = polcyclo(ordmax, vt);
  gel(vT,ordmax) = Qab_trace_init(ordmax, ordmax, P, P);
  deghecke = myeulerphiu(ordmax);
  res = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
  {
    GEN T, v, vchi = ZV_to_zv(gel(Lvchi,j));
    GEN chi, chin = char_normalize_zv(vchi, cycn);
    long o, vnum, k0, degrel;
    v = bnrchartwist2conrey(chin, cycn, bnrconreyN, kroconreyN);
    o = itou(Q_denom(v));
    T = gel(vT, o);
    if (!T) gel(vT,o) = T = Qab_trace_init(ordmax, o, P, polcyclo(o,vt));
    chi = conreydenormalize(znN, v);
    vnum = itou(znconreyexp(znN, chi));
    chi = ZV_to_zv(znconreychar(znN,chi));
    degrel = deghecke / degpol(gel(T,1));
    k0 = zv_cyc_minimize(cycN, chi, coprimes_zv(o));
    vnum = Fl_powu(vnum, k0, N);
    /* encodes degrel forms: jdeg = 0..degrel-1 */
    gel(res,j) = mkvec3(mkvecsmalln(5, N, k0 % o, vnum, D, degrel),
                        id, mkvec3(cycn,chin,T));
  }
  return res;
}

static long
is_cond(long D, long n)
{
  if (D > 0) return n != 4 || (D&7L) == 1;
  return n != 2 && n != 3 && (n != 4 || (D&7L)!=1);
}
/* Append to v all dihedral weight 1 forms coming from D, if fundamental.
 * level in [l1, l2] */
static void
append_dihedral(GEN v, long D, long l1, long l2)
{
  long Da = labs(D), no, i, numi, ct, min, max;
  GEN bnf, con, vI, resall, arch1, arch2;
  pari_sp av;

  /* min <= Nf <= max */
  max = l2 / Da;
  if (l1 == l2)
  { /* assume Da | l2 */
    min = max;
    if (D > 0 && min < 3) return;
  }
  else /* assume l1 < l2 */
    min = (l1 + Da-1)/Da;
  if (!sisfundamental(D)) return;

  av = avma;
  bnf = dihan_bnf(D);
  con = nf2_get_conj(bnf_get_nf(bnf));
  vI = ideallist(bnf, max);
  numi = 0; for (i = min; i <= max; i++) numi += lg(gel(vI, i)) - 1;
  if (D > 0)
  {
    numi <<= 1;
    arch1 = mkvec2(gen_1,gen_0);
    arch2 = mkvec2(gen_0,gen_1);
  }
  else
    arch1 = arch2 = NULL;
  resall = cgetg(numi+1, t_VEC); ct = 1;
  for (no = min; no <= max; no++) if (is_cond(D, no))
  {
    long N = Da*no, lc, lI;
    GEN I = gel(vI, no), znN = znstar0(utoipos(N), 1), conreyN, kroconreyN;

    conreyN = znstar_get_conreygen(znN); lc = lg(conreyN);
    kroconreyN = cgetg(lc, t_VECSMALL);
    for (i = 1; i < lc; i++) kroconreyN[i] = krosi(D, gel(conreyN, i));
    lI = lg(I);
    for (i = 1; i < lI; i++)
    {
      GEN id = gel(I, i), idcon, z;
      long j;
      if (typ(id) == t_INT) continue;
      idcon = galoisapply(bnf, con, id);
      for (j = i; j < lI; j++)
        if (gequal(idcon, gel(I, j))) { gel(I, j) = gen_0; break; }
      if (D < 0)
      {
        GEN conk = i == j ? con : NULL;
        z = mfdihedralcommon(bnf, id, znN, kroconreyN, N, D, conk);
        if (z) gel(resall, ct++) = z;
      }
      else
      {
        GEN ide;
        ide = mkvec2(id, arch1);
        z = mfdihedralcommon(bnf, ide, znN, kroconreyN, N, D, NULL);
        if (z) gel(resall, ct++) = z;
        if (gequal(idcon,id)) continue;
        ide = mkvec2(id, arch2);
        z = mfdihedralcommon(bnf, ide, znN, kroconreyN, N, D, NULL);
        if (z) gel(resall, ct++) = z;
      }
    }
  }
  if (ct == 1) set_avma(av);
  else
  {
    setlg(resall, ct);
    vectrunc_append(v, gerepilecopy(av, shallowconcat1(resall)));
  }
}

static long
di_N(GEN a) { return gel(a,1)[1]; }
static GEN
mfdihedral(long N)
{
  GEN D = mydivisorsu(N), res = vectrunc_init(2*N);
  long j, l = lg(D);
  for (j = 2; j < l; j++)
  { /* skip d = 1 */
    long d = D[j];
    if (d == 2) continue;
    append_dihedral(res, -d, N,N);
    if (d >= 5 && D[l-j] >= 3) append_dihedral(res, d, N,N); /* Nf >= 3 */
  }
  if (lg(res) > 1) res = shallowconcat1(res);
  return res;
}
/* All primitive dihedral weight 1 forms of leven in [1, N], N > 1 */
static GEN
mfdihedralall(long N)
{
  GEN res = vectrunc_init(2*N), z;
  long D, ct, i;

  for (D = -3; D >= -N; D--) append_dihedral(res, D, 1,N);
  /* Nf >= 3 (GTM 193, prop 3.3.18) */
  for (D = N / 3; D >= 5; D--) append_dihedral(res, D, 1,N);
  ct = lg(res);
  if (ct > 1)
  { /* sort wrt N */
    res = shallowconcat1(res);
    res = vecpermute(res, indexvecsort(res, mkvecsmall(1)));
    ct = lg(res);
  }
  z = const_vec(N, cgetg(1,t_VEC));
  for (i = 1; i < ct;)
  { /* regroup result sharing the same N */
    long n = di_N(gel(res,i)), j = i+1, k;
    GEN v;
    while (j < ct && di_N(gel(res,j)) == n) j++;
    gel(z, n) = v = cgetg(j-i+1, t_VEC);
    for (k = 1; i < j; k++,i++) gel(v,k) = gel(res,i);
  }
  return z;
}

/* return [vF, index], where vecpermute(vF,index) generates dihedral forms
 * for character CHI */
static GEN
mfdihedralnew_i(long N, GEN CHI, GEN SP)
{
  GEN bnf, Tinit, Pm, vf, M, V, NK;
  long Dold, d, ordw, i, SB, c, l, k0, k1, chino, chinoorig, lv;

  lv = lg(SP); if (lv == 1) return NULL;
  CHI = mfcharinduce(CHI,N);
  ordw = mfcharorder(CHI);
  chinoorig = mfcharno(CHI);
  k0 = mfconreyminimize(CHI);
  chino = Fl_powu(chinoorig, k0, N);
  k1 = Fl_inv(k0 % ordw, ordw);
  V = cgetg(lv, t_VEC);
  d = 0;
  for (i = l = 1; i < lv; i++)
  {
    GEN sp = gel(SP,i), T = gel(sp,1);
    if (T[3] != chino) continue;
    d += T[5];
    if (k1 != 1)
    {
      GEN t = leafcopy(T);
      t[3] = chinoorig;
      t[2] = (t[2]*k1) % ordw;
      sp = mkvec4(t, gel(sp,2), gel(sp,3), gel(sp,4));
    }
    gel(V, l++) = sp;
  }
  setlg(V, l); /* dihedral forms of level N and character CHI */
  if (l == 1) return NULL;

  SB = mfsturmNk(N,1) + 1;
  M = cgetg(d+1, t_MAT);
  vf = cgetg(d+1, t_VEC);
  NK = mkNK(N, 1, CHI);
  bnf = NULL; Dold = 0;
  for (i = c = 1; i < l; i++)
  { /* T = [N, k0, conreyno, D, degrel] */
    GEN bnr, Vi = gel(V,i), T = gel(Vi,1), id = gel(Vi,2), w = gel(Vi,3);
    long jdeg, k0i = T[2], D = T[4], degrel = T[5];

    if (D != Dold) { Dold = D; bnf = dihan_bnf(D); }
    bnr = dihan_bnr(bnf, id);
    for (jdeg = 0; jdeg < degrel; jdeg++,c++)
    {
      GEN k0j = mkvecsmall2(k0i, jdeg), an = dihan(bnr, w, k0j, ordw, SB);
      settyp(an, t_COL); gel(M,c) = an;
      gel(vf,c) = tag3(t_MF_DIHEDRAL, NK, bnr, w, k0j);
    }
  }
  Tinit = gmael3(V,1,3,3); Pm = gel(Tinit,1);
  V = QabM_indexrank(M, degpol(Pm)==1? NULL: Pm, ordw);
  return mkvec2(vf,gel(V,2));
}
static long
mfdihedralnewdim(long N, GEN CHI, GEN SP)
{
  pari_sp av = avma;
  GEN S = mfdihedralnew_i(N, CHI, SP);
  return gc_long(av, S? lg(gel(S,2))-1: 0);
}
static GEN
mfdihedralnew(long N, GEN CHI, GEN SP)
{
  pari_sp av = avma;
  GEN S = mfdihedralnew_i(N, CHI, SP);
  if (!S) { set_avma(av); return cgetg(1, t_VEC); }
  return vecpermute(gel(S,1), gel(S,2));
}

static long
mfdihedralcuspdim(long N, GEN CHI, GEN vSP)
{
  pari_sp av = avma;
  GEN D, CHIP;
  long F, i, lD, dim;

  CHIP = mfchartoprimitive(CHI, &F);
  D = mydivisorsu(N/F); lD = lg(D);
  dim = mfdihedralnewdim(N, CHI, gel(vSP,N)); /* d = 1 */
  for (i = 2; i < lD; i++)
  {
    long d = D[i], a = mfdihedralnewdim(N/d, CHIP, gel(vSP, N/d));
    if (a) dim += a * mynumdivu(d);
  }
  return gc_long(av,dim);
}

static GEN
mfbdall(GEN E, long N)
{
  GEN v, D = mydivisorsu(N);
  long i, j, nD = lg(D) - 1, nE = lg(E) - 1;
  v = cgetg(nD*nE + 1, t_VEC);
  for (j = 1; j <= nE; j++)
  {
    GEN Ej = gel(E, j);
    for (i = 0; i < nD; i++) gel(v, i*nE + j) = mfbd_i(Ej, D[i+1]);
  }
  return v;
}
static GEN
mfdihedralcusp(long N, GEN CHI, GEN vSP)
{
  pari_sp av = avma;
  GEN D, CHIP, z;
  long F, i, lD;

  CHIP = mfchartoprimitive(CHI, &F);
  D = mydivisorsu(N/F); lD = lg(D);
  z = cgetg(lD, t_VEC);
  gel(z,1) = mfdihedralnew(N, CHI, gel(vSP,N));
  for (i = 2; i < lD; i++) /* skip 1 */
  {
    GEN LF = mfdihedralnew(N / D[i], CHIP, gel(vSP, N / D[i]));
    gel(z,i) = mfbdall(LF, D[i]);
  }
  return gerepilecopy(av, shallowconcat1(z));
}

/* used to decide between ratlift and comatrix for ZM_inv; ratlift is better
 * when N has many divisors */
static int
abundant(ulong N) { return mynumdivu(N) >= 8; }

/* CHI an mfchar */
static int
cmp_ord(void *E, GEN a, GEN b)
{
  GEN chia = MF_get_CHI(a), chib = MF_get_CHI(b);
  (void)E; return cmpii(gmfcharorder(chia), gmfcharorder(chib));
}
/* mfinit structure.
-- mf[1] contains [N,k,CHI,space],
-- mf[2] contains vector of closures of Eisenstein series, empty if not
   full space.
-- mf[3] contains vector of closures, so #mf[3] = dimension of cusp/new space.
-- mf[4] contains the corresponding indices: either j for T(j)tf if newspace,
   or [M,j,d] for B(d)T(j)tf_M if cuspspace or oldspace.
-- mf[5] contains the matrix M of first coefficients of basis, never cleaned.
 * NK is either [N,k] or [N,k,CHI].
 * mfinit does not do the splitting, only the basis generation. */

/* Set flraw to 1 if do not need mf[5]: no mftobasis etc..., only the
   expansions of the basis elements are needed. */

static GEN
mfinit_Nkchi(long N, long k, GEN CHI, long space, long flraw)
{
  GEN M = NULL, mf = NULL, mf1 = mkvec4(utoi(N), stoi(k), CHI, utoi(space));
  long sb = mfsturmNk(N, k);
  if (k < 0 || badchar(N, k, CHI)) return mfEMPTY(mf1);
  if (k == 0 || space == mf_EISEN) /*nothing*/;
  else if (k == 1)
  {
    switch (space)
    {
      case mf_NEW:
      case mf_FULL:
      case mf_CUSP: mf = mf1init(N, CHI, NULL, get_vDIH(N,NULL), space, flraw);
                    break;
      case mf_OLD: pari_err_IMPL("mfinit in weight 1 for old space");
      default: pari_err_FLAG("mfinit");
    }
  }
  else /* k >= 2 */
  {
    long ord = mfcharorder(CHI);
    GEN z = NULL, P = (ord <= 2)? NULL: mfcharpol(CHI);
    cachenew_t cache;
    switch(space)
    {
      case mf_NEW:
        mf = mfnewinit(N, k, CHI, &cache, 1);
        if (mf && !flraw) { M = MF_get_M(mf); z = MF_get_Mindex(mf); }
        break;
      case mf_OLD:
      case mf_CUSP:
      case mf_FULL:
        if (!(mf = mfinitcusp(N, k, CHI, &cache, space))) break;
        if (!flraw)
        {
          M = bhnmat_extend(M, sb+1, 1, MF_get_S(mf), &cache);
          if (space != mf_FULL) gel(mf,5) = mfcleanCHI(M, CHI, abundant(N));
        }
        dbg_cachenew(&cache); break;
      default: pari_err_FLAG("mfinit");
    }
    if (z) gel(mf,5) = mfclean2(M, z, P, ord);
  }
  if (!mf) mf = mfEMPTY(mf1);
  else
  {
    gel(mf,1) = mf1;
    if (flraw) gel(mf,5) = zerovec(3);
  }
  if (!space_is_cusp(space))
  {
    GEN E = mfeisensteinbasis(N, k, CHI);
    gel(mf,2) = E;
    if (!flraw)
    {
      if (M)
        M = shallowconcat(mfvectomat(E, sb+1, 1), M);
      else
        M = mfcoefs_mf(mf, sb+1, 1);
      gel(mf,5) = mfcleanCHI(M, CHI, abundant(N));
    }
  }
  return mf;
}

/* mfinit for k = nk/dk */
static GEN
mfinit_Nndkchi(long N, long nk, long dk, GEN CHI, long space, long flraw)
{ return (dk == 2)? mf2init_Nkchi(N, nk >> 1, CHI, space, flraw)
                  : mfinit_Nkchi(N, nk, CHI, space, flraw); }
static GEN
mfinit_i(GEN NK, long space)
{
  GEN CHI, mf;
  long N, k, dk, joker;
  if (checkmf_i(NK))
  {
    N = mf_get_N(NK);
    Qtoss(mf_get_gk(NK), &k, &dk);
    CHI = mf_get_CHI(NK);
  }
  else if ((mf = checkMF_i(NK)))
  {
    long s = MF_get_space(mf);
    if (s == space) return mf;
    Qtoss(MF_get_gk(mf), &k, &dk);
    if (dk == 1 && k > 1 && space == mf_NEW && (s == mf_CUSP || s == mf_FULL))
      return mfinittonew(mf);
    N = MF_get_N(mf);
    CHI = MF_get_CHI(mf);
  }
  else
    checkNK2(NK, &N, &k, &dk, &CHI, 1);
  joker = !CHI || typ(CHI) == t_COL;
  if (joker)
  {
    GEN mf, vCHI = CHI;
    long i, j, l;
    if (CHI && lg(CHI) == 1) return cgetg(1,t_VEC);
    if (k < 0) return mfEMPTYall(N, uutoQ(k,dk), CHI, space);
    if (k == 1 && dk == 1 && space != mf_EISEN)
    {
      GEN TMP, vSP, gN, gs;
      pari_timer tt;
      if (space != mf_CUSP && space != mf_NEW)
        pari_err_IMPL("mfinit([N,1,wildcard], space != cusp or new space)");
      if (wt1empty(N)) return mfEMPTYall(N, gen_1, CHI, space);
      vCHI = mf1chars(N,vCHI);
      l = lg(vCHI); mf = cgetg(l, t_VEC); if (l == 1) return mf;
      TMP = mf1_pre(N); vSP = get_vDIH(N, NULL);
      gN = utoipos(N); gs = utoi(space);
      if (DEBUGLEVEL) timer_start(&tt);
      for (i = j = 1; i < l; i++)
      {
        pari_sp av = avma;
        GEN c = gel(vCHI,i), z = mf1init(N, c, TMP, vSP, space, 0);
        if (z) z = gerepilecopy(av, z);
        else
        {
          set_avma(av);
          if (CHI) z = mfEMPTY(mkvec4(gN,gen_1,c,gs));
        }
        if (z) gel(mf, j++) = z;
        if (DEBUGLEVEL)
          timer_printf(&tt, "mf1basis: character %ld / %ld (order = %ld)",
                       i, l-1, mfcharorder(c));
      }
    }
    else
    {
      vCHI = mfchars(N,k,dk,vCHI);
      l = lg(vCHI); mf = cgetg(l, t_VEC);
      for (i = j = 1; i < l; i++)
      {
        pari_sp av = avma;
        GEN v = mfinit_Nndkchi(N, k, dk, gel(vCHI,i), space, 0);
        if (MF_get_dim(v) || CHI) gel(mf, j++) = v; else set_avma(av);
      }
    }
    setlg(mf,j);
    if (!CHI) gen_sort_inplace(mf, NULL, &cmp_ord, NULL);
    return mf;
  }
  return mfinit_Nndkchi(N, k, dk, CHI, space, 0);
}
GEN
mfinit(GEN NK, long space)
{
  pari_sp av = avma;
  return gerepilecopy(av, mfinit_i(NK, space));
}

/* UTILITY FUNCTIONS */
static void
cusp_canon(GEN cusp, long N, long *pA, long *pC)
{
  pari_sp av = avma;
  long A, C, tc, cg;
  if (N <= 0) pari_err_DOMAIN("mfcuspwidth","N","<=",gen_0,stoi(N));
  if (!cusp || (tc = typ(cusp)) == t_INFINITY) { *pA = 1; *pC = N; return; }
  if (tc != t_INT && tc != t_FRAC) pari_err_TYPE("checkcusp", cusp);
  Qtoss(cusp, &A,&C);
  if (N % C)
  {
    ulong uC;
    long u = Fl_invgen((C-1)%N + 1, N, &uC);
    A = Fl_mul(A, u, N);
    C = (long)uC;
  }
  cg = ugcd(C, N/C);
  while (ugcd(A, N) > 1) A += cg;
  *pA = A % N; *pC = C; set_avma(av);
}
static long
mfcuspcanon_width(long N, long C)
{ return (!C || C == N)? 1 : N / ugcd(N, Fl_sqr(umodsu(C,N),N)); }
/* v = [a,c] a ZC, width of cusp (a:c) */
static long
mfZC_width(long N, GEN v)
{
  ulong C = umodiu(gel(v,2), N);
  return (C == 0)? 1: N / ugcd(N, Fl_sqr(C,N));
}
long
mfcuspwidth(GEN gN, GEN cusp)
{
  long N = 0, A, C;
  GEN mf;
  if (typ(gN) == t_INT) N = itos(gN);
  else if ((mf = checkMF_i(gN))) N = MF_get_N(mf);
  else pari_err_TYPE("mfcuspwidth", gN);
  cusp_canon(cusp, N, &A, &C);
  return mfcuspcanon_width(N, C);
}

/* Q a t_INT */
static GEN
findq(GEN al, GEN Q)
{
  long n;
  if (typ(al) == t_FRAC && cmpii(gel(al,2), Q) <= 0)
    return mkvec(mkvec2(gel(al,1), gel(al,2)));
  n = 1 + (long)ceil(2.0781*gtodouble(glog(Q, LOWDEFAULTPREC)));
  return contfracpnqn(gboundcf(al,n), n);
}
static GEN
findqga(long N, GEN z)
{
  GEN Q, LDC, CK = NULL, DK = NULL, ma, x, y = imag_i(z);
  long j, l;
  if (gcmpgs(gmulsg(2*N, y), 1) >= 0) return NULL;
  x = real_i(z);
  Q = ground(ginv(gsqrt(gmulsg(N, y), LOWDEFAULTPREC)));
  LDC = findq(gmulsg(-N,x), Q);
  ma = gen_1; l = lg(LDC);
  for (j = 1; j < l; j++)
  {
    GEN D, DC = gel(LDC,j), C1 = gel(DC,2);
    if (cmpii(C1,Q) > 0) break;
    D = gel(DC,1);
    if (ugcdiu(D,N) == 1)
    {
      GEN C = mului(N, C1), den;
      den = gadd(gsqr(gmul(C,y)), gsqr(gadd(D, gmul(C,x))));
      if (gcmp(den, ma) < 0) { ma = den; CK = C; DK = D; }
    }
  }
  return DK? mkvec2(CK, DK): NULL;
}

static long
valNC2(GEN P, GEN E, long e)
{
  long i, d = 1, l = lg(P);
  for (i = 1; i < l; i++)
  {
    long v = u_lval(e, P[i]) << 1;
    if (v == E[i] + 1) v--;
    d *= upowuu(P[i], v);
  }
  return d;
}

static GEN
findqganew(long N, GEN z)
{
  GEN MI, DI, x = real_i(z), y = imag_i(z), Ck = gen_0, Dk = gen_1, fa, P, E;
  long i;
  MI = uutoQ(1,N);
  DI = mydivisorsu(mysqrtu(N));
  fa = myfactoru(N); P = gel(fa,1); E = gel(fa,2);
  for (i = 1; i < lg(DI); i++)
  {
    long e = DI[i], g;
    GEN U, C, D, m;
    (void)cxredsl2(gmulsg(e, z), &U);
    C = gcoeff(U,2,1); if (!signe(C)) continue;
    D = gcoeff(U,2,2);
    g = ugcdiu(D,e);
    if (g > 1) { C = muliu(C,e/g); D = diviuexact(D,g); } else C = muliu(C,e);
    m = gadd(gsqr(gadd(gmul(C, x), D)), gsqr(gmul(C, y)));
    m = gdivgu(m, valNC2(P, E, e));
    if (gcmp(m, MI) < 0) { MI = m; Ck = C; Dk = D; }
  }
  return signe(Ck)? mkvec2(Ck, Dk): NULL;
}

/* Return z' and U = [a,b;c,d] \in SL_2(Z), z' = U*z,
 * Im(z')/width(U.oo) > sqrt(3)/(2N). Set *pczd = c*z+d */
static GEN
cxredga0N(long N, GEN z, GEN *pU, GEN *pczd, long flag)
{
  GEN v = NULL, A, B, C, D;
  long e;
  if (N == 1) return cxredsl2_i(z, pU, pczd);
  e = gexpo(gel(z,2));
  if (e < 0) z = gprec_wensure(z, precision(z) + nbits2extraprec(-e));
  v = flag? findqganew(N,z): findqga(N,z);
  if (!v) { *pU = matid(2); *pczd = gen_1; return z; }
  C = gel(v,1);
  D = gel(v,2);
  if (!is_pm1(bezout(C,D, &B,&A))) pari_err_BUG("cxredga0N [gcd > 1]");
  B = negi(B);
  *pU = mkmat2(mkcol2(A,C), mkcol2(B,D));
  *pczd = gadd(gmul(C,z), D);
  return gdiv(gadd(gmul(A,z), B), *pczd);
}

static GEN
lfunthetaall(GEN b, GEN vL, GEN t, long bitprec)
{
  long i, l = lg(vL);
  GEN v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN T, L = gel(vL,i), a0 = gel(L,1), ldata = gel(L,2);
    GEN van = gel(ldata_get_an(ldata),2);
    if (lg(van) == 1)
    {
      T = gmul(b, a0);
      if (isexactzero(T)) { GEN z = real_0_bit(-bitprec); T = mkcomplex(z,z); }
    }
    else
    {
      T = gmul2n(lfuntheta(ldata, t, 0, bitprec), -1);
      T = gmul(b, gadd(a0, T));
    }
    gel(v,i) = T;
  }
  return l == 2? gel(v,1): v;
}

/* P in ZX, irreducible */
static GEN
ZX_roots(GEN P, long prec)
{
  long d = degpol(P);
  if (d == 1) return mkvec(gen_0);
  if (d == 2 && isint1(gel(P,2)) && isintzero(gel(P,3)) && isint1(gel(P,4)))
    return mkvec2(powIs(3), gen_I()); /* order as polroots */
  return (ZX_sturm_irred(P) == d)? ZX_realroots_irred(P, prec)
                                 : QX_complex_roots(P, prec);
}
/* initializations for RgX_RgV_eval / RgC_embed */
static GEN
rootspowers(GEN v)
{
  long i, l = lg(v);
  GEN w = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(w,i) = gpowers(gel(v,i), l-2);
  return w;
}
/* mf embeddings attached to Q(chi)/(T), chi attached to cyclotomic P */
static GEN
getembed(GEN P, GEN T, GEN zcyclo, long prec)
{
  long i, l;
  GEN v;
  if (degpol(P) == 1) P = NULL; /* mfcharpol for quadratic char */
  if (degpol(T) == 1) T = NULL; /* dim 1 orbit */
  if (T && P)
  { /* K(y) / (T(y)), K = Q(t)/(P) cyclotomic */
    GEN vr = RgX_is_ZX(T)? ZX_roots(T,prec): roots(RgX_embed1(T,zcyclo), prec);
    v = rootspowers(vr); l = lg(v);
    for (i = 1; i < l; i++) gel(v,i) = mkcol3(P,zcyclo,gel(v,i));
  }
  else if (T)
  { /* Q(y) / (T(y)), T noncyclotomic */
    GEN vr = ZX_roots(T, prec);
    v = rootspowers(vr); l = lg(v);
    for (i = 1; i < l; i++) gel(v,i) = mkcol2(T, gel(v,i));
  }
  else /* cyclotomic or rational */
    v = mkvec(P? mkvec2(P, zcyclo): cgetg(1,t_VEC));
  return v;
}
static GEN
grootsof1_CHI(GEN CHI, long prec)
{ return grootsof1(mfcharorder(CHI), prec); }
/* return the [Q(F):Q(chi)] embeddings of F */
static GEN
mfgetembed(GEN F, long prec)
{
  GEN T = mf_get_field(F), CHI = mf_get_CHI(F), P = mfcharpol(CHI);
  return getembed(P, T, grootsof1_CHI(CHI, prec), prec);
}
static GEN
mfchiembed(GEN mf, long prec)
{
  GEN CHI = MF_get_CHI(mf), P = mfcharpol(CHI);
  return getembed(P, pol_x(0), grootsof1_CHI(CHI, prec), prec);
}
/* mfgetembed for the successive eigenforms in MF_get_newforms */
static GEN
mfeigenembed(GEN mf, long prec)
{
  GEN vP = MF_get_fields(mf), vF = MF_get_newforms(mf);
  GEN zcyclo, vE, CHI = MF_get_CHI(mf), P = mfcharpol(CHI);
  long i, l = lg(vP);
  vF = Q_remove_denom(liftpol_shallow(vF), NULL);
  prec += nbits2extraprec(gexpo(vF));
  zcyclo = grootsof1_CHI(CHI, prec);
  vE = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(vE,i) = getembed(P, gel(vP,i), zcyclo, prec);
  return vE;
}

static int
checkPv(GEN P, GEN v)
{ return typ(P) == t_POL && is_vec_t(typ(v)) && lg(v)-1 >= degpol(P); }
static int
checkemb_i(GEN E)
{
  long t = typ(E), l = lg(E);
  if (t == t_VEC) return l == 1 || (l == 3 && checkPv(gel(E,1), gel(E,2)));
  if (t != t_COL) return 0;
  if (l == 3) return checkPv(gel(E,1), gel(E,2));
  return l == 4 && is_vec_t(typ(gel(E,2))) && checkPv(gel(E,1), gel(E,3));
}
static GEN
anyembed(GEN v, GEN E)
{
  switch(typ(v))
  {
    case t_VEC: case t_COL: return mfvecembed(E, v);
    case t_MAT: return mfmatembed(E, v);
  }
  return mfembed(E, v);
}
GEN
mfembed0(GEN E, GEN v, long prec)
{
  pari_sp av = avma;
  GEN mf, vE = NULL;
  if (checkmf_i(E)) vE = mfgetembed(E, prec);
  else if ((mf = checkMF_i(E))) vE = mfchiembed(mf, prec);
  if (vE)
  {
    long i, l = lg(vE);
    GEN w;
    if (!v) return gerepilecopy(av, l == 2? gel(vE,1): vE);
    w = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(w,i) = anyembed(v, gel(vE,i));
    return gerepilecopy(av, l == 2? gel(w,1): w);
  }
  if (!checkemb_i(E) || !v) pari_err_TYPE("mfembed", E);
  return gerepilecopy(av, anyembed(v,E));
}

/* dummy lfun create for theta evaluation */
static GEN
mfthetaancreate(GEN van, GEN N, GEN k)
{
  GEN L = zerovec(6);
  gel(L,1) = lfuntag(t_LFUN_GENERIC, van);
  gel(L,3) = mkvec2(gen_0, gen_1);
  gel(L,4) = k;
  gel(L,5) = N; return L;
}
/* destroy van and prepare to evaluate theta(sigma(van)), for all sigma in
 * embeddings vector vE */
static GEN
van_embedall(GEN van, GEN vE, GEN gN, GEN gk)
{
  GEN a0 = gel(van,1), vL;
  long i, lE = lg(vE), l = lg(van);
  van++; van[0] = evaltyp(t_VEC) | evallg(l-1); /* remove a0 */
  vL = cgetg(lE, t_VEC);
  for (i = 1; i < lE; i++)
  {
    GEN E = gel(vE,i), v = mfvecembed(E, van);
    gel(vL,i) = mkvec2(mfembed(E,a0), mfthetaancreate(v, gN, gk));
  }
  return vL;
}

static int
cusp_AC(GEN cusp, long *A, long *C)
{
  switch(typ(cusp))
  {
    case t_INFINITY: *A = 1; *C = 0; break;
    case t_INT:  *A = itos(cusp); *C = 1; break;
    case t_FRAC: *A = itos(gel(cusp, 1)); *C = itos(gel(cusp, 2)); break;
    case t_REAL: case t_COMPLEX:
      *A = 0; *C = 0;
      if (gsigne(imag_i(cusp)) <= 0)
        pari_err_DOMAIN("mfeval","imag(tau)","<=",gen_0,cusp);
      return 0;
    default: pari_err_TYPE("cusp_AC", cusp);
  }
  return 1;
}
static GEN
cusp2mat(long A, long C)
{ long B, D;
  cbezout(A, C, &D, &B);
  return mkmat22s(A, -B, C, D);
}
static GEN
mkS(void) { return mkmat22s(0,-1,1,0); }

/* if t is a cusp, return F(t), else NULL */
static GEN
evalcusp(GEN mf, GEN F, GEN t, long prec)
{
  long A, C;
  GEN R;
  if (!cusp_AC(t, &A,&C)) return NULL;
  if (C % mf_get_N(F) == 0) return gel(mfcoefs_i(F, 0, 1), 1);
  R = mfgaexpansion(mf, F, cusp2mat(A,C), 0, prec);
  return gequal0(gel(R,1))? gmael(R,3,1): gen_0;
}
/* Evaluate an mf closure numerically, i.e., in the usual sense, either for a
 * single tau or a vector of tau; for each, return a vector of results
 * corresponding to all complex embeddings of F. If flag is nonzero, allow
 * replacing F by F | gamma to increase imag(gamma^(-1).tau) [ expensive if
 * MF_EISENSPACE not present ] */
static GEN
mfeval_i(GEN mf, GEN F, GEN vtau, long flag, long bitprec)
{
  GEN L0, vL, vb, sqN, vczd, vTAU, vs, van, vE;
  long N = MF_get_N(mf), N0, ta, lv, i, prec = nbits2prec(bitprec);
  GEN gN = utoipos(N), gk = mf_get_gk(F), gk1 = gsubgs(gk,1), vgk;
  long flscal = 0;

  /* gen_0 is ignored, second component assumes Ramanujan-Petersson in
   * 1/2-integer weight */
  vgk = mkvec2(gen_0, mfiscuspidal(mf,F)? gmul2n(gk1,-1): gk1);
  ta = typ(vtau);
  if (!is_vec_t(ta)) { flscal = 1; vtau = mkvec(vtau); ta = t_VEC; }
  lv = lg(vtau);
  sqN = sqrtr_abs(utor(N, prec));
  vs = const_vec(lv-1, NULL);
  vb = const_vec(lv-1, NULL);
  vL = cgetg(lv, t_VEC);
  vTAU = cgetg(lv, t_VEC);
  vczd = cgetg(lv, t_VEC);
  L0 = mfthetaancreate(NULL, gN, vgk); /* only for thetacost */
  vE = mfgetembed(F, prec);
  N0 = 0;
  for (i = 1; i < lv; i++)
  {
    GEN z = gel(vtau,i), tau, U;
    long w, n;

    gel(vs,i) = evalcusp(mf, F, z, prec);
    if (gel(vs,i)) continue;
    tau = cxredga0N(N, z, &U, &gel(vczd,i), flag);
    if (!flag) w = 0; else { w = mfZC_width(N, gel(U,1)); tau = gdivgu(tau,w); }
    gel(vTAU,i) = mulcxmI(gmul(tau, sqN));
    n = lfunthetacost(L0, real_i(gel(vTAU,i)), 0, bitprec);
    if (N0 < n) N0 = n;
    if (flag)
    {
      GEN A, al, v = mfslashexpansion(mf, F, ZM_inv(U,NULL), n, 0, &A, prec);
      gel(vL,i) = van_embedall(v, vE, gN, vgk);
      al = gel(A,1);
      if (!gequal0(al))
        gel(vb,i) = gexp(gmul(gmul(gmulsg(w,al),PiI2(prec)), tau), prec);
    }
  }
  if (!flag)
  {
    van = mfcoefs_i(F, N0, 1);
    vL = const_vec(lv-1, van_embedall(van, vE, gN, vgk));
  }
  for (i = 1; i < lv; i++)
  {
    GEN T;
    if (gel(vs,i)) continue;
    T = gpow(gel(vczd,i), gneg(gk), prec);
    if (flag && gel(vb,i)) T = gmul(T, gel(vb,i));
    gel(vs,i) = lfunthetaall(T, gel(vL,i), gel(vTAU,i), bitprec);
  }
  return flscal? gel(vs,1): vs;
}

static long
mfistrivial(GEN F)
{
  switch(mf_get_type(F))
  {
    case t_MF_CONST: return lg(gel(F,2)) == 1;
    case t_MF_LINEAR: case t_MF_LINEAR_BHN: return gequal0(gel(F,3));
    default: return 0;
  }
}

static long
mf_same_k(GEN mf, GEN f) { return gequal(MF_get_gk(mf), mf_get_gk(f)); }
static long
mf_same_CHI(GEN mf, GEN f)
{
  GEN F1, F2, chi1, chi2, CHI1 = MF_get_CHI(mf), CHI2 = mf_get_CHI(f);
  /* are the primitive chars attached to CHI1 and CHI2 equal ? */
  F1 = znconreyconductor(gel(CHI1,1), gel(CHI1,2), &chi1);
  if (typ(F1) == t_VEC) F1 = gel(F1,1);
  F2 = znconreyconductor(gel(CHI2,1), gel(CHI2,2), &chi2);
  if (typ(F2) == t_VEC) F2 = gel(F2,1);
  return equalii(F1,F2) && ZV_equal(chi1,chi2);
}
/* check k and CHI rigorously, but not coefficients nor N */
static long
mfisinspace_i(GEN mf, GEN F)
{
  return mfistrivial(F) || (mf_same_k(mf,F) && mf_same_CHI(mf,F));
}
static void
err_space(GEN F)
{ pari_err_DOMAIN("mftobasis", "form", "does not belong to",
                  strtoGENstr("space"), F); }

static long
mfcheapeisen(GEN mf)
{
  long k, L, N = MF_get_N(mf);
  GEN P;
  if (N <= 70) return 1;
  k = itos(gceil(MF_get_gk(mf)));
  if (odd(k)) k--;
  switch (k)
  {
    case 2:  L = 190; break;
    case 4:  L = 162; break;
    case 6:
    case 8:  L = 88; break;
    case 10: L = 78; break;
    default: L = 66; break;
  }
  P = gel(myfactoru(N), 1);
  return P[lg(P)-1] <= L;
}

static GEN
myimag_i(GEN tau)
{
  long tc = typ(tau);
  if (tc == t_INFINITY || tc == t_INT || tc == t_FRAC)
    return gen_1;
  if (tc == t_VEC)
  {
    long ltau, i;
    GEN z = cgetg_copy(tau, &ltau);
    for (i=1; i<ltau; i++) gel(z,i) = myimag_i(gel(tau,i));
    return z;
  }
  return imag_i(tau);
}

static GEN
mintau(GEN vtau)
{
  if (!is_vec_t(typ(vtau))) return myimag_i(vtau);
  return (lg(vtau) == 1)? gen_1: vecmin(myimag_i(vtau));
}

/* initialization for mfgaexpansion: what does not depend on cusp */
static GEN
mf_eisendec(GEN mf, GEN F, long prec)
{
  GEN B = liftpol_shallow(mfeisensteindec(mf, F)), v = variables_vecsmall(B);
  GEN Mvecj = obj_check(mf, MF_EISENSPACE);
  long l = lg(v), i, ord;
  if (lg(Mvecj) < 5) Mvecj = gel(Mvecj,1);
  ord = itou(gel(Mvecj,4));
  for (i = 1; i < l; i++)
    if (v[i] != 1)
    {
      GEN d;
      long e;
      B = Q_remove_denom(B, &d);
      e = gexpo(B);
      if (e > 0) prec += nbits2prec(e);
      B = gsubst(B, v[i], rootsof1u_cx(ord, prec));
      if (d) B = gdiv(B, d);
      break;
    }
  return B;
}

GEN
mfeval(GEN mf0, GEN F, GEN vtau, long bitprec)
{
  pari_sp av = avma;
  long flnew = 1;
  GEN mf = checkMF_i(mf0);
  if (!mf) pari_err_TYPE("mfeval", mf0);
  if (!checkmf_i(F)) pari_err_TYPE("mfeval", F);
  if (!mfisinspace_i(mf, F)) err_space(F);
  if (!obj_check(mf, MF_EISENSPACE)) flnew = mfcheapeisen(mf);
  if (flnew && gcmpgs(gmulsg(2*MF_get_N(mf), mintau(vtau)), 1) >= 0) flnew = 0;
  return gerepilecopy(av, mfeval_i(mf, F, vtau, flnew, bitprec));
}

static long
val(GEN v, long bit)
{
  long c, l = lg(v);
  for (c = 1; c < l; c++)
    if (gexpo(gel(v,c)) > -bit) return c-1;
  return -1;
}
GEN
mfcuspval(GEN mf, GEN F, GEN cusp, long bitprec)
{
  pari_sp av = avma;
  long lvE, w, N, sb, n, A, C, prec = nbits2prec(bitprec);
  GEN ga, gk, vE;
  mf = checkMF(mf);
  if (!checkmf_i(F)) pari_err_TYPE("mfcuspval",F);
  N = MF_get_N(mf);
  cusp_canon(cusp, N, &A, &C);
  gk = mf_get_gk(F);
  if (typ(gk) != t_INT)
  {
    GEN FT = mfmultheta(F), mf2 = obj_checkbuild(mf, MF_MF2INIT, &mf2init);
    GEN r = mfcuspval(mf2, FT, cusp, bitprec);
    if ((C & 3L) == 2)
    {
      GEN z = uutoQ(1,4);
      r = gsub(r, typ(r) == t_VEC? const_vec(lg(r)-1, z): z);
    }
    return gerepileupto(av, r);
  }
  vE = mfgetembed(F, prec);
  lvE = lg(vE);
  w = mfcuspcanon_width(N, C);
  sb = w * mfsturmNk(N, itos(gk));
  ga = cusp2mat(A,C);
  for (n = 8;; n = minss(sb, n << 1))
  {
    GEN R = mfgaexpansion(mf, F, ga, n, prec), res = liftpol_shallow(gel(R,3));
    GEN v = cgetg(lvE-1, t_VECSMALL);
    long j, ok = 1;
    res = RgC_embedall(res, vE);
    for (j = 1; j < lvE; j++)
    {
      v[j] = val(gel(res,j), bitprec/2);
      if (v[j] < 0) ok = 0;
    }
    if (ok)
    {
      res = cgetg(lvE, t_VEC);
      for (j = 1; j < lvE; j++) gel(res,j) = gadd(gel(R,1), uutoQ(v[j], w));
      return gerepilecopy(av, lvE==2? gel(res,1): res);
    }
    if (n == sb) return lvE==2? mkoo(): const_vec(lvE-1, mkoo()); /* 0 */
  }
}

long
mfiscuspidal(GEN mf, GEN F)
{
  pari_sp av = avma;
  GEN mf2;
  if (space_is_cusp(MF_get_space(mf))) return 1;
  if (typ(mf_get_gk(F)) == t_INT)
  {
    GEN v = mftobasis(mf,F,0), vE = vecslice(v, 1, lg(MF_get_E(mf))-1);
    return gc_long(av, gequal0(vE));
  }
  if (!gequal0(mfak_i(F, 0))) return 0;
  mf2 = obj_checkbuild(mf, MF_MF2INIT, &mf2init);
  return mfiscuspidal(mf2, mfmultheta(F));
}

/* F = vector of newforms in mftobasis format */
static GEN
mffrickeeigen_i(GEN mf, GEN F, GEN vE, long prec)
{
  GEN M, Z, L0, gN = MF_get_gN(mf), gk = MF_get_gk(mf);
  long N0, i, lM, bit = prec2nbits(prec), k = itou(gk);
  long LIM = 5; /* Sturm bound is enough */

  L0 = mfthetaancreate(NULL, gN, gk); /* only for thetacost */
START:
  N0 = lfunthetacost(L0, gen_1, LIM, bit);
  M = mfcoefs_mf(mf, N0, 1);
  lM = lg(F);
  Z = cgetg(lM, t_VEC);
  for (i = 1; i < lM; i++)
  { /* expansion of D * F[i] */
    GEN D, z, van = RgM_RgC_mul(M, Q_remove_denom(gel(F,i), &D));
    GEN L = van_embedall(van, gel(vE,i), gN, gk);
    long l = lg(L), j, bit_add = D? expi(D): 0;
    gel(Z,i) = z = cgetg(l, t_VEC);
    for (j = 1; j < l; j++)
    {
      GEN v, C, C0;
      long m, e;
      for (m = 0; m <= LIM; m++)
      {
        v = lfuntheta(gmael(L,j,2), gen_1, m, bit);
        if (gexpo(v) > bit_add - bit/2) break;
      }
      if (m > LIM) { LIM <<= 1; goto START; }
      C = mulcxpowIs(gdiv(v,conj_i(v)), 2*m - k);
      C0 = grndtoi(C, &e); if (e < 5-bit_accuracy(precision(C))) C = C0;
      gel(z,j) = C;
    }
  }
  return Z;
}
static GEN
mffrickeeigen(GEN mf, GEN vE, long prec)
{
  GEN D = obj_check(mf, MF_FRICKE);
  if (D) { long p = gprecision(D); if (!p || p >= prec) return D; }
  D = mffrickeeigen_i(mf, MF_get_newforms(mf), vE, prec);
  return obj_insert(mf, MF_FRICKE, D);
}

/* integral weight, new space for primitive quadratic character CHIP;
 * MF = vector of embedded eigenforms coefs on mfbasis, by orbit.
 * Assume N > Q > 1 and (Q,f(CHIP)) = 1 */
static GEN
mfatkineigenquad(GEN mf, GEN CHIP, long Q, GEN MF, long bitprec)
{
  GEN L0, la2, S, F, vP, tau, wtau, Z, va, vb, den, coe, sqrtQ, sqrtN;
  GEN M, gN, gk = MF_get_gk(mf);
  long N0, x, yq, i, j, lF, dim, muQ, prec = nbits2prec(bitprec);
  long N = MF_get_N(mf), k = itos(gk), NQ = N / Q;

  /* Q coprime to FC */
  F = MF_get_newforms(mf);
  vP = MF_get_fields(mf);
  lF = lg(F);
  Z = cgetg(lF, t_VEC);
  S = MF_get_S(mf); dim = lg(S) - 1;
  muQ = mymoebiusu(Q);
  if (muQ)
  {
    GEN SQ = cgetg(dim+1,t_VEC), Qk = gpow(stoi(Q), sstoQ(k-2, 2), prec);
    long i, bit2 = bitprec >> 1;
    for (j = 1; j <= dim; j++) gel(SQ,j) = mfak_i(gel(S,j), Q);
    for (i = 1; i < lF; i++)
    {
      GEN S = RgV_dotproduct(gel(F,i), SQ), T = gel(vP,i);
      long e;
      if (degpol(T) > 1 && typ(S) != t_POLMOD) S = gmodulo(S, T);
      S = grndtoi(gdiv(conjvec(S, prec), Qk), &e);
      if (e > -bit2) pari_err_PREC("mfatkineigenquad");
      if (muQ == -1) S = gneg(S);
      gel(Z,i) = S;
    }
    return Z;
  }
  la2 = mfchareval(CHIP, Q); /* 1 or -1 */
  (void)cbezout(Q, NQ, &x, &yq);
  sqrtQ = sqrtr_abs(utor(Q,prec));
  tau = mkcomplex(gadd(sstoQ(-1, NQ), uutoQ(1, 1000)),
                  divru(sqrtQ, N));
  den = gaddgs(gmulsg(NQ, tau), 1);
  wtau = gdiv(gsub(gmulsg(x, tau), sstoQ(yq, Q)), den);
  coe = gpowgs(gmul(sqrtQ, den), k);

  sqrtN = sqrtr_abs(utor(N,prec));
  tau  = mulcxmI(gmul(tau,  sqrtN));
  wtau = mulcxmI(gmul(wtau, sqrtN));
  gN = utoipos(N);
  L0 = mfthetaancreate(NULL, gN, gk); /* only for thetacost */
  N0 = maxss(lfunthetacost(L0,real_i(tau), 0,bitprec),
             lfunthetacost(L0,real_i(wtau),0,bitprec));
  M = mfcoefs_mf(mf, N0, 1);
  va = cgetg(dim+1, t_VEC);
  vb = cgetg(dim+1, t_VEC);
  for (j = 1; j <= dim; j++)
  {
    GEN L, v = vecslice(gel(M,j), 2, N0+1); /* remove a0 */
    settyp(v, t_VEC); L = mfthetaancreate(v, gN, gk);
    gel(va,j) = lfuntheta(L, tau,0,bitprec);
    gel(vb,j) = lfuntheta(L,wtau,0,bitprec);
  }
  for (i = 1; i < lF; i++)
  {
    GEN z, FE = gel(MF,i);
    long l = lg(FE);
    z = cgetg(l, t_VEC);
    for (j = 1; j < l; j++)
    {
      GEN f = gel(FE,j), a = RgV_dotproduct(va,f), b = RgV_dotproduct(vb,f);
      GEN la = ground( gdiv(b, gmul(a,coe)) );
      if (!gequal(gsqr(la), la2)) pari_err_PREC("mfatkineigenquad");
      if (typ(la) == t_INT)
      {
        if (j != 1) pari_err_BUG("mfatkineigenquad");
        z = const_vec(l-1, la); break;
      }
      gel(z,j) = la;
    }
    gel(Z,i) = z;
  }
  return Z;
}

static GEN
myusqrt(ulong a, long prec)
{
  if (a == 1UL) return gen_1;
  if (uissquareall(a, &a)) return utoipos(a);
  return sqrtr_abs(utor(a, prec));
}
/* Assume mf is a nontrivial new space, rational primitive character CHIP
 * and (Q,FC) = 1 */
static GEN
mfatkinmatnewquad(GEN mf, GEN CHIP, long Q, long flag, long PREC)
{
  GEN cM, M, D, MF, den, vE, F = MF_get_newforms(mf);
  long i, c, e, prec, bitprec, lF = lg(F), N = MF_get_N(mf), k = MF_get_k(mf);

  if (Q == 1) return mkvec4(gen_0, matid(MF_get_dim(mf)), gen_1, mf);
  den = gel(MF_get_Minv(mf), 2);
  bitprec = expi(den) + 64;
  if (!flag) bitprec = maxss(bitprec, prec2nbits(PREC));

START:
  prec = nbits2prec(bitprec);
  vE = mfeigenembed(mf, prec);
  M = cgetg(lF, t_VEC);
  for (i = 1; i < lF; i++) gel(M,i) = RgC_embedall(gel(F,i), gel(vE,i));
  if (Q != N)
  {
    D = mfatkineigenquad(mf, CHIP, Q, M, bitprec);
    c = odd(k)? Q: 1;
  }
  else
  {
    D = mffrickeeigen(mf, vE, prec);
    c = mfcharmodulus(CHIP); if (odd(k)) c = -Q/c;
  }
  D = shallowconcat1(D);
  if (vec_isconst(D)) { MF = diagonal_shallow(D); flag = 0; }
  else
  {
    M = shallowconcat1(M);
    MF = RgM_mul(matmuldiagonal(M,D), ginv(M));
  }
  if (!flag) return mkvec4(gen_0, MF, gen_1, mf);

  if (c > 0)
    cM = myusqrt(c, PREC);
  else
  {
    MF = imag_i(MF); c = -c;
    cM = mkcomplex(gen_0, myusqrt(c,PREC));
  }
  if (c != 1) MF = RgM_Rg_mul(MF, myusqrt(c,prec));
  MF = grndtoi(RgM_Rg_mul(MF,den), &e);
  if (e > -32) { bitprec <<= 1; goto START; }
  MF = RgM_Rg_div(MF, den);
  if (is_rational_t(typ(cM)) && !isint1(cM))
  { MF = RgM_Rg_div(MF, cM); cM = gen_1; }
  return mkvec4(gen_0, MF, cM, mf);
}

/* let CHI mod N, Q || N, return \bar{CHI_Q} * CHI_{N/Q} */
static GEN
mfcharAL(GEN CHI, long Q)
{
  GEN G = gel(CHI,1), c = gel(CHI,2), cycc, d, P, E, F;
  long l = lg(c), N = mfcharmodulus(CHI), i;
  if (N == Q) return mfcharconj(CHI);
  if (N == 1) return CHI;
  CHI = leafcopy(CHI);
  gel(CHI,2) = d = leafcopy(c);
  F = znstar_get_faN(G);
  P = gel(F,1);
  E = gel(F,2);
  cycc = znstar_get_conreycyc(G);
  if (!odd(Q) && equaliu(gel(P,1), 2) && E[1] >= 3)
    gel(d,2) = Fp_neg(gel(d,2), gel(cycc,2));
  else for (i = 1; i < l; i++)
    if (!umodui(Q, gel(P,i))) gel(d,i) = Fp_neg(gel(d,i), gel(cycc,i));
  return CHI;
}
static long
atkin_get_NQ(long N, long Q, const char *f)
{
  long NQ = N / Q;
  if (N % Q) pari_err_DOMAIN(f,"N % Q","!=",gen_0,utoi(Q));
  if (ugcd(NQ, Q) > 1) pari_err_DOMAIN(f,"gcd(Q,N/Q)","!=",gen_1,utoi(Q));
  return NQ;
}

/* transform mf to new_NEW if possible */
static GEN
MF_set_new(GEN mf)
{
  GEN vMjd, vj, gk = MF_get_gk(mf);
  long l, j;
  if (MF_get_space(mf) != mf_CUSP
      || typ(gk) != t_INT || itou(gk) == 1) return mf;
  vMjd = MFcusp_get_vMjd(mf); l = lg(vMjd);
  if (l > 1 && gel(vMjd,1)[1] != MF_get_N(mf)) return mf; /* oldspace != 0 */
  mf = shallowcopy(mf);
  gel(mf,1) = shallowcopy(gel(mf,1));
  MF_set_space(mf, mf_NEW);
  vj = cgetg(l, t_VECSMALL);
  for (j = 1; j < l; j++) vj[j] = gel(vMjd, j)[2];
  gel(mf,4) = vj; return mf;
}

/* if flag = 1, rationalize, else don't */
static GEN
mfatkininit_i(GEN mf, long Q, long flag, long prec)
{
  GEN M, B, C, CHI, CHIAL, G, chi, P, z, g, mfB, s, Mindex, Minv;
  long j, l, lim, ord, FC, NQ, cQ, nk, dk, N = MF_get_N(mf);

  B = MF_get_basis(mf); l = lg(B);
  M = cgetg(l, t_MAT); if (l == 1) return mkvec4(gen_0,M,gen_1,mf);
  Qtoss(MF_get_gk(mf), &nk,&dk);
  Q = labs(Q);
  NQ = atkin_get_NQ(N, Q, "mfatkininit");
  CHI = MF_get_CHI(mf);
  CHI = mfchartoprimitive(CHI, &FC);
  ord = mfcharorder(CHI);
  mf = MF_set_new(mf);
  if (MF_get_space(mf) == mf_NEW && ord <= 2 && NQ % FC == 0 && dk == 1)
    return mfatkinmatnewquad(mf, CHI, Q, flag, prec);
  /* now flag != 0 */
  G   = gel(CHI,1);
  chi = gel(CHI,2);
  if (Q == N) { g = mkmat22s(0, -1, N, 0); cQ = NQ; } /* Fricke */
  else
  {
    GEN F, gQP = utoi(ugcd(Q, FC));
    long t, v;
    chi = znchardecompose(G, chi, gQP);
    F = znconreyconductor(G, chi, &chi);
    G = znstar0(F,1);
    (void)cbezout(Q, NQ, &t, &v);
    g = mkmat22s(Q*t, 1, -N*v, Q);
    cQ = -NQ*v;
  }
  C = s = gen_1;
  /* N.B. G,chi are G_Q,chi_Q [primitive] at this point */
  if (lg(chi) != 1) C = ginv( znchargauss(G, chi, gen_1, prec2nbits(prec)) );
  if (dk == 1)
  { if (odd(nk)) s = myusqrt(Q,prec); }
  else
  {
    long r = nk >> 1; /* k-1/2 */
    s = gpow(utoipos(Q), mkfracss(odd(r)? 1: 3, 4), prec);
    if (odd(cQ))
    {
      long t = r + ((cQ-1) >> 1);
      s = mkcomplex(s, odd(t)? gneg(s): s);
    }
  }
  if (!isint1(s)) C = gmul(C, s);
  CHIAL = mfcharAL(CHI, Q);
  if (dk == 2)
  {
    ulong q = odd(Q)? Q << 2: Q, Nq = ulcm(q, mfcharmodulus(CHIAL));
    CHIAL = induceN(Nq, CHIAL);
    CHIAL = mfcharmul(CHIAL, induce(gel(CHIAL,1), utoipos(q)));
  }
  CHIAL = mfchartoprimitive(CHIAL,NULL);
  mfB = gequal(CHIAL,CHI)? mf: mfinit_Nndkchi(N,nk,dk,CHIAL,MF_get_space(mf),0);
  Mindex = MF_get_Mindex(mfB);
  Minv = MF_get_Minv(mfB);
  P = z = NULL;
  if (ord > 2) { P = mfcharpol(CHI); z = rootsof1u_cx(ord, prec); }
  lim = maxss(mfsturm(mfB), mfsturm(mf)) + 1;
  for (j = 1; j < l; j++)
  {
    GEN v = mfslashexpansion(mf, gel(B,j), g, lim, 0, NULL, prec+EXTRAPREC64);
    long junk;
    if (!isint1(C)) v = RgV_Rg_mul(v, C);
    v = bestapprnf(v, P, z, prec);
    v = vecpermute_partial(v, Mindex, &junk);
    v = Minv_RgC_mul(Minv, v); /* cf mftobasis_i */
    gel(M, j) = v;
  }
  if (is_rational_t(typ(C)) && !gequal1(C)) { M = gdiv(M, C); C = gen_1; }
  if (mfB == mf) mfB = gen_0;
  return mkvec4(mfB, M, C, mf);
}
GEN
mfatkininit(GEN mf, long Q, long prec)
{
  pari_sp av = avma;
  mf = checkMF(mf); return gerepilecopy(av, mfatkininit_i(mf, Q, 1, prec));
}
static void
checkmfa(GEN z)
{
  if (typ(z) != t_VEC || lg(z) != 5 || typ(gel(z,2)) != t_MAT
      || !checkMF_i(gel(z,4))
      || (!isintzero(gel(z,1)) && !checkMF_i(gel(z,1))))
    pari_err_TYPE("mfatkin [please apply mfatkininit()]",z);
}

/* Apply atkin Q to closure F */
GEN
mfatkin(GEN mfa, GEN F)
{
  pari_sp av = avma;
  GEN z, mfB, MQ, mf;
  checkmfa(mfa);
  mfB= gel(mfa,1);
  MQ = gel(mfa,2);
  mf = gel(mfa,4);
  if (typ(mfB) == t_INT) mfB = mf;
  z = RgM_RgC_mul(MQ, mftobasis_i(mf,F));
  return gerepileupto(av, mflinear(mfB, z));
}

GEN
mfatkineigenvalues(GEN mf, long Q, long prec)
{
  pari_sp av = avma;
  GEN vF, L, CHI, M, mfatk, C, MQ, vE, mfB;
  long N, NQ, l, i;

  mf = checkMF(mf); N = MF_get_N(mf);
  vF = MF_get_newforms(mf); l = lg(vF);
  /* N.B. k is integral */
  if (l == 1) { set_avma(av); return cgetg(1, t_VEC); }
  L = cgetg(l, t_VEC);
  if (Q == 1)
  {
    GEN vP = MF_get_fields(mf);
    for (i = 1; i < l; i++) gel(L,i) = const_vec(degpol(gel(vP,i)), gen_1);
    return L;
  }
  vE = mfeigenembed(mf,prec);
  if (Q == N) return gerepileupto(av, mffrickeeigen(mf, vE, prec));
  Q = labs(Q);
  NQ = atkin_get_NQ(N, Q, "mfatkineigenvalues"); /* != 1 */
  mfatk = mfatkininit(mf, Q, prec);
  mfB= gel(mfatk,1); if (typ(mfB) != t_VEC) mfB = mf;
  MQ = gel(mfatk,2);
  C  = gel(mfatk,3);
  M = row(mfcoefs_mf(mfB,1,1), 2); /* vec of a_1(b_i) for mfbasis functions */
  for (i = 1; i < l; i++)
  {
    GEN c = RgV_dotproduct(RgM_RgC_mul(MQ,gel(vF,i)), M); /* C * eigen_i */
    gel(L,i) = Rg_embedall_i(c, gel(vE,i));
  }
  if (!gequal1(C)) L = gdiv(L, C);
  CHI = MF_get_CHI(mf);
  if (mfcharorder(CHI) <= 2 && NQ % mfcharconductor(CHI) == 0) L = ground(L);
  return gerepilecopy(av, L);
}

/* expand B_d V, keeping same length */
static GEN
bdexpand(GEN V, long d)
{
  GEN W;
  long N, n;
  if (d == 1) return V;
  N = lg(V)-1; W = zerovec(N);
  for (n = 0; n <= (N-1)/d; n++) gel(W, n*d+1) = gel(V, n+1);
  return W;
}
/* expand B_d V, increasing length up to lim */
static GEN
bdexpandall(GEN V, long d, long lim)
{
  GEN W;
  long N, n;
  if (d == 1) return V;
  N = lg(V)-1; W = zerovec(lim);
  for (n = 0; n <= N-1 && n*d <= lim; n++) gel(W, n*d+1) = gel(V, n+1);
  return W;
}

static void
parse_vecj(GEN T, GEN *E1, GEN *E2)
{
  if (lg(T)==3) { *E1 = gel(T,1); *E2 = gel(T,2); }
  else { *E1 = T; *E2 = NULL; }
}

/* g in M_2(Z) ? */
static int
check_M2Z(GEN g)
{  return typ(g) == t_MAT && lg(g) == 3 && lgcols(g) == 3 && RgM_is_ZM(g); }
/* g in SL_2(Z) ? */
static int
check_SL2Z(GEN g) { return check_M2Z(g) && equali1(ZM_det(g)); }

static GEN
mfcharcxeval(GEN CHI, long n, long prec)
{
  ulong ord, N = mfcharmodulus(CHI);
  GEN ordg;
  if (N == 1) return gen_1;
  if (ugcd(N, labs(n)) > 1) return gen_0;
  ordg = gmfcharorder(CHI);
  ord = itou(ordg);
  return rootsof1q_cx(znchareval_i(CHI,n,ordg), ord, prec);
}

static GEN
RgV_shift(GEN V, GEN gn)
{
  long i, n, l;
  GEN W;
  if (typ(gn) != t_INT) pari_err_BUG("RgV_shift [n not integral]");
  n = itos(gn);
  if (n < 0) pari_err_BUG("RgV_shift [n negative]");
  if (!n) return V;
  W = cgetg_copy(V, &l); if (n > l-1) n = l-1;
  for (i=1; i <= n; i++) gel(W,i) = gen_0;
  for (    ; i < l; i++) gel(W,i) = gel(V, i-n);
  return W;
}
static GEN
hash_eisengacx(hashtable *H, void *E, long w, GEN ga, long n, long prec)
{
  ulong h = H->hash(E);
  hashentry *e = hash_search2(H, E, h);
  GEN v;
  if (e) v = (GEN)e->val;
  else
  {
    v = mfeisensteingacx((GEN)E, w, ga, n, prec);
    hash_insert2(H, E, (void*)v, h);
  }
  return v;
}
static GEN
vecj_expand(GEN B, hashtable *H, long w, GEN ga, long n, long prec)
{
  GEN E1, E2, v;
  parse_vecj(B, &E1, &E2);
  v = hash_eisengacx(H, (void*)E1, w, ga, n, prec);
  if (E2)
  {
    GEN u = hash_eisengacx(H, (void*)E2, w, ga, n, prec);
    GEN a = gadd(gel(v,1), gel(u,1));
    GEN b = RgV_mul_RgXn(gel(v,2), gel(u,2));
    v = mkvec2(a,b);
  }
  return v;
}
static GEN
shift_M(GEN M, GEN Valpha, long w)
{
  long i, l = lg(Valpha);
  GEN almin = vecmin(Valpha);
  for (i = 1; i < l; i++)
  {
    GEN alpha = gel(Valpha, i), gsh = gmulsg(w, gsub(alpha,almin));
    gel(M,i) = RgV_shift(gel(M,i), gsh);
  }
  return almin;
}
static GEN mfeisensteinspaceinit(GEN NK);
#if 0
/* ga in M_2^+(Z)), n >= 0 */
static GEN
mfgaexpansion_init(GEN mf, GEN ga, long n, long prec)
{
  GEN M, Mvecj, vecj, almin, Valpha;
  long i, w, l, N = MF_get_N(mf), c = itos(gcoeff(ga,2,1));
  hashtable *H;

  if (c % N == 0)
  { /* ga in G_0(N), trivial case; w = 1 */
    GEN chid = mfcharcxeval(MF_get_CHI(mf), itos(gcoeff(ga,2,2)), prec);
    return mkvec2(chid, utoi(n));
  }

  Mvecj = obj_checkbuild(mf, MF_EISENSPACE, &mfeisensteinspaceinit);
  if (lg(Mvecj) < 5) pari_err_IMPL("mfgaexpansion_init in this case");
  w = mfcuspcanon_width(N, c);
  vecj = gel(Mvecj, 3);
  l = lg(vecj);
  M = cgetg(l, t_VEC);
  Valpha = cgetg(l, t_VEC);
  H = hash_create(l, (ulong(*)(void*))&hash_GEN,
                     (int(*)(void*,void*))&gidentical, 1);
  for (i = 1; i < l; i++)
  {
    GEN v = vecj_expand(gel(vecj,i), H, w, ga, n, prec);
    gel(Valpha,i) = gel(v,1);
    gel(M,i) = gel(v,2);
  }
  almin = shift_M(M, Valpha, w);
  return mkvec3(almin, utoi(w), M);
}
/* half-integer weight not supported; vF = [F,eisendec(F)].
 * Minit = mfgaexpansion_init(mf, ga, n, prec) */
static GEN
mfgaexpansion_with_init(GEN Minit, GEN vF)
{
  GEN v;
  if (lg(Minit) == 3)
  { /* ga in G_0(N) */
    GEN chid = gel(Minit,1), gn = gel(Minit,2);
    v = mfcoefs_i(gel(vF,1), itou(gn), 1);
    v = mkvec3(gen_0, gen_1, RgV_Rg_mul(v,chid));
  }
  else
  {
    GEN V = RgM_RgC_mul(gel(Minit,3), gel(vF,2));
    v = mkvec3(gel(Minit,1), gel(Minit,2), V);
  }
  return v;
}
#endif

/* B = mfeisensteindec(F) already embedded, ga in M_2^+(Z)), n >= 0 */
static GEN
mfgaexpansion_i(GEN mf, GEN B0, GEN ga, long n, long prec)
{
  GEN M, Mvecj, vecj, almin, Valpha, B, E = NULL;
  long i, j, w, nw, l, N = MF_get_N(mf), bit = prec2nbits(prec) / 2;
  hashtable *H;

  Mvecj = obj_check(mf, MF_EISENSPACE);
  if (lg(Mvecj) < 5) { E = gel(Mvecj, 2); Mvecj = gel(Mvecj, 1); }
  vecj = gel(Mvecj, 3);
  l = lg(vecj);
  B = cgetg(l, t_COL);
  M = cgetg(l, t_VEC);
  Valpha = cgetg(l, t_VEC);
  w = mfZC_width(N, gel(ga,1));
  nw = E ? n + w : n;
  H = hash_create(l, (ulong(*)(void*))&hash_GEN,
                     (int(*)(void*,void*))&gidentical, 1);
  for (i = j = 1; i < l; i++)
  {
    GEN v;
    if (gequal0(gel(B0,i))) continue;
    v = vecj_expand(gel(vecj,i), H, w, ga, nw, prec);
    gel(B,j) = gel(B0,i);
    gel(Valpha,j) = gel(v,1);
    gel(M,j) = gel(v,2); j++;
  }
  setlg(Valpha, j);
  setlg(B, j);
  setlg(M, j); l = j;
  if (l == 1) return mkvec3(gen_0, utoi(w), zerovec(n+1));
  almin = shift_M(M, Valpha, w);
  B = RgM_RgC_mul(M, B); l = lg(B);
  for (i = 1; i < l; i++)
    if (gexpo(gel(B,i)) < -bit) gel(B,i) = gen_0;
  settyp(B, t_VEC);
  if (E)
  {
    GEN v, e;
    long ell = 0, vB, ve;
    for (i = 1; i < l; i++)
      if (!gequal0(gel(B,i))) break;
    vB = i-1;
    v = hash_eisengacx(H, (void*)E, w, ga, n + vB, prec);
    e = gel(v,2); l = lg(e);
    for (i = 1; i < l; i++)
      if (!gequal0(gel(e,i))) break;
    ve = i-1;
    almin = gsub(almin, gel(v,1));
    if (gsigne(almin) < 0)
    {
      GEN gell = gceil(gmulsg(-w, almin));
      ell = itos(gell);
      almin = gadd(almin, gdivgu(gell, w));
      if (nw < ell) pari_err_IMPL("alpha < 0 in mfgaexpansion");
    }
    if (ve) { ell += ve; e = vecslice(e, ve+1, l-1); }
    B = vecslice(B, ell + 1, minss(n + ell + 1, lg(B)-1));
    B = RgV_div_RgXn(B, e);
  }
  return mkvec3(almin, utoi(w), B);
}

/* Theta multiplier: assume 4 | C, (C,D)=1 */
static GEN
mfthetamultiplier(GEN C, GEN D)
{
  long s = kronecker(C, D);
  if (Mod4(D) == 1) return s > 0 ? gen_1: gen_m1;
  return s > 0? powIs(3): gen_I();
}
/* theta | [*,*;C,D] defined over Q(i) [else over Q] */
static int
mfthetaI(long C, long D) { return odd(C) || (D & 3) == 3; }
/* (theta | M) [0..n], assume (C,D) = 1 */
static GEN
mfthetaexpansion(GEN M, long n)
{
  GEN w, s, al, sla, E, V = zerovec(n+1), C = gcoeff(M,2,1), D = gcoeff(M,2,2);
  long lim, la, f, C4 = Mod4(C);
  switch (C4)
  {
    case 0: al = gen_0; w = gen_1;
      s = mfthetamultiplier(C,D);
      lim = usqrt(n); gel(V, 1) = s;
      s = gmul2n(s, 1);
      for (f = 1; f <= lim; f++) gel(V, f*f + 1) = s;
      break;
    case 2: al = uutoQ(1,4); w = gen_1;
      E = subii(C, shifti(D,1)); /* (E, D) = 1 */
      s = gmul2n(mfthetamultiplier(E, D), 1);
      if ((!signe(E) && equalim1(D)) || (signe(E) > 0 && signe(C) < 0))
        s = gneg(s);
      lim = (usqrt(n << 2) - 1) >> 1;
      for (f = 0; f <= lim; f++) gel(V, f*(f+1) + 1) = s;
      break;
    default: al = gen_0; w = utoipos(4);
      la = (-Mod4(D)*C4) & 3L;
      E = negi(addii(D, mului(la, C)));
      s = mfthetamultiplier(E, C); /* (E,C) = 1 */
      if (signe(C) < 0 && signe(E) >= 0) s = gneg(s);
      s = gsub(s, mulcxI(s));
      sla = gmul(s, powIs(-la));
      lim = usqrt(n); gel(V, 1) = gmul2n(s, -1);
      for (f = 1; f <= lim; f++) gel(V, f*f + 1) = odd(f) ? sla : s;
      break;
  }
  return mkvec3(al, w, V);
}

/* F 1/2 integral weight */
static GEN
mf2gaexpansion(GEN mf2, GEN F, GEN ga, long n, long prec)
{
  GEN FT = mfmultheta(F), mf = obj_checkbuild(mf2, MF_MF2INIT, &mf2init);
  GEN res, V1, Tres, V2, al, V, gsh, C = gcoeff(ga,2,1);
  long w2, N = MF_get_N(mf), w = mfcuspcanon_width(N, umodiu(C,N));
  long ext = (Mod4(C) != 2)? 0: (w+3) >> 2;
  long prec2 = prec + nbits2extraprec((long)M_PI/(2*M_LN2)*sqrt(n + ext));
  res = mfgaexpansion(mf, FT, ga, n + ext, prec2);
  Tres = mfthetaexpansion(ga, n + ext);
  V1 = gel(res,3);
  V2 = gel(Tres,3);
  al = gsub(gel(res,1), gel(Tres,1));
  w2 = itos(gel(Tres,2));
  if (w != itos(gel(res,2)) || w % w2)
    pari_err_BUG("mf2gaexpansion [incorrect w2 or w]");
  if (w2 != w) V2 = bdexpand(V2, w/w2);
  V = RgV_div_RgXn(V1, V2);
  gsh = gfloor(gmulsg(w, al));
  if (!gequal0(gsh))
  {
    al = gsub(al, gdivgu(gsh, w));
    if (gsigne(gsh) > 0)
    {
      V = RgV_shift(V, gsh);
      V = vecslice(V, 1, n + 1);
    }
    else
    {
      long sh = -itos(gsh), i;
      if (sh > ext) pari_err_BUG("mf2gaexpansion [incorrect sh]");
      for (i = 1; i <= sh; i++)
        if (!gequal0(gel(V,i))) pari_err_BUG("mf2gaexpansion [sh too large]");
      V = vecslice(V, sh+1, n + sh+1);
    }
  }
  obj_free(mf); return mkvec3(al, stoi(w), gprec_wtrunc(V, prec));
}

static GEN
mfgaexpansionatkin(GEN mf, GEN F, GEN C, GEN D, long Q, long n, long prec)
{
  GEN mfa = mfatkininit_i(mf, Q, 0, prec), MQ = gel(mfa,2);
  long i, FC, k = MF_get_k(mf);
  GEN x, v, V, z, s, CHI = mfchartoprimitive(MF_get_CHI(mf), &FC);

  /* V = mfcoefs(F | w_Q, n), can't use mfatkin because MQ nonrational */
  V = RgM_RgC_mul(mfcoefs_mf(mf,n,1), RgM_RgC_mul(MQ, mftobasis_i(mf,F)));
  (void)bezout(utoipos(Q), C, &x, &v);
  s = mfchareval(CHI, (umodiu(x, FC) * umodiu(D, FC)) % FC);
  s = gdiv(s, gpow(utoipos(Q), uutoQ(k,2), prec));
  V = RgV_Rg_mul(V, s);
  z = rootsof1powinit(umodiu(D,Q)*umodiu(v,Q) % Q, Q, prec);
  for (i = 1; i <= n+1; i++) gel(V,i) = gmul(gel(V,i), rootsof1pow(z, i-1));
  return mkvec3(gen_0, utoipos(Q), V);
}

static long
inveis_extraprec(long N, GEN ga, GEN Mvecj, long n)
{
  long e, w = mfZC_width(N, gel(ga,1));
  GEN f, E = gel(Mvecj,2), v = mfeisensteingacx(E, w, ga, n, DEFAULTPREC);
  v = gel(v,2);
  f = RgV_to_RgX(v,0); n -= RgX_valrem(f, &f);
  e = gexpo(RgXn_inv(f, n+1));
  return (e > 0)? nbits2extraprec(e): 0;
}
/* allow F of the form [F, mf_eisendec(F)]~ */
static GEN
mfgaexpansion(GEN mf, GEN F, GEN ga, long n, long prec)
{
  GEN v, EF = NULL, res, Mvecj, c, d;
  long precnew, N;

  if (n < 0) pari_err_DOMAIN("mfgaexpansion", "n", "<", gen_0, stoi(n));
  if (typ(F) == t_COL && lg(F) == 3) { EF = gel(F,2); F = gel(F,1); }
  if (!checkmf_i(F)) pari_err_TYPE("mfgaexpansion", F);
  if (!check_SL2Z(ga)) pari_err_TYPE("mfgaexpansion",ga);
  if (typ(mf_get_gk(F)) != t_INT) return mf2gaexpansion(mf, F, ga, n, prec);
  c = gcoeff(ga,2,1);
  d = gcoeff(ga,2,2);
  N = MF_get_N(mf);
  if (!umodiu(c, mf_get_N(F)))
  { /* trivial case: ga in Gamma_0(N) */
    long w = mfcuspcanon_width(N, umodiu(c,N));
    GEN CHI = mf_get_CHI(F);
    GEN chid = mfcharcxeval(CHI, umodiu(d,mfcharmodulus(CHI)), prec);
    v = mfcoefs_i(F, n/w, 1); if (!isint1(chid)) v = RgV_Rg_mul(v,chid);
    return mkvec3(gen_0, stoi(w), bdexpandall(v,w,n+1));
  }
  mf = MF_set_new(mf);
  if (MF_get_space(mf) == mf_NEW)
  {
    long cN = umodiu(c,N), g = ugcd(cN,N), Q = N/g;
    GEN CHI = MF_get_CHI(mf);
    if (ugcd(cN, Q)==1 && mfcharorder(CHI) <= 2
                       && g % mfcharconductor(CHI) == 0
                       && degpol(mf_get_field(F)) == 1)
      return mfgaexpansionatkin(mf, F, c, d, Q, n, prec);
  }
  Mvecj = obj_checkbuild(mf, MF_EISENSPACE, &mfeisensteinspaceinit);
  precnew = prec;
  if (lg(Mvecj) < 5) precnew += inveis_extraprec(N, ga, Mvecj, n);
  if (!EF) EF = mf_eisendec(mf, F, precnew);
  res = mfgaexpansion_i(mf, EF, ga, n, precnew);
  return precnew == prec ? res : gprec_wtrunc(res, prec);
}

/* parity = -1 or +1 */
static GEN
findd(long N, long parity)
{
  GEN L, D = mydivisorsu(N);
  long i, j, l = lg(D);
  L = cgetg(l, t_VEC);
  for (i = j = 1; i < l; i++)
  {
    long d = D[i];
    if (parity == -1) d = -d;
    if (sisfundamental(d)) gel(L,j++) = stoi(d);
  }
  setlg(L,j); return L;
}
/* does ND contain a divisor of N ? */
static int
seenD(long N, GEN ND)
{
  long j, l = lg(ND);
  for (j = 1; j < l; j++)
    if (N % ND[j] == 0) return 1;
  return 0;
}
static GEN
search_levels(GEN vN, const char *f)
{
  switch(typ(vN))
  {
    case t_INT: vN = mkvecsmall(itos(vN)); break;
    case t_VEC: case t_COL: vN = ZV_to_zv(vN); break;
    case t_VECSMALL: vN = leafcopy(vN); break;
    default: pari_err_TYPE(f, vN);
  }
  vecsmall_sort(vN); return vN;
}
GEN
mfsearch(GEN NK, GEN V, long space)
{
  pari_sp av = avma;
  GEN F, gk, NbyD, vN;
  long n, nk, dk, parity, nV, i, lvN;

  if (typ(NK) != t_VEC || lg(NK) != 3) pari_err_TYPE("mfsearch", NK);
  gk = gel(NK,2);
  if (typ(gmul2n(gk, 1)) != t_INT) pari_err_TYPE("mfsearch [k]", gk);
  switch(typ(V))
  {
    case t_VEC: V = shallowtrans(V);
    case t_COL: break;
    default: pari_err_TYPE("mfsearch [V]", V);
  }
  vN = search_levels(gel(NK,1), "mfsearch [N]");
  if (gequal0(V)) { set_avma(av); retmkvec(mftrivial()); }
  lvN = lg(vN);

  Qtoss(gk, &nk,&dk);
  parity = (dk == 1 && odd(nk)) ? -1 : 1;
  nV = lg(V)-2;
  F = cgetg(1, t_VEC);
  NbyD = const_vec(vN[lvN-1], cgetg(1,t_VECSMALL));
  for (n = 1; n < lvN; n++)
  {
    long N = vN[n];
    GEN L;
    if (N <= 0 || (dk == 2 && (N & 3))) continue;
    L = findd(N, parity);
    for (i = 1; i < lg(L); i++)
    {
      GEN mf, M, CO, gD = gel(L,i);
      GEN *ND = (GEN*)NbyD + itou(gD); /* points to NbyD[|D|] */

      if (seenD(N, *ND)) continue;
      mf = mfinit_Nndkchi(N, nk, dk, get_mfchar(gD), space, 1);
      M = mfcoefs_mf(mf, nV, 1);
      CO = inverseimage(M, V); if (lg(CO) == 1) continue;

      F = vec_append(F, mflinear(mf,CO));
      *ND = vecsmall_append(*ND, N); /* add to NbyD[|D|] */
    }
  }
  return gerepilecopy(av, F);
}

static GEN
search_from_split(GEN mf, GEN vap, GEN vlp)
{
  pari_sp av = avma;
  long lvlp = lg(vlp), j, jv, l1;
  GEN v, NK, S1, S, M = NULL;

  S1 = gel(split_i(mf, 1, 0), 1); /* rational newforms */
  l1 = lg(S1);
  if (l1 == 1) return gc_NULL(av);
  v = cgetg(l1, t_VEC);
  S = MF_get_S(mf);
  NK = mf_get_NK(gel(S,1));
  if (lvlp > 1) M = rowpermute(mfcoefs_mf(mf, vlp[lvlp-1], 1), vlp);
  for (j = jv = 1; j < l1; j++)
  {
    GEN vF = gel(S1,j);
    long t;
    for (t = lvlp-1; t > 0; t--)
    { /* lhs = vlp[j]-th coefficient of eigenform */
      GEN rhs = gel(vap,t), lhs = RgMrow_RgC_mul(M, vF, t);
      if (!gequal(lhs, rhs)) break;
    }
    if (!t) gel(v,jv++) = mflinear_i(NK,S,vF);
  }
  if (jv == 1) return gc_NULL(av);
  setlg(v,jv); return v;
}
GEN
mfeigensearch(GEN NK, GEN AP)
{
  pari_sp av = avma;
  GEN k, vN, vap, vlp, vres = cgetg(1, t_VEC), D;
  long n, lvN, i, l, even;

  if (!AP) l = 1;
  else
  {
    l = lg(AP);
    if (typ(AP) != t_VEC) pari_err_TYPE("mfeigensearch",AP);
  }
  vap = cgetg(l, t_VEC);
  vlp = cgetg(l, t_VECSMALL);
  if (l > 1)
  {
    GEN perm = indexvecsort(AP, mkvecsmall(1));
    for (i = 1; i < l; i++)
    {
      GEN v = gel(AP,perm[i]), gp, ap;
      if (typ(v) != t_VEC || lg(v) != 3) pari_err_TYPE("mfeigensearch", AP);
      gp = gel(v,1);
      ap = gel(v,2);
      if (typ(gp) != t_INT || (typ(ap) != t_INT && typ(ap) != t_INTMOD))
        pari_err_TYPE("mfeigensearch", AP);
      gel(vap,i) = ap;
      vlp[i] = itos(gp)+1; if (vlp[i] < 0) pari_err_TYPE("mfeigensearch", AP);
    }
  }
  l = lg(NK);
  if (typ(NK) != t_VEC || l != 3) pari_err_TYPE("mfeigensearch",NK);
  k = gel(NK,2);
  vN = search_levels(gel(NK,1), "mfeigensearch [N]");
  lvN = lg(vN);
  vecsmall_sort(vlp);
  even = !mpodd(k);
  for (n = 1; n < lvN; n++)
  {
    pari_sp av2 = avma;
    GEN mf, L;
    long N = vN[n];
    if (even) D = gen_1;
    else
    {
      long r = (N&3L);
      if (r == 1 || r == 2) continue;
      D = stoi( corediscs(-N, NULL) ); /* < 0 */
    }
    mf = mfinit_i(mkvec3(utoipos(N), k, D), mf_NEW);
    L = search_from_split(mf, vap, vlp);
    if (L) vres = shallowconcat(vres, L); else set_avma(av2);
  }
  return gerepilecopy(av, vres);
}

/* tf_{N,k}(n) */
static GEN
mfnewtracecache(long N, long k, long n, cachenew_t *cache)
{
  GEN C = NULL, S;
  long lcache;
  if (!n) return gen_0;
  S = gel(cache->vnew,N);
  lcache = lg(S);
  if (n < lcache) C = gel(S, n);
  if (C) cache->newHIT++;
  else C = mfnewtrace_i(N,k,n,cache);
  cache->newTOTAL++;
  if (n < lcache) gel(S,n) = C;
  return C;
}

static long
mfdim_Nkchi(long N, long k, GEN CHI, long space)
{
  if (k < 0 || badchar(N,k,CHI)) return 0;
  if (k == 0)
    return mfcharistrivial(CHI) && !space_is_cusp(space)? 1: 0;
  switch(space)
  {
    case mf_NEW: return mfnewdim(N,k,CHI);
    case mf_CUSP:return mfcuspdim(N,k,CHI);
    case mf_OLD: return mfolddim(N,k,CHI);
    case mf_FULL:return mffulldim(N,k,CHI);
    case mf_EISEN: return mfeisensteindim(N,k,CHI);
    default: pari_err_FLAG("mfdim");
  }
  return 0;/*LCOV_EXCL_LINE*/
}
static long
mf1dimsum(long N, long space)
{
  switch(space)
  {
    case mf_NEW:  return mf1newdimsum(N);
    case mf_CUSP: return mf1cuspdimsum(N);
    case mf_OLD:  return mf1olddimsum(N);
  }
  pari_err_FLAG("mfdim");
  return 0; /*LCOV_EXCL_LINE*/
}
/* mfdim for k = nk/dk */
static long
mfdim_Nndkchi(long N, long nk, long dk, GEN CHI, long space)
{ return (dk == 2)? mf2dim_Nkchi(N, nk >> 1, CHI, space)
                  : mfdim_Nkchi(N, nk, CHI, space); }
/* FIXME: use direct dim Gamma1(N) formula, don't compute individual spaces */
static long
mfkdimsum(long N, long k, long dk, long space)
{
  GEN w = mfchars(N, k, dk, NULL);
  long i, j, D = 0, l = lg(w);
  for (i = j = 1; i < l; i++)
  {
    GEN CHI = gel(w,i);
    long d = mfdim_Nndkchi(N,k,dk,CHI,space);
    if (d) D += d * myeulerphiu(mfcharorder(CHI));
  }
  return D;
}
static GEN
mf1dims(long N, GEN vCHI, long space)
{
  GEN D = NULL;
  switch(space)
  {
    case mf_NEW: D = mf1newdimall(N, vCHI); break;
    case mf_CUSP:D = mf1cuspdimall(N, vCHI); break;
    case mf_OLD: D = mf1olddimall(N, vCHI); break;
    default: pari_err_FLAG("mfdim");
  }
  return D;
}
static GEN
mfkdims(long N, long k, long dk, GEN vCHI, long space)
{
  GEN D, w = mfchars(N, k, dk, vCHI);
  long i, j, l = lg(w);
  D = cgetg(l, t_VEC);
  for (i = j = 1; i < l; i++)
  {
    GEN CHI = gel(w,i);
    long d = mfdim_Nndkchi(N,k,dk,CHI,space);
    if (vCHI)
      gel(D, j++) = mkvec2s(d, 0);
    else if (d)
      gel(D, j++) = fmt_dim(CHI, d, 0);
  }
  setlg(D,j); return D;
}
GEN
mfdim(GEN NK, long space)
{
  pari_sp av = avma;
  long N, k, dk, joker;
  GEN CHI, mf;
  if ((mf = checkMF_i(NK))) return utoi(MF_get_dim(mf));
  checkNK2(NK, &N, &k, &dk, &CHI, 2);
  if (!CHI) joker = 1;
  else
    switch(typ(CHI))
    {
      case t_INT: joker = 2; break;
      case t_COL: joker = 3; break;
      default: joker = 0; break;
    }
  if (joker)
  {
    long d;
    GEN D;
    if (k < 0) switch(joker)
    {
      case 1: return cgetg(1,t_VEC);
      case 2: return gen_0;
      case 3: return mfdim0all(CHI);
    }
    if (k == 0)
    {
      if (space_is_cusp(space)) switch(joker)
      {
        case 1: return cgetg(1,t_VEC);
        case 2: return gen_0;
        case 3: return mfdim0all(CHI);
      }
      switch(joker)
      {
        long i, l;
        case 1: retmkvec(fmt_dim(mfchartrivial(),0,0));
        case 2: return gen_1;
        case 3: l = lg(CHI); D = cgetg(l,t_VEC);
                for (i = 1; i < l; i++)
                {
                  long t = mfcharistrivial(gel(CHI,i));
                  gel(D,i) = mkvec2(t? gen_1: gen_0, gen_0);
                }
                return D;
      }
    }
    if (dk == 1 && k == 1 && space != mf_EISEN)
    {
      long fix = 0, space0 = space;
      if (space == mf_FULL) space = mf_CUSP; /* remove Eisenstein part */
      if (joker == 2)
      {
        d = mf1dimsum(N, space);
        if (space0 == mf_FULL) d += mfkdimsum(N,k,dk,mf_EISEN);/*add it back*/
        return gc_utoi(av, d);
      }
      /* must initialize explicitly: trivial spaces for E_k/S_k differ */
      if (space0 == mf_FULL)
      {
        if (!CHI) fix = 1; /* must remove 0 spaces */
        CHI = mfchars(N, k, dk, CHI);
      }
      D = mf1dims(N, CHI, space);
      if (space0 == mf_FULL)
      {
        GEN D2 = mfkdims(N, k, dk, CHI, mf_EISEN);
        D = merge_dims(D, D2, fix? CHI: NULL);
      }
    }
    else
    {
      if (joker==2) { d = mfkdimsum(N,k,dk,space); return gc_utoi(av,d); }
      D = mfkdims(N, k, dk, CHI, space);
    }
    if (!CHI) return gerepileupto(av, vecsort(D, mkvecsmall(1)));
    return gerepilecopy(av, D);
  }
  return utoi( mfdim_Nndkchi(N, k, dk, CHI, space) );
}

GEN
mfbasis(GEN NK, long space)
{
  pari_sp av = avma;
  long N, k, dk;
  GEN mf, CHI;
  if ((mf = checkMF_i(NK))) return gconcat(gel(mf,2), gel(mf,3));
  checkNK2(NK, &N, &k, &dk, &CHI, 0);
  if (dk == 2) return gerepilecopy(av, mf2basis(N, k>>1, CHI, NULL, space));
  mf = mfinit_Nkchi(N, k, CHI, space, 1);
  return gerepilecopy(av, MF_get_basis(mf));
}

static GEN
deg1ser_shallow(GEN a1, GEN a0, long v, long e)
{ return RgX_to_ser(deg1pol_shallow(a1, a0, v), e+2); }
/* r / x + O(1) */
static GEN
simple_pole(GEN r)
{
  GEN S = deg1ser_shallow(gen_0, r, 0, 1);
  setvalser(S, -1); return S;
}

/* F form, E embedding; mfa = mfatkininit or root number (eigenform case) */
static GEN
mflfuncreate(GEN mfa, GEN F, GEN E, GEN N, GEN gk)
{
  GEN LF = cgetg(8,t_VEC), polar = cgetg(1,t_COL), eps;
  long k = itou(gk);
  gel(LF,1) = lfuntag(t_LFUN_MFCLOS, mkvec3(F,E,gen_1));
  if (typ(mfa) != t_VEC)
    eps = mfa; /* cuspidal eigenform: root number; no poles */
  else
  { /* mfatkininit */
    GEN a0, b0, vF, vG, G = NULL;
    GEN M = gel(mfa,2), C = gel(mfa,3), mf = gel(mfa,4);
    M = gdiv(mfmatembed(E, M), C);
    vF = mfvecembed(E, mftobasis_i(mf, F));
    vG = RgM_RgC_mul(M, vF);
    if (gequal(vF,vG)) eps = gen_1;
    else if (gequal(vF,gneg(vG))) eps = gen_m1;
    else
    { /* not self-dual */
      eps = NULL;
      G = mfatkin(mfa, F);
      gel(LF,2) = lfuntag(t_LFUN_MFCLOS, mkvec3(G,E,ginv(C)));
      gel(LF,6) = powIs(k);
    }
    /* polar part */
    a0 = mfembed(E, mfcoef(F,0));
    b0 = eps? gmul(eps,a0): gdiv(mfembed(E, mfcoef(G,0)), C);
    if (!gequal0(b0))
    {
      b0 = mulcxpowIs(gmul2n(b0,1), k);
      polar = vec_append(polar, mkvec2(gk, simple_pole(b0)));
    }
    if (!gequal0(a0))
    {
      a0 = gneg(gmul2n(a0,1));
      polar = vec_append(polar, mkvec2(gen_0, simple_pole(a0)));
    }
  }
  if (eps) /* self-dual */
  {
    gel(LF,2) = mfcharorder(mf_get_CHI(F)) <= 2? gen_0: gen_1;
    gel(LF,6) = mulcxpowIs(eps,k);
  }
  gel(LF,3) = mkvec2(gen_0, gen_1);
  gel(LF,4) = gk;
  gel(LF,5) = N;
  if (lg(polar) == 1) setlg(LF,7); else gel(LF,7) = polar;
  return LF;
}
static GEN
mflfuncreateall(long sd, GEN mfa, GEN F, GEN vE, GEN gN, GEN gk)
{
  long i, l = lg(vE);
  GEN L = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(L,i) = mflfuncreate(sd? gel(mfa,i): mfa, F, gel(vE,i), gN, gk);
  return L;
}
GEN
lfunmf(GEN mf, GEN F, long bitprec)
{
  pari_sp av = avma;
  long i, l, prec = nbits2prec(bitprec);
  GEN L, gk, gN;
  mf = checkMF(mf);
  gk = MF_get_gk(mf);
  gN = MF_get_gN(mf);
  if (typ(gk)!=t_INT) pari_err_IMPL("half-integral weight");
  if (F)
  {
    GEN v;
    long s = MF_get_space(mf);
    if (!checkmf_i(F)) pari_err_TYPE("lfunmf", F);
    if (!mfisinspace_i(mf, F)) err_space(F);
    L = NULL;
    if ((s == mf_NEW || s == mf_CUSP || s == mf_FULL)
        && gequal(mfcoefs_i(F,1,1), mkvec2(gen_0,gen_1)))
    { /* check if eigenform */
      GEN vP, vF, b = mftobasis_i(mf, F);
      long lF, d = degpol(mf_get_field(F));
      v = mfsplit(mf, d, 0);
      vF = gel(v,1);
      vP = gel(v,2); lF = lg(vF);
      for (i = 1; i < lF; i++)
        if (degpol(gel(vP,i)) == d && gequal(gel(vF,i), b))
        {
          GEN vE = mfgetembed(F, prec);
          GEN Z = mffrickeeigen_i(mf, mkvec(b), mkvec(vE), prec);
          L = mflfuncreateall(1, gel(Z,1), F, vE, gN, gk);
          break;
        }
    }
    if (!L)
    { /* not an eigenform: costly general case */
      GEN mfa = mfatkininit_i(mf, itou(gN), 1, prec);
      L = mflfuncreateall(0,mfa, F, mfgetembed(F,prec), gN, gk);
    }
    if (lg(L) == 2) L = gel(L,1);
  }
  else
  {
    GEN M = mfeigenbasis(mf), vE = mfeigenembed(mf, prec);
    GEN v = mffrickeeigen(mf, vE, prec);
    l = lg(vE); L = cgetg(l, t_VEC);
    for (i = 1; i < l; i++)
      gel(L,i) = mflfuncreateall(1,gel(v,i), gel(M,i), gel(vE,i), gN, gk);
  }
  return gerepilecopy(av, L);
}

GEN
mffromell(GEN E)
{
  pari_sp av = avma;
  GEN mf, F, z, v, S;
  long N, i, l;

  checkell(E);
  if (ell_get_type(E) != t_ELL_Q) pari_err_TYPE("mfffromell [E not over Q]", E);
  N = itos(ellQ_get_N(E));
  mf = mfinit_i(mkvec2(utoi(N), gen_2), mf_NEW);
  v = split_i(mf, 1, 0);
  S = gel(v,1); l = lg(S); /* rational newforms */
  F = tag(t_MF_ELL, mkNK(N,2,mfchartrivial()), E);
  z = mftobasis_i(mf, F);
  for(i = 1; i < l; i++)
    if (gequal(z, gel(S,i))) break;
  if (i == l) pari_err_BUG("mffromell [E is not modular]");
  return gerepilecopy(av, mkvec3(mf, F, z));
}

/* returns -1 if not, degree otherwise */
long
polishomogeneous(GEN P)
{
  long i, D, l;
  if (typ(P) != t_POL) return 0;
  D = -1; l = lg(P);
  for (i = 2; i < l; i++)
  {
    GEN c = gel(P,i);
    long d;
    if (gequal0(c)) continue;
    d = polishomogeneous(c);
    if (d < 0) return -1;
    if (D < 0) D = d + i-2; else if (D != d + i-2) return -1;
  }
  return D;
}

/* M a pp((Gram q)^(-1)) ZM; P a homogeneous t_POL, is P spherical ? */
static int
RgX_isspherical(GEN M, GEN P)
{
  pari_sp av = avma;
  GEN S, v = variables_vecsmall(P);
  long i, j, l = lg(v);
  if (l > lg(M)) pari_err(e_MISC, "too many variables in mffromqf");
  S = gen_0;
  for (j = 1; j < l; j++)
  {
    GEN Mj = gel(M, j), Pj = deriv(P, v[j]);
    for (i = 1; i <= j; i++)
    {
      GEN c = gel(Mj, i);
      if (!signe(c)) continue;
      if (i != j) c = shifti(c, 1);
      S = gadd(S, gmul(c, deriv(Pj, v[i])));
    }
  }
  return gc_bool(av, gequal0(S));
}

static GEN
c_QFsimple_i(long n, GEN Q, GEN P)
{
  GEN V, v = qfrep0(Q, utoi(n), 1);
  long i, l = lg(v);
  V = cgetg(l+1, t_VEC);
  if (!P || equali1(P))
  {
    gel(V,1) = gen_1;
    for (i = 2; i <= l; i++) gel(V,i) = utoi(v[i-1] << 1);
  }
  else
  {
    gel(V,1) = gcopy(P);
    for (i = 2; i <= l; i++) gel(V,i) = gmulgu(P, v[i-1] << 1);
  }
  return V;
}

/* v a t_VECSMALL of variable numbers, lg(r) >= lg(v), r is a vector of
 * scalars [not involving any variable in v] */
static GEN
gsubstvec_i(GEN e, GEN v, GEN r)
{
  long i, l = lg(v);
  for(i = 1; i < l; i++) e = gsubst(e, v[i], gel(r,i));
  return e;
}
static GEN
c_QF_i(long n, GEN Q, GEN P)
{
  pari_sp av = avma;
  GEN V, v, va;
  long i, l;
  if (!P || typ(P) != t_POL) return gerepileupto(av, c_QFsimple_i(n, Q, P));
  v = gel(minim(Q, utoi(2*n), NULL), 3);
  va = variables_vecsmall(P);
  V = zerovec(n + 1); l = lg(v);
  for (i = 1; i < l; i++)
  {
    pari_sp av = avma;
    GEN X = gel(v,i);
    long c = (itos(qfeval(Q, X)) >> 1) + 1;
    gel(V, c) = gerepileupto(av, gadd(gel(V, c), gsubstvec_i(P, va, X)));
  }
  return gmul2n(V, 1);
}

GEN
mffromqf(GEN Q, GEN P)
{
  pari_sp av = avma;
  GEN G, Qi, F, D, N, mf, v, gk, chi;
  long m, d, space;
  if (typ(Q) != t_MAT) pari_err_TYPE("mffromqf", Q);
  if (!RgM_is_ZM(Q) || !qfiseven(Q))
    pari_err_TYPE("mffromqf [not integral or even]", Q);
  m = lg(Q)-1;
  Qi = ZM_inv(Q, &N);
  if (!qfiseven(Qi)) N = shifti(N, 1);
  d = 0;
  if (!P || gequal1(P)) P = NULL;
  else
  {
    P = simplify_shallow(P);
    if (typ(P) == t_POL)
    {
      d = polishomogeneous(P);
      if (d < 0) pari_err_TYPE("mffromqf [not homogeneous t_POL]", P);
      if (!RgX_isspherical(Qi, P))
        pari_err_TYPE("mffromqf [not a spherical t_POL]", P);
    }
  }
  gk = uutoQ(m + 2*d, 2);
  D = ZM_det(Q);
  if (!odd(m)) { if ((m & 3) == 2) D = negi(D); } else D = shifti(D, 1);
  space = d > 0 ? mf_CUSP : mf_FULL;
  G = znstar0(N,1);
  chi = mkvec2(G, znchar_quad(G,D));
  mf = mfinit(mkvec3(N, gk, chi), space);
  if (odd(d))
  {
    F = mftrivial();
    v = zerocol(MF_get_dim(mf));
  }
  else
  {
    F = c_QF_i(mfsturm(mf), Q, P);
    v = mftobasis_i(mf, F);
    F = mflinear(mf, v);
  }
  return gerepilecopy(av, mkvec3(mf, F, v));
}

/***********************************************************************/
/*                          Eisenstein Series                          */
/***********************************************************************/
/* \sigma_{k-1}(\chi,n) */
static GEN
sigchi(long k, GEN CHI, long n)
{
  pari_sp av = avma;
  GEN S = gen_1, D = mydivisorsu(u_ppo(n,mfcharmodulus(CHI)));
  long i, l = lg(D), ord = mfcharorder(CHI), vt = varn(mfcharpol(CHI));
  for (i = 2; i < l; i++) /* skip D[1] = 1 */
  {
    long d = D[i], a = mfcharevalord(CHI, d, ord);
    S = gadd(S, Qab_Czeta(a, ord, powuu(d, k-1), vt));
  }
  return gerepileupto(av,S);
}

/* write n = n0*n1*n2, (n0,N1*N2) = 1, n1 | N1^oo, n2 | N2^oo;
 * return NULL if (n,N1,N2) > 1, else return factoru(n0) */
static GEN
sigchi2_dec(long n, long N1, long N2, long *pn1, long *pn2)
{
  GEN P0, E0, P, E, fa = myfactoru(n);
  long i, j, l;
  *pn1 = 1;
  *pn2 = 1;
  if (N1 == 1 && N2 == 1) return fa;
  P = gel(fa,1); l = lg(P);
  E = gel(fa,2);
  P0 = cgetg(l, t_VECSMALL);
  E0 = cgetg(l, t_VECSMALL);
  for (i = j = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (N1 % p == 0)
    {
      if (N2 % p == 0) return NULL;
      *pn1 *= upowuu(p,e);
    }
    else if (N2 % p == 0)
      *pn2 *= upowuu(p,e);
    else { P0[j] = p; E0[j] = e; j++; }
  }
  setlg(P0, j);
  setlg(E0, j); return mkvec2(P0,E0);
}

/* sigma_{k-1}(\chi_1,\chi_2,n), ord multiple of lcm(ord(CHI1),ord(CHI2)) */
static GEN
sigchi2(long k, GEN CHI1, GEN CHI2, long n, long ord)
{
  pari_sp av = avma;
  GEN S, D;
  long i, l, n1, n2, vt, N1 = mfcharmodulus(CHI1), N2 = mfcharmodulus(CHI2);
  D = sigchi2_dec(n, N1, N2, &n1, &n2); if (!D) return gc_const(av, gen_0);
  D = divisorsu_fact(D); l = lg(D);
  vt = varn(mfcharpol(CHI1));
  for (i = 1, S = gen_0; i < l; i++)
  { /* S += d^(k-1)*chi1(d)*chi2(n/d) */
    long a, d = n2*D[i], nd = n1*D[l-i]; /* (d,N1)=1; (n/d,N2) = 1 */
    a = mfcharevalord(CHI1, d, ord) + mfcharevalord(CHI2, nd, ord);
    if (a >= ord) a -= ord;
    S = gadd(S, Qab_Czeta(a, ord, powuu(d, k-1), vt));
  }
  return gerepileupto(av, S);
}

/**************************************************************************/
/**           Dirichlet characters with precomputed values               **/
/**************************************************************************/
/* CHI mfchar */
static GEN
mfcharcxinit(GEN CHI, long prec)
{
  GEN G = gel(CHI,1), chi = gel(CHI,2), z, V;
  GEN v = ncharvecexpo(G, znconrey_normalized(G,chi));
  long n, l = lg(v), o = mfcharorder(CHI);
  V = cgetg(l, t_VEC);
  z = grootsof1(o, prec); /* Mod(t, Phi_o(t)) -> e(1/o) */
  for (n = 1; n < l; n++) gel(V,n) = v[n] < 0? gen_0: gel(z, v[n]+1);
  return mkvecn(6, G, chi, gmfcharorder(CHI), v, V, mfcharpol(CHI));
}
/* v a "CHIvec" */
static long
CHIvec_N(GEN v) { return itou(znstar_get_N(gel(v,1))); }
static GEN
CHIvec_CHI(GEN v)
{ return mkvec4(gel(v,1), gel(v,2), gel(v,3), gel(v,6)); }
/* character order */
static long
CHIvec_ord(GEN v) { return itou(gel(v,3)); }
/* character exponents, i.e. t such that chi(n) = e(t) */
static GEN
CHIvec_expo(GEN v) { return gel(v,4); }
/* character values chi(n) */
static GEN
CHIvec_val(GEN v) { return gel(v,5); }
/* CHI(n) */
static GEN
mychareval(GEN v, long n)
{
  long N = CHIvec_N(v), ind = n%N;
  if (ind <= 0) ind += N;
  return gel(CHIvec_val(v), ind);
}
/* return c such that CHI(n) = e(c / ordz) or -1 if (n,N) > 1 */
static long
mycharexpo(GEN v, long n)
{
  long N = CHIvec_N(v), ind = n%N;
  if (ind <= 0) ind += N;
  return CHIvec_expo(v)[ind];
}
/* faster than mfcharparity */
static long
CHIvec_parity(GEN v) { return mycharexpo(v,-1) ? -1: 1; }
/**************************************************************************/

static ulong
sigchi2_Fl(long k, GEN CHI1vec, GEN CHI2vec, long n, GEN vz, ulong p)
{
  pari_sp av = avma;
  long ordz = lg(vz)-2, i, l, n1, n2;
  ulong S = 0;
  GEN D = sigchi2_dec(n, CHIvec_N(CHI1vec), CHIvec_N(CHI2vec), &n1, &n2);
  if (!D) return gc_ulong(av,S);
  D = divisorsu_fact(D);
  l = lg(D);
  for (i = 1; i < l; i++)
  { /* S += d^(k-1)*chi1(d)*chi2(n/d) */
    long a, d = n2*D[i], nd = n1*D[l-i]; /* (d,N1)=1, (n/d,N2)=1 */
    a = mycharexpo(CHI2vec, nd) + mycharexpo(CHI1vec, d);
    if (a >= ordz) a -= ordz;
    S = Fl_add(S, Qab_Czeta_Fl(a, vz, Fl_powu(d,k-1,p), p), p);
  }
  return gc_ulong(av,S);
}

/**********************************************************************/
/* Fourier expansions of Eisenstein series                            */
/**********************************************************************/
/* L(CHI_t,0) / 2, CHI_t(n) = CHI(n)(t/n) as a character modulo N*t,
 * order(CHI) | ord != 0 */
static GEN
charLFwt1(long N, GEN CHI, long ord, long t)
{
  GEN S;
  long r, vt;

  if (N == 1 && t == 1) return mkfrac(gen_m1,stoi(4));
  S = gen_0; vt = varn(mfcharpol(CHI));
  for (r = 1; r < N; r++)
  { /* S += r*chi(r) */
    long a, c;
    if (ugcd(N,r) != 1) continue;
    a = mfcharevalord(CHI,r,ord);
    c = (t != 1 && kross(t, r) < 0)? -r: r;
    S = gadd(S, Qab_Czeta(a, ord, stoi(c), vt));
  }
  return gdivgs(S, -2*N);
}
/* L(CHI,0) / 2, mod p */
static ulong
charLFwt1_Fl(GEN CHIvec, GEN vz, ulong p)
{
  long r, m = CHIvec_N(CHIvec);
  ulong S;
  if (m == 1) return Rg_to_Fl(mkfrac(gen_m1,stoi(4)), p);
  S = 0;
  for (r = 1; r < m; r++)
  { /* S += r*chi(r) */
    long a = mycharexpo(CHIvec,r);
    if (a < 0) continue;
    S = Fl_add(S, Qab_Czeta_Fl(a, vz, r, p), p);
  }
  return Fl_div(Fl_neg(S,p), 2*m, p);
}
/* L(CHI_t,1-k) / 2, CHI_t(n) = CHI(n) * (t/n), order(CHI) | ord != 0;
 * assume conductor of CHI_t divides N */
static GEN
charLFwtk(long N, long k, GEN CHI, long ord, long t)
{
  GEN S, P, dS;
  long r, vt;

  if (k == 1) return charLFwt1(N, CHI, ord, t);
  if (N == 1 && t == 1) return gdivgs(bernfrac(k),-2*k);
  vt = varn(mfcharpol(CHI));
  P = bern_init(N, k, &dS);
  dS = mul_denom(dS, stoi(-2*N*k));
  for (r = 1, S = gen_0; r < N; r++)
  { /* S += P(r)*chi(r) */
    long a;
    GEN C;
    if (ugcd(r,N) != 1) continue;
    a = mfcharevalord(CHI,r,ord);
    C = ZX_Z_eval(P, utoi(r));
    if (t != 1 && kross(t, r) < 0) C = gneg(C);
    S = gadd(S, Qab_Czeta(a, ord, C, vt));
  }
  return gdiv(S, dS);
}
/* L(CHI,1-k) / 2, mod p */
static ulong
charLFwtk_Fl(long k, GEN CHIvec, GEN vz, ulong p)
{
  GEN P, dS;
  long r, m;
  ulong S, d;
  if (k == 1) return charLFwt1_Fl(CHIvec, vz, p);
  m = CHIvec_N(CHIvec);
  if (m == 1) return Rg_to_Fl(gdivgs(bernfrac(k),-2*k), p);
  P = ZX_to_Flx(bern_init(m, k, &dS), p);
  for (r = 1, S = 0; r < m; r++)
  { /* S += P(r)*chi(r) */
    long a = mycharexpo(CHIvec,r);
    if (a < 0) continue;
    S = Fl_add(S, Qab_Czeta_Fl(a, vz, Flx_eval(P,r,p), p), p);
  }
  d = (2 * k * m) % p; if (dS) d = Fl_mul(d, umodiu(dS, p), p);
  return Fl_div(Fl_neg(S,p), d, p);
}

static GEN
mfeisenstein2_0(long k, GEN CHI1, GEN CHI2, long ord)
{
  long N1 = mfcharmodulus(CHI1), N2 = mfcharmodulus(CHI2);
  if (k == 1 && N1 == 1) return charLFwtk(N2, 1, CHI2, ord, 1);
  if (N2 == 1) return charLFwtk(N1, k, CHI1, ord, 1);
  return gen_0;
}
static ulong
mfeisenstein2_0_Fl(long k, GEN CHI1vec, GEN CHI2vec, GEN vz, ulong p)
{
  if (k == 1 && CHIvec_N(CHI1vec) == 1)
    return charLFwtk_Fl(k, CHI2vec, vz, p);
  else if (CHIvec_N(CHI2vec) == 1)
    return charLFwtk_Fl(k, CHI1vec, vz, p);
  else return 0;
}
static GEN
NK_eisen2(long k, GEN CHI1, GEN CHI2, long ord)
{
  long o, N = mfcharmodulus(CHI1)*mfcharmodulus(CHI2);
  GEN CHI = mfcharmul(CHI1, CHI2);
  o = mfcharorder(CHI);
  if ((ord & 3) == 2) ord >>= 1;
  if ((o & 3) == 2) o >>= 1;
  if (ord != o) pari_err_IMPL("mfeisenstein for these characters");
  return mkNK(N, k, CHI);
}
static GEN
mfeisenstein_i(long k, GEN CHI1, GEN CHI2)
{
  long s = 1, ord, vt;
  GEN E0, NK, vchi, T;
  if (CHI2) { CHI2 = get_mfchar(CHI2); if (mfcharparity(CHI2) < 0) s = -s; }
  if (CHI1) { CHI1 = get_mfchar(CHI1); if (mfcharparity(CHI1) < 0) s = -s; }
  if (s != m1pk(k)) return mftrivial();
  if (!CHI1) CHI1 = mfchartrivial();
  if (!CHI2)
  { /* E_k(chi1) */
    vt = varn(mfcharpol(CHI1));
    ord = mfcharorder(CHI1);
    NK = mkNK(mfcharmodulus(CHI1), k, CHI1);
    E0 = charLFwtk(mfcharmodulus(CHI1), k, CHI1, ord, 1);
    vchi = mkvec3(E0, mkvec(mfcharpol(CHI1)), CHI1);
    return tag(t_MF_EISEN, NK, vchi);
  }
  /* E_k(chi1,chi2) */
  vt = varn(mfcharpol(CHI1));
  ord = ulcm(mfcharorder(CHI1), mfcharorder(CHI2));
  NK = NK_eisen2(k, CHI1, CHI2, ord);
  E0 = mfeisenstein2_0(k, CHI1, CHI2, ord);
  T = mkvec(polcyclo(ord, vt));
  vchi = mkvec4(E0, T, CHI1, CHI2);
  return tag2(t_MF_EISEN, NK, vchi, mkvecsmall2(ord,0));
}
GEN
mfeisenstein(long k, GEN CHI1, GEN CHI2)
{
  pari_sp av = avma;
  if (k < 1) pari_err_DOMAIN("mfeisenstein", "k", "<", gen_1, stoi(k));
  return gerepilecopy(av, mfeisenstein_i(k, CHI1, CHI2));
}

static GEN
mfeisenstein2all(long N0, GEN NK, long k, GEN CHI1, GEN CHI2, GEN T, long o)
{
  GEN E, E0 = mfeisenstein2_0(k, CHI1,CHI2, o), vchi = mkvec4(E0, T, CHI1,CHI2);
  long j, d = (lg(T)==4)? itou(gmael(T,3,1)): 1;
  E = cgetg(d+1, t_VEC);
  for (j=1; j<=d; j++) gel(E,j) = tag2(t_MF_EISEN, NK,vchi,mkvecsmall2(o,j-1));
  return mfbdall(E, N0 / mf_get_N(gel(E,1)));
}

/* list of characters on G = (Z/NZ)^*, v[i] = NULL if (i,N) > 1, else
 * the conductor of Conrey label i, [conductor, primitive char].
 * Trivial chi (label 1) comes first */
static GEN
zncharsG(GEN G)
{
  long i, l, N = itou(znstar_get_N(G));
  GEN vCHI, V;
  if (N == 1) return mkvec2(gen_1,cgetg(1,t_COL));
  vCHI = const_vec(N,NULL);
  V = cyc2elts(znstar_get_conreycyc(G));
  l = lg(V);
  for (i = 1; i < l; i++)
  {
    GEN chi0, chi = zc_to_ZC(gel(V,i)), n, F;
    F = znconreyconductor(G, chi, &chi0);
    if (typ(F) != t_INT) F = gel(F,1);
    n = znconreyexp(G, chi);
    gel(vCHI, itos(n)) = mkvec2(chi0, F);
  }
  return vCHI;
}

/* CHI primitive, f(CHI) | N. Return pairs (CHI1,CHI2) both primitive
 * such that f(CHI1)*f(CHI2) | N and CHI1 * CHI2 = CHI;
 * if k = 1, CHI1 is even; if k = 2, omit (1,1) if CHI = 1 */
static GEN
mfeisensteinbasis_i(long N0, long k, GEN CHI)
{
  GEN G = gel(CHI,1), chi = gel(CHI,2), vT = const_vec(myeulerphiu(N0), NULL);
  GEN CHI0, GN, chiN, Lchi, LG, V, RES, NK, T, C = mfcharpol(CHI);
  long i, j, l, n, n1, N, ord = mfcharorder(CHI);
  long F = mfcharmodulus(CHI), vt = varn(mfcharpol(CHI));

  CHI0 = (F == 1)? CHI: mfchartrivial();
  j = 1; RES = cgetg(N0+1, t_VEC);
  T = gel(vT,ord) = Qab_trace_init(ord, ord, C, C);
  if (F != 1 || k != 2)
  { /* N1 = 1 */
    NK = mkNK(F, k, CHI);
    gel(RES, j++) = mfeisenstein2all(N0, NK, k, CHI0, CHI, T, ord);
    if (F != 1 && k != 1)
      gel(RES, j++) = mfeisenstein2all(N0, NK, k, CHI, CHI0, T, ord);
  }
  if (N0 == 1) { setlg(RES,j); return RES; }
  GN = G; chiN = chi;
  if (F == N0) N = N0;
  else
  {
    GEN faN = myfactoru(N0), P = gel(faN,1), E = gel(faN,2);
    long lP = lg(P);
    for (i = N = 1; i < lP; i++)
    {
      long p = P[i];
      N *= upowuu(p, maxuu(E[i]/2, z_lval(F,p)));
    }
    if ((N & 3) == 2) N >>= 1;
    if (N == 1) { setlg(RES,j); return RES; }
    if (F != N)
    {
      GN = znstar0(utoipos(N),1);
      chiN = zncharinduce(G, chi, GN);
    }
  }
  LG = const_vec(N, NULL); /* LG[d] = znstar(d,1) or NULL */
  gel(LG,1) = gel(CHI0,1);
  gel(LG,F) = G;
  gel(LG,N) = GN;
  Lchi = coprimes_zv(N);
  n = itou(znconreyexp(GN,chiN));
  V = zncharsG(GN); l = lg(V);
  for (n1 = 2; n1 < l; n1++) /* skip 1 (trivial char) */
  {
    GEN v = gel(V,n1), w, chi1, chi2, G1, G2, CHI1, CHI2;
    long N12, N1, N2, no, o12, t, m;
    if (!Lchi[n1] || n1 == n) continue; /* skip trivial chi2 */
    chi1 = gel(v,1); N1 = itou(gel(v,2)); /* conductor of chi1 */
    w = gel(V, Fl_div(n,n1,N));
    chi2 = gel(w,1); N2 = itou(gel(w,2)); /* conductor of chi2 */
    N12 = N1 * N2;
    if (N0 % N12) continue;

    G1 = gel(LG,N1); if (!G1) gel(LG,N1) = G1 = znstar0(utoipos(N1), 1);
    if (k == 1 && zncharisodd(G1,chi1)) continue;
    G2 = gel(LG,N2); if (!G2) gel(LG,N2) = G2 = znstar0(utoipos(N2), 1);
    CHI1 = mfcharGL(G1, chi1);
    CHI2 = mfcharGL(G2, chi2);
    o12 = ulcm(mfcharorder(CHI1), mfcharorder(CHI2));
    /* remove Galois orbit: same trace */
    no = Fl_powu(n1, ord, N);
    for (t = 1+ord, m = n1; t <= o12; t += ord)
    { /* m <-> CHI1^t, if t in Gal(Q(chi1,chi2)/Q), omit (CHI1^t,CHI2^t) */
      m = Fl_mul(m, no, N); if (!m) break;
      if (ugcd(t, o12) == 1) Lchi[m] = 0;
    }
    T = gel(vT,o12);
    if (!T) T = gel(vT,o12) = Qab_trace_init(o12, ord, polcyclo(o12,vt), C);
    NK = mkNK(N12, k, CHI);
    gel(RES, j++) = mfeisenstein2all(N0, NK, k, CHI1, CHI2, T, o12);
  }
  setlg(RES,j); return RES;
}

static GEN
mfbd_E2(GEN E2, long d, GEN CHI)
{
  GEN E2d = mfbd_i(E2, d);
  GEN F = mkvec2(E2, E2d), L = mkvec2(gen_1, utoineg(d));
  /* cannot use mflinear_i: E2 and E2d do not have the same level */
  return tag3(t_MF_LINEAR, mkNK(d,2,CHI), F, L, gen_1);
}
/* C-basis of E_k(Gamma_0(N),chi). If k = 1, the first basis element must not
 * vanish at oo [used in mf1basis]. Here E_1(CHI), whose q^0 coefficient
 * does not vanish (since L(CHI,0) does not) *if* CHI is not trivial; which
 * must be the case in weight 1.
 *
 * (k>=3): In weight k >= 3, basis is B(d) E(CHI1,(CHI/CHI1)_prim), where
 * CHI1 is primitive modulo N1, and if N2 is the conductor of CHI/CHI1
 * then d*N1*N2 | N.
 * (k=2): In weight k=2, same if CHI is nontrivial. If CHI is trivial, must
 * not take CHI1 trivial, and must add E_2(tau)-dE_2(d tau)), where
 * d|N, d > 1.
 * (k=1): In weight k=1, same as k >= 3 except that we restrict to CHI1 even */
static GEN
mfeisensteinbasis(long N, long k, GEN CHI)
{
  long i, F;
  GEN L;
  if (badchar(N, k, CHI)) return cgetg(1, t_VEC);
  if (k == 0) return mfcharistrivial(CHI)? mkvec(mf1()): cgetg(1, t_VEC);
  CHI = mfchartoprimitive(CHI, &F);
  L = mfeisensteinbasis_i(N, k, CHI);
  if (F == 1 && k == 2)
  {
    GEN v, E2 = mfeisenstein(2, NULL, NULL), D = mydivisorsu(N);
    long nD = lg(D)-1;
    v = cgetg(nD, t_VEC); L = vec_append(L,v);
    for (i = 1; i < nD; i++) gel(v,i) = mfbd_E2(E2, D[i+1], CHI);
  }
  return lg(L) == 1? L: shallowconcat1(L);
}

static GEN
not_in_space(GEN F, long flag)
{
  if (!flag) err_space(F);
  return cgetg(1, t_COL);
}
/* when flag set, no error */
GEN
mftobasis(GEN mf, GEN F, long flag)
{
  pari_sp av2, av = avma;
  GEN G, v, y, gk;
  long N, B, ismf = checkmf_i(F);

  mf = checkMF(mf);
  if (ismf)
  {
    if (mfistrivial(F)) return zerocol(MF_get_dim(mf));
    if (!mf_same_k(mf, F) || !mf_same_CHI(mf, F)) return not_in_space(F, flag);
  }
  N = MF_get_N(mf);
  gk = MF_get_gk(mf);
  if (ismf)
  {
    long NF = mf_get_N(F);
    B = maxuu(mfsturmNgk(NF,gk), mfsturmNgk(N,gk)) + 1;
    v = mfcoefs_i(F,B,1);
  }
  else
  {
    B = mfsturmNgk(N, gk) + 1;
    switch(typ(F))
    { /* F(0),...,F(lg(v)-2) */
      case t_SER: v = sertocol(F); settyp(v,t_VEC); break;
      case t_VEC: v = F; break;
      case t_COL: v = shallowtrans(F); break;
      default: pari_err_TYPE("mftobasis",F);
               v = NULL;/*LCOV_EXCL_LINE*/
    }
    if (flag) B = minss(B, lg(v)-2);
  }
  y = mftobasis_i(mf, v);
  if (typ(y) == t_VEC)
  {
    if (flag) return gerepilecopy(av, y);
    pari_err(e_MISC, "not enough coefficients in mftobasis");
  }
  av2 = avma;
  if (MF_get_space(mf) == mf_FULL || mfsturm(mf)+1 == B) return y;
  G = mflinear(mf, y);
  if (!gequal(v, mfcoefs_i(G, lg(v)-2,1))) y = NULL;
  if (!y) { set_avma(av); return not_in_space(F, flag); }
  set_avma(av2); return gerepileupto(av, y);
}

/* assume N > 0; first cusp is always 0 */
static GEN
mfcusps_i(long N)
{
  long i, c, l;
  GEN D, v;

  if (N == 1) return mkvec(gen_0);
  D = mydivisorsu(N); l = lg(D); /* left on stack */
  c = mfnumcuspsu_fact(myfactoru(N));
  v = cgetg(c + 1, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    long C = D[i], NC = D[l-i], lima = ugcd(C, NC), A0, A;
    for (A0 = 0; A0 < lima; A0++)
      if (ugcd(A0, lima) == 1)
      {
        A = A0; while (ugcd(A,C) > 1) A += lima;
        gel(v, c++) = uutoQ(A, C);
      }
  }
  return v;
}
/* List of cusps of Gamma_0(N) */
GEN
mfcusps(GEN gN)
{
  long N;
  GEN mf;
  if (typ(gN) == t_INT) N = itos(gN);
  else if ((mf = checkMF_i(gN))) N = MF_get_N(mf);
  else { pari_err_TYPE("mfcusps", gN); N = 0; }
  if (N <= 0) pari_err_DOMAIN("mfcusps", "N", "<=", gen_0, stoi(N));
  return mfcusps_i(N);
}

long
mfcuspisregular(GEN NK, GEN cusp)
{
  long v, N, dk, nk, t, o;
  GEN mf, CHI, go, A, C, g, c, d;
  if ((mf = checkMF_i(NK)))
  {
    GEN gk = MF_get_gk(mf);
    N = MF_get_N(mf);
    CHI = MF_get_CHI(mf);
    Qtoss(gk, &nk, &dk);
  }
  else
    checkNK2(NK, &N, &nk, &dk, &CHI, 0);
  if (typ(cusp) == t_INFINITY) return 1;
  if (typ(cusp) == t_FRAC) { A = gel(cusp,1); C = gel(cusp,2); }
  else { A = cusp; C = gen_1; }
  g = diviuexact(mului(N,C), ugcd(N, Fl_sqr(umodiu(C,N), N)));
  c = mulii(negi(C),g);
  d = addiu(mulii(A,g), 1);
  if (!CHI) return 1;
  go = gmfcharorder(CHI);
  v = vali(go); if (v < 2) go = shifti(go, 2-v);
  t = itou( znchareval(gel(CHI,1), gel(CHI,2), d, go) );
  if (dk == 1) return t == 0;
  o = itou(go);
  if (kronecker(c,d) < 0) t = Fl_add(t, o/2, o);
  if (Mod4(d) == 1) return t == 0;
  t = Fl_sub(t, Fl_mul(o/4, nk, o), o);
  return t == 0;
}

/* Some useful closures */

/* sum_{d|n} d^k */
static GEN
mysumdivku(ulong n, ulong k)
{
  GEN fa = myfactoru(n);
  return k == 1? usumdiv_fact(fa): usumdivk_fact(fa,k);
}
static GEN
c_Ek(long n, long d, GEN F)
{
  GEN E = cgetg(n + 2, t_VEC), C = gel(F,2);
  long i, k = mf_get_k(F);
  gel (E, 1) = gen_1;
  for (i = 1; i <= n; i++)
  {
    pari_sp av = avma;
    gel(E, i+1) = gerepileupto(av, gmul(C, mysumdivku(i*d, k-1)));
  }
  return E;
}

GEN
mfEk(long k)
{
  pari_sp av = avma;
  GEN E0, NK;
  if (k < 0 || odd(k)) pari_err_TYPE("mfEk [incorrect k]", stoi(k));
  if (!k) return mf1();
  E0 = gdivsg(-2*k, bernfrac(k));
  NK = mkNK(1,k,mfchartrivial());
  return gerepilecopy(av, tag(t_MF_Ek, NK, E0));
}

GEN
mfDelta(void)
{
  pari_sp av = avma;
  return gerepilecopy(av, tag0(t_MF_DELTA, mkNK(1,12,mfchartrivial())));
}

GEN
mfTheta(GEN psi)
{
  pari_sp av = avma;
  GEN N, gk, psi2;
  long par;
  if (!psi) { psi = mfchartrivial(); N = utoipos(4); par = 1; }
  else
  {
    long FC;
    psi = get_mfchar(psi);
    FC = mfcharconductor(psi);
    if (mfcharmodulus(psi) != FC)
      pari_err_TYPE("mfTheta [nonprimitive character]", psi);
    par = mfcharparity(psi);
    N = shifti(sqru(FC),2);
  }
  if (par > 0) { gk = ghalf; psi2 = psi; }
  else { gk = gsubsg(2, ghalf); psi2 = mfcharmul(psi, get_mfchar(stoi(-4))); }
  return gerepilecopy(av, tag(t_MF_THETA, mkgNK(N, gk, psi2, pol_x(1)), psi));
}

/* Output 0 if not desired eta product: if flag=0 (default) require
 * holomorphic at cusps. If flag set, accept meromorphic, but sill in some
 * modular function space */
GEN
mffrometaquo(GEN eta, long flag)
{
  pari_sp av = avma;
  GEN NK, N, k, BR, P;
  long v, cusp = 0;
  if (!etaquotype(&eta, &N,&k,&P, &v, NULL, flag? NULL: &cusp) || cusp < 0)
    return gc_const(av, gen_0);
  if (lg(gel(eta,1)) == 1) { set_avma(av); return mf1(); }
  BR = mkvec2(ZV_to_zv(gel(eta,1)), ZV_to_zv(gel(eta,2)));
  if (v < 0) v = 0;
  NK = mkgNK(N, k, get_mfchar(P), pol_x(1));
  return gerepilecopy(av, tag2(t_MF_ETAQUO, NK, BR, utoi(v)));
}

/* Q^(-r) */
static GEN
RgXn_negpow(GEN Q, long r, long L)
{
  if (r < 0) r = -r; else Q = RgXn_inv_i(Q, L);
  if (r != 1) Q = RgXn_powu_i(Q, r, L);
  return Q;
}
/* flag same as in mffrometaquo: if set, accept meromorphic. */
static GEN
mfisetaquo_i(GEN F, long flag)
{
  GEN gk, P, E, M, S, G, CHI, v, w;
  long b, l, L, N, vS, m, j;
  const long bextra = 10;

  if (!checkmf_i(F)) pari_err_TYPE("mfisetaquo",F);
  CHI = mf_get_CHI(F); if (mfcharorder(CHI) > 2) return NULL;
  N = mf_get_N(F);
  gk = mf_get_gk(F);
  b = mfsturmNgk(N, gk);
  L = maxss(N, b) + bextra;
  S = mfcoefs_i(F, L, 1);
  if (!RgV_is_ZV(S)) return NULL;
  for (vS = 1; vS <= L+1; vS++)
    if (signe(gel(S,vS))) break;
  vS--;
  if (vS >= bextra - 1) { L += vS; S = mfcoefs_i(F, L, 1); }
  if (vS) { S = vecslice(S, vS+1, L+1); L -= vS; }
  S = RgV_to_RgX(S, 0); l = lg(S)-2;
  P = cgetg(l, t_COL);
  E = cgetg(l, t_COL); w = v = gen_0; /* w = weight, v = valuation */
  for (m = j = 1; m+2 < lg(S); m++)
  {
    GEN c = gel(S,m+2);
    long r;
    if (is_bigint(c)) return NULL;
    r = -itos(c);
    if (r)
    {
      S = ZXn_mul(S, RgXn_negpow(eta_ZXn(m, L), r, L), L);
      gel(P,j) = utoipos(m);
      gel(E,j) = stoi(r);
      v = addmuliu(v, gel(E,j), m);
      w = addis(w, r);
      j++;
    }
  }
  if (!equalii(w, gmul2n(gk, 1)) || (!flag && !equalii(v, muluu(24,vS))))
    return NULL;
  setlg(P, j);
  setlg(E, j); M = mkmat2(P, E); G = mffrometaquo(M, flag);
  return (typ(G) != t_INT
          && (mfsturmmf(G) <= b + bextra || mfisequal(F, G, b)))? M: NULL;
}
GEN
mfisetaquo(GEN F, long flag)
{
  pari_sp av = avma;
  GEN M = mfisetaquo_i(F, flag);
  return M? gerepilecopy(av, M): gc_const(av, gen_0);
}

#if 0
/* number of primitive characters modulo N */
static ulong
numprimchars(ulong N)
{
  GEN fa, P, E;
  long i, l;
  ulong n;
  if ((N & 3) == 2) return 0;
  fa = myfactoru(N);
  P = gel(fa,1); l = lg(P);
  E = gel(fa,2);
  for (i = n = 1; i < l; i++)
  {
    ulong p = P[i], e = E[i];
    if (e == 2) n *= p-2; else n *= (p-1)*(p-1)*upowuu(p,e-2);
  }
  return n;
}
#endif

/* Space generated by products of two Eisenstein series */

static int
cmp_small_priority(void *E, GEN a, GEN b)
{
  GEN prio = (GEN)E;
  return cmpss(prio[(long)a], prio[(long)b]);
}
static long
znstar_get_expo(GEN G) { return itou(cyc_get_expo(znstar_get_cyc(G))); }

/* Return [vchi, bymod, vG]:
 * vG[f] = znstar(f,1) for f a conductor of (at least) a char mod N; else NULL
 * bymod[f] = vecsmall of conrey indexes of chars modulo f | N; else NULL
 * vchi[n] = a list of CHIvec [G0,chi0,o,ncharvecexpo(G0,nchi0),...]:
 *   chi0 = primitive char attached to Conrey Mod(n,N)
 * (resp. NULL if (n,N) > 1) */
static GEN
charsmodN(long N)
{
  GEN D, G, prio, phio, dummy = cgetg(1,t_VEC);
  GEN vP, vG = const_vec(N,NULL), vCHI  = const_vec(N,NULL);
  GEN bymod = const_vec(N,NULL);
  long pn, i, l, vt = fetch_user_var("t");
  D = mydivisorsu(N); l = lg(D);
  for (i = 1; i < l; i++)
    gel(bymod, D[i]) = vecsmalltrunc_init(myeulerphiu(D[i])+1);
  gel(vG,N) = G = znstar0(utoipos(N),1);
  pn = znstar_get_expo(G);  /* exponent(Z/NZ)^* */
  vP = const_vec(pn,NULL);
  for (i = 1; i <= N; i++)
  {
    GEN P, gF, G0, chi0, nchi0, chi, v, go;
    long j, F, o;
    if (ugcd(i,N) != 1) continue;
    chi = znconreylog(G, utoipos(i));
    gF = znconreyconductor(G, chi, &chi0);
    F = (typ(gF) == t_INT)? itou(gF): itou(gel(gF,1));
    G0 = gel(vG, F); if (!G0) G0 = gel(vG,F) = znstar0(gF, 1);
    nchi0 = znconreylog_normalize(G0,chi0);
    go = gel(nchi0,1); o = itou(go); /* order(chi0) */
    v = ncharvecexpo(G0, nchi0);
    if (!equaliu(go, pn)) v = zv_z_mul(v, pn / o);
    P = gel(vP, o); if (!P) P = gel(vP,o) = polcyclo(o,vt);
    /* mfcharcxinit with dummy complex powers */
    gel(vCHI,i) = mkvecn(6, G0, chi0, go, v, dummy, P);
    D = mydivisorsu(N / F); l = lg(D);
    for (j = 1; j < l; j++) vecsmalltrunc_append(gel(bymod, F*D[j]), i);
  }
  phio = zero_zv(pn); l = lg(vCHI); prio = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN CHI = gel(vCHI,i);
    long o;
    if (!CHI) continue;
    o = CHIvec_ord(CHI);
    if (!phio[o]) phio[o] = myeulerphiu(o);
    prio[i] = phio[o];
  }
  l = lg(bymod);
  /* sort characters by increasing value of phi(order) */
  for (i = 1; i < l; i++)
  {
    GEN z = gel(bymod,i);
    if (z) gen_sort_inplace(z, (void*)prio, &cmp_small_priority, NULL);
  }
  return mkvec3(vCHI, bymod, vG);
}

static GEN
mfeisenstein2pure(long k, GEN CHI1, GEN CHI2, long ord, GEN P, long lim)
{
  GEN c, V = cgetg(lim+2, t_COL);
  long n;
  c = mfeisenstein2_0(k, CHI1, CHI2, ord);
  if (P) c = grem(c, P);
  gel(V,1) = c;
  for (n=1; n <= lim; n++)
  {
    c = sigchi2(k, CHI1, CHI2, n, ord);
    if (P) c = grem(c, P);
    gel(V,n+1) = c;
  }
  return V;
}
static GEN
mfeisenstein2pure_Fl(long k, GEN CHI1vec, GEN CHI2vec, GEN vz, ulong p, long lim)
{
  GEN V = cgetg(lim+2, t_VECSMALL);
  long n;
  V[1] = mfeisenstein2_0_Fl(k, CHI1vec, CHI2vec, vz, p);
  for (n=1; n <= lim; n++) V[n+1] = sigchi2_Fl(k, CHI1vec, CHI2vec, n, vz, p);
  return V;
}

static GEN
getcolswt2(GEN M, GEN D, ulong p)
{
  GEN R, v = gel(M,1);
  long i, l = lg(M) - 1;
  R = cgetg(l, t_MAT); /* skip D[1] = 1 */
  for (i = 1; i < l; i++)
  {
    GEN w = Flv_Fl_mul(gel(M,i+1), D[i+1], p);
    gel(R,i) = Flv_sub(v, w, p);
  }
  return R;
}
static GEN
expandbd(GEN V, long d)
{
  long L, n, nd;
  GEN W;
  if (d == 1) return V;
  L = lg(V)-1; W = zerocol(L); /* nd = n/d */
  for (n = nd = 0; n < L; n += d, nd++) gel(W, n+1) = gel(V, nd+1);
  return W;
}
static GEN
expandbd_Fl(GEN V, long d)
{
  long L, n, nd;
  GEN W;
  if (d == 1) return V;
  L = lg(V)-1; W = zero_Flv(L); /* nd = n/d */
  for (n = nd = 0; n < L; n += d, nd++) W[n+1] = V[nd+1];
  return W;
}
static void
getcols_i(GEN *pM, GEN *pvj, GEN gk, GEN CHI1vec, GEN CHI2vec, long NN1, GEN vz,
          ulong p, long lim)
{
  GEN CHI1 = CHIvec_CHI(CHI1vec), CHI2 = CHIvec_CHI(CHI2vec);
  long N2 = CHIvec_N(CHI2vec);
  GEN vj, M, D = mydivisorsu(NN1/N2);
  long i, l = lg(D), k = gk[2];
  GEN V = mfeisenstein2pure_Fl(k, CHI1vec, CHI2vec, vz, p, lim);
  M = cgetg(l, t_MAT);
  for (i = 1; i < l; i++) gel(M,i) = expandbd_Fl(V, D[i]);
  if (k == 2 && N2 == 1 && CHIvec_N(CHI1vec) == 1)
  {
    M = getcolswt2(M, D, p); l--;
    D = vecslice(D, 2, l);
  }
  *pM = M;
  *pvj = vj = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(vj,i) = mkvec4(gk, CHI1, CHI2, utoipos(D[i]));
}

/* find all CHI1, CHI2 mod N such that CHI1*CHI2 = CHI, f(CHI1)*f(CHI2) | N.
 * set M = mfcoefs(B_e E(CHI1,CHI2), lim), vj = [e,i1,i2] */
static void
getcols(GEN *pM, GEN *pv, long k, long nCHI, GEN allN, GEN vz, ulong p,
        long lim)
{
  GEN vCHI = gel(allN,1), gk = utoi(k);
  GEN M = cgetg(1,t_MAT), v = cgetg(1,t_VEC);
  long i1, N = lg(vCHI)-1;
  for (i1 = 1; i1 <= N; i1++)
  {
    GEN CHI1vec = gel(vCHI, i1), CHI2vec, M1, v1;
    long NN1, i2;
    if (!CHI1vec) continue;
    if (k == 1 && CHIvec_parity(CHI1vec) == -1) continue;
    NN1 = N/CHIvec_N(CHI1vec); /* N/f(chi1) */;
    i2 = Fl_div(nCHI,i1, N);
    if (!i2) i2 = 1;
    CHI2vec = gel(vCHI,i2);
    if (NN1 % CHIvec_N(CHI2vec)) continue; /* f(chi1)f(chi2) | N ? */
    getcols_i(&M1, &v1, gk, CHI1vec, CHI2vec, NN1, vz, p, lim);
    M = shallowconcat(M, M1);
    v = shallowconcat(v, v1);
  }
  *pM = M;
  *pv = v;
}

static void
update_Mj(GEN *M, GEN *vecj, GEN *pz, ulong p)
{
  GEN perm;
  *pz = Flm_indexrank(*M, p); perm = gel(*pz,2);
  *M = vecpermute(*M, perm);
  *vecj = vecpermute(*vecj, perm);
}
static int
getcolsgen(long dim, GEN *pM, GEN *pvj, GEN *pz, long k, long ell, long nCHI,
           GEN allN, GEN vz, ulong p, long lim)
{
  GEN vCHI = gel(allN,1), bymod = gel(allN,2), gell = utoi(ell);
  long i1, N = lg(vCHI)-1;
  long L = lim+1;
  if (lg(*pvj)-1 >= dim) update_Mj(pM, pvj, pz, p);
  if (lg(*pvj)-1 == dim) return 1;
  for (i1 = 1; i1 <= N; i1++)
  {
    GEN CHI1vec = gel(vCHI, i1), T;
    long par1, j, l, N1, NN1;

    if (!CHI1vec) continue;
    par1 = CHIvec_parity(CHI1vec);
    if (ell == 1 && par1 == -1) continue;
    if (odd(ell)) par1 = -par1;
    N1 = CHIvec_N(CHI1vec);
    NN1 = N/N1;
    T = gel(bymod, NN1); l = lg(T);
    for (j = 1; j < l; j++)
    {
      long i2 = T[j], l1, l2, j1, s, nC;
      GEN M, M1, M2, vj, vj1, vj2, CHI2vec = gel(vCHI, i2);
      if (CHIvec_parity(CHI2vec) != par1) continue;
      nC = Fl_div(nCHI, Fl_mul(i1,i2,N), N);
      getcols(&M2, &vj2, k-ell, nC, allN, vz, p, lim);
      l2 = lg(M2); if (l2 == 1) continue;
      getcols_i(&M1, &vj1, gell, CHI1vec, CHI2vec, NN1, vz, p, lim);
      l1 = lg(M1);
      M1 = Flm_to_FlxV(M1, 0);
      M2 = Flm_to_FlxV(M2, 0);
      M  = cgetg((l1-1)*(l2-1) + 1, t_MAT);
      vj = cgetg((l1-1)*(l2-1) + 1, t_VEC);
      for (j1 = s = 1; j1 < l1; j1++)
      {
        GEN E = gel(M1,j1), v = gel(vj1,j1);
        long j2;
        for (j2 = 1; j2 < l2; j2++, s++)
        {
          GEN c = Flx_to_Flv(Flxn_mul(E, gel(M2,j2), L, p), L);
          gel(M,s) = c;
          gel(vj,s) = mkvec2(v, gel(vj2,j2));
        }
      }
      *pM = shallowconcat(*pM, M);
      *pvj = shallowconcat(*pvj, vj);
      if (lg(*pvj)-1 >= dim) update_Mj(pM, pvj, pz, p);
      if (lg(*pvj)-1 == dim) return 1;
    }
  }
  if (ell == 1)
  {
    update_Mj(pM, pvj, pz, p);
    return (lg(*pvj)-1 == dim);
  }
  return 0;
}

static GEN
mkF2bd(long d, long lim)
{
  GEN V = zerovec(lim + 1);
  long n;
  gel(V, 1) = sstoQ(-1, 24);
  for (n = 1; n <= lim/d; n++) gel(V, n*d + 1) = mysumdivku(n, 1);
  return V;
}

static GEN
mkeisen(GEN E, long ord, GEN P, long lim)
{
  long k = itou(gel(E,1)), e = itou(gel(E,4));
  GEN CHI1 = gel(E,2), CHI2 = gel(E,3);
  if (k == 2 && mfcharistrivial(CHI1) && mfcharistrivial(CHI2))
    return gsub(mkF2bd(1,lim), gmulgu(mkF2bd(e,lim), e));
  else
  {
    GEN V = mfeisenstein2pure(k, CHI1, CHI2, ord, P, lim);
    return expandbd(V, e);
  }
}
static GEN
mkM(GEN vj, long pn, GEN P, long lim)
{
  long j, l = lg(vj), L = lim+1;
  GEN M = cgetg(l, t_MAT);
  for (j = 1; j < l; j++)
  {
    GEN E1, E2;
    parse_vecj(gel(vj,j), &E1,&E2);
    E1 = RgV_to_RgX(mkeisen(E1, pn, P, lim), 0);
    if (E2)
    {
      E2 = RgV_to_RgX(mkeisen(E2, pn, P, lim), 0);
      E1 = RgXn_mul(E1, E2, L);
    }
    E1 = RgX_to_RgC(E1, L);
    if (P && E2) E1 = RgXQV_red(E1, P);
    gel(M,j) = E1;
  }
  return M;
}

/* assume N > 2 */
static GEN
mffindeisen1(long N)
{
  GEN G = znstar0(utoipos(N), 1), L = chargalois(G, NULL), chi0 = NULL;
  long j, m = N, l = lg(L);
  for (j = 1; j < l; j++)
  {
    GEN chi = gel(L,j);
    long r = myeulerphiu(itou(zncharorder(G,chi)));
    if (r >= m) continue;
    chi = znconreyfromchar(G, chi);
    if (zncharisodd(G,chi)) { m = r; chi0 = chi; if (r == 1) break; }
  }
  if (!chi0) pari_err_BUG("mffindeisen1 [no Eisenstein series found]");
  chi0 = znchartoprimitive(G,chi0);
  return mfcharGL(gel(chi0,1), gel(chi0,2));
}

static GEN
mfeisensteinspaceinit_i(long N, long k, GEN CHI)
{
  GEN M, Minv, vj, vG, GN, allN, P, vz, z = NULL;
  long nCHI, lim, ell, ord, dim = mffulldim(N, k, CHI);
  ulong r, p;

  if (!dim) retmkvec3(cgetg(1,t_VECSMALL),
                      mkvec2(cgetg(1,t_MAT),gen_1),cgetg(1,t_VEC));
  lim = mfsturmNk(N, k) + 1;
  allN = charsmodN(N);
  vG = gel(allN,3);
  GN = gel(vG,N);
  ord = znstar_get_expo(GN);
  P = ord <= 2? NULL: polcyclo(ord, varn(mfcharpol(CHI)));
  CHI = induce(GN, CHI); /* lift CHI mod N before mfcharno*/
  nCHI = mfcharno(CHI);
  r = QabM_init(ord, &p);
  vz = Fl_powers(r, ord, p);
  getcols(&M, &vj, k, nCHI, allN, vz, p, lim);
  for (ell = k>>1; ell >= 1; ell--)
    if (getcolsgen(dim, &M, &vj, &z, k, ell, nCHI, allN, vz, p, lim)) break;
  if (!z) update_Mj(&M, &vj, &z, p);
  if (lg(vj) - 1 < dim) return NULL;
  M = mkM(vj, ord, P, lim);
  Minv = QabM_Minv(rowpermute(M, gel(z,1)), P, ord);
  return mkvec4(gel(z,1), Minv, vj, utoi(ord));
}
/* true mf */
static GEN
mfeisensteinspaceinit(GEN mf)
{
  pari_sp av = avma;
  GEN z, CHI = MF_get_CHI(mf);
  long N = MF_get_N(mf), k = MF_get_k(mf);
  if (!CHI) CHI = mfchartrivial();
  z = mfeisensteinspaceinit_i(N, k, CHI);
  if (!z)
  {
    GEN E, CHIN = mffindeisen1(N), CHI0 = mfchartrivial();
    z = mfeisensteinspaceinit_i(N, k+1, mfcharmul(CHI, CHIN));
    if (z) E = mkvec4(gen_1, CHI0, CHIN, gen_1);
    else
    {
      z = mfeisensteinspaceinit_i(N, k+2, CHI);
      E = mkvec4(gen_2, CHI0, CHI0, utoipos(N));
    }
    z = mkvec2(z, E);
  }
  return gerepilecopy(av, z);
}

/* decomposition of modular form on eisenspace */
static GEN
mfeisensteindec(GEN mf, GEN F)
{
  pari_sp av = avma;
  GEN M, Mindex, Mvecj, V, B, CHI;
  long o, ord;

  Mvecj = obj_checkbuild(mf, MF_EISENSPACE, &mfeisensteinspaceinit);
  if (lg(Mvecj) < 5)
  {
    GEN E, e = gel(Mvecj,2), gkE = gel(e,1);
    long dE = itou(gel(e,4));
    Mvecj = gel(Mvecj,1);
    E = mfeisenstein(itou(gkE), NULL, gel(e,3));
    if (dE != 1) E = mfbd_E2(E, dE, gel(e,2)); /* here k = 2 */
    F = mfmul(F, E);
  }
  M = gel(Mvecj, 2);
  if (lg(M) == 1) return cgetg(1, t_VEC);
  Mindex = gel(Mvecj, 1);
  ord = itou(gel(Mvecj,4));
  V = mfcoefs(F, Mindex[lg(Mindex)-1]-1, 1); settyp(V, t_COL);
  CHI = mf_get_CHI(F);
  o = mfcharorder(CHI);
  if (o > 2 && o != ord)
  { /* convert Mod(.,polcyclo(o)) to Mod(., polcyclo(N)) for o | N,
     * o and N both != 2 (mod 4) */
    GEN z, P = gel(M,4); /* polcyclo(ord) */
    long vt = varn(P);
    z = gmodulo(pol_xn(ord/o, vt), P);
    if (ord % o) pari_err_TYPE("mfeisensteindec", V);
    V = gsubst(liftpol_shallow(V), vt, z);
  }
  B = Minv_RgC_mul(M, vecpermute(V, Mindex));
  return gerepileupto(av, B);
}

/*********************************************************************/
/*                        END EISENSPACE                             */
/*********************************************************************/

static GEN
sertocol2(GEN S, long l)
{
  GEN C = cgetg(l + 2, t_COL);
  long i;
  for (i = 0; i <= l; i++) gel(C, i+1) = polcoef_i(S, i, -1);
  return C;
}

/* Compute polynomial P0 such that F=E4^(k/4)P0(E6/E4^(3/2)). */
static GEN
mfcanfindp0(GEN F, long k)
{
  pari_sp ltop = avma;
  GEN E4, E6, V, V1, Q, W, res, M, B;
  long l, j;
  l = k/6 + 2;
  V = mfcoefsser(F,l);
  E4 = mfcoefsser(mfEk(4),l);
  E6 = mfcoefsser(mfEk(6),l);
  V1 = gdiv(V, gpow(E4, uutoQ(k,4), 0));
  Q = gdiv(E6, gpow(E4, uutoQ(3,2), 0));
  W = gpowers(Q, l - 1);
  M = cgetg(l + 1, t_MAT);
  for (j = 1; j <= l; j++) gel(M,j) = sertocol2(gel(W,j), l);
  B = sertocol2(V1, l);
  res = inverseimage(M, B);
  if (lg(res) == 1) err_space(F);
  return gerepilecopy(ltop, gtopolyrev(res, 0));
}

/* Compute the first n+1 Taylor coeffs at tau=I of a modular form
 * on SL_2(Z). */
GEN
mftaylor(GEN F, long n, long flreal, long prec)
{
  pari_sp ltop = avma;
  GEN P0, Pm1 = gen_0, v;
  GEN X2 = mkpoln(3, ghalf,gen_0,gneg(ghalf)); /* (x^2-1) / 2 */
  long k, m;
  if (!checkmf_i(F)) pari_err_TYPE("mftaylor",F);
  k = mf_get_k(F);
  if (mf_get_N(F) != 1 || k < 0) pari_err_IMPL("mftaylor for this form");
  P0 = mfcanfindp0(F, k);
  v = cgetg(n+2, t_VEC); gel(v, 1) = RgX_coeff(P0,0);
  for (m = 0; m < n; m++)
  {
    GEN P1 = gdivgu(gmulsg(-(k + 2*m), RgX_shift(P0,1)), 12);
    P1 = gadd(P1, gmul(X2, RgX_deriv(P0)));
    if (m) P1 = gsub(P1, gdivgu(gmulsg(m*(m+k-1), Pm1), 144));
    Pm1 = P0; P0 = P1;
    gel(v, m+2) = RgX_coeff(P0, 0);
  }
  if (flreal)
  {
    GEN pi2 = Pi2n(1, prec), pim4 = gmulsg(-2, pi2), VPC;
    GEN C = gmulsg(3, gdiv(gpowgs(ggamma(uutoQ(1,4), prec), 8), gpowgs(pi2, 6)));
    /* E_4(i): */
    GEN facn = gen_1;
    VPC = gpowers(gmul(pim4, gsqrt(C, prec)), n);
    C = gpow(C, uutoQ(k,4), prec);
    for (m = 0; m <= n; m++)
    {
      gel(v, m+1) = gdiv(gmul(C, gmul(gel(v, m+1), gel(VPC, m+1))), facn);
      facn = gmulgu(facn, m+1);
    }
  }
  return gerepilecopy(ltop, v);
}

#if 0
/* To be used in mfeigensearch() */
GEN
mfreadratfile()
{
  GEN eqn;
  pariFILE *F = pari_fopengz("rateigen300.gp");
  eqn = gp_readvec_stream(F->file);
  pari_fclose(F);
  return eqn;
}
#endif
 /*****************************************************************/
/*           EISENSTEIN CUSPS: COMPLEX DIRECTLY: one F_k         */
/*****************************************************************/

/* CHIvec = charinit(CHI); data = [N1g/g1,N2g/g2,g1/g,g2/g,C/g1,C/g2,
 * (N1g/g1)^{-1},(N2g/g2)^{-1}] */

/* nm = n/m;
 * z1 = powers of \z_{C/g}^{(Ae/g)^{-1}},
 * z2 = powers of \z_N^{A^{-1}(g1g2/C)}]
 * N.B. : we compute value and conjugate at the end, so it is (Ae/g)^{-1}
 * and not -(Ae/g)^{-1} */
static GEN
eiscnm(long nm, long m, GEN CHI1vec, GEN CHI2vec, GEN data, GEN z1)
{
  long Cg1 = data[5], s10 = (nm*data[7]) % Cg1, r10 = (nm - data[1]*s10) / Cg1;
  long Cg2 = data[6], s20 = (m *data[8]) % Cg2, r20 = (m  - data[2]*s20) / Cg2;
  long j1, r1, s1;
  GEN T = gen_0;
  for (j1 = 0, r1 = r10, s1 = s10; j1 < data[3]; j1++, r1 -= data[1], s1 += Cg1)
  {
    GEN c1 = mychareval(CHI1vec, r1);
    if (!gequal0(c1))
    {
      long j2, r2, s2;
      GEN S = gen_0;
      for (j2 = 0, r2 = r20, s2 = s20; j2 < data[4]; j2++, r2 -= data[2], s2 += Cg2)
      {
        GEN c2 = mychareval(CHI2vec, r2);
        if (!gequal0(c2)) S = gadd(S, gmul(c2, rootsof1pow(z1, s1*s2)));
      }
      T = gadd(T, gmul(c1, S));
    }
  }
  return conj_i(T);
}

static GEN
fg1g2n(long n, long k, GEN CHI1vec, GEN CHI2vec, GEN data, GEN z1, GEN z2)
{
  pari_sp av = avma;
  GEN S = gen_0, D = mydivisorsu(n);
  long i, l = lg(D);
  for (i = 1; i < l; i++)
  {
    long m = D[i], nm = D[l-i]; /* n/m */
    GEN u = eiscnm( nm,  m, CHI1vec, CHI2vec, data, z1);
    GEN v = eiscnm(-nm, -m, CHI1vec, CHI2vec, data, z1);
    GEN w = odd(k) ? gsub(u, v) : gadd(u, v);
    S = gadd(S, gmul(powuu(m, k-1), w));
  }
  return gerepileupto(av, gmul(S, rootsof1pow(z2, n)));
}

static GEN
gausssumcx(GEN CHIvec, long prec)
{
  GEN z, S, V;
  long m, N = CHIvec_N(CHIvec);
  if (N == 1) return gen_1;
  V = CHIvec_val(CHIvec);
  z = rootsof1u_cx(N, prec);
  S = gmul(z, gel(V, N));
  for (m = N-1; m >= 1; m--) S = gmul(z, gadd(gel(V, m), S));
  return S;
}

/* Computation of Q_k(\z_N^s) as a polynomial in \z_N^s. FIXME: explicit
 * formula ? */
static GEN
mfqk(long k, long N)
{
  GEN X, P, ZI, Q, Xm1, invden;
  long i;
  ZI = gdivgu(RgX_shift_shallow(RgV_to_RgX(identity_ZV(N-1), 0), 1), N);
  if (k == 1) return ZI;
  P = gsubgs(pol_xn(N,0), 1);
  invden = RgXQ_powu(ZI, k, P);
  X = pol_x(0); Q = gneg(X); Xm1 = gsubgs(X, 1);
  for (i = 2; i < k; i++)
    Q = RgX_shift_shallow(ZX_add(gmul(Xm1, ZX_deriv(Q)), gmulsg(-i, Q)), 1);
  return RgXQ_mul(Q, invden, P);
}

/* CHI mfchar; M is a multiple of the conductor of CHI, but is NOT
 * necessarily its modulus */
static GEN
mfskcx(long k, GEN CHI, long M, long prec)
{
  GEN S, CHIvec, P;
  long F, m, i, l;
  CHI = mfchartoprimitive(CHI, &F);
  CHIvec = mfcharcxinit(CHI, prec);
  if (F == 1) S = gdivgu(bernfrac(k), k);
  else
  {
    GEN Q = mfqk(k, F), V = CHIvec_val(CHIvec);
    S = gmul(gel(V, F), RgX_coeff(Q, 0));
    for (m = 1; m < F; m++) S = gadd(S, gmul(gel(V, m), RgX_coeff(Q, m)));
    S = conj_i(S);
  }
  /* prime divisors of M not dividing f(chi) */
  P = gel(myfactoru(u_ppo(M/F,F)), 1); l = lg(P);
  for (i = 1; i < l; i++)
  {
    long p = P[i];
    S = gmul(S, gsubsg(1, gdiv(mychareval(CHIvec, p), powuu(p, k))));
  }
  return gmul(gmul(gausssumcx(CHIvec, prec), S), powuu(M/F, k));
}

static GEN
f00_i(long k, GEN CHI1vec, GEN CHI2vec, GEN G2, GEN S, long prec)
{
  GEN c, a;
  long N1 = CHIvec_N(CHI1vec), N2 = CHIvec_N(CHI2vec);
  if (S[2] != N1) return gen_0;
  c = mychareval(CHI1vec, S[3]);
  if (isintzero(c)) return gen_0;
  a = mfskcx(k, mfchardiv(CHIvec_CHI(CHI2vec), CHIvec_CHI(CHI1vec)), N1*N2, prec);
  a = gmul(a, conj_i(gmul(c,G2)));
  return gdiv(a, mulsi(-N2, powuu(S[1], k-1)));
}

static GEN
f00(long k, GEN CHI1vec,GEN CHI2vec, GEN G1,GEN G2, GEN data, long prec)
{
  GEN T1, T2;
  T2 = f00_i(k, CHI1vec, CHI2vec, G2, data, prec);
  if (k > 1) return T2;
  T1 = f00_i(k, CHI2vec, CHI1vec, G1, data, prec);
  return gadd(T1, T2);
}

/* ga in SL_2(Z), find beta [a,b;c,d] in Gamma_0(N) and mu in Z such that
 * beta * ga * T^u = [A',B';C',D'] with C' | N and N | B', C' > 0 */
static void
mfgatogap(GEN ga, long N, long *pA, long *pC, long *pD, long *pd, long *pmu)
{
  GEN A = gcoeff(ga,1,1), B = gcoeff(ga,1,2);
  GEN C = gcoeff(ga,2,1), D = gcoeff(ga,2,2), a, b, c, d;
  long t, Ap, Cp, B1, D1, mu;
  Cp = itou(bezout(muliu(A,N), C, &c, &d)); /* divides N */
  t = 0;
  if (Cp > 1)
  { /* (d, N/Cp) = 1, find t such that (d - t*(A*N/Cp), N) = 1 */
    long dN = umodiu(d,Cp), Q = (N/Cp * umodiu(A,Cp)) % Cp;
    while (ugcd(dN, Cp) > 1) { t++; dN = Fl_sub(dN, Q, Cp); }
  }
  if (t)
  {
    c = addii(c, mului(t, diviuexact(C,Cp)));
    d = subii(d, mului(t, muliu(A, N/Cp))); /* (d,N) = 1 */
  }
  D1 = umodiu(mulii(d,D), N);
  (void)bezout(d, mulis(c,-N), &a, &b); /* = 1 */
  t = 0; Ap = umodiu(addii(mulii(a,A), mulii(b,C)), N); /* (Ap,Cp) = 1 */
  while (ugcd(Ap, N) > 1) { t++; Ap = Fl_add(Ap, Cp, N); }
  B1 = umodiu(a,N)*umodiu(B,N) + umodiu(b,N)*umodiu(D,N) + t*D1;
  B1 %= N;
  *pmu = mu = Fl_neg(Fl_div(B1, Ap, N), N);
  /* A', D' and d only needed modulo N */
  *pd = umodiu(d, N);
  *pA = Ap;
  *pC = Cp; *pD = (D1 + Cp*mu) % N;
}

#if 0
/* CHI is a mfchar, return alpha(CHI) */
static long
mfalchi(GEN CHI, long AN, long cg)
{
  GEN G = gel(CHI,1), chi = gel(CHI,2), go = gmfcharorder(CHI);
  long o = itou(go), a = itos( znchareval(G, chi, stoi(1 + AN/cg), go) );
  if (a < 0 || (cg * a) % o) pari_err_BUG("mfalchi");
  return (cg * a) / o;
}
#endif
/* return A such that CHI1(t) * CHI2(t) = e(A) or NULL if (t,N1*N2) > 1 */
static GEN
mfcharmuleval(GEN CHI1vec, GEN CHI2vec, long t)
{
  long a1 = mycharexpo(CHI1vec, t), o1 = CHIvec_ord(CHI1vec);
  long a2 = mycharexpo(CHI2vec, t), o2 = CHIvec_ord(CHI2vec);;
  if (a1 < 0 || a2 < 0) return NULL;
  return sstoQ(a1*o2 + a2*o1, o1*o2);
}
static GEN
mfcharmulcxeval(GEN CHI1vec, GEN CHI2vec, long t, long prec)
{
  GEN A = mfcharmuleval(CHI1vec, CHI2vec, t);
  long n, d;
  if (!A) return gen_0;
  Qtoss(A, &n,&d); return rootsof1q_cx(n, d, prec);
}
/* alpha(CHI1 * CHI2) */
static long
mfalchi2(GEN CHI1vec, GEN CHI2vec, long AN, long cg)
{
  GEN A = mfcharmuleval(CHI1vec, CHI2vec, 1 + AN/cg);
  long a;
  if (!A) pari_err_BUG("mfalchi2");
  A = gmulsg(cg, A);
  if (typ(A) != t_INT) pari_err_BUG("mfalchi2");
  a = itos(A) % cg; if (a < 0) a += cg;
  return a;
}

/* return g = (a,b), set u >= 0 s.t. g = a * u (mod b) */
static long
mybezout(long a, long b, long *pu)
{
  long junk, g = cbezout(a, b, pu, &junk);
  if (*pu < 0) *pu += b/g;
  return g;
}

/* E = [k, CHI1,CHI2, e], CHI1 and CHI2 primitive mfchars such that,
 * CHI1(-1)*CHI2(-1) = (-1)^k; expansion of (B_e (E_k(CHI1,CHI2))) | ga.
 * w is the width for the space of the calling function. */
static GEN
mfeisensteingacx(GEN E, long w, GEN ga, long lim, long prec)
{
  GEN CHI1vec, CHI2vec, CHI1 = gel(E,2), CHI2 = gel(E,3), v, S, ALPHA;
  GEN G1, G2, z1, z2, data;
  long k = itou(gel(E,1)), e = itou(gel(E,4));
  long N1 = mfcharmodulus(CHI1);
  long N2 = mfcharmodulus(CHI2), N = e * N1 * N2;
  long NsurC, cg, wN, A, C, Ai, d, mu, alchi, na, da;
  long eg, g, gH, U, u0, u1, u2, Aig, H, m, n, t, Cg, NC1, NC2;

  mfgatogap(ga, N, &A, &C, &Ai, &d, &mu);
  CHI1vec = mfcharcxinit(CHI1, prec);
  CHI2vec = mfcharcxinit(CHI2, prec);
  NsurC = N/C; cg  = ugcd(C, NsurC); wN = NsurC / cg;
  if (w%wN) pari_err_BUG("mfeisensteingacx [wN does not divide w]");
  alchi = mfalchi2(CHI1vec, CHI2vec, A*N, cg);
  ALPHA = sstoQ(alchi, NsurC);

  g = mybezout(A*e, C, &u0); Cg = C/g; eg = e/g;
  NC1 = mybezout(N1, Cg, &u1);
  NC2 = mybezout(N2, Cg, &u2);
  H = (NC1*NC2*g)/Cg;
  Aig = (Ai*H)%N; if (Aig < 0) Aig += N;
  z1 = rootsof1powinit(u0, Cg, prec);
  z2 = rootsof1powinit(Aig, N, prec);
  data = mkvecsmalln(8, N1/NC1, N2/NC2, NC1, NC2, Cg/NC1, Cg/NC2, u1, u2);
  v = zerovec(lim + 1);
  /* need n*H = alchi (mod cg) */
  gH = mybezout(H, cg, &U);
  if (gH > 1)
  {
    if (alchi % gH) return mkvec2(gen_0, v);
    alchi /= gH; cg /= gH; H /= gH;
  }
  G1 = gausssumcx(CHI1vec, prec);
  G2 = gausssumcx(CHI2vec, prec);
  if (!alchi)
    gel(v,1) = f00(k, CHI1vec,CHI2vec,G1,G2, mkvecsmall3(NC2,Cg,A*eg), prec);
  n = Fl_mul(alchi,U,cg); if (!n) n = cg;
  m = (n*H - alchi) / cg; /* positive, exact division */
  for (; m <= lim; n+=cg, m+=H)
    gel(v, m+1) = fg1g2n(n, k, CHI1vec, CHI2vec, data, z1,z2);
  t = (2*e)/g; if (odd(k)) t = -t;
  v = gdiv(v, gmul(conj_i(gmul(G1,G2)), mulsi(t, powuu(eg*N2/NC2, k-1))));
  if (k == 2 && N1 == 1 && N2 == 1) v = gsub(mkF2bd(wN,lim), gmulsg(e,v));

  Qtoss(ALPHA, &na,&da);
  S = conj_i( mfcharmulcxeval(CHI1vec,CHI2vec,d,prec) ); /* CHI(1/d) */
  if (wN > 1)
  {
    GEN z = rootsof1powinit(-mu, wN, prec);
    long i, l = lg(v);
    for (i = 1; i < l; i++) gel(v,i) = gmul(gel(v,i), rootsof1pow(z,i-1));
  }
  v = RgV_Rg_mul(v, gmul(S, rootsof1q_cx(-mu*na, da, prec)));
  return mkvec2(ALPHA, bdexpand(v, w/wN));
}

/*****************************************************************/
/*                       END EISENSTEIN CUSPS                    */
/*****************************************************************/

static GEN
mfchisimpl(GEN CHI)
{
  GEN G, chi;
  if (typ(CHI) == t_INT) return CHI;
  G = gel(CHI, 1); chi = gel(CHI, 2);
  switch(mfcharorder(CHI))
  {
    case 1: chi = gen_1; break;
    case 2: chi = znchartokronecker(G,chi,1); break;
    default:chi = mkintmod(znconreyexp(G,chi), znstar_get_N(G)); break;
  }
  return chi;
}

GEN
mfparams(GEN F)
{
  pari_sp av = avma;
  GEN z, mf, CHI;
  if ((mf = checkMF_i(F)))
  {
    long N = MF_get_N(mf);
    GEN gk = MF_get_gk(mf);
    CHI = MF_get_CHI(mf);
    z = mkvec5(utoi(N), gk, CHI, utoi(MF_get_space(mf)), mfcharpol(CHI));
  }
  else
  {
    if (!checkmf_i(F)) pari_err_TYPE("mfparams", F);
    z = vec_append(mf_get_NK(F), mfcharpol(mf_get_CHI(F)));
  }
  gel(z,3) = mfchisimpl(gel(z,3));
  return gerepilecopy(av, z);
}

GEN
mfisCM(GEN F)
{
  pari_sp av = avma;
  forprime_t S;
  GEN D, v;
  long N, k, lD, sb, p, i;
  if (!checkmf_i(F)) pari_err_TYPE("mfisCM", F);
  N = mf_get_N(F);
  k = mf_get_k(F); if (N < 0 || k < 0) pari_err_IMPL("mfisCM for this F");
  D = mfunram(N, -1);
  lD = lg(D);
  sb = maxss(mfsturmNk(N, k), 4*N);
  v = mfcoefs_i(F, sb, 1);
  u_forprime_init(&S, 2, sb);
  while ((p = u_forprime_next(&S)))
  {
    GEN ap = gel(v, p+1);
    if (!gequal0(ap))
      for (i = 1; i < lD; i++)
        if (kross(D[i], p) == -1) { D = vecsplice(D, i); lD--; }
  }
  if (lD == 1) return gc_const(av, gen_0);
  if (lD == 2) return gc_stoi(av, D[1]);
  if (k > 1) pari_err_BUG("mfisCM");
  return gerepileupto(av, zv_to_ZV(D));
}

static long
mfspace_i(GEN mf, GEN F)
{
  GEN v, vF, gk;
  long n, nE, i, l, s, N;

  mf = checkMF(mf); s = MF_get_space(mf);
  if (!F) return s;
  if (!checkmf_i(F)) pari_err_TYPE("mfspace",F);
  v = mftobasis(mf, F, 1);
  n = lg(v)-1; if (!n) return -1;
  nE = lg(MF_get_E(mf))-1;
  switch(s)
  {
    case mf_NEW: case mf_OLD: case mf_EISEN: return s;
    case mf_FULL:
      if (mf_get_type(F) == t_MF_THETA) return mf_EISEN;
      if (!gequal0(vecslice(v,1,nE)))
        return gequal0(vecslice(v,nE+1,n))? mf_EISEN: mf_FULL;
  }
  /* mf is mf_CUSP or mf_FULL, F a cusp form */
  gk = mf_get_gk(F);
  if (typ(gk) == t_FRAC || equali1(gk)) return mf_CUSP;
  vF = mftonew_i(mf, vecslice(v, nE+1, n), &N);
  if (N != MF_get_N(mf)) return mf_OLD;
  l = lg(vF);
  for (i = 1; i < l; i++)
    if (itos(gmael(vF,i,1)) != N) return mf_CUSP;
  return mf_NEW;
}
long
mfspace(GEN mf, GEN F)
{ pari_sp av = avma; return gc_long(av, mfspace_i(mf,F)); }
static GEN
lfunfindchi(GEN ldata, GEN van, long prec)
{
  GEN gN = ldata_get_conductor(ldata), gk = ldata_get_k(ldata);
  GEN G = znstar0(gN,1), cyc = znstar_get_conreycyc(G), L, go, vz;
  long N = itou(gN), odd = typ(gk) == t_INT && mpodd(gk);
  long i, j, o, l, B0 = 2, B = lg(van)-1, bit = 10 - prec2nbits(prec);

  /* if van is integral, chi must be trivial */
  if (typ(van) == t_VECSMALL) return mfcharGL(G, zerocol(lg(cyc)-1));
  L = cyc2elts(cyc); l = lg(L);
  for (i = j = 1; i < l; i++)
  {
    GEN chi = zc_to_ZC(gel(L,i));
    if (zncharisodd(G,chi) == odd) gel(L,j++) = mfcharGL(G,chi);
  }
  setlg(L,j); l = j;
  if (l <= 2) return gel(L,1);
  o = znstar_get_expo(G); go = utoi(o);
  vz = grootsof1(o, prec);
  for (;;)
  {
    long n;
    for (n = B0; n <= B; n++)
    {
      GEN an, r;
      long j;
      if (ugcd(n, N) != 1) continue;
      an = gel(van,n); if (gexpo(an) < bit) continue;
      r = gdiv(an, conj_i(an));
      for (i = 1; i < l; i++)
      {
        GEN CHI = gel(L,i);
        if (gexpo(gsub(r, gel(vz, znchareval_i(CHI,n,go)+1))) > bit)
          gel(L,i) = NULL;
      }
      for (i = j = 1; i < l; i++)
        if (gel(L,i)) gel(L,j++) = gel(L,i);
      l = j; setlg(L,l);
      if (l == 2) return gel(L,1);
    }
    B0 = B+1; B <<= 1;
    van = ldata_vecan(ldata_get_an(ldata), B, prec);
  }
}

GEN
mffromlfun(GEN L, long prec)
{
  pari_sp av = avma;
  GEN ldata = lfunmisc_to_ldata_shallow(L), Vga = ldata_get_gammavec(ldata);
  GEN van, a0, CHI, NK, gk = ldata_get_k(ldata);
  long N, space;
  if (!gequal(Vga, mkvec2(gen_0, gen_1))) pari_err_TYPE("mffromlfun", L);
  N = itou(ldata_get_conductor(ldata));
  van = ldata_vecan(ldata_get_an(ldata), mfsturmNgk(N,gk) + 2, prec);
  CHI = lfunfindchi(ldata, van, prec);
  if (typ(van) != t_VEC) van = vecsmall_to_vec_inplace(van);
  space = (lg(ldata) == 7)? mf_CUSP: mf_FULL;
  a0 = (space == mf_CUSP)? gen_0: gneg(lfun(L, gen_0, prec2nbits(prec)));
  NK = mkvec3(utoi(N), gk, mfchisimpl(CHI));
  return gerepilecopy(av, mkvec3(NK, utoi(space), shallowconcat(a0, van)));
}
/*******************************************************************/
/*                                                                 */
/*                       HALF-INTEGRAL WEIGHT                      */
/*                                                                 */
/*******************************************************************/
/* We use the prefix mf2; k represents the weight -1/2, so e.g.
   k = 2 is weight 5/2. N is the level, so 4\mid N, and CHI is the
   character, always even. */

static long
lamCO(long r, long s, long p)
{
  if ((s << 1) <= r)
  {
    long rp = r >> 1;
    if (odd(r)) return upowuu(p, rp) << 1;
    else return (p + 1)*upowuu(p, rp - 1);
  }
  else return upowuu(p, r - s) << 1;
}

static int
condC(GEN faN, GEN valF)
{
  GEN P = gel(faN, 1), E = gel(faN, 2);
  long l = lg(P), i;
  for (i = 1; i < l; i++)
    if ((P[i] & 3L) == 3)
    {
      long r = E[i];
      if (odd(r) || r < (valF[i] << 1)) return 1;
    }
  return 0;
}

/* returns 2*zetaCO; weight is k + 1/2 */
static long
zeta2CO(GEN faN, GEN valF, long r2, long s2, long k)
{
  if (r2 >= 4) return lamCO(r2, s2, 2) << 1;
  if (r2 == 3) return 6;
  if (condC(faN, valF)) return 4;
  if (odd(k)) return s2 ? 3 : 5; else return s2 ? 5: 3;
}

/* returns 4 times last term in formula */
static long
dim22(long N, long F, long k)
{
  pari_sp av = avma;
  GEN vF, faN = myfactoru(N), P = gel(faN, 1), E = gel(faN, 2);
  long i, D, l = lg(P);
  vF = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) vF[i] = u_lval(F, P[i]);
  D = zeta2CO(faN, vF, E[1], vF[1], k);
  for (i = 2; i < l; i++) D *= lamCO(E[i], vF[i], P[i]);
  return gc_long(av,D);
}

/* PSI not necessarily primitive, of conductor F */
static int
charistotallyeven(GEN PSI, long F)
{
  pari_sp av = avma;
  GEN P = gel(myfactoru(F), 1);
  GEN G = gel(PSI,1), psi = gel(PSI,2);
  long i;
  for (i = 1; i < lg(P); i++)
  {
    GEN psip = znchardecompose(G, psi, utoipos(P[i]));
    if (zncharisodd(G, psip)) return gc_bool(av,0);
  }
  return gc_bool(av,1);
}

static GEN
get_PSI(GEN CHI, long t)
{
  long r = t & 3L, t2 = (r == 2 || r == 3) ? t << 2 : t;
  return mfcharmul_i(CHI, induce(gel(CHI,1), utoipos(t2)));
}
/* space = mf_CUSP, mf_EISEN or mf_FULL, weight k + 1/2 */
static long
mf2dimwt12(long N, GEN CHI, long space)
{
  pari_sp av = avma;
  GEN D = mydivisorsu(N >> 2);
  long i, l = lg(D), dim3 = 0, dim4 = 0;

  CHI = induceN(N, CHI);
  for (i = 1; i < l; i++)
  {
    long rp, t = D[i], Mt = D[l-i];
    GEN PSI = get_PSI(CHI,t);
    rp = mfcharconductor(PSI);
    if (Mt % (rp*rp) == 0) { dim4++; if (charistotallyeven(PSI,rp)) dim3++; }
  }
  set_avma(av);
  switch (space)
  {
    case mf_CUSP: return dim4 - dim3;
    case mf_EISEN:return dim3;
    case mf_FULL: return dim4;
  }
  return 0; /*LCOV_EXCL_LINE*/
}

static long
mf2dimwt32(long N, GEN CHI, long F, long space)
{
  long D;
  switch(space)
  {
    case mf_CUSP: D = mypsiu(N) - 6*dim22(N, F, 1);
      if (D%24) pari_err_BUG("mfdim");
      return D/24 + mf2dimwt12(N, CHI, 4);
    case mf_FULL: D = mypsiu(N) + 6*dim22(N, F, 0);
      if (D%24) pari_err_BUG("mfdim");
      return D/24 + mf2dimwt12(N, CHI, 1);
    case mf_EISEN: D = dim22(N, F, 0) + dim22(N, F, 1);
      if (D & 3L) pari_err_BUG("mfdim");
      return (D >> 2) - mf2dimwt12(N, CHI, 3);
  }
  return 0; /*LCOV_EXCL_LINE*/
}

/* F = conductor(CHI), weight k = r+1/2 */
static long
checkmf2(long N, long r, GEN CHI, long F, long space)
{
  switch(space)
  {
    case mf_FULL: case mf_CUSP: case mf_EISEN: break;
    case mf_NEW: case mf_OLD:
      pari_err_TYPE("half-integral weight [new/old spaces]", utoi(space));
    default:
      pari_err_TYPE("half-integral weight [incorrect space]",utoi(space));
  }
  if (N & 3L)
    pari_err_DOMAIN("half-integral weight", "N % 4", "!=", gen_0, stoi(N));
  return r >= 0 && mfcharparity(CHI) == 1 && N % F == 0;
}

/* weight k = r + 1/2 */
static long
mf2dim_Nkchi(long N, long r, GEN CHI, ulong space)
{
  long D, D2, F = mfcharconductor(CHI);
  if (!checkmf2(N, r, CHI, F, space)) return 0;
  if (r == 0) return mf2dimwt12(N, CHI, space);
  if (r == 1) return mf2dimwt32(N, CHI, F, space);
  if (space == mf_EISEN)
  {
    D = dim22(N, F, r) + dim22(N, F, 1-r);
    if (D & 3L) pari_err_BUG("mfdim");
    return D >> 2;
  }
  D2 = space == mf_FULL? dim22(N, F, 1-r): -dim22(N, F, r);
  D = (2*r-1)*mypsiu(N) + 6*D2;
  if (D%24) pari_err_BUG("mfdim");
  return D/24;
}

/* weight k=r+1/2 */
static GEN
mf2init_Nkchi(long N, long r, GEN CHI, long space, long flraw)
{
  GEN CHI1, Minv, Minvmat, B, M, gk = gaddsg(r,ghalf);
  GEN mf1 = mkvec4(utoi(N),gk,CHI,utoi(space));
  long L;
  if (!checkmf2(N, r, CHI, mfcharconductor(CHI), space)) return mfEMPTY(mf1);
  if (space==mf_EISEN) pari_err_IMPL("half-integral weight Eisenstein space");
  L = mfsturmNgk(N, gk) + 1;
  B = mf2basis(N, r, CHI, &CHI1, space);
  M = mflineardivtomat(N,B,L); /* defined modulo T = charpol(CHI) */
  if (flraw) M = mkvec3(gen_0,gen_0,M);
  else
  {
    long o1 = mfcharorder(CHI1), o = mfcharorder(CHI);
    M = mfcleanCHI(M, CHI, 0);
    Minv = gel(M,2);
    Minvmat = RgM_Minv_mul(NULL, Minv); /* mod T */
    if (o1 != o)
    {
      GEN tr = Qab_trace_init(o, o1, mfcharpol(CHI), mfcharpol(CHI1));
      Minvmat = QabM_tracerel(tr, 0, Minvmat);
    }
    /* Minvmat mod T1 = charpol(CHI1) */
    B = vecmflineardiv_linear(B, Minvmat);
    gel(M,3) = RgM_Minv_mul(gel(M,3), Minv);
    gel(M,2) = mkMinv(matid(lg(B)-1), NULL,NULL,NULL);
  }
  return mkmf(mf1, cgetg(1,t_VEC), B, gen_0, M);
}

/**************************************************************************/
/*                          Kohnen + space                                */
/**************************************************************************/

static GEN
mfkohnenbasis_i(GEN mf, GEN CHI, long eps, long sb)
{
  GEN M = mfcoefs_mf(mf, sb, 1), p, P;
  long c, i, n = mfcharorder(CHI), l = sb + 2;
  p = cgetg(l, t_VECSMALL);
  /* keep the a_n, n = (2 or 2+eps) mod 4 */
  for (i = 3, c = 1; i < l; i+=4) p[c++] = i;
  for (i = 3+eps;    i < l; i+=4) p[c++] = i;
  P = n <= 2? NULL: mfcharpol(CHI);
  setlg(p, c);
  return QabM_ker(rowpermute(M, p), P, n);
}
GEN
mfkohnenbasis(GEN mf)
{
  pari_sp av = avma;
  GEN gk, CHI, CHIP, K;
  long N4, r, eps, sb;
  mf = checkMF(mf);
  if (MF_get_space(mf) != mf_CUSP)
    pari_err_TYPE("mfkohnenbasis [not a cuspidal space", mf);
  if (!MF_get_dim(mf)) return cgetg(1, t_MAT);
  N4 = MF_get_N(mf) >> 2; gk = MF_get_gk(mf); CHI = MF_get_CHI(mf);
  if (typ(gk) == t_INT) pari_err_TYPE("mfkohnenbasis", gk);
  r = MF_get_r(mf);
  CHIP = mfcharchiliftprim(CHI, N4);
  eps = CHIP==CHI? 1: -1;
  if (odd(r)) eps = -eps;
  if (uissquarefree(N4))
  {
    long d = mfdim_Nkchi(N4, 2*r, mfcharpow(CHI, gen_2), mf_CUSP);
    sb = mfsturmNgk(N4 << 2, gk) + 1;
    K = mfkohnenbasis_i(mf, CHIP, eps, sb);
    if (lg(K) - 1 == d) return gerepilecopy(av, K);
  }
  sb = mfsturmNgk(N4 << 4, gk) + 1;
  K = mfkohnenbasis_i(mf, CHIP, eps, sb);
  return gerepilecopy(av, K);
}

static GEN
get_Shimura(GEN mf, GEN CHI, GEN vB, long t)
{
  long N = MF_get_N(mf), r = MF_get_k(mf) >> 1;
  long i, d = MF_get_dim(mf), sb = mfsturm_mf(mf);
  GEN a = cgetg(d+1, t_MAT);
  for (i = 1; i <= d; i++)
  {
    pari_sp av = avma;
    GEN f = c_deflate(sb*sb, t, gel(vB,i));
    f = mftobasis_i(mf, RgV_shimura(f, sb, t, N, r, CHI));
    gel(a,i) = gerepileupto(av, f);
  }
  return a;
}
static long
QabM_rank(GEN M, GEN P, long n)
{
  GEN z = QabM_indexrank(M, P, n);
  return lg(gel(z,2))-1;
}
/* discard T[*i] */
static void
discard_Ti(GEN T, long *i, long *lt)
{
  long j, l = *lt-1;
  for (j = *i; j < l; j++) T[j] = T[j+1];
  (*i)--; *lt = l;
}
/* return [mf3, bijection, mfkohnenbasis, codeshi] */
static GEN
mfkohnenbijection_i(GEN mf)
{
  GEN CHI = MF_get_CHI(mf), K = mfkohnenbasis(mf);
  GEN mres, dMi, Mi, M, C, vB, mf3, SHI, T, P;
  long N4 = MF_get_N(mf)>>2, r = MF_get_r(mf), dK = lg(K) - 1;
  long i, c, n, oldr, lt, ltold, sb3, t, limt;
  const long MAXlt = 100;

  mf3 = mfinit_Nkchi(N4, r<<1, mfcharpow(CHI,gen_2), mf_CUSP, 0);
  if (MF_get_dim(mf3) != dK)
    pari_err_BUG("mfkohnenbijection [different dimensions]");
  if (!dK) return mkvec4(mf3, cgetg(1, t_MAT), K, cgetg(1, t_VEC));
  CHI = mfcharchiliftprim(CHI, N4);
  n = mfcharorder(CHI);
  P = n<=2? NULL: mfcharpol(CHI);
  SHI = cgetg(MAXlt, t_COL);
  T = cgetg(MAXlt, t_VECSMALL);
  sb3 = mfsturm_mf(mf3);
  limt = 6; oldr = 0; vB = C = M = NULL;
  for (t = lt = ltold = 1; lt < MAXlt; t++)
  {
    pari_sp av;
    if (!uissquarefree(t)) continue;
    T[lt++] = t; if (t <= limt) continue;
    av = avma;
    if (vB) gunclone(vB);
    /* could improve the rest but 99% of running time is spent here */
    vB = gclone( RgM_mul(mfcoefs_mf(mf, t*sb3*sb3, 1), K) );
    set_avma(av);
    for (i = ltold; i < lt; i++)
    {
      pari_sp av;
      long r;
      M = get_Shimura(mf3, CHI, vB, T[i]);
      r = QabM_rank(M, P, n); if (!r) { discard_Ti(T, &i, &lt); continue; }
      gel(SHI, i) = M; setlg(SHI, i+1);
      if (r >= dK) { C = vecsmall_ei(dK, i); goto DONE; }
      if (i == 1) { oldr = r; continue; }
      av = avma; M = shallowmatconcat(SHI);
      r = QabM_rank(M, P, n); /* >= rank(sum C[j] SHI[j]), probably sharp */
      if (r >= dK)
      {
        M = RgV_sum(SHI);
        if (QabM_rank(M, P, n) >= dK) { C = const_vecsmall(dK, 1); goto DONE; }
        C = random_Flv(dK, 16);
        M = RgV_zc_mul(SHI, C);
        if (QabM_rank(M, P, n) >= dK) goto DONE;
      }
      else if (r == oldr) discard_Ti(T, &i, &lt);
      oldr = r; set_avma(av);
    }
    limt *= 2; ltold = lt;
  }
  pari_err_BUG("mfkohnenbijection");
DONE:
  gunclone(vB); lt = lg(SHI);
  Mi = QabM_pseudoinv(M,P,n, NULL,&dMi); Mi = RgM_Rg_div(Mi,dMi);
  mres = cgetg(lt, t_VEC);
  for (i = c = 1; i < lt; i++)
    if (C[i]) gel(mres,c++) = mkvec2s(T[i], C[i]);
  setlg(mres,c); return mkvec4(mf3, Mi, K, mres);
}
GEN
mfkohnenbijection(GEN mf)
{
  pari_sp av = avma;
  long N;
  mf = checkMF(mf); N = MF_get_N(mf);
  if (!uissquarefree(N >> 2))
    pari_err_TYPE("mfkohnenbijection [N/4 not squarefree]", utoi(N));
  if (MF_get_space(mf) != mf_CUSP || MF_get_r(mf) == 0 || !mfshimura_space_cusp(mf))
    pari_err_TYPE("mfkohnenbijection [incorrect mf for Kohnen]", mf);
  return gerepilecopy(av, mfkohnenbijection_i(mf));
}

static int
checkbij_i(GEN b)
{
  return typ(b) == t_VEC && lg(b) == 5 && checkMF_i(gel(b,1))
         && typ(gel(b,2)) == t_MAT
         && typ(gel(b,3)) == t_MAT
         && typ(gel(b,4)) == t_VEC;
}

/* bij is the output of mfkohnenbijection */
GEN
mfkohneneigenbasis(GEN mf, GEN bij)
{
  pari_sp av = avma;
  GEN mf3, mf30, B, KM, M, k;
  long r, i, l, N4;
  mf = checkMF(mf);
  if (!checkbij_i(bij))
    pari_err_TYPE("mfkohneneigenbasis [bijection]", bij);
  if (MF_get_space(mf) != mf_CUSP)
    pari_err_TYPE("mfkohneneigenbasis [not a cuspidal space]", mf);
  if (!MF_get_dim(mf))
    retmkvec3(cgetg(1, t_VEC), cgetg(1, t_VEC), cgetg(1, t_VEC));
  N4 = MF_get_N(mf) >> 2; k = MF_get_gk(mf);
  if (typ(k) == t_INT) pari_err_TYPE("mfkohneneigenbasis", k);
  if (!uissquarefree(N4))
    pari_err_TYPE("mfkohneneigenbasis [N not squarefree]", utoipos(N4));
  r = MF_get_r(mf);
  KM = RgM_mul(gel(bij,3), gel(bij,2));
  mf3 = gel(bij,1);
  mf30 = mfinit_Nkchi(N4, 2*r, MF_get_CHI(mf3), mf_NEW, 0);
  B = mfcoefs_mf(mf30, mfsturm_mf(mf3), 1); l = lg(B);
  M = cgetg(l, t_MAT);
  for (i=1; i<l; i++) gel(M,i) = RgM_RgC_mul(KM, mftobasis_i(mf3, gel(B,i)));
  return gerepilecopy(av, mkvec3(mf30, M, RgM_mul(M, MF_get_newforms(mf30))));
}
/*************************** End Kohnen ************************************/
/***************************************************************************/

static GEN desc(GEN F);
static GEN
desc_mfeisen(GEN F)
{
  GEN R, gk = mf_get_gk(F);
  if (typ(gk) == t_FRAC)
    R = gsprintf("H_{%Ps}", gk);
  else
  {
    GEN vchi = gel(F, 2), CHI = mfchisimpl(gel(vchi, 3));
    long k = itou(gk);
    if (lg(vchi) < 5) R = gsprintf("F_%ld(%Ps)", k, CHI);
    else
    {
      GEN CHI2 = mfchisimpl(gel(vchi, 4));
      R = gsprintf("F_%ld(%Ps, %Ps)", k, CHI, CHI2);
    }
  }
  return R;
}
static GEN
desc_hecke(GEN F)
{
  long n, N;
  GEN D = gel(F,2);
  if (typ(D) == t_VECSMALL) { N = D[3]; n = D[1]; }
  else { GEN nN = gel(D,2); n = nN[1]; N = nN[2]; } /* half integer */
  return gsprintf("T_%ld(%ld)(%Ps)", N, n, desc(gel(F,3)));
}
static GEN
desc_linear(GEN FLD, GEN dL)
{
  GEN F = gel(FLD,2), L = gel(FLD,3), R = strtoGENstr("LIN([");
  long n = lg(F) - 1, i;
  for (i = 1; i <= n; i++)
  {
    R = shallowconcat(R, desc(gel(F,i))); if (i == n) break;
    R = shallowconcat(R, strtoGENstr(", "));
  }
  return shallowconcat(R, gsprintf("], %Ps)", gdiv(L, dL)));
}
static GEN
desc_dihedral(GEN F)
{
  GEN bnr = gel(F,2), D = nf_get_disc(bnr_get_nf(bnr)), f = bnr_get_mod(bnr);
  GEN cyc = bnr_get_cyc(bnr);
  GEN w = gel(F,3), chin = zv_to_ZV(gel(w,2)), o = utoi(gel(w,1)[1]);
  GEN chi = char_denormalize(cyc, o, chin);
  if (lg(gel(f,2)) == 1) f = gel(f,1);
  return gsprintf("DIH(%Ps, %Ps, %Ps, %Ps)",D,f,cyc,chi);
}

static void
unpack0(GEN *U)
{ if (U) *U = mkvec2(cgetg(1, t_VEC), cgetg(1, t_VEC)); }
static void
unpack2(GEN F, GEN *U)
{ if (U) *U = mkvec2(mkvec2(gel(F,2), gel(F,3)), cgetg(1, t_VEC)); }
static void
unpack23(GEN F, GEN *U)
{ if (U) *U = mkvec2(mkvec(gel(F,2)), mkvec(gel(F,3))); }
static GEN
desc_i(GEN F, GEN *U)
{
  switch(mf_get_type(F))
  {
    case t_MF_CONST: unpack0(U); return gsprintf("CONST(%Ps)", gel(F,2));
    case t_MF_EISEN: unpack0(U); return desc_mfeisen(F);
    case t_MF_Ek: unpack0(U); return gsprintf("E_%ld", mf_get_k(F));
    case t_MF_DELTA: unpack0(U); return gsprintf("DELTA");
    case t_MF_THETA: unpack0(U);
      return gsprintf("THETA(%Ps)", mfchisimpl(gel(F,2)));
    case t_MF_ETAQUO: unpack0(U);
      return gsprintf("ETAQUO(%Ps, %Ps)", gel(F,2), gel(F,3));
    case t_MF_ELL: unpack0(U);
      return gsprintf("ELL(%Ps)", vecslice(gel(F,2), 1, 5));
    case t_MF_TRACE: unpack0(U); return gsprintf("TR(%Ps)", mfparams(F));
    case t_MF_NEWTRACE: unpack0(U); return gsprintf("TR^new(%Ps)", mfparams(F));
    case t_MF_DIHEDRAL: unpack0(U); return desc_dihedral(F);
    case t_MF_MUL: unpack2(F, U);
      return gsprintf("MUL(%Ps, %Ps)", desc(gel(F,2)), desc(gel(F,3)));
    case t_MF_DIV: unpack2(F, U);
      return gsprintf("DIV(%Ps, %Ps)", desc(gel(F,2)), desc(gel(F,3)));
    case t_MF_POW: unpack23(F, U);
      return gsprintf("POW(%Ps, %ld)", desc(gel(F,2)), itos(gel(F,3)));
    case t_MF_SHIFT: unpack23(F, U);
      return gsprintf("SHIFT(%Ps, %ld)", desc(gel(F,2)), itos(gel(F,3)));
    case t_MF_DERIV: unpack23(F, U);
      return gsprintf("DER^%ld(%Ps)", itos(gel(F,3)), desc(gel(F,2)));
    case t_MF_DERIVE2: unpack23(F, U);
      return gsprintf("DERE2^%ld(%Ps)", itos(gel(F,3)), desc(gel(F,2)));
    case t_MF_TWIST: unpack23(F, U);
      return gsprintf("TWIST(%Ps, %Ps)", desc(gel(F,2)), gel(F,3));
    case t_MF_BD: unpack23(F, U);
      return gsprintf("B(%ld)(%Ps)", itou(gel(F,3)), desc(gel(F,2)));
    case t_MF_BRACKET:
      if (U) *U = mkvec2(mkvec2(gel(F,2), gel(F,3)), mkvec(gel(F,4)));
      return gsprintf("MULRC_%ld(%Ps, %Ps)", itos(gel(F,4)), desc(gel(F,2)), desc(gel(F,3)));
    case t_MF_LINEAR_BHN:
    case t_MF_LINEAR:
      if (U) *U = mkvec2(gel(F,2), mkvec(gdiv(gel(F,3), gel(F,4))));
      return desc_linear(F,gel(F,4));
    case t_MF_HECKE:
      if (U) *U = mkvec2(mkvec(gel(F,3)), mkvec(stoi(gel(F,2)[1])));
      return desc_hecke(F);
    default: pari_err_TYPE("mfdescribe",F);
    return NULL;/*LCOV_EXCL_LINE*/
  }
}
static GEN
desc(GEN F) { return desc_i(F, NULL); }
GEN
mfdescribe(GEN F, GEN *U)
{
  pari_sp av = avma;
  GEN mf;
  if ((mf = checkMF_i(F)))
  {
    const char *f = NULL;
    switch (MF_get_space(mf))
    {
      case mf_NEW:  f = "S_%Ps^new(G_0(%ld, %Ps))"; break;
      case mf_CUSP: f = "S_%Ps(G_0(%ld, %Ps))"; break;
      case mf_OLD:  f = "S_%Ps^old(G_0(%ld, %Ps))"; break;
      case mf_EISEN:f = "E_%Ps(G_0(%ld, %Ps))"; break;
      case mf_FULL: f = "M_%Ps(G_0(%ld, %Ps))"; break;
    }
    if (U) *U = cgetg(1, t_VEC);
    return gsprintf(f, MF_get_gk(mf), MF_get_N(mf), mfchisimpl(MF_get_CHI(mf)));
  }
  if (!checkmf_i(F)) pari_err_TYPE("mfdescribe", F);
  F = desc_i(F, U); return gc_all(av, U ? 2: 1, &F, U);
}

/***********************************************************************/
/*               Eisenstein series H_r of weight r+1/2                 */
/***********************************************************************/
/* radical(u_ppo(g,q)) */
static long
u_pporad(long g, long q)
{
  GEN F = myfactoru(g), P = gel(F,1);
  long i, l, n;
  if (q == 1) return zv_prod(P);
  l = lg(P);
  for (i = n = 1; i < l; i++)
  {
    long p = P[i];
    if (q % p) n *= p;
  }
  return n;
}
static void
c_F2TH4(long n, GEN *pF2, GEN *pTH4)
{
  GEN v = mfcoefs_i(mfEk(2), n, 1), v2 = bdexpand(v,2), v4 = bdexpand(v,4);
  GEN F2 = gdivgs(ZC_add(ZC_sub(v, ZC_z_mul(v2,3)), ZC_z_mul(v4,2)), -24);
  GEN TH4 = gdivgs(ZC_sub(v, ZC_z_mul(v4,4)), -3);
  settyp(F2,t_VEC); *pF2 = F2;
  settyp(TH4,t_VEC);*pTH4= TH4;
}
/* r > 0, N >= 0 */
static GEN
mfEHcoef(long r, long N)
{
  long D0, f, i, l, s;
  GEN S, Df;

  if (r == 1) return hclassno(utoi(N));
  if (N == 0) return gdivgs(bernfrac(2*r), -2*r);
  s = N & 3L;
  if (odd(r))
  {
    if (s == 2 || s == 1) return gen_0;
    D0 = mycoredisc2neg(N,&f);
  }
  else
  {
    if (s == 2 || s == 3) return gen_0;
    D0 = mycoredisc2pos(N,&f);
  }
  Df = mydivisorsu(u_pporad(f, D0)); l = lg(Df);
  S = gen_0;
  for (i = 1; i < l; i++)
  {
    long d = Df[i], s = mymoebiusu(d)*kross(D0, d); /* != 0 */
    GEN c = gmul(powuu(d, r-1), mysumdivku(f/d, 2*r-1));
    S = s > 0? addii(S, c): subii(S, c);
  }
  return gmul(lfunquadneg_naive(D0, r), S);
}
static GEN
mfEHmat(long lim, long r)
{
  long j, l, d = r/2;
  GEN f2, th4, th3, v, vth4, vf2;
  c_F2TH4(lim, &f2, &th4);
  f2 =  RgV_to_ser(f2, 0, lim+3);
  th4 = RgV_to_ser(th4, 0, lim+3);
  th3 = RgV_to_ser(c_theta(lim, 1, mfchartrivial()), 0, lim+3);
  if (odd(r)) th3 = gpowgs(th3, 3);
  vth4 = gpowers(th4, d);
  vf2 = gpowers0(f2, d, th3); /* th3 f2^j */
  l = d+2; v = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
    gel(v, j) = ser2rfrac_i(gmul(gel(vth4, l-j), gel(vf2, j)));
  return RgXV_to_RgM(v, lim);
}
static GEN
Hfind(long r, GEN *pden)
{
  long lim = (r/2)+3, i;
  GEN res, M, B;

  if (r <= 0) pari_err_DOMAIN("mfEH", "r", "<=", gen_0, stoi(r));
  M = mfEHmat(lim, r);
  B = cgetg(lim+1, t_COL);
  for (i = 1; i <= lim; i++) gel(B, i) = mfEHcoef(r, i-1);
  res = QM_gauss(M, B);
  if (lg(res) == 1) pari_err_BUG("mfEH");
  return Q_remove_denom(res,pden);
}
GEN
mfEH(GEN gk)
{
  pari_sp av = avma;
  GEN v, d, NK, gr = gsub(gk, ghalf);
  long r;
  if (typ(gr) != t_INT) pari_err_TYPE("mfEH", gk);
  r = itos(gr);
  switch (r)
  {
    case 1: v=cgetg(1,t_VEC); d=gen_1; break;
    case 2: v=mkvec2s(1,-20); d=utoipos(120); break;
    case 3: v=mkvec2s(-1,14); d=utoipos(252); break;
    case 4: v=mkvec3s(1,-16,16); d=utoipos(240); break;
    case 5: v=mkvec3s(-1,22,-88); d=utoipos(132); break;
    case 6: v=mkvec4s(691,-18096,110136,-4160); d=utoipos(32760); break;
    case 7: v=mkvec4s(-1,30,-240,224); d=utoipos(12); break;
    default: v = Hfind(r, &d); break;
  }
  NK = mkgNK(utoipos(4), gaddgs(ghalf,r), mfchartrivial(), pol_x(1));
  return gerepilecopy(av, tag(t_MF_EISEN, NK, mkvec2(v,d)));
}

/**********************************************************/
/*             T(f^2) for half-integral weight            */
/**********************************************************/

/* T_p^2 V, p2 = p^2, c1 = chi(p) (-1/p)^r p^(r-1), c2 = chi(p^2)*p^(2r-1) */
static GEN
tp2apply(GEN V, long p, long p2, GEN c1, GEN c2)
{
  long lw = (lg(V) - 2)/p2 + 1, m, n;
  GEN a0 = gel(V,1), W = cgetg(lw + 1, t_VEC);

  gel(W,1) = gequal0(a0)? gen_0: gmul(a0, gaddsg(1, c2));
  for (n = 1; n < lw; n++)
  {
    GEN c = gel(V, p2*n + 1);
    if (n%p) c = gadd(c, gmulsg(kross(n,p), gmul(gel(V,n+1), c1)));
    gel(W, n+1) = c; /* a(p^2*n) + c1 * (n/p) a(n) */
  }
  for (m = 1, n = p2; n < lw; m++, n += p2)
    gel(W, n+1) = gadd(gel(W,n+1), gmul(gel(V,m+1), c2));
  return W;
}

/* T_{p^{2e}} V; can derecursify [Purkait, Hecke operators in half-integral
 * weight, Prop 4.3], not worth it */
static GEN
tp2eapply(GEN V, long p, long p2, long e, GEN q, GEN c1, GEN c2)
{
  GEN V4 = NULL;
  if (e > 1)
  {
    V4 = vecslice(V, 1, (lg(V) - 2)/(p2*p2) + 1);
    V = tp2eapply(V, p, p2, e-1, q, c1, c2);
  }
  V = tp2apply(V, p, p2, c1, c2);
  if (e > 1)
    V = gsub(V, (e == 2)? gmul(q, V4)
                        : gmul(c2, tp2eapply(V4, p, p2, e-2, q, c1, c2)));
  return V;
}
/* weight k = r+1/2 */
static GEN
RgV_heckef2(long n, long d, GEN V, GEN F, GEN DATA)
{
  GEN CHI = mf_get_CHI(F), fa = gel(DATA,1), P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), r = mf_get_r(F), s4 = odd(r)? -4: 4, k2m2 = (r<<1)-1;
  if (typ(V) == t_COL) V = shallowtrans(V);
  for (i = 1; i < l; i++)
  { /* p does not divide N */
    long p = P[i], e = E[i], p2 = p*p;
    GEN c1, c2, a, b, q = NULL, C = mfchareval(CHI,p), C2 = gsqr(C);
    a = r? powuu(p,r-1): mkfrac(gen_1,utoipos(p)); /* p^(r-1) = p^(k-3/2) */
    b = r? mulii(powuu(p,r), a): a; /* p^(2r-1) = p^(2k-2) */
    c1 = gmul(C, gmulsg(kross(s4,p),a));
    c2 = gmul(C2, b);
    if (e > 1)
    {
      q = r? powuu(p,k2m2): a;
      if (e == 2) q = gmul(q, uutoQ(p+1,p)); /* special case T_{p^4} */
      q = gmul(C2, q); /* chi(p^2) [ p^(2k-2) or (p+1)p^(2k-3) ] */
    }
    V = tp2eapply(V, p, p2, e, q, c1, c2);
  }
  return c_deflate(n, d, V);
}

static GEN
GL2toSL2(GEN g, GEN *abd)
{
  GEN A, B, C, D, u, v, a, b, d, q;
  g = Q_primpart(g);
  if (!check_M2Z(g)) pari_err_TYPE("GL2toSL2", g);
  A = gcoeff(g,1,1); B = gcoeff(g,1,2);
  C = gcoeff(g,2,1); D = gcoeff(g,2,2);
  a = bezout(A, C, &u, &v);
  if (!equali1(a)) { A = diviiexact(A,a); C = diviiexact(C,a); }
  d = subii(mulii(A,D), mulii(B,C));
  if (signe(d) <= 0) pari_err_TYPE("GL2toSL2",g);
  q = dvmdii(addii(mulii(u,B), mulii(v,D)), d, &b);
  *abd = (equali1(a) && equali1(d))? NULL: mkvec3(a, b, d);
  return mkmat22(A, subii(mulii(q,A), v), C, addii(mulii(q,C), u));
}

static GEN
Rg_approx(GEN t, long bit)
{
  GEN a = real_i(t), b = imag_i(t);
  long e1 = gexpo(a), e2 = gexpo(b);
  if (e2 < -bit) { t = e1 < -bit? gen_0: a; }
  else if (e1 < -bit) t = gmul(b, gen_I());
  return t;
}
static GEN
RgV_approx(GEN v, long bit)
{
  long i, l = lg(v);
  GEN w = cgetg_copy(v, &l);
  for (i = 1; i < lg(v); i++) gel(w,i) = Rg_approx(gel(v,i), bit);
  return w;
}
/* m != 2 (mod 4), D t_INT; V has "denominator" D, recognize in Q(zeta_m) */
static GEN
bestapprnf2(GEN V, long m, GEN D, long prec)
{
  long i, j, f, vt = fetch_user_var("t"), bit = prec2nbits_mul(prec, 0.8);
  GEN Tinit, Vl, H, Pf, P = polcyclo(m, vt);

  V = liftpol_shallow(V);
  V = gmul(RgV_approx(V, bit), D);
  V = bestapprnf(V, P, NULL, prec);
  Vl = liftpol_shallow(V);
  H = coprimes_zv(m);
  for (i = 2; i < m; i++)
  {
    if (H[i] != 1) continue;
    if (!gequal(Vl, vecGalois(Vl, i, P, m))) H[i] = 0;
    else for (j = i; j < m; j *= i) H[i] = 3;
  }
  f = znstar_conductor_bits(Flv_to_F2v(H));
  if (f == 1) return gdiv(V, D);
  if (f == m) return gmodulo(gdiv(V, D), P);
  Pf = polcyclo(f, vt);
  Tinit = Qab_trace_init(m, f, P, Pf);
  Vl = QabV_tracerel(Tinit, 0, Vl);
  return gmodulo(gdiv(Vl, muliu(D, degpol(P)/degpol(Pf))), Pf);
}

/* f | ga expansion; [f, mf_eisendec(f)]~ allowed */
GEN
mfslashexpansion(GEN mf, GEN f, GEN ga, long n, long flrat, GEN *params, long prec)
{
  pari_sp av = avma;
  GEN a, b, d, res, al, V, M, ad, abd, gk, A, awd = NULL;
  long i, w;

  mf = checkMF(mf);
  gk = MF_get_gk(mf);
  M = GL2toSL2(ga, &abd);
  if (abd) { a = gel(abd,1); b = gel(abd,2); d = gel(abd,3); }
  else { a = d = gen_1; b = gen_0; }
  ad = gdiv(a,d);
  res = mfgaexpansion(mf, f, M, n, prec);
  al = gel(res,1);
  w = itou(gel(res,2));
  V = gel(res,3);
  if (flrat)
  {
    GEN CHI = MF_get_CHI(mf);
    long N = MF_get_N(mf), F = mfcharconductor(CHI);
    long ord = mfcharorder(CHI), k, deg;
    long B = umodiu(gcoeff(M,1,2), N);
    long C = umodiu(gcoeff(M,2,1), N);
    long D = umodiu(gcoeff(M,2,2), N);
    long CD = (C * D) % N, BC = (B * C) % F;
    GEN CV, t;
    /* weight of f * Theta in 1/2-integral weight */
    k = typ(gk) == t_INT? (long) itou(gk): MF_get_r(mf)+1;
    CV = odd(k) ? powuu(N, k - 1) : powuu(N, k >> 1);
    deg = ulcm(ulcm(ord, N/ugcd(N,CD)), F/ugcd(F,BC));
    if ((deg & 3) == 2) deg >>= 1;
    if (typ(gk) != t_INT && odd(deg) && mfthetaI(C,D)) deg <<= 2;
    V = bestapprnf2(V, deg, CV, prec);
    if (abd && !signe(b))
    { /* can [a,0; 0,d] be simplified to id ? */
      long nk, dk; Qtoss(gk, &nk, &dk);
      if (ispower(ad, utoipos(2*dk), &t)) /* t^(2*dk) = a/d or t = NULL */
      {
        V = RgV_Rg_mul(V, powiu(t,nk));
        awd = gdiv(a, muliu(d,w));
      }
    }
  }
  else if (abd)
  { /* ga = M * [a,b;0,d] * rational, F := f | M = q^al * \sum V[j] q^(j/w) */
    GEN u, t = NULL, wd = muliu(d,w);
    /* a > 0, 0 <= b < d; f | ga = (a/d)^(k/2) * F(tau + b/d) */
    if (signe(b))
    {
      long ns, ds;
      GEN z;
      Qtoss(gdiv(b, wd), &ns, &ds); z = rootsof1powinit(ns, ds, prec);
      for (i = 1; i <= n+1; i++) gel(V,i) = gmul(gel(V,i), rootsof1pow(z, i-1));
      if (!gequal0(al)) t = gexp(gmul(PiI2(prec), gmul(al, gdiv(b,d))), prec);
    }
    awd = gdiv(a, wd);
    u = gpow(ad, gmul2n(gk,-1), prec);
    t = t? gmul(t, u): u;
    V = RgV_Rg_mul(V, t);
  }
  if (!awd) A = mkmat22(a, b, gen_0, d);
  else
  { /* rescale and update w from [a,0; 0,d] */
    long ns;
    Qtoss(awd, &ns, &w); /* update w */
    V = bdexpand(V, ns);
    if (!gequal0(al))
    {
      GEN adal = gmul(ad, al), sh = gfloor(adal);
      al = gsub(adal, sh);
      V = RgV_shift(V, sh);
    }
    A = matid(2);
  }
  if (params) *params = mkvec3(al, utoipos(w), A);
  return gc_all(av,params?2:1,&V,params);
}

/**************************************************************/
/*         Alternative method for 1/2-integral weight         */
/**************************************************************/
static GEN
mf2basis(long N, long r, GEN CHI, GEN *pCHI1, long space)
{
  GEN CHI1, CHI2, mf1, mf2, B1, B2, BT, M1, M2, M, M2i, T, Th, v, den;
  long sb, N2, o1, o2, k1 = r + 1;

  if (odd(k1))
  {
    CHI1 = mfcharmul(CHI, get_mfchar(stoi(-4)));
    CHI2 = mfcharmul(CHI, get_mfchar(stoi(-8)));
  }
  else
  {
    CHI1 = CHI;
    CHI2 = mfcharmul(CHI, get_mfchar(utoi(8)));
  }
  mf1 = mfinit_Nkchi(N, k1, CHI1, space, 1);
  if (pCHI1) *pCHI1 = CHI1;
  B1 = MF_get_basis(mf1); if (lg(B1) == 1) return cgetg(1,t_VEC);
  N2 = ulcm(8, N);
  mf2 = mfinit_Nkchi(N2, k1, CHI2, space, 1);
  B2 = MF_get_basis(mf2); if (lg(B2) == 1) return cgetg(1,t_VEC);
  sb = mfsturmNgk(N2, gaddsg(k1, ghalf));
  M1 = mfcoefs_mf(mf1, sb, 1);
  M2 = mfcoefs_mf(mf2, sb, 1);
  Th = mfTheta(NULL);
  BT = mfcoefs_i(Th, sb, 1);
  M1 = mfmatsermul(M1, RgV_to_RgX(expandbd(BT,2),0));
  M2 = mfmatsermul(M2, RgV_to_RgX(BT,0));
  o1= mfcharorder(CHI1);
  T = (o1 <= 2)? NULL: mfcharpol(CHI1);
  if (o1 > 2) M1 = liftpol_shallow(M1);
  o2= mfcharorder(CHI2);
  if (T)
  {
    if (o2 == o1) M2 = liftpol_shallow(M2);
    else
    {
      GEN tr = Qab_trace_init(o2, o1, mfcharpol(CHI2), mfcharpol(CHI1));
      M2 = QabM_tracerel(tr, 0, M2);
    }
  }
  /* now everything is defined mod T = mfcharpol(CHI1) */
  M2i = QabM_pseudoinv_i(M2, T, o1, &v, &den);
  M = RgM_mul(M2i, rowpermute(M1, gel(v,1)));
  M = RgM_mul(M2, M);
  M1 = RgM_Rg_mul(M1, den);
  M = RgM_sub(M1, M); if (T) M = RgXQM_red(M, T);
  return vecmflineardiv0(B1, QabM_ker(M, T, o1), Th);
}

/*******************************************************************/
/*                         Integration                             */
/*******************************************************************/
static GEN
vanembed(GEN F, GEN v, long prec)
{
  GEN CHI = mf_get_CHI(F);
  long o = mfcharorder(CHI);
  if (o > 2 || degpol(mf_get_field(F)) > 1) v = liftpol_shallow(v);
  if (o > 2) v = gsubst(v, varn(mfcharpol(CHI)), rootsof1u_cx(o, prec));
  return v;
}

static long
mfperiod_prelim_double(double t0, long k, long bitprec)
{
  double nlim, c = 2*M_PI*t0;
  nlim = ceil(bitprec * M_LN2 / c);
  c -= (k - 1)/(2*nlim); if (c < 1) c = 1.;
  nlim += ceil((0.7 + (k-1)/2*log(nlim))/c);
  return (long)nlim;
}
static long
mfperiod_prelim(GEN t0, long k, long bitprec)
{ return mfperiod_prelim_double(gtodouble(t0), k, bitprec); }

/* (-X)^(k-2) * P(-1/X) = (-1)^{k-2} P|_{k-2} S */
static GEN
RgX_act_S(GEN P, long k)
{
  P = RgX_unscale(RgX_recipspec_shallow(P+2, lgpol(P), k-1), gen_m1);
  setvarn(P,0); return P;
}
static int
RgX_act_typ(GEN P, long k)
{
  switch(typ(P))
  {
    case t_RFRAC: return t_RFRAC;
    case t_POL:
      if (varn(P) == 0)
      {
        long d = degpol(P);
        if (d > k-2) return t_RFRAC;
        if (d) return t_POL;
      }
  }
  return 0;
}
static GEN
act_S(GEN P, long k)
{
  GEN X;
  switch(RgX_act_typ(P, k))
  {
    case t_RFRAC:
      X = gneg(pol_x(0));
      return gmul(gsubst(P, 0, ginv(X)), gpowgs(X, k - 2));
    case t_POL: return RgX_act_S(P, k);
  }
  return P;
}

static GEN
AX_B(GEN M)
{ GEN A = gcoeff(M,1,1), B = gcoeff(M,1,2); return deg1pol_shallow(A,B,0); }
static GEN
CX_D(GEN M)
{ GEN C = gcoeff(M,2,1), D = gcoeff(M,2,2); return deg1pol_shallow(C,D,0); }

/* P|_{2-k}M = (CX+D)^{k-2}P((AX+B)/(CX+D)) */
static GEN
RgX_act_gen(GEN P, GEN M, long k)
{
  GEN S = gen_0, PCD, PAB;
  long i;
  PCD = gpowers(CX_D(M), k-2);
  PAB = gpowers(AX_B(M), k-2);
  for (i = 0; i <= k-2; i++)
  {
    GEN t = RgX_coeff(P, i);
    if (!gequal0(t)) S = gadd(S, gmul(t, gmul(gel(PCD, k-i-1), gel(PAB, i+1))));
  }
  return S;
}
static GEN
act_GL2(GEN P, GEN M, long k)
{
  switch(RgX_act_typ(P, k))
  {
    case t_RFRAC:
    {
      GEN AB = AX_B(M), CD = CX_D(M);
      return gmul(gsubst(P, 0, gdiv(AB, CD)), gpowgs(CD, k - 2));
    }
    case t_POL: return RgX_act_gen(P, M, k);
  }
  return P;
}
static GEN
vecact_GL2(GEN x, GEN M, long k)
{ pari_APPLY_same(act_GL2(gel(x,i), M, k)); }

static GEN
normalizeapprox(GEN R, long bit)
{
  GEN X = gen_1, Q;
  long i, l;
  if (is_vec_t(typ(R)))
  {
    Q = cgetg_copy(R, &l);
    for (i = 1; i < l; i++)
    {
      pari_sp av = avma;
      gel(Q,i) = gerepileupto(av, normalizeapprox(gel(R,i), bit));
    }
    return Q;
  }
  if (typ(R) == t_RFRAC && varn(gel(R,2)) == 0) { X = gel(R,2); R = gel(R,1); }
  if (typ(R) != t_POL || varn(R) != 0) return gdiv(R, X);
  Q = cgetg_copy(R, &l); Q[1] = R[1];
  for (i = 2; i < l; ++i) gel(Q,i) = Rg_approx(gel(R,i),bit);
  Q = normalizepol_lg(Q,l); return gdiv(Q, X);
}

/* make sure T is a t_POL in variable 0 */
static GEN
toRgX0(GEN T)
{ return typ(T) == t_POL && varn(T) == 0? T: scalarpol_shallow(T,0); }

/* integrate by summing  nlim+1 terms of van [may be < lg(van)]
 * van can be an expansion with vector coefficients
 * \int_A^oo \sum_n van[n] * q^(n/w + al) * P(z-A) dz, q = e(z) */
static GEN
intAoo(GEN van, long nlim, GEN al, long w, GEN P, GEN A, long k, long prec)
{
  GEN alw, P1, piI2A, q, S, van0;
  long n, vz = varn(gel(P,2));

  if (nlim < 1) nlim = 1;
  alw = gmulsg(w, al);
  P1 = RgX_translate(P, gneg(A));
  piI2A = gmul(PiI2n(1, prec), A);
  q = gexp(gdivgu(piI2A, w), prec);
  S = gen_0;
  for (n = nlim; n >= 1; n--)
  {
    GEN t = gsubst(P1, vz, gdivsg(w, gaddsg(n, alw)));
    S = gadd(gmul(gel(van, n+1), t), gmul(q, S));
  }
  S = gmul(q, S);
  van0 = gel(van, 1);
  if (!gequal0(al))
  {
    S = gadd(S, gmul(gsubst(P1, vz, ginv(al)), van0));
    S = gmul(S, gexp(gmul(piI2A, al), prec));
  }
  else if (!gequal0(van0))
    S = gsub(S, gdivgu(gmul(van0, gpowgs(gsub(pol_x(0), A), k-1)), k-1));
  if (is_vec_t(typ(S)))
  {
    long j, l = lg(S);
    for (j = 1; j < l; j++) gel(S,j) = toRgX0(gel(S,j));
  }
  else
    S = toRgX0(S);
  return gneg(S);
}

/* \sum_{j <= k} X^j * (Y / (2I\pi))^{k+1-j} k! / j! */
static GEN
get_P(long k, long v, long prec)
{
  GEN a, S = cgetg(k + 1, t_POL), u = invr(Pi2n(1, prec+EXTRAPREC64));
  long j, K = k-2;
  S[1] = evalsigne(1)|evalvarn(0); a = u;
  gel(S,K+2) = monomial(mulcxpowIs(a,3), 1, v); /* j = K */
  for(j = K-1; j >= 0; j--)
  {
    a = mulrr(mulru(a,j+1), u);
    gel(S,j+2) = monomial(mulcxpowIs(a,3*(K+1-j)), K+1-j, v);
  }
  return S;
}

static GEN
getw1w2(long N, GEN ga)
{ return mkvecsmall2(mfZC_width(N, gel(ga,1)),
                     mfZC_width(N, gel(ga,2))); }

static GEN
intAoowithvanall(GEN mf, GEN vanall, GEN P, GEN cosets, long bitprec)
{
  GEN vvan = gel(vanall,1), vaw = gel(vanall,2), W1W2, resall;
  long prec = nbits2prec(bitprec), N, k, lco, j;

  N = MF_get_N(mf); k = MF_get_k(mf);
  lco = lg(cosets);
  W1W2 = cgetg(lco, t_VEC); resall = cgetg(lco, t_VEC);
  for (j = 1; j < lco; j++) gel(W1W2,j) = getw1w2(N, gel(cosets, j));
  for (j = 1; j < lco; j++)
  {
    GEN w1w2j = gel(W1W2,j), alj, M, VAN, RES, AR, Q;
    long jq, c, w1, w2, w;
    if (!w1w2j) continue;
    alj = gel(vaw,j);
    w1 = w1w2j[1]; Q = cgetg(lco, t_VECSMALL);
    w2 = w1w2j[2]; M = cgetg(lco, t_COL);
    for (c = 1, jq = j; jq < lco; jq++)
    {
      GEN W = gel(W1W2, jq);
      if (jq == j || (W && gequal(W, w1w2j) && gequal(gel(vaw, jq), alj)))
      {
        Q[c] = jq; gel(W1W2, jq) = NULL;
        gel(M, c) = gel(vvan, jq); c++;
      }
    }
    setlg(M,c); VAN = shallowmatconcat(M);
    AR = mkcomplex(gen_0, sqrtr_abs(divru(utor(w1, prec+EXTRAPREC64), w2)));
    w = itos(gel(alj,2));
    RES = intAoo(VAN, lg(VAN)-2, gel(alj,1),w, P, AR, k, prec);
    for (jq = 1; jq < c; jq++) gel(resall, Q[jq]) = gel(RES, jq);
  }
  return resall;
}

GEN
mftobasisES(GEN mf, GEN F)
{
  GEN v = mftobasis(mf, F, 0);
  long nE = lg(MF_get_E(mf))-1;
  return mkvec2(vecslice(v,1,nE), vecslice(v,nE+1,lg(v)-1));
}

static long
wt1mulcond(GEN F, long D, long space)
{
  GEN E = mfeisenstein_i(1, mfchartrivial(), get_mfchar(stoi(D))), mf;
  F = mfmul(F, E);
  mf = mfinit_Nkchi(mf_get_N(F), mf_get_k(F), mf_get_CHI(F), space, 0);
  return mfconductor(mf, F);
}
static int
wt1newlevel(long N)
{
  GEN P = gel(myfactoru(N),1);
  long l = lg(P), i;
  for (i = 1; i < l; i++)
    if (!wt1empty(N/P[i])) return 0;
  return 1;
}
long
mfconductor(GEN mf, GEN F)
{
  pari_sp av = avma;
  GEN gk;
  long space, N, M;

  mf = checkMF(mf);
  if (!checkmf_i(F)) pari_err_TYPE("mfconductor",F);
  if (mfistrivial(F)) return 1;
  space = MF_get_space(mf);
  if (space == mf_NEW) return mf_get_N(F);
  gk = MF_get_gk(mf);
  if (isint1(gk))
  {
    N = mf_get_N(F);
    if (!wt1newlevel(N))
    {
      long s = space_is_cusp(space)? mf_CUSP: mf_FULL;
      N = ugcd(N, wt1mulcond(F,-3,s));
      if (!wt1newlevel(N)) N = ugcd(N, wt1mulcond(F,-4,s));
    }
    return gc_long(av,N);
  }
  if (typ(gk) != t_INT)
  {
    F = mfmultheta(F);
    mf = obj_checkbuild(mf, MF_MF2INIT, &mf2init); /* mf_FULL */
  }
  N = 1;
  if (space_is_cusp(space))
  {
    F = mftobasis_i(mf, F);
    if (typ(gk) != t_INT) F = vecslice(F, lg(MF_get_E(mf)), lg(F) - 1);
  }
  else
  {
    GEN EF = mftobasisES(mf, F), vE = gel(EF,1), B = MF_get_E(mf);
    long i, l = lg(B);
    for (i = 1; i < l; i++)
      if (!gequal0(gel(vE,i))) N = ulcm(N, mf_get_N(gel(B, i)));
    F = gel(EF,2);
  }
  (void)mftonew_i(mf, F, &M); /* M = conductor of cuspidal part */
  return gc_long(av, ulcm(M, N));
}

static GEN
fs_get_MF(GEN fs) { return gel(fs,1); }
static GEN
fs_get_vES(GEN fs) { return gel(fs,2); }
static GEN
fs_get_pols(GEN fs) { return gel(fs,3); }
static GEN
fs_get_cosets(GEN fs) { return gel(fs,4); }
static long
fs_get_bitprec(GEN fs) { return itou(gel(fs,5)); }
static GEN
fs_get_vE(GEN fs) { return gel(fs,6); }
static GEN
fs_get_EF(GEN fs) { return gel(fs,7); }
static GEN
fs_get_expan(GEN fs) { return gel(fs,8); }
static GEN
fs_set_expan(GEN fs, GEN vanall)
{ GEN f = shallowcopy(fs); gel(f,8) = vanall; return f; }
static int
mfs_checkmf(GEN fs, GEN mf)
{ GEN mfF = fs_get_MF(fs); return gequal(gel(mfF,1), gel(mf,1)); }
static long
checkfs_i(GEN v)
{ return typ(v) == t_VEC && lg(v) == 9 && checkMF_i(fs_get_MF(v))
    && typ(fs_get_vES(v)) == t_VEC
    && typ(fs_get_pols(v)) == t_VEC
    && typ(fs_get_cosets(v)) == t_VEC
    && typ(fs_get_vE(v)) == t_VEC
    && lg(fs_get_pols(v)) == lg(fs_get_cosets(v))
    && typ(fs_get_expan(v)) == t_VEC
    && lg(fs_get_expan(v)) == 3
    && lg(gel(fs_get_expan(v), 1)) == lg(fs_get_cosets(v))
    && typ(gel(v,5)) == t_INT; }
GEN
checkMF_i(GEN mf)
{
  long l = lg(mf);
  GEN v;
  if (typ(mf) != t_VEC) return NULL;
  if (l == 9) return checkMF_i(fs_get_MF(mf));
  if (l != 7) return NULL;
  v = gel(mf,1);
  if (typ(v) != t_VEC || lg(v) != 5) return NULL;
  return (typ(gel(v,1)) == t_INT
         && typ(gmul2n(gel(v,2), 1)) == t_INT
         && typ(gel(v,3)) == t_VEC
         && typ(gel(v,4)) == t_INT)? mf: NULL; }
GEN
checkMF(GEN T)
{
  GEN mf = checkMF_i(T);
  if (!mf) pari_err_TYPE("checkMF [please use mfinit]", T);
  return mf;
}

/* c,d >= 0; c * Nc = N, find coset whose image in P1(Z/NZ) ~ (c, d + k(N/c)) */
static GEN
coset_complete(long c, long d, long Nc)
{
  long a, b;
  while (ugcd(c, d) > 1) d += Nc;
  (void)cbezout(c, d, &b, &a);
  return mkmat22s(a, -b, c, d);
}
/* right cosets of $\G_0(N)$: $\G=\bigsqcup_j \G_0(N)\ga_j$. */
/* We choose them with c\mid N and d mod N/c, not the reverse */
GEN
mfcosets(GEN gN)
{
  pari_sp av = avma;
  GEN V, D, mf;
  long l, i, ct, N = 0;
  if (typ(gN) == t_INT) N = itos(gN);
  else if ((mf = checkMF_i(gN))) N = MF_get_N(mf);
  else pari_err_TYPE("mfcosets", gN);
  if (N <= 0) pari_err_DOMAIN("mfcosets", "N", "<=", gen_0, stoi(N));
  V = cgetg(mypsiu(N) + 1, t_VEC);
  D = mydivisorsu(N); l = lg(D);
  for (i = ct = 1; i < l; i++)
  {
    long d, c = D[i], Nc = D[l-i], e = ugcd(Nc, c);
    for (d = 0; d < Nc; d++)
      if (ugcd(d,e) == 1) gel(V, ct++) = coset_complete(c, d, Nc);
  }
  return gerepilecopy(av, V);
}
static int
cmp_coset(void *E, GEN A, GEN B)
{
  ulong N = (ulong)E, Nc, c = itou(gcoeff(A,2,1));
  int r = cmpuu(c, itou(gcoeff(B,2,1)));
  if (r) return r;
  Nc = N / c;
  return cmpuu(umodiu(gcoeff(A,2,2), Nc), umodiu(gcoeff(B,2,2), Nc));
}
/* M in SL_2(Z) */
static long
mftocoset_i(ulong N, GEN M, GEN cosets)
{
  pari_sp av = avma;
  long A = itos(gcoeff(M,1,1)), c, u, v, Nc, i;
  long C = itos(gcoeff(M,2,1)), D = itos(gcoeff(M,2,2));
  GEN ga;
  c = cbezout(N*A, C, &u, &v); Nc = N/c;
  ga = coset_complete(c, umodsu(v*D, Nc), Nc);
  i = gen_search(cosets, ga, (void*)N, &cmp_coset);
  if (i < 0) pari_err_BUG("mftocoset [no coset found]");
  return gc_long(av,i);
}
/* (U * V^(-1))[2,2] mod N, assuming V in SL2(Z) */
static long
SL2_div_D(ulong N, GEN U, GEN V)
{
  long c = umodiu(gcoeff(U,2,1), N), d = umodiu(gcoeff(U,2,2), N);
  long a2 = umodiu(gcoeff(V,1,1), N), b2 = umodiu(gcoeff(V,1,2), N);
  return (a2*d - b2*c) % (long)N;
}
static long
mftocoset_iD(ulong N, GEN M, GEN cosets, long *D)
{
  long i = mftocoset_i(N, M, cosets);
  *D = SL2_div_D(N, M, gel(cosets,i)); return i;
}
GEN
mftocoset(ulong N, GEN M, GEN cosets)
{
  long i;
  if (!check_SL2Z(M)) pari_err_TYPE("mftocoset",M);
  i = mftocoset_i(N, M, cosets);
  retmkvec2(gdiv(M,gel(cosets,i)), utoipos(i));
}

static long
getnlim2(long N, long w1, long w2, long nlim, long k, long bitprec)
{
  if (w2 == N) return nlim;
  return mfperiod_prelim_double(1./sqrt((double)w1*w2), k, bitprec + 32);
}

/* g * S, g 2x2 */
static GEN
ZM_mulS(GEN g)
{ return mkmat2(gel(g,2), ZC_neg(gel(g,1))); }
/* g * T, g 2x2 */
static GEN
ZM_mulT(GEN g)
{ return mkmat2(gel(g,1), ZC_add(gel(g,2), gel(g,1))); }
/* g * T^(-1), g 2x2 */
static GEN
ZM_mulTi(GEN g)
{ return mkmat2(gel(g,1), ZC_sub(gel(g,2), gel(g,1))); }

/* Compute all slashexpansions for all cosets */
static GEN
mfgaexpansionall(GEN mf, GEN FE, GEN cosets, double height, long prec)
{
  GEN CHI = MF_get_CHI(mf), vres, vresaw;
  long lco, j, k = MF_get_k(mf), N = MF_get_N(mf), bitprec = prec2nbits(prec);

  lco = lg(cosets);
  vres = const_vec(lco-1, NULL);
  vresaw = cgetg(lco, t_VEC);
  for (j = 1; j < lco; j++) if (!gel(vres,j))
  {
    GEN ga = gel(cosets, j), van, aw, al, z, gai;
    long w1 = mfZC_width(N, gel(ga,1));
    long w2 = mfZC_width(N, gel(ga,2));
    long nlim, nlim2, daw, da, na, i;
    double sqNinvdbl = height ? height/w1 : 1./sqrt((double)w1*N);
    nlim = mfperiod_prelim_double(sqNinvdbl, k, bitprec + 32);
    van = mfslashexpansion(mf, FE, ga, nlim, 0, &aw, prec + EXTRAPREC64);
    van = vanembed(gel(FE, 1), van, prec + EXTRAPREC64);
    al = gel(aw, 1);
    nlim2 = height? nlim: getnlim2(N, w1, w2, nlim, k, bitprec);
    gel(vres, j) = vecslice(van, 1, nlim2 + 1);
    gel(vresaw, j) = aw;
    Qtoss(al, &na, &da); daw = da*w1;
    z = rootsof1powinit(1, daw, prec + EXTRAPREC64);
    gai = ga;
    for (i = 1; i < w1; i++)
    {
      GEN V, coe;
      long Di, n, ind, w2, s = ((i*na) % da) * w1, t = i*da;
      gai = ZM_mulT(gai);
      ind = mftocoset_iD(N, gai, cosets, &Di);
      w2 = mfZC_width(N, gel(gel(cosets,ind), 2));
      nlim2 = height? nlim: getnlim2(N, w1, w2, nlim, k, bitprec);
      gel(vresaw, ind) = aw;
      V = cgetg(nlim2 + 2, t_VEC);
      for (n = 0; n <= nlim2; n++, s = Fl_add(s, t, daw))
        gel(V, n+1) = gmul(gel(van, n+1), rootsof1pow(z, s));
      coe = mfcharcxeval(CHI, Di, prec + EXTRAPREC64);
      if (!gequal1(coe)) V = RgV_Rg_mul(V, conj_i(coe));
      gel(vres, ind) = V;
    }
  }
  return mkvec2(vres, vresaw);
}

/* Compute all period pols of F|_k\ga_j, vF = mftobasis(F_S) */
static GEN
mfperiodpols_i(GEN mf, GEN FE, GEN cosets, GEN *pvan, long bit)
{
  long N, i, prec = nbits2prec(bit), k = MF_get_k(mf);
  GEN vP, P, CHI, intall = gen_0;

  *pvan = gen_0;
  if (k == 0 && gequal0(gel(FE,2)))
    return cosets? const_vec(lg(cosets)-1, pol_0(0)): pol_0(0);
  N = MF_get_N(mf);
  CHI = MF_get_CHI(mf);
  P = get_P(k, fetch_var(), prec);
  if (!cosets)
  { /* ga = id */
    long nlim, PREC = prec + EXTRAPREC64;
    GEN F = gel(FE,1), sqNinv = invr(sqrtr_abs(utor(N, PREC))); /* A/w */
    GEN AR, v, van, T1, T2;

    nlim = mfperiod_prelim(sqNinv, k, bit + 32);
    /* F|id: al = 0, w = 1 */
    v = mfcoefs_i(F, nlim, 1);
    van = vanembed(F, v, PREC);
    AR = mkcomplex(gen_0, sqNinv);
    T1 = intAoo(van, nlim, gen_0,1, P, AR, k, prec);
    if (N == 1) T2 = T1;
    else
    { /* F|S: al = 0, w = N */
      v = mfgaexpansion(mf, FE, mkS(), nlim, PREC);
      van = vanembed(F, gel(v,3), PREC);
      AR = mkcomplex(gen_0, mulur(N,sqNinv));
      T2 = intAoo(van, nlim, gen_0,N, P, AR, k, prec);
    }
    T1 = gsub(T1, act_S(T2, k));
    T1 = normalizeapprox(T1, bit-20);
    vP = gprec_wtrunc(T1, prec);
  }
  else
  {
    long lco = lg(cosets);
    GEN vanall = mfgaexpansionall(mf, FE, cosets, 0, prec);
    *pvan = vanall;
    intall = intAoowithvanall(mf, vanall, P, cosets, bit);
    vP = const_vec(lco-1, NULL);
    for (i = 1; i < lco; i++)
    {
      GEN P, P1, P2, c, ga = gel(cosets, i);
      long iS, DS;
      if (gel(vP,i)) continue;
      P1 = gel(intall, i);
      iS = mftocoset_iD(N, ZM_mulS(ga), cosets, &DS);
      c = mfcharcxeval(CHI, DS, prec + EXTRAPREC64);
      P2 = gel(intall, iS);

      P = act_S(isint1(c)? P2: gmul(c, P2), k);
      P = normalizeapprox(gsub(P1, P), bit-20);
      gel(vP,i) = gprec_wtrunc(P, prec);
      if (iS == i) continue;

      P = act_S(isint1(c)? P1: gmul(conj_i(c), P1), k);
      if (!odd(k)) P = gneg(P);
      P = normalizeapprox(gadd(P, P2), bit-20);
      gel(vP,iS) = gprec_wtrunc(P, prec);
    }
  }
  delete_var(); return vP;
}

/* when cosets = NULL, return a "fake" symbol containing only fs(oo->0) */
static GEN
mfsymbol_i(GEN mf, GEN F, GEN cosets, long bit)
{
  GEN FE, van, vP, vE, Mvecj, vES = mftobasisES(mf,F);
  long precnew, prec = nbits2prec(bit), k = MF_get_k(mf);
  vE = mfgetembed(F, prec);
  Mvecj = obj_checkbuild(mf, MF_EISENSPACE, &mfeisensteinspaceinit);
  if (lg(Mvecj) >= 5) precnew = prec;
  else
  {
    long N = MF_get_N(mf), n = mfperiod_prelim_double(1/(double)N, k, bit + 32);
    precnew = prec + inveis_extraprec(N, mkS(), Mvecj, n);
  }
  FE = mkcol2(F, mf_eisendec(mf,F,precnew));
  vP = mfperiodpols_i(mf, FE, cosets, &van, bit);
  return mkvecn(8, mf, vES, vP, cosets, utoi(bit), vE, FE, van);
}

static GEN
fs2_get_cusps(GEN f) { return gel(f,3); }
static GEN
fs2_get_MF(GEN f) { return gel(f,1); }
static GEN
fs2_get_W(GEN f) { return gel(f,2); }
static GEN
fs2_get_F(GEN f) { return gel(f,4); }
static long
fs2_get_bitprec(GEN f) { return itou(gel(f,5)); }
static GEN
fs2_get_al0(GEN f) { return gel(f,6); }
static GEN
fs2_get_den(GEN f) { return gel(f,7); }
static int
checkfs2_i(GEN f)
{
  GEN W, C, F, al0;
  long l;
  if (typ(f) != t_VEC || lg(f) != 8 || typ(gel(f,5)) != t_INT) return 0;
  C = fs2_get_cusps(f); l = lg(C);
  W = fs2_get_W(f);
  F = fs2_get_F(f);
  al0 = fs2_get_al0(f);
  return checkMF_i(fs2_get_MF(f))
      && typ(W) == t_VEC && typ(F) == t_VEC && typ(al0) == t_VECSMALL
      && lg(W) == l && lg(F) == l && lg(al0) == l;
}
static GEN fs2_init(GEN mf, GEN F, long bit);
GEN
mfsymbol(GEN mf, GEN F, long bit)
{
  pari_sp av = avma;
  GEN cosets = NULL;
  if (!F)
  {
    F = mf;
    if (!checkmf_i(F)) pari_err_TYPE("mfsymbol", F);
    mf = mfinit_i(F, mf_FULL);
  }
  else if (!checkmf_i(F)) pari_err_TYPE("mfsymbol", F);
  if (checkfs2_i(mf)) return fs2_init(mf, F, bit);
  if (checkfs_i(mf))
  {
    cosets = fs_get_cosets(mf);
    mf = fs_get_MF(mf);
  }
  else if (checkMF_i(mf))
  {
    GEN gk = MF_get_gk(mf);
    if (typ(gk) != t_INT || equali1(gk)) return fs2_init(mf, F, bit);
    if (signe(gk) <= 0) pari_err_TYPE("mfsymbol [k <= 0]", mf);
    cosets = mfcosets(MF_get_gN(mf));
  }
  else pari_err_TYPE("mfsymbol",mf);
  return gerepilecopy(av, mfsymbol_i(mf, F, cosets, bit));
}

static GEN
RgX_by_parity(GEN P, long odd)
{
  long i, l = lg(P);
  GEN Q;
  if (l < 4) return odd ? pol_x(0): P;
  Q = cgetg(l, t_POL); Q[1] = P[1];
  for (i = odd? 2: 3; i < l; i += 2) gel(Q,i) = gen_0;
  for (i = odd? 3: 2; i < l; i += 2) gel(Q,i) = gel(P,i);
  return normalizepol_lg(Q, l);
}
/* flag 0: period polynomial of F, >0 or <0 with corresponding parity */
GEN
mfperiodpol(GEN mf0, GEN F, long flag, long bit)
{
  pari_sp av = avma;
  GEN pol, mf = checkMF_i(mf0);
  if (!mf) pari_err_TYPE("mfperiodpol",mf0);
  if (checkfs_i(F))
  {
    GEN mfpols = fs_get_pols(F);
    if (!mfs_checkmf(F, mf)) pari_err_TYPE("mfperiodpol [different mf]",F);
    pol = veclast(mfpols); /* trivial coset is last */
  }
  else
  {
    GEN gk = MF_get_gk(mf);
    if (typ(gk) != t_INT) pari_err_TYPE("mfperiodpol [half-integral k]", mf);
    if (equali1(gk)) pari_err_TYPE("mfperiodpol [k = 1]", mf);
    F = mfsymbol_i(mf, F, NULL, bit);
    pol = fs_get_pols(F);
  }
  if (flag) pol = RgX_by_parity(pol, flag < 0);
  return gerepilecopy(av, RgX_embedall(pol, fs_get_vE(F)));
}

static int
mfs_iscusp(GEN mfs) { return gequal0(gmael(mfs,2,1)); }
/* given cusps s1 and s2 (rationals or oo)
 * compute $\int_{s1}^{s2}(X-\tau)^{k-2}F|_k\ga_j(\tau)\,d\tau$ */
/* If flag = 1, do not give an error message if divergent, but
   give the rational function as result. */

static GEN
col2cusp(GEN v)
{
  GEN A, C;
  if (lg(v) != 3 || !RgV_is_ZV(v)) pari_err_TYPE("col2cusp",v);
  A = gel(v,1);
  C = gel(v,2);
  if (gequal0(C))
  {
    if (gequal0(A)) pari_err_TYPE("mfsymboleval", mkvec2(A, C));
    return mkoo();
  }
  return gdiv(A, C);
}
/* g.oo */
static GEN
mat2cusp(GEN g) { return col2cusp(gel(g,1)); }

static GEN
pathmattovec(GEN path)
{ return mkvec2(col2cusp(gel(path,1)), col2cusp(gel(path,2))); }

static void
get_mf_F(GEN fs, GEN *mf, GEN *F)
{
  if (lg(fs) == 3) { *mf = gel(fs,1); *F = gel(fs,2); }
  else { *mf = fs_get_MF(fs); *F = NULL; }
}
static GEN
mfgetvan(GEN fs, GEN ga, GEN *pal, long nlim, long prec)
{
  GEN van, mf, F, W;
  long PREC;
  get_mf_F(fs, &mf, &F);
  if (!F)
  {
    GEN vanall = fs_get_expan(fs), cosets = fs_get_cosets(fs);
    long D, jga = mftocoset_iD(MF_get_N(mf), ga, cosets, &D);
    van = gmael(vanall, 1, jga);
    W   = gmael(vanall, 2, jga);
    if (lg(van) >= nlim + 2)
    {
      GEN z = mfcharcxeval(MF_get_CHI(mf), D, prec);
      if (!gequal1(z)) van = RgV_Rg_mul(van, z);
      *pal = gel(W,1); return van;
    }
    F = gel(fs_get_EF(fs), 1);
  }
  PREC = prec + EXTRAPREC64;
  van = mfslashexpansion(mf, F, ga, nlim, 0, &W, PREC);
  van = vanembed(F, van, PREC);
  *pal = gel(W,1); return van;
}
/* Computation of int_A^oo (f | ga)(t)(X-t)^{k-2} dt, assuming convergence;
 * fs is either a symbol or a triple [mf,F,bitprec]. A != oo and im(A) > 0 */
static GEN
intAoo0(GEN fs, GEN A, GEN ga, GEN P, long bit)
{
  long nlim, N, k, w, prec = nbits2prec(bit);
  GEN van, mf, F, al;
  get_mf_F(fs, &mf,&F); N = MF_get_N(mf); k = MF_get_k(mf);
  w = mfZC_width(N, gel(ga,1));
  nlim = mfperiod_prelim(gdivgu(imag_i(A), w), k, bit + 32);
  van = mfgetvan(fs, ga, &al, nlim, prec);
  return intAoo(van, nlim, al,w, P, A, k, prec);
}

/* fs symbol, naive summation, A != oo, im(A) > 0 and B = oo or im(B) > 0 */
static GEN
mfsymboleval_direct(GEN fs, GEN path, GEN ga, GEN P)
{
  GEN A, B, van, S, al, mf = fs_get_MF(fs);
  long w, nlimA, nlimB = 0, N = MF_get_N(mf), k = MF_get_k(mf);
  long bit = fs_get_bitprec(fs), prec = nbits2prec(bit);

  A = gel(path, 1);
  B = gel(path, 2); if (typ(B) == t_INFINITY) B = NULL;
  w = mfZC_width(N, gel(ga,1));
  nlimA = mfperiod_prelim(gdivgu(imag_i(A),w), k, bit + 32);
  if (B) nlimB = mfperiod_prelim(gdivgu(imag_i(B),w), k, bit + 32);
  van = mfgetvan(fs, ga, &al, maxss(nlimA,nlimB), prec);
  S = intAoo(van, nlimA, al,w, P, A, k, prec);
  if (B) S = gsub(S, intAoo(van, nlimB, al,w, P, B, k, prec));
  return RgX_embedall(S, fs_get_vE(fs));
}

/* Computation of int_A^oo (f | ga)(t)(X-t)^{k-2} dt, assuming convergence;
 * fs is either a symbol or a pair [mf,F]. */
static GEN
mfsymbolevalpartial(GEN fs, GEN A, GEN ga, long bit)
{
  GEN Y, F, S, P, mf;
  long N, k, w, prec = nbits2prec(bit);

  get_mf_F(fs, &mf, &F);
  N = MF_get_N(mf); w = mfZC_width(N, gel(ga,1));
  k = MF_get_k(mf);
  Y = gdivgu(imag_i(A), w);
  P = get_P(k, fetch_var(), prec);
  if (lg(fs) != 3 && gtodouble(Y)*(2*N) < 1)
  { /* true symbol + low imaginary part: use GL_2 action to improve */
    GEN U, ga2, czd, A2 = cxredga0N(N, A, &U, &czd, 1);
    GEN vE = fs_get_vE(fs);
    ga2 = ZM_mul(ga, ZM_inv(U, NULL));
    S = RgX_embedall(intAoo0(fs, A2, ga2, P, bit), vE);
    S = gsub(S, mfsymboleval(fs, mkvec2(mat2cusp(U), mkoo()), ga2, bit));
    S = typ(S) == t_VEC? vecact_GL2(S, U, k): act_GL2(S, U, k);
  }
  else
  {
    S = intAoo0(fs, A, ga, P, bit);
    S = RgX_embedall(S, F? mfgetembed(F,prec): fs_get_vE(fs));
  }
  delete_var(); return normalizeapprox(S, bit-20);
}

static GEN
actal(GEN x, GEN vabd)
{
  if (typ(x) == t_INFINITY) return x;
  return gdiv(gadd(gmul(gel(vabd,1), x), gel(vabd,2)), gel(vabd,3));
}

static GEN
unact(GEN z, GEN vabd, long k, long prec)
{
  GEN res = gsubst(z, 0, actal(pol_x(0), vabd));
  GEN CO = gpow(gdiv(gel(vabd,3), gel(vabd,1)), sstoQ(k-2, 2), prec);
  return gmul(CO, res);
}

GEN
mfsymboleval(GEN fs, GEN path, GEN ga, long bitprec)
{
  pari_sp av = avma;
  GEN tau, V, LM, S, CHI, mfpols, cosets, al, be, mf, F, vabd = NULL;
  long D, B, m, u, v, a, b, c, d, j, k, N, prec, tsc1, tsc2;

  if (checkfs_i(fs))
  {
    get_mf_F(fs, &mf, &F);
    bitprec = minss(bitprec, fs_get_bitprec(fs));
  }
  else
  {
    if (checkfs2_i(fs)) pari_err_TYPE("mfsymboleval [need integral k > 1]",fs);
    if (typ(fs) != t_VEC || lg(fs) != 3) pari_err_TYPE("mfsymboleval",fs);
    get_mf_F(fs, &mf, &F);
    mf = checkMF_i(mf);
    if (!mf ||!checkmf_i(F)) pari_err_TYPE("mfsymboleval",fs);
  }
  if (lg(path) != 3) pari_err_TYPE("mfsymboleval",path);
  if (typ(path) == t_MAT) path = pathmattovec(path);
  if (typ(path) != t_VEC) pari_err_TYPE("mfsymboleval",path);
  al = gel(path,1);
  be = gel(path,2);
  ga = ga? GL2toSL2(ga, &vabd): matid(2);
  if (vabd)
  {
    al = actal(al, vabd);
    be = actal(be, vabd); path = mkvec2(al, be);
  }
  tsc1 = cusp_AC(al, &a, &c);
  tsc2 = cusp_AC(be, &b, &d);
  prec = nbits2prec(bitprec);
  k = MF_get_k(mf);
  if (!tsc1)
  {
    GEN z2, z = mfsymbolevalpartial(fs, al, ga, bitprec);
    if (tsc2)
      z2 = d? mfsymboleval(fs, mkvec2(be, mkoo()), ga, bitprec): gen_0;
    else
      z2 = mfsymbolevalpartial(fs, be, ga, bitprec);
    z = gsub(z, z2);
    if (vabd) z = unact(z, vabd, k, prec);
    return gerepileupto(av, normalizeapprox(z, bitprec-20));
  }
  else if (!tsc2)
  {
    GEN z = mfsymbolevalpartial(fs, be, ga, bitprec);
    if (c) z = gsub(mfsymboleval(fs, mkvec2(al, mkoo()), ga, bitprec), z);
    if (vabd) z = unact(z, vabd, k, prec);
    return gerepileupto(av, normalizeapprox(z, bitprec-20));
  }
  if (F) pari_err_TYPE("mfsymboleval", fs);
  D = a*d-b*c;
  if (!D) { set_avma(av); return RgX_embedall(gen_0, fs_get_vE(fs)); }
  mfpols = fs_get_pols(fs);
  cosets = fs_get_cosets(fs);
  CHI = MF_get_CHI(mf); N = MF_get_N(mf);
  cbezout(a, c, &u, &v); B = u*b + v*d; tau = mkmat22s(a, -v, c, u);
  V = gcf(sstoQ(B, D));
  LM = shallowconcat(mkcol2(gen_1, gen_0), contfracpnqn(V, lg(V)));
  S = gen_0; m = lg(LM) - 2;
  for (j = 0; j < m; j++)
  {
    GEN M, P;
    long D, iN;
    M = mkmat2(gel(LM, j+2), gel(LM, j+1));
    if (!odd(j)) gel(M,1) = ZC_neg(gel(M,1));
    M = ZM_mul(tau, M);
    iN = mftocoset_iD(N, ZM_mul(ga, M), cosets, &D);
    P = gmul(gel(mfpols,iN), mfcharcxeval(CHI,D,prec));
    S = gadd(S, act_GL2(P, ZM_inv(M, NULL), k));
  }
  if (typ(S) == t_RFRAC)
  {
    GEN R, S1, co;
    gel(S,2) = primitive_part(gel(S,2), &co);
    if (co) gel(S,1) = gdiv(gel(S,1), gtofp(co,prec));
    S1 = poldivrem(gel(S,1), gel(S,2), &R);
    if (gexpo(R) < -bitprec + 20) S = S1;
  }
  if (vabd) S = unact(S, vabd, k, prec);
  S = RgX_embedall(S, fs_get_vE(fs));
  return gerepileupto(av, normalizeapprox(S, bitprec-20));
}

/* v a scalar or t_POL; set *pw = a if expo(a) > E for some coefficient;
 * take the 'a' with largest exponent */
static void
improve(GEN v, GEN *pw, long *E)
{
  if (typ(v) != t_POL)
  {
    long e = gexpo(v);
    if (e > *E) { *E = e; *pw = v; }
  }
  else
  {
    long j, l = lg(v);
    for (j = 2; j < l; j++) improve(gel(v,j), pw, E);
  }
}
static GEN
polabstorel(GEN rnfeq, GEN T)
{
  long i, l;
  GEN U;
  if (typ(T) != t_POL) return T;
  U = cgetg_copy(T, &l); U[1] = T[1];
  for (i = 2; i < l; i++) gel(U,i) = eltabstorel(rnfeq, gel(T,i));
  return U;
}
static GEN
bestapprnfrel(GEN x, GEN polabs, GEN roabs, GEN rnfeq, long prec)
{
  x = bestapprnf(x, polabs, roabs, prec);
  if (rnfeq) x = polabstorel(rnfeq, liftpol_shallow(x));
  return x;
}
/* v vector of polynomials polynomial in C[X] (possibly scalar).
 * Set *w = coeff with largest exponent and return T / *w, rationalized */
static GEN
normal(GEN v, GEN polabs, GEN roabs, GEN rnfeq, GEN *w, long prec)
{
  long i, l = lg(v), E = -(long)HIGHEXPOBIT;
  GEN dv;
  for (i = 1; i < l; i++) improve(gel(v,i), w, &E);
  v = RgV_Rg_mul(v, ginv(*w));
  for (i = 1; i < l; i++)
    gel(v,i) = bestapprnfrel(gel(v,i), polabs,roabs,rnfeq,prec);
  v = Q_primitive_part(v,&dv);
  if (dv) *w = gmul(*w,dv);
  return v;
}

static GEN mfpetersson_i(GEN FS, GEN GS);

GEN
mfmanin(GEN FS, long bitprec)
{
  pari_sp av = avma;
  GEN mf, M, vp, vm, cosets, CHI, vpp, vmm, f, T, P, vE, polabs, roabs, rnfeq;
  GEN pet;
  long N, k, lco, i, prec, lvE;

  if (!checkfs_i(FS))
  {
    if (checkfs2_i(FS)) pari_err_TYPE("mfmanin [need integral k > 1]",FS);
    pari_err_TYPE("mfmanin",FS);
  }
  if (!mfs_iscusp(FS)) pari_err_TYPE("mfmanin [noncuspidal]",FS);
  mf = fs_get_MF(FS);
  vp = fs_get_pols(FS);
  cosets = fs_get_cosets(FS);
  bitprec = fs_get_bitprec(FS);
  N = MF_get_N(mf); k = MF_get_k(mf); CHI = MF_get_CHI(mf);
  lco = lg(cosets); vm = cgetg(lco, t_VEC);
  prec = nbits2prec(bitprec);
  for (i = 1; i < lco; i++)
  {
    GEN g = gel(cosets, i), c;
    long A = itos(gcoeff(g,1,1)), B = itos(gcoeff(g,1,2));
    long C = itos(gcoeff(g,2,1)), D = itos(gcoeff(g,2,2));
    long Dbar, ibar = mftocoset_iD(N, mkmat22s(-B,-A,D,C), cosets, &Dbar);

    c = mfcharcxeval(CHI, Dbar, prec); if (odd(k)) c = gneg(c);
    T = RgX_Rg_mul(gel(vp,ibar), c);
    if (typ(T) == t_POL && varn(T) == 0) T = RgX_recip(T);
    gel(vm,i) = T;
  }
  vpp = gadd(vp,vm);
  vmm = gsub(vp,vm);

  vE = fs_get_vE(FS); lvE = lg(vE);
  f = gel(fs_get_EF(FS), 1);
  P = mf_get_field(f); if (degpol(P) == 1) P = NULL;
  T = mfcharpol(CHI);  if (degpol(T) == 1) T = NULL;
  if (T && P)
  {
    rnfeq = nf_rnfeqsimple(T, P);
    polabs = gel(rnfeq,1);
    roabs = gel(QX_complex_roots(polabs,prec), 1);
  }
  else
  {
    rnfeq = roabs = NULL;
    polabs = P? P: T;
  }
  pet = mfpetersson_i(FS, NULL);
  M = cgetg(lvE, t_VEC);
  for (i = 1; i < lvE; i++)
  {
    GEN p, m, wp, wm, petdiag, r, E = gel(vE,i);
    p = normal(RgXV_embed(vpp, E), polabs, roabs, rnfeq, &wp, prec);
    m = normal(RgXV_embed(vmm, E), polabs, roabs, rnfeq, &wm, prec);
    petdiag = typ(pet)==t_MAT? gcoeff(pet,i,i): pet;
    r = gdiv(imag_i(gmul(wp, conj_i(wm))), petdiag);
    r = bestapprnfrel(r, polabs, roabs, rnfeq, prec);
    gel(M,i) = mkvec2(mkvec2(p,m), mkvec3(wp,wm,r));
  }
  return gerepilecopy(av, lvE == 2? gel(M,1): M);
}

/* flag = 0: full, flag = +1 or -1, odd/even */
/* Basis of period polynomials in level 1. */
GEN
mfperiodpolbasis(long k, long flag)
{
  pari_sp av = avma;
  long i, j, n = k - 2;
  GEN M, C, v;
  if (k <= 4) return cgetg(1,t_VEC);
  M = cgetg(k, t_MAT);
  C = matpascal(n);
  if (!flag)
    for (j = 0; j <= n; j++)
    {
      gel(M, j+1) = v = cgetg(k, t_COL);
      for (i = 0; i <= j; i++) gel(v, i+1) = gcoeff(C, j+1, i+1);
      for (; i <= n; i++) gel(v, i+1) = gcoeff(C, n-j+1, i-j+1);
    }
  else
    for (j = 0; j <= n; j++)
    {
      gel(M, j+1) = v = cgetg(k, t_COL);
      for (i = 0; i <= n; i++)
      {
        GEN a = i < j ? gcoeff(C, j+1, i+1) : gen_0;
        if (i + j >= n)
        {
          GEN b = gcoeff(C, j+1, i+j-n+1);
          a = flag < 0 ? addii(a,b) : subii(a,b);
        }
        gel(v, i+1) = a;
      }
    }
  return gerepilecopy(av, RgM_to_RgXV(ZM_ker(M), 0));
}

static int
zero_at_cusp(GEN mf, GEN F, GEN c)
{
  GEN v = evalcusp(mf, F, c, LOWDEFAULTPREC);
  return (gequal0(v) || gexpo(v) <= -62);
}
/* Compute list E of j such that F|_k g_j vanishes at oo: return [E, s(E)] */
static void
mffvanish(GEN mf, GEN F, GEN G, GEN cosets, GEN *pres, GEN *press)
{
  long j, lc = lg(cosets), N = MF_get_N(mf);
  GEN v, vs;
  *pres = v  = zero_zv(lc-1);
  *press= vs = zero_zv(lc-1);
  for (j = 1; j < lc; j++)
  {
    GEN ga = gel(cosets,j), c = mat2cusp(ga);
    if (zero_at_cusp(mf, F, c))
      v[j] = vs[ mftocoset_i(N, ZM_mulS(ga), cosets) ] = 1;
    else if (!zero_at_cusp(mf, G, c))
      pari_err_IMPL("divergent Petersson product");
  }
}
static GEN
Haberland(GEN PF, GEN PG, GEN vEF, GEN vEG, long k)
{
  GEN S = gen_0, vC = vecbinomial(k-2); /* vC[n+1] = (-1)^n binom(k-2,n) */
  long n, j, l = lg(PG);
  for (n = 2; n < k; n+=2) gel(vC,n) = negi(gel(vC,n));
  for (j = 1; j < l; j++)
  {
    GEN PFj = gel(PF,j), PGj = gel(PG,j);
    for (n = 0; n <= k-2; n++)
    {
      GEN a = RgX_coeff(PGj, k-2-n), b = RgX_coeff(PFj, n);
      a = Rg_embedall(a, vEG);
      b = Rg_embedall(b, vEF);
      a = conj_i(a); if (typ(a) == t_VEC) settyp(a, t_COL);
      /* a*b = scalar or t_VEC or t_COL or t_MAT */
      S = gadd(S, gdiv(gmul(a,b), gel(vC,n+1)));
    }
  }
  S = mulcxpowIs(gmul2n(S, 1-k), 1+k);
  return vEF==vEG? real_i(S): S;
}
/* F1S, F2S both symbols, same mf */
static GEN
mfpeterssonnoncusp(GEN F1S, GEN F2S)
{
  pari_sp av = avma;
  GEN mf, F1, F2, GF1, GF2, P2, cosets, vE1, vE2, FE1, FE2, P;
  GEN I, IP1, RHO, RHOP1, INF, res, ress;
  const double height = sqrt(3.)/2;
  long k, r, j, bitprec, prec;

  mf = fs_get_MF(F1S);
  FE1 = fs_get_EF(F1S); F1 = gel(FE1, 1);
  FE2 = fs_get_EF(F2S); F2 = gel(FE2, 1);
  cosets = fs_get_cosets(F1S);
  bitprec = minuu(fs_get_bitprec(F1S), fs_get_bitprec(F2S));
  prec = nbits2prec(bitprec);
  F1S = fs_set_expan(F1S, mfgaexpansionall(mf, FE1, cosets, height, prec));
  if (F2S != F1S)
    F2S = fs_set_expan(F2S, mfgaexpansionall(mf, FE2, cosets, height, prec));
  k = MF_get_k(mf); r = lg(cosets) - 1;
  vE1 = fs_get_vE(F1S);
  vE2 = fs_get_vE(F2S);
  I = gen_I();
  IP1 = mkcomplex(gen_1,gen_1);
  RHO = rootsof1u_cx(3, prec+EXTRAPREC64);
  RHOP1 = gaddsg(1, RHO);
  INF = mkoo();
  mffvanish(mf, F1, F2, cosets, &res, &ress);
  P2 = fs_get_pols(F2S);
  GF1 = cgetg(r+1, t_VEC);
  GF2 = cgetg(r+1, t_VEC); P = get_P(k, fetch_var(), prec);
  for (j = 1; j <= r; j++)
  {
    GEN g = gel(cosets,j);
    if (res[j]) {
      gel(GF1,j) = mfsymboleval_direct(F1S, mkvec2(RHOP1,INF), g, P);
      gel(GF2,j) = mfsymboleval_direct(F2S, mkvec2(I,IP1), g, P);
    } else if (ress[j]) {
      gel(GF1,j) = mfsymboleval_direct(F1S, mkvec2(RHOP1,RHO), g, P);
      gel(GF2,j) = mfsymboleval_direct(F2S, mkvec2(I,INF), g, P);
    } else {
      gel(GF1,j) = mfsymboleval_direct(F1S, mkvec2(RHO,I), g, P);
      gel(GF2,j) = gneg(gel(P2,j)); /* - symboleval(F2S, [0,oo] */
    }
  }
  delete_var();
  return gerepileupto(av, gdivgu(Haberland(GF1,GF2, vE1,vE2, k), r));
}

/* Petersson product of F and G, given by mfsymbol's [k > 1 integral] */
static GEN
mfpetersson_i(GEN FS, GEN GS)
{
  pari_sp av = avma;
  GEN mf, ESF, ESG, PF, PG, PH, CHI, cosets, vEF, vEG;
  long k, r, j, N, bitprec, prec;

  if (!checkfs_i(FS)) pari_err_TYPE("mfpetersson",FS);
  mf = fs_get_MF(FS);
  ESF = fs_get_vES(FS);
  if (!GS) GS = FS;
  else
  {
    if (!checkfs_i(GS)) pari_err_TYPE("mfpetersson",GS);
    if (!mfs_checkmf(GS, mf))
      pari_err_TYPE("mfpetersson [different mf]", mkvec2(FS,GS));
  }
  ESG = fs_get_vES(GS);
  if (!gequal0(gel(ESF,1)) && !gequal0(gel(ESG,1)))
    return mfpeterssonnoncusp(FS, GS);
  if (gequal0(gel(ESF,2)) || gequal0(gel(ESG,2))) return gc_const(av, gen_0);
  N = MF_get_N(mf);
  k = MF_get_k(mf);
  CHI = MF_get_CHI(mf);
  PF = fs_get_pols(FS); vEF = fs_get_vE(FS);
  PG = fs_get_pols(GS); vEG = fs_get_vE(GS);
  cosets = fs_get_cosets(FS);
  bitprec = minuu(fs_get_bitprec(FS), fs_get_bitprec(GS));
  prec = nbits2prec(bitprec);
  r = lg(PG)-1;
  PH = cgetg(r+1, t_VEC);
  for (j = 1; j <= r; j++)
  {
    GEN ga = gel(cosets,j), PGj1, PGjm1;
    long iT, D;
    iT = mftocoset_iD(N, ZM_mulTi(ga), cosets, &D);
    PGj1 = RgX_translate(gel(PG, iT), gen_1);
    PGj1 = RgX_Rg_mul(PGj1, mfcharcxeval(CHI, D, prec));
    iT = mftocoset_iD(N, ZM_mulT(ga), cosets, &D);
    PGjm1 = RgX_translate(gel(PG,iT), gen_m1);
    PGjm1 = RgX_Rg_mul(PGjm1, mfcharcxeval(CHI, D, prec));
    gel(PH,j) = gsub(PGj1, PGjm1);
  }
  return gerepileupto(av, gdivgu(Haberland(PF, PH, vEF, vEG, k), 6*r));
}

/****************************************************************/
/*           Petersson products using Nelson-Collins            */
/****************************************************************/
/* Compute W(k,z) = sum_{m >= 1} (mz)^{k-1}(mzK_{k-2}(mz)-K_{k-1}(mz))
 * for z>0 and absolute accuracy < 2^{-B}.
 * K_k(x) ~ (Pi/(2x))^{1/2} e^{-x} */

static void
Wparams(GEN *ph, long *pN, long k, double x, long prec)
{
  double B = prec2nbits(prec) + 10;
  double C = B + k*log(x)/M_LN2 + 1, D = C*M_LN2 + 2.065;
  double F = 2 * ((C - 1) * M_LN2 + log(gtodouble(mpfact(k)))) / x;
  double T = log(F) * (1 + 2*k/x/F), PI2 = M_PI*M_PI;
  *pN = (long)ceil((T/PI2) * (D + log(D/PI2)));
  *ph = gprec_w(dbltor(T / *pN), prec);
}

static void
Wcoshall(GEN *pCH, GEN *pCHK, GEN *pCHK1, GEN h, long N, long k, long prec)
{
  GEN CH, CHK, CHK1, z = gexp(h, prec);
  GEN PO = gpowers(z, N), POK1 = gpowers(gpowgs(z, k-1), N);
  GEN E = ginv(gel(PO, N + 1)); /* exp(-hN) */
  GEN E1 = ginv(gel(POK1, N + 1)); /* exp(-(k-1)h) */
  long j;
  *pCH = CH = cgetg(N+1, t_VEC);
  *pCHK = CHK = cgetg(N+1, t_VEC);
  *pCHK1 = CHK1 = cgetg(N+1, t_VEC);
  for (j = 1; j <= N; j++)
  {
    GEN eh = gel(PO, j+1), emh = gmul(gel(PO, N-j+1), E); /* e^{jh}, e^{-jh} */
    GEN ek1h = gel(POK1, j+1), ek1mh = gmul(gel(POK1, N-j+1), E1);
    gel(CH, j) = gmul2n(gadd(eh, emh), -1); /* cosh(jh) */
    gel(CHK1,j) = gmul2n(gadd(ek1h, ek1mh), -1); /* cosh((k-1)jh) */
    gel(CHK, j) = gmul2n(gadd(gmul(eh, ek1h), gmul(emh, ek1mh)), -1);
  }
}

/* computing W(k,x) via integral */
static GEN
Wint(long k, GEN vP, GEN x, long prec)
{
  GEN P, P1, S1, S, h, CH, CHK, CHK1;
  long N, j;
  Wparams(&h, &N, k, gtodouble(x), prec);
  Wcoshall(&CH, &CHK, &CHK1, h, N, k, prec);
  P = gel(vP, k+1); P1 = gel(vP, k); S = S1 = NULL;
  for (j = 0; j <= N; j++)
  {
    GEN eh = gexp(j? gmul(x, gel(CH, j)): x, prec);
    GEN eh1 = gsubgs(eh, 1), eh1k = gpowgs(eh1, k), t1, t;
    t = gdiv(poleval(P, eh), gmul(eh1, eh1k));
    t1 = gdiv(poleval(P1, eh), eh1k);
    if (j)
    {
      S = gadd(S, gmul(t, gel(CHK, j)));
      S1 = gadd(S1, gmul(t1, gel(CHK1, j)));
    }
    else
    {
      S = gmul2n(t, -1);
      S1 = gmul2n(t1, -1);
    }
  }
  return gmul(gmul(h, gpowgs(x, k-1)), gsub(gmul(x, S), gmulsg(2*k-1, S1)));
}

static GEN
get_vP(long k)
{
  GEN P, v = cgetg(k+2, t_VEC), Q = deg1pol_shallow(gen_1,gen_m1,0);
  long j;
  gel(v,1) = gen_1;
  gel(v,2) = P = pol_x(0);
  for (j = 2; j <= k; j++)
    gel(v,j+1) = P = RgX_shift_shallow(gsub(gmulsg(j, P),
                                            gmul(Q, ZX_deriv(P))), 1);
  return v;
}
/* vector of (-1)^j(1/(exp(x)-1))^(j) [x = z] * z^j for 0<=j<=r */
static GEN
VS(long r, GEN z, GEN V, long prec)
{
  GEN e = gexp(z, prec), c = ginv(gsubgs(e,1));
  GEN T = gpowers0(gmul(c, z), r, c);
  long j;
  V = gsubst(V, 0, e);
  for (j = 1; j <= r + 1; j++) gel(V,j) = gmul(gel(V,j), gel(T,j));
  return V;
}

/* U(r,x)=sum_{m >= 1} (mx)^k K_k(mx), k = r+1/2 */
static GEN
Unelson(long r, GEN V)
{
  GEN S = gel(V,r+1), C = gen_1; /* (r+j)! / j! / (r-j)! */
  long j;
  if (!r) return S;
  for (j = 1; j <= r; j++)
  {
    C = gdivgu(gmulgu(C, (r+j)*(r-j+1)), j);
    S = gadd(S, gmul2n(gmul(C, gel(V, r-j+1)), -j));
  }
  return S;
}
/* W(r+1/2,z) / sqrt(Pi/2) */
static GEN
Wint2(long r, GEN vP, GEN z, long prec)
{
  GEN R, V = VS(r, z, vP, prec);
  R = Unelson(r, V);
  if (r) R = gsub(R, gmulsg(2*r, Unelson(r-1, V)));
  return R;
}
typedef GEN(*Wfun_t)(long, GEN, GEN, long);
static GEN
WfromZ(GEN Z, GEN vP, GEN gkm1, Wfun_t W, long k, GEN pi4, long prec)
{
  pari_sp av = avma;
  GEN Zk = gpow(Z, gkm1, prec), z = gmul(pi4, gsqrt(Z,prec));
  return gerepileupto(av, gdiv(W(k, vP, z, prec), Zk));
}
/* mf a true mf or an fs2 */
static GEN
fs2_init(GEN mf, GEN F, long bit)
{
  pari_sp av = avma;
  long i, l, lim, N, k, k2, prec = nbits2prec(bit);
  GEN DEN, cusps, tab, gk, gkm1, W0, vW, vVW, vVF, vP, al0;
  GEN vE = mfgetembed(F, prec), pi4 = Pi2n(2, prec);
  Wfun_t Wf;

  if (lg(mf) == 7)
  {
    vW = cusps = NULL; /* true mf */
    DEN = tab = NULL; /* -Wall */
  }
  else
  { /* mf already an fs2, reset if its precision is too low */
    vW = (fs2_get_bitprec(mf) < bit)? NULL: fs2_get_W(mf);
    cusps = fs2_get_cusps(mf);
    DEN = fs2_get_den(mf);
    mf = fs2_get_MF(mf);
  }
  N = MF_get_N(mf);
  gk = MF_get_gk(mf); gkm1 = gsubgs(gk, 1);
  k2 = itos(gmul2n(gk,1));
  Wf = odd(k2)? Wint2: Wint;
  k = k2 >> 1; vP = get_vP(k);
  if (vW) lim = (lg(gel(vW,1)) - 2) / N; /* vW[1] attached to cusp 0, width N */
  else
  { /* true mf */
    double B = (bit + 10)*M_LN2;
    double L = (B + k2*log(B)/2 + k2*k2*log(B)/(4*B)) / (4*M_PI);
    long n, Lw;
    lim = ((long)ceil(L*L));
    Lw = N*lim;
    tab = cgetg(Lw+1,t_VEC);
    for (n = 1; n <= Lw; n++)
      gel(tab,n) = WfromZ(uutoQ(n,N), vP, gkm1, Wf, k, pi4, prec);
    if (!cusps) cusps = mfcusps_i(N);
    DEN = gmul2n(gmulgu(gpow(Pi2n(3, prec), gkm1, prec), mypsiu(N)), -2);
    if (odd(k2)) DEN = gdiv(DEN, sqrtr_abs(Pi2n(-1,prec)));
  }
  l = lg(cusps);
  vVF = cgetg(l, t_VEC);
  vVW = cgetg(l, t_VEC);
  al0 = cgetg(l, t_VECSMALL);
  W0 = k2==1? ginv(pi4): gen_0;
  for (i = 1; i < l; i++)
  {
    long A, C, w, wi, Lw, n;
    GEN VF, W, paramsF, al;
    (void)cusp_AC(gel(cusps,i), &A,&C);
    wi = ugcd(N, C*C); w = N / wi; Lw = w * lim;
    VF = mfslashexpansion(mf, F, cusp2mat(A,C), Lw, 0, &paramsF, prec);
    /* paramsF[2] = w */
    al = gel(paramsF, 1); if (gequal0(al)) al = NULL;
    for (n = 0; n <= Lw; n++)
    {
      GEN a = gel(VF,n+1);
      gel(VF,n+1) = gequal0(a)? gen_0: Rg_embedall(a, vE);
    }
    if (vW)
      W = gel(vW, i);
    else
    {
      W = cgetg(Lw+2, t_VEC);
      for (n = 0; n <= Lw; n++)
        gel(W, n+1) = al? WfromZ(gadd(al,uutoQ(n,w)),vP,gkm1,Wf,k,pi4, prec)
                        : (n? gel(tab, n * wi): W0);
    }
    al0[i] = !al;
    gel(vVF, i) = VF;
    gel(vVW, i) = W;
  }
  if (k2 <= 1) al0 = zero_zv(l-1); /* no need to test for convergence */
  return gerepilecopy(av, mkvecn(7, mf,vVW,cusps,vVF,utoipos(bit),al0,DEN));
}

static GEN
mfpetersson2(GEN Fs, GEN Gs)
{
  pari_sp av = avma;
  GEN VC, RES, vF, vG, vW = fs2_get_W(Fs), al0 = fs2_get_al0(Fs);
  long N = MF_get_N(fs2_get_MF(Fs)), j, lC;

  VC = fs2_get_cusps(Fs); lC = lg(VC);
  vF = fs2_get_F(Fs);
  vG = Gs? fs2_get_F(Gs): vF;
  RES = gen_0;
  for (j = 1; j < lC; j++)
  {
    GEN W = gel(vW,j), VF = gel(vF,j), VG = gel(vG,j), T = gen_0;
    long A, C, w, n, L = lg(W);
    pari_sp av = avma;
    (void)cusp_AC(gel(VC,j), &A,&C); w = N/ugcd(N, C*C);
    if (al0[j] && !isintzero(gel(VF,1)) && !isintzero(gel(VG,1)))
      pari_err_IMPL("divergent Petersson product");
    for (n = 1; n < L; n++)
    {
      GEN b = gel(VF,n), a = gel(VG,n);
      if (!isintzero(a) && !isintzero(b))
      {
        T = gadd(T, gmul(gel(W,n), gmul(conj_i(a),b)));
        if (gc_needed(av,2)) T = gerepileupto(av,T);
      }
    }
    if (w != 1) T = gmulgu(T,w);
    RES = gerepileupto(av, gadd(RES, T));
  }
  if (!Gs) RES = real_i(RES);
  return gerepileupto(av, gdiv(RES, fs2_get_den(Fs)));
}

static long
symbol_type(GEN F)
{
  if (checkfs_i(F)) return 1;
  if (checkfs2_i(F)) return 2;
  return 0;
}
static int
symbol_same_mf(GEN F, GEN G) { return gequal(gmael(F,1,1), gmael(G,1,1)); }
GEN
mfpetersson(GEN F, GEN G)
{
  long tF = symbol_type(F);
  if (!tF) pari_err_TYPE("mfpetersson",F);
  if (G)
  {
    long tG = symbol_type(G);
    if (!tG) pari_err_TYPE("mfpetersson",F);
    if (tF != tG || !symbol_same_mf(F,G))
      pari_err_TYPE("mfpetersson [incompatible symbols]", mkvec2(F,G));
  }
  return (tF == 1)? mfpetersson_i(F, G): mfpetersson2(F, G);
}

/****************************************************************/
/*         projective Galois representation, weight 1           */
/****************************************************************/
static void
moreorders(long N, GEN CHI, GEN F, GEN *pP, GEN *pO, ulong *bound)
{
  pari_sp av = avma;
  forprime_t iter;
  ulong a = *bound+1, b = 2*(*bound), p;
  long i = 1;
  GEN P, O, V = mfcoefs_i(F, b, 1);
  *bound = b;
  P = cgetg(b-a+2, t_VECSMALL);
  O = cgetg(b-a+2, t_VECSMALL);
  u_forprime_init(&iter, a, b);
  while((p = u_forprime_next(&iter))) if (N % p)
  {
    O[i] = mffindrootof1(V, p, CHI);
    P[i++] = p;
  }
  setlg(P, i); *pP = shallowconcat(*pP, P);
  setlg(O, i); *pO = shallowconcat(*pO, O);
  gerepileall(av, 2, pP, pO);
}

static GEN
search_abelian(GEN nf, long n, long k, GEN N, GEN CHI, GEN F,
               GEN *pP, GEN *pO, ulong *bound, long prec)
{
  pari_sp av = avma;
  GEN bnr, cond, H, cyc, gn, T, Bquo, P, E;
  long sN = itos(N), r1 = nf_get_r1(nf), i, j, d;

  cond = idealfactor(nf, N);
  P = gel(cond,1);
  E = gel(cond,2);
  for (i = j = 1; i < lg(P); i++)
  {
    GEN pr = gel(P,i), Ej = gen_1;
    long p = itos(pr_get_p(pr));
    if (p == n)
    {
      long e = pr_get_e(pr); /* 1 + [e*p/(p-1)] */
      Ej = utoipos(1 + (e*p) / (p-1));
    }
    else
    {
      long f = pr_get_f(pr);
      if (Fl_powu(p % n, f, n) != 1) continue;
    }
    gel(P,j) = pr;
    gel(E,j) = Ej; j++;
  }
  setlg(P,j);
  setlg(E,j);
  cond = mkvec2(cond, const_vec(r1, gen_1));
  bnr = Buchraymod(Buchall(nf, nf_FORCE, prec), cond, nf_INIT, utoipos(n));
  cyc = bnr_get_cyc(bnr);
  d = lg(cyc)-1;
  H = zv_diagonal(ZV_to_Flv(cyc, n));
  gn = utoi(n);
  for (i = 1;;)
  {
    for(j = 2; i < lg(*pO); i++)
    {
      long o, q = (*pP)[i];
      GEN pr = idealprimedec_galois(nf, stoi(q));
      o = ((*pO)[i] / pr_get_f(pr)) % n;
      if (o)
      {
        GEN v = ZV_to_Flv(isprincipalray(bnr, pr), n);
        H = vec_append(H, Flv_Fl_mul(v, o, n));
      }
    }
    H = Flm_image(H, n); if (lg(cyc)-lg(H) <= k) break;
    moreorders(sN, CHI, F, pP, pO, bound);
  }
  H = hnfmodid(shallowconcat(zm_to_ZM(H), diagonal_shallow(cyc)), gn);

  Bquo = cgetg(k+1, t_MAT);
  for (i = j = 1; i <= d; i++)
    if (!equali1(gcoeff(H,i,i))) gel(Bquo,j++) = col_ei(d,i);

  for (i = 1, T = NULL; i<=k; i++)
  {
    GEN Hi = hnfmodid(shallowconcat(H, vecsplice(Bquo,i)), gn);
    GEN pol = rnfkummer(bnr, Hi, prec);
    T = T? nfcompositum(nf, T, pol, 2): pol;
  }
  T = rnfequation(nf, T); return gc_all(av, 3, &T, pP, pO);
}

static GEN
search_solvable(GEN LG, GEN mf, GEN F, long prec)
{
  GEN N = MF_get_gN(mf), CHI = MF_get_CHI(mf), pol, O, P, nf, Nfa;
  long i, l = lg(LG), v = fetch_var();
  ulong bound = 1;
  O = cgetg(1, t_VECSMALL); /* projective order of rho(Frob_p) */
  P = cgetg(1, t_VECSMALL);
  Nfa = Z_factor(N);
  pol = pol_x(v);
  for (i = 1; i < l; i++)
  { /* n prime, find a (Z/nZ)^k - extension */
    GEN G = gel(LG,i);
    long n = G[1], k = G[2];
    nf = nfinitred(mkvec2(pol,Nfa), prec);
    pol = search_abelian(nf, n, k, N, CHI, F, &P, &O, &bound, prec);
    setvarn(pol,v);
  }
  delete_var(); setvarn(pol,0); return pol;
}

static GEN
search_A5(GEN mf, GEN F)
{
  GEN CHI = MF_get_CHI(mf), O, P, L;
  long N = MF_get_N(mf), i, j, lL, nd, r;
  ulong bound = 1;
  r = radicalu(N);
  L = veccond_to_A5(zv_z_mul(divisorsu(N/r),r), 2); lL = lg(L); nd = lL-1;
  if (nd == 1) return gmael(L,1,1);
  O = cgetg(1, t_VECSMALL); /* projective order of rho(Frob_p) */
  P = cgetg(1, t_VECSMALL);
  for(i = 1; nd > 1; )
  {
    long l;
    moreorders(N, CHI, F, &P, &O, &bound);
    l = lg(P);
    for (   ; i < l; i++)
    {
      ulong p = P[i], f = O[i];
      for (j = 1; j < lL; j++)
        if (gel(L,j))
        {
          GEN FE = ZpX_primedec(gmael(L,j,1), utoi(p)), F = gel(FE,1);
          long nF = lg(F)-1;
          if (!equaliu(gel(F, nF), f)) { gel(L,j) = NULL; nd--; }
        }
      if (nd <= 1) break;
    }
  }
  for (j = 1; j < lL; j++)
    if (gel(L,j)) return gmael(L,j,1);
  return NULL;
}

GEN
mfgaloisprojrep(GEN mf, GEN F, long prec)
{
  pari_sp av = avma;
  GEN LG = NULL;
  if (!checkMF_i(mf) && !checkmf_i(F)) pari_err_TYPE("mfgaloisrep", F);
  switch( itos(mfgaloistype(mf,F)) )
  {
    case 0: case -12:
      LG = mkvec2(mkvecsmall2(3,1), mkvecsmall2(2,2)); break;
    case -24:
      LG = mkvec3(mkvecsmall2(2,1), mkvecsmall2(3,1), mkvecsmall2(2,2)); break;
    case -60: return gerepilecopy(av, search_A5(mf, F));
    default: pari_err_IMPL("mfgaloisprojrep for types D_n");
  }
  return gerepilecopy(av, search_solvable(LG, mf, F, prec));
}
