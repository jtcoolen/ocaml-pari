/* Copyright (C) 2020  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_nflist

/* Code s: if s >= 0 number of complex embeddings; s = -1: all signatures;
 * s = -2: all signatures separated.
 * Known groups: C1 = 1T1, C2 = 2T1, C3 = 3T1, S3 = 3T2, C4 = 4T1, V4 = 4T2,
 * D4 = 4T3, A4 = 4T4, S4 = 4T5, C5 = 5T1, D5 = 5T2, F5 = M20 = 5T3, A5 = 5T4,
 * C6 = 6T1, S36 = D66 = 6T2, D612 = 6T3, A46 = 6T4, S3C3 = 6T5,
 * A462 = 6T6, S46P = 6T7, S46M = 6T8, C32C4 = 6T10, S462 = 6T11,
 * A56 = PSL25 = 6T12, C32D4 = 6T13,
 * C7 = 7T1, D7 = 7T2, M21 = 7T3, M42 = 7T4, C9 = 9T1, C3C3 = 9T2,
 * CL and DL for L prime. A5cond is A5 ordered by twist-minimal conductor.
 *
 * For each group G:
 * - makeG(GEN N, GEN F, long s, long prec):
 * fields of given Galois group G, absolute discriminant N, signature s, and
 * auxiliary field (possibly NULL) F.
 * - makeGvec(GEN X, GEN Xinf, GEN F, long s, long prec):
 * fields of given Galois group G, absolute discriminant between Xinf and X,
 * signature s, and auxiliary field (possibly NULL) F.
 * - makeGresolvent(GEN pol, long flag)
 * - makeGsome(long s, long n, long prec): find n fields of given Galois group
 * G and signature s, not necessarily the smallest; n is assumed small. Useful
 * only when makeG/makeGvec take a long time. */

/* Relations between discriminants:
* D = disc of quadratic subfield or resolvent; Dk = disc subfield of degree k.
* f,g are integers, not necessarily conductors.
*
* C_ell: f^(ell-1), conductor f; ell odd prime
* D_ell: (Df^2)^((ell-1)/2), (f) = cond. over quad. subfield; ell odd prime
* C4: D^3f^2 (D sum of 2 squares); "conductor" Df
* V4: D^2f^2 (D for any of the 3 max. subfields); conductor lcm(D1,D2)
* D4: D^2f; "conductor" Df
* A4: Df^2, D = g^2
* S4: Df^2
* F5 = M20: Df^4 or 25D f^4 iff 125|D (what f's ?)
* C6: D^3D3^2/gcd(D,D3)^2; conductor lcm(D,D3)
* D6 = D6(12): D^3 D3^2 / gcd(D,D3)^2 * (1 or 4)
* S3(6) = D6(6): D^3 f^4 = D3^2 D
* A4(6), S4(6)+: D3 * D4
* S4(6)-: D*D3*D4 / gcd(g,D)^2*(1,4,16); D4 = D3*g^2, D3 = Df^2, D4 = D3g^2
*         D*D3*D4 = D^3f^4g^2 = D3^3(g/f)^2 = D3^2*D*g^2=D4^3/(g^4f^2)
*         Disc = D3^2 D * (g * (1,2,4) / gcd(g,D))^2
*         16 iff v(D2)=2 and v(g)=2 or 3.
*         4 iff v(g)=1 or (v(D2)=2 and v(g)=4)
*         1 or 4 if v(D2)=3 and v(g)=4.
* C32C4: D D4 f^2
* S32: disc of 2 S3 subfields: D1F1^2, D2F2^2, D1, D2 fund. disc
*      disc = (D1D2)^3 / gcd(D1,D2)^4 * lcm(F1,F2)^2 * g^2
* M21: D^2f^6 (D = g^2) or 7^4 D^2 f^6 iff 49|D
* M42: Df^6 or 7^2 Df^6 iff 7^4|D or 7^4 Df^6 iff 7^5|D
* C9: D^4f^6
* C3xC3: lcm(D3, D3')^3
* D9: D^4 g^2 f^6 (disc subic subfield = Dg^2)
*
* Minimimal discriminants for each group, by s
* C1: [1]
* C2: [5, 3]
* C3: [49, 0]
* S3: [148, 23]
* C4: [1125, 0, 125]
* V4: [1600, 0, 144]
* D4: [725, 275, 117]
* A4: [26569, 0, 3136]
* S4: [1957, 283, 229]
* C5: [14641, 0, 0]
* D5: [160801, 0, 2209]
* F5: [2382032, 0, 35152]
* A5: [3104644, 0, 18496]
* C6: [300125, 0, 0, 16807]
* S36: [810448, 0, 0, 12167]
* D612: [2738000, 0, 66125, 14283]
* A46: [25969216, 0, 153664, 0]
* S3C3: [722000, 0, 0, 9747]
* A462: [434581, 103243, 31213, 0]
* S46+: [3356224, 0, 33856, 0]
* S46-: [7495014493, 0, 3241792, 778688]
* S32: [27848000, 0, 242000, 309123]
* C32C4: [55130625, 0, 525625, 0]
* S462: [1387029, 309123, 28037, 10051]
* C7: [594823321, 0, 0]
* D7: [192100033, 0, 0, 357911]
* M21: [1817487424, 0, 0, 0]
* M42: [12431698517, 0, 0, 38014691]
* C9: [16983563041, 0, 0, 0, 0]
* C3C3: [62523502209, 0, 0, 0, 0]
* D9: [1624709678881, 0, 0, 0, 775511104]
* C11: [41426511213649, 0, 0, 0, 0]
* D11: [3670285774226257, 0, 0, 0, 0, 129891985607] */

/* FIXME: export */
static long
RgVV_nb(GEN v)
{
  long i, l = lg(v), n = 0;
  for (i = 1; i < l; i++) n += lg(gel(v,i)) - 1;
  return n;
}
/* FIXME: export */
static GEN
gtoset_shallow(GEN x)
{
  GEN p = gen_indexsort_uniq(x, (void*)&cmp_universal, cmp_nodata);
  return vecpermute(x, p);
}

static GEN
nflist_parapply(const char *s, GEN v, GEN w)
{
  GEN L;
  if (DEBUGLEVEL>=3) err_printf("%s: ",s);
  L = gen_parapply_percent(snm_closure(is_entry(s), v), w, DEBUGLEVEL>=3);
  if (DEBUGLEVEL>=3) err_printf("done\n");
  return L;
}

/**************************************************************************/
/*                        Utility functions                               */
/**************************************************************************/
static long
divissquareall(GEN x, GEN y, GEN *f)
{ GEN r, q = dvmdii(x, y, &r); return r == gen_0 && Z_issquareall(q,f); }
static long
divissquare(GEN x, GEN y) { return divissquareall(x, y, NULL); }
static long
divispowerall(GEN x, GEN y, ulong k, GEN *f)
{ GEN r, q = dvmdii(x, y, &r); return r == gen_0 && Z_ispowerall(q,k,f); }
/* x / y if y | x, else NULL */
static GEN
divide(GEN x, GEN y)
{ GEN r, q = dvmdii(x, y, &r); return r == gen_0? q: NULL; }

/* ceil(X^(1/n)) */
static long
ceilsqrtn(GEN X, long n)
{
  pari_sp av = avma;
  ulong x = itou(sqrtnint(X, n));
  if (cmpii(powuu(x, n), X) < 0) x++;
  return gc_long(av, x);
}
static long
ceilsqrt(GEN X)
{
  pari_sp av = avma;
  GEN r;
  ulong x = itou(sqrtremi(X, &r));
  return gc_long(av, r==gen_0? x: x+1);
}
static GEN
gceilsqrtn(GEN X, long n)
{
  GEN x = sqrtnint(X, n);
  if (cmpii(powiu(x, n), X) < 0) x = addiu(x, 1);
  return x;
}
/* assume X >= 0 or n odd */
static long
sceilsqrtn(long X, long n)
{
  ulong x, Xa;
  if (!X) return 0;
  Xa = labs(X); x = usqrtn(Xa, n);
  if (X > 0 && upowuu(x, n) != Xa) x++;
  return X > 0? (long)x: -(long)x;
}
/* ceil((X/Y)^1/n)*/
static long
ceilsqrtndiv(GEN X, GEN Y, long n)
{
  pari_sp av = avma;
  ulong x = itou(sqrtnint(divii(X, Y), n));
  if (cmpii(mulii(powuu(x, n), Y), X) < 0) x++;
  return gc_long(av, x);
}
long
ceilsqrtdiv(GEN X, GEN Y)
{
  pari_sp av = avma;
  GEN r, q = dvmdii(X, Y, &r);
  ulong x = itou((r == gen_0)? sqrtremi(q, &r): sqrti(q));
  return gc_long(av, r==gen_0? x: x+1);
}
static GEN
gceilsqrtdiv(GEN X, GEN Y)
{
  GEN r, q = dvmdii(X, Y, &r);
  q = (r == gen_0)? sqrtremi(q, &r): sqrti(q);
  return r == gen_0? q: addiu(q, 1);
}
static GEN
gfloorsqrtdiv(GEN X, GEN Y) { return sqrti(divii(X, Y)); }
/* floor(X^(1/n)) */
static long
floorsqrtn(GEN X, long n)
{ pari_sp av = avma; return gc_long(av, itou(sqrtnint(X, n))); }
static long
floorsqrt(GEN X)
{ pari_sp av = avma; return gc_long(av, itou(sqrti(X))); }
/* floor((X/Y)^(1/n)) */
static long
floorsqrtndiv(GEN X, GEN Y, long n)
{ pari_sp av = avma; return gc_long(av, itou(sqrtnint(divii(X,Y), n))); }
static long
floorsqrtdiv(GEN X, GEN Y)
{ pari_sp av = avma; return gc_long(av, itou(gfloorsqrtdiv(X, Y))); }
static GEN
ceildiv(GEN X, GEN Y)
{
  GEN r, q = dvmdii(X, Y, & r);
  return (r == gen_0)? q: addiu(q, 1);
}

static GEN
nfY(GEN T)
{ T = shallowcopy(T); setvarn(T,1); return nfinit(T, DEFAULTPREC); }
static GEN
bnfY(GEN T)
{ T = shallowcopy(T); setvarn(T,1); return Buchall(T, nf_FORCE, DEFAULTPREC); }
static GEN
bnf_get_disc(GEN b) { return nf_get_disc(bnf_get_nf(b)); }

/* Compute n s.t. d | n <=> d^k | N. Return [n, factor(n)] */
static GEN
cored(GEN N, long k)
{
  GEN fa = Z_factor(N), P = gel(fa,1), E = gel(fa,2), n = gen_1;
  long i, c, l = lg(P);

  for (i = c = 1; i < l; i++)
  {
    long e = itou(gel(E,i));
    if (e >= k)
    {
      e /= k; n = mulii(n, powiu(gel(P,i), e));
      gel(P,c) = gel(P,i); gel(E,c) = utoipos(e); c++;
    }
  }
  setlg(P,c); setlg(E,c); return mkvec2(n, fa);
}

/* return D = nfdisc(T), set d = coredisc */
static GEN
nfcoredisc(GEN T, GEN *pd)
{
  GEN D = nfdiscfactors(T), d = core(D); /* d = core(|D|) */
  D = gel(D,1); if (signe(D) < 0) togglesign_safe(&d);
  if (Mod4(d) != 1) d = shifti(d,2); /* = coredisc(D) */
  *pd = d; return D;
}
static GEN
nfcoredisc2(GEN T, GEN *pd, GEN *pf)
{
  GEN D = nfcoredisc(T, pd);
  if (pf) *pf = sqrti(diviiexact(D, *pd));
  return D;
}

/* \prod {pr | ell} pr */
static GEN
getpell(GEN nf, long ell, long *pteell)
{
  GEN P = idealprimedec(nf, utoipos(ell));
  *pteell = pr_get_e(gel(P,1)); return idealfactorback(nf, P, NULL, 0);
}

static void
checkfield_i(GEN F, long d)
{ if (F && degpol(F) != d) pari_err_TYPE("nflist", F); }
static GEN
checkfield(GEN F, long d) { checkfield_i(F, d); return nfdisc(F); }

static long
pol2s(GEN T) { return (degpol(T) - ZX_sturm_irred(T)) >> 1; }

static GEN
sturmseparate(GEN V, long s, long deg)
{
  GEN w, C;
  long l, ls , i;

  if (s != -2) return V;
  l = lg(V); ls = (deg >> 1) + 2;
  w = cgetg(ls, t_VEC);
  C = cgetg(ls, t_VECSMALL);
  for (i = 1; i < ls; i++) { gel(w, i) = cgetg(l, t_VEC); C[i] = 1; }
  for (i = 1; i < l; i++)
  {
    long k = pol2s(gel(V, i)) + 1;
    gmael(w, k, C[k]++) = gel(V, i);
  }
  for (i = 1; i < ls; i++) setlg(gel(w, i), C[i]);
  return w;
}

/* fa = factorization of positive integer N. Are +N and/or -N fundamental ? */
static void
fa_is_fundamental_pm(GEN N, GEN fa, long s, int *p, int *m)
{
  GEN P = gel(fa,1), E = gel(fa,2);
  long l = lg(P), i;
  ulong r, r4;

  if (l == 1) { *m = 0; *p = (s <= 0); return; }
  r = Mod16(N); r4 = r & 3UL;
  if (!r || r4 == 2) { *p = *m = 0; return; } /* v_2 > 3 or N=2 mod 4 */
  *p = (s <= 0);
  *m = s? 1: 0;
  if (odd(r))
  {
    if (r4 == 1) { *m = 0; if (!*p) return; }
    else         { *p = 0; if (!*m) return; }
    i = 1;
  }
  else
  { /* P[1] = 2 => 4 | N */
    if (r == 4)       { *p = 0; if (!*m) return; }
    else if (r == 12) { *m = 0; if (!*p) return; }
    i = 2;
  }
  for (; i < l; i++)
    if (itou(gel(E,i)) > 1) { *p = *m = 0; return; }
}
/* if flag is set assume the odd part of N is squarefree */
static void
uis_fundamental_pm_i(ulong N, long s, int *p, int *m, long flag)
{
  ulong r, r4;

  if (N == 1UL) { *m = 0; *p = (s <= 0); return; }
  r = N & 15UL; r4 = r & 3UL;
  if (!r || r4 == 2) { *p = *m = 0; return; } /* v_2 > 3 or N=2 mod 4 */
  *p = (s <= 0);
  *m = s? 1: 0;
  if (odd(r))
  {
    if (r4 == 1) { *m = 0; if (!*p) return; }
    else         { *p = 0; if (!*m) return; }
  }
  else
  { /* P[1] = 2 => 4 | N */
    if (r == 4)       { *p = 0; if (!*m) return; }
    else if (r == 12) { *m = 0; if (!*p) return; }
    N >>= (r == 8? 3: 2); /* odd part */
  }
  if (!flag && !uissquarefree(N)) { *p = *m = 0; }
}
static void
uis_fundamental_pm(ulong N, long s, int *p, int *m)
{ uis_fundamental_pm_i(N, s, p, m, 0); }

static void
is_fundamental_pm(GEN N, long s, int *p, int *m)
{
  ulong r, r4;

  if (lgefint(N) == 3) { uis_fundamental_pm(N[2], s, p, m); return; }
  r = Mod16(N); r4 = r & 3UL;
  if (!r || r4 == 2) { *p = *m = 0; return; } /* v_2 > 3 or N=2 mod 4 */
  *p = (s <= 0);
  *m = s? 1: 0;
  if (odd(r))
  {
    if (r4 == 1) { *m = 0; if (!*p) return; }
    else         { *p = 0; if (!*m) return; }
  }
  else
  { /* P[1] = 2 => 4 | N */
    if (r == 4)       { *p = 0; if (!*m) return; }
    else if (r == 12) { *m = 0; if (!*p) return; }
    N = shifti(N, r == 8? -3: -2); /* odd part */
  }
  if (!Z_issquarefree(N)) { *p = *m = 0; }
}
static GEN
fund_pm(GEN N, int p, int m)
{
  if (p && m) return mkvec2(N, negi(N));
  if (p) return mkvec(N);
  if (m) return mkvec(negi(N));
  return NULL;
}
static GEN
ufund_pm(ulong N, int p, int m)
{
  if (p && m) return mkvec2(utoipos(N), utoineg(N));
  if (p) return mkvec(utoipos(N));
  if (m) return mkvec(utoineg(N));
  return NULL;
}

static GEN
divisorsdisc(GEN N, long s)
{
  GEN D, V;
  long l, c = 1, i;

  if (typ(N) == t_VEC)
  { /* [n, factor(n)]; assume n > 0 */
    GEN n = gel(N,1), fa = gel(N,2);
    if (Mod4(n) == 2) N = mkvec2(shifti(n,-1), rowsplice(fa, 1));
  }
  else
    if (Mod4(N) == 2) N = shifti(N, -1);
  D = divisors_factored(N); l = lg(D);
  V = cgetg(2 * l - 1, t_VEC);
  for (i = 2; i < l; i++)
  {
    GEN d = gel(D, i);
    int p, m;
    fa_is_fundamental_pm(gel(d,1), gel(d,2), s, &p, &m);
    if (p) gel(V, c++) = gel(d,1);
    if (m) gel(V, c++) = negi(gel(d,1));
  }
  setlg(V, c); return V;
}

static int
usum2sq(ulong m)
{
  pari_sp av = avma;
  GEN fa, P, E;
  long i, v2 = vals(m);
  if (v2)
  {
    if (v2 != 3) return 0;
    m >>= 3;
  }
  if ((m & 3L) != 1) return 0;
  fa = factoru(m); P = gel(fa, 1); E = gel(fa, 2);
  for (i = 1; i < lg(P); i++)
    if (E[i] >= 2 || (P[i] & 3L) == 3) { set_avma(av); return 0; }
  set_avma(av); return 1;
}
static int
sum2sq(GEN m)
{
  pari_sp av = avma;
  GEN fa, P, E;
  long i, v2;
  if (lgefint(m) == 3) return usum2sq(m[2]);
  v2 = vali(m);
  if (v2)
  {
    if (v2 != 3) return 0;
    m = shifti(m, -3);
  }
  if (Mod4(m) != 1) return 0;
  fa = Z_factor(m); P = gel(fa, 1); E = gel(fa, 2);
  for (i = 1; i < lg(P); i++)
    if (!equali1(gel(E,i)) || Mod4(gel(P,i)) == 3) { set_avma(av); return 0; }
  set_avma(av); return 1;
}

static int
ok_int(GEN d, GEN X, GEN Xinf)
{ return (abscmpii(d, X) <= 0 && abscmpii(d, Xinf) >= 0); }
static int
ok_intu(GEN d, ulong X, ulong Xinf)
{ return (abscmpiu(d, X) <= 0 && abscmpiu(d, Xinf) >= 0); }

static int
ok_disc(GEN d, GEN X, GEN Xinf)
{
  if (!Xinf) return absequalii(d, X);
  return ok_int(d, X, Xinf);
}

/* G cyclic galoisinit */
static GEN
cyclicgalois(GEN bnr, GEN G, long *o)
{
  GEN g = galoispermtopol(G, gel(gal_get_gen(G), 1));
  *o = gal_get_orders(G)[1];
  return bnrautmatrix(bnr, g); /* order o */
}
/* Cl_f / H cyclic of prime order, return i s.t bnr.cyc[i] is generator */
static long
cyclicprimegen(GEN H)
{
  long i, l = lg(H);
  for (i = 1; i < l; i++) if (!is_pm1(gcoeff(H,i,i))) return i;
  return -1;/*LCOV_EXCL_LINE*/
}
/* k/Q cyclic and M the bnrautmatrix for the generator s of its Galois group
 * (action on bnr = Cl_f(k)). vH a vector of congruence subgroups for bnr,
 * attached to abelian extensions K/k of prime degree, assumed to be Galois
 * over Q [sf = f and sH = H]. Filter out the H corresponding to K/Q abelian */
static void
nonabelianfilter(GEN vH, GEN M)
{
  long i, c, l = lg(vH);
  for (i = c = 1; i < l; i++)
  {
    GEN v, H = gel(vH,i);
    long k = cyclicprimegen(H);
    v = shallowcopy(gel(M,k));
    gel(v,k) = subiu(gel(v,k), 1);
    if (!hnf_invimage(H, v)) gel(vH, c++) = H;
  }
  setlg(vH, c);
}


/* bnf attached to K. Cyclic extensions L/K of degree d and exact conductor
 * F; if F = [F,Finf]~, check that Finf | conductor | F;
 * check that |disc L/Q| in [Xinf,X] if not NULL. If G != NULL,
 * then K/Q = <s> is cyclic, we assume s.F = F and
 * G = [galoisinit(bnf), flag], with flag > 0 (resp. 0) to insist L be
 * Galois / Q (resp. not Galois). If flag = 2, insist that L/Q is non abelian.
 * In the non-Galois case, keep only one among isomorphic extensions attached
 * to sigma.H; sigma in Gal(K/Q). For simplicity assume the base is cyclic;
 * will extend it later if needed. */
static GEN
mybnrclassfield_X(GEN bnf, GEN F, long d, GEN X, GEN Xinf, GEN G)
{
  GEN gd = utoipos(d), Finf = NULL, bnr, L;
  long i, j, c, l;

  if (typ(F) == t_COL) { Finf = gel(F,1); F = gel(F,2); }
  bnr = bnrinitmod(bnf, F, 0, gd);
  L = subgrouplist0(bnr, mkvec(gd), Finf? 1: 0); l = lg(L);
  if (Finf)
  {
    GEN Fi = idealinv(bnr, Finf);
    for (i = c = 1; i < l; i++)
    { /* for now assume that F and Finf are finite */
      GEN f = gel(bnrconductor_raw(bnr, gel(L,i)), 1);
      if (equali1(Q_denom(idealmul(bnr, f, Fi)))) gel(L,c++) = gel(L,i);
    }
    setlg(L, c); l = c;
  }
  if (l == 1) return L;
  if (!uisprime(d))
  {
    for (i = j = 1; i < l; i++)
      if (lg(smithclean(ZM_snf(gel(L,i)))) == 2) gel(L,j++) = gel(L,i);
    setlg(L, l = j); if (l == 1) return L;
  }
  if (G)
  {
    GEN M;
    long o, gal = itou(gel(G,2));
    if (l == 2)
    {  /* => L[1] is fixed: must be Galois */
      if (!gal) { setlg(L,1); return L; }
      if (gal == 2)
      {
        M = cyclicgalois(bnr, gel(G,1), &o);
        nonabelianfilter(L, M);
      }
    }
    else
    {
      M = cyclicgalois(bnr, gel(G,1), &o); /* assume cyclic for now */
      if (gal)
      {
        for (i = j = 1; i < l; i++)
        {
          GEN H = gel(L,i);
          if (ZM_equal(bnrgaloisapply(bnr, M, H), H)) gel(L,j++) = H;
        }
        setlg(L, l = j);
        if (gal == 2) nonabelianfilter(L, M);
      }
      else
      {
        for (i = 1; i < l; i++)
        {
          GEN H = gel(L,i), K = bnrgaloisapply(bnr, M, H);
          long k;

          /* \sigma H = H <=> Galois : delete */
          if (ZM_equal(K, H)) { L = vecsplice(L,i--); l--; continue; }
          /* else delete the rest of Galois orbit */
          for (j = 1; j < o; j++)
          {
            for (k = i+1; k < l; k++)
              if (ZM_equal(K, gel(L,k))) { L = vecsplice(L,k); l--; break; }
            if (j != o-1) K = bnrgaloisapply(bnr, M, K);
          }
        }
      }
    }
    if ((l = lg(L)) == 1) return L;
  }
  if (X)
  {
    for (i = j = 1; i < l; i++)
    {
      GEN D = gel(bnrdisc(bnr, gel(L,i), 0), 3);
      if (ok_disc(D, X, Xinf)) gel(L,j++) = gel(L,i);
    }
    setlg(L, j); if (j == 1) return L;
  }
  return shallowconcat1(bnrclassfield(bnr, L, 0, DEFAULTPREC));
}
static GEN
mybnrclassfield_N(GEN bnf, GEN F, GEN N, long d)
{ return mybnrclassfield_X(bnf, F, d, N, NULL, NULL); }
static GEN
mybnrclassfield(GEN bnf, GEN F, long d)
{ return mybnrclassfield_X(bnf, F, d, NULL, NULL, NULL); }

/* N > 1 */
static int
checkcondell_i(GEN N, long ell, GEN D2, GEN *pP)
{
  GEN fa, P, E;
  long l, i, e;

  if (typ(N) == t_VEC)
  {
    fa = gel(N,2); P = gel(fa, 1); E = gel(fa, 2);
    i = ZV_search(P, utoipos(ell));
    if (!i) e = 0;
    else
    {
      e = itou(gel(E,i)); if (e != 2) return 0;
      P = vecsplice(P, i);
      E = vecsplice(E, i);
    }
  }
  else
  {
    e = Z_lvalrem(N, ell, &N);
    if (e != 0 && e != 2) return 0;
    fa = Z_factor(N); P = gel(fa, 1); E = gel(fa, 2);
  }
  l = lg(P);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(P,i);
    long r;
    if (!equaliu(gel(E,i), 1)) return 0;
    r = umodiu(p, ell);
    if (!D2) { if (r != 1) return 0; }
    else
    {
      r -= kronecker(D2, p);
      if (r && r != ell) return 0;
    }
  }
  *pP = P; return 1;
}
/* ell odd prime, N potential conductor for C_ell field, *pP contains
 * the prime divisors of N different from ell */
static int
checkcondCL(GEN N, long ell, GEN *pP)
{ GEN n = typ(N) == t_VEC? gel(N, 1): N;
  return odd(Mod4(n)) && !equali1(n) && checkcondell_i(N, ell, NULL, pP); }
/* D2 fundamental discriminant, ell odd prime, N potential conductor for
 * D_ell field over Q(sqrt(D2)) */
static int
checkcondDL(GEN D2, GEN N, long ell, GEN *pP)
{
  ulong N4;
  if (!umodiu(D2, ell))
  {
    long v = Z_lvalrem(N, ell, &N);
    if (v && v > 2) return 0;
  }
  if (equali1(N)) { *pP = cgetg(1,t_VEC); return 1; }
  N4 = Mod4(N);
  return N4 && (N4 != 2 || ell == 3) && checkcondell_i(N, ell, D2, pP);
}

static GEN
myshallowconcat1(GEN V)
{
  if (lg(V) == 1) return V;
  return shallowconcat1(V);
}

static GEN
_nfsubfields(GEN pol, long d) { return nfsubfields0(pol, d, 1); }
static GEN
_nfsubfields1(GEN pol, long d) { return gel(_nfsubfields(pol, d), 1); }
static GEN
mynfsubfields(GEN pol, long d)
{
  GEN V = _nfsubfields(pol, d), W;
  long l = lg(V), i;
  W = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(W,i) = polredabs(gel(V,i));
  return W;
}
static GEN
mynfsubfield(GEN pol, long d)
{
  if (d == 2 && (degpol(pol) & 3) == 2)
    return quadpoly_i(quaddisc(ZX_disc(pol)));
  return polredabs(gel(_nfsubfields(pol, d), 1));
}

/* global checks to be done:
-- in nflist: if s > deg / 2, return empty.
-- in nfresolvent: check polynomial of correct degree.
*/

/***************************************************************/

static GEN
makeC1(GEN N, GEN field, long s)
{
  checkfield_i(field, 1);
  if (!equali1(N)) return NULL;
  return mkvec(s != -2? pol_x(0): mkvec(pol_x(0)));
}
static GEN
makeC1resolvent(long flag)
{ return odd(flag)? mkvec2(pol_x(0), gen_1): pol_x(0); }
static GEN
makeC1vec(GEN Xinf, GEN field, long s) { return makeC1(Xinf, field, s); }

/**********************************************************************/
/*                                 C2                                 */
/**********************************************************************/
static GEN
makeC2(GEN N, GEN field, long s)
{
  GEN V = NULL;
  long l, i;
  int p, m;

  checkfield_i(field, 1);
  if (equali1(N) || Mod4(N) == 2) return NULL;
  is_fundamental_pm(N, s, &p, &m);
  if (!(V = fund_pm(N, p, m))) return NULL;
  l = lg(V);
  for (i = 1; i < l; i++) gel(V, i) = quadpoly_i(gel(V, i));
  return sturmseparate(V, s, 2);
}

static GEN
makeC2resolvent(GEN pol, long flag)
{ return odd(flag)? mkvec2(pol_x(0), absi_shallow(nfdisc(pol))): pol_x(0); }

static GEN
makeC2vec(GEN X, GEN Xinf, GEN field, long s)
{
  long M, cv, cw, l = itou(subii(X, Xinf)) + 1;
  GEN v, w;

  checkfield_i(field, 1);
  v = (s <= 0)? cgetg(l, t_VEC): NULL;
  w = s? cgetg(l, t_VEC): NULL;
  for (M = equali1(Xinf)? 2: 1, cv = cw = 1; M < l; M++)
  {
    GEN N = addiu(Xinf, M);
    int p, m;
    is_fundamental_pm(N, s, &p, &m);
    if (p) gel(v, cv++) = quadpoly_i(N);
    if (m) gel(w, cw++) = quadpoly_i(negi(N));
  }
  if (cv == 1 && cw == 1) return NULL;
  switch (s)
  {
    case 0:  setlg(v, cv); return v;
    case 1:  setlg(w, cw); return w;
    case -1: setlg(v, cv); setlg(w, cw); return shallowconcat(v, w);
    default: setlg(v, cv); setlg(w, cw); return mkvec2(v, w);
  }
}

/**********************************************************************/
/*                                 C3                                 */
/**********************************************************************/
/* \prod x[i]^e[i], e[i] in {0,1} */
static GEN
eltlist2(GEN nf, GEN x)
{
  long i, j, c, l = lg(x);
  GEN v;
  if (l == 1) return mkvec(gen_1);
  v = cgetg((1 << (l-1))+1, t_VEC);
  gel(v,1) = gen_1;
  gel(v,2) = gel(x,1);
  for (i = c = 2; i < l; i++, c <<= 1)
    for (j = 1; j <= c; j++) gel(v, c + j) = nfmul(nf, gel(v,j), gel(x,i));
  return v;
}
/* { x[1][1] * \prod_i>=2 x[i][e_i], (e) in {1,2}^(#x-1)} */
static GEN
mullist2(GEN x)
{
  long i, j, c, l = lg(x);
  GEN v;
  if (l == 2) return mkvec(gmael(x,1,1));
  v = cgetg((1 << (l-2))+1, t_VEC);
  gel(v,1) = gel(v,2) = gmael(x,1,1);
  for (i = 2, c = 1; i < l; i++, c <<= 1)
    for (j = 1; j <= c; j++)
    {
      gel(v, c + j) = gmul(gel(v, j), gmael(x,i,2));
      gel(v, j) = gmul(gel(v, j), gmael(x,i,1));
    }
  return v;
}

static GEN
makepolC3(GEN n, GEN u, long fl3)
{
  GEN T = cgetg(6, t_POL), n3, nu27;
  T[1] = evalsigne(1) | evalvarn(0);
  gel(T, 5) = gen_1;
  gel(T, 4) = fl3 ? gen_m1 : gen_0;
  if (!fl3)
  { n3 = divis(n, -3); nu27 = mulii(n, u); }
  else
  { n3 = divis(subiu(n, 1), -3); nu27 = addiu(mulii(n, subiu(u, 3)), 1); }
  gel(T, 3) = n3;
  gel(T, 2) = divis(nu27, -27); return T;
}

static GEN
decp(GEN Q, GEN t, GEN p)
{
  GEN u, v, z;
  if (equaliu(p, 3)) { u = utoineg(3); v = utoipos(3); }
  else
  {
    GEN uv = qfbsolve(Q, shifti(p, 2), 2);
    u = gel(uv,1); if (umodiu(u, 3) == 1) togglesign(u);
    v = muliu(gel(uv,2), 3); if (signe(v) < 0) togglesign(v);
  }
  z = gadd(gmul(v, t), shifti(subii(u, v), -1));
  return mkvec2(z, conj_i(z));
}

static int
checkcondC3(GEN n, GEN *pP)
{
  GEN fa = NULL, P, E;
  long l, i, n27;

  *pP = NULL;
  if (typ(n) == t_VEC) { fa = gel(n,2); n = gel(n,1); }
  if (cmpiu(n, 7) < 0 || !mpodd(n)) return 0;
  n27 = umodiu(n, 27);
  switch(n27 % 3)
  {
    case 2: return 0;
    case 1: i = 1; break;
    default: i = 2; if (n27 != 9 && n27 != 18) return 0;
  }
  if (!fa) fa = Z_factor(n);
  P = gel(fa, 1); E = gel(fa, 2); l = lg(P);
  for (; i < l; i++)
    if (umodiu(gel(P,i), 3) != 1 || !equali1(gel(E,i))) return 0;
  *pP = P; return 1;
}

static GEN
makeC3_i(GEN sqN, GEN P)
{
  GEN v, t, Q = mkqfb(gen_1, gen_0, utoipos(27), utoineg(108));
  long i, j, l, n = lg(P)-1, fl3 = umodiu(gel(P,1), 3);

  t = quadgen0(utoineg(3), 1); v = cgetg(n+1, t_VEC);
  for (i = 1; i <= n; i++) gel(v,i) = decp(Q, t, gel(P,i));
  v = mullist2(v); l = lg(v);
  for (j = 1; j < l; j++) gel(v,j) = makepolC3(sqN, gtrace(gel(v,j)), fl3);
  return v;
}
/* makeC3(f^2, 0) */
static GEN
makeC3_f(GEN f)
{
  GEN P;
  return checkcondC3(f, &P)? makeC3_i(f, P): cgetg(1, t_VEC);
}
static GEN
vecs(long ns, GEN x)
{ GEN v = const_vec(ns, cgetg(1,t_VEC)); gel(v,1) = x; return v; }
static GEN
vecs14(GEN x, GEN y) { GEN v = cgetg(1,t_VEC); return mkvec4(x,v,v,y); }

static GEN
makeC3(GEN N, GEN field, long s)
{
  GEN v, f, P;

  checkfield_i(field, 1);
  if (s > 0 || cmpiu(N, 49) < 0 || !Z_issquareall(N, &f)
      || !checkcondC3(f, &P)) return NULL;
  v = makeC3_i(f, P); return s == -2 ? vecs(2, v): v;
}

static GEN
makeC3resolvent(GEN pol, long flag)
{ return odd(flag)? mkvec2(pol_x(0), sqrti(nfdisc(pol))): pol_x(0); }

GEN
nflist_C3_worker(GEN gv, GEN T)
{
  long v = itos(gv), sX = T[1], sXinf = T[2], c, r, u;
  long v227 = 27 * v * v, limu = usqrt((sX << 2) - v227);
  GEN V = cgetg(limu + 2, t_VEC);

  if (odd(limu - v)) limu--; /* make sure u = v (mod 2) */
  for (u = -limu, r = smodss(u, 9), c = 1; u <= limu; u += 2, r += 2)
  {
    if (r >= 9) r -= 9; /* r = u % 9 */
    if (r == 2 || r == 5 || r == 6 || r == 8) /* u = 2 (mod 3) or 6 (mod 9) */
    {
      long e;
      if (ugcd(labs(u), v) > 2) continue;
      e = (u * u + v227) >> 2; /* conductor, disc = e^2 */
      if (e < sXinf) continue;
      if (r == 6) e /= 9; /* 9 | e */
      if (!uissquarefree(e)) continue;
      gel(V, c++) = r==6? mkvecsmall4(1, 0, -3 * e, -e * u / 3)
                        : mkvecsmall4(1, -1, (1-e) / 3, -(1 + e * (u-3)) / 27 );
    }
  }
  setlg(V, c); return V;
}

static GEN
zvV_to_ZXV(GEN v)
{
  long i, l = lg(v);
  GEN w = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(w,i) = gtopoly(gel(v,i), 0);
  return w;
}
static GEN
C3vec(GEN V, long s)
{
  if (s != -2) return zvV_to_ZXV(V);
  retmkvec2(zvV_to_ZXV(V), cgetg(1,t_VEC));
}

/* t a C3 t_VECSMALL generated by C3_worker. Return its conductor f */
static long
uC3pol_f(GEN t) { return - t[2] - 3 * t[3]; }
/* t a C3 t_POL = gtopoly(C3_worker t_VECSMALL) */
static GEN
C3pol_f(GEN t) { return subii(mulsi(-3, gel(t,3)), gel(t,4)); }
/* C3vec for discriminant f^2, f in [sX,sXinf] */
static GEN
C3vec_F(long sX, long sXinf, GEN *pF)
{
  GEN v, F, perm, T = mkvecsmall2(sX, sXinf);
  long i, l, lim = usqrt((sX << 2) / 27);
  v = nflist_parapply("_nflist_C3_worker", mkvec(T), identity_ZV(lim));
  v = myshallowconcat1(v); l = lg(v); if (l == 1) return NULL;
  F = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) F[i] = uC3pol_f(gel(v,i));
  perm = vecsmall_indexsort(F);
  if (pF) *pF = vecsmallpermute(F, perm);
  return vecpermute(v, perm);
}
static GEN
makeC3vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v;
  checkfield_i(field, 1);
  if (s > 0 || !(v = C3vec_F(floorsqrt(X), ceilsqrt(Xinf), NULL))) return NULL;
  return C3vec(v, s);
}

/**********************************************************************/
/*                                 S3                                 */
/**********************************************************************/
/* Quadratic resolvent field. */

static GEN makeDL(long ell, GEN N, GEN field, long s);
static GEN makeDLvec(long ell, GEN X, GEN Xinf, GEN field, long s);

/* Cubic programs from KB and HC */
#define min(a, b) ((a) >= (b) ? b : a)
#define max(a, b) ((a) >= (b) ? a : b)

static GEN
checkU(long a, long b, long c, long d, long P, long Q, long R, long D)
{
  long t, f = cgcd(cgcd(P, Q), R);
  GEN F;

  if (odd(f)) { long e = D & 15L; if (e == 0 || e == 12) return NULL; }
  else if ((D & 7L) == 0) return NULL;
  if (f % 3 == 0)
  {
    if ((a % 9 == 0) || (a % 3 && (d % 9 == 0))) return NULL;
    if ((a % 3) && (d % 3))
    {
      long e = (a - d) % 3 ? - 1 : 1;
      if ((a + c - e * (b + d)) % 9 == 0) return NULL;
    }
    if (!uissquarefree(f / 9)) return NULL;
  }
  else if (D % 27 == 0 || !uissquarefree(f)) return NULL;
  t = labs(D) / (f * f); t >>= vals(t); while (t % 3 == 0) t /= 3;
  if (cgcd(t, f) > 1 || !uissquarefree(t)) return NULL;
  F = cgetg(6, t_POL); F[1] = evalsigne(1)|evalvarn(0);
  gel(F,2) = stoi(d * a * a);
  gel(F,3) = stoi(c * a);
  gel(F,4) = stoi(b);
  gel(F,5) = gen_1; return F;
}

/* ceil(m/d), assume d != 0 */
static long
sceildiv(long m, long d)
{
  long q;
  if (d == 1) return m;
  if (!m) return 0;
  if (d < 0) { d = -d; m = -m; }
  if (m < 0) return -((-m) / d);
  q = m / d; return m%d? q+1: q;
}
/* floor(m/d), assume d != 0 */
static long
sfloordiv(long m, long d)
{
  long q;
  if (d == 1) return m;
  if (!m) return 0;
  if (d < 0) { d = -d; m = -m; }
  if (m > 0) return m / d;
  q = -((-m) / d); return (-m)%d? q-1: q;
}

GEN
nflist_S3R_worker(GEN ga, GEN S)
{
  long a = itos(ga), a3 = 3 * a, a9 = 9 * a, b, c, d, ct = 1;
  long x = S[1], xinf = S[2], sqx = S[3], cplus = S[4], cminus = S[5];
  long cmin = S[6], Dmin = S[7], Dsup = S[8], bsup = S[9], binf = S[10];
  long csupa = usqrtn(cplus / a, 3), cinfa = sceilsqrtn(sceildiv(cminus, a), 3);
  long dsupa = Dsup / a, dinfa = sceildiv(Dmin, a);
  GEN RET = cgetg(x / 3, t_VEC);

  for (b = binf; b <= bsup; b++)
  {
    long cinf = cinfa, csup = csupa, dinfb = dinfa, dsupb = dsupa;
    long bb = b * b, b3 = 3 * b, gcdab = cgcd(a, b);
    if (b)
    {
      long bbb = bb * b, sqxb = sqx / labs(b), m, M;
      if (b < 0)
      {
        cinf = -sqxb; csup = -1;
        M = sfloordiv(cminus,bbb);
        m = sceildiv(cplus, bbb);
      }
      else
      {
        cinf = cmin; csup = minss(csup, sqxb);
        M = cplus / bbb;
        m = sceildiv(cminus, bbb);
      }
      dsupb = minss(dsupb, M);
      dinfb = maxss(dinfb, m); cinf = maxss(cinfa, cinf);
    }
    for (c = cinf; c <= csup; c++)
    {
      long dsup, dinf, gcdabc = cgcd(gcdab, c);
      long bc = b * c, cc = c * c, P = bb - a3 * c;
      dsup = minss(dsupb, sfloordiv(bc, a9)); /* Q >= 0 */
      /* bc-9ad <= 4x / 3c^2 */
      dinf = c? maxss(dinfb, sceildiv(bc - ((4 * x) / (cc * 3)), a9)): dinfb;
      for (d = dinf; d <= dsup; d++)
      {
        long Q, R, D, DF;
        GEN F;
        if (cgcd(gcdabc, d) > 1) continue;
        Q = bc - a9 * d; if (Q < 0 || Q > P) continue;
        if (Q == 0 && b <= 0) continue;
        R = cc - b3 * d; if (P > R) continue;
        D = 4 * P * R - Q * Q; DF = D / 3; if (DF > x || DF < xinf) continue;
        if (P == Q && (Q == R || labs(b) >= labs(3 * a - b))) continue;
        if (P == R && (a > labs(d) || (a == labs(d) && labs(b) >= labs(c))))
          continue;
        if ((F = checkU(a, b, c, d, P, Q, R, D))) gel(RET, ct++) = F;
      }
    }
  }
  setlg(RET, ct); return RET;
}

/* x >= xinf >= 1 */
static GEN
cubicreal(long x, long xinf)
{
  double sqx, sqx4, sq13, sq3x;
  long A, bsup, binf, cmin, cplus, cminus, Dmin, Dsup;
  GEN V, S;

  if (x < 148) return NULL;
  sqx = sqrt((double)x); sq3x = sqrt((double)(3 * x)); sqx4 = sqrt(sqx);
  sq13 = sqrt(13.);
  cplus = ((-35 + 13 * sq13) * x) / 216;
  cminus = ceil((-(35 + 13 * sq13) * x) / 216);
  cmin = ceil(-sq3x / 4);
  Dmin = ceil(-4./27 * sqx);
  Dsup = sq3x / 36;
  A = floor(sqx4 * 2. / sqrt(27));
  bsup = floor(sqx4 * 2. / sqrt(3));
  binf = ceil(-sqx4);
  S = mkvecsmalln(10, x, xinf, (long)sqx, cplus, cminus, cmin, Dmin, Dsup,
                  bsup, binf);
  V = nflist_parapply("_nflist_S3R_worker", mkvec(S), identity_ZV(A));
  V = myshallowconcat1(V); return lg(V) == 1? NULL: V;
}

GEN
nflist_S3I_worker(GEN ga, GEN S)
{
  long a = itos(ga), a3 = a * 3, a9 = a * 9, b, c, d, ct = 1;
  long x = S[1], xinf = S[2], cplus = S[3], Dsup = S[4], limb = S[5];
  long x4 = x * 4, csupa = usqrtn(cplus / a, 3), dsupa = Dsup / a;
  GEN RET = cgetg(x, t_VEC);

  for (b = 0; b <= limb; b++)
  {
    long b3 = b * 3, bb = b * b, gcdab = cgcd(a, b);
    long apb = a + b, amb = a - b;
    long dsupb = b? minuu(dsupa, cplus / (bb * b)): dsupa;
    long csup = b? min(csupa, 4 * Dsup / b): csupa;
    for (c = -csup; c <= csup; c++)
    {
      long dsup = dsupb, dinf = b? -dsupb: 1, gcdabc = cgcd(gcdab, c);
      long bc = b * c, cc = c * c, P = bb - a3 * c;
      if (c)
      { /* c^2|bc-9ad| <= 4x */
        long t = x4 / cc;
        dsup = minss(dsup, sfloordiv(bc + t, a));
        dinf = maxss(dinf, sceildiv(bc - t, a));
      }
      dinf = maxss(dinf, sceildiv(-amb * (amb + c) + 1, a));
      dsup = minss(dsup, (apb * (apb + c) - 1) / a);
      for (d = dinf; d <= dsup; d++)
      {
        GEN F;
        long Q, R, D, DF;
        if (!d || cgcd(gcdabc, d) > 1) continue;
        if (d * (d - b) + a * (c - a) <= 0) continue;
        Q = bc - a9 * d;
        R = cc - b3 * d; D = 4 * P * R - Q * Q; DF = D / 3;
        if (DF > -xinf || DF < -x) continue;
        if ((F = checkU(a, b, c, d, P, Q, R, D))) gel(RET, ct++) = F;
      }
    }
  }
  setlg(RET, ct); return RET;
}

static GEN
cubicimag(long x, long xinf)
{
  double sqx, sqx4;
  long lima, limb, Dsup, cplus;
  GEN V, S;

  if (x < 31) return NULL;
  sqx = sqrt((double)x / 27); sqx4 = sqrt(sqx);
  cplus = (11 + 5 * sqrt(5.)) / 8 * x;
  Dsup = 3 * sqx;
  lima = 2 * sqx4;
  limb = sqrt(3.) * 2 * sqx4;
  S = mkvecsmall5(x, xinf, cplus, Dsup, limb);
  V = nflist_parapply("_nflist_S3I_worker", mkvec(S), identity_ZV(lima));
  V = myshallowconcat1(V); return lg(V) == 1? NULL: V;
}

static GEN
makeS3resolvent(GEN T, long flag)
{
  GEN P, d, f = NULL;
  (void)nfcoredisc2(T, &d, odd(flag)? &f: NULL);
  P = quadpoly_i(d); return f? mkvec2(P, f): P;
}

static GEN
makeS3vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN R, I;
  long x, xinf;

  if (field) return makeDLvec(3, X, Xinf, field, s);
  x = itos(X); xinf = itos(Xinf);
  R = (s <= 0)? cubicreal(x, xinf): NULL;
  I = s? cubicimag(x, xinf): NULL;
  switch (s)
  {
    case 0: return R;
    case 1: return I;
    case -1: return R? (I? shallowconcat(R, I): R): I;
    default: if (!R && !I) return NULL; /* -2 */
             return mkvec2(R? R: cgetg(1,t_VEC), I? I: cgetg(1,t_VEC));
  }
}

/**********************************************************************/
/*                                 C4                                 */
/**********************************************************************/

static GEN
makepolC4(GEN S, GEN T)
{
  GEN V = cgetg(7, t_POL);
  V[1] = evalsigne(1)|evalvarn(0);
  gel(V, 6) = gen_1;
  gel(V, 5) = gen_0;
  gel(V, 4) = S;
  gel(V, 3) = gen_0;
  gel(V, 2) = T; return V;
}

static GEN
C4qfbsolve(GEN Q, GEN D)
{
  GEN v = qfbsolve(Q, D, 1), w;
  long i, c, n = lg(v) - 1;

  w = cgetg(2 * n + 1, t_VEC);
  for (i = c = 1; i <= n; i++)
  {
    GEN BC = gel(v, i), B = gel(BC,1), C = gel(BC,2);
    gel(w, c++) = absi_shallow(B);
    if (!absequalii(B, C)) gel(w, c++) = absi_shallow(C);
  }
  setlg(w, c); return gtoset_shallow(w);
}

/* D squarefree in [D,factor(D)] form, D = B^2 + C^2,
 * A*(odd part of D) = n2 = prod_{odd p | n} p, v2 = v2(n) */
static GEN
polsubC4_D(GEN Q, GEN A, GEN Dfa, GEN n2, long v2, long s, long fli)
{
  GEN v, S, mS, AD, A2D, D = gel(Dfa,1), vB = C4qfbsolve(Q, Dfa);
  long i, c, l = lg(vB), A4 = Mod4(A); /* 1 or 3 */

  AD = mpodd(D)? n2: shifti(n2, 1);
  A2D = mulii(A, AD);
  S = mulsi(-2, AD); mS = negi(S);
  v = cgetg(2 * l - 1, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    GEN B = gel(vB, i), T;
    long B4 = Mod4(B);
    int p = (s <= 0), m = !!s;
    if (v2 <= 2 && odd(B4)) continue;
    if (!v2)
    { if (((A4 + B4) & 3) == 1) m = 0; else p = 0; }
    else if (fli)
    {
      if (v2 == 3)
      { if (!odd(B4)) continue; }
      else if (v2 == 2)
      { if (((A4 + B4) & 3) == 1) p = 0; else m = 0; }
    }
    if (!p && !m) continue;
    T = mulii(A2D, subii(D, sqri(B)));
    if (p) gel(v, c++) = makepolC4(S, T);
    if (m) gel(v, c++) = makepolC4(mS, T);
  }
  setlg(v, c); return v;
}


/* vector of distinct primes -> squarefree famat */
static GEN
P2fa(GEN P) { return mkmat2(P, const_col(lg(P)-1, gen_1)); }
/* vector of distinct primes -> [factorback, P2fa] */
static GEN
P2Nfa(GEN P) { return mkvec2(ZV_prod(P), P2fa(P)); }
/* P = prime divisors of f different from ell; nf = Q or quadratic */
static GEN
Pell2prfa(GEN nf, GEN P, long ell, GEN f)
{
  long v = Z_lval(f, ell);
  if (v) P = ZV_sort_shallow(vec_append(P, utoipos(ell)));
  P = nf_pV_to_prV(nf, P); settyp(P, t_COL); P = P2fa(P);
  if (v)
  { /* add pr^{2e} for all pr | ell */
    long i, l = lg(gel(P,1));
    for (i = 1; i < l; i++)
    {
      GEN pr = gcoeff(P,i,1);
      if (equaliu(pr_get_p(pr), ell)) gcoeff(P,i,2) = utoipos(v * pr_get_e(pr));
    }
  }
  return P;
}
static int
ZV_is_1(GEN x, long i0)
{
  long i, l = lg(x);
  for (i = i0; i < l; i++) if (!equali1(gel(x,i))) return 0;
  return 1;
}
static int
zv_is_1(GEN x, long i0)
{
  long i, l = lg(x);
  for (i = i0; i < l; i++) if (x[i] != 1) return 0;
  return 1;
}

/* n > 0, D sqfree, sum2sq(odd(D)? D: 4*D) is true */
static GEN
polsubcycloC4_i(GEN n, long s, long fli, GEN D)
{
  GEN fa = NULL, P, Q, v, n2;
  long v2;

  if (typ(n) == t_VEC) { fa = gel(n,2); n = gel(n,1); }
  if (s == 1 || equali1(n)) return NULL;
  /* s = -1, 0 or 2 */
  v2 = vali(n); if (fli && (v2 == 1 || v2 > 4)) return NULL;
  if (!fa) fa = Z_factor(n);
  P = gel(fa,1);
  if (fli && !ZV_is_1(gel(fa,2), v2? 2: 1)) return NULL;
  n2 = ZV_prod(v2? vecsplice(P, 1): P); /* odd part of rad(n) */
  Q = mkqfb(gen_1, gen_0, gen_1, utoineg(4));
  if (D)
  {
    GEN A, PD, LD;
    if (fli && mpodd(D) == (v2 == 4)) return NULL;
    if (!(A = divide(n2, mpodd(D) ? D : gmul2n(D, -1)))) return NULL;
    (void)Z_smoothen(D, P, &PD, &LD);
    D = mkvec2(D, mkmat2(PD, LD));
    v = polsubC4_D(Q, A, D, n2, v2, s, fli);
  }
  else
  {
    long c, i, lv, l = lg(P);
    GEN M2 = NULL;
    c = (v2 && v2 < 4)? 2:  1; /* leave 2 in P if 16 | n */
    if (c == 2) M2 = mkmat2(mkcol(gen_2),mkcol(gen_1));
    for (i = v2? 2: 1; i < l; i++) /* odd prime divisors of n */
      if (Mod4(gel(P,i)) == 1) gel(P, c++) = gel(P,i);
    setlg(P, c);
    v = divisors_factored(P2Nfa(P)); lv = lg(v);
    for (i = c = 1; i < lv; i++)
    {
      GEN A, D = gel(v,i), d = gel(D,1);
      if (M2) /* replace (odd) D by 2*D */
      {
        gel(D,1) = shifti(d,1);
        gel(D,2) = famat_mul(M2, gel(D,2));
      } else if (i == 1) continue; /* ommit D = 1 */
      A = diviiexact(n2, mpodd(d)? d: shifti(d,-1));
      gel(v,c++) = polsubC4_D(Q, A, D, n2, v2, s, fli);
    }
    if (c == 1) return NULL;
    setlg(v, c); v = shallowconcat1(v);
  }
  return v;
}
static GEN
polsubcycloC4(GEN n, long s)
{
  long i, l, c;
  GEN D = divisors_factored(n);
  l = lg(D);
  for (i = 2, c = 1; i < l; i++)
  {
    GEN v = polsubcycloC4_i(gel(D,i), s, 1, NULL);
    if (v) gel(D,c++) = v;
  }
  setlg(D, c); return myshallowconcat1(D);
}

/* x^2 + a */
static GEN
X2p(GEN a) { return deg2pol_shallow(gen_1, gen_0, a, 0); }
/* x^2 - a */
static GEN
X2m(GEN a) { return deg2pol_shallow(gen_1, gen_0, negi(a), 0); }
/* y^2 - a */
static GEN
Y2m(GEN a) { return deg2pol_shallow(gen_1, gen_0, negi(a), 1); }

static GEN
makeC4(GEN N, GEN field, long s)
{
  GEN D;
  long i, c;

  if (s == 1) return NULL;
  if (field)
  {
    GEN d = checkfield(field, 2);
    if (signe(d) < 0 || !divissquare(N, powiu(d,3))) return NULL;
    D = mkvec(d);
  }
  else D = divisorsdisc(cored(N, 3), 0);
  for (i = c = 1; i < lg(D); i++)
  {
    GEN cond, v, d = gel(D, i);
    if (sum2sq(d) && Z_issquareall(divii(N, powiu(d, 3)), &cond)
        && (v = polsubcycloC4_i(mulii(d,cond),s,1, mpodd(d)? d: shifti(d,-2))))
          gel(D, c++) = v;
  }
  if (c == 1) return NULL;
  setlg(D, c); return sturmseparate(myshallowconcat1(D), s, 4);
}

static GEN
condrel_i(GEN P, GEN pol)
{
  GEN bnf = bnfY(P), T = gcoeff(nffactor(bnf, pol), 1, 1);
  GEN f = gel(rnfconductor0(bnf, T, 2), 1);
  GEN id = gel(f, 1), arch = gel(f, 2), co = gcoeff(id, 1, 1);
  if (ZM_isscalar(id, co)) id = co;
  return mkvec2(P, gequal0(arch) ? id : mkvec2(id, arch));
}
static GEN
condrel(GEN P, GEN pol, long flag)
{ return odd(flag)? condrel_i(P, pol): P; }
static GEN
condrel_dummy(GEN P, long flag)
{ return odd(flag)? mkvec2(P, gen_1): P; }
static GEN
condrelresolvent(GEN pol, long d, long flag)
{ return condrel(mynfsubfield(pol, d), pol, flag); }


static GEN
makeC4resolvent(GEN pol, long flag)
{
  GEN d; (void)nfcoredisc(pol, &d);
  return condrel(quadpoly_i(d), pol, flag);
}

static GEN
C4vec(GEN X, GEN Xinf, GEN m, long s)
{
  GEN v, M, inf, m3 = powiu(m, 3), limf = gfloorsqrtdiv(X, m3);
  long l, n, c;
  pari_sp av;
  inf = cmpiu(Xinf, 500) >= 0? gceilsqrtdiv(Xinf, m3): gen_1;
  l = itos(subii(limf, inf)) + 2;
  M = mpodd(m)? m: shifti(m, -2); av = avma;

  v = const_vec(l-1, cgetg(1,t_VEC));
  for (n = c = 1; n < l; n++)
  {
    GEN w, cond = addui(n-1, inf);
    if ((w = polsubcycloC4_i(mulii(m, cond), s, 1, M))) gel(v, c++) = w;
    if ((n & 0xfff) == 0 && gc_needed(av, 3))
    { /* let parisizemax handle some of it */
      if (DEBUGMEM>1) pari_warn(warnmem,"C4vec, n = %ld/%ld", n, l-1);
      v = gerepilecopy(av, v);
    }
  }
  setlg(v, c); return myshallowconcat1(v);
}

GEN
nflist_C4vec_worker(GEN m, GEN X, GEN Xinf, GEN gs)
{
  pari_sp av = avma;
  return gerepilecopy(av, C4vec(X, Xinf, m, itos(gs)));
}

static GEN
makeC4vec_i(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v;
  long limD = floorsqrtn(X,3), m, c, snew = s == -2 ? -1 : s;
  if (s == 1) return NULL;
  if (field)
  {
    GEN gm = checkfield(field, 2);
    return sum2sq(gm)? C4vec(X, Xinf, gm, snew): NULL;
  }
  v = cgetg(limD >> 1, t_VEC);
  for (m = 5, c = 1; m <= limD; m += odd(m) ? 3 : 1)
    if (usum2sq(m)) gel(v, c++) = utoipos(m);
  setlg(v, c);
  v = nflist_parapply("_nflist_C4vec_worker", mkvec3(X, Xinf, stoi(snew)), v);
  return myshallowconcat1(v);
}
static GEN
makeC4vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v = makeC4vec_i(X, Xinf, field, s);
  return v? sturmseparate(v, s, 4): NULL;
}

/**********************************************************************/
/*                                 V4                                 */
/**********************************************************************/

static GEN
makeV4(GEN N, GEN field, long s)
{
  GEN V, R;
  long lV, i1, i2, c = 1;
  if (s == 1) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 2);
    if (signe(D) < 0) pari_err_TYPE("makeV4 [real quadratic subfield]", field);
    V = mkvec(D);
  }
  else V = divisorsdisc(N, -1);
  lV = lg(V); R = cgetg((lV - 1) * (lV - 2) >> 1, t_VEC);
  for (i1 = 1; i1 < lV; i1++)
  {
    GEN V2, D1 = gel(V, i1);
    if (s == 0 && signe(D1) < 0) continue;
    if (cmpii(sqri(D1), N) > 0) continue;
    V2 = divisorsdisc(diviiexact(N, absi_shallow(D1)), -1);
    for (i2 = 1; i2 < lg(V2); i2++)
    {
      GEN D2 = gel(V2, i2), D3, D12;
      if (s == 0 && signe(D2) < 0) continue;
      if (s > 0 && signe(D1) > 0 && signe(D2) > 0) continue;
      if ((!field && cmpii(D1, D2) >= 0) || equalii(D1, D2)) continue;
      D12 = mulii(D1, D2); D3 = coredisc(D12);
      if (cmpii(D2, D3) < 0 && !equalii(D1, D3)
          && absequalii(mulii(D12, D3), N))
        gel(R, c++) = mkpoln(5, gen_1, gen_0, mulsi(-2, addii(D1, D2)),
                              gen_0, sqri(subii(D1, D2)));
    }
  }
  if (c == 1) return NULL;
  setlg(R, c); return sturmseparate(R, s, 4);
}

static GEN
makeV4resolvent(GEN pol, long flag)
{
  GEN P, V = mynfsubfields(pol, 2);
  long i;
  if (lg(V) != 4) pari_err_BUG("makeV4resolvent");
  if (flag >= 2)
  {
    if (flag == 2) return V;
    return mkvec3(condrel_i(gel(V, 1), pol),
                  condrel_i(gel(V, 2), pol),
                  condrel_i(gel(V, 3), pol));
  }
  for (i = 1; i <= 3; i++) { P = gel(V, i); if (signe(ZX_disc(P)) > 0) break; }
  return condrel(P, pol, flag);
}

static GEN
polV4(long d1, long d2)
{ return mkpoln(5, gen_1, gen_0, mulss(-2, d1+d2), gen_0, sqrs(d1-d2)); }

GEN
nflist_V4_worker(GEN D1, GEN X, GEN Xinf, GEN gs)
{
  pari_sp av = avma, av2;
  GEN V, W;
  long d2a, e1 = signe(D1), d1 = itos(D1), d1a = labs(d1);
  long limg, limg2, s2 = -1, s = itos(gs);
  long limD2 = itos(sqrti(divis(X, d1a)));
  long limQ = floorsqrtdiv(X, sqru(d1a));

  limg2 = limg = usqrt(d1a);
  if (!odd(d1a))
  { /* limg2 = sqrt(d1a * 4), to be used when d2 is also even */
    long r = d1a - limg*limg;
    limg2 *= 2; if (r >= limg) limg2++;
  }

  if (s == 2 && e1 > 0) s2 = 1; /* forbid d2 > 0 */
  else if (!s) s2 = 0; /* forbid d2 < 0 */
  W = vectrunc_init(2 * limD2);
  V = e1 < 0? W: vectrunc_init(2 * limD2); av2 = avma;
  for (d2a = d1a; d2a <= limD2; d2a++, set_avma(av2))
  {
    long g, d2ag, LIMg;
    GEN D3, d1d2a, d3;
    int p, m;
    if (odd(d2a)) LIMg = limg;
    else
    {
      if ((d2a & 3) == 2 || !(d2a & 15)) continue; /* v2(d2) = 1 or >= 4 */
      LIMg = limg2;
    }
    g = ugcd(d2a, d1a); if (g > LIMg) continue;
    d2ag = d2a / g; if (d2ag > limQ) continue;
    uis_fundamental_pm(d2a, s2, &p, &m);
    if (!p && !m) continue;
    d3 = muluu(d1a / g, d2ag); d1d2a = muluu(d1a, d2a);
    if (p)
    { /* D2 = d2a is fundamental */
      setsigne(d3, e1);
      D3 = Mod4(d3) > 1? shifti(d3, 2): d3; /* now D3 = coredisc(D1*D2) */
      if (abscmpiu(D3, d2a) > 0 && ok_int(mulii(d1d2a, D3), X, Xinf))
      { vectrunc_append(V,  polV4(d1, d2a)); av2 = avma; }
    }
    if (m)
    { /* D2 = - d2a is fundamental */
      int fl;
      setsigne(d3, -e1);
      D3 = Mod4(d3) > 1? shifti(d3, 2): d3; /* now D3 = coredisc(D1*D2) */
      fl = abscmpiu(D3, d2a);
      if (fl < 0 || (!fl && e1 > 0)) continue;
      if (ok_int(mulii(d1d2a, D3), X, Xinf))
      { set_avma(av2); vectrunc_append(W, polV4(d1, -d2a)); av2 = avma; }
    }
  }
  return gerepilecopy(av, mkvec2(e1 < 0? cgetg(1, t_VEC): V, W));
}

static GEN
Sextract(GEN v, long ind)
{
  long j, l;
  GEN w = cgetg_copy(v, &l);
  for (j = 1; j < l; j++) gel(w, j) = gmael(v, j, ind);
  return myshallowconcat1(w);
}
static GEN
makeV4vec(GEN X, GEN Xinf, GEN field, long s)
{
  long s2, d, dinf, dsup, l, c;
  GEN v;

  if (s == 1) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 2), DSQ = sqri(D);
    if (signe(D) < 0) pari_err_TYPE("makeV4 [real quadratic subfield]", field);
    if (cmpii(DSQ, X) > 0) return NULL;
    dinf = itos(D); dsup = dinf; l = 2; s2 = 0;
  }
  else
  { dinf = 3; dsup = floorsqrtn(X,3); l = dsup << 1; s2 = s? -1: 0; }
  v = cgetg(l, t_VEC); c = 1;
  for (d = dinf; d <= dsup; d++)
  {
    int p, m;
    uis_fundamental_pm(d, s2, &p, &m);
    if (m) gel(v, c++) = utoineg(d);
    if (p) gel(v, c++) = utoipos(d);
  }
  setlg(v, c);
  v = nflist_parapply("_nflist_V4_worker", mkvec3(X, Xinf, stoi(s)), v);
  switch (s)
  {
    case 0: return Sextract(v,1);
    case 2: return Sextract(v,2);
    case -1: return shallowconcat(Sextract(v,1), Sextract(v,2));
    default: return mkvec3(Sextract(v,1), cgetg(1, t_VEC), Sextract(v,2));
  }
}

/**********************************************************************/
/*                                 D4                                 */
/**********************************************************************/
static GEN
archD40() { return mkvec(cgetg(1, t_VECSMALL)); }

static GEN
archD41() { return mkvec2(mkvecsmall(2), mkvecsmall(1)); }

static GEN
archD42() { return mkvec(mkvecsmall2(1, 2)); }

static GEN
getarchD4(long s)
{
  switch (s)
  {
    case 0: return archD40();
    case 1: return archD41();
    case 2: return archD42();
    default: return shallowconcat1(mkvec3(archD40(), archD41(), archD42()));
  }
  return gen_0;
}

/* x = [N, a;0, m] quadratic ideal in HNF, apply quadratic automorphism */
static GEN
aut2(GEN x, long oddD)
{
  GEN N = gcoeff(x,1,1), t = subii(N, gcoeff(x,1,2)), m = gcoeff(x,2,2);
  if (oddD) t = addii(t, m);
  return mkmat2(gel(x,1), mkcol2(modii(t, N), m));
}
/* I a vector of quadratic ideals of same norm */
static GEN
authI(GEN nf, GEN I, GEN *pstable, GEN D)
{
  long l = lg(I), i, oddD;
  GEN v, w;

  if (l == 1) { *pstable = NULL; return I; }
  if (l == 2) { *pstable = mkvecsmall(1); return I; }
  if (l == 3) { *pstable = mkvecsmall2(0,0); gel(I,2) = NULL; return I; }
  v = w = shallowcopy(I); *pstable = zero_zv(l-1); oddD = mpodd(D);
  if (typ(gcoeff(gel(I,1), 1, 1)) != t_INT) /* vector of factorizations */
  {
    w = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(w,i) = idealfactorback(nf,gel(v,i),NULL,0);
  }

  for (i = 1; i < l; i++)
  {
    GEN a = gel(w, i), b;
    long j;
    if (!a) continue;
    b = aut2(a, oddD);
    if (ZM_equal(b, a)) { (*pstable)[i] = 1; continue; }
    for (j = i + 1; j < l; j++)
      if (ZM_equal(b, gel(w,j))) { gel(v,j) = gel(w,j) = NULL; break;}
    if (j == l) pari_err_BUG("makeD4 [conjugate not found]");
  }
  return v;
}

/* kronecker(D, cond) != -1, Arch a vector of arch in t_VECSMALL form */
static GEN
polD4onecond(GEN bnf, GEN G, GEN D, GEN I, GEN Arch)
{
  GEN stable, v0, v1, v2;
  long j, k, m, l, lA, ok = 0, r1 = signe(D) > 0? 2: 0;

  v0 = v1 = v2 = cgetg(1,t_VEC);
  I = authI(bnf, I, &stable, D); l = lg(I); lA = lg(Arch);
  for (j = 1; j < l; j++)
  {
    GEN id = gel(I, j);
    if (!id) continue;
    for (k = 1; k < lA; k++)
    {
      GEN arch = gel(Arch, k), R = NULL;
      long st = lg(arch)-1, lR;
      if (stable[j])
      {
        if (st == 1 && arch[1] == 1) continue;
        if (st != 1) R = mybnrclassfield_X(bnf, mkvec2(id,arch), 2,NULL,NULL,G);
      }
      if (!R) R = mybnrclassfield(bnf, mkvec2(id, arch), 2);
      lR = lg(R); if (lR == 1) continue;
      ok = 1;
      for (m = 1; m < lR; m++)
      {
        GEN P = rnfequation(bnf, gel(R, m));
        if (st == 0 && r1) v0 = vec_append(v0, P);
        else if (st == 1) v1 = vec_append(v1, P);
        else v2 = vec_append(v2, P);
      }
    }
  }
  return ok? mkvec3(v0, v1, v2): NULL;
}

static GEN
makeD4(GEN N, GEN field, long s)
{
  pari_sp av2;
  GEN vD, v, v0, v1, v2, archempty, listarch = getarchD4(s);
  long l, i;

  if (field)
  {
    GEN D = checkfield(field, 2);
    if ((signe(D) < 0 && (s == 0 || s == 1)) || !dvdii(N, sqri(D)))
      return NULL;
    vD = mkvec(D);
  }
  else
    vD = divisorsdisc(cored(N,2), (s == 0 || s == 1)? 0 : -1);
  archempty = mkvec(cgetg(1, t_VECSMALL));
  l = lg(vD); av2 = avma;
  v0 = const_vec(l-1, cgetg(1,t_VEC));
  v1 = const_vec(l-1, cgetg(1,t_VEC));
  v2 = const_vec(l-1, cgetg(1,t_VEC));
  for (i = 1; i < l; i++)
  {
    GEN bnf, G, I, Arch, RET, D = gel(vD, i);
    pari_sp av3 = avma;
    long cond = itou(divii(N, sqri(D)));

    set_avma(av3);
    if (kroiu(D, cond) == -1) continue;
    bnf = Buchall(Y2m(D), nf_FORCE, DEFAULTPREC);
    I = ideals_by_norm(bnf_get_nf(bnf), utoipos(cond));
    Arch = signe(D) > 0 ? listarch : archempty;
    /* restrict to fields which are not Galois over Q [eliminate V4/C4] */
    G = s != 1? mkvec2(galoisinit(bnf, NULL), gen_0): NULL;
    if (!(RET = polD4onecond(bnf, G, D, I, Arch)))
    { set_avma(av3); continue; }
    gel(v0,i) = gel(RET,1);
    gel(v1,i) = gel(RET,2);
    gel(v2,i) = gel(RET,3);
    if (gc_needed(av2, 2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"makeD4");
      gerepileall(av2, 3, &v0,&v1,&v2);
    }
  }
  if      (s == 0) v = myshallowconcat1(v0);
  else if (s == 1) v = myshallowconcat1(v1);
  else if (s == 2) v = myshallowconcat1(v2);
  else
  {
    v0 = myshallowconcat1(v0);
    v1 = myshallowconcat1(v1);
    v2 = myshallowconcat1(v2); v = mkvec3(v0, v1, v2);
    if (s == -1) v = myshallowconcat1(v);
  }
  return v;
}

GEN
nflist_D4_worker(GEN D, GEN X, GEN Xinf, GEN listarch)
{
  pari_sp av = avma, av2;
  GEN bnf, G, vI, v0, v1, v2, Arch, D2 = sqri(D);
  long c0, c1, c2, cond, l = itos(divii(X, D2)) + 1;
  long lmin = itos(ceildiv(Xinf, D2));

  bnf = Buchall(Y2m(D), nf_FORCE, DEFAULTPREC);
  vI = ideallist(bnf, l-1);
  Arch = signe(D) > 0 ? listarch : mkvec(cgetg(1,t_VECSMALL));
  G = lg(Arch) != 3? mkvec2(galoisinit(bnf, NULL), gen_0): NULL;
  av2 = avma;
  v0 = const_vec(l-1, cgetg(1,t_VEC));
  v1 = const_vec(l-1, cgetg(1,t_VEC));
  v2 = const_vec(l-1, cgetg(1,t_VEC)); c0 = c1 = c2 = 1;
  for (cond = lmin; cond < l; cond++)
  {
    pari_sp av3 = avma;
    GEN R, R1, R2, R3;
    if (kroiu(D, cond) == -1) continue;
    if (!(R = polD4onecond(bnf, G, D, gel(vI, cond), Arch)))
    { set_avma(av3); continue; }
    R1 = gel(R,1); if (lg(R1) > 1) gel(v0, c0++) = R1;
    R2 = gel(R,2); if (lg(R2) > 1) gel(v1, c1++) = R2;
    R3 = gel(R,3); if (lg(R3) > 1) gel(v2, c2++) = R3;
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"makeD4vec, cond = %ld/%ld",cond,l-1);
      gerepileall(av2, 3, &v0,&v1,&v2);
    }
  }
  setlg(v0,c0); v0 = myshallowconcat1(v0);
  setlg(v1,c1); v1 = myshallowconcat1(v1);
  setlg(v2,c2); v2 = myshallowconcat1(v2);
  return gerepilecopy(av, mkvec3(v0, v1, v2));
}

static GEN
makeD4vec(GEN X, GEN Xinf, GEN field, long s)
{
  long s2, limdinf, limdsup, c, da;
  GEN v, D;

  if (field)
  {
    GEN D = checkfield(field, 2);
    if (cmpii(sqri(D), X) > 0) return NULL;
    limdsup = limdinf = labs(itos(D));
    s2 = signe(D) < 0? 1: 0;
  }
  else
  {
    limdinf = 3; limdsup = itou(sqrti(X));
    s2 = (s == 0 || s == 1) ? 0 : -1;
  }
  D = cgetg(2 * limdsup + 1, t_VEC); c = 1;
  for (da = limdinf; da <= limdsup; da++)
  {
    int p, m;
    uis_fundamental_pm(da, s2, &p, &m);
    if (p) gel(D, c++) = utoipos(da);
    if (m) gel(D, c++) = utoineg(da);
  }
  setlg(D, c);
  v = nflist_parapply("_nflist_D4_worker", mkvec3(X, Xinf, getarchD4(s)), D);
  if (s >= 0) v = Sextract(v,s+1);
  else
  {
    v = mkvec3(Sextract(v,1), Sextract(v,2), Sextract(v,3));
    if (s == -1) v = shallowconcat1(v);
  }
  return v;
}

/**********************************************************************/
/*                              A4 and S4                             */
/**********************************************************************/
/* FIXME: export */
static GEN
to_principal_unit(GEN nf, GEN x, GEN pr, GEN sprk)
{
  if (pr_get_f(pr) != 1)
  {
    GEN prk = gel(sprk,3);
    x = nfpowmodideal(nf, x, gmael(sprk,5,1), prk);
  }
  return x;
}
static long
ZV_iseven(GEN zlog)
{
  long i, l = lg(zlog);
  for (i = 1; i < l; i++)
    if (mpodd(gel(zlog,i))) return 0;
  return 1;
}

/* x^2 = t (mod bid) solvable ? */
static int
issolvable(GEN nf, GEN t, GEN sprk)
{
  GEN pr = sprk_get_pr(sprk);
  (void)nfvalrem(nf, t, pr, &t);
  t = to_principal_unit(nf, t, pr, sprk);
  return ZV_iseven(sprk_log_prk1(nf, t, sprk));
}

/* true nf, cubic field */
static GEN
makeGid(GEN nf)
{
  GEN P = idealprimedec(nf, gen_2), Sprk;
  GEN bid4 = cgetg(1, t_VEC), bid6 = cgetg(1, t_VEC);
  long l = lg(P), i, parity = mpodd(nf_get_disc(nf));
  if (l == 3)
  {
    if (parity) /* ensure f(P[1]/2) = 2 */
    { swap(gel(P,1), gel(P,2)); }
    else /* make sure e(P[1]/2) = 2 */
    { if (pr_get_e(gel(P,1)) == 1) swap(gel(P,1), gel(P,2)); }
  }
  Sprk = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(Sprk, i) = log_prk_init(nf, gel(P,i), 2, gen_2);
  if (!parity)
  {
    bid4 = log_prk_init(nf, gel(P,1), 4, gen_2);
    if (l == 2) bid6 = log_prk_init(nf, gel(P,1), 6, gen_2);
  }
  return mkvecn(3, Sprk, bid4, bid6);
}

static long
r2(GEN v)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++) if (mpodd(gel(v,i))) break;
  return i - 1;
}

/* list of virtual units whose norm is a square */
static GEN
makevunits(GEN bnf)
{
  GEN cyc = bnf_get_cyc(bnf), G = bnf_get_gen(bnf), nf = bnf_get_nf(bnf), v;
  long rc = r2(cyc), l, i;

  v = cgetg(rc + 1, t_VEC);
  for (i = 1; i <= rc; i++)
  {
    GEN g = idealpows(nf, gel(G, i), itos(gel(cyc, i)) >> 1);
    g = idealsqr(nf, idealred(nf, g));
    g = bnfisprincipal0(bnf, g, nf_GEN | nf_FORCE);
    gel(v, i) = gel(g, 2);
  }
  v = shallowconcat(v, bnf_get_fu(bnf)); l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN u = gel(v,i);
    if (signe(nfnorm(nf, u)) < 0) gel(v,i) = gneg(u); /*norm is a square now*/
  }
  return eltlist2(nf, v);
}

static GEN
A4clean3(GEN v, long c)
{
  if (c)
  {
    GEN w = gtoset_shallow(vecslice(v, 1, c)); /* #w /= 3 */
    if (c != lg(v)-1) w = shallowconcat(w, vecslice(v, c+1, lg(v)-1));
    v = w;
  }
  return v;
}

/* nf cubic field, possible local factors for ideals of square norm p^2 */
static GEN
cubictypedec(GEN nf, GEN p)
{
  GEN P = idealprimedec(nf, p);
  switch (lg(P))
  {
    case 2: return NULL;
    case 3: if (pr_get_f(gel(P,2)) == 2)
              return mkvec(idealhnf_shallow(nf, gel(P,2)));
            return mkvec(idealmul(nf, gel(P, 1), gel(P, 2)));
    default: return mkvec3(idealmul(nf, gel(P, 1), gel(P, 2)),
                           idealmul(nf, gel(P, 2), gel(P, 3)),
                           idealmul(nf, gel(P, 3), gel(P, 1)));
  }
}

/* x = gen_1 or in HNF */
static int
oddnorm(GEN x) { return typ(x) == t_INT || mpodd(gcoeff(x,1,1)); }

/* return k < 4, s.t 2^(2k) = discriminant of quadratic extension over cubic */
static long
quadcubpow(GEN bnf, GEN Gid, GEN ideal, GEN a)
{
  GEN nf = bnf_get_nf(bnf), Sprk = gel(Gid,1);
  long t;
  if (mpodd(nf_get_disc(nf)))
    switch (lg(Sprk))
    { /* 2 = P3, P2*P1, P1*P1*P1 */
      case 2: return issolvable(nf, a, gel(Sprk,1))? 0: 3; /* ideal is odd */
      case 3:
        t = issolvable(nf, a, gel(Sprk,2))? 2: 3;
        if (oddnorm(ideal) && issolvable(nf, a, gel(Sprk,1))) t -= 2;
        return t;
      default:
        t = 3;
        if (oddnorm(ideal))
        {
          if (issolvable(nf,a,gel(Sprk,1))) t--;
          if (issolvable(nf,a,gel(Sprk,2))) t--;
          if (issolvable(nf,a,gel(Sprk,3))) t--;
        }
        else
        { /* exactly 2 of the 3 primes divide ideal, test solvability by 3rd */
          if (!idealval(nf,ideal,sprk_get_pr(gel(Sprk,1))))
          { if (issolvable(nf,a,gel(Sprk,1))) t--; }
          else if (!idealval(nf,ideal,sprk_get_pr(gel(Sprk,2))))
          { if (issolvable(nf,a,gel(Sprk,2))) t--; }
          else
          { if (issolvable(nf,a,gel(Sprk,3))) t--; }
        }
        return t;
    }
  if (lg(Sprk) == 3)
  { /* 2 = P1^2 P2 */
    if (!oddnorm(ideal)) return 3;
    t = issolvable(nf, a, gel(Sprk,2))? 2: 3;
    if (issolvable(nf, a, gel(Gid,2))) t -= 2; /* solvable mod P1^4 */
    else if (issolvable(nf, a, gel(Sprk,1))) t--; /* solvable mod P1^2 */
    return t;
  }
  else
  { /* 2 = P1^3, ideal must be odd */
    if (issolvable(nf, a, gel(Gid,3)))  return 0; /* solvable mod pr^6 */
    if (issolvable(nf, a, gel(Gid,2)))  return 1; /* solvable mod pr^4 */
    if (issolvable(nf, a, gel(Sprk,1))) return 2; /* solvable mod pr^2 */
    return 3;
  }
}

/* idealfactorback for { I } x W[1] x ... assuming all W[i] have d entries */
static GEN
idlist(GEN nf, GEN I, GEN W)
{
  long i, j, d, l = lg(W);
  GEN v, w;
  if (l == 1) return mkvec(I? I: gen_1);
  w = gel(W,1); d = lg(w)-1;
  if (!I) v = w;
  else
  {
    v = cgetg(d + 1, t_VEC);
    for (j = 1; j <= d; j++) gel(v,j) = idealmul(nf, I, gel(w,j));
  }
  for (i = 2; i < l; i++)
  {
    long nv = lg(v)-1, c, k;
    GEN V = cgetg(d*nv + 1, t_VEC);
    w = gel(W,i);
    for (j = c = 1; j <= nv; j++)
      for (k = 1; k <= d; k++) gel(V, c++) = idealmul(nf, gel(v,j), gel(w,k));
    v = V;
  }
  return v;
}

static GEN
issquareclass(GEN bnf, GEN x, long rc)
{
  GEN v, d, cyc = bnf_get_cyc(bnf), nf = bnf_get_nf(bnf);
  GEN e = isprincipal(bnf, x);
  long l = lg(cyc), j;

  v = cgetg(l, t_VEC);
  for (j = 1; j <= rc; j++)
  {
    if (mpodd(gel(e,j))) return NULL;
    gel(v, j) = subii(gel(cyc,j), gel(e,j));
  }
  for (; j < l; j++)
  {
    GEN t = subii(gel(cyc,j), gel(e,j));
    if (mpodd(t)) t = addii(t, gel(cyc,j));
    gel(v, j) = t;
  }
  /* all exponents are even */
  x = isprincipalfact(bnf, x, bnf_get_gen(bnf), v, nf_GENMAT | nf_FORCE);
  /* reduce generator mod squares */
  x = gel(x,2); x = nffactorback(nf, gel(x,1), ZV_to_Flv(gel(x,2), 2));
  x = nfmul(nf, x, nfsqr(nf, idealredmodpower(nf, x, 2, 0)));
  x = Q_remove_denom(x, &d); return d? gmul(x, d): x;
}
/* v a vector of ideals of same norm */
static GEN
S4makeidclass(GEN bnf, GEN v, long rc)
{
  long j, c, l = lg(v);
  GEN w = cgetg(l, t_VEC);
  for (j = c = 1; j < l; j++)
  {
    GEN N, a, I = gel(v,j);
    if (typ(I) == t_INT) a = gen_1;
    else
    {
      if (!(a = issquareclass(bnf, I, rc))) continue; /* I^2 = (a)(mod K^*)^2 */
      N = nfnorm(bnf,a);
      if (signe(N) < 0) a = gneg(a);/* Norm(a)=|N| is a square */
    }
    gel(w, c++) = mkvec2(I, a);
  }
  setlg(w, c); return w;
}

/* L squarefree outside of 2, v2(L) <= 4, P = prime divisors of L.
 * Write L = N*2^v2, v2 <= 4, N odd sqfree.
 * List of squarefree ideals A of norm N^2 (and 4N^2 if v2 > 1) */
static GEN
S4makeid(GEN bnf, long isA4, long v2, GEN P)
{
  GEN V, V3, v, id, w, d2, nf = bnf_get_nf(bnf);
  long c, c3, i, k, l, n2, rc;

  l = lg(P); V = cgetg(l, t_VEC); V3 = cgetg(l, t_VEC);
  for (i = v2? 2: 1, c = c3 = 1; i < l; i++)
  {
    GEN p = gel(P, i), d = cubictypedec(nf, p);
    if (!d) return NULL;
    if (lg(d) == 4) gel(V3, c3++) = d; else gel(V,c++) = gel(d, 1);
  }
  d2 = (v2 > 1)? cubictypedec(nf, gen_2): NULL;
  if (isA4)
  { /* choose representative in C3 orbit */
    if (c3 > 1) gel(V, c++) = gmael(V3, --c3, 1);
    else if (d2 && lg(d2) == 4) d2 = mkvec(gel(d2,1));
  }
  setlg(V,c); setlg(V3,c3);
  id = c > 1? idealfactorback(nf, V, NULL, 0): NULL;
  v = idlist(nf, id, V3); rc = r2(bnf_get_cyc(bnf));
  if (!d2) return mkvec(S4makeidclass(bnf, v, rc));
  n2 = lg(d2)-1; l = lg(v); w = cgetg(n2 * (l-1) + 1, t_VEC);
  for (i = c = 1; i < l; i++)
    for (k = 1; k <= n2; k++) gel(w, c++) = idealmul(nf, gel(v,i), gel(d2,k));
  return mkvec2(v2 == 4? cgetg(1,t_VEC): S4makeidclass(bnf, v, rc),
                S4makeidclass(bnf, w, rc));
}

static int
checkS4data(GEN x)
{ return lg(x) == 6 && typ(gel(x, 5)) == t_VECSMALL; }

static GEN
S4data(GEN pol, long s)
{
  GEN bnf, nf, lvunit, Gid, sgnu;
  long isA4;

  if (checkS4data(pol)) return pol;
  bnf = Buchall(pol, nf_FORCE, DEFAULTPREC);
  nf = bnf_get_nf(bnf); Gid = makeGid(nf);
  lvunit = makevunits(bnf); isA4 = Z_issquare(nf_get_disc(nf));
  sgnu = (s != -1 && nf_get_r1(nf) == 3)? nfsign(nf, lvunit): gen_0;
  return mkvecn(5, bnf, lvunit, Gid, sgnu, mkvecsmall(isA4));
}
static GEN
S4_get_disc(GEN S) { return nf_get_disc(bnf_get_nf(gel(S,1))); }

static int
cmp2(void *E,GEN x,GEN y)
{ (void)E; return  signe(gel(x,2))==0 ? 1: signe(gel(y,2))==0 ? -1: cmpii(gel(x,2), gel(y,2)); }

/* Find quartic A4 or S4-extensions of Q with resolvent pol and square root of
 * norm of relative discriminant = L; disc(K/Q) = L^2 nfdisc(pol).
 * Here s = -1 or (0, 1, 2) */
static GEN
makeA4S4(GEN pol, GEN L, long s)
{
  GEN DATA = S4data(pol, s), bnf = gel(DATA, 1), nf = bnf_get_nf(bnf);
  GEN lvunit = gel(DATA, 2), Gid = gel(DATA, 3), sgnunit = gel(DATA, 4);
  GEN sgnal0 = NULL, vI, V, P, L2;
  long nu, l, c, c1, i, j, k, v2, r1 = nf_get_r1(nf), isA4 = gel(DATA, 5)[1];

  if (s != -1 && ((r1 == 1 && s != 1) || (r1 == 3 && s == 1))) return NULL;
  if (typ(L) == t_VEC)
  { P = gel(L,1); L = gel(L,2); v2 = vali(L); }
  else
  {
    GEN fa;
    v2 = vali(L); if (v2 > 4) return NULL;
    fa = Z_factor(L); if (!ZV_is_1(gel(fa,2), v2? 2: 1)) return NULL;
    P = gel(fa,1);
  }
  L2 = v2? shifti(L, -v2): L;
  vI = S4makeid(bnf, isA4, v2, P); if (!vI) return NULL;
  l = lg(vI); nu = lg(lvunit) - 1;
  V = cgetg(RgVV_nb(vI) * nu + 1, t_VEC); c = 1; c1 = 0;
  for (k = 1; k < l; k++) /* l = 2 or 3 */
  {
    GEN I = gel(vI, k);
    int norm1 = k == 1 && equali1(L2);
    for (j = 1; j < lg(I); j++)
    {
      GEN ideal = gmael(I, j, 1), al0 = gmael(I, j, 2);
      if (s != -1 && r1 == 3) sgnal0 = nfsign(nf, al0);
      for (i = norm1? 2 : 1; i <= nu; i++)
      {
        GEN T, a1, a2, a3, a;
        if (sgnal0 && !!s == zv_equal0(Flv_add(sgnal0, gel(sgnunit,i), 2)))
          continue;
        a = nfmul(nf, al0, gel(lvunit, i));
        if (v2 != quadcubpow(bnf, Gid, ideal, a) + (k == 1? 0 : 1)) continue;
        if (isA4 && norm1) c1++;
        a = nf_to_scalar_or_alg(nf, a);
        T = QXQ_charpoly(a, nf_get_pol(nf), 0);
        a1 = gel(T,4); a2 = gel(T,3); a3 = negi(gel(T,2));
        T = mkpoln(5, gen_1, gen_0, shifti(a1,1), mulsi(-8,sqrti(a3)),
                   subii(sqri(a1), shifti(a2, 2)));
        gel(V, c++) = isA4? polredabs(T): T;
      }
    }
  }
  if (c == 1) return NULL;
  setlg(V, c); return isA4? A4clean3(V, c1): V;
}

/* A4 fields of square root discriminant = N and "signature" s */
static GEN
makeA4_i(GEN N2, GEN field, long s)
{
  GEN N, v;
  if (s == 1 || !Z_issquareall(N2, &N)) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 3), d, cond;
    if (!Z_issquareall(D, &d) || !(cond = divide(N,d))
        || !(v = makeA4S4(field, cond, s))) return NULL;
  }
  else
  {
    GEN D = divisors(N);
    long i, cv, l = lg(D);
    v = cgetg(l, t_VEC);
    for (i = cv = 1; i < l; i++)
    {
      GEN w, m = gel(D, i), n = gel(D, l-i), C = makeC3_f(m);
      long j, c, lC = lg(C);
      for (j = c = 1; j < lC; j++)
        if ((w = makeA4S4(gel(C,j), n, s))) gel(C,c++) = w;
      if (c == 1) continue;
      setlg(C,c); gel(v,cv++) = shallowconcat1(C);
    }
    setlg(v,cv); v = myshallowconcat1(v);
  }
  return v;
}
static GEN
makeA4(GEN N, GEN field, long s)
{
  GEN v = makeA4_i(N, field, maxss(s, -1));
  return v? sturmseparate(v, s, 4): NULL;
}

static GEN
makeS4_i(GEN N, GEN field, long s)
{
  GEN v;
  if (field)
  {
    GEN q, f, D = checkfield(field, 3);
    if (!(q = divide(N, D))) return NULL;
    setsigne(q, s == 2? -signe(q): 1);
    if (!Z_issquareall(q, &f) || !(v = makeA4S4(field, f, s))) return NULL;
  }
  else
  {
    GEN f, M = divisors(N);
    long i, cv, l = lg(M);
    v = cgetg(l, t_VEC); if (!odd(s)) s = 0;
    for (i = cv = 1; i < l; i++)
      if (Z_issquareall(gel(M, l-i), &f))
      {
        GEN w, D;
        long j, c, lD;
        if (!(D = makeDL(3, gel(M,i), NULL, s))) continue;
        lD = lg(D);
        for (j = c = 1; j < lD; j++)
          if ((w = makeA4S4(gel(D,j), f, s))) gel(D, c++) = w;
        if (c == 1) continue;
        setlg(D, c); gel(v, cv++) = shallowconcat1(D);
      }
    if (cv == 1) return NULL;
    setlg(v,cv); v = shallowconcat1(v);
  }
  return v;
}
static GEN
makeS4(GEN N, GEN field, long s)
{
  GEN v = makeS4_i(N, field, maxss(s, -1));
  return v? sturmseparate(v, s, 4): NULL;
}

static long
gal_get_order(GEN G) { return degpol(gal_get_pol(G)); }

/* P is monic */
static GEN
makeA4S4resolvent(GEN P, long flag)
{
  GEN R, a0 = gel(P,2), a1 = gel(P,3), a2 = gel(P,4), a3 = gel(P,5);
  GEN b0 = subii(mulii(a0, subii(shifti(a2,2), sqri(a3))), sqri(a1));
  GEN b1 = subii(mulii(a3, a1), shifti(a0,2));
  R = mkpoln(4, gen_1, negi(a2), b1, b0); setvarn(R, varn(P));
  R = polredabs(R);
  return flag? mkvec2(R, sqrti(divii(nfdisc(P), nfdisc(R)))): R;
}

static GEN
A4S4_fa(GEN DATA, GEN fa, long cond, long s)
{
  pari_sp av = avma;
  GEN w, P = gel(fa,1), E = gel(fa,2);
  if (odd(cond))
  { if (!zv_is_1(E, 1)) return gc_NULL(av); }
  else
    if (E[1] > 4 || !zv_is_1(E, 2)) return gc_NULL(av);
  if (!(w = makeA4S4(DATA, mkvec2(Flv_to_ZV(P), utoipos(cond)), s)))
    return gc_NULL(av);
  return gerepilecopy(av, w);
}
static GEN
nflist_A4S4_worker_i(GEN P3, GEN X, GEN Xinf, long s)
{
  GEN v, w, F, DATA = S4data(P3, s), D3 = absi_shallow(S4_get_disc(DATA));
  long i, c, f, linf, limf = floorsqrtdiv(X, D3);

  linf = cmpii(Xinf, shifti(D3, 2)) >= 0? ceilsqrtdiv(Xinf, D3): 1;
  v = cgetg(limf - linf + 2, t_VEC);
  F = vecfactoru_i(linf, limf);
  for (f = linf, i = c = 1; f <= limf; f++, i++)
    if ((w = A4S4_fa(DATA, gel(F,i), f, s))) gel(v, c++) = w;
  setlg(v, c); return myshallowconcat1(v);
}
GEN
nflist_A4S4_worker(GEN P3, GEN X, GEN Xinf, GEN gs)
{
  pari_sp av = avma;
  return gerepilecopy(av, nflist_A4S4_worker_i(P3, X, Xinf, gs[1]));
}

static GEN
makeA4S4vec(long A4, GEN X, GEN Xinf, GEN field, long s)
{
  long snew = s == -2? -1: s;
  GEN v;

  if (field)
  {
    GEN D = checkfield(field, 3);
    long sD = signe(D);
    if (A4 != Z_issquare(D) || abscmpii(D, X) > 0 ||
        (sD > 0 && snew == 1) || (sD < 0 && !odd(snew)))
          return NULL;
    v = nflist_A4S4_worker_i(field, X, Xinf, snew);
  }
  else
  {
    v = A4? makeC3vec(X, gen_1, NULL, 0)
          : makeS3vec(X, gen_1, NULL, odd(snew)? snew: 0);
    if (!v) return NULL;
    v = nflist_parapply("_nflist_A4S4_worker",
                             mkvec3(X,Xinf,mkvecsmall(snew)), v);
    v = myshallowconcat1(v);
  }
  return sturmseparate(v, s, 4);
}

/**********************************************************************/
/*                                 C5                                 */
/**********************************************************************/
/* elements in B have the same norm */
static void
C5cleanB(GEN nf, GEN aut, GEN B)
{
  long l = lg(B), c, i, j, k;
  GEN W = const_vecsmall(l - 1, 1);
  for (i = c = 1; i < l; i++)
  {
    GEN bi, d;
    if (!W[i]) continue;
    gel(B, c++) = gel(B,i);
    bi = Q_remove_denom(nfinv(nf, gel(B,i)), &d); /*1/b = bi / d */
    for (j = 1; j <= 3; j++)
    {
      bi = galoisapply(nf, aut, bi);
      for (k = i + 1; k < l; k++)
      {
        GEN a;
        if (!W[k]) continue;
        a = nfmuli(nf, bi, gel(B,k)); /* bi/d * B[k] has norm 1 or -1 */
        if (absequalii(content(a), d)) { W[k] = 0; break; }
      }
    }
  }
  setlg(B, c);
}

static GEN
makepolC5(GEN nf, GEN e, GEN b, GEN aut)
{
  GEN b1 = galoisapply(nf, aut, b), t1 = nfmuli(nf, b, b1);
  GEN b2 = galoisapply(nf, aut, b1);
  GEN t2 = nfmuli(nf, t1, nfmuli(nf, b1, b2));
  GEN v = cgetg(8, t_POL);
  v[1] = evalsigne(1) | evalvarn(0);
  gel(v, 7) = gen_1;
  gel(v, 6) = gen_0;
  gel(v, 5) = mulsi(-10, e);
  gel(v, 4) = mulsi(-5, mulii(e, nftrace(nf, t1)));
  gel(v, 3) = mului(5, mulii(e, subii(e, nftrace(nf,t2))));
  gel(v, 2) = mulii(negi(e), nftrace(nf, nfmuli(nf, t1, t2)));
  if (umodiu(e, 5)) v = ZX_translate(v, gen_m1);
  return ZX_Z_divexact(ZX_z_unscale(v, 5), utoipos(3125));
}

/* b a pr-unit; multiply by a unit so that z u = 1 (mod pr5^2) */
static GEN
C5prim(GEN nf, GEN pr5, GEN z, GEN eps, GEN b)
{
  GEN pol = nf_get_pol(nf);
  long k, j;
  if (typ(b) != t_POL) b = scalarpol_shallow(b, varn(pol));
  for (j = 0; j <= 1; j++)
  {
    GEN g = j ? b : ZXQ_mul(b, eps, pol);
    for (k = 0; k <= 9; k++)
    {
      if (idealval(nf, gsubgs(g, 1), pr5) > 1) return g;
      if (k < 9) g = ZXQ_mul(g, z, pol);
    }
  }
  pari_err_BUG("C5prim");
  return NULL; /* LCOV_EXCL_LINE */
}

static GEN
C5bnf()
{
  GEN bnf = Buchall(polcyclo(5,1), nf_FORCE, DEFAULTPREC), nf = bnf_get_nf(bnf);
  GEN aut = poltobasis(nf, pol_xn(2, 1));
  GEN p5 = idealprimedec_galois(nf, utoipos(5));
  return mkvec3(bnf, aut, p5);
}

static GEN
polsubcycloC5_i(GEN N, GEN T)
{
  GEN bnf, nf, pol, B, aut, z, eps, p5, N5, P;
  long fl5, i, l, v;

  if (!checkcondCL(N, 5, &P)) return NULL;
  if (typ(N) == t_VEC) N = gel(N,1);
  if (!T) T = C5bnf();
  bnf = gel(T, 1); nf = bnf_get_nf(bnf); pol = nf_get_pol(nf);
  aut = gel(T, 2);
  p5 = gel(T, 3); v = varn(pol);
  z = monomial(gen_m1, 1, v); /* tu */
  eps = deg1pol_shallow(gen_1, gen_1, v); /* fu */
  N5 = divis_rem(N, 25, &fl5); if (fl5) N5 = N; /* fl5 is set if 5 \nmid N */
  N5 = mkvec2(N5, P2fa(P));
  B = bnfisintnorm(bnf, N5); l = lg(B);
  for (i = 1; i < l; i++) gel(B, i) = C5prim(nf, p5, z, eps, gel(B, i));
  if (fl5)
  {
    B = matalgtobasis(nf, B);
    C5cleanB(nf, aut, B);
  }
  else
  {
    GEN b5 =  mkpoln(4, gen_m1, gen_1, gen_1, gen_m1); /* norm 25 */
    setvarn(b5, v); B = matalgtobasis(nf, RgXQV_RgXQ_mul(B, b5, pol));
  }
  for (i = 1; i < l; i++) gel(B, i) = makepolC5(nf, N, gel(B, i), aut);
  return B;
}

static GEN
makeC5(GEN N, GEN field, long s)
{
  GEN sqN, v;
  checkfield_i(field, 1);
  if (s > 0 || !Z_ispowerall(N, 4, &sqN)
     || !(v = polsubcycloC5_i(sqN, NULL))) return NULL;
  return s == -2? vecs(3,v): v;
}

GEN
nflist_C5_worker(GEN N, GEN T)
{
  pari_sp av = avma;
  GEN v = polsubcycloC5_i(N, T);
  if (!v) { set_avma(av); return cgetg(1, t_VEC); }
  return gerepilecopy(av, v);
}

static GEN
makeC5vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, F, bnfC5;
  long x, xinf, i, l;

  checkfield_i(field, 1); if (s > 0) return NULL;
  xinf = ceilsqrtn(Xinf, 4);
  x = floorsqrtn(X, 4); bnfC5 = C5bnf();
  if (!odd(xinf)) xinf++;
  if (!odd(x)) x--;
  F = vecfactoroddu_i(xinf, x); l = lg(F);
  for (i = 1; i < l; i++)
    gel(F,i) = mkvec2(utoipos(xinf + ((i - 1) << 1)), zm_to_ZM(gel(F,i)));
  v = nflist_parapply("_nflist_C5_worker", mkvec(bnfC5), F);
  v = myshallowconcat1(v); return s == -2? vecs(3, v): v;
}

/**********************************************************************/
/*                            CL (ell prime)                          */
/**********************************************************************/
/* polredabs iff |nfdisc(pol)| = N */
static GEN
ZX_red_disc(GEN pol, GEN N)
{
  GEN d, B = nfbasis(mkvec2(pol, utoipos(500000)), &d);
  return absequalii(d, N)? polredabs(mkvec2(pol,B)): NULL;
}
/* polredabs iff Xinf <= |nfdisc(pol)| <= X */
static GEN
ZX_red_disc2(GEN pol, GEN Xinf, GEN X)
{
  GEN d, B = nfbasis(mkvec2(pol, utoipos(500000)), &d);
  if (abscmpii(d, X) > 0 || abscmpii(d, Xinf) < 0) return NULL;
  return polredabs(mkvec2(pol,B));
}

/* make CL(f^(ell-1), 0) */
static GEN
makeCL_f(long ell, GEN F)
{
  GEN bnf, P, f = typ(F) == t_VEC? gel(F,1): F;
  if (!checkcondCL(F, ell, &P)) return cgetg(1,t_VEC);
  bnf = bnfY(pol_x(1));
  P = Pell2prfa(bnf_get_nf(bnf), P, ell, f);
  return mybnrclassfield(bnf, P, ell);
}
/* ell odd prime */
static GEN
makeCL(long ell, GEN N, GEN field, long s)
{
  GEN F, v;
  checkfield_i(field, 1);
  if (s > 0 || !Z_ispowerall(N, ell-1, &F)) return NULL;
  v = makeCL_f(ell, F); return s != -2? v: vecs((ell-1)/2, v);
}

static GEN
makeCLresolvent(long ell, GEN pol, long flag)
{
  if (!odd(flag)) return pol_x(0);
  return mkvec2(pol_x(0), sqrtnint(checkfield(pol, ell), ell-1));
}

static GEN
RgXV_polred(GEN x)
{ pari_APPLY_same(polredabs(gel(x,i))); }

GEN
nflist_CL_worker(GEN f, GEN bnf, GEN gell)
{
  pari_sp av = avma;
  return gerepileupto(av, RgXV_polred(mybnrclassfield(bnf, f, gell[1])));
}

static GEN
makeCLvec(long ell, GEN X, GEN Xinf, GEN field, long s)
{
  long em1 = ell - 1, x, xinf, f;
  GEN v, bnf, F;

  checkfield_i(field, 1); if (s > 0) return NULL;
  xinf = ceilsqrtn(Xinf, em1);
  x = floorsqrtn(X, em1); bnf = bnfY(pol_x(1));
  F = cgetg(x - xinf + 2, t_VEC);
  for (f = xinf; f <= x; f++)
    gel(F, f - xinf + 1) = utoipos(f);
  v = nflist_parapply("_nflist_CL_worker", mkvec2(bnf, mkvecsmall(ell)), F);
  v = myshallowconcat1(v); return s == -2? vecs(em1>>1, v): v;
}

/**********************************************************************/
/*                            DL (ell prime)                          */
/**********************************************************************/
/* For metacyclic groups; assume G is Galois and non-abelian */
static GEN
getpol(GEN nf, GEN T)
{
  GEN G = galoisinit(rnfequation(nf, T), NULL);
  return galoisfixedfield(G, vecsplice(gal_get_gen(G), 1), 1, 0);
}

static GEN
makeDL(long ell, GEN N, GEN field, long s)
{
  GEN v, vD, F = N;
  long i, l, c, si = 0, pow = (ell - 1) >> 1;

  if (s > 0 && s != pow) return NULL;
  if (ell != 3 && !Z_ispowerall(N, pow, &F)) return NULL;
  if (field)
  {
    GEN q, D = checkfield(field, 2);
    si = signe(D);
    if ((s > 0 && si > 0) || (!s && si < 0)) return NULL;
    D = absi_shallow(D);
    if (!(q = divide(F, D))) return NULL;
    vD = mkvec2(q, D);
  }
  else vD = divisors(F);
  l = lg(vD); v = cgetg(2 * l, t_VEC);
  for (i = 2, c = 1; i < l; i++) /* omit 1 */
  {
    GEN LD, f, M = gel(vD, i);
    int p, m;
    long j;
    if (!Z_issquareall(gel(vD, l-i), &f)) continue;
    is_fundamental_pm(M, s, &p, &m);
    if (si < 0) p = 0;
    if (si > 0) m = 0;
    if (!(LD = fund_pm(M, p, m))) continue;
    for (j = 1; j < lg(LD); j++)
    {
      GEN D = gel(LD, j), R, bnf, P, G, pol;
      long k, lR;
      if (!checkcondDL(D, f, ell, &P)) continue;
      pol = Y2m(gel(LD,j)); bnf = bnfY(pol);
      G = mkvec2(galoisinit(pol,NULL), gen_2);
      P = Pell2prfa(bnf_get_nf(bnf), P, ell, f);
      R = mybnrclassfield_X(bnf, P, ell, NULL, NULL, G);
      lR = lg(R); if (lR == 1) continue;
      /* L/Q degree ell subfield of R; d(L) = F^pow, F = D f^2 */
      for (k = 1; k < lR; k++) gel(R,k) = polredabs(getpol(bnf, gel(R,k)));
      gel(v, c++) = R;
    }
  }
  if (c == 1) return NULL;
  setlg(v, c); return sturmseparate(myshallowconcat1(v), s, ell);
}
/* ell >= 5 prime */
static GEN
makeDLresolvent(long ell, GEN pol, long flag)
{
  GEN Dpow = checkfield(pol, ell), D, DF, F;
  long d4, pow = (ell - 1) >> 1, si = signe(Dpow);
  D = si > 0 ? sqrtnint(Dpow, pow) : negi(sqrtnint(negi(Dpow), pow));
  d4 = Mod4(D);
  if (d4 == 3 || (d4 == 0 && si > 0 && pol2s(pol))) D = negi(D);
  else if (d4 == 2) D = shifti(D, 2);
  DF = coredisc2(D); D = quadpoly_i(gel(DF, 1)); F = gel(DF, 2);
  return flag? mkvec2(D, F): D;
}

GEN
nflist_DL_worker(GEN P2, GEN X1pow, GEN X0pow, GEN X2, GEN Xinf2, GEN gell)
{
  pari_sp av = avma;
  GEN X, Xinf, G, D, Da, V, bnf = bnfY(P2), nf = bnf_get_nf(bnf);
  long f, c, limf, linf, ell = gell[1];

  G = mkvec2(galoisinit(nf_get_pol(nf),NULL), gen_2);
  D = bnf_get_disc(bnf);
  Da = absi_shallow(D);
  limf = floorsqrtdiv(X1pow, Da);
  linf = cmpii(X0pow, shifti(Da, 2)) >= 0? ceilsqrtdiv(X0pow, Da): 1;
  V = cgetg(limf + 1, t_VEC);
  /* K/Q degree l with D_l Galois closure L/Q and k/Q quadratic resolvent
   * Then d_k = D, d_K = D^(l-1)/2 f^(l-1), d_L = D^l f^(2l-2).
   * Want d_K in [Xinf,X], i.e. d_L in D [Xinf^2,X^2] */
  Xinf = mulii(Da, Xinf2); X = mulii(Da, X2);
  for (f = linf, c = 1; f <= limf; f++)
  {
    pari_sp av2 = avma;
    GEN P, R, F = utoipos(f);
    long lR, k;
    if (!checkcondDL(D, F, ell, &P)) { set_avma(av2); continue; }
    P = Pell2prfa(nf, P, ell, F);
    R = mybnrclassfield_X(bnf, P, ell, X, Xinf, G);
    lR = lg(R); if (lR == 1) { set_avma(av2); continue; }
    for (k = 1; k < lR; k++) gel(R,k) = polredabs(getpol(bnf, gel(R,k)));
    gel(V, c++) = R;
  }
  setlg(V,c); return gerepilecopy(av, myshallowconcat1(V));
}

static GEN
makeDLvec(long ell, GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, X1pow, X0pow, V2;
  long pow = (ell - 1) >> 1;

  checkfield_i(field, 2); if (s > 0 && s != pow) return NULL;
  if (s == pow) s = 1;
  X1pow = sqrtnint(X, pow);
  X0pow = gceilsqrtn(Xinf, pow);
  V2 = field? mkvec(field): makeC2vec(X1pow, gen_1, NULL, s == -2? -1: s);
  if (!V2) return NULL;
  v = nflist_parapply("_nflist_DL_worker", mkvec5(X1pow, X0pow, sqri(X),
                           sqri(Xinf), mkvecsmall(ell)), V2);
  return sturmseparate(myshallowconcat1(v), s, ell);
}
/**********************************************************************/
/*                                 D9                                 */
/**********************************************************************/
/* disc = D^4 g^2 f^6 (quad. subfield: D, cubic subfield: D g^2) */
static GEN
makeD9(GEN N, GEN field, long s)
{
  GEN v, LD, D, D4;
  long i, j, si;

  if ((s > 0 && s != 4) || !Z_issquare(N)) return NULL;
  if (field)
  {
    D = checkfield(field, 2); D4 = powiu(D, 4);
    si = signe(D);
    if ((s > 0 && si > 0) || (s == 0 && si < 0) || !dvdii(N, D4)) return NULL;
    LD = mkvec(field);
  }
  else
  {
    GEN t = divisorsdisc(cored(N, 4), s);
    long l = lg(t);
    LD = cgetg(l, t_VEC);
    for (j = 1; j < l; j++) gel(LD, j) = quadpoly_i(gel(t, j));
  }
  v = cgetg(1, t_VEC);
  for (i = 1; i < lg(LD); i++)
  {
    GEN bnf = bnfY(gel(LD, i)), Q, F, G;
    G = mkvec2(galoisinit(bnf, NULL), gen_2);
    D4 = powiu(bnf_get_disc(bnf), 4); Q = divii(N, D4);
    F = divisors(cored(Q, 6));
    for (j = 1; j < lg(F); j++)
    {
      GEN R = mybnrclassfield_X(bnf, gel(F,j), 9, NULL, NULL, G);
      long k;
      for (k = 1; k < lg(R); k++)
      {
        GEN pol = getpol(bnf, gel(R, k));
        if (pol && (pol = ZX_red_disc(pol, N))) v = shallowconcat(v, pol);
      }
    }
  }
  return sturmseparate(v, s, 9);
}

GEN
nflist_D9_worker(GEN P2, GEN X, GEN Xinf)
{
  pari_sp av = avma;
  GEN v, bnf = bnfY(P2), D2 = bnf_get_disc(bnf);
  GEN G = mkvec2(galoisinit(bnf, NULL), gen_2);
  long l, f, c;

  l = floorsqrtndiv(X, powiu(D2, 4), 6) + 1;
  v = cgetg(l, t_VEC); c = 1;
  for (f = 1; f < l; f++)
  {
    GEN R = mybnrclassfield_X(bnf, utoipos(f), 9, NULL, NULL, G);
    long k, ci, lR = lg(R);
    for (k = ci = 1; k < lR; k++)
    {
      GEN pol = getpol(bnf, gel(R, k));
      if ((pol = ZX_red_disc2(pol, Xinf, X))) gel(R, ci++) = pol;
    }
    if (ci > 1) { setlg(R, ci); gel(v, c++) = R; }
  }
  setlg(v,c); return gerepilecopy(av, myshallowconcat1(v));
}

static GEN
makeD9resolvent(GEN G, long flag)
{
  GEN R = polredabs(galoisfixedfield(G, vecsplice(gal_get_gen(G), 2), 1, 0));
  return condrel(R, gal_get_pol(G), flag);
}

static GEN
makeD9vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN X1pow, V2, v;

  checkfield_i(field,2); if (s > 0 && s != 4) return NULL;
  if (s == 4) s = 1;
  X1pow = sqrtnint(X, 4);
  V2 = field? mkvec(field): makeC2vec(X1pow, gen_1, NULL, s == -2? -1: s);
  if (!V2) return NULL;
  v = nflist_parapply("_nflist_D9_worker", mkvec2(X, Xinf), V2);
  return sturmseparate(myshallowconcat1(v), s, 9);
}
/**********************************************************************/
/* Metacyclic C_a \rtimes C_ell groups with ell prime and a | ell - 1 */
/*                includes F5 = M20, M21, and M42                     */
/**********************************************************************/
/* C_a resolvent field. */
static GEN nfmakenum(long n, long t, GEN N, GEN field, long s);
static GEN nfmakevecnum(long n, long t, GEN X, GEN Xinf, GEN field, long s);

static GEN
MgenF(long ell, GEN d, GEN Fn, long *vell)
{
  GEN F;
  if (umodiu(d, ell)) *vell = 0;
  else
  {
    *vell = Z_lval(Fn, ell) % (ell - 1);
    if (*vell) Fn = diviiexact(Fn, powuu(ell, *vell));
  }
  return Z_ispowerall(Fn, ell - 1, &F)? F: NULL;
}

static int
okgal(GEN P, GEN g)
{
  GEN G = polgalois(P, DEFAULTPREC);
  return equaliu(gel(G,1), g[1]) && equalis(gel(G,2), g[2])
                                 && equaliu(gel(G,3), g[3]);
}
static int
okgal1(GEN P, long d)
{ GEN G = polgalois(P, DEFAULTPREC); return equaliu(gel(G,1), d); }
static int
okgal2(GEN P, long d, long p)
{
  GEN G = polgalois(P, DEFAULTPREC);
  return equaliu(gel(G,1), d) && equalis(gel(G,2), p);
}
static int
ok_s(GEN P, long s) { return s < 0 || pol2s(P) == s; }

/* a | ell - 1, (Z/aZ)^* cyclic, F^(ell-1)*D^((ell-1)/a) */
static GEN
makeMgen(long ell, long a, GEN N, GEN field, long s)
{
  GEN v, Fn, F;
  long i, lv, c, vell, deg = ell * a, drel = (ell - 1) / a;

  if (field)
  {
    GEN d = absi_shallow(checkfield(field, a));
    Fn = gdiv(N, powiu(d, drel));
    if (typ(Fn) != t_INT || !(F = MgenF(ell, d, Fn, &vell))) return NULL;
    v = mkvec(mkvec3(mkvec(field), F, utoi(vell)));
  }
  else
  {
    long s2 = maxss(s, -1);
    v = divisors(cored(N, drel)); lv = lg(v);
    for (i = c = 1; i < lv; i++)
    {
      GEN R, d = gel(v, i); Fn = diviiexact(N, powiu(d, drel));
      if ((F = MgenF(ell, d, Fn, &vell))
           && (R = nfmakenum(a, 1, d, NULL, s2))) /* C_a, disc d */
        gel(v, c++) = mkvec3(R, F, utoi(vell));
    }
    setlg(v, c);
  }
  lv = lg(v);
  for (i = 1; i < lv; i++)
  {
    GEN T = gel(v, i), R = gel(T, 1), F0 = gel(T, 2);
    long vell = itou(gel(T, 3)), lR = lg(R), j;
    for (j = c = 1; j < lR; j++)
    {
      GEN nf = nfY(gel(R,j)), F = F0, K, G;
      long k, ck, l;
      if (vell)
      { /* ell ramified in nf */
        long eell, q;
        GEN pell = getpell(nf, ell, &eell);
        q = (ell - 1) / eell; if (vell % q) continue;
        F = idealmul(nf, F, idealpows(nf, pell, vell / q));
      }
      G = mkvec2(galoisinit(nf, NULL), gen_2);
      K = mybnrclassfield_X(Buchall(nf, nf_FORCE, DEFAULTPREC),
                            F, ell, NULL, NULL, G);
      l = lg(K);
      for (k = ck = 1; k < l; k++)
      {
        GEN q = getpol(nf, gel(K, k));
        if ((deg == 21 || okgal1(q, deg)) && /* automatic for M21;FIXME */
            (q = ZX_red_disc(q, N))) gel(K, ck++) = q;
      }
      if (ck > 1) { setlg(K, ck); gel(R,c++) = K; }
    }
    setlg(R, c); gel(v, i) = myshallowconcat1(R);
  }
  return sturmseparate(gtoset_shallow(myshallowconcat1(v)), s, ell);
}

/* (ell,a) = (5,4), (7,3) or (7,6) */
static GEN
makeMgenresolvent(long ell, long a, GEN pol, long flag)
{
  GEN Dpow = checkfield(pol, ell), G, R, DR, F2, nf, pell, F;

  G = galoissplittinginit(pol, utoipos(a*ell));
  if (gal_get_order(G) != a * ell) pari_err_BUG("nfresolvent [Galois group]");
  R = polredabs(galoisfixedfield(G, vecsplice(gal_get_gen(G), 2), 1, 0));
  if (!flag) return R;
  DR = nfdisc(R);
  if (ell == 5 && a == 4)
  {
    F2 = sqrti(divii(Dpow, DR));
    if (!Z_issquareall(F2, &F))
    {
      long e;
      F2 = divis(F2, 5);
      if (!Z_issquareall(F2, &F)) pari_err_BUG("nfresolvent [F5]");
      nf = nfinit(R, DEFAULTPREC); pell = getpell(nf, 5, &e);
      if (e == 4) pell = idealsqr(nf, pell);
      F = idealmul(nf, F, pell);
    }
  }
  else
  { /* ell == 7 && (a == 3 || a == 6) */
    long v;
    if (a == 3) DR = sqri(DR);
    if (!Z_issquareall(divii(Dpow, DR), &F2))
      pari_err_BUG("nfresolvent [M21/M42]");
    /* F2 = F^3 or 7F^3 or 7^2F^3 */
    v = Z_lval(F2, 7) % 3;
    if (v) F2 = divii(F2, powuu(7, v));
    if (!Z_ispowerall(F2, 3, &F)) pari_err_BUG("nfresolvent [M21/M42]");
    if (v)
    {
      long e;
      nf = nfinit(R, DEFAULTPREC); pell = getpell(nf, 7, &e);
      if (e == 6) v *= 2;
      F = idealmul(nf, F, idealpows(nf, pell, v));
    }
  }
  return mkvec2(R, F);
}

GEN
nflist_Mgen_worker(GEN field, GEN X, GEN Xinf, GEN T)
{
  pari_sp av = avma;
  GEN v, Fn, pell, lpow, bnf = bnfY(field), D = bnf_get_disc(bnf);
  GEN G = mkvec2(galoisinit(bnf, NULL), gen_2);
  long ell = T[1], drel = T[2], deg = T[3], c, e;
  long vd = Z_lval(D, ell), limf, f;

  Fn = divii(X, drel == 1 ? absi_shallow(D) : sqri(D));
  limf = floorsqrtn(Fn, ell - 1);
  pell = getpell(bnf, ell, &e); /* e | a */
  lpow = powuu(ell, (ell - 1) / e);
  v = cgetg(limf + 1, t_VEC);
  for (f = c = 1; f <= limf; f++)
  {
    GEN F = utoipos(f), K;
    long k, ci, lK;

    if (vd)
    {
      GEN fn = powuu(f, ell - 1);
      long imax = minss(e - 1, logint(divii(Fn, fn), lpow));
      F = mkcol2(F, gmulgu(idealpows(bnf, pell, imax), f));
    }
    K = mybnrclassfield_X(bnf, F, ell, NULL, NULL, G); lK = lg(K);
    for (k = ci = 1; k < lK; k++)
    {
      GEN q = getpol(bnf, gel(K, k));
      if (degpol(q) == ell && (deg == 21 || okgal1(q, deg)) && /* FIXME */
          (q = ZX_red_disc2(q, Xinf, X))) gel(K, ci++) = q;
    }
    if (ci > 1) { setlg(K, ci); gel(v, c++) = K; }
  }
  setlg(v, c); v = gtoset_shallow(myshallowconcat1(v));
  return gerepilecopy(av, v);
}

/* (a,ell) = (3,7), (4,5) or (6,7) */
static GEN
makeMgenvec(long ell, long a, GEN X, GEN Xinf, GEN field, long s)
{
  GEN L, v, T;
  long drel = (ell - 1) / a;

  if (field)
  {
    if (degpol(field) != a || !okgal2(field, a, a==3? 1: -1))
      pari_err_TYPE("makeMgenvec [field]", field);
    L = mkvec(field);
  }
  else L = nfmakevecnum(a, 1, drel == 1? X: sqrti(X), gen_1, NULL, maxss(s,-1));
  if (!L) return NULL;
  T = mkvecsmall3(ell, drel, ell * a);
  v = nflist_parapply("_nflist_Mgen_worker", mkvec3(X, Xinf, T), L);
  return sturmseparate(myshallowconcat1(v), s, ell);
}

/**********************************************************************/
/*                        A5 by table lookup                          */
/**********************************************************************/
/* V a vector of [T, n] sorted wrt t_INT n. Return elts with Xinf <= n <= X.
 * If flag = 0 return only the T's. */
static GEN
vecslicebyX(GEN V, GEN Xinf, GEN X, long flag)
{
  long l = lg(V), i = 1, c;
  GEN W;
  if (cmpii(Xinf,gmael(V,1,2)) > 0) /* frequent special case */
  {
    i = gen_search(V, mkvec2(NULL,Xinf), NULL, &cmp2);
    if (i > 0) /* found in list, rewind to first occurence */
    { while (i > 1 && equalii(gmael(V, i-1, 2), Xinf)) i--; }
    else /* not in list */
      i = -i;
  }
  W = cgetg(l, t_VEC);
  for (c = 1; i < l; i++)
  {
    GEN C = gmael(V, i, 2), x;
    if (isintzero(C)) /* marker for incomplete slice */
    {
      GEN B = gmael(V, i-1, 2);
      if (equalii(B, X)) break;
      pari_err_DOMAIN("nflist(A5)", "sqrt(N)", ">", B, X);
    }
    if (cmpii(C, X) > 0) break;
    x = RgV_to_RgX(gmael(V, i, 1), 0);
    gel(W, c++) = flag ? mkvec2(x, gmael(V, i, 2)): x;
  }
  setlg(W, c); return W;
}

/* assume 1 <= t < 1000, 1 <= 2s <= n < 100 */
static GEN
nflistfile(const char *suf, long n, long t, long s, long u)
{
  pariFILE *F;
  GEN z;
  char *f = stack_sprintf("%s/nflistdata/%ld/%ld/%ld%s/%ld",
                          pari_datadir, n, t, s, suf, u);
  F = pari_fopengz(f);
  if (!F) pari_err_FILE("nflistdata file",f);
  z = gp_readvec_stream(F->file); pari_fclose(F); return z;
}

static GEN
A5file(const char *suf, long s, long u) { return nflistfile(suf, 5, 4, s, u); }

/* If flag = 0 return only the T's. */
static GEN
vecsliceA5all(const char *suf, long s, ulong sl, GEN Xinf, GEN X, long flag)
{
  long i, l;
  GEN V;
  ulong  uinf = itou(divis(Xinf, sl));
  ulong  usup = itou(divis(X, sl));
  l = usup-uinf+2;
  V = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(V, i) = vecslicebyX(A5file(suf, s, uinf+i-1), Xinf, X, flag);
  return shallowconcat1(V);
}

static GEN
vecsliceA5(long s, GEN Xinf, GEN X, long flag)
{
  return vecsliceA5all("", s, 100000, Xinf, X, flag);
}

static GEN
vecsliceA5cond(long s, GEN Xinf, GEN X, long flag)
{
  return vecsliceA5all("cond", s, 100000, Xinf, X, flag);
}

static GEN
A5vec(GEN X, GEN Xinf, long s, long fl)
{
  GEN L1, L5;
  const char *suf = fl? "cond": "";

  L1 = L5 = NULL;
  if (s <= 0) L5 = vecsliceA5all(suf, 0, 100000, Xinf, X, fl);
  if (s)      L1 = vecsliceA5all(suf, 2, 100000, Xinf, X, fl);
  switch (s)
  {
    case 2: return L1;
    case 0: return L5;
    case -1:
      return shallowconcat(L1, L5);
    default:
      return mkvec3(L5, cgetg(1, t_VEC), L1);
  }
}
static GEN
makeA5_i(GEN N, long s, long fl)
{ return s == 1 ? NULL: A5vec(N, N, s, fl); }
static GEN
makeA5(GEN N, long s)
{
  GEN rN;
  if (!Z_issquareall(N, &rN)) return NULL;
  return makeA5_i(rN, s, 0);
}
static GEN
makeA5cond(GEN N, long s) { return makeA5_i(N, s, 1); }

/* D a sorted t_VECSMALL of conductors; return all [T, d] with d = D[i]
 * for some i and Gal(T) = A5 with s complex places */
GEN
veccond_to_A5(GEN D, long s)
{
  pari_sp av = avma;
  long l, j, lD = lg(D), c = 1;
  GEN W, V = vecsliceA5cond(s, utoi(D[1]), utoi(D[lD-1]), 1);
  l = lg(V);
  W = cgetg(lD, t_VEC);
  for (j = 1; j < lD; j++)
  {
    GEN Xinf = utoi(D[j]);
    long i = gen_search(V, mkvec2(NULL, Xinf), NULL, &cmp2);
    if (i > 0) /* found in list, rewind to first occurence */
    {
      long ii;
      while (i > 1 && equalii(gmael(V, i-1, 2), Xinf)) i--;
      for (ii = i; ii < l && equaliu(gmael(V,ii,2),D[j]); ii++);
      gel(W, c++) = vecslice(V, i, ii-1);
    }
  }
  setlg(W, c); return gerepilecopy(av, shallowconcat1(W));
}

/* Sextic resolvent of A5 field */
static GEN
makeA5resolvent(GEN pol, long flag)
{
  GEN R = cgetg(9, t_POL), D = ZX_disc(pol), c, d, e, f, v;
  GEN c2, d2, e2, c4, df;
  pol = RgX_translate(pol, gdivgs(gel(pol, 6), -5));
  c = gdivgu(gel(pol, 5), 10);
  d = gdivgu(gel(pol, 4), 10);
  e = gdivgu(gel(pol, 3), 5);
  f = gel(pol, 2);
  c2 = gsqr(c); c4 = gsqr(c2); d2 = gsqr(d); e2 = gsqr(e);
  df = gmul(d, f);
  R[1] = evalsigne(1)|evalvarn(0);
  gel(R, 8) = gen_1;
  gel(R, 7) = gen_0;
  gel(R, 6) = gmulsg(-25, gadd(e, gmulsg(3, c2)));
  gel(R, 5) = gen_0;

  v = cgetg(6, t_VEC);
  gel(v, 1) = gmulsg(15, c4);
  gel(v, 2) = gmulsg(8, gmul(c, d2));
  gel(v, 3) = gmulsg(-2, gmul(c2, e));
  gel(v, 4) = gmulsg(3, e2);
  gel(v, 5) = gmulsg(-2, df);
  gel(R, 4) = gmulsg(125, vecsum(v));
  gel(R, 3) = sqrti(D);

  v = cgetg(11, t_VEC);
  gel(v, 1) = gmulsg(-25, gmul(c2, c4));
  gel(v, 2) = gmulsg(-40, gmul(gmul(c, c2), d2));
  gel(v, 3) = gmulsg(-16, gsqr(d2));
  gel(v, 4) = gmulsg(35, gmul(c4, e));
  gel(v, 5) = gmulsg(28, gmul(c, gmul(d2, e)));
  gel(v, 6) = gmulsg(-11, gsqr(gmul(c, e)));
  gel(v, 7) = gmul(e, e2);
  gel(v, 8) = gmulsg(-2, gmul(c2, df));
  gel(v, 9) = gmulsg(-2, gmul(e, df));
  gel(v, 10) = gmul(c, gsqr(f));
  gel(R, 2) = gmulsg(625, vecsum(v));
  R = polredabs(R);
  return odd(flag)? mkvec2(R, gen_1): R;
}

/* For now field ignored. */
static GEN
makeA5vec_i(GEN X, GEN Xinf, GEN field, long s, long fl)
{
  (void)field; if (s == 1) return NULL;
  return A5vec(X, Xinf, s, fl);
}

static GEN
makeA5vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN rX = sqrti(X), sXinf, rXinf = sqrtremi(Xinf, &sXinf);
  if (signe(sXinf)) rXinf = addiu(rXinf, 1);
  return makeA5vec_i(rX, rXinf, field, s, 0);
}

static GEN
makeA5condvec(GEN X, GEN Xinf, GEN field, long s)
{ return makeA5vec_i(X, Xinf, field, s, 1); }

static GEN
makeA56vec_i(GEN V, GEN X, GEN Xinf)
{
  long l = lg(V), i, c;
  GEN W = cgetg(l, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    GEN pol = makeA5resolvent(gel(V, i), 0), D = nfdisc(pol);
    if (cmpii(D, X) <= 0 && cmpii(D, Xinf) >= 0) gel(W, c++) = pol;
  }
  setlg(W, c); return W;
}

static GEN
makeA56vec(GEN X, GEN Xinf, long s)
{
  GEN v;
  if (s == 1 || s == 3 || !(v = makeA5vec(X, Xinf, NULL, s))) return NULL;
  if (s != -2) return makeA56vec_i(v, X, Xinf);
  return mkvec3(makeA56vec_i(gel(v, 1), X, Xinf), cgetg(1, t_VEC),
                makeA56vec_i(gel(v, 3), X, Xinf));
}
static GEN
makeA56(GEN N, long s) { return makeA56vec(N, N, s); }

/* Stupid for now */
static GEN
makeA56resolvent(GEN pol, long flag)
{
  GEN D6 = sqrti(nfdisc(pol)), LD = divisors(D6);
  long i, s;
  pol = polredabs(pol);
  s = pol2s(pol)? 2: 0;
  for (i = 1; i < lg(LD); i++)
  {
    GEN D5 = gel(LD,i);
    if (dvdii(sqri(D5), D6))
    {
      GEN L = vecsliceA5(s, D5, D5, 0);
      long j;
      for (j = 1; j < lg(L); j++)
      {
        GEN P = gel(L, j);
        if (ZX_equal(makeA5resolvent(P, 0), pol))
          return odd(flag)? mkvec2(P, gen_1): P;
      }
    }
  }
  pari_err_BUG("nfresolvent [A56 resolvent not found]");
  return NULL; /* LCOV_EXCL_LINE */
}

/**********************************************************************/
/*                                 C6                                 */
/**********************************************************************/

static GEN
makepol6(GEN P3, GEN P2) { return polcompositum0(P3, P2, 2); }
static GEN
makepol6abs(GEN P3, GEN P2) { return polredabs(makepol6(P3, P2)); }

static GEN
makeC6(GEN N, GEN field, long s)
{
  GEN R, D, d3 = NULL;
  long i, j, lD, s2, c;

  if (s == 1 || s == 2) return NULL;
  if (!field) D = divisorsdisc(cored(N, 3), s);
  else
  {
    if (degpol(field) == 2)
    {
      GEN D2 = nfdisc(field);
      long si = signe(D2);
      if ((s == 3 && si > 0) || (s == 0 && si < 0)
          || !divissquare(N, powiu(D2,3))) return NULL;
      D = mkvec(D2);
    }
    else
    {
      GEN q, D3 = checkfield(field, 3);
      if (!Z_issquareall(D3, &d3)) pari_err_TYPE("makeC6 [field]", field);
      if (!(q = divide(N, sqri(D3)))) return NULL;
      D = divisorsdisc(cored(gcdii(N, powiu(q,3)), 3), s);
    }
  }
  s2 = maxss(s, -1); if (s2 == 3) s2 = 1;
  lD = lg(D); R = cgetg(lD, t_VEC);
  for (i = c = 1; i < lD; i++)
  {
    GEN R0, D2 = gel(D, i), D2a = absi_shallow(D2);
    GEN M = diviiexact(N, powiu(D2a, 3)), F, L, V2;
    long l, l2;
    if (!Z_issquareall(M, &F)) continue;
    if (d3) { L = mkvec(mkvec(field)); l = 2; }
    else
    {
      long k;
      L = divisors(cored(mulii(F, D2a), 2)); l = lg(L);
      for (j = k = 1; j < l; j ++)
      {
        GEN C = makeC3_f(gel(L, j));
        if (lg(C) > 1) gel(L, k++) = C;
      }
      setlg(L, k); l = k; if (l == 1) continue;
    }
    V2 = makeC2(D2a, NULL, s2); l2 = lg(V2);
    R0 = cgetg(l, t_VEC);
    for (j = 1; j < l; j++)
    {
      GEN R3, C3 = gel(L, j);
      long i2, c3, i3, l3 = lg(C3);

      R3 = cgetg(l2 * l3, t_VEC);
      for (i3 = c3 = 1; i3 < l3; i3++)
      {
        GEN P3 = gel(C3, i3);
        for (i2 = 1; i2 < l2; i2++)
        {
          GEN P6 = makepol6(P3, gel(V2, i2));
          if (absequalii(nfdisc(P6), N)) gel(R3, c3++) = P6;
        }
      }
      setlg(R3, c3); gel(R0, j) = R3;
    }
    gel(R, c++) = shallowconcat1(R0);
  }
  setlg(R,c); return sturmseparate(myshallowconcat1(R), s, 6);
}

static GEN
makeC6resolvent(GEN pol, long flag)
{
  GEN V, R3, R = mynfsubfield(pol, 2);
  R3 = (flag >= 2)? mynfsubfield(pol, 3): NULL;
  switch (flag)
  {
    case 0: V = R; break;
    case 1: V = condrel_i(R, pol); break;
    case 2: V = mkvec2(R, R3); break;
    default:V = mkvec2(condrel_i(R, pol), condrel_i(R3, pol)); break;
  }
  return V;
}

/* assume the odd part of M is squarefree, disc is OK */
static void
C6fill(long M, GEN P3, long s, GEN vp,GEN vm)
{
  int p, m;
  uis_fundamental_pm_i(M, s, &p, &m, 1);
  if (p) vectrunc_append(vp, makepol6(P3, X2p(utoineg(M))));
  if (m) vectrunc_append(vm, makepol6(P3, X2p(utoipos(M))));
}

GEN
nflist_C6_worker(GEN P3, GEN X, GEN Xinf, GEN M, GEN T)
{
  pari_sp av = avma;
  GEN D3, f, D32, vp, vm, G, Ginf;
  long i, limD2, l = lg(M), s = T[1];

  if (typ(P3)==t_VEC) { f = gel(P3,2); P3 = gel(P3,1);  } else f = C3pol_f(P3);
  D3 = sqri(f); D32 = sqri(D3); G = divii(X, D32); Ginf = ceildiv(Xinf, D32);
  limD2 = cmpiu(G, T[2]) < 0 ? itou(G) : T[2];

  /* D3 = f^2 is odd, gcd(M,D3) = gcd(M,f); disc = D3^2 / (D3,M)^2 * M^3 */
  vp = vectrunc_init(limD2);
  vm = vectrunc_init(limD2);
  for (i = 1; i < l; i++)
  {
    long m = M[i];
    GEN g;
    if (!odd(m)) continue;
    if (m > limD2) break;
    g = muliu(sqru(m / ugcdiu(f, m)), m);
    if (m != 1 && ok_int(g, G, Ginf)) C6fill(m, P3, s, vp, vm);
    if ((m << 2) <= limD2 && ok_int(shifti(g,6), G, Ginf))
      C6fill(m << 2, P3, s, vp, vm);
    if ((m << 3) <= limD2 && ok_int(shifti(g,9), G, Ginf))
      C6fill(m << 3, P3, s, vp, vm);
  }
  return gerepilecopy(av, mkvec2(vp, vm));
}

static GEN
makeC6vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN T, v, M;

  if (s == 1 || s == 2) return NULL;
  if (field)
  {
    GEN D, f;
    if (degpol(field) == 2)
    {
      long si, m, i, c, l;
      GEN F;
      D = nfdisc(field); si = signe(D);
      if (cmpii(powiu(D, 3), X) > 0 || (s == 3 && si > 0)
          || (s == 0 && si < 0)) return NULL;
      m = itou(D); v = C3vec_F(floorsqrtdiv(X,D), 1, &F); l = lg(v);
      for (i = c = 1; i < l; i++)
      {
        long f = F[i]; /* conductor */
        GEN g = muliu(sqru(m / ugcd(f, m)), m);
        if (ok_int(mulii(powuu(f, 4), g), X, Xinf))
          gel(v, c++) = makepol6(gtopoly(gel(v,i), 0), field);
      }
      setlg(v, c);
      if (s == -2) v = si > 0? vecs14(v, cgetg(1,t_VEC)): vecs(4, v);
      return v;
    }
    D = checkfield(field, 3);
    if (!Z_issquareall(D, &f)) pari_err_TYPE("makeC6 [field]", field);
    if (cmpii(sqri(D), X) > 0) return NULL;
    v = mkvec(mkvec2(field, f));
  }
  else if (!(v = makeC3vec(sqrti(divis(X, 3)), gen_1, NULL, 0))) return NULL;
  T = mkvecsmall2(s, floorsqrtn(X, 3));
  M = vecsquarefreeu(1, T[2]);
  v = nflist_parapply("_nflist_C6_worker", mkvec4(X, Xinf, M, T), v);
  switch (s)
  {
    case -1: return shallowconcat(Sextract(v,1), Sextract(v,2));
    case -2: return vecs14(Sextract(v,1), Sextract(v,2)); /* -2 */
    default: return Sextract(v, s? 2: 1);
  }
}

/**********************************************************************/
/*                             S36 = D66                              */
/**********************************************************************/
static GEN
makeS36(GEN N, GEN field, long s)
{
  GEN vD, P, vp, vm;
  long i, l, cp, cm;
  if (s == 1 || s == 2) return NULL;
  if (s == 3) s = 1;
  if (field)
  {
    long sf = s != -1? pol2s(field): 0/*dummy*/;
    if (s >= 0 && s != sf) return NULL;
    if (degpol(field) == 3)
    {
      GEN d, D = nfcoredisc(field, &d);
      if (!absequalii(mulii(sqri(D), d), N)) return NULL;
      P = mkvec(makepol6abs(field, X2m(d)));
      if (s == -2) { P = vecs(4, P); if (sf) swap(gel(P,1), gel(P,4)); }
      return P;
    }
    else
    {
      GEN D2 = checkfield(field, 2);
      if (!divispowerall(N,  powiu(absi_shallow(D2),3), 4, NULL)) return NULL;
      vD = mkvec(D2);
    }
  }
  else vD = divisorsdisc(cored(N, 3), s);
  l = lg(vD);
  vp = cgetg(l, t_VEC);
  vm = cgetg(l, t_VEC);
  for (i = cp = cm = 1; i < l; i++)
  {
    GEN F, w, P2, D = gel(vD, i), Da = absi_shallow(D);
    long lw, j, s2 = signe(D) > 0? 0: 1;
    if (!Z_ispowerall(divii(N, powiu(Da, 3)), 4, &F)) continue;
    P2 = X2m(D); if (!(w = makeDL(3, mulii(Da, sqri(F)), P2, s2))) continue;
    lw = lg(w);
    for (j = 1; j < lw; j++) gel(w, j) = makepol6abs(gel(w, j), P2);
    if (signe(D) < 0) gel(vm, cm++) = w; else gel(vp, cp++) = w;
  }
  setlg(vp, cp); vp = myshallowconcat1(vp);
  setlg(vm, cm); vm = myshallowconcat1(vm);
  return s == -2? vecs14(vp, vm): shallowconcat(vp, vm);
}

static GEN
makeS36resolvent(GEN pol, long flag)
{
  GEN R2, V, S = mynfsubfields(pol, 3);
  if (flag < 2) return condrel(gel(S,1), pol, flag);
  R2 = mynfsubfield(pol, 2);
  if (flag == 2)
    V = vec_append(S, R2);
  else
    V = mkvec4(condrel_i(gel(S,1), pol), condrel_i(gel(S,2), pol),
               condrel_i(gel(S,3), pol), condrel_i(R2, pol));
  return V;
}

GEN
nflist_S36_worker(GEN pol, GEN X, GEN Xinf)
{
  GEN d, D = nfcoredisc(pol, &d);
  if (ok_int(mulii(sqri(D), d), X, Xinf)) return makepol6(pol, X2m(d));
  return gen_0;
}

static GEN
parselectS36(GEN v, GEN X, GEN Xinf)
{
  GEN w = nflist_parapply("_nflist_S36_worker", mkvec2(X, Xinf), v);
  long l = lg(w), i, c;

  for (i = c = 1; i < l; i++)
  {
    GEN t = gel(w, i);
    if (typ(t) == t_POL) gel(w, c++) = t;
  }
  setlg(w, c); return w;
}

static GEN
makeS36vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v;

  if (s == 1 || s == 2) return NULL;
  if (s == 3) s = 1;
  if (field)
  {
    if (degpol(field) == 3)
    {
      GEN d, D = nfcoredisc(field,&d);
      long ss = signe(D) < 0? 1: 0;
      if (s >= 0 && s != ss) return NULL;
      if (abscmpii(mulii(sqri(D), d), X) > 0) return NULL;
      v = mkvec(field);
    }
    else
    {
      GEN D2a = absi_shallow(checkfield(field, 2)), D2a3 = powiu(D2a, 3), RES;
      long Fsup, Finf, F, c;
      if ((s >= 0 && s != pol2s(field)) || cmpii(D2a3, X) > 0) return NULL;
      Fsup = floorsqrtndiv(X, D2a3, 4);
      Finf = ceilsqrtndiv(Xinf, D2a3, 4);
      RES = cgetg(Fsup + 1, t_VEC);
      for (F = Finf, c = 1; F <= Fsup; F++)
      {
        pari_sp av = avma;
        GEN w, N = mulii(powuu(F, 4), D2a3);
        if (!(w = makeS36(N, field, s))) set_avma(av);
        else gel(RES, c++) = gerepilecopy(av, w);
      }
      setlg(RES,c); return myshallowconcat1(RES);
    }
  }
  else
    if (!(v = makeS3vec(sqrti(divis(X, 3)), gen_1, NULL, s))) return NULL;
  if (s != -2) return parselectS36(v, X, Xinf);
  return mkvec4(parselectS36(gel(v,1), X, Xinf), cgetg(1, t_VEC),
                cgetg(1, t_VEC), parselectS36(gel(v,2), X, Xinf));
}
/**********************************************************************/
/*                              D612                                  */
/**********************************************************************/
static void
gets2s3(long s, long *s2, long *s3)
{
  switch (s)
  {
    case 0: *s2 = *s3 = 0; break;
    case 2: *s2 = 0; *s3 = 1; break;
    case 3: *s2 = 1; *s3 = -1; break;
    default: *s2 = *s3 = -1; break;
  }
}

static GEN makeD612vec(GEN X, GEN Xinf, GEN field, long s);
static GEN
makeD612(GEN N, GEN field, long s)
{
  long i, j, l, c3, s2, s3;
  GEN v;

  if (s == 1) return NULL;
  gets2s3(s, &s2, &s3);
  if (field)
  {
    GEN D2;
    long si;
    if (degpol(field) == 3) return makeD612vec(N,N,field,s);
    D2 = checkfield(field, 2); si = signe(D2);
    if ((si == 1 && s2 > 0) || (si == -1 && !s2)
        || !divissquare(N, powiu(D2,3))) return NULL;
    v = mkvec(D2);
  }
  else v = divisorsdisc(cored(N, 3), s2);
  l = lg(v);
  for (i = c3 = 1; i < l; i++)
  {
    GEN D2 = gel(v, i), D2a = absi_shallow(D2), M = divii(N, powiu(D2a, 3));
    GEN P2, F = gel(core2(M), 2), L = divisors(mulii(F, D2a));
    long c2, lL = lg(L);
    if (lL == 1) continue;
    P2 = quadpoly_i(D2);
    for (j = c2 = 1; j < lL; j++)
    {
      GEN w, D3 = gel(L, j);
      long k, c, lw;
      if (Mod4(D3) == 2 || !dvdii(F, divii(D3, gcdii(D2a, D3)))
          || !(w = makeDL(3, D3, NULL, s3))) continue;
      lw = lg(w);
      for (k = c = 1; k < lw; k++)
      {
        GEN P3 = gel(w, k), P6, d;
        (void)nfcoredisc(P3, &d); if (equalii(d, D2)) continue;
        if ((P6 = ZX_red_disc(makepol6(P3, P2), N))) gel(w, c++) = P6;
      }
      if (c > 1) { setlg(w, c); gel(L, c2++) = w; }
    }
    if (c2 > 1) { setlg(L, c2); gel(v, c3++) = shallowconcat1(L); }
  }
  setlg(v, c3); return sturmseparate(myshallowconcat1(v), s, 6);
}

static GEN
makeD612resolvent(GEN pol, long flag)
{
  GEN R3, R = mynfsubfield(pol, 2);
  if (flag < 2) return condrel(R, pol, flag);
  R3 = mynfsubfield(pol, 3);
  if (flag == 3) { R = condrel_i(R, pol); R3 = condrel_i(R3, pol); }
  return mkvec2(R, R3);
}

GEN
nflist_D612_worker(GEN P3, GEN X, GEN Xinf, GEN limd2s2)
{
  pari_sp av = avma;
  GEN v, D2, D3 = nfcoredisc(P3, &D2), D32 = sqri(D3), Q = divii(X, D32);
  long limD2 = limd2s2[1], s2 = limd2s2[2];
  long c, M, limD = cmpis(Q, limD2) < 0 ? itos(Q) : limD2;
  v = cgetg(2 * limD + 1, t_VEC);
  for (M = 3, c = 1; M <= limD; M++)
  {
    GEN N, LD = cgetg(1, t_VEC);
    long g, i;
    int p, m;
    uis_fundamental_pm(M, s2, &p, &m);
    if (absequaliu(D2, M))
    { if (signe(D2) > 0) p = 0; else m = 0; }
    if (!(LD = ufund_pm(M, p, m))) continue;
    g = ugcdiu(D3, M);
    N = mulii(D32, muliu(sqru(M/g), M));
    if (cmpii(N, X) <= 0 && cmpii(shifti(N, 2), Xinf) >= 0)
    {
      long l = lg(LD);
      for (i = 1; i < l; i++)
      {
        GEN P = makepol6(P3, X2m(gel(LD,i)));
        if (odd(g)) gel(v, c++) = polredabs(P);
        else if ((P = ZX_red_disc2(P, Xinf, X))) gel(v, c++) = P;
      }
    }
  }
  setlg(v, c); return gerepilecopy(av, v);
}

static GEN
makeD612vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, T;
  long s2, s3;

  if (s == 1) return NULL;
  v = NULL; gets2s3(s, &s2, &s3);
  if (field)
  {
    if (degpol(field) == 3)
    {
      GEN D3 = nfdisc(field);
      long si = signe(D3);
      if ((si > 0 && s2 > 0) || (si < 0 && !s2)
          || cmpii(sqri(D3), X) > 0) return NULL;
      v = mkvec(field);
    }
    else
    {
      GEN D2a = absi_shallow(checkfield(field, 2));
      long l, j, c;
      if (!(v = makeS3vec(sqrti(divii(X, D2a)), gen_1, NULL, s3))) return NULL;
      l = lg(v);
      for (j = c = 1; j < l; j++)
      {
        GEN P = makepol6(gel(v, j), field);
        if ((P = ZX_red_disc2(P, Xinf, X))) gel(v, c++) = P;
      }
      setlg(v, c); return sturmseparate(v, s, 6);
    }
  }
  else if (!(v = makeS3vec(sqrti(X), gen_1, NULL, s3))) return NULL;
  T = mkvecsmall2(floorsqrtn(X, 3), s2);
  v = nflist_parapply("_nflist_D612_worker", mkvec3(X, Xinf, T), v);
  return sturmseparate(myshallowconcat1(v), s, 6);
}

/**********************************************************************/
/*                          A46 and S46P                              */
/**********************************************************************/

/* A46, S46P, in place */
static GEN
makeS46Ppols(long card, GEN v)
{
  long l = lg(v), i;
  GEN d = utoipos(card);
  for (i = 1; i < l; i++)
  {
    GEN G = galoissplittinginit(gel(v,i), d), g = gal_get_gen(G);
    GEN p = (card == 12)? gel(g, 1): mkvec2(gel(g, 1), gel(g, 4));
    gel(v,i) = polredabs(galoisfixedfield(G, p, 1, 0));
  }
  return v;
}
/* S46M, in place */
static GEN
makeS46Mpols(GEN v, GEN X, GEN Xinf)
{
  long l = lg(v), i, c;
  GEN d = utoipos(24);
  for (i = c = 1; i < l; i++)
  {
    GEN G = galoissplittinginit(gel(v,i), d), g = gal_get_gen(G);
    GEN p = perm_mul(gel(g, 4), gel(g, 2));
    p = galoisfixedfield(G, p, 1, 0);
    p = Xinf? ZX_red_disc2(p, Xinf, X): ZX_red_disc(p, X);
    if (p) gel(v, c++) = p;
  }
  setlg(v, c); return v;
}

static GEN
makeA46(GEN N, GEN field, long s)
{
  GEN n, v, D;
  long i, l, c;

  if (s== 1 || s==3 || !Z_issquareall(N, &n)) return NULL;
  if (field)
  {
    GEN t, q, D = checkfield(field, 3);
    if (!Z_issquare(D)
        || !(q = divide(n, D)) || !(t = makeA4S4(field, q, s))) return NULL;
    return makeS46Ppols(12, t);
  }
  D = divisors(gel(core2(n), 2));
  l = lg(D); v = cgetg(l, t_VEC);
  for (i = 2, c = 1; i < l; i++)
  {
    GEN t, q, g3 = gel(D,i), C = makeC3_f(g3);
    long j, l = lg(C);
    if (l == 1) continue;
    q = diviiexact(n, sqri(g3));
    for (j = 1; j < l; j++)
      if ((t = makeA4S4(gel(C,j), q, s))) gel(v, c++) = makeS46Ppols(12,t);
  }
  setlg(v,c); return sturmseparate(myshallowconcat1(v), s, 6);
}

static GEN
makeS46P(GEN N, GEN field, long s)
{
  GEN n, v, D;
  long i, snew, l, c;

  if (s==1 || s==3 || !Z_issquareall(N, &n)) return NULL;
  /* s = -2, -1, 0, 2 */
  if (field)
  {
    GEN D3 = checkfield(field, 3), f, t;
    if (Z_issquare(D3) || !dvdii(n, D3)) return NULL;
    snew = s == 2 && signe(D3) < 0 ? 1 : s;
    f = divii(n, absi_shallow(D3));
    if (!(t = makeA4S4(field, f, snew))) return NULL;
    return makeS46Ppols(24, t);
  }
  D = divisors(n); l = lg(D); v = cgetg(l, t_VEC);
  for (i = 2, c = 1; i < l; i++)
  {
    GEN f, P, D3a = gel(D,i);
    long c3, j, lv3;
    if (!(P =  makeDL(3, D3a, NULL, s? -1: 0))) continue;
    f = gel(D, l-i); lv3 = lg(P);
    for (j = c3 = 1; j < lv3; j++)
    {
      GEN T, P3 = gel(P,j);
      long snew = (s == 2 && signe(ZX_disc(P3)) == -1) ? 1 : s;
      if ((T = makeA4S4(P3, f, snew))) gel(P,c3++) = T;
    }
    if (c3 == 1) continue;
    setlg(P, c3); gel(v, c++) = makeS46Ppols(24, shallowconcat1(P));
  }
  setlg(v,c); return sturmseparate(myshallowconcat1(v), s, 6);
}

GEN
nflist_A46S46P_worker(GEN P3, GEN Xinf, GEN sqX, GEN cards)
{
  pari_sp av = avma;
  long card = cards[1], s = cards[2];
  GEN w, F, V, DATA = S4data(P3, s), D3 = S4_get_disc(DATA);
  GEN D3a = absi_shallow(D3);
  long limf = itos(divii(sqX, D3a)), linf = 1, snew, f, i, c;

  if (cmpii(Xinf, sqri(shifti(D3a, 2))) >= 0)
    linf = ceilsqrtdiv(Xinf, sqri(D3));
  snew = s == 2 && signe(D3) < 0 ? 1 : s;
  V = cgetg(limf, t_VEC);
  F = vecfactoru_i(linf, limf);
  for (f = linf, i = c = 1; f <= limf; f++, i++)
    if ((w = A4S4_fa(DATA, gel(F,i), f, snew)))
      gel(V, c++) = makeS46Ppols(card, w);
  setlg(V,c); V = myshallowconcat1(V);
  return gerepilecopy(av, V);
}

static GEN
makeA46S46Pvec(long card, GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, sqX, T;

  if (s == 1 || s == 3) return NULL;
  sqX = sqrti(X);
  if (field)
  {
    GEN D = checkfield(field, 3);
    long fl = Z_issquare(D);
    if ((card == 12 && !fl) || (card == 24 && fl)) return  NULL;
    v = mkvec(field);
  }
  else
    v = card == 12? makeC3vec(sqX, gen_1, NULL, 0)
                  : makeS3vec(sqX, gen_1, NULL, s? -1: 0);
  if (!v) return NULL;
  T = mkvec3(Xinf, sqX, mkvecsmall2(card, s == -2? -1: s));
  v = nflist_parapply("_nflist_A46S46P_worker", T, v);
  return sturmseparate(myshallowconcat1(v), s, 6);
}

/**********************************************************************/
/*                              S46M                                  */
/**********************************************************************/
static GEN
glco46M(GEN F, GEN D2a)
{
  GEN C, F0, D = divisors(D2a);
  long k, i, c, l = lg(D), klim = vali(D2a)? minss(2, vali(F)): 0;
  /* could restrict divisors to multiples of (D2,F)/2^klim */

  F0 = klim? shifti(F, -klim): F;
  C = cgetg((klim+1) * (l-1) + 1, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    GEN g = gcdii(F, gel(D,l-i));
    long v = vali(g);
    if (v) g = shifti(g, -v);
    if (!is_pm1(g) || v > klim) continue;
    /* (F,D[l-i]) = 2^v; if v <= k <= klim, add F*D[i]>>k */
    gel(C, c++) = g = mulii(F0, gel(D,i));
    for (k = v; k < klim; k++) gel(C, c++) = g = shifti(g, 1);
  }
  setlg(C, c); return C;
}

static GEN
doA4S4(GEN field, GEN C, long s)
{
  long l = lg(C), i, c;
  GEN w, v = cgetg(l, t_VEC);
  for (i = c = 1; i < l; i++)
    if ((w = makeA4S4(field, gel(C,i), s))) gel(v, c++) = w;
  setlg(v,c); return myshallowconcat1(v);
}

static GEN
makeS46M(GEN N, GEN field, long s)
{
  GEN v, D, LC, F;
  long i, c, l, snew;

  if (s == 1) return NULL;
  snew = s == 3 ? 1 : maxss(s, -1);
  if (field)
  {
    GEN D3, D2, D2a, t, Dpow;
    checkfield_i(field, 3); D3 = nfcoredisc(field, &D2); D2a = absi_shallow(D2);
    Dpow = mulii(D2a, sqri(D3));
    if ((signe(D3) < 0 && (s == 0 || s == 2))
        || (signe(D3) > 0 && (s == 3 || Z_issquare(D3)))
        || !divissquareall(N, Dpow, &F)) return NULL;
    LC = glco46M(F, D2a);
    t = doA4S4(field, LC, snew); return makeS46Mpols(t, N, NULL);
  }
  D = divisorsdisc(cored(N, 3), snew);
  l = lg(D); v = cgetg(l*l, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    GEN D2 = gel(D, i), D2a = absi_shallow(D2);
    GEN NSD2 = divii(N, powiu(D2a, 3)), NSD4, F;
    long j;
    if (!Z_issquareall(NSD2, &NSD4)) continue;
    F = divisors(cored(NSD2, 4));
    for (j = 1; j < lg(F); j++)
    {
      GEN f2 = sqri(gel(F, j)), P;
      long k, lP;
      if (!(P = makeDL(3, mulii(D2a, f2), NULL, minss(snew, 1)))) continue;
      lP = lg(P); F = divii(NSD4, f2); LC = glco46M(F, D2a);
      for (k = 1; k < lP; k++) gel(P,k) = doA4S4(gel(P,k), LC, snew);
      gel(v, c++) = makeS46Mpols(shallowconcat1(P), N, NULL);
    }
  }
  if (c == 1) return NULL;
  setlg(v,c); return sturmseparate(gtoset_shallow(shallowconcat1(v)), s, 6);
}

GEN
nflist_S46M_worker(GEN P3, GEN X, GEN Xinf, GEN gs)
{
  pari_sp av = avma;
  long s = gs[1], snew = s == 3 ? 1 : s;
  GEN V, DATA = S4data(P3, s), D3 = S4_get_disc(DATA);
  GEN D2a = absi_shallow(coredisc(D3));
  long lim = floorsqrtdiv(X, mulii(sqri(D3), D2a)), f, c;

  V = cgetg(lim + 1, t_VEC);;
  for (f = 1, c = 1; f <= lim; f++)
  {
    GEN C = glco46M(utoipos(f), D2a), t = doA4S4(DATA, C, snew);
    gel(V, c++) = makeS46Mpols(t, X, Xinf);
  }
  setlg(V,c); V = myshallowconcat1(V);
  return gerepileupto(av, gtoset(V));
}

static GEN
makeS46Mvec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v;

  if (s == 1) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 3);
    if (Z_issquare(D)) return NULL;
    v = mkvec(field);
  }
  else
  {
    long s3 = s == 3? 1: (s < 0? -1: 0), l2, i, c;
    GEN v2 = makeC2vec(sqrtnint(X,3), gen_1, NULL, s3);
    if (!v2) return NULL;
    l2 = lg(v2); v = cgetg(l2, t_VEC);
    for (i = c = 1; i < l2; i++)
    {
      GEN w, T = gel(v2, i), D2a = absi_shallow(nfdisc(T));
      if ((w = makeS3vec(sqrti(divii(X, D2a)), gen_1, T, s3))) gel(v, c++) = w;
    }
    setlg(v,c); v = myshallowconcat1(v);
  }
  v = nflist_parapply("_nflist_S46M_worker",
                           mkvec3(X, Xinf, mkvecsmall(s == -2? -1: s)), v);
  return sturmseparate(myshallowconcat1(v), s, 6);
}

/************************************************************************/
/*                                  A462                                */
/************************************************************************/
static GEN
arch0() { return mkvec(mkvec3(gen_0, gen_0, gen_0)); }
static GEN
arch1g() { return mkvec(mkvec3(gen_1, gen_0, gen_0)); }
static GEN
arch1() { return mkvec3(mkvec3(gen_1, gen_0, gen_0),
                        mkvec3(gen_0, gen_1, gen_0),
                        mkvec3(gen_0, gen_0, gen_1)); }
static GEN
arch2g() { return mkvec(mkvec3(gen_0, gen_1, gen_1)); }
static GEN
arch2() { return mkvec3(mkvec3(gen_0, gen_1, gen_1),
                        mkvec3(gen_1, gen_0, gen_1),
                        mkvec3(gen_1, gen_1, gen_0)); }
static GEN
arch3() { return mkvec(mkvec3(gen_1, gen_1, gen_1)); }

static GEN
archA462(long s)
{
  switch (s)
  {
    case 0: return arch0();
    case 1: return arch1g();
    case 2: return arch2g();
    default: return shallowconcat1(mkvec3(arch0(),arch1g(),arch2g()));
  }
}

static int
stable_arch(GEN v)
{
  long i, l = lg(v);
  GEN x = gel(v,1);
  for (i = 2; i < l; i++) if (!equalii(x, gel(v,i))) return 0;
  return 1;
}
/* nf cyclic of prime degree, return a generator of */
static GEN
cycfindaut(GEN nf)
{
  GEN A = galoisconj(nf, NULL);
  return nfgaloismatrix(nf, gel(A, gequalX(gel(A,1))? 2 : 1));
}

static int
isprM(GEN x)
{ return typ(x) == t_MAT && lg(x) == 3; }
static GEN
doA462(GEN bnf, GEN L, GEN Arch, GEN aut, GEN G, GEN GAL)
{
  pari_sp av = avma;
  long c, k, i, m, lA = lg(Arch), l = lg(L);
  int stable0;
  GEN v;
  if (l == 1) return NULL;
  v = cgetg((lA-1) * (l-1) + 1, t_VEC);
  stable0 = !isprM(gel(L,l-1)); /* not implemented for prM */
  for (i = c = 1; i < lA; i++)
  {
    GEN arch = gel(Arch, i);
    int stable = stable0 && stable_arch(arch);
    for (k = 1; k < l; k++)
    {
      GEN R, id = gel(L,k), F = mkvec2(id, arch);
      long cR, lR;
      if (stable && ZM_equal(nfgaloismatrixapply(bnf, aut, id), id))
        R = mybnrclassfield_X(bnf, F, 2, NULL, NULL, G);
      else
        R = mybnrclassfield(bnf, F, 2);
      lR = lg(R);
      for (m = cR = 1; m < lR; m++)
      {
        GEN P = rnfequation(bnf, gel(R, m));
        if (okgal(P, GAL)) gel(R, cR++) = polredabs(P);
      }
      if (cR > 1) { setlg(R, cR); gel(v, c++) = R; }
    }
  }
  if (c == 1) { set_avma(av); return NULL; }
  setlg(v, c); return gtoset_shallow(shallowconcat1(v));
}
static GEN
makeA462(GEN N, GEN field, long s)
{
  GEN v, L, Arch, GAL;
  long i, c, l;

  if (s == 3) return NULL;
  Arch = archA462(s);
  GAL = mkvecsmall3(24, -1, 2);
  if (field)
  {
    GEN D3 = checkfield(field, 3);
    if (!Z_issquare(D3) || !dvdii(N, sqri(D3))) return NULL;
    L = mkvec(field);
  }
  else
  {
    GEN LD = divisors(cored(N, 4));
    L = cgetg(1, t_VEC);
    for (i = 1; i < lg(LD); i++)
    {
      GEN t = makeC3_f(gel(LD,i));
      if (lg(t) > 1) L = shallowconcat(L, t);
    }
  }
  l = lg(L); v = cgetg(l, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    GEN bnf = bnfY(gel(L,i)), nf = bnf_get_nf(bnf), aut = cycfindaut(bnf);
    GEN T, I = ideals_by_norm(nf, divii(N, sqri(nf_get_disc(nf))));
    GEN G = mkvec2(galoisinit(nf, NULL), gen_0);
    if ((T = doA462(bnf, I, Arch, aut, G, GAL))) gel(v, c++) = T;
  }
  if (c == 1) return NULL;
  setlg(v, c); return sturmseparate(shallowconcat1(v), s, 6);
}

GEN
nflist_A462_worker(GEN P3, GEN X, GEN Xinf, GEN Arch, GEN GAL)
{
  pari_sp av = avma;
  GEN bnf = bnfY(P3), aut = cycfindaut(bnf), v, t;
  GEN G = mkvec2(galoisinit(bnf, NULL), gen_0), D2 = sqri(bnf_get_disc(bnf));
  long c, l, j, lim = itos(divii(X, D2)), liminf = itos(ceildiv(Xinf, D2));

  v = ideallist(bnf, lim); l = lg(v);
  for (c = 1, j = liminf; j < l; j++)
    if ((t = doA462(bnf, gel(v,j), Arch, aut, G, GAL))) gel(v,c++) = t;
  setlg(v, c); return gerepilecopy(av, myshallowconcat1(v));
}
static GEN
makeA462vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, GAL;

  if (s == 3) return NULL;
  if (field)
  {
    GEN D3 = checkfield(field, 3);
    if (!Z_issquare(D3) || cmpii(sqri(D3), X) > 0) return NULL;
    v = mkvec(field);
  }
  else if (!(v = makeC3vec(sqrti(X), gen_1, NULL, 0))) return NULL;
  GAL = mkvecsmall3(24, -1, 2);
  v = nflist_parapply("_nflist_A462_worker", mkvec4(X, Xinf, archA462(s), GAL), v);
  return sturmseparate(myshallowconcat1(v), s, 6);
}

/************************************************************************/
/*                                  S3C3                                */
/************************************************************************/

static int
isok3(ulong N)
{
  GEN fa, P, E;
  long v = u_lvalrem(N, 3, &N), i, l;
  if (v == 1 || v >= 4) return 0;
  fa = factoru(N); P = gel(fa, 1); E = gel(fa, 2); l = lg(P);
  for (i = 1; i < l; i++)
    if (P[i] % 3 == 1) { if (E[i] != 1) return 0; }
    else               { if (E[i] != 2) return 0; }
  return 1;
}

static GEN
makeS3C3(GEN N, GEN field, long s)
{
  GEN v, LD, cond;
  long s2, i;

  if (s == 1 || s == 2) return NULL;
  s2 = s == 3 ? 1 : s;
  if (field)
  {
    GEN D = checkfield(field, 2);
    if (!divissquareall(N, powiu(absi_shallow(D), 3), &cond)) return NULL;
    LD = mkvec(D);
  }
  else LD = divisorsdisc(cored(N, 3), s2);
  v = cgetg(1, t_VEC);
  for (i = 1; i < lg(LD); i++)
  {
    GEN L, bnf, nf, D = gel(LD, i);
    long j, k;
    if (!divissquareall(N, powiu(absi_shallow(D), 3), &cond)) continue;
    bnf = bnfY(Y2m(D)); nf = bnf_get_nf(bnf);
    L = ideals_by_norm(nf, cond);
    for (j = 1; j < lg(L); j++)
    {
      GEN R = mybnrclassfield_N(bnf, gel(L,j), N, 3);
      for (k = 1; k < lg(R); k++)
      {
        GEN P = rnfequation(nf, gel(R, k));
        if (okgal1(P, 18)) v = vec_append(v, polredabs(P));
      }
    }
  }
  return sturmseparate(gtoset_shallow(v), s, 6);
}

GEN
nflist_S3C3_worker(GEN D2, GEN X, GEN Xinf)
{
  pari_sp av = avma;
  GEN bnf = bnfY(Y2m(D2)), nf = bnf_get_nf(bnf), aut = cycfindaut(nf);
  GEN G = mkvec2(galoisinit(bnf, NULL), gen_0);
  long f, c, limf = floorsqrtdiv(X, powuu(itou(D2), 3));
  GEN v = ideallist0(nf, limf, 4 | 8);

  for (f = c = 1; f <= limf; f++)
  {
    pari_sp av2;
    long j, k, cL;
    GEN L;

    if (!isok3(f)) continue;
    av2 = avma; L = gel(v, f);
    for (j = cL = 1; j < lg(L); j++)
    {
      pari_sp av3 = avma;
      long stable = gequal(gel(L,j), nfgaloismatrixapply(nf, aut, gel(L,j)));
      GEN R = mybnrclassfield_X(bnf, gel(L,j), 3, X, Xinf, stable? G: NULL);
      long lR = lg(R), cR;
      for (k = cR = 1; k < lR; k++)
      {
        GEN P = rnfequation(nf, gel(R, k));
        if (okgal1(P, 18)) gel(R, cR++) = polredabs(P);
      }
      if (cR == 1) { set_avma(av3); continue; }
      setlg(R, cR); gel(L, cL++) = R;
    }
    if (cL == 1) { set_avma(av2); continue; }
    setlg(L, cL); gel(v, c++) = shallowconcat1(L);
  }
  setlg(v, c); return gerepilecopy(av, gtoset_shallow(myshallowconcat1(v)));
}

static GEN
makeS3C3vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v;

  if (s == 1 || s == 2) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 2);
    v = mkvec(D);
  }
  else
  {
    long lim = floorsqrtn(X, 3), Da, c;
    v = cgetg(2 * lim + 1, t_VEC);
    for (Da = 3, c = 1; Da <= lim; Da++)
    {
      int p, m;
      uis_fundamental_pm(Da, s, &p, &m);
      if (p) gel(v, c++) = utoipos(Da);
      if (m) gel(v, c++) = utoineg(Da);
    }
    if (c == 1) return NULL;
    setlg(v, c);
  }
  v = nflist_parapply("_nflist_S3C3_worker", mkvec2(X, Xinf), v);
  return sturmseparate(myshallowconcat1(v), s, 6);
}

/************************************************************************/
/*                               S462                                   */
/************************************************************************/

static GEN
archS4621(long s)
{
  switch(s)
  {
    case 0: case 1: return cgetg(1, t_VEC);
    case 2: retmkvec(mkvec(gen_0));
    case 3: retmkvec(mkvec(gen_1));
    default:retmkvec2(mkvec(gen_0), mkvec(gen_1));
  }
}

static GEN
archS4623(long s)
{
  switch (s)
  {
    case 0: return arch0();
    case 1: return arch1();
    case 2: return arch2();
    case 3: return arch3();
    default:return shallowconcat1(mkvec4(arch0(),arch1(),arch2(),arch3()));
  }
}

static GEN
makeS462(GEN N, GEN field, long s)
{
  GEN RES = cgetg(1, t_VEC), L, listarch1, listarch3, GAL;
  long i, j, l, m;
  listarch1 = archS4621(s); listarch3 = archS4623(s);
  GAL = mkvecsmall3(48, -1, 1);
  if (field)
  {
    GEN d = checkfield(field, 3);
    if (Z_issquare(d) || !dvdii(N, sqri(d))) return NULL;
    L = mkvec(field);
  }
  else
  {
    GEN T;
    long c;
    L =  divisors(cored(N, 2));
    for (i = c = 1; i < lg(L); i++)
      if ((T = makeDL(3, gel(L,i), NULL, (s == 0 || s == 1) ? 0 : -1)))
        gel(L, c++) = T;
    if (c == 1) return NULL;
    setlg(L, c); L = shallowconcat1(L);
  }
  for (i = 1; i < lg(L); i++)
  {
    GEN bnf = bnfY(gel(L,i)), nf = bnf_get_nf(bnf);
    GEN I = ideals_by_norm(nf, divii(N, sqri(nf_get_disc(nf))));
    GEN Arch = nf_get_r1(nf) == 1 ? listarch1 : listarch3;
    for (j = 1; j < lg(I); j++)
    {
      GEN id = gel(I, j);
      for (l = 1; l < lg(Arch); l++)
      {
        GEN R = mybnrclassfield(bnf, mkvec2(id, gel(Arch, l)), 2);
        for (m = 1; m < lg(R); m++)
        {
          GEN P = rnfequation(bnf, gel(R, m));
          if (okgal(P, GAL) && (P = ZX_red_disc(P, N))) RES = vec_append(RES, P);
        }
      }
    }
  }
  return sturmseparate(gtoset_shallow(RES), s, 6);
}

GEN
nflist_S462_worker(GEN P3, GEN X, GEN Xinf, GEN vArch, GEN GAL)
{
  pari_sp av = avma;
  GEN bnf = bnfY(P3), nf = bnf_get_nf(bnf), D2 = sqri(nf_get_disc(nf));
  long limf = itos(divii(X, D2)), liminf = itos(ceildiv(Xinf, D2));
  long r1 = nf_get_r1(nf), c, j, k, l, m;
  GEN v, vI = ideallist(bnf, limf), Arch = gel(vArch, r1 == 1? 1 : 2);

  v = cgetg(limf + 1, t_VEC);
  for (c = 1, j = liminf; j <= limf; j++)
  {
    GEN I = gel(vI, j), REU = cgetg(1, t_VEC);
    for (k = 1; k < lg(I); k++)
    {
      GEN id = gel(I, k);
      for (l = 1; l < lg(Arch); l++)
      {
        GEN R = mybnrclassfield(bnf, mkvec2(id, gel(Arch, l)), 2);
        for (m = 1; m < lg(R); m++)
        {
          GEN P = rnfequation(bnf, gel(R, m));
          if (okgal(P, GAL)) REU = vec_append(REU, polredabs(P));
        }
      }
    }
    if (lg(REU) > 1) gel(v, c++) = REU;
  }
  setlg(v,c); v = myshallowconcat1(v);
  return gerepilecopy(av, gtoset_shallow(v));
}
static GEN
makeS462vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, T, GAL;

  if (field)
  {
    GEN D3 = checkfield(field, 3);
    long si = signe(D3);
    if (Z_issquare(D3) || (si < 0 && (s == 0 || s == 1))) return NULL;
    v = mkvec(field);
  }
  else if (!(v = makeS3vec(sqrti(X), gen_1, NULL, (s==0 || s==1)? 0: -1)))
    return NULL;
  GAL = mkvecsmall3(48, -1, 1);
  T = mkvec4(X, Xinf, mkvec2(archS4621(s), archS4623(s)), GAL);
  v = nflist_parapply("_nflist_S462_worker", T, v);
  return sturmseparate(myshallowconcat1(v), s, 6);
}
/************************************************************************/
/*                              C32C4                                   */
/************************************************************************/
static GEN
doC32C4_i(GEN bnf, GEN L, GEN GAL)
{
  long i, l = lg(L);
  GEN v;
  if (l == 1) return L;
  v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN w = cgetg(1, t_VEC), R = mybnrclassfield(bnf, gel(L,i), 3);
    long j, lR = lg(R);
    for (j = 1; j < lR; j++)
    {
      GEN P12 = rnfequation(bnf, gel(R, j)), S = _nfsubfields(P12, 6);
      long k, lS = lg(S);
      for (k = 1; k < lS; k++)
      {
        GEN P = gel(S,k);
        if (okgal(P, GAL)) w = vec_append(w, polredabs(P));
      }
    }
    gel(v,i) = gtoset_shallow(w);
  }
  return myshallowconcat1(v);
}
static GEN
doC32C4(GEN N, GEN P4, GEN GAL)
{
  GEN nf, bnf, F, F2, D4 = nfdisc(P4), D2 = nfdisc(_nfsubfields1(P4, 2));
  if (!(F2 = divide(N, mulii(D2,D4))) || !Z_issquareall(F2, &F)) return NULL;
  bnf = bnfY(P4); nf = bnf_get_nf(bnf);
  return doC32C4_i(bnf, ideals_by_norm(nf, F2), GAL);
}
static GEN
makeC32C4_i(GEN N, GEN field, long s)
{
  GEN GAL = mkvecsmall3(36, 1, 1), v, w, C;
  long c, i, j, l;
  if (!Z_issquare(N) || s == 1 || s == 3) return NULL;
  if (field)
  {
    checkfield_i(field, 4);
    return (okgal2(field,4,-1) && ok_s(field, s))? doC32C4(N, field, GAL): NULL;
  }
  v = divisors(N); l = lg(v);
  for (i = c = 1; i < l; i++)
  {
    long cw, lC;
    if (!(C = makeC4(gel(v, i), NULL, maxss(s, -1)))) continue;
    lC = lg(C);
    for (j = cw = 1; j < lC; j++)
      if ((w = doC32C4(N, gel(C,j), GAL))) gel(C,cw++) = w;
    if (cw > 1) { setlg(C, cw); gel(v, c++) = shallowconcat1(C); }
  }
  setlg(v, c); return myshallowconcat1(v);
}
static GEN
makeC32C4(GEN N, GEN field, long s)
{
  GEN v = makeC32C4_i(N, field, s);
  return v? sturmseparate(v, s, 6): NULL;
}

static GEN
makeC32C4resolvent(GEN pol, long flag)
{
  GEN P12 = polredabs(gel(compositum(pol, pol), 2));
  return condrel(mynfsubfield(P12,4), P12, flag);
}

/* ideals of square norm < lim^2 */
static GEN
ideallistsquare(GEN bnf, long lim)
{
  pari_sp av = avma;
  GEN nf = bnf_get_nf(bnf), V, Z, F;
  long d = nf_get_degree(nf), lim2 = lim * lim, p;
  forprime_t T;

  if (lim <= 0) return cgetg(1, t_VEC);
  V = const_vec(lim, cgetg(1, t_VEC)); gel(V, 1) = mkvec(trivial_fact());
  u_forprime_init(&T, 2, lim);
  F = cgetg(d+1, t_VECSMALL);
  Z = cgetg(d+1, t_VECSMALL);
  while ((p = u_forprime_next(&T)))
  {
    long lv, i, llp = ulogint(lim2, p), tot, m;
    GEN P = idealprimedec_limit_f(nf, utoipos(p), llp);
    GEN W = shallowcopy(V);
    lv = lg(P);
    for (i = tot = 1; i < lv; i++)
    {
      F[i] = pr_get_f(gel(P,i));
      Z[i] = llp / F[i] + 1; tot *= Z[i];
    }
    for (m = 1; m < tot; m++)
    {
      GEN v = cgetg(lv, t_VECSMALL);
      long n = m, S = 0;
      for (i = 1; i < lv; i++) { v[i] = n % Z[i]; n /= Z[i]; S += v[i] * F[i]; }
      if (!odd(S) && S <= llp)
      {
        GEN id = famat_remove_trivial(mkvec2(P, zc_to_ZC(v)));
        long j, pS = upowuu(p, S >> 1);
        for (j = 1; j <= lim / pS; j++)
        {
          GEN vs = shallowcopy(gel(V, j));
          long k, l = lg(vs);
          for (k = 1; k < l; k++) gel(vs, k) = famat_mul(gel(vs, k), id);
          gel(W, pS * j) = shallowconcat(gel(W, pS * j), vs);
        }
      }
    }
    V = W;
  }
  return gerepilecopy(av, V);
}

GEN
nflist_C32C4_worker(GEN P4, GEN X, GEN Xinf, GEN GAL)
{
  pari_sp av = avma;
  GEN bnf = bnfY(P4), D4 = bnf_get_disc(bnf), D2 = nfdisc(_nfsubfields1(P4, 2));
  GEN vI, v, w, D4D2 = mulii(D4, D2);
  long f, c, limf = floorsqrtdiv(X, D4D2), liminf = ceilsqrtdiv(Xinf, D4D2);

  vI = ideallistsquare(bnf, limf); v = cgetg(limf + 1, t_VEC);
  for (c = 1, f = liminf; f <= limf; f++)
    if ((w = doC32C4_i(bnf,  gel(vI, f), GAL))) gel(v, c++) = w;
  setlg(v,c); return gerepilecopy(av, gtoset_shallow(myshallowconcat1(v)));
}
static GEN
makeC32C4vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v, L, GAL;

  if (s == 1 || s == 3) return NULL;
  GAL = mkvecsmall3(36, 1, 1);
  if (field)
  {
    checkfield_i(field, 4);
    if (!okgal2(field, 4, -1) || !ok_s(field, s)) return NULL;
    L = mkvec(field);
  }
  else L = makeC4vec(divis(X, 5), gen_1, NULL, s == -2? -1: s);
  v = nflist_parapply("_nflist_C32C4_worker", mkvec3(X, Xinf, GAL), L);
  return sturmseparate(myshallowconcat1(v), s, 6);
}
/************************************************************************/
/*                                 C9                                   */
/************************************************************************/

static GEN
bnrcfC9(GEN bnf, GEN P, GEN F)
{
  GEN v, cond = F, vec9 = mkvec(utoipos(9)), nf = bnf_get_nf(bnf);
  long i, l, c, lP = lg(P);
  for (i = 1; i < lP; i++)
  {
    GEN p = gel(P, i), pr = idealprimedec_galois(nf, p);
    if (equaliu(p, 3)) pr = idealsqr(nf, pr);
    cond = idealmul(nf, cond, pr);
  }
  v = mybnrclassfield(bnf, cond, 3);
  l = lg(v); if (l == 1) return v;
  for (i = c = 1; i < l; i++)
  {
    GEN P = rnfequation(nf, gel(v,i)), G = galoisinit(P, NULL);
    if (typ(G) != t_INT && gequal(galoisisabelian(G, 2), vec9))
      gel(v, c++) = polredabs(P);
  }
  setlg(v, c); return gtoset_shallow(v);
}

static GEN
makeC9(GEN N, GEN field, long s)
{
  GEN v, D, F;
  long i, lD;

  if (s > 0) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 3), d, P;
    if (!Z_issquareall(D, &d)
        || !divispowerall(N, powiu(D,4), 6, &F)) return NULL;
    P = gel(Z_factor(d), 1);
    return bnrcfC9(bnfY(field), P, F);
  }
  v = cgetg(1, t_VEC);
  D = divisors(cored(N, 8)); lD = lg(D);
  for (i = 2; i < lD; i++)
  {
    GEN v3, P, d = gel(D,i);
    long j, l3;
    if (!Z_ispowerall(divii(N, powiu(d, 8)), 6, &F)
        || !checkcondC3(d, &P)) continue;
    v3 = makeC3_i(d, P); l3 = lg(v3);
    for (j = 1; j < l3; j++)
      v = shallowconcat(v, bnrcfC9(bnfY(gel(v3,j)), P, F));
  }
  return s == -2? vecs(5, v): v;
}

GEN
nflist_C9_worker(GEN T, GEN X, GEN Xinf)
{
  pari_sp av = avma;
  GEN bnf = bnfY(T), D3 = bnf_get_disc(bnf), D34 = powiu(D3, 4);
  GEN sqD = sqrti(D3), P = gel(Z_factor(sqD), 1), v;
  long fl = umodiu(D3, 3) == 0;
  long limf = floorsqrtndiv(X, D34, 6), f, c;
  long limi = ceilsqrtndiv(Xinf, D34, 6);

  v = cgetg(limf + 1, t_VEC); c = 1;
  for (f = limi; f <= limf; f++)
  {
    GEN t;
    if (fl) { long r = f % 9; if (r != 3 && r != 6) continue; }
    t = bnrcfC9(bnf, P, utoipos(f));
    if (lg(t) > 1) gel(v, c++) = t;
  }
  if (c == 1) { set_avma(av); return cgetg(1, t_VEC); }
  setlg(v,c); return gerepilecopy(av, myshallowconcat1(v));
}

static GEN
makeC9vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN v;
  if (s > 0) return NULL;
  if (field)
  {
    GEN D = checkfield(field, 3);
    if (!Z_issquare(D) || cmpii(powiu(D,4), X) > 0) return NULL;
    v = mkvec(field);
  }
  else if (!(v = makeC3vec(sqrtnint(X, 4), gen_1, NULL, 0))) return NULL;
  v = nflist_parapply("_nflist_C9_worker", mkvec2(X, Xinf), v);
  v = myshallowconcat1(v);
  return (s == -2)? vecs(5, v): v;
}
/************************************************************************/
/*                                C3xC3                                 */
/************************************************************************/

static GEN
makeC3C3(GEN N, GEN field, long s)
{
  GEN D, v, f, L;
  long i, j, l, c;

  if (s > 0 || !Z_ispowerall(N, 6, &f)) return NULL;
  D = divisors(f); l = lg(D);
  if (field)
  {
    GEN d = checkfield(field, 3), g;
    if (!Z_issquareall(d, &g) || !dvdii(f, g)) return NULL;
    v = cgetg(l, t_VEC);
    for (i = 2, c = 1; i < l; i++)
    {
      GEN t, g3 = gel(D, i);
      long lt;
      if (equalii(g3, g) || !equalii(lcmii(g,g3), f)) continue;
      t = makeC3_f(g3); lt = lg(t); if (lt == 1) continue;
      for (j = 1; j < lt; j++)
        gel(t,j) = polredabs(polcompositum0(field, gel(t,j), 2));
      gel(v, c++) = t;
    }
    setlg(v, c); return gtoset_shallow(myshallowconcat1(v));
  }
  L = const_vec(l-1, NULL);
  v = cgetg(l * (l-1) / 2 + 1, t_VEC);
  for (i = c = 1; i < l; i++)
  {
    GEN g = gel(D,i);
    for (j = i; j < l; j++)
      if (equalii(lcmii(g, gel(D,j)), f))
      {
        GEN Li, Lj, w;
        long li, lj, a, b, cw;
        if (!gel(L,i)) gel(L,i) = makeC3_f(g);
        if (!gel(L,j)) gel(L,j) = makeC3_f(gel(D,j));
        Li = gel(L,i); li = lg(Li);
        Lj = gel(L,j); lj = lg(Lj); w = cgetg(li * lj, t_VEC);
        for (a = cw = 1; a < li; a++)
          for (b = i == j? a+1: 1; b < lj; b++)
            gel(w, cw++) = polredabs(polcompositum0(gel(Li,a), gel(Lj,b), 2));
        setlg(w, cw); gel(v, c++) = w;
      }
  }
  setlg(v, c); v = gtoset_shallow(myshallowconcat1(v));
  return s == -2? vecs(5, v): v;
}

static GEN
makeC3C3resolvent(GEN pol, long flag)
{
  GEN V = mynfsubfields(pol, 3);
  if (lg(V) != 5) pari_err_BUG("makeC3C3resolvent");
  if (flag < 2) return condrel(gel(V,1), pol, flag);
  if (flag == 2) return V;
  return mkvec4(condrel_i(gel(V,1), pol),
                condrel_i(gel(V,2), pol),
                condrel_i(gel(V,3), pol),
                condrel_i(gel(V,4), pol));
}

/* x, y > 0 */
static GEN
lcmiu(GEN x, ulong y) { return muliu(x, y / ugcd(umodiu(x,y), y)); }
static GEN
lcmuu(ulong x, ulong y) { return muluu(x, y / ugcd(x, y)); }

GEN
nflist_C3C3_worker(GEN gi, GEN w, GEN F, GEN X)
{
  pari_sp av = avma;
  long c, j, i = itos(gi), l = lg(w), f = F[i], x = X[1], xinf = X[2];
  GEN P3 = gel(w, i), v = cgetg(l, t_VEC);
  for (j = i + 1, c = 1; j < l; j++)
    if (ok_intu(lcmuu(f, F[j]), x, xinf))
      gel(v, c++) = polredabs(polcompositum0(P3, gel(w, j), 2));
  setlg(v, c); return gerepilecopy(av, v);
}

static GEN
makeC3C3vec(GEN X, GEN Xinf, GEN field, long s)
{
  GEN F, v, v3;
  long j, l, x, xinf;

  if (s > 0) return NULL;
  x = floorsqrtn(X, 6);
  v3 = C3vec_F(x, 1, &F); if (!v3) return NULL;
  v3 = zvV_to_ZXV(v3); l = lg(v3);
  v = cgetg((l - 1) * l / 2 + 1, t_VEC);
  xinf = ceilsqrtn(Xinf, 6);
  if (field)
  {
    GEN F3, D3 = checkfield(field, 3);
    long c;
    if (!Z_issquareall(D3, &F3)) return NULL;
    for (j = c = 1; j < l; j++)
      if (ok_intu(lcmiu(F3, F[j]), x, xinf) && !ZX_equal(gel(v3,j), field))
        gel(v, c++) = polredabs(polcompositum0(field, gel(v3,j), 2));
    setlg(v, c);
  }
  else
  {
    GEN T = mkvec3(v3, F, mkvecsmall2(x,xinf));
    v = nflist_parapply("_nflist_C3C3_worker", T, identity_ZV(l-1));
    v = myshallowconcat1(v);
  }
  v = gtoset_shallow(v); return s == -2? vecs(5, v): v;
}

/************************************************************************/
/*                                S32                                   */
/************************************************************************/

static GEN
makepolS32(GEN P1, GEN P2)
{
  GEN G = galoissplittinginit(polcompositum0(P1, P2, 2), utoipos(36));
  GEN vH = galoissubgroups(G), g = mkvec2(gal_get_gen(G), gal_get_orders(G));
  long i, l = lg(vH);
  for (i = 1; i < l; i++)
  {
    GEN H = gel(vH, i);
    if (group_order(H) == 6 && !group_isabelian(H) /*S3*/
        && group_subgroup_is_faithful(g, H))
      return polredabs(galoisfixedfield(G, H, 1, 0));
  }
  return NULL; /*LCOV_EXCL_LINE*/
}

static GEN
extractS3cond(GEN V3, GEN sqX, GEN field, long s)
{
  GEN v, v2 = NULL;
  long l = lg(V3), c, c2, i;

  v = cgetg(l, t_VEC);
  if (s == 3) v2 = cgetg(l, t_VEC);
  for (i = c = c2 = 1; i < l; i++)
  {
    GEN pol = gel(V3, i), D, F, DF;
    (void)nfcoredisc2(pol, &D, &F); DF = mulii(D, F);
    if (abscmpii(DF, sqX) <= 0)
    {
      GEN ind = field || s == 3 ? gen_0 : utoipos(c);
      GEN V = mkvecn(5, pol, F, mulii(sqri(DF), D), D, ind);
      if (s != 3 || signe(D) > 0) gel(v, c++) = V; else gel(v2, c2++) = V;
    }
  }
  setlg(v, c); if (s != 3) return v;
  setlg(v2, c2); return mkvec2(v, v2);
}

static GEN makeS32common(GEN V3, GEN X, GEN Xinf, GEN field, long s);
static GEN
makeS32(GEN N, GEN field, long s)
{
  long s3, i, c, l;
  GEN v, t;

  if (s == 1) return NULL;
  s3 = -1; if (s == 0) s3 = 0; if (s == 2) s3 = 1;
  v = divisors(N); l = lg(v);
  for (i = 2, c = 1; i < l; i++)
    if ((t = makeDL(3, gel(v, i), NULL, s3))) gel(v,c++) = t;
  setlg(v,c); return makeS32common(myshallowconcat1(v), N, N, field, s);
}

static GEN
group_add_elt(GEN H, GEN g, long r)
{ return mkvec2(vec_append(gel(H,1),g), vecsmall_append(gel(H,2), r)); }

static GEN
makeS32resolvent(GEN pol, long flag)
{
  GEN w, g1, g2, H1, H2, G = galoissplittinginit(pol, utoipos(36));
  GEN v = galoissubgroups(G), g = galois_group(G);
  long i, c, l = lg(v);
  for (i = c = 1; i < l; i++)
  {
    GEN H = gel(v,i);
    if (group_order(H) == 6 && group_subgroup_isnormal(g,H)) gel(v, c++) = H;
  }
  H1 = gel(v,1); g1 = gel(H1,1);
  H2 = gel(v,2); g2 = gel(H2,1); /* G = H1 x H2, Hi ~ S3 */
  H1 = group_add_elt(H1, gel(g2,2), 2);
  H2 = group_add_elt(H2, gel(g1,2), 2);
  w = condrel_dummy(galoisfixedfield(G,H1,1,0), flag);
  if (flag >= 2) w = mkvec2(w, condrel_dummy(galoisfixedfield(G,H2,1,0),flag));
  return w;
}

/* s = 0: real, real; s = 1 imp; s = 2: imag, imag; s = 3: real, imag. */
GEN
nflist_S32_worker(GEN S1, GEN X, GEN Xinf, GEN w, GEN gs)
{
  pari_sp av = avma;
  GEN pol1 = gel(S1, 1), F1 = gel(S1, 2), A1 = gel(S1, 3), D1 = gel(S1, 4), v;
  long c, j, l = lg(w), i = itos(gel(S1, 5)), s = gs[1];
  v = cgetg(l, t_VEC);
  for (j = s == 3 ? 1 : i + 1, c = 1; j < l; j++)
  {
    GEN S2 = gel(w,j), F2 = gel(S2,2), A2 = gel(S2,3), D2 = gel(S2,4), Q,P;
    if (equalii(D2, D1)) continue;
    P = mulii(sqri(gcdii(D1,D2)), gcdii(F1,F2)); /* usually 1 */
    Q = diviiexact(mulii(A1, A2), sqri(P)); if (abscmpii(Q, X) > 0) continue;
    P = makepolS32(pol1, gel(S2,1));
    if (ok_int(nfdisc(P), X, Xinf)) gel(v, c++) = P;
  }
  setlg(v, c); return gerepilecopy(av, v);
}

static GEN
makeS32common(GEN v, GEN X, GEN Xinf, GEN field, long s)
{
  GEN v1, v2;
  v = extractS3cond(v, sqrti(X), field, s);
  if (field)
  {
    GEN D, F;
    long si;
    checkfield_i(field, 3); nfcoredisc2(field, &D, &F); si = signe(D);
    if ((si > 0 && s == 2) || (si < 0 && s == 0) || equali1(D)) return NULL;
    v2 = mkvec(mkvecn(5, field, F, mulii(sqri(F), powiu(D, 3)), D, gen_0));
    if (s != 3) v1 = v; else v1 = gel(v, si > 0 ? 2 : 1);
  }
  else
    if (s != 3) v1 = v2 = v; else { v1 = gel(v, 1); v2 = gel(v, 2); }
  v = nflist_parapply("_nflist_S32_worker", mkvec4(X, Xinf, v2, mkvecsmall(s)), v1);
  return sturmseparate(gtoset_shallow(myshallowconcat1(v)), s, 6);
}

static GEN
makeS32vec(GEN X, GEN Xinf, GEN field, long s)
{
  long s3 = -1;
  GEN v;

  if (s == 1) return NULL;
  if (s == 0) s3 = 0; else if (s == 2) s3 = 1;
  if (!(v = makeS3vec(divis(X, s? 3: 5), gen_1, NULL, s3))) return NULL;
  return makeS32common(v, X, Xinf, field, s);
}

/************************************************************************/
/*                            C32:D4                                    */
/************************************************************************/

static GEN
makeC32D4resolvent(GEN pol, long flag) { return makeC32C4resolvent(pol, flag); }
static int
cyc_is_trivial(GEN c) { return lg(c) == 1 || equali1(gel(c,1)); }

static GEN
C32D4pol(GEN bnf, GEN id)
{
  GEN v, g3 = utoipos(3), bnr = bnrinitmod(bnf, id, 0, g3);
  long l, i, c;

  if (cyc_is_trivial(bnr_get_cyc(bnr))) return NULL;
  v = bnrclassfield(bnr, g3, 0, DEFAULTPREC);
  if (typ(v) == t_POL) v = mkvec(v);
  l = lg(v);
  for (i = c = 1; i < l; i++)
  {
    GEN Q = rnfequation0(bnf, gel(v, i), 0);
    Q = _nfsubfields(Q, 6);
    if (lg(Q) > 1)
    {
      Q = polredabs(gel(Q, 1));
      if (okgal1(Q, 72)) gel(v, c++) = Q;
    }
  }
  if (c == 1) return NULL;
  setlg(v, c); return v;
}

static GEN
bigdisc(GEN P) { return mulii(nfdisc(P), nfdisc(_nfsubfields1(P, 2))); }

/* v,w = factorization matrices (for ideals of the same norm), do we have
 * aut(v) = w ? */
static int
prMconj(GEN nf, GEN v, GEN w, GEN aut)
{
  GEN P = gel(v,1), E = gel(v,2), Q = gel(w,1), F = gel(w,2);
  long i, j, l = lg(P);
  if (lg(Q) != l) return 0;
  if (!ZV_equal(ZV_sort_shallow(E), ZV_sort_shallow(F))) return 0;
  Q = shallowcopy(Q);
  for (i = 1; i < l; i++)
  {
    GEN pr = gel(P,i), p = pr_get_p(pr), e = gel(E,i), pi = pr_get_gen(pr);
    long ep = pr_get_e(pr), fp = pr_get_f(pr);
    pi = nfgaloismatrixapply(nf, aut, pi);
    for (j = 1; j < l; j++)
    {
      GEN qr = gel(Q,j);
      if (!qr) continue;
      if (pr_get_f(qr) == fp && pr_get_e(qr) == ep
          && equalii(gel(F,j), e) && equalii(pr_get_p(qr), p)
          && nfval(nf, pi, qr)) { gel(Q,j) = NULL; break; }
    }
  }
  return 1;
}

GEN
nflist_C32D4_worker(GEN P, GEN X, GEN Xinf, GEN gs)
{
  pari_sp av = avma;
  GEN bd = bigdisc(P), RES = cgetg(1, t_VEC), L, bnf, nf, aut;
  long s = itos(gs), lim, j;

  if (absi_cmp(bd, X) > 0) { set_avma(av); return cgetg(1, t_VEC); }
  bnf = bnfY(P); nf = bnf_get_nf(bnf); aut = cycfindaut(nf);
  lim = itos(divii(X, absi_shallow(bd)));
  L = ideallistsquare(bnf, lim);
  for (j = 1; j <= lim; j++)
  {
    GEN v = gel(L, j);
    long k, lv = lg(v);
    for (k = 1; k < lv; k++)
    {
      GEN R, vk = gel(v, k);
      long m, n, c, lR;
      if (!vk || !(R = C32D4pol(bnf, vk))) continue;
      lR = lg(R);
      for (m = c = 1; m < lR; m++)
      {
        GEN Z = gel(R, m);
        if (ok_s(Z, s) && ok_int(nfdisc(Z), X, Xinf)) gel(R, c++) = Z;
      }
      if (c > 1) { setlg(R, c); RES = shallowconcat(RES, R); }
      for (n = k + 1; n < lv; n++)
        if (gel(v,n) && prMconj(nf, vk, gel(v,n), aut)) {gel(v,n)=NULL; break;}
    }
  }
  return gerepilecopy(av, RES);
}

static GEN
makeC32D4vec(GEN X, GEN Xinf, GEN field, long s)
{
  long s4;
  GEN v;

  if (s == -2) s4 = -1; else if (s == 3) s4 = 2; else s4 = s;
  if (field)
  {
    checkfield_i(field, 4);
    if (!okgal1(field, 8) || !ok_s(field,s4)) return NULL;
    v = mkvec(field);
  }
  else v = makeD4vec(X, gen_1, NULL, s4);
  v = nflist_parapply("_nflist_C32D4_worker", mkvec3(X, Xinf, stoi(s)), v);
  return sturmseparate(gtoset_shallow(myshallowconcat1(v)), s, 6);
}

static GEN
makeC32D4(GEN N, GEN field, long s)
{
  long s4, i, lv;
  GEN v;
  if (s == -2) s4 = -1; else if (s == 3) s4 = 2; else s4 = s;
  if (field)
  {
    GEN D = checkfield(field, 4);
    if (!okgal1(field, 8) || !ok_s(field,s4) || !dvdii(N, D)) return NULL;
    v = mkvec(field);
  }
  else
  {
    long c;
    GEN C;
    v = divisors(absi(N)); lv = lg(v);
    for (i = c = 1; i < lv; i++)
      if ((C = makeD4(gel(v, i), NULL, s4)))
      {
        long j, cC, lC = lg(C);
        for (j = cC = 1; j < lC; j++)
          if (dvdii(N, bigdisc(gel(C,j)))) gel(C, cC++) = gel(C,j);
        if (cC > 1) { setlg(C, cC); gel(v, c++) = C; }
      }
    if (c == 1) return NULL;
    setlg(v, c); v = shallowconcat1(v);
  }
  lv = lg(v);
  for (i = 1; i < lv; i++)
    gel(v,i) = nflist_C32D4_worker(gel(v,i), N, N, stoi(s));
  return sturmseparate(gtoset_shallow(myshallowconcat1(v)), s, 6);
}

/************************************************************************/
/*                         Global Programs                              */
/************************************************************************/
static long
grouptranslate(const char *g, long *t, int QT)
{
  long ell;
  char r;

  if (QT)
  {
    r = *g; ell = itos( strtoi(g + 1) );
    if (ell < 0) return 0;
    if (r == 'A') { *t = -2; return ell; }
    if (r == 'S') { *t = -1; return ell; }
    if (!strcmp(g, "C3")) { *t = -2; return ell; }
  }
  if (!strcmp(g, "C1")) { *t = 1; return 1; }
  if (!strcmp(g, "C2") || !strcmp(g, "D2")) { *t = 1; return 2; }
  if (!strcmp(g, "C3")) { *t = 1; return 3; }
  if (!strcmp(g, "S3") || !strcmp(g,"D3")) { *t = 2; return 3; }
  if (!strcmp(g, "C4")) { *t = 1; return 4; }
  if (!strcmp(g, "V4")) { *t = 2; return 4; }
  if (!strcmp(g, "D4")) { *t = 3; return 4; }
  if (!strcmp(g, "A4")) { *t = 4; return 4; }
  if (!strcmp(g, "S4")) { *t = 5; return 4; }
  if (!strcmp(g, "C5")) { *t = 1; return 5; }
  if (!strcmp(g, "D5"))  { *t = 2; return 5; }
  if (!strcmp(g, "F5") || !strcmp(g, "M20"))  { *t = 3; return 5; }
  if (!strcmp(g, "A5")) { *t = 4; return 5; }
  if (!strcmp(g, "A5cond")) { *t = 9; return 5; }
  if (!strcmp(g, "C6")) { *t = 1; return 6; }
  if (!strcmp(g, "D6")) { *t = 2; return 6; }
  if (!strcmp(g, "C7")) { *t = 1; return 7; }
  if (!strcmp(g, "D7")) { *t = 2; return 7; }
  if (!strcmp(g, "M21")) { *t = 3; return 7; }
  if (!strcmp(g, "M42")) { *t = 4; return 7; }
  if (!strcmp(g, "C9")) { *t = 1; return 9; }
  if (!strcmp(g, "D9")) { *t = 3; return 9; }
  if (QT)
  {
    if (!strcmp(g, "C8")) { *t = 1; return 8; }
    if (!strcmp(g, "D8")) { *t = 2; return 8; }
    if (!strcmp(g, "C10")) { *t = 1; return 10; }
    if (!strcmp(g, "D10")) { *t = 3; return 10; }
    if (!strcmp(g, "C11")) { *t = 1; return 11; }
    if (!strcmp(g, "D11")) { *t = 2; return 11; }
  }
  r = *g; ell = itos( strtoi(g + 1) );
  if (ell >= 8 && uisprime(ell))
  {
    if (r == 'C') { *t = 1; return ell; }
    if (r == 'D') { *t = 2; return ell; }
  }
  *t = 0; return 0;
}
static long
group_nTk(GEN g, long *t, int QT)
{
  long L[] = { 0, /* https://oeis.org/A002106 */
  1,1,2,5,5,16,7,50,34,45,8,301,9,63,104,1954,10,
  983,8,1117,164,59,7,25000,211,96,2392,1854,8,5712,
  12,2801324,162,115,407,121279,11,76,306,315842,10,
  9491,10,2113,10923,56,6 };
  long N = numberof(L), n, k;

  if (lg(g) != 3 || !RgV_is_ZV(g)) { *t = 0; return 0; }
  n = itos(gel(g,1)); if (n <= 0) return 0;
  if (n >= N) pari_err_IMPL(stack_sprintf("group nTk with n > %ld", N-1));
  *t = k = itos(gel(g,2));
  if (k <= 0 || k > L[n])
  {
    char *s;
    s = stack_sprintf("incorrect group %ldTk with k = %ld not in [1,%ld]",
                      n, k, L[n]);
    pari_err(e_MISC, s);
  }
  if (!QT)
  {
    if (n <= 9)
    {
      long v[] = { 0, 1, 1, 2, 5, 4, 13, 4, 0, 3 };
      return k <= v[n]? n: 0;
    }
    return (uisprime(n) && k <= 2)? n: 0;
  }
  if (n <= 2) *t = -2; /* An */
  else if (k == L[n]) *t = -1; /* Sn */
  else if (k == L[n]-1) *t = -2; /* An */
  return n;
}

static int
okfield(GEN F) { return typ(F) == t_POL && RgX_is_ZX(F) && ZX_is_irred(F); }
static GEN
nfmakenum(long n, long t, GEN N, GEN field, long s)
{
  GEN v = NULL;
  switch(100 * n + t)
  {
    case 101: return makeC1(N, field, s);
    case 201: return makeC2(N, field, s);
    case 301: return makeC3(N, field, s);
    case 302: return makeDL(3, N, field, s);
    case 401: return makeC4(N, field, s);
    case 402: return makeV4(N, field, s);
    case 403: return makeD4(N, field, s);
    case 404: return makeA4(N, field, s);
    case 405: return makeS4(N, field, s);
    case 501: return makeC5(N, field, s);
    case 502: return makeDL(5, N, field, s);
    case 503: return makeMgen(5, 4, N, field, s); /*F5*/
    case 504: return makeA5(N, s);
    case 509: return makeA5cond(N, s);
    case 601: return makeC6(N, field, s);
    case 602: return makeS36(N, field, s);
    case 603: return makeD612(N, field, s);
    case 604: return makeA46(N, field, s);
    case 605: return makeS3C3(N, field, s);
    case 606: return makeA462(N, field, s);
    case 607: return makeS46P(N, field, s);
    case 608: return makeS46M(N, field, s);
    case 609: return makeS32(N, field, s);
    case 610: return makeC32C4(N, field, s);
    case 611: return makeS462(N, field, s);
    case 612: return makeA56(N, s);
    case 613: return makeC32D4(N, field, s);
    case 701: return makeCL(7, N, field, s);
    case 702: return makeDL(7, N, field, s);
    case 703: return makeMgen(7, 3, N, field, s); /*M21*/
    case 704: return makeMgen(7, 6, N, field, s);
    case 901: return makeC9(N, field, s);
    case 902: return makeC3C3(N, field, s);
    case 903: return makeD9(N, field, s);
  }
  if (!v && uisprime(n)) switch(t)
  {
    case 1: return makeCL(n, N, field, s);
    case 2: return makeDL(n, N, field, s);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}
/* deg(pol) < 8 */
static GEN
nfresolvent_small(GEN pol, long flag)
{
  long deg = degpol(pol), dP, s;
  GEN G;
  if (deg == 1) return makeC1resolvent(flag);
  if (deg == 2) return makeC2resolvent(pol, flag);
  G = polgalois(pol, DEFAULTPREC);
  dP = itos(gel(G,1));
  if (deg == 3) return dP == 3? makeC3resolvent(pol, flag)
                              : makeS3resolvent(pol, flag);
  s = itos(gel(G,2));
  if (deg == 4)
  {
    if (dP == 4) return s == -1? makeC4resolvent(pol, flag)
                               : makeV4resolvent(pol, flag);
    if (dP == 8) return condrelresolvent(pol, 2, flag); /*D4*/
    return makeA4S4resolvent(pol, flag);
  }
  if (deg == 5)
  {
    if (dP == 5)  return makeCLresolvent(5, pol, flag);
    if (dP == 10) return makeDLresolvent(5, pol, flag);
    if (dP == 20) return makeMgenresolvent(5, 4, pol, flag); /*F5*/
    if (dP == 60) return makeA5resolvent(pol, flag);
  }
  if (deg == 6)
  {
    if (dP == 6 && s == -1)
    { /* works both with new_galois_format set or unset */
      long k = itos(gel(G,3));
      return k == 1? makeC6resolvent(pol, flag)
                   : makeS36resolvent(pol, flag);
    }
    if (dP == 12) return s == -1? makeD612resolvent(pol, flag)
                                : condrelresolvent(pol,3,flag); /*A46*/
    if (dP == 18) return condrelresolvent(pol,2,flag); /*S3C3*/
    if (dP == 24) return condrelresolvent(pol,3,flag); /*S46P,S46M,A462*/
    if (dP == 36) return (s == 1)? makeC32C4resolvent(pol, flag)
                                 : makeS32resolvent(pol, flag);
    if (dP == 48) return condrelresolvent(pol,3,flag); /*S462*/
    if (dP == 60) return makeA56resolvent(pol,flag);
    if (dP == 72) return makeC32D4resolvent(pol, flag);
  }
  if (deg == 7)
  {
    if (dP == 7)  return makeCLresolvent(7, pol, flag);
    if (dP == 14) return makeDLresolvent(7, pol, flag);
    if (dP == 21) return makeMgenresolvent(7, 3, pol, flag); /*M21*/
    if (dP == 42) return makeMgenresolvent(7, 6, pol, flag); /*M42*/
  }
  return gen_0;
}

static GEN
nfresolvent_i(GEN pol, long flag)
{
  long d;
  GEN G;

  if (!okfield(pol)) pari_err_TYPE("nfresolvent", pol);
  if (flag < 0 || flag > 3) pari_err_FLAG("nfresolvent");
  d = degpol(pol);
  if (d < 8) return nfresolvent_small(pol, flag);
  if (d != 9 && !uisprime(d)) return gen_0;
  G = galoisinit(pol, NULL);
  if (typ(G) != t_INT)
  {
    if (d == 9)
    {
      long n = lg(gal_get_gen(G))-1;
      return n == 1? condrelresolvent(pol,3,flag) /*C9*/
                   : makeC3C3resolvent(pol, flag); /*C3xC3*/
    }
    return makeCLresolvent(d, pol, flag);
  }
  G = galoissplittinginit(pol, utoipos(2*d));
  if (gal_get_order(G) != 2*d) return gen_0;
  return d == 9? makeD9resolvent(G, flag): makeDLresolvent(d, pol, flag);
}
GEN
nfresolvent(GEN pol, long flag)
{ pari_sp av = avma; return gerepilecopy(av, nfresolvent_i(pol, flag)); }

/* 1 <= Xinf <= X */
static GEN
nfmakevecnum(long n, long t, GEN X, GEN Xinf, GEN field, long s)
{
  switch(n * 100 + t)
  {
    case 101: return makeC1vec(Xinf, field, s);
    case 201: return makeC2vec(X, Xinf, field, s);
    case 301: return makeC3vec(X, Xinf, field, s);
    case 302: return makeS3vec(X, Xinf, field, s);
    case 401: return makeC4vec(X, Xinf, field, s);
    case 402: return makeV4vec(X, Xinf, field, s);
    case 403: return makeD4vec(X, Xinf, field, s);
    case 404: return makeA4S4vec(1, X, Xinf, field, s);
    case 405: return makeA4S4vec(0, X, Xinf, field, s);
    case 501: return makeC5vec(X, Xinf, field, s);
    case 502: return makeDLvec(5, X, Xinf, field, s);
    case 503: return makeMgenvec(5, 4, X, Xinf, field, s); /*F5*/
    case 504: return makeA5vec(X, Xinf, field, s);
    case 509: return makeA5condvec(X, Xinf, field, s);
    case 601: return makeC6vec(X, Xinf, field, s);
    case 602: return makeS36vec(X, Xinf, field, s);
    case 603: return makeD612vec(X, Xinf, field, s);
    case 604: return makeA46S46Pvec(12, X, Xinf, field, s);/*A46S*/
    case 605: return makeS3C3vec(X, Xinf, field, s);
    case 606: return makeA462vec(X, Xinf, field, s);
    case 607: return makeA46S46Pvec(24, X, Xinf, field, s); /*S46P*/
    case 608: return makeS46Mvec(X, Xinf, field, s);
    case 609: return makeS32vec(X, Xinf, field, s);
    case 610: return makeC32C4vec(X, Xinf, field, s);
    case 611: return makeS462vec(X, Xinf, field, s);
    case 612: return makeA56vec(X, Xinf, s);
    case 613: return makeC32D4vec(X, Xinf, field, s);
    case 701: return makeCLvec(7, X, Xinf, field, s);
    case 702: return makeDLvec(7, X, Xinf, field, s);
    case 703: return makeMgenvec(7, 3, X, Xinf, field, s); /*M21*/
    case 704: return makeMgenvec(7, 6, X, Xinf, field, s); /*M41*/
    case 901: return makeC9vec(X, Xinf, field, s);
    case 902: return makeC3C3vec(X, Xinf, field, s);
    case 903: return makeD9vec(X, Xinf, field, s);
  }
  if (uisprime(n)) switch(t)
  {
    case 1: return makeCLvec(n, X, Xinf, field, s);
    case 2: return makeDLvec(n, X, Xinf, field, s);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

/* s > -2 */
static GEN
nfmakesomehard(long n, long t, long s)
{
  pari_sp av = avma;
  long i;
  for (i = 1;; i++, set_avma(av))
  {
    GEN v = nfmakevecnum(n, t, int2n(18 + 2*i), gen_1, NULL, s);
    if (v && lg(v) > 2) return v;
  }
}
static long
minlim(GEN v)
{
  long i, m = LONG_MAX;
  if (!v) return m;
  for (i = lg(v)-1; i; i--) if (v[i] && m > v[i]) m = v[i];
  return m;
}
static GEN
nfmakesome(long n, long t, long s)
{
  GEN v = NULL;
  long lim, flag = 0;
  switch(n * 100 + t)
  {
    case 101: v = mkvecsmall(1); break;
    case 201: v = mkvecsmall2(33, 24); break;
    case 301: v = mkvecsmall2(3969, 0); break;
    case 302: v = mkvecsmall2(568, 108); break;
    case 401: v = mkvecsmall3(35152, 0, 44217); break;
    case 402: v = mkvecsmall3(14400, 0, 1225); break;
    case 403: v = mkvecsmall3(5125, 1375, 549); break;
    case 404: v = mkvecsmall3(270400, 0, 29241); break;
    case 405: v = mkvecsmall3(8468, 976, 1076); break;
    case 501: v = mkvecsmall3(1073283121, 0, 0); break;
    case 502: v = mkvecsmall3(4330561, 0, 51529); break;
    case 503: v = mkvecsmall3(LONG_MAX, 0, 253125); break;
    case 504: v = mkvecsmall3(11812969, 0, 149769); break;
    case 509: v = mkvecsmall3(5105, 0, 992); break;
    case 601: v = mkvecsmall4(4148928, 0, 0, 2250423); break;
    case 602: v = mkvecsmall4(32166277, 0, 0, 273375); break;
    case 603: v = mkvecsmall4(9045125, 0, 242000, 86528); break;
    case 604: v = mkvecsmall4(125238481, 0, 4439449, 0); break;
    case 605: v = mkvecsmall4(7442000, 0, 0, 143883); break;
    case 606: v = mkvecsmall4(2115281, 419904, 373977, 0); break;
    case 607: v = mkvecsmall4(12730624, 0, 118336, 0); break;
    case 608: v = mkvecsmall4(183250432, 0, 440711081, 13144256); break;
    case 609: v = mkvecsmall4(LONG_MAX, 0, 1382400, 1494108); break;
    case 610: v = mkvecsmall4(765905625, 0, 4950625, 0); break;
    case 611: v = mkvecsmall4(5695040, 941872, 57661, 37479); break;
    case 612: v = mkvecsmall4(185313769, 0, 1907161, 0); break;
    case 613: v = mkvecsmall4(LONG_MAX, 221875, 87625, 44496); break;
    case 701: v = mkvecsmall4(LONG_MAX, 0, 0, 0); break;
    case 702: v = mkvecsmall4(LONG_MAX, 0, 0, 80062991); break;
    case 703: v = mkvecsmall4(LONG_MAX, 0, 0, 0); break;
    case 704: v = mkvecsmall4(LONG_MAX, 0, 0, LONG_MAX); break;
    case 901: v = mkvecsmall5(LONG_MAX, 0, 0, 0, 0); break;
    case 902: v = mkvecsmall5(LONG_MAX, 0, 0, 0, 0); break;
    case 903: v = mkvecsmall5(LONG_MAX, 0, 0, 0, LONG_MAX); break;
  }
  if (!v) flag = uisprime(n) && t <= 2? t: 0;
  if (s == -2)
  {
    long i, l = (n >> 1) + 2;
    GEN W = cgetg(l, t_VEC);
    for (i = 1; i < l; i++)
    {
      GEN w = NULL;
      if (!v)
      { if (i == 1 || (i == l-1 && flag == 2)) w = nfmakesomehard(n, t, i-1); }
      else if (v[i] == LONG_MAX)
        w = nfmakesomehard(n, t, i-1);
      else if (v[i])
        w = nfmakevecnum(n, t, utoipos(v[i]), gen_1, NULL, i-1);
      gel(W, i) = w? w: cgetg(1, t_VEC);
    }
    return W;
  }
  else if (s == -1)
    lim = minlim(v);
  else
  {
    lim = v[s + 1];
    if (!lim) return cgetg(1, t_VEC);
  }
  if (lim == LONG_MAX) return nfmakesomehard(n, t, s);
  return nfmakevecnum(n, t, utoipos(lim), gen_1, NULL, s);
}

GEN
nflist(GEN GP, GEN N, long s, GEN field)
{
  pari_sp av = avma;
  GEN v, X, Xinf;
  long n = 0, t = 0, tp = typ(GP);
  long QT = N && typ(N) == t_POL;

  if (s < -2) pari_err_DOMAIN("nflist", "s", "<", gen_m2, stoi(s));
  if (field && !okfield(field)) pari_err_TYPE("nflist", field);
  switch(tp)
  {
    case t_STR: n = grouptranslate(GSTR(GP), &t, QT); break;
    case t_VEC: n = group_nTk(GP, &t, QT); break;
  }
  if (!n)
  {
    const char *s =
    "unsupported group (%Ps). Use one of\n\
  \"C1\"=[1,1];\n\
  \"C2\"=[2,1];\n\
  \"C3\"=[3,1], \"S3\"=[3,2];\n\
  \"C4\"=[4,1], \"V4\"=[4,2], \"D4\"=[4,3], \"A4\"=[4,4], \"S4\"=[4,5];\n\
  \"C5\"=[5,1], \"D5\"=[5,2], \"F5\"=\"M20\"=[5,3], \"A5\"=[5,4];\n\
  \"C6\"=[6,1], \"D6\"=[6,2], [6,3], [6,4],..., [6,13];\n\
  \"C7\"=[7,1], \"D7\"=[7,2], \"M21\"=[7,3], \"M42\"=[7,4];\n\
  \"C9\"=[9,1], [9,2], \"D9\"=[9,3].\"\n\
  Also supported are \"Cp\"=[p,1] and \"Dp\"=[p,2] for any odd prime p";
    pari_err(e_MISC, s, GP);
  }
  if (QT) return gerepilecopy(av, nflistQT(n, t, varn(N)));
  if (s > (n >> 1)) return cgetg(1, t_VEC);
  if (!N) return gerepilecopy(av, nfmakesome(n, t, s));
  switch(typ(N))
  {
    case t_INT: X = Xinf = N; break;
    case t_VEC: case t_COL:
      if (lg(N) == 3) { Xinf = gel(N,1); X = gel(N,2); break; }
    default: pari_err_TYPE("nflist", N);
      Xinf = X = NULL;/*LCOV_EXCL_LINE*/
  }
  if (typ(X) != t_INT)
  {
    X = gfloor(X);
    if (typ(X) != t_INT) pari_err_TYPE("nflist", N);
  }
  if (typ(Xinf) != t_INT)
  {
    Xinf = gceil(Xinf);
    if (typ(Xinf) != t_INT) pari_err_TYPE("nflist", N);
  }
  if (signe(Xinf) <= 0)
  {
    if (signe(Xinf) < 0) pari_err_DOMAIN("nflist", "Xinf", "<=", gen_0, Xinf);
    Xinf = gen_1;
  }
  if (signe(X) < 0) pari_err_DOMAIN("nflist", "X", "<=", gen_0, X);
  switch(cmpii(Xinf, X))
  {
    case 1: v = NULL; break;
    case 0: v = nfmakenum(n, t, X, field, s); break;
    default: v = nfmakevecnum(n, t, X, Xinf, field, s);
  }
  if (!v)
  {
    set_avma(av); if (s != -2) return cgetg(1,t_VEC);
    retconst_vec((n>>1) + 1, cgetg(1,t_VEC));
  }
  return gerepilecopy(av, v);
}

/*****************************************************************/
/*                          Polsubcyclo                          */
/*****************************************************************/
/* auxiliary functions assume that trivial impossibilities for s or n
 * are already handled in caller */
static GEN
polsubcycloC2(GEN n, long s)
{
  GEN V = divisorsdisc(n, s), W;
  long l = lg(V), i;
  W = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(W, i) = quadpoly_i(gel(V, i));
  return W;
}
static GEN
polsubcycloC2_i(GEN n, long s)
{
  long l, i;
  GEN V;
  int p, m;
  if (typ(n) == t_VEC)
  {
    fa_is_fundamental_pm(gel(n,1), gel(n,2), s, &p, &m);
    n = gel(n,1);
  }
  else
    is_fundamental_pm(n, s, &p, &m);
  if (!(V = fund_pm(n, p, m))) return NULL;
  l = lg(V);
  for (i = 1; i < l; i++) gel(V, i) = quadpoly_i(gel(V, i));
  return V;
}

static GEN
polsubcycloC3_i(GEN n)
{ GEN P; return checkcondC3(n, &P)? makeC3_i(typ(n) == t_VEC? gel(n,1): n, P)
                                  : NULL; }
/* Cyclic cubic subfields of Q(zeta_n). */
static GEN
polsubcycloC3(GEN n)
{
  long i, l, c;
  GEN D = divisors_factored(n);
  l = lg(D);
  for (i = 2, c = 1; i < l; i++)
  {
    GEN v = polsubcycloC3_i(gel(D,i));
    if (v) gel(D,c++) = v;
  }
  setlg(D, c); return myshallowconcat1(D);
}

static GEN
makeV4pairssimple(GEN D, GEN P, GEN f)
{
  long l = lg(D), n = l-1, i, j, c;
  GEN R = cgetg((n-1) * n / 2 + 1, t_VEC);
  for (i = c = 1; i < n; i++)
  {
    GEN Di = gel(D,i);
    for (j = i + 1; j < l; j++)
    {
      if (f && !equalii(lcmii(Di, gel(D,j)), f)) continue;
      gel(R, c++) = polcompositum0(gel(P,i), gel(P,j), 2);
    }
  }
  setlg(R,c); return R;
}
static GEN
makeV4pairs(GEN D, GEN P, GEN f)
{
  long l = lg(D), n = l-1, i, j, c;
  GEN V = cgetg(l, t_VEC), R = cgetg((n-1) * n / 2 + 1, t_VEC);

  for (i = 1; i < l; i++) gel(V, i) = const_vecsmall(n, 1);
  for (i = c = 1; i < n; i++)
  {
    GEN C = gel(V,i);
    for (j = i + 1; j < l; j++)
      if (C[j])
      { /* Di, Dj fundamental discs */
        GEN d, Di = gel(D,i), Dj = gel(D,j), g = gcdii(Di, Dj);
        long k;
        if (!is_pm1(g)) { Di = diviiexact(Di, g); Dj = diviiexact(Dj, g); }
        d = mulii(Di, Dj); if (f && !equalii(f, mulii(d, g))) continue;
        if (Mod4(d) > 1) d = shifti(d, 2);
        k = vecsearch(D, d, NULL); /* d = coredisc(Di*Dj), j < k */
        C[k] = gel(V, j)[k] = 0;
        gel(R, c++) = polcompositum0(gel(P,i), gel(P,j), 2);
      }
  }
  setlg(R, c); return R;
}
static GEN
polsubcycloV4_i(GEN V, long s, GEN n)
{
  long i, l = lg(V);
  GEN P = cgetg(l, t_VEC);
  if (s <= 0) ZV_sort_inplace(V); /* for vecsearch */
  for (i = 1; i < l; i++) gel(P,i) = quadpoly_i(gel(V,i));
  return (s <= 0)? makeV4pairs(V, P, n): makeV4pairssimple(V, P, n);
}

static GEN
polsubcycloC5(GEN n)
{
  GEN v, D = divisors_factored(n), T = C5bnf();
  long i, c, l = lg(D);
  for (i = 2, c = 1; i < l; i++)
    if ((v = polsubcycloC5_i(gel(D,i), T))) gel(D,c++) = v;
  setlg(D, c); return myshallowconcat1(D);
}

/* ell odd prime */
static GEN
makeCLall(long ell, GEN F)
{
  GEN D = divisors(F);
  long i, l = lg(D);
  for (i = 1; i < l; i++) gel(D,i) = makeCL_f(ell, gel(D,i));
  return shallowconcat1(D);
}

static GEN
polsubcycloC6(GEN n, long s)
{
  GEN v3 = polsubcycloC3(n), v2, R;
  long n3 = lg(v3) - 1, n2, i, j, c;
  if (!n3) return v3;
  v2 = polsubcycloC2(n, s); n2 = lg(v2) - 1;
  if (!n2) return NULL;
  R = cgetg(n2 * n3 + 1, t_VEC);
  for (i = c = 1; i <= n3; i++)
  {
    GEN p3 = gel(v3, i);
    for (j = 1; j <= n2; j++)
      gel(R, c++) = polcompositum0(p3, gel(v2,j), 2);
  }
  return R;
}

static GEN
polsubcycloC6_i(GEN n, long s)
{
  GEN D = divisors_factored(n), R;
  long l = lg(D), i, j, c, L = 2 * (l-1) * omega(n);

  if (typ(n) == t_VEC) n = gel(n,1);
  R = cgetg(L + 1, t_VEC); c = 1;
  for (i = 2; i < l; i++)
  {
    GEN d = gel(D, i), V2 = polsubcycloC2_i(d, s);
    long l2;
    if (!V2) continue;
    l2 = lg(V2);
    if (typ(d) == t_VEC) d = gel(d,1);
    for (j = 1; j < l; j++)
    {
      GEN V3, e = gel(D, j);
      long l3, i3;
      if (!equalii(lcmii(d, typ(e) == t_VEC? gel(e,1): e), n)) continue;
      V3 = polsubcycloC3_i(e); if (!V3) continue;
      l3 = lg(V3);
      for (i3 = 1; i3 < l3; i3++)
      {
        GEN p3 = gel(V3, i3);
        long i2;
        for (i2 = 1; i2 < l2; i2++)
          gel(R, c++) = polcompositum0(p3, gel(V2,i2), 2);
      }
    }
  }
  setlg(R, c); return R;
}

/* fli = 1 for conductor n, else all subfields of Q(zeta_n) */
static GEN
polsubcyclofast_i(GEN n, long ell, long s, long fli)
{
  GEN N, fa = check_arith_pos(n, "polsubcyclofast");

  if (fa && typ(n) != t_VEC) n = mkvec2(factorback(fa), fa);
  /* n either t_INT or [N, factor(N)] */
  if (ell <= 0 && ell != -4)
    pari_err_DOMAIN("polsubcyclofast", "d", "<=", gen_0, stoi(ell));
  /* translate wrt r2 for compatibility with nflist functions */
  if (!s) s = odd(ell)? 0: -1;
  else if (s == 1) s = 0;
  else if (s ==-1)
  {
    if (odd(ell)) return NULL;
    s = labs(ell) >> 1;
  }
  else pari_err_FLAG("polsubcyclo");
  N = fa? gel(n, 1): n;
  if (Mod4(N) == 2)
  {
    if (fli) return NULL;
    N = shifti(N, -1);
    if (fa)
    { /* remove 2^1 */
      GEN P = vecsplice(gel(fa,1), 1), E = vecsplice(gel(fa,2), 1);
      n = mkvec2(N, mkmat2(P, E));
    }
  }
  if (ell == 1)
  {
    if (fli && !equali1(N)) return NULL;
    retmkvec(pol_x(0));
  }
  if (equali1(N)) return NULL;
  if (ell == -4) return polsubcycloV4_i(divisorsdisc(n,s), s, fli? N: NULL);
  if (ell >= 7) return fli? makeCLall(ell,n): makeCL_f(ell,n);
  switch(ell)
  {
    case 2: return fli? polsubcycloC2_i(n, s): polsubcycloC2(n, s);
    case 3: return fli? polsubcycloC3_i(n): polsubcycloC3(n);
    case 4: return fli? polsubcycloC4_i(n, s, fli, NULL): polsubcycloC4(n, s);
    case 5: return fli? polsubcycloC5_i(n, NULL): polsubcycloC5(n);
    case 6: return fli? polsubcycloC6_i(n, s): polsubcycloC6(n, s);
  }
  return NULL; /* LCOV_EXCL_LINE */
}
GEN
polsubcyclofast(GEN n, long ell, long s, long fli)
{
  pari_sp av = avma;
  GEN v = polsubcyclofast_i(n, ell, s, fli);
  if (!v) { set_avma(av); return cgetg(1, t_VEC); }
  return gerepilecopy(av, v);
}
