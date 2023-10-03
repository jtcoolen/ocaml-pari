/* Copyright (C) 2018  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/********************************************************************/
/**       L-functions: values at integers of L-functions           **/
/**             of primitive quadratic characters                  **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

static GEN
RCpol(long k, long t, GEN c)
{
  GEN P = cgetg(k+3, t_POL);
  long l;

  gel(P,k+2) = c;
  for(l = 0; l < k; l++)
  {
    c = diviiexact(mulii(c, muluu(2*k-1 - 2*l, k-l)), mulss(l+1, 2*l-t));
    gel(P,k-l+1) = c;
  }
  P[1] = evalsigne(1) | evalvarn(0); return P;
}
static GEN
vecRCpol(long r, long d)
{
  long k, K = d - 1, t = 2*r - 3;
  GEN v = cgetg(d + 1, t_VEC), c = int2n(2*K);
  for (k = 0; k <= K; k++)
  { /* c = 2^(2K) binomial(n/2,k), an integer */
    gel(v,k+1) = RCpol(k, t, c);
    if (k == K) break;
    c = diviuexact(muliu(c, t - 2*k), 2*k + 2);
  }
  return v;
}
static GEN
euler_sumdiv(GEN q, long v)
{
  GEN u = addui(1, q);
  for (; v > 1; v--) u = addui(1, mulii(q, u));
  return u;
}

/* [p^{k-1},p^{k-3},...,p^{k-2(d-1)-1}] * (s/p), s = 1 or -1 */
static GEN
vpowp(long k, long d, long p, long s)
{
  GEN v = cgetg(d + 1, t_VEC), p2 = sqru(p);
  long j;
  gel(v, d) = powuu(p, k - 2*d + 1);
  if (s == -1 && (p & 3L) == 3) togglesign_safe(&gel(v,d));
  for (j = d-1; j >= 1; j--) gel(v, j) = mulii(p2, gel(v, j+1));
  return v;
}
static GEN
usumdivk_0_all(long k, long d)
{
  GEN v = cgetg(d + 1, t_COL);
  long j;
  constbern(k >> 1);
  for (j = 1; j <= d; j++)
  {
    long n = k + 2 - 2*j;
    gel(v,j) = gdivgs(bernfrac(n), - (n << 1));
  }
  return v;
}
static GEN
usumdivk_fact_all(GEN fa, long k, long d)
{
  GEN res, P, E, pow;
  long i, j, l;
  res = cgetg(d + 1, t_COL);
  P = gel(fa, 1); l = lg(P);
  E = gel(fa, 2); pow = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(pow, i) = vpowp(k, d, P[i], 1);
  for (j = 1; j <= d; j++)
  {
    GEN v = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(v,i) = euler_sumdiv(gmael(pow,i,j), E[i]);
    gel(res, j) = ZV_prod(v);
  }
  return res;
}

/* Hadamard product */
static GEN
RgV_mul(GEN a, GEN b)
{
  long j, l = lg(a);
  GEN v = cgetg(l, t_COL);
  for (j = 1; j < l; j++) gel(v,j) = gmul(gel(a,j), gel(b,j));
  return v;
}
static GEN
RgV_multwist(GEN a, GEN P, long k, long dim, long d, long v2, long N4)
{
  GEN v = cgetg(dim+1, t_COL);
  long j;
  for (j = 1; j <= d; j++)
  {
    GEN z;
    gel(v,j) = z = gmul(gel(a,j), gel(P,j));
    if (j + d <= dim)
    {
      if (N4 == 3) z = negi(z);
      if (v2) z = shifti(z, (k - 2*j + 1)*v2);
      gel(v, j + d) = z;
    }
  }
  return v;
}

/* r = k - 2*j, 0<=j<d, factor s=an+b, 0<=s<lim. Check if n starts at 0 or 1
 * P(D,(an+b)^2), (D-s^2)/N = (D-b^2)/N - 2abn/N - a^2n^2/N and guarantee
 *  N | D-b^2, N | 2ab, and N | a^2 (except N=8, D odd):
 * N=4: a=2, b=0,1\equiv D: D = 0,1 mod 4.
 * N=8: a=4, b=2 if D/4 odd, 0 if D/4 even: D = 0 mod 4 or 1 mod 8
 * N=12: a=6, b=3 if D odd, 0 if D even: D = 0,1 mod 4
 * N=-12: a=6, b=5,1 if D odd, 4,2 if D even: D = 0,1 mod 4
 * N=16: a=8, b=7,1 if D = 1 mod 16, 5,3 if D = 9 mod 16: D = 1 mod 8 */
/* Cost: O( sqrt(D)/a d^3 log(D) ) */
static GEN
sigsum(long k, long d, long a, long b, long D, long N, GEN vs, GEN vP)
{
  pari_sp av;
  GEN S, keep0 = NULL, vPD = RgXV_rescale(vP, stoi(D));
  long D2, n, c1, c2, s, lim = usqrt(labs(D));

  D2 = (D - b*b)/N; c1 = (2*a*b)/N; c2 = (a*a)/N;
  av = avma; S = zerocol(d);
  for (s = b, n = 0; s <= lim; s += a, n++)
  {
    long Ds = c2 ? D2 - n*(c2*n + c1) : D2 - ((n*(n+1)) >> 1);
    GEN v, P = gsubst(vPD, 0, utoi(s*s));
    if (vs)
      v = gel(vs, Ds+1);
    else
      v = Ds? usumdivk_fact_all(factoru(Ds), k, d)
            : usumdivk_0_all(k,d);
    v = RgV_mul(v, P);
    if (!s) keep0 = gclone(v); else S = gadd(S, v);
    if (gc_needed(av, 1)) S = gerepileupto(av, S);
  }
  S = gmul2n(S, 1);
  if (keep0) { S = gadd(S, keep0); gunclone(keep0); }
  return S;
}

static GEN
sigsum4(long k, long d, long D, GEN vs, GEN vP)
{ return sigsum(k, d, 2, odd(D), D, 4, vs, vP); }

/* D != 5 (mod 8) */
static GEN
sigsum8(long k, long d, long D, GEN vs, GEN vP)
{
  if (D&1L) return gmul2n(sigsum(k, d, 2, 1, D, 8, vs, vP), -1);
  return sigsum(k, d, 4, 2*odd(D >> 2), D, 8, vs, vP);
}

/* D = 0 (mod 3) */
static GEN
sigsum12(long k, long d, long D, GEN vs, GEN vP)
{ return sigsum(k, d, 6, 3*odd(D), D, 12, vs, vP); }

/* D = 1 (mod 3) */
static GEN
sigsumm12(long k, long d, long D, GEN vs, GEN vP)
{
  long fl = odd(D);
  GEN res = sigsum(k, d, 6, 4 + fl, D, 12, vs, vP);
  res = gadd(res, sigsum(k, d, 6, 2 - fl, D, 12, vs, vP));
  return gmul2n(res, -1);
}

/* D = 1 (mod 8) */
static GEN
sigsum16(long k, long d, long D, GEN vs, GEN vP)
{
  long fl = (D&15L) == 1;
  GEN res = sigsum(k, d, 8, 5 + 2*fl, D, 16, vs, vP);
  return gadd(res, sigsum(k, d, 8, 3 - 2*fl, D, 16, vs, vP));
}

/* N = 4 (as above), 8 (factor (1+(D/2))), 12 (factor (1+(D/3))),
   16 (only D=1 mod 8). */
static GEN
Dpos(long d, long N, long B)
{
  GEN vD = cgetg(maxss(B, d) + 1, t_VECSMALL);
  long D, step, c;
  switch(N)
  {
    case 4:  D = 5;  step = 1; break;
    case 8:  D = 8;  step = 4; break;
    case 12: D = 12; step = 3; break;
    case 16: D = 17; step = 8; break;
    default: D = 13; step = 3; break; /* -12 */
  }
  for (c = 1; c <= d || D <= B; D += step)
    if (sisfundamental(D)) vD[c++] = D;
  setlg(vD, c); return vD;
}

typedef GEN (*SIGMA_F)(long,long,long,GEN,GEN);
static SIGMA_F
get_S_even(long N)
{
  switch(N) {
    case 4: return sigsum4;
    case 8: return sigsum8;
    case 12:return sigsum12;
    case 16:return sigsum16;
    default:return sigsumm12; /* -12 */
  }
}

static GEN
mfDcoefs(GEN F, GEN vD, long d)
{
  long l = lg(vD), i;
  GEN v = mfcoefs(F, vD[l-1], d), w = cgetg(l, t_COL);
  if (d == 4)
    for (i = 1; i < l; i++) gel(w, i) = gel(v, (vD[i]>>2)+1);
  else
    for (i = 1; i < l; i++) gel(w, i) = gel(v, vD[i]+1);
  return w;
}

static GEN
myinverseimage(GEN M, GEN R, GEN *pden)
{
  GEN c = Q_remove_denom(QM_gauss_i(M, R, 1), pden);/* M*res / den = R */
  if (!c) pari_err_BUG("theta brackets");
  return c;
}

static GEN Lfeq(long D, long k);
static GEN
Hcol(GEN k, long r, GEN vD, long d, long N2)
{
  long i, l = lg(vD);
  GEN v;
  if (r < 5)
  {
    v = mfDcoefs(mfEH(k),vD,d);
    for (i = 1; i < l; i++)
      if (N2 != 1 && vD[i] % N2) gel(v,i) = gmul2n(gel(v,i), 1);
    return v;
  }
  v = cgetg(l, t_COL);
  for (i = 1; i < l; i++)
  {
    pari_sp av = avma;
    GEN c = Lfeq(odd(r)? -vD[i]: vD[i], r); /* fundamental */
    if (N2 != 1 && vD[i] % N2) c = gmul2n(c, 1);
    gel(v, i) = gerepileupto(av, c);
  }
  return v;
}

/***********************************************************/
/*   Modular form method using Half-Integral Weight forms  */
/*                      Case D > 0                         */
/***********************************************************/
static long
dimeven(long r, long N)
{
  switch(N)
  {
    case 4:  return r / 6 + 1;
    case 12: return r / 3 + 1;
    default: return r / 4 + 1;
  }
}
static long
muleven(long N) { return (N == 4)? 1: 2; }

/* L(\chi_D, 1-r) for D > 0 and r > 0 even. */
static GEN
modulareven(long D, long r, long N0)
{
  long B, d, i, l, N = labs(N0);
  GEN V, vs, R, M, C, den, L, vP, vD, k = uutoQ(2*r+1, 2);
  SIGMA_F S = get_S_even(N0);

  d = dimeven(r, N);
  B = muleven(N) * mfsturmNgk(N, k);
  vD = Dpos(d, N0, B);
  vP = vecRCpol(r, d);
  l = lg(vD); B = vD[l-1] / N; V = vecfactoru_i(1, B);
  vs = cgetg(B+2, t_VEC); gel(vs,1) = usumdivk_0_all(r, d);
  for (i = 1; i <= B; i++) gel(vs, i+1) = usumdivk_fact_all(gel(V,i), r, d);
  M = cgetg(l, t_MAT);
  for (i = 1; i < l; i++)
  {
    pari_sp av = avma;
    gel(M,i) = gerepileupto(av, S(r, d, vD[i], vs, vP));
  }
  M = shallowtrans(M);
  if (r == 2*d)
  { /* r = 2 or (r = 4 and N = 4) */
    GEN v = mfDcoefs(mfderiv(mfTheta(NULL), d+1), vD, 1);
    gel(M, d) = gadd(gel(M, d), gdivgu(v, N*(2*d - 1)));
  }
  R = Hcol(k, r, vD, 1, (N == 8 || N0 == 12)? N >> 2: 1);
  /* Cost is O(d^2) * bitsize(result) ~ O(d^3.8) [heuristic] */
  C = myinverseimage(M, R, &den);

  /* Cost: O( sqrt(D)/c d^3 log(D) ), c from findNeven */
  L = RgV_dotproduct(C, S(r, lg(C)-1, D, NULL, vP));
  return den? gdiv(L, den): L;
}

/***********************************************************/
/*   Modular form method using Half-Integral Weight forms  */
/*                      Case D < 0                         */
/***********************************************************/

static long
dimodd(long r, long kro, long N)
{
  switch(N)
  {
    case 1: switch (kro)
    {
      case -1:return (r + 3) >> 2;
      case 0: return (r + 2)/3;
      case 1: return (r + 1) >> 2;
    }
    case 3: return kro? (r + 1) >> 1: ((r << 1) + 2)/3;
    case 5: switch (kro)
    {
      case -1:return (3*r + 2) >> 2;
      case 0: return r;
      case 1: return (3*r - 1) >> 2;
    }
    case 6: return kro == 1 ? (r + 1) >> 1 : r;
    default: return r;
  }
}

static GEN
Dneg(long n, long kro, long d, long N)
{
  GEN vD = cgetg(maxss(n, d) + 1, t_VECSMALL);
  long D, c, step, N2 = odd(N)? N: N>> 1;
  switch(kro)
  {
    case -1: D = -3; step = 8; break;
    case 1:  D = -7; step = 8; break;
    default: D = -8; step = 4; break;
  }
  for (c = 1; D >= -n || c <= d; D -= step)
    if (kross(-D, N2) != -1 && sisfundamental(D)) vD[c++] = -D;
  setlg(vD, c); return vD;
}

static GEN
div4(GEN V)
{
  long l = lg(V), i;
  GEN W = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) W[i] = V[i] >> 2;
  return W;
}

static GEN
usumdivktwist_fact_all(GEN fa, long k, long d)
{
  GEN V, P, E, pow, res = cgetg(d + 1, t_VEC);
  long i, j, l;

  P = gel(fa, 1); l = lg(P);
  E = gel(fa, 2);
  if (l > 1 && P[1] == 2) { l--; P++; E++; } /* odd part */
  pow = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(pow, i) = vpowp(k, d, P[i], -1);
  V = cgetg(l, t_VEC);
  for (j = 1; j <= d; j++)
  {
    for (i = 1; i < l; i++) gel(V,i) = euler_sumdiv(gmael(pow,i,j), E[i]);
    gel(res, j) = ZV_prod(V);
  }
  return res;
}

static long
mulodd(long N, long kro)
{
  if (N == 1 || N == 2) return 1;
  if (kro != 1) return kro? 5: 7;
  if (N == 3) return 4;
  if (N == 5) return 5;
  return 2;
}

/* Cost: O( sqrt(D)/a d^3 log(D) ) */
static GEN
sigsumtwist(long k, long dim, long a, long b, long Da, long N, GEN vs, GEN vP)
{
  GEN vPD, S = zerocol(dim), keep0 = NULL;
  long D2, n, c1, c2, s, lim = usqrt(Da), d;
  pari_sp av;

  if (N > 2 && kross(Da, N == 6 ? 3 : N) == -1) return S;
  d = (dim + 1) >> 1;
  vPD = RgXV_rescale(vP, stoi(Da));
  D2 = (Da - b*b)/N; c1 = (2*a*b)/N; c2 = (a*a)/N;
  av = avma;
  for (s = b, n = 0; s <= lim; s += a, n++)
  {
    long v2, D4, Ds2, Ds = D2 - n*(c2*n + c1); /* (Da - s^2) / N */
    GEN v, P;
    if (!Ds) continue;
    v2 = vals(Ds); Ds2 = Ds >> v2; D4 = Ds2 & 3L; /* (Ds/2^oo) mod 4 */
    if (vs)
      v = gel(vs, Ds+1);
    else
      v = usumdivktwist_fact_all(factoru(Ds2), k, d);
    P = gsubst(vPD, 0, utoi(s*s));
    v = RgV_multwist(v, P, k, dim, d, v2, D4);
    if (!s) keep0 = gclone(v); else S = gadd(S, v);
    if (gc_needed(av, 1)) S = gerepileupto(av, S);
  }
  S = gmul2n(S, 1);
  if (keep0) { S = gadd(S, keep0); gunclone(keep0); }
  return gmul2n(S, -2*(d-1));
}

/* Da = |D|; [sum sigma_r^(1)(Da-s^2), sum sigma_r^(2)(Da-s^2)], N = 1 */
static GEN
sigsumtwist11(long k, long dim, long Da, long N, GEN vs, GEN vP)
{ return sigsumtwist(k, dim, 1, 0, Da, N, vs, vP); }

/* Da = |D| or |D|/4 */
/* [sum sigma_r^(1)((Da-s^2)/N), sum sigma_r^(2)((Da-s^2)/N)] */
/* Case N|Da; N not necessarily prime. */
static GEN
sigsumtwist12p0(long k, long dim, long Da, long N, GEN vs, GEN vP)
{ return sigsumtwist(k, dim, N, 0, Da, N, vs, vP); }

/* [sum sigma_r^(1)((Da-s^2)/p), sum sigma_r^(2)((Da-s^2)/p)] */
/* Case p\nmid Da */
/* p = 3: s = +-1 mod 3;
 * p = 5: s = +-1 mod 5 if Da = 1 mod 5, s = +-2 mod 5 if Da = 2 mod 5;
 * p = 7: s=+-1, +-2, +-3 if Da=1,4,2 mod 7;
 * p = 6: s=+-1, +-2, +-3 if Da=1,4,3 mod 6 */
static GEN
sigsumtwist12pt(long k, long dim, long Da, long N, GEN vs, GEN vP)
{
  long t = Da%N, e = 0;
  GEN res;
  if (t == 1) e = 1;
  else if (t == 4) e = 2;
  else if (t == 2 || t == 3) e = 3;
  res = sigsumtwist(k, dim, N, N-e, Da, N, vs, vP);
  if (N-e != e) res = gadd(res, sigsumtwist(k, dim, N, e, Da, N, vs, vP));
  return res;
}

static GEN
sigsumtwist12_6(long r, long dim, long Da, long N, GEN vs, GEN vP)
{
  if (Da%12 == 6) return sigsumtwist12p0(r, dim, Da, N, vs, vP);
  return sigsumtwist12pt(r, dim, Da, N, vs, vP);
}
static GEN
sigsumtwist12_N(long r, long dim, long Da, long N, GEN vs, GEN vP)
{
  if (Da%N == 0) return sigsumtwist12p0(r, dim, Da, N, vs, vP);
  return sigsumtwist12pt(r, dim, Da, N, vs, vP);
}

typedef GEN (*SIGMA_Fodd)(long,long,long,long,GEN,GEN);
static SIGMA_Fodd
get_S_odd(long N)
{
  if (N == 1) return sigsumtwist11;
  if (N == 6) return sigsumtwist12_6;
  return sigsumtwist12_N;
}

/* L(\chi_D, 1-r) for D < 0 and r > 0 odd. */
static GEN
modularodd(long D, long r, long N0)
{
  long B, d, i, l, dim, kro = kross(D, 2), Da = labs(D), N = labs(N0);
  GEN V, vs, R, M, C, den, L, vP, vD, vD4, k = uutoQ(2*r+1, 2);
  SIGMA_Fodd S = get_S_odd(N);

  dim = dimodd(r, kro, N); d = (dim + 1) >> 1;
  vP = vecRCpol(r, d);
  B = mulodd(N, kro) * mfsturmNgk(4*N, k);
  vD = Dneg(B, kro, dim + 5, N);
  vD4 = kro ? vD : div4(vD);
  l = lg(vD); B = vD4[l-1] / N; V = vecfactoru_i(1, B);
  vs = cgetg(B+2, t_VEC); gel(vs,1) = NULL; /* unused */
  for (i = 1; i <= B; i++) gel(vs,i+1) = usumdivktwist_fact_all(gel(V,i), r, d);
  M = cgetg(l, t_MAT);
  for (i = 1; i < l; i++)
  {
    pari_sp av = avma;
    gel(M,i) = gerepileupto(av, S(r, dim, vD4[i], N, vs, vP));
  }
  M = shallowtrans(M);
  R = Hcol(k, r, vD, kro? 1: 4, odd(N)? N: N >>1);
  /* Cost O(d^2) * bitsize(result) ~ O(d^3.7) [heuristic] */
  C = myinverseimage(M, R, &den);

  if (!kro) Da >>= 2;
  /* Cost: O( sqrt(D)/c d^3 log(D) ), c from findNodd */
  L = RgV_dotproduct(C, S(r, lg(C)-1, Da, N, NULL, vP));
  if (N0 < 0 && (N0 != -6 || Da%3)) den = den? shifti(den,1): gen_2;
  return den? gdiv(L, den): L;
}

/********************************************************/
/*        Using the Full Functional Equation            */
/********************************************************/
/* prod_p (1 - (D/p)p^(-k))
 * Cost O( D/log(D) (k log(kD))^mu ), mu = multiplication exponent */
static GEN
Linv(long D, long k, ulong den)
{
  pari_sp av;
  long s, bit, lim, Da = labs(D), prec;
  double km = k - 1, B = (k-0.5) * log(km*Da/17.079) + 12; /* 17.079 ~ 2Pi e */
  forprime_t iter;
  ulong p;
  GEN P, Q;
  if (den) B += log((double)den);
  bit = maxss((long)(B * k)/(M_LN2 * km), 32) + 32;
  prec = nbits2prec(bit);
  lim = (long)exp( (B-log(km)) / km ); /* ~ D / (2Pi e) */
  u_forprime_init(&iter, 3, lim); av = avma;
  s = kross(D, 2);
  if (!s) P = real_1(prec);
  else
  {
    Q = real2n(-k, nbits2prec(bit - k));
    P = (s == 1)? subir(gen_1, Q): addir(gen_1, Q);
  }
  while ((p = u_forprime_next(&iter)))
  {
    long bitnew;
    GEN Q;
    s = kross(D, p); if (!s) continue;
    bitnew = (long)(bit - k * log2(p));
    Q = divrr(P, rpowuu(p, k, nbits2prec(maxss(64, bitnew))));
    P = s == 1? subrr(P, Q): addrr(P, Q);
    if (gc_needed(av,1)) P = gerepileuptoleaf(av, P);
  }
  return P;
}

static GEN
myround(GEN z, ulong d)
{
  long e;
  if (d) z = mulru(z, d);
  z = grndtoi(z, &e); if (e >= -4) pari_err_BUG("lfunquad");
  return d? Qdiviu(z, d): z;
}

/* D != 1, k > 2; L(\chi_D, 1-k) using func. eq. */
static GEN
Lfeq(long D, long k)
{
  GEN z, res;
  long Da, prec, den = 0;

  if ((D > 0 && odd(k)) || (D < 0 && !odd(k))) return gen_0;
  Da = labs(D);
  if (Da & 3)
  {
    long d = (Da - 1) >> 1, kd = k / d;
    if (odd(kd) && !(k % d) && uisprime(Da)) den = kd * Da;
  }
  else if (Da == 4) den = 2;
  z = Linv(D, k, den); prec = lg(z);
  z = mulrr(z, powrs(divru(Pi2n(1, prec), Da), k));
  if (Da != 4) { z = mulrr(z, sqrtr_abs(utor(Da,prec))); shiftr_inplace(z,-1); }
  res = divrr(mpfactr(k-1, prec), z);
  if (odd(k/2)) togglesign(res);
  return myround(res, den);
}

/* heuristic */
static long
usefeq(long D, long k, double c)
{
  if (k == 2) return 0;
  if (D < 0) { k = 2*k; D = -D; }
  return sqrt(D*c) <= k;
}

static long
findNeven(long D, double *c)
{
  long r = D%3;
  if (!r) { *c = 3; return 12; }
  if ((D&7L) == 1) { *c = 2; return 16; }
  if (!odd(D)) { *c = 2; return 8; }
  if (r == 1) { *c = 1.5; return -12; }
  *c = 1; return 4;
}
static long
findNodd(long D, long k, double *c)
{
  long Dmod8 = D&7L, r;
  if (log(k) > 0.7 * log((double)-D)) { *c = 1; return odd(D)? 2: 1; }
  if (D%7 == 0 && Dmod8 == 5) { *c = 3.5; return 7; }
  if (D%6 == 0) { *c = 3; return 6; }
  if (D%5 == 0) { *c = 2.5; return 5; }
  if (D%3 == 0) { *c = 1.5; return 3; }
  if (Dmod8 == 5)
  {
    r = smodss(D, 7);
    if (r!=1 && r!=2 && r!=4) { *c = 7./6; return -7; }
  }
  if (smodss(D, 3) != 1 && !odd(D)) { *c = 1.5; return -6; }
  r = smodss(D, 5); if (r != 2 && r != 3) { *c = 5./4; return -5; }
  *c = 1; return 2;
}

/* k <= 0 */
static GEN
lfunquadneg_i(long D, long k)
{
  double c;
  long N;

  if (D == 1) return k == 0 ? gneg(ghalf) : gdivgs(bernfrac(1-k), k-1);
  if (!sisfundamental(D)) pari_err_TYPE("lfunquad [D not fundamental]",stoi(D));
  if (k == 0) return D < 0? hclassno(stoi(-D)): gen_0;
  if ((D > 0 && !odd(k)) || (D < 0 && odd(k))) return gen_0;
  if (D == -4) return gmul2n(eulerfrac(-k), -1);
  k = 1 - k;
  N = D < 0? findNodd(D, k, &c): findNeven(D, &c);
  if (usefeq(D, k, c)) return Lfeq(D, k);
  return D < 0? modularodd(D,k,N): modulareven(D,k,N);
}
/* need k <= 0 and D fundamental */
GEN
lfunquadneg(long D, long k)
{ pari_sp av = avma; return gerepileupto(av, lfunquadneg_i(D, k)); }
