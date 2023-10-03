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

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_bern

/********************************************************************/
/**                                                                **/
/**                     BERNOULLI NUMBERS B_2k                     **/
/**                                                                **/
/********************************************************************/

/* D = divisorsu(n). Return a/b = \sum_{p-1 | 2n: p prime} 1/p
 * B_2k + a/b in Z [Clausen-von Staudt] */
static GEN
fracB2k(GEN D)
{
  GEN a = utoipos(5), b = utoipos(6); /* 1/2 + 1/3 */
  long i, l = lg(D);
  for (i = 2; i < l; i++) /* skip 1 */
  {
    ulong p = 2*D[i] + 1; /* a/b += 1/p */
    if (uisprime(p)) { a = addii(muliu(a,p), b); b = muliu(b,p); }
  }
  return mkfrac(a,b);
}
/* precision needed to compute B_k for all k <= N */
long
bernbitprec(long N)
{ /* 1.612086 ~ log(8Pi) / 2 */
  const double log2PI = 1.83787706641;
  double logN = log((double)N);
  double t = (N + 4) * logN - N*(1 + log2PI) + 1.612086;
  return (long)ceil(t / M_LN2) + 10;
}
static long
bernprec(long N) { return nbits2prec(bernbitprec(N)); }
/* \sum_{k > M} k^(-n) <= M^(1-n) / (n-1) < 2^-b */
static long
zetamaxpow(long n)
{
  long M = (long)ceil(n / (2 * M_PI * M_E));
  return M | 1; /* make it odd */
}
/* v * zeta(k) using r precomputed odd powers */
static GEN
bern_zeta(GEN v, long k, GEN pow, long r, long p)
{
  GEN z, s = gel(pow, r);
  long j;
  for (j = r - 2; j >= 3; j -= 2) s = addii(s, gel(pow,j));
  z = s = itor(s, nbits2prec(p));
  shiftr_inplace(s, -p); /* zeta(k)(1 - 2^(-k)) - 1*/
  s = addrs(s, 1); shiftr_inplace(s, -k);
  /* divide by 1 - 2^(-k): s + s/2^k + s/2^(2k) + ... */
  for (; k < p; k <<= 1) s = addrr(s, shiftr(s, -k));
  return addrr(v, mulrr(v, addrr(z, s)));
}
/* z * j^2 */
static GEN
muliu2(GEN z, ulong j)
{ return (j | HIGHMASK)? mulii(z, sqru(j)): muliu(z, j*j); }
/* 1 <= m <= n, set y[1] = B_{2m}, ... y[n-m+1] = B_{2n} in Q */
static void
bernset(GEN *y, long m, long n)
{
  long i, j, k, p, prec, r, N = n << 1; /* up to B_N */
  GEN u, b, v, t;
  p = bernbitprec(N); prec = nbits2prec(p);
  u = sqrr(Pi2n(1, prec)); /* (2Pi)^2 */
  v = divrr(mpfactr(N, prec), powru(u, n)); shiftr_inplace(v,1);
  r = zetamaxpow(N);
  t = cgetg(r+1, t_VEC); b = int2n(p); /* fixed point */
  for (j = 3; j <= r; j += 2)
  {
    GEN z = cgeti(prec);
    pari_sp av2 = avma;
    affii(divii(b, powuu(j, N)), z);
    gel(t,j) = z; set_avma(av2);
  }
  y += n - m;
  for (i = n, k = N;; i--)
  { /* set B_n, k = 2i */
    pari_sp av2 = avma;
    GEN z = fracB2k(divisorsu(i)), B = bern_zeta(v, k, t, r, p);
    long j;
    /* B = v * zeta(k), v = 2*k! / (2Pi)^k */
    if (!odd(i)) setsigne(B, -1); /* B ~ B_n */
    B = roundr(addrr(B, fractor(z,LOWDEFAULTPREC))); /* B - z = B_n */
    *y-- = gclone(gsub(B, z));
    if (i == m) break;
    affrr(divrunextu(mulrr(v,u), k-1), v);
    for (j = r; j >= 3; j -= 2) affii(muliu2(gel(t,j), j), gel(t,j));
    set_avma(av2); k -= 2;
    if (((N - k) & 0x7f) == 0x7e)
    { /* reduce precision if possible */
      long p2 = p, prec2 = prec;
      p = bernbitprec(k); prec = nbits2prec(p); if (prec2 == prec) continue;
      setprec(v, prec); r = zetamaxpow(k);
      for (j = 3; j <= r; j += 2) affii(shifti(gel(t,j), p - p2), gel(t,j));
      set_avma(av2);
    }
  }
}
/* need B[2..2*nb] as t_INT or t_FRAC */
void
constbern(long nb)
{
  const pari_sp av = avma;
  long i, l;
  GEN B;
  pari_timer T;

  l = bernzone? lg(bernzone): 0;
  if (l > nb) return;

  nb = maxss(nb, l + 127);
  B = cgetg_block(nb+1, t_VEC);
  if (bernzone)
  { for (i = 1; i < l; i++) gel(B,i) = gel(bernzone,i); }
  else
  {
    gel(B,1) = gclone(mkfracss(1,6));
    gel(B,2) = gclone(mkfracss(-1,30));
    gel(B,3) = gclone(mkfracss(1,42));
    gel(B,4) = gclone(mkfracss(-1,30));
    gel(B,5) = gclone(mkfracss(5,66));
    gel(B,6) = gclone(mkfracss(-691,2730));
    gel(B,7) = gclone(mkfracss(7,6));
    gel(B,8) = gclone(mkfracss(-3617,510));
    gel(B,9) = gclone(mkfracss(43867,798));
    gel(B,10)= gclone(mkfracss(-174611,330));
    gel(B,11)= gclone(mkfracss(854513,138));
    gel(B,12)= gclone(mkfracss(-236364091,2730));
    gel(B,13)= gclone(mkfracss(8553103,6)); /* B_26 */
    l = 14;
  }
  set_avma(av);
  if (DEBUGLEVEL) {
    err_printf("caching Bernoulli numbers 2*%ld to 2*%ld\n", l, nb);
    timer_start(&T);
  }
  bernset((GEN*)B + l, l, nb);
  if (DEBUGLEVEL) timer_printf(&T, "Bernoulli");
  swap(B, bernzone); guncloneNULL(B);
  set_avma(av);
#if 0
  if (nb > 200000)
#endif
  {
    const ulong p = 4294967291UL;
    long n = 2 * nb + 2;
    GEN t = const_vecsmall(n+1, 1);
    t[1] = evalvarn(0); t[2] = 0;
    t = Flx_shift(Flx_invLaplace(t, p), -1); /* t = (exp(x)-1)/x */
    t = Flx_Laplace(Flxn_inv(t, n, p), p);
    for (i = 1; i <= nb; i++)
      if (Rg_to_Fl(bernfrac(2*i), p) != uel(t,2*i+2))
      {
        gunclone(bernzone); bernzone = NULL;
        pari_err_BUG(stack_sprintf("B_{2*%ld}", i));
      }
    set_avma(av);
  }
}
/* Obsolete, kept for backward compatibility */
void
mpbern(long n, long prec) { (void)prec; constbern(n); }

/* assume n even > 0, if iz != NULL, assume iz = 1/zeta(n) */
static GEN
bernreal_using_zeta(long n, long prec)
{
  GEN pi2 = Pi2n(1, prec+EXTRAPREC64);
  GEN iz = inv_szeta_euler(n, prec);
  GEN z = divrr(mpfactr(n, prec), mulrr(powru(pi2, n), iz));
  shiftr_inplace(z, 1); /* 2 * n! * zeta(n) / (2Pi)^n */
  if ((n & 3) == 0) setsigne(z, -1);
  return z;
}
/* assume n even > 0, B = NULL or good approximation to B_n */
static GEN
bernfrac_i(long n, GEN B)
{
  GEN z = fracB2k(divisorsu(n >> 1));
  if (!B) B = bernreal_using_zeta(n, bernprec(n));
  B = roundr( gadd(B, fractor(z,LOWDEFAULTPREC)) );
  return gsub(B, z);
}
GEN
bernfrac(long n)
{
  pari_sp av;
  long k;
  if (n <= 1)
  {
    if (n < 0) pari_err_DOMAIN("bernfrac", "index", "<", gen_0, stoi(n));
    return n? mkfrac(gen_m1,gen_2): gen_1;
  }
  if (odd(n)) return gen_0;
  k = n >> 1;
  if (!bernzone) constbern(0);
  if (bernzone && k < lg(bernzone)) return gel(bernzone, k);
  av = avma;
  return gerepileupto(av, bernfrac_i(n, NULL));
}
GEN
bernvec(long n)
{
  long i, l;
  GEN y;
  if (n < 0) return cgetg(1, t_VEC);
  constbern(n);
  l = n+2; y = cgetg(l, t_VEC); gel(y,1) = gen_1;
  for (i = 2; i < l; i++) gel(y,i) = gel(bernzone,i-1);
  return y;
}

/* x := pol_x(v); B_k(x) = \sum_{i=0}^k binomial(k, i) B_i x^{k-i} */
static GEN
bernpol_i(long k, long v)
{
  GEN B, C;
  long i;
  if (v < 0) v = 0;
  constbern(k >> 1); /* cache B_2, ..., B_2[k/2] */
  C = vecbinomial(k);
  B = cgetg(k + 3, t_POL);
  for (i = 0; i <= k; ++i) gel(B, k-i+2) = gmul(gel(C,i+1), bernfrac(i));
  B[1] = evalsigne(1) | evalvarn(v);
  return B;
}
GEN
bernpol(long k, long v)
{
  pari_sp av = avma;
  if (k < 0) pari_err_DOMAIN("bernpol", "index", "<", gen_0, stoi(k));
  return gerepileupto(av, bernpol_i(k, v));
}
GEN
bernpol_eval(long k, GEN x)
{
  pari_sp av = avma;
  GEN B;
  if (!x) return bernpol(k, 0);
  if (gequalX(x)) return bernpol(k, varn(x));
  if (k < 0) pari_err_DOMAIN("bernpol", "index", "<", gen_0, stoi(k));
  B = bernpol_i(k, fetch_var_higher());
  return gerepileupto(av, poleval(B, x));
}

/* x := pol_x(v); return 1^e + ... + x^e = x^e + (B_{e+1}(x) - B_{e+1})/(e+1) */
static GEN
faulhaber(long e, long v)
{
  GEN B;
  if (e == 0) return pol_x(v);
  B = RgX_integ(bernpol_i(e, v)); /* (B_{e+1}(x) - B_{e+1}) / (e+1) */
  gel(B,e+2) = gaddgs(gel(B,e+2), 1); /* add x^e, in place */
  return B;
}
/* sum_v T(v), T a polynomial expression in v */
GEN
sumformal(GEN T, long v)
{
  pari_sp av = avma, av2;
  long i, t, d;
  GEN R;

  T = simplify_shallow(T);
  t = typ(T);
  if (is_scalar_t(t))
    return gerepileupto(av, monomialcopy(T, 1, v < 0? 0: v));
  if (t != t_POL) pari_err_TYPE("sumformal [not a t_POL]", T);
  if (v < 0) v = varn(T);
  av2 = avma;
  R = gen_0;
  d = poldegree(T,v);
  for (i = d; i >= 0; i--)
  {
    GEN c = polcoef_i(T, i, v);
    if (gequal0(c)) continue;
    R = gadd(R, gmul(c, faulhaber(i, v)));
    if (gc_needed(av2,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"sumformal, i = %ld/%ld", i,d);
      R = gerepileupto(av2, R);
    }
  }
  return gerepileupto(av, R);
}

/* 1/zeta(n) using Euler product. Assume n > 0. */
GEN
inv_szeta_euler(long n, long prec)
{
  long bit = prec2nbits(prec);
  GEN z, res;
  pari_sp av, av2;
  double A, D, lba;
  ulong p, lim;
  forprime_t S;

  if (n > bit) return real_1(prec);

  lba = prec2nbits_mul(prec, M_LN2);
  D = exp((lba - log((double)(n-1))) / (n-1));
  lim = 1 + (ulong)ceil(D);
  if (lim < 3) return subir(gen_1,real2n(-n,prec));
  res = cgetr(prec); av = avma; incrprec(prec);

  (void)u_forprime_init(&S, 3, lim);
  av2 = avma; A = n / M_LN2; z = subir(gen_1, real2n(-n, prec));
  while ((p = u_forprime_next(&S)))
  {
    long l = bit - (long)floor(A * log((double)p));
    GEN h;

    if (l < BITS_IN_LONG) l = BITS_IN_LONG;
    l = minss(prec, nbits2prec(l));
    h = divrr(z, rpowuu(p, (ulong)n, l));
    z = subrr(z, h);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"inv_szeta_euler, p = %lu/%lu", p,lim);
      z = gerepileuptoleaf(av2, z);
    }
  }
  affrr(z, res); set_avma(av); return res;
}

/* Return B_n */
GEN
bernreal(long n, long prec)
{
  pari_sp av;
  GEN B;
  long p, k;
  if (n < 0) pari_err_DOMAIN("bernreal", "index", "<", gen_0, stoi(n));
  if (n == 0) return real_1(prec);
  if (n == 1) return real_m2n(-1,prec); /*-1/2*/
  if (odd(n)) return real_0(prec);

  k = n >> 1;
  if (!bernzone) constbern(0);
  if (k < lg(bernzone)) return fractor(gel(bernzone,k), prec);
  p = bernprec(n); av = avma;
  B = bernreal_using_zeta(n, minss(p, prec));
  if (p < prec) B = fractor(bernfrac_i(n, B), prec);
  return gerepileuptoleaf(av, B);
}

GEN
eulerpol(long k, long v)
{
  pari_sp av = avma;
  GEN B, E;
  if (k < 0) pari_err_DOMAIN("eulerpol", "index", "<", gen_0, stoi(k));
  k++; B = bernpol_i(k, v);
  E = RgX_Rg_mul(RgX_sub(B, RgX_rescale(B, gen_2)), uutoQ(2,k));
  return gerepileupto(av, E);
}

/*******************************************************************/
/**                      HARMONIC NUMBERS                         **/
/*******************************************************************/
/* 1/a + ... + 1/(b-1); a < b <= 2^(BIL-1) */
static GEN
hrec(ulong a, ulong b)
{
  ulong m;
  switch(b - a)
  {
    case 1: retmkfrac(gen_1, utoipos(a));
    case 2: if (a < 65536) retmkfrac(utoipos(2*a + 1), utoipos(a * a + a));
      retmkfrac(utoipos(2*a + 1), muluu(a, a+1));
  }
  m = (a + b) >> 1;
  return gadd(hrec(a, m), hrec(m, b));
}
/* exact Harmonic number H_n, n < 2^(BIL-1).
 * Could use H_n = sum_k 2^(-k) H^odd_{n \ 2^k} */
GEN
harmonic(ulong n)
{
  pari_sp av = avma;
  return n? gerepileupto(av, hrec(1, n+1)): gen_0;
}

/* 1/a^k + ... + 1/(b-1)^k; a < b */
static GEN
hreck(ulong a, ulong b, ulong k)
{
  ulong m;
  switch(b - a)
  {
    GEN x, y;
    case 1: retmkfrac(gen_1, powuu(a, k));
    case 2:
      x = powuu(a, k); y = powuu(a + 1, k);
      retmkfrac(addii(x, y), mulii(x, y));
  }
  m = (a + b) >> 1;
  return gadd(hreck(a, m, k), hreck(m, b, k));
}
GEN
harmonic0(ulong n, GEN k)
{
  pari_sp av = avma;
  ulong r;
  if (!n) return gen_0;
  if (n & HIGHBIT) pari_err_OVERFLOW("harmonic");
  if (!k) return harmonic(n);
  if (typ(k) != t_INT) pari_err_TYPE("harmonic", k);
  if (signe(k) < 0)
  {
    GEN H = poleval(faulhaber(-itos(k), 0), utoipos(n));
    return gerepileuptoint(av, H);
  }
  r = itou(k);
  if (!r) return utoipos(n);
  if (r == 1) return harmonic(n);
  return gerepileupto(av, hreck(1, n+1, r));
}


/**************************************************************/
/*                      Euler numbers                         */
/**************************************************************/

/* precision needed to compute E_k for all k <= N */
static long
eulerbitprec(long N)
{ /* 1.1605 ~ log(32/Pi) / 2 */
  const double logPIS2 = 0.4515827;
  double t = (N + 1) * log((double)N) - N*(1 + logPIS2) + 1.1605;
  return (long)ceil(t / M_LN2) + 10;
}
static long
eulerprec(long N) { return nbits2prec(eulerbitprec(N)); }

/* sum_{k > M, k odd} (-1)^((k-1)/2)k^(-n) < M^(-n) < 2^-b */
static long
lfun4maxpow(long n)
{
  long M = (long)ceil(2 * n / (M_E * M_PI));
  return M | 1; /* make it odd */
}

/* lfun4(k) using r precomputed odd powers */
static GEN
euler_lfun4(GEN v, GEN pow, long r, long p)
{
  GEN s = ((r & 3L) == 1)? gel(pow, r): negi(gel(pow, r));
  long j;
  for (j = r - 2; j >= 3; j -= 2)
    s = ((j & 3L) == 1)? addii(s, gel(pow,j)): subii(s, gel(pow,j));
  s = mulri(v, s); shiftr_inplace(s, -p);
  return addrr(v, s);
}

/* 1 <= m <= n, set y[1] = E_{2m}, ... y[n-m+1] = E_{2n} in Z */
static void
eulerset(GEN *y, long m, long n)
{
  long i, j, k, p, prec, r, N = n << 1, N1 = N + 1; /* up to E_N */
  GEN b, u, v, t;
  p = eulerbitprec(N); prec = nbits2prec(p);
  u = sqrr(Pi2n(-1, prec)); /* (Pi/2)^2 */
  v = divrr(mpfactr(N, prec), mulrr(powru(u, n), Pi2n(-2,prec)));
  r = lfun4maxpow(N1);
  t = cgetg(r+1, t_VEC); b = int2n(p); /* fixed point */
  for (j = 3; j <= r; j += 2)
  {
    GEN z = cgeti(prec);
    pari_sp av2 = avma;
    affii(divii(b, powuu(j, N+1)), z);
    gel(t,j) = z; set_avma(av2);
  }
  y += n - m;
  for (i = n, k = N1;; i--)
  { /* set E_n, k = 2i + 1 */
    pari_sp av2 = avma;
    GEN E = euler_lfun4(v, t, r, p);
    long j;
    /* E = v * lfun4(k), v = (4/Pi)*k! / (Pi/2)^k */
    E = roundr(E); if (odd(i)) setsigne(E, -1); /* E ~ E_n */
    *y-- = gclone(E);
    if (i == m) break;
    affrr(divrunextu(mulrr(v,u), k-2), v);
    for (j = r; j >= 3; j -= 2) affii(muliu2(gel(t,j), j), gel(t,j));
    set_avma(av2); k -= 2;
    if (((N1 - k) & 0x7f) == 0x7e)
    { /* reduce precision if possible */
      long p2 = p, prec2 = prec;
      p = eulerbitprec(k); prec = nbits2prec(p); if (prec2 == prec) continue;
      setprec(v, prec); r = lfun4maxpow(k);
      for (j = 3; j <= r; j += 2) affii(shifti(gel(t,j), p - p2), gel(t,j));
      set_avma(av2);
    }
  }
}

/* need E[2..2*nb] as t_INT */
static void
constreuler(long nb)
{
  const pari_sp av = avma;
  long i, l;
  GEN E;
  pari_timer T;

  l = eulerzone? lg(eulerzone): 0;
  if (l > nb) return;

  nb = maxss(nb, l + 127);
  E = cgetg_block(nb+1, t_VEC);
  if (eulerzone)
  { for (i = 1; i < l; i++) gel(E,i) = gel(eulerzone,i); }
  else
  {
    gel(E,1) = gclone(stoi(-1));
    gel(E,2) = gclone(stoi(5));
    gel(E,3) = gclone(stoi(-61));
    gel(E,4) = gclone(stoi(1385));
    gel(E,5) = gclone(stoi(-50521));
    gel(E,6) = gclone(stoi(2702765));
    gel(E,7) = gclone(stoi(-199360981));
    l = 8;
  }
  set_avma(av);
  if (DEBUGLEVEL) {
    err_printf("caching Euler numbers 2*%ld to 2*%ld\n", l, nb);
    timer_start(&T);
  }
  eulerset((GEN*)E + l, l, nb);
  if (DEBUGLEVEL) timer_printf(&T, "Euler");
  swap(E, eulerzone); guncloneNULL(E);
  set_avma(av);
}

/* 1/lfun(-4,n) using Euler product. Assume n > 0. */
static GEN
inv_lfun4(long n, long prec)
{
  long bit = prec2nbits(prec);
  GEN z, res;
  pari_sp av, av2;
  double A;
  ulong p, lim;
  forprime_t S;

  if (n > bit) return real_1(prec);

  lim = 1 + (ulong)ceil(exp2((double)bit / n));
  res = cgetr(prec); av = avma; incrprec(prec);

  (void)u_forprime_init(&S, 3, lim);
  av2 = avma; A = n / M_LN2; z = real_1(prec);
  while ((p = u_forprime_next(&S)))
  {
    long l = bit - (long)floor(A * log((double)p));
    GEN h;

    if (l < BITS_IN_LONG) l = BITS_IN_LONG;
    l = minss(prec, nbits2prec(l));
    h = rpowuu(p, (ulong)n, l); if ((p & 3UL) == 1) setsigne(h, -1);
    z = addrr(z, divrr(z, h)); /* z *= 1 - chi_{-4}(p) / p^n */
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"inv_lfun4, p = %lu/%lu", p,lim);
      z = gerepileuptoleaf(av2, z);
    }
  }
  affrr(z, res); set_avma(av); return res;
}
/* assume n even > 0, E_n = (-1)^(n/2) (4/Pi) n! lfun4(n+1) / (Pi/2)^n */
static GEN
eulerreal_using_lfun4(long n, long prec)
{
  GEN pisur2 = Pi2n(-1, prec+EXTRAPREC64);
  GEN iz = inv_lfun4(n+1, prec);
  GEN z = divrr(mpfactr(n, prec), mulrr(powru(pisur2, n+1), iz));
  if ((n & 3L) == 2) setsigne(z, -1);
  shiftr_inplace(z, 1); return z;
}
/* Euler numbers: 1, 0, -1, 0, 5, 0, -61,... */
GEN
eulerfrac(long n)
{
  pari_sp av;
  long k;
  GEN E;
  if (n <= 0)
  {
    if (n < 0) pari_err_DOMAIN("eulerfrac", "index", "<", gen_0, stoi(n));
    return gen_1;
  }
  if (odd(n)) return gen_0;
  k = n >> 1;
  if (!eulerzone) constreuler(0);
  if (eulerzone && k < lg(eulerzone)) return gel(eulerzone, k);
  av = avma; E = eulerreal_using_lfun4(n, eulerprec(n));
  return gerepileuptoleaf(av, roundr(E));
}
GEN
eulervec(long n)
{
  long i, l;
  GEN y;
  if (n < 0) return cgetg(1, t_VEC);
  constreuler(n);
  l = n+2; y = cgetg(l, t_VEC); gel(y,1) = gen_1;
  for (i = 2; i < l; i++) gel(y,i) = gel(eulerzone,i-1);
  return y;
}

/* Return E_n */
GEN
eulerreal(long n, long prec)
{
  pari_sp av = avma;
  GEN B;
  long p, k;
  if (n < 0) pari_err_DOMAIN("eulerreal", "index", "<", gen_0, stoi(n));
  if (n == 0) return real_1(prec);
  if (odd(n)) return real_0(prec);

  k = n >> 1;
  if (!eulerzone) constreuler(0);
  if (k < lg(eulerzone)) return itor(gel(eulerzone,k), prec);
  p = eulerprec(n);
  B = eulerreal_using_lfun4(n, minss(p, prec));
  if (p < prec) B = itor(roundr(B), prec);
  return gerepileuptoleaf(av, B);
}
