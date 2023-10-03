/* Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/*******************************************************************/
/*                                                                 */
/*                ROOTS OF COMPLEX POLYNOMIALS                     */
/*  (original code contributed by Xavier Gourdon, INRIA RR 1852)   */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_polroots

static const double pariINFINITY = 1./0.;

static long
isvalidcoeff(GEN x)
{
  switch (typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC: return 1;
    case t_COMPLEX: return isvalidcoeff(gel(x,1)) && isvalidcoeff(gel(x,2));
  }
  return 0;
}

static void
checkvalidpol(GEN p, const char *f)
{
  long i,n = lg(p);
  for (i=2; i<n; i++)
    if (!isvalidcoeff(gel(p,i))) pari_err_TYPE(f, gel(p,i));
}

/********************************************************************/
/**                                                                **/
/**                   FAST ARITHMETIC over Z[i]                    **/
/**                                                                **/
/********************************************************************/

static GEN
ZX_to_ZiX(GEN Pr, GEN Pi)
{
  long i, lr = lg(Pr), li = lg(Pi), l = maxss(lr, li), m = minss(lr, li);
  GEN P = cgetg(l, t_POL);
  P[1] = Pr[1];
  for(i = 2; i < m; i++)
    gel(P,i) = signe(gel(Pi,i)) ? mkcomplex(gel(Pr,i), gel(Pi,i))
                                : gel(Pr,i);
  for(     ; i < lr; i++)
    gel(P,i) = gel(Pr, i);
  for(     ; i < li; i++)
    gel(P,i) = mkcomplex(gen_0, gel(Pi, i));
  return normalizepol_lg(P, l);
}

static GEN
ZiX_sqr(GEN P)
{
  pari_sp av = avma;
  GEN Pr2, Pi2, Qr, Qi;
  GEN Pr = real_i(P), Pi = imag_i(P);
  if (signe(Pi)==0) return gerepileupto(av, ZX_sqr(Pr));
  if (signe(Pr)==0) return gerepileupto(av, ZX_neg(ZX_sqr(Pi)));
  Pr2 = ZX_sqr(Pr); Pi2 = ZX_sqr(Pi);
  Qr = ZX_sub(Pr2, Pi2);
  if (degpol(Pr)==degpol(Pi))
    Qi = ZX_sub(ZX_sqr(ZX_add(Pr, Pi)), ZX_add(Pr2, Pi2));
  else
    Qi = ZX_shifti(ZX_mul(Pr, Pi), 1);
  return gerepilecopy(av, ZX_to_ZiX(Qr, Qi));
}

static GEN
graeffe(GEN p)
{
  pari_sp av = avma;
  GEN p0, p1, s0, s1;
  long n = degpol(p);

  if (!n) return RgX_copy(p);
  RgX_even_odd(p, &p0, &p1);
  /* p = p0(x^2) + x p1(x^2) */
  s0 = ZiX_sqr(p0);
  s1 = ZiX_sqr(p1);
  return gerepileupto(av, RgX_sub(s0, RgX_shift_shallow(s1,1)));
}

GEN
ZX_graeffe(GEN p)
{
  pari_sp av = avma;
  GEN p0, p1, s0, s1;
  long n = degpol(p);

  if (!n) return ZX_copy(p);
  RgX_even_odd(p, &p0, &p1);
  /* p = p0(x^2) + x p1(x^2) */
  s0 = ZX_sqr(p0);
  s1 = ZX_sqr(p1);
  return gerepileupto(av, ZX_sub(s0, RgX_shift_shallow(s1,1)));
}
GEN
polgraeffe(GEN p)
{
  pari_sp av = avma;
  GEN p0, p1, s0, s1;
  long n = degpol(p);

  if (typ(p) != t_POL) pari_err_TYPE("polgraeffe",p);
  n = degpol(p);
  if (!n) return gcopy(p);
  RgX_even_odd(p, &p0, &p1);
  /* p = p0(x^2) + x p1(x^2) */
  s0 = RgX_sqr(p0);
  s1 = RgX_sqr(p1);
  return gerepileupto(av, RgX_sub(s0, RgX_shift_shallow(s1,1)));
}

/********************************************************************/
/**                                                                **/
/**                       MODULUS OF ROOTS                         **/
/**                                                                **/
/********************************************************************/

/* Quick approximation to log2(|x|); first define y s.t. |y-x| < 2^-32 then
 * return y rounded to 2 ulp. In particular, if result < 2^21, absolute error
 * is bounded by 2^-31. If result > 2^21, it is correct to 2 ulp */
static double
mydbllog2i(GEN x)
{
#ifdef LONG_IS_64BIT
  const double W = 1/(4294967296. * 4294967296.); /* 2^-64 */
#else
  const double W = 1/4294967296.; /*2^-32*/
#endif
  GEN m;
  long lx = lgefint(x);
  double l;
  if (lx == 2) return -pariINFINITY;
  m = int_MSW(x);
  l = (double)(ulong)*m;
  if (lx == 3) return log2(l);
  l += ((double)(ulong)*int_precW(m)) * W;
  /* at least m = min(53,BIL) bits are correct in the mantissa, thus log2
   * is correct with error < log(1 + 2^-m) ~ 2^-m. Adding the correct
   * exponent BIL(lx-3) causes 1ulp further round-off error */
  return log2(l) + (double)(BITS_IN_LONG*(lx-3));
}

/* return log(|x|) or -pariINFINITY */
static double
mydbllogr(GEN x) {
  if (!signe(x)) return -pariINFINITY;
  return M_LN2*dbllog2r(x);
}

/* return log2(|x|) or -pariINFINITY */
static double
mydbllog2r(GEN x) {
  if (!signe(x)) return -pariINFINITY;
  return dbllog2r(x);
}
double
dbllog2(GEN z)
{
  double x, y;
  switch(typ(z))
  {
    case t_INT: return mydbllog2i(z);
    case t_FRAC: return mydbllog2i(gel(z,1))-mydbllog2i(gel(z,2));
    case t_REAL: return mydbllog2r(z);
    default: /*t_COMPLEX*/
      x = dbllog2(gel(z,1));
      y = dbllog2(gel(z,2));
      if (x == -pariINFINITY) return y;
      if (y == -pariINFINITY) return x;
      if (fabs(x-y) > 10) return maxdd(x,y);
      return x + 0.5*log2(1 + exp2(2*(y-x)));
  }
}
static GEN /* beware overflow */
dblexp(double x) { return fabs(x) < 100.? dbltor(exp(x)): mpexp(dbltor(x)); }

/* find s such that  A_h <= 2^s <= 2 A_i  for one h and all i < n = deg(p),
 * with  A_i := (binom(n,i) lc(p) / p_i) ^ 1/(n-i), and  p = sum p_i X^i */
static long
findpower(GEN p)
{
  double x, L, mins = pariINFINITY;
  long n = degpol(p),i;

  L = dbllog2(gel(p,n+2)); /* log2(lc * binom(n,i)) */
  for (i=n-1; i>=0; i--)
  {
    L += log2((double)(i+1) / (double)(n-i));
    x = dbllog2(gel(p,i+2));
    if (x != -pariINFINITY)
    {
      double s = (L - x) / (double)(n-i);
      if (s < mins) mins = s;
    }
  }
  i = (long)ceil(mins);
  if (i - mins > 1 - 1e-12) i--;
  return i;
}

/* returns the exponent for logmodulus(), from the Newton diagram */
static long
newton_polygon(GEN p, long k)
{
  pari_sp av = avma;
  long n = degpol(p), i, j, h, l, *vertex = (long*)new_chunk(n+1);
  double *L = (double*)stack_malloc_align((n+1)*sizeof(double), sizeof(double));

  /* vertex[i] = 1 if i a vertex of convex hull, 0 otherwise */
  for (i=0; i<=n; i++) { L[i] = dbllog2(gel(p,2+i)); vertex[i] = 0; }
  vertex[0] = 1; /* sentinel */
  for (i=0; i < n; i=h)
  {
    double slope;
    h = i+1;
    while (L[i] == -pariINFINITY) { vertex[h] = 1; i = h; h = i+1; }
    slope = L[h] - L[i];
    for (j = i+2; j<=n; j++) if (L[j] != -pariINFINITY)
    {
      double pij = (L[j] - L[i])/(double)(j - i);
      if (slope < pij) { slope = pij; h = j; }
    }
    vertex[h] = 1;
  }
  h = k;   while (!vertex[h]) h++;
  l = k-1; while (!vertex[l]) l--;
  set_avma(av);
  return (long)floor((L[h]-L[l])/(double)(h-l) + 0.5);
}

/* change z into z*2^e, where z is real or complex of real */
static void
myshiftrc(GEN z, long e)
{
  if (typ(z)==t_COMPLEX)
  {
    if (signe(gel(z,1))) shiftr_inplace(gel(z,1), e);
    if (signe(gel(z,2))) shiftr_inplace(gel(z,2), e);
  }
  else
    if (signe(z)) shiftr_inplace(z, e);
}

/* return z*2^e, where z is integer or complex of integer (destroy z) */
static GEN
myshiftic(GEN z, long e)
{
  if (typ(z)==t_COMPLEX)
  {
    gel(z,1) = signe(gel(z,1))? mpshift(gel(z,1),e): gen_0;
    gel(z,2) = mpshift(gel(z,2),e);
    return z;
  }
  return signe(z)? mpshift(z,e): gen_0;
}

static GEN
RgX_gtofp_bit(GEN q, long bit) { return RgX_gtofp(q, nbits2prec(bit)); }

static GEN
mygprecrc(GEN x, long prec, long e)
{
  GEN y;
  switch(typ(x))
  {
    case t_REAL:
      if (!signe(x)) return real_0_bit(e);
      return realprec(x) == prec? x: rtor(x, prec);
    case t_COMPLEX:
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = mygprecrc(gel(x,1),prec,e);
      gel(y,2) = mygprecrc(gel(x,2),prec,e);
      return y;
    default: return x;
  }
}

/* gprec behaves badly with the zero for polynomials.
The second parameter in mygprec is the precision in base 2 */
static GEN
mygprec(GEN x, long bit)
{
  long lx, i, e, prec;
  GEN y;

  if (bit < 0) bit = 0; /* should rarely happen */
  e = gexpo(x) - bit;
  prec = nbits2prec(bit);
  switch(typ(x))
  {
    case t_POL:
      y = cgetg_copy(x, &lx); y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = mygprecrc(gel(x,i),prec,e);
      break;

    default: y = mygprecrc(x,prec,e);
  }
  return y;
}

/* normalize a polynomial p, that is change it with coefficients in Z[i],
after making product by 2^shift */
static GEN
pol_to_gaussint(GEN p, long shift)
{
  long i, l = lg(p);
  GEN q = cgetg(l, t_POL); q[1] = p[1];
  for (i=2; i<l; i++) gel(q,i) = gtrunc2n(gel(p,i), shift);
  return q;
}

/* returns a polynomial q in Z[i][x] keeping bit bits of p */
static GEN
eval_rel_pol(GEN p, long bit)
{
  long i;
  for (i = 2; i < lg(p); i++)
    if (gequal0(gel(p,i))) gel(p,i) = gen_0; /* bad behavior of gexpo */
  return pol_to_gaussint(p, bit-gexpo(p)+1);
}

/* returns p(R*x)/R^n (in R or R[i]), R = exp(lrho), bit bits of precision */
static GEN
homothetie(GEN p, double lrho, long bit)
{
  GEN q, r, t, iR;
  long n = degpol(p), i;

  iR = mygprec(dblexp(-lrho),bit);
  q = mygprec(p, bit);
  r = cgetg(n+3,t_POL); r[1] = p[1];
  t = iR; r[n+2] = q[n+2];
  for (i=n-1; i>0; i--)
  {
    gel(r,i+2) = gmul(t, gel(q,i+2));
    t = mulrr(t, iR);
  }
  gel(r,2) = gmul(t, gel(q,2)); return r;
}

/* change q in 2^(n*e) p(x*2^(-e)), n=deg(q)  [ ~as above with R = 2^-e ]*/
static void
homothetie2n(GEN p, long e)
{
  if (e)
  {
    long i,n = lg(p)-1;
    for (i=2; i<=n; i++) myshiftrc(gel(p,i), (n-i)*e);
  }
}

/* return 2^f * 2^(n*e) p(x*2^(-e)), n=deg(q) */
static void
homothetie_gauss(GEN p, long e, long f)
{
  if (e || f)
  {
    long i, n = lg(p)-1;
    for (i=2; i<=n; i++) gel(p,i) = myshiftic(gel(p,i), f+(n-i)*e);
  }
}

/* Lower bound on the modulus of the largest root z_0
 * k is set to an upper bound for #{z roots, |z-z_0| < eps} */
static double
lower_bound(GEN p, long *k, double eps)
{
  long n = degpol(p), i, j;
  pari_sp ltop = avma;
  GEN a, s, S, ilc;
  double r, R, rho;

  if (n < 4) { *k = n; return 0.; }
  S = cgetg(5,t_VEC);
  a = cgetg(5,t_VEC); ilc = gdiv(real_1(DEFAULTPREC), gel(p,n+2));
  for (i=1; i<=4; i++) gel(a,i) = gmul(ilc,gel(p,n+2-i));
  /* i = 1 split out from next loop for efficiency and initialization */
  s = gel(a,1);
  gel(S,1) = gneg(s); /* Newton sum S_i */
  rho = r = gtodouble(gabs(s,3));
  R = r / n;
  for (i=2; i<=4; i++)
  {
    s = gmulsg(i,gel(a,i));
    for (j=1; j<i; j++) s = gadd(s, gmul(gel(S,j),gel(a,i-j)));
    gel(S,i) = gneg(s); /* Newton sum S_i */
    r = gtodouble(gabs(s,3));
    if (r > 0.)
    {
      r = exp(log(r/n) / (double)i);
      if (r > R) R = r;
    }
  }
  if (R > 0. && eps < 1.2)
    *k = (long)floor((rho/R + n) / (1 + exp(-eps)*cos(eps)));
  else
    *k = n;
  return gc_double(ltop, R);
}

/* return R such that exp(R - tau) <= rho_n(P) <= exp(R + tau)
 * P(0) != 0 and P non constant */
static double
logmax_modulus(GEN p, double tau)
{
  GEN r, q, aux, gunr;
  pari_sp av, ltop = avma;
  long i,k,n=degpol(p),nn,bit,M,e;
  double rho,eps, tau2 = (tau > 3.0)? 0.5: tau/6.;

  r = cgeti(BIGDEFAULTPREC);
  av = avma;

  eps = - 1/log(1.5*tau2); /* > 0 */
  bit = (long) ((double) n*log2(1./tau2)+3*log2((double) n))+1;
  gunr = real_1_bit(bit+2*n);
  aux = gdiv(gunr, gel(p,2+n));
  q = RgX_Rg_mul(p, aux); gel(q,2+n) = gunr;
  e = findpower(q);
  homothetie2n(q,e);
  affsi(e, r);
  q = pol_to_gaussint(q, bit);
  M = (long) (log2( log(4.*n) / (2*tau2) )) + 2;
  nn = n;
  for (i=0,e=0;;)
  { /* nn = deg(q) */
    rho = lower_bound(q, &k, eps);
    if (rho > exp2(-(double)e)) e = (long)-floor(log2(rho));
    affii(shifti(addis(r,e), 1), r);
    if (++i == M) break;

    bit = (long) ((double)k * log2(1./tau2) +
                     (double)(nn-k)*log2(1./eps) + 3*log2((double)nn)) + 1;
    homothetie_gauss(q, e, bit-(long)floor(dbllog2(gel(q,2+nn))+0.5));
    nn -= RgX_valrem(q, &q);
    q = gerepileupto(av, graeffe(q));
    tau2 *= 1.5; if (tau2 > 0.9) tau2 = 0.5;
    eps = -1/log(tau2); /* > 0 */
    e = findpower(q);
  }
  if (!signe(r)) return gc_double(ltop,0.);
  r = itor(r, DEFAULTPREC); shiftr_inplace(r, -M);
  return gc_double(ltop, -rtodbl(r) * M_LN2); /* -log(2) sum e_i 2^-i */
}

static GEN
RgX_normalize1(GEN x)
{
  long i, n = lg(x)-1;
  GEN y;
  for (i = n; i > 1; i--)
    if (!gequal0( gel(x,i) )) break;
  if (i == n) return x;
  pari_warn(warner,"normalizing a polynomial with 0 leading term");
  if (i == 1) pari_err_ROOTS0("roots");
  y = cgetg(i+1, t_POL); y[1] = x[1];
  for (; i > 1; i--) gel(y,i) = gel(x,i);
  return y;
}

static GEN
polrootsbound_i(GEN P, double TAU)
{
  pari_sp av = avma;
  double d;
  (void)RgX_valrem_inexact(P,&P);
  P = RgX_normalize1(P);
  switch(degpol(P))
  {
    case -1: pari_err_ROOTS0("roots");
    case 0:  set_avma(av); return gen_0;
  }
  d = logmax_modulus(P, TAU) + TAU;
  /* not dblexp: result differs on ARM emulator */
  return gerepileuptoleaf(av, mpexp(dbltor(d)));
}
GEN
polrootsbound(GEN P, GEN tau)
{
  if (typ(P) != t_POL) pari_err_TYPE("polrootsbound",P);
  checkvalidpol(P, "polrootsbound");
  return polrootsbound_i(P, tau? gtodouble(tau): 0.01);
}

/* log of modulus of the smallest root of p, with relative error tau */
static double
logmin_modulus(GEN p, double tau)
{
  pari_sp av = avma;
  if (gequal0(gel(p,2))) return -pariINFINITY;
  return gc_double(av, - logmax_modulus(RgX_recip_i(p),tau));
}

/* return the log of the k-th modulus (ascending order) of p, rel. error tau*/
static double
logmodulus(GEN p, long k, double tau)
{
  GEN q;
  long i, kk = k, imax, n = degpol(p), nn, bit, e;
  pari_sp av, ltop=avma;
  double r, tau2 = tau/6;

  bit = (long)(n * (2. + log2(3.*n/tau2)));
  av = avma;
  q = gprec_w(p, nbits2prec(bit));
  q = RgX_gtofp_bit(q, bit);
  e = newton_polygon(q,k);
  r = (double)e;
  homothetie2n(q,e);
  imax = (long)(log2(3./tau) + log2(log(4.*n)))+1;
  for (i=1; i<imax; i++)
  {
    q = eval_rel_pol(q,bit);
    kk -= RgX_valrem(q, &q);
    nn = degpol(q);

    q = gerepileupto(av, graeffe(q));
    e = newton_polygon(q,kk);
    r += e / exp2((double)i);
    q = RgX_gtofp_bit(q, bit);
    homothetie2n(q,e);

    tau2 *= 1.5; if (tau2 > 1.) tau2 = 1.;
    bit = 1 + (long)(nn*(2. + log2(3.*nn/tau2)));
  }
  return gc_double(ltop, -r * M_LN2);
}

/* return the log of the k-th modulus r_k of p, rel. error tau, knowing that
 * rmin < r_k < rmax. This information helps because we may reduce precision
 * quicker */
static double
logpre_modulus(GEN p, long k, double tau, double lrmin, double lrmax)
{
  GEN q;
  long n = degpol(p), i, imax, imax2, bit;
  pari_sp ltop = avma, av;
  double lrho, aux, tau2 = tau/6.;

  aux = (lrmax - lrmin) / 2. + 4*tau2;
  imax = (long) log2(log((double)n)/ aux);
  if (imax <= 0) return logmodulus(p,k,tau);

  lrho  = (lrmin + lrmax) / 2;
  av = avma;
  bit = (long)(n*(2. + aux / M_LN2 - log2(tau2)));
  q = homothetie(p, lrho, bit);
  imax2 = (long)(log2(3./tau * log(4.*n))) + 1;
  if (imax > imax2) imax = imax2;

  for (i=0; i<imax; i++)
  {
    q = eval_rel_pol(q,bit);
    q = gerepileupto(av, graeffe(q));
    aux = 2*aux + 2*tau2;
    tau2 *= 1.5;
    bit = (long)(n*(2. + aux / M_LN2 - log2(1-exp(-tau2))));
    q = RgX_gtofp_bit(q, bit);
  }
  aux = exp2((double)imax);
  return gc_double(ltop, lrho + logmodulus(q,k, aux*tau/3.) / aux);
}

static double
ind_maxlog2(GEN q)
{
  long i, k = -1;
  double L = - pariINFINITY;
  for (i=0; i<=degpol(q); i++)
  {
    double d = dbllog2(gel(q,2+i));
    if (d > L) { L = d; k = i; }
  }
  return k;
}

/* Returns k such that r_k e^(-tau) < R < r_{k+1} e^tau.
 * Assume that l <= k <= n-l */
static long
dual_modulus(GEN p, double lrho, double tau, long l)
{
  long i, imax, delta_k = 0, n = degpol(p), nn, v2, v, bit, ll = l;
  double tau2 = tau * 7./8.;
  pari_sp av = avma;
  GEN q;

  bit = 6*n - 5*l + (long)(n*(-log2(tau2) + tau2 * 8./7.));
  q = homothetie(p, lrho, bit);
  imax = (long)(log(log(2.*n)/tau2)/log(7./4.)+1);

  for (i=0; i<imax; i++)
  {
    q = eval_rel_pol(q,bit); v2 = n - degpol(q);
    v = RgX_valrem(q, &q);
    ll -= maxss(v, v2); if (ll < 0) ll = 0;

    nn = degpol(q); delta_k += v;
    if (!nn) return delta_k;

    q = gerepileupto(av, graeffe(q));
    tau2 *= 7./4.;
    bit = 6*nn - 5*ll + (long)(nn*(-log2(tau2) + tau2 * 8./7.));
  }
  return gc_long(av, delta_k + (long)ind_maxlog2(q));
}

/********************************************************************/
/**                                                                **/
/**              FACTORS THROUGH CIRCLE INTEGRATION                **/
/**                                                                **/
/********************************************************************/
/* l power of 2, W[step*j] = w_j; set f[j] = p(w_j)
 * if inv, w_j = exp(2IPi*j/l), else exp(-2IPi*j/l) */

static void
fft2(GEN W, GEN p, GEN f, long step, long l)
{
  pari_sp av;
  long i, s1, l1, step2;

  if (l == 2)
  {
    gel(f,0) = gadd(gel(p,0), gel(p,step));
    gel(f,1) = gsub(gel(p,0), gel(p,step)); return;
  }
  av = avma;
  l1 = l>>1; step2 = step<<1;
  fft2(W,p,          f,   step2,l1);
  fft2(W,p+step,     f+l1,step2,l1);
  for (i = s1 = 0; i < l1; i++, s1 += step)
  {
    GEN f0 = gel(f,i);
    GEN f1 = gmul(gel(W,s1), gel(f,i+l1));
    gel(f,i)    = gadd(f0, f1);
    gel(f,i+l1) = gsub(f0, f1);
  }
  gerepilecoeffs(av, f, l);
}

static void
fft(GEN W, GEN p, GEN f, long step, long l, long inv)
{
  pari_sp av;
  long i, s1, l1, l2, l3, step4;
  GEN f1, f2, f3, f02;

  if (l == 2)
  {
    gel(f,0) = gadd(gel(p,0), gel(p,step));
    gel(f,1) = gsub(gel(p,0), gel(p,step)); return;
  }
  av = avma;
  if (l == 4)
  {
    pari_sp av2;
    f1 = gadd(gel(p,0), gel(p,step<<1));
    f2 = gsub(gel(p,0), gel(p,step<<1));
    f3 = gadd(gel(p,step), gel(p,3*step));
    f02 = gsub(gel(p,step), gel(p,3*step));
    f02 = inv? mulcxI(f02): mulcxmI(f02);
    av2 = avma;
    gel(f,0) = gadd(f1, f3);
    gel(f,1) = gadd(f2, f02);
    gel(f,2) = gsub(f1, f3);
    gel(f,3) = gsub(f2, f02);
    gerepileallsp(av,av2,4,&gel(f,0),&gel(f,1),&gel(f,2),&gel(f,3));
    return;
  }
  l1 = l>>2; l2 = 2*l1; l3 = l1+l2; step4 = step<<2;
  fft(W,p,          f,   step4,l1,inv);
  fft(W,p+step,     f+l1,step4,l1,inv);
  fft(W,p+(step<<1),f+l2,step4,l1,inv);
  fft(W,p+3*step,   f+l3,step4,l1,inv);
  for (i = s1 = 0; i < l1; i++, s1 += step)
  {
    long s2 = s1 << 1, s3 = s1 + s2;
    GEN g02, g13, f13;
    f1 = gmul(gel(W,s1), gel(f,i+l1));
    f2 = gmul(gel(W,s2), gel(f,i+l2));
    f3 = gmul(gel(W,s3), gel(f,i+l3));

    f02 = gadd(gel(f,i),f2);
    g02 = gsub(gel(f,i),f2);
    f13 = gadd(f1,f3);
    g13 = gsub(f1,f3); g13 = inv? mulcxI(g13): mulcxmI(g13);

    gel(f,i)    = gadd(f02, f13);
    gel(f,i+l1) = gadd(g02, g13);
    gel(f,i+l2) = gsub(f02, f13);
    gel(f,i+l3) = gsub(g02, g13);
  }
  gerepilecoeffs(av, f, l);
}

#define code(t1,t2) ((t1 << 6) | t2)

static GEN
FFT_i(GEN W, GEN x)
{
  long i, l = lg(W), n = lg(x), tx = typ(x), tw, pa;
  GEN y, z, p, pol;
  if (l==1 || ((l-1) & (l-2))) pari_err_DIM("fft");
  tw = RgV_type(W, &p, &pol, &pa);
  if (tx == t_POL) { x++; n--; }
  else if (!is_vec_t(tx)) pari_err_TYPE("fft",x);
  if (n > l) pari_err_DIM("fft");
  if (n < l) {
    z = cgetg(l, t_VECSMALL); /* cf stackdummy */
    for (i = 1; i < n; i++) gel(z,i) = gel(x,i);
    for (     ; i < l; i++) gel(z,i) = gen_0;
  }
  else z = x;
  if (l == 2) return mkveccopy(gel(z,1));
  y = cgetg(l, t_VEC);
  if (tw==code(t_COMPLEX,t_INT) || tw==code(t_COMPLEX,t_REAL))
  {
    long inv = (l >= 5 && signe(imag_i(gel(W,1+(l>>2))))==1) ? 1 : 0;
    fft(W+1, z+1, y+1, 1, l-1, inv);
  } else
    fft2(W+1, z+1, y+1, 1, l-1);
  return y;
}

#undef code

GEN
FFT(GEN W, GEN x)
{
  if (!is_vec_t(typ(W))) pari_err_TYPE("fft",W);
  return FFT_i(W, x);
}

GEN
FFTinv(GEN W, GEN x)
{
  long l = lg(W), i;
  GEN w;
  if (!is_vec_t(typ(W))) pari_err_TYPE("fft",W);
  if (l==1 || ((l-1) & (l-2))) pari_err_DIM("fft");
  w = cgetg(l, t_VECSMALL); /* cf stackdummy */
  gel(w,1) = gel(W,1); /* w = gconj(W), faster */
  for (i = 2; i < l; i++) gel(w, i) = gel(W, l-i+1);
  return FFT_i(w, x);
}

/* returns 1 if p has only real coefficients, 0 else */
static int
isreal(GEN p)
{
  long i;
  for (i = lg(p)-1; i > 1; i--)
    if (typ(gel(p,i)) == t_COMPLEX) return 0;
  return 1;
}

/* x non complex */
static GEN
abs_update_r(GEN x, double *mu) {
  GEN y = gtofp(x, DEFAULTPREC);
  double ly = mydbllogr(y); if (ly < *mu) *mu = ly;
  setabssign(y); return y;
}
/* return |x|, low accuracy. Set *mu = min(log(y), *mu) */
static GEN
abs_update(GEN x, double *mu) {
  GEN y, xr, yr;
  double ly;
  if (typ(x) != t_COMPLEX) return abs_update_r(x, mu);
  xr = gel(x,1);
  yr = gel(x,2);
  if (gequal0(xr)) return abs_update_r(yr,mu);
  if (gequal0(yr)) return abs_update_r(xr,mu);
  /* have to treat 0 specially: 0E-10 + 1e-20 = 0E-10 */
  xr = gtofp(xr, DEFAULTPREC);
  yr = gtofp(yr, DEFAULTPREC);
  y = sqrtr(addrr(sqrr(xr), sqrr(yr)));
  ly = mydbllogr(y); if (ly < *mu) *mu = ly;
  return y;
}

static void
initdft(GEN *Omega, GEN *prim, long N, long Lmax, long bit)
{
  long prec = nbits2prec(bit);
  *Omega = grootsof1(Lmax, prec) + 1;
  *prim = rootsof1u_cx(N, prec);
}

static void
parameters(GEN p, long *LMAX, double *mu, double *gamma,
           int polreal, double param, double param2)
{
  GEN q, pc, Omega, A, RU, prim, g, TWO;
  long n = degpol(p), bit, NN, K, i, j, Lmax;
  pari_sp av2, av = avma;

  bit = gexpo(p) + (long)param2+8;
  Lmax = 4; while (Lmax <= n) Lmax <<= 1;
  NN = (long)(param*3.14)+1; if (NN < Lmax) NN = Lmax;
  K = NN/Lmax; if (K & 1) K++;
  NN = Lmax*K;
  if (polreal) K = K/2+1;

  initdft(&Omega, &prim, NN, Lmax, bit);
  q = mygprec(p,bit) + 2;
  A = cgetg(Lmax+1,t_VEC); A++;
  pc= cgetg(Lmax+1,t_VEC); pc++;
  for (i=0; i <= n; i++) gel(pc,i)= gel(q,i);
  for (   ; i<Lmax; i++) gel(pc,i) = gen_0;

  *mu = pariINFINITY;
  g = real_0_bit(-bit);
  TWO = real2n(1, DEFAULTPREC);
  av2 = avma;
  RU = gen_1;
  for (i=0; i<K; i++)
  {
    if (i) {
      GEN z = RU;
      for (j=1; j<n; j++)
      {
        gel(pc,j) = gmul(gel(q,j),z);
        z = gmul(z,RU); /* RU = prim^i, z=prim^(ij) */
      }
      gel(pc,n) = gmul(gel(q,n),z);
    }

    fft(Omega,pc,A,1,Lmax,1);
    if (polreal && i>0 && i<K-1)
      for (j=0; j<Lmax; j++) g = addrr(g, divrr(TWO, abs_update(gel(A,j),mu)));
    else
      for (j=0; j<Lmax; j++) g = addrr(g, invr(abs_update(gel(A,j),mu)));
    RU = gmul(RU, prim);
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"parameters");
      gerepileall(av2,2, &g,&RU);
    }
  }
  *gamma = mydbllog2r(divru(g,NN));
  *LMAX = Lmax; set_avma(av);
}

/* NN is a multiple of Lmax */
static void
dft(GEN p, long k, long NN, long Lmax, long bit, GEN F, GEN H, long polreal)
{
  GEN Omega, q, qd, pc, pd, A, B, C, RU, aux, U, W, prim, prim2;
  long n = degpol(p), i, j, K;
  pari_sp ltop;

  initdft(&Omega, &prim, NN, Lmax, bit);
  RU = cgetg(n+2,t_VEC) + 1;

  K = NN/Lmax; if (polreal) K = K/2+1;
  q = mygprec(p,bit);
  qd = RgX_deriv(q);

  A = cgetg(Lmax+1,t_VEC); A++;
  B = cgetg(Lmax+1,t_VEC); B++;
  C = cgetg(Lmax+1,t_VEC); C++;
  pc = cgetg(Lmax+1,t_VEC); pc++;
  pd = cgetg(Lmax+1,t_VEC); pd++;
  gel(pc,0) = gel(q,2);  for (i=n+1; i<Lmax; i++) gel(pc,i) = gen_0;
  gel(pd,0) = gel(qd,2); for (i=n;   i<Lmax; i++) gel(pd,i) = gen_0;

  ltop = avma;
  W = cgetg(k+1,t_VEC);
  U = cgetg(k+1,t_VEC);
  for (i=1; i<=k; i++) gel(W,i) = gel(U,i) = gen_0;

  gel(RU,0) = gen_1;
  prim2 = gen_1;
  for (i=0; i<K; i++)
  {
    gel(RU,1) = prim2;
    for (j=1; j<n; j++) gel(RU,j+1) = gmul(gel(RU,j),prim2);
    /* RU[j] = prim^(ij)= prim2^j */

    for (j=1; j<n; j++) gel(pd,j) = gmul(gel(qd,j+2),gel(RU,j));
    fft(Omega,pd,A,1,Lmax,1);
    for (j=1; j<=n; j++) gel(pc,j) = gmul(gel(q,j+2),gel(RU,j));
    fft(Omega,pc,B,1,Lmax,1);
    for (j=0; j<Lmax; j++) gel(C,j) = ginv(gel(B,j));
    for (j=0; j<Lmax; j++) gel(B,j) = gmul(gel(A,j),gel(C,j));
    fft(Omega,B,A,1,Lmax,1);
    fft(Omega,C,B,1,Lmax,1);

    if (polreal) /* p has real coefficients */
    {
      if (i>0 && i<K-1)
      {
        for (j=1; j<=k; j++)
        {
          gel(W,j) = gadd(gel(W,j), gshift(mulreal(gel(A,j+1),gel(RU,j+1)),1));
          gel(U,j) = gadd(gel(U,j), gshift(mulreal(gel(B,j),gel(RU,j)),1));
        }
      }
      else
      {
        for (j=1; j<=k; j++)
        {
          gel(W,j) = gadd(gel(W,j), mulreal(gel(A,j+1),gel(RU,j+1)));
          gel(U,j) = gadd(gel(U,j), mulreal(gel(B,j),gel(RU,j)));
        }
      }
    }
    else
    {
      for (j=1; j<=k; j++)
      {
        gel(W,j) = gadd(gel(W,j), gmul(gel(A,j+1),gel(RU,j+1)));
        gel(U,j) = gadd(gel(U,j), gmul(gel(B,j),gel(RU,j)));
      }
    }
    prim2 = gmul(prim2,prim);
    gerepileall(ltop,3, &W,&U,&prim2);
  }

  for (i=1; i<=k; i++)
  {
    aux=gel(W,i);
    for (j=1; j<i; j++) aux = gadd(aux, gmul(gel(W,i-j),gel(F,k+2-j)));
    gel(F,k+2-i) = gdivgs(aux,-i*NN);
  }
  for (i=0; i<k; i++)
  {
    aux=gel(U,k-i);
    for (j=1+i; j<k; j++) aux = gadd(aux,gmul(gel(F,2+j),gel(U,j-i)));
    gel(H,i+2) = gdivgu(aux,NN);
  }
}

#define NEWTON_MAX 10
static GEN
refine_H(GEN F, GEN G, GEN HH, long bit, long Sbit)
{
  GEN H = HH, D, aux;
  pari_sp ltop = avma;
  long error, i, bit1, bit2;

  D = Rg_RgX_sub(gen_1, RgX_rem(RgX_mul(H,G),F)); error = gexpo(D);
  bit2 = bit + Sbit;
  for (i=0; error>-bit && i<NEWTON_MAX && error<=0; i++)
  {
    if (gc_needed(ltop,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"refine_H");
      gerepileall(ltop,2, &D,&H);
    }
    bit1 = -error + Sbit;
    aux = RgX_mul(mygprec(H,bit1), mygprec(D,bit1));
    aux = RgX_rem(mygprec(aux,bit1), mygprec(F,bit1));

    bit1 = -error*2 + Sbit; if (bit1 > bit2) bit1 = bit2;
    H = RgX_add(mygprec(H,bit1), aux);
    D = Rg_RgX_sub(gen_1, RgX_rem(RgX_mul(H,G),F));
    error = gexpo(D); if (error < -bit1) error = -bit1;
  }
  if (error > -bit/2) return NULL; /* FAIL */
  return gerepilecopy(ltop,H);
}

/* return 0 if fails, 1 else */
static long
refine_F(GEN p, GEN *F, GEN *G, GEN H, long bit, double gamma)
{
  GEN f0, FF, GG, r, HH = H;
  long error, i, bit1 = 0, bit2, Sbit, Sbit2,  enh, normF, normG, n = degpol(p);
  pari_sp av = avma;

  FF = *F; GG = RgX_divrem(p, FF, &r);
  error = gexpo(r); if (error <= -bit) error = 1-bit;
  normF = gexpo(FF);
  normG = gexpo(GG);
  enh = gexpo(H); if (enh < 0) enh = 0;
  Sbit = normF + 2*normG + enh + (long)(4.*log2((double)n)+gamma) + 1;
  Sbit2 = enh + 2*(normF+normG) + (long)(2.*gamma+5.*log2((double)n)) + 1;
  bit2 = bit + Sbit;
  for (i=0; error>-bit && i<NEWTON_MAX && error<=0; i++)
  {
    if (bit1 == bit2 && i >= 2) { Sbit += n; Sbit2 += n; bit2 += n; }
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"refine_F");
      gerepileall(av,4, &FF,&GG,&r,&HH);
    }

    bit1 = -error + Sbit2;
    HH = refine_H(mygprec(FF,bit1), mygprec(GG,bit1), mygprec(HH,bit1),
                  1-error, Sbit2);
    if (!HH) return 0; /* FAIL */

    bit1 = -error + Sbit;
    r = RgX_mul(mygprec(HH,bit1), mygprec(r,bit1));
    f0 = RgX_rem(mygprec(r,bit1), mygprec(FF,bit1));

    bit1 = -2*error + Sbit; if (bit1 > bit2) bit1 = bit2;
    FF = gadd(mygprec(FF,bit1),f0);

    bit1 = -3*error + Sbit; if (bit1 > bit2) bit1 = bit2;
    GG = RgX_divrem(mygprec(p,bit1), mygprec(FF,bit1), &r);
    error = gexpo(r); if (error < -bit1) error = -bit1;
  }
  if (error>-bit) return 0; /* FAIL */
  *F = FF; *G = GG; return 1;
}

/* returns F and G from the unit circle U such that |p-FG|<2^(-bit) |cd|,
where cd is the leading coefficient of p */
static void
split_fromU(GEN p, long k, double delta, long bit,
            GEN *F, GEN *G, double param, double param2)
{
  GEN pp, FF, GG, H;
  long n = degpol(p), NN, bit2, Lmax;
  int polreal = isreal(p);
  pari_sp ltop;
  double mu, gamma;

  pp = gdiv(p, gel(p,2+n));
  parameters(pp, &Lmax,&mu,&gamma, polreal,param,param2);

  H  = cgetg(k+2,t_POL); H[1] = p[1];
  FF = cgetg(k+3,t_POL); FF[1]= p[1];
  gel(FF,k+2) = gen_1;

  NN = (long)(0.5/delta); NN |= 1; if (NN < 2) NN = 2;
  NN *= Lmax; ltop = avma;
  for(;;)
  {
    bit2 = (long)(((double)NN*delta-mu)/M_LN2) + gexpo(pp) + 8;
    dft(pp, k, NN, Lmax, bit2, FF, H, polreal);
    if (refine_F(pp,&FF,&GG,H,bit,gamma)) break;
    NN <<= 1; set_avma(ltop);
  }
  *G = gmul(GG,gel(p,2+n)); *F = FF;
}

static void
optimize_split(GEN p, long k, double delta, long bit,
            GEN *F, GEN *G, double param, double param2)
{
  long n = degpol(p);
  GEN FF, GG;

  if (k <= n/2)
    split_fromU(p,k,delta,bit,F,G,param,param2);
  else
  {
    split_fromU(RgX_recip_i(p),n-k,delta,bit,&FF,&GG,param,param2);
    *F = RgX_recip_i(GG);
    *G = RgX_recip_i(FF);
  }
}

/********************************************************************/
/**                                                                **/
/**               SEARCH FOR SEPARATING CIRCLE                     **/
/**                                                                **/
/********************************************************************/

/* return p(2^e*x) *2^(-n*e) */
static void
scalepol2n(GEN p, long e)
{
  long i,n=lg(p)-1;
  for (i=2; i<=n; i++) gel(p,i) = gmul2n(gel(p,i),(i-n)*e);
}

/* returns p(x/R)*R^n; assume R is at the correct accuracy */
static GEN
scalepol(GEN p, GEN R, long bit)
{ return RgX_rescale(mygprec(p, bit), R); }

/* return (conj(a)X-1)^n * p[ (X-a) / (conj(a)X-1) ] */
static GEN
conformal_basecase(GEN p, GEN a)
{
  GEN z, r, ma, ca;
  long i, n = degpol(p);
  pari_sp av;

  if (n <= 0) return p;
  ma = gneg(a); ca = conj_i(a);
  av = avma;
  z = deg1pol_shallow(ca, gen_m1, 0);
  r = scalarpol_shallow(gel(p,2+n), 0);
  for (i=n-1; ; i--)
  {
    r = RgX_addmulXn_shallow(r, gmul(ma,r), 1); /* r *= (X - a) */
    r = gadd(r, gmul(z, gel(p,2+i)));
    if (i == 0) return gerepileupto(av, r);
    z = RgX_addmulXn_shallow(gmul(z,ca), gneg(z), 1); /* z *= conj(a)X - 1 */
    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"conformal_pol (%ld/%ld)",n-i, n);
      gerepileall(av,2, &r,&z);
    }
  }
}
static GEN
conformal_pol(GEN p, GEN a)
{
  pari_sp av = avma;
  long d, nR, n = degpol(p), v;
  GEN Q, R, S, T;
  if (n < 35) return conformal_basecase(p, a);
  d = (n+1) >> 1; v = varn(p);
  Q = RgX_shift_shallow(p, -d);
  R = RgXn_red_shallow(p, d);
  Q = conformal_pol(Q, a);
  R = conformal_pol(R, a);
  S = gpowgs(deg1pol_shallow(gen_1, gneg(a), v), d);
  T = RgX_recip_i(S);
  if (typ(a) == t_COMPLEX) T = gconj(T);
  if (odd(d)) T = RgX_neg(T);
  /* S = (X - a)^d, T = (conj(a) X - 1)^d */
  nR = n - degpol(R) - d; /* >= 0 */
  if (nR) T = RgX_mul(T, gpowgs(deg1pol_shallow(gconj(a), gen_m1, v), nR));
  return gerepileupto(av, RgX_add(RgX_mul(Q, S), RgX_mul(R, T)));
}

static const double UNDEF = -100000.;

static double
logradius(double *radii, GEN p, long k, double aux, double *delta)
{
  long i, n = degpol(p);
  double lrho, lrmin, lrmax;
  if (k > 1)
  {
    i = k-1; while (i>0 && radii[i] == UNDEF) i--;
    lrmin = logpre_modulus(p,k,aux, radii[i], radii[k]);
  }
  else /* k=1 */
    lrmin = logmin_modulus(p,aux);
  radii[k] = lrmin;

  if (k+1<n)
  {
    i = k+2; while (i<=n && radii[i] == UNDEF) i++;
    lrmax = logpre_modulus(p,k+1,aux, radii[k+1], radii[i]);
  }
  else /* k+1=n */
    lrmax = logmax_modulus(p,aux);
  radii[k+1] = lrmax;

  lrho = radii[k];
  for (i=k-1; i>=1; i--)
  {
    if (radii[i] == UNDEF || radii[i] > lrho)
      radii[i] = lrho;
    else
      lrho = radii[i];
  }
  lrho = radii[k+1];
  for (i=k+1; i<=n; i++)
  {
    if (radii[i] == UNDEF || radii[i] < lrho)
      radii[i] = lrho;
    else
      lrho = radii[i];
  }
  *delta = (lrmax - lrmin) / 2;
  if (*delta > 1.) *delta = 1.;
  return (lrmin + lrmax) / 2;
}

static void
update_radius(long n, double *radii, double lrho, double *par, double *par2)
{
  double t, param = 0., param2 = 0.;
  long i;
  for (i=1; i<=n; i++)
  {
    radii[i] -= lrho;
    t = fabs(rtodbl( invr(subsr(1, dblexp(radii[i]))) ));
    param += t; if (t > 1.) param2 += log2(t);
  }
  *par = param; *par2 = param2;
}

/* apply the conformal mapping then split from U */
static void
conformal_mapping(double *radii, GEN ctr, GEN p, long k, long bit,
                  double aux, GEN *F,GEN *G)
{
  long bit2, n = degpol(p), i;
  pari_sp ltop = avma, av;
  GEN q, FF, GG, a, R;
  double lrho, delta, param, param2;
  /* n * (2.*log2(2.732)+log2(1.5)) + 1 */
  bit2 = bit + (long)(n*3.4848775) + 1;
  a = sqrtr_abs( utor(3, 2*MEDDEFAULTPREC - 2) );
  a = divrs(a, -6);
  a = gmul(mygprec(a,bit2), mygprec(ctr,bit2)); /* a = -ctr/2sqrt(3) */

  av = avma;
  q = conformal_pol(mygprec(p,bit2), a);
  for (i=1; i<=n; i++)
    if (radii[i] != UNDEF) /* update array radii */
    {
      pari_sp av2 = avma;
      GEN t, r = dblexp(radii[i]), r2 = sqrr(r);
      /* 2(r^2 - 1) / (r^2 - 3(r-1)) */
      t = divrr(shiftr((subrs(r2,1)),1), subrr(r2, mulur(3,subrs(r,1))));
      radii[i] = mydbllogr(addsr(1,t)) / 2;
      set_avma(av2);
    }
  lrho = logradius(radii, q,k,aux/10., &delta);
  update_radius(n, radii, lrho, &param, &param2);

  bit2 += (long)(n * fabs(lrho)/M_LN2 + 1.);
  R = mygprec(dblexp(-lrho), bit2);
  q = scalepol(q,R,bit2);
  gerepileall(av,2, &q,&R);

  optimize_split(q,k,delta,bit2,&FF,&GG,param,param2);
  bit2 += n; R = invr(R);
  FF = scalepol(FF,R,bit2);
  GG = scalepol(GG,R,bit2);

  a = mygprec(a,bit2);
  FF = conformal_pol(FF,a);
  GG = conformal_pol(GG,a);

  a = invr(subsr(1, gnorm(a)));
  FF = RgX_Rg_mul(FF, powru(a,k));
  GG = RgX_Rg_mul(GG, powru(a,n-k));

  *F = mygprec(FF,bit+n);
  *G = mygprec(GG,bit+n); gerepileall(ltop,2, F,G);
}

/* split p, this time without scaling. returns in F and G two polynomials
 * such that |p-FG|< 2^(-bit)|p| */
static void
split_2(GEN p, long bit, GEN ctr, double thickness, GEN *F, GEN *G)
{
  GEN q, FF, GG, R;
  double aux, delta, param, param2;
  long n = degpol(p), i, j, k, bit2;
  double lrmin, lrmax, lrho, *radii;

  radii = (double*) stack_malloc_align((n+1) * sizeof(double), sizeof(double));

  for (i=2; i<n; i++) radii[i] = UNDEF;
  aux = thickness/(double)(4 * n);
  lrmin = logmin_modulus(p, aux);
  lrmax = logmax_modulus(p, aux);
  radii[1] = lrmin;
  radii[n] = lrmax;
  i = 1; j = n;
  lrho = (lrmin + lrmax) / 2;
  k = dual_modulus(p, lrho, aux, 1);
  if (5*k < n || (n < 2*k && 5*k < 4*n))
    { lrmax = lrho; j=k+1; radii[j] = lrho; }
  else
    { lrmin = lrho; i=k;   radii[i] = lrho; }
  while (j > i+1)
  {
    if (i+j == n+1)
      lrho = (lrmin + lrmax) / 2;
    else
    {
      double kappa = 2. - log(1. + minss(i,n-j)) / log(1. + minss(j,n-i));
      if (i+j < n+1) lrho = lrmax * kappa + lrmin;
      else           lrho = lrmin * kappa + lrmax;
      lrho /= 1+kappa;
    }
    aux = (lrmax - lrmin) / (4*(j-i));
    k = dual_modulus(p, lrho, aux, minss(i,n+1-j));
    if (k-i < j-k-1 || (k-i == j-k-1 && 2*k > n))
      { lrmax = lrho; j=k+1; radii[j] = lrho - aux; }
    else
      { lrmin = lrho; i=k;   radii[i] = lrho + aux; }
  }
  aux = lrmax - lrmin;

  if (ctr)
  {
    lrho = (lrmax + lrmin) / 2;
    for (i=1; i<=n; i++)
      if (radii[i] != UNDEF) radii[i] -= lrho;

    bit2 = bit + (long)(n * fabs(lrho)/M_LN2 + 1.);
    R = mygprec(dblexp(-lrho), bit2);
    q = scalepol(p,R,bit2);
    conformal_mapping(radii, ctr, q, k, bit2, aux, &FF, &GG);
  }
  else
  {
    lrho = logradius(radii, p, k, aux/10., &delta);
    update_radius(n, radii, lrho, &param, &param2);

    bit2 = bit + (long)(n * fabs(lrho)/M_LN2 + 1.);
    R = mygprec(dblexp(-lrho), bit2);
    q = scalepol(p,R,bit2);
    optimize_split(q, k, delta, bit2, &FF, &GG, param, param2);
  }
  bit  += n;
  bit2 += n; R = invr(mygprec(R,bit2));
  *F = mygprec(scalepol(FF,R,bit2), bit);
  *G = mygprec(scalepol(GG,R,bit2), bit);
}

/* procedure corresponding to steps 5,6,.. page 44 in RR n. 1852 */
/* put in F and G two polynomial such that |p-FG|<2^(-bit)|p|
 * where the maximum modulus of the roots of p is <=1.
 * Assume sum of roots is 0. */
static void
split_1(GEN p, long bit, GEN *F, GEN *G)
{
  long i, imax, n = degpol(p), polreal = isreal(p), ep = gexpo(p), bit2 = bit+n;
  GEN ctr, q, qq, FF, GG, v, gr, r, newq;
  double lrmin, lrmax, lthick;
  const double LOG3 = 1.098613;

  lrmax = logmax_modulus(p, 0.01);
  gr = mygprec(dblexp(-lrmax), bit2);
  q = scalepol(p,gr,bit2);

  bit2 = bit + gexpo(q) - ep + (long)((double)n*2.*log2(3.)+1);
  v = cgetg(5,t_VEC);
  gel(v,1) = gen_2;
  gel(v,2) = gen_m2;
  gel(v,3) = mkcomplex(gen_0, gel(v,1));
  gel(v,4) = mkcomplex(gen_0, gel(v,2));
  q = mygprec(q,bit2); lthick = 0;
  newq = ctr = NULL; /* -Wall */
  imax = polreal? 3: 4;
  for (i=1; i<=imax; i++)
  {
    qq = RgX_translate(q, gel(v,i));
    lrmin = logmin_modulus(qq,0.05);
    if (LOG3 > lrmin + lthick)
    {
      double lquo = logmax_modulus(qq,0.05) - lrmin;
      if (lquo > lthick) { lthick = lquo; newq = qq; ctr = gel(v,i); }
    }
    if (lthick > M_LN2) break;
    if (polreal && i==2 && lthick > LOG3 - M_LN2) break;
  }
  bit2 = bit + gexpo(newq) - ep + (long)(n*LOG3/M_LN2 + 1);
  split_2(newq, bit2, ctr, lthick, &FF, &GG);
  r = gneg(mygprec(ctr,bit2));
  FF = RgX_translate(FF,r);
  GG = RgX_translate(GG,r);

  gr = invr(gr); bit2 = bit - ep + gexpo(FF)+gexpo(GG);
  *F = scalepol(FF,gr,bit2);
  *G = scalepol(GG,gr,bit2);
}

/* put in F and G two polynomials such that |P-FG|<2^(-bit)|P|,
where the maximum modulus of the roots of p is < 0.5 */
static int
split_0_2(GEN p, long bit, GEN *F, GEN *G)
{
  GEN q, b;
  long n = degpol(p), k, bit2, eq;
  double aux0 = dbllog2(gel(p,n+2)); /* != -oo */
  double aux1 = dbllog2(gel(p,n+1)), aux;

  if (aux1 == -pariINFINITY) /* p1 = 0 */
    aux = 0;
  else
  {
    aux = aux1 - aux0; /* log2(p1/p0) */
    /* beware double overflow */
    if (aux >= 0 && (aux > 1e4 || exp2(aux) > 2.5*n)) return 0;
    aux = (aux < -300)? 0.: n*log2(1 + exp2(aux)/(double)n);
  }
  bit2 = bit+1 + (long)(log2((double)n) + aux);
  q = mygprec(p,bit2);
  if (aux1 == -pariINFINITY) b = NULL;
  else
  {
    b = gdivgs(gdiv(gel(q,n+1),gel(q,n+2)),-n);
    q = RgX_translate(q,b);
  }
  gel(q,n+1) = gen_0; eq = gexpo(q);
  k = 0;
  while (k <= n/2 && (- gexpo(gel(q,k+2)) > bit2 + 2*(n-k) + eq
                      || gequal0(gel(q,k+2)))) k++;
  if (k > 0)
  {
    if (k > n/2) k = n/2;
    bit2 += k<<1;
    *F = pol_xn(k, 0);
    *G = RgX_shift_shallow(q, -k);
  }
  else
  {
    split_1(q,bit2,F,G);
    bit2 = bit + gexpo(*F) + gexpo(*G) - gexpo(p) + (long)aux+1;
    *F = mygprec(*F,bit2);
  }
  *G = mygprec(*G,bit2);
  if (b)
  {
    GEN mb = mygprec(gneg(b), bit2);
    *F = RgX_translate(*F, mb);
    *G = RgX_translate(*G, mb);
  }
  return 1;
}

/* put in F and G two polynomials such that |P-FG|<2^(-bit)|P|.
 * Assume max_modulus(p) < 2 */
static void
split_0_1(GEN p, long bit, GEN *F, GEN *G)
{
  GEN FF, GG;
  long n, bit2, normp;

  if  (split_0_2(p,bit,F,G)) return;

  normp = gexpo(p);
  scalepol2n(p,2); /* p := 4^(-n) p(4*x) */
  n = degpol(p); bit2 = bit + 2*n + gexpo(p) - normp;
  split_1(mygprec(p,bit2), bit2,&FF,&GG);
  scalepol2n(FF,-2);
  scalepol2n(GG,-2); bit2 = bit + gexpo(FF) + gexpo(GG) - normp;
  *F = mygprec(FF,bit2);
  *G = mygprec(GG,bit2);
}

/* put in F and G two polynomials such that |P-FG|<2^(-bit)|P| */
static void
split_0(GEN p, long bit, GEN *F, GEN *G)
{
  const double LOG1_9 = 0.6418539;
  long n = degpol(p), k = 0;
  GEN q;

  while (gexpo(gel(p,k+2)) < -bit && k <= n/2) k++;
  if (k > 0)
  {
    if (k > n/2) k = n/2;
    *F = pol_xn(k, 0);
    *G = RgX_shift_shallow(p, -k);
  }
  else
  {
    double lr = logmax_modulus(p, 0.05);
    if (lr < LOG1_9) split_0_1(p, bit, F, G);
    else
    {
      q = RgX_recip_i(p);
      lr = logmax_modulus(q,0.05);
      if (lr < LOG1_9)
      {
        split_0_1(q, bit, F, G);
        *F = RgX_recip_i(*F);
        *G = RgX_recip_i(*G);
      }
      else
        split_2(p,bit,NULL, 1.2837,F,G);
    }
  }
}

/********************************************************************/
/**                                                                **/
/**                ERROR ESTIMATE FOR THE ROOTS                    **/
/**                                                                **/
/********************************************************************/

static GEN
root_error(long n, long k, GEN roots_pol, long err, GEN shatzle)
{
  GEN rho, d, eps, epsbis, eps2, aux, rap = NULL;
  long i, j;

  d = cgetg(n+1,t_VEC);
  for (i=1; i<=n; i++)
  {
    if (i!=k)
    {
      aux = gsub(gel(roots_pol,i), gel(roots_pol,k));
      gel(d,i) = gabs(mygprec(aux,31), DEFAULTPREC);
    }
  }
  rho = gabs(mygprec(gel(roots_pol,k),31), DEFAULTPREC);
  if (expo(rho) < 0) rho = real_1(DEFAULTPREC);
  eps = mulrr(rho, shatzle);
  aux = shiftr(powru(rho,n), err);

  for (j=1; j<=2 || (j<=5 && cmprr(rap, dbltor(1.2)) > 0); j++)
  {
    GEN prod = NULL; /* 1. */
    long m = n;
    epsbis = mulrr(eps, dbltor(1.25));
    for (i=1; i<=n; i++)
    {
      if (i != k && cmprr(gel(d,i),epsbis) > 0)
      {
        GEN dif = subrr(gel(d,i),eps);
        prod = prod? mulrr(prod, dif): dif;
        m--;
      }
    }
    eps2 = prod? divrr(aux, prod): aux;
    if (m > 1) eps2 = sqrtnr(shiftr(eps2, 2*m-2), m);
    rap = divrr(eps,eps2); eps = eps2;
  }
  return eps;
}

/* round a complex or real number x to an absolute value of 2^(-bit) */
static GEN
mygprec_absolute(GEN x, long bit)
{
  long e;
  GEN y;

  switch(typ(x))
  {
    case t_REAL:
      e = expo(x) + bit;
      return (e <= 0 || !signe(x))? real_0_bit(-bit): rtor(x, nbits2prec(e));
    case t_COMPLEX:
      if (gexpo(gel(x,2)) < -bit) return mygprec_absolute(gel(x,1),bit);
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = mygprec_absolute(gel(x,1),bit);
      gel(y,2) = mygprec_absolute(gel(x,2),bit);
      return y;
    default: return x;
  }
}

static long
a_posteriori_errors(GEN p, GEN roots_pol, long err)
{
  long i, n = degpol(p), e_max = -(long)EXPOBITS;
  GEN sigma, shatzle;

  err += (long)log2((double)n) + 1;
  if (err > -2) return 0;
  sigma = real2n(-err, LOWDEFAULTPREC);
  /*  2 / ((s - 1)^(1/n) - 1) */
  shatzle = divur(2, subrs(sqrtnr(subrs(sigma,1),n), 1));
  for (i=1; i<=n; i++)
  {
    pari_sp av = avma;
    GEN x = root_error(n,i,roots_pol,err,shatzle);
    long e = gexpo(x);
    set_avma(av); if (e > e_max) e_max = e;
    gel(roots_pol,i) = mygprec_absolute(gel(roots_pol,i), -e);
  }
  return e_max;
}

/********************************************************************/
/**                                                                **/
/**                           MAIN                                 **/
/**                                                                **/
/********************************************************************/
static GEN
append_clone(GEN r, GEN a) { a = gclone(a); vectrunc_append(r, a); return a; }

/* put roots in placeholder roots_pol so that |P - L_1...L_n| < 2^(-bit)|P|
 * returns prod (x-roots_pol[i]) */
static GEN
split_complete(GEN p, long bit, GEN roots_pol)
{
  long n = degpol(p);
  pari_sp ltop;
  GEN p1, F, G, a, b, m1, m2;

  if (n == 1)
  {
    a = gneg_i(gdiv(gel(p,2), gel(p,3)));
    (void)append_clone(roots_pol,a); return p;
  }
  ltop = avma;
  if (n == 2)
  {
    F = gsub(gsqr(gel(p,3)), gmul2n(gmul(gel(p,2),gel(p,4)), 2));
    F = gsqrt(F, nbits2prec(bit));
    p1 = ginv(gmul2n(gel(p,4),1));
    a = gneg_i(gmul(gadd(F,gel(p,3)), p1));
    b =        gmul(gsub(F,gel(p,3)), p1);
    a = append_clone(roots_pol,a);
    b = append_clone(roots_pol,b); set_avma(ltop);
    a = mygprec(a, 3*bit);
    b = mygprec(b, 3*bit);
    return gmul(gel(p,4), mkpoln(3, gen_1, gneg(gadd(a,b)), gmul(a,b)));
  }
  split_0(p,bit,&F,&G);
  m1 = split_complete(F,bit,roots_pol);
  m2 = split_complete(G,bit,roots_pol);
  return gerepileupto(ltop, gmul(m1,m2));
}

static GEN
quicktofp(GEN x)
{
  const long prec = DEFAULTPREC;
  switch(typ(x))
  {
    case t_INT: return itor(x, prec);
    case t_REAL: return rtor(x, prec);
    case t_FRAC: return fractor(x, prec);
    case t_COMPLEX: {
      GEN a = gel(x,1), b = gel(x,2);
      /* avoid problem with 0, e.g. x = 0 + I*1e-100. We don't want |x| = 0. */
      if (isintzero(a)) return cxcompotor(b, prec);
      if (isintzero(b)) return cxcompotor(a, prec);
      a = cxcompotor(a, prec);
      b = cxcompotor(b, prec); return sqrtr(addrr(sqrr(a), sqrr(b)));
    }
    default: pari_err_TYPE("quicktofp",x);
      return NULL;/*LCOV_EXCL_LINE*/
  }

}

/* bound log_2 |largest root of p| (Fujiwara's bound) */
double
fujiwara_bound(GEN p)
{
  pari_sp av = avma;
  long i, n = degpol(p);
  GEN cc;
  double loglc, Lmax;

  if (n <= 0) pari_err_CONSTPOL("fujiwara_bound");
  loglc = mydbllog2r( quicktofp(gel(p,n+2)) ); /* log_2 |lc(p)| */
  cc = gel(p, 2);
  if (gequal0(cc))
    Lmax = -pariINFINITY-1;
  else
    Lmax = (mydbllog2r(quicktofp(cc)) - loglc - 1) / n;
  for (i = 1; i < n; i++)
  {
    GEN y = gel(p,i+2);
    double L;
    if (gequal0(y)) continue;
    L = (mydbllog2r(quicktofp(y)) - loglc) / (n-i);
    if (L > Lmax) Lmax = L;
  }
  return gc_double(av, Lmax+1);
}

/* Fujiwara's bound, real roots. Based on the following remark: if
 *   p = x^n + sum a_i x^i and q = x^n + sum min(a_i,0)x^i
 * then for all x >= 0, p(x) >= q(x). Thus any bound for the (positive) roots
 * of q is a bound for the positive roots of p. */
double
fujiwara_bound_real(GEN p, long sign)
{
  pari_sp av = avma;
  GEN x;
  long n = degpol(p), i, signodd, signeven;
  if (n <= 0) pari_err_CONSTPOL("fujiwara_bound");
  x = shallowcopy(p);
  if (gsigne(gel(x, n+2)) > 0)
  { signeven = 1; signodd = sign; }
  else
  { signeven = -1; signodd = -sign; }
  for (i = 0; i < n; i++)
  {
    if ((n - i) % 2)
    { if (gsigne(gel(x, i+2)) == signodd ) gel(x, i+2) = gen_0; }
    else
    { if (gsigne(gel(x, i+2)) == signeven) gel(x, i+2) = gen_0; }
  }
  return gc_double(av, fujiwara_bound(x));
}

static GEN
mygprecrc_special(GEN x, long prec, long e)
{
  GEN y;
  switch(typ(x))
  {
    case t_REAL:
      if (!signe(x)) return real_0_bit(minss(e, expo(x)));
      return (prec > realprec(x))? rtor(x, prec): x;
    case t_COMPLEX:
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = mygprecrc_special(gel(x,1),prec,e);
      gel(y,2) = mygprecrc_special(gel(x,2),prec,e);
      return y;
    default: return x;
  }
}

/* like mygprec but keep at least the same precision as before */
static GEN
mygprec_special(GEN x, long bit)
{
  long lx, i, e, prec;
  GEN y;

  if (bit < 0) bit = 0; /* should not happen */
  e = gexpo(x) - bit;
  prec = nbits2prec(bit);
  switch(typ(x))
  {
    case t_POL:
      y = cgetg_copy(x, &lx); y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = mygprecrc_special(gel(x,i),prec,e);
      break;

    default: y = mygprecrc_special(x,prec,e);
  }
  return y;
}

static GEN
fix_roots1(GEN R)
{
  long i, l = lg(R);
  GEN v = cgetg(l, t_VEC);
  for (i=1; i < l; i++) { GEN r = gel(R,i); gel(v,i) = gcopy(r); gunclone(r); }
  return v;
}
static GEN
fix_roots(GEN R, long h, long bit)
{
  long i, j, c, n, prec;
  GEN v, Z, gh;

  if (h == 1) return fix_roots1(R);
  prec = nbits2prec(bit); Z = grootsof1(h, prec); gh = utoipos(h);
  n = lg(R)-1; v = cgetg(h*n + 1, t_VEC);
  for (c = i = 1; i <= n; i++)
  {
    GEN s, r = gel(R,i);
    s = (h == 2)? gsqrt(r, prec): gsqrtn(r, gh, NULL, prec);
    for (j = 1; j <= h; j++) gel(v, c++) = gmul(s, gel(Z,j));
    gunclone(r);
  }
  return v;
}

static GEN
all_roots(GEN p, long bit)
{
  long bit2, i, e, h, n = degpol(p), elc = gexpo(leading_coeff(p));
  GEN q, R, m, pd = RgX_deflate_max(p, &h);
  double fb = fujiwara_bound(pd);
  pari_sp av;

  if (fb < 0) fb = 0;
  bit2 = bit + maxss(gexpo(p), 0) + (long)ceil(log2(n / h) + 2 * fb);
  for (av = avma, i = 1, e = 0;; i++, set_avma(av))
  {
    R = vectrunc_init(n+1);
    bit2 += e + (n << i);
    q = RgX_gtofp_bit(mygprec(pd,bit2), bit2);
    q[1] = evalsigne(1)|evalvarn(0);
    m = split_complete(q, bit2, R);
    R = fix_roots(R, h, bit2);
    q = mygprec_special(pd,bit2);
    q[1] = evalsigne(1)|evalvarn(0);
    e = gexpo(RgX_sub(q, m)) - elc + (long)log2((double)n) + 1;
    if (e < 0)
    {
      if (e < -2*bit2) e = -2*bit2; /* avoid e = -oo */
      e = bit + a_posteriori_errors(p, R, e);
      if (e < 0) return R;
    }
    if (DEBUGLEVEL)
      err_printf("all_roots: restarting, i = %ld, e = %ld\n", i,e);
  }
}

INLINE int
isexactscalar(GEN x) { long tx = typ(x); return is_rational_t(tx); }

static int
isexactpol(GEN p)
{
  long i,n = degpol(p);
  for (i=0; i<=n; i++)
    if (!isexactscalar(gel(p,i+2))) return 0;
  return 1;
}

/* p(0) != 0 [for efficiency] */
static GEN
solve_exact_pol(GEN p, long bit)
{
  long i, j, k, m, n = degpol(p), iroots = 0;
  GEN ex, factors, v = zerovec(n);

  factors = ZX_squff(Q_primpart(p), &ex);
  for (i=1; i<lg(factors); i++)
  {
    GEN roots_fact = all_roots(gel(factors,i), bit);
    n = degpol(gel(factors,i));
    m = ex[i];
    for (j=1; j<=n; j++)
      for (k=1; k<=m; k++) v[++iroots] = roots_fact[j];
  }
  return v;
}

/* return the roots of p with absolute error bit */
static GEN
roots_com(GEN q, long bit)
{
  GEN L, p;
  long v = RgX_valrem_inexact(q, &p);
  int ex = isexactpol(p);
  if (!ex) p = RgX_normalize1(p);
  if (lg(p) == 3)
    L = cgetg(1,t_VEC); /* constant polynomial */
  else
    L = ex? solve_exact_pol(p,bit): all_roots(p,bit);
  if (v)
  {
    GEN M, z, t = gel(q,2);
    long i, x, y, l, n;

    if (isrationalzero(t)) x = -bit;
    else
    {
      n = gexpo(t);
      x = n / v; l = degpol(q);
      for (i = v; i <= l; i++)
      {
        t  = gel(q,i+2);
        if (isrationalzero(t)) continue;
        y = (n - gexpo(t)) / i;
        if (y < x) x = y;
      }
    }
    z = real_0_bit(x); l = v + lg(L);
    M = cgetg(l, t_VEC); L -= v;
    for (i = 1; i <= v; i++) gel(M,i) = z;
    for (     ; i <  l; i++) gel(M,i) = gel(L,i);
    L = M;
  }
  return L;
}

static GEN
tocomplex(GEN x, long l, long bit)
{
  GEN y;
  if (typ(x) == t_COMPLEX)
  {
    if (signe(gel(x,1))) return mygprecrc(x, l, -bit);
    x = gel(x,2);
    y = cgetg(3,t_COMPLEX);
    gel(y,1) = real_0_bit(-bit);
    gel(y,2) = mygprecrc(x, l, -bit);
  }
  else
  {
    y = cgetg(3,t_COMPLEX);
    gel(y,1) = mygprecrc(x, l, -bit);
    gel(y,2) = real_0_bit(-bit);
  }
  return y;
}

/* x,y are t_COMPLEX of t_REALs or t_REAL, compare wrt |Im x| - |Im y|,
 * then Re x - Re y, up to 2^-e absolute error */
static int
cmp_complex_appr(void *E, GEN x, GEN y)
{
  long e = (long)E;
  GEN z, xi, yi, xr, yr;
  long sz, sxi, syi;
  if (typ(x) == t_COMPLEX) { xr = gel(x,1); xi = gel(x,2); sxi = signe(xi); }
  else { xr = x; xi = NULL; sxi = 0; }
  if (typ(y) == t_COMPLEX) { yr = gel(y,1); yi = gel(y,2); syi = signe(yi); }
  else { yr = y; yi = NULL; syi = 0; }
  /* Compare absolute values of imaginary parts */
  if (!sxi)
  {
    if (syi && expo(yi) >= e) return -1;
    /* |Im x| ~ |Im y| ~ 0 */
  }
  else if (!syi)
  {
    if (sxi && expo(xi) >= e) return 1;
    /* |Im x| ~ |Im y| ~ 0 */
  }
  else
  {
    z = addrr_sign(xi, 1, yi, -1); sz = signe(z);
    if (sz && expo(z) >= e) return (int)sz;
  }
  /* |Im x| ~ |Im y|, sort according to real parts */
  z = subrr(xr, yr); sz = signe(z);
  if (sz && expo(z) >= e) return (int)sz;
  /* Re x ~ Re y. Place negative imaginary part before positive */
  return (int) (sxi - syi);
}

static GEN
clean_roots(GEN L, long l, long bit, long clean)
{
  long i, n = lg(L), ex = 5 - bit;
  GEN res = cgetg(n,t_COL);
  for (i=1; i<n; i++)
  {
    GEN c = gel(L,i);
    if (clean && isrealappr(c,ex))
    {
      if (typ(c) == t_COMPLEX) c = gel(c,1);
      c = mygprecrc(c, l, -bit);
    }
    else
      c = tocomplex(c, l, bit);
    gel(res,i) = c;
  }
  gen_sort_inplace(res, (void*)ex, &cmp_complex_appr, NULL);
  return res;
}

/* the vector of roots of p, with absolute error 2^(- prec2nbits(l)) */
static GEN
roots_aux(GEN p, long l, long clean)
{
  pari_sp av = avma;
  long bit;
  GEN L;

  if (typ(p) != t_POL)
  {
    if (gequal0(p)) pari_err_ROOTS0("roots");
    if (!isvalidcoeff(p)) pari_err_TYPE("roots",p);
    return cgetg(1,t_COL); /* constant polynomial */
  }
  if (!signe(p)) pari_err_ROOTS0("roots");
  checkvalidpol(p,"roots");
  if (lg(p) == 3) return cgetg(1,t_COL); /* constant polynomial */
  if (l < LOWDEFAULTPREC) l = LOWDEFAULTPREC;
  bit = prec2nbits(l);
  L = roots_com(p, bit);
  return gerepilecopy(av, clean_roots(L, l, bit, clean));
}
GEN
roots(GEN p, long l) { return roots_aux(p,l, 0); }
/* clean up roots. If root is real replace it by its real part */
GEN
cleanroots(GEN p, long l) { return roots_aux(p,l, 1); }

/* private variant of conjvec. Allow non rational coefficients, shallow
 * function. */
GEN
polmod_to_embed(GEN x, long prec)
{
  GEN v, T = gel(x,1), A = gel(x,2);
  long i, l;
  if (typ(A) != t_POL || varn(A) != varn(T))
  {
    checkvalidpol(T,"polmod_to_embed");
    return const_col(degpol(T), A);
  }
  v = cleanroots(T,prec); l = lg(v);
  for (i=1; i<l; i++) gel(v,i) = poleval(A,gel(v,i));
  return v;
}

GEN
QX_complex_roots(GEN p, long l)
{
  pari_sp av = avma;
  long bit, v;
  GEN L;

  if (!signe(p)) pari_err_ROOTS0("QX_complex_roots");
  if (lg(p) == 3) return cgetg(1,t_COL); /* constant polynomial */
  if (l < LOWDEFAULTPREC) l = LOWDEFAULTPREC;
  bit = prec2nbits(l);
  v = RgX_valrem(p, &p);
  L = lg(p) > 3? all_roots(Q_primpart(p), bit): cgetg(1,t_COL);
  if (v) L = shallowconcat(const_vec(v, real_0_bit(-bit)), L);
  return gerepilecopy(av, clean_roots(L, l, bit, 1));
}

/********************************************************************/
/**                                                                **/
/**                REAL ROOTS OF INTEGER POLYNOMIAL                **/
/**                                                                **/
/********************************************************************/

/* Count sign changes in the coefficients of (x+1)^deg(P)*P(1/(x+1)), P
 * has no rational root. The inversion is implicit (we take coefficients
 * backwards). */
static long
X2XP1(GEN P, GEN *Premapped)
{
  const pari_sp av = avma;
  GEN v = shallowcopy(P);
  long i, j, nb, s, dP = degpol(P), vlim = dP+2;

  for (j = 2; j < vlim; j++) gel(v, j+1) = addii(gel(v, j), gel(v, j+1));
  s = -signe(gel(v, vlim));
  vlim--; nb = 0;
  for (i = 1; i < dP; i++)
  {
    long s2 = -signe(gel(v, 2));
    int flag = (s2 == s);
    for (j = 2; j < vlim; j++)
    {
      gel(v, j+1) = addii(gel(v, j), gel(v, j+1));
      if (flag) flag = (s2 != signe(gel(v, j+1)));
    }
    if (s == signe(gel(v, vlim)))
    {
      if (++nb >= 2) return gc_long(av,2);
      s = -s;
    }
    /* if flag is set there will be no further sign changes */
    if (flag && (!Premapped || !nb)) return gc_long(av, nb);
    vlim--;
    if (gc_needed(av, 3))
    {
      if (DEBUGMEM>1) pari_warn(warnmem, "X2XP1, i = %ld/%ld", i, dP-1);
      if (!Premapped) setlg(v, vlim + 2);
      v = gerepilecopy(av, v);
    }
  }
  if (vlim >= 2 && s == signe(gel(v, vlim))) nb++;
  if (Premapped && nb == 1) *Premapped = v; else set_avma(av);
  return nb;
}

static long
_intervalcmp(GEN x, GEN y)
{
  if (typ(x) == t_VEC) x = gel(x, 1);
  if (typ(y) == t_VEC) y = gel(y, 1);
  return gcmp(x, y);
}

static GEN
_gen_nored(void *E, GEN x) { (void)E; return x; }
static GEN
_mp_add(void *E, GEN x, GEN y) { (void)E; return mpadd(x, y); }
static GEN
_mp_sub(void *E, GEN x, GEN y) { (void)E; return mpsub(x, y); }
static GEN
_mp_mul(void *E, GEN x, GEN y) { (void)E; return mpmul(x, y); }
static GEN
_mp_sqr(void *E, GEN x) { (void)E; return mpsqr(x); }
static GEN
_gen_one(void *E) { (void)E; return gen_1; }
static GEN
_gen_zero(void *E) { (void)E; return gen_0; }

static struct bb_algebra mp_algebra = { _gen_nored, _mp_add, _mp_sub,
                         _mp_mul, _mp_sqr, _gen_one, _gen_zero };

static GEN
_mp_cmul(void *E, GEN P, long a, GEN x) {(void)E; return mpmul(gel(P,a+2), x);}

/* Split the polynom P in two parts, whose coeffs have constant sign:
 * P(X) = X^D*Pp + Pm. Also compute the two parts of the derivative of P,
 * Pprimem = Pm', Pprimep = X*Pp'+ D*Pp => P' = X^(D-1)*Pprimep + Pprimem;
 * Pprimep[i] = (i+D) Pp[i]. Return D */
static long
split_pols(GEN P, GEN *pPp, GEN *pPm, GEN *pPprimep, GEN *pPprimem)
{
  long i, D, dP = degpol(P), s0 = signe(gel(P,2));
  GEN Pp, Pm, Pprimep, Pprimem;
  for(i=1; i <= dP; i++)
    if (signe(gel(P, i+2)) == -s0) break;
  D = i;
  Pm = cgetg(D + 2, t_POL);
  Pprimem = cgetg(D + 1, t_POL);
  Pp = cgetg(dP-D + 3, t_POL);
  Pprimep = cgetg(dP-D + 3, t_POL);
  Pm[1] = Pp[1] = Pprimem[1] = Pprimep[1] = P[1];
  for(i=0; i < D; i++)
  {
    GEN c = gel(P, i+2);
    gel(Pm, i+2) = c;
    if (i) gel(Pprimem, i+1) = mului(i, c);
  }
  for(; i <= dP; i++)
  {
    GEN c = gel(P, i+2);
    gel(Pp, i+2-D) = c;
    gel(Pprimep, i+2-D) = mului(i, c);
  }
  *pPm = normalizepol_lg(Pm, D+2);
  *pPprimem = normalizepol_lg(Pprimem, D+1);
  *pPp = normalizepol_lg(Pp, dP-D+3);
  *pPprimep = normalizepol_lg(Pprimep, dP-D+3);
  return dP - degpol(*pPp);
}

static GEN
bkeval_single_power(long d, GEN V)
{
  long mp = lg(V) - 2;
  if (d > mp) return gmul(gpowgs(gel(V, mp+1), d/mp), gel(V, (d%mp)+1));
  return gel(V, d+1);
}

static GEN
splitpoleval(GEN Pp, GEN Pm, GEN pows, long D, long bitprec)
{
  GEN vp = gen_bkeval_powers(Pp, degpol(Pp), pows, NULL, &mp_algebra, _mp_cmul);
  GEN vm = gen_bkeval_powers(Pm, degpol(Pm), pows, NULL, &mp_algebra, _mp_cmul);
  GEN xa = bkeval_single_power(D, pows);
  GEN r;
  if (!signe(vp)) return vm;
  vp = gmul(vp, xa);
  r = gadd(vp, vm);
  if (gexpo(vp) - (signe(r)? gexpo(r): 0) > prec2nbits(realprec(vp)) - bitprec)
    return NULL;
  return r;
}

/* optimized Cauchy bound for P = X^D*Pp + Pm, D > deg(Pm) */
static GEN
splitcauchy(GEN Pp, GEN Pm, long prec)
{
  GEN S = gel(Pp,2), A = gel(Pm,2);
  long i, lPm = lg(Pm), lPp = lg(Pp);
  for (i=3; i < lPm; i++) { GEN c = gel(Pm,i); if (abscmpii(A, c) < 0) A = c; }
  for (i=3; i < lPp; i++) S = addii(S, gel(Pp, i));
  return subsr(1, rdivii(A, S, prec)); /* 1 + |Pm|_oo / |Pp|_1 */
}

static GEN
ZX_deg1root(GEN P, long prec)
{
  GEN a = gel(P,3), b = gel(P,2);
  if (is_pm1(a))
  {
    b = itor(b, prec); if (signe(a) > 0) togglesign(b);
    return b;
  }
  return rdivii(negi(b), a, prec);
}

/* Newton for polynom P, P(0)!=0, with unique sign change => one root in ]0,oo[
 * P' has also at most one zero there */
static GEN
polsolve(GEN P, long bitprec)
{
  pari_sp av;
  GEN Pp, Pm, Pprimep, Pprimem, Pprime, Pprime2, ra, rb, rc, Pc;
  long dP = degpol(P), prec = nbits2prec(bitprec);
  long expoold, iter, D, rt, s0, bitaddprec, cprec, PREC;

  if (dP == 1) return ZX_deg1root(P, prec);
  Pprime = ZX_deriv(P);
  Pprime2 = ZX_deriv(Pprime);
  bitaddprec = 1 + 2*expu(dP); PREC = prec + nbits2prec(bitaddprec);
  D = split_pols(P, &Pp, &Pm, &Pprimep, &Pprimem); /* P = X^D*Pp + Pm */
  s0 = signe(gel(P, 2));
  rt = maxss(D, brent_kung_optpow(maxss(degpol(Pp), degpol(Pm)), 2, 1));
  rb = splitcauchy(Pp, Pm, DEFAULTPREC);
  for (cprec = DEFAULTPREC, expoold = LONG_MAX;;)
  {
    GEN pows = gen_powers(rb, rt, 1, NULL, _mp_sqr, _mp_mul, _gen_one);
    Pc = splitpoleval(Pp, Pm, pows, D, bitaddprec);
    if (!Pc) { cprec += EXTRAPREC64; rb = rtor(rb, cprec); continue; }
    if (signe(Pc) != s0) break;
    shiftr_inplace(rb,1);
  }
  for (iter = 0, ra = NULL;;)
  {
    GEN wdth;
    iter++;
    if (ra)
      rc = shiftr(addrr(ra, rb), -1);
    else
      rc = shiftr(rb, -1);
    for(;;)
    {
      GEN pows = gen_powers(rc, rt, 1, NULL, _mp_sqr, _mp_mul, _gen_one);
      Pc = splitpoleval(Pp, Pm, pows, D, bitaddprec+2);
      if (Pc) break;
      cprec += EXTRAPREC64;
      rc = rtor(rc, cprec);
    }
    if (signe(Pc) == s0)
      ra = rc;
    else
      rb = rc;
    if (!ra) continue;
    wdth = subrr(rb, ra);
    if (!(iter % 8))
    {
      GEN m1 = poleval(Pprime, ra), M2;
      if (signe(m1) == s0) continue;
      M2 = poleval(Pprime2, rb);
      if (abscmprr(gmul(M2, wdth), shiftr(m1, 1)) > 0) continue;
      break;
    }
    else if (gexpo(wdth) <= -bitprec)
      break;
  }
  rc = rb; av = avma;
  for(;; rc = gerepileuptoleaf(av, rc))
  {
    long exponew;
    GEN Ppc, dist, rcold = rc;
    GEN pows = gen_powers(rc, rt, 1, NULL, _mp_sqr, _mp_mul, _gen_one);
    Ppc = splitpoleval(Pprimep, Pprimem, pows, D-1, bitaddprec+4);
    if (Ppc) Pc = splitpoleval(Pp, Pm, pows, D, bitaddprec+4);
    if (!Ppc || !Pc)
    {
      if (cprec >= PREC)
        cprec += EXTRAPREC64;
      else
        cprec = minss(2*cprec, PREC);
      rc = rtor(rc, cprec); continue; /* backtrack one step */
    }
    dist = typ(Ppc) == t_REAL? divrr(Pc, Ppc): divri(Pc, Ppc);
    rc = subrr(rc, dist);
    if (cmprr(ra, rc) > 0 || cmprr(rb, rc) < 0)
    {
      if (cprec >= PREC) break;
      cprec = minss(2*cprec, PREC);
      rc = rtor(rcold, cprec); continue; /* backtrack one step */
    }
    if (expoold == LONG_MAX) { expoold = expo(dist); continue; }
    exponew = expo(dist);
    if (exponew < -bitprec - 1)
    {
      if (cprec >= PREC) break;
      cprec = minss(2*cprec, PREC);
      rc = rtor(rc, cprec); continue;
    }
    if (exponew > expoold - 2)
    {
      if (cprec >= PREC) break;
      expoold = LONG_MAX;
      cprec = minss(2*cprec, PREC);
      rc = rtor(rc, cprec); continue;
    }
    expoold = exponew;
  }
  return rtor(rc, prec);
}

/* Return primpart(P(x / 2)) */
static GEN
ZX_rescale2prim(GEN P)
{
  long i, l = lg(P), v, n;
  GEN Q;
  if (l==2) return pol_0(varn(P));
  Q = cgetg(l,t_POL); v = vali(gel(P,l-1));
  for (i = l-2, n = 1; v > n && i >= 2; i--, n++)
    v = minss(v, vali(gel(P,i)) + n);
  gel(Q,l-1) = v? shifti(gel(P,l-1), -v): gel(P,l-1);
  for (i = l-2, n = 1-v; i >= 2; i--, n++)
    gel(Q,i) = shifti(gel(P,i), n);
  Q[1] = P[1]; return Q;
}

/* assume Q0 has no rational root */
static GEN
usp(GEN Q0, long flag, long bitprec)
{
  const pari_sp av = avma;
  GEN Qremapped, Q, c, Lc, Lk, sol;
  GEN *pQremapped = flag == 1? &Qremapped: NULL;
  const long prec = nbits2prec(bitprec), deg = degpol(Q0);
  long listsize = 64, nbr = 0, nb_todo, ind, indf, i, k, nb;

  sol = zerocol(deg);
  Lc = zerovec(listsize);
  Lk = cgetg(listsize+1, t_VECSMALL);
  k = Lk[1] = 0;
  ind = 1; indf = 2;
  Q = Q0;
  c = gen_0;
  nb_todo = 1;
  while (nb_todo)
  {
    GEN nc = gel(Lc, ind);
    pari_sp av2;
    if (Lk[ind] == k + 1)
    {
      Q = Q0 = ZX_rescale2prim(Q0);
      c = gen_0;
    }
    if (!equalii(nc, c)) Q = ZX_translate(Q, subii(nc, c));
    av2 = avma;
    k = Lk[ind];
    ind++;
    c = nc;
    nb_todo--;
    nb = X2XP1(Q, pQremapped);

    if (nb == 1)
    { /* exactly one root */
      GEN s = gen_0;
      if (flag == 0)
      {
        s = mkvec2(gmul2n(c,-k), gmul2n(addiu(c,1),-k));
        s = gerepilecopy(av2, s);
      }
      else if (flag == 1) /* Caveat: Qremapped is the reciprocal polynomial */
      {
        s = polsolve(*pQremapped, bitprec+1);
        s = addir(c, divrr(s, addsr(1, s)));
        shiftr_inplace(s, -k);
        if (realprec(s) != prec) s = rtor(s, prec);
        s = gerepileupto(av2, s);
      }
      else set_avma(av2);
      gel(sol, ++nbr) = s;
    }
    else if (nb)
    { /* unknown, add two nodes to refine */
      if (indf + 2 > listsize)
      {
        if (ind>1)
        {
          for (i = ind; i < indf; i++)
          {
            gel(Lc, i-ind+1) = gel(Lc, i);
            Lk[i-ind+1] = Lk[i];
          }
          indf -= ind-1;
          ind = 1;
        }
        if (indf + 2 > listsize)
        {
          listsize *= 2;
          Lc = vec_lengthen(Lc, listsize);
          Lk = vecsmall_lengthen(Lk, listsize);
        }
        for (i = indf; i <= listsize; i++) gel(Lc, i) = gen_0;
      }
      gel(Lc, indf) = nc = shifti(c, 1);
      gel(Lc, indf + 1) = addiu(nc, 1);
      Lk[indf] = Lk[indf + 1] = k + 1;
      indf += 2;
      nb_todo += 2;
    }
    if (gc_needed(av, 2))
    {
      gerepileall(av, 6, &Q0, &Q, &c, &Lc, &Lk, &sol);
      if (DEBUGMEM > 1) pari_warn(warnmem, "ZX_Uspensky", avma);
    }
  }
  setlg(sol, nbr+1);
  return gerepilecopy(av, sol);
}

static GEN
ZX_Uspensky_equal_yes(GEN a, long flag, long bit)
{
  if (flag == 2) return gen_1;
  if (flag == 1 && typ(a) != t_REAL)
  {
    if (typ(a) == t_INT && !signe(a))
      a = real_0_bit(bit);
    else
      a = gtofp(a, nbits2prec(bit));
  }
  return mkcol(a);
}
static GEN
ZX_Uspensky_no(long flag)
{ return flag <= 1 ? cgetg(1, t_COL) : gen_0; }
/* ZX_Uspensky(P, [a,a], flag) */
static GEN
ZX_Uspensky_equal(GEN P, GEN a, long flag, long bit)
{
  if (typ(a) != t_INFINITY && gequal0(poleval(P, a)))
    return ZX_Uspensky_equal_yes(a, flag, bit);
  else
    return ZX_Uspensky_no(flag);
}
static int
sol_ok(GEN r, GEN a, GEN b) { return gcmp(a, r) <= 0 && gcmp(r, b) <= 0; }

/* P a ZX without real double roots; better if primitive and squarefree but
 * caller should ensure that. If flag & 4 assume that P has no rational root
 * (modest speedup) */
GEN
ZX_Uspensky(GEN P, GEN ab, long flag, long bitprec)
{
  pari_sp av = avma;
  GEN a, b, res, sol;
  double fb;
  long l, nbz, deg;

  if (ab)
  {
    if (typ(ab) == t_VEC)
    {
      if (lg(ab) != 3) pari_err_DIM("ZX_Uspensky");
      a = gel(ab, 1);
      b = gel(ab, 2);
    }
    else
    {
      a = ab;
      b = mkoo();
    }
  }
  else
  {
    a = mkmoo();
    b = mkoo();
  }
  if (flag & 4)
  {
    if (gcmp(a, b) >= 0) { set_avma(av); return ZX_Uspensky_no(flag); }
    flag &= ~4;
    sol = cgetg(1, t_COL);
  }
  else
  {
    switch (gcmp(a, b))
    {
      case 1: set_avma(av); return ZX_Uspensky_no(flag);
      case 0: return gerepilecopy(av, ZX_Uspensky_equal(P, a, flag, bitprec));
    }
    sol = nfrootsQ(P);
  }
  nbz = 0; l = lg(sol);
  if (l > 1)
  {
    long i, j;
    P = RgX_div(P, roots_to_pol(sol, varn(P)));
    if (!RgV_is_ZV(sol)) P = Q_primpart(P);
    for (i = j = 1; i < l; i++)
      if (sol_ok(gel(sol,i), a, b)) gel(sol,j++) = gel(sol,i);
    setlg(sol, j);
    if (flag == 2) { nbz = j-1; sol = utoi(nbz); }
    else if (flag == 1) sol = RgC_gtofp(sol, nbits2prec(bitprec));
  }
  else if (flag == 2) sol = gen_0;
  deg = degpol(P);
  if (deg == 0) return gerepilecopy(av, sol);
  if (typ(a) == t_INFINITY && typ(b) != t_INFINITY && gsigne(b))
  {
    fb = fujiwara_bound_real(P, -1);
    if (fb <= -pariINFINITY) a = gen_0;
    else if (fb < 0) a = gen_m1;
    else a = negi(int2n((long)ceil(fb)));
  }
  if (typ(b) == t_INFINITY && typ(a) != t_INFINITY && gsigne(a))
  {
    fb = fujiwara_bound_real(P, 1);
    if (fb <= -pariINFINITY) b = gen_0;
    else if (fb < 0) b = gen_1;
    else b = int2n((long)ceil(fb));
  }
  if (typ(a) != t_INFINITY && typ(b) != t_INFINITY)
  {
    GEN d, ad, bd, diff;
    long i;
    /* can occur if one of a,b was initially a t_INFINITY */
    if (gequal(a,b)) return gerepilecopy(av, sol);
    d = lcmii(Q_denom(a), Q_denom(b));
    if (is_pm1(d)) { d = NULL; ad = a; bd = b; }
    else
    { P = ZX_rescale(P, d); ad = gmul(a, d); bd = gmul(b, d); }
    diff = subii(bd, ad);
    P = ZX_affine(P, diff, ad);
    res = usp(P, flag, bitprec);
    if (flag <= 1)
    {
      for (i = 1; i < lg(res); i++)
      {
        GEN z = gmul(diff, gel(res, i));
        if (typ(z) == t_VEC)
        {
          gel(z, 1) = gadd(ad, gel(z, 1));
          gel(z, 2) = gadd(ad, gel(z, 2));
        }
        else
          z = gadd(ad, z);
        if (d) z = gdiv(z, d);
        gel(res, i) = z;
      }
      sol = shallowconcat(sol, res);
    }
    else
      nbz += lg(res) - 1;
  }
  if (typ(b) == t_INFINITY && (fb=fujiwara_bound_real(P, 1)) > -pariINFINITY)
  {
    long bp = maxss((long)ceil(fb), 0);
    res = usp(ZX_unscale2n(P, bp), flag, bitprec);
    if (flag <= 1)
      sol = shallowconcat(sol, gmul2n(res, bp));
    else
      nbz += lg(res)-1;
  }
  if (typ(a) == t_INFINITY && (fb=fujiwara_bound_real(P,-1)) > -pariINFINITY)
  {
    long i, bm = maxss((long)ceil(fb), 0);
    res = usp(ZX_unscale2n(ZX_z_unscale(P, -1), bm), flag, bitprec);
    if (flag <= 1)
    {
      for (i = 1; i < lg(res); i++)
      {
        GEN z = gneg(gmul2n(gel(res, i), bm));
        if (typ(z) == t_VEC) swap(gel(z, 1), gel(z, 2));
        gel(res, i) = z;
      }
      sol = shallowconcat(res, sol);
    }
    else
      nbz += lg(res)-1;
  }
  if (flag >= 2) return utoi(nbz);
  if (flag)
    sol = sort(sol);
  else
    sol = gen_sort(sol, (void *)_intervalcmp, cmp_nodata);
  return gerepileupto(av, sol);
}

/* x a scalar */
static GEN
rootsdeg0(GEN x)
{
  if (!is_real_t(typ(x))) pari_err_TYPE("realroots",x);
  if (gequal0(x)) pari_err_ROOTS0("realroots");
  return cgetg(1,t_COL); /* constant polynomial */
}
static void
checkbound(GEN a)
{
  switch(typ(a))
  {
    case t_INT: case t_FRAC: case t_INFINITY: break;
    default: pari_err_TYPE("polrealroots", a);
  }
}
static GEN
check_ab(GEN ab)
{
  GEN a, b;
  if (!ab) return NULL;
  if (typ(ab) != t_VEC || lg(ab) != 3) pari_err_TYPE("polrootsreal",ab);
  a = gel(ab,1); checkbound(a);
  b = gel(ab,2); checkbound(b);
  if (typ(a) == t_INFINITY && inf_get_sign(a) < 0 &&
      typ(b) == t_INFINITY && inf_get_sign(b) > 0) ab = NULL;
  return ab;
}
/* e^(1/h) assuming the h-th root is real, beware that sqrtnr assumes e >= 0 */
static GEN
_sqrtnr(GEN e, long h)
{
  long s;
  GEN r;
  if (h == 2) return sqrtr(e);
  s = signe(e); setsigne(e, 1); /* e < 0 is possible, implies h is odd */
  r = sqrtnr(e, h); if (s < 0) setsigne(r, -1);
  return r;
}
GEN
realroots(GEN P, GEN ab, long prec)
{
  pari_sp av = avma;
  GEN sol = NULL, fa, ex;
  long i, j, v, l;

  ab = check_ab(ab);
  if (typ(P) != t_POL) return rootsdeg0(P);
  switch(degpol(P))
  {
    case -1: return rootsdeg0(gen_0);
    case 0: return rootsdeg0(gel(P,2));
  }
  if (!RgX_is_ZX(P)) P = RgX_rescale_to_int(P);
  v = ZX_valrem(Q_primpart(P), &P);
  fa = ZX_squff(P, &ex); l = lg(fa); sol = cgetg(l + 1, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN Pi = gel(fa, i), soli, soli2;
    long n, h;
    if (ab) h = 1; else Pi = ZX_deflate_max(Pi, &h);
    soli = ZX_Uspensky(Pi, odd(h)? ab: gen_0, 1, prec2nbits(prec));
    n = lg(soli); soli2 = odd(h)? NULL: cgetg(n, t_COL);
    for (j = 1; j < n; j++)
    {
      GEN r = gel(soli, j); /* != 0 */
      if (typ(r) != t_REAL) gel(soli, j) = r = gtofp(r, prec);
      if (h > 1)
      {
        gel(soli, j) = r = _sqrtnr(r, h);
        if (soli2) gel(soli2, j) = negr(r);
      }
    }
    if (soli2) soli = shallowconcat(soli, soli2);
    if (ex[i] > 1) soli = shallowconcat1( const_vec(ex[i], soli) );
    gel(sol, i) = soli;
  }
  if (v && (!ab || (gsigne(gel(ab,1)) <= 0 && gsigne(gel(ab,2)) >= 0)))
    gel(sol, i++) = const_col(v, real_0(prec));
  setlg(sol, i); if (i == 1) { set_avma(av); return cgetg(1,t_COL); }
  return gerepileupto(av, sort(shallowconcat1(sol)));
}
GEN
ZX_realroots_irred(GEN P, long prec)
{
  long dP = degpol(P), j, n, h;
  GEN sol, sol2;
  pari_sp av;
  if (dP == 1) retmkvec(ZX_deg1root(P, prec));
  av = avma; P = ZX_deflate_max(P, &h);
  if (h == dP)
  {
    GEN r = _sqrtnr(ZX_deg1root(P, prec), h);
    return gerepilecopy(av, odd(h)? mkvec(r): mkvec2(negr(r), r));
  }
  sol = ZX_Uspensky(P, odd(h)? NULL: gen_0, 1 | 4, prec2nbits(prec));
  n = lg(sol); sol2 = odd(h)? NULL: cgetg(n, t_COL);
  for (j = 1; j < n; j++)
  {
    GEN r = gel(sol, j);
    if (typ(r) != t_REAL) gel(sol, j) = r = gtofp(r, prec);
    if (h > 1)
    {
      gel(sol, j) = r = _sqrtnr(r, h);
      if (sol2) gel(sol2, j) = negr(r);
    }
  }
  if (sol2) sol = shallowconcat(sol, sol2);
  return gerepileupto(av, sort(sol));
}

static long
ZX_sturm_i(GEN P, long flag)
{
  pari_sp av;
  long h, r, dP = degpol(P);
  if (dP == 1) return 1;
  av = avma; P = ZX_deflate_max(P, &h);
  if (h == dP)
  { /* now deg P = 1 */
    if (odd(h))
      r = 1;
    else
      r = (signe(gel(P,2)) != signe(gel(P,3)))? 2: 0;
    return gc_long(av, r);
  }
  if (odd(h))
    r = itou(ZX_Uspensky(P, NULL, flag, 0));
  else
    r = 2*itou(ZX_Uspensky(P, gen_0, flag, 0));
  return gc_long(av,r);
}
/* P nonconstant, squarefree ZX */
long
ZX_sturmpart(GEN P, GEN ab)
{
  pari_sp av = avma;
  if (!check_ab(ab)) return ZX_sturm(P);
  return gc_long(av, itou(ZX_Uspensky(P, ab, 2, 0)));
}
/* P nonconstant, squarefree ZX */
long
ZX_sturm(GEN P) { return ZX_sturm_i(P, 2); }
/* P irreducible ZX */
long
ZX_sturm_irred(GEN P) { return ZX_sturm_i(P, 2 + 4); }
