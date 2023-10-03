/* Copyright (C) 2021  The PARI group.

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

/***********************************************************************/
/**                 LAMBERT's W_K FUNCTIONS                           **/
/***********************************************************************/
/* roughly follows Veberic, https://arxiv.org/abs/1003.1628 */

static double
dblL1L2(double L1)
{
  double L2 = log(-L1), LI = 1 / L1, N2, N3, N4, N5;
  N2 = (L2-2.)/2.; N3 = (6.+L2*(-9.+2.*L2))/6.;
  N4 = (-12.+L2*(36.+L2*(-22.+3*L2)))/12.;
  N5 = (60.+L2*(-300.+L2*(350.+L2*(-125.+12*L2))))/60.;
  return L1-L2+L2*LI*(1+LI*(N2+LI*(N3+LI*(N4+LI*N5))));
}

/* rough approximation to W0(a > -1/e), < 1% relative error */
double
dbllambertW0(double a)
{
  if (a < -0.2583)
  {
    const double c2 = -1./3, c3 = 11./72, c4 = -43./540, c5 = 769./17280;
    double p = sqrt(2 * (M_E * a + 1));
    if (a < -0.3243) return -1+p*(1+p*(c2+p*c3));
    return -1+p*(1+p*(c2+p*(c3+p*(c4+p*c5))));
  }
  else
  {
    double Wd = log(1.+a);
    Wd *= (1.-log(Wd/a))/(1.+Wd);
    if (a < 0.6482 && a > -0.1838) return Wd;
    return Wd*(1.-log(Wd/a))/(1.+Wd);
  }
}
/* uniform approximation to W0, at least 15 bits. */
static double
dbllambertW0init(double a)
{
  if (a < -0.323581)
  {
    const double c2 = 1./3., c3 = 11./72., c4 = 43./540., c5 = 769./17280.;
    const double c6 = 221./8505., c7 = 680863./43545600.;
    const double c8 = 1963./204120., c9 = 226287557./37623398400.;
    double p = M_E * a + 1;
    if (p <= 0) return -1;
    p = -sqrt(2 * p);
    return -(1.+p*(1.+p*(c2+p*(c3+p*(c4+p*(c5+p*(c6+p*(c7+p*(c8+p*c9)))))))));
  }
  if (a < 0.145469)
  {
    const double a1 = 5.931375, a2 = 11.392205, a3 = 7.338883, a4 = 0.653449;
    const double b1 = 6.931373, b2 = 16.823494, b3 = 16.430723, b4 = 5.115235;
    double n = 1.+a*(a1+a*(a2+a*(a3+a*a4)));
    double d = 1.+a*(b1+a*(b2+a*(b3+a*b4)));
    return a * n / d;
  }
  if (a < 8.706658)
  {
    const double a1 = 2.445053, a2 = 1.343664, a3 = 0.148440, a4 = 0.000804;
    const double b1 = 3.444708, b2 = 3.292489, b3 = 0.916460, b4 = 0.053068;
    double n = 1.+a*(a1+a*(a2+a*(a3+a*a4)));
    double d = 1.+a*(b1+a*(b2+a*(b3+a*b4)));
    return a * n / d;
  }
  else
  {
    double w = log(1.+a);
    w *= (1.-log(w/a)) / (1.+w);
    return w * (1.-log(w/a)) / (1.+w);
  }
}

/* rough approximation to W_{-1}(0 > a > -1/e), < 1% relative error */
double
dbllambertW_1(double a)
{
  if (a < -0.2464)
  {
    const double c2 = -1./3, c3 = 11./72, c4 = -43./540, c5 = 769./17280;
    double p = -sqrt(2 * (M_E * a + 1));
    if (a < -0.3243) return -1+p*(1+p*(c2+p*c3));
    return -1+p*(1+p*(c2+p*(c3+p*(c4+p*c5))));
  }
  else
  {
    double Wd;
    a = -a; Wd = -log(a);
    Wd *= (1.-log(Wd/a))/(1.-Wd);
    if (a < 0.0056) return -Wd;
    return -Wd*(1.-log(Wd/a))/(1.-Wd);
  }
}
/* uniform approximation to W_{-1}, at least 15 bits. */
static double
dbllambertW_1init(double a)
{
  if (a < -0.302985)
  {
    const double c2 = 1./3., c3 = 11./72., c4 = 43./540., c5 = 769./17280.;
    const double c6 = 221./8505., c7 = 680863./43545600.;
    const double c8 = 1963./204120., c9 = 226287557./37623398400.;
    double p = M_E * a + 1;
    if (p <= 0) return -1;
    p = sqrt(2 * p);
    return -(1.+p*(1.+p*(c2+p*(c3+p*(c4+p*(c5+p*(c6+p*(c7+p*(c8+p*c9)))))))));
  }
  if (a <= -0.051012)
  {
    const double a0 = -7.814176, a1 = 253.888101, a2 = 657.949317;
    const double b1 = -60.439587, b2 = 99.985670, b3 = 682.607399;
    const double b4 = 962.178439, b5 = 1477.934128;
    double n = a0+a*(a1+a*a2);
    double d = 1+a*(b1+a*(b2+a*(b3+a*(b4+a*b5))));
    return n / d;
  }
  return dblL1L2(log(-a));
}

/* uniform approximation to more than 46 bits, 50 bits away from -1/e;
 * branch = -1 or 0 */
static double
dbllambertWfritsch(GEN ga, int branch)
{
  double a, z, w1, q, w;
  if (expo(ga) >= 0x3fe)
  { /* branch = 0 */
    double w = dbllog2(ga) * M_LN2; /* ~ log(1+a) ~ log a */
    return w * (1.+w-log(w)) / (1.+w);
  }
  a = rtodbl(ga);
  w = branch? dbllambertW_1init(a): dbllambertW0init(a);
  if (w == -1.|| w == 0.) return w;
  z = log(a / w) - w; w1 = 1. + w;
  q = 2. * w1 * (w1 + (2./3.) * z);
  return w * (1 + (z / w1) * (q - z) / (q - 2 * z));
}

static double
dbllambertWhalleyspec(double loga)
{
  double w = dblL1L2(loga);
  for(;;)
  {
    double n = w + log(-w) - loga, d = 1 - w, r = n / (d + n / d);
    w *= 1 - r; if (r < 2.e-15) return w;
  }
}
/* k = 0 or -1. */
static GEN
lambertW(GEN z, long k, long prec)
{
  pari_sp av = avma;
  long bit = prec2nbits(prec), L = -(bit / 3 + 10), ct = 0, p, pb;
  double wd;
  GEN w, vp;

  if (gequal0(z) && !k) return real_0(prec);
  z = gtofp(z, prec);
  if (k == -1)
  {
    long e = expo(z);
    if (signe(z) >= 0) pari_err_DOMAIN("lambertw", "z", ">", gen_0, z);
    wd = e < -512? dbllambertWhalleyspec(dbllog2(z) * M_LN2)
                 : dbllambertWfritsch(z, -1);
  }
  else
    wd = dbllambertWfritsch(z, 0);
  if (fabs(wd + 1) < 1e-5)
  {
    long prec2 = prec + EXTRAPREC64;
    GEN Z = rtor(z, prec2);
    GEN t = addrs(mulrr(Z, gexp(gen_1, prec2)), 1);
    if (signe(t) <= 0) { set_avma(av); return real_m1(prec); }
    if (realprec(t) < prec)
    {
      prec2 += prec - realprec(t);
      Z = rtor(z, prec2);
      t = addrs(mulrr(Z, gexp(gen_1, prec2)), 1);
    }
    t = sqrtr(shiftr(t, 1));
    w = gprec_w(k == -1? subsr(-1, t) : subrs(t, 1), prec);
    p = prec - 2; vp = NULL;
  }
  else
  { /* away from -1/e: can reduce accuracy and self-correct */
    w = wd == 0.? z: dbltor(wd);
    vp = cgetg(30, t_VECSMALL); ct = 0; pb = bit;
    while (pb > BITS_IN_LONG * 3/4)
    { vp[++ct] = (pb + BITS_IN_LONG-1) >> TWOPOTBITS_IN_LONG; pb = (pb + 2) / 3; }
    p = vp[ct]; w = gprec_w(w, p + 2);
  }
  if ((k == -1 && (bit < 192 || bit > 640)) || (k == 0 && bit > 1024))
  {
    for(;;)
    {
      GEN t, ew, n, d;
      ew = mplog(divrr(w, z)); n = addrr(w, ew); d = addrs(w, 1);
      t = divrr(n, shiftr(d, 1));
      w = mulrr(w, subsr(1, divrr(n, addrr(d, t))));
      if (p >= prec-2 && expo(n) - expo(d) - expo(w) <= L) break;
      if (vp) { if (--ct) p = vp[ct]; w = gprec_w(w, ct? p + 2: prec); }
    }
  }
  else
  {
    for(;;)
    {
      GEN t, ew, wew, n, d;
      ew = mpexp(w); wew = mulrr(w, ew); n = subrr(wew, z); d = addrr(ew, wew);
      t = divrr(mulrr(addrs(w, 2), n), shiftr(addrs(w, 1), 1));
      w = subrr(w, divrr(n, subrr(d, t)));
      if (p >= prec-2 && expo(n) - expo(d) - expo(w) <= L) break;
      if (vp) { if (--ct) p = vp[ct]; w = gprec_w(w, ct? p + 2: prec); }
    }
  }
  return gerepileupto(av, w);
}

/*********************************************************************/
/*                       Complex branches                            */
/*********************************************************************/

/* x *= (1 - (x + log(x) - L) / (x + 1)); L = log(z) + 2IPi * k */
static GEN
lamaux(GEN x, GEN L, long *pe, long prec)
{
  GEN n = gsub(gadd(x, glog(x, prec)), L);
  if (pe) *pe = maxss(4, -gexpo(n));
  if (gequal0(imag_i(n))) n = real_i(n);
  return gmul(x, gsubsg(1, gdiv(n, gaddsg(1, x))));
}

/* Complex branches, experimental */
static GEN
lambertWC(GEN z, long branch, long prec)
{
  pari_sp av = avma;
  GEN w, pii2k, zl, lzl, L, Lz;
  long bit0, si, j, fl = 0, lim = 6, lp = DEFAULTPREC, bit = prec2nbits(prec);

  si = gsigne(imag_i(z)); if (!si) z = real_i(z);
  pii2k = gmulsg(branch, PiI2(lp));
  zl = gtofp(z, lp); lzl = glog(zl, lp);
  /* From here */
  if (branch == 0 || branch * si < 0
      || (si == 0 && gsigne(z) < 0 && branch == -1))
  {
    GEN lnzl1 = gaddsg(1, glog(gneg(zl), lp));
    if (si == 0) si = gsigne(lnzl1);
    if ((branch == 0 || branch * si < 0) && gexpo(lnzl1) < -1)
    { /* close to -1/e */
      w = gaddsg(1, gmul(z, gexp(gen_1, prec)));
      w = gprec_wtrunc(w, lp);
      w = gsqrt(gmul2n(w, 1), lp);
      w = branch * si < 0? gsubsg(-1, w): gaddsg(-1, w);
      lim = 10; fl = 1;
    }
    if (branch == 0 && !fl && gexpo(lzl) < 0) { w = zl; fl = 1; }
  }
  if (!fl)
  {
    if (branch)
    {
      GEN lr = glog(pii2k, lp);
      w = gadd(gsub(gadd(pii2k, lzl), lr), gdiv(gsub(lr, lzl), pii2k));
    }
    else
    {
      GEN p = gaddsg(1, gmul(z, gexp(gen_1, lp)));
      w = gexpo(p) > 0? lzl: gaddgs(gsqrt(p, lp), -1);
    }
  }
  /* to here: heuristic */
  L = gadd(lzl, pii2k);
  for (j = 1; j < lim; j++) w = lamaux(w, L, NULL, lp);
  Lz = NULL;
  if (branch == 0 || branch == -1)
  {
    Lz = glog(z, prec);
    if (branch == -1)
    {
      long flag = 1;
      if (!si && signe(z) <= 0 && signe(addrs(Lz, 1))) flag = 0;
      if (flag) Lz = gsub(Lz, PiI2(prec));
    }
  }
  w = lamaux(w, L, &bit0, lp);
  while (bit0 < bit || (Lz && gexpo(gsub(gadd(w, glog(w, prec)), Lz)) > 16-bit))
  {
    long p = nbits2prec(bit0 <<= 1);
    L = gadd(gmulsg(branch, PiI2(p)), glog(gprec_w(z, p), p));
    w = lamaux(gprec_w(w, p), L, NULL, p);
  }
  return gerepilecopy(av, gprec_w(w, nbits2prec(bit)));
}

/* exp(t (1 + O(t^n))), n >= 0 */
static GEN
serexp0(long v, long n)
{
  GEN y = cgetg(n+3, t_SER), t;
  long i;
  y[1] = evalsigne(1) | evalvarn(v) | evalvalser(0);
  gel(y,2) = gen_1; if (!n) return y;
  gel(y,3) = gen_1; if (n == 1) return y;
  for (i=2, t = gen_2; i < n; i++, t = muliu(t,i)) gel(y,i+2) = mkfrac(gen_1,t);
  gel(y,i+2) = mkfrac(gen_1,t); return y;
}

/* series expansion of W at -1/e */
static GEN
Wbra(long N)
{
  GEN v = cgetg(N + 2, t_VEC);
  long n;
  gel(v, 1) = gen_m1;
  gel(v, 2) = gen_1;
  for (n = 2; n <= N; n++)
  {
    GEN t = gel(v,n), a = gen_0;
    long k, K = (n - 1) >> 1;
    for (k = 1; k <= K; k++) t = gadd(t, gmul2n(gel(v,n-2*k), -k));
    for (k = 2; k < n; k++) a = gadd(a, gmul(gel(v,k+1), gel(v,n+2-k)));
    gel(v,n+1) = gsub(gdivgs(t, -n-1), gmul2n(a, -1));
  }
  return RgV_to_RgX(v, 0);
}

static GEN
reverse(GEN y)
{
  GEN z = ser2rfrac_i(y);
  long l = lg(z);
  return RgX_to_ser(RgXn_reverse(z, l-2), l-1);
}
static GEN
serlambertW(GEN y, long branch, long prec)
{
  long n, vy, val, v;
  GEN t = NULL;

  if (!signe(y)) return gcopy(y);
  v = valser(y);
  if (v < 0) pari_err_DOMAIN("lambertw","valuation", "<", gen_0, y);
  if (v > 0 && branch)
    pari_err_DOMAIN("lambertw [k != 0]", "x", "~", gen_0, y);
  vy = varn(y); n = lg(y)-3;
  for (val = 1; val < n; val++)
    if (!gequal0(polcoef_i(y, val, vy))) break;
  if (v)
  {
    t = serexp0(vy, n / val);
    setvalser(t, 1); t = reverse(t); /* rev(x exp(x)) */
  }
  else
  {
    GEN y0 = gel(y,2), x = glambertW(y0, branch, prec);
    if (val > n) return scalarser(x, vy, n+1);
    y = serchop0(y);
    if (gequalm1(x))
    { /* y0 ~ -1/e, branch = 0 or -1 */
      GEN p = gmul(shiftr(gexp(gen_1,prec), 1), y);
      if (odd(val)) pari_err(e_MISC, "odd valuation at branch point");
      p = gsqrt(p, prec); if (odd(branch)) p = gneg(p);
      n -= val >> 1;
      t = RgXn_eval(Wbra(n), ser2rfrac_i(p), n);
      return gtoser(t, varn(t), lg(p));
    }
    t = serexp0(vy, n / val);
    /* (x + t) exp(x + t) = (y0 + t y0/x) * exp(t) */
    t = gmul(deg1pol_shallow(gdiv(y0,x), y0, vy), t);
    t = gadd(x, reverse(serchop0(t)));
  }
  return normalizeser(gsubst(t, vy, y));
}

static GEN
lambertp(GEN x)
{
  pari_sp av = avma;
  long k;
  GEN y;

  if (gequal0(x)) return gcopy(x);
  if (!valp(x)) { x = leafcopy(x); setvalp(x, 1); }
  k = Qp_exp_prec(x);
  if (k < 0) return NULL;
  y = gpowgs(cvstop2(k, x), k - 1);
  for (k--; k; k--)
    y = gsub(gpowgs(cvstop2(k, x), k - 1), gdivgu(gmul(x, y), k + 1));
  return gerepileupto(av, gmul(x, y));
}

/* y a t_REAL */
static int
useC(GEN y, long k)
{
  if (signe(y) > 0 || (k && k != -1)) return k ? 1: 0;
  return gsigne(addsr(1, logr_abs(y))) > 0;
}
static GEN
glambertW_i(void *E, GEN y, long prec)
{
  pari_sp av;
  long k = (long)E, p;
  GEN z;
  if (gequal0(y))
  {
    if (k) pari_err_DOMAIN("glambertW","argument","",gen_0,y);
    return gcopy(y);
  }
  switch(typ(y))
  {
    case t_REAL:
      p = minss(prec, realprec(y));
      return useC(y, k)? lambertWC(y, k, p): lambertW(y, k, p);
    case t_PADIC: z = lambertp(y);
      if (!z) pari_err_DOMAIN("glambertW(t_PADIC)","argument","",gen_0,y);
      return z;
    case t_COMPLEX:
      p = precision(y);
      return lambertWC(y, k, p? p: prec);
    default:
      av = avma; if (!(z = toser_i(y))) break;
      return gerepileupto(av, serlambertW(z, k, prec));
  }
  return trans_evalgen("lambert", E, glambertW_i, y, prec);
}

GEN
glambertW(GEN y, long k, long prec)
{
  return glambertW_i((void*)k, y, prec);
}
GEN
mplambertW(GEN y, long prec) { return lambertW(y, 0, prec); }

/*********************************************************************/
/*                        Application                                */
/*********************************************************************/
/* Solve x - a * log(x) = b with a > 0 and b >= a * (1 - log(a)). */
GEN
mplambertx_logx(GEN a, GEN b, long bit)
{
  pari_sp av = avma;
  GEN e = gexp(gneg(gdiv(b, a)), nbits2prec(bit));
  return gerepileupto(av, gmul(gneg(a), lambertW(gneg(gdiv(e, a)), -1, bit)));
}
/* Special case a = 1, b = log(y): solve e^x / x = y with y >= exp(1). */
GEN
mplambertX(GEN y, long bit)
{
  pari_sp av = avma;
  return gerepileupto(av, gneg(lambertW(gneg(ginv(y)), -1, bit)));
}

/* Solve x * log(x) - a * x = b; if b < 0, assume a >= 1 + log |b|. */
GEN
mplambertxlogx_x(GEN a, GEN b, long bit)
{
  pari_sp av = avma;
  long s = gsigne(b);
  GEN e;
  if (!s) return gen_0;
  e = gexp(gneg(a), nbits2prec(bit));
  return gerepileupto(av, gdiv(b, lambertW(gmul(b, e), s > 0? 0: -1, bit)));
}
