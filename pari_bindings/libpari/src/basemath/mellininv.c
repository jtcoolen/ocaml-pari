/* Copyright (C) 2015  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_gammamellininv

/*******************************************************************/
/*               Computation of inverse Mellin                     */
/*               transforms of gamma products.                     */
/*******************************************************************/
/* Handle complex Vga whose sum is real */
static GEN
sumVga(GEN Vga) { return real_i(vecsum(Vga)); }

/* ac != 0 */
static double
lemma526_i(double ac, double c, double t, double B)
{
  double D = -B/ac; /* sgn(t) = sgn(a) = - sgn(D) */
  if (D <= 0)
  {
    if (D > -100)
    {
      D = -exp(D) / t;
      if (D < - 1/M_E) return 0;
      D = dbllambertW_1(D);
    }
    else
    { /* avoid underflow, use asymptotic expansion */
      double U = D - log(t);
      D = U - log(-U);
    }
    return pow(maxdd(t, -t * D), c);
  }
  else
  {
    if (D < 100)
      D = dbllambertW0(-exp(D) / t);
    else
    { /* avoid overflow, use asymptotic expansion */
      double U = D - log(-t);
      D = U - log(U);
    }
    return pow(-t * D, c);
  }
}
/* b > 0, c > 0; solve x^a exp(-b x^(1/c)) < e^(-B) for x >= 0 */
double
dbllemma526(double a, double b, double c, double B)
{
  double ac;
  if (!a) return B <= 0? 0: pow(B/b, c);
  ac = a*c; if (B < 0) B = 1e-9;
  return lemma526_i(ac, c, ac/b, B);
}
/* Same, special case b/c = 2Pi, the only one needed: for c = d/2 */
double
dblcoro526(double a, double c, double B)
{
  if (!a) return B <= 0? 0: pow(B/(2*M_PI*c), c);
  if (B < 0) B = 1e-9;
  return lemma526_i(a*c, c, a/(2*M_PI), B);
}

static const double MELLININV_CUTOFF = 121.; /* C*C */

/* x real */
static GEN
RMOD2(GEN x) { return gsub(x, gmul2n(gdiventgs(x,2), 1)); }
/* x real or complex, return canonical representative for x mod 2Z */
static GEN
MOD2(GEN x)
{ return typ(x) == t_COMPLEX? mkcomplex(RMOD2(gel(x,1)), gel(x,2)): RMOD2(x); }
static GEN
RgV_MOD2(GEN x)
{ pari_APPLY_same(MOD2(gel(x,i))); }

/* classes of poles of the gamma factor mod 2Z, sorted by increasing
 * Re(s) mod 2 (in [0,2[).*/
static GEN
gammapoles(GEN Vga, long *pdV, long bit)
{
  long i, m, emax, l = lg(Vga);
  GEN P, B = RgV_MOD2(Vga), V = cgetg(l, t_VEC);
  P = gen_indexsort(B, (void*)lexcmp, cmp_nodata);
  for (i = m = 1; i < l;)
  {
    GEN u = gel(B, P[i]);
    long k;
    for(k = i+1; k < l; k++)
    {
      GEN v = gsub(u, gel(B, P[k]));
      if (!gequal0(v) && (!isinexactreal(v) || gexpo(v) > -bit)) break;
    }
    gel(V, m++) = vecslice(P,i,k-1);
    i = k;
  }
  setlg(V, m); emax = 0;
  for (i = 1; i < m; i++)
  {
    long j, e = 0, li = lg(gel(V,i))-1;
    GEN b = gel(B, gel(V,i)[1]);
    for (j = 1; j < m; j++)
      if (j != i) e -= gexpo(gsub(gel(B, gel(V,j)[1]), b));
    emax = maxss(emax, e * li);
  }
  for (i = 1; i < m; i++) gel(V,i) = vecpermute(Vga, gel(V,i));
  *pdV = emax; return V;
}

static GEN
sercoeff(GEN x, long n, long prec)
{
  long N = n - valser(x);
  return (N < 0)? gen_0: gprec_wtrunc(gel(x, N+2), prec);
}

/* prod_i Gamma(s/2 + (m+LA[i])/2), set t *= prod_i (s/2 + (m+LA[i])/2) */
static GEN
get_gamma(GEN *pt, GEN LA, GEN m, int round, long precdl, long prec)
{
  long i, l = lg(LA);
  GEN pr = NULL, t = *pt;
  for (i = 1; i < l; i++)
  {
    GEN u, g, a = gmul2n(gadd(m, gel(LA,i)), -1);
    if (round) a = ground(a);
    u = deg1pol_shallow(ghalf, a, 0);
    g = ggamma(RgX_to_ser(u, precdl), prec);
    pr = pr? gmul(pr, g): g;
    t = t? gmul(t, u): u;
  }
  *pt = t; return pr;
}
/* generalized power series expansion of inverse Mellin around x = 0;
 * m-th derivative */
static GEN
Kderivsmallinit(GEN ldata, GEN Vga, long m, long bit)
{
  const double C2 = MELLININV_CUTOFF;
  long prec2, N, j, l, dLA, limn, d = lg(Vga)-1;
  GEN piA, LA, L, M, mat;

  LA = gammapoles(Vga, &dLA, bit); N = lg(LA)-1;
  prec2 = nbits2prec(dLA + bit * (1 + M_PI*d/C2));
#if BITS_IN_LONG == 32
  prec2 += prec2 & 1L;
#endif
  if (ldata) Vga = ldata_get_gammavec(ldata_newprec(ldata, prec2));
  L = cgetg(N+1, t_VECSMALL);
  M = cgetg(N+1, t_VEC);
  mat = cgetg(N+1, t_VEC);
  limn = ceil(2*M_LN2*bit / (d * dbllambertW0(C2/(M_PI*M_E*d))));
  l = limn + 2;
  for (j = 1; j <= N; j++)
  {
    GEN S = gel(LA,j);
    GEN C, c, mj, G = NULL, t = NULL, tj = NULL;
    long i, k, n, jj, lj = L[j] = lg(S)-1, precdl = lj+3;

    gel(M,j) = mj = gsubsg(2, gel(S, vecindexmin(real_i(S))));
    for (jj = 1; jj <= N; jj++)
    {
      GEN g;
      if (jj == j) /* poles come from this class only */
        g = get_gamma(&tj, gel(LA,jj), mj, 1, precdl, prec2);
      else
        g = get_gamma(&t, gel(LA,jj), mj, 0, precdl, prec2);
      G = G? gmul(G, g): g;
    }
    c = cgetg(limn+2,t_COL); gel(c,1) = G;
    for (n=1; n <= limn; n++)
    {
      GEN A = utoineg(2*n), T = RgX_translate(tj, A);
      /* T = exact polynomial, may vanish at 0 (=> pole in c[n+1]) */
      if (t) T = RgX_mul(T, RgX_translate(t, A)); /* no pole here */
      gel(c,n+1) = gdiv(gel(c,n), T);
    }
    gel(mat, j) = C = cgetg(lj+1, t_COL);
    for (k = 1; k <= lj; k++)
    {
      GEN L = cgetg(l, t_POL);
      for (n = 2; n < l; n++) gel(L,n) = sercoeff(gel(c,n), -k, prec2);
      L[1] = evalsigne(1)|evalvarn(0); gel(C,k) = L;
    }
    /* C[k] = \sum_n c_{j,k} t^n =: C_k(t) in Dokchitser's Algo 3.3
     * m-th derivative of t^(-M+2) sum_k (-ln t)^k/k! C_k(t^2) */
    if (m)
    {
      mj = gsubgs(mj, 2);
      for (i = 1; i <= m; i++, mj = gaddgs(mj,1))
        for (k = 1; k <= lj; k++)
        {
          GEN c = gel(C,k), d = RgX_shift_shallow(gmul2n(RgX_deriv(c),1), 1);
          c = RgX_Rg_mul(c, mj);
          if (k < lj) c = RgX_add(c, gel(C,k+1));
          gel(C,k) = RgX_sub(d, c);
        }
      gel(M,j) = gaddgs(mj,2);
    }
    for (k = 1; k <= lj; k++)
    {
      GEN c = gel(C,k);
      if (k > 2) c = RgX_Rg_div(c, mpfact(k-1));
      gel(C,k) = RgX_to_RgC(c, lgpol(c));
    }
  }
  /* Algo 3.3: * \phi^(m)(t) = sum_j t^m_j sum_k (-ln t)^k mat[j,k](t^2) */
  piA = gsubsg(m*d, sumVga(Vga));
  if (!gequal0(piA)) piA = powPis(gmul2n(piA,-1), prec2);
  return mkvec5(L, RgV_neg(M), mat, mkvecsmall(prec2), piA);
}

/* Evaluate a vector considered as a polynomial using Horner. */
static GEN
evalvec(GEN vec, long N, GEN u)
{
  GEN S = gen_0;
  long n;
  N = minss(N, lg(vec)-1);
  for (n = N; n >= 1; n--) S = gmul(u, gadd(gel(vec,n), S));
  return S;
}

/* gammamellininvinit accessors */
static double
get_tmax(long bitprec)
{ return (M_LN2 / MELLININV_CUTOFF) * bitprec ; }
static GEN
GMi_get_Vga(GEN K) { return gel(K,2); }
static long
GMi_get_degree(GEN K) { return lg(gel(K,2))-1; }
static long
GMi_get_m(GEN K) { return itos( gel(K,3) ); }
static GEN /* [lj,mj,mat,prec2], Kderivsmall only */
GMi_get_VS(GEN K) { return gel(K,4); }
/* K[5] = [Ms,cd,A2], Kderivlarge only */
static long/*Kderivlarge*/
GMi_get_status(GEN K) { return itos(gmael3(K,5,1,2)); }
static GEN/*Kderivlarge*/
GMi_get_M(GEN K) { return gmael3(K,5,1,1); }
static GEN/*Kderivlarge*/
GMi_get_cd(GEN K) { return gmael(K,5,2); }
static GEN/*Kderivlarge*/
GMi_get_A2(GEN K) { return gmael(K,5,3); }

static double
GMi_get_tmax(GEN K, long bitprec)
{ return (typ(GMi_get_VS(K)) == t_INT)? -1.0 : get_tmax(bitprec); }

/* Compute m-th derivative of inverse Mellin at x by generalized power series
 * around x = 0; x2d = x^(2/d), x is possibly NULL (don't bother about
 * complex branches). Assume |x|^(2/d) <= tmax = M_LN2*bitprec/MELLININV_CUTOFF*/
static GEN
Kderivsmall(GEN K, GEN x, GEN x2d, long bitprec)
{
  GEN VS = GMi_get_VS(K), L = gel(VS,1), M = gel(VS,2), mat = gel(VS,3);
  GEN d2, Lx, x2, S, pi, piA = gel(VS,5);
  long j, k, prec = gel(VS,4)[1], d = GMi_get_degree(K), limn, N = lg(L)-1;
  double xd, Wd, Ed = M_LN2*bitprec / d;

  xd = maxdd(M_PI*dblmodulus(x2d), 1E-13); /* pi |x|^2/d unless x tiny */
  if (xd > Ed) pari_err_BUG("Kderivsmall (x2d too large)");
  /* Lemma 5.2.6 (2), a = 1 + log(Pi x^(2/d)) = log(e / xd),
   * B = log(2)*bitprec / d = Ed */
  Wd = dbllambertW0( Ed / (M_E*xd) ); /* solution of w exp(w) = B exp(-a)*/
  limn = (long) ceil(2*Ed/Wd);
  pi = mppi(prec);
  d2 = gdivsg(d,gen_2);
  if (x)
    x = gmul(gtofp(x,prec), gpow(pi,d2,prec));
  else
    x = gpow(gmul(gtofp(x2d,prec),pi), d2, prec);
  /* at this stage, x has been replaced by pi^(d/2) x */
  x2 = gsqr(x);
  Lx = gpowers(gneg(glog(x,prec)), vecsmall_max(L));
  S = gen_0;
  for (j = 1; j <= N; ++j)
  {
    long lj = L[j];
    GEN s = gen_0;
    for (k = 1; k <= lj; k++)
      s = gadd(s, gmul(gel(Lx,k), evalvec(gmael(mat,j,k), limn, x2)));
    S = gadd(S, gmul(gpow(x, gel(M,j), prec), s));
  }
  if (!gequal0(piA)) S = gmul(S, piA);
  return S;
}

/* In Klarge, we conpute K(t) as (asymptotic) * F(z), where F ~ 1 is given by
 * a continued fraction and z = Pi t^(2/d). If we take 2n terms in F (n terms
 * in Euler form), F_n(z) - F(z) is experimentally in exp(- C sqrt(n*z))
 * where C ~ 8 for d > 2 [HEURISTIC] and C = 4 (theorem) for d = 1 or d = 2
 * and vga = [0,1]. For e^(-E) absolute error, we want
 *   exp(-C sqrt(nz)) < e^-(E+a), where a ~ ln(asymptotic)
 * i.e.  2n > (E+a)^2 / t^(2/d) * 2/(C^2 Pi); C^2*Pi/2 ~ 100.5 ~ 101
 *
 * In fact, this model becomes wrong for z large: we use instead
 *
 *   exp(- sqrt(D * nz/log(z+1))) < e^-(E+a),
 * i.e.  2n > (E+a)^2 * log(1 + Pi t^(2/d))/ t^(2/d) * 2/(D Pi); */
static double
get_D(long d) { return d <= 2 ? 157. : 180.; }
/* if (abs), absolute error rather than relative */
static void
Kderivlarge_optim(GEN K, int abs, GEN t2d, double cd, long *pbitprec, long *pnlim)
{
  GEN A2 = GMi_get_A2(K);
  long bitprec = *pbitprec, d = GMi_get_degree(K);
  const double D = get_D(d), td = dblmodulus(t2d);
  double a, rtd, E = M_LN2*bitprec;

  /* t = 0 can happen with finite continued fraction or easyvga */
  if (!td) { *pnlim = 0; return; }
  rtd = (typ(t2d) == t_COMPLEX)? gtodouble(gel(t2d,1)): td;
  /* A2/2 = A, log(td) = (2/d)*log t */
  a = d*gtodouble(A2)*log2(td)/2 - (M_PI/M_LN2)*d*rtd + log2(cd);/*log2 K(t)~a*/
  /* if bitprec <= 0, caller should return K(t) ~ 0 */
  bitprec += 64;
  if (abs)
  {
    bitprec += ceil(a);
    if (a <= -65) E = M_LN2*bitprec; /* guarantees E <= initial E */
  }
  *pbitprec = bitprec;
  *pnlim = ceil(E*E * log2(1+M_PI*td) / (D*td));
}

/* Compute m-th derivative of inverse Mellin at t by continued fraction of
 * asymptotic expansion; t2d = t^(2/d). If t is NULL, "lfun" mode: don't
 * bother about complex branches + use absolute (rather than relative)
 * accuracy */
static GEN
Kderivlarge(GEN K, GEN t, GEN t2d, long bitprec0)
{
  GEN tdA, P, S, pi, z;
  const long d = GMi_get_degree(K);
  GEN M = GMi_get_M(K), cd = GMi_get_cd(K), A2 = GMi_get_A2(K);
  long prec, nlim, status = GMi_get_status(K), m = GMi_get_m(K), bitprec = bitprec0;

  Kderivlarge_optim(K, !t, t2d, gtodouble(cd), &bitprec, &nlim);
  if (bitprec <= 0) return gen_0;
  prec = nbits2prec(bitprec);
  t2d = gtofp(t2d, prec);
  if (t)
    tdA = gpow(t, gdivgu(A2,d), prec);
  else
    tdA = gpow(t2d, gdivgu(A2,2), prec);
  pi = mppi(prec); z = gmul(pi, t2d);
  P = gmul(gmul(cd, tdA), gexp(gmulsg(-d, z), prec));
  if (m) P = gmul(P, gpowgs(mulsr(-2, pi), m));
  if (status == 2) /* finite continued fraction */
    S = (lg(M) == 2)? gel(M,1): poleval(M, ginv(z));
  else
  {
    S = contfraceval_inv(M, z, nlim/2);
    if (DEBUGLEVEL>3)
    {
      GEN S0 = contfraceval_inv(M, z, minss(nlim/2 + 1, minss(lg(gel(M, 1)) - 1, lg(gel(M, 2)))));
      long e = gexpo(gmul(P, gsub(S,S0)));
      if (-e < bitprec0)
        err_printf("Kderivlarge: e = %ld, bit = %ld\n",e,bitprec0);
    }
    if (status == 1) S = gmul(S, gsubsg(1, ginv(gmul(z, pi))));
  }
  return gmul(P, S);
}

/* Dokchitser's coefficients used for asymptotic expansion of inverse Mellin
 * 2 <= p <= min(n+1, d), c = 2n-p+1; sh = (sh(x)/x)^(d-p) */
static GEN
vp(long p, long c, GEN SMd, GEN sh)
{
  GEN s, ve = cgetg(p+2, t_VEC);
  long m, j, k;

  gel(ve,1) = gen_1; gel(ve,2) = utoipos(c);
  for (j = 2; j <= p; j++) gel(ve,j+1) = gdivgs(gmulgu(gel(ve,j), c), j);
  s = gel(SMd, 1);
  for (m = 1; m <= p; m++)
  {
    GEN t, c = gel(SMd, m+1);
    if (gequal0(c)) continue;
    t = gel(ve, m+1);
    for (k = 1; k <= m/2; k++)
      t = gadd(t, gmul(gel(ve, m-2*k+1), RgX_coeff(sh, k)));
    s = gadd(s, gmul(c, t));
  }
  return s;
}

static GEN
get_SM(GEN Vga)
{
  long k, m, d = lg(Vga)-1;
  GEN pol, nS1, SM, C, t = vecsum(Vga);

  pol = roots_to_pol(gmulgs(Vga, -d), 0); /* deg(pol) = d */
  SM = cgetg(d+2, t_VEC); gel(SM,1) = gen_1;
  if (gequal0(t))
  { /* shortcut */
    for (m = 1; m <= d; m++) gel(SM,m+1) = gel(pol,d+2-m);
    return SM;
  }
  nS1 = gpowers(gneg(t), d); C = matpascal(d);
  for (m = 1; m <= d; m++)
  {
    pari_sp av = avma;
    GEN s = gmul(gel(nS1, m+1), gcoeff(C, d+1, m+1));
    for (k = 1; k <= m; k++)
    {
      GEN e = gmul(gel(nS1, m-k+1), gcoeff(C, d-k+1, m-k+1));
      s = gadd(s, gmul(e, RgX_coeff(pol, d-k)));
    }
    gel(SM, m+1) = gerepileupto(av, s);
  }
  return SM;
}

static GEN
get_SMd(GEN Vga)
{
  GEN M, SM = get_SM(Vga);
  long p, m, d = lg(Vga)-1;

  M = cgetg(d, t_VEC);
  for (p = 2; p <= d; p++)
  {
    GEN a = gen_1, c;
    long D = d - p;
    gel(M, p-1) = c = cgetg(p+2, t_COL);
    gel(c, 1) = gel(SM, p+1);
    for (m = 1; m <= p; m++)
    {
      a = muliu(a, D + m);
      gel(c, m+1) = gmul(gel(SM, p-m+1), a);
    }
  }
  return M;
}

/* Asymptotic expansion of inverse Mellin, to length nlimmax. Set status = 0
 * (regular), 1 (one Hankel determinant vanishes => contfracinit will fail)
 * or 2 (same as 1, but asymptotic expansion is finite!)
 *
 * If status = 2, the asymptotic expansion is finite so return only
 * the necessary number of terms nlim <= nlimmax + d. */
static GEN
Klargeinit(GEN Vga, long nlimmax, long *status, long prec)
{
  long d = lg(Vga) - 1, p, n, cnt;
  GEN M, SMd, se, vsinh, vd;

  if (Vgaeasytheta(Vga)) { *status = 2; return mkvec(gen_1); }
  /* d >= 2 */
  *status = 0; prec += prec >> 1;
  SMd = get_SMd(Vga);
  se = gsinh(RgX_to_ser(pol_x(0), d+2), 0); setvalser(se,0);
  se = gdeflate(se, 0, 2); /* se(x^2) = sinh(x)/x */
  vsinh = gpowers(se, d);
  vd = gpowers(utoipos(2*d), d);
  M = cgetg(nlimmax + d + 1, t_VEC); gel(M,1) = gen_1;
  for (n = 2, cnt = 0; n <= nlimmax || cnt; n++)
  {
    pari_sp av = avma;
    long ld = minss(d, n);
    GEN s = gen_0;
    for (p = 2; p <= ld; p++)
    {
      GEN z = vp(p, 2*n-1-p, gel(SMd, p-1), gel(vsinh, d-p+1));
      s = gadd(s, gmul(gdiv(z, gel(vd, p+1)), gel(M, n+1-p)));
    }
    if (prec && !isinexact(s)) s = gtofp(s, prec);
    gel(M,n) = s = gerepileupto(av, gdivgs(s, 1-n));
    if (gequal0(s))
    {
      cnt++; *status = 1;
      if (cnt >= d-1) { *status = 2; n -= d-2; break; }
    }
    else
    {
      if (n >= nlimmax) { n++; break; }
      cnt = 0;
    }
  }
  setlg(M, n); return M;
}

/* remove trailing zeros from vector. */
static void
stripzeros(GEN M)
{
  long i;
  for(i = lg(M)-1; i >= 1; --i)
    if (!gequal0(gel(M, i))) break;
  setlg(M, i+1);
}

/* Asymptotic expansion of the m-th derivative of inverse Mellin, to length
 * nlimmax. If status = 2, the asymptotic expansion is finite so return only
 * the necessary number of terms nlim <= nlimmax + d. */
static GEN
gammamellininvasymp_i(GEN Vga, long nlimmax, long m, long *status, long prec)
{
  GEN M, A, Aadd;
  long d, i, nlim, n;

  M = Klargeinit(Vga, nlimmax, status, prec);
  if (!m) return M;
  d = lg(Vga)-1;
  /* half the exponent of t in asymptotic expansion. */
  A = gdivgu(gaddsg(1-d, sumVga(Vga)), 2*d);
  if (*status == 2) M = shallowconcat(M, zerovec(m));
  nlim = lg(M)-1;
  Aadd = sstoQ(2-d, 2*d); /* (1/d) - (1/2) */
  for (i = 1; i <= m; i++, A = gadd(A,Aadd))
    for (n = nlim-1; n >= 1; --n)
      gel(M, n+1) = gsub(gel(M, n+1),
                         gmul(gel(M, n), gsub(A, uutoQ(n-1, d))));
  stripzeros(M); return M;
}

INLINE int
RgV_is_CV(GEN x)
{
  long i;
  for (i = lg(x)-1; i > 0; i--)
  {
    long t = typ(gel(x,i));
    if (!is_real_t(t) && t!= t_COMPLEX) return 0;
  }
  return 1;
}

static GEN
get_Vga(GEN x, GEN *ldata)
{
  if (typ(x)==t_VEC && RgV_is_CV(x)) { *ldata = NULL; return x; }
  *ldata = lfunmisc_to_ldata_shallow_i(x);
  if (*ldata) x = ldata_get_gammavec(*ldata);
  return x;
}
GEN
gammamellininvasymp(GEN Vga, long nlim, long m)
{
  pari_sp av = avma;
  long status;
  GEN ldata;
  Vga = get_Vga(Vga, &ldata);
  if (!is_vec_t(typ(Vga)) || lg(Vga) == 1)
    pari_err_TYPE("gammamellininvasymp",Vga);
  return gerepilecopy(av, gammamellininvasymp_i(Vga, nlim, m, &status, 0));
}

/* Does the continued fraction of the asymptotic expansion M at oo of inverse
 * Mellin transform attached to Vga have zero Hankel determinants ? */
static long
ishankelspec(GEN Vga, GEN M)
{
  long status, i, d = lg(Vga)-1;

  if (d == 5 || d == 7)
  { /* known bad cases: a x 5 and a x 7 */
    GEN v1 = gel(Vga, 1);
    for (i = 2; i <= d; ++i)
      if (!gequal(gel(Vga,i), v1)) break;
    if (i > d) return 1;
  }
  else if (d==10 || d==14)
  { /* [ a x 5 , (a+1) x 5 ] */
    long d2 = d>>1;
    long s0 = 1, s1 = 0, sm1 = 0;
    GEN v1 = gel(Vga, 1);
    for (i = 2; i <= d; i++)
    {
      GEN s = gsub(gel(Vga,i),v1);
      if (gequal0(s)) s0++;
      else if(gequal1(s)) s1++;
      else if(gequalm1(s)) sm1++;
    }
    if (s0==d2 && (s1==d2 || sm1==d2)) return 1;
  }
  status = 0;
  /* Heuristic: if 6 first terms in contfracinit don't fail, assume it's OK */
  pari_CATCH(e_INV) { status = 1; }
  pari_TRY { contfracinit(M, minss(lg(M)-2,6)); }
  pari_ENDCATCH; return status;
}

/* Initialize data for computing m-th derivative of inverse Mellin */
GEN
gammamellininvinit(GEN Vga, long m, long bitprec)
{
  const double C2 = MELLININV_CUTOFF;
  pari_sp ltop = avma;
  GEN A2, M, VS, VL, cd, ldata;
  long nlimmax, status, d, prec = nbits2prec((4*bitprec)/3);
  double E = M_LN2*bitprec, tmax = get_tmax(bitprec); /* = E/C2 */

  if (m < 0)
    pari_err_DOMAIN("gammamellininvinit", "derivation order", "<", gen_0, stoi(m));
  Vga = get_Vga(Vga, &ldata); d = lg(Vga)-1;
  if (!is_vec_t(typ(Vga)) || !d) pari_err_TYPE("gammamellininvinit",Vga);
  nlimmax = ceil(E * log2(1+M_PI*tmax) * C2 / get_D(d));
  A2 = gaddsg(m*(2-d) + 1-d, sumVga(Vga));
  cd = (d <= 2)? gen_2: gsqrt(gdivgu(int2n(d+1), d), nbits2prec(bitprec));
  /* if in Klarge, we have |t| > tmax = E/C2, thus nlim < E*C2/D. */
  M = gammamellininvasymp_i(Vga, nlimmax, m, &status, prec);
  if (status == 2)
  {
    tmax = -1.; /* only use Klarge */
    VS = gen_0;
  }
  else
  {
    VS = Kderivsmallinit(ldata, Vga, m, bitprec);
    if (status == 0 && ishankelspec(Vga, M)) status = 1;
    if (status == 1)
    { /* a Hankel determinant vanishes => contfracinit is undefined.
         So compute K(t) / (1 - 1/(pi^2*t)) instead of K(t)*/
      GEN t = ginv(mppi(prec));
      long i;
      for (i = 2; i < lg(M); ++i)
        gel(M, i) = gadd(gel(M, i), gmul(gel(M, i-1), t));
    }
    else
      M = RgC_gtofp(M, prec); /* convert from rationals to t_REAL: faster */
    M = contfracinit(M, lg(M)-2);
  }
  VL = mkvec3(mkvec2(M, stoi(status)), cd, A2);
  return gerepilecopy(ltop, mkvec5(dbltor(tmax), Vga, stoi(m), VS, VL));
}

/* Compute m-th derivative of inverse Mellin at s2d = s^(d/2) using
 * initialization data. Use Taylor expansion at 0 for |s2d| < tmax, and
 * asymptotic expansion at oo otherwise. WARNING: assume that accuracy
 * has been increased according to tmax by the CALLING program. */
static GEN
gammamellininvrt_i(GEN K, GEN s, GEN s2d, long bit)
{
  if (dblmodulus(s2d) < GMi_get_tmax(K, bit))
    return Kderivsmall(K, s, s2d, bit);
  else
    return Kderivlarge(K, s, s2d, bit);
}
GEN
gammamellininvrt(GEN K, GEN s2d, long bit)
{ return gammamellininvrt_i(K, NULL, s2d, bit); }

/* Compute inverse Mellin at s. K from gammamellininv OR a Vga, in which
 * case the initialization data is computed. */
GEN
gammamellininv(GEN K, GEN s, long m, long bitprec)
{
  pari_sp av = avma;
  GEN s2d;
  long d;

  if (!is_vec_t(typ(K)) || lg(K) != 6 || !is_vec_t(typ(GMi_get_Vga(K))))
    K = gammamellininvinit(K, m, bitprec);
  d = GMi_get_degree(K);
  s2d = gpow(s, gdivgu(gen_2, d), nbits2prec(bitprec));
  return gerepileupto(av, gammamellininvrt_i(K, s, s2d, bitprec));
}
