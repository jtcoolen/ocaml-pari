/* Copyright (C) 2000-2005  The PARI group.

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
/*******************************************************************/
/*                                                                 */
/*         QUADRATIC POLYNOMIAL ASSOCIATED TO A DISCRIMINANT       */
/*                                                                 */
/*******************************************************************/

void
check_quaddisc(GEN x, long *s, long *pr, const char *f)
{
  long r;
  if (typ(x) != t_INT) pari_err_TYPE(f,x);
  *s = signe(x);
  if (Z_issquare(x)) pari_err_DOMAIN(f,"issquare(disc)","=", gen_1,x);
  r = mod4(x); if (*s < 0 && r) r = 4 - r;
  if (r > 1) pari_err_DOMAIN(f,"disc % 4",">", gen_1,x);
  if (pr) *pr = r;
}
void
check_quaddisc_real(GEN x, long *r, const char *f)
{
  long sx; check_quaddisc(x, &sx, r, f);
  if (sx < 0) pari_err_DOMAIN(f, "disc","<",gen_0,x);
}
void
check_quaddisc_imag(GEN x, long *r, const char *f)
{
  long sx; check_quaddisc(x, &sx, r, f);
  if (sx > 0) pari_err_DOMAIN(f, "disc",">",gen_0,x);
}

/* X^2 + b X + c is the canonical quadratic t_POL of discriminant D.
 * Dodd is nonzero iff D is odd */
static void
quadpoly_bc(GEN D, long Dodd, GEN *b, GEN *c)
{
  if (Dodd)
  {
    pari_sp av = avma;
    *b = gen_m1;
    *c = gerepileuptoint(av, shifti(subui(1,D), -2));
  }
  else
  {
    *b = gen_0;
    *c = shifti(D,-2); togglesign(*c);
  }
}
/* X^2 - X - (D-1)/4 or X^2 - D/4 */
static GEN
quadpoly_ii(GEN D, long Dmod4)
{
  GEN b, c, y = cgetg(5,t_POL);
  y[1] = evalsigne(1) | evalvarn(0);
  quadpoly_bc(D, Dmod4, &b,&c);
  gel(y,2) = c;
  gel(y,3) = b;
  gel(y,4) = gen_1; return y;
}
GEN
quadpoly(GEN D)
{
  long s, Dmod4;
  check_quaddisc(D, &s, &Dmod4, "quadpoly");
  return quadpoly_ii(D, Dmod4);
}
GEN /* no checks */
quadpoly_i(GEN D) { return quadpoly_ii(D, Mod4(D)); }

GEN
quadpoly0(GEN x, long v)
{
  GEN T = quadpoly(x);
  if (v > 0) setvarn(T, v);
  return T;
}

GEN
quadgen(GEN x)
{ retmkquad(quadpoly(x), gen_0, gen_1); }

GEN
quadgen0(GEN x, long v)
{
  if (v==-1) v = fetch_user_var("w");
  retmkquad(quadpoly0(x, v), gen_0, gen_1);
}

/***********************************************************************/
/**                                                                   **/
/**                      BINARY QUADRATIC FORMS                       **/
/**                                                                   **/
/***********************************************************************/
static int
is_qfi(GEN q) { return typ(q)==t_QFB && qfb_is_qfi(q); }

static GEN
check_qfbext(const char *fun, GEN x)
{
  long t = typ(x);
  if (t == t_QFB) return x;
  if (t == t_VEC && lg(x)==3)
  {
    GEN q = gel(x,1);
    if (!is_qfi(q) && typ(gel(x,2))==t_REAL) return q;
  }
  pari_err_TYPE(fun, x);
  return NULL;/* LCOV_EXCL_LINE */
}

static GEN
qfb3(GEN x, GEN y, GEN z)
{ retmkqfb(icopy(x), icopy(y), icopy(z), qfb_disc3(x,y,z)); }

static int
qfb_equal(GEN x, GEN y)
{
  return equalii(gel(x,1),gel(y,1))
      && equalii(gel(x,2),gel(y,2))
      && equalii(gel(x,3),gel(y,3));
}

/* valid for t_QFB, qfr3, qfr5; shallow */
static GEN
qfb_inv(GEN x) {
  GEN z = shallowcopy(x);
  gel(z,2) = negi(gel(z,2));
  return z;
}
/* valid for t_QFB, gerepile-safe */
static GEN
qfbinv(GEN x)
{ retmkqfb(icopy(gel(x,1)),negi(gel(x,2)),icopy(gel(x,3)), icopy(gel(x,4))); }

GEN
Qfb0(GEN a, GEN b, GEN c)
{
  GEN q, D;
  if (!b)
  {
    if (c) pari_err_TYPE("Qfb",c);
    if (typ(a) == t_VEC && lg(a) == 4)
    { b = gel(a,2); c = gel(a,3); a = gel(a,1); }
    else if (typ(a) == t_POL && degpol(a) == 2)
    { b = gel(a,3); c = gel(a,2); a = gel(a,4); }
    else if (typ(a) == t_MAT && lg(a)==3 && lgcols(a)==3)
    {
      b = gadd(gcoeff(a,2,1), gcoeff(a,1,2));
      c = gcoeff(a,2,2); a = gcoeff(a,1,1);
    }
    else
      pari_err_TYPE("Qfb",a);
  }
  else if (!c)
    pari_err_TYPE("Qfb",b);
  if (typ(a)!=t_INT) pari_err_TYPE("Qfb",a);
  if (typ(b)!=t_INT) pari_err_TYPE("Qfb",b);
  if (typ(c)!=t_INT) pari_err_TYPE("Qfb",c);
  q = qfb3(a, b, c); D = qfb_disc(q);
  if (signe(D) < 0)
  { if (signe(a) < 0) pari_err_IMPL("negative definite t_QFB"); }
  else if (Z_issquare(D)) pari_err_DOMAIN("Qfb","issquare(disc)","=", gen_1,q);
  return q;
}

/***********************************************************************/
/**                                                                   **/
/**                         Reduction                                 **/
/**                                                                   **/
/***********************************************************************/

/* assume a > 0. Write b = q*2a + r, with -a < r <= a */
static GEN
dvmdii_round(GEN b, GEN a, GEN *r)
{
  GEN a2 = shifti(a, 1), q = dvmdii(b, a2, r);
  if (signe(b) >= 0) {
    if (abscmpii(*r, a) > 0) { q = addiu(q, 1); *r = subii(*r, a2); }
  } else { /* r <= 0 */
    if (abscmpii(*r, a) >= 0){ q = subiu(q, 1); *r = addii(*r, a2); }
  }
  return q;
}
/* Assume 0 < a <= LONG_MAX. Ensure no overflow */
static long
dvmdsu_round(long b, ulong a, long *r)
{
  ulong a2 = a << 1, q, ub, ur;
  if (b >= 0) {
    ub = b;
    q = ub / a2;
    ur = ub % a2;
    if (ur > a) { ur -= a; q++; *r = (long)ur; *r -= (long)a; }
    else *r = (long)ur;
    return (long)q;
  } else { /* r <= 0 */
    ub = (ulong)-b; /* |b| */
    q = ub / a2;
    ur = ub % a2;
    if (ur >= a) { ur -= a; q++; *r = (long)ur; *r = (long)a - *r; }
    else *r = -(long)ur;
    return -(long)q;
  }
}
/* reduce b mod 2*a. Update b,c */
static void
REDB(GEN a, GEN *b, GEN *c)
{
  GEN r, q = dvmdii_round(*b, a, &r);
  if (!signe(q)) return;
  *c = subii(*c, mulii(q, shifti(addii(*b, r),-1)));
  *b = r;
}
/* Assume a > 0. Reduce b mod 2*a. Update b,c */
static void
sREDB(ulong a, long *b, ulong *c)
{
  long r, q;
  ulong uz;
  if (a > LONG_MAX) return; /* b already reduced */
  q = dvmdsu_round(*b, a, &r);
  if (q == 0) return;
  /* Final (a,r,c2) satisfies |r| <= |b| hence c2 <= c, c2 = c - q*z,
   * where z = (b+r) / 2, representable as long, has the same sign as q. */
  if (*b < 0)
  { /* uz = -z >= 0, q < 0 */
    if (r >= 0) /* different signs=>no overflow, exact division */
      uz = (ulong)-((*b + r)>>1);
    else
    {
      ulong ub = (ulong)-*b, ur = (ulong)-r;
      uz = (ub + ur) >> 1;
    }
    *c -= (-q) * uz; /* c -= qz */
  }
  else
  { /* uz = z >= 0, q > 0 */
    if (r <= 0)
      uz = (*b + r)>>1;
    else
    {
      ulong ub = (ulong)*b, ur = (ulong)r;
      uz = ((ub + ur) >> 1);
    }
    *c -= q * uz; /* c -= qz */
  }
  *b = r;
}
static void
REDBU(GEN a, GEN *b, GEN *c, GEN u1, GEN *u2)
{ /* REDB(a,b,c) */
  GEN r, q = dvmdii_round(*b, a, &r);
  *c = subii(*c, mulii(q, shifti(addii(*b, r),-1)));
  *b = r;
  *u2 = subii(*u2, mulii(q, u1));
}

/* q t_QFB, return reduced representative and set base change U in Sl2(Z) */
GEN
redimagsl2(GEN q, GEN *U)
{
  pari_sp av = avma;
  GEN z, u1,u2,v1,v2,Q;
  GEN a = gel(q,1), b = gel(q,2), c = gel(q,3);
  long cmp;
  u1 = gen_1; u2 = gen_0;
  cmp = abscmpii(a, b);
  if (cmp < 0)
    REDBU(a,&b,&c, u1,&u2);
  else if (cmp == 0 && signe(b) < 0)
  { /* b = -a */
    b = negi(b);
    u2 = gen_1;
  }
  for(;;)
  {
    cmp = abscmpii(a, c); if (cmp <= 0) break;
    swap(a,c); b = negi(b);
    z = u1; u1 = u2; u2 = negi(z);
    REDBU(a,&b,&c, u1,&u2);
    if (gc_needed(av, 1)) {
      if (DEBUGMEM>1) pari_warn(warnmem, "redimagsl2");
      gerepileall(av, 5, &a,&b,&c, &u1,&u2);
    }
  }
  if (cmp == 0 && signe(b) < 0)
  {
    b = negi(b);
    z = u1; u1 = u2; u2 = negi(z);
  }
  /* Let q = (A,B,C). q o [u1,u2; v1,v2] = Q implies
   * [v1,v2] = (1/C) [(b-B)/2 u1 - a u2, c u1 - (b+B)/2 u2] */
  z = shifti(subii(b, gel(q,2)), -1);
  v1 = subii(mulii(z, u1), mulii(a, u2)); v1 = diviiexact(v1, gel(q,3));
  z = subii(z, b);
  v2 = addii(mulii(z, u2), mulii(c, u1)); v2 = diviiexact(v2, gel(q,3));
  *U = mkmat2(mkcol2(u1,v1), mkcol2(u2,v2));
  Q = lg(q)==5 ? mkqfb(a,b,c,gel(q,4)): mkvec3(a,b,c);
  return gc_all(av, 2, &Q, U);
}

static GEN
setq_b0(ulong a, ulong c, GEN D)
{ retmkqfb(utoipos(a), gen_0, utoipos(c), icopy(D)); }
/* assume |sb| = 1 */
static GEN
setq(ulong a, ulong b, ulong c, long sb, GEN D)
{ retmkqfb(utoipos(a), sb==1? utoipos(b): utoineg(b), utoipos(c), icopy(D)); }
/* 0 < a, c < 2^BIL, b = 0 */
static GEN
redimag_1_b0(ulong a, ulong c, GEN D)
{ return (a <= c)? setq_b0(a, c, D): setq_b0(c, a, D); }

/* 0 < a, c < 2^BIL: single word affair */
static GEN
redimag_1(pari_sp av, GEN a, GEN b, GEN c, GEN D)
{
  ulong ua, ub, uc;
  long sb;
  for(;;)
  { /* at most twice */
    long lb = lgefint(b); /* <= 3 after first loop */
    if (lb == 2) return redimag_1_b0(a[2],c[2], D);
    if (lb == 3 && uel(b,2) <= (ulong)LONG_MAX) break;
    REDB(a,&b,&c);
    if (uel(a,2) <= uel(c,2))
    { /* lg(b) <= 3 but may be too large for itos */
      long s = signe(b);
      set_avma(av);
      if (!s) return redimag_1_b0(a[2], c[2], D);
      if (a[2] == c[2]) s = 1;
      return setq(a[2], b[2], c[2], s, D);
    }
    swap(a,c); b = negi(b);
  }
  /* b != 0 */
  set_avma(av);
  ua = a[2];
  ub = sb = b[2]; if (signe(b) < 0) sb = -sb;
  uc = c[2];
  if (ua < ub)
    sREDB(ua, &sb, &uc);
  else if (ua == ub && sb < 0) sb = (long)ub;
  while(ua > uc)
  {
    lswap(ua,uc); sb = -sb;
    sREDB(ua, &sb, &uc);
  }
  if (!sb) return setq_b0(ua, uc, D);
  else
  {
    long s = 1;
    if (sb < 0)
    {
      sb = -sb;
      if (ua != uc) s = -1;
    }
    return setq(ua, sb, uc, s, D);
  }
}

static GEN
redimag_av(pari_sp av, GEN q)
{
  GEN a = gel(q,1), b = gel(q,2), c = gel(q,3), D = gel(q,4);
  long cmp, lc = lgefint(c);

  if (lgefint(a) == 3 && lc == 3) return redimag_1(av, a, b, c, D);
  cmp = abscmpii(a, b);
  if (cmp < 0)
    REDB(a,&b,&c);
  else if (cmp == 0 && signe(b) < 0)
    b = negi(b);
  for(;;)
  {
    cmp = abscmpii(a, c); if (cmp <= 0) break;
    lc = lgefint(a); /* lg(future c): we swap a & c next */
    if (lc == 3) return redimag_1(av, a, b, c, D);
    swap(a,c); b = negi(b); /* apply rho */
    REDB(a,&b,&c);
    if (gc_needed(av, 2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"redimag, lc = %ld", lc);
      gerepileall(av, 3, &a,&b,&c);
    }
  }
  if (cmp == 0 && signe(b) < 0) b = negi(b);
  return gerepilecopy(av, mkqfb(a, b, c, D));
}
static GEN
redimag(GEN q) { return redimag_av(avma, q); }

static GEN
rhoimag(GEN x)
{
  pari_sp av = avma;
  GEN a = gel(x,1), b = gel(x,2), c = gel(x,3);
  int fl = abscmpii(a, c);
  if (fl <= 0)
  {
    int fg = abscmpii(a, b);
    if (fg >= 0)
    {
      x = gcopy(x);
      if ((!fl || !fg) && signe(gel(x,2)) < 0) setsigne(gel(x,2), 1);
      return x;
    }
  }
  swap(a,c); b = negi(b);
  REDB(a, &b, &c);
  return gerepilecopy(av, mkqfb(a,b,c, qfb_disc(x)));
}

/* qfr3 / qfr5 */

/* t_QFB are unusable: D, sqrtD, isqrtD are recomputed all the time and the
 * logarithmic Shanks's distance is costly and hard to control.
 * qfr3 / qfr5 routines take a container of t_INTs (e.g a t_VEC) as argument,
 * at least 3 (resp. 5) components [it is a feature that they do not check the
 * precise type or length of the input]. They return a vector of length 3
 * (resp. 5). A qfr3 [a,b,c] contains the form coeffs, in a qfr5 [a,b,c, e,d]
 * the t_INT e is a binary exponent, d a t_REAL, coding the distance in
 * multiplicative form: the true distance is obtained from qfr5_dist.
 * All other qfr routines are obsolete (inefficient) wrappers */

/* static functions are not stack-clean. Unless mentionned otherwise, public
 * functions are. */

#define EMAX 22
static void
fix_expo(GEN x)
{
  if (expo(gel(x,5)) >= (1L << EMAX)) {
    gel(x,4) = addiu(gel(x,4), 1);
    shiftr_inplace(gel(x,5), - (1L << EMAX));
  }
}

/* (1/2) log (d * 2^{e * 2^EMAX}). Not stack clean if e != 0 */
GEN
qfr5_dist(GEN e, GEN d, long prec)
{
  GEN t = logr_abs(d);
  if (signe(e)) {
    GEN u = mulir(e, mplog2(prec));
    shiftr_inplace(u, EMAX); t = addrr(t, u);
  }
  shiftr_inplace(t, -1); return t;
}

static void
rho_get_BC(GEN *B, GEN *C, GEN b, GEN c, struct qfr_data *S)
{
  GEN t, u;
  u = shifti(c,1);
  t = (abscmpii(S->isqrtD,c) >= 0)? S->isqrtD: c;
  u = remii(addii_sign(t,1, b,signe(b)), u);
  *B = addii_sign(t, 1, u, -signe(u)); /* |t| - (|t|+b) % |2c| */
  if (*B == gen_0)
  { u = shifti(S->D, -2); setsigne(u, -1); }
  else
    u = shifti(addii_sign(sqri(*B),1, S->D,-1), -2);
  *C = diviiexact(u, c); /* = (B^2-D)/4c */
}
/* Not stack-clean */
GEN
qfr3_rho(GEN x, struct qfr_data *S)
{
  GEN B, C, b = gel(x,2), c = gel(x,3);
  rho_get_BC(&B, &C, b, c, S);
  return mkvec3(c,B,C);
}
/* Not stack-clean */
GEN
qfr5_rho(GEN x, struct qfr_data *S)
{
  GEN B, C, y, b = gel(x,2), c = gel(x,3);
  long sb = signe(b);

  rho_get_BC(&B, &C, b, c, S);
  y = mkvec5(c,B,C, gel(x,4), gel(x,5));
  if (sb) {
    GEN t = subii(sqri(b), S->D);
    if (sb < 0)
      t = divir(t, sqrr(subir(b,S->sqrtD)));
    else
      t = divri(sqrr(addir(b,S->sqrtD)), t);
    /* t = (b + sqrt(D)) / (b - sqrt(D)), evaluated stably */
    gel(y,5) = mulrr(t, gel(y,5)); fix_expo(y);
  }
  return y;
}

/* Not stack-clean */
GEN
qfr_to_qfr5(GEN x, long prec)
{ return mkvec5(gel(x,1),gel(x,2),gel(x,3),gen_0,real_1(prec)); }

/* d0 = initial distance, x = [a,b,c, expo(d), d], d = exp(2*distance) */
GEN
qfr5_to_qfr(GEN x, GEN D, GEN d0)
{
  if (d0)
  {
    GEN n = gel(x,4), d = absr(gel(x,5));
    if (signe(n))
    {
      n = addis(shifti(n, EMAX), expo(d));
      setexpo(d, 0); d = logr_abs(d);
      if (signe(n)) d = addrr(d, mulir(n, mplog2(lg(d0))));
      shiftr_inplace(d, -1);
      d0 = addrr(d0, d);
    }
    else if (!gequal1(d)) /* avoid loss of precision */
    {
      d = logr_abs(d);
      shiftr_inplace(d, -1);
      d0 = addrr(d0, d);
    }
  }
  x = qfr3_to_qfr(x, D);
  return d0 ? mkvec2(x,d0): x;
}

/* Not stack-clean */
GEN
qfr3_to_qfr(GEN x, GEN d) { retmkqfb(gel(x,1), gel(x,2), gel(x,3), d); }

static int
ab_isreduced(GEN a, GEN b, GEN isqrtD)
{
  if (signe(b) > 0 && abscmpii(b, isqrtD) <= 0)
  {
    GEN t = addii_sign(isqrtD,1, shifti(a,1),-1);
    long l = abscmpii(b, t); /* compare |b| and |floor(sqrt(D)) - |2a|| */
    if (l > 0 || (l == 0 && signe(t) < 0)) return 1;
  }
  return 0;
}

INLINE int
qfr_isreduced(GEN x, GEN isqrtD)
{
  return ab_isreduced(gel(x,1),gel(x,2),isqrtD);
}

/* Not stack-clean */
GEN
qfr5_red(GEN x, struct qfr_data *S)
{
  pari_sp av = avma;
  while (!qfr_isreduced(x, S->isqrtD))
  {
    x = qfr5_rho(x, S);
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"qfr5_red");
      x = gerepilecopy(av, x);
    }
  }
  return x;
}
/* Not stack-clean */
GEN
qfr3_red(GEN x, struct qfr_data *S)
{
  pari_sp av = avma;
  while (!qfr_isreduced(x, S->isqrtD))
  {
    x = qfr3_rho(x, S);
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"qfr3_red");
      x = gerepilecopy(av, x);
    }
  }
  return x;
}

void
qfr_data_init(GEN D, long prec, struct qfr_data *S)
{
  S->D = D;
  S->sqrtD = sqrtr(itor(S->D,prec));
  S->isqrtD = truncr(S->sqrtD);
}

static GEN
qfr5_init(GEN x, GEN d, struct qfr_data *S)
{
  long prec = realprec(d), l = -expo(d);
  if (l < BITS_IN_LONG) l = BITS_IN_LONG;
  prec = maxss(prec, nbits2prec(l));
  S->D = qfb_disc(x);
  x = qfr_to_qfr5(x,prec);
  if (!S->sqrtD) S->sqrtD = sqrtr(itor(S->D,prec));
  else if (typ(S->sqrtD) != t_REAL) pari_err_TYPE("qfr_init",S->sqrtD);

  if (!S->isqrtD)
  {
    pari_sp av=avma;
    long e;
    S->isqrtD = gcvtoi(S->sqrtD,&e);
    if (e>-2) { set_avma(av); S->isqrtD = sqrti(S->D); }
  }
  else if (typ(S->isqrtD) != t_INT) pari_err_TYPE("qfr_init",S->isqrtD);
  return x;
}
static GEN
qfr3_init(GEN x, struct qfr_data *S)
{
  S->D = qfb_disc(x);
  if (!S->isqrtD) S->isqrtD = sqrti(S->D);
  else if (typ(S->isqrtD) != t_INT) pari_err_TYPE("qfr_init",S->isqrtD);
  return x;
}

#define qf_NOD  2
#define qf_STEP 1

static GEN
redreal_i(GEN x, long flag, GEN isqrtD, GEN sqrtD)
{
  struct qfr_data S;
  GEN d = NULL, y;
  if (typ(x)==t_VEC) { d = gel(x,2); x = gel(x,1); } else flag |= qf_NOD;
  S.sqrtD = sqrtD;
  S.isqrtD = isqrtD;
  y = (flag & qf_NOD)? qfr3_init(x, &S): qfr5_init(x, d, &S);
  switch(flag) {
    case 0:              y = qfr5_red(y,&S); break;
    case qf_NOD:         y = qfr3_red(y,&S); break;
    case qf_STEP:        y = qfr5_rho(y,&S); break;
    case qf_STEP|qf_NOD: y = qfr3_rho(y,&S); break;
    default: pari_err_FLAG("qfbred");
  }
  return qfr5_to_qfr(y, qfb_disc(x), d);
}
static GEN
redreal(GEN x) { return redreal_i(x,0,NULL,NULL); }

GEN
qfbred0(GEN x, long flag, GEN isqrtD, GEN sqrtD)
{
  pari_sp av;
  GEN q = check_qfbext("qfbred",x);
  if (qfb_is_qfi(q)) return (flag & qf_STEP)? rhoimag(x): redimag(x);
  if (typ(x)==t_QFB) flag |= qf_NOD;
  else               flag &= ~qf_NOD;
  av = avma;
  return gerepilecopy(av, redreal_i(x,flag,isqrtD,sqrtD));
}
/* t_QFB */
GEN
qfbred_i(GEN x) { return qfb_is_qfi(x)? redimag(x): redreal(x); }
GEN
qfbred(GEN x) { return qfbred0(x, 0, NULL, NULL); }

static void
_rhorealsl2(GEN *pa, GEN *pb, GEN *pc, GEN *pu1, GEN *pu2, GEN *pv1,
            GEN *pv2, GEN d, GEN rd)
{
  GEN C = mpabs_shallow(*pc), t = addii(*pb, gmax_shallow(rd,C));
  GEN r, q = truedvmdii(t, shifti(C,1), &r);
  *pb = subii(t, addii(r, *pb));
  *pa = *pc;
  *pc = diviiexact(subii(sqri(*pb), d), shifti(*pa, 2));
  if (signe(*pa) < 0) togglesign(q);
  r = *pu1; *pu1 = *pv1; *pv1 = subii(mulii(q, *pv1), r);
  r = *pu2; *pu2 = *pv2; *pv2 = subii(mulii(q, *pv2), r);
}

static GEN
rhorealsl2(GEN A, GEN rd)
{
  GEN V = gel(A,1), M = gel(A,2);
  GEN a = gel(V,1), b = gel(V,2), c = gel(V,3), d = qfb_disc(V);
  GEN u1 = gcoeff(M,1,1), v1 = gcoeff(M,1,2);
  GEN u2 = gcoeff(M,2,1), v2 = gcoeff(M,2,2);
  _rhorealsl2(&a,&b,&c, &u1,&u2,&v1,&v2, d, rd);
  return mkvec2(mkqfb(a,b,c,d), mkmat22(u1,v1,u2,v2));
}

static GEN
redrealsl2(GEN V, GEN rd)
{
  pari_sp av = avma;
  GEN u1 = gen_1, u2 = gen_0, v1 = gen_0, v2 = gen_1;
  GEN a = gel(V,1), b = gel(V,2), c = gel(V,3), d = qfb_disc(V);
  while (!ab_isreduced(a,b,rd))
  {
    _rhorealsl2(&a,&b,&c, &u1,&u2,&v1,&v2, d, rd);
    if (gc_needed(av, 1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"redrealsl2");
      gerepileall(av, 7, &a,&b,&c,&u1,&u2,&v1,&v2);
    }
  }
  return gerepilecopy(av, mkvec2(mkqfb(a,b,c,d), mkmat22(u1,v1,u2,v2)));
}

GEN
qfbredsl2(GEN q, GEN isD)
{
  pari_sp av;
  if (typ(q) != t_QFB) pari_err_TYPE("qfbredsl2",q);
  if (qfb_is_qfi(q))
  {
    GEN v = cgetg(3,t_VEC);
    if (isD) pari_err_TYPE("qfbredsl2", isD);
    gel(v,1) = redimagsl2(q, &gel(v,2)); return v;
  }
  av = avma;
  if (!isD) isD = sqrti(qfb_disc(q));
  else if (typ(isD) != t_INT) pari_err_TYPE("qfbredsl2",isD);
  return gerepileupto(av, redrealsl2(q, isD));
}



/***********************************************************************/
/**                                                                   **/
/**                         Composition                               **/
/**                                                                   **/
/***********************************************************************/

static void
qfb_sqr(GEN z, GEN x)
{
  GEN c, d1, x2, v1, v2, c3, m, p1, r;

  d1 = bezout(gel(x,2),gel(x,1),&x2, NULL); /* usually 1 */
  c = gel(x,3);
  m = mulii(c,x2);
  if (equali1(d1))
    v1 = v2 = gel(x,1);
  else
  {
    v1 = diviiexact(gel(x,1),d1);
    v2 = mulii(v1, gcdii(d1,c)); /* = v1 iff x primitive */
    c = mulii(c, d1);
  }
  togglesign(m);
  r = modii(m,v2);
  p1 = mulii(r, v1);
  c3 = addii(c, mulii(r,addii(gel(x,2),p1)));
  gel(z,1) = mulii(v1,v2);
  gel(z,2) = addii(gel(x,2), shifti(p1,1));
  gel(z,3) = diviiexact(c3,v2);
}
/* z <- x * y */
static void
qfb_comp(GEN z, GEN x, GEN y)
{
  GEN n, c, d, y1, v1, v2, c3, m, p1, r;

  if (x == y) { qfb_sqr(z,x); return; }
  n = shifti(subii(gel(y,2),gel(x,2)), -1);
  v1 = gel(x,1);
  v2 = gel(y,1);
  c  = gel(y,3);
  d = bezout(v2,v1,&y1,NULL);
  if (equali1(d))
    m = mulii(y1,n);
  else
  {
    GEN s = subii(gel(y,2), n);
    GEN x2, y2, d1 = bezout(s,d,&x2,&y2); /* x2 s + y2 (x1 v1 + y1 v2) = d1 */
    if (!equali1(d1))
    {
      v1 = diviiexact(v1,d1);
      v2 = diviiexact(v2,d1); /* gcd = 1 iff x or y primitive */
      v1 = mulii(v1, gcdii(c,gcdii(gel(x,3),gcdii(d1,n))));
      c = mulii(c, d1);
    }
    m = addii(mulii(mulii(y1,y2),n), mulii(gel(y,3),x2));
  }
  togglesign(m);
  r = modii(m, v1);
  p1 = mulii(r, v2);
  c3 = addii(c, mulii(r,addii(gel(y,2),p1)));
  gel(z,1) = mulii(v1,v2);
  gel(z,2) = addii(gel(y,2), shifti(p1,1));
  gel(z,3) = diviiexact(c3,v1);
}

/* not meant to be efficient */
static GEN
qfb_comp_gen(GEN x, GEN y)
{
  GEN d1 = qfb_disc(x), d2 = qfb_disc(y);
  GEN a1 = gel(x,1), b1 = gel(x,2), c1 = gel(x,3), n1;
  GEN a2 = gel(y,1), b2 = gel(y,2), c2 = gel(y,3), n2;
  GEN cx = content(x), cy = content(y), A, B, C, D, U, m, m2;

  if (!is_pm1(cx))
  {
    a1 = diviiexact(a1, cx); b1 = diviiexact(b1, cx);
    c1 = diviiexact(c1, cx); d1 = diviiexact(d1, sqri(cx));
  }
  if (!is_pm1(cy))
  {
    a2 = diviiexact(a2, cy); c2 = diviiexact(c2, cy);
    b2 = diviiexact(b2, cy); d2 = diviiexact(d2, sqri(cy));
  }
  D = gcdii(d1, d2); if (signe(d1) < 0) setsigne(D, -1);
  if (!Z_issquareall(diviiexact(d1, D), &n1) ||
      !Z_issquareall(diviiexact(d2, D), &n2)) return NULL;
  A = mulii(a1, n2);
  B = mulii(a2, n1);
  C = shifti(addii(mulii(b1, n2), mulii(b2, n1)), -1);
  U = ZV_extgcd(mkvec3(A, B, C));
  m = gel(U,1); U = gmael(U,2,3);
  A = mulii(diviiexact(mulii(a1,b2),m), gel(U,1));
  B = mulii(diviiexact(mulii(a2,b1),m), gel(U,2));
  C = addii(mulii(b1,b2), mulii(D, mulii(n1,n2)));
  C = mulii(diviiexact(shifti(C,-1), m), gel(U,3));
  B = addii(A, addii(B, C));
  m2 = sqri(m);
  A = diviiexact(mulii(a1, a2), m2);
  C = diviiexact(shifti(subii(sqri(B),D), -2), A);
  cx = mulii(cx, cy);
  if (!is_pm1(cx))
  {
    A = mulii(A, cx); B = mulii(B, cx);
    C = mulii(C, cx); D = mulii(D, sqri(cx));
  }
  return mkqfb(A, B, C, D);
}

static GEN redimag_av(pari_sp av, GEN q);
static GEN
qficomp0(GEN x, GEN y, int raw)
{
  pari_sp av = avma;
  GEN z = cgetg(5,t_QFB);
  gel(z,4) = gel(x,4);
  qfb_comp(z, x,y);
  if (raw) return gerepilecopy(av,z);
  return redimag_av(av, z);
}
static GEN redreal(GEN x);
static GEN
qfrcomp0(GEN x, GEN y, int raw)
{
  pari_sp av = avma;
  GEN dx = NULL, dy = NULL;
  GEN z = cgetg(5,t_QFB);
  if (typ(x)==t_VEC) { dx = gel(x,2); x = gel(x,1); }
  if (typ(y)==t_VEC) { dy = gel(y,2); y = gel(y,1); }
  gel(z,4) = gel(x,4);
  qfb_comp(z, x,y);
  if (dx) z = mkvec2(z, dy? addrr(dx, dy): dx); else if (dy) z = mkvec2(z, dy);
  if (!raw) z = redreal(z);
  return gerepilecopy(av, z);
}
/* same discriminant, no distance, no checks */
GEN
qfbcomp_i(GEN x, GEN y)
{ return qfb_is_qfi(x)? qficomp0(x,y,0): qfrcomp0(x,y,0); }
GEN
qfbcomp(GEN x, GEN y)
{
  GEN qx = check_qfbext("qfbcomp", x);
  GEN qy = check_qfbext("qfbcomp", y);
  if (!equalii(gel(qx,4),gel(qy,4)))
  {
    pari_sp av = avma;
    GEN z = qfb_comp_gen(qx, qy);
    if (typ(x) == t_VEC || typ(y) == t_VEC)
      pari_err_IMPL("Shanks's distance in general composition");
    if (!z) pari_err_OP("*",x,y);
    return gerepileupto(av, qfbred(z));
  }
  return qfb_is_qfi(qx)? qficomp0(x,y,0): qfrcomp0(x,y,0);
}
/* same discriminant, no distance, no checks */
GEN
qfbcompraw_i(GEN x, GEN y)
{ return qfb_is_qfi(x)? qficomp0(x,y,1): qfrcomp0(x,y,1); }
GEN
qfbcompraw(GEN x, GEN y)
{
  GEN qx = check_qfbext("qfbcompraw", x);
  GEN qy = check_qfbext("qfbcompraw", y);
  if (!equalii(gel(qx,4),gel(qy,4)))
  {
    pari_sp av = avma;
    GEN z = qfb_comp_gen(qx, qy);
    if (typ(x) == t_VEC || typ(y) == t_VEC)
      pari_err_IMPL("Shanks's distance in general composition");
    if (!z) pari_err_OP("qfbcompraw",x,y);
    return gerepilecopy(av, z);
  }
  if (!equalii(gel(qx,4),gel(qy,4))) pari_err_OP("qfbcompraw",x,y);
  return qfb_is_qfi(qx)? qficomp0(x,y,1): qfrcomp0(x,y,1);
}

static GEN
qfisqr0(GEN x, long raw)
{
  pari_sp av = avma;
  GEN z = cgetg(5,t_QFB);
  gel(z,4) = gel(x,4);
  qfb_sqr(z,x);
  if (raw) return gerepilecopy(av,z);
  return redimag_av(av, z);
}
static GEN
qfrsqr0(GEN x, long raw)
{
  pari_sp av = avma;
  GEN dx = NULL, z = cgetg(5,t_QFB);
  if (typ(x) == t_VEC) { dx = gel(x,2); x = gel(x,1); }
  gel(z,4) = gel(x,4); qfb_sqr(z,x);
  if (dx) z = mkvec2(z, shiftr(dx,1));
  if (!raw) z = redreal(z);
  return gerepilecopy(av, z);
}
/* same discriminant, no distance, no checks */
GEN
qfbsqr_i(GEN x)
{ return qfb_is_qfi(x)? qfisqr0(x,0): qfrsqr0(x,0); }
GEN
qfbsqr(GEN x)
{
  GEN qx = check_qfbext("qfbsqr", x);
  return qfb_is_qfi(qx)? qfisqr0(x,0): qfrsqr0(x,0);
}

static GEN
qfr_1_by_disc(GEN D)
{
  GEN y, r, s;
  check_quaddisc_real(D, NULL, "qfr_1_by_disc");
  y = cgetg(5,t_QFB);
  s = sqrtremi(D, &r); togglesign(r); /* s^2 - r = D */
  if (mpodd(r))
  {
    s = subiu(s,1);
    r = subii(r, addiu(shifti(s, 1), 1));
    r = shifti(r, -2); set_avma((pari_sp)y); s = icopy(s);
  }
  else
  { r = shifti(r, -2); set_avma((pari_sp)s); }
  gel(y,1) = gen_1;
  gel(y,2) = s;
  gel(y,3) = icopy(r);
  gel(y,4) = icopy(D); return y;
}

static GEN
qfr_disc(GEN x)
{ return qfb_disc(typ(x)==t_VEC ? gel(x,1): x); }

static GEN
qfr_1(GEN x)
{ return qfr_1_by_disc(qfr_disc(x)); }

static void
qfr_1_fill(GEN y, struct qfr_data *S)
{
  pari_sp av = avma;
  GEN y2 = S->isqrtD;
  gel(y,1) = gen_1;
  if (mod2(S->D) != mod2(y2)) y2 = subiu(y,1);
  gel(y,2) = y2; av = avma;
  gel(y,3) = gerepileuptoint(av, shifti(subii(sqri(y2), S->D),-2));
}
static GEN
qfr5_1(struct qfr_data *S, long prec)
{
  GEN y = cgetg(6, t_VEC);
  qfr_1_fill(y, S);
  gel(y,4) = gen_0;
  gel(y,5) = real_1(prec); return y;
}
static GEN
qfr3_1(struct qfr_data *S)
{
  GEN y = cgetg(4, t_VEC);
  qfr_1_fill(y, S); return y;
}

/* Assume D < 0 is the discriminant of a t_QFB */
static GEN
qfi_1_by_disc(GEN D)
{
  GEN b,c, y = cgetg(5,t_QFB);
  quadpoly_bc(D, mod2(D), &b,&c);
  if (b == gen_m1) b = gen_1;
  gel(y,1) = gen_1;
  gel(y,2) = b;
  gel(y,3) = c;
  gel(y,4) = icopy(D); return y;
}
static GEN
qfi_1(GEN x)
{
  if (typ(x) != t_QFB) pari_err_TYPE("qfi_1",x);
  return qfi_1_by_disc(qfb_disc(x));
}

GEN
qfb_1(GEN x) { return qfb_is_qfi(x) ? qfi_1(x): qfr_1(x); }

static GEN
_qfimul(void *E, GEN x, GEN y) { (void) E; return qficomp0(x,y,0); }
static GEN
_qfisqr(void *E, GEN x) { (void) E; return qficomp0(x,x,0); }
static GEN
_qfimulraw(void *E, GEN x, GEN y) { (void) E; return qficomp0(x,y,1); }
static GEN
_qfisqrraw(void *E, GEN x) { (void) E; return qficomp0(x,x,1); }

static GEN
qfipowraw(GEN x, long n)
{
  pari_sp av = avma;
  GEN y;
  if (!n) return qfi_1(x);
  if (n== 1) return gcopy(x);
  if (n==-1) { x = gcopy(x); togglesign(gel(x,2)); return x; }
  if (n < 0) x = qfb_inv(x);
  y = gen_powu(x, labs(n), NULL, &_qfisqrraw, &_qfimulraw);
  return gerepilecopy(av,y);
}

static GEN
qfipow(GEN x, GEN n)
{
  pari_sp av = avma;
  GEN y;
  long s = signe(n);
  if (!s) return qfi_1(x);
  if (s < 0) x = qfb_inv(x);
  y = gen_pow(qfbred_i(x), n, NULL, &_qfisqr, &_qfimul);
  return gerepilecopy(av,y);
}

static long
parteucl(GEN L, GEN *d, GEN *v3, GEN *v, GEN *v2)
{
  long z;
  *v = gen_0; *v2 = gen_1;
  for (z=0; abscmpii(*v3,L) > 0; z++)
  {
    GEN t3, t2 = subii(*v, mulii(truedvmdii(*d,*v3,&t3),*v2));
    *v = *v2; *d = *v3; *v2 = t2; *v3 = t3;
  }
  return z;
}

/* composition: Shanks' NUCOMP & NUDUPL */
/* L = floor((|d|/4)^(1/4)) */
GEN
nucomp(GEN x, GEN y, GEN L)
{
  pari_sp av = avma;
  long z;
  GEN a, a1, a2, b2, b, d, d1, g, n, p1, q1, q2, s, u, u1, v, v2, v3, Q;

  if (x==y) return nudupl(x,L);
  if (!is_qfi(x)) pari_err_TYPE("nucomp",x);
  if (!is_qfi(y)) pari_err_TYPE("nucomp",y);

  if (abscmpii(gel(x,1),gel(y,1)) < 0) swap(x, y);
  s = shifti(addii(gel(x,2),gel(y,2)), -1);
  n = subii(gel(y,2), s);
  a1 = gel(x,1);
  a2 = gel(y,1); d = bezout(a2,a1,&u,&v);
  if (equali1(d)) { a = negi(mulii(u,n)); d1 = d; }
  else if (dvdii(s,d)) /* d | s */
  {
    a = negi(mulii(u,n)); d1 = d;
    a1 = diviiexact(a1, d1);
    a2 = diviiexact(a2, d1);
    s = diviiexact(s, d1);
  }
  else
  {
    GEN p2, l;
    d1 = bezout(s,d,&u1,NULL);
    if (!equali1(d1))
    {
      a1 = diviiexact(a1,d1);
      a2 = diviiexact(a2,d1);
      s = diviiexact(s,d1);
      d = diviiexact(d,d1);
    }
    p1 = remii(gel(x,3),d);
    p2 = remii(gel(y,3),d);
    l = modii(mulii(negi(u1), addii(mulii(u,p1),mulii(v,p2))), d);
    a = subii(mulii(l,diviiexact(a1,d)), mulii(u,diviiexact(n,d)));
  }
  a = modii(a,a1); p1 = subii(a,a1); if (abscmpii(a,p1) > 0) a = p1;
  d = a1; v3 = a; z = parteucl(L, &d,&v3, &v,&v2);
  Q = cgetg(5,t_QFB);
  if (!z) {
    g = diviiexact(addii(mulii(v3,s),gel(y,3)), d);
    b = a2;
    b2 = gel(y,2);
    v2 = d1;
    gel(Q,1) = mulii(d,b);
  } else {
    GEN e, q3, q4;
    if (z&1) { v3 = negi(v3); v2 = negi(v2); }
    b = diviiexact(addii(mulii(a2,d), mulii(n,v)), a1);
    e = diviiexact(addii(mulii(s,d),mulii(gel(y,3),v)), a1);
    q3 = mulii(e,v2);
    q4 = subii(q3,s);
    b2 = addii(q3,q4);
    g = diviiexact(q4,v);
    if (!equali1(d1)) { v2 = mulii(d1,v2); v = mulii(d1,v); b2 = mulii(d1,b2); }
    gel(Q,1) = addii(mulii(d,b), mulii(e,v));
  }
  q1 = mulii(b, v3);
  q2 = addii(q1,n);
  gel(Q,2) = addii(b2, z? addii(q1,q2): shifti(q1, 1));
  gel(Q,3) = addii(mulii(v3,diviiexact(q2,d)), mulii(g,v2));
  gel(Q,4) = gel(x,4);
  return redimag_av(av, Q);
}

GEN
nudupl(GEN x, GEN L)
{
  pari_sp av = avma;
  long z;
  GEN u, v, d, d1, p1, a, b, c, a2, b2, c2, Q, v2, v3, g;

  if (!is_qfi(x)) pari_err_TYPE("nudupl",x);
  a = gel(x,1);
  b = gel(x,2);
  d1 = bezout(b,a, &u,NULL);
  if (!equali1(d1))
  {
    a = diviiexact(a, d1);
    b = diviiexact(b, d1);
  }
  c = modii(negi(mulii(u,gel(x,3))), a);
  p1 = subii(c,a); if (abscmpii(c,p1) > 0) c = p1;
  d = a; v3 = c; z = parteucl(L, &d,&v3, &v,&v2);
  a2 = sqri(d);
  c2 = sqri(v3);
  Q = cgetg(5,t_QFB);
  if (!z) {
    g = diviiexact(addii(mulii(v3,b),gel(x,3)), d);
    b2 = gel(x,2);
    v2 = d1;
    gel(Q,1) = a2;
  } else {
    GEN e;
    if (z&1) { v = negi(v); d = negi(d); }
    e = diviiexact(addii(mulii(gel(x,3),v), mulii(b,d)), a);
    g = diviiexact(subii(mulii(e,v2), b), v);
    b2 = addii(mulii(e,v2), mulii(v,g));
    if (!equali1(d1)) { b2 = mulii(d1,b2); v = mulii(d1,v); v2 = mulii(d1,v2); }
    gel(Q,1) = addii(a2, mulii(e,v));
  }
  gel(Q,2) = addii(b2, subii(sqri(addii(d,v3)), addii(a2,c2)));
  gel(Q,3) = addii(c2, mulii(g,v2));
  gel(Q,4) = gel(x,4);
  return redimag_av(av, Q);
}

static GEN
mul_nucomp(void *l, GEN x, GEN y) { return nucomp(x, y, (GEN)l); }
static GEN
mul_nudupl(void *l, GEN x) { return nudupl(x, (GEN)l); }
GEN
nupow(GEN x, GEN n, GEN L)
{
  pari_sp av;
  GEN y, D;

  if (typ(n) != t_INT) pari_err_TYPE("nupow",n);
  if (!is_qfi(x)) pari_err_TYPE("nupow",x);
  if (gequal1(n)) return gcopy(x);
  av = avma;
  D = qfb_disc(x);
  y = qfi_1_by_disc(D);
  if (!signe(n)) return y;
  if (!L) L = sqrtnint(absi_shallow(D), 4);
  y = gen_pow_i(x, n, (void*)L, &mul_nudupl, &mul_nucomp);
  if (signe(n) < 0
  && !absequalii(gel(y,1),gel(y,2))
  && !absequalii(gel(y,1),gel(y,3))) togglesign(gel(y,2));
  return gerepilecopy(av, y);
}

/* Not stack-clean */
GEN
qfr5_compraw(GEN x, GEN y)
{
  GEN z = cgetg(6,t_VEC); qfb_comp(z,x,y);
  if (x == y)
  {
    gel(z,4) = shifti(gel(x,4),1);
    gel(z,5) = sqrr(gel(x,5));
  }
  else
  {
    gel(z,4) = addii(gel(x,4),gel(y,4));
    gel(z,5) = mulrr(gel(x,5),gel(y,5));
  }
  fix_expo(z); return z;
}
GEN
qfr5_comp(GEN x, GEN y, struct qfr_data *S)
{ return qfr5_red(qfr5_compraw(x, y), S); }
/* Not stack-clean */
GEN
qfr3_compraw(GEN x, GEN y)
{
  GEN z = cgetg(4,t_VEC); qfb_comp(z,x,y);
  return z;
}
GEN
qfr3_comp(GEN x, GEN y, struct qfr_data *S)
{ return qfr3_red(qfr3_compraw(x,y), S); }

/* m > 0. Not stack-clean */
static GEN
qfr5_powraw(GEN x, long m)
{
  GEN y = NULL;
  for (; m; m >>= 1)
  {
    if (m&1) y = y? qfr5_compraw(y,x): x;
    if (m == 1) break;
    x = qfr5_compraw(x,x);
  }
  return y;
}

/* return x^n. Not stack-clean */
GEN
qfr5_pow(GEN x, GEN n, struct qfr_data *S)
{
  GEN y = NULL;
  long i, m, s = signe(n);
  if (!s) return qfr5_1(S, lg(gel(x,5)));
  if (s < 0) x = qfb_inv(x);
  for (i=lgefint(n)-1; i>1; i--)
  {
    m = n[i];
    for (; m; m>>=1)
    {
      if (m&1) y = y? qfr5_comp(y,x,S): x;
      if (m == 1 && i == 2) break;
      x = qfr5_comp(x,x,S);
    }
  }
  return y;
}
/* m > 0; return x^m. Not stack-clean */
static GEN
qfr3_powraw(GEN x, long m)
{
  GEN y = NULL;
  for (; m; m>>=1)
  {
    if (m&1) y = y? qfr3_compraw(y,x): x;
    if (m == 1) break;
    x = qfr3_compraw(x,x);
  }
  return y;
}
/* return x^n. Not stack-clean */
GEN
qfr3_pow(GEN x, GEN n, struct qfr_data *S)
{
  GEN y = NULL;
  long i, m, s = signe(n);
  if (!s) return qfr3_1(S);
  if (s < 0) x = qfb_inv(x);
  for (i=lgefint(n)-1; i>1; i--)
  {
    m = n[i];
    for (; m; m>>=1)
    {
      if (m&1) y = y? qfr3_comp(y,x,S): x;
      if (m == 1 && i == 2) break;
      x = qfr3_comp(x,x,S);
    }
  }
  return y;
}

static GEN
qfrinvraw(GEN x)
{
  if (typ(x) == t_VEC) retmkvec2(qfbinv(gel(x,1)), negr(gel(x,2)));
 return qfbinv(x);
}
static GEN
qfrpowraw(GEN x, long n)
{
  struct qfr_data S = { NULL, NULL, NULL };
  pari_sp av = avma;
  if (n==1) return gcopy(x);
  if (n==-1) return qfrinvraw(x);
  if (typ(x)==t_QFB)
  {
    GEN D = qfb_disc(x);
    if (!n) return qfr_1(x);
    if (n < 0) { x = qfb_inv(x); n = -n; }
    x = qfr3_powraw(x, n);
    x = qfr3_to_qfr(x, D);
  }
  else
  {
    GEN d0 = gel(x,2);
    x = gel(x,1);
    if (!n) retmkvec2(qfr_1(x), real_0(precision(d0)));
    if (n < 0) { x = qfb_inv(x); n = -n; }
    x = qfr5_init(x, d0, &S);
    if (labs(n) != 1) x = qfr5_powraw(x, n);
    x = qfr5_to_qfr(x, S.D, mulrs(d0,n));
  }
  return gerepilecopy(av, x);
}
static GEN
qfrpow(GEN x, GEN n)
{
  struct qfr_data S = { NULL, NULL, NULL };
  long s = signe(n);
  pari_sp av = avma;
  if (typ(x)==t_QFB)
  {
    if (!s) return qfr_1(x);
    if (s < 0) x = qfb_inv(x);
    x = qfr3_init(x, &S);
    x = is_pm1(n)? qfr3_red(x, &S): qfr3_pow(x, n, &S);
    x = qfr3_to_qfr(x, S.D);
  }
  else
  {
    GEN d0 = gel(x,2);
    x = gel(x,1);
    if (!s) retmkvec2(qfr_1(x), real_0(precision(d0)));
    if (s < 0) x = qfb_inv(x);
    x = qfr5_init(x, d0, &S);
    x = is_pm1(n)? qfr5_red(x, &S): qfr5_pow(x, n, &S);
    x = qfr5_to_qfr(x, S.D, mulri(d0,n));
  }
  return gerepilecopy(av, x);
}
GEN
qfbpowraw(GEN x, long n)
{
  GEN q = check_qfbext("qfbpowraw",x);
  return qfb_is_qfi(q)? qfipowraw(x,n): qfrpowraw(x,n);
}
/* same discriminant, no distance, no checks */
GEN
qfbpow_i(GEN x, GEN n) { return qfb_is_qfi(x)? qfipow(x,n): qfrpow(x,n); }
GEN
qfbpow(GEN x, GEN n)
{
  GEN q = check_qfbext("qfbpow",x);
  return qfb_is_qfi(q)? qfipow(x,n): qfrpow(x,n);
}
GEN
qfbpows(GEN x, long n)
{
  long N[] = { evaltyp(t_INT) | _evallg(3), 0, 0};
  affsi(n, N); return qfbpow(x, N);
}

/* Prime forms attached to prime ideals of degree 1 */

/* assume x != 0 a t_INT, p > 0
 * Return a t_QFB, but discriminant sign is not checked: can be used for
 * real forms as well */
GEN
primeform_u(GEN x, ulong p)
{
  GEN c, y = cgetg(5, t_QFB);
  pari_sp av = avma;
  ulong b;
  long s;

  s = mod8(x); if (signe(x) < 0 && s) s = 8-s;
  /* 2 or 3 mod 4 */
  if (s & 2) pari_err_DOMAIN("primeform", "disc % 4", ">",gen_1, x);
  if (p == 2) {
    switch(s) {
      case 0: b = 0; break;
      case 1: b = 1; break;
      case 4: b = 2; break;
      default: pari_err_SQRTN("primeform", mkintmod(x,utoi(p)) );
               b = 0; /* -Wall */
    }
    c = shifti(subsi(s,x), -3);
  } else {
    b = Fl_sqrt(umodiu(x,p), p);
    if (b == ~0UL) pari_err_SQRTN("primeform", mkintmod(x,utoi(p)) );
    /* mod(b) != mod2(x) ? */
    if ((b ^ s) & 1) b = p - b;
    c = diviuexact(shifti(subii(sqru(b), x), -2), p);
  }
  gel(y,3) = gerepileuptoint(av, c);
  gel(y,4) = icopy(x);
  gel(y,2) = utoi(b);
  gel(y,1) = utoipos(p); return y;
}

/* special case: p = 1 return unit form */
GEN
primeform(GEN x, GEN p)
{
  const char *f = "primeform";
  pari_sp av;
  long s, sx = signe(x), sp = signe(p);
  GEN y, b, absp;

  if (typ(x) != t_INT) pari_err_TYPE(f,x);
  if (typ(p) != t_INT) pari_err_TYPE(f,p);
  if (!sp) pari_err_DOMAIN(f,"p","=",gen_0,p);
  if (!sx) pari_err_DOMAIN(f,"D","=",gen_0,x);
  if (lgefint(p) == 3)
  {
    ulong pp = p[2];
    if (pp == 1) {
      if (sx < 0) {
        long r;
        if (sp < 0) pari_err_IMPL("negative definite t_QFB");
        r = mod4(x);
        if (r && r != 3) pari_err_DOMAIN(f,"disc % 4",">", gen_1,x);
        return qfi_1_by_disc(x);
      }
      y = qfr_1_by_disc(x);
      if (sp < 0) { gel(y,1) = gen_m1; togglesign(gel(y,3)); }
      return y;
    }
    y = primeform_u(x, pp);
    if (sx < 0) {
      if (sp < 0) pari_err_IMPL("negative definite t_QFB");
      return y;
    }
    if (sp < 0) { togglesign(gel(y,1)); togglesign(gel(y,3)); }
    return gcopy( qfr3_to_qfr(y, x) );
  }
  s = mod8(x);
  if (sx < 0)
  {
    if (sp < 0) pari_err_IMPL("negative definite t_QFB");
    if (s) s = 8-s;
  }
  y = cgetg(5, t_QFB);
  /* 2 or 3 mod 4 */
  if (s & 2) pari_err_DOMAIN(f, "disc % 4", ">",gen_1, x);
  absp = absi_shallow(p); av = avma;
  b = Fp_sqrt(x, absp); if (!b) pari_err_SQRTN(f, mkintmod(x,absp));
  s &= 1; /* s = x mod 2 */
  /* mod(b) != mod2(x) ? [Warning: we may have b == 0] */
  if ((!signe(b) && s) || mod2(b) != s) b = gerepileuptoint(av, subii(absp,b));

  av = avma;
  gel(y,3) = gerepileuptoint(av, diviiexact(shifti(subii(sqri(b), x), -2), p));
  gel(y,4) = icopy(x);
  gel(y,2) = b;
  gel(y,1) = icopy(p);
  return y;
}

static GEN
normforms(GEN D, GEN fa)
{
  long i, j, k, lB, aN, sa;
  GEN a, L, V, B, N, N2;
  int D_odd = mpodd(D);
  a = typ(fa) == t_INT ? fa: typ(fa) == t_VEC? gel(fa,1): factorback(fa);
  sa = signe(a);
  if (sa==0 || (signe(D)<0 && sa<0)) return NULL;
  V = D_odd? Zn_quad_roots(fa, gen_1, shifti(subsi(1, D), -2))
           : Zn_quad_roots(fa, gen_0, negi(shifti(D, -2)));
  if (!V) return NULL;
  N = gel(V,1); B = gel(V,2); lB = lg(B);
  N2 = shifti(N,1);
  aN = itou(diviiexact(a, N)); /* |a|/N */
  L = cgetg((lB-1)*aN+1, t_VEC);
  for (k = 1, i = 1; i < lB; i++)
  {
    GEN b = shifti(gel(B,i), 1), c, C;
    if (D_odd) b = addiu(b , 1);
    c = diviiexact(shifti(subii(sqri(b), D), -2), a);
    for (j = 0;; b = addii(b, N2))
    {
      gel(L, k++) = mkqfb(a, b, c, D);
      if (++j == aN) break;
      C = addii(b, N); if (aN > 1) C = diviuexact(C, aN);
      c = sa > 0? addii(c, C): subii(c, C);
    }
  }
  return L;
}

/* Let M and N in SL2(Z), return (N*M^-1)[,1] */
static GEN
SL2_div_mul_e1(GEN N, GEN M)
{
  GEN b = gcoeff(M,2,1), d = gcoeff(M,2,2);
  GEN A = mulii(gcoeff(N,1,1), d), B = mulii(gcoeff(N,1,2), b);
  GEN C = mulii(gcoeff(N,2,1), d), D = mulii(gcoeff(N,2,2), b);
  retmkvec2(subii(A,B), subii(C,D));
}
static GEN
qfisolve_normform(GEN Q, GEN P)
{
  GEN a = gel(Q,1), N = gel(Q,2);
  GEN M, b = redimagsl2(P, &M);
  if (!qfb_equal(a,b)) return NULL;
  return SL2_div_mul_e1(N,M);
}

/* Test equality modulo GL2 of two reduced forms */
static int
GL2_qfb_equal(GEN a, GEN b)
{
  return equalii(gel(a,1),gel(b,1))
   && absequalii(gel(a,2),gel(b,2))
   &&    equalii(gel(a,3),gel(b,3));
}

/* Q(u,v) = p; if s < 0 return that solution; else the set of all solutions */
static GEN
allsols(GEN Q, long s, GEN u, GEN v)
{
  GEN w = mkvec2(u, v), b;
  if (signe(v) < 0) { u = negi(u); v = negi(v); } /* normalize for v >= 0 */
  w = mkvec2(u, v); if (s < 0) return w;
  if (!s) return mkvec(w);
  b = gel(Q,2); /* sum of the 2 solutions (if they exist) is -bv / a */
  if (signe(b))
  { /* something to check */
    GEN r, t;
    t = dvmdii(mulii(b, v), gel(Q,1), &r);
    if (signe(r)) return mkvec(w);
    u = addii(u, t);
  }
  return mkvec2(w, mkvec2(negi(u), v));
}
static GEN
qfisolvep_all(GEN Q, GEN p, long all)
{
  GEN R, U, V, M, N, x, q, D = qfb_disc(Q);
  long s = kronecker(D, p);

  if (s < 0) return NULL;
  if (!all) s = -1; /* to indicate we want a single solution */
  /* Solutions iff a class of maximal ideal above p is the class of Q;
   * Two solutions iff (s > 0 and the class has order > 2), else one */
  if (!signe(gel(Q,2)))
  { /* if principal form, use faster cornacchia */
    GEN a = gel(Q,1), c = gel(Q,3);
    if (equali1(a))
    {
      if (!cornacchia(c, p, &M,&N)) return NULL;
      return allsols(Q, s, M, N);
    }
    if (equali1(c))
    {
      if (!cornacchia(a, p, &M,&N)) return NULL;
      return allsols(Q, s, N, M);
    }
  }
  R = redimagsl2(Q, &U);
  if (equali1(gel(R,1)))
  { /* principal form */
    if (!signe(gel(R,2)))
    {
      if (!cornacchia(gel(R,3), p, &M,&N)) return NULL;
      x = mkvec2(M,N);
    }
    else
    { /* x^2 + xy + ((1-D)/4)y^2 = p <==> (2x + y)^2 - D y^2 = 4p */
      if (!cornacchia2(negi(D), p, &M, &N)) return NULL;
      x = subii(M,N); if (mpodd(x)) return NULL;
      x = mkvec2(shifti(x,-1), N);
    }
    x = ZM_ZC_mul(U, x); x[0] = evaltyp(t_VEC) | _evallg(3); /* transpose */
    return allsols(Q, s, gel(x,1), gel(x,2));
  }
  q = redimagsl2(primeform(D, p), &V);
  if (!GL2_qfb_equal(R,q)) return NULL;
  if (signe(gel(R,2)) != signe(gel(q,2))) gcoeff(V,2,1) = negi(gcoeff(V,2,1));
  x = SL2_div_mul_e1(U,V); return allsols(Q, s, gel(x,1), gel(x,2));
}
GEN
qfisolvep(GEN Q, GEN p)
{
  pari_sp av = avma;
  GEN x = qfisolvep_all(Q, p, 0);
  return x? gerepilecopy(av, x): gc_const(av, gen_0);
}

static GEN
qfrsolve_normform(GEN N, GEN Ps, GEN rd)
{
  pari_sp av = avma, btop;
  GEN M = N, P = redrealsl2(Ps, rd), Q = P;

  btop = avma;
  for(;;)
  {
    if (qfb_equal(gel(M,1), gel(P,1)))
      return gerepileupto(av, SL2_div_mul_e1(gel(M,2),gel(P,2)));
    if (qfb_equal(gel(N,1), gel(Q,1)))
      return gerepileupto(av, SL2_div_mul_e1(gel(N,2),gel(Q,2)));
    M = rhorealsl2(M, rd);
    if (qfb_equal(gel(M,1), gel(N,1))) return gc_NULL(av);
    Q = rhorealsl2(Q, rd);
    if (qfb_equal(gel(P,1), gel(Q,1))) return gc_NULL(av);
    if (gc_needed(btop, 1)) gerepileall(btop, 2, &M, &Q);
  }
}

GEN
qfrsolvep(GEN Q, GEN p)
{
  pari_sp av = avma;
  GEN N, x, rd, d = qfb_disc(Q);

  if (kronecker(d, p) < 0) return gc_const(av, gen_0);
  rd = sqrti(d);
  N = redrealsl2(Q, rd);
  x = qfrsolve_normform(N, primeform(d, p), rd);
  return x? gerepileupto(av, x): gc_const(av, gen_0);
}

static GEN
known_prime(GEN v)
{
  GEN p, e, fa = check_arith_all(v, "qfbsolve");
  if (!fa) return BPSW_psp(v)? v: NULL;
  if (lg(gel(fa,1)) != 2) return NULL;
  p = gcoeff(fa,1,1);
  e = gcoeff(fa,1,2);
  return (equali1(e) && !is_pm1(p) && signe(p) > 0)? p: NULL;
}
static GEN
qfsolve_normform(GEN Q, GEN f, GEN rd)
{ return rd? qfrsolve_normform(Q, f, rd): qfisolve_normform(Q, f); }
static GEN
qfbsolve_primitive_i(GEN Q, GEN rd, GEN *Qr, GEN fa, long all)
{
  GEN x, W, F, p;
  long i, j, l;
  if (!rd && (p = known_prime(fa))) return qfisolvep_all(Q, p, all);
  F = normforms(qfb_disc(Q), fa);
  if (!F) return NULL;
  if (!*Qr) *Qr = qfbredsl2(Q, rd);
  l = lg(F); W = all? cgetg(l, t_VEC): NULL;
  for (j = i = 1; i < l; i++)
    if ((x = qfsolve_normform(*Qr, gel(F,i), rd)))
    {
      if (!all) return x;
      gel(W,j++) = x;
    }
  if (j == 1) return NULL;
  setlg(W,j); return W;
}

static GEN
qfb_initrd(GEN Q) { GEN d = qfb_disc(Q); return signe(d) > 0? sqrti(d): NULL; }
static GEN
qfbsolve_primitive(GEN Q, GEN fa, long all)
{
  GEN x, Qr = NULL, rdQ = qfb_initrd(Q);
  x = qfbsolve_primitive_i(Q, rdQ, &Qr, fa, all);
  if (!x) return cgetg(1, t_VEC);
  return x;
}

/* f / g^2 */
static GEN
famat_divsqr(GEN f, GEN g)
{ return famat_reduce(famat_div_shallow(f, famat_pows_shallow(g,2))); }
static GEN
qfbsolve_all(GEN Q, GEN n, long all)
{
  GEN W, Qr = NULL, fa = factorint(n, 0), rdQ = qfb_initrd(Q);
  GEN D = divisors_factored(mkmat2(gel(fa,1), gshift(gel(fa,2),-1)));
  long i, j, l = lg(D);
  W = all? cgetg(l, t_VEC): NULL;
  for (i = j = 1; i < l; i++)
  {
    GEN w, d = gel(D,i), FA = i == 1? fa: famat_divsqr(fa, gel(d,2));
    if ((w = qfbsolve_primitive_i(Q, rdQ, &Qr, FA, all)))
    {
      if (i != 1) w = RgV_Rg_mul(w, gel(d,1));
      if (!all) return w;
      gel(W,j++) = w;
    }
  }
  if (j == 1) return cgetg(1, t_VEC);
  setlg(W,j); return shallowconcat1(W);
}

GEN
qfbsolve(GEN Q, GEN n, long fl)
{
  pari_sp av = avma;
  if (typ(Q) != t_QFB) pari_err_TYPE("qfbsolve",Q);
  if (fl < 0 || fl > 3) pari_err_FLAG("qfbsolve");
  return gerepilecopy(av, (fl & 2)? qfbsolve_all(Q, n, fl & 1)
                                  : qfbsolve_primitive(Q, n, fl & 1));
}

/* 1 if there exists x,y such that x^2 + dy^2 = p, 0 otherwise;
 * Assume d > 0 and p is prime */
long
cornacchia(GEN d, GEN p, GEN *px, GEN *py)
{
  pari_sp av = avma;
  GEN b, c, r;

  *px = *py = gen_0;
  b = subii(p, d);
  if (signe(b) < 0) return gc_long(av,0);
  if (signe(b) == 0) { *py = gen_1; return gc_long(av,1); }
  b = Fp_sqrt(b, p); /* sqrt(-d) */
  if (!b) return gc_long(av,0);
  b = gmael(halfgcdii(p, b), 2, 2);
  c = dvmdii(subii(p, sqri(b)), d, &r);
  if (r != gen_0 || !Z_issquareall(c, &c)) return gc_long(av,0);
  set_avma(av);
  *px = icopy(b);
  *py = icopy(c); return 1;
}

static GEN
lastqi(GEN Q)
{
  GEN s = gcoeff(Q,1,1), q = gcoeff(Q,1,2), p = absi_shallow(gcoeff(Q,2,2));
  if (!signe(q)) return gen_0;
  if (!signe(s)) return p;
  if (is_pm1(q)) return subiu(p,1);
  return divii(p, absi_shallow(q));
}

static long
cornacchia2_i(long av, GEN d, GEN p, GEN b, GEN px4, GEN *px, GEN *py)
{
  GEN M, Q, V, c, r, b2;
  if (!signe(b)) { /* d = p,2p,3p,4p */
    set_avma(av);
    if (absequalii(d, px4)){ *py = gen_1; return 1; }
    if (absequalii(d, p))  { *py = gen_2; return 1; }
    return 0;
  }
  if (mod2(b) != mod2(d)) b = subii(p,b);
  M = halfgcdii(shifti(p,1), b); Q = gel(M,1); V = gel(M, 2);
  b = addii(mulii(gel(V,1), lastqi(Q)), gel(V,2));
  b2 = sqri(b);
  if (cmpii(b2,px4) > 0)
  {
    b = gel(V,1); b2 = sqri(b);
    if (cmpii(b2,px4) > 0) { b = gel(V,2); b2 = sqri(b); }
  }
  c = dvmdii(subii(px4, b2), d, &r);
  if (r != gen_0 || !Z_issquareall(c, &c)) return gc_long(av,0);
  set_avma(av);
  *px = icopy(b);
  *py = icopy(c); return 1;
}

/* 1 if there exists x,y such that x^2 + dy^2 = 4p, 0 otherwise;
 * Assume d > 0 is congruent to 0 or 3 mod 4 and p is prime */
long
cornacchia2(GEN d, GEN p, GEN *px, GEN *py)
{
  pari_sp av = avma;
  GEN b, p4 = shifti(p,2);

  *px = *py = gen_0;
  if (abscmpii(p4, d) < 0) return gc_long(av,0);
  if (absequaliu(p, 2))
  {
    set_avma(av);
    switch (itou_or_0(d)) {
      case 4: *px = gen_2; break;
      case 7: *px = gen_1; break;
      default: return 0;
    }
    *py = gen_1; return 1;
  }
  b = Fp_sqrt(negi(d), p);
  if (!b) return gc_long(av,0);
  return cornacchia2_i(av, d, p, b, p4, px, py);
}

/* 1 if there exists x,y such that x^2 + dy^2 = 4p [p prime], 0 otherwise */
long
cornacchia2_sqrt(GEN d, GEN p, GEN b, GEN *px, GEN *py)
{
  pari_sp av = avma;
  GEN p4 = shifti(p,2);
  *px = *py = gen_0;
  if (abscmpii(p4, d) < 0) return gc_long(av,0);
  return cornacchia2_i(av, d, p, b, p4, px, py);
}

GEN
qfbcornacchia(GEN d, GEN p)
{
  pari_sp av = avma;
  GEN x, y;
  if (typ(d) != t_INT || signe(d) <= 0) pari_err_TYPE("qfbcornacchia", d);
  if (typ(p) != t_INT || cmpiu(p, 2) < 0) pari_err_TYPE("qfbcornacchia", p);
  if (mod4(p)? cornacchia(d, p, &x, &y): cornacchia2(d, shifti(p, -2), &x, &y))
    return gerepilecopy(av, mkvec2(x, y));
  set_avma(av); return cgetg(1, t_VEC);
}
