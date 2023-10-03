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

/********************************************************************/
/**                                                                **/
/**                       ELLIPTIC CURVES                          **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_ell

#undef coordch

/* Transforms a curve E into short Weierstrass form E' modulo p.
   Returns a vector, the first two entries of which are a4' and a6'.
   The third entry is a vector describing the isomorphism E' \to E.
*/

static ulong
Fl_c4_to_a4(ulong c4, ulong p)
{ return Fl_neg(Fl_mul(c4, 27, p), p); }
static ulong
Fl_c6_to_a6(ulong c6, ulong p)
{ return Fl_neg(Fl_mul(c6, 54, p), p); }
static void
Fl_c4c6_to_a4a6(ulong c4, ulong c6, ulong p, ulong *a4, ulong *a6)
{
  *a4 = Fl_c4_to_a4(c4, p);
  *a6 = Fl_c6_to_a6(c6, p);
}
static GEN
c4_to_a4(GEN c4, GEN p)
{ return Fp_neg(Fp_mulu(c4, 27, p), p); }
static void
c4c6_to_a4a6(GEN c4, GEN c6, GEN p, GEN *a4, GEN *a6)
{
  *a4 = c4_to_a4(c4, p);
  *a6 = Fp_neg(Fp_mulu(c6, 54, p), p);
}
static GEN
Fq_c4_to_a4(GEN c4, GEN T, GEN p)
{ return Fq_neg(Fq_mulu(c4, 27, T,p), T,p); }
static void
Fq_c4c6_to_a4a6(GEN c4, GEN c6, GEN T, GEN p, GEN *a4, GEN *a6)
{
  *a4 = Fq_c4_to_a4(c4, T,p);
  *a6 = Fq_neg(Fq_mulu(c6, 54, T,p), T,p);
}
static void
ell_to_a4a6(GEN E, GEN p, GEN *a4, GEN *a6)
{
  GEN c4 = Rg_to_Fp(ell_get_c4(E),p);
  GEN c6 = Rg_to_Fp(ell_get_c6(E),p);
  c4c6_to_a4a6(c4, c6, p, a4, a6);
}
static void
Fl_ell_to_a4a6(GEN E, ulong p, ulong *a4, ulong *a6)
{
  ulong c4 = Rg_to_Fl(ell_get_c4(E),p);
  ulong c6 = Rg_to_Fl(ell_get_c6(E),p);
  Fl_c4c6_to_a4a6(c4, c6, p, a4, a6);
}

/* [6,3b2,3a1,108a3] */
static GEN
a4a6_ch(GEN E, GEN p)
{
  GEN a1 = Rg_to_Fp(ell_get_a1(E),p);
  GEN a3 = Rg_to_Fp(ell_get_a3(E),p);
  GEN b2 = Rg_to_Fp(ell_get_b2(E),p);
  retmkvec4(modsi(6,p),Fp_mulu(b2,3,p),Fp_mulu(a1,3,p),Fp_mulu(a3,108,p));
}
static GEN
a4a6_ch_Fl(GEN E, ulong p)
{
  ulong a1 = Rg_to_Fl(ell_get_a1(E),p);
  ulong a3 = Rg_to_Fl(ell_get_a3(E),p);
  ulong b2 = Rg_to_Fl(ell_get_b2(E),p);
  return mkvecsmall4(6 % p,Fl_mul(b2,3,p),Fl_mul(a1,3,p),Fl_mul(a3,108,p));
}

static GEN
ell_to_a4a6_bc(GEN E, GEN p)
{
  GEN A4, A6;
  ell_to_a4a6(E, p, &A4, &A6);
  retmkvec3(A4, A6, a4a6_ch(E,p));
}
GEN
point_to_a4a6(GEN E, GEN P, GEN p, GEN *pa4)
{
  GEN c4 = Rg_to_Fp(ell_get_c4(E),p);
  *pa4 = c4_to_a4(c4, p);
  return FpE_changepointinv(RgV_to_FpV(P,p), a4a6_ch(E,p), p);
}
GEN
point_to_a4a6_Fl(GEN E, GEN P, ulong p, ulong *pa4)
{
  ulong c4 = Rg_to_Fl(ell_get_c4(E),p);
  *pa4 = Fl_c4_to_a4(c4, p);
  return Fle_changepointinv(RgV_to_Flv(P,p), a4a6_ch_Fl(E,p), p);
}

/* shallow basistoalg */
static GEN
nftoalg(GEN nf, GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_FRAC: case t_POLMOD: return x;
    default: return basistoalg(nf, x);
  }
}

void
checkellpt(GEN z)
{
  if (typ(z)!=t_VEC) pari_err_TYPE("checkellpt", z);
  switch(lg(z))
  {
    case 3: break;
    case 2: if (isintzero(gel(z,1))) break;
    /* fall through */
    default: pari_err_TYPE("checkellpt", z);
  }
}
void
checkell5(GEN E)
{
  long l = lg(E);
  if (typ(E)!=t_VEC || (l != 17 && l != 6)) pari_err_TYPE("checkell5",E);
}
void
checkell(GEN E)
{ if (!checkell_i(E)) pari_err_TYPE("checkell",E); }
void
checkellisog(GEN v)
{ if (typ(v)!=t_VEC || lg(v) != 4) pari_err_TYPE("checkellisog",v); }

void
checkell_Q(GEN E)
{
  if (!checkell_i(E) || ell_get_type(E)!=t_ELL_Q)
    pari_err_TYPE("checkell over Q",E);
}

void
checkell_Qp(GEN E)
{
  if (!checkell_i(E) || ell_get_type(E)!=t_ELL_Qp)
    pari_err_TYPE("checkell over Qp",E);
}

static int
ell_over_Fq(GEN E)
{
  long t = ell_get_type(E);
  return t==t_ELL_Fp || t==t_ELL_Fq;
}

void
checkell_Fq(GEN E)
{
  if (!checkell_i(E) || !ell_over_Fq(E)) pari_err_TYPE("checkell over Fq", E);
}

GEN
ellff_get_p(GEN E)
{
  GEN fg = ellff_get_field(E);
  return typ(fg)==t_INT? fg: FF_p_i(fg);
}

int
ell_is_integral(GEN E)
{
  return typ(ell_get_a1(E)) == t_INT
      && typ(ell_get_a2(E)) == t_INT
      && typ(ell_get_a3(E)) == t_INT
      && typ(ell_get_a4(E)) == t_INT
      && typ(ell_get_a6(E)) == t_INT;
}

static void
checkcoordch(GEN z)
{ if (typ(z)!=t_VEC || lg(z) != 5) pari_err_TYPE("checkcoordch",z); }

/* 4 X^3 + b2 X^2 + 2b4 X + b6 */
GEN
ec_bmodel(GEN e, long v)
{
  GEN b2 = ell_get_b2(e), b6 = ell_get_b6(e), b42 = gmul2n(ell_get_b4(e),1);
  GEN P = mkpoln(4, utoipos(4), b2, b42, b6);
  setvarn(P, v); return P;
}

/* X^4 - b4*X^2 - 2b6*X - b8 */
GEN
ec_phi2(GEN e, long v)
{
  GEN b4 = ell_get_b4(e), b6 = ell_get_b6(e), b8 = ell_get_b8(e);
  GEN P = mkpoln(5, gen_1, gen_0, gneg(b4), gmul2n(gneg(b6),1), gneg(b8));
  setvarn(P, v); return P;
}

static int
invcmp(void *E, GEN x, GEN y) { (void)E; return -gcmp(x,y); }

/* prec = working precision, prec0 = target precision */
static GEN
doellR_roots_i(GEN e, long prec, long prec0)
{
  GEN d1, d2, d3, e1, e2, e3, R = cleanroots(ec_bmodel(e,0), prec);
  long s = ellR_get_sign(e);
  if (s > 0)
  { /* sort 3 real roots in decreasing order */
    R = real_i(R);
    gen_sort_inplace(R, NULL, &invcmp, NULL);
    e1 = gel(R,1); e2 = gel(R,2); e3 = gel(R,3);
    d3 = subrr(e1,e2);
    d1 = subrr(e2,e3);
    d2 = subrr(e1,e3);
    if (realprec(d3) < prec0 || realprec(d1) < prec0) return NULL;
  } else {
    e1 = gel(R,1); e2 = gel(R,2); e3 = gel(R,3);
    if (s < 0)
    { /* make sure e1 is real, imag(e2) > 0 and imag(e3) < 0 */
      e1 = real_i(e1);
      if (signe(gel(e2,2)) < 0) swap(e2, e3);
      d1 = mkcomplex(gen_0, gsub(gel(e2,2),gel(e3,2)));
    }
    else
      d1 = gsub(e2,e3);
    d3 = gsub(e1,e2);
    d2 = gsub(e1,e3);
    if (precision(d1) < prec0
        || precision(d2) < prec0
        || precision(d3) < prec0) return NULL;
  }
  return mkcol6(e1,e2,e3,d1,d2,d3);
}
static GEN
doellR_roots(GEN e, long prec0)
{
  long p;
  for (p = prec0;; p = precdbl(p))
  {
    GEN v = doellR_roots_i(e, p, prec0);
    if (v) return v;
    if (DEBUGLEVEL) pari_warn(warnprec,"doellR_roots", p);
  }
}
static GEN
ellR_root(GEN e, long prec) { return gel(ellR_roots(e,prec),1); }

/* Given E and the x-coordinate of a point Q = [xQ, yQ], return
 *   f(xQ) = xQ^3 + E.a2 * xQ^2 + E.a4 * xQ + E.a6
 * where E is given by y^2 + h(x)y = f(x). */
GEN
ec_f_evalx(GEN E, GEN x)
{
  pari_sp av = avma;
  GEN z;
  z = gadd(ell_get_a2(E),x);
  z = gadd(ell_get_a4(E), gmul(x,z));
  z = gadd(ell_get_a6(E), gmul(x,z));
  return gerepileupto(av, z); /* ((x + E.a2) * x + E.a4) * x + E.a6 */
}

/* a1 x + a3 */
GEN
ec_h_evalx(GEN e, GEN x)
{
  GEN a1 = ell_get_a1(e);
  GEN a3 = ell_get_a3(e);
  return gadd(a3, gmul(x,a1));
}
static GEN
Zec_h_evalx(GEN e, GEN x)
{
  GEN a1 = ell_get_a1(e);
  GEN a3 = ell_get_a3(e);
  return signe(a1)? addii(a3, mulii(x, a1)): a3;
}
/* y^2 + a1 xy + a3 y = y^2 + h(x)y */
static GEN
ec_LHS_evalQ(GEN e, GEN Q)
{
  GEN x = gel(Q,1), y = gel(Q,2);
  return gmul(y, gadd(y, ec_h_evalx(e,x)));
}

/* Given E and a point Q = [xQ, yQ], return
 *   3 * xQ^2 + 2 * E.a2 * xQ + E.a4 - E.a1 * yQ.
 * which is the derivative of the curve equation
 *   f(x) - (y^2 + h(x)y) = 0
 * wrt x evaluated at Q */
GEN
ec_dFdx_evalQ(GEN E, GEN Q)
{
  pari_sp av = avma;
  GEN x = gel(Q,1), y = gel(Q,2);
  GEN a1 = ell_get_a1(E);
  GEN a2 = ell_get_a2(E);
  GEN a4 = ell_get_a4(E);
  GEN tmp = gmul(gadd(gmulsg(3L,x), gmul2n(a2,1)), x);
  return gerepileupto(av, gadd(tmp, gsub(a4, gmul(a1, y))));
}

/* 2y + a1 x + a3 = -ec_dFdy_evalQ */
GEN
ec_dmFdy_evalQ(GEN e, GEN Q)
{
  GEN x = gel(Q,1), y = gel(Q,2);
  return gadd(ec_h_evalx(e,x), gmul2n(y,1));
}
/* Given E and a point Q = [xQ, yQ], return
 *  -(2 * yQ + E.a1 * xQ + E.a3)
 * which is the derivative of the curve equation
 *  f(x) - (y^2 + h(x)y) = 0
 * wrt y evaluated at Q */
GEN
ec_dFdy_evalQ(GEN E, GEN Q)
{
  pari_sp av = avma;
  return gerepileupto(av, gneg(ec_dmFdy_evalQ(E,Q)));
}

/* Given E and a point Q = [xQ, yQ], return
 *   4 xQ^3 + E.b2 xQ^2 + 2 E.b4 xQ + E.b6
 * which is the 2-division polynomial of E evaluated at Q */
GEN
ec_2divpol_evalx(GEN E, GEN x)
{
  pari_sp av = avma;
  GEN b2 = ell_get_b2(E), x4 = gmul2n(x,2), t1, t2;
  GEN b42 = gmul2n(ell_get_b4(E), 1);
  GEN b6 = ell_get_b6(E);
  if (ell_get_type(E) == t_ELL_NF && typ(x)==t_COL)
  {
    GEN nf = ellnf_get_nf(E);
    t1 = nfmul(nf, nfadd(nf, x4, b2), x);
    t2 = nfadd(nf, t1, b42);
    t2 = nfadd(nf, nfmul(nf, t2, x), b6);
    t2 = nftoalg(nf, t2);
  }
  else
  {
    t1 = gmul(gadd(x4, b2), x);
    t2 = gadd(t1, b42);
    t2 = gadd(gmul(t2, x), b6);
  }
  return gerepileupto(av, t2);
}

/* Given E and a point Q = [xQ, yQ], return
 *   3 xQ^4 + E.b2 xQ^3 + 3 E.b4 xQ^2 + 3*E.b6 xQ + E.b8
 * which is the 3-division polynomial of E evaluated at Q */
GEN
ec_3divpol_evalx(GEN E, GEN x)
{
  pari_sp av = avma;
  GEN b2 = ell_get_b2(E);
  GEN b4 = ell_get_b4(E);
  GEN b6 = ell_get_b6(E);
  GEN b8 = ell_get_b8(E);
  GEN x2 = gsqr(x);
  GEN t1 = gadd(gadd(gmulsg(3L, x2), gmul(b2, x)), gmulsg(3L, b4));
  GEN t2 = gadd(gmul(gmulsg(3L, b6), x), b8);
  return gerepileupto(av, gadd(gmul(t1, x2), t2));
}

/* Given E and a point Q = [xQ, yQ], return
 *   6 xQ^2 + E.b2 xQ + E.b4
 * which, if f is the curve equation, is 2 dfdx - E.a1 dfdy evaluated at Q */
GEN
ec_half_deriv_2divpol(GEN E, long v)
{ return deg2pol_shallow(utoi(6), ell_get_b2(E), ell_get_b4(E), v); }

GEN
ec_half_deriv_2divpol_evalx(GEN E, GEN x)
{
  pari_sp av = avma;
  GEN b2 = ell_get_b2(E);
  GEN b4 = ell_get_b4(E);
  GEN res = gadd(gmul(gadd(gmulsg(6L, x), b2), x), b4);
  return gerepileupto(av, res);
}

/* Return the characteristic of the ring over which E is defined. */
GEN
ellbasechar(GEN E)
{
  pari_sp av = avma;
  GEN D = ell_get_disc(E);
  return gerepileuptoint(av, characteristic(D));
}

/* Initialize basic elliptic struct y[1..12] for initsmall
 * (do not include j to allow for singular Weistrass model)
 * Also allocate room for n dynamic members. */
static GEN
initsmall46(GEN a4, GEN a6, long n)
{
  GEN y = obj_init(15, n);
  gel(y,1) = gen_0;
  gel(y,2) = gen_0;
  gel(y,3) = gen_0;
  gel(y,4) = a4;
  gel(y,5) = a6;
  gel(y,6) = gen_0;
  gel(y,7) = gmul2n(a4,1);
  gel(y,8) = gmul2n(a6,2);
  gel(y,9) = gneg(gsqr(a4));
  gel(y,10)= gmulgs(a4,-48);
  gel(y,11)= gmulgs(a6,-864);
  gel(y,12)= gadd(gmul(gmulgs(a4,-64), gsqr(a4)), gmulsg(-432,gsqr(a6)));
  gel(y,16) = zerovec(n); return y;
}
/* [a1,a2,a3,a4,a6] */
static GEN
initsmall5(GEN x, long n)
{
  GEN a1 = gel(x,1), a2 = gel(x,2), a3 = gel(x,3);
  GEN a4 = gel(x,4), a6 = gel(x,5);
  GEN y, b2, b4, b6, b8, c4, c6, D, a11, a13, a33, b22;
  if (gequal0(a1) && gequal0(a2) && gequal0(a3)) return initsmall46(a4, a6, n);
  a11= gsqr(a1);
  b2 = gadd(a11, gmul2n(a2,2));
  a13= gmul(a1, a3);
  b4 = gadd(a13, gmul2n(a4,1));
  a33= gsqr(a3);
  b6 = gadd(a33, gmul2n(a6,2));
  b8 = gsub(gadd(gmul(a11,a6), gmul(b6, a2)), gmul(a4, gadd(a4,a13)));
  b22= gsqr(b2);
  c4 = gadd(b22, gmulsg(-24,b4));
  c6 = gadd(gmul(b2,gsub(gmulsg(36,b4),b22)), gmulsg(-216,b6));
  D  = gsub(gmul(b4, gadd(gmulsg(9,gmul(b2,b6)),gmulsg(-8,gsqr(b4)))),
            gadd(gmul(b22,b8),gmulsg(27,gsqr(b6))));
  y = obj_init(15, n);
  gel(y,1) = a1;
  gel(y,2) = a2;
  gel(y,3) = a3;
  gel(y,4) = a4;
  gel(y,5) = a6;
  gel(y,6) = b2; /* a1^2 + 4a2 */
  gel(y,7) = b4; /* a1 a3 + 2a4 */
  gel(y,8) = b6; /* a3^2 + 4 a6 */
  gel(y,9) = b8; /* a1^2 a6 + 4a6 a2 + a2 a3^2 - a4(a4 + a1 a3) */
  gel(y,10)= c4; /* b2^2 - 24 b4 */
  gel(y,11)= c6; /* 36 b2 b4 - b2^3 - 216 b6 */
  gel(y,12)= D;
  gel(y,16) = zerovec(n); return y;
}

static GEN
get_j(GEN c4, GEN D)
{
  GEN g, d, c;
  if (typ(D) != t_POL || typ(c4) != t_POL || varn(D) != varn(c4))
    return gdiv(gmul(gsqr(c4),c4), D);
  /* c4^3 / D, simplifying incrementally */
  g = RgX_gcd(D, c4);
  if (degpol(g) == 0) return gred_rfrac_simple(gmul(gsqr(c4),c4), D);
  c = RgX_div(c4, g);
  D = RgX_div(D, g);
  g = RgX_gcd(D,c4);
  if (degpol(g) == 0) return gred_rfrac_simple(gmul(gsqr(c4),c), D);
  D = RgX_div(D, g); d = RgX_div(c4, g);
  g = RgX_gcd(D,c4);
  if (degpol(g)) { D = RgX_div(D, g); c4 = RgX_div(c4, g); }
  return gred_rfrac_simple(gmul(gmul(c4, d),c), D);
}

/* return basic elliptic struct y[1..13], y[14] (domain type) and y[15]
 * (domain-specific data) are left uninitialized, from x[1], ..., x[5].
 * Also allocate room for n dynamic members (actually stored in the last
 * component y[16])*/
static GEN
initsmall(GEN x, long n)
{
  GEN y, D;

  switch(lg(x))
  {
    case 2: y = initsmall5(ellfromj(gel(x,1)), n); break;
    case 3: y = initsmall46(gel(x,1), gel(x,2), n); break;
    case 6:
    case 17: y = initsmall5(x, n); break;
    default:
      pari_err_TYPE("ellxxx [not an elliptic curve (ell5)]",x);
      return NULL; /* LCOV_EXCL_LINE */
  }
  D = ell_get_disc(y); if (gequal0(D)) return NULL;
  gel(y,13) = get_j(ell_get_c4(y), D); return y;
}
void
ellprint(GEN e)
{
  pari_sp av = avma;
  long vx, vy;
  GEN z;
  checkell5(e);
  vx = fetch_var(); name_var(vx, "X");
  vy = fetch_var(); name_var(vy, "Y"); z = mkvec2(pol_x(vx), pol_x(vy));
  err_printf("%Ps - (%Ps)\n", ec_LHS_evalQ(e, z), ec_f_evalx(e, pol_x(vx)));
  (void)delete_var();
  (void)delete_var(); set_avma(av);
}

/* compute a,b such that E1: y^2 = x(x-a)(x-b) ~ E */
static GEN
doellR_ab(GEN E, long prec)
{
  GEN b2 = ell_get_b2(E), R = ellR_roots(E, prec);
  GEN e1 = gel(R,1), d2 = gel(R,5), d3 =  gel(R,6), a, b, t;

  t = gmul2n(gadd(mulur(12,e1), b2), -4); /* = (12 e1 + b2) / 16 */
  if (ellR_get_sign(E) > 0)
    b = mulrr(d3,d2);
  else
    b = cxnorm(d3);
  b = sqrtr(b); /* = sqrt( (e1 - e2)(e1 - e3) ) */
  if (gsigne(t) > 0) togglesign(b);
  a = gsub(gmul2n(b,-1),t);
  return mkvec2(a, b);
}
GEN
ellR_ab(GEN E, long prec)
{ return obj_checkbuild_realprec(E, R_AB, &doellR_ab, prec); }

/* q a t_REAL*/
static long
real_prec(GEN q)
{ return signe(q)? realprec(q): LONG_MAX; }
/* q a t_PADIC */
static long
padic_prec(GEN q)
{ return signe(gel(q,4))? precp(q)+valp(q): valp(q); }

/* check whether moduli are consistent */
static void
chk_p(GEN p, GEN p2)
{ if (!equalii(p, p2)) pari_err_MODULUS("ellinit", p,p2); }

static int
fix_nftype(GEN *pp)
{
  switch(nftyp(*pp))
  {
    case typ_NF: case typ_BNF: break;
    case typ_BNR:*pp = bnr_get_bnf(*pp); break;
    default: return 0;
  }
  return 1;
}
static long
base_ring(GEN x, GEN *pp, long *prec)
{
  long i, e = *prec, ep = LONG_MAX, imax = minss(lg(x), 6);
  GEN p = NULL, pol = NULL;
  long t = t_FRAC;
  if (*pp) switch(t = typ(*pp))
  {
    case t_INT:
      if (cmpis(*pp,2) < 0) { t = t_FRAC; p = NULL; break; }
      p = *pp;
      t = t_INTMOD;
      break;
    case t_INTMOD:
      p = gel(*pp, 1);
      break;
    case t_REAL:
      e = real_prec(*pp);
      p = NULL;
      break;
    case t_PADIC:
      ep = padic_prec(*pp);
      p = gel(*pp, 2);
      break;
    case t_FFELT:
      p = *pp;
      break;
    case t_VEC:
      t = t_VEC; p = *pp;
      if (fix_nftype(&p)) break;
    default:
      pari_err_TYPE("elliptic curve base_ring", *pp);
      return 0;
  }
  if (t==t_VEC) pol = nf_get_pol(checknf(p));
  /* Possible cases:
   * t = t_VEC (p an nf or bnf)
   * t = t_FFELT (p t_FFELT)
   * t = t_INTMOD (p a prime)
   * t = t_PADIC (p a prime, ep = padic prec)
   * t = t_REAL (p = NULL, e = real prec)
   * t = t_FRAC (p = NULL) */
  for (i = 1; i < imax; i++)
  {
    GEN p2, q = gel(x,i);
    switch(typ(q)) {
      case t_PADIC:
        p2 = gel(q,2);
        switch(t)
        {
          case t_FRAC:  t = t_PADIC; p = p2; break;
          case t_PADIC: chk_p(p,p2); break;
          default: pari_err_TYPE("elliptic curve base_ring", x);
        }
        ep = minss(ep, padic_prec(q));
        break;
      case t_INTMOD:
        p2 = gel(q,1);
        switch(t)
        {
          case t_FRAC:  t = t_INTMOD; p = p2; break;
          case t_FFELT: chk_p(FF_p_i(p),p2); break;
          case t_INTMOD:chk_p(p,p2); break;
          default: pari_err_TYPE("elliptic curve base_ring", x);
        }
        break;
      case t_FFELT:
        switch(t)
        {
          case t_INTMOD: chk_p(p, FF_p_i(q)); /* fall through */
          case t_FRAC:   t = t_FFELT; p = q; break;
          case t_FFELT:
            if (!FF_samefield(p,q) && FF_f(q)>1) pari_err_MODULUS("ellinit", p,q);
            break;
          default: pari_err_TYPE("elliptic curve base_ring", x);
        }
        break;

      case t_INT: case t_FRAC: break;
      case t_REAL:
        switch(t)
        {
          case t_REAL: e = minss(e, real_prec(q)); break;
          case t_FRAC: e = real_prec(q); t = t_REAL; break;
          default: pari_err_TYPE("elliptic curve base_ring", x);
        }
        break;
      case t_POLMOD:
        if (pol && !RgX_equal(pol, gel(q,1)))
          pari_err_MODULUS("ellinit",gel(q,1), pol);
      case t_COL:
      case t_POL:
        if (t == t_VEC) break;
      default: /* base ring too general */
        return t_COMPLEX;
    }
  }
  *pp = p; *prec = (t == t_PADIC)? ep: e; return t;
}

/* s = 0 complex, else real;
 * if (s = 2) set s = sign(D), else accept s as is */
static GEN
ellinit_Rg(GEN x, long s, long prec)
{
  GEN y;
  if (lg(x) > 6) switch(ell_get_type(x))
  {
    case t_ELL_Rg:
    case t_ELL_Q: break;
    default: pari_err_TYPE("elliptic curve base_ring", x);
  }
  if (!(y = initsmall(x, 4))) return NULL;
  if (s == 2) s = gsigne(ell_get_disc(y));
  gel(y,14) = mkvecsmall(t_ELL_Rg);
  gel(y,15) = mkvec(mkvecsmall2(prec2nbits(prec), s));
  return y;
}

static GEN
ellinit_Qp(GEN x, GEN p, long prec)
{
  GEN y;
  if (lg(x) > 6)
  {
    switch(ell_get_type(x))
    { /* sanity checks */
      case t_ELL_Q: break;
      case t_ELL_Qp: chk_p(ellQp_get_p(x), p); break;
      default: pari_err_TYPE("elliptic curve base_ring", x);
    }
    x = vecslice(x,1,5);
  }
  x = QpV_to_QV(x); /* make entries rational */
  if (!(y = initsmall(x, 2))) return NULL;
  gel(y,14) = mkvecsmall(t_ELL_Qp);
  gel(y,15) = mkvec(zeropadic(p, prec));
  return y;
}

static GEN
ellinit_Q(GEN x, long prec)
{
  GEN y;
  long s;
  if (!(y = initsmall(x, 8))) return NULL;
  s = gsigne( ell_get_disc(y) );
  gel(y,14) = mkvecsmall(t_ELL_Q);
  gel(y,15) = mkvec(mkvecsmall2(prec2nbits(prec), s));
  return y;
}

static GEN
nfVtoalg(GEN nf, GEN x) { pari_APPLY_same(nftoalg(nf,gel(x,i))); }

static GEN
ellinit_nf(GEN x, GEN p)
{
  GEN y, nf;
  if (lg(x) > 6) x = vecslice(x,1,5);
  nf = checknf(p);
  x = nfVtoalg(nf, x);
  if (!(y = initsmall(x, 5))) return NULL;
  gel(y,14) = mkvecsmall(t_ELL_NF);
  gel(y,15) = mkvec(p);
  return y;
}

/* FF_ellinit allows singular cubic, return NULL in that case */
static GEN
FF_ellinit_ns(GEN x, GEN fg)
{
  x = FF_ellinit(x,fg);
  return FF_equal0(ell_get_disc(x))? NULL: x;
}

static GEN
to_mod(GEN x, GEN p) { return mkintmod(Rg_to_Fp(x,p), p); }
static GEN
ellinit_Fp(GEN x, GEN p)
{
  long i;
  GEN y, disc;
  if (lg(x) > 6) switch(ell_get_type(x))
  {
    case t_ELL_Q: break;
    case t_ELL_Fp: chk_p(ellff_get_p(x),p); break;
    case t_ELL_Qp: chk_p(ellQp_get_p(x),p); break;
    default: pari_err_TYPE("elliptic curve base_ring", x);
  }
  if (lg(x) == 2) x = ellfromj(to_mod(gel(x,1), p));
  if (!(y = initsmall(x, 4))) return NULL;
  /* ell_to_a4a6_bc does not handle p<=3 */
  if (abscmpiu(p, 3) <= 0) return FF_ellinit_ns(y,p_to_FF(p,0));
  disc = Rg_to_Fp(ell_get_disc(y),p);
  if (!signe(disc)) return NULL;
  for(i = 1; i <= 13; i++) gel(y,i) = to_mod(gel(y,i),p);
  gel(y,14) = mkvecsmall(t_ELL_Fp);
  gel(y,15) = mkvec2(p, ell_to_a4a6_bc(y, p));
  return y;
}

static GEN
ellinit_Fq(GEN x, GEN fg)
{
  GEN y;
  if (lg(x) == 2)
  {
    GEN j = gel(x,1);
    if (typ(j) != t_FFELT) j = Fq_to_FF(j, fg);
    x = ellfromj(j);
  }
  if (!(y = initsmall(x, 4))) return NULL;
  return FF_ellinit_ns(y,fg);
}

static GEN
ellnf_to_Fq(GEN nf, GEN x, GEN P, GEN *pp, GEN *pT)
{
  GEN e = vecslice(x,1,5);
  GEN p, modP;
  if (get_modpr(P))
  { /* modpr accept */
    modP = P;
    p = modpr_get_p(modP);
  }
  else
  { /* pr, initialize modpr */
    GEN d = Q_denom(e);
    p = pr_get_p(P);
    modP = dvdii(d,p)? nfmodprinit(nf,P): zkmodprinit(nf,P);
  }
  *pp = p;
  *pT = modpr_get_T(modP);
  return nfV_to_FqV(e, nf, modP);
}
static GEN
ellinit_nf_to_Fq(GEN nf, GEN E, GEN P)
{
  GEN T,p;
  E = ellnf_to_Fq(nf, E, P, &p, &T);
  return T? ellinit_Fq(E,Tp_to_FF(T,p)): ellinit_Fp(E,p);
}

static GEN
ellinit_i(GEN x, GEN D, long prec)
{
  GEN y;

  switch(typ(x))
  {
    case t_STR: x = gel(ellsearchcurve(x),2); break;
    case t_VEC:
      switch(lg(x))
      {
        case 2: case 3: case 6: case 17: break;
        default: pari_err_TYPE("ellinit",x);
      }
      break;
    default: pari_err_TYPE("ellinit",x);
  }
  if (D && get_prid(D))
  {
    if (lg(x) == 6 || ell_get_type(x) != t_ELL_NF) pari_err_TYPE("ellinit",x);
    return ellinit_nf_to_Fq(ellnf_get_nf(x), x, D);
  }
  switch (base_ring(x, &D, &prec))
  {
  case t_PADIC:
    y = ellinit_Qp(x, D, prec);
    break;
  case t_INTMOD:
    y = ellinit_Fp(x, D);
    break;
  case t_FFELT:
    y = ellinit_Fq(x, D);
    break;
  case t_FRAC:
    y = ellinit_Q(x, prec);
    break;
  case t_REAL:
    y = ellinit_Rg(x, 2, prec);
    break;
  case t_VEC:
    y = ellinit_nf(x, D);
    break;
  default:
    y = ellinit_Rg(x, 0, prec);
  }
  return y;
}
GEN
ellinit(GEN x, GEN D, long prec)
{
  pari_sp av = avma;
  GEN y = ellinit_i(x, D, prec);
  if (!y) { set_avma(av); return cgetg(1,t_VEC); }
  return gerepilecopy(av,y);
}

/********************************************************************/
/**                                                                **/
/**                     COORDINATE CHANGE                          **/
/**  Apply [u,r,s,t]. All coordch_* functions update E[1..14] only **/
/**  and copy E[15] (type-specific data), E[16] (dynamic data)     **/
/**  verbatim                                                      **/
/**                                                                **/
/********************************************************************/
/* [1,0,0,0] */
static GEN
init_ch(void) { return mkvec4(gen_1,gen_0,gen_0,gen_0); }
static int
is_trivial_change(GEN v)
{
  GEN u, r, s, t;
  if (typ(v) == t_INT) return 1;
  u = gel(v,1); r = gel(v,2); s = gel(v,3); t = gel(v,4);
  return isint1(u) && isintzero(r) && isintzero(s) && isintzero(t);
}

/* Accumulate the effects of variable changes w o v, where
 * w = [u,r,s,t], *vtotal = v = [U,R,S,T]. No assumption on types */
static void
gcomposev(GEN *vtotal, GEN w)
{
  GEN v = *vtotal;
  GEN U2, U, R, S, T, u, r, s, t;

  if (!v || typ(v) == t_INT) { *vtotal = w; return; }
  U = gel(v,1); R = gel(v,2); S = gel(v,3); T = gel(v,4);
  u = gel(w,1); r = gel(w,2); s = gel(w,3); t = gel(w,4);
  U2 = NULL;
  if (!gequal0(r))
  {
    GEN rU2;
    U2 = gsqr(U); rU2 = gmul(U2, r);
    R = gadd(R, rU2);
    T = gadd(T, gmul(S, rU2));
  }
  if (!gequal0(s)) S = gadd(S, gmul(U, s));
  if (!gequal0(t))
  {
    if (!U2) U2 = gsqr(U);
    T = gadd(T, gmul(gmul(U,U2), t));
  }
  gel(v,1) = gmul(U, u);
  gel(v,2) = R;
  gel(v,3) = S;
  gel(v,4) = T;
}

/* [u,r,s,t]^-1 = [ 1/u,-r/u^2,-s/u, (rs-t)/u^3 ] */
GEN
ellchangeinvert(GEN w)
{
  GEN u,r,s,t, u2,u3, U,R,S,T;
  if (typ(w) == t_INT) return w;
  u = gel(w,1);
  r = gel(w,2);
  s = gel(w,3);
  t = gel(w,4);
  u2 = gsqr(u); u3 = gmul(u2,u);
  U = ginv(u);
  R = gdiv(gneg(r), u2);
  S = gdiv(gneg(s), u);
  T = gdiv(gsub(gmul(r,s), t), u3);
  return mkvec4(U,R,S,T);
}

static GEN
ell_to_nfell10(GEN e)
{
  long i;
  GEN nf = ellnf_get_nf(e);
  GEN y = cgetg(11,t_VEC);
  for(i=1; i<=10; i++)
    gel(y, i) = nf_to_scalar_or_basis(nf, gel(e, i));
  return y;
}

/* apply [u^(-1),0,0,0] */
static GEN
nf_coordch_uinv(GEN nf, GEN e, GEN u)
{
  GEN y, u2, u3, u4, u6, u8;
  long lx;
  if (gequal1(u)) return e;
  y = cgetg_copy(e, &lx);
  u2 = nfsqr(nf,u); u3 = nfmul(nf,u,u2); u4 = nfsqr(nf,u2);
  u6 = nfsqr(nf,u3); u8 = nfsqr(nf,u4);
  gel(y,1) = nfmul(nf,ell_get_a1(e),  u);
  gel(y,2) = nfmul(nf,ell_get_a2(e), u2);
  gel(y,3) = nfmul(nf,ell_get_a3(e), u3);
  gel(y,4) = nfmul(nf,ell_get_a4(e), u4);
  gel(y,5) = nfmul(nf,ell_get_a6(e), u6);
  if (lx == 6) return y;
  gel(y,6) = nfmul(nf,ell_get_b2(e), u2);
  gel(y,7) = nfmul(nf,ell_get_b4(e), u4);
  gel(y,8) = nfmul(nf,ell_get_b6(e), u6);
  gel(y,9) = nfmul(nf,ell_get_b8(e), u8);
  return y;
}
/* apply [1,r,0,0] */
static GEN
nf_coordch_r(GEN nf, GEN e, GEN r)
{
  GEN a2, a4, b4, b6, y, p1, r2, b2r, rx3;
  long lx;
  if (gequal0(r)) return e;
  y = cgetg_copy(e, &lx);
  a2 = ell_get_a2(e); a4 = ell_get_a4(e);
  rx3 = gmulsg(3,r);

  gel(y,1) = ell_get_a1(e);
  /* A2 = a2 + 3r */
  gel(y,2) = nfadd(nf,a2,rx3);
  /* A3 = a1 r + a3 */
  gel(y,3) = nfadd(nf,ell_get_a3(e), nfmul(nf,ell_get_a1(e),r));
  /* A4 = 3r^2 + 2a2 r + a4 */
  gel(y,4) = nfadd(nf,a4, nfmul(nf,r,nfadd(nf,gmul2n(a2,1),rx3)));
  /* A6 = r^3 + a2 r^2 + a4 r + a6 */
  gel(y,5) = nfadd(nf,ell_get_a6(e),nfmul(nf,r,nfadd(nf, a4, nfmul(nf,r,nfadd(nf,a2, r)))));
  if (lx == 6) return y;

  b4 = ell_get_b4(e);
  b6 = ell_get_b6(e);
  /* B2 = 12r + b2 */
  gel(y,6) = nfadd(nf,ell_get_b2(e),gmul2n(rx3,2));
  b2r = nfmul(nf,r, ell_get_b2(e));
  r2 = nfsqr(nf,r);
  /* B4 = 6r^2 + b2 r + b4 */
  gel(y,7) = nfadd(nf,b4,nfadd(nf,b2r, gmulsg(6,r2)));
  /* B6 = 4r^3 + 2b2 r^2 + 2b4 r + b6 */
  gel(y,8) = nfadd(nf,b6,nfmul(nf,r,nfadd(nf,gmul2n(b4,1), nfadd(nf,b2r,gmul2n(r2,2)))));
  /* B8 = 3r^4 + b2 r^3 + 3b4 r^2 + 3b6 r + b8 */
  p1 = nfadd(nf,gmulsg(3,b4),nfadd(nf,b2r, gmulsg(3,r2)));
  gel(y,9) = nfadd(nf,ell_get_b8(e), nfmul(nf,r,nfadd(nf,gmulsg(3,b6), nfmul(nf,r,p1))));
  return y;
}

static GEN
nf_coordch_s(GEN nf, GEN e, GEN s)
{
  GEN a1, y;
  if (gequal0(s)) return e;
  a1 = ell_get_a1(e);
  y = leafcopy(e);

  /* A1 = a1 + 2s */
  gel(y,1) = nfadd(nf,a1,gmul2n(s,1));
  /* A2 = a2 - (a1 s + s^2) */
  gel(y,2) = nfsub(nf,ell_get_a2(e),nfmul(nf,s,nfadd(nf,a1,s)));
  /* A4 = a4 - s a3 */
  gel(y,4) = nfsub(nf,ell_get_a4(e),nfmul(nf,s,ell_get_a3(e)));
  return y;
}
/* apply [1,0,0,t] */
static GEN
nf_coordch_t(GEN nf, GEN e, GEN t)
{
  GEN a1, a3, y;
  if (gequal0(t)) return e;
  a1 = ell_get_a1(e); a3 = ell_get_a3(e);
  y = leafcopy(e);
  /* A3 = 2t + a3 */
  gel(y,3) = nfadd(nf,a3, gmul2n(t,1));
  /* A4 = a4 - a1 t */
  gel(y,4) = nfsub(nf,ell_get_a4(e), nfmul(nf,t,a1));
  /* A6 = a6 - t(t + a3) */
  gel(y,5) = nfsub(nf,ell_get_a6(e), nfmul(nf,t,nfadd(nf,t, a3)));
  return y;
}

/* apply [1,0,s,t] */
static GEN
nf_coordch_st(GEN nf, GEN e, GEN s, GEN t)
{
  GEN y, a1, a3;
  if (gequal0(s)) return nf_coordch_t(nf, e, t);
  if (gequal0(t)) return nf_coordch_s(nf, e, s);
  a1 = ell_get_a1(e); a3 = ell_get_a3(e);
  y = leafcopy(e);
  /* A1 = a1 + 2s */
  gel(y,1) = nfadd(nf,a1,gmul2n(s,1));
  /* A2 = a2 - (a1 s + s^2) */
  gel(y,2) = nfsub(nf,ell_get_a2(e),nfmul(nf,s,nfadd(nf,a1,s)));
  /* A3 = 2t + a3 */
  gel(y,3) = nfadd(nf,a3,gmul2n(t,1));
  /* A4 = a4 - (a1 t + s (2t + a3)) */
  gel(y,4) = nfsub(nf,ell_get_a4(e),nfadd(nf,nfmul(nf,t,a1),nfmul(nf,s,gel(y,3))));
  /* A6 = a6 - t(t + a3) */
  gel(y,5) = nfsub(nf,ell_get_a6(e), nfmul(nf,t,nfadd(nf,t, a3)));
  return y;
}

static GEN
nf_coordch_rt(GEN nf, GEN e, GEN r, GEN t)
{
  e = nf_coordch_r(nf, e, r);
  return nf_coordch_t(nf, e, t);
}

/* apply [1,r,s,t] */
static GEN
nf_coordch_rst(GEN nf, GEN e, GEN r, GEN s, GEN t)
{
  e = nf_coordch_r(nf, e, r);
  return nf_coordch_st(nf, e, s, t);
}
/* apply w = [u,r,s,t] */
static GEN
nf_coordch(GEN nf, GEN e, GEN w)
{
  if (typ(w) == t_INT) return e;
  e = nf_coordch_rst(nf, e, gel(w,2), gel(w,3), gel(w,4));
  return nf_coordch_uinv(nf, e, nfinv(nf, gel(w,1)));
}

/* apply [u^(-1),0,0,0] */
static GEN
coordch_uinv(GEN e, GEN u)
{
  GEN y, u2, u3, u4, u6, u12, D, c4, c6;
  long lx;
  if (gequal1(u)) return e;
  y = cgetg_copy(e, &lx);
  u2 = gsqr(u); u3 = gmul(u,u2); u4 = gsqr(u2); u6 = gsqr(u3);
  gel(y,1) = gmul(ell_get_a1(e),  u);
  gel(y,2) = gmul(ell_get_a2(e), u2);
  gel(y,3) = gmul(ell_get_a3(e), u3);
  gel(y,4) = gmul(ell_get_a4(e), u4);
  gel(y,5) = gmul(ell_get_a6(e), u6);
  if (lx == 6) return y;
  gel(y,6) = gmul(ell_get_b2(e), u2);
  gel(y,7) = gmul(ell_get_b4(e), u4);
  gel(y,8) = gmul(ell_get_b6(e), u6);
  gel(y,9) = gmul(ell_get_b8(e), gsqr(u4));
  u12 = gsqr(u6);
  D = ell_get_disc(e);
  c4 = ell_get_c4(e);
  c6 = ell_get_c6(e);
  c4 = gmul(c4, u4);
  c6 = gmul(c6, u6);
  D = gmul(D, u12);
  gel(y,10)= c4;
  gel(y,11)= c6;
  gel(y,12)= D;
  gel(y,13)= ell_get_j(e);
  gel(y,14)= gel(e,14);
  gel(y,15)= gel(e,15);
  gel(y,16)= gel(e,16);
  return y;
}
/* apply [1,r,0,0] */
static GEN
coordch_r(GEN e, GEN r)
{
  GEN a2, b4, b6, y, p1, r2, b2r, rx3;
  if (gequal0(r)) return e;
  y = leafcopy(e);
  a2 = ell_get_a2(e);
  rx3 = gmulsg(3,r);

  /* A2 = a2 + 3r */
  gel(y,2) = gadd(a2,rx3);
  /* A3 = a1 r + a3 */
  gel(y,3) = ec_h_evalx(e,r);
  /* A4 = 3r^2 + 2a2 r + a4 */
  gel(y,4) = gadd(ell_get_a4(e), gmul(r,gadd(gmul2n(a2,1),rx3)));
  /* A6 = r^3 + a2 r^2 + a4 r + a6 */
  gel(y,5) = ec_f_evalx(e,r);
  if (lg(y) == 6) return y;

  b4 = ell_get_b4(e);
  b6 = ell_get_b6(e);
  /* B2 = 12r + b2 */
  gel(y,6) = gadd(ell_get_b2(e),gmul2n(rx3,2));
  b2r = gmul(r, ell_get_b2(e));
  r2 = gsqr(r);
  /* B4 = 6r^2 + b2 r + b4 */
  gel(y,7) = gadd(b4,gadd(b2r, gmulsg(6,r2)));
  /* B6 = 4r^3 + 2b2 r^2 + 2b4 r + b6 */
  gel(y,8) = gadd(b6,gmul(r,gadd(gmul2n(b4,1), gadd(b2r,gmul2n(r2,2)))));
  /* B8 = 3r^4 + b2 r^3 + 3b4 r^2 + 3b6 r + b8 */
  p1 = gadd(gmulsg(3,b4),gadd(b2r, gmulsg(3,r2)));
  gel(y,9) = gadd(ell_get_b8(e), gmul(r,gadd(gmulsg(3,b6), gmul(r,p1))));
  return y;
}
/* apply [1,0,s,0] */
static GEN
coordch_s(GEN e, GEN s)
{
  GEN a1, y;
  if (gequal0(s)) return e;
  a1 = ell_get_a1(e);
  y = leafcopy(e);

  /* A1 = a1 + 2s */
  gel(y,1) = gadd(a1,gmul2n(s,1));
  /* A2 = a2 - (a1 s + s^2) */
  gel(y,2) = gsub(ell_get_a2(e),gmul(s,gadd(a1,s)));
  /* A4 = a4 - s a3 */
  gel(y,4) = gsub(ell_get_a4(e),gmul(s,ell_get_a3(e)));
  return y;
}
/* apply [1,0,0,t] */
static GEN
coordch_t(GEN e, GEN t)
{
  GEN a1, a3, y;
  if (gequal0(t)) return e;
  a1 = ell_get_a1(e); a3 = ell_get_a3(e);
  y = leafcopy(e);
  /* A3 = 2t + a3 */
  gel(y,3) = gadd(a3, gmul2n(t,1));
  /* A4 = a4 - a1 t */
  gel(y,4) = gsub(ell_get_a4(e), gmul(t,a1));
  /* A6 = a6 - t(t + a3) */
  gel(y,5) = gsub(ell_get_a6(e), gmul(t,gadd(t, a3)));
  return y;
}
/* apply [1,0,s,t] */
static GEN
coordch_st(GEN e, GEN s, GEN t)
{
  GEN y, a1, a3;
  if (gequal0(s)) return coordch_t(e, t);
  if (gequal0(t)) return coordch_s(e, s);
  a1 = ell_get_a1(e); a3 = ell_get_a3(e);
  y = leafcopy(e);
  /* A1 = a1 + 2s */
  gel(y,1) = gadd(a1,gmul2n(s,1));
  /* A2 = a2 - (a1 s + s^2) */
  gel(y,2) = gsub(ell_get_a2(e),gmul(s,gadd(a1,s)));
  /* A3 = 2t + a3 */
  gel(y,3) = gadd(a3,gmul2n(t,1));
  /* A4 = a4 - (a1 t + s (2t + a3)) */
  gel(y,4) = gsub(ell_get_a4(e),gadd(gmul(t,a1),gmul(s,gel(y,3))));
  /* A6 = a6 - t(t + a3) */
  gel(y,5) = gsub(ell_get_a6(e), gmul(t,gadd(t, a3)));
  return y;
}
/* apply [1,r,s,t] */
static GEN
coordch_rst(GEN e, GEN r, GEN s, GEN t)
{
  e = coordch_r(e, r);
  return coordch_st(e, s, t);
}
/* apply w = [u,r,s,t] */
static GEN
coordch(GEN e, GEN w)
{
  if (typ(w) == t_INT) return e;
  e = coordch_rst(e, gel(w,2), gel(w,3), gel(w,4));
  return coordch_uinv(e, ginv(gel(w,1)));
}

/* the ch_* routines update E[14] (type), E[15] (type specific data), E[16]
 * (dynamic data) */
static GEN
ch_Qp(GEN E, GEN e, GEN w)
{
  GEN S, p = ellQp_get_zero(E), u2 = NULL, u = gel(w,1), r = gel(w,2);
  long prec = valp(p);
  if (base_ring(E, &p, &prec) != t_PADIC) return ellinit(E, p, prec);
  if ((S = obj_check(e, Qp_ROOT)))
  {
    if (!u2) u2 = gsqr(u);
    obj_insert_shallow(E, Qp_ROOT, gdiv(gsub(S, r), u2));
  }
  if ((S = obj_check(e, Qp_TATE)))
  {
    GEN U2 = gel(S,1), U = gel(S,2), Q = gel(S,3), AB = gel(S,4), L = gel(S,5);
    if (!u2) u2 = gsqr(u);
    U2 = gmul(U2, u2);
    U = gmul(U, u);
    AB = gdiv(AB, u2);
    obj_insert_shallow(E, Qp_TATE, mkvec5(U2,U,Q,AB,L));
  }
  return E;
}

/* common to Q and Rg */
static GEN
ch_R(GEN E, GEN e, GEN w)
{
  GEN S, u = gel(w,1), r = gel(w,2);
  if ((S = obj_check(e, R_PERIODS)))
    obj_insert_shallow(E, R_PERIODS, gmul(S, u));
  if ((S = obj_check(e, R_ETA)))
    obj_insert_shallow(E, R_ETA, gmul(S, u));
  if ((S = obj_check(e, R_ROOTS)))
  {
    GEN ro = cgetg(4, t_VEC), u2 = gsqr(u);
    long i;
    for (i = 1; i <= 3; i++) gel(ro,i) = gdiv(gsub(gel(S,i), r), u2);
    obj_insert_shallow(E, R_ROOTS, ro);
  }
  return E;
}

static GEN
ch_Rg(GEN E, GEN e, GEN w)
{
  GEN p = NULL;
  long prec = ellR_get_prec(E);
  if (base_ring(E, &p, &prec) != t_REAL) return ellinit(E, p, prec);
  ch_R(E, e, w); return E;
}

static GEN
ch_Q(GEN E, GEN e, GEN w)
{
  long prec = ellR_get_prec(E);
  GEN S, v = NULL, p = NULL;
  if (base_ring(E, &p, &prec) != t_FRAC) return ellinit(E, p, prec);
  ch_R(E, e, w);
  if ((S = obj_check(e, Q_GROUPGEN)))
    S = obj_insert_shallow(E, Q_GROUPGEN, ellchangepoint(S, w));
  if ((S = obj_check(e, Q_MINIMALMODEL)))
  {
    if (lg(S) == 2)
    { /* model was minimal */
      if (!is_trivial_change(w)) /* no longer minimal */
        S = mkvec3(gel(S,1), ellchangeinvert(w), e);
      (void)obj_insert_shallow(E, Q_MINIMALMODEL, S);
    }
    else
    {
      v = gel(S,2);
      if (gequal(v, w) || (is_trivial_change(v) && is_trivial_change(w)))
        S = mkvec(gel(S,1)); /* now minimal */
      else
      {
        w = ellchangeinvert(w);
        gcomposev(&w, v); v = w;
        S = leafcopy(S); /* don't modify S in place: would corrupt e */
        gel(S,2) = v;
      }
      (void)obj_insert_shallow(E, Q_MINIMALMODEL, S);
    }
  }
  if ((S = obj_check(e, Q_GLOBALRED)))
    S = obj_insert_shallow(E, Q_GLOBALRED, S);
  if ((S = obj_check(e, Q_ROOTNO)))
    S = obj_insert_shallow(E, Q_ROOTNO, S);
  return E;
}

static void
ch_FF(GEN E, GEN e, GEN w)
{
  GEN S;
  if ((S = obj_check(e, FF_CARD)))
    S = obj_insert_shallow(E, FF_CARD, S);
  if ((S = obj_check(e, FF_GROUP)))
    S = obj_insert_shallow(E, FF_GROUP, S);
  if ((S = obj_check(e, FF_GROUPGEN)))
    S = obj_insert_shallow(E, FF_GROUPGEN, ellchangepoint(S, w));
  if ((S = obj_check(e, FF_O)))
    S = obj_insert_shallow(E, FF_O, S);
}

/* FF_CARD, FF_GROUP, FF_O are invariant */
static GEN
ch_Fp(GEN E, GEN e, GEN w)
{
  long prec = 0;
  GEN p = ellff_get_field(E);
  if (base_ring(E, &p, &prec) != t_INTMOD) return ellinit(E, p, prec);
  gel(E,15) = mkvec2(p, ell_to_a4a6_bc(E, p));
  ch_FF(E, e, w); return E;
}
static GEN
ch_Fq(GEN E, GEN e, GEN w)
{
  long prec = 0;
  GEN p = ellff_get_field(E);
  if (base_ring(E, &p, &prec) != t_FFELT) return ellinit(E, p, prec);
  gel(E,15) = FF_elldata(E, p);
  ch_FF(E, e, w); return E;
}

static void
ell_reset(GEN E)
{ gel(E,16) = zerovec(lg(gel(E,16))-1); }

GEN
ellchangecurve(GEN e, GEN w)
{
  pari_sp av = avma;
  GEN E;
  checkell5(e);
  if (equali1(w)) return gcopy(e);
  checkcoordch(w);
  E = coordch(leafcopy(e), w);
  if (lg(E) != 6)
  {
    ell_reset(E);
    switch(ell_get_type(E))
    {
      case t_ELL_Qp: E = ch_Qp(E,e,w); break;
      case t_ELL_Fp: E = ch_Fp(E,e,w); break;
      case t_ELL_Fq: E = ch_Fq(E,e,w); break;
      case t_ELL_Q:  E = ch_Q(E,e,w);  break;
      case t_ELL_Rg: E = ch_Rg(E,e,w); break;
    }
  }
  return gerepilecopy(av, E);
}

/* v o= [1,r,0,0] */
static void
nf_compose_r(GEN nf, GEN *vtotal, GEN *e, GEN r)
{
  GEN v = *vtotal;
  GEN U2, R, S, T;
  if (gequal0(r)) return;
  *e = nf_coordch_r(nf, *e,r);
  U2 = nfsqr(nf,gel(v,1)); R = gel(v,2); S = gel(v, 3); T = gel(v, 4);
  gel(v,2) = nfadd(nf,R, nfmul(nf,U2, r));
  gel(v,4) = nfadd(nf,T, nfmul(nf,U2, nfmul(nf,S, r)));
}
/* v o= [1,0,s,0]; never used for s = 0 */
static void
nf_compose_s(GEN nf, GEN *vtotal, GEN *e, GEN s)
{
  GEN v = *vtotal;
  GEN U, S;
  *e = nf_coordch_s(nf,*e,s);
  U = gel(v,1); S = gel(v,3);
  gel(v,3) = nfadd(nf, S, nfmul(nf, U, s));
}
/* v o= [1,0,0,t] */
static void
nf_compose_t(GEN nf ,GEN *vtotal, GEN *e, GEN t)
{
  GEN v = *vtotal;
  GEN U3, U, T;
  if (gequal0(t)) return;
  *e = nf_coordch_t(nf,*e,t);
  U = gel(v,1); U3 = nfmul(nf,U, nfsqr(nf,U)); T = gel(v,4);
  gel(v,4) = nfadd(nf,T, nfmul(nf,U3, t));
}
/* v o= [1,r,0,t] */
static void
nf_compose_rt(GEN nf, GEN *vtotal, GEN *e, GEN r, GEN t)
{
  GEN v = *vtotal;
  GEN U2, U, R, S, T;
  if (gequal0(t)) { nf_compose_r(nf, vtotal, e, r); return; }
  *e = nf_coordch_rt(nf,*e,r,t);
  U = gel(v,1); R = gel(v,2); S = gel(v,3); T = gel(v,4);
  U2 = nfsqr(nf,U);
  gel(v,2) = nfadd(nf,R, nfmul(nf,U2, r));
  gel(v,4) = nfadd(nf,T, nfmul(nf,U2, nfadd(nf,nfmul(nf,U, t), nfmul(nf,S, r))));
}
/* v o= [1,0,s,t] */
static void
nf_compose_st(GEN nf, GEN *vtotal, GEN *e, GEN s, GEN t)
{
  GEN v = *vtotal;
  GEN U3, U, S, T;
  if (gequal0(s)) { nf_compose_t(nf, vtotal, e, t); return; }
  if (gequal0(t)) { nf_compose_s(nf, vtotal, e, s); return; }
  *e = nf_coordch_st(nf, *e,s,t);
  U = gel(v,1); U3 = nfmul(nf,U,nfsqr(nf,U)); S = gel(v,3); T = gel(v,4);
  gel(v,3) = nfadd(nf, S, nfmul(nf,U, s));
  gel(v,4) = nfadd(nf, T, nfmul(nf,U3, t));
}

/* v o= [u,0,0,0] */
static void
nf_compose_u(GEN nf, GEN *vtotal, GEN *e, GEN u, GEN uinv)
{
  GEN v = *vtotal;
  *e = nf_coordch_uinv(nf, *e,uinv); gel(v,1) = nfmul(nf,gel(v,1), u);
}

/* X = (x-r)/u^2
 * Y = (y - s(x-r) - t) / u^3 */
static GEN
ellchangepoint0(GEN P, GEN v2, GEN v3, GEN r, GEN s, GEN t)
{
  GEN a, x, y;
  long l = lg(P);
  if (typ(P) != t_VEC || l == 1 || l > 4) pari_err_TYPE("ellchangepoint", P);
  if (ell_is_inf(P)) return P;
  x = gel(P,1); y = gel(P,2); a = gsub(x,r);
  retmkvec2(gmul(v2, a), gmul(v3, gsub(y, gadd(gmul(s,a),t))));
}

GEN
ellchangepoint(GEN x, GEN ch)
{
  GEN y, v, v2, v3, r, s, t, u;
  long tx, i, lx = lg(x);
  pari_sp av = avma;

  if (typ(x) != t_VEC) pari_err_TYPE("ellchangepoint",x);
  if (equali1(ch)) return gcopy(x);
  checkcoordch(ch);
  if (lx == 1) return cgetg(1, t_VEC);
  u = gel(ch,1); r = gel(ch,2); s = gel(ch,3); t = gel(ch,4);
  v = ginv(u); v2 = gsqr(v); v3 = gmul(v,v2);
  tx = typ(gel(x,1));
  if (is_matvec_t(tx))
  {
    y = cgetg(lx,tx);
    for (i=1; i<lx; i++)
      gel(y,i) = ellchangepoint0(gel(x,i),v2,v3,r,s,t);
  }
  else
    y = ellchangepoint0(x,v2,v3,r,s,t);
  return gerepilecopy(av,y);
}

/* x = u^2*X + r
 * y = u^3*Y + s*u^2*X + t */
static GEN
ellchangepointinv0(GEN P, GEN u2, GEN u3, GEN r, GEN s, GEN t)
{
  GEN a, X, Y;
  if (ell_is_inf(P)) return P;
  X = gel(P,1); Y = gel(P,2); a = gmul(u2,X);
  return mkvec2(gadd(a, r), gadd(gmul(u3, Y), gadd(gmul(s, a), t)));
}
GEN
ellchangepointinv(GEN x, GEN ch)
{
  GEN y, u, r, s, t, u2, u3;
  long tx, i, lx = lg(x);
  pari_sp av = avma;

  if (typ(x) != t_VEC) pari_err_TYPE("ellchangepointinv",x);
  if (equali1(ch)) return gcopy(x);
  checkcoordch(ch);
  if (lx == 1) return cgetg(1, t_VEC);
  u = gel(ch,1); r = gel(ch,2); s = gel(ch,3); t = gel(ch,4);
  u2 = gsqr(u); u3 = gmul(u,u2);
  tx = typ(gel(x,1));
  if (is_matvec_t(tx))
  {
    y = cgetg(lx,tx);
    for (i=1; i<lx; i++)
      gel(y,i) = ellchangepointinv0(gel(x,i),u2,u3,r,s,t);
  }
  else
    y = ellchangepointinv0(x,u2,u3,r,s,t);
  return gerepilecopy(av,y);
}

static GEN
elltwist_card(GEN h, GEN q) { return subii(shifti(addiu(q, 1), 1), h); }
GEN
elltwist(GEN E, GEN P)
{
  pari_sp av = avma;
  GEN a1, a2, a3, a4, a6, a, b, c, ac, D, D2, V, DOM = NULL;
  long prec = DEFAULTPREC, isell = (lg(E) == 17);

  if (typ(E) != t_VEC) pari_err_TYPE("elltwist",E);
  if (isell) switch(ell_get_type(E))
  {
    case t_ELL_Q:
    case t_ELL_Rg: prec = ellR_get_prec(E); break;
  }
  if (!P)
  {
    GEN Et, S, a4, a6, e, fg, q;
    if (!isell)
    { /* Could avoid this ellinit. Don't bother. */
      e = E; E = ellinit_i(E, NULL, prec);
      if (!E) pari_err_TYPE("elltwist", e);
    }
    switch (ell_get_type(E))
    {
      case t_ELL_Fp:
        q = ellff_get_field(E);
        e = ellff_get_a4a6(E);
        Fp_elltwist(gel(e,1), gel(e,2), q, &a4, &a6);
        Et = ellinit_Fp(mkvec2(a4,a6), q); break;
      case t_ELL_Fq:
        fg = ellff_get_field(E); q = FF_q(fg);
        Et = ellinit_Fq(FF_elltwist(E), fg); break;
      default: pari_err_TYPE("elltwist [missing P]", E);
        return NULL;/*LCOV_EXCL_LINE*/
    }
    if ((S = obj_check(E, FF_CARD)))
      obj_insert_shallow(Et, FF_CARD, elltwist_card(S, q));
    return gerepilecopy(av, Et);
  }
  if (isell && ell_get_type(E) == t_ELL_NF)
    if (!(DOM = ellnf_get_bnf(E))) DOM = ellnf_get_nf(E);
  if (typ(P) == t_INT)
  {
    if (equali1(P)) return ellinit(E, DOM, prec);
    P = quadpoly(P);
  }
  else
  {
    if (typ(P) != t_POL) pari_err_TYPE("elltwist",P);
    if (degpol(P) != 2 )
      pari_err_DOMAIN("elltwist", "degree(P)", "!=", gen_2, P);
  }
  switch(lg(E))
  {
    case 3:
      a1 = a2 = a3 = gen_0;
      a4 = gel(E,1);
      a6 = gel(E,2); break;
    case 6: case 17:
      a1 = ell_get_a1(E);
      a2 = ell_get_a2(E);
      a3 = ell_get_a3(E);
      a4 = ell_get_a4(E);
      a6 = ell_get_a6(E); break;
    default:
      pari_err_TYPE("ellxxx [not an elliptic curve (ell5)]",E);
      return NULL; /* LCOV_EXCL_LINE */
  }
  a = gel(P,4); b = gel(P,3); c = gel(P,2); ac = gmul(a, c);
  D = gsub(gsqr(b), gmulsg(4, ac)); D2 = gsqr(D);
  if (gequal0(a1) && gequal0(a2) && gequal0(a3))
    V = mkvec2(gmul(a4, D2), gmul(gmul(a6, D), D2));
  else
  {
    GEN a3D = gmul(a3, D);
    V = cgetg(6, t_VEC);
    gel(V,1) = gmul(a1, b);
    gel(V,2) = gsub(gmul(a2, D), gmul(gsqr(a1), ac));
    gel(V,3) = gmul(a3D, b);
    gel(V,4) = gsub(gmul(a4, D2), gmul(gmul(gmulsg(2, a3D), a1), ac));
    gel(V,5) = gmul(gsub(gmul(a6, D), gmul(gsqr(a3), ac)), D2);
  }
  E = ellinit_i(V, DOM, prec);
  if (!E) pari_err_TYPE("elltwist", V);
  return gerepilecopy(av, E);
}

/********************************************************************/
/**                      E/Q: MINIMAL TWIST                        **/
/**      Cf Ian Connell, Elliptic Curve Handbook, chap. 5          **/
/**                http://www.math.mcgill.ca/connell/              **/
/********************************************************************/

static long
safe_Z_lval(GEN n, ulong p)
{ return signe(n)==0? -1: Z_lval(n, p); }

/* Twist by d2 = 1,-4,-8,8, to get minimal discriminant at 2 after
 * ellminimalmodel / ellQ_minimalu; assume vg = min(3*v4,2*v6,vD) >= 6.
 * If nontrivial, v(d2) = 2 or 3 and let t = [(vg+6v(d2))/12].
 * Good case if reduction in ellQ_minimalu i.e. t = 2 (v4 < 6 or v6 < 9 or
 * vD < 18) or 3 and "d--" does not occur. Minimal model => t = 3 iff
 * v4 = 6, v6 = 9 and vD >= 18. Total net effect is
 *   v4 += 2v(d2) - 4t, v6 += 3v(d2) - 6t, vD += 6 v(d2) - 12t
 * After rescaling in ellQ_minimalu (c4 >>= 4t, c6 >>= 6t) we need
 *   c6 % 4 = 3 OR  (v4 >= 4 AND (v6 >= 5 or c6 % 32 = 8)) */
static long
twist2(GEN c4, GEN c6, GEN disc, long vg)
{ /* v4 >= 4, v6 >= 3 (and c6 = 0,8 mod 32). After twist + minimization,
   * either same condition OR v(C4) = 0, C6 = 0,3 mod 4 */
  long v4, v6, vD;

  if (vg == 18) /* v4=6, v6=9, vD>=18; only case with t = 3 */
    return (umodi2n(c6, 11)>>9) == 1 ? -8: 8; /* need C6 % 4 = 3 */

  /* 100 = oo, any number >= 8 would do */
  v4 = signe(c4)? vali(c4): 100; if (v4 == 5) return 1;
  /* 100 = oo, any number > 9 would do */
  v6 = signe(c6)? vali(c6): 100; if (v6 == 7) return 1;

  /* handle case v(DISC) = 0 or v(C4) = 0 after twist, only case with d2 = -4 */
  if (vg == 12 && ((v4==4 && v6==6) || (v4>=8 && v6==9))) return -4;

  /* Now, d2 = 1 OR v(d2) = 3, t = 2 => v4 -= 2, v6 -= 3, vD -= 6 */
  if (v4 < 6 || v6 < 6) return 1; /* v(C4) >= 4, v(C6) >= 3 */
  vD = vali(disc);
  if (v6==6 && vD==6 && (umodi2n(c6,8)>>6) == 1) return 8; /* C6 % 32 = 8 */
  return -8;
}

/* Return D such that E_D has minimal discriminant. It also has minimal
 * conductor in Z[1/2] */
GEN
ellminimaltwist(GEN e)
{
  pari_sp av = avma;
  GEN c4, c6, disc, g, N, M, F, E, D = gen_1;
  long i, lF;
  checkell_Q(e);
  E = ellminimalmodel(e, NULL);
  c4 = ell_get_c4(E);
  c6 = ell_get_c6(E);
  disc = ell_get_disc(E);
  g = gcdii(disc, sqri(c6));
  ellQ_get_Nfa(E, &N, &M);
  F = gel(M, 1); lF = lg(F);
  /* on twist by d, (c4,c6,D,g) -> (d^2 c4, d^3 c6, d^6 D, d^6 g),
   * then apply ellQ_minimalu(). Since model is minimal, v(g) < 12 unless p=3
   * and v(g) < 14 or p = 2 and v(g) <= 18 */
  for(i = 1; i < lF; i++)
  {
    GEN p = gel(F, i);
    long vg = Z_pval(g,p), d2;
    if (vg < 6) continue;
    /* twist by fund. discriminant d2; in ellQ_minimalu,
     * we have v(g) = vg + 6*v(d2) */
    switch(itou_or_0(p))
    {
      default: /* p > 3, 6 <= v(g) < 12 => v(D) -= 6 */
        D = mulii(D, (mod4(p)==1)? p: negi(p));
        break;
      case 3: /* bad case: v(final_c6) = 2 => no reduction; else v(D) -= 6 */
        if (safe_Z_lval(c6,3) != 5) D = mulis(D, -3);
        break;
      case 2:
        d2 = twist2(c4,c6,disc,vg);
        if (d2 != 1) D = mulis(D,d2);
        break;
    }
  }
  obj_free(E);
  return gerepileuptoleaf(av, D);
}

/*
Reference:
William A. Stein and Mark Watkins
A Database of Elliptic Curves-First Report
ANTS 5
<http://modular.math.washington.edu/papers/stein-watkins/ants.pdf>
*/
static GEN localred_23(GEN e, long p);
GEN
ellminimaltwistcond(GEN e)
{
  pari_sp av = avma;
  GEN D = ellminimaltwist(e), eD = elltwist(e, D);
  GEN R = localred_23(ellintegralmodel_i(eD,NULL), 2);
  long f = itos(gel(R,1)), v = vali(D);
  if (f==4) D = negi(v==3 ? D: shifti(D, v==0? 2: -2));
  else if (f==6)
  {
    long s, t;
    if (v < 3) s = v==0? 3: 1;
    else
    {
      t = (v==3 && mod32(D) == 8)? 1: -1;
      s = signe(D)==t ? -3: -1;
    }
    D = shifti(D, s);
  }
  return gerepileuptoleaf(av, D);
}

GEN
ellminimaltwist0(GEN e, long flag)
{
  switch(flag)
  {
    case 0: return ellminimaltwist(e);
    case 1: return ellminimaltwistcond(e);
  }
  pari_err_FLAG("ellminimaltwist");
  return NULL; /* LCOV_EXCL_LINE */
}

static long
ellexpo(GEN E)
{
  long i, f, e = -(long)HIGHEXPOBIT;
  for (i=1; i<=5; i++)
  {
    f = gexpo(gel(E,i));
    if (f > e) e = f;
  }
  return e;
}


static int
oncurve_exact(GEN e, GEN z)
{
  pari_sp av = avma;
  GEN A = ec_LHS_evalQ(e,z), B = ec_f_evalx(e,gel(z,1));
  return gc_bool(av, gequal(A, B));
}
/* Exactness of lhs and rhs in the following depends in nonobvious ways
 * on the coeffs of the curve as well as on the components of the point z.
 * Thus if e is exact, with a1==0, and z has exact y coordinate only, the
 * lhs will be exact but the rhs won't. */
int
oncurve(GEN e, GEN z)
{
  GEN LHS, RHS, x;
  long pl, pr, ex, expx;
  pari_sp av;

  checkellpt(z); if (ell_is_inf(z)) return 1; /* oo */
  if (ell_get_type(e) == t_ELL_NF) z = nfVtoalg(ellnf_get_nf(e), z);
  av = avma;
  LHS = ec_LHS_evalQ(e,z);
  RHS = ec_f_evalx(e,gel(z,1)); x = gsub(LHS,RHS);
  if (gequal0(x)) return gc_bool(av,1);
  pl = precision(LHS);
  pr = precision(RHS);
  if (!pl && !pr) return gc_bool(av,0); /* both of LHS, RHS are exact */
  /* at least one of LHS,RHS is inexact */
  ex = pr? gexpo(RHS): gexpo(LHS); /* don't take exponent of exact 0 */
  if (!pr || (pl && pl < pr)) pr = pl; /* min among nonzero elts of {pl,pr} */
  expx = gexpo(x);
  pr = (expx < ex - prec2nbits(pr) + 15
     || expx < ellexpo(e) - prec2nbits(pr) + 5);
  return gc_bool(av,pr);
}

GEN
ellisoncurve(GEN e, GEN x)
{
  long i, tx = typ(x), lx;

  checkell(e);
  if (!is_vec_t(tx)) pari_err_TYPE("ellisoncurve [point]", x);
  lx = lg(x); if (lx==1) return cgetg(1,tx);
  tx = typ(gel(x,1));
  if (is_vec_t(tx))
  {
    GEN z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = ellisoncurve(e,gel(x,i));
    return z;
  }
  return oncurve(e, x)? gen_1: gen_0;
}

/* y1 = y2 or -LHS0-y2 */
static GEN
slope_samex(GEN e, GEN x, GEN y1, GEN y2)
{
  GEN dy,dx;
  if (y1 != y2)
  {
    int eq;
    if (precision(y1) || precision(y2))
      eq = (gexpo(gadd(ec_h_evalx(e,x),gadd(y1,y2))) >= gexpo(y1));
    else
      eq = gequal(y1,y2);
    if (!eq) return NULL;
  }
  dx = ec_dmFdy_evalQ(e,mkvec2(x,y1));
  if (gequal0(dx)) return NULL;
  dy = gadd(gsub(ell_get_a4(e),gmul(ell_get_a1(e),y1)),
            gmul(x,gadd(gmul2n(ell_get_a2(e),1),gmulsg(3,x))));
  return gdiv(dy,dx);
}

GEN
elladd(GEN e, GEN z1, GEN z2)
{
  GEN s, z, x, y, x1, x2, y1, y2;
  pari_sp av = avma;

  checkell(e); checkellpt(z1); checkellpt(z2);
  if (ell_is_inf(z1)) return gcopy(z2);
  if (ell_is_inf(z2)) return gcopy(z1);

  x1 = gel(z1,1); y1 = gel(z1,2);
  x2 = gel(z2,1); y2 = gel(z2,2);
  if (ell_get_type(e) == t_ELL_NF)
  {
    GEN nf = ellnf_get_nf(e);
    x1 = nftoalg(nf, x1);
    x2 = nftoalg(nf, x2);
    y1 = nftoalg(nf, y1);
    y2 = nftoalg(nf, y2);
  }
  if (cx_approx_equal(x1,x2))
  {
    s = slope_samex(e, x1, y1, y2);
    if (!s) { set_avma(av); return ellinf(); }
  }
  else
    s = gdiv(gsub(y2,y1), gsub(x2,x1));
  x = gsub(gmul(s,gadd(s,ell_get_a1(e))), gadd(gadd(x1,x2),ell_get_a2(e)));
  y = gadd(gadd(y1, ec_h_evalx(e,x)), gmul(s,gsub(x,x1)));
  z = cgetg(3,t_VEC);
  gel(z,1) = gcopy(x);
  gel(z,2) = gneg(y); return gerepileupto(av, z);
}

static GEN
ellneg_i(GEN e, GEN z)
{
  GEN t, x, y;
  if (ell_is_inf(z)) return z;
  x = gel(z,1);
  y = gel(z,2);
  if (ell_get_type(e) == t_ELL_NF)
  {
    GEN nf = ellnf_get_nf(e);
    x = nftoalg(nf,x);
    y = nftoalg(nf,y);
  }
  t = cgetg(3,t_VEC);
  gel(t,1) = x;
  gel(t,2) = gneg_i(gadd(y, ec_h_evalx(e,x)));
  return t;
}

GEN
ellneg(GEN e, GEN z)
{
  pari_sp av;
  GEN t, y;
  checkell(e); checkellpt(z);
  if (ell_is_inf(z)) return z;
  t = cgetg(3,t_VEC);
  gel(t,1) = gcopy(gel(z,1));
  av = avma;
  y = gneg(gadd(gel(z,2), ec_h_evalx(e,gel(z,1))));
  gel(t,2) = gerepileupto(av, y);
  return t;
}

GEN
ellsub(GEN e, GEN z1, GEN z2)
{
  pari_sp av = avma;
  checkell(e); checkellpt(z2);
  return gerepileupto(av, elladd(e, z1, ellneg_i(e,z2)));
}

/* E an ell, x a scalar */
static GEN
ellordinate_i(GEN E, GEN x, long prec)
{
  pari_sp av = avma;
  GEN a, b, D, d, y, p, nf = NULL;

  if (ell_get_type(E) == t_ELL_NF)
  {
    nf = ellnf_get_nf(E);
    x = nftoalg(nf,x);
  }
  a = ec_f_evalx(E,x);
  b = ec_h_evalx(E,x);
  D = gadd(gsqr(b), gmul2n(a,2));
  /* solve y*(y+b) = a */
  if (gequal0(D)) {
    if (ell_get_type(E) == t_ELL_Fq && absequaliu(ellff_get_p(E),2))
      retmkvec( FF_sqrt(a) );
    b = gneg_i(b); y = cgetg(2,t_VEC);
    gel(y,1) = gmul2n(b,-1);
    return gerepileupto(av,y);
  }
  /* D != 0 */
  switch(ell_get_type(E))
  {
    case t_ELL_Fp: /* imply p!=2 */
      p = ellff_get_p(E);
      D = gel(D,2);
      if (kronecker(D, p) < 0) { set_avma(av); return cgetg(1,t_VEC); }
      d = Fp_sqrt(D, p);
      break;
    case t_ELL_Fq:
      if (absequaliu(ellff_get_p(E),2))
      {
        GEN F = FFX_roots(mkpoln(3, gen_1, b, a), D);
        if (lg(F) == 1) { set_avma(av); return cgetg(1,t_VEC); }
        return gerepileupto(av, F);
      }
      if (!FF_issquareall(D,&d)) { set_avma(av); return cgetg(1,t_VEC); }
      break;
    case t_ELL_Q:
      if (typ(x) == t_COMPLEX) { d = gsqrt(D, prec); break; }
      if (!issquareall(D,&d)) { set_avma(av); return cgetg(1,t_VEC); }
      break;

    case t_ELL_NF:
    {
      GEN T = mkpoln(3, gen_1, gen_0, gneg(D));
      setvarn(T, fetch_var_higher());
      d = nfroots(nf, T);
      delete_var();
      if (lg(d) == 1) { set_avma(av); return cgetg(1, t_VEC); }
      d = gel(d,1);
      break;
    }

    case t_ELL_Qp:
      p = ellQp_get_p(E);
      D = cvtop(D, p, ellQp_get_prec(E));
      if (!issquare(D)) { set_avma(av); return cgetg(1,t_VEC); }
      d = Qp_sqrt(D);
      break;

    default:
      d = gsqrt(D,prec);
  }
  a = gsub(d,b); y = cgetg(3,t_VEC);
  gel(y,1) = gmul2n(a, -1);
  gel(y,2) = gsub(gel(y,1),d);
  return gerepileupto(av,y);
}

GEN
ellordinate(GEN e, GEN x, long prec)
{
  checkell(e);
  if (is_matvec_t(typ(x)))
  {
    long i, lx;
    GEN v = cgetg_copy(x, &lx);
    for (i=1; i<lx; i++) gel(v,i) = ellordinate(e,gel(x,i),prec);
    return v;
  }
  return ellordinate_i(e, x, prec);
}

GEN
ellrandom(GEN E)
{
  GEN fg;
  checkell_Fq(E);
  fg = ellff_get_field(E);
  if (typ(fg)==t_FFELT)
    return FF_ellrandom(E);
  else
  {
    pari_sp av = avma;
    GEN p = fg, e = ellff_get_a4a6(E);
    GEN P = random_FpE(gel(e,1),gel(e,2),p);
    P = FpE_to_mod(FpE_changepoint(P,gel(e,3),p),p);
    return gerepileupto(av, P);
  }
}

/* n t_QUAD or t_COMPLEX, P != [0] */
static GEN
ellmul_CM(GEN e, GEN P, GEN n)
{
  GEN p1p, q1p, x, y, p0, p1, q0, q1, z1, z2, grdx, b2ov12, N = gnorm(n);
  long ln, vn;

  if (typ(N) != t_INT)
    pari_err_TYPE("ellmul (non integral CM exponent)",N);
  ln = itos_or_0(shifti(addiu(N, 1UL), 3));
  if (!ln) pari_err_OVERFLOW("ellmul_CM [norm too large]");
  vn = ((ln>>1)-4)>>2;
  z1 = ellwpseries(e, 0, ln);
  z2 = ser_unscale(z1, n);
  p0 = gen_0; p1 = gen_1;
  q0 = gen_1; q1 = gen_0;
  do
  {
    GEN p2,q2, ss = gen_0;
    do
    {
      long ep = (-valser(z2)) >> 1;
      ss = gadd(ss, gmul(gel(z2,2), pol_xnall(ep, 0)));
      z2 = gsub(z2, gmul(gel(z2,2), gpowgs(z1, ep)));
    }
    while (valser(z2) <= 0);
    p2 = gadd(p0, gmul(ss,p1)); p0 = p1; p1 = p2;
    q2 = gadd(q0, gmul(ss,q1)); q0 = q1; q1 = q2;
    if (!signe(z2)) break;
    z2 = ginv(z2);
  }
  while (degpol(p1) < vn);
  if (degpol(p1) > vn || signe(z2))
    pari_err_TYPE("ellmul [not a complex multiplication]", n);
  q1p = RgX_deriv(q1);
  b2ov12 = gdivgu(ell_get_b2(e), 12);
  grdx = gadd(gel(P,1), b2ov12); /* x(P) + b2/12 */
  q1 = poleval(q1, grdx);
  if (gequal0(q1)) return ellinf();

  p1p = RgX_deriv(p1);
  p1 = poleval(p1, grdx);
  p1p = poleval(p1p, grdx);
  q1p = poleval(q1p, grdx);

  x = gdiv(p1,q1);
  y = gdiv(gsub(gmul(p1p,q1), gmul(p1,q1p)), gmul(n,gsqr(q1)));
  x = gsub(x, b2ov12);
  y = gsub( gmul(ec_dmFdy_evalQ(e,P), y), ec_h_evalx(e,x));
  return mkvec2(x, gmul2n(y,-1));
}

static GEN
_sqr(void *e, GEN x) { return elladd((GEN)e, x, x); }
static GEN
_mul(void *e, GEN x, GEN y) { return elladd((GEN)e, x, y); }

static GEN
ellffmul(GEN E, GEN P, GEN n)
{
  GEN fg = ellff_get_field(E);
  if (typ(fg)==t_FFELT)
    return FF_ellmul(E, P, n);
  else
  {
    pari_sp av = avma;
    GEN p = fg, e = ellff_get_a4a6(E), Q;
    GEN Pp = FpE_changepointinv(RgE_to_FpE(P, p), gel(e,3), p);
    GEN Qp = FpE_mul(Pp, n, gel(e,1), p);
    Q = FpE_to_mod(FpE_changepoint(Qp, gel(e,3), p), p);
    return gerepileupto(av, Q);
  }
}
/* [n] z, n integral */
static GEN
ellmul_Z(GEN e, GEN z, GEN n)
{
  long s;
  if (ell_is_inf(z)) return ellinf();
  if (ell_over_Fq(e)) return ellffmul(e,z,n);
  s = signe(n);
  if (!s) return ellinf();
  if (s < 0) z = ellneg_i(e,z);
  if (is_pm1(n)) return z;
  return gen_pow(z, n, (void*)e, &_sqr, &_mul);
}

/* x a t_REAL, try to round it to an integer */
enum { OK, LOW_PREC, NO };
static long
myroundr(GEN *px)
{
  GEN x = *px;
  long e;
  if (bit_prec(x) - expo(x) < 5) return LOW_PREC;
  *px = grndtoi(x, &e);
  if (e >= -5) return NO;
  return OK;
}

/* E has CM by Q, t_COMPLEX or t_QUAD. Return q such that E has CM by Q/q
 * or gen_1 (couldn't find q > 1)
 * or NULL (doesn't have CM by Q) */
static GEN
CM_factor(GEN E, GEN Q)
{
  GEN w, tau, D, v, x, y, F, dF, q, r, fk, fkb, fkc;
  long prec;

  if (ell_get_type(E) != t_ELL_Q) return gen_1;
  switch(typ(Q))
  {
    case t_COMPLEX:
      D = utoineg(4);
      v = gel(Q,2);
      break;
    case t_QUAD:
      D = quad_disc(Q);
      v = gel(Q,3);
      break;
    default:
      return NULL; /*-Wall*/
  }
  /* disc Q = v^2 D, D < 0 fundamental */
  w = ellR_omega(E, DEFAULTPREC + nbits2nlong(expi(D)));
  tau = gdiv(gel(w,2), gel(w,1));
  prec = precision(tau);
  /* disc tau = -4 k^2 (Im tau)^2 for some integral k
   * Assuming that E has CM by Q, then disc Q / disc tau = f^2 is a square.
   * Compute f*k */
  x = gel(tau,1);
  y = gel(tau,2); /* tau = x + Iy */
  fk = gmul(gdiv(v, gmul2n(y, 1)), sqrtr_abs(itor(D, prec)));
  switch(myroundr(&fk))
  {
    case NO: return NULL;
    case LOW_PREC: return gen_1;
  }
  fk = absi_shallow(fk);

  fkb = gmul(fk, gmul2n(x,1));
  switch(myroundr(&fkb))
  {
    case NO: return NULL;
    case LOW_PREC: return gen_1;
  }

  fkc = gmul(fk, cxnorm(tau));
  switch(myroundr(&fkc))
  {
    case NO: return NULL;
    case LOW_PREC: return gen_1;
  }

  /* tau is a root of fk (X^2 - b X + c) \in Z[X],  */
  F = Q_primpart(mkvec3(fk, fkb, fkc));
  dF = qfb_disc3(gel(F,1), gel(F,2), gel(F,3));
  /* = disc tau, E has CM by orders of disc dF q^2, all q */
  q = dvmdii(dF, D, &r);
  if (r != gen_0 || !Z_issquareall(q, &q)) return NULL;
  /* disc(Q) = disc(tau) (v / q)^2 */
  v = dvmdii(absi_shallow(v), q, &r);
  if (r != gen_0) return NULL;
  return is_pm1(v)? gen_1: v; /* E has CM by Q/q: [Q] = [q] o [Q/q] */
}

/* [a + w] z, a integral, w pure imaginary */
static GEN
ellmul_CM_aux(GEN e, GEN z, GEN a, GEN w)
{
  GEN A, B, q;
  if (typ(a) != t_INT) pari_err_TYPE("ellmul_Z",a);
  q = CM_factor(e, w);
  if (!q) pari_err_TYPE("ellmul [not a complex multiplication]",w);
  if (q != gen_1) w = gdiv(w, q);
  /* compute [a + q w] z, z has CM by w */
  if (typ(w) == t_QUAD && is_pm1(gel(gel(w,1), 3)))
  { /* replace w by w - u, u in Z, so that N(w-u) is minimal
     * N(w - u) = N w - Tr w u + u^2, minimal for u = Tr w / 2 */
    GEN u = gtrace(w);
    if (typ(u) != t_INT) pari_err_TYPE("ellmul_CM",w);
    u = shifti(u, -1);
    if (signe(u))
    {
      w = gsub(w, u);
      a = addii(a, mulii(q,u));
    }
    /* [a + w]z = [(a + qu)] z + [q] [(w - u)] z */
  }
  A = ellmul_Z(e,z,a);
  B = ellmul_CM(e,z,w);
  if (q != gen_1) B = ellmul_Z(e, B, q);
  return elladd(e, A, B);
}
GEN
ellmul(GEN e, GEN z, GEN n)
{
  pari_sp av = avma;

  checkell(e); checkellpt(z);
  if (ell_is_inf(z)) return ellinf();
  switch(typ(n))
  {
    case t_INT: return gerepilecopy(av, ellmul_Z(e,z,n));
    case t_QUAD: {
      GEN pol = gel(n,1), a = gel(n,2), b = gel(n,3);
      if (signe(gel(pol,2)) < 0) pari_err_TYPE("ellmul_CM",n); /* disc > 0 ? */
      return gerepileupto(av, ellmul_CM_aux(e,z,a,mkquad(pol, gen_0,b)));
    }
    case t_COMPLEX: {
      GEN a = gel(n,1), b = gel(n,2);
      return gerepileupto(av, ellmul_CM_aux(e,z,a,mkcomplex(gen_0,b)));
    }
  }
  pari_err_TYPE("ellmul (non integral, non CM exponent)",n);
  return NULL; /* LCOV_EXCL_LINE */
}

/********************************************************************/
/**                                                                **/
/**                       Periods                                  **/
/**                                                                **/
/********************************************************************/

/* References:
  The complex AGM, periods of elliptic curves over C and complex elliptic logarithms
  John E. Cremona, Thotsaphon Thongjunthug, arXiv:1011.0914
*/

static GEN
ellomega_agm(GEN a, GEN b, GEN c, long prec)
{
  GEN pi = mppi(prec), mIpi = mkcomplex(gen_0, negr(pi));
  GEN Mac = agm(a,c,prec), Mbc = agm(b,c,prec);
  retmkvec2(gdiv(pi, Mac), gdiv(mIpi, Mbc));
}

static GEN
ellomega_cx(GEN E, long prec)
{
  pari_sp av = avma;
  GEN roots = ellR_roots(E, prec + EXTRAPREC64);
  GEN d1=gel(roots,4), d2=gel(roots,5), d3=gel(roots,6);
  GEN a = gsqrt(d3,prec), b = gsqrt(d1,prec), c = gsqrt(d2,prec);
  return gerepileupto(av, ellomega_agm(a,b,c,prec));
}

/* return [w1,w2] for E / R; w1 > 0 is real.
 * If e.disc > 0, w2 = -I r; else w2 = w1/2 - I r, for some real r > 0.
 * => tau = w1/w2 is in upper half plane */
static GEN
doellR_omega(GEN E, long prec)
{
  pari_sp av = avma;
  GEN roots, d2, z, a, b, c;
  if (ellR_get_sign(E) >= 0) return ellomega_cx(E,prec);
  roots = ellR_roots(E,prec + EXTRAPREC64);
  d2 = gel(roots,5);
  z = gsqrt(d2,prec); /* imag(e1-e3) > 0, so that b > 0*/
  a = gel(z,1); /* >= 0 */
  b = gel(z,2);
  c = gabs(z, prec);
  z = ellomega_agm(a,b,c,prec);
  return gerepilecopy(av, mkvec2(gel(z,1),gmul2n(gadd(gel(z,1),gel(z,2)),-1)));
}
static GEN
doellR_eta(GEN E, long prec)
{ GEN w = ellR_omega(E, prec + EXTRAPREC64); return elleta(w, prec); }

GEN
ellR_omega(GEN E, long prec)
{ return obj_checkbuild_realprec(E, R_PERIODS, &doellR_omega, prec); }
GEN
ellR_eta(GEN E, long prec)
{ return obj_checkbuild_realprec(E, R_ETA, &doellR_eta, prec); }
GEN
ellR_roots(GEN E, long prec)
{ return obj_checkbuild_realprec(E, R_ROOTS, &doellR_roots, prec); }

GEN
ellR_area(GEN E, long prec)
{
  pari_sp av = avma;
  GEN w, w1, w2, a,b,c,d;
  w = ellR_omega(E, prec);
  w1 = gel(w,1); a = real_i(w1); b = imag_i(w1);
  w2 = gel(w,2); c = real_i(w2); d = imag_i(w2);
  return gerepileupto(av, gabs(gsub(gmul(a,d),gmul(b,c)), prec));
}

/********************************************************************/
/**                                                                **/
/**                       ELLIPTIC FUNCTIONS                       **/
/**                                                                **/
/********************************************************************/
/* P = [x,0] is 2-torsion on y^2 = g(x). Return w1/2, (w1+w2)/2, or w2/2
 * depending on whether x is closest to e1,e2, or e3, the 3 complex root of g */
static GEN
zell_closest_0(GEN om, GEN x, GEN ro)
{
  GEN e1 = gel(ro,1), e2 = gel(ro,2), e3 = gel(ro,3);
  GEN d1 = gnorm(gsub(x,e1));
  GEN d2 = gnorm(gsub(x,e2));
  GEN d3 = gnorm(gsub(x,e3));
  GEN z = gel(om,2);
  if (gcmp(d1, d2) <= 0)
  { if (gcmp(d1, d3) <= 0) z = gel(om,1); }
  else
  { if (gcmp(d2, d3)<=0) z = gadd(gel(om,1),gel(om,2)); }
  return gmul2n(z, -1);
}

static GEN
zellcx(GEN E, GEN P, long prec)
{
  GEN R = ellR_roots(E, prec+EXTRAPREC64);
  GEN x0 = gel(P,1), y0 = ec_dmFdy_evalQ(E,P);
  if (gequal0(y0))
    return zell_closest_0(ellomega_cx(E,prec),x0,R);
  else
  {
    GEN e2 = gel(R,2), e3 = gel(R,3), d2 = gel(R,5), d3 = gel(R,6);
    GEN a = gsqrt(d2,prec), b = gsqrt(d3,prec);
    GEN r = gsqrt(gdiv(gsub(x0,e3), gsub(x0,e2)),prec);
    GEN t = gdiv(gneg(y0), gmul2n(gmul(r,gsub(x0,e2)),1));
    GEN ar = real_i(a), br = real_i(b), ai = imag_i(a), bi = imag_i(b);
    /* |a+b| < |a-b| */
    if (gcmp(gmul(ar,br), gneg(gmul(ai,bi))) < 0) b = gneg(b);
    return zellagmcx(a,b,r,t,prec);
  }
}

/* Assume E/R, disc E < 0, and P \in E(R) ==> z \in R */
static GEN
zellrealneg(GEN E, GEN P, long prec)
{
  GEN x0 = gel(P,1), y0 = ec_dmFdy_evalQ(E,P);
  if (gequal0(y0)) return gmul2n(gel(ellR_omega(E,prec),1),-1);
  else
  {
    GEN R = ellR_roots(E, prec+EXTRAPREC64);
    GEN d2 = gel(R,5), e3 = gel(R,3);
    GEN a = gsqrt(d2,prec);
    GEN z = gsqrt(gsub(x0,e3), prec);
    GEN ar = real_i(a), zr = real_i(z), ai = imag_i(a), zi = imag_i(z);
    GEN t = gdiv(gneg(y0), gmul2n(gnorm(z),1));
    GEN r2 = ginv(gsqrt(gaddsg(1,gdiv(gmul(ai,zi),gmul(ar,zr))),prec));
    return zellagmcx(ar,gabs(a,prec),r2,gmul(t,r2),prec);
  }
}

/* Assume E/R, disc E > 0, and P \in E(R) */
static GEN
zellrealpos(GEN E, GEN P, long prec)
{
  GEN R = ellR_roots(E, prec+EXTRAPREC64);
  GEN d2,d3,e1,e2,e3, a,b, x0 = gel(P,1), y0 = ec_dmFdy_evalQ(E,P);
  if (gequal0(y0)) return zell_closest_0(ellR_omega(E,prec), x0,R);
  e1 = gel(R,1);
  e2 = gel(R,2);
  e3 = gel(R,3);
  d2 = gel(R,5);
  d3 = gel(R,6);
  a = gsqrt(d2,prec);
  b = gsqrt(d3,prec);
  if (gcmp(x0,e1)>0) {
    GEN r = gsqrt(gdiv(gsub(x0,e3), gsub(x0,e2)),prec);
    GEN t = gdiv(gneg(y0), gmul2n(gmul(r,gsub(x0,e2)),1));
    return zellagmcx(a,b,r,t,prec);
  } else {
    GEN om = ellR_omega(E,prec);
    GEN r = gdiv(a,gsqrt(gsub(e1,x0),prec));
    GEN t = gdiv(gmul(r,y0),gmul2n(gsub(x0,e3),1));
    return gsub(zellagmcx(a,b,r,t,prec),gmul2n(gel(om,2),-1));
  }
}

static void
ellQp_P2t_err(GEN E, GEN z)
{
  if (typ(ellQp_u(E,1)) == t_POLMOD)
    pari_err_IMPL("ellpointtoz when u not in Qp");
  pari_err_DOMAIN("ellpointtoz", "point", "not on", strtoGENstr("E"),z);
}
static GEN
get_r0(GEN E, long prec)
{
  GEN b2 = ell_get_b2(E), e1 = ellQp_root(E, prec);
  return gadd(e1,gmul2n(b2,-2));
}
static GEN
ellQp_P2t(GEN E, GEN P, long prec)
{
  pari_sp av = avma;
  GEN a, b, ab, c0, r0, ar, r, x, delta, x1, y1, t, u, q;
  long vq, vt, Q, R;
  if (ell_is_inf(P)) return gen_1;
  ab = ellQp_ab(E, prec); a = gel(ab,1); b = gel(ab,2);
  u = ellQp_u(E, prec);
  q = ellQp_q(E, prec);
  x = gel(P,1);
  r0 = get_r0(E, prec);
  c0 = gadd(x, gmul2n(r0,-1));
  if (typ(c0) != t_PADIC || !is_scalar_t(typ(gel(P,2))))
    pari_err_TYPE("ellpointtoz",P);
  r = gsub(a,b);
  ar = gmul(a, r);
  if (gequal0(c0))
  {
    x1 = Qp_sqrt(gneg(ar));
    if (!x1) ellQp_P2t_err(E,P);
  }
  else
  {
    delta = gdiv(ar, gsqr(c0));
    t = Qp_sqrt(gsubsg(1,gmul2n(delta,2)));
    if (!t) ellQp_P2t_err(E,P);
    x1 = gmul(gmul2n(c0,-1), gaddsg(1,t));
  }
  y1 = gsubsg(1, gdiv(ar, gsqr(x1)));
  if (gequal0(y1))
  {
    y1 = Qp_sqrt(gmul(x1, gmul(gadd(x1, a), gadd(x1, r))));
    if (!y1) ellQp_P2t_err(E,P);
  }
  else
    y1 = gdiv(gmul2n(ec_dmFdy_evalQ(E,P), -1), y1);
  Qp_descending_Landen(ellQp_AGM(E,prec), &x1,&y1);

  t = gmul(u, gmul2n(y1,1)); /* 2u y_oo */
  t = gdiv(gsub(t, x1), gadd(t, x1));
  /* Reduce mod q^Z: we want 0 <= v(t) < v(q) */
  if (typ(t) == t_PADIC)
    vt = valp(t);
  else
    vt = valp(gnorm(t)) / 2; /* v(t) = v(Nt) / (e*f) */
  vq = valp(q); /* > 0 */
  Q = vt / vq; R = vt % vq; if (R < 0) Q--;
  if (Q) t = gdiv(t, gpowgs(q,Q));
  if (padicprec_relative(t) > prec) t = gprec(t, prec);
  return gerepileupto(av, t);
}

static GEN
ellQp_t2P(GEN E, GEN t, long prec)
{
  pari_sp av = avma;
  GEN AB, A, R, x0,x1, y0,y1, u, u2, r0, s0, ar;
  long v;
  if (gequal1(t)) return ellinf();

  AB = ellQp_AGM(E,prec); A = gel(AB,1); R = gel(AB,3); v = itos(gel(AB,4));
  u = ellQp_u(E,prec);
  u2= ellQp_u2(E,prec);
  x1 = gdiv(t, gmul(u2, gsqr(gsubsg(1,t))));
  y1 = gdiv(gmul(x1,gaddsg(1,t)), gmul(gmul2n(u,1),gsubsg(1,t)));
  Qp_ascending_Landen(AB, &x1,&y1);
  r0 = get_r0(E, prec);

  ar = gmul(gel(A,1), gel(R,1)); setvalp(ar, valp(ar)+v);
  x0 = gsub(gadd(x1, gdiv(ar, x1)), gmul2n(r0,-1));
  s0 = gmul2n(ec_h_evalx(E, x0), -1);
  y0 = gsub(gmul(y1, gsubsg(1, gdiv(ar,gsqr(x1)))), s0);
  return gerepilecopy(av, mkvec2(x0,y0));
}

static GEN
zell_i(GEN e, GEN z, long prec)
{
  GEN t;
  long s;
  (void)ellR_omega(e, prec); /* type checking */
  if (ell_is_inf(z)) return gen_0;
  s = ellR_get_sign(e);
  if (s && typ(gel(z,1))!=t_COMPLEX && typ(gel(z,2))!=t_COMPLEX)
    t = (s < 0)? zellrealneg(e,z,prec): zellrealpos(e,z,prec);
  else
    t = zellcx(e,z,prec);
  return t;
}
static GEN ellnfembed(GEN E, long prec);
static GEN ellpointnfembed(GEN E, GEN P, long prec);
static void ellnfembed_free(GEN L);
GEN
zell(GEN E, GEN P, long prec)
{
  pari_sp av = avma;
  checkell(E); checkellpt(P);
  switch(ell_get_type(E))
  {
    case t_ELL_Qp:
      prec = minss(ellQp_get_prec(E), padicprec_relative(P));
      return ellQp_P2t(E, P, prec);
    case t_ELL_NF:
    {
      GEN Ee = ellnfembed(E, prec), Pe = ellpointnfembed(E, P, prec);
      long i, l = lg(Pe);
      for (i = 1; i < l; i++) gel(Pe,i) = zell_i(gel(Ee,i), gel(Pe,i), prec);
      ellnfembed_free(Ee); return gerepilecopy(av, Pe);
    }
    case t_ELL_Q: break;
    case t_ELL_Rg: break;
    default: pari_err_TYPE("ellpointtoz", E);
  }
  return gerepileupto(av, zell_i(E, P, prec));
}

enum period_type { t_PER_W, t_PER_WETA, t_PER_ELL };
/* normalization / argument reduction for ellptic functions */
typedef struct {
  enum period_type type;
  GEN in; /* original input */
  GEN w1,w2,tau; /* original basis for L = <w1,w2> = w2 <1,tau> */
  GEN W1,W2,Tau; /* new basis for L = <W1,W2> = W2 <1,tau> */
  GEN a,b,c,d; /* t_INT; tau in F = h/Sl2, tau = g.t, g=[a,b;c,d] in SL(2,Z) */
  GEN z,Z; /* z/w2 defined mod <1,tau>, Z = z/w2 + x*tau+y reduced mod <1,tau>*/
  GEN x,y; /* t_INT */
  int swap; /* 1 if we swapped w1 and w2 */
  int some_q_is_real; /* exp(2iPi g.tau) for some g \in SL(2,Z) */
  int some_z_is_real; /* z + xw1 + yw2 is real for some x,y \in Z */
  int some_z_is_pure_imag; /* z + xw1 + yw2 in i*R */
  int q_is_real; /* exp(2iPi tau) \in R */
  int abs_u_is_1; /* |exp(2iPi Z)| = 1 */
  long prec; /* precision(Z) */
  long prec0; /* required precision for result */
} ellred_t;

/* compute g in SL_2(Z), g.t is in the usual
   fundamental domain. Internal function no check, no garbage. */
static void
set_gamma(GEN *pt, GEN *pa, GEN *pb, GEN *pc, GEN *pd)
{
  GEN a, b, c, d, t, t0 = *pt, run = dbltor(1. - 1e-8);
  long e = gexpo(gel(t0,2));
  if (e < 0) t0 = gprec_wensure(t0, precision(t0)+nbits2extraprec(-e));
  t = t0;
  a = d = gen_1;
  b = c = gen_0;
  for(;;)
  {
    GEN m, n = ground(gel(t,1));
    if (signe(n))
    { /* apply T^n */
      t = gsub(t,n);
      a = subii(a, mulii(n,c));
      b = subii(b, mulii(n,d));
    }
    m = cxnorm(t); if (gcmp(m,run) > 0) break;
    t = gneg_i(gdiv(conj_i(t), m)); /* apply S */
    togglesign_safe(&c); swap(a,c);
    togglesign_safe(&d); swap(b,d);
  }
  if (e < 0 && (signe(b) || signe(c))) *pt = t0;
  *pa = a; *pb = b; *pc = c; *pd = d;
}
/* Im z > 0. Return U.z in PSl2(Z)'s standard fundamental domain.
 * Set *pU to U. */
GEN
cxredsl2_i(GEN z, GEN *pU, GEN *czd)
{
  GEN a,b,c,d;
  set_gamma(&z, &a, &b, &c, &d);
  *pU = mkmat2(mkcol2(a,c), mkcol2(b,d));
  *czd = gadd(gmul(c,z), d);
  return gdiv(gadd(gmul(a,z), b), *czd);
}
GEN
cxredsl2(GEN t, GEN *pU)
{
  pari_sp av = avma;
  GEN czd;
  t = cxredsl2_i(t, pU, &czd);
  return gc_all(av, 2, &t, pU);
}

/* swap w1, w2 so that Im(t := w1/w2) > 0. Set tau = representative of t in
 * the standard fundamental domain, and g in Sl_2, such that tau = g.t */
static void
red_modSL2(ellred_t *T, long prec)
{
  long s, p;
  T->tau = gdiv(T->w1,T->w2);
  if (isintzero(real_i(T->tau))) T->some_q_is_real = 1;
  s = gsigne(imag_i(T->tau));
  if (!s) pari_err_DOMAIN("elliptic function", "det(w1,w2)", "=", gen_0,
                          mkvec2(T->w1,T->w2));
  T->swap = (s < 0);
  if (T->swap) { swap(T->w1, T->w2); T->tau = ginv(T->tau); }
  p = precision(T->tau); T->prec0 = p? p: prec;
  set_gamma(&T->tau, &T->a, &T->b, &T->c, &T->d);
  /* update lattice */
  p = precision(T->tau);
  if (p)
  {
    T->w1 = gprec_wensure(T->w1, p);
    T->w2 = gprec_wensure(T->w2, p);
  }
  T->W1 = gadd(gmul(T->a,T->w1), gmul(T->b,T->w2));
  T->W2 = gadd(gmul(T->c,T->w1), gmul(T->d,T->w2));
  T->Tau = gdiv(T->W1, T->W2);
  if (isintzero(real_i(T->Tau))) T->some_q_is_real = T->q_is_real = 1;
  p = precision(T->Tau); T->prec = p? p: prec;
}
/* is z real or pure imaginary ? */
static void
check_complex(GEN z, int *real, int *imag)
{
  if (typ(z) != t_COMPLEX)      { *real = 1; *imag = 0; }
  else if (isintzero(gel(z,1))) { *real = 0; *imag = 1; }
  else *real = *imag = 0;
}
static void
reduce_z(GEN z, ellred_t *T)
{
  GEN x, Z;
  long p, e;
  switch(typ(z))
  {
    case t_INT: case t_REAL: case t_FRAC: case t_COMPLEX: break;
    case t_QUAD:
      z = isexactzero(gel(z,2))? gel(z,1): quadtofp(z, T->prec);
      break;
    default: pari_err_TYPE("reduction mod 2-dim lattice (reduce_z)", z);
  }
  Z = gdiv(z, T->W2);
  T->z = z;
  x = gdiv(imag_i(Z), imag_i(T->Tau));
  T->x = grndtoi(x, &e); /* |Im(Z - x*Tau)| <= Im(Tau)/2 */
  /* Avoid Im(Z) << 0; take 0 <= Im(Z - x*Tau) < Im(Tau) instead.
   * Leave round when Im(Z - x*Tau) ~ 0 to allow detecting Z in <1,Tau>
   * at the end */
  if (e > -10) T->x = gfloor(x);
  if (signe(T->x)) Z = gsub(Z, gmul(T->x,T->Tau));
  T->y = ground(real_i(Z));/* |Re(Z - y)| <= 1/2 */
  if (signe(T->y)) Z = gsub(Z, T->y);
  T->abs_u_is_1 = (typ(Z) != t_COMPLEX);
  /* Z = - y - x tau + z/W2, x,y integers */
  check_complex(z, &(T->some_z_is_real), &(T->some_z_is_pure_imag));
  if (!T->some_z_is_real && !T->some_z_is_pure_imag)
  {
    int W2real, W2imag;
    check_complex(T->W2,&W2real,&W2imag);
    if (W2real)
      check_complex(Z, &(T->some_z_is_real), &(T->some_z_is_pure_imag));
    else if (W2imag)
      check_complex(Z, &(T->some_z_is_pure_imag), &(T->some_z_is_real));
  }
  p = precision(Z);
  if (gequal0(Z) || (p && gexpo(Z) < 5 - prec2nbits(p))) Z = NULL; /*z in L*/
  if (p && p < T->prec) T->prec = p;
  T->Z = Z;
}
/* return x.eta1 + y.eta2 */
static GEN
eta_correction(ellred_t *T, GEN eta)
{
  GEN y1 = NULL, y2 = NULL;
  if (signe(T->x)) y1 = gmul(T->x, gel(eta,1));
  if (signe(T->y)) y2 = gmul(T->y, gel(eta,2));
  if (!y1) return y2? y2: gen_0;
  return y2? gadd(y1, y2): y1;
}
/* e is either
 * - [w1,w2]
 * - [[w1,w2],[eta1,eta2]]
 * - an ellinit structure */
static void
compute_periods(ellred_t *T, GEN z, long prec)
{
  GEN w, e;
  T->q_is_real = 0;
  T->some_q_is_real = 0;
  switch(T->type)
  {
    case t_PER_ELL:
    {
      long pr, p = prec;
      if (z && (pr = precision(z))) p = pr;
      e = T->in;
      w = ellR_omega(e, p);
      T->some_q_is_real = T->q_is_real = 1;
      break;
    }
    case t_PER_W:
      w = T->in; break;
    default: /*t_PER_WETA*/
      w = gel(T->in,1); break;
  }
  T->w1 = gel(w,1);
  T->w2 = gel(w,2);
  red_modSL2(T, prec);
  if (z) reduce_z(z, T);
}
static int
check_periods(GEN e, ellred_t *T)
{
  GEN w1;
  if (typ(e) != t_VEC) return 0;
  T->in = e;
  switch(lg(e))
  {
    case 17:
      T->type = t_PER_ELL;
      break;
    case 3:
      w1 = gel(e,1);
      if (typ(w1) != t_VEC)
        T->type = t_PER_W;
      else
      {
        if (lg(w1) != 3) return 0;
        T->type = t_PER_WETA;
      }
      break;
    default: return 0;
  }
  return 1;
}
static int
get_periods(GEN e, GEN z, ellred_t *T, long prec)
{
  if (!check_periods(e, T)) return 0;
  compute_periods(T, z, prec); return 1;
}

/* 2iPi/x, more efficient when x pure imaginary */
static GEN
PiI2div(GEN x, long prec) { return gdiv(Pi2n(1, prec), mulcxmI(x)); }
/* (2iPi/W2)^k E_k(W1/W2) */
static GEN
_elleisnum(ellred_t *T, long k)
{
  GEN z = gmul(cxEk(T->Tau, k, T->prec), gpowgs(PiI2div(T->W2, T->prec), k));
  return cxtoreal(z);
}

/* Return (2iPi)^k E_k(L) = (2iPi/w2)^k E_k(tau), with L = <w1,w2>, k > 0 even
 * E_k(tau) = 1 + 2/zeta(1-k) * sum(n>=1, n^(k-1) q^n/(1-q^n))
 * If flag is != 0 and k=4 or 6, compute g2 = E4/12 or g3 = -E6/216 resp. */
GEN
elleisnum(GEN om, long k, long flag, long prec)
{
  pari_sp av = avma;
  GEN y;
  ellred_t T;

  if (k<=0) pari_err_DOMAIN("elleisnum", "k", "<=", gen_0, stoi(k));
  if (k&1) pari_err_DOMAIN("elleisnum", "k % 2", "!=", gen_0, stoi(k));
  if (!get_periods(om, NULL, &T, prec)) pari_err_TYPE("elleisnum",om);
  y = _elleisnum(&T, k);
  if (k==2 && signe(T.c))
  {
    GEN a = gmul(Pi2n(1,T.prec), mului(12, T.c));
    y = gsub(y, mulcxI(gdiv(a, gmul(T.w2, T.W2))));
  }
  else if (k==4 && flag) y = gdivgu(y,  12);
  else if (k==6 && flag) y = gdivgs(y,-216);
  return gerepilecopy(av, gprec_wtrunc(y, T.prec0));
}

/* return quasi-periods attached to [T->W1,T->W2] */
static GEN
_elleta(ellred_t *T)
{
  GEN y1, y2, e2 = gdivgs(_elleisnum(T,2), -12);
  y2 = gmul(T->W2, e2);
  y1 = gsub(gmul(T->W1,e2), PiI2div(T->W2, T->prec));
  retmkvec2(y1, y2);
}

/* compute eta1, eta2 */
GEN
elleta(GEN om, long prec)
{
  pari_sp av = avma;
  GEN y1, y2, E2, pi;
  ellred_t T;

  if (!check_periods(om, &T))
  {
    pari_err_TYPE("elleta",om);
    return NULL;/*LCOV_EXCL_LINE*/
  }
  if (T.type == t_PER_ELL) return ellR_eta(om, prec);

  compute_periods(&T, NULL, prec);
  prec = T.prec;
  pi = mppi(prec);
  E2 = cxEk(T.Tau, 2, prec); /* E_2(Tau) */
  if (signe(T.c))
  {
    GEN u = gdiv(T.w2, T.W2);
    /* E2 := u^2 E2 + 6iuc/pi = E_2(tau) */
    E2 = gadd(gmul(gsqr(u), E2), mulcxI(gdiv(gmul(mului(6,T.c), u), pi)));
  }
  y2 = gdiv(gmul(E2, sqrr(pi)), gmulsg(3, T.w2));
  if (T.swap)
  {
    y1 = y2;
    y2 = gadd(gmul(T.tau,y1), PiI2div(T.w2, prec));
  }
  else
    y1 = gsub(gmul(T.tau,y2), PiI2div(T.w2, prec));
  switch(typ(T.w1))
  {
    case t_INT: case t_FRAC: case t_REAL:
      y1 = real_i(y1);
  }
  return gerepilecopy(av, mkvec2(y1,y2));
}
GEN
ellperiods(GEN w, long flag, long prec)
{
  pari_sp av = avma;
  ellred_t T;
  if (!get_periods(w, NULL, &T, prec)) pari_err_TYPE("ellperiods",w);
  switch(flag)
  {
    case 0: return gerepilecopy(av, mkvec2(T.W1, T.W2));
    case 1: return gerepilecopy(av, mkvec2(mkvec2(T.W1, T.W2), _elleta(&T)));
    default: pari_err_FLAG("ellperiods");
             return NULL;/*LCOV_EXCL_LINE*/
  }
}

/* 2Pi Im(z)/log(2) */
static double
get_toadd(GEN z) { return (2*M_PI/M_LN2)*gtodouble(imag_i(z)); }

/* computes the numerical value of wp(z | L), L = om1 Z + om2 Z
 * return NULL if z in L.  If flall=1, compute also wp' */
static GEN
ellwpnum_all(GEN e, GEN z, long flall, long prec)
{
  long toadd;
  pari_sp av = avma, av1;
  GEN q, u, y, yp, u1, u2, qn;
  ellred_t T;
  int simple_case;

  if (!get_periods(e, z, &T, prec)) pari_err_TYPE("ellwp",e);
  if (!T.Z) return NULL;
  prec = T.prec;

  /* Now L,Z normalized to <1,tau>. Z in fund. domain of <1, tau> */
  q = expIPiC(gmul2n(T.Tau,1), prec);
  u = expIPiC(gmul2n(T.Z,1), prec);
  u1 = gsubsg(1,u);
  u2 = gsqr(u1); /* (1-u)^2 = -4u sin^2(Pi Z) */
  if (gequal0(gnorm(u2))) return NULL; /* possible if loss of accuracy */
  y = gdiv(u,u2); /* -1/4(sin^2(Pi Z)) */
  if (T.abs_u_is_1) y = real_i(y);
  simple_case = T.abs_u_is_1 && T.q_is_real;
  y = gadd(mkfrac(gen_1, utoipos(12)), y);
  yp = flall? gen_0: NULL;
  toadd = (long)ceil(get_toadd(T.Z));

  av1 = avma; qn = q;
  for(;;)
  { /* y += u q^n [ 1/(1-q^n u)^2 + 1/(q^n-u)^2 ] - 2q^n /(1-q^n)^2 */
    /* analogous formula for yp */
    GEN yadd, ypadd = NULL;
    GEN qnu = gmul(qn,u); /* q^n u */
    GEN a = gsubsg(1,qnu);/* 1 - q^n u */
    GEN a2 = gsqr(a);     /* (1 - q^n u)^2 */
    if (yp) ypadd = gdiv(gaddsg(1,qnu),gmul(a,a2));
    if (simple_case) /* conj(u) = 1/u: formula simplifies */
      yadd = gmul2n(real_i(gdiv(u,a2)), 1);
    else
    {
      GEN b = gsub(qn,u);/* q^n - u */
      GEN b2 = gsqr(b);  /* (q^n - u)^2 */
      yadd = gmul(u, gadd(ginv(a2),ginv(b2)));
      if (yp) ypadd = gadd(ypadd, gdiv(gadd(qn,u),gmul(b,b2)));
    }
    yadd = gsub(yadd, gmul2n(ginv(gsqr(gsubsg(1,qn))), 1));
    y = gadd(y, gmul(qn,yadd));
    if (yp) yp = gadd(yp, gmul(qn,ypadd));

    qn = gmul(q,qn);
    if (gexpo(qn) <= - prec2nbits(prec) - 5 - toadd) break;
    if (gc_needed(av1,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ellwp");
      gerepileall(av1, flall? 3: 2, &y, &qn, &yp);
    }
  }
  if (yp)
  {
    if (simple_case) yp = gsub(yp, conj_i(gmul(yp,gsqr(u))));
    yp = gadd(yp, gdiv(gaddsg(1,u), gmul(u1,u2)));
  }

  u1 = PiI2div(T.W2, prec);
  u2 = gsqr(u1);
  y = gmul(u2,y); /* y *= (2i pi / w2)^2 */
  if (T.some_q_is_real && (T.some_z_is_real || T.some_z_is_pure_imag))
    y = real_i(y);
  if (yp)
  {
    yp = gmul(u, gmul(gmul(u1,u2),yp));/* yp *= u (2i pi / w2)^3 */
    if (T.some_q_is_real)
    {
      if (T.some_z_is_real) yp = real_i(yp);
      else if (T.some_z_is_pure_imag) yp = mkcomplex(gen_0, imag_i(yp));
    }
    y = mkvec2(y, yp);
  }
  return gerepilecopy(av, gprec_wtrunc(y, T.prec0));
}
static GEN
ellwpseries_aux(GEN c4, GEN c6, long v, long PRECDL)
{
  long i, k, l;
  pari_sp av;
  GEN _1, t, res = cgetg(PRECDL+2,t_SER), *P = (GEN*)(res + 2);

  res[1] = evalsigne(1) | _evalvalser(-2) | evalvarn(v);
  if (!PRECDL) { setsigne(res,0); return res; }

  for (i=1; i<PRECDL; i+=2) P[i]= gen_0;
  _1 = Rg_get_1(c4);
  switch(PRECDL)
  {
    default:P[6] = gdivgu(c6,6048);
    case 6:
    case 5: P[4] = gdivgu(c4, 240);
    case 4:
    case 3: P[2] = gmul(_1,gen_0);
    case 2:
    case 1: P[0] = _1;
  }
  if (PRECDL <= 8) return res;
  av = avma;
  P[8] = gerepileupto(av, gdivgu(gsqr(P[4]), 3));
  for (k=5; (k<<1) < PRECDL; k++)
  {
    av = avma;
    t = gmul(P[4], P[(k-2)<<1]);
    for (l=3; (l<<1) < k; l++) t = gadd(t, gmul(P[l<<1], P[(k-l)<<1]));
    t = gmul2n(t, 1);
    if ((k & 1) == 0) t = gadd(gsqr(P[k]), t);
    if (k % 3 == 2)
      t = gdivgu(gmulsg(3, t), (k-3)*(2*k+1));
    else /* same value, more efficient */
      t = gdivgu(t, ((k-3)*(2*k+1)) / 3);
    P[k<<1] = gerepileupto(av, t);
  }
  return res;
}

static int
get_c4c6(GEN w, GEN *c4, GEN *c6, long prec)
{
  if (typ(w) == t_VEC) switch(lg(w))
  {
    case 17:
      *c4 = ell_get_c4(w);
      *c6 = ell_get_c6(w);
      return 1;
    case 3:
    {
      ellred_t T;
      if (!get_periods(w,NULL,&T, prec)) break;
      *c4 = _elleisnum(&T, 4);
      *c6 = gneg(_elleisnum(&T, 6));
      return 1;
    }
  }
  *c4 = *c6 = NULL;
  return 0;
}

GEN
ellwpseries(GEN e, long v, long PRECDL)
{
  GEN c4, c6;
  checkell(e);
  c4 = ell_get_c4(e);
  c6 = ell_get_c6(e); return ellwpseries_aux(c4,c6,v,PRECDL);
}

GEN
ellwp(GEN w, GEN z, long prec)
{ return ellwp0(w,z,0,prec); }

GEN
ellwp0(GEN w, GEN z, long flag, long prec)
{
  pari_sp av = avma;
  GEN y;

  if (flag && flag != 1) pari_err_FLAG("ellwp");
  if (!z) z = pol_x(0);
  y = toser_i(z);
  if (y)
  {
    long vy = varn(y), v = valser(y);
    GEN P, Q, c4,c6;
    if (!get_c4c6(w,&c4,&c6,prec)) pari_err_TYPE("ellwp",w);
    if (v <= 0) pari_err(e_IMPL,"ellwp(t_SER) away from 0");
    if (gequal0(y)) {
      set_avma(av);
      if (!flag) return zeroser(vy, -2*v);
      retmkvec2(zeroser(vy, -2*v), zeroser(vy, -3*v));
    }
    P = ellwpseries_aux(c4,c6, vy, lg(y)-2);
    Q = gsubst(P, varn(P), y);
    if (!flag)
      return gerepileupto(av, Q);
    else
    {
      GEN R = mkvec2(Q, gdiv(derivser(Q), derivser(y)));
      return gerepilecopy(av, R);
    }
  }
  y = ellwpnum_all(w,z,flag,prec);
  if (!y) pari_err_DOMAIN("ellwp", "argument","=", gen_0,z);
  return gerepileupto(av, y);
}

GEN
ellzeta(GEN w, GEN z, long prec0)
{
  long prec;
  pari_sp av = avma;
  GEN pi2, q, y, et = NULL;
  ellred_t T;

  if (!z) z = pol_x(0);
  y = toser_i(z);
  if (y)
  {
    long vy = varn(y), v = valser(y);
    GEN P, Q, c4,c6;
    if (!get_c4c6(w,&c4,&c6,prec0)) pari_err_TYPE("ellzeta",w);
    if (v <= 0) pari_err(e_IMPL,"ellzeta(t_SER) away from 0");
    if (gequal0(y)) { set_avma(av); return zeroser(vy, -v); }
    P = ellwpseries_aux(c4,c6, vy, lg(y)-2);
    P = integser(gneg(P)); /* \zeta' = - \wp*/
    Q = gsubst(P, varn(P), y);
    return gerepileupto(av, Q);
  }
  if (!get_periods(w, z, &T, prec0)) pari_err_TYPE("ellzeta", w);
  if (!T.Z) pari_err_DOMAIN("ellzeta", "z", "=", gen_0, z);
  prec = T.prec;
  if (signe(T.x) || signe(T.y)) et = eta_correction(&T, _elleta(&T));

  pi2 = Pi2n(1, prec);
  q = expIPiC(gmul2n(T.Tau,1), prec);
  y = mulcxI(gmul(cxEk(T.Tau,2,prec), gmul(T.Z,divrs(pi2,-12))));
  if (!T.abs_u_is_1 || (!gequal(T.Z,ghalf) && !gequal(T.Z,gneg(ghalf))))
  { /* else u = -1 and this vanishes */
    long toadd = (long)ceil(get_toadd(T.Z));
    GEN qn, u, v, S = gen_0;
    pari_sp av1;
    u = expIPiC(gmul2n(T.Z,1), prec);
    v = gadd(ghalf, ginv(gsubgs(u, 1)));
    if (T.abs_u_is_1) gel(v,1) = gen_0; /*v = (u+1)/2(u-1), pure imaginary*/
    y = gadd(y, v);
    /* add sum_n q^n ( u/(u*q^n - 1) + 1/(u - q^n) )
     *     = (u^2 - 1) sum_n q^n / (uq^n - 1)(u - q^n) */
    av1 = avma;
    for (qn = q;;)
    {
      S = gadd(S, gdiv(qn, gmul(gsubgs(gmul(qn,u),1), gsub(u,qn))));
      qn = gmul(q,qn);
      if (gexpo(qn) <= - prec2nbits(prec) - 5 - toadd) break;
      if (gc_needed(av1,1))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"ellzeta");
        gerepileall(av1,2, &S,&qn);
      }
    }
    y = gadd(y, gmul(gsubgs(gsqr(u),1), S));
  }
  y = mulcxI(gmul(gdiv(pi2,T.W2), y));
  if (T.some_q_is_real)
  {
    if (T.some_z_is_real)
    {
      if (!et || typ(et) != t_COMPLEX) y = real_i(y);
    }
    else if (T.some_z_is_pure_imag)
    {
      if (!et || (typ(et) == t_COMPLEX && isintzero(gel(et,1))))
        gel(y,1) = gen_0;
    }
  }
  if (et) y = gadd(y, et);
  return gerepilecopy(av, gprec_wtrunc(y, T.prec0));
}

/* if flag=0, return ellsigma, otherwise return log(ellsigma) */
GEN
ellsigma(GEN w, GEN z, long flag, long prec0)
{
  long toadd, prec, n;
  pari_sp av = avma, av1;
  GEN u, urn, urninv, z0, pi, pi2, q, q8, qn2, qn, y, y1, uinv, et, etnew;
  ellred_t T;

  if (flag < 0 || flag > 1) pari_err_FLAG("ellsigma");

  if (!z) z = pol_x(0);
  y = toser_i(z);
  if (y)
  {
    long vy = varn(y), v = valser(y);
    GEN P, Q, c4,c6;
    if (!get_c4c6(w,&c4,&c6,prec0)) pari_err_TYPE("ellsigma",w);
    if (v <= 0) pari_err_IMPL("ellsigma(t_SER) away from 0");
    if (flag) pari_err_TYPE("log(ellsigma)",y);
    if (gequal0(y)) { set_avma(av); return zeroser(vy, -v); }
    P = ellwpseries_aux(c4,c6, vy, lg(y)-2);
    P = integser(gneg(P)); /* \zeta' = - \wp*/
    /* (log \sigma)' = \zeta; remove log-singularity first */
    P = integser(serchop0(P));
    P = gexp(P, prec0);
    setvalser(P, valser(P)+1);
    Q = gsubst(P, varn(P), y);
    return gerepileupto(av, Q);
  }
  if (!get_periods(w, z, &T, prec0)) pari_err_TYPE("ellsigma",w);
  if (!T.Z)
  {
    if (!flag) return gen_0;
    pari_err_DOMAIN("log(ellsigma)", "argument","=",gen_0,z);
  }
  prec = T.prec;
  pi2 = Pi2n(1,prec);
  pi  = mppi(prec);

  urninv = uinv = NULL;
  if (typ(T.Z) == t_FRAC && equaliu(gel(T.Z,2), 2) && equalim1(gel(T.Z,1)))
  {
    toadd = 0;
    urn = mkcomplex(gen_0, gen_m1); /* Z = -1/2 => urn = -I */
    u = gen_1;
  }
  else
  {
    toadd = (long)ceil(fabs( get_toadd(T.Z) ));
    urn = expIPiC(T.Z, prec); /* exp(i Pi Z) */
    u = gneg_i(gsqr(urn));
    if (!T.abs_u_is_1) { urninv = ginv(urn); uinv = gneg_i(gsqr(urninv)); }
  }
  q8 = expIPiC(gmul2n(T.Tau, -2), prec);
  q = gpowgs(q8,8); av1 = avma;
  y = gen_0; qn = q; qn2 = gen_1;
  for(n=0;;n++)
  { /* qn = q^(n+1), qn2 = q^(n(n+1)/2), urn = u^((n+1)/2)
     * if |u| = 1, will multiply by 2*I at the end ! */
    y = gadd(y, gmul(qn2, uinv? gsub(urn,urninv): imag_i(urn)));
    qn2 = gmul(qn,qn2);
    if (gexpo(qn2) + n*toadd <= - prec2nbits(prec) - 5) break;
    qn  = gmul(q,qn);
    urn = gmul(urn,u);
    if (uinv) urninv = gmul(urninv,uinv);
    if (gc_needed(av1,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ellsigma");
      gerepileall(av1,urninv? 5: 4, &y,&qn,&qn2,&urn,&urninv);
    }
  }
  y = gmul(y, gdiv(q8, gmul(pi2, gpowgs(trueeta(T.Tau,prec),3))));
  y = gmul(y, T.abs_u_is_1? gmul2n(T.W2,1): mulcxmI(T.W2));

  et = _elleta(&T);
  z0 = gmul(T.Z,T.W2);
  y1 = gadd(z0, gmul2n(gadd(gmul(T.x,T.W1), gmul(T.y,T.W2)),-1));
  etnew = gmul(eta_correction(&T, et), y1);
  y1 = gadd(etnew, gmul2n(gmul(gmul(T.Z,z0),gel(et,2)),-1));
  if (flag)
  {
    y = gadd(y1, glog(y,prec));
    if (mpodd(T.x) || mpodd(T.y)) y = gadd(y, mulcxI(pi));
    /* log(real number): im(y) = 0 or Pi */
    if (T.some_q_is_real && isintzero(imag_i(z)) && gexpo(imag_i(y)) < 1)
      y = real_i(y);
  }
  else
  {
    y = gmul(y, gexp(y1,prec));
    if (mpodd(T.x) || mpodd(T.y)) y = gneg_i(y);
    if (T.some_q_is_real)
    {
      int re, cx;
      check_complex(z,&re,&cx);
      if (re) y = real_i(y);
      else if (cx && typ(y) == t_COMPLEX) gel(y,1) = gen_0;
    }
  }
  return gerepilecopy(av, gprec_wtrunc(y, T.prec0));
}

GEN
pointell(GEN e, GEN z, long prec)
{
  pari_sp av = avma;
  GEN v;

  checkell(e);
  if (ell_get_type(e) == t_ELL_Qp)
  {
    prec = minss(ellQp_get_prec(e), padicprec_relative(z));
    return ellQp_t2P(e, z, prec);
  }
  v = ellwpnum_all(e,z,1,prec);
  if (!v) { set_avma(av); return ellinf(); }
  gel(v,1) = gsub(gel(v,1), gdivgu(ell_get_b2(e),12));
  gel(v,2) = gmul2n(gsub(gel(v,2), ec_h_evalx(e,gel(v,1))),-1);
  return gerepilecopy(av, v);
}

/********************************************************************/
/**                                                                **/
/**                 Tate's algorithm e (cf Anvers IV)              **/
/**               Kodaira types, global minimal model              **/
/**                                                                **/
/********************************************************************/
/* structure to hold incremental computation of standard minimal model/Q */
typedef struct {
  long a1; /*{0,1}*/
  long a2; /*{-1,0,1}*/
  long a3; /*{0,1}*/
  long b2; /* centermod(-c6, 12), in [-5,6] */
  GEN u, u2, u3, u4, u6;
  GEN a4, a6, b4, b6, b8, c4, c6, D;
} ellmin_t;

/* u from [u,r,s,t] */
static void
min_set_u(ellmin_t *M, GEN u)
{
  M->u = u;
  if (is_pm1(u))
    M->u2 = M->u3 = M->u4 = M->u6 = gen_1;
  else
  {
    M->u2 = sqri(u);
    M->u3 = mulii(M->u2, u);
    M->u4 = sqri(M->u2);
    M->u6 = sqri(M->u3);
  }
}
/* E = original curve */
static void
min_set_c(ellmin_t *M, GEN E)
{
  GEN c4 = ell_get_c4(E), c6 = ell_get_c6(E);
  if (!is_pm1(M->u4)) {
    c4 = diviiexact(c4, M->u4);
    c6 = diviiexact(c6, M->u6);
  }
  M->c4 = c4;
  M->c6 = c6;
}
static void
min_set_D(ellmin_t *M, GEN E)
{
  GEN D = ell_get_disc(E);
  if (!is_pm1(M->u6)) D = diviiexact(D, sqri(M->u6));
  M->D = D;
}
static void
min_set_b(ellmin_t *M)
{
  long b22, b2;
  M->b2 = b2 = Fl_center(12 - umodiu(M->c6,12), 12, 6);
  b22 = b2 * b2; /* in [0,36] */
  M->b4 = diviuexact(subui(b22, M->c4), 24);
  M->b6 = diviuexact(subii(mulsi(b2, subiu(mului(36,M->b4),b22)), M->c6), 216);
}
static void
min_set_a(ellmin_t *M)
{
  long a1, a2, a3, a13, b2 = M->b2;
  GEN b4 = M->b4, b6 = M->b6;
  if (odd(b2))
  {
    a1 = 1;
    a2 = (b2 - 1) >> 2;
  }
  else
  {
    a1 = 0;
    a2 = b2 >> 2;
  }
  M->a1 = a1;
  M->a2 = a2;
  M->a3 = a3 = Mod2(b6)? 1: 0;
  a13 = a1 & a3; /* a1 * a3 */
  M->a4 = shifti(subiu(b4, a13), -1);
  M->a6 = shifti(subiu(b6, a3), -2);
}
static void
min_set_all(ellmin_t *M, GEN E, GEN u)
{
  min_set_u(M, u);
  min_set_c(M, E);
  min_set_D(M, E);
  min_set_b(M);
  min_set_a(M);
}
static GEN
min_to_ell(ellmin_t *M, GEN E)
{
  GEN b8, y = obj_init(15, 8);
  long a11, a13;
  gel(y,1) = M->a1? gen_1: gen_0;
  gel(y,2) = stoi(M->a2);
  gel(y,3) = M->a3? gen_1: gen_0;
  gel(y,4) = M->a4;
  gel(y,5) = M->a6;
  gel(y,6) = stoi(M->b2);
  gel(y,7) = M->b4;
  gel(y,8) = M->b6;
  a11 = M->a1;
  a13 = M->a1 & M->a3;
  b8 = subii(addii(mului(a11,M->a6), mulis(M->b6, M->a2)),
             mulii(M->a4, addiu(M->a4,a13)));
  gel(y,9) = b8; /* a1^2 a6 + 4a6 a2 + a2 a3^2 - a4(a4 + a1 a3) */
  gel(y,10)= M->c4;
  gel(y,11)= M->c6;
  gel(y,12)= M->D;
  gel(y,13)= gel(E,13);
  gel(y,14)= gel(E,14);
  gel(y,15)= gel(E,15);
  return y;
}
static GEN
min_get_v(ellmin_t *M, GEN E)
{
  GEN r, s, t;
  r = diviuexact(subii(mulis(M->u2,M->b2), ell_get_b2(E)), 12);
  s = shifti(subii(M->a1? M->u: gen_0, ell_get_a1(E)), -1);
  t = shifti(subii(M->a3? M->u3: gen_0, Zec_h_evalx(E,r)), -1);
  return mkvec4(M->u,r,s,t);
}

/* return v_p(u), where [u,r,s,t] is the variable change to minimal model */
static long
get_vp_u_small(GEN E, ulong p, long *pv6, long *pvD)
{
  GEN c6 = ell_get_c6(E);
  long d, v6, vD = Z_lval(ell_get_disc(E), p);
  if (!signe(c6))
  {
    d = vD / 12;
    if (d)
    {
      if (p == 2)
      {
        GEN c4 = ell_get_c4(E);
        long a = Mod16( shifti(c4, -4*d) );
        if (a) d--;
      }
      if (d) vD -= 12*d; /* non minimal model */
    }
    v6 = 12; /* +oo */
  }
  else
  {
    v6 = Z_lval(c6,p);
    d = minss(2*v6, vD) / 12;
    if (d) {
      if (p == 2) {
        GEN c4 = ell_get_c4(E);
        long a = Mod16( shifti(c4, -4*d) );
        long b = Mod32( shifti(c6, -6*d) );
        if ((b & 3) != 3 && (a || (b && b!=8))) d--;
      } else if (p == 3) {
        if (v6 == 6*d+2) d--;
      }
      if (d) { v6 -= 6*d; vD -= 12*d; } /* non minimal model */
    }
  }
  *pv6 = v6; *pvD = vD; return d;
}
static long
get_vp_u(GEN E, GEN p, long *pv6, long *pvD)
{
  GEN c6;
  long d, v6, vD;
  if (lgefint(p) == 3) return get_vp_u_small(E, p[2], pv6, pvD);
  c6 = ell_get_c6(E);
  vD = Z_pval(ell_get_disc(E), p);
  if (!signe(c6))
  {
    d = vD / 12;
    if (d) vD -= 12*d; /* non minimal model */
    v6 = 12; /* +oo */
  }
  else
  {
    v6 = Z_pval(c6,p);
    d = minss(2*v6, vD) / 12;
    if (d) { v6 -= 6*d; vD -= 12*d; } /* non minimal model */
  }
  *pv6 = v6; *pvD = vD; return d;
}

/* Given an integral elliptic curve in ellinit form, and a prime p, returns the
  type of the fiber at p of the Neron model, as well as the change of variables
  in the form [f, kod, v, c].

  * The integer f is the conductor's exponent.

  * The integer kod is the Kodaira type using the following notation:
    II , III , IV  -->  2, 3, 4
    I0  -->  1
    Inu --> 4+nu for nu > 0
  A '*' negates the code (e.g I* --> -2)

  * v is a quadruple [u, r, s, t] yielding a minimal model

  * c is the Tamagawa number.

  Uses Tate's algorithm (Anvers IV). Given the remarks at the bottom of
  page 46, the "long" algorithm is used for p = 2,3 only. */
static GEN
localred_result(long f, long kod, long c, GEN v)
{
  GEN z = cgetg(5, t_VEC);
  gel(z,1) = stoi(f);
  gel(z,2) = stoi(kod);
  gel(z,3) = gcopy(v);
  gel(z,4) = stoi(c); return z;
}
static GEN
localredbug(GEN p, const char *s)
{
  if (BPSW_psp(p)) pari_err_BUG(s);
  pari_err_PRIME("localred",p);
  return NULL; /* LCOV_EXCL_LINE */
}

/* v_p( denom(j(E)) ) >= 0 */
static long
j_pval(GEN E, GEN p) { return Z_pval(Q_denom(ell_get_j(E)), p); }

/* p > 3, e integral */
static GEN
localred_p(GEN e, GEN p)
{
  long k, f, kod, c, nuj, nuD, nu6;
  GEN p2, v, tri, c4, c6, D = ell_get_disc(e);

  c4 = ell_get_c4(e);
  c6 = ell_get_c6(e);
  nuj = j_pval(e, p);
  nuD = Z_pval(D, p);
  k = get_vp_u(e, p, &nu6, &nuD);
  if (!k) v = init_ch();
  else
  { /* model not minimal */
    ellmin_t M;
    min_set_all(&M, e, powiu(p,k));
    v = min_get_v(&M, e);
    c4 = M.c4; c6 = M.c6; D = M.D;
  }

  if (nuj > 0) switch(nuD - nuj)
  {
    case 0: f = 1; kod = 4+nuj; /* Inu */
      switch(kronecker(negi(c6),p))
      {
        case  1: c = nuD; break;
        case -1: c = odd(nuD)? 1: 2; break;
        default: return localredbug(p,"localred (p | c6)");
      }
      break;
    case 6:
    {
      GEN d = Fp_red(diviiexact(D, powiu(p, 6+nuj)), p);
      if (nuj & 1) d = Fp_mul(d, diviiexact(c6, powiu(p,3)), p);
      f = 2; kod = -4-nuj; c = 3 + kronecker(d, p); /* Inu* */
      break;
    }
    default: return localredbug(p,"localred (nu_D - nu_j != 0,6)");
  }
  else switch(nuD)
  {
    case  0: f = 0; kod = 1; c = 1; break; /* I0, regular */
    case  2: f = 2; kod = 2; c = 1; break; /* II   */
    case  3: f = 2; kod = 3; c = 2; break; /* III  */
    case  4: f = 2; kod = 4; /* IV   */
      c = 2 + krosi(-6,p) * kronecker(diviiexact(c6,sqri(p)), p);
      break;
    case  6: f = 2; kod = -1; /* I0*  */
      p2 = sqri(p);
      /* x^3 - 3c4/p^2 x - 2c6/p^3 */
      tri = mkpoln(4, gen_1, gen_0,
                            negi(mului(3, diviiexact(c4, p2))),
                            negi(shifti(diviiexact(c6, mulii(p2,p)), 1)));
      c = 1 + FpX_nbroots(tri, p);
      break;
    case  8: f = 2; kod = -4; /* IV*  */
      c = 2 + krosi(-6,p) * kronecker(diviiexact(c6, sqri(sqri(p))), p);
      break;
    case  9: f = 2; kod = -3; c = 2; break; /* III* */
    case 10: f = 2; kod = -2; c = 1; break; /* II*  */
    default: return localredbug(p,"localred");
  }
  return localred_result(f, kod, c, v);
}

/* return a_{ k,l } in Tate's notation, pl = p^l */
static ulong
aux(GEN ak, ulong q, ulong pl)
{ return umodiu(ak, q) / pl; }

static ulong
aux2(GEN ak, ulong p, GEN pl)
{ pari_sp av = avma; return gc_ulong(av, umodiu(diviiexact(ak, pl), p)); }

/* number of distinct roots of X^3 + aX^2 + bX + c modulo p = 2 or 3
 * assume a,b,c in {0, 1} [ p = 2 ] or {0, 1, 2} [ p = 3 ]
 * if there's a multiple root, put it in *mult */
static long
numroots3(long a, long b, long c, long p, long *mult)
{
  if (p == 2)
  {
    if (odd(c + a * b)) return 3;
    *mult = b; return odd(a + b)? 2: 1;
  }
  /* p = 3 */
  if (!a) { *mult = -c; return b? 3: 1; }
  *mult = a * b;
  if (b == 2)
    return (a + c) == 3 ? 2 : 3;
  else
    return c ? 3 : 2;
}

/* same for aX^2 +bX + c */
static long
numroots2(long a, long b, long c, long p, long *mult)
{
  if (p == 2) { *mult = c; return odd(b)? 2: 1; }
  /* p = 3 */
  *mult = a * b; return (b * b - a * c) % 3 ? 2 : 1;
}

/* p = 2 or 3 */
static GEN
localred_23(GEN e, long p)
{
  long c, nu, nu6, nuD, r, s, t;
  long k, theroot, p2, p3, p4, p5, a21, a42, a63, a32, a64;
  GEN v;

  k = get_vp_u_small(e, p, &nu6, &nuD);
  if (!k) v = init_ch();
  else
  {
    ellmin_t M;
    min_set_all(&M, e, powuu(p, k));
    v = min_get_v(&M, e);
    e = min_to_ell(&M, e);
  }
  /* model is minimal */
  nuD = Z_lval(ell_get_disc(e), (ulong)p);
  if (!nuD) return localred_result(0, 1, 1, v); /* I0 */
  if (p == 2) { p2 = 4; p3 = 8;  p4 = 16; p5 = 32; }
  else        { p2 = 9; p3 = 27; p4 = 81; p5 =243; }

  if (umodiu(ell_get_b2(e), p)) /* p \nmid b2 */
  {
    if (umodiu(negi(ell_get_c6(e)), p == 2 ? 8 : 3) == 1)
      c = nuD;
    else
      c = odd(nuD)? 1: 2;
    return localred_result(1, 4 + nuD, c, v); /* Inu */
  }
  if (p == 2)
  {
    r = umodiu(ell_get_a4(e), 2);
    s = umodiu(ell_get_a2(e), 2);
    t = umodiu(ell_get_a6(e), 2);
    if (r) { t = (s + t) & 1; s = (s + 1) & 1; }
  }
  else /* p == 3 */
  {
    r = - umodiu(ell_get_b6(e), 3);
    s = umodiu(ell_get_a1(e), 3);
    t = umodiu(ell_get_a3(e), 3);
    if (s) { t  = (t + r*s) % 3; if (t < 0) t += 3; }
  }
  /* p | (a1, a2, a3, a4, a6) */
  if (r || s || t) e = coordch_rst(e, stoi(r), stoi(s), stoi(t));
  if (umodiu(ell_get_a6(e), p2))
    return localred_result(nuD, 2, 1, v); /* II */
  if (umodiu(ell_get_b8(e), p3))
    return localred_result(nuD - 1, 3, 2, v); /* III */
  if (umodiu(ell_get_b6(e), p3))
  {
    if (umodiu(ell_get_b6(e), (p==2)? 32: 27) == (ulong)p2)
      c = 3;
    else
      c = 1;
    return localred_result(nuD - 2, 4, c, v); /* IV */
  }

  if (umodiu(ell_get_a6(e), p3))
    e = coordch_t(e, p == 2? gen_2: modis(ell_get_a3(e), 9));
      /* p | a1, a2; p^2  | a3, a4; p^3 | a6 */
  a21 = aux(ell_get_a2(e), p2, p);
  a42 = aux(ell_get_a4(e), p3, p2);
  a63 = aux(ell_get_a6(e), p4, p3);
  switch (numroots3(a21, a42, a63, p, &theroot))
  {
    case 3:
      c = a63 ? 1: 2;
      if (p == 2)
        c += ((a21 + a42 + a63) & 1);
      else {
        if (((1 + a21 + a42 + a63) % 3) == 0) c++;
        if (((1 - a21 + a42 - a63) % 3) == 0) c++;
      }
      return localred_result(nuD - 4, -1, c, v); /* I0* */
    case 2:
    { /* compute nu */
      GEN pk, pk1, p2k;
      long al, be, ga;
      if (theroot) e = coordch_r(e, stoi(theroot * p));
          /* p | a1; p^2  | a2, a3; p^3 | a4; p^4 | a6 */
      nu = 1;
      pk  = utoipos(p2);
      p2k = utoipos(p4);
      for(;;)
      {
        be =  aux2(ell_get_a3(e), p, pk);
        ga = -aux2(ell_get_a6(e), p, p2k);
        al = 1;
        if (numroots2(al, be, ga, p, &theroot) == 2) break;
        if (theroot) e = coordch_t(e, mulsi(theroot,pk));
        pk1 = pk;
        pk  = mului(p, pk);
        p2k = mului(p, p2k); nu++;

        al = a21;
        be = aux2(ell_get_a4(e), p, pk);
        ga = aux2(ell_get_a6(e), p, p2k);
        if (numroots2(al, be, ga, p, &theroot) == 2) break;
        if (theroot) e = coordch_r(e, mulsi(theroot, pk1));
        p2k = mului(p, p2k); nu++;
      }
      if (p == 2)
        c = odd(ga)? 2: 4;
      else
        c = 3 + kross(be * be - al * ga, 3);
      return localred_result(nuD - 4 - nu, -4 - nu, c, v); /* Inu* */
    }
    case 1:
      if (theroot) e = coordch_r(e, stoi(theroot*p));
          /* p | a1; p^2  | a2, a3; p^3 | a4; p^4 | a6 */
      a32 = aux(ell_get_a3(e), p3, p2);
      a64 = aux(ell_get_a6(e), p5, p4);
      if (numroots2(1, a32, -a64, p, &theroot) == 2)
      {
        if (p == 2)
          c = 3 - 2 * a64;
        else
          c = 2 + kross(a32 * a32 + a64, 3);
        return localred_result(nuD - 6, -4, c, v); /* IV* */
      }
      if (theroot) e = coordch_t(e, stoi(theroot*p2));
          /* p | a1; p^2 | a2; p^3 | a3, a4; p^5 | a6 */
      if (umodiu(ell_get_a4(e), p4))
        return localred_result(nuD - 7, -3, 2, v); /* III* */

      /* p^6 \nmid a6, otherwise wouldn't be minimal */
      return localred_result(nuD - 8, -2, 1, v); /* II* */
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/* e is integral */
static GEN
localred(GEN e, GEN p)
{
  if (abscmpiu(p, 3) > 0) /* p != 2,3 */
    return localred_p(e,p);
  else
  {
    long l = itos(p);
    if (l < 2) pari_err_PRIME("localred",p);
    return localred_23(e, l);
  }
}

/* Given J an ideal in HNF coprime to 2 and z algebraic integer,
 * return b algebraic integer such that z + 2b in  J */
static GEN
approx_mod2(GEN J, GEN z)
{
  GEN b = z;
  long i;
  if (typ(b) == t_INT)
  {
    if (mpodd(b)) b = addii(b, gcoeff(J,1,1));
    return shifti(negi(b),-1);
  }
  for (i = lg(J)-1; i >= 1; i--)
  {
    if (mpodd(gel(b,i))) b = ZC_add(b, gel(J,i));
  }
  return gshift(ZC_neg(b), -1);
}

/* Given J an ideal in HNF coprime to 3 and z algebraic integer,
 * return b algebraic integer such that z + 3b in  J */
static GEN
approx_mod3(GEN J, GEN z)
{
  GEN b = z;
  long i;
  if (typ(b) == t_INT)
  {
    long s = smodis(b,3);
    if (s)
    {
      GEN Jz = gcoeff(J,1,1);
      if (smodis(Jz, 3) == s)
        b = subii(b, Jz);
      else
        b = addii(b, Jz);
    }
    return diviiexact(b, stoi(-3));
  }
  for (i = lg(J)-1; i >= 1; i--)
  {
    long s = smodis(gel(b,i), 3);
    if (!s) continue;
    if (smodis(gcoeff(J,i,i), 3) == s)
      b = ZC_sub(b, gel(J,i));
    else
      b = ZC_add(b, gel(J,i));
  }
  return ZC_Z_divexact(b, stoi(-3));
}

/* return a such that v_P(a) = -1, integral elsewhere */
static GEN
get_piinv(GEN P)
{
  GEN z = pr_get_tau(P);
  if (typ(z) == t_MAT) z = gel(z,1);
  return gdiv(z, pr_get_p(P));
}
/* pi = local uniformizer, pv = 1/pi */
static void
get_uniformizers(GEN nf, GEN P, GEN *pi, GEN *pv)
{
  if (pr_is_inert(P))
  {
    *pi = pr_get_p(P);
    *pv = mkfrac(gen_1, *pi);
  }
  else
  {
    *pv = get_piinv(P);
    *pi = nfinv(nf, *pv);
  }
}
/* x^2+E.a1*x-E.a2 */
static GEN
pola1a2(GEN e, GEN nf, GEN modP)
{
  GEN a1 = nf_to_Fq(nf, ell_get_a1(e), modP);
  GEN a2 = nf_to_Fq(nf, ell_get_a2(e), modP);
  return mkpoln(3, gen_1, a1, gneg(a2));
}

/* x^2+E.a3*pv3*x-E.a6*pv6 */
static GEN
pola3a6(GEN e, GEN nf, GEN modP, GEN pv3, GEN pv6)
{
  GEN a3 = nf_to_Fq(nf, nfmul(nf, ell_get_a3(e), pv3), modP);
  GEN a6 = nf_to_Fq(nf, nfmul(nf, ell_get_a6(e), pv6), modP);
  return mkpoln(3, gen_1, a3, gneg(a6));
}

/* E.a2*pv2*x^2 + E.a4*pv4*x + E.a6*pv6 */

static GEN
pola2a4a6(GEN e, GEN nf, GEN modP, GEN pv2, GEN pv4, GEN pv6)
{
  GEN a2 = nf_to_Fq(nf, nfmul(nf, ell_get_a2(e), pv2), modP);
  GEN a4 = nf_to_Fq(nf, nfmul(nf, ell_get_a4(e), pv4), modP);
  GEN a6 = nf_to_Fq(nf, nfmul(nf, ell_get_a6(e), pv6), modP);
  return mkpoln(3, a2, a4, a6);
}

static GEN
pol2sqrt_23(GEN modP, GEN Q)
{
  GEN p = modpr_get_p(modP), T = modpr_get_T(modP);
  GEN r = absequaliu(p,2) ? gel(Q,2): gel(Q,3);
  if (!gequal1(gel(Q,4))) r = Fq_div(r, gel(Q,4), T, p);
  if (absequaliu(p,2)) r = Fq_sqrt(r,T,p);
  return Fq_to_nf(r, modP);
}

static GEN
nflocalred_section7(GEN e, GEN nf, GEN modP, GEN pi, GEN pv, long vD, GEN ch)
{
  GEN p = modpr_get_p(modP), T = modpr_get_T(modP);
  GEN pi3 = nfsqr(nf,pi);
  GEN pv3 = nfsqr(nf,pv), pv4 = nfmul(nf,pv,pv3), pv6 = nfsqr(nf,pv3);
  long n = 1;
  while(1)
  {
    GEN Q = pola3a6(e, nf, modP, pv3, pv6);
    GEN gama;
    if (FqX_is_squarefree(Q, T, p))
    {
      long nr = FqX_nbroots(Q,T,p);
      return localred_result(vD-n-4,-4-n,nr+2,ch);
    }
    gama = pol2sqrt_23(modP, Q);
    nf_compose_t(nf, &ch, &e, nfmul(nf, gama,pi3));
    pv6 = nfmul(nf,pv,pv6); n++;
    Q = pola2a4a6(e, nf, modP, pv, pv4, pv6);
    if (FqX_is_squarefree(Q, T, p))
    {
      long nr = FqX_nbroots(Q,T,p);
      return localred_result(vD-n-4,-4-n,nr+2,ch);
    }
    gama = pol2sqrt_23(modP, Q);
    nf_compose_r(nf, &ch, &e, nfmul(nf, gama, pi3));
    pi3 = nfmul(nf,pi, pi3);
    pv3 = pv4; pv4 = nfmul(nf,pv,pv4); pv6 = nfmul(nf,pv,pv6); n++;
  }
}

/* Tate algorithm, following J.H. Silverman GTM 151, chapt. IV, algo 9.4 */
/* Dedicated to John Tate for his kind words */

static GEN
nflocalred_23(GEN nf, GEN e, GEN D, GEN P, long *ap)
{
  GEN T, p, modP;
  long vD;
  GEN ch, pv, pv2, pv4, pi, pol;
  modP = nf_to_Fq_init(nf,&P,&T,&p);
  get_uniformizers(nf,P, &pi, &pv);
  ch = init_ch();
  vD = nfval(nf,D,P);
  *ap = 0;
  while(1)
  {
    if (vD==0)
      return localred_result(0,1,1,ch);
    else
    {
      GEN a1 = nf_to_Fq(nf, ell_get_a1(e), modP);
      GEN a2 = nf_to_Fq(nf, ell_get_a2(e), modP);
      GEN a3 = nf_to_Fq(nf, ell_get_a3(e), modP);
      GEN a4 = nf_to_Fq(nf, ell_get_a4(e), modP);
      GEN a6 = nf_to_Fq(nf, ell_get_a6(e), modP);
      GEN x0, y0;
      if (absequaliu(p,2))
      {
        GEN x02, y02;
        if (signe(a1))
        {
          x0 = Fq_div(a3, a1, T, p);
          x02 = Fq_sqr(x0,T,p);
          y02 = Fq_add(Fq_mul(x02,Fq_add(x0,a2,T,p),T,p),Fq_add(Fq_mul(a4,x0,T,p),a6,T,p),T,p);
        }
        else
        {
          x0 = Fq_sqrt(a4, T, p);
          y02 = Fq_add(Fq_mul(a4,a2,T,p),a6,T,p);
        }
        y0 = Fq_sqrt(y02,T,p);
      }
      else
      {
        GEN a12 = Fq_add(Fq_sqr(a1,T,p),a2,T,p);
        if (signe(a12))
          x0 = Fq_div(Fq_sub(a4,Fq_mul(a3,a1,T,p),T,p),a12,T,p);
        else
          x0 = Fq_sqrtn(Fq_neg(Fq_add(Fq_sqr(a3,T,p),a6,T,p),T,p),p,T,p,NULL);
        y0 = Fq_add(Fq_mul(a1, x0, T, p), a3, T, p);
      }
      x0 = Fq_to_nf(x0, modP);
      y0 = Fq_to_nf(y0, modP);
      nf_compose_rt(nf, &ch, &e, x0, y0);
    }
    /* 2 */
    {
      GEN b2 = nf_to_Fq(nf, ell_get_b2(e), modP);
      if (signe(b2) != 0)
      {
        GEN Q = pola1a2(e, nf, modP);
        long nr = FqX_nbroots(Q, T, p);
        if (nr==2) { *ap =  1; return localred_result(1,vD+4,vD,ch); /* Inu */ }
        else       { *ap = -1; return localred_result(1,vD+4,odd(vD)?1:2,ch);  }
      }
    }
    /* 3 */
    {
      long va6 = nfval(nf,ell_get_a6(e),P);
      if (va6 <= 1) return localred_result(vD,2,1,ch); /* II */
    }
    /* 4 */
    {
      long vb8 = nfval(nf,ell_get_b8(e),P);
      if (vb8 <= 2) return localred_result(vD-1,3,2,ch);/* III */
    }
    /* 5 */
    pv2 = nfsqr(nf,pv);
    {
      long vb6 = nfval(nf,ell_get_b6(e),P);
      if (vb6<=2)
      {
        GEN Q = pola3a6(e, nf, modP, pv, pv2);
        long nr = FqX_nbroots(Q,T,p);
        return localred_result(vD-2,4,1+nr,ch);/* IV */
      }
    }
    /* 6 */
    {
      GEN pv3 = nfmul(nf,pv, pv2);
      GEN alpha = pol2sqrt_23(modP, pola1a2(e, nf, modP));
      GEN beta  = pol2sqrt_23(modP, pola3a6(e, nf, modP, pv, pv2));
      GEN po2, E, F, mr;
      long i, lE;
      nf_compose_st(nf, &ch, &e, alpha, nfmul(nf, beta, pi));
      po2 = pola2a4a6(e, nf, modP, pv, pv2, pv3);
      if (signe(po2)) /* po2 = 0 is frequent when nonminimal */
      {
        pol = RgX_add(pol_xn(3,0), po2);
        F = FqX_factor(pol, T, p); E = gel(F,2);
        lE = lg(E);
        if (E[1] == 1 && (lE == 2 || E[2] == 1))
        { /* T squarefree, degree pattern is (3), (12) or (111) */
          long c; /* 1 + number of roots */
          switch(lE)
          {
            case 2: c = 1; break;
            case 3: c = 2; break;
            default: c = 4; break;
          }
          return localred_result(vD-4,-1,c,ch);/* I0* */
        }
      /* 7 */
        i = (lE == 2 || E[1] == 2)? 1: 2; /* index of multiple root */
        mr = constant_coeff(gmael(F,1,i)); /* - multiple root */
        if (!gequal0(mr))
        { /* not so frequent */
          GEN gama = Fq_to_nf(Fq_neg(mr, T, p), modP);
          nf_compose_r(nf, &ch, &e, nfmul(nf, gama,pi));
        }
        if (lE == 3)
          return nflocalred_section7(e, nf, modP, pi, pv, vD, ch); /* Inu* */
      }
    }
    pv4 = nfsqr(nf,pv2);
    pol = pola3a6(e, nf, modP, pv2, pv4);
    /*  8 */
    if (FqX_is_squarefree(pol,T,p))
    {
      long nr = FqX_nbroots(pol, T, p);
      return localred_result(vD-6,-4,1+nr,ch); /* IV* */
    }
    /*  9 */
    {
      GEN alpha = pol2sqrt_23(modP, pol);
      nf_compose_t(nf, &ch, &e, nfmul(nf, alpha, nfsqr(nf,pi)));
      if (nfval(nf, ell_get_a4(e), P) == 3)
        return localred_result(vD-7,-3,2,ch); /* III* */
    }
    /* 10 */
    if (nfval(nf, ell_get_a6(e), P) == 5)
      return localred_result(vD-8,-2,1,ch); /* II* */
    /* 11 */
    nf_compose_u(nf, &ch, &e, pi, pv);
    vD -= 12;
  }
}

/* Local reduction, residual characteristic >= 5. E/nf, P prid
* Output: f, kod, [u,r,s,t], c */
static GEN
nflocalred_p(GEN e, GEN P)
{
  GEN nf = ellnf_get_nf(e), T,p, modP = nf_to_Fq_init(nf,&P,&T,&p);
  long c, f, vD, nuj, kod, m;
  GEN ch, c4, c6, D, z, pi, piinv;

  c4 = ell_get_c4(e);
  c6 = ell_get_c6(e);
  D = ell_get_disc(e);
  vD = nfval(nf,D,P);
  nuj = nfval(nf,ell_get_j(e),P);
  nuj = nuj >= 0? 0: -nuj; /* v_P(denom(j)) */
  m = (vD - nuj)/12;
  get_uniformizers(nf,P, &pi, &piinv);

  if(m <= 0) ch = init_ch();
  else
  { /* model not minimal */
    GEN r,s,t, a1,a2,a3, u,ui,ui2,ui4,ui6,ui12;
    u = nfpow_u(nf,pi,m);
    ui = nfpow_u(nf,piinv,m);
    ui2 = nfsqr(nf,ui);
    ui4 = nfsqr(nf,ui2);
    ui6 = nfmul(nf,ui2,ui4);
    ui12 = nfsqr(nf,ui6);
    c4 = nfmul(nf,c4,ui4);
    c6 = nfmul(nf,c6,ui6);
    D = nfmul(nf,D,ui12);  vD -= 12*m;
    a1 = nf_to_scalar_or_basis(nf, ell_get_a1(e));
    a2 = nf_to_scalar_or_basis(nf, ell_get_a2(e));
    a3 = nf_to_scalar_or_basis(nf, ell_get_a3(e));
    s = approx_mod2(idealpow(nf,P,stoi(m)),   a1);
    r = gsub(a2, nfmul(nf,s,gadd(a1,s)));
    r = approx_mod3(idealpow(nf,P,stoi(2*m)), r);
    t = gadd(a3, nfmul(nf,r,a1));
    t = approx_mod2(idealpow(nf,P,stoi(3*m)), t);
    ch = mkvec4(u,r,s,t);
  }

  kod = 1; c = 1;
  /* minimal at P */
  if (nuj > 0)
  { /* v(j) < 0 */
    if (vD == nuj)
    { /* v(c4) = v(c6) = 0, multiplicative reduction */
      f = 1; kod = 4+vD;
      z = Fq_neg(nf_to_Fq(nf,c6,modP), T,p);
      if (Fq_issquare(z,T,p))
        c = vD;/* split */
      else
        c = odd(vD)?1 : 2; /* nonsplit */
    }
    else
    { /* v(c4) = 2, v(c6) = 3, potentially multiplicative */
      GEN Du;
      f = 2; kod = 2-vD;
      (void)nfvalrem(nf, D, P, &Du);
      z = nf_to_Fq(nf, Du, modP);
      if(odd(vD))
      {
        GEN c6u;
        (void)nfvalrem(nf, c6, P, &c6u);
        c6u = nf_to_Fq(nf, c6u, modP);
        z = Fq_mul(z, c6u, T,p);
      }
      c = Fq_issquare(z,T,p)? 4: 2;
    }
  }
  else
  { /* v(j) >= 0 */
    f = vD? 2: 0;
    switch(vD)
    {
      GEN piinv2, piinv3, piinv4, w;
      case 0: kod = 1; c = 1; break;
      case 2: kod = 2; c = 1; break;
      case 3: kod = 3; c = 2; break;
      case 4: kod = 4;
        z = nfmul(nf,c6,nfsqr(nf,piinv));
        z = nf_to_Fq(nf, z, modP);
        z = Fq_Fp_mul(z,stoi(-6),T,p);
        c = Fq_issquare(z,T,p)? 3: 1;
        break;
      case 6: kod = -1;
        piinv2 = nfsqr(nf,piinv);
        piinv3 = nfmul(nf,piinv,piinv2);
        z = nfmul(nf,c4,piinv2); z = nf_to_Fq(nf, z, modP);
        z = Fq_Fp_mul(z,stoi(-3), T,p);
        w = nfmul(nf,c6,piinv3); w = nf_to_Fq(nf, w, modP);
        w = Fq_Fp_mul(w,gen_m2, T,p);
        c = 1 + FqX_nbroots(mkpoln(4, gen_1,gen_0,z,w), T,p);
        break;
      case 8: kod = -4;
        piinv4 = nfpow_u(nf,piinv,4);
        z = nfmul(nf,c6,piinv4); z = nf_to_Fq(nf, z, modP);
        z = Fq_Fp_mul(z,stoi(-6),T,p);
        c = Fq_issquare(z,T,p)? 3: 1;
        break;
      case 9: kod = -3; c = 2; break;
      case 10: kod = -2; c = 1; break;
    }
  }
  return localred_result(f,kod,c,ch);
}
/* E is integral */
static GEN
nflocalred(GEN E, GEN pr)
{
  GEN p = pr_get_p(pr);
  if (abscmpiu(p, 3) <= 0)
  {
    long i, ap, vu;
    GEN nf = ellnf_get_nf(E), e = ell_to_nfell10(E), D = ell_get_disc(E);
    GEN q = nflocalred_23(nf,e,D,pr,&ap), v = gel(q,3), u = gel(v,1);
    gel(q,3) = v;
    /* do nothing if already minimal or equation was not pr-integral */
    vu = nfval(nf, u, pr);
    if (vu > 0)
    { /* remove denominators in r,s,t on nf.zk */
      GEN D, r = gel(v,2), s = gel(v,3), t = gel(v,4);
      D = Q_denom(mkvec3(r, s, t));
      if (!equali1(D))
      { /* Beware: D may not be coprime to pr */
        GEN a;
        (void)nfvalrem(nf, D, pr, &D);
        /* a in D/p^oo, = 1 mod (u^6) locally */
        a = idealaddtoone_i(nf, D, idealpows(nf, pr, 6*vu));
        gel(v,2) = nfmul(nf, r, a);
        gel(v,3) = nfmul(nf, s, a);
        gel(v,4) = nfmul(nf, t, a);
      }
    }
    for(i=1; i <= 4; i++) gel(v,i) = nftoalg(nf, gel(v,i));
    return q;
  }
  return nflocalred_p(E,pr);
}

static GEN
checkellp(GEN *pE, GEN p, GEN *pv, const char *s)
{
  GEN q, E = *pE;
  long tE;
  checkell(E); tE = ell_get_type(E);
  if (pv) *pv = NULL;
  if (p) switch(typ(p))
  {
    case t_INT:
      if (cmpis(p, 2) < 0) pari_err_DOMAIN(s,"p", "<", gen_2, p);
      break;
    case t_VEC:
      q = get_prid(p);
      if (q && tE == t_ELL_NF)
      {
        *pE = ellintegralmodel_i(E, pv);
        return q;
      }
    default: pari_err_TYPE(s,p);
  }
  switch(tE)
  {
    case t_ELL_Fp:
    case t_ELL_Fq: q = ellff_get_p(E); break;
    case t_ELL_Qp: q = ellQp_get_p(E); break;
    case t_ELL_Q: if (p) { q = p; p = NULL; break; }
    default:
      pari_err_TYPE(stack_strcat(s," [can't determine p]"), E);
      return NULL;/*LCOV_EXCL_LINE*/
  }
  if (p && !equalii(p, q)) pari_err_MODULUS(s, p,q);
  if (tE == t_ELL_Q || tE == t_ELL_Qp || tE == t_ELL_NF)
    *pE = ellintegralmodel_i(E, pv);
  return q;
}

GEN
elllocalred(GEN E, GEN p)
{
  pari_sp av = avma;
  GEN v, q;
  checkell(E);
  p = checkellp(&E, p, &v, "elllocalred");
  switch(ell_get_type(E))
  {
    case t_ELL_Qp:
    case t_ELL_Q:  q = localred(E, p); break;
    case t_ELL_NF: q = nflocalred(E, p); break;
    default: pari_err_TYPE("elllocalred", E);
      return NULL;/*LCOV_EXCL_LINE*/
  }
  if (v)
  { /* compose local change of variables with v */
    GEN u = gel(v,1), w = gel(q,3);
    if (is_trivial_change(w))
      gel(q,3) = mkvec4(u,gen_0,gen_0,gen_0);
    else
      gel(w,1) = gmul(u, gel(w,1));
  }
  return gerepilecopy(av, q);
}

/* typ(c) = t_INT or t_FRAC */
static GEN
handle_Q(GEN c, GEN *pd)
{
  *pd = (typ(c) == t_INT)? NULL: gel(c,2);
  return c;
}
static GEN
handle_coeff(GEN nf, GEN c, GEN *pd)
{
  *pd = NULL;
  switch(typ(c))
  {
    case t_INT: *pd = NULL; return c;
    case t_FRAC: *pd = gel(c,2); return c;
    case t_POL: case t_POLMOD: case t_COL:
      if (nf)
      {
        c = nf_to_scalar_or_basis(nf,c);
        return handle_Q(Q_content(c), pd);
      }
    default: pari_err_TYPE("ellintegralmodel",c);
      return NULL;/*LCOV_EXCL_LINE*/
  }
}
/* Return an integral model for e / nf, Q. Set v = NULL (already integral)
 * or the variable change [u,0,0,0], u = 1/t, t > 1 integer making e integral */
GEN
ellintegralmodel_i(GEN e, GEN *pv)
{
  GEN a, t, u, L, nf;
  long i, l, k;

  if (pv) *pv = NULL;
  /* t_ELL_Qp is also possible */
  nf = (ell_get_type(e) == t_ELL_NF)?ellnf_get_nf(e): NULL;
  L = cgetg(1, t_VEC); a = cgetg(6, t_VEC);
  for (i = 1; i < 6; i++)
  {
    GEN d;
    gel(a,i) = handle_coeff(nf, gel(e,i), &d);
    if (d) /* partial factorization of denominator */
      L = shallowconcat(L, gel(Z_factor_limit(d, 0),1));
  }
  /* a = [a1, a2, a3, a4, a6] */
  l = lg(L); if (l == 1) return e;
  L = ZV_sort_uniq_shallow(L);
  l = lg(L);

  t = gen_1;
  for (k = 1; k < l; k++)
  {
    GEN p = gel(L,k);
    long n = 0, m;
    for (i = 1; i < 6; i++)
      if (!gequal0(gel(a,i)))
      {
        long r = (i == 5)? 6: i; /* a5 is missing */
        m = r * n + Q_pval(gel(a,i), p);
        while (m < 0) { n++; m += r; }
      }
    t = mulii(t, powiu(p, n));
  }
  u = ginv(t);
  if (pv) *pv = mkvec4(u,gen_0,gen_0,gen_0);
  return coordch_uinv(e, t);
}
GEN
ellintegralmodel(GEN e, GEN *pv)
{
  pari_sp av = avma;
  checkell(e);
  switch(ell_get_type(e))
  {
    case t_ELL_Q:
    case t_ELL_Qp:
    case t_ELL_NF: break;
    default: pari_err_TYPE("ellintegralmodel",e);
  }
  e = ellintegralmodel_i(e, pv);
  if (pv && *pv) return gc_all(av, 2, &e, pv);
  e = gerepilecopy(av, e);
  if (pv) *pv = init_ch();
  return e;
}

/* return an integral model with a1 = a3 = 0 */
GEN
ellintegralbmodel(GEN e, GEN *pv)
{
  pari_sp av = avma;
  GEN f, a1, a3;

  checkell(e); f = ellintegralmodel_i(e, pv);
  a1 = ell_get_a1(f);
  a3 = ell_get_a3(f);
  if (!signe(a1) && !signe(a3))
  { if (!*pv) *pv = init_ch(); }
  else
  {
    GEN v = mkvec4(mpodd(a1) || mpodd(a3) ? ghalf : gen_1,
                   gen_0, gdivgs(a1,-2), gdivgs(a3,-2));
    gcomposev(pv, v); f = coordch(f, v);
  }
  if (f != e) ell_reset(f);
  return gc_all(av, 2, &f, pv);
}

static long
F2_card(ulong a1, ulong a2, ulong a3, ulong a4, ulong a6)
{
  long N = 1; /* oo */
  if (!a3) N ++; /* x = 0, y=0 or 1 */
  else if (!a6) N += 2; /* x = 0, y arbitrary */
  if ((a3 ^ a1) == 0) N++; /* x = 1, y = 0 or 1 */
  else if (a2 ^ a4 ^ a6) N += 2; /* x = 1, y arbitrary */
  return N;
}
static long
F3_card(ulong b2, ulong b4, ulong b6)
{
  ulong Po = 1+2*b4, Pe = b2+b6;
  /* kro(x,3)+1 = (x+1)%3, N = 4 + sum(kro) = 1+ sum(1+kro) */
  return 1+(b6+1)%3+(Po+Pe+1)%3+(2*Po+Pe+1)%3;
}
static long
cardmod2(GEN e)
{ /* solve y(1 + a1x + a3) = x (1 + a2 + a4) + a6 */
  ulong a1 = Rg_to_F2(ell_get_a1(e));
  ulong a2 = Rg_to_F2(ell_get_a2(e));
  ulong a3 = Rg_to_F2(ell_get_a3(e));
  ulong a4 = Rg_to_F2(ell_get_a4(e));
  ulong a6 = Rg_to_F2(ell_get_a6(e));
  return F2_card(a1,a2,a3,a4,a6);
}
static long
cardmod3(GEN e)
{
  ulong b2 = Rg_to_Fl(ell_get_b2(e), 3);
  ulong b4 = Rg_to_Fl(ell_get_b4(e), 3);
  ulong b6 = Rg_to_Fl(ell_get_b6(e), 3);
  return F3_card(b2,b4,b6);
}

static ulong
ZtoF2(GEN x) { return (ulong)mpodd(x); }

/* complete local reduction at 2, u = 2^d */
static void
min_set_2(ellmin_t *M, GEN E, long d)
{
  min_set_u(M, int2n(d));
  min_set_c(M, E);
  min_set_b(M);
  min_set_a(M);
}
/* local reduction at 3, u = 3^d, don't compute the a_i */
static void
min_set_3(ellmin_t *M, GEN E, long d)
{
  min_set_u(M, powuu(3, d));
  min_set_c(M, E);
  min_set_b(M);
}

static long
ellQap_u(GEN E, ulong p, int *good_red)
{
  long vc6, vD, d = get_vp_u_small(E, p, &vc6, &vD);
  if (vD) /* bad reduction */
  {
    GEN c6;
    long s;
    *good_red = 0;
    if (vc6) return 0;
    c6 = ell_get_c6(E);
    if (d) c6 = diviiexact(c6, powuu(p, 6*d));
    s = kroiu(c6,p);
    if ((p & 3) == 3) s = -s;
    return s;
  }
  *good_red = 1;
  if (p == 2)
  {
    ellmin_t M;
    if (!d) return 3 - cardmod2(E);
    min_set_2(&M, E, d);
    return 3 - F2_card(M.a1, M.a2 & 1, M.a3, ZtoF2(M.a4), ZtoF2(M.a6));
  }
  else if (p == 3)
  {
    ellmin_t M;
    if (!d) return 4 - cardmod3(E);
    min_set_3(&M, E, d);
    return 4 - F3_card(M.b2, umodiu(M.b4,3), umodiu(M.b6,3));
  }
  else
  {
    ellmin_t M;
    GEN a4, a6, pp = utoipos(p);
    min_set_u(&M, powuu(p,d));
    min_set_c(&M, E);
    c4c6_to_a4a6(M.c4, M.c6, pp, &a4,&a6);
    return itos( subui(p+1, Fp_ellcard(a4, a6, pp)) );
  }
}

static GEN
ellQap(GEN E, GEN p, int *good_red)
{
  GEN a4,a6, c4, c6, D;
  long vc6, vD, d;
  if (lgefint(p) == 3) return stoi( ellQap_u(E, p[2], good_red) );
  c6 = ell_get_c6(E);
  D = ell_get_disc(E);
  vc6 = Z_pval(c6,p); vD = Z_pval(D,p);
  d = minss(2*vc6, vD) / 12;
  if (d) { vc6 -= 6*d; vD -= 12*d; } /* non minimal model */
  if (vD) /* bad reduction */
  {
    long s;
    *good_red = 0;
    if (vc6) return gen_0;
    if (d) c6 = diviiexact(c6, powiu(p, 6*d));
    s = kronecker(c6,p);
    if (mod4(p) == 3) s = -s;
    return s < 0? gen_m1: gen_1;
  }
  *good_red = 1;
  c4 = ell_get_c4(E);
  if (d)
  {
    GEN u2 = powiu(p, 2*d), u4 = sqri(u2), u6 = mulii(u2,u4);
    c4 = diviiexact(c4, u4);
    c6 = diviiexact(c6, u6);
  }
  c4c6_to_a4a6(c4, c6, p, &a4,&a6);
  return subii(addiu(p,1), Fp_ellcard(a4, a6, p));
}

static GEN
doellcard(GEN E)
{
  GEN fg = ellff_get_field(E);
  if (typ(fg)==t_FFELT)
    return FF_ellcard(E);
  else
  {
    GEN e = ellff_get_a4a6(E);
    return Fp_ellcard(gel(e,1),gel(e,2),fg);
  }
}

static GEN
ellnfap(GEN E, GEN P, int *good_red)
{
  GEN a4,a6, card, nf = ellnf_get_nf(E);
  GEN T,p, modP = nf_to_Fq_init(nf,&P,&T,&p);
  if (abscmpiu(p, 3) <= 0)
  {
    long ap;
    GEN nf = ellnf_get_nf(E), e = ell_to_nfell10(E), D = ell_get_disc(E);
    GEN L = nflocalred_23(nf, e,D,P,&ap), kod = gel(L,2);
    if (!equali1(kod)) { *good_red = 0; return stoi(ap); }
    *good_red = 1;
    E = nf_coordch(nf, vecslice(e,1,5), gel(L,3));
    E = ellinit_nf_to_Fq(nf, E, modP);
    card = FF_ellcard(E);
  }
  else
  {
    GEN c6 = ell_get_c6(E), c4 = ell_get_c4(E);
    long vD = nfval(nf, ell_get_disc(E), P);
    if (vD)
    {
      GEN c6new;
      long d, vc6 = nfvalrem(nf,c6,P, &c6new);
      d = ((vc6 == LONG_MAX)? vD: minss(vD,2*vc6)) / 12;
      if (vD > 12*d)
      { /* bad reduction */
        *good_red = 0;
        if (vc6 != 6*d) return gen_0;
        c6 = nf_to_Fq(nf, c6new, modP);
        return Fq_issquare(gneg(c6),T,p)? gen_1: gen_m1;
      }
      if (d)
      { /* model not minimal at P */
        GEN piinv = get_piinv(P);
        GEN ui2 = nfpow(nf, piinv, stoi(2*d));
        GEN ui4 = nfsqr(nf, ui2);
        GEN ui6 = nfmul(nf, ui2, ui4);
        c4 = nfmul(nf, c4, ui4);
        c6 = nfmul(nf, c6, ui6);
      }
    }
    *good_red = 1;
    c4 = nf_to_Fq(nf, c4, modP);
    c6 = nf_to_Fq(nf, c6, modP);
    Fq_c4c6_to_a4a6(c4, c6, T,p, &a4,&a6);
    card = T? FpXQ_ellcard(Fq_to_FpXQ(a4,T,p),Fq_to_FpXQ(a6,T,p),T,p)
            : Fp_ellcard(a4,a6,p);
  }
  return subii(addiu(pr_norm(P),1), card);
}

/* a, b not both 0; sorted list of primes dividing gcd(a,b), using coprime
 * basis */
static GEN
Z_gcd_primes(GEN a, GEN b)
{
  GEN P;
  if (!signe(a))
    P = gel(absZ_factor(b), 1);
  else if (!signe(b))
    P = gel(absZ_factor(a), 1);
  else
  {
    GEN A, B, v = Z_ppio(a,b), d = gel(v,1); /* = gcd(a,b) */
    long k, l;
    if (is_pm1(d)) return cgetg(1, t_COL);
    A = gel(v,2); /* gcd(a, b^oo) */
    B = diviiexact(b, Z_ppo(b, d)); /* gcd(b, a^oo) */
    /* d = gcd(A,B) */
    P = Z_cba(A, B); /* use coprime basis to help as much as possible */
    l = lg(P);
    for (k = 1; k < l; k++) gel(P,k) = gel(Z_factor(gel(P,k)), 1);
    P = shallowconcat1(P);
    ZV_sort_inplace(P);
  }
  settyp(P, t_VEC); return P;
}
/* E/Q, integral model, Laska-Kraus-Connell algorithm. Set *pDP to a list
 * of known prime divisors of minimal discriminant */
static GEN
ellQ_minimalu(GEN E, GEN *pDP)
{
  pari_sp av;
  GEN D = ell_get_disc(E);
  GEN c4 = ell_get_c4(E);
  GEN c6 = ell_get_c6(E), g, u, P, DP;
  long l, k;

  P = Z_gcd_primes(c4, c6);
  l = lg(P); if (l == 1) { if(pDP) *pDP = P; return gen_1; }
  DP = coltrunc_init(l);
  av = avma;
  g = gcdii(sqri(c6), D);
  u = gen_1;
  for (k = 1; k < l; k++)
  {
    GEN p = gel(P, k);
    long vg = Z_pval(g, p), d = vg / 12, r = vg % 12;
    if (d) switch(itou_or_0(p))
    {
      case 2:
      {
        long a, b;
        a = Mod16( shifti(c4, -4*d) );
        b = Mod32( shifti(c6, -6*d) );
        if ((b & 3) != 3 && (a || (b && b!=8))) { d--; r += 12; }
        break;
      }
      case 3:
        if (safe_Z_lval(c6,3) == 6*d+2) { d--; r += 12; }
        break;
    }
    if (r) vectrunc_append(DP, p);
    if (d) u = mulii(u, powiu(p, d));
  }
  if (pDP) *pDP = DP;
  return gerepileuptoint(av, u);
}

/* Ensure a1 and a3 are 2-restricted and a2 is 3-restricted */
static GEN
nfrestrict23(GEN nf, GEN E)
{
  GEN a1 = nf_to_scalar_or_basis(nf, ell_get_a1(E)), A1, A2, A3, r, s, t;
  GEN a2 = nf_to_scalar_or_basis(nf, ell_get_a2(E));
  GEN a3 = nf_to_scalar_or_basis(nf, ell_get_a3(E));

  A1 = gmodgs(a1,2);
  s = gshift(gsub(A1,a1), -1);
  s = lift_if_rational(basistoalg(nf, s));
  A2 = nfsub(nf, a2, nfmul(nf,s, nfadd(nf,a1,s)));
  r = gdivgu(gsub(gmodgs(A2,3), A2), 3);
  r = lift_if_rational(basistoalg(nf, r));
  A3 = nfadd(nf, a3, nfmul(nf,r,A1));
  t = nfadd(nf, nfmul(nf, r,s), gshift(gsub(gmodgs(A3,2), A3), -1));
  t = lift_if_rational(basistoalg(nf, t));
  return mkvec4(gen_1, r, s, t);
}

static GEN
zk_capZ(GEN nf, GEN x)
{
  GEN mx = zk_scalar_or_multable(nf, x);
  return (typ(mx) == t_INT)? mx: zkmultable_capZ(mx);
}
static GEN
ellnf_c4c6_primes(GEN E)
{
  GEN nf = ellnf_get_nf(E);
  GEN c4Z = zk_capZ(nf, ell_get_c4(E));
  GEN c6Z = zk_capZ(nf, ell_get_c6(E));
  return Z_gcd_primes(c4Z, c6Z); /* primes dividing (c4,c6) \cap Z */
}
static GEN
ellnf_D_primes(GEN E)
{
  GEN nf = ellnf_get_nf(E);
  GEN P = ellnf_c4c6_primes(E);
  GEN DZ = zk_capZ(nf, ell_get_disc(E));
  long k, l = lg(P);
  for (k = 1; k < l; k++) (void)Z_pvalrem(DZ, gel(P,k), &DZ);
  if (!is_pm1(DZ))
  {
    GEN Q = gel(absZ_factor(DZ),1);
    settyp(Q, t_VEC); P = shallowconcat(P, Q); ZV_sort_inplace(P);
  }
  return P;
}

/* convert vector of localreds to NF_MINIMALPRIMES */
static GEN
Q_to_minimalprimes(GEN nf, GEN P, GEN Q)
{
  GEN L, Lr, Ls, Lt, U;
  long k, l = lg(P);
  Lr = vectrunc_init(l);
  Ls = vectrunc_init(l);
  Lt = vectrunc_init(l);
  L = vectrunc_init(l); settyp(L,t_COL);
  U = vectrunc_init(l); settyp(U,t_COL);
  for (k = 1; k < l; k++)
  {
    GEN pr = gel(P, k), q = gel(Q, k), v, u;
    long vu;
    v = gel(q,3);
    u = gel(v,1);
    vu = nfval(nf, u, pr);
    if (!vu) continue;
    vectrunc_append(Lr, gel(v,2));
    vectrunc_append(Ls, gel(v,3));
    vectrunc_append(Lt, gel(v,4));
    vectrunc_append(L, pr);
    vectrunc_append(U, stoi(vu));
  }
  return mkvec5(L, U, Lr, Ls, Lt);
}
/* E integral */
static GEN
ellminimalprimes(GEN E)
{
  GEN S, nf, c4, c6, P, Q;
  long j, k, l;

  if ((S = obj_check(E, NF_MINIMALPRIMES))) return S;
  nf = ellnf_get_nf(E);
  c4 = nf_to_scalar_or_basis(nf, ell_get_c4(E));
  c6 = nf_to_scalar_or_basis(nf, ell_get_c6(E));
  if (typ(c4) == t_INT) c4 = NULL;
  if (typ(c6) == t_INT) c6 = NULL;
  P = nf_pV_to_prV(nf, ellnf_c4c6_primes(E));
  Q = cgetg_copy(P, &l);
  for (k = j = 1; k < l; k++)
  {
    GEN pr = gel(P, k);
    if (c4 && !ZC_prdvd(c4,pr)) continue;
    if (c6 && !ZC_prdvd(c6,pr)) continue;
    gel(Q,j) = nflocalred(E, pr); /* pr | (c4,c6) */
    gel(P,j++) = pr;
  }
  setlg(P,j); setlg(Q,j);
  return obj_insert(E, NF_MINIMALPRIMES, Q_to_minimalprimes(nf,P,Q));
}
static GEN
ellnf_minimalnormu(GEN E0)
{
  GEN E, S, L, U, P, v, Nu = NULL, nf = ellnf_get_nf(E0);
  long i, l;
  E = ellintegralmodel_i(E0, &v);
  S = ellminimalprimes(E);
  L = gel(S,1);
  U = gel(S,2);
  if (v) Nu = idealnorm(nf, gel(v,1));
  P = cgetg_copy(L, &l);
  for (i = 1; i < l; i++) gel(P,i) = pr_norm(gel(L,i));
  P = factorback2(P, U);
  if (Nu) P = gmul(Nu, P);
  return P;
}
/* E integral model; return change of variable to minimal model (t_VEC)
 * or (nontrivial) Weierstrass class (t_COL) */
static GEN
bnf_get_v(GEN E)
{
  GEN bnf = ellnf_get_bnf(E);
  GEN nf, L, Lr, Ls, Lt, F, C, U, R, S, T;

  if (!bnf) pari_err_TYPE("ellminimalmodel (need a bnf)", ellnf_get_nf(E));
  S = ellminimalprimes(E);
  L = gel(S,1);
  U = gel(S,2);
  Lr = gel(S,3);
  Ls = gel(S,4);
  Lt = gel(S,5);
  F = isprincipalfact(bnf, NULL, L, U, nf_GEN);
  if (!gequal0(gel(F,1))) return gel(F,1);
  nf = bnf_get_nf(bnf);
  C = idealchinese(nf, mkmat2(L, ZC_z_mul(U,6)), NULL);
  U = basistoalg(nf, gel(F,2));
  R = basistoalg(nf, idealchinese(nf, C, Lr));
  S = basistoalg(nf, idealchinese(nf, C, Ls));
  T = basistoalg(nf, idealchinese(nf, C, Lt));
  return lift_if_rational(mkvec4(U,R,S,T));
}

GEN
ellminimaldisc(GEN E)
{
  pari_sp av = avma;
  checkell(E);
  switch(ell_get_type(E))
  {
    case t_ELL_Q:
      E = ellminimalmodel(E,NULL);
      return gerepileuptoint(av, absi_shallow(ell_get_disc(E)));
    case t_ELL_NF:
    {
      GEN nf = ellnf_get_nf(E), S, L, U, D;
      E = ellintegralmodel_i(E,NULL);
      S = ellminimalprimes(E);
      L = gel(S,1);
      U = ZC_z_mul(gel(S,2), 12);
      D = idealfactorback(nf, L, U, 0);
      return gerepileupto(av, idealdiv(nf, ell_get_disc(E), D));
    }
    default: pari_err_TYPE("ellminimaldisc (E / number field)", E);
             return NULL; /*LCOV_EXCL_LINE*/
  }
}

/* update Q_MINIMALMODEL entry in E, but don't update type-specific data on
 * ellminimalmodel(E) */
static GEN
ellminimalmodel_i(GEN E, GEN *ptv, GEN *pS)
{
  GEN S, y, e, v, v0, u, DP;
  ellmin_t M;
  if ((S = obj_check(E, Q_MINIMALMODEL)))
  {
    if (lg(S) != 2)
    {
      E = gel(S,3);
      v = gel(S,2);
    }
    else
      v = init_ch();
    if (ptv) *ptv = v;
    if (pS) *pS = S;
    return gcopy(E);
  }
  e = ellintegralmodel_i(E, &v0);
  u = ellQ_minimalu(e, &DP);
  min_set_all(&M, e, u);
  v = min_get_v(&M, e);
  y = min_to_ell(&M, e);
  if (v0) { gcomposev(&v0, v); v = v0; }
  if (is_trivial_change(v))
  {
    v = init_ch();
    S = mkvec(DP);
  }
  else
    S = mkvec3(DP, v, y);
  obj_insert(E, Q_MINIMALMODEL, S);
  if (pS) *pS = S;
  if (ptv) *ptv = v; return y;
}

static GEN
ellQminimalmodel(GEN E, GEN *ptv)
{
  pari_sp av = avma;
  GEN S, DP, v, y = ellminimalmodel_i(E, &v, &S);
  if (!is_trivial_change(v)) ch_Q(y, E, v);
  DP = gel(S,1);
  obj_insert_shallow(y, Q_MINIMALMODEL, mkvec(DP));
  if (!ptv) return gerepilecopy(av, y);
  *ptv = v; return gc_all(av, 2, &y, ptv);
}

static GEN
ellnfminimalmodel_i(GEN E, GEN *ptv)
{
  GEN S, y, v, v2;
  if ((S = obj_check(E, NF_MINIMALMODEL)))
  {
    switch(lg(S))
    {
      case 1: v = init_ch(); break;
      case 2: v = NULL; E = gel(S,1); break;
      default: E = gel(S,2); v = gel(S,1); break;
    }
    *ptv = v;
    return gcopy(E);
  }
  *ptv = NULL;
  y = ellintegralmodel_i(E, &v);
  v2 = bnf_get_v(y);
  if (typ(v2) == t_COL)
  {
    obj_insert(E, NF_MINIMALMODEL, mkvec(v2));
    return v2; /* nontrivial Weierstrass class */
  }
  y = coordch(y, v2);
  gcomposev(&v, v2);
  v2 = nfrestrict23(ellnf_get_nf(E), y);
  y = coordch(y, v2);
  /* copy to avoid inserting twice in y = E */
  y = obj_reinit(y);
  gcomposev(&v, v2);
  if (is_trivial_change(v))
  {
    v = init_ch();
    S = cgetg(1,t_VEC);
  }
  else
  {
    v = lift_if_rational(v);
    S = mkvec2(v, y);
  }
  obj_insert(E, NF_MINIMALMODEL, S);
  *ptv = v; return y;
}
static GEN
ellnfminimalmodel(GEN E, GEN *ptv)
{
  pari_sp av = avma;
  GEN v, y = ellnfminimalmodel_i(E, &v);
  if (v) obj_insert_shallow(y, NF_MINIMALMODEL, cgetg(1,t_VEC));
  if (!v || !ptv) return gerepilecopy(av, y);
  *ptv = v; return gc_all(av, 2, &y, ptv);
}
GEN
ellminimalmodel(GEN E, GEN *ptv)
{
  checkell(E);
  switch(ell_get_type(E))
  {
    case t_ELL_Q: return ellQminimalmodel(E, ptv);
    case t_ELL_NF: return ellnfminimalmodel(E, ptv);
    default: pari_err_TYPE("ellminimalmodel (E / number field)", E);
             return NULL; /*LCOV_EXCL_LINE*/
  }
}

/* return a model minimal among b models */
GEN
ellminimalbmodel(GEN e, GEN *pv)
{
  pari_sp av = avma;
  GEN f, a1, a3;

  checkell(e); f = ellminimalmodel(e, pv);
  a1 = ell_get_a1(f);
  a3 = ell_get_a3(f);
  if (!signe(a1) && !signe(a3))
  { if (!*pv) *pv = init_ch(); }
  else
  {
    GEN v = mkvec4(mpodd(a1) || mpodd(a3) ? ghalf : gen_1,
                   gen_0, gdivgs(a1,-2), gdivgs(a3,-2));
    gcomposev(pv, v); f = coordch(f, v);
  }
  if (f != e) ell_reset(f);
  return gc_all(av, 2, &f, pv);
}

/* Reduction of a rational curve E to its standard minimal model, don't
 * update type-dependant components.
 * Set v = [u, r, s, t] = change of variable E -> minimal model, with u > 0
 * Set gr = [N, [u,r,s,t], c, fa, L], where
 *   N = arithmetic conductor of E
 *   c = product of the local Tamagawa numbers cp
 *   fa = factorization of N
 *   L = list of localred(E,p) for p | N. */
static GEN
ellQ_globalred(GEN e)
{
  long k, l, iN;
  GEN S, c, E, L, P, NP, NE, D;

  E = ellminimalmodel_i(e, NULL, &S);
  P = gel(S,1); l = lg(P); /* some known prime divisors of D */
  D  = ell_get_disc(E);
  for (k = 1; k < l; k++) (void)Z_pvalrem(D, gel(P,k), &D);
  if (!is_pm1(D))
  {
    P = shallowconcat(P, gel(absZ_factor(D),1));
    ZV_sort_inplace(P);
  }
  l = lg(P); c = gen_1;
  iN = 1;
  NP = cgetg(l, t_COL);
  NE = cgetg(l, t_COL);
  L = cgetg(l, t_VEC);
  for (k = 1; k < l; k++)
  {
    GEN p = gel(P,k), q = localred(E, p), ex = gel(q,1);
    if (!signe(ex)) continue;
    gel(NP, iN) = p;
    gel(NE, iN) = ex;
    gel(L, iN) = q; iN++;
    gel(q,3) = gen_0; /*delete variable change*/
    c = mulii(c, gel(q,4));
  }
  setlg(L, iN);
  setlg(NP, iN);
  setlg(NE, iN);
  return mkvec4(factorback2(NP,NE), c, mkmat2(NP,NE), L);
}
static GEN
ellglobalred_i(GEN E)
{ return obj_checkbuild(E, Q_GLOBALRED, &ellQ_globalred); }

static GEN
Q_to_globalred(GEN nf, GEN P, GEN Q, GEN v)
{
  GEN c, L, NP, NE;
  long j, k, l = lg(P);
  c = gen_1;
  NP = cgetg(l, t_COL);
  NE = cgetg(l, t_COL);
  L = cgetg(l, t_VEC);
  for (k = j = 1; k < l; k++)
  {
    GEN p = gel(P,k), q = gel(Q,k), ex;
    ex = gel(q,1);
    if (!signe(ex)) continue;
    gel(NP, j) = p;
    gel(NE, j) = ex;
    gel(L, j) = q; j++;
    c = mulii(c, gel(q,4));
  }
  setlg(L, j); setlg(NP, j); setlg(NE, j);
  return mkvec5(idealfactorback(nf,NP,NE,0), v, c, mkmat2(NP,NE), L);
}

static GEN
ellnfglobalred(GEN E0)
{
  GEN E, P, Q, D, nf, v;
  long j, k, l;

  E = ellintegralmodel_i(E0, &v);
  if (!v) v = init_ch();
  nf = ellnf_get_nf(E);
  P = nf_pV_to_prV(nf, ellnf_D_primes(E));
  D = nf_to_scalar_or_basis(nf, ell_get_disc(E));
  if (typ(D) == t_INT) D = NULL;
  Q = cgetg_copy(P, &l);
  for (k = j = 1; k < l; k++)
  {
    GEN p = gel(P,k);
    if (D && !ZC_prdvd(D, p)) continue;
    gel(Q,j) = nflocalred(E, p);
    gel(P,j++) = p;
  }
  setlg(P,j); setlg(Q,j);
  if (!obj_check(E0, NF_MINIMALPRIMES))
    (void)obj_insert(E0, NF_MINIMALPRIMES, Q_to_minimalprimes(nf,P,Q));
  return Q_to_globalred(nf,P,Q,v);
}

GEN
ellglobalred(GEN E)
{
  pari_sp av = avma;
  GEN S, gr, v;
  checkell(E);
  switch(ell_get_type(E))
  {
    default: pari_err_TYPE("ellglobalred",E);
    case t_ELL_Q:
      gr = ellglobalred_i(E);
      S = obj_check(E, Q_MINIMALMODEL);
      v = (lg(S) == 2)? init_ch(): gel(S,2);
      v = mkvec5(gel(gr,1), v, gel(gr,2),gel(gr,3),gel(gr,4));
      break;
    case t_ELL_NF:
      v = obj_checkbuild(E, NF_GLOBALRED, &ellnfglobalred);
      break;
  }
  return gerepilecopy(av, v);
}

static GEN doellrootno(GEN e);
/* Return E = ellminimalmodel(e), but only update E[1..14].
 * insert MINIMALMODEL, GLOBALRED, ROOTNO in both e (regular insertion)
 * and E (shallow insert) */
GEN
ellanal_globalred(GEN e, GEN *ch)
{
  GEN E, S, v = NULL;
  checkell_Q(e);
  if (!(S = obj_check(e, Q_MINIMALMODEL)))
  {
    E = ellminimalmodel_i(e, &v, &S);
    obj_insert_shallow(E, Q_MINIMALMODEL, mkvec(gel(S,1)));
  }
  else if (lg(S) == 2) /* trivial change */
    E = e;
  else
  {
    v = gel(S,2);
    E = gcopy(gel(S,3));
    obj_insert_shallow(E, Q_MINIMALMODEL, mkvec(gel(S,1)));
  }
  if (ch) *ch = v;
  S = ellglobalred_i(e);
  if (E != e) obj_insert_shallow(E, Q_GLOBALRED, S);
  S = obj_check(e, Q_ROOTNO);
  if (!S)
  {
    S = doellrootno(E);
    obj_insert(e, Q_ROOTNO, S); /* insert in e */
  }
  if (E != e) obj_insert_shallow(E, Q_ROOTNO, S); /* ... and in E */
  return E;
}

static long
nb_real_components(GEN E) { return gsigne(ell_get_disc(E)) > 0? 2: 1; }
/* E minimal, \Omega_E^s in "La constante de Manin et le degre modulaire
 * d'une courbe elliptique" */
GEN
ellQtwist_bsdperiod(GEN E, long s)
{
  GEN w = ellR_omega(E,DEFAULTPREC);
  if (s == 1)
    w = gel(w,1);
  else if (nb_real_components(E) == 2)
    w = gneg(gel(w,2));
  else
    w = mkcomplex(gen_0, gneg(gmul2n(imag_i(gel(w,2)), 1)));
  return w;
}

static GEN
ellQ_tamagawa(GEN e)
{
  GEN red = ellglobalred(e), tam = gel(red,3);
  return muliu(tam, nb_real_components(e));
}

static GEN
ellnf_tamagawa(GEN e)
{
  GEN red = ellglobalred(e), tam = gel(red,3);
  GEN nf  = ellnf_get_nf(e), s = nfsign(nf, ell_get_disc(e));
  long r1, r2;
  nf_get_sign(nf, &r1, &r2);
  return shifti(tam, r2 + r1 - hammingweight(s));
}

GEN
elltamagawa(GEN E)
{
  pari_sp av = avma;
  GEN v;
  checkell(E);
  switch(ell_get_type(E))
  {
    default: pari_err_TYPE("elltamagawa",E);
    case t_ELL_Q:  v = ellQ_tamagawa(E);  break;
    case t_ELL_NF: v = ellnf_tamagawa(E); break;
  }
  return gerepileuptoint(av, v);
}

static GEN
ellnf_get_nf_prec(GEN E, long prec)
{
  GEN S, nf = ellnf_get_nf(E);
  if (nf_get_prec(nf) >= prec) return nf;
  if ((S = obj_check(E, NF_NF)) && nf_get_prec(S) >= prec) return S;
  return obj_insert(E, NF_NF, nfnewprec_shallow(nf, prec));
}
/* true nf, use nf prec */
static GEN
nfembedall(GEN nf, GEN x)
{
  long r1, r2;
  GEN cx;
  nf_get_sign(nf,&r1,&r2);
  x = nf_to_scalar_or_basis(nf,x);
  if (typ(x) != t_COL) return const_vec(r1+r2, x);
  x = Q_primitive_part(x, &cx);
  x = RgM_RgC_mul(nf_get_M(nf), x);
  if (cx) x = RgC_Rg_mul(x,cx);
  return x;
}
static long
nfembed_extraprec(GEN x)
{ long e = gexpo(x); return (e < 8)? 0: nbits2extraprec(e); }
static GEN
ellnfembed(GEN E, long prec)
{
  GEN E0, nf = ellnf_get_nf(E), Eb = cgetg(6,t_VEC), e = cgetg(6,t_VEC), L, sD;
  long prec0, r1, r2, n, i;

  nf_get_sign(nf, &r1, &r2); n = r1+r2;
  E0 = RgC_to_nfC(nf, vecslice(E,1,5));
  prec0 = prec + EXTRAPREC64;
  /* need accuracy 3b for bmodel to ensure roots are correct to b bits */
  prec += (prec0-2)*3 + nfembed_extraprec(E0);
  L =  cgetg(n+1, t_VEC);
  sD = nfeltsign(nf, ell_get_disc(E), identity_perm(r1));
  for(;;)
  {
    nf = ellnf_get_nf_prec(E, prec);
    for (i=1; i<=5; i++) gel(Eb,i) = nfembedall(nf,gel(E0,i));
    for (i=1; i<=n; i++)
    {
      GEN Ei, r;
      long j;
      for (j=1; j<=5; j++) gel(e,j) = gmael(Eb,j,i);
      gel(L,i) = Ei = ellinit_Rg(e, i<=r1? signe(gel(sD,i)): 0, prec);
      if (!Ei) break;
      r = doellR_roots_i(Ei, prec, prec0);
      if (!r) break;
    }
    if (i > n) return L;
    prec = precdbl(prec);
    if (DEBUGLEVEL>1) pari_warn(warnprec,"ellnfembed", prec);
  }
}

static GEN
ellpointnfembed(GEN E, GEN P, long prec)
{
  GEN nf = ellnf_get_nf(E), Px, Py, L;
  long i, l;
  P = RgC_to_nfC(nf, P);
  prec += nfembed_extraprec(P);
  nf = ellnf_get_nf_prec(E, prec);
  Px = nfembedall(nf, gel(P,1));
  Py = nfembedall(nf, gel(P,2));
  l = lg(Px); L =  cgetg(l, t_VEC);
  for(i = 1; i < l; i++) gel(L,i) = mkvec2(gel(Px,i), gel(Py,i));
  return L;
}

static void
ellnfembed_free(GEN L)
{
  long i, l = lg(L);
  for(i = 1; i < l; i++) obj_free(gel(L,i));
}

static GEN
ellnf_vec_wrap(GEN (*fun)(GEN, long), GEN E, long prec)
{
  pari_sp av = avma;
  GEN V = ellnfembed(E, prec);
  long i, l = lg(V);
  GEN P = cgetg(l, t_VEC);
  for(i=1; i<l; i++) gel(P,i) = fun(gel(V,i), prec);
  ellnfembed_free(V);
  return gerepilecopy(av, P);
}

GEN
ellnf_vecarea(GEN E, long prec)
{ return ellnf_vec_wrap(&ellR_area, E, prec); }

GEN
ellnf_veceta(GEN E, long prec)
{ return ellnf_vec_wrap(&ellR_eta, E, prec); }

GEN
ellnf_vecomega(GEN E, long prec)
{ return ellnf_vec_wrap(&ellR_omega, E, prec); }

static GEN
ellnf_bsdperiod(GEN E, long prec)
{
  pari_sp av = avma;
  GEN Eb = ellnfembed(E, prec), per = gtofp(ellnf_minimalnormu(E), prec);
  long i, l = lg(Eb), r1 = nf_get_r1(ellnf_get_nf(E));
  for(i = 1; i < l; i++)
  {
    GEN e = gel(Eb, i);
    GEN pi = (i <= r1)? gel(ellR_omega(e, prec),1): ellR_area(e, prec);
    per = mulrr(per, pi);
  }
  ellnfembed_free(Eb);
  return gerepileuptoleaf(av, per);
}
static GEN
ellnf_adelicvolume(GEN E, long prec)
{
  GEN t = ellnf_tamagawa(E);
  return gmul(t, ellnf_bsdperiod(E, prec));
}

static GEN
ellnf_bsd(GEN E, long prec)
{
  GEN v = ellnf_adelicvolume(E, prec);
  GEN tor = gel(elltors(E),1);
  GEN D = itor(nf_get_disc(ellnf_get_nf(E)), prec);
  return divrr(divri(v, sqri(tor)), sqrtr_abs(D));
}

static GEN
ellQ_bsd(GEN E, long prec)
{
  GEN per = gel(ellR_omega(E, prec),1);
  GEN tam = ellQ_tamagawa(E);
  GEN tor = gel(elltors(E),1);
  GEN S = obj_check(E, Q_MINIMALMODEL);
  if (lg(S) != 2)
  { /* switch to minimal model if needed */
    GEN v = gel(S,2), u = gel(v,1);
    per = gmul(per,u);
  }
  return divri(mulri(per,tam), sqri(tor));
}

GEN
ellbsd(GEN E, long prec)
{
  pari_sp av = avma;
  GEN v;
  checkell(E);
  switch(ell_get_type(E))
  {
    default: pari_err_TYPE("ellbsd",E);
    case t_ELL_Q:  v = ellQ_bsd(E, prec);  break;
    case t_ELL_NF: v = ellnf_bsd(E, prec); break;
  }
  return gerepileupto(av, v);
}

static GEN
QE_to_ZJ(GEN P)
{
  if (ell_is_inf(P))
    return mkvec3(gen_1, gen_1, gen_0);
  else
  {
    pari_sp av = avma;
    GEN D1 = denom(gel(P,1)), D2 = denom(gel(P,2));
    GEN R  = diviiexact(D2, gcdii(D1,D2));
    GEN R2 = sqri(R), R3 = mulii(R2, R);
    GEN Q1 = gmul(gel(P,1),R2);
    GEN Q2 = gmul(gel(P,2),R3);
    GEN Z  = denom(mkvec2(Q1, Q2));
    GEN Z2 = sqri(Z), Z3 = mulii(Z, Z2);
    return gerepilecopy(av, mkvec3(gmul(Q1, Z2), gmul(Q2, Z3), mulii(Z, R)));
  }
}

static GEN
QEV_to_ZJV(GEN x)
{ pari_APPLY_same(QE_to_ZJ(gel(x,i))) }

static GEN
FljV_changepointinv_pre(GEN x, GEN a4a6, ulong p, ulong pi)
{
  pari_APPLY_same(Flj_changepointinv_pre(gel(x,i), a4a6, p, pi))
}

static GEN
ellQ_factorback_filter(GEN A, GEN P, GEN *pQ)
{
  long i, j, k, l = lg(A);
  GEN B, Q;
  for (i = k = 1; i < l; i++)
    if (!ell_is_inf(gel(A,i))) k++;
  if (k == 1 || k == l) { *pQ = P; return A; }
  B = cgetg(k, t_VEC);
  Q = cgetg(k, typ(P));
  for (i = j = 1; i < l; i++)
    if (!ell_is_inf(gel(A,i))) { gel(B,j) = gel(A,i); Q[j] = P[i]; j++; }
  *pQ = Q; return B;
}

static GEN
ellQ_factorback_chinese(GEN A, GEN P, GEN *mod)
{
  GEN Q, B = ellQ_factorback_filter(A, P, &Q);
  return ncV_chinese_center(B, Q, mod);
}

static GEN
ellQ_factorback1(GEN A, GEN L, ulong l, GEN E, long CM, ulong p)
{
  pari_sp av = avma;
  ulong pi = get_Fl_red(p);
  GEN c4 = ell_get_c4(E);
  ulong a4 = Fl_c4_to_a4(Rg_to_Fl(c4, p), p);
  GEN a4a6, a, Hp;
  ulong d = 1;
  if (l != 1)
  {
    GEN c6 = ell_get_c6(E);
    ulong a6 = Fl_c6_to_a6(Rg_to_Fl(c6, p), p);
    ulong c = p + 1 - Fl_elltrace_CM(CM, a4, a6, p);
    d = Fl_invsafe(l % c, c);
    if (!d) return NULL;
  }
  a4a6 = a4a6_ch_Fl(E,p);
  a = FljV_changepointinv_pre(A, a4a6, p, pi);
  Hp = FljV_factorback_pre(a, L, a4, p, pi);
  if (d != 1)
    Hp = Flj_mulu_pre(Hp, d, a4, p, pi);
  Hp = Flj_to_Fle_pre(Hp, p, pi);
  Hp = Fle_changepoint(Hp, a4a6, p);
  return gerepileuptoleaf(av, Hp);
}

static GEN
ellQ_factorback_slice(GEN A, GEN L, ulong l, GEN E, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1;
  long CM = ellQ_get_CM(E);
  GEN H, T, B, Q;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN Hp = ellQ_factorback1(ZM_to_Flm(A, p), L, l, E, CM, p);
    if (!Hp) { *mod = gen_1; return ellinf(); }
    *mod = utoi(p);
    return Flv_to_ZV(Hp);
  }
  T = ZV_producttree(P);
  A = ZM_nv_mod_tree(A, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    gel(H,i) = ellQ_factorback1(gel(A,i), L, l, E, CM, uel(P,i));
    if (gel(H,i)==NULL) { gel(H,i) = ellinf(); uel(P,i) = 1; }
  }
  B = ellQ_factorback_filter(H, P, &Q);
  if (lg(Q) != lg(P)) T = ZV_producttree(Q);
  H = ncV_chinese_center_tree(B, Q, T, ZV_chinesetree(Q,T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ellQ_factorback_worker(GEN P, GEN E, GEN A, GEN L, ulong l)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = ellQ_factorback_slice(A, L, l, E, P, &gel(V,2));
  return V;
}

/* If a single non-zero entry, equal to 1, return its index. Else 0 */
static long
ZV_is_ei(GEN v)
{
  long i, ei = 0;
  for (i = lg(v)-1; i; i--)
    if (signe(gel(v,i)))
    {
      if (ei || !equali1(gel(v,i))) return 0;
      ei = i;
    }
  return ei;
}

/* A vector of points, L a ZV, return (sum L[i]*A[i]) / l;
 * h is the canonical height of result. Assume the result is NOT
 * torsion */
static GEN
ellQ_factorback(GEN E, GEN A, GEN L, ulong l, GEN h, long prec)
{
  pari_sp av = avma;
  GEN hn, D, worker, mod = gen_1, H = NULL;
  forprime_t S;
  ulong bound = 1;

  if (l == 1)
  {
    long i = ZV_is_ei(L);
    if (i) return gel(A,i);
  }
  hn = l==1 ? NULL: hnaive_max(E, h);
  D = ell_get_disc(E);
  worker = snm_closure(is_entry("_ellQ_factorback_worker"),
                       mkvec4(E, QEV_to_ZJV(A), L, utoi(l)));
  if (l==1)
    init_modular_big(&S);
  else
    init_modular_small(&S);
  for (bound = 1;; bound <<= 1)
  {
    GEN amax, r;
    gen_inccrt("ellQ_factorback", worker, D, bound, 0,
               &S, &H, &mod, ellQ_factorback_chinese, NULL);
    amax = sqrti(shifti(mod,-2));
    if (!ell_is_inf(H) && (r = FpC_ratlift(H, mod, amax, amax, NULL))
                       && oncurve_exact(E,r))
    {
      GEN g;
      settyp(r,t_VEC); g = ellheight(E,r,prec);
      if (signe(g) && expo(subrs(divrr(g,h),1))<-prec2nbits(prec)/2)
        return gerepileupto(av, r);
    }
    if (hn && gcmpsg(expi(mod)>>2,hn) > 0) return gc_NULL(av);
  }
}

GEN
ellQ_genreduce(GEN E, GEN G, GEN M, long prec)
{
  pari_sp av = avma;
  long i, j, l = lg(G);
  GEN L, V = cgetg(l, t_VEC);

  if (!M) M = ellheightmatrix(E, G, prec);
  L = lllgram(M); l = lg(L); /* can decrease */
  for (i = j = 1; i < l; i++)
  {
    GEN Li = gel(L, i), h = qfeval(M, Li);
    if (expo(h) > -prec2nbits(prec)/2)
      gel(V,j++) = ellQ_factorback(E, G, Li, 1, h, prec);
  }
  setlg(V, j); return gerepilecopy(av, V);
}

static long
ellQ_isdivisible_test(forprime_t *S, GEN E, long CM, GEN P, ulong l, long nb)
{
  GEN D = ell_get_disc(E);
  pari_sp av = avma;
  long m;
  for (m = 1; m <= nb; set_avma(av))
  {
    ulong o, a4, a6, p = u_forprime_next(S);
    if (dvdiu(D, p)) continue;
    Fl_ell_to_a4a6(E, p, &a4, &a6);
    o = p+1 - Fl_elltrace_CM(CM, a4, a6, p);
    if (o % l == 0)
    {
      ulong pi = get_Fl_red(p);
      GEN a4a6 = a4a6_ch_Fl(E,p);
      GEN Q = Flj_changepointinv_pre(ZV_to_Flv(P, p), a4a6, p, pi);
      GEN R = Flj_mulu_pre(Q, o/l, a4, p, pi);
      if (uel(R, 3) != 0) return 0;
      m++;
    }
  }
  return 1;
}

/* Assume l prime to 210 */
GEN
ellQ_isdivisible(GEN E, GEN P, ulong l)
{
  pari_sp av = avma;
  GEN worker, mod = gen_1, H = NULL, D = ell_get_disc(E), PJ = QE_to_ZJ(P);
  forprime_t S, U;
  long CM = ellQ_get_CM(E);
  ulong bound;

  u_forprime_init(&U, l+1, ULONG_MAX);
  if (!ellQ_isdivisible_test(&U, E, CM, PJ, l, 10)) return gc_NULL(av);
  worker = snm_closure(is_entry("_ellQ_factorback_worker"),
                       mkvec4(E, mkvec(PJ), mkvecs(1), utoi(l)));
  init_modular_small(&S);
  for (bound = 1;; bound <<= 1)
  {
    GEN amax, r;
    gen_inccrt("ellQ_factorback", worker, D, bound, 0,
               &S, &H, &mod, ellQ_factorback_chinese, NULL);
    amax = sqrti(shifti(mod,-2));
    if (!ell_is_inf(H) && (r = FpC_ratlift(H, mod, amax, amax, NULL))
                       && oncurve_exact(E,r))
    {
      settyp(r,t_VEC);
      if (gequal(ellmul(E,r,utoi(l)), P)) return gerepileupto(av, r);
      if (!ellQ_isdivisible_test(&U, E, CM, PJ, l, 10)) return gc_NULL(av);
    }
  }
}

/********************************************************************/
/**                                                                **/
/**           ROOT NUMBER (after Halberstadt at p = 2,3)           **/
/**                                                                **/
/********************************************************************/
/* x a t_INT */
static long
val_aux(GEN x, long p, long pk, long *u)
{
  long v;
  GEN z;
  if (!signe(x)) { *u = 0; return 12; }
  v = Z_lvalrem(x,p,&z);
  *u = umodiu(z,pk); return v;
}
static void
val_init(GEN e, long p, long pk,
         long *v4, long *u, long *v6, long *v, long *vD, long *d1)
{
  GEN c4 = ell_get_c4(e), c6 = ell_get_c6(e), D = ell_get_disc(e);
  pari_sp av = avma;
  *v4 = val_aux(c4, p,pk, u);
  *v6 = val_aux(c6, p,pk, v);
  *vD = val_aux(D , p,pk, d1); set_avma(av);
}

static long
kod_23(GEN e, long p)
{
  GEN S, nv;
  if ((S = obj_check(e, Q_GLOBALRED)))
  {
    GEN NP = gmael(S,3,1), L = gel(S,4);
    nv = absequaliu(gel(NP,1), p)? gel(L,1): gel(L,2); /* localred(p) */
  }
  else
    nv = localred_23(e, p);
  return itos(gel(nv,2));
}

/* v(c4), v(c6), v(D) for minimal model, +oo is coded by 12 */
static long
neron_2(long v4, long v6, long vD, long kod)
{
  if (kod > 4) return 1;
  switch(kod)
  {
    case 1: return (v6>0) ? 2 : 1;
    case 2:
      if (vD==4) return 1;
      else
      {
        if (vD==7) return 3;
        else return v4==4 ? 2 : 4;
      }
    case 3:
      switch(vD)
      {
        case 6: return 3;
        case 8: return 4;
        case 9: return 5;
        default: return v4==5 ? 2 : 1;
      }
    case 4: return v4>4 ? 2 : 1;
    case -1:
      switch(vD)
      {
        case 9: return 2;
        case 10: return 4;
        default: return v4>4 ? 3 : 1;
      }
    case -2:
      switch(vD)
      {
        case 12: return 2;
        case 14: return 3;
        default: return 1;
      }
    case -3:
      switch(vD)
      {
        case 12: return 2;
        case 14: return 3;
        case 15: return 4;
        default: return 1;
      }
    case -4: return v6==7 ? 2 : 1;
    case -5: return (v6==7 || v4==6) ? 2 : 1;
    case -6:
      switch(vD)
      {
        case 12: return 2;
        case 13: return 3;
        default: return v4==6 ? 2 : 1;
      }
    case -7: return (vD==12 || v4==6) ? 2 : 1;
    default: return v4==6 ? 2 : 1;
  }
}
/* p = 3; v(c4), v(c6), v(D) for minimal model, +oo is coded by 12 */
static long
neron_3(long v4, long v6, long vD, long kod)
{
  if (labs(kod) > 4) return 1;
  switch(kod)
  {
    case -1: case 1: return odd(v4)? 2: 1;
    case -3: case 3: return (2*v6>vD+3)? 2: 1;
    case -4: case 2:
      switch (vD%6)
      {
        case 4: return 3;
        case 5: return 4;
        default: return v6%3==1 ? 2 : 1;
      }
    default: /* kod = -2 et 4 */
      switch (vD%6)
      {
        case 0: return 2;
        case 1: return 3;
        default: return 1;
      }
  }
}

static long
ellrootno_2(GEN e)
{
  long n2, kod, u, v, x1, y1, D1, vD, v4, v6;
  long d = get_vp_u_small(e, 2, &v6, &vD);

  if (!vD) return 1;
  if (d) { /* not minimal */
    ellmin_t M;
    min_set_2(&M, e, d);
    min_set_D(&M, e);
    e = min_to_ell(&M, e);
  }
  val_init(e, 2,64,&v4,&u, &v6,&v, &vD,&D1);
  kod = kod_23(e,2);
  n2 = neron_2(v4,v6,vD, kod);
  if (kod>=5)
  {
    long a2, a3;
    a2 = ZtoF2(ell_get_a2(e));
    a3 = ZtoF2(ell_get_a3(e));
    return odd(a2 + a3) ? 1 : -1;
  }
  if (kod<-9) return (n2==2) ? -kross(-1,v) : -1;
  x1 = u+v+v;
  switch(kod)
  {
    case 1: return 1;
    case 2:
      switch(n2)
      {
        case 1:
          switch(v4)
          {
            case 4: return kross(-1,u);
            case 5: return 1;
            default: return -1;
          }
        case 2: return (v6==7) ? 1 : -1;
        case 3: return (v%8==5 || (u*v)%8==5) ? 1 : -1;
        case 4: if (v4>5) return kross(-1,v);
          return (v4==5) ? -kross(-1,u) : -1;
      }
    case 3:
      switch(n2)
      {
        case 1: return -kross(2,u*v);
        case 2: return -kross(2,v);
        case 3: y1 = (u - (v << (v6-5))) & 15;
          return (y1==7 || y1==11) ? 1 : -1;
        case 4: return (v%8==3 || (2*u+v)%8==7) ? 1 : -1;
        case 5: return v6==8 ? kross(2,x1) : kross(-2,u);
      }
    case -1:
      switch(n2)
      {
        case 1: return -kross(2,x1);
        case 2: return (v%8==7) || (x1%32==11) ? 1 : -1;
        case 3: return v4==6 ? 1 : -1;
        case 4: if (v4>6) return kross(-1,v);
          return v4==6 ? -kross(-1,u*v) : -1;
      }
    case -2: return n2==1 ? kross(-2,v) : kross(-1,v);
    case -3:
      switch(n2)
      {
        case 1: y1=(u-2*v)%64; if (y1<0) y1+=64;
          return (y1==3) || (y1==19) ? 1 : -1;
        case 2: return kross(2*kross(-1,u),v);
        case 3: return -kross(-1,u)*kross(-2*kross(-1,u),u*v);
        case 4: return v6==11 ? kross(-2,x1) : -kross(-2,u);
      }
    case -5:
      if (n2==1) return x1%32==23 ? 1 : -1;
      else return -kross(2,2*u+v);
    case -6:
      switch(n2)
      {
        case 1: return 1;
        case 2: return v6==10 ? 1 : -1;
        case 3: return (u%16==11) || ((u+4*v)%16==3) ? 1 : -1;
      }
    case -7:
      if (n2==1) return 1;
      else
      {
        y1 = (u + (v << (v6-8))) & 15;
        if (v6==10) return (y1==9 || y1==13) ? 1 : -1;
        else return (y1==9 || y1==5) ? 1 : -1;
      }
    case -8: return n2==2 ? kross(-1,v*D1) : -1;
    case -9: return n2==2 ? -kross(-1,D1) : -1;
    default: return -1;
  }
}

static long
ellrootno_3(GEN e)
{
  long n2, kod, u, v, D1, r6, K4, K6, vD, v4, v6;
  long d = get_vp_u_small(e, 3, &v6, &vD);

  if (!vD) return 1;
  if (d) { /* not minimal */
    ellmin_t M;
    min_set_3(&M, e, d);
    min_set_a(&M);
    min_set_D(&M, e);
    e = min_to_ell(&M, e);
  }
  val_init(e, 3,81, &v4,&u, &v6,&v, &vD,&D1);
  kod = kod_23(e,3);
  K6 = kross(v,3); if (kod>4) return K6;
  n2 = neron_3(v4,v6,vD,kod);
  r6 = v%9; K4 = kross(u,3);
  switch(kod)
  {
    case 1: case 3: case -3: return 1;
    case 2:
      switch(n2)
      {
        case 1: return (r6==4 || r6>6) ? 1 : -1;
        case 2: return -K4*K6;
        case 3: return 1;
        case 4: return -K6;
      }
    case 4:
      switch(n2)
      {
        case 1: return K6*kross(D1,3);
        case 2: return -K4;
        case 3: return -K6;
      }
    case -2: return n2==2 ? 1 : K6;
    case -4:
      switch(n2)
      {
        case 1:
          if (v4==4) return (r6==4 || r6==8) ? 1 : -1;
          else return (r6==1 || r6==2) ? 1 : -1;
        case 2: return -K6;
        case 3: return (r6==2 || r6==7) ? 1 : -1;
        case 4: return K6;
      }
    default: return -1;
  }
}

/* p > 3. Don't assume that e is minimal or even integral at p */
static long
ellrootno_p(GEN e, GEN p)
{
  long nuj, nuD, nu;
  GEN D = ell_get_disc(e);
  long ep, z;

  nuD = Q_pval(D, p);
  if (!nuD) return 1;
  nuj = j_pval(e, p);
  nu = (nuD - nuj) % 12;
  if (nu == 0)
  {
    GEN c6;
    long d, vg;
    if (!nuj) return 1; /* good reduction */
   /* p || N */
    c6 = ell_get_c6(e); /* != 0 */
    vg = minss(2*Q_pval(c6, p), nuD);
    d = vg / 12;
    if (d)
    {
      GEN q = powiu(p,6*d);
      c6 = (typ(c6) == t_INT)? diviiexact(c6, q): gdiv(c6, q);
    }
    if (typ(c6) != t_INT) c6 = Rg_to_Fp(c6,p);
    /* c6 in minimal model */
    return -kronecker(negi(c6), p);
  }
  if (nuj) return krosi(-1,p);
  ep = 12 / ugcd(12, nu);
  if (ep==4) z = 2; else z = odd(ep)? 3: 1;
  return krosi(-z, p);
}

static GEN
doellrootno(GEN e)
{
  GEN V, P, S = ellglobalred_i(e);
  long i, l, s = -1;

  V = obj_check(e, Q_MINIMALMODEL);
  if (lg(V) != 2) e = gel(V,3);
  P = gmael(S,3,1); l = lg(P);
  V = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(P,i);
    long t;
    switch(itou_or_0(p))
    {
      case 2: t = ellrootno_2(e); break;
      case 3: t = ellrootno_3(e); break;
      default:t = ellrootno_p(e, p);
    }
    V[i] = t; if (t < 0) s = -s;
  }
  return mkvec2(stoi(s), V);
}

/* local epsilon factor at p (over Q), including p=0 for the infinite place.
 * Global if p==1 or NULL. */
static long
ellQ_rootno(GEN e, GEN p)
{
  pari_sp av = avma;
  GEN S;
  long s;
  if (!p || isint1(p)) return ellrootno_global(e);
  if (!signe(p)) return -1; /* local factor at infinity */
  if ( (S = obj_check(e, Q_ROOTNO)) )
  {
    GEN T = obj_check(e, Q_GLOBALRED), NP = gmael(T,3,1);
    long i = ZV_search(NP, p);
    if (i) { GEN V = gel(S,2); return V[i]; }
    return 1;
  }
  switch(itou_or_0(p))
  {
    case 2:
      e = ellintegralmodel_i(e, NULL);
      s = ellrootno_2(e); break;
    case 3:
      e = ellintegralmodel_i(e, NULL);
      s = ellrootno_3(e); break;
    default:
      s = ellrootno_p(e,p); break;
  }
  return gc_long(av, s);
}

/* global root number over number field
 * Root numbers and parity of ranks of elliptic curves, Tim and Vladimir Dokchitser
 * https://arxiv.org/abs/0906.1815
 */

static GEN
ellrnfup(GEN rnf, GEN E, long prec)
{
  long i;
  GEN Eb = cgetg(6, t_VEC);
  for(i=1; i<=5; i++)
    gel(Eb, i) = rnfeltup(rnf,gel(E, i));
  return ellinit_nf(Eb, rnf_build_nfabs(rnf, prec));
}

static GEN
ellnf2isog(GEN E, GEN z)
{
  long v = fetch_var_higher();
  GEN S = deg1pol(gen_1, gneg(z), v);
  GEN E2 = ellisogeny(E, S, 1, -1, -1);
  delete_var();
  return ellinit_nf(E2, ellnf_get_nf(E));
}

static GEN
ellnf_reladelicvolume(GEN E, GEN P, GEN z, long prec)
{
  pari_sp av = avma;
  GEN nf = ellnf_get_nf(E);
  GEN rnf = rnfinit0(nf, P, 1);
  GEN Et = ellrnfup(rnf, E, prec);
  GEN E2 = ellnf2isog(Et, rnfeltreltoabs(rnf, z));
  GEN c1 = ellnf_adelicvolume(Et, prec), c2 = ellnf_adelicvolume(E2, prec);
  obj_free(rnf); obj_free(Et); obj_free(E2);
  return gerepilecopy(av, mkvec2(c1,c2));
}

static long
rootnovalp(GEN z, ulong p, long prec)
{ return mpodd(ground(gdiv(glog(z, prec), glog(utoi(p),prec)))); }

static long
ellnf_rootno_global(GEN E)
{
  pari_sp av = avma;
  GEN nf = ellnf_get_nf(E);
  long prec = nf_get_prec(nf);
  long v, var = fetch_var_higher();
  GEN F;
  E = ellintegralmodel_i(E, NULL);
  F = nfroots(nf, ec_bmodel(E, var));
  if (lg(F)>1)
  {
    GEN Et = ellnf2isog(E, gel(F,1));
    GEN cK = ellnf_adelicvolume(E, prec), cKt = ellnf_adelicvolume(Et, prec);
    obj_free(Et);
    v = rootnovalp(divrr(cK,cKt), 2, prec);
  } else
  {
    GEN D = deg2pol_shallow(gen_1, gen_0, gneg(ell_get_disc(E)), var);
    GEN P = RgX_divs(RgX_rescale(ec_bmodel(E, var), utoi(4)), 4);
    GEN c = ellnf_reladelicvolume(E, P, gmul2n(pol_x(var),-2), prec);
    GEN cL = gel(c,1), cLt = gel(c,2);
    GEN F = nfroots(nf, D);
    if (lg(F)>1)
      v = rootnovalp(divrr(cL,cLt), 2, prec);
    else
    {
      GEN cK = ellnf_adelicvolume(E, prec);
      GEN cp = nfcompositum(nf, P, D, 3);
      GEN cc = ellnf_reladelicvolume(E, gel(cp,1), gmul2n(gel(cp,2),-2), prec);
      GEN cF = gel(cc,1), cFt = gel(cc,2);
      GEN rnf = rnfinit0(nf,D,1);
      GEN Et = ellrnfup(rnf, E, prec);
      GEN cKv = ellnf_adelicvolume(Et, prec);
      long v2 = rootnovalp(divrr(gmul(cL,cF),gmul(cLt,cFt)), 2, prec);
      long v3 = rootnovalp(divrr(gmul(cF,gsqr(cK)),gmul(cKv,gsqr(cL))), 3, prec);
      obj_free(rnf); obj_free(Et);
      v = odd(v2+v3);
    }
  }
  delete_var();
  return gc_long(av, v? -1: 1);
}

static GEN
doellnfrootno(GEN e)
{ return stoi(ellnf_rootno_global(e)); }

long
ellrootno_global(GEN e)
{
  pari_sp av = avma;
  GEN S;
  switch(ell_get_type(e))
  {
    case t_ELL_Q:
      S = gel(obj_checkbuild(e, Q_ROOTNO, &doellrootno),1);
      break;
    case t_ELL_NF:
      S = obj_checkbuild(e, NF_ROOTNO, &doellnfrootno);
      break;
    default:
      pari_err_TYPE("ellrootno", e); return 0; /*LCOV_EXCL_LINE*/
  }
  return gc_long(av, itos(S));
}

long
ellrootno(GEN e, GEN p)
{
  checkell(e);
  if (p && typ(p) != t_INT) pari_err_TYPE("ellrootno", p);
  if (p && signe(p) < 0) pari_err_PRIME("ellrootno",p);
  switch(ell_get_type(e))
  {
    case t_ELL_Q:
      return ellQ_rootno(e, p);
    default: pari_err_TYPE("ellrootno", e);
    case t_ELL_NF:
      if (p) pari_err_IMPL("local root number for number fields");
      return ellrootno_global(e);
  }
}

/********************************************************************/
/**                                                                **/
/**                       TRACE OF FROBENIUS                       **/
/**                                                                **/
/********************************************************************/

/* assume p does not divide disc E */
long
ellap_CM_fast(GEN E, ulong p, long CM)
{
  ulong a4, a6;
  if (p == 2) return 3 - cardmod2(E);
  if (p == 3) return 4 - cardmod3(E);
  Fl_ell_to_a4a6(E, p, &a4, &a6);
  return Fl_elltrace_CM(CM, a4, a6, p);
}

static void
checkell_int(GEN e)
{
  checkell_Q(e);
  if (typ(ell_get_a1(e)) != t_INT ||
      typ(ell_get_a2(e)) != t_INT ||
      typ(ell_get_a3(e)) != t_INT ||
      typ(ell_get_a4(e)) != t_INT ||
      typ(ell_get_a6(e)) != t_INT) pari_err_TYPE("ellanQ [not an integral model]",e);
}

long
ellQ_get_CM(GEN e)
{
  GEN j = ell_get_j(e);
  if (typ(j) != t_INT) return 0;
  if (is_bigint(j))
  {
#ifndef LONG_IS_64BIT
    if (signe(j) < 0)
    {
      pari_sp av = avma;
      if (absequalii(j, uu32toi(0x22UL,0x45ae8000UL))) return gc_long(av,-67);
      if (absequalii(j, uu32toi(0x03a4b862,0xc4b40000UL))) return gc_long(av,-163);
    }
#endif
    return 0;
  }
  switch(signe(j))
  {
    default: return -3; /* j = 0 */
    case 1:
      switch(j[2])
      {
      case 1728: return -4;
      case 8000: return -8;
      case 54000: return -12;
      case 287496: return -16;
      case 16581375: return -28;
      default: return 0;
      }
    case -1:
      switch(j[2]) {
      case 3375: return -7;
      case 32768: return -11;
      case 884736: return -19;
      case 12288000: return -27;
      case 884736000: return -43;
#ifdef LONG_IS_64BIT
      case 147197952000L: return -67;
      case 262537412640768000L: return -163;
#endif
      default: return 0;
    }
  }
}

static long
ellnf_get_CM(GEN E)
{
  long av = avma;
  GEN j = ell_get_j(E), nf = ellnf_get_nf(E);
  GEN P = minpoly(basistoalg(nf, j), 0);
  return gc_long(av, polisclass(P));
}

long
elliscm(GEN E)
{
  checkell(E);
  switch(ell_get_type(E))
  {
    case t_ELL_Q:  return ellQ_get_CM(E);
    case t_ELL_NF: return ellnf_get_CM(E);
    default: pari_err_TYPE("elliscm", E);
             return 0; /*LCOV_EXCL_LINE*/
  }
}

/* E/Q or Qp, return cardinality including the (possible) ramified point */
static GEN
ellcard_ram(GEN E, GEN p, int *good_red)
{
  GEN a4, a6, D = Rg_to_Fp(ell_get_disc(E), p);
  if (!signe(D))
  {
    pari_sp av = avma;
    GEN ap = ellQap(E, p, good_red);
    return gerepileuptoint(av, subii(addiu(p,1), ap));
  }
  *good_red = 1;
  if (absequaliu(p,2)) return utoi(cardmod2(E));
  if (absequaliu(p,3)) return utoi(cardmod3(E));
  ell_to_a4a6(E,p,&a4,&a6);
  return Fp_ellcard(a4, a6, p);
}


/* bad reduction at p */
static void
sievep_bad(ulong p, GEN an, ulong n)
{
  ulong m, N;
  switch (an[p]) /* (-c6/p) */
  {
    case -1: /* nonsplit */
      N = n/p;
      for (m=2; m<=N; m++)
        if (an[m] != LONG_MAX) an[m*p] = -an[m];
      break;
    case 0: /* additive */
      for (m=2*p; m<=n; m+=p) an[m] = 0;
      break;
    case 1: /* split */
      N = n/p;
      for (m=2; m<=N; m++)
        if (an[m] != LONG_MAX) an[m*p] = an[m];
      break;
  }
}
/* good reduction at p */
static void
sievep_good(ulong p, GEN an, ulong n, ulong SQRTn)
{
  const long ap = an[p];
  ulong m;
  if (p <= SQRTn) {
    ulong pk, oldpk = 1;
    for (pk=p; pk <= n; oldpk=pk, pk *= p)
    {
      if (pk != p) an[pk] = ap * an[oldpk] - p * an[oldpk/p];
      for (m = n/pk; m > 1; m--)
        if (an[m] != LONG_MAX && m%p) an[m*pk] = an[m] * an[pk];
    }
  } else {
    for (m = n/p; m > 1; m--)
      if (an[m] != LONG_MAX) an[m*p] = ap * an[m];
  }
}
static void
sievep(ulong p, GEN an, ulong n, ulong SQRTn, int good_red)
{
  if (good_red)
    sievep_good(p, an, n, SQRTn);
  else
    sievep_bad(p, an, n);
}

static long
ellan_get_ap(ulong p, int *good_red, int CM, GEN e)
{
  if (!umodiu(ell_get_disc(e),p)) /* p|D, bad reduction or nonminimal model */
    return ellQap_u(e, p, good_red);
  else /* good reduction */
  {
    *good_red = 1;
    return ellap_CM_fast(e, p, CM);
  }
}
GEN
ellanQ_zv(GEN e, long n0)
{
  pari_sp av;
  ulong p, SQRTn, n = (ulong)n0;
  GEN an;
  int CM;

  if (n0 <= 0) return cgetg(1,t_VEC);
  if (n >= LGBITS)
    pari_err_IMPL( stack_sprintf("ellan for n >= %lu", LGBITS) );
  e = ellintegralmodel_i(e,NULL);
  SQRTn = usqrt(n);
  CM = ellQ_get_CM(e);

  an = const_vecsmall(n, LONG_MAX);
  an[1] = 1; av = avma;
  for (p=2; p<=n; p++)
  {
    int good_red;
    if (an[p] != LONG_MAX) continue; /* p not prime */
    an[p] = ellan_get_ap(p, &good_red, CM, e);
    sievep(p, an, n, SQRTn, good_red);
  }
  set_avma(av); return an;
}

static GEN
ellQ_eulerf(GEN e, GEN p)
{
  int good_red;
  GEN card = ellcard_ram(e, p, &good_red);
  GEN ap = subii(addiu(p, 1), card);
  if (good_red)
    return mkrfrac(gen_1,deg2pol_shallow(p, gneg(ap), gen_1, 0));
  if (!signe(ap)) return pol_1(0);
  return mkrfrac(gen_1,deg1pol_shallow(negi(ap), gen_1,0));
}

static GEN
ellanQ(GEN e, long N)
{ return vecsmall_to_vec_inplace(ellanQ_zv(e,N)); }

static GEN
ellnflocal(GEN E, GEN p, long n)
{
  pari_sp av = avma;
  GEN nf = ellnf_get_nf(E);
  GEN LP = idealprimedec_limit_f(nf, p, n ? n-1: nf_get_degree(nf)), T = NULL;
  long l = lg(LP), i;
  for (i = 1; i < l; i++)
  {
    int goodred;
    GEN P = gel(LP,i), T2;
    GEN ap = ellnfap(E, P, &goodred);
    long f = pr_get_f(P);
    if (goodred)
      T2 = mkpoln(3, pr_norm(P), negi(ap), gen_1);
    else
    {
      if (!signe(ap)) continue;
      T2 = deg1pol_shallow(negi(ap), gen_1, 0);
    }
    if (f > 1) T2 = RgX_inflate(T2, f);
    T = T? ZX_mul(T, T2): T2;
  }
  if (!T) { set_avma(av); return pol_1(0); }
  if (n==0) return gerepilecopy(av, mkrfrac(gen_1,T));
  return gerepileupto(av, RgXn_inv_i(T, n));
}

GEN
direllnf_worker(GEN P, ulong X, GEN E)
{
  pari_sp av = avma;
  long i, l = lg(P);
  GEN W = cgetg(l, t_VEC);
  for(i = 1; i < l; i++)
  {
    ulong p = uel(P,i);
    long d = ulogint(X, p) + 1; /* minimal d such that p^d > X */
    gel(W,i) = ellnflocal(E, utoi(uel(P,i)), d);
  }
  return gerepilecopy(av, mkvec2(P,W));
}

static GEN
ellnfan(GEN E, long N)
{
  GEN worker = snm_closure(is_entry("_direllnf_worker"), mkvec(E));
  return pardireuler(worker, gen_2, stoi(N), NULL, NULL);
}

GEN
elleulerf(GEN E, GEN p)
{
  checkell(E);
  switch(ell_get_type(E))
  {
    case t_ELL_Q: return ellQ_eulerf(E, p);
    case t_ELL_NF: return ellnflocal(E, p, 0);
    default:
      pari_err_TYPE("elleulerf",E);
      return NULL; /*LCOV_EXCL_LINE*/
  }
}

GEN
ellan(GEN E, long N)
{
  checkell(E);
  switch(ell_get_type(E))
  {
    case t_ELL_Q: return ellanQ(E, N);
    case t_ELL_NF: return ellnfan(E, N);
    default:
      pari_err_TYPE("ellan",E);
      return NULL; /*LCOV_EXCL_LINE*/
  }
}

static GEN
apk_good(GEN ap, GEN p, long e)
{
  GEN u, v, w;
  long j;
  if (e == 1) return ap;
  u = ap;
  w = subii(sqri(ap), p);
  for (j=3; j<=e; j++)
  {
    v = u; u = w;
    w = subii(mulii(ap,u), mulii(p,v));
  }
  return w;
}

GEN
akell(GEN e, GEN n)
{
  long i, j, s;
  pari_sp av = avma;
  GEN fa, P, E, D, u, y;

  checkell_int(e);
  if (typ(n) != t_INT) pari_err_TYPE("akell",n);
  if (signe(n)<= 0) return gen_0;
  if (gequal1(n)) return gen_1;
  D = ell_get_disc(e);
  u = Z_ppo(n, D);
  y = gen_1;
  s = 1;
  if (!equalii(u, n))
  { /* bad reduction at primes dividing n/u */
    fa = Z_factor(diviiexact(n, u));
    P = gel(fa,1);
    E = gel(fa,2);
    for (i=1; i<lg(P); i++)
    {
      GEN p = gel(P,i);
      long ex = itos(gel(E,i));
      int good_red;
      GEN ap = ellQap(e,p,&good_red);
      if (good_red) { y = mulii(y, apk_good(ap, p, ex)); continue; }
      j = signe(ap);
      if (!j) { set_avma(av); return gen_0; }
      if (odd(ex) && j < 0) s = -s;
    }
  }
  if (s < 0) y = negi(y);
  fa = Z_factor(u);
  P = gel(fa,1);
  E = gel(fa,2);
  for (i=1; i<lg(P); i++)
  { /* good reduction */
    GEN p = gel(P,i);
    GEN ap = ellap(e,p);
    y = mulii(y, apk_good(ap, p, itos(gel(E,i))));
  }
  return gerepileuptoint(av,y);
}

GEN
ellQ_get_N(GEN e)
{ GEN v = ellglobalred_i(e); return gel(v,1); }
void
ellQ_get_Nfa(GEN e, GEN *N, GEN *faN)
{ GEN v = ellglobalred_i(e); *N = gel(v,1); *faN = gel(v,3); }

GEN
elllseries(GEN e, GEN s, GEN A, long prec)
{
  pari_sp av = avma, av1;
  ulong l, n;
  long eps, flun;
  GEN z, cg, v, cga, cgb, s2, K, gs, N;

  if (!A) A = gen_1;
  else
  {
    if (gsigne(A)<=0)
      pari_err_DOMAIN("elllseries", "cut-off point", "<=", gen_0,A);
    if (gcmpgs(A,1) < 0) A = ginv(A);
  }
  if (isint(s, &s) && signe(s) <= 0) { set_avma(av); return gen_0; }
  flun = gequal1(A) && gequal1(s);
  checkell_Q(e);
  e = ellanal_globalred(e, NULL);
  N = ellQ_get_N(e);
  eps = ellrootno_global(e);
  if (flun && eps < 0) { set_avma(av); return real_0(prec); }

  gs = ggamma(s, prec);
  cg = divrr(Pi2n(1, prec), gsqrt(N,prec));
  cga = gmul(cg, A);
  cgb = gdiv(cg, A);
  l = (ulong)((prec2nbits_mul(prec, M_LN2) +
              fabs(gtodouble(real_i(s))-1.) * log(rtodbl(cga)))
            / rtodbl(cgb) + 1);
  if ((long)l < 1) l = 1;
  v = ellanQ_zv(e, minss(l,LGBITS-1));
  s2 = K = NULL; /* gcc -Wall */
  if (!flun) { s2 = gsubsg(2,s); K = gpow(cg, gsubgs(gmul2n(s,1),2),prec); }
  z = gen_0;
  av1 = avma;
  for (n = 1; n <= l; n++)
  {
    GEN p1, an, gn = utoipos(n), ns;
    an = ((ulong)n<LGBITS)? stoi(v[n]): akell(e,gn);
    if (!signe(an)) continue;

    ns = gpow(gn,s,prec);
    p1 = gdiv(incgam0(s,mulur(n,cga),gs,prec), ns);
    if (flun)
      p1 = gmul2n(p1, 1);
    else
    {
      GEN p2 = gdiv(gmul(gmul(K,ns), incgam(s2,mulur(n,cgb),prec)), sqru(n));
      if (eps < 0) p2 = gneg_i(p2);
      p1 = gadd(p1, p2);
    }
    z = gadd(z, gmul(p1, an));
    if (gc_needed(av1,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lseriesell");
      z = gerepilecopy(av1,z);
    }
  }
  return gerepileupto(av, gdiv(z,gs));
}

/********************************************************************/
/**                                                                **/
/**                       CANONICAL HEIGHT                         **/
/**                                                                **/
/********************************************************************/

static GEN
ellnf_volume(GEN e, long prec)
{
  GEN V = ellnf_vecarea(e,prec);
  long i, r1 = nf_get_r1(ellnf_get_nf(e)), l = lg(V);
  GEN r = gen_1;
  for(i=1; i <= r1; i++) r = gmul(r, gel(V,i));
  for(   ; i < l  ; i++) r = gmul(r, gsqr(gel(V,i)));
  return r;
}

/* The function follows
<https://publications.ias.edu/sites/default/files/Number52.pdf>
<https://resnumtheor.springeropen.com/track/pdf/10.1007/s40993-017-0077-7>
*/

static GEN
ellheightfaltings(GEN e, long prec)
{
  GEN h;
  long d;
  pari_sp av = avma;
  checkell(e);
  switch(ell_get_type(e))
  {
    case t_ELL_Q:
      d = 1; e = ellintegralmodel_i(e,NULL);
      h = gmul(gsqr(ellQ_minimalu(e,NULL)), ellR_area(e, prec));
      break;
    case t_ELL_NF:
      d = nf_get_degree(ellnf_get_nf(e));
      h = gmul(gsqr(ellnf_minimalnormu(e)), ellnf_volume(e, prec));
      break;
    default:
      pari_err_TYPE("ellheight", e);
      return NULL; /*LCOV_EXCL_LINE*/
  }
  return gerepileupto(av, gdivgs(logr_abs(h), -2*d));
}

static GEN
Q_numer(GEN x) { return typ(x) == t_INT? x: gel(x,1); }

/* one root of X^2 - t X + c */
static GEN
quad_root(GEN t, GEN c, long prec)
{
  return gmul2n(gadd(t, gsqrt(gsub(gsqr(t), gmul2n(c,2)),prec)), -1);
}

/* exp( h_oo(z) ), assume z on neutral component.
 * If flag, return exp(4 h_oo(z)) instead */
static GEN
exphellagm(GEN e, GEN z, int flag, long prec)
{
  GEN x_a, ab, a, b, e1, r, V = cgetg(1, t_VEC), x = gel(z,1);
  long n, ex = 5-prec2nbits(prec), p = prec+EXTRAPREC64;

  if (typ(x) == t_REAL && realprec(x) < p) x = gprec_w(x, p);
  ab = ellR_ab(e, p);
  a = gel(ab, 1);
  b = gel(ab, 2);
  e1= gel(obj_check(e,R_ROOTS), 1); /* use maximal accuracy, don't truncate */
  x = gsub(x, e1);
  x = quad_root(gadd(x,b), gmul(a,x), prec);

  x_a = gsub(x, a);
  if (gsigne(a) > 0) { GEN a0=a; x = gsub(x, b); a = gneg(b); b = gsub(a0, b); }
  a = gsqrt(gneg(a), prec);
  b = gsqrt(gneg(b), prec);
  /* compute height on isogenous curve E1 ~ E0 */
  for(n=0;; n++)
  {
    GEN p1, p2, ab, a0 = a;
    a = gmul2n(gadd(a0,b), -1);
    r = gsub(a, a0);
    if (gequal0(r) || gexpo(r) < ex) break;
    ab = gmul(a0, b);
    b = gsqrt(ab, prec);

    p1 = gmul2n(gsub(x, ab), -1);
    p2 = gsqr(a);
    x = gadd(p1, gsqrt(gadd(gsqr(p1), gmul(x, p2)), prec));
    V = shallowconcat(V, gadd(x, p2));
  }
  if (n) {
    x = gel(V,n);
    while (--n > 0) x = gdiv(gsqr(x), gel(V,n));
  } else
    x = gadd(x, gsqr(a));
  /* height on E1 is log(x)/2. Go back to E0 */
  return flag? gsqr(gdiv(gsqr(x), x_a)): gdiv(x, sqrtr(mpabs_shallow(x_a)));
}
/* is P \in E(R)^0, the neutral component ? */
static int
ellR_on_neutral(GEN E, GEN P, long prec)
{
  GEN x = gel(P,1), e1 = ellR_root(E, prec);
  return gcmp(x, e1) >= 0;
}

/* hoo + 1/2 log(den(x)) */
static GEN
hoo_aux(GEN E, GEN z, GEN d, long prec)
{
  pari_sp av = avma;
  GEN h;
  if (!ellR_on_neutral(E, z, prec))
  {
    GEN eh = exphellagm(E, elladd(E, z,z), 0, prec);
    /* h_oo(2P) = 4h_oo(P) - log |2y + a1x + a3| */
    h = gmul(eh, gabs(ec_dmFdy_evalQ(E, z), prec));
  }
  else
    h = exphellagm(E, z, 1, prec);
  if (!is_pm1(d)) h = gmul(h, sqri(d));
  return gerepileuptoleaf(av, gmul2n(mplog(h), -2));
}
GEN
ellheightoo(GEN E, GEN z, long prec) { return hoo_aux(E, z, gen_1, prec); }

/* Formula from Silverman GTM 151 Theorem 3.2 page 466 */
static GEN
ellheight_C(GEN E, GEN P, long prec)
{
  pari_sp av = avma;
  GEN z = zell(E, P, prec);
  GEN per = ellperiods(E, 1, prec);
  GEN w = gel(per,1), w1 = gel(w,1), w2 = gel(w, 2), w1c = conj_i(w1);
  GEN e = gel(per,2), e1 = gel(e,1), e2 = gel(e, 2);
  GEN D = gsub(gmul(w1, conj_i(w2)),gmul(w1c, w2));
  GEN b = gdiv(gsub(gmul(w1, conj_i(z)),gmul(w1c, z)), D);
  GEN a = gdiv(gsub(z, gmul(b, w2)), w1);
  GEN eta = gadd(gmul(a, e1), gmul(b, e2));
  GEN r = gmul2n(real_i(gmul(z, eta)), -1);
  GEN l = real_i(ellsigma(per, z, 1, prec));
  return gerepileupto(av, gsub(r, l));
}

static GEN
_hell(GEN E, GEN p, long n, GEN P)
{ return p? ellpadicheight(E,p,n, P): ellheight(E,P,n); }
static GEN
ellheightpairing(GEN E, GEN p, long n, GEN P, GEN Q)
{
  pari_sp av = avma;
  GEN a = _hell(E,p,n, elladd(E,P,Q));
  GEN b = _hell(E,p,n, ellsub(E,P,Q));
  return gerepileupto(av, gmul2n(gsub(a,b), -2));
}
GEN
ellheight0(GEN e, GEN a, GEN b, long n)
{
  if (!a)
  {
    if (b) pari_err(e_MISC, "cannot omit P and set Q");
    return ellheightfaltings(e,n);
  }
  return b? ellheightpairing(e,NULL,n, a,b): ellheight(e,a,n);
}
GEN
ellpadicheight0(GEN e, GEN p, long n, GEN P, GEN Q)
{ return Q? ellheightpairing(e,p,n, P,Q): ellpadicheight(e,p,n, P); }

static GEN
ellnf_localheight(GEN e, GEN P, GEN pr)
{
  long v1, v2, vD, vu;
  GEN nf = ellnf_get_nf(e);
  GEN lr = nflocalred(e,pr);
  GEN k = gel(lr, 2), urst = gel(lr, 3), u = gel(urst, 1);
  GEN E = ellchangecurve(e, urst);
  GEN Q = ellchangepoint(P, urst);
  GEN v;
  vu = nfval(nf, u, pr);
  v1 = nfval(nf, ec_dFdx_evalQ(E, Q), pr);
  v2 = nfval(nf, ec_dmFdy_evalQ(E, Q), pr);
  vD = nfval(nf, ell_get_disc(E), pr); /* >= 0 */
  if (v1<0)
    vu = 0;
  if (v1<=0 || v2<=0)
    v = gen_0;
  else if (cmpis(k,5) >= 0)
  {
    GEN a = uutoQ(minss(2*v2,vD), 2*vD);
    v = gmul(gsub(gsqr(a),a), uutoQ(vD,2));
  }
  else
  {
    long v3 = nfval(nf, ec_3divpol_evalx(E, gel(Q,1)), pr);
    v = (v2<LONG_MAX && v3>=3*v2) ? sstoQ(-v2,3):
                                    sstoQ(-v3,8);
  }
  return gsubgs(v,vu);
}

static GEN
ellnf_height(GEN E, GEN P, long prec)
{
  pari_sp av = avma;
  GEN x, nf, disc, d, F, Ee, Pe, s;
  long i, n, l, r1;
  if (signe(ellorder(E, P, NULL))) return gen_0;
  x = gel(P,1);
  if (gequal0(ec_2divpol_evalx(E, x))) { set_avma(av); return gen_0; }
  nf = ellnf_get_nf(E); r1 = nf_get_r1(nf);
  disc = ell_get_disc(E);
  d = idealnorm(nf, gel(idealnumden(nf, x), 2));
  F = gel(idealfactor(nf, disc), 1);
  Ee = ellnfembed(E, prec);
  Pe = ellpointnfembed(E, P, prec);
  n = lg(Ee); l = lg(F);
  s = gmul2n(glog(d, prec), -1);
  for (i=1; i<=r1; i++)
    s = gadd(s, ellheightoo(gel(Ee, i), gel(Pe, i), prec));
  for (   ; i<n; i++)
    s = gadd(s, gmul2n(ellheight_C(gel(Ee, i), gel(Pe, i), prec), 1));
  for (i=1; i<l; i++)
  {
    GEN pr = gel(F,i), p = pr_get_p(pr);
    long f = pr_get_f(pr);
    GEN lam = ellnf_localheight(E, P, pr);
    s = gadd(s, gmul(lam, mulrs(glog(p, prec), f)));
  }
  return gerepileupto(av, gmul2n(s, 1));
}

static GEN
ellQ_height(GEN e, GEN a, long prec)
{
  long i, lx, newell = 0;
  pari_sp av;
  GEN Lp, x, y, z, phi2, psi2, psi3;
  GEN v, S, b2, b4, b6, b8, a1, a2, a4, c4, D;

  if (!RgV_is_QV(a)) pari_err_TYPE("ellheight [not a rational point]",a);
  if (ellorder_Q(e, a)) return gen_0;
  av = avma;
  if ((S = obj_check(e, Q_MINIMALMODEL)))
  { /* switch to minimal model if needed */
    if (lg(S) != 2)
    {
      v = gel(S,2);
      e = gel(S,3);
      a = ellchangepoint(a, v);
    }
  }
  else
  {
    newell = 1;
    e = ellminimalmodel_i(e, &v, NULL);
    a = ellchangepoint(a, v);
  }
  if (!oncurve(e,a))
    pari_err_DOMAIN("ellheight", "point", "not on", strtoGENstr("E"),a);
  psi2 = Q_numer(ec_dmFdy_evalQ(e,a));
  if (!signe(psi2)) { set_avma(av); return gen_0; }
  x = gel(a,1);
  y = gel(a,2);
  b2 = ell_get_b2(e);
  b4 = ell_get_b4(e);
  b6 = ell_get_b6(e);
  b8 = ell_get_b8(e);
  psi3 = Q_numer( /* b8 + 3x b6 + 3x^2 b4 + x^3 b2 + 3 x^4 */
    poleval(mkvec5(b8, mului(3,b6), mului(3,b4), b2, utoipos(3)), x)
  );
  if (!signe(psi3)) { set_avma(av); return gen_0; }
  a1 = ell_get_a1(e);
  a2 = ell_get_a2(e);
  a4 = ell_get_a4(e);
  phi2 = Q_numer( /* a4 + 2a2 x + 3x^2 - y a1*/
    poleval(mkvec3(gsub(a4,gmul(a1,y)), shifti(a2,1), utoipos(3)), x)
  );
  c4 = ell_get_c4(e);
  D = ell_get_disc(e);
  z = hoo_aux(e,a,Q_denom(x),prec);  /* hoo(a) + log(den(x))/2 */
  Lp = gel(Z_factor(gcdii(psi2,phi2)),1);
  lx = lg(Lp);
  for (i=1; i<lx; i++)
  {
    GEN p = gel(Lp,i);
    long u, v, n, n2;
    if (!dvdii(c4,p))
    { /* p \nmid c4 */
      long N = Z_pval(D,p);
      if (!N) continue;
      n2 = Z_pval(psi2,p); n = n2<<1;
      if (n > N) n = N;
      u = n * ((N<<1) - n);
      v = N << 3;
    }
    else
    {
      n2 = Z_pval(psi2, p);
      n  = Z_pval(psi3, p);
      if (n >= 3*n2) { u = n2; v = 3; } else { u = n; v = 8; }
    }
    /* z -= u log(p) / v */
    z = gsub(z, divru(mulur(u, logr_abs(itor(p,prec))), v));
  }
  if (newell) obj_free(e);
  return gerepileupto(av, gmul2n(z, 1));
}

GEN
ellheight(GEN e, GEN a, long prec)
{
  checkell(e); checkellpt(a);
  switch(ell_get_type(e))
  {
    case t_ELL_Q:
      return ellQ_height(e, a, prec);
    default: pari_err_TYPE("ellheight", e);
    case t_ELL_NF:
      return ellnf_height(e, a, prec);
  }
}

GEN
ellpadicheightmatrix(GEN e, GEN p, long n, GEN x)
{
  GEN D, A, B;
  long lx = lg(x), i, j;
  pari_sp av = avma;

  if (!is_vec_t(typ(x))) pari_err_TYPE("ellheightmatrix",x);
  D = cgetg(lx,t_VEC);
  A = cgetg(lx,t_MAT);
  B = cgetg(lx,t_MAT);
  for (i=1; i<lx; i++)
  {
    gel(D,i) = _hell(e,p,n, gel(x,i));
    gel(A,i) = cgetg(lx,t_COL);
    gel(B,i) = cgetg(lx,t_COL); /*unused if p = NULL */
  }
  for (i=1; i<lx; i++)
  {
    GEN h = gel(D,i);
    if (p)
    {
      gcoeff(A,i,i) = gel(h,1);
      gcoeff(B,i,i) = gel(h,2);
    }
    else
      gcoeff(A,i,i) = h;
    for (j=i+1; j<lx; j++)
    {
      h = _hell(e,p,n, elladd(e,gel(x,i),gel(x,j)));
      h = gmul2n(gsub(h, gadd(gel(D,i),gel(D,j))), -1);
      if (p)
      {
        gcoeff(A,j,i) = gcoeff(A,i,j) = gel(h,1);
        gcoeff(B,j,i) = gcoeff(B,i,j) = gel(h,2);
      }
      else
        gcoeff(A,j,i) = gcoeff(A,i,j) = h;
    }
  }
  return gerepilecopy(av, p? mkvec2(A,B): A);
}
GEN
ellheightmatrix(GEN E, GEN x, long n)
{ return ellpadicheightmatrix(E,NULL,n, x); }

/* Q an actual point, P a point or vector/matrix of points */
static GEN
bilhell_i(GEN E, GEN P, GEN Q, long n)
{
  GEN y;
  long i, l = lg(P);
  if (l==1) return cgetg(1,typ(P));
  if (!is_matvec_t( typ(gel(P,1)) )) return ellheight0(E,P,Q,n);
  y = cgetg(l, typ(P));
  for (i=1; i<l; i++) gel(y,i) = bilhell_i(E,gel(P,i),Q,n);
  return y;
}
GEN
bilhell(GEN E, GEN P, GEN Q, long n)
{
  long t1 = typ(P), t2 = typ(Q);
  if (!is_matvec_t(t1)) pari_err_TYPE("ellbil",P);
  if (!is_matvec_t(t2)) pari_err_TYPE("ellbil",Q);
  if (lg(P)==1) return cgetg(1,t1);
  if (lg(Q)==1) return cgetg(1,t2);
  t2 = typ(gel(Q,1));
  if (is_matvec_t(t2))
  {
    t1 = typ(gel(P,1));
    if (is_matvec_t(t1)) pari_err_TYPE("bilhell",P);
    return bilhell_i(E,Q,P,n);
  }
  return bilhell_i(E,P,Q,n);
}
/********************************************************************/
/**                                                                **/
/**                    Modular Parametrization                     **/
/**                                                                **/
/********************************************************************/
/* t*x^v (1 + O(x)), t != 0 */
static GEN
triv_ser(GEN t, long v)
{
  GEN s = cgetg(3,t_SER);
  s[1] = evalsigne(1) | _evalvalser(v) | evalvarn(0);
  gel(s,2) = t; return s;
}

GEN
elltaniyama(GEN e, long prec)
{
  GEN x, w, c, d, X, C, b2, b4;
  long n, m;
  pari_sp av = avma;

  checkell_Q(e);
  if (prec < 0) pari_err_DOMAIN("elltaniyama","precision","<",gen_0,stoi(prec));
  if (!prec) retmkvec2(triv_ser(gen_1,-2), triv_ser(gen_m1,-3));

  x = cgetg(prec+3,t_SER);
  x[1] = evalsigne(1) | _evalvalser(-2) | evalvarn(0);
  d = ginv(RgV_to_ser(ellanQ(e,prec+1), 0, prec+3)); setvalser(d,-1);
  /* 2y(q) + a1x + a3 = d qx'(q). Solve for x(q),y(q):
   * 4y^2 = 4x^3 + b2 x^2 + 2b4 x + b6 */
  c = gsqr(d);
  /* solve 4x^3 + b2 x^2 + 2b4 x + b6 = c (q x'(q))^2; c = 1/q^2 + O(1/q)
   * Take derivative then divide by 2x':
   *  b2 x + b4 = (1/2) (q c')(q x') + c q (q x')' - 6x^2.
   * Write X[i] = coeff(x, q^i), C[i] = coeff(c, q^i), we obtain for all n
   *  ((n+1)(n+2)-12) X[n+2] =  b2 X[n] + b4 delta_{n = 0}
   *   + 6    \sum_{m = -1}^{n+1} X[m] X[n-m]
   *   - (1/2)\sum_{m = -2}^{n+1} (n+m) m C[n-m]X[m].
   * */
  C = c+4;
  X = x+4;
  gel(X,-2) = gen_1;
  gel(X,-1) = gmul2n(gel(C,-1), -1); /* n = -3, X[-1] = C[-1] / 2 */
  b2 = ell_get_b2(e);
  b4 = ell_get_b4(e);
  for (n=-2; n <= prec-4; n++)
  {
    pari_sp av2 = avma;
    GEN s1, s2, s3;
    if (n != 2)
    {
      s3 = gmul(b2, gel(X,n));
      if (!n) s3 = gadd(s3, b4);
      s2 = gen_0;
      for (m=-2; m<=n+1; m++)
        if (m) s2 = gadd(s2, gmulsg(m*(n+m), gmul(gel(X,m),gel(C,n-m))));
      s2 = gmul2n(s2,-1);
      s1 = gen_0;
      for (m=-1; m+m < n; m++) s1 = gadd(s1, gmul(gel(X,m),gel(X,n-m)));
      s1 = gmul2n(s1, 1);
      if (m+m==n) s1 = gadd(s1, gsqr(gel(X,m)));
      /* ( (n+1)(n+2) - 12 ) X[n+2] = (6 s1 + s3 - s2) */
      s1 = gdivgs(gsub(gadd(gmulsg(6,s1),s3),s2), (n+2)*(n+1)-12);
    }
    else
    {
      GEN b6 = ell_get_b6(e);
      GEN U = cgetg(9, t_SER);
      U[1] = evalsigne(1) | _evalvalser(-2) | evalvarn(0);
      gel(U,2) = gel(x,2);
      gel(U,3) = gel(x,3);
      gel(U,4) = gel(x,4);
      gel(U,5) = gel(x,5);
      gel(U,6) = gel(x,6);
      gel(U,7) = gel(x,7);
      gel(U,8) = gen_0; /* defined mod q^5 */
      /* write x = U + x_4 q^4 + O(q^5) and expand original equation */
      w = derivser(U); setvalser(w,-2); /* q X' */
      /* 4X^3 + b2 U^2 + 2b4 U + b6 */
      s1 = gadd(b6, gmul(U, gadd(gmul2n(b4,1), gmul(U,gadd(b2,gmul2n(U,2))))));
      /* s2 = (qX')^2 - (4X^3 + b2 U^2 + 2b4 U + b6) = 28 x_4 + O(q) */
      s2 = gsub(gmul(c,gsqr(w)), s1);
      s1 = signe(s2)? gdivgu(gel(s2,2), 28): gen_0; /* = x_4 */
    }
    gel(X,n+2) = gerepileupto(av2, s1);
  }
  w = gmul(d,derivser(x)); setvalser(w, valser(w)+1);
  w = gsub(w, ec_h_evalx(e,x));
  c = cgetg(3,t_VEC);
  gel(c,1) = gcopy(x);
  gel(c,2) = gmul2n(w,-1); return gerepileupto(av, c);
}

/********************************************************************/
/**                                                                **/
/**                       TORSION POINTS (over Q)                  **/
/**                                                                **/
/********************************************************************/
static GEN
doellff_get_o(GEN E)
{
  GEN G = ellff_get_group(E), d = (lg(G) == 1)? gen_1: gel(G,1);
  return mkvec2(d, Z_factor(d));
}
GEN
ellff_get_o(GEN E)
{ return obj_checkbuild(E, FF_O, &doellff_get_o); }

static void
RgE2_Fp_init(GEN E, GEN *pP, GEN *pQ, GEN *a4, GEN p)
{
  GEN e = ellff_get_a4a6(E);
  *a4 = gel(e, 1);
  *pP = FpE_changepointinv(RgE_to_FpE(*pP,p), gel(e,3), p);
  *pQ = FpE_changepointinv(RgE_to_FpE(*pQ,p), gel(e,3), p);
}
GEN
elllog(GEN E, GEN a, GEN g, GEN o)
{
  pari_sp av = avma;
  GEN p;
  checkell_Fq(E); checkellpt(a); checkellpt(g);
  p = ellff_get_field(E);
  if (!o) o = ellff_get_o(E);
  if (typ(p)==t_FFELT) return FF_elllog(E, a, g, o);
  else
  {
    GEN a4;
    RgE2_Fp_init(E, &a, &g, &a4, p);
    return gerepileuptoint(av, FpE_log(a, g, o, a4, p));
  }
}

GEN
ellweilpairing(GEN E, GEN P, GEN Q, GEN m)
{
  GEN p;
  checkell_Fq(E); checkellpt(P); checkellpt(Q);
  if (typ(m)!=t_INT) pari_err_TYPE("ellweilpairing",m);
  p = ellff_get_field(E);
  if (typ(p)==t_FFELT) return FF_ellweilpairing(E, P, Q, m);
  else
  {
    pari_sp av = avma;
    GEN w, a4;
    RgE2_Fp_init(E, &P, &Q, &a4, p);
    w = FpE_weilpairing(P, Q, m, a4, p);
    return gerepileupto(av, Fp_to_mod(w, p));
  }
}

GEN
elltatepairing(GEN E, GEN P, GEN Q, GEN m)
{
  GEN p;
  checkell_Fq(E); checkellpt(P); checkellpt(Q);
  if (typ(m)!=t_INT) pari_err_TYPE("elltatepairing",m);
  p = ellff_get_field(E);
  if (typ(p)==t_FFELT) return FF_elltatepairing(E, P, Q, m);
  else
  {
    pari_sp av = avma;
    GEN t, a4;
    RgE2_Fp_init(E, &P, &Q, &a4, p);
    t = FpE_tatepairing(P, Q, m, a4, p);
    return gerepileupto(av, Fp_to_mod(t, p));
  }
}

GEN
ellap(GEN E, GEN p)
{
  pari_sp av = avma;
  GEN q, card;
  int goodred;
  p = checkellp(&E, p, NULL, "ellap");
  switch(ell_get_type(E))
  {
  case t_ELL_Fp:
    q = p; card = ellff_get_card(E);
    break;
  case t_ELL_Fq:
    q = FF_q(ellff_get_field(E)); card = ellff_get_card(E);
    break;
  case t_ELL_Qp:
  case t_ELL_Q:
    q = p; card = ellcard_ram(E, p, &goodred);
    break;
  case t_ELL_NF:
    return ellnfap(E, p, &goodred);
  default:
    pari_err_TYPE("ellap",E);
    return NULL; /*LCOV_EXCL_LINE*/
  }
  return gerepileuptoint(av, subii(addiu(q,1), card));
}

/* N.B. q > minq, then the list of potential orders in ellsea will not contain
 * an ambiguity => oo-loop. E.g. ellsea(ellinit([1,519],523)) */
GEN
ellsea(GEN E, long smallfact)
{
  const ulong minq = 523;
  checkell_Fq(E);
  switch(ell_get_type(E))
  {
  case t_ELL_Fp:
    {
      GEN p = ellff_get_field(E), e = ellff_get_a4a6(E);
      if (abscmpiu(p, minq) <= 0) return Fp_ellcard(gel(e,1), gel(e,2), p);
      return Fp_ellcard_SEA(gel(e,1), gel(e,2), p, smallfact);
    }
  case t_ELL_Fq:
    {
      GEN fg = ellff_get_field(E);
      if (abscmpiu(FF_p_i(fg), 7) <= 0 || abscmpiu(FF_q(fg), minq) <= 0)
        return FF_ellcard(E);
      return FF_ellcard_SEA(E, smallfact);
    }
  }
  return NULL; /*LCOV_EXCL_LINE*/
}

GEN
ellff_get_card(GEN E)
{ return obj_checkbuild(E, FF_CARD, &doellcard); }

GEN
ellcard(GEN E, GEN p)
{
  p = checkellp(&E, p, NULL, "ellcard");
  switch(ell_get_type(E))
  {
  case t_ELL_Fp: case t_ELL_Fq:
    return icopy(ellff_get_card(E));
  case t_ELL_Qp:
  case t_ELL_Q:
    {
      pari_sp av = avma;
      int goodred;
      GEN N = ellcard_ram(E, p, &goodred);
      if (!goodred) N = subiu(N, 1); /* remove singular point */
      return gerepileuptoint(av, N);
    }
  case t_ELL_NF:
    {
      pari_sp av = avma;
      int goodred;
      GEN N = subii(pr_norm(p), ellnfap(E, p, &goodred));
      if (goodred) N = addiu(N, 1);
      return gerepileuptoint(av, N);
    }
  default:
    pari_err_TYPE("ellcard",E);
    return NULL; /*LCOV_EXCL_LINE*/
  }
}

/* assume model is p-minimal */
static GEN
ellgroup_m(GEN E, GEN p, GEN *pm)
{
  GEN a4, a6, N = ellcard(E, p); /* #E^ns(Fp) */
  *pm = gen_1;
  if (equali1(N)) return cgetg(1,t_VEC);
  if (absequaliu(p, 2)) return mkvec(N);
  if (absequaliu(p, 3))
  { /* The only possible noncyclic group is [2,2] which happens 9 times */
    ulong b2, b4, b6;
    if (!absequaliu(N, 4)) return mkvec(N);
    /* If the group is not cyclic, T = 4x^3 + b2 x^2 + 2b4 x + b6
     * must have 3 roots else 1 root. Test T(0) = T(1) = 0 mod 3 */
    b6 = Rg_to_Fl(ell_get_b6(E), 3);
    if (b6) return mkvec(N);
    /* b6 = T(0) = 0 mod 3. Test T(1) */
    b2 = Rg_to_Fl(ell_get_b2(E), 3);
    b4 = Rg_to_Fl(ell_get_b4(E), 3);
    if ((1 + b2 + (b4<<1)) % 3) return mkvec(N);
    return mkvec2s(2, 2);
  } /* Now assume p > 3 */
  ell_to_a4a6(E, p, &a4,&a6);
  return Fp_ellgroup(a4,a6,N,p, pm);
}

static GEN
doellGm(GEN E)
{
  GEN fg = ellff_get_field(E);
  GEN m, G = (typ(fg) == t_FFELT)? FF_ellgroup(E, &m): ellgroup_m(E, fg, &m);
  return mkvec2(G, m);
}
static GEN
ellff_Gm(GEN E)
{ return obj_checkbuild(E, FF_GROUP, &doellGm); }
GEN
ellff_get_group(GEN E) { return gel(ellff_Gm(E), 1); }
GEN
ellff_get_m(GEN E) { return gel(ellff_Gm(E), 2); }
GEN
ellff_get_D(GEN E)
{
  GEN G = ellff_get_group(E), o = ellff_get_o(E);
  switch(lg(G))
  {
    case 1: return G;
    case 2: return mkvec(o);
    default: return mkvec2(o, gel(G,2));
  }
}

/* E / Fp */
static GEN
doellgens(GEN E)
{
  GEN fg = ellff_get_field(E);
  if (typ(fg)==t_FFELT)
    return FF_ellgens(E);
  else
  {
    GEN F, p = fg, e = ellff_get_a4a6(E);
    F = Fp_ellgens(gel(e,1),gel(e,2),gel(e,3), ellff_get_D(E),ellff_get_m(E),p);
    return FpVV_to_mod(F,p);
  }
}

GEN
ellff_get_gens(GEN E)
{ return obj_checkbuild(E, FF_GROUPGEN, &doellgens); }

GEN
ellgroup(GEN E, GEN p)
{
  pari_sp av = avma;
  GEN m, G;
  p = checkellp(&E,p, NULL, "ellgroup");
  switch(ell_get_type(E))
  {
    case t_ELL_Fp:
    case t_ELL_Fq: G = ellff_get_group(E); break;
    case t_ELL_Qp:
    case t_ELL_Q:
      if (Z_pval(Q_numer(ell_get_disc(E)), p))
      {
        GEN Q = localred(E,p), kod = gel(Q,2);
        E = ellchangecurve(E, gel(Q,3));
        if (!equali1(kod)) { G = mkvec(ellcard(E,p)); break; }
      }
      G = ellgroup_m(E,p,&m); break;
    case t_ELL_NF:
      if (nfval(ellnf_get_nf(E), ell_get_disc(E), p))
      {
        GEN Q = nflocalred(E,p), kod = gel(Q,2);
        E = ellchangecurve(E, gel(Q,3));
        if (!equali1(kod)) { G = mkvec(ellcard(E,p)); break; }
      }
      E = ellinit(E, p, 0);
      G = ellff_get_group(E);
      G = gcopy(G); obj_free(E); break;
    default:
      pari_err_TYPE("ellgroup", E);
      return NULL;/*LCOV_EXCL_LINE*/
  }
  return gerepilecopy(av, G);
}

GEN
ellgroup0(GEN E, GEN p, long flag)
{
  pari_sp av = avma;
  long tE, freeE = 0;
  GEN G;
  if (flag==0) return ellgroup(E, p);
  if (flag!=1) pari_err_FLAG("ellgroup");
  checkell(E); tE = ell_get_type(E);
  if (tE != t_ELL_Fp && tE != t_ELL_Fq)
  {
    GEN Q = elllocalred(E, p), v = gel(Q,3), u = gel(v,1), kod = gel(Q,2);
    long vu;
    switch(tE)
    {
      case t_ELL_Qp: p = ellQp_get_p(E);/*fall through*/
      case t_ELL_Q:  vu = Q_pval(u, p); break;
      case t_ELL_NF: vu = nfval(ellnf_get_nf(E), u, p); break;
      default: pari_err_TYPE("ellgroup", E); vu = 0;
    }
    if (vu) pari_err_TYPE("ellgroup [not a p-minimal curve]",E);
    if (!equali1(kod)) /* bad reduction */
    {
      GEN Ep, T = NULL, q = p, ap = ellap(E,p);
      if (typ(p) == t_INT)
      {
        long i;
        Ep = obj_init(15, 4);
        for (i = 1; i <= 12; i++) gel(Ep,i) = gel(E,i);
      }
      else
      {
        q = pr_norm(p);
        Ep = initsmall5(ellnf_to_Fq(ellnf_get_nf(E), E, p, &p, &T), 4);
      }
      E = FF_ellinit(Ep, Tp_to_FF(T, p)); /* singular curve */
      gel(E,14) = mkvecsmall(t_ELL_Fq);
      obj_insert(E, FF_CARD, subii(q, ap));
    }
    else
      E = ellinit(E, p, 0);
    freeE = 1;
  }
  G = mkvec3(ellff_get_card(E), ellff_get_group(E), ellff_get_gens(E));
  if (!freeE) return gerepilecopy(av, G);
  G = gcopy(G); obj_free(E); return gerepileupto(av, G);
}

GEN
ellgenerators(GEN E)
{
  checkell(E);
  switch(ell_get_type(E))
  {
    case t_ELL_Q:
      return obj_checkbuild(E, Q_GROUPGEN, &elldatagenerators);
    case t_ELL_Fp: case t_ELL_Fq:
      return gcopy(ellff_get_gens(E));
    default:
      pari_err_TYPE("ellgenerators",E);
      return NULL;/*LCOV_EXCL_LINE*/
  }
}

/* char != 2,3, j != 0, 1728 */
static GEN
ellfromj_simple(GEN j)
{
  pari_sp av = avma;
  GEN k = gsubsg(1728,j), kj = gmul(k, j), k2j = gmul(kj, k);
  GEN E = zerovec(5);
  gel(E,4) = gmulsg(3,kj);
  gel(E,5) = gmulsg(2,k2j); return gerepileupto(av, E);
}
GEN
ellfromj(GEN j)
{
  GEN T = NULL, p = typ(j)==t_FFELT? FF_p_i(j): NULL;
  /* trick: use j^0 to get 1 in the proper base field */
  if ((p || (Rg_is_FpXQ(j,&T,&p) && p)) && lgefint(p) == 3) switch(p[2])
  {
    case 2:
      if (gequal0(j))
        retmkvec5(gen_0,gen_0, gpowgs(j,0), gen_0,gen_0);
      else
        retmkvec5(gpowgs(j,0),gen_0,gen_0, gen_0,ginv(j));
    case 3:
      if (gequal0(j))
        retmkvec5(gen_0,gen_0,gen_0, gpowgs(j,0), gen_0);
      else
      {
        GEN E = zerovec(5);
        pari_sp av = avma;
        gel(E,5) = gerepileupto(av, gneg(gsqr(j)));
        gel(E,2) = gcopy(j);
        return E;
      }
  }
  if (gequal0(j)) retmkvec5(gen_0,gen_0,gen_0,gen_0, gpowgs(j,0));
  if (gequalgs(j,1728)) retmkvec5(gen_0,gen_0,gen_0, gpowgs(j,0), gen_0);
  return ellfromj_simple(j);
}

/********************************************************************/
/**                                                                **/
/**                       IS SUPERSINGULAR                         **/
/**                                                                **/
/********************************************************************/

int
elljissupersingular(GEN x)
{
  pari_sp av = avma;
  int res;

  if (typ(x) == t_INTMOD) {
    GEN p = gel(x, 1);
    GEN j = gel(x, 2);
    res = Fp_elljissupersingular(j, p);
  } else if (typ(x) == t_FFELT) {
    GEN j = FF_to_FpXQ_i(x);
    GEN p = FF_p_i(x);
    GEN T = FF_mod(x);
    res = FpXQ_elljissupersingular(j, T, p);
  } else {
    pari_err_TYPE("elljissupersingular", x);
    return 0; /*LCOV_EXCL_LINE*/
  }
  set_avma(av);
  return res;
}

int
ellissupersingular(GEN E, GEN p)
{
  pari_sp av;
  GEN j;
  if (typ(E)!=t_VEC && !p) return elljissupersingular(E);
  p = checkellp(&E, p, NULL, "ellissupersingular");
  j = ell_get_j(E);
  switch(ell_get_type(E))
  {
  case t_ELL_Fp:
  case t_ELL_Fq:
    return elljissupersingular(j);
  case t_ELL_Qp:
  case t_ELL_Q:
    if (typ(j)==t_FRAC && dvdii(gel(j,2), p)) return 0;
    av = avma;
    return gc_bool(av, Fp_elljissupersingular(Rg_to_Fp(j,p), p));
  case t_ELL_NF:
    {
      GEN modP, T, nf = ellnf_get_nf(E), pr = p;
      int res;
      av = avma;
      j = nf_to_scalar_or_basis(nf, j);
      if (dvdii(Q_denom(j), pr_get_p(pr)))
      {
        if (typ(j) == t_FRAC || nfval(nf, j, pr) < 0) return 0;
        modP = nf_to_Fq_init(nf,&pr,&T,&p);
      }
      else
        modP = zk_to_Fq_init(nf,&pr,&T,&p);
      j = nf_to_Fq(nf, j, modP);
      if (typ(j) == t_INT)
        res = Fp_elljissupersingular(j, p);
      else
        res = FpXQ_elljissupersingular(j, T, p);
      return gc_bool(av, res);
    }
  default:
    pari_err_TYPE("ellissupersingular",E);
  }
  return 0; /*LCOV_EXCL_LINE*/
}

GEN
ellsupersingularj(GEN a)
{
  pari_sp av = avma;
  GEN r, T, p;
  long d;
  switch(typ(a))
  {
    case t_INT:
      p = a;
      if (Z_issquare(p)) pari_err_PRIME("ellsupersingularj", p);
      T = init_Fq(p, 2, fetch_user_var("w"));
      d = 2;
      break;
    case t_FFELT:
      p = FF_p_i(a); T = FF_mod(a); d = degpol(T);
      if (!odd(d))
      {
        if (d != 2)
          T = init_Fq(p, 2, varn(T));
        break;
      }
    default: /* FALL THROUGH */
      pari_err_TYPE("ellsupersingular", a);
      return NULL; /* LCOV_EXCL_LINE */
  }
  r = Fq_to_FF(ellsupersingularj_FpXQ(T, p), Tp_to_FF(T, p));
  if (d != 2)
    r = ffmap(ffembed(r, a), r);
  return gerepilecopy(av, r);
}

/* n <= 4, N is the characteristic of the base ring or NULL (char 0) */
static GEN
elldivpol4(GEN e, GEN N, long n, long v)
{
  GEN b2,b4,b6,b8, res;
  if (n==0) return pol_0(v);
  if (n<=2) return N? scalarpol_shallow(mkintmod(gen_1,N),v): pol_1(v);
  b2 = ell_get_b2(e); b4 = ell_get_b4(e);
  b6 = ell_get_b6(e); b8 = ell_get_b8(e);
  if (n==3)
    res = mkpoln(5, N? modsi(3,N): utoi(3),b2,gmulsg(3,b4),gmulsg(3,b6),b8);
  else
  {
    GEN b10 = gsub(gmul(b2, b8), gmul(b4, b6));
    GEN b12 = gsub(gmul(b8, b4), gsqr(b6));
    res = mkpoln(7, N? modsi(2,N): gen_2,b2,gmulsg(5,b4),gmulsg(10,b6),gmulsg(10,b8),b10,b12);
  }
  setvarn(res, v); return res;
}

/* T = (2y + a1x + a3)^4 modulo the curve equation. Store elldivpol(e,n,v)
 * in t[n]. N is the caracteristic of the base ring or NULL (char 0) */
static GEN
elldivpol0(GEN e, GEN t, GEN N, GEN T, long n, long v)
{
  GEN ret;
  long m = n/2;
  if (gel(t,n)) return gel(t,n);
  if (n<=4) ret = elldivpol4(e, N, n, v);
  else if (odd(n))
  {
    GEN t1 = RgX_mul(elldivpol0(e,t,N,T,m+2,v),
                     gpowgs(elldivpol0(e,t,N,T,m,v),3));
    GEN t2 = RgX_mul(elldivpol0(e,t,N,T,m-1,v),
                     gpowgs(elldivpol0(e,t,N,T,m+1,v),3));
    if (odd(m))/*f_{4l+3} = f_{2l+3}f_{2l+1}^3 - T f_{2l}f_{2l+2}^3, m=2l+1*/
      ret = RgX_sub(t1, RgX_mul(T,t2));
    else       /*f_{4l+1} = T f_{2l+2}f_{2l}^3 - f_{2l-1}f_{2l+1}^3, m=2l*/
      ret = RgX_sub(RgX_mul(T,t1), t2);
  }
  else
  { /* f_2m = f_m(f_{m+2}f_{m-1}^2 - f_{m-2}f_{m+1}^2) */
    GEN t1 = RgX_mul(elldivpol0(e,t,N,T,m+2,v),
                     RgX_sqr(elldivpol0(e,t,N,T,m-1,v)));
    GEN t2 = RgX_mul(elldivpol0(e,t,N,T,m-2,v),
                     RgX_sqr(elldivpol0(e,t,N,T,m+1,v)));
    ret = RgX_mul(elldivpol0(e,t,N,T,m,v), RgX_sub(t1,t2));
  }
  gel(t,n) = ret;
  return ret;
}

GEN
elldivpol(GEN e, long n0, long v)
{
  pari_sp av = avma;
  GEN f, D, N;
  long n = labs(n0);

  checkell(e); D = ell_get_disc(e);
  if (v < 0) v = 0;
  if (varncmp(gvar(D), v) <= 0) pari_err_PRIORITY("elldivpol", e, "<=", v);
  N = characteristic(D); if (!signe(N)) N = NULL;
  if (n==1 || n==3)
    f = elldivpol4(e, N, n, v);
  else
  {
    GEN d2 = ec_bmodel(e,v); /* (2y + a1x + a3)^2 mod E */
    if (N && !mod2(N)) { gel(d2,5) = modsi(4,N); d2 = normalizepol(d2); }
    if (n <= 4)
      f = elldivpol4(e, N, n, v);
    else
      f = elldivpol0(e, const_vec(n,NULL), N,RgX_sqr(d2), n, v);
    if (n%2==0) f = RgX_mul(f, d2);
  }
  if (n0 < 0) return gerepileupto(av, RgX_neg(f));
  return gerepilecopy(av, f);
}

/* return [phi_n, (psi_n)^2] such that x[nP] = phi_n / (psi_n)^2 */
GEN
ellxn(GEN e, long n, long v)
{
  pari_sp av = avma;
  GEN d2, D, N, A, B;
  checkell(e); D = ell_get_disc(e);
  if (v==-1) v = 0;
  if (varncmp(gvar(D), v) <= 0) pari_err_PRIORITY("elldivpol", e, "<=", v);
  N = characteristic(D);
  if (!signe(N)) N = NULL;
  if (n < 0) n = -n;
  d2 = ec_bmodel(e,v); /* (2y + a1x + 3)^2 mod E */
  if (N && !mod2(N)) { gel(d2,5) = modsi(4,N); d2 = normalizepol(d2); }
  if (n == 0)
  {
    A = pol_0(v);
    B = pol_0(v);
  }
  else if (n == 1)
  {
    A = pol_1(v);
    B = pol_x(v);
  }
  else if (n == 2)
  {
    A = d2;
    B = ec_phi2(e, v);
  }
  else
  {
    GEN t = const_vec(n+1,NULL), T = RgX_sqr(d2);
    GEN f = elldivpol0(e, t, N, T, n, v); /* f_n / d2^(n odd)*/
    GEN g = elldivpol0(e, t, N, T, n-1, v); /* f_{n-1} / d2^(n even) */
    GEN h = elldivpol0(e, t, N, T, n+1, v); /* f_{n+1} / d2^(n even) */
    GEN f2 = RgX_sqr(f), u = RgX_mul(g,h);
    if (!odd(n))
      A = RgX_mul(f2, d2);
    else
    { A = f2; u = RgX_mul(u,d2); }
    /* A = psi_n^2, u = psi_{n-1} psi_{n+1} */
    B = RgX_sub(RgX_shift(A,1), u);
  }
  return gerepilecopy(av, mkvec2(B,A));
}

/* l and p primes; p = 1 mod l; return an element of order l in (Z/pZ)^* */
static ulong
ltors_Fl(ulong l, ulong p)
{
  ulong x, y, r = (p-1)/l;
  for (x = 2;; x++) { y = Fl_powu(x, r, p); if (y != 1) return y; }
}

/* Assume that l|o but p!=1 [l] so r_l E(F_p) = 1 */
static void
FljV_vecsat_Siksek(GEN E, GEN P, ulong o, ulong l, ulong a4, ulong a6, ulong p,
                   GEN S, long *m)
{
  long i, n = lg(P)-1;
  GEN a4a6, g, F, v = zero_zv(n);
  pari_sp av = avma;
  ulong q = o / l;

  F = mkmat2(mkcols(l), mkcols(1));
  a4a6 = a4a6_ch_Fl(E, p);
  g = gel(Fl_ellptors(l, o, a4, a6, p), 1);
  for (i=1; i <= n; i++)
  {
    GEN Q = Fle_changepointinv(Flj_to_Fle(gel(P,i), p), a4a6, p);
    if (!ell_is_inf(Q))
      v[i] = itou(Fle_log(Fle_mulu(Q, q, a4, p), g, F, a4, p));
  }
  gel(S,(*m)++) = v;
  set_avma(av);
}

/* Assume that l|o and p=1 [l] so r_l E(F_p) = 1 or 2 */
static void
FljV_vecsat_Prickett(GEN E, GEN P, ulong o, ulong l, ulong a4, ulong a6,
                     ulong p, GEN S, long *m)
{
  long i, n = lg(P)-1;
  GEN a4a6, G, G1, G2, v = zero_zv(n), w = zero_zv(n);
  ulong g = ltors_Fl(l, p), q = (p-1)/l;
  pari_sp av = avma;

  a4a6 = a4a6_ch_Fl(E, p);
  G = Fl_ellptors(l, o, a4, a6, p);
  G1 = gel(G,1);
  G2 = lg(G)==3 ? gel(G, 2): NULL;
  for (i = 1; i <= n; i++)
  {
    GEN Q = Fle_changepointinv(Flj_to_Fle(gel(P,i), p), a4a6, p);
    if (!ell_is_inf(Q))
    {
      ulong u = Fl_powu(Fle_tatepairing(G1, Q, l, a4, p), q, p);
      v[i] = Fl_log(u, g, l, p);
      if (G2)
      {
        ulong u = Fl_powu(Fle_tatepairing(G2, Q, l, a4, p), q, p);
        w[i] = Fl_log(u, g, l, p);
      }
    }
  }
  gel(S,(*m)++) = v;
  if (G2 && *m < lg(S)) gel(S,(*m)++) = w;
  set_avma(av);
}

static void
FljV_vecsat(GEN E, GEN P, ulong o, ulong l, ulong a4, ulong a6, ulong p,
            GEN S, long *m)
{
  P = ZM_to_Flm(P, p);
  if (p % l == 1)
    FljV_vecsat_Prickett(E, P, o, l, a4, a6, p, S, m);
  else
    FljV_vecsat_Siksek(E, P, o, l, a4, a6, p, S, m);
}

/* P a vector of points in E(Q), return a linear map M from the abelian group
 * they generate to Z/lZ; sum x[i] P[i] is l-divisible => x M = 0 */
static GEN
ellsatp_mat(hashtable *h, GEN E, long CM, GEN P, ulong l, long nb)
{
  long m = 1;
  GEN D = ell_get_disc(E), M = cgetg(nb+1, t_MAT);
  forprime_t S;

  P = QEV_to_ZJV(P);
  (void)u_forprime_init(&S, 5, ULONG_MAX);
  while (m <= nb)
  {
    ulong a4, a6, p = u_forprime_next(&S);
    long o;
    if (dvdiu(D, p)) continue;
    Fl_ell_to_a4a6(E, p, &a4, &a6);
    if (!hash_haskey_long(h, (void*)p, &o))
    {
      o = p+1 - Fl_elltrace_CM(CM, a4, a6, p);
      hash_insert_long(h,(void*)p, o);
    }
    if (o % l == 0) FljV_vecsat(E, P, o, l, a4, a6, p, M, &m);
  }
  return M;
}

INLINE long
Flv_firstnonzero(GEN v)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++)
    if (v[i]) break;
  return i;
}

/* update M in place */
static GEN
ellsatp(hashtable *hh, GEN E, long CM, GEN T, GEN H, GEN M, ulong l, GEN *xl,
        long vxl, long nb, long prec)
{
  GEN P = T ? shallowconcat(H, T): H;
  GEN S = ellsatp_mat(hh, E, CM, P, l, nb); /* fill hh */
  pari_sp av = avma;
  GEN K = Flm_ker(Flm_transpose(S), l);
  long i, lK = lg(K), nH = lg(H)-1;

  if (lK==1) return gc_NULL(av);
  if (DEBUGLEVEL >= 3)
    err_printf("ellsat: potential factor %lu, dim Ker = %ld\n",l,lK-1);
  /* Mazur bound for torsion of isogenous curves */
  if (!*xl && l <= 7) *xl = ellxn(E, l, vxl);
  for (i = 1; i < lK; i++)
  {
    GEN ki = gel(K,i), Ki, h, R;
    long f = Flv_firstnonzero(ki);

    /* for T != NULL: avoid solving for [p]Q = R when R is p-torsion */
    if (f > nH) continue;
    if (ki[f] != 1) ki = Flv_Fl_div(ki, ki[f], l);
    Ki = zv_to_ZV(Flv_center(ki, l, l >> 1));
    h = qfeval(M, T? vecslice(Ki, 1, nH): Ki);
    if (*xl)
    {
      GEN Q = ellQ_factorback(E, P, Ki, 1, h, prec);
      if (ellisdivisible(E, Q, *xl, &R)) h = gdiv(h, sqru(l)); else R = NULL;
    }
    else
    {
      h = gdiv(h, sqru(l));
      R = ellQ_factorback(E, P, Ki, l, h, prec);
    }
    if (DEBUGLEVEL >= 2)
      err_printf("ellsat: %s divisible by %lu\n", R? "": "not", l);
    if (!R)
    {
      if (lK == 2) break;
      return l > 7? gc_const(av,H): H; /* fail: return and retry */
    }
    gcoeff(M, f, f) = h;
    for (i = 1; i <= nH; i++)
      if (i != f) gcoeff(M, f, i) = gdivgu(RgV_dotproduct(gel(M,i), Ki), l);
    for (i = 1; i <= nH; i++) gcoeff(M, i, f) = gcoeff(M, f, i);
    gel(H,f) = R; return H; /* found l-divisible point: return new lattice */
  }
  return gc_NULL(av); /* l-saturated */
}

static GEN
ellQ_saturation(GEN E, GEN P, long B, long prec)
{
  forprime_t S;
  GEN M = ellheightmatrix(E, P, prec);
  long CM = ellQ_get_CM(E), w = fetch_var_higher();
  hashtable h;
  ulong p;

  hash_init_ulong(&h, 16, 1);
  (void)u_forprime_init(&S, 2, B);
  P = leafcopy(P); /* modified in place by ellsatp */
  while((p = u_forprime_next(&S)))
  {
    long nb = lg(P)-1 + 25 / log2(p) - 1; /* error ~ 2^{-25} */
    GEN xp = NULL, T = gel(elltors_psylow(E, p), 3);
    if (lg(T)==1) T = NULL;
    while (1)
    {
      GEN Q = ellsatp(&h, E, CM, T, P, M, p, &xp, w, nb, prec);
      if (!Q) break;
      nb += lg(P)-1;
      P = Q;
    }
  }
  return ellQ_genreduce(E, P, M, prec);
}

GEN
ellsaturation(GEN E, GEN P, long B, long prec)
{
  pari_sp av = avma;
  GEN urst;

  if (lg(P) == 1) return cgetg(1, t_VEC);
  E = ellminimalmodel(E, &urst);
  if (is_trivial_change(urst)) urst = NULL;
  if (urst) P = ellchangepoint(P, urst);
  P = ellQ_saturation(E, P, B, prec);
  if (urst) P = ellchangepoint(P, ellchangeinvert(urst));
  obj_free(E); return gerepilecopy(av, P);
}

static GEN
to_RgX(GEN P, long vx)
{ return typ(P)==t_POL && varn(P) == vx ? P: scalarpol_shallow(P, vx); }
GEN
elltrace(GEN E, GEN P)
{
  pari_sp av = avma;
  GEN xP,yP, T, lxP, Q, LP, M, K, U,V,R, xQ,yQ;
  long v, n, i, j, d;

  checkell(E);
  checkellpt(P);
  if (ell_is_inf(P)) return gcopy(P); /* P == oo */
  /* More checks */
  xP = gel(P,1); if (typ(xP)!=t_POLMOD) pari_err_TYPE("elltrace",xP);
  yP = gel(P,2); if (typ(yP)!=t_POLMOD) pari_err_TYPE("elltrace",yP);
  T = gel(xP,1); v = varn(T); n = degpol(T);
  if (!gequal(gel(yP,1),T)) pari_err_MODULUS("elltrace",xP,yP);
  /* Trivial cases */
  if (n == 1) { set_avma(av); return gcopy(P); }
  lxP = to_RgX(gel(xP,2), v);
  if (!degpol(lxP))
  {
    GEN lyP = to_RgX(gel(yP,2), v);
    if (degpol(lyP)){ set_avma(av); retmkvec(gen_0); }
    P = mkvec2(gel(lxP,2), gel(lyP,2));
    return gerepileupto(av, ellmul(E, P, utoipos(n)));
  }
  /* Strategy: look for a function with divisor equal to
   *   [P_1] + ... + [P_n] + [-Tr(P)] - (n+1)[0]. */
  LP = cgetg(n+2,t_VEC); /* basis of the Riemann-Roch space evaluated at P */
  gel(LP,1) = gen_1;
  gel(LP,2) = xP;
  gel(LP,3) = yP;
  for (i = 4; i <= n+1; i++) gel(LP,i) = gmul(gel(LP,i-2), xP);
  M = cgetg(n+2,t_MAT); /* functions defined over K vanishing at P */
  for (j = 1; j <= n+1; j++)
  {
    GEN t = to_RgX(liftpol_shallow(gel(LP,j)), v);
    for (i = 1; i <= n; i++) gel(M,j) = RgX_to_RgC(t, n);
  }
  K = gel(ker(M),1);
  /* Coords on 1,x,y,x^2,xy,.. of function f of smallest degree vanishing at P
   * div f = [P_1] + ... + [P_d] + [-Tr(P)] - (d+1)[0]
   * with deg(K(P)) = d+1 if Tr(P) != 0; = d otherwise; f = U(x) + y*V(x) */
  U = cgetg((n+1)/2+3,t_POL);
  V = cgetg((n-2)/2+3,t_POL); U[1] = V[1] = evalvarn(0);
  gel(U,2) = gel(K,1); /* Coef of 1 */
  for(i = 1; 2*i <= n+1; i++) gel(U,i+2) = gel(K,2*i); /* Coef of x^i */
  for(i = 0; 2*i+3 <= n+1; i++) gel(V,i+2) = gel(K,2*i+3); /* Coef of x^i*y */
  U = normalizepol(U); V = normalizepol(V);
  /* f does not depend on y, so trace = oo */
  if (signe(V)==0) { set_avma(av); retmkvec(gen_0); }
  /* Plug y = -U(x)/V(x) into Weierstrass equation:
   * 0 = ((x^3+a2x^2+a4x+a6)*V + (a1x+a3)*U)*V - U^2 */
  R = mkpoln(4, gen_1, ell_get_a2(E), ell_get_a4(E), ell_get_a6(E));
  R = gmul(R, V);
  R = gadd(R, gmul(U, mkpoln(2,ell_get_a1(E),ell_get_a3(E))));
  R = gmul(R, V);
  R = gsub(R, gsqr(U));
  /* Discard Galois orbit of P */
  R = RgX_div(R, minpoly(xP,0));
  /* What is left is either constant -> return oo, or deg 1 -> nontrivial trace. */
  if(degpol(R)==0) { set_avma(av); retmkvec(gen_0); }
  /* Recover the trace */
  xQ = gneg(gdiv(gel(R,2), gel(R,3)));
  yQ = gneg(gdiv(poleval(U, xQ), poleval(V, xQ)));
  Q = mkvec2(xQ, yQ);
  /* So far, we have computed -Tr(P) over the extension K(P)/K
   * we still need to compute [L:K(P)] */
  d = 0;
  for (i = n+1; i > 0; i--) if (!gequal0(gel(K,i))) { d = i; break; }
  return gerepileupto(av, ellmul(E, Q, stoi(-n / (d-1))));
}
