/* Copyright (C) 2012  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_ellcard

/* Not so fast arithmetic with points over elliptic curves over Fq,
small characteristic. */

/***********************************************************************/
/**                                                                   **/
/**                              FlxqE                                **/
/**                                                                   **/
/***********************************************************************/
/* These functions deal with point over elliptic curves over Fq defined
 * by an equation of the form y^2=x^3+a4*x+a6. Most of the time a6 is omitted
 * since it can be recovered from any point on the curve. */

GEN
RgE_to_FlxqE(GEN x, GEN T, ulong p)
{
  if (ell_is_inf(x)) return x;
  retmkvec2(Rg_to_Flxq(gel(x,1),T,p), Rg_to_Flxq(gel(x,2),T,p));
}

GEN
FlxqE_changepoint(GEN x, GEN ch, GEN T, ulong p)
{
  pari_sp av = avma;
  GEN p1, p2, z, u, r, s, t, v, v2, v3;
  ulong pi;
  if (ell_is_inf(x)) return x;
  pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  u = gel(ch,1); r = gel(ch,2);
  s = gel(ch,3); t = gel(ch,4);
  v = Flxq_inv_pre(u, T, p, pi);
  v2 = Flxq_sqr_pre(v, T, p, pi);
  v3 = Flxq_mul_pre(v,v2, T, p, pi);
  p1 = Flx_sub(gel(x,1), r, p);
  p2 = Flx_sub(gel(x,2), Flx_add(Flxq_mul_pre(s, p1, T, p, pi),t, p), p);
  z = cgetg(3,t_VEC);
  gel(z,1) = Flxq_mul_pre(v2, p1, T, p, pi);
  gel(z,2) = Flxq_mul_pre(v3, p2, T, p, pi);
  return gerepileupto(av, z);
}

GEN
FlxqE_changepointinv(GEN x, GEN ch, GEN T, ulong p)
{
  pari_sp av = avma;
  GEN p1, p2, u, r, s, t, X, Y, u2, u3, u2X, z;
  ulong pi;
  if (ell_is_inf(x)) return x;
  pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  X = gel(x,1); Y = gel(x,2);
  u = gel(ch,1); r = gel(ch,2);
  s = gel(ch,3); t = gel(ch,4);
  u2 = Flxq_sqr_pre(u, T, p, pi);
  u3 = Flxq_mul_pre(u,u2, T, p, pi);
  u2X = Flxq_mul_pre(u2,X, T, p, pi);
  p1 = Flxq_mul_pre(u3,Y, T, p, pi);
  p2 = Flx_add(Flxq_mul_pre(s,u2X, T, p, pi), t, p);
  z = cgetg(3, t_VEC);
  gel(z,1) = Flx_add(u2X, r, p);
  gel(z,2) = Flx_add(p1, p2, p);
  return gerepileupto(av, z);
}

static GEN
nonsquare_Flxq(GEN T, ulong p)
{
  pari_sp av = avma;
  long n = degpol(T), vs = T[1];
  GEN a;
  if (odd(n))
    return mkvecsmall2(vs, nonsquare_Fl(p));
  do
  {
    set_avma(av);
    a = random_Flx(n, vs, p);
  } while (Flxq_issquare(a, T, p));
  return a;
}

void
Flxq_elltwist(GEN a, GEN a6, GEN T, ulong p, GEN *pt_a, GEN *pt_a6)
{
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  GEN d = nonsquare_Flxq(T, p);
  GEN d2 = Flxq_sqr_pre(d, T, p, pi), d3 = Flxq_mul_pre(d2, d, T, p, pi);
  if (typ(a)==t_VECSMALL)
  {
    *pt_a  = Flxq_mul_pre(a,  d2, T, p, pi);
    *pt_a6 = Flxq_mul_pre(a6, d3, T, p, pi);
  } else
  {
    *pt_a  = mkvec(Flxq_mul_pre(gel(a,1), d, T, p, pi));
    *pt_a6 = Flxq_mul_pre(a6, d3, T, p, pi);
  }
}

static GEN
FlxqE_dbl_slope(GEN P, GEN a4, GEN T, ulong p, ulong pi, GEN *ps)
{
  GEN x, y, Q, s;
  if (ell_is_inf(P) || !lgpol(gel(P,2))) return ellinf();
  x = gel(P,1); y = gel(P,2);
  if (p==3UL)
    s = typ(a4)==t_VEC? Flxq_div_pre(Flxq_mul_pre(x, gel(a4,1), T,p,pi), y, T,p,pi)
                      : Flxq_div_pre(a4, Flx_neg(y, p), T,p,pi);
  else
  {
    GEN sx = Flx_add(Flx_triple(Flxq_sqr_pre(x, T, p, pi), p), a4, p);
    s = Flxq_div_pre(sx, Flx_double(y, p), T, p, pi);
  }
  Q = cgetg(3,t_VEC);
  gel(Q,1) = Flx_sub(Flxq_sqr_pre(s, T, p, pi), Flx_double(x, p), p);
  if (typ(a4)==t_VEC) gel(Q, 1) = Flx_sub(gel(Q,1), gel(a4,1), p);
  gel(Q,2) = Flx_sub(Flxq_mul_pre(s, Flx_sub(x, gel(Q,1), p), T, p, pi),
                     y, p);
  if (ps) *ps = s;
  return Q;
}

GEN
FlxqE_dbl(GEN P, GEN a4, GEN T, ulong p)
{
  pari_sp av = avma;
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return gerepileupto(av, FlxqE_dbl_slope(P,a4, T, p, pi, NULL));
}

static GEN
FlxqE_add_slope(GEN P, GEN Q, GEN a4, GEN T, ulong p, ulong pi, GEN *ps)
{
  GEN Px, Py, Qx, Qy, R, s;
  if (ell_is_inf(P)) return Q;
  if (ell_is_inf(Q)) return P;
  Px = gel(P,1); Py = gel(P,2);
  Qx = gel(Q,1); Qy = gel(Q,2);
  if (Flx_equal(Px, Qx))
  {
    if (Flx_equal(Py, Qy))
      return FlxqE_dbl_slope(P, a4, T, p, pi, ps);
    else
      return ellinf();
  }
  s = Flxq_div_pre(Flx_sub(Py, Qy, p), Flx_sub(Px, Qx, p), T, p, pi);
  R = cgetg(3,t_VEC);
  gel(R,1) = Flx_sub(Flx_sub(Flxq_sqr_pre(s, T, p, pi), Px, p), Qx, p);
  if (typ(a4)==t_VEC) gel(R,1) = Flx_sub(gel(R,1), gel(a4,1), p);
  gel(R,2) = Flx_sub(Flxq_mul_pre(s, Flx_sub(Px, gel(R,1), p), T, p, pi),
                     Py, p);
  if (ps) *ps = s;
  return R;
}

GEN
FlxqE_add(GEN P, GEN Q, GEN a4, GEN T, ulong p)
{
  pari_sp av = avma;
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return gerepileupto(av, FlxqE_add_slope(P,Q,a4, T,p,pi, NULL));
}

static GEN
FlxqE_neg_i(GEN P, ulong p)
{
  if (ell_is_inf(P)) return P;
  return mkvec2(gel(P,1), Flx_neg(gel(P,2), p));
}

GEN
FlxqE_neg(GEN P, GEN T, ulong p)
{
  (void) T;
  if (ell_is_inf(P)) return ellinf();
  return mkvec2(gcopy(gel(P,1)), Flx_neg(gel(P,2), p));
}

GEN
FlxqE_sub(GEN P, GEN Q, GEN a4, GEN T, ulong p)
{
  pari_sp av = avma;
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return gerepileupto(av, FlxqE_add_slope(P, FlxqE_neg_i(Q, p), a4, T,p,pi, NULL));
}

struct _FlxqE
{
  GEN a4, a6, T;
  ulong p, pi;
};

static GEN
_FlxqE_dbl(void *E, GEN P)
{
  struct _FlxqE *e = (struct _FlxqE *) E;
  return FlxqE_dbl_slope(P, e->a4, e->T, e->p, e->pi, NULL);
}

static GEN
_FlxqE_add(void *E, GEN P, GEN Q)
{
  struct _FlxqE *e = (struct _FlxqE *) E;
  return FlxqE_add_slope(P, Q, e->a4, e->T, e->p, e->pi, NULL);
}

static GEN
_FlxqE_sub(void *E, GEN P, GEN Q)
{
  struct _FlxqE *e = (struct _FlxqE *) E;
  return FlxqE_add_slope(P, FlxqE_neg_i(Q,e->p), e->a4, e->T,e->p,e->pi, NULL);
}

static GEN
_FlxqE_mul(void *E, GEN P, GEN n)
{
  pari_sp av = avma;
  struct _FlxqE *e=(struct _FlxqE *) E;
  long s = signe(n);
  if (!s || ell_is_inf(P)) return ellinf();
  if (s < 0) P = FlxqE_neg(P, e->T, e->p);
  if (is_pm1(n)) return s>0? gcopy(P): P;
  return gerepilecopy(av, gen_pow_i(P, n, e, &_FlxqE_dbl, &_FlxqE_add));
}

GEN
FlxqE_mul(GEN P, GEN n, GEN a4, GEN T, ulong p)
{
  struct _FlxqE E;
  E.a4= a4; E.T = T; E.p = p; E.pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return _FlxqE_mul(&E, P, n);
}

/* 3*x^2+2*a2*x = -a2*x, and a2!=0 */

/* Finds a random nonsingular point on E */
static GEN
random_F3xqE(GEN a2, GEN a6, GEN T)
{
  pari_sp ltop = avma;
  GEN x, y, rhs;
  const ulong p = 3;
  do
  {
    set_avma(ltop);
    x   = random_Flx(get_Flx_degree(T),get_Flx_var(T),p);
    rhs = Flx_add(Flxq_mul(Flxq_sqr(x, T, p), Flx_add(x, a2, p), T, p), a6, p);
  } while ((!lgpol(rhs) && !lgpol(x)) || !Flxq_issquare(rhs, T, p));
  y = Flxq_sqrt(rhs, T, p);
  if (!y) pari_err_PRIME("random_F3xqE", T);
  return gerepilecopy(ltop, mkvec2(x, y));
}

/* Finds a random nonsingular point on E */
static GEN
random_FlxqE_pre(GEN a4, GEN a6, GEN T, ulong p, ulong pi)
{
  pari_sp ltop = avma;
  GEN x, x2, y, rhs;
  if (typ(a4)==t_VEC) return random_F3xqE(gel(a4,1), a6, T);
  do
  {
    set_avma(ltop);
    x   = random_Flx(get_Flx_degree(T),get_Flx_var(T),p);
    x2  = Flxq_sqr_pre(x, T, p, pi); /*  x^3+a4*x+a6 = x*(x^2+a4)+a6  */
    rhs = Flx_add(Flxq_mul_pre(x, Flx_add(x2, a4, p), T, p, pi), a6, p);
  } while ((!lgpol(rhs) && !lgpol(Flx_add(Flx_triple(x2, p), a4, p)))
          || !Flxq_issquare(rhs, T, p));
  y = Flxq_sqrt(rhs, T, p);
  if (!y) pari_err_PRIME("random_FlxqE", T);
  return gerepilecopy(ltop, mkvec2(x, y));
}
GEN
random_FlxqE(GEN a4, GEN a6, GEN T, ulong p)
{ return random_FlxqE_pre(a4, a6, T, p, SMALL_ULONG(p)? 0: get_Fl_red(p)); }

static GEN
_FlxqE_rand(void *E)
{
  struct _FlxqE *e=(struct _FlxqE *) E;
  return random_FlxqE_pre(e->a4, e->a6, e->T, e->p, e->pi);
}

static const struct bb_group FlxqE_group={_FlxqE_add,_FlxqE_mul,_FlxqE_rand,hash_GEN,zvV_equal,ell_is_inf, NULL};

const struct bb_group *
get_FlxqE_group(void ** pt_E, GEN a4, GEN a6, GEN T, ulong p)
{
  struct _FlxqE *e = (struct _FlxqE *) stack_malloc(sizeof(struct _FlxqE));
  e->a4 = a4; e->a6 = a6;
  e->pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  e->p = p;
  e->T = Flx_get_red_pre(T, p, e->pi);
  *pt_E = (void *)e; return &FlxqE_group;
}

GEN
FlxqE_order(GEN z, GEN o, GEN a4, GEN T, ulong p)
{
  pari_sp av = avma;
  struct _FlxqE e;
  e.a4 = a4; e.T = T; e.p = p; e.pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return gerepileuptoint(av, gen_order(z, o, (void*)&e, &FlxqE_group));
}

GEN
FlxqE_log(GEN a, GEN b, GEN o, GEN a4, GEN T, ulong p)
{
  pari_sp av = avma;
  struct _FlxqE e;
  e.a4 = a4; e.T = T; e.p = p; e.pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return gerepileuptoint(av, gen_PH_log(a, b, o, (void*)&e, &FlxqE_group));
}

/***********************************************************************/
/**                            Pairings                               **/
/***********************************************************************/
/* Derived from APIP by Jerome Milan, 2012 */
static GEN
FlxqE_vert(GEN P, GEN Q, GEN a4, GEN T, ulong p, ulong pi)
{
  long vT = get_Flx_var(T);
  GEN df;
  if (ell_is_inf(P)) return pol1_Flx(vT);
  if (!Flx_equal(gel(Q,1), gel(P,1))) return Flx_sub(gel(Q,1), gel(P,1), p);
  if (lgpol(gel(P,2))!=0) return pol1_Flx(vT);
  df = typ(a4)==t_VEC ? Flxq_mul_pre(gel(P,1), Flx_double(gel(a4,1), p), T,p,pi)
                      : a4;
  return Flxq_inv_pre(Flx_add(Flx_triple(Flxq_sqr_pre(gel(P,1), T,p, pi), p),
                              df, p), T, p, pi);
}

static GEN
FlxqE_Miller_line(GEN R, GEN Q, GEN slope, GEN a4, GEN T, ulong p, ulong pi)
{
  long vT = get_Flx_var(T);
  GEN x = gel(Q,1), y = gel(Q,2);
  GEN tmp1 = Flx_sub(x, gel(R,1), p);
  GEN tmp2 = Flx_add(Flxq_mul_pre(tmp1, slope, T, p, pi), gel(R,2), p);
  if (!Flx_equal(y, tmp2)) return Flx_sub(y, tmp2, p);
  if (lgpol(y) == 0) return pol1_Flx(vT);
  else
  {
    GEN s1, s2, a2 = typ(a4)==t_VEC ? gel(a4,1): NULL;
    GEN y2i = Flxq_inv_pre(Flx_mulu(y, 2, p), T, p, pi);
    GEN df = a2 ? Flxq_mul_pre(x, Flx_mulu(a2, 2, p), T, p, pi): a4;
    GEN x3, ddf;
    s1 = Flxq_mul_pre(Flx_add(Flx_triple(Flxq_sqr_pre(x, T, p, pi), p), df, p), y2i, T, p, pi);
    if (!Flx_equal(s1, slope)) return Flx_sub(s1, slope, p);
    x3 = Flx_triple(x, p);
    ddf = a2 ? Flx_add(x3, a2, p): x3;
    s2 = Flxq_mul_pre(Flx_sub(ddf, Flxq_sqr_pre(s1, T,p,pi), p), y2i, T,p,pi);
    return lgpol(s2)!=0 ? s2: y2i;
  }
}

/* Computes the equation of the line tangent to R and returns its
 * evaluation at the point Q. Also doubles the point R. */
static GEN
FlxqE_tangent_update(GEN R, GEN Q, GEN a4, GEN T, ulong p, ulong pi, GEN *pt_R)
{
  if (ell_is_inf(R))
  {
    *pt_R = ellinf();
    return pol1_Flx(get_Flx_var(T));
  }
  else if (!lgpol(gel(R,2)))
  {
    *pt_R = ellinf();
    return FlxqE_vert(R, Q, a4, T, p, pi);
  } else {
    GEN slope;
    *pt_R = FlxqE_dbl_slope(R, a4, T, p, pi, &slope);
    return FlxqE_Miller_line(R, Q, slope, a4, T, p, pi);
  }
}

/* Computes the equation of the line through R and P, and returns its
 * evaluation at the point Q. Also adds P to the point R. */
static GEN
FlxqE_chord_update(GEN R, GEN P, GEN Q, GEN a4, GEN T, ulong p, ulong pi, GEN *pt_R)
{
  if (ell_is_inf(R))
  {
    *pt_R = gcopy(P);
    return FlxqE_vert(P, Q, a4, T, p, pi);
  }
  else if (ell_is_inf(P))
  {
    *pt_R = gcopy(R);
    return FlxqE_vert(R, Q, a4, T, p, pi);
  }
  else if (Flx_equal(gel(P, 1), gel(R, 1)))
  {
    if (Flx_equal(gel(P, 2), gel(R, 2)))
      return FlxqE_tangent_update(R, Q, a4, T, p, pi, pt_R);
    else
    {
      *pt_R = ellinf();
      return FlxqE_vert(R, Q, a4, T, p, pi);
    }
  } else {
    GEN slope;
    *pt_R = FlxqE_add_slope(P, R, a4, T, p, pi, &slope);
    return FlxqE_Miller_line(R, Q, slope, a4, T, p, pi);
  }
}

struct _FlxqE_miller
{
  ulong p, pi;
  GEN T, a4, P;
};

static GEN
FlxqE_Miller_dbl(void* E, GEN d)
{
  struct _FlxqE_miller *m = (struct _FlxqE_miller *)E;
  ulong p = m->p, pi = m->pi;
  GEN T = m->T, a4 = m->a4, P = m->P;
  GEN v, line, point = gel(d,3);
  GEN N = Flxq_sqr_pre(gel(d,1), T, p, pi);
  GEN D = Flxq_sqr_pre(gel(d,2), T, p, pi);
  line = FlxqE_tangent_update(point, P, a4, T, p, pi, &point);
  N  = Flxq_mul_pre(N, line, T, p, pi);
  v = FlxqE_vert(point, P, a4, T, p, pi);
  D = Flxq_mul_pre(D, v, T, p, pi); return mkvec3(N, D, point);
}

static GEN
FlxqE_Miller_add(void* E, GEN va, GEN vb)
{
  struct _FlxqE_miller *m = (struct _FlxqE_miller *)E;
  ulong p = m->p, pi = m->pi;
  GEN T = m->T, a4 = m->a4, P = m->P;
  GEN v, line, point;
  GEN na = gel(va,1), da = gel(va,2), pa = gel(va,3);
  GEN nb = gel(vb,1), db = gel(vb,2), pb = gel(vb,3);
  GEN N = Flxq_mul_pre(na, nb, T, p, pi);
  GEN D = Flxq_mul_pre(da, db, T, p, pi);
  line = FlxqE_chord_update(pa, pb, P, a4, T, p, pi, &point);
  N  = Flxq_mul_pre(N, line, T, p, pi);
  v = FlxqE_vert(point, P, a4, T, p, pi);
  D = Flxq_mul_pre(D, v, T, p, pi); return mkvec3(N, D, point);
}

/* Returns the Miller function f_{m, Q} evaluated at the point P using
 * the standard Miller algorithm. */
static GEN
FlxqE_Miller(GEN Q, GEN P, GEN m, GEN a4, GEN T, ulong p, ulong pi)
{
  pari_sp av = avma;
  struct _FlxqE_miller d;
  GEN v, N, D, g1;

  d.a4 = a4; d.T = T; d.p = p; d.P = P; d.pi = pi;
  g1 = pol1_Flx(get_Flx_var(T));
  v = gen_pow_i(mkvec3(g1,g1,Q), m, (void*)&d,
                FlxqE_Miller_dbl, FlxqE_Miller_add);
  N = gel(v,1); D = gel(v,2);
  return gerepileupto(av, Flxq_div_pre(N, D, T, p, pi));
}

GEN
FlxqE_weilpairing_pre(GEN P, GEN Q, GEN m, GEN a4, GEN T, ulong p, ulong pi)
{
  pari_sp av = avma;
  GEN N, D, w;
  if (ell_is_inf(P) || ell_is_inf(Q)
    || (Flx_equal(gel(P,1),gel(Q,1)) && Flx_equal(gel(P,2),gel(Q,2))))
    return pol1_Flx(get_Flx_var(T));
  N = FlxqE_Miller(P, Q, m, a4, T, p, pi);
  D = FlxqE_Miller(Q, P, m, a4, T, p, pi);
  w = Flxq_div_pre(N, D, T, p, pi); if (mpodd(m)) w = Flx_neg(w, p);
  return gerepileupto(av, w);
}
GEN
FlxqE_weilpairing(GEN P, GEN Q, GEN m, GEN a4, GEN T, ulong p)
{ return FlxqE_weilpairing_pre(P,Q,m,a4,T,p, SMALL_ULONG(p)?0:get_Fl_red(p)); }

GEN
FlxqE_tatepairing(GEN P, GEN Q, GEN m, GEN a4, GEN T, ulong p)
{
  if (ell_is_inf(P) || ell_is_inf(Q)) return pol1_Flx(get_Flx_var(T));
  return FlxqE_Miller(P, Q, m, a4, T, p, SMALL_ULONG(p)? 0: get_Fl_red(p));
}

static GEN
_FlxqE_pairorder(void *E, GEN P, GEN Q, GEN m, GEN F)
{
  struct _FlxqE *e = (struct _FlxqE *) E;
  return  Flxq_order(FlxqE_weilpairing_pre(P,Q,m,e->a4,e->T,e->p,e->pi), F, e->T, e->p);
}

GEN
Flxq_ellgroup(GEN a4, GEN a6, GEN N, GEN T, ulong p, GEN *pt_m)
{
  struct _FlxqE e;
  GEN q = powuu(p, get_Flx_degree(T));
  e.a4=a4; e.a6=a6; e.T=T; e.p=p; e.pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  return gen_ellgroup(N, subiu(q,1), pt_m, (void*)&e, &FlxqE_group, _FlxqE_pairorder);
}

GEN
Flxq_ellgens(GEN a4, GEN a6, GEN ch, GEN D, GEN m, GEN T, ulong p)
{
  GEN P;
  pari_sp av = avma;
  struct _FlxqE e;
  e.a4=a4; e.a6=a6; e.T=T; e.p=p; e.pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  switch(lg(D)-1)
  {
  case 0:
    return cgetg(1,t_VEC);
  case 1:
    P = gen_gener(gel(D,1), (void*)&e, &FlxqE_group);
    P = mkvec(FlxqE_changepoint(P, ch, T, p));
    break;
  default:
    P = gen_ellgens(gel(D,1), gel(D,2), m, (void*)&e, &FlxqE_group, _FlxqE_pairorder);
    gel(P,1) = FlxqE_changepoint(gel(P,1), ch, T, p);
    gel(P,2) = FlxqE_changepoint(gel(P,2), ch, T, p);
    break;
  }
  return gerepilecopy(av, P);
}
/***********************************************************************/
/**                          Point counting                           **/
/***********************************************************************/

/* assume a and n are coprime */
static GEN
RgX_circular_shallow(GEN P, long a, long n)
{
  long i, l = lgpol(P);
  GEN Q = cgetg(2+n,t_POL);
  Q[1] = P[1];
  for(i=0; i<l; i++)
    gel(Q,2+(i*a)%n) = gel(P,2+i);
  for(   ; i<n; i++)
    gel(Q,2+(i*a)%n) = gen_0;
  return normalizepol_lg(Q,2+n);
}

static GEN
ZpXQ_frob_cyc(GEN x, GEN T, GEN q, ulong p)
{
  long n = get_FpX_degree(T);
  return FpX_rem(RgX_circular_shallow(x,p,n+1), T, q);
}

static GEN
ZpXQ_frob(GEN x, GEN Xm, GEN T, GEN q, ulong p)
{
  if (lg(Xm)==1)
    return ZpXQ_frob_cyc(x, T, q, p);
  else
  {
    long n = get_FpX_degree(T);
    GEN V = RgX_blocks(RgX_inflate(x, p), n, p);
    GEN W = ZXV_dotproduct(V, Xm);
    return FpX_rem(W, T, q);
  }
}

struct _lift_lin
{
  ulong p, pi;
  GEN sqx, Tp, ai, Xm;
};

static GEN
_lift_invl(void *E, GEN x)
{
  struct _lift_lin *d = (struct _lift_lin *) E;
  GEN T = d->Tp;
  ulong p = d->p, pi = d->pi;
  GEN xai = Flxq_mul_pre(ZX_to_Flx(x, p), d->ai, T, p, pi);
  return Flx_to_ZX(Flxq_lroot_fast_pre(xai, d->sqx, T, p, pi));
}
static GEN
_lift_lin(void *E, GEN F, GEN x2, GEN q)
{
  struct _lift_lin *d = (struct _lift_lin *) E;
  pari_sp av = avma;
  GEN T = gel(F,3), Xm = gel(F,4);
  GEN y2  = ZpXQ_frob(x2, Xm, T, q, d->p);
  GEN lin = FpX_add(ZX_mul(gel(F,1), y2), ZX_mul(gel(F,2), x2), q);
  return gerepileupto(av, FpX_rem(lin, T, q));
}

static GEN
FpM_FpXV_bilinear(GEN P, GEN X, GEN Y, GEN p)
{
   pari_sp av = avma;
   GEN s =  ZX_mul(FpXV_FpC_mul(X,gel(P,1),p),gel(Y,1));
   long i, l = lg(P);
   for(i=2; i<l; i++)
     s = ZX_add(s, ZX_mul(FpXV_FpC_mul(X,gel(P,i),p),gel(Y,i)));
   return gerepileupto(av, FpX_red(s, p));
}

static GEN
FpM_FpXQV_bilinear(GEN P, GEN X, GEN Y, GEN T, GEN p)
{ return FpX_rem(FpM_FpXV_bilinear(P,X,Y,p),T,p); }

static GEN
FpXC_powderiv(GEN M, GEN p)
{
  long i, l;
  long v = varn(gel(M,2));
  GEN m = cgetg_copy(M, &l);
  gel(m,1) = pol_0(v);
  gel(m,2) = pol_1(v);
  for(i=2; i<l-1; i++)
    gel(m,i+1) = FpX_Fp_mul(gel(M,i),utoi(i), p);
  return m;
}

struct _lift_iso
{
  GEN phi, Xm, T, sqx, Tp;
  ulong p, pi;
};

static GEN
_lift_iter(void *E, GEN x2, GEN q)
{
  struct _lift_iso *d = (struct _lift_iso *) E;
  ulong p = d->p;
  long n = lg(d->phi)-2;
  GEN TN = FpXT_red(d->T, q), XN = FpXV_red(d->Xm, q);
  GEN y2 = ZpXQ_frob(x2, XN, TN, q, p);
  GEN xp = FpXQ_powers(x2, n, TN, q);
  GEN yp = FpXQ_powers(y2, n, TN, q);
  GEN V  = FpM_FpXQV_bilinear(d->phi,xp,yp,TN,q);
  return mkvec3(V,xp,yp);
}

static GEN
_lift_invd(void *E, GEN V, GEN v, GEN qM, long M)
{
  struct _lift_iso *d = (struct _lift_iso *) E;
  struct _lift_lin e;
  ulong p = d->p, pi = d->pi;
  GEN TM = FpXT_red(d->T, qM), XM = FpXV_red(d->Xm, qM);
  GEN xp = FpXV_red(gel(v,2), qM);
  GEN yp = FpXV_red(gel(v,3), qM);
  GEN Dx = FpM_FpXQV_bilinear(d->phi, FpXC_powderiv(xp, qM), yp, TM, qM);
  GEN Dy = FpM_FpXQV_bilinear(d->phi, xp, FpXC_powderiv(yp, qM), TM, qM);
  GEN F = mkvec4(Dy, Dx, TM, XM);
  e.ai = Flxq_inv_pre(ZX_to_Flx(Dy,p),d->Tp, p, pi);
  e.sqx = d->sqx; e.Tp = d->Tp; e.p=p; e.pi=pi; e.Xm = XM;
  return gen_ZpX_Dixon(F,V,qM,utoipos(p),M,(void*) &e, _lift_lin, _lift_invl);
}

static GEN
lift_isogeny(GEN phi, GEN x0, long n, GEN Xm, GEN T, GEN sqx, GEN Tp,
  ulong p, ulong pi)
{
  struct _lift_iso d;
  d.phi = phi; d.Xm = Xm; d.T = T;
  d.sqx = sqx; d.Tp = Tp; d.p = p; d.pi = pi;
  return gen_ZpX_Newton(x0, utoipos(p), n,(void*)&d, _lift_iter, _lift_invd);
}

static GEN
getc2(GEN act, GEN X, GEN T, GEN q, ulong p, long N)
{
  GEN A1 = RgV_to_RgX(gel(act,1),0), A2 =  RgV_to_RgX(gel(act,2),0);
  long n = brent_kung_optpow(maxss(degpol(A1),degpol(A2)),2,1);
  GEN xp = FpXQ_powers(X,n,T,q);
  GEN P  = FpX_FpXQV_eval(A1, xp, T, q);
  GEN Q  = FpX_FpXQV_eval(A2, xp, T, q);
  return ZpXQ_div(P, Q, T, q, utoipos(p), N);
}

struct _ZpXQ_norm
{
  long n;
  GEN T, p;
};

static GEN
ZpXQ_norm_mul(void *E, GEN x, GEN y)
{
  struct _ZpXQ_norm *D = (struct _ZpXQ_norm*)E;
  GEN P = gel(x,1), Q = gel(y,1);
  long a = mael(x,2,1), b = mael(y,2,1);
  retmkvec2(FpXQ_mul(P,ZpXQ_frob_cyc(Q, D->T, D->p, a), D->T, D->p),
            mkvecsmall((a*b)%D->n));
}
static GEN
ZpXQ_norm_sqr(void *E, GEN x) { return ZpXQ_norm_mul(E, x, x); }

/* Assume T = Phi_(n) and n prime */
GEN
ZpXQ_norm_pcyc(GEN x, GEN T, GEN q, GEN p)
{
  GEN z;
  struct _ZpXQ_norm D;
  long d = get_FpX_degree(T);
  D.T = T; D.p = q; D.n = d+1;
  if (d==1) return ZX_copy(x);
  z = mkvec2(x,mkvecsmall(p[2]));
  z = gen_powu_i(z,d,(void*)&D,ZpXQ_norm_sqr,ZpXQ_norm_mul);
  return gmael(z,1,2);
}

/* Assume T = Phi_(n) and n prime */
static GEN
ZpXQ_sqrtnorm_pcyc(GEN x, GEN T, GEN q, GEN p, long e)
{
  GEN z = ZpXQ_norm_pcyc(x, T, q, p);
  return Zp_sqrtlift(z,Fp_sqrt(z,p),p,e);
}

/* Assume a = 1 [p], return the square root of the norm */
static GEN
ZpXQ_sqrtnorm(GEN a, GEN T, GEN q, GEN p, long e)
{
  GEN s = Fp_div(FpXQ_trace(ZpXQ_log(a, T, p, e), T, q), gen_2, q);
  return modii(gel(Qp_exp(cvtop(s, p, e-1)),4), q);
}

struct _teich_lin
{
  ulong p, pi;
  GEN sqx, Tp;
  long m;
};

static GEN
_teich_invl(void *E, GEN x)
{
  struct _teich_lin *d = (struct _teich_lin *) E;
  ulong p = d->p, pi = d->pi;
  return Flx_to_ZX(Flxq_lroot_fast_pre(ZX_to_Flx(x,p), d->sqx, d->Tp, p, pi));
}

static GEN
_teich_lin(void *E, GEN F, GEN x2, GEN q)
{
  struct _teich_lin *d = (struct _teich_lin *) E;
  pari_sp av = avma;
  GEN T = gel(F,2), Xm = gel(F,3);
  GEN y2  = ZpXQ_frob(x2, Xm, T, q, d->p);
  GEN lin = FpX_sub(y2, ZX_mulu(ZX_mul(gel(F,1), x2), d->p), q);
  return gerepileupto(av, FpX_rem(lin, T, q));
}

struct _teich_iso
{
  GEN Xm, T, sqx, Tp;
  ulong p, pi;
};

static GEN
_teich_iter(void *E, GEN x2, GEN q)
{
  struct _teich_iso *d = (struct _teich_iso *) E;
  ulong p = d->p;
  GEN TN = FpXT_red(d->T, q), XN = FpXV_red(d->Xm, q);
  GEN y2 = ZpXQ_frob(x2, XN, TN, q, d->p);
  GEN x1 = FpXQ_powu(x2, p-1, TN, q);
  GEN xp = FpXQ_mul(x2, x1, TN, q);
  GEN V = FpX_sub(y2,xp,q);
  return mkvec2(V,x1);
}

static GEN
_teich_invd(void *E, GEN V, GEN v, GEN qM, long M)
{
  struct _teich_iso *d = (struct _teich_iso *) E;
  struct _teich_lin e;
  ulong p = d->p;
  GEN TM = FpXT_red(d->T, qM), XM = FpXV_red(d->Xm, qM);
  GEN x1 = FpX_red(gel(v,2), qM);
  GEN F = mkvec3(x1, TM, XM);
  e.sqx = d->sqx; e.Tp = d->Tp; e.p = p; e.pi = d->pi;
  return gen_ZpX_Dixon(F,V,qM,utoipos(p),M,(void*) &e, _teich_lin, _teich_invl);
}

static GEN
Teichmuller_lift(GEN x, GEN Xm, GEN T, GEN sqx, GEN Tp, ulong p, ulong pi,
  long N)
{
  struct _teich_iso d;
  d.Xm = Xm; d.T = T; d.sqx = sqx; d.Tp = Tp; d.p = p; d.pi = pi;
  return gen_ZpX_Newton(x,utoipos(p), N,(void*)&d, _teich_iter, _teich_invd);
}

static GEN
get_norm(GEN a4, GEN a6, GEN T, ulong p, ulong pi, long N)
{
  long sv=T[1];
  GEN a;
  if (p==3) a = gel(a4,1);
  else
  {
    GEN P = mkpoln(4, pol1_Flx(sv), pol0_Flx(sv), a4, a6);
    a = gel(FlxqX_powu_pre(P, p>>1, T,p,pi), 2+p-1);
  }
  return Zp_sqrtnlift(gen_1,subss(p,1),utoi(Flxq_norm(a,T,p)),utoipos(p), N);
}

static GEN
fill_pols(long n, const long *v, long m, const long *vn,
          const long *vd, GEN *act)
{
  long i, j;
  long d = upowuu(n,12/(n-1));
  GEN N, D, M = zeromatcopy(n+1,n+1);
  gmael(M,1,n+1) = gen_1;
  for (i = 2; i <= n+1; i++)
    for (j = i-1; j <= n; j++) gmael(M,i,j) = mulis(powuu(d,i-2), v[j-i+1]);
  N = cgetg(m+1,t_COL);
  D = cgetg(m+1,t_COL);
  for(i = 1; i <= m; i++)
  {
    gel(N,i) = stoi(*vn++);
    gel(D,i) = stoi(*vd++);
  }
  *act = mkmat2(N,D); return M;
}

/*
  These polynomials were extracted from the ECHIDNA databases
  available at <http://echidna.maths.usyd.edu.au/echidna/>
  and computed by David R. Kohel.
  Return the matrix of the modular polynomial, set act to the parametrization,
  and set dj to the opposite of the supersingular j-invariant.
*/
static GEN
get_Kohel_polynomials(ulong p, GEN *act, long *dj)
{
  const long mat3[] = {-1,-36,-270};
  const long num3[] = {1,-483,-21141,-59049};
  const long den3[] = {1,261, 4347, -6561};
  const long mat5[] = {-1,-30,-315,-1300,-1575};
  const long num5[] = {-1,490,20620,158750,78125};
  const long den5[] = {-1,-254,-4124,-12250,3125};
  const long mat7[] = {-1,-28,-322,-1904,-5915,-8624,-4018};
  const long num7[] = {1,-485,-24058,-343833,-2021642,-4353013,-823543};
  const long den7[] = {1,259,5894,49119,168406,166355,-16807};
  const long mat13[]= {-1,-26,-325,-2548,-13832,-54340,-157118,-333580,-509366,
                       -534820,-354536,-124852,-15145};
  const long num13[]= {1,-487,-24056,-391463,-3396483,-18047328,-61622301,
                       -133245853,-168395656,-95422301,-4826809};
  const long den13[]= {1,257,5896,60649,364629,1388256,3396483,5089019,4065464,
                       1069939,-28561};
  switch(p)
  {
  case 3:
    *dj = 0;
    return fill_pols(3,mat3,4,num3,den3,act);
  case 5:
    *dj = 0;
    return fill_pols(5,mat5,5,num5,den5,act);
  case 7:
    *dj = 1;
    return fill_pols(7,mat7,7,num7,den7,act);
  case 13:
    *dj = 8;
    return fill_pols(13,mat13,11,num13,den13,act);
  }
  *dj=0; *act = NULL; return NULL; /* LCOV_EXCL_LINE */
}

long
zx_is_pcyc(GEN T)
{
  long i, n = degpol(T);
  if (!uisprime(n+1)) return 0;
  for (i = 0; i <= n; i++)
    if (T[i+2] != 1UL) return 0;
  return 1;
}

static GEN
Flxq_ellcard_Kohel(GEN a4, GEN a6, GEN T, ulong p)
{
  pari_sp av = avma, av2;
  pari_timer ti;
  long n = get_Flx_degree(T), N = (n+4)/2, dj;
  GEN q = powuu(p, N);
  GEN T2, Xm, s1, c2, t, lr, S1, sqx, Nc2, Np;
  GEN act, phi = get_Kohel_polynomials(p, &act, &dj);
  long ispcyc = zx_is_pcyc(get_Flx_mod(T));
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  timer_start(&ti);
  if (!ispcyc)
  {
    T2 = Flx_Teichmuller(get_Flx_mod(T),p,N);
    if (DEBUGLEVEL) timer_printf(&ti,"Teich");
  } else
    T2 = Flx_to_ZX(get_Flx_mod(T));

  T2 = FpX_get_red(T2, q); T = ZXT_to_FlxT(T2, p);
  av2 = avma;
  if (DEBUGLEVEL) timer_printf(&ti,"Barrett");
  if (!ispcyc)
  {
    Xm = FpXQ_powers(pol_xn(n,get_FpX_var(T2)),p-1,T2,q);
    if (DEBUGLEVEL) timer_printf(&ti,"Xm");
  } else
    Xm = cgetg(1,t_VEC);
  s1 = Flxq_inv_pre(Flx_Fl_add(Flxq_ellj(a4,a6,T,p),dj, p),T,p,pi);
  lr = Flxq_lroot_pre(polx_Flx(get_Flx_var(T)), T,p,pi);
  sqx = Flxq_powers_pre(lr, p-1, T, p, pi);
  S1 = lift_isogeny(phi, Flx_to_ZX(s1), N, Xm, T2, sqx, T,p,pi);
  if (DEBUGLEVEL) timer_printf(&ti,"Lift isogeny");
  c2 = getc2(act, S1, T2, q, p, N);
  if (DEBUGLEVEL) timer_printf(&ti,"c^2");
  if (p>3 && !ispcyc)
  {
    GEN c2p = Flx_to_ZX(Flxq_inv_pre(ZX_to_Flx(c2,p),T,p,pi));
    GEN tc2 = Teichmuller_lift(c2p,Xm, T2,sqx,T,p,pi,N);
    if (DEBUGLEVEL) timer_printf(&ti,"Teichmuller/Fq");
    c2 = FpX_rem(FpX_mul(tc2,c2,q),T2,q);
  }
  c2 = gerepileupto(av2, c2);
  if (DEBUGLEVEL) timer_printf(&ti,"tc2");
  Nc2 = (ispcyc? ZpXQ_sqrtnorm_pcyc: ZpXQ_sqrtnorm)(c2, T2, q, utoipos(p), N);
  if (DEBUGLEVEL) timer_printf(&ti,"Norm");
  Np = get_norm(a4,a6,T,p,pi,N);
  if (p>3 && ispcyc)
  {
    GEN Ncpi =  utoi(Fl_inv(umodiu(Nc2,p), p));
    GEN tNc2 = Zp_sqrtnlift(gen_1, subss(p,1), Ncpi, utoipos(p),N);
    if (DEBUGLEVEL) timer_printf(&ti,"Teichmuller/Fp");
    Nc2 = Fp_mul(Nc2,tNc2,q);
  }
  t = Fp_center_i(Fp_mul(Nc2,Np,q),q,shifti(q,-1));
  return gerepileupto(av, subii(addiu(powuu(p,n),1),t));
}

/* Use Damien Robert's method */
static GEN
get_trace_Robert(GEN J, GEN phi, GEN Xm, GEN T, GEN q, ulong p, long e)
{
  long n = lg(phi)-2;
  GEN K = ZpXQ_frob(J, Xm, T, q, p);
  GEN Jp = FpXQ_powers(J, n, T, q);
  GEN Kp = FpXQ_powers(K, n, T, q);
  GEN Jd = FpXC_powderiv(Jp, q);
  GEN Kd = FpXC_powderiv(Kp, q);
  GEN Dx = FpM_FpXQV_bilinear(phi, Kd, Jp, T, q);
  GEN Dy = FpM_FpXQV_bilinear(phi, Kp, Jd, T, q);
  GEN C = ZpXQ_inv(ZX_divuexact(Dy, p), T, utoi(p), e);
  return FpX_neg(FpXQ_mul(Dx, C, T, q), q);
}

/* in p^2, so p is tiny */
static GEN
Flxq_ellcard_Harley(GEN a4, GEN a6, GEN T, ulong p)
{
  pari_sp av = avma, av2;
  pari_timer ti;
  long n = get_Flx_degree(T), N = (n+5)/2;
  GEN pp = utoipos(p), q = powuu(p, N);
  GEN T2, j, t, phi, J1, sqx, Xm, c2, tc2, c2p, Nc2, Np;
  long ispcyc = zx_is_pcyc(get_Flx_mod(T));
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p); /* = 0 here */
  timer_start(&ti);
  if (!ispcyc)
  {
    T2 = Flx_Teichmuller(get_Flx_mod(T),p,N);
    if (DEBUGLEVEL) timer_printf(&ti,"Teich");
  } else
    T2 = Flx_to_ZX(get_Flx_mod(T));
  T2 = FpX_get_red(T2, q); T = ZXT_to_FlxT(T2, p);
  av2 = avma;
  if (DEBUGLEVEL) timer_printf(&ti,"Barrett");
  if (!ispcyc)
  {
    Xm = FpXQ_powers(pol_xn(n,get_FpX_var(T2)),p-1,T2,q);
    if (DEBUGLEVEL) timer_printf(&ti,"Xm");
  } else
    Xm = cgetg(1,t_VEC);
  j = Flxq_ellj(a4,a6,T,p);
  sqx = Flxq_powers_pre(Flxq_lroot_pre(polx_Flx(T[1]), T,p,pi), p-1, T,p,pi);
  phi = polmodular_ZM(p, 0);
  if (DEBUGLEVEL) timer_printf(&ti,"phi");
  J1 = lift_isogeny(phi, Flx_to_ZX(j), N, Xm, T2,sqx,T,p,pi);
  if (DEBUGLEVEL) timer_printf(&ti,"Lift isogeny");
  c2 = get_trace_Robert(J1, phi, Xm, T2, q, p, N);
  q = diviuexact(q,p); N--;
  if (DEBUGLEVEL) timer_printf(&ti,"c^2");
  if (!ispcyc)
  {
    c2p = Flx_to_ZX(Flxq_inv_pre(ZX_to_Flx(c2,p),T,p,pi));
    tc2 = Teichmuller_lift(c2p,Xm, T2,sqx,T,p,pi,N);
    if (DEBUGLEVEL) timer_printf(&ti,"teichmuller");
    c2 = FpX_rem(FpX_mul(tc2,c2,q),T2,q);
  }
  c2 = gerepileupto(av2, c2);
  q = powuu(p, N);
  Nc2 = (ispcyc? ZpXQ_sqrtnorm_pcyc: ZpXQ_sqrtnorm)(c2, T2, q, pp, N);
  if (DEBUGLEVEL) timer_printf(&ti,"Norm");
  Np = get_norm(a4,a6,T,p,pi,N);
  if (ispcyc)
  {
    GEN Ncpi = utoi(Fl_inv(umodiu(Nc2,p), p));
    GEN tNc2 = Zp_sqrtnlift(gen_1, subss(p,1), Ncpi, pp, N);
    if (DEBUGLEVEL) timer_printf(&ti,"Teichmuller/Fp");
    Nc2 = Fp_mul(Nc2,tNc2,q);
  }
  t = Fp_center_i(Fp_mul(Nc2,Np,q),q,shifti(q,-1));
  return gerepileupto(av, subii(addiu(powuu(p,n),1),t));
}

/***************************************************************************/
/*                          Shanks-Mestre                                  */
/***************************************************************************/

/* Return the lift of a (mod b), which is closest to h */
static GEN
closest_lift(GEN a, GEN b, GEN h)
{ return addii(a, mulii(b, diviiround(subii(h,a), b))); }

/* find multiple of order of f using Baby Step/Giant Step, f^h close to 1,
 * order lies in an interval of size <= 'bound' and known mod B */
static GEN
_FlxqE_order_multiple(void *E, GEN f, GEN h, GEN bound, GEN B)
{
  pari_sp av = avma, av1;
  pari_timer Ti;
  long i, s = ceilsqrtdiv(bound, B) >> 1;
  GEN P, F, tx, ti, fg, fh;

  P = fh = _FlxqE_mul(E, f, h);
  if (DEBUGLEVEL >= 6) timer_start(&Ti);
  if (ell_is_inf(fh)) return h;
  F = _FlxqE_mul(E, f, B);
  if (s < 3)
  { /* we're nearly done: naive search */
    GEN Q = P;
    for (i=1;; i++)
    {
      P = _FlxqE_add(E, P, F); /* h.f + i.F */
      if (ell_is_inf(P)) return gerepileupto(av, addii(h, mului(i,B)));
      Q = _FlxqE_sub(E, Q, F); /* h.f - i.F */
      if (ell_is_inf(Q)) return gerepileupto(av, subii(h, mului(i,B)));
    }
  }
  tx = cgetg(s+1,t_VECSMALL); av1 = avma;
  for (i=1; i<=s; i++)
  { /* baby steps */
    tx[i] = hash_GEN(gel(P, 1));
    P = _FlxqE_add(E, P, F); /* h.f + i.F */
    if (ell_is_inf(P)) return gerepileupto(av, addii(h, mului(i,B)));
    if (gc_needed(av1,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"[Flxq_ellcard] baby steps, i=%ld",i);
      P = gerepileupto(av1,P);
    }
  }
  if (DEBUGLEVEL >= 6) timer_printf(&Ti,"[Flxq_ellcard] baby steps, s = %ld",s);
  /* giant steps: fg = s.F */
  fg = gerepileupto(av1, _FlxqE_sub(E, P, fh));
  if (ell_is_inf(fg)) return gerepileupto(av, mului(s,B));
  ti = vecsmall_indexsort(tx); /* = permutation sorting tx */
  tx = perm_mul(tx,ti);
  if (DEBUGLEVEL >= 6) timer_printf(&Ti, "[Flxq_ellcard] sorting");
  av1 = avma;
  for (P=fg, i=1; ; i++)
  {
    long k = hash_GEN(gel(P,1)), r = zv_search(tx, k);
    if (r)
    {
      while (r && tx[r] == k) r--;
      for (r++; r <= s && tx[r] == k; r++)
      {
        long j = ti[r]-1;
        GEN Q = _FlxqE_add(E, _FlxqE_mul(E, F, stoi(j)), fh);
        if (DEBUGLEVEL >= 6)
          timer_printf(&Ti, "[Flxq_ellcard] giant steps, i = %ld",i);
        if (Flx_equal(gel(P,1), gel(Q,1)))
        {
          if (Flx_equal(gel(P,2), gel(Q,2))) i = -i;
          return gerepileupto(av,addii(h, mulii(addis(mulss(s,i), j), B)));
        }
      }
    }
    P = _FlxqE_add(E, P, fg);
    if (gc_needed(av1,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"[Flxq_ellcard] giants steps, i=%ld",i);
      P = gerepileupto(av1,P);
    }
  }
}
static GEN
_FlxqE_order(void *E, GEN f, GEN h, GEN bound, GEN B)
{
  GEN o = _FlxqE_order_multiple(E, f, h, bound, B);
  return gen_order(f, o, E, &FlxqE_group);
}

static void
Flx_next(GEN t, ulong p)
{
  long i;
  for(i=2;;i++)
    if (uel(t,i)==p-1) t[i]=0; else { t[i]++; break; }
}

static void
Flx_renormalize_ip(GEN x, long lx)
{
  long i;
  for (i = lx-1; i>=2; i--)
    if (x[i]) break;
  setlg(x, i+1);
}

static ulong
F3xq_ellcard_naive(GEN a2, GEN a6, GEN T)
{
  pari_sp av = avma;
  long i, d = get_Flx_degree(T), lx = d+2;
  long q = upowuu(3, d), a;
  GEN x = zero_zv(lx); x[1] = get_Flx_var(T);
  for(a=1, i=0; i<q; i++)
  {
    GEN rhs;
    Flx_renormalize_ip(x, lx);
    rhs = Flx_add(Flxq_mul(Flxq_sqr(x, T, 3), Flx_add(x, a2, 3), T, 3), a6, 3);
    if (!lgpol(rhs)) a++; else if (Flxq_issquare(rhs, T, 3)) a+=2;
    Flx_next(x, 3);
  }
  set_avma(av); return a;
}

/* p^deg(T) is tiny */
static ulong
Flxq_ellcard_naive(GEN a4, GEN a6, GEN T, ulong p)
{
  pari_sp av = avma;
  long i, d = get_Flx_degree(T), lx = d+2;
  long q = upowuu(p, d), a;
  GEN x = zero_zv(lx); x[1] = get_Flx_var(T);
  for(a = 1, i = 0; i < q; i++)
  {
    GEN x2, rhs;
    Flx_renormalize_ip(x, lx);
    x2  = Flxq_sqr_pre(x, T, p, 0);
    rhs = Flx_add(Flxq_mul_pre(x, Flx_add(x2, a4, p), T, p, 0), a6, p);
    if (!lgpol(rhs)) a++; else if (Flxq_issquare(rhs,T,p)) a += 2;
    Flx_next(x,p);
  }
  set_avma(av); return a;
}

static long
Flxq_kronecker(GEN x, GEN T, ulong p)
{
  pari_sp av;
  if (lgpol(x) == 0) return 0;
  av = avma; return gc_long(av, krouu(Flxq_norm(x, T, p), p));
}

/* Find x such that kronecker(u = x^3+a4x+a6, p) is KRO.
 * Return point [x*u,u^2] on E (KRO=1) / E^twist (KRO=-1) */
static GEN
Flxq_ellpoint(long KRO, GEN a4, GEN a6, GEN T, ulong p, ulong pi)
{
  long v = get_Flx_var(T), n = get_Flx_degree(T);
  for(;;)
  {
    GEN x = random_Flx(n, v, p), x2 = Flxq_sqr_pre(x,T,p,pi);
    GEN u = Flx_add(a6, Flxq_mul_pre(Flx_add(a4, x2, p), x, T,p, pi), p);
    if (Flxq_kronecker(u,T,p) == KRO)
      return mkvec2(Flxq_mul_pre(u,x, T,p,pi), Flxq_sqr_pre(u, T,p,pi));
  }
}

static GEN
Flxq_ellcard_Shanks(GEN a4, GEN a6, GEN q, GEN T, ulong p)
{
  pari_sp av = avma;
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  long v = get_Flx_var(T), KRO = -1;
  GEN h,f, A, B;
  GEN q1p = addiu(q,1), q2p = shifti(q1p, 1);
  GEN bound = addiu(sqrti(gmul2n(q,4)), 1); /* ceil( 4sqrt(q) ) */
  struct _FlxqE e;
  e.p = p; e.pi = pi; e.T = Flx_get_red_pre(T, p, pi);
  /* once #E(Flxq) is known mod B >= bound, it is determined */
  switch(FlxqX_nbroots(mkpoln(4, pol1_Flx(v), pol0_Flx(v), a4, a6), T, p))
  { /* how many 2-torsion points ? */
  case 3:  A = gen_0; B = utoipos(4); break;
  case 1:  A = gen_0; B = gen_2; break;
  default: A = gen_1; B = gen_2; break; /* 0 */
  }
  for(;;)
  {
    h = closest_lift(A, B, q1p);
    /* [ux, u^2] is on E_u: y^2 = x^3 + c4 u^2 x + c6 u^3
     * E_u isomorphic to E (resp. E') iff KRO = 1 (resp. -1)
     * #E(F_p) = p+1 - a_p, #E'(F_p) = p+1 + a_p
     *
     * #E_u(Flxq) = A (mod B),  h is close to #E_u(Flxq) */
    KRO = -KRO;
    f = Flxq_ellpoint(KRO, a4,a6, T,p,pi);
    e.a4 = Flxq_mul_pre(a4, gel(f,2), T,p,pi); /* a4 for E_u */
    h = _FlxqE_order((void*)&e, f, h, bound, B);
    /* h | #E_u(Flxq) = A (mod B) */
    A = Z_chinese_all(A, gen_0, B, h, &B);
    if (cmpii(B, bound) >= 0) break;
    /* not done, update A mod B for the _next_ curve, isomorphic to
     * the quadratic twist of this one */
    A = remii(subii(q2p,A), B); /* #E(Fq)+#E'(Fq) = 2q+2 */
  }
  h = closest_lift(A, B, q1p);
  return gerepileuptoint(av, KRO == 1? h: subii(q2p,h));
}

static GEN
F3xq_ellcard(GEN a2, GEN a6, GEN T)
{
  long n = get_Flx_degree(T);
  if (n <= 2)
    return utoi(F3xq_ellcard_naive(a2, a6, T));
  else
  {
    GEN q1 = addiu(powuu(3, get_Flx_degree(T)), 1), t;
    GEN a = Flxq_div(a6,Flxq_powu(a2,3,T,3),T,3);
    if (Flx_equal1(Flxq_powu(a, 8, T, 3)))
    {
      GEN P = Flxq_minpoly(a,T,3);
      long dP = degpol(P); /* dP <= 2 */
      ulong q = upowuu(3,dP);
      GEN A2 = pol1_Flx(P[1]), A6 = Flx_rem(polx_Flx(P[1]), P, 3);
      long tP = q + 1 - F3xq_ellcard_naive(A2, A6, P);
      t = elltrace_extension(stoi(tP), n/dP, utoi(q));
      if (umodiu(t, 3)!=1) t = negi(t);
      return Flx_equal1(a2) || Flxq_issquare(a2,T,3) ? subii(q1,t): addii(q1,t);
    }
    else return Flxq_ellcard_Kohel(mkvec(a2), a6, T, 3);
  }
}

static GEN
Flxq_ellcard_Satoh(GEN a4, GEN a6, GEN j, GEN T, ulong p)
{
  long n = get_Flx_degree(T);
  if (n <= 2)
    return utoi(Flxq_ellcard_naive(a4, a6, T, p));
  else
  {
    GEN jp = Flxq_powu(j, p, T, p);
    GEN s = Flx_add(j, jp, p);
    if (degpol(s) <= 0)
    { /* it is assumed j not in F_p */
      GEN m = Flxq_mul(j, jp, T, p);
      if (degpol(m) <= 0)
      {
        GEN q = sqru(p);
        GEN q1 = addiu(powuu(p, get_Flx_degree(T)), 1);
        GEN sk = Flx_Fl_add(Flx_neg(j, p), 1728%p, p);
        GEN sA4 = Flx_triple(Flxq_mul(sk, j, T, p), p);
        GEN u = Flxq_div(a4, sA4, T, p);
        ulong ns = lgpol(s) ? Fl_neg(s[2], p): 0UL;
        GEN P = mkvecsmall4(T[1], m[2], ns, 1L);
        GEN A4, A6, t, tP;
        Flxq_ellj_to_a4a6(polx_Flx(T[1]), P, p, &A4, &A6);
        tP = addis(q, 1 - Flxq_ellcard_naive(A4, A6, P, p));
        t = elltrace_extension(tP, n>>1, q);
        return Flxq_is2npower(u, 2, T, p) ? subii(q1,t): addii(q1,t);
      }
    }
    if (p<=7 || p==13) return Flxq_ellcard_Kohel(a4, a6, T, p);
    else return Flxq_ellcard_Harley(a4, a6, T, p);
  }
}

static GEN
Flxq_ellcard_Kedlaya(GEN a4, GEN a6, GEN T, ulong p)
{
  pari_sp av = avma;
  GEN H = mkpoln(4, gen_1, gen_0, Flx_to_ZX(a4), Flx_to_ZX(a6));
  GEN Tp = Flx_to_ZX(get_Flx_mod(T));
  long n = degpol(Tp), e = ((p < 16 ? n+1: n)>>1)+1;
  GEN M = ZlXQX_hyperellpadicfrobenius(H, Tp, p, e);
  GEN N = ZpXQM_prodFrobenius(M, Tp, utoipos(p), e);
  GEN q = powuu(p, e);
  GEN tp = Fq_add(gcoeff(N,1,1), gcoeff(N,2,2), Tp, q);
  GEN t = Fp_center_i(typ(tp)==t_INT ? tp: leading_coeff(tp), q, shifti(q,-1));
  return gerepileupto(av, subii(addiu(powuu(p, n), 1), t));
}

GEN
Flxq_ellj(GEN a4, GEN a6, GEN T, ulong p)
{
  pari_sp av=avma;
  if (p==3)
  {
    GEN J;
    if (typ(a4)!=t_VEC) return pol0_Flx(get_Flx_var(T));
    J = Flxq_div(Flxq_powu(gel(a4,1),3, T, p),Flx_neg(a6,p), T, p);
    return gerepileuptoleaf(av, J);
  }
  else
  {
    pari_sp av=avma;
    GEN a43 = Flxq_mul(a4,Flxq_sqr(a4,T,p),T,p);
    GEN a62 = Flxq_sqr(a6,T,p);
    GEN num = Flx_mulu(a43,6912,p);
    GEN den = Flx_add(Flx_mulu(a43,4,p),Flx_mulu(a62,27,p),p);
    return gerepileuptoleaf(av, Flxq_div(num, den, T, p));
  }
}

void
Flxq_ellj_to_a4a6(GEN j, GEN T, ulong p, GEN *pt_a4, GEN *pt_a6)
{
  ulong zagier = 1728 % p;
  if (lgpol(j)==0)
    { *pt_a4 = pol0_Flx(T[1]); *pt_a6 =pol1_Flx(T[1]); }
  else if (lgpol(j)==1 && uel(j,2) == zagier)
    { *pt_a4 = pol1_Flx(T[1]); *pt_a6 =pol0_Flx(T[1]); }
  else
  {
    GEN k = Flx_Fl_add(Flx_neg(j, p), zagier, p);
    GEN kj = Flxq_mul(k, j, T, p);
    GEN k2j = Flxq_mul(kj, k, T, p);
    *pt_a4 = Flx_triple(kj, p);
    *pt_a6 = Flx_double(k2j, p);
  }
}

static GEN
F3xq_ellcardj(GEN a4, GEN a6, GEN T, GEN q, long n)
{
  const ulong p = 3;
  ulong t;
  GEN q1 = addiu(q,1);
  GEN na4 = Flx_neg(a4,p), ra4;
  if (!Flxq_issquare(na4,T,p))
    return q1;
  ra4 = Flxq_sqrt(na4,T,p);
  t = Flxq_trace(Flxq_div(a6,Flxq_mul(na4,ra4,T,p),T,p),T,p);
  if (n%2==1)
  {
    GEN q3;
    if (t==0) return q1;
    q3 = powuu(p,(n+1)>>1);
    return (t==1)^(n%4==1) ? subii(q1,q3): addii(q1,q3);
  }
  else
  {
    GEN q22, q2 = powuu(p,n>>1);
    GEN W = Flxq_pow(a4,shifti(q,-2),T,p);
    long s = (W[2]==1)^(n%4==2);
    if (t!=0) return s ? addii(q1,q2): subii(q1, q2);
    q22 = shifti(q2,1);
    return s ? subii(q1,q22):  addii(q1, q22);
  }
}

static GEN
Flxq_ellcardj(GEN a4, GEN a6, ulong j, GEN T, GEN q, ulong p, long n)
{
  GEN q1 = addiu(q,1);
  if (j==0)
  {
    ulong w;
    GEN W, t, N;
    if (umodiu(q,6)!=1) return q1;
    N = Fp_ffellcard(gen_0,gen_1,q,n,utoipos(p));
    t = subii(q1, N);
    W = Flxq_pow(a6,diviuexact(shifti(q,-1), 3),T,p);
    if (degpol(W)>0) /*p=5 mod 6*/
      return Flx_equal1(Flxq_powu(W,3,T,p)) ? addii(q1,shifti(t,-1)):
                                              subii(q1,shifti(t,-1));
    w = W[2];
    if (w==1)   return N;
    if (w==p-1) return addii(q1,t);
    else /*p=1 mod 6*/
    {
      GEN u = shifti(t,-1), v = sqrtint(diviuexact(subii(q,sqri(u)),3));
      GEN a = addii(u,v), b = shifti(v,1);
      if (Fl_powu(w,3,p)==1)
      {
        if (Fl_add(umodiu(a,p),Fl_mul(w,umodiu(b,p),p),p)==0)
          return subii(q1,subii(shifti(b,1),a));
        else
          return addii(q1,addii(a,b));
      }
      else
      {
        if (Fl_sub(umodiu(a,p),Fl_mul(w,umodiu(b,p),p),p)==0)
          return subii(q1,subii(a,shifti(b,1)));
        else
          return subii(q1,addii(a,b));
      }
    }
  } else if (j==1728%p)
  {
    ulong w;
    GEN W, N, t;
    if (mod4(q)==3) return q1;
    W = Flxq_pow(a4,shifti(q,-2),T,p);
    if (degpol(W)>0) return q1; /*p=3 mod 4*/
    w = W[2];
    N = Fp_ffellcard(gen_1,gen_0,q,n,utoipos(p));
    if(w==1) return N;
    t = subii(q1, N);
    if(w==p-1) return addii(q1, t);
    else /*p=1 mod 4*/
    {
      GEN u = shifti(t,-1), v = sqrtint(subii(q,sqri(u)));
      if (Fl_add(umodiu(u,p),Fl_mul(w,umodiu(v,p),p),p)==0)
        return subii(q1,shifti(v,1));
      else
        return addii(q1,shifti(v,1));
    }
  } else
  {
    ulong g = Fl_div(j, Fl_sub(1728%p, j, p), p);
    GEN N = Fp_ffellcard(utoi(Fl_triple(g,p)),utoi(Fl_double(g,p)),q,n,utoipos(p));
    GEN l = Flxq_mul(Flx_triple(a6,p),Flx_double(a4,p),T,p);
    if (Flxq_issquare(l,T,p)) return N;
    return subii(shifti(q1,1),N);
  }
}

static GEN
Flxq_ffellcard(GEN a4, GEN a6, GEN M, GEN q, GEN T, ulong p, long n)
{
  long m = degpol(M);
  GEN j = polx_Flx(M[1]);
  GEN g = Flxq_div(j, mkvecsmall3(M[1],1728%p,p-1), M, p);
  GEN N = Flxq_ellcard(Flx_triple(g, p), Flx_double(g, p), M, p);
  GEN qm =  powuu(p, m), q1 = addiu(q, 1), qm1 = addiu(qm, 1);
  GEN l = Flxq_mul(Flx_triple(a6,p), Flx_double(a4,p), T, p);
  GEN te = elltrace_extension(subii(qm1, N), n/m, qm);
  return Flxq_issquare(l,T,p) ? subii(q1, te): addii(q1, te);
}

static GEN
Flxq_ellcard_i(GEN a4, GEN a6, GEN T, ulong p)
{
  long n = get_Flx_degree(T);
  GEN J, M, q = powuu(p,  n);
  if (typ(a4)==t_VEC)
    return F3xq_ellcard(gel(a4,1), a6, T);
  if (p==3)
    return F3xq_ellcardj(a4, a6, T, q, n);
  if (degpol(a4)<=0 && degpol(a6)<=0)
    return Fp_ffellcard(utoi(Flx_eval(a4,0,p)),utoi(Flx_eval(a6,0,p)),q,n,utoipos(p));
  J = Flxq_ellj(a4,a6,T,p);
  if (degpol(J)<=0)
    return Flxq_ellcardj(a4,a6,lgpol(J)?J[2]:0,T,q,p,n);
  M = Flxq_minpoly(J, T, p);
  if (degpol(M) < n)
    return Flxq_ffellcard(a4, a6, M, q, T, p, n);
  if (p <= 7)
    return Flxq_ellcard_Satoh(a4, a6, J, T, p);
  if (cmpis(q,100)<0)
    return utoi(Flxq_ellcard_naive(a4, a6, T, p));
  if (p == 13 || (7*p <= (ulong)10*n && (BITS_IN_LONG==64 || p <= 103)))
    return Flxq_ellcard_Satoh(a4, a6, J, T, p);
  if (p <= (ulong)2*n)
    return Flxq_ellcard_Kedlaya(a4, a6, T, p);
  if (expi(q)<=62)
    return Flxq_ellcard_Shanks(a4, a6, q, T, p);
  else
    return Fq_ellcard_SEA(Flx_to_ZX(a4),Flx_to_ZX(a6),q,Flx_to_ZX(T),utoipos(p),0);
}

GEN
Flxq_ellcard(GEN a4, GEN a6, GEN T, ulong p)
{
  pari_sp av = avma;
  return gerepileuptoint(av, Flxq_ellcard_i(a4, a6, T, p));
}

static long
Fl_ellj_trace(ulong j, ulong p)
{
  ulong a4, a6;
  Fl_ellj_to_a4a6(j, p, &a4, &a6);
  return Fl_elltrace(a4, a6, p);
}

/* Given ordinary E/Fq, a prime ell, and the height of the ell-volcano
 * containing j(E) (= v_ell(conductor of Z[pi_E]) returns the height of j(E)
 * on its ell-volcano (= v_ell(conductor of the order End(E)). */
static long
Fl_ellheightabovefloor(ulong j, long ell, long e, ulong p)
{
  pari_sp av = avma;
  GEN Xp, G, phi, phix, j0, j1;
  long h, i, nj1;
  if (e==0) return 0;
  if (j==0 || j==1728%p) return e;
  phi = ZXX_to_FlxX(polmodular_ZXX(ell, 0, 0, 1), p, 1);
  phix = FlxY_evalx(phi, j, p);
  Xp = Flx_Frobenius(phix, p);
  G  = Flx_gcd(Flx_sub(Xp, polx_Flx(0), p), phix, p);
  nj1 = degpol(G);
  if (nj1 < ell) return 0;
  if (e==1 || nj1 != ell+1) return e;
  j1 = Flx_roots(G, p);
  nj1 = lg(j1)-1;
  if (nj1 < 3) return 0;
  j0 = mkvecsmall3(j,j,j);
  for (h = 1; ; h++)
    for(i = 1; i <= 3; i++)
    {
      GEN P = Flx_div_by_X_x(FlxY_evalx(phi, uel(j1,i), p), uel(j0,i), p, NULL);
      GEN r = Flx_roots(P, p);
      if (lg(r) == 1) return gc_long(av, h);
      j0[i] = j1[i];
      j1[i] = r[1];
    }
}

/* Given an ordinary elliptic curve E/Fp and an integer h, returns
 * D = disc(End(E)) assuming h(D) = h, using the approach sketched in
 * Remark 13. If the algorithm returns 0 it has proved that h(D) != h, but it
 * is under no obligation to do so and is allowed to return any value when the
 * assumption h(d) = h is false. */
static long
Fl_end13(ulong j, ulong h, ulong p)
{
  ulong D0, v, h0;
  long i, lL, lc, lD, nc;
  GEN D, DF, cs, L, vP, vE;
  ulong t = Fl_ellj_trace(j, p);

  D0 = coredisc2u_fact(factoru(4*p-t*t), -1, &vP, &vE);
  h0 = itou(classno(stoi(-D0)));
  if (h % h0 != 0) return 0;
  h /= h0;
  D = divisorsu_fact_factored(mkmat2(vP,vE));
  DF = gel(D,2); D = gel(D,1);
  lD = lg(D); v = D[lD-1];
  cs = cgetg(lD,t_VECSMALL); nc = 0;
  for (i = 1; i < lD; i++)
  {
    GEN F = gel(DF,i);
    ulong w = uquadclassnoF_fact(D0, -1, gel(F,1), gel(F,2));
    if (w == h) uel(cs,++nc) = v / uel(D,i);
  }
  if (nc==0) return 0;
  if (nc==1) { v /= uel(cs,1); return -D0*v*v; }
  L = cgetg(nc+1, t_VEC);
  for (i = 1; i <= nc; i++) gel(L,i) = gel(factoru(uel(cs,i)), 1);
  L = vecsmall_uniq(shallowconcat1(L));
  lL = lg(L); lc = nc+1;
  for (i = 1; i < lL; i++)
  {
    ulong ell = L[i];
    long k, e = Fl_ellheightabovefloor(j, ell, z_lval(v,ell), p);
    for (k = 1; k < lc; k++)
      if(cs[k] && z_lval(cs[k], ell) != e) { cs[k] = 0; nc--; }
    if (nc==0) return 0;
    if (nc==1)
    {
      for (k = 1; k < lc; k++)
        if (cs[k]) { v /= uel(cs,k); return -D0*v*v; }
    }
  }
  return 0;
}

INLINE int
RgX_is_monic_ZX(GEN pol)
{ return RgX_is_ZX(pol) && ZX_is_monic(pol); }

long
polisclass(GEN H)
{
  pari_sp av = avma, btop;
  long h = degpol(H), hl, i, pmin, vH = varn(H), vh;
  double lmin;
  ulong p;
  GEN h2list;
  forprime_t T;

  if (typ(H)!= t_POL) pari_err_TYPE("polsisclass",H);
  if (h <= 0 || !RgX_is_monic_ZX(H)) return 0;
  vh = vals(h);
  h2list = cgetg(vh+2, t_VECSMALL); hl = 1;
  for (i = 0; i <= vh; i++)
  {
    ulong d = 1UL<<i;
    if (((d-h)&1)==0) h2list[hl++] = d;
  }
  setlg(h2list, hl);
  lmin = h * (log(log(h+2))+2);
  pmin = 33 * ceil(lmin*lmin);
  u_forprime_init(&T, pmin, ULONG_MAX);
  btop = avma;
  while((p = u_forprime_next(&T)))
  {
    ulong r;
    long D, nroots;
    GEN Xp, G, Hp = ZX_to_Flx(H,p);
    if (!Flx_is_squarefree(Hp, p)) { set_avma(btop); continue; }
    Xp = Flx_Frobenius(Hp, p);
    G  = Flx_gcd(Flx_sub(Xp, polx_Flx(evalvarn(vH)), p), Hp, p);
    nroots = degpol(G);
    if (nroots==0) { set_avma(btop); continue; }
    if (nroots < h && !zv_search(h2list,nroots)) return gc_long(av, 0);
    r = Flx_oneroot(G, p);
    if (Fp_elljissupersingular(utoi(r), utoi(p))) { set_avma(btop); continue; }
    D = Fl_end13(r, h, p);
    if (D && gequal(H, polclass(stoi(D), 0, vH))) return gc_long(av, D);
    return gc_long(av, 0);
  }
  pari_err_BUG("polisclass");
  return 0; /* LCOV_EXCL_LINE */
}
