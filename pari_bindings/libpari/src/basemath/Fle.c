/* Copyright (C) 2014  The PARI group.

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

/* Not so fast arithmetic with points over elliptic curves over Fl */

/***********************************************************************/
/**                              Flj                                  **/
/***********************************************************************/
/* Jacobian coordinates: we represent a projective point (x:y:z) on E by
 * [z*x, z^2*y, z]. Not the fastest representation available for the problem,
 * but easy to implement and up to 60% faster than the school-book method. */

/* Cost: 1M + 8S + 1*a + 10add + 1*8 + 2*2 + 1*3.
 * http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl */
INLINE void
Flj_dbl_indir_pre(GEN P, GEN Q, ulong a4, ulong p, ulong pi)
{
  ulong X1, Y1, Z1;
  ulong XX, YY, YYYY, ZZ, S, M, T;

  X1 = P[1]; Y1 = P[2]; Z1 = P[3];

  if (Z1 == 0) { Q[1] = X1; Q[2] = Y1; Q[3] = Z1; return; }

  XX = Fl_sqr_pre(X1, p, pi);
  YY = Fl_sqr_pre(Y1, p, pi);
  YYYY = Fl_sqr_pre(YY, p, pi);
  ZZ = Fl_sqr_pre(Z1, p, pi);
  S = Fl_double(Fl_sub(Fl_sqr_pre(Fl_add(X1, YY, p), p, pi),
                       Fl_add(XX, YYYY, p), p), p);
  M = Fl_addmul_pre(Fl_triple(XX, p), a4, Fl_sqr_pre(ZZ, p, pi), p, pi);
  T = Fl_sub(Fl_sqr_pre(M, p, pi), Fl_double(S, p), p);
  Q[1] = T;
  Q[2] = Fl_sub(Fl_mul_pre(M, Fl_sub(S, T, p), p, pi),
                Fl_double(Fl_double(Fl_double(YYYY, p), p), p), p);
  Q[3] = Fl_sub(Fl_sqr_pre(Fl_add(Y1, Z1, p), p, pi),
                Fl_add(YY, ZZ, p), p);
}

INLINE void
Flj_dbl_pre_inplace(GEN P, ulong a4, ulong p, ulong pi)
{
  Flj_dbl_indir_pre(P, P, a4, p, pi);
}

GEN
Flj_dbl_pre(GEN P, ulong a4, ulong p, ulong pi)
{
  GEN Q = cgetg(4, t_VECSMALL);
  Flj_dbl_indir_pre(P, Q, a4, p, pi);
  return Q;
}

/* Cost: 11M + 5S + 9add + 4*2.
 * http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl */
INLINE void
Flj_add_indir_pre(GEN P, GEN Q, GEN R, ulong a4, ulong p, ulong pi)
{
  ulong X1, Y1, Z1, X2, Y2, Z2;
  ulong Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, r, V, W;
  X1 = P[1]; Y1 = P[2]; Z1 = P[3];
  X2 = Q[1]; Y2 = Q[2]; Z2 = Q[3];

  if (Z2 == 0) { R[1] = X1; R[2] = Y1; R[3] = Z1; return; }
  if (Z1 == 0) { R[1] = X2; R[2] = Y2; R[3] = Z2; return; }

  Z1Z1 = Fl_sqr_pre(Z1, p, pi);
  Z2Z2 = Fl_sqr_pre(Z2, p, pi);
  U1 = Fl_mul_pre(X1, Z2Z2, p, pi);
  U2 = Fl_mul_pre(X2, Z1Z1, p, pi);
  S1 = Fl_mul_pre(Y1, Fl_mul_pre(Z2, Z2Z2, p, pi), p, pi);
  S2 = Fl_mul_pre(Y2, Fl_mul_pre(Z1, Z1Z1, p, pi), p, pi);
  H = Fl_sub(U2, U1, p);
  r = Fl_double(Fl_sub(S2, S1, p), p);

  if (H == 0) {
    if (r == 0) Flj_dbl_indir_pre(P, R, a4, p, pi); /* P = Q */
    else { R[1] = R[2] = 1; R[3] = 0; } /* P = -Q */
    return;
  }
  I = Fl_sqr_pre(Fl_double(H, p), p, pi);
  J = Fl_mul_pre(H, I, p, pi);
  V = Fl_mul_pre(U1, I, p, pi);
  W = Fl_sub(Fl_sqr_pre(r, p, pi), Fl_add(J, Fl_double(V, p), p), p);
  R[1] = W;
  R[2] = Fl_sub(Fl_mul_pre(r, Fl_sub(V, W, p), p, pi),
                Fl_double(Fl_mul_pre(S1, J, p, pi), p), p);
  R[3] = Fl_mul_pre(Fl_sub(Fl_sqr_pre(Fl_add(Z1, Z2, p), p, pi),
                           Fl_add(Z1Z1, Z2Z2, p), p), H, p, pi);
}

INLINE void
Flj_add_pre_inplace(GEN P, GEN Q, ulong a4, ulong p, ulong pi)
{ Flj_add_indir_pre(P, Q, P, a4, p, pi); }

GEN
Flj_add_pre(GEN P, GEN Q, ulong a4, ulong p, ulong pi)
{
  GEN R = cgetg(4, t_VECSMALL);
  Flj_add_indir_pre(P, Q, R, a4, p, pi);
  return R;
}

GEN
Flj_neg(GEN Q, ulong p)
{ return mkvecsmall3(Q[1], Fl_neg(Q[2], p), Q[3]); }

typedef struct {
  ulong pbits, nbits;  /* Positive bits and negative bits */
  ulong lnzb; /* Leading nonzero bit */
} naf_t;

/* Return the signed binary representation (i.e. the Non-Adjacent Form
 * in base 2) of a; a = x.pbits - x.nbits (+ 2^BILif < 0; this
 * exceptional case can happen if a > 2^(BIL-1)) */
static void
naf_repr(naf_t *x, ulong a)
{
  ulong pbits = 0, nbits = 0, c0 = 0, c1, a0;
  long t, i;

  for (i = 0; a; a >>= 1, ++i) {
    a0 = a & 1;
    c1 = (c0 + a0 + ((a & 2) >> 1)) >> 1;
    t = c0 + a0 - (c1 << 1);
    if (t < 0)
      nbits |= (1UL << i);
    else if (t > 0)
      pbits |= (1UL << i);
    c0 = c1;
  }
  c1 = c0 >> 1;
  t = c0 - (c1 << 1);
  /* since a >= 0, we have t >= 0; if i = BIL, pbits (virtually) overflows;
   * that leading overflowed bit is implied and not recorded in pbits */
  if (t > 0 && i != BITS_IN_LONG) pbits |= (1UL << i);
  x->pbits = pbits;
  x->nbits = nbits;
  x->lnzb = t? i-2: i-3;
}

/* Standard left-to-right signed double-and-add to compute [n]P. */
static GEN
Flj_mulu_pre_naf(GEN P, ulong n, ulong a4, ulong p, ulong pi, const naf_t *x)
{
  GEN R, Pi;
  ulong pbits, nbits;
  ulong m;

  if (n == 0) return mkvecsmall3(1, 1, 0);
  if (n == 1) return Flv_copy(P);

  R = Flj_dbl_pre(P, a4, p, pi);
  if (n == 2) return R;

  pbits = x->pbits;
  nbits = x->nbits;
  Pi = nbits? Flj_neg(P, p): NULL;
  m = (1UL << x->lnzb);
  for ( ; m; m >>= 1) {
    Flj_dbl_pre_inplace(R, a4, p, pi);
    if (m & pbits)
      Flj_add_pre_inplace(R, P, a4, p, pi);
    else if (m & nbits)
      Flj_add_pre_inplace(R, Pi, a4, p, pi);
  }
  return gc_const((pari_sp)R, R);
}

GEN
Flj_mulu_pre(GEN P, ulong n, ulong a4, ulong p, ulong pi)
{
  naf_t x; naf_repr(&x, n);
  return Flj_mulu_pre_naf(P, n, a4, p, pi, &x);
}

struct _Flj { ulong a4, p, pi; };

static GEN
_Flj_add(void *E, GEN P, GEN Q)
{
  struct _Flj *ell=(struct _Flj *) E;
  return Flj_add_pre(P, Q, ell->a4, ell->p, ell->pi);
}

static GEN
_Flj_mul(void *E, GEN P, GEN n)
{
  struct _Flj *ell = (struct _Flj *) E;
  long s = signe(n);
  GEN Q;
  if (s==0) return mkvecsmall3(1, 1, 0);
  Q = Flj_mulu_pre(P, itou(n), ell->a4, ell->p, ell->pi);
  return s>0 ? Q : Flj_neg(Q, ell->p);
}
static GEN
_Flj_one(void *E)
{ (void) E; return mkvecsmall3(1, 1, 0); }

GEN
FljV_factorback_pre(GEN P, GEN L, ulong a4, ulong p, ulong pi)
{
  struct _Flj E;
  E.a4 = a4; E.p = p; E.pi = pi;
  return gen_factorback(P, L, (void*)&E, &_Flj_add, &_Flj_mul, &_Flj_one);
}

ulong
Flj_order_ufact(GEN P, ulong n, GEN fa, ulong a4, ulong p, ulong pi)
{
  pari_sp av = avma;
  GEN T = gel(fa,1), E = gel(fa,2);
  long i, l = lg(T);
  ulong res = 1;

  for (i = 1; i < l; i++, set_avma(av))
  {
    ulong j, t = T[i], e = E[i];
    GEN b = P;
    naf_t x; naf_repr(&x, t);
    if (l != 2) b = Flj_mulu_pre(b, n / upowuu(t,e), a4, p, pi);
    for (j = 0; j < e && b[3]; j++) b = Flj_mulu_pre_naf(b, t, a4, p, pi, &x);
    if (b[3]) return 0;
    res *= upowuu(t, j);
  }
  return res;
}

GEN
Fle_to_Flj(GEN P)
{ return ell_is_inf(P) ? mkvecsmall3(1UL, 1UL, 0UL):
                         mkvecsmall3(P[1], P[2], 1UL); }

GEN
Flj_to_Fle(GEN P, ulong p)
{
  if (P[3] == 0) return ellinf();
  else
  {
    ulong Z = Fl_inv(P[3], p);
    ulong Z2 = Fl_sqr(Z, p);
    ulong X3 = Fl_mul(P[1], Z2, p);
    ulong Y3 = Fl_mul(P[2], Fl_mul(Z, Z2, p), p);
    return mkvecsmall2(X3, Y3);
  }
}

GEN
Flj_to_Fle_pre(GEN P, ulong p, ulong pi)
{
  if (P[3] == 0) return ellinf();
  else
  {
    ulong Z = Fl_inv(P[3], p);
    ulong Z2 = Fl_sqr_pre(Z, p, pi);
    ulong X3 = Fl_mul_pre(P[1], Z2, p, pi);
    ulong Y3 = Fl_mul_pre(P[2], Fl_mul_pre(Z, Z2, p, pi), p, pi);
    return mkvecsmall2(X3, Y3);
  }
}

INLINE void
random_Fle_pre_indir(ulong a4, ulong a6, ulong p, ulong pi,
                     ulong *pt_x, ulong *pt_y)
{
  ulong x, x2, y, rhs;
  do
  {
    x   = random_Fl(p); /*  x^3+a4*x+a6 = x*(x^2+a4)+a6  */
    x2  = Fl_sqr_pre(x, p, pi);
    rhs = Fl_addmul_pre(a6, x, Fl_add(x2, a4, p), p, pi);
  } while ((!rhs && !Fl_add(Fl_triple(x2,p),a4,p)) || krouu(rhs, p) < 0);
  y = Fl_sqrt_pre(rhs, p, pi);
  *pt_x = x; *pt_y = y;
}

GEN
random_Flj_pre(ulong a4, ulong a6, ulong p, ulong pi)
{
  ulong x, y;
  random_Fle_pre_indir(a4, a6, p, pi, &x, &y);
  return mkvecsmall3(x, y, 1);
}

GEN
Flj_changepointinv_pre(GEN P, GEN ch, ulong p, ulong pi)
{
  ulong c, u, r, s, t, u2, u3, z2;
  ulong x  = uel(P,1), y = uel(P,2), z = uel(P,3);
  GEN w;
  if (z == 0) return Flv_copy(P);
  u = ch[1]; r = ch[2];
  s = ch[3]; t = ch[4];
  u2 = Fl_sqr_pre(u, p, pi); u3 = Fl_mul_pre(u, u2, p, pi);
  c = Fl_mul_pre(u2, x, p, pi);
  z2 = Fl_sqr_pre(z, p, pi);
  w = cgetg(4, t_VECSMALL);
  uel(w,1) = Fl_add(c, Fl_mul_pre(r, z2, p, pi), p);
  uel(w,2) = Fl_add(Fl_mul_pre(u3 ,y, p, pi),
                    Fl_mul_pre(z, Fl_add(Fl_mul_pre(s,c,p,pi),
                                         Fl_mul_pre(z2,t,p,pi), p), p, pi), p);
  uel(w,3) = z;
  return w;
}

/***********************************************************************/
/**                              Fle                                  **/
/***********************************************************************/
GEN
Fle_changepoint(GEN P, GEN ch, ulong p)
{
  ulong c, u, r, s, t, v, v2, v3;
  GEN z;
  if (ell_is_inf(P)) return ellinf();
  u = ch[1]; r = ch[2];
  s = ch[3]; t = ch[4];
  v = Fl_inv(u, p); v2 = Fl_sqr(v,p); v3 = Fl_mul(v,v2,p);
  c = Fl_sub(uel(P,1),r,p);
  z = cgetg(3,t_VECSMALL);
  z[1] = Fl_mul(v2, c, p);
  z[2] = Fl_mul(v3, Fl_sub(uel(P,2), Fl_add(Fl_mul(s,c, p),t, p),p),p);
  return z;
}

GEN
Fle_changepointinv(GEN P, GEN ch, ulong p)
{
  ulong c, u, r, s, t, u2, u3;
  GEN z;
  if (ell_is_inf(P)) return ellinf();
  u = ch[1]; r = ch[2];
  s = ch[3]; t = ch[4];
  u2 = Fl_sqr(u, p); u3 = Fl_mul(u,u2,p);
  c = Fl_mul(u2,uel(P,1), p);
  z = cgetg(3, t_VECSMALL);
  z[1] = Fl_add(c,r,p);
  z[2] = Fl_add(Fl_mul(u3,uel(P,2),p), Fl_add(Fl_mul(s,c,p), t, p), p);
  return z;
}
static GEN
Fle_dbl_slope(GEN P, ulong a4, ulong p, ulong *slope)
{
  ulong x, y, Qx, Qy;
  if (ell_is_inf(P) || !P[2]) return ellinf();
  x = P[1]; y = P[2];
  *slope = Fl_div(Fl_add(Fl_triple(Fl_sqr(x,p), p), a4, p),
                  Fl_double(y, p), p);
  Qx = Fl_sub(Fl_sqr(*slope, p), Fl_double(x, p), p);
  Qy = Fl_sub(Fl_mul(*slope, Fl_sub(x, Qx, p), p), y, p);
  return mkvecsmall2(Qx, Qy);
}

GEN
Fle_dbl(GEN P, ulong a4, ulong p)
{
  ulong slope;
  return Fle_dbl_slope(P,a4,p,&slope);
}

static GEN
Fle_add_slope(GEN P, GEN Q, ulong a4, ulong p, ulong *slope)
{
  ulong Px, Py, Qx, Qy, Rx, Ry;
  if (ell_is_inf(P)) return Q;
  if (ell_is_inf(Q)) return P;
  Px = P[1]; Py = P[2];
  Qx = Q[1]; Qy = Q[2];
  if (Px==Qx) return Py==Qy ? Fle_dbl_slope(P, a4, p, slope): ellinf();
  *slope = Fl_div(Fl_sub(Py, Qy, p), Fl_sub(Px, Qx, p), p);
  Rx = Fl_sub(Fl_sub(Fl_sqr(*slope, p), Px, p), Qx, p);
  Ry = Fl_sub(Fl_mul(*slope, Fl_sub(Px, Rx, p), p), Py, p);
  return mkvecsmall2(Rx, Ry);
}

GEN
Fle_add(GEN P, GEN Q, ulong a4, ulong p)
{
  ulong slope;
  return Fle_add_slope(P,Q,a4,p,&slope);
}

static GEN
Fle_neg(GEN P, ulong p)
{
  if (ell_is_inf(P)) return P;
  return mkvecsmall2(P[1], Fl_neg(P[2], p));
}

GEN
Fle_sub(GEN P, GEN Q, ulong a4, ulong p)
{
  pari_sp av = avma;
  ulong slope;
  return gerepileupto(av, Fle_add_slope(P, Fle_neg(Q, p), a4, p, &slope));
}

struct _Fle { ulong a4, a6, p; };

static GEN
_Fle_dbl(void *E, GEN P)
{
  struct _Fle *ell = (struct _Fle *) E;
  return Fle_dbl(P, ell->a4, ell->p);
}

static GEN
_Fle_add(void *E, GEN P, GEN Q)
{
  struct _Fle *ell=(struct _Fle *) E;
  return Fle_add(P, Q, ell->a4, ell->p);
}

GEN
Fle_mulu(GEN P, ulong n, ulong a4, ulong p)
{
  ulong pi;
  if (!n || ell_is_inf(P)) return ellinf();
  if (n==1) return zv_copy(P);
  if (n==2) return Fle_dbl(P, a4, p);
  pi = get_Fl_red(p);
  return Flj_to_Fle_pre(Flj_mulu_pre(Fle_to_Flj(P), n, a4, p, pi), p, pi);
}

static GEN
_Fle_mul(void *E, GEN P, GEN n)
{
  pari_sp av = avma;
  struct _Fle *e=(struct _Fle *) E;
  long s = signe(n);
  GEN Q;
  if (!s || ell_is_inf(P)) return ellinf();
  if (s < 0) P = Fle_neg(P, e->p);
  if (is_pm1(n)) return s > 0? zv_copy(P): P;
  Q = (lgefint(n)==3) ? Fle_mulu(P, uel(n,2), e->a4, e->p):
                        gen_pow(P, n, (void*)e, &_Fle_dbl, &_Fle_add);
  return s > 0? Q: gerepileuptoleaf(av, Q);
}

GEN
Fle_mul(GEN P, GEN n, ulong a4, ulong p)
{
  struct _Fle E;
  E.a4 = a4; E.p = p;
  return _Fle_mul(&E, P, n);
}

/* Finds a random nonsingular point on E */
GEN
random_Fle_pre(ulong a4, ulong a6, ulong p, ulong pi)
{
  ulong x, y;
  random_Fle_pre_indir(a4, a6, p, pi, &x, &y);
  return mkvecsmall2(x, y);
}

GEN
random_Fle(ulong a4, ulong a6, ulong p)
{ return random_Fle_pre(a4, a6, p, get_Fl_red(p)); }

static GEN
_Fle_rand(void *E)
{
  struct _Fle *e=(struct _Fle *) E;
  return random_Fle(e->a4, e->a6, e->p);
}

static const struct bb_group Fle_group={_Fle_add,_Fle_mul,_Fle_rand,hash_GEN,zv_equal,ell_is_inf,NULL};

GEN
Fle_order(GEN z, GEN o, ulong a4, ulong p)
{
  pari_sp av = avma;
  struct _Fle e;
  e.a4=a4;
  e.p=p;
  return gerepileuptoint(av, gen_order(z, o, (void*)&e, &Fle_group));
}

GEN
Fle_log(GEN a, GEN b, GEN o, ulong a4, ulong p)
{
  pari_sp av = avma;
  struct _Fle e;
  e.a4=a4;
  e.p=p;
  return gerepileuptoint(av, gen_PH_log(a, b, o, (void*)&e, &Fle_group));
}

ulong
Fl_ellj(ulong a4, ulong a6, ulong p)
{
  if (SMALL_ULONG(p))
  { /* a43 = 4 a4^3 */
    ulong a43 = Fl_double(Fl_double(Fl_mul(a4, Fl_sqr(a4, p), p), p), p);
    /* a62 = 27 a6^2 */
    ulong a62 = Fl_mul(Fl_sqr(a6, p), 27 % p, p);
    ulong z1 = Fl_mul(a43, 1728 % p, p);
    ulong z2 = Fl_add(a43, a62, p);
    return Fl_div(z1, z2, p);
  }
  return Fl_ellj_pre(a4, a6, p, get_Fl_red(p));
}

void
Fl_ellj_to_a4a6(ulong j, ulong p, ulong *pt_a4, ulong *pt_a6)
{
  ulong zagier = 1728 % p;
  if (j == 0)           { *pt_a4 = 0; *pt_a6 =1; }
  else if (j == zagier) { *pt_a4 = 1; *pt_a6 =0; }
  else
  {
    ulong k = Fl_sub(zagier, j, p);
    ulong kj = Fl_mul(k, j, p);
    ulong k2j = Fl_mul(kj, k, p);
    *pt_a4 = Fl_triple(kj, p);
    *pt_a6 = Fl_double(k2j, p);
  }
}

ulong
Fl_elldisc_pre(ulong a4, ulong a6, ulong p, ulong pi)
{ /* D = -(4A^3 + 27B^2) */
  ulong t1, t2;
  t1 = Fl_mul_pre(a4, Fl_sqr_pre(a4, p, pi), p, pi);
  t1 = Fl_double(Fl_double(t1, p), p);
  t2 = Fl_mul_pre(27 % p, Fl_sqr_pre(a6, p, pi), p, pi);
  return Fl_neg(Fl_add(t1, t2, p), p);
}

ulong
Fl_elldisc(ulong a4, ulong a6, ulong p)
{
  if (SMALL_ULONG(p))
  { /* D = -(4A^3 + 27B^2) */
    ulong t1, t2;
    t1 = Fl_mul(a4, Fl_sqr(a4, p), p);
    t1 = Fl_double(Fl_double(t1, p), p);
    t2 = Fl_mul(27 % p, Fl_sqr(a6, p), p);
    return Fl_neg(Fl_add(t1, t2, p), p);
  }
  return Fl_elldisc_pre(a4, a6, p, get_Fl_red(p));
}

void
Fl_elltwist_disc(ulong a4, ulong a6, ulong D, ulong p, ulong *pa4, ulong *pa6)
{
  ulong D2 = Fl_sqr(D, p);
  *pa4 = Fl_mul(a4, D2, p);
  *pa6 = Fl_mul(a6, Fl_mul(D, D2, p), p);
}

void
Fl_elltwist(ulong a4, ulong a6, ulong p, ulong *pt_a4, ulong *pt_a6)
{ Fl_elltwist_disc(a4, a6, nonsquare_Fl(p), p, pt_a4, pt_a6); }

static void
Fle_dbl_sinv_pre_inplace(GEN P, ulong a4, ulong sinv, ulong p, ulong pi)
{
  ulong x, y, slope;
  if (uel(P,1)==p) return;
  if (!P[2]) { P[1] = p; return; }
  x = P[1]; y = P[2];
  slope = Fl_mul_pre(Fl_add(Fl_triple(Fl_sqr_pre(x, p, pi), p), a4, p),
                sinv, p, pi);
  P[1] = Fl_sub(Fl_sqr_pre(slope, p, pi), Fl_double(x, p), p);
  P[2] = Fl_sub(Fl_mul_pre(slope, Fl_sub(x, P[1], p), p, pi), y, p);
}

static void
Fle_add_sinv_pre_inplace(GEN P, GEN Q, ulong a4, ulong sinv, ulong p, ulong pi)
{
  ulong Px, Py, Qx, Qy, slope;
  if (uel(P,1)==p) { P[1] = Q[1]; P[2] = Q[2]; }
  if (ell_is_inf(Q)) return;
  Px = P[1]; Py = P[2];
  Qx = Q[1]; Qy = Q[2];
  if (Px==Qx)
  {
    if (Py==Qy) Fle_dbl_sinv_pre_inplace(P, a4, sinv, p, pi);
    else P[1] = p;
    return;
  }
  slope = Fl_mul_pre(Fl_sub(Py, Qy, p), sinv, p, pi);
  P[1] = Fl_sub(Fl_sub(Fl_sqr_pre(slope, p, pi), Px, p), Qx, p);
  P[2] = Fl_sub(Fl_mul_pre(slope, Fl_sub(Px, P[1], p), p, pi), Py, p);
}

static void
Fle_sub_sinv_pre_inplace(GEN P, GEN Q, ulong a4, ulong sinv, ulong p, ulong pi)
{
  ulong Px, Py, Qx, Qy, slope;
  if (uel(P,1)==p) { P[1] = Q[1]; P[2] = Fl_neg(Q[2], p); }
  if (ell_is_inf(Q)) return;
  Px = P[1]; Py = P[2];
  Qx = Q[1]; Qy = Q[2];
  if (Px==Qx)
  {
    if (Py==Qy) P[1] = p;
    else
      Fle_dbl_sinv_pre_inplace(P, a4, sinv, p, pi);
    return;
  }
  slope = Fl_mul_pre(Fl_add(Py, Qy, p), sinv, p, pi);
  P[1] = Fl_sub(Fl_sub(Fl_sqr_pre(slope, p, pi), Px, p), Qx, p);
  P[2] = Fl_sub(Fl_mul_pre(slope, Fl_sub(Px, P[1], p), p, pi), Py, p);
}

static long
skipzero(long n) { return n ? n:1; }

void
FleV_add_pre_inplace(GEN P, GEN Q, GEN a4, ulong p, ulong pi)
{
  long i, l=lg(a4);
  GEN sinv = cgetg(l, t_VECSMALL);
  for(i=1; i<l; i++)
    uel(sinv,i) = umael(P,i,1)==p? 1: skipzero(Fl_sub(mael(P,i,1), mael(Q,i,1), p));
  Flv_inv_pre_inplace(sinv, p, pi);
  for (i=1; i<l; i++)
    Fle_add_sinv_pre_inplace(gel(P,i), gel(Q,i), uel(a4,i), uel(sinv,i), p, pi);
}

void
FleV_sub_pre_inplace(GEN P, GEN Q, GEN a4, ulong p, ulong pi)
{
  long i, l=lg(a4);
  GEN sinv = cgetg(l, t_VECSMALL);
  for(i=1; i<l; i++)
    uel(sinv,i) = umael(P,i,1)==p? 1: skipzero(Fl_sub(mael(P,i,1), mael(Q,i,1), p));
  Flv_inv_pre_inplace(sinv, p, pi);
  for (i=1; i<l; i++)
    Fle_sub_sinv_pre_inplace(gel(P,i), gel(Q,i), uel(a4,i), uel(sinv,i), p, pi);
}

void
FleV_dbl_pre_inplace(GEN P, GEN a4, ulong p, ulong pi)
{
  long i, l=lg(a4);
  GEN sinv = cgetg(l, t_VECSMALL);
  for(i=1; i<l; i++)
    uel(sinv,i) = umael(P,i,1)==p? 1: skipzero(Fl_double(umael(P,i,2), p));
  Flv_inv_pre_inplace(sinv, p, pi);
  for(i=1; i<l; i++)
    Fle_dbl_sinv_pre_inplace(gel(P,i), uel(a4,i), uel(sinv,i), p, pi);
}

static void
FleV_mulu_pre_naf_inplace(GEN P, ulong n, GEN a4, ulong p, ulong pi, const naf_t *x)
{
  pari_sp av = avma;
  ulong pbits, nbits, m;
  GEN R;
  if (n == 1) return;

  R = P; P = gcopy(P);
  FleV_dbl_pre_inplace(R, a4, p, pi);
  if (n == 2) return;

  pbits = x->pbits;
  nbits = x->nbits;
  m = (1UL << x->lnzb);
  for ( ; m; m >>= 1) {
    FleV_dbl_pre_inplace(R, a4, p, pi);
    if (m & pbits)
      FleV_add_pre_inplace(R, P, a4, p, pi);
    else if (m & nbits)
      FleV_sub_pre_inplace(R, P, a4, p, pi);
  }
  set_avma(av);
}

void
FleV_mulu_pre_inplace(GEN P, ulong n, GEN a4, ulong p, ulong pi)
{
  naf_t x; naf_repr(&x, n);
  FleV_mulu_pre_naf_inplace(P, n, a4, p, pi, &x);
}

/***********************************************************************/
/**                            Pairings                               **/
/**               Derived from APIP by Jerome Milan, 2012             **/
/***********************************************************************/
static ulong
Fle_vert(GEN P, GEN Q, ulong a4, ulong p)
{
  if (ell_is_inf(P))
    return 1;
  if (uel(Q, 1) != uel(P, 1))
    return Fl_sub(uel(Q, 1), uel(P, 1), p);
  if (uel(P,2) != 0 ) return 1;
  return Fl_inv(Fl_add(Fl_triple(Fl_sqr(uel(P,1),p), p), a4, p), p);
}

static ulong
Fle_Miller_line(GEN R, GEN Q, ulong slope, ulong a4, ulong p)
{
  ulong x = uel(Q, 1), y = uel(Q, 2);
  ulong tmp1 = Fl_sub(x, uel(R, 1), p);
  ulong tmp2 = Fl_add(Fl_mul(tmp1, slope, p), uel(R,2), p);
  if (y != tmp2)
    return Fl_sub(y, tmp2, p);
  if (y == 0)
    return 1;
  else
  {
    ulong s1, s2;
    ulong y2i = Fl_inv(Fl_double(y, p), p);
    s1 = Fl_mul(Fl_add(Fl_triple(Fl_sqr(x, p), p), a4, p), y2i, p);
    if (s1 != slope)
      return Fl_sub(s1, slope, p);
    s2 = Fl_mul(Fl_sub(Fl_triple(x, p), Fl_sqr(s1, p), p), y2i, p);
    return s2 ? s2: y2i;
  }
}

/* Computes the equation of the line tangent to R and returns its
 * evaluation at the point Q. Also doubles the point R. */
static ulong
Fle_tangent_update(GEN R, GEN Q, ulong a4, ulong p, GEN *pt_R)
{
  if (ell_is_inf(R)) { *pt_R = ellinf(); return 1; }
  else if (uel(R,2) == 0) { *pt_R = ellinf(); return Fle_vert(R, Q, a4, p); }
  else
  {
    ulong slope;
    *pt_R = Fle_dbl_slope(R, a4, p, &slope);
    return Fle_Miller_line(R, Q, slope, a4, p);
  }
}

/* Computes the equation of the line through R and P, and returns its
 * evaluation at the point Q. Also adds P to the point R. */
static ulong
Fle_chord_update(GEN R, GEN P, GEN Q, ulong a4, ulong p, GEN *pt_R)
{
  if (ell_is_inf(R)) { *pt_R = P; return Fle_vert(P, Q, a4, p); }
  else if (ell_is_inf(P)) { *pt_R = R; return Fle_vert(R, Q, a4, p); }
  else if (uel(P, 1) == uel(R, 1))
  {
    if (uel(P, 2) == uel(R, 2)) return Fle_tangent_update(R, Q, a4, p, pt_R);
    else { *pt_R = ellinf(); return Fle_vert(R, Q, a4, p); }
  }
  else
  {
    ulong slope;
    *pt_R = Fle_add_slope(P, R, a4, p, &slope);
    return Fle_Miller_line(R, Q, slope, a4, p);
  }
}

struct _Fle_miller { ulong p, a4; GEN P; };
static GEN
Fle_Miller_dbl(void* E, GEN d)
{
  struct _Fle_miller *m = (struct _Fle_miller *)E;
  ulong p = m->p, a4 = m->a4;
  GEN P = m->P;
  ulong v, line;
  ulong N = Fl_sqr(umael(d,1,1), p);
  ulong D = Fl_sqr(umael(d,1,2), p);
  GEN point = gel(d,2);
  line = Fle_tangent_update(point, P, a4, p, &point);
  N  = Fl_mul(N, line, p);
  v = Fle_vert(point, P, a4, p);
  D = Fl_mul(D, v, p); return mkvec2(mkvecsmall2(N, D), point);
}
static GEN
Fle_Miller_add(void* E, GEN va, GEN vb)
{
  struct _Fle_miller *m = (struct _Fle_miller *)E;
  ulong p = m->p, a4= m->a4;
  GEN P = m->P, point;
  ulong v, line;
  ulong na = umael(va,1,1), da = umael(va,1,2);
  ulong nb = umael(vb,1,1), db = umael(vb,1,2);
  GEN pa = gel(va,2), pb = gel(vb,2);
  ulong N = Fl_mul(na, nb, p);
  ulong D = Fl_mul(da, db, p);
  line = Fle_chord_update(pa, pb, P, a4, p, &point);
  N = Fl_mul(N, line, p);
  v = Fle_vert(point, P, a4, p);
  D = Fl_mul(D, v, p); return mkvec2(mkvecsmall2(N, D), point);
}

/* Returns the Miller function f_{m, Q} evaluated at the point P using
 * the standard Miller algorithm. */
static ulong
Fle_Miller(GEN Q, GEN P, ulong m, ulong a4, ulong p)
{
  pari_sp av = avma;
  struct _Fle_miller d;
  GEN v;
  ulong N, D;

  d.a4 = a4; d.p = p; d.P = P;
  v = gen_powu_i(mkvec2(mkvecsmall2(1,1), Q), m, (void*)&d,
                Fle_Miller_dbl, Fle_Miller_add);
  N = umael(v,1,1); D = umael(v,1,2);
  return gc_ulong(av, Fl_div(N, D, p));
}

ulong
Fle_weilpairing(GEN P, GEN Q, ulong m, ulong a4, ulong p)
{
  pari_sp ltop = avma;
  ulong N, D, w;
  if (ell_is_inf(P) || ell_is_inf(Q) || zv_equal(P,Q)) return 1;
  N = Fle_Miller(P, Q, m, a4, p);
  D = Fle_Miller(Q, P, m, a4, p);
  w = Fl_div(N, D, p);
  if (odd(m)) w  = Fl_neg(w, p);
  return gc_ulong(ltop, w);
}

ulong
Fle_tatepairing(GEN P, GEN Q, ulong m, ulong a4, ulong p)
{
  if (ell_is_inf(P) || ell_is_inf(Q)) return 1;
  return Fle_Miller(P, Q, m, a4, p);
}

GEN
Fl_ellptors(ulong l, ulong N, ulong a4, ulong a6, ulong p)
{
  long v = z_lval(N, l);
  ulong N0, N1;
  GEN F;
  if (v==0) return cgetg(1,t_VEC);
  N0 = upowuu(l, v); N1 = N/N0;
  F = mkmat2(mkcols(l), mkcols(v));
  while(1)
  {
    pari_sp av2 = avma;
    GEN P, Q;
    ulong d, s, t;

    P = Fle_mulu(random_Fle(a4, a6, p), N1, a4, p);
    Q = Fle_mulu(random_Fle(a4, a6, p), N1, a4, p);
    s = itou(Fle_order(P, F, a4, p));
    t = itou(Fle_order(Q, F, a4, p));
    if (s < t) { swap(P,Q); lswap(s,t); }
    if (s == N0) retmkvec(Fle_mulu(P, s/l, a4, p));
    d = Fl_order(Fle_weilpairing(P, Q, s, a4, p), s, p);
    if (s*d == N0)
      retmkvec2(Fle_mulu(P, s/l, a4, p), Fle_mulu(Q, t/l, a4, p));
    set_avma(av2);
  }
}
