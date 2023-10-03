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

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_arith

/*********************************************************************/
/**                                                                 **/
/**                    FUNDAMENTAL DISCRIMINANTS                    **/
/**                                                                 **/
/*********************************************************************/
static long
fa_isfundamental(GEN F)
{
  GEN P = gel(F,1), E = gel(F,2);
  long i, s, l = lg(P);

  if (l == 1) return 1;
  s = signe(gel(P,1)); /* = signe(x) */
  if (!s) return 0;
  if (s < 0) { l--; P = vecslice(P,2,l); E = vecslice(E,2,l); }
  if (l == 1) return 0;
  if (!absequaliu(gel(P,1), 2))
    i = 1; /* need x = 1 mod 4 */
  else
  {
    i = 2;
    switch(itou(gel(E,1)))
    {
      case 2: s = -s; break; /* need x/4 = 3 mod 4 */
      case 3: s = 0; break; /* no condition mod 4 */
      default: return 0;
    }
  }
  for(; i < l; i++)
  {
    if (!equali1(gel(E,i))) return 0;
    if (s && Mod4(gel(P,i)) == 3) s = -s;
  }
  return s >= 0;
}
long
isfundamental(GEN x)
{
  if (typ(x) != t_INT)
  {
    pari_sp av = avma;
    long v = fa_isfundamental(check_arith_all(x,"isfundamental"));
    return gc_long(av,v);
  }
  return Z_isfundamental(x);
}

/* x fundamental ? */
long
uposisfundamental(ulong x)
{
  ulong r = x & 15; /* x mod 16 */
  if (!r) return 0;
  switch(r & 3)
  { /* x mod 4 */
    case 0: return (r == 4)? 0: uissquarefree(x >> 2);
    case 1: return uissquarefree(x);
    default: return 0;
  }
}
/* -x fundamental ? */
long
unegisfundamental(ulong x)
{
  ulong r = x & 15; /* x mod 16 */
  if (!r) return 0;
  switch(r & 3)
  { /* x mod 4 */
    case 0: return (r == 12)? 0: uissquarefree(x >> 2);
    case 3: return uissquarefree(x);
    default: return 0;
  }
}
long
sisfundamental(long x)
{ return x < 0? unegisfundamental((ulong)(-x)): uposisfundamental(x); }

long
Z_isfundamental(GEN x)
{
  long r;
  switch(lgefint(x))
  {
    case 2: return 0;
    case 3: return signe(x) < 0? unegisfundamental(x[2])
                               : uposisfundamental(x[2]);
  }
  r = mod16(x);
  if (!r) return 0;
  if ((r & 3) == 0)
  {
    pari_sp av;
    r >>= 2; /* |x|/4 mod 4 */
    if (signe(x) < 0) r = 4-r;
    if (r == 1) return 0;
    av = avma;
    r = Z_issquarefree( shifti(x,-2) );
    return gc_long(av, r);
  }
  r &= 3; /* |x| mod 4 */
  if (signe(x) < 0) r = 4-r;
  return (r==1) ? Z_issquarefree(x) : 0;
}

static GEN
fa_quaddisc(GEN f)
{
  GEN P = gel(f,1), E = gel(f,2), s = gen_1;
  long i, l = lg(P);
  for (i = 1; i < l; i++) /* possibly including -1 */
    if (mpodd(gel(E,i))) s = mulii(s, gel(P,i));
  if (Mod4(s) > 1) s = shifti(s,2);
  return s;
}

GEN
quaddisc(GEN x)
{
  const pari_sp av = avma;
  long tx = typ(x);
  GEN F;
  if (is_rational_t(tx)) F = factor(x);
  else
  {
    F = check_arith_all(x,"quaddisc");
    if (tx == t_VEC && typ(gel(x,1)) == t_INT
                    && Z_issquarefree_fact(gel(x,2)))
    {
      x = gel(x,1);
      return (Mod4(x) > 1)? shifti(x, 2): icopy(x);
    }
  }
  return gerepileuptoint(av, fa_quaddisc(F));
}


/***********************************************************************/
/**                                                                   **/
/**         FUNDAMENTAL UNIT AND REGULATOR (QUADRATIC FIELDS)         **/
/**                                                                   **/
/***********************************************************************/
/* replace f by f * [u,1; 1,0] */
static void
update_f(GEN f, GEN u)
{
  GEN a = gcoeff(f,1,1), b = gcoeff(f,1,2);
  GEN c = gcoeff(f,2,1), d = gcoeff(f,2,2);
  gcoeff(f,1,1) = addmulii(b, u,a); gcoeff(f,1,2) = a;
  gcoeff(f,2,1) = addmulii(d, u,c); gcoeff(f,2,2) = c;
}

/* f is a vector of matrices and i an index whose bits give the non-zero
 * entries; the product of the non zero entries is the actual result.
 * if i odd, f[1] may be an int: implicitely represent [f[1],1;1,0] */
static long
update_fm(GEN f, GEN a, long i)
{
#ifdef LONG_IS_64BIT
  const long LIM = 10;
#else
  const long LIM = 18;
#endif
  pari_sp av = avma;
  long k, v;
  GEN u;
  if (!odd(i)) { gel(f,1) = a; return i+1; }
  u = gel(f, 1);
  if (typ(u) == t_INT) /* [u,1;1,0] * [a,1;1,0] */
  { gel(f,1) = mkmat22(addiu(mulii(a, u), 1), u, a, gen_1); return i; }
  update_f(u, a); if (lgefint(gcoeff(u,1,1)) < LIM) return i;
  v = vals(i+1); gel(f,1) = gen_0;
  for (k = 1; k < v; k++) { u = ZM2_mul(gel(f,k+1), u); gel(f,k+1) = gen_0; }
  if (v != 1) u = gerepileupto(av, u);
  gel(f,v+1) = u; return i+1;
}
/* \prod f[j]; if first only return column 1 */
static GEN
prod_fm(GEN f, long i, long first)
{
  long k, v = vals(i) + 1;
  GEN u = gel(f, v);
  /* i a power of 2: f[1] can't be a t_INT */
  if (!(i >>= v)) return first? gel(u,1): u;
  for (k = v+1; i; i >>= 1, k++)
    if (odd(i))
    {
      GEN w = gel(f,k);
      switch(typ(u))
      {
        case t_INT: update_f(w, u);
          u = first? gel(w,1): w; break;
        case t_COL: /* implies 'first' */
          u = ZM_ZC_mul(w, u); break;
        default: /* t_MAT */
          u = first? ZM_ZC_mul(w, gel(u,1)): ZM2_mul(w, u); break;
      }
    }
  return u;
}

GEN
quadunit0(GEN x, long v)
{
  GEN y = quadunit(x);
  if (v==-1) v = fetch_user_var("w");
  setvarn(gel(y,1), v); return y;
}

struct uimod { GEN N, T; };
static GEN
ui_pow(void *E, GEN x, GEN n)
{ struct uimod *S = (struct uimod*)E; return FpXQ_pow(x, n, S->T, S->N); }
static int
ui_equal1(GEN x) { return degpol(x) < 1; }
static const struct bb_group
ui_group={ NULL,ui_pow,NULL,NULL,NULL,ui_equal1,NULL};

static void
quadunit_uvmod(GEN D, GEN d, GEN N, GEN *pu, GEN *pv)
{
  GEN u1, u2, v1, v2, p, q, q1, u, v;
  int m = mpodd(D), first = 1;
  pari_sp av = avma;
  p = (mpodd(d) == m)? d: subiu(d, 1);
  u1 = negi(p); u2 = gen_2;
  v1 = gen_1; v2 = gen_0; q = gen_2;
  q1 = shifti(subii(D, sqri(p)), -1);
  for(;;)
  {
    GEN r, A = dvmdii(addii(p, d), q, &r), p1 = p, t;
    p = subii(d, r);
    if (equalii(p1, p) && !first)
    { /* even period */
      u = addmulii(sqri(u2), D, sqri(v2));
      v = shifti(mulii(u2,v2), 1);
      break;
    }
    first = 0;
    t = Fp_addmul(u1, A, u2, N); u1 = u2; u2 = t;
    t = Fp_addmul(v1, A, v2, N); v1 = v2; v2 = t;
    t = q; q = submulii(q1, A, subii(p, p1)); q1 = t;
    if (equalii(q, t))
    { /* odd period */
      u = addmulii(mulii(u1,u2), D, mulii(v1,v2));
      v = addmulii(mulii(u1,v2), u2, v1);
      break;
    }
    if (gc_needed(av, 2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"quadunit_uvmod");
      gerepileall(av, 7, &p, &u1,&u2,&v1,&v2, &q,&q1);
    }
  }
  *pu = modii(u, N);
  *pv = modii(v, N); if (m) *pu = Fp_sub(*pu, *pv, N);
}
/* fundamental unit is u + vx mod quadpoly(D); always called with D
 * fundamental and relatively small but would work in all cases. Should be
 * called whenever the fundamental unit is so "small" that asymptotically
 * fast multiplication is not used in the continued fraction loop */
static void
quadunit_uv_basecase(GEN D, GEN *pu, GEN *pv)
{
  GEN u1, u2, v1, v2, p, q, q1, u, v, a, b, c, d = sqrtremi(D, &a);
  int m = mpodd(D);
  long first = 1;

  p = d; q1 = shifti(a, -1); q = gen_2;
  if (mpodd(d) != m) { p = subiu(d,1); q1 = addii(q1,d); } /* q1 = (D-p^2)/2 */
  u1 = gen_2; u2 = p;
  v1 = gen_0; v2 = gen_1;
  for(;;)
  {
    GEN t = q;
    if (first) { first = 0; q = q1; }
    else
    {
      GEN r, A = dvmdii(addii(p, d), q, &r), p1 = p;
      p = subii(d, r);
      if (equalii(p1, p)) /* even period */
      { a = sqri(u2); b = sqri(v2); c = sqri(addii(u2, v2)); break; }
      r = addmulii(u1, A, u2); u1 = u2; u2 = r;
      r = addmulii(v1, A, v2); v1 = v2; v2 = r;
      q = submulii(q1, A, subii(p, p1));
    }
    q1 = t;
    if (equalii(q, t))
    { /* odd period */
      a = mulii(u1, u2); b = mulii(v1, v2);
      c = mulii(addii(u1, v1), addii(u2, v2)); break;
    }
  }
  u = diviiexact(addmulii(a, D, b), q);
  v = diviiexact(subii(c, addii(a, b)), q);
  if (m == 1) u = subii(u, v);
  *pu = shifti(u, -1); *pv = v;
}

/* D > 0, d = sqrti(D) */
static GEN
quadunit_q(GEN D, GEN d, long *pN)
{
  pari_sp av = avma;
  GEN p, q, q1;
  long first = 1;
  p = (Mod2(d) == Mod2(D))? d: subiu(d, 1);
  q = gen_2;
  q1 = shifti(subii(D, sqri(p)), -1);
  for(;;)
  {
    GEN r, A = dvmdii(addii(p, d), q, &r), p1 = p, t;
    p = subii(d, r);
    if (!first && equalii(p1, p)) { *pN = 1; return q; } /* even period */
    first = 0;
    t = q; q = submulii(q1, A, subii(p, p1)); q1 = t;
    if (equalii(q, t)) { *pN = -1; return q; } /* odd period */
    if (gc_needed(av, 2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"quadunitnorm");
      gerepileall(av, 3, &p, &q, &q1);
    }
  }
}
/* fundamental unit mod N */
static GEN
quadunit_mod(GEN D, GEN N)
{
  GEN q, u, v, d = sqrti(D);
  pari_sp av = avma;
  long s;
  q = gerepileuptoint(av, quadunit_q(D, d, &s));
  if (mpodd(N) && equali1(gcdii(q, N)))
  {
    quadunit_uvmod(D, d, N, &u, &v);
    q = Fp_inv(shifti(q, 1), N);
    u = Fp_mul(u, q, N);
    v = Fp_mul(v, q, N); v = modii(shifti(v, 1), N);
  }
  else
  {
    GEN M = shifti(mulii(q, N), 1);
    quadunit_uvmod(D, d, M, &u, &v);
    u = diviiexact(u, q);
    v = modii(diviiexact(v, q), N); u = shifti(u,-1);
  }
  return deg1pol_shallow(v, u, 0);
}

/* f \prod_{p|f}  [ 1 - (D/p) p^-1 ] = \prod_{p^e||f} p^(e-1) [ p - (D/p) ] */
static GEN
quadclassnoEuler_fact(GEN D, GEN P, GEN E)
{
  long i, l = lg(P);
  GEN H;
  if (typ(E) != t_VECSMALL) E = vec_to_vecsmall(E);
  for (i = 1, H = gen_1; i < l; i++)
  {
    GEN p = gel(P,i);
    long e = E[i], s = kronecker(D,p);
    if (!s)
      H = mulii(H, e == 1? p: powiu(p, e));
    else
    {
      H = mulii(H, subis(p, s));
      if (e >= 2) H = mulii(H, e == 2? p: powiu(p,e-1));
    }
  }
  return H;
}

/* D > 0; y mod (N,T) congruent to fundamental unit of maximal order and
 * disc D. Return unit index of order of conductor N */
static GEN
quadunitindex_ii(GEN D, GEN N, GEN F, GEN y, GEN T)
{
  GEN H = quadclassnoEuler_fact(D, gel(F,1), gel(F,2));
  GEN P, E, a = Z_smoothen(H, gel(F,1), &P, &E), faH = mkmat2(P, E);
  struct uimod S;

  if (a) faH = merge_factor(Z_factor(a), faH,(void*)&cmpii,cmp_nodata);
  /* multiple of unit index, in [H, factor(H)] format */
  S.N = N; S.T = FpX_red(T, N);
  return gen_order(y, mkvec2(H,faH), (void*)&S, &ui_group);
}
static GEN
quadunitindex_i(GEN D, GEN N, GEN F)
{ return quadunitindex_ii(D, N, F, quadunit_mod(D, N), quadpoly_i(D)); }
GEN
quadunitindex(GEN D, GEN N)
{
  pari_sp av = avma;
  long r, s;
  GEN F;
  check_quaddisc(D, &s, &r, "quadunitindex");
  if ((F = check_arith_pos(N,"quadunitindex")))
    N = typ(N) == t_VEC? gel(N,1): factorback(F);
  if (equali1(N)) return gen_1;
  if (s < 0) switch(itos_or_0(D)) {
    case -3: return utoipos(3);
    case -4: return utoipos(2);
    default: return gen_1;
  }
  return gerepileuptoint(av, quadunitindex_i(D, N, F? F: Z_factor(N)));
}

/* fundamental unit is u + vx mod quadpoly(D); always called with D
 * fundamental but would work in all cases. Same algorithm as basecase,
 * except we compute the product of elementary matrices with a product tree */
static void
quadunit_uv(GEN D, GEN *pu, GEN *pv)
{
  GEN a, b, c, u, v, p, q, q1, f, d = sqrtremi(D, &a);
  pari_sp av = avma;
  long i = 0;
  int m = mpodd(D);

  p = d; q1 = shifti(a, -1); q = gen_2;
  if (mpodd(d) != m) { p = subiu(d,1); q1 = addii(q1,d); } /* q1 = (D-p^2)/2 */
  f = zerovec(2 + (expi(D)>>1));
  gel(f,1) = mkmat22(p, gen_2, gen_1, gen_0);
  for(;;)
  {
    GEN t = q, u1,u2, v1,v2;
    if (!i) { i = 1; q = q1; }
    else
    {
      GEN r, A = dvmdii(addii(p, d), q, &r), p1 = p;
      p = subii(d, r);
      if (equalii(p1, p))
      { /* even period */
        f = prod_fm(f, i, 1); u2 = gel(f,1); v2 = gel(f,2);
        a = sqri(u2); b = sqri(v2); c = sqri(addii(u2, v2)); break;
      }
      i = update_fm(f, A, i);
      q = submulii(q1, A, subii(p, p1));
    }
    q1 = t;
    if (equalii(q, t))
    { /* odd period */
      f = prod_fm(f, i, 0);
      u2 = gcoeff(f,1,1); u1 = gcoeff(f,1,2); a = mulii(u1, u2);
      v2 = gcoeff(f,2,1); v1 = gcoeff(f,2,2); b = mulii(v1, v2);
      c = mulii(addii(u1, v1), addii(u2, v2)); break;
    }
    if (gc_needed(av, 2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"quadunit (%ld)", i);
      gerepileall(av, 4, &p, &f, &q,&q1);
    }
  }
  u = diviiexact(addmulii(a, D, b), q);
  v = diviiexact(subii(c, addii(a, b)), q);
  if (m == 1) u = subii(u, v);
  *pu = shifti(u, -1); *pv = v;
}
GEN
quadunit(GEN D0)
{
  pari_sp av = avma;
  GEN P, E, D, u, v;
  long s = signe(D0);
  /* check_quaddisc_real omitting test for squares */
  if (typ(D0) != t_INT) pari_err_TYPE("quadunit", D0);
  if (s <= 0) pari_err_DOMAIN("quadunit", "disc","<=",gen_0,D0);
  if (mod4(D0) > 1) pari_err_DOMAIN("quadunit","disc % 4",">", gen_1,D0);
  D = coredisc2_fact(Z_factor(D0), s, &P, &E);
  /* test for squares done here for free */
  if (equali1(D)) pari_err_DOMAIN("quadunit","issquare(disc)","=", gen_1,D0);
  if (cmpiu(D, 2000000) < 0)
    quadunit_uv_basecase(D, &u, &v);
  else
    quadunit_uv(D, &u, &v);
  if (lg(P) != 1)
  { /* non-trivial conductor N > 1 */
    GEN N = factorback2(P,E), qD = quadpoly_i(D);
    GEN n, y = deg1pol_shallow(v, u, 0); /* maximal order fund unit */
    n = quadunitindex_ii(D, N, mkvec2(P,E), FpX_red(y,N), qD); /* unit index */
    y = ZXQ_powu(y, itou(n), qD); /* fund unit of order of conductor N */
    v = gel(y,3); u = gel(y,2); /* u + v w_D */
    if (mpodd(D))
    { /* w_D = (1+sqrt(D))/2 */
      if (mpodd(D0))
      { /* w_D0 = (1 + N sqrt(D)) / 2 */
        GEN v0 = v;
        v = diviiexact(v, N);
        u = addii(u, shifti(subii(v0, v), -1));
      }
      else
      { /* w_D0 = N sqrt(D)/2, N is even */
        v = shifti(v, -1);
        u = addii(u, v);
        v = diviiexact(v, shifti(N,-1));
      }
    }
    else /* w_D = sqrt(D), w_D0 = N sqrt(D) */
      v = diviiexact(v, N);
  }
  return gerepilecopy(av, mkquad(quadpoly_i(D0), u, v));
}
long
quadunitnorm(GEN D)
{
  pari_sp av = avma;
  long s, r;
  check_quaddisc(D, &s, &r, "quadunitnorm");
  if (s < 0) return 1;
  (void)quadunit_q(D, sqrti(D), &s); return gc_long(av, s);
}

GEN
quadregulator(GEN x, long prec)
{
  pari_sp av = avma, av2;
  GEN R, rsqd, u, v, sqd;
  long r, e;

  check_quaddisc_real(x, &r, "quadregulator");
  sqd = sqrti(x);
  rsqd = gsqrt(x,prec); av2 = avma;
  e = 0; R = real2n(1, prec); u = utoi(r); v = gen_2;
  for(;;)
  {
    GEN u1 = subii(mulii(divii(addii(u,sqd),v), v), u);
    GEN v1 = divii(subii(x,sqri(u1)),v);
    if (equalii(v,v1)) { R = mulrr(sqrr(R), divri(addir(u1,rsqd),v)); break; }
    if (equalii(u,u1)) { R = sqrr(R); break; }
    R = mulrr(R, divri(addir(u1,rsqd),v));
    e += expo(R); setexpo(R,0);
    u = u1; v = v1;
    if (e & ~EXPOBITS) pari_err_OVERFLOW("quadregulator [exponent]");
    if (gc_needed(av2,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"quadregulator");
      gerepileall(av2,3, &R,&u,&v);
    }
  }
  R = divri(R, v); e = 2*e - 1;
  /* avoid loss of accuracy */
  if (!((e + expo(R)) & ~EXPOBITS)) { setexpo(R, e + expo(R)); e = 0; }
  R = logr_abs(R);
  if (e) R = addrr(R, mulsr(e, mplog2(prec)));
  return gerepileuptoleaf(av, R);
}

/*************************************************************************/
/**                                                                     **/
/**                            CLASS NUMBER                             **/
/**                                                                     **/
/*************************************************************************/

int
qfb_equal1(GEN f) { return equali1(gel(f,1)); }

static GEN qfi_pow(void *E, GEN f, GEN n)
{ return E? nupow(f,n,(GEN)E): qfbpow_i(f,n); }
static GEN qfi_comp(void *E, GEN f, GEN g)
{ return E? nucomp(f,g,(GEN)E): qfbcomp_i(f,g); }
static const struct bb_group qfi_group={ qfi_comp,qfi_pow,NULL,hash_GEN,
                                         gidentical,qfb_equal1,NULL};

GEN
qfi_order(GEN q, GEN o)
{ return gen_order(q, o, NULL, &qfi_group); }

GEN
qfi_log(GEN a, GEN g, GEN o)
{ return gen_PH_log(a, g, o, NULL, &qfi_group); }

GEN
qfi_Shanks(GEN a, GEN g, long n)
{
  pari_sp av = avma;
  GEN T, X;
  long rt_n, c;

  a = qfbred_i(a);
  g = qfbred_i(g);

  rt_n = sqrt((double)n);
  c = n / rt_n;
  c = (c * rt_n < n + 1) ? c + 1 : c;

  T = gen_Shanks_init(g, rt_n, NULL, &qfi_group);
  X = gen_Shanks(T, a, c, NULL, &qfi_group);
  return X? gerepileuptoint(av, X): gc_NULL(av);
}

GEN
qfbclassno0(GEN x,long flag)
{
  switch(flag)
  {
    case 0: return map_proto_G(classno,x);
    case 1: return map_proto_G(classno2,x);
    default: pari_err_FLAG("qfbclassno");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/* f^h = 1, return order(f). Set *pfao to its factorization */
static GEN
find_order(void *E, GEN f, GEN h, GEN *pfao)
{
  GEN v = gen_factored_order(f, h, E, &qfi_group);
  *pfao = gel(v,2); return gel(v,1);
}

static int
ok_q(GEN q, GEN h, GEN d2, long r2)
{
  if (d2)
  {
    if (r2 <= 2 && !mpodd(q)) return 0;
    return is_pm1(Z_ppo(q,d2));
  }
  else
  {
    if (r2 <= 1 && !mpodd(q)) return 0;
    return is_pm1(Z_ppo(q,h));
  }
}

/* a,b given by their factorizations. Return factorization of lcm(a,b).
 * Set A,B such that A*B = lcm(a, b), (A,B)=1, A|a, B|b */
static GEN
split_lcm(GEN a, GEN Fa, GEN b, GEN Fb, GEN *pA, GEN *pB)
{
  GEN P = ZC_union_shallow(gel(Fa,1), gel(Fb,1));
  GEN A = gen_1, B = gen_1;
  long i, l = lg(P);
  GEN E = cgetg(l, t_COL);
  for (i=1; i<l; i++)
  {
    GEN p = gel(P,i);
    long va = Z_pval(a,p);
    long vb = Z_pval(b,p);
    if (va < vb)
    {
      B = mulii(B,powiu(p,vb));
      gel(E,i) = utoi(vb);
    }
    else
    {
      A = mulii(A,powiu(p,va));
      gel(E,i) = utoi(va);
    }
  }
  *pA = A;
  *pB = B; return mkmat2(P,E);
}

/* g1 has order d1, f has order o, replace g1 by an element of order lcm(d1,o)*/
static void
update_g1(GEN *pg1, GEN *pd1, GEN *pfad1, GEN f, GEN o, GEN fao)
{
  GEN A,B, g1 = *pg1, d1 = *pd1;
  *pfad1 = split_lcm(d1,*pfad1, o,fao, &A,&B);
  *pg1 = gmul(qfbpow_i(g1, diviiexact(d1,A)),  qfbpow_i(f, diviiexact(o,B)));
  *pd1 = mulii(A,B); /* g1 has order d1 <- lcm(d1,o) */
}

/* Let s = 1 or -1; D = s * d; assume Df^2 fits in an ulong
 * Return  f / [O_{Df^2}^*:O_D^*] * \prod_{p|f}  [ 1 - (D/p) p^-1 ]
 * The Euler product is \prod_{p^e||f} p^(e-1) [ p - (D/p) ] */
ulong
uquadclassnoF_fact(ulong d, long s, GEN P, GEN E)
{
  long i, l = lg(P);
  ulong H = 1;
  for (i = 1; i < l; i++)
  {
    ulong p = P[i], e = E[i];
    long D = (long)(p == 2? d & 7: d % p), a;
    if (s < 0) D = -D;
    a = kross(D,p);
    if (!a)
      H *= upowuu(p, e);
    else
    {
      H *= p - a;
      if (e >= 2) H *= upowuu(p, e-1);
    }
  }
  if (l == 1) return H;
  if (s < 0)
  {
    switch(d)
    { /* divide by [ O_K^* : O^* ] */
      case 4: H >>= 1; break;
      case 3: H /= 3; break;
    }
  }
  else
  {
    GEN fa = mkmat2(zc_to_ZC(P), zc_to_ZC(E));
    H /= itou(quadunitindex_i(utoipos(d), factorback(fa), fa));
  }
  return H;
}
GEN
quadclassnoF_fact(GEN D, GEN P, GEN E)
{
  GEN H = quadclassnoEuler_fact(D, P, E);
  if (lg(P) == 1) return H;
  if (signe(D) < 0)
  {
    switch(itou_or_0(D))
    { /* divide by [ O_K^* : O^* ] */
      case 4: H = shifti(H,-1); break;
      case 3: H = diviuexact(H,3); break;
    }
  }
  else
  {
    GEN fa = mkmat2(P, E);
    H = diviiexact(H, quadunitindex_i(D, factorback(fa), fa));
  }
  return H;
}

static ulong
quadclassnoF_u(ulong x, long s, ulong *pD)
{
  pari_sp av = avma;
  GEN P, E;
  ulong D = coredisc2u_fact(factoru(x), s, &P, &E);
  long H = uquadclassnoF_fact(D, s, P, E);
  *pD = D; return gc_ulong(av, H);
}
ulong
unegquadclassnoF(ulong x, ulong *pD) { return quadclassnoF_u(x, -1, pD); }
ulong
uposquadclassnoF(ulong x, ulong *pD) { return quadclassnoF_u(x, 1, pD); }

/* *pD = coredisc(x), *pR = regulator (x > 0) or NULL */
GEN
quadclassnoF(GEN x, GEN *pD)
{
  GEN D, P, E;
  if (lgefint(x) == 3)
  {
    long s = signe(x);
    ulong d, h = quadclassnoF_u(x[2], s, &d);
    if (pD) *pD = s > 0? utoipos(d): utoineg(d);
    return utoipos(h);
  }
  D = coredisc2_fact(absZ_factor(x), signe(x), &P, &E);
  if (pD) *pD = D;
  return quadclassnoF_fact(D, P, E);
}

static long
two_rank(GEN x)
{
  GEN p = gel(absZ_factor(x),1);
  long l = lg(p)-1;
#if 0 /* positive disc not needed */
  if (signe(x) > 0)
  {
    long i;
    for (i=1; i<=l; i++)
      if (mod4(gel(p,i)) == 3) { l--; break; }
  }
#endif
  return l-1;
}

static GEN
sqr_primeform(GEN x, ulong p) { return qfbsqr_i(primeform_u(x, p)); }
/* return a set of forms hopefully generating Cl(K)^2; set L ~ L(chi_D,1) */
static GEN
get_forms(GEN D, GEN *pL)
{
  const long MAXFORM = 20;
  GEN L, sqrtD = gsqrt(absi_shallow(D),DEFAULTPREC);
  GEN forms = vectrunc_init(MAXFORM+1);
  long s, nforms = 0;
  ulong p;
  forprime_t S;
  L = mulrr(divrr(sqrtD,mppi(DEFAULTPREC)), dbltor(1.005));/*overshoot by 0.5%*/
  s = itos_or_0( truncr(shiftr(sqrtr(sqrtD), 1)) );
  if (!s) pari_err_OVERFLOW("classno [discriminant too large]");
  if      (s < 10)   s = 200;
  else if (s < 20)   s = 1000;
  else if (s < 5000) s = 5000;
  u_forprime_init(&S, 2, s);
  while ( (p = u_forprime_next(&S)) )
  {
    long d, k = kroiu(D,p);
    pari_sp av2;
    if (!k) continue;
    if (k > 0)
    {
      if (++nforms < MAXFORM) vectrunc_append(forms, sqr_primeform(D,p));
      d = p-1;
    }
    else
      d = p+1;
    av2 = avma; affrr(divru(mulur(p,L),d), L); set_avma(av2);
  }
  *pL = L; return forms;
}

/* h ~ #G, return o = order of f, set fao = its factorization */
static  GEN
Shanks_order(void *E, GEN f, GEN h, GEN *pfao)
{
  long s = minss(itos(sqrti(h)), 10000);
  GEN T = gen_Shanks_init(f, s, E, &qfi_group);
  GEN v = gen_Shanks(T, ginv(f), ULONG_MAX, E, &qfi_group);
  return find_order(E, f, addiu(v,1), pfao);
}

/* if g = 1 in  G/<f> ? */
static int
equal1(void *E, GEN T, ulong N, GEN g)
{ return !!gen_Shanks(T, g, N, E, &qfi_group); }

/* Order of 'a' in G/<f>, T = gen_Shanks_init(f,n), order(f) < n*N
 * FIXME: should be gen_order, but equal1 has the wrong prototype */
static GEN
relative_order(void *E, GEN a, GEN o, ulong N,  GEN T)
{
  pari_sp av = avma;
  long i, l;
  GEN m;

  m = get_arith_ZZM(o);
  if (!m) pari_err_TYPE("gen_order [missing order]",a);
  o = gel(m,1);
  m = gel(m,2); l = lgcols(m);
  for (i = l-1; i; i--)
  {
    GEN t, y, p = gcoeff(m,i,1);
    long j, e = itos(gcoeff(m,i,2));
    if (l == 2) {
      t = gen_1;
      y = a;
    } else {
      t = diviiexact(o, powiu(p,e));
      y = powgi(a, t);
    }
    if (equal1(E, T,N,y)) o = t;
    else {
      for (j = 1; j < e; j++)
      {
        y = powgi(y, p);
        if (equal1(E, T,N,y)) break;
      }
      if (j < e) {
        if (j > 1) p = powiu(p, j);
        o = mulii(t, p);
      }
    }
  }
  return gerepilecopy(av, o);
}

/* h(x) for x<0 using Baby Step/Giant Step.
 * Assumes G is not too far from being cyclic.
 *
 * Compute G^2 instead of G so as to kill most of the noncyclicity */
GEN
classno(GEN x)
{
  pari_sp av = avma;
  long r2, k, s, i, l;
  GEN forms, hin, Hf, D, g1, d1, d2, q, L, fad1, order_bound;
  void *E;

  if (signe(x) >= 0) return classno2(x);

  check_quaddisc(x, &s, &k, "classno");
  if (abscmpiu(x,12) <= 0) return gen_1;

  Hf = quadclassnoF(x, &D);
  if (abscmpiu(D,12) <= 0) return gerepilecopy(av, Hf);
  forms =  get_forms(D, &L);
  r2 = two_rank(D);
  hin = roundr(shiftr(L, -r2)); /* rough approximation for #G, G = Cl(K)^2 */

  l = lg(forms);
  order_bound = const_vec(l-1, NULL);
  E = expi(D) > 60? (void*)sqrtnint(shifti(absi_shallow(D),-2),4): NULL;
  g1 = gel(forms,1);
  gel(order_bound,1) = d1 = Shanks_order(E, g1, hin, &fad1);
  q = diviiround(hin, d1); /* approximate order of G/<g1> */
  d2 = NULL; /* not computed yet */
  if (is_pm1(q)) goto END;
  for (i=2; i < l; i++)
  {
    GEN o, fao, a, F, fd, f = gel(forms,i);
    fd = qfbpow_i(f, d1); if (is_pm1(gel(fd,1))) continue;
    F = qfbpow_i(fd, q);
    a = gel(F,1);
    o = is_pm1(a)? find_order(E, fd, q, &fao): Shanks_order(E, fd, q, &fao);
    /* f^(d1 q) = 1 */
    fao = merge_factor(fad1,fao, (void*)&cmpii, &cmp_nodata);
    o = find_order(E, f, fao, &fao);
    gel(order_bound,i) = o;
    /* o = order of f, fao = factor(o) */
    update_g1(&g1,&d1,&fad1, f,o,fao);
    q = diviiround(hin, d1);
    if (is_pm1(q)) goto END;
  }
  /* very probably d1 = expo(Cl^2(K)), q ~ #Cl^2(K) / d1 */
  if (expi(q) > 3)
  { /* q large: compute d2, 2nd elt divisor */
    ulong N, n = 2*itou(sqrti(d1));
    GEN D = d1, T = gen_Shanks_init(g1, n, E, &qfi_group);
    d2 = gen_1;
    N = itou( gceil(gdivgu(d1,n)) ); /* order(g1) <= n*N */
    for (i = 1; i < l; i++)
    {
      GEN d, f = gel(forms,i), B = gel(order_bound,i);
      if (!B) B = find_order(E, f, fad1, /*junk*/&d);
      f = qfbpow_i(f,d2);
      if (equal1(E, T,N,f)) continue;
      B = gdiv(B,d2); if (typ(B) == t_FRAC) B = gel(B,1);
      /* f^B = 1 */
      d = relative_order(E, f, B, N,T);
      d2= mulii(d,d2);
      D = mulii(d1,d2);
      q = diviiround(hin,D);
      if (is_pm1(q)) { d1 = D; goto END; }
    }
    /* very probably, d2 is the 2nd elementary divisor */
    d1 = D; /* product of first two elt divisors */
  }
  /* impose q | d2^oo (d1^oo if d2 not computed), and compatible with known
   * 2-rank */
  if (!ok_q(q,d1,d2,r2))
  {
    GEN q0 = q;
    long d;
    if (cmpii(mulii(q,d1), hin) < 0)
    { /* try q = q0+1,-1,+2,-2 */
      d = 1;
      do { q = addis(q0,d); d = d>0? -d: 1-d; } while(!ok_q(q,d1,d2,r2));
    }
    else
    { /* q0-1,+1,-2,+2  */
      d = -1;
      do { q = addis(q0,d); d = d<0? -d: -1-d; } while(!ok_q(q,d1,d2,r2));
    }
  }
  d1 = mulii(d1,q);

END:
  return gerepileuptoint(av, shifti(mulii(d1,Hf), r2));
}

/* use Euler products */
GEN
classno2(GEN x)
{
  pari_sp av = avma;
  const long prec = DEFAULTPREC;
  long n, i, s;
  GEN p1, p2, S, p4, p5, p7, Hf, Pi, logd, sqrtd, D, half, reg = NULL;

  check_quaddisc(x, &s, NULL, "classno2");
  if (s < 0 && abscmpiu(x,12) <= 0) return gen_1;

  Hf = quadclassnoF(x, &D);
  if (s < 0 && abscmpiu(D,12) <= 0) return gerepilecopy(av, Hf); /* |D| < 12*/

  Pi = mppi(prec);
  sqrtd = sqrtr_abs(itor(D, prec));
  logd = logr_abs(sqrtd); shiftr_inplace(logd,1);
  p1 = sqrtr_abs(divrr(mulir(D,logd), gmul2n(Pi,1)));
  if (s > 0)
  {
    GEN invlogd = invr(logd);
    reg = quadregulator(D, prec);
    p2 = subsr(1, shiftr(mulrr(logr_abs(reg),invlogd),1));
    if (cmprr(sqrr(p2), shiftr(invlogd,1)) >= 0) p1 = mulrr(p2,p1);
  }
  n = itos_or_0( mptrunc(p1) );
  if (!n) pari_err_OVERFLOW("classno [discriminant too large]");

  p4 = divri(Pi, D); setsigne(p4, 1);
  p7 = invr(sqrtr_abs(Pi));
  half = real2n(-1, prec);
  if (s > 0)
  { /* i = 1, shortcut */
    p1 = sqrtd;
    p5 = subsr(1, mulrr(p7,incgamc(half,p4,prec)));
    S = addrr(mulrr(p1,p5), eint1(p4,prec));
    for (i=2; i<=n; i++)
    {
      long k = kroiu(D,i); if (!k) continue;
      p2 = mulir(sqru(i), p4);
      p5 = subsr(1, mulrr(p7,incgamc(half,p2,prec)));
      p5 = addrr(divru(mulrr(p1,p5),i), eint1(p2,prec));
      S = (k>0)? addrr(S,p5): subrr(S,p5);
    }
    S = shiftr(divrr(S,reg),-1);
  }
  else
  { /* i = 1, shortcut */
    p1 = gdiv(sqrtd, Pi);
    p5 = subsr(1, mulrr(p7,incgamc(half,p4,prec)));
    S = addrr(p5, divrr(p1, mpexp(p4)));
    for (i=2; i<=n; i++)
    {
      long k = kroiu(D,i); if (!k) continue;
      p2 = mulir(sqru(i), p4);
      p5 = subsr(1, mulrr(p7,incgamc(half,p2,prec)));
      p5 = addrr(p5, divrr(p1, mulur(i, mpexp(p2))));
      S = (k>0)? addrr(S,p5): subrr(S,p5);
    }
  }
  return gerepileuptoint(av, mulii(Hf, roundr(S)));
}

/* 1 + q + ... + q^v, v > 0 */
static GEN
geomsumu(ulong q, long v)
{
  GEN u = utoipos(1+q);
  for (; v > 1; v--) u = addui(1, mului(q, u));
  return u;
}
static GEN
geomsum(GEN q, long v)
{
  GEN u;
  if (lgefint(q) == 3) return geomsumu(q[2], v);
  u = addiu(q,1);
  for (; v > 1; v--) u = addui(1, mulii(q, u));
  return u;
}

/* 1+p+...+p^(e-1), e >= 1; assuming result fits in an ulong */
static ulong
usumpow(ulong p, long e)
{
  ulong q;
  long i;
  if (p == 2) return (1UL << e) - 1; /* also OK if e = BITS_IN_LONG */
  e--; for (i = 1, q = p + 1; i < e; i++) q = p*q + 1;
  return q;
}
long
uhclassnoF_fact(GEN faF, long D)
{
  GEN P = gel(faF,1), E = gel(faF,2);
  long i, t, l = lg(P);
  for (i = t = 1; i < l; i++)
  {
    long p = P[i], e = E[i], s = kross(D,p);
    if (e == 1) { t *= 1 + p - s; continue; }
    if (s == 1) { t *= upowuu(p,e); continue; }
    t *= 1 + usumpow(p, e) * (p - s);
  }
  return t;
}
/* Hurwitz(D F^2)/ Hurwitz(D)
 * = \sum_{f|F}  f \prod_{p|f} (1-kro(D/p)/p)
 * = \prod_{p^e || F} (1 + (p^e-1) / (p-1) * (p-kro(D/p))) */
GEN
hclassnoF_fact(GEN P, GEN E, GEN D)
{
  GEN H;
  long i, l = lg(P);
  if (l == 1) return gen_1;
  for (i = 1, H = NULL; i < l; i++)
  {
    GEN t, p = gel(P,i);
    long e = E[i], s = kronecker(D,p);
    if (e == 1) t = addiu(p, 1-s);
    else if (s == 1) t = powiu(p, e);
    else t = addui(1, mulii(subis(p, s), geomsum(p, e-1)));
    H = H? mulii(H,t): t;
  }
  return H;
}
static GEN
hclassno6_large(GEN x)
{
  GEN H = NULL, P, E, D = coredisc2_fact(absZ_factor(x), -1, &P, &E);
  long l = lg(P);

  if (l > 1 && lgefint(x) == 3)
  { /* F != 1, second chance */
    ulong h = hclassno6u_no_cache(x[2]);
    if (h) H = utoipos(h);
  }
  if (!H)
  {
    H = quadclassno(D);
    switch(itou_or_0(D))
    {
      case 3: H = shifti(H,1);break;
      case 4: H = muliu(H,3); break;
      default:H = muliu(H,6); break;
    }
  }
  return mulii(H, hclassnoF_fact(P, E, D));
}

/* x > 0, x = 0,3 (mod 4). Return 6*hclassno(x), an integer */
GEN
hclassno6(GEN x)
{
  ulong d = itou_or_0(x);
  if (d)
  { /* create cache if d small */
    ulong h = d < 500000 ? hclassno6u(d): hclassno6u_no_cache(d);
    if (h) return utoipos(h);
  }
  return hclassno6_large(x);
}

GEN
hclassno(GEN x)
{
  long a, s;
  if (typ(x) != t_INT) pari_err_TYPE("hclassno",x);
  s = signe(x);
  if (s < 0) return gen_0;
  if (!s) return gdivgs(gen_1, -12);
  a = mod4(x); if (a == 1 || a == 2) return gen_0;
  return gdivgu(hclassno6(x), 6);
}

/* return [N',v]; v contains all x mod N' s.t. x^2 + B x + C = 0 modulo N */
GEN
Zn_quad_roots(GEN N, GEN B, GEN C)
{
  pari_sp av = avma;
  GEN fa = NULL, D, w, v, P, E, F0, Q0, F, mF, A, Q, T, R, Np, N4;
  long l, i, j, ct;

  if ((fa = check_arith_non0(N,"Zn_quad_roots")))
  {
    N = typ(N) == t_VEC? gel(N,1): factorback(N);
    fa = clean_Z_factor(fa);
  }
  N = absi_shallow(N);
  N4 = shifti(N,2);
  D = modii(subii(sqri(B), shifti(C,2)), N4);
  if (!signe(D))
  { /* (x + B/2)^2 = 0 (mod N), D = B^2-4C = 0 (4N) => B even */
    if (!fa) fa = Z_factor(N);
    P = gel(fa,1);
    E = ZV_to_zv(gel(fa,2));
    l = lg(P);
    for (i = 1; i < l; i++) E[i] = (E[i]+1) >> 1;
    Np = factorback2(P, E); /* x = -B mod N' */
    B = shifti(B,-1);
    return gerepilecopy(av, mkvec2(Np, mkvec(Fp_neg(B,Np))));
  }
  if (!fa)
    fa = Z_factor(N4);
  else  /* convert to factorization of N4 = 4*N */
    fa = famat_reduce(famat_mulpows_shallow(fa, gen_2, 2));
  P = gel(fa,1); l = lg(P);
  E = ZV_to_zv(gel(fa,2));
  F = cgetg(l, t_VEC);
  mF= cgetg(l, t_VEC); F0 = gen_0;
  Q = cgetg(l, t_VEC); Q0 = gen_1;
  for (i = j = 1, ct = 0; i < l; i++)
  {
    GEN p = gel(P,i), q, f, mf, D0;
    long t2, s = E[i], t = Z_pvalrem(D, p, &D0), d = s - t;
    if (d <= 0)
    {
      q = powiu(p, (s+1)>>1);
      Q0 = mulii(Q0, q); continue;
    }
    /* d > 0 */
    if (odd(t)) return NULL;
    t2 = t >> 1;
    if (i > 1)
    { /* p > 2 */
      if (kronecker(D0, p) == -1) return NULL;
      q = powiu(p, s - t2);
      f = Zp_sqrt(D0, p, d);
      if (!f) return NULL; /* p was not actually prime... */
      if (t2) f = mulii(powiu(p,t2), f);
      mf = Fp_neg(f, q);
    }
    else
    { /* p = 2 */
      if (d <= 3)
      {
        if (d == 3 && Mod8(D0) != 1) return NULL;
        if (d == 2 && Mod4(D0) != 1) return NULL;
        Q0 = int2n(1+t2); F0 = NULL; continue;
      }
      if (Mod8(D0) != 1) return NULL;
      q = int2n(d - 1 + t2);
      f = Z2_sqrt(D0, d);
      if (t2) f = shifti(f, t2);
      mf = Fp_neg(f, q);
    }
    gel(Q,j) = q;
    gel(F,j) = f;
    gel(mF,j)= mf; j++;
  }
  setlg(Q,j);
  setlg(F,j);
  setlg(mF,j);
  if (is_pm1(Q0)) A = leafcopy(F);
  else
  { /* append the fixed congruence (F0 mod Q0) */
    if (!F0) F0 = shifti(Q0,-1);
    A = shallowconcat(F, F0);
    Q = shallowconcat(Q, Q0);
  }
  ct = 1 << (j-1);
  T = ZV_producttree(Q);
  R = ZV_chinesetree(Q,T);
  Np = gmael(T, lg(T)-1, 1);
  B = modii(B, Np);
  if (!signe(B)) B = NULL;
  Np = shifti(Np, -1); /* N' = (\prod_i Q[i]) / 2 */
  w = cgetg(3, t_VEC);
  gel(w,1) = icopy(Np);
  gel(w,2) = v = cgetg(ct+1, t_VEC);
  l = lg(F);
  for (j = 1; j <= ct; j++)
  {
    pari_sp av2 = avma;
    long m = j - 1;
    GEN u;
    for (i = 1; i < l; i++)
    {
      gel(A,i) = (m&1L)? gel(mF,i): gel(F,i);
      m >>= 1;
    }
    u = ZV_chinese_tree(A,Q,T,R); /* u mod N' st u^2 = B^2-4C modulo 4N */
    if (B) u = subii(u,B);
    gel(v,j) = gerepileuptoint(av2, modii(shifti(u,-1), Np));
  }
  return gerepileupto(av, w);
}
