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

/***********************************************************************/
/**                                                                   **/
/**               ARITHMETIC OPERATIONS ON POLYNOMIALS                **/
/**                         (third part)                              **/
/**                                                                   **/
/***********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_pol

/************************************************************************
 **                                                                    **
 **                      Ring membership                               **
 **                                                                    **
 ************************************************************************/
struct charact {
  GEN q;
  int isprime;
};
static void
char_update_prime(struct charact *S, GEN p)
{
  if (!S->isprime) { S->isprime = 1; S->q = p; }
  if (!equalii(p, S->q)) pari_err_MODULUS("characteristic", S->q, p);
}
static void
char_update_int(struct charact *S, GEN n)
{
  if (S->isprime)
  {
    if (dvdii(n, S->q)) return;
    pari_err_MODULUS("characteristic", S->q, n);
  }
  S->q = gcdii(S->q, n);
}
static void
charact(struct charact *S, GEN x)
{
  const long tx = typ(x);
  long i, l;
  switch(tx)
  {
    case t_INTMOD:char_update_int(S, gel(x,1)); break;
    case t_FFELT: char_update_prime(S, gel(x,4)); break;
    case t_COMPLEX: case t_QUAD:
    case t_POLMOD: case t_POL: case t_SER: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      l = lg(x);
      for (i=lontyp[tx]; i < l; i++) charact(S,gel(x,i));
      break;
    case t_LIST:
      x = list_data(x);
      if (x) charact(S, x);
      break;
  }
}
static void
charact_res(struct charact *S, GEN x)
{
  const long tx = typ(x);
  long i, l;
  switch(tx)
  {
    case t_INTMOD:char_update_int(S, gel(x,1)); break;
    case t_FFELT: char_update_prime(S, gel(x,4)); break;
    case t_PADIC: char_update_prime(S, gel(x,2)); break;
    case t_COMPLEX: case t_QUAD:
    case t_POLMOD: case t_POL: case t_SER: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      l = lg(x);
      for (i=lontyp[tx]; i < l; i++) charact_res(S,gel(x,i));
      break;
    case t_LIST:
      x = list_data(x);
      if (x) charact_res(S, x);
      break;
  }
}
GEN
characteristic(GEN x)
{
  struct charact S;
  S.q = gen_0; S.isprime = 0;
  charact(&S, x); return S.q;
}
GEN
residual_characteristic(GEN x)
{
  struct charact S;
  S.q = gen_0; S.isprime = 0;
  charact_res(&S, x); return S.q;
}

int
Rg_is_Fp(GEN x, GEN *pp)
{
  GEN mod;
  switch(typ(x))
  {
  case t_INTMOD:
    mod = gel(x,1);
    if (!*pp) *pp = mod;
    else if (mod != *pp && !equalii(mod, *pp))
    {
      if (DEBUGLEVEL) pari_warn(warner,"different moduli in Rg_is_Fp");
      return 0;
    }
    return 1;
  case t_INT:
    return 1;
  default: return 0;
  }
}

int
RgX_is_FpX(GEN x, GEN *pp)
{
  long i, lx = lg(x);
  for (i=2; i<lx; i++)
    if (!Rg_is_Fp(gel(x, i), pp))
      return 0;
  return 1;
}

int
RgV_is_FpV(GEN x, GEN *pp)
{
  long i, lx = lg(x);
  for (i=1; i<lx; i++)
    if (!Rg_is_Fp(gel(x,i), pp)) return 0;
  return 1;
}

int
RgM_is_FpM(GEN x, GEN *pp)
{
  long i, lx = lg(x);
  for (i=1; i<lx; i++)
    if (!RgV_is_FpV(gel(x, i), pp)) return 0;
  return 1;
}

int
Rg_is_FpXQ(GEN x, GEN *pT, GEN *pp)
{
  GEN pol, mod, p;
  switch(typ(x))
  {
  case t_INTMOD:
    return Rg_is_Fp(x, pp);
  case t_INT:
    return 1;
  case t_POL:
    return RgX_is_FpX(x, pp);
  case t_FFELT:
    mod = x; p = FF_p_i(x);
    if (!*pp) *pp = p;
    if (!*pT) *pT = mod;
    else if (typ(*pT)!=t_FFELT || !FF_samefield(*pT,mod))
    {
      if (DEBUGLEVEL) pari_warn(warner,"different moduli in Rg_is_FpXQ");
      return 0;
    }
    return 1;
  case t_POLMOD:
    mod = gel(x,1); pol = gel(x, 2);
    if (!RgX_is_FpX(mod, pp)) return 0;
    if (typ(pol)==t_POL)
    {
      if (!RgX_is_FpX(pol, pp)) return 0;
    }
    else if (!Rg_is_Fp(pol, pp)) return 0;
    if (!*pT) *pT = mod;
    else if (mod != *pT && !gequal(mod, *pT))
    {
      if (DEBUGLEVEL) pari_warn(warner,"different moduli in Rg_is_FpXQ");
      return 0;
    }
    return 1;

  default: return 0;
  }
}

int
RgX_is_FpXQX(GEN x, GEN *pT, GEN *pp)
{
  long i, lx = lg(x);
  for (i = 2; i < lx; i++)
    if (!Rg_is_FpXQ(gel(x,i), pT, pp)) return 0;
  return 1;
}

/************************************************************************
 **                                                                    **
 **                      Ring conversion                               **
 **                                                                    **
 ************************************************************************/

/* p > 0 a t_INT, return lift(x * Mod(1,p)).
 * If x is an INTMOD, assume modulus is a multiple of p. */
GEN
Rg_to_Fp(GEN x, GEN p)
{
  if (lgefint(p) == 3) return utoi(Rg_to_Fl(x, uel(p,2)));
  switch(typ(x))
  {
    case t_INT: return modii(x, p);
    case t_FRAC: {
      pari_sp av = avma;
      GEN z = modii(gel(x,1), p);
      if (z == gen_0) return gen_0;
      return gerepileuptoint(av, remii(mulii(z, Fp_inv(gel(x,2), p)), p));
    }
    case t_PADIC: return padic_to_Fp(x, p);
    case t_INTMOD: {
      GEN q = gel(x,1), a = gel(x,2);
      if (equalii(q, p)) return icopy(a);
      if (!dvdii(q,p)) pari_err_MODULUS("Rg_to_Fp", q, p);
      return remii(a, p);
    }
    default: pari_err_TYPE("Rg_to_Fp",x);
      return NULL; /* LCOV_EXCL_LINE */
  }
}
/* If x is a POLMOD, assume modulus is a multiple of T. */
GEN
Rg_to_FpXQ(GEN x, GEN T, GEN p)
{
  long ta, tx = typ(x), v = get_FpX_var(T);
  GEN a, b;
  if (is_const_t(tx))
  {
    if (tx == t_FFELT)
    {
      GEN z = FF_to_FpXQ(x);
      setvarn(z, v);
      return z;
    }
    return scalar_ZX(degpol(T)? Rg_to_Fp(x, p): gen_0, v);
  }
  switch(tx)
  {
    case t_POLMOD:
      b = gel(x,1);
      a = gel(x,2); ta = typ(a);
      if (is_const_t(ta))
        return scalar_ZX(degpol(T)? Rg_to_Fp(a, p): gen_0, v);
      b = RgX_to_FpX(b, p); if (varn(b) != v) break;
      a = RgX_to_FpX(a, p);
      if (ZX_equal(b,get_FpX_mod(T)) || signe(FpX_rem(b,T,p))==0)
        return FpX_rem(a, T, p);
      break;
    case t_POL:
      if (varn(x) != v) break;
      return FpX_rem(RgX_to_FpX(x,p), T, p);
    case t_RFRAC:
      a = Rg_to_FpXQ(gel(x,1), T,p);
      b = Rg_to_FpXQ(gel(x,2), T,p);
      return FpXQ_div(a,b, T,p);
  }
  pari_err_TYPE("Rg_to_FpXQ",x);
  return NULL; /* LCOV_EXCL_LINE */
}
GEN
RgX_to_FpX(GEN x, GEN p)
{
  long i, l;
  GEN z = cgetg_copy(x, &l); z[1] = x[1];
  for (i = 2; i < l; i++) gel(z,i) = Rg_to_Fp(gel(x,i), p);
  return FpX_renormalize(z, l);
}

GEN
RgV_to_FpV(GEN x, GEN p)
{ pari_APPLY_type(t_VEC, Rg_to_Fp(gel(x,i), p)) }

GEN
RgC_to_FpC(GEN x, GEN p)
{ pari_APPLY_type(t_COL, Rg_to_Fp(gel(x,i), p)) }

GEN
RgM_to_FpM(GEN x, GEN p)
{ pari_APPLY_same(RgC_to_FpC(gel(x,i), p)) }

GEN
RgV_to_Flv(GEN x, ulong p)
{ pari_APPLY_ulong(Rg_to_Fl(gel(x,i), p)) }

GEN
RgM_to_Flm(GEN x, ulong p)
{ pari_APPLY_same(RgV_to_Flv(gel(x,i), p)) }

GEN
RgX_to_FpXQX(GEN x, GEN T, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_POL); z[1] = x[1];
  for (i = 2; i < l; i++) gel(z,i) = Rg_to_FpXQ(gel(x,i), T,p);
  return FpXQX_renormalize(z, l);
}
GEN
RgX_to_FqX(GEN x, GEN T, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_POL); z[1] = x[1];
  if (T)
    for (i = 2; i < l; i++) gel(z,i) = Rg_to_FpXQ(gel(x,i), T, p);
  else
    for (i = 2; i < l; i++) gel(z,i) = Rg_to_Fp(gel(x,i), p);
  return FpXQX_renormalize(z, l);
}

GEN
RgC_to_FqC(GEN x, GEN T, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_COL);
  if (T)
    for (i = 1; i < l; i++) gel(z,i) = Rg_to_FpXQ(gel(x,i), T, p);
  else
    for (i = 1; i < l; i++) gel(z,i) = Rg_to_Fp(gel(x,i), p);
  return z;
}

GEN
RgM_to_FqM(GEN x, GEN T, GEN p)
{ pari_APPLY_same(RgC_to_FqC(gel(x, i), T, p)) }

/* lg(V) > 1 */
GEN
FpXV_FpC_mul(GEN V, GEN W, GEN p)
{
  pari_sp av = avma;
  long i, l = lg(V);
  GEN z = ZX_Z_mul(gel(V,1),gel(W,1));
  for(i=2; i<l; i++)
  {
    z = ZX_add(z, ZX_Z_mul(gel(V,i),gel(W,i)));
    if ((i & 7) == 0) z = gerepileupto(av, z);
  }
  return gerepileupto(av, FpX_red(z,p));
}

GEN
FqX_Fq_add(GEN y, GEN x, GEN T, GEN p)
{
  long i, lz = lg(y);
  GEN z;
  if (!T) return FpX_Fp_add(y, x, p);
  if (lz == 2) return scalarpol(x, varn(y));
  z = cgetg(lz,t_POL); z[1] = y[1];
  gel(z,2) = Fq_add(gel(y,2),x, T, p);
  if (lz == 3) z = FpXX_renormalize(z,lz);
  else
    for(i=3;i<lz;i++) gel(z,i) = gcopy(gel(y,i));
  return z;
}

GEN
FqX_Fq_sub(GEN y, GEN x, GEN T, GEN p)
{
  long i, lz = lg(y);
  GEN z;
  if (!T) return FpX_Fp_sub(y, x, p);
  if (lz == 2) return scalarpol(x, varn(y));
  z = cgetg(lz,t_POL); z[1] = y[1];
  gel(z,2) = Fq_sub(gel(y,2), x, T, p);
  if (lz == 3) z = FpXX_renormalize(z,lz);
  else
    for(i=3;i<lz;i++) gel(z,i) = gcopy(gel(y,i));
  return z;
}

GEN
FqX_Fq_mul_to_monic(GEN P, GEN U, GEN T, GEN p)
{
  long i, lP;
  GEN res = cgetg_copy(P, &lP); res[1] = P[1];
  for(i=2; i<lP-1; i++) gel(res,i) = Fq_mul(U,gel(P,i), T,p);
  gel(res,lP-1) = gen_1; return res;
}

GEN
FpXQX_normalize(GEN z, GEN T, GEN p)
{
  GEN lc;
  if (lg(z) == 2) return z;
  lc = leading_coeff(z);
  if (typ(lc) == t_POL)
  {
    if (lg(lc) > 3) /* nonconstant */
      return FqX_Fq_mul_to_monic(z, Fq_inv(lc,T,p), T,p);
    /* constant */
    lc = gel(lc,2);
    z = shallowcopy(z);
    gel(z, lg(z)-1) = lc;
  }
  /* lc a t_INT */
  if (equali1(lc)) return z;
  return FqX_Fq_mul_to_monic(z, Fp_inv(lc,p), T,p);
}

GEN
FqX_eval(GEN x, GEN y, GEN T, GEN p)
{
  pari_sp av;
  GEN p1, r;
  long j, i=lg(x)-1;
  if (i<=2)
    return (i==2)? Fq_red(gel(x,2), T, p): gen_0;
  av=avma; p1=gel(x,i);
  /* specific attention to sparse polynomials (see poleval)*/
  /*You've guessed it! It's a copy-paste(tm)*/
  for (i--; i>=2; i=j-1)
  {
    for (j=i; !signe(gel(x,j)); j--)
      if (j==2)
      {
        if (i!=j) y = Fq_pow(y,utoipos(i-j+1), T, p);
        return gerepileupto(av, Fq_mul(p1,y, T, p));
      }
    r = (i==j)? y: Fq_pow(y, utoipos(i-j+1), T, p);
    p1 = Fq_add(Fq_mul(p1,r,T,p), gel(x,j), T, p);
  }
  return gerepileupto(av, p1);
}

GEN
FqXY_evalx(GEN Q, GEN x, GEN T, GEN p)
{
  long i, lb = lg(Q);
  GEN z;
  if (!T) return FpXY_evalx(Q, x, p);
  z = cgetg(lb, t_POL); z[1] = Q[1];
  for (i=2; i<lb; i++)
  {
    GEN q = gel(Q,i);
    gel(z,i) = typ(q) == t_INT? modii(q,p): FqX_eval(q, x, T, p);
  }
  return FpXQX_renormalize(z, lb);
}

/* Q an FpXY, evaluate at (X,Y) = (x,y) */
GEN
FqXY_eval(GEN Q, GEN y, GEN x, GEN T, GEN p)
{
  pari_sp av = avma;
  if (!T) return FpXY_eval(Q, y, x, p);
  return gerepileupto(av, FqX_eval(FqXY_evalx(Q, x, T, p), y, T, p));
}

/* a X^d */
GEN
monomial(GEN a, long d, long v)
{
  long i, n;
  GEN P;
  if (d < 0) {
    if (isrationalzero(a)) return pol_0(v);
    retmkrfrac(a, pol_xn(-d, v));
  }
  if (gequal0(a))
  {
    if (isexactzero(a)) return scalarpol_shallow(a,v);
    n = d+2; P = cgetg(n+1, t_POL);
    P[1] = evalsigne(0) | evalvarn(v);
  }
  else
  {
    n = d+2; P = cgetg(n+1, t_POL);
    P[1] = evalsigne(1) | evalvarn(v);
  }
  for (i = 2; i < n; i++) gel(P,i) = gen_0;
  gel(P,i) = a; return P;
}
GEN
monomialcopy(GEN a, long d, long v)
{
  long i, n;
  GEN P;
  if (d < 0) {
    if (isrationalzero(a)) return pol_0(v);
    retmkrfrac(gcopy(a), pol_xn(-d, v));
  }
  if (gequal0(a))
  {
    if (isexactzero(a)) return scalarpol(a,v);
    n = d+2; P = cgetg(n+1, t_POL);
    P[1] = evalsigne(0) | evalvarn(v);
  }
  else
  {
    n = d+2; P = cgetg(n+1, t_POL);
    P[1] = evalsigne(1) | evalvarn(v);
  }
  for (i = 2; i < n; i++) gel(P,i) = gen_0;
  gel(P,i) = gcopy(a); return P;
}
GEN
pol_x_powers(long N, long v)
{
  GEN L = cgetg(N+1,t_VEC);
  long i;
  for (i=1; i<=N; i++) gel(L,i) = pol_xn(i-1, v);
  return L;
}

GEN
FqXQ_powers(GEN x, long l, GEN S, GEN T, GEN p)
{
  return T ? FpXQXQ_powers(x, l, S, T, p): FpXQ_powers(x, l, S, p);
}

GEN
FqXQ_matrix_pow(GEN y, long n, long m, GEN S, GEN T, GEN p)
{
  return T ? FpXQXQ_matrix_pow(y, n, m, S, T, p): FpXQ_matrix_pow(y, n, m, S, p);
}

/*******************************************************************/
/*                                                                 */
/*                             Fq                                  */
/*                                                                 */
/*******************************************************************/

GEN
Fq_add(GEN x, GEN y, GEN T/*unused*/, GEN p)
{
  (void)T;
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return Fp_add(x,y,p);
    case 1: return FpX_Fp_add(x,y,p);
    case 2: return FpX_Fp_add(y,x,p);
    case 3: return FpX_add(x,y,p);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
Fq_sub(GEN x, GEN y, GEN T/*unused*/, GEN p)
{
  (void)T;
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return Fp_sub(x,y,p);
    case 1: return FpX_Fp_sub(x,y,p);
    case 2: return Fp_FpX_sub(x,y,p);
    case 3: return FpX_sub(x,y,p);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
Fq_neg(GEN x, GEN T/*unused*/, GEN p)
{
  (void)T;
  return (typ(x)==t_POL)? FpX_neg(x,p): Fp_neg(x,p);
}

GEN
Fq_halve(GEN x, GEN T/*unused*/, GEN p)
{
  (void)T;
  return (typ(x)==t_POL)? FpX_halve(x,p): Fp_halve(x,p);
}

/* If T==NULL do not reduce*/
GEN
Fq_mul(GEN x, GEN y, GEN T, GEN p)
{
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return Fp_mul(x,y,p);
    case 1: return FpX_Fp_mul(x,y,p);
    case 2: return FpX_Fp_mul(y,x,p);
    case 3: if (T) return FpXQ_mul(x,y,T,p);
            else return FpX_mul(x,y,p);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

/* If T==NULL do not reduce*/
GEN
Fq_mulu(GEN x, ulong y, /*unused*/GEN T, GEN p)
{
  (void) T;
  return typ(x)==t_POL ? FpX_Fp_mul(x,utoi(y),p): Fp_mulu(x, y, p);
}

/* y t_INT */
GEN
Fq_Fp_mul(GEN x, GEN y, GEN T/*unused*/, GEN p)
{
  (void)T;
  return (typ(x) == t_POL)? FpX_Fp_mul(x,y,p)
                          : Fp_mul(x,y,p);
}
/* If T==NULL do not reduce*/
GEN
Fq_sqr(GEN x, GEN T, GEN p)
{
  if (typ(x) == t_POL)
  {
    if (T) return FpXQ_sqr(x,T,p);
    else return FpX_sqr(x,p);
  }
  else
    return Fp_sqr(x,p);
}

GEN
Fq_neg_inv(GEN x, GEN T, GEN p)
{
  if (typ(x) == t_INT) return Fp_inv(Fp_neg(x,p),p);
  return FpXQ_inv(FpX_neg(x,p),T,p);
}

GEN
Fq_invsafe(GEN x, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_invsafe(x,p);
  return FpXQ_invsafe(x,pol,p);
}

GEN
Fq_inv(GEN x, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_inv(x,p);
  return FpXQ_inv(x,pol,p);
}

GEN
Fq_div(GEN x, GEN y, GEN pol, GEN p)
{
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return Fp_div(x,y,p);
    case 1: return FpX_Fp_div(x,y,p);
    case 2: return FpX_Fp_mul(FpXQ_inv(y,pol,p),x,p);
    case 3: return FpXQ_div(x,y,pol,p);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
Fq_pow(GEN x, GEN n, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_pow(x,n,p);
  return FpXQ_pow(x,n,pol,p);
}

GEN
Fq_powu(GEN x, ulong n, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_powu(x,n,p);
  return FpXQ_powu(x,n,pol,p);
}

GEN
Fq_sqrt(GEN x, GEN T, GEN p)
{
  if (typ(x) == t_INT)
  {
    if (!T || odd(get_FpX_degree(T))) return Fp_sqrt(x,p);
    x = scalarpol_shallow(x, get_FpX_var(T));
  }
  return FpXQ_sqrt(x,T,p);
}
GEN
Fq_sqrtn(GEN x, GEN n, GEN T, GEN p, GEN *zeta)
{
  if (typ(x) == t_INT)
  {
    long d;
    if (!T) return Fp_sqrtn(x,n,p,zeta);
    d = get_FpX_degree(T);
    if (ugcdiu(n,d) == 1)
    {
      if (!zeta) return Fp_sqrtn(x,n,p,NULL);
      /* gcd(n,p-1)=gcd(n,q-1): same number of solutions in Fp and F_q */
      if (equalii(gcdii(subiu(p,1),n), gcdii(subiu(Fp_powu(p,d,n), 1), n)))
        return Fp_sqrtn(x,n,p,zeta);
    }
    x = scalarpol(x, get_FpX_var(T)); /* left on stack */
  }
  return FpXQ_sqrtn(x,n,T,p,zeta);
}

struct _Fq_field
{
  GEN T, p;
};

static GEN
_Fq_red(void *E, GEN x)
{ struct _Fq_field *s = (struct _Fq_field *)E;
  return Fq_red(x, s->T, s->p);
}

static GEN
_Fq_add(void *E, GEN x, GEN y)
{
  (void) E;
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return addii(x,y);
    case 1: return ZX_Z_add(x,y);
    case 2: return ZX_Z_add(y,x);
    default: return ZX_add(x,y);
  }
}

static GEN
_Fq_neg(void *E, GEN x) { (void) E; return typ(x)==t_POL?ZX_neg(x):negi(x); }

static GEN
_Fq_mul(void *E, GEN x, GEN y)
{
  (void) E;
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return mulii(x,y);
    case 1: return ZX_Z_mul(x,y);
    case 2: return ZX_Z_mul(y,x);
    default: return ZX_mul(x,y);
  }
}

static GEN
_Fq_inv(void *E, GEN x)
{ struct _Fq_field *s = (struct _Fq_field *)E;
  return Fq_inv(x,s->T,s->p);
}

static int
_Fq_equal0(GEN x) { return signe(x)==0; }

static GEN
_Fq_s(void *E, long x) { (void) E; return stoi(x); }

static const struct bb_field Fq_field={_Fq_red,_Fq_add,_Fq_mul,_Fq_neg,
                                       _Fq_inv,_Fq_equal0,_Fq_s};

const struct bb_field *get_Fq_field(void **E, GEN T, GEN p)
{
  if (!T)
    return get_Fp_field(E, p);
  else
  {
    GEN z = new_chunk(sizeof(struct _Fq_field));
    struct _Fq_field *e = (struct _Fq_field *) z;
    e->T = T; e->p  = p; *E = (void*)e;
    return &Fq_field;
  }
}

/*******************************************************************/
/*                                                                 */
/*                             Fq[X]                               */
/*                                                                 */
/*******************************************************************/
/* P(X + c) */
GEN
FpX_translate(GEN P, GEN c, GEN p)
{
  pari_sp av = avma;
  GEN Q, *R;
  long i, k, n;

  if (!signe(P) || !signe(c)) return ZX_copy(P);
  Q = leafcopy(P);
  R = (GEN*)(Q+2); n = degpol(P);
  for (i=1; i<=n; i++)
  {
    for (k=n-i; k<n; k++)
      R[k] = Fp_add(R[k], Fp_mul(c, R[k+1], p), p);

    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"FpX_translate, i = %ld/%ld", i,n);
      Q = gerepilecopy(av, Q); R = (GEN*)Q+2;
    }
  }
  return gerepilecopy(av, FpX_renormalize(Q, lg(Q)));
}
/* P(X + c), c an Fq */
GEN
FqX_translate(GEN P, GEN c, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN Q, *R;
  long i, k, n;

  /* signe works for t_(INT|POL) */
  if (!signe(P) || !signe(c)) return RgX_copy(P);
  Q = leafcopy(P);
  R = (GEN*)(Q+2); n = degpol(P);
  for (i=1; i<=n; i++)
  {
    for (k=n-i; k<n; k++)
      R[k] = Fq_add(R[k], Fq_mul(c, R[k+1], T, p), T, p);

    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"FqX_translate, i = %ld/%ld", i,n);
      Q = gerepilecopy(av, Q); R = (GEN*)Q+2;
    }
  }
  return gerepilecopy(av, FpXQX_renormalize(Q, lg(Q)));
}

GEN
FqV_roots_to_pol(GEN V, GEN T, GEN p, long v)
{
  pari_sp ltop = avma;
  long k;
  GEN W;
  if (lgefint(p) == 3)
  {
    ulong pp = p[2];
    GEN Tl = ZX_to_Flx(T, pp);
    GEN Vl = ZXC_to_FlxC(V, pp, get_Flx_var(Tl));
    Tl = FlxqV_roots_to_pol(Vl, Tl, pp, v);
    return gerepileupto(ltop, FlxX_to_ZXX(Tl));
  }
  W = cgetg(lg(V),t_VEC);
  for(k=1; k < lg(V); k++)
    gel(W,k) = deg1pol_shallow(gen_1,Fq_neg(gel(V,k),T,p),v);
  return gerepileupto(ltop, FpXQXV_prod(W, T, p));
}

GEN
FqV_red(GEN x, GEN T, GEN p)
{ pari_APPLY_same(Fq_red(gel(x,i), T, p)) }

GEN
FqC_add(GEN x, GEN y, GEN T, GEN p)
{
  if (!T) return FpC_add(x, y, p);
  pari_APPLY_type(t_COL, Fq_add(gel(x,i), gel(y,i), T, p))
}

GEN
FqC_sub(GEN x, GEN y, GEN T, GEN p)
{
  if (!T) return FpC_sub(x, y, p);
  pari_APPLY_type(t_COL, Fq_sub(gel(x,i), gel(y,i), T, p))
}

GEN
FqC_Fq_mul(GEN x, GEN y, GEN T, GEN p)
{
  if (!T) return FpC_Fp_mul(x, y, p);
  pari_APPLY_type(t_COL, Fq_mul(gel(x,i),y,T,p))
}

GEN
FqC_FqV_mul(GEN x, GEN y, GEN T, GEN p)
{
  long i,j, lx=lg(x), ly=lg(y);
  GEN z;
  if (ly==1) return cgetg(1,t_MAT);
  z = cgetg(ly,t_MAT);
  for (j=1; j < ly; j++)
  {
    GEN zj = cgetg(lx,t_COL);
    for (i=1; i<lx; i++) gel(zj,i) = Fq_mul(gel(x,i),gel(y,j), T, p);
    gel(z, j) = zj;
  }
  return z;
}

GEN
FpXC_center(GEN x, GEN p, GEN pov2)
{ pari_APPLY_type(t_COL, FpX_center(gel(x,i), p, pov2)) }

GEN
FpXM_center(GEN x, GEN p, GEN pov2)
{ pari_APPLY_same(FpXC_center(gel(x,i), p, pov2)) }

/*******************************************************************/
/*                                                                 */
/*                          GENERIC CRT                            */
/*                                                                 */
/*******************************************************************/
static GEN
primelist(forprime_t *S, long n, GEN dB)
{
  GEN P = cgetg(n+1, t_VECSMALL);
  long i = 1;
  ulong p;
  while (i <= n && (p = u_forprime_next(S)))
    if (!dB || umodiu(dB, p)) P[i++] = p;
  return P;
}

void
gen_inccrt_i(const char *str, GEN worker, GEN dB, long n, long mmin,
             forprime_t *S, GEN *pH, GEN *pmod, GEN crt(GEN, GEN, GEN*),
             GEN center(GEN, GEN, GEN))
{
  long m = mmin? minss(mmin, n): usqrt(n);
  GEN  H, P, mod;
  pari_timer ti;
  if (DEBUGLEVEL > 4)
  {
    timer_start(&ti);
    err_printf("%s: nb primes: %ld\n",str, n);
  }
  if (m == 1)
  {
    GEN P = primelist(S, n, dB);
    GEN done = closure_callgen1(worker, P);
    H = gel(done,1);
    mod = gel(done,2);
    if (!*pH && center) H = center(H, mod, shifti(mod,-1));
    if (DEBUGLEVEL>4) timer_printf(&ti,"%s: modular", str);
  }
  else
  {
    long i, s = (n+m-1)/m, r = m - (m*s-n), di = 0;
    struct pari_mt pt;
    long pending = 0;
    H = cgetg(m+1, t_VEC); P = cgetg(m+1, t_VEC);
    mt_queue_start_lim(&pt, worker, m);
    for (i=1; i<=m || pending; i++)
    {
      GEN done;
      GEN pr = i <= m ? mkvec(primelist(S, i<=r ? s: s-1, dB)): NULL;
      mt_queue_submit(&pt, i, pr);
      done = mt_queue_get(&pt, NULL, &pending);
      if (done)
      {
        di++;
        gel(H, di) = gel(done,1);
        gel(P, di) = gel(done,2);
        if (DEBUGLEVEL>5) err_printf("%ld%% ",100*di/m);
      }
    }
    mt_queue_end(&pt);
    if (DEBUGLEVEL>5) err_printf("\n");
    if (DEBUGLEVEL>4) timer_printf(&ti,"%s: modular", str);
    H = crt(H, P, &mod);
    if (DEBUGLEVEL>4) timer_printf(&ti,"%s: chinese", str);
  }
  if (*pH) H = crt(mkvec2(*pH, H), mkvec2(*pmod, mod), &mod);
  *pH = H; *pmod = mod;
}
void
gen_inccrt(const char *str, GEN worker, GEN dB, long n, long mmin,
           forprime_t *S, GEN *pH, GEN *pmod, GEN crt(GEN, GEN, GEN*),
           GEN center(GEN, GEN, GEN))
{
  pari_sp av = avma;
  gen_inccrt_i(str, worker, dB, n, mmin, S, pH, pmod, crt, center);
  gerepileall(av, 2, pH, pmod);
}

GEN
gen_crt(const char *str, GEN worker, forprime_t *S, GEN dB, ulong bound, long mmin, GEN *pmod,
        GEN crt(GEN, GEN, GEN*), GEN center(GEN, GEN, GEN))
{
  GEN mod = gen_1, H = NULL;
  ulong e;

  bound++;
  while (bound > (e = expi(mod)))
  {
    long n = (bound - e) / expu(S->p) + 1;
    gen_inccrt(str, worker, dB, n, mmin, S, &H, &mod, crt, center);
  }
  if (pmod) *pmod = mod;
  return H;
}

/*******************************************************************/
/*                                                                 */
/*                          MODULAR GCD                            */
/*                                                                 */
/*******************************************************************/
/* return z = a mod q, b mod p (p,q) = 1; qinv = 1/q mod p; a in ]-q,q] */
static GEN
Fl_chinese_coprime(GEN a, ulong b, GEN q, ulong p, ulong qinv, GEN pq, GEN pq2)
{
  ulong d, amod = umodiu(a, p);
  pari_sp av = avma;
  GEN ax;

  if (b == amod) return NULL;
  d = Fl_mul(Fl_sub(b, amod, p), qinv, p); /* != 0 */
  if (d >= 1 + (p>>1))
    ax = subii(a, mului(p-d, q));
  else
  {
    ax = addii(a, mului(d, q)); /* in ]0, pq[ assuming a in ]-q,q[ */
    if (cmpii(ax,pq2) > 0) ax = subii(ax,pq);
  }
  return gerepileuptoint(av, ax);
}
GEN
Z_init_CRT(ulong Hp, ulong p) { return stoi(Fl_center(Hp, p, p>>1)); }
GEN
ZX_init_CRT(GEN Hp, ulong p, long v)
{
  long i, l = lg(Hp), lim = (long)(p>>1);
  GEN H = cgetg(l, t_POL);
  H[1] = evalsigne(1) | evalvarn(v);
  for (i=2; i<l; i++)
    gel(H,i) = stoi(Fl_center(Hp[i], p, lim));
  return ZX_renormalize(H,l);
}

GEN
ZM_init_CRT(GEN Hp, ulong p)
{
  long i,j, m, l = lg(Hp), lim = (long)(p>>1);
  GEN c, cp, H = cgetg(l, t_MAT);
  if (l==1) return H;
  m = lgcols(Hp);
  for (j=1; j<l; j++)
  {
    cp = gel(Hp,j);
    c = cgetg(m, t_COL);
    gel(H,j) = c;
    for (i=1; i<m; i++) gel(c,i) = stoi(Fl_center(cp[i],p, lim));
  }
  return H;
}

int
Z_incremental_CRT(GEN *H, ulong Hp, GEN *ptq, ulong p)
{
  GEN h, q = *ptq, qp = muliu(q,p);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  int stable = 1;
  h = Fl_chinese_coprime(*H,Hp,q,p,qinv,qp,shifti(qp,-1));
  if (h) { *H = h; stable = 0; }
  *ptq = qp; return stable;
}

static int
ZX_incremental_CRT_raw(GEN *ptH, GEN Hp, GEN q, GEN qp, ulong p)
{
  GEN H = *ptH, h, qp2 = shifti(qp,-1);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  long i, l = lg(H), lp = lg(Hp);
  int stable = 1;

  if (l < lp)
  { /* degree increases */
    GEN x = cgetg(lp, t_POL);
    for (i=1; i<l; i++)  x[i] = H[i];
    for (   ; i<lp; i++) gel(x,i) = gen_0;
    *ptH = H = x;
    stable = 0;
  } else if (l > lp)
  { /* degree decreases */
    GEN x = cgetg(l, t_VECSMALL);
    for (i=1; i<lp; i++)  x[i] = Hp[i];
    for (   ; i<l; i++) x[i] = 0;
    Hp = x; lp = l;
  }
  for (i=2; i<lp; i++)
  {
    h = Fl_chinese_coprime(gel(H,i),Hp[i],q,p,qinv,qp,qp2);
    if (h) { gel(H,i) = h; stable = 0; }
  }
  (void)ZX_renormalize(H,lp);
  return stable;
}

int
ZX_incremental_CRT(GEN *ptH, GEN Hp, GEN *ptq, ulong p)
{
  GEN q = *ptq, qp = muliu(q,p);
  int stable = ZX_incremental_CRT_raw(ptH, Hp, q, qp, p);
  *ptq = qp; return stable;
}

int
ZM_incremental_CRT(GEN *pH, GEN Hp, GEN *ptq, ulong p)
{
  GEN h, H = *pH, q = *ptq, qp = muliu(q, p), qp2 = shifti(qp,-1);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  long i,j, l = lg(H), m = lgcols(H);
  int stable = 1;
  for (j=1; j<l; j++)
    for (i=1; i<m; i++)
    {
      h = Fl_chinese_coprime(gcoeff(H,i,j), coeff(Hp,i,j),q,p,qinv,qp,qp2);
      if (h) { gcoeff(H,i,j) = h; stable = 0; }
    }
  *ptq = qp; return stable;
}

GEN
ZXM_init_CRT(GEN Hp, long deg, ulong p)
{
  long i, j, k;
  GEN H;
  long m, l = lg(Hp), lim = (long)(p>>1), n;
  H = cgetg(l, t_MAT);
  if (l==1) return H;
  m = lgcols(Hp);
  n = deg + 3;
  for (j=1; j<l; j++)
  {
    GEN cp = gel(Hp,j);
    GEN c = cgetg(m, t_COL);
    gel(H,j) = c;
    for (i=1; i<m; i++)
    {
      GEN dp = gel(cp, i);
      long l = lg(dp);
      GEN d = cgetg(n, t_POL);
      gel(c, i) = d;
      d[1] = dp[1] | evalsigne(1);
      for (k=2; k<l; k++)
        gel(d,k) = stoi(Fl_center(dp[k], p, lim));
      for (   ; k<n; k++)
        gel(d,k) = gen_0;
    }
  }
  return H;
}

int
ZXM_incremental_CRT(GEN *pH, GEN Hp, GEN *ptq, ulong p)
{
  GEN v, H = *pH, q = *ptq, qp = muliu(q, p), qp2 = shifti(qp,-1);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  long i,j,k, l = lg(H), m = lgcols(H), n = lg(gmael(H,1,1));
  int stable = 1;
  for (j=1; j<l; j++)
    for (i=1; i<m; i++)
    {
      GEN h = gmael(H,j,i), hp = gmael(Hp,j,i);
      long lh = lg(hp);
      for (k=2; k<lh; k++)
      {
        v = Fl_chinese_coprime(gel(h,k),uel(hp,k),q,p,qinv,qp,qp2);
        if (v) { gel(h,k) = v; stable = 0; }
      }
      for (; k<n; k++)
      {
        v = Fl_chinese_coprime(gel(h,k),0,q,p,qinv,qp,qp2);
        if (v) { gel(h,k) = v; stable = 0; }
      }
    }
  *ptq = qp; return stable;
}

/* record the degrees of Euclidean remainders (make them as large as
 * possible : smaller values correspond to a degenerate sequence) */
static void
Flx_resultant_set_dglist(GEN a, GEN b, GEN dglist, ulong p)
{
  long da,db,dc, ind;
  pari_sp av = avma;

  if (lgpol(a)==0 || lgpol(b)==0) return;
  da = degpol(a);
  db = degpol(b);
  if (db > da)
  { swapspec(a,b, da,db); }
  else if (!da) return;
  ind = 0;
  while (db)
  {
    GEN c = Flx_rem(a,b, p);
    a = b; b = c; dc = degpol(c);
    if (dc < 0) break;

    ind++;
    if (dc > dglist[ind]) dglist[ind] = dc;
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"Flx_resultant_all");
      gerepileall(av, 2, &a,&b);
    }
    db = dc; /* = degpol(b) */
  }
  if (ind+1 > lg(dglist)) setlg(dglist,ind+1);
  set_avma(av);
}
/* assuming the PRS finishes on a degree 1 polynomial C0 + C1X, with
 * "generic" degree sequence as given by dglist, set *Ci and return
 * resultant(a,b). Modular version of Collins's subresultant */
static ulong
Flx_resultant_all(GEN a, GEN b, long *C0, long *C1, GEN dglist, ulong p)
{
  long da,db,dc, ind;
  ulong lb, res, g = 1UL, h = 1UL, ca = 1UL, cb = 1UL;
  int s = 1;
  pari_sp av = avma;

  *C0 = 1; *C1 = 0;
  if (lgpol(a)==0 || lgpol(b)==0) return 0;
  da = degpol(a);
  db = degpol(b);
  if (db > da)
  {
    swapspec(a,b, da,db);
    if (both_odd(da,db)) s = -s;
  }
  else if (!da) return 1; /* = a[2] ^ db, since 0 <= db <= da = 0 */
  ind = 0;
  while (db)
  { /* sub-resultant algo., applied to ca * a and cb * b, ca,cb scalars,
     * da = deg a, db = deg b */
    GEN c = Flx_rem(a,b, p);
    long delta = da - db;

    if (both_odd(da,db)) s = -s;
    lb = Fl_mul(b[db+2], cb, p);
    a = b; b = c; dc = degpol(c);
    ind++;
    if (dc != dglist[ind]) return gc_ulong(av,0); /* degenerates */
    if (g == h)
    { /* frequent */
      ulong cc = Fl_mul(ca, Fl_powu(Fl_div(lb,g,p), delta+1, p), p);
      ca = cb;
      cb = cc;
    }
    else
    {
      ulong cc = Fl_mul(ca, Fl_powu(lb, delta+1, p), p);
      ulong ghdelta = Fl_mul(g, Fl_powu(h, delta, p), p);
      ca = cb;
      cb = Fl_div(cc, ghdelta, p);
    }
    da = db; /* = degpol(a) */
    db = dc; /* = degpol(b) */

    g = lb;
    if (delta == 1)
      h = g; /* frequent */
    else
      h = Fl_mul(h, Fl_powu(Fl_div(g,h,p), delta, p), p);

    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"Flx_resultant_all");
      gerepileall(av, 2, &a,&b);
    }
  }
  if (da > 1) return 0; /* Failure */
  /* last nonconstant polynomial has degree 1 */
  *C0 = Fl_mul(ca, a[2], p);
  *C1 = Fl_mul(ca, a[3], p);
  res = Fl_mul(cb, b[2], p);
  if (s == -1) res = p - res;
  return gc_ulong(av,res);
}

/* Q a vector of polynomials representing B in Fp[X][Y], evaluate at X = x,
 * Return 0 in case of degree drop. */
static GEN
FlxY_evalx_drop(GEN Q, ulong x, ulong p)
{
  GEN z;
  long i, lb = lg(Q);
  ulong leadz = Flx_eval(leading_coeff(Q), x, p);
  long vs=mael(Q,2,1);
  if (!leadz) return zero_Flx(vs);

  z = cgetg(lb, t_VECSMALL); z[1] = vs;
  for (i=2; i<lb-1; i++) z[i] = Flx_eval(gel(Q,i), x, p);
  z[i] = leadz; return z;
}

GEN
FpXY_FpXQ_evaly(GEN Q, GEN y, GEN T, GEN p, long vx)
{
  pari_sp av = avma;
  long i, lb = lg(Q);
  GEN z;
  if (lb == 2) return pol_0(vx);
  z = gel(Q, lb-1);
  if (lb == 3 || !signe(y)) return typ(z)==t_INT? scalar_ZX(z, vx): ZX_copy(z);

  if (typ(z) == t_INT) z = scalar_ZX_shallow(z, vx);
  for (i=lb-2; i>=2; i--)
  {
    GEN c = gel(Q,i);
    z = FqX_Fq_mul(z, y, T, p);
    z = typ(c) == t_INT? FqX_Fq_add(z,c,T,p): FqX_add(z,c,T,p);
  }
  return gerepileupto(av, z);
}

static GEN
ZX_norml1(GEN x)
{
  long i, l = lg(x);
  GEN s;

  if (l == 2) return gen_0;
  s = gel(x, l-1); /* != 0 */
  for (i = l-2; i > 1; i--) {
    GEN xi = gel(x,i);
    if (!signe(xi)) continue;
    s = addii_sign(s,1, xi,1);
  }
  return s;
}
/* x >= 0, y != 0, return x + |y| */
static GEN
addii_abs(GEN x, GEN y)
{
  if (!signe(x)) return absi_shallow(y);
  return addii_sign(x,1, y,1);
}

/* x a ZX, return sum_{i >= k} |x[i]| binomial(i, k) */
static GEN
ZX_norml1_1(GEN x, long k)
{
  long i, d = degpol(x);
  GEN s, C; /* = binomial(i, k) */

  if (!d || k > d) return gen_0;
  s = absi_shallow(gel(x, k+2)); /* may be 0 */
  C = gen_1;
  for (i = k+1; i <= d; i++) {
    GEN xi = gel(x,i+2);
    if (k) C = diviuexact(muliu(C, i), i-k);
    if (signe(xi)) s = addii_abs(s, mulii(C, xi));
  }
  return s;
}
/* x has non-negative real coefficients */
static GEN
RgX_norml1_1(GEN x, long k)
{
  long i, d = degpol(x);
  GEN s, C; /* = binomial(i, k) */

  if (!d || k > d) return gen_0;
  s = gel(x, k+2); /* may be 0 */
  C = gen_1;
  for (i = k+1; i <= d; i++) {
    GEN xi = gel(x,i+2);
    if (k) C = diviuexact(muliu(C, i), i-k);
    if (!gequal0(xi)) s = gadd(s, gmul(C, xi));
  }
  return s;
}

/* N_2(A)^2 */
static GEN
sqrN2(GEN A, long prec)
{
  pari_sp av = avma;
  long i, l = lg(A);
  GEN a = gen_0;
  for (i = 2; i < l; i++)
  {
    a = gadd(a, gabs(gnorm(gel(A,i)), prec));
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"RgX_RgXY_ResBound i = %ld",i);
      a = gerepileupto(av, a);
    }
  }
  return a;
}
/* Interpolate at roots of 1 and use Hadamard bound for univariate resultant:
 *   bound = N_2(A)^degpol B N_2(B)^degpol(A),  where
 *     N_2(A) = sqrt(sum (N_1(Ai))^2)
 * Return e such that Res(A, B) < 2^e */
static GEN
RgX_RgXY_ResBound(GEN A, GEN B, long prec)
{
  pari_sp av = avma;
  GEN b = gen_0, bnd;
  long i, lB = lg(B);
  for (i=2; i<lB; i++)
  {
    GEN t = gel(B,i);
    if (typ(t) == t_POL) t = gnorml1(t, prec);
    b = gadd(b, gabs(gsqr(t), prec));
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"RgX_RgXY_ResBound i = %ld",i);
      b = gerepileupto(av, b);
    }
  }
  bnd = gsqrt(gmul(gpowgs(sqrN2(A,prec), degpol(B)),
                   gpowgs(b, degpol(A))), prec);
  return gerepileupto(av, bnd);
}
/* A,B in C[X] return RgX_RgXY_ResBound(A, B(x+y)) */
static GEN
RgX_RgXY_ResBound_1(GEN A, GEN B, long prec)
{
  pari_sp av = avma, av2;
  GEN b = gen_0, bnd;
  long i, lB = lg(B);
  B = shallowcopy(B);
  for (i=2; i<lB; i++) gel(B,i) = gabs(gel(B,i), prec);
  av2 = avma;
  for (i=2; i<lB; i++)
  {
    b = gadd(b, gsqr(RgX_norml1_1(B, i-2)));
    if (gc_needed(av2,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"RgX_RgXY_ResBound i = %ld",i);
      b = gerepileupto(av2, b);
    }
  }
  bnd = gsqrt(gmul(gpowgs(sqrN2(A,prec), degpol(B)),
                   gpowgs(b, degpol(A))), prec);
  return gerepileupto(av, bnd);
}

/* log2 N_2(A)^2 */
static double
log2N2(GEN A)
{
  pari_sp av = avma;
  long i, l = lg(A);
  GEN a = gen_0;
  for (i=2; i < l; i++)
  {
    a = addii(a, sqri(gel(A,i)));
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZX_ZXY_ResBound i = %ld",i);
      a = gerepileupto(av, a);
    }
  }
  return gc_double(av, dbllog2(a));
}
/* Interpolate at roots of 1 and use Hadamard bound for univariate resultant:
 *   bound = N_2(A)^degpol B N_2(B)^degpol(A),  where
 *     N_2(A) = sqrt(sum (N_1(Ai))^2)
 * Return e such that Res(A, B) < 2^e */
ulong
ZX_ZXY_ResBound(GEN A, GEN B, GEN dB)
{
  pari_sp av = avma;
  GEN b = gen_0;
  long i, lB = lg(B);
  double logb;
  for (i=2; i<lB; i++)
  {
    GEN t = gel(B,i);
    if (typ(t) == t_POL) t = ZX_norml1(t);
    b = addii(b, sqri(t));
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZX_ZXY_ResBound i = %ld",i);
      b = gerepileupto(av, b);
    }
  }
  logb = dbllog2(b); if (dB) logb -= 2 * dbllog2(dB);
  i = (long)((degpol(B) * log2N2(A) + degpol(A) * logb) / 2);
  return gc_ulong(av, (i <= 0)? 1: 1 + (ulong)i);
}
/* A,B ZX. Return ZX_ZXY_ResBound(A(x), B(x+y)) */
static ulong
ZX_ZXY_ResBound_1(GEN A, GEN B)
{
  pari_sp av = avma;
  GEN b = gen_0;
  long i, lB = lg(B);
  for (i=2; i<lB; i++)
  {
    b = addii(b, sqri(ZX_norml1_1(B, i-2)));
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZX_ZXY_ResBound i = %ld",i);
      b = gerepileupto(av, b);
    }
  }
  i = (long)((degpol(B) * log2N2(A) + degpol(A) * dbllog2(b)) / 2);
  return gc_ulong(av, (i <= 0)? 1: 1 + (ulong)i);
}
/* special case B = A' */
static ulong
ZX_discbound(GEN A)
{
  pari_sp av = avma;
  GEN a = gen_0, b = gen_0;
  long i , lA = lg(A), dA = degpol(A);
  double loga, logb;
  for (i = 2; i < lA; i++)
  {
    GEN c = sqri(gel(A,i));
    a = addii(a, c);
    if (i > 2) b = addii(b, mulii(c, sqru(i-2)));
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZX_discbound i = %ld",i);
      gerepileall(av, 2, &a, &b);
    }
  }
  loga = dbllog2(a);
  logb = dbllog2(b); set_avma(av);
  i = (long)(((dA-1) * loga + dA * logb) / 2);
  return (i <= 0)? 1: 1 + (ulong)i;
}

/* return Res(a(Y), b(n,Y)) over Fp. la = leading_coeff(a) [for efficiency] */
static ulong
Flx_FlxY_eval_resultant(GEN a, GEN b, ulong n, ulong p, ulong pi, ulong la)
{
  GEN ev = FlxY_evalx_pre(b, n, p, pi);
  long drop = lg(b) - lg(ev);
  ulong r = Flx_resultant_pre(a, ev, p, pi);
  if (drop && la != 1) r = Fl_mul(r, Fl_powu_pre(la, drop, p, pi), p);
  return r;
}
static GEN
FpX_FpXY_eval_resultant(GEN a, GEN b, GEN n, GEN p, GEN la, long db, long vX)
{
  GEN ev = FpXY_evaly(b, n, p, vX);
  long drop = db-degpol(ev);
  GEN r = FpX_resultant(a, ev, p);
  if (drop && !gequal1(la)) r = Fp_mul(r, Fp_powu(la, drop,p),p);
  return r;
}

/* assume dres := deg(Res_X(a,b), Y) <= deg(a,X) * deg(b,Y) < p */
/* Return a Fly */
static GEN
Flx_FlxY_resultant_polint(GEN a, GEN b, ulong p, ulong pi, long dres, long sx)
{
  long i;
  ulong n, la = Flx_lead(a);
  GEN  x = cgetg(dres+2, t_VECSMALL);
  GEN  y = cgetg(dres+2, t_VECSMALL);
 /* Evaluate at dres+ 1 points: 0 (if dres even) and +/- n, so that P_n(X) =
  * P_{-n}(-X), where P_i is Lagrange polynomial: P_i(j) = delta_{i,j} */
  for (i=0,n = 1; i < dres; n++)
  {
    x[++i] = n;   y[i] = Flx_FlxY_eval_resultant(a,b, x[i], p,pi,la);
    x[++i] = p-n; y[i] = Flx_FlxY_eval_resultant(a,b, x[i], p,pi,la);
  }
  if (i == dres)
  {
    x[++i] = 0;   y[i] = Flx_FlxY_eval_resultant(a,b, x[i], p,pi,la);
  }
  return Flv_polint(x,y, p, sx);
}

static GEN
FlxX_pseudorem(GEN x, GEN y, ulong p, ulong pi)
{
  long vx = varn(x), dx, dy, dz, i, lx, dp;
  pari_sp av = avma, av2;

  if (!signe(y)) pari_err_INV("FlxX_pseudorem",y);
  (void)new_chunk(2);
  dx=degpol(x); x = RgX_recip_i(x)+2;
  dy=degpol(y); y = RgX_recip_i(y)+2; dz=dx-dy; dp = dz+1;
  av2 = avma;
  for (;;)
  {
    gel(x,0) = Flx_neg(gel(x,0), p); dp--;
    for (i=1; i<=dy; i++)
      gel(x,i) = Flx_add( Flx_mul_pre(gel(y,0), gel(x,i), p, pi),
                          Flx_mul_pre(gel(x,0), gel(y,i), p, pi), p );
    for (   ; i<=dx; i++)
      gel(x,i) = Flx_mul_pre(gel(y,0), gel(x,i), p, pi);
    do { x++; dx--; } while (dx >= 0 && lg(gel(x,0))==2);
    if (dx < dy) break;
    if (gc_needed(av2,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"FlxX_pseudorem dx = %ld >= %ld",dx,dy);
      gerepilecoeffs(av2,x,dx+1);
    }
  }
  if (dx < 0) return zero_Flx(0);
  lx = dx+3; x -= 2;
  x[0]=evaltyp(t_POL) | evallg(lx);
  x[1]=evalsigne(1) | evalvarn(vx);
  x = RgX_recip_i(x);
  if (dp)
  { /* multiply by y[0]^dp   [beware dummy vars from FpX_FpXY_resultant] */
    GEN t = Flx_powu_pre(gel(y,0), dp, p, pi);
    for (i=2; i<lx; i++) gel(x,i) = Flx_mul_pre(gel(x,i), t, p, pi);
  }
  return gerepilecopy(av, x);
}

/* return a Flx */
GEN
FlxX_resultant(GEN u, GEN v, ulong p, long sx)
{
  pari_sp av = avma, av2;
  long degq, dx, dy, du, dv, dr, signh;
  ulong pi;
  GEN z, g, h, r, p1;

  dx = degpol(u); dy = degpol(v); signh = 1;
  if (dx < dy)
  {
    swap(u,v); lswap(dx,dy);
    if (both_odd(dx, dy)) signh = -signh;
  }
  if (dy < 0) return zero_Flx(sx);
  pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  if (dy==0) return gerepileupto(av, Flx_powu_pre(gel(v,2),dx,p,pi));

  g = h = pol1_Flx(sx); av2 = avma;
  for(;;)
  {
    r = FlxX_pseudorem(u,v,p,pi); dr = lg(r);
    if (dr == 2) { set_avma(av); return zero_Flx(sx); }
    du = degpol(u); dv = degpol(v); degq = du-dv;
    u = v; p1 = g; g = leading_coeff(u);
    switch(degq)
    {
      case 0: break;
      case 1:
        p1 = Flx_mul_pre(h,p1, p, pi); h = g; break;
      default:
        p1 = Flx_mul_pre(Flx_powu_pre(h,degq,p,pi), p1, p, pi);
        h = Flx_div_pre(Flx_powu_pre(g,degq,p,pi),
                        Flx_powu_pre(h,degq-1,p,pi), p, pi);
    }
    if (both_odd(du,dv)) signh = -signh;
    v = FlxY_Flx_div(r, p1, p);
    if (dr==3) break;
    if (gc_needed(av2,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"FlxX_resultant, dr = %ld",dr);
      gerepileall(av2,4, &u, &v, &g, &h);
    }
  }
  z = gel(v,2);
  if (dv > 1) z = Flx_div_pre(Flx_powu_pre(z,dv,p,pi),
                              Flx_powu_pre(h,dv-1,p,pi), p, pi);
  if (signh < 0) z = Flx_neg(z,p);
  return gerepileupto(av, z);
}

/* Warning:
 * This function switches between valid and invalid variable ordering*/

static GEN
FlxY_to_FlyX(GEN b, long sv)
{
  long i, n=-1;
  long sw = b[1]&VARNBITS;
  for(i=2;i<lg(b);i++) n = maxss(n,lgpol(gel(b,i)));
  return Flm_to_FlxX(Flm_transpose(FlxX_to_Flm(b,n)),sv,sw);
}

/* Return a Fly*/
GEN
Flx_FlxY_resultant(GEN a, GEN b, ulong p)
{
  pari_sp ltop=avma;
  long dres = degpol(a)*degpol(b);
  long sx=a[1], sy=b[1]&VARNBITS;
  GEN z;
  b = FlxY_to_FlyX(b,sx);
  if ((ulong)dres >= p)
    z = FlxX_resultant(Fly_to_FlxY(a, sy), b, p, sx);
  else
  {
    ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
    z = Flx_FlxY_resultant_polint(a, b, p, pi, (ulong)dres, sy);
  }
  return gerepileupto(ltop,z);
}

/* return a t_POL (in variable v >= 0) whose coeffs are the coeffs of b,
 * in variable v. This is an incorrect PARI object if initially varn(b) << v.
 * We could return a vector of coeffs, but it is convenient to have degpol()
 * and friends available. Even in that case, it will behave nicely with all
 * functions treating a polynomial as a vector of coeffs (eg poleval).
 * FOR INTERNAL USE! */
GEN
swap_vars(GEN b0, long v)
{
  long i, n = RgX_degree(b0, v);
  GEN b, x;
  if (n < 0) return pol_0(v);
  b = cgetg(n+3, t_POL); x = b + 2;
  b[1] = evalsigne(1) | evalvarn(v);
  for (i=0; i<=n; i++) gel(x,i) = polcoef_i(b0, i, v);
  return b;
}

/* assume varn(b) << varn(a) */
/* return a FpY*/
GEN
FpX_FpXY_resultant(GEN a, GEN b, GEN p)
{
  long i,n,dres, db, vY = varn(b), vX = varn(a);
  GEN la,x,y;

  if (lgefint(p) == 3)
  {
    ulong pp = uel(p,2);
    b = ZXX_to_FlxX(b, pp, vX);
    a = ZX_to_Flx(a, pp);
    x = Flx_FlxY_resultant(a, b, pp);
    return Flx_to_ZX(x);
  }
  db = RgXY_degreex(b);
  dres = degpol(a)*degpol(b);
  la = leading_coeff(a);
  x = cgetg(dres+2, t_VEC);
  y = cgetg(dres+2, t_VEC);
 /* Evaluate at dres+ 1 points: 0 (if dres even) and +/- n, so that P_n(X) =
  * P_{-n}(-X), where P_i is Lagrange polynomial: P_i(j) = delta_{i,j} */
  for (i=0,n = 1; i < dres; n++)
  {
    gel(x,++i) = utoipos(n);
    gel(y,i) = FpX_FpXY_eval_resultant(a,b,gel(x,i),p,la,db,vY);
    gel(x,++i) = subiu(p,n);
    gel(y,i) = FpX_FpXY_eval_resultant(a,b,gel(x,i),p,la,db,vY);
  }
  if (i == dres)
  {
    gel(x,++i) = gen_0;
    gel(y,i) = FpX_FpXY_eval_resultant(a,b, gel(x,i), p,la,db,vY);
  }
  return FpV_polint(x,y, p, vY);
}

GEN
FpX_composedsum(GEN P, GEN Q, GEN p)
{
  pari_sp av = avma;
  if (lgefint(p)==3)
  {
    ulong pp = p[2];
    GEN z = Flx_composedsum(ZX_to_Flx(P, pp), ZX_to_Flx(Q, pp), pp);
    return gerepileupto(av, Flx_to_ZX(z));
  }
  else
  {
    long n = 1+ degpol(P)*degpol(Q);
    GEN Pl = FpX_invLaplace(FpX_Newton(P,n,p), p);
    GEN Ql = FpX_invLaplace(FpX_Newton(Q,n,p), p);
    GEN L = FpX_Laplace(FpXn_mul(Pl, Ql, n, p), p);
    GEN lead = Fp_mul(Fp_powu(leading_coeff(P),degpol(Q), p),
        Fp_powu(leading_coeff(Q),degpol(P), p), p);
    GEN R = FpX_fromNewton(L, p);
    return gerepileupto(av, FpX_Fp_mul(R, lead, p));
  }
}

GEN
FpX_composedprod(GEN P, GEN Q, GEN p)
{
  pari_sp av = avma;
  if (lgefint(p)==3)
  {
    ulong pp = p[2];
    GEN z = Flx_composedprod(ZX_to_Flx(P, pp), ZX_to_Flx(Q, pp), pp);
    return gerepileupto(av, Flx_to_ZX(z));
  }
  else
  {
    long n = 1+ degpol(P)*degpol(Q);
    GEN L = FpX_convol(FpX_Newton(P,n,p), FpX_Newton(Q,n,p), p);
    return gerepileupto(av,FpX_fromNewton(L, p));
  }
}

static GEN
_FpX_composedsum(void *E, GEN a, GEN b)
{ return FpX_composedsum(a,b, (GEN)E); }

GEN
FpXV_composedsum(GEN V, GEN p)
{
  if (lgefint(p)==3)
  {
    ulong pp = p[2];
    return Flx_to_ZX(FlxV_composedsum(ZXV_to_FlxV(V, pp), pp));
  }
  return gen_product(V, (void *)p, &_FpX_composedsum);
}

/* 0, 1, -1, 2, -2, ... */
#define next_lambda(a) (a>0 ? -a : 1-a)

/* Assume A in Z[Y], B in Q[Y][X], B squarefree in (Q[Y]/(A))[X] and
 * Res_Y(A, B) in Z[X]. Find a small lambda (start from *lambda, use
 * next_lambda successively) such that C(X) = Res_Y(A(Y), B(X + lambda Y))
 * is squarefree, reset *lambda to the chosen value and return C. Set LERS to
 * the Last nonconstant polynomial in the Euclidean Remainder Sequence */
static GEN
ZX_ZXY_resultant_LERS(GEN A, GEN B0, long *plambda, GEN *LERS)
{
  ulong bound, dp;
  pari_sp av = avma, av2 = 0;
  long lambda = *plambda, degA = degpol(A), dres = degA*degpol(B0);
  long stable, checksqfree, i,n, cnt, degB;
  long v, vX = varn(B0), vY = varn(A); /* vY < vX */
  GEN x, y, dglist, B, q, a, b, ev, H, H0, H1, Hp, H0p, H1p, C0, C1;
  forprime_t S;

  if (degA == 1)
  {
    GEN a1 = gel(A,3), a0 = gel(A,2);
    B = lambda? RgX_translate(B0, monomial(stoi(lambda), 1, vY)): B0;
    H = gsubst(B, vY, gdiv(gneg(a0),a1));
   if (!equali1(a1)) H = RgX_Rg_mul(H, powiu(a1, poldegree(B,vY)));
    *LERS = mkvec2(scalarpol_shallow(a0,vX), scalarpol_shallow(a1,vX));
    return gc_all(av, 2, &H, LERS);
  }

  dglist = Hp = H0p = H1p = C0 = C1 = NULL; /* gcc -Wall */
  C0 = cgetg(dres+2, t_VECSMALL);
  C1 = cgetg(dres+2, t_VECSMALL);
  dglist = cgetg(dres+1, t_VECSMALL);
  x = cgetg(dres+2, t_VECSMALL);
  y = cgetg(dres+2, t_VECSMALL);
  B0 = leafcopy(B0);
  A = leafcopy(A);
  B = B0;
  v = fetch_var_higher(); setvarn(A,v);
  /* make sure p large enough */
INIT:
  /* always except the first time */
  if (av2) { set_avma(av2); lambda = next_lambda(lambda); }
  if (lambda) B = RgX_translate(B0, monomial(stoi(lambda), 1, vY));
  B = swap_vars(B, vY); setvarn(B,v);
  /* B0(lambda v + x, v) */
  if (DEBUGLEVEL>4) err_printf("Trying lambda = %ld\n", lambda);
  av2 = avma;

  if (degA <= 3)
  { /* sub-resultant faster for small degrees */
    H = RgX_resultant_all(A,B,&q);
    if (typ(q) != t_POL || degpol(q)!=1) goto INIT;
    H0 = gel(q,2);
    if (typ(H0) == t_POL) setvarn(H0,vX); else H0 = scalarpol(H0,vX);
    H1 = gel(q,3);
    if (typ(H1) == t_POL) setvarn(H1,vX); else H1 = scalarpol(H1,vX);
    if (!ZX_is_squarefree(H)) goto INIT;
    goto END;
  }

  H = H0 = H1 = NULL;
  degB = degpol(B);
  bound = ZX_ZXY_ResBound(A, B, NULL);
  if (DEBUGLEVEL>4) err_printf("bound for resultant coeffs: 2^%ld\n",bound);
  dp = 1;
  init_modular_big(&S);
  for(cnt = 0, checksqfree = 1;;)
  {
    ulong p = u_forprime_next(&S);
    GEN Hi;
    a = ZX_to_Flx(A, p);
    b = ZXX_to_FlxX(B, p, varn(A));
    if (degpol(a) < degA || degpol(b) < degB) continue; /* p | lc(A)lc(B) */
    if (checksqfree)
    { /* find degree list for generic Euclidean Remainder Sequence */
      long goal = minss(degpol(a), degpol(b)); /* longest possible */
      for (n=1; n <= goal; n++) dglist[n] = 0;
      setlg(dglist, 1);
      for (n=0; n <= dres; n++)
      {
        ev = FlxY_evalx_drop(b, n, p);
        Flx_resultant_set_dglist(a, ev, dglist, p);
        if (lg(dglist)-1 == goal) break;
      }
      /* last pol in ERS has degree > 1 ? */
      goal = lg(dglist)-1;
      if (degpol(B) == 1) { if (!goal) goto INIT; }
      else
      {
        if (goal <= 1) goto INIT;
        if (dglist[goal] != 0 || dglist[goal-1] != 1) goto INIT;
      }
      if (DEBUGLEVEL>4)
        err_printf("Degree list for ERS (trials: %ld) = %Ps\n",n+1,dglist);
    }

    for (i=0,n = 0; i <= dres; n++)
    {
      ev = FlxY_evalx_drop(b, n, p);
      x[++i] = n; y[i] = Flx_resultant_all(a, ev, C0+i, C1+i, dglist, p);
      if (!C1[i]) i--; /* C1(i) = 0. No way to recover C0(i) */
    }
    Hi = Flv_Flm_polint(x, mkvec3(y,C0,C1), p, 0);
    Hp = gel(Hi,1); H0p = gel(Hi,2); H1p = gel(Hi,3);
    if (!H && degpol(Hp) != dres) continue;
    if (dp != 1) Hp = Flx_Fl_mul(Hp, Fl_powu(Fl_inv(dp,p), degA, p), p);
    if (checksqfree) {
      if (!Flx_is_squarefree(Hp, p)) goto INIT;
      if (DEBUGLEVEL>4) err_printf("Final lambda = %ld\n", lambda);
      checksqfree = 0;
    }

    if (!H)
    { /* initialize */
      q = utoipos(p); stable = 0;
      H = ZX_init_CRT(Hp, p,vX);
      H0= ZX_init_CRT(H0p, p,vX);
      H1= ZX_init_CRT(H1p, p,vX);
    }
    else
    {
      GEN qp = muliu(q,p);
      stable  = ZX_incremental_CRT_raw(&H, Hp, q,qp, p)
              & ZX_incremental_CRT_raw(&H0,H0p, q,qp, p)
              & ZX_incremental_CRT_raw(&H1,H1p, q,qp, p);
      q = qp;
    }
    /* could make it probabilistic for H ? [e.g if stable twice, etc]
     * Probabilistic anyway for H0, H1 */
    if (DEBUGLEVEL>5 && (stable ||  ++cnt==100))
    { cnt=0; err_printf("%ld%%%s ",100*expi(q)/bound,stable?"s":""); }
    if (stable && (ulong)expi(q) >= bound) break; /* DONE */
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"ZX_ZXY_rnfequation");
      gerepileall(av2, 4, &H, &q, &H0, &H1);
    }
  }
END:
  if (DEBUGLEVEL>5) err_printf(" done\n");
  setvarn(H, vX); (void)delete_var();
  *LERS = mkvec2(H0,H1);
  *plambda = lambda; return gc_all(av, 2, &H, LERS);
}

GEN
ZX_ZXY_resultant_all(GEN A, GEN B, long *plambda, GEN *LERS)
{
  if (LERS)
  {
    if (!plambda)
      pari_err_BUG("ZX_ZXY_resultant_all [LERS != NULL needs lambda]");
    return ZX_ZXY_resultant_LERS(A, B, plambda, LERS);
  }
  return ZX_ZXY_rnfequation(A, B, plambda);
}

/* If lambda = NULL, return caract(Mod(A, T)), T,A in Z[X].
 * Otherwise find a small lambda such that caract (Mod(A + lambda X, T)) is
 * squarefree */
GEN
ZXQ_charpoly_sqf(GEN A, GEN T, long *lambda, long v)
{
  pari_sp av = avma;
  GEN R, a;
  long dA;
  int delvar;

  if (v < 0) v = 0;
  switch (typ(A))
  {
    case t_POL: dA = degpol(A); if (dA > 0) break;
      A = constant_coeff(A);
    default:
      if (lambda) { A = scalar_ZX_shallow(A,varn(T)); dA = 0; break;}
      return gerepileupto(av, gpowgs(gsub(pol_x(v), A), degpol(T)));
  }
  delvar = 0;
  if (varn(T) == 0)
  {
    long v0 = fetch_var(); delvar = 1;
    T = leafcopy(T); setvarn(T,v0);
    A = leafcopy(A); setvarn(A,v0);
  }
  R = ZX_ZXY_rnfequation(T, deg1pol_shallow(gen_1, gneg_i(A), 0), lambda);
  if (delvar) (void)delete_var();
  setvarn(R, v); a = leading_coeff(T);
  if (!gequal1(a)) R = gdiv(R, powiu(a, dA));
  return gerepileupto(av, R);
}

/* charpoly(Mod(A,T)), A may be in Q[X], but assume T and result are integral */
GEN
ZXQ_charpoly(GEN A, GEN T, long v)
{
  return (degpol(T) < 16) ? RgXQ_charpoly_i(A,T,v): ZXQ_charpoly_sqf(A,T, NULL, v);
}

GEN
QXQ_charpoly(GEN A, GEN T, long v)
{
  pari_sp av = avma;
  GEN den, B = Q_remove_denom(A, &den);
  GEN P = ZXQ_charpoly(B, T, v);
  return gerepilecopy(av, den ? RgX_rescale(P, ginv(den)): P);
}

static ulong
ZX_resultant_prime(GEN a, GEN b, GEN dB, long degA, long degB, ulong p)
{
  long dropa = degA - degpol(a), dropb = degB - degpol(b);
  ulong H, dp;
  if (dropa && dropb) return 0; /* p | lc(A), p | lc(B) */
  H = Flx_resultant(a, b, p);
  if (dropa)
  { /* multiply by ((-1)^deg B lc(B))^(deg A - deg a) */
    ulong c = b[degB+2]; /* lc(B) */
    if (odd(degB)) c = p - c;
    c = Fl_powu(c, dropa, p);
    if (c != 1) H = Fl_mul(H, c, p);
  }
  else if (dropb)
  { /* multiply by lc(A)^(deg B - deg b) */
    ulong c = a[degA+2]; /* lc(A) */
    c = Fl_powu(c, dropb, p);
    if (c != 1) H = Fl_mul(H, c, p);
  }
  dp = dB ? umodiu(dB, p): 1;
  if (dp != 1) H = Fl_mul(H, Fl_powu(Fl_inv(dp,p), degA, p), p);
  return H;
}

/* If B=NULL, assume B=A' */
static GEN
ZX_resultant_slice(GEN A, GEN B, GEN dB, GEN P, GEN *mod)
{
  pari_sp av = avma, av2;
  long degA, degB, i, n = lg(P)-1;
  GEN H, T;

  degA = degpol(A);
  degB = B? degpol(B): degA - 1;
  if (n == 1)
  {
    ulong Hp, p = uel(P,1);
    GEN a = ZX_to_Flx(A, p), b = B? ZX_to_Flx(B, p): Flx_deriv(a, p);
    Hp = ZX_resultant_prime(a, b, dB, degA, degB, p);
    set_avma(av); *mod = utoipos(p); return utoi(Hp);
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  if (B) B = ZX_nv_mod_tree(B, P, T);
  H = cgetg(n+1, t_VECSMALL); av2 = avma;
  for(i=1; i <= n; i++, set_avma(av2))
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = B? gel(B,i): Flx_deriv(a, p);
    H[i] = ZX_resultant_prime(a, b, dB, degA, degB, p);
  }
  H = ZV_chinese_tree(H, P, T, ZV_chinesetree(P,T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZX_resultant_worker(GEN P, GEN A, GEN B, GEN dB)
{
  GEN V = cgetg(3, t_VEC);
  if (typ(B) == t_INT) B = NULL;
  if (!signe(dB)) dB = NULL;
  gel(V,1) = ZX_resultant_slice(A, B, dB, P, &gel(V,2));
  return V;
}

/* Compute Res(A, B/dB) in Z, assuming A,B in Z[X], dB in Z or NULL (= 1)
 * If B=NULL, take B = A' and assume deg A > 1 and 'bound' is set */
GEN
ZX_resultant_all(GEN A, GEN B, GEN dB, ulong bound)
{
  pari_sp av = avma;
  forprime_t S;
  GEN  H, worker;
  if (B)
  {
    long a = degpol(A), b = degpol(B);
    if (a < 0 || b < 0) return gen_0;
    if (!a) return powiu(gel(A,2), b);
    if (!b) return powiu(gel(B,2), a);
    if (!bound) bound = ZX_ZXY_ResBound(A, B, dB);
  }
  worker = snm_closure(is_entry("_ZX_resultant_worker"),
                       mkvec3(A, B? B: gen_0, dB? dB: gen_0));
  init_modular_big(&S);
  H = gen_crt("ZX_resultant_all", worker, &S, dB, bound, 0, NULL,
              ZV_chinese_center, Fp_center);
  return gerepileuptoint(av, H);
}

/* A0 and B0 in Q[X] */
GEN
QX_resultant(GEN A0, GEN B0)
{
  GEN s, a, b, A, B;
  pari_sp av = avma;

  A = Q_primitive_part(A0, &a);
  B = Q_primitive_part(B0, &b);
  s = ZX_resultant(A, B);
  if (!signe(s)) { set_avma(av); return gen_0; }
  if (a) s = gmul(s, gpowgs(a,degpol(B)));
  if (b) s = gmul(s, gpowgs(b,degpol(A)));
  return gerepileupto(av, s);
}

GEN
ZX_resultant(GEN A, GEN B) { return ZX_resultant_all(A,B,NULL,0); }

GEN
QXQ_intnorm(GEN A, GEN B)
{
  GEN c, n, R, lB;
  long dA = degpol(A), dB = degpol(B);
  pari_sp av = avma;
  if (dA < 0) return gen_0;
  A = Q_primitive_part(A, &c);
  if (!c || typ(c) == t_INT) {
    n = c;
    R = ZX_resultant(B, A);
  } else {
    n = gel(c,1);
    R = ZX_resultant_all(B, A, gel(c,2), 0);
  }
  if (n && !equali1(n)) R = mulii(R, powiu(n, dB));
  lB = leading_coeff(B);
  if (!equali1(lB)) R = diviiexact(R, powiu(lB, dA));
  return gerepileuptoint(av, R);
}

GEN
QXQ_norm(GEN A, GEN B)
{
  GEN c, R, lB;
  long dA = degpol(A), dB = degpol(B);
  pari_sp av = avma;
  if (dA < 0) return gen_0;
  A = Q_primitive_part(A, &c);
  R = ZX_resultant(B, A);
  if (c) R = gmul(R, gpowgs(c, dB));
  lB = leading_coeff(B);
  if (!equali1(lB)) R = gdiv(R, gpowgs(lB, dA));
  return gerepileupto(av, R);
}

/* assume x has integral coefficients */
GEN
ZX_disc_all(GEN x, ulong bound)
{
  pari_sp av = avma;
  long s, d = degpol(x);
  GEN l, R;

  if (d <= 1) return d == 1? gen_1: gen_0;
  s = (d & 2) ? -1: 1;
  l = leading_coeff(x);
  if (!bound) bound = ZX_discbound(x);
  R = ZX_resultant_all(x, NULL, NULL, bound);
  if (is_pm1(l))
  { if (signe(l) < 0) s = -s; }
  else
    R = diviiexact(R,l);
  if (s == -1) togglesign_safe(&R);
  return gerepileuptoint(av,R);
}

GEN
ZX_disc(GEN x) { return ZX_disc_all(x,0); }

static GEN
ZXQX_resultant_prime(GEN a, GEN b, GEN dB, long degA, long degB, GEN T, ulong p)
{
  long dropa = degA - degpol(a), dropb = degB - degpol(b);
  GEN H, dp;
  if (dropa && dropb) return pol0_Flx(T[1]); /* p | lc(A), p | lc(B) */
  H = FlxqX_saferesultant(a, b, T, p);
  if (!H) return NULL;
  if (dropa)
  { /* multiply by ((-1)^deg B lc(B))^(deg A - deg a) */
    GEN c = gel(b,degB+2); /* lc(B) */
    if (odd(degB)) c = Flx_neg(c, p);
    c = Flxq_powu(c, dropa, T, p);
    if (!Flx_equal1(c)) H = Flxq_mul(H, c, T, p);
  }
  else if (dropb)
  { /* multiply by lc(A)^(deg B - deg b) */
    GEN c = gel(a,degA+2); /* lc(A) */
    c = Flxq_powu(c, dropb, T, p);
    if (!Flx_equal1(c)) H = Flxq_mul(H, c, T, p);
  }
  dp = dB ? ZX_to_Flx(dB, p): pol1_Flx(T[1]);
  if (!Flx_equal1(dp))
  {
    GEN idp = Flxq_invsafe(dp, T, p);
    if (!idp) return NULL;
    H = Flxq_mul(H, Flxq_powu(idp, degA, T, p), T, p);
  }
  return H;
}

/* If B=NULL, assume B=A' */
static GEN
ZXQX_resultant_slice(GEN A, GEN B, GEN U, GEN dB, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long degA, degB, i, n = lg(P)-1;
  GEN H, T;
  long v = varn(U), redo = 0;

  degA = degpol(A);
  degB = B? degpol(B): degA - 1;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZXX_to_FlxX(A, p, v), b = B? ZXX_to_FlxX(B, p, v): FlxX_deriv(a, p);
    GEN u = ZX_to_Flx(U, p);
    GEN Hp = ZXQX_resultant_prime(a, b, dB, degA, degB, u, p);
    if (!Hp) { set_avma(av); *mod = gen_1; return pol_0(v); }
    Hp = gerepileupto(av, Flx_to_ZX(Hp)); *mod = utoipos(p); return Hp;
  }
  T = ZV_producttree(P);
  A = ZXX_nv_mod_tree(A, P, T, v);
  if (B) B = ZXX_nv_mod_tree(B, P, T, v);
  U = ZX_nv_mod_tree(U, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = B? gel(B,i): FlxX_deriv(a, p), u = gel(U, i);
    GEN h = ZXQX_resultant_prime(a, b, dB, degA, degB, u, p);
    if (!h)
    {
      gel(H,i) = pol_0(v);
      P[i] = 1; redo = 1;
    }
    else
      gel(H,i) = h;
  }
  if (redo) T = ZV_producttree(P);
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZXQX_resultant_worker(GEN P, GEN A, GEN B, GEN T, GEN dB)
{
  GEN V = cgetg(3, t_VEC);
  if (isintzero(B)) B = NULL;
  if (!signe(dB)) dB = NULL;
  gel(V,1) = ZXQX_resultant_slice(A, B, T, dB, P, &gel(V,2));
  return V;
}

static ulong
ZXQX_resultant_bound_i(GEN nf, GEN A, GEN B, GEN (*f)(GEN,GEN,long))
{
  pari_sp av = avma;
  GEN r, M = nf_L2_bound(nf, NULL, &r);
  long v = nf_get_varn(nf), i, l = lg(r);
  GEN a = cgetg(l, t_COL);
  for (i = 1; i < l; i++)
    gel(a, i) = f(gsubst(A, v, gel(r,i)), gsubst(B, v, gel(r,i)), DEFAULTPREC);
  return gc_ulong(av, (ulong) dbllog2(gmul(M,RgC_fpnorml2(a, DEFAULTPREC))));
}
static ulong
ZXQX_resultant_bound(GEN nf, GEN A, GEN B)
{ return ZXQX_resultant_bound_i(nf, A, B, &RgX_RgXY_ResBound); }

/* Compute Res(A, B/dB) in Z[X]/T, assuming A,B in Z[X,Y], dB in Z or NULL (= 1)
 * If B=NULL, take B = A' and assume deg A > 1 */
static GEN
ZXQX_resultant_all(GEN A, GEN B, GEN T, GEN dB, ulong bound)
{
  pari_sp av = avma;
  forprime_t S;
  GEN  H, worker;
  if (B)
  {
    long a = degpol(A), b = degpol(B);
    if (a < 0 || b < 0) return gen_0;
    if (!a) return gpowgs(gel(A,2), b);
    if (!b) return gpowgs(gel(B,2), a);
  } else
    if (!bound) B = RgX_deriv(A);
  if (!bound) bound = ZXQX_resultant_bound(nfinit(T, DEFAULTPREC), A, B);
  worker = snm_closure(is_entry("_ZXQX_resultant_worker"),
                       mkvec4(A, B? B: gen_0, T, dB? dB: gen_0));
  init_modular_big(&S);
  H = gen_crt("ZXQX_resultant_all", worker, &S, dB, bound, 0, NULL,
              nxV_chinese_center, FpX_center);
  if (DEBUGLEVEL)
    err_printf("ZXQX_resultant_all: a priori bound: %lu, a posteriori: %lu\n",
               bound, expi(gsupnorm(H, DEFAULTPREC)));
  return gerepileupto(av, H);
}

GEN
nfX_resultant(GEN nf, GEN x, GEN y)
{
  pari_sp av = avma;
  GEN cx, cy, D, T = nf_get_pol(nf);
  ulong bound;
  long d = degpol(x), v = varn(T);
  if (d <= 1) return d == 1? pol_1(v): pol_0(v);
  x = Q_primitive_part(x, &cx);
  y = Q_primitive_part(y, &cy);
  bound = ZXQX_resultant_bound(nf, x, y);
  D = ZXQX_resultant_all(x, y, T, NULL, bound);
  if (cx) D = gmul(D, gpowgs(cx, degpol(y)));
  if (cy) D = gmul(D, gpowgs(cy, degpol(x)));
  return gerepileupto(av, D);
}

static GEN
to_ZX(GEN a, long v) { return typ(a)==t_INT? scalarpol(a,v): a; }

static GEN
ZXQX_disc_all(GEN x, GEN T, ulong bound)
{
  pari_sp av = avma;
  long s, d = degpol(x), v = varn(T);
  GEN l, R;

  if (d <= 1) return d == 1? pol_1(v): pol_0(v);
  s = (d & 2) ? -1: 1;
  l = leading_coeff(x);
  R = ZXQX_resultant_all(x, NULL, T, NULL, bound);
  if (!gequal1(l)) R = QXQ_div(R, to_ZX(l,v), T);
  if (s == -1) R = RgX_neg(R);
  return gerepileupto(av, R);
}

GEN
QX_disc(GEN x)
{
  pari_sp av = avma;
  GEN c, d = ZX_disc( Q_primitive_part(x, &c) );
  if (c) d = gmul(d, gpowgs(c, 2*degpol(x) - 2));
  return gerepileupto(av, d);
}

GEN
nfX_disc(GEN nf, GEN x)
{
  pari_sp av = avma;
  GEN c, D, T = nf_get_pol(nf);
  ulong bound;
  long d = degpol(x), v = varn(T);
  if (d <= 1) return d == 1? pol_1(v): pol_0(v);
  x = Q_primitive_part(x, &c);
  bound = ZXQX_resultant_bound(nf, x, RgX_deriv(x));
  D = ZXQX_disc_all(x, T, bound);
  if (c) D = gmul(D, gpowgs(c, 2*d - 2));
  return gerepileupto(av, D);
}

GEN
QXQ_mul(GEN x, GEN y, GEN T)
{
  GEN dx, nx = Q_primitive_part(x, &dx);
  GEN dy, ny = Q_primitive_part(y, &dy);
  GEN z = ZXQ_mul(nx, ny, T);
  if (dx || dy)
  {
    GEN d = dx ? dy ? gmul(dx, dy): dx : dy;
    if (!gequal1(d)) z = ZX_Q_mul(z, d);
  }
  return z;
}

GEN
QXQ_sqr(GEN x, GEN T)
{
  GEN dx, nx = Q_primitive_part(x, &dx);
  GEN z = ZXQ_sqr(nx, T);
  if (dx)
    z = ZX_Q_mul(z, gsqr(dx));
  return z;
}

static GEN
QXQ_inv_slice(GEN A, GEN B, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1, v = varn(A), redo = 0;
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZX_to_Flx(A, p), b = ZX_to_Flx(B, p);
    GEN U = Flxq_invsafe(a, b, p);
    if (!U)
    {
      set_avma(av);
      *mod = gen_1; return pol_0(v);
    }
    H = gerepilecopy(av, Flx_to_ZX(U));
    *mod = utoipos(p); return H;
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  B = ZX_nv_mod_tree(B, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i);
    GEN U = Flxq_invsafe(a, b, p);
    if (!U)
    {
      gel(H,i) = pol_0(v);
      P[i] = 1; redo = 1;
    }
    else
      gel(H,i) = U;
  }
  if (redo) T = ZV_producttree(P);
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
QXQ_inv_worker(GEN P, GEN A, GEN B)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = QXQ_inv_slice(A, B, P, &gel(V,2));
  return V;
}

/* lift(1 / Mod(A,B)). B a ZX, A a scalar or a QX */
GEN
QXQ_inv(GEN A, GEN B)
{
  GEN D, Ap, Bp;
  ulong pp;
  pari_sp av2, av = avma;
  forprime_t S;
  GEN worker, U, H = NULL, mod = gen_1;
  pari_timer ti;
  long k, dA, dB;
  if (is_scalar_t(typ(A))) return scalarpol(ginv(A), varn(B));
  /* A a QX, B a ZX */
  A = Q_primitive_part(A, &D);
  dA = degpol(A); dB= degpol(B);
  /* A, B in Z[X] */
  init_modular_small(&S);
  do {
    pp = u_forprime_next(&S);
    Ap = ZX_to_Flx(A, pp);
    Bp = ZX_to_Flx(B, pp);
  } while (degpol(Ap) != dA || degpol(Bp) != dB);
  if (degpol(Flx_gcd(Ap, Bp, pp)) != 0 && degpol(ZX_gcd(A,B))!=0)
    pari_err_INV("QXQ_inv",mkpolmod(A,B));
  worker = snm_closure(is_entry("_QXQ_inv_worker"), mkvec2(A, B));
  av2 = avma;
  for (k = 1; ;k *= 2)
  {
    GEN res, b, N, den;
    gen_inccrt_i("QXQ_inv", worker, NULL, (k+1)>>1, 0, &S, &H, &mod,
                 nxV_chinese_center, FpX_center);
    gerepileall(av2, 2, &H, &mod);
    b = sqrti(shifti(mod,-1));
    if (DEBUGLEVEL>5) timer_start(&ti);
    U = FpX_ratlift(H, mod, b, b, NULL);
    if (DEBUGLEVEL>5) timer_printf(&ti,"QXQ_inv: ratlift");
    if (!U) continue;
    N = Q_remove_denom(U, &den); if (!den) den = gen_1;
    res = Flx_rem(Flx_Fl_sub(Flx_mul(Ap, ZX_to_Flx(N,pp), pp),
                  umodiu(den, pp), pp), Bp, pp);
    if (degpol(res) >= 0) continue;
    res = ZX_Z_sub(ZX_mul(A, N), den);
    res = ZX_is_monic(B) ? ZX_rem(res, B): RgX_pseudorem(res, B);
    if (DEBUGLEVEL>5) timer_printf(&ti,"QXQ_inv: final check");
    if (degpol(res)<0)
    {
      if (D) U = RgX_Rg_div(U, D);
      return gerepilecopy(av, U);
    }
  }
}

static GEN
QXQ_div_slice(GEN A, GEN B, GEN C, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1, v = varn(A), redo = 0;
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZX_to_Flx(A, p), b = ZX_to_Flx(B, p), c = ZX_to_Flx(C, p);
    GEN bi = Flxq_invsafe(b, c, p), U;
    if (!bi)
    {
      set_avma(av);
      *mod = gen_1; return pol_0(v);
    }
    U = Flxq_mul(a, bi, c, p);
    H = gerepilecopy(av, Flx_to_ZX(U));
    *mod = utoipos(p); return H;
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  B = ZX_nv_mod_tree(B, P, T);
  C = ZX_nv_mod_tree(C, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i), c = gel(C, i);
    GEN bi = Flxq_invsafe(b, c, p);
    if (!bi)
    {
      gel(H,i) = pol_0(v);
      P[i] = 1; redo = 1;
    }
    else
      gel(H,i) = Flxq_mul(a, bi, c, p);
  }
  if (redo) T = ZV_producttree(P);
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
QXQ_div_worker(GEN P, GEN A, GEN B, GEN C)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = QXQ_div_slice(A, B, C, P, &gel(V,2));
  return V;
}

/* lift(Mod(A/B, C)). C a ZX, A, B a scalar or a QX */
GEN
QXQ_div(GEN A, GEN B, GEN C)
{
  GEN DA, DB, Ap, Bp, Cp;
  ulong pp;
  pari_sp av2, av = avma;
  forprime_t S;
  GEN worker, U, H = NULL, mod = gen_1;
  pari_timer ti;
  long k, dA, dB, dC;
  if (is_scalar_t(typ(A))) return scalarpol(ginv(A), varn(B));
  /* A a QX, B a ZX */
  A = Q_primitive_part(A, &DA);
  B = Q_primitive_part(B, &DB);
  dA = degpol(A); dB = degpol(B); dC = degpol(C);
  /* A, B in Z[X] */
  init_modular_small(&S);
  do {
    pp = u_forprime_next(&S);
    Ap = ZX_to_Flx(A, pp);
    Bp = ZX_to_Flx(B, pp);
    Cp = ZX_to_Flx(C, pp);
  } while (degpol(Ap) != dA || degpol(Bp) != dB || degpol(Cp) != dC);
  if (degpol(Flx_gcd(Bp, Cp, pp)) != 0 && degpol(ZX_gcd(B,C))!=0)
    pari_err_INV("QXQ_div",mkpolmod(B,C));
  worker = snm_closure(is_entry("_QXQ_div_worker"), mkvec3(A, B, C));
  av2 = avma;
  for (k = 1; ;k *= 2)
  {
    GEN res, b, N, den;
    gen_inccrt_i("QXQ_div", worker, NULL, (k+1)>>1, 0, &S, &H, &mod,
                 nxV_chinese_center, FpX_center);
    gerepileall(av2, 2, &H, &mod);
    b = sqrti(shifti(mod,-1));
    if (DEBUGLEVEL>5) timer_start(&ti);
    U = FpX_ratlift(H, mod, b, b, NULL);
    if (DEBUGLEVEL>5) timer_printf(&ti,"QXQ_div: ratlift");
    if (!U) continue;
    N = Q_remove_denom(U, &den); if (!den) den = gen_1;
    res = Flx_rem(Flx_sub(Flx_mul(Bp, ZX_to_Flx(N,pp), pp),
                          Flx_Fl_mul(Ap, umodiu(den, pp), pp), pp), Cp, pp);
    if (degpol(res) >= 0) continue;
    res = ZX_sub(ZX_mul(B, N), ZX_Z_mul(A,den));
    res = ZX_is_monic(C) ? ZX_rem(res, C): RgX_pseudorem(res, C);
    if (DEBUGLEVEL>5) timer_printf(&ti,"QXQ_div: final check");
    if (degpol(res)<0)
    {
      if (DA && DB) U = RgX_Rg_mul(U, gdiv(DA,DB));
      else if (DA) U = RgX_Rg_mul(U, DA);
      else if (DB) U = RgX_Rg_div(U, DB);
      return gerepilecopy(av, U);
    }
  }
}

/************************************************************************
 *                                                                      *
 *                           ZXQ_minpoly                                *
 *                                                                      *
 ************************************************************************/

static GEN
ZXQ_minpoly_slice(GEN A, GEN B, long d, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1, v = evalvarn(varn(B));
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZX_to_Flx(A, p), b = ZX_to_Flx(B, p);
    GEN Hp = Flxq_minpoly(a, b, p);
    if (degpol(Hp) != d) { p = 1; Hp = pol0_Flx(v); }
    H = gerepileupto(av, Flx_to_ZX(Hp));
    *mod = utoipos(p); return H;
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  B = ZX_nv_mod_tree(B, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i);
    GEN m = Flxq_minpoly(a, b, p);
    if (degpol(m) != d) { P[i] = 1; m = pol0_Flx(v); }
    gel(H, i) = m;
  }
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZXQ_minpoly_worker(GEN P, GEN A, GEN B, long d)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = ZXQ_minpoly_slice(A, B, d, P, &gel(V,2));
  return V;
}

GEN
ZXQ_minpoly(GEN A, GEN B, long d, ulong bound)
{
  pari_sp av = avma;
  GEN worker, H, dB;
  forprime_t S;
  B = Q_remove_denom(B, &dB);
  worker = strtoclosure("_ZXQ_minpoly_worker", 3, A, B, stoi(d));
  init_modular_big(&S);
  H = gen_crt("ZXQ_minpoly", worker, &S, dB, bound, 0, NULL,
               nxV_chinese_center, FpX_center_i);
  return gerepilecopy(av, H);
}

/************************************************************************
 *                                                                      *
 *                   ZX_ZXY_resultant                                   *
 *                                                                      *
 ************************************************************************/

static GEN
ZX_ZXY_resultant_prime(GEN a, GEN b, ulong dp, ulong p,
                       long degA, long degB, long dres, long sX)
{
  long dropa = degA - degpol(a), dropb = degB - degpol(b);
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  GEN Hp = Flx_FlxY_resultant_polint(a, b, p, pi, dres, sX);
  if (dropa && dropb)
    Hp = zero_Flx(sX);
  else {
    if (dropa)
    { /* multiply by ((-1)^deg B lc(B))^(deg A - deg a) */
      GEN c = gel(b,degB+2); /* lc(B) */
      if (odd(degB)) c = Flx_neg(c, p);
      if (!Flx_equal1(c)) {
        c = Flx_powu_pre(c, dropa, p, pi);
        if (!Flx_equal1(c)) Hp = Flx_mul_pre(Hp, c, p, pi);
      }
    }
    else if (dropb)
    { /* multiply by lc(A)^(deg B - deg b) */
      ulong c = uel(a, degA+2); /* lc(A) */
      c = Fl_powu(c, dropb, p);
      if (c != 1) Hp = Flx_Fl_mul_pre(Hp, c, p, pi);
    }
  }
  if (dp != 1) Hp = Flx_Fl_mul_pre(Hp, Fl_powu_pre(Fl_inv(dp,p), degA, p, pi), p, pi);
  return Hp;
}

static GEN
ZX_ZXY_resultant_slice(GEN A, GEN B, GEN dB, long degA, long degB, long dres,
                       GEN P, GEN *mod, long sX, long vY)
{
  pari_sp av = avma;
  long i, n = lg(P)-1;
  GEN H, T, D;
  if (n == 1)
  {
    ulong p = uel(P,1);
    ulong dp = dB ? umodiu(dB, p): 1;
    GEN a = ZX_to_Flx(A, p), b = ZXX_to_FlxX(B, p, vY);
    GEN Hp = ZX_ZXY_resultant_prime(a, b, dp, p, degA, degB, dres, sX);
    H = gerepileupto(av, Flx_to_ZX(Hp));
    *mod = utoipos(p); return H;
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  B = ZXX_nv_mod_tree(B, P, T, vY);
  D = dB ? Z_ZV_mod_tree(dB, P, T): NULL;
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i);
    ulong dp = D ? uel(D, i): 1;
    gel(H,i) = ZX_ZXY_resultant_prime(a, b, dp, p, degA, degB, dres, sX);
  }
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZX_ZXY_resultant_worker(GEN P, GEN A, GEN B, GEN dB, GEN v)
{
  GEN V = cgetg(3, t_VEC);
  if (isintzero(dB)) dB = NULL;
  gel(V,1) = ZX_ZXY_resultant_slice(A, B, dB, v[1], v[2], v[3], P, &gel(V,2), v[4], v[5]);
  return V;
}

GEN
ZX_ZXY_resultant(GEN A, GEN B)
{
  pari_sp av = avma;
  forprime_t S;
  ulong bound;
  long v = fetch_var_higher();
  long degA = degpol(A), degB, dres = degA * degpol(B);
  long vX = varn(B), vY = varn(A); /* assume vY has lower priority */
  long sX = evalvarn(vX);
  GEN worker, H, dB;
  B = Q_remove_denom(B, &dB);
  if (!dB) B = leafcopy(B);
  A = leafcopy(A); setvarn(A,v);
  B = swap_vars(B, vY); setvarn(B,v); degB = degpol(B);
  bound = ZX_ZXY_ResBound(A, B, dB);
  if (DEBUGLEVEL>4) err_printf("bound for resultant coeffs: 2^%ld\n",bound);
  worker = snm_closure(is_entry("_ZX_ZXY_resultant_worker"),
                       mkvec4(A, B, dB? dB: gen_0,
                              mkvecsmall5(degA, degB, dres, sX, vY)));
  init_modular_big(&S);
  H = gen_crt("ZX_ZXY_resultant_all", worker, &S, dB, bound, 0, NULL,
               nxV_chinese_center, FpX_center_i);
  setvarn(H, vX); (void)delete_var();
  return gerepilecopy(av, H);
}

static long
ZX_ZXY_rnfequation_lambda(GEN A, GEN B0, long lambda)
{
  pari_sp av = avma;
  long degA = degpol(A), degB, dres = degA*degpol(B0);
  long v = fetch_var_higher();
  long vX = varn(B0), vY = varn(A); /* assume vY has lower priority */
  long sX = evalvarn(vX);
  GEN dB, B, a, b, Hp;
  forprime_t S;

  B0 = Q_remove_denom(B0, &dB);
  if (!dB) B0 = leafcopy(B0);
  A = leafcopy(A);
  B = B0;
  setvarn(A,v);
INIT:
  if (lambda) B = RgX_translate(B0, monomial(stoi(lambda), 1, vY));
  B = swap_vars(B, vY); setvarn(B,v);
  /* B0(lambda v + x, v) */
  if (DEBUGLEVEL>4) err_printf("Trying lambda = %ld\n", lambda);

  degB = degpol(B);
  init_modular_big(&S);
  while (1)
  {
    ulong p = u_forprime_next(&S);
    ulong dp = dB ? umodiu(dB, p): 1;
    if (!dp) continue;
    a = ZX_to_Flx(A, p);
    b = ZXX_to_FlxX(B, p, v);
    Hp = ZX_ZXY_resultant_prime(a, b, dp, p, degA, degB, dres, sX);
    if (degpol(Hp) != dres) continue;
    if (dp != 1) Hp = Flx_Fl_mul(Hp, Fl_powu(Fl_inv(dp,p), degA, p), p);
    if (!Flx_is_squarefree(Hp, p)) { lambda = next_lambda(lambda); goto INIT; }
    if (DEBUGLEVEL>4) err_printf("Final lambda = %ld\n", lambda);
    (void)delete_var(); return gc_long(av,lambda);
  }
}

GEN
ZX_ZXY_rnfequation(GEN A, GEN B, long *lambda)
{
  if (lambda)
  {
    *lambda = ZX_ZXY_rnfequation_lambda(A, B, *lambda);
    if (*lambda) B = RgX_translate(B, monomial(stoi(*lambda), 1, varn(A)));
  }
  return ZX_ZXY_resultant(A,B);
}

static GEN
ZX_composedsum_slice(GEN A, GEN B, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1;
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZX_to_Flx(A, p), b = ZX_to_Flx(B, p);
    GEN Hp = Flx_composedsum(a, b, p);
    H = gerepileupto(av, Flx_to_ZX(Hp));
    *mod = utoipos(p); return H;
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  B = ZX_nv_mod_tree(B, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i);
    gel(H,i) = Flx_composedsum(a, b, p);
  }
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZX_composedsum_worker(GEN P, GEN A, GEN B)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = ZX_composedsum_slice(A, B, P, &gel(V,2));
  return V;
}

static GEN
ZX_composedsum_i(GEN A, GEN B, GEN lead)
{
  pari_sp av = avma;
  forprime_t S;
  ulong bound;
  GEN H, worker, mod;
  if (degpol(A) < degpol(B)) swap(A, B);
  if (!lead) lead  = mulii(leading_coeff(A),leading_coeff(B));
  bound = ZX_ZXY_ResBound_1(A, B);
  worker = snm_closure(is_entry("_ZX_composedsum_worker"), mkvec2(A,B));
  init_modular_big(&S);
  H = gen_crt("ZX_composedsum", worker, &S, lead, bound, 0, &mod,
              nxV_chinese_center, FpX_center);
  return gerepileupto(av, H);
}

static long
ZX_compositum_lambda(GEN A, GEN B, GEN lead, long lambda)
{
  pari_sp av = avma;
  forprime_t S;
  ulong p;
  init_modular_big(&S);
  p = u_forprime_next(&S);
  while (1)
  {
    GEN Hp, a;
    if (DEBUGLEVEL>4) err_printf("Trying lambda = %ld\n", lambda);
    if (lead && dvdiu(lead,p)) { p = u_forprime_next(&S); continue; }
    a = ZX_to_Flx(ZX_rescale(A, stoi(-lambda)), p);
    Hp = Flx_composedsum(a, ZX_to_Flx(B, p), p);
    if (!Flx_is_squarefree(Hp, p)) { lambda = next_lambda(lambda); continue; }
    if (DEBUGLEVEL>4) err_printf("Final lambda = %ld\n", lambda);
    return gc_long(av, lambda);
  }
}

GEN
ZX_compositum(GEN A, GEN B, long *lambda)
{
  GEN lead  = mulii(leading_coeff(A),leading_coeff(B));
  if (lambda)
  {
    *lambda = ZX_compositum_lambda(A, B, lead, *lambda);
    A = ZX_rescale(A, stoi(-*lambda));
  }
  return ZX_composedsum_i(A, B, lead);
}

GEN
ZX_composedsum(GEN A, GEN B)
{ return ZX_composedsum_i(A, B, NULL); }

static GEN
ZXQX_composedsum_slice(GEN A, GEN B, GEN C, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1, dC = degpol(C), v = varn(C);
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZXX_to_FlxX(A, p, v), b = ZXX_to_FlxX(B, p, v);
    GEN c = ZX_to_Flx(C, p);
    GEN Hp = FlxX_to_Flm(FlxqX_composedsum(a, b, c, p), dC);
    H = gerepileupto(av, Flm_to_ZM(Hp));
    *mod = utoipos(p); return H;
  }
  T = ZV_producttree(P);
  A = ZXX_nv_mod_tree(A, P, T, v);
  B = ZXX_nv_mod_tree(B, P, T, v);
  C = ZX_nv_mod_tree(C, P, T);
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i), c = gel(C,i);
    gel(H,i) = FlxX_to_Flm(FlxqX_composedsum(a, b, c, p), dC);
  }
  H = nmV_chinese_center_tree_seq(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZXQX_composedsum_worker(GEN P, GEN A, GEN B, GEN C)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = ZXQX_composedsum_slice(A, B, C, P, &gel(V,2));
  return V;
}

static GEN
ZXQX_composedsum(GEN A, GEN B, GEN T, ulong bound)
{
  pari_sp av = avma;
  forprime_t S;
  GEN H, worker, mod;
  GEN lead = mulii(Q_content(leading_coeff(A)), Q_content(leading_coeff(B)));
  worker = snm_closure(is_entry("_ZXQX_composedsum_worker")
                      , mkvec3(A,B,T));
  init_modular_big(&S);
  H = gen_crt("ZXQX_composedsum", worker, &S, lead, bound, 0, &mod,
              nmV_chinese_center, FpM_center);
  if (DEBUGLEVEL > 4)
    err_printf("nfcompositum: a priori bound: %lu, a posteriori: %lu\n",
               bound, expi(gsupnorm(H, DEFAULTPREC)));
  return gerepilecopy(av, RgM_to_RgXX(H, varn(A), varn(T)));
}

static long
ZXQX_composedsum_bound(GEN nf, GEN A, GEN B)
{ return ZXQX_resultant_bound_i(nf, A, B, &RgX_RgXY_ResBound_1); }

GEN
nf_direct_compositum(GEN nf, GEN A, GEN B)
{
  ulong bnd = ZXQX_composedsum_bound(nf, A, B);
  return ZXQX_composedsum(A, B, nf_get_pol(nf), bnd);
}

/************************************************************************
 *                                                                      *
 *                   IRREDUCIBLE POLYNOMIAL / Fp                        *
 *                                                                      *
 ************************************************************************/

/* irreducible (unitary) polynomial of degree n over Fp */
GEN
ffinit_rand(GEN p,long n)
{
  for(;;) {
    pari_sp av = avma;
    GEN pol = ZX_add(pol_xn(n, 0), random_FpX(n-1,0, p));
    if (FpX_is_irred(pol, p)) return pol;
    set_avma(av);
  }
}

/* return an extension of degree 2^l of F_2, assume l > 0
 * Not stack clean. */
static GEN
ffinit_Artin_Schreier_2(long l)
{
  GEN Q, T, S;
  long i, v;

  if (l == 1) return mkvecsmall4(0,1,1,1); /*x^2 + x + 1*/
  v = fetch_var_higher();
  S = mkvecsmall5(0, 0, 0, 1, 1); /* y(y^2 + y) */
  Q = mkpoln(3, pol1_Flx(0), pol1_Flx(0), S); /* x^2 + x + y(y^2+y) */
  setvarn(Q, v);

  /* x^4+x+1, irred over F_2, minimal polynomial of a root of Q */
  T = mkvecsmalln(6,evalvarn(v),1UL,1UL,0UL,0UL,1UL);
  /* Q = x^2 + x + a(y) irred. over K = F2[y] / (T(y))
   * ==> x^2 + x + a(y) b irred. over K for any root b of Q
   * ==> x^2 + x + (b^2+b)b */
  for (i=2; i<l; i++) T = Flx_FlxY_resultant(T, Q, 2); /* minpoly(b) / F2*/
  (void)delete_var(); T[1] = 0; return T;
}

/* return an extension of degree p^l of F_p, assume l > 0
 * Not stack clean. */
GEN
ffinit_Artin_Schreier(ulong p, long l)
{
  long i, v;
  GEN Q, R, S, T, xp;
  if (p==2) return ffinit_Artin_Schreier_2(l);
  xp = polxn_Flx(p,0); /* x^p */
  T = Flx_sub(xp, mkvecsmall3(0,1,1),p); /* x^p - x - 1 */
  if (l == 1) return T;

  v = evalvarn(fetch_var_higher());
  xp[1] = v;
  R = Flx_sub(polxn_Flx(2*p-1,0), polxn_Flx(p,0),p);
  S = Flx_sub(xp, polx_Flx(0), p);
  Q = FlxX_Flx_sub(Flx_to_FlxX(S, v), R, p); /* x^p - x - (y^(2p-1)-y^p) */
  for (i = 2; i <= l; ++i) T = Flx_FlxY_resultant(T, Q, p);
  (void)delete_var(); T[1] = 0; return T;
}

static long
flinit_check(ulong p, long n, long l)
{
  ulong q;
  if (!uisprime(n)) return 0;
  q = p % n; if (!q) return 0;
  return ugcd((n-1)/Fl_order(q, n-1, n), l) == 1;
}

static GEN
flinit(ulong p, long l)
{
  ulong n = 1+l;
  while (!flinit_check(p,n,l)) n += l;
  if (DEBUGLEVEL>=4) err_printf("FFInit: using polsubcyclo(%ld, %ld)\n",n,l);
  return ZX_to_Flx(polsubcyclo(n,l,0), p);
}

static GEN
ffinit_fact_Flx(ulong p, long n)
{
  GEN P, F = factoru_pow(n), Fp = gel(F,1), Fe = gel(F,2), Fm = gel(F,3);
  long i, l = lg(Fm);
  P = cgetg(l, t_VEC);
  for (i = 1; i < l; ++i)
    gel(P,i) = p==uel(Fp,i) ?
                 ffinit_Artin_Schreier(uel(Fp,i), Fe[i])
               : flinit(p, uel(Fm,i));
  return FlxV_composedsum(P, p);
}

static GEN
init_Flxq_i(ulong p, long n, long sv)
{
  GEN P;
  if (n == 1) return polx_Flx(sv);
  if (flinit_check(p, n+1, n))
  {
    P = const_vecsmall(n+2,1);
    P[1] = sv; return P;
  }
  P = ffinit_fact_Flx(p,n);
  P[1] = sv; return P;
}

GEN
init_Flxq(ulong p, long n, long v)
{
  pari_sp av = avma;
  return gerepileupto(av, init_Flxq_i(p, n, v));
}

/* check if polsubcyclo(n,l,0) is irreducible modulo p */
static long
fpinit_check(GEN p, long n, long l)
{
  ulong q;
  if (!uisprime(n)) return 0;
  q = umodiu(p,n); if (!q) return 0;
  return ugcd((n-1)/Fl_order(q, n-1, n), l) == 1;
}

/* let k=2 if p%4==1, and k=4 else and assume k*p does not divide l.
 * Return an irreducible polynomial of degree l over F_p.
 * Variant of Adleman and Lenstra "Finding irreducible polynomials over
 * finite fields", ACM, 1986 (5) 350--355.
 * Not stack clean */
static GEN
fpinit(GEN p, long l)
{
  ulong n = 1+l;
  while (!fpinit_check(p,n,l)) n += l;
  if (DEBUGLEVEL>=4) err_printf("FFInit: using polsubcyclo(%ld, %ld)\n",n,l);
  return FpX_red(polsubcyclo(n,l,0),p);
}

static GEN
ffinit_fact(GEN p, long n)
{
  GEN P, F = factoru_pow(n), Fp = gel(F,1), Fe = gel(F,2), Fm = gel(F,3);
  long i, l = lg(Fm);
  P = cgetg(l, t_VEC);
  for (i = 1; i < l; ++i)
    gel(P,i) = absequaliu(p, Fp[i]) ?
                 Flx_to_ZX(ffinit_Artin_Schreier(Fp[i], Fe[i]))
               : fpinit(p, Fm[i]);
  return FpXV_composedsum(P, p);
}

static GEN
init_Fq_i(GEN p, long n, long v)
{
  GEN P;
  if (n <= 0) pari_err_DOMAIN("ffinit", "degree", "<=", gen_0, stoi(n));
  if (typ(p) != t_INT) pari_err_TYPE("ffinit",p);
  if (cmpiu(p, 2) < 0) pari_err_PRIME("ffinit",p);
  if (v < 0) v = 0;
  if (n == 1) return pol_x(v);
  if (lgefint(p) == 3)
    return Flx_to_ZX(init_Flxq_i(p[2], n, evalvarn(v)));
  if (fpinit_check(p, n+1, n)) return polcyclo(n+1, v);
  P = ffinit_fact(p,n);
  setvarn(P, v); return P;
}
GEN
init_Fq(GEN p, long n, long v)
{
  pari_sp av = avma;
  return gerepileupto(av, init_Fq_i(p, n, v));
}
GEN
ffinit(GEN p, long n, long v)
{
  pari_sp av = avma;
  return gerepileupto(av, FpX_to_mod(init_Fq_i(p, n, v), p));
}

GEN
ffnbirred(GEN p, long n)
{
  pari_sp av = avma;
  GEN s = powiu(p,n), F = factoru(n), D = divisorsu_moebius(gel(F, 1));
  long j, l = lg(D);
  for (j = 2; j < l; j++) /* skip d = 1 */
  {
    long md = D[j]; /* mu(d) * d, d squarefree */
    GEN pd = powiu(p, n / labs(md)); /* p^{n/d} */
    s = md > 0? addii(s, pd): subii(s,pd);
  }
  return gerepileuptoint(av, diviuexact(s, n));
}

GEN
ffsumnbirred(GEN p, long n)
{
  pari_sp av = avma, av2;
  GEN q, t = p, v = vecfactoru_i(1, n);
  long i;
  q = cgetg(n+1,t_VEC); gel(q,1) = p;
  for (i=2; i<=n; i++) gel(q,i) = mulii(gel(q,i-1), p);
  av2 = avma;
  for (i=2; i<=n; i++)
  {
    GEN s = gel(q,i), F = gel(v,i), D = divisorsu_moebius(gel(F,1));
    long j, l = lg(D);
    for (j = 2; j < l; j++) /* skip 1 */
    {
      long md = D[j];
      GEN pd = gel(q, i / labs(md)); /* p^{i/d} */
      s = md > 0? addii(s, pd): subii(s, pd);
    }
    t = gerepileuptoint(av2, addii(t, diviuexact(s, i)));
  }
  return gerepileuptoint(av, t);
}

GEN
ffnbirred0(GEN p, long n, long flag)
{
  if (typ(p) != t_INT) pari_err_TYPE("ffnbirred", p);
  if (n <= 0) pari_err_DOMAIN("ffnbirred", "degree", "<=", gen_0, stoi(n));
  switch(flag)
  {
    case 0: return ffnbirred(p, n);
    case 1: return ffsumnbirred(p, n);
  }
  pari_err_FLAG("ffnbirred");
  return NULL; /* LCOV_EXCL_LINE */
}

static void
checkmap(GEN m, const char *s)
{
  if (typ(m)!=t_VEC || lg(m)!=3 || typ(gel(m,1))!=t_FFELT)
    pari_err_TYPE(s,m);
}

GEN
ffembed(GEN a, GEN b)
{
  pari_sp av = avma;
  GEN p, Ta, Tb, g, r = NULL;
  if (typ(a)!=t_FFELT) pari_err_TYPE("ffembed",a);
  if (typ(b)!=t_FFELT) pari_err_TYPE("ffembed",b);
  p = FF_p_i(a); g = FF_gen(a);
  if (!equalii(p, FF_p_i(b))) pari_err_MODULUS("ffembed",a,b);
  Ta = FF_mod(a);
  Tb = FF_mod(b);
  if (degpol(Tb)%degpol(Ta)!=0)
    pari_err_DOMAIN("ffembed",GENtostr_raw(a),"is not a subfield of",b,a);
  r = gel(FFX_roots(Ta, b), 1);
  return gerepilecopy(av, mkvec2(g,r));
}

GEN
ffextend(GEN a, GEN P, long v)
{
  pari_sp av = avma;
  long n;
  GEN p, T, R, g, m;
  if (typ(a)!=t_FFELT) pari_err_TYPE("ffextend",a);
  T = a; p = FF_p_i(a);
  if (typ(P)!=t_POL || !RgX_is_FpXQX(P,&T,&p)) pari_err_TYPE("ffextend", P);
  if (!FF_samefield(a, T)) pari_err_MODULUS("ffextend",a,T);
  if (v < 0) v = varn(P);
  n = FF_f(T) * degpol(P); R = ffinit(p, n, v); g = ffgen(R, v);
  m = ffembed(a, g);
  R = FFX_roots(ffmap(m, P),g);
  return gerepilecopy(av, mkvec2(gel(R,1), m));
}

GEN
fffrobenius(GEN a, long n)
{
  if (typ(a)!=t_FFELT) pari_err_TYPE("fffrobenius",a);
  retmkvec2(FF_gen(a), FF_Frobenius(a, n));
}

GEN
ffinvmap(GEN m)
{
  pari_sp av = avma;
  long i, l;
  GEN T, F, a, g, r, f = NULL;
  checkmap(m, "ffinvmap");
  a = gel(m,1); r = gel(m,2);
  if (typ(r) != t_FFELT)
   pari_err_TYPE("ffinvmap", m);
  g = FF_gen(a);
  T = FF_mod(r);
  F = gel(FFX_factor(T, a), 1);
  l = lg(F);
  for(i=1; i<l; i++)
  {
    GEN s = FFX_rem(FF_to_FpXQ_i(r), gel(F, i), a);
    if (degpol(s)==0 && gequal(constant_coeff(s),g)) { f = gel(F, i); break; }
  }
  if (f==NULL) pari_err_TYPE("ffinvmap", m);
  if (degpol(f)==1) f = FF_neg_i(gel(f,2));
  return gerepilecopy(av, mkvec2(FF_gen(r),f));
}

static GEN
ffpartmapimage(const char *s, GEN r)
{
   GEN a = NULL, p = NULL;
   if (typ(r)==t_POL && degpol(r) >= 1
      && RgX_is_FpXQX(r,&a,&p) && a && typ(a)==t_FFELT) return a;
   pari_err_TYPE(s, r);
   return NULL; /* LCOV_EXCL_LINE */
}

static GEN
ffeltmap_i(GEN m, GEN x)
{
   GEN r = gel(m,2);
   if (!FF_samefield(x, gel(m,1)))
     pari_err_DOMAIN("ffmap","m","domain does not contain", x, r);
   if (typ(r)==t_FFELT)
     return FF_map(r, x);
   else
     return FFX_preimage(x, r, ffpartmapimage("ffmap", r));
}

static GEN
ffmap_i(GEN m, GEN x)
{
  GEN y;
  long i, lx, tx = typ(x);
  switch(tx)
  {
    case t_FFELT:
      return ffeltmap_i(m, x);
    case t_POL: case t_RFRAC: case t_SER:
    case t_VEC: case t_COL: case t_MAT:
      y = cgetg_copy(x, &lx);
      for (i=1; i<lontyp[tx]; i++) y[i] = x[1];
      for (i=lontyp[tx]; i<lx; i++)
      {
        GEN yi = ffmap_i(m, gel(x,i));
        if (!yi) return NULL;
        gel(y,i) = yi;
      }
      return y;
  }
  return gcopy(x);
}

GEN
ffmap(GEN m, GEN x)
{
  pari_sp ltop = avma;
  GEN y;
  checkmap(m, "ffmap");
  y = ffmap_i(m, x);
  if (y) return y;
  set_avma(ltop); return cgetg(1,t_VEC);
}

static GEN
ffeltmaprel_i(GEN m, GEN x)
{
   GEN g = gel(m,1), r = gel(m,2);
   if (!FF_samefield(x, g))
     pari_err_DOMAIN("ffmap","m","domain does not contain", x, r);
   if (typ(r)==t_FFELT)
     retmkpolmod(FF_map(r, x), pol_x(FF_var(g)));
   else
     retmkpolmod(FFX_preimagerel(x, r, ffpartmapimage("ffmap", r)), gcopy(r));
}

static GEN
ffmaprel_i(GEN m, GEN x)
{
  GEN y;
  long i, lx, tx = typ(x);
  switch(tx)
  {
    case t_FFELT:
      return ffeltmaprel_i(m, x);
    case t_POL: case t_RFRAC: case t_SER:
    case t_VEC: case t_COL: case t_MAT:
      y = cgetg_copy(x, &lx);
      for (i=1; i<lontyp[tx]; i++) y[i] = x[1];
      for (i=lontyp[tx]; i<lx; i++)
        gel(y,i) = ffmaprel_i(m, gel(x,i));
      return y;
  }
  return gcopy(x);
}

GEN
ffmaprel(GEN m, GEN x)
{
  checkmap(m, "ffmaprel");
  return ffmaprel_i(m, x);
}

static void
err_compo(GEN m, GEN n)
{ pari_err_DOMAIN("ffcompomap","m","domain does not contain codomain of",n,m); }

GEN
ffcompomap(GEN m, GEN n)
{
  pari_sp av = avma;
  GEN g = gel(n,1), r, m2, n2;
  checkmap(m, "ffcompomap");
  checkmap(n, "ffcompomap");
  m2 = gel(m,2); n2 = gel(n,2);
  switch((typ(m2)==t_POL)|((typ(n2)==t_POL)<<1))
  {
    case 0:
      if (!FF_samefield(gel(m,1),n2)) err_compo(m,n);
      r = FF_map(gel(m,2), n2);
      break;
    case 2:
      r = ffmap_i(m, n2);
      if (lg(r) == 1) err_compo(m,n);
      break;
    case 1:
      r = ffeltmap_i(m, n2);
      if (!r)
      {
        GEN a, A, R, M;
        long dm, dn;
        a = ffpartmapimage("ffcompomap",m2);
        A = FF_to_FpXQ_i(FF_neg(n2));
        setvarn(A, 1);
        R = deg1pol(gen_1, A, 0);
        setvarn(R, 0);
        M = gcopy(m2);
        setvarn(M, 1);
        r = polresultant0(R, M, 1, 0);
        dm = FF_f(gel(m,1)); dn = FF_f(gel(n,1));
        if (dm % dn || !FFX_ispower(r, dm/dn, a, &r)) err_compo(m,n);
        setvarn(r, varn(FF_mod(g)));
      }
      break;
    case 3:
    {
      GEN M, R, T, p, a;
      a = ffpartmapimage("ffcompomap",n2);
      if (!FF_samefield(a, gel(m,1))) err_compo(m,n);
      p = FF_p_i(gel(n,1));
      T = FF_mod(gel(n,1));
      setvarn(T, 1);
      R = RgX_to_FpXQX(n2,T,p);
      setvarn(R, 0);
      M = gcopy(m2);
      setvarn(M, 1);
      r = polresultant0(R, M, 1, 0);
      setvarn(r, varn(n2));
    }
  }
  return gerepilecopy(av, mkvec2(g,r));
}
