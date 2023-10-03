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

/**************************************************************/
/*                                                            */
/*                        NUMBER FIELDS                       */
/*                                                            */
/**************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_nf

int new_galois_format = 0;

/* v a t_VEC, lg(v) = 13, sanity check for true rnf */
static int
v13checkrnf(GEN v)
{ return typ(gel(v,6)) == t_VEC; }
static int
rawcheckbnf(GEN v) { return typ(v)==t_VEC && lg(v)==11; }
static int
rawchecknf(GEN v) { return typ(v)==t_VEC && lg(v)==10; }
/* v a t_VEC, lg(v) = 11, sanity check for true bnf */
static int
v11checkbnf(GEN v) { return rawchecknf(bnf_get_nf(v)); }
/* v a t_VEC, lg(v) = 13, sanity check for true nf and true bnf */
static int
v13checkgchar(GEN v) { return rawchecknf(gchar_get_nf(v)) && rawcheckbnf(gchar_get_bnf(v)); }
/* v a t_VEC, lg(v) = 10, sanity check for true nf */
static int
v10checknf(GEN v) { return typ(gel(v,1))==t_POL; }
/* v a t_VEC, lg(v) = 9, sanity check for true gal */
static int
v9checkgal(GEN v)
{ GEN x = gel(v,2); return typ(x) == t_VEC && lg(x) == 4; }

int
checkrnf_i(GEN rnf)
{ return (typ(rnf)==t_VEC && lg(rnf)==13 && v13checkrnf(rnf)); }

void
checkrnf(GEN rnf)
{ if (!checkrnf_i(rnf)) pari_err_TYPE("checkrnf",rnf); }

GEN
checkbnf_i(GEN X)
{
  if (typ(X) == t_VEC)
    switch (lg(X))
    {
      case 11:
        if (typ(gel(X,6)) != t_INT) return NULL; /* pre-2.2.4 format */
        if (lg(gel(X,10)) != 4) return NULL; /* pre-2.8.1 format */
        return X;
      case 7: return checkbnf_i(bnr_get_bnf(X));
    }
  return NULL;
}

GEN
checknf_i(GEN X)
{
  if (typ(X)==t_VEC)
    switch(lg(X))
    {
      case 10: return X;
      case 11: return checknf_i(bnf_get_nf(X));
      case 7:  return checknf_i(bnr_get_bnf(X));
      case 3: if (typ(gel(X,2)) == t_POLMOD) return checknf_i(gel(X,1));
    }
  return NULL;
}

GEN
checkbnf(GEN x)
{
  GEN bnf = checkbnf_i(x);
  if (!bnf) pari_err_TYPE("checkbnf [please apply bnfinit()]",x);
  return bnf;
}

GEN
checknf(GEN x)
{
  GEN nf = checknf_i(x);
  if (!nf) pari_err_TYPE("checknf [please apply nfinit()]",x);
  return nf;
}

GEN
checkbnr_i(GEN bnr)
{
  if (typ(bnr)!=t_VEC || lg(bnr)!=7 || !checkbnf_i(bnr_get_bnf(bnr)))
    return NULL;
  return bnr;
}
void
checkbnr(GEN bnr)
{
  if (!checkbnr_i(bnr))
    pari_err_TYPE("checkbnr [please apply bnrinit()]",bnr);
}

void
checksqmat(GEN x, long N)
{
  if (typ(x)!=t_MAT) pari_err_TYPE("checksqmat",x);
  if (lg(x) == 1 || lgcols(x) != N+1) pari_err_DIM("checksqmat");
}

GEN
checkbid_i(GEN bid)
{
  GEN f;
  if (typ(bid)!=t_VEC || lg(bid)!=6 || typ(bid_get_U(bid)) != t_VEC)
    return NULL;
  f = bid_get_mod(bid);
  if (typ(f)!=t_VEC || lg(f)!=3) return NULL;
  return bid;
}
void
checkbid(GEN bid)
{
  if (!checkbid_i(bid)) pari_err_TYPE("checkbid",bid);
}
void
checkabgrp(GEN v)
{
  if (typ(v) == t_VEC) switch(lg(v))
  {
    case 4: if (typ(gel(v,3)) != t_VEC) break;
    case 3: if (typ(gel(v,2)) != t_VEC) break;
            if (typ(gel(v,1)) != t_INT) break;
            return;/*OK*/
    default: break;
  }
  pari_err_TYPE("checkabgrp",v);
}

GEN
checknfelt_mod(GEN nf, GEN x, const char *s)
{
  GEN T = gel(x,1), a = gel(x,2), Tnf = nf_get_pol(nf);
  if (!RgX_equal_var(T, Tnf)) pari_err_MODULUS(s, T, Tnf);
  return a;
}

int
check_ZKmodule_i(GEN M)
{
  return (typ(M) ==t_VEC && lg(M) >= 3
          && typ(gel(M,1)) == t_MAT
          && typ(gel(M,2)) == t_VEC
          && lgcols(M) == lg(gel(M,2)));
}
void
check_ZKmodule(GEN M, const char *s)
{ if (!check_ZKmodule_i(M)) pari_err_TYPE(s, M); }

static long
typv6(GEN x)
{
  if (typ(gel(x,1)) == t_VEC && lg(gel(x,3)) == 3)
  {
    GEN t = gel(x,3);
    if (typ(t) != t_VEC) return typ_NULL;
    t = gel(x,5);
    switch(typ(gel(x,5)))
    {
      case t_VEC: return typ_BID;
      case t_MAT: return typ_BIDZ;
      default: return typ_NULL;
    }
  }
  if (typ(gel(x,2)) == t_COL && typ(gel(x,3)) == t_INT) return typ_PRID;
  return typ_NULL;
}

GEN
get_bnf(GEN x, long *t)
{
  switch(typ(x))
  {
    case t_POL: *t = typ_POL;  return NULL;
    case t_QUAD: *t = typ_Q  ; return NULL;
    case t_VEC:
      switch(lg(x))
      {
        case 5: if (typ(gel(x,1)) != t_INT) break;
                *t = typ_QUA; return NULL;
        case 6: *t = typv6(x); return NULL;
        case 7:  *t = typ_BNR;
          x = bnr_get_bnf(x);
          if (!rawcheckbnf(x)) break;
          return x;
        case 9:
          if (!v9checkgal(x)) break;
          *t = typ_GAL; return NULL;
        case 10:
          if (!v10checknf(x)) break;
          *t = typ_NF; return NULL;
        case 11:
          if (!v11checkbnf(x)) break;
          *t = typ_BNF; return x;
        case 13:
          if (v13checkgchar(x)) { *t = typ_GCHAR; return gchar_get_bnf(x); }
          if (!v13checkrnf(x)) break;
          *t = typ_RNF; return NULL;
        case 17: *t = typ_ELL; return NULL;
      }
      break;
    case t_COL:
      if (get_prid(x)) { *t = typ_MODPR; return NULL; }
      break;
  }
  *t = typ_NULL; return NULL;
}

GEN
get_nf(GEN x, long *t)
{
  switch(typ(x))
  {
    case t_POL : *t = typ_POL; return NULL;
    case t_QUAD: *t = typ_Q  ; return NULL;
    case t_VEC:
      switch(lg(x))
      {
        case 3:
          if (typ(gel(x,2)) != t_POLMOD) break;
          return get_nf(gel(x,1),t);
        case 5:
          if (typ(gel(x,1)) != t_INT) break;
          *t = typ_QUA; return NULL;
        case 6: *t = typv6(x); return NULL;
        case 7:
          x = bnr_get_bnf(x);
          if (!rawcheckbnf(x) || !rawchecknf(x = bnf_get_nf(x))) break;
          *t = typ_BNR; return x;
        case 9:
          if (!v9checkgal(x)) break;
          *t = typ_GAL; return NULL;
        case 10:
          if (!v10checknf(x)) break;
          *t = typ_NF; return x;
        case 11:
          if (!rawchecknf(x = bnf_get_nf(x))) break;
          *t = typ_BNF; return x;
        case 13:
          if (v13checkgchar(x)) { *t = typ_GCHAR; return gchar_get_nf(x); }
          if (!v13checkrnf(x)) break;
          *t = typ_RNF; return NULL;
        case 17: *t = typ_ELL; return NULL;
      }
      break;
    case t_QFB: *t = typ_QFB; return NULL;
    case t_COL:
      if (get_prid(x)) { *t = typ_MODPR; return NULL; }
      break;
  }
  *t = typ_NULL; return NULL;
}

long
nftyp(GEN x)
{
  switch(typ(x))
  {
    case t_POL : return typ_POL;
    case t_QUAD: return typ_Q;
    case t_VEC:
      switch(lg(x))
      {
        case 13:
          if (v13checkgchar(x)) return typ_GCHAR;
          if (!v13checkrnf(x)) break;
          return typ_RNF;
        case 10:
          if (!v10checknf(x)) break;
          return typ_NF;
        case 11:
          if (!v11checkbnf(x)) break;
          return typ_BNF;
        case 7:
          x = bnr_get_bnf(x);
          if (!rawcheckbnf(x) || !v11checkbnf(x)) break;
          return typ_BNR;
        case 6:
          return typv6(x);
        case 9:
          if (!v9checkgal(x)) break;
          return typ_GAL;
        case 17: return typ_ELL;
      }
  }
  return typ_NULL;
}

/*************************************************************************/
/**                                                                     **/
/**                           GALOIS GROUP                              **/
/**                                                                     **/
/*************************************************************************/

GEN
tschirnhaus(GEN x)
{
  pari_sp av = avma, av2;
  long a, v = varn(x);
  GEN u, y = cgetg(5,t_POL);

  if (typ(x)!=t_POL) pari_err_TYPE("tschirnhaus",x);
  if (lg(x) < 4) pari_err_CONSTPOL("tschirnhaus");
  if (v) { u = leafcopy(x); setvarn(u,0); x=u; }
  y[1] = evalsigne(1)|evalvarn(0);
  do
  {
    a = random_bits(2); if (a==0) a  = 1; gel(y,4) = stoi(a);
    a = random_bits(3); if (a>=4) a -= 8; gel(y,3) = stoi(a);
    a = random_bits(3); if (a>=4) a -= 8; gel(y,2) = stoi(a);
    u = RgXQ_charpoly(y,x,v); av2 = avma;
  }
  while (degpol(RgX_gcd(u,RgX_deriv(u)))); /* while u not separable */
  if (DEBUGLEVEL>1)
    err_printf("Tschirnhaus transform. New pol: %Ps",u);
  set_avma(av2); return gerepileupto(av,u);
}

/* Assume pol in Z[X], monic of degree n. Find L in Z such that
 * POL = L^(-n) pol(L x) is monic in Z[X]. Return POL and set *ptk = L.
 * No GC. */
GEN
ZX_Z_normalize(GEN pol, GEN *ptk)
{
  long i,j, sk, n = degpol(pol); /* > 0 */
  GEN k, fa, P, E, a, POL;

  if (ptk) *ptk = gen_1;
  if (!n) return pol;
  a = pol + 2; k = gel(a,n-1); /* a[i] = coeff of degree i */
  for (i = n-2; i >= 0; i--)
  {
    k = gcdii(k, gel(a,i));
    if (is_pm1(k)) return pol;
  }
  sk = signe(k);
  if (!sk) return pol; /* monomial! */
  fa = absZ_factor_limit(k, 0); k = gen_1;
  P = gel(fa,1);
  E = gel(fa,2);
  POL = leafcopy(pol); a = POL+2;
  for (i = lg(P)-1; i > 0; i--)
  {
    GEN p = gel(P,i), pv, pvj;
    long vmin = itos(gel(E,i));
    /* find v_p(k) = min floor( v_p(a[i]) / (n-i)) */
    for (j=n-1; j>=0; j--)
    {
      long v;
      if (!signe(gel(a,j))) continue;
      v = Z_pval(gel(a,j), p) / (n - j);
      if (v < vmin) vmin = v;
    }
    if (!vmin) continue;
    pvj = pv = powiu(p,vmin); k = mulii(k, pv);
    /* a[j] /= p^(v*(n-j)) */
    for (j=n-1; j>=0; j--)
    {
      if (j < n-1) pvj = mulii(pvj, pv);
      gel(a,j) = diviiexact(gel(a,j), pvj);
    }
  }
  if (ptk) *ptk = k;
  return POL;
}

/* Assume pol != 0 in Z[X]. Find C in Q, L in Z such that POL = C pol(x/L) monic
 * in Z[X]. Return POL and set *pL = L. Wasteful (but correct) if pol is not
 * primitive: better if caller used Q_primpart already. No GC. */
GEN
ZX_primitive_to_monic(GEN pol, GEN *pL)
{
  long i,j, n = degpol(pol);
  GEN lc = leading_coeff(pol), L, fa, P, E, a, POL;

  if (is_pm1(lc))
  {
    if (pL) *pL = gen_1;
    return signe(lc) < 0? ZX_neg(pol): pol;
  }
  if (signe(lc) < 0)
    POL = ZX_neg(pol);
  else
    POL = leafcopy(pol);
  a = POL+2; lc = gel(a,n);
  fa = absZ_factor_limit(lc,0); L = gen_1;
  P = gel(fa,1);
  E = gel(fa,2);
  for (i = lg(P)-1; i > 0; i--)
  {
    GEN p = gel(P,i), pk, pku;
    long v, j0, e = itos(gel(E,i)), k = e/n, d = k*n - e;

    if (d < 0) { k++; d += n; }
    /* k = ceil(e[i] / n); find d, k such that  p^d pol(x / p^k) monic */
    for (j=n-1; j>0; j--)
    {
      if (!signe(gel(a,j))) continue;
      v = Z_pval(gel(a,j), p);
      while (v + d < k * j) { k++; d += n; }
    }
    pk = powiu(p,k); j0 = d/k;
    L = mulii(L, pk);

    pku = powiu(p,d - k*j0);
    /* a[j] *= p^(d - kj) */
    for (j=j0; j>=0; j--)
    {
      if (j < j0) pku = mulii(pku, pk);
      gel(a,j) = mulii(gel(a,j), pku);
    }
    j0++;
    pku = powiu(p,k*j0 - d);
    /* a[j] /= p^(kj - d) */
    for (j=j0; j<=n; j++)
    {
      if (j > j0) pku = mulii(pku, pk);
      gel(a,j) = diviiexact(gel(a,j), pku);
    }
  }
  if (pL) *pL = L;
  return POL;
}
/* Assume pol != 0 in Z[X]. Find C,L in Q such that POL = C pol(x/L)
 * monic in Z[X]. Return POL and set *pL = L.
 * Wasteful (but correct) if pol is not primitive: better if caller used
 * Q_primpart already. No GC. */
GEN
ZX_Q_normalize(GEN pol, GEN *pL)
{
  GEN lc, POL = ZX_primitive_to_monic(pol, &lc);
  POL = ZX_Z_normalize(POL, pL);
  if (pL) *pL = gdiv(lc, *pL);
  return POL;
}

GEN
ZX_Q_mul(GEN A, GEN z)
{
  pari_sp av = avma;
  long i, l = lg(A);
  GEN d, n, Ad, B, u;
  if (typ(z)==t_INT) return ZX_Z_mul(A,z);
  n = gel(z, 1); d = gel(z, 2);
  Ad = RgX_to_RgC(FpX_red(A, d), l-2);
  u = gcdii(d, FpV_factorback(Ad, NULL, d));
  B = cgetg(l, t_POL);
  B[1] = A[1];
  if (equali1(u))
  {
    for(i=2; i<l; i++)
      gel(B, i) = mkfrac(mulii(n, gel(A,i)), d);
  } else
  {
    for(i=2; i<l; i++)
    {
      GEN di = gcdii(gel(Ad, i-1), u);
      GEN ni = mulii(n, diviiexact(gel(A,i), di));
      if (equalii(d, di))
        gel(B, i) = ni;
      else
        gel(B, i) = mkfrac(ni, diviiexact(d, di));
    }
  }
  return gerepilecopy(av, B);
}

/* T != 0 in Z[x], returns a monic polynomial U in Z[x] generating the
 * same field: there exist C in Q, L in Z such that U(x) = C T(x/L).
 * Set *L = NULL if L = 1, and to L otherwise. No garbage collecting. */
GEN
ZX_to_monic(GEN T, GEN *L)
{
  GEN lc = leading_coeff(T);
  if (is_pm1(lc)) { *L = gen_1; return signe(lc) > 0? T: ZX_neg(T); }
  return ZX_primitive_to_monic(Q_primpart(T), L);
}

GEN
poltomonic(GEN T, GEN *L)
{
  pari_sp av = avma;
  if (typ(T) != t_POL || !RgX_is_QX(T)) pari_err_TYPE("poltomonic", T);
  if (degpol(T) < 0) pari_err_ROOTS0("poltomonic");
  T = ZX_Q_normalize(Q_primpart(T), L); return gc_all(av, L? 2: 1, &T, L);
}

GEN
ZXX_Q_mul(GEN A, GEN z)
{
  long i, l;
  GEN B;
  if (typ(z)==t_INT) return ZXX_Z_mul(A,z);
  B = cgetg_copy(A, &l);
  B[1] = A[1];
  for (i=2; i<l; i++)
  {
    GEN Ai = gel(A,i);
    gel(B,i) = typ(Ai)==t_POL ? ZX_Q_mul(Ai, z): gmul(Ai, z);
  }
  return B;
}

/* Evaluate pol in s using nfelt arithmetic and Horner rule */
GEN
nfpoleval(GEN nf, GEN pol, GEN s)
{
  pari_sp av=avma;
  long i=lg(pol)-1;
  GEN res;
  if (i==1) return gen_0;
  res = nf_to_scalar_or_basis(nf, gel(pol,i));
  for (i-- ; i>=2; i--)
    res = nfadd(nf, nfmul(nf, s, res), gel(pol,i));
  return gerepileupto(av, res);
}

static GEN
QX_table_nfpoleval(GEN nf, GEN pol, GEN m)
{
  pari_sp av = avma;
  long i = lg(pol)-1;
  GEN res, den;
  if (i==1) return gen_0;
  pol = Q_remove_denom(pol, &den);
  res = scalarcol_shallow(gel(pol,i), nf_get_degree(nf));
  for (i-- ; i>=2; i--)
    res = ZC_Z_add(ZM_ZC_mul(m, res), gel(pol,i));
  if (den) res = RgC_Rg_div(res, den);
  return gerepileupto(av, res);
}

GEN
FpX_FpC_nfpoleval(GEN nf, GEN pol, GEN a, GEN p)
{
  pari_sp av=avma;
  long i=lg(pol)-1, n=nf_get_degree(nf);
  GEN res, Ma;
  if (i==1) return zerocol(n);
  Ma = FpM_red(zk_multable(nf, a), p);
  res = scalarcol(gel(pol,i),n);
  for (i-- ; i>=2; i--)
  {
    res = FpM_FpC_mul(Ma, res, p);
    gel(res,1) = Fp_add(gel(res,1), gel(pol,i), p);
  }
  return gerepileupto(av, res);
}

/* compute s(x), not stack clean */
static GEN
ZC_galoisapply(GEN nf, GEN s, GEN x)
{
  x = nf_to_scalar_or_alg(nf, x);
  if (typ(x) != t_POL) return scalarcol(x, nf_get_degree(nf));
  return QX_table_nfpoleval(nf, x, zk_multable(nf, s));
}

/* true nf; S = automorphism in basis form, return an FpC = S(z) mod p */
GEN
zk_galoisapplymod(GEN nf, GEN z, GEN S, GEN p)
{
  GEN den, pe, pe1, denpe, R;

  z = nf_to_scalar_or_alg(nf, z);
  if (typ(z) != t_POL) return z;
  if (gequalX(z)) return FpC_red(S, p); /* common, e.g. modpr_genFq */
  z = Q_remove_denom(z,&den);
  denpe = pe = NULL;
  pe1 = p;
  if (den)
  {
    ulong e = Z_pvalrem(den, p, &den);
    if (e) { pe = powiu(p, e); pe1 = mulii(pe, p); }
    denpe = Zp_inv(den, p, e+1);
  }
  R = FpX_FpC_nfpoleval(nf, FpX_red(z, pe1), FpC_red(S, pe1), pe1);
  if (denpe) R = FpC_Fp_mul(R, denpe, pe1);
  if (pe) R = gdivexact(R, pe);
  return R;
}

/* true nf */
static GEN
pr_make(GEN nf, GEN p, GEN u, GEN e, GEN f)
{
  GEN t = FpM_deplin(zk_multable(nf, u), p);
  t = zk_scalar_or_multable(nf, t);
  return mkvec5(p, u, e, f, t);
}
static GEN
pr_galoisapply(GEN nf, GEN pr, GEN aut)
{
  GEN p = pr_get_p(pr), u = zk_galoisapplymod(nf, pr_get_gen(pr), aut, p);
  return pr_make(nf, p, u, gel(pr,3), gel(pr,4));
}
static GEN
pr_galoismatrixapply(GEN nf, GEN pr, GEN M)
{
  GEN p = pr_get_p(pr), u = FpC_red(ZM_ZC_mul(M, pr_get_gen(pr)), p);
  return pr_make(nf, p, u, gel(pr,3), gel(pr,4));
}

static GEN
vecgaloisapply(GEN nf, GEN aut, GEN x)
{ pari_APPLY_same(galoisapply(nf, aut, gel(x,i))); }
static GEN
vecgaloismatrixapply(GEN nf, GEN aut, GEN x)
{ pari_APPLY_same(nfgaloismatrixapply(nf, aut, gel(x,i))); }

/* x: famat or standard algebraic number, aut automorphism in ZC form
 * simplified from general galoisapply */
static GEN
elt_galoisapply(GEN nf, GEN aut, GEN x)
{
  pari_sp av = avma;
  switch(typ(x))
  {
    case t_INT:  return icopy(x);
    case t_FRAC: return gcopy(x);
    case t_POLMOD: x = gel(x,2); /* fall through */
    case t_POL: {
      GEN y = basistoalg(nf, ZC_galoisapply(nf, aut, x));
      return gerepileupto(av,y);
    }
    case t_COL:
      return gerepileupto(av, ZC_galoisapply(nf, aut, x));
    case t_MAT:
      switch(lg(x)) {
        case 1: return cgetg(1, t_MAT);
        case 3: retmkmat2(vecgaloisapply(nf,aut,gel(x,1)), ZC_copy(gel(x,2)));
      }
  }
  pari_err_TYPE("galoisapply",x);
  return NULL; /* LCOV_EXCL_LINE */
}
/* M automorphism in matrix form */
static GEN
elt_galoismatrixapply(GEN nf, GEN M, GEN x)
{
  if (typ(x) == t_MAT)
    switch(lg(x)) {
      case 1: return cgetg(1, t_MAT);
      case 3: retmkmat2(vecgaloismatrixapply(nf,M,gel(x,1)), ZC_copy(gel(x,2)));
    }
  return nfgaloismatrixapply(nf, M, x);
}

GEN
galoisapply(GEN nf, GEN aut, GEN x)
{
  pari_sp av = avma;
  long lx;
  GEN y;

  nf = checknf(nf);
  switch(typ(x))
  {
    case t_INT:  return icopy(x);
    case t_FRAC: return gcopy(x);

    case t_POLMOD: x = gel(x,2); /* fall through */
    case t_POL:
      aut = algtobasis(nf, aut);
      y = basistoalg(nf, ZC_galoisapply(nf, aut, x));
      return gerepileupto(av,y);

    case t_VEC:
      aut = algtobasis(nf, aut);
      switch(lg(x))
      {
        case 6:
          if (pr_is_inert(x)) { set_avma(av); return gcopy(x); }
          return gerepilecopy(av, pr_galoisapply(nf, x, aut));
        case 3: y = cgetg(3,t_VEC);
          gel(y,1) = galoisapply(nf, aut, gel(x,1));
          gel(y,2) = elt_galoisapply(nf, aut, gel(x,2));
          return gerepileupto(av, y);
      }
      break;

    case t_COL:
      aut = algtobasis(nf, aut);
      return gerepileupto(av, ZC_galoisapply(nf, aut, x));

    case t_MAT: /* ideal */
      lx = lg(x); if (lx==1) return cgetg(1,t_MAT);
      if (nbrows(x) != nf_get_degree(nf)) break;
      y = RgM_mul(nfgaloismatrix(nf,aut), x);
      return gerepileupto(av, idealhnf_shallow(nf,y));
  }
  pari_err_TYPE("galoisapply",x);
  return NULL; /* LCOV_EXCL_LINE */
}

/* M automorphism in galoismatrix form */
GEN
nfgaloismatrixapply(GEN nf, GEN M, GEN x)
{
  pari_sp av = avma;
  long lx;
  GEN y;

  nf = checknf(nf);
  switch(typ(x))
  {
    case t_INT:  return icopy(x);
    case t_FRAC: return gcopy(x);

    case t_POLMOD: x = gel(x,2); /* fall through */
    case t_POL:
      x = algtobasis(nf, x);
      return gerepileupto(av, basistoalg(nf, RgM_RgC_mul(M, x)));

    case t_VEC:
      switch(lg(x))
      {
        case 6:
          if (pr_is_inert(x)) { set_avma(av); return gcopy(x); }
          return gerepilecopy(av, pr_galoismatrixapply(nf, x, M));
        case 3: y = cgetg(3,t_VEC);
          gel(y,1) = nfgaloismatrixapply(nf, M, gel(x,1));
          gel(y,2) = elt_galoismatrixapply(nf, M, gel(x,2));
          return gerepileupto(av, y);
      }
      break;

    case t_COL: return RgM_RgC_mul(M, x);

    case t_MAT: /* ideal */
      lx = lg(x); if (lx==1) return cgetg(1,t_MAT);
      if (nbrows(x) != nf_get_degree(nf)) break;
      return gerepileupto(av, idealhnf_shallow(nf,RgM_mul(M, x)));
  }
  pari_err_TYPE("galoisapply",x);
  return NULL; /* LCOV_EXCL_LINE */
}

/* compute action of automorphism s on nf.zk */
GEN
nfgaloismatrix(GEN nf, GEN s)
{
  pari_sp av2, av = avma;
  GEN zk, D, M, H, m;
  long k, n;

  nf = checknf(nf);
  zk = nf_get_zkprimpart(nf); n = lg(zk)-1;
  M = cgetg(n+1, t_MAT);
  gel(M,1) = col_ei(n, 1); /* s(1) = 1 */
  if (n == 1) return M;
  av2 = avma;
  if (typ(s) != t_COL) s = algtobasis(nf, s);
  D = nf_get_zkden(nf);
  H = RgV_to_RgM(zk, n);
  if (n == 2)
  {
    GEN t = gel(H,2); /* D * s(w_2) */
    t = ZC_Z_add(ZC_Z_mul(s, gel(t,2)), gel(t,1));
    gel(M,2) = gerepileupto(av2, gdiv(t, D));
    return M;
  }
  m = zk_multable(nf, s);
  gel(M,2) = s; /* M[,k] = s(x^(k-1)) */
  for (k = 3; k <= n; k++) gel(M,k) = ZM_ZC_mul(m, gel(M,k-1));
  M = ZM_mul(M, H);
  if (!equali1(D)) M = ZM_Z_divexact(M, D);
  return gerepileupto(av, M);
}

static GEN
get_aut(GEN nf, GEN gal, GEN aut, GEN g)
{
  return aut ? gel(aut, g[1]): poltobasis(nf, galoispermtopol(gal, g));
}

static GEN
idealquasifrob(GEN nf, GEN gal, GEN grp, GEN pr, GEN subg, GEN *S, GEN aut)
{
  pari_sp av = avma;
  long i, n = nf_get_degree(nf), f = pr_get_f(pr);
  GEN pi = pr_get_gen(pr);
  for (i=1; i<=n; i++)
  {
    GEN g = gel(grp,i);
    if ((!subg && perm_orderu(g) == (ulong)f)
      || (subg && perm_relorder(g, subg)==f))
    {
      *S = get_aut(nf, gal, aut, g);
      if (ZC_prdvd(ZC_galoisapply(nf, *S, pi), pr)) return g;
      set_avma(av);
    }
  }
  pari_err_BUG("idealquasifrob [Frobenius not found]");
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
nfgaloispermtobasis(GEN nf, GEN gal)
{
  GEN grp = gal_get_group(gal);
  long i, n = lg(grp)-1;
  GEN aut = cgetg(n+1, t_VEC);
  for(i=1; i<=n; i++)
  {
    pari_sp av = avma;
    GEN g = gel(grp, i);
    GEN vec = poltobasis(nf, galoispermtopol(gal, g));
    gel(aut, g[1]) = gerepileupto(av, vec);
  }
  return aut;
}

static void
gal_check_pol(const char *f, GEN x, GEN y)
{ if (!RgX_equal_var(x,y)) pari_err_MODULUS(f,x,y); }

/* true nf */
GEN
idealfrobenius_aut(GEN nf, GEN gal, GEN pr, GEN aut)
{
  pari_sp av = avma;
  GEN S=NULL, g=NULL; /*-Wall*/
  GEN T, p, a, b, modpr;
  long f, n, s;
  f = pr_get_f(pr); n = nf_get_degree(nf);
  if (f==1) { set_avma(av); return identity_perm(n); }
  g = idealquasifrob(nf, gal, gal_get_group(gal), pr, NULL, &S, aut);
  if (f==2) return gerepileuptoleaf(av, g);
  modpr = zk_to_Fq_init(nf,&pr,&T,&p);
  a = pol_x(nf_get_varn(nf));
  b = nf_to_Fq(nf, zk_galoisapplymod(nf, modpr_genFq(modpr), S, p), modpr);
  for (s = 1; s < f-1; s++)
  {
    a = Fq_pow(a, p, T, p);
    if (ZX_equal(a, b)) break;
  }
  g = perm_powu(g, Fl_inv(s, f));
  return gerepileupto(av, g);
}

GEN
idealfrobenius(GEN nf, GEN gal, GEN pr)
{
  nf = checknf(nf);
  checkgal(gal);
  checkprid(pr);
  gal_check_pol("idealfrobenius",nf_get_pol(nf),gal_get_pol(gal));
  if (pr_get_e(pr)>1) pari_err_DOMAIN("idealfrobenius","pr.e", ">", gen_1,pr);
  return idealfrobenius_aut(nf, gal, pr, NULL);
}

/* true nf */
GEN
idealramfrobenius_aut(GEN nf, GEN gal, GEN pr, GEN ram, GEN aut)
{
  pari_sp av = avma;
  GEN S=NULL, g=NULL; /*-Wall*/
  GEN T, p, a, b, modpr;
  GEN isog, deco;
  long f, n, s;
  f = pr_get_f(pr); n = nf_get_degree(nf);
  if (f==1) { set_avma(av); return identity_perm(n); }
  modpr = zk_to_Fq_init(nf,&pr,&T,&p);
  deco = group_elts(gel(ram,1), nf_get_degree(nf));
  isog = group_set(gel(ram,2),  nf_get_degree(nf));
  g = idealquasifrob(nf, gal, deco, pr, isog, &S, aut);
  a = pol_x(nf_get_varn(nf));
  b = nf_to_Fq(nf, zk_galoisapplymod(nf, modpr_genFq(modpr), S, p), modpr);
  for (s=0; !ZX_equal(a, b); s++)
    a = Fq_pow(a, p, T, p);
  g = perm_powu(g, Fl_inv(s, f));
  return gerepileupto(av, g);
}

GEN
idealramfrobenius(GEN nf, GEN gal, GEN pr, GEN ram)
{
  return idealramfrobenius_aut(nf, gal, pr, ram, NULL);
}

static GEN
idealinertiagroup(GEN nf, GEN gal, GEN aut, GEN pr)
{
  long i, n = nf_get_degree(nf);
  GEN p, T, modpr = zk_to_Fq_init(nf,&pr,&T,&p);
  GEN b = modpr_genFq(modpr);
  long e = pr_get_e(pr), coprime = ugcd(e, pr_get_f(pr)) == 1;
  GEN grp = gal_get_group(gal), pi = pr_get_gen(pr);
  pari_sp ltop = avma;
  for (i=1; i<=n; i++)
  {
    GEN iso = gel(grp,i);
    if (perm_orderu(iso) == (ulong)e)
    {
      GEN S = get_aut(nf, gal, aut, iso);
      if (ZC_prdvd(ZC_galoisapply(nf, S, pi), pr)
          && (coprime || gequalX(nf_to_Fq(nf, galoisapply(nf,S,b), modpr))))
          return iso;
      set_avma(ltop);
    }
  }
  pari_err_BUG("idealinertiagroup [no isotropic element]");
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN
idealramgroupstame(GEN nf, GEN gal, GEN aut, GEN pr)
{
  pari_sp av = avma;
  GEN iso, frob, giso, isog, S, res;
  long e = pr_get_e(pr), f = pr_get_f(pr);
  GEN grp = gal_get_group(gal);
  if (e == 1)
  {
    if (f==1)
      return cgetg(1,t_VEC);
    frob = idealquasifrob(nf, gal, grp, pr, NULL, &S, aut);
    set_avma(av);
    res = cgetg(2, t_VEC);
    gel(res, 1) = cyclicgroup(frob, f);
    return res;
  }
  res = cgetg(3, t_VEC);
  av = avma;
  iso = idealinertiagroup(nf, gal, aut, pr);
  set_avma(av);
  giso = cyclicgroup(iso, e);
  gel(res, 2) = giso;
  if (f==1)
  {
    gel(res, 1) = giso;
    return res;
  }
  av = avma;
  isog = group_set(giso, nf_get_degree(nf));
  frob = idealquasifrob(nf, gal, grp, pr, isog, &S, aut);
  set_avma(av);
  gel(res, 1) = dicyclicgroup(iso,frob,e,f);
  return res;
}

/* true nf, p | e */
static GEN
idealramgroupswild(GEN nf, GEN gal, GEN aut, GEN pr)
{
  pari_sp av2, av = avma;
  GEN p, T, idx, g, gbas, pi, pibas, Dpi, modpr = zk_to_Fq_init(nf,&pr,&T,&p);
  long bound, i, vDpi, vDg, n = nf_get_degree(nf);
  long e = pr_get_e(pr);
  long f = pr_get_f(pr);
  ulong nt,rorder;
  GEN pg, ppi, grp = gal_get_group(gal);

  /* G_i = {s: v(s(pi) - pi) > i} trivial for i > bound;
   * v_pr(Diff) = sum_{i = 0}^{bound} (#G_i - 1) >= e-1 + bound*(p-1)*/
  bound = (idealval(nf, nf_get_diff(nf), pr) - (e-1)) / (itou(p)-1);
  (void) u_pvalrem(n,p,&nt);
  rorder = e*f*(n/nt);
  idx = const_vecsmall(n,-1);
  pg = NULL;
  vDg = 0;
  if (f == 1)
    g = gbas = NULL;
  else
  {
    GEN Dg;
    g = nf_to_scalar_or_alg(nf, modpr_genFq(modpr));
    if (!gequalX(g)) /* p | nf.index */
    {
      g = Q_remove_denom(g, &Dg);
      vDg = Z_pval(Dg,p);
      pg = powiu(p, vDg + 1);
      g = FpX_red(g, pg);
    }
    gbas = nf_to_scalar_or_basis(nf, g);
  }
  pi = nf_to_scalar_or_alg(nf, pr_get_gen(pr));
  pi = Q_remove_denom(pi, &Dpi);
  vDpi = Dpi ? Z_pval(Dpi, p): 0;
  ppi = powiu(p, vDpi + (bound + e)/e);
  pi = FpX_red(pi, ppi);
  pibas = nf_to_scalar_or_basis(nf, pi);
  av2 = avma;
  for (i = 2; i <= n; i++)
  {
    GEN S, Spi, piso, iso = gel(grp, i);
    long j, o, ix = iso[1];
    if (idx[ix] >= 0 || rorder % (o = (long)perm_orderu(iso))) continue;

    piso = iso;
    S = get_aut(nf, gal, aut, iso);
    Spi = FpX_FpC_nfpoleval(nf, pi, FpC_red(S, ppi), ppi);
    /* Computation made knowing that the valuation is <= bound + 1. Correct
     * to maximal value if reduction mod ppi altered this */
    idx[ix] = minss(bound+1, idealval(nf, gsub(Spi,pibas), pr) - e*vDpi);
    if (idx[ix] == 0) idx[ix] = -1;
    else if (g)
    {
      GEN Sg = pg? FpX_FpC_nfpoleval(nf, g, FpC_red(S, pg), pg): S;
      if (vDg)
      { if (nfval(nf, gsub(Sg, gbas), pr) - e*vDg <= 0) idx[ix] = 0; }
      else /* same, more efficient */
      { if (!ZC_prdvd(gsub(Sg, gbas), pr)) idx[ix] = 0; }
    }
    for (j = 2; j < o; j++)
    {
      piso = perm_mul(piso,iso);
      if (ugcd(j,o)==1) idx[ piso[1] ] = idx[ix];
    }
    set_avma(av2);
  }
  return gerepileuptoleaf(av, idx);
}

GEN
idealramgroups_aut(GEN nf, GEN gal, GEN pr, GEN aut)
{
  pari_sp av = avma;
  GEN tbl, idx, res, set, sub;
  long i, j, e, n, maxm, p;
  ulong et;
  nf = checknf(nf);
  checkgal(gal);
  checkprid(pr);
  gal_check_pol("idealramgroups",nf_get_pol(nf),gal_get_pol(gal));
  e = pr_get_e(pr); n = nf_get_degree(nf);
  p = itos(pr_get_p(pr));
  if (e%p) return idealramgroupstame(nf, gal, aut, pr);
  (void) u_lvalrem(e,p,&et);
  idx = idealramgroupswild(nf, gal, aut, pr);
  sub = group_subgroups(galois_group(gal));
  tbl = subgroups_tableset(sub, n);
  maxm = vecsmall_max(idx)+1;
  res = cgetg(maxm+1,t_VEC);
  set = zero_F2v(n); F2v_set(set,1);
  for(i=maxm; i>0; i--)
  {
    long ix;
    for(j=1;j<=n;j++)
      if (idx[j]==i-1)
        F2v_set(set,j);
    ix = tableset_find_index(tbl, set);
    if (ix==0) pari_err_BUG("idealramgroups");
    gel(res,i) = gel(sub, ix);
  }
  return gerepilecopy(av, res);
}

GEN
idealramgroups(GEN nf, GEN gal, GEN pr)
{
  return idealramgroups_aut(nf, gal, pr, NULL);
}

/* x = relative polynomial nf = absolute nf, bnf = absolute bnf */
GEN
get_bnfpol(GEN x, GEN *bnf, GEN *nf)
{
  *bnf = checkbnf_i(x);
  *nf  = checknf_i(x);
  if (*nf) x = nf_get_pol(*nf);
  if (typ(x) != t_POL) pari_err_TYPE("get_bnfpol",x);
  return x;
}

GEN
get_nfpol(GEN x, GEN *nf)
{
  if (typ(x) == t_POL) { *nf = NULL; return x; }
  *nf = checknf(x); return nf_get_pol(*nf);
}

static GEN
incl_disc(GEN nfa, GEN a, int nolocal)
{
  GEN d;
  if (nfa) return nf_get_disc(nfa);
  if (nolocal) return NULL;
  d = ZX_disc(a);
  if (!signe(d)) pari_err_IRREDPOL("nfisincl",a);
  return d;
}

static int
badp(GEN fa, GEN db, long q)
{
  GEN P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P);
  for (i = 1; i < l; i++)
    if (mod2(gel(E,i)) && !dvdii(db, powiu(gel(P,i),q))) return 1;
  return 0;
}

/* is isomorphism / inclusion (a \subset b) compatible with what we know about
 * basic invariants ? (degree, signature, discriminant); test for isomorphism
 * if fliso is set and for inclusion otherwise */
static int
tests_OK(GEN a, GEN nfa, GEN b, GEN nfb, long fliso)
{
  GEN da2, da, db, fa, P, U;
  long i, l, q, m = degpol(a), n = degpol(b);

  if (m <= 0) pari_err_IRREDPOL("nfisincl",a);
  if (n <= 0) pari_err_IRREDPOL("nfisincl",b);
  q = n / m; /* relative degree */
  if (fliso) { if (n != m) return 0; } else { if (n % m) return 0; }
  if (m == 1) return 1;

  /*local test expensive if n^2 >> m^4 <=> q = n/m >> m */
  db = incl_disc(nfb, b, q > m);
  da = db? incl_disc(nfa, a, 0): NULL;
  if (nfa && nfb) /* both nf structures available */
  {
    long r1a = nf_get_r1(nfa), r1b = nf_get_r1(nfb);
    return fliso ? (r1a == r1b && equalii(da, db))
                 : (r1b <= r1a * q && dvdii(db, powiu(da, q)));
  }
  if (!db) return 1;
  if (fliso) return issquare(gdiv(da,db));

  if (nfa)
  {
    P = nf_get_ramified_primes(nfa); l = lg(P);
    for (i = 1; i < l; i++)
      if (Z_pval(db, gel(P,i)) < q * Z_pval(da, gel(P,i))) return 0;
    return 1;
  }
  else if (nfb)
  {
    P = nf_get_ramified_primes(nfb); l = lg(P);
    for (i = 1; i < l; i++)
    {
      GEN p = gel(P,i);
      long va = Z_pval(nfdisc(mkvec2(a, mkvec(p))), p);
      if (va && Z_pval(db, gel(P,i)) < va * q) return 0;
    }
    return 1;
  }
  /* da = dK A^2, db = dL B^2, dL = dK^q * N(D)
   * da = da1 * da2, da2 maximal s.t. (da2, db) = 1: let p a prime divisor of
   * da2 then p \nmid da1 * dK and p | A => v_p(da) = v_p(da2) is even */
  da2 = Z_ppo(da, db);
  if (!is_pm1(da2))
  { /* replace da by da1 all of whose prime divisors divide db */
    da2 = absi_shallow(da2);
    if (!Z_issquare(da2)) return 0;
    da = diviiexact(da, da2);
  }
  if (is_pm1(da)) return 1;
  fa = absZ_factor_limit_strict(da, 0, &U);
  if (badp(fa, db, q)) return 0;
  if (U && mod2(gel(U,2)) && expi(gel(U,1)) < 150)
  { /* cofactor is small, finish */
    fa = absZ_factor(gel(U,1));
    if (badp(fa, db, q)) return 0;
  }
  return 1;
}

GEN
nfisisom(GEN a, GEN b)
{
  pari_sp av = avma;
  long i, va, vb, lx;
  GEN nfa, nfb, y, la, lb;
  int newvar, sw = 0;

  a = get_nfpol(a, &nfa);
  b = get_nfpol(b, &nfb);
  if (!nfa) { a = Q_primpart(a); RgX_check_ZX(a, "nfisisom"); }
  if (!nfb) { b = Q_primpart(b); RgX_check_ZX(b, "nfisisom"); }
  if (ZX_equal(a, b))
  {
    y = galoisconj(nfb? nfb: b, NULL); settyp(y, t_VEC);
    return gerepilecopy(av,y);
  }
  if (nfa && !nfb) { swap(a,b); nfb = nfa; nfa = NULL; sw = 1; }
  if (!tests_OK(a, nfa, b, nfb, 1)) { set_avma(av); return gen_0; }

  if (nfb) lb = gen_1; else nfb = b = ZX_Q_normalize(b,&lb);
  if (nfa) la = gen_1; else nfa = a = ZX_Q_normalize(a,&la);
  va = varn(a); vb = varn(b); newvar = (varncmp(vb,va) <= 0);
  if (newvar) { a = leafcopy(a); setvarn(a, fetch_var_higher()); }
  y = lift_shallow(nfroots(nfb,a));
  if (newvar) (void)delete_var();
  lx = lg(y); if (lx==1) { set_avma(av); return gen_0; }
  if (sw) { vb = va; b = leafcopy(b); setvarn(b, vb); }
  for (i=1; i<lx; i++)
  {
    GEN t = gel(y,i);
    if (typ(t) == t_POL) setvarn(t, vb); else t = scalarpol(t, vb);
    if (lb != gen_1) t = RgX_unscale(t, lb);
    if (la != gen_1) t = RgX_Rg_div(t, la);
    gel(y,i) = sw? RgXQ_reverse(t, b): t;
  }
  settyp(y, t_VEC); return gerepilecopy(av,y);
}

static GEN
partmap_reverse(GEN a, GEN b, GEN t, GEN la, GEN lb, long v)
{
  pari_sp av = avma;
  GEN rnf = rnfequation2(a, t), z;
  if (!RgX_equal(b, gel(rnf,1)))
    { setvarn(b,v); pari_err_IRREDPOL("nfisincl", b); }
  z = liftpol_shallow(gel(rnf, 2));
  setvarn(z, v);
  if (!isint1(lb)) z = RgX_unscale(z, lb);
  if (!isint1(la)) z = RgX_Rg_div(z, la);
  return gerepilecopy(av, z);
}

static GEN
partmap_reverse_frac(GEN a, GEN b, GEN t, GEN la, GEN lb, long v)
{
  pari_sp av = avma;
  long k = 0;
  GEN N, D, G, L, de;
  GEN C = ZX_ZXY_resultant_all(a, Q_remove_denom(t,&de), &k, &L);
  if (k || degpol(b) != degpol(C))
    { setvarn(b,v); pari_err_IRREDPOL("nfisincl", b); }
  L =  Q_primpart(L);
  N = gel(L,1); if (!signe(N)) { set_avma(av); return pol_0(v); }
  D = gel(L,2);
  N = RgX_neg(N); setvarn(N, v); setvarn(D, v);
  G = QX_gcd(N,D);
  if (degpol(G)) { N = RgX_div(N,G); D = RgX_div(D,G); }
  if (!isint1(lb)) { N = RgX_unscale(N, lb); D = RgX_unscale(D, lb); }
  if (!isint1(la)) D = RgX_Rg_mul(D, la);
  return gerepilecopy(av, mkrfrac(N,D));
}

GEN
partmap_reverse_frac_worker(GEN t, GEN a, GEN b, GEN la, GEN lb, long v)
{ return partmap_reverse_frac(a, b, t, la, lb, v); }

static GEN
nfisincl_from_fact(GEN a, long da, GEN b, GEN la, GEN lb, long vb, GEN y,
                   long flag)
{
  long i, k, l = lg(y), db = degpol(b), d = db / da;
  GEN x = cgetg(l, t_VEC);
  for (i= k = 1; i < l; i++)
  {
    GEN t = gel(y,i);
    if (degpol(t) != d) continue;
    gel(x, k++) = partmap_reverse(a, b, t, la, lb, vb);
    if (flag) return gel(x,1);
  }
  if (k==1) return gen_0;
  setlg(x, k);
  gen_sort_inplace(x, (void*)&cmp_RgX, &cmp_nodata, NULL);
  return x;
}

static GEN
nfisincl_from_fact_frac(GEN a, GEN b, GEN la, GEN lb, long vb, GEN y)
{
  long i, k, l = lg(y), d = degpol(b) / degpol(a);
  GEN worker, x = cgetg(l, t_VEC);
  for (i = k = 1; i < l; i++)
  {
    GEN t = gel(y,i);
    if (degpol(t) != d) continue;
    gel(x, k++) = t;
  }
  if (k==1) return gen_0;
  worker = snm_closure(is_entry("_partmap_reverse_frac_worker"),
                       mkvec5(a,b,la,lb,stoi(vb)));
  setlg(x, k); return gen_parapply(worker, x);
}

GEN
nfisincl0(GEN fa, GEN fb, long flag)
{
  pari_sp av = avma;
  long vb;
  GEN a, b, nfa, nfb, x, y, la, lb;
  int newvar;
  if (flag < 0 || flag > 3) pari_err_FLAG("nfisincl");

  a = get_nfpol(fa, &nfa);
  b = get_nfpol(fb, &nfb);
  if (!nfa) { a = Q_primpart(a); RgX_check_ZX(a, "nsisincl"); }
  if (!nfb) { b = Q_primpart(b); RgX_check_ZX(b, "nsisincl"); }
  if (ZX_equal(a, b) && flag<=1)
  {
    if (flag==1)
    {
      x = pol_x(varn(b));
      return degpol(b) > 1 ? x: RgX_rem(x,b);
    }
    x = galoisconj(fb, NULL); settyp(x, t_VEC);
    return gerepilecopy(av,x);
  }
  if (flag==0 && !tests_OK(a, nfa, b, nfb, 0)) { set_avma(av); return gen_0; }

  if (nfb) lb = gen_1; else nfb = b = ZX_Q_normalize(b,&lb);
  if (nfa) la = gen_1; else nfa = a = ZX_Q_normalize(a,&la);
  vb = varn(b); newvar = (varncmp(varn(a),vb) <= 0);
  if (newvar) { b = leafcopy(b); setvarn(b, fetch_var_higher()); }
  y = lift_shallow(gel(nffactor(nfa,b),1));
  if (flag==2)
    x = nfisincl_from_fact_frac(a, b, la, lb, vb, y);
  else
    x = nfisincl_from_fact(nfa, degpol(a), b, la, lb, vb, y, flag);
  if (newvar) (void)delete_var();
  return gerepilecopy(av,x);
}

GEN
nfisincl(GEN fa, GEN fb) { return nfisincl0(fa, fb, 0); }

static GEN
RgF_to_Flxq(GEN F, GEN T, ulong p)
{
  GEN N, D, iD;
  if (typ(F)==t_POL) return RgX_to_Flx(F, p);
  N = RgX_to_Flx(gel(F,1), p); D = RgX_to_Flx(gel(F,2), p);
  iD = Flxq_invsafe(D, T, p);
  if (!iD) return NULL;
  return Flxq_mul(N, iD, T, p);
}

#define pari_APPLY_abort(EXPR)\
  { \
    long i, _l; \
    GEN _y = cgetg_copy(x, &_l);\
    for (i=1; i<_l; i++) \
    { GEN _z = EXPR;\
      if (!_z) return _z;\
       gel(_y,i) = _z;\
    } return _y;\
  }

static GEN
RgFV_to_FlxqV(GEN x, GEN T, ulong p)
{ pari_APPLY_abort(RgF_to_Flxq(gel(x,i), T, p)) }

static GEN
nfsplitting_auto(GEN g, GEN R)
{
  pari_sp av;
  forprime_t T;
  long i, d = degpol(g);
  ulong p;
  GEN P, K, N, G, q, den = Q_denom(R), Rp, Gp;
  u_forprime_init(&T, d*maxss(expu(d)-3, 2), ULONG_MAX);
  av = avma;
  for(;; set_avma(av))
  {
    p = u_forprime_next(&T);
    if (dvdiu(den,p)) continue;
    Gp = ZX_to_Flx(g, p);
    if (!Flx_is_totally_split(Gp, p)) continue;
    P = Flx_roots(Gp, p);
    Rp = RgFV_to_FlxqV(R, Gp, p);
    if (Rp) break;
    if (DEBUGLEVEL) err_printf("nfsplitting_auto: bad p : %lu\n",p);
  }
  if (d == 1) return mkvec3(g, mkcol(gel(Rp,1)), utoi(p));
  K = Flm_Flc_invimage(FlxV_to_Flm(Rp, d), vecsmall_ei(d, 2), p);
  if (!K) pari_err_BUG("nfsplitting_auto");
  N = Flm_transpose(FlxV_Flv_multieval(Rp, P, p));
  q = perm_inv(vecsmall_indexsort(gel(N,1)));
  G = cgetg(d+1, t_COL);
  for (i=1; i<=d; i++)
  {
    GEN r = perm_mul(vecsmall_indexsort(gel(N,i)), q);
    gel(G,i) = FlxV_Flc_mul(Rp, vecpermute(K, r), p);
  }
  return mkvec3(g, G, utoi(p));
}

static GEN
nfsplitting_composite(GEN P)
{
  GEN F = gel(ZX_factor(P), 1), Q = NULL;
  long i, n = lg(F)-1;
  for (i = 1; i <= n; i++)
  {
    GEN Fi = gel(F, i);
    if (degpol(Fi) == 1) continue;
    Q = Q ? veclast(compositum(Q, Fi)): Fi;
  }
  return Q ? Q: pol_x(varn(P));
}
GEN
nfsplitting0(GEN T0, GEN D, long flag)
{
  pari_sp av = avma;
  long d, Ds, v;
  GEN T, F, K, N = NULL, lT = NULL;
  if (flag < 0 || flag > 3) pari_err_FLAG("nfsplitting");
  T = T0 = get_nfpol(T0, &K);
  if (!K)
  {
    GEN c;
    if (typ(T) != t_POL) pari_err_TYPE("nfsplitting",T);
    T = Q_primitive_part(T, &c);
    lT = leading_coeff(T); if (isint1(lT)) lT = NULL;
    if (flag && (c || lT)) pari_err_TYPE("nfsplitting", T0);
    RgX_check_ZX(T,"nfsplitting");
  }
  T = nfsplitting_composite(T);
  if (flag && !ZX_equal(T, T0)) pari_err_IRREDPOL("nfsplitting", T0);
  d = degpol(T); v = varn(T);
  if (d <= 1 && !flag) { set_avma(av); return pol_x(v); }
  if (!K) {
    if (lT) T = polredbest(T,0);
    K = T;
  }
  if (D)
  { if (typ(D) != t_INT || signe(D) < 1) pari_err_TYPE("nfsplitting",D); }
  else if (d <= 7 ||
           (d <= 11 && pari_is_dir(stack_strcat(pari_datadir, "/galdata"))))
    D = gel(polgalois(T,DEFAULTPREC), 1);
  else
    D = mpfact(d);
  Ds = itos_or_0(D);
  T = leafcopy(T); setvarn(T, fetch_var_higher());
  for(F = T;;)
  {
    GEN P = gel(nffactor(K, F), 1), Q = veclast(P);
    if (degpol(gel(P,1)) == degpol(Q))
    {
      if (!flag) break;
      P = liftall_shallow(P);
      if (flag==1)
        N = nfisincl_from_fact(K, d, F, gen_1, gen_1, v, P, flag);
      else
        N = nfisincl_from_fact_frac(T0, F, gen_1, gen_1, v, P);
      break;
    }
    F = rnfequation(K,Q);
    if (degpol(F) == Ds && !flag) break;
  }
  if (umodiu(D,degpol(F)))
  {
    char *sD = itostr(D);
    pari_warn(warner,stack_strcat("ignoring incorrect degree bound ",sD));
  }
  setvarn(F, v); (void)delete_var();
  if (flag) F = flag == 3? nfsplitting_auto(F, N): mkvec2(F, N);
  return gerepilecopy(av, F);
}

GEN
nfsplitting(GEN T, GEN D) { return nfsplitting0(T, D, 0); }

/*************************************************************************/
/**                                                                     **/
/**                               INITALG                               **/
/**                                                                     **/
/*************************************************************************/
typedef struct {
  GEN T;
  GEN ro; /* roots of T */
  long r1;
  GEN basden;
  long prec;
  long extraprec; /* possibly -1 = irrelevant or not computed */
  GEN M, G; /* possibly NULL = irrelevant or not computed */
} nffp_t;

static GEN
get_roots(GEN x, long r1, long prec)
{
  long i, ru;
  GEN z;
  if (typ(x) != t_POL)
  {
    z = leafcopy(x);
    ru = (lg(z)-1 + r1) >> 1;
  }
  else
  {
    long n = degpol(x);
    z = (r1 == n)? ZX_realroots_irred(x, prec): QX_complex_roots(x,prec);
    ru = (n+r1)>>1;
  }
  for (i=r1+1; i<=ru; i++) gel(z,i) = gel(z, (i<<1)-r1);
  z[0]=evaltyp(t_VEC)|evallg(ru+1); return z;
}

GEN
nf_get_allroots(GEN nf)
{
  return embed_roots(nf_get_roots(nf), nf_get_r1(nf));
}

/* For internal use. compute trace(x mod pol), sym=polsym(pol,deg(pol)-1) */
static GEN
quicktrace(GEN x, GEN sym)
{
  GEN p1 = gen_0;
  long i;

  if (typ(x) != t_POL) return gmul(x, gel(sym,1));
  if (signe(x))
  {
    sym--;
    for (i=lg(x)-1; i>1; i--)
      p1 = gadd(p1, gmul(gel(x,i),gel(sym,i)));
  }
  return p1;
}

static GEN
get_Tr(GEN mul, GEN x, GEN basden)
{
  GEN t, bas = gel(basden,1), den = gel(basden,2);
  long i, j, n = lg(bas)-1;
  GEN T = cgetg(n+1,t_MAT), TW = cgetg(n+1,t_COL), sym = polsym(x, n-1);

  gel(TW,1) = utoipos(n);
  for (i=2; i<=n; i++)
  {
    t = quicktrace(gel(bas,i), sym);
    if (den && gel(den,i)) t = diviiexact(t,gel(den,i));
    gel(TW,i) = t; /* tr(w[i]) */
  }
  gel(T,1) = TW;
  for (i=2; i<=n; i++)
  {
    gel(T,i) = cgetg(n+1,t_COL); gcoeff(T,1,i) = gel(TW,i);
    for (j=2; j<=i; j++) /* Tr(W[i]W[j]) */
      gcoeff(T,i,j) = gcoeff(T,j,i) = ZV_dotproduct(gel(mul,j+(i-1)*n), TW);
  }
  return T;
}

/* return [bas[i]*denom(bas[i]), denom(bas[i])], denom 1 is given as NULL */
static GEN
get_bas_den(GEN bas)
{
  GEN b,d,den, dbas = leafcopy(bas);
  long i, l = lg(bas);
  int power = 1;
  den = cgetg(l,t_VEC);
  for (i=1; i<l; i++)
  {
    b = Q_remove_denom(gel(bas,i), &d);
    gel(dbas,i) = b;
    gel(den,i) = d; if (d) power = 0;
  }
  if (power) den = NULL; /* power basis */
  return mkvec2(dbas, den);
}

/* return multiplication table for S->basis */
static GEN
nf_multable(nfmaxord_t *S, GEN invbas)
{
  GEN T = S->T, w = gel(S->basden,1), den = gel(S->basden,2);
  long i,j, n = degpol(T);
  GEN mul = cgetg(n*n+1,t_MAT);

  /* i = 1 split for efficiency, assume w[1] = 1 */
  for (j=1; j<=n; j++)
    gel(mul,j) = gel(mul,1+(j-1)*n) = col_ei(n, j);
  for (i=2; i<=n; i++)
    for (j=i; j<=n; j++)
    {
      pari_sp av = avma;
      GEN z = (i == j)? ZXQ_sqr(gel(w,i), T): ZXQ_mul(gel(w,i),gel(w,j), T);
      z = ZM_ZX_mul(invbas, z); /* integral column */
      if (den)
      {
        GEN d = mul_denom(gel(den,i), gel(den,j));
        if (d) z = ZC_Z_divexact(z, d);
      }
      gel(mul,j+(i-1)*n) = gel(mul,i+(j-1)*n) = gerepileupto(av,z);
    }
  return mul;
}

/* as get_Tr, mul_table not precomputed */
static GEN
make_Tr(nfmaxord_t *S)
{
  GEN T = S->T, w = gel(S->basden,1), den = gel(S->basden,2);
  long i,j, n = degpol(T);
  GEN c, t, d, M = cgetg(n+1,t_MAT), sym = polsym(T, n-1);

  /* W[i] = w[i]/den[i]; assume W[1] = 1, case i = 1 split for efficiency */
  c = cgetg(n+1,t_COL); gel(M,1) = c;
  gel(c, 1) = utoipos(n);
  for (j=2; j<=n; j++)
  {
    pari_sp av = avma;
    t = quicktrace(gel(w,j), sym);
    if (den)
    {
      d = gel(den,j);
      if (d) t = diviiexact(t, d);
    }
    gel(c,j) = gerepileuptoint(av, t);
  }
  for (i=2; i<=n; i++)
  {
    c = cgetg(n+1,t_COL); gel(M,i) = c;
    for (j=1; j<i ; j++) gel(c,j) = gcoeff(M,i,j);
    for (   ; j<=n; j++)
    {
      pari_sp av = avma;
      t = (i == j)? ZXQ_sqr(gel(w,i), T): ZXQ_mul(gel(w,i),gel(w,j), T);
      t = quicktrace(t, sym);
      if (den)
      {
        d = mul_denom(gel(den,i),gel(den,j));
        if (d) t = diviiexact(t, d);
      }
      gel(c,j) = gerepileuptoint(av, t); /* Tr (W[i]W[j]) */
    }
  }
  return M;
}

/* [bas[i]/den[i]]= integer basis. roo = real part of the roots */
static void
make_M(nffp_t *F, int trunc)
{
  GEN bas = gel(F->basden,1), den = gel(F->basden,2), ro = F->ro;
  GEN m, d, M;
  long i, j, l = lg(ro), n = lg(bas);
  M = cgetg(n,t_MAT);
  gel(M,1) = const_col(l-1, gen_1); /* bas[1] = 1 */
  for (j=2; j<n; j++) gel(M,j) = cgetg(l,t_COL);
  for (i=1; i<l; i++)
  {
    GEN r = gel(ro,i), ri;
    ri = (gexpo(r) > 1)? ginv(r): NULL;
    for (j=2; j<n; j++) gcoeff(M,i,j) = RgX_cxeval(gel(bas,j), r, ri);
  }
  if (den)
    for (j=2; j<n; j++)
    {
      d = gel(den,j); if (!d) continue;
      m = gel(M,j);
      for (i=1; i<l; i++) gel(m,i) = gdiv(gel(m,i), d);
    }

  if (trunc && gprecision(M) > F->prec)
  {
    M     = gprec_w(M, F->prec);
    F->ro = gprec_w(ro,F->prec);
  }
  F->M = M;
}

/* return G real such that G~ * G = T_2 */
static void
make_G(nffp_t *F)
{
  GEN G, M = F->M;
  long i, j, k, r1 = F->r1, l = lg(M);

  if (r1 == l-1) { F->G = M; return; }
  G = cgetg(l, t_MAT);
  for (j = 1; j < l; j++)
  {
    GEN g, m = gel(M,j);
    gel(G,j) = g = cgetg(l, t_COL);
    for (k = i = 1; i <= r1; i++) gel(g,k++) = gel(m,i);
    for (     ; k < l; i++)
    {
      GEN r = gel(m,i);
      if (typ(r) == t_COMPLEX)
      {
        GEN a = gel(r,1), b = gel(r,2);
        gel(g,k++) = mpadd(a, b);
        gel(g,k++) = mpsub(a, b);
      }
      else
      {
        gel(g,k++) = r;
        gel(g,k++) = r;
      }
    }
  }
  F->G = G;
}

static long
prec_fix(long prec)
{
#ifndef LONG_IS_64BIT
  /* make sure that default accuracy is the same on 32/64bit */
  if (odd(prec)) prec++;
#endif
  return prec;
}
static void
make_M_G(nffp_t *F, int trunc)
{
  long n, eBD, prec;
  if (F->extraprec < 0)
  { /* not initialized yet; compute roots so that absolute accuracy
     * of M & G >= prec */
    double er;
    n = degpol(F->T);
    eBD = 1 + gexpo(gel(F->basden,1));
    er  = F->ro? (1+gexpo(F->ro)): fujiwara_bound(F->T);
    if (er < 0) er = 0;
    F->extraprec = nbits2extraprec(n*er + eBD + log2(n));
  }
  prec = prec_fix(F->prec + F->extraprec);
  if (!F->ro || gprecision(gel(F->ro,1)) < prec)
    F->ro = get_roots(F->T, F->r1, prec);

  make_M(F, trunc);
  make_G(F);
}

static void
nffp_init(nffp_t *F, nfmaxord_t *S, long prec)
{
  F->T  = S->T;
  F->r1 = S->r1;
  F->basden = S->basden;
  F->ro = NULL;
  F->extraprec = -1;
  F->prec = prec;
}

/* let bas a t_VEC of QX giving a Z-basis of O_K. Return the index of the
 * basis. Assume bas[1] = 1 and that the leading coefficient of elements
 * of bas are of the form 1/b for a t_INT b */
static GEN
get_nfindex(GEN bas)
{
  pari_sp av = avma;
  long n = lg(bas)-1, i;
  GEN D, d, mat;

  /* assume bas[1] = 1 */
  D = gel(bas,1);
  if (! is_pm1(simplify_shallow(D))) pari_err_TYPE("get_nfindex", D);
  D = gen_1;
  for (i = 2; i <= n; i++)
  { /* after nfbasis, basis is upper triangular! */
    GEN B = gel(bas,i), lc;
    if (degpol(B) != i-1) break;
    lc = gel(B, i+1);
    switch (typ(lc))
    {
      case t_INT: continue;
      case t_FRAC: if (is_pm1(gel(lc,1)) ) {D = mulii(D, gel(lc,2)); continue;}
      default: pari_err_TYPE("get_nfindex", B);
    }
  }
  if (i <= n)
  { /* not triangular after all */
    bas = vecslice(bas,i,n);
    bas = Q_remove_denom(bas, &d);
    if (!d) return D;
    mat = RgV_to_RgM(bas, n);
    mat = rowslice(mat, i,n);
    D = mulii(D, diviiexact(powiu(d, n-i+1), absi_shallow(ZM_det(mat))));
  }
  return gerepileuptoint(av, D);
}
/* make sure all components of S are initialized */
static void
nfmaxord_complete(nfmaxord_t *S)
{
  if (!S->dT) S->dT = ZX_disc(S->T);
  if (!S->index)
  {
    if (S->dK) /* fast */
      S->index = sqrti( diviiexact(S->dT, S->dK) );
    else
      S->index = get_nfindex(S->basis);
  }
  if (!S->dK) S->dK = diviiexact(S->dT, sqri(S->index));
  if (S->r1 < 0) S->r1 = ZX_sturm_irred(S->T);
  if (!S->basden) S->basden = get_bas_den(S->basis);
}

GEN
nfmaxord_to_nf(nfmaxord_t *S, GEN ro, long prec)
{
  GEN nf = cgetg(10,t_VEC);
  GEN T = S->T, Tr, D, w, A, dA, MDI, mat = cgetg(9,t_VEC);
  long n = degpol(T);
  nffp_t F;
  nfmaxord_complete(S);
  nffp_init(&F,S,prec);
  F.ro = ro;
  make_M_G(&F, 0);

  gel(nf,1) = S->T;
  gel(nf,2) = mkvec2s(S->r1, (n - S->r1)>>1);
  gel(nf,3) = S->dK;
  gel(nf,4) = S->index;
  gel(nf,5) = mat;
  if (gprecision(gel(F.ro,1)) > prec) F.ro = gprec_wtrunc(F.ro, prec);
  gel(nf,6) = F.ro;
  w = S->basis;
  if (!is_pm1(S->index)) w = Q_remove_denom(w, NULL);
  gel(nf,7) = w;
  gel(nf,8) = ZM_inv(RgV_to_RgM(w,n), NULL);
  gel(nf,9) = nf_multable(S, nf_get_invzk(nf));
  gel(mat,1) = F.M;
  gel(mat,2) = F.G;

  Tr = get_Tr(gel(nf,9), T, S->basden);
  gel(mat,6) = A = ZM_inv(Tr, &dA); /* dA T^-1, primitive */
  A = ZM_hnfmodid(A, dA);
  /* CAVEAT: nf is not complete yet, but the fields needed for
   * idealtwoelt, zk_scalar_or_multable and idealinv are present ! */
  MDI = idealtwoelt(nf, A);
  gel(MDI,2) = zk_scalar_or_multable(nf, gel(MDI,2));
  gel(mat,7) = MDI;
  if (is_pm1(S->index))
  { /* principal ideal (T'), whose norm is |dK| */
    D = zk_scalar_or_multable(nf, ZX_deriv(T));
    if (typ(D) == t_MAT) D = ZM_hnfmod(D, absi_shallow(S->dK));
  }
  else
  {
    GEN c = diviiexact(dA, gcoeff(A,1,1));
    D = idealHNF_inv_Z(nf, A); /* (A\cap Z) / A */
    if (!is_pm1(c)) D = ZM_Z_mul(D, c);
  }
  gel(mat,3) = RM_round_maxrank(F.G);
  gel(mat,4) = Tr;
  gel(mat,5) = D;
  w = S->dKP; if (!w) { w = gel(absZ_factor(S->dK), 1); settyp(w, t_VEC); }
  gel(mat,8) = w; return nf;
}

static GEN
primes_certify(GEN dK, GEN dKP)
{
  long i, l = lg(dKP);
  GEN v, w, D = dK;
  v = vectrunc_init(l);
  w = vectrunc_init(l);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(dKP,i);
    vectrunc_append(isprime(p)? w: v, p);
    (void)Z_pvalrem(D, p, &D);
  }
  if (!is_pm1(D))
  {
    if (signe(D) < 0) D = negi(D);
    vectrunc_append(isprime(D)? w: v, D);
  }
  return mkvec2(v,w);
}
GEN
nfcertify(GEN nf)
{
  pari_sp av = avma;
  GEN vw;
  nf = checknf(nf);
  vw = primes_certify(nf_get_disc(nf), nf_get_ramified_primes(nf));
  return gerepilecopy(av, gel(vw,1));
}

/* set *pro to roots of S->T */
static GEN
get_red_G(nfmaxord_t *S, GEN *pro)
{
  pari_sp av = avma;
  GEN G, u, u0 = NULL;
  long prec, n = degpol(S->T);
  nffp_t F;

  prec = nbits2prec(n+32);
  nffp_init(&F, S, prec);
  for (;;)
  {
    F.prec = prec; make_M_G(&F, 0); G = F.G;
    if (u0) G = RgM_mul(G, u0);
    if (DEBUGLEVEL)
      err_printf("get_red_G: starting LLL, prec = %ld (%ld + %ld)\n",
                  prec + F.extraprec, prec, F.extraprec);
    if ((u = lllfp(G, 0.99, LLL_KEEP_FIRST)))
    {
      if (lg(u)-1 == n) break;
      /* singular ==> loss of accuracy */
      if (u0) u0 = gerepileupto(av, RgM_mul(u0,u));
      else    u0 = gerepilecopy(av, u);
    }
    prec = precdbl(prec) + nbits2extraprec(gexpo(u0));
    F.ro = NULL;
    if (DEBUGLEVEL) pari_warn(warnprec,"get_red_G", prec);
  }
  if (u0) u = RgM_mul(u0,u);
  *pro = F.ro; return u;
}

/* Compute an LLL-reduced basis for the integer basis of nf(T).
 * set *pro = roots of x if computed [NULL if not computed] */
static void
set_LLL_basis(nfmaxord_t *S, GEN *pro, long flag, double DELTA)
{
  GEN B = S->basis;
  long N = degpol(S->T);
  if (S->r1 < 0)
  {
    S->r1 = ZX_sturm_irred(S->T);
    if (odd(N - S->r1)) pari_err_IRREDPOL("set_LLL_basis", S->T);
  }
  if (!S->basden) S->basden = get_bas_den(B);
  *pro = NULL; if (flag & nf_NOLLL) return;
  if (S->r1 == N) {
    pari_sp av = avma;
    GEN u = ZM_lll(make_Tr(S), DELTA, LLL_GRAM|LLL_KEEP_FIRST|LLL_IM);
    B = gerepileupto(av, RgV_RgM_mul(B, u));
  }
  else
    B = RgV_RgM_mul(B, get_red_G(S, pro));
  S->basis = B;
  S->basden = get_bas_den(B);
}

/* = 1 iff |a| > |b| or equality and a > 0 */
static int
cmpii_polred(GEN a, GEN b)
{
  int fl = abscmpii(a, b);
  long sa, sb;
  if (fl) return fl;
  sa = signe(a);
  sb = signe(b);
  if (sa == sb) return 0;
  return sa == 1? 1: -1;
}
static int
ZX_cmp(GEN x, GEN y)
{  return gen_cmp_RgX((void*)cmpii_polred, x, y); }
/* current best: ZX x of discriminant *dx, is ZX y better than x ?
 * (if so update *dx); both x and y are monic */
static int
ZX_is_better(GEN y, GEN x, GEN *dx)
{
  pari_sp av;
  int cmp;
  GEN d;
  if (!*dx) *dx = ZX_disc(x);
  av = avma; d = ZX_disc(y);
  cmp = abscmpii(d, *dx);
  if (cmp < 0) { *dx = d; return 1; }
  return gc_bool(av, cmp? 0: (ZX_cmp(y, x) < 0));
}

static void polredbest_aux(nfmaxord_t *S, GEN *pro, GEN *px, GEN *pdx, GEN *pa);
/* Seek a simpler, polynomial pol defining the same number field as
 * x (assumed to be monic at this point) */
static GEN
nfpolred(nfmaxord_t *S, GEN *pro)
{
  GEN x = S->T, dx, b, rev;
  long n = degpol(x), v = varn(x);

  if (n == 1) {
    S->T = pol_x(v);
    *pro = NULL;
    return scalarpol_shallow(negi(gel(x,2)), v);
  }
  polredbest_aux(S, pro, &x, &dx, &b);
  if (x == S->T) return NULL; /* no improvement */
  if (DEBUGLEVEL>1) err_printf("xbest = %Ps\n",x);

  /* update T */
  rev = QXQ_reverse(b, S->T);
  S->basis = QXV_QXQ_eval(S->basis, rev, x);
  S->index = sqrti( diviiexact(dx,S->dK) );
  S->basden = get_bas_den(S->basis);
  S->dT = dx;
  S->T = x;
  *pro = NULL; /* reset */
  return rev;
}

/* Either nf type or ZX or [monic ZX, data], where data is either an integral
 * basis (deprecated), or listP data (nfbasis input format) to specify
 * a set of primes at with the basis order must be maximal.
 * 1) nf type (or unrecognized): return t_VEC
 * 2) ZX or [ZX, listP]: return t_POL
 * 3) [ZX, order basis]: return 0 (deprecated)
 * incorrect: return -1 */
static long
nf_input_type(GEN x)
{
  GEN T, V, DKP = NULL;
  long i, d, v;
  switch(typ(x))
  {
    case t_POL: return t_POL;
    case t_VEC:
      switch(lg(x))
      {
        case 4: DKP = gel(x,3);
        case 3: break;
        default: return t_VEC; /* nf or incorrect */
      }
      T = gel(x,1); V = gel(x,2);
      if (typ(T) != t_POL) return -1;
      switch(typ(V))
      {
        case t_INT: case t_MAT: return t_POL;
        case t_VEC: case t_COL:
          if (RgV_is_ZV(V)) return t_POL;
          break;
        default: return -1;
      }
      d = degpol(T); v = varn(T);
      if (d<1 || !RgX_is_ZX(T) || !isint1(gel(T,d+2)) || lg(V)-1!=d) return -1;
      for (i = 1; i <= d; i++)
      { /* check integer basis */
        GEN c = gel(V,i);
        switch(typ(c))
        {
          case t_INT: break;
          case t_POL: if (varn(c) == v && RgX_is_QX(c) && degpol(c) < d) break;
          /* fall through */
          default: return -1;
        }
      }
      if (DKP && (typ(DKP) != t_VEC || !RgV_is_ZV(DKP))) return -1;
      return 0;
  }
  return t_VEC; /* nf or incorrect */
}

/* cater for obsolete nf_PARTIALFACT flag */
static void
nfinit_basic_partial(nfmaxord_t *S, GEN T)
{
  if (typ(T) == t_POL) { nfmaxord(S, mkvec2(T,utoipos(500000)), 0); }
  else nfinit_basic(S, T);
}
static void
nfinit_basic_flag(nfmaxord_t *S, GEN x, long flag)
{
  if (flag & nf_PARTIALFACT)
    nfinit_basic_partial(S, x);
  else
    nfinit_basic(S, x);
}

/* true nf */
static GEN
nf_basden(GEN nf)
{
  GEN zkD = nf_get_zkprimpart(nf), D = nf_get_zkden(nf);
  D = equali1(D)? NULL: const_vec(lg(zkD)-1, D);
  return mkvec2(zkD, D);
}
void
nfinit_basic(nfmaxord_t *S, GEN T)
{
  switch (nf_input_type(T))
  {
    case t_POL: nfmaxord(S, T, 0); return;
    case t_VEC:
    { /* nf, bnf, bnr */
      GEN nf = checknf(T);
      S->T = S->T0 = nf_get_pol(nf);
      S->basis = nf_get_zk(nf); /* probably useless */
      S->basden = nf_basden(nf);
      S->index = nf_get_index(nf);
      S->dK = nf_get_disc(nf);
      S->dKP = nf_get_ramified_primes(nf);
      S->dT = mulii(S->dK, sqri(S->index));
      S->r1 = nf_get_r1(nf); break;
    }
    case 0: /* monic integral polynomial + integer basis (+ ramified primes)*/
      S->T = S->T0 = gel(T,1);
      S->basis = gel(T,2);
      S->basden = NULL;
      S->index = NULL;
      S->dK = NULL;
      S->dKP = NULL;
      if (lg(T) == 4)
      {
        GEN v = gel(T,3); if (typ(v) == t_COL) v = shallowtrans(v);
        S->dKP = v;
      }
      S->dT = NULL;
      S->r1 = -1; break;
    default: /* -1 */
      pari_err_TYPE("nfinit_basic", T);
  }
  S->dTP = S->dTE = S->dKE = NULL;
  S->unscale = gen_1;
}

GEN
nfinit_complete(nfmaxord_t *S, long flag, long prec)
{
  GEN nf, unscale = S->unscale, rev = NULL;

  if (!ZX_is_irred(S->T)) pari_err_IRREDPOL("nfinit",S->T);
  if (!(flag & nf_RED) && !ZX_is_monic(S->T0))
  {
    pari_warn(warner,"nonmonic polynomial. Result of the form [nf,c]");
    flag |= nf_RED | nf_ORIG;
  }
  if (!(flag & nf_RED) && !isint1(unscale))
  { /* implies lc(x0) = 1 and L := 1/unscale is integral */
    long d = degpol(S->T0);
    GEN L = ginv(unscale); /* x = L^(-deg(x)) x0(L X) */
    GEN f = powiu(L, (d*(d-1)) >> 1);
    S->T = S->T0; /* restore original user-supplied x0, unscale data */
    S->unscale = gen_1;
    S->dT    = gmul(S->dT, sqri(f));
    S->basis = RgXV_unscale(S->basis, unscale);
    S->index = gmul(S->index, f);
  }
  nfmaxord_complete(S); /* more expensive after set_LLL_basis */
  if (flag & nf_RED)
  {
    GEN ro;
    /* lie to polred: more efficient to update *after* modreverse, than to
     * unscale in the polred subsystem */
    S->unscale = gen_1;
    rev = nfpolred(S, &ro);
    nf = nfmaxord_to_nf(S, ro, prec);
    S->unscale = unscale; /* restore */
  }
  else
  {
    GEN ro; set_LLL_basis(S, &ro, flag, 0.99);
    nf = nfmaxord_to_nf(S, ro, prec);
  }
  if (flag & nf_ORIG)
  {
    if (!rev)
    { /* no improvement */
      long v = varn(S->T);
      rev = degpol(S->T) == 1? pol_0(v): pol_x(v);
    }
    if (!isint1(unscale)) rev = RgX_Rg_div(rev, unscale);
    nf = mkvec2(nf, mkpolmod(rev, S->T));
  }
  return nf;
}
/* Initialize the number field defined by the polynomial x (in variable v)
 * flag & nf_RED: try a polred first.
 * flag & nf_ORIG: return [nfinit(x), Mod(a,red)], where
 *    Mod(a,red) = Mod(v,x) (i.e return the base change). */
GEN
nfinit0(GEN x, long flag,long prec)
{
  const pari_sp av = avma;
  nfmaxord_t S;
  if (flag < 0 || flag > 7) pari_err_FLAG("nfinit");
  if (checkrnf_i(x)) return rnf_build_nfabs(x, prec);
  nfinit_basic(&S, x);
  return gerepilecopy(av, nfinit_complete(&S, flag, prec));
}
GEN
nfinitred(GEN x, long prec)  { return nfinit0(x, nf_RED, prec); }
GEN
nfinitred2(GEN x, long prec) { return nfinit0(x, nf_RED|nf_ORIG, prec); }
GEN
nfinit(GEN x, long prec)     { return nfinit0(x, 0, prec); }

/* assume x a bnr/bnf/nf */
long
nf_get_prec(GEN x)
{
  GEN nf = checknf(x), ro = nf_get_roots(nf);
  return (typ(ro)==t_VEC)? precision(gel(ro,1)): DEFAULTPREC;
}

/* true nf */
GEN
nfnewprec_shallow(GEN nf, long prec)
{
  GEN m, NF = leafcopy(nf);
  nffp_t F;

  F.T  = nf_get_pol(nf);
  F.ro = NULL;
  F.r1 = nf_get_r1(nf);
  F.basden = nf_basden(nf);
  F.extraprec = -1;
  F.prec = prec; make_M_G(&F, 0);
  gel(NF,5) = m = leafcopy(gel(NF,5));
  gel(m,1) = F.M;
  gel(m,2) = F.G;
  gel(NF,6) = F.ro; return NF;
}

GEN
nfnewprec(GEN nf, long prec)
{
  GEN z;
  switch(nftyp(nf))
  {
    default: pari_err_TYPE("nfnewprec", nf);
    case typ_BNF: z = bnfnewprec(nf,prec); break;
    case typ_BNR: z = bnrnewprec(nf,prec); break;
    case typ_NF: {
      pari_sp av = avma;
      z = gerepilecopy(av, nfnewprec_shallow(checknf(nf), prec));
      break;
    }
  }
  return z;
}

/********************************************************************/
/**                                                                **/
/**                           POLRED                               **/
/**                                                                **/
/********************************************************************/
GEN
embednorm_T2(GEN x, long r1)
{
  pari_sp av = avma;
  GEN p = RgV_sumpart(x, r1);
  GEN q = RgV_sumpart2(x,r1+1, lg(x)-1);
  if (q != gen_0) p = gadd(p, gmul2n(q,1));
  return avma == av? gcopy(p): gerepileupto(av, p);
}

/* simplified version of gnorm for scalar, noncomplex inputs, without GC */
static GEN
real_norm(GEN x)
{
  switch(typ(x))
  {
    case t_INT:  return sqri(x);
    case t_REAL: return sqrr(x);
    case t_FRAC: return sqrfrac(x);
  }
  pari_err_TYPE("real_norm", x);
  return NULL;/*LCOV_EXCL_LINE*/
}
/* simplified version of gnorm, without GC */
static GEN
complex_norm(GEN x)
{
  return typ(x) == t_COMPLEX? cxnorm(x): real_norm(x);
}
/* return T2(x), argument r1 needed in case x has components whose type
 * is unexpected, e.g. all of them t_INT for embed(gen_1) */
GEN
embed_T2(GEN x, long r1)
{
  pari_sp av = avma;
  long i, l = lg(x);
  GEN c, s = NULL, t = NULL;
  if (typ(gel(x,1)) == t_INT) return muliu(gel(x,1), 2*(l-1)-r1);
  for (i = 1; i <= r1; i++)
  {
    c = real_norm(gel(x,i));
    s = s? gadd(s, c): c;
  }
  for (; i < l; i++)
  {
    c = complex_norm(gel(x,i));
    t = t? gadd(t, c): c;
  }
  if (t) { t = gmul2n(t,1); s = s? gadd(s,t): t; }
  return gerepileupto(av, s);
}
/* return N(x) */
GEN
embed_norm(GEN x, long r1)
{
  pari_sp av = avma;
  long i, l = lg(x);
  GEN c, s = NULL, t = NULL;
  if (typ(gel(x,1)) == t_INT) return powiu(gel(x,1), 2*(l-1)-r1);
  for (i = 1; i <= r1; i++)
  {
    c = gel(x,i);
    s = s? gmul(s, c): c;
  }
  for (; i < l; i++)
  {
    c = complex_norm(gel(x,i));
    t = t? gmul(t, c): c;
  }
  if (t) s = s? gmul(s,t): t;
  return gerepileupto(av, s);
}

typedef struct {
  long r1, v, prec;
  GEN ZKembed; /* embeddings of fincke-pohst-reduced Zk basis */
  GEN u; /* matrix giving fincke-pohst-reduced Zk basis */
  GEN M; /* embeddings of initial (LLL-reduced) Zk basis */
  GEN bound; /* T2 norm of the polynomial defining nf */
  long expo_best_disc; /* expo(disc(x)), best generator so far */
} CG_data;

/* characteristic pol of x (given by embeddings) */
static GEN
get_pol(CG_data *d, GEN x)
{
  long e;
  GEN g = grndtoi(roots_to_pol_r1(x, d->v, d->r1), &e);
  return (e > -5)? NULL: g;
}

/* characteristic pol of x (given as vector on (w_i)) */
static GEN
get_polchar(CG_data *d, GEN x)
{ return get_pol(d, RgM_RgC_mul(d->ZKembed,x)); }

/* Choose a canonical polynomial in the pair { Pmin_a, Pmin_{-a} }, i.e.
 * { z(X), (-1)^(deg z) z(-Z) } and keeping the smallest wrt cmpii_polred
 * Either leave z alone (return 1) or set z <- (-1)^n z(-X). In place. */
int
ZX_canon_neg(GEN z)
{
  long i, s;
  for (i = lg(z)-2; i >= 2; i -= 2)
  { /* examine the odd (resp. even) part of z if deg(z) even (resp. odd). */
    s = signe(gel(z,i));
    if (!s) continue;
    /* non trivial */
    if (s < 0) break; /* z(X) < (-1)^n z(-X) */

    for (; i>=2; i-=2) gel(z,i) = negi(gel(z,i));
    return 1;
  }
  return 0;
}
/* return a defining polynomial for Q(alpha), v = embeddings of alpha.
 * Return NULL on failure: discriminant too large or non primitive */
static GEN
try_polmin(CG_data *d, nfmaxord_t *S, GEN v, long flag, GEN *ai)
{
  const long best = flag & nf_ABSOLUTE;
  long ed;
  pari_sp av = avma;
  GEN g;
  if (best)
  {
    ed = expo(embed_disc(v, d->r1, LOWDEFAULTPREC));
    set_avma(av); if (d->expo_best_disc < ed) return NULL;
  }
  else
    ed = 0;
  g = get_pol(d, v);
  /* accuracy too low, compute algebraically */
  if (!g) { set_avma(av); g = ZXQ_charpoly(*ai, S->T, varn(S->T)); }
  g = ZX_radical(g);
  if (best && degpol(g) != degpol(S->T)) return gc_NULL(av);
  g = gerepilecopy(av, g);
  d->expo_best_disc = ed;
  if (flag & nf_ORIG)
  {
    if (ZX_canon_neg(g)) *ai = RgX_neg(*ai);
    if (!isint1(S->unscale)) *ai = RgX_unscale(*ai, S->unscale);
  }
  else
    (void)ZX_canon_neg(g);
  if (DEBUGLEVEL>3) err_printf("polred: generator %Ps\n", g);
  return g;
}

/* does x generate the correct field ? */
static GEN
chk_gen(void *data, GEN x)
{
  pari_sp av = avma, av1;
  GEN h, g = get_polchar((CG_data*)data,x);
  if (!g) pari_err_PREC("chk_gen");
  av1 = avma;
  h = ZX_gcd(g, ZX_deriv(g));
  if (degpol(h)) return gc_NULL(av);
  if (DEBUGLEVEL>3) err_printf("  generator: %Ps\n",g);
  set_avma(av1); return gerepileupto(av, g);
}

static long
chk_gen_prec(long N, long bit)
{ return prec_fix(nbits2prec(10 + (long)log2((double)N) + bit)); }

/* v = [P,A] two vectors (of ZX and ZV resp.) of same length; remove duplicate
 * polynomials in P, updating A, in place. Among elements having the same
 * characteristic pol, choose the smallest according to ZV_abscmp */
static void
remove_duplicates(GEN v)
{
  GEN x, a, P = gel(v,1), A = gel(v,2);
  long k, i, l = lg(P);
  pari_sp av = avma;

  if (l < 2) return;
  (void)sort_factor_pol(mkvec2(P, A), cmpii);
  x = gel(P,1); a = gel(A,1);
  for  (k=1,i=2; i<l; i++)
    if (ZX_equal(gel(P,i), x))
    {
      if (ZV_abscmp(gel(A,i), a) < 0) a = gel(A,i);
    }
    else
    {
      gel(A,k) = a;
      gel(P,k) = x;
      k++;
      x = gel(P,i); a = gel(A,i);
    }
  l = k+1;
  gel(A,k) = a; setlg(A,l);
  gel(P,k) = x; setlg(P,l); set_avma(av);
}

static void
polred_init(nfmaxord_t *S, nffp_t *F, CG_data *d)
{
  long e, prec, n = degpol(S->T);
  double log2rho;
  GEN ro;
  set_LLL_basis(S, &ro, 0, 0.9999);
  /* || polchar ||_oo < 2^e ~ 2 (n * rho)^n, rho = max modulus of root */
  log2rho = ro ? (double)gexpo(ro): fujiwara_bound(S->T);
  e = n * (long)(log2rho + log2((double)n)) + 1;
  if (e < 0) e = 0; /* can occur if n = 1 */
  prec = chk_gen_prec(n, e);
  nffp_init(F,S,prec);
  F->ro = ro;
  make_M_G(F, 1);

  d->v = varn(S->T);
  d->expo_best_disc = -1;
  d->ZKembed = NULL;
  d->M = NULL;
  d->u = NULL;
  d->r1= S->r1;
}
static GEN
findmindisc(GEN y)
{
  GEN x = gel(y,1), dx = NULL;
  long i, l = lg(y);
  for (i = 2; i < l; i++)
  {
    GEN yi = gel(y,i);
    if (ZX_is_better(yi,x,&dx)) x = yi;
  }
  return x;
}
/* filter [y,b] from polred_aux: keep a single polynomial of degree n in y
 * [ the best wrt discriminant ordering ], but keep all imprimitive
 * polynomials */
static void
filter(GEN y, GEN b, long n)
{
  GEN x, a, dx;
  long i, k = 1, l = lg(y);
  a = x = dx = NULL;
  for (i = 1; i < l; i++)
  {
    GEN yi = gel(y,i), ai = gel(b,i);
    if (degpol(yi) == n)
    {
      if (!dx) dx = ZX_disc(yi); else if (!ZX_is_better(yi,x,&dx)) continue;
      x = yi; a = ai; continue;
    }
    gel(y,k) = yi;
    gel(b,k) = ai; k++;
  }
  if (dx)
  {
    gel(y,k) = x;
    gel(b,k) = a; k++;
  }
  setlg(y, k);
  setlg(b, k);
}

static GEN
polred_aux(nfmaxord_t *S, GEN *pro, long flag)
{ /* only keep polynomials of max degree and best discriminant */
  const long best = flag & nf_ABSOLUTE;
  const long orig = flag & nf_ORIG;
  GEN M, b, y, x = S->T;
  long maxi, i, j, k, v = varn(x), n = lg(S->basis)-1;
  nffp_t F;
  CG_data d;

  if (n == 1)
  {
    if (!best)
    {
      GEN X = pol_x(v);
      return orig? mkmat2(mkcol(X),mkcol(gen_1)): mkvec(X);
    }
    else
      return orig? trivial_fact(): cgetg(1,t_VEC);
  }

  polred_init(S, &F, &d);
  if (pro) *pro = F.ro;
  M = F.M;
  if (best)
  {
    if (!S->dT) S->dT = ZX_disc(S->T);
    d.expo_best_disc = expi(S->dT);
  }

  /* n + 2 sum_{1 <= i <= n} n-i = n + n(n-1) = n*n */
  y = cgetg(n*n + 1, t_VEC);
  b = cgetg(n*n + 1, t_COL);
  k = 1;
  if (!best) { gel(y,1) = pol_x(v); gel(b,1) = gen_0; k++; }
  for (i = 2; i <= n; i++)
  {
    GEN ch, ai;
    ai = gel(S->basis,i);
    ch = try_polmin(&d, S, gel(M,i), flag, &ai);
    if (ch) { gel(y,k) = ch; gel(b,k) = ai; k++; }
  }
  maxi = minss(n, 3);
  for (i = 1; i <= maxi; i++)
    for (j = i+1; j <= n; j++)
    {
      GEN ch, ai, v;
      ai = gadd(gel(S->basis,i), gel(S->basis,j));
      v = RgV_add(gel(M,i), gel(M,j));
      /* defining polynomial for Q(w_i+w_j) */
      ch = try_polmin(&d, S, v, flag, &ai);
      if (ch) { gel(y,k) = ch; gel(b,k) = ai; k++; }

      ai = gsub(gel(S->basis,i), gel(S->basis,j));
      v = RgV_sub(gel(M,i), gel(M,j));
      /* defining polynomial for Q(w_i-w_j) */
      ch = try_polmin(&d, S, v, flag, &ai);
      if (ch) { gel(y,k) = ch; gel(b,k) = ai; k++; }
    }
  setlg(y, k);
  setlg(b, k); filter(y, b, n);
  if (!orig) return gen_sort_uniq(y, (void*)cmpii, &gen_cmp_RgX);
  settyp(y, t_COL);
  (void)sort_factor_pol(mkmat2(y, b), cmpii);
  return mkmat2(b, y);
}

/* FIXME: obsolete */
static GEN
Polred(GEN x, long flag, GEN fa)
{
  pari_sp av = avma;
  nfmaxord_t S;
  if (fa)
    nfinit_basic(&S, mkvec2(x,fa));
  else
    nfinit_basic_flag(&S, x, flag);
  return gerepilecopy(av, polred_aux(&S, NULL, flag));
}

/* finds "best" polynomial in polred_aux list, defaulting to S->T if none of
 * them is primitive. *px a ZX, characteristic polynomial of Mod(*pb,S->T),
 * *pdx its discriminant if pdx != NULL. Set *pro = polroots(S->T) */
static void
polredbest_aux(nfmaxord_t *S, GEN *pro, GEN *px, GEN *pdx, GEN *pb)
{
  GEN y, dx, x = S->T; /* default value */
  long i, l;
  y = polred_aux(S, pro, pb? nf_ORIG|nf_ABSOLUTE: nf_ABSOLUTE);
  dx = S->dT;
  if (pb)
  {
    GEN a, b = deg1pol_shallow(S->unscale, gen_0, varn(x));
    a = gel(y,1); l = lg(a);
    y = gel(y,2);
    for (i=1; i<l; i++)
    {
      GEN yi = gel(y,i);
      if (ZX_is_better(yi,x,&dx)) { x = yi; b = gel(a,i); }
    }
    *pb = b;
  }
  else
  {
    l = lg(y);
    for (i=1; i<l; i++)
    {
      GEN yi = gel(y,i);
      if (ZX_is_better(yi,x,&dx)) x = yi;
    }
  }
  if (pdx) { if (!dx) dx = ZX_disc(x); *pdx = dx; }
  *px = x;
}
static GEN
polredbest_i(GEN T, long flag)
{
  nfmaxord_t S;
  GEN a;
  nfinit_basic_partial(&S, T);
  polredbest_aux(&S, NULL, &T, NULL, flag? &a: NULL);
  if (flag == 2)
    T = mkvec2(T, a);
  else if (flag == 1)
  {
    GEN b = (S.T0 == T)? pol_x(varn(T)): QXQ_reverse(a, S.T0);
    /* charpoly(Mod(a,T0)) = T; charpoly(Mod(b,T)) = S.x */
    if (degpol(T) == 1) b = grem(b,T);
    T = mkvec2(T, mkpolmod(b,T));
  }
  return T;
}
GEN
polredbest(GEN T, long flag)
{
  pari_sp av = avma;
  if (flag < 0 || flag > 1) pari_err_FLAG("polredbest");
  return gerepilecopy(av, polredbest_i(T, flag));
}
/* DEPRECATED: backward compatibility */
GEN
polred0(GEN x, long flag, GEN fa)
{
  long fl = 0;
  if (flag & 1) fl |= nf_PARTIALFACT;
  if (flag & 2) fl |= nf_ORIG;
  return Polred(x, fl, fa);
}

GEN
polredord(GEN x)
{
  pari_sp av = avma;
  GEN v, lt;
  long i, n, vx;

  if (typ(x) != t_POL) pari_err_TYPE("polredord",x);
  x = Q_primpart(x); RgX_check_ZX(x,"polredord");
  n = degpol(x); if (n <= 0) pari_err_CONSTPOL("polredord");
  if (n == 1) return gerepilecopy(av, mkvec(x));
  lt = leading_coeff(x); vx = varn(x);
  if (is_pm1(lt))
  {
    if (signe(lt) < 0) x = ZX_neg(x);
    v = pol_x_powers(n, vx);
  }
  else
  { GEN L;
    /* basis for Dedekind order */
    v = cgetg(n+1, t_VEC);
    gel(v,1) = scalarpol_shallow(lt, vx);
    for (i = 2; i <= n; i++)
      gel(v,i) = RgX_Rg_add(RgX_mulXn(gel(v,i-1), 1), gel(x,n+3-i));
    gel(v,1) = pol_1(vx);
    x = ZX_Q_normalize(x, &L);
    v = gsubst(v, vx, monomial(ginv(L),1,vx));
    for (i=2; i <= n; i++)
      if (Q_denom(gel(v,i)) == gen_1) gel(v,i) = pol_xn(i-1, vx);
  }
  return gerepileupto(av, polred(mkvec2(x, v)));
}

GEN
polred(GEN x) { return Polred(x, 0, NULL); }
GEN
smallpolred(GEN x) { return Polred(x, nf_PARTIALFACT, NULL); }
GEN
factoredpolred(GEN x, GEN fa) { return Polred(x, 0, fa); }
GEN
polred2(GEN x) { return Polred(x, nf_ORIG, NULL); }
GEN
smallpolred2(GEN x) { return Polred(x, nf_PARTIALFACT|nf_ORIG, NULL); }
GEN
factoredpolred2(GEN x, GEN fa) { return Polred(x, nf_PARTIALFACT, fa); }

/********************************************************************/
/**                                                                **/
/**                           POLREDABS                            **/
/**                                                                **/
/********************************************************************/
/* set V[k] := matrix of multiplication by nk.zk[k] */
static GEN
set_mulid(GEN V, GEN M, GEN Mi, long r1, long r2, long N, long k)
{
  GEN v, Mk = cgetg(N+1, t_MAT);
  long i, e;
  for (i = 1; i < k; i++) gel(Mk,i) = gmael(V, i, k);
  for (     ; i <=N; i++)
  {
    v = vecmul(gel(M,k), gel(M,i));
    v = RgM_RgC_mul(Mi, split_realimag(v, r1, r2));
    gel(Mk,i) = grndtoi(v, &e);
    if (e > -5) return NULL;
  }
  gel(V,k) = Mk; return Mk;
}

static GEN
ZM_image_shallow(GEN M, long *pr)
{
  long j, k, r;
  GEN y, d = ZM_pivots(M, &k);
  r = lg(M)-1 - k;
  y = cgetg(r+1,t_MAT);
  for (j=k=1; j<=r; k++)
    if (d[k]) gel(y,j++) = gel(M,k);
  *pr = r; return y;
}

/* U = base change matrix, R = Cholesky form of the quadratic form [matrix
 * Q from algo 2.7.6] */
static GEN
chk_gen_init(FP_chk_fun *chk, GEN R, GEN U)
{
  CG_data *d = (CG_data*)chk->data;
  GEN P, V, D, inv, bound, S, M;
  long N = lg(U)-1, r1 = d->r1, r2 = (N-r1)>>1;
  long i, j, prec, firstprim = 0, skipfirst = 0;
  pari_sp av;

  d->u = U;
  d->ZKembed = M = RgM_mul(d->M, U);

  av = avma; bound = d->bound;
  D = cgetg(N+1, t_VECSMALL);
  for (i = 1; i <= N; i++)
  {
    pari_sp av2 = avma;
    P = get_pol(d, gel(M,i));
    if (!P) pari_err_PREC("chk_gen_init");
    P = gerepilecopy(av2, ZX_radical(P));
    D[i] = degpol(P);
    if (D[i] == N)
    { /* primitive element */
      GEN B = embed_T2(gel(M,i), r1);
      if (!firstprim) firstprim = i; /* index of first primitive element */
      if (DEBUGLEVEL>2) err_printf("chk_gen_init: generator %Ps\n",P);
      if (gcmp(B,bound) < 0) bound = gerepileuptoleaf(av2, B);
    }
    else
    {
      if (DEBUGLEVEL>2) err_printf("chk_gen_init: subfield %Ps\n",P);
      if (firstprim)
      { /* cycle basis vectors so that primitive elements come last */
        GEN u = d->u, e = M;
        GEN te = gel(e,i), tu = gel(u,i), tR = gel(R,i);
        long tS = D[i];
        for (j = i; j > firstprim; j--)
        {
          u[j] = u[j-1];
          e[j] = e[j-1];
          R[j] = R[j-1];
          D[j] = D[j-1];
        }
        gel(u,firstprim) = tu;
        gel(e,firstprim) = te;
        gel(R,firstprim) = tR;
        D[firstprim] = tS; firstprim++;
      }
    }
  }
  if (!firstprim)
  { /* try (a little) to find primitive elements to improve bound */
    GEN x = cgetg(N+1, t_VECSMALL);
    if (DEBUGLEVEL>1)
      err_printf("chk_gen_init: difficult field, trying random elements\n");
    for (i = 0; i < 10; i++)
    {
      GEN e, B;
      for (j = 1; j <= N; j++) x[j] = (long)random_Fl(7) - 3;
      e = RgM_zc_mul(M, x);
      B = embed_T2(e, r1);
      if (gcmp(B,bound) >= 0) continue;
      P = get_pol(d, e); if (!P) pari_err_PREC( "chk_gen_init");
      if (!ZX_is_squarefree(P)) continue;
      if (DEBUGLEVEL>2) err_printf("chk_gen_init: generator %Ps\n",P);
      bound = B ;
    }
  }

  if (firstprim != 1)
  {
    inv = ginv( split_realimag(M, r1, r2) ); /*TODO: use QR?*/
    V = gel(inv,1);
    for (i = 2; i <= r1+r2; i++) V = gadd(V, gel(inv,i));
    /* V corresponds to 1_Z */
    V = grndtoi(V, &j);
    if (j > -5) pari_err_BUG("precision too low in chk_gen_init");
    S = mkmat(V); /* 1 */

    V = cgetg(N+1, t_VEC);
    for (i = 1; i <= N; i++,skipfirst++)
    { /* S = Q-basis of subfield generated by nf.zk[1..i-1] */
      GEN Mx, M2;
      long j, k, h, rkM, dP = D[i];

      if (dP == N) break; /* primitive */
      Mx = set_mulid(V, M, inv, r1, r2, N, i);
      if (!Mx) break; /* prec. problem. Stop */
      if (dP == 1) continue;
      rkM = lg(S)-1;
      M2 = cgetg(N+1, t_MAT); /* we will add to S the elts of M2 */
      gel(M2,1) = col_ei(N, i); /* nf.zk[i] */
      k = 2;
      for (h = 1; h < dP; h++)
      {
        long r; /* add to M2 the elts of S * nf.zk[i]  */
        for (j = 1; j <= rkM; j++) gel(M2,k++) = ZM_ZC_mul(Mx, gel(S,j));
        setlg(M2, k); k = 1;
        S = ZM_image_shallow(shallowconcat(S,M2), &r);
        if (r == rkM) break;
        if (r > rkM)
        {
          rkM = r;
          if (rkM == N) break;
        }
      }
      if (rkM == N) break;
      /* Q(w[1],...,w[i-1]) is a strict subfield of nf */
    }
  }
  /* x_1,...,x_skipfirst generate a strict subfield [unless N=skipfirst=1] */
  chk->skipfirst = skipfirst;
  if (DEBUGLEVEL>2) err_printf("chk_gen_init: skipfirst = %ld\n",skipfirst);

  /* should be DEF + gexpo( max_k C^n_k (bound/k)^(k/2) ) */
  bound = gerepileuptoleaf(av, bound);
  prec = chk_gen_prec(N, (gexpo(bound)*N)/2);
  if (DEBUGLEVEL)
    err_printf("chk_gen_init: new prec = %ld (initially %ld)\n", prec, d->prec);
  if (prec > d->prec) pari_err_BUG("polredabs (precision problem)");
  if (prec < d->prec) d->ZKembed = gprec_w(M, prec);
  return bound;
}

static GEN
polredabs_i(GEN x, nfmaxord_t *S, GEN *u, long flag)
{
  FP_chk_fun chk = { &chk_gen, &chk_gen_init, NULL, NULL, 0 };
  nffp_t F;
  CG_data d;
  GEN v, y, a;
  long i, l;

  nfinit_basic_flag(S, x, flag);
  x = S->T0;
  if (degpol(x) == 1)
  {
    long vx = varn(x);
    *u = NULL;
    return mkvec2(mkvec( pol_x(vx) ),
                  mkvec( deg1pol_shallow(gen_1, negi(gel(S->T,2)), vx) ));
  }
  chk.data = (void*)&d;
  polred_init(S, &F, &d);
  d.bound = embed_T2(F.ro, d.r1);
  if (realprec(d.bound) > F.prec) d.bound = rtor(d.bound, F.prec);
  for (;;)
  {
    GEN R = R_from_QR(F.G, F.prec);
    if (R)
    {
      d.prec = F.prec;
      d.M    = F.M;
      v = fincke_pohst(mkvec(R),NULL,-1, 0, &chk);
      if (v) break;
    }
    F.prec = precdbl(F.prec);
    F.ro = NULL;
    make_M_G(&F, 1);
    if (DEBUGLEVEL) pari_warn(warnprec,"polredabs0",F.prec);
  }
  y = gel(v,1);
  a = gel(v,2); l = lg(a);
  for (i = 1; i < l; i++) /* normalize wrt z -> -z */
    if (ZX_canon_neg(gel(y,i)) && (flag & (nf_ORIG|nf_RAW)))
      gel(a,i) = ZC_neg(gel(a,i));
  *u = d.u; return v;
}

GEN
polredabs0(GEN x, long flag)
{
  pari_sp av = avma;
  GEN Y, A, u, v;
  nfmaxord_t S;
  long i, l;

  v = polredabs_i(x, &S, &u, flag);
  remove_duplicates(v);
  Y = gel(v,1);
  A = gel(v,2);
  l = lg(A); if (l == 1) pari_err_BUG("polredabs (missing vector)");
  if (DEBUGLEVEL) err_printf("Found %ld minimal polynomials.\n",l-1);
  if (!(flag & nf_ALL))
  {
    GEN y = findmindisc(Y);
    for (i = 1; i < l; i++)
      if (ZX_equal(gel(Y,i), y)) break;
    Y = mkvec(gel(Y,i));
    A = mkvec(gel(A,i)); l = 2;
  }
  if (flag & (nf_RAW|nf_ORIG)) for (i = 1; i < l; i++)
  {
    GEN y = gel(Y,i), a = gel(A,i);
    if (u) a = RgV_RgC_mul(S.basis, ZM_ZC_mul(u, a));
    if (flag & nf_ORIG)
    {
      a = QXQ_reverse(a, S.T);
      if (!isint1(S.unscale)) a = gdiv(a, S.unscale); /* not RgX_Rg_div */
      a = mkpolmod(a,y);
    }
    gel(Y,i) = mkvec2(y, a);
  }
  return gerepilecopy(av, (flag & nf_ALL)? Y: gel(Y,1));
}

GEN
polredabsall(GEN x, long flun) { return polredabs0(x, flun | nf_ALL); }
GEN
polredabs(GEN x) { return polredabs0(x,0); }
GEN
polredabs2(GEN x) { return polredabs0(x,nf_ORIG); }

/* relative polredabs/best. Returns relative polynomial by default (flag = 0)
 * flag & nf_ORIG: + element (base change)
 * flag & nf_ABSOLUTE: absolute polynomial */
static GEN
rnfpolred_i(GEN nf, GEN R, long flag, long best)
{
  const char *f = best? "rnfpolredbest": "rnfpolredabs";
  const long abs = ((flag & nf_ORIG) && (flag & nf_ABSOLUTE));
  GEN listP = NULL, red, pol, A, P, T, rnfeq;
  pari_sp av = avma;

  if (typ(R) == t_VEC) {
    if (lg(R) != 3) pari_err_TYPE(f,R);
    listP = gel(R,2);
    R = gel(R,1);
  }
  if (typ(R) != t_POL) pari_err_TYPE(f,R);
  nf = checknf(nf);
  T = nf_get_pol(nf);
  R = RgX_nffix(f, T, R, 0);
  if (best || (flag & nf_PARTIALFACT))
  {
    rnfeq = abs? nf_rnfeq(nf, R): nf_rnfeqsimple(nf, R);
    pol = gel(rnfeq,1);
    if (listP) pol = mkvec2(pol, listP);
    red = best? polredbest_i(pol, abs? 1: 2)
              : polredabs0(pol, (abs? nf_ORIG: nf_RAW)|nf_PARTIALFACT);
    P = gel(red,1);
    A = gel(red,2);
  }
  else
  {
    nfmaxord_t S;
    GEN rnf, u, v, y, a;
    long i, j, l;
    pari_timer ti;
    if (DEBUGLEVEL>1) timer_start(&ti);
    rnf = rnfinit(nf, R);
    rnfeq = rnf_get_map(rnf);
    pol = rnf_zkabs(rnf);
    if (DEBUGLEVEL>1) timer_printf(&ti, "absolute basis");
    v = polredabs_i(pol, &S, &u, nf_ORIG);
    pol = gel(pol,1);
    y = gel(v,1); P = findmindisc(y);
    a = gel(v,2);
    l = lg(y); A = cgetg(l, t_VEC);
    for (i = j = 1; i < l; i++)
      if (ZX_equal(gel(y,i),P))
      {
        GEN t = gel(a,i);
        if (u) t = RgV_RgC_mul(S.basis, ZM_ZC_mul(u,t));
        gel(A,j++) = t;
      }
    setlg(A,j); /* mod(A[i], pol) are all roots of P in Q[X]/(pol) */
  }
  if (DEBUGLEVEL>1) err_printf("reduced absolute generator: %Ps\n",P);
  if (flag & nf_ABSOLUTE)
  {
    if (flag & nf_ORIG)
    {
      GEN a = gel(rnfeq,2); /* Mod(a,pol) root of T */
      GEN k = gel(rnfeq,3); /* Mod(variable(R),R) + k*a root of pol */
      if (typ(A) == t_VEC) A = gel(A,1); /* any root will do */
      a = RgX_RgXQ_eval(a, lift_shallow(A), P); /* Mod(a, P) root of T */
      P = mkvec3(P, mkpolmod(a,P), gsub(A, gmul(k,a)));
    }
    return gerepilecopy(av, P);
  }
  if (typ(A) != t_VEC)
  {
    A = eltabstorel_lift(rnfeq, A);
    P = lift_if_rational( RgXQ_charpoly(A, R, varn(R)) );
  }
  else
  { /* canonical factor */
    long i, l = lg(A), v = varn(R);
    GEN besta = NULL;
    for (i = 1; i < l; i++)
    {
      GEN a = eltabstorel_lift(rnfeq, gel(A,i));
      GEN p = lift_if_rational( RgXQ_charpoly(a, R, v) );
      if (i == 1 || cmp_universal(p, P) < 0) { P = p; besta = a; }
    }
    A = besta;
  }
  if (flag & nf_ORIG) P = mkvec2(P, mkpolmod(RgXQ_reverse(A,R),P));
  return gerepilecopy(av, P);
}
GEN
rnfpolredabs(GEN nf, GEN R, long flag)
{ return rnfpolred_i(nf,R,flag, 0); }
GEN
rnfpolredbest(GEN nf, GEN R, long flag)
{
  if (flag < 0 || flag > 3) pari_err_FLAG("rnfpolredbest");
  return rnfpolred_i(nf,R,flag, 1);
}
