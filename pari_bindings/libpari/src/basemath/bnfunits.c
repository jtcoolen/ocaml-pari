/* Copyright (C) 2020  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_bnf

/* if x a famat, assume it is a true unit (very costly to check even that
 * it's an algebraic integer) */
GEN
bnfisunit(GEN bnf, GEN x)
{
  long tx = typ(x), i, r1, RU, e, n, prec;
  pari_sp av = avma;
  GEN t, logunit, ex, nf, pi2_sur_w, rx, emb, arg;

  bnf = checkbnf(bnf); nf = bnf_get_nf(bnf);
  RU = lg(bnf_get_logfu(bnf));
  n = bnf_get_tuN(bnf); /* # { roots of 1 } */
  if (tx == t_MAT)
  { /* famat, assumed OK */
    if (lg(x) != 3) pari_err_TYPE("bnfisunit [not a factorization]", x);
  } else {
    x = nf_to_scalar_or_basis(nf,x);
    if (typ(x) != t_COL)
    { /* rational unit ? */
      GEN v;
      long s;
      if (typ(x) != t_INT || !is_pm1(x)) return cgetg(1,t_COL);
      s = signe(x); set_avma(av); v = zerocol(RU);
      gel(v,RU) = utoi((s > 0)? 0: n>>1);
      return v;
    }
    if (!isint1(Q_denom(x))) { set_avma(av); return cgetg(1,t_COL); }
  }

  r1 = nf_get_r1(nf);
  prec = nf_get_prec(nf);
  for (i = 1;; i++)
  {
    GEN rlog;
    nf = bnf_get_nf(bnf);
    logunit = bnf_get_logfu(bnf);
    rlog = real_i(logunit);
    rx = nflogembed(nf,x,&emb, prec);
    if (rx)
    {
      GEN logN = RgV_sum(rx); /* log(Nx), should be ~ 0 */
      if (gexpo(logN) > -20)
      { /* precision problem ? */
        if (typ(logN) != t_REAL) { set_avma(av); return cgetg(1,t_COL); } /*no*/
        if (i == 1 && typ(x) != t_MAT &&
            !is_pm1(nfnorm(nf, x))) { set_avma(av); return cgetg(1, t_COL); }
      }
      else
      {
        ex = RU == 1? cgetg(1,t_COL)
                    : RgM_solve(rlog, rx); /* ~ fundamental units exponents */
        if (ex) { ex = grndtoi(ex, &e); if (e < -4) break; }
      }
    }
    if (i == 1)
      prec = nbits2prec(gexpo(x) + 128);
    else
    {
      if (i > 4) pari_err_PREC("bnfisunit");
      prec = precdbl(prec);
    }
    if (DEBUGLEVEL) pari_warn(warnprec,"bnfisunit",prec);
    bnf = bnfnewprec_shallow(bnf, prec);
  }
  /* choose a large embedding => small relative error */
  for (i = 1; i < RU; i++)
    if (signe(gel(rx,i)) > -1) break;
  if (RU == 1) t = gen_0;
  else
  {
    t = imag_i( row_i(logunit,i, 1,RU-1) );
    t = RgV_dotproduct(t, ex);
    if (i > r1) t = gmul2n(t, -1);
  }
  if (typ(emb) != t_MAT) arg = garg(gel(emb,i), prec);
  else
  {
    GEN p = gel(emb,1), e = gel(emb,2);
    long j, l = lg(p);
    arg = NULL;
    for (j = 1; j < l; j++)
    {
      GEN a = gmul(gel(e,j), garg(gel(gel(p,j),i), prec));
      arg = arg? gadd(arg, a): a;
    }
  }
  t = gsub(arg, t); /* = arg(the missing root of 1) */
  pi2_sur_w = divru(mppi(prec), n>>1); /* 2pi / n */
  e = umodiu(roundr(divrr(t, pi2_sur_w)), n);
  if (n > 2)
  {
    GEN z = algtobasis(nf, bnf_get_tuU(bnf)); /* primitive root of 1 */
    GEN ro = RgV_dotproduct(row(nf_get_M(nf), i), z);
    GEN p2 = roundr(divrr(garg(ro, prec), pi2_sur_w));
    e *= Fl_inv(umodiu(p2,n), n);
    e %= n;
  }
  gel(ex,RU) = utoi(e); setlg(ex, RU+1); return gerepilecopy(av, ex);
}

/* split M a square ZM in HNF as [H, B; 0, Id], H in HNF without 1-eigenvalue */
static GEN
hnfsplit(GEN M, GEN *pB)
{
  long i, l = lg(M);
  for (i = l-1; i; i--)
    if (!equali1(gcoeff(M,i,i))) break;
  if (!i) { *pB = zeromat(0, l-1); return cgetg(1, t_MAT); }
  *pB = matslice(M, 1, i, i+1, l-1); return matslice(M, 1, i, 1, i);
}

/* S a list of prime ideal in idealprimedec format. If pH != NULL, set it to
 * the HNF of the S-class group and return bnfsunit, else return bnfunits */
static GEN
bnfsunit_i(GEN bnf, GEN S, GEN *pH, GEN *pA, GEN *pden)
{
  long FLAG, i, lS = lg(S);
  GEN M, U1, U2, U, V, H, Sunit, B, g;

  if (!is_vec_t(typ(S))) pari_err_TYPE("bnfsunit",S);
  bnf = checkbnf(bnf);
  if (lS == 1)
  {
    *pA = cgetg(1,t_MAT);
    *pden = gen_1; return cgetg(1,t_VEC);
  }
  M = cgetg(lS,t_MAT); /* relation matrix for the S class group */
  g = cgetg(lS,t_MAT); /* principal part */
  FLAG = pH ? 0: nf_GENMAT;
  for (i = 1; i < lS; i++)
  {
    GEN pr = gel(S,i);
    checkprid(pr);
    if (pH)
      gel(M,i) = isprincipal(bnf, pr);
    else
    {
      GEN v = bnfisprincipal0(bnf, pr, FLAG);
      gel(M,i) = gel(v,1);
      gel(g,i) = gel(v,2);
    }
  }
  /* S class group and S units, use ZM_hnflll to get small 'U' */
  M = shallowconcat(M, diagonal_shallow(bnf_get_cyc(bnf)));
  H = ZM_hnflll(M, &U, 1); setlg(U, lS); if (pH) *pH = H;
  U1 = rowslice(U,1, lS-1);
  U2 = rowslice(U,lS, lg(M)-1); /* (M | cyc) [U1; U2] = 0 */
  H = ZM_hnflll(U1, pH? NULL: &V, 0);
 /* U1 = upper left corner of U, invertible. S * U1 = principal ideals
  * whose generators generate the S-units */
  H = hnfsplit(H, &B);
 /*                     [ H B  ]            [ H^-1   - H^-1 B ]
  * U1 * V = HNF(U1) =  [ 0 Id ], inverse = [  0         Id   ]
  * S * HNF(U1) = integral generators for S-units = Sunit */
  Sunit = cgetg(lS, t_VEC);
  if (pH)
  {
    long nH = lg(H) - 1;
    FLAG = nf_FORCE | nf_GEN;
    for (i = 1; i < lS; i++)
    {
      GEN C = NULL, E;
      if (i <= nH) E = gel(H,i); else { C = gel(S,i), E = gel(B,i-nH); }
      gel(Sunit,i) = gel(isprincipalfact(bnf, C, S, E, FLAG), 2);
    }
  }
  else
  {
    GEN cycgen = bnf_build_cycgen(bnf);
    U1 = ZM_mul(U1, V);
    U2 = ZM_mul(U2, V);
    FLAG = nf_FORCE | nf_GENMAT;
    for (i = 1; i < lS; i++)
    {
      GEN a = famatV_factorback(g, gel(U1,i));
      GEN b = famatV_factorback(cycgen, ZC_neg(gel(U2,i)));
      gel(Sunit,i) = famat_reduce(famat_mul(a, b));
    }
  }
  H = ZM_inv(H, pden);
  *pA = shallowconcat(H, ZM_neg(ZM_mul(H,B))); /* top inverse * den */
  return Sunit;
}
GEN
bnfsunit(GEN bnf,GEN S,long prec)
{
  pari_sp av = avma;
  long i, l = lg(S);
  GEN v, R, nf, A, den, U, cl, H = NULL;
  bnf = checkbnf(bnf); nf = bnf_get_nf(bnf);
  v = cgetg(7, t_VEC);
  gel(v,1) = U = bnfsunit_i(bnf, S, &H, &A, &den);
  gel(v,2) = mkvec2(A, den);
  gel(v,3) = cgetg(1,t_VEC); /* dummy */
  R = bnf_get_reg(bnf);
  cl = bnf_get_clgp(bnf);
  if (l != 1)
  {
    GEN u,A, G = bnf_get_gen(bnf), D = ZM_snf_group(H,NULL,&u), h = ZV_prod(D);
    long lD = lg(D);
    A = cgetg(lD, t_VEC);
    for(i = 1; i < lD; i++) gel(A,i) = idealfactorback(nf, G, gel(u,i), 1);
    cl = mkvec3(h, D, A);
    R = mpmul(R, h);
    for (i = 1; i < l; i++)
    {
      GEN pr = gel(S,i), p = pr_get_p(pr);
      long f = pr_get_f(pr);
      R = mpmul(R, logr_abs(itor(p,prec)));
      if (f != 1) R = mulru(R, f);
      gel(U,i) = nf_to_scalar_or_alg(nf, gel(U,i));
    }
  }
  gel(v,4) = R;
  gel(v,5) = cl;
  gel(v,6) = S; return gerepilecopy(av, v);
}
GEN
bnfunits(GEN bnf, GEN S)
{
  pari_sp av = avma;
  GEN A, den, U, fu, tu;
  bnf = checkbnf(bnf);
  U = bnfsunit_i(bnf, S? S: cgetg(1,t_VEC), NULL, &A, &den);
  if (!S) S = cgetg(1,t_VEC);
  fu = bnf_compactfu(bnf);
  if (!fu)
  {
    long i, l;
    fu = bnf_has_fu(bnf); if (!fu) bnf_build_units(bnf);
    fu = shallowcopy(fu); l = lg(fu);
    for (i = 1; i < l; i++) gel(fu,i) = to_famat_shallow(gel(fu,i),gen_1);
  }
  tu = nf_to_scalar_or_basis(bnf_get_nf(bnf), bnf_get_tuU(bnf));
  U = shallowconcat(U, vec_append(fu, to_famat_shallow(tu,gen_1)));
  return gerepilecopy(av, mkvec4(U, S, A, den));
}
GEN
sunits_mod_units(GEN bnf, GEN S)
{
  pari_sp av = avma;
  GEN A, den;
  bnf = checkbnf(bnf);
  return gerepilecopy(av, bnfsunit_i(bnf, S, NULL, &A, &den));
}

/* v_S(x), x in famat form */
static GEN
sunit_famat_val(GEN nf, GEN S, GEN x)
{
  long i, l = lg(S);
  GEN v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(v,i) = famat_nfvalrem(nf, x, gel(S,i), NULL);
  return v;
}
/* v_S(x), x algebraic number */
static GEN
sunit_val(GEN nf, GEN S, GEN x, GEN N)
{
  long i, l = lg(S);
  GEN v = zero_zv(l-1), N0 = N;
  for (i = 1; i < l; i++)
  {
    GEN P = gel(S,i), p = pr_get_p(P);
    if (dvdii(N, p)) { v[i] = nfval(nf,x,P); (void)Z_pvalrem(N0, p, &N0); }
  }
  return is_pm1(N0)? v: NULL;
}

/* if *px a famat, assume it's an S-unit */
static GEN
make_unit(GEN nf, GEN U, GEN *px)
{
  GEN den, v, w, A, gen = gel(U,1), S = gel(U,2), x = *px;
  long cH, i, l = lg(S);

  if (l == 1) return cgetg(1, t_COL);
  A = gel(U,3); den = gel(U,4);
  cH = nbrows(A);
  if (typ(x) == t_MAT && lg(x) == 3)
  {
    w = sunit_famat_val(nf, S, x); /* x = S v */
    v = ZM_ZC_mul(A, w);
    w += cH; w[0] = evaltyp(t_COL) | evallg(lg(A) - cH);
  }
  else
  {
    GEN N;
    x = nf_to_scalar_or_basis(nf,x);
    switch(typ(x))
    {
      case t_INT:  N = x; if (!signe(N)) return NULL; break;
      case t_FRAC: N = mulii(gel(x,1),gel(x,2)); break;
      default: { GEN d = Q_denom(x); N = mulii(idealnorm(nf,gmul(x,d)), d); }
    }
    /* relevant primes divide N */
    if (is_pm1(N)) return zerocol(l-1);
    w = sunit_val(nf, S, x, N);
    if (!w) return NULL;
    v = ZM_zc_mul(A, w);
    w += cH; w[0] = evaltyp(t_VECSMALL) | evallg(lg(A) - cH);
    w = zc_to_ZC(w);
  }
  if (!is_pm1(den)) for (i = 1; i <= cH; i++)
  {
    GEN r;
    gel(v,i) = dvmdii(gel(v,i), den, &r);
    if (r != gen_0) return NULL;
  }
  v = shallowconcat(v, w); /* append bottom of w (= [0 Id] part) */
  for (i = 1; i < l; i++)
  {
    GEN e = gel(v,i);
    if (signe(e)) x = famat_mulpow_shallow(x, gel(gen,i), negi(e));
  }
  *px = x; return v;
}

static GEN
bnfissunit_i(GEN bnf, GEN x, GEN U)
{
  GEN w, nf, v = NULL;
  bnf = checkbnf(bnf);
  nf = bnf_get_nf(bnf);
  if ( (w = make_unit(nf, U, &x)) ) v = bnfisunit(bnf, x);
  if (!v || lg(v) == 1) return NULL;
  return mkvec2(v, w);
}
static int
checkU(GEN U)
{
  if (typ(U) != t_VEC || lg(U) != 5) return 0;
  return typ(gel(U,1)) == t_VEC && is_vec_t(typ(gel(U,2)))
         && typ(gel(U,4))==t_INT;
}
GEN
bnfisunit0(GEN bnf, GEN x, GEN U)
{
  pari_sp av = avma;
  GEN z;
  if (!U) return bnfisunit(bnf, x);
  if (!checkU(U)) pari_err_TYPE("bnfisunit",U);
  z = bnfissunit_i(bnf, x, U);
  if (!z) { set_avma(av); return cgetg(1,t_COL); }
  return gerepilecopy(av, shallowconcat(gel(z,2), gel(z,1)));
}

/* OBSOLETE */
static int
checkbnfS_i(GEN v)
{
  GEN S, g, w;
  if (typ(v) != t_VEC || lg(v) != 7) return 0;
  g = gel(v,1); w = gel(v,2); S = gel(v,6);
  if (typ(g) != t_VEC || !is_vec_t(typ(S)) || lg(S) != lg(g)) return 0;
  return typ(w) == t_VEC && lg(w) == 3;
}
/* OBSOLETE */
GEN
bnfissunit(GEN bnf, GEN bnfS, GEN x)
{
  pari_sp av = avma;
  GEN z, U;
  if (!checkbnfS_i(bnfS)) pari_err_TYPE("bnfissunit",bnfS);
  U = mkvec4(gel(bnfS,1), gel(bnfS,6), gmael(bnfS,2,1), gmael(bnfS,2,2));
  z = bnfissunit_i(bnf, x, U);
  if (!z) { set_avma(av); return cgetg(1,t_COL); }
  return gerepilecopy(av, shallowconcat(gel(z,1), gel(z,2)));
}
