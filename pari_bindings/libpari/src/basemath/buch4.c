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

/*******************************************************************/
/*                        RNFISNORM                                */
/*     (Adapted from Denis Simon's original implementation)        */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

static void
p_append(GEN p, hashtable *H, hashtable *H2)
{
  ulong h = H->hash(p);
  hashentry *e = hash_search2(H, (void*)p, h);
  if (!e)
  {
    hash_insert2(H, (void*)p, NULL, h);
    if (H2) hash_insert2(H2, (void*)p, NULL, h);
  }
}

/* N a t_INT */
static void
Zfa_append(GEN N, hashtable *H, hashtable *H2)
{
  if (!is_pm1(N))
  {
    GEN v = gel(absZ_factor(N),1);
    long i, l = lg(v);
    for (i=1; i<l; i++) p_append(gel(v,i), H, H2);
  }
}
/* N a t_INT or t_FRAC or ideal in HNF*/
static void
fa_append(GEN N, hashtable *H, hashtable *H2)
{
  switch(typ(N))
  {
    case t_INT:
      Zfa_append(N,H,H2);
      break;
    case t_FRAC:
      Zfa_append(gel(N,1),H,H2);
      Zfa_append(gel(N,2),H,H2);
      break;
    default: /*t_MAT*/
      Zfa_append(gcoeff(N,1,1),H,H2);
      break;
  }
}

/* apply lift(rnfeltup) to all coeffs, without rnf structure */
static GEN
nfX_eltup(GEN nf, GEN rnfeq, GEN x)
{
  long i, l;
  GEN y = cgetg_copy(x, &l), zknf = nf_nfzk(nf, rnfeq);
  for (i=2; i<l; i++) gel(y,i) = nfeltup(nf, gel(x,i), zknf);
  y[1] = x[1]; return y;
}

static hashtable *
hash_create_INT(ulong s)
{ return hash_create(s, (ulong(*)(void*))&hash_GEN,
                        (int(*)(void*,void*))&equalii, 1); }
GEN
rnfisnorminit(GEN T, GEN R, int galois)
{
  pari_sp av = avma;
  long i, l, dR;
  GEN S, gen, cyc, bnf, nf, nfabs, rnfeq, bnfabs, k, polabs;
  GEN y = cgetg(9, t_VEC);
  hashtable *H = hash_create_INT(100UL);

  if (galois < 0 || galois > 2) pari_err_FLAG("rnfisnorminit");
  T = get_bnfpol(T, &bnf, &nf);
  if (!bnf) bnf = Buchall(nf? nf: T, nf_FORCE, DEFAULTPREC);
  if (!nf) nf = bnf_get_nf(bnf);

  R = get_bnfpol(R, &bnfabs, &nfabs);
  if (!gequal1(leading_coeff(R))) pari_err_IMPL("non monic relative equation");
  dR = degpol(R);
  if (dR <= 2) galois = 1;

  R = RgX_nffix("rnfisnorminit", T, R, 1);
  if (nf_get_degree(nf) == 1) /* over Q */
    rnfeq = mkvec5(R,gen_0,gen_0,T,R);
  else if (galois == 2) /* needs eltup+abstorel */
    rnfeq = nf_rnfeq(nf, R);
  else /* needs abstorel */
    rnfeq = nf_rnfeqsimple(nf, R);
  polabs = gel(rnfeq,1);
  k = gel(rnfeq,3);
  if (!bnfabs || !gequal0(k))
    bnfabs = Buchall(polabs, nf_FORCE, nf_get_prec(nf));
  if (!nfabs) nfabs = bnf_get_nf(bnfabs);

  if (galois == 2)
  {
    GEN P = polabs==R? leafcopy(R): nfX_eltup(nf, rnfeq, R);
    setvarn(P, fetch_var_higher());
    galois = !!nfroots_if_split(&nfabs, P);
    (void)delete_var();
  }

  cyc = bnf_get_cyc(bnfabs);
  gen = bnf_get_gen(bnfabs); l = lg(cyc);
  for(i=1; i<l; i++)
  {
    GEN g = gel(gen,i);
    if (ugcdiu(gel(cyc,i), dR) == 1) break;
    Zfa_append(gcoeff(g,1,1), H, NULL);
  }
  if (!galois)
  {
    GEN Ndiscrel = diviiexact(nf_get_disc(nfabs), powiu(nf_get_disc(nf), dR));
    Zfa_append(Ndiscrel, H, NULL);
  }
  S = hash_keys(H); settyp(S,t_VEC);
  gel(y,1) = bnf;
  gel(y,2) = bnfabs;
  gel(y,3) = R;
  gel(y,4) = rnfeq;
  gel(y,5) = S;
  gel(y,6) = nf_pV_to_prV(nf, S);
  gel(y,7) = nf_pV_to_prV(nfabs, S);
  gel(y,8) = stoi(galois); return gerepilecopy(av, y);
}

/* T as output by rnfisnorminit
 * if flag=0 assume extension is Galois (==> answer is unconditional)
 * if flag>0 add to S all primes dividing p <= flag
 * if flag<0 add to S all primes dividing abs(flag)

 * answer is a vector v = [a,b] such that
 * x = N(a)*b and x is a norm iff b = 1  [assuming S large enough] */
GEN
rnfisnorm(GEN T, GEN x, long flag)
{
  pari_sp av = avma;
  GEN bnf, rel, R, rnfeq, nfpol;
  GEN nf, aux, H, U, Y, M, A, bnfS, sunitrel, futu, S, S1, S2, Sx;
  long L, i, itu;
  hashtable *H0, *H2;
  if (typ(T) != t_VEC || lg(T) != 9)
    pari_err_TYPE("rnfisnorm [please apply rnfisnorminit()]", T);
  bnf = gel(T,1);
  rel = gel(T,2);
  bnf = checkbnf(bnf);
  rel = checkbnf(rel);
  nf = bnf_get_nf(bnf);
  x = nf_to_scalar_or_alg(nf,x);
  if (gequal0(x)) { set_avma(av); return mkvec2(gen_0, gen_1); }
  if (gequal1(x)) { set_avma(av); return mkvec2(gen_1, gen_1); }
  R = gel(T,3);
  rnfeq = gel(T,4);
  if (gequalm1(x) && odd(degpol(R)))
  { set_avma(av); return mkvec2(gen_m1, gen_1); }

  /* build set T of ideals involved in the solutions */
  nfpol = nf_get_pol(nf);
  S = gel(T,5);
  H0 = hash_create_INT(100UL);
  H2 = hash_create_INT(100UL);
  L = lg(S);
  for (i = 1; i < L; i++) p_append(gel(S,i),H0,NULL);
  S1 = gel(T,6);
  S2 = gel(T,7);
  if (flag > 0)
  {
    forprime_t T;
    ulong p;
    u_forprime_init(&T, 2, flag);
    while ((p = u_forprime_next(&T))) p_append(utoipos(p), H0,H2);
  }
  else if (flag < 0)
    Zfa_append(utoipos(-flag),H0,H2);
  /* overkill: prime ideals dividing x would be enough */
  A = idealnumden(nf, x);
  fa_append(gel(A,1), H0,H2);
  fa_append(gel(A,2), H0,H2);
  Sx = hash_keys(H2); L = lg(Sx);
  if (L > 1)
  { /* new primes */
    settyp(Sx, t_VEC);
    S1 = shallowconcat(S1, nf_pV_to_prV(nf, Sx));
    S2 = shallowconcat(S2, nf_pV_to_prV(rel, Sx));
  }

  /* computation on T-units */
  futu = shallowconcat(bnf_get_fu(rel), bnf_get_tuU(rel));
  bnfS = bnfsunit(bnf,S1,LOWDEFAULTPREC);
  sunitrel = shallowconcat(futu, gel(bnfsunit(rel,S2,LOWDEFAULTPREC), 1));

  A = lift_shallow(bnfissunit(bnf,bnfS,x));
  L = lg(sunitrel);
  itu = lg(nf_get_roots(nf))-1; /* index of torsion unit in bnfsunit(nf) output */
  M = cgetg(L+1,t_MAT);
  for (i=1; i<L; i++)
  {
    GEN u = eltabstorel(rnfeq, gel(sunitrel,i));
    gel(sunitrel,i) = u;
    u = bnfissunit(bnf,bnfS, gnorm(u));
    if (lg(u) == 1) pari_err_BUG("rnfisnorm");
    gel(u,itu) = lift_shallow(gel(u,itu)); /* lift root of 1 part */
    gel(M,i) = u;
  }
  aux = zerocol(lg(A)-1); gel(aux,itu) = utoipos( bnf_get_tuN(rel) );
  gel(M,L) = aux;
  H = ZM_hnfall(M, &U, 2);
  Y = RgM_RgC_mul(U, inverseimage(H,A));
  /* Y: sols of MY = A over Q */
  setlg(Y, L);
  aux = factorback2(sunitrel, gfloor(Y));
  x = mkpolmod(x,nfpol);
  if (!gequal1(aux)) x = gdiv(x, gnorm(aux));
  x = lift_if_rational(x);
  if (typ(aux) == t_POLMOD && degpol(nfpol) == 1)
    gel(aux,2) = lift_if_rational(gel(aux,2));
  return gerepilecopy(av, mkvec2(aux, x));
}

GEN
bnfisnorm(GEN bnf, GEN x, long flag)
{
  pari_sp av = avma;
  GEN T = rnfisnorminit(pol_x(fetch_var()), bnf, flag == 0? 1: 2);
  GEN r = rnfisnorm(T, x, flag == 1? 0: flag);
  (void)delete_var();
  return gerepileupto(av,r);
}
