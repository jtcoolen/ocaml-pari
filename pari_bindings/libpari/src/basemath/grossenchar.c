/* Copyright (C) 2000  The PARI group.

   This file is part of the PARI/GP package.

   PARI/GP is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation. It is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY WHATSOEVER.

   Check the License for details. You should have received a copy of it, along
   with the package; see the file 'COPYING'. If not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_gchar

static GEN gchari_eval(GEN gc, GEN chi, GEN x, long flag, GEN logchi, GEN s0, long prec);

/*********************************************************************/
/**                                                                 **/
/**                    HECKE GROSSENCHARACTERS                      **/
/**         contributed by Pascal Molin and Aurel Page (2022)       **/
/**                                                                 **/
/*********************************************************************/

/*
  characters can be represented by:
   - t_COL of coordinates on the snf basis (mostly for gp use): prefix gchar_
   - t_VEC of coordinates on the internal basis: prefix gchari_
   - t_VEC of R/Z components (logs): prefix gcharlog_

   see gchar_internal for SNF -> internal
   and gchari_duallog for internal -> R/Z components
*/

/*
localstar: represent (Z_F/m)^* . {+-1}^r + generators of U_{i-1}(p)/U_i
structure:
- [ sprk for p^k | m ] , size np
- [ Ufil_p for p | m ], size np
- m_oo, t_VECSMALL of size nci <= r1 (indices of real places)

where Ufil_p = [ Mat([gen[j], t_COL of size ncp]_j) ]_{1<=i<=k}
*/

#define GC_LENGTH   12
#define LOCS_LENGTH 4

static GEN
compute_Lcyc(GEN Lsprk, GEN moo)
{
  long i, l = lg(Lsprk), len = l+lg(moo)-1;
  GEN Lcyc = cgetg(len,t_VEC);
  for (i = 1; i < l; i++)   gel(Lcyc,i) = sprk_get_cyc(gel(Lsprk,i));
  for (     ; i < len; i++) gel(Lcyc,i) = mkvec(gen_2);
  return Lcyc;
}

static long
sprk_get_ncp(GEN sprk) { return lg(sprk_get_cyc(sprk)) - 1; }

/* true nf; modulus = [ factor(m_f), m_oo ] */
static GEN
localstar(GEN nf, GEN famod, GEN moo)
{
  GEN Lcyc, cyc, Lsprk, Lgenfil, P = gel(famod,1), E = gel(famod,2);
  long i, l = lg(P);

  Lsprk = cgetg(l, t_VEC);
  Lgenfil = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    long e, k = itos(gel(E,i));
    GEN v, sprk = log_prk_init(nf, gel(P,i), k, NULL);
    gel(Lsprk,i) = sprk;
    gel(Lgenfil,i) = v = cgetg(k+1, t_VEC);
    /* log on sprk of generators of U_{e-1}/U_e(pr) */
    gel(v, 1) = col_ei(sprk_get_ncp(sprk), 1);
    for (e = 2; e <= k; e++) gel(v, e) = sprk_log_gen_pr2(nf, sprk, e);
  }
  Lcyc = compute_Lcyc(Lsprk, moo);
  if (lg(Lcyc) > 1)
    cyc = shallowconcat1(Lcyc);
  else
    cyc = cgetg(1, t_VEC);
  return mkvec4(cyc, Lsprk, Lgenfil, mkvec2(famod,moo));
}

/* (nv * log|x^sigma|/norm, arg(x^sigma))/2*Pi
 * substract norm so that we project to the hyperplane
 * H : sum n_s x_s = 0 */
static GEN
nfembedlog(GEN *pnf, GEN x, long prec)
{
  pari_sp av = avma;
  GEN logs, cxlogs, nf = *pnf;
  long k, r1, r2, n, extrabit, extranfbit = 0, nfprec, nfprec0, logprec;

  nfprec0 = nf_get_prec(nf);
  nf_get_sign(nf, &r1, &r2);
  n = r1 + 2*r2;
  logprec = prec;
  extrabit = expu(n) + gexpo(nf_get_M(nf)) + 10;
  if (typ(x) == t_MAT)
  {
    long l = lg(gel(x,2));
    if (l > 1)
    {
      extrabit += expu(l) + gexpo(gel(x,2));
      extranfbit = gexpo(gel(x,1));
    }
  }
  else
  {
    x = nf_to_scalar_or_basis(nf,x);
    extranfbit = gexpo(x);
  }
  if (DEBUGLEVEL>3)
    err_printf("  nfembedlog: prec=%d extrabit=%d nfprec=%d extralogprec=%d\n",
               prec, nbits2extraprec(extrabit + extranfbit), nfprec0,
               nbits2extraprec(extrabit));
  nfprec = prec + nbits2extraprec(extrabit + extranfbit);
  logprec = prec + nbits2extraprec(extrabit);
  if (nfprec0 < nfprec)
  {
    if (DEBUGLEVEL)
      err_printf("  nfembedlog: increasing prec %d -> %d\n", nfprec0, nfprec);
    *pnf = nf = nfnewprec_shallow(nf, nfprec);
    av = avma;
  }
  if (!(cxlogs = nf_cxlog(nf, x, logprec))) return gc_NULL(av);
  if (!(cxlogs = nf_cxlog_normalize(nf, cxlogs, logprec))) return gc_NULL(av);
  logs = cgetg(n+1,t_COL);
  for (k = 1; k <= r1+r2; k++) gel(logs,k) = real_i(gel(cxlogs,k));
  for (     ; k <= n; k++) gel(logs,k) = gmul2n(imag_i(gel(cxlogs,k-r2)), -1);
  extrabit = gexpo(logs);
  if (extrabit < 0) extrabit = 0;
  prec += nbits2extraprec(extrabit);
  return gerepileupto(av, gdiv(logs, Pi2n(1,prec)));
}

static GEN
gchar_Sval(GEN nf, GEN S, GEN x)
{
  GEN res = cgetg(lg(S),t_COL);
  long i;
  if (typ(x)==t_MAT)
    for (i=1; i<lg(S); i++)
      gel(res, i) = famat_nfvalrem(nf, x, gel(S,i), NULL);
  else
    for (i=1; i<lg(S); i++)
      gel(res, i) = stoi(nfval(nf, x, gel(S,i)));
  return res;
}

/* true nf; log_prk(x*pi_pr^{-v_pr(x)}), sign(sigma(x)) */
static GEN
gchar_logm(GEN nf, GEN locs, GEN x)
{
  GEN moo, loga, Lsprk = locs_get_Lsprk(locs);
  long i, l = lg(Lsprk);

  moo = locs_get_m_infty(locs);
  if (typ(x) != t_MAT) x = to_famat_shallow(x, gen_1);
  loga = cgetg(l+1, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN sprk = gel(Lsprk, i), pr = sprk_get_pr(sprk);
    GEN g = vec_append(gel(x,1), pr_get_gen(pr));
    GEN v = famat_nfvalrem(nf, x, pr, NULL);
    GEN e = vec_append(gel(x,2), gneg(v));
    gel(loga, i) = famat_zlog_pr(nf, g, e, sprk, NULL);
  }
  gel(loga, i) = zc_to_ZC( nfsign_arch(nf, x, moo) );
  return shallowconcat1(loga);
}

static GEN
gchar_nflog(GEN *pnf, GEN zm, GEN S, GEN x, long prec)
{
  GEN emb = nfembedlog(pnf, x, prec);
  if (!emb) return NULL;
  return shallowconcat1(mkvec3(gchar_Sval(*pnf,S,x),
                               gchar_logm(*pnf,zm,x), emb));
}

/*******************************************************************/
/*                                                                 */
/*                        CHARACTER GROUP                          */
/*                                                                 */
/*******************************************************************/

/* embed v in vector of length size, shifted to the right */
static GEN
embedcol(GEN v, long size, long shift)
{
  long k;
  GEN w = zerocol(size);
  for (k = 1; k < lg(v); k++) gel(w, shift+k) = gel(v,k);
  return w;
}

/* write exact rationals from real approximation, in place. */
static void
shallow_clean_rat(GEN v, long k0, long k1, GEN den, long prec)
{
  long k, e, bit = -prec2nbits(prec)/2;
  for (k = k0; k <= k1; k++)
  {
    GEN t = gel(v,k);
    if (den) t = gmul(t, den);
    t = grndtoi(t, &e);
    if (DEBUGLEVEL>1) err_printf("[%Ps*%Ps=%Ps..e=%ld|prec=%ld]\n",
                                 gel(v,k), den? den: gen_1, t, e, prec);
    if (e > bit)
      pari_err_BUG("gcharinit, non rational entry"); /*LCOV_EXCL_LINE*/
    gel(v, k) = den? gdiv(t, den): t;
  }
}

/* FIXME: export ? */
static GEN
rowreverse(GEN m)
{
  long i, l;
  GEN v;
  if (lg(m) == 1) return m;
  l = lgcols(m); v = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++) v[i] = l - i;
  return rowpermute(m, v);
}

/* lower-left hnf on subblock m[r0+1..r0+nr, c0+1..c0+nc]
 * return base change matrix U */
static GEN
hnf_block(GEN m, long r0, long nr, long c0, long nc)
{
  GEN block, u, uu;
  long nm = lg(m)-1, k;
  pari_sp av = avma;

  block = matslice(m, r0+1, r0+nr, c0+1, c0+nc);
  block = Q_remove_denom(block, NULL);
  block = rowreverse(block);

  (void)ZM_hnfall(block, &u, 1);
  vecreverse_inplace(u); uu = matid(nm); /* embed in matid */
  for (k = 1; k <= nc; k++) gel(uu,c0+k) = embedcol(gel(u,k),nm,c0);
  return gerepilecopy(av, uu);
}

/* (matrix, starting row, nb rows, starting column, nb columns) */
static GEN
lll_block(GEN m, long r0, long nr, long c0, long nc)
{
  GEN block, u, uu;
  long nm = lg(m)-1, k;
  pari_sp av = avma;

  block = matslice(m, r0+1, r0+nr, c0+1, c0+nc);
  u = lll(block); if (lg(u) <= nc) return NULL;
  vecreverse_inplace(u); uu = matid(nm); /* embed in matid */
  for (k = 1; k <= nc; k++) gel(uu,c0+k) = embedcol(gel(u,k),nm,c0);
  return gerepilecopy(av, uu);
}

/* insert a column at index i, no copy */
static GEN
shallowmatinsert(GEN m, GEN x, long i)
{
  long k, n = lg(m);
  GEN mm = cgetg(n+1,t_MAT);
  for (k=1; k < i; k++) gel(mm, k) = gel(m, k);
  gel(mm, i) = x;
  for (k=i; k < n; k++) gel(mm, k+1) = gel(m, k);
  return mm;
}

static GEN
vec_v0(long n, long n0, long r1, long r2)
{
  long k;
  GEN C = zerocol(n);
  for (k = 1; k <= r1; k++) gel(C, n0++) = gen_1;
  for (k = 1; k <= r2; k++) gel(C, n0++) = gen_2;
  return C;
}

/* select cm embeddings; return a matrix */
/* TODO detect if precision was insufficient */
static GEN
cm_select(GEN nf, GEN cm, long prec)
{
  GEN emb, keys, v, m_sel, imag_emb;
  long nalg, d_cm, r_cm, c, i, j, r2 = nf_get_r2(nf);
  pari_sp av;

  d_cm = degpol(gel(cm, 1)); /* degree of the cm field; even */
  nalg = d_cm / 2; /* nb of clusters */
  r_cm = nf_get_degree(nf) / d_cm; /* nb by cluster; nalg * r_cm = r2 */
  m_sel = zeromatcopy(nalg, r2); /* selection matrix */
  av = avma;
  /* group complex embeddings */
  emb = nfeltembed(nf, gel(cm, 2), NULL, prec);
  /* sort */
  imag_emb = imag_i(emb);
  keys = gadd(gmul(mppi(prec), real_i(emb)), gabs(imag_emb, prec));
  v = indexsort(keys);

  for (j = c = 1; c <= nalg; c++)
  {
    int ref = gsigne(gel(imag_emb, v[j]));
    gcoeff(m_sel, c, v[j]) = gen_1;
    j++;
    for (i = 2; i <= r_cm; i++)
    {
      int s = gsigne(gel(imag_emb, v[j]));
      gcoeff(m_sel, c, v[j]) = (s == ref) ? gen_1 : gen_m1;
      j++;
    }
  }
  return gc_const(av, m_sel);
}

static GEN gchar_hnfreduce_shallow(GEN gc, GEN cm);
static void gchar_snfbasis_shallow(GEN gc, GEN rel);
static void gcharmat_tinverse(GEN gc, GEN m, long prec);
static GEN gcharmatnewprec_shallow(GEN gc, long mprec);

/* return a set S of prime ideals such that Cl_S(K) = 1 */
static GEN
bestS(GEN bnf)
{
  GEN v, S, hw, hv = bnf_get_no(bnf), DL, dl;
  long i, lS;
  ulong l;
  forprime_t P;

  if (equali1(hv)) return mkvec2(cgetg(1,t_VEC), cgetg(1,t_VEC));
  v = diagonal_shallow(bnf_get_cyc(bnf));
  S = cgetg(expi(hv)+2, t_VEC); lS = 1;
  DL = cgetg(expi(hv)+2, t_VEC);
  u_forprime_init(&P,2,ULONG_MAX);
  while ((l = u_forprime_next(&P)))
  {
    pari_sp av = avma;
    GEN w, Sl = idealprimedec(bnf, utoi(l));
    long nSl = lg(Sl)-1;
    for (i = 1; i < nSl; i++) /* remove one prime ideal */
    {
      dl = isprincipal(bnf, gel(Sl,i));
      w = ZM_hnf(shallowconcat(v, dl));
      hw = ZM_det(w);
      if (cmpii(hw, hv) < 0)
      {
        gel(DL,lS) = dl;
        gel(S,lS++) = gel(Sl,i);
        hv = hw; v = w; av = avma;
        if (equali1(hv)) { setlg(S, lS); setlg(DL, lS); return mkvec2(S,DL); }
      }
    }
    set_avma(av);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN
gcharDLdata(GEN bnf, GEN S, GEN DL)
{
  GEN M, h, Minv, Lalpha, t, dl, alpha, gen, cyc = bnf_get_cyc(bnf);
  long i;
  M = shallowmatconcat(DL);
  h = bnf_get_no(bnf);
  gen = bnf_get_gen(bnf);

  /* compute right inverse of M modulo cyc */
  M = shallowtrans(M);
  M = shallowmatconcat(mkcol2(M,diagonal_shallow(cyc)));
  Minv = matinvmod(M,h);
  Minv = vecslice(Minv,1,lg(Minv)-lg(cyc));
  Minv = shallowtrans(Minv);

  Lalpha = cgetg(lg(Minv),t_VEC);
  for (i=1; i<lg(Minv); i++)
  {
    /* gen[i] = (alpha) * prod_j S[j]^Minv[j,i] */
    t = isprincipalfact(bnf, gel(gen,i), S, gneg(gel(Minv,i)), nf_GENMAT);
    dl = gel(t, 1); alpha = gel(t, 2);
    if (!gequal0(dl)) pari_err_BUG("gcharDLdata (non-principal ideal)");
    gel(Lalpha,i) = alpha;
  }
  return mkvec2(Minv, Lalpha);
}

/* compute basis of characters; gc[1] generating family, as rows */
GEN
gcharinit(GEN bnf, GEN mod, long prec)
{
  pari_sp av = avma;
  GEN nf, zm, zmcyc, S, DLdata, sfu, logx;
  GEN fa2, archp, z, C, gc, cm, cyc, rel, U, Ui, m, m_inv, m0, u0;
  long n, k, r1, r2, ns, nc, nu, nm, order;
  long evalprec = prec, nfprec, mprec, embprec;

  prec = evalprec + EXTRAPREC64; /* default 64 extra bits */

  /* note on precision:

     - evalprec = precision requested for evaluation

     - prec = precision available
            = (relative) precision of parameters = m_inv
            = (relative) precision of evaluation of small chars
              if no cancellation

     - nfprec = internal nf precision used for
       the embedding matrix m

     In the structure we store [evalprec,prec,nfprec]

     When evaluating chi(x) at evalprec,
     we need prec >= max(evalprec + exponent(chi), nfprec+exponent(x))
     where exponent(x) is the exponent of the number field element alpha
     obtained after principalisation of x.

     If prec is not sufficient, we call gcharnewprec.

     To mitigate precision increase, we always initialize the structure
     with 64 extra bits.

     Final remark: a gchar struct may have inconsistent values
     for prec and nfprec, for example if it has been saved and loaded at
     default prec : one should call gcharnewprec then.
  */

  if (!checkbnf_i(bnf))
  {
    nfprec = prec;
    bnf = bnfinit0(bnf, 1, NULL, nfprec);
    nf = shallowcopy(bnf_get_nf(bnf));
  }
  else
  {
    GEN fu = bnf_get_sunits(bnf);
    if (!fu) fu = bnf_get_fu(bnf); /* impose fundamental units */
    nf = shallowcopy(bnf_get_nf(bnf));
    nfprec = nf_get_prec(nf);
  }

  /* Dirichlet group + make sure mod contains archimedean places */
  mod = check_mod_factored(nf,mod,NULL,&fa2,&archp,NULL);
  sort_factor(fa2, (void*)&cmp_prime_ideal, &cmp_nodata);
  zm = localstar(nf, fa2, archp);
  zmcyc = locs_get_cyc(zm);

  /* set of primes S and valuations of generators */
  S = bestS(bnf);
  DLdata = gel(S,2);
  S = gel(S,1);
  DLdata = gcharDLdata(bnf, S, DLdata);

  nf_get_sign(nf, &r1, &r2);
  n = r1+2*r2;
  ns = lg(S) - 1;
  nu = r1+r2-1 + ns;
  nc = lg(zmcyc) - 1;
  nm = ns+nc+n; /* number of parameters = ns + nc + r1 + r2 + r2 */

  /* units and S-units */
  sfu = gel(bnfunits(bnf,S), 1);
  sfu = vec_shorten(sfu, nu); /* remove torsion */

  /* root of unity */
  order = bnf_get_tuN(bnf);
  z = bnf_get_tuU(bnf);

  /* maximal cm subfield */
  cm = nfsubfieldscm(nf, 0);

  /*
   Now compute matrix of parameters,
   until we obtain the right precision
   FIXME: make sure generators, units, and discrete logs
          do not depend on precision.

   m0 is the matrix of units embeddings
   u  is the HNF base change, m = m0*u

   subsequent steps may lead to precision increase, we put everything in gc
   struct and modify it in place.

     A) sets m0
     B) sets U, cyc, rel, U and Ui
     C) sets m_inv
  */

  /* A) make big matrix m0 of embeddings of units */

  if (DEBUGLEVEL>2) err_printf("start matrix m\n");
  m = cgetg(nm + 1, t_MAT);
  mprec = nbits2prec(nm+10) + EXTRAPREC64;
  embprec = mprec;
  for(;;)
  {
    for (k = 1; k <= nu; k++)
    { /* Lambda_S (S-units) then Lambda_f, fund. units */
      logx = gchar_nflog(&nf, zm, S, gel(sfu,k), embprec);
      if (!logx) break;
      gel(m, k) = logx;
    }
    if (k > nu) break;
    if (DEBUGLEVEL) err_printf("gcharinit: increasing embprec %d -> %d\n",
                               embprec, precdbl(embprec));
    embprec = precdbl(embprec);
  }
  for (k = 1; k <= nc; k++) /* Gamma, structure of (Z/m)* */
  {
    C = zerocol(nm);
    gel(C, ns+k) = gel(zmcyc, k);
    gel(m, nu+k) = C;
  }
  /* zeta, root of unity */
  gel(m, nu+nc+1) = gchar_nflog(&nf, zm, S, z, mprec);
  shallow_clean_rat(gel(m, nu+nc+1), 1, nm, stoi(order), mprec);
  for (k = 1; k <= r2; k++) /* embed Z^r_2 */
  {
    C = zerocol(nm);
    gel(C, ns+nc+r1+r2+k) = gen_1;
    gel(m, nu+nc+1+k) = C;
  }
  if (DEBUGLEVEL>1) err_printf("matrix m = %Ps\n", m);

  m0 = m;
  u0 = gen_0;
  rel = U = Ui = gen_0;
  cyc = gen_0;
  m_inv = gen_0;

  /* only m and m_inv depend on prec and are recomputed under gcharnewprec. */
  gc = mkvecn(GC_LENGTH,
              m_inv, /* internal basis, characters as rows */
              bnf,
              nf,
              zm,    /* Zk/mod, nc components */
              S,     /* generators of clgp, ns components */
              DLdata,
              sfu,
              mkvec2(mkvecsmall3(evalprec,prec,nfprec),
                     mkvecsmall3(0,0,0)), /* ntors, nfree, nalg */
              cyc, /* reduced components */
              mkvec3(rel, U, Ui), /* internal / SNF base change */
              m0,                 /* embeddings of units */
              u0);                /* m_inv = (m0 u0)~^-1 */

  /* B) do HNF reductions + LLL (may increase precision) */
  m = gchar_hnfreduce_shallow(gc, cm);

  /* C) compute snf basis of torsion subgroup */
  rel = shallowtrans(matslice(m, 1, ns+nc, 1, ns+nc));
  gchar_snfbasis_shallow(gc, rel);

  /* D) transpose inverse m_inv = (m0*u)~^-1 (may increase precision) */
  gcharmat_tinverse(gc, m, prec);
  return gerepilecopy(av, gc);
}

/* b) do HNF reductions + LLL, keep base change u0 */
static GEN
gchar_hnfreduce_shallow(GEN gc, GEN cm)
{
  GEN bnf = gchar_get_bnf(gc), nf = gchar_get_nf(gc), u, u0, m;
  long order, r1, r2, ns, nc, n, nu, nm, nalg = 0, mprec;

  nf_get_sign(nf, &r1, &r2);
  n = r1 + 2*r2;
  nu = r1 + r2 - 1;
  ns = gchar_get_ns(gc);
  nc = gchar_get_nc(gc);
  nm = ns+nc+n; /* ns + nc + r1 + r2 + r2 */
  order = 2*bnf_get_tuN(bnf);
  u0 = matid(nm);
  m = shallowcopy(gchar_get_m0(gc)); /* keep m0 unchanged */
  mprec = gprecision(m);
  if (DEBUGLEVEL>1) err_printf("matrix m = %Ps\n", m);
  if (nc)
  { /* keep steps 1&2 to make sure we have zeta_m */
    u = hnf_block(m, ns,nc, ns+nu,nc+1);
    u0 = ZM_mul(u0, u); m = RgM_ZM_mul(m, u);
    if (DEBUGLEVEL>2) err_printf("step 1 -> %Ps\n", m);
    u = hnf_block(m, ns,nc, ns,nu+nc);
    u0 = ZM_mul(u0, u); m = RgM_ZM_mul(m, u);
    if (DEBUGLEVEL>2) err_printf("step 2 -> %Ps\n", m);
  }
  if (r2)
  {
    u = hnf_block(m, nm-r2,r2, nm-r2-1,r2+1);
    u0 = ZM_mul(u0, u); m = RgM_ZM_mul(m, u);
    if (DEBUGLEVEL>2) err_printf("step 3 -> %Ps\n", m);
  }
  /* remove last column */
  setlg(u0, nm); setlg(m, nm);
  if (DEBUGLEVEL>2) err_printf("remove last col -> %Ps\n", m);

  if (!gequal0(cm))
  {
    GEN v, Nargs;
    long bit = - prec2nbits(mprec) + 16 + expu(order);
    /* reduce on Norm arguments */
    v = cm_select(nf, cm, gchar_get_nfprec(gc));
    if (DEBUGLEVEL>2) err_printf("cm_select -> %Ps\n", v);
    nalg = nbrows(v);
    gchar_set_u0(gc, u0);
    for(;;)
    {
      long e, emax, i;
      Nargs = gmul(v, rowslice(m, nm-r2+1, nm));
      if (DEBUGLEVEL>2) err_printf("Nargs -> %Ps\n", Nargs);
      emax = bit-1;
      for (i = ns+nc+1; i < lg(Nargs); i++)
      {
        gel(Nargs,i) = grndtoi(gmulgs(gel(Nargs,i), order), &e);
        emax = maxss(emax,e);
      }
      if (emax < bit) break;
      if (DEBUGLEVEL>1) err_printf("cm select: doubling prec\n");
      mprec = precdbl(mprec); m = gcharmatnewprec_shallow(gc, mprec);
    }
    if (DEBUGLEVEL>2) err_printf("rounded Nargs -> %Ps\n", Nargs);
    u = hnf_block(Nargs, 0, nalg, ns+nc, n-1);
    u0 = ZM_mul(u0, u); m = RgM_ZM_mul(m, u);
    if (DEBUGLEVEL>2) err_printf("after cm reduction -> %Ps\n", m);
  }

  /* apply LLL on Lambda_m, may need to increase prec */
  gchar_set_nalg(gc, nalg);
  gchar_set_u0(gc, u0);

  /* TODO factor these two LLL reduction codes in a function? */
  if (nalg > 0)
  {
    GEN u = NULL;
    while (1)
    {
      if ((u = lll_block(m, ns+nc, n, ns+nc, nalg))) break;
      mprec = precdbl(mprec); m = gcharmatnewprec_shallow(gc, mprec);
    }
    u0 = ZM_mul(u0, u); m = RgM_ZM_mul(m, u);
    if (DEBUGLEVEL>1) err_printf("after LLL reduction (CM block) -> %Ps\n", m);
  }
  gchar_set_u0(gc, u0);

  if (nu > 0)
  {
    GEN u = NULL;
    while (1)
    {
      if ((u = lll_block(m, ns+nc, n, ns+nc+nalg, n-1-nalg))) break;
      mprec = precdbl(mprec); m = gcharmatnewprec_shallow(gc, mprec);
    }
    u0 = ZM_mul(u0, u); m = RgM_ZM_mul(m, u);
    if (DEBUGLEVEL>1) err_printf("after LLL reduction (trans block) -> %Ps\n", m);
  }
  gchar_set_u0(gc, u0); return m;
}

/* convert to snf basis of torsion + Z^(r1+2*r2-1) */
static void
gchar_snfbasis_shallow(GEN gc, GEN rel)
{
  long n, r1, r2, lU, lUi;
  GEN nf, cyc, U, Ui;

  nf = gchar_get_nf(gc);
  nf_get_sign(nf, &r1, &r2);
  n = r1+2*r2;

  rel = ZM_hnf(rel);
  if (DEBUGLEVEL>1) err_printf("relations after hnf: %Ps\n", rel);
  cyc = ZM_snf_group(rel, &U, &Ui);
  if (lg(cyc)==1)
  {
    cyc = zerovec(n-1);
    U = shallowconcat(zeromat(n-1,lg(rel)-1),matid(n-1));
    Ui = shallowtrans(U);
  }
  else if (n!=1)
  {
    cyc = shallowconcat(cyc, zerovec(n-1));
    U = shallowmatconcat(mkmat22(U,gen_0,gen_0,matid(n-1)));
    Ui = shallowmatconcat(mkmat22(Ui,gen_0,gen_0,matid(n-1)));
  }

  rel = shallowtrans(rel);
  lU = lg(U);
  lUi = lg(Ui);
  U = shallowtrans(U);
  if (lg(U)==1 && lUi>1) U = zeromat(0,lUi-1);
  Ui = shallowtrans(Ui);
  if (lg(Ui)==1 && lU>1) Ui = zeromat(0,lU-1);

  gchar_set_nfree(gc, n-1);
  gchar_set_ntors(gc, (lg(cyc)-1) - (n-1));
  gchar_set_cyc(gc, shallowconcat(cyc, real_0(gchar_get_prec(gc))));
  gchar_set_HUUi(gc, rel, U, Ui);
}

static long
mextraprec(GEN m_inv, GEN m, GEN gc)
{
  return nbits2extraprec(2*maxss(gexpo(m_inv),1) + expu(lg(m))
          + gexpo(gchar_get_u0(gc)) + 10);
}

/* c) transpose inverse + clean rationals.
 * prec = target prec for m^-1,
 * mprec = prec of m */
static void
gcharmat_tinverse(GEN gc, GEN m, long prec)
{
  GEN m_inv;
  long k, r1, r2, ns, nc, nalg, nm, mprec, bitprec = prec2nbits(prec);

  nf_get_sign(gchar_get_nf(gc), &r1, &r2);
  ns = gchar_get_ns(gc);
  nc = gchar_get_nc(gc);
  nalg = gchar_get_nalg(gc);
  nm = ns + nc + r1 + 2*r2;
  if (lg(m)==1) { gchar_set_basis(gc,zeromat(0,nm)); return; }

  mprec = gprecision(m); /* possibly 0, if m is exact */
  while (1)
  {
    long targetmprec = 0;
    GEN v0, mm;
    /* insert at column ns+nc+r1+r2, or last column if cm */
    v0 = vec_v0(nm, ns+nc+1, r1, r2);
    mm = shallowmatinsert(m, v0, nalg? nm: nm-r2);
    if (DEBUGLEVEL>1) err_printf("add v0 -> %Ps\n", mm);
    mm = shallowtrans(mm);
    m_inv = RgM_inv(mm); /* invert matrix, may need to increase prec */
    if (m_inv)
    {
      if (DEBUGLEVEL>1) err_printf("inverse: %Ps\n",m_inv);
      m_inv = vecsplice(m_inv, nalg? nm: nm-r2); /* remove v0 */
      if (DEBUGLEVEL>1) err_printf("v0 removed: %Ps\n", m_inv);
      m_inv = shallowtrans(m_inv);
      if (!mprec) break;
      /* enough precision? */
      /* |B - A^(-1)| << |B|.|Id-B*A| */
      if (gexpo(m_inv) + gexpo(gsub(RgM_mul(m_inv, m), gen_1)) + expu(lg(m))
          <= -bitprec)
      {
        /* |A^(-1) - (A+H)^(-1)| << |H|.|A^(-1)|^2 */
        targetmprec = prec + mextraprec(m_inv,m,gc);
        if (mprec >= targetmprec) break;
      }
      else targetmprec = 0;
    }
    mprec = maxss(precdbl(mprec), targetmprec);
    if (mprec < DEFAULTPREC) mprec = DEFAULTPREC;
    m = gcharmatnewprec_shallow(gc, mprec); /* m0 * u0 to higher prec */
  }
  /* clean rationals */
  if (nc)
  { /* reduce mod exponent of the group, TODO reduce on each component */
    GEN zmcyc = locs_get_cyc(gchar_get_zm(gc));
    GEN e = ZV_lcm(zmcyc);
    for (k = 1; k <= nc; k++)
      shallow_clean_rat(gel(m_inv, ns+k), 1, nm - 1, /*zmcyc[k]*/e, prec);
  }
  if (r2)
  {
    for (k = 1; k <= r2; k++)
      shallow_clean_rat(gel(m_inv,nm-k+1), 1, nm-1, NULL, prec);
  }
  if (nalg)
  {
    long i, j, e;
    for (i = 1; i<=r1+r2; i++)
      for (j = 1; j <= nalg; j++)
      {
        e = gexpo(gcoeff(m_inv, ns+nc+j, ns+nc+i));
        if (e > -20) /* TODO is this bound too permissive? */
          pari_err_BUG("gcharinit (nonzero entry)"); /*LCOV_EXCL_LINE*/
        gcoeff(m_inv, ns+nc+j, ns+nc+i) = gen_0;
      }
  }
  if (DEBUGLEVEL>1) err_printf("cyc/cm cleaned: %Ps", shallowtrans(m_inv));
  /* normalize characters, parameters mod Z */
  for (k = 1; k <= ns+nc; k++) gel(m_inv, k) = gfrac(gel(m_inv, k));
  /* increase relative prec of real values */
  gchar_set_basis(gc, gprec_w(m_inv, prec));
}

/* make sure same argument determination is chosen */
static void
same_arg(GEN v1, GEN v2, long s1, long s2)
{
  long i1, i2, l = lg(v1);
  for (i1 = s1, i2 = s2; i1 < l; i1++,i2++)
  {
    GEN d = grndtoi(gsub(gel(v2,i2),gel(v1,i1)), NULL);
    if (signe(d)) gel(v1,i1) = gadd(gel(v1,i1), d);
  }
}

static void
vaffect_shallow(GEN x, long i0, GEN y)
{
  long i;
  for (i = 1; i < lg(y); i++)
    gel(x, i+i0) = gel(y, i);
}

/* recompute matrix with precision increased */
/* u0 the base change, returns m0 * u0 */
/* mprec: requested precision for m0 */
static GEN
gcharmatnewprec_shallow(GEN gc, long mprec)
{
  GEN nf, m0, u0, sfu;
  long k, l, ns, nc, r1, r2, nfprec, embprec;

  ns = gchar_get_ns(gc);
  nc = gchar_get_nc(gc);
  nf = gchar_get_nf(gc);
  sfu = gchar_get_sfu(gc);
  nf_get_sign(nf, &r1, &r2);
  nfprec = nf_get_prec(gchar_get_nf(gc));

  m0 = gchar_get_m0(gc);
  u0 = gchar_get_u0(gc);

  if (DEBUGLEVEL>3) err_printf("gcharmatnewprec_shallow mprec=%d nfprec=%d\n", mprec, nfprec);

  embprec = mprec; l = lg(sfu);
  while(1)
  { /* recompute the nfembedlogs of s-units and fundamental units */
    for (k = 1; k < l; k++) /* ns s-units, then r1+r2-1 fundamental units */
    {
      GEN emb = nfembedlog(&nf, gel(sfu,k), embprec);
      if (!emb) break;
      same_arg(emb, gel(m0,k),  r1+r2, ns+nc+r1+r2);
      vaffect_shallow(gel(m0, k), ns+nc, emb);
    }
    if (k == l) break;
    if (DEBUGLEVEL)
      err_printf("gcharmatnewprec_shallow: increasing embprec %d -> %d\n",
                 embprec, precdbl(embprec));
    embprec = precdbl(embprec);
  }
  gchar_set_nf(gc, nf);
  gchar_set_nfprec(gc, nfprec);
  return RgM_ZM_mul(m0, u0);
}

static void _check_gchar_group(GEN gc, long flag);
static void
check_gchar_group(GEN gc) { _check_gchar_group(gc, 0); }

/* increase prec if needed. newprec: requested precision for m_inv */
static GEN
gcharnewprec_i(GEN gc, long newprec)
{
  long prec, prec0, nfprec, nfprec0, mprec;
  GEN gc2 = shallowcopy(gc);

  _check_gchar_group(gc2, 1); /* ignore illegal prec */
  prec = gchar_get_prec(gc2);
  nfprec = gchar_get_nfprec(gc2);

  if (newprec > prec)
  { /* increase precision */
    if (DEBUGLEVEL) pari_warn(warnprec,"gcharnewprec",newprec);
    nfprec += newprec - prec;
    prec = newprec;
    gchar_copy_precs(gc, gc2);
    gchar_set_prec(gc2, prec);
    gchar_set_nfprec(gc2, nfprec);
  }

  nfprec0 = nf_get_prec(gchar_get_nf(gc2));
  if (nfprec0 && nfprec > nfprec0)
  {
    if (DEBUGLEVEL) pari_warn(warnprec,"gcharnewprec (nf)", nfprec);
    gchar_set_nf(gc2, nfnewprec_shallow(gchar_get_nf(gc2), nfprec));
  }

  prec0 = gprecision(gchar_get_basis(gc2));
  if (prec0 && prec > prec0)
  {
    GEN cyc, m, m0 = gchar_get_m0(gc);
    if (DEBUGLEVEL) pari_warn(warnprec,"gcharnewprec (minv)", prec);
    gchar_set_m0(gc2, shallowcopy(m0));
    mprec = prec + mextraprec(gchar_get_basis(gc), m0, gc);
    m = gcharmatnewprec_shallow(gc2, mprec);
    if (DEBUGLEVEL>2) err_printf("m0*u0 recomputed -> %Ps\n", m);
    gcharmat_tinverse(gc2, m, prec);
    cyc = shallowcopy(gchar_get_cyc(gc2));
    gel(cyc, lg(cyc)-1) = real_0(prec);
    gchar_set_cyc(gc2, cyc);
  }
  return gc2;
}

/* newprec: requested evalprec */
GEN
gcharnewprec(GEN gc, long newprec)
{
  pari_sp av = avma;
  GEN gc2 = gcharnewprec_i(gc, newprec + EXTRAPREC64);
  gchar_set_evalprec(gc2, newprec);
  return gerepilecopy(av, gc2);
}

static void
check_localstar(GEN x)
{
  if (typ(x) != t_VEC || lg(x) != LOCS_LENGTH + 1)
    pari_err_TYPE("char group", x);
  /* FIXME: check components once settled */
}

int
is_gchar_group(GEN gc)
{
  return (typ(gc) == t_VEC
      &&  lg(gc) == GC_LENGTH + 1
      &&  typ(gel(gc, 8)) == t_VEC
      &&  lg(gel(gc, 8)) == 3
      &&  typ(gmael(gc, 8, 1))  == t_VECSMALL
      &&  typ(gmael(gc, 8, 2))  == t_VECSMALL
      &&  (checkbnf_i(gchar_get_bnf(gc)) != NULL)
      &&  (checknf_i(gchar_get_nf(gc)) != NULL));
}

/* validates the structure format.
 * Raise error if inconsistent precision, unless flag=1. */
static void
_check_gchar_group(GEN gc, long flag)
{
  GEN b, bnf, nf;
  if (typ(gc) != t_VEC || lg(gc) != GC_LENGTH + 1)
    pari_err_TYPE("char group", gc);
  check_localstar(gchar_get_zm(gc));
  if (typ(gchar_get_loccyc(gc)) != t_VEC)
    pari_err_TYPE("gchar group (loccyc)", gc);
  b = gchar_get_basis(gc);
  if (typ(b) != t_MAT) pari_err_TYPE("gchar group (basis)", gc);
  bnf = gchar_get_bnf(gc); checkbnf(bnf);
  nf = gchar_get_nf(gc); checknf(nf);
  if (!gequal(nf_get_pol(nf),nf_get_pol(bnf_get_nf(bnf))))
    pari_err_TYPE("gchar group (nf != bnf.nf)", gc);
  if (typ(gel(gc,8)) != t_VEC ||typ(gmael(gc,8,1)) != t_VECSMALL)
    pari_err_TYPE("gchar group (gc[8])", gc);
  if (!flag)
  {
    long prec0 = gprecision(b), nfprec0 = nf_get_prec(nf);
    if ((prec0 && gchar_get_prec(gc) > prec0)
        || (nfprec0 && gchar_get_nfprec(gc) > nfprec0))
      pari_err_PREC("gchar group, please call gcharnewprec");
  }
}

/* basis of algebraic characters + cyc vector */
static GEN
gchar_algebraic_basis(GEN gc)
{
  long ntors, nfree, nc, nm, r2, nalg, n0, k;
  GEN basis, args, m, w, normchar, alg_basis, tors_basis;

  /* in snf basis */
  ntors = gchar_get_ntors(gc); /* number of torsion generators */
  nfree = gchar_get_nfree(gc);
  nc = ntors + nfree;
  /* in internal basis */
  n0 = gchar_get_ns(gc) + gchar_get_nc(gc); /* last index of torsion char */
  r2 = gchar_get_r2(gc);
  nm = gchar_get_nm(gc);
  /* in both */
  nalg = gchar_get_nalg(gc); /* number of generators of free algebraic chars */

  /* finite order characters have weight 0 */
  tors_basis = matid(ntors);

  /* the norm is an algebraic character */
  normchar = gneg(col_ei(nc+1,nc+1));

  /* add sublattice of algebraic */

  if (!nalg)
  {
    if (DEBUGLEVEL>2) err_printf("nalg=0\n");
    return shallowmatconcat(mkvec2(tors_basis,normchar));
  }

  /* block of k_s parameters of free algebraic */
  args = matslice(gchar_get_basis(gc),n0+1,n0+nalg,nm-r2+1,nm);

  if (r2 == 1)
  {
    /* no parity condition */
    if (DEBUGLEVEL>2) err_printf("r2 = 1 -> args = %Ps\n", args);
    alg_basis = matid(nalg);
    w = gel(args,1);
  }
  else
  {
    /* parity condition: w + k_s = 0 mod 2 for all s,
       ie solve x.K constant mod 2, ie solve
       x.K' = 0 mod 2 where K' = [ C-C0 ] (substract first column)
     */
    /* select block k_s in char parameters and */
    if (DEBUGLEVEL>2) err_printf("block ks -> %Ps\n", args);
    m = cgetg(r2, t_MAT);
    for (k = 1; k < r2; k++)
      gel(m,k) = gsub(gel(args,k+1),gel(args,1));
    if (DEBUGLEVEL>2) err_printf("block ks' -> %Ps", m);
    alg_basis = shallowtrans(gel(matsolvemod(shallowtrans(m),gen_2,gen_0,1),2));
    if (DEBUGLEVEL>2) err_printf("alg_basis -> %Ps\n", alg_basis);
    w = gmul(alg_basis, gel(args,1));
    if (DEBUGLEVEL>2) err_printf("w -> %Ps\n", w);
  }
  /* add weight to infinite order characters, at position nc+1 */
  w = gneg(gdivgs(gmodgs(w, 2), 2));
  alg_basis = shallowmatconcat(mkvec3(alg_basis, zerovec(nfree-nalg), w));

  if (lg(tors_basis)>1)
    basis = shallowmatconcat(mkmat22(tors_basis, gen_0, gen_0, alg_basis));
  else
    basis = alg_basis;
  return shallowmatconcat(mkvec2(shallowtrans(basis),normchar));
}
static GEN
gchar_algebraicnormtype(GEN gc, GEN type)
{
  GEN w = NULL, w1, v;
  long i, r1, r2, n;
  r1 = gchar_get_r1(gc);
  r2 = gchar_get_r2(gc);
  for (i=1; i<=r1; i++)
  {
    w1 = addii(gmael(type,i,1),gmael(type,i,2));
    if (!w) w = w1;
    else if (!equalii(w,w1)) return NULL;
  }
  for (i=r1+1; i<=r1+r2; i++)
  {
    w1 = gmael(type,i,1);
    if (!w) w = w1;
    else if (!equalii(w,w1)) return NULL;
    if (!equalii(w,gmael(type,i,2))) return NULL;
  }
  n = lg(gchar_get_cyc(gc))-1;
  v = zerocol(n);
  gel(v,n) = negi(w);
  return mkvec(v);
}
static GEN
gchar_algebraicoftype(GEN gc, GEN type)
{
  long i, r1, r2, nalg, n0, nm;
  GEN p, q, w, k, matk, chi;
  /* in internal basis */
  n0 = gchar_get_ns(gc) + gchar_get_nc(gc); /* last index of torsion chars */
  r1 = gchar_get_r1(gc);
  r2 = gchar_get_r2(gc);
  nm = gchar_get_nm(gc);
  /* in both */
  nalg = gchar_get_nalg(gc); /* number of generators of free algebraic chars */

  if (typ(type)!=t_VEC || lg(type) != r1+r2+1)
    pari_err_TYPE("gcharalgebraic (1)", type);
  for (i = 1; i<=r1+r2; i++)
    if (typ(gel(type,i)) != t_VEC
        ||lg(gel(type,i)) != 3
        ||typ(gmael(type,i,1)) != t_INT
        ||typ(gmael(type,i,2)) != t_INT)
      pari_err_TYPE("gcharalgebraic (2)", type);

  chi = gchar_algebraicnormtype(gc, type);
  if (chi) return chi;

  if (!nalg) return NULL;
  k = cgetg(r2+1,t_VEC);
  p = gmael(type, 1, 1);
  q = gmael(type, 1, 2); w = addii(p, q);
  gel(k, 1) = subii(q, p);
  for (i = 2; i <= r2; i++)
  {
    p = gmael(type, i, 1);
    q = gmael(type, i, 2);
    if (!equalii(w, addii(p, q))) return NULL;
    gel(k, i) = subii(q, p);
  }
  /* block of k_s parameters of free algebraic */
  matk = matslice(gchar_get_basis(gc),n0+1,n0+nalg,nm-r2+1,nm);
  chi = inverseimage(shallowtrans(matk),shallowtrans(k));
  if (lg(chi) == 1) return NULL;
  for (i=1; i<lg(chi); i++) if (typ(gel(chi,i)) != t_INT) return NULL;
  chi = mkvec4(zerocol(gchar_get_ntors(gc)), chi,
               zerocol(gchar_get_nfree(gc) - nalg), gmul2n(negi(w),-1));
  return mkvec(shallowconcat1(chi));
}

GEN
gcharalgebraic(GEN gc, GEN type)
{
  pari_sp av = avma;
  GEN b;
  check_gchar_group(gc);
  b = type? gchar_algebraicoftype(gc, type): gchar_algebraic_basis(gc);
  if (!b) { set_avma(av); return cgetg(1, t_VEC); }
  return gerepilecopy(av, b);
}

/*********************************************************************/
/*                                                                   */
/*                          CHARACTERS                               */
/*                                                                   */
/*********************************************************************/
/* return chi without (optional) norm component; set *s = to the latter  */
static GEN
check_gchar_i(GEN chi, long l, GEN *s)
{
  long i, n = lg(chi)-1;
  if (n == l)
  { /* norm component */
    if (s)
    {
      *s = gel(chi,l);
      switch(typ(*s))
      {
        case t_INT:
        case t_FRAC:
        case t_REAL:
        case t_COMPLEX: break;
        default: pari_err_TYPE("check_gchar_i [s]", *s);
      }
    }
    chi = vecslice(chi, 1, l-1);
  }
  else
  { /* missing norm component */
    if (n != l-1) pari_err_DIM("check_gchar_i [chi]");
    if (s) *s = gen_0;
  }
  for (i = 1; i < l; i++) if (typ(gel(chi,i))!=t_INT)
    pari_err_TYPE("check_gchar_i [coefficient]", gel(chi,i));
  return chi;
}

static GEN
check_gchar(GEN gc, GEN chi, GEN *s)
{
  if (typ(chi)!=t_COL) pari_err_TYPE("check_gchar [chi]", chi);
  return check_gchar_i(chi, lg(gchar_get_cyc(gc))-1, s);
}

static long
safelgcols(GEN m)
{
  return lg(m)==1 ? 1 : lg(gel(m,1));
}

static GEN
check_gchari(GEN gc, GEN chi, GEN *s)
{
  if (typ(chi)!=t_VEC) pari_err_TYPE("check_gchari [chi]", chi);
  return check_gchar_i(chi, safelgcols(gchar_get_basis(gc)), s);
}

/* from coordinates on snf basis, return coordinates on internal basis, set
 * s to the norm component */
static GEN
gchar_internal(GEN gc, GEN chi, GEN *s)
{
  chi = check_gchar(gc, chi, s);
  return ZV_ZM_mul(chi, gchar_get_Ui(gc));
}

static GEN
RgV_frac_inplace(GEN v, long n)
{
  long i;
  for (i = 1; i <= n; i++) gel(v,i) = gfrac(gel(v,i));
  return v;
}
/* from internal basis form, return the R/Z components */
static GEN
gchari_duallog(GEN gc, GEN chi)
{
  chi = RgV_RgM_mul(chi, gchar_get_basis(gc)); /* take exponents mod Z */
  return RgV_frac_inplace(chi, gchar_get_ns(gc) + gchar_get_nc(gc));
}

/* chip has ncp = #zm[1][i].cyc components */
static GEN
conductor_expo_pr(GEN gens_fil, GEN chip)
{
  long i;
  for (i = lg(gens_fil) - 1; i > 0; i--)
  {
    GEN gens = gel(gens_fil, i);
    long j;
    for (j = 1; j < lg(gens); j++)
    {
      GEN v = gmul(chip, gel(gens, j));
      if (denom_i(v) != gen_1) return stoi(i);
    }
  }
  return gen_0;
}

/* assume chi in log form */
static GEN
gcharlog_conductor_f(GEN gc, GEN chi, GEN *faN)
{
  GEN zm, P, E, Lsprk, ufil;
  long i, l, ic;

  if (gchar_get_nc(gc) == 0)
  {
    if (faN) *faN = trivial_fact();
    return gen_1;
  }
  zm = gchar_get_zm(gc);
  Lsprk = locs_get_Lsprk(zm);
  ufil = locs_get_Lgenfil(zm);
  P = gel(locs_get_famod(zm), 1);
  l = lg(Lsprk); E = cgetg(l, t_VEC);
  for (i = 1, ic = gchar_get_ns(gc); i < l ; i++)
  {
    GEN gens = gel(ufil, i), chip;
    long ncp = sprk_get_ncp(gel(Lsprk,i));

    chip = vecslice(chi, ic + 1, ic + ncp);
    gel(E, i) = conductor_expo_pr(gens, chip);
    ic += ncp;
  }
  if (faN) *faN = famat_remove_trivial(mkmat2(P,E));
  return idealfactorback(gchar_get_nf(gc), P, E, 0); /* red = 0 */
}

/* ={sigma} s.t. k_sigma = 1 */
static GEN
gcharlog_conductor_oo(GEN gc, GEN chi)
{
  long noo, i, n0 = gchar_get_ns(gc) + gchar_get_nc(gc);
  GEN k_real, chi_oo, moo, den;

  moo = locs_get_m_infty(gchar_get_zm(gc));
  noo = lg(moo)-1;
  k_real = vecslice(chi, n0-noo+1, n0);
  chi_oo = zerovec(gchar_get_r1(gc));
  for (i=1; i<=noo; i++)
  {
    den = Q_denom(gel(k_real,i));
    if (den && !equali1(den)) gel(chi_oo, moo[i]) = gen_1;
  }
  return chi_oo;
}

static GEN
gchari_conductor(GEN gc, GEN chi)
{
  chi = gchari_duallog(gc, chi);
  return mkvec2(gcharlog_conductor_f(gc, chi, NULL),
                gcharlog_conductor_oo(gc, chi));
}
GEN
gchar_conductor(GEN gc, GEN chi)
{
  pari_sp av = avma;
  check_gchar_group(gc);
  return gerepilecopy(av, gchari_conductor(gc, gchar_internal(gc, chi, NULL)));
}

int
gcharisalgebraic(GEN gc, GEN chi, GEN *pq)
{
  long i, ntors, nfree, n0, nalg, r1, r2;
  GEN w, chii, v;
  pari_sp av = avma;

  check_gchar_group(gc);
  /* in snf basis */
  ntors = gchar_get_ntors(gc);
  nfree = gchar_get_nfree(gc);
  /* in internal basis */
  r1 = gchar_get_r1(gc);
  r2 = gchar_get_r2(gc);
  n0 = gchar_get_nm(gc) - r2; /* last index before k_s */
  /* in both */
  nalg = gchar_get_nalg(gc); /* number of generators of free algebraic chars */

  chii = gchar_internal(gc, chi, &w);
  /* check component are on algebraic generators */
  for (i = ntors+nalg+1; i <= ntors+nfree; i++)
    if (!gequal0(gel(chi,i))) return gc_bool(av, 0);
  chii = gchari_duallog(gc, chii);

  if (r1) /* not totally complex: finite order * integral power of norm */
  {
    if (typ(w) != t_INT) return gc_bool(av, 0);
    w = negi(w);
    if (pq)
    { /* set the infinity type */
      v = cgetg(r1+r2+1, t_VEC);
      for (i = 1; i <= r1; i++) gel(v, i) = mkvec2(w, gen_0);
      for (  ; i <= r1+r2; i++) gel(v, i) = mkvec2(w, w);
    }
  }
  else /* totally complex */
  {
    int w2;
    w = gneg(gmul2n(w, 1));
    if (typ(w) != t_INT) return gc_bool(av, 0);
    w2 = mpodd(w);
    for (i = 1; i <= r2; i++) /* need k_s + w = 0 mod 2 for all s */
      if (mpodd(gel(chii, n0 + i)) != w2) return gc_bool(av, 0);
    if (pq)
    { /* set the infinity type */
      v = cgetg(r2+1, t_VEC);
      for (i = 1; i <= r2; i++)
      {
        GEN p = gmul2n(subii(w, gel(chii, n0+i)), -1);
        gel(v, i) = mkvec2(p, subii(w, p));
      }
    }
  }
  if (pq) { *pq = gerepilecopy(av, v); av = avma; }
  return gc_bool(av, 1);
}

GEN
gcharlocal(GEN gc, GEN chi, GEN v, long prec, GEN* pbid)
{
  pari_sp av = avma;
  GEN nf = gchar_get_nf(gc), s, chiv, logchi;

  check_gchar_group(gc);
  chi = gchar_internal(gc, chi, &s);
  logchi = gchari_duallog(gc, chi);
  if (typ(v) == t_INT) /* v infinite */
  {
    long i, r1, r2, tau = itos(v), n0 = gchar_get_ns(gc) + gchar_get_nc(gc);
    GEN phi, k;
    nf_get_sign(nf, &r1, &r2);
    if (tau <= 0)
      pari_err_DOMAIN("gcharlocal [index of an infinite place]", "v", "<=", gen_0, v);
    if (tau > r1+r2)
      pari_err_DOMAIN("gcharlocal [index of an infinite place]", "v", ">", stoi(r1+r2), v);
    if (r1+r2 == 1) phi = gen_0;
    else            phi = gel(logchi, n0 + tau);
    if (tau <= r1) /* v real */
    {
      GEN moo = gel(gchar_get_mod(gc),2);
      i = zv_search(moo, tau);
      if (i==0) k = gen_0;
      else
      {
        k = gel(logchi, n0 - lg(moo) + 1 + i); /* 0 or 1/2 */
        if (!gequal0(k)) k = gen_1;
      }
    }
    else /* v complex */
      k = gel(logchi, n0 + r2 + tau);
    if (s) phi = gsub(phi, mulcxI(s));
    chiv = mkvec2(k, phi);
  }
  else /* v finite */
  {
    GEN P = gchar_get_modP(gc);
    long iv;
    checkprid(v);
    iv = gen_search(P, v, (void*)cmp_prime_ideal, cmp_nodata);
    chiv = gchari_eval(gc, NULL, v, 0, logchi, s, prec);
    if (iv <= 0) chiv = mkvec(chiv);
    else
    {
      GEN cyc, w, Lsprk, bid, chip = NULL;
      long i, ic, l = lg(P);
      Lsprk = locs_get_Lsprk(gchar_get_zm(gc));
      for (i = 1, ic = gchar_get_ns(gc); i < l; i++)
      {
        long ncp = sprk_get_ncp(gel(Lsprk,i));
        if (i == iv) { chip = vecslice(logchi, ic + 1, ic + ncp); break; }
        ic += ncp;
      }
      if (!chip) pari_err_BUG("gcharlocal (chip not found)");
      /* TODO store bid instead of recomputing? */
      bid = sprk_to_bid(nf, gel(Lsprk,i), nf_INIT);
      cyc = bid_get_cyc(bid);
      w = RgV_RgM_mul(chip, gel(bid_get_U(bid),1));
      for (i = 1; i < lg(w); i++)
        gel(w,i) = modii(gmul(gel(w,i), gel(cyc,i)), gel(cyc,i));
      chiv = vec_append(w, chiv);
      if (pbid) { *pbid = bid; gerepileall(av, 2, &chiv, pbid); return chiv; }
    }
  }
  return gerepilecopy(av, chiv);
}


/*******************************************************************/
/*                                                                 */
/*                EVALUATION OF CHARACTERS                         */
/*                                                                 */
/*******************************************************************/
GEN
gcharduallog(GEN gc, GEN chi)
{
  pari_sp av = avma;
  GEN logchi, s;
  check_gchar_group(gc);
  logchi = gchari_duallog(gc, gchar_internal(gc, chi, &s));
  return gerepilecopy(av, shallowconcat1(mkcol2(logchi,s)));
}

static GEN
gcharisprincipal(GEN gc, GEN x, GEN *palpha)
{
  GEN t, v, bnf = gchar_get_bnf(gc), DLdata = gchar_get_DLdata(gc);
  t = bnfisprincipal0(bnf, x, nf_GENMAT | nf_FORCE); v = gel(t, 1);
  *palpha = famat_reduce(famat_mul(nffactorback(bnf, gel(DLdata,2), v),
                                   gel(t, 2)));
  return ZM_ZC_mul(gel(DLdata,1), v);
}

/* complete log of ideal; if logchi != NULL make sure precision is
 * sufficient to evaluate gcharlog_eval_raw to attain 'prec' */
static GEN
gchar_log(GEN gc, GEN x, GEN logchi, long prec)
{
  GEN zm, v, alpha, arch_log = NULL, zm_log, nf;

  nf = gchar_get_nf(gc);
  zm = gchar_get_zm(gc);
  v = gcharisprincipal(gc, x, &alpha); /* alpha a GENMAT */
  if (DEBUGLEVEL>2) err_printf("v %Ps\n", v);
  zm_log = gchar_logm(nf, zm, alpha);
  if (DEBUGLEVEL>2) err_printf("zm_log(alpha) %Ps\n", zm_log);

  if (logchi)
  { /* check if precision is sufficient, take care of gexpo = -infty */
    long e, bit = expu(nf_get_degree(nf) + lg(zm_log)-1) + 3;
    e = gexpo(logchi); if (e > 0) bit += e;
    e = gexpo(gel(alpha,1)); if (e > 0) bit += e;
    prec += nbits2extraprec(bit);
  }
  for(;;)
  {
    arch_log = nfembedlog(&nf, alpha, prec);
    if (arch_log) break;
    prec = precdbl(prec);
  }
  if (DEBUGLEVEL>2) err_printf("arch log %Ps\n", arch_log);
  return shallowconcat1(mkvec3(v, gneg(zm_log), gneg(arch_log)));
}

/* gp version, with norm component */
GEN
gcharlog(GEN gc, GEN x, long prec)
{
  pari_sp av = avma;
  GEN norm;

  check_gchar_group(gc);
  norm = idealnorm(gchar_get_nf(gc), x);
  norm = mkcomplex(gen_0, gdiv(glog(norm,prec), Pi2n(1,prec)));
  return gerepilecopy(av, vec_append(gchar_log(gc, x, NULL, prec), norm));
}

static GEN
gcharlog_eval_raw(GEN logchi, GEN logx)
{ GEN v = RgV_dotproduct(logchi, logx); return gsub(v, ground(v)); }

/* if flag = 1 -> value in C, flag = 0 -> value in R/Z
 * chi in chari format without norm component (given in s)
 * if logchi is provided, assume it has enough precision
 * TODO check that logchi argument is indeed used correctly by callers */
static GEN
gchari_eval(GEN gc, GEN chi, GEN x, long flag, GEN logchi, GEN s, long prec)
{
  GEN z, norm, logx;
  if (!logchi)
  {
    long e, precgc = prec, prec0 = gchar_get_prec(gc);
    logchi = gchari_duallog(gc, chi);
    e = gexpo(logchi); if (e > 0) precgc += nbits2extraprec(e);
    if (precgc > prec0)
    {
      gc = gcharnewprec(gc, precgc);
      logchi = gchari_duallog(gc, chi);
    }
  }
  logx = gchar_log(gc, x, logchi, prec);
  norm = gequal0(s)? NULL: idealnorm(gchar_get_nf(gc), x);
  z = gcharlog_eval_raw(logchi, logx);
  if (flag)
  {
    z = expIPiC(gmul2n(z, 1), prec);
    if (norm) z = gmul(z, gpow(norm, gneg(s), prec));
  }
  else if (norm)
    z = gadd(z, mulcxI(gdiv(gmul(s, glog(norm,prec)), Pi2n(1,prec))));
  if (DEBUGLEVEL>1) err_printf("char value %Ps\n", z);
  return z;
}

GEN
gchareval(GEN gc, GEN chi, GEN x, long flag)
{
  GEN s;
  long prec;
  pari_sp av = avma;
  check_gchar_group(gc);
  prec = gchar_get_evalprec(gc);
  chi = gchar_internal(gc, chi, &s);
  return gerepileupto(av, gchari_eval(gc, chi, x, flag, NULL, s, prec));
}

/*******************************************************************/
/*                                                                 */
/*                         IDENTIFICATION                          */
/*                                                                 */
/*******************************************************************/
static GEN
col_2ei(long n, long i) { GEN e = zerocol(n); gel(e,i) = gen_2; return e; }
static GEN
gchar_identify_init(GEN gc, GEN Lv, long prec)
{
  GEN M, cyc, mult, Lpr, Lk1, Lphi1, Lk2, Llog, eps, U, P, nf, moo, vlogchi;
  long beps, r1, r2, d, nmiss, n, npr, nk1, nchi, s, i, j, l, dim, n0, ncol;

  check_gchar_group(gc);
  n0 = gchar_get_ns(gc) + gchar_get_nc(gc); /* last index of torsion char */
  cyc = gchar_get_cyc(gc);
  P = gchar_get_modP(gc);
  moo = gel(gchar_get_mod(gc), 2);
  nf = gchar_get_nf(gc); nf_get_sign(nf, &r1, &r2);
  nchi = lg(cyc)-2; /* ignore norm */
  d = r1 + 2*r2;
  mult = (nchi >= d)? gel(cyc,1): gen_1;
  s = (8*prec2nbits(prec))/10; mult = shifti(mult, s);
  beps = -(7*s) / 16; eps = real2n(beps, prec);
  l = lg(Lv);
  if (lg(gen_sort_uniq(Lv, (void*)cmp_universal, cmp_nodata)) != l)
    pari_err(e_MISC, "components of Lv must be distinct");
  /* log of prime ideals */
  Llog = cgetg(l, t_VEC);
  /* index in Lchiv corresponding to the places */
  Lpr = cgetg(l, t_VECSMALL);
  Lk1 = cgetg(l, t_VECSMALL);
  Lphi1 = zero_zv(r1);
  Lk2 = zero_zv(r2); nmiss = d;
  for (i = 1, npr = nk1 = 0; i < l; i++)
    if (typ(gel(Lv,i)) == t_INT)
    {
      long v = itos(gel(Lv,i));
      if (v <= 0) pari_err_DOMAIN("gcharidentify", "v", "<=", gen_0, stoi(v));
      if (v > r1+r2)
        pari_err_DOMAIN("gcharidentify", "v", ">", stoi(r1+r2), stoi(v));
      if (v <= r1)
      { /* don't put in k1 if not in conductor (but keep as phi) */
        if (zv_search(moo, v)) { nk1++; Lk1[nk1] = i; }
        Lphi1[v] = i; nmiss--;
      }
      else
      {
        Lk2[v-r1] = i; nmiss -= 2;
      }
    }
    else
    {
      GEN pr = gel(Lv,i);
      long lP = lg(P);
      if (idealtyp(&pr, NULL) != id_PRIME)
        pari_err_TYPE("gcharidentify [ideal]", pr);
      for (j = 1; j < lP; j++)
        if (pr_equal(pr, gel(P,j)))
          pari_err_COPRIME("gcharidentify", pr, gel(gchar_get_mod(gc), 1));
      npr++;
      Lpr[npr] = i;
      gel(Llog,npr) = gchar_log(gc, pr, NULL, prec);
    }
  setlg(Llog, npr+1); setlg(Lpr, npr+1); setlg(Lk1, nk1+1);
  /* build matrix M */
  n = npr + nk1; dim = n + d; ncol = n + 1 + nchi;
  M = cgetg(ncol+1, t_MAT);
  for (j = 1; j <= npr; j++) gel(M,j) = col_ei(dim, j);
  for (;  j <= n; j++) gel(M,j) = col_2ei(dim, j);
  gel(M,j) = zerocol(dim);
  for (i = n+1; i <= n+r1+r2; i++) gcoeff(M,i,j) = eps;
  vlogchi = RgM_mul(gchar_get_Ui(gc), gchar_get_basis(gc));
  for (j=1; j<=nchi; j++)
  {
    GEN logchi = RgV_frac_inplace(row(vlogchi, j), n0); /* duallog(e_j) */
    GEN chi_oo = gcharlog_conductor_oo(gc, logchi), C = cgetg(dim+1, t_COL);
    for (i=1; i<=npr; i++) gel(C,i) = gcharlog_eval_raw(logchi, gel(Llog,i));
    for (i=1; i<=nk1; i++) gel(C,npr+i) = gel(chi_oo, itos(gel(Lv,Lk1[i])));
    for (i=1; i<=d; i++) gel(C,n+i) = gel(logchi, n0 + i);
    gel(M,n+1+j) = C;
  }
  for (i = 1; i <= r1; i++)
    if (!Lphi1[i])
    {
      long a = n + i;
      for (j = 1; j <= ncol; j++) gcoeff(M,a,j) = gmul2n(gcoeff(M,a,j), beps);
    }
  for (i = 1; i <= r2; i++)
    if (!Lk2[i])
    {
      long a = n + r1 + i, b = a + r2;
      for (j = 1; j <= ncol; j++)
      {
        gcoeff(M,a,j) = gmul2n(gcoeff(M,a,j), beps);
        gcoeff(M,b,j) = gmul2n(gcoeff(M,b,j), beps);
      }
    }
  /* rescale and truncate M, then LLL-reduce M */
  M = gtrunc(RgM_Rg_mul(M, mult));
  U = ZM_lll(M, 0.99, LLL_IM);
  M = ZM_mul(M, U);
  U = rowslice(U, n + 2, n + 1 + nchi);
  return mkvecn(9, M, U, Lpr, Lk1, Lphi1, Lk2, mult, Lv,
                   mkvecsmall3(prec, nmiss, beps));
}

/* TODO return the distance between the character found and the conditions? */
static GEN
gchar_identify_i(GEN gc, GEN idinit, GEN Lchiv)
{
  GEN M, U, Lpr, Lk1, Lphi1, Lk2, mult, cyc, y, x, sumphi, Lv, Norm, nf;
  long i, l, r1, r2, beps, npr, nk1, n, nmiss, nnorm, prec;
  M = gel(idinit,1);
  U = gel(idinit,2);
  Lpr = gel(idinit,3);
  Lk1 = gel(idinit,4);
  Lphi1 = gel(idinit,5);
  Lk2 = gel(idinit,6);
  mult = gel(idinit,7);
  Lv = gel(idinit,8);
  prec = gel(idinit,9)[1];
  nmiss = gel(idinit,9)[2];
  beps = gel(idinit,9)[3];
  npr = lg(Lpr)-1;
  nk1 = lg(Lk1)-1; n = npr + nk1;
  cyc = gchar_get_cyc(gc);
  nf = gchar_get_nf(gc);
  nf_get_sign(nf, &r1, &r2);
  nnorm = 0;
  Norm = gen_0;

  l = lg(Lchiv); Lchiv = shallowcopy(Lchiv);
  if (lg(Lv) != l) pari_err_DIM("gcharidentify [#Lv != #Lchiv]");
  for (i = 1; i < l; i++)
  {
    GEN t = gel(Lchiv,i), u;
    if (typ(gel(Lv,i)) != t_INT)
    {
      if (typ(t) == t_VEC) /* value at last component */
          gel(Lchiv,i) = t = gel(t, lg(t)-1);
      if (typ(t) == t_COMPLEX)
      {
        nnorm++; /* 2 Pi Im(theta) / log N(pr) */
        Norm = gadd(Norm, gdiv(gmul(Pi2n(1,prec), gel(t,2)),
                               glog(idealnorm(nf,gel(Lv,i)),prec)));
        gel(Lchiv,i) = t = gel(t,1);
      }
      if (!is_real_t(typ(t)))
        pari_err_TYPE("gcharidentify [character value: should be real or complex]", t);
    }
    else
    {
      if (typ(t) != t_VEC)
        pari_err_TYPE("gcharidentify [character at infinite place: should be a t_VEC]", t);
      if (lg(t) != 3)
        pari_err_DIM("gcharidentify [character at infinite place: should have two components]");
      if (typ(gel(t,1)) != t_INT)
        pari_err_TYPE("gcharidentify [k parameter at infinite place: should be a t_INT]", gel(t,1));
      u = gel(t,2);
      if (typ(u) == t_COMPLEX)
      {
        nnorm++;
        Norm = gsub(Norm, gel(u,2)); u = gel(u,1);
        gel(Lchiv, i) = mkvec2(gel(t,1), u);
      }
      if (!is_real_t(typ(u)))
        pari_err_TYPE("gcharidentify [phi parameter at infinite place: should be real or complex]", u);
    }
  }

  /* construct vector */
  y = zerocol(n + r1 + 2*r2); sumphi = gen_0;
  for (i=1; i<=npr; i++) gel(y,i) = gel(Lchiv, Lpr[i]);
  for (i=1; i<=nk1; i++) gel(y,npr+i) = gmael(Lchiv,Lk1[i],1);
  for (i=1; i<=r1; i++)
    if (Lphi1[i])
    {
      gel(y, n+i) = x =  gmael(Lchiv,Lphi1[i],2);
      sumphi = gadd(sumphi, x);
    }
  for (i=1; i<=r2; i++)
    if (Lk2[i])
    {
      long a = n + r1 + i;
      gel(y, a + r2) = gmael(Lchiv,Lk2[i],1);
      gel(y, a) = x =  gmael(Lchiv,Lk2[i],2);
      sumphi = gadd(sumphi, gshift(x,1));
    }
  if (nmiss)
  {
    sumphi = gmul2n(gdivgs(sumphi, -nmiss), beps);
    for (i = 1; i <= r1; i++) if (!Lphi1[i]) gel(y, n + i) = sumphi;
    for (i = 1; i <= r2; i++) if (!Lk2[i])   gel(y, n + r1+i) = sumphi;
  }
  y = gtrunc(RgC_Rg_mul(y, mult));

  /* find approximation */
  x = ZM_ZC_mul(U, RgM_Babai(M, y));
  for (i = 1; i < lg(cyc)-1; i++) /* ignore norm */
    if (signe(gel(cyc,i))) gel(x,i) = modii(gel(x,i), gel(cyc,i));
  if (nnorm > 0) x = vec_append(x, gdivgu(Norm, lg(Lv)-1));
  return x;
}

/* TODO export the init interface */
GEN
gchar_identify(GEN gc, GEN Lv, GEN Lchiv, long prec)
{
  pari_sp av = avma;
  GEN idinit = gchar_identify_init(gc, Lv, prec);
  return gerepilecopy(av, gchar_identify_i(gc,idinit,Lchiv));
}

/*******************************************************************/
/*                                                                 */
/*                          L FUNCTION                             */
/*                                                                 */
/*******************************************************************/

/* TODO: could merge with vecan_chigen */

/* return x + yz; y != 0; z = 0,1 "often"; x = 0 "often" */
static GEN
gaddmul(GEN x, GEN y, GEN z)
{
  pari_sp av;
  if (typ(z) == t_INT)
  {
    if (!signe(z)) return x;
    if (equali1(z)) return gadd(x,y);
  }
  if (isintzero(x)) return gmul(y,z);
  av = avma;
  return gerepileupto(av, gadd(x, gmul(y,z)));
}

GEN
vecan_gchar(GEN an, long n, long prec)
{
  forprime_t iter;
  GEN gc = gel(an,1), chi = gel(an,2), P = gel(an,3), PZ = gel(an,4);
  GEN BOUND = stoi(n), v = vec_ei(n, 1);
  GEN gp = cgetipos(3), nf, chilog, s;
  ulong p;

  /* prec increase: 1/n*log(N(pmax)) < log(pmax) */
  if (DEBUGLEVEL > 1)
    err_printf("vecan_gchar: need extra prec %ld\n", nbits2extraprec(expu(n)));
  gc = gcharnewprec(gc, prec + nbits2extraprec(expu(n)));
  chilog = gchari_duallog(gc, check_gchari(gc, chi, &s));

  nf = gchar_get_nf(gc);
  /* FIXME: when many log of many primes are computed:
     - bnfisprincipal may not be improved
     - however we can precompute the logs of generators
       for principal part
     - when galois, should compute one ideal by orbit.
     - when real, clear imaginary part
   */

  u_forprime_init(&iter, 2, n);
  while ((p = u_forprime_next(&iter)))
  {
    GEN L;
    long j;
    int check = !umodiu(PZ,p);
    gp[2] = p;
    L = idealprimedec_limit_norm(nf, gp, BOUND);
    for (j = 1; j < lg(L); j++)
    {
      GEN pr = gel(L, j), ch;
      pari_sp av;
      ulong k, q;
      if (check && gen_search(P, pr, (void*)cmp_prime_ideal, cmp_nodata) > 0)
        continue;
      /* TODO: extract code and use precom sprk? */
      av = avma;
      ch = gchari_eval(gc, chi, pr, 1, chilog, gen_0, prec);
      ch = gerepileupto(av, ch);
      q = upr_norm(pr);
      gel(v, q) = gadd(gel(v, q), ch);
      for (k = 2*q; k <= (ulong)n; k += q)
        gel(v, k) = gaddmul(gel(v, k), ch, gel(v, k/q));
    }
  }
  /* weight, could consider shifting s at eval instead */
  if (!gequal0(s))
  {
    GEN ns = dirpowers(n, gneg(s), prec);
    long j;
    for (j = 1; j <= n; j++)
      if (gel(v,j) != gen_0) gel(v, j) = gmul(gel(v,j),gel(ns,j));
  }
  return v;
}

GEN
eulerf_gchar(GEN an, GEN p, long prec)
{
  GEN gc = gel(an,1), chi = gel(an,2), P = gel(an,3), PZ = gel(an,4);
  GEN f, L, nf, chilog, s;
  int check;
  long j, l;

  /* prec increase: 1/n*log(N(pmax)) < log(pmax) */
  if (DEBUGLEVEL > 1)
    err_printf("vecan_gchar: need extra prec %ld\n", nbits2extraprec(expi(p)));
  gc = gcharnewprec(gc, prec + nbits2extraprec(expi(p)));
  chilog = gchari_duallog(gc, check_gchari(gc, chi, &s));

  nf = gchar_get_nf(gc);
  f = pol_1(0);
  check = dvdii(PZ, p);
  L = idealprimedec(nf, p); l = lg(L);
  for (j = 1; j < l; j++)
  {
    GEN pr = gel(L, j), ch;
    if (check && gen_search(P, pr, (void*)cmp_prime_ideal, cmp_nodata) > 0)
      continue;
    ch =  gchari_eval(gc, chi, pr, 1, chilog, s, prec);
    f = gmul(f, gsub(gen_1, monomial(ch, pr_get_f(pr), 0)));
  }
  return mkrfrac(gen_1,f);
}

static GEN
cleanup_vga(GEN vga, long prec)
{
  GEN ind;
  long bitprec, i, l;
  if (!prec) return vga; /* already exact */
  bitprec = bit_accuracy(prec);
  vga = shallowcopy(vga); l = lg(vga);
  for (i = 1; i < l; i++)
  {
    GEN z = gel(vga,i);
    if (typ(z) != t_COMPLEX) continue;
    if (gexpo(gel(z,2)) < -bitprec+20) gel(vga,i) = gel(z,1);
  }
  ind = indexsort(imag_i(vga));
  for (i = 2; i < l; i++)
  {
    GEN z = gel(vga,ind[i]), t;
    if (typ(z) != t_COMPLEX) continue;
    t = imag_i(gel(vga, ind[i-1]));
    if (gexpo(gsub(gel(z,2), t)) < -bitprec+20)
      gel(vga, ind[i]) = mkcomplex(gel(z,1), t);
   }
  for (i = 1; i < l; i++)
  {
    GEN z = gel(vga,i);
    if (typ(z) != t_COMPLEX) continue;
    gel(vga, i) = mkcomplex(gel(z,1), bestappr(gel(z,2), int2n(bitprec/2)));
  }
  return vga;
}

/* TODO: move to lfunutils, use lfunzeta and lfunzetak */
GEN
gchari_lfun(GEN gc, GEN chi, GEN s0)
{
  GEN nf, chilog, s, cond_f, cond_oo, vga_r, vga_c, chiw;
  GEN v_phi, v_arg, sig, k, NN, faN, P;
  long ns, nc, nm, r1, r2;

  nf = gchar_get_nf(gc);
  ns = gchar_get_ns(gc);
  nc = gchar_get_nc(gc);
  nm = gchar_get_nm(gc);
  nf_get_sign(nf, &r1, &r2);
  chi = check_gchari(gc, chi, &s);
  chilog = gchari_duallog(gc, chi);
  s = gadd(s0,s); chiw = shallowconcat(chi, s);
  if (!gequal0(imag_i(s)))
    pari_err_IMPL("lfun for gchar with imaginary norm component");
  cond_f =  gcharlog_conductor_f(gc, chilog, &faN);
  P = gel(faN, 1); /* prime ideals dividing cond(chi) */
  cond_oo =  gcharlog_conductor_oo(gc, chilog);

  NN = mulii(idealnorm(nf, cond_f), absi_shallow(nf_get_disc(nf)));
  if (equali1(NN)) return lfunshift(lfuncreate(gen_1), gneg(s), 0,
      prec2nbits(gchar_get_evalprec(gc)));
  if (ZV_equal0(chi)) return lfunshift(lfuncreate(nf), gneg(s), 0,
      prec2nbits(gchar_get_evalprec(gc)));

  /* vga_r = vector(r1,k,I*c[ns+nc+k]-s + cond_oo[k]);
   * vga_c = vector(r2,k,I*c[ns+nc+r1+k]+abs(c[ns+nc+r1+r2+k])/2-s) */
  v_phi = gmul(vecslice(chilog, ns+nc+1, ns+nc+r1+r2), gen_I());
  v_arg = gdivgs(gabs(vecslice(chilog,ns+nc+r1+r2+1,nm),BITS_IN_LONG), 2);
  vga_r = gadd(vecslice(v_phi, 1, r1), cond_oo);
  vga_c = gadd(vecslice(v_phi, r1+1, r1+r2), v_arg);
  sig = shallowconcat1(mkvec3(vga_r,vga_c,gadd(vga_c,const_vec(r2,gen_1))));
  /* TODO: remove cleanup when gammamellinv takes ldata*/
  sig = cleanup_vga(sig, gchar_get_prec(gc));
  k = gen_1;
  if (!gequal0(s))
  {
    long j;
    for (j = 1; j < lg(sig); j++) gel(sig, j) = gadd(gel(sig, j), s);
    k = gsub(k, gmulgs(s,2));
  }

#define tag(x,t)  mkvec2(mkvecsmall(t),x)
  return mkvecn(6, tag(mkvec4(gc, chiw, P, prV_lcm_capZ(P)), t_LFUN_HECKE),
                gen_1, sig, k, NN, gen_0);
}

GEN
lfungchar(GEN gc, GEN chi)
{
  pari_sp av = avma;
  GEN s;
  check_gchar_group(gc);
  chi = gchar_internal(gc, chi, &s);
  return gerepilecopy(av, gchari_lfun(gc, chi, s));
}
