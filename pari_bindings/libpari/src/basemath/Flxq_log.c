/* Copyright (C) 2013 The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_fflog

/* Let [ be the following order on Fp: 0 [ p-1 [ 1 [ p-2 [ 2 .. [ p\2
and [[ the lexicographic extension of [ to Fp[T]. Compute the
isomorphism (Fp[X], [[) -> (N,<) on P */

static long
Flx_cindex(GEN P, ulong p)
{
  long d = degpol(P), i;
  ulong s = 0, p2 = (p-1)>>1;
  for (i = 0; i <= d; ++i)
  {
    ulong x = P[d-i+2];
    if (x<=p2) x = 2*x; else x = 1+2*(p-1-x);
    s = p*s+x;
  }
  return s;
}

/* Compute the polynomial immediately after t for the [[ order */

static void
Flx_cnext(GEN t, ulong p)
{
  long i;
  long p2 = p>>1;
  for(i=2;;i++)
    if (t[i]==p2)
      t[i]=0;
    else
    {
      t[i] = t[i]<p2 ? p-1-t[i]: p-t[i];
      break;
    }
}

static int
has_deg1_auto(GEN T, ulong p, ulong pi)
{
  long i, n = degpol(T);
  GEN a = polx_Flx(get_Flx_var(T));
  for (i=1; i<n; i++)
  {
    a = Flxq_powu_pre(a, p, T, p, pi);
    if (degpol(a)==1) return 1;
  }
  return 0;
}

static void
smallirred_Flx_next(GEN a, long p, ulong pi)
{
  do
  {
    long i;
    for(i=2;;i++)
      if (++a[i]==p) a[i]=0;
      else break;
  } while (!Flx_is_irred(a, p) || has_deg1_auto(a,p,pi) );
}

/* Avoid automorphisms of degree 1 */
static GEN
smallirred_Flx(long p, ulong n, long sv, ulong pi)
{
  GEN a = zero_zv(n+2);
  a[1] = sv; a[3] = 1; a[n+2] = 1;
  smallirred_Flx_next(a, p, pi);
  return a;
}

struct Flxq_log_rel
{
  long nbrel;
  GEN rel;
  long nb;
  long r, off, nbmax, nbexp;
  ulong nbtest;
};

static GEN
cindex_Flx(long c, long d, ulong p, long v)
{
  GEN P = cgetg(d+3, t_VECSMALL);
  long i;
  P[1] = v;
  for (i = 0; i <= d; ++i)
  {
    ulong x = c%p;
    P[i+2] = (x&1) ? p-1-(x>>1) : x>>1;
    c/=p;
  }
  return Flx_renormalize(P, d+3);
}

static GEN
factorel(GEN h, ulong p)
{
  GEN F = Flx_factor(h, p);
  GEN F1 = gel(F, 1), F2 = gel(F, 2);
  long i, l1 = lg(F1)-1;
  GEN p2 = cgetg(l1+1, t_VECSMALL);
  GEN e2 = cgetg(l1+1, t_VECSMALL);
  for (i = 1; i <= l1; ++i)
  {
    p2[i] = Flx_cindex(gel(F1, i), p);
    e2[i] = F2[i];
  }
  return mkmat2(p2, e2);
}

static long
Flx_addifsmooth3(pari_sp *av, struct Flxq_log_rel *r, GEN h, long u, long v, long w, ulong p)
{
  long off = r->off;
  r->nbtest++;
  if (Flx_is_smooth(h, r->r, p))
  {
    GEN z = factorel(h, p);
    if (v<0)
      z = mkmat2(vecsmall_append(gel(z,1),off+u),vecsmall_append(gel(z,2),-1));
    else
      z = famatsmall_reduce(mkmat2(
            vecsmall_concat(gel(z,1),mkvecsmall3(off+u,off+v,off+w)),
            vecsmall_concat(gel(z,2),mkvecsmall3(-1,-1,-1))));
    gel(r->rel,++r->nbrel) = gerepilecopy(*av,z);
    if (DEBUGLEVEL && (r->nbrel&511UL)==0)
      err_printf("%ld%% ",r->nbrel*100/r->nbexp);
    *av = avma;
  } else set_avma(*av);
  return r->nbrel==r->nb || r->nbrel==r->nbmax;
}

static void
Flx_renormalize_inplace(GEN x, long lx)
{
  long i;
  for (i = lx-1; i>1; i--)
    if (x[i]) break;
  setlg(x, i+1);
}

/*
   Let T*X^e=C^3-R
   a+b+c = 0
   (C+a)*(C+b)*(C+c) = C^3+ (a*b+a*c+b*c)*C+a*b*c
   = R + (a*b+a*c+b*c)*C+a*b*c
   = R + (a*b-c^2)*C+a*b*c
 */
static void
Flxq_log_cubic(struct Flxq_log_rel *r, GEN C, GEN R, ulong p, ulong pi)
{
  long l = lg(C);
  GEN a = zero_zv(l); /*We allocate one extra word to catch overflow*/
  GEN b = zero_zv(l);
  pari_sp av = avma;
  long i,j,k;
  for(i=0; ; i++, Flx_cnext(a, p))
  {
    Flx_renormalize_inplace(a, l+1);
    r->nb++;
    if (Flx_addifsmooth3(&av, r, Flx_add(a, C, p), i, -1, -1, p)) return;
    for(j=2; j<=l; j++) b[j] = 0;
    for(j=0; j<=i; j++, Flx_cnext(b, p))
    {
      GEN h,c;
      GEN pab,pabc,pabc2;
      Flx_renormalize_inplace(b, l+1);
      c = Flx_neg(Flx_add(a,b,p),p);
      k = Flx_cindex(c, p);
      if (k > j) continue;
      pab  = Flx_mul_pre(a, b, p, pi);
      pabc = Flx_mul_pre(pab,c,p,pi);
      pabc2= Flx_sub(pab,Flx_sqr_pre(c,p,pi),p);
      h = Flx_add(R,Flx_add(Flx_mul_pre(C,pabc2,p,pi),pabc,p), p);
      h = Flx_normalize(h, p);
      if (Flx_addifsmooth3(&av, r, h, i, j, k, p)) return;
    }
  }
}

static GEN
Flxq_log_find_rel(GEN b, long r, GEN T, ulong p, ulong pi, GEN *g, long *e)
{
  pari_sp av = avma;
  while (1)
  {
    GEN M, z;
    *g = Flxq_mul_pre(*g, b, T, p, pi); (*e)++;
    M = Flx_halfgcd_all_pre(*g,T,p,pi,&z,NULL);
    if (Flx_is_smooth_pre(gcoeff(M,1,1), r, p, pi)
     && Flx_is_smooth_pre(z, r, p, pi))
    {
      GEN F = factorel(z, p);
      GEN G = factorel(gcoeff(M,1,1), p);
      GEN rel = mkmat2(vecsmall_concat(gel(F, 1),gel(G, 1)),
          vecsmall_concat(gel(F, 2),zv_neg(gel(G, 2))));
      return gc_all(av,2,&rel,g);
    }
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"Flxq_log_find_rel");
      *g = gerepilecopy(av, *g);
    }
  }
}

/* Generalised Odlyzko formulae ( EUROCRYPT '84, LNCS 209, pp. 224-314, 1985. ) */
/* Return the number of monic, k smooth, degree n polynomials for k=1..r */
static GEN
smoothness_vec(ulong p, long r, long n)
{
  long i,j,k;
  GEN R = cgetg(r+1, t_VEC), pp = utoipos(p);
  GEN V = cgetg(n+1, t_VEC);
  for (j = 1; j <= n; ++j)
    gel(V, j) =  binomialuu(p+j-1,j);
  gel(R, 1) = gel(V, n);
  for (k = 2; k <= r; ++k)
  {
    GEN W = cgetg(n+1, t_VEC);
    GEN Ik = ffnbirred(pp, k);
    for (j = 1; j <= n; ++j)
    {
      long l = j/k;
      GEN s = gen_0;
      pari_sp av2 = avma;
      if (l*k == j)
      {
        s = binomial(addiu(Ik,l-1), l);
        l--;
      }
      for (i = 0; i <= l; ++i)
        s = addii(s, mulii(gel(V, j-k*i), binomial(addis(Ik,i-1), i)));
      gel(W, j) = gerepileuptoint(av2, s);
    }
    V = W;
    gel(R, k) = gel(V, n);
  }
  return R;
}

/* Solve N^2*pr/6 + N*prC = N+fb
   N^2*pr/6 + N*(prC-1) -fb = 0
 */

static GEN
smooth_cost(GEN fb, GEN pr, GEN prC)
{
  GEN a = gdivgu(pr,6);
  GEN b = gsubgs(prC,1);
  GEN c = gneg(fb);
  GEN vD = gsqrt(gsub(gsqr(b),gmul2n(gmul(a,c),2)),BIGDEFAULTPREC);
  return ceil_safe(gdiv(gsub(vD,b),gmul2n(a,1)));
}

/* Return best choice of r.
   We loop over d until there is sufficiently many triples (a,b,c) (a+b+c=0)
   of degree <=d with respect to the probability of smoothness of (a*b-c^2)*C
 */

static GEN
smooth_best(long p, long n, long *pt_r, long *pt_nb)
{
  pari_sp av = avma, av2;
  GEN bestc = NULL, pp = utoipos(p);
  long bestr = 0, bestFB = 0;
  long r,d, dC = (n+2)/3;
  for (r = 1; r < dC; ++r)
  {
    GEN fb = ffsumnbirred(pp, r);
    GEN smoothC = smoothness_vec(p,r,dC);
    GEN prC = gdiv(gel(smoothC,r), powuu(p,dC));
    ulong rels = 0;
    av2 = avma;
    for(d=0; d<dC && rels < ULONG_MAX; d++)
    {
      GEN c;
      long dt = dC+2*d;
      GEN smooth = smoothness_vec(p,r,dt);
      GEN pr = gdiv(gel(smooth,r), powuu(p,dt));
      GEN FB = addii(fb,powuu(p,d));
      GEN N = smooth_cost(subiu(FB,rels),pr,prC);
      GEN Nmax = powuu(p,d+1);
      if (gcmp(N,Nmax) >= 0)
      {
        rels = itou_or_0(addui(rels, gceil(gmul(gdivgu(sqri(Nmax),6),pr))));
        if (!rels) rels = ULONG_MAX;
        set_avma(av2);
        continue;
      }
      c = gdivgu(addii(powuu(p,2*d),sqri(N)),6);
      FB = addii(FB,N);
      if ((!bestc || gcmp(gmul2n(c,r), gmul2n(bestc,bestr)) < 0))
      {
        if (DEBUGLEVEL)
          err_printf("r=%ld d=%ld fb=%Ps early rels=%lu P=%.5Pe -> C=%.5Pe \n",
                      r, dt, FB, rels, pr, c);
        bestc = c;
        bestr = r;
        bestFB = itos_or_0(FB);
      }
      break;
    }
  }
  *pt_r=bestr;
  *pt_nb=bestFB;
  return bestc ? gerepileupto(av, gceil(bestc)): NULL;
}

static GEN
check_kernel(long r, GEN M, long nbi, long nbrow, GEN T, ulong p, ulong pi, GEN m)
{
  pari_sp av = avma;
  long N = 3*upowuu(p, r);
  GEN K = FpMs_leftkernel_elt(M, nbrow, m);
  long i, f=0, tbs;
  long lm = lgefint(m), u=1;
  GEN tab, g;
  GEN q = powuu(p,degpol(T));
  GEN idx = diviiexact(subiu(q,1),m);
  pari_timer ti;
  if (DEBUGLEVEL) timer_start(&ti);
  while (signe(gel(K,u))==0)
    u++;
  K = FpC_Fp_mul(K, Fp_inv(gel(K, u), m), m);
  g = Flxq_pow_pre(cindex_Flx(u, r, p, T[1]), idx, T, p, pi);
  tbs = maxss(1, expu(nbi/expi(m)));
  tab = Flxq_pow_init_pre(g, q, tbs, T, p, pi);
  setlg(K, N);
  for (i=1; i<N; i++)
  {
    GEN k = gel(K,i);
    pari_sp av = avma;
    long t = signe(k) && Flx_equal(Flxq_pow_table_pre(tab, k, T, p, pi),
                                   Flxq_pow_pre(cindex_Flx(i,r,p,T[1]), idx, T, p, pi));
    set_avma(av);
    if (!t)
      gel(K,i) = cgetineg(lm);
    else
      f++;
  }
  if (DEBUGLEVEL) timer_printf(&ti,"found %ld/%ld logs", f, nbi);
  if (f < maxss(3,maxss(p/2,nbi/p))) return NULL; /* Not enough logs found */
  return gerepilecopy(av, K);
}

static GEN
Flxq_log_rec(GEN W, GEN a, long r, GEN T, ulong p, ulong pi, GEN m)
{
  long AV = 0, u = 1;
  GEN g = a, b;
  pari_timer ti;
  while (!equali1(gel(W,u)))
   u++;
  b = cindex_Flx(u, r, p, T[1]);
  while(1)
  {
    long i, l;
    GEN V, F, E, Ao;
    timer_start(&ti);
    V = Flxq_log_find_rel(b, r, T, p, pi, &g, &AV);
    if (DEBUGLEVEL>1) timer_printf(&ti,"%ld-smooth element",r);
    F = gel(V,1); E = gel(V,2);
    l = lg(F);
    Ao = gen_0;
    for(i=1; i<l; i++)
    {
      GEN R = gel(W,F[i]);
      if (signe(R)<=0)
        break;
      Ao = Fp_add(Ao, mulis(R, E[i]), m);
    }
    if (i==l) return subis(Ao,AV);
  }
}

static int
Flxq_log_use_index_cubic(GEN m, GEN T0, ulong p)
{
  pari_sp av = avma;
  long n = get_Flx_degree(T0), r, nb;
  GEN cost = smooth_best(p, n, &r, &nb);
  GEN cost_rho = sqrti(shifti(m,2));
  int use = (cost && gcmp(cost,cost_rho)<0);
  set_avma(av);
  return use;
}

static GEN
Flxq_log_index_cubic(GEN a0, GEN b0, GEN m, GEN T0, ulong p)
{
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  long n = get_Flx_degree(T0), r, nb;
  pari_sp av = avma;
  struct Flxq_log_rel rel;
  long nbi;
  GEN W, M, S, T, a, b, Ao, Bo, e, C, R;
  pari_timer ti;
  GEN cost = smooth_best(p, n, &r, &nb);
  GEN cost_rho = sqrti(shifti(m,2));
  if (!cost || gcmp(cost,cost_rho)>=0) return gc_NULL(av);
  nbi = itos(ffsumnbirred(stoi(p), r));
  if (DEBUGLEVEL)
  {
    err_printf("Size FB=%ld, looking for %ld relations, %Ps tests needed\n", nbi, nb,cost);
    timer_start(&ti);
  }
  T = smallirred_Flx(p,n,get_Flx_var(T0), pi);
  for(;;)
  {
    S = Flx_ffisom(T0,T,p);
    a = Flx_Flxq_eval_pre(a0, S, T, p, pi);
    b = Flx_Flxq_eval_pre(b0, S, T, p, pi);
    C = Flx_shift(pol1_Flx(get_Flx_var(T)), (n+2)/3);
    R = Flxq_powu_pre(C,3,T,p,pi);
    if (DEBUGLEVEL)
      timer_printf(&ti," model change: %Ps",Flx_to_ZX(T));
    rel.nbmax=2*nb;
    M = cgetg(rel.nbmax+1, t_VEC);
    rel.rel = M;
    rel.nbrel = 0; rel.r = r; rel.off = 3*upowuu(p,r);
    rel.nb = nbi; rel.nbexp = nb; rel.nbtest=0;
    Flxq_log_cubic(&rel, C, R, p, pi);
    setlg(M,1+rel.nbrel);
    if (DEBUGLEVEL)
    {
      err_printf("\n");
      timer_printf(&ti," %ld relations, %ld generators (%ld tests)",rel.nbrel,rel.nb,rel.nbtest);
    }
    W = check_kernel(r, M, nbi, rel.off + rel.nb - nbi, T, p, pi, m);
    if (W) break;
    if (DEBUGLEVEL) timer_start(&ti);
    smallirred_Flx_next(T,p, pi);
  }
  if (DEBUGLEVEL) timer_start(&ti);
  Ao = Flxq_log_rec(W, a, r, T, p, pi, m);
  if (DEBUGLEVEL) timer_printf(&ti,"smooth element");
  Bo = Flxq_log_rec(W, b, r, T, p, pi, m);
  if (DEBUGLEVEL) timer_printf(&ti,"smooth generator");
  e = Fp_div(Ao, Bo, m);
  if (!Flx_equal(Flxq_pow_pre(b0, e, T0, p, pi), a0)) pari_err_BUG("Flxq_log");
  return gerepileupto(av, e);
}

INLINE GEN Flx_frob(GEN u, ulong p) { return Flx_inflate(u, p); }

static GEN
rel_Coppersmith(long r, GEN u, GEN v, long h, GEN R, long d, ulong p, ulong pi)
{
  GEN a, b, F, G, M;
  if (degpol(Flx_gcd_pre(u,v,p,pi))) return NULL;
  a = Flx_add(Flx_shift(u, h), v, p);
  if (lgpol(a)==0 || !Flx_is_smooth_pre(a, r, p, pi)) return NULL;
  b = Flx_add(Flx_mul_pre(R, Flx_frob(u, p), p, pi),
              Flx_shift(Flx_frob(v, p),d), p);
  if (!Flx_is_smooth_pre(b, r, p, pi)) return NULL;
  F = factorel(a, p); G = factorel(b, p);
  M = mkmat2(vecsmall_concat(gel(F, 1), vecsmall_append(gel(G, 1), 2*p)),
             vecsmall_concat(zv_z_mul(gel(F, 2),p), vecsmall_append(zv_neg(gel(G, 2)),d)));
  return famatsmall_reduce(M);
}

GEN
Flxq_log_Coppersmith_worker(GEN u, long i, GEN V, GEN R)
{
  long r = V[1], h = V[2], d = V[3], p = V[4], pi = V[5], dT = V[6];
  pari_sp ltop = avma;
  GEN v = zero_zv(dT+2);
  GEN L = cgetg(2*i+1, t_VEC);
  pari_sp av = avma;
  long j;
  long nbtest=0, rel = 1;
  ulong lu = Flx_lead(u), lv;
  for (j=1; j<=i; j++)
  {
    GEN z;
    Flx_cnext(v, p);
    Flx_renormalize_inplace(v, dT+2);
    lv = Flx_lead(v);
    set_avma(av);
    if (lu != 1 && lv != 1) continue;
    if (degpol(Flx_gcd_pre(u, v, p, pi))!=0) continue;
    if (lu==1)
    {
      z = rel_Coppersmith(r, u, v, h, R, d, p, pi);
      nbtest++;
      if (z) { gel(L, rel++) = z; av = avma; }
    }
    if (i==j) continue;
    if (lv==1)
    {
      z = rel_Coppersmith(r, v, u, h, R, d, p, pi);
      nbtest++;
      if (z) { gel(L, rel++) = z; av = avma; }
    }
  }
  setlg(L,rel);
  return gerepilecopy(ltop, mkvec2(stoi(nbtest), L));
}

static GEN
Flxq_log_Coppersmith(long nbrel, long r, GEN T, ulong p, ulong pi)
{
  pari_sp av;
  long dT = degpol(T);
  long h = dT/p, d = dT-(h*p);
  GEN R = Flx_sub(Flx_shift(pol1_Flx(T[1]), dT), T, p);
  GEN u = zero_zv(dT+2);
  GEN done;
  long nbtest = 0, rel = 0;
  GEN M = cgetg(nbrel+1, t_VEC);
  long i = 1;
  GEN worker = snm_closure(is_entry("_Flxq_log_Coppersmith_worker"),
               mkvec2(mkvecsmalln(6, r,h,d,p,pi,dT), R));
  struct pari_mt pt;
  long running, pending = 0, stop=0;
  if (DEBUGLEVEL) err_printf("Coppersmith (R = %ld): ",degpol(R));
  mt_queue_start(&pt, worker);
  av = avma;
  while ((running = !stop) || pending)
  {
    GEN L;
    long l, j;
    Flx_cnext(u, p);
    Flx_renormalize_inplace(u, dT+2);
    mt_queue_submit(&pt, 0, running ? mkvec2(u, stoi(i)): NULL);
    done = mt_queue_get(&pt, NULL, &pending);
    if (!done) continue;
    L = gel(done, 2); nbtest += itos(gel(done,1));
    l = lg(L);
    if (l > 1)
    {
      for (j=1; j<l; j++)
      {
        if (rel>nbrel) break;
        gel(M,++rel) = gel(L,j);
        if (DEBUGLEVEL && (rel&511UL)==0)
          err_printf("%ld%%[%ld] ",rel*100/nbrel,i);
      }
      av = avma;
    }
    else set_avma(av);
    if (rel>nbrel) stop = 1;
    i++;
  }
  mt_queue_end(&pt);
  if (DEBUGLEVEL) err_printf(": %ld tests\n", nbtest);
  return M;
}

static GEN Flxq_log_Coppersmith_d(GEN W, GEN g, long r, GEN T, ulong p, ulong pi, GEN mo);

static GEN
Flxq_log_from_rel(GEN W, GEN rel, long r, GEN T, ulong p, ulong pi, GEN m)
{
  pari_sp av = avma;
  GEN F = gel(rel,1), E = gel(rel,2), o = gen_0;
  long i, l = lg(F);
  for(i=1; i<l; i++)
  {
    GEN R = gel(W, F[i]);
    if (signe(R)==0) /* Already failed */
      return NULL;
    else if (signe(R)<0) /* Not yet tested */
    {
      setsigne(gel(W,F[i]),0);
      R = Flxq_log_Coppersmith_d(W, cindex_Flx(F[i],r,p,T[1]), r, T, p, pi, m);
      if (!R) return NULL;
    }
    o = Fp_add(o, mulis(R, E[i]), m);
  }
  return gerepileuptoint(av, o);
}

static GEN
Flxq_log_Coppersmith_d(GEN W, GEN g, long r, GEN T, ulong p, ulong pi, GEN mo)
{
  pari_sp av = avma, av2;
  long dg = degpol(g), k = r-1, m = maxss((dg-k)/2,0);
  long i, j, l = dg-m, N;
  GEN v = cgetg(k+m+1,t_MAT);
  long dT = degpol(T);
  long h = dT/p, d = dT-h*p;
  GEN R = Flx_rem_pre(Flx_shift(pol1_Flx(T[1]), dT), T, p, pi);
  GEN z = Flx_rem_pre(Flx_shift(pol1_Flx(T[1]), h), g, p, pi);
  for(i=1; i<=k+m; i++)
  {
    gel(v,i) = Flx_to_Flv(Flx_shift(z,-l),m);
    z = Flx_rem_pre(Flx_shift(z,1),g,p,pi);
  }
  v = Flm_ker(v,p);
  for(i=1; i<=k; i++)
    gel(v,i) = Flv_to_Flx(gel(v,i),T[1]);
  N = upowuu(p,k);
  av2 = avma;
  for (i=1; i<N; i++)
  {
    GEN p0,q,qh,a,b;
    ulong el = i;
    set_avma(av2);
    q = pol0_Flx(T[1]);
    for (j=1; j<=k; j++)
    {
      ulong r = el % p;
      el /= p;
      if (r) q = Flx_add(q, Flx_Fl_mul(gel(v,j), r, p), p);
    }
    qh = Flx_shift(q, h);
    p0 = Flx_rem_pre(qh, g, p, pi);
    b = Flx_sub(Flx_mul_pre(R, Flx_frob(q, p), p, pi),
                Flx_shift(Flx_frob(p0, p), d), p);
    if (lgpol(b)==0 || !Flx_is_smooth_pre(b, r, p, pi)) continue;
    a = Flx_div_pre(Flx_sub(qh, p0, p), g, p, pi);
    if (degpol(Flx_gcd_pre(a, q, p, pi)) && degpol(Flx_gcd_pre(a, p0, p, pi)))
      continue;
    if (!(lgpol(a)==0 || !Flx_is_smooth_pre(a, r, p, pi)))
    {
      GEN F = factorel(b, p);
      GEN G = factorel(a, p);
      GEN FG = vecsmall_concat(vecsmall_append(gel(F, 1), 2*p), gel(G, 1));
      GEN E  = vecsmall_concat(vecsmall_append(gel(F, 2), -d),
          zv_z_mul(gel(G, 2),-p));
      GEN R  = famatsmall_reduce(mkmat2(FG, E));
      GEN l  = Flxq_log_from_rel(W, R, r, T, p, pi, mo);
      if (!l) continue;
      l = Fp_divu(l,p,mo);
      if (dg <= r)
      {
        long idx = Flx_cindex(g, p);
        affii(l, gel(W, idx));
        if (DEBUGLEVEL>1) err_printf("Found %lu\n", idx);
      }
      return gerepileuptoint(av, l);
    }
  }
  set_avma(av);
  return NULL;
}

static GEN
Flxq_log_Coppersmith_rec(GEN W, long r2, GEN a, long r, GEN T, ulong p, ulong pi, GEN m)
{
  GEN b = polx_Flx(T[1]);
  long AV = 0;
  GEN g = a, bad = pol0_Flx(T[1]);
  pari_timer ti;
  while(1)
  {
    long i, l;
    GEN V, F, E, Ao;
    timer_start(&ti);
    V = Flxq_log_find_rel(b, r2, T, p, pi, &g, &AV);
    if (DEBUGLEVEL>1) timer_printf(&ti,"%ld-smooth element",r2);
    F = gel(V,1); E = gel(V,2);
    l = lg(F);
    Ao = gen_0;
    for(i=1; i<l; i++)
    {
      GEN Fi = cindex_Flx(F[i], r2, p, T[1]);
      GEN R;
      if (degpol(Fi) <= r)
      {
        if (signe(gel(W,F[i]))==0)
          break;
        else if (signe(gel(W,F[i]))<0)
        {
          setsigne(gel(W,F[i]),0);
          R = Flxq_log_Coppersmith_d(W,Fi,r,T,p,pi,m);
        } else
          R = gel(W,F[i]);
      }
      else
      {
        if (Flx_equal(Fi,bad)) break;
        R = Flxq_log_Coppersmith_d(W,Fi,r,T,p,pi,m);
        if (!R) bad = Fi;
      }
      if (!R) break;
      Ao = Fp_add(Ao, mulis(R, E[i]), m);
    }
    if (i==l) return subis(Ao,AV);
  }
}

static GEN
Flxq_log_index_Coppersmith(GEN a0, GEN b0, GEN m, GEN T0, ulong p)
{
  ulong pi = SMALL_ULONG(p)? 0: get_Fl_red(p);
  pari_sp av = avma;
  GEN  M, S, a, b, Ao=NULL, Bo=NULL, W, e;
  pari_timer ti;
  double rf = p ==3 ? 1.2 : .9;
  long n = degpol(T0), r = (long) sqrt(n*rf);
  GEN T;
  long r2 = 3*r/2;
  long nbi = itos(ffsumnbirred(utoipos(p), r)), nbrel=nbi*5/4;
  if (DEBUGLEVEL)
  {
    err_printf("Coppersmith: Parameters r=%ld r2=%ld\n", r,r2);
    err_printf("Coppersmith: Size FB=%ld rel. needed=%ld\n", nbi, nbrel);
    timer_start(&ti);
  }
  T = smallirred_Flx(p,n,get_Flx_var(T0), pi);
  S = Flx_ffisom(T0,T,p);
  a = Flx_Flxq_eval_pre(a0, S, T, p, pi);
  b = Flx_Flxq_eval_pre(b0, S, T, p, pi);
  if (DEBUGLEVEL) timer_printf(&ti,"model change");
  M = Flxq_log_Coppersmith(nbrel, r, T, p, pi);
  if (DEBUGLEVEL) timer_printf(&ti,"relations");
  W = check_kernel(r, M, nbi, 3*upowuu(p,r), T, p, pi, m);
  timer_start(&ti);
  Ao = Flxq_log_Coppersmith_rec(W, r2, a, r, T, p, pi, m);
  if (DEBUGLEVEL) timer_printf(&ti,"smooth element");
  Bo = Flxq_log_Coppersmith_rec(W, r2, b, r, T, p, pi, m);
  if (DEBUGLEVEL) timer_printf(&ti,"smooth generator");
  e = Fp_div(Ao, Bo, m);
  if (!Flx_equal(Flxq_pow_pre(b0,e,T0,p,pi), a0)) pari_err_BUG("Flxq_log");
  return gerepileupto(av, e);
}

GEN
Flxq_log_index(GEN a, GEN b, GEN m, GEN T, ulong p)
{
  long d = get_Flx_degree(T);
  if (p==3 || (p==5 && d>41))
    return Flxq_log_index_Coppersmith(a, b, m, T, p);
  else    return Flxq_log_index_cubic(a, b, m, T, p);
}

int
Flxq_log_use_index(GEN m, GEN T, ulong p)
{
  long d = get_Flx_degree(T);
  if (p==3 || (p==5 && d>41))
    return 1;
  else if (d<=4 || d==6)
    return 0;
  else
    return Flxq_log_use_index_cubic(m, T, p);
}
