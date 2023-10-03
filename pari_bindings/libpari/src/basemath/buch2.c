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

#define DEBUGLEVEL DEBUGLEVEL_bnf

/*******************************************************************/
/*                                                                 */
/*         CLASS GROUP AND REGULATOR (McCURLEY, BUCHMANN)          */
/*                    GENERAL NUMBER FIELDS                        */
/*                                                                 */
/*******************************************************************/
/* get_random_ideal */
static const long RANDOM_BITS = 4;
/* Buchall */
static const long RELSUP = 5;
static const long FAIL_DIVISOR = 32;
static const long MINFAIL = 10;
/* small_norm */
static const long BNF_RELPID = 4;
static const long maxtry_FACT = 500;
/* rnd_rel */
static const long RND_REL_RELPID = 1;
/* random relations */
static const long MINSFB = 3;
static const long SFB_MAX = 3;
static const long DEPSIZESFBMULT = 16;
static const long DEPSFBDIV = 10;
/* add_rel_i */
static const ulong mod_p = 27449UL;
/* be_honest */
static const long maxtry_HONEST = 50;

typedef struct FACT {
    long pr, ex;
} FACT;

typedef struct subFB_t {
  GEN subFB;
  struct subFB_t *old;
} subFB_t;

/* a factor base contains only noninert primes
 * KC = # of P in factor base (p <= n, NP <= n2)
 * KC2= # of P assumed to generate class group (NP <= n2)
 *
 * KCZ = # of rational primes under ideals counted by KC
 * KCZ2= same for KC2 */

typedef struct FB_t {
  GEN FB; /* FB[i] = i-th rational prime used in factor base */
  GEN LP; /* vector of all prime ideals in FB, by increasing norm */
  GEN LV; /* LV[p] = vector of P|p, NP <= n2
            * isclone() is set for LV[p] iff all P|p are in FB
            * LV[i], i not prime or i > n2, is undefined! */
  GEN iLP; /* iLP[p] = i such that LV[p] = [LP[i],...] */
  GEN L_jid; /* indexes of "useful" prime ideals for rnd_rel */
  long KC, KCZ, KCZ2;
  GEN prodZ; /* product of the primes in KCZ*/
  GEN subFB; /* LP o subFB =  part of FB used to build random relations */
  int sfb_chg; /* need to change subFB ? */
  GEN perm; /* permutation of LP used to represent relations [updated by
               hnfspec/hnfadd: dense rows come first] */
  GEN idealperm; /* permutation of ideals under field automorphisms */
  GEN minidx; /* minidx[i] min ideal in orbit of LP[i] under field autom */
  subFB_t *allsubFB; /* all subFB's used */
  GEN embperm; /* permutations of the complex embeddings */
  long MAXDEPSIZESFB; /* # trials before increasing subFB */
  long MAXDEPSFB; /* MAXDEPSIZESFB / DEPSFBDIV, # trials befor rotating subFB */
  double ballvol;
} FB_t;

enum { sfb_CHANGE = 1, sfb_INCREASE = 2 };

typedef struct REL_t {
  GEN R; /* relation vector as t_VECSMALL; clone */
  long nz; /* index of first nonzero elt in R (hash) */
  GEN m; /* pseudo-minimum yielding the relation; clone */
  long relorig; /* relation this one is an image of */
  long relaut; /* automorphim used to compute this relation from the original */
  GEN emb; /* archimedean embeddings */
  GEN junk[2]; /*make sure sizeof(struct) is a power of two.*/
} REL_t;

typedef struct RELCACHE_t {
  REL_t *chk; /* last checkpoint */
  REL_t *base; /* first rel found */
  REL_t *last; /* last rel found so far */
  REL_t *end; /* target for last relation. base <= last <= end */
  size_t len; /* number of rels pre-allocated in base */
  long relsup; /* how many linearly dependent relations we allow */
  GEN basis; /* mod p basis (generating family actually) */
  ulong missing; /* missing vectors in generating family above */
} RELCACHE_t;

typedef struct FP_t {
  double **q;
  GEN x;
  double *y;
  double *z;
  double *v;
} FP_t;

typedef struct RNDREL_t {
  long jid;
  GEN ex;
} RNDREL_t;

static void
wr_rel(GEN e)
{
  long i, l = lg(e);
  for (i = 1; i < l; i++)
    if (e[i]) err_printf("%ld^%ld ",i,e[i]);
}
static void
dbg_newrel(RELCACHE_t *cache)
{
  if (DEBUGLEVEL > 1)
  {
    err_printf("\n++++ cglob = %ld\nrel = ", cache->last - cache->base);
    wr_rel(cache->last->R);
    err_printf("\n");
  }
  else
    err_printf("%ld ", cache->last - cache->base);
}

static void
delete_cache(RELCACHE_t *M)
{
  REL_t *rel;
  for (rel = M->base+1; rel <= M->last; rel++)
  {
    gunclone(rel->R);
    if (rel->m) gunclone(rel->m);
  }
  pari_free((void*)M->base); M->base = NULL;
}

static void
delete_FB(FB_t *F)
{
  subFB_t *s, *sold;
  for (s = F->allsubFB; s; s = sold) { sold = s->old; pari_free(s); }
  gunclone(F->minidx);
  gunclone(F->idealperm);
}

static void
reallocate(RELCACHE_t *M, long len)
{
  M->len = len;
  if (!M->base)
    M->base = (REL_t*)pari_malloc((len+1) * sizeof(REL_t));
  else
  {
    size_t l = M->last - M->base, c = M->chk - M->base, e = M->end - M->base;
    pari_realloc_ip((void**)&M->base, (len+1) * sizeof(REL_t));
    M->last = M->base + l;
    M->chk  = M->base + c;
    M->end  = M->base + e;
  }
}

#define pr_get_smallp(pr) gel(pr,1)[2]

/* don't take P|p all other Q|p are already there */
static int
bad_subFB(FB_t *F, long t)
{
  GEN LP, P = gel(F->LP,t);
  long p = pr_get_smallp(P);
  LP = gel(F->LV,p);
  return (isclone(LP) && t == F->iLP[p] + lg(LP)-1);
}

static void
assign_subFB(FB_t *F, GEN yes, long iyes)
{
  long i, lv = sizeof(subFB_t) + iyes*sizeof(long); /* for struct + GEN */
  subFB_t *s = (subFB_t *)pari_malloc(lv);
  s->subFB = (GEN)&s[1];
  s->old = F->allsubFB; F->allsubFB = s;
  for (i = 0; i < iyes; i++) s->subFB[i] = yes[i];
  F->subFB = s->subFB;
  F->MAXDEPSIZESFB = (iyes-1) * DEPSIZESFBMULT;
  F->MAXDEPSFB = F->MAXDEPSIZESFB / DEPSFBDIV;
}

/* Determine the permutation of the ideals made by each field automorphism */
static GEN
FB_aut_perm(FB_t *F, GEN auts, GEN cyclic)
{
  long i, j, m, KC = F->KC, nauts = lg(auts)-1;
  GEN minidx, perm = zero_Flm_copy(KC, nauts);

  if (!nauts) { F->minidx = gclone(identity_zv(KC)); return cgetg(1,t_MAT); }
  minidx = zero_Flv(KC);
  for (m = 1; m < lg(cyclic); m++)
  {
    GEN thiscyc = gel(cyclic, m);
    long k0 = thiscyc[1];
    GEN aut = gel(auts, k0), permk0 = gel(perm, k0), ppermk;
    i = 1;
    while (i <= KC)
    {
      pari_sp av2 = avma;
      GEN seen = zero_Flv(KC), P = gel(F->LP, i);
      long imin = i, p, f, l;
      p = pr_get_smallp(P);
      f = pr_get_f(P);
      do
      {
        if (++i > KC) break;
        P = gel(F->LP, i);
      }
      while (p == pr_get_smallp(P) && f == pr_get_f(P));
      for (j = imin; j < i; j++)
      {
        GEN img = ZM_ZC_mul(aut, pr_get_gen(gel(F->LP, j)));
        for (l = imin; l < i; l++)
          if (!seen[l] && ZC_prdvd(img, gel(F->LP, l)))
          {
            seen[l] = 1; permk0[j] = l; break;
          }
      }
      set_avma(av2);
    }
    for (ppermk = permk0, i = 2; i < lg(thiscyc); i++)
    {
      GEN permk = gel(perm, thiscyc[i]);
      for (j = 1; j <= KC; j++) permk[j] = permk0[ppermk[j]];
      ppermk = permk;
    }
  }
  for (j = 1; j <= KC; j++)
  {
    if (minidx[j]) continue;
    minidx[j] = j;
    for (i = 1; i <= nauts; i++) minidx[coeff(perm, j, i)] = j;
  }
  F->minidx = gclone(minidx); return perm;
}

/* set subFB.
 * Fill F->perm (if != NULL): primes ideals sorted by increasing norm (except
 * the ones in subFB come first [dense rows for hnfspec]) */
static void
subFBgen(FB_t *F, GEN auts, GEN cyclic, double PROD, long minsFB)
{
  GEN y, perm, yes, no;
  long i, j, k, iyes, ino, lv = F->KC + 1;
  double prod;
  pari_sp av;

  F->LP   = cgetg(lv, t_VEC);
  F->L_jid = F->perm = cgetg(lv, t_VECSMALL);
  av = avma;
  y = cgetg(lv,t_COL); /* Norm P */
  for (k=0, i=1; i <= F->KCZ; i++)
  {
    GEN LP = gel(F->LV,F->FB[i]);
    long l = lg(LP);
    for (j = 1; j < l; j++)
    {
      GEN P = gel(LP,j);
      k++;
      gel(y,k) = pr_norm(P);
      gel(F->LP,k) = P;
    }
  }
  /* perm sorts LP by increasing norm */
  perm = indexsort(y);
  no  = cgetg(lv, t_VECSMALL); ino  = 1;
  yes = cgetg(lv, t_VECSMALL); iyes = 1;
  prod = 1.0;
  for (i = 1; i < lv; i++)
  {
    long t = perm[i];
    if (bad_subFB(F, t)) { no[ino++] = t; continue; }

    yes[iyes++] = t;
    prod *= (double)itos(gel(y,t));
    if (iyes > minsFB && prod > PROD) break;
  }
  setlg(yes, iyes);
  for (j=1; j<iyes; j++)     F->perm[j] = yes[j];
  for (i=1; i<ino; i++, j++) F->perm[j] =  no[i];
  for (   ; j<lv; j++)       F->perm[j] =  perm[j];
  F->allsubFB = NULL;
  F->idealperm = gclone(FB_aut_perm(F, auts, cyclic));
  if (iyes) assign_subFB(F, yes, iyes);
  set_avma(av);
}
static int
subFB_change(FB_t *F)
{
  long i, iyes, minsFB, lv = F->KC + 1, l = lg(F->subFB)-1;
  pari_sp av = avma;
  GEN yes, L_jid = F->L_jid, present = zero_zv(lv-1);

  switch (F->sfb_chg)
  {
    case sfb_INCREASE: minsFB = l + 1; break;
    default: minsFB = l; break;
  }

  yes = cgetg(minsFB+1, t_VECSMALL); iyes = 1;
  if (L_jid)
  {
    for (i = 1; i < lg(L_jid); i++)
    {
      long l = L_jid[i];
      yes[iyes++] = l;
      present[l] = 1;
      if (iyes > minsFB) break;
    }
  }
  else i = 1;
  if (iyes <= minsFB)
  {
    for ( ; i < lv; i++)
    {
      long l = F->perm[i];
      if (present[l]) continue;
      yes[iyes++] = l;
      if (iyes > minsFB) break;
    }
    if (i == lv) return 0;
  }
  if (zv_equal(F->subFB, yes))
  {
    if (DEBUGLEVEL) err_printf("\n*** NOT Changing sub factor base\n");
  }
  else
  {
    if (DEBUGLEVEL) err_printf("\n*** Changing sub factor base\n");
    assign_subFB(F, yes, iyes);
  }
  F->sfb_chg = 0; return gc_bool(av, 1);
}

/* make sure enough room to store n more relations */
static void
pre_allocate(RELCACHE_t *cache, size_t n)
{
  size_t len = (cache->last - cache->base) + n;
  if (len >= cache->len) reallocate(cache, len << 1);
}

void
init_GRHcheck(GRHcheck_t *S, long N, long R1, double LOGD)
{
  const double c1 = M_PI*M_PI/2;
  const double c2 = 3.663862376709;
  const double c3 = 3.801387092431; /* Euler + log(8*Pi)*/
  S->clone = 0;
  S->cN = R1*c2 + N*c1;
  S->cD = LOGD - N*c3 - R1*M_PI/2;
  S->maxprimes = 16000; /* sufficient for LIMC=176081*/
  S->primes = (GRHprime_t*)pari_malloc(S->maxprimes*sizeof(*S->primes));
  S->nprimes = 0;
  S->limp = 0;
  u_forprime_init(&S->P, 2, ULONG_MAX);
}

void
free_GRHcheck(GRHcheck_t *S)
{
  if (S->clone)
  {
    long i = S->nprimes;
    GRHprime_t *pr;
    for (pr = S->primes, i = S->nprimes; i > 0; pr++, i--) gunclone(pr->dec);
  }
  pari_free(S->primes);
}

int
GRHok(GRHcheck_t *S, double L, double SA, double SB)
{
  return (S->cD + (S->cN + 2*SB) / L - 2*SA < -1e-8);
}

/* Return factorization pattern of p: [f,n], where n[i] primes of
 * residue degree f[i] */
static GEN
get_fs(GEN nf, GEN P, GEN index, ulong p)
{
  long j, k, f, n, l;
  GEN fs, ns;

  if (umodiu(index, p))
  { /* easy case: p does not divide index */
    GEN F = Flx_degfact(ZX_to_Flx(P,p), p);
    fs = gel(F,1); l = lg(fs);
  }
  else
  {
    GEN F = idealprimedec(nf, utoipos(p));
    l = lg(F);
    fs = cgetg(l, t_VECSMALL);
    for (j = 1; j < l; j++) fs[j] = pr_get_f(gel(F,j));
  }
  ns = cgetg(l, t_VECSMALL);
  f = fs[1]; n = 1;
  for (j = 2, k = 1; j < l; j++)
    if (fs[j] == f)
      n++;
    else
    {
      ns[k] = n; fs[k] = f; k++;
      f = fs[j]; n = 1;
    }
  ns[k] = n; fs[k] = f; k++;
  setlg(fs, k);
  setlg(ns, k); return mkvec2(fs,ns);
}

/* cache data for all rational primes up to the LIM */
static void
cache_prime_dec(GRHcheck_t *S, ulong LIM, GEN nf)
{
  pari_sp av = avma;
  GRHprime_t *pr;
  GEN index, P;
  double nb;

  if (S->limp >= LIM) return;
  S->clone = 1;
  nb = primepi_upper_bound((double)LIM); /* #{p <= LIM} <= nb */
  GRH_ensure(S, nb+1); /* room for one extra prime */
  P = nf_get_pol(nf);
  index = nf_get_index(nf);
  for (pr = S->primes + S->nprimes;;)
  {
    ulong p = u_forprime_next(&(S->P));
    pr->p = p;
    pr->logp = log((double)p);
    pr->dec = gclone(get_fs(nf, P, index, p));
    S->nprimes++;
    pr++;
    set_avma(av);
    /* store up to nextprime(LIM) included */
    if (p >= LIM) { S->limp = p; break; }
  }
}

static double
tailresback(long R1, long R2, double rK, long C, double C2, double C3, double r1K, double r2K, double logC, double logC2, double logC3)
{
  const double  rQ = 1.83787706641;
  const double r1Q = 1.98505372441;
  const double r2Q = 1.07991541347;
  return fabs((R1+R2-1)*(12*logC3+4*logC2-9*logC-6)/(2*C*logC3)
         + (rK-rQ)*(6*logC2 + 5*logC + 2)/(C*logC3)
         - R2*(6*logC2+11*logC+6)/(C2*logC2)
         - 2*(r1K-r1Q)*(3*logC2 + 4*logC + 2)/(C2*logC3)
         + (R1+R2-1)*(12*logC3+40*logC2+45*logC+18)/(6*C3*logC3)
         + (r2K-r2Q)*(2*logC2 + 3*logC + 2)/(C3*logC3));
}

static double
tailres(long R1, long R2, double al2K, double rKm, double rKM, double r1Km,
        double r1KM, double r2Km, double r2KM, double C, long i)
{
  /* C >= 3*2^i, lower bound for eint1(log(C)/2) */
  /* for(i=0,30,print(eint1(log(3*2^i)/2))) */
  static double tab[] = {
    0.50409264803,
    0.26205336997,
    0.14815491171,
    0.08770540561,
    0.05347651832,
    0.03328934284,
    0.02104510690,
    0.01346475900,
    0.00869778586,
    0.00566279855,
    0.00371111950,
    0.00244567837,
    0.00161948049,
    0.00107686891,
    0.00071868750,
    0.00048119961,
    0.00032312188,
    0.00021753772,
    0.00014679818,
    9.9272855581E-5,
    6.7263969995E-5,
    4.5656812967E-5,
    3.1041124593E-5,
    2.1136011590E-5,
    1.4411645381E-5,
    9.8393304088E-6,
    6.7257395409E-6,
    4.6025878272E-6,
    3.1529719271E-6,
    2.1620490021E-6,
    1.4839266071E-6
  };
  const double logC = log(C), logC2 = logC*logC, logC3 = logC*logC2;
  const double C2 = C*C, C3 = C*C2;
  double E1 = i >30? 0: tab[i];
  return al2K*((33*logC2+22*logC+8)/(8*logC3*sqrt(C))+15*E1/16)
    + maxdd(tailresback(rKm,r1KM,r2Km, C,C2,C3,R1,R2,logC,logC2,logC3),
            tailresback(rKM,r1Km,r2KM, C,C2,C3,R1,R2,logC,logC2,logC3))/2
    + ((R1+R2-1)*4*C+R2)*(C2+6*logC)/(4*C2*C2*logC2);
}

static long
primeneeded(long N, long R1, long R2, double LOGD)
{
  const double lim = 0.25; /* should be log(2)/2 == 0.34657... */
  const double al2K =  0.3526*LOGD - 0.8212*N + 4.5007;
  const double  rKm = -1.0155*LOGD + 2.1042*N - 8.3419;
  const double  rKM = -0.5   *LOGD + 1.2076*N + 1;
  const double r1Km = -       LOGD + 1.4150*N;
  const double r1KM = -       LOGD + 1.9851*N;
  const double r2Km = -       LOGD + 0.9151*N;
  const double r2KM = -       LOGD + 1.0800*N;
  long Cmin = 3, Cmax = 3, i = 0;
  while (tailres(R1, R2, al2K, rKm, rKM, r1Km, r1KM, r2Km, r2KM, Cmax, i) > lim)
  {
    Cmin = Cmax;
    Cmax *= 2;
    i++;
  }
  i--;
  while (Cmax - Cmin > 1)
  {
    long t = (Cmin + Cmax)/2;
    if (tailres(R1, R2, al2K, rKm, rKM, r1Km, r1KM, r2Km, r2KM, t, i) > lim)
      Cmin = t;
    else
      Cmax = t;
  }
  return Cmax;
}

/* ~ 1 / Res(s = 1, zeta_K) */
static GEN
compute_invres(GRHcheck_t *S, long LIMC)
{
  pari_sp av = avma;
  double loginvres = 0.;
  GRHprime_t *pr;
  long i;
  double logLIMC = log((double)LIMC);
  double logLIMC2 = logLIMC*logLIMC, denc;
  double c0, c1, c2;
  denc = 1/(pow((double)LIMC, 3.) * logLIMC * logLIMC2);
  c2 = (    logLIMC2 + 3 * logLIMC / 2 + 1) * denc;
  denc *= LIMC;
  c1 = (3 * logLIMC2 + 4 * logLIMC     + 2) * denc;
  denc *= LIMC;
  c0 = (3 * logLIMC2 + 5 * logLIMC / 2 + 1) * denc;
  for (pr = S->primes, i = S->nprimes; i > 0; pr++, i--)
  {
    GEN dec, fs, ns;
    long addpsi;
    double addpsi1, addpsi2;
    double logp = pr->logp, NPk;
    long j, k, limp = logLIMC/logp;
    ulong p = pr->p, p2 = p*p;
    if (limp < 1) break;
    dec = pr->dec;
    fs = gel(dec, 1); ns = gel(dec, 2);
    loginvres += 1./p;
    /* NB: limp = 1 nearly always and limp > 2 for very few primes */
    for (k=2, NPk = p; k <= limp; k++) { NPk *= p; loginvres += 1/(k * NPk); }
    addpsi = limp;
    addpsi1 = p *(pow((double)p , (double)limp)-1)/(p -1);
    addpsi2 = p2*(pow((double)p2, (double)limp)-1)/(p2-1);
    j = lg(fs);
    while (--j > 0)
    {
      long f, nb, kmax;
      double NP, NP2, addinvres;
      f = fs[j]; if (f > limp) continue;
      nb = ns[j];
      NP = pow((double)p, (double)f);
      addinvres = 1/NP;
      kmax = limp / f;
      for (k=2, NPk = NP; k <= kmax; k++) { NPk *= NP; addinvres += 1/(k*NPk); }
      NP2 = NP*NP;
      loginvres -= nb * addinvres;
      addpsi -= nb * f * kmax;
      addpsi1 -= nb*(f*NP *(pow(NP ,(double)kmax)-1)/(NP -1));
      addpsi2 -= nb*(f*NP2*(pow(NP2,(double)kmax)-1)/(NP2-1));
    }
    loginvres -= (addpsi*c0 - addpsi1*c1 + addpsi2*c2)*logp;
  }
  return gerepileuptoleaf(av, mpexp(dbltor(loginvres)));
}

static long
nthideal(GRHcheck_t *S, GEN nf, long n)
{
  pari_sp av = avma;
  GEN P = nf_get_pol(nf);
  ulong p = 0, *vecN = (ulong*)const_vecsmall(n, LONG_MAX);
  long i, N = poldegree(P, -1);
  for (i = 0; ; i++)
  {
    GRHprime_t *pr;
    GEN fs;
    cache_prime_dec(S, p+1, nf);
    pr = S->primes + i;
    fs = gel(pr->dec, 1);
    p = pr->p;
    if (fs[1] != N)
    {
      GEN ns = gel(pr->dec, 2);
      long k, l, j = lg(fs);
      while (--j > 0)
      {
        ulong NP = upowuu(p, fs[j]);
        long nf;
        if (!NP) continue;
        for (k = 1; k <= n; k++) if (vecN[k] > NP) break;
        if (k > n) continue;
        /* vecN[k] <= NP */
        nf = ns[j]; /*#{primes of norme NP} = nf, insert them here*/
        for (l = k+nf; l <= n; l++) vecN[l] = vecN[l-nf];
        for (l = 0; l < nf && k+l <= n; l++) vecN[k+l] = NP;
        while (l <= k) vecN[l++] = NP;
      }
    }
    if (p > vecN[n]) break;
  }
  return gc_long(av, vecN[n]);
}

/* volume of unit ball in R^n: \pi^{n/2} / \Gamma(n/2 + 1) */
static double
ballvol(long n)
{
  double v = odd(n)? 2: 1;
  for (; n > 1; n -= 2) v *= (2 * M_PI) / n;
  return v;
}

/* Compute FB, LV, iLP + KC*. Reset perm
 * C2: bound for norm of tested prime ideals (includes be_honest())
 * C1: bound for p, such that P|p (NP <= C2) used to build relations */
static void
FBgen(FB_t *F, GEN nf, long N, ulong C1, ulong C2, GRHcheck_t *S)
{
  GRHprime_t *pr;
  long i, ip;
  GEN prim;
  const double L = log((double)C2 + 0.5);

  cache_prime_dec(S, C2, nf);
  pr = S->primes;
  F->sfb_chg = 0;
  F->FB  = cgetg(C2+1, t_VECSMALL);
  F->iLP = cgetg(C2+1, t_VECSMALL);
  F->LV = zerovec(C2);

  prim = icopy(gen_1);
  i = ip = 0;
  F->KC = F->KCZ = 0;
  for (;; pr++) /* p <= C2 */
  {
    ulong p = pr->p;
    long k, l, m;
    GEN LP, nb, f;

    if (!F->KC && p > C1) { F->KCZ = i; F->KC = ip; }
    if (p > C2) break;

    if (DEBUGLEVEL>1) err_printf(" %ld",p);

    f = gel(pr->dec, 1); nb = gel(pr->dec, 2);
    if (f[1] == N)
    {
      if (p == C2) break;
      continue; /* p inert */
    }
    l = (long)(L/pr->logp); /* p^f <= C2  <=> f <= l */
    for (k=0, m=1; m < lg(f) && f[m]<=l; m++) k += nb[m];
    if (!k)
    { /* too inert to appear in FB */
      if (p == C2) break;
      continue;
    }
    prim[2] = p; LP = idealprimedec_limit_f(nf,prim, l);
    /* keep noninert ideals with Norm <= C2 */
    if (m == lg(f)) setisclone(LP); /* flag it: all prime divisors in FB */
    F->FB[++i]= p;
    gel(F->LV,p) = LP;
    F->iLP[p] = ip; ip += k;
    if (p == C2) break;
  }
  if (!F->KC) { F->KCZ = i; F->KC = ip; }
  /* Note F->KC > 0 otherwise GRHchk is false */
  setlg(F->FB, F->KCZ+1); F->KCZ2 = i;
  F->prodZ = zv_prod_Z(F->FB);
  if (DEBUGLEVEL>1)
  {
    err_printf("\n");
    if (DEBUGLEVEL>6)
    {
      err_printf("########## FACTORBASE ##########\n\n");
      err_printf("KC2=%ld, KC=%ld, KCZ=%ld, KCZ2=%ld\n",
                  ip, F->KC, F->KCZ, F->KCZ2);
      for (i=1; i<=F->KCZ; i++) err_printf("++ LV[%ld] = %Ps",i,gel(F->LV,F->FB[i]));
    }
  }
  F->perm = NULL; F->L_jid = NULL;
  F->ballvol = ballvol(nf_get_degree(nf));
}

static int
GRHchk(GEN nf, GRHcheck_t *S, ulong LIMC)
{
  double logC = log((double)LIMC), SA = 0, SB = 0;
  GRHprime_t *pr = S->primes;

  cache_prime_dec(S, LIMC, nf);
  for (pr = S->primes;; pr++)
  {
    ulong p = pr->p;
    GEN dec, fs, ns;
    double logCslogp;
    long j;

    if (p > LIMC) break;
    dec = pr->dec; fs = gel(dec, 1); ns = gel(dec,2);
    logCslogp = logC/pr->logp;
    for (j = 1; j < lg(fs); j++)
    {
      long f = fs[j], M, nb;
      double logNP, q, A, B;
      if (f > logCslogp) break;
      logNP = f * pr->logp;
      q = 1/sqrt((double)upowuu(p, f));
      A = logNP * q; B = logNP * A; M = (long)(logCslogp/f);
      if (M > 1)
      {
        double inv1_q = 1 / (1-q);
        A *= (1 - pow(q, (double)M)) * inv1_q;
        B *= (1 - pow(q, (double)M)*(M+1 - M*q)) * inv1_q * inv1_q;
      }
      nb = ns[j];
      SA += nb * A;
      SB += nb * B;
    }
    if (p == LIMC) break;
  }
  return GRHok(S, logC, SA, SB);
}

/*  SMOOTH IDEALS */
static void
store(long i, long e, FACT *fact)
{
  ++fact[0].pr;
  fact[fact[0].pr].pr = i; /* index */
  fact[fact[0].pr].ex = e; /* exponent */
}

/* divide out x by all P|p, where x as in can_factor().  k = v_p(Nx) */
static int
divide_p_elt(GEN LP, long ip, long k, GEN m, FACT *fact)
{
  long j, l = lg(LP);
  for (j=1; j<l; j++)
  {
    GEN P = gel(LP,j);
    long v = ZC_nfval(m, P);
    if (!v) continue;
    store(ip + j, v, fact); /* v = v_P(m) > 0 */
    k -= v * pr_get_f(P);
    if (!k) return 1;
  }
  return 0;
}
static int
divide_p_id(GEN LP, long ip, long k, GEN nf, GEN I, FACT *fact)
{
  long j, l = lg(LP);
  for (j=1; j<l; j++)
  {
    GEN P = gel(LP,j);
    long v = idealval(nf,I, P);
    if (!v) continue;
    store(ip + j, v, fact); /* v = v_P(I) > 0 */
    k -= v * pr_get_f(P);
    if (!k) return 1;
  }
  return 0;
}
static int
divide_p_quo(GEN LP, long ip, long k, GEN nf, GEN I, GEN m, FACT *fact)
{
  long j, l = lg(LP);
  for (j=1; j<l; j++)
  {
    GEN P = gel(LP,j);
    long v = ZC_nfval(m, P);
    if (!v) continue;
    v -= idealval(nf,I, P);
    if (!v) continue;
    store(ip + j, v, fact); /* v = v_P(m / I) > 0 */
    k -= v * pr_get_f(P);
    if (!k) return 1;
  }
  return 0;
}

/* |*N| != 0 is the norm of a primitive ideal, in particular not divisible by
 * any inert prime. Is |*N| a smooth rational integer wrt F ?
 */
static int
Z_issmooth_prod(GEN N, GEN P)
{
  P = gcdii(P,N);
  while (!is_pm1(P))
  {
    N = diviiexact(N, P);
    P = gcdii(N, P);
  }
  return is_pm1(N);
}

static int
divide_p(FB_t *F, long p, long k, GEN nf, GEN I, GEN m, FACT *fact)
{
  GEN LP = gel(F->LV,p);
  long ip = F->iLP[p];
  if (!m) return divide_p_id (LP,ip,k,nf,I,fact);
  if (!I) return divide_p_elt(LP,ip,k,m,fact);
  return divide_p_quo(LP,ip,k,nf,I,m,fact);
}

/* Let x = m if I == NULL,
 *         I if m == NULL,
 *         m/I otherwise.
 * Can we factor the integral primitive ideal x ? |N| = Norm x > 0 */
static long
can_factor(FB_t *F, GEN nf, GEN I, GEN m, GEN N, FACT *fact)
{
  GEN f, p, e;
  long i, l;
  fact[0].pr = 0;
  if (is_pm1(N)) return 1;
  if (!Z_issmooth_prod(N, F->prodZ)) return 0;
  f = absZ_factor(N); p = gel(f,1); e = gel(f,2); l = lg(p);
  for (i = 1; i < l; i++)
    if (!divide_p(F, itou(gel(p,i)), itou(gel(e,i)), nf, I, m, fact))
    {
      if (DEBUGLEVEL > 1) err_printf(".");
      return 0;
    }
  return 1;
}

/* can we factor m/I ? [m in I from idealpseudomin_nonscalar], NI = norm I */
static long
factorgen(FB_t *F, GEN nf, GEN I, GEN NI, GEN m, FACT *fact)
{
  long e, r1 = nf_get_r1(nf);
  GEN M = nf_get_M(nf);
  GEN N = divri(embed_norm(RgM_RgC_mul(M,m), r1), NI); /* ~ N(m/I) */
  N = grndtoi(N, &e);
  if (e > -32)
  {
    if (DEBUGLEVEL > 1) err_printf("+");
    return 0;
  }
  return can_factor(F, nf, I, m, N, fact);
}

/*  FUNDAMENTAL UNITS */

/* a, y real. Return  (Re(x) + a) + I * (Im(x) % y) */
static GEN
addRe_modIm(GEN x, GEN a, GEN y, GEN iy)
{
  GEN z;
  if (typ(x) == t_COMPLEX)
  {
    GEN re, im = modRr_i(gel(x,2), y, iy);
    if (!im) return NULL;
    re = gadd(gel(x,1), a);
    z = gequal0(im)? re: mkcomplex(re, im);
  }
  else
    z = gadd(x, a);
  return z;
}
static GEN
modIm(GEN x, GEN y, GEN iy)
{
  if (typ(x) == t_COMPLEX)
  {
    GEN im = modRr_i(gel(x,2), y, iy);
    if (!im) return NULL;
    x = gequal0(im)? gel(x,1): mkcomplex(gel(x,1), im);
  }
  return x;
}

/* clean archimedean components. ipi = 2^n / pi (n arbitrary); its
 * exponent may be modified */
static GEN
cleanarch(GEN x, long N, GEN ipi, long prec)
{
  long i, l, R1, RU;
  GEN s, y = cgetg_copy(x, &l);

  if (!ipi) ipi = invr(mppi(prec));
  if (typ(x) == t_MAT)
  {
    for (i = 1; i < l; i++)
      if (!(gel(y,i) = cleanarch(gel(x,i), N, ipi, prec))) return NULL;
    return y;
  }
  RU = l-1; R1 = (RU<<1) - N;
  s = gdivgs(RgV_sum(real_i(x)), -N); /* -log |norm(x)| / N */
  i = 1;
  if (R1)
  {
    GEN pi2 = Pi2n(1, prec);
    setexpo(ipi, -3); /* 1/(2pi) */
    for (; i <= R1; i++)
      if (!(gel(y,i) = addRe_modIm(gel(x,i), s, pi2, ipi))) return NULL;
  }
  if (i <= RU)
  {
    GEN pi4 = Pi2n(2, prec), s2 = gmul2n(s, 1);
    setexpo(ipi, -4); /* 1/(4pi) */
    for (; i <= RU; i++)
      if (!(gel(y,i) = addRe_modIm(gel(x,i), s2, pi4, ipi))) return NULL;
  }
  return y;
}
GEN
nf_cxlog_normalize(GEN nf, GEN x, long prec)
{
  long N = nf_get_degree(nf);
  return cleanarch(x, N, NULL, prec);
}

/* clean unit archimedean components. ipi = 2^n / pi (n arbitrary); its
 * exponent may be modified */
static GEN
cleanarchunit(GEN x, long N, GEN ipi, long prec)
{
  long i, l, R1, RU;
  GEN y = cgetg_copy(x, &l);

  if (!ipi) ipi = invr(mppi(prec));
  if (typ(x) == t_MAT)
  {
    for (i = 1; i < l; i++)
      if (!(gel(y,i) = cleanarchunit(gel(x,i), N, ipi, prec))) return NULL;
    return y;
  }
  if (gexpo(RgV_sum(real_i(x))) > -10) return NULL;
  RU = l-1; R1 = (RU<<1) - N;
  i = 1;
  if (R1)
  {
    GEN pi2 = Pi2n(1, prec);
    setexpo(ipi, -3); /* 1/(2pi) */
    for (; i <= R1; i++)
      if (!(gel(y,i) = modIm(gel(x,i), pi2, ipi))) return NULL;
  }
  if (i <= RU)
  {
    GEN pi4 = Pi2n(2, prec);
    setexpo(ipi, -4); /* 1/(4pi) */
    for (; i <= RU; i++)
      if (!(gel(y,i) = modIm(gel(x,i), pi4, ipi))) return NULL;
  }
  return y;
}

static GEN
not_given(long reason)
{
  if (DEBUGLEVEL)
    switch(reason)
    {
      case fupb_LARGE:
        pari_warn(warner,"fundamental units too large, not given");
        break;
      case fupb_PRECI:
        pari_warn(warner,"insufficient precision for fundamental units, not given");
        break;
    }
  return NULL;
}

/* check whether exp(x) will 1) get too big (real(x) large), 2) require
 * large accuracy for argument reduction (imag(x) large) */
static long
expbitprec(GEN x, long *e)
{
  GEN re, im;
  if (typ(x) != t_COMPLEX) re = x;
  else
  {
    im = gel(x,2); *e = maxss(*e, expo(im) + 5 - bit_prec(im));
    re = gel(x,1);
  }
  return (expo(re) <= 20);

}
static long
RgC_expbitprec(GEN x)
{
  long l = lg(x), i, e = - (long)HIGHEXPOBIT;
  for (i = 1; i < l; i++)
    if (!expbitprec(gel(x,i), &e)) return LONG_MAX;
  return e;
}
static long
RgM_expbitprec(GEN x)
{
  long i, j, I, J, e = - (long)HIGHEXPOBIT;
  RgM_dimensions(x, &I,&J);
  for (j = 1; j <= J; j++)
    for (i = 1; i <= I; i++)
      if (!expbitprec(gcoeff(x,i,j), &e)) return LONG_MAX;
  return e;
}

static GEN
FlxqX_chinese_unit(GEN X, GEN U, GEN invzk, GEN D, GEN T, ulong p)
{
  long i, lU = lg(U), lX = lg(X), d = lg(invzk)-1;
  GEN M = cgetg(lU, t_MAT);
  if (D)
  {
    D = Flv_inv(D, p);
    for (i = 1; i < lX; i++)
      if (uel(D, i) != 1)
        gel(X,i) = Flx_Fl_mul(gel(X,i), uel(D,i), p);
  }
  for (i = 1; i < lU; i++)
  {
    GEN H = FlxqV_factorback(X, gel(U, i), T, p);
    gel(M, i) = Flm_Flc_mul(invzk, Flx_to_Flv(H, d), p);
  }
  return M;
}

static GEN
chinese_unit_slice(GEN A, GEN U, GEN B, GEN D, GEN C, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1, v = varn(C);
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1);
    GEN a = ZXV_to_FlxV(A, p), b = ZM_to_Flm(B, p), c = ZX_to_Flx(C, p);
    GEN d = D ? ZV_to_Flv(D, p): NULL;
    GEN Hp = FlxqX_chinese_unit(a, U, b, d, c, p);
    H = gerepileupto(av, Flm_to_ZM(Hp));
    *mod = utoi(p);
    return H;
  }
  T = ZV_producttree(P);
  A = ZXC_nv_mod_tree(A, P, T, v);
  B = ZM_nv_mod_tree(B, P, T);
  D = D ? ZV_nv_mod_tree(D, P, T): NULL;
  C = ZX_nv_mod_tree(C, P, T);

  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i), c = gel(C,i), d = D ? gel(D,i): NULL;
    gel(H,i) = FlxqX_chinese_unit(a, U, b, d, c, p);
  }
  H = nmV_chinese_center_tree_seq(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
chinese_unit_worker(GEN P, GEN A, GEN U, GEN B, GEN D, GEN C)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = chinese_unit_slice(A, U, B, isintzero(D) ? NULL: D, C, P, &gel(V,2));
  return V;
}

/* Let x = \prod X[i]^E[i] = u, return u.
 * If dX != NULL, X[i] = nX[i] / dX[i] where nX[i] is a ZX, dX[i] in Z */
static GEN
chinese_unit(GEN nf, GEN nX, GEN dX, GEN U, ulong bnd)
{
  pari_sp av = avma;
  GEN f = nf_get_index(nf), T = nf_get_pol(nf), invzk = nf_get_invzk(nf);
  GEN H, mod;
  forprime_t S;
  GEN worker = snm_closure(is_entry("_chinese_unit_worker"),
               mkcol5(nX, U, invzk, dX? dX: gen_0, T));
  init_modular_big(&S);
  H = gen_crt("chinese_units", worker, &S, f, bnd, 0, &mod, nmV_chinese_center, FpM_center);
  settyp(H, t_VEC); return gerepilecopy(av, H);
}

/* *pE a ZM */
static void
ZM_remove_unused(GEN *pE, GEN *pX)
{
  long j, k, l = lg(*pX);
  GEN E = *pE, v = cgetg(l, t_VECSMALL);
  for (j = k = 1; j < l; j++)
    if (!ZMrow_equal0(E, j)) v[k++] = j;
  if (k < l)
  {
    setlg(v, k);
    *pX = vecpermute(*pX,v);
    *pE = rowpermute(E,v);
  }
}

/* s = -log|norm(x)|/N */
static GEN
fixarch(GEN x, GEN s, long R1)
{
  long i, l;
  GEN y = cgetg_copy(x, &l);
  for (i = 1; i <= R1; i++) gel(y,i) = gadd(s, gel(x,i));
  for (     ; i <   l; i++) gel(y,i) = gadd(s, gmul2n(gel(x,i),-1));
  return y;
}

static GEN
getfu(GEN nf, GEN *ptA, GEN *ptU, long prec)
{
  GEN U, y, matep, A, T = nf_get_pol(nf), M = nf_get_M(nf);
  long e, j, R1, RU, N = degpol(T);

  R1 = nf_get_r1(nf); RU = (N+R1) >> 1;
  if (RU == 1) return cgetg(1,t_VEC);

  A = *ptA;
  matep = cgetg(RU,t_MAT);
  for (j = 1; j < RU; j++)
  {
    GEN Aj = gel(A,j), s = gdivgs(RgV_sum(real_i(Aj)), -N);
    gel(matep,j) = fixarch(Aj, s, R1);
  }
  U = lll(real_i(matep));
  if (lg(U) < RU) return not_given(fupb_PRECI);
  if (ptU) { *ptU = U; *ptA = A = RgM_ZM_mul(A,U); }
  y = RgM_ZM_mul(matep,U);
  e = RgM_expbitprec(y);
  if (e >= 0) return not_given(e == LONG_MAX? fupb_LARGE: fupb_PRECI);
  if (prec <= 0) prec = gprecision(A);
  y = RgM_solve_realimag(M, gexp(y,prec));
  if (!y) return not_given(fupb_PRECI);
  y = grndtoi(y, &e); if (e >= 0) return not_given(fupb_PRECI);
  settyp(y, t_VEC);

  if (!ptU) *ptA = A = RgM_ZM_mul(A, U);
  for (j = 1; j < RU; j++)
  { /* y[i] are hopefully unit generators. Normalize: smallest T2 norm */
    GEN u = gel(y,j), v = zk_inv(nf, u);
    if (!v || !is_pm1(Q_denom(v)) || ZV_isscalar(u))
      return not_given(fupb_PRECI);
    if (gcmp(RgC_fpnorml2(v,DEFAULTPREC), RgC_fpnorml2(u,DEFAULTPREC)) < 0)
    {
      gel(A,j) = RgC_neg(gel(A,j));
      if (ptU) gel(U,j) = ZC_neg(gel(U,j));
      u = v;
    }
    gel(y,j) = nf_to_scalar_or_alg(nf, u);
  }
  return y;
}

static void
err_units() { pari_err_PREC("makeunits [cannot get units, use bnfinit(,1)]"); }

/* bound for log2 |sigma(u)|, sigma complex embedding, u fundamental unit
 * attached to bnf_get_logfu */
static double
log2fubound(GEN bnf)
{
  GEN LU = bnf_get_logfu(bnf);
  long i, j, l = lg(LU), r1 = nf_get_r1(bnf_get_nf(bnf));
  double e = 0.0;
  for (j = 1; j < l; j++)
  {
    GEN u = gel(LU,j);
    for (i = 1; i <= r1; i++)
    {
      GEN E = real_i(gel(u,i));
      e = maxdd(e, gtodouble(E));
    }
    for (     ; i <= l; i++)
    {
      GEN E = real_i(gel(u,i));
      e = maxdd(e, gtodouble(E) / 2);
    }
  }
  return e / M_LN2;
}
/* bound for log2(|RgM_solve_realimag(M, y)|_oo / |y|_oo)*/
static double
log2Mbound(GEN nf)
{
  GEN G = nf_get_G(nf), D = nf_get_disc(nf);
  long r2 = nf_get_r2(nf), l = lg(G), i;
  double e, d = dbllog2(D)/2 - r2 * M_LN2; /* log2 |det(split_realimag(M))| */
  e = log2(nf_get_degree(nf));
  for (i = 2; i < l; i++) e += dbllog2(gnorml2(gel(G,i))); /* Hadamard bound */
  return e / 2 - d;
}

static GEN
vec_chinese_units(GEN bnf)
{
  GEN nf = bnf_get_nf(bnf), SUnits = bnf_get_sunits(bnf);
  double bnd = ceil(log2Mbound(nf) + log2fubound(bnf));
  GEN X, dX, Y, U, f = nf_get_index(nf);
  long j, l, v = nf_get_varn(nf);
  if (!SUnits) err_units(); /* no compact units */
  Y = gel(SUnits,1);
  U = gel(SUnits,2);
  ZM_remove_unused(&U, &Y); l = lg(Y); X = cgetg(l, t_VEC);
  if (is_pm1(f)) f = dX = NULL; else dX = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
  {
    GEN t = nf_to_scalar_or_alg(nf, gel(Y,j));
    if (f)
    {
      GEN den;
      t = Q_remove_denom(t, &den);
      gel(dX,j) = den ? den: gen_1;
    }
    gel(X,j) = typ(t) == t_INT? scalarpol_shallow(t,v): t;
  }
  if (dblexpo(bnd) >= BITS_IN_LONG)
    pari_err_OVERFLOW("vec_chinese_units [units too large]");
  return chinese_unit(nf, X, dX, U, (ulong)bnd);
}

static GEN
makeunits(GEN bnf)
{
  GEN nf = bnf_get_nf(bnf), fu = bnf_get_fu_nocheck(bnf);
  GEN tu = nf_to_scalar_or_basis(nf, bnf_get_tuU(bnf));
  fu = (typ(fu) == t_MAT)? vec_chinese_units(bnf): matalgtobasis(nf, fu);
  return vec_prepend(fu, tu);
}

/*******************************************************************/
/*                                                                 */
/*           PRINCIPAL IDEAL ALGORITHM (DISCRETE LOG)              */
/*                                                                 */
/*******************************************************************/

/* G: prime ideals, E: vector of nonnegative exponents.
 * C = possible extra prime (^1) or NULL
 * Return Norm (product) */
static GEN
get_norm_fact_primes(GEN G, GEN E, GEN C)
{
  pari_sp av=avma;
  GEN N = gen_1, P, p;
  long i, c = lg(E);
  for (i=1; i<c; i++)
  {
    GEN ex = gel(E,i);
    long s = signe(ex);
    if (!s) continue;

    P = gel(G,i); p = pr_get_p(P);
    N = mulii(N, powii(p, mului(pr_get_f(P), ex)));
  }
  if (C) N = mulii(N, pr_norm(C));
  return gerepileuptoint(av, N);
}

/* gen: HNF ideals */
static GEN
get_norm_fact(GEN gen, GEN ex, GEN *pd)
{
  long i, c = lg(ex);
  GEN d,N,I,e,n,ne,de;
  d = N = gen_1;
  for (i=1; i<c; i++)
    if (signe(gel(ex,i)))
    {
      I = gel(gen,i); e = gel(ex,i); n = ZM_det_triangular(I);
      ne = powii(n,e);
      de = equalii(n, gcoeff(I,1,1))? ne: powii(gcoeff(I,1,1), e);
      N = mulii(N, ne);
      d = mulii(d, de);
    }
  *pd = d; return N;
}

static GEN
get_pr_lists(GEN FB, long N, int list_pr)
{
  GEN pr, L;
  long i, l = lg(FB), p, pmax;

  pmax = 0;
  for (i=1; i<l; i++)
  {
    pr = gel(FB,i); p = pr_get_smallp(pr);
    if (p > pmax) pmax = p;
  }
  L = const_vec(pmax, NULL);
  if (list_pr)
  {
    for (i=1; i<l; i++)
    {
      pr = gel(FB,i); p = pr_get_smallp(pr);
      if (!L[p]) gel(L,p) = vectrunc_init(N+1);
      vectrunc_append(gel(L,p), pr);
    }
    for (p=1; p<=pmax; p++)
      if (L[p]) gen_sort_inplace(gel(L,p), (void*)&cmp_prime_over_p,
                                 &cmp_nodata, NULL);
  }
  else
  {
    for (i=1; i<l; i++)
    {
      pr = gel(FB,i); p = pr_get_smallp(pr);
      if (!L[p]) gel(L,p) = vecsmalltrunc_init(N+1);
      vecsmalltrunc_append(gel(L,p), i);
    }
  }
  return L;
}

/* recover FB, LV, iLP, KCZ from Vbase */
static GEN
recover_partFB(FB_t *F, GEN Vbase, long N)
{
  GEN FB, LV, iLP, L = get_pr_lists(Vbase, N, 0);
  long l = lg(L), p, ip, i;

  i = ip = 0;
  FB = cgetg(l, t_VECSMALL);
  iLP= cgetg(l, t_VECSMALL);
  LV = cgetg(l, t_VEC);
  for (p = 2; p < l; p++)
  {
    if (!L[p]) continue;
    FB[++i] = p;
    gel(LV,p) = vecpermute(Vbase, gel(L,p));
    iLP[p]= ip; ip += lg(gel(L,p))-1;
  }
  F->KCZ = i;
  F->KC = ip;
  F->FB = FB; setlg(FB, i+1);
  F->prodZ = zv_prod_Z(F->FB);
  F->LV = LV;
  F->iLP= iLP; return L;
}

/* add v^e to factorization */
static void
add_to_fact(long v, long e, FACT *fact)
{
  long i, l = fact[0].pr;
  for (i=1; i<=l && fact[i].pr < v; i++)/*empty*/;
  if (i <= l && fact[i].pr == v) fact[i].ex += e; else store(v, e, fact);
}
static void
inv_fact(FACT *fact)
{
  long i, l = fact[0].pr;
  for (i=1; i<=l; i++) fact[i].ex = -fact[i].ex;
}

/* L (small) list of primes above the same p including pr. Return pr index */
static int
pr_index(GEN L, GEN pr)
{
  long j, l = lg(L);
  GEN al = pr_get_gen(pr);
  for (j=1; j<l; j++)
    if (ZV_equal(al, pr_get_gen(gel(L,j)))) return j;
  pari_err_BUG("codeprime");
  return 0; /* LCOV_EXCL_LINE */
}

static long
Vbase_to_FB(FB_t *F, GEN pr)
{
  long p = pr_get_smallp(pr);
  return F->iLP[p] + pr_index(gel(F->LV,p), pr);
}

/* x, y 2 extended ideals whose first component is an integral HNF and second
 * a famat */
static GEN
idealHNF_mulred(GEN nf, GEN x, GEN y)
{
  GEN A = idealHNF_mul(nf, gel(x,1), gel(y,1));
  GEN F = famat_mul_shallow(gel(x,2), gel(y,2));
  return idealred(nf, mkvec2(A, F));
}
/* idealred(x * pr^n), n > 0 is small, x extended ideal. Reduction in order to
 * avoid prec pb: don't let id become too large as lgsub increases */
static GEN
idealmulpowprime2(GEN nf, GEN x, GEN pr, ulong n)
{
  GEN A = idealmulpowprime(nf, gel(x,1), pr, utoipos(n));
  return mkvec2(A, gel(x,2));
}
static GEN
init_famat(GEN x) { return mkvec2(x, trivial_fact()); }
/* optimized idealfactorback + reduction; z = init_famat() */
static GEN
genback(GEN z, GEN nf, GEN P, GEN E)
{
  long i, l = lg(E);
  GEN I = NULL;
  for (i = 1; i < l; i++)
    if (signe(gel(E,i)))
    {
      GEN J;
      gel(z,1) = gel(P,i);
      J = idealpowred(nf, z, gel(E,i));
      I = I? idealHNF_mulred(nf, I, J): J;
    }
  return I; /* != NULL since a generator */
}

/* return famat y (principal ideal) such that y / x is smooth [wrt Vbase] */
static GEN
SPLIT(FB_t *F, GEN nf, GEN x, GEN Vbase, FACT *fact)
{
  GEN vecG, ex, Ly, y, x0, Nx = ZM_det_triangular(x);
  long nbtest_lim, nbtest, i, j, k, ru, lgsub;
  pari_sp av;

  /* try without reduction if x is small */
  if (gexpo(gcoeff(x,1,1)) < 100 &&
      can_factor(F, nf, x, NULL, Nx, fact)) return NULL;

  av = avma;
  Ly = idealpseudominvec(x, nf_get_roundG(nf));
  for(k=1; k<lg(Ly); k++)
  {
    y = gel(Ly,k);
    if (factorgen(F, nf, x, Nx, y, fact)) return y;
  }
  set_avma(av);

  /* reduce in various directions */
  ru = lg(nf_get_roots(nf));
  vecG = cgetg(ru, t_VEC);
  for (j=1; j<ru; j++)
  {
    gel(vecG,j) = nf_get_Gtwist1(nf, j);
    av = avma;
    Ly = idealpseudominvec(x, gel(vecG,j));
    for(k=1; k<lg(Ly); k++)
    {
      y = gel(Ly,k);
      if (factorgen(F, nf, x, Nx, y, fact)) return y;
    }
    set_avma(av);
  }

  /* tough case, multiply by random products */
  lgsub = 3;
  ex = cgetg(lgsub, t_VECSMALL);
  x0 = init_famat(x);
  nbtest = 1; nbtest_lim = 4;
  for(;;)
  {
    GEN Ired, I, NI, id = x0;
    av = avma;
    if (DEBUGLEVEL>2) err_printf("# ideals tried = %ld\n",nbtest);
    for (i=1; i<lgsub; i++)
    {
      ex[i] = random_bits(RANDOM_BITS);
      if (ex[i]) id = idealmulpowprime2(nf, id, gel(Vbase,i), ex[i]);
    }
    if (id == x0) continue;
    /* I^(-1) * \prod Vbase[i]^ex[i] = (id[2]) / x */

    I = gel(id,1); NI = ZM_det_triangular(I);
    if (can_factor(F, nf, I, NULL, NI, fact))
    {
      inv_fact(fact); /* I^(-1) */
      for (i=1; i<lgsub; i++)
        if (ex[i]) add_to_fact(Vbase_to_FB(F,gel(Vbase,i)), ex[i], fact);
      return gel(id,2);
    }
    Ired = ru == 2? I: ZM_lll(I, 0.99, LLL_INPLACE);
    for (j=1; j<ru; j++)
    {
      pari_sp av2 = avma;
      Ly = idealpseudominvec(Ired, gel(vecG,j));
      for (k=1; k < lg(Ly); k++)
      {
        y = gel(Ly,k);
        if (factorgen(F, nf, I, NI, y, fact))
        {
          for (i=1; i<lgsub; i++)
            if (ex[i]) add_to_fact(Vbase_to_FB(F,gel(Vbase,i)), ex[i], fact);
          return famat_mul_shallow(gel(id,2), y);
        }
      }
      set_avma(av2);
    }
    set_avma(av);
    if (++nbtest > nbtest_lim)
    {
      nbtest = 0;
      if (++lgsub < minss(8, lg(Vbase)-1))
      {
        nbtest_lim <<= 1;
        ex = cgetg(lgsub, t_VECSMALL);
      }
      else nbtest_lim = LONG_MAX; /* don't increase further */
      if (DEBUGLEVEL>2) err_printf("SPLIT: increasing factor base [%ld]\n",lgsub);
    }
  }
}

INLINE GEN
bnf_get_W(GEN bnf) { return gel(bnf,1); }
INLINE GEN
bnf_get_B(GEN bnf) { return gel(bnf,2); }
INLINE GEN
bnf_get_C(GEN bnf) { return gel(bnf,4); }
INLINE GEN
bnf_get_vbase(GEN bnf) { return gel(bnf,5); }
INLINE GEN
bnf_get_Ur(GEN bnf) { return gmael(bnf,9,1); }
INLINE GEN
bnf_get_ga(GEN bnf) { return gmael(bnf,9,2); }
INLINE GEN
bnf_get_GD(GEN bnf) { return gmael(bnf,9,3); }

/* Return y (as an elt of K or a t_MAT representing an elt in Z[K])
 * such that x / (y) is smooth and store the exponents of  its factorization
 * on g_W and g_B in Wex / Bex; return NULL for y = 1 */
static GEN
split_ideal(GEN bnf, GEN x, GEN *pWex, GEN *pBex)
{
  GEN L, y, Vbase = bnf_get_vbase(bnf);
  GEN Wex, W  = bnf_get_W(bnf);
  GEN Bex, B  = bnf_get_B(bnf);
  long p, j, i, l, nW, nB;
  FACT *fact;
  FB_t F;

  L = recover_partFB(&F, Vbase, lg(x)-1);
  fact = (FACT*)stack_malloc((F.KC+1)*sizeof(FACT));
  y = SPLIT(&F, bnf_get_nf(bnf), x, Vbase, fact);
  nW = lg(W)-1; *pWex = Wex = zero_zv(nW);
  nB = lg(B)-1; *pBex = Bex = zero_zv(nB); l = lg(F.FB);
  p = j = 0; /* -Wall */
  for (i = 1; i <= fact[0].pr; i++)
  { /* decode index C = ip+j --> (p,j) */
    long a, b, t, C = fact[i].pr;
    for (t = 1; t < l; t++)
    {
      long q = F.FB[t], k = C - F.iLP[q];
      if (k <= 0) break;
      p = q;
      j = k;
    }
    a = gel(L, p)[j];
    b = a - nW;
    if (b <= 0) Wex[a] = y? -fact[i].ex: fact[i].ex;
    else        Bex[b] = y? -fact[i].ex: fact[i].ex;
  }
  return y;
}

GEN
init_red_mod_units(GEN bnf, long prec)
{
  GEN s = gen_0, p1,s1,mat, logfu = bnf_get_logfu(bnf);
  long i,j, RU = lg(logfu);

  if (RU == 1) return NULL;
  mat = cgetg(RU,t_MAT);
  for (j=1; j<RU; j++)
  {
    p1 = cgetg(RU+1,t_COL); gel(mat,j) = p1;
    s1 = gen_0;
    for (i=1; i<RU; i++)
    {
      gel(p1,i) = real_i(gcoeff(logfu,i,j));
      s1 = mpadd(s1, mpsqr(gel(p1,i)));
    }
    gel(p1,RU) = gen_0; if (mpcmp(s1,s) > 0) s = s1;
  }
  s = gsqrt(gmul2n(s,RU),prec);
  if (expo(s) < 27) s = utoipos(1UL << 27);
  return mkvec2(mat, s);
}

/* z computed above. Return unit exponents that would reduce col (arch) */
GEN
red_mod_units(GEN col, GEN z)
{
  long i,RU;
  GEN x,mat,N2;

  if (!z) return NULL;
  mat= gel(z,1);
  N2 = gel(z,2);
  RU = lg(mat); x = cgetg(RU+1,t_COL);
  for (i=1; i<RU; i++) gel(x,i) = real_i(gel(col,i));
  gel(x,RU) = N2;
  x = lll(shallowconcat(mat,x));
  if (typ(x) != t_MAT || lg(x) <= RU) return NULL;
  x = gel(x,RU);
  if (signe(gel(x,RU)) < 0) x = gneg_i(x);
  if (!gequal1(gel(x,RU))) pari_err_BUG("red_mod_units");
  setlg(x,RU); return x;
}

static GEN
add(GEN a, GEN t) { return a = a? RgC_add(a,t): t; }

/* [x] archimedian components, A column vector. return [x] A */
static GEN
act_arch(GEN A, GEN x)
{
  GEN a;
  long i,l = lg(A), tA = typ(A);
  if (tA == t_MAT)
  { /* assume lg(x) >= l */
    a = cgetg(l, t_MAT);
    for (i=1; i<l; i++) gel(a,i) = act_arch(gel(A,i), x);
    return a;
  }
  if (l==1) return cgetg(1, t_COL);
  a = NULL;
  if (tA == t_VECSMALL)
  {
    for (i=1; i<l; i++)
    {
      long c = A[i];
      if (c) a = add(a, gmulsg(c, gel(x,i)));
    }
  }
  else
  { /* A a t_COL of t_INT. Assume lg(A)==lg(x) */
    for (i=1; i<l; i++)
    {
      GEN c = gel(A,i);
      if (signe(c)) a = add(a, gmul(c, gel(x,i)));
    }
  }
  return a? a: zerocol(lgcols(x)-1);
}
/* act_arch(matdiagonal(v), x) */
static GEN
diagact_arch(GEN v, GEN x)
{
  long i, l = lg(v);
  GEN a = cgetg(l, t_MAT);
  for (i = 1; i < l; i++) gel(a,i) = gmul(gel(x,i), gel(v,i));
  return a;
}

static long
prec_arch(GEN bnf)
{
  GEN a = bnf_get_C(bnf);
  long i, l = lg(a), prec;

  for (i=1; i<l; i++)
    if ( (prec = gprecision(gel(a,i))) ) return prec;
  return DEFAULTPREC;
}

static long
needed_bitprec(GEN x)
{
  long i, e = 0, l = lg(x);
  for (i = 1; i < l; i++)
  {
    GEN c = gel(x,i);
    long f = gexpo(c) - prec2nbits(gprecision(c));
    if (f > e) e = f;
  }
  return e;
}

/* col = archimedian components of x, Nx its norm, dx a multiple of its
 * denominator. Return x or NULL (fail) */
GEN
isprincipalarch(GEN bnf, GEN col, GEN kNx, GEN e, GEN dx, long *pe)
{
  GEN nf, x, y, logfu, s, M;
  long N, prec = gprecision(col);
  bnf = checkbnf(bnf); nf = bnf_get_nf(bnf); M = nf_get_M(nf);
  if (!prec) prec = prec_arch(bnf);
  *pe = 128;
  logfu = bnf_get_logfu(bnf);
  N = nf_get_degree(nf);
  if (!(col = cleanarch(col,N,NULL,prec))) return NULL;
  if (lg(col) > 2)
  { /* reduce mod units */
    GEN u, z = init_red_mod_units(bnf,prec);
    if (!(u = red_mod_units(col,z))) return NULL;
    col = RgC_add(col, RgM_RgC_mul(logfu, u));
    if (!(col = cleanarch(col,N,NULL,prec))) return NULL;
  }
  s = divru(mulir(e, glog(kNx,prec)), N);
  col = fixarch(col, s, nf_get_r1(nf));
  if (RgC_expbitprec(col) >= 0) return NULL;
  col = gexp(col, prec);
  /* d.alpha such that x = alpha \prod gj^ej */
  x = RgM_solve_realimag(M,col); if (!x) return NULL;
  x = RgC_Rg_mul(x, dx);
  y = grndtoi(x, pe);
  if (*pe > -5) { *pe = needed_bitprec(x); return NULL; }
  return RgC_Rg_div(y, dx);
}

/* y = C \prod g[i]^e[i] ? */
static int
fact_ok(GEN nf, GEN y, GEN C, GEN g, GEN e)
{
  pari_sp av = avma;
  long i, c = lg(e);
  GEN z = C? C: gen_1;
  for (i=1; i<c; i++)
    if (signe(gel(e,i))) z = idealmul(nf, z, idealpow(nf, gel(g,i), gel(e,i)));
  if (typ(z) != t_MAT) z = idealhnf_shallow(nf,z);
  if (typ(y) != t_MAT) y = idealhnf_shallow(nf,y);
  return gc_bool(av, ZM_equal(y,z));
}
static GEN
ZV_divrem(GEN A, GEN B, GEN *pR)
{
  long i, l = lg(A);
  GEN Q = cgetg(l, t_COL), R = cgetg(l, t_COL);
  for (i = 1; i < l; i++) gel(Q,i) = truedvmdii(gel(A,i), gel(B,i), &gel(R,i));
  *pR = R; return Q;
}

static GEN
Ur_ZC_mul(GEN bnf, GEN v)
{
  GEN w, U = bnf_get_Ur(bnf);
  long i, l = lg(bnf_get_cyc(bnf)); /* may be < lgcols(U) */

  w = cgetg(l, t_COL);
  for (i = 1; i < l; i++) gel(w,i) = ZMrow_ZC_mul(U, v, i);
  return w;
}

static GEN
ZV_mul(GEN x, GEN y)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_COL);
  for (i = 1; i < l; i++) gel(z,i) = mulii(gel(x,i), gel(y,i));
  return z;
}
static int
dump_gen(GEN SUnits, GEN x, long flag)
{
  GEN d;
  long e;
  if (!(flag & nf_GENMAT) || !SUnits) return 0;
  e = gexpo(gel(SUnits,2)); if (e > 64) return 0; /* U large */
  x = Q_remove_denom(x, &d);
  return (d && expi(d) > 32) || gexpo(x) > 32;
}

/* assume x in HNF; cf class_group_gen for notations. Return NULL iff
 * flag & nf_FORCE and computation of principal ideal generator fails */
static GEN
isprincipalall(GEN bnf, GEN x, long *pprec, long flag)
{
  GEN xar, Wex, Bex, gen, xc, col, A, Q, R, UA, SUnits;
  GEN C = bnf_get_C(bnf), nf = bnf_get_nf(bnf), cyc = bnf_get_cyc(bnf);
  long nB, nW, e;

  if (lg(cyc) == 1 && !(flag & (nf_GEN|nf_GENMAT|nf_GEN_IF_PRINCIPAL)))
    return cgetg(1,t_COL);
  if (lg(x) == 2)
  { /* nf = Q */
    col = gel(x,1);
    if (flag & nf_GENMAT) col = to_famat_shallow(col, gen_1);
    return (flag & nf_GEN_IF_PRINCIPAL)? col: mkvec2(cgetg(1,t_COL), col);
  }

  x = Q_primitive_part(x, &xc);
  if (equali1(gcoeff(x,1,1))) /* trivial ideal */
  {
    R = zerocol(lg(cyc)-1);
    if (!(flag & (nf_GEN|nf_GENMAT|nf_GEN_IF_PRINCIPAL))) return R;
    if (flag & nf_GEN_IF_PRINCIPAL)
      return scalarcol_shallow(xc? xc: gen_1, nf_get_degree(nf));
    if (flag & nf_GENMAT)
      col = xc? to_famat_shallow(xc, gen_1): trivial_fact();
    else
      col = scalarcol_shallow(xc? xc: gen_1, nf_get_degree(nf));
    return mkvec2(R, col);
  }
  xar = split_ideal(bnf, x, &Wex, &Bex);
  /* x = g_W Wex + g_B Bex + [xar] = g_W (Wex - B*Bex) + [xar] + [C_B]Bex */
  A = zc_to_ZC(Wex); nB = lg(Bex)-1;
  if (nB) A = ZC_sub(A, ZM_zc_mul(bnf_get_B(bnf), Bex));
  UA = Ur_ZC_mul(bnf, A);
  Q = ZV_divrem(UA, cyc, &R);
  /* g_W (Wex - B*Bex) = G Ur A - [ga]A = G R + [GD]Q - [ga]A
   * Finally: x = G R + [xar] + [C_B]Bex + [GD]Q - [ga]A */
  if (!(flag & (nf_GEN|nf_GENMAT|nf_GEN_IF_PRINCIPAL))) return R;
  if ((flag & nf_GEN_IF_PRINCIPAL) && !ZV_equal0(R)) return gen_0;

  nW = lg(Wex)-1;
  gen = bnf_get_gen(bnf);
  col = NULL;
  SUnits = bnf_get_sunits(bnf);
  if (lg(R) == 1
      || abscmpiu(gel(R,vecindexmax(R)), 4 * bit_accuracy(*pprec)) < 0)
  { /* q = N (x / prod gj^ej) = N(alpha), denom(alpha) | d */
    GEN d, q = gdiv(ZM_det_triangular(x), get_norm_fact(gen, R, &d));
    col = xar? nf_cxlog(nf, xar, *pprec): NULL;
    if (nB) col = add(col, act_arch(Bex, nW? vecslice(C,nW+1,lg(C)-1): C));
    if (nW) col = add(col, RgC_sub(act_arch(Q, bnf_get_GD(bnf)),
                                   act_arch(A, bnf_get_ga(bnf))));
    col = isprincipalarch(bnf, col, q, gen_1, d, &e);
    if (col && (dump_gen(SUnits, col, flag)
                || !fact_ok(nf,x, col,gen,R))) col = NULL;
  }
  if (!col && (flag & nf_GENMAT))
  {
    if (SUnits)
    {
      GEN X = gel(SUnits,1), U = gel(SUnits,2), C = gel(SUnits,3);
      GEN v = gel(bnf,9), Ge = gel(v,4), M1 = gel(v,5), M2 = gel(v,6);
      GEN z = NULL, F = NULL;
      if (nB)
      {
        GEN C2 = nW? vecslice(C, nW+1, lg(C)-1): C;
        z = ZM_zc_mul(C2, Bex);
      }
      if (nW)
      { /* [GD]Q - [ga]A = ([X]M1 - [Ge]D) Q - ([X]M2 - [Ge]Ur) A */
        GEN C1 = vecslice(C, 1, nW);
        GEN v = ZC_sub(ZM_ZC_mul(M1,Q), ZM_ZC_mul(M2,A));
        z = add(z, ZM_ZC_mul(C1, v));
        F = famat_reduce(famatV_factorback(Ge, ZC_sub(UA, ZV_mul(cyc,Q))));
        if (lgcols(F) == 1) F = NULL;
      }
      /* reduce modulo units and Q^* */
      if (lg(U) != 1) z = ZC_sub(z, ZM_ZC_mul(U, RgM_Babai(U,z)));
      col = mkmat2(X, z);
      if (F) col = famat_mul_shallow(col, F);
      col = famat_remove_trivial(col);
      if (xar) col = famat_mul_shallow(col, xar);
    }
    else if (!ZV_equal0(R))
    { /* in case isprincipalfact calls bnfinit() due to prec trouble...*/
      GEN y = isprincipalfact(bnf, x, gen, ZC_neg(R), flag);
      if (typ(y) != t_VEC) return y;
      col = gel(y,2);
    }
  }
  if (col)
  { /* add back missing content */
    if (xc) col = (typ(col)==t_MAT)? famat_mul_shallow(col,xc)
                                   : RgC_Rg_mul(col,xc);
    if (typ(col) != t_MAT && lg(col) != 1 && (flag & nf_GENMAT))
      col = to_famat_shallow(col, gen_1);
  }
  else
  {
    if (e < 0) e = 0;
    *pprec += nbits2extraprec(e + 128);
    if (flag & nf_FORCE)
    {
      if (DEBUGLEVEL)
        pari_warn(warner,"precision too low for generators, e = %ld",e);
      return NULL;
    }
    pari_warn(warner,"precision too low for generators, not given");
    col = cgetg(1, t_COL);
  }
  return (flag & nf_GEN_IF_PRINCIPAL)? col: mkvec2(R, col);
}

static GEN
triv_gen(GEN bnf, GEN x, long flag)
{
  pari_sp av = avma;
  GEN nf = bnf_get_nf(bnf);
  long c;
  if (flag & nf_GEN_IF_PRINCIPAL)
  {
    if (!(flag & nf_GENMAT)) return algtobasis(nf,x);
    x = nf_to_scalar_or_basis(nf,x);
    if (typ(x) == t_INT && is_pm1(x)) return trivial_fact();
    return gerepilecopy(av, to_famat_shallow(x, gen_1));
  }
  c = lg(bnf_get_cyc(bnf)) - 1;
  if (flag & nf_GENMAT)
    retmkvec2(zerocol(c), to_famat_shallow(algtobasis(nf,x), gen_1));
  if (flag & nf_GEN)
    retmkvec2(zerocol(c), algtobasis(nf,x));
  return zerocol(c);
}

GEN
bnfisprincipal0(GEN bnf,GEN x,long flag)
{
  pari_sp av = avma;
  GEN c, nf;
  long pr;

  bnf = checkbnf(bnf);
  nf = bnf_get_nf(bnf);
  switch( idealtyp(&x, NULL) )
  {
    case id_PRINCIPAL:
      if (gequal0(x)) pari_err_DOMAIN("bnfisprincipal","ideal","=",gen_0,x);
      return triv_gen(bnf, x, flag);
    case id_PRIME:
      if (pr_is_inert(x)) return triv_gen(bnf, pr_get_p(x), flag);
      x = pr_hnf(nf, x);
      break;
    case id_MAT:
      if (lg(x)==1) pari_err_DOMAIN("bnfisprincipal","ideal","=",gen_0,x);
      if (nf_get_degree(nf) != lg(x)-1)
        pari_err_TYPE("idealtyp [dimension != degree]", x);
  }
  pr = prec_arch(bnf); /* precision of unit matrix */
  c = getrand();
  for (;;)
  {
    pari_sp av1 = avma;
    GEN y = isprincipalall(bnf,x,&pr,flag);
    if (y) return gerepilecopy(av, y);

    if (DEBUGLEVEL) pari_warn(warnprec,"isprincipal",pr);
    set_avma(av1); bnf = bnfnewprec_shallow(bnf,pr); setrand(c);
  }
}
GEN
isprincipal(GEN bnf,GEN x) { return bnfisprincipal0(bnf,x,0); }

/* FIXME: OBSOLETE */
GEN
isprincipalgen(GEN bnf,GEN x)
{ return bnfisprincipal0(bnf,x,nf_GEN); }
GEN
isprincipalforce(GEN bnf,GEN x)
{ return bnfisprincipal0(bnf,x,nf_FORCE); }
GEN
isprincipalgenforce(GEN bnf,GEN x)
{ return bnfisprincipal0(bnf,x,nf_GEN | nf_FORCE); }

/* lg(u) > 1 */
static int
RgV_is1(GEN u) { return isint1(gel(u,1)) && RgV_isscalar(u); }
static GEN
add_principal_part(GEN nf, GEN u, GEN v, long flag)
{
  if (flag & nf_GENMAT)
    return (typ(u) == t_COL && RgV_is1(u))? v: famat_mul_shallow(v,u);
  else
    return nfmul(nf, v, u);
}

#if 0
/* compute C prod P[i]^e[i],  e[i] >=0 for all i. C may be NULL (omitted)
 * e destroyed ! */
static GEN
expand(GEN nf, GEN C, GEN P, GEN e)
{
  long i, l = lg(e), done = 1;
  GEN id = C;
  for (i=1; i<l; i++)
  {
    GEN ei = gel(e,i);
    if (signe(ei))
    {
      if (mod2(ei)) id = id? idealmul(nf, id, gel(P,i)): gel(P,i);
      ei = shifti(ei,-1);
      if (signe(ei)) done = 0;
      gel(e,i) = ei;
    }
  }
  if (id != C) id = idealred(nf, id);
  if (done) return id;
  return idealmulred(nf, id, idealsqr(nf, expand(nf,id,P,e)));
}
/* C is an extended ideal, possibly with C[1] = NULL */
static GEN
expandext(GEN nf, GEN C, GEN P, GEN e)
{
  long i, l = lg(e), done = 1;
  GEN A = gel(C,1);
  for (i=1; i<l; i++)
  {
    GEN ei = gel(e,i);
    if (signe(ei))
    {
      if (mod2(ei)) A = A? idealmul(nf, A, gel(P,i)): gel(P,i);
      ei = shifti(ei,-1);
      if (signe(ei)) done = 0;
      gel(e,i) = ei;
    }
  }
  if (A == gel(C,1))
    A = C;
  else
    A = idealred(nf, mkvec2(A, gel(C,2)));
  if (done) return A;
  return idealmulred(nf, A, idealsqr(nf, expand(nf,A,P,e)));
}
#endif

static GEN
expand(GEN nf, GEN C, GEN P, GEN e)
{
  long i, l = lg(e);
  GEN B, A = C;
  for (i=1; i<l; i++) /* compute prod P[i]^e[i] */
    if (signe(gel(e,i)))
    {
      B = idealpowred(nf, gel(P,i), gel(e,i));
      A = A? idealmulred(nf,A,B): B;
    }
  return A;
}
static GEN
expandext(GEN nf, GEN C, GEN P, GEN e)
{
  long i, l = lg(e);
  GEN B, A = gel(C,1), C1 = A;
  for (i=1; i<l; i++) /* compute prod P[i]^e[i] */
    if (signe(gel(e,i)))
    {
      gel(C,1) = gel(P,i);
      B = idealpowred(nf, C, gel(e,i));
      A = A? idealmulred(nf,A,B): B;
    }
  return A == C1? C: A;
}

/* isprincipal for C * \prod P[i]^e[i] (C omitted if NULL) */
GEN
isprincipalfact(GEN bnf, GEN C, GEN P, GEN e, long flag)
{
  const long gen = flag & (nf_GEN|nf_GENMAT|nf_GEN_IF_PRINCIPAL);
  long prec;
  pari_sp av = avma;
  GEN C0, Cext, c, id, nf = bnf_get_nf(bnf);

  if (gen)
  {
    Cext = (flag & nf_GENMAT)? trivial_fact()
                             : mkpolmod(gen_1,nf_get_pol(nf));
    C0 = mkvec2(C, Cext);
    id = expandext(nf, C0, P, e);
  } else {
    Cext = NULL;
    C0 = C;
    id = expand(nf, C, P, e);
  }
  if (id == C0) /* e = 0 */
  {
    if (!C) return bnfisprincipal0(bnf, gen_1, flag);
    switch(typ(C))
    {
      case t_INT: case t_FRAC: case t_POL: case t_POLMOD: case t_COL:
        return triv_gen(bnf, C, flag);
    }
    C = idealhnf_shallow(nf,C);
  }
  else
  {
    if (gen) { C = gel(id,1); Cext = gel(id,2); } else C = id;
  }
  prec = prec_arch(bnf);
  c = getrand();
  for (;;)
  {
    pari_sp av1 = avma;
    GEN y = isprincipalall(bnf, C, &prec, flag);
    if (y)
    {
      if (flag & nf_GEN_IF_PRINCIPAL)
      {
        if (typ(y) == t_INT) return gc_NULL(av);
        y = add_principal_part(nf, y, Cext, flag);
      }
      else
      {
        GEN u = gel(y,2);
        if (!gen || typ(y) != t_VEC) return gerepileupto(av,y);
        if (lg(u) != 1) gel(y,2) = add_principal_part(nf, u, Cext, flag);
      }
      return gerepilecopy(av, y);
    }
    if (DEBUGLEVEL) pari_warn(warnprec,"isprincipal",prec);
    set_avma(av1); bnf = bnfnewprec_shallow(bnf,prec); setrand(c);
  }
}
GEN
isprincipalfact_or_fail(GEN bnf, GEN C, GEN P, GEN e)
{
  const long flag = nf_GENMAT|nf_FORCE;
  long prec;
  pari_sp av = avma;
  GEN u, y, id, C0, Cext, nf = bnf_get_nf(bnf);

  Cext = trivial_fact();
  C0 = mkvec2(C, Cext);
  id = expandext(nf, C0, P, e);
  if (id == C0) /* e = 0 */
    C = idealhnf_shallow(nf,C);
  else {
    C = gel(id,1); Cext = gel(id,2);
  }
  prec = prec_arch(bnf);
  y = isprincipalall(bnf, C, &prec, flag);
  if (!y) return gc_utoipos(av, prec);
  u = gel(y,2);
  if (lg(u) != 1) gel(y,2) = add_principal_part(nf, u, Cext, flag);
  return gerepilecopy(av, y);
}

GEN
nfsign_from_logarch(GEN LA, GEN invpi, GEN archp)
{
  long l = lg(archp), i;
  GEN y = cgetg(l, t_VECSMALL);
  pari_sp av = avma;

  for (i=1; i<l; i++)
  {
    GEN c = ground( gmul(imag_i(gel(LA,archp[i])), invpi) );
    y[i] = mpodd(c)? 1: 0;
  }
  set_avma(av); return y;
}

GEN
nfsign_tu(GEN bnf, GEN archp)
{
  long n;
  if (bnf_get_tuN(bnf) != 2) return cgetg(1, t_VECSMALL);
  n = archp? lg(archp) - 1: nf_get_r1(bnf_get_nf(bnf));
  return const_vecsmall(n, 1);
}
GEN
nfsign_fu(GEN bnf, GEN archp)
{
  GEN invpi, y, A = bnf_get_logfu(bnf), nf = bnf_get_nf(bnf);
  long j = 1, RU = lg(A);

  if (!archp) archp = identity_perm( nf_get_r1(nf) );
  invpi = invr( mppi(nf_get_prec(nf)) );
  y = cgetg(RU,t_MAT);
  for (j = 1; j < RU; j++)
    gel(y,j) = nfsign_from_logarch(gel(A,j), invpi, archp);
  return y;
}
GEN
nfsign_units(GEN bnf, GEN archp, int add_zu)
{
  GEN sfu = nfsign_fu(bnf, archp);
  return add_zu? vec_prepend(sfu, nfsign_tu(bnf, archp)): sfu;
}

/* obsolete */
GEN
signunits(GEN bnf)
{
  pari_sp av;
  GEN S, y, nf;
  long i, j, r1, r2;

  bnf = checkbnf(bnf); nf = bnf_get_nf(bnf);
  nf_get_sign(nf, &r1,&r2);
  S = zeromatcopy(r1, r1+r2-1); av = avma;
  y = nfsign_fu(bnf, NULL);
  for (j = 1; j < lg(y); j++)
  {
    GEN Sj = gel(S,j), yj = gel(y,j);
    for (i = 1; i <= r1; i++) gel(Sj,i) = yj[i]? gen_m1: gen_1;
  }
  set_avma(av); return S;
}

static GEN
get_log_embed(REL_t *rel, GEN M, long RU, long R1, long prec)
{
  GEN arch, C, z = rel->m;
  long i;
  arch = typ(z) == t_COL? RgM_RgC_mul(M, z): const_col(nbrows(M), z);
  C = cgetg(RU+1, t_COL); arch = glog(arch, prec);
  for (i=1; i<=R1; i++) gel(C,i) = gel(arch,i);
  for (   ; i<=RU; i++) gel(C,i) = gmul2n(gel(arch,i), 1);
  return C;
}
static GEN
rel_embed(REL_t *rel, FB_t *F, GEN embs, long ind, GEN M, long RU, long R1,
          long prec)
{
  GEN C, D, perm;
  long i, n;
  if (!rel->relaut) return get_log_embed(rel, M, RU, R1, prec);
  /* image of another relation by automorphism */
  C = gel(embs, ind - rel->relorig);
  perm = gel(F->embperm, rel->relaut);
  D = cgetg_copy(C, &n);
  for (i = 1; i < n; i++)
  {
    long v = perm[i];
    gel(D,i) = (v > 0)? gel(C,v): conj_i(gel(C,-v));
  }
  return D;
}
static GEN
get_embs(FB_t *F, RELCACHE_t *cache, GEN nf, GEN embs, long PREC)
{
  long ru, j, k, l = cache->last - cache->chk + 1, r1 = nf_get_r1(nf);
  GEN M = nf_get_M(nf), nembs = cgetg(cache->last - cache->base+1, t_MAT);
  REL_t *rel;

  for (k = 1; k <= cache->chk - cache->base; k++) gel(nembs,k) = gel(embs,k);
  embs = nembs; ru = nbrows(M);
  for (j=1,rel = cache->chk + 1; j < l; rel++,j++,k++)
    gel(embs,k) = rel_embed(rel, F, embs, k, M, ru, r1, PREC);
  return embs;
}
static void
set_rel_alpha(REL_t *rel, GEN auts, GEN vA, long ind)
{
  GEN u;
  if (!rel->relaut)
    u = rel->m;
  else
    u = ZM_ZC_mul(gel(auts, rel->relaut), gel(vA, ind - rel->relorig));
  gel(vA, ind) = u;
}
static GEN
set_fact(FB_t *F, FACT *fact, GEN e, long *pnz)
{
  long n = fact[0].pr;
  GEN c = zero_Flv(F->KC);
  if (!n) /* trivial factorization */
    *pnz = F->KC+1;
  else
  {
    long i, nz = minss(fact[1].pr, fact[n].pr);
    for (i = 1; i <= n; i++) c[fact[i].pr] = fact[i].ex;
    if (e)
    {
      long l = lg(e);
      for (i = 1; i < l; i++)
        if (e[i]) { long v = F->subFB[i]; c[v] += e[i]; if (v < nz) nz = v; }
    }
    *pnz = nz;
  }
  return c;
}

/* Is cols already in the cache ? bs = index of first non zero coeff in cols
 * General check for colinearity useless since exceedingly rare */
static int
already_known(RELCACHE_t *cache, long bs, GEN cols)
{
  REL_t *r;
  long l = lg(cols);
  for (r = cache->last; r > cache->base; r--)
    if (bs == r->nz)
    {
      GEN coll = r->R;
      long b = bs;
      while (b < l && cols[b] == coll[b]) b++;
      if (b == l) return 1;
    }
  return 0;
}

/* Add relation R to cache, nz = index of first non zero coeff in R.
 * If relation is a linear combination of the previous ones, return 0.
 * Otherwise, update basis and return > 0. Compute mod p (much faster)
 * so some kernel vector might not be genuine. */
static int
add_rel_i(RELCACHE_t *cache, GEN R, long nz, GEN m, long orig, long aut, REL_t **relp, long in_rnd_rel)
{
  long i, k, n = lg(R)-1;

  if (nz == n+1) { k = 0; goto ADD_REL; }
  if (already_known(cache, nz, R)) return -1;
  if (cache->last >= cache->base + cache->len) return 0;
  if (DEBUGLEVEL>6)
  {
    err_printf("adding vector = %Ps\n",R);
    err_printf("generators =\n%Ps\n", cache->basis);
  }
  if (cache->missing)
  {
    GEN a = leafcopy(R), basis = cache->basis;
    k = lg(a);
    do --k; while (!a[k]);
    while (k)
    {
      GEN c = gel(basis, k);
      if (c[k])
      {
        long ak = a[k];
        for (i=1; i < k; i++) if (c[i]) a[i] = (a[i] + ak*(mod_p-c[i])) % mod_p;
        a[k] = 0;
        do --k; while (!a[k]); /* k cannot go below 0: codeword is a sentinel */
      }
      else
      {
        ulong invak = Fl_inv(uel(a,k), mod_p);
        /* Cleanup a */
        for (i = k; i-- > 1; )
        {
          long j, ai = a[i];
          c = gel(basis, i);
          if (!ai || !c[i]) continue;
          ai = mod_p-ai;
          for (j = 1; j < i; j++) if (c[j]) a[j] = (a[j] + ai*c[j]) % mod_p;
          a[i] = 0;
        }
        /* Insert a/a[k] as k-th column */
        c = gel(basis, k);
        for (i = 1; i<k; i++) if (a[i]) c[i] = (a[i] * invak) % mod_p;
        c[k] = 1; a = c;
        /* Cleanup above k */
        for (i = k+1; i<n; i++)
        {
          long j, ck;
          c = gel(basis, i);
          ck = c[k];
          if (!ck) continue;
          ck = mod_p-ck;
          for (j = 1; j < k; j++) if (a[j]) c[j] = (c[j] + ck*a[j]) % mod_p;
          c[k] = 0;
        }
        cache->missing--;
        break;
      }
    }
  }
  else
    k = (cache->last - cache->base) + 1;
  if (k || cache->relsup > 0 || (m && in_rnd_rel))
  {
    REL_t *rel;

ADD_REL:
    rel = ++cache->last;
    if (!k && cache->relsup && nz < n+1)
    {
      cache->relsup--;
      k = (rel - cache->base) + cache->missing;
    }
    rel->R  = gclone(R);
    rel->m  =  m ? gclone(m) : NULL;
    rel->nz = nz;
    if (aut)
    {
      rel->relorig = (rel - cache->base) - orig;
      rel->relaut = aut;
    }
    else
      rel->relaut = 0;
    if (relp) *relp = rel;
    if (DEBUGLEVEL) dbg_newrel(cache);
  }
  return k;
}

static int
add_rel(RELCACHE_t *cache, FB_t *F, GEN R, long nz, GEN m, long in_rnd_rel)
{
  REL_t *rel;
  long k, l, reln;
  const long lauts = lg(F->idealperm), KC = F->KC;

  k = add_rel_i(cache, R, nz, m, 0, 0, &rel, in_rnd_rel);
  if (k > 0 && typ(m) != t_INT)
  {
    GEN Rl = cgetg(KC+1, t_VECSMALL);
    reln = rel - cache->base;
    for (l = 1; l < lauts; l++)
    {
      GEN perml = gel(F->idealperm, l);
      long i, nzl = perml[nz];

      for (i = 1; i <= KC; i++) Rl[i] = 0;
      for (i = nz; i <= KC; i++)
        if (R[i])
        {
          long v = perml[i];

          if (v < nzl) nzl = v;
          Rl[v] = R[i];
        }
      (void)add_rel_i(cache, Rl, nzl, NULL, reln, l, NULL, in_rnd_rel);
    }
  }
  return k;
}

INLINE void
step(GEN x, double *y, GEN inc, long k)
{
  if (!y[k])
    x[k]++; /* leading coeff > 0 */
  else
  {
    long i = inc[k];
    x[k] += i;
    inc[k] = (i > 0)? -1-i: 1-i;
  }
}

static double
Fincke_Pohst_bound(double T, GEN r)
{
  pari_sp av = avma;
  GEN zT = dbltor(T * T), p = gmael(r,1,1), B = real_1(DEFAULTPREC);
  long i, n = lg(r)-1;
  for (i = 2; i <= n; i++)
  {
    p = gmul(p, gmael(r,i,i));
    B = sqrtnr(gmul(zT,p), i);
    if (i == n || cmprr(B, gmael(r,i+1,i+1)) < 0) break;
  }
  return gc_double(av, rtodbl(B));
}

INLINE long
Fincke_Pohst_ideal(RELCACHE_t *cache, FB_t *F, GEN nf, GEN M, GEN I,
    GEN NI, FACT *fact, long Nrelid, FP_t *fp, RNDREL_t *rr, long prec,
    long *Nsmall, long *Nfact)
{
  pari_sp av;
  const long N = nf_get_degree(nf), R1 = nf_get_r1(nf);
  GEN G = nf_get_G(nf), G0 = nf_get_roundG(nf), r, u, gx, inc, ideal;
  double BOUND, B1, B2;
  long j, k, skipfirst, relid=0, try_factor=0;

  inc = const_vecsmall(N, 1);
  u = ZM_lll(ZM_mul(G0, I), 0.99, LLL_IM);
  ideal = ZM_mul(I,u); /* approximate T2-LLL reduction */
  r = gaussred_from_QR(RgM_mul(G, ideal), prec); /* Cholesky for T2 | ideal */
  if (!r) pari_err_BUG("small_norm (precision too low)");

  for (k=1; k<=N; k++)
  {
    if (!gisdouble(gcoeff(r,k,k),&(fp->v[k]))) return 0;
    for (j=1; j<k; j++) if (!gisdouble(gcoeff(r,j,k),&(fp->q[j][k]))) return 0;
    if (DEBUGLEVEL>3) err_printf("v[%ld]=%.4g ",k,fp->v[k]);
  }
  B1 = fp->v[1]; /* T2(ideal[1]) */
  B2 = fp->v[2] + B1 * fp->q[1][2] * fp->q[1][2]; /* T2(ideal[2]) */
  skipfirst = ZV_isscalar(gel(ideal,1));
  BOUND = maxdd(2*B2, Fincke_Pohst_bound(4 * maxtry_FACT / F->ballvol, r));
  if (DEBUGLEVEL>1)
  {
    if (DEBUGLEVEL>3) err_printf("\n");
    err_printf("BOUND = %.4g\n",BOUND);
  }

  k = N; fp->y[N] = fp->z[N] = 0; fp->x[N] = 0;
  for (av = avma;; set_avma(av), step(fp->x,fp->y,inc,k))
  {
    GEN R;
    long nz;
    do
    { /* look for primitive element of small norm, cf minim00 */
      int fl = 0;
      double p;
      if (k > 1)
      {
        long l = k-1;
        fp->z[l] = 0;
        for (j=k; j<=N; j++) fp->z[l] += fp->q[l][j]*fp->x[j];
        p = (double)fp->x[k] + fp->z[k];
        fp->y[l] = fp->y[k] + p*p*fp->v[k];
        if (l <= skipfirst && !fp->y[1]) fl = 1;
        fp->x[l] = (long)floor(-fp->z[l] + 0.5);
        k = l;
      }
      for(;; step(fp->x,fp->y,inc,k))
      {
        if (!fl)
        {
          p = (double)fp->x[k] + fp->z[k];
          if (fp->y[k] + p*p*fp->v[k] <= BOUND) break;

          step(fp->x,fp->y,inc,k);

          p = (double)fp->x[k] + fp->z[k];
          if (fp->y[k] + p*p*fp->v[k] <= BOUND) break;
        }
        fl = 0; inc[k] = 1;
        if (++k > N) goto END_Fincke_Pohst_ideal;
      }
    } while (k > 1);

    /* element complete */
    if (zv_content(fp->x) !=1) continue; /* not primitive */
    gx = ZM_zc_mul(ideal,fp->x);
    if (ZV_isscalar(gx)) continue;
    if (++try_factor > maxtry_FACT) break;

    if (!Nrelid)
    {
      if (!factorgen(F,nf,I,NI,gx,fact)) continue;
      return 1;
    }
    else if (rr)
    {
      if (!factorgen(F,nf,I,NI,gx,fact)) continue;
      add_to_fact(rr->jid, 1, fact);
    }
    else
    {
      GEN Nx, xembed = RgM_RgC_mul(M, gx);
      long e;
      if (Nsmall) (*Nsmall)++;
      Nx = grndtoi(embed_norm(xembed, R1), &e);
      if (e >= 0) {
        if (DEBUGLEVEL > 1) err_printf("+");
        continue;
      }
      if (!can_factor(F, nf, NULL, gx, Nx, fact)) continue;
    }

    /* smooth element */
    R = set_fact(F, fact, rr ? rr->ex : NULL, &nz);
    /* make sure we get maximal rank first, then allow all relations */
    if (add_rel(cache, F, R, nz, gx, rr ? 1 : 0) <= 0)
    { /* probably Q-dependent from previous ones: forget it */
      if (DEBUGLEVEL>1) err_printf("*");
      if (DEBUGLEVEL && Nfact && rr) (*Nfact)++;
      continue;
    }
    if (DEBUGLEVEL && Nfact) (*Nfact)++;
    if (cache->last >= cache->end) return 1; /* we have enough */
    if (++relid == Nrelid) break;
  }
  END_Fincke_Pohst_ideal:
  return 0;
}

static void
small_norm(RELCACHE_t *cache, FB_t *F, GEN nf, long Nrelid, GEN M,
           FACT *fact, GEN p0)
{
  const long prec = nf_get_prec(nf);
  FP_t fp;
  pari_sp av;
  GEN L_jid = F->L_jid, Np0 = NULL;
  long Nsmall, Nfact, n = lg(L_jid);
  pari_timer T;

  if (DEBUGLEVEL)
  {
    timer_start(&T);
    err_printf("#### Look for %ld relations in %ld ideals (small_norm)\n",
               cache->end - cache->last, lg(L_jid)-1);
    if (p0) err_printf("Look in p0 = %Ps\n", vecslice(p0,1,4));
  }
  Nsmall = Nfact = 0;
  minim_alloc(lg(M), &fp.q, &fp.x, &fp.y, &fp.z, &fp.v);
  if (p0)
  {
    GEN n = pr_norm(p0);
    ulong e = maxuu(1,logint0(sqri(pr_norm(veclast(F->LP))), n, NULL));
    p0 = idealpow(nf, p0, utoi(e));
    Np0 = powiu(n,e);
  }
  for (av = avma; --n; set_avma(av))
  {
    long j = L_jid[n];
    GEN id = gel(F->LP, j), Nid;
    if (DEBUGLEVEL>1)
      err_printf("\n*** Ideal no %ld: %Ps\n", j, vecslice(id,1,4));
    if (p0)
    { Nid = mulii(Np0, pr_norm(id)); id = idealmul(nf, p0, id); }
    else
    { Nid = pr_norm(id); id = pr_hnf(nf, id);}
    if (Fincke_Pohst_ideal(cache, F, nf, M, id, Nid, fact, Nrelid, &fp,
                           NULL, prec, &Nsmall, &Nfact)) break;
  }
  if (DEBUGLEVEL && Nsmall)
  {
    if (DEBUGLEVEL == 1)
    { if (Nfact) err_printf("\n"); }
    else
      err_printf("  \nnb. fact./nb. small norm = %ld/%ld = %.3f\n",
                  Nfact,Nsmall,((double)Nfact)/Nsmall);
    if (timer_get(&T)>1) timer_printf(&T,"small_norm");
  }
}

static GEN
get_random_ideal(FB_t *F, GEN nf, GEN ex)
{
  long i, l = lg(ex);
  for (;;)
  {
    GEN I = NULL;
    for (i = 1; i < l; i++)
      if ((ex[i] = random_bits(RANDOM_BITS)))
      {
        GEN pr = gel(F->LP, F->subFB[i]), e = utoipos(ex[i]);
        I = I? idealmulpowprime(nf, I, pr, e): idealpow(nf, pr, e);
      }
    if (I && !ZM_isscalar(I,NULL)) return I; /* != (n)Z_K */
  }
}

static void
rnd_rel(RELCACHE_t *cache, FB_t *F, GEN nf, FACT *fact)
{
  pari_timer T;
  GEN L_jid = F->L_jid, M = nf_get_M(nf), R, NR;
  long i, l = lg(L_jid), prec = nf_get_prec(nf), Nfact = 0;
  RNDREL_t rr;
  FP_t fp;
  pari_sp av;

  if (DEBUGLEVEL) {
    timer_start(&T);
    err_printf("#### Look for %ld relations in %ld ideals (rnd_rel)\n",
               cache->end - cache->last, l-1);
  }
  rr.ex = cgetg(lg(F->subFB), t_VECSMALL);
  R = get_random_ideal(F, nf, rr.ex); /* random product from subFB */
  NR = ZM_det_triangular(R);
  minim_alloc(lg(M), &fp.q, &fp.x, &fp.y, &fp.z, &fp.v);
  for (av = avma, i = 1; i < l; i++, set_avma(av))
  { /* try P[j] * base */
    long j = L_jid[i];
    GEN P = gel(F->LP, j), Nid = mulii(NR, pr_norm(P));
    if (DEBUGLEVEL>1) err_printf("\n*** Ideal %ld: %Ps\n", j, vecslice(P,1,4));
    rr.jid = j;
    if (Fincke_Pohst_ideal(cache, F, nf, M, idealHNF_mul(nf, R, P), Nid, fact,
                           RND_REL_RELPID, &fp, &rr, prec, NULL, &Nfact)) break;
  }
  if (DEBUGLEVEL)
  {
    if (Nfact) err_printf("\n");
    if (timer_get(&T)>=0) timer_printf(&T,"rnd_rel");
  }
}

static GEN
automorphism_perms(GEN M, GEN auts, GEN cyclic, long r1, long r2, long N)
{
  long L = lgcols(M), lauts = lg(auts), lcyc = lg(cyclic), i, j, l, m;
  GEN Mt, perms = cgetg(lauts, t_VEC);
  pari_sp av;

  for (l = 1; l < lauts; l++) gel(perms, l) = cgetg(L, t_VECSMALL);
  av = avma;
  Mt = shallowtrans(gprec_w(M, LOWDEFAULTPREC));
  Mt = shallowconcat(Mt, conj_i(vecslice(Mt, r1+1, r1+r2)));
  for (l = 1; l < lcyc; l++)
  {
    GEN thiscyc = gel(cyclic, l), thisperm, perm, prev, Nt;
    long k = thiscyc[1];

    Nt = RgM_mul(shallowtrans(gel(auts, k)), Mt);
    perm = gel(perms, k);
    for (i = 1; i < L; i++)
    {
      GEN v = gel(Nt, i), minD;
      minD = gnorml2(gsub(v, gel(Mt, 1)));
      perm[i] = 1;
      for (j = 2; j <= N; j++)
      {
        GEN D = gnorml2(gsub(v, gel(Mt, j)));
        if (gcmp(D, minD) < 0) { minD = D; perm[i] = j >= L ? r2-j : j; }
      }
    }
    for (prev = perm, m = 2; m < lg(thiscyc); m++, prev = thisperm)
    {
      thisperm = gel(perms, thiscyc[m]);
      for (i = 1; i < L; i++)
      {
        long pp = labs(prev[i]);
        thisperm[i] = prev[i] < 0 ? -perm[pp] : perm[pp];
      }
    }
  }
  set_avma(av); return perms;
}

/* Determine the field automorphisms as matrices on the integral basis */
static GEN
automorphism_matrices(GEN nf, GEN *cycp)
{
  pari_sp av = avma;
  GEN auts = galoisconj(nf, NULL), mats, cyclic, cyclicidx;
  long nauts = lg(auts)-1, i, j, k, l;

  cyclic = cgetg(nauts+1, t_VEC);
  cyclicidx = zero_Flv(nauts);
  for (l = 1; l <= nauts; l++)
  {
    GEN aut = gel(auts, l);
    if (gequalX(aut)) { swap(gel(auts, l), gel(auts, nauts)); break; }
  }
  /* trivial automorphism is last */
  for (l = 1; l <= nauts; l++) gel(auts, l) = algtobasis(nf, gel(auts, l));
  /* Compute maximal cyclic subgroups */
  for (l = nauts; --l > 0; ) if (!cyclicidx[l])
  {
    GEN elt = gel(auts, l), aut = elt, cyc = cgetg(nauts+1, t_VECSMALL);
    cyc[1] = cyclicidx[l] = l; j = 1;
    do
    {
      elt = galoisapply(nf, elt, aut);
      for (k = 1; k <= nauts; k++) if (gequal(elt, gel(auts, k))) break;
      cyclicidx[k] = l; cyc[++j] = k;
    }
    while (k != nauts);
    setlg(cyc, j);
    gel(cyclic, l) = cyc;
  }
  for (i = j = 1; i < nauts; i++)
    if (cyclicidx[i] == i) cyclic[j++] = cyclic[i];
  setlg(cyclic, j);
  mats = cgetg(nauts, t_VEC);
  while (--j > 0)
  {
    GEN cyc = gel(cyclic, j);
    long id = cyc[1];
    GEN M, Mi, aut = gel(auts, id);

    gel(mats, id) = Mi = M = nfgaloismatrix(nf, aut);
    for (i = 2; i < lg(cyc); i++) gel(mats, cyc[i]) = Mi = ZM_mul(Mi, M);
  }
  gerepileall(av, 2, &mats, &cyclic);
  if (cycp) *cycp = cyclic;
  return mats;
}

/* vP a list of maximal ideals above the same p from idealprimedec: f(P/p) is
 * increasing; 1 <= j <= #vP; orbit a zc of length <= #vP; auts a vector of
 * automorphisms in ZM form.
 * Set orbit[i] = 1 for all vP[i], i >= j, in the orbit of pr = vP[j] wrt auts.
 * N.B.1 orbit need not be initialized to 0: useful to incrementally run
 * through successive orbits
 * N.B.2 i >= j, so primes with index < j will be missed; run incrementally
 * starting from j = 1 ! */
static void
pr_orbit_fill(GEN orbit, GEN auts, GEN vP, long j)
{
  GEN pr = gel(vP,j), gen = pr_get_gen(pr);
  long i, l = lg(auts), J = lg(orbit), f = pr_get_f(pr);
  orbit[j] = 1;
  for (i = 1; i < l; i++)
  {
    GEN g = ZM_ZC_mul(gel(auts,i), gen);
    long k;
    for (k = j+1; k < J; k++)
    {
      GEN prk = gel(vP,k);
      if (pr_get_f(prk) > f) break; /* f(P[k]) increases with k */
      /* don't check that e matches: (almost) always 1 ! */
      if (!orbit[k] && ZC_prdvd(g, prk)) { orbit[k] = 1; break; }
    }
  }
}
/* remark: F->KCZ changes if be_honest() fails */
static int
be_honest(FB_t *F, GEN nf, GEN auts, FACT *fact)
{
  long i, iz, nbtest;
  long lgsub = lg(F->subFB), KCZ0 = F->KCZ;
  long N = nf_get_degree(nf), prec = nf_get_prec(nf);
  GEN M = nf_get_M(nf);
  FP_t fp;
  pari_sp av;

  if (DEBUGLEVEL) {
    err_printf("Be honest for %ld primes from %ld to %ld\n", F->KCZ2 - F->KCZ,
               F->FB[ F->KCZ+1 ], F->FB[ F->KCZ2 ]);
  }
  minim_alloc(N+1, &fp.q, &fp.x, &fp.y, &fp.z, &fp.v);
  if (lg(auts) == 1) auts = NULL;
  av = avma;
  for (iz=F->KCZ+1; iz<=F->KCZ2; iz++, set_avma(av))
  {
    long p = F->FB[iz];
    GEN pr_orbit, P = gel(F->LV,p);
    long j, J = lg(P); /* > 1 */
    /* the P|p, NP > C2 are assumed in subgroup generated by FB + last P
     * with NP <= C2 is unramified --> check all but last */
    if (pr_get_e(gel(P,J-1)) == 1) J--;
    if (J == 1) continue;
    if (DEBUGLEVEL>1) err_printf("%ld ", p);
    pr_orbit = auts? zero_zv(J-1): NULL;
    for (j = 1; j < J; j++)
    {
      GEN Nid, id, id0;
      if (pr_orbit)
      {
        if (pr_orbit[j]) continue;
        /* discard all primes in automorphism orbit simultaneously */
        pr_orbit_fill(pr_orbit, auts, P, j);
      }
      id = id0 = pr_hnf(nf,gel(P,j));
      Nid = pr_norm(gel(P,j));
      for (nbtest=0;;)
      {
        if (Fincke_Pohst_ideal(NULL, F, nf, M, id, Nid, fact, 0, &fp,
                               NULL, prec, NULL, NULL)) break;
        if (++nbtest > maxtry_HONEST)
        {
          if (DEBUGLEVEL)
            pari_warn(warner,"be_honest() failure on prime %Ps\n", gel(P,j));
          return 0;
        }
        /* occurs at most once in the whole function */
        for (i = 1, id = id0; i < lgsub; i++)
        {
          long ex = random_bits(RANDOM_BITS);
          if (ex)
          {
            GEN pr = gel(F->LP, F->subFB[i]);
            id = idealmulpowprime(nf, id, pr, utoipos(ex));
          }
        }
        if (!equali1(gcoeff(id,N,N))) id = Q_primpart(id);
        if (expi(gcoeff(id,1,1)) > 100) id = idealred(nf, id);
        Nid = ZM_det_triangular(id);
      }
    }
    F->KCZ++; /* SUCCESS, "enlarge" factorbase */
  }
  F->KCZ = KCZ0; return gc_bool(av,1);
}

/* all primes with N(P) <= BOUND factor on factorbase ? */
void
bnftestprimes(GEN bnf, GEN BOUND)
{
  pari_sp av0 = avma, av;
  ulong count = 0;
  GEN auts, p, nf = bnf_get_nf(bnf), Vbase = bnf_get_vbase(bnf);
  GEN fb = gen_sort_shallow(Vbase, (void*)&cmp_prime_ideal, cmp_nodata);
  ulong pmax = pr_get_smallp(veclast(fb)); /*largest p in factorbase*/
  forprime_t S;
  FACT *fact;
  FB_t F;

  (void)recover_partFB(&F, Vbase, nf_get_degree(nf));
  fact = (FACT*)stack_malloc((F.KC+1)*sizeof(FACT));
  forprime_init(&S, gen_2, BOUND);
  auts = automorphism_matrices(nf, NULL);
  if (lg(auts) == 1) auts = NULL;
  av = avma;
  while (( p = forprime_next(&S) ))
  {
    GEN pr_orbit, vP;
    long j, J;

    if (DEBUGLEVEL == 1 && ++count > 1000)
    {
      err_printf("passing p = %Ps / %Ps\n", p, BOUND);
      count = 0;
    }
    set_avma(av);
    vP = idealprimedec_limit_norm(nf, p, BOUND);
    J = lg(vP);
    /* if last is unramified, all P|p in subgroup generated by FB: skip last */
    if (J > 1 && pr_get_e(gel(vP,J-1)) == 1) J--;
    if (J == 1) continue;
    if (DEBUGLEVEL>1) err_printf("*** p = %Ps\n",p);
    pr_orbit = auts? zero_zv(J-1): NULL;
    for (j = 1; j < J; j++)
    {
      GEN P = gel(vP,j);
      long k = 0;
      if (pr_orbit)
      {
        if (pr_orbit[j]) continue;
        /* discard all primes in automorphism orbit simultaneously */
        pr_orbit_fill(pr_orbit, auts, vP, j);
      }
      if (abscmpiu(p, pmax) > 0 || !(k = tablesearch(fb, P, &cmp_prime_ideal)))
        (void)SPLIT(&F, nf, pr_hnf(nf,P), Vbase, fact);
      if (DEBUGLEVEL>1)
      {
        err_printf("  Testing P = %Ps\n",P);
        if (k) err_printf("    #%ld in factor base\n",k);
        else err_printf("    is %Ps\n", isprincipal(bnf,P));
      }
    }
  }
  set_avma(av0);
}

/* A t_MAT of complex floats, in fact reals. Extract a submatrix B
 * whose columns are definitely nonzero, i.e. gexpo(A[j]) >= -2
 *
 * If possible precision problem (t_REAL 0 with large exponent), set
 * *precpb to 1 */
static GEN
clean_cols(GEN A, int *precpb)
{
  long l = lg(A), h, i, j, k;
  GEN B;
  *precpb = 0;
  if (l == 1) return A;
  h = lgcols(A);;
  B = cgetg(l, t_MAT);
  for (i = k = 1; i < l; i++)
  {
    GEN Ai = gel(A,i);
    int non0 = 0;
    for (j = 1; j < h; j++)
    {
      GEN c = gel(Ai,j);
      if (gexpo(c) >= -2)
      {
        if (gequal0(c)) *precpb = 1; else non0 = 1;
      }
    }
    if (non0) gel(B, k++) = Ai;
  }
  setlg(B, k); return B;
}

static long
compute_multiple_of_R_pivot(GEN X, GEN x0/*unused*/, long ix, GEN c)
{
  GEN x = gel(X,ix);
  long i, k = 0, ex = - (long)HIGHEXPOBIT, lx = lg(x);
  (void)x0;
  for (i=1; i<lx; i++)
    if (!c[i] && !gequal0(gel(x,i)))
    {
      long e = gexpo(gel(x,i));
      if (e > ex) { ex = e; k = i; }
    }
  return (k && ex > -32)? k: lx;
}

/* Ar = (log |sigma_i(u_j)|) for units (u_j) found so far;
 * RU = R1+R2 = target rank for unit matrix, after adding [1 x r1, 2 x r2];
 * N = field degree, need = unit rank defect;
 * L = NULL (prec problem) or B^(-1) * A with approximate rational entries
 * (as t_REAL), B a submatrix of A, with (probably) maximal rank RU */
static GEN
compute_multiple_of_R(GEN Ar, long RU, long N, long *pneed, long *bit, GEN *ptL)
{
  GEN T, d, mdet, Im_mdet, kR, L;
  long i, j, r, R1 = 2*RU - N;
  int precpb;
  pari_sp av = avma;

  if (RU == 1) { *ptL = zeromat(0, lg(Ar)-1); return gen_1; }

  if (DEBUGLEVEL) err_printf("\n#### Computing regulator multiple\n");
  mdet = clean_cols(Ar, &precpb);
  /* will cause precision to increase on later failure, but we may succeed! */
  *ptL = precpb? NULL: gen_1;
  T = cgetg(RU+1,t_COL);
  for (i=1; i<=R1; i++) gel(T,i) = gen_1;
  for (   ; i<=RU; i++) gel(T,i) = gen_2;
  mdet = shallowconcat(T, mdet); /* det(Span(mdet)) = N * R */

  /* could be using indexrank(), but need custom "get_pivot" function */
  d = RgM_pivots(mdet, NULL, &r, &compute_multiple_of_R_pivot);
  /* # of independent columns = target rank ? */
  if (lg(mdet)-1 - r != RU)
  {
    if (DEBUGLEVEL)
      err_printf("Units matrix target rank = %ld < %ld\n",lg(mdet)-1 - r, RU);
    *pneed = RU - (lg(mdet)-1-r); return gc_NULL(av);
  }

  Im_mdet = cgetg(RU+1, t_MAT); /* extract independent columns */
  /* N.B: d[1] = 1, corresponding to T above */
  gel(Im_mdet, 1) = T;
  for (i = j = 2; i <= RU; j++)
    if (d[j]) gel(Im_mdet, i++) = gel(mdet,j);

  /* integral multiple of R: the cols we picked form a Q-basis, they have an
   * index in the full lattice. First column is T */
  kR = divru(det2(Im_mdet), N);
  /* R > 0.2 uniformly */
  if (!signe(kR) || expo(kR) < -3)
  {
    if (DEBUGLEVEL) err_printf("Regulator is zero.\n");
    *pneed = 0; return gc_NULL(av);
  }
  d = det2(rowslice(vecslice(Im_mdet, 2, RU), 2, RU));
  setabssign(d); setabssign(kR);
  if (gexpo(gsub(d,kR)) - gexpo(d) > -20) { *ptL = NULL; return gc_NULL(av); }
  L = RgM_inv(Im_mdet);
  /* estimate # of correct bits in result */
  if (!L || (*bit = -gexpo(RgM_Rg_sub_shallow(RgM_mul(L,Im_mdet), gen_1))) < 16)
  { *ptL = NULL; return gc_NULL(av); }

  *ptL = RgM_mul(rowslice(L,2,RU), Ar); /* approximate rational entries */
  return gc_all(av,2, &kR, ptL);
}

/* leave small integer n as is, convert huge n to t_REAL (for readability) */
static GEN
i2print(GEN n)
{ return lgefint(n) <= DEFAULTPREC? n: itor(n,LOWDEFAULTPREC); }

static long
bad_check(GEN c)
{
  long ec = gexpo(c);
  if (DEBUGLEVEL) err_printf("\n ***** check = %.28Pg\n",c);
  /* safe check for c < 0.75 : avoid underflow in gtodouble() */
  if (ec < -1 || (ec == -1 && gtodouble(c) < 0.75)) return fupb_PRECI;
  /* safe check for c > 1.3 : avoid overflow */
  if (ec > 0 || (ec == 0 && gtodouble(c) > 1.3)) return fupb_RELAT;
  return fupb_NONE;
}
/* Input:
 * lambda = approximate rational entries: coords of units found so far on a
 * sublattice of maximal rank (sublambda)
 * *ptkR = regulator of sublambda = multiple of regulator of lambda
 * Compute R = true regulator of lambda.
 *
 * If c := Rz ~ 1, by Dirichlet's formula, then lambda is the full group of
 * units AND the full set of relations for the class group has been computed.
 * In fact z is a very rough approximation and we only expect 0.75 < Rz < 1.3
 *
 * Output: *ptkR = R, *ptL = numerator(units) (in terms of lambda) */
static long
compute_R(GEN lambda, GEN z, GEN *ptL, GEN *ptkR)
{
  pari_sp av = avma;
  long bit, r, reason, RU = lg(lambda) == 1? 1: lgcols(lambda);
  GEN L, H, D, den, R, c;

  *ptL = NULL;
  if (RU == 1) { *ptkR = gen_1; *ptL = lambda; return bad_check(z); }
  D = gmul2n(mpmul(*ptkR,z), 1); /* bound for denom(lambda) */
  if (expo(D) < 0 && rtodbl(D) < 0.95) return fupb_PRECI;
  L = bestappr(lambda,D);
  if (lg(L) == 1)
  {
    if (DEBUGLEVEL) err_printf("truncation error in bestappr\n");
    return fupb_PRECI;
  }
  den = Q_denom(L);
  if (mpcmp(den,D) > 0)
  {
    if (DEBUGLEVEL) err_printf("D = %Ps\nden = %Ps\n",D, i2print(den));
    return fupb_PRECI;
  }
  bit = -gexpo(gsub(L, lambda)); /* input accuracy */
  L = Q_muli_to_int(L, den);
  if (gexpo(L) + expi(den) > bit - 32)
  {
    if (DEBUGLEVEL) err_printf("dubious bestappr; den = %Ps\n", i2print(den));
    return fupb_PRECI;
  }
  H = ZM_hnf(L); r = lg(H)-1;
  if (!r || r != nbrows(H))
    R = gen_0; /* wrong rank */
  else
    R = gmul(*ptkR, gdiv(ZM_det_triangular(H), powiu(den, r)));
  /* R = tentative regulator; regulator > 0.2 uniformly */
  if (gexpo(R) < -3) {
    if (DEBUGLEVEL) err_printf("\n#### Tentative regulator: %.28Pg\n", R);
    return gc_long(av, fupb_PRECI);
  }
  c = gmul(R,z); /* should be n (= 1 if we are done) */
  if (DEBUGLEVEL) err_printf("\n#### Tentative regulator: %.28Pg\n", R);
  if ((reason = bad_check(c))) return gc_long(av, reason);
  *ptkR = R; *ptL = L; return fupb_NONE;
}
static GEN
get_clg2(GEN cyc, GEN Ga, GEN C, GEN Ur, GEN Ge, GEN M1, GEN M2)
{
  GEN GD = gsub(act_arch(M1, C), diagact_arch(cyc, Ga));
  GEN ga = gsub(act_arch(M2, C), act_arch(Ur, Ga));
  return mkvecn(6, Ur, ga, GD, Ge, M1, M2);
}
/* compute class group (clg1) + data for isprincipal (clg2) */
static GEN
class_group_gen(GEN nf,GEN W,GEN C,GEN Vbase,long prec, GEN *pclg2)
{
  GEN M1, M2, z, G, Ga, Ge, cyc, X, Y, D, U, V, Ur, Ui, Uir;
  long j, l;

  D = ZM_snfall(W,&U,&V); /* UWV=D, D diagonal, G = g Ui (G=new gens, g=old) */
  Ui = ZM_inv(U, NULL);
  l = lg(D); cyc = cgetg(l, t_VEC); /* elementary divisors */
  for (j = 1; j < l; j++)
  {
    gel(cyc,j) = gcoeff(D,j,j); /* strip useless components */
    if (is_pm1(gel(cyc,j))) break;
  }
  l = j;
  Ur  = ZM_hnfdivrem(U, D, &Y);
  Uir = ZM_hnfdivrem(Ui,W, &X);
 /* {x} = logarithmic embedding of x (arch. component)
  * NB: [J,z] = idealred(I) --> I = y J, with {y} = - z
  * G = g Uir - {Ga},  Uir = Ui + WX
  * g = G Ur  - {ga},  Ur  = U + DY */
  G = cgetg(l,t_VEC);
  Ga= cgetg(l,t_MAT);
  Ge= cgetg(l,t_COL);
  z = init_famat(NULL);
  for (j = 1; j < l; j++)
  {
    GEN I = genback(z, nf, Vbase, gel(Uir,j));
    gel(G,j) = gel(I,1); /* generator, order cyc[j] */
    gel(Ge,j)= gel(I,2);
    gel(Ga,j)= nf_cxlog(nf, gel(I,2), prec);
    if (!gel(Ga,j)) pari_err_PREC("class_group_gen");
  }
  /* {ga} = {GD}Y + G U - g = {GD}Y - {Ga} U + gW X U
                            = gW (X Ur + V Y) - {Ga}Ur */
  M2 = ZM_add(ZM_mul(X,Ur), ZM_mul(V,Y));
  setlg(cyc,l); setlg(V,l); setlg(D,l);
  /* G D =: {GD} = g (Ui + W X) D - {Ga}D = g W (V + X D) - {Ga}D
   * NB: Ui D = W V. gW is given by (first l-1 cols of) C */
  M1 = ZM_add(V, ZM_mul(X,D));
  *pclg2 = get_clg2(cyc, Ga, C, Ur, Ge, M1, M2);
  return mkvec3(ZV_prod(cyc), cyc, G);
}

/* compute principal ideals corresponding to (gen[i]^cyc[i]) */
static GEN
makecycgen(GEN bnf)
{
  GEN cyc = bnf_get_cyc(bnf), gen = bnf_get_gen(bnf), nf = bnf_get_nf(bnf);
  GEN h, y, GD = bnf_get_GD(bnf), W = bnf_get_W(bnf); /* HNF */
  GEN Sunits = bnf_get_sunits(bnf);
  GEN X = Sunits? gel(Sunits,1): NULL, C = Sunits? gel(Sunits,3): NULL;
  long e, i, l;

  if (DEBUGLEVEL) pari_warn(warner,"completing bnf (building cycgen)");
  h = cgetg_copy(gen, &l);
  for (i = 1; i < l; i++)
  {
    GEN gi = gel(gen,i), ci = gel(cyc,i);
    if (X && equalii(ci, gcoeff(W,i,i)))
    {
      long j;
      for (j = i+1; j < l; j++)
        if (signe(gcoeff(W,i,j))) break;
      if (j == i) { gel(h,i) = mkmat2(X, gel(C,i)); continue; }
    }
    if (abscmpiu(ci, 5) < 0)
    {
      GEN N = ZM_det_triangular(gi);
      y = isprincipalarch(bnf,gel(GD,i), N, ci, gen_1, &e);
      if (y && fact_ok(nf,y,NULL,mkvec(gi),mkvec(ci)))
      {
        gel(h,i) = to_famat_shallow(y,gen_1);
        continue;
      }
    }
    y = isprincipalfact(bnf, NULL, mkvec(gi), mkvec(ci), nf_GENMAT|nf_FORCE);
    gel(h,i) = gel(y,2);
  }
  return h;
}

static GEN
get_y(GEN bnf, GEN W, GEN B, GEN C, GEN pFB, long j)
{
  GEN y, nf  = bnf_get_nf(bnf);
  long e, lW = lg(W)-1;
  GEN ex = (j<=lW)? gel(W,j): gel(B,j-lW);
  GEN P = (j<=lW)? NULL: gel(pFB,j);
  if (C)
  { /* archimedean embeddings known: cheap trial */
    GEN Nx = get_norm_fact_primes(pFB, ex, P);
    y = isprincipalarch(bnf,gel(C,j), Nx,gen_1, gen_1, &e);
    if (y && fact_ok(nf,y,P,pFB,ex)) return y;
  }
  y = isprincipalfact_or_fail(bnf, P, pFB, ex);
  return typ(y) == t_INT? y: gel(y,2);
}
/* compute principal ideals corresponding to bnf relations */
static GEN
makematal(GEN bnf)
{
  GEN W = bnf_get_W(bnf), B = bnf_get_B(bnf), C = bnf_get_C(bnf);
  GEN pFB, ma, retry;
  long lma, j, prec = 0;

  if (DEBUGLEVEL) pari_warn(warner,"completing bnf (building matal)");
  lma=lg(W)+lg(B)-1;
  pFB = bnf_get_vbase(bnf);
  ma = cgetg(lma,t_VEC);
  retry = vecsmalltrunc_init(lma);
  for (j=lma-1; j>0; j--)
  {
    pari_sp av = avma;
    GEN y = get_y(bnf, W, B, C, pFB, j);
    if (typ(y) == t_INT)
    {
      long E = itos(y);
      if (DEBUGLEVEL>1) err_printf("\n%ld done later at prec %ld\n",j,E);
      set_avma(av);
      vecsmalltrunc_append(retry, j);
      if (E > prec) prec = E;
    }
    else
    {
      if (DEBUGLEVEL>1) err_printf("%ld ",j);
      gel(ma,j) = gerepileupto(av,y);
    }
  }
  if (prec)
  {
    long k, l = lg(retry);
    GEN y, nf = bnf_get_nf(bnf);
    if (DEBUGLEVEL) pari_warn(warnprec,"makematal",prec);
    nf = nfnewprec_shallow(nf,prec);
    bnf = Buchall(nf, nf_FORCE, prec);
    if (DEBUGLEVEL) err_printf("makematal, adding missing entries:");
    for (k=1; k<l; k++)
    {
      pari_sp av = avma;
      long j = retry[k];
      y = get_y(bnf,W,B,NULL, pFB, j);
      if (typ(y) == t_INT) pari_err_PREC("makematal");
      if (DEBUGLEVEL>1) err_printf("%ld ",j);
      gel(ma,j) = gerepileupto(av,y);
    }
  }
  if (DEBUGLEVEL>1) err_printf("\n");
  return ma;
}

enum { MATAL = 1, CYCGEN, UNITS };
GEN
bnf_build_cycgen(GEN bnf)
{ return obj_checkbuild(bnf, CYCGEN, &makecycgen); }
GEN
bnf_build_matalpha(GEN bnf)
{ return obj_checkbuild(bnf, MATAL, &makematal); }
GEN
bnf_build_units(GEN bnf)
{ return obj_checkbuild(bnf, UNITS, &makeunits); }

/* return fu in compact form if available; in terms of a fixed basis
 * of S-units */
GEN
bnf_compactfu_mat(GEN bnf)
{
  GEN X, U, SUnits = bnf_get_sunits(bnf);
  if (!SUnits) return NULL;
  X = gel(SUnits,1);
  U = gel(SUnits,2); ZM_remove_unused(&U, &X);
  return mkvec2(X, U);
}
/* return fu in compact form if available; individually as famat */
GEN
bnf_compactfu(GEN bnf)
{
  GEN fu, X, U, SUnits = bnf_get_sunits(bnf);
  long i, l;
  if (!SUnits) return NULL;
  X = gel(SUnits,1);
  U = gel(SUnits,2); l = lg(U); fu = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(fu,i) = famat_remove_trivial(mkmat2(X, gel(U,i)));
  return fu;
}
/* return expanded fu if available */
GEN
bnf_has_fu(GEN bnf)
{
  GEN fu = obj_check(bnf, UNITS);
  if (fu) return vecsplice(fu, 1);
  fu = bnf_get_fu_nocheck(bnf);
  return (typ(fu) == t_MAT)? NULL: fu;
}
/* return expanded fu if available; build if cheap */
GEN
bnf_build_cheapfu(GEN bnf)
{
  GEN fu, SUnits;
  if ((fu = bnf_has_fu(bnf))) return fu;
  if ((SUnits = bnf_get_sunits(bnf)))
  {
    pari_sp av = avma;
    long e = gexpo(real_i(bnf_get_logfu(bnf)));
    set_avma(av); if (e < 13) return vecsplice(bnf_build_units(bnf), 1);
  }
  return NULL;
}

static GEN
get_regulator(GEN A)
{
  pari_sp av = avma;
  GEN R;

  if (lg(A) == 1) return gen_1;
  R = det( rowslice(real_i(A), 1, lgcols(A)-2) );
  setabssign(R); return gerepileuptoleaf(av, R);
}

/* return corrected archimedian components for elts of x (vector)
 * (= log(sigma_i(x)) - log(|Nx|) / [K:Q]) */
static GEN
get_archclean(GEN nf, GEN x, long prec, int units)
{
  long k, N, l = lg(x);
  GEN M = cgetg(l, t_MAT);

  if (l == 1) return M;
  N = nf_get_degree(nf);
  for (k = 1; k < l; k++)
  {
    pari_sp av = avma;
    GEN c = nf_cxlog(nf, gel(x,k), prec);
    if (!c || (!units && !(c = cleanarch(c, N, NULL,prec)))) return NULL;
    gel(M,k) = gerepilecopy(av, c);
  }
  return M;
}
static void
Sunits_archclean(GEN nf, GEN Sunits, GEN *pmun, GEN *pC, long prec)
{
  GEN ipi, M, X = gel(Sunits,1), U = gel(Sunits,2), G = gel(Sunits,3);
  long k, N = nf_get_degree(nf), l = lg(X);

  M = cgetg(l, t_MAT);
  for (k = 1; k < l; k++)
    if (!(gel(M,k) = nf_cxlog(nf, gel(X,k), prec))) return;
  ipi = invr(mppi(prec));
  *pmun = cleanarch(RgM_ZM_mul(M, U), N, ipi, prec); /* not cleanarchunit ! */
  if (*pmun) *pC = cleanarch(RgM_ZM_mul(M, G), N, ipi, prec);
}

GEN
bnfnewprec_shallow(GEN bnf, long prec)
{
  GEN nf0 = bnf_get_nf(bnf), nf, v, fu, matal, y, A, C;
  GEN Sunits = bnf_get_sunits(bnf), Ur, Ga, Ge, M1, M2;
  long r1, r2, prec0 = prec;

  nf_get_sign(nf0, &r1, &r2);
  if (Sunits)
  {
    fu = matal = NULL;
    prec += nbits2extraprec(gexpo(Sunits));
  }
  else
  {
    fu = bnf_build_units(bnf);
    fu = vecslice(fu, 2, lg(fu)-1);
    if (r1 + r2 > 1) {
      long e = gexpo(bnf_get_logfu(bnf)) + 1 - TWOPOTBITS_IN_LONG;
      if (e >= 0) prec += nbits2extraprec(e);
    }
    matal = bnf_build_matalpha(bnf);
  }

  if (DEBUGLEVEL && prec0 != prec) pari_warn(warnprec,"bnfnewprec",prec);
  for(C = NULL;;)
  {
    pari_sp av = avma;
    nf = nfnewprec_shallow(nf0,prec);
    if (Sunits)
      Sunits_archclean(nf, Sunits, &A, &C, prec);
    else
    {
      A = get_archclean(nf, fu, prec, 1);
      if (A) C = get_archclean(nf, matal, prec, 0);
    }
    if (C) break;
    set_avma(av); prec = precdbl(prec);
    if (DEBUGLEVEL) pari_warn(warnprec,"bnfnewprec(extra)",prec);
  }
  y = leafcopy(bnf);
  gel(y,3) = A;
  gel(y,4) = C;
  gel(y,7) = nf;
  gel(y,8) = v = leafcopy(gel(bnf,8));
  gel(v,2) = get_regulator(A);
  v = gel(bnf,9);
  if (lg(v) < 7) pari_err_TYPE("bnfnewprec [obsolete bnf format]", bnf);
  Ur = gel(v,1);
  Ge = gel(v,4);
  Ga = nfV_cxlog(nf, Ge, prec);
  M1 = gel(v,5);
  M2 = gel(v,6);
  gel(y,9) = get_clg2(bnf_get_cyc(bnf), Ga, C, Ur, Ge, M1, M2);
  return y;
}
GEN
bnfnewprec(GEN bnf, long prec)
{
  pari_sp av = avma;
  return gerepilecopy(av, bnfnewprec_shallow(checkbnf(bnf), prec));
}

GEN
bnrnewprec_shallow(GEN bnr, long prec)
{
  GEN y = cgetg(7,t_VEC);
  long i;
  gel(y,1) = bnfnewprec_shallow(bnr_get_bnf(bnr), prec);
  for (i=2; i<7; i++) gel(y,i) = gel(bnr,i);
  return y;
}
GEN
bnrnewprec(GEN bnr, long prec)
{
  GEN y = cgetg(7,t_VEC);
  long i;
  checkbnr(bnr);
  gel(y,1) = bnfnewprec(bnr_get_bnf(bnr), prec);
  for (i=2; i<7; i++) gel(y,i) = gcopy(gel(bnr,i));
  return y;
}

static GEN
buchall_end(GEN nf,GEN res, GEN clg2, GEN W, GEN B, GEN A, GEN C,GEN Vbase)
{
  GEN z = obj_init(9, 3);
  gel(z,1) = W;
  gel(z,2) = B;
  gel(z,3) = A;
  gel(z,4) = C;
  gel(z,5) = Vbase;
  gel(z,6) = gen_0;
  gel(z,7) = nf;
  gel(z,8) = res;
  gel(z,9) = clg2;
  return z;
}

GEN
bnfinit0(GEN P, long flag, GEN data, long prec)
{
  double c1 = 0., c2 = 0.;
  long fl, relpid = BNF_RELPID;

  if (data)
  {
    long lx = lg(data);
    if (typ(data) != t_VEC || lx > 5) pari_err_TYPE("bnfinit",data);
    switch(lx)
    {
      case 4: relpid = itos(gel(data,3));
      case 3: c2 = gtodouble(gel(data,2));
      case 2: c1 = gtodouble(gel(data,1));
    }
  }
  switch(flag)
  {
    case 2:
    case 0: fl = 0; break;
    case 1: fl = nf_FORCE; break;
    default: pari_err_FLAG("bnfinit");
      return NULL; /* LCOV_EXCL_LINE */
  }
  return Buchall_param(P, c1, c2, relpid, fl, prec);
}
GEN
Buchall(GEN P, long flag, long prec)
{ return Buchall_param(P, 0., 0., BNF_RELPID, flag & nf_FORCE, prec); }

static GEN
Buchall_deg1(GEN nf)
{
  GEN v = cgetg(1,t_VEC), m = cgetg(1,t_MAT);
  GEN res, W, A, B, C, Vbase = cgetg(1,t_COL);
  GEN fu = v, R = gen_1, zu = mkvec2(gen_2, gen_m1);
  GEN clg1 = mkvec3(gen_1,v,v), clg2 = mkvecn(6, m,m,m,v,m,m);

  W = A = B = C = m; res = mkvec5(clg1, R, gen_1, zu, fu);
  return buchall_end(nf,res,clg2,W,B,A,C,Vbase);
}

/* return (small set of) indices of columns generating the same lattice as x.
 * Assume HNF(x) is inexpensive (few rows, many columns).
 * Dichotomy approach since interesting columns may be at the very end */
GEN
extract_full_lattice(GEN x)
{
  long dj, j, k, l = lg(x);
  GEN h, h2, H, v;

  if (l < 200) return NULL; /* not worth it */

  v = vecsmalltrunc_init(l);
  H = ZM_hnf(x);
  h = cgetg(1, t_MAT);
  dj = 1;
  for (j = 1; j < l; )
  {
    pari_sp av = avma;
    long lv = lg(v);

    for (k = 0; k < dj; k++) v[lv+k] = j+k;
    setlg(v, lv + dj);
    h2 = ZM_hnf(vecpermute(x, v));
    if (ZM_equal(h, h2))
    { /* these dj columns can be eliminated */
      set_avma(av); setlg(v, lv);
      j += dj;
      if (j >= l) break;
      dj <<= 1;
      if (j + dj >= l) { dj = (l - j) >> 1; if (!dj) dj = 1; }
    }
    else if (dj > 1)
    { /* at least one interesting column, try with first half of this set */
      set_avma(av); setlg(v, lv);
      dj >>= 1; /* > 0 */
    }
    else
    { /* this column should be kept */
      if (ZM_equal(h2, H)) break;
      h = h2; j++;
    }
  }
  return v;
}

static void
init_rel(RELCACHE_t *cache, FB_t *F, long add_need)
{
  const long n = F->KC + add_need; /* expected # of needed relations */
  long i, j, k, p;
  GEN c, P;
  GEN R;

  if (DEBUGLEVEL) err_printf("KCZ = %ld, KC = %ld, n = %ld\n", F->KCZ,F->KC,n);
  reallocate(cache, 10*n + 50); /* make room for lots of relations */
  cache->chk = cache->base;
  cache->end = cache->base + n;
  cache->relsup = add_need;
  cache->last = cache->base;
  cache->missing = lg(cache->basis) - 1;
  for (i = 1; i <= F->KCZ; i++)
  { /* trivial relations (p) = prod P^e */
    p = F->FB[i]; P = gel(F->LV,p);
    if (!isclone(P)) continue;

    /* all prime divisors in FB */
    c = zero_Flv(F->KC); k = F->iLP[p];
    R = c; c += k;
    for (j = lg(P)-1; j; j--) c[j] = pr_get_e(gel(P,j));
    add_rel(cache, F, R, k+1, pr_get_p(gel(P,1)), 0);
  }
}

/* Let z = \zeta_n in nf. List of not-obviously-dependent generators for
 * cyclotomic units modulo torsion in Q(z) [independent when n a prime power]:
 * - z^a - 1,  n/(a,n) not a prime power, a \nmid n unless a=1,  1 <= a < n/2
 * - (Z^a - 1)/(Z - 1),  p^k || n, Z = z^{n/p^k}, (p,a) = 1, 1 < a <= (p^k-1)/2
 */
GEN
nfcyclotomicunits(GEN nf, GEN zu)
{
  long n = itos(gel(zu, 1)), n2, lP, i, a;
  GEN z, fa, P, E, L, mz, powz;
  if (n <= 6) return cgetg(1, t_VEC);

  z = algtobasis(nf,gel(zu, 2));
  if ((n & 3) == 2) { n = n >> 1; z = ZC_neg(z); } /* ensure n != 2 (mod 4) */
  n2 = n/2;
  mz = zk_multable(nf, z); /* multiplication by z */
  powz = cgetg(n2, t_VEC); gel(powz,1) = z;
  for (i = 2; i < n2; i++) gel(powz,i) = ZM_ZC_mul(mz, gel(powz,i-1));
  /* powz[i] = z^i */

  L = vectrunc_init(n);
  fa = factoru(n);
  P = gel(fa,1); lP = lg(P);
  E = gel(fa,2);
  for (i = 1; i < lP; i++)
  { /* second kind */
    long p = P[i], k = E[i], pk = upowuu(p,k), pk2 = (pk-1) / 2;
    GEN u = gen_1;
    for (a = 2; a <= pk2; a++)
    {
      u = nfadd(nf, u, gel(powz, (n/pk) * (a-1))); /* = (Z^a-1)/(Z-1) */
      if (a % p) vectrunc_append(L, u);
    }
  }
  if (lP > 2) for (a = 1; a < n2; a++)
  { /* first kind, when n not a prime power */
    ulong p;
    if (a > 1 && (n % a == 0 || uisprimepower(n/ugcd(a,n), &p))) continue;
    vectrunc_append(L, nfadd(nf, gel(powz, a), gen_m1));
  }
  return L;
}
static void
add_cyclotomic_units(GEN nf, GEN zu, RELCACHE_t *cache, FB_t *F)
{
  pari_sp av = avma;
  GEN L = nfcyclotomicunits(nf, zu);
  long i, l = lg(L);
  if (l > 1)
  {
    GEN R = zero_Flv(F->KC);
    for(i = 1; i < l; i++) add_rel(cache, F, R, F->KC+1, gel(L,i), 0);
  }
  set_avma(av);
}

static GEN
trim_list(FB_t *F)
{
  pari_sp av = avma;
  GEN v, L_jid = F->L_jid, minidx = F->minidx, present = zero_Flv(F->KC);
  long i, j, imax = minss(lg(L_jid), F->KC + 1);

  v = cgetg(imax, t_VECSMALL);
  for (i = j = 1; i < imax; i++)
  {
    long k = minidx[ L_jid[i] ];
    if (!present[k]) { v[j++] = L_jid[i]; present[k] = 1; }
  }
  setlg(v, j); return gerepileuptoleaf(av, v);
}

static void
try_elt(RELCACHE_t *cache, FB_t *F, GEN nf, GEN x, FACT *fact)
{
  pari_sp av = avma;
  GEN R, Nx;
  long nz, tx = typ(x);

  if (tx == t_INT || tx == t_FRAC) return;
  if (tx != t_COL) x = algtobasis(nf, x);
  if (RgV_isscalar(x)) return;
  x = Q_primpart(x);
  Nx = nfnorm(nf, x);
  if (!can_factor(F, nf, NULL, x, Nx, fact)) return;

  /* smooth element */
  R = set_fact(F, fact, NULL, &nz);
  /* make sure we get maximal rank first, then allow all relations */
  (void) add_rel(cache, F, R, nz, x, 0);
  set_avma(av);
}

static void
matenlarge(GEN C, long h)
{
  GEN _0 = zerocol(h);
  long i;
  for (i = lg(C); --i; ) gel(C,i) = shallowconcat(gel(C,i), _0);
}

/* E = floating point embeddings */
static GEN
matbotidembs(RELCACHE_t *cache, GEN E)
{
  long w = cache->last - cache->chk, h = cache->last - cache->base;
  long j, d = h - w, hE = nbrows(E);
  GEN y = cgetg(w+1,t_MAT), _0 = zerocol(h);
  for (j = 1; j <= w; j++)
  {
    GEN c = shallowconcat(gel(E,j), _0);
    if (d + j >= 1) gel(c, d + j + hE) = gen_1;
    gel(y,j) = c;
  }
  return y;
}
static GEN
matbotid(RELCACHE_t *cache)
{
  long w = cache->last - cache->chk, h = cache->last - cache->base;
  long j, d = h - w;
  GEN y = cgetg(w+1,t_MAT);
  for (j = 1; j <= w; j++)
  {
    GEN c = zerocol(h);
    if (d + j >= 1) gel(c, d + j) = gen_1;
    gel(y,j) = c;
  }
  return y;
}

static long
myprecdbl(long prec, GEN C)
{
  long p = prec2nbits(prec) < 1280? precdbl(prec): (long)(prec * 1.5);
  if (C) p = maxss(p, minss(3*p, prec + nbits2extraprec(gexpo(C))));
  return p;
}

static GEN
_nfnewprec(GEN nf, long prec, long *isclone)
{
  GEN NF = gclone(nfnewprec_shallow(nf, prec));
  if (*isclone) gunclone(nf);
  *isclone = 1; return NF;
}

/* Nrelid = nb relations per ideal, possibly 0. If flag is set, keep data in
 * algebraic form. */
GEN
Buchall_param(GEN P, double cbach, double cbach2, long Nrelid, long flag, long prec)
{
  pari_timer T;
  pari_sp av0 = avma, av, av2;
  long PREC, N, R1, R2, RU, low, high, LIMC0, LIMC, LIMC2, LIMCMAX, zc, i;
  long LIMres, bit = 0, flag_nfinit = 0;
  long nreldep, sfb_trials, need, old_need, precdouble = 0, TRIES = 0;
  long nfisclone = 0;
  long done_small, small_fail, fail_limit, squash_index, small_norm_prec;
  double LOGD, LOGD2, lim;
  GEN computed = NULL, fu = NULL, zu, nf, M_sn, D, A, W, R, h, Ce, PERM;
  GEN small_multiplier, auts, cyclic, embs, SUnits;
  GEN res, L, invhr, B, C, lambda, dep, clg1, clg2, Vbase;
  const char *precpb = NULL;
  REL_t *old_cache = NULL;
  nfmaxord_t nfT;
  RELCACHE_t cache;
  FB_t F;
  GRHcheck_t GRHcheck;
  FACT *fact;

  if (DEBUGLEVEL) timer_start(&T);
  P = get_nfpol(P, &nf);
  if (nf)
    D = nf_get_disc(nf);
  else
  {
    nfinit_basic(&nfT, P);
    D = nfT.dK;
    if (!ZX_is_monic(nfT.T0))
    {
      pari_warn(warner,"nonmonic polynomial in bnfinit, using polredbest");
      flag_nfinit = nf_RED;
    }
  }
  PREC = maxss(DEFAULTPREC, prec);
  N = degpol(P);
  if (N <= 1)
  {
    if (!nf) nf = nfinit_complete(&nfT, flag_nfinit, PREC);
    return gerepilecopy(av0, Buchall_deg1(nf));
  }
  D = absi_shallow(D);
  LOGD = dbllog2(D) * M_LN2;
  LOGD2 = LOGD*LOGD;
  LIMCMAX = (long)(4.*LOGD2);
  /* In small_norm, LLL reduction produces v0 in I such that
   *     T2(v0) <= (4/3)^((n-1)/2) NI^(2/n) disc(K)^(1/n)
   * We consider v with T2(v) <= BMULT * T2(v0)
   * Hence Nv <= ((4/3)^((n-1)/2) * BMULT / n)^(n/2) NI sqrt(disc(K)).
   * NI <= LIMCMAX^2 */
  if (nf) PREC = maxss(PREC, nf_get_prec(nf));
  PREC = maxss(PREC, nbits2prec((long)(LOGD2 * 0.02) + N*N));
  if (DEBUGLEVEL) err_printf("PREC = %ld\n", PREC);
  small_norm_prec = nbits2prec( BITS_IN_LONG +
    (N/2. * ((N-1)/2.*log(4./3) + log(8/(double)N))
     + 2*log((double) LIMCMAX) + LOGD/2) / M_LN2 ); /*enough to compute norms*/
  if (small_norm_prec > PREC) PREC = small_norm_prec;
  if (!nf)
    nf = nfinit_complete(&nfT, flag_nfinit, PREC);
  else if (nf_get_prec(nf) < PREC)
    nf = nfnewprec_shallow(nf, PREC);
  M_sn = nf_get_M(nf);
  if (PREC > small_norm_prec) M_sn = gprec_w(M_sn, small_norm_prec);

  zu = nfrootsof1(nf);
  gel(zu,2) = nf_to_scalar_or_alg(nf, gel(zu,2));

  nf_get_sign(nf, &R1, &R2); RU = R1+R2;
  auts = automorphism_matrices(nf, &cyclic);
  F.embperm = automorphism_perms(nf_get_M(nf), auts, cyclic, R1, R2, N);
  if (DEBUGLEVEL)
  {
    timer_printf(&T, "nfinit & nfrootsof1");
    err_printf("%s bnf: R1 = %ld, R2 = %ld\nD = %Ps\n",
               flag? "Algebraic": "Floating point", R1,R2, D);
  }
  if (LOGD < 20.)
  { /* tiny disc, Minkowski may be smaller than Bach */
    lim = exp(-N + R2 * log(4/M_PI) + LOGD/2) * sqrt(2*M_PI*N);
    if (lim < 3) lim = 3;
  }
  else /* to be ignored */
    lim = -1;
  if (cbach > 12.) {
    if (cbach2 < cbach) cbach2 = cbach;
    cbach = 12.;
  }
  if (cbach < 0.)
    pari_err_DOMAIN("Buchall","Bach constant","<",gen_0,dbltor(cbach));

  cache.base = NULL; F.subFB = NULL; F.LP = NULL; SUnits = Ce = NULL;
  init_GRHcheck(&GRHcheck, N, R1, LOGD);
  high = low = LIMC0 = maxss((long)(cbach2*LOGD2), 1);
  while (!GRHchk(nf, &GRHcheck, high)) { low = high; high *= 2; }
  while (high - low > 1)
  {
    long test = (low+high)/2;
    if (GRHchk(nf, &GRHcheck, test)) high = test; else low = test;
  }
  LIMC2 = (high == LIMC0+1 && GRHchk(nf, &GRHcheck, LIMC0))? LIMC0: high;
  if (LIMC2 > LIMCMAX) LIMC2 = LIMCMAX;
  /* Assuming GRH, {P, NP <= LIMC2} generate Cl(K) */
  if (DEBUGLEVEL) err_printf("LIMC2 = %ld\n", LIMC2);
  LIMC0 = (long)(cbach*LOGD2); /* initial value for LIMC */
  LIMC = cbach? LIMC0: LIMC2; /* use {P, NP <= LIMC} as a factorbase */
  LIMC = maxss(LIMC, nthideal(&GRHcheck, nf, N));
  if (DEBUGLEVEL) timer_printf(&T, "computing Bach constant");
  LIMres = primeneeded(N, R1, R2, LOGD);
  cache_prime_dec(&GRHcheck, LIMres, nf);
  /* invhr ~ 2^r1 (2pi)^r2 / sqrt(D) w * Res(zeta_K, s=1) = 1 / hR */
  invhr = gmul(gdiv(gmul2n(powru(mppi(DEFAULTPREC), R2), RU),
              mulri(gsqrt(D,DEFAULTPREC),gel(zu,1))),
              compute_invres(&GRHcheck, LIMres));
  if (DEBUGLEVEL) timer_printf(&T, "computing inverse of hR");
  av = avma;

START:
  if (DEBUGLEVEL) timer_start(&T);
  if (TRIES) LIMC = bnf_increase_LIMC(LIMC,LIMCMAX);
  if (DEBUGLEVEL && LIMC > LIMC0)
    err_printf("%s*** Bach constant: %f\n", TRIES?"\n":"", LIMC/LOGD2);
  if (cache.base)
  {
    REL_t *rel;
    for (i = 1, rel = cache.base + 1; rel < cache.last; rel++)
      if (rel->m) i++;
    computed = cgetg(i, t_VEC);
    for (i = 1, rel = cache.base + 1; rel < cache.last; rel++)
      if (rel->m) gel(computed, i++) = rel->m;
    computed = gclone(computed); delete_cache(&cache);
  }
  TRIES++; set_avma(av);
  if (F.LP) delete_FB(&F);
  if (LIMC2 < LIMC) LIMC2 = LIMC;
  if (DEBUGLEVEL) { err_printf("LIMC = %ld, LIMC2 = %ld\n",LIMC,LIMC2); }

  FBgen(&F, nf, N, LIMC, LIMC2, &GRHcheck);
  if (!F.KC) goto START;
  av = avma;
  subFBgen(&F,auts,cyclic,lim < 0? LIMC2: mindd(lim,LIMC2),MINSFB);
  if (lg(F.subFB) == 1) goto START;
  if (DEBUGLEVEL)
    timer_printf(&T, "factorbase (#subFB = %ld) and ideal permutations",
                     lg(F.subFB)-1);

  fact = (FACT*)stack_malloc((F.KC+1)*sizeof(FACT));
  PERM = leafcopy(F.perm); /* to be restored in case of precision increase */
  cache.basis = zero_Flm_copy(F.KC,F.KC);
  small_multiplier = zero_Flv(F.KC);
  done_small = small_fail = squash_index = zc = sfb_trials = nreldep = 0;
  fail_limit = F.KC + 1;
  W = A = R = NULL;
  av2 = avma;
  init_rel(&cache, &F, RELSUP + RU-1);
  old_need = need = cache.end - cache.last;
  add_cyclotomic_units(nf, zu, &cache, &F);
  if (DEBUGLEVEL) err_printf("\n");
  cache.end = cache.last + need;

  if (computed)
  {
    for (i = 1; i < lg(computed); i++)
      try_elt(&cache, &F, nf, gel(computed, i), fact);
    gunclone(computed);
    if (DEBUGLEVEL && i > 1)
      timer_printf(&T, "including already computed relations");
    need = 0;
  }

  do
  {
    GEN Ar, C0;
    do
    {
      pari_sp av4 = avma;
      if (need > 0)
      {
        long oneed = cache.end - cache.last;
        /* Test below can be true if small_norm did not find enough linearly
         * dependent relations */
        if (need < oneed) need = oneed;
        pre_allocate(&cache, need+lg(auts)-1+(R ? lg(W)-1 : 0));
        cache.end = cache.last + need;
        F.L_jid = trim_list(&F);
      }
      if (need > 0 && Nrelid > 0 && (done_small <= F.KC+1 || A) &&
          small_fail <= fail_limit &&
          cache.last < cache.base + 2*F.KC+2*RU+RELSUP /* heuristic */)
      {
        long j, k, LIE = (R && lg(W) > 1 && (done_small % 2));
        REL_t *last = cache.last;
        pari_sp av3 = avma;
        GEN p0;
        if (LIE)
        { /* We have full rank for class group and unit. The following tries to
           * improve the prime group lattice by looking for relations involving
           * the primes generating the class group. */
          long n = lg(W)-1; /* need n relations to squash the class group */
          F.L_jid = vecslice(F.perm, 1, n);
          cache.end = cache.last + n;
          /* Lie to the add_rel subsystem: pretend we miss relations involving
           * the primes generating the class group (and only those). */
          cache.missing = n;
          for ( ; n > 0; n--) mael(cache.basis, F.perm[n], F.perm[n]) = 0;
        }
        j = done_small % (F.KC+1);
        if (j == 0) p0 = NULL;
        else
        {
          p0 = gel(F.LP, j);
          if (!A)
          { /* Prevent considering both P_iP_j and P_jP_i in small_norm */
            /* Not all elements end up in F.L_jid (eliminated by hnfspec/add or
             * by trim_list): keep track of which ideals are being considered
             * at each run. */
            long mj = small_multiplier[j];
            for (i = k = 1; i < lg(F.L_jid); i++)
              if (F.L_jid[i] > mj)
              {
                small_multiplier[F.L_jid[i]] = j;
                F.L_jid[k++] = F.L_jid[i];
              }
            setlg(F.L_jid, k);
          }
        }
        if (lg(F.L_jid) > 1)
          small_norm(&cache, &F, nf, Nrelid, M_sn, fact, p0);
        F.L_jid = F.perm; set_avma(av3);
        if (!A && cache.last != last) small_fail = 0; else small_fail++;
        if (LIE)
        { /* restore add_rel subsystem: undo above lie */
          long n = lg(W) - 1;
          for ( ; n > 0; n--) mael(cache.basis, F.perm[n], F.perm[n]) = 1;
          cache.missing = 0;
        }
        cache.end = cache.last;
        done_small++;
        need = F.sfb_chg = 0;
      }
      if (need > 0)
      { /* Random relations */
        if (++nreldep > F.MAXDEPSIZESFB) {
          if (++sfb_trials > SFB_MAX && LIMC < LIMCMAX/2) goto START;
          F.sfb_chg = sfb_INCREASE;
          nreldep = 0;
        }
        else if (!(nreldep % F.MAXDEPSFB))
          F.sfb_chg = sfb_CHANGE;
        if (F.sfb_chg && !subFB_change(&F)) goto START;
        rnd_rel(&cache, &F, nf, fact);
        F.L_jid = F.perm;
      }
      if (DEBUGLEVEL) timer_start(&T);
      if (precpb)
      {
        REL_t *rel;
        if (DEBUGLEVEL)
        {
          char str[64]; sprintf(str,"Buchall_param (%s)",precpb);
          pari_warn(warnprec,str,PREC);
        }
        nf = _nfnewprec(nf, PREC, &nfisclone);
        precdouble++; precpb = NULL;

        if (flag)
        { /* recompute embs only, no need to redo HNF */
          long j, le = lg(embs), lC = lg(C);
          GEN E, M = nf_get_M(nf);
          set_avma(av4);
          for (rel = cache.base+1, i = 1; i < le; i++,rel++)
            gel(embs,i) = rel_embed(rel, &F, embs, i, M, RU, R1, PREC);
          E = RgM_ZM_mul(embs, rowslice(C, RU+1, nbrows(C)));
          for (j = 1; j < lC; j++)
            for (i = 1; i <= RU; i++) gcoeff(C,i,j) = gcoeff(E,i,j);
          av4 = avma;
        }
        else
        { /* recompute embs + HNF */
          for(i = 1; i < lg(PERM); i++) F.perm[i] = PERM[i];
          cache.chk = cache.base;
          W = NULL;
        }
        if (DEBUGLEVEL) timer_printf(&T, "increasing accuracy");
      }
      set_avma(av4);
      if (cache.chk != cache.last)
      { /* Reduce relation matrices */
        long l = cache.last - cache.chk + 1, j;
        GEN mat = cgetg(l, t_MAT);
        REL_t *rel;

        for (j=1,rel = cache.chk + 1; j < l; rel++,j++) gel(mat,j) = rel->R;
        if (!flag || W)
        {
          embs = get_embs(&F, &cache, nf, embs, PREC);
          if (DEBUGLEVEL && timer_get(&T) > 1)
            timer_printf(&T, "floating point embeddings");
        }
        if (!W)
        { /* never reduced before */
          C = flag? matbotid(&cache): embs;
          W = hnfspec_i(mat, F.perm, &dep, &B, &C, F.subFB ? lg(F.subFB)-1:0);
          if (DEBUGLEVEL)
            timer_printf(&T, "hnfspec [%ld x %ld]", lg(F.perm)-1, l-1);
          if (flag)
          {
            PREC += nbits2extraprec(gexpo(C));
            if (nf_get_prec(nf) < PREC) nf = _nfnewprec(nf, PREC, &nfisclone);
            embs = get_embs(&F, &cache, nf, embs, PREC);
            C = vconcat(RgM_ZM_mul(embs, C), C);
          }
          if (DEBUGLEVEL)
            timer_printf(&T, "hnfspec floating points");
        }
        else
        {
          long k = lg(embs);
          GEN E = vecslice(embs, k-l+1,k-1);
          if (flag)
          {
            E = matbotidembs(&cache, E);
            matenlarge(C, cache.last - cache.chk);
          }
          W = hnfadd_i(W, F.perm, &dep, &B, &C, mat, E);
          if (DEBUGLEVEL)
            timer_printf(&T, "hnfadd (%ld + %ld)", l-1, lg(dep)-1);
        }
        gerepileall(av2, 5, &W,&C,&B,&dep,&embs);
        cache.chk = cache.last;
      }
      else if (!W)
      {
        need = old_need;
        F.L_jid = vecslice(F.perm, 1, need);
        continue;
      }
      need = F.KC - (lg(W)-1) - (lg(B)-1);
      if (!need && cache.missing)
      { /* The test above will never be true except if 27449|class number.
         * Ensure that if we have maximal rank for the ideal lattice, then
         * cache.missing == 0. */
        for (i = 1; cache.missing; i++)
          if (!mael(cache.basis, i, i))
          {
            long j;
            cache.missing--; mael(cache.basis, i, i) = 1;
            for (j = i+1; j <= F.KC; j++) mael(cache.basis, j, i) = 0;
          }
      }
      zc = (lg(C)-1) - (lg(B)-1) - (lg(W)-1);
      if (RU-1-zc > 0) need = minss(need + RU-1-zc, F.KC); /* for units */
      if (need)
      { /* dependent rows */
        F.L_jid = vecslice(F.perm, 1, need);
        vecsmall_sort(F.L_jid);
        if (need != old_need) { nreldep = 0; old_need = need; }
      }
      else
      { /* If the relation lattice is too small, check will be > 1 and we will
         * do a new run of small_norm/rnd_rel asking for 1 relation. This often
         * gives a relation involving L_jid[1]. We rotate the first element of
         * L_jid in order to increase the probability of finding relations that
         * increases the lattice. */
        long j, n = lg(W) - 1;
        if (n > 1 && squash_index % n)
        {
          F.L_jid = leafcopy(F.perm);
          for (j = 1; j <= n; j++)
            F.L_jid[j] = F.perm[1 + (j + squash_index - 1) % n];
        }
        else
          F.L_jid = F.perm;
        squash_index++;
      }
    }
    while (need);

    if (!A)
    {
      small_fail = old_need = 0;
      fail_limit = maxss(F.KC / FAIL_DIVISOR, MINFAIL);
    }
    A = vecslice(C, 1, zc); /* cols corresponding to units */
    if (flag) A = rowslice(A, 1, RU);
    Ar = real_i(A);
    R = compute_multiple_of_R(Ar, RU, N, &need, &bit, &lambda);
    if (need < old_need) small_fail = 0;
#if 0 /* A good idea if we are indeed stuck but needs tuning */
    /* we have computed way more relations than should be necessary */
    if (TRIES < 3 && LIMC < LIMCMAX / 8 &&
                     cache.last - cache.base > 10 * F.KC) goto START;
#endif
    old_need = need;
    if (!lambda)
    { precpb = "bestappr"; PREC = myprecdbl(PREC, flag? C: NULL); continue; }
    if (!R)
    { /* not full rank for units */
      if (!need)
      { precpb = "regulator"; PREC = myprecdbl(PREC, flag? C: NULL); }
      continue;
    }
    if (cache.last==old_cache) { need=1; continue; }
    old_cache = cache.last;
    h = ZM_det_triangular(W);
    if (DEBUGLEVEL) err_printf("\n#### Tentative class number: %Ps\n", h);
    i = compute_R(lambda, mulir(h,invhr), &L, &R);
    if (DEBUGLEVEL)
    {
      err_printf("\n");
      timer_printf(&T, "computing regulator and check");
    }
    switch(i)
    {
      case fupb_RELAT:
        need = 1; /* not enough relations */
        continue;
      case fupb_PRECI: /* prec problem unless we cheat on Bach constant */
        if ((precdouble&7) == 7 && LIMC <= LIMCMAX/2) goto START;
        precpb = "compute_R"; PREC = myprecdbl(PREC, flag? C: NULL);
        continue;
    }
    /* DONE */

    if (F.KCZ2 > F.KCZ)
    {
      if (F.sfb_chg && !subFB_change(&F)) goto START;
      if (!be_honest(&F, nf, auts, fact)) goto START;
      if (DEBUGLEVEL) timer_printf(&T, "to be honest");
    }
    F.KCZ2 = 0; /* be honest only once */

    /* fundamental units */
    {
      GEN AU, CU, U, v = extract_full_lattice(L); /* L may be large */
      CU = NULL;
      if (v) { A = vecpermute(A, v); L = vecpermute(L, v); }
      /* arch. components of fund. units */
      U = ZM_lll(L, 0.99, LLL_IM);
      U = ZM_mul(U, lll(RgM_ZM_mul(real_i(A), U)));
      if (DEBUGLEVEL) timer_printf(&T, "units LLL");
      AU = RgM_ZM_mul(A, U);
      A = cleanarchunit(AU, N, NULL, PREC);
      if (!A || lg(A) < RU || expo(gsub(get_regulator(A), R)) > -1)
      {
        long add = nbits2extraprec( gexpo(AU) + 64 ) - gprecision(AU);
        long t = maxss((PREC-2) * 0.15, add);
        if (!A && DEBUGLEVEL) err_printf("### Incorrect units lognorm");
        precpb = "cleanarch"; PREC += maxss(t, EXTRAPREC64); continue;
      }
      if (flag)
      {
        long l = lgcols(C) - RU;
        REL_t *rel;
        SUnits = cgetg(l, t_COL);
        for (rel = cache.base+1, i = 1; i < l; i++,rel++)
          set_rel_alpha(rel, auts, SUnits, i);
        if (RU > 1)
        {
          GEN c = v? vecpermute(C,v): vecslice(C,1,zc);
          CU = ZM_mul(rowslice(c, RU+1, nbrows(c)), U);
        }
      }
      if (DEBUGLEVEL) err_printf("\n#### Computing fundamental units\n");
      fu = getfu(nf, &A, CU? &U: NULL, PREC);
      CU = CU? ZM_mul(CU, U): cgetg(1, t_MAT);
      if (DEBUGLEVEL) timer_printf(&T, "getfu");
      Ce = vecslice(C, zc+1, lg(C)-1);
      if (flag) SUnits = mkvec4(SUnits, CU, rowslice(Ce, RU+1, nbrows(Ce)),
                                utoipos(LIMC));
    }
    /* class group generators */
    if (flag) Ce = rowslice(Ce, 1, RU);
    C0 = Ce; Ce = cleanarch(Ce, N, NULL, PREC);
    if (!Ce) {
      long add = nbits2extraprec( gexpo(C0) + 64 ) - gprecision(C0);
      precpb = "cleanarch"; PREC += maxss(add, 1);
    }
    if (DEBUGLEVEL) timer_printf(&T, "cleanarch");
  } while (need || precpb);

  Vbase = vecpermute(F.LP, F.perm);
  if (!fu) fu = cgetg(1, t_MAT);
  if (!SUnits) SUnits = gen_1;
  clg1 = class_group_gen(nf,W,Ce,Vbase,PREC, &clg2);
  res = mkvec5(clg1, R, SUnits, zu, fu);
  res = buchall_end(nf,res,clg2,W,B,A,Ce,Vbase);
  delete_FB(&F);
  res = gerepilecopy(av0, res);
  if (flag) obj_insert_shallow(res, MATAL, cgetg(1,t_VEC));
  if (nfisclone) gunclone(nf);
  delete_cache(&cache);
  free_GRHcheck(&GRHcheck);
  return res;
}
