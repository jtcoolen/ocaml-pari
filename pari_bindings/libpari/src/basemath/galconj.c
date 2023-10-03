/* Copyright (C) 2000-2003  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_galois

/*************************************************************************/
/**                                                                     **/
/**                           GALOIS CONJUGATES                         **/
/**                                                                     **/
/*************************************************************************/

static int
is2sparse(GEN x)
{
  long i, l = lg(x);
  if (odd(l-3)) return 0;
  for(i=3; i<l; i+=2)
    if (signe(gel(x,i))) return 0;
  return 1;
}

static GEN
galoisconj1(GEN nf)
{
  GEN x = get_nfpol(nf, &nf), f = nf? nf : x, y, z;
  long i, lz, v = varn(x), nbmax;
  pari_sp av = avma;
  RgX_check_ZX(x, "nfgaloisconj");
  nbmax = numberofconjugates(x, 2);
  if (nbmax==1) retmkcol(pol_x(v));
  if (nbmax==2 && is2sparse(x))
  {
    GEN res = cgetg(3,t_COL);
    gel(res,1) = deg1pol_shallow(gen_m1, gen_0, v);
    gel(res,2) = pol_x(v);
    return res;
  }
  x = leafcopy(x); setvarn(x, fetch_var_higher());
  z = nfroots(f, x); lz = lg(z);
  y = cgetg(lz, t_COL);
  for (i = 1; i < lz; i++)
  {
    GEN t = lift(gel(z,i));
    if (typ(t) == t_POL) setvarn(t, v);
    gel(y,i) = t;
  }
  (void)delete_var();
  return gerepileupto(av, y);
}

/*************************************************************************/
/**                                                                     **/
/**                           GALOISCONJ4                               **/
/**                                                                     **/
/*************************************************************************/
/*DEBUGLEVEL:
  1: timing
  2: outline
  4: complete outline
  6: detail
  7: memory
  9: complete detail
*/
struct galois_borne {
  GEN  l;
  long valsol;
  long valabs;
  GEN  bornesol;
  GEN  ladicsol;
  GEN  ladicabs;
  GEN  dis;
};
struct galois_lift {
  GEN  T;
  GEN  den;
  GEN  p;
  GEN  L;
  GEN  Lden;
  long e;
  GEN  Q;
  GEN  TQ;
  struct galois_borne *gb;
};
struct galois_testlift {
  long n;
  long f;
  long g;
  GEN  bezoutcoeff;
  GEN  pauto;
  GEN  C;
  GEN  Cd;
};
struct galois_test { /* data for permutation test */
  GEN order; /* order of tests pour galois_test_perm */
  GEN borne, lborne; /* coefficient bounds */
  GEN ladic;
  GEN PV; /* NULL or vector of test matrices (Vmatrix) */
  GEN TM; /* transpose of M */
  GEN L; /* p-adic roots, known mod ladic */
  GEN M; /* vandermonde inverse */
};
/* result of the study of Frobenius degrees */
enum ga_code {ga_all_normal=1,ga_ext_2=2,ga_non_wss=4,
              ga_all_nilpotent=8,ga_easy=16};
struct galois_analysis {
  long p; /* prime to be lifted */
  long deg; /* degree of the lift */
  long ord;
  long l; /* l: prime number such that T is totally split mod l */
  long p4;
  long group;
};
struct galois_frobenius {
  long p;
  long fp;
  long deg;
  GEN Tmod;
  GEN psi;
};

/* #r = r1 + r2 */
GEN
embed_roots(GEN ro, long r1)
{
  long r2 = lg(ro)-1-r1;
  GEN L;
  if (!r2) L = ro;
  else
  {
    long i,j, N = r1+2*r2;
    L = cgetg(N+1, t_VEC);
    for (i = 1; i <= r1; i++) gel(L,i) = gel(ro,i);
    for (j = i; j <= N; i++)
    {
      GEN z = gel(ro,i);
      gel(L,j++) = z;
      gel(L,j++) = mkcomplex(gel(z,1), gneg(gel(z,2)));
    }
  }
  return L;
}
GEN
embed_disc(GEN z, long r1, long prec)
{
  pari_sp av = avma;
  GEN t = real_1(prec);
  long i, j, n = lg(z)-1, r2 = n-r1;
  for (i = 1; i < r1; i++)
  {
    GEN zi = gel(z,i);
    for (j = i+1; j <= r1; j++) t = gmul(t, gsub(zi, gel(z,j)));
  }
  for (j = r1+1; j <= n; j++)
  {
    GEN zj = gel(z,j), a = gel(zj,1), b = gel(zj,2), b2 = gsqr(b);
    for (i = 1; i <= r1; i++)
    {
      GEN zi = gel(z,i);
      t = gmul(t, gadd(gsqr(gsub(zi, a)), b2));
    }
    t = gmul(t, b);
  }
  if (r2) t = gmul2n(t, r2);
  if (r2 > 1)
  {
    GEN T = real_1(prec);
    for (i = r1+1; i < n; i++)
    {
      GEN zi = gel(z,i), a = gel(zi,1), b = gel(zi,2);
      for (j = i+1; j <= n; j++)
      {
        GEN zj = gel(z,j), c = gel(zj,1), d = gel(zj,2);
        GEN f = gsqr(gsub(a,c)), g = gsqr(gsub(b,d)), h = gsqr(gadd(b,d));
        T = gmul(T, gmul(gadd(f,g), gadd(f,h)));
      }
    }
    t = gmul(t, T);
  }
  t = gsqr(t);
  if (odd(r2)) t = gneg(t);
  return gerepileupto(av, t);
}

/* Compute bound for the coefficients of automorphisms.
 * T a ZX, den a t_INT denominator or NULL */
GEN
initgaloisborne(GEN T, GEN den, long prec, GEN *pL, GEN *pprep, GEN *pD)
{
  GEN L, prep, nf, r;
  pari_timer ti;

  if (DEBUGLEVEL>=4) timer_start(&ti);
  T = get_nfpol(T, &nf);
  r = nf ? nf_get_roots(nf) : NULL;
  if (nf &&  precision(gel(r, 1)) >= prec)
    L = embed_roots(r, nf_get_r1(nf));
  else
    L = QX_complex_roots(T, prec);
  if (DEBUGLEVEL>=4) timer_printf(&ti,"roots");
  prep = vandermondeinverseinit(L);
  if (!den || pD)
  {
    GEN res = RgV_prod(gabs(prep,prec));
    GEN D = ZX_disc_all(T, 1 + gexpo(res)); /* +1 for safety */
    if (pD) *pD = D;
    if (!den) den = indexpartial(T,D);
  }
  if (pprep) *pprep = prep;
  *pL = L; return den;
}

/* ||| M ||| with respect to || x ||_oo, M t_MAT */
GEN
matrixnorm(GEN M, long prec)
{
  long i,j,m, l = lg(M);
  GEN B = real_0(prec);

  if (l == 1) return B;
  m = lgcols(M);
  for (i = 1; i < m; i++)
  {
    GEN z = gabs(gcoeff(M,i,1), prec);
    for (j = 2; j < l; j++) z = gadd(z, gabs(gcoeff(M,i,j), prec));
    if (gcmp(z, B) > 0) B = z;
  }
  return B;
}

static GEN
galoisborne(GEN T, GEN dn, struct galois_borne *gb, long d)
{
  pari_sp ltop, av2;
  GEN borne, borneroots, bornetrace, borneabs;
  long prec;
  GEN L, M, prep, den;
  pari_timer ti;
  const long step=3;

  prec = nbits2prec(bit_accuracy(ZX_max_lg(T)));
  den = initgaloisborne(T,dn,prec, &L,&prep,&gb->dis);
  if (!dn) dn = den;
  ltop = avma;
  if (DEBUGLEVEL>=4) timer_start(&ti);
  M = vandermondeinverse(L, RgX_gtofp(T, prec), den, prep);
  if (DEBUGLEVEL>=4) timer_printf(&ti,"vandermondeinverse");
  borne = matrixnorm(M, prec);
  borneroots = gsupnorm(L, prec); /*t_REAL*/
  bornetrace = mulur((2*step)*degpol(T)/d,
                     powru(borneroots, minss(degpol(T), step)));
  borneroots = ceil_safe(gmul(borne, borneroots));
  borneabs = ceil_safe(gmax_shallow(gmul(borne, bornetrace),
                                    powru(bornetrace, d)));
  av2 = avma;
  /*We use d-1 test, so we must overlift to 2^BITS_IN_LONG*/
  gb->valsol = logint(shifti(borneroots,2+BITS_IN_LONG), gb->l) + 1;
  gb->valabs = logint(shifti(borneabs,2), gb->l) + 1;
  gb->valabs = maxss(gb->valsol, gb->valabs);
  if (DEBUGLEVEL >= 4)
    err_printf("GaloisConj: val1=%ld val2=%ld\n", gb->valsol, gb->valabs);
  set_avma(av2);
  gb->bornesol = gerepileuptoint(ltop, shifti(borneroots,1));
  if (DEBUGLEVEL >= 9)
    err_printf("GaloisConj: Bound %Ps\n",borneroots);
  gb->ladicsol = powiu(gb->l, gb->valsol);
  gb->ladicabs = powiu(gb->l, gb->valabs);
  return dn;
}

static GEN
makeLden(GEN L,GEN den, struct galois_borne *gb)
{ return FpC_Fp_mul(L, den, gb->ladicsol); }

/* Initialize the galois_lift structure */
static void
initlift(GEN T, GEN den, ulong p, GEN L, GEN Lden, struct galois_borne *gb, struct galois_lift *gl)
{
  pari_sp av;
  long e;
  gl->gb = gb;
  gl->T = T;
  gl->den = is_pm1(den)? gen_1: den;
  gl->p = utoipos(p);
  gl->L = L;
  gl->Lden = Lden;
  av = avma;
  e = logint(shifti(gb->bornesol, 2+BITS_IN_LONG), gl->p) + 1;
  set_avma(av);
  if (e < 2) e = 2;
  gl->e = e;
  gl->Q = powuu(p, e);
  gl->TQ = FpX_red(T,gl->Q);
}

/* Check whether f is (with high probability) a solution and compute its
 * permutation */
static int
poltopermtest(GEN f, struct galois_lift *gl, GEN pf)
{
  pari_sp av;
  GEN fx, fp, B = gl->gb->bornesol;
  long i, j, ll;
  for (i = 2; i < lg(f); i++)
    if (abscmpii(gel(f,i),B) > 0)
    {
      if (DEBUGLEVEL>=4) err_printf("GaloisConj: Solution too large.\n");
      if (DEBUGLEVEL>=8) err_printf("f=%Ps\n borne=%Ps\n",f,B);
      return 0;
    }
  ll = lg(gl->L);
  fp = const_vecsmall(ll-1, 1); /* left on stack */
  av = avma;
  for (i = 1; i < ll; i++, set_avma(av))
  {
    fx = FpX_eval(f, gel(gl->L,i), gl->gb->ladicsol);
    for (j = 1; j < ll; j++)
      if (fp[j] && equalii(fx, gel(gl->Lden,j))) { pf[i]=j; fp[j]=0; break; }
    if (j == ll) return 0;
  }
  return 1;
}

static long
galoisfrobeniustest(GEN aut, struct galois_lift *gl, GEN frob)
{
  pari_sp av = avma;
  GEN tlift = aut;
  if (gl->den != gen_1) tlift = FpX_Fp_mul(tlift, gl->den, gl->Q);
  tlift = FpX_center_i(tlift, gl->Q, shifti(gl->Q,-1));
  return gc_long(av, poltopermtest(tlift, gl, frob));
}

static GEN
monoratlift(void *E, GEN S, GEN q)
{
  pari_sp ltop = avma;
  struct galois_lift *gl = (struct galois_lift *) E;
  GEN qm1 = sqrti(shifti(q,-2)), N = gl->Q;
  GEN tlift = FpX_ratlift(S, q, qm1, qm1, gl->den);
  if (tlift)
  {
    pari_sp ltop = avma;
    GEN frob = cgetg(lg(gl->L), t_VECSMALL);
    if(DEBUGLEVEL>=4)
      err_printf("MonomorphismLift: trying early solution %Ps\n",tlift);
    if (gl->den != gen_1)
      tlift = FpX_Fp_mul(FpX_red(Q_muli_to_int(tlift, gl->den), N),
                         Fp_inv(gl->den, N), N);
    if (galoisfrobeniustest(tlift, gl, frob))
    {
      if(DEBUGLEVEL>=4) err_printf("MonomorphismLift: true early solution.\n");
      return gerepilecopy(ltop, tlift);
    }
    if(DEBUGLEVEL>=4) err_printf("MonomorphismLift: false early solution.\n");
  }
  set_avma(ltop);
  return NULL;
}

static GEN
monomorphismratlift(GEN P, GEN S, struct galois_lift *gl)
{
  pari_timer ti;
  if (DEBUGLEVEL >= 1) timer_start(&ti);
  S = ZpX_ZpXQ_liftroot_ea(P,S,gl->T,gl->p, gl->e, (void*)gl, monoratlift);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "monomorphismlift()");
  return S;
}

/* Let T be a polynomial in Z[X] , p a prime number, S in Fp[X]/(T) so
 * that T(S)=0 [p,T]. Lift S in S_0 so that T(S_0)=0 [T,p^e]
 * Unclean stack */
static GEN
automorphismlift(GEN S, struct galois_lift *gl)
{
  return monomorphismratlift(gl->T, S, gl);
}

static GEN
galoisdolift(struct galois_lift *gl)
{
  pari_sp av = avma;
  GEN Tp = FpX_red(gl->T, gl->p);
  GEN S = FpX_Frobenius(Tp, gl->p);
  return gerepileupto(av, automorphismlift(S, gl));
}

static GEN
galoisdoliftn(struct galois_lift *gl, long e)
{
  pari_sp av = avma;
  GEN Tp = FpX_red(gl->T, gl->p);
  GEN S = FpXQ_autpow(FpX_Frobenius(Tp, gl->p), e, Tp, gl->p);
  return gerepileupto(av, automorphismlift(S, gl));
}

static ulong
findpsi(GEN D, ulong pstart, GEN P, GEN S, long o, GEN *Tmod, GEN *Tpsi)
{
  forprime_t iter;
  ulong p;
  long n = degpol(P), i, j, g = n/o;
  GEN psi = cgetg(g+1, t_VECSMALL);
  u_forprime_init(&iter, pstart, ULONG_MAX);
  while ((p = u_forprime_next(&iter)))
  {
    GEN F, Sp;
    long gp = 0;
    if (smodis(D, p) == 0)
      continue;
    F = gel(Flx_factor(ZX_to_Flx(P, p), p), 1);
    if (lg(F)-1 != g) continue;
    Sp = RgX_to_Flx(S, p);
    for (j = 1; j <= g; j++)
    {
      GEN Fj = gel(F, j);
      GEN Sj = Flx_rem(Sp, Fj, p);
      GEN A = Flxq_autpowers(Flx_Frobenius(Fj, p), o,  Fj, p);
      for (i = 1; i <= o; i++)
        if (gequal(Sj, gel(A,i+1)))
        {
          psi[j] = i; break;
        }
      if (i > o) break;
      if (gp==0 && i==1) gp=j;
    }
    if (gp && j > g)
    {
      /* Normalize result so that psi[l]=1 */
      if (gp!=1)
      {
        swap(gel(F,1),gel(F,gp));
        lswap(uel(psi,1),uel(psi,gp));
      }
      *Tpsi = Flv_Fl_div(psi,psi[g],o);
      *Tmod = FlxV_to_ZXV(F);
      return p;
    }
  }
  return 0;
}

static void
inittestlift(GEN plift, GEN Tmod, struct galois_lift *gl,
             struct galois_testlift *gt)
{
  pari_timer ti;
  gt->n = lg(gl->L) - 1;
  gt->g = lg(Tmod) - 1;
  gt->f = gt->n / gt->g;
  gt->bezoutcoeff = bezout_lift_fact(gl->T, Tmod, gl->p, gl->e);
  if (DEBUGLEVEL >= 2) timer_start(&ti);
  gt->pauto = FpXQ_autpowers(plift, gt->f-1, gl->TQ, gl->Q);
  if (DEBUGLEVEL >= 2) timer_printf(&ti, "Frobenius power");
}

/* Explanation of the intheadlong technique:
 * Let C be a bound, B = BITS_IN_LONG, M > C*2^B a modulus and 0 <= a_i < M for
 * i=1,...,n where n < 2^B. We want to test if there exist k,l, |k| < C < M/2^B,
 * such that sum a_i = k + l*M
 * We write a_i*2^B/M = b_i+c_i with b_i integer and 0<=c_i<1, so that
 *   sum b_i - l*2^B = k*2^B/M - sum c_i
 * Since -1 < k*2^B/M < 1 and 0<=c_i<1, it follows that
 *   -n-1 < sum b_i - l*2^B < 1  i.e.  -n <= sum b_i -l*2^B <= 0
 * So we compute z = - sum b_i [mod 2^B] and check if 0 <= z <= n. */

/* Assume 0 <= x < mod. */
static ulong
intheadlong(GEN x, GEN mod)
{
  pari_sp av = avma;
  long res = (long) itou(divii(shifti(x,BITS_IN_LONG),mod));
  return gc_long(av,res);
}
static GEN
vecheadlong(GEN W, GEN mod)
{
  long i, l = lg(W);
  GEN V = cgetg(l, t_VECSMALL);
  for(i=1; i<l; i++) V[i] = intheadlong(gel(W,i), mod);
  return V;
}
static GEN
matheadlong(GEN W, GEN mod)
{
  long i, l = lg(W);
  GEN V = cgetg(l,t_MAT);
  for(i=1; i<l; i++) gel(V,i) = vecheadlong(gel(W,i), mod);
  return V;
}
static ulong
polheadlong(GEN P, long n, GEN mod)
{
  return (lg(P)>n+2)? intheadlong(gel(P,n+2),mod): 0;
}

#define headlongisint(Z,N) (-(ulong)(Z)<=(ulong)(N))

static long
frobeniusliftall(GEN sg, long el, GEN *psi, struct galois_lift *gl,
                 struct galois_testlift *gt, GEN frob)
{
  pari_sp av, ltop2, ltop = avma;
  long i,j,k, c = lg(sg)-1, n = lg(gl->L)-1, m = gt->g, d = m / c;
  GEN pf, u, v, C, Cd, SG, cache;
  long N1, N2, R1, Ni, ord = gt->f, c_idx = gt->g-1;
  ulong headcache;
  long hop = 0;
  GEN NN, NQ;
  pari_timer ti;

  *psi = pf = cgetg(m, t_VECSMALL);
  ltop2 = avma;
  NN = diviiexact(mpfact(m), mului(c, powiu(mpfact(d), c)));
  if (DEBUGLEVEL >= 4)
    err_printf("GaloisConj: I will try %Ps permutations\n", NN);
  N1=10000000;
  NQ=divis_rem(NN,N1,&R1);
  if (abscmpiu(NQ,1000000000)>0)
  {
    pari_warn(warner,"Combinatorics too hard : would need %Ps tests!\n"
        "I will skip it, but it may induce an infinite loop",NN);
    *psi = NULL; return gc_long(ltop,0);
  }
  N2=itos(NQ); if(!N2) N1=R1;
  if (DEBUGLEVEL>=4) timer_start(&ti);
  set_avma(ltop2);
  C = gt->C;
  Cd= gt->Cd;
  v = FpXQ_mul(gel(gt->pauto, 1+el%ord), gel(gt->bezoutcoeff, m),gl->TQ,gl->Q);
  if (gl->den != gen_1) v = FpX_Fp_mul(v, gl->den, gl->Q);
  SG = cgetg(lg(sg),t_VECSMALL);
  for(i=1; i<lg(SG); i++) SG[i] = (el*sg[i])%ord + 1;
  cache = cgetg(m+1,t_VECSMALL); cache[m] = polheadlong(v,1,gl->Q);
  headcache = polheadlong(v,2,gl->Q);
  for (i = 1; i < m; i++) pf[i] = 1 + i/d;
  av = avma;
  for (Ni = 0, i = 0; ;i++)
  {
    for (j = c_idx ; j > 0; j--)
    {
      long h = SG[pf[j]];
      if (!mael(C,h,j))
      {
        pari_sp av3 = avma;
        GEN r = FpXQ_mul(gel(gt->pauto,h), gel(gt->bezoutcoeff,j),gl->TQ,gl->Q);
        if (gl->den != gen_1) r = FpX_Fp_mul(r, gl->den, gl->Q);
        gmael(C,h,j) = gclone(r);
        mael(Cd,h,j) = polheadlong(r,1,gl->Q);
        set_avma(av3);
      }
      uel(cache,j) = uel(cache,j+1)+umael(Cd,h,j);
    }
    if (headlongisint(uel(cache,1),n))
    {
      ulong head = headcache;
      for (j = 1; j < m; j++) head += polheadlong(gmael(C,SG[pf[j]],j),2,gl->Q);
      if (headlongisint(head,n))
      {
        u = v;
        for (j = 1; j < m; j++) u = ZX_add(u, gmael(C,SG[pf[j]],j));
        u = FpX_center_i(FpX_red(u, gl->Q), gl->Q, shifti(gl->Q,-1));
        if (poltopermtest(u, gl, frob))
        {
          if (DEBUGLEVEL >= 4)
          {
            timer_printf(&ti, "");
            err_printf("GaloisConj: %d hops on %Ps tests\n",hop,addis(mulss(Ni,N1),i));
          }
          return gc_long(ltop2,1);
        }
        if (DEBUGLEVEL >= 4) err_printf("M");
      }
      else hop++;
    }
    if (DEBUGLEVEL >= 4 && i % maxss(N1/20, 1) == 0)
      timer_printf(&ti, "GaloisConj:Testing %Ps", addis(mulss(Ni,N1),i));
    set_avma(av);
    if (i == N1-1)
    {
      if (Ni==N2-1) N1 = R1;
      if (Ni==N2) break;
      Ni++; i = 0;
      if (DEBUGLEVEL>=4) timer_start(&ti);
    }
    for (j = 2; j < m && pf[j-1] >= pf[j]; j++)
      /*empty*/; /* to kill clang Warning */
    for (k = 1; k < j-k && pf[k] != pf[j-k]; k++) { lswap(pf[k], pf[j-k]); }
    for (k = j - 1; pf[k] >= pf[j]; k--)
      /*empty*/;
    lswap(pf[j], pf[k]); c_idx = j;
  }
  if (DEBUGLEVEL>=4) err_printf("GaloisConj: not found, %d hops \n",hop);
  *psi = NULL; return gc_long(ltop,0);
}

/* Compute the test matrix for the i-th line of V. Clone. */
static GEN
Vmatrix(long i, struct galois_test *td)
{
  pari_sp av = avma;
  GEN m = gclone( matheadlong(FpC_FpV_mul(td->L, row(td->M,i), td->ladic), td->ladic));
  set_avma(av); return m;
}

/* Initialize galois_test */
static void
inittest(GEN L, GEN M, GEN borne, GEN ladic, struct galois_test *td)
{
  long i, n = lg(L)-1;
  GEN p = cgetg(n+1, t_VECSMALL);
  if (DEBUGLEVEL >= 8) err_printf("GaloisConj: Init Test\n");
  td->order = p;
  for (i = 1; i <= n-2; i++) p[i] = i+2;
  p[n-1] = 1; p[n] = 2;
  td->borne = borne;
  td->lborne = subii(ladic, borne);
  td->ladic = ladic;
  td->L = L;
  td->M = M;
  td->TM = shallowtrans(M);
  td->PV = zero_zv(n);
  gel(td->PV, 2) = Vmatrix(2, td);
}

/* Free clones stored inside galois_test */
static void
freetest(struct galois_test *td)
{
  long i;
  for (i = 1; i < lg(td->PV); i++)
    if (td->PV[i]) { gunclone(gel(td->PV,i)); td->PV[i] = 0; }
}

/* Check if the integer P seen as a p-adic number is close to an integer less
 * than td->borne in absolute value */
static long
padicisint(GEN P, struct galois_test *td)
{
  pari_sp ltop = avma;
  GEN U  = modii(P, td->ladic);
  long r = cmpii(U, td->borne) <= 0 || cmpii(U, td->lborne) >= 0;
  return gc_long(ltop, r);
}

/* Check if the permutation pf is valid according to td.
 * If not, update td to make subsequent test faster (hopefully) */
static long
galois_test_perm(struct galois_test *td, GEN pf)
{
  pari_sp av = avma;
  long i, j, n = lg(td->L)-1;
  GEN V, P = NULL;
  for (i = 1; i < n; i++)
  {
    long ord = td->order[i];
    GEN PW = gel(td->PV, ord);
    if (PW)
    {
      ulong head = umael(PW,1,pf[1]);
      for (j = 2; j <= n; j++) head += umael(PW,j,pf[j]);
      if (!headlongisint(head,n)) break;
    } else
    {
      if (!P) P = vecpermute(td->L, pf);
      V = FpV_dotproduct(gel(td->TM,ord), P, td->ladic);
      if (!padicisint(V, td)) {
        gel(td->PV, ord) = Vmatrix(ord, td);
        if (DEBUGLEVEL >= 4) err_printf("M");
        break;
      }
    }
  }
  if (i == n) return gc_long(av,1);
  if (DEBUGLEVEL >= 4) err_printf("%d.", i);
  if (i > 1)
  {
    long z = td->order[i];
    for (j = i; j > 1; j--) td->order[j] = td->order[j-1];
    td->order[1] = z;
    if (DEBUGLEVEL >= 8) err_printf("%Ps", td->order);
  }
  return gc_long(av,0);
}
/*Compute a*b/c when a*b will overflow*/
static long
muldiv(long a,long b,long c)
{
  return (long)((double)a*(double)b/c);
}

/* F = cycle decomposition of sigma,
 * B = cycle decomposition of cl(tau).
 * Check all permutations pf who can possibly correspond to tau, such that
 * tau*sigma*tau^-1 = sigma^s and tau^d = sigma^t, where d = ord cl(tau)
 * x: vector of choices,
 * G: vector allowing linear access to elts of F.
 * Choices multiple of e are not changed. */
static GEN
testpermutation(GEN F, GEN B, GEN x, long s, long e, long cut,
                struct galois_test *td)
{
  pari_sp av, avm = avma;
  long a, b, c, d, n, p1, p2, p3, p4, p5, p6, l1, l2, N1, N2, R1;
  long i, j, cx, hop = 0, start = 0;
  GEN pf, ar, G, W, NN, NQ;
  pari_timer ti;
  if (DEBUGLEVEL >= 1) timer_start(&ti);
  a = lg(F)-1; b = lg(gel(F,1))-1;
  c = lg(B)-1; d = lg(gel(B,1))-1;
  n = a*b;
  s = (b+s) % b;
  pf = cgetg(n+1, t_VECSMALL);
  av = avma;
  ar = cgetg(a+2, t_VECSMALL); ar[a+1]=0;
  G  = cgetg(a+1, t_VECSMALL);
  W  = gel(td->PV, td->order[n]);
  for (cx=1, i=1, j=1; cx <= a; cx++, i++)
  {
    gel(G,cx) = gel(F, coeff(B,i,j));
    if (i == d) { i = 0; j++; }
  }
  NN = divis(powuu(b, c * (d - d/e)), cut);
  if (DEBUGLEVEL>=4) err_printf("GaloisConj: I will try %Ps permutations\n", NN);
  N1 = 1000000;
  NQ = divis_rem(NN,N1,&R1);
  if (abscmpiu(NQ,100000000)>0)
  {
    set_avma(avm);
    pari_warn(warner,"Combinatorics too hard: would need %Ps tests!\n"
                     "I'll skip it but you will get a partial result...",NN);
    return identity_perm(n);
  }
  N2 = itos(NQ);
  for (l2 = 0; l2 <= N2; l2++)
  {
    long nbiter = (l2<N2) ? N1: R1;
    if (DEBUGLEVEL >= 2 && N2) err_printf("%d%% ", muldiv(l2,100,N2));
    for (l1 = 0; l1 < nbiter; l1++)
    {
      if (start)
      {
        for (i=1, j=e; i < a;)
        {
          if ((++(x[i])) != b) break;
          x[i++] = 0;
          if (i == j) { i++; j += e; }
        }
      }
      else { start=1; i = a-1; }
      /* intheadlong test: overflow in + is OK, we compute mod 2^BIL */
      for (p1 = i+1, p5 = p1%d - 1 ; p1 >= 1; p1--, p5--) /* p5 = (p1%d) - 1 */
      {
        GEN G1, G6;
        ulong V = 0;
        if (p5 == - 1) { p5 = d - 1; p6 = p1 + 1 - d; } else p6 = p1 + 1;
        G1 = gel(G,p1); G6 = gel(G,p6);
        p4 = p5 ? x[p1-1] : 0;
        for (p2 = 1+p4, p3 = 1 + x[p1]; p2 <= b; p2++)
        {
          V += umael(W,uel(G6,p3),uel(G1,p2));
          p3 += s; if (p3 > b) p3 -= b;
        }
        p3 = 1 + x[p1] - s; if (p3 <= 0) p3 += b;
        for (p2 = p4; p2 >= 1; p2--)
        {
          V += umael(W,uel(G6,p3),uel(G1,p2));
          p3 -= s; if (p3 <= 0) p3 += b;
        }
        uel(ar,p1) = uel(ar,p1+1) + V;
      }
      if (!headlongisint(uel(ar,1),n)) continue;

      /* intheadlong succeeds. Full computation */
      for (p1=1, p5=d; p1 <= a; p1++, p5++)
      {
        if (p5 == d) { p5 = 0; p4 = 0; } else p4 = x[p1-1];
        if (p5 == d-1) p6 = p1+1-d; else p6 = p1+1;
        for (p2 = 1+p4, p3 = 1 + x[p1]; p2 <= b; p2++)
        {
          pf[mael(G,p1,p2)] = mael(G,p6,p3);
          p3 += s; if (p3 > b) p3 -= b;
        }
        p3 = 1 + x[p1] - s; if (p3 <= 0) p3 += b;
        for (p2 = p4; p2 >= 1; p2--)
        {
          pf[mael(G,p1,p2)] = mael(G,p6,p3);
          p3 -= s; if (p3 <= 0) p3 += b;
        }
      }
      if (galois_test_perm(td, pf))
      {
        if (DEBUGLEVEL >= 1)
        {
          GEN nb = addis(mulss(l2,N1),l1);
          timer_printf(&ti, "testpermutation(%Ps)", nb);
          if (DEBUGLEVEL >= 2 && hop)
            err_printf("GaloisConj: %d hop over %Ps iterations\n", hop, nb);
        }
        set_avma(av); return pf;
      }
      hop++;
    }
  }
  if (DEBUGLEVEL >= 1)
  {
    timer_printf(&ti, "testpermutation(%Ps)", NN);
    if (DEBUGLEVEL >= 2 && hop)
      err_printf("GaloisConj: %d hop over %Ps iterations\n", hop, NN);
  }
  return gc_NULL(avm);
}

/* List of subgroups of (Z/mZ)^* whose order divide o, and return the list
 * of their elements, sorted by increasing order */
static GEN
listznstarelts(long m, long o)
{
  pari_sp av = avma;
  GEN L, zn, zns;
  long i, phi, ind, l;
  if (m == 2) retmkvec(mkvecsmall(1));
  zn = znstar(stoi(m));
  phi = itos(gel(zn,1));
  o = ugcd(o, phi); /* do we impose this on input ? */
  zns = znstar_small(zn);
  L = cgetg(o+1, t_VEC);
  for (i=1,ind = phi; ind; ind -= phi/o, i++) /* by *decreasing* exact index */
    gel(L,i) = subgrouplist(gel(zn,2), mkvec(utoipos(ind)));
  L = shallowconcat1(L); l = lg(L);
  for (i = 1; i < l; i++) gel(L,i) = znstar_hnf_elts(zns, gel(L,i));
  return gerepilecopy(av, L);
}

/* A sympol is a symmetric polynomial
 *
 * Currently sympol are couple of t_VECSMALL [v,w]
 * v[1]...v[k], w[1]...w[k]  represent the polynomial sum(i=1,k,v[i]*s_w[i])
 * where s_i(X_1,...,X_n) = sum(j=1,n,X_j^i) */

static GEN
Flm_newtonsum(GEN M, ulong e, ulong p)
{
  long f = lg(M), g = lg(gel(M,1)), i, j;
  GEN NS = cgetg(f, t_VECSMALL);
  for(i=1; i<f; i++)
  {
    ulong s = 0;
    GEN Mi = gel(M,i);
    for(j = 1; j < g; j++)
      s = Fl_add(s, Fl_powu(uel(Mi,j), e, p), p);
    uel(NS,i) = s;
  }
  return NS;
}

static GEN
Flv_sympol_eval(GEN v, GEN NS, ulong p)
{
  pari_sp av = avma;
  long i, l = lg(v);
  GEN S = Flv_Fl_mul(gel(NS,1), uel(v,1), p);
  for (i=2; i<l; i++)
    if (v[i]) S = Flv_add(S, Flv_Fl_mul(gel(NS,i), uel(v,i), p), p);
  return gerepileuptoleaf(av, S);
}

static GEN
sympol_eval_newtonsum(long e, GEN O, GEN mod)
{
  long f = lg(O), g = lg(gel(O,1)), i, j;
  GEN PL = cgetg(f, t_COL);
  for(i=1; i<f; i++)
  {
    pari_sp av = avma;
    GEN s = gen_0;
    for(j=1; j<g; j++) s = addii(s, Fp_powu(gmael(O,i,j), e, mod));
    gel(PL,i) = gerepileuptoint(av, remii(s,mod));
  }
  return PL;
}

static GEN
sympol_eval(GEN sym, GEN O, GEN mod)
{
  pari_sp av = avma;
  long i;
  GEN v = gel(sym,1), w = gel(sym,2);
  GEN S = gen_0;
  for (i=1; i<lg(v); i++)
    if (v[i]) S = gadd(S, gmulsg(v[i],  sympol_eval_newtonsum(w[i], O, mod)));
  return gerepileupto(av, S);
}

/* Let sigma be an automorphism of L (as a polynomial with rational coefs)
 * Let 'sym' be a symmetric polynomial defining alpha in L.
 * We have alpha = sym(x,sigma(x),,,sigma^(g-1)(x)). Compute alpha mod p */
static GEN
sympol_aut_evalmod(GEN sym, long g, GEN sigma, GEN Tp, GEN p)
{
  pari_sp ltop=avma;
  long i, j, npows = brent_kung_optpow(degpol(Tp)-1, g-1, 1);
  GEN s, f, pows, v = zv_to_ZV(gel(sym,1)), w = zv_to_ZV(gel(sym,2));
  sigma = RgX_to_FpX(sigma, p);
  pows  = FpXQ_powers(sigma,npows,Tp,p);
  f = pol_x(varn(sigma));
  s = pol_0(varn(sigma));
  for(i=1; i<=g;i++)
  {
    if (i > 1) f = FpX_FpXQV_eval(f,pows,Tp,p);
    for(j=1; j<lg(v); j++)
      s = FpX_add(s, FpX_Fp_mul(FpXQ_pow(f,gel(w,j),Tp,p),gel(v,j),p),p);
  }
  return gerepileupto(ltop, s);
}

/* Let Sp be as computed with sympol_aut_evalmod
 * Let Tmod be the factorisation of T mod p.
 * Return the factorisation of the minimal polynomial of S mod p w.r.t. Tmod */
static GEN
fixedfieldfactmod(GEN Sp, GEN p, GEN Tmod)
{
  long i, l = lg(Tmod);
  GEN F = cgetg(l,t_VEC);
  for(i=1; i<l; i++)
  {
    GEN Ti = gel(Tmod,i);
    gel(F,i) = FpXQ_minpoly(FpX_rem(Sp,Ti,p), Ti,p);
  }
  return F;
}

static GEN
fixedfieldsurmer(ulong l, GEN NS, GEN W)
{
  const long step=3;
  long i, j, n = lg(W)-1, m = 1L<<((n-1)<<1);
  GEN sym = cgetg(n+1,t_VECSMALL);
  for (j=1;j<n;j++) sym[j] = step;
  sym[n] = 0;
  if (DEBUGLEVEL>=4) err_printf("FixedField: Weight: %Ps\n",W);
  for (i=0; i<m; i++)
  {
    pari_sp av = avma;
    GEN L;
    for (j=1; sym[j]==step; j++) sym[j]=0;
    sym[j]++;
    if (DEBUGLEVEL>=6) err_printf("FixedField: Sym: %Ps\n",sym);
    L = Flv_sympol_eval(sym, NS, l);
    if (!vecsmall_is1to1(L)) { set_avma(av); continue; }
    return mkvec2(sym,W);
  }
  return NULL;
}

/*Check whether the line of NS are pair-wise distinct.*/
static long
sympol_is1to1_lg(GEN NS, long n)
{
  long i, j, k, l = lgcols(NS);
  for (i=1; i<l; i++)
    for(j=i+1; j<l; j++)
    {
      for(k=1; k<n; k++)
        if (mael(NS,k,j)!=mael(NS,k,i)) break;
      if (k>=n) return 0;
    }
  return 1;
}

/* Let O a set of orbits of roots (see fixedfieldorbits) modulo mod,
 * l | mod and p two prime numbers. Return a vector [sym,s,P] where:
 * sym is a sympol, s is the set of images of sym on O and
 * P is the polynomial with roots s. */
static GEN
fixedfieldsympol(GEN O, ulong l)
{
  pari_sp ltop=avma;
  const long n=(BITS_IN_LONG>>1)-1;
  GEN NS = cgetg(n+1,t_MAT), sym = NULL, W = cgetg(n+1,t_VECSMALL);
  long i, e=1;
  if (DEBUGLEVEL>=4)
    err_printf("FixedField: Size: %ldx%ld\n",lg(O)-1,lg(gel(O,1))-1);
  O = ZM_to_Flm(O,l);
  for (i=1; !sym && i<=n; i++)
  {
    GEN L = Flm_newtonsum(O, e++, l);
    if (lg(O)>2)
      while (vecsmall_isconst(L)) L = Flm_newtonsum(O, e++, l);
    W[i] = e-1; gel(NS,i) = L;
    if (sympol_is1to1_lg(NS,i+1))
      sym = fixedfieldsurmer(l,NS,vecsmall_shorten(W,i));
  }
  if (!sym) pari_err_BUG("fixedfieldsympol [p too small]");
  if (DEBUGLEVEL>=2) err_printf("FixedField: Found: %Ps\n",gel(sym,1));
  return gerepilecopy(ltop,sym);
}

/* Let O a set of orbits as indices and L the corresponding roots.
 * Return the set of orbits as roots. */
static GEN
fixedfieldorbits(GEN O, GEN L)
{
  GEN S = cgetg(lg(O), t_MAT);
  long i;
  for (i = 1; i < lg(O); i++) gel(S,i) = vecpermute(L, gel(O,i));
  return S;
}

static GEN
fixedfieldinclusion(GEN O, GEN PL)
{
  long i, j, f = lg(O)-1, g = lg(gel(O,1))-1;
  GEN S = cgetg(f*g + 1, t_COL);
  for (i = 1; i <= f; i++)
  {
    GEN Oi = gel(O,i);
    for (j = 1; j <= g; j++) gel(S, Oi[j]) = gel(PL, i);
  }
  return S;
}

/* Polynomial attached to a vector of conjugates. Not stack clean */
static GEN
vectopol(GEN v, GEN M, GEN den , GEN mod, GEN mod2, long x)
{
  long l = lg(v)+1, i;
  GEN z = cgetg(l,t_POL);
  z[1] = evalsigne(1)|evalvarn(x);
  for (i=2; i<l; i++)
    gel(z,i) = gdiv(centermodii(ZMrow_ZC_mul(M,v,i-1), mod, mod2), den);
  return normalizepol_lg(z, l);
}

/* Polynomial associate to a permutation of the roots. Not stack clean */
static GEN
permtopol(GEN p, GEN L, GEN M, GEN den, GEN mod, GEN mod2, long x)
{
  if (lg(p) != lg(L)) pari_err_TYPE("permtopol [permutation]", p);
  return vectopol(vecpermute(L,p), M, den, mod, mod2, x);
}

static GEN
galoisvecpermtopol(GEN gal, GEN vec, GEN mod, GEN mod2)
{
  long i, l = lg(vec);
  long v = varn(gal_get_pol(gal));
  GEN L = gal_get_roots(gal);
  GEN M = gal_get_invvdm(gal);
  GEN P = cgetg(l, t_MAT);
  for (i=1; i<l; i++)
    gel(P, i) = vecpermute(L,gel(vec,i));
  P = RgM_to_RgXV(FpM_center(FpM_mul(M, P, mod), mod, mod2), v);
  return gdiv(P, gal_get_den(gal));
}

static void
notgalois(long p, struct galois_analysis *ga)
{
  if (DEBUGLEVEL >= 2) err_printf("GaloisAnalysis:non Galois for p=%ld\n", p);
  ga->p = p;
  ga->deg = 0;
}

/*Gather information about the group*/
static long
init_group(long n, long np, GEN Fp, GEN Fe, long *porder)
{
  const long prim_nonwss_orders[] = { 48,56,60,72,75,80,196,200,216 };
  long i, phi_order = 1, order = 1, group = 0;
  ulong p;

 /* non-WSS groups of this order? */
  for (i=0; i < (long)numberof(prim_nonwss_orders); i++)
    if (n % prim_nonwss_orders[i] == 0) { group |= ga_non_wss; break; }
  if (np == 2 && Fp[2] == 3 && Fe[2] == 1 && Fe[1] > 2) group |= ga_ext_2;

  for (i = np; i > 0; i--)
  {
    long p = Fp[i];
    if (phi_order % p == 0) { group |= ga_all_normal; break; }
    order *= p; phi_order *= p-1;
    if (Fe[i] > 1) break;
  }
  if (uisprimepower(n, &p) || n == 135) group |= ga_all_nilpotent;
  if (n <= 104) group |= ga_easy; /* no need to use polynomial algo */
  *porder = order; return group;
}

/*is a "better" than b ? (if so, update karma) */
static int
improves(long a, long b, long plift, long p, long n, long *karma)
{
  if (!plift || a > b) { *karma = ugcd(p-1,n); return 1; }
  if (a == b) {
    long k = ugcd(p-1,n);
    if (k > *karma) { *karma = k; return 1; }
  }
  return 0; /* worse */
}

/* return 0 if not galois or not wss */
static int
galoisanalysis(GEN T, struct galois_analysis *ga, long calcul_l, GEN bad)
{
  pari_sp ltop = avma, av;
  long group, linf, n, p, i, karma = 0;
  GEN F, Fp, Fe, Fpe, O;
  long np, order, plift, nbmax, nbtest, deg;
  pari_timer ti;
  forprime_t S;
  if (DEBUGLEVEL >= 1) timer_start(&ti);
  n = degpol(T);
  O = zero_zv(n);
  F = factoru_pow(n);
  Fp = gel(F,1); np = lg(Fp)-1;
  Fe = gel(F,2);
  Fpe= gel(F,3);
  group = init_group(n, np, Fp, Fe, &order);

  /*Now we study the orders of the Frobenius elements*/
  deg = Fp[np]; /* largest prime | n */
  plift = 0;
  nbtest = 0;
  nbmax = 8+(n>>1);
  u_forprime_init(&S, n*maxss(expu(n)-3, 2), ULONG_MAX);
  av = avma;
  while (!plift || (nbtest < nbmax && (nbtest <=8 || order < (n>>1)))
                || ((n == 24 || n==36) && O[6] == 0 && O[4] == 0)
                || ((group&ga_non_wss) && order == Fp[np]))
  {
    long d, o, norm_o = 1;
    GEN D, Tp;

    if ((group&ga_non_wss) && nbtest >= 3*nbmax) break; /* in all cases */
    nbtest++; set_avma(av);
    p = u_forprime_next(&S);
    if (!p) pari_err_OVERFLOW("galoisanalysis [ran out of primes]");
    if (bad && dvdiu(bad, p)) continue;
    Tp = ZX_to_Flx(T,p);
    if (!Flx_is_squarefree(Tp,p)) { if (!--nbtest) nbtest = 1; continue; }

    D = Flx_nbfact_by_degree(Tp, &d, p);
    o = n / d; /* d factors, all should have degree o */
    if (D[o] != d) { notgalois(p, ga); return gc_bool(ltop,0); }

    if (!O[o]) O[o] = p;
    if (o % deg) goto ga_end; /* NB: deg > 1 */
    if ((group&ga_all_normal) && o < order) goto ga_end;

    /*Frob_p has order o > 1, find a power which generates a normal subgroup*/
    if (o * Fp[1] >= n)
      norm_o = o; /*subgroups of smallest index are normal*/
    else
    {
      for (i = np; i > 0; i--)
      {
        if (o % Fpe[i]) break;
        norm_o *= Fpe[i];
      }
    }
    /* Frob_p^(o/norm_o) generates a normal subgroup of order norm_o */
    if (norm_o != 1)
    {
      if (!(group&ga_all_normal) || o > order)
        karma = ugcd(p-1,n);
      else if (!improves(norm_o, deg, plift,p,n, &karma)) goto ga_end;
      /* karma0=0, deg0<=norm_o -> the first improves() returns 1 */
      deg = norm_o; group |= ga_all_normal; /* STORE */
    }
    else if (group&ga_all_normal) goto ga_end;
    else if (!improves(o, order, plift,p,n, &karma)) goto ga_end;

    order = o; plift = p; /* STORE */
    ga_end:
    if (DEBUGLEVEL >= 5)
      err_printf("GaloisAnalysis:Nbtest=%ld,p=%ld,o=%ld,n_o=%d,best p=%ld,ord=%ld,k=%ld\n", nbtest, p, o, norm_o, plift, order,karma);
  }
  /* To avoid looping on non-WSS group.
   * TODO: check for large groups. Would it be better to disable this check if
   * we are in a good case (ga_all_normal && !(ga_ext_2) (e.g. 60)) ?*/
  ga->p = plift;
  if (!plift || ((group&ga_non_wss) && order == Fp[np]))
  {
    if (DEBUGLEVEL)
      pari_warn(warner,"Galois group probably not weakly super solvable");
    return 0;
  }
  linf = 2*n*usqrt(n);
  if (calcul_l && O[1] <= linf)
  {
    pari_sp av2;
    forprime_t S2;
    ulong p;
    u_forprime_init(&S2, linf+1,ULONG_MAX);
    av2 = avma;
    while ((p = u_forprime_next(&S2)))
    { /*find a totally split prime l > linf*/
      GEN Tp = ZX_to_Flx(T, p);
      long nb = Flx_nbroots(Tp, p);
      if (nb == n) { O[1] = p; break; }
      if (nb && Flx_is_squarefree(Tp,p)) { notgalois(p,ga); return gc_bool(ltop,0); }
      set_avma(av2);
    }
    if (!p) pari_err_OVERFLOW("galoisanalysis [ran out of primes]");
  }
  ga->group = group;
  ga->deg = deg;
  ga->ord = order;
  ga->l  = O[1];
  ga->p4 = n >= 4 ? O[4] : 0;
  if (DEBUGLEVEL >= 4)
    err_printf("GaloisAnalysis:p=%ld l=%ld group=%ld deg=%ld ord=%ld\n",
               plift, O[1], group, deg, order);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "galoisanalysis()");
  return gc_bool(ltop,1);
}

static GEN
a4galoisgen(struct galois_test *td)
{
  const long n = 12;
  pari_sp ltop = avma, av, av2;
  long i, j, k, N, hop = 0;
  GEN MT, O,O1,O2,O3, ar, mt, t, u, res, orb, pft, pfu, pfv;

  res = cgetg(3, t_VEC);
  pft = cgetg(n+1, t_VECSMALL);
  pfu = cgetg(n+1, t_VECSMALL);
  pfv = cgetg(n+1, t_VECSMALL);
  gel(res,1) = mkvec3(pft,pfu,pfv);
  gel(res,2) = mkvecsmall3(2,2,3);
  av = avma;
  ar = cgetg(5, t_VECSMALL);
  mt = gel(td->PV, td->order[n]);
  t = identity_perm(n) + 1; /* Sorry for this hack */
  u = cgetg(n+1, t_VECSMALL) + 1; /* too lazy to correct */
  MT = cgetg(n+1, t_MAT);
  for (j = 1; j <= n; j++) gel(MT,j) = cgetg(n+1, t_VECSMALL);
  for (j = 1; j <= n; j++)
    for (i = 1; i < j; i++)
      ucoeff(MT,i,j) = ucoeff(MT,j,i) = ucoeff(mt,i,j)+ucoeff(mt,j,i);
  /* MT(i,i) unused */

  av2 = avma;
  /* N = itos(gdiv(mpfact(n), mpfact(n >> 1))) >> (n >> 1); */
  /* n = 2k = 12; N = (2k)! / (k! * 2^k) = 10395 */
  N = 10395;
  if (DEBUGLEVEL>=4) err_printf("A4GaloisConj: will test %ld permutations\n", N);
  uel(ar,4) = umael(MT,11,12);
  uel(ar,3) = uel(ar,4) + umael(MT,9,10);
  uel(ar,2) = uel(ar,3) + umael(MT,7,8);
  uel(ar,1) = uel(ar,2) + umael(MT,5,6);
  for (i = 0; i < N; i++)
  {
    long g;
    if (i)
    {
      long a, x = i, y = 1;
      do { y += 2; a = x%y; x = x/y; } while (!a);
      switch (y)
      {
      case 3:
        lswap(t[2], t[2-a]);
        break;
      case 5:
        x = t[0]; t[0] = t[2]; t[2] = t[1]; t[1] = x;
        lswap(t[4], t[4-a]);
        uel(ar,1) = uel(ar,2) + umael(MT,t[4],t[5]);
        break;
      case 7:
        x = t[0]; t[0] = t[4]; t[4] = t[3]; t[3] = t[1]; t[1] = t[2]; t[2] = x;
        lswap(t[6], t[6-a]);
        uel(ar,2) = uel(ar,3) + umael(MT,t[6],t[7]);
        uel(ar,1) = uel(ar,2) + umael(MT,t[4],t[5]);
        break;
      case 9:
        x = t[0]; t[0] = t[6]; t[6] = t[5]; t[5] = t[3]; t[3] = x;
        lswap(t[1], t[4]);
        lswap(t[8], t[8-a]);
        uel(ar,3) = uel(ar,4) + umael(MT,t[8],t[9]);
        uel(ar,2) = uel(ar,3) + umael(MT,t[6],t[7]);
        uel(ar,1) = uel(ar,2) + umael(MT,t[4],t[5]);
        break;
      case 11:
        x = t[0]; t[0] = t[8]; t[8] = t[7]; t[7] = t[5]; t[5] = t[1];
        t[1] = t[6]; t[6] = t[3]; t[3] = t[2]; t[2] = t[4]; t[4] = x;
        lswap(t[10], t[10-a]);
        uel(ar,4) = umael(MT,t[10],t[11]);
        uel(ar,3) = uel(ar,4) + umael(MT,t[8],t[9]);
        uel(ar,2) = uel(ar,3) + umael(MT,t[6],t[7]);
        uel(ar,1) = uel(ar,2) + umael(MT,t[4],t[5]);
      }
    }
    g = uel(ar,1)+umael(MT,t[0],t[1])+umael(MT,t[2],t[3]);
    if (headlongisint(g,n))
    {
      for (k = 0; k < n; k += 2)
      {
        pft[t[k]] = t[k+1];
        pft[t[k+1]] = t[k];
      }
      if (galois_test_perm(td, pft)) break;
      hop++;
    }
    set_avma(av2);
  }
  if (DEBUGLEVEL >= 1 && hop)
    err_printf("A4GaloisConj: %ld hop over %ld iterations\n", hop, N);
  if (i == N) return gc_NULL(ltop);
  /* N = itos(gdiv(mpfact(n >> 1), mpfact(n >> 2))) >> 1; */
  N = 60;
  if (DEBUGLEVEL >= 4) err_printf("A4GaloisConj: sigma=%Ps \n", pft);
  for (k = 0; k < n; k += 4)
  {
    u[k+3] = t[k+3];
    u[k+2] = t[k+1];
    u[k+1] = t[k+2];
    u[k]   = t[k];
  }
  for (i = 0; i < N; i++)
  {
    ulong g = 0;
    if (i)
    {
      long a, x = i, y = -2;
      do { y += 4; a = x%y; x = x/y; } while (!a);
      lswap(u[0],u[2]);
      switch (y)
      {
      case 2:
        break;
      case 6:
        lswap(u[4],u[6]);
        if (!(a & 1))
        {
          a = 4 - (a>>1);
          lswap(u[6], u[a]);
          lswap(u[4], u[a-2]);
        }
        break;
      case 10:
        x = u[6];
        u[6] = u[3];
        u[3] = u[2];
        u[2] = u[4];
        u[4] = u[1];
        u[1] = u[0];
        u[0] = x;
        if (a >= 3) a += 2;
        a = 8 - a;
        lswap(u[10],u[a]);
        lswap(u[8], u[a-2]);
        break;
      }
    }
    for (k = 0; k < n; k += 2) g += mael(MT,u[k],u[k+1]);
    if (headlongisint(g,n))
    {
      for (k = 0; k < n; k += 2)
      {
        pfu[u[k]] = u[k+1];
        pfu[u[k+1]] = u[k];
      }
      if (galois_test_perm(td, pfu)) break;
      hop++;
    }
    set_avma(av2);
  }
  if (i == N) return gc_NULL(ltop);
  if (DEBUGLEVEL >= 1 && hop)
    err_printf("A4GaloisConj: %ld hop over %ld iterations\n", hop, N);
  if (DEBUGLEVEL >= 4) err_printf("A4GaloisConj: tau=%Ps \n", pfu);
  set_avma(av2);
  orb = mkvec2(pft,pfu);
  O = vecperm_orbits(orb, 12);
  if (DEBUGLEVEL >= 4) {
    err_printf("A4GaloisConj: orb=%Ps\n", orb);
    err_printf("A4GaloisConj: O=%Ps \n", O);
  }
  av2 = avma;
  O1 = gel(O,1); O2 = gel(O,2); O3 = gel(O,3);
  for (j = 0; j < 2; j++)
  {
    pfv[O1[1]] = O2[1];
    pfv[O1[2]] = O2[3+j];
    pfv[O1[3]] = O2[4 - (j << 1)];
    pfv[O1[4]] = O2[2+j];
    for (i = 0; i < 4; i++)
    {
      ulong g = 0;
      switch (i)
      {
      case 0: break;
      case 1: lswap(O3[1], O3[2]); lswap(O3[3], O3[4]); break;
      case 2: lswap(O3[1], O3[4]); lswap(O3[2], O3[3]); break;
      case 3: lswap(O3[1], O3[2]); lswap(O3[3], O3[4]); break;
      }
      pfv[O2[1]]          = O3[1];
      pfv[O2[3+j]]        = O3[4-j];
      pfv[O2[4 - (j<<1)]] = O3[2 + (j<<1)];
      pfv[O2[2+j]]        = O3[3-j];
      pfv[O3[1]]          = O1[1];
      pfv[O3[4-j]]        = O1[2];
      pfv[O3[2 + (j<<1)]] = O1[3];
      pfv[O3[3-j]]        = O1[4];
      for (k = 1; k <= n; k++) g += mael(mt,k,pfv[k]);
      if (headlongisint(g,n) && galois_test_perm(td, pfv))
      {
        set_avma(av);
        if (DEBUGLEVEL >= 1)
          err_printf("A4GaloisConj: %ld hop over %d iterations max\n",
                     hop, 10395 + 68);
        return res;
      }
      hop++; set_avma(av2);
    }
  }
  return gc_NULL(ltop);
}

/* S4 */
static GEN
s4makelift(GEN u, struct galois_lift *gl)
{ return FpXQ_powers(u, degpol(gl->T)-1, gl->TQ, gl->Q); }

static long
s4test(GEN u, GEN liftpow, struct galois_lift *gl, GEN phi)
{
  pari_sp av = avma;
  GEN res, Q, Q2;
  long bl, i, d = lg(u)-2;
  pari_timer ti;
  if (DEBUGLEVEL >= 6) timer_start(&ti);
  if (!d) return 0;
  Q = gl->Q; Q2 = shifti(Q,-1);
  res = gel(u,2);
  for (i = 2; i <= d; i++)
    if (lg(gel(liftpow,i))>2)
      res = addii(res, mulii(gmael(liftpow,i,2), gel(u,i+1)));
  res = remii(res,Q);
  if (gl->den != gen_1) res = mulii(res, gl->den);
  res = centermodii(res, Q,Q2);
  if (abscmpii(res, gl->gb->bornesol) > 0) return gc_long(av,0);
  res = scalar_ZX_shallow(gel(u,2),varn(u));
  for (i = 2; i <= d ; i++)
    if (lg(gel(liftpow,i))>2)
      res = ZX_add(res, ZX_Z_mul(gel(liftpow,i), gel(u,i+1)));
  res = FpX_red(res, Q);
  if (gl->den != gen_1) res = FpX_Fp_mul(res, gl->den, Q);
  res = FpX_center_i(res, Q, shifti(Q,-1));
  bl = poltopermtest(res, gl, phi);
  if (DEBUGLEVEL >= 6) timer_printf(&ti, "s4test()");
  return gc_long(av,bl);
}

static GEN
s4releveauto(GEN M, GEN B, GEN T, GEN p,long a1,long a2,long a3,long a4,long a5,long a6)
{
  GEN       F = ZX_mul(gel(M,a1),gel(B,a2));
  F = ZX_add(F, ZX_mul(gel(M,a2),gel(B,a1)));
  F = ZX_add(F, ZX_mul(gel(M,a3),gel(B,a4)));
  F = ZX_add(F, ZX_mul(gel(M,a4),gel(B,a3)));
  F = ZX_add(F, ZX_mul(gel(M,a5),gel(B,a6)));
  F = ZX_add(F, ZX_mul(gel(M,a6),gel(B,a5)));
  return FpXQ_red(F, T, p);
}

static GEN
lincomb(GEN B, long a, long b, long j)
{
  long k = (-j) & 3;
  return ZX_add(gmael(B,a,j+1), gmael(B,b,k+1));
}

static GEN
FpXV_ffisom(GEN V, GEN p)
{
  pari_sp av = avma;
  long i, j, l = lg(V);
  GEN S = cgetg(l, t_VEC), Si = cgetg(l, t_VEC), M;
  for (i = 1; i < l; i++)
  {
    gel(S,i) = FpX_ffisom(gel(V,1), gel(V,i), p);
    gel(Si,i) = FpXQ_ffisom_inv(gel(S,i), gel(V,i), p);
  }
  M = cgetg(l, t_MAT);
  for (j = 1; j < l; j++)
    gel(M,j) = FpXC_FpXQ_eval(Si, gel(S,j), gel(V,j), p);
  return gerepileupto(av, M);
}

static GEN
mkliftpow(GEN x, GEN T, GEN p, struct galois_lift *gl)
{ pari_APPLY_same(automorphismlift(FpXV_chinese(gel(x,i), T, p, NULL), gl)) }

#define rot3(x,y,z) {long _t=x; x=y; y=z; z=_t;}
#define rot4(x,y,z,t) {long _u=x; x=y; y=z; z=t; t=_u;}

/* FIXME: could use the intheadlong technique */
static GEN
s4galoisgen(struct galois_lift *gl)
{
  const long n = 24;
  struct galois_testlift gt;
  pari_sp av, ltop2, ltop = avma;
  long i, j;
  GEN sigma, tau, phi, res, r1,r2,r3,r4, pj, p = gl->p, Q = gl->Q, TQ = gl->TQ;
  GEN sg, Tp, Tmod, misom, B, Bcoeff, liftpow, liftp, aut;

  res = cgetg(3, t_VEC);
  r1 = cgetg(n+1, t_VECSMALL);
  r2 = cgetg(n+1, t_VECSMALL);
  r3 = cgetg(n+1, t_VECSMALL);
  r4 = cgetg(n+1, t_VECSMALL);
  gel(res,1)= mkvec4(r1,r2,r3,r4);
  gel(res,2) = mkvecsmall4(2,2,3,2);
  ltop2 = avma;
  sg = identity_perm(6);
  pj = zero_zv(6);
  sigma = cgetg(n+1, t_VECSMALL);
  tau = cgetg(n+1, t_VECSMALL);
  phi = cgetg(n+1, t_VECSMALL);
  Tp = FpX_red(gl->T,p);
  Tmod = gel(FpX_factor(Tp,p), 1);
  misom = FpXV_ffisom(Tmod, p);
  aut = galoisdolift(gl);
  inittestlift(aut, Tmod, gl, &gt);
  B = FqC_FqV_mul(gt.pauto, gt.bezoutcoeff, gl->TQ, Q);
  Bcoeff = gt.bezoutcoeff;
  liftp = mkliftpow(shallowtrans(misom), Tmod, p, gl);
  av = avma;
  for (i = 0; i < 3; i++, set_avma(av))
  {
    pari_sp av1, av2, av3;
    GEN u, u1, u2, u3;
    long j1, j2, j3;
    if (i)
    {
      if (i == 1) { lswap(sg[2],sg[3]); }
      else        { lswap(sg[1],sg[3]); }
    }
    u = s4releveauto(liftp,Bcoeff,TQ,Q,sg[1],sg[2],sg[3],sg[4],sg[5],sg[6]);
    liftpow = s4makelift(u, gl);
    av1 = avma;
    for (j1 = 0; j1 < 4; j1++, set_avma(av1))
    {
      u1 = lincomb(B,sg[5],sg[6],j1);
      av2 = avma;
      for (j2 = 0; j2 < 4; j2++, set_avma(av2))
      {
        u2 = lincomb(B,sg[3],sg[4],j2);
        u2 = FpX_add(u1, u2, Q); av3 = avma;
        for (j3 = 0; j3 < 4; j3++, set_avma(av3))
        {
          u3 = lincomb(B,sg[1],sg[2],j3);
          u3 = FpX_add(u2, u3, Q);
          if (DEBUGLEVEL >= 4)
            err_printf("S4GaloisConj: Testing %d/3:%d/4:%d/4:%d/4:%Ps\n",
                       i,j1,j2,j3, sg);
          if (s4test(u3, liftpow, gl, sigma))
          {
            pj[1] = j3;
            pj[2] = j2;
            pj[3] = j1; goto suites4;
          }
        }
      }
    }
  }
  return gc_NULL(ltop);
suites4:
  if (DEBUGLEVEL >= 4) err_printf("S4GaloisConj: sigma=%Ps\n", sigma);
  if (DEBUGLEVEL >= 4) err_printf("S4GaloisConj: pj=%Ps\n", pj);
  set_avma(av);
  for (j = 1; j <= 3; j++)
  {
    pari_sp av2, av3;
    GEN u;
    long w, l;
    rot3(sg[1], sg[3], sg[5])
    rot3(sg[2], sg[4], sg[6])
    rot3(pj[1], pj[2], pj[3])
    for (l = 0; l < 2; l++, set_avma(av))
    {
      u = s4releveauto(liftp,Bcoeff,TQ,Q,sg[1],sg[3],sg[2],sg[4],sg[5],sg[6]);
      liftpow = s4makelift(u, gl);
      av2 = avma;
      for (w = 0; w < 4; w += 2, set_avma(av2))
      {
        GEN uu;
        pj[6] = (w + pj[3]) & 3;
        uu = lincomb(B, sg[5], sg[6], pj[6]);
        uu = FpX_red(uu, Q);
        av3 = avma;
        for (i = 0; i < 4; i++, set_avma(av3))
        {
          GEN u;
          pj[4] = i;
          pj[5] = (i + pj[2] - pj[1]) & 3;
          if (DEBUGLEVEL >= 4)
            err_printf("S4GaloisConj: Testing %d/3:%d/2:%d/2:%d/4:%Ps:%Ps\n",
                       j-1, w >> 1, l, i, sg, pj);
          u = ZX_add(lincomb(B,sg[1],sg[3],pj[4]),
                     lincomb(B,sg[2],sg[4],pj[5]));
          u = FpX_add(uu,u,Q);
          if (s4test(u, liftpow, gl, tau)) goto suites4_2;
        }
      }
      lswap(sg[3], sg[4]);
      pj[2] = (-pj[2]) & 3;
    }
  }
  return gc_NULL(ltop);
suites4_2:
  set_avma(av);
  {
    long abc = (pj[1] + pj[2] + pj[3]) & 3;
    long abcdef = ((abc + pj[4] + pj[5] - pj[6]) & 3) >> 1;
    GEN u;
    pari_sp av2;
    u = s4releveauto(liftp,Bcoeff,TQ,Q,sg[1],sg[4],sg[2],sg[5],sg[3],sg[6]);
    liftpow = s4makelift(u, gl);
    av2 = avma;
    for (j = 0; j < 8; j++)
    {
      long h, g, i;
      h = j & 3;
      g = (abcdef + ((j & 4) >> 1)) & 3;
      i = (h + abc - g) & 3;
      u = ZX_add(lincomb(B,sg[1],sg[4], g), lincomb(B,sg[2],sg[5], h));
      u = FpX_add(u, lincomb(B,sg[3],sg[6], i),Q);
      if (DEBUGLEVEL >= 4)
        err_printf("S4GaloisConj: Testing %d/8 %d:%d:%d\n", j,g,h,i);
      if (s4test(u, liftpow, gl, phi)) break;
      set_avma(av2);
    }
  }
  if (j == 8) return gc_NULL(ltop);
  for (i = 1; i <= n; i++)
  {
    r1[i] = sigma[tau[i]];
    r2[i] = phi[sigma[tau[phi[i]]]];
    r3[i] = phi[sigma[i]];
    r4[i] = sigma[i];
  }
  set_avma(ltop2); return res;
}

static GEN
f36releveauto2(GEN Bl, GEN T, GEN p,GEN a)
{
  GEN      F = gmael(Bl,a[1],a[1]);
  F = ZX_add(F,gmael(Bl,a[2],a[3]));
  F = ZX_add(F,gmael(Bl,a[3],a[2]));
  F = ZX_add(F,gmael(Bl,a[4],a[5]));
  F = ZX_add(F,gmael(Bl,a[5],a[4]));
  F = ZX_add(F,gmael(Bl,a[6],a[7]));
  F = ZX_add(F,gmael(Bl,a[7],a[6]));
  F = ZX_add(F,gmael(Bl,a[8],a[9]));
  F = ZX_add(F,gmael(Bl,a[9],a[8]));
  return FpXQ_red(F, T, p);
}

static GEN
f36releveauto4(GEN Bl, GEN T, GEN p,GEN a)
{
  GEN      F = gmael(Bl,a[1],a[1]);
  F = ZX_add(F,gmael(Bl,a[2],a[3]));
  F = ZX_add(F,gmael(Bl,a[3],a[4]));
  F = ZX_add(F,gmael(Bl,a[4],a[5]));
  F = ZX_add(F,gmael(Bl,a[5],a[2]));
  F = ZX_add(F,gmael(Bl,a[6],a[7]));
  F = ZX_add(F,gmael(Bl,a[7],a[8]));
  F = ZX_add(F,gmael(Bl,a[8],a[9]));
  F = ZX_add(F,gmael(Bl,a[9],a[6]));
  return FpXQ_red(F, T, p);
}

static GEN
f36galoisgen(struct galois_lift *gl)
{
  const long n = 36;
  struct galois_testlift gt;
  pari_sp av, ltop2, ltop = avma;
  long i;
  GEN sigma, tau, rho, res, r1,r2,r3, pj, pk, p = gl->p, Q = gl->Q, TQ = gl->TQ;
  GEN sg, s4, sp,  Tp, Tmod, misom, Bcoeff, liftpow, aut, liftp, B, Bl, tam;
  res = cgetg(3, t_VEC);
  r1 = cgetg(n+1, t_VECSMALL);
  r2 = cgetg(n+1, t_VECSMALL);
  r3 = cgetg(n+1, t_VECSMALL);
  gel(res,1)= mkvec3(r1,r2,r3);
  gel(res,2) = mkvecsmall3(3,3,4);
  ltop2 = avma;
  sg = identity_perm(9);
  s4 = identity_perm(9);
  sp = identity_perm(9);
  pj = zero_zv(4);
  pk = zero_zv(2);
  sigma = cgetg(n+1, t_VECSMALL);
  tau = r3;
  rho = cgetg(n+1, t_VECSMALL);
  Tp = FpX_red(gl->T,p);
  Tmod = gel(FpX_factor(Tp,p), 1);
  misom = FpXV_ffisom(Tmod, p);
  aut = galoisdolift(gl);
  inittestlift(aut, Tmod, gl, &gt);
  Bcoeff = gt.bezoutcoeff;
  B = FqC_FqV_mul(gt.pauto, Bcoeff, gl->TQ, gl->Q);
  liftp = mkliftpow(shallowtrans(misom), Tmod, p, gl);
  Bl = FqC_FqV_mul(liftp,Bcoeff, gl->TQ, gl->Q);
  av = avma;
  for (i = 0; i < 105; i++, set_avma(av))
  {
    pari_sp av0, av1, av2, av3;
    GEN u0, u1, u2, u3;
    long j0, j1, j2, j3, s;
    if (i)
    {
      rot3(sg[7],sg[8],sg[9])
      if (i%3==0)
      {
        s=sg[5]; sg[5]=sg[6]; sg[6]=sg[7]; sg[7]=sg[8]; sg[8]=sg[9]; sg[9]=s;
        if (i%15==0)
        {
          s=sg[3]; sg[3]=sg[4]; sg[4]=sg[5];
          sg[5]=sg[6]; sg[6]=sg[7]; sg[7]=sg[8]; sg[8]=sg[9]; sg[9]=s;
        }
      }
    }
    liftpow = s4makelift(f36releveauto2(Bl, TQ, Q, sg), gl);
    av0 = avma;
    for (j0 = 0; j0 < 4; j0++, set_avma(av0))
    {
      u0 = lincomb(B,sg[8],sg[9],j0);
      u0 = FpX_add(u0, gmael(B,sg[1],3), Q); av1 = avma;
      for (j1 = 0; j1 < 4; j1++, set_avma(av1))
      {
        u1 = lincomb(B,sg[6],sg[7],j1);
        u1 = FpX_add(u0, u1, Q); av2 = avma;
        for (j2 = 0; j2 < 4; j2++, set_avma(av2))
        {
          u2 = lincomb(B,sg[4],sg[5],j2);
          u2 = FpX_add(u1, u2, Q); av3 = avma;
          for (j3 = 0; j3 < 4; j3++, set_avma(av3))
          {
            u3 = lincomb(B,sg[2],sg[3],j3);
            u3 = FpX_add(u2, u3, Q);
            if (s4test(u3, liftpow, gl, sigma))
            {
              pj[1] = j3;
              pj[2] = j2;
              pj[3] = j1;
              pj[4] = j0; goto suitef36;
            }
          }
        }
      }
    }
  }
  return gc_NULL(ltop);
suitef36:
  s4[1]=sg[1]; s4[2]=sg[2]; s4[4]=sg[3];
  s4[3]=sg[4]; s4[5]=sg[5]; s4[6]=sg[6];
  s4[8]=sg[7]; s4[7]=sg[8]; s4[9]=sg[9];
  for (i = 0; i < 12; i++, set_avma(av))
  {
    pari_sp av0, av1;
    GEN u0, u1;
    long j0, j1;
    if (i)
    {
      lswap(s4[3],s4[5]); pj[2] = (-pj[2])&3;
      if (odd(i)) { lswap(s4[7],s4[9]); pj[4]=(-pj[4])&3; }
      if (i%4==0)
      {
        rot3(s4[3],s4[6],s4[7]);
        rot3(s4[5],s4[8],s4[9]);
        rot3(pj[2],pj[3],pj[4]);
      }
    }
    liftpow = s4makelift(f36releveauto4(Bl, TQ, Q, s4), gl);
    av0 = avma;
    for (j0 = 0; j0 < 4; j0++, set_avma(av0))
    {
      u0 = FpX_add(gmael(B,s4[1],2), gmael(B,s4[2],1+j0),Q);
      u0 = FpX_add(u0, gmael(B,s4[3],1+smodss(pj[2]-j0,4)),Q);
      u0 = FpX_add(u0, gmael(B,s4[4],1+smodss(j0-pj[1]-pj[2],4)),Q);
      u0 = FpX_add(u0, gmael(B,s4[5],1+smodss(pj[1]-j0,4)),Q);
      av1 = avma;
      for (j1 = 0; j1 < 4; j1++, set_avma(av1))
      {
        u1 = FpX_add(u0, gmael(B,s4[6],1+j1),Q);
        u1 = FpX_add(u1, gmael(B,s4[7],1+smodss(pj[4]-j1,4)),Q);
        u1 = FpX_add(u1, gmael(B,s4[8],1+smodss(j1-pj[3]-pj[4],4)),Q);
        u1 = FpX_add(u1, gmael(B,s4[9],1+smodss(pj[3]-j1,4)),Q);
        if (s4test(u1, liftpow, gl, tau))
        {
          pk[1] = j0;
          pk[2] = j1; goto suitef36_2;
        }
      }
    }
  }
  return gc_NULL(ltop);
suitef36_2:
  sp[1]=s4[9]; sp[2]=s4[1]; sp[3]=s4[2];
  sp[4]=s4[7]; sp[5]=s4[3]; sp[6]=s4[8];
  sp[8]=s4[4]; sp[7]=s4[5]; sp[9]=s4[6];
  for (i = 0; i < 4; i++, set_avma(av))
  {
    const int w[4][6]={{0,0,1,3,0,2},{1,0,2,1,1,2},{3,3,2,0,3,1},{0,1,3,0,0,3}};
    pari_sp av0, av1, av2;
    GEN u0, u1, u2;
    long j0, j1,j2,j3,j4,j5;
    if (i)
    {
      rot4(sp[3],sp[5],sp[8],sp[7])
      pk[1]=(-pk[1])&3;
    }
    liftpow = s4makelift(f36releveauto4(Bl,TQ,Q,sp), gl);
    av0 = avma;
    for (j0 = 0; j0 < 4; j0++, set_avma(av0))
    {
      u0 = FpX_add(gmael(B,sp[1],2), gmael(B,sp[2],1+j0),Q);
      av1 = avma;
      for (j1 = 0; j1 < 4; j1++, set_avma(av1))
      {
        u1 = FpX_add(u0, gmael(B,sp[3],1+j1),Q);
        j3 = (-pk[1]-pj[3]+j0+j1-w[i][0]*pj[1]-w[i][3]*pj[2])&3;
        u1 = FpX_add(u1, gmael(B,sp[6],1+j3),Q);
        j5 = (-pk[1]+2*j0+2*j1-w[i][2]*pj[1]-w[i][5]*pj[2])&3;
        u1 = FpX_add(u1, gmael(B,sp[8],1+j5),Q);
        av2 = avma;
        for (j2 = 0; j2 < 4; j2++, set_avma(av2))
        {
          u2 = FpX_add(u1, gmael(B,sp[4],1+j2),Q);
          u2 = FpX_add(u2, gmael(B,sp[5],1+smodss(-j0-j1-j2,4)),Q);
          j4 = (-pk[1]-pk[2]+pj[3]+pj[4]-j2-w[i][1]*pj[1]-w[i][4]*pj[2])&3;
          u2 = FpX_add(u2, gmael(B,sp[7],1+j4),Q);
          u2 = FpX_add(u2, gmael(B,sp[9],1+smodss(-j3-j4-j5,4)),Q);
          if (s4test(u2, liftpow, gl, rho))
            goto suitef36_3;
        }
      }
    }
  }
  return gc_NULL(ltop);
suitef36_3:
  tam = perm_inv(tau);
  for (i = 1; i <= n; i++)
  {
    r1[tau[i]] = rho[i];
    r2[i] = tam[rho[i]];
  }
  set_avma(ltop2); return res;
}

/* return a vecvecsmall */
static GEN
galoisfindgroups(GEN lo, GEN sg, long f)
{
  pari_sp ltop = avma;
  long i, j, k;
  GEN V = cgetg(lg(lo), t_VEC);
  for(j=1,i=1; i<lg(lo); i++)
  {
    pari_sp av = avma;
    GEN loi = gel(lo,i), W = cgetg(lg(loi),t_VECSMALL);
    for (k=1; k<lg(loi); k++) W[k] = loi[k] % f;
    W = vecsmall_uniq(W);
    if (zv_equal(W, sg)) gel(V,j++) = loi;
    set_avma(av);
  }
  setlg(V,j); return gerepilecopy(ltop,V);
}

static GEN
galoismakepsi(long g, GEN sg, GEN pf)
{
  GEN psi=cgetg(g+1,t_VECSMALL);
  long i;
  for (i = 1; i < g; i++) psi[i] = sg[pf[i]];
  psi[g] = sg[1]; return psi;
}

static GEN
galoisfrobeniuslift_nilp(GEN T, GEN den, GEN L,  GEN Lden,
    struct galois_frobenius *gf,  struct galois_borne *gb)
{
  pari_sp ltop=avma, av2;
  struct galois_lift gl;
  long i, k, deg = 1, g = lg(gf->Tmod)-1;
  GEN F,Fp,Fe, aut, frob, res = cgetg(lg(L), t_VECSMALL);
  gf->psi = const_vecsmall(g,1);
  av2 = avma;
  initlift(T, den, gf->p, L, Lden, gb, &gl);
  if (DEBUGLEVEL >= 4)
    err_printf("GaloisConj: p=%ld e=%ld deg=%ld fp=%ld\n",
                            gf->p, gl.e, deg, gf->fp);
  aut = galoisdolift(&gl);
  if (galoisfrobeniustest(aut,&gl,res))
  {
    set_avma(av2); gf->deg = gf->fp; return res;
  }

  F =factoru(gf->fp);
  Fp = gel(F,1);
  Fe = gel(F,2);
  frob = cgetg(lg(L), t_VECSMALL);
  for(k = lg(Fp)-1; k>=1; k--)
  {
    pari_sp btop=avma;
    GEN fres=NULL;
    long el = gf->fp, dg = 1, dgf = 1, e, pr;
    for(e=1; e<=Fe[k]; e++)
    {
      dg *= Fp[k]; el /= Fp[k];
      if (DEBUGLEVEL>=4) err_printf("Trying degre %d.\n",dg);
      if (el==1) break;
      aut = galoisdoliftn(&gl, el);
      if (!galoisfrobeniustest(aut,&gl,frob))
        break;
      dgf = dg; fres = gcopy(frob);
    }
    if (dgf == 1) { set_avma(btop); continue; }
    pr = deg*dgf;
    if (deg == 1)
    {
      for(i=1;i<lg(res);i++) res[i]=fres[i];
    }
    else
    {
      GEN cp = perm_mul(res,fres);
      for(i=1;i<lg(res);i++) res[i] = cp[i];
    }
    deg = pr; set_avma(btop);
  }
  if (DEBUGLEVEL>=4 && res) err_printf("Best lift: %d\n",deg);
  if (deg==1) return gc_NULL(ltop);
  else
  {
    set_avma(av2); gf->deg = deg; return res;
  }
}

static GEN
galoisfrobeniuslift(GEN T, GEN den, GEN L,  GEN Lden,
    struct galois_frobenius *gf,  struct galois_borne *gb)
{
  pari_sp ltop=avma, av2;
  struct galois_testlift gt;
  struct galois_lift gl;
  long i, j, k, n = lg(L)-1, deg = 1, g = lg(gf->Tmod)-1;
  GEN F,Fp,Fe, aut, frob, res = cgetg(lg(L), t_VECSMALL);
  gf->psi = const_vecsmall(g,1);
  av2 = avma;
  initlift(T, den, gf->p, L, Lden, gb, &gl);
  if (DEBUGLEVEL >= 4)
    err_printf("GaloisConj: p=%ld e=%ld deg=%ld fp=%ld\n",
                            gf->p, gl.e, deg, gf->fp);
  aut = galoisdolift(&gl);
  if (galoisfrobeniustest(aut,&gl,res))
  {
    set_avma(av2); gf->deg = gf->fp; return res;
  }
  inittestlift(aut,gf->Tmod, &gl, &gt);
  gt.C = cgetg(gf->fp+1,t_VEC);
  gt.Cd= cgetg(gf->fp+1,t_VEC);
  for (i = 1; i <= gf->fp; i++) {
    gel(gt.C,i)  = zero_zv(gt.g);
    gel(gt.Cd,i) = zero_zv(gt.g);
  }

  F =factoru(gf->fp);
  Fp = gel(F,1);
  Fe = gel(F,2);
  frob = cgetg(lg(L), t_VECSMALL);
  for(k=lg(Fp)-1;k>=1;k--)
  {
    pari_sp btop=avma;
    GEN psi=NULL, fres=NULL, sg = identity_perm(1);
    long el=gf->fp, dg=1, dgf=1, e, pr;
    for(e=1; e<=Fe[k]; e++)
    {
      GEN lo, pf;
      long l;
      dg *= Fp[k]; el /= Fp[k];
      if (DEBUGLEVEL>=4) err_printf("Trying degre %d.\n",dg);
      if (galoisfrobeniustest(gel(gt.pauto,el+1),&gl,frob))
      {
        psi = const_vecsmall(g,1);
        dgf = dg; fres = leafcopy(frob); continue;
      }
      lo = listznstarelts(dg, n / gf->fp);
      if (e!=1) lo = galoisfindgroups(lo, sg, dgf);
      if (DEBUGLEVEL>=4) err_printf("Galoisconj:Subgroups list:%Ps\n", lo);
      for (l = 1; l < lg(lo); l++)
        if (lg(gel(lo,l))>2 && frobeniusliftall(gel(lo,l), el, &pf,&gl,&gt, frob))
        {
          sg  = leafcopy(gel(lo,l));
          psi = galoismakepsi(g,sg,pf);
          dgf = dg; fres = leafcopy(frob); break;
        }
      if (l == lg(lo)) break;
    }
    if (dgf == 1) { set_avma(btop); continue; }
    pr = deg*dgf;
    if (deg == 1)
    {
      for(i=1;i<lg(res);i++) res[i]=fres[i];
      for(i=1;i<lg(psi);i++) gf->psi[i]=psi[i];
    }
    else
    {
      GEN cp = perm_mul(res,fres);
      for(i=1;i<lg(res);i++) res[i] = cp[i];
      for(i=1;i<lg(psi);i++) gf->psi[i] = (dgf*gf->psi[i] + deg*psi[i]) % pr;
    }
    deg = pr; set_avma(btop);
  }
  for (i = 1; i <= gf->fp; i++)
    for (j = 1; j <= gt.g; j++) guncloneNULL(gmael(gt.C,i,j));
  if (DEBUGLEVEL>=4 && res) err_printf("Best lift: %d\n",deg);
  if (deg==1) return gc_NULL(ltop);
  else
  {
    /* Normalize result so that psi[g]=1 */
    ulong im = Fl_inv(gf->psi[g], deg);
    GEN cp = perm_powu(res, im);
    for(i=1;i<lg(res);i++) res[i] = cp[i];
    for(i=1;i<lg(gf->psi);i++) gf->psi[i] = Fl_mul(im, gf->psi[i], deg);
    set_avma(av2); gf->deg = deg; return res;
  }
}

/* return NULL if not Galois */
static GEN
galoisfindfrobenius(GEN T, GEN L, GEN den, GEN bad, struct galois_frobenius *gf,
    struct galois_borne *gb, const struct galois_analysis *ga)
{
  pari_sp ltop = avma, av;
  long Try = 0, n = degpol(T), deg, gmask = (ga->group&ga_ext_2)? 3: 1;
  GEN frob, Lden = makeLden(L,den,gb);
  long is_nilpotent = ga->group&ga_all_nilpotent;
  forprime_t S;

  u_forprime_init(&S, ga->p, ULONG_MAX);
  av = avma;
  deg = gf->deg = ga->deg;
  while ((gf->p = u_forprime_next(&S)))
  {
    pari_sp lbot;
    GEN Ti, Tp;
    long nb, d;
    set_avma(av);
    Tp = ZX_to_Flx(T, gf->p);
    if (!Flx_is_squarefree(Tp, gf->p)) continue;
    if (bad && dvdiu(bad, gf->p)) continue;
    Ti = gel(Flx_factor(Tp, gf->p), 1);
    nb = lg(Ti)-1; d = degpol(gel(Ti,1));
    if (nb > 1 && degpol(gel(Ti,nb)) != d) return gc_NULL(ltop);
    if (((gmask&1)==0 || d % deg) && ((gmask&2)==0 || odd(d))) continue;
    if (DEBUGLEVEL >= 1) err_printf("GaloisConj: Trying p=%ld\n", gf->p);
    FlxV_to_ZXV_inplace(Ti);
    gf->fp = d;
    gf->Tmod = Ti; lbot = avma;
    if (is_nilpotent)
      frob = galoisfrobeniuslift_nilp(T, den, L, Lden, gf, gb);
    else
      frob = galoisfrobeniuslift(T, den, L, Lden, gf, gb);
    if (frob)
    {
      GEN *gptr[3];
      gf->Tmod = gcopy(Ti);
      gptr[0]=&gf->Tmod; gptr[1]=&gf->psi; gptr[2]=&frob;
      gerepilemanysp(ltop,lbot,gptr,3); return frob;
    }
    if (is_nilpotent) continue;
    if ((ga->group&ga_all_normal) && d % deg == 0) gmask &= ~1;
    /* The first prime degree is always divisible by deg, so we don't
     * have to worry about ext_2 being used before regular supersolvable*/
    if (!gmask) return gc_NULL(ltop);
    if ((ga->group&ga_non_wss) && ++Try > ((3*n)>>1))
    {
      if (DEBUGLEVEL)
        pari_warn(warner,"Galois group probably not weakly super solvable");
      return NULL;
    }
  }
  if (!gf->p) pari_err_OVERFLOW("galoisfindfrobenius [ran out of primes]");
  return NULL;
}

/* compute g such that tau(Pmod[#])= tau(Pmod[g]) */

static long
get_image(GEN tau, GEN P, GEN Pmod, GEN p)
{
  pari_sp av = avma;
  long g, gp = lg(Pmod)-1;
  tau = RgX_to_FpX(tau, p);
  tau = FpX_FpXQ_eval(gel(Pmod, gp), tau, P, p);
  tau = FpX_normalize(FpX_gcd(P, tau, p), p);
  for (g = 1; g <= gp; g++)
    if (ZX_equal(tau, gel(Pmod,g))) return gc_long(av,g);
  return gc_long(av,0);
}

static GEN
gg_get_std(GEN G)
{
  return !G ? NULL: lg(G)==3 ? G: mkvec2(gel(G,1),gmael(G,5,1));
}

static GEN galoisgen(GEN T, GEN L, GEN M, GEN den, GEN bad, struct galois_borne *gb,
          const struct galois_analysis *ga);

static GEN
galoisgenfixedfield(GEN Tp, GEN Pmod, GEN PL, GEN P, GEN ip, GEN bad, struct galois_borne *gb)
{
  GEN  Pden, PM;
  GEN  tau, PG, Pg;
  long g, lP;
  long x = varn(Tp);
  GEN Pp = FpX_red(P, ip);
  if (DEBUGLEVEL>=6)
    err_printf("GaloisConj: Fixed field %Ps\n",P);
  if (degpol(P)==2 && !bad)
  {
    PG=cgetg(3,t_VEC);
    gel(PG,1) = mkvec( mkvecsmall2(2,1) );
    gel(PG,2) = mkvecsmall(2);
    tau = deg1pol_shallow(gen_m1, negi(gel(P,3)), x);
    g = get_image(tau, Pp, Pmod, ip);
    if (!g) return NULL;
    Pg = mkvecsmall(g);
  }
  else
  {
    struct galois_analysis Pga;
    struct galois_borne Pgb;
    GEN mod, mod2;
    long j;
    if (!galoisanalysis(P, &Pga, 0, NULL)) return NULL;
    if (bad) Pga.group &= ~ga_easy;
    Pgb.l = gb->l;
    Pden = galoisborne(P, NULL, &Pgb, degpol(P));

    if (Pgb.valabs > gb->valabs)
    {
      if (DEBUGLEVEL>=4)
        err_printf("GaloisConj: increase prec of p-adic roots of %ld.\n"
            ,Pgb.valabs-gb->valabs);
      PL = ZpX_liftroots(P,PL,gb->l,Pgb.valabs);
    }
    else if (Pgb.valabs < gb->valabs)
      PL = FpC_red(PL, Pgb.ladicabs);
    PM = FpV_invVandermonde(PL, Pden, Pgb.ladicabs);
    PG = galoisgen(P, PL, PM, Pden, bad ? lcmii(Pgb.dis, bad): NULL, &Pgb, &Pga);
    if (!PG) return NULL;
    lP = lg(gel(PG,1));
    mod = Pgb.ladicabs; mod2 = shifti(mod, -1);
    Pg = cgetg(lP, t_VECSMALL);
    for (j = 1; j < lP; j++)
    {
      pari_sp btop=avma;
      tau = permtopol(gmael(PG,1,j), PL, PM, Pden, mod, mod2, x);
      g = get_image(tau, Pp, Pmod, ip);
      if (!g) return NULL;
      Pg[j] = g;
      set_avma(btop);
    }
  }
  return mkvec2(PG,Pg);
}

static GEN
galoisgenfixedfield0(GEN O, GEN L, GEN sigma, GEN T, GEN bad, GEN *pt_V,
                     struct galois_frobenius *gf, struct galois_borne *gb)
{
  pari_sp btop = avma;
  long vT = varn(T);
  GEN mod = gb->ladicabs, mod2 = shifti(gb->ladicabs,-1);
  GEN OL, sym, P, PL, p, Tp, Sp, Pmod, PG;
  OL = fixedfieldorbits(O,L);
  sym  = fixedfieldsympol(OL, itou(gb->l));
  PL = sympol_eval(sym, OL, mod);
  P = FpX_center_i(FpV_roots_to_pol(PL, mod, vT), mod, mod2);
  if (!FpX_is_squarefree(P,utoipos(gf->p)))
  {
    GEN badp = lcmii(bad? bad: gb->dis, ZX_disc(P));
    gf->p  = findpsi(badp, gf->p, T, sigma, gf->deg, &gf->Tmod, &gf->psi);
  }
  p  = utoipos(gf->p);
  Tp = FpX_red(T,p);
  Sp = sympol_aut_evalmod(sym, gf->deg, sigma, Tp, p);
  Pmod = fixedfieldfactmod(Sp, p, gf->Tmod);
  PG = galoisgenfixedfield(Tp, Pmod, PL, P, p, bad, gb);
  if (PG == NULL) return NULL;
  if (DEBUGLEVEL >= 4)
    err_printf("GaloisConj: Back to Earth:%Ps\n", gg_get_std(gel(PG,1)));
  if (pt_V) *pt_V = mkvec3(sym, PL, P);
  return gc_all(btop, pt_V ? 4: 3, &PG, &gf->Tmod, &gf->psi, pt_V);
}

/* Let sigma^m=1, tau*sigma*tau^-1=sigma^s. Return n = sum_{0<=k<e,0} s^k mod m
 * so that (sigma*tau)^e = sigma^n*tau^e. N.B. n*(1-s) = 1-s^e mod m,
 * unfortunately (1-s) may not invertible mod m */
static long
stpow(long s, long e, long m)
{
  long i, n = 1;
  for (i = 1; i < e; i++) n = (1 + n * s) % m;
  return n;
}

static GEN
wpow(long s, long m, long e, long n)
{
  GEN   w = cgetg(n+1,t_VECSMALL);
  long si = s;
  long i;
  w[1] = 1;
  for(i=2; i<=n; i++) w[i] = w[i-1]*e;
  for(i=n; i>=1; i--)
  {
    si = Fl_powu(si,e,m);
    w[i] = Fl_mul(s-1, stpow(si, w[i], m), m);
  }
  return w;
}

static GEN
galoisgenliftauto(GEN O, GEN gj, long s, long n, struct galois_test *td)
{
  pari_sp av = avma;
  long sr, k;
  long deg = lg(gel(O,1))-1;
  GEN  X  = cgetg(lg(O), t_VECSMALL);
  GEN  oX = cgetg(lg(O), t_VECSMALL);
  GEN  B  = perm_cycles(gj);
  long oj = lg(gel(B,1)) - 1;
  GEN  F  = factoru(oj);
  GEN  Fp = gel(F,1);
  GEN  Fe = gel(F,2);
  GEN  pf = identity_perm(n);
  if (DEBUGLEVEL >= 6)
    err_printf("GaloisConj: %Ps of relative order %d\n", gj, oj);
  for (k=lg(Fp)-1; k>=1; k--)
  {
    long f, dg = 1, el = oj, osel = 1, a = 0;
    long p  = Fp[k], e  = Fe[k], op = oj / upowuu(p,e);
    long i;
    GEN  pf1 = NULL, w, wg, Be = cgetg(e+1,t_VEC);
    gel(Be,e) = cyc_pow(B, op);
    for(i=e-1; i>=1; i--) gel(Be,i) = cyc_pow(gel(Be,i+1), p);
    w = wpow(Fl_powu(s,op,deg),deg,p,e);
    wg = cgetg(e+2,t_VECSMALL);
    wg[e+1] = deg;
    for (i=e; i>=1; i--) wg[i] = ugcd(wg[i+1], w[i]);
    for (i=1; i<lg(O); i++) oX[i] = 0;
    for (f=1; f<=e; f++)
    {
      long sel, t;
      GEN Bel = gel(Be,f);
      dg *= p; el /= p;
      sel = Fl_powu(s,el,deg);
      if (DEBUGLEVEL >= 6) err_printf("GaloisConj: B=%Ps\n", Bel);
      sr  = ugcd(stpow(sel,p,deg),deg);
      if (DEBUGLEVEL >= 6)
        err_printf("GaloisConj: exp %d: s=%ld [%ld] a=%ld w=%ld wg=%ld sr=%ld\n",
            dg, sel, deg, a, w[f], wg[f+1], sr);
      for (t = 0; t < sr; t++)
        if ((a+t*w[f])%wg[f+1]==0)
        {
          long i, j, k, st;
          for (i = 1; i < lg(X); i++) X[i] = 0;
          for (i = 0; i < lg(X)-1; i+=dg)
            for (j = 1, k = p, st = t; k <= dg; j++, k += p)
            {
              X[k+i] = (oX[j+i] + st)%deg;
              st = (t + st*osel)%deg;
            }
          pf1 = testpermutation(O, Bel, X, sel, p, sr, td);
          if (pf1) break;
        }
      if (!pf1) return NULL;
      for (i=1; i<lg(O); i++) oX[i] = X[i];
      osel = sel; a = (a+t*w[f])%deg;
    }
    pf = perm_mul(pf, perm_powu(pf1, el));
  }
  return gerepileuptoleaf(av, pf);
}

static GEN
FlxV_Flx_gcd(GEN x, GEN T, ulong p)
{ pari_APPLY_same(Flx_normalize(Flx_gcd(gel(x,i),T,p),p)) }

static GEN
Flx_FlxV_minpolymod(GEN y, GEN x, ulong p)
{ pari_APPLY_same(Flxq_minpoly(Flx_rem(y, gel(x,i), p), gel(x,i), p)) }

static GEN
FlxV_minpolymod(GEN x, GEN y, ulong p)
{ pari_APPLY_same(Flx_FlxV_minpolymod(gel(x,i), y, p)) }

static GEN
factperm(GEN x)
{
  pari_APPLY_same(gen_indexsort(gel(x,i), (void*)cmp_Flx, cmp_nodata))
}

/* compute (prod p_i^e_i)(1) */

static long
permprodeval(GEN p, GEN e, long s)
{
  long i, j, l = lg(p);
  for (i=l-1; i>=1; i--)
  {
    GEN pi = gel(p,i);
    long ei = uel(e,i);
    for(j = 1; j <= ei; j++)
      s = uel(pi, s);
  }
  return s;
}

static GEN
pc_to_perm(GEN pc, GEN gen, long n)
{
  long i, l = lg(pc);
  GEN s = identity_perm(n);
  for (i=1; i<l; i++)
    s = perm_mul(gel(gen,pc[i]),s);
  return s;
}

static GEN
genorbit(GEN ordH, GEN permfact_Hp, long fr, long n, long k, long j)
{
  pari_sp av = avma;
  long l = lg(gel(permfact_Hp,1))-1, no = 1, b, i;
  GEN W = zero_zv(l);
  GEN orb = cgetg(l+1, t_VECSMALL);
  GEN gen = cgetg(l+1, t_VEC);
  GEN E = cgetg(k+1, t_VECSMALL);
  for(b = 0; b < n; b++)
  {
    long bb = b, s;
    for(i = 1; i <= k; i++)
    {
      uel(E,i) = bb % uel(ordH,i);
      bb /= uel(ordH,i);
    }
    if (E[j]) continue;
    s = permprodeval(permfact_Hp, E, fr);
    if (s>lg(W)-1) pari_err_BUG("W1");
    if (W[s]) continue;
    W[s] = 1;
    if (no > l) pari_err_BUG("genorbit");
    uel(orb,no) = s;
    gel(gen,no) = zv_copy(E);
    no++;
  }
  if(no<l) pari_err_BUG("genorbit");
  return gerepilecopy(av, mkvec2(orb,gen));
}

INLINE GEN br_get(GEN br, long i, long j) { return gmael(br,j,i-j); }
static GEN pcgrp_get_ord(GEN G) { return gel(G,1); }
static GEN pcgrp_get_pow(GEN G) { return gel(G,2); }
static GEN pcgrp_get_br(GEN G)  { return gel(G,3); }

static GEN
cyclic_pc(long n)
{
  return mkvec3(mkvecsmall(n),mkvec(cgetg(1,t_VECSMALL)), mkvec(cgetg(1,t_VEC)));
}

static GEN
pc_normalize(GEN g, GEN G)
{
  long i, l = lg(g)-1, o = 1;
  GEN ord = pcgrp_get_ord(G), pw = pcgrp_get_pow(G), br = pcgrp_get_br(G);
  for (i = 1; i < l; i++)
  {
    if (g[i] == g[i+1])
    {
      if (++o == ord[g[i]])
      {
        GEN v = vecsmall_concat(vecslice(g,1,i-o+1),gel(pw,g[i]));
        GEN w = vecsmall_concat(v,vecslice(g,i+2,l));
        return pc_normalize(w, G);
      }
    }
    else if (g[i] > g[i+1])
    {
      GEN v = vecsmall_concat(vecslice(g,1,i-1), br_get(br,g[i],g[i+1]));
      GEN w = vecsmall_concat(mkvecsmall2(g[i+1],g[i]),vecslice(g,i+2,l));
      v = vecsmall_concat(v, w);
      return pc_normalize(v, G);
    }
    else o = 1;
  }
  return g;
}

static GEN
pc_inv(GEN g, GEN G)
{
  long i, l = lg(g);
  GEN ord = pcgrp_get_ord(G), pw  = pcgrp_get_pow(G);
  GEN v = cgetg(l, t_VEC);
  if (l==1) return v;
  for(i = 1; i < l; i++)
  {
    ulong gi = uel(g,i);
    gel(v,l-i) = vecsmall_concat(pc_inv(gel(pw, gi), G),
                                 const_vecsmall(uel(ord,gi)-1,gi));
  }
  return pc_normalize(shallowconcat1(v), G);
}

static GEN
pc_mul(GEN g, GEN h, GEN G)
{
  return pc_normalize(vecsmall_concat(g,h), G);
}

static GEN
pc_bracket(GEN g, GEN h, GEN G)
{
  GEN gh = pc_mul(g, h, G);
  GEN hg = pc_mul(h, g, G);
  long i, l1 = lg(gh), l2 = lg(hg), lm = minss(l1,l2);
  for (i = 1; i < lm; i++)
    if (gh[l1-i] != hg[l2-i]) break;
  return pc_mul(vecsmall_shorten(gh,l1-i), pc_inv(vecsmall_shorten(hg,l2-i), G), G);
}

static GEN
pc_exp(GEN v)
{
  long i, l = lg(v);
  GEN w = cgetg(l, t_VEC);
  if (l==1) return w;
  for (i = 1; i < l; i++)
    gel(w,i) = const_vecsmall(v[i], i+1);
  return shallowconcat1(w);
}
static GEN
vecsmall_increase(GEN x)
{ pari_APPLY_ulong(x[i]+1) }

static GEN
vecvecsmall_increase(GEN x)
{ pari_APPLY_same(vecsmall_increase(gel(x,i))) }

static GEN
pcgrp_lift(GEN G, long deg)
{
  GEN ord = pcgrp_get_ord(G), pw  = pcgrp_get_pow(G), br = pcgrp_get_br(G);
  long i, l = lg(br);
  GEN Ord = vecsmall_prepend(ord, deg);
  GEN Pw = vec_prepend(vecvecsmall_increase(pw), cgetg(1,t_VECSMALL));
  GEN Br = cgetg(l+1, t_VEC);
  gel(Br,1) = const_vec(l-1, cgetg(1, t_VECSMALL));
  for (i = 1; i < l; i++)
    gel(Br,i+1) = vecvecsmall_increase(gel(br, i));
  return mkvec3(Ord, Pw, Br);
}

static GEN
brl_add(GEN x, GEN a)
{
  pari_APPLY_same(vecsmall_concat(const_vecsmall(uel(a,i),1),gel(x,i)))
}

static void
pcgrp_insert(GEN G, long j, GEN a)
{
  GEN pw  = pcgrp_get_pow(G), br = pcgrp_get_br(G);
  gel(pw,j) = vecsmall_concat(gel(a,1),gel(pw, j));
  gel(br,j) = brl_add(gel(br, j), gel(a,2));
}

static long
getfr(GEN f, GEN h)
{
  long i, l = lg(f);
  for (i = 1; i < l; i++)
    if (zv_equal(gel(f,i), h)) return i;
  pari_err_BUG("galoisinit");
  return 0;
}

static long
get_pow(GEN pf, ulong o, GEN pw, GEN gen)
{
  long i, n  = lg(pf)-1;
  GEN p1 = perm_powu(pf, o);
  GEN p2 = pc_to_perm(pw, gen, n);
  for(i = 0; ; i++)
  {
    if (zv_equal(p1, p2)) break;
    p2 = perm_mul(gel(gen,1), p2);
  }
  return i;
}

struct galois_perm
{
  GEN L;
  GEN M;
  GEN den;
  GEN mod, mod2;
  long x;
  GEN cache;
};

static void
galoisperm_init(struct galois_perm *gp, GEN L, GEN M, GEN den, GEN mod, GEN mod2, long x)
{
  gp->L = L;
  gp->M = M;
  gp->den = den;
  gp->mod = mod;
  gp->mod2 = mod2;
  gp->x = x;
  gp->cache = zerovec(lg(L)-1);
}

static void
galoisperm_free(struct galois_perm *gp)
{
  long i, l = lg(gp->cache);
  for (i=1; i<l; i++)
    if (!isintzero(gel(gp->cache,i)))
      gunclone(gel(gp->cache,i));
}

static GEN
permtoaut(GEN p, struct galois_perm *gp)
{
  pari_sp av = avma;
  if (isintzero(gel(gp->cache,p[1])))
  {
    GEN pol = permtopol(p, gp->L, gp->M, gp->den, gp->mod, gp->mod2, gp->x);
    gel(gp->cache,p[1]) = gclone(pol);
  }
  set_avma(av);
  return gel(gp->cache,p[1]);
}

static GEN
pc_evalcache(GEN W, GEN u, GEN sp, GEN T, GEN p, struct galois_perm *gp)
{
  GEN v;
  long ns = sp[1];
  if (!isintzero(gel(W,ns))) return gel(W,ns);
  v = RgX_to_FpX(permtoaut(sp, gp), p);
  gel(W,ns) = FpX_FpXQV_eval(v, u, T, p);
  return gel(W,ns);
}

static ulong
findp(GEN D, GEN P, GEN S, long o, GEN *Tmod)
{
  forprime_t iter;
  ulong p;
  long n = degpol(P);
  u_forprime_init(&iter, n*maxss(expu(n)-3, 2), ULONG_MAX);
  while ((p = u_forprime_next(&iter)))
  {
    GEN F, F1, Sp;
    if (smodis(D, p) == 0)
      continue;
    F = gel(Flx_factor(ZX_to_Flx(P, p), p), 1);
    F1 = gel(F,1);
    if (degpol(F1) != o)
      continue;
    Sp = RgX_to_Flx(S, p);
    if (gequal(Flx_rem(Sp, F1, p), Flx_Frobenius(F1, p)))
    {
      *Tmod = FlxV_to_ZXV(F);
      return p;
    }
  }
  return 0;
}

static GEN
nilp_froblift(GEN genG, GEN autH, long j, GEN pcgrp,
  GEN idp, GEN incl, GEN H, struct galois_lift *gl, struct galois_perm *gp)
{
  pari_sp av = avma;
  GEN T = gl->T, p = gl->p, pe = gl->Q;
  ulong pp = itou(p);
  long e   = gl->e;
  GEN pf   = cgetg(lg(gl->L), t_VECSMALL);
  GEN Tp   = ZX_to_Flx(T, pp);
  GEN Hp   = ZX_to_Flx(H, pp);
  GEN ord = pcgrp_get_ord(pcgrp);
  GEN pcp = gel(pcgrp_get_pow(pcgrp),j+1);
  long o  = uel(ord,1);
  GEN ordH = vecslice(ord,2,lg(ord)-1);
  long n = zv_prod(ordH), k = lg(ordH)-1, l = k-j, m = upowuu(o, l), v = varn(T);
  GEN factTp = gel(Flx_factor(Tp, pp), 1);
  long fp = degpol(gel(factTp, 1));
  GEN frobp = Flxq_autpow(Flx_Frobenius(Tp, pp), fp-1, Tp, pp);
  GEN frob = ZpX_ZpXQ_liftroot(T, Flx_to_ZX(frobp), T, p, e);
  if (galoisfrobeniustest(frob, gl, pf))
  {
    GEN pfi = perm_inv(pf);
    long d = get_pow(pfi, uel(ord,j+1), pcp, genG);
    return mkvec3(pfi, mkvec2(const_vecsmall(d,1),zero_zv(l+1)), gel(factTp, 1));
  }
  else
  {
    GEN frobG = FpXQ_powers(frob, usqrt(degpol(T)), T, pe);
    GEN autHp = RgXV_to_FlxV(autH,pp);
    GEN inclp = RgX_to_Flx(incl,pp);
    GEN factHp = gel(Flx_factor(Hp, pp),1);
    long fr = getfr(factHp, idp);
    GEN minHp  = FlxV_minpolymod(autHp, factHp, pp);
    GEN permfact_Hp = factperm(minHp);
    GEN permfact_Gp = FlxV_Flx_gcd(FlxC_Flxq_eval(factHp, inclp, Tp, pp), Tp, pp);
    GEN bezout_Gpe = bezout_lift_fact(T, FlxV_to_ZXV(permfact_Gp), p, e);
    GEN id = gmael(Flx_factor(gel(permfact_Gp, fr),pp),1,1);
    GEN orbgen = genorbit(ordH, permfact_Hp, fr, n, k, j);
    GEN orb = gel(orbgen,1), gen = gel(orbgen,2);
    long nborb = lg(orb)-1;
    GEN A = cgetg(l+1, t_VECSMALL);
    GEN W = zerovec(lg(gl->L)-1);
    GEN U = zeromatcopy(nborb,degpol(T));
    GEN br = pcgrp_get_br(pcgrp), brj = gcopy(gel(br, j+1));
    GEN Ui = cgetg(nborb+1, t_VEC);
    long a, b, i;
    for(a = 0; a < m; a++)
    {
      pari_timer ti;
      pari_sp av2;
      GEN B = pol_0(v);
      long aa = a;
      if (DEBUGLEVEL>=4) timer_start(&ti);
      for(i = 1; i <= l; i++)
      {
        uel(A,i) = aa % o;
        aa /= o;
      }
      gel(br,j+1) = brl_add(brj, A);
      for(b = 1; b <= nborb; b++)
      {
        GEN br = pc_bracket(pc_exp(gel(gen,b)), mkvecsmall(j+1), pcgrp);
        GEN sp = pc_to_perm(br, genG, degpol(T));
        long u = sp[1];
        long s = permprodeval(permfact_Hp, gel(gen,b), fr);
        if (isintzero(gmael(U,u,s)))
        {
          GEN Ub = pc_evalcache(W, frobG, sp, T, pe, gp);
          gmael(U,u,s) = FpXQ_mul(Ub, gel(bezout_Gpe,s), T, pe);
        }
        gel(Ui, b) = gmael(U,u,s);
      }
      av2 = avma;
      for(b = 1; b <= nborb; b++)
        B = FpX_add(B, gel(Ui,b), pe);
      if (DEBUGLEVEL >= 4) timer_printf(&ti,"Testing candidate %ld",a);
      if (galoisfrobeniustest(B, gl, pf))
      {
        GEN pfi = perm_inv(pf);
        long d = get_pow(pfi, uel(ord,j+1), pcp, genG);
        gel(br,j+1) = brj;
        return gerepilecopy(av,mkvec3(pfi,mkvec2(const_vecsmall(d,1),A),id));
      }
      set_avma(av2);
    }
    return gc_NULL(av);
  }
}

static GEN
galoisgenlift_nilp(GEN PG, GEN O, GEN V, GEN T, GEN frob, GEN sigma,
  struct galois_borne *gb, struct galois_frobenius *gf, struct galois_perm *gp)
{
  long j, n = degpol(T), deg = gf->deg;
  ulong p = gf->p;
  GEN L = gp->L, M =  gp->M, den = gp->den;
  GEN S = fixedfieldinclusion(O, gel(V,2));
  GEN incl = vectopol(S, M, den, gb->ladicabs, shifti(gb->ladicabs,-1), varn(T));
  GEN H = gel(V,3);
  GEN PG1 = gmael(PG, 1, 1);
  GEN PG2 = gmael(PG, 1, 2);
  GEN PG3 = gmael(PG, 1, 3);
  GEN PG4 = gmael(PG, 1, 4);
  long lP = lg(PG1);
  GEN PG5 = pcgrp_lift(gmael(PG, 1, 5), deg);
  GEN res = cgetg(6, t_VEC), res1, res2, res3;
  gel(res,1) = res1 = cgetg(lP + 1, t_VEC);
  gel(res,2) = res2 = cgetg(lP + 1, t_VEC);
  gel(res,3) = res3 = cgetg(lP + 1, t_VEC);
  gel(res,4) = vecsmall_prepend(PG4, p);
  gel(res,5) = PG5;
  gel(res1, 1) = frob;
  gel(res2, 1) = ZX_to_Flx(gel(gf->Tmod,1), p);
  gel(res3, 1) = sigma;
  for (j = 1; j < lP; j++)
  {
    struct galois_lift gl;
    GEN Lden = makeLden(L,den,gb);
    GEN pf;
    initlift(T, den, uel(PG4,j), L, Lden, gb, &gl);
    pf = nilp_froblift(vecslice(res1,1,j), PG3, j, PG5, gel(PG2,j), incl, H, &gl, gp);
    if (!pf) return NULL;
    if (DEBUGLEVEL>=2)
      err_printf("found: %ld/%ld: %Ps: %Ps\n", n, j+1, gel(pf,2),gel(pf,1));
    pcgrp_insert(PG5, j+1, gel(pf,2));
    gel(res1, j+1) = gel(pf,1);
    gel(res2, j+1) = gel(pf,3);
    gel(res3, j+1) = gcopy(permtoaut(gel(pf,1), gp));
  }
  if (DEBUGLEVEL >= 4) err_printf("GaloisConj: Fini!\n");
  return res;
}

static GEN
galoisgenlift(GEN PG, GEN Pg, GEN O, GEN L, GEN M, GEN frob,
              struct galois_borne *gb, struct galois_frobenius *gf)
{
  struct galois_test td;
  GEN res, res1;
  GEN PG1 = gel(PG, 1), PG2 = gel(PG, 2);
  long lP = lg(PG1), j, n = lg(L)-1;
  inittest(L, M, gb->bornesol, gb->ladicsol, &td);
  res = cgetg(3, t_VEC);
  gel(res,1) = res1 = cgetg(lP + 1, t_VEC);
  gel(res,2) = vecsmall_prepend(PG2, gf->deg);
  gel(res1, 1) = vecsmall_copy(frob);
  for (j = 1; j < lP; j++)
  {
    GEN pf = galoisgenliftauto(O, gel(PG1, j), gf->psi[Pg[j]], n, &td);
    if (!pf) { freetest(&td); return NULL; }
    gel(res1, j+1) = pf;
  }
  if (DEBUGLEVEL >= 4) err_printf("GaloisConj: Fini!\n");
  freetest(&td);
  return res;
}

static ulong
psi_order(GEN psi, ulong d)
{
  long i, l = lg(psi);
  ulong s = 1;
  for (i=1; i<l; i++)
    s = clcm(s, d/cgcd(uel(psi, i)-1, d));
  return s;
}

static GEN
galoisgen(GEN T, GEN L, GEN M, GEN den, GEN bad, struct galois_borne *gb,
          const struct galois_analysis *ga)
{
  struct galois_test td;
  struct galois_frobenius gf, ogf;
  pari_sp ltop = avma;
  long x, n = degpol(T), is_central;
  ulong po;
  GEN sigma, res, frob, O, PG, V, ofrob = NULL;

  if (!ga->deg) return NULL;
  x = varn(T);
  if (DEBUGLEVEL >= 9) err_printf("GaloisConj: denominator:%Ps\n", den);
  if (n == 12 && ga->ord==3 && !ga->p4)
  { /* A4 is very probable: test it first */
    pari_sp av = avma;
    if (DEBUGLEVEL >= 4) err_printf("GaloisConj: Testing A4 first\n");
    inittest(L, M, gb->bornesol, gb->ladicsol, &td);
    PG = a4galoisgen(&td);
    freetest(&td);
    if (PG) return gerepileupto(ltop, PG);
    set_avma(av);
  }
  if (n == 24 && ga->ord==3 && ga->p4)
  { /* S4 is very probable: test it first */
    pari_sp av = avma;
    struct galois_lift gl;
    if (DEBUGLEVEL >= 4) err_printf("GaloisConj: Testing S4 first\n");
    initlift(T, den, ga->p4, L, makeLden(L,den,gb), gb, &gl);
    PG = s4galoisgen(&gl);
    if (PG) return gerepileupto(ltop, PG);
    set_avma(av);
  }
  if (n == 36 && ga->ord==3 && ga->p4)
  { /* F36 is very probable: test it first */
    pari_sp av = avma;
    struct galois_lift gl;
    if (DEBUGLEVEL >= 4) err_printf("GaloisConj: Testing 3x3:4 first (p=%ld)\n",ga->p4);
    initlift(T, den, ga->p4, L, makeLden(L,den,gb), gb, &gl);
    PG = f36galoisgen(&gl);
    if (PG) return gerepileupto(ltop, PG);
    set_avma(av);
  }
  frob = galoisfindfrobenius(T, L, den, bad, &gf, gb, ga);
  if (!frob) return gc_NULL(ltop);
  po = psi_order(gf.psi, gf.deg);
  if (!(ga->group&ga_easy) && po < (ulong) gf.deg && gf.deg/radicalu(gf.deg)%po == 0)
  {
    is_central = 1;
    if (!bad) bad = gb->dis;
    if (po > 1)
    {
      ofrob = frob; ogf = gf;
      frob = perm_powu(frob, po);
      gf.deg /= po;
    }
  } else is_central = 0;
  sigma = permtopol(frob, L, M, den, gb->ladicabs, shifti(gb->ladicabs,-1), x);
  if (is_central && gf.fp != gf.deg)
  { gf.p = findp(bad, T, sigma, gf.deg, &gf.Tmod); gf.fp = gf.deg;
    gf.psi = const_vecsmall(lg(gf.Tmod)-1, 1);
  }
  if (gf.deg == n)        /* cyclic */
  {
    GEN Tp = ZX_to_Flx(gel(gf.Tmod,1), gf.p);
    res = mkvec5(mkvec(frob), mkvec(Tp), mkvec(sigma), mkvecsmall(gf.p), cyclic_pc(n));
    return gerepilecopy(ltop, res);
  }
  O = perm_cycles(frob);
  if (DEBUGLEVEL >= 9) err_printf("GaloisConj: Frobenius:%Ps\n", sigma);
  PG = galoisgenfixedfield0(O, L, sigma, T, is_central ? bad: NULL,
                                            is_central ? &V:  NULL, &gf, gb);
  if (PG == NULL) return gc_NULL(ltop);
  if (is_central && lg(gel(PG,1))!=3)
  {
    struct galois_perm gp;
    galoisperm_init(&gp, L, M, den, gb->ladicabs, shifti(gb->ladicabs,-1), varn(T));
    res = galoisgenlift_nilp(PG, O, V, T, frob, sigma, gb, &gf, &gp);
    galoisperm_free(&gp);
  }
  else
  {
    if (is_central && po > 1)
    { /* backtrack powering of frob */
      frob = ofrob; gf = ogf;
      O = perm_cycles(ofrob);
      sigma = permtopol(ofrob, L, M, den, gb->ladicabs, shifti(gb->ladicabs,-1), x);
      PG = galoisgenfixedfield0(O, L, sigma, T, NULL, NULL, &gf, gb);
      if (PG == NULL) return gc_NULL(ltop);
    }
    res = galoisgenlift(gg_get_std(gel(PG,1)), gel(PG,2), O, L, M, frob, gb, &gf);
  }
  if (!res) return gc_NULL(ltop);
  return gerepilecopy(ltop, res);
}

/* T = polcyclo(N) */
static GEN
conjcyclo(GEN T, long N)
{
  pari_sp av = avma;
  long i, k = 1, d = eulerphiu(N), v = varn(T);
  GEN L = cgetg(d+1,t_COL);
  for (i=1; i<=N; i++)
    if (ugcd(i, N)==1)
    {
      GEN s = pol_xn(i, v);
      if (i >= d) s = ZX_rem(s, T);
      gel(L,k++) = s;
    }
  return gerepileupto(av, gen_sort(L, (void*)&gcmp, &gen_cmp_RgX));
}

static GEN
aut_to_groupelts(GEN aut, GEN L, ulong p)
{
  pari_sp av = avma;
  long i, d = lg(aut)-1;
  GEN P = ZV_to_Flv(L, p);
  GEN N = FlxV_Flv_multieval(aut, P, p);
  GEN q = perm_inv(vecsmall_indexsort(P));
  GEN G = cgetg(d+1, t_VEC);
  for (i=1; i<=d; i++)
    gel(G,i) = perm_mul(vecsmall_indexsort(gel(N,i)), q);
  return gerepilecopy(av, vecvecsmall_sort_shallow(G));
}

static ulong
galois_find_totally_split(GEN P, GEN Q)
{
  pari_sp av = avma;
  forprime_t iter;
  ulong p;
  long n = degpol(P);
  u_forprime_init(&iter, n*maxss(expu(n)-3, 2), ULONG_MAX);
  while ((p = u_forprime_next(&iter)))
  {
    if (Flx_is_totally_split(ZX_to_Flx(P, p), p)
       && (!Q || Flx_is_squarefree(ZX_to_Flx(Q, p), p)))
      return gc_ulong(av, p);
    set_avma(av);
  }
  return 0;
}

GEN
galoisinitfromaut(GEN T, GEN aut, ulong l)
{
  pari_sp ltop = avma;
  GEN nf, A, G, L, M, grp, den=NULL;
  struct galois_analysis ga;
  struct galois_borne gb;
  long n;
  pari_timer ti;

  T = get_nfpol(T, &nf);
  n = degpol(T);
  if (nf)
  { if (!den) den = nf_get_zkden(nf); }
  else
  {
    if (n <= 0) pari_err_IRREDPOL("galoisinit",T);
    RgX_check_ZX(T, "galoisinit");
    if (!ZX_is_squarefree(T))
      pari_err_DOMAIN("galoisinit","issquarefree(pol)","=",gen_0,T);
    if (!gequal1(gel(T,n+2))) pari_err_IMPL("galoisinit(nonmonic)");
  }
  if (lg(aut)-1 != n)
    return gen_0;
  ga.l = l? l: galois_find_totally_split(T, NULL);
  if (!l) aut = RgXV_to_FlxV(aut, ga.l);
  gb.l = utoipos(ga.l);
  if (DEBUGLEVEL >= 1) timer_start(&ti);
  den = galoisborne(T, den, &gb, degpol(T));
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "galoisborne()");
  L = ZpX_roots(T, gb.l, gb.valabs);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "ZpX_roots");
  M = FpV_invVandermonde(L, den, gb.ladicabs);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "FpV_invVandermonde()");
  A = aut_to_groupelts(aut, L, ga.l);
  G = groupelts_to_group(A);
  if (!G) G = trivialgroup();
  else A = group_elts(G,n);
  grp = cgetg(9, t_VEC);
  gel(grp,1) = T;
  gel(grp,2) = mkvec3(utoipos(ga.l), utoipos(gb.valabs), gb.ladicabs);
  gel(grp,3) = L;
  gel(grp,4) = M;
  gel(grp,5) = den;
  gel(grp,6) = A;
  gel(grp,7) = gel(G,1);
  gel(grp,8) = gel(G,2);
  return gerepilecopy(ltop, grp);
}

GEN
galoissplittinginit(GEN T, GEN D)
{
  pari_sp av = avma;
  GEN R = nfsplitting0(T, D, 3), P = gel(R,1), aut = gel(R,2);
  ulong p = itou(gel(R,3));
  return gerepileupto(av, galoisinitfromaut(P, aut, p));
}

/* T: polynomial or nf, den multiple of common denominator of solutions or
 * NULL (unknown). If T is nf, and den unknown, use den = denom(nf.zk) */
static GEN
galoisconj4_main(GEN T, GEN den, long flag)
{
  pari_sp ltop = avma;
  GEN nf, G, L, M, aut, grp;
  struct galois_analysis ga;
  struct galois_borne gb;
  long n;
  pari_timer ti;

  T = get_nfpol(T, &nf);
  n = poliscyclo(T);
  if (n) return flag? galoiscyclo(n, varn(T)): conjcyclo(T, n);
  n = degpol(T);
  if (nf)
  { if (!den) den = nf_get_zkden(nf); }
  else
  {
    if (n <= 0) pari_err_IRREDPOL("galoisinit",T);
    RgX_check_ZX(T, "galoisinit");
    if (!ZX_is_squarefree(T))
      pari_err_DOMAIN("galoisinit","issquarefree(pol)","=",gen_0,T);
    if (!gequal1(gel(T,n+2))) pari_err_IMPL("galoisinit(nonmonic)");
  }
  if (n == 1)
  {
    if (!flag) { G = cgetg(2, t_COL); gel(G,1) = pol_x(varn(T)); return G;}
    ga.l = 3;
    ga.deg = 1;
    den = gen_1;
  }
  else if (!galoisanalysis(T, &ga, 1, NULL)) return gc_NULL(ltop);

  if (den)
  {
    if (typ(den) != t_INT) pari_err_TYPE("galoisconj4 [2nd arg integer]", den);
    den = absi_shallow(den);
  }
  gb.l = utoipos(ga.l);
  if (DEBUGLEVEL >= 1) timer_start(&ti);
  den = galoisborne(T, den, &gb, degpol(T));
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "galoisborne()");
  L = ZpX_roots(T, gb.l, gb.valabs);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "ZpX_roots");
  M = FpV_invVandermonde(L, den, gb.ladicabs);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "FpV_invVandermonde()");
  if (n == 1)
  {
    G = cgetg(3, t_VEC);
    gel(G,1) = cgetg(1, t_VEC);
    gel(G,2) = cgetg(1, t_VECSMALL);
  }
  else
    G = gg_get_std(galoisgen(T, L, M, den, NULL, &gb, &ga));
  if (DEBUGLEVEL >= 6) err_printf("GaloisConj: %Ps\n", G);
  if (!G) return gc_NULL(ltop);
  if (DEBUGLEVEL >= 1) timer_start(&ti);
  grp = cgetg(9, t_VEC);
  gel(grp,1) = T;
  gel(grp,2) = mkvec3(utoipos(ga.l), utoipos(gb.valabs), gb.ladicabs);
  gel(grp,3) = L;
  gel(grp,4) = M;
  gel(grp,5) = den;
  gel(grp,6) = group_elts(G,n);
  gel(grp,7) = gel(G,1);
  gel(grp,8) = gel(G,2);
  if (flag) return gerepilecopy(ltop, grp);
  aut = galoisvecpermtopol(grp, gal_get_group(grp), gb.ladicabs, shifti(gb.ladicabs,-1));
  settyp(aut, t_COL);
  if (DEBUGLEVEL >= 1) timer_printf(&ti, "Computation of polynomials");
  return gerepileupto(ltop, gen_sort(aut, (void*)&gcmp, &gen_cmp_RgX));
}

/* Heuristic computation of #Aut(T), pinit = first prime to be tested */
long
numberofconjugates(GEN T, long pinit)
{
  pari_sp av = avma;
  long c, nbtest, nbmax, n = degpol(T);
  ulong p;
  forprime_t S;

  if (n == 1) return 1;
  nbmax = (n < 10)? 20: (n<<1) + 1;
  nbtest = 0;
#if 0
  c = ZX_sturm(T); c = ugcd(c, n-c); /* too costly: finite primes are cheaper */
#else
  c = n;
#endif
  u_forprime_init(&S, pinit, ULONG_MAX);
  while((p = u_forprime_next(&S)))
  {
    GEN L, Tp = ZX_to_Flx(T,p);
    long i, nb;
    if (!Flx_is_squarefree(Tp, p)) continue;
    /* unramified */
    nbtest++;
    L = Flx_nbfact_by_degree(Tp, &nb, p); /* L[i] = #factors of degree i */
    if (L[n/nb] == nb) {
      if (c == n && nbtest > 10) break; /* probably Galois */
    }
    else
    {
      c = ugcd(c, L[1]);
      for (i = 2; i <= n; i++)
        if (L[i]) { c = ugcd(c, L[i]*i); if (c == 1) break; }
      if (c == 1) break;
    }
    if (nbtest == nbmax) break;
    if (DEBUGLEVEL >= 6)
      err_printf("NumberOfConjugates [%ld]:c=%ld,p=%ld\n", nbtest,c,p);
    set_avma(av);
  }
  if (DEBUGLEVEL >= 2) err_printf("NumberOfConjugates:c=%ld,p=%ld\n", c, p);
  return gc_long(av,c);
}
static GEN
galoisconj4(GEN nf, GEN d)
{
  pari_sp av = avma;
  GEN G, T;
  G = galoisconj4_main(nf, d, 0);
  if (G) return G; /* Success */
  set_avma(av); T = get_nfpol(nf, &nf);
  G = cgetg(2, t_COL); gel(G,1) = pol_x(varn(T)); return G; /* Fail */

}

/* d multiplicative bound for the automorphism's denominators */
static GEN
galoisconj_monic(GEN nf, GEN d)
{
  pari_sp av = avma;
  GEN G, NF, T = get_nfpol(nf,&NF);
  if (degpol(T) == 2)
  { /* fast shortcut */
    GEN b = gel(T,3);
    long v = varn(T);
    G = cgetg(3, t_COL);
    gel(G,1) = deg1pol_shallow(gen_m1, negi(b), v);
    gel(G,2) = pol_x(v);
    return G;
  }
  G = galoisconj4_main(nf, d, 0);
  if (G) return G; /* Success */
  set_avma(av); return galoisconj1(nf);
}

GEN
galoisconj(GEN nf, GEN d)
{
  pari_sp av;
  GEN NF, S, L, T = get_nfpol(nf,&NF);
  if (NF) return galoisconj_monic(NF, d);
  RgX_check_QX(T, "galoisconj");
  av = avma;
  T = Q_primpart(T);
  if (ZX_is_monic(T)) return galoisconj_monic(T, d);
  S = galoisconj_monic(poltomonic(T,&L), NULL);
  return gerepileupto(av, gdiv(RgXV_unscale(S, L),L));
}

/* FIXME: obsolete, use galoisconj(nf, d) directly */
GEN
galoisconj0(GEN nf, long flag, GEN d, long prec)
{
  (void)prec;
  switch(flag) {
    case 2:
    case 0: return galoisconj(nf, d);
    case 1: return galoisconj1(nf);
    case 4: return galoisconj4(nf, d);
  }
  pari_err_FLAG("nfgaloisconj");
  return NULL; /*LCOV_EXCL_LINE*/
}

/******************************************************************************/
/* Galois theory related algorithms                                           */
/******************************************************************************/
GEN
checkgal(GEN gal)
{
  if (typ(gal) == t_POL) pari_err_TYPE("checkgal [apply galoisinit first]",gal);
  if (typ(gal) != t_VEC || lg(gal) != 9) pari_err_TYPE("checkgal",gal);
  return gal;
}

GEN
galoisinit(GEN nf, GEN den)
{
  GEN G;
  if (is_vec_t(typ(nf)) && lg(nf)==3 && is_vec_t(typ(gel(nf,2))))
    return galoisinitfromaut(gel(nf,1), gel(nf,2), 0);
  G = galoisconj4_main(nf, den, 1);
  return G? G: gen_0;
}

static GEN
galoispermtopol_i(GEN gal, GEN perm, GEN mod, GEN mod2)
{
  switch (typ(perm))
  {
    case t_VECSMALL:
      return permtopol(perm, gal_get_roots(gal), gal_get_invvdm(gal),
                             gal_get_den(gal), mod, mod2,
                             varn(gal_get_pol(gal)));
    case t_VEC: case t_COL: case t_MAT:
      return galoisvecpermtopol(gal, perm, mod, mod2);
  }
  pari_err_TYPE("galoispermtopol", perm);
  return NULL; /* LCOV_EXCL_LINE */
}

GEN
galoispermtopol(GEN gal, GEN perm)
{
  pari_sp av = avma;
  GEN mod, mod2;
  gal = checkgal(gal);
  mod = gal_get_mod(gal);
  mod2 = shifti(mod,-1);
  return gerepilecopy(av, galoispermtopol_i(gal, perm, mod, mod2));
}

GEN
galoiscosets(GEN O, GEN perm)
{
  long i, j, k, u, f, l = lg(O);
  GEN RC, C = cgetg(l,t_VECSMALL), o = gel(O,1);
  pari_sp av = avma;
  f = lg(o); u = o[1]; RC = zero_zv(lg(perm)-1);
  for(i=1,j=1; j<l; i++)
  {
    GEN p = gel(perm,i);
    if (RC[ p[u] ]) continue;
    for(k=1; k<f; k++) RC[ p[ o[k] ] ] = 1;
    C[j++] = i;
  }
  set_avma(av); return C;
}

static GEN
fixedfieldfactor(GEN L, GEN O, GEN perm, GEN M, GEN den, GEN mod, GEN mod2,
                 long x,long y)
{
  pari_sp ltop = avma;
  long i, j, k, l = lg(O), lo = lg(gel(O,1));
  GEN V, res, cosets = galoiscosets(O,perm), F = cgetg(lo+1,t_COL);

  gel(F, lo) = gen_1;
  if (DEBUGLEVEL>=4) err_printf("GaloisFixedField:cosets=%Ps \n",cosets);
  if (DEBUGLEVEL>=6) err_printf("GaloisFixedField:den=%Ps mod=%Ps \n",den,mod);
  V = cgetg(l,t_COL); res = cgetg(l,t_VEC);
  for (i = 1; i < l; i++)
  {
    pari_sp av = avma;
    GEN G = cgetg(l,t_VEC), Lp = vecpermute(L, gel(perm, cosets[i]));
    for (k = 1; k < l; k++)
      gel(G,k) = FpV_roots_to_pol(vecpermute(Lp, gel(O,k)), mod, x);
    for (j = 1; j < lo; j++)
    {
      for(k = 1; k < l; k++) gel(V,k) = gmael(G,k,j+1);
      gel(F,j) = vectopol(V, M, den, mod, mod2, y);
    }
    gel(res,i) = gerepileupto(av,gtopolyrev(F,x));
  }
  return gerepileupto(ltop,res);
}

static void
chk_perm(GEN perm, long n)
{
  if (typ(perm) != t_VECSMALL || lg(perm)!=n+1)
    pari_err_TYPE("galoisfixedfield", perm);
}

static int
is_group(GEN g)
{
  if (typ(g) == t_VEC && lg(g) == 3)
  {
    GEN a = gel(g,1), o = gel(g,2);
    return typ(a)==t_VEC && typ(o)==t_VECSMALL && lg(a) == lg(o);
  }
  return 0;
}

GEN
galoisfixedfield(GEN gal, GEN perm, long flag, long y)
{
  pari_sp ltop = avma;
  GEN T, L, P, S, PL, O, res, mod, mod2, OL, sym;
  long vT, n, i;
  if (flag<0 || flag>2) pari_err_FLAG("galoisfixedfield");
  gal = checkgal(gal); T = gal_get_pol(gal);
  vT = varn(T);
  L = gal_get_roots(gal); n = lg(L)-1;
  mod = gal_get_mod(gal);
  if (typ(perm) == t_VEC)
  {
    if (is_group(perm)) perm = gel(perm, 1);
    for (i = 1; i < lg(perm); i++) chk_perm(gel(perm,i), n);
    O = vecperm_orbits(perm, n);
  }
  else
  {
    chk_perm(perm, n);
    O = perm_cycles(perm);
  }
  mod2 = shifti(mod,-1);
  OL = fixedfieldorbits(O, L);
  sym = fixedfieldsympol(OL, itou(gal_get_p(gal)));
  PL = sympol_eval(sym, OL, mod);
  P = FpX_center_i(FpV_roots_to_pol(PL, mod, vT), mod, mod2);
  if (flag==1) return gerepilecopy(ltop,P);
  S = fixedfieldinclusion(O, PL);
  S = vectopol(S, gal_get_invvdm(gal), gal_get_den(gal), mod, mod2, vT);
  if (flag==0)
    res = cgetg(3, t_VEC);
  else
  {
    GEN PM, Pden;
    struct galois_borne Pgb;
    long val = itos(gal_get_e(gal));
    Pgb.l = gal_get_p(gal);
    Pden = galoisborne(P, NULL, &Pgb, degpol(T)/degpol(P));
    if (Pgb.valabs > val)
    {
      if (DEBUGLEVEL>=4)
        err_printf("GaloisConj: increase p-adic prec by %ld.\n", Pgb.valabs-val);
      PL = ZpX_liftroots(P, PL, Pgb.l, Pgb.valabs);
      L  = ZpX_liftroots(T, L, Pgb.l, Pgb.valabs);
      mod = Pgb.ladicabs; mod2 = shifti(mod,-1);
    }
    PM = FpV_invVandermonde(PL, Pden, mod);
    if (y < 0) y = 1;
    if (varncmp(y, vT) <= 0)
      pari_err_PRIORITY("galoisfixedfield", T, "<=", y);
    setvarn(P, y);
    res = cgetg(4, t_VEC);
    gel(res,3) = fixedfieldfactor(L,O,gal_get_group(gal), PM,Pden,mod,mod2,vT,y);
  }
  gel(res,1) = gcopy(P);
  gel(res,2) = gmodulo(S, T);
  return gerepileupto(ltop, res);
}

/* gal a galois group output the underlying wss group */
GEN
galois_group(GEN gal) { return mkvec2(gal_get_gen(gal), gal_get_orders(gal)); }

GEN
checkgroup(GEN g, GEN *S)
{
  if (is_group(g)) { *S = NULL; return g; }
  g  = checkgal(g);
  *S = gal_get_group(g); return galois_group(g);
}

GEN
checkgroupelts(GEN G)
{
  long i, n;
  if (typ(G)!=t_VEC) pari_err_TYPE("checkgroupelts", G);
  if (is_group(G))
  { /* subgroup of S_n */
    if (lg(gel(G,1))==1) return mkvec(mkvecsmall(1));
    return group_elts(G, group_domain(G));
  }
  if (lg(G)==9 && typ(gel(G,1))==t_POL)
    return gal_get_group(G); /* galoisinit */
  /* vector of permutations ? */
  n = lg(G)-1;
  if (n==0) pari_err_DIM("checkgroupelts");
  for (i = 1; i <= n; i++)
  {
    if (typ(gel(G,i)) != t_VECSMALL)
      pari_err_TYPE("checkgroupelts (element)", gel(G,i));
    if (lg(gel(G,i)) != lg(gel(G,1)))
      pari_err_DIM("checkgroupelts [length of permutations]");
  }
  return G;
}

GEN
galoisisabelian(GEN gal, long flag)
{
  pari_sp av = avma;
  GEN S, G = checkgroup(gal,&S);
  if (!group_isabelian(G)) { set_avma(av); return gen_0; }
  switch(flag)
  {
    case 0: return gerepileupto(av, group_abelianHNF(G,S));
    case 1: set_avma(av); return gen_1;
    case 2: return gerepileupto(av, group_abelianSNF(G,S));
    default: pari_err_FLAG("galoisisabelian");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

long
galoisisnormal(GEN gal, GEN sub)
{
  pari_sp av = avma;
  GEN S, G = checkgroup(gal, &S), H = checkgroup(sub, &S);
  long res = group_subgroup_isnormal(G, H);
  set_avma(av);
  return res;
}

static GEN
conjclasses_count(GEN conj, long nb)
{
  long i, l = lg(conj);
  GEN c = zero_zv(nb);
  for (i = 1; i < l; i++) c[conj[i]]++;
  return c;
}
GEN
galoisconjclasses(GEN G)
{
  pari_sp av = avma;
  GEN c, e, cc = group_to_cc(G);
  GEN elts = gel(cc,1), conj = gel(cc,2), repr = gel(cc,3);
  long i, l = lg(conj), lc = lg(repr);
  c = conjclasses_count(conj, lc-1);
  e = cgetg(lc, t_VEC);
  for (i = 1; i < lc; i++) gel(e,i) = cgetg(c[i]+1, t_VEC);
  for (i = 1; i < l; i++)
  {
    long ci = conj[i];
    gmael(e, ci, c[ci]) = gel(elts, i);
    c[ci]--;
  }
  return gerepilecopy(av, e);
}

static GEN
groupelts_to_group_or_elts(GEN elts)
{
  GEN G = groupelts_to_group(elts);
  return G ? G: gcopy(elts);
}

static GEN
vec_groupelts_to_group_or_elts(GEN x)
{ pari_APPLY_same(groupelts_to_group_or_elts(gel(x,i))) }

GEN
galoissubgroups(GEN gal)
{
  pari_sp av = avma;
  GEN S, G = checkgroup(gal,&S);
  if (lg(gel(G,1))==1 && lg(S)>2)
    return gerepileupto(av,
      vec_groupelts_to_group_or_elts(groupelts_solvablesubgroups(S)));
  return gerepileupto(av, group_subgroups(G));
}

GEN
galoissubfields(GEN G, long flag, long v)
{
  pari_sp av = avma;
  GEN L = galoissubgroups(G);
  long i, l = lg(L);
  GEN S = cgetg(l, t_VEC);
  for (i = 1; i < l; ++i) gel(S,i) = galoisfixedfield(G, gmael(L,i,1), flag, v);
  return gerepileupto(av, S);
}

GEN
galoisexport(GEN gal, long format)
{
  pari_sp av = avma;
  GEN S, G = checkgroup(gal,&S);
  return gerepileupto(av, group_export(G,format));
}

GEN
galoisidentify(GEN gal)
{
  pari_sp av = avma;
  GEN S, G = checkgroup(gal,&S);
  long idx = group_ident(G,S), card = S ? lg(S)-1: group_order(G);
  set_avma(av); return mkvec2s(card, idx);
}

/* index of conjugacy class containing g */
static long
cc_id(GEN cc, GEN g)
{
  GEN conj = gel(cc,2);
  long k = signe(gel(cc,4))? g[1]: vecvecsmall_search(gel(cc,1), g);
  return conj[k];
}

static GEN
Qevproj_RgX(GEN c, long d, GEN pro)
{ return RgV_to_RgX(Qevproj_down(RgX_to_RgC(c,d), pro), varn(c)); }
/* c in Z[X] / (X^o-1), To = polcyclo(o), T = polcyclo(expo), e = expo/o
 * return c(X^e) mod T as an element of Z[X] / (To) */
static GEN
chival(GEN c, GEN T, GEN To, long e, GEN pro, long phie)
{
  c = ZX_rem(c, To);
  if (e != 1) c = ZX_rem(RgX_inflate(c,e), T);
  if (pro) c = Qevproj_RgX(c, phie, pro);
  return c;
}
/* chi(g^l) = sum_{k=0}^{o-1} a_k zeta_o^{l*k} for all l;
* => a_k = 1/o sum_{l=0}^{o-1} chi(g^l) zeta_o^{-k*l}. Assume o > 1 */
static GEN
chiFT(GEN cp, GEN jg, GEN vze, long e, long o, ulong p, ulong pov2)
{
  const long var = 1;
  ulong oinv = Fl_inv(o,p);
  long k, l;
  GEN c = cgetg(o+2, t_POL);
  for (k = 0; k < o; k++)
  {
    ulong a = 0;
    for (l=0; l<o; l++)
    {
      ulong z = vze[Fl_mul(k,l,o)*e + 1];/* zeta_o^{-k*l} */
      a = Fl_add(a, Fl_mul(uel(cp,jg[l+1]), z, p), p);
    }
    gel(c,k+2) = stoi(Fl_center(Fl_mul(a,oinv,p), p, pov2)); /* a_k */
  }
  c[1] = evalvarn(var) | evalsigne(1); return ZX_renormalize(c,o+2);
}
static GEN
cc_chartable(GEN cc)
{
  GEN al, elts, rep, ctp, ct, dec, id, vjg, H, vord, operm;
  long i, j, k, f, l, expo, lcl, n;
  ulong p, pov2;

  elts = gel(cc,1); n = lg(elts)-1;
  if (n == 1) return mkvec2(mkmat(mkcol(gen_1)), gen_1);
  rep = gel(cc,3);
  lcl = lg(rep);
  vjg = cgetg(lcl, t_VEC);
  vord = cgetg(lcl,t_VECSMALL);
  id = identity_perm(lg(gel(elts,1))-1);
  expo = 1;
  for(j=1;j<lcl;j++)
  {
    GEN jg, h = id, g = gel(elts,rep[j]);
    long o;
    vord[j] = o = perm_orderu(g);
    expo = ulcm(expo, o);
    gel(vjg,j) = jg = cgetg(o+1,t_VECSMALL);
    for (l=1; l<=o; l++)
    {
      jg[l] = cc_id(cc, h); /* index of conjugacy class of g^(l-1) */
      if (l < o) h = perm_mul(h, g);
    }
  }
  /* would sort conjugacy classes by inc. order */
  operm = vecsmall_indexsort(vord);

  /* expo > 1, exponent of G */
  p = unextprime(2*n+1);
  while (p%expo != 1) p = unextprime(p+1);
  /* compute character table modulo p: idempotents of Z(KG) */
  al = conjclasses_algcenter(cc, utoipos(p));
  dec = algsimpledec_ss(al,1);
  ctp = cgetg(lcl,t_VEC);
  for(i=1; i<lcl; i++)
  {
    GEN e = ZV_to_Flv(gmael3(dec,i,3,1), p); /*(1/n)[(dim chi)chi(g): g in G]*/
    ulong d = usqrt(Fl_mul(e[1], n, p)); /* = chi(1) <= sqrt(n) < sqrt(p) */
    gel(ctp,i) = Flv_Fl_mul(e,Fl_div(n,d,p), p); /*[chi(g): g in G]*/
  }
  /* Find minimal f such that table is defined over Q(zeta(f)): the conductor
   * of the class field Q(\zeta_e)^H defined by subgroup
   * H = { k in (Z/e)^*: g^k ~ g, for all g } */
  H = coprimes_zv(expo);
  for (k = 2; k < expo; k++)
  {
    if (!H[k]) continue;
    for (j = 2; j < lcl; j++) /* skip g ~ 1 */
      if (umael(vjg,j,(k % vord[j])+1) != umael(vjg,j,2)) { H[k] = 0; break; }
  }
  f = znstar_conductor_bits(Flv_to_F2v(H));
  /* lift character table to Z[zeta_f] */
  pov2 = p>>1;
  ct = cgetg(lcl, t_MAT);
  if (f == 1)
  { /* rational representation */
    for (j=1; j<lcl; j++) gel(ct,j) = cgetg(lcl,t_COL);
    for(j=1; j<lcl; j++)
    {
      GEN jg = gel(vjg,j); /* jg[l+1] = class of g^l */
      long t = lg(jg) > 2? jg[2]: jg[1];
      for(i=1; i<lcl; i++)
      {
        GEN cp = gel(ctp,i); /* cp[i] = chi(g_i) mod \P */
        gcoeff(ct,j,i) = stoi(Fl_center(cp[t], p, pov2));
      }
    }
  }
  else
  {
    const long var = 1;
    ulong ze = Fl_powu(pgener_Fl(p),(p-1)/expo, p); /* seen as zeta_e^(-1) */
    GEN vze = Fl_powers(ze, expo-1, p); /* vze[i] = ze^(i-1) */
    GEN vzeZX = const_vec(p, gen_0);
    GEN T = polcyclo(expo, var), vT = const_vec(expo,NULL), pro = NULL;
    long phie = degpol(T), id1 = gel(vjg,1)[1]; /* index of 1_G, always 1 ? */
    gel(vT, expo) = T;
    if (f != expo)
    {
      long phif = eulerphiu(f);
      GEN zf = ZX_rem(pol_xn(expo/f,var), T), zfj = zf;
      GEN M = cgetg(phif+1, t_MAT);
      gel(M,1) = col_ei(phie,1);
      for (j = 2; j <= phif; j++)
      {
        gel(M,j) = RgX_to_RgC(zfj, phie);
        if (j < phif) zfj = ZX_rem(ZX_mul(zfj, zf), T);
      }
      pro = Qevproj_init(M);
    }
    gel(vzeZX,1) = pol_1(var);
    for (i = 2; i <= expo; i++)
    {
      GEN t = ZX_rem(pol_xn(expo-(i-1), var), T);
      if (pro) t = Qevproj_RgX(t, phie, pro);
      gel(vzeZX, vze[i]) = t;
    }
    for(i=1; i<lcl; i++)
    { /* loop over characters */
      GEN cp = gel(ctp,i), C, cj; /* cp[j] = chi(g_j) mod \P */
      long dim = cp[id1];
      gel(ct, i) = C = const_col(lcl-1, NULL);
      gel(C,operm[1]) = utoi(dim); /* chi(1_G) */
      for (j=lcl-1; j > 1; j--)
      { /* loop over conjugacy classes, decreasing order: skip 1_G */
        long e, jperm = operm[j], o = vord[jperm];
        GEN To, jg = gel(vjg,jperm); /* jg[l+1] = class of g^l */

        if (gel(C, jperm)) continue; /* done already */
        if (dim == 1) { gel(C, jperm) = gel(vzeZX, cp[jg[2]]); continue; }
        e = expo / o;
        cj = chiFT(cp, jg, vze, e, o, p, pov2);
        To = gel(vT, o); if (!To) To = gel(vT,o) = polcyclo(o, var);
        gel(C, jperm) = chival(cj, T, To, e, pro, phie);
        for (k = 2; k < o; k++)
        {
          GEN ck = RgX_inflate(cj, k); /* chi(g^k) */
          gel(C, jg[k+1]) = chival(ck, T, To, e, pro, phie);
        }
      }
    }
  }
  ct = gen_sort_shallow(ct,(void*)cmp_universal,cmp_nodata);
  i = 1; while (!vec_isconst(gel(ct,i))) i++;
  if (i > 1) swap(gel(ct,1), gel(ct,i));
  return mkvec2(ct, utoipos(f));
}
GEN
galoischartable(GEN gal)
{
  pari_sp av = avma;
  GEN cc = group_to_cc(gal);
  return gerepilecopy(av, cc_chartable(cc));
}

static void
checkgaloischar(GEN ch, GEN repr)
{
  if (gvar(ch) == 0) pari_err_PRIORITY("galoischarpoly",ch,"=",0);
  if (!is_vec_t(typ(ch))) pari_err_TYPE("galoischarpoly", ch);
  if (lg(repr) != lg(ch)) pari_err_DIM("galoischarpoly");
}

static long
galoischar_dim(GEN ch)
{
  pari_sp av = avma;
  long d = gtos(simplify_shallow(lift_shallow(gel(ch,1))));
  return gc_long(av,d);
}

static GEN
galoischar_aut_charpoly(GEN cc, GEN ch, GEN p, long d)
{
  GEN q = p, V = cgetg(d+2, t_POL);
  long i;
  V[1] = evalsigne(1)|evalvarn(0);
  for (i = 1; i <= d; i++)
  {
    gel(V,i+1) = gel(ch, cc_id(cc,q));
    if (i < d) q = perm_mul(q, p);
  }
  return liftpol_shallow(RgXn_expint(RgX_neg(V),d+1));
}

static GEN
galoischar_charpoly(GEN cc, GEN ch, long o)
{
  GEN chm, V, elts = gel(cc,1), repr = gel(cc,3);
  long i, d, l = lg(ch), v = gvar(ch);
  checkgaloischar(ch, repr);
  chm = v < 0 ? ch: gmodulo(ch, polcyclo(o, v));
  V = cgetg(l, t_COL); d = galoischar_dim(ch);
  for (i = 1; i < l; i++)
    gel(V,i) = galoischar_aut_charpoly(cc, chm, gel(elts,repr[i]), d);
  return V;
}

GEN
galoischarpoly(GEN gal, GEN ch, long o)
{
  pari_sp av = avma;
  GEN cc = group_to_cc(gal);
  return gerepilecopy(av, galoischar_charpoly(cc, ch, o));
}

static GEN
cc_char_det(GEN cc, GEN ch, long o)
{
  long i, l = lg(ch), d = galoischar_dim(ch);
  GEN V = galoischar_charpoly(cc, ch, o);
  for (i = 1; i < l; i++) gel(V,i) = leading_coeff(gel(V,i));
  return odd(d)? gneg(V): V;
}

GEN
galoischardet(GEN gal, GEN ch, long o)
{
  pari_sp av = avma;
  GEN cc = group_to_cc(gal);
  return gerepilecopy(av, cc_char_det(cc, ch, o));
}
