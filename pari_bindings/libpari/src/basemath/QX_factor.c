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

#define DEBUGLEVEL DEBUGLEVEL_factor

/* x,y two ZX, y non constant. Return q = x/y if y divides x in Z[X] and NULL
 * otherwise. If not NULL, B is a t_INT upper bound for ||q||_oo. */
static GEN
ZX_divides_i(GEN x, GEN y, GEN B)
{
  long dx, dy, dz, i, j;
  pari_sp av;
  GEN z,p1,y_lead;

  dy=degpol(y);
  dx=degpol(x);
  dz=dx-dy; if (dz<0) return NULL;
  z=cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;
  y_lead = gel(y,dy);
  if (equali1(y_lead)) y_lead = NULL;

  p1 = gel(x,dx);
  if (y_lead) {
    GEN r;
    p1 = dvmdii(p1,y_lead, &r);
    if (r != gen_0) return NULL;
  }
  else p1 = icopy(p1);
  gel(z,dz) = p1;
  for (i=dx-1; i>=dy; i--)
  {
    av = avma; p1 = gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    if (y_lead) {
      GEN r;
      p1 = dvmdii(p1,y_lead, &r);
      if (r != gen_0) return NULL;
    }
    if (B && abscmpii(p1, B) > 0) return NULL;
    p1 = gerepileuptoint(av, p1);
    gel(z,i-dy) = p1;
  }
  av = avma;
  for (; i >= 0; i--)
  {
    p1 = gel(x,i);
    /* we always enter this loop at least once */
    for (j=0; j<=i && j<=dz; j++)
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    if (signe(p1)) return NULL;
    set_avma(av);
  }
  return z - 2;
}
static GEN
ZX_divides(GEN x, GEN y) { return ZX_divides_i(x,y,NULL); }

#if 0
/* cf Beauzamy et al: upper bound for
 *      lc(x) * [2^(5/8) / pi^(3/8)] e^(1/4n) 2^(n/2) sqrt([x]_2)/ n^(3/8)
 * where [x]_2 = sqrt(\sum_i=0^n x[i]^2 / binomial(n,i)). One factor has
 * all coeffs less than then bound */
static GEN
two_factor_bound(GEN x)
{
  long i, j, n = lg(x) - 3;
  pari_sp av = avma;
  GEN *invbin, c, r = cgetr(3), z;

  x += 2; invbin = (GEN*)new_chunk(n+1);
  z = real_1(LOWDEFAULTPREC); /* invbin[i] = 1 / binomial(n, i) */
  for (i=0,j=n; j >= i; i++,j--)
  {
    invbin[i] = invbin[j] = z;
    z = divru(mulru(z, i+1), n-i);
  }
  z = invbin[0]; /* = 1 */
  for (i=0; i<=n; i++)
  {
    c = gel(x,i); if (!signe(c)) continue;
    affir(c, r);
    z = addrr(z, mulrr(sqrr(r), invbin[i]));
  }
  z = shiftr(sqrtr(z), n);
  z = divrr(z, dbltor(pow((double)n, 0.75)));
  z = roundr_safe(sqrtr(z));
  z = mulii(z, absi_shallow(gel(x,n)));
  return gerepileuptoint(av, shifti(z, 1));
}
#endif

/* A | S ==> |a_i| <= binom(d-1, i-1) || S ||_2 + binom(d-1, i) lc(S) */
static GEN
Mignotte_bound(GEN S)
{
  long i, d = degpol(S);
  GEN C, N2, t, binlS, lS = leading_coeff(S), bin = vecbinomial(d-1);

  N2 = sqrtr(RgX_fpnorml2(S,DEFAULTPREC));
  binlS = is_pm1(lS)? bin: ZC_Z_mul(bin, lS);

  /* i = 0 */
  C = gel(binlS,1);
  /* i = d */
  t = N2; if (gcmp(C, t) < 0) C = t;
  for (i = 1; i < d; i++)
  {
    t = addri(mulir(gel(bin,i), N2), gel(binlS,i+1));
    if (mpcmp(C, t) < 0) C = t;
  }
  return C;
}
/* A | S ==> |a_i|^2 <= 3^{3/2 + d} / (4 \pi d) [P]_2^2,
 * where [P]_2 is Bombieri's 2-norm */
static GEN
Beauzamy_bound(GEN S)
{
  const long prec = DEFAULTPREC;
  long i, d = degpol(S);
  GEN bin, lS, s, C;
  bin = vecbinomial(d);

  s = real_0(prec);
  for (i=0; i<=d; i++)
  {
    GEN c = gel(S,i+2);
    if (gequal0(c)) continue;
    /* s += P_i^2 / binomial(d,i) */
    s = addrr(s, divri(itor(sqri(c), prec), gel(bin,i+1)));
  }
  /* s = [S]_2^2 */
  C = powruhalf(utor(3,prec), 3 + 2*d); /* 3^{3/2 + d} */
  C = divrr(mulrr(C, s), mulur(4*d, mppi(prec)));
  lS = absi_shallow(leading_coeff(S));
  return mulir(lS, sqrtr(C));
}

static GEN
factor_bound(GEN S)
{
  pari_sp av = avma;
  GEN a = Mignotte_bound(S);
  GEN b = Beauzamy_bound(S);
  if (DEBUGLEVEL>2)
  {
    err_printf("Mignotte bound: %Ps\n",a);
    err_printf("Beauzamy bound: %Ps\n",b);
  }
  return gerepileupto(av, ceil_safe(gmin_shallow(a, b)));
}

/* Naive recombination of modular factors: combine up to maxK modular
 * factors, degree <= klim
 *
 * target = polynomial we want to factor
 * famod = array of modular factors.  Product should be congruent to
 * target/lc(target) modulo p^a
 * For true factors: S1,S2 <= p^b, with b <= a and p^(b-a) < 2^31 */
static GEN
cmbf(GEN pol, GEN famod, GEN bound, GEN p, long a, long b,
     long klim, long *pmaxK, int *done)
{
  long K = 1, cnt = 1, i,j,k, curdeg, lfamod = lg(famod)-1;
  ulong spa_b, spa_bs2, Sbound;
  GEN lc, lcpol, pa = powiu(p,a), pas2 = shifti(pa,-1);
  GEN trace1   = cgetg(lfamod+1, t_VECSMALL);
  GEN trace2   = cgetg(lfamod+1, t_VECSMALL);
  GEN ind      = cgetg(lfamod+1, t_VECSMALL);
  GEN deg      = cgetg(lfamod+1, t_VECSMALL);
  GEN degsofar = cgetg(lfamod+1, t_VECSMALL);
  GEN listmod  = cgetg(lfamod+1, t_VEC);
  GEN fa       = cgetg(lfamod+1, t_VEC);

  *pmaxK = cmbf_maxK(lfamod);
  lc = absi_shallow(leading_coeff(pol));
  if (equali1(lc)) lc = NULL;
  lcpol = lc? ZX_Z_mul(pol, lc): pol;

  {
    GEN pa_b,pa_bs2,pb, lc2 = lc? sqri(lc): NULL;

    pa_b = powiu(p, a-b); /* < 2^31 */
    pa_bs2 = shifti(pa_b,-1);
    pb= powiu(p, b);
    for (i=1; i <= lfamod; i++)
    {
      GEN T1,T2, P = gel(famod,i);
      long d = degpol(P);

      deg[i] = d; P += 2;
      T1 = gel(P,d-1);/* = - S_1 */
      T2 = sqri(T1);
      if (d > 1) T2 = subii(T2, shifti(gel(P,d-2),1));
      T2 = modii(T2, pa); /* = S_2 Newton sum */
      if (lc)
      {
        T1 = Fp_mul(lc, T1, pa);
        T2 = Fp_mul(lc2,T2, pa);
      }
      uel(trace1,i) = itou(diviiround(T1, pb));
      uel(trace2,i) = itou(diviiround(T2, pb));
    }
    spa_b   = uel(pa_b,2); /* < 2^31 */
    spa_bs2 = uel(pa_bs2,2); /* < 2^31 */
  }
  degsofar[0] = 0; /* sentinel */

  /* ind runs through strictly increasing sequences of length K,
   * 1 <= ind[i] <= lfamod */
nextK:
  if (K > *pmaxK || 2*K > lfamod) goto END;
  if (DEBUGLEVEL > 3)
    err_printf("\n### K = %d, %Ps combinations\n", K,binomial(utoipos(lfamod), K));
  setlg(ind, K+1); ind[1] = 1;
  Sbound = (ulong) ((K+1)>>1);
  i = 1; curdeg = deg[ind[1]];
  for(;;)
  { /* try all combinations of K factors */
    for (j = i; j < K; j++)
    {
      degsofar[j] = curdeg;
      ind[j+1] = ind[j]+1; curdeg += deg[ind[j+1]];
    }
    if (curdeg <= klim) /* trial divide */
    {
      GEN y, q, list;
      pari_sp av;
      ulong t;

      /* d - 1 test */
      for (t=uel(trace1,ind[1]),i=2; i<=K; i++)
        t = Fl_add(t, uel(trace1,ind[i]), spa_b);
      if (t > spa_bs2) t = spa_b - t;
      if (t > Sbound)
      {
        if (DEBUGLEVEL>6) err_printf(".");
        goto NEXT;
      }
      /* d - 2 test */
      for (t=uel(trace2,ind[1]),i=2; i<=K; i++)
        t = Fl_add(t, uel(trace2,ind[i]), spa_b);
      if (t > spa_bs2) t = spa_b - t;
      if (t > Sbound)
      {
        if (DEBUGLEVEL>6) err_printf("|");
        goto NEXT;
      }

      av = avma;
      /* check trailing coeff */
      y = lc;
      for (i=1; i<=K; i++)
      {
        GEN q = constant_coeff(gel(famod,ind[i]));
        if (y) q = mulii(y, q);
        y = centermodii(q, pa, pas2);
      }
      if (!signe(y) || !dvdii(constant_coeff(lcpol), y))
      {
        if (DEBUGLEVEL>3) err_printf("T");
        set_avma(av); goto NEXT;
      }
      y = lc; /* full computation */
      for (i=1; i<=K; i++)
      {
        GEN q = gel(famod,ind[i]);
        if (y) q = gmul(y, q);
        y = centermod_i(q, pa, pas2);
      }

      /* y is the candidate factor */
      if (! (q = ZX_divides_i(lcpol,y,bound)) )
      {
        if (DEBUGLEVEL>3) err_printf("*");
        set_avma(av); goto NEXT;
      }
      /* found a factor */
      list = cgetg(K+1, t_VEC);
      gel(listmod,cnt) = list;
      for (i=1; i<=K; i++) list[i] = famod[ind[i]];

      y = Q_primpart(y);
      gel(fa,cnt++) = y;
      /* fix up pol */
      pol = q;
      if (lc) pol = Q_div_to_int(pol, leading_coeff(y));
      for (i=j=k=1; i <= lfamod; i++)
      { /* remove used factors */
        if (j <= K && i == ind[j]) j++;
        else
        {
          gel(famod,k) = gel(famod,i);
          uel(trace1,k) = uel(trace1,i);
          uel(trace2,k) = uel(trace2,i);
          deg[k] = deg[i]; k++;
        }
      }
      lfamod -= K;
      *pmaxK = cmbf_maxK(lfamod);
      if (lfamod < 2*K) goto END;
      i = 1; curdeg = deg[ind[1]];
      bound = factor_bound(pol);
      if (lc) lc = absi_shallow(leading_coeff(pol));
      lcpol = lc? ZX_Z_mul(pol, lc): pol;
      if (DEBUGLEVEL>3)
        err_printf("\nfound factor %Ps\nremaining modular factor(s): %ld\n",
                   y, lfamod);
      continue;
    }

NEXT:
    for (i = K+1;;)
    {
      if (--i == 0) { K++; goto nextK; }
      if (++ind[i] <= lfamod - K + i)
      {
        curdeg = degsofar[i-1] + deg[ind[i]];
        if (curdeg <= klim) break;
      }
    }
  }
END:
  *done = 1;
  if (degpol(pol) > 0)
  { /* leftover factor */
    if (signe(leading_coeff(pol)) < 0) pol = ZX_neg(pol);
    if (lfamod >= 2*K) *done = 0;

    setlg(famod, lfamod+1);
    gel(listmod,cnt) = leafcopy(famod);
    gel(fa,cnt++) = pol;
  }
  if (DEBUGLEVEL>6) err_printf("\n");
  setlg(listmod, cnt);
  setlg(fa, cnt); return mkvec2(fa, listmod);
}

/* recombination of modular factors: van Hoeij's algorithm */

/* Q in Z[X], return Q(2^n) */
static GEN
shifteval(GEN Q, long n)
{
  pari_sp av = avma;
  long i, l = lg(Q);
  GEN s;

  if (!signe(Q)) return gen_0;
  s = gel(Q,l-1);
  for (i = l-2; i > 1; i--)
  {
    s = addii(gel(Q,i), shifti(s, n));
    if (gc_needed(av,1)) s = gerepileuptoint(av, s);
  }
  return s;
}

/* return integer y such that all |a| <= y if P(a) = 0 */
static GEN
root_bound(GEN P0)
{
  GEN Q = leafcopy(P0), lP = absi_shallow(leading_coeff(Q)), x,y,z;
  long k, d = degpol(Q);

  /* P0 = lP x^d + Q, deg Q < d */
  Q = normalizepol_lg(Q, d+2);
  for (k=lg(Q)-1; k>1; k--) gel(Q,k) = absi_shallow(gel(Q,k));
  k = (long)(fujiwara_bound(P0));
  for (  ; k >= 0; k--)
  {
    pari_sp av = avma;
    /* y = 2^k; Q(y) >= lP y^d ? */
    if (cmpii(shifteval(Q,k), shifti(lP, d*k)) >= 0) break;
    set_avma(av);
  }
  if (k < 0) k = 0;
  y = int2n(k+1);
  if (d > 2000) return y; /* likely to be expensive, don't bother */
  x = int2n(k);
  for(k=0; ; k++)
  {
    z = shifti(addii(x,y), -1);
    if (equalii(x,z) || k > 5) break;
    if (cmpii(ZX_Z_eval(Q,z), mulii(lP, powiu(z, d))) < 0)
      y = z;
    else
      x = z;
  }
  return y;
}

GEN
chk_factors_get(GEN lt, GEN famod, GEN c, GEN T, GEN N)
{
  long i = 1, j, l = lg(famod);
  GEN V = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
    if (signe(gel(c,j))) gel(V,i++) = gel(famod,j);
  if (lt && i > 1) gel(V,1) = RgX_Rg_mul(gel(V,1), lt);
  setlg(V, i);
  return T? FpXQXV_prod(V, T, N): FpXV_prod(V,N);
}

static GEN
chk_factors(GEN P, GEN M_L, GEN bound, GEN famod, GEN pa)
{
  long i, r;
  GEN pol = P, list, piv, y, ltpol, lt, paov2;

  piv = ZM_hnf_knapsack(M_L);
  if (!piv) return NULL;
  if (DEBUGLEVEL>7) err_printf("ZM_hnf_knapsack output:\n%Ps\n",piv);

  r  = lg(piv)-1;
  list = cgetg(r+1, t_VEC);
  lt = absi_shallow(leading_coeff(pol));
  if (equali1(lt)) lt = NULL;
  ltpol = lt? ZX_Z_mul(pol, lt): pol;
  paov2 = shifti(pa,-1);
  for (i = 1;;)
  {
    if (DEBUGLEVEL) err_printf("LLL_cmbf: checking factor %ld\n",i);
    y = chk_factors_get(lt, famod, gel(piv,i), NULL, pa);
    y = FpX_center_i(y, pa, paov2);
    if (! (pol = ZX_divides_i(ltpol,y,bound)) ) return NULL;
    if (lt) y = Q_primpart(y);
    gel(list,i) = y;
    if (++i >= r) break;

    if (lt)
    {
      pol = ZX_Z_divexact(pol, leading_coeff(y));
      lt = absi_shallow(leading_coeff(pol));
      ltpol = ZX_Z_mul(pol, lt);
    }
    else
      ltpol = pol;
  }
  y = Q_primpart(pol);
  gel(list,i) = y; return list;
}

GEN
LLL_check_progress(GEN Bnorm, long n0, GEN m, int final, long *ti_LLL)
{
  GEN norm, u;
  long i, R;
  pari_timer T;

  if (DEBUGLEVEL>2) timer_start(&T);
  u = ZM_lll_norms(m, final? 0.999: 0.75, LLL_INPLACE, &norm);
  if (DEBUGLEVEL>2) *ti_LLL += timer_delay(&T);
  for (R=lg(m)-1; R > 0; R--)
    if (cmprr(gel(norm,R), Bnorm) < 0) break;
  for (i=1; i<=R; i++) setlg(u[i], n0+1);
  if (R <= 1)
  {
    if (!R) pari_err_BUG("LLL_cmbf [no factor]");
    return NULL; /* irreducible */
  }
  setlg(u, R+1); return u;
}

static ulong
next2pow(ulong a)
{
  ulong b = 1;
  while (b < a) b <<= 1;
  return b;
}

/* Recombination phase of Berlekamp-Zassenhaus algorithm using a variant of
 * van Hoeij's knapsack
 *
 * P = squarefree in Z[X].
 * famod = array of (lifted) modular factors mod p^a
 * bound = Mignotte bound for the size of divisors of P (for the sup norm)
 * previously recombined all set of factors with less than rec elts */
static GEN
LLL_cmbf(GEN P, GEN famod, GEN p, GEN pa, GEN bound, long a, long rec)
{
  const long N0 = 1; /* # of traces added at each step */
  double BitPerFactor = 0.4; /* nb bits in p^(a-b) / modular factor */
  long i,j,tmax,n0,C, dP = degpol(P);
  double logp = log((double)itos(p)), LOGp2 = M_LN2/logp;
  double b0 = log((double)dP*2) / logp, logBr;
  GEN lP, Br, Bnorm, Tra, T2, TT, CM_L, m, list, ZERO;
  pari_sp av, av2;
  long ti_LLL = 0, ti_CF  = 0;

  lP = absi_shallow(leading_coeff(P));
  if (equali1(lP)) lP = NULL;
  Br = root_bound(P);
  if (lP) Br = mulii(lP, Br);
  logBr = gtodouble(glog(Br, DEFAULTPREC)) / logp;

  n0 = lg(famod) - 1;
  C = (long)ceil( sqrt(N0 * n0 / 4.) ); /* > 1 */
  Bnorm = dbltor(n0 * (C*C + N0*n0/4.) * 1.00001);
  ZERO = zeromat(n0, N0);

  av = avma;
  TT = cgetg(n0+1, t_VEC);
  Tra  = cgetg(n0+1, t_MAT);
  for (i=1; i<=n0; i++)
  {
    TT[i]  = 0;
    gel(Tra,i) = cgetg(N0+1, t_COL);
  }
  CM_L = scalarmat_s(C, n0);
  /* tmax = current number of traces used (and computed so far) */
  for (tmax = 0;; tmax += N0)
  {
    long b, bmin, bgood, delta, tnew = tmax + N0, r = lg(CM_L)-1;
    GEN M_L, q, CM_Lp, oldCM_L;
    int first = 1;
    pari_timer ti2, TI;

    bmin = (long)ceil(b0 + tnew*logBr);
    if (DEBUGLEVEL>2)
      err_printf("\nLLL_cmbf: %ld potential factors (tmax = %ld, bmin = %ld)\n",
                 r, tmax, bmin);

    /* compute Newton sums (possibly relifting first) */
    if (a <= bmin)
    {
      a = (long)ceil(bmin + 3*N0*logBr) + 1; /* enough for 3 more rounds */
      a = (long)next2pow((ulong)a);

      pa = powiu(p,a);
      famod = ZpX_liftfact(P, famod, pa, p, a);
      for (i=1; i<=n0; i++) TT[i] = 0;
    }
    for (i=1; i<=n0; i++)
    {
      GEN p1 = gel(Tra,i);
      GEN p2 = polsym_gen(gel(famod,i), gel(TT,i), tnew, NULL, pa);
      gel(TT,i) = p2;
      p2 += 1+tmax; /* ignore traces number 0...tmax */
      for (j=1; j<=N0; j++) gel(p1,j) = gel(p2,j);
      if (lP)
      { /* make Newton sums integral */
        GEN lPpow = powiu(lP, tmax);
        for (j=1; j<=N0; j++)
        {
          lPpow = mulii(lPpow,lP);
          gel(p1,j) = mulii(gel(p1,j), lPpow);
        }
      }
    }

    /* compute truncation parameter */
    if (DEBUGLEVEL>2) { timer_start(&ti2); timer_start(&TI); }
    oldCM_L = CM_L;
    av2 = avma;
    delta = b = 0; /* -Wall */
AGAIN:
    M_L = Q_div_to_int(CM_L, utoipos(C));
    T2 = centermod( ZM_mul(Tra, M_L), pa );
    if (first)
    { /* initialize lattice, using few p-adic digits for traces */
      double t = gexpo(T2) - maxdd(32.0, BitPerFactor*r);
      bgood = (long) (t * LOGp2);
      b = maxss(bmin, bgood);
      delta = a - b;
    }
    else
    { /* add more p-adic digits and continue reduction */
      long b0 = (long)(gexpo(T2) * LOGp2);
      if (b0 < b) b = b0;
      b = maxss(b-delta, bmin);
      if (b - delta/2 < bmin) b = bmin; /* near there. Go all the way */
    }

    q = powiu(p, b);
    m = vconcat( CM_L, gdivround(T2, q) );
    if (first)
    {
      GEN P1 = scalarmat(powiu(p, a-b), N0);
      first = 0;
      m = shallowconcat( m, vconcat(ZERO, P1) );
      /*     [ C M_L        0     ]
       * m = [                    ]   square matrix
       *     [  T2'  p^(a-b) I_N0 ]   T2' = Tra * M_L  truncated
       */
    }

    CM_L = LLL_check_progress(Bnorm, n0, m, b == bmin, /*dbg:*/ &ti_LLL);
    if (DEBUGLEVEL>2)
      err_printf("LLL_cmbf: (a,b) =%4ld,%4ld; r =%3ld -->%3ld, time = %ld\n",
                 a,b, lg(m)-1, CM_L? lg(CM_L)-1: 1, timer_delay(&TI));
    if (!CM_L) { list = mkvec(P); break; }
    if (b > bmin)
    {
      CM_L = gerepilecopy(av2, CM_L);
      goto AGAIN;
    }
    if (DEBUGLEVEL>2) timer_printf(&ti2, "for this block of traces");

    i = lg(CM_L) - 1;
    if (i == r && ZM_equal(CM_L, oldCM_L))
    {
      CM_L = oldCM_L;
      set_avma(av2); continue;
    }

    CM_Lp = FpM_image(CM_L, utoipos(27449)); /* inexpensive test */
    if (lg(CM_Lp) != lg(CM_L))
    {
      if (DEBUGLEVEL>2) err_printf("LLL_cmbf: rank decrease\n");
      CM_L = ZM_hnf(CM_L);
    }

    if (i <= r && i*rec < n0)
    {
      pari_timer ti;
      if (DEBUGLEVEL>2) timer_start(&ti);
      list = chk_factors(P, Q_div_to_int(CM_L,utoipos(C)), bound, famod, pa);
      if (DEBUGLEVEL>2) ti_CF += timer_delay(&ti);
      if (list) break;
      if (DEBUGLEVEL>2) err_printf("LLL_cmbf: chk_factors failed");
    }
    CM_L = gerepilecopy(av2, CM_L);
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"LLL_cmbf");
      gerepileall(av, 5, &CM_L, &TT, &Tra, &famod, &pa);
    }
  }
  if (DEBUGLEVEL>2)
    err_printf("* Time LLL: %ld\n* Time Check Factor: %ld\n",ti_LLL,ti_CF);
  return list;
}

/* Find a,b minimal such that A < q^a, B < q^b, 1 << q^(a-b) < 2^31 */
static int
cmbf_precs(GEN q, GEN A, GEN B, long *pta, long *ptb, GEN *qa, GEN *qb)
{
  long a,b,amin,d = (long)(31 * M_LN2/gtodouble(glog(q,DEFAULTPREC)) - 1e-5);
  int fl = 0;

  b = logintall(B, q, qb) + 1;
  *qb = mulii(*qb, q);
  amin = b + d;
  if (gcmp(powiu(q, amin), A) <= 0)
  {
    a = logintall(A, q, qa) + 1;
    *qa = mulii(*qa, q);
    b = a - d; *qb = powiu(q, b);
  }
  else
  { /* not enough room */
    a = amin;  *qa = powiu(q, a);
    fl = 1;
  }
  if (DEBUGLEVEL > 3) {
    err_printf("S_2   bound: %Ps^%ld\n", q,b);
    err_printf("coeff bound: %Ps^%ld\n", q,a);
  }
  *pta = a;
  *ptb = b; return fl;
}

/* use van Hoeij's knapsack algorithm */
static GEN
combine_factors(GEN target, GEN famod, GEN p, long klim)
{
  GEN la, B, A, res, L, pa, pb, listmod;
  long a,b, l, maxK, n = degpol(target);
  int done;
  pari_timer T;

  A = factor_bound(target);

  la = absi_shallow(leading_coeff(target));
  B = mului(n, sqri(mulii(la, root_bound(target)))); /* = bound for S_2 */

  (void)cmbf_precs(p, A, B, &a, &b, &pa, &pb);

  if (DEBUGLEVEL>2) timer_start(&T);
  famod = ZpX_liftfact(target, famod, pa, p, a);
  if (DEBUGLEVEL>2) timer_printf(&T, "Hensel lift (mod %Ps^%ld)", p,a);
  L = cmbf(target, famod, A, p, a, b, klim, &maxK, &done);
  if (DEBUGLEVEL>2) timer_printf(&T, "Naive recombination");

  res     = gel(L,1);
  listmod = gel(L,2); l = lg(listmod)-1;
  famod = gel(listmod,l);
  if (maxK > 0 && lg(famod)-1 > 2*maxK)
  {
    if (l!=1) A = factor_bound(gel(res,l));
    if (DEBUGLEVEL > 4) err_printf("last factor still to be checked\n");
    L = LLL_cmbf(gel(res,l), famod, p, pa, A, a, maxK);
    if (DEBUGLEVEL>2) timer_printf(&T,"Knapsack");
    /* remove last elt, possibly unfactored. Add all new ones. */
    setlg(res, l); res = shallowconcat(res, L);
  }
  return res;
}

/* Assume 'a' a squarefree ZX; return 0 if no root (fl=1) / irreducible (fl=0).
 * Otherwise return prime p such that a mod p has fewest roots / factors */
static ulong
pick_prime(GEN a, long fl, pari_timer *T)
{
  pari_sp av = avma, av1;
  const long MAXNP = 7, da = degpol(a);
  long nmax = da+1, np;
  ulong chosenp = 0;
  GEN lead = gel(a,da+2);
  forprime_t S;
  if (equali1(lead)) lead = NULL;
  u_forprime_init(&S, 2, ULONG_MAX);
  av1 = avma;
  for (np = 0; np < MAXNP; set_avma(av1))
  {
    ulong p = u_forprime_next(&S);
    long nfacp;
    GEN z;

    if (!p) pari_err_OVERFLOW("DDF [out of small primes]");
    if (lead && !umodiu(lead,p)) continue;
    z = ZX_to_Flx(a, p);
    if (!Flx_is_squarefree(z, p)) continue;

    if (fl==1)
    {
      nfacp = Flx_nbroots(z, p);
      if (!nfacp) { chosenp = 0; break; } /* no root */
    }
    else if(fl==0)
    {
      nfacp = Flx_nbfact(z, p);
      if (nfacp == 1) { chosenp = 0; break; } /* irreducible */
    } else
    {
      GEN f = gel(Flx_degfact(z, p),1);
      nfacp = lg(f)-1;
      if (f[1] > fl) { chosenp = 0; break; } /* no small factors */
    }
    if (DEBUGLEVEL>4)
      err_printf("...tried prime %3lu (%-3ld %s). Time = %ld\n",
                  p, nfacp, fl==1? "roots": "factors", timer_delay(T));
    if (nfacp < nmax)
    {
      nmax = nfacp; chosenp = p;
      if (da > 100 && nmax < 5) break; /* large degree, few factors. Enough */
    }
    np++;
  }
  return gc_ulong(av, chosenp);
}

/* Assume A a squarefree ZX; return the vector of its rational roots */
static GEN
DDF_roots(GEN A)
{
  GEN p, lc, lcpol, z, pe, pes2, bound;
  long i, m, e, lz;
  ulong pp;
  pari_sp av;
  pari_timer T;

  if (DEBUGLEVEL>2) timer_start(&T);
  pp = pick_prime(A, 1, &T);
  if (!pp) return cgetg(1,t_COL); /* no root */
  p = utoipos(pp);
  lc = leading_coeff(A);
  if (is_pm1(lc))
  { lc = NULL; lcpol = A; }
  else
  { lc = absi_shallow(lc); lcpol = ZX_Z_mul(A, lc); }
  bound = root_bound(A); if (lc) bound = mulii(lc, bound);
  e = logintall(addiu(shifti(bound, 1), 1), p, &pe) + 1;
  pe = mulii(pe, p);
  pes2 = shifti(pe, -1);
  if (DEBUGLEVEL>2) timer_printf(&T, "Root bound");
  av = avma;
  z = ZpX_roots(A, p, e); lz = lg(z);
  z = deg1_from_roots(z, varn(A));
  if (DEBUGLEVEL>2) timer_printf(&T, "Hensel lift (mod %lu^%ld)", pp,e);
  for (m=1, i=1; i < lz; i++)
  {
    GEN q, r, y = gel(z,i);
    if (lc) y = ZX_Z_mul(y, lc);
    y = centermod_i(y, pe, pes2);
    if (! (q = ZX_divides(lcpol, y)) ) continue;

    lcpol = q;
    r = negi( constant_coeff(y) );
    if (lc) {
      r = gdiv(r,lc);
      lcpol = Q_primpart(lcpol);
      lc = absi_shallow( leading_coeff(lcpol) );
      if (is_pm1(lc)) lc = NULL; else lcpol = ZX_Z_mul(lcpol, lc);
    }
    gel(z,m++) = r;
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"DDF_roots, m = %ld", m);
      gerepileall(av, lc? 3:2, &z, &lcpol, &lc);
    }
  }
  if (DEBUGLEVEL>2) timer_printf(&T, "Recombination");
  setlg(z, m); return z;
}

/* Assume a squarefree ZX, deg(a) > 0, return rational factors.
 * In fact, a(0) != 0 but we don't use this
 * if dmax>0, Only look for factor of degree at most dmax */
GEN
ZX_DDF_max(GEN a, long dmax)
{
  GEN ap, prime, famod, z;
  long ti = 0;
  ulong p = 0;
  pari_sp av = avma;
  pari_timer T, T2;

  if (DEBUGLEVEL>2) { timer_start(&T); timer_start(&T2); }
  p = pick_prime(a, dmax, &T2);
  if (!p) return mkvec(a);
  prime = utoipos(p);
  ap = Flx_normalize(ZX_to_Flx(a, p), p);
  famod = gel(Flx_factor(ap, p), 1);
  if (DEBUGLEVEL>2)
  {
    if (DEBUGLEVEL>4) timer_printf(&T2, "splitting mod p = %lu", p);
    ti = timer_delay(&T);
    err_printf("Time setup: %ld\n", ti);
  }
  z = combine_factors(a, FlxV_to_ZXV(famod), prime, degpol(a)-1);
  if (DEBUGLEVEL>2)
    err_printf("Total Time: %ld\n===========\n", ti + timer_delay(&T));
  return gerepilecopy(av, z);
}

/* Distinct Degree Factorization (deflating first)
 * Assume x squarefree, degree(x) > 0, x(0) != 0 */
GEN
ZX_DDF(GEN x)
{
  GEN L;
  long m;
  x = ZX_deflate_max(x, &m);
  L = ZX_DDF_max(x,0);
  if (m > 1)
  {
    GEN e, v, fa = factoru(m);
    long i,j,k, l;

    e = gel(fa,2); k = 0;
    fa= gel(fa,1); l = lg(fa);
    for (i=1; i<l; i++) k += e[i];
    v = cgetg(k+1, t_VECSMALL); k = 1;
    for (i=1; i<l; i++)
      for (j=1; j<=e[i]; j++) v[k++] = fa[i];
    for (k--; k; k--)
    {
      GEN L2 = cgetg(1,t_VEC);
      for (i=1; i < lg(L); i++)
              L2 = shallowconcat(L2, ZX_DDF_max(RgX_inflate(gel(L,i), v[k]),0));
      L = L2;
    }
  }
  return L;
}

/* SquareFree Factorization in Z[X] (char 0 is enough, if ZX_gcd -> RgX_gcd)
 * f = prod Q[i]^E[i], E[1] < E[2] < ..., and Q[i] squarefree and coprime.
 * Return Q, set *pE = E. For efficiency, caller should have used ZX_valrem
 * so that f(0) != 0 */
GEN
ZX_squff(GEN f, GEN *pE)
{
  GEN T, V, P, E;
  long i, k, n = 1 + degpol(f);

  if (signe(leading_coeff(f)) < 0) f = ZX_neg(f);
  E = cgetg(n, t_VECSMALL);
  P = cgetg(n, t_COL);
  T = ZX_gcd_all(f, ZX_deriv(f), &V);
  for (k = i = 1;; k++)
  {
    GEN W = ZX_gcd_all(T,V, &T); /* V and W are squarefree */
    long dW = degpol(W), dV = degpol(V);
    /* f = prod_i T_i^{e_i}
     * W = prod_{i: e_i > k} T_i,
     * V = prod_{i: e_i >= k} T_i,
     * T = prod_{i: e_i > k} T_i^{e_i - k} */
    if (!dW)
    {
      if (dV) { gel(P,i) = Q_primpart(V); E[i] = k; i++; }
      break;
    }
    if (dW == dV)
    {
      GEN U;
      while ( (U = ZX_divides(T, V)) ) { k++; T = U; }
    }
    else
    {
      gel(P,i) = Q_primpart(RgX_div(V,W));
      E[i] = k; i++; V = W;
    }
  }
  setlg(P,i);
  setlg(E,i); *pE = E; return P;
}

static GEN
fact_from_DDF(GEN Q, GEN E, long n)
{
  GEN v,w, y = cgetg(3, t_MAT);
  long i,j,k, l = lg(Q);

  v = cgetg(n+1, t_COL); gel(y,1) = v;
  w = cgetg(n+1, t_COL); gel(y,2) = w;
  for (k = i = 1; i < l; i++)
  {
    GEN L = gel(Q,i), e = utoipos(E[i]);
    long J = lg(L);
    for (j = 1; j < J; j++,k++)
    {
      gel(v,k) = ZX_copy(gel(L,j));
      gel(w,k) = e;
    }
  }
  return y;
}

/* Factor T in Z[x] */
static GEN
ZX_factor_i(GEN T)
{
  GEN Q, E, y;
  long n, i, l, v;

  if (!signe(T)) return prime_fact(T);
  v = ZX_valrem(T, &T);
  Q = ZX_squff(T, &E); l = lg(Q);
  for (i = 1, n = 0; i < l; i++)
  {
    gel(Q,i) = ZX_DDF(gel(Q,i));
    n += lg(gel(Q,i)) - 1;
  }
  if (v)
  {
    Q = vec_append(Q, mkvec(pol_x(varn(T))));
    E = vecsmall_append(E, v); n++;
  }
  y = fact_from_DDF(Q, E, n);
  return sort_factor_pol(y, cmpii);
}
GEN
ZX_factor(GEN x)
{
  pari_sp av = avma;
  return gerepileupto(av, ZX_factor_i(x));
}
GEN
QX_factor(GEN x)
{
  pari_sp av = avma;
  return gerepileupto(av, ZX_factor_i(Q_primpart(x)));
}

long
ZX_is_irred(GEN x)
{
  pari_sp av = avma;
  long l = lg(x);
  GEN y;
  if (l <= 3) return 0; /* degree < 1 */
  if (l == 4) return 1; /* degree 1 */
  if (ZX_val(x)) return 0;
  if (!ZX_is_squarefree(x)) return 0;
  y = ZX_DDF(x); set_avma(av);
  return (lg(y) == 2);
}

GEN
nfrootsQ(GEN x)
{
  pari_sp av = avma;
  GEN z;
  long val;

  if (typ(x)!=t_POL) pari_err_TYPE("nfrootsQ",x);
  if (!signe(x)) pari_err_ROOTS0("nfrootsQ");
  x = Q_primpart(x);
  RgX_check_ZX(x,"nfrootsQ");
  val = ZX_valrem(x, &x);
  z = DDF_roots( ZX_radical(x) );
  if (val) z = vec_append(z, gen_0);
  return gerepileupto(av, sort(z));
}

/************************************************************************
 *                   GCD OVER Z[X] / Q[X]                               *
 ************************************************************************/
int
ZX_is_squarefree(GEN x)
{
  pari_sp av = avma;
  GEN d;
  long m;
  if (lg(x) == 2) return 0;
  m = ZX_deflate_order(x);
  if (m > 1)
  {
    if (!signe(gel(x,2))) return 0;
    x = RgX_deflate(x, m);
  }
  d = ZX_gcd(x,ZX_deriv(x));
  return gc_bool(av, lg(d) == 3);
}

static int
ZX_gcd_filter(GEN *pt_A, GEN *pt_P)
{
  GEN A = *pt_A, P = *pt_P;
  long i, j, l = lg(A), n = 1, d = degpol(gel(A,1));
  GEN B, Q;
  for (i=2; i<l; i++)
  {
    long di = degpol(gel(A,i));
    if (di==d) n++;
    else if (d > di)
    { n=1; d = di; }
  }
  if (n == l-1)
    return 0;
  B = cgetg(n+1, t_VEC);
  Q = cgetg(n+1, typ(P));
  for (i=1, j=1; i<l; i++)
  {
    if (degpol(gel(A,i))==d)
    {
      gel(B,j) = gel(A,i);
      Q[j] = P[i];
      j++;
    }
  }
  *pt_A = B; *pt_P = Q; return 1;
}

static GEN
ZX_gcd_Flx(GEN a, GEN b, ulong g, ulong p)
{
  GEN H = Flx_gcd(a, b, p);
  if (!g)
    return Flx_normalize(H, p);
  else
  {
    ulong t = Fl_mul(g, Fl_inv(Flx_lead(H), p), p);
    return Flx_Fl_mul(H, t, p);
  }
}

static GEN
ZX_gcd_slice(GEN A, GEN B, GEN g, GEN P, GEN *mod)
{
  pari_sp av = avma;
  long i, n = lg(P)-1;
  GEN H, T;
  if (n == 1)
  {
    ulong p = uel(P,1), gp = g ? umodiu(g, p): 0;
    GEN a = ZX_to_Flx(A, p), b = ZX_to_Flx(B, p);
    GEN Hp = ZX_gcd_Flx(a, b, gp, p);
    H = gerepileupto(av, Flx_to_ZX(Hp));
    *mod = utoi(p);
    return H;
  }
  T = ZV_producttree(P);
  A = ZX_nv_mod_tree(A, P, T);
  B = ZX_nv_mod_tree(B, P, T);
  g = g ?  Z_ZV_mod_tree(g, P, T): NULL;
  H = cgetg(n+1, t_VEC);
  for(i=1; i <= n; i++)
  {
    ulong p = P[i];
    GEN a = gel(A,i), b = gel(B,i);
    gel(H,i) = ZX_gcd_Flx(a, b, g? g[i]: 0, p);
  }
  if (ZX_gcd_filter(&H, &P))
    T = ZV_producttree(P);
  H = nxV_chinese_center_tree(H, P, T, ZV_chinesetree(P, T));
  *mod = gmael(T, lg(T)-1, 1); return gc_all(av, 2, &H, mod);
}

GEN
ZX_gcd_worker(GEN P, GEN A, GEN B, GEN g)
{
  GEN V = cgetg(3, t_VEC);
  gel(V,1) = ZX_gcd_slice(A, B, equali1(g)? NULL: g , P, &gel(V,2));
  return V;
}

static GEN
ZX_gcd_chinese(GEN A, GEN P, GEN *mod)
{
  ZX_gcd_filter(&A, &P);
  return nxV_chinese_center(A, P, mod);
}

GEN
ZX_gcd_all(GEN A, GEN B, GEN *Anew)
{
  pari_sp av = avma;
  long k, valH, valA, valB, vA = varn(A), dA = degpol(A), dB = degpol(B);
  GEN worker, c, cA, cB, g, Ag, Bg, H = NULL, mod = gen_1, R;
  GEN Ap, Bp, Hp;
  forprime_t S;
  ulong pp;
  if (dA < 0) { if (Anew) *Anew = pol_0(vA); return ZX_copy(B); }
  if (dB < 0) { if (Anew) *Anew = pol_1(vA); return ZX_copy(A); }
  A = Q_primitive_part(A, &cA);
  B = Q_primitive_part(B, &cB);
  valA = ZX_valrem(A, &A); dA -= valA;
  valB = ZX_valrem(B, &B); dB -= valB;
  valH = minss(valA, valB);
  valA -= valH; /* valuation(Anew) */
  c = (cA && cB)? gcdii(cA, cB): NULL; /* content(gcd) */
  if (!dA || !dB)
  {
    if (Anew) *Anew = RgX_shift_shallow(A, valA);
    return monomial(c? c: gen_1, valH, vA);
  }
  g = gcdii(leading_coeff(A), leading_coeff(B)); /* multiple of lead(gcd) */
  if (is_pm1(g)) {
    g = NULL;
    Ag = A;
    Bg = B;
  } else {
    Ag = ZX_Z_mul(A,g);
    Bg = ZX_Z_mul(B,g);
  }
  init_modular_big(&S);
  do {
    pp = u_forprime_next(&S);
    Ap = ZX_to_Flx(Ag, pp);
    Bp = ZX_to_Flx(Bg, pp);
  } while (degpol(Ap) != dA || degpol(Bp) != dB);
  if (degpol(Flx_gcd(Ap, Bp, pp)) == 0)
  {
    if (Anew) *Anew = RgX_shift_shallow(A, valA);
    return monomial(c? c: gen_1, valH, vA);
  }
  worker = snm_closure(is_entry("_ZX_gcd_worker"), mkvec3(A, B, g? g: gen_1));
  av = avma;
  for (k = 1; ;k *= 2)
  {
    gen_inccrt_i("ZX_gcd", worker, g, (k+1)>>1, 0, &S, &H, &mod, ZX_gcd_chinese, NULL);
    gerepileall(av, 2, &H, &mod);
    Hp = ZX_to_Flx(H, pp);
    if (lgpol(Flx_rem(Ap, Hp, pp)) || lgpol(Flx_rem(Bp, Hp, pp))) continue;
    if (!ZX_divides(Bg, H)) continue;
    R = ZX_divides(Ag, H);
    if (R) break;
  }
  /* lead(H) = g */
  if (g) H = Q_primpart(H);
  if (c) H = ZX_Z_mul(H,c);
  if (DEBUGLEVEL>5) err_printf("done\n");
  if (Anew) *Anew = RgX_shift_shallow(R, valA);
  return valH? RgX_shift_shallow(H, valH): H;
}

#if 0
/* ceil( || p ||_oo / lc(p) ) */
static GEN
maxnorm(GEN p)
{
  long i, n = degpol(p), av = avma;
  GEN x, m = gen_0;

  p += 2;
  for (i=0; i<n; i++)
  {
    x = gel(p,i);
    if (abscmpii(x,m) > 0) m = x;
  }
  m = divii(m, gel(p,n));
  return gerepileuptoint(av, addiu(absi_shallow(m),1));
}
#endif

GEN
ZX_gcd(GEN A, GEN B)
{
  pari_sp av = avma;
  return gerepilecopy(av, ZX_gcd_all(A,B,NULL));
}

GEN
ZX_radical(GEN A) { GEN B; (void)ZX_gcd_all(A,ZX_deriv(A),&B); return B; }

static GEN
_gcd(GEN a, GEN b)
{
  if (!a) a = gen_1;
  if (!b) b = gen_1;
  return Q_gcd(a,b);
}
/* A0 and B0 in Q[X] */
GEN
QX_gcd(GEN A0, GEN B0)
{
  GEN a, b, D;
  pari_sp av = avma, av2;

  D = ZX_gcd(Q_primitive_part(A0, &a), Q_primitive_part(B0, &b));
  av2 = avma; a = _gcd(a,b);
  if (isint1(a)) set_avma(av2); else D = ZX_Q_mul(D, a);
  return gerepileupto(av, D);
}

/*****************************************************************************
 * Variants of the Bradford-Davenport algorithm: look for cyclotomic         *
 * factors, and decide whether a ZX is cyclotomic or a product of cyclotomic *
 *****************************************************************************/
/* f of degree 1, return a cyclotomic factor (Phi_1 or Phi_2) or NULL */
static GEN
BD_deg1(GEN f)
{
  GEN a = gel(f,3), b = gel(f,2); /* f = ax + b */
  if (!absequalii(a,b)) return NULL;
  return polcyclo((signe(a) == signe(b))? 2: 1, varn(f));
}

/* f a squarefree ZX; not divisible by any Phi_n, n even */
static GEN
BD_odd(GEN f)
{
  while(degpol(f) > 1)
  {
    GEN f1 = ZX_graeffe(f); /* contain all cyclotomic divisors of f */
    if (ZX_equal(f1, f)) return f; /* product of cyclotomics */
    f = ZX_gcd(f, f1);
  }
  if (degpol(f) == 1) return BD_deg1(f);
  return NULL; /* no cyclotomic divisor */
}

static GEN
myconcat(GEN v, GEN x)
{
  if (typ(x) != t_VEC) x = mkvec(x);
  if (!v) return x;
  return shallowconcat(v, x);
}

/* Bradford-Davenport algorithm.
 * f a squarefree ZX of degree > 0, return NULL or a vector of coprime
 * cyclotomic factors of f [ possibly reducible ] */
static GEN
BD(GEN f)
{
  GEN G = NULL, Gs = NULL, Gp = NULL, Gi = NULL;
  GEN fs2, fp, f2, f1, fe, fo, fe1, fo1;
  RgX_even_odd(f, &fe, &fo);
  fe1 = ZX_eval1(fe);
  fo1 = ZX_eval1(fo);
  if (absequalii(fe1, fo1)) /* f(1) = 0 or f(-1) = 0 */
  {
    long i, v = varn(f);
    if (!signe(fe1))
      G = mkvec2(polcyclo(1, v), polcyclo(2, v)); /* both 0 */
    else if (signe(fe1) == signe(fo1))
      G = mkvec(polcyclo(2, v)); /*f(-1) = 0*/
    else
      G = mkvec(polcyclo(1, v)); /*f(1) = 0*/
    for (i = lg(G)-1; i; i--) f = RgX_div(f, gel(G,i));
  }
  /* f no longer divisible by Phi_1 or Phi_2 */
  if (degpol(f) <= 1) return G;
  f1 = ZX_graeffe(f); /* has at most square factors */
  if (ZX_equal(f1, f)) return myconcat(G,f); /* f = product of Phi_n, n odd */

  fs2 = ZX_gcd_all(f1, ZX_deriv(f1), &f2); /* fs2 squarefree */
  if (degpol(fs2))
  { /* fs contains all Phi_n | f, 4 | n; and only those */
    /* In that case, Graeffe(Phi_n) = Phi_{n/2}^2, and Phi_n = Phi_{n/2}(x^2) */
    GEN fs = RgX_inflate(fs2, 2);
    (void)ZX_gcd_all(f, fs, &f); /* remove those Phi_n | f, 4 | n */
    Gs = BD(fs2);
    if (Gs)
    {
      long i;
      for (i = lg(Gs)-1; i; i--) gel(Gs,i) = RgX_inflate(gel(Gs,i), 2);
      /* prod Gs[i] is the product of all Phi_n | f, 4 | n */
      G = myconcat(G, Gs);
    }
    /* f2 = f1 / fs2 */
    f1 = RgX_div(f2, fs2); /* f1 / fs2^2 */
  }
  fp = ZX_gcd(f, f1); /* contains all Phi_n | f, n > 1 odd; and only those */
  if (degpol(fp))
  {
    Gp = BD_odd(fp);
    /* Gp is the product of all Phi_n | f, n odd */
    if (Gp) G = myconcat(G, Gp);
    f = RgX_div(f, fp);
  }
  if (degpol(f))
  { /* contains all Phi_n originally dividing f, n = 2 mod 4, n > 2;
     * and only those
     * In that case, Graeffe(Phi_n) = Phi_{n/2}, and Phi_n = Phi_{n/2}(-x) */
    Gi = BD_odd(ZX_z_unscale(f, -1));
    if (Gi)
    { /* N.B. Phi_2 does not divide f */
      Gi = ZX_z_unscale(Gi, -1);
      /* Gi is the product of all Phi_n | f, n = 2 mod 4 */
      G = myconcat(G, Gi);
    }
  }
  return G;
}

/* Let f be a nonzero QX, return the (squarefree) product of cyclotomic
 * divisors of f */
GEN
polcyclofactors(GEN f)
{
  pari_sp av = avma;
  if (typ(f) != t_POL || !signe(f)) pari_err_TYPE("polcyclofactors",f);
  (void)RgX_valrem(f, &f);
  f = Q_primpart(f);
  RgX_check_ZX(f,"polcyclofactors");
  if (degpol(f))
  {
    f = BD(ZX_radical(f));
    if (f) return gerepilecopy(av, f);
  }
  set_avma(av); return cgetg(1,t_VEC);
}

/* list of all squarefree odd x such that phi(x) = n, P^-(x) > m. Unsorted */
static GEN
invphi(ulong n, ulong m)
{
  GEN C, D;
  long l, i;
  if (n == 1) return mkvecsmall(1);
  D = divisorsu(n); l = lg(D);
  C = cgetg(1, t_VECSMALL);
  for (i = 2; i < l; i++) /* skip 1 */
  {
    ulong d = D[i], p;
    if (d < m) continue;
    p = d + 1; if (!uisprime(p)) continue;
    C = vecsmall_concat(C, zv_z_mul(invphi(D[l-i], p), p));
  }
  return C;
}

long
poliscyclo(GEN f)
{
  const ulong p = 2147483647; /* prime */
  pari_sp av;
  long i, n, e, l;
  ulong f3, fm3;
  GEN D, fp, _3;
  if (typ(f) != t_POL) pari_err_TYPE("poliscyclo", f);
  n = degpol(f);
  if (n <= 0 || !RgX_is_ZX(f)) return 0;
  if (!equali1(gel(f,n+2)) || !is_pm1(gel(f,2))) return 0;
  if (n == 1) return signe(gel(f,2)) > 0? 2: 1;
  av = avma;
  f = ZX_deflate_max(f, &e); if (e != 1) n = degpol(f);
  D = invphi(n, 1); /* squareefree odd d s.t. phi(d) = n */
  l = lg(D); _3 = gmodulss(3, p);
  fp = ZX_to_Flx(f, p);
  f3 = Flx_eval(fp, 3, p);
  fm3 = Flx_eval(fp, p-3, p);
  /* f(x^e) is cyclotomic (= Phi_{de}) iff f = Phi_d, where all prime dividing
   * e also divide d. */
  for (i = 1; i < l; i++)
  {
    long d = D[i]; /* squarefree odd */
    if (odd(e))
    {
      if (e == 1 || u_ppo(e, d) == 1)
      { /* early abort: check whether f(3) = Phi_d(3) or Phi_2d(3) = Phi_d(-3)
         * mod p before checking in Z. N.B. phi(d) and value at 3 mod p
         * determine Phi_d for all d <= 10^7 */
        ulong F3 = Rg_to_Fl(polcyclo_eval(d, _3), p);
        if (F3 == f3 && ZX_equal(f, polcyclo(d, varn(f))))
          return gc_long(av, d * e);
        if (F3 == fm3 && ZX_equal(f, polcyclo(2*d, varn(f))))
          return gc_long(av, 2* d * e);
      }
    }
    else
    {
      if (u_ppo(e, 2*d) == 1)
      { /* early abort: check whether f(3) = Phi_2d(3) mod p */
        ulong F3 = Rg_to_Fl(polcyclo_eval(2*d, _3), p);
        if (F3 == f3 && ZX_equal(f, polcyclo(2*d, varn(f))))
          return gc_long(av, 2* d * e);
      }
    }
  }
  return gc_long(av, 0);
}

long
poliscycloprod(GEN f)
{
  pari_sp av = avma;
  long i, d = degpol(f);
  if (typ(f) != t_POL) pari_err_TYPE("poliscycloprod",f);
  if (!RgX_is_ZX(f)) return 0;
  if (!ZX_is_monic(f) || !is_pm1(constant_coeff(f))) return 0;
  if (d < 2) return (d == 1);
  if ( degpol(ZX_gcd_all(f, ZX_deriv(f), &f)) )
  {
    d = degpol(f);
    if (d == 1) return 1;
  }
  f = BD(f); if (!f) return 0;
  for (i = lg(f)-1; i; i--) d -= degpol(gel(f,i));
  return gc_long(av, d == 0);
}
