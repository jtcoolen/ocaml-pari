/* Copyright (C) 2000-2004  The PARI group.

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
/*                                                                 */
/*               SUBFIELDS OF A NUMBER FIELD                       */
/*   J. Klueners and M. Pohst, J. Symb. Comp. (1996), vol. 11      */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_nfsubfields

typedef struct _poldata {
  GEN pol;
  GEN dis; /* |disc(pol)| */
  GEN roo; /* roots(pol) */
  GEN den; /* multiple of index(pol) */
} poldata;
typedef struct _primedata {
  GEN p;  /* prime */
  GEN pol; /* pol mod p, squarefree */
  GEN ff; /* factorization of pol mod p */
  GEN Z; /* cycle structure of the above [ Frobenius orbits ] */
  long lcm; /* lcm of the above */
  GEN T;  /* ffinit(p, lcm) */

  GEN fk;      /* factorization of pol over F_(p^lcm) */
  GEN firstroot; /* *[i] = index of first root of fk[i] */
  GEN interp;    /* *[i] = interpolation polynomial for fk[i]
                  * [= 1 on the first root firstroot[i], 0 on the others] */
  GEN bezoutC; /* Bezout coefficients attached to the ff[i] */
  GEN Trk;     /* used to compute traces (cf poltrace) */
} primedata;
typedef struct _blockdata {
  poldata *PD; /* data depending from pol */
  primedata *S;/* data depending from pol, p */
  GEN DATA; /* data depending from pol, p, degree, # translations [updated] */
  long N; /* deg(PD.pol) */
  long d; /* subfield degree */
  long size;/* block degree = N/d */
  long fl;
} blockdata;

static GEN print_block_system(blockdata *B, GEN Y, GEN BS);
static GEN test_block(blockdata *B, GEN L, GEN D);

/* COMBINATORIAL PART: generate potential block systems */

#define BIL 32 /* for 64bit machines also */
/* Computation of potential block systems of given size d attached to a
 * rational prime p: give a row vector of row vectors containing the
 * potential block systems of imprimitivity; a potential block system is a
 * vector of row vectors (enumeration of the roots). */
static GEN
calc_block(blockdata *B, GEN Z, GEN Y, GEN SB)
{
  long r = lg(Z), lK, i, j, t, tp, T, u, nn, lnon, lY;
  GEN K, n, non, pn, pnon, e, Yp, Zp, Zpp;
  pari_sp av0 = avma;

  if (DEBUGLEVEL>3)
  {
    err_printf("lg(Z) = %ld, lg(Y) = %ld\n", r,lg(Y));
    if (DEBUGLEVEL > 5)
    {
      err_printf("Z = %Ps\n",Z);
      err_printf("Y = %Ps\n",Y);
    }
  }
  lnon = minss(BIL, r);
  e    = new_chunk(BIL);
  n    = new_chunk(r);
  non  = new_chunk(lnon);
  pnon = new_chunk(lnon);
  pn   = new_chunk(lnon);

  Zp   = cgetg(lnon,t_VEC);
  Zpp  = cgetg(lnon,t_VEC); nn = 0;
  for (i=1; i<r; i++) { n[i] = lg(gel(Z,i))-1; nn += n[i]; }
  lY = lg(Y); Yp = cgetg(lY+1,t_VEC);
  for (j=1; j<lY; j++) gel(Yp,j) = gel(Y,j);

  {
    pari_sp av = avma;
    long k = nn / B->size;
    for (j = 1; j < r; j++)
      if (n[j] % k) break;
    if (j == r)
    {
      gel(Yp,lY) = Z;
      SB = print_block_system(B, Yp, SB);
      set_avma(av);
    }
  }
  gel(Yp,lY) = Zp;

  K = divisorsu(n[1]); lK = lg(K);
  for (i=1; i<lK; i++)
  {
    long ngcd = n[1], k = K[i], dk = B->size*k, lpn = 0;
    for (j=2; j<r; j++)
      if (n[j]%k == 0)
      {
        if (++lpn >= BIL) pari_err_OVERFLOW("calc_block");
        pn[lpn] = n[j]; pnon[lpn] = j;
        ngcd = ugcd(ngcd, n[j]);
      }
    if (dk % ngcd) continue;
    T = 1L<<lpn;
    if (lpn == r-2)
    {
      T--; /* done already above --> print_block_system */
      if (!T) continue;
    }

    if (dk == n[1])
    { /* empty subset, t = 0. Split out for clarity */
      Zp[1] = Z[1]; setlg(Zp, 2);
      for (u=1,j=2; j<r; j++) Zpp[u++] = Z[j];
      setlg(Zpp, u);
      SB = calc_block(B, Zpp, Yp, SB);
    }

    for (t = 1; t < T; t++)
    { /* loop through all nonempty subsets of [1..lpn] */
      for (nn=n[1],tp=t, u=1; u<=lpn; u++,tp>>=1)
      {
        if (tp&1) { nn += pn[u]; e[u] = 1; } else e[u] = 0;
      }
      if (dk != nn) continue;

      for (j=1; j<r; j++) non[j]=0;
      Zp[1] = Z[1];
      for (u=2,j=1; j<=lpn; j++)
        if (e[j]) { Zp[u] = Z[pnon[j]]; non[pnon[j]] = 1; u++; }
      setlg(Zp, u);
      for (u=1,j=2; j<r; j++)
        if (!non[j]) Zpp[u++] = Z[j];
      setlg(Zpp, u);
      SB = calc_block(B, Zpp, Yp, SB);
    }
  }
  return gc_const(av0, SB);
}

/* product of permutations. Put the result in perm1. */
static void
perm_mul_i(GEN perm1, GEN perm2)
{
  long i, N = lg(perm1);
  pari_sp av = avma;
  GEN perm = new_chunk(N);
  for (i=1; i<N; i++) perm[i] = perm1[perm2[i]];
  for (i=1; i<N; i++) perm1[i]= perm[i];
  set_avma(av);
}

/* cy is a cycle; compute cy^l as a permutation */
static GEN
cycle_power_to_perm(GEN perm,GEN cy,long l)
{
  long lp,i,j,b, N = lg(perm), lcy = lg(cy)-1;

  lp = l % lcy;
  for (i=1; i<N; i++) perm[i] = i;
  if (lp)
  {
    pari_sp av = avma;
    GEN p1 = new_chunk(N);
    b = cy[1];
    for (i=1; i<lcy; i++) b = (perm[b] = cy[i+1]);
    perm[b] = cy[1];
    for (i=1; i<N; i++) p1[i] = perm[i];

    for (j=2; j<=lp; j++) perm_mul_i(perm,p1);
    set_avma(av);
  }
  return perm;
}

/* image du block system D par la permutation perm */
static GEN
im_block_by_perm(GEN D,GEN perm)
{
  long i, lb = lg(D);
  GEN Dn = cgetg(lb,t_VEC);
  for (i=1; i<lb; i++) gel(Dn,i) = vecsmallpermute(perm, gel(D,i));
  return Dn;
}

static void
append(GEN D, GEN a)
{
  long i,l = lg(D), m = lg(a);
  GEN x = D + (l-1);
  for (i=1; i<m; i++) gel(x,i) = gel(a,i);
  setlg(D, l+m-1);
}

static GEN
print_block_system(blockdata *B, GEN Y, GEN SB)
{
  long i, j, l, ll, lp, u, v, ns, r = lg(Y), N = B->N;
  long *k, *n, **e, *t;
  GEN D, De, Z, cyperm, perm, VOID = cgetg(1, t_VECSMALL);

  if (DEBUGLEVEL>5) err_printf("Y = %Ps\n",Y);
  n = new_chunk(N+1);
  D = vectrunc_init(N+1);
  t = new_chunk(r+1);
  k = new_chunk(r+1);
  Z = cgetg(r+1, t_VEC);
  for (ns=0,i=1; i<r; i++)
  {
    GEN Yi = gel(Y,i);
    long ki = 0, si = lg(Yi)-1;

    for (j=1; j<=si; j++) { n[j] = lg(gel(Yi,j))-1; ki += n[j]; }
    ki /= B->size;
    De = cgetg(ki+1,t_VEC);
    for (j=1; j<=ki; j++) gel(De,j) = VOID;
    for (j=1; j<=si; j++)
    {
      GEN cy = gel(Yi,j);
      for (l=1,lp=0; l<=n[j]; l++)
      {
        lp++; if (lp > ki) lp = 1;
        gel(De,lp) = vecsmall_append(gel(De,lp), cy[l]);
      }
    }
    append(D, De);
    if (si>1 && ki>1)
    {
      GEN p1 = cgetg(si,t_VEC);
      for (j=2; j<=si; j++) p1[j-1] = Yi[j];
      ns++;
      t[ns] = si-1;
      k[ns] = ki-1;
      gel(Z,ns) = p1;
    }
  }
  if (DEBUGLEVEL>2) err_printf("\nns = %ld\n",ns);
  if (!ns) return test_block(B, SB, D);

  setlg(Z, ns+1);
  e = (long**)new_chunk(ns+1);
  for (i=1; i<=ns; i++)
  {
    e[i] = new_chunk(t[i]+1);
    for (j=1; j<=t[i]; j++) e[i][j] = 0;
  }
  cyperm= cgetg(N+1,t_VECSMALL);
  perm  = cgetg(N+1,t_VECSMALL); i = ns;
  do
  {
    pari_sp av = avma;
    for (u=1; u<=N; u++) perm[u] = u;
    for (u=1; u<=ns; u++)
      for (v=1; v<=t[u]; v++)
        perm_mul_i(perm, cycle_power_to_perm(cyperm, gmael(Z,u,v), e[u][v]));
    SB = test_block(B, SB, im_block_by_perm(D,perm));
    set_avma(av);

    /* i = 1..ns, j = 1..t[i], e[i][j] loop through 0..k[i].
     * TODO: flatten to 1-dimensional loop */
    if (++e[ns][t[ns]] > k[ns])
    {
      j = t[ns]-1;
      while (j>=1 && e[ns][j] == k[ns]) j--;
      if (j >= 1) { e[ns][j]++; for (l=j+1; l<=t[ns]; l++) e[ns][l] = 0; }
      else
      {
        i = ns-1;
        while (i>=1)
        {
          j = t[i];
          while (j>=1 && e[i][j] == k[i]) j--;
          if (j<1) i--;
          else
          {
            e[i][j]++;
            for (l=j+1; l<=t[i]; l++) e[i][l] = 0;
            for (ll=i+1; ll<=ns; ll++)
              for (l=1; l<=t[ll]; l++) e[ll][l] = 0;
            break;
          }
        }
      }
    }
  }
  while (i > 0);
  return SB;
}

/* ALGEBRAIC PART: test potential block systems */

static GEN
polsimplify(GEN x)
{
  long i,lx = lg(x);
  for (i=2; i<lx; i++)
    if (typ(gel(x,i)) == t_POL) gel(x,i) = constant_coeff(gel(x,i));
  return x;
}

/* return 0 if |g[i]| > M[i] for some i; 1 otherwise */
static long
ok_coeffs(GEN g,GEN M)
{
  long i, lg = lg(g)-1; /* g is monic, and cst term is ok */
  for (i=3; i<lg; i++)
    if (abscmpii(gel(g,i), gel(M,i)) > 0) return 0;
  return 1;
}

/* assume x in Fq, return Tr_{Fq/Fp}(x) as a t_INT */
static GEN
trace(GEN x, GEN Trq, GEN p)
{
  long i, l;
  GEN s;
  if (typ(x) == t_INT) return Fp_mul(x, gel(Trq,1), p);
  l = lg(x)-1; if (l == 1) return gen_0;
  x++; s = mulii(gel(x,1), gel(Trq,1));
  for (i=2; i<l; i++)
    s = addii(s, mulii(gel(x,i), gel(Trq,i)));
  return modii(s, p);
}

/* assume x in Fq[X], return Tr_{Fq[X]/Fp[X]}(x), varn(X) = 0 */
static GEN
poltrace(GEN x, GEN Trq, GEN p)
{
  long i,l;
  GEN y;
  if (typ(x) == t_INT || varn(x) != 0) return trace(x, Trq, p);
  y = cgetg_copy(x, &l); y[1] = x[1];
  for (i=2; i<l; i++) gel(y,i) = trace(gel(x,i),Trq,p);
  return normalizepol_lg(y, l);
}

/* Find h in Fp[X] such that h(a[i]) = listdelta[i] for all modular factors
 * ff[i], where a[i] is a fixed root of ff[i] in Fq = Z[Y]/(p,T) [namely the
 * first one in FpX_factorff_irred output]. Let f = ff[i], A the given root,
 * then h mod f is Tr_Fq/Fp ( h(A) f(X)/(X-A)f'(A) ), most of the expression
 * being precomputed. The complete h is recovered via chinese remaindering */
static GEN
chinese_retrieve_pol(GEN DATA, primedata *S, GEN listdelta)
{
  GEN interp, bezoutC, h, p = S->p, pol = FpX_red(gel(DATA,1), p);
  long i, l;
  interp = gel(DATA,9);
  bezoutC= gel(DATA,6);

  h = NULL; l = lg(interp);
  for (i=1; i<l; i++)
  { /* h(firstroot[i]) = listdelta[i] */
    GEN t = FqX_Fq_mul(gel(interp,i), gel(listdelta,i), S->T, p);
    t = poltrace(t, gel(S->Trk,i), p);
    t = FpX_mul(t, gel(bezoutC,i), p);
    h = h? FpX_add(h,t,p): t;
  }
  return FpX_rem(h, pol, p);
}

/* g in Z[X] potentially defines a subfield of Q[X]/f. It is a subfield iff A
 * (cf subfield) was a block system; then there
 * exists h in Q[X] such that f | g o h. listdelta determines h s.t f | g o h
 * in Fp[X] (cf chinese_retrieve_pol). Try to lift it; den is a
 * multiplicative bound for denominator of lift. */
static GEN
embedding(GEN g, GEN DATA, primedata *S, GEN den, GEN listdelta)
{
  GEN TR, w0_Q, w0, w1_Q, w1, wpow, h0, gp, T, q2, q, maxp, a, p = S->p;
  long rt;
  pari_sp av;

  T   = gel(DATA,1); rt = brent_kung_optpow(degpol(T), 4, 3);
  maxp= gel(DATA,7);
  gp = RgX_deriv(g); av = avma;
  w0 = chinese_retrieve_pol(DATA, S, listdelta);
  w0_Q = centermod(gmul(w0,den), p);
  h0 = FpXQ_inv(FpX_FpXQ_eval(gp,w0, T,p), T,p); /* = 1/g'(w0) mod (T,p) */
  wpow = NULL; q = sqri(p);
  for(;;)
  {/* Given g,w0,h0 in Z[x], s.t. h0.g'(w0) = 1 and g(w0) = 0 mod (T,p), find
    * [w1,h1] satisfying the same conditions mod p^2, [w1,h1] = [w0,h0] (mod p)
    * (cf. Dixon: J. Austral. Math. Soc., Series A, vol.49, 1990, p.445) */
    if (DEBUGLEVEL>1)
      err_printf("lifting embedding mod p^k = %Ps^%ld\n",S->p, Z_pval(q,S->p));

    /* w1 := w0 - h0 g(w0) mod (T,q) */
    if (wpow) a = FpX_FpXQV_eval(g,wpow, T,q);
    else      a = FpX_FpXQ_eval(g,w0, T,q); /* first time */
    /* now, a = 0 (p) */
    a = FpXQ_mul(ZX_neg(h0), ZX_Z_divexact(a, p), T,p);
    w1 = ZX_add(w0, ZX_Z_mul(a, p));

    w1_Q = centermod(ZX_Z_mul(w1, remii(den,q)), q);
    if (ZX_equal(w1_Q, w0_Q))
    {
      GEN G = is_pm1(den)? g: RgX_rescale(g,den);
      if (gequal0(RgX_RgXQ_eval(G, w1_Q, T))) break;
    }
    else if (cmpii(q,maxp) > 0)
    {
      GEN G = is_pm1(den)? g: RgX_rescale(g,den);
      if (gequal0(RgX_RgXQ_eval(G, w1_Q, T))) break;
      if (DEBUGLEVEL) err_printf("coeff too big for embedding\n");
      return NULL;
    }
    gerepileall(av, 5, &w1,&h0,&w1_Q,&q,&p);
    q2 = sqri(q);
    wpow = FpXQ_powers(w1, rt, T, q2);
    /* h0 := h0 * (2 - h0 g'(w1)) mod (T,q)
     *     = h0 + h0 * (1 - h0 g'(w1)) */
    a = FpXQ_mul(ZX_neg(h0), FpX_FpXQV_eval(gp, FpXV_red(wpow,q),T,q), T,q);
    a = ZX_Z_add_shallow(a, gen_1); /* 1 - h0 g'(w1) = 0 (p) */
    a = FpXQ_mul(h0, ZX_Z_divexact(a, p), T,p);
    h0 = ZX_add(h0, ZX_Z_mul(a, p));
    w0 = w1; w0_Q = w1_Q; p = q; q = q2;
  }
  TR = gel(DATA,5);
  if (!gequal0(TR)) w1_Q = RgX_translate(w1_Q, TR);
  return gdiv(w1_Q,den);
}

/* return U list of polynomials s.t U[i] = 1 mod fk[i] and 0 mod fk[j] for all
 * other j */
static GEN
get_bezout(GEN pol, GEN fk, GEN p)
{
  long i, l = lg(fk);
  GEN A, B, d, u, v, U = cgetg(l, t_VEC);
  for (i=1; i<l; i++)
  {
    A = gel(fk,i);
    B = FpX_div(pol, A, p);
    d = FpX_extgcd(A,B,p, &u, &v);
    if (degpol(d) > 0) pari_err_COPRIME("get_bezout",A,B);
    d = constant_coeff(d);
    if (!gequal1(d)) v = FpX_Fp_div(v, d, p);
    gel(U,i) = FpX_mul(B,v, p);
  }
  return U;
}

static GEN
init_traces(GEN ff, GEN T, GEN p)
{
  long N = degpol(T),i,j,k, r = lg(ff);
  GEN Frob = FpX_matFrobenius(T,p);
  GEN y,p1,p2,Trk,pow,pow1;

  k = degpol(gel(ff,r-1)); /* largest degree in modular factorization */
  pow = cgetg(k+1, t_VEC);
  gel(pow,1) = gen_0; /* dummy */
  gel(pow,2) = Frob;
  pow1= cgetg(k+1, t_VEC); /* 1st line */
  for (i=3; i<=k; i++)
    gel(pow,i) = FpM_mul(gel(pow,i-1), Frob, p);
  gel(pow1,1) = gen_0; /* dummy */
  for (i=2; i<=k; i++)
  {
    p1 = cgetg(N+1, t_VEC);
    gel(pow1,i) = p1; p2 = gel(pow,i);
    for (j=1; j<=N; j++) gel(p1,j) = gcoeff(p2,1,j);
  }

  /* Trk[i] = line 1 of x -> x + x^p + ... + x^{p^(i-1)} */
  Trk = pow; /* re-use (destroy) pow */
  gel(Trk,1) = vec_ei(N,1);
  for (i=2; i<=k; i++)
    gel(Trk,i) = gadd(gel(Trk,i-1), gel(pow1,i));
  y = cgetg(r, t_VEC);
  for (i=1; i<r; i++) y[i] = Trk[degpol(gel(ff,i))];
  return y;
}

static void
init_primedata(primedata *S)
{
  long i, j, l, lff = lg(S->ff), v = fetch_var(), N = degpol(S->pol);
  GEN T, p = S->p;

  if (S->lcm == degpol(gel(S->ff,lff-1)))
  {
    T = leafcopy(gel(S->ff,lff-1));
    setvarn(T, v);
  }
  else
    T = init_Fq(p, S->lcm, v);
  S->T = T;
  S->firstroot = cgetg(lff, t_VECSMALL);
  S->interp = cgetg(lff, t_VEC);
  S->fk = cgetg(N+1, t_VEC);
  for (l=1,j=1; j<lff; j++)
  { /* compute roots and fix ordering (Frobenius cycles) */
    GEN F = gel(S->ff, j), deg1 = FpX_factorff_irred(F, T,p);
    GEN H = gel(deg1,1), a = Fq_neg(constant_coeff(H), T,p);
    GEN Q = FqX_div(F, H, T,p);
    GEN q = Fq_inv(FqX_eval(Q, a, T,p), T,p);
    gel(S->interp,j) = FqX_Fq_mul(Q, q, T,p); /* = 1 at a, 0 at other roots */
    S->firstroot[j] = l;
    for (i=1; i<lg(deg1); i++,l++) gel(S->fk, l) = gel(deg1, i);
  }
  S->Trk     = init_traces(S->ff, T,p);
  S->bezoutC = get_bezout(S->pol, S->ff, p);
}

static int
choose_prime(primedata *S, GEN pol)
{
  long i, j, k, r, lcm, oldr, K, N = degpol(pol);
  ulong p, pp;
  GEN Z, ff, n, oldn;
  pari_sp av;
  forprime_t T;

  u_forprime_init(&T, (N*N) >> 2, ULONG_MAX);
  oldr = S->lcm = LONG_MAX;
  S->ff = oldn = NULL; pp = 0; /* gcc -Wall */
  av = avma; K = N + 10;
  for(k = 1; k < K || !S->ff; k++,set_avma(av))
  {
    GEN Tp;
    if (k > 5 * N) return 0;
    do
    {
      p = u_forprime_next(&T);
      Tp = ZX_to_Flx(pol, p);
    }
    while (!Flx_is_squarefree(Tp, p));
    ff = gel(Flx_factor(Tp, p), 1);
    r = lg(ff)-1;
    if (r == N || r >= BIL) continue;

    n = cgetg(r+1, t_VECSMALL); lcm = n[1] = degpol(gel(ff,1));
    for (j=2; j<=r; j++) { n[j] = degpol(gel(ff,j)); lcm = ulcm(lcm, n[j]); }
    if (r > oldr || (r == oldr && (lcm <= S->lcm || S->lcm > 2*N)))
      continue;
    if (DEBUGLEVEL) err_printf("p = %lu,\tlcm = %ld,\torbits: %Ps\n",p,lcm,n);

    pp = p;
    oldr = r;
    oldn = n;
    S->ff = ff;
    S->lcm = lcm; if (r == 1) break;
    av = avma;
  }
  if (oldr > 6) return 0;
  if (DEBUGLEVEL) err_printf("Chosen prime: p = %ld\n", pp);
  FlxV_to_ZXV_inplace(S->ff);
  S->p  = utoipos(pp);
  S->pol = FpX_red(pol, S->p); init_primedata(S);
  n = oldn; r = lg(n); S->Z = Z = cgetg(r,t_VEC);
  for (k=0,i=1; i<r; i++)
  {
    GEN t = cgetg(n[i]+1, t_VECSMALL); gel(Z,i) = t;
    for (j=1; j<=n[i]; j++) t[j] = ++k;
  }
  return 1;
}

/* maxroot t_REAL */
static GEN
bound_for_coeff(long m, GEN rr, GEN *maxroot)
{
  long i,r1, lrr=lg(rr);
  GEN p1,b1,b2,B,M, C = matpascal(m-1);

  for (r1=1; r1 < lrr; r1++)
    if (typ(gel(rr,r1)) != t_REAL) break;
  r1--;

  rr = gabs(rr,0); *maxroot = vecmax(rr);
  for (i=1; i<lrr; i++)
    if (gcmp(gel(rr,i), gen_1) < 0) gel(rr,i) = gen_1;
  for (b1=gen_1,i=1; i<=r1; i++) b1 = gmul(b1, gel(rr,i));
  for (b2=gen_1    ; i<lrr; i++) b2 = gmul(b2, gel(rr,i));
  B = gmul(b1, gsqr(b2)); /* Mahler measure */
  M = cgetg(m+2, t_VEC); gel(M,1) = gel(M,2) = gen_0; /* unused */
  for (i=1; i<m; i++)
  {
    p1 = gadd(gmul(gcoeff(C, m, i+1), B),/* binom(m-1, i)   */
              gcoeff(C, m, i));          /* binom(m-1, i-1) */
    gel(M,i+2) = ceil_safe(p1);
  }
  return M;
}

/* d = requested degree for subfield. Return DATA, valid for given pol, S and d
 * If DATA != NULL, translate pol [ --> pol(X+1) ] and update DATA
 * 1: polynomial pol
 * 2: p^e (for Hensel lifts) such that p^e > max(M),
 * 3: Hensel lift to precision p^e of DATA[4]
 * 4: roots of pol in F_(p^S->lcm),
 * 5: number of polynomial changes (translations)
 * 6: Bezout coefficients attached to the S->ff[i]
 * 7: Hadamard bound for coefficients of h(x) such that g o h = 0 mod pol.
 * 8: bound M for polynomials defining subfields x PD->den
 * 9: *[i] = interpolation polynomial for S->ff[i] [= 1 on the first root
      S->firstroot[i], 0 on the others] */
static void
compute_data(blockdata *B)
{
  GEN ffL, roo, pe, p1, p2, fk, fhk, MM, maxroot, pol;
  primedata *S = B->S;
  GEN p = S->p, T = S->T, ff = S->ff, DATA = B->DATA;
  long i, j, l, e, N, lff = lg(ff);

  if (DEBUGLEVEL>1) err_printf("Entering compute_data()\n\n");
  pol = B->PD->pol; N = degpol(pol);
  roo = B->PD->roo;
  if (DATA)
  {
    GEN Xm1 = gsub(pol_x(varn(pol)), gen_1);
    GEN TR = addiu(gel(DATA,5), 1);
    GEN mTR = negi(TR), interp, bezoutC;

    if (DEBUGLEVEL>1) err_printf("... update (translate) an existing DATA\n\n");

    gel(DATA,5) = TR;
    pol = RgX_translate(gel(DATA,1), gen_m1);
    p1 = cgetg_copy(roo, &l);
    for (i=1; i<l; i++) gel(p1,i) = gadd(TR, gel(roo,i));
    roo = p1;

    fk = gel(DATA,4); l = lg(fk);
    for (i=1; i<l; i++) gel(fk,i) = gsub(Xm1, gel(fk,i));

    bezoutC = gel(DATA,6); l = lg(bezoutC);
    interp  = gel(DATA,9);
    for (i=1; i<l; i++)
    {
      if (degpol(gel(interp,i)) > 0) /* do not turn pol_1(0) into gen_1 */
      {
        p1 = RgX_translate(gel(interp,i), gen_m1);
        gel(interp,i) = FpXX_red(p1, p);
      }
      if (degpol(gel(bezoutC,i)) > 0)
      {
        p1 = RgX_translate(gel(bezoutC,i), gen_m1);
        gel(bezoutC,i) = FpXX_red(p1, p);
      }
    }
    ff = cgetg(lff, t_VEC); /* copy, do not overwrite! */
    for (i=1; i<lff; i++)
      gel(ff,i) = FpX_red(RgX_translate(gel(S->ff,i), mTR), p);
  }
  else
  {
    DATA = cgetg(10,t_VEC);
    fk = S->fk;
    gel(DATA,5) = gen_0;
    gel(DATA,6) = leafcopy(S->bezoutC);
    gel(DATA,9) = leafcopy(S->interp);
  }
  gel(DATA,1) = pol;
  MM = gmul2n(bound_for_coeff(B->d, roo, &maxroot), 1);
  gel(DATA,8) = MM;
  e = logintall(shifti(vecmax(MM),20), p, &pe); /* overlift 2^20 [d-1 test] */
  gel(DATA,2) = pe;
  gel(DATA,4) = roots_from_deg1(fk);

  /* compute fhk = ZpX_liftfact(pol,fk,T,p,e,pe) in 2 steps
   * 1) lift in Zp to precision p^e */
  ffL = ZpX_liftfact(pol, ff, pe, p, e);
  fhk = NULL;
  for (l=i=1; i<lff; i++)
  { /* 2) lift factorization of ff[i] in Qp[X] / T */
    GEN F, L = gel(ffL,i);
    long di = degpol(L);
    F = cgetg(di+1, t_VEC);
    for (j=1; j<=di; j++) F[j] = fk[l++];
    L = ZqX_liftfact(L, F, T, pe, p, e);
    fhk = fhk? shallowconcat(fhk, L): L;
  }
  gel(DATA,3) = roots_from_deg1(fhk);

  p1 = mulur(N, powruhalf(utor(N-1,DEFAULTPREC), N-1));
  p2 = powru(maxroot, B->size + N*(N-1)/2);
  p1 = divrr(mulrr(p1,p2), gsqrt(B->PD->dis,DEFAULTPREC));
  gel(DATA,7) = mulii(shifti(ceil_safe(p1), 1), B->PD->den);

  if (DEBUGLEVEL>1) {
    err_printf("f = %Ps\n",DATA[1]);
    err_printf("p = %Ps, lift to p^%ld\n", p, e);
    err_printf("2 * Hadamard bound * ind = %Ps\n",DATA[7]);
    err_printf("2 * M = %Ps\n",DATA[8]);
  }
  if (B->DATA) { DATA = gclone(DATA); if (isclone(B->DATA)) gunclone(B->DATA); }
  B->DATA = DATA;
}

/* g = polynomial, h = embedding. Return [[g,h]] */
static GEN
_subfield(GEN g, GEN h) { return mkvec(mkvec2(g,h)); }

/* Return a subfield, gen_0 [ change p ] or NULL [ not a subfield ] */
static GEN
subfield(GEN A, blockdata *B)
{
  long N, i, j, d, lf, m = lg(A)-1;
  GEN M, pe, pol, fhk, g, e, d_1_term, delta, listdelta, whichdelta;
  GEN T = B->S->T, p = B->S->p, firstroot = B->S->firstroot;

  pol= gel(B->DATA,1); N = degpol(pol); d = N/m; /* m | N */
  pe = gel(B->DATA,2);
  fhk= gel(B->DATA,3);
  M  = gel(B->DATA,8);

  delta = cgetg(m+1,t_VEC);
  whichdelta = cgetg(N+1, t_VECSMALL);
  d_1_term = gen_0;
  for (i=1; i<=m; i++)
  {
    GEN Ai = gel(A,i), p1 = gel(fhk,Ai[1]);
    for (j=2; j<=d; j++)
      p1 = Fq_mul(p1, gel(fhk,Ai[j]), T, pe);
    gel(delta,i) = p1;
    if (DEBUGLEVEL>5) err_printf("delta[%ld] = %Ps\n",i,p1);
    /* g = prod (X - delta[i])
     * if g o h = 0 (pol), we'll have h(Ai[j]) = delta[i] for all j */
    /* fk[k] belongs to block number whichdelta[k] */
    for (j=1; j<=d; j++) whichdelta[Ai[j]] = i;
    if (typ(p1) == t_POL) p1 = constant_coeff(p1);
    d_1_term = addii(d_1_term, p1);
  }
  d_1_term = centermod(d_1_term, pe); /* Tr(g) */
  if (abscmpii(d_1_term, gel(M,3)) > 0) {
    if (DEBUGLEVEL>1) err_printf("d-1 test failed\n");
    return NULL;
  }
  g = FqV_roots_to_pol(delta, T, pe, 0);
  g = centermod(polsimplify(g), pe); /* assume g in Z[X] */
  if (!ok_coeffs(g,M)) {
    if (DEBUGLEVEL>2) err_printf("pol. found = %Ps\n",g);
    if (DEBUGLEVEL>1) err_printf("coeff too big for pol g(x)\n");
    return NULL;
  }
  if (!FpX_is_squarefree(g, p)) {
    if (DEBUGLEVEL>2) err_printf("pol. found = %Ps\n",g);
    if (DEBUGLEVEL>1) err_printf("changing f(x): p divides disc(g)\n");
    compute_data(B);
    return subfield(A, B);
  }

  lf = lg(firstroot); listdelta = cgetg(lf, t_VEC);
  for (i=1; i<lf; i++) listdelta[i] = delta[whichdelta[firstroot[i]]];
  if (DEBUGLEVEL) err_printf("candidate = %Ps\n", g);
  e = embedding(g, B->DATA, B->S, B->PD->den, listdelta);
  if (!e) return NULL;
  if (DEBUGLEVEL) err_printf("... OK!\n");
  return B->fl==1? mkvec(g):_subfield(g, e);
}

/* L list of current subfields, test whether potential block D is a block,
 * if so, append corresponding subfield */
static GEN
test_block(blockdata *B, GEN L, GEN D)
{
  pari_sp av = avma;
  GEN sub = subfield(D, B);
  if (sub) {
    GEN old = L;
    L = gclone( L? shallowconcat(L, sub): sub );
    guncloneNULL(old);
  }
  return gc_const(av,L);
}

/* subfields of degree d */
static GEN
subfields_of_given_degree(blockdata *B)
{
  pari_sp av = avma;
  GEN L;

  if (DEBUGLEVEL) err_printf("\n* Look for subfields of degree %ld\n\n", B->d);
  B->DATA = NULL; compute_data(B);
  L = calc_block(B, B->S->Z, cgetg(1,t_VEC), NULL);
  if (DEBUGLEVEL>9)
    err_printf("\nSubfields of degree %ld: %Ps\n", B->d, L? L: cgetg(1,t_VEC));
  if (isclone(B->DATA)) gunclone(B->DATA);
  return gc_const(av,L);
}

static void
setvarn2(GEN t, long v) { setvarn(gel(t,1),v); setvarn(gel(t,2),v); }
static GEN
fix_var(GEN x, long v, long fl)
{
  long i, l = lg(x);
  if (!v) return x;
  if (fl)
    for (i = 1; i < l; i++) setvarn(gel(x,i), v);
  else
    for (i = 1; i < l; i++) setvarn2(gel(x,i), v);
  return x;
}

static void
subfields_poldata(GEN nf, GEN T, poldata *PD)
{
  GEN L, dis;

  PD->pol = T;
  if (nf)
  {
    PD->den = nf_get_zkden(nf);
    PD->roo = nf_get_roots(nf);
    PD->dis = mulii(absi_shallow(nf_get_disc(nf)), sqri(nf_get_index(nf)));
  }
  else
  {
    PD->den = initgaloisborne(T,NULL,nbits2prec(bit_accuracy(ZX_max_lg(T))), &L,NULL,&dis);
    PD->roo = L;
    PD->dis = absi_shallow(dis);
  }
}

static GEN nfsubfields_fa(GEN nf, long d, long fl);
static GEN
subfieldsall(GEN nf0, long fl)
{
  pari_sp av = avma;
  long N, ld, i, v;
  GEN nf, G, T, dg, LSB, NLSB;
  poldata PD;
  primedata S;
  blockdata B;

  /* much easier if nf is Galois (WSS) */
  G = galoisinit(nf0, NULL);
  T = get_nfpol(nf0, &nf);
  if (G != gen_0)
  {
    GEN L, S;
    long l;
    L = lift_shallow( galoissubfields(G, fl, varn(T)) );
    l = lg(L); S = cgetg(l, t_VECSMALL);
    for (i=1; i<l; i++) S[i] = lg(fl==1? gel(L,i): gmael(L,i,1));
    return gerepilecopy(av, vecpermute(L, vecsmall_indexsort(S)));
  }
  v = varn(T); N = degpol(T);
  dg = divisorsu(N); ld = lg(dg)-1;
  LSB = fl==1 ? mkvec(pol_x(v)): _subfield(pol_x(v), pol_0(v));
  if (ld <= 2)
  {
    if (ld == 2)
      LSB = shallowconcat(LSB, fl==1? mkvec(T): _subfield(T, pol_x(v)));
    return gerepilecopy(av, LSB);
  }
  if (varn(T) != 0) { T = leafcopy(T); setvarn(T, 0); }
  if (!choose_prime(&S, T)) { set_avma(av); return nfsubfields_fa(nf0, 0, fl); }
  subfields_poldata(nf, T, &PD);

  if (DEBUGLEVEL) err_printf("\n***** Entering subfields\n\npol = %Ps\n",T);
  B.PD = &PD;
  B.S  = &S;
  B.N  = N;
  B.fl = fl;
  for (i=ld-1; i>1; i--)
  {
    B.size  = dg[i];
    B.d = N / B.size;
    NLSB = subfields_of_given_degree(&B);
    if (NLSB) { LSB = gconcat(LSB, NLSB); gunclone(NLSB); }
  }
  (void)delete_var(); /* from init_primedata() */
  LSB = shallowconcat(LSB, fl==1? mkvec(T):_subfield(T, pol_x(0)));
  if (DEBUGLEVEL) err_printf("\n***** Leaving subfields\n\n");
  return fix_var(gerepilecopy(av, LSB), v, fl);
}

GEN
nfsubfields0(GEN nf0, long d, long fl)
{
  pari_sp av = avma;
  long N, v0;
  GEN nf, LSB, T, G;
  poldata PD;
  primedata S;
  blockdata B;
  if (fl<0 || fl>1) pari_err_FLAG("nfsubfields");
  if (typ(nf0)==t_VEC && lg(nf0)==3) return nfsubfields_fa(nf0, d, fl);
  if (!d) return subfieldsall(nf0, fl);

  /* treat trivial cases */
  T = get_nfpol(nf0, &nf); v0 = varn(T); N = degpol(T);
  RgX_check_ZX(T,"nfsubfields");
  if (d == N)
    return gerepilecopy(av, fl==1 ? mkvec(T) : _subfield(T, pol_x(v0)));
  if (d == 1)
    return gerepilecopy(av, fl==1 ? mkvec(pol_x(v0)) : _subfield(pol_x(v0), zeropol(v0)));
  if (d < 1 || d > N || N % d) return cgetg(1,t_VEC);

  /* much easier if nf is Galois (WSS) */
  G = galoisinit(nf0, NULL);
  if (G != gen_0)
  { /* Bingo */
    GEN L = galoissubgroups(G), F;
    long k,i, l = lg(L), o = N/d;
    F = cgetg(l, t_VEC);
    k = 1;
    for (i=1; i<l; i++)
    {
      GEN H = gel(L,i);
      if (group_order(H) == o)
        gel(F,k++) = lift_shallow(galoisfixedfield(G, gel(H,1), fl, v0));
    }
    setlg(F, k);
    return gerepilecopy(av, F);
  }
  if (varn(T) != 0) { T = leafcopy(T); setvarn(T, 0); }
  if (!choose_prime(&S, T)) { set_avma(av); return nfsubfields_fa(nf0, d, fl); }
  subfields_poldata(nf, T, &PD);
  B.PD = &PD;
  B.S  = &S;
  B.N  = N;
  B.d  = d;
  B.size = N/d;
  B.fl = fl;
  LSB = subfields_of_given_degree(&B);
  (void)delete_var(); /* from init_primedata */
  set_avma(av);
  if (!LSB) return cgetg(1, t_VEC);
  G = gcopy(LSB); gunclone(LSB);
  return fix_var(G, v0, fl);
}

GEN
nfsubfields(GEN nf0, long d)
{ return nfsubfields0(nf0, d, 0); }

/******************************/
/*                            */
/*    Maximal CM subfield     */
/*     Aurel Page (2019)      */
/*                            */
/******************************/

/* ero: maximum exponent+1 of roots of pol */
static GEN
try_subfield_generator(GEN pol, GEN v, long e, long p, long ero, long fl)
{
  GEN a, P, Q;
  long d, bound, i, B, bi, ed;

  a = gtopolyrev(v, varn(pol));
  P = Flxq_charpoly(ZX_to_Flx(a,p), ZX_to_Flx(pol,p), p);
  Flx_ispower(P, e, p, &Q);
  if (!Flx_is_squarefree(Q,p)) return NULL;
  d = degpol(pol)/e;
  B = 0;
  for (i=1; i<lg(v); i++)
  {
    bi = (i-1)*ero + expi(gel(v,i));
    if (bi > B) B = bi;
  }
  ed = expu(d);
  B += ed+1;
  bound = 0;
  for (i=0; 2*i<=d; i++)
  {
    if (!i) bi = d*B;
    else    bi = (d-i)*B + i*(3+ed-expu(i));
    if (bi > bound) bound = bi;
  }
  Q = ZXQ_minpoly(a,pol,d,bound);
  return fl==1? Q: mkvec2(Q, a);
}

/* subfield sub of nf of degree d assuming:
   - V is contained in sub
   - V is not contained in a proper subfield of sub
   ero: maximum exponent+1 of roots of pol
   output as nfsubfields:
   - pair [g,h] where g absolute equation for the  subfield and h expresses
   - one of the roots of g in terms of the generator of nf
*/
static GEN
subfield_generator(GEN pol, GEN V, long d, long ero, long fl)
{
  long p, i, e, vp = varn(pol);
  GEN a = NULL, v = cgetg(lg(V),t_COL), B;

  if (d==1) return fl ? pol_x(vp): mkvec2(pol_x(vp), pol_0(vp));
  e = degpol(pol)/d;
  p = 1009;
  for (i=1; i<lg(V); i++)
  {
    a = try_subfield_generator(pol, gel(V,i), e, p, ero, fl);
    if (a) return a;
    p = unextprime(p+1);
  }
  B = stoi(10);
  while(1)
  {
    for (i=1; i<lg(v); i++) gel(v,i) = randomi(B);
    a = try_subfield_generator(pol, QM_QC_mul(V,v), e, p, ero, fl);
    if (a) return a;
    p = unextprime(p+1);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN
RgXY_to_RgC(GEN P, long dx, long dy)
{
  GEN res, c;
  long i, j, k, d = degpol(P);
  if (d > dy) pari_err_BUG("RgXY_to_RgC [incorrect degree]");
  res = cgetg((dx+1)*(dy+1)+1, t_COL);
  k = 1;
  for (i=0; i<=d; i++)
  {
    c = gel(P,i+2);
    if (typ(c)==t_POL)
    {
      long dc = degpol(c);
      if (dc > dx) pari_err_BUG("RgXY_to_RgC [incorrect degree]");
      for (j=0; j<=dc; j++)
        gel(res,k++) = gel(c,j+2);
    } else
    {
      gel(res,k++) = c; j=1;
    }
    for (  ; j<=dx; j++)
      gel(res,k++) = gen_0;
  }
  for(  ; i<=dy; i++)
    for (j=0; j<=dx; j++)
      gel(res,k++) = gen_0;
  return res;
}

/* lambda: t_VEC of t_INT; 0 means ignore this factor */
static GEN
twoembequation(GEN pol, GEN fa, GEN lambda)
{
  GEN m, vpolx, poly;
  long i,j, lfa = lg(fa), dx = degpol(pol);
  long vx = varn(pol), vy = varn(gel(fa,1)); /* vx < vy ! */

  if (varncmp(vx,vy) <= 0) pari_err_BUG("twoembequation [incorrect variable priorities]");

  lambda = shallowcopy(lambda);
  fa = shallowcopy(fa);
  j = 1;
  for (i=1; i<lfa; i++)
    if (signe(gel(lambda,i)))
    {
      gel(lambda,j) = negi(gel(lambda,i));
      gel(fa,j) = gel(fa,i);
      j++;
    }
  setlg(lambda, j);
  setlg(fa, j); lfa = j;

  vpolx = ZXQ_powers(pol_x(vx),dx-1,pol);
  m = cgetg(dx+1, t_MAT);
  for (j=1; j <= dx; j++)
    gel(m,j) = cgetg(lfa, t_COL);
  for(i=1; i<lfa; i++)
  {
    long dy = degpol(gel(fa,i));
    poly = pol_1(vy);
    for (j=1; j <= dx; j++)
    {
      gcoeff(m,i,j) = RgXY_to_RgC(gadd(ZX_Z_mul(gel(vpolx,j),gel(lambda,i)),poly), dx, dy);
      poly = RgXQX_rem(RgX_shift(poly,1), gel(fa,i), pol);
    }
  }
  for(j=1; j<=dx; j++) gel(m,j) = shallowconcat1(gel(m,j));
  return QM_ker(m);
}

static void
subfields_cleanup(GEN* nf, GEN* pol, long* n, GEN* fa)
{
  *fa = NULL;
  if (typ(*nf) != t_VEC && typ(*nf) != t_POL) pari_err_TYPE("subfields_cleanup", *nf);
  if (typ(*nf) == t_VEC && lg(*nf) == 3)
  {
    *fa = gel(*nf,2);
    *nf = gel(*nf,1);
    if (typ(*fa)!=t_MAT || lg(*fa)!=3)
      pari_err_TYPE("subfields_cleanup [fa should be a factorisation matrix]", *fa);
  }
  if (typ(*nf) == t_POL)
  {
    *pol = *nf;
    *nf = NULL;
    if (!RgX_is_ZX(*pol)) pari_err_TYPE("subfields_cleanup [not integral]", *pol);
    if (!equali1(leading_coeff(*pol))) pari_err_TYPE("subfields_cleanup [not monic]", *pol);
    *n = degpol(*pol);
    if (*n<=0) pari_err_TYPE("subfields_cleanup [constant polynomial]", *pol);
  }
  else
  {
    *nf = checknf(*nf);
    *pol = nf_get_pol(*nf);
    *n = degpol(*pol);
  }
  if(*fa)
  {
    long v = varn(*pol);
    GEN o = gcoeff(*fa,1,1);
    if (varncmp(varn(o),v) >= 0) pari_err_PRIORITY("nfsubfields_fa", o, "<=", v);
  }
}

static GEN
rootsuptoconj(GEN pol, long prec)
{
  GEN ro;
  long n, i;
  ro = roots(pol,prec);
  n = lg(ro)-1;
  for (i=1; i<=n/2; i++)
    gel(ro,i) = gel(ro,2*i-1);
  setlg(ro,n/2+1);
  return ro;
}
static GEN
cmsubfield_get_roots(GEN pol, GEN nf, long n, long* r2, long *prec)
{
  GEN ro;
  if (nf)
  {
    if (nf_get_r1(nf)) return NULL;
    *r2 = nf_get_r2(nf);
    *prec = nf_get_prec(nf);
    ro = nf_get_roots(nf);
  }
  else
  {
    if (n%2 || sturm(pol)) return NULL;
    *r2 = n/2;
    *prec = MEDDEFAULTPREC;
    ro = rootsuptoconj(pol, *prec);
  }
  return ro;
}

static GEN
subfields_get_fa(GEN pol, GEN nf, GEN fa)
{
  if (!fa)
  {
    GEN poly = shallowcopy(pol);
    setvarn(poly, fetch_var_higher());
    fa = nffactor(nf? nf: pol, poly);
  }
  return liftpol_shallow(gel(fa,1));
}

static long
subfields_get_ero(GEN pol, GEN nf)
{
  return 1 + gexpo(nf? nf_get_roots(nf):
                       QX_complex_roots(pol, LOWDEFAULTPREC));
}

static GEN
try_imag(GEN x, GEN c, GEN pol, long v, ulong p, GEN emb, GEN galpol, long fl)
{
  GEN a = Q_primpart(RgX_sub(RgX_RgXQ_eval(x,c,pol),x));
  if (Flx_is_squarefree(Flxq_charpoly(ZX_to_Flx(a,p),ZX_to_Flx(pol,p),p),p))
  {
    pol = ZXQ_charpoly(a, pol, v);
    return fl ? pol : mkvec2(pol, RgX_RgXQ_eval(a, emb, galpol));
  }
  return NULL;
}

static GEN
galoissubfieldcm(GEN G, long fl)
{
  pari_sp av = avma;
  GEN c, H, elts, g, Hset, c2, gene, sub, pol, emb, a, galpol, B, b;
  long n, i, j, nH, ind, v, d;
  ulong p = 1009;

  galpol = gal_get_pol(G);
  n = degpol(galpol);
  v = varn(galpol);
  c = galois_get_conj(G);
  /* compute the list of c*g*c*g^(-1) : product of all pairs of conjugations
   * maximal CM subfield is the field fixed by those elements, if c does not
   * belong to the group they generate */
  checkgroup(G, &elts);
  elts = gen_sort_shallow(elts,(void*)vecsmall_lexcmp,cmp_nodata);
  H = vecsmall_ei(n,1); /* indices of elements of H */
  Hset = zero_F2v(n);
  F2v_set(Hset,1);
  nH = 1;
  for (i=2; i<=n; i++)
  {
    g = gel(elts,i);
    c2 = perm_mul(c,perm_conj(g,c));
    if (!F2v_coeff(Hset,c2[1]))
    {
      nH++;
      H[nH] = c2[1];
      F2v_set(Hset,c2[1]);
    }
  }
  /* group generated */
  gene = gcopy(H);
  setlg(gene,nH+1);
  i = 1; /* last element that has been multiplied by the generators */
  while (i < nH)
  {
    for (j=1; j<lg(gene); j++)
    {
      g = gel(elts,gene[j]);
      ind = g[H[i]]; /* index of the product */
      if (!F2v_coeff(Hset,ind))
      {
        nH++;
        if (ind==c[1] || 2*nH>n) return gc_const(av, gen_0);
        H[nH] = ind;
        F2v_set(Hset,ind);
      }
    }
    i++;
  }
  H = cgetg(lg(gene), t_VEC);
  for (i=1; i<lg(H); i++)
    gel(H,i) = gel(elts,gene[i]);
  sub = galoisfixedfield(G, H, 0, -1);

  /* compute a totally imaginary generator */
  pol = gel(sub,1);
  emb = liftpol_shallow(gel(sub,2));
  d = degpol(pol);
  if (!(ZX_deflate_order(pol)%2) && sturm(RgX_deflate(pol,2))==d/2)
  {
    setvarn(pol,v);
    return fl==1 ? pol: mkvec2(pol,emb);
  }

  /* compute action of c on the subfield from that on the large field */
  c = galoispermtopol(G,c);
  if (d<n)
  {
    GEN M = cgetg(d+1,t_MAT), contc, contM;
    gel(M,1) = col_ei(n,1); a = pol_1(v);
    for (i=2; i<=d; i++)
    {
      a = RgX_rem(QX_mul(a,emb), galpol);
      gel(M,i) = RgX_to_RgC(a,n);
    }
    c = RgX_RgXQ_eval(emb,c,galpol);
    c = Q_primitive_part(c,&contc);
    c = RgX_to_RgC(c,n);
    M = Q_primitive_part(M,&contM);
    c = RgM_RgC_invimage(M,c);
    if (contc)
    {
      if (contM) contc = gdiv(contc,contM);
      c = RgV_Rg_mul(c, contc);
    }
    else if (contM) c = RgV_Rg_mul(c, ginv(contM));
    c = RgV_to_RgX(c, v);
  }

  /* search for a generator of the form c(b)-b */
  for (i=1; i<d; i++)
  {
    a = try_imag(pol_xn(i,v),c,pol,v,p,emb,galpol,fl);
    if (a) return a;
    p = unextprime(p+1);
  }
  B = stoi(10);
  b = pol_xn(d-1,v);
  while(1)
  {
    for (i=2; i<lg(b); i++) gel(b,i) = randomi(B);
    a = try_imag(b,c,pol,v,p,emb,galpol,fl);
    if (a) return a;
    p = unextprime(p+1);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN
quadsubfieldcm(GEN pol, long fl)
{
  GEN a = gel(pol,3), b = gel(pol,2), d, P;
  long v = varn(pol);
  if (mpodd(a))
  { b = mului(4, b); d = gen_2; }
  else
  { a = divis(a,2);  d = gen_1; }
  P = deg2pol_shallow(gen_1, gen_0, subii(b, sqri(a)), v);
  return fl==1 ? P: mkvec2(P, deg1pol_shallow(d,a,v));
}

GEN
nfsubfieldscm(GEN nf, long fl)
{
  pari_sp av = avma;
  GEN fa, lambda, V, res, ro, a, aa, ev, minev, pol, G;
  long i, j, n, r2, minj=0, prec, emax, emin, e, precbound, ero;

  subfields_cleanup(&nf, &pol, &n, &fa);
  ro = cmsubfield_get_roots(pol, nf, n, &r2, &prec);
  if (!ro) return gc_const(av, gen_0);
  /* now r2 == 2*n */

  if (n==2) return gerepilecopy(av, quadsubfieldcm(pol, fl));
  G = galoisinit(nf? nf: pol, NULL);
  if (G != gen_0) return gerepilecopy(av, galoissubfieldcm(G, fl));

  ero = 0;
  for (i=1; i<lg(ro); i++)
  {
    e = 1+gexpo(gel(ro,i));
    if (e > ero) ero = e;
  }
  ero++;
  fa = subfields_get_fa(pol, nf, fa);

  emax = 1;
  emin = -1;
  for (i=1; i<lg(ro); i++)
    for (j=i+1; j<lg(ro); j++)
    {
      e = gexpo(gsub(gel(ro,i),gel(ro,j)));
      if (e > emax) emax = e;
      if (e < emin) emin = e;
    }
  precbound = n*(emax-emin) + gexpo(fa) + n*n + 5;
  precbound = 3 + precbound/BITS_IN_LONG;
  if (prec < precbound)
  {
    prec = precbound;
    ro = rootsuptoconj(pol, prec);
  }

  lambda = zerovec(lg(fa)-1);
  for (i=1; i<=r2; i++)
  {
    a = gel(ro,i);
    aa = conj_i(a);
    for (j=1; j<lg(fa); j++)
    {
      ev = cxnorm(poleval(poleval(gel(fa,j),aa),a));
      if (j==1 || cmprr(minev,ev)>0) { minj = j; minev = ev; }
    }
    gel(lambda,minj) = gen_m1;
  }

  V = twoembequation(pol, fa, lambda);
  if (lg(V)==1) { delete_var(); return gc_const(av, gen_0); }
  res = subfield_generator(pol, V, 2*(lg(V)-1), ero, fl);
  delete_var();
  return gerepilecopy(av, res);
}

static int
field_is_contained(GEN V, GEN W, int strict)
{
  GEN VW;
  ulong p = 1073741827;
  /* distinct overfield must have different dimension */
  if (strict && lg(V) == lg(W)) return 0;
  /* dimension of overfield must be multiple */
  if ((lg(W)-1) % (lg(V)-1)) return 0;
  VW = shallowconcat(V,W);
  if (Flm_rank(ZM_to_Flm(VW,p),p) > lg(W)-1) return 0;
  return ZM_rank(VW) == lg(W)-1;
}

/***********************************************/
/*                                             */
/*    Maximal, generating, all subfields       */
/*             Aurel Page (2019)               */
/*     after van Hoeij, Klueners, Novocin      */
/*  Journal of Symbolic Computation 52 (2013)  */
/*                                             */
/***********************************************/

const long subf_MAXIMAL = 1; /* return the maximal subfields */
const long subf_GENERATING = 2; /* return the generating subfields */
static GEN
maxgen_subfields(GEN pol, GEN fa, long flag)
{
  pari_sp av = avma;
  GEN principal, ismax, isgene, Lmax = NULL, Lgene, res, V, W, W1;
  long i, i2, j, flmax, flgene, nbmax = 0, nbgene = 0;

  if (!flag) return cgetg(1,t_VEC);
  flmax = (flag & subf_MAXIMAL)!=0;
  flgene = (flag & subf_GENERATING)!=0;

  /* compute principal subfields */
  principal = cgetg(lg(fa),t_VEC);
  for (i=1; i<lg(fa); i++)
    gel(principal,i) = twoembequation(pol, fa, vec_ei(lg(fa)-1,i));
  principal = gen_sort_uniq(principal, (void*)&cmp_universal, &cmp_nodata);
  /* remove nf and duplicates (sort_uniq possibly not enough) */
  i2 = 1;
  for (i=1; i<lg(principal)-1; i++)
  {
    long dup = 0;
    V = gel(principal,i);
    j = i2-1;
    while (j > 0 && lg(gel(principal,j)) == lg(V))
    {
      if (field_is_contained(gel(principal,j),V,0)) { dup=1; break; }
      j--;
    }
    if (!dup) gel(principal,i2++) = V;
  }
  setlg(principal, i2);

  /* a subfield is generating iff all overfields contain the first overfield */
  ismax = cgetg(lg(principal),t_VECSMALL);
  isgene = cgetg(lg(principal),t_VECSMALL);
  for (i=1; i<lg(principal); i++)
  {
    V = gel(principal,i);
    ismax[i] = flmax;
    isgene[i] = flgene;
    W1 = NULL; /* intersection of strict overfields */
    for (j=i+1; j<lg(principal); j++)
    {
      W = gel(principal,j);
      if (!field_is_contained(V,W,1)) continue;
      ismax[i] = 0;
      if (!flgene) break;
      if (!W1) { W1 = W; continue; }
      if (!field_is_contained(W1,W,1))
      {
        W1 = intersect(W1,W);
        if (lg(W1)==lg(V)) { isgene[i]=0; break; }
      }
    }
  }

  for (i=1; i<lg(principal); i++)
  {
    nbmax += ismax[i];
    nbgene += isgene[i];
  }

  if (flmax)
  {
    Lmax = cgetg(nbmax+1, t_VEC);
    j=1;
    for (i=1; i<lg(principal); i++)
      if (ismax[i]) gel(Lmax,j++) = gel(principal,i);
  }

  if (flgene)
  {
    Lgene = cgetg(nbgene+1, t_VEC);
    j=1;
    for (i=1; i<lg(principal); i++)
      if (isgene[i]) gel(Lgene,j++) = gel(principal,i);
  }

  if (!flgene) res = Lmax;
  else if (!flmax) res = Lgene;
  else res = mkvec2(Lmax,Lgene);
  return gerepilecopy(av, res);
}

GEN
nfsubfieldsmax(GEN nf, long fl)
{
  pari_sp av = avma;
  GEN pol, fa, Lmax, V;
  long n, i, ero;

  subfields_cleanup(&nf, &pol, &n, &fa);
  if (n==1) { set_avma(av); return cgetg(1,t_VEC); }
  if (uisprime(n))
    return gerepilecopy(av, fl==1 ? mkvec(pol_x(varn(pol)))
      : mkvec(mkvec2(pol_x(varn(pol)),gen_0)));
  ero = subfields_get_ero(pol, nf);
  fa = subfields_get_fa(pol, nf, fa);
  Lmax = maxgen_subfields(pol, fa, subf_MAXIMAL);
  for (i=1; i<lg(Lmax); i++)
  {
    V = gel(Lmax,i);
    gel(Lmax,i) = subfield_generator(pol, V, lg(V)-1, ero, fl);
  }
  delete_var();
  return gerepilecopy(av, Lmax);
}

static void
heap_climb(GEN* H, long i)
{
  long j;
  if (i==1) return;
  j = i/2;
  if (cmp_universal(gel(*H,i),gel(*H,j)) > 0)
  {
    swap(gel(*H,i), gel(*H,j));
    return heap_climb(H,j);
  }
}

static void
heap_push(GEN* H, long *len, GEN x)
{
  if (*len+1 == lg(*H))
  {
    GEN H2 = zerovec(2*(*len));
    long i;
    for(i=1; i<lg(*H); i++)
      gel(H2,i) = gel(*H,i);
    *H = H2;
  }
  (*len)++;
  gel(*H,*len) = x;
  return heap_climb(H,*len);
}

static void
heap_descend(GEN* H, long len, long i)
{
  long maxi = i, j = 2*i;
  if (j > len) return;
  if (cmp_universal(gel(*H,j),gel(*H,i)) > 0) maxi = j;
  j++;
  if (j<=len && cmp_universal(gel(*H,j),gel(*H,maxi))>0) maxi = j;
  if (maxi == i) return;
  swap(gel(*H,i), gel(*H,maxi));
  return heap_descend(H,len,maxi);
}

static void
heap_pop(GEN *H, long *len, GEN* top)
{
  *top = gel(*H,1);
  gel(*H,1) = gel(*H,*len);
  (*len)--;
  return heap_descend(H,*len,1);
};

static GEN
nfsubfields_fa(GEN nf, long d, long fl)
{
  pari_sp av = avma;
  GEN pol, fa, gene, res, res2, H, V, v, W, w, data;
  long n, r, i, j, nres, len, s, newfield, ero, vp;

  subfields_cleanup(&nf, &pol, &n, &fa); vp = varn(pol);
  if (d && (d<1 || d>n || n%d)) return gerepilecopy(av, cgetg(1,t_VEC));
  if (!d && uisprime(n)) return gerepilecopy(av,
    fl==1 ? mkvec2( pol_x(varn(pol)), pol)
          : mkvec2( mkvec2(pol_x(vp),pol_0(vp)), mkvec2(pol,pol_x(vp))));
  if (n==1 || d==1) return gerepilecopy(av,
    fl==1 ? mkvec(pol_x(varn(pol))): _subfield(pol_x(vp),pol_0(vp)));
  if (d==n) return gerepilecopy(av,
    fl==1 ? mkvec(pol): _subfield(pol,pol_x(vp)));
  ero = subfields_get_ero(pol, nf);
  fa = subfields_get_fa(pol, nf, fa);
  gene = maxgen_subfields(pol, fa, subf_GENERATING);

  if (d)
  {
    /* keep only generating subfields of degree a multiple of d */
    j=1;
    for (i=1; i<lg(gene); i++)
      if ((lg(gel(gene,i))-1) % d == 0)
      {
        gel(gene,j) = gel(gene,i);
        j++;
      }
    setlg(gene,j);
  }
  r = lg(gene)-1;

  res = zerovec(10);
  nres = 0;
  H = zerovec(10);
  gel(H,1) = mkvec3(matid(n),zero_F2v(r),mkvecsmall(0));
  len = 1;

  while (len>0)
  {
    heap_pop(&H, &len, &data);
    V = gel(data,1);
    v = gel(data,2);
    s = gel(data,3)[1];
    for (i=s+1; i<=r; i++)
      if (!F2v_coeff(v,i))
      {
        W = vec_Q_primpart(intersect(V, gel(gene,i)));
        w = F2v_copy(v);
        F2v_set(w, i);
        newfield = 1;
        for (j=1; j<=r; j++)
          if (!F2v_coeff(w,j) && field_is_contained(W,gel(gene,j),1))
          {
            if (j<i) { newfield = 0; break; }
            F2v_set(w,j);
          }
        if (newfield && (!d || (lg(W)-1)%d==0)) heap_push(&H, &len, mkvec3(W,w,mkvecsmall(i)));
      }

    if (!d || lg(V)-1==d)
    {
      nres++;
      if (nres == lg(res))
      {
        res2 = zerovec(2*lg(res));
        for(j=1; j<lg(res); j++) gel(res2,j) = gel(res,j);
        res = res2;
      }
      gel(res,nres) = subfield_generator(pol, V, lg(V)-1, ero, fl);
    }
  }
  setlg(res,nres+1);
  vecreverse_inplace(res);

  delete_var();
  return gerepilecopy(av, res);
}
