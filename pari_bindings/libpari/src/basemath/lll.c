/* Copyright (C) 2008  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_qflll

static double
pari_rint(double a)
{
#ifdef HAS_RINT
  return rint(a);
#else
  const double two_to_52 = 4.5035996273704960e+15;
  double fa = fabs(a);
  double r = two_to_52 + fa;
  if (fa >= two_to_52) {
    r = a;
  } else {
    r = r - two_to_52;
    if (a < 0) r = -r;
  }
  return r;
#endif
}

/* default quality ratio for LLL */
static const double LLLDFT = 0.99;

/* assume flag & (LLL_KER|LLL_IM|LLL_ALL). LLL_INPLACE implies LLL_IM */
static GEN
lll_trivial(GEN x, long flag)
{
  if (lg(x) == 1)
  { /* dim x = 0 */
    if (! (flag & LLL_ALL)) return cgetg(1,t_MAT);
    retmkvec2(cgetg(1,t_MAT), cgetg(1,t_MAT));
  }
  /* dim x = 1 */
  if (gequal0(gel(x,1)))
  {
    if (flag & LLL_KER) return matid(1);
    if (flag & (LLL_IM|LLL_INPLACE)) return cgetg(1,t_MAT);
    retmkvec2(matid(1), cgetg(1,t_MAT));
  }
  if (flag & LLL_INPLACE) return gcopy(x);
  if (flag & LLL_KER) return cgetg(1,t_MAT);
  if (flag & LLL_IM)  return matid(1);
  retmkvec2(cgetg(1,t_MAT), (flag & LLL_GRAM)? gcopy(x): matid(1));
}

/* vecslice(x,#x-k,#x) in place. Works for t_MAT, t_VEC/t_COL */
static GEN
vectail_inplace(GEN x, long k)
{
  if (!k) return x;
  x[k] = ((ulong)x[0] & ~LGBITS) | evallg(lg(x) - k);
  return x + k;
}

/* k = dim Kernel */
static GEN
lll_finish(GEN h, long k, long flag)
{
  GEN g;
  if (!(flag & (LLL_IM|LLL_KER|LLL_ALL|LLL_INPLACE))) return h;
  if (flag & (LLL_IM|LLL_INPLACE)) return vectail_inplace(h, k);
  if (flag & LLL_KER) { setlg(h,k+1); return h; }
  g = vecslice(h,1,k); /* done first: vectail_inplace kills h */
  return mkvec2(g, vectail_inplace(h, k));
}

/* y * z * 2^e, e >= 0; y,z t_INT */
INLINE GEN
mulshift(GEN y, GEN z, long e)
{
  long ly = lgefint(y), lz;
  pari_sp av;
  GEN t;
  if (ly == 2) return gen_0;
  lz = lgefint(z);
  av = avma; (void)new_chunk(ly+lz+nbits2lg(e)); /* HACK */
  t = mulii(z, y);
  set_avma(av); return shifti(t, e);
}

/* x - y * z * 2^e, e >= 0; x,y,z t_INT */
INLINE GEN
submulshift(GEN x, GEN y, GEN z, long e)
{
  long lx = lgefint(x), ly, lz;
  pari_sp av;
  GEN t;
  if (!e) return submulii(x, y, z);
  if (lx == 2) { t = mulshift(y, z, e); togglesign(t); return t; }
  ly = lgefint(y);
  if (ly == 2) return icopy(x);
  lz = lgefint(z);
  av = avma; (void)new_chunk(lx+ly+lz+nbits2lg(e)); /* HACK */
  t = shifti(mulii(z, y), e);
  set_avma(av); return subii(x, t);
}
static void
subzi(GEN *a, GEN b)
{
  pari_sp av = avma;
  b = subii(*a, b);
  if (lgefint(b)<=lg(*a) && isonstack(*a)) { affii(b,*a); set_avma(av); }
  else *a = b;
}

static void
addzi(GEN *a, GEN b)
{
  pari_sp av = avma;
  b = addii(*a, b);
  if (lgefint(b)<=lg(*a) && isonstack(*a)) { affii(b,*a); set_avma(av); }
  else *a = b;
}

/* x - u*y * 2^e */
INLINE GEN
submuliu2n(GEN x, GEN y, ulong u, long e)
{
  pari_sp av;
  long ly = lgefint(y);
  if (ly == 2) return x;
  av = avma;
  (void)new_chunk(3+ly+lgefint(x)+nbits2lg(e)); /* HACK */
  y = shifti(mului(u,y), e);
  set_avma(av); return subii(x, y);
}
/* *x -= u*y * 2^e */
INLINE void
submulzu2n(GEN *x, GEN y, ulong u, long e)
{
  pari_sp av;
  long ly = lgefint(y);
  if (ly == 2) return;
  av = avma;
  (void)new_chunk(3+ly+lgefint(*x)+nbits2lg(e)); /* HACK */
  y = shifti(mului(u,y), e);
  set_avma(av); return subzi(x, y);
}

/* x + u*y * 2^e */
INLINE GEN
addmuliu2n(GEN x, GEN y, ulong u, long e)
{
  pari_sp av;
  long ly = lgefint(y);
  if (ly == 2) return x;
  av = avma;
  (void)new_chunk(3+ly+lgefint(x)+nbits2lg(e)); /* HACK */
  y = shifti(mului(u,y), e);
  set_avma(av); return addii(x, y);
}

/* *x += u*y * 2^e */
INLINE void
addmulzu2n(GEN *x, GEN y, ulong u, long e)
{
  pari_sp av;
  long ly = lgefint(y);
  if (ly == 2) return;
  av = avma;
  (void)new_chunk(3+ly+lgefint(*x)+nbits2lg(e)); /* HACK */
  y = shifti(mului(u,y), e);
  set_avma(av); return addzi(x, y);
}

/* n < 10; gerepileall supporting &NULL arguments. Maybe rename and export ? */
INLINE void
gc_lll(pari_sp av, int n, ...)
{
  int i, j;
  GEN *gptr[10];
  size_t s;
  va_list a; va_start(a, n);
  for (i=j=0; i<n; i++)
  {
    GEN *x = va_arg(a,GEN*);
    if (*x) { gptr[j++] = x; *x = (GEN)copy_bin(*x); }
  }
  va_end(a); set_avma(av);
  for (--j; j>=0; j--) *gptr[j] = bin_copy((GENbin*)*gptr[j]);
  s = pari_mainstack->top - pari_mainstack->bot;
  /* size of saved objects ~ stacksize / 4 => overflow */
  if (av - avma > (s >> 2))
  {
    size_t t = avma - pari_mainstack->bot;
    av = avma; new_chunk((s + t) / sizeof(long)); set_avma(av); /* double */
  }
}

/********************************************************************/
/**                                                                **/
/**                   FPLLL (adapted from D. Stehle's code)        **/
/**                                                                **/
/********************************************************************/
/* Babai* and fplll* are a conversion to libpari API and data types
   of fplll-1.3 by Damien Stehle'.

  Copyright 2005, 2006 Damien Stehle'.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program implements ideas from the paper "Floating-point LLL Revisited",
  by Phong Nguyen and Damien Stehle', in the Proceedings of Eurocrypt'2005,
  Springer-Verlag; and was partly inspired by Shoup's NTL library:
  http://www.shoup.net/ntl/ */

/* x t_REAL, |x| >= 1/2. Test whether |x| <= 3/2 */
static int
absrsmall2(GEN x)
{
  long e = expo(x), l, i;
  if (e < 0) return 1;
  if (e > 0 || (ulong)x[2] > (3UL << (BITS_IN_LONG-2))) return 0;
  /* line above assumes l > 2. OK since x != 0 */
  l = lg(x); for (i = 3; i < l; i++) if (x[i]) return 0;
  return 1;
}
/* x t_REAL; test whether |x| <= 1/2 */
static int
absrsmall(GEN x)
{
  long e, l, i;
  if (!signe(x)) return 1;
  e = expo(x); if (e < -1) return 1;
  if (e > -1 || (ulong)x[2] > HIGHBIT) return 0;
  l = lg(x); for (i = 3; i < l; i++) if (x[i]) return 0;
  return 1;
}

static void
rotate(GEN A, long k2, long k)
{
  long i;
  GEN B = gel(A,k2);
  for (i = k2; i > k; i--) gel(A,i) = gel(A,i-1);
  gel(A,k) = B;
}

/************************* FAST version (double) ************************/
#define dmael(x,i,j) ((x)[i][j])
#define del(x,i) ((x)[i])

static double *
cget_dblvec(long d)
{ return (double*) stack_malloc_align(d*sizeof(double), sizeof(double)); }

static double **
cget_dblmat(long d) { return (double **) cgetg(d, t_VECSMALL); }

static double
itodbl_exp(GEN x, long *e)
{
  pari_sp av = avma;
  GEN r = itor(x,DEFAULTPREC);
  *e = expo(r); setexpo(r,0);
  return gc_double(av, rtodbl(r));
}

static double
dbldotproduct(double *x, double *y, long n)
{
  long i;
  double sum = del(x,1) * del(y,1);
  for (i=2; i<=n; i++) sum += del(x,i) * del(y,i);
  return sum;
}

static double
dbldotsquare(double *x, long n)
{
  long i;
  double sum = del(x,1) * del(x,1);
  for (i=2; i<=n; i++) sum += del(x,i) * del(x,i);
  return sum;
}

static long
set_line(double *appv, GEN v, long n)
{
  long i, maxexp = 0;
  pari_sp av = avma;
  GEN e = cgetg(n+1, t_VECSMALL);
  for (i = 1; i <= n; i++)
  {
    del(appv,i) = itodbl_exp(gel(v,i), e+i);
    if (e[i] > maxexp) maxexp = e[i];
  }
  for (i = 1; i <= n; i++) del(appv,i) = ldexp(del(appv,i), e[i]-maxexp);
  set_avma(av); return maxexp;
}

static void
dblrotate(double **A, long k2, long k)
{
  long i;
  double *B = del(A,k2);
  for (i = k2; i > k; i--) del(A,i) = del(A,i-1);
  del(A,k) = B;
}
/* update G[kappa][i] from appB */
static void
setG_fast(double **appB, long n, double **G, long kappa, long a, long b)
{ long i;
  for (i = a; i <= b; i++)
    dmael(G,kappa,i) = dbldotproduct(del(appB,kappa), del(appB,i), n);
}
/* update G[i][kappa] from appB */
static void
setG2_fast(double **appB, long n, double **G, long kappa, long a, long b)
{ long i;
  for (i = a; i <= b; i++)
    dmael(G,i,kappa) = dbldotproduct(del(appB,kappa), del(appB,i), n);
}
const long EX0 = -2; /* uninitialized; any value less than expo(0.51) = -1 */

#ifdef LONG_IS_64BIT
typedef long s64;
#define addmuliu64_inplace addmuliu_inplace
#define submuliu64_inplace submuliu_inplace
#define submuliu642n submuliu2n
#define addmuliu642n addmuliu2n
#else
typedef long long s64;
typedef unsigned long long u64;

INLINE GEN
u64toi(u64 x)
{
  GEN y;
  ulong h;
  if (!x) return gen_0;
  h = x>>32;
  if (!h) return utoipos(x);
  y = cgetipos(4);
  *int_LSW(y) = x&0xFFFFFFFF;
  *int_MSW(y) = x>>32;
  return y;
}

INLINE GEN
u64toineg(u64 x)
{
  GEN y;
  ulong h;
  if (!x) return gen_0;
  h = x>>32;
  if (!h) return utoineg(x);
  y = cgetineg(4);
  *int_LSW(y) = x&0xFFFFFFFF;
  *int_MSW(y) = x>>32;
  return y;
}
INLINE GEN
addmuliu64_inplace(GEN x, GEN y, u64 u) { return addmulii(x, y, u64toi(u)); }

INLINE GEN
submuliu64_inplace(GEN x, GEN y, u64 u) { return submulii(x, y, u64toi(u)); }

INLINE GEN
addmuliu642n(GEN x, GEN y, u64 u, long e) { return submulshift(x, y, u64toineg(u), e); }

INLINE GEN
submuliu642n(GEN x, GEN y, u64 u, long e) { return submulshift(x, y, u64toi(u), e); }

#endif

/* Babai's Nearest Plane algorithm (iterative); see Babai() */
static int
Babai_fast(pari_sp av, long kappa, GEN *pB, GEN *pU, double **mu, double **r,
           double *s, double **appB, GEN expoB, double **G,
           long a, long zeros, long maxG, double eta)
{
  GEN B = *pB, U = *pU;
  const long n = nbrows(B), d = U ? lg(U)-1: 0;
  long k, aa = (a > zeros)? a : zeros+1;
  long emaxmu = EX0, emax2mu = EX0;
  s64 xx;
  int did_something = 0;
  /* N.B: we set d = 0 (resp. n = 0) to avoid updating U (resp. B) */

  for (;;) {
    int go_on = 0;
    long i, j, emax3mu = emax2mu;

    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Babai[1], a=%ld", aa);
      gc_lll(av,2,&B,&U);
    }
    /* Step2: compute the GSO for stage kappa */
    emax2mu = emaxmu; emaxmu = EX0;
    for (j=aa; j<kappa; j++)
    {
      double g = dmael(G,kappa,j);
      for (k = zeros+1; k < j; k++) g -= dmael(mu,j,k) * dmael(r,kappa,k);
      dmael(r,kappa,j) = g;
      dmael(mu,kappa,j) = dmael(r,kappa,j) / dmael(r,j,j);
      emaxmu = maxss(emaxmu, expoB[kappa]-expoB[j]);
    }
    /* maxmu doesn't decrease fast enough */
    if (emax3mu != EX0 && emax3mu <= emax2mu + 5) {*pB = B; *pU = U; return 1;}

    for (j=kappa-1; j>zeros; j--)
    {
      double tmp = fabs(ldexp (dmael(mu,kappa,j), expoB[kappa]-expoB[j]));
      if (tmp>eta) { go_on = 1; break; }
    }

    /* Step3--5: compute the X_j's  */
    if (go_on)
      for (j=kappa-1; j>zeros; j--)
      { /* The code below seemingly handles U = NULL, but in this case d = 0 */
        int e = expoB[j] - expoB[kappa];
        double tmp = ldexp(dmael(mu,kappa,j), -e), atmp = fabs(tmp);
        /* tmp = Inf is allowed */
        if (atmp <= .5) continue; /* size-reduced */
        if (gc_needed(av,2))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"Babai[2], a=%ld, j=%ld", aa,j);
          gc_lll(av,2,&B,&U);
        }
        did_something = 1;
        /* we consider separately the case |X| = 1 */
        if (atmp <= 1.5)
        {
          if (dmael(mu,kappa,j) > 0) { /* in this case, X = 1 */
            for (k=zeros+1; k<j; k++)
              dmael(mu,kappa,k) -= ldexp(dmael(mu,j,k), e);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = subii(gmael(B,kappa,i), gmael(B,j,i));
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = subii(gmael(U,kappa,i), gmael(U,j,i));
          } else { /* otherwise X = -1 */
            for (k=zeros+1; k<j; k++)
              dmael(mu,kappa,k) += ldexp(dmael(mu,j,k), e);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = addii(gmael(B,kappa,i), gmael(B,j,i));
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = addii(gmael(U,kappa,i), gmael(U,j,i));
          }
          continue;
        }
        /* we have |X| >= 2 */
        if (atmp < 9007199254740992.)
        {
          tmp = pari_rint(tmp);
          for (k=zeros+1; k<j; k++)
            dmael(mu,kappa,k) -= ldexp(tmp * dmael(mu,j,k), e);
          xx = (s64) tmp;
          if (xx > 0) /* = xx */
          {
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = submuliu64_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = submuliu64_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
          }
          else /* = -xx */
          {
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = addmuliu64_inplace(gmael(B,kappa,i), gmael(B,j,i), -xx);
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = addmuliu64_inplace(gmael(U,kappa,i), gmael(U,j,i), -xx);
          }
        }
        else
        {
          int E;
          xx = (s64) ldexp(frexp(dmael(mu,kappa,j), &E), 53);
          E -= e + 53;
          if (E <= 0)
          {
            xx = xx << -E;
            for (k=zeros+1; k<j; k++)
              dmael(mu,kappa,k) -= ldexp(((double)xx) * dmael(mu,j,k), e);
            if (xx > 0) /* = xx */
            {
              for (i=1; i<=n; i++)
                gmael(B,kappa,i) = submuliu64_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
              for (i=1; i<=d; i++)
                gmael(U,kappa,i) = submuliu64_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
            }
            else /* = -xx */
            {
              for (i=1; i<=n; i++)
                gmael(B,kappa,i) = addmuliu64_inplace(gmael(B,kappa,i), gmael(B,j,i), -xx);
              for (i=1; i<=d; i++)
                gmael(U,kappa,i) = addmuliu64_inplace(gmael(U,kappa,i), gmael(U,j,i), -xx);
            }
          } else
          {
            for (k=zeros+1; k<j; k++)
              dmael(mu,kappa,k) -= ldexp(((double)xx) * dmael(mu,j,k), E + e);
            if (xx > 0) /* = xx */
            {
              for (i=1; i<=n; i++)
                gmael(B,kappa,i) = submuliu642n(gmael(B,kappa,i), gmael(B,j,i), xx, E);
              for (i=1; i<=d; i++)
                gmael(U,kappa,i) = submuliu642n(gmael(U,kappa,i), gmael(U,j,i), xx, E);
            }
            else /* = -xx */
            {
              for (i=1; i<=n; i++)
                gmael(B,kappa,i) = addmuliu642n(gmael(B,kappa,i), gmael(B,j,i), -xx, E);
              for (i=1; i<=d; i++)
                gmael(U,kappa,i) = addmuliu642n(gmael(U,kappa,i), gmael(U,j,i), -xx, E);
            }
          }
        }
      }
    if (!go_on) break; /* Anything happened? */
    expoB[kappa] = set_line(del(appB,kappa), gel(B,kappa), n);
    setG_fast(appB, n, G, kappa, zeros+1, kappa-1);
    aa = zeros+1;
  }
  if (did_something) setG2_fast(appB, n, G, kappa, kappa, maxG);

  del(s,zeros+1) = dmael(G,kappa,kappa);
  /* the last s[kappa-1]=r[kappa][kappa] is computed only if kappa increases */
  for (k=zeros+1; k<=kappa-2; k++)
    del(s,k+1) = del(s,k) - dmael(mu,kappa,k)*dmael(r,kappa,k);
  *pB = B; *pU = U; return 0;
}

static void
update_alpha(GEN alpha, long kappa, long kappa2, long kappamax)
{
  long i;
  for (i = kappa; i < kappa2; i++)
    if (kappa <= alpha[i]) alpha[i] = kappa;
  for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];
  for (i = kappa2+1; i <= kappamax; i++)
    if (kappa < alpha[i]) alpha[i] = kappa;
  alpha[kappa] = kappa;
}
static void
rotateG(GEN G, long kappa2, long kappa, long maxG, GEN Gtmp)
{
  long i, j;
  for (i=1; i<=kappa2; i++) gel(Gtmp,i) = gmael(G,kappa2,i);
  for (   ; i<=maxG; i++)   gel(Gtmp,i) = gmael(G,i,kappa2);
  for (i=kappa2; i>kappa; i--)
    {
      for (j=1; j<kappa; j++) gmael(G,i,j) = gmael(G,i-1,j);
      gmael(G,i,kappa) = gel(Gtmp,i-1);
      for (j=kappa+1; j<=i; j++) gmael(G,i,j) = gmael(G,i-1,j-1);
      for (j=kappa2+1; j<=maxG; j++) gmael(G,j,i) = gmael(G,j,i-1);
    }
  for (i=1; i<kappa; i++) gmael(G,kappa,i) = gel(Gtmp,i);
  gmael(G,kappa,kappa) = gel(Gtmp,kappa2);
  for (i=kappa2+1; i<=maxG; i++) gmael(G,i,kappa) = gel(Gtmp,i);
}
static void
rotateG_fast(double **G, long kappa2, long kappa, long maxG, double *Gtmp)
{
  long i, j;
  for (i=1; i<=kappa2; i++) del(Gtmp,i) = dmael(G,kappa2,i);
  for (   ; i<=maxG; i++) del(Gtmp,i) = dmael(G,i,kappa2);
  for (i=kappa2; i>kappa; i--)
  {
    for (j=1; j<kappa; j++) dmael(G,i,j) = dmael(G,i-1,j);
    dmael(G,i,kappa) = del(Gtmp,i-1);
    for (j=kappa+1; j<=i; j++) dmael(G,i,j) = dmael(G,i-1,j-1);
    for (j=kappa2+1; j<=maxG; j++) dmael(G,j,i) = dmael(G,j,i-1);
  }
  for (i=1; i<kappa; i++) dmael(G,kappa,i) = del(Gtmp,i);
  dmael(G,kappa,kappa) = del(Gtmp,kappa2);
  for (i=kappa2+1; i<=maxG; i++) dmael(G,i,kappa) = del(Gtmp,i);
}

/* LLL-reduces (B,U) in place [apply base change transforms to B and U].
 * Gram matrix, and GSO performed on matrices of 'double'.
 * If (keepfirst), never swap with first vector.
 * Return -1 on failure, else zeros = dim Kernel (>= 0) */
static long
fplll_fast(GEN *pB, GEN *pU, double delta, double eta, long keepfirst)
{
  pari_sp av;
  long kappa, kappa2, d, n, i, j, zeros, kappamax, maxG;
  double **mu, **r, *s, tmp, *Gtmp, **G, **appB;
  GEN alpha, expoB, B = *pB, U;
  long cnt = 0;

  d = lg(B)-1;
  n = nbrows(B);
  U = *pU; /* NULL if inplace */

  G = cget_dblmat(d+1);
  appB = cget_dblmat(d+1);
  mu = cget_dblmat(d+1);
  r  = cget_dblmat(d+1);
  s  = cget_dblvec(d+1);
  for (j = 1; j <= d; j++)
  {
    del(mu,j) = cget_dblvec(d+1);
    del(r,j) = cget_dblvec(d+1);
    del(appB,j) = cget_dblvec(n+1);
    del(G,j) = cget_dblvec(d+1);
    for (i=1; i<=d; i++) dmael(G,j,i) = 0.;
  }
  expoB = cgetg(d+1, t_VECSMALL);
  for (i=1; i<=d; i++) expoB[i] = set_line(del(appB,i), gel(B,i), n);
  Gtmp = cget_dblvec(d+1);
  alpha = cgetg(d+1, t_VECSMALL);
  av = avma;

  /* Step2: Initializing the main loop */
  kappamax = 1;
  i = 1;
  maxG = d; /* later updated to kappamax */

  do {
    dmael(G,i,i) = dbldotsquare(del(appB,i),n);
  } while (dmael(G,i,i) <= 0 && (++i <=d));
  zeros = i-1; /* all vectors B[i] with i <= zeros are zero vectors */
  kappa = i;
  if (zeros < d) dmael(r,zeros+1,zeros+1) = dmael(G,zeros+1,zeros+1);
  for (i=zeros+1; i<=d; i++) alpha[i]=1;
  while (++kappa <= d)
  {
    if (kappa > kappamax)
    {
      if (DEBUGLEVEL>=4) err_printf("K%ld ",kappa);
      maxG = kappamax = kappa;
      setG_fast(appB, n, G, kappa, zeros+1, kappa);
    }
    /* Step3: Call to the Babai algorithm, mu,r,s updated in place */
    if (Babai_fast(av, kappa, &B,&U, mu,r,s, appB, expoB, G, alpha[kappa],
                   zeros, maxG, eta)) { *pB=B; *pU=U; return -1; }

    tmp = ldexp(r[kappa-1][kappa-1] * delta, 2*(expoB[kappa-1]-expoB[kappa]));
    if ((keepfirst && kappa == 2) || tmp <= del(s,kappa-1))
    { /* Step4: Success of Lovasz's condition */
      alpha[kappa] = kappa;
      tmp = dmael(mu,kappa,kappa-1) * dmael(r,kappa,kappa-1);
      dmael(r,kappa,kappa) = del(s,kappa-1)- tmp;
      continue;
    }
    /* Step5: Find the right insertion index kappa, kappa2 = initial kappa */
    if (DEBUGLEVEL>=4 && kappa==kappamax && del(s,kappa-1)!=0)
      if (++cnt > 20) { cnt = 0; err_printf("(%ld) ", 2*expoB[1] + dblexpo(del(s,1))); }
    kappa2 = kappa;
    do {
      kappa--;
      if (kappa<zeros+2 + (keepfirst ? 1: 0)) break;
      tmp = dmael(r,kappa-1,kappa-1) * delta;
      tmp = ldexp(tmp, 2*(expoB[kappa-1]-expoB[kappa2]));
    } while (del(s,kappa-1) <= tmp);
    update_alpha(alpha, kappa, kappa2, kappamax);

    /* Step6: Update the mu's and r's */
    dblrotate(mu,kappa2,kappa);
    dblrotate(r,kappa2,kappa);
    dmael(r,kappa,kappa) = del(s,kappa);

    /* Step7: Update B, appB, U, G */
    rotate(B,kappa2,kappa);
    dblrotate(appB,kappa2,kappa);
    if (U) rotate(U,kappa2,kappa);
    rotate(expoB,kappa2,kappa);
    rotateG_fast(G,kappa2,kappa, maxG, Gtmp);

    /* Step8: Prepare the next loop iteration */
    if (kappa == zeros+1 && dmael(G,kappa,kappa)<= 0)
    {
      zeros++; kappa++;
      dmael(G,kappa,kappa) = dbldotsquare(del(appB,kappa),n);
      dmael(r,kappa,kappa) = dmael(G,kappa,kappa);
    }
  }
  *pB = B; *pU = U; return zeros;
}

/***************** HEURISTIC version (reduced precision) ****************/
static GEN
realsqrdotproduct(GEN x)
{
  long i, l = lg(x);
  GEN z = sqrr(gel(x,1));
  for (i=2; i<l; i++) z = addrr(z, sqrr(gel(x,i)));
  return z;
}
/* x, y non-empty vector of t_REALs, same length */
static GEN
realdotproduct(GEN x, GEN y)
{
  long i, l;
  GEN z;
  if (x == y) return realsqrdotproduct(x);
  l = lg(x); z = mulrr(gel(x,1),gel(y,1));
  for (i=2; i<l; i++) z = addrr(z, mulrr(gel(x,i), gel(y,i)));
  return z;
}
static void
setG_heuristic(GEN appB, GEN G, long kappa, long a, long b)
{ pari_sp av = avma;
  long i;
  for (i = a; i <= b; i++)
    affrr(realdotproduct(gel(appB,kappa),gel(appB,i)), gmael(G,kappa,i));
  set_avma(av);
}
static void
setG2_heuristic(GEN appB, GEN G, long kappa, long a, long b)
{ pari_sp av = avma;
  long i;
  for (i = a; i <= b; i++)
    affrr(realdotproduct(gel(appB,kappa),gel(appB,i)), gmael(G,i,kappa));
  set_avma(av);
}

/* approximate t_REAL x as m * 2^e, where |m| < 2^bit */
static GEN
truncexpo(GEN x, long bit, long *e)
{
  *e = expo(x) + 1 - bit;
  if (*e >= 0) return mantissa2nr(x, 0);
  *e = 0; return roundr_safe(x);
}
/* Babai's Nearest Plane algorithm (iterative); see Babai() */
static int
Babai_heuristic(pari_sp av, long kappa, GEN *pB, GEN *pU, GEN mu, GEN r, GEN s,
                GEN appB, GEN G, long a, long zeros, long maxG,
                GEN eta, long prec)
{
  GEN B = *pB, U = *pU;
  const long n = nbrows(B), d = U ? lg(U)-1: 0, bit = prec2nbits(prec);
  long k, aa = (a > zeros)? a : zeros+1;
  int did_something = 0;
  long emaxmu = EX0, emax2mu = EX0;
  /* N.B: we set d = 0 (resp. n = 0) to avoid updating U (resp. B) */

  for (;;) {
    int go_on = 0;
    long i, j, emax3mu = emax2mu;

    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Babai[1], a=%ld", aa);
      gc_lll(av,2,&B,&U);
    }
    /* Step2: compute the GSO for stage kappa */
    emax2mu = emaxmu; emaxmu = EX0;
    for (j=aa; j<kappa; j++)
    {
      pari_sp btop = avma;
      GEN g = gmael(G,kappa,j);
      for (k = zeros+1; k<j; k++)
        g = subrr(g, mulrr(gmael(mu,j,k), gmael(r,kappa,k)));
      affrr(g, gmael(r,kappa,j));
      affrr(divrr(gmael(r,kappa,j), gmael(r,j,j)), gmael(mu,kappa,j));
      emaxmu = maxss(emaxmu, expo(gmael(mu,kappa,j)));
      set_avma(btop);
    }
    if (emax3mu != EX0 && emax3mu <= emax2mu + 5)
    { *pB = B; *pU = U; return 1; }

    for (j=kappa-1; j>zeros; j--)
      if (abscmprr(gmael(mu,kappa,j), eta) > 0) { go_on = 1; break; }

    /* Step3--5: compute the X_j's  */
    if (go_on)
      for (j=kappa-1; j>zeros; j--)
      { /* The code below seemingly handles U = NULL, but in this case d = 0 */
        pari_sp btop;
        GEN tmp = gmael(mu,kappa,j);
        if (absrsmall(tmp)) continue; /* size-reduced */

        if (gc_needed(av,2))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"Babai[2], a=%ld, j=%ld", aa,j);
          gc_lll(av,2,&B,&U);
        }
        btop = avma; did_something = 1;
        /* we consider separately the case |X| = 1 */
        if (absrsmall2(tmp))
        {
          if (signe(tmp) > 0) { /* in this case, X = 1 */
            for (k=zeros+1; k<j; k++)
              affrr(subrr(gmael(mu,kappa,k), gmael(mu,j,k)), gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = subii(gmael(B,kappa,i), gmael(B,j,i));
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = subii(gmael(U,kappa,i), gmael(U,j,i));
          } else { /* otherwise X = -1 */
            for (k=zeros+1; k<j; k++)
              affrr(addrr(gmael(mu,kappa,k), gmael(mu,j,k)), gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = addii(gmael(B,kappa,i), gmael(B,j,i));
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = addii(gmael(U,kappa,i),gmael(U,j,i));
          }
          continue;
        }
        /* we have |X| >= 2 */
        if (expo(tmp) < BITS_IN_LONG)
        {
          ulong xx = roundr_safe(tmp)[2]; /* X fits in an ulong */
          if (signe(tmp) > 0) /* = xx */
          {
            for (k=zeros+1; k<j; k++)
              affrr(subrr(gmael(mu,kappa,k), mulur(xx, gmael(mu,j,k))),
                  gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = submuliu_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = submuliu_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
          }
          else /* = -xx */
          {
            for (k=zeros+1; k<j; k++)
              affrr(addrr(gmael(mu,kappa,k), mulur(xx, gmael(mu,j,k))),
                  gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = addmuliu_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = addmuliu_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
          }
        }
        else
        {
          long e;
          GEN X = truncexpo(tmp, bit, &e); /* tmp ~ X * 2^e */
          btop = avma;
          for (k=zeros+1; k<j; k++)
          {
            GEN x = mulir(X, gmael(mu,j,k));
            if (e) shiftr_inplace(x, e);
            affrr(subrr(gmael(mu,kappa,k), x), gmael(mu,kappa,k));
          }
          set_avma(btop);
          for (i=1; i<=n; i++)
            gmael(B,kappa,i) = submulshift(gmael(B,kappa,i), gmael(B,j,i), X, e);
          for (i=1; i<=d; i++)
            gmael(U,kappa,i) = submulshift(gmael(U,kappa,i), gmael(U,j,i), X, e);
        }
      }
    if (!go_on) break; /* Anything happened? */
    for (i=1 ; i<=n; i++) affir(gmael(B,kappa,i), gmael(appB,kappa,i));
    setG_heuristic(appB, G, kappa, zeros+1, kappa-1);
    aa = zeros+1;
  }
  if (did_something) setG2_heuristic(appB, G, kappa, kappa, maxG);
  affrr(gmael(G,kappa,kappa), gel(s,zeros+1));
  /* the last s[kappa-1]=r[kappa][kappa] is computed only if kappa increases */
  av = avma;
  for (k=zeros+1; k<=kappa-2; k++)
    affrr(subrr(gel(s,k), mulrr(gmael(mu,kappa,k), gmael(r,kappa,k))),
          gel(s,k+1));
  *pB = B; *pU = U; return gc_bool(av, 0);
}

static GEN
ZC_to_RC(GEN x, long prec)
{ pari_APPLY_type(t_COL,itor(gel(x,i),prec)) }

static GEN
ZM_to_RM(GEN x, long prec)
{ pari_APPLY_same(ZC_to_RC(gel(x,i),prec)) }

/* LLL-reduces (B,U) in place [apply base change transforms to B and U].
 * Gram matrix made of t_REAL at precision prec2, performe GSO at prec.
 * If (keepfirst), never swap with first vector.
 * Return -1 on failure, else zeros = dim Kernel (>= 0) */
static long
fplll_heuristic(GEN *pB, GEN *pU, double DELTA, double ETA, long keepfirst,
                long prec, long prec2)
{
  pari_sp av, av2;
  long kappa, kappa2, d, i, j, zeros, kappamax, maxG;
  GEN mu, r, s, tmp, Gtmp, alpha, G, appB, B = *pB, U;
  GEN delta = dbltor(DELTA), eta = dbltor(ETA);
  long cnt = 0;

  d = lg(B)-1;
  U = *pU; /* NULL if inplace */

  G = cgetg(d+1, t_MAT);
  mu = cgetg(d+1, t_MAT);
  r  = cgetg(d+1, t_MAT);
  s  = cgetg(d+1, t_VEC);
  appB = ZM_to_RM(B, prec2);
  for (j = 1; j <= d; j++)
  {
    GEN M = cgetg(d+1, t_COL), R = cgetg(d+1, t_COL), S = cgetg(d+1, t_COL);
    gel(mu,j)= M;
    gel(r,j) = R;
    gel(G,j) = S;
    gel(s,j) = cgetr(prec);
    for (i = 1; i <= d; i++)
    {
      gel(R,i) = cgetr(prec);
      gel(M,i) = cgetr(prec);
      gel(S,i) = cgetr(prec2);
    }
  }
  Gtmp = cgetg(d+1, t_VEC);
  alpha = cgetg(d+1, t_VECSMALL);
  av = avma;

  /* Step2: Initializing the main loop */
  kappamax = 1;
  i = 1;
  maxG = d; /* later updated to kappamax */

  do {
    affrr(RgV_dotsquare(gel(appB,i)), gmael(G,i,i));
  } while (signe(gmael(G,i,i)) == 0 && (++i <=d));
  zeros = i-1; /* all vectors B[i] with i <= zeros are zero vectors */
  kappa = i;
  if (zeros < d) affrr(gmael(G,zeros+1,zeros+1), gmael(r,zeros+1,zeros+1));
  for (i=zeros+1; i<=d; i++) alpha[i]=1;

  while (++kappa <= d)
  {
    if (kappa > kappamax)
    {
      if (DEBUGLEVEL>=4) err_printf("K%ld ",kappa);
      maxG = kappamax = kappa;
      setG_heuristic(appB, G, kappa, zeros+1, kappa);
    }
    /* Step3: Call to the Babai algorithm, mu,r,s updated in place */
    if (Babai_heuristic(av, kappa, &B,&U, mu,r,s, appB, G, alpha[kappa], zeros,
                        maxG, eta, prec)) { *pB = B; *pU = U; return -1; }
    av2 = avma;
    if ((keepfirst && kappa == 2) ||
        cmprr(mulrr(gmael(r,kappa-1,kappa-1), delta), gel(s,kappa-1)) <= 0)
    { /* Step4: Success of Lovasz's condition */
      alpha[kappa] = kappa;
      tmp = mulrr(gmael(mu,kappa,kappa-1), gmael(r,kappa,kappa-1));
      affrr(subrr(gel(s,kappa-1), tmp), gmael(r,kappa,kappa));
      set_avma(av2); continue;
    }
    /* Step5: Find the right insertion index kappa, kappa2 = initial kappa */
    if (DEBUGLEVEL>=4 && kappa==kappamax && signe(gel(s,kappa-1)))
      if (++cnt > 20) { cnt = 0; err_printf("(%ld) ", expo(gel(s,1))); }
    kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros+2 + (keepfirst ? 1: 0)) break;
      tmp = mulrr(gmael(r,kappa-1,kappa-1), delta);
    } while (cmprr(gel(s,kappa-1), tmp) <= 0 );
    set_avma(av2);
    update_alpha(alpha, kappa, kappa2, kappamax);

    /* Step6: Update the mu's and r's */
    rotate(mu,kappa2,kappa);
    rotate(r,kappa2,kappa);
    affrr(gel(s,kappa), gmael(r,kappa,kappa));

    /* Step7: Update B, appB, U, G */
    rotate(B,kappa2,kappa);
    rotate(appB,kappa2,kappa);
    if (U) rotate(U,kappa2,kappa);
    rotateG(G,kappa2,kappa, maxG, Gtmp);

    /* Step8: Prepare the next loop iteration */
    if (kappa == zeros+1 && !signe(gmael(G,kappa,kappa)))
    {
      zeros++; kappa++;
      affrr(RgV_dotsquare(gel(appB,kappa)), gmael(G,kappa,kappa));
      affrr(gmael(G,kappa,kappa), gmael(r,kappa,kappa));
    }
  }
  *pB=B; *pU=U; return zeros;
}

/************************* PROVED version (t_INT) ***********************/
/* dpe inspired by dpe.h by Patrick Pelissier, Paul Zimmermann
 * https://gforge.inria.fr/projects/dpe/
 */

typedef struct
{
  double d;  /* significand */
  long e; /* exponent */
} dpe_t;

#define Dmael(x,i,j) (&((x)[i][j]))
#define Del(x,i) (&((x)[i]))

static void
dperotate(dpe_t **A, long k2, long k)
{
  long i;
  dpe_t *B = A[k2];
  for (i = k2; i > k; i--) A[i] = A[i-1];
  A[k] = B;
}

static void
dpe_normalize0(dpe_t *x)
{
  int e;
  x->d = frexp(x->d, &e);
  x->e += e;
}

static void
dpe_normalize(dpe_t *x)
{
  if (x->d == 0.0)
    x->e = -LONG_MAX;
  else
    dpe_normalize0(x);
}

static GEN
dpetor(dpe_t *x)
{
  GEN r = dbltor(x->d);
  if (signe(r)==0) return r;
  setexpo(r, x->e-1);
  return r;
}

static void
affdpe(dpe_t *y, dpe_t *x)
{
  x->d = y->d;
  x->e = y->e;
}

static void
affidpe(GEN y, dpe_t *x)
{
  pari_sp av = avma;
  GEN r = itor(y, DEFAULTPREC);
  x->e = expo(r)+1;
  setexpo(r,-1);
  x->d = rtodbl(r);
  set_avma(av);
}

static void
affdbldpe(double y, dpe_t *x)
{
  x->d = (double)y;
  x->e = 0;
  dpe_normalize(x);
}

static void
dpe_mulz(dpe_t *x, dpe_t *y, dpe_t *z)
{
  z->d = x->d * y->d;
  if (z->d == 0.0)
    z->e = -LONG_MAX;
  else
  {
    z->e = x->e + y->e;
    dpe_normalize0(z);
  }
}

static void
dpe_divz(dpe_t *x, dpe_t *y, dpe_t *z)
{
  z->d = x->d / y->d;
  if (z->d == 0.0)
    z->e = -LONG_MAX;
  else
  {
    z->e = x->e - y->e;
    dpe_normalize0(z);
  }
}

static void
dpe_negz(dpe_t *y, dpe_t *x)
{
  x->d = - y->d;
  x->e = y->e;
}

static void
dpe_addz(dpe_t *y, dpe_t *z, dpe_t *x)
{
  if (y->e > z->e + 53)
    affdpe(y, x);
  else if (z->e > y->e + 53)
    affdpe(z, x);
  else
  {
    long d = y->e - z->e;

    if (d >= 0)
    {
      x->d = y->d + ldexp(z->d, -d);
      x->e  = y->e;
    }
    else
    {
      x->d = z->d + ldexp(y->d, d);
      x->e = z->e;
    }
    dpe_normalize(x);
  }
}
static void
dpe_subz(dpe_t *y, dpe_t *z, dpe_t *x)
{
  if (y->e > z->e + 53)
    affdpe(y, x);
  else if (z->e > y->e + 53)
    dpe_negz(z, x);
  else
  {
    long d = y->e - z->e;

    if (d >= 0)
    {
      x->d = y->d - ldexp(z->d, -d);
      x->e = y->e;
    }
    else
    {
      x->d = ldexp(y->d, d) - z->d;
      x->e = z->e;
    }
    dpe_normalize(x);
  }
}

static void
dpe_muluz(dpe_t *y, ulong t, dpe_t *x)
{
  x->d = y->d * (double)t;
  x->e = y->e;
  dpe_normalize(x);
}

static void
dpe_addmuluz(dpe_t *y,  dpe_t *z, ulong t, dpe_t *x)
{
  dpe_t tmp;
  dpe_muluz(z, t, &tmp);
  dpe_addz(y, &tmp, x);
}

static void
dpe_submuluz(dpe_t *y,  dpe_t *z, ulong t, dpe_t *x)
{
  dpe_t tmp;
  dpe_muluz(z, t, &tmp);
  dpe_subz(y, &tmp, x);
}

static void
dpe_submulz(dpe_t *y,  dpe_t *z, dpe_t *t, dpe_t *x)
{
  dpe_t tmp;
  dpe_mulz(z, t, &tmp);
  dpe_subz(y, &tmp, x);
}

static int
dpe_cmp(dpe_t *x, dpe_t *y)
{
  int sx = x->d < 0. ? -1: x->d > 0.;
  int sy = y->d < 0. ? -1: y->d > 0.;
  int d  = sx - sy;

  if (d != 0)
    return d;
  else if (x->e > y->e)
    return (sx > 0) ? 1 : -1;
  else if (y->e > x->e)
    return (sx > 0) ? -1 : 1;
  else
    return (x->d < y->d) ? -1 : (x->d > y->d);
}

static int
dpe_abscmp(dpe_t *x, dpe_t *y)
{
  if (x->e > y->e)
    return 1;
  else if (y->e > x->e)
    return -1;
  else
    return (fabs(x->d) < fabs(y->d)) ? -1 : (fabs(x->d) > fabs(y->d));
}

static int
dpe_abssmall(dpe_t *x)
{
  return (x->e <= 0) || (x->e == 1 && fabs(x->d) <= .75);
}

static int
dpe_cmpmul(dpe_t *x, dpe_t *y, dpe_t *z)
{
  dpe_t t;
  dpe_mulz(x,y,&t);
  return dpe_cmp(&t, z);
}

static dpe_t *
cget_dpevec(long d)
{ return (dpe_t*) stack_malloc_align(d*sizeof(dpe_t), sizeof(dpe_t)); }

static dpe_t **
cget_dpemat(long d) { return (dpe_t **) cgetg(d, t_VECSMALL); }

static GEN
dpeM_diagonal_shallow(dpe_t **m, long d)
{
  long i;
  GEN y = cgetg(d+1,t_VEC);
  for (i=1; i<=d; i++) gel(y, i) = dpetor(Dmael(m,i,i));
  return y;
}

static void
affii_or_copy_gc(pari_sp av, GEN x, GEN *y)
{
  long l = lg(*y);
  if (lgefint(x) <= l && isonstack(*y))
  {
    affii(x,*y);
    set_avma(av);
  }
  else
    *y = gerepileuptoint(av, x);
}

/* *x -= u*y */
INLINE void
submulziu(GEN *x, GEN y, ulong u)
{
  pari_sp av;
  long ly = lgefint(y);
  if (ly == 2) return;
  av = avma;
  (void)new_chunk(3+ly+lgefint(*x)); /* HACK */
  y = mului(u,y);
  set_avma(av); subzi(x, y);
}

/* *x += u*y */
INLINE void
addmulziu(GEN *x, GEN y, ulong u)
{
  pari_sp av;
  long ly = lgefint(y);
  if (ly == 2) return;
  av = avma;
  (void)new_chunk(3+ly+lgefint(*x)); /* HACK */
  y = mului(u,y);
  set_avma(av); addzi(x, y);
}

/************************** PROVED version (dpe) *************************/

/* Babai's Nearest Plane algorithm (iterative).
 * Size-reduces b_kappa using mu_{i,j} and r_{i,j} for j<=i <kappa
 * Update B[,kappa]; compute mu_{kappa,j}, r_{kappa,j} for j<=kappa and s[kappa]
 * mu, r, s updated in place (affrr). Return 1 on failure, else 0. */
static int
Babai_dpe(pari_sp av, long kappa, GEN *pG, GEN *pB, GEN *pU, dpe_t **mu, dpe_t **r, dpe_t *s,
      long a, long zeros, long maxG, dpe_t *eta)
{
  GEN G = *pG, B = *pB, U = *pU, ztmp;
  long k, d, n, aa = a > zeros? a: zeros+1;
  long emaxmu = EX0, emax2mu = EX0;
  /* N.B: we set d = 0 (resp. n = 0) to avoid updating U (resp. B) */
  d = U? lg(U)-1: 0;
  n = B? nbrows(B): 0;
  for (;;) {
    int go_on = 0;
    long i, j, emax3mu = emax2mu;

    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Babai[1], a=%ld", aa);
      gc_lll(av,3,&G,&B,&U);
    }
    /* Step2: compute the GSO for stage kappa */
    emax2mu = emaxmu; emaxmu = EX0;
    for (j=aa; j<kappa; j++)
    {
      dpe_t g;
      affidpe(gmael(G,kappa,j), &g);
      for (k = zeros+1; k < j; k++)
        dpe_submulz(&g, Dmael(mu,j,k), Dmael(r,kappa,k), &g);
      affdpe(&g, Dmael(r,kappa,j));
      dpe_divz(Dmael(r,kappa,j), Dmael(r,j,j), Dmael(mu,kappa,j));
      emaxmu = maxss(emaxmu, Dmael(mu,kappa,j)->e);
    }
    if (emax3mu != EX0 && emax3mu <= emax2mu + 5) /* precision too low */
    { *pG = G; *pB = B; *pU = U; return 1; }

    for (j=kappa-1; j>zeros; j--)
      if (dpe_abscmp(Dmael(mu,kappa,j), eta) > 0) { go_on = 1; break; }

    /* Step3--5: compute the X_j's  */
    if (go_on)
      for (j=kappa-1; j>zeros; j--)
      {
        pari_sp btop;
        dpe_t *tmp = Dmael(mu,kappa,j);
        if (tmp->e < 0) continue; /* (essentially) size-reduced */

        if (gc_needed(av,2))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"Babai[2], a=%ld, j=%ld", aa,j);
          gc_lll(av,3,&G,&B,&U);
        }
        /* we consider separately the case |X| = 1 */
        if (dpe_abssmall(tmp))
        {
          if (tmp->d > 0) { /* in this case, X = 1 */
            for (k=zeros+1; k<j; k++)
              dpe_subz(Dmael(mu,kappa,k), Dmael(mu,j,k), Dmael(mu,kappa,k));
            for (i=1; i<=n; i++)
              subzi(&gmael(B,kappa,i), gmael(B,j,i));
            for (i=1; i<=d; i++)
              subzi(&gmael(U,kappa,i), gmael(U,j,i));
            btop = avma;
            ztmp = subii(gmael(G,j,j), shifti(gmael(G,kappa,j), 1));
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            affii_or_copy_gc(btop, ztmp, &gmael(G,kappa,kappa));
            for (i=1; i<=j; i++)
              subzi(&gmael(G,kappa,i), gmael(G,j,i));
            for (i=j+1; i<kappa; i++)
              subzi(&gmael(G,kappa,i), gmael(G,i,j));
            for (i=kappa+1; i<=maxG; i++)
              subzi(&gmael(G,i,kappa), gmael(G,i,j));
          } else { /* otherwise X = -1 */
            for (k=zeros+1; k<j; k++)
              dpe_addz(Dmael(mu,kappa,k), Dmael(mu,j,k), Dmael(mu,kappa,k));
            for (i=1; i<=n; i++)
              addzi(&gmael(B,kappa,i),gmael(B,j,i));
            for (i=1; i<=d; i++)
              addzi(&gmael(U,kappa,i),gmael(U,j,i));
            btop = avma;
            ztmp = addii(gmael(G,j,j), shifti(gmael(G,kappa,j), 1));
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            affii_or_copy_gc(btop, ztmp, &gmael(G,kappa,kappa));
            for (i=1; i<=j; i++)
              addzi(&gmael(G,kappa,i), gmael(G,j,i));
            for (i=j+1; i<kappa; i++)
              addzi(&gmael(G,kappa,i), gmael(G,i,j));
            for (i=kappa+1; i<=maxG; i++)
              addzi(&gmael(G,i,kappa), gmael(G,i,j));
          }
          continue;
        }
        /* we have |X| >= 2 */
        if (tmp->e < BITS_IN_LONG-1)
        {
          if (tmp->d > 0)
          {
            ulong xx = (ulong) pari_rint(ldexp(tmp->d, tmp->e)); /* X fits in an ulong */
            for (k=zeros+1; k<j; k++)
              dpe_submuluz(Dmael(mu,kappa,k), Dmael(mu,j,k), xx, Dmael(mu,kappa,k));
            for (i=1; i<=n; i++)
              submulziu(&gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              submulziu(&gmael(U,kappa,i), gmael(U,j,i), xx);
            btop = avma;
            ztmp = submuliu2n(mulii(gmael(G,j,j), sqru(xx)), gmael(G,kappa,j), xx, 1);
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            affii_or_copy_gc(btop, ztmp, &gmael(G,kappa,kappa));
            for (i=1; i<=j; i++)
              submulziu(&gmael(G,kappa,i), gmael(G,j,i), xx);
            for (i=j+1; i<kappa; i++)
              submulziu(&gmael(G,kappa,i), gmael(G,i,j), xx);
            for (i=kappa+1; i<=maxG; i++)
              submulziu(&gmael(G,i,kappa), gmael(G,i,j), xx);
          }
          else
          {
            ulong xx = (ulong) pari_rint(ldexp(-tmp->d, tmp->e)); /* X fits in an ulong */
            for (k=zeros+1; k<j; k++)
              dpe_addmuluz(Dmael(mu,kappa,k), Dmael(mu,j,k), xx, Dmael(mu,kappa,k));
            for (i=1; i<=n; i++)
              addmulziu(&gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              addmulziu(&gmael(U,kappa,i), gmael(U,j,i), xx);
            btop = avma;
            ztmp = addmuliu2n(mulii(gmael(G,j,j), sqru(xx)), gmael(G,kappa,j), xx, 1);
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            affii_or_copy_gc(btop, ztmp, &gmael(G,kappa,kappa));
            for (i=1; i<=j; i++)
              addmulziu(&gmael(G,kappa,i), gmael(G,j,i), xx);
            for (i=j+1; i<kappa; i++)
              addmulziu(&gmael(G,kappa,i), gmael(G,i,j), xx);
            for (i=kappa+1; i<=maxG; i++)
              addmulziu(&gmael(G,i,kappa), gmael(G,i,j), xx);
          }
        }
        else
        {
          long e = tmp->e - BITS_IN_LONG + 1;
          if (tmp->d > 0)
          {
            ulong xx = (ulong) pari_rint(ldexp(tmp->d, BITS_IN_LONG - 1));
            for (k=zeros+1; k<j; k++)
            {
              dpe_t x;
              dpe_muluz(Dmael(mu,j,k), xx, &x);
              x.e += e;
              dpe_subz(Dmael(mu,kappa,k), &x, Dmael(mu,kappa,k));
            }
            for (i=1; i<=n; i++)
              submulzu2n(&gmael(B,kappa,i), gmael(B,j,i), xx, e);
            for (i=1; i<=d; i++)
              submulzu2n(&gmael(U,kappa,i), gmael(U,j,i), xx, e);
            btop = avma;
            ztmp = submuliu2n(mulshift(gmael(G,j,j), sqru(xx), 2*e),
                gmael(G,kappa,j), xx, e+1);
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            affii_or_copy_gc(btop, ztmp, &gmael(G,kappa,kappa));
            for (i=1; i<=j; i++)
              submulzu2n(&gmael(G,kappa,i), gmael(G,j,i), xx, e);
            for (   ; i<kappa; i++)
              submulzu2n(&gmael(G,kappa,i), gmael(G,i,j), xx, e);
            for (i=kappa+1; i<=maxG; i++)
              submulzu2n(&gmael(G,i,kappa), gmael(G,i,j), xx, e);
          } else
          {
            ulong xx = (ulong) pari_rint(ldexp(-tmp->d, BITS_IN_LONG - 1));
            for (k=zeros+1; k<j; k++)
            {
              dpe_t x;
              dpe_muluz(Dmael(mu,j,k), xx, &x);
              x.e += e;
              dpe_addz(Dmael(mu,kappa,k), &x, Dmael(mu,kappa,k));
            }
            for (i=1; i<=n; i++)
              addmulzu2n(&gmael(B,kappa,i), gmael(B,j,i), xx, e);
            for (i=1; i<=d; i++)
              addmulzu2n(&gmael(U,kappa,i), gmael(U,j,i), xx, e);
            btop = avma;
            ztmp = addmuliu2n(mulshift(gmael(G,j,j), sqru(xx), 2*e),
                gmael(G,kappa,j), xx, e+1);
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            affii_or_copy_gc(btop, ztmp, &gmael(G,kappa,kappa));
            for (i=1; i<=j; i++)
              addmulzu2n(&gmael(G,kappa,i), gmael(G,j,i), xx, e);
            for (   ; i<kappa; i++)
              addmulzu2n(&gmael(G,kappa,i), gmael(G,i,j), xx, e);
            for (i=kappa+1; i<=maxG; i++)
              addmulzu2n(&gmael(G,i,kappa), gmael(G,i,j), xx, e);
          }
        }
      }
    if (!go_on) break; /* Anything happened? */
    aa = zeros+1;
  }

  affidpe(gmael(G,kappa,kappa), Del(s,zeros+1));
  /* the last s[kappa-1]=r[kappa][kappa] is computed only if kappa increases */
  for (k=zeros+1; k<=kappa-2; k++)
    dpe_submulz(Del(s,k), Dmael(mu,kappa,k), Dmael(r,kappa,k), Del(s,k+1));
  *pG = G; *pB = B; *pU = U; return 0;
}

/* G integral Gram matrix, LLL-reduces (G,B,U) in place [apply base change
 * transforms to B and U]. If (keepfirst), never swap with first vector.
 * If G = NULL, we compute the Gram matrix incrementally.
 * Return -1 on failure, else zeros = dim Kernel (>= 0) */
static long
fplll_dpe(GEN *pG, GEN *pB, GEN *pU, GEN *pr, double DELTA, double ETA,
      long keepfirst)
{
  pari_sp av;
  GEN Gtmp, alpha, G = *pG, B = *pB, U = *pU;
  long d, maxG, kappa, kappa2, i, j, zeros, kappamax, incgram = !G, cnt = 0;
  dpe_t delta, eta, **mu, **r, *s;
  affdbldpe(DELTA,&delta);
  affdbldpe(ETA,&eta);

  if (incgram)
  { /* incremental Gram matrix */
    maxG = 2; d = lg(B)-1;
    G = zeromatcopy(d, d);
  }
  else
    maxG = d = lg(G)-1;

  mu = cget_dpemat(d+1);
  r  = cget_dpemat(d+1);
  s  = cget_dpevec(d+1);
  for (j = 1; j <= d; j++)
  {
    mu[j]= cget_dpevec(d+1);
    r[j] = cget_dpevec(d+1);
  }
  Gtmp = cgetg(d+1, t_VEC);
  alpha = cgetg(d+1, t_VECSMALL);
  av = avma;

  /* Step2: Initializing the main loop */
  kappamax = 1;
  i = 1;
  do {
    if (incgram) gmael(G,i,i) = ZV_dotsquare(gel(B,i));
    affidpe(gmael(G,i,i), Dmael(r,i,i));
  } while (!signe(gmael(G,i,i)) && ++i <= d);
  zeros = i-1; /* all basis vectors b_i with i <= zeros are zero vectors */
  kappa = i;
  for (i=zeros+1; i<=d; i++) alpha[i]=1;

  while (++kappa <= d)
  {
    if (kappa > kappamax)
    {
      if (DEBUGLEVEL>=4) err_printf("K%ld ",kappa);
      kappamax = kappa;
      if (incgram)
      {
        for (i=zeros+1; i<=kappa; i++)
          gmael(G,kappa,i) = ZV_dotproduct(gel(B,kappa), gel(B,i));
        maxG = kappamax;
      }
    }
    /* Step3: Call to the Babai algorithm, mu,r,s updated in place */
    if (Babai_dpe(av, kappa, &G,&B,&U, mu,r,s, alpha[kappa], zeros, maxG, &eta))
    { *pG = incgram? NULL: G; *pB = B; *pU = U; return -1; }
    if ((keepfirst && kappa == 2) ||
        dpe_cmpmul(Dmael(r,kappa-1,kappa-1), &delta, Del(s,kappa-1)) <= 0)
    { /* Step4: Success of Lovasz's condition */
      alpha[kappa] = kappa;
      dpe_submulz(Del(s,kappa-1), Dmael(mu,kappa,kappa-1), Dmael(r,kappa,kappa-1), Dmael(r,kappa,kappa));
      continue;
    }
    /* Step5: Find the right insertion index kappa, kappa2 = initial kappa */
    if (DEBUGLEVEL>=4 && kappa==kappamax && Del(s,kappa-1)->d)
      if (++cnt > 20) { cnt = 0; err_printf("(%ld) ", Del(s,1)->e-1); }
    kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros+2 + (keepfirst ? 1: 0)) break;
    } while (dpe_cmpmul(Dmael(r,kappa-1,kappa-1), &delta, Del(s,kappa-1)) >= 0);
    update_alpha(alpha, kappa, kappa2, kappamax);

    /* Step6: Update the mu's and r's */
    dperotate(mu, kappa2, kappa);
    dperotate(r, kappa2, kappa);
    affdpe(Del(s,kappa), Dmael(r,kappa,kappa));

    /* Step7: Update G, B, U */
    if (U) rotate(U, kappa2, kappa);
    if (B) rotate(B, kappa2, kappa);
    rotateG(G,kappa2,kappa, maxG, Gtmp);

    /* Step8: Prepare the next loop iteration */
    if (kappa == zeros+1 && !signe(gmael(G,kappa,kappa)))
    {
      zeros++; kappa++;
      affidpe(gmael(G,kappa,kappa), Dmael(r,kappa,kappa));
    }
  }
  if (pr) *pr = dpeM_diagonal_shallow(r,d);
  *pG = G; *pB = B; *pU = U; return zeros; /* success */
}


/************************** PROVED version (t_INT) *************************/

/* Babai's Nearest Plane algorithm (iterative).
 * Size-reduces b_kappa using mu_{i,j} and r_{i,j} for j<=i <kappa
 * Update B[,kappa]; compute mu_{kappa,j}, r_{kappa,j} for j<=kappa and s[kappa]
 * mu, r, s updated in place (affrr). Return 1 on failure, else 0. */
static int
Babai(pari_sp av, long kappa, GEN *pG, GEN *pB, GEN *pU, GEN mu, GEN r, GEN s,
      long a, long zeros, long maxG, GEN eta, long prec)
{
  GEN G = *pG, B = *pB, U = *pU, ztmp;
  long k, aa = a > zeros? a: zeros+1;
  const long n = B? nbrows(B): 0, d = U ? lg(U)-1: 0, bit = prec2nbits(prec);
  long emaxmu = EX0, emax2mu = EX0;
  /* N.B: we set d = 0 (resp. n = 0) to avoid updating U (resp. B) */

  for (;;) {
    int go_on = 0;
    long i, j, emax3mu = emax2mu;

    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Babai[1], a=%ld", aa);
      gc_lll(av,3,&G,&B,&U);
    }
    /* Step2: compute the GSO for stage kappa */
    emax2mu = emaxmu; emaxmu = EX0;
    for (j=aa; j<kappa; j++)
    {
      pari_sp btop = avma;
      GEN g = gmael(G,kappa,j);
      for (k = zeros+1; k < j; k++)
        g = mpsub(g, mulrr(gmael(mu,j,k), gmael(r,kappa,k)));
      mpaff(g, gmael(r,kappa,j));
      affrr(divrr(gmael(r,kappa,j), gmael(r,j,j)), gmael(mu,kappa,j));
      emaxmu = maxss(emaxmu, expo(gmael(mu,kappa,j)));
      set_avma(btop);
    }
    if (emax3mu != EX0 && emax3mu <= emax2mu + 5) /* precision too low */
    { *pG = G; *pB = B; *pU = U; return 1; }

    for (j=kappa-1; j>zeros; j--)
      if (abscmprr(gmael(mu,kappa,j), eta) > 0) { go_on = 1; break; }

    /* Step3--5: compute the X_j's  */
    if (go_on)
      for (j=kappa-1; j>zeros; j--)
      {
        pari_sp btop;
        GEN tmp = gmael(mu,kappa,j);
        if (absrsmall(tmp)) continue; /* size-reduced */

        if (gc_needed(av,2))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"Babai[2], a=%ld, j=%ld", aa,j);
          gc_lll(av,3,&G,&B,&U);
        }
        btop = avma;
        /* we consider separately the case |X| = 1 */
        if (absrsmall2(tmp))
        {
          if (signe(tmp) > 0) { /* in this case, X = 1 */
            for (k=zeros+1; k<j; k++)
              affrr(subrr(gmael(mu,kappa,k), gmael(mu,j,k)), gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = subii(gmael(B,kappa,i), gmael(B,j,i));
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = subii(gmael(U,kappa,i), gmael(U,j,i));
            btop = avma;
            ztmp = subii(gmael(G,j,j), shifti(gmael(G,kappa,j), 1));
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
            for (i=1; i<=j; i++)
              gmael(G,kappa,i) = subii(gmael(G,kappa,i), gmael(G,j,i));
            for (i=j+1; i<kappa; i++)
              gmael(G,kappa,i) = subii(gmael(G,kappa,i), gmael(G,i,j));
            for (i=kappa+1; i<=maxG; i++)
              gmael(G,i,kappa) = subii(gmael(G,i,kappa), gmael(G,i,j));
          } else { /* otherwise X = -1 */
            for (k=zeros+1; k<j; k++)
              affrr(addrr(gmael(mu,kappa,k), gmael(mu,j,k)), gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = addii(gmael(B,kappa,i),gmael(B,j,i));
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = addii(gmael(U,kappa,i),gmael(U,j,i));
            btop = avma;
            ztmp = addii(gmael(G,j,j), shifti(gmael(G,kappa,j), 1));
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
            for (i=1; i<=j; i++)
              gmael(G,kappa,i) = addii(gmael(G,kappa,i), gmael(G,j,i));
            for (i=j+1; i<kappa; i++)
              gmael(G,kappa,i) = addii(gmael(G,kappa,i), gmael(G,i,j));
            for (i=kappa+1; i<=maxG; i++)
              gmael(G,i,kappa) = addii(gmael(G,i,kappa), gmael(G,i,j));
          }
          continue;
        }
        /* we have |X| >= 2 */
        if (expo(tmp) < BITS_IN_LONG)
        {
          ulong xx = roundr_safe(tmp)[2]; /* X fits in an ulong */
          if (signe(tmp) > 0) /* = xx */
          {
            for (k=zeros+1; k<j; k++)
              affrr(subrr(gmael(mu,kappa,k), mulur(xx, gmael(mu,j,k))),
                  gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = submuliu_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = submuliu_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
            btop = avma;
            ztmp = submuliu2n(mulii(gmael(G,j,j), sqru(xx)), gmael(G,kappa,j), xx, 1);
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
            for (i=1; i<=j; i++)
              gmael(G,kappa,i) = submuliu_inplace(gmael(G,kappa,i), gmael(G,j,i), xx);
            for (i=j+1; i<kappa; i++)
              gmael(G,kappa,i) = submuliu_inplace(gmael(G,kappa,i), gmael(G,i,j), xx);
            for (i=kappa+1; i<=maxG; i++)
              gmael(G,i,kappa) = submuliu_inplace(gmael(G,i,kappa), gmael(G,i,j), xx);
          }
          else /* = -xx */
          {
            for (k=zeros+1; k<j; k++)
              affrr(addrr(gmael(mu,kappa,k), mulur(xx, gmael(mu,j,k))),
                  gmael(mu,kappa,k));
            set_avma(btop);
            for (i=1; i<=n; i++)
              gmael(B,kappa,i) = addmuliu_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
            for (i=1; i<=d; i++)
              gmael(U,kappa,i) = addmuliu_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
            btop = avma;
            ztmp = addmuliu2n(mulii(gmael(G,j,j), sqru(xx)), gmael(G,kappa,j), xx, 1);
            ztmp = addii(gmael(G,kappa,kappa), ztmp);
            gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
            for (i=1; i<=j; i++)
              gmael(G,kappa,i) = addmuliu_inplace(gmael(G,kappa,i), gmael(G,j,i), xx);
            for (i=j+1; i<kappa; i++)
              gmael(G,kappa,i) = addmuliu_inplace(gmael(G,kappa,i), gmael(G,i,j), xx);
            for (i=kappa+1; i<=maxG; i++)
              gmael(G,i,kappa) = addmuliu_inplace(gmael(G,i,kappa), gmael(G,i,j), xx);
          }
        }
        else
        {
          long e;
          GEN X = truncexpo(tmp, bit, &e); /* tmp ~ X * 2^e */
          btop = avma;
          for (k=zeros+1; k<j; k++)
          {
            GEN x = mulir(X, gmael(mu,j,k));
            if (e) shiftr_inplace(x, e);
            affrr(subrr(gmael(mu,kappa,k), x), gmael(mu,kappa,k));
          }
          set_avma(btop);
          for (i=1; i<=n; i++)
            gmael(B,kappa,i) = submulshift(gmael(B,kappa,i), gmael(B,j,i), X, e);
          for (i=1; i<=d; i++)
            gmael(U,kappa,i) = submulshift(gmael(U,kappa,i), gmael(U,j,i), X, e);
          btop = avma;
          ztmp = submulshift(mulshift(gmael(G,j,j), sqri(X), 2*e),
              gmael(G,kappa,j), X, e+1);
          ztmp = addii(gmael(G,kappa,kappa), ztmp);
          gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
          for (i=1; i<=j; i++)
            gmael(G,kappa,i) = submulshift(gmael(G,kappa,i), gmael(G,j,i), X, e);
          for (   ; i<kappa; i++)
            gmael(G,kappa,i) = submulshift(gmael(G,kappa,i), gmael(G,i,j), X, e);
          for (i=kappa+1; i<=maxG; i++)
            gmael(G,i,kappa) = submulshift(gmael(G,i,kappa), gmael(G,i,j), X, e);
        }
      }
    if (!go_on) break; /* Anything happened? */
    aa = zeros+1;
  }

  affir(gmael(G,kappa,kappa), gel(s,zeros+1));
  /* the last s[kappa-1]=r[kappa][kappa] is computed only if kappa increases */
  av = avma;
  for (k=zeros+1; k<=kappa-2; k++)
    affrr(subrr(gel(s,k), mulrr(gmael(mu,kappa,k), gmael(r,kappa,k))),
          gel(s,k+1));
  *pG = G; *pB = B; *pU = U; return gc_bool(av, 0);
}

/* G integral Gram matrix, LLL-reduces (G,B,U) in place [apply base change
 * transforms to B and U]. If (keepfirst), never swap with first vector.
 * If G = NULL, we compute the Gram matrix incrementally.
 * Return -1 on failure, else zeros = dim Kernel (>= 0) */
static long
fplll(GEN *pG, GEN *pB, GEN *pU, GEN *pr, double DELTA, double ETA,
      long keepfirst, long prec)
{
  pari_sp av, av2;
  GEN mu, r, s, tmp, Gtmp, alpha, G = *pG, B = *pB, U = *pU;
  GEN delta = dbltor(DELTA), eta = dbltor(ETA);
  long d, maxG, kappa, kappa2, i, j, zeros, kappamax, incgram = !G, cnt = 0;

  if (incgram)
  { /* incremental Gram matrix */
    maxG = 2; d = lg(B)-1;
    G = zeromatcopy(d, d);
  }
  else
    maxG = d = lg(G)-1;

  mu = cgetg(d+1, t_MAT);
  r  = cgetg(d+1, t_MAT);
  s  = cgetg(d+1, t_VEC);
  for (j = 1; j <= d; j++)
  {
    GEN M = cgetg(d+1, t_COL), R = cgetg(d+1, t_COL);
    gel(mu,j)= M;
    gel(r,j) = R;
    gel(s,j) = cgetr(prec);
    for (i = 1; i <= d; i++)
    {
      gel(R,i) = cgetr(prec);
      gel(M,i) = cgetr(prec);
    }
  }
  Gtmp = cgetg(d+1, t_VEC);
  alpha = cgetg(d+1, t_VECSMALL);
  av = avma;

  /* Step2: Initializing the main loop */
  kappamax = 1;
  i = 1;
  do {
    if (incgram) gmael(G,i,i) = ZV_dotsquare(gel(B,i));
    affir(gmael(G,i,i), gmael(r,i,i));
  } while (!signe(gmael(G,i,i)) && ++i <= d);
  zeros = i-1; /* all basis vectors b_i with i <= zeros are zero vectors */
  kappa = i;
  for (i=zeros+1; i<=d; i++) alpha[i]=1;

  while (++kappa <= d)
  {
    if (kappa > kappamax)
    {
      if (DEBUGLEVEL>=4) err_printf("K%ld ",kappa);
      kappamax = kappa;
      if (incgram)
      {
        for (i=zeros+1; i<=kappa; i++)
          gmael(G,kappa,i) = ZV_dotproduct(gel(B,kappa), gel(B,i));
        maxG = kappamax;
      }
    }
    /* Step3: Call to the Babai algorithm, mu,r,s updated in place */
    if (Babai(av, kappa, &G,&B,&U, mu,r,s, alpha[kappa], zeros, maxG, eta, prec))
    { *pG = incgram? NULL: G; *pB = B; *pU = U; return -1; }
    av2 = avma;
    if ((keepfirst && kappa == 2) ||
        cmprr(mulrr(gmael(r,kappa-1,kappa-1), delta), gel(s,kappa-1)) <= 0)
    { /* Step4: Success of Lovasz's condition */
      alpha[kappa] = kappa;
      tmp = mulrr(gmael(mu,kappa,kappa-1), gmael(r,kappa,kappa-1));
      affrr(subrr(gel(s,kappa-1), tmp), gmael(r,kappa,kappa));
      set_avma(av2); continue;
    }
    /* Step5: Find the right insertion index kappa, kappa2 = initial kappa */
    if (DEBUGLEVEL>=4 && kappa==kappamax && signe(gel(s,kappa-1)))
      if (++cnt > 20) { cnt = 0; err_printf("(%ld) ", expo(gel(s,1))); }
    kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros+2 + (keepfirst ? 1: 0)) break;
      tmp = mulrr(gmael(r,kappa-1,kappa-1), delta);
    } while (cmprr(gel(s,kappa-1), tmp) <= 0);
    set_avma(av2);
    update_alpha(alpha, kappa, kappa2, kappamax);

    /* Step6: Update the mu's and r's */
    rotate(mu, kappa2, kappa);
    rotate(r, kappa2, kappa);
    affrr(gel(s,kappa), gmael(r,kappa,kappa));

    /* Step7: Update G, B, U */
    if (U) rotate(U, kappa2, kappa);
    if (B) rotate(B, kappa2, kappa);
    rotateG(G,kappa2,kappa, maxG, Gtmp);

    /* Step8: Prepare the next loop iteration */
    if (kappa == zeros+1 && !signe(gmael(G,kappa,kappa)))
    {
      zeros++; kappa++;
      affir(gmael(G,kappa,kappa), gmael(r,kappa,kappa));
    }
  }
  if (pr) *pr = RgM_diagonal_shallow(r);
  *pG = G; *pB = B; *pU = U; return zeros; /* success */
}

/* Assume x a ZM, if pN != NULL, set it to Gram-Schmidt (squared) norms
 * The following modes are supported:
 * - flag & LLL_INPLACE: x a lattice basis, return x*U
 * - flag & LLL_GRAM: x a Gram matrix / else x a lattice basis; return
 *     LLL base change matrix U [LLL_IM]
 *     kernel basis [LLL_KER, nonreduced]
 *     both [LLL_ALL] */
GEN
ZM_lll_norms(GEN x, double DELTA, long flag, GEN *pN)
{
  pari_sp av = avma;
  const double ETA = 0.51;
  const long keepfirst = flag & LLL_KEEP_FIRST;
  long p, zeros = -1, n = lg(x)-1;
  GEN G, B, U;
  pari_timer T;
  if (n <= 1) return lll_trivial(x, flag);
  if (nbrows(x)==0)
  {
    if (flag & LLL_KER) return matid(n);
    if (flag & (LLL_INPLACE|LLL_IM)) return cgetg(1,t_MAT);
    retmkvec2(matid(n), cgetg(1,t_MAT));
  }
  x = gcopy(x);
  if (flag & LLL_GRAM)
  { G = x; B = NULL; U = matid(n); }
  else
  { G = NULL; B = x; U = (flag & LLL_INPLACE)? NULL: matid(n); }
  if(DEBUGLEVEL>=4) timer_start(&T);
  if (!(flag & LLL_GRAM))
  {
    long t;
    if(DEBUGLEVEL>=4)
      err_printf("Entering L^2 (double): LLL-parameters (%.3f,%.3f)\n",
                 DELTA,ETA);
    zeros = fplll_fast(&B, &U, DELTA, ETA, keepfirst);
    if (DEBUGLEVEL>=4) timer_printf(&T, zeros < 0? "LLL (failed)": "LLL");
    for (p = DEFAULTPREC, t = 0; zeros < 0 && t < 1 ; p += EXTRAPREC64, t++)
    {
      if (DEBUGLEVEL>=4)
        err_printf("Entering L^2 (heuristic): LLL-parameters (%.3f,%.3f), prec = %d/%d\n", DELTA, ETA, p, p);
      zeros = fplll_heuristic(&B, &U, DELTA, ETA, keepfirst, p, p);
      gc_lll(av, 2, &B, &U);
      if (DEBUGLEVEL>=4) timer_printf(&T, zeros < 0? "LLL (failed)": "LLL");
    }
  }
  if(DEBUGLEVEL>=4)
    err_printf("Entering L^2 (dpe): LLL-parameters (%.3f,%.3f)\n", DELTA,ETA);
  zeros = fplll_dpe(&G, &B, &U, pN, DELTA, ETA, keepfirst);
  if (DEBUGLEVEL>=4) timer_printf(&T, zeros < 0? "LLL (failed)": "LLL");
  if (zeros < 0)
    for (p = DEFAULTPREC;; p += EXTRAPREC64)
    {
      if (DEBUGLEVEL>=4)
        err_printf("Entering L^2: LLL-parameters (%.3f,%.3f), prec = %d\n",
                    DELTA,ETA, p);
      zeros = fplll(&G, &B, &U, pN, DELTA, ETA, keepfirst, p);
      if (DEBUGLEVEL>=4) timer_printf(&T, zeros < 0? "LLL (failed)": "LLL");
      if (zeros >= 0) break;
      gc_lll(av, 3, &G, &B, &U);
    }
  return lll_finish(U? U: B, zeros, flag);
}

/********************************************************************/
/**                                                                **/
/**                        LLL OVER K[X]                           **/
/**                                                                **/
/********************************************************************/
static int
pslg(GEN x)
{
  long tx;
  if (gequal0(x)) return 2;
  tx = typ(x); return is_scalar_t(tx)? 3: lg(x);
}

static int
REDgen(long k, long l, GEN h, GEN L, GEN B)
{
  GEN q, u = gcoeff(L,k,l);
  long i;

  if (pslg(u) < pslg(B)) return 0;

  q = gneg(gdeuc(u,B));
  gel(h,k) = gadd(gel(h,k), gmul(q,gel(h,l)));
  for (i=1; i<l; i++) gcoeff(L,k,i) = gadd(gcoeff(L,k,i), gmul(q,gcoeff(L,l,i)));
  gcoeff(L,k,l) = gadd(gcoeff(L,k,l), gmul(q,B)); return 1;
}

static int
do_SWAPgen(GEN h, GEN L, GEN B, long k, GEN fl, int *flc)
{
  GEN p1, la, la2, Bk;
  long ps1, ps2, i, j, lx;

  if (!fl[k-1]) return 0;

  la = gcoeff(L,k,k-1); la2 = gsqr(la);
  Bk = gel(B,k);
  if (fl[k])
  {
    GEN q = gadd(la2, gmul(gel(B,k-1),gel(B,k+1)));
    ps1 = pslg(gsqr(Bk));
    ps2 = pslg(q);
    if (ps1 <= ps2 && (ps1 < ps2 || !*flc)) return 0;
    *flc = (ps1 != ps2);
    gel(B,k) = gdiv(q, Bk);
  }

  swap(gel(h,k-1), gel(h,k)); lx = lg(L);
  for (j=1; j<k-1; j++) swap(gcoeff(L,k-1,j), gcoeff(L,k,j));
  if (fl[k])
  {
    for (i=k+1; i<lx; i++)
    {
      GEN t = gcoeff(L,i,k);
      p1 = gsub(gmul(gel(B,k+1),gcoeff(L,i,k-1)), gmul(la,t));
      gcoeff(L,i,k) = gdiv(p1, Bk);
      p1 = gadd(gmul(la,gcoeff(L,i,k-1)), gmul(gel(B,k-1),t));
      gcoeff(L,i,k-1) = gdiv(p1, Bk);
    }
  }
  else if (!gequal0(la))
  {
    p1 = gdiv(la2, Bk);
    gel(B,k+1) = gel(B,k) = p1;
    for (i=k+2; i<=lx; i++) gel(B,i) = gdiv(gmul(p1,gel(B,i)),Bk);
    for (i=k+1; i<lx; i++)
      gcoeff(L,i,k-1) = gdiv(gmul(la,gcoeff(L,i,k-1)), Bk);
    for (j=k+1; j<lx-1; j++)
      for (i=j+1; i<lx; i++)
        gcoeff(L,i,j) = gdiv(gmul(p1,gcoeff(L,i,j)), Bk);
  }
  else
  {
    gcoeff(L,k,k-1) = gen_0;
    for (i=k+1; i<lx; i++)
    {
      gcoeff(L,i,k) = gcoeff(L,i,k-1);
      gcoeff(L,i,k-1) = gen_0;
    }
    gel(B,k) = gel(B,k-1); fl[k] = 1; fl[k-1] = 0;
  }
  return 1;
}

static void
incrementalGSgen(GEN x, GEN L, GEN B, long k, GEN fl)
{
  GEN u = NULL; /* gcc -Wall */
  long i, j;
  for (j = 1; j <= k; j++)
    if (j==k || fl[j])
    {
      u = gcoeff(x,k,j);
      if (!is_extscalar_t(typ(u))) pari_err_TYPE("incrementalGSgen",u);
      for (i=1; i<j; i++)
        if (fl[i])
        {
          u = gsub(gmul(gel(B,i+1),u), gmul(gcoeff(L,k,i),gcoeff(L,j,i)));
          u = gdiv(u, gel(B,i));
        }
      gcoeff(L,k,j) = u;
    }
  if (gequal0(u)) gel(B,k+1) = gel(B,k);
  else
  {
    gel(B,k+1) = gcoeff(L,k,k); gcoeff(L,k,k) = gen_1; fl[k] = 1;
  }
}

static GEN
lllgramallgen(GEN x, long flag)
{
  long lx = lg(x), i, j, k, l, n;
  pari_sp av;
  GEN B, L, h, fl;
  int flc;

  n = lx-1; if (n<=1) return lll_trivial(x,flag);
  if (lgcols(x) != lx) pari_err_DIM("lllgramallgen");

  fl = cgetg(lx, t_VECSMALL);

  av = avma;
  B = scalarcol_shallow(gen_1, lx);
  L = cgetg(lx,t_MAT);
  for (j=1; j<lx; j++) { gel(L,j) = zerocol(n); fl[j] = 0; }

  h = matid(n);
  for (i=1; i<lx; i++)
    incrementalGSgen(x, L, B, i, fl);
  flc = 0;
  for(k=2;;)
  {
    if (REDgen(k, k-1, h, L, gel(B,k))) flc = 1;
    if (do_SWAPgen(h, L, B, k, fl, &flc)) { if (k > 2) k--; }
    else
    {
      for (l=k-2; l>=1; l--)
        if (REDgen(k, l, h, L, gel(B,l+1))) flc = 1;
      if (++k > n) break;
    }
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lllgramallgen");
      gerepileall(av,3,&B,&L,&h);
    }
  }
  k=1; while (k<lx && !fl[k]) k++;
  return lll_finish(h,k-1,flag);
}

static int
RgM_square(GEN x) { long l = lg(x); return l == 1 || l == lgcols(x); }
static GEN
lllallgen(GEN x, long flag)
{
  pari_sp av = avma;
  if (!(flag & LLL_GRAM)) x = gram_matrix(x);
  else if (!RgM_square(x)) pari_err_DIM("qflllgram");
  return gerepilecopy(av, lllgramallgen(x, flag));
}
GEN
lllgen(GEN x) { return lllallgen(x, LLL_IM); }
GEN
lllkerimgen(GEN x) { return lllallgen(x, LLL_ALL); }
GEN
lllgramgen(GEN x)  { return lllallgen(x, LLL_IM|LLL_GRAM); }
GEN
lllgramkerimgen(GEN x)  { return lllallgen(x, LLL_ALL|LLL_GRAM); }

static GEN
lllall(GEN x, long flag)
{ pari_sp av = avma; return gerepilecopy(av, ZM_lll(x, LLLDFT, flag)); }
GEN
lllint(GEN x) { return lllall(x, LLL_IM); }
GEN
lllkerim(GEN x) { return lllall(x, LLL_ALL); }
GEN
lllgramint(GEN x)
{ if (!RgM_square(x)) pari_err_DIM("qflllgram");
  return lllall(x, LLL_IM | LLL_GRAM); }
GEN
lllgramkerim(GEN x)
{ if (!RgM_square(x)) pari_err_DIM("qflllgram");
  return lllall(x, LLL_ALL | LLL_GRAM); }

GEN
lllfp(GEN x, double D, long flag)
{
  long n = lg(x)-1;
  pari_sp av = avma;
  GEN h;
  if (n <= 1) return lll_trivial(x,flag);
  if ((flag & LLL_GRAM) && !RgM_square(x)) pari_err_DIM("qflllgram");
  h = ZM_lll(RgM_rescale_to_int(x), D, flag);
  return gerepilecopy(av, h);
}

GEN
lllgram(GEN x) { return lllfp(x,LLLDFT,LLL_GRAM|LLL_IM); }
GEN
lll(GEN x) { return lllfp(x,LLLDFT,LLL_IM); }

GEN
qflll0(GEN x, long flag)
{
  if (typ(x) != t_MAT) pari_err_TYPE("qflll",x);
  switch(flag)
  {
    case 0: return lll(x);
    case 1: RgM_check_ZM(x,"qflll"); return lllint(x);
    case 2: RgM_check_ZM(x,"qflll"); return lllintpartial(x);
    case 3: RgM_check_ZM(x,"qflll"); return lllall(x, LLL_INPLACE);
    case 4: RgM_check_ZM(x,"qflll"); return lllkerim(x);
    case 5: return lllkerimgen(x);
    case 8: return lllgen(x);
    default: pari_err_FLAG("qflll");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

GEN
qflllgram0(GEN x, long flag)
{
  if (typ(x) != t_MAT) pari_err_TYPE("qflllgram",x);
  switch(flag)
  {
    case 0: return lllgram(x);
    case 1: RgM_check_ZM(x,"qflllgram"); return lllgramint(x);
    case 4: RgM_check_ZM(x,"qflllgram"); return lllgramkerim(x);
    case 5: return lllgramkerimgen(x);
    case 8: return lllgramgen(x);
    default: pari_err_FLAG("qflllgram");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/********************************************************************/
/**                                                                **/
/**                   INTEGRAL KERNEL (LLL REDUCED)                **/
/**                                                                **/
/********************************************************************/
static GEN
kerint0(GEN M)
{
  /* return ZM_lll(M, LLLDFT, LLL_KER); */
  GEN U, H = ZM_hnflll(M,&U,1);
  long d = lg(M)-lg(H);
  if (!d) return cgetg(1, t_MAT);
  return ZM_lll(vecslice(U,1,d), LLLDFT, LLL_INPLACE);
}
GEN
kerint(GEN M)
{
  pari_sp av = avma;
  return gerepilecopy(av, kerint0(M));
}
/* OBSOLETE: use kerint */
GEN
matkerint0(GEN M, long flag)
{
  pari_sp av = avma;
  if (typ(M) != t_MAT) pari_err_TYPE("matkerint",M);
  M = Q_primpart(M);
  RgM_check_ZM(M, "kerint");
  switch(flag)
  {
    case 0:
    case 1: return gerepilecopy(av, kerint0(M));
    default: pari_err_FLAG("matkerint");
  }
  return NULL; /* LCOV_EXCL_LINE */
}
