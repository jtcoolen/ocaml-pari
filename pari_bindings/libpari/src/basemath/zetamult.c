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

/********************************************************************/
/**                                                                **/
/**                         MULTIPLE ZETA VALUES                   **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_zetamult

/********************************************************************/
/**                           CONVERSIONS                          **/
/********************************************************************/
/* vecsmall to binary */
static long
fd(GEN evec, long a, long b)
{
  long i, s = 0;
  for (i = a; i <= b; i++) s = evec[i] | (s << 1);
  return s;
}
/* 2^(b-a+1) + fd(evec) */
static long
fd1(GEN evec, long a, long b)
{
  long i, s = 1;
  for (i = a; i <= b; i++) s = evec[i] | (s << 1);
  return s;
}

/* m > 0 */
static GEN
mtoevec(GEN m)
{
  GEN e = vecsmall_append(binary_zv(m), 1);
  e[1] = 0; return e;
}
static GEN
etoindex(GEN evec) { return utoipos(fd1(evec, 2, lg(evec)-2)); }

/* dual of evec[1..l-1] */
static GEN
revslice(GEN evec, long l)
{
  GEN res = cgetg(l, t_VECSMALL);
  long i;
  for (i = 1; i < l; i++) res[i] = 1 - evec[l-i];
  return res;
}

/* N.B. evec[ne] = 1 */
static GEN
etoa(GEN evec)
{
  long c = 0, cold = 0, i = 1, l = lg(evec);
  GEN avec = cgetg(l, t_VECSMALL);
  while (++c < l)
    if (evec[c] == 1) { avec[i++] = c - cold; cold = c; }
  setlg(avec, i); return avec;
}

static GEN
atoe(GEN avec)
{
  long i, j, l = lg(avec), n = zv_sum(avec);
  GEN evec = zero_zv(n);
  for (i = 1, j = 0; i < l; i++) { long a = avec[i]; j += a; evec[j] = 1; }
  return evec;
}


/* Conversions: types are evec, avec, m (if evec=0y1, m=(1y)_2).
   fl is respectively 0, 1, 2. Type of a is autodetected. */
static GEN
zetamultconvert_i(GEN a, long fl)
{
  long i, l;
  if (fl < 0 || fl > 2) pari_err_FLAG("zetamultconvert");
  switch(typ(a))
  {
    case t_INT:
      if (signe(a) <= 0) pari_err_TYPE("zetamultconvert",a);
      switch (fl)
      {
        case 0: a = mtoevec(a); break;
        case 1: a = etoa(mtoevec(a)); break;
        case 2: a = icopy(a); break;
      }
      break;
    case t_VEC: case t_COL: case t_VECSMALL:
      a = gtovecsmall(a);
      l = lg(a);
      if (a[1] == 0)
      {
        if (!a[l-1]) pari_err_TYPE("zetamultconvert", a);
        for (i = 1; i < l; i++)
          if (a[i] & ~1UL) pari_err_TYPE("zetamultconvert", a);
        switch (fl)
        {
          case 1: a = etoa(a); break;
          case 2: a = etoindex(a);
        }
      }
      else
      {
        if (a[1] < 2) pari_err_TYPE("zetamultconvert", a);
        for (i = 2; i < l; i++)
          if (a[i] <= 0) pari_err_TYPE("zetamultconvert", a);
        switch (fl)
        {
          case 0: a = atoe(a); break;
          case 2: a = etoindex(atoe(a));
        }
      }
      break;
    default: pari_err_TYPE("zetamultconvert", a);
  }
  return a;
}
GEN
zetamultconvert(GEN a, long fl)
{
  pari_sp av = avma;
  return gerepileuptoleaf(av, zetamultconvert_i(a, fl));
}

GEN
zetamultdual(GEN s)
{
  pari_sp av = avma;
  GEN e = zetamultconvert_i(s, 0);
  return gerepileuptoleaf(av, etoa(revslice(e, lg(e))));
}

/********************************************************************/
/**                  AVEC -> LIST OF STAR AVECS                    **/
/********************************************************************/
/* all star avecs corresponding to given t_VECSMALL avec */
static GEN
allstar(GEN avec)
{
  long i, la = lg(avec), K = 1 << (la - 2);
  GEN w = cgetg(K + 1, t_VEC);

  gel(w, 1) = avec;
  for (i = 2; i < la; i++)
  {
    long L = 1 << (i - 2), j;
    for (j = 1; j <= L; j++)
    {
      GEN u = gel(w,j), v;
      long k, l = lg(u) - 1, ind = l - la + i;
      gel(w, L + j) = v = cgetg(l, t_VECSMALL);
      for (k = 1; k < ind; k++) v[k] = u[k];
      v[ind] = u[ind] + u[ind + 1];
      for (k = ind + 1; k < l; k++) v[k] = u[k+1];
    }
  }
  return w;
}
/* same for multipolylogs */
static GEN
allstar2(GEN avec, GEN zvec)
{
  long la = lg(avec), i, K = 1 << (la - 2);
  GEN W = cgetg(K + 1, t_VEC), Z = cgetg(K + 1, t_VEC);

  gel(W, 1) = avec;
  gel(Z, 1) = zvec = gtovec(zvec);
  for (i = 2; i < la; i++)
  {
    long L = 1 << (i - 2), j;
    for (j = 1; j <= L; j++)
    {
      GEN u = gel(W, j), w, y = gel(Z, j), z;
      long l = lg(u) - 1, ind = l - la + i, k;
      w = cgetg(l, t_VECSMALL);
      z = cgetg(l, t_VEC);
      for (k = 1; k < ind; k++) { w[k] = u[k]; gel(z, k) = gel(y, k); }
      w[ind] = u[ind] + u[ind + 1];
      gel(z, ind) = gmul(gel(y, ind), gel(y, ind + 1));
      for (k = ind + 1; k < l; k++) { w[k] = u[k+1]; gel(z, k) = gel(y, k+1); }
      gel(W, L + j) = w;
      gel(Z, L + j) = z;
    }
  }
  return mkvec2(W, Z);
}

/**************************************************************/
/*              ZAGIER & RADCHENKO'S ALGORITHM                */
/**************************************************************/
/* accuracy 2^(-b); s << (b/log b)^2, l << b/sqrt(log b) */
static void
zparams(long *s, long *l, long b)
{
  double D = b * LOG10_2, E = 3*D/2 / log(3*D);
  *s = (long)floor(E*E);
  *l = (long)floor(sqrt(*s * log((double)*s)/2.));
}

static GEN
zetamult_Zagier(GEN avec, long bit, long prec)
{
  pari_sp av;
  GEN ze, z = NULL, b;
  long h, r, n, s, l, t = lg(avec) - 1;

  zparams(&s, &l, bit);
  ze= cgetg(s + 1, t_VEC);
  b = cgetg(l + 1, t_VEC);
  for (r = 1; r <= s; r++) gel(ze,r) = cgetr(prec);
  for (r = 1; r <= l; r++) { gel(b,r) = cgetr(prec); affur(0,gel(b,r)); }
  affur(1, gel(b,1)); av = avma;
  for (r = 1, h = -1; r <= t; r++)
  {
    long m, j = avec[r];
    GEN q;

    h += j - 1; z = gen_0;
    q = h? invr(itor(powuu(s,h), prec)): real_1(prec);
    for (m = 1; m <= l; m++)
    {
      pari_sp av2;
      GEN S = gel(b, m), C;
      q = divru(q, s); av2 = avma;
      C = binomialuu(h + m, m);
      for (n = 1; n < m; n++)
      { /* C = binom(h+m, m-n+1) */
        S = gsub(S, mulir(C, gel(b, n)));
        C = diviuexact(muliu(C, m-n+1), h+n);
      }
      affrr(divru(S, h+m), gel(b,m)); set_avma(av2);
      z = gadd(z, gmul(gel(b, m), q));
    }
    for (m = s; m >= 1; m--)
    {
      GEN z1 = r == 1? ginv(powuu(m,j)): gdiv(gel(ze, m), powuu(m,j));
      z1 = gadd(z, z1);
      affrr(z, gel(ze, m)); z = z1;
    }
    set_avma(av);
  }
  return z;
}

/* Compute t-mzvs; slower than Zagier's code for t=0 */
static GEN
zetamult_interpolate2_i(GEN avec, GEN t, long prec)
{
  pari_sp av;
  GEN a, b, ze, _1;
  long i, j, n, s, l;

  zparams(&s, &l, prec2nbits(prec));
  n = lg(avec) - 1;
  a = zeromatcopy(n + 1, l);
  b = zerovec(l + 1);
  for (i = 1; i <= n; i++)
  {
    long h = avec[n + 1 - i];
    if (i == 1) gel(b, h) = gen_1;
    else
      for (j = l + 1; j >= 1; j--)
      {
        if (j <= h) gel(b, j) = gen_0;
        else gel(b, j) = gadd(gcoeff(a, i, j-h), gmul(t, gel(b, j-h)));
      }
    gcoeff(a, i+1, 1) = gel(b, 2); /* j = 1 */
    for (j = 2; j <= l; j++)
    { /* b[j+1] - sum_{0 <= u < j-1} binom(j,u) a[i+1,u+1]*/
      pari_sp av = avma;
      GEN S = gel(b, j + 1);
      S = gsub(S, gcoeff(a, i+1, 1)); /* u = 0 */
      if (j > 2) S = gsub(S, gmulgu(gcoeff(a, i+1, 2), j)); /* u = 1 */
      if (j >= 4)
      {
        GEN C = utoipos(j*(j-1) / 2);
        long u, U = (j-1)/2;
        for (u = 2; u <= U; u++)
        { /* C = binom(j, u) = binom(j, j-u) */
          GEN A = gadd(gcoeff(a, i+1, u+1), gcoeff(a, i+1, j-u+1));
          S = gsub(S, gmul(C, A));
          C = diviuexact(muliu(C, j-u), u+1);
        }
        if (!odd(j)) S = gsub(S, gmul(C, gcoeff(a,i+1, j/2+1)));
      }
      gcoeff(a, i+1, j) = gerepileupto(av, gdivgu(S, j));
    }
  }
  _1 = real_1(prec + EXTRAPREC64); av = avma;
  ze = cgetg(n+1, t_VEC);
  for (i = 1; i <= n; i++)
  {
    GEN S = gdivgu(gcoeff(a, n+2-i, 1), s), sj = utoipos(s);
    for (j = 2; j <= l; j++)
    {
      sj = muliu(sj, s); /* = s^j */
      S = gadd(S, gdiv(gcoeff(a, n+2-i, j), sj));
    }
    gel(ze, i) = S;
  }
  for (i = s; i >= 1; i--)
  {
    GEN b0 = divri(_1, powuu(i, avec[n])), z;
    z = gel(ze,n); gel(ze,n) = gadd(z, b0);
    for (j = n-1; j >= 1; j--)
    {
      b0 = gdiv(gadd(gmul(t, b0), z), powuu(i, avec[j]));
      z = gel(ze,j); gel(ze,j) = gadd(gel(ze,j), b0);
    }
    if (gc_needed(av, 1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"zetamult: i = %ld", i);
      ze = gerepilecopy(av, ze);
    }
  }
  return gel(ze, 1);
}

/********************************************************************/
/**                      AKHILESH ALGORITHM                        **/
/********************************************************************/
/* a t_VECSMALL, upper bound for -log2(zeta(a)) */
static long
log2zeta_bound(GEN a)
{ return ceil(-dbllog2(zetamult_Zagier(a, 32, LOWDEFAULTPREC))); }
/* ibin[n+1] = 1 / binom(2n, n) as a t_REAL */
static void
get_ibin(GEN *pibin, GEN *pibin1, long N, long prec)
{
  GEN ibin, ibin1;
  long n;
  *pibin = ibin = cgetg(N + 2, t_VEC);
  *pibin1= ibin1= cgetg(N + 2, t_VEC);
  gel(ibin,1) = gel(ibin1,1) = gen_0; /* unused */
  gel(ibin,2) = gel(ibin1,2) = real2n(-1,prec);
  for (n = 2; n <= N; n++)
  {
    gel(ibin, n+1) = divru(mulru(gel(ibin, n), n), 4*n-2);
    gel(ibin1, n+1) = divru(gel(ibin, n+1), n);
  }
}
/**************************************************************/
/*                         ALL MZV's                          */
/**************************************************************/

/* Generalization to Multiple Polylogarithms.
The basic function is polylogmult which takes two arguments
avec, as above, and zvec, the vector of z_j, and the result
is $\sum_{n_1>n_2>...>n_r>0}z_1^{n_1}...z_r^{n_r}/(n_1^a_1...n_r^{a_r})$. */

/* Given admissible evec = xe_2....e_{k-1}y, (k>=2), compute a,b,v such that
evec = x{1}_{a-1}v{0}_{b-1}y with v empty or admissible.
Input: vector w=evec
Output: v=wmid, wini, wfin */
static void
findabvgen(GEN evec, long _0, long _1, GEN *pwmid, GEN *pwini, GEN *pwfin,
           long *pa, long *pb)
{
  long s = lg(evec) - 1, m, a, b, j, x = evec[1], y = evec[s];
  GEN wmid, wini, wfin;
  if (s == 2)
  {
    *pwmid = cgetg(1, t_VECSMALL);
    *pwini = mkvecsmall(x);
    *pwfin = mkvecsmall(y);
    *pa = *pb = 1; return;
  }
  a = s - 1;
  for (j = 1; j <= s - 2; j++) if (evec[j + 1] != _1) { a = j; break; }
  *pa = a;
  b = s - 1;
  for (j = s - 2; j >= 1; j--) if (evec[j + 1] != _0) { b = s - 1 - j; break; }
  *pb = b;

  *pwmid = wmid = a+b < s? vecslice(evec, a+1, s-b): cgetg(1, t_VECSMALL);
  m = lg(wmid) - 1;
  *pwini = wini = cgetg(a + m + 1, t_VECSMALL);
  wini[1] = x; for (j = 2; j <= a; j++) wini[j] = _1;
  for (; j <= a + m; j++) wini[j] = wmid[j-a];
  *pwfin = wfin = cgetg(b + m + 1, t_VECSMALL);
  for (j = 1; j <= m; j++) wfin[j] = wmid[j];
  for (; j < b + m; j++) wfin[j] = _0;
  wfin[j] = y;
}

/* y != 0,1 */
static GEN
filllg1(GEN ibin1, GEN r1, GEN y, long N, long prec)
{
  GEN v, y1 = gsubgs(gmulsg(2, y), 1), y3 = gmul(y, gsubsg(1, y));
  long n;

  v = cgetg(N + 2, t_VEC); gel(v, N + 1) = gen_0;
  if (gcmpgs(gnorm(y3),1) > 0)
  {
    GEN y2 = gdiv(r1, y3);
    for (n = N; n >= 1; n--)
    {
      pari_sp av2 = avma;
      GEN z = gmul(y2, gsub(gel(v, n+1), gmul(y1, gel(ibin1, n+1))));
      gel(v, n) = gerepileupto(av2, z);
    }
  }
  else
  {
    pari_sp av0 = avma;
    gel(v, 1) = gerepileupto(av0, glog(gdiv(y, gsubgs(y, 1)), prec));
    for (n = 1; n < N; n++)
    {
      pari_sp av2 = avma;
      GEN z = gadd(gmul(y3, gel(v, n)), gmul(y1, gel(ibin1, n+1)));
      gel(v, n + 1) = gerepileupto(av2, z);
    }
  }
  return v;
}
static GEN
fillrec(hashtable *H, GEN evec, long _0, long _1, GEN X, GEN pab, GEN r1,
        long N)
{
  long n, a, b, s, x0;
  GEN xy1, x, y, r, wmid, wini, wfin, mid, ini, fin;
  hashentry *ep = hash_search(H, evec);

  if (ep) return (GEN)ep->val;
  findabvgen(evec, _0, _1, &wmid, &wini, &wfin, &a, &b);
  x = gel(X, evec[1]); s = lg(evec)-1; /* > 1 */
  y = gel(X, evec[s]);
  mid = fillrec(H, wmid, _0, _1, X, pab, r1, N);
  ini = fillrec(H, wini, _0, _1, X, pab, r1, N);
  fin = fillrec(H, wfin, _0, _1, X, pab, r1, N);
  if (evec[1] == _0) { x0 = 1; xy1 = gdiv(r1, y); }
  else { x0 = 0; xy1 = gdiv(r1, gmul(gsubsg(1, x), y)); }
  r = cgetg(N+2, t_VEC); gel(r, N+1) = gen_0;
  for (n = N; n > 1; n--)
  {
    pari_sp av = avma;
    GEN t = gmul(gel(ini, n+1), gmael(pab, n, a));
    GEN u = gadd(gmul(gel(fin, n+1), gmael(pab, n, b)), gel(mid, n+1));
    GEN v = gdiv(x0? gadd(t, u): gsub(t, u), gmael(pab, n, a+b));
    gel(r, n) = gerepileupto(av, gmul(xy1, gadd(gel(r, n+1), v)));
  }
  { /* n = 1 */
    pari_sp av = avma;
    GEN t = gel(ini, 2), u = gadd(gel(fin, 2), gel(mid, 2));
    GEN v = x0? gadd(t, u): gsub(t, u);
    gel(r,1) = gerepileupto(av, gmul(xy1, gadd(gel(r,2), v)));
  }
  hash_insert(H, (void*)evec, (void*)r); return r;
}

static GEN
aztoe(GEN avec, GEN zvec, long prec)
{
  GEN y, E, u = subsr(1, real2n(10-prec2nbits(prec), LOWDEFAULTPREC));
  long i, l = lg(avec);

  E = cgetg(l, t_VEC); if (l == 1) return E;
  y = gen_1;
  for (i = 1; i < l; i++)
  {
    long a = avec[i];
    GEN e, zi = gel(zvec, i);
    if (a <= 0 || (a == 1 && i == 1 && gequal1(zi)))
      pari_err_TYPE("polylogmult [divergent]", avec);
    if (gequal0(zi)) return NULL;
    gel(E, i) = e = zerovec(a);
    gel(e, a) = y = gdiv(y, zi);
    if (gcmp(gnorm(y), u) < 0) pari_err_TYPE("polylogmult [divergent]", zvec);
  }
  return shallowconcat1(E);
}

/***********************************************************/
/* Special case of zvec = [1,1,...], i.e., zetamult.       */
/***********************************************************/
static void
findabvgens(GEN evec, GEN *pwmid, GEN *pwini, GEN *pwfin, long *pa, long *pb)
{
  GEN wmid, wini, wfin;
  long s = lg(evec) - 1, a, b, j, m;
  if (s == 2)
  {
    *pwmid = cgetg(1, t_VECSMALL);
    *pwini = mkvecsmall(0);
    *pwfin = mkvecsmall(1);
    *pa = *pb = 1; return;
  }
  a = s - 1;
  for (j = 1; j <= s - 2; j++) if (!evec[j + 1]) { a = j; break; }
  *pa = a;
  b = s - 1;
  for (j = s - 2; j >= 1; j--) if (evec[j + 1]) { b = s - 1 - j; break; }
  *pb = b;

  *pwmid = wmid = a+b < s? vecslice(evec, a+1, s-b): cgetg(1, t_VECSMALL);
  m = lg(wmid) - 1;
  *pwini = wini = cgetg(a + m + 1, t_VECSMALL);
  wini[1] = 0; for (j = 2; j <= a; j++) wini[j] = 1;
  for (; j <= a + m; j++) wini[j] = wmid[j-a];
  *pwfin = wfin = cgetg(b + m + 1, t_VECSMALL);
  for (j = 1; j <= m; j++) wfin[j] = wmid[j];
  for (; j < b + m; j++) wfin[j] = 0;
  wfin[j] = 1;
}
static GEN
fillrecs(hashtable *H, GEN evec, GEN pab, long N, long prec)
{
  long n, a, b;
  GEN r, wmid, wini, wfin, mid, ini, fin;
  hashentry *ep = hash_search(H, evec);

  if (ep) return (GEN)ep->val;
  findabvgens(evec, &wmid, &wini, &wfin, &a, &b);
  mid = fillrecs(H, wmid, pab, N, prec);
  ini = fillrecs(H, wini, pab, N, prec);
  fin = fillrecs(H, wfin, pab, N, prec);
  r = cgetg(N + 2, t_VEC); gel(r, N+1) = gen_0;
  for (n = N; n > 1; n--)
  {
    GEN z = cgetr(prec);
    pari_sp av = avma;
    GEN t = mpmul(gel(ini, n+1), gmael(pab, n, a));
    GEN u = mpadd(mpmul(gel(fin, n+1), gmael(pab, n, b)), gel(mid,n+1));
    GEN v = mpdiv(mpadd(t, u), gmael(pab, n, a+b));
    mpaff(mpadd(gel(r, n+1), v), z); set_avma(av); gel(r,n) = z;
  }
  { /* n = 1 */
    GEN z = cgetr(prec);
    pari_sp av = avma;
    GEN t = gel(ini,2), u = mpadd(gel(fin,2), gel(mid,2)), v = mpadd(t, u);
    mpaff(mpadd(gel(r, 2), v), z); set_avma(av); gel(r,1) = z;
  }
  hash_insert(H, (void*)evec, (void*)r); return r;
}
/* [n, ..., n^k] */
static GEN
powersu(ulong n, long k)
{
  GEN v, gn = utoipos(n);
  long i, l = k+1;
  v = cgetg(l, t_VEC); gel(v,1) = gn;
  for (i = 2; i < l; i++) gel(v,i) = muliu(gel(v,i-1), n);
  return v;
}
/* n^a = pab[n][a] */
static GEN
get_pab(long N, long k)
{
  GEN v = cgetg(N+1, t_VEC);
  long j;
  gel(v, 1) = gen_0; /* not needed */
  for (j = 2; j <= N; j++) gel(v, j) = powersu(j, k);
  return v;
}
static hashtable *
zetamult_hash(long _0, long _1, GEN ibin, GEN ibin1)
{
  hashtable *H = hash_create(4096, (ulong(*)(void*))&hash_zv,
                                   (int(*)(void*,void*))&zv_equal, 1);
  hash_insert(H, (void*)cgetg(1, t_VECSMALL), (void*)ibin);
  hash_insert(H, (void*)mkvecsmall(_0), (void*)ibin1);
  hash_insert(H, (void*)mkvecsmall(_1), (void*)ibin1); return H;
}
/* Akhilesh recursive algorithm, #a > 1;
 * e t_VECSMALL, prec final precision, bit required bitprecision */
static GEN
zetamult_Akhilesh(GEN e, long bit, long prec)
{
  long k = lg(e) - 1, N = 1 + bit/2, prec2 = nbits2prec(bit);
  GEN r, pab, ibin, ibin1;
  hashtable *H;

  get_ibin(&ibin, &ibin1, N, prec2);
  pab = get_pab(N, k);
  H = zetamult_hash(0, 1, ibin, ibin1);
  r = fillrecs(H, e, pab, lg(pab)-1, prec2);
  if (DEBUGLEVEL) err_printf("polylogmult: k = %ld, %ld nodes\n", k, H->nb);
  return gprec_wtrunc(gel(r,1), prec);
}

/* lump together close entries + round to 1 entries that are ~ 1 */
static GEN
vec_round(GEN V, long b)
{
  GEN v = shallowcopy(V), w = shallowcopy(v);
  long i, j, l = lg(v);
  for (i = 1; i < l; i++)
  {
    long e;
    if (!gel(v,i)) continue;
    if (gexpo(gsubgs(gel(v,i), 1)) < b) gel(w,i) = gel(v,i) = gen_1;
    e = gexpo(gel(v,i));
    for (j = i+1; j < l; j++)
      if (gel(v,j) && gexpo(gsub(gel(v,i), gel(v,j))) - e < b)
      {
        gel(v,j) = NULL;
        gel(w,j) = gel(w,i);
      }
  }
  return w;
}

/* evec t_VEC */
static GEN
zetamultevec(GEN evec, long prec)
{
  pari_sp av = avma;
  double *x, *y, z = 0;
  long b, i, j, l, lH, bitprec, prec2, N, _0 = 0, _1 = 0, k = lg(evec) - 1;
  GEN r1, r, pab, ibin, ibin1, X, Evec, v;
  hashtable *H;

  if (k == 0) return gen_1;
  evec = vec_round(evec, 3 - prec2nbits(prec));
  v = vec_equiv(evec); l = lg(v);
  Evec = cgetg(k+1, t_VECSMALL);
  X = cgetg(l + 2, t_VEC); /* our alphabet */
  for (i = lH = 1; i < l; i++)
  {
    GEN vi = gel(v,i), xi = gel(evec, vi[1]);
    long li = lg(vi);
    gel(X,i) = xi;
    if (!_0 && gequal0(xi)) _0 = i;
    else if (!_1 && gequal1(xi)) _1 = i;
    for (j = 1; j < li; j++) Evec[vi[j]] = i;
  }
  /* add 0,1 if needed */
  if (!_0) { gel(X, i) = gen_0; _0 = i++; }
  if (!_1) { gel(X, i) = gen_1; _1 = i++; }
  l = i; setlg(X, l);
  av = avma;
  x = (double*)stack_malloc_align(l * sizeof(double), sizeof(double));
  y = (double*)stack_malloc_align(l * sizeof(double), sizeof(double));
  for (j = 1; j < l; j++)
  {
    GEN t = gel(X,j);
    x[j] = (j == _1)? 0: -dbllog2(gsubsg(1, t));
    y[j] = (j == _0)? 0: -dbllog2(t);
  }
  for (i = 1; i < l; i++)
    for (j = i+1; j < l; j++) z = maxdd(z, x[i] + y[j]);
  b = 0;
  for (i = 1; i < l; i++)
  {
    GEN t = real_i(gel(X,i));
    long e = -gexpo(gsubgs(gmul2n(t,1), 1));
    b = maxss(b, e);
  }
  set_avma(av);
  if (z >= 2) pari_err_IMPL("polylogmult in this range");
  bitprec = prec2nbits(prec) + 64*(1 + (k >> 5));
  N = 1 + bitprec / (2 - z);
  bitprec += z * N;
  prec2 = nbits2prec(bitprec + b);
  X = gprec_wensure(X, prec2);
  get_ibin(&ibin, &ibin1, N, prec2);
  pab = get_pab(N, k);
  H = zetamult_hash(_0, _1, ibin, ibin1);
  r1 = real_1(prec2);
  for (i = 1; i < l; i++)
    if (i != _0 && i != _1)
      hash_insert(H, mkvecsmall(i), filllg1(ibin1, r1, gel(X,i), N, prec2));
  r = fillrec(H, Evec, _0, _1, X, pab, r1, N);
  if (DEBUGLEVEL) err_printf("polylogmult: k = %ld, %ld nodes\n", k, H->nb);
  return gprec_wtrunc(gel(r,1), prec);
}

/* a t_VECSMALL */
static GEN
zetamult_i(GEN a, long prec)
{
  long r = lg(a)-1, k, bit;
  if (r == 0) return gen_1;
  if (r == 1) return szeta(a[1], prec);
  bit = prec2nbits(prec);
  if (bit <= 128)
    return zetamult_Zagier(a, bit, prec + EXTRAPREC64);
  k = zv_sum(a);
  if (((double)r) / (k*k) * bit / log((double)10*bit) < 0.5)
    return zetamult_Zagier(a, bit, prec + EXTRAPREC64);
  bit += maxss(log2zeta_bound(a), 64);
  return zetamult_Akhilesh(atoe(a), bit, prec);
}
GEN
zetamult(GEN s, long prec)
{
  pari_sp av0 = avma, av;
  GEN z, avec, r = cgetr(prec);

  av = avma; avec = zetamultconvert_i(s,1);
  if (lg(avec) == 1) return gc_const(av0, gen_1);
  z = zetamult_i(avec, prec); affrr(z, r); return gc_const(av, r);
}

/* If star = NULL: MZV, otherwise Yamamoto interpolation (MZSV for t=1) */
GEN
zetamult_interpolate(GEN s, GEN t, long prec)
{
  pari_sp av = avma, av2;
  long i, k, l, la, bit;
  GEN avec, v, V;

  if (lg(s) == 1) return gen_1;
  if (!t) return zetamult(s, prec);
  avec = zetamultconvert_i(s, 1); k = zv_sum(avec);
  bit = prec2nbits(prec);
  if (bit <= 128 || k > 20 || (bit >> k) < 4)
    return zetamult_interpolate2_i(vecsmall_reverse(avec), t, prec);
  v = allstar(avec); l = lg(v); la = lg(avec);
  V = cgetg(la, t_VEC);
  for (i = 1; i < la; i++)
  { gel(V,i) = cgetr(prec + EXTRAPREC64); affur(0, gel(V,i)); }
  av2 = avma;
  for (i = 1; i < l; i++, set_avma(av2))
  {
    GEN a = gel(v,i); /* avec */
    long n = lg(a)-1; /* > 0 */
    affrr(addrr(gel(V,n), zetamult_i(a, prec)), gel(V,n));
  }
  return gerepileupto(av, poleval(vecreverse(V),t));
}

GEN
polylogmult_interpolate(GEN a, GEN z, GEN t, long prec)
{
  pari_sp av = avma;
  GEN V, avec, A, AZ, Z;
  long i, la, l;
  switch(typ(a))
  {
    case t_VEC:
    case t_COL: a = gtovecsmall(a); break;
    case t_VECSMALL: break;
    default: pari_err_TYPE("polylogmult", a);
             return NULL;/*LCOV_EXCL_LINE*/
  }
  if (!z) return zetamult_interpolate(a, t, prec);
  switch (typ(z))
  {
    case t_INT: case t_FRAC: case t_REAL: case t_COMPLEX:
      z = mkvec(z); break;
    case t_VEC: case t_COL: break;
    case t_VECSMALL: z = zv_to_ZV(z); break;
    default: pari_err_TYPE("polylogmult [z]", z);
  }
  if (lg(z) != lg(a))
    pari_err_TYPE("polylogmult [#s != #z]", mkvec2(a,z));
  if (!t) return gerepilecopy(av, zetamultevec(aztoe(a,z,prec), prec));
  avec = zetamultconvert_i(a, 1); la = lg(avec);
  AZ = allstar2(avec, z);
  A = gel(AZ, 1); l = lg(A);
  Z = gel(AZ, 2); V = zerovec(la-1);
  for (i = 1; i < l; i++)
  {
    pari_sp av2 = avma;
    GEN ai = gel(A,i), e = aztoe(ai, gel(Z,i), prec);
    long n = lg(ai)-1; /* > 0 */
    gel(V,n) = gerepileupto(av2, gadd(gel(V,n), zetamultevec(e, prec)));
  }
  return gerepileupto(av, poleval(vecreverse(V),t));
}

GEN
polylogmult(GEN a, GEN z, long prec)
{
  return polylogmult_interpolate(a, z, NULL, prec);
}

/**************************************************************/
/*                           ALL MZV's                        */
/**************************************************************/

/* Given admissible evec w = 0e_2....e_{k-1}1, compute a,b,v such that
 * w=0{1}_{b-1}v{0}_{a-1}1 with v empty or admissible.
 * Input: binary vector evec */
static void
findabv(GEN w, long *pa, long *pb, long *pminit, long *pmmid, long *pmfin)
{
  long le = lg(w) - 2;
  if (le == 0)
  {
    *pa = 1;
    *pb = 1;
    *pminit = 2;
    *pmfin = 2;
    *pmmid = 1;
  }
  else
  {
    long a, b, j, lv;
    for (j = 1; j <= le; j++)
      if (!w[j+1]) break;
    *pb = b = j;
    for (j = le; j >= 1; j--)
      if (w[j+1]) break;
    *pa = a = le + 1 - j;
    lv = le + 2 - a - b;
    if (lv > 0)
    {
      long v = fd(w, b + 1, le - a + 2), u = v + (1 << (lv-1));
      *pminit = (((1 << b) - 1) << (lv - 1)) + (v/2) + 2;
      *pmfin = (u << (a - 1)) + 2;
      *pmmid = (u >> 1) + 2;
    }
    else
    {
      *pminit = (1 << (b - 1)) + 1;
      *pmfin = (a == 1) ? 2 : (1 << (a - 2)) + 2;
      *pmmid = 1;
    }
  }
}

/* Returns L:
* L[1] contains zeta(emptyset)_{n-1,n-1},
* L[2] contains zeta({0})_{n-1,n-1}=zeta({1})_{n-1,n-1} for n >= 2,
* L[m+2][n] : 1 <= m < 2^{k-2}, 1 <= n <= N + 1
* contains zeta(w)_{n-1,n-1}, w corresponding to m,n
* L[m+2] : 2^{k-2} <= m < 2^{k-1} contains zeta(w), w corresponding to m
(code: w=0y1 iff m=1y). */
static GEN
fillL(long k, long bitprec)
{
  long N = 1 + bitprec/2, prec = nbits2prec(bitprec);
  long s, j, n, m, K = 1 << (k - 1), K2 = K/2;
  GEN p1, p2, pab = get_pab(N, k), L = cgetg(K + 2, t_VEC);

  get_ibin(&gel(L,1), &gel(L,2), N, prec);
  for (m = 1; m < K2; m++)
  {
    gel(L, m+2) = p1 = cgetg(N+1, t_VEC);
    for (n = 1; n < N; n++) gel(p1, n) = cgetr(prec);
    gel(p1, n) = gen_0;
  }
  for (m = K2; m < K; m++) gel(L, m+2) = utor(0, prec);
  for (s = 2; s <= k; s++)
  { /* Assume length evec < s filled */
    /* If evec = 0e_2...e_{s-1}1 then m = (1e_2...e_{s-1})_2 */
    GEN w = cgetg(s, t_VECSMALL);
    long M = 1 << (s - 2);
    pari_sp av = avma;
    for (m = M; m < 2*M; m++)
    {
      GEN pinit, pfin, pmid;
      long comp, a, b, mbar, minit, mfin, mmid, mc;
      p1 = gel(L, m + 2);
      for (j = s - 1, mc = m, mbar = 1; j >= 2; j--, mc >>= 1)
      {
        w[j] = mc & 1;
        mbar = (1 - w[j]) | (mbar << 1);
      }
      /* m, mbar are dual; handle smallest, copy the other */
      comp = mbar - m; if (comp < 0) continue; /* m > mbar */
      if (comp)
      {
        p2 = gel(L, mbar + 2);
        setisclone(p2); /* flag as dual */
      }
      else
        p2 = NULL; /* no copy needed if m = mbar */
      findabv(w, &a,&b,&minit,&mmid,&mfin);
      pinit= gel(L, minit);
      pfin = gel(L, mfin);
      pmid = gel(L, mmid);
      for (n = N-1; n > 1; n--, set_avma(av))
      {
        GEN t = mpmul(gel(pinit,n+1), gmael(pab, n, b));
        GEN u = mpmul(gel(pfin, n+1), gmael(pab, n, a));
        GEN v = gel(pmid, n+1), S = s < k ? gel(p1, n+1): p1;
        S = mpadd(S, mpdiv(mpadd(mpadd(t, u), v), gmael(pab, n, a+b)));
        mpaff(S, s < k ? gel(p1, n) : p1);
        if (p2 && s < k) mpaff(S, gel(p2, n));
      }
      { /* n = 1: same formula simplifies */
        GEN t = gel(pinit,2), u = gel(pfin,2), v = gel(pmid,2);
        GEN S = s < k ? gel(p1,2): p1;
        S = mpadd(S, mpadd(mpadd(t, u), v));
        mpaff(S, s < k ? gel(p1,1) : p1);
        if (p2 && s < k) mpaff(S, gel(p2, 1));
        set_avma(av);
      }
      if (p2 && s == k) mpaff(p1, p2);
    }
  }
  return L;
}

/* bit 1 of flag unset: full, otherwise up to duality (~ half)
 * bit 2 of flag unset: all <= k, otherwise only k
 * half: 2^(k-3)+ delta_{k even} * 2^(k/2-2), sum = 2^(k-2)+2^(floor(k/2)-1)-1
 * full: 2^(k-2); sum = 2^(k-1)-1 */
static GEN
zetamultall_i(long k, long flag, long prec)
{
  GEN res, ind, L = fillL(k, prec2nbits(prec) + 32);
  long m, K2 = 1 << (k-2), n = lg(L) - 1, m0 = (flag & 4L) ? K2 : 1;

  if (!(flag & 2L))
  {
    res = cgetg(n - m0, t_VEC);
    ind = cgetg(n - m0, t_VECSMALL);
    for (m = m0; m < n - 1; m++)
    {
      gel(res, m - m0 + 1) = m < K2 ? gmael(L, m + 2, 1) : gel(L, m + 2);
      ind[m - m0 + 1] = m;
    }
  }
  else
  { /* up to duality */
    long nres, c;
    if (k == 2) nres = 1;
    else if (!(flag & 2L))
      nres = (1 << (k - 2)) + (1 << ((k/2) - 1)) - 1;
    else
      nres = (1 << (k - 1));
    res = cgetg(nres + 1, t_VEC);
    ind = cgetg(nres + 1, t_VECSMALL);
    for (m = m0, c = 1; m < n - 1; m++)
    {
      GEN z = gel(L,m+2);
      if (isclone(z)) continue; /* dual */
      if (m < K2) z = gel(z,1);
      gel(res, c) = z;
      ind[c] = m; c++;
    }
    setlg(res, c);
    setlg(ind, c);
  }
  return mkvec2(res, ind);
}

/* fd(e, 2, lg(e)-2), e = atoe(avec) */
static long
atom(GEN avec)
{
  long i, m, l = lg(avec);
  if (l < 3) return 0;
  m = 1; /* avec[1] != 0 */
  for (i = 2; i < l-1; i++) m = (m << avec[i]) + 1;
  return m << (avec[i]-1);
}
static long
atoind(GEN avec, long flag)
{ return atom(avec) + (flag? 1: (1 << (zv_sum(avec) - 2))); }
/* If flag is unset, L has all k1 <= k, otherwise only k */
static GEN
zetamultstar_i(GEN L, GEN avec, long flag)
{
  GEN s = allstar(avec), S = gen_0;
  long i, l = lg(s);
  for (i = 1; i < l; i++) S = gadd(S, gel(L, atoind(gel(s,i), flag)));
  return S;
}

/* bit 0: notstar/star
 * bit 1: full/half (ignored if star, always full)
 * bit 2: all <= k / only k
 * bit 3: without / with index */
GEN
zetamultall(long k, long flag, long prec)
{
  pari_sp av = avma;
  GEN Lind, L, res;
  long K, k1, ct, fl;

  if (flag < 0 || flag > 15) pari_err_FLAG("zetamultall");
  if (k < 1) pari_err_DOMAIN("zetamultall", "k", "<", gen_1, stoi(k));
  if (k == 1) return cgetg(1, t_VEC);
  if (k >= 64) pari_err_OVERFLOW("zetamultall");
  if (!(flag & 1L))
  { /* not star */
    Lind = zetamultall_i(k, flag, prec);
    res = (flag & 8L)? Lind : gel(Lind, 1);
    return gerepilecopy(av, res);
  }
  /* star */
  fl = flag & 4L; /* 4 if k, else 0 (all <= k) */
  Lind = gerepilecopy(av, zetamultall_i(k, fl, prec)); /* full */
  L = gel(Lind, 1);
  K = 1 << (k - 2);
  res = cgetg(fl? K+1: 2*K, t_VEC);
  for (ct = 1, k1 = fl? k: 2; k1 <= k; k1++)
  {
    GEN w = cgetg(k1 + 1, t_VECSMALL);
    long M = 1 << (k1 - 1), m;
    for (m = 1; m <= M; m += 2)
    {
      pari_sp av = avma;
      long j, mc = m;
      for (j = k1; j >= 1; j--) { w[j] = mc & 1; mc >>= 1; }
      gel(res, ct++) = gerepileupto(av, zetamultstar_i(L, etoa(w), fl));
    }
  }
  if (flag & 8L) res = mkvec2(res, gel(Lind, 2));
  return gerepilecopy(av, res);
}
