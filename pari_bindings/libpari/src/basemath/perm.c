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

/*************************************************************************/
/**                                                                     **/
/**                   Routines for handling VEC/COL                     **/
/**                                                                     **/
/*************************************************************************/
int
vec_isconst(GEN v)
{
  long i, l = lg(v);
  GEN w;
  if (l==1) return 1;
  w = gel(v,1);
  for(i=2; i<l; i++)
    if (!gequal(gel(v,i), w)) return 0;
  return 1;
}

int
vecsmall_isconst(GEN v)
{
  long i, l = lg(v);
  ulong w;
  if (l==1) return 1;
  w = uel(v,1);
  for(i=2; i<l; i++)
    if (uel(v,i) != w) return 0;
  return 1;
}

/* Check if all the elements of v are different.
 * Use a quadratic algorithm. Could be done in n*log(n) by sorting. */
int
vec_is1to1(GEN v)
{
  long i, j, l = lg(v);
  for (i=1; i<l; i++)
  {
    GEN w = gel(v,i);
    for(j=i+1; j<l; j++)
      if (gequal(gel(v,j), w)) return 0;
  }
  return 1;
}

GEN
vec_insert(GEN v, long n, GEN x)
{
  long i, l=lg(v);
  GEN V = cgetg(l+1,t_VEC);
  for(i=1; i<n; i++) gel(V,i) = gel(v,i);
  gel(V,n) = x;
  for(i=n+1; i<=l; i++) gel(V,i) = gel(v,i-1);
  return V;
}
/*************************************************************************/
/**                                                                     **/
/**                   Routines for handling VECSMALL                    **/
/**                                                                     **/
/*************************************************************************/
/* Sort v[0]...v[n-1] and put result in w[0]...w[n-1].
 * We accept v==w. w must be allocated. */
static void
vecsmall_sortspec(GEN v, long n, GEN w)
{
  pari_sp ltop=avma;
  long nx=n>>1, ny=n-nx;
  long m, ix, iy;
  GEN x, y;
  if (n<=2)
  {
    if (n==1)
      w[0]=v[0];
    else if (n==2)
    {
      long v0=v[0], v1=v[1];
      if (v0<=v1) { w[0]=v0; w[1]=v1; }
      else        { w[0]=v1; w[1]=v0; }
    }
    return;
  }
  x=new_chunk(nx); y=new_chunk(ny);
  vecsmall_sortspec(v,nx,x);
  vecsmall_sortspec(v+nx,ny,y);
  for (m=0, ix=0, iy=0; ix<nx && iy<ny; )
    if (x[ix]<=y[iy])
      w[m++]=x[ix++];
    else
      w[m++]=y[iy++];
  for(;ix<nx;) w[m++]=x[ix++];
  for(;iy<ny;) w[m++]=y[iy++];
  set_avma(ltop);
}

static long
vecsmall_sort_max(GEN v)
{
  long i, l = lg(v), max = -1;
  for (i = 1; i < l; i++)
    if (v[i] > max) { max = v[i]; if (max >= l) return -1; }
    else if (v[i] < 0) return -1;
  return max;
}
/* assume 0 <= v[i] <= M. In place. */
void
vecsmall_counting_sort(GEN v, long M)
{
  pari_sp av;
  long i, j, k, l;
  GEN T;
  if (M == 0) return;
  av = avma; T = new_chunk(M + 1); l = lg(v);
  for (i = 0; i <= M; i++) T[i] = 0;
  for (i = 1; i < l; i++) T[v[i]]++; /* T[j] is # keys = j */
  for (j = 0, k = 1; j <= M; j++)
    for (i = 1; i <= T[j]; i++) v[k++] = j;
  set_avma(av);
}
/* not GC-clean, suitable for gerepileupto */
GEN
vecsmall_counting_uniq(GEN v, long M)
{
  long i, k, l = lg(v);
  GEN T, U;
  if (l == 1) return cgetg(1, t_VECSMALL);
  if (M == 0) return mkvecsmall(0);
  if (l == 2) return leafcopy(v);
  U = new_chunk(M + 2);
  T = U+1; /* allows to rewrite result over T also if T[0] = 1 */
  for (i = 0; i <= M; i++) T[i] = 0;
  for (i = 1; i < l; i++) T[v[i]] = 1;
  for (i = 0, k = 1; i <= M; i++)
    if (T[i]) U[k++] = i;
  U[0] = evaltyp(t_VECSMALL) | _evallg(k); return U;
}
GEN
vecsmall_counting_indexsort(GEN v, long M)
{
  pari_sp av;
  long i, l = lg(v);
  GEN T, p;
  if (M == 0 || l <= 2) return identity_zv(l - 1);
  p = cgetg(l, t_VECSMALL); av = avma; T = new_chunk(M + 1);
  for (i = 0; i <= M; i++) T[i] = 0;
  for (i = 1; i < l; i++) T[v[i]]++; /* T[j] is # keys = j */
  for (i = 1; i <= M; i++) T[i] += T[i-1]; /* T[j] is # keys <= j */
  for (i = l-1; i > 0; i--) { p[T[v[i]]] = i; T[v[i]]--; }
  return gc_const(av, p);
}

/* in place sort */
void
vecsmall_sort(GEN v)
{
  long n = lg(v) - 1, max;
  if (n <= 1) return;
  if ((max = vecsmall_sort_max(v)) >= 0)
    vecsmall_counting_sort(v, max);
  else
    vecsmall_sortspec(v+1, n, v+1);
}

/* cf gen_sortspec */
static GEN
vecsmall_indexsortspec(GEN v, long n)
{
  long nx, ny, m, ix, iy;
  GEN x, y, w;
  switch(n)
  {
    case 1: return mkvecsmall(1);
    case 2: return (v[1] <= v[2])? mkvecsmall2(1,2): mkvecsmall2(2,1);
    case 3:
      if (v[1] <= v[2]) {
        if (v[2] <= v[3]) return mkvecsmall3(1,2,3);
        return (v[1] <= v[3])? mkvecsmall3(1,3,2)
                             : mkvecsmall3(3,1,2);
      } else {
        if (v[1] <= v[3]) return mkvecsmall3(2,1,3);
        return (v[2] <= v[3])? mkvecsmall3(2,3,1)
                             : mkvecsmall3(3,2,1);
      }
  }
  nx = n>>1; ny = n-nx;
  w = cgetg(n+1,t_VECSMALL);
  x = vecsmall_indexsortspec(v,nx);
  y = vecsmall_indexsortspec(v+nx,ny);
  for (m=1, ix=1, iy=1; ix<=nx && iy<=ny; )
    if (v[x[ix]] <= v[y[iy]+nx])
      w[m++] = x[ix++];
    else
      w[m++] = y[iy++]+nx;
  for(;ix<=nx;) w[m++] = x[ix++];
  for(;iy<=ny;) w[m++] = y[iy++]+nx;
  set_avma((pari_sp)w); return w;
}

/*indirect sort.*/
GEN
vecsmall_indexsort(GEN v)
{
  long n = lg(v) - 1, max;
  if (n==0) return cgetg(1, t_VECSMALL);
  if ((max = vecsmall_sort_max(v)) >= 0)
    return vecsmall_counting_indexsort(v, max);
  else
    return vecsmall_indexsortspec(v,n);
}

/* assume V sorted */
GEN
vecsmall_uniq_sorted(GEN v)
{
  long i, j, l;
  GEN w = cgetg_copy(v, &l);
  if (l == 1) return w;
  w[1] = v[1];
  for(i = j = 2; i < l; i++)
    if (v[i] != w[j-1]) w[j++] = v[i];
  stackdummy((pari_sp)(w + l), (pari_sp)(w + j));
  setlg(w, j); return w;
}

GEN
vecsmall_uniq(GEN v)
{
  pari_sp av = avma;
  long max;
  if ((max = vecsmall_sort_max(v)) >= 0)
    v = vecsmall_counting_uniq(v, max);
  else
  { v = zv_copy(v); vecsmall_sort(v); v = vecsmall_uniq_sorted(v); }
  return gerepileuptoleaf(av, v);
}

/* assume x sorted */
long
vecsmall_duplicate_sorted(GEN x)
{
  long i,k,l=lg(x);
  if (l==1) return 0;
  for (k=x[1],i=2; i<l; k=x[i++])
    if (x[i] == k) return i;
  return 0;
}

long
vecsmall_duplicate(GEN x)
{
  pari_sp av=avma;
  GEN p=vecsmall_indexsort(x);
  long k,i,r=0,l=lg(x);
  if (l==1) return 0;
  for (k=x[p[1]],i=2; i<l; k=x[p[i++]])
    if (x[p[i]] == k) { r=p[i]; break; }
  set_avma(av);
  return r;
}

static int
vecsmall_is1to1spec(GEN v, long n, GEN w)
{
  pari_sp ltop=avma;
  long nx=n>>1, ny=n-nx;
  long m, ix, iy;
  GEN x, y;
  if (n<=2)
  {
    if (n==1)
      w[0]=v[0];
    else if (n==2)
    {
      long v0=v[0], v1=v[1];
      if (v0==v1) return 0;
      else if (v0<v1) { w[0]=v0; w[1]=v1; }
      else            { w[0]=v1; w[1]=v0; }
    }
    return 1;
  }
  x = new_chunk(nx);
  if (!vecsmall_is1to1spec(v,nx,x))    return 0;
  y = new_chunk(ny);
  if (!vecsmall_is1to1spec(v+nx,ny,y)) return 0;
  for (m=0, ix=0, iy=0; ix<nx && iy<ny; )
    if (x[ix]==y[iy]) return 0;
    else if (x[ix]<y[iy])
      w[m++]=x[ix++];
    else
      w[m++]=y[iy++];
  for(;ix<nx;) w[m++]=x[ix++];
  for(;iy<ny;) w[m++]=y[iy++];
  set_avma(ltop);
  return 1;
}

int
vecsmall_is1to1(GEN V)
{
  pari_sp av = avma;
  long l;
  GEN W = cgetg_copy(V, &l);
  if (l <= 2) return 1;
  return gc_bool(av, vecsmall_is1to1spec(V+1,l,W+1));
}

/*************************************************************************/
/**                                                                     **/
/**             Routines for handling vectors of VECSMALL               **/
/**                                                                     **/
/*************************************************************************/

GEN
vecvecsmall_sort(GEN x)
{ return gen_sort(x, (void*)&vecsmall_lexcmp, cmp_nodata); }
GEN
vecvecsmall_sort_shallow(GEN x)
{ return gen_sort_shallow(x, (void*)&vecsmall_lexcmp, cmp_nodata); }

void
vecvecsmall_sort_inplace(GEN x, GEN *perm)
{ gen_sort_inplace(x, (void*)&vecsmall_lexcmp, cmp_nodata, perm); }

GEN
vecvecsmall_sort_uniq(GEN x)
{ return gen_sort_uniq(x, (void*)&vecsmall_lexcmp, cmp_nodata); }

GEN
vecvecsmall_indexsort(GEN x)
{ return gen_indexsort(x, (void*)&vecsmall_lexcmp, cmp_nodata); }

long
vecvecsmall_search(GEN x, GEN y)
{ return gen_search(x,y,(void*)vecsmall_prefixcmp, cmp_nodata); }

/* assume x non empty */
long
vecvecsmall_max(GEN x)
{
  long i, l = lg(x), m = vecsmall_max(gel(x,1));
  for (i = 2; i < l; i++)
  {
    long t = vecsmall_max(gel(x,i));
    if (t > m) m = t;
  }
  return m;
}

/*************************************************************************/
/**                                                                     **/
/**                  Routines for handling permutations                 **/
/**                                                                     **/
/*************************************************************************/

/* Permutations may be given by
 * perm (VECSMALL): a bijection from 1...n to 1...n i-->perm[i]
 * cyc (VEC of VECSMALL): a product of disjoint cycles. */

/* Multiply (compose) two permutations, putting the result in the second one. */
static void
perm_mul_inplace2(GEN s, GEN t)
{
  long i, l = lg(s);
  for (i = 1; i < l; i++) t[i] = s[t[i]];
}

GEN
vecperm_extendschreier(GEN C, GEN v, long n)
{
  pari_sp av = avma;
  long mj, lv = lg(v), m = 1, mtested = 1;
  GEN bit = const_vecsmall(n, 0);
  GEN cy = cgetg(n+1, t_VECSMALL);
  GEN sh = const_vec(n, gen_0);
  for(mj=1; mj<=n; mj++)
  {
    if (isintzero(gel(C,mj))) continue;
    gel(sh,mj) = gcopy(gel(C,mj));
    if (bit[mj]) continue;
    cy[m++] = mj;
    bit[mj] = 1;
    for(;;)
    {
      long o, mold = m;
      for (o = 1; o < lv; o++)
      {
        GEN vo = gel(v,o);
        long p;
        for (p = mtested; p < mold; p++) /* m increases! */
        {
          long j = vo[ cy[p] ];
          if (!bit[j])
          {
            gel(sh,j) = perm_mul(vo, gel(sh, cy[p]));
            cy[m++] = j;
          }
          bit[j] = 1;
        }
      }
      mtested = mold;
      if (m == mold) break;
    }
  }
  return gerepileupto(av, sh);
}

/* Orbits of the subgroup generated by v on {1,..,n} */
static GEN
vecperm_orbits_i(GEN v, long n)
{
  long mj = 1, lv = lg(v), k, l;
  GEN cycle = cgetg(n+1, t_VEC), bit = const_vecsmall(n, 0);
  for (k = 1, l = 1; k <= n;)
  {
    pari_sp ltop = avma;
    long m = 1;
    GEN cy = cgetg(n+1, t_VECSMALL);
    for (  ; bit[mj]; mj++) /*empty*/;
    k++; cy[m++] = mj;
    bit[mj++] = 1;
    for(;;)
    {
      long o, mold = m;
      for (o = 1; o < lv; o++)
      {
        GEN vo = gel(v,o);
        long p;
        for (p = 1; p < m; p++) /* m increases! */
        {
          long j = vo[ cy[p] ];
          if (!bit[j]) cy[m++] = j;
          bit[j] = 1;
        }
      }
      if (m == mold) break;
      k += m - mold;
    }
    setlg(cy, m);
    gel(cycle,l++) = gerepileuptoleaf(ltop, cy);
  }
  setlg(cycle, l); return cycle;
}
/* memory clean version */
GEN
vecperm_orbits(GEN v, long n)
{
  pari_sp av = avma;
  return gerepilecopy(av, vecperm_orbits_i(v, n));
}

static int
isperm(GEN v)
{
  pari_sp av = avma;
  long i, n = lg(v)-1;
  GEN w;
  if (typ(v) != t_VECSMALL) return 0;
  w = zero_zv(n);
  for (i=1; i<=n; i++)
  {
    long d = v[i];
    if (d < 1 || d > n || w[d]) return gc_bool(av,0);
    w[d] = 1;
  }
  return gc_bool(av,1);
}

/* Compute the cyclic decomposition of a permutation */
GEN
perm_cycles(GEN v)
{
  pari_sp av = avma;
  return gerepilecopy(av, vecperm_orbits_i(mkvec(v), lg(v)-1));
}

GEN
permcycles(GEN v)
{
  if (!isperm(v)) pari_err_TYPE("permcycles",v);
  return perm_cycles(v);
}

/* Output the order of p */
ulong
perm_orderu(GEN v)
{
  pari_sp av = avma;
  GEN c = vecperm_orbits_i(mkvec(v), lg(v)-1);
  long i, d;
  for(i=1, d=1; i<lg(c); i++) d = ulcm(d, lg(gel(c,i))-1);
  return gc_ulong(av,d);
}

static GEN
_domul(void *data, GEN x, GEN y)
{
  GEN (*mul)(GEN,GEN)=(GEN (*)(GEN,GEN)) data;
  return mul(x,y);
}

/* Output the order of p */
GEN
perm_order(GEN v)
{
  pari_sp av = avma;
  GEN c = vecperm_orbits_i(mkvec(v), lg(v)-1);
  long i, l = lg(c);
  GEN V = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(V,i) = utoi(lg(gel(c,i))-1);
  return gerepileuptoint(av, gen_product(V, (void *)lcmii, _domul));
}

GEN
permorder(GEN v)
{
  if (!isperm(v)) pari_err_TYPE("permorder",v);
  return perm_order(v);
}

/* sign of a permutation */
long
perm_sign(GEN v)
{
  pari_sp av = avma;
  GEN c = vecperm_orbits_i(mkvec(v), lg(v)-1);
  long i, l = lg(c), s = 1;
  for (i = 1; i < l; i++)
    if (odd(lg(gel(c, i)))) s = -s;
  return gc_long(av,s);
}

long
permsign(GEN v)
{
  if (!isperm(v)) pari_err_TYPE("permsign",v);
  return perm_sign(v);
}

GEN
Z_to_perm(long n, GEN x)
{
  pari_sp av;
  ulong i, r;
  GEN v = cgetg(n+1, t_VECSMALL);
  if (n==0) return v;
  uel(v,n) = 1; av = avma;
  if (signe(x) <= 0) x = modii(x, mpfact(n));
  for (r=n-1; r>=1; r--)
  {
    ulong a;
    x = absdiviu_rem(x, n+1-r, &a);
    for (i=r+1; i<=(ulong)n; i++)
      if (uel(v,i) > a) uel(v,i)++;
    uel(v,r) = a+1;
  }
  set_avma(av); return v;
}
GEN
numtoperm(long n, GEN x)
{
  if (n < 0) pari_err_DOMAIN("numtoperm", "n", "<", gen_0, stoi(n));
  if (typ(x) != t_INT) pari_err_TYPE("numtoperm",x);
  return Z_to_perm(n, x);
}

/* destroys v */
static GEN
perm_to_Z_inplace(GEN v)
{
  long l = lg(v), i, r;
  GEN x = gen_0;
  if (!isperm(v)) return NULL;
  for (i = 1; i < l; i++)
  {
    long vi = v[i];
    if (vi <= 0) return NULL;
    x = i==1 ? utoi(vi-1): addiu(muliu(x,l-i), vi-1);
    for (r = i+1; r < l; r++)
      if (v[r] > vi) v[r]--;
  }
  return x;
}
GEN
perm_to_Z(GEN v)
{
  pari_sp av = avma;
  GEN x = perm_to_Z_inplace(leafcopy(v));
  if (!x) pari_err_TYPE("permtonum",v);
  return gerepileuptoint(av, x);
}
GEN
permtonum(GEN p)
{
  pari_sp av = avma;
  GEN v, x;
  switch(typ(p))
  {
    case t_VECSMALL: return perm_to_Z(p);
    case t_VEC: case t_COL:
      if (RgV_is_ZV(p)) { v = ZV_to_zv(p); break; }
    default: pari_err_TYPE("permtonum",p);
      return NULL;/*LCOV_EXCL_LINE*/
  }
  x = perm_to_Z_inplace(v);
  if (!x) pari_err_TYPE("permtonum",p);
  return gerepileuptoint(av, x);
}

GEN
cyc_pow(GEN cyc, long exp)
{
  long i, j, k, l, r;
  GEN c;
  for (r = j = 1; j < lg(cyc); j++)
  {
    long n = lg(gel(cyc,j)) - 1;
    r += cgcd(n, exp);
  }
  c = cgetg(r, t_VEC);
  for (r = j = 1; j < lg(cyc); j++)
  {
    GEN v = gel(cyc,j);
    long n = lg(v) - 1, e = umodsu(exp,n), g = (long)ugcd(n, e), m = n / g;
    for (i = 0; i < g; i++)
    {
      GEN p = cgetg(m+1, t_VECSMALL);
      gel(c,r++) = p;
      for (k = 1, l = i; k <= m; k++)
      {
        p[k] = v[l+1];
        l += e; if (l >= n) l -= n;
      }
    }
  }
  return c;
}

/* Compute the power of a permutation given by product of cycles
 * Ouput a perm, not a cyc */
GEN
cyc_pow_perm(GEN cyc, long exp)
{
  long e, j, k, l, n;
  GEN p;
  for (n = 0, j = 1; j < lg(cyc); j++) n += lg(gel(cyc,j))-1;
  p = cgetg(n + 1, t_VECSMALL);
  for (j = 1; j < lg(cyc); j++)
  {
    GEN v = gel(cyc,j);
    n = lg(v) - 1; e = umodsu(exp, n);
    for (k = 1, l = e; k <= n; k++)
    {
      p[v[k]] = v[l+1];
      if (++l == n) l = 0;
    }
  }
  return p;
}

GEN
perm_pow(GEN perm, GEN exp)
{
  long i, r = lg(perm)-1;
  GEN p = zero_zv(r);
  pari_sp av = avma;
  GEN v = cgetg(r+1, t_VECSMALL);
  for (i=1; i<=r; i++)
  {
    long e, n, k, l;
    if (p[i]) continue;
    v[1] = i;
    for (n=1, k=perm[i]; k!=i; k=perm[k], n++) v[n+1] = k;
    e = umodiu(exp, n);
    for (k = 1, l = e; k <= n; k++)
    {
      p[v[k]] = v[l+1];
      if (++l == n) l = 0;
    }
  }
  set_avma(av); return p;
}

GEN
perm_powu(GEN perm, ulong exp)
{
  ulong i, r = lg(perm)-1;
  GEN p = zero_zv(r);
  pari_sp av = avma;
  GEN v = cgetg(r+1, t_VECSMALL);
  for (i=1; i<=r; i++)
  {
    ulong e, n, k, l;
    if (p[i]) continue;
    v[1] = i;
    for (n=1, k=perm[i]; k!=i; k=perm[k], n++) v[n+1] = k;
    e = exp % n;
    for (k = 1, l = e; k <= n; k++)
    {
      p[v[k]] = v[l+1];
      if (++l == n) l = 0;
    }
  }
  set_avma(av); return p;
}

GEN
perm_to_GAP(GEN p)
{
  pari_sp ltop=avma;
  GEN gap;
  GEN x;
  long i;
  long nb, c=0;
  char *s;
  long sz;
  long lp=lg(p)-1;
  if (typ(p) != t_VECSMALL)  pari_err_TYPE("perm_to_GAP",p);
  x = perm_cycles(p);
  sz = (long) ((bfffo(lp)+1) * LOG10_2 + 1);
  /*Dry run*/
  for (i = 1, nb = 1; i < lg(x); ++i)
  {
    GEN z = gel(x,i);
    long lz = lg(z)-1;
    nb += 1+lz*(sz+2);
  }
  nb++;
  /*Real run*/
  gap = cgetg(nchar2nlong(nb) + 1, t_STR);
  s = GSTR(gap);
  for (i = 1; i < lg(x); ++i)
  {
    long j;
    GEN z = gel(x,i);
    if (lg(z) > 2)
    {
      s[c++] = '(';
      for (j = 1; j < lg(z); ++j)
      {
        if (j > 1)
        {
          s[c++] = ','; s[c++] = ' ';
        }
        sprintf(s+c,"%ld",z[j]);
        while(s[c++]) /* empty */;
        c--;
      }
      s[c++] = ')';
    }
  }
  if (!c) { s[c++]='('; s[c++]=')'; }
  s[c] = '\0';
  return gerepileupto(ltop,gap);
}

int
perm_commute(GEN s, GEN t)
{
  long i, l = lg(t);
  for (i = 1; i < l; i++)
    if (t[ s[i] ] != s[ t[i] ]) return 0;
  return 1;
}

/*************************************************************************/
/**                                                                     **/
/**                  Routines for handling groups                       **/
/**                                                                     **/
/*************************************************************************/
/* A Group is a t_VEC [gen,orders]
 * gen (vecvecsmall): list of generators given by permutations
 * orders (vecsmall): relatives orders of generators. */
INLINE GEN grp_get_gen(GEN G) { return gel(G,1); }
INLINE GEN grp_get_ord(GEN G) { return gel(G,2); }

/* A Quotient Group is a t_VEC [gen,coset]
 * gen (vecvecsmall): coset generators
 * coset (vecsmall): gen[coset[p[1]]] generate the p-coset.
 */
INLINE GEN quo_get_gen(GEN C) { return gel(C,1); }
INLINE GEN quo_get_coset(GEN C) { return gel(C,2); }

static GEN
trivialsubgroups(void)
{ GEN L = cgetg(2, t_VEC); gel(L,1) = trivialgroup(); return L; }

/* Compute the order of p modulo the group given by a set */
long
perm_relorder(GEN p, GEN set)
{
  pari_sp ltop = avma;
  long n = 1, q = p[1];
  while (!F2v_coeff(set,q)) { q = p[q]; n++; }
  return gc_long(ltop,n);
}

GEN
perm_generate(GEN S, GEN H, long o)
{
  long i, n = lg(H)-1;
  GEN L = cgetg(n*o + 1, t_VEC);
  for(i=1; i<=n;     i++) gel(L,i) = vecsmall_copy(gel(H,i));
  for(   ; i <= n*o; i++) gel(L,i) = perm_mul(gel(L,i-n), S);
  return L;
}

/*Return the order (cardinality) of a group */
long
group_order(GEN G)
{
  return zv_prod(grp_get_ord(G));
}

/* G being a subgroup of S_n, output n */
long
group_domain(GEN G)
{
  GEN gen = grp_get_gen(G);
  if (lg(gen) < 2) pari_err_DOMAIN("group_domain", "#G", "=", gen_1,G);
  return lg(gel(gen,1)) - 1;
}

/*Left coset of g mod G: gG*/
GEN
group_leftcoset(GEN G, GEN g)
{
  GEN gen = grp_get_gen(G), ord = grp_get_ord(G);
  GEN res = cgetg(group_order(G)+1, t_VEC);
  long i, j, k;
  gel(res,1) = vecsmall_copy(g);
  k = 1;
  for (i = 1; i < lg(gen); i++)
  {
    long c = k * (ord[i] - 1);
    for (j = 1; j <= c; j++) gel(res,++k) = perm_mul(gel(res,j), gel(gen,i));
  }
  return res;
}
/*Right coset of g mod G: Gg*/
GEN
group_rightcoset(GEN G, GEN g)
{
  GEN gen = grp_get_gen(G), ord = grp_get_ord(G);
  GEN res = cgetg(group_order(G)+1, t_VEC);
  long i, j, k;
  gel(res,1) = vecsmall_copy(g);
  k = 1;
  for (i = 1; i < lg(gen); i++)
  {
    long c = k * (ord[i] - 1);
    for (j = 1; j <= c; j++) gel(res,++k) = perm_mul(gel(gen,i), gel(res,j));
  }
  return res;
}
/*Elements of a group from the generators, cf group_leftcoset*/
GEN
group_elts(GEN G, long n)
{
  GEN gen = grp_get_gen(G), ord = grp_get_ord(G);
  GEN res = cgetg(group_order(G)+1, t_VEC);
  long i, j, k;
  gel(res,1) = identity_perm(n);
  k = 1;
  for (i = 1; i < lg(gen); i++)
  {
    long c = k * (ord[i] - 1);
    /* j = 1, use res[1] = identity */
    gel(res,++k) = vecsmall_copy(gel(gen,i));
    for (j = 2; j <= c; j++) gel(res,++k) = perm_mul(gel(res,j), gel(gen,i));
  }
  return res;
}

GEN
groupelts_conj_set(GEN elts, GEN p)
{
  long i, j, l = lg(elts), n = lg(p)-1;
  GEN res = zero_F2v(n);
  for(j = 1; j < n; j++)
    if (p[j]==1) break;
  for(i = 1; i < l; i++)
    F2v_set(res, p[mael(elts,i,j)]);
  return res;
}

GEN
groupelts_set(GEN elts, long n)
{
  GEN res = zero_F2v(n);
  long i, l = lg(elts);
  for(i=1; i<l; i++)
    F2v_set(res,mael(elts,i,1));
  return res;
}

/*Elements of a group from the generators, returned as a set (bitmap)*/
GEN
group_set(GEN G, long n)
{
  GEN res = zero_F2v(n);
  pari_sp av = avma;
  GEN elts = group_elts(G, n);
  long i, l = lg(elts);
  for(i=1; i<l; i++)
    F2v_set(res,mael(elts,i,1));
  set_avma(av);
  return res;
}

static int
sgcmp(GEN a, GEN b) { return vecsmall_lexcmp(gel(a,1),gel(b,1)); }

GEN
subgroups_tableset(GEN S, long n)
{
  long i, l = lg(S);
  GEN  v = cgetg(l, t_VEC);
  for(i=1; i<l; i++)
    gel(v,i) = mkvec2(group_set(gel(S,i), n), mkvecsmall(i));
  gen_sort_inplace(v,(void*)sgcmp,cmp_nodata, NULL);
  return v;
}

long
tableset_find_index(GEN tbl, GEN set)
{
  long i = tablesearch(tbl,mkvec2(set,mkvecsmall(0)),sgcmp);
  if (!i) return 0;
  return mael3(tbl,i,2,1);
}

GEN
trivialgroup(void) { retmkvec2(cgetg(1,t_VEC), cgetg(1,t_VECSMALL)); }

/*Cyclic group generated by g of order s*/
GEN
cyclicgroup(GEN g, long s)
{ retmkvec2(mkvec( vecsmall_copy(g) ), mkvecsmall(s)); }

/*Return the group generated by g1,g2 of relative orders s1,s2*/
GEN
dicyclicgroup(GEN g1, GEN g2, long s1, long s2)
{ retmkvec2( mkvec2(vecsmall_copy(g1), vecsmall_copy(g2)),
             mkvecsmall2(s1, s2) ); }

/* return the quotient map G --> G/H */
/*The ouput is [gen,hash]*/
/* gen (vecvecsmall): coset generators
 * coset (vecsmall): vecsmall of coset number) */
GEN
groupelts_quotient(GEN elt, GEN H)
{
  pari_sp ltop = avma;
  GEN  p2, p3;
  long i, j, a = 1;
  long n = lg(gel(elt,1))-1, o = group_order(H);
  GEN  el;
  long le = lg(elt)-1;
  GEN used = zero_F2v(le+1);
  long l = le/o;
  p2 = cgetg(l+1, t_VEC);
  p3 = zero_zv(n);
  el = zero_zv(n);
  for (i = 1; i<=le; i++)
    el[mael(elt,i,1)]=i;
  for (i = 1; i <= l; ++i)
  {
    GEN V;
    while(F2v_coeff(used,a)) a++;
    V = group_leftcoset(H,gel(elt,a));
    gel(p2,i) = gel(V,1);
    for(j=1;j<lg(V);j++)
    {
      long b = el[mael(V,j,1)];
      if (b==0) pari_err_IMPL("group_quotient for a non-WSS group");
      F2v_set(used,b);
    }
    for (j = 1; j <= o; j++)
      p3[mael(V, j, 1)] = i;
  }
  return gerepilecopy(ltop,mkvec2(p2,p3));
}

GEN
group_quotient(GEN G, GEN H)
{
  return groupelts_quotient(group_elts(G, group_domain(G)), H);
}

/*Compute the image of a permutation by a quotient map.*/
GEN
quotient_perm(GEN C, GEN p)
{
  GEN gen = quo_get_gen(C);
  GEN coset = quo_get_coset(C);
  long j, l = lg(gen);
  GEN p3 = cgetg(l, t_VECSMALL);
  for (j = 1; j < l; ++j)
  {
    p3[j] = coset[p[mael(gen,j,1)]];
    if (p3[j]==0) pari_err_IMPL("quotient_perm for a non-WSS group");
  }
  return p3;
}

/* H is a subgroup of G, C is the quotient map G --> G/H
 *
 * Lift a subgroup S of G/H to a subgroup of G containing H */
GEN
quotient_subgroup_lift(GEN C, GEN H, GEN S)
{
  GEN genH = grp_get_gen(H);
  GEN genS = grp_get_gen(S);
  GEN genC = quo_get_gen(C);
  long l1 = lg(genH)-1;
  long l2 = lg(genS)-1, j;
  GEN p1 = cgetg(3, t_VEC), L = cgetg(l1+l2+1, t_VEC);
  for (j = 1; j <= l1; ++j) gel(L,j) = gel(genH,j);
  for (j = 1; j <= l2; ++j) gel(L,l1+j) = gel(genC, mael(genS,j,1));
  gel(p1,1) = L;
  gel(p1,2) = vecsmall_concat(grp_get_ord(H), grp_get_ord(S));
  return p1;
}

/* Let G a group and C a quotient map G --> G/H
 * Assume H is normal, return the group G/H */
GEN
quotient_group(GEN C, GEN G)
{
  pari_sp ltop = avma;
  GEN Qgen, Qord, Qelt, Qset, Q;
  GEN Cgen = quo_get_gen(C);
  GEN Ggen = grp_get_gen(G);
  long i,j, n = lg(Cgen)-1, l = lg(Ggen);
  Qord = cgetg(l, t_VECSMALL);
  Qgen = cgetg(l, t_VEC);
  Qelt = mkvec(identity_perm(n));
  Qset = groupelts_set(Qelt, n);
  for (i = 1, j = 1; i < l; ++i)
  {
    GEN  g = quotient_perm(C, gel(Ggen,i));
    long o = perm_relorder(g, Qset);
    gel(Qgen,j) = g;
    Qord[j] = o;
    if (o != 1)
    {
      Qelt = perm_generate(g, Qelt, o);
      Qset = groupelts_set(Qelt, n);
      j++;
    }
  }
  setlg(Qgen,j);
  setlg(Qord,j); Q = mkvec2(Qgen, Qord);
  return gerepilecopy(ltop,Q);
}

GEN
quotient_groupelts(GEN C)
{
  GEN G = quo_get_gen(C);
  long i, l = lg(G);
  GEN Q = cgetg(l, t_VEC);
  for (i = 1; i < l; ++i)
    gel(Q,i) = quotient_perm(C, gel(G,i));
  return Q;
}

/* Return 1 if g normalizes N, 0 otherwise */
long
group_perm_normalize(GEN N, GEN g)
{
  pari_sp ltop = avma;
  long r = gequal(vecvecsmall_sort_shallow(group_leftcoset(N, g)),
                  vecvecsmall_sort_shallow(group_rightcoset(N, g)));
  return gc_long(ltop, r);
}

/* L is a list of subgroups, C is a coset and r a relative order.*/
static GEN
liftlistsubgroups(GEN L, GEN C, long r)
{
  pari_sp ltop = avma;
  long c = lg(C)-1, l = lg(L)-1, n = lg(gel(C,1))-1, i, k;
  GEN R;
  if (!l) return cgetg(1,t_VEC);
  R = cgetg(l*c+1, t_VEC);
  for (i = 1, k = 1; i <= l; ++i)
  {
    GEN S = gel(L,i), Selt = group_set(S,n);
    GEN gen = grp_get_gen(S);
    GEN ord = grp_get_ord(S);
    long j;
    for (j = 1; j <= c; ++j)
    {
      GEN p = gel(C,j);
      if (perm_relorder(p, Selt) == r && group_perm_normalize(S, p))
        gel(R,k++) = mkvec2(vec_append(gen, p),
                            vecsmall_append(ord, r));
    }
  }
  setlg(R, k);
  return gerepilecopy(ltop, R);
}

/* H is a normal subgroup, C is the quotient map G -->G/H,
 * S is a subgroup of G/H, and G is embedded in Sym(l)
 * Return all the subgroups K of G such that
 * S= K mod H and K inter H={1} */
static GEN
liftsubgroup(GEN C, GEN H, GEN S)
{
  pari_sp ltop = avma;
  GEN V = trivialsubgroups();
  GEN Sgen = grp_get_gen(S);
  GEN Sord = grp_get_ord(S);
  GEN Cgen = quo_get_gen(C);
  long n = lg(Sgen), i;
  for (i = 1; i < n; ++i)
  { /*loop over generators of S*/
    GEN W = group_leftcoset(H, gel(Cgen, mael(Sgen, i, 1)));
    V = liftlistsubgroups(V, W, Sord[i]);
  }
  return gerepilecopy(ltop,V);
}

/* 1:A4, 2:S4, 3:F36, 0: other */
long
group_isA4S4(GEN G)
{
  GEN elt = grp_get_gen(G);
  GEN ord = grp_get_ord(G);
  long n = lg(ord);
  if (n != 4 && n != 5) return 0;
  if (n==4 && ord[1]==3 && ord[2]==3 && ord[3]==4)
  {
    long i;
    GEN p = gel(elt,1), q = gel(elt,2), r = gel(elt,3);
    for(i=1; i<=36; i++)
      if (p[r[i]]!=r[q[i]]) return 0;
    return 3;
  }
  if (ord[1]!=2 || ord[2]!=2 || ord[3]!=3) return 0;
  if (perm_commute(gel(elt,1),gel(elt,3))) return 0;
  if (n==4) return 1;
  if (ord[4]!=2) return 0;
  if (perm_commute(gel(elt,3),gel(elt,4))) return 0;
  return 2;
}
/* compute all the subgroups of a group G */
GEN
group_subgroups(GEN G)
{
  pari_sp ltop = avma;
  GEN p1, H, C, Q, M, sg1, sg2, sg3;
  GEN gen = grp_get_gen(G);
  GEN ord = grp_get_ord(G);
  long lM, i, j, n = lg(gen);
  long t;
  if (n == 1) return trivialsubgroups();
  t = group_isA4S4(G);
  if (t == 3)
  {
    GEN H = mkvec2(mkvec3(gel(gen,1), gel(gen,2), perm_sqr(gel(gen,3))),
                   mkvecsmall3(3, 3, 2));
    GEN S = group_subgroups(H);
    GEN V = cgetg(11,t_VEC);
    gel(V,1) = cyclicgroup(gel(gen,3),4);
    for (i=2; i<10; i++)
      gel(V,i) = cyclicgroup(perm_mul(gmael3(V,i-1,1,1),gel(gen,i%3==1 ? 2:1)),4);
    gel(V,10) = G;
    return gerepilecopy(ltop,shallowconcat(S,V));
  }
  else if (t)
  {
    GEN s = gel(gen,1);       /*s = (1,2)(3,4) */
    GEN t = gel(gen,2);       /*t = (1,3)(2,4) */
    GEN st = perm_mul(s, t); /*st = (1,4)(2,3) */
    H = dicyclicgroup(s, t, 2, 2);
    /* sg3 is the list of subgroups intersecting only partially with H*/
    sg3 = cgetg((n==4)?4: 10, t_VEC);
    gel(sg3,1) = cyclicgroup(s, 2);
    gel(sg3,2) = cyclicgroup(t, 2);
    gel(sg3,3) = cyclicgroup(st, 2);
    if (n==5)
    {
      GEN u = gel(gen,3);
      GEN v = gel(gen,4), w, u2;
      if (zv_equal(perm_conj(u,s), t)) /*u=(2,3,4)*/
        u2 = perm_sqr(u);
      else
      {
        u2 = u;
        u = perm_sqr(u);
      }
      if (perm_orderu(v)==2)
      {
        if (!perm_commute(s,v)) /*v=(1,2)*/
        {
          v = perm_conj(u,v);
          if (!perm_commute(s,v)) v = perm_conj(u,v);
        }
        w = perm_mul(v,t); /*w=(1,4,2,3)*/
      }
      else
      {
        w = v;
        if (!zv_equal(perm_sqr(w), s)) /*w=(1,4,2,3)*/
        {
          w = perm_conj(u,w);
          if (!zv_equal(perm_sqr(w), s)) w = perm_conj(u,w);
        }
        v = perm_mul(w,t); /*v=(1,2)*/
      }
      gel(sg3,4) = dicyclicgroup(s,v,2,2);
      gel(sg3,5) = dicyclicgroup(t,perm_conj(u,v),2,2);
      gel(sg3,6) = dicyclicgroup(st,perm_conj(u2,v),2,2);
      gel(sg3,7) = dicyclicgroup(s,w,2,2);
      gel(sg3,8) = dicyclicgroup(t,perm_conj(u,w),2,2);
      gel(sg3,9) = dicyclicgroup(st,perm_conj(u2,w),2,2);
    }
  }
  else
  {
    ulong osig = mael(factoru(ord[1]), 1, 1);
    GEN sig = perm_powu(gel(gen,1), ord[1]/osig);
    H = cyclicgroup(sig,osig);
    sg3 = NULL;
  }
  C = group_quotient(G,H);
  Q = quotient_group(C,G);
  M = group_subgroups(Q); lM = lg(M);
  /* sg1 is the list of subgroups containing H*/
  sg1 = cgetg(lM, t_VEC);
  for (i = 1; i < lM; ++i) gel(sg1,i) = quotient_subgroup_lift(C,H,gel(M,i));
  /*sg2 is a list of lists of subgroups not intersecting with H*/
  sg2 = cgetg(lM, t_VEC);
  /* Loop over all subgroups of G/H */
  for (j = 1; j < lM; ++j) gel(sg2,j) = liftsubgroup(C, H, gel(M,j));
  p1 = gconcat(sg1, shallowconcat1(sg2));
  if (sg3)
  {
    p1 = gconcat(p1, sg3);
    if (n==5) /*ensure that the D4 subgroups of S4 are in supersolvable format*/
      for(j = 3; j <= 5; j++)
      {
        GEN c = gmael(p1,j,1);
        if (!perm_commute(gel(c,1),gel(c,3)))
        {
          if (perm_commute(gel(c,2),gel(c,3))) { swap(gel(c,1), gel(c,2)); }
          else
            perm_mul_inplace2(gel(c,2), gel(c,1));
        }
      }
  }
  return gerepileupto(ltop,p1);
}

/*return 1 if G is abelian, else 0*/
long
group_isabelian(GEN G)
{
  GEN g = grp_get_gen(G);
  long i, j, n = lg(g);
  for(i=2; i<n; i++)
    for(j=1; j<i; j++)
      if (!perm_commute(gel(g,i), gel(g,j))) return 0;
  return 1;
}

/*If G is abelian, return its HNF matrix*/
GEN
group_abelianHNF(GEN G, GEN S)
{
  GEN M, g = grp_get_gen(G), o = grp_get_ord(G);
  long i, j, k, n = lg(g);
  if (!group_isabelian(G)) return NULL;
  if (n==1) return cgetg(1,t_MAT);
  if (!S) S = group_elts(G, group_domain(G));
  M = cgetg(n,t_MAT);
  for(i=1; i<n; i++)
  {
    GEN P, C = cgetg(n,t_COL);
    pari_sp av = avma;
    gel(M,i) = C;
    P = perm_inv(perm_powu(gel(g,i), o[i]));
    for(j=1; j<lg(S); j++)
      if (zv_equal(P, gel(S,j))) break;
    set_avma(av);
    if (j==lg(S)) pari_err_BUG("galoisisabelian [inconsistent group]");
    j--;
    for(k=1; k<i; k++)
    {
      long q = j / o[k];
      gel(C,k) = stoi(j - q*o[k]);
      j = q;
    }
    gel(C,k) = stoi(o[i]);
    for (k++; k<n; k++) gel(C,k) = gen_0;
  }
  return M;
}

/*If G is abelian, return its abstract SNF matrix*/
GEN
group_abelianSNF(GEN G, GEN L)
{
  pari_sp ltop = avma;
  GEN H = group_abelianHNF(G,L);
  if (!H) return NULL;
  return gerepileupto(ltop, smithclean( ZM_snf(H) ));
}

GEN
abelian_group(GEN v)
{
  long card = zv_prod(v), i, d = 1, l = lg(v);
  GEN G = cgetg(3,t_VEC), gen = cgetg(l,t_VEC);
  gel(G,1) = gen;
  gel(G,2) = vecsmall_copy(v);
  for(i=1; i<l; i++)
  {
    GEN p = cgetg(card+1, t_VECSMALL);
    long o = v[i], u = d*(o-1), j, k, l;
    gel(gen, i) = p;
    /* The following loop is over-optimized. Remember that I wrote it for
     * testpermutation. Something has survived... BA */
    for(j=1;j<=card;)
    {
      for(k=1;k<o;k++)
        for(l=1;l<=d; l++,j++) p[j] = j+d;
      for (l=1; l<=d; l++,j++) p[j] = j-u;
    }
    d += u;
  }
  return G;
}

static long
groupelts_subgroup_isnormal(GEN G, GEN H)
{
  long i, n = lg(G);
  for(i = 1; i < n; i++)
    if (!group_perm_normalize(H, gel(G,i))) return 0;
  return 1;
}

/*return 1 if H is a normal subgroup of G*/
long
group_subgroup_isnormal(GEN G, GEN H)
{
  if (lg(grp_get_gen(H)) > 1 && group_domain(G) != group_domain(H))
    pari_err_DOMAIN("group_subgroup_isnormal","domain(H)","!=",
                    strtoGENstr("domain(G)"), H);
  return groupelts_subgroup_isnormal(grp_get_gen(G), H);
}

static GEN
group_subgroup_kernel_set(GEN G, GEN H)
{
  pari_sp av;
  GEN g = grp_get_gen(G);
  long i, n = lg(g);
  GEN S, elts;
  long d = group_domain(G);
  if (lg(grp_get_gen(H)) > 1 && group_domain(G) != group_domain(H))
    pari_err_DOMAIN("group_subgroup_isnormal","domain(H)","!=",
                    strtoGENstr("domain(G)"), H);
  elts = group_elts(H,d);
  S = groupelts_set(elts, d);
  av = avma;
  for(i=1; i<n; i++)
  {
    F2v_and_inplace(S, groupelts_conj_set(elts,gel(g,i)));
    set_avma(av);
  }
  return S;
}

int
group_subgroup_is_faithful(GEN G, GEN H)
{
  pari_sp av = avma;
  GEN K = group_subgroup_kernel_set(G,H);
  F2v_clear(K,1);
  return gc_long(av, F2v_equal0(K));
}

long
groupelts_exponent(GEN elts)
{
  long i, n = lg(elts)-1, expo = 1;
  for(i=1; i<=n; i++) expo = ulcm(expo, perm_orderu(gel(elts,i)));
  return expo;
}

GEN
groupelts_center(GEN S)
{
  pari_sp ltop = avma;
  long i, j, n = lg(S)-1, l = n;
  GEN V, elts = zero_F2v(n+1);
  for(i=1; i<=n; i++)
  {
    if (F2v_coeff(elts,i)) { l--;  continue; }
    for(j=1; j<=n; j++)
      if (!perm_commute(gel(S,i),gel(S,j)))
      {
        F2v_set(elts,i);
        F2v_set(elts,j); l--; break;
      }
  }
  V = cgetg(l+1,t_VEC);
  for (i=1, j=1; i<=n ;i++)
    if (!F2v_coeff(elts,i)) gel(V,j++) = vecsmall_copy(gel(S,i));
  return gerepileupto(ltop,V);
}

GEN
groupelts_conjclasses(GEN elts, long *pnbcl)
{
  long i, j, cl = 0, n = lg(elts)-1;
  GEN c = const_vecsmall(n,0);
  pari_sp av = avma;
  for (i=1; i<=n; i++)
  {
    GEN g = gel(elts,i);
    if (c[i]) continue;
    c[i] = ++cl;
    for(j=1; j<=n; j++)
      if (j != i)
      {
        GEN h = perm_conj(gel(elts,j), g);
        long i2 = gen_search(elts,h,(void*)&vecsmall_lexcmp,&cmp_nodata);
        c[i2] = cl; set_avma(av);
      }
  }
  if (pnbcl) *pnbcl = cl;
  return c;
}

GEN
conjclasses_repr(GEN conj, long nb)
{
  long i, l = lg(conj);
  GEN e = const_vecsmall(nb, 0);
  for(i=1; i<l; i++)
  {
    long ci = conj[i];
    if (!e[ci]) e[ci] = i;
  }
  return e;
}

/* elts of G sorted wrt vecsmall_lexcmp order: g in G is determined by g[1]
 * so sort by increasing g[1] */
static GEN
galois_elts_sorted(GEN gal)
{
  long i, l;
  GEN elts = gal_get_group(gal), v = cgetg_copy(elts, &l);
  for (i = 1; i < l; i++) { GEN g = gel(elts,i); gel(v, g[1]) = g; }
  return v;
}
GEN
group_to_cc(GEN G)
{
  GEN elts = checkgroupelts(G), z = cgetg(5,t_VEC);
  long n, flag = 1;
  if (typ(gel(G,1)) == t_POL)
    elts = galois_elts_sorted(G); /* galoisinit */
  else
  {
    long i, l = lg(elts);
    elts = gen_sort_shallow(elts,(void*)vecsmall_lexcmp,cmp_nodata);
    for (i = 1; i < l; i++)
      if (gel(elts,i)[1] != i) { flag = 0; break; }
  }
  gel(z,1) = elts;
  gel(z,2) = groupelts_conjclasses(elts,&n);
  gel(z,3) = conjclasses_repr(gel(z,2),n);
  gel(z,4) = utoi(flag); return z;
}

/* S a list of generators */
GEN
groupelts_abelian_group(GEN S)
{
  pari_sp ltop = avma;
  GEN Qgen, Qord, Qelt;
  long i, j, n = lg(gel(S,1))-1, l = lg(S);
  Qord = cgetg(l, t_VECSMALL);
  Qgen = cgetg(l, t_VEC);
  Qelt = mkvec(identity_perm(n));
  for (i = 1, j = 1; i < l; ++i)
  {
    GEN  g = gel(S,i);
    long o = perm_relorder(g, groupelts_set(Qelt, n));
    gel(Qgen,j) = g;
    Qord[j] = o;
    if (o != 1) { Qelt = perm_generate(g, Qelt, o); j++; }
  }
  setlg(Qgen,j);
  setlg(Qord,j);
  return gerepilecopy(ltop, mkvec2(Qgen, Qord));
}

GEN
group_export_GAP(GEN G)
{
  pari_sp av = avma;
  GEN s, comma, g = grp_get_gen(G);
  long i, k, l = lg(g);
  if (l == 1) return strtoGENstr("Group(())");
  s = cgetg(2*l, t_VEC);
  comma = strtoGENstr(", ");
  gel(s,1) = strtoGENstr("Group(");
  for (i=1, k=2; i < l; ++i)
  {
    if (i > 1) gel(s,k++) = comma;
    gel(s,k++) = perm_to_GAP(gel(g,i));
  }
  gel(s,k++) = strtoGENstr(")");
  return gerepilecopy(av, shallowconcat1(s));
}

GEN
group_export_MAGMA(GEN G)
{
  pari_sp av = avma;
  GEN s, comma, g = grp_get_gen(G);
  long i, k, l = lg(g);
  if (l == 1) return strtoGENstr("PermutationGroup<1|>");
  s = cgetg(2*l, t_VEC);
  comma = strtoGENstr(", ");
  gel(s,1) = gsprintf("PermutationGroup<%ld|",group_domain(G));
  for (i=1, k=2; i < l; ++i)
  {
    if (i > 1) gel(s,k++) = comma;
    gel(s,k++) = GENtoGENstr( vecsmall_to_vec(gel(g,i)) );
  }
  gel(s,k++) = strtoGENstr(">");
  return gerepilecopy(av, shallowconcat1(s));
}

GEN
group_export(GEN G, long format)
{
  switch(format)
  {
  case 0: return group_export_GAP(G);
  case 1: return group_export_MAGMA(G);
  }
  pari_err_FLAG("galoisexport");
  return NULL; /*-Wall*/
}

static GEN
groupelts_cyclic_subgroups(GEN G)
{
  pari_sp av = avma;
  long i, j, n = lg(G)-1;
  GEN elts, f, gen, ord;
  if (n==1) return cgetg(1,t_VEC);
  elts = zero_F2v(lg(gel(G,1))-1);
  gen = cgetg(n+1, t_VECSMALL);
  ord = cgetg(n+1, t_VECSMALL);
  for (i=1, j=1; i<=n; i++)
  {
    long k = 1, o, c = 0;
    GEN p = gel(G, i);
    if (F2v_coeff(elts, p[1])) continue;
    o = perm_orderu(p);
    gen[j] = i; ord[j] = o; j++;
    do
    {
      if (cgcd(o, ++c)==1) F2v_set(elts, p[k]);
      k = p[k];
    } while (k!=1);
  }
  setlg(gen, j);
  setlg(ord, j);
  f = vecsmall_indexsort(ord);
  return gerepilecopy(av, mkvec2(vecpermute(gen, f), vecpermute(ord, f)));
}

GEN
groupelts_to_group(GEN G)
{
  pari_sp av = avma;
  GEN L, cyc, ord;
  long i, l, n = lg(G)-1;
  if (n==1) return trivialgroup();
  L = groupelts_cyclic_subgroups(G);
  cyc = gel(L,1); ord = gel(L,2);
  l = lg(cyc);
  for (i = l-1; i >= 2; i--)
  {
    GEN p = gel(G,cyc[i]);
    long o = ord[i];
    GEN H = cyclicgroup(p, o);
    if (o == n) return gerepileupto(av, H);
    if (groupelts_subgroup_isnormal(G, H))
    {
      GEN C = groupelts_quotient(G, H);
      GEN Q = quotient_groupelts(C);
      GEN R = groupelts_to_group(Q);
      if (!R) return gc_NULL(av);
      return gerepilecopy(av, quotient_subgroup_lift(C, H, R));
    }
  }
  if (n==12 && l==9 && ord[2]==2 && ord[3]==2 && ord[5]==3)
    return gerepilecopy(av,
      mkvec2(mkvec3(gel(G,cyc[2]), gel(G,cyc[3]), gel(G,cyc[5])), mkvecsmall3(2,2,3)));
  if (n==24 && l==18 && ord[11]==3 && ord[15]==4 && ord[16]==4)
  {
    GEN t21 = perm_sqr(gel(G,cyc[15]));
    GEN t22 = perm_sqr(gel(G,cyc[16]));
    GEN s = perm_mul(t22, gel(G,cyc[15]));
    return gerepilecopy(av,
      mkvec2(mkvec4(t21,t22, gel(G,cyc[11]), s), mkvecsmall4(2,2,3,2)));
  }
  if (n==36 && l==24 && ord[11]==3 && ord[15]==4)
  {
    GEN t1 = gel(G,cyc[11]), t3 = gel(G,cyc[15]);
    return gerepilecopy(av,
      mkvec2(mkvec3(perm_conj(t3, t1), t1, t3), mkvecsmall3(3,3,4)));
  }
  return gc_NULL(av);
}

static GEN
subg_get_gen(GEN subg) {  return gel(subg, 1); }

static GEN
subg_get_set(GEN subg) {  return gel(subg, 2); }

static GEN
groupelt_subg_normalize(GEN elt, GEN subg, GEN cyc)
{
  GEN gen = subg_get_gen(subg), set =  subg_get_set(subg);
  long i, j, u, n = lg(elt)-1, lgen = lg(gen);
  GEN b = F2v_copy(cyc), res = zero_F2v(n);
  for(i = 1; i <= n; i++)
  {
    GEN g;
    if (!F2v_coeff(b, i)) continue;
    g = gel(elt,i);
    for(u=1; u<=n; u++)
      if (g[u]==1) break;
    for(j=1; j<lgen; j++)
    {
      GEN h = gel(elt,gen[j]);
      if (!F2v_coeff(set,g[h[u]])) break;
    }
    if (j < lgen) continue;
    F2v_set(res,i);
    for(j=1; j <= n; j++)
      if (F2v_coeff(set, j))
        F2v_clear(b,g[gel(elt,j)[1]]);
  }
  return res;
}

static GEN
triv_subg(GEN elt)
{
  GEN v = cgetg(3, t_VEC);
  gel(v,1) = cgetg(1,t_VECSMALL);
  gel(v,2) = zero_F2v(lg(elt)-1);
  F2v_set(gel(v,2),1);
  return v;
}

static GEN
subg_extend(GEN U, long e, long o, GEN elt)
{
  long i, j, n = lg(elt)-1;
  GEN g = gel(elt, e);
  GEN gen = vecsmall_append(subg_get_gen(U), e);
  GEN set = subg_get_set(U);
  GEN Vset = zv_copy(set);
  for(i = 1; i <= n; i++)
    if (F2v_coeff(set, i))
    {
      long h = gel(elt, i)[1];
      for(j = 1; j < o; j++)
      {
        h = g[h];
        F2v_set(Vset, h);
      }
    }
  return mkvec2(gen, Vset);
}

static GEN
cyclic_subg(long e, long o, GEN elt)
{
  long j, n = lg(elt)-1, h = 1;
  GEN g = gel(elt, e);
  GEN gen = mkvecsmall(e);
  GEN set = zero_F2v(n);
  F2v_set(set,1);
  for(j = 1; j < o; j++)
  {
    h = g[h];
    F2v_set(set, h);
  }
  return mkvec2(gen, set);
}

static GEN
groupelts_to_regular(GEN elt)
{
  long i, j, n = lg(elt)-1;
  GEN V = cgetg(n+1,t_VEC);
  for (i=1; i<=n; i++)
  {
    pari_sp av = avma;
    GEN g = gel(elt, i);
    GEN W = cgetg(n+1,t_VEC);
    for(j=1; j<=n; j++)
      gel(W,j) = perm_mul(g, gel(elt,j));
    gel(V, i) = gerepileuptoleaf(av,vecvecsmall_indexsort(W));
  }
  vecvecsmall_sort_inplace(V, NULL);
  return V;
}

static long
groupelts_pow(GEN elt, long j, long n)
{
  GEN g = gel(elt,j);
  long i, h = 1;
  for (i=1; i<=n; i++)
    h = g[h];
  return h;
}

static GEN
groupelts_cyclic_primepow(GEN elt, GEN *pt_pr, GEN *pt_po)
{
  GEN R = groupelts_cyclic_subgroups(elt);
  GEN gen = gel(R,1), ord = gel(R,2);
  long i, n = lg(elt)-1, l = lg(gen);
  GEN set = zero_F2v(n);
  GEN pr  = zero_Flv(n);
  GEN po  = zero_Flv(n);
  for (i = 1; i < l; i++)
  {
    long h = gen[i];
    ulong p;
    if (uisprimepower(ord[i], &p))
    {
      F2v_set(set, h);
      uel(pr,h) = p;
      po[h] = groupelts_pow(elt, h, p);
    }
  }
  *pt_pr = pr; *pt_po = po;
  return set;
}

static GEN
all_cyclic_subg(GEN pr, GEN po, GEN elt)
{
  long i, n = lg(pr)-1, m = 0, k = 1;
  GEN W;
  for (i=1; i <= n; i++)
    m += po[i]==1;
  W = cgetg(m+1, t_VEC);
  for (i=1; i <= n; i++)
    if (po[i]==1)
      gel(W, k++) = cyclic_subg(i, pr[i], elt);
  return W;
}

static GEN
groupelts_subgroups_raw(GEN elts)
{
  pari_sp av = avma;
  GEN elt = groupelts_to_regular(elts);
  GEN pr, po, cyc = groupelts_cyclic_primepow(elt, &pr, &po);
  long n = lg(elt)-1;
  long i, j, nS = 1;
  GEN S, L;
  S = cgetg(1+bigomegau(n)+1, t_VEC);
  gel(S, nS++) = mkvec(triv_subg(elt));
  gel(S, nS++) = L = all_cyclic_subg(pr, po, elt);
  if (DEBUGLEVEL) err_printf("subgroups: level %ld: %ld\n",nS,lg(L)-1);
  while (lg(L) > 1)
  {
    pari_sp av2 = avma;
    long nW = 1, lL = lg(L);
    long ng = n;
    GEN W = cgetg(1+ng, t_VEC);
    for (i=1; i<lL; i++)
    {
      GEN U = gel(L, i), set = subg_get_set(U);
      GEN G = groupelt_subg_normalize(elt, U, cyc);
      for (j=1; j<nW; j++)
      {
        GEN Wj = subg_get_set(gel(W, j));
        if (F2v_subset(set, Wj))
          F2v_negimply_inplace(G, Wj);
      }
      for (j=1; j<=n; j++)
        if(F2v_coeff(G,j))
        {
          long p = pr[j];
          if (F2v_coeff(set, j)) continue;
          if (F2v_coeff(set, po[j]))
          {
            GEN U2 = subg_extend(U, j, p, elt);
            F2v_negimply_inplace(G, subg_get_set(U2));
            if (nW > ng) { ng<<=1; W = vec_lengthen(W, ng); }
            gel(W, nW++) = U2;
          }
        }
    }
    setlg(W, nW);
    if (DEBUGLEVEL) err_printf("subgroups: level %ld: %ld\n",nS,nW-1);
    L = W;
    if (nW > 1) gel(S, nS++) = L = gerepilecopy(av2, W);
  }
  setlg(S, nS);
  return gerepilecopy(av, shallowconcat1(S));
}

static GEN
set_groupelts(GEN S, GEN x)
{
  long i, n = F2v_hamming(x), k=1, m = x[1];
  GEN v = cgetg(n+1, t_VEC);
  for (i=1; i<=m; i++)
    if (F2v_coeff(x,i))
      gel(v,k++) = gel(S,i);
  return v;
}

static GEN
subg_to_elts(GEN S, GEN x)
{ pari_APPLY_type(t_VEC, set_groupelts(S, gmael(x,i,2))); }

GEN
groupelts_solvablesubgroups(GEN G)
{
  pari_sp av = avma;
  GEN S = vecvecsmall_sort(checkgroupelts(G));
  GEN L = groupelts_subgroups_raw(S);
  return gerepilecopy(av, subg_to_elts(S, L));
}
