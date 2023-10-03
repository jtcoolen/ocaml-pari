/* Copyright (C) 2015  The PARI group.

This file is part of the PARI package.

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

#define tvalue(i)  gmael(t,(i),1)
#define tleft(i)   mael3(t,(i),2,1)
#define tright(i)  mael3(t,(i),2,2)
#define theight(i) mael3(t,(i),2,3)

static GEN
treesearch(GEN T, GEN x)
{
  long i = 1;
  GEN t = list_data(T);
  if (!t || lg(t)==1) return NULL;
  while (i)
  {
    long c = cmp_universal(x, gel(tvalue(i),1));
    if (!c) return tvalue(i);
    i = c < 0 ? tleft(i): tright(i);
  }
  return NULL;
}

static long
treeparent_r(GEN t, GEN x, long i, long parent)
{
  long c;
  if (i==0) return parent;
  c = cmp_universal(x, gel(tvalue(i),1));
  if (c < 0)
    return treeparent_r(t,x,tleft(i),i);
  else if (c > 0)
    return treeparent_r(t,x,tright(i),i);
  else
    return parent;
}

static void
treekeys(GEN t, long i, GEN V, long *n)
{
  if (i==0) return;
  treekeys(t, tleft(i), V, n);
  gel(V, ++*n) = gel(tvalue(i),1);
  treekeys(t, tright(i), V, n);
}

GEN
mapdomain_shallow(GEN T)
{
  GEN V, t = list_data(T);
  long n = 0;
  if (!t || lg(t)==1) return cgetg(1, t_VEC);
  V = cgetg(lg(t), t_VEC); treekeys(t, 1, V, &n); return V;
}

static void
treemat(GEN t, long i, GEN V, long *n)
{
  if (i==0) return;
  treemat(t, tleft(i), V, n);
  ++*n;
  gmael(V, 1, *n) = gel(tvalue(i), 1);
  gmael(V, 2, *n) = gel(tvalue(i), 2);
  treemat(t, tright(i), V, n);
}

GEN
maptomat_shallow(GEN T)
{
  GEN V, t = list_data(T);
  long n = 0;
  if (!t || lg(t)==1) return cgetg(1, t_MAT);
  V = cgetg(3, t_MAT);
  gel(V,1) = cgetg(lg(t), t_COL);
  gel(V,2) = cgetg(lg(t), t_COL);
  treemat(t, 1, V, &n); return V;
}

static void
treemap_i_r(GEN t, long i, long a, long c, GEN p, GEN M)
{
  long b = (a+c)>>1;
  GEN x = mkvec2(gcopy(gmael(M, 1, p[b])), gcopy(gmael(M, 2, p[b])));
  if (a == c)
    gel(t, i) = mkvec2(x, mkvecsmall3(0, 0, 1));
  else if (a+1 == c)
  {
    treemap_i_r(t, i+1, a+1, c, p, M);
    gel(t, i) = mkvec2(x, mkvecsmall3(0, i+1, theight(i+1) + 1));
  }
  else
  {
    long l = i+1, r = l + b - a, h;
    treemap_i_r(t, l, a, b-1, p, M);
    treemap_i_r(t, r, b+1, c, p, M);
    h = maxss(theight(l), theight(r))+1;
    gel(t, i) = mkvec2(x, mkvecsmall3(l, r, h));
  }
}

static void
treemap_i(GEN t, GEN p, GEN M) { treemap_i_r(t, 1, 1, lg(p)-1, p, M); }

#define value(i)  gmael(list_data(T),(i),1)
#define left(i)   mael3(list_data(T),(i),2,1)
#define right(i)  mael3(list_data(T),(i),2,2)
#define height(i) mael3(list_data(T),(i),2,3)

static long
treeheight(GEN T, long i) { return i? height(i): 0; }

static void
change_leaf(GEN T, GEN x, long p)
{
  pari_sp av = avma;
  listput(T, mkvec2(x, gmael(list_data(T), p, 2)), p);
  set_avma(av);
}

static long
new_leaf(GEN T, GEN x)
{
  pari_sp av = avma;
  listput(T, mkvec2(x, mkvecsmall3(0,0,1)), 0);
  return gc_long(av, lg(list_data(T))-1);
}

static void
fix_height(GEN T, long x)
{ height(x) = maxss(treeheight(T,left(x)), treeheight(T,right(x)))+1; }
static long
treebalance(GEN T, long i)
{ return i ? treeheight(T,left(i)) - treeheight(T,right(i)): 0; }

static long
rotright(GEN T, long y)
{
  long x = left(y), t = right(x);
  right(x) = y;
  left(y)  = t;
  fix_height(T, y);
  fix_height(T, x);
  return x;
}

static long
rotleft(GEN T, long x)
{
  long y = right(x), t = left(y);
  left(y)  = x;
  right(x) = t;
  fix_height(T, x);
  fix_height(T, y);
  return y;
}

static long
treeinsert_r(GEN T, GEN x, long i, long *d)
{
  long b, c;
  if (i==0 || !list_data(T) || lg(list_data(T))==1) return new_leaf(T, x);
  c = cmp_universal(gel(x,1), gel(value(i),1));
  if (c < 0)
  {
    long s = treeinsert_r(T, x, left(i), d);
    if (s < 0) return s;
    left(i) = s;
  }
  else if (c > 0)
  {
    long s = treeinsert_r(T, x, right(i), d);
    if (s < 0) return s;
    right(i) = s;
  }
  else return -i;
  fix_height(T, i);
  b = treebalance(T, i);
  if (b > 1)
  {
    if (*d > 0) left(i) = rotleft(T, left(i));
    return rotright(T, i);
  }
  if (b < -1)
  {
    if (*d < 0) right(i) = rotright(T, right(i));
    return rotleft(T, i);
  }
  *d = c; return i;
}

static long
treeinsert(GEN T, GEN x)
{
  long c = 0, r = treeinsert_r(T, x, 1, &c);
  GEN d;
  if (r < 0) return -r;
  if (r == 1) return 0;
  d = list_data(T);
  /* By convention we want the root to be 1 */
  swap(gel(d,1), gel(d,r));
  if (left(1) == 1) left(1) = r;
  else if (right(1) == 1) right(1) = r;
  else pari_err_BUG("treeadd");
  return 0;
}

static long
treedelete_r(GEN T, GEN x, long i, long *dead)
{
  long b, c;
  if (i==0 || !list_data(T) || lg(list_data(T))==1) return -1;
  c = cmp_universal(x, gel(value(i),1));
  if (c < 0)
  {
    long s = treedelete_r(T, x, left(i), dead);
    if (s < 0) return s;
    left(i) = s;
  }
  else if (c > 0)
  {
    long s = treedelete_r(T, x, right(i), dead);
    if (s < 0) return s;
    right(i) = s;
  }
  else
  {
    *dead = i;
    if (left(i)==0 && right(i)==0) return 0;
    else if (left(i)==0) return right(i);
    else if (right(i)==0) return left(i);
    else
    {
      GEN v, d = list_data(T);
      long j = right(i);
      while (left(j)) j = left(j);
      v = gel(value(j), 1);
      right(i) = treedelete_r(T, v, right(i), dead);
      swap(gel(d,i), gel(d,j));
      lswap(left(i),left(j));
      lswap(right(i),right(j));
      lswap(height(i),height(j));
    }
  }
  fix_height(T, i);
  b = treebalance(T, i);
  if (b > 1 && treebalance(T, left(i)) >= 0) return rotright(T, i);
  if (b > 1 && treebalance(T, left(i)) < 0)
  { left(i) = rotleft(T, left(i)); return rotright(T, i); }
  if (b < -1 && treebalance(T, right(i)) <= 0) return rotleft(T,i);
  if (b < -1 && treebalance(T, right(i)) > 0)
  { right(i) = rotright(T, right(i)); return rotleft(T, i); }
  return i;
}

static long
treedelete(GEN T, GEN x)
{
  long dead, l, r = treedelete_r(T, x, 1, &dead);
  GEN d;
  if (r < 0) return 0;
  d = list_data(T); /* != NULL and nonempty */
  if (r > 1)
  { /* By convention we want the root to be 1 */
    swap(gel(d,1), gel(d,r));
    if (left(1) == 1) left(1) = r;
    else if (right(1) == 1) right(1) = r;
    else dead = r;
  }
  /* We want the dead to be last */
  l = lg(d)-1;
  if (dead != l)
  {
    long p = treeparent_r(d, gel(value(l),1), 1, 0);
    if (left(p) == l) left(p) = dead;
    else if (right(p) == l) right(p) = dead;
    else pari_err_BUG("treedelete2");
    swap(gel(d, dead),gel(d, l));
  }
  listpop(T,0); return 1;
}

static int
ismap(GEN T) { return typ(T) == t_LIST && list_typ(T) == t_LIST_MAP; }

void
mapput(GEN T, GEN a, GEN b)
{
  pari_sp av = avma;
  GEN p = mkvec2(a, b);
  long i;
  if (!ismap(T)) pari_err_TYPE("mapput",T);
  i = treeinsert(T, p); if (i) change_leaf(T, p, i);
  set_avma(av);
}

void
mapdelete(GEN T, GEN a)
{
  pari_sp av = avma;
  long s;
  if (!ismap(T)) pari_err_TYPE("mapdelete",T);
  s = treedelete(T, a); set_avma(av);
  if (!s) pari_err_COMPONENT("mapdelete", "not in", strtoGENstr("map"), a);
}

GEN
mapget(GEN T, GEN a)
{
  GEN x;
  if (!ismap(T)) pari_err_TYPE("mapget",T);
  x = treesearch(T, a);
  if (!x) pari_err_COMPONENT("mapget", "not in", strtoGENstr("map"), a);
  return gcopy(gel(x, 2));
}

int
mapisdefined(GEN T, GEN a, GEN *pt_z)
{
  GEN x;
  if (!ismap(T)) pari_err_TYPE("mapisdefined",T);
  x = treesearch(T, a); if (!x) return 0;
  if (pt_z) *pt_z = gcopy(gel(x, 2));
  return 1;
}

GEN
mapdomain(GEN T)
{
  long i, l;
  GEN V;
  if (!ismap(T)) pari_err_TYPE("mapdomain",T);
  V = mapdomain_shallow(T); l = lg(V);
  for (i = 1; i < l; i++) gel(V,i) = gcopy(gel(V,i));
  return V;
}

GEN
maptomat(GEN T)
{
  long i, l;
  GEN V;
  if (!ismap(T)) pari_err_TYPE("maptomat",T);
  V = maptomat_shallow(T); if (lg(V) == 1) return V;
  l = lgcols(V);
  for (i = 1; i < l; i++)
  {
    gcoeff(V,i,1) = gcopy(gcoeff(V,i,1));
    gcoeff(V,i,2) = gcopy(gcoeff(V,i,2));
  }
  return V;
}

GEN
gtomap(GEN x)
{
  if (!x) return mkmap();
  switch(typ(x))
  {
  case t_MAT:
    {
      long l = lg(x);
      GEN M, p;
      if (l == 1 || lgcols(x)==1) return mkmap();
      if (l != 3) pari_err_TYPE("Map",x);
      p = gen_indexsort_uniq(gel(x,1),(void*)&cmp_universal, cmp_nodata);
      l = lgcols(x);
      if (lg(p) != l)
        pari_err_DOMAIN("Map","x","is not",strtoGENstr("one-to-one"),x);
      M = cgetg(3, t_LIST);
      M[1] = evaltyp(t_LIST_MAP); /* do not set list_nmax! */
      list_data(M) = cgetg(l, t_VEC);
      treemap_i(list_data(M), p, x);
      return M;
    }
  default:
    pari_err_TYPE("Map",x);
  }
  return NULL; /* LCOV_EXCL_LINE */
}
