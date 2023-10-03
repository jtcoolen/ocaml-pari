/* Copyright (C) 2000, 2012  The PARI group.

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

static long
conginlist(GEN L, GEN g, void *E, long (*in)(void *, GEN ))
{
  pari_sp av = avma;
  long i, l = lg(L);
  GEN gi = ginv(g);
  for (i = 1; i < l; i++)
    if (in(E, gmul(gel(L,i), gi))) break;
  return gc_long(av, i);
}

static GEN
normalise(GEN M)
{
  long sd = signe(gcoeff(M,2,2));
  if (sd < 0 || (!sd && signe(gcoeff(M,1,2)) < 0)) M = ZM_neg(M);
  return M;
}

static void
filln(GEN V, long n, long a, long c)
{
  long i, j;
  for (j = a + 1, i = 1; i < n; i++)
  { /* j != a (mod n) */
    gel(V,i) = mkvecsmall2(c, j);
    if (++j > n) j = 1;
  }
}
/* set v[k+1..k+n-1] or (k == l) append to v; 0 <= a < n */
static GEN
vec_insertn(GEN v, long n, long k, long a, long c)
{
  long i, j, l = lg(v), L = l + n-1;
  GEN V = cgetg(L, t_VEC);
  if (k == l)
  {
    for (i = 1; i < l; i++) gel(V,i) = gel(v,i);
    filln(V + i-1, n, a, c);
  }
  else
  {
    for (i = 1; i <= k; i++) gel(V,i) = gel(v,i);
    filln(V + i-1, n, a, c);
    i += n - 1;
    for (j = k + 1; j < l; j++) gel(V,i++) = gel(v,j);
  }
  return V;
}
/* append the [c,L[i]], i=1..#L to v */
static GEN
vec_appendL(GEN v, GEN L, long c)
{
  long i, j, lv, l = lg(L);
  GEN w;
  if (l == 1) return v;
  lv = lg(v); w = cgetg(lv + l -1, typ(v));
  for (i = 1; i < lv; i++) gel(w,i) = gel(v,i);
  for (j = 1; j < l; i++, j++) gel(w,i) = mkvecsmall2(c, L[j]);
  return w;
}
#define newcoset(g, k, a) \
{ \
  long _c = lg(C); \
  C = vec_append(C, g); \
  M = vec_append(M, zero_zv(n)); \
  L3= vec_appendL(L3, list3, _c); \
  L = vec_appendL(L, list, _c); \
  B = vec_insertn(B, n, k, a % n, _c); \
}

static long
_isin2(GEN L, long m, long a)
{
  pari_sp av = avma;
  long k = RgV_isin(L, mkvecsmall2(m,a));
  return gc_long(av, k? k: lg(L));
}
static void
get2(GEN x, long *a, long *b) { *a = x[1]; *b = x[2]; }

static GEN
denval(GEN g)
{
  GEN a = gcoeff(g,1,1), c = gcoeff(g,2,1);
  return signe(c)? denom_i(gdiv(a,c)): gen_0;
}
/* M * S, S = [0,1;-1,0] */
static GEN
mulS(GEN g)
{
  GEN a = gcoeff(g,1,1), b = gcoeff(g,1,2);
  GEN c = gcoeff(g,2,1), d = gcoeff(g,2,2);
  retmkmat22(negi(b), a, negi(d), c);
}
/* remove extra scales and reduce ast to involution */
static GEN
rectify(GEN V, GEN ast, GEN gam)
{
  long n = lg(V)-1, n1, i, def, m, dec;
  GEN V1, a1, g1, d, inj;
  pari_sp av;

  for(i = 1, def = 0; i <= n; i++)
    if (ast[ast[i]] != i) def++;
  def /= 3;

  if (!def) return mkvec3(V, ast, gam);
  n1 = n + def;
  g1 = cgetg(n1+1, t_VEC);
  V1 = cgetg(n1+1, t_VEC);
  a1 = cgetg(n1+1, t_VECSMALL);
  d = cgetg(def+1, t_VECSMALL);
  av = avma;
  for (i = m = 1; i <= n; i++)
  {
    long i2 = ast[i], i3 = ast[i2];
    if (i2 > i && i3 > i)
    {
      GEN d1 = denval(ZM_mul(gel(gam,i),  gel(V,ast[i])));
      GEN d2 = denval(ZM_mul(gel(gam,i2), gel(V,ast[i2])));
      GEN d3 = denval(ZM_mul(gel(gam,i3), gel(V,ast[i3])));
      if (cmpii(d1,d2) <= 0)
        d[m++] = cmpii(d1,d3) <= 0? i: i3;
      else
        d[m++] = cmpii(d2,d3) <= 0? i2: i3;
    }
  }
  set_avma(av); inj = zero_zv(n);
  for (i = 1; i <= def; i++) inj[d[i]] = 1;
  for (i = 1, dec = 0; i <= n; i++) { dec += inj[i]; inj[i] = i + dec; }
  for (i = 1; i <= n; i++)
    if (ast[ast[i]] == i)
    {
      gel(g1, inj[i]) = gel(gam,i);
      gel(V1, inj[i]) = gel(V,i);
      a1[inj[i]] = inj[ast[i]];
    }
  for (i = 1; i <= def; i++)
  {
    long a = d[i], b = ast[a], c = ast[b];
    GEN igc;

    gel(V1, inj[b]) = gel(V, b);
    gel(g1, inj[b]) = normalise(SL2_inv_shallow(gel(gam,a)));
    a1[inj[b]] = inj[a]-1;

    gel(V1, inj[c]) = gel(V, c);
    gel(g1, inj[c]) = gel(gam, c);
    a1[inj[c]] = inj[a];

    gel(V1, inj[a]-1) = normalise(ZM_mul(gel(gam,a), mulS(gel(V,b))));
    gel(g1, inj[a]-1) = gel(gam, a);
    a1[inj[a]-1] = inj[b];

    igc = SL2_inv_shallow(gel(gam,c));
    gel(V1, inj[a]) = normalise(ZM_mul(igc, mulS(gel(V,c))));
    gel(g1, inj[a]) = normalise(igc);
    a1[inj[a]] = inj[c];
  }
  return mkvec3(V1, a1, g1);
}
static GEN
vecpop(GEN v)
{
  long l = lg(v);
  *v++ = evaltyp(t_VEC)|_evallg(1); /* stackdummy */
  *v = evaltyp(t_VEC)|_evallg(l-1);
  return v;
}

GEN
msfarey(GEN F, void *E, long (*in)(void *, GEN), GEN *pCM)
{
  pari_sp av = avma, av2, av3;
  GEN V = gel(F,1), ast = gel(F,2), gam = gel(F,3), V2, ast2, gam2;
  GEN C, M, L3, L, B, g, list3, list, perm, v2;
  long n = lg(gam)-1, i, k, m, a, l, c, c3;

  list = cgetg(n+1, t_VECSMALL);
  list3 = cgetg(n+1, t_VECSMALL);
  for (i = c = c3 = 1; i <= n; i++)
  {
    long t;
    if (ast[i] == i)
      t = !isintzero(gtrace(gel(gam,i)));
    else
      t = ast[ast[i]] != i;
    if (t) list3[c3++] = i; else list[c++] = i;
  }
  setlg(list, c); setlg(list3, c3);
  if (typ(ast) == t_VEC) ast = ZV_to_zv(ast);
  av2 = avma;
  C = M = L = L3 = cgetg(1, t_VEC);
  B = mkvec(mkvecsmall2(1,1));
  newcoset(matid(2),1,1);
  while(lg(L)-1 + lg(L3)-1)
  {
    while(lg(L3)-1)
    {
      get2(gel(L3,1), &m,&a); L3 = vecpop(L3);
      av3 = avma;
      g = ZM_mul(gel(C,m), gel(gam,a));
      k = conginlist(C, g, E, in);
      gel(M,m)[a] = k;
      if (k < lg(C)) set_avma(av3);
      else
      {
        k = _isin2(B, m, a);
        newcoset(g, k, ast[a]);
        newcoset(ZM_mul(g,gel(gam,ast[a])), k+n-1, ast[ast[a]]);
        B = vecsplice(B, k);
      }
    }
    get2(gel(L,1), &m,&a); L = vecpop(L);
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"msfarey, #L = %ld", lg(L)-1);
      gerepileall(av2, 4, &C, &M, &L, &B); L3 = cgetg(1, t_VEC);
    }
    av3 = avma;
    g = ZM_mul(gel(C,m), gel(gam,a));
    k = conginlist(C, g, E, in);
    gel(M,m)[a] = k; /* class of C[m]*gam[a] */
    if (k < lg(C)) set_avma(av3);
    else
    {
      k = _isin2(B, m, a);
      newcoset(g,k,ast[a]);
      B = vecsplice(B,k);
    }
  }
  vecvecsmall_sort_inplace(B, &perm);
  l = lg(B);
  V2 = cgetg(l, t_VEC);
  gam2 = cgetg(l, t_VEC);
  ast2 = cgetg(l, t_VECSMALL);
  v2 = cgetg(3, t_VECSMALL);
  for (i = 1; i < l; i++)
  {
    long r, j = perm[i];
    GEN ig;
    get2(gel(B,i), &m,&a);
    r = gel(M,m)[a]; ig = SL2_inv_shallow(gel(C,r));
    gel(V2, j) = normalise(ZM_mul(gel(C,m), gel(V,a)));
    gel(gam2, j) = normalise(ZM_mul(ZM_mul(gel(C,m), gel(gam,a)), ig));
    v2[1] = r; v2[2] = ast[a]; k = vecvecsmall_search(B,v2);
    if (k < 0)
      pari_err(e_MISC, "msfarey: H is not a subgroup of PSL_2(Z)");
    ast2[j] = perm[k];
  }
  F = rectify(V2, ast2, gam2);
  if (pCM) *pCM = mkvec2(C,M);
  return gc_all(av, pCM? 2: 1, &F, pCM);
}

GEN
mscosets(GEN G, void *E, long (*in)(void *, GEN))
{
  pari_sp av = avma;
  GEN g, L, M;
  long n = lg(G)-1, i, m, k;
  g = gel(G,1);
  L = mkvec(typ(g) == t_VECSMALL? identity_perm(lg(g)-1): gdiv(g,g));
  M = mkvec(zero_zv(n));
  for (m = 1; m < lg(L); m++)
    for (i = 1; i <= n; i++)
    {
      g = gmul(gel(L,m), gel(G,i));
      mael(M, m, i) = k = conginlist(L, g, E, in);
      if (k > lg(L)-1) { L = vec_append(L,g); M = vec_append(M, zero_zv(n)); }
      if (gc_needed(av,2))
      {
        if (DEBUGMEM>1) pari_warn(warnmem,"mscosets, #L = %ld", lg(L)-1);
        gerepileall(av, 2, &M, &L);
      }
    }
  return gerepilecopy(av, mkvec2(L, M));
}

int
checkfarey_i(GEN F)
{
  GEN V, ast, gam;
  if (typ(F) != t_VEC || lg(F) < 4) return 0;
  V   = gel(F,1);
  ast = gel(F,2);
  gam = gel(F,3);
  if (typ(V) != t_VEC
      || (typ(ast) != t_VECSMALL && (typ(ast) != t_VEC || !RgV_is_ZV(ast)))
      || typ(gam) != t_VEC
      || lg(V) != lg(ast) || lg(ast) != lg(gam)) return 0;
  return 1;
}
static int
check_inH(GEN inH)
{
  return (typ(inH) == t_CLOSURE && closure_arity(inH) == 1
          && !closure_is_variadic(inH));
}
GEN
msfarey0(GEN F, GEN code, GEN *pCM)
{
  if (!checkfarey_i(F)) pari_err_TYPE("msfarey", F);
  if (!check_inH(code)) pari_err_TYPE("msfarey", code);
  return msfarey(F, (void*)code, gp_callbool, pCM);
}
GEN
mscosets0(GEN V, GEN code)
{
  if (typ(V) != t_VEC) pari_err_TYPE("mscosets", V);
  if (!check_inH(code)) pari_err_TYPE("mscosets", code);
  if (lg(V) == 1) pari_err_TYPE("mscosets [trivial group]", V);
  return mscosets(V, (void*)code, gp_callbool);
}
