/* Copyright (C) 2012-2019 The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_mat

/***********************************************************************/
/**                                                                   **/
/**                             F2v                                   **/
/**                                                                   **/
/***********************************************************************/
/* F2v objects are defined as follows:
   An F2v is a t_VECSMALL:
   v[0] = codeword
   v[1] = number of components
   x[2] = a_0...a_31 x[3] = a_32..a_63, etc.  on 32bit
   x[2] = a_0...a_63 x[3] = a_64..a_127, etc. on 64bit

   where the a_i are bits.
*/

int
F2v_equal0(GEN V)
{
  long l = lg(V);
  while (--l > 1)
    if (V[l]) return 0;
  return 1;
}

GEN
F2c_to_ZC(GEN x)
{
  long l = x[1]+1, lx = lg(x);
  GEN  z = cgetg(l, t_COL);
  long i, j, k;
  for (i=2, k=1; i<lx; i++)
    for (j=0; j<BITS_IN_LONG && k<l; j++,k++)
      gel(z,k) = (x[i]&(1UL<<j))? gen_1: gen_0;
  return z;
}
GEN
F2c_to_mod(GEN x)
{
  long l = x[1]+1, lx = lg(x);
  GEN  z = cgetg(l, t_COL);
  GEN _0 = mkintmod(gen_0,gen_2);
  GEN _1 = mkintmod(gen_1,gen_2);
  long i, j, k;
  for (i=2, k=1; i<lx; i++)
    for (j=0; j<BITS_IN_LONG && k<l; j++,k++)
      gel(z,k) = (x[i]&(1UL<<j))? _1: _0;
  return z;
}

/* x[a..b], a <= b */
GEN
F2v_slice(GEN x, long a, long b)
{
  long i,j,k, l = b-a+1;
  GEN z = cgetg(nbits2lg(l), t_VECSMALL);
  z[1] = l;
  for(i=a,k=1,j=BITS_IN_LONG; i<=b; i++,j++)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    if (F2v_coeff(x,i)) z[k] |= 1UL<<j;
  }
  return z;
}
/* x[a..b,], a <= b */
GEN
F2m_rowslice(GEN x, long a, long b)
{ pari_APPLY_same(F2v_slice(gel(x,i),a,b)) }

GEN
F2v_to_Flv(GEN x)
{
  long l = x[1]+1, lx = lg(x);
  GEN  z = cgetg(l, t_VECSMALL);
  long i,j,k;
  for (i=2, k=1; i<lx; i++)
    for (j=0; j<BITS_IN_LONG && k<l; j++,k++)
      z[k] = (x[i]>>j)&1UL;
  return z;
}

GEN
F2m_to_ZM(GEN x) { pari_APPLY_same(F2c_to_ZC(gel(x,i))) }
GEN
F2m_to_mod(GEN x) { pari_APPLY_same(F2c_to_mod(gel(x,i))) }
GEN
F2m_to_Flm(GEN x) { pari_APPLY_same(F2v_to_Flv(gel(x,i))) }

GEN
RgV_F2v_extract_shallow(GEN V, GEN x)
{
  long n = F2v_hamming(x), m = 1;
  long l = x[1]+1, lx = lg(x);
  GEN W = cgetg(n+1, t_VEC);
  long i,j,k;
  for (i=2, k=1; i<lx; i++)
    for (j=0; j<BITS_IN_LONG && k<l; j++,k++)
      if ((x[i]>>j)&1UL)
        gel(W, m++) = gel(V,k);
  return W;
}

GEN
ZV_to_F2v(GEN x)
{
  long i, j, k, l = lg(x)-1;
  GEN z = cgetg(nbits2lg(l), t_VECSMALL);
  z[1] = l;
  for(i=1,k=1,j=BITS_IN_LONG; i<=l; i++,j++)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    if (mpodd(gel(x,i))) z[k] |= 1UL<<j;
  }
  return z;
}

GEN
RgV_to_F2v(GEN x)
{
  long l = lg(x)-1;
  GEN z = cgetg(nbits2lg(l), t_VECSMALL);
  long i,j,k;
  z[1] = l;
  for(i=1,k=1,j=BITS_IN_LONG; i<=l; i++,j++)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    if (Rg_to_F2(gel(x,i))) z[k] |= 1UL<<j;
  }
  return z;
}

GEN
Flv_to_F2v(GEN x)
{
  long l = lg(x)-1;
  GEN z = cgetg(nbits2lg(l), t_VECSMALL);
  long i,j,k;
  z[1] = l;
  for(i=1,k=1,j=BITS_IN_LONG; i<=l; i++,j++)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    if (x[i]&1L) z[k] |= 1UL<<j;
  }
  return z;
}

GEN
ZM_to_F2m(GEN x) { pari_APPLY_same(ZV_to_F2v(gel(x,i))) }
GEN
RgM_to_F2m(GEN x) { pari_APPLY_same(RgV_to_F2v(gel(x,i))) }
GEN
Flm_to_F2m(GEN x) { pari_APPLY_same(Flv_to_F2v(gel(x,i))) }

GEN
const_F2v(long m)
{
  long i, l = nbits2lg(m);
  GEN c = cgetg(l, t_VECSMALL);
  c[1] = m;
  for (i = 2; i < l; i++) uel(c,i) = -1UL;
  if (remsBIL(m)) uel(c,l-1) = (1UL<<remsBIL(m))-1UL;
  return c;
}

/* Allow lg(y)<lg(x) */
void
F2v_add_inplace(GEN x, GEN y)
{
  long n = lg(y);
  long r = (n-2)&7L, q = n-r, i;
  for (i = 2; i < q; i += 8)
  {
    x[  i] ^= y[  i]; x[1+i] ^= y[1+i]; x[2+i] ^= y[2+i]; x[3+i] ^= y[3+i];
    x[4+i] ^= y[4+i]; x[5+i] ^= y[5+i]; x[6+i] ^= y[6+i]; x[7+i] ^= y[7+i];
  }
  switch (r)
  {
    case 7: x[i] ^= y[i]; i++; case 6: x[i] ^= y[i]; i++;
    case 5: x[i] ^= y[i]; i++; case 4: x[i] ^= y[i]; i++;
    case 3: x[i] ^= y[i]; i++; case 2: x[i] ^= y[i]; i++;
    case 1: x[i] ^= y[i]; i++;
  }
}

/* Allow lg(y)<lg(x) */
void
F2v_and_inplace(GEN x, GEN y)
{
  long n = lg(y);
  long r = (n-2)&7L, q = n-r, i;
  for (i = 2; i < q; i += 8)
  {
    x[  i] &= y[  i]; x[1+i] &= y[1+i]; x[2+i] &= y[2+i]; x[3+i] &= y[3+i];
    x[4+i] &= y[4+i]; x[5+i] &= y[5+i]; x[6+i] &= y[6+i]; x[7+i] &= y[7+i];
  }
  switch (r)
  {
    case 7: x[i] &= y[i]; i++; case 6: x[i] &= y[i]; i++;
    case 5: x[i] &= y[i]; i++; case 4: x[i] &= y[i]; i++;
    case 3: x[i] &= y[i]; i++; case 2: x[i] &= y[i]; i++;
    case 1: x[i] &= y[i]; i++;
  }
}

/* Allow lg(y)<lg(x) */
void
F2v_or_inplace(GEN x, GEN y)
{
  long n = lg(y);
  long r = (n-2)&7L, q = n-r, i;
  for (i = 2; i < q; i += 8)
  {
    x[  i] |= y[  i]; x[1+i] |= y[1+i]; x[2+i] |= y[2+i]; x[3+i] |= y[3+i];
    x[4+i] |= y[4+i]; x[5+i] |= y[5+i]; x[6+i] |= y[6+i]; x[7+i] |= y[7+i];
  }
  switch (r)
  {
    case 7: x[i] |= y[i]; i++; case 6: x[i] |= y[i]; i++;
    case 5: x[i] |= y[i]; i++; case 4: x[i] |= y[i]; i++;
    case 3: x[i] |= y[i]; i++; case 2: x[i] |= y[i]; i++;
    case 1: x[i] |= y[i]; i++;
  }
}

/* Allow lg(y)<lg(x) */
void
F2v_negimply_inplace(GEN x, GEN y)
{
  long n = lg(y);
  long r = (n-2)&7L, q = n-r, i;
  for (i = 2; i < q; i += 8)
  {
    x[  i] &= ~y[  i]; x[1+i] &= ~y[1+i]; x[2+i] &= ~y[2+i]; x[3+i] &= ~y[3+i];
    x[4+i] &= ~y[4+i]; x[5+i] &= ~y[5+i]; x[6+i] &= ~y[6+i]; x[7+i] &= ~y[7+i];
  }
  switch (r)
  {
    case 7: x[i] &= ~y[i]; i++; case 6: x[i] &= ~y[i]; i++;
    case 5: x[i] &= ~y[i]; i++; case 4: x[i] &= ~y[i]; i++;
    case 3: x[i] &= ~y[i]; i++; case 2: x[i] &= ~y[i]; i++;
    case 1: x[i] &= ~y[i]; i++;
  }
}

ulong
F2v_dotproduct(GEN x, GEN y)
{
  long i, lx = lg(x);
  ulong c;
  if (lx <= 2) return 0;
  c = uel(x,2) & uel(y,2);
  for (i=3; i<lx; i++) c ^= uel(x,i) & uel(y,i);
#ifdef LONG_IS_64BIT
  c ^= c >> 32;
#endif
  c ^= c >> 16;
  c ^= c >>  8;
  c ^= c >>  4;
  c ^= c >>  2;
  c ^= c >>  1;
  return c & 1;
}

ulong
F2v_hamming(GEN H)
{
  long i, n=0, l=lg(H);
  for (i=2; i<l; i++) n += hammingl(uel(H,i));
  return n;
}

int
F2v_subset(GEN x, GEN y)
{
  long i, n = lg(y);
  for (i = 2; i < n; i ++)
    if ((x[i] & y[i]) != x[i]) return 0;
  return 1;
}

GEN
matid_F2m(long n)
{
  GEN y = cgetg(n+1,t_MAT);
  long i;
  if (n < 0) pari_err_DOMAIN("matid_F2m", "dimension","<",gen_0,stoi(n));
  for (i=1; i<=n; i++) { gel(y,i) = zero_F2v(n); F2v_set(gel(y,i),i); }
  return y;
}

GEN
F2m_row(GEN x, long j)
{
  long i, l = lg(x);
  GEN V = zero_F2v(l-1);
  for(i = 1; i < l; i++)
    if (F2m_coeff(x,j,i))
      F2v_set(V,i);
  return V;
}

GEN
F2m_transpose(GEN x)
{
  long i, dx, lx = lg(x);
  GEN y;
  if (lx == 1) return cgetg(1,t_MAT);
  dx = coeff(x,1,1); y = cgetg(dx+1,t_MAT);
  for (i=1; i<=dx; i++) gel(y,i) = F2m_row(x,i);
  return y;
}

INLINE GEN
F2m_F2c_mul_i(GEN x, GEN y, long lx, long l)
{
  long j;
  GEN z = NULL;

  for (j=1; j<lx; j++)
  {
    if (!F2v_coeff(y,j)) continue;
    if (!z) z = vecsmall_copy(gel(x,j));
    else F2v_add_inplace(z,gel(x,j));
  }
  if (!z) z = zero_F2v(l);
  return z;
}

GEN
F2m_mul(GEN x, GEN y)
{
  long i,j,l,lx=lg(x), ly=lg(y);
  GEN z;
  if (ly==1) return cgetg(1,t_MAT);
  z = cgetg(ly,t_MAT);
  if (lx==1)
  {
    for (i=1; i<ly; i++) gel(z,i) = mkvecsmall(0);
    return z;
  }
  l = coeff(x,1,1);
  for (j=1; j<ly; j++) gel(z,j) = F2m_F2c_mul_i(x, gel(y,j), lx, l);
  return z;
}

GEN
F2m_F2c_mul(GEN x, GEN y)
{
  long l, lx = lg(x);
  if (lx==1) return cgetg(1,t_VECSMALL);
  l = coeff(x,1,1);
  return F2m_F2c_mul_i(x, y, lx, l);
}

static GEN
_F2m_mul(void *data, GEN x, GEN y) { (void) data; return F2m_mul(x,y); }
static GEN
_F2m_sqr(void *data, GEN x) { (void) data; return F2m_mul(x,x); }
GEN
F2m_powu(GEN x, ulong n)
{
  if (!n) return matid(lg(x)-1);
  return gen_powu(x, n,NULL, &_F2m_sqr, &_F2m_mul);
}

static long
F2v_find_nonzero(GEN x0, GEN mask0, long m)
{
  ulong *x = (ulong *)x0+2, *mask = (ulong *)mask0+2, e;
  long i, l = lg(x0)-2;
  for (i = 0; i < l; i++)
  {
    e = *x++ & *mask++;
    if (e) return i*BITS_IN_LONG+vals(e)+1;
  }
  return m+1;
}

/* in place, destroy x */
GEN
F2m_ker_sp(GEN x, long deplin)
{
  GEN y, c, d;
  long i, j, k, r, m, n;

  n = lg(x)-1;
  m = mael(x,1,1); r=0;

  d = cgetg(n+1, t_VECSMALL);
  c = const_F2v(m);
  for (k=1; k<=n; k++)
  {
    GEN xk = gel(x,k);
    j = F2v_find_nonzero(xk, c, m);
    if (j>m)
    {
      if (deplin) {
        GEN v = zero_F2v(n);
        for (i=1; i<k; i++)
          if (F2v_coeff(xk, d[i])) F2v_set(v, i);
        F2v_set(v, k); return v;
      }
      r++; d[k] = 0;
    }
    else
    {
      F2v_clear(c,j); d[k] = j;
      F2v_clear(xk, j);
      for (i=k+1; i<=n; i++)
      {
        GEN xi = gel(x,i);
        if (F2v_coeff(xi,j)) F2v_add_inplace(xi, xk);
      }
      F2v_set(xk, j);
    }
  }
  if (deplin) return NULL;

  y = zero_F2m_copy(n,r);
  for (j=k=1; j<=r; j++,k++)
  {
    GEN C = gel(y,j); while (d[k]) k++;
    for (i=1; i<k; i++)
      if (d[i] && F2m_coeff(x,d[i],k)) F2v_set(C,i);
    F2v_set(C, k);
  }
  return y;
}

/* not memory clean */
GEN
F2m_ker(GEN x) { return F2m_ker_sp(F2m_copy(x), 0); }
GEN
F2m_deplin(GEN x) { return F2m_ker_sp(F2m_copy(x), 1); }

ulong
F2m_det_sp(GEN x) { return !F2m_ker_sp(x, 1); }

ulong
F2m_det(GEN x)
{
  pari_sp av = avma;
  ulong d = F2m_det_sp(F2m_copy(x));
  return gc_ulong(av, d);
}

/* Destroy x */
GEN
F2m_gauss_pivot(GEN x, long *rr)
{
  GEN c, d;
  long i, j, k, r, m, n;

  n = lg(x)-1; if (!n) { *rr=0; return NULL; }
  m = mael(x,1,1); r=0;

  d = cgetg(n+1, t_VECSMALL);
  c = const_F2v(m);
  for (k=1; k<=n; k++)
  {
    GEN xk = gel(x,k);
    j = F2v_find_nonzero(xk, c, m);
    if (j>m) { r++; d[k] = 0; }
    else
    {
      F2v_clear(c,j); d[k] = j;
      for (i=k+1; i<=n; i++)
      {
        GEN xi = gel(x,i);
        if (F2v_coeff(xi,j)) F2v_add_inplace(xi, xk);
      }
    }
  }

  *rr = r; return gc_const((pari_sp)d, d);
}

long
F2m_rank(GEN x)
{
  pari_sp av = avma;
  long r;
  (void)F2m_gauss_pivot(F2m_copy(x),&r);
  return gc_long(av, lg(x)-1 - r);
}

static GEN
F2m_inv_upper_1_ind(GEN A, long index)
{
  pari_sp av = avma;
  long n = lg(A)-1, i = index, j;
  GEN u = const_vecsmall(n, 0);
  u[i] = 1;
  for (i--; i>0; i--)
  {
    ulong m = F2m_coeff(A,i,i+1) & uel(u,i+1); /* j = i+1 */
    for (j=i+2; j<=n; j++) m ^= F2m_coeff(A,i,j) & uel(u,j);
    u[i] = m & 1;
  }
  return gerepileuptoleaf(av, Flv_to_F2v(u));
}
static GEN
F2m_inv_upper_1(GEN A)
{
  long i, l;
  GEN B = cgetg_copy(A, &l);
  for (i = 1; i < l; i++) gel(B,i) = F2m_inv_upper_1_ind(A, i);
  return B;
}

static GEN
F2_get_col(GEN b, GEN d, long li, long aco)
{
  long i, l = nbits2lg(aco);
  GEN u = cgetg(l, t_VECSMALL);
  u[1] = aco;
  for (i = 1; i <= li; i++)
    if (d[i]) /* d[i] can still be 0 if li > aco */
    {
      if (F2v_coeff(b, i))
        F2v_set(u, d[i]);
      else
        F2v_clear(u, d[i]);
    }
  return u;
}

/* destroy a, b */
GEN
F2m_gauss_sp(GEN a, GEN b)
{
  long i, j, k, l, li, bco, aco = lg(a)-1;
  GEN u, d;

  if (!aco) return cgetg(1,t_MAT);
  li = gel(a,1)[1];
  d = zero_Flv(li);
  bco = lg(b)-1;
  for (i=1; i<=aco; i++)
  {
    GEN ai = vecsmall_copy(gel(a,i));
    if (!d[i] && F2v_coeff(ai, i))
      k = i;
    else
      for (k = 1; k <= li; k++) if (!d[k] && F2v_coeff(ai,k)) break;
    /* found a pivot on row k */
    if (k > li) return NULL;
    d[k] = i;

    /* Clear k-th row but column-wise instead of line-wise */
    /*  a_ij -= a_i1*a1j/a_11
       line-wise grouping:  L_j -= a_1j/a_11*L_1
       column-wise:         C_i -= a_i1/a_11*C_1
    */
    F2v_clear(ai,k);
    for (l=1; l<=aco; l++)
    {
      GEN al = gel(a,l);
      if (F2v_coeff(al,k)) F2v_add_inplace(al,ai);
    }
    for (l=1; l<=bco; l++)
    {
      GEN bl = gel(b,l);
      if (F2v_coeff(bl,k)) F2v_add_inplace(bl,ai);
    }
  }
  u = cgetg(bco+1,t_MAT);
  for (j = 1; j <= bco; j++) gel(u,j) = F2_get_col(gel(b,j), d, li, aco);
  return u;
}

GEN
F2m_gauss(GEN a, GEN b)
{
  pari_sp av = avma;
  if (lg(a) == 1) return cgetg(1,t_MAT);
  return gerepileupto(av, F2m_gauss_sp(F2m_copy(a), F2m_copy(b)));
}
GEN
F2m_F2c_gauss(GEN a, GEN b)
{
  pari_sp av = avma;
  GEN z = F2m_gauss(a, mkmat(b));
  if (!z) return gc_NULL(av);
  if (lg(z) == 1) { set_avma(av); return cgetg(1,t_VECSMALL); }
  return gerepileuptoleaf(av, gel(z,1));
}

GEN
F2m_inv(GEN a)
{
  pari_sp av = avma;
  if (lg(a) == 1) return cgetg(1,t_MAT);
  return gerepileupto(av, F2m_gauss_sp(F2m_copy(a), matid_F2m(gel(a,1)[1])));
}

GEN
F2m_invimage_i(GEN A, GEN B)
{
  GEN d, x, X, Y;
  long i, j, nY, nA = lg(A)-1, nB = lg(B)-1;
  x = F2m_ker_sp(shallowconcat(A, B), 0);
  /* AX = BY, Y in strict upper echelon form with pivots = 1.
   * We must find T such that Y T = Id_nB then X T = Z. This exists iff
   * Y has at least nB columns and full rank */
  nY = lg(x)-1;
  if (nY < nB) return NULL;

  /* implicitly: Y = rowslice(x, nA+1, nA+nB), nB rows */
  d = cgetg(nB+1, t_VECSMALL);
  for (i = nB, j = nY; i >= 1; i--, j--)
  {
    for (; j>=1; j--)
      if (F2m_coeff(x,nA+i,j)) { d[i] = j; break; } /* Y[i,j] */
    if (!j) return NULL;
  }
  x = vecpermute(x, d);

  X = F2m_rowslice(x, 1, nA);
  Y = F2m_rowslice(x, nA+1, nA+nB);
  return F2m_mul(X, F2m_inv_upper_1(Y));
}
GEN
F2m_invimage(GEN A, GEN B)
{
  pari_sp av = avma;
  GEN X = F2m_invimage_i(A,B);
  if (!X) return gc_NULL(av);
  return gerepileupto(av, X);
}

GEN
F2m_F2c_invimage(GEN A, GEN y)
{
  pari_sp av = avma;
  long i, l = lg(A);
  GEN M, x;

  if (l==1) return NULL;
  if (lg(y) != lgcols(A)) pari_err_DIM("F2m_F2c_invimage");
  M = cgetg(l+1,t_MAT);
  for (i=1; i<l; i++) gel(M,i) = gel(A,i);
  gel(M,l) = y; M = F2m_ker(M);
  i = lg(M)-1; if (!i) return gc_NULL(av);

  x = gel(M,i);
  if (!F2v_coeff(x,l)) return gc_NULL(av);
  F2v_clear(x, x[1]); x[1]--; /* remove last coord */
  return gerepileuptoleaf(av, x);
}

/*  Block Lanczos algorithm for kernel of sparse matrix (F2Ms)
    Based on lanczos.cpp by Jason Papadopoulos
    <https://github.com/sagemath/FlintQS/blob/master/src/lanczos.cpp>
    Copyright Jason Papadopoulos 2006
    Released under the GNU General Public License v2 or later version.
*/

/* F2Ms are vector of vecsmall which represents nonzero entries of columns
 * F2w are vecsmall whoses entries are columns of a n x BIL F2m
 * F2wB are F2w in the special case where dim = BIL.
 */

#define BIL BITS_IN_LONG

static GEN
F2w_transpose_F2m(GEN x)
{
  long i, j, l = lg(x)-1;
  GEN z = cgetg(BIL+1, t_MAT);
  for (j = 1; j <= BIL; j++)
    gel(z,j) = zero_F2v(l);
  for (i = 1; i <= l; i++)
  {
    ulong xi = uel(x,i);
    for(j = 1; j <= BIL; j++)
      if (xi&(1UL<<(j-1)))
        F2v_set(gel(z, j), i);
  }
  return z;
}

static GEN
F2wB_mul(GEN a, GEN b)
{
  long i, j;
  GEN c = cgetg(BIL+1, t_VECSMALL);
  for (i = 1; i <= BIL; i++)
  {
    ulong s = 0, ai = a[i];
    for (j = 1; ai; j++, ai>>=1)
      if (ai & 1)
        s ^= b[j];
    c[i] = s;
  }
  return c;
}

static void
precompute_F2w_F2wB(GEN x, GEN c)
{
  ulong z, xk;
  ulong i, j, k, index;
  x++; c++;
  for (j = 0; j < (BIL>>3); j++)
  {
    for (i = 0; i < 256; i++)
    {
      k = 0;
      index = i;
      z = 0;
      while (index) {
        xk = x[k];
        if (index & 1)
          z ^= xk;
        index >>= 1;
        k++;
      }
      c[i] = z;
    }
    x += 8; c += 256;
  }
}

static void
F2w_F2wB_mul_add_inplace(GEN v, GEN x, GEN y)
{
  long i, n = lg(y)-1;
  ulong word;
  GEN c = cgetg(1+(BIL<<5), t_VECSMALL);
  precompute_F2w_F2wB(x, c);
  c++;
  for (i = 1; i <= n; i++)
  {
    word = v[i];
    y[i] ^=  c[ 0*256 + ((word>> 0) & 0xff) ]
           ^ c[ 1*256 + ((word>> 8) & 0xff) ]
           ^ c[ 2*256 + ((word>>16) & 0xff) ]
           ^ c[ 3*256 + ((word>>24) & 0xff) ]
#ifdef LONG_IS_64BIT
           ^ c[ 4*256 + ((word>>32) & 0xff) ]
           ^ c[ 5*256 + ((word>>40) & 0xff) ]
           ^ c[ 6*256 + ((word>>48) & 0xff) ]
           ^ c[ 7*256 + ((word>>56)       ) ]
#endif
           ;
  }
}

/* Return x*y~, which is a F2wB */
static GEN
F2w_transmul(GEN x, GEN y)
{
  long i, j, n = lg(x)-1;
  GEN z = zero_zv(BIL);
  pari_sp av = avma;
  GEN c = zero_zv(BIL<<5) + 1;
  GEN xy = z + 1;

  for (i = 1; i <= n; i++)
  {
    ulong xi = x[i];
    ulong yi = y[i];
    c[ 0*256 + ( xi        & 0xff) ] ^= yi;
    c[ 1*256 + ((xi >>  8) & 0xff) ] ^= yi;
    c[ 2*256 + ((xi >> 16) & 0xff) ] ^= yi;
    c[ 3*256 + ((xi >> 24) & 0xff) ] ^= yi;
#ifdef LONG_IS_64BIT
    c[ 4*256 + ((xi >> 32) & 0xff) ] ^= yi;
    c[ 5*256 + ((xi >> 40) & 0xff) ] ^= yi;
    c[ 6*256 + ((xi >> 48) & 0xff) ] ^= yi;
    c[ 7*256 + ((xi >> 56)       ) ] ^= yi;
#endif
  }
  for(i = 0; i < 8; i++)
  {
    ulong a0 = 0, a1 = 0, a2 = 0, a3 = 0;
#ifdef LONG_IS_64BIT
    ulong a4 = 0, a5 = 0, a6 = 0, a7 = 0;
#endif
    for (j = 0; j < 256; j++) {
      if ((j >> i) & 1) {
        a0 ^= c[0*256 + j];
        a1 ^= c[1*256 + j];
        a2 ^= c[2*256 + j];
        a3 ^= c[3*256 + j];
#ifdef LONG_IS_64BIT
        a4 ^= c[4*256 + j];
        a5 ^= c[5*256 + j];
        a6 ^= c[6*256 + j];
        a7 ^= c[7*256 + j];
#endif
      }
    }
    xy[ 0] = a0; xy[ 8] = a1; xy[16] = a2; xy[24] = a3;
#ifdef LONG_IS_64BIT
    xy[32] = a4; xy[40] = a5; xy[48] = a6; xy[56] = a7;
#endif
    xy++;
  }
  return gc_const(av, z);
}

static GEN
identity_F2wB(void)
{
  long i;
  GEN M = cgetg(BIL+1, t_VECSMALL);
  for (i = 1; i <= BIL; i++)
    uel(M,i) = 1UL<<(i-1);
  return M;
}

static GEN
find_nonsingular_sub(GEN t, GEN last_s, GEN *pt_s)
{
  long i, j, dim = 0;
  ulong mask, row_i, row_j;
  long last_dim = lg(last_s)-1;
  GEN s = cgetg(BIL+1, t_VECSMALL);
  GEN M1 = identity_F2wB();
  pari_sp av = avma;
  GEN cols = cgetg(BIL+1, t_VECSMALL);
  GEN M0 = zv_copy(t);

  mask = 0;
  for (i = 1; i <= last_dim; i++)
  {
    cols[BIL + 1 - i] = last_s[i];
    mask |= 1UL<<(last_s[i]-1);
  }
  for (i = j = 1; i <= BIL; i++)
    if (!(mask & (1UL<<(i-1))))
      cols[j++] = i;

  /* compute the inverse of t[][] */

  for (i = 1; i <= BIL; i++)
  {
    mask = 1UL<<(cols[i]-1);
    row_i = cols[i];
    for (j = i; j <= BIL; j++)
    {
      row_j = cols[j];
      if (uel(M0,row_j) & mask)
      {
        swap(gel(M0, row_j), gel(M0, row_i));
        swap(gel(M1, row_j), gel(M1, row_i));
        break;
      }
    }
    if (j <= BIL)
    {
      for (j = 1; j <= BIL; j++)
      {
        row_j = cols[j];
        if (row_i != row_j && (M0[row_j] & mask))
        {
          uel(M0,row_j) ^= uel(M0,row_i);
          uel(M1,row_j) ^= uel(M1,row_i);
        }
      }
      s[++dim] = cols[i];
      continue;
    }
    for (j = i; j <= BIL; j++)
    {
      row_j = cols[j];
      if (uel(M1,row_j) & mask)
      {
        swap(gel(M0, row_j), gel(M0, row_i));
        swap(gel(M1, row_j), gel(M1, row_i));
        break;
      }
    }
    if (j > BIL) return NULL;
    for (j = 1; j <= BIL; j++)
    {
      row_j = cols[j];
      if (row_i != row_j && (M1[row_j] & mask))
      {
        uel(M0,row_j) ^= uel(M0,row_i);
        uel(M1,row_j) ^= uel(M1,row_i);
      }
    }
    M0[row_i] = M1[row_i] = 0;
  }
  mask = 0;
  for (i = 1; i <= dim; i++)
    mask |= 1UL<<(s[i]-1);
  for (i = 1; i <= last_dim; i++)
    mask |= 1UL<<(last_s[i]-1);
  if (mask != (ulong)(-1))
    return NULL; /* Failure */
  setlg(s, dim+1);
  set_avma(av);
  *pt_s = s;
  return M1;
}

/* Compute x * A~ */
static GEN
F2w_F2Ms_transmul(GEN x, GEN A, long nbrow)
{
  long i, j, l = lg(A);
  GEN b = zero_zv(nbrow);
  for (i = 1; i < l; i++)
  {
    GEN c = gel(A,i);
    long lc = lg(c);
    ulong xi = x[i];
    for (j = 1; j < lc; j++)
      b[c[j]] ^= xi;
  }
  return b;
}

/* Compute x * A */
static GEN
F2w_F2Ms_mul(GEN x, GEN A)
{
  long i, j, l = lg(A);
  GEN b = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; i++)
  {
    GEN c = gel(A,i);
    long lc = lg(c);
    ulong s = 0;
    for (j = 1; j < lc; j++)
      s ^= x[c[j]];
    b[i] = s;
  }
  return b;
}

static void
F2wB_addid_inplace(GEN f)
{
  long i;
  for (i = 1; i <= BIL; i++)
    uel(f,i) ^= 1UL<<(i-1);
}

static void
F2w_mask_inplace(GEN f, ulong m)
{
  long i, l = lg(f);
  for (i = 1; i < l; i++)
    uel(f,i) &= m;
}

static GEN
block_lanczos(GEN B, ulong nbrow)
{
  pari_sp av = avma, av2;
  GEN v0, v1, v2, vnext, x, w;
  GEN winv0, winv1, winv2;
  GEN vt_a_v0, vt_a_v1, vt_a2_v0, vt_a2_v1;
  GEN d, e, f, f2, s0;
  long i, iter;
  long n = lg(B)-1;
  long dim0;
  ulong mask0, mask1;
  v1 = zero_zv(n);
  v2 = zero_zv(n);
  vt_a_v1 = zero_zv(BIL);
  vt_a2_v1 = zero_zv(BIL);
  winv1 = zero_zv(BIL);
  winv2 = zero_zv(BIL);
  s0 = identity_zv(BIL);
  mask1 = (ulong)(-1);

  x = random_zv(n);
  w = F2w_F2Ms_mul(F2w_F2Ms_transmul(x, B, nbrow), B);
  v0 = w;
  av2 = avma;
  for (iter=1;;iter++)
  {
    vnext = F2w_F2Ms_mul(F2w_F2Ms_transmul(v0, B, nbrow), B);
    vt_a_v0  = F2w_transmul(v0, vnext);
    if (zv_equal0(vt_a_v0)) break; /* success */
    vt_a2_v0 = F2w_transmul(vnext, vnext);
    winv0 = find_nonsingular_sub(vt_a_v0, s0, &s0);
    if (!winv0) return gc_NULL(av); /* failure */
    dim0 = lg(s0)-1;
    mask0 = 0;
    for (i = 1; i <= dim0; i++)
      mask0 |= 1UL<<(s0[i]-1);
    d = cgetg(BIL+1, t_VECSMALL);
    for (i = 1; i <= BIL; i++)
      d[i] = (vt_a2_v0[i] & mask0) ^ vt_a_v0[i];

    d = F2wB_mul(winv0, d);
    F2wB_addid_inplace(d);
    e = F2wB_mul(winv1, vt_a_v0);
    F2w_mask_inplace(e, mask0);
    f = F2wB_mul(vt_a_v1, winv1);
    F2wB_addid_inplace(f);
    f = F2wB_mul(winv2, f);
    f2 = cgetg(BIL+1, t_VECSMALL);
    for (i = 1; i <= BIL; i++)
      f2[i] = ((vt_a2_v1[i] & mask1) ^ vt_a_v1[i]) & mask0;

    f = F2wB_mul(f, f2);
    F2w_mask_inplace(vnext, mask0);
    F2w_F2wB_mul_add_inplace(v0, d, vnext);
    F2w_F2wB_mul_add_inplace(v1, e, vnext);
    F2w_F2wB_mul_add_inplace(v2, f, vnext);
    d = F2wB_mul(winv0, F2w_transmul(v0, w));
    F2w_F2wB_mul_add_inplace(v0, d, x);
    v2 = v1; v1 = v0; v0 = vnext;
    winv2 = winv1; winv1 = winv0;
    vt_a_v1 = vt_a_v0;
    vt_a2_v1 = vt_a2_v0;
    mask1 = mask0;
    gerepileall(av2, 9, &x, &s0, &v0, &v1, &v2,
                        &winv1, &winv2, &vt_a_v1, &vt_a2_v1);
  }
  if (DEBUGLEVEL >= 5)
    err_printf("Lanczos halted after %ld iterations\n", iter);
  v1 = F2w_F2Ms_transmul(x, B, nbrow);
  v2 = F2w_F2Ms_transmul(v0, B, nbrow);
  x  = shallowconcat(F2w_transpose_F2m(x), F2w_transpose_F2m(v0));
  v1 = shallowconcat(F2w_transpose_F2m(v1), F2w_transpose_F2m(v2));
  s0 = gel(F2m_indexrank(x), 2);
  x = shallowextract(x, s0);
  v1 = shallowextract(v1, s0);
  return F2m_mul(x, F2m_ker(v1));
}

static GEN
F2v_inflate(GEN v, GEN p, long n)
{
  long i, l = lg(p)-1;
  GEN w = zero_F2v(n);
  for (i=1; i<=l; i++)
    if (F2v_coeff(v,i))
      F2v_set(w, p[i]);
  return w;
}

static GEN
F2m_inflate(GEN x, GEN p, long n)
{ pari_APPLY_same(F2v_inflate(gel(x,i), p, n)) }

GEN
F2Ms_ker(GEN M, long nbrow)
{
  pari_sp av = avma;
  long nbcol = lg(M)-1;
  GEN Mp, R, Rp, p;
  if (nbrow <= 640)
    return gerepileupto(av, F2m_ker(F2Ms_to_F2m(M, nbrow)));
  p = F2Ms_colelim(M, nbrow);
  Mp = vecpermute(M, p);
  do
  {
    R = block_lanczos(Mp, nbrow);
  } while(!R);
  Rp = F2m_inflate(R, p, nbcol);
  return gerepilecopy(av, Rp);
}

GEN
F2m_to_F2Ms(GEN M)
{
  long ncol = lg(M)-1;
  GEN B = cgetg(ncol+1, t_MAT);
  long i, j, k;
  for(i = 1; i <= ncol; i++)
  {
    GEN D, V = gel(M,i);
    long n = F2v_hamming(V), l = V[1];
    D = cgetg(n+1, t_VECSMALL);
    for (j=1, k=1; j<=l; j++)
      if( F2v_coeff(V,j))
        D[k++] = j;
    gel(B, i) = D;
  }
  return B;
}

GEN
F2Ms_to_F2m(GEN M, long nrow)
{
  long i, j, l = lg(M);
  GEN B = cgetg(l, t_MAT);
  for(i = 1; i < l; i++)
  {
    GEN Bi = zero_F2v(nrow), Mi = gel(M,i);
    long l = lg(Mi);
    for (j = 1; j < l; j++)
      F2v_set(Bi, Mi[j]);
    gel(B, i) = Bi;
  }
  return B;
}
