/* Copyright (C) 2021 The PARI group.

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

GEN
zero_F3v(long m)
{
  long l = nbits2nlong(2*m);
  GEN v  = const_vecsmall(l+1, 0);
  v[1] = m;
  return v;
}

GEN
zero_F3m_copy(long m, long n)
{
  long i;
  GEN M = cgetg(n+1, t_MAT);
  for (i = 1; i <= n; i++)
    gel(M,i)= zero_F3v(m);
  return M;
}
#define TRITS_IN_LONG (BITS_IN_LONG>>1)
#define TRITS_MASK (ULONG_MAX/3UL)
#define TWOPOTTRITS_IN_LONG (TWOPOTBITS_IN_LONG-1)

ulong
F3v_coeff(GEN x,long v)
{
  long pos = (v-1)>>TWOPOTTRITS_IN_LONG;
  long r = (v-1)&(TRITS_IN_LONG-1);
  ulong u=(ulong)x[2+pos];
  return (u>>(2*r))&3UL;
}

void
F3v_clear(GEN x, long v)
{
  long pos = (v-1)>>TWOPOTTRITS_IN_LONG;
  long r = (v-1)&(TRITS_IN_LONG-1);
  ulong *u=(ulong*)&x[2+pos];
  *u&=~(3UL<<(2*r));
}

void
F3v_set(GEN x, long v, ulong n)
{
  long pos = (v-1)>>TWOPOTTRITS_IN_LONG;
  long r = (v-1)&(TRITS_IN_LONG-1);
  ulong *u=(ulong*)&x[2+pos];
  *u&=~(3UL<<(2*r));
  *u|=(n<<(2*r));
}

INLINE void
F3v_setneg(GEN x, long v)
{
  long pos = (v-1)>>TWOPOTTRITS_IN_LONG;
  long r = (v-1)&(TRITS_IN_LONG-1);
  ulong *u=(ulong*)&x[2+pos];
  if ((*u>>(2*r))&3UL)
    *u^=(3UL<<(2*r));
}

INLINE void
F3m_setneg(GEN x, long a, long b) { F3v_setneg(gel(x,b), a); }

static ulong
bitswap(ulong a)
{
  const ulong m  = TRITS_MASK;
  return ((a&m)<<1)|((a>>1)&m);
}

static ulong
F3_add(ulong a, ulong b)
{
  ulong c = a^b^bitswap(a&b);
  return c&~bitswap(c);
}

static ulong
F3_sub(ulong a, ulong b)
{
  ulong bi = bitswap(b);
  ulong c = a^bi^bitswap(a&bi);
  return c&~bitswap(c);
}

/* Allow lg(y)<lg(x) */
static void
F3v_add_inplace(GEN x, GEN y)
{
  long n = lg(y);
  long i;
  for (i = 2; i < n; i++)
    x[i] = F3_add(x[i], y[i]);
}

/* Allow lg(y)<lg(x) */
static void
F3v_sub_inplace(GEN x, GEN y)
{
  long n = lg(y);
  long i;
  for (i = 2; i < n; i++)
    x[i] = F3_sub(x[i], y[i]);
}

GEN
Flv_to_F3v(GEN x)
{
  long l = lg(x)-1;
  GEN z = cgetg(nbits2lg(2*l), t_VECSMALL);
  long i,j,k;
  z[1] = l;
  for(i=1,k=1,j=BITS_IN_LONG; i<=l; i++,j+=2)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    z[k] |= (uel(x,i)%3)<<j;
  }
  return z;
}

GEN
Flm_to_F3m(GEN x) { pari_APPLY_same(Flv_to_F3v(gel(x,i))) }

GEN
ZV_to_F3v(GEN x)
{
  long l = lg(x)-1;
  GEN z = cgetg(nbits2lg(2*l), t_VECSMALL);
  long i,j,k;
  z[1] = l;
  for(i=1,k=1,j=BITS_IN_LONG; i<=l; i++,j+=2)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    z[k] |= umodiu(gel(x,i),3)<<j;
  }
  return z;
}

GEN
ZM_to_F3m(GEN x) { pari_APPLY_same(ZV_to_F3v(gel(x,i))) }

GEN
RgV_to_F3v(GEN x)
{
  long l = lg(x)-1;
  GEN z = cgetg(nbits2lg(2*l), t_VECSMALL);
  long i,j,k;
  z[1] = l;
  for(i=1,k=1,j=BITS_IN_LONG; i<=l; i++,j+=2)
  {
    if (j==BITS_IN_LONG) { j=0; z[++k]=0; }
    z[k] |= Rg_to_Fl(gel(x,i),3)<<j;
  }
  return z;
}

GEN
RgM_to_F3m(GEN x) { pari_APPLY_same(RgV_to_F3v(gel(x,i))) }

GEN
F3v_to_Flv(GEN x)
{
  long l = x[1]+1, i, j, k;
  GEN z = cgetg(l, t_VECSMALL);
  for (i=2,k=1; i<lg(x); i++)
    for (j=0; j<BITS_IN_LONG && k<l; j+=2,k++)
      z[k] = (uel(x,i)>>j)&3UL;
  return z;
}
GEN
F3c_to_ZC(GEN x)
{
  long l = x[1]+1, i, j, k;
  GEN z = cgetg(l, t_COL);
  for (i=2,k=1; i<lg(x); i++)
    for (j=0; j<BITS_IN_LONG && k<l; j+=2,k++)
      switch((uel(x,i)>>j)&3UL)
      {
      case 0: gel(z,k) = gen_0; break;
      case 1: gel(z,k) = gen_1; break;
      default:gel(z,k) = gen_2; break;
      }
  return z;
}
GEN
F3c_to_mod(GEN x)
{
  long l = x[1]+1, i, j, k;
  GEN z = cgetg(l, t_COL), N = utoipos(3);
  GEN _0 = mkintmod(gen_0, N);
  GEN _1 = mkintmod(gen_1, N);
  GEN _2 = mkintmod(gen_2, N);
  for (i=2,k=1; i<lg(x); i++)
    for (j=0; j<BITS_IN_LONG && k<l; j+=2,k++)
      switch((uel(x,i)>>j)&3UL)
      {
      case 0: gel(z,k) = _0; break;
      case 1: gel(z,k) = _1; break;
      default: gel(z,k)= _2; break;
      }
  return z;
}

GEN
F3m_to_ZM(GEN x) { pari_APPLY_same(F3c_to_ZC(gel(x,i))) }
GEN
F3m_to_mod(GEN x) { pari_APPLY_same(F3c_to_mod(gel(x,i))) }
GEN
F3m_to_Flm(GEN x) { pari_APPLY_same(F3v_to_Flv(gel(x,i))) }

/* in place, destroy x */
GEN
F3m_ker_sp(GEN x, long deplin)
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
    for (j=1; j<=m; j++)
      if (F2v_coeff(c,j) && F3m_coeff(x,j,k)) break;
    if (j>m)
    {
      if (deplin) {
        GEN v = zero_F3v(n);
        for (i=1; i<k; i++) F3v_set(v, i, F3v_coeff(xk, d[i]));
        F3v_set(v, k, 1); return v;
      }
      r++; d[k] = 0;
    }
    else
    {
      ulong xkj = F3v_coeff(xk,j);
      F3v_clear(xk, j);
      F2v_clear(c,j); d[k] = j;
      for (i=k+1; i<=n; i++)
      {
        GEN xi = gel(x,i);
        ulong u = F3v_coeff(xi,j);
        if (u)
        {
          if (u==xkj) F3v_sub_inplace(xi, xk);
          else        F3v_add_inplace(xi, xk);
        }
      }
      F3v_set(xk, j, 2);
      if (xkj==1)
        for (i=k+1; i<=n; i++) F3m_setneg(x,j,i);
    }
  }
  if (deplin) return NULL;
  y = zero_F3m_copy(n,r);
  for (j=k=1; j<=r; j++,k++)
  {
    GEN C = gel(y,j);
    while (d[k]) k++;
    for (i=1; i<k; i++)
      if (d[i]) F3v_set(C,i,F3m_coeff(x,d[i],k));
    F3v_set(C, k, 1);
  }
  return y;
}

GEN
F3m_ker(GEN x) { return F3m_ker_sp(F3m_copy(x), 0); }

INLINE GEN
F3m_F3c_mul_i(GEN x, GEN y, long lx, long l)
{
  long j;
  GEN z = zero_F3v(l);

  for (j=1; j<lx; j++)
  {
    ulong c = F3v_coeff(y,j);
    if (!c) continue;
    if (c==1)
      F3v_add_inplace(z,gel(x,j));
    else
      F3v_sub_inplace(z,gel(x,j));
  }
  return z;
}

GEN
F3m_mul(GEN x, GEN y)
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
  for (j=1; j<ly; j++) gel(z,j) = F3m_F3c_mul_i(x, gel(y,j), lx, l);
  return z;
}

GEN
F3m_row(GEN x, long j)
{
  long i, l = lg(x);
  GEN V = zero_F3v(l-1);
  for(i = 1; i < l; i++) F3v_set(V, i, F3m_coeff(x,j,i));
  return V;
}

GEN
F3m_transpose(GEN x)
{
  long i, l;
  GEN y;
  if (lg(x) == 1) return cgetg(1,t_MAT);
  l = coeff(x,1,1) + 1; y = cgetg(l, t_MAT);
  for (i = 1; i < l; i++) gel(y,i) = F3m_row(x,i);
  return y;
}
