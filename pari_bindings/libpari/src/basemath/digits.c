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

/*********************************************************************/
/**                    DIGITS / SUM OF DIGITS                       **/
/*********************************************************************/
#include "pari.h"
#include "paripriv.h"

/* set v[i] = 1 iff B^i is needed in the digits_dac algorithm */
static void
set_vexp(GEN v, long l)
{
  long m;
  if (v[l]) return;
  v[l] = 1; m = l>>1;
  set_vexp(v, m);
  set_vexp(v, l-m);
}

/* return all needed B^i for DAC algorithm, for lz digits */
static GEN
get_vB(GEN T, long lz, void *E, struct bb_ring *r)
{
  GEN vB, vexp = const_vecsmall(lz, 0);
  long i, l = (lz+1) >> 1;
  vexp[1] = 1;
  vexp[2] = 1;
  set_vexp(vexp, lz);
  vB = zerovec(lz); /* unneeded entries remain = 0 */
  gel(vB, 1) = T;
  for (i = 2; i <= l; i++)
    if (vexp[i])
    {
      long j = i >> 1;
      GEN B2j = r->sqr(E, gel(vB,j));
      gel(vB,i) = odd(i)? r->mul(E, B2j, T): B2j;
    }
  return vB;
}

static void
gen_digits_dac(GEN x, GEN vB, long l, GEN *z,
               void *E, GEN div(void *E, GEN a, GEN b, GEN *r))
{
  GEN q, r;
  long m = l>>1;
  if (l==1) { *z=x; return; }
  q = div(E, x, gel(vB,m), &r);
  gen_digits_dac(r, vB, m, z, E, div);
  gen_digits_dac(q, vB, l-m, z+m, E, div);
}

static GEN
gen_fromdigits_dac(GEN x, GEN vB, long i, long l, void *E,
                   GEN add(void *E, GEN a, GEN b),
                   GEN mul(void *E, GEN a, GEN b))
{
  GEN a, b;
  long m = l>>1;
  if (l==1) return gel(x,i);
  a = gen_fromdigits_dac(x, vB, i, m, E, add, mul);
  b = gen_fromdigits_dac(x, vB, i+m, l-m, E, add, mul);
  return add(E, a, mul(E, b, gel(vB, m)));
}

static GEN
gen_digits_i(GEN x, GEN B, long n, void *E, struct bb_ring *r,
                          GEN (*div)(void *E, GEN x, GEN y, GEN *r))
{
  GEN z, vB;
  if (n==1) retmkvec(gcopy(x));
  vB = get_vB(B, n, E, r);
  z = cgetg(n+1, t_VEC);
  gen_digits_dac(x, vB, n, (GEN*)(z+1), E, div);
  return z;
}

GEN
gen_digits(GEN x, GEN B, long n, void *E, struct bb_ring *r,
                          GEN (*div)(void *E, GEN x, GEN y, GEN *r))
{
  pari_sp av = avma;
  return gerepilecopy(av, gen_digits_i(x, B, n, E, r, div));
}

GEN
gen_fromdigits(GEN x, GEN B, void *E, struct bb_ring *r)
{
  pari_sp av = avma;
  long n = lg(x)-1;
  GEN vB = get_vB(B, n, E, r);
  GEN z = gen_fromdigits_dac(x, vB, 1, n, E, r->add, r->mul);
  return gerepilecopy(av, z);
}

static GEN
_addii(void *data /* ignored */, GEN x, GEN y)
{ (void)data; return addii(x,y); }
static GEN
_sqri(void *data /* ignored */, GEN x) { (void)data; return sqri(x); }
static GEN
_mulii(void *data /* ignored */, GEN x, GEN y)
 { (void)data; return mulii(x,y); }
static GEN
_dvmdii(void *data /* ignored */, GEN x, GEN y, GEN *r)
 { (void)data; return dvmdii(x,y,r); }

static struct bb_ring Z_ring = { _addii, _mulii, _sqri };

/* does not affect stack unless B = NULL */
static GEN
check_basis(GEN B)
{
  if (!B) return utoipos(10);
  if (typ(B)!=t_INT) pari_err_TYPE("digits",B);
  if (abscmpiu(B,2)<0) pari_err_DOMAIN("digits","B","<",gen_2,B);
  return B;
}

/* x has l digits in base B, write them to z[0..l-1], vB[i] = B^i */
static void
digits_dacsmall(GEN x, GEN vB, long l, ulong* z)
{
  pari_sp av = avma;
  GEN q,r;
  long m;
  if (l==1) { *z=itou(x); return; }
  m=l>>1;
  q = dvmdii(x, gel(vB,m), &r);
  digits_dacsmall(q,vB,l-m,z);
  digits_dacsmall(r,vB,m,z+l-m);
  set_avma(av);
}

/* x t_INT */
static GEN
digits_i(GEN x, GEN B)
{
  pari_sp av = avma;
  long lz;
  GEN z;
  B = check_basis(B);
  if (signe(B) < 0) pari_err_DOMAIN("digits","B","<",gen_0,B);
  if (!signe(x))       {set_avma(av); return cgetg(1,t_VEC); }
  if (abscmpii(x,B)<0) {set_avma(av); retmkvec(absi(x)); }
  if (Z_ispow2(B))
  {
    long k = expi(B);
    if (k == 1) return binaire(x);
    if (k >= BITS_IN_LONG) return binary_2k(x, k);
    (void)new_chunk(4*(expi(x) + 2)); /* HACK */
    z = binary_2k_nv(x, k);
    set_avma(av); return Flv_to_ZV(z);
  }
  x = absi_shallow(x);
  lz = logint(x,B) + 1;
  if (lgefint(B) > 3)
  {
    z = gerepileupto(av, gen_digits_i(x, B, lz, NULL, &Z_ring, _dvmdii));
    vecreverse_inplace(z); return z;
  }
  else
  {
    GEN vB = get_vB(B, lz, NULL, &Z_ring);
    (void)new_chunk(3*lz); /* HACK */
    z = zero_zv(lz);
    digits_dacsmall(x,vB,lz,(ulong*)(z+1));
    set_avma(av); return Flv_to_ZV(z);
  }
}
GEN
digits(GEN x, GEN B)
{
  pari_sp av = avma;
  long v = 0;
  if (typ(x) == t_INT) return digits_i(x, B);
  if (typ(x) != t_PADIC || (v = valp(x)) < 0 || (B && !gequal(B, gel(x,2))))
    pari_err_TYPE("digits",x);
  if (!signe(gel(x, 4))) return cgetg(1, t_VEC);
  x = digits_i(gel(x, 4), gel(x, 2));
  vecreverse_inplace(x);
  if (!v) return x;
  return gerepileupto(av, concat(zerovec(v), x));
}

static GEN
fromdigitsu_dac(GEN x, GEN vB, long i, long l)
{
  GEN a, b;
  long m = l>>1;
  if (l==1) return utoi(uel(x,i));
  if (l==2) return addui(uel(x,i), mului(uel(x,i+1), gel(vB, m)));
  a = fromdigitsu_dac(x, vB, i, m);
  b = fromdigitsu_dac(x, vB, i+m, l-m);
  return addii(a, mulii(b, gel(vB, m)));
}

static GEN
fromdigitsu_i(GEN x, GEN B)
{
  long n = lg(x)-1;
  GEN vB;
  if (n == 0) return gen_0;
  vB = get_vB(B, n, NULL, &Z_ring);
  return fromdigitsu_dac(x, vB, 1, n);
}
GEN
fromdigitsu(GEN x, GEN B)
{ pari_sp av = avma; return gerepileuptoint(av, fromdigitsu_i(x, B)); }

static int
ZV_in_range(GEN v, GEN B)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN vi = gel(v, i);
    if (signe(vi) < 0 || cmpii(vi, B) >= 0) return 0;
  }
  return 1;
}
static int
zv_nonnegative(GEN v)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++) if (v[i] < 0) return 0;
  return 1;
}

GEN
fromdigits(GEN x, GEN B)
{
  pari_sp av = avma;
  long tx = typ(x);
  if (tx == t_VECSMALL)
  {
    if (lg(x)==1) return gen_0;
    if (zv_nonnegative(x))
    {
      B = check_basis(B); x = vecsmall_reverse(x);
      return gerepileuptoint(av, fromdigitsu_i(x, B));
    }
    x = zv_to_ZV(x);
  }
  else if (!is_vec_t(tx) || !RgV_is_ZV(x)) pari_err_TYPE("fromdigits",x);
  if (lg(x) == 1) return gen_0;
  B = check_basis(B);
  if (Z_ispow2(B) && ZV_in_range(x, B)) return fromdigits_2k(x, expi(B));
  x = vecreverse(x);
  return gerepileuptoint(av, gen_fromdigits(x, B, NULL, &Z_ring));
}

static const ulong digsum[] ={
  0,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,11,3,4,5,6,7,8,
  9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,
  12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,
  12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,11,3,4,5,6,7,8,
  9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,
  12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,
  12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11,3,
  4,5,6,7,8,9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,
  9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,
  9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,13,14,15,
  16,17,18,19,20,3,4,5,6,7,8,9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,
  11,12,13,14,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,
  12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,
  19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,19,20,21,4,5,6,7,8,9,
  10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,
  12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,
  11,12,13,14,15,16,17,18,19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,
  18,19,20,21,13,14,15,16,17,18,19,20,21,22,5,6,7,8,9,10,11,12,13,14,6,7,8,9,
  10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,
  10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,13,14,15,16,
  17,18,19,20,12,13,14,15,16,17,18,19,20,21,13,14,15,16,17,18,19,20,21,22,14,
  15,16,17,18,19,20,21,22,23,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,
  15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,11,12,13,
  14,15,16,17,18,19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,19,20,
  21,13,14,15,16,17,18,19,20,21,22,14,15,16,17,18,19,20,21,22,23,15,16,17,18,
  19,20,21,22,23,24,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,
  10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,13,14,15,16,
  17,18,19,20,12,13,14,15,16,17,18,19,20,21,13,14,15,16,17,18,19,20,21,22,14,
  15,16,17,18,19,20,21,22,23,15,16,17,18,19,20,21,22,23,24,16,17,18,19,20,21,
  22,23,24,25,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,11,
  12,13,14,15,16,17,18,19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,
  19,20,21,13,14,15,16,17,18,19,20,21,22,14,15,16,17,18,19,20,21,22,23,15,16,
  17,18,19,20,21,22,23,24,16,17,18,19,20,21,22,23,24,25,17,18,19,20,21,22,23,
  24,25,26,9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,
  13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,19,20,21,13,14,15,16,17,18,19,
  20,21,22,14,15,16,17,18,19,20,21,22,23,15,16,17,18,19,20,21,22,23,24,16,17,
  18,19,20,21,22,23,24,25,17,18,19,20,21,22,23,24,25,26,18,19,20,21,22,23,24,
  25,26,27
};

ulong
sumdigitsu(ulong n)
{
  ulong s = 0;
  while (n) { s += digsum[n % 1000]; n /= 1000; }
  return s;
}

/* res=array of 9-digits integers, return sum_{0 <= i < l} sumdigits(res[i]) */
static ulong
sumdigits_block(ulong *res, long l)
{
  ulong s = sumdigitsu(*--res);
  while (--l > 0) s += sumdigitsu(*--res);
  return s;
}

GEN
sumdigits(GEN n)
{
  const long L = (long)(ULONG_MAX / 81);
  pari_sp av = avma;
  ulong *res;
  long l;

  if (typ(n) != t_INT) pari_err_TYPE("sumdigits", n);
  switch(lgefint(n))
  {
    case 2: return gen_0;
    case 3: return utoipos(sumdigitsu(n[2]));
  }
  res = convi(n, &l);
  if (l < L)
  {
    ulong s = sumdigits_block(res, l);
    return gc_utoipos(av, s);
  }
  else /* Huge. Overflows ulong */
  {
    GEN S = gen_0;
    while (l > L)
    {
      S = addiu(S, sumdigits_block(res, L));
      res += L; l -= L;
    }
    if (l)
      S = addiu(S, sumdigits_block(res, l));
    return gerepileuptoint(av, S);
  }
}

GEN
sumdigits0(GEN x, GEN B)
{
  pari_sp av = avma;
  GEN z;
  long lz;

  if (!B) return sumdigits(x);
  if (typ(x) != t_INT) pari_err_TYPE("sumdigits", x);
  B = check_basis(B);
  if (Z_ispow2(B))
  {
    long k = expi(B);
    if (k == 1) return gc_utoi(av,hammingweight(x));
    if (k < BITS_IN_LONG)
    {
      GEN z = binary_2k_nv(x, k);
      if (lg(z)-1 > 1L<<(BITS_IN_LONG-k)) /* may overflow */
        return gerepileuptoint(av, ZV_sum(Flv_to_ZV(z)));
      return gc_utoi(av,zv_sum(z));
    }
    return gerepileuptoint(av, ZV_sum(binary_2k(x, k)));
  }
  if (!signe(x))       { set_avma(av); return gen_0; }
  if (abscmpii(x,B)<0) { set_avma(av); return absi(x); }
  if (absequaliu(B,10))   { set_avma(av); return sumdigits(x); }
  x = absi_shallow(x); lz = logint(x,B) + 1;
  z = gen_digits_i(x, B, lz, NULL, &Z_ring, _dvmdii);
  return gerepileuptoint(av, ZV_sum(z));
}
