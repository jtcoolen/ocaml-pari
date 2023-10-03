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

#define DEBUGLEVEL DEBUGLEVEL_hensel

/* Assume n > 0. We want to go to accuracy n, starting from accuracy 1, using
 * a quadratically convergent algorithm. Goal: 9 -> 1,2,3,5,9 instead of
 * 1,2,4,8,9 (sequence of accuracies).
 *
 * Let a0 = 1, a1 = 2, a2, ... ak = n, the sequence of accuracies. To obtain
 * it, work backwards:
 *   a(k) = n, a(i-1) = (a(i) + 1) \ 2,
 * but we do not want to store a(i) explicitly, even as a t_VECSMALL, since
 * this would leave an object on the stack. We store a(i) implicitly in a
 * MASK: let a(0) = 1, if the i-bit of MASK is set, set a(i+1) = 2 a(i) - 1,
 * and 2a(i) otherwise.
 *
 * In fact, we do something a little more complicated to simplify the
 * function interface and avoid returning k and MASK separately: we return
 * MASK + 2^(k+1), so the highest bit of the mask indicates the length of the
 * sequence, and the following ones are as above. */

ulong
quadratic_prec_mask(long n)
{
  long a = n, i;
  ulong mask = 0;
  for(i = 1;; i++, mask <<= 1)
  {
    mask |= (a&1); a = (a+1)>>1;
    if (a==1) return mask | (1UL << i);
  }
}

/***********************************************************************/
/**                                                                   **/
/**                            Zp                                     **/
/**                                                                   **/
/***********************************************************************/

static GEN
Zp_divlift(GEN b, GEN a, GEN x, GEN p, long n)
{
  pari_sp ltop = avma, av;
  ulong mask;
  GEN q = p;
  if (n == 1) return gcopy(x);
  mask = quadratic_prec_mask(n);
  av = avma;
  while (mask > 1)
  {
    GEN v, q2 = q;
    q = sqri(q);
    if (mask & 1UL) q = diviiexact(q,p);
    mask >>= 1;
    if (mask > 1 || !b)
    {
      v = Fp_sub(Fp_mul(x, modii(a, q), q), gen_1, q);
      x = Fp_sub(x, Fp_mul(v, x, q), q);
    }
    else
    {
      GEN y = Fp_mul(x, b, q), yt = modii(y, q2);
      v = Fp_sub(Fp_mul(x, modii(a, q), q), gen_1, q);
      x = Fp_sub(y, Fp_mul(v, yt, q), q);
    }
    if (gc_needed(av, 1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"gen_Zp_Newton");
      gerepileall(av, 2, &x, &q);
    }
  }
  return gerepileupto(ltop, x);
}

GEN
Zp_invlift(GEN a, GEN x, GEN p, long e)
{ return Zp_divlift(NULL, a, x, p, e); }

GEN
Zp_inv(GEN a, GEN p, long e)
{
  pari_sp av=avma;
  GEN ai;
  if (lgefint(p)==3)
  {
    ulong pp = p[2];
    ai = utoi(Fl_inv(umodiu(a,pp), pp));
  } else
    ai = Fp_inv(modii(a, p), p);
  return gerepileupto(av, Zp_invlift(a, ai, p, e));
}

GEN
Zp_div(GEN b, GEN a, GEN p, long e)
{
  pari_sp av=avma;
  GEN ai;
  if (lgefint(p)==3)
  {
    ulong pp = p[2];
    ai = utoi(Fl_inv(umodiu(a,pp), pp));
  } else
    ai = Fp_inv(modii(a, p), p);
  return gerepileupto(av, Zp_divlift(b, a, ai, p, e));
}

static GEN
mul2n(void *E, GEN x, GEN y) { return remi2n(mulii(x, y), (long)E); }
static GEN
sqr2n(void *E, GEN x) { return remi2n(sqri(x), (long)E); }

/* a^n mod 2^e using remi2n (result has the same sign as a) */
static GEN
Fp_pow2n(GEN a, GEN n, long e)
{ return gen_pow(a, n, (void*)e, &sqr2n, &mul2n); }

/* Same as ZpX_liftroot for the polynomial X^n-b*/
GEN
Zp_sqrtnlift(GEN b, GEN n, GEN a, GEN p, long e)
{
  pari_sp av = avma;
  int nis2, pis2;
  GEN q, w, n_1;
  ulong mask;

  if (e == 1) return icopy(a);
  nis2 = equaliu(n, 2); n_1 = nis2? NULL: subiu(n,1);
  pis2 = equaliu(p, 2);
  mask = quadratic_prec_mask(e);
  w = nis2 ? shifti(a,1): Fp_mul(n, Fp_pow(a,n_1,p), p);
  w = Fp_inv(w, p);
  q = p; /* q = p^e; use e instead of q iff p = 2 */
  e = 1;
  for(;;)
  {
    if (pis2)
    {
      e <<= 1; if (mask & 1) e--;
      mask >>= 1;
      /* a -= w (a^n - b) */
      a = remi2n(subii(a, mulii(w, subii(Fp_pow2n(a, n, e), b))), e);
      if (mask == 1) break;
      /* w += w - w^2 n a^(n-1)*/
      w = subii(shifti(w,1), remi2n(mulii(remi2n(sqri(w), e),
                                          mulii(n, Fp_pow2n(a, n_1, e))), e));
      continue;
    }
    q = sqri(q); if (mask & 1) q = diviiexact(q, p);
    mask >>= 1;
    if (lgefint(q) == 3 && lgefint(n) == 3)
    {
      ulong Q = uel(q,2), N = uel(n,2);
      ulong A = umodiu(a, Q);
      ulong B = umodiu(b, Q);
      ulong W = umodiu(w, Q);

      A = Fl_sub(A, Fl_mul(W, Fl_sub(Fl_powu(A,N,Q), B, Q), Q), Q);
      a = utoi(A);
      if (mask == 1) break;
      if (nis2)
        W = Fl_double(Fl_sub(W, Fl_mul(Fl_sqr(W,Q), A, Q), Q), Q);
      else
        W = Fl_sub(Fl_double(W,Q),
                   Fl_mul(Fl_sqr(W,Q), Fl_mul(N,Fl_powu(A, N-1, Q), Q), Q), Q);
      w = utoi(W);
    }
    else
    {
      /* a -= w (a^n - b) */
      a = modii(subii(a, mulii(w, subii(Fp_pow(a,n,q), b))), q);
      if (mask == 1) break;
      /* w += w - w^2 n a^(n-1)*/
      if (nis2)
        w = shifti(subii(w, Fp_mul(Fp_sqr(w,q), a, q)), 1);
      else
        w = subii(shifti(w,1), Fp_mul(Fp_sqr(w,q),
                                      mulii(n, Fp_pow(a,n_1,q)), q));
    }
  }
  if (pis2 && signe(a) < 0) a = addii(a, int2n(e));
  return gerepileuptoint(av, a);
}

/* ZpX_liftroot for the polynomial X^2-b */
GEN
Zp_sqrtlift(GEN b, GEN a, GEN p, long e)
{ return Zp_sqrtnlift(b, gen_2, a, p, e); }

GEN
Zp_sqrt(GEN x, GEN p, long e)
{
  pari_sp av;
  GEN z;
  if (absequaliu(p,2)) return Z2_sqrt(x,e);
  av = avma; z = Fp_sqrt(Fp_red(x, p), p);
  if (!z) return NULL;
  if (e > 1) z = Zp_sqrtlift(x, z, p, e);
  return gerepileuptoint(av, z);
}

/* p-adic logarithm of a = 1 mod p (adapted from a C program of Xavier Caruso)
 * Algorithm:
 * 1. raise a at the power p^(v-1)  in order to make it closer to 1
 * 2. write the new a as a product
 *      1 / a = (1 - a_0*p^v) (1 - a_1*p^(2*v) (1 - a_2*p^(4*v) ...
 *    with 0 <= a_i < p^(v*2^i).
 * 3. compute each log(1 - a_i*p^(v*2^i)) using Taylor expansion
 *    and binary splitting */
GEN
Zp_log(GEN a, GEN p, ulong e)
{
  pari_sp av = avma;
  ulong v, N, Np, trunc, pp = itou_or_0(p);
  GEN pe, pv, trunc_mod, num, den, ans = gen_0;

  if (equali1(a)) return ans; /* very frequent! */
  /* First make the argument closer to 1 by raising it to the p^(v-1) */
  v = pp? ulogint(e, pp): 0;  /* v here is v-1 */
  pe = powiu(p,e); pv = powiu(p,v);
  a = Fp_pow(a, pv, mulii(pe, pv));
  e += v;

  /* Where do we truncate the Taylor expansion */
  N = e + v; N /= ++v; /* note the ++v */
  Np = N;
  while(1)
  {
    ulong e = Np;
    if (pp) e += ulogint(N, pp) / v;
    if (e == N) break;
    N = e;
  }

  num = cgetg(N+1, t_INT);
  den = cgetg(N+1, t_INT);
  trunc = v << 1;
  trunc_mod = powiu(p, trunc);
  while(1)
  { /* compute f = 1 - a_i*p^((v+1)*2^i); trunc_mod = p^((v+1)*2^(i+1)) */
    GEN f = modii(a, trunc_mod);
    if (!equali1(f))
    {
      ulong i, step = 1;
      GEN h, hpow;

      f = subui(2, f);
      a = mulii(a, f);

      /* compute the Taylor expansion of log(f), over Q for now */
      for (i = 1; i <= N; i++) { gel(num,i) = gen_1; gel(den,i) = utoipos(i); }
      hpow = h = subui(1, f); /* h = a_i*p^(2^i) */
      while(1)
      {
        for (i = 1; i <= N - step; i += step << 1)
        {
          GEN t = mulii(mulii(hpow, gel(num,i+step)), gel(den,i));
          gel(num,i) = mulii(gel(num,i), gel(den,i+step));
          gel(num,i) = addii(gel(num,i), t);
          gel(den,i) = mulii(gel(den,i), gel(den,i+step));
        }
        step <<= 1; if (step >= N) break;
        hpow = sqri(hpow);
      }

      if (pp)
      { /* simplify the fraction */
        GEN d = powuu(pp, factorial_lval(N, pp));
        gel(num,1) = diviiexact(gel(num,1), d);
        gel(den,1) = diviiexact(gel(den,1), d);
      }

      h = diviiexact(h, pv);
      ans = addii(ans, mulii(mulii(gel(num,1), h), Zp_inv(gel(den,1), p, e)));
    }
    if (trunc > e) break;
    trunc_mod = sqri(trunc_mod); trunc <<= 1; N >>= 1;
  }
  return gerepileuptoint(av, modii(ans, pe));
}

/* p-adic exponential of a = 0 (mod 2p)
 * 1. write a as a sum a = a_0*p + a_1*p^2 + a_2*p^4 + ...
 *    with 0 <= a_i < p^(2^i).
 * 2. compute exp(a_i*p^(2^i)) using Taylor expansion and binary splitting */
GEN
Zp_exp(GEN a, GEN p, ulong e)
{
  pari_sp av = avma;
  ulong trunc, N = e, pp = itou_or_0(p);
  GEN num, den, trunc_mod = NULL, denominator = gen_1, ans = gen_1;
  GEN pe = powiu(p, e);
  int pis2 = pp == 2;

  /* Where do we truncate the Taylor expansion */
  if (!pis2) N += sdivsi(e, subis(p,2));
  num = cgetg(N+2, t_VEC);
  den = cgetg(N+2, t_VEC);
  if (pis2) trunc = 4;
  else
  {
    trunc = 2;
    trunc_mod = sqri(p);
  }
  while(1)
  {
    GEN h, hpow, f = pis2? remi2n(a, trunc): modii(a, trunc_mod);
    a = subii(a, f);
    if (signe(f))
    {
      ulong step = 1, i;
      /* Taylor expansion of exp(f), over Q for now */
      gel(num,1) = gen_1;
      gel(den,1) = gen_1;
      for (i = 2; i <= N+1; i++)
      {
        gel(num,i) = gen_1;
        gel(den,i) = utoipos(i-1);
      }
      hpow = h = f;
      while(1)
      {
        for (i = 1; i <= N+1 - step; i += step << 1) {
          gel(num,i) = mulii(gel(num,i), gel(den,i+step));
          gel(num,i) = addii(gel(num,i), mulii(hpow, gel(num,i+step)));
          gel(den,i) = mulii(gel(den,i), gel(den,i+step));
        }
        step <<= 1; if (step > N) break;
        hpow = sqri(hpow);
      }

      /* We simplify the fraction */
      if (pp)
      {
        GEN d = powuu(pp, factorial_lval(N, pp));
        gel(num,1) = diviiexact(gel(num,1), d);
        gel(den,1) = diviiexact(gel(den,1), d);
      }
      /* We add this contribution to exp(f) */
      ans = Fp_mul(ans, gel(num,1), pe);
      denominator = Fp_mul(denominator, gel(den,1), pe);
    }

    if (trunc > e) break;
    if (!pis2) trunc_mod = sqri(trunc_mod);
    trunc <<= 1; N >>= 1;
  }
  return gerepileuptoint(av, Zp_div(ans, denominator, p, e));
}

/***********************************************************************/
/**                                                                   **/
/**       QUADRATIC HENSEL LIFT (adapted from V. Shoup's NTL)         **/
/**                                                                   **/
/***********************************************************************/
/* Setup for divide/conquer quadratic Hensel lift
 *  a = set of k t_POL in Z[X] = factors over Fp (T=NULL) or Fp[Y]/(T)
 *  V = set of products of factors built as follows
 *  1) V[1..k] = initial a
 *  2) iterate:
 *      append to V the two smallest factors (minimal degree) in a, remove them
 *      from a and replace them by their product [net loss for a = 1 factor]
 *
 * W = bezout coeffs W[i]V[i] + W[i+1]V[i+1] = 1
 *
 * link[i] = -j if V[i] = a[j]
 *            j if V[i] = V[j] * V[j+1]
 * Arrays (link, V, W) pre-allocated for 2k - 2 elements */
static void
BuildTree(GEN link, GEN V, GEN W, GEN a, GEN T, GEN p)
{
  long k = lg(a)-1;
  long i, j, s, minp, mind;

  for (i=1; i<=k; i++) { gel(V,i) = gel(a,i); link[i] = -i; }

  for (j=1; j <= 2*k-5; j+=2,i++)
  {
    minp = j;
    mind = degpol(gel(V,j));
    for (s=j+1; s<i; s++)
      if (degpol(gel(V,s)) < mind) { minp = s; mind = degpol(gel(V,s)); }

    swap(gel(V,j), gel(V,minp));
    lswap(link[j], link[minp]);

    minp = j+1;
    mind = degpol(gel(V,j+1));
    for (s=j+2; s<i; s++)
      if (degpol(gel(V,s)) < mind) { minp = s; mind = degpol(gel(V,s)); }

    swap(gel(V,j+1), gel(V,minp));
    lswap(link[j+1], link[minp]);

    gel(V,i) = FqX_mul(gel(V,j), gel(V,j+1), T, p);
    link[i] = j;
  }

  for (j=1; j <= 2*k-3; j+=2)
  {
    GEN d, u, v;
    d = FqX_extgcd(gel(V,j), gel(V,j+1), T, p, &u, &v);
    if (degpol(d) > 0) pari_err_COPRIME("BuildTree", gel(V,j), gel(V,j+1));
    d = gel(d,2);
    if (!gequal1(d))
    {
      if (typ(d)==t_POL)
      {
        d = FpXQ_inv(d, T, p);
        u = FqX_Fq_mul(u, d, T, p);
        v = FqX_Fq_mul(v, d, T, p);
      }
      else
      {
        d = Fp_inv(d, p);
        u = FqX_Fp_mul(u, d, T,p);
        v = FqX_Fp_mul(v, d, T,p);
      }
    }
    gel(W,j) = u;
    gel(W,j+1) = v;
  }
}

/* au + bv = 1 (p0), ab = f (p0). Lift mod p1 = p0 pd (<= p0^2).
 * If noinv is set, don't lift the inverses u and v */
static void
ZpX_HenselLift(GEN V, GEN W, long j, GEN f, GEN pd, GEN p0, GEN p1, int noinv)
{
  pari_sp av = avma;
  long space = lg(f) * lgefint(p1);
  GEN a2, b2, g, z, s, t;
  GEN a = gel(V,j), b = gel(V,j+1);
  GEN u = gel(W,j), v = gel(W,j+1);

  (void)new_chunk(space); /* HACK */
  g = ZX_sub(f, ZX_mul(a,b));
  g = ZX_Z_divexact(g, p0);
  g = FpX_red(g, pd);
  z = FpX_mul(v,g, pd);
  t = FpX_divrem(z,a, pd, &s);
  t = ZX_add(ZX_mul(u,g), ZX_mul(t,b));
  t = FpX_red(t, pd);
  t = ZX_Z_mul(t,p0);
  s = ZX_Z_mul(s,p0);
  set_avma(av);
  a2 = ZX_add(a,s);
  b2 = ZX_add(b,t);

  /* already reduced mod p1 = pd p0 */
  gel(V,j)   = a2;
  gel(V,j+1) = b2;
  if (noinv) return;

  av = avma;
  (void)new_chunk(space); /* HACK */
  g = ZX_add(ZX_mul(u,a2), ZX_mul(v,b2));
  g = Z_ZX_sub(gen_1, g);
  g = ZX_Z_divexact(g, p0);
  g = FpX_red(g, pd);
  z = FpX_mul(v,g, pd);
  t = FpX_divrem(z,a, pd, &s);
  t = ZX_add(ZX_mul(u,g), ZX_mul(t,b));
  t = FpX_red(t, pd);
  t = ZX_Z_mul(t,p0);
  s = ZX_Z_mul(s,p0);
  set_avma(av);
  gel(W,j)   = ZX_add(u,t);
  gel(W,j+1) = ZX_add(v,s);
}

static void
ZpXQ_HenselLift(GEN V, GEN W, long j, GEN f, GEN Td, GEN T1, GEN pd, GEN p0, GEN p1, int noinv)
{
  pari_sp av = avma;
  const long n = degpol(T1), vT = varn(T1);
  long space = lg(f) * lgefint(p1) * lg(T1);
  GEN a2, b2, g, z, s, t;
  GEN a = gel(V,j), b = gel(V,j+1);
  GEN u = gel(W,j), v = gel(W,j+1);

  (void)new_chunk(space); /* HACK */
  g = RgX_sub(f, Kronecker_to_ZXX(ZXX_mul_Kronecker(a,b,n), n, vT));
  g = FpXQX_red(g, T1, p1);
  g = RgX_Rg_divexact(g, p0);
  z = FpXQX_mul(v,g, Td,pd);
  t = FpXQX_divrem(z,a, Td,pd, &s);
  t = ZX_add(ZXX_mul_Kronecker(u,g,n), ZXX_mul_Kronecker(t,b,n));
  t = Kronecker_to_ZXX(t, n, vT);
  t = FpXQX_red(t, Td, pd);
  t = RgX_Rg_mul(t,p0);
  s = RgX_Rg_mul(s,p0);
  set_avma(av);

  a2 = RgX_add(a,s);
  b2 = RgX_add(b,t);
  /* already reduced mod p1 = pd p0 */
  gel(V,j)   = a2;
  gel(V,j+1) = b2;
  if (noinv) return;

  av = avma;
  (void)new_chunk(space); /* HACK */
  g = ZX_add(ZXX_mul_Kronecker(u,a2,n), ZXX_mul_Kronecker(v,b2,n));
  g = Kronecker_to_ZXX(g, n, vT);
  g = Rg_RgX_sub(gen_1, g);
  g = FpXQX_red(g, T1, p1);
  g = RgX_Rg_divexact(g, p0);
  z = FpXQX_mul(v,g, Td,pd);
  t = FpXQX_divrem(z,a, Td,pd, &s);
  t = ZX_add(ZXX_mul_Kronecker(u,g,n), ZXX_mul_Kronecker(t,b,n));
  t = Kronecker_to_ZXX(t, n, vT);
  t = FpXQX_red(t, Td, pd);
  t = RgX_Rg_mul(t,p0);
  s = RgX_Rg_mul(s,p0);
  set_avma(av);
  gel(W,j)   = RgX_add(u,t);
  gel(W,j+1) = RgX_add(v,s);
}

/* v list of factors, w list of inverses.  f = v[j] v[j+1]
 * Lift v[j] and v[j+1] mod p0 pd (possibly mod T), then all their divisors */
static void
ZpX_RecTreeLift(GEN link, GEN v, GEN w, GEN pd, GEN p0, GEN p1,
                GEN f, long j, int noinv)
{
  if (j < 0) return;
  ZpX_HenselLift(v, w, j, f, pd, p0,p1, noinv);
  ZpX_RecTreeLift(link, v, w, pd, p0,p1, gel(v,j)  , link[j  ], noinv);
  ZpX_RecTreeLift(link, v, w, pd, p0,p1, gel(v,j+1), link[j+1], noinv);
}
static void
ZpXQ_RecTreeLift(GEN link, GEN v, GEN w, GEN Td, GEN T1, GEN pd, GEN p0, GEN p1,
                 GEN f, long j, int noinv)
{
  if (j < 0) return;
  ZpXQ_HenselLift(v, w, j, f, Td,T1, pd, p0,p1, noinv);
  ZpXQ_RecTreeLift(link, v, w, Td,T1, pd, p0,p1, gel(v,j)  , link[j  ], noinv);
  ZpXQ_RecTreeLift(link, v, w, Td,T1, pd, p0,p1, gel(v,j+1), link[j+1], noinv);
}

/* Lift to precision p^e0.
 * a = modular factors of f mod (p,T) [possibly T=NULL]
 *  OR a TreeLift structure [e, link, v, w]: go on lifting
 * flag = 0: standard.
 * flag = 1: return TreeLift structure */
static GEN
MultiLift(GEN f, GEN a, GEN T, GEN p, long e0, long flag)
{
  long i, eold, e, k = lg(a) - 1;
  GEN E, v, w, link, penew, Tnew;
  ulong mask;
  pari_timer Ti;

  if (k < 2) pari_err_DOMAIN("MultiLift", "#(modular factors)", "<", gen_2,a);
  if (e0 < 1) pari_err_DOMAIN("MultiLift", "precision", "<", gen_1,stoi(e0));
  if (e0 == 1) return a;

  if (DEBUGLEVEL > 3) timer_start(&Ti);
  if (typ(gel(a,1)) == t_INT)
  { /* a = TreeLift structure */
    e = itos(gel(a,1));
    link = gel(a,2);
    v    = gel(a,3);
    w    = gel(a,4);
  }
  else
  {
    e = 1;
    v = cgetg(2*k-2 + 1, t_VEC);
    w = cgetg(2*k-2 + 1, t_VEC);
    link=cgetg(2*k-2 + 1, t_VECSMALL);
    BuildTree(link, v, w, a, T? FpX_red(T,p): NULL, p);
    if (DEBUGLEVEL > 3) timer_printf(&Ti, "building tree");
  }
  mask = quadratic_prec_mask(e0);
  eold = 1;
  penew = NULL;
  Tnew = NULL;
  if (DEBUGLEVEL > 3) err_printf("lifting to prec %ld\n", e0);
  while (mask > 1)
  {
    long enew = eold << 1;
    if (mask & 1) enew--;
    mask >>= 1;
    if (enew >= e) { /* mask == 1: last iteration */
      GEN peold = penew? penew: powiu(p, eold);
      GEN Td = NULL, pd;
      long d = enew - eold; /* = eold or eold-1 */
      /* lift from p^eold to p^enew */
      pd = (d == eold)? peold: diviiexact(peold, p); /* p^d */
      penew = mulii(peold,pd);
      if (T) {
        if (Tnew)
          Td = (d == eold)? Tnew: FpX_red(Tnew,pd);
        else
          Td = FpX_red(T, peold);
        Tnew = FpX_red(T, penew);
        ZpXQ_RecTreeLift(link, v, w, Td, Tnew, pd, peold, penew, f, lg(v)-2,
                         (flag == 0 && mask == 1));
      }
      else
        ZpX_RecTreeLift(link, v, w, pd, peold, penew, f, lg(v)-2,
                        (flag == 0 && mask == 1));
      if (DEBUGLEVEL > 3) timer_printf(&Ti, "reaching prec %ld", enew);
    }
    eold = enew;
  }

  if (flag)
    E = mkvec4(utoipos(e0), link, v, w);
  else
  {
    E = cgetg(k+1, t_VEC);
    for (i = 1; i <= 2*k-2; i++)
    {
      long t = link[i];
      if (t < 0) gel(E,-t) = gel(v,i);
    }
  }
  return E;
}

/* Q list of (coprime, monic) factors of pol mod (T,p). Lift mod p^e = pe.
 * T may be NULL */
GEN
ZpX_liftfact(GEN pol, GEN Q, GEN pe, GEN p, long e)
{
  pari_sp av = avma;
  pol = FpX_normalize(pol, pe);
  if (lg(Q) == 2) return mkvec(pol);
  return gerepilecopy(av, MultiLift(pol, Q, NULL, p, e, 0));
}

GEN
ZpXQX_liftfact(GEN pol, GEN Q, GEN T, GEN pe, GEN p, long e)
{
  pari_sp av = avma;
  pol = FpXQX_normalize(pol, T, pe);
  if (lg(Q) == 2) return mkvec(pol);
  return gerepilecopy(av, MultiLift(pol, Q, T, p, e, 0));
}

GEN
ZqX_liftfact(GEN f, GEN a, GEN T, GEN pe, GEN p, long e)
{ return T ? ZpXQX_liftfact(f, a, T, pe, p, e): ZpX_liftfact(f, a, pe, p, e); }
GEN
ZqX_liftroot(GEN f, GEN a, GEN T, GEN p, long e)
{ return T ? ZpXQX_liftroot(f, a,T , p, e): ZpX_liftroot(f, a, p, e); }

/* U = NULL treated as 1 */
static void
BezoutPropagate(GEN link, GEN v, GEN w, GEN pe, GEN U, GEN f, long j)
{
  GEN Q, R;
  if (j < 0) return;

  Q = FpX_mul(gel(v,j), gel(w,j), pe);
  if (U)
  {
    Q = FpXQ_mul(Q, U, f, pe);
    R = FpX_sub(U, Q, pe);
  }
  else
    R = Fp_FpX_sub(gen_1, Q, pe);
  gel(w,j+1) = Q; /*  0 mod U v[j],  1 mod (1-U) v[j+1] */
  gel(w,j) = R; /*  1 mod U v[j],  0 mod (1-U) v[j+1] */
  BezoutPropagate(link, v, w, pe, R, f, link[j  ]);
  BezoutPropagate(link, v, w, pe, Q, f, link[j+1]);
}

/* as above, but return the Bezout coefficients for the lifted modular factors
 *   U[i] = 1 mod Qlift[i]
 *          0 mod Qlift[j], j != i */
GEN
bezout_lift_fact(GEN pol, GEN Q, GEN p, long e)
{
  pari_sp av = avma;
  GEN E, link, v, w, pe;
  long i, k = lg(Q)-1;
  if (k == 1) retmkvec(pol_1(varn(pol)));
  pe = powiu(p, e);
  pol = FpX_normalize(pol, pe);
  E = MultiLift(pol, Q, NULL, p, e, 1);
  link = gel(E,2);
  v    = gel(E,3);
  w    = gel(E,4);
  BezoutPropagate(link, v, w, pe, NULL, pol, lg(v)-2);
  E = cgetg(k+1, t_VEC);
  for (i = 1; i <= 2*k-2; i++)
  {
    long t = link[i];
    if (t < 0) E[-t] = w[i];
  }
  return gerepilecopy(av, E);
}

/* Front-end for ZpX_liftfact:
   lift the factorization of pol mod p given by L to p^N (if possible) */
GEN
polhensellift(GEN pol, GEN L, GEN Tp, long N)
{
  GEN T, p;
  long i, l;
  pari_sp av = avma;
  void (*chk)(GEN, const char*);

  if (typ(pol) != t_POL) pari_err_TYPE("polhensellift",pol);
  RgX_check_ZXX(pol, "polhensellift");
  if (!is_vec_t(typ(L)) || lg(L) < 3) pari_err_TYPE("polhensellift",L);
  if (N < 1) pari_err_DOMAIN("polhensellift", "precision", "<", gen_1,stoi(N));
  if (!ff_parse_Tp(Tp, &T, &p, 0)) pari_err_TYPE("polhensellift",Tp);
  chk = T? RgX_check_ZXX: RgX_check_ZX;
  l = lg(L); L = leafcopy(L);
  for (i = 1; i < l; i++)
  {
    GEN q = gel(L,i);
    if (typ(q) != t_POL) gel(L,i) = scalar_ZX_shallow(q, varn(pol));
    else chk(q, "polhensellift");
  }
  return gerepilecopy(av, ZqX_liftfact(pol, L, T, powiu(p,N), p, N));
}

static GEN
FqV_roots_from_deg1(GEN x, GEN T, GEN p)
{
  long i,l = lg(x);
  GEN r = cgetg(l,t_COL);
  for (i=1; i<l; i++) { GEN P = gel(x,i); gel(r,i) = Fq_neg(gel(P,2), T, p); }
  return r;
}

static GEN
ZpX_liftroots_full(GEN f, GEN S, GEN q, GEN p, long e)
{
  pari_sp av = avma;
  GEN y = ZpX_liftfact(f, deg1_from_roots(S, varn(f)), q, p, e);
  return gerepileupto(av, FqV_roots_from_deg1(y, NULL, q));
}

GEN
ZpX_roots(GEN F, GEN p, long e)
{
  pari_sp av = avma;
  GEN q = powiu(p, e);
  GEN f = FpX_normalize(F, p);
  GEN g = FpX_normalize(FpX_split_part(f, p), p);
  GEN S;
  if (degpol(g) < degpol(f))
  {
    GEN h = FpX_div(f, g, p);
    F = gel(ZpX_liftfact(F, mkvec2(g, h), q, p, e), 1);
  }
  S = FpX_roots(g, p);
  return gerepileupto(av, ZpX_liftroots_full(F, S, q, p, e));
}

static GEN
ZpXQX_liftroots_full(GEN f, GEN S, GEN T, GEN q, GEN p, long e)
{
  pari_sp av = avma;
  GEN y = ZpXQX_liftfact(f, deg1_from_roots(S, varn(f)), T, q, p, e);
  return gerepileupto(av, FqV_roots_from_deg1(y, T, q));
}

GEN
ZpXQX_roots(GEN F, GEN T, GEN p, long e)
{
  pari_sp av = avma;
  GEN q = powiu(p, e);
  GEN f = FpXQX_normalize(F, T, p);
  GEN g = FpXQX_normalize(FpXQX_split_part(f, T, p), T, p);
  GEN S;
  if (degpol(g) < degpol(f))
  {
    GEN h = FpXQX_div(f, g, T, p);
    F = gel(ZpXQX_liftfact(F, mkvec2(g, h), T, q, p, e), 1);
  }
  S = FpXQX_roots(g, T, p);
  return gerepileupto(av, ZpXQX_liftroots_full(F, S, T, q, p, e));
}

GEN
ZqX_roots(GEN F, GEN T, GEN p, long e)
{
  return T ? ZpXQX_roots(F, T, p, e): ZpX_roots(F, p, e);
}

/* SPEC:
 * p is a t_INT > 1, e >= 1
 * f is a ZX with leading term prime to p.
 * a is a simple root mod l for all l|p.
 * Return roots of f mod p^e, as integers (implicitly mod p^e)
 * STANDARD USE: p is a prime power */
GEN
ZpX_liftroot(GEN f, GEN a, GEN p, long e)
{
  pari_sp av = avma;
  GEN q = p, fr, W;
  ulong mask;

  a = modii(a,q);
  if (e == 1) return a;
  mask = quadratic_prec_mask(e);
  fr = FpX_red(f,q);
  W = Fp_inv(FpX_eval(ZX_deriv(fr), a, q), q); /* 1/f'(a) mod p */
  for(;;)
  {
    q = sqri(q);
    if (mask & 1) q = diviiexact(q, p);
    mask >>= 1;
    fr = FpX_red(f,q);
    a = Fp_sub(a, Fp_mul(W, FpX_eval(fr, a,q), q), q);
    if (mask == 1) return gerepileuptoint(av, a);
    W = Fp_sub(shifti(W,1), Fp_mul(Fp_sqr(W,q), FpX_eval(ZX_deriv(fr),a,q), q), q);
  }
}

GEN
ZpX_liftroots(GEN f, GEN S, GEN p, long e)
{
  long i, n = lg(S)-1, d = degpol(f);
  GEN r;
  if (n == d) return ZpX_liftroots_full(f, S, powiu(p, e), p, e);
  r = cgetg(n+1, typ(S));
  for (i=1; i <= n; i++)
    gel(r,i) = ZpX_liftroot(f, gel(S,i), p, e);
  return r;
}

GEN
ZpXQX_liftroot_vald(GEN f, GEN a, long v, GEN T, GEN p, long e)
{
  pari_sp av = avma, av2;
  GEN pv = p, q, W, df, Tq, fr, dfr;
  ulong mask;
  a = Fq_red(a, T, p);
  if (e <= v+1) return a;
  df = RgX_deriv(f);
  if (v) { pv = powiu(p,v); df = ZXX_Z_divexact(df, pv); }
  mask = quadratic_prec_mask(e-v);
  Tq = FpXT_red(T, p); dfr = FpXQX_red(df, Tq, p);
  W = Fq_inv(FqX_eval(dfr, a, Tq, p), Tq, p); /* 1/f'(a) mod (T,p) */
  q = p;
  av2 = avma;
  for (;;)
  {
    GEN u, fa, qv, q2v, q2, Tq2;
    q2 = q; q = sqri(q);
    if (mask & 1) q = diviiexact(q,p);
    mask >>= 1;
    if (v) { qv = mulii(q, pv); q2v = mulii(q2, pv); }
    else { qv = q; q2v = q2; }
    Tq2 = FpXT_red(T, q2v); Tq = FpXT_red(T, qv);
    fr = FpXQX_red(f, Tq, qv);
    fa = FqX_eval(fr, a, Tq, qv);
    fa = typ(fa)==t_INT? diviiexact(fa,q2v): ZX_Z_divexact(fa, q2v);
    a = Fq_sub(a, Fq_mul(Fq_mul(W,fa,Tq2,q2v),q2, Tq,qv), Tq, qv);
    if (mask == 1) return gerepileupto(av, a);
    dfr = FpXQX_red(df, Tq, q);
    u = Fq_sub(Fq_mul(W,FqX_eval(dfr,a,Tq,q),Tq,q),gen_1,Tq,q);
    u = typ(u)==t_INT? diviiexact(u,q2): ZX_Z_divexact(u,q2);
    W = Fq_sub(W, Fq_mul(Fq_mul(u,W,Tq2,q2),q2, Tq,q),Tq,q);
    if (gc_needed(av2,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZpXQX_liftroot, e = %ld", e);
      gerepileall(av2, 3, &a, &W, &q);
    }
  }
}

GEN
ZpXQX_liftroot(GEN f, GEN a, GEN T, GEN p, long e) { return ZpXQX_liftroot_vald(f,a,0,T,p,e); }

GEN
ZpXQX_liftroots(GEN f, GEN S, GEN T, GEN p, long e)
{
  long i, n = lg(S)-1, d = degpol(f);
  GEN r;
  if (n == d) return ZpXQX_liftroots_full(f, S, T, powiu(p, e), p, e);
  r = cgetg(n+1, typ(S));
  for (i=1; i <= n; i++)
    gel(r,i) = ZpXQX_liftroot(f, gel(S,i), T, p, e);
  return r;
}

/* Compute (x-1)/(x+1)/p^k */
static GEN
ZpXQ_log_to_ath(GEN x, long k, GEN T, GEN p, long e, GEN pe)
{
  pari_sp av = avma;
  long vT = get_FpX_var(T);
  GEN bn, bdi;
  GEN bd = ZX_Z_add(x, gen_1);
  if (absequaliu(p,2)) /*For p=2, we need to simplify by 2*/
  {
    bn = ZX_shifti(x,-(k+1));
    bdi= ZpXQ_invlift(ZX_shifti(bd ,-1), pol_1(vT), T, p, e);
  }
  else
  {
    bn = ZX_Z_divexact(ZX_Z_sub(x, gen_1),powiu(p,k));
    bdi= ZpXQ_invlift(bd, scalarpol(Fp_inv(gen_2,p),vT), T, p, e);
  }
  return gerepileupto(av, FpXQ_mul(bn, bdi, T, pe));
}

/* Assume p odd, a = 1 [p], return log(a) */
GEN
ZpXQ_log(GEN a, GEN T, GEN p, long N)
{
  pari_sp av = avma;
  pari_timer ti;
  long is2 = absequaliu(p,2);
  ulong pp = is2 ? 0: itou_or_0(p);
  double lp = is2 ? 1: pp ? log2(pp): expi(p);
  long k = maxss(1 , (long) .5+pow((double)(N>>1)/(lp*lp), 1./3));
  GEN ak, s, b, pol;
  long e = is2 ? N-1: N;
  long i, l = (e-2)/(2*(k+is2));
  GEN pe = powiu(p,e);
  GEN TNk, pNk = powiu(p,N+k);
  if( DEBUGLEVEL>=3) timer_start(&ti);
  TNk = FpX_get_red(get_FpX_mod(T), pNk);
  ak = FpXQ_pow(a, powiu(p,k), TNk, pNk);
  if( DEBUGLEVEL>=3) timer_printf(&ti,"FpXQ_pow(%ld)",k);
  b = ZpXQ_log_to_ath(ak, k, T, p, e, pe);
  if( DEBUGLEVEL>=3) timer_printf(&ti,"ZpXQ_log_to_ath");
  pol= cgetg(l+3,t_POL);
  pol[1] = evalsigne(1)|evalvarn(0);
  for(i=0; i<=l; i++)
  {
    GEN g;
    ulong z = 2*i+1;
    if (pp)
    {
      long w = u_lvalrem(z, pp, &z);
      g = powuu(pp,2*i*k-w);
    }
    else g = powiu(p,2*i*k);
    gel(pol,i+2) = Fp_divu(g, z,pe);
  }
  if( DEBUGLEVEL>=3) timer_printf(&ti,"pol(%ld)",l);
  s = FpX_FpXQ_eval(pol, FpXQ_sqr(b, T, pe), T,  pe);
  if( DEBUGLEVEL>=3) timer_printf(&ti,"FpX_FpXQ_eval");
  s = ZX_shifti(FpXQ_mul(b, s, T, pe), 1);
  return gerepileupto(av, is2? s: FpX_red(s, pe));
}

/***********************************************************************/
/**                                                                   **/
/**                 Generic quadratic hensel lift over Zp[X]          **/
/**                                                                   **/
/***********************************************************************/
/* q = p^N */
GEN
gen_ZpM_Dixon(GEN F, GEN V, GEN q, GEN p, long N, void *E,
                            GEN lin(void *E, GEN F, GEN d, GEN q),
                            GEN invl(void *E, GEN d))
{
  pari_sp av = avma;
  long N2, M;
  GEN VN2, V2, VM, bil;
  GEN q2, qM;
  V = FpM_red(V, q);
  if (N == 1) return invl(E, V);
  N2 = (N + 1)>>1; M = N - N2;
  F = FpM_red(F, q);
  qM = powiu(p, M);
  q2 = M == N2? qM: mulii(qM, p);
  /* q2 = p^N2, qM = p^M, q = q2 * qM */
  VN2 = gen_ZpM_Dixon(F, V, q2, p, N2, E, lin, invl);
  bil = lin(E, F, VN2, q);
  V2 = ZM_Z_divexact(ZM_sub(V, bil), q2);
  VM = gen_ZpM_Dixon(F, V2, qM, p, M, E, lin, invl);
  return gerepileupto(av, FpM_red(ZM_add(VN2, ZM_Z_mul(VM, q2)), q));
}

GEN
gen_ZpX_Dixon(GEN F, GEN V, GEN q, GEN p, long N, void *E,
                            GEN lin(void *E, GEN F, GEN d, GEN q),
                            GEN invl(void *E, GEN d))
{
  pari_sp av = avma;
  long N2, M;
  GEN VN2, V2, VM, bil;
  GEN q2, qM;
  V = FpX_red(V, q);
  if (N == 1) return invl(E, V);
  N2 = (N + 1)>>1; M = N - N2;
  F = FpXT_red(F, q);
  qM = powiu(p, M);
  q2 = M == N2? qM: mulii(qM, p);
  /* q2 = p^N2, qM = p^M, q = q2 * qM */
  VN2 = gen_ZpX_Dixon(F, V, q2, p, N2, E, lin, invl);
  bil = lin(E, F, VN2, q);
  V2 = ZX_Z_divexact(ZX_sub(V, bil), q2);
  VM = gen_ZpX_Dixon(F, V2, qM, p, M, E, lin, invl);
  return gerepileupto(av, FpX_red(ZX_add(VN2, ZX_Z_mul(VM, q2)), q));
}

GEN
gen_ZpM_Newton(GEN x, GEN p, long n, void *E,
                      GEN eval(void *E, GEN f, GEN q),
                      GEN invd(void *E, GEN V, GEN v, GEN q, long M))
{
  pari_sp ltop = avma, av;
  long N = 1, N2, M;
  long mask;
  GEN q = p;
  if (n == 1) return gcopy(x);
  mask = quadratic_prec_mask(n);
  av = avma;
  while (mask > 1)
  {
    GEN qM, q2, v, V;
    N2 = N; N <<= 1;
    q2 = q;
    if (mask&1UL) { /* can never happen when q2 = p */
      N--; M = N2-1;
      qM = diviiexact(q2,p); /* > 1 */
      q = mulii(qM,q2);
    } else {
      M = N2;
      qM = q2;
      q = sqri(q2);
    }
    /* q2 = p^N2, qM = p^M, q = p^N = q2 * qM */
    mask >>= 1;
    v = eval(E, x, q);
    V = ZM_Z_divexact(gel(v,1), q2);
    x = FpM_sub(x, ZM_Z_mul(invd(E, V, v, qM, M), q2), q);
    if (gc_needed(av, 1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"gen_ZpM_Newton");
      gerepileall(av, 2, &x, &q);
    }
  }
  return gerepileupto(ltop, x);
}

static GEN
_ZpM_invd(void *E, GEN V, GEN v, GEN q, long M/*unused*/)
{
  (void)E; (void)M;
  return FpM_mul(V, gel(v,2), q);
}

static GEN
_ZpM_eval(void *E, GEN x, GEN q)
{
  GEN f = RgM_Rg_sub_shallow(FpM_mul(x, FpM_red((GEN) E, q), q), gen_1);
  return mkvec2(f, x);
}

GEN
ZpM_invlift(GEN M, GEN C, GEN p, long n)
{
  return gen_ZpM_Newton(C, p, n, (void *)M, _ZpM_eval, _ZpM_invd);
}

GEN
gen_ZpX_Newton(GEN x, GEN p, long n, void *E,
                      GEN eval(void *E, GEN f, GEN q),
                      GEN invd(void *E, GEN V, GEN v, GEN q, long M))
{
  pari_sp ltop = avma, av;
  long N = 1, N2, M;
  long mask;
  GEN q = p;
  if (n == 1) return gcopy(x);
  mask = quadratic_prec_mask(n);
  av = avma;
  while (mask > 1)
  {
    GEN qM, q2, v, V;
    N2 = N; N <<= 1;
    q2 = q;
    if (mask&1UL) { /* can never happen when q2 = p */
      N--; M = N2-1;
      qM = diviiexact(q2,p); /* > 1 */
      q = mulii(qM,q2);
    } else {
      M = N2;
      qM = q2;
      q = sqri(q2);
    }
    /* q2 = p^N2, qM = p^M, q = p^N = q2 * qM */
    mask >>= 1;
    v = eval(E, x, q);
    V = ZX_Z_divexact(gel(v,1), q2);
    x = FpX_sub(x, ZX_Z_mul(invd(E, V, v, qM, M), q2), q);
    if (gc_needed(av, 1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"gen_ZpX_Newton");
      gerepileall(av, 2, &x, &q);
    }
  }
  return gerepileupto(ltop, x);
}

struct _ZpXQ_inv
{
  GEN T, a, p ,n;
};

static GEN
_inv_invd(void *E, GEN V, GEN v, GEN q, long M/*unused*/)
{
  struct _ZpXQ_inv *d = (struct _ZpXQ_inv *) E;
  GEN Tq = FpXT_red(d->T, q);
  (void)M;
  return FpXQ_mul(V, gel(v,2), Tq, q);
}

static GEN
_inv_eval(void *E, GEN x, GEN q)
{
  struct _ZpXQ_inv *d = (struct _ZpXQ_inv *) E;
  GEN Tq = FpXT_red(d->T, q);
  GEN f = FpX_Fp_sub(FpXQ_mul(x, FpX_red(d->a, q), Tq, q), gen_1, q);
  return mkvec2(f, x);
}

GEN
ZpXQ_invlift(GEN a, GEN x, GEN T, GEN p, long e)
{
  struct _ZpXQ_inv d;
  d.a = a; d.T = T; d.p = p;
  return gen_ZpX_Newton(x, p, e, &d, _inv_eval, _inv_invd);
}

GEN
ZpXQ_inv(GEN a, GEN T, GEN p, long e)
{
  pari_sp av=avma;
  GEN ai;
  if (lgefint(p)==3)
  {
    ulong pp = p[2];
    ai = Flx_to_ZX(Flxq_inv(ZX_to_Flx(a,pp), ZXT_to_FlxT(T, pp), pp));
  } else
    ai = FpXQ_inv(FpX_red(a,p), FpXT_red(T,p),p);
  return gerepileupto(av, ZpXQ_invlift(a, ai, T, p, e));
}

GEN
ZpXQ_div(GEN a, GEN b, GEN T, GEN q, GEN p, long e)
{
  return FpXQ_mul(a, ZpXQ_inv(b, T, p, e), T, q);
}

GEN
ZpXQX_divrem(GEN x, GEN Sp, GEN T, GEN q, GEN p, long e, GEN *pr)
{
  pari_sp av = avma;
  GEN S = get_FpXQX_mod(Sp);
  GEN b = leading_coeff(S), bi;
  GEN S2, Q;
  if (typ(b)==t_INT) return FpXQX_divrem(x, Sp, T, q, pr);
  bi = ZpXQ_inv(b, T, p, e);
  S2 = FqX_Fq_mul_to_monic(S, bi, T, q);
  Q = FpXQX_divrem(x, S2, T, q, pr);
  if (pr==ONLY_DIVIDES && !Q) { set_avma(av); return NULL; }
  if (pr==ONLY_REM || pr==ONLY_DIVIDES) return gerepileupto(av, Q);
  Q = FpXQX_FpXQ_mul(Q, bi, T, q);
  return gc_all(av, 2, &Q, pr);
}

GEN
ZpXQX_digits(GEN x, GEN B, GEN T, GEN q, GEN p, long e)
{
  pari_sp av = avma;
  GEN b = leading_coeff(B), bi;
  GEN B2, P, V, W;
  long i, lV;
  if (typ(b)==t_INT) return FpXQX_digits(x, B, T, q);
  bi = ZpXQ_inv(b, T, p, e);
  B2 = FqX_Fq_mul_to_monic(B, bi, T, q);
  V = FpXQX_digits(x, B2, T, q);
  lV = lg(V)-1;
  P = FpXQ_powers(bi, lV-1, T, q);
  W = cgetg(lV+1, t_VEC);
  for(i=1; i<=lV; i++)
    gel(W, i) = FpXQX_FpXQ_mul(gel(V,i), gel(P, i), T, q);
  return gerepileupto(av, W);
}

struct _ZpXQ_sqrtn
{
  GEN T, a, n, ai;
};

static GEN
_sqrtn_invd(void *E, GEN V, GEN v, GEN q, long M)
{
  struct _ZpXQ_sqrtn *d = (struct _ZpXQ_sqrtn *) E;
  GEN Tq = FpX_red(d->T, q), aiq = FpX_red(d->ai, q);
  (void)M;
  return FpXQ_mul(FpXQ_mul(V, gel(v,2), Tq, q), aiq, Tq, q);
}

static GEN
_sqrtn_eval(void *E, GEN x, GEN q)
{
  struct _ZpXQ_sqrtn *d = (struct _ZpXQ_sqrtn *) E;
  GEN Tq = FpX_red(d->T, q);
  GEN f = FpX_sub(FpXQ_pow(x, d->n, Tq, q), d->a, q);
  return mkvec2(f, x);
}

GEN
ZpXQ_sqrtnlift(GEN a, GEN n, GEN x, GEN T, GEN p, long e)
{
  struct _ZpXQ_sqrtn d;
  d.a = a; d.T = T; d.n = n;
  d.ai = ZpXQ_inv(ZX_Z_mul(a, n),T,p,(e+1)>>1);
  return gen_ZpX_Newton(x, p, e, &d, _sqrtn_eval, _sqrtn_invd);
}

static GEN
to_ZX(GEN a, long v) { return typ(a)==t_INT? scalarpol_shallow(a,v): a; }

GEN
Zq_sqrtnlift(GEN a, GEN n, GEN x, GEN T, GEN p, long e)
{
  return T? ZpXQ_sqrtnlift(to_ZX(a,varn(T)), n, to_ZX(x,varn(T)), T, p, e)
          : Zp_sqrtnlift(a, n, x, p, e);
}

GEN
ZpXQ_sqrt(GEN a, GEN T, GEN p, long e)
{
  pari_sp av = avma;
  GEN z = FpXQ_sqrt(FpX_red(a, p), T, p);
  if (!z) return NULL;
  if (e <= 1) return gerepileupto(av, z);
  return gerepileupto(av, ZpXQ_sqrtnlift(a, gen_2, z, T, p, e));
}

GEN
ZpX_ZpXQ_liftroot_ea(GEN P, GEN S, GEN T, GEN p, long n, void *E,
                     GEN early(void *E, GEN x, GEN q))
{
  pari_sp ltop = avma, av;
  long N, r;
  long mask;
  GEN q2, q, W, Q, Tq2, Tq, Pq;
  pari_timer ti;
  T = FpX_get_red(T, powiu(p, n));
  if (n == 1) return gcopy(S);
  mask = quadratic_prec_mask(n);
  av = avma;
  q2 = p; q = sqri(p); mask >>= 1; N = 2;
  if (DEBUGLEVEL > 3) timer_start(&ti);
  Tq = FpXT_red(T,q);
  Tq2 = FpXT_red(Tq,q2);
  Pq = FpX_red(P,q);
  W = FpXQ_inv(FpX_FpXQ_eval(FpX_deriv(P,q2), S, Tq2, q2), Tq2, q2);
  Q  = ZX_Z_divexact(FpX_FpXQ_eval(Pq, S, Tq, q), q2);
  r = brent_kung_optpow(degpol(P), 4, 3);
  if (DEBUGLEVEL > 3)
    err_printf("ZpX_ZpXQ_liftroot: lifting to prec %ld\n",n);
  for (;;)
  {
    GEN H, Sq, Wq, Spow, dP, qq, Pqq, Tqq;
    H  = FpXQ_mul(W, Q, Tq2, q2);
    Sq = FpX_sub(S, ZX_Z_mul(H, q2), q);
    if (DEBUGLEVEL > 3)
      timer_printf(&ti,"ZpX_ZpXQ_liftroot: reaching prec %ld",N);
    if (mask==1)
      return gerepileupto(ltop, Sq);
    if (early)
    {
      GEN Se = early(E, Sq, q);
      if (Se) return gerepileupto(ltop, Se);
    }
    qq = sqri(q); N <<= 1;
    if (mask&1UL) { qq = diviiexact(qq, p); N--; }
    mask >>= 1;
    Pqq  = FpX_red(P, qq);
    Tqq  = FpXT_red(T, qq);
    Spow = FpXQ_powers(Sq, r, Tqq, qq);
    Q  = ZX_Z_divexact(FpX_FpXQV_eval(Pqq, Spow, Tqq, qq), q);
    dP = FpX_FpXQV_eval(FpX_deriv(Pq, q), FpXV_red(Spow, q), Tq, q);
    Wq = ZX_Z_divexact(FpX_Fp_sub(FpXQ_mul(W, dP, Tq, q), gen_1, q), q2);
    Wq = ZX_Z_mul(FpXQ_mul(W, Wq, Tq2, q2), q2);
    Wq = FpX_sub(W, Wq, q);
    S = Sq; W = Wq; q2 = q; q = qq; Tq2 = Tq; Tq = Tqq; Pq = Pqq;
    if (gc_needed(av, 1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZpX_ZpXQ_liftroot");
      gerepileall(av, 8, &S, &W, &Q, &Tq2, &Tq, &Pq, &q, &q2);
    }
  }
}

GEN
ZpX_ZpXQ_liftroot(GEN P, GEN S, GEN T, GEN p, long n)
{
  return ZpX_ZpXQ_liftroot_ea(P, S, T, p, n, NULL, NULL);
}

GEN
ZpX_Frobenius(GEN T, GEN p, long e)
{
  return ZpX_ZpXQ_liftroot(get_FpX_mod(T), FpX_Frobenius(T, p), T, p, e);
}

GEN
ZpXQM_prodFrobenius(GEN M, GEN T, GEN p, long e)
{
  pari_sp av = avma;
  GEN xp = ZpX_Frobenius(T, p, e);
  GEN z = FpXQM_autsum(mkvec2(xp, M), get_FpX_degree(T), T, powiu(p,e));
  return gerepilecopy(av, gel(z,2));
}

GEN
ZpXQX_ZpXQXQ_liftroot(GEN P, GEN S, GEN U, GEN T, GEN p, long n)
{
  pari_sp ltop = avma, av;
  long N, r;
  long mask;
  GEN  qn, q2, q, W, Q, Tq2, Tq, Pq, Uq, Uq2;
  pari_timer ti;
  qn = powiu(p, n);
  T = FpX_get_red(T, qn);
  U = FpXQX_get_red(U, T, qn);
  if (n == 1) return gcopy(S);
  mask = quadratic_prec_mask(n);
  av = avma;
  q2 = p; q = sqri(p); mask >>= 1; N = 2;
  if (DEBUGLEVEL > 3) timer_start(&ti);
  Tq = FpXT_red(T,q);
  Uq = FpXQXT_red(U, Tq, q);
  Tq2 = FpXT_red(Tq, q2);
  Uq2 = FpXQXT_red(U, Tq2, q2);
  Pq = FpXQX_red(P, Tq, q);
  W = FpXQXQ_inv(FpXQX_FpXQXQ_eval(FpXX_deriv(P,q2), S, Uq2, Tq2, q2), Uq2, Tq2, q2);
  Q  = ZXX_Z_divexact(FpXQX_FpXQXQ_eval(Pq, S, Uq, Tq, q), q2);
  r = brent_kung_optpow(degpol(P), 4, 3);
  if (DEBUGLEVEL > 3)
    err_printf("ZpXQX_ZpXQXQ_liftroot: lifting to prec %ld\n",n);
  for (;;)
  {
    GEN H, Sq, Wq, Spow, dP, qq, Pqq, Tqq, Uqq;
    H  = FpXQXQ_mul(W, Q, Uq2, Tq2, q2);
    Sq = FpXX_sub(S, ZXX_Z_mul(H, q2), q);
    if (DEBUGLEVEL > 3)
      timer_printf(&ti,"ZpXQX_ZpXQXQ_liftroot: reaching prec %ld",N);
    if (mask==1)
      return gerepileupto(ltop, Sq);
    qq = sqri(q); N <<= 1;
    if (mask&1UL) { qq = diviiexact(qq, p); N--; }
    mask >>= 1;
    Tqq  = FpXT_red(T, qq);
    Uqq  = FpXQXT_red(U, Tqq, qq);
    Pqq  = FpXQX_red(P, Tqq, qq);
    Spow = FpXQXQ_powers(Sq, r, Uqq, Tqq, qq);
    Q  = ZXX_Z_divexact(FpXQX_FpXQXQV_eval(Pqq, Spow, Uqq, Tqq, qq), q);
    dP = FpXQX_FpXQXQV_eval(FpXX_deriv(Pq, q), FpXQXV_red(Spow, Tq, q), Uq, Tq, q);
    Wq = ZXX_Z_divexact(gsub(FpXQXQ_mul(W, dP, Uq, Tq, q), gen_1), q2);
    Wq = ZXX_Z_mul(FpXQXQ_mul(W, Wq, Uq2, Tq2, q2), q2);
    Wq = FpXX_sub(W, Wq, q);
    S = Sq; W = Wq; q2 = q; q = qq; Tq2 = Tq; Tq = Tqq; Uq2 = Uq; Uq = Uqq;  Pq = Pqq;
    if (gc_needed(av, 1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ZpXQX_ZpXQXQ_liftroot");
      gerepileall(av, 10, &S, &W, &Q, &Uq2, &Uq, &Tq2, &Tq, &Pq, &q, &q2);
    }
  }
}

GEN
ZqX_ZqXQ_liftroot(GEN f, GEN a, GEN P, GEN T, GEN p, long e)
{ return T ? ZpXQX_ZpXQXQ_liftroot(f, a, P, T , p, e): ZpX_ZpXQ_liftroot(f, a, P, p, e); }

/* Canonical lift of polynomial */

static GEN _can_invl(void *E, GEN V) {(void) E; return V; }

static GEN _can_lin(void *E, GEN F, GEN V, GEN q)
{
  GEN v = RgX_splitting(V, 3);
  (void) E;
  return FpX_sub(V,ZXV_dotproduct(v, F), q);
}

static GEN
_can_iter(void *E, GEN f, GEN q)
{
  GEN h = RgX_splitting(f,3);
  GEN h1s = ZX_sqr(gel(h,1)), h2s = ZX_sqr(gel(h,2)), h3s = ZX_sqr(gel(h,3));
  GEN h12 = ZX_mul(gel(h,1), gel(h,2));
  GEN h13 = ZX_mul(gel(h,1), gel(h,3));
  GEN h23 = ZX_mul(gel(h,2), gel(h,3));
  GEN h1c = ZX_mul(gel(h,1), h1s);
  GEN h3c = ZX_mul(gel(h,3), h3s);
  GEN th = ZX_mul(ZX_sub(h2s,ZX_mulu(h13,3)),gel(h,2));
  GEN y = FpX_sub(f,ZX_add(RgX_shift_shallow(h3c,2),ZX_add(RgX_shift_shallow(th,1),h1c)),q);
  (void) E;
  return mkvecn(7,y,h1s,h2s,h3s,h12,h13,h23);
}

static GEN
_can_invd(void *E, GEN V, GEN v, GEN qM, long M)
{
  GEN h1s=gel(v,2), h2s=gel(v,3), h3s=gel(v,4);
  GEN h12=gel(v,5), h13=gel(v,6), h23=gel(v,7);
  GEN F = mkvec3(ZX_sub(h1s,RgX_shift_shallow(h23,1)),RgX_shift_shallow(ZX_sub(h2s,h13),1),
                 ZX_sub(RgX_shift_shallow(h3s,2),RgX_shift_shallow(h12,1)));
  (void)E;
  return gen_ZpX_Dixon(ZXV_Z_mul(F, utoi(3)), V, qM, utoi(3), M, NULL,
                                                 _can_lin, _can_invl);
}

static GEN
F3x_frobeniuslift(GEN P, long n)
{ return gen_ZpX_Newton(Flx_to_ZX(P),utoi(3), n, NULL, _can_iter, _can_invd); }

static GEN _can5_invl(void *E, GEN V) {(void) E; return V; }

static GEN _can5_lin(void *E, GEN F, GEN V, GEN q)
{
  ulong p = *(ulong*)E;
  GEN v = RgX_splitting(V, p);
  return FpX_sub(V,ZXV_dotproduct(v, F), q);
}

/* P(X,t) -> P(X*t^n,t) mod (t^p-1) */
static GEN
_shift(GEN P, long n, ulong p, long v)
{
  long i, l=lg(P);
  GEN r = cgetg(l,t_POL); r[1] = P[1];
  for(i=2;i<l;i++)
  {
    long s = n*(i-2)%p;
    GEN ci = gel(P,i);
    if (typ(ci)==t_INT)
      gel(r,i) = monomial(ci, s, v);
    else
      gel(r,i) = RgX_rotate_shallow(ci, s, p);
  }
  return FpXX_renormalize(r, l);
}

struct _can_mul
{
  GEN T, q;
  ulong p;
};

static GEN
_can5_mul(void *E, GEN A, GEN B)
{
  struct _can_mul *d = (struct _can_mul *)E;
  GEN a = gel(A,1), b = gel(B,1);
  long n = itos(gel(A,2));
  GEN bn = _shift(b, n, d->p, get_FpX_var(d->T));
  GEN c = FpXQX_mul(a, bn, d->T, d->q);
  return mkvec2(c, addii(gel(A,2), gel(B,2)));
}

static GEN
_can5_sqr(void *E, GEN A)
{
  return _can5_mul(E,A,A);
}

static GEN
_can5_iter(void *E, GEN f, GEN q)
{
  pari_sp av = avma;
  struct _can_mul D;
  ulong p = *(ulong*)E;
  long i, vT = fetch_var();
  GEN N, P, d, V, fs;
  D.q = q; D.T = ZX_Z_sub(pol_xn(p,vT),gen_1);
  D.p = p;
  fs = mkvec2(_shift(f, 1, p, vT), gen_1);
  N = gel(gen_powu_i(fs,p-1,(void*)&D,_can5_sqr,_can5_mul),1);
  N = ZXX_evalx0(FpXQX_red(N,polcyclo(p,vT),q));
  P = FpX_mul(N,f,q);
  P = RgX_deflate(P, p);
  d = RgX_splitting(N, p);
  V = cgetg(p+1,t_VEC);
  gel(V,1) = ZX_mulu(gel(d,1), p);
  for(i=2; i<= (long)p; i++)
    gel(V,i) = ZX_mulu(RgX_shift_shallow(gel(d,p+2-i), 1), p);
  (void)delete_var(); return gerepilecopy(av, mkvec2(ZX_sub(f,P),V));
}

static GEN
_can5_invd(void *E, GEN H, GEN v, GEN qM, long M)
{
  ulong p = *(long*)E;
  return gen_ZpX_Dixon(gel(v,2), H, qM, utoipos(p), M, E, _can5_lin, _can5_invl);
}

GEN
Flx_Teichmuller(GEN P, ulong p, long n)
{
  return p==3 ? F3x_frobeniuslift(P,n):
         gen_ZpX_Newton(Flx_to_ZX(P),utoipos(p), n, &p, _can5_iter, _can5_invd);
}

GEN
polteichmuller(GEN P, ulong p, long n)
{
  pari_sp av = avma;
  GEN q = NULL;
  if (typ(P)!=t_POL || !RgX_is_FpX(P,&q)) pari_err_TYPE("polteichmuller",P);
  if (q && !equaliu(q,p)) pari_err_MODULUS("polteichmuller",q,utoi(p));
  if (n <= 0)
    pari_err_DOMAIN("polteichmuller", "precision", "<=",gen_0,stoi(n));
  return gerepileupto(av, p==2 ? F2x_Teichmuller(RgX_to_F2x(P), n)
                               : Flx_Teichmuller(RgX_to_Flx(P, p), p, n));
}
