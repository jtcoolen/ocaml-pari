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

/******************************************************************/
/*                                                                */
/*                 RAMANUJAN's TAU FUNCTION                       */
/*                                                                */
/******************************************************************/
/* 4|N > 0, not fundamental at 2; 6 * Hurwitz class number in level 2,
 * equal to 6*(H(N)+2H(N/4)), H=qfbhclassno */
static GEN
Hspec(GEN N)
{
  long v2 = Z_lvalrem(N, 2, &N), v2f = v2 >> 1;
  GEN t;
  if (odd(v2)) { v2f--; N = shifti(N,3); }
  else if (mod4(N)!=3) { v2f--; N = shifti(N,2); }
  /* N fundamental at 2, v2f = v2(f) s.t. N = f^2 D, D fundamental */
  t = addui(3, muliu(subiu(int2n(v2f+1), 3), 2 - kroiu(N,2)));
  return mulii(t, hclassno6(N));
}

static GEN
tauprime_i(ulong t, GEN p2_7, GEN p_9, GEN p, ulong tin)
{
  GEN h, a, t2 = sqru(t), D = shifti(subii(p, t2), 2); /* 4(p-t^2) */
  /* t mod 2 != tin <=> D not fundamental at 2 */
  h = ((t&1UL) == tin)? hclassno6(D): Hspec(D);
  a = mulii(powiu(t2,3), addii(p2_7, mulii(t2, subii(shifti(t2,2), p_9))));
  return mulii(a, h);
}
GEN
ramanujantau_worker(GEN T, GEN p2_7, GEN p_9, GEN p)
{
  ulong tin = mod4(p) == 3? 1: 0;
  long i, l = lg(T);
  GEN s = gen_0;
  for (i = 1; i < l; i++) s = addii(s, tauprime_i(T[i], p2_7, p_9, p, tin));
  return s;
}

/* B <- {a + k * m : k = 0, ..., (b-a)/m)} */
static void
arithprogset(GEN B, ulong a, ulong b, ulong m)
{
  long k;
  for (k = 1; a <= b; a += m, k++) B[k] = a;
  setlg(B, k);
}
/* sum_{1 <= t <= N} f(t), f has integer values */
static GEN
parsum_u(ulong N, GEN worker)
{
  long a, r, pending = 0, m = usqrt(N);
  struct pari_mt pt;
  GEN v, s = gen_0;
  pari_sp av;

  mt_queue_start_lim(&pt, worker, m);
  v = mkvec(cgetg(m + 2, t_VECSMALL)); av = avma;
  for (r = 1, a = 1; r <= m || pending; r++)
  {
    long workid;
    GEN done;
    if (r <= m) { arithprogset(gel(v,1), a, N, m); a++; }
    mt_queue_submit(&pt, 0, r <= m? v: NULL);
    done = mt_queue_get(&pt, &workid, &pending);
    if (done) s = gerepileuptoint(av, addii(s, done));
  }
  mt_queue_end(&pt); return s;
}

static int
tau_parallel(GEN n) { return mt_nbthreads() > 1 && expi(n) > 18; }

/* Ramanujan tau function for p prime */
static GEN
tauprime(GEN p)
{
  pari_sp av = avma;
  GEN s, p2, p2_7, p_9, T;
  ulong lim;

  if (absequaliu(p, 2)) return utoineg(24);
  /* p > 2 */
  p2 = sqri(p); p2_7 = mului(7, p2); p_9 = mului(9, p);
  lim = itou(sqrtint(p));
  if (tau_parallel(p))
  {
    GEN worker = snm_closure(is_entry("_ramanujantau_worker"),
                             mkvec3(p2_7, p_9, p));
    s = parsum_u(lim, worker);
  }
  else
  {
    pari_sp av2 = avma;
    ulong tin = mod4(p) == 3? 1: 0, t;
    s = gen_0;
    for (t = 1; t <= lim; t++)
    {
      s = addii(s, tauprime_i(t, p2_7, p_9, p, tin));
      if (!(t & 255)) s = gerepileuptoint(av2, s);
    }
  }
  /* 28p^3 - 28p^2 - 90p - 35 */
  T = subii(shifti(mulii(p2_7, subiu(p,1)), 2), addiu(mului(90,p), 35));
  s = shifti(diviuexact(s, 3), 6);
  return gerepileuptoint(av, subii(mulii(mulii(p2, p), T), addui(1, s)));
}

static GEN
taugen_n_i(ulong t, GEN G, GEN n4)
{
  GEN t2 = sqru(t);
  return mulii(mfrhopol_eval(G, t2), hclassno6(subii(n4, t2)));
}
GEN
taugen_n_worker(GEN T, GEN G, GEN n4)
{
  long i, l = lg(T);
  GEN s = gen_0;
  for (i = 1; i < l; i++) s = addii(s, taugen_n_i(T[i], G, n4));
  return s;
}

static GEN
taugen_n(GEN n, GEN G)
{
  GEN S, r, n4 = shifti(n, 2);
  ulong t, lim = itou(sqrtremi(n4, &r));

  if (r == gen_0) lim--;
  G = ZX_unscale(G, n);
  if (tau_parallel(n))
  {
    GEN worker = snm_closure(is_entry("_taugen_n_worker"), mkvec2(G, n4));
    S = parsum_u(lim, worker);
  }
  else
  {
    pari_sp av2 = avma;
    S = gen_0;
    for (t = 1; t <= lim; t++)
    {
      S = addii(S, taugen_n_i(t, G, n4));
      if (!(t & 255)) S = gerepileuptoint(av2, S);
    }
  }
  S = addii(shifti(S,1), mulii(leading_coeff(G), hclassno6(n4)));
  return gdivgu(S, 12);
}

/* ell != 12 */
static GEN
newtrace(GEN fan, GEN n, long ell)
{
  pari_sp av = avma;
  GEN D = divisors(fan), G = mfrhopol(ell-2), T = taugen_n(n, G);
  long i, l = lg(D);

  for (i = 1; i < l; i++)
  {
    GEN d = gel(D, i), q;
    long c = cmpii(sqri(d), n);
    if (c > 0) break;
    q = powiu(d, ell-1);
    if (c < 0) T = gadd(T, q);
    else /* d^2 = n */
    {
      T = gadd(T, gmul2n(q, -1));
      T = gsub(T, gdivgu(mulii(diviiexact(q,d), mfrhopol_eval(G, utoipos(4))), 12));
      break;
    }
  }
  return gerepileuptoint(av, negi(T));
}

#ifdef DEBUG
static void
checkellcong(GEN T, GEN n, long ell)
{
  long V[] = { 0, 691, 0, 3617, 43867, 283*617, 131*593, 0, 657931 };
  if (typ(T) != t_INT
      || umodiu(subii(T, sumdivk(n, ell-1)), V[ell / 2 - 5]))
    pari_err_BUG("ramanujantau");
}
#endif

/* Ramanujan tau function for weights ell = 12, 16, 18, 20, 22, 26,
 * return 0 for <= 0 */
GEN
ramanujantau(GEN n, long ell)
{
  pari_sp av = avma;
  GEN T, P, E, G, F = check_arith_all(n, "ramanujantau");
  long j, lP;

  if (ell < 12 || ell == 14 || odd(ell)) return gen_0;
  if (!F)
  {
    if (signe(n) <= 0) return gen_0;
    F = Z_factor(n); P = gel(F,1);
  }
  else
  {
    P = gel(F,1);
    if (lg(P) == 1 || signe(gel(P,1)) <= 0) return gen_0;
    n = typ(n) == t_VEC? gel(n,1): NULL;
  }
  if (ell > 26 || ell == 24) return newtrace(F, n? n: factorback(F), ell);
  /* dimension 1: tau is multiplicative */
  E = gel(F,2); lP = lg(P); T = gen_1;
  G = ell == 12? NULL: mfrhopol(ell - 2);
  for (j = 1; j < lP; j++)
  {
    GEN p = gel(P,j), q = powiu(p, ell-1), t, t1, t0 = gen_1;
    long k, e = itou(gel(E,j));
    t1 = t = G? subsi(-1, taugen_n(p, G))
              : tauprime(p);
    for (k = 1; k < e; k++)
    {
      GEN t2 = subii(mulii(t, t1), mulii(q, t0));
      t0 = t1; t1 = t2;
    }
    T = mulii(T, t1);
  }
#ifdef DEBUG
  checkellcong(T, n, ell);
#endif
  return gerepileuptoint(av, T);
}
