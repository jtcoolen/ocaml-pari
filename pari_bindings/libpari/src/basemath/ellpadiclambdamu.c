/* Copyright (C) 2011  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

#include "pari.h"
#include "paripriv.h"

/** Batch p-adic logarithms **/
/* a/b mod q */
static ulong
divmodulo(ulong a, ulong b, ulong p, ulong q)
{
  long v = u_lvalrem(b, p, &b);
  if (v) a /= upowuu(p,v);
  /* a/b is now a p-integer */
  return Fl_div(a, b, q);
}
/* to compute log_p(a) mod q = p^n, p < 2^31 */
static GEN
initQplog(long p, ulong q, long n)
{
  long i, nn, nt;
  ulong a, P;
  GEN C;
  for(nn = n, nt = n + 1; nn >= p; nn /= p) nt++;
  if (p == 2)
  {
    P = q - 8;
    while (3 * (nt - 1) > u_lval(nt-1, p) + n) nt--;
  }
  else
  {
    P = q - p;
    while (nt > u_lval(nt-1, p) + n) nt--;
  }
  C = cgetg(nt, t_VECSMALL);
  /* [ P^(k-1) / k, k = 1..nt ] / (p - 1) */
  for(i = 1, a = 1; i < nt; i++, a = Fl_mul(a, P, q))
    C[i] = divmodulo(a, i * (p-1), p, q);
  return C;
}
/* compute log_p(a) / p (resp. log_2(a) / 4) mod p^n, 'a' a p-unit,
 * C = initQplog(p,q,n), q = p^n */
static ulong
logp(GEN C, ulong a, ulong p, ulong q, ulong pn)
{
  long i, b, z, n = lg(C)-1;
  a %= q;
  /* Euclidean quotient (a^2 = 1 mod 8, a^(p-1) = 1 mod p) */
  if (p == 2)
    b = Fl_sqr(a, q << 1) >> 3;
  else
    b = Fl_powu(a, p-1, q) / p;
  z = Fl_mul(b, C[n], pn);
  for (i = n-1; i > 0; i--) z = Fl_mul(z + C[i], b, pn);
  return z;
}

/* p-adic integration against d mu_E, mod p^n, m > 0. D fundamental,
 * W a msinit for k = 2, xpm an eigensymbol.
 * Assume |D p^n| < MAX_LONG. Much cheaper than oms at low precision.
 * Return 1/2 * correct value [half loop over a]. Assume xpm belongs to
 * (-1)^i (D/-1) - part */
static GEN
polPn(GEN W, GEN xpm, long p, long D, long R, long n)
{
  pari_sp av = avma, av2;
  long N = (p == 2)? n + 2: n + 1;
  ulong aD = labs(D), pn = upowuu(p, n), q = upowuu(p, N), ic = 0;
  GEN v, Dq = muluu(aD, q), nc = icopy(gen_1), c = mkfrac(nc, Dq);
  GEN C = n? initQplog(p, q, N): NULL;
  GEN tab = R? ZV_to_Flv(teichmullerinit(p, N), q): NULL;
  long a, A = itou(shifti(Dq,-1));

  if (n) ic = Fl_inv(logp(C, 1 + (p == 2? 4: p), p, q, pn), pn);
  v = zero_zv(pn); av2 = avma;
  for (a = 1; a <= A; a++, set_avma(av2))
  {
    GEN X;
    long s, ap = a % p;
    ulong x, j;
    if (!ap || !(s = kross(D,a))) continue;
    nc[2] = (long)a;
    X = mseval2_ooQ(W, xpm, c); /* xpm(a / Dq) */
    x = umodiu(X, q); if (!x) continue;
    /* log(a) / log(c) */
    j = n? Fl_mul(logp(C, a, p, q, pn), ic, pn): 0;
    if (R) x = Fl_mul(x, tab[Fl_powu(ap, R, p)], q);
    if (s < 0) x = Fl_neg(x, q);
    v[j + 1] = Fl_add(v[j + 1], x, q);
  }
  v = Flv_to_Flx(v, evalvarn(0));
  v = zlx_translate1(v, p, N);
  return gerepileupto(av, Flx_to_ZX(v));
}
/* return lambda, assuming mu = 0 [else infinite loop] */
static long
lambda_ss(GEN W, GEN xpm, long v, long p, long D, long R, long n)
{
  for (;; n += 2)
  {
    GEN P = polPn(W, xpm, p, D, R, n);
    long mu;
    if (!signe(P)) continue;
    mu = ZX_lvalrem(P, p, &P) + v;
    if (!mu)
    {
      long M = upowuu(p,n);
      if (odd(n)) M -= p; else M--;
      M /= p + 1;
      return Flx_val(ZX_to_Flx(P, p)) - M;
    }
  }
}
/* return lambda, assuming mu = 0 [else infinite loop] */
static long
lambda_ord(GEN W, GEN xpm, long v, long p, long D, long R, GEN ap)
{
  GEN P1, P0 = polPn(W, xpm, p, D, R, 0);
  long n;
  for (n = 1;; n++, P0 = P1)
  {
    long mu, pn = upowuu(p,n);
    GEN P, xi, alpha, Q = utoipos(pn);

    P1 = polPn(W, xpm, p, D, R, n);
    alpha = mspadic_unit_eigenvalue(ap, 2, utoipos(p), n+1);
    alpha = padic_to_Q(ginv(alpha));
    xi = FpX_translate(polcyclo(pn, 0), gen_1, Q);
    P = ZX_sub(P1, ZX_Z_mul(FpX_mul(P0, xi, Q), alpha)); /* mod p^n */

    if (!signe(P) || n + v <= 0) continue;
    mu = ZX_lvalrem(P, p, &P) + v;
    if (!mu) return Flx_val(ZX_to_Flx(P, p));
  }
}
GEN
ellpadiclambdamu(GEN E, long p, long D, long R)
{
  pari_sp av = avma;
  long vC, s, muadd = 0;
  GEN W, xpm, C, ap;

  if (!sisfundamental(D))
    pari_err_DOMAIN("ellpadiclambdamu", "isfundamental(D)","=", gen_0, stoi(D));
  s = D < 0? -1: 1;
  if (odd(R)) s = -s;

  ap = ellap(E, utoi(p));
  if (ell_get_type(E) != t_ELL_Q)
    pari_err_TYPE("ellpadiclambdamu", E);
  if (!umodiu(ap, p))
  {
    if (Z_lval(ellQ_get_N(E), p) >= 2)
      pari_err_IMPL("additive reduction in ellpadiclambdamu");
    ap = NULL; /* supersingular */
  }
  if (ap)
  { /* ordinary */
    GEN v = ellisomat(E, p, 1), vE = gel(v,1), M = gel(v,2);
    if (lg(M) != 2) /* some p-isogeny */
    {
      long i, imax = 0, l = lg(vE);
      GEN wmax = NULL, vw = cgetg(l, t_VEC);
      for (i = 1; i < l; i++)
      {
        GEN w, e = ellinit(gel(vE,i), gen_1, 0), em = ellminimalmodel(e, NULL);
        gel(vE,i) = em; obj_free(e);
        w = ellQtwist_bsdperiod(gel(vE,i), s);
        if (s < 0) w = gel(w,2);
        gel(vw,i) = w;
        /* w is a positive real number, either Omega^+ or Omega^- / I */
        if (!imax || gcmp(w, wmax) > 0) { imax = i; wmax = w; }
      }
      if (imax != 1) muadd = Z_lval(ground(gdiv(gel(vw,imax), gel(vw,1))), p);
      for (i = 1; i < l; i++) obj_free(gel(vE,i));
      E = gel(vE, imax);
    }
  }

  W = msfromell(E, s);
  xpm = gel(W,2); W = gel(W,1);
  xpm = Q_primitive_part(xpm, &C);
  vC = C? Q_pval(C, utoipos(p)): 0;
  if (p == 2) vC++; /* due to half-loop in polPn */
  if (vC > 0) pari_err_BUG("ellpadiclambdamu [mu > 0]");

  if (ap)
  { /* ordinary */
    long L = lambda_ord(W, xpm, vC, p, D, R, ap);
    set_avma(av); return mkvec2s(L, muadd);
  }
  else
  {
    long Lp = lambda_ss(W, xpm, vC, p, D, R, 0);
    long Lm = lambda_ss(W, xpm, vC, p, D, R, 1);
    set_avma(av); retmkvec2(mkvec2s(Lp, Lm), zerovec(2));
  }
}
