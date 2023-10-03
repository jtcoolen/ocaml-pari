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
/*********************************************************************/
/**                      PERFECT POWERS                             **/
/*********************************************************************/
#define DEBUGLEVEL DEBUGLEVEL_arith

/*********************************************************************/
/**                     INTEGRAL LOGARITHM                          **/
/*********************************************************************/
/* y > 1 and B > 0 integers. Return e such that y^e <= B < y^(e+1), i.e
 * e = floor(log_y B). Set *ptq = y^e if non-NULL */
long
ulogintall(ulong B, ulong y, ulong *ptq)
{
  ulong r, r2;
  long e;

  if (y == 2)
  {
    long eB = expu(B); /* 2^eB <= B < 2^(eB + 1) */
    if (ptq) *ptq = 1UL << eB;
    return eB;
  }
  r = y, r2 = 1UL;
  for (e=1;; e++)
  { /* here, r = y^e, r2 = y^(e-1) */
    if (r >= B)
    {
      if (r != B) { e--; r = r2; }
      if (ptq) *ptq = r;
      return e;
    }
    r2 = r;
    r = umuluu_or_0(y, r);
    if (!r)
    {
      if (ptq) *ptq = r2;
      return e;
    }
  }
}

/* y > 1 and B > 0 integers. Return e such that y^e <= B < y^(e+1), i.e
 * e = floor(log_y B). Set *ptq = y^e if non-NULL */
long
logintall(GEN B, GEN y, GEN *ptq)
{
  pari_sp av;
  long ey, e, emax, i, eB = expi(B); /* 2^eB <= B < 2^(eB + 1) */
  GEN q, pow2;

  if (lgefint(B) == 3)
  {
    ulong q;
    if (lgefint(y) > 3)
    {
      if (ptq) *ptq = gen_1;
      return 0;
    }
    if (!ptq) return ulogintall(B[2], y[2], NULL);
    e = ulogintall(B[2], y[2], &q);
    *ptq = utoi(q); return e;
  }
  if (equaliu(y,2))
  {
    if (ptq) *ptq = int2n(eB);
    return eB;
  }
  av = avma;
  ey = expi(y);
  /* eB/(ey+1) - 1 < e <= eB/ey */
  emax = eB/ey;
  if (emax <= 13) /* e small, be naive */
  {
    GEN r = y, r2 = gen_1;
    for (e=1;; e++)
    { /* here, r = y^e, r2 = y^(e-1) */
      long fl = cmpii(r, B);
      if (fl >= 0)
      {
        if (fl) { e--; cgiv(r); r = r2; }
        if (ptq) *ptq = gerepileuptoint(av, r); else set_avma(av);
        return e;
      }
      r2 = r; r = mulii(r,y);
    }
  }
  /* e >= 13 ey / (ey+1) >= 6.5 */

  /* binary splitting: compute bits of e one by one */
  /* compute pow2[i] = y^(2^i) [i < crude upper bound for log_2 log_y(B)] */
  pow2 = new_chunk((long)log2(eB)+2);
  gel(pow2,0) = y;
  for (i=0, q=y;; )
  {
    GEN r = gel(pow2,i); /* r = y^2^i */
    long fl = cmpii(r,B);
    if (!fl)
    {
      e = 1L<<i;
      if (ptq) *ptq = gerepileuptoint(av, r); else set_avma(av);
      return e;
    }
    if (fl > 0) { i--; break; }
    q = r;
    if (1L<<(i+1) > emax) break;
    gel(pow2,++i) = sqri(q);
  }

  for (e = 1L<<i;;)
  { /* y^e = q < B < r = q * y^(2^i) */
    pari_sp av2 = avma;
    long fl;
    GEN r;
    if (--i < 0) break;
    r = mulii(q, gel(pow2,i));
    fl = cmpii(r, B);
    if (fl > 0) set_avma(av2);
    else
    {
      e += (1L<<i);
      q = r;
      if (!fl) break; /* B = r */
    }
  }
  if (ptq) *ptq = gerepileuptoint(av, q); else set_avma(av);
  return e;
}

long
logint0(GEN B, GEN y, GEN *ptq)
{
  const char *f = "logint";
  if (typ(y) != t_INT) pari_err_TYPE(f,y);
  if (cmpis(y, 2) < 0) pari_err_DOMAIN(f, "b" ,"<=", gen_1, y);
  if (typ(B) != t_INT)
  {
    pari_sp av = avma;
    long a;
    if (typ(B) == t_REAL)
    {
      long e, p;
      if (cmprs(B, 1) < 1) pari_err_DOMAIN(f, "x", "<", gen_1, B);
      e = expo(B); if (e < 0) return 0;
      if (equaliu(y, 2)) return e;
      if (expu(e) < 50)
      {
        a = floor(dbllog2(B) / dbllog2(y));
        if (ptq) *ptq = powiu(y, a);
        return a;
      }
      /* play safe */
      p = lg(B);
      if (nbits2lg(e+1) > p)
      { /* try to avoid precision loss in truncation */
        if (p > DEFAULTPREC) { p = DEFAULTPREC; B = rtor(B, p); }
        a = itos(floorr(divrr(logr_abs(B), logr_abs(itor(y, p)))));
        set_avma(av); if (ptq) *ptq = powiu(y, a);
        return a;
      }
      a = logintall(truncr(B), y, ptq);
    }
    else
    {
      GEN b = gfloor(B);
      if (typ(b) != t_INT) pari_err_TYPE(f,B);
      if (signe(b) <= 0) pari_err_DOMAIN(f, "x", "<", gen_1, B);
      a = logintall(b, y, ptq);
    }
    if (!ptq) return gc_long(av, a);
    *ptq = gerepileuptoint(av, *ptq); return a;
  }
  if (signe(B) <= 0) pari_err_DOMAIN(f, "x" ,"<=", gen_0, B);
  return logintall(B,y,ptq);
}

/*********************************************************************/
/**                     INTEGRAL SQUARE ROOT                        **/
/*********************************************************************/
GEN
sqrtint(GEN a)
{
  if (typ(a) != t_INT)
  {
    pari_sp av = avma;
    if (typ(a) == t_REAL)
    {
      long e;
      switch(signe(a))
      {
        case 0: return gen_0;
        case -1: pari_err_DOMAIN("sqrtint", "argument", "<", gen_0,a);
      }
      e = expo(a); if (e < 0) return gen_0;
      if (nbits2lg(e+1) > lg(a))
        a = floorr(sqrtr(a)); /* try to avoid precision loss in truncation */
      else
        a = sqrti(truncr(a));
    }
    else
    {
      GEN b = gfloor(a);
      if (typ(b) != t_INT) pari_err_TYPE("sqrtint",a);
      if (signe(b) < 0) pari_err_DOMAIN("sqrtint", "argument", "<", gen_0,a);
      a = sqrti(b);
    }
    return gerepileuptoleaf(av, a);
  }
  switch (signe(a))
  {
    case 1: return sqrti(a);
    case 0: return gen_0;
    default: pari_err_DOMAIN("sqrtint", "argument", "<", gen_0,a);
  }
  return NULL; /* LCOV_EXCL_LINE */
}
GEN
sqrtint0(GEN a, GEN *r)
{
  if (!r) return sqrtint(a);
  if (typ(a) != t_INT)
  {
    GEN b = sqrtint(a);
    pari_sp av = avma;
    *r = gerepileupto(av, gsub(a, sqri(b))); return b;
  }
  switch (signe(a))
  {
    case 1: return sqrtremi(a, r);
    case 0: *r = gen_0; return gen_0;
    default: pari_err_DOMAIN("sqrtint", "argument", "<", gen_0,a);
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/*********************************************************************/
/**                      PERFECT SQUARE                             **/
/*********************************************************************/
static int
squaremod(ulong A)
{
  const int squaremod64[]={
    1,1,0,0,1,0,0,0,0,1, 0,0,0,0,0,0,1,1,0,0, 0,0,0,0,0,1,0,0,0,0,
    0,0,0,1,0,0,1,0,0,0, 0,1,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,1,0,0, 0,0,0,0};
  const int squaremod63[]={
    1,1,0,0,1,0,0,1,0,1, 0,0,0,0,0,0,1,0,1,0, 0,0,1,0,0,1,0,0,1,0,
    0,0,0,0,0,0,1,1,0,0, 0,0,0,1,0,0,1,0,0,1, 0,0,0,0,0,0,0,0,1,0, 0,0,0};
  const int squaremod65[]={
    1,1,0,0,1,0,0,0,0,1, 1,0,0,0,1,0,1,0,0,0, 0,0,0,0,0,1,1,0,0,1,
    1,0,0,0,0,1,1,0,0,1, 1,0,0,0,0,0,0,0,0,1, 0,1,0,0,0,1,1,0,0,0, 0,1,0,0,1};
  const int squaremod11[]={1,1,0,1,1,1,0,0,0,1, 0};
  return (squaremod64[A & 0x3fUL]
    && squaremod63[A % 63UL]
    && squaremod65[A % 65UL]
    && squaremod11[A % 11UL]);
}

/* emulate Z_issquareall on single-word integers */
long
uissquareall(ulong A, ulong *sqrtA)
{
  if (!A) { *sqrtA = 0; return 1; }
  if (squaremod(A))
  {
    ulong a = usqrt(A);
    if (a * a == A) { *sqrtA = a; return 1; }
  }
  return 0;
}
long
uissquare(ulong A)
{
  if (!A) return 1;
  if (squaremod(A)) { ulong a = usqrt(A); if (a * a == A) return 1; }
  return 0;
}

long
Z_issquareall(GEN x, GEN *pt)
{
  pari_sp av;
  GEN y, r;

  switch(signe(x))
  {
    case -1: return 0;
    case 0: if (pt) *pt=gen_0; return 1;
  }
  if (lgefint(x) == 3)
  {
    ulong u = uel(x,2), a;
    if (!pt) return uissquare(u);
    if (!uissquareall(u, &a)) return 0;
    *pt = utoipos(a); return 1;
  }
  if (!squaremod(umodiu(x, 64*63*65*11))) return 0;
  av = avma; y = sqrtremi(x, &r);
  if (r != gen_0) return gc_long(av,0);
  if (pt) { *pt = y; set_avma((pari_sp)y); } else set_avma(av);
  return 1;
}

/* a t_INT, p prime */
long
Zp_issquare(GEN a, GEN p)
{
  long v;
  GEN ap;

  if (!signe(a) || equali1(a)) return 1;
  v = Z_pvalrem(a, p, &ap);
  if (v&1) return 0;
  return absequaliu(p, 2)? Mod8(ap) == 1
                         : kronecker(ap,p) == 1;
}

static long
polissquareall(GEN x, GEN *pt)
{
  pari_sp av;
  long v;
  GEN y, a, b, p;

  if (!signe(x))
  {
    if (pt) *pt = gcopy(x);
    return 1;
  }
  if (odd(degpol(x))) return 0; /* odd degree */
  av = avma;
  v = RgX_valrem(x, &x);
  if (v & 1) return gc_long(av,0);
  a = gel(x,2); /* test constant coeff */
  if (!pt)
  { if (!issquare(a)) return gc_long(av,0); }
  else
  { if (!issquareall(a,&b)) return gc_long(av,0); }
  if (!degpol(x)) { /* constant polynomial */
    if (!pt) return gc_long(av,1);
    y = scalarpol(b, varn(x)); goto END;
  }
  p = characteristic(x);
  if (signe(p) && !mod2(p))
  {
    long i, lx;
    if (!absequaliu(p,2)) pari_err_IMPL("issquare for even characteristic != 2");
    x = gmul(x, mkintmod(gen_1, gen_2));
    lx = lg(x);
    if ((lx-3) & 1) return gc_long(av,0);
    for (i = 3; i < lx; i+=2)
      if (!gequal0(gel(x,i))) return gc_long(av,0);
    if (pt) {
      y = cgetg((lx+3) / 2, t_POL);
      for (i = 2; i < lx; i+=2)
        if (!issquareall(gel(x,i), &gel(y, (i+2)>>1))) return gc_long(av,0);
      y[1] = evalsigne(1) | evalvarn(varn(x));
      goto END;
    } else {
      for (i = 2; i < lx; i+=2)
        if (!issquare(gel(x,i))) return gc_long(av,0);
      return gc_long(av,1);
    }
  }
  else
  {
    long m = 1;
    x = RgX_Rg_div(x,a);
    /* a(x^m) = B^2 => B = b(x^m) provided a(0) != 0 */
    if (!signe(p)) x = RgX_deflate_max(x,&m);
    y = ser2rfrac_i(gsqrt(RgX_to_ser(x,lg(x)-1),0));
    if (!RgX_equal(RgX_sqr(y), x)) return gc_long(av,0);
    if (!pt) return gc_long(av,1);
    if (!gequal1(a)) y = gmul(b, y);
    if (m != 1) y = RgX_inflate(y,m);
  }
END:
  if (v) y = RgX_shift_shallow(y, v>>1);
  *pt = gerepilecopy(av, y); return 1;
}

/* b unit mod p */
static int
Up_ispower(GEN b, GEN K, GEN p, long d, GEN *pt)
{
  if (d == 1)
  { /* mod p: faster */
    if (!Fp_ispower(b, K, p)) return 0;
    if (pt) *pt = Fp_sqrtn(b, K, p, NULL);
  }
  else
  { /* mod p^{2 +} */
    if (!ispower(cvtop(b, p, d), K, pt)) return 0;
    if (pt) *pt = gtrunc(*pt);
  }
  return 1;
}

/* We're studying whether a mod (q*p^e) is a K-th power, (q,p) = 1.
 * Decide mod p^e, then reduce a mod q unless q = NULL. */
static int
handle_pe(GEN *pa, GEN q, GEN L, GEN K, GEN p, long e)
{
  GEN t, A;
  long v = Z_pvalrem(*pa, p, &A), d = e - v;
  if (d <= 0) t = gen_0;
  else
  {
    ulong r;
    v = uabsdivui_rem(v, K, &r);
    if (r || !Up_ispower(A, K, p, d, L? &t: NULL)) return 0;
    if (L && v) t = mulii(t, powiu(p, v));
  }
  if (q) *pa = modii(*pa, q);
  if (L) vectrunc_append(L, mkintmod(t, powiu(p, e)));
  return 1;
}
long
Zn_ispower(GEN a, GEN q, GEN K, GEN *pt)
{
  GEN L, N;
  pari_sp av;
  long e, i, l;
  ulong pp;
  forprime_t S;

  if (!signe(a))
  {
    if (pt) {
      GEN t = cgetg(3, t_INTMOD);
      gel(t,1) = icopy(q); gel(t,2) = gen_0; *pt = t;
    }
    return 1;
  }
  /* a != 0 */
  av = avma;

  if (typ(q) != t_INT) /* integer factorization */
  {
    GEN P = gel(q,1), E = gel(q,2);
    l = lg(P);
    L = pt? vectrunc_init(l): NULL;
    for (i = 1; i < l; i++)
    {
      GEN p = gel(P,i);
      long e = itos(gel(E,i));
      if (!handle_pe(&a, NULL, L, K, p, e)) return gc_long(av,0);
    }
    goto END;
  }
  if (!mod2(K)
      && kronecker(a, shifti(q,-vali(q))) == -1) return gc_long(av,0);
  L = pt? vectrunc_init(expi(q)+1): NULL;
  u_forprime_init(&S, 2, tridiv_bound(q));
  while ((pp = u_forprime_next(&S)))
  {
    int stop;
    e = Z_lvalrem_stop(&q, pp, &stop);
    if (!e) continue;
    if (!handle_pe(&a, q, L, K, utoipos(pp), e)) return gc_long(av,0);
    if (stop)
    {
      if (!is_pm1(q) && !handle_pe(&a, q, L, K, q, 1)) return gc_long(av,0);
      goto END;
    }
  }
  l = lg(primetab);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(primetab,i);
    e = Z_pvalrem(q, p, &q);
    if (!e) continue;
    if (!handle_pe(&a, q, L, K, p, e)) return gc_long(av,0);
    if (is_pm1(q)) goto END;
  }
  N = gcdii(a,q);
  if (!is_pm1(N))
  {
    if (ifac_isprime(N))
    {
      e = Z_pvalrem(q, N, &q);
      if (!handle_pe(&a, q, L, K, N, e)) return gc_long(av,0);
    }
    else
    {
      GEN part = ifac_start(N, 0);
      for(;;)
      {
        long e;
        GEN p;
        if (!ifac_next(&part, &p, &e)) break;
        e = Z_pvalrem(q, p, &q);
        if (!handle_pe(&a, q, L, K, p, e)) return gc_long(av,0);
      }
    }
  }
  if (!is_pm1(q))
  {
    if (ifac_isprime(q))
    {
      if (!handle_pe(&a, q, L, K, q, 1)) return gc_long(av,0);
    }
    else
    {
      GEN part = ifac_start(q, 0);
      for(;;)
      {
        long e;
        GEN p;
        if (!ifac_next(&part, &p, &e)) break;
        if (!handle_pe(&a, q, L, K, p, e)) return gc_long(av,0);
      }
    }
  }
END:
  if (pt) *pt = gerepileupto(av, chinese1_coprime_Z(L));
  return 1;
}

static long
polmodispower(GEN x, GEN K, GEN *pt)
{
  pari_sp av = avma;
  GEN p = NULL, T = NULL;
  if (Rg_is_FpXQ(x, &T,&p) && p)
  {
    x = liftall_shallow(x);
    if (T) T = liftall_shallow(T);
    if (!Fq_ispower(x, K, T, p)) return gc_long(av,0);
    if (!pt) return gc_long(av,1);
    x = Fq_sqrtn(x, K, T,p, NULL);
    if (typ(x) == t_INT)
      x = Fp_to_mod(x,p);
    else
      x = mkpolmod(FpX_to_mod(x,p), FpX_to_mod(T,p));
    *pt = gerepilecopy(av, x); return 1;
  }
  pari_err_IMPL("ispower for general t_POLMOD");
  return 0;
}
static long
rfracispower(GEN x, GEN K, GEN *pt)
{
  pari_sp av = avma;
  GEN n = gel(x,1), d = gel(x,2);
  long v = -RgX_valrem(d, &d), vx = varn(d);
  if (typ(n) == t_POL && varn(n) == vx) v += RgX_valrem(n, &n);
  if (!dvdsi(v, K)) return gc_long(av, 0);
  if (lg(d) >= 3)
  {
    GEN a = gel(d,2); /* constant term */
    if (!gequal1(a)) { d = RgX_Rg_div(d, a); n = gdiv(n, a); }
  }
  if (!ispower(d, K, pt? &d: NULL) || !ispower(n, K, pt? &n: NULL))
    return gc_long(av, 0);
  if (!pt) return gc_long(av, 1);
  x = gdiv(n, d);
  if (v) x = gmul(x, monomial(gen_1, v / itos(K), vx));
  *pt = gerepileupto(av, x); return 1;
}
long
issquareall(GEN x, GEN *pt)
{
  long tx = typ(x);
  GEN F;
  pari_sp av;

  if (!pt) return issquare(x);
  switch(tx)
  {
    case t_INT: return Z_issquareall(x, pt);
    case t_FRAC: av = avma;
      F = cgetg(3, t_FRAC);
      if (   !Z_issquareall(gel(x,1), &gel(F,1))
          || !Z_issquareall(gel(x,2), &gel(F,2))) return gc_long(av,0);
      *pt = F; return 1;

    case t_POLMOD:
      return polmodispower(x, gen_2, pt);
    case t_POL: return polissquareall(x,pt);
    case t_RFRAC: return rfracispower(x, gen_2, pt);

    case t_REAL: case t_COMPLEX: case t_PADIC: case t_SER:
      if (!issquare(x)) return 0;
      *pt = gsqrt(x, DEFAULTPREC); return 1;

    case t_INTMOD:
      return Zn_ispower(gel(x,2), gel(x,1), gen_2, pt);

    case t_FFELT: return FF_issquareall(x, pt);

  }
  pari_err_TYPE("issquareall",x);
  return 0; /* LCOV_EXCL_LINE */
}

long
issquare(GEN x)
{
  GEN a, p;
  long v;

  switch(typ(x))
  {
    case t_INT:
      return Z_issquare(x);

    case t_REAL:
      return (signe(x)>=0);

    case t_INTMOD:
      return Zn_ispower(gel(x,2), gel(x,1), gen_2, NULL);

    case t_FRAC:
      return Z_issquare(gel(x,1)) && Z_issquare(gel(x,2));

    case t_FFELT: return FF_issquareall(x, NULL);

    case t_COMPLEX:
      return 1;

    case t_PADIC:
      a = gel(x,4); if (!signe(a)) return 1;
      if (valp(x)&1) return 0;
      p = gel(x,2);
      if (!absequaliu(p, 2)) return (kronecker(a,p) != -1);

      v = precp(x); /* here p=2, a is odd */
      if ((v>=3 && mod8(a) != 1 ) ||
          (v==2 && mod4(a) != 1)) return 0;
      return 1;

    case t_POLMOD:
      return polmodispower(x, gen_2, NULL);

    case t_POL:
      return polissquareall(x,NULL);

    case t_SER:
      if (!signe(x)) return 1;
      if (valser(x)&1) return 0;
      return issquare(gel(x,2));

    case t_RFRAC:
      return rfracispower(x, gen_2, NULL);
  }
  pari_err_TYPE("issquare",x);
  return 0; /* LCOV_EXCL_LINE */
}
GEN gissquare(GEN x) { return issquare(x)? gen_1: gen_0; }
GEN gissquareall(GEN x, GEN *pt) { return issquareall(x,pt)? gen_1: gen_0; }

long
ispolygonal(GEN x, GEN S, GEN *N)
{
  pari_sp av = avma;
  GEN D, d, n;
  if (typ(x) != t_INT) pari_err_TYPE("ispolygonal", x);
  if (typ(S) != t_INT) pari_err_TYPE("ispolygonal", S);
  if (abscmpiu(S,3) < 0) pari_err_DOMAIN("ispolygonal","s","<", utoipos(3),S);
  if (signe(x) < 0) return 0;
  if (signe(x) == 0) { if (N) *N = gen_0; return 1; }
  if (is_pm1(x)) { if (N) *N = gen_1; return 1; }
  /* n = (sqrt( (8s - 16) x + (s-4)^2 ) + s - 4) / 2(s - 2) */
  if (abscmpiu(S, 1<<16) < 0) /* common case ! */
  {
    ulong s = S[2], r;
    if (s == 4) return Z_issquareall(x, N);
    if (s == 3)
      D = addiu(shifti(x, 3), 1);
    else
      D = addiu(mului(8*s - 16, x), (s-4)*(s-4));
    if (!Z_issquareall(D, &d)) return gc_long(av,0);
    if (s == 3)
      d = subiu(d, 1);
    else
      d = addiu(d, s - 4);
    n = absdiviu_rem(d, 2*s - 4, &r);
    if (r) return gc_long(av,0);
  }
  else
  {
    GEN r, S_2 = subiu(S,2), S_4 = subiu(S,4);
    D = addii(mulii(shifti(S_2,3), x), sqri(S_4));
    if (!Z_issquareall(D, &d)) return gc_long(av,0);
    d = addii(d, S_4);
    n = dvmdii(shifti(d,-1), S_2, &r);
    if (r != gen_0) return gc_long(av,0);
  }
  if (N) *N = gerepileuptoint(av, n); else set_avma(av);
  return 1;
}

/*********************************************************************/
/**                        PERFECT POWER                            **/
/*********************************************************************/
static long
polispower(GEN x, GEN K, GEN *pt)
{
  pari_sp av;
  long v, d, k = itos(K);
  GEN y, a, b;
  GEN T = NULL, p = NULL;

  if (!signe(x))
  {
    if (pt) *pt = gcopy(x);
    return 1;
  }
  d = degpol(x);
  if (d % k) return 0; /* degree not multiple of k */
  av = avma;
  if (RgX_is_FpXQX(x, &T, &p) && p)
  { /* over Fq */
    if (T && typ(T) == t_FFELT)
    {
      if (!FFX_ispower(x, k, T, pt)) return gc_long(av,0);
      return 1;
    }
    x = RgX_to_FqX(x,T,p);
    if (!FqX_ispower(x, k, T,p, pt)) return gc_long(av,0);
    if (pt) *pt = gerepileupto(av, FqX_to_mod(*pt, T, p));
    return 1;
  }
  v = RgX_valrem(x, &x);
  if (v % k) return 0;
  v /= k;
  a = gel(x,2); b = NULL;
  if (!ispower(a, K, &b)) return gc_long(av,0);
  if (d)
  {
    GEN p = characteristic(x);
    a = leading_coeff(x);
    if (!ispower(a, K, &b)) return gc_long(av,0);
    x = RgX_normalize(x);
    if (signe(p) && cmpii(p,K) <= 0)
      pari_err_IMPL("ispower(general t_POL) in small characteristic");
    y = gtrunc(gsqrtn(RgX_to_ser(x,lg(x)), K, NULL, 0));
    if (!RgX_equal(powgi(y, K), x)) return gc_long(av,0);
  }
  else
    y = pol_1(varn(x));
  if (pt)
  {
    if (!gequal1(a))
    {
      if (!b) b = gsqrtn(a, K, NULL, DEFAULTPREC);
      y = gmul(b,y);
    }
    if (v) y = RgX_shift_shallow(y, v);
    *pt = gerepilecopy(av, y);
  }
  else set_avma(av);
  return 1;
}

long
Z_ispowerall(GEN x, ulong k, GEN *pt)
{
  long s = signe(x);
  ulong mask;
  if (!s) { if (pt) *pt = gen_0; return 1; }
  if (s > 0) {
    if (k == 2) return Z_issquareall(x, pt);
    if (k == 3) { mask = 1; return !!is_357_power(x, pt, &mask); }
    if (k == 5) { mask = 2; return !!is_357_power(x, pt, &mask); }
    if (k == 7) { mask = 4; return !!is_357_power(x, pt, &mask); }
    return is_kth_power(x, k, pt);
  }
  if (!odd(k)) return 0;
  if (Z_ispowerall(absi_shallow(x), k, pt))
  {
    if (pt) *pt = negi(*pt);
    return 1;
  };
  return 0;
}

/* is x a K-th power mod p ? Assume p prime. */
int
Fp_ispower(GEN x, GEN K, GEN p)
{
  pari_sp av = avma;
  GEN p_1;
  x = modii(x, p);
  if (!signe(x) || equali1(x)) return gc_bool(av,1);
  /* implies p > 2 */
  p_1 = subiu(p,1);
  K = gcdii(K, p_1);
  if (absequaliu(K, 2)) return gc_bool(av, kronecker(x,p) > 0);
  x = Fp_pow(x, diviiexact(p_1,K), p);
  return gc_bool(av, equali1(x));
}

/* x unit defined modulo 2^e, e > 0, p prime */
static int
U2_issquare(GEN x, long e)
{
  long r = signe(x)>=0?mod8(x):8-mod8(x);
  if (e==1) return 1;
  if (e==2) return (r&3L) == 1;
  return r == 1;
}
/* x unit defined modulo p^e, e > 0, p prime */
static int
Up_issquare(GEN x, GEN p, long e)
{ return (absequaliu(p,2))? U2_issquare(x, e): kronecker(x,p)==1; }

long
Zn_issquare(GEN d, GEN fn)
{
  long j, np;
  if (typ(d) != t_INT) pari_err_TYPE("Zn_issquare",d);
  if (typ(fn) == t_INT) return Zn_ispower(d, fn, gen_2, NULL);
  /* integer factorization */
  np = nbrows(fn);
  for (j = 1; j <= np; ++j)
  {
    GEN  r, p = gcoeff(fn, j, 1);
    long e = itos(gcoeff(fn, j, 2));
    long v = Z_pvalrem(d,p,&r);
    if (v < e && (odd(v) || !Up_issquare(r, p, e-v))) return 0;
  }
  return 1;
}

static long
Qp_ispower(GEN x, GEN K, GEN *pt)
{
  pari_sp av = avma;
  GEN z = Qp_sqrtn(x, K, NULL);
  if (!z) return gc_long(av,0);
  if (pt) *pt = z;
  return 1;
}

long
ispower(GEN x, GEN K, GEN *pt)
{
  GEN z;

  if (!K) return gisanypower(x, pt);
  if (typ(K) != t_INT) pari_err_TYPE("ispower",K);
  if (signe(K) <= 0) pari_err_DOMAIN("ispower","exponent","<=",gen_0,K);
  if (equali1(K)) { if (pt) *pt = gcopy(x); return 1; }
  switch(typ(x)) {
    case t_INT:
      if (lgefint(K) != 3) return 0;
      return Z_ispowerall(x, itou(K), pt);
    case t_FRAC:
    {
      GEN a = gel(x,1), b = gel(x,2);
      ulong k;
      if (lgefint(K) != 3) return 0;
      k = itou(K);
      if (pt) {
        z = cgetg(3, t_FRAC);
        if (Z_ispowerall(a, k, &a) && Z_ispowerall(b, k, &b)) {
          *pt = z; gel(z,1) = a; gel(z,2) = b; return 1;
        }
        set_avma((pari_sp)(z + 3)); return 0;
      }
      return Z_ispower(a, k) && Z_ispower(b, k);
    }
    case t_INTMOD:
      return Zn_ispower(gel(x,2), gel(x,1), K, pt);
    case t_FFELT:
      return FF_ispower(x, K, pt);

    case t_PADIC:
      return Qp_ispower(x, K, pt);
    case t_POLMOD:
      return polmodispower(x, K, pt);
    case t_POL:
      return polispower(x, K, pt);
    case t_RFRAC:
      return rfracispower(x, K, pt);
    case t_REAL:
      if (signe(x) < 0 && !mpodd(K)) return 0;
    case t_COMPLEX:
      if (pt) *pt = gsqrtn(x, K, NULL, DEFAULTPREC);
      return 1;

    case t_SER:
      if (signe(x) && (!dvdsi(valser(x), K) || !ispower(gel(x,2), K, NULL)))
        return 0;
      if (pt) *pt = gsqrtn(x, K, NULL, DEFAULTPREC);
      return 1;
  }
  pari_err_TYPE("ispower",x);
  return 0; /* LCOV_EXCL_LINE */
}

long
gisanypower(GEN x, GEN *pty)
{
  long tx = typ(x);
  ulong k, h;
  if (tx == t_INT) return Z_isanypower(x, pty);
  if (tx == t_FRAC)
  {
    pari_sp av = avma;
    GEN fa, P, E, a = gel(x,1), b = gel(x,2);
    long i, j, p, e;
    int sw = (abscmpii(a, b) > 0);

    if (sw) swap(a, b);
    k = Z_isanypower(a, pty? &a: NULL);
    if (!k)
    { /* a = -1,1 or not a pure power */
      if (!is_pm1(a)) return gc_long(av,0);
      if (signe(a) < 0) b = negi(b);
      k = Z_isanypower(b, pty? &b: NULL);
      if (!k || !pty) return gc_long(av,k);
      *pty = gerepileupto(av, ginv(b));
      return k;
    }
    fa = factoru(k);
    P = gel(fa,1);
    E = gel(fa,2); h = k;
    for (i = lg(P) - 1; i > 0; i--)
    {
      p = P[i];
      e = E[i];
      for (j = 0; j < e; j++)
        if (!is_kth_power(b, p, &b)) break;
      if (j < e) k /= upowuu(p, e - j);
    }
    if (k == 1) return gc_long(av,0);
    if (!pty) return gc_long(av,k);
    if (k != h) a = powiu(a, h/k);
    *pty = gerepilecopy(av, mkfrac(a, b));
    return k;
  }
  pari_err_TYPE("gisanypower", x);
  return 0; /* LCOV_EXCL_LINE */
}

/* v_p(x) = e != 0 for some p; return ispower(x,,&x), updating x.
 * No need to optimize for 2,3,5,7 powers (done before) */
static long
split_exponent(ulong e, GEN *x)
{
  GEN fa, P, E;
  long i, j, l, k = 1;
  if (e == 1) return 1;
  fa = factoru(e);
  P = gel(fa,1);
  E = gel(fa,2); l = lg(P);
  for (i = 1; i < l; i++)
  {
    ulong p = P[i];
    for (j = 0; j < E[i]; j++)
    {
      GEN y;
      if (!is_kth_power(*x, p, &y)) break;
      k *= p; *x = y;
    }
  }
  return k;
}

static long
Z_isanypower_nosmalldiv(GEN *px)
{ /* any prime divisor of x is > 102 */
  const double LOG2_103 = 6.6865; /* lower bound for log_2(103) */
  const double LOG103 = 4.6347; /* lower bound for log(103) */
  forprime_t T;
  ulong mask = 7, e2;
  long k, ex;
  GEN y, x = *px;

  k = 1;
  while (Z_issquareall(x, &y)) { k <<= 1; x = y; }
  while ( (ex = is_357_power(x, &y, &mask)) ) { k *= ex; x = y; }
  e2 = (ulong)((expi(x) + 1) / LOG2_103); /* >= log_103 (x) */
  if (u_forprime_init(&T, 11, e2))
  {
    GEN logx = NULL;
    const ulong Q = 30011; /* prime */
    ulong p, xmodQ;
    double dlogx = 0;
    /* cut off at x^(1/p) ~ 2^30 bits which seems to be about optimum;
     * for large p the modular checks are no longer competitively fast */
    while ( (ex = is_pth_power(x, &y, &T, 30)) )
    {
      k *= ex; x = y;
      e2 = (ulong)((expi(x) + 1) / LOG2_103);
      u_forprime_restrict(&T, e2);
    }
    if (DEBUGLEVEL>4)
      err_printf("Z_isanypower: now k=%ld, x=%ld-bit\n", k, expi(x)+1);
    xmodQ = umodiu(x, Q);
    /* test Q | x, just in case */
    if (!xmodQ) { *px = x; return k * split_exponent(Z_lval(x,Q), px); }
    /* x^(1/p) < 2^31 */
    p = T.p;
    if (p <= e2)
    {
      logx = logr_abs( itor(x, DEFAULTPREC) );
      dlogx = rtodbl(logx);
      e2 = (ulong)(dlogx / LOG103); /* >= log_103(x) */
    }
    while (p && p <= e2)
    { /* is x a p-th power ? By computing y = round(x^(1/p)).
       * Check whether y^p = x, first mod Q, then exactly. */
      pari_sp av = avma;
      long e;
      GEN logy = divru(logx, p), y = grndtoi(mpexp(logy), &e);
      ulong ymodQ = umodiu(y,Q);
      if (e >= -10 || Fl_powu(ymodQ, p % (Q-1), Q) != xmodQ
                   || !equalii(powiu(y, p), x)) set_avma(av);
      else
      {
        k *= p; x = y; xmodQ = ymodQ; logx = logy; dlogx /= p;
        e2 = (ulong)(dlogx / LOG103); /* >= log_103(x) */
        u_forprime_restrict(&T, e2);
        continue; /* if success, retry same p */
      }
      p = u_forprime_next(&T);
    }
  }
  *px = x; return k;
}

static ulong tinyprimes[] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
  157, 163, 167, 173, 179, 181, 191, 193, 197, 199
};

/* disregard the sign of x, caller will take care of x < 0 */
static long
Z_isanypower_aux(GEN x, GEN *pty)
{
  long ex, v, i, l, k;
  GEN y, P, E;
  ulong mask, e = 0;

  if (abscmpii(x, gen_2) < 0) return 0; /* -1,0,1 */

  if (signe(x) < 0) x = negi(x);
  k = l = 1;
  P = cgetg(26 + 1, t_VECSMALL);
  E = cgetg(26 + 1, t_VECSMALL);
  /* trial division */
  for(i = 0; i < 26; i++)
  {
    ulong p = tinyprimes[i];
    int stop;
    v = Z_lvalrem_stop(&x, p, &stop);
    if (v)
    {
      P[l] = p;
      E[l] = v; l++;
      e = ugcd(e, v); if (e == 1) goto END;
    }
    if (stop) {
      if (is_pm1(x)) k = e;
      goto END;
    }
  }

  if (e)
  { /* Bingo. Result divides e */
    long v3, v5, v7;
    ulong e2 = e;
    v = u_lvalrem(e2, 2, &e2);
    if (v)
    {
      for (i = 0; i < v; i++)
      {
        if (!Z_issquareall(x, &y)) break;
        k <<= 1; x = y;
      }
    }
    mask = 0;
    v3 = u_lvalrem(e2, 3, &e2); if (v3) mask = 1;
    v5 = u_lvalrem(e2, 5, &e2); if (v5) mask |= 2;
    v7 = u_lvalrem(e2, 7, &e2); if (v7) mask |= 4;
    while ( (ex = is_357_power(x, &y, &mask)) ) {
      x = y;
      switch(ex)
      {
        case 3: k *= 3; if (--v3 == 0) mask &= ~1; break;
        case 5: k *= 5; if (--v5 == 0) mask &= ~2; break;
        case 7: k *= 7; if (--v7 == 0) mask &= ~4; break;
      }
    }
    k *= split_exponent(e2, &x);
  }
  else
    k = Z_isanypower_nosmalldiv(&x);
END:
  if (pty && k != 1)
  {
    if (e)
    { /* add missing small factors */
      y = powuu(P[1], E[1] / k);
      for (i = 2; i < l; i++) y = mulii(y, powuu(P[i], E[i] / k));
      x = equali1(x)? y: mulii(x,y);
    }
    *pty = x;
  }
  return k == 1? 0: k;
}

long
Z_isanypower(GEN x, GEN *pty)
{
  pari_sp av = avma;
  long k = Z_isanypower_aux(x, pty);
  if (!k) return gc_long(av,0);
  if (signe(x) < 0)
  {
    long v = vals(k);
    if (v)
    {
      k >>= v;
      if (k == 1) return gc_long(av,0);
      if (!pty) return gc_long(av,k);
      *pty = gerepileuptoint(av, powiu(*pty, 1<<v));
      togglesign(*pty); return k;
    }
    if (pty) togglesign_safe(pty);
  }
  if (!pty) return gc_long(av, k);
  *pty = gerepilecopy(av, *pty); return k;
}

/* Faster than expi(n) == vali(n) or hamming(n) == 1 even for single-word
 * values. If all you have is a word, you can just use n & !(n & (n-1)). */
long
Z_ispow2(GEN n)
{
  GEN xp;
  long i, l;
  ulong u;
  if (signe(n) != 1) return 0;
  xp = int_LSW(n); u = *xp; l = lgefint(n);
  for (i = 3; i < l; ++i)
  {
    if (u) return 0;
    xp = int_nextW(xp); u = *xp;
  }
  return !(u & (u-1));
}

static long
isprimepower_i(GEN n, GEN *pt, long flag)
{
  pari_sp av = avma;
  long i, v;

  if (typ(n) != t_INT) pari_err_TYPE("isprimepower", n);
  if (signe(n) <= 0) return 0;

  if (lgefint(n) == 3)
  {
    ulong p;
    v = uisprimepower(n[2], &p);
    if (v)
    {
      if (pt) *pt = utoipos(p);
      return v;
    }
    return 0;
  }
  for (i = 0; i < 26; i++)
  {
    ulong p = tinyprimes[i];
    v = Z_lvalrem(n, p, &n);
    if (v)
    {
      set_avma(av);
      if (!is_pm1(n)) return 0;
      if (pt) *pt = utoipos(p);
      return v;
    }
  }
  /* p | n => p >= 103 */
  v = Z_isanypower_nosmalldiv(&n); /* expensive */
  if (!(flag? isprime(n): BPSW_psp(n))) return gc_long(av,0);
  if (pt) *pt = gerepilecopy(av, n); else set_avma(av);
  return v;
}
long
isprimepower(GEN n, GEN *pt) { return isprimepower_i(n,pt,1); }
long
ispseudoprimepower(GEN n, GEN *pt) { return isprimepower_i(n,pt,0); }

long
uisprimepower(ulong n, ulong *pp)
{ /* We must have CUTOFF^11 >= ULONG_MAX and CUTOFF^3 < ULONG_MAX.
   * Tests suggest that 200-300 is the best range for 64-bit platforms. */
  const ulong CUTOFF = 200UL;
  const long TINYCUTOFF = 46;  /* tinyprimes[45] = 199 */
  const ulong CUTOFF3 = CUTOFF*CUTOFF*CUTOFF;
#ifdef LONG_IS_64BIT
  /* primes preceeding the appropriate root of ULONG_MAX. */
  const ulong ROOT9 = 137;
  const ulong ROOT8 = 251;
  const ulong ROOT7 = 563;
  const ulong ROOT5 = 7129;
  const ulong ROOT4 = 65521;
#else
  const ulong ROOT9 = 11;
  const ulong ROOT8 = 13;
  const ulong ROOT7 = 23;
  const ulong ROOT5 = 83;
  const ulong ROOT4 = 251;
#endif
  ulong mask;
  long v, i;
  int e;
  if (n < 2) return 0;
  if (!odd(n)) {
    if (n & (n-1)) return 0;
    *pp = 2; return vals(n);
  }
  if (n < 8) { *pp = n; return 1; } /* 3,5,7 */
  for (i = 1/*skip p=2*/; i < TINYCUTOFF; i++)
  {
    ulong p = tinyprimes[i];
    if (n % p == 0)
    {
      v = u_lvalrem(n, p, &n);
      if (n == 1) { *pp = p; return v; }
      return 0;
    }
  }
  /* p | n => p >= CUTOFF */

  if (n < CUTOFF3)
  {
    if (n < CUTOFF*CUTOFF || uisprime_101(n)) { *pp = n; return 1; }
    if (uissquareall(n, &n)) { *pp = n; return 2; }
    return 0;
  }

  /* Check for squares, fourth powers, and eighth powers as appropriate. */
  v = 1;
  if (uissquareall(n, &n)) {
    v <<= 1;
    if (CUTOFF <= ROOT4 && uissquareall(n, &n)) {
      v <<= 1;
      if (CUTOFF <= ROOT8 && uissquareall(n, &n)) v <<= 1;
    }
  }

  if (CUTOFF > ROOT5) mask = 1;
  else
  {
    const ulong CUTOFF5 = CUTOFF3*CUTOFF*CUTOFF;
    if (n < CUTOFF5) mask = 1; else mask = 3;
    if (CUTOFF <= ROOT7)
    {
      const ulong CUTOFF7 = CUTOFF5*CUTOFF*CUTOFF;
      if (n >= CUTOFF7) mask = 7;
    }
  }

  if (CUTOFF <= ROOT9 && (e = uis_357_power(n, &n, &mask))) { v *= e; mask=1; }
  if ((e = uis_357_power(n, &n, &mask))) v *= e;

  if (uisprime_101(n)) { *pp = n; return v; }
  return 0;
}

