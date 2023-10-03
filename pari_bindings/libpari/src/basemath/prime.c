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

#define DEBUGLEVEL DEBUGLEVEL_isprime

/*********************************************************************/
/**                                                                 **/
/**               PSEUDO PRIMALITY (MILLER-RABIN)                   **/
/**                                                                 **/
/*********************************************************************/
typedef struct {
  GEN n, sqrt1, sqrt2, t1, t;
  long r1;
} MR_Jaeschke_t;

static void
init_MR_Jaeschke(MR_Jaeschke_t *S, GEN n)
{
  S->n = n = absi_shallow(n);
  S->t = subiu(n,1);
  S->r1 = vali(S->t);
  S->t1 = shifti(S->t, -S->r1);
  S->sqrt1 = cgeti(lg(n)); S->sqrt1[1] = evalsigne(0)|evallgefint(2);
  S->sqrt2 = cgeti(lg(n)); S->sqrt2[1] = evalsigne(0)|evallgefint(2);
}

/* is n strong pseudo-prime for base a ? 'End matching' (check for square
 * roots of -1): if ends do mismatch, then we have factored n, and this
 * information should be made available to the factoring machinery. But so
 * exceedingly rare... besides we use BSPW now. */
static int
ispsp(MR_Jaeschke_t *S, ulong a)
{
  pari_sp av = avma;
  GEN c = Fp_pow(utoipos(a), S->t1, S->n);
  long r;

  if (is_pm1(c) || equalii(S->t, c)) return 1;
  /* go fishing for -1, not for 1 (saves one squaring) */
  for (r = S->r1 - 1; r; r--) /* r1 - 1 squarings */
  {
    GEN c2 = c;
    c = remii(sqri(c), S->n);
    if (equalii(S->t, c))
    {
      if (!signe(S->sqrt1))
      {
        affii(subii(S->n, c2), S->sqrt2);
        affii(c2, S->sqrt1); return 1;
      }
      /* saw one earlier: too many sqrt(-1)s mod n ? */
      return equalii(c2, S->sqrt1) || equalii(c2, S->sqrt2);
    }
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ispsp, r = %ld", r);
      c = gerepileuptoint(av, c);
    }
  }
  return 0;
}

/* is n > 0 strong pseudo-prime for base 2 ? Only used when lgefint(n) > 3,
 * so don't test */
static int
is2psp(GEN n)
{
  GEN c, t = subiu(n, 1);
  pari_sp av = avma;
  long e = vali(t);

  c = Fp_pow(gen_2, shifti(t, -e), n);
  if (is_pm1(c) || equalii(t, c)) return 1;
  while (--e)
  { /* go fishing for -1, not for 1 (e - 1 squaring) */
    c = remii(sqri(c), n);
    if (equalii(t, c)) return 1;
    /* can return 0 if (c == 1) but very infrequent */
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"is2psp, r = %ld", e);
      c = gerepileuptoint(av, c);
    }
  }
  return 0;
}
static int
uispsp_pre(ulong a, ulong n, ulong ni)
{
  ulong c, t = n - 1;
  long e = vals(t);

  c = Fl_powu_pre(a, t >> e, n, ni);
  if (c == 1 || c == t) return 1;
  while (--e)
  { /* go fishing for -1, not for 1 (saves one squaring) */
    c = Fl_sqr_pre(c, n, ni);
    if (c == t) return 1;
    /* can return 0 if (c == 1) but very infrequent */
  }
  return 0;
}
int
uispsp(ulong a, ulong n)
{
  ulong c, t;
  long e;

  if (n & HIGHMASK) return uispsp_pre(a, n, get_Fl_red(n));
  t = n - 1;
  e = vals(t);
  c = Fl_powu(a, t >> e, n);
  if (c == 1 || c == t) return 1;
  while (--e)
  { /* go fishing for -1, not for 1 (e - 1 squaring) */
    c = Fl_sqr(c, n);
    if (c == t) return 1;
    /* can return 0 if (c == 1) but very infrequent */
  }
  return 0;
}
int
uis2psp(ulong n) { return uispsp(2, n); }

/* Miller-Rabin test for k random bases */
long
millerrabin(GEN n, long k)
{
  pari_sp av2, av = avma;
  ulong r;
  long i;
  MR_Jaeschke_t S;

  if (typ(n) != t_INT) pari_err_TYPE("millerrabin",n);
  if (signe(n) <= 0) return 0;
  /* If |n| <= 3, check if n = +- 1 */
  if (lgefint(n) == 3 && uel(n,2) <= 3) return uel(n,2) != 1;

  if (!mod2(n)) return 0;
  init_MR_Jaeschke(&S, n); av2 = avma;
  for (i = 1; i <= k; i++)
  {
    do r = umodui(pari_rand(), n); while (!r);
    if (DEBUGLEVEL > 4) err_printf("Miller-Rabin: testing base %ld\n", r);
    if (!ispsp(&S, r)) return gc_long(av,0);
    set_avma(av2);
  }
  return gc_long(av,1);
}

GEN
gispseudoprime(GEN x, long flag)
{ return flag? map_proto_lGL(millerrabin, x, flag): map_proto_lG(BPSW_psp,x); }

long
ispseudoprime(GEN x, long flag)
{ return flag? millerrabin(x, flag): BPSW_psp(x); }

int
MR_Jaeschke(GEN n)
{
  MR_Jaeschke_t S;
  pari_sp av;

  if (lgefint(n) == 3) return uisprime(uel(n,2));
  if (!mod2(n)) return 0;
  av = avma; init_MR_Jaeschke(&S, n);
  return gc_int(av, ispsp(&S, 31) && ispsp(&S, 73));
}

/*********************************************************************/
/**                                                                 **/
/**                      PSEUDO PRIMALITY (LUCAS)                   **/
/**                                                                 **/
/*********************************************************************/
/* compute n-th term of Lucas sequence modulo N.
 * v_{k+2} = P v_{k+1} - v_k, v_0 = 2, v_1 = P.
 * Assume n > 0 */
static GEN
LucasMod(GEN n, ulong P, GEN N)
{
  pari_sp av = avma;
  GEN nd = int_MSW(n);
  ulong m = *nd;
  long i, j;
  GEN v = utoipos(P), v1 = utoipos(P*P - 2);

  if (m == 1)
    j = 0;
  else
  {
    j = 1+bfffo(m); /* < BIL */
    m <<= j; j = BITS_IN_LONG - j;
  }
  for (i=lgefint(n)-2;;) /* cf. leftright_pow */
  {
    for (; j; m<<=1,j--)
    { /* v = v_k, v1 = v_{k+1} */
      if (m&HIGHBIT)
      { /* set v = v_{2k+1}, v1 = v_{2k+2} */
        v = subiu(mulii(v,v1), P);
        v1= subiu(sqri(v1), 2);
      }
      else
      {/* set v = v_{2k}, v1 = v_{2k+1} */
        v1= subiu(mulii(v,v1), P);
        v = subiu(sqri(v), 2);
      }
      v = modii(v, N);
      v1= modii(v1,N);
      if (gc_needed(av,1))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"LucasMod");
        gerepileall(av, 2, &v,&v1);
      }
    }
    if (--i == 0) return v;
    j = BITS_IN_LONG;
    nd=int_precW(nd);
    m = *nd;
  }
}
/* compute n-th term of Lucas sequence modulo N.
 * v_{k+2} = P v_{k+1} - v_k, v_0 = 2, v_1 = P.
 * Assume n > 0 */
static ulong
u_LucasMod_pre(ulong n, ulong P, ulong N, ulong NI)
{
  ulong v, v1, m;
  long j;

  if (n == 1) return P;
  j = 1 + bfffo(n); /* < BIL */
  v = P; v1 = P*P - 2;
  m = n<<j; j = BITS_IN_LONG - j;
  for (; j; m<<=1,j--)
  { /* v = v_k, v1 = v_{k+1} */
    if (m & HIGHBIT)
    { /* set v = v_{2k+1}, v1 = v_{2k+2} */
      v = Fl_sub(Fl_mul_pre(v,v1,N,NI), P, N);
      v1= Fl_sub(Fl_sqr_pre(v1,N,NI), 2UL, N);
    }
    else
    {/* set v = v_{2k}, v1 = v_{2k+1} */
      v1= Fl_sub(Fl_mul_pre(v,v1,N,NI),P, N);
      v = Fl_sub(Fl_sqr_pre(v,N,NI), 2UL, N);
    }
  }
  return v;
}

static ulong
get_disc(ulong n)
{
  ulong b;
  long i;
  for (b = 3, i = 0;; b += 2, i++)
  {
    ulong c = b*b - 4; /* = 1 mod 4 */
    if (krouu(n % c, c) < 0) break;
    if (i == 64 && uissquareall(n, &c)) return 0; /* oo loop if N = m^2 */
  }
  return b;
}
static int
uislucaspsp_pre(ulong n, ulong ni)
{
  long i, v;
  ulong b, z, m = n + 1;

  if (!m) return 0; /* neither 2^32-1 nor 2^64-1 are Lucas-pp */
  b = get_disc(n); if (!b) return 0;
  v = vals(m); m >>= v;
  z = u_LucasMod_pre(m, b, n, ni);
  if (z == 2 || z == n-2) return 1;
  for (i=1; i<v; i++)
  {
    if (!z) return 1;
    z = Fl_sub(Fl_sqr_pre(z,n,ni), 2UL, n);
    if (z == 2) return 0;
  }
  return 0;
}
/* never called; no need to optimize */
int
uislucaspsp(ulong n)
{ return uislucaspsp_pre(n, get_Fl_red(n)); }

/* N > 3. Caller should check that N is not a square first (taken care of here,
 * but inefficient) */
static int
islucaspsp(GEN N)
{
  pari_sp av = avma;
  GEN m, z;
  long i, v;
  ulong b;

  for (b=3;; b+=2)
  {
    ulong c = b*b - 4; /* = 1 mod 4 */
    if (b == 129 && Z_issquare(N)) return 0; /* avoid oo loop if N = m^2 */
    if (krouu(umodiu(N,c), c) < 0) break;
  }
  m = addiu(N,1); v = vali(m); m = shifti(m,-v);
  z = LucasMod(m, b, N);
  if (absequaliu(z, 2)) return 1;
  if (equalii(z, subiu(N,2))) return 1;
  for (i=1; i<v; i++)
  {
    if (!signe(z)) return 1;
    z = modii(subiu(sqri(z), 2), N);
    if (absequaliu(z, 2)) return 0;
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"islucaspsp");
      z = gerepileupto(av, z);
    }
  }
  return 0;
}

/* composite strong 2-pseudoprime < 1016801 whose prime divisors are > 101.
 * All have a prime divisor <= 661 */
static int
is_2_prp_101(ulong n)
{
  switch(n) {
  case 42799:
  case 49141:
  case 88357:
  case 90751:
  case 104653:
  case 130561:
  case 196093:
  case 220729:
  case 253241:
  case 256999:
  case 271951:
  case 280601:
  case 357761:
  case 390937:
  case 458989:
  case 486737:
  case 489997:
  case 514447:
  case 580337:
  case 741751:
  case 838861:
  case 873181:
  case 877099:
  case 916327:
  case 976873:
  case 983401: return 1;
  } return 0;
}

static int
_uispsp(ulong a, long n) { a %= n; return !a || uispsp(a, n); }
static int
_uisprime(ulong n)
{
#ifdef LONG_IS_64BIT
  ulong ni;
  if (n < 341531)
    return _uispsp(9345883071009581737UL, n);
  if (n < 1050535501)
    return _uispsp(336781006125UL, n)
       && _uispsp(9639812373923155UL, n);
  if (n < 350269456337)
    return _uispsp(4230279247111683200UL, n)
        && _uispsp(14694767155120705706UL, n)
        && _uispsp(16641139526367750375UL, n);
  /* n & HIGHMASK */
  ni = get_Fl_red(n);
  return (uispsp_pre(2, n, ni) && uislucaspsp_pre(n,ni));
#else
  if (n < 360018361) return _uispsp(1143370UL, n) && _uispsp(2350307676UL, n);
  return uispsp(15, n) && uispsp(176006322UL, n) && _uispsp(4221622697UL, n);
#endif
}

int
uisprime(ulong n)
{
  if (n < 103)
    switch(n)
    {
      case 2:
      case 3:
      case 5:
      case 7:
      case 11:
      case 13:
      case 17:
      case 19:
      case 23:
      case 29:
      case 31:
      case 37:
      case 41:
      case 43:
      case 47:
      case 53:
      case 59:
      case 61:
      case 67:
      case 71:
      case 73:
      case 79:
      case 83:
      case 89:
      case 97:
      case 101: return 1;
      default: return 0;
    }
  /* gcd-extraction is much slower */
  return odd(n) && n % 3 && n % 5 && n % 7 && n % 11 && n % 13 && n % 17
                && n % 19 && n % 23 && n % 29 && n % 31 && n % 37 && n % 41
                && (n < 1849 || _uisprime(n));
}

/* assume no prime divisor <= 101 */
int
uisprime_101(ulong n)
{
  if (n < 1016801) return n < 10609? 1: (uispsp(2, n) && !is_2_prp_101(n));
  return _uisprime(n);
}

/* assume no prime divisor <= 661 */
int
uisprime_661(ulong n)
{
  if (n < 1016801) return n < 452929? 1: uispsp(2, n);
  return _uisprime(n);
}

static int
iu_coprime(GEN N, ulong u) { return ugcd(u, umodiu(N,u)) == 1; }
long
BPSW_psp(GEN N)
{
  pari_sp av;
  if (typ(N) != t_INT) pari_err_TYPE("BPSW_psp",N);
  if (signe(N) <= 0) return 0;
  if (lgefint(N) == 3) return uisprime(uel(N,2));
  if (!mod2(N)) return 0;
#ifdef LONG_IS_64BIT
  /* 16294579238595022365 = 3*5*7*11*13*17*19*23*29*31*37*41*43*47*53
   *  7145393598349078859 = 59*61*67*71*73*79*83*89*97*101 */
  if (!iu_coprime(N, 16294579238595022365UL) ||
      !iu_coprime(N,  7145393598349078859UL)) return 0;
#else
  /* 4127218095 = 3*5*7*11*13*17*19*23*37
   * 3948078067 = 29*31*41*43*47*53
   * 4269855901 = 59*83*89*97*101
   * 1673450759 = 61*67*71*73*79 */
  if (!iu_coprime(N, 4127218095UL) ||
      !iu_coprime(N, 3948078067UL) ||
      !iu_coprime(N, 1673450759UL) ||
      !iu_coprime(N, 4269855901UL)) return 0;
#endif
  /* no prime divisor < 103 */
  av = avma;
  return gc_long(av, is2psp(N) && islucaspsp(N));
}

/* can we write n = x^k ? Assume N has no prime divisor <= 2^14, else may
 * miss som powers. Not memory clean */
long
isanypower_nosmalldiv(GEN N, GEN *px)
{
  GEN x = N, y;
  ulong mask = 7;
  long ex, k = 1;
  forprime_t T;
  while (Z_issquareall(x, &y)) { k <<= 1; x = y; }
  while ( (ex = is_357_power(x, &y, &mask)) ) { k *= ex; x = y; }
  (void)u_forprime_init(&T, 11, ULONG_MAX);
  /* stop when x^(1/k) < 2^14 */
  while ( (ex = is_pth_power(x, &y, &T, 15)) ) { k *= ex; x = y; }
  *px = x; return k;
}

/* no prime divisor <= 2^14 (> 661) */
long
BPSW_psp_nosmalldiv(GEN N)
{
  pari_sp av;
  long l = lgefint(N);

  if (l == 3) return uisprime_661(uel(N,2));
  av = avma;
  /* N large: test for pure power, rarely succeeds, but requires < 1% of
   * compositeness test times */
  if (bit_accuracy(l) > 512 && isanypower_nosmalldiv(N,&N) != 1)
    return gc_long(av,0);
  N = absi_shallow(N);
  return gc_long(av, is2psp(N) && islucaspsp(N));
}

/***********************************************************************/
/**                                                                   **/
/**                       Pocklington-Lehmer                          **/
/**                        P-1 primality test                         **/
/**                                                                   **/
/***********************************************************************/
/* Assume x BPSW pseudoprime. Check whether it's small enough to be certified
 * prime (< 2^64). Reference for strong 2-pseudoprimes:
 *   http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html */
static int
BPSW_isprime_small(GEN x)
{
  long l = lgefint(x);
#ifdef LONG_IS_64BIT
  return (l == 3);
#else
  return (l <= 4);
#endif
}

/* Assume N > 1, p^e || N-1, p prime. Find a witness a(p) such that
 *   a^(N-1) = 1 (mod N)
 *   a^(N-1)/p - 1 invertible mod N.
 * Proves that any divisor of N is 1 mod p^e. Return 0 if N is composite */
static ulong
pl831(GEN N, GEN p)
{
  GEN b, c, g, Nmunp = diviiexact(subiu(N,1), p);
  pari_sp av = avma;
  ulong a;
  for(a = 2;; a++, set_avma(av))
  {
    b = Fp_pow(utoipos(a), Nmunp, N);
    if (!equali1(b)) break;
  }
  c = Fp_pow(b,p,N);
  g = gcdii(subiu(b,1), N); /* 0 < g < N */
  return (equali1(c) && equali1(g))? a: 0;
}

/* Brillhart, Lehmer, Selfridge test (Crandall & Pomerance, Th 4.1.5)
 * N^(1/3) <= f fully factored, f | N-1. If pl831(p) is true for
 * any prime divisor p of f, then any divisor of N is 1 mod f.
 * In that case return 1 iff N is prime */
static int
BLS_test(GEN N, GEN f)
{
  GEN c1, c2, r, q;
  q = dvmdii(N, f, &r);
  if (!is_pm1(r)) return 0;
  c2 = dvmdii(q, f, &c1);
  /* N = 1 + f c1 + f^2 c2, 0 <= c_i < f; check whether it is of the form
   * (1 + fa)(1 + fb) */
  return ! Z_issquare(subii(sqri(c1), shifti(c2,2)));
}

/* BPSW_psp(N) && !BPSW_isprime_small(N). Decide between Pocklington-Lehmer
 * and APRCL/ECPP. Return a vector of (small) primes such that PL-witnesses
 * guarantee the primality of N. Return NULL if PL is likely too expensive.
 * Return gen_0 if BLS test finds N to be composite */
static GEN
BPSW_try_PL(GEN N)
{
  ulong B = minuu(1UL<<19, maxprime());
  GEN E, p, U, F, N_1 = subiu(N,1);
  GEN fa = absZ_factor_limit_strict(N_1, B, &U), P = gel(fa,1);

  if (!U) return P; /* N-1 fully factored */
  p = gel(U,1);
  E = gel(fa,2);
  U = powii(p, gel(U,2)); /* unfactored part of N-1 */
  F = (lg(P) == 2)? powii(gel(P,1), gel(E,1)): diviiexact(N_1,  U);

  /* N-1 = F U, F factored, U composite, (U,F) = 1 */
  if (cmpii(F, U) > 0) return P; /* 1/2-smooth */
  if (cmpii(sqri(F), U) > 0) return BLS_test(N,F)? P: gen_0; /* 1/3-smooth */
  return NULL; /* not smooth enough */
}

static GEN isprimePL(GEN N);
/* F is a vector whose entries are primes. For each of them, find a PL
 * witness. Return 0 if caller lied and F contains a composite */
static long
PL_certify(GEN N, GEN F)
{
  long i, l = lg(F);
  for(i = 1; i < l; i++)
    if (! pl831(N, gel(F,i))) return 0;
  return 1;
}
/* F is a vector whose entries are *believed* to be primes (BPSW_psp).
 * For each of them, recording a witness and recursive primality certificate */
static GEN
PL_certificate(GEN N, GEN F)
{
  long i, l = lg(F);
  GEN C;
  if (BPSW_isprime_small(N)) return N;
  C = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(F,i), C0;
    ulong w;

    if (BPSW_isprime_small(p)) { gel(C,i) = p; continue; }
    w = pl831(N,p); if (!w) return gen_0;
    C0 = isprimePL(p);
    if (isintzero(C0))
    { /* composite in prime factorisation ! */
      err_printf("Not a prime: %Ps", p);
      pari_err_BUG("PL_certificate [false prime number]");
    }
    gel(C,i) = mkvec3(p,utoipos(w), C0);
  }
  return mkvec2(N, C);
}
/* M a t_MAT */
static int
PL_isvalid(GEN v)
{
  GEN C, F, F2, N, N1, U;
  long i, l;
  switch(typ(v))
  {
    case t_INT: return BPSW_isprime_small(v) && BPSW_psp(v);
    case t_VEC: if (lg(v) == 3) break;/*fall through */
    default: return 0;
  }
  N = gel(v,1);
  C = gel(v,2);
  if (typ(N) != t_INT || signe(N) <= 0 || typ(C) != t_VEC) return 0;
  U = N1 = subiu(N,1);
  l = lg(C);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(C,i), a = NULL, C0 = NULL, ap;
    long vp;
    if (typ(p) != t_INT)
    {
      if (lg(p) != 4) return 0;
      a = gel(p,2); C0 = gel(p,3); p = gel(p,1);
      if (typ(p) != t_INT || typ(a) != t_INT || !PL_isvalid(C0)) return 0;
    }
    vp = Z_pvalrem(U, p, &U); if (!vp) return 0;
    if (!a)
    {
      if (!BPSW_isprime_small(p)) return 0;
      continue;
    }
    if (!equalii(gel(C0,1), p)) return 0;
    ap = Fp_pow(a, diviiexact(N1,p), N);
    if (!equali1(gcdii(subiu(ap,1), N)) || !equali1(Fp_pow(ap, p, N))) return 0;
  }
  F = diviiexact(N1, U); /* factored part of N-1 */
  F2= sqri(F);
  if (cmpii(F2, N) > 0) return 1;
  if (cmpii(mulii(F2,F), N) <= 0) return 0;
  return BLS_test(N,F);
}

/* Assume N is a strong BPSW pseudoprime, Pocklington-Lehmer primality proof.
 * Return gen_0 (nonprime), N (small prime), matrix (large prime)
 *
 * The matrix has 3 columns, [a,b,c] with
 * a[i] prime factor of N-1,
 * b[i] witness for a[i] as in pl831
 * c[i] check_prime(a[i]) */
static GEN
isprimePL(GEN N)
{
  GEN cbrtN, N_1, F, f;
  if (BPSW_isprime_small(N)) return N;
  cbrtN = sqrtnint(N,3);
  N_1 = subiu(N,1);
  F = Z_factor_until(N_1, sqri(cbrtN));
  f = factorback(F); /* factored part of N-1, f^3 > N */
  if (DEBUGLEVEL>3)
  {
    GEN r = divri(itor(f,LOWDEFAULTPREC), N);
    err_printf("Pocklington-Lehmer: proving primality of N = %Ps\n", N);
    err_printf("Pocklington-Lehmer: N-1 factored up to %Ps! (%.3Ps%%)\n", f, r);
  }
  /* if N-1 is only N^(1/3)-smooth, BLS test */
  if (!equalii(f,N_1) && cmpii(sqri(f),N) <= 0 && !BLS_test(N,f))
    return gen_0; /* Failed, N is composite */
  F = gel(F,1); settyp(F, t_VEC);
  return PL_certificate(N, F);
}

/* assume N a BPSW pseudoprime, in particular, it is odd > 2. Prove N prime */
long
BPSW_isprime(GEN N)
{
  pari_sp av;
  long t;
  GEN P;
  if (BPSW_isprime_small(N)) return 1;
  av = avma; P = BPSW_try_PL(N);
  if (!P) /* not smooth enough */
    t = expi(N) < 768? isprimeAPRCL(N): isprimeECPP(N);
  else
    t = (typ(P) == t_INT)? 0: PL_certify(N,P);
  return gc_long(av,t);
}

static long
_isprimePL(GEN x)
{
  pari_sp av = avma;
  if (!BPSW_psp(x)) return 0;
  return gc_long(av, !isintzero(isprimePL(x)));
}
GEN
gisprime(GEN x, long flag)
{
  switch (flag)
  {
    case 0: return map_proto_lG(isprime,x);
    case 1: return map_proto_lG(_isprimePL,x);
    case 2: return map_proto_lG(isprimeAPRCL,x);
    case 3: return map_proto_lG(isprimeECPP,x);
  }
  pari_err_FLAG("gisprime");
  return NULL;/*LCOV_EXCL_LINE*/
}

long
isprime(GEN x) { return BPSW_psp(x) && BPSW_isprime(x); }

GEN
primecert0(GEN x, long flag, long stopat)
{
  if ((flag || typ(x) == t_INT) && !BPSW_psp(x)) return gen_0;
  switch(flag)
  {
    case 0: return ecpp0(x, stopat);
    case 1: { pari_sp av = avma; return gerepilecopy(av, isprimePL(x)); }
  }
  pari_err_FLAG("primecert");
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
primecert(GEN x, long flag)
{ return primecert0(x, flag, 0); }

enum { c_VOID = 0, c_ECPP, c_N1 };
static long
cert_type(GEN c)
{
  switch(typ(c))
  {
    case t_INT: return c_ECPP;
    case t_VEC:
      if (lg(c) == 3 && typ(gel(c,1)) == t_INT) return c_N1;
      return c_ECPP;
  }
  return c_VOID;
}

long
primecertisvalid(GEN c)
{
  switch(typ(c))
  {
    case t_INT: return BPSW_isprime_small(c) && BPSW_psp(c);
    case t_VEC:
      return cert_type(c) == c_ECPP? ecppisvalid(c): PL_isvalid(c);
  }
  return 0;
}

static long
check_ecppcertentry(GEN c)
{
  GEN v;
  long i,l = lg(c);
  if (typ(c)!=t_VEC || l!=6) return 0;
  for(i=1; i<=4; i++)
    if (typ(gel(c,i))!=t_INT) return 0;
  v = gel(c,5);
  if(typ(v)!=t_VEC) return 0;
  for(i=1; i<=2; i++)
    if (typ(gel(v,i))!=t_INT) return 0;
  return 1;
}

long
check_ecppcert(GEN c)
{
  long i, l;
  switch(typ(c))
  {
    case t_INT: return signe(c) >= 0;
    case t_VEC: break;
    default: return 0;
  }
  l = lg(c); if (l == 1) return 0;
  for (i = 1; i < l; i++)
    if (check_ecppcertentry(gel(c,i))==0) return 0;
  return 1;
}

GEN
primecertexport(GEN c, long flag)
{
  if (cert_type(c) != c_ECPP) pari_err_IMPL("N-1 certificate");
  if (!check_ecppcert(c))
    pari_err_TYPE("primecertexport - invalid certificate", c);
  return ecppexport(c, flag);
}

/***********************************************************************/
/**                                                                   **/
/**                          PRIME NUMBERS                            **/
/**                                                                   **/
/***********************************************************************/

static struct {
  ulong p;
  long n;
} prime_table[] = {
  {           0,          0},
  {        7919,       1000},
  {       17389,       2000},
  {       27449,       3000},
  {       37813,       4000},
  {       48611,       5000},
  {       59359,       6000},
  {       70657,       7000},
  {       81799,       8000},
  {       93179,       9000},
  {      104729,      10000},
  {      224737,      20000},
  {      350377,      30000},
  {      479909,      40000},
  {      611953,      50000},
  {      746773,      60000},
  {      882377,      70000},
  {     1020379,      80000},
  {     1159523,      90000},
  {     1299709,     100000},
  {     2750159,     200000},
  {     7368787,     500000},
  {    15485863,    1000000},
  {    32452843,    2000000},
  {    86028121,    5000000},
  {   179424673,   10000000},
  {   373587883,   20000000},
  {   982451653,   50000000},
  {  2038074743,  100000000},
  {  4000000483UL,189961831},
  {  4222234741UL,200000000},
#if BITS_IN_LONG == 64
  { 5336500537UL,   250000000L},
  { 6461335109UL,   300000000L},
  { 7594955549UL,   350000000L},
  { 8736028057UL,   400000000L},
  { 9883692017UL,   450000000L},
  { 11037271757UL,  500000000L},
  { 13359555403UL,  600000000L},
  { 15699342107UL,  700000000L},
  { 18054236957UL,  800000000L},
  { 20422213579UL,  900000000L},
  { 22801763489UL, 1000000000L},
  { 47055833459UL, 2000000000L},
  { 71856445751UL, 3000000000L},
  { 97011687217UL, 4000000000L},
  {122430513841UL, 5000000000L},
  {148059109201UL, 6000000000L},
  {173862636221UL, 7000000000L},
  {200000000507UL, 8007105083L},
  {225898512559UL, 9000000000L},
  {252097800623UL,10000000000L},
  {384489816343UL,15000000000L},
  {518649879439UL,20000000000L},
  {654124187867UL,25000000000L},
  {790645490053UL,30000000000L},
  {928037044463UL,35000000000L},
 {1066173339601UL,40000000000L},
 {1344326694119UL,50000000000L},
 {1624571841097UL,60000000000L},
 {1906555030411UL,70000000000L},
 {2190026988349UL,80000000000L},
 {2474799787573UL,90000000000L},
 {2760727302517UL,100000000000L}
#endif
};
static const int prime_table_len = numberof(prime_table);

/* find prime closest to n in prime_table. */
static long
prime_table_closest_p(ulong n)
{
  long i;
  for (i = 1; i < prime_table_len; i++)
  {
    ulong p = prime_table[i].p;
    if (p > n)
    {
      ulong u = n - prime_table[i-1].p;
      if (p - n > u) i--;
      break;
    }
  }
  if (i == prime_table_len) i = prime_table_len - 1;
  return i;
}

/* return the n-th successor of prime p > 2; n > 0 */
static GEN
prime_successor(ulong p, ulong n)
{
  GEN Q = utoipos(p), P = NULL;
  ulong i;
RESET:
  for(;;)
  {
    forprime_t S;
    GEN L = dbltor(4*n * log(gtodouble(Q) + 1)), R = addii(Q, ceil_safe(L));

    forprime_init(&S, addiu(Q, 1), R);
    Q = R;
    for (i = 1; i <= n; i++)
    {
      P = forprime_next(&S);
      if (!P) { n -= i-1; goto RESET; }
    }
    return P;
  }
}
/* find the N-th prime */
static GEN
prime_table_find_n(ulong N)
{
  byteptr d;
  ulong n, p, maxp = maxprime();
  long i;
  for (i = 1; i < prime_table_len; i++)
  {
    n = prime_table[i].n;
    if (n > N)
    {
      ulong u = N - prime_table[i-1].n;
      if (n - N > u) i--;
      break;
    }
  }
  if (i == prime_table_len) i = prime_table_len - 1;
  p = prime_table[i].p;
  n = prime_table[i].n;
  if (n > N && p > maxp)
  {
    i--;
    p = prime_table[i].p;
    n = prime_table[i].n;
  }
  /* if beyond prime table, then n <= N */
  d = diffptr + n;
  if (n > N)
  {
    n -= N;
    do { n--; PREC_PRIME_VIADIFF(p,d); } while (n) ;
  }
  else if (n < N)
  {
    ulong maxpN = maxprimeN();
    if (N >= maxpN)
    {
      if (N == maxpN) return utoipos(maxp);
      if (p < maxp) { p = maxp; n = maxpN; }
      return prime_successor(p, N - n);
    }
    /* can find prime(N) in table */
    n = N - n;
    do { n--; NEXT_PRIME_VIADIFF(p,d); } while (n) ;
  }
  return utoipos(p);
}

ulong
uprime(long N)
{
  pari_sp av = avma;
  GEN p;
  if (N <= 0) pari_err_DOMAIN("prime", "n", "<=",gen_0, stoi(N));
  p = prime_table_find_n(N);
  if (lgefint(p) != 3) pari_err_OVERFLOW("uprime");
  return gc_ulong(av, p[2]);
}
GEN
prime(long N)
{
  pari_sp av = avma;
  GEN p;
  if (N <= 0) pari_err_DOMAIN("prime", "n", "<=",gen_0, stoi(N));
  new_chunk(4); /*HACK*/
  p = prime_table_find_n(N);
  set_avma(av); return icopy(p);
}

static void
prime_interval(GEN N, GEN *pa, GEN *pb, GEN *pd)
{
  GEN a, b, d;
  switch(typ(N))
  {
    case t_INT:
      a = gen_2;
      b = subiu(N,1); /* between 2 and N-1 */
      d = subiu(N,2);
      if (signe(d) <= 0)
        pari_err_DOMAIN("randomprime","N", "<=", gen_2, N);
      break;
    case t_VEC:
      if (lg(N) != 3) pari_err_TYPE("randomprime",N);
      a = gel(N,1);
      b = gel(N,2);
      if (gcmp(b, a) < 0)
        pari_err_DOMAIN("randomprime","b-a", "<", gen_0, mkvec2(a,b));
      if (typ(a) != t_INT)
      {
        a = gceil(a);
        if (typ(a) != t_INT) pari_err_TYPE("randomprime",a);
      }
      if (typ(b) != t_INT)
      {
        b = gfloor(b);
        if (typ(b) != t_INT) pari_err_TYPE("randomprime",b);
      }
      if (cmpiu(a, 2) < 0)
      {
        a = gen_2;
        d = subiu(b,1);
      }
      else
        d = addiu(subii(b,a), 1);
      if (signe(d) <= 0)
        pari_err_DOMAIN("randomprime","floor(b) - max(ceil(a),2)", "<",
                        gen_0, mkvec2(a,b));
      break;
    default:
      pari_err_TYPE("randomprime", N);
      a = b = d = NULL; /*LCOV_EXCL_LINE*/
  }
  *pa = a; *pb = b; *pd = d;
}

/* random b-bit prime */
GEN
randomprime(GEN N)
{
  pari_sp av = avma, av2;
  GEN a, b, d;
  if (!N)
    for(;;)
    {
      ulong p = random_bits(31);
      if (uisprime(p)) return utoipos(p);
    }
  prime_interval(N, &a, &b, &d); av2 = avma;
  for (;;)
  {
    GEN p = addii(a, randomi(d));
    if (BPSW_psp(p)) return gerepileuptoint(av, p);
    set_avma(av2);
  }
}
GEN
randomprime0(GEN N, GEN q)
{
  pari_sp av = avma, av2;
  GEN C, D, a, b, d, r;
  if (!q) return randomprime(N);
  switch(typ(q))
  {
    case t_INT: C = gen_1; D = q; break;
    case t_INTMOD: C = gel(q,2); D = gel(q,1); break;
    default:
      pari_err_TYPE("randomprime", q);
      return NULL;/*LCOV_EXCL_LINE*/
  }
  if (!N) N = int2n(31);
  prime_interval(N, &a, &b, &d);
  r = modii(subii(C, a), D);
  if (signe(r)) { a = addii(a, r); d = subii(d, r); }
  if (!equali1(gcdii(C,D)))
  {
    if (isprime(a)) return gerepilecopy(av, a);
    pari_err_COPRIME("randomprime", C, D);
  }
  d = divii(d, D); if (!signe(d)) d = gen_1;
  av2 = avma;
  for (;;)
  {
    GEN p = addii(a, mulii(D, randomi(d)));
    if (BPSW_psp(p)) return gerepileuptoint(av, p);
    set_avma(av2);
  }
  return NULL;
}

/* set *pp = nextprime(a) = p
 *     *pd so that NEXT_PRIME_VIADIFF(d, p) = nextprime(p+1)
 *     *pn so that p = the n-th prime
 * error if nextprime(a) is out of primetable bounds */
void
prime_table_next_p(ulong a, byteptr *pd, ulong *pp, ulong *pn)
{
  byteptr d;
  ulong p, n, maxp = maxprime();
  long i = prime_table_closest_p(a);
  p = prime_table[i].p;
  if (p > a && p > maxp)
  {
    i--;
    p = prime_table[i].p;
  }
  /* if beyond prime table, then p <= a */
  n = prime_table[i].n;
  d = diffptr + n;
  if (p < a)
  {
    if (a > maxp) pari_err_MAXPRIME(a);
    do { n++; NEXT_PRIME_VIADIFF(p,d); } while (p < a);
  }
  else if (p != a)
  {
    do { n--; PREC_PRIME_VIADIFF(p,d); } while (p > a) ;
    if (p < a) { NEXT_PRIME_VIADIFF(p,d); n++; }
  }
  *pn = n;
  *pp = p;
  *pd = d;
}

ulong
uprimepi(ulong a)
{
  ulong p, n, maxp = maxprime();
  if (a <= maxp)
  {
    byteptr d;
    prime_table_next_p(a, &d, &p, &n);
    return p == a? n: n-1;
  }
  else
  {
    long i = prime_table_closest_p(a);
    forprime_t S;
    p = prime_table[i].p;
    if (p > a)
    {
      i--;
      p = prime_table[i].p;
    }
    /* p = largest prime in table <= a */
    n = prime_table[i].n;
    (void)u_forprime_init(&S, p+1, a);
    for (; p; n++) p = u_forprime_next(&S);
    return n-1;
  }
}

GEN
primepi(GEN x)
{
  pari_sp av = avma;
  GEN pp, nn, N = typ(x) == t_INT? x: gfloor(x);
  forprime_t S;
  ulong n, p;
  long i;
  if (typ(N) != t_INT) pari_err_TYPE("primepi",N);
  if (signe(N) <= 0) return gen_0;
  if (lgefint(N) == 3) return gc_utoi(av, uprimepi(N[2]));
  i = prime_table_len-1;
  p = prime_table[i].p;
  n = prime_table[i].n;
  (void)forprime_init(&S, utoipos(p+1), N);
  nn = setloop(utoipos(n));
  pp = gen_0;
  for (; pp; incloop(nn)) pp = forprime_next(&S);
  return gerepileuptoint(av, subiu(nn,1));
}

/* pi(x) < x/log x * (1 + 1/log x + 2.51/log^2 x)), x>=355991 [ Dusart ]
 * pi(x) < x/(log x - 1.1), x >= 60184 [ Dusart ]
 * ? \p9
 * ? M = 0; for(x = 4, 60184, M = max(M, log(x) - x/primepi(x))); M
 * %1 = 1.11196252 */
double
primepi_upper_bound(double x)
{
  if (x >= 355991)
  {
    double L = 1/log(x);
    return x * L * (1 + L + 2.51*L*L);
  }
  if (x >= 60184) return x / (log(x) - 1.1);
  if (x < 5) return 2; /* don't bother */
  return x / (log(x) - 1.111963);
}
/* pi(x) > x/log x (1 + 1/log x), x >= 599 [ Dusart ]
 * pi(x) > x / (log x + 2), x >= 55 [ Rosser ] */
double
primepi_lower_bound(double x)
{
  if (x >= 599)
  {
    double L = 1/log(x);
    return x * L * (1 + L);
  }
  if (x < 55) return 0; /* don't bother */
  return x / (log(x) + 2.);
}
GEN
gprimepi_upper_bound(GEN x)
{
  pari_sp av = avma;
  double L;
  if (typ(x) != t_INT) x = gfloor(x);
  if (expi(x) <= 1022)
  {
    set_avma(av);
    return dbltor(primepi_upper_bound(gtodouble(x)));
  }
  x = itor(x, LOWDEFAULTPREC);
  L = 1 / rtodbl(logr_abs(x));
  x = mulrr(x, dbltor(L * (1 + L + 2.51*L*L)));
  return gerepileuptoleaf(av, x);
}
GEN
gprimepi_lower_bound(GEN x)
{
  pari_sp av = avma;
  double L;
  if (typ(x) != t_INT) x = gfloor(x);
  if (abscmpiu(x, 55) <= 0) return gen_0;
  if (expi(x) <= 1022)
  {
    set_avma(av);
    return dbltor(primepi_lower_bound(gtodouble(x)));
  }
  x = itor(x, LOWDEFAULTPREC);
  L = 1 / rtodbl(logr_abs(x));
  x = mulrr(x, dbltor(L * (1 + L)));
  return gerepileuptoleaf(av, x);
}

GEN
primes(long n)
{
  forprime_t S;
  long i;
  GEN y;
  if (n <= 0) return cgetg(1, t_VEC);
  y = cgetg(n+1, t_VEC);
  (void)new_chunk(3*n); /*HACK*/
  u_forprime_init(&S, 2, (ulong)n > maxprimeN()? ULONG_MAX: maxprime());
  set_avma((pari_sp)y);
  for (i = 1; i <= n; i++) gel(y, i) = utoipos( u_forprime_next(&S) );
  return y;
}
GEN
primes_zv(long n)
{
  forprime_t S;
  long i;
  GEN y;
  if (n <= 0) return cgetg(1, t_VECSMALL);
  y = cgetg(n+1,t_VECSMALL);
  /* don't initialize sieve if all fits in primetable */
  u_forprime_init(&S, 2, (ulong)n > maxprimeN()? ULONG_MAX: maxprime());
  for (i = 1; i <= n; i++) y[i] =  u_forprime_next(&S);
  set_avma((pari_sp)y); return y;
}
GEN
primes0(GEN N)
{
  switch(typ(N))
  {
    case t_INT: return primes(itos(N));
    case t_VEC:
      if (lg(N) == 3) return primes_interval(gel(N,1),gel(N,2));
  }
  pari_err_TYPE("primes", N);
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
primes_interval(GEN a, GEN b)
{
  pari_sp av = avma;
  forprime_t S;
  long i, n;
  GEN y, d, p;
  if (typ(a) != t_INT)
  {
    a = gceil(a);
    if (typ(a) != t_INT) pari_err_TYPE("primes_interval",a);
  }
  if (typ(b) != t_INT)
  {
    b = gfloor(b);
    if (typ(b) != t_INT) pari_err_TYPE("primes_interval",b);
  }
  if (signe(a) < 0) a = gen_2;
  d = subii(b, a);
  if (signe(d) < 0 || signe(b) <= 0) { set_avma(av); return cgetg(1, t_VEC); }
  if (lgefint(b) == 3)
  {
    set_avma(av);
    y = primes_interval_zv(itou(a), itou(b));
    n = lg(y); settyp(y, t_VEC);
    for (i = 1; i < n; i++) gel(y,i) = utoipos(y[i]);
    return y;
  }
  /* at most d+1 primes in [a,b]. If d large, try better bound to lower
   * memory use */
  if (abscmpiu(d,100000) > 0)
  {
    GEN D = gsub(gprimepi_upper_bound(b), gprimepi_lower_bound(a));
    D = ceil_safe(D);
    if (cmpii(D, d) < 0) d = D;
  }
  n = itos(d)+1;
  forprime_init(&S, a, b);
  y = cgetg(n+1, t_VEC); i = 1;
  while ((p = forprime_next(&S))) gel(y, i++) = icopy(p);
  setlg(y, i); return gerepileupto(av, y);
}

/* a <= b, at most d primes in [a,b]. Return them */
static GEN
primes_interval_i(ulong a, ulong b, ulong d)
{
  ulong p, i = 1, n = d + 1;
  forprime_t S;
  GEN y = cgetg(n+1, t_VECSMALL);
  pari_sp av = avma;
  u_forprime_init(&S, a, b);
  while ((p = u_forprime_next(&S))) y[i++] = p;
  set_avma(av); setlg(y, i); stackdummy((pari_sp)(y + i), (pari_sp)(y + n+1));
  return y;
}
GEN
primes_interval_zv(ulong a, ulong b)
{
  ulong d;
  if (a < 3) return primes_upto_zv(b);
  if (b < a) return cgetg(1, t_VECSMALL);
  d = b - a;
  if (d > 100000UL)
  {
    ulong D = (ulong)ceil(primepi_upper_bound(b)-primepi_lower_bound(a));
    if (D < d) d = D;
  }
  return primes_interval_i(a, b, d);
}
GEN
primes_upto_zv(ulong b)
{
  ulong d;
  if (b < 2) return cgetg(1, t_VECSMALL);
  d = (b > 100000UL)? (ulong)primepi_upper_bound(b): b;
  return primes_interval_i(2, b, d);
}

/***********************************************************************/
/**                                                                   **/
/**                       PRIVATE PRIME TABLE                         **/
/**                                                                   **/
/***********************************************************************/

void
pari_set_primetab(GEN global_primetab)
{
  if (global_primetab)
  {
    long i, l = lg(global_primetab);
    primetab = cgetg_block(l, t_VEC);
    for (i = 1; i < l; i++)
      gel(primetab,i) = gclone(gel(global_primetab,i));
  } else primetab = cgetg_block(1, t_VEC);
}

/* delete dummy NULL entries */
static void
cleanprimetab(GEN T)
{
  long i,j, l = lg(T);
  for (i = j = 1; i < l; i++)
    if (T[i]) T[j++] = T[i];
  setlg(T,j);
}
/* remove p from T */
static void
rmprime(GEN T, GEN p)
{
  long i;
  if (typ(p) != t_INT) pari_err_TYPE("removeprimes",p);
  i = ZV_search(T, p);
  if (!i)
    pari_err_DOMAIN("removeprimes","prime","not in",
                    strtoGENstr("primetable"), p);
  gunclone(gel(T,i)); gel(T,i) = NULL;
  cleanprimetab(T);
}

/*stolen from ZV_union_shallow() : clone entries from y */
static GEN
addp_union(GEN x, GEN y)
{
  long i, j, k, lx = lg(x), ly = lg(y);
  GEN z = cgetg(lx + ly - 1, t_VEC);
  i = j = k = 1;
  while (i<lx && j<ly)
  {
    int s = cmpii(gel(x,i), gel(y,j));
    if (s < 0)
      gel(z,k++) = gel(x,i++);
    else if (s > 0)
      gel(z,k++) = gclone(gel(y,j++));
    else {
      gel(z,k++) = gel(x,i++);
      j++;
    }
  }
  while (i<lx) gel(z,k++) = gel(x,i++);
  while (j<ly) gel(z,k++) = gclone(gel(y,j++));
  setlg(z, k); return z;
}

/* p is NULL, or a single element or a row vector with "primes" to add to
 * prime table. */
static GEN
addp(GEN *T, GEN p)
{
  pari_sp av = avma;
  long i, l;
  GEN v;

  if (!p || lg(p) == 1) return *T;
  if (!is_vec_t(typ(p))) p = mkvec(p);

  RgV_check_ZV(p, "addprimes");
  v = gen_indexsort_uniq(p, (void*)&cmpii, &cmp_nodata);
  p = vecpermute(p, v);
  if (abscmpiu(gel(p,1), 2) < 0) pari_err_DOMAIN("addprimes", "p", "<", gen_2,p);
  p = addp_union(*T, p);
  l = lg(p);
  if (l != lg(*T))
  {
    GEN old = *T, t = cgetg_block(l, t_VEC);
    for (i = 1; i < l; i++) gel(t,i) = gel(p,i);
    *T = t; gunclone(old);
  }
  set_avma(av); return *T;
}
GEN
addprimes(GEN p) { return addp(&primetab, p); }

static GEN
rmprimes(GEN T, GEN prime)
{
  long i,tx;

  if (!prime) return T;
  tx = typ(prime);
  if (is_vec_t(tx))
  {
    if (prime == T)
    {
      for (i=1; i < lg(prime); i++) gunclone(gel(prime,i));
      setlg(prime, 1);
    }
    else
    {
      for (i=1; i < lg(prime); i++) rmprime(T, gel(prime,i));
    }
    return T;
  }
  rmprime(T, prime); return T;
}
GEN
removeprimes(GEN prime) { return rmprimes(primetab, prime); }
