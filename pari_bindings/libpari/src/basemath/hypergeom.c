/* Copyright (C) 2017  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/********************************************************************/
/**                                                                **/
/**                   HYPERGEOMETRIC FUNCTIONS                     **/
/**                                                                **/
/********************************************************************/

#include "pari.h"
#include "paripriv.h"

static GEN
F10(GEN a, GEN z, long prec)
{ return gpow(gsubsg(1,z), gneg(a), prec); }

static int
isnegint2(GEN a, long *pa)
{
  GEN b;
  if (!gequal0(imag_i(a))) return 0;
  a = real_i(a); if (gsigne(a) > 0) return 0;
  b = ground(a); if (!gequal(a, b)) return 0;
  if (pa) *pa = -itos(b);
  return 1;
}
static int
isnegint(GEN a) { return isnegint2(a, NULL); }
static int
isnegint_approx(GEN a, long bit)
{
  GEN b;
  if (gexpo(imag_i(a)) > -bit) return 0;
  a = real_i(a); if (gsigne(a) > 0) return 0;
  b = ground(a); return gexpo(gsub(a, b)) < -bit;
}
static int
is0(GEN a, long bit) { return gequal0(a) || gexpo(a) < -bit; }
static int
islong(GEN a, long *m, long prec)
{
  *m = itos(ground(real_i(a)));
  if (is0(gsubgs(a, *m), prec2nbits(prec) - 5)) return 1;
  return 0;
}

static GEN
F01(GEN a, GEN z, long prec)
{
  GEN A, B, al, sz;
  if (is0(z, prec2nbits(prec)-5)) return real_1(prec);
  sz = gsqrt(z, prec); al = gsubgs(a, 1);
  A = gmul(ggamma(a, prec), gpow(sz, gneg(al), prec));
  B = ibessel(al, gmul2n(sz,1), prec);
  return isexactzero(imag_i(z))? mulreal(A,B): gmul(A,B);
}

/* Airy functions [Ai,Bi] */
static GEN
airy_i(GEN x, long prec)
{
  long bit = prec2nbits(prec), tx = typ(x), prec2;
  GEN a, b, A, B, z, z2;
  if (!is_scalar_t(tx)) pari_err_TYPE("airy",x);
  if (is0(x, bit))
  {
    GEN s = sqrtnr_abs(utor(3,prec), 6), s3 = powrs(s,3), s4 = mulrr(s,s3);
    A = invr(mulrr(s4, ggamma(uutoQ(2,3), prec)));
    B = mulrr(A, s3); return mkvec2(A, B);
  }
  prec2 = prec + EXTRAPREC64;
  x = gprec_wensure(x, prec2);
  z = gsqrt(gpowgs(x,3), prec2); z2 = gdivgu(gmul2n(z,1),3);
  if (is_real_t(tx) && gsigne(x) > 0)
    a = b = gsqrt(x, prec2); /* expression simplifies */
  else
  {
    a = gsqrtn(z, utoipos(3), NULL, prec2);
    b = gdiv(x, a);
  }
  a = gmul(a, ibessel(sstoQ(-1,3),z2, prec));
  b = gmul(b, ibessel(uutoQ(1,3), z2, prec));
  if (isexactzero(imag_i(x))) { a = real_i(a); b = real_i(b); }
  A = gdivgu(gsub(a,b), 3);
  B = gdiv(gadd(a,b), sqrtr_abs(utor(3, prec)));

  bit -= gexpo(a) + 16;
  if (!is0(A, bit) && !is0(B, bit)) return mkvec2(A, B);
  prec = precdbl(prec); x = gprec_wensure(x, prec); return airy_i(x, prec);
}
GEN
airy(GEN z, long prec)
{ pari_sp av = avma; return gerepilecopy(av, airy_i(z, prec)); }

/* Gamma(a)*Gamma(b) */
static GEN
mulgamma2(GEN a, GEN b, long prec)
{ return gmul(ggamma(a, prec), ggamma(b, prec)); }
/* Gamma(a)/Gamma(b) */
static GEN
divgamma2(GEN a, GEN b, long prec)
{ return gdiv(ggamma(a, prec), ggamma(b, prec)); }
/* Gamma(v[1])*Gamma(v[2]) */
static GEN
mulgammav2(GEN v, long prec)
{ return mulgamma2(gel(v,1), gel(v,2), prec); }

/***********************************************************************/
/**                 CONFLUENT HYPERGEOMETRIC U(a,b,z)                 **/
/***********************************************************************/
static GEN Ftaylor(GEN N, GEN D, GEN z, long prec);
/* b not integral; use 1F1; prec0 is precision we really want */
static GEN
hyperu_F11(GEN a, GEN b, GEN z, long prec0, long prec)
{
  GEN S1, S2, b1 = gsubsg(1, b), ab1 = gadd(a, b1);
  if (isnegint(ab1)) S1 = gen_0;
  else
  {
    GEN tmp = Ftaylor(mkvec(a), mkvec(b), z, prec);
    S1 = gmul(divgamma2(b1, ab1, prec), tmp);
  }
  if (isnegint(a)) S2 = gen_0;
  else
  {
    GEN tmp = Ftaylor(mkvec(ab1), mkvec(gaddsg(1, b1)), z, prec);
    S2 = gmul(divgamma2(gneg(b1), a, prec), tmp);
    S2 = gmul(S2, gpow(z, b1, prec));
  }
  S1 = gadd(S1, S2);
  if (gexpo(S1)-gexpo(S2) >= prec2nbits(prec0) - prec2nbits(prec))
    return S1;
  prec = precdbl(prec);
  a = gprec_wensure(a, prec);
  b = gprec_wensure(b, prec);
  z = gprec_wensure(z, prec);
  return hyperu_F11(a, b, z, prec0, prec);
}
/* one branch of this must assume x > 0 (a,b complex); see Temme, The
 * numerical computation of the confluent hypergeometric function U(a,b,z),
 * Numer. Math. 41 (1983), no. 1, 63-82. */
static GEN
hyperu_i(GEN a, GEN b, GEN x, long prec)
{
  GEN u, S, P, T, zf, a1;
  long k, n, bit, l, bigx;

  if (gequal0(imag_i(x))) x = real_i(x);
  l = precision(x); if (!l) l = prec;
  prec = l;
  a1 = gaddsg(1, gsub(a,b));
  P = gmul(a1, a);
  S = gadd(a1, a);
  n = (long)(prec2nbits_mul(l, M_LN2) + M_PI*sqrt(dblmodulus(P)));
  bigx = dbllog2(x) >= log2((double)n);
  if (!bigx && (!isint(b,&b) || typ(x) == t_COMPLEX || gsigne(x) <= 0))
  {
    if (typ(b) == t_INT)
    {
      bit = prec2nbits(l); l += l-2;
      b = gadd(b, real2n(-bit, l));
      a = gprec_wensure(a, l);
      x = gprec_wensure(x, l);
    }
    return hyperu_F11(a, b, x, l, l);
  }
  bit = prec2nbits(l)-1;
  l += EXTRAPREC64;
  T = gadd(gadd(P, gmulsg(n-1, S)), sqru(n-1));
  x = gtofp(x, l);
  if (!bigx)
  { /* this part only works if x is real and positive; only used with b t_INT */
    pari_sp av2 = avma;
    GEN q, v, c, s = real_1(l), t = real_0(l);
    for (k = n-1; k >= 0; k--)
    { /* T = (a+k)*(a1+k) = a*a1 + k(a+a1) + k^2 = previous(T) - S - 2k + 1 */
      GEN p1 = gdiv(T, mulss(-n, k+1));
      s = gprec_wtrunc(gaddgs(gmul(p1,s), 1), l);
      t = gprec_wtrunc(gadd(gmul(p1,t), gaddgs(a,k)), l);
      if (!k) break;
      T = gsubgs(gsub(T, S), 2*k-1);
      if (gc_needed(av2,3)) gerepileall(av2, 3, &s,&t,&T);
    }
    q = utor(n, l);
    zf = gpow(utoi(n), gneg_i(a), l);
    u = gprec_wensure(gmul(zf, s), l);
    v = gprec_wensure(gmul(zf, gdivgs(t,-n)), l);
    for(;;)
    {
      GEN p1, e, f, d = real_1(l), qmb = gsub(q,b);
      pari_sp av3;
      c = divur(5,q); if (expo(c) >= -1) c = real2n(-1, l);
      p1 = subsr(1, divrr(x,q)); if (cmprr(c,p1) > 0) c = p1;
      togglesign(c); av3 = avma;
      e = u;
      f = v;
      for(k = 1;; k++)
      {
        GEN w = gadd(gmul(gaddgs(a,k-1),u), gmul(gaddgs(qmb,1-k),v));
        u = gmul(divru(q,k),v);
        v = gdivgu(w, k);
        d = mulrr(d, c);
        e = gadd(e, gmul(d,u));
        f = gadd(f, p1 = gmul(d,v));
        if (gequal0(p1) || gexpo(p1) - gexpo(f) <= 1-prec2nbits(precision(p1)))
          break;
        if (gc_needed(av3,3)) gerepileall(av3,5,&u,&v,&d,&e,&f);
      }
      u = e;
      v = f;
      q = mulrr(q, addrs(c,1));
      if (expo(x) - expo(subrr(q,x)) >= bit) break;
      gerepileall(av2, 3, &u,&v,&q);
    }
  }
  else
  { /* this part works for large complex x */
    GEN zz = gneg_i(ginv(x)), s = gen_1;
    zf = gpow(x, gneg_i(a), l);
    for (k = n-1; k >= 0; k--)
    {
      s = gaddsg(1, gmul(gmul(T, gdivgu(zz,k+1)), s));
      if (!k) break;
      T = gsubgs(gsub(T, S), 2*k-1);
    }
    u = gmul(s, zf);
  }
  return gprec_wtrunc(u, prec);
}
GEN
hyperu(GEN a, GEN b, GEN x, long prec)
{ pari_sp av = avma; return gerepilecopy(av, hyperu_i(a,b,x,prec)); }

static GEN
mkendpt(GEN z, GEN a)
{
  a = real_i(a);
  if (gcmpgs(a,-1) <= 0) pari_err_IMPL("hypergeom for these parameters");
  return (gcmpgs(a,1) >= 0 || gequal0(a))? z: mkvec2(z, a);
}

/*z != 0 */
static GEN
F20(GEN a, GEN b, GEN z, long prec)
{
  GEN U;
  z = gneg_i(z); U = hyperu_i(a, gadd(a, gsubsg(1, b)), ginv(z), prec);
  return gmul(U, gpow(z, gneg(a), prec));
}

static GEN F21(GEN a, GEN b, GEN c, GEN z, long prec);

static GEN
fF31(void *E, GEN t)
{
  pari_sp av = avma;
  GEN a1 = gel(E,1), b = gel(E,2), c = gel(E,3), d = gel(E,4), z = gel(E,5);
  long prec = precision(t);
  return gerepileupto(av, gmul(gmul(gexp(gneg(t), prec), gpow(t, a1, prec)),
                               F21(b, c, d, gmul(t, z), prec)));
}
/* F31(a,b,c;d,z) = \int_0^oo\exp(-t)t^{a-1}F_{21}(b,c;d,tz) / gamma(a) */
static GEN
F31(GEN a, GEN b, GEN c, GEN d, GEN z, long prec)
{
  GEN p1, p2, a1;
  if (gcmp(real_i(a), real_i(b)) < 0) swap(a,b);
  if (gcmp(real_i(a), real_i(c)) < 0) swap(a,c);
  if (gsigne(real_i(a)) <= 0) pari_err_IMPL("F31 with a, b, and c <= 0");
  a1 = gsubgs(a,1);
  p1 = mkendpt(gen_0, a1);
  p2 = mkvec2(mkoo(), gen_1);
  return gdiv(intnum(mkvecn(5, a1, b, c, d, z), fF31, p1, p2, NULL, prec),
              ggamma(a, prec));
}

/* F32(a,b,c;d,e;z)=\int_0^1,x^{c-1}(1-x)^{e-c-1}F21(a,b,d,x*z)*gamma(e)/(gamma(c)*gamma(e-c)) */
static GEN
fF32(void *E, GEN x)
{
  pari_sp av = avma;
  GEN c1 = gel(E,1), ec = gel(E,2), a = gel(E,3), b = gel(E,4);
  GEN d = gel(E,5), z = gel(E,6), T;
  long prec = precision(x);
  T = F21(a, b, d, gmul(x, z), prec);
  if (!gequal0(c1)) T = gmul(T, gpow(x,c1,prec));
  if (!gequal0(ec)) T = gmul(T, gpow(gsubsg(1,x),ec,prec));
  return gerepileupto(av,T);
}

static GEN
myint32(GEN a, GEN b, GEN c, GEN d, GEN e, GEN z, long prec)
{
  GEN c1 = gsubgs(c,1), c2 = gsub(e, gaddgs(c,1));
  GEN p0 = mkendpt(gen_0, c1);
  GEN p1 = mkendpt(gen_1, c2);
  return intnum(mkvecn(6, c1, c2, a,b,d, z), fF32, p0, p1, NULL, prec);
}

static void
check_hyp1(GEN x)
{
  if (gsigne(real_i(x)) <= 0)
    pari_err_DOMAIN("hypergeom","real(vecsum(D)-vecsum(N))", "<=", gen_0, x);
}

/* 0 < re(z) < e ? */
static int
ok_F32(GEN z, GEN e)
{ GEN x = real_i(z); return gsigne(x) > 0 && gcmp(e, x) > 0; }

/* |z| very close to 1 but z != 1 */
static GEN
F32(GEN N, GEN D, GEN z, long prec)
{
  GEN tmp, a,b,c, d,e, re;
  a = gel(N,1); d = gel(D,1);
  b = gel(N,2); e = gel(D,2);
  c = gel(N,3);
  if (gcmp(real_i(e), real_i(d)) < 0) swap(e,d);
  re = real_i(e);
  if (!ok_F32(c, re))
  {
    if (ok_F32(b, re))  {swap(b,c);}
    else if (ok_F32(a, re)) {swap(a,c);}
    else pari_err_IMPL("3F2 for these arguments");
  }
  tmp = gdiv(ggamma(e, prec), mulgamma2(c, gsub(e,c), prec));
  return gmul(tmp, myint32(a, b, c, d, e, z, prec));
}

static GEN
poch(GEN a, long n, long prec)
{
  GEN S = real_1(prec);
  long j;
  for (j = 0; j < n; j++) S = gmul(S, gaddsg(j, a));
  return S;
}
static GEN
vpoch(GEN a, long n)
{
  GEN v = cgetg(n+1, t_VEC);
  long j;
  gel(v,1) = a;
  for (j = 1; j < n; j++) gel(v,j+1) = gmul(gel(v,j), gaddsg(j, a));
  return v;
}

static GEN
Npoch(GEN a, long n) { return gnorm(poch(a, n, LOWDEFAULTPREC)); }
static GEN
Npochden(GEN a, long n)
{
  GEN L = Npoch(a, n), r = ground(real_i(a));
  if (signe(r) <= 0)
  {
    GEN t = gnorm(gsub(a, r));
    if (gcmpgs(t,1) > 0) L = gmul(L, t);
  }
  return L;
}

/* |x + z|^2 */
static GEN
normpol2(GEN z)
{
  GEN a = real_i(z), b = imag_i(z);
  return deg2pol_shallow(gen_1, gmul2n(a,1), gadd(gsqr(a),gsqr(b)), 0);
}
/* \prod |x + v[i]|^2 */
static GEN
vnormpol2(GEN v)
{
  long i, l = lg(v);
  GEN P;
  if (l == 1) return pol_1(0);
  P = normpol2(gel(v,1));
  for (i = 2; i < l; i++) P = RgX_mul(P, normpol2(gel(v,i)));
  return P;
}

static long
precFtaylor(GEN N, GEN D, GEN z, long *pmi)
{
  GEN v, ma, P = vnormpol2(D), Q = vnormpol2(N), Nz = gnorm(z);
  double wma, logNz = (gexpo(Nz) < -27)? -27: dbllog2(Nz) / 2;
  long pr = LOWDEFAULTPREC, prec = precision(z);
  long lN = lg(N), lD = lg(D), mi, j, i, lv;

  P = RgX_shift_shallow(P, 2);
  /* avoid almost cancellation of leading coeff if |z| ~ 1 */
  if (!prec || fabs(logNz) > 1e-38) Q = RgX_Rg_mul(Q,Nz);
  for (j = 1, ma = NULL; j < lN; ++j)
  {
    GEN Nj = gel(N,j);
    if (isint(Nj,&Nj) && signe(Nj) <= 0
        && (!ma || abscmpii(ma, Nj) < 0)) ma = Nj;
  }
  /* use less sharp fujiwara_bound_real(,1) ? */
  v = ground(realroots(gsub(P,Q), mkvec2(gen_0,mkoo()), pr));
  v = ZV_to_zv(v); lv = lg(v);
  if (ma)
  {
    long sma = is_bigint(ma)? LONG_MAX: maxss(labs(itos(ma))-1, 1);
    for (i = 1; i < lv; ++i) v[i] = maxss(minss(sma, v[i]), 1);
  }
  for (i = 1, wma = 0., mi = 0; i < lv; ++i)
  {
    GEN t1 = gen_1, t2 = gen_1;
    long u = v[i];
    double t;
    mi = maxss(mi, u);
    for (j = 1; j < lN; j++) t1 = gmul(t1, Npoch(gel(N,j), u));
    for (j = 1; j < lD; j++) t2 = gmul(t2, Npochden(gel(D,j), u));
    t = dbllog2(gdiv(t1,t2)) / 2 + u * logNz - dbllog2(mpfactr(u,pr));
    wma = maxdd(wma, t); /* t ~ log2 | N(u)/D(u) z^u/u! | */
  }
  /* make up for exponential decrease in exp() */
  if (gsigne(real_i(z)) < 0) wma -= gtodouble(real_i(z)) / M_LN2;
  *pmi = mi; return ceil(wma/BITS_IN_LONG) + 1;
}

static GEN
Ftaylor(GEN N, GEN D, GEN z, long prec)
{
  pari_sp av;
  GEN C, S;
  long i, j, ct, lN = lg(N), lD = lg(D), pradd, mi, bitmin, tol;
  pradd = precFtaylor(N, D, z, &mi);
  if (pradd > 0)
  {
    prec += pradd;
    N = gprec_wensure(N, prec);
    D = gprec_wensure(D, prec);
    z = gprec_wensure(z, prec);
  }
  bitmin = -(prec2nbits(prec) + 10);
  S = C = real_1(prec); ct = 0; j = 0; tol = 0;
  av = avma;
  for(;;)
  {
    GEN a = gen_1, b = gen_1;
    for (i = 1; i < lN; i++) a = gmul(a, gaddsg(j, gel(N,i)));
    for (i = 1; i < lD; i++) b = gmul(b, gaddsg(j, gel(D,i)));
    C = gmul(C, gmul(gdiv(a, gmulsg(j+1, b)), z));
    if (gequal0(C)) break;
    if (j > mi) tol = gequal0(S)? 0: gexpo(C) - gexpo(S);
    S = gadd(S, C); ++j;
    if (j > mi)
    { if (tol > bitmin) ct = 0; else if (++ct >= lN+lD-2) break; }
    if (gc_needed(av, 1)) gerepileall(av, 2, &S, &C);
  }
  return S;
}

static GEN
bind(GEN a, GEN b, GEN c, long ind)
{
  switch(ind)
  {
    case 1: case 2: return gsub(c, b);
    case 5: case 6: return gsub(gaddsg(1, a), c);
    default: return b;
  }
}
static GEN
cind(GEN a, GEN b, GEN c, long ind)
{
  switch(ind)
  {
    case 1: case 6: return gaddsg(1, gsub(a,b));
    case 4: case 5: return gsub(gaddsg(1, gadd(a,b)), c);
    default: return c;
  }
}
static GEN
zind(GEN z, long ind)
{
  switch(ind)
  {
    case 1: return ginv(gsubsg(1, z));
    case 2: return gdiv(z, gsubgs(z, 1));
    case 4: return gsubsg(1, z);
    case 5: return gsubsg(1, ginv(z));
    case 6: return ginv(z);
  }
  return z;
}

/* z not 0 or 1, c not a nonpositive integer */
static long
F21ind(GEN a, GEN b, GEN c, GEN z, long bit)
{
  GEN v = const_vec(6, mkoo());
  long ind = 0, B = bit - 5;
  const long LD = LOWDEFAULTPREC;
  if (!isnegint_approx(cind(a,b,c, 1),B)) gel(v,1) = gabs(zind(z,1), LD);
  gel(v,2) = gabs(zind(z,2), LD);
  gel(v,3) = gabs(z, LD);
  if (!isnegint_approx(cind(a,b,c, 4),B)) gel(v,4) = gabs(zind(z,4), LD);
  if (!isnegint_approx(cind(a,b,c, 5),B)) gel(v,5) = gabs(zind(z,5), LD);
  if (!isnegint_approx(cind(a,b,c, 6),B)) gel(v,6) = gabs(zind(z,6), LD);
  ind = vecindexmin(v); /* |znew| <= 1; close to 1 ? */
  return (gexpo(gsubgs(gel(v,ind),1)) > -maxss(bit / 4, 32))? -ind: ind;
}
static GEN
mul4(GEN a, GEN b, GEN c, GEN d) { return gmul(a,gmul(b, gmul(c, d))); }
static GEN
mul3(GEN a, GEN b, GEN c) { return gmul(a,gmul(b, c)); }

/* (1 - zt)^a t^b (1-t)^c */
static GEN
fF212(void *E, GEN t)
{
  GEN z = gel(E,1), a = gel(E,2), b = gel(E,3), c = gel(E,4);
  GEN u = gsubsg(1, gmul(z, t));
  long prec = precision(t);
  return mul3(gpow(u, a, prec), gpow(t, b, prec), gpow(gsubsg(1,t), c, prec));
}

/* (1 - zt)^a T(1-zt) t^b (1-t)^c */
static GEN
fF21neg2(void *E, GEN t)
{
  GEN z = gel(E,1), a = gel(E,2), b = gel(E,3), c = gel(E,4), T = gel(E,5);
  GEN u = gsubsg(1, gmul(z, t));
  long prec = precision(t);
  return mul4(gsubst(T, 0, u), gpow(u, a, prec), gpow(t, b, prec),
              gpow(gsubsg(1,t), c, prec));
}

/* N >= 1 */
static GEN
F21lam(long N, GEN a, GEN c)
{
  long i;
  GEN C = vecbinomial(N), S = cgetg(N+2, t_VEC);
  GEN vb = vpoch(gsub(c,a), N), va = vpoch(a, N);
  gel(S,1) = gel(va,N);
  for (i = 1; i < N; i++) gel(S,i+1) = mul3(gel(C,i+1), gel(vb,i), gel(va,N-i));
  gel(S,i+1) = gel(vb,N); return RgV_to_RgX(S,0);
}

/* F(-m,b; c; z), m >= 0 */
static GEN
F21finitetaylor(long m, GEN b, GEN c, GEN z, long prec)
{
  pari_sp av;
  GEN C, S;
  long j, ct, pradd, mi, bitmin, mb;
  if (isnegint2(b, &mb) && mb < m) { b = stoi(-m); m = mb; }
  pradd = precFtaylor(mkvec2(stoi(-m), b), mkvec(c), z, &mi);
  if (pradd > 0)
  {
    prec += pradd;
    b = gprec_wensure(b, prec);
    c = gprec_wensure(c, prec);
    z = gprec_wensure(z, prec);
  }
  bitmin = -(prec2nbits(prec) + 10);
  C = real_1(prec); S = C; ct = 0;
  av = avma;
  for(j = 0; j < m; ++j)
  {
    C = gmul(C, gdiv(gmulsg(j-m, gaddsg(j, b)), gmulsg(j+1, gaddsg(j, c))));
    C = gmul(C, z);
    if (j > mi && !gequal0(S))
    { if (gexpo(C) - gexpo(S) > bitmin) ct = 0; else if (++ct == 3) break; }
    S = gadd(S, C);
    if (gc_needed(av, 1)) gerepileall(av, 2, &S, &C);
  }
  return S;
}

/* c not a nonpositive integer */
static GEN
F21finite_i(long m, GEN b, GEN c, GEN z, GEN B, GEN C, GEN coe, long prec)
{
  return mul3(poch(B, m, prec), gdiv(gpowgs(coe, m), poch(C, m, prec)),
              F21finitetaylor(m, b, c, z, prec));
}

/* F(-m,b; c; z), m >= 0; c not a nonpositive integer */
static GEN
F21finite(long m, GEN b, GEN c, GEN z, long prec)
{
  GEN a = stoi(-m), b1 = b, c1 = c, z1;
  long ind = F21ind(a, b, c, z, prec2nbits(prec)), inda = labs(ind);
  z1 = zind(z, inda);
  if (ind < 0)
  {
    b1 = bind(a, b, c, inda);
    c1 = cind(a, b, c, inda); /* not a nonpositive integer */
  }
  switch (inda)
  {
    case 1: return F21finite_i(m, b1, c1, z1, b, c, gsubsg(1,z), prec);
    case 2: return gmul(gpowgs(gsubsg(1,z), m),
                        F21finitetaylor(m, b1, c, z1, prec));
    case 4: return F21finite_i(m, b1, c1, z1, gsub(c,b), c, gen_1, prec);
    case 5: return F21finite_i(m, b1, c1, z1, gsub(c,b), c, z, prec);
    case 6: return F21finite_i(m, b1, c1, z1, b, c, gneg(z), prec);
    default:return F21finitetaylor(m, b1, c1, z, prec);
  }
}

/**********************************************************/

static GEN
multgam(GEN a, GEN b, GEN c, GEN d, long prec)
{
  if (isnegint(c) || isnegint(d)) return gen_0;
  return gdiv(mulgamma2(a, b, prec), mulgamma2(c, d, prec));
}

static GEN
intnumsplit(void *E, GEN (*f)(void*, GEN), GEN a, GEN b, GEN z, long prec)
{
  if (!z) return intnum(E, f, a, b, NULL, prec);
  return gadd(intnum(E, f, a, z, NULL, prec),
              intnum(E, f, z, b, NULL, prec));
}
/* z != 1 */
static GEN
myint21(void *E, GEN (*f)(void*, GEN), long prec)
{
  GEN z = gel(E,1), a = real_i(gel(E,2)), b = gel(E,3), c = gel(E,4);
  GEN pz = NULL, p0 = mkendpt(gen_0, b), p1 = mkendpt(gen_1, c);
  if (gcmpgs(a, 1) <= 0 && is0(imag_i(z), 10))
  {
    GEN r;
    pz = ginv(z); r = real_i(pz);
    if (gsigne(r) <= 0 || gcmp(r, gen_1) >= 0) pz = NULL;
  }
  if (pz) pz = mkendpt(pz,a);
  else if (gcmpgs(a,-1) <= 0) prec += ((gexpo(a)+1)>>1) * EXTRAPREC64;
  return intnumsplit(E, f, p0, p1, pz, prec);
}

/* Algorithm used for F21(a,b;c;z)
Basic transforms:
  1: (c-b,1+a-b,1/(1-z))
  2: (c-b,c,z/(z-1))
  3: (b,c,z)
  4: (b,b-c+a+1,1-z)
  5: (1+a-c,b-c+a+1,1-1/z)
  6: (1+a-c,1+a-b,1/z)

F21: calls F21_i and increase prec if too much cancellation
F21_i: c is not a non-positive integer
- z ~ 0 or 1: return special value
- if a, b, c-b or c-a a non-positive integer: use F21finite
- compute index, value of z
   if |z| < 1-epsilon return F21taylorind
   if Re(b)<=0, swap and/or recurse
   so may assume Re(b)>0 and Re(a)>=Re(b) and integrate.

F21finite:
- compute index, value of z
- call F21finitetaylor

F21ind: find best index (1 to 6, -1 to -6 if |z| < 1-epsilon)
F21finitetaylor: a or b in Z_{<=0}; calls precFtaylor

F21taylorind: in case 2, may lose accuracy, possible bug.
- calls F21taylor[1456] or F21taylor

F21taylor: calls Ftaylor / gamma: may lose accuracy

FBaux1: F21taylor twice
FBaux2: F21finitelim + F21taylorlim
F21taylor[45]: if c-(a+b) integer, calls FBaux1, else calls FBaux2
F21taylor[16]: if b-a integer, calls FBaux1, else calls FBaux2

F21taylorlim: calls precFtaylor then compute
F21finitelim: direct */
static GEN F21taylorind(GEN a, GEN b, GEN c, GEN z, long ind, long prec);
/* c not a nonpositive integer */
static GEN
F21_i(GEN a, GEN b, GEN c, GEN z, long prec)
{
  GEN res;
  long m, ind, prec2, bitprec = prec2nbits(prec);
  if (is0(imag_i(z), bitprec)) z = real_i(z);
  if (is0(z, bitprec)) return real_1(prec);
  if (gequal1(z))
  {
    GEN x = gsub(c, gadd(a, b)); check_hyp1(x);
    return multgam(c, x, gsub(c,a), gsub(c,b), prec);
  }
  if (isnegint2(b, &m)) return F21finite(m, a, c, z, prec);
  if (isnegint2(a, &m)) return F21finite(m, b, c, z, prec);
  if (isnegint(gsub(c, b))) swap(a, b);
  if (isnegint2(gsub(c, a), &m))
  {
    GEN x = gpow(gsubsg(1, z), gneg(gaddsg(m, b)), prec);
    return gmul(x, F21finite(m, gsub(c, b), c, z, prec));
  }
  /* Here a, b, c, c-a, c-b are not nonpositive integers */
  ind = F21ind(a, b, c, z, bitprec);
  prec2 = prec + EXTRAPREC64;
  a = gprec_wensure(a,prec2);
  b = gprec_wensure(b,prec2);
  c = gprec_wensure(c,prec2);
  z = gprec_wensure(z,prec2);
  if (ind < 0) return gmul(ggamma(c, prec), F21taylorind(a,b,c, z, ind, prec));
  if (gsigne(real_i(b)) <= 0)
  {
    if (gsigne(real_i(a)) <= 0)
    {
      GEN p1,p2;
      if (gcmp(real_i(b), real_i(a)) < 0) swap(a,b);
      /* FIXME: solve recursion as below with F21auxpol */
      p1 = gmul(gsubsg(1, z), F21_i(a, gaddsg(1,b), c, z, prec));
      p2 = gmul(gmul(gsubsg(1, gdiv(a,c)), z),
                F21_i(a, gaddsg(1,b), gaddsg(1,c), z, prec));
      return gadd(p1, p2);
    }
    swap(a,b);
  }
  if (gcmp(real_i(a), real_i(b)) < 0 && gsigne(real_i(a)) > 0) swap(a,b);
  /* Here real(b) > 0 and either real(a) <= 0 or real(a) > real(b) */
  if (gcmp(real_i(c), real_i(b)) <= 0)
  {
    long N = 1 + itos(gfloor(gsub(real_i(b),real_i(c)))); /* >= 1 */
    GEN T = F21lam(N, a, c), c0 = c;
    void *E;
    c = gaddsg(N,c);
    E = (void*)mkvec5(z, gsubsg(-N,a), gsubgs(b,1), gsubgs(gsub(c,b),1), T);
    res = gdiv(myint21(E, fF21neg2, prec2), poch(c0, N, prec));
  }
  else
  {
    void *E = (void*)mkvec4(z, gneg(a), gsubgs(b,1), gsubgs(gsub(c,b),1));
    res = myint21(E, fF212, prec2);
  }
  return gmul(multgam(gen_1, c, b, gsub(c,b), prec), res);
}

/* c not a non-positive integer */
static GEN
F21(GEN a, GEN b, GEN c, GEN z, long prec)
{
  GEN res = F21_i(a, b, c, z, prec);
  long ex = labs(gexpo(res)), bitprec = prec2nbits(prec);
  if (ex > bitprec)
  {
    prec = nbits2prec(ex + bitprec);
    res = F21_i(gprec_wensure(a,prec), gprec_wensure(b,prec),
                gprec_wensure(c,prec), gprec_wensure(z,prec), prec);
  }
  return res;
}

static GEN
F21taylor(GEN a, GEN b, GEN c, GEN z, long prec)
{
  pari_sp av = avma;
  GEN r = gdiv(Ftaylor(mkvec2(a,b), mkvec(c), z, prec), ggamma(c, prec));
  return gerepileupto(av, r);
}

static GEN
F21taylorlim(GEN N, long m, GEN z, GEN Z, long ind, long prec)
{
  pari_sp av;
  GEN C, P, S, a, b;
  long j, ct, pradd, mi, fl, bitmin, tol, si = (ind == 5 || ind == 6)? -1: 1;
  pradd = precFtaylor(N, mkvec(stoi(m + 1)), z, &mi);
  if (pradd)
  {
    prec += pradd;
    N = gprec_wensure(N, prec);
    z = gprec_wensure(z, prec);
    Z = gprec_wensure(Z, prec);
  }
  av = avma; a = gel(N,1); b = gel(N,2);
  bitmin = -(prec2nbits(prec) + 10);
  P = glog(Z, prec); if (ind == 4 || ind == 5) P = gneg(P);
  P = gadd(P, gsub(gpsi(stoi(m+1), prec), mpeuler(prec)));
  P = gsub(P, gadd(gpsi(a, prec), gpsi(si == -1? gsubsg(1, b): b, prec)));
  C = real_1(prec); ct = 0; tol = 0; S = P;
  for(j = 0, fl = 1;;)
  {
    GEN v1 = gaddgs(a, j), v2 = gaddgs(b, j);
    long J = (j+1) * (j+1+m);
    C = gmul(C, gdivgs(gmul(z, v1), J));
    if (gequal0(v2)) fl = 0; else C = gmul(C, v2);
    if (j > mi) tol = gequal0(S) ? 0 : gexpo(C) - gexpo(S);
    if (fl)
    {
      P = gadd(P, gsub(uutoQ(2*j+2+m, J), gadd(ginv(v1), ginv(v2))));
      S = gadd(S, gmul(C, P));
    }
    else
      S = (si == 1)? gadd(S, C): gsub(S, C);
    if (++j > mi) { if (tol > bitmin) ct = 0; else if (++ct == 3) break; }
    if (gc_needed(av, 1)) gerepileall(av, 3, &S, &C, &P);
  }
  return gdiv(S, mpfact(m));
}

/* N = [a,b]; (m-1)! sum_{1 <= k < m} (-z)^k (a)_k (b)_k / k! (m-k)!*/
static GEN
F21finitelim(GEN N, long m, GEN z, long prec)
{
  GEN C, S, a, b;
  long j;
  if (!m) return gen_0;
  a = gel(N,1);
  b = gel(N,2);
  S = C = real_1(prec + EXTRAPREC64);
  for (j = 1; j < m; j++)
  {
    GEN v1 = gaddsg(j-1, a), v2 = gaddsg(j-1, b);
    C = gdivgs(gmul(C, gmul(gmul(v1, v2), z)), j*(j-m));
    S = gadd(S, C);
  }
  return gmul(S, mpfact(m-1));
}

static GEN
OK_gadd(GEN x, GEN y, long prec0, long *pprec,
        GEN *d1,GEN *d2,GEN *d3,GEN *d4,GEN *d5,GEN *d6,GEN *d7,GEN *d8)
{
  long prec = *pprec;
  GEN z = gadd(x,y);
  if (!gequal0(z) && gexpo(z)-gexpo(x) >= prec2nbits(prec0) - prec2nbits(prec))
    return z;
  *pprec = prec = precdbl(prec);
  *d1 = gprec_wensure(*d1, prec); *d2 = gprec_wensure(*d2, prec);
  *d3 = gprec_wensure(*d3, prec); *d4 = gprec_wensure(*d4, prec);
  *d5 = gprec_wensure(*d5, prec); *d6 = gprec_wensure(*d6, prec);
  *d7 = gprec_wensure(*d7, prec); *d8 = gprec_wensure(*d8, prec);
  return NULL;
}

static GEN
FBaux1(GEN v1, GEN g1, GEN c1, GEN v2, GEN g2, GEN c2, GEN z, GEN bma,
       long prec0, long prec)
{
  GEN pi = mppi(prec);
  for (;;)
  {
    GEN t1 = gdiv(c1, mulgammav2(g1, prec)), r1;
    GEN t2 = gdiv(c2, mulgammav2(g2, prec)), r2, F;
    r1 = gmul(t1, F21taylor(gel(v1,1), gel(v1,2), gel(v1,3), z, prec));
    r2 = gmul(t2, F21taylor(gel(v2,1), gel(v2,2), gel(v2,3), z, prec));
    F = OK_gadd(r1, r2, prec0, &prec, &c1,&c2, &g1,&g2, &v1,&v2, &z,&bma);
    if (F) return gmul(F, gdiv(pi, gsin(gmul(pi, bma), prec)));
  }
}

static GEN
FBaux2(GEN v1, GEN g1, GEN c1, long m, GEN z1, GEN c2, GEN g2, GEN v2, GEN z2,
       GEN Z2, long ind, long prec)
{
  GEN t1 = gdiv(c1, mulgammav2(g1, prec)), r1;
  GEN t2 = gdiv(c2, mulgammav2(g2, prec)), r2;
  r1 = gmul(t1, F21finitelim(v1, m, z1, prec));
  r2 = gmul(t2, F21taylorlim(v2, m, z2, Z2, ind, prec));
  return gadd(r1, r2);
}

/* 1 / (1-z) */
static GEN
F21taylor1(GEN a, GEN b, GEN c, GEN z, long prec)
{
  GEN bma = gsub(b, a), coe1, coe2, z1, g1, g2, v1, v2, Z;
  long m;
  if (!islong(bma,&m, prec))
  {
    GEN b1, c1, e1, c2;
    b1 = gsub(c, b);
    c1 = gsubsg(1, bma);
    e1 = gsub(c, a);
    coe1 = gpow(gsubsg(1, z), gneg(a), prec);
    c2 = gaddgs(bma, 1);
    coe2 = gneg(gpow(gsubsg(1, z), gneg(b), prec));
    z1 = ginv(gsubsg(1, z));
    return FBaux1(mkvec3(a,b1,c1), mkvec2(b,e1), coe1, mkvec3(b,e1,c2),
                  mkvec2(a,b1), coe2, z1, bma, prec, prec);
  }
  if (m < 0) { swap(a,b); m = -m; }
  Z = gsubsg(1, z);
  coe1 = gpow(Z, gneg(a), prec);
  v2 = g1 = mkvec2(gaddgs(a, m), gsub(c, a));
  v1 = mkvec2(a, gsub(c, gaddgs(a, m)));
  z1 = ginv(Z);
  coe2 = gmul(coe1, gpowgs(z1, m)); if (m & 1L) coe2 = gneg(coe2);
  g2 = mkvec2(a, gsub(c, gaddgs(a, m))); /* 15.8.9 */
  return FBaux2(v1, g1, coe1, m, z1, coe2, g2, v2, z1, Z, 1, prec);
}

/* 1 - z */
static GEN
F21taylor4(GEN a, GEN b, GEN c, GEN z, long prec)
{
  GEN bma = gsub(c, gadd(a, b)), coe2, z1, z2, g1, g2, v1, v2;
  long m;
  if (!islong(bma,&m, prec))
  {
    GEN c1, a2, b2, c2;
    c1 = gsubsg(1, bma);
    a2 = gsub(c, a);
    b2 = gsub(c, b);
    c2 = gaddsg(1, bma);
    z1 = gsubsg(1, z); coe2 = gneg(gpow(z1, bma, prec));
    return FBaux1(mkvec3(a,b,c1), mkvec2(a2,b2), gen_1, mkvec3(a2,b2,c2),
                  mkvec2(a,b), coe2, z1, bma, prec, prec);
  }
  if (m < 0)
  {
    GEN F = F21taylor4(gaddgs(a,m), gaddgs(b,m), gaddgs(gadd(a,b), m), z, prec);
    return gmul(gpowgs(gsubsg(1,z), m), F);
  }
  v2 = g1 = mkvec2(gaddgs(a,m), gaddgs(b,m));
  v1 = g2 = mkvec2(a, b);
  z1 = gsubgs(z, 1);
  z2 = gneg(z1); coe2 = gpowgs(z1, m); /* 15.8.10 */
  return FBaux2(v1, g1, gen_1, m, z1, coe2, g2, v2, z2, z2, 4, prec);
}

/* 1 - 1/z */
static GEN
F21taylor5(GEN a, GEN b, GEN c, GEN z, long prec)
{
  GEN bma = gsub(c, gadd(a,b)), tmp, coe1, coe2, z1, g1, g2, v1, v2, z2;
  long m;
  if (!islong(bma,&m, prec))
  {
    GEN b1, c1, d1, e1, b2, c2;
    d1 = gsub(c, a);
    b1 = gsubsg(1, d1);
    c1 = gsubsg(1, bma);
    e1 = gsub(c, b);
    b2 = gsubsg(1, a);
    c2 = gaddsg(1, bma);
    coe1 = gpow(z, gneg(a), prec);
    coe2 = gneg(gmul(gpow(gsubsg(1, z), bma, prec), gpow(z, gneg(d1), prec)));
    z1 = gsubsg(1, ginv(z));
    return FBaux1(mkvec3(a,b1,c1), mkvec2(d1,e1), coe1,
                  mkvec3(d1,b2,c2), mkvec2(a, b), coe2, z1, bma, prec, prec);
  }
  /* c - (a + b) ~ m */
  if (m < 0)
  {
    tmp = F21taylor5(gaddgs(a,m), gaddgs(b,m), c, z, prec);
    return gmul(gpowgs(gsubsg(1, z), m), tmp);
  }
  g1 = mkvec2(gaddgs(a,m), gaddgs(b,m));
  v1 = mkvec2(a, gsubsg(1-m, b));
  v2 = mkvec2(gaddgs(a,m), gsubsg(1,b));
  z1 = gsubgs(ginv(z), 1);
  z2 = gneg(z1);
  g2 = mkvec2(a, b);
  coe1 = gpow(z, gneg(a), prec);
  coe2 = gmul(coe1, gpowgs(z2, m)); /* 15.8.11 */
  return FBaux2(v1, g1, coe1, m, z1, coe2, g2, v2, z2, z1, 5, prec);
}

/* 1 / z */
static GEN
F21taylor6(GEN a, GEN b, GEN c, GEN z, long prec)
{
  GEN bma = gsub(b, a), cma, am, coe1, coe2, z1, g1, g2, v1, v2, z2, Z;
  long m;
  if (!islong(bma,&m, prec))
  {
    GEN e1, e2, b1, b2, c1, c2;
    b1 = gaddgs(gsub(a,c), 1);
    c1 = gsubsg(1, bma);
    e1 = gsub(c,a);
    b2 = gaddgs(gsub(b,c), 1);
    c2 = gaddgs(bma, 1);
    e2 = gsub(c, b);
    coe1 = gpow(gneg(z), gneg(a), prec);
    coe2 = gneg(gpow(gneg(z), gneg(b), prec));
    z1 = ginv(z);
    return FBaux1(mkvec3(a,b1,c1), mkvec2(b,e1), coe1,
                  mkvec3(b,b2,c2), mkvec2(a,e2), coe2, z1, bma, prec, prec);
  }
  /* b - a ~ m */
  if (m < 0) { swap(a,b); m = -m; }
  cma = gsub(c, a); am = gaddgs(a,m);
  Z = gneg(z);
  coe1 = gpow(Z, gneg(a), prec);
  coe2 = gdiv(coe1, gpowgs(z, m));
  g1 = mkvec2(am, cma);
  v1 = mkvec2(a, gsubsg(1, cma));
  g2 = mkvec2(a, gsubgs(cma, m));
  v2 = mkvec2(am, gsubsg(m+1, cma));
  z2 = ginv(z);
  z1 = gneg(z2); /* 15.8.8 */
  return FBaux2(v1, g1, coe1, m, z1, coe2, g2, v2, z2, Z, 6, prec);
}

/* (new b, new c, new z): given by bind, cind, zind
 * case 1: (c-b,1+a-b,1/(1-z))
 * case 2: (c-b,c,z/(z-1))
 * case 3: (b,c,z)
 * case 4: (b,b-c+a+1,1-z)
 * case 5: (1+a-c,b-c+a+1,1-1/z)
 * case 6: (1+a-c,1+a-b,1/z) */
static GEN
F21taylorind(GEN a, GEN b, GEN c, GEN z, long ind, long prec)
{
  GEN res;
  switch (labs(ind))
  {
    case 1: res = F21taylor1(a, b, c, z, prec); break;
    case 2: res = F21taylor(a, gsub(c, b), c, gdiv(z, gsubgs(z, 1)), prec);
            res = gmul(res, gpow(gsubsg(1,z), gneg(a), prec)); break;
    case 3: res = F21taylor(a, b, c, z, prec); break;
    case 4: res = F21taylor4(a, b, c, z, prec); break;
    case 5: res = F21taylor5(a, b, c, z, prec); break;
    default:res = F21taylor6(a, b, c, z, prec); break;
  }
  return gprec_wtrunc(res, prec);
}

static long
hypersimplify(GEN *pn, GEN *pd)
{
  GEN n = *pn, d = *pd;
  long j, ld = lg(d), ln = lg(n);
  for (j = 1; j < ld; j++)
  {
    GEN t = gel(d, j);
    long k;
    for (k = 1; k < ln; k++)
      if (gequal(t, gel(n, k)))
      {
        *pd = vecsplice(d, j);
        *pn = vecsplice(n, k); return hypersimplify(pd, pn) + 1;
      }
  }
  return 0;
}

static GEN
f_pochall(void *E, GEN n)
{
  GEN S, a, tmp, N = gel(E,1), D = gel(E,2);
  long j, prec = itou(gel(E,3));
  S = gen_0;
  for (j = 1; j < lg(N); j++)
  {
    a = gel(N, j);
    tmp = gsub(glngamma(gadd(n, a), prec), glngamma(a, prec));
    S = gadd(S, tmp);
  }
  for (j = 1; j < lg(D); j++)
  {
    a = gel(D, j);
    tmp = gsub(glngamma(gadd(n, a), prec), glngamma(a, prec));
    S = gsub(S, tmp);
  }
  return gexp(gsub(S, glngamma(gaddsg(1, n), prec)), prec);
}
static GEN
f_pochall_alt(void *E, GEN n)
{ GEN z = f_pochall(E,n); return mpodd(n)? gneg(z): z; }
/* z = \pm1 */
static GEN
sumz(GEN N, GEN D, long z, long prec)
{
  void *E = (void*)mkvec3(N, D, utoi(prec));
  GEN tab, be;
  if (z == -1) return sumalt(E, f_pochall_alt, gen_0, prec);
  be = gsub(vecsum(D), vecsum(N)); check_hyp1(be);
  tab = sumnummonieninit(be, NULL, gen_0, prec);
  return sumnummonien(E, f_pochall, gen_0, tab, prec);
}

static GEN
hypergeom_arg(GEN x)
{
  if (!x) return cgetg(1,t_VEC);
  return (typ(x) == t_VEC)? x: mkvec(x);
}

/* [x[i]+k, i=1..#v] */
static GEN
RgV_z_add(GEN x, long k)
{
  if (!k) return x;
  pari_APPLY_same(gaddgs(gel(x,i), k));
}
static GEN
vp(GEN a, GEN p, GEN dft)
{
  long v = gvaluation(a, p);
  return v < 0? stoi(v): dft;
}
static GEN
Qp_hypergeom(GEN N, GEN D, GEN z)
{
  pari_sp av = avma;
  GEN r, S = gen_1, R = gen_1, p = gel(z, 2), dft = ginv(subis(p, 1));
  long l, i, prec = precp(z) + valp(z) + 1;

  r = gsub(stoi(valp(z)), dft);
  l = lg(N); for (i = 1; i < l; i++) r = gadd(r, vp(gel(N,i), p, dft));
  l = lg(D); for (i = 1; i < l; i++) r = gsub(r, vp(gel(D,i), p, dft));
  if (gsigne(r) <= 0) pari_err(e_MISC, "divergent p-adic hypergeometric sum");
  l = itou(gceil(gdivsg(prec, r)));
  for (i = 1; i <= l; i++)
  {
    GEN u = vecprod(RgV_z_add(N, i-1)), v = vecprod(RgV_z_add(D, i-1));
    R = gmul(R, gmul(z, gdiv(u, gmulsg(i, v))));
    S = gadd(S, R);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"hypergeom, i = %ld / %ld", i,l);
      gerepileall(av, 2, &R, &S);
    }
  }
  return S;
}

/* assume is_scalar_t(typ(z)) */
static GEN
hypergeom_i(GEN N, GEN D, GEN z, long prec)
{
  long nN, nD;
  if (typ(z) == t_PADIC) return Qp_hypergeom(N, D, z);
  if (gequal0(z)) return gen_1;
  nN = lg(N) - 1;
  nD = lg(D) - 1;
  if (nD >= (nN? nN: 2)) return Ftaylor(N, D, z, prec);
  if (nD == nN - 1 && nN >= 3)
  {
    GEN d = gsubsg(1, gabs(z,LOWDEFAULTPREC));
    long ed = gexpo(d);
    /* z in unit disc but "away" from unit circle */
    if (gsigne(d) > 0 && ed > -prec2nbits(prec)/4
        && (nN != 3 || ed > -15)) /* For 3F2 we can use integral */
      return Ftaylor(N, D, z, prec);
    if (gequal1(z))  return sumz(N, D, 1, prec);
    if (gequalm1(z)) return sumz(N, D,-1, prec);
  }
  switch (nN)
  {
    case 0:
      if (nD == 0) return gexp(z, prec);
      if (nD == 1) return F01(gel(D,1), z, prec);
    case 1: return F10(gel(N, 1), z, prec);
    case 2:
      if (nD == 0) return F20(gel(N,1), gel(N,2), z, prec);
      if (nD == 1) return F21(gel(N,1), gel(N,2), gel(D,1), z, prec);
    case 3:
      if (nD == 0) break;
      if (nD == 1) return F31(gel(N,1), gel(N,2), gel(N,3), gel(D,1), z, prec);
      if (nD == 2) return F32(N, D, z, prec);
  }
  pari_err_IMPL("this hypergeometric function");
  return NULL; /*LCOV_EXCL_LINE*/
}

static GEN
serhypergeom(GEN N, GEN D, GEN y, long prec)
{
  GEN Di, Ni, S = gen_1, R = gen_1, y0 = NULL;
  long v, l, i;
  pari_sp av;
  if (!signe(y)) return gadd(gen_1, y);
  v = valser(y); l = lg(y);
  if (v < 0) pari_err_DOMAIN("hypergeom","valuation", "<", gen_0, y);
  if (!v)
  {
    y0 = gel(y, 2);
    if (!is_scalar_t(typ(y0))) pari_err_TYPE("hypergeom",y);
    y = serchop0(y); l = 3 + (l - 3) / valser(y);
    S = hypergeom(N, D, y0, prec);
  }
  av = avma; Ni = N; Di = D;
  for (i = 1; i < l; i++)
  {
    R = gmul(R, gmul(y, gdiv(vecprod(Ni), gmulsg(i, vecprod(Di)))));
    /* have to offset by 1 when y0 != NULL; keep Ni,Di to avoid recomputing
     * in next loop */
    Ni = RgV_z_add(N, i);
    Di = RgV_z_add(D, i);
    S = gadd(S, y0? gmul(R, hypergeom_i(Ni, Di, y0, prec)): R);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"hypergeom, i = %ld / %ld", i,l-1);
      gerepileall(av, 4, &S, &R, &Ni, &Di);
    }
  }
  return S;
}

GEN
hypergeom(GEN N, GEN D, GEN y, long prec)
{
  pari_sp av = avma;
  GEN z;
  long j, n;
  N = hypergeom_arg(N);
  D = hypergeom_arg(D);
  n = hypersimplify(&N, &D);
  for (j = 1; j < lg(D); j++)
    if (isnegint(gel(D,j)))
      pari_err_DOMAIN("hypergeom", stack_sprintf("b[%ld]", j + n),
                      "<=", gen_0, gel(D,j));
  if (is_scalar_t(typ(y)))
    return gerepilecopy(av, hypergeom_i(N, D, y, prec));
  if (!(z = toser_i(y))) pari_err_TYPE("hypergeom", y);
  return gerepileupto(av, serhypergeom(N, D, z, prec));
}
