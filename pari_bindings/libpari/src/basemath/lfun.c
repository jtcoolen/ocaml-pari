/* Copyright (C) 2015  The PARI group.

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
/**                       L-functions                              **/
/**                                                                **/
/********************************************************************/

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_lfun

/*******************************************************************/
/*  Accessors                                                      */
/*******************************************************************/

static GEN
mysercoeff(GEN x, long n)
{
  long N = n - valser(x);
  return (N < 0)? gen_0: gel(x, N+2);
}

long
ldata_get_type(GEN ldata) { return mael3(ldata, 1, 1, 1); }

GEN
ldata_get_an(GEN ldata) { return gel(ldata, 1); }

GEN
ldata_get_dual(GEN ldata) { return gel(ldata, 2); }

long
ldata_isreal(GEN ldata) { return isintzero(gel(ldata, 2)); }

GEN
ldata_get_gammavec(GEN ldata) { return gel(ldata, 3); }

long
ldata_get_degree(GEN ldata) { return lg(gel(ldata, 3))-1; }

GEN
ldata_get_k(GEN ldata)
{
  GEN w = gel(ldata,4);
  if (typ(w) == t_VEC) w = gel(w,1);
  return w;
}

/* a_n = O(n^{k1 + epsilon}) */
GEN
ldata_get_k1(GEN ldata)
{
  GEN w = gel(ldata,4);
  if (typ(w) == t_VEC) return gel(w,2);
  /* by default, assume that k1 = k-1 and even (k-1)/2 for entire functions */
  w = gaddgs(w,-1);
  return ldata_get_residue(ldata)? w: gmul2n(w, -1);
}

/* a_n = O(n^{k1 + epsilon}) */
static double
ldata_get_k1_dbl(GEN ldata)
{
  GEN w = gel(ldata,4);
  double k;
  if (typ(w) == t_VEC) return gtodouble(gel(w,2));
  /* by default, assume that k1 = k-1 and even (k-1)/2 for entire functions */
  k = gtodouble(w);
  return ldata_get_residue(ldata)? k-1: (k-1)/2.;
}

GEN
ldata_get_conductor(GEN ldata) { return gel(ldata, 5); }

GEN
ldata_get_rootno(GEN ldata) { return gel(ldata, 6); }

GEN
ldata_get_residue(GEN ldata) { return lg(ldata) == 7 ? NULL: gel(ldata, 7); }

long
linit_get_type(GEN linit) { return mael(linit, 1, 1); }

GEN
linit_get_ldata(GEN linit) { return gel(linit, 2); }

GEN
linit_get_tech(GEN linit) { return gel(linit, 3); }

long
is_linit(GEN data)
{
  return lg(data) == 4 && typ(data) == t_VEC
                       && typ(gel(data, 1)) == t_VECSMALL;
}

GEN
lfun_get_step(GEN tech) { return gmael(tech, 2, 1);}

GEN
lfun_get_pol(GEN tech) { return gmael(tech, 2, 2);}

GEN
lfun_get_Residue(GEN tech) { return gmael(tech, 2, 3);}

GEN
lfun_get_k2(GEN tech) { return gmael(tech, 3, 1);}

GEN
lfun_get_w2(GEN tech) { return gmael(tech, 3, 2);}

GEN
lfun_get_expot(GEN tech) { return gmael(tech, 3, 3);}

GEN
lfun_get_factgammavec(GEN tech) { return gmael(tech, 3, 4); }

/* Handle complex Vga whose sum is real */
static GEN
sumVga(GEN Vga) { return real_i(vecsum(Vga)); }
/* sum_i max (Im v[i],0) */
static double
sumVgaimpos(GEN v)
{
  double d = 0.;
  long i, l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN c = imag_i(gel(v,i));
    if (gsigne(c) > 0) d += gtodouble(c);
  }
  return d;
}

static long
vgaell(GEN Vga)
{
  if (lg(Vga) == 3)
  { GEN c = gsub(gel(Vga,1), gel(Vga,2)); return gequal1(c) || gequalm1(c); }
  return 0;
}
int
Vgaeasytheta(GEN Vga) { return lg(Vga)-1 == 1 || vgaell(Vga); }
/* return b(n) := a(n) * n^c, when Vgaeasytheta(Vga) is set */
static GEN
antwist(GEN an, GEN Vga, long prec)
{
  long l, i;
  GEN b, c = vecmin(Vga);
  if (gequal0(c)) return an;
  l = lg(an); b = cgetg(l, t_VEC);
  if (gequal1(c))
  {
    if (typ(an) == t_VECSMALL)
      for (i = 1; i < l; i++) gel(b,i) = mulss(an[i], i);
    else
      for (i = 1; i < l; i++) gel(b,i) = gmulgu(gel(an,i), i);
  }
  else
  {
    GEN v = vecpowug(l-1, c, prec);
    if (typ(an) == t_VECSMALL)
      for (i = 1; i < l; i++) gel(b,i) = gmulsg(an[i], gel(v,i));
    else
      for (i = 1; i < l; i++) gel(b,i) = gmul(gel(an,i), gel(v,i));
  }
  return b;
}

static GEN
theta_dual(GEN theta, GEN bn)
{
  if (typ(bn)==t_INT) return NULL;
  else
  {
    GEN thetad = shallowcopy(theta), ldata = linit_get_ldata(theta);
    GEN Vga = ldata_get_gammavec(ldata);
    GEN tech = shallowcopy(linit_get_tech(theta));
    GEN an = theta_get_an(tech);
    long prec = nbits2prec(theta_get_bitprec(tech));
    GEN vb = ldata_vecan(bn, lg(an)-1, prec);
    if (!theta_get_m(tech) && Vgaeasytheta(Vga)) vb = antwist(vb, Vga, prec);
    gel(tech,1) = vb;
    gel(thetad,3) = tech; return thetad;
  }
}

static GEN
domain_get_dom(GEN domain)  { return gel(domain,1); }
static long
domain_get_der(GEN domain)  { return mael2(domain, 2, 1); }
static long
domain_get_bitprec(GEN domain)  { return mael2(domain, 2, 2); }
GEN
lfun_get_domain(GEN tech) { return gel(tech,1); }
long
lfun_get_bitprec(GEN tech){ return domain_get_bitprec(lfun_get_domain(tech)); }
GEN
lfun_get_dom(GEN tech) { return domain_get_dom(lfun_get_domain(tech)); }

GEN
lfunprod_get_fact(GEN tech)  { return gel(tech, 2); }

GEN
theta_get_an(GEN tdata)      { return gel(tdata, 1);}
GEN
theta_get_K(GEN tdata)       { return gel(tdata, 2);}
GEN
theta_get_R(GEN tdata)       { return gel(tdata, 3);}
long
theta_get_bitprec(GEN tdata) { return itos(gel(tdata, 4));}
long
theta_get_m(GEN tdata)       { return itos(gel(tdata, 5));}
GEN
theta_get_tdom(GEN tdata)    { return gel(tdata, 6);}
GEN
theta_get_isqrtN(GEN tdata)  { return gel(tdata, 7);}

/*******************************************************************/
/*  Helper functions related to Gamma products                     */
/*******************************************************************/
/* x != 0 */
static int
serisscalar(GEN x)
{
  long i;
  if (valser(x)) return 0;
  for (i = lg(x)-1; i > 3; i--) if (!gequal0(gel(x,i))) return 0;
  return 1;
}

/* return -itos(s) >= 0 if scalar s is (approximately) equal to a nonpositive
 * integer, and -1 otherwise */
static long
isnegint(GEN s)
{
  GEN r = ground(real_i(s));
  if (signe(r) <= 0 && gequal(s, r)) return -itos(r);
  return -1;
}
/* if s = a + O(x^n), a <= 0 integer, replace by a + b*x^n + O(x^(n+1)) */
static GEN
serextendifnegint(GEN s, GEN b, long *ext)
{
  if (!signe(s) || (serisscalar(s) && isnegint(gel(s,2)) >= 0))
  {
    long l = lg(s);
    GEN t = cgetg(l+1, t_SER);
    gel(t, l) = b; while (--l > 1) gel(t,l) = gel(s,l);
    if (gequal0(gel(t,2))) gel(t,2) = gen_0;
    t[1] = s[1]; s = normalizeser(t); *ext = 1;
  }
  return s;
}

/* r/x + O(1), r != 0 */
static GEN
serpole(GEN r)
{
  GEN s = cgetg(3, t_SER);
  s[1] = evalsigne(1)|evalvalser(-1)|evalvarn(0);
  gel(s,2) = r; return s;
}
/* a0 +  a1 x + O(x^e), e >= 0 */
static GEN
deg1ser_shallow(GEN a1, GEN a0, long v, long e)
{ return RgX_to_ser(deg1pol_shallow(a1, a0, v), e+2); }

/* pi^(-s/2) Gamma(s/2) */
static GEN
gamma_R(GEN s, long *ext, long prec)
{
  GEN s2 = gmul2n(s, -1);
  long ms;

  if (typ(s) == t_SER)
    s2 = serextendifnegint(s2, ghalf, ext);
  else if ((ms = isnegint(s2)) >= 0)
  {
    GEN r = gmul(powPis(stoi(ms),prec), gdivsg(odd(ms)? -2: 2, mpfact(ms)));
    return serpole(r);
  }
  return gdiv(ggamma(s2,prec), powPis(s2,prec));
}
/* gamma_R(s)gamma_R(s+1) = 2 (2pi)^(-s) Gamma(s) */
static GEN
gamma_C(GEN s, long *ext, long prec)
{
  long ms;
  if (typ(s) == t_SER)
    s = serextendifnegint(s, gen_1, ext);
  else if ((ms = isnegint(s)) >= 0)
  {
    GEN r = gmul(pow2Pis(stoi(ms),prec), gdivsg(odd(ms)? -2: 2, mpfact(ms)));
    return serpole(r);
  }
  return gmul2n(gdiv(ggamma(s,prec), pow2Pis(s,prec)), 1);
}

static GEN
gammafrac(GEN r, long d)
{
  long i, l = labs(d) + 1, j = (d > 0)? 0: 2*d;
  GEN T, v = cgetg(l, t_COL);
  for (i = 1; i < l; i++, j += 2)
    gel(v,i) = deg1pol_shallow(gen_1, gaddgs(r, j), 0);
  T = RgV_prod(v); return d > 0? T: mkrfrac(gen_1, T);
}

/*
GR(s)=Pi^-(s/2)*gamma(s/2);
GC(s)=2*(2*Pi)^-s*gamma(s)
gdirect(F,s)=prod(i=1,#F,GR(s+F[i]))
gfact(F,s)=
{ my([R,A,B]=gammafactor(F), [a,e]=A, [b,f]=B, p=poldegree(R));
  subst(R,x,s) * (2*Pi)^-p * prod(i=1,#a,GR(s+a[i])^e[i])
                           * prod(i=1,#b,GC(s+b[i])^f[i]); }
*/
static GEN
gammafactor(GEN Vga)
{
  long i, r, c, l = lg(Vga);
  GEN v, P, a, b, e, f, E, F = cgetg(l, t_VEC), R = gen_1;
  for (i = 1; i < l; ++i)
  {
    GEN a = gel(Vga,i), r = gmul2n(real_i(a), -1);
    long q = itos(gfloor(r)); /* [Re a/2] */
    r = gmul2n(gsubgs(r, q), 1);
    gel(F,i) = gequal0(imag_i(a)) ? r : mkcomplex(r, imag_i(a)); /* 2{Re a/2} + I*(Im a) */
    if (q) R = gmul(R, gammafrac(gel(F,i), q));
  }
  F = vec_reduce(F, &E); l = lg(E);
  v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
      gel(v,i) = mkvec2(gsub(gel(F,i),gfloor(real_i(gel(F,i)))), stoi(E[i]));
  gen_sort_inplace(v, (void*)cmp_universal, cmp_nodata, &P);
  a = cgetg(l, t_VEC); e = cgetg(l, t_VECSMALL);
  b = cgetg(l, t_VEC); f = cgetg(l, t_VECSMALL);
  for (i = r = c = 1; i < l;)
    if (i==l-1 || cmp_universal(gel(v,i), gel(v,i+1)))
    { gel(a, r) = gel(F, P[i]); e[r++] = E[P[i]]; i++; }
    else
    { gel(b, c) = gel(F, P[i]); f[c++] = E[P[i]]; i+=2; }
  setlg(a, r); setlg(e, r);
  setlg(b, c); setlg(f, c); return mkvec3(R, mkvec2(a,e), mkvec2(b,f));
}

static GEN
polgammaeval(GEN F, GEN s)
{
  GEN r = poleval(F, s);
  if (typ(s) != t_SER && gequal0(r))
  { /* here typ(F) = t_POL */
    long e;
    for (e = 1;; e++)
    {
      F = RgX_deriv(F); r = poleval(F,s);
      if (!gequal0(r)) break;
    }
    if (e > 1) r = gdiv(r, mpfact(e));
    r = serpole(r); setvalser(r, e);
  }
  return r;
}
static long
rfrac_degree(GEN R)
{
  GEN a = gel(R,1), b = gel(R,2);
  return ((typ(a) == t_POL)? degpol(a): 0) - degpol(b);
}
static GEN
fracgammaeval(GEN F, GEN s, long prec)
{
  GEN R = gel(F,1);
  long d;
  switch(typ(R))
  {
    case t_POL:
      d = degpol(R);
      R = polgammaeval(R, s); break;
    case t_RFRAC:
      d = rfrac_degree(R);
      R = gdiv(polgammaeval(gel(R,1), s), polgammaeval(gel(R,2), s)); break;
    default: return R;
  }
  return gmul(R, powrs(Pi2n(1,prec), -d));
}

static GEN
gammafactproduct(GEN F, GEN s, long *ext, long prec)
{
  pari_sp av = avma;
  GEN R = gel(F,2), Rw = gel(R,1), Re = gel(R,2);
  GEN C = gel(F,3), Cw = gel(C,1), Ce = gel(C,2), z = fracgammaeval(F,s,prec);
  long i, lR = lg(Rw), lC = lg(Cw);
  *ext = 0;
  for (i = 1; i < lR; i++)
    z = gmul(z, gpowgs(gamma_R(gadd(s,gel(Rw, i)), ext, prec), Re[i]));
  for (i = 1; i < lC; i++)
    z = gmul(z, gpowgs(gamma_C(gadd(s,gel(Cw, i)), ext, prec), Ce[i]));
  return gerepileupto(av, z);
}

static int
gammaordinary(GEN Vga, GEN s)
{
  long i, d = lg(Vga)-1;
  for (i = 1; i <= d; i++)
  {
    GEN z = gadd(s, gel(Vga,i));
    long e;
    if (gexpo(imag_i(z)) < -10)
    {
      z = real_i(z);
      if (gsigne(z) <= 0) { (void)grndtoi(z, &e); if (e < -10) return 0; }
    }
  }
  return 1;
}

/* Exponent A of t in asymptotic expansion; K(t) ~ C t^A exp(-pi d t^(2/d)).
 * suma = vecsum(Vga)*/
static double
gammavec_expo(long d, double suma) { return (1 - d + suma) / d; }

/*******************************************************************/
/*       First part: computations only involving Theta(t)          */
/*******************************************************************/

static void
get_cone(GEN t, double *r, double *a)
{
  const long prec = LOWDEFAULTPREC;
  if (typ(t) == t_COMPLEX)
  {
    t  = gprec_w(t, prec);
    *r = gtodouble(gabs(t, prec));
    *a = fabs(gtodouble(garg(t, prec)));
  }
  else
  {
    *r = fabs(gtodouble(t));
    *a = 0.;
  }
  if (!*r && !*a) pari_err_DOMAIN("lfunthetainit","t","=",gen_0,t);
}
/* slightly larger cone than necessary, to avoid round-off problems */
static void
get_cone_fuzz(GEN t, double *r, double *a)
{ get_cone(t, r, a); *r -= 1e-10; if (*a) *a += 1e-10; }

/* Initialization m-th Theta derivative. tdom is either
 * - [rho,alpha]: assume |t| >= rho and |arg(t)| <= alpha
 * - a positive real scalar: assume t real, t >= tdom;
 * - a complex number t: compute at t;
 * N is the conductor (either the true one from ldata or a guess from
 * lfunconductor) */
long
lfunthetacost(GEN ldata, GEN tdom, long m, long bitprec)
{
  pari_sp av = avma;
  GEN Vga = ldata_get_gammavec(ldata);
  long d = lg(Vga)-1;
  double k1 = maxdd(ldata_get_k1_dbl(ldata), 0.);
  double c = d/2., a, A, B, logC, al, rho, T;
  double N = gtodouble(ldata_get_conductor(ldata));

  if (!N) pari_err_TYPE("lfunthetaneed [missing conductor]", ldata);
  if (typ(tdom) == t_VEC && lg(tdom) == 3)
  {
    rho= gtodouble(gel(tdom,1));
    al = gtodouble(gel(tdom,2));
  }
  else
    get_cone_fuzz(tdom, &rho, &al);
  A = gammavec_expo(d, gtodouble(sumVga(Vga))); set_avma(av);
  a = (A+k1+1) + (m-1)/c;
  if (fabs(a) < 1e-10) a = 0.;
  logC = c*M_LN2 - log(c)/2;
  /* +1: fudge factor */
  B = M_LN2*bitprec+logC+m*log(2*M_PI) + 1 + (k1+1)*log(N)/2 - (k1+m+1)*log(rho);
  if (al)
  { /* t = rho e^(i*al), T^(1/c) = Re(t^(1/c)) > 0, T = rho cos^c(al/c) */
    double z = cos(al/c);
    if (z <= 0)
      pari_err_DOMAIN("lfunthetaneed", "arg t", ">", dbltor(c*M_PI/2), tdom);
    T = (d == 2 && typ(tdom) != t_VEC)? gtodouble(real_i(tdom)): rho*pow(z,c);
    B -= log(z) * (c * (k1+A+1) + m);
  }
  else
    T = rho;
  if (B <= 0) return 0;
  A = floor(0.9 + dblcoro526(a,c,B) / T * sqrt(N));
  if (dblexpo(A) >= BITS_IN_LONG-1) pari_err_OVERFLOW("lfunthetacost");
  return (long)A;
}
long
lfunthetacost0(GEN L, GEN tdom, long m, long bitprec)
{
  long n;
  if (is_linit(L) && linit_get_type(L)==t_LDESC_THETA)
  {
    GEN tech = linit_get_tech(L);
    n = lg(theta_get_an(tech))-1;
  }
  else
  {
    pari_sp av = avma;
    GEN ldata = lfunmisc_to_ldata_shallow(L);
    n = lfunthetacost(ldata, tdom? tdom: gen_1, m, bitprec);
    set_avma(av);
  }
  return n;
}

static long
fracgammadegree(GEN FVga)
{ GEN F = gel(FVga,1); return (typ(F)==t_RFRAC)? degpol(gel(F,2)): 0; }

/* Poles of a L-function can be represented in the following ways:
 * 1) Nothing (ldata has only 6 components, ldata_get_residue = NULL).
 * 2) a complex number (single pole at s = k with given residue, unknown if 0).
 * 3) A vector (possibly empty) of 2-component vectors [a, ra], where a is the
 * pole, ra a t_SER: its Taylor expansion at a. A t_VEC encodes the polar
 * part of L, a t_COL, the polar part of Lambda */

/* 'a' a complex number (pole), 'r' the polar part of L at 'a';
 * return 'R' the polar part of Lambda at 'a' */
static GEN
rtoR(GEN a, GEN r, GEN FVga, GEN N, long prec)
{
  long v = lg(r)-2, d = fracgammadegree(FVga), ext;
  GEN Na, as = deg1ser_shallow(gen_1, a, varn(r), v);
  Na = gpow(N, gdivgu(as, 2), prec);
  /* make up for a possible loss of accuracy */
  if (d) as = deg1ser_shallow(gen_1, a, varn(r), v + d);
  return gmul(gmul(r, Na), gammafactproduct(FVga, as, &ext, prec));
}

/* assume r in normalized form: t_VEC of pairs [be,re] */
GEN
lfunrtopoles(GEN r)
{
  long j, l = lg(r);
  GEN v = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
  {
    GEN rj = gel(r,j), a = gel(rj,1);
    gel(v,j) = a;
  }
  gen_sort_inplace(v, (void*)&cmp_universal, cmp_nodata, NULL);
  return v;
}

/* r / x + O(1) */
static GEN
simple_pole(GEN r)
{ return isintzero(r)? gen_0: serpole(r); }
static GEN
normalize_simple_pole(GEN r, GEN k)
{
  long tx = typ(r);
  if (is_vec_t(tx)) return r;
  if (!is_scalar_t(tx)) pari_err_TYPE("lfunrootres [poles]", r);
  return mkvec(mkvec2(k, simple_pole(r)));
}
/* normalize the description of a polar part */
static GEN
normalizepoles(GEN r, GEN k)
{
  long iv, j, l;
  GEN v;
  if (!is_vec_t(typ(r))) return normalize_simple_pole(r, k);
  v = cgetg_copy(r, &l);
  for (j = iv = 1; j < l; j++)
  {
    GEN rj = gel(r,j), a = gel(rj,1), ra = gel(rj,2);
    if (!is_scalar_t(typ(a)) || typ(ra) != t_SER)
      pari_err_TYPE("lfunrootres [poles]",r);
    if (valser(ra) >= 0) continue;
    gel(v,iv++) = rj;
  }
  setlg(v, iv); return v;
}
static int
residues_known(GEN r)
{
  long i, l = lg(r);
  if (isintzero(r)) return 0;
  if (!is_vec_t(typ(r))) return 1;
  for (i = 1; i < l; i++)
  {
    GEN ri = gel(r,i);
    if (!is_vec_t(typ(ri)) || lg(ri)!=3)
      pari_err_TYPE("lfunrootres [poles]",r);
    if (isintzero(gel(ri, 2))) return 0;
  }
  return 1;
}

/* Compute R's from r's (r = Taylor devts of L(s), R of Lambda(s)).
 * 'r/eno' passed to override the one from ldata  */
static GEN
lfunrtoR_i(GEN ldata, GEN r, GEN eno, long prec)
{
  GEN Vga = ldata_get_gammavec(ldata), N = ldata_get_conductor(ldata);
  GEN R, vr, FVga;
  pari_sp av = avma;
  long lr, j, jR;
  GEN k = ldata_get_k(ldata);

  if (!r || isintzero(eno) || !residues_known(r))
    return gen_0;
  r = normalizepoles(r, k);
  if (typ(r) == t_COL) return gerepilecopy(av, r);
  if (typ(ldata_get_dual(ldata)) != t_INT)
    pari_err(e_MISC,"please give the Taylor development of Lambda");
  vr = lfunrtopoles(r); lr = lg(vr);
  FVga = gammafactor(Vga);
  R = cgetg(2*lr, t_COL);
  for (j = jR = 1; j < lr; j++)
  {
    GEN rj = gel(r,j), a = gel(rj,1), ra = gel(rj,2);
    GEN Ra = rtoR(a, ra, FVga, N, prec);
    GEN b = gsub(k, conj_i(a));
    if (lg(Ra)-2 < -valser(Ra))
      pari_err(e_MISC,
        "please give more terms in L function's Taylor development at %Ps", a);
    gel(R,jR++) = mkvec2(a, Ra);
    if (!tablesearch(vr, b, (int (*)(GEN,GEN))&cmp_universal))
    {
      GEN mX = gneg(pol_x(varn(Ra)));
      GEN Rb = gmul(eno, gsubst(conj_i(Ra), varn(Ra), mX));
      gel(R,jR++) = mkvec2(b, Rb);
    }
  }
  setlg(R, jR); return gerepilecopy(av, R);
}
static GEN
lfunrtoR_eno(GEN ldata, GEN eno, long prec)
{ return lfunrtoR_i(ldata, ldata_get_residue(ldata), eno, prec); }
static GEN
lfunrtoR(GEN ldata, long prec)
{ return lfunrtoR_eno(ldata, ldata_get_rootno(ldata), prec); }

static long
prec_fix(long prec)
{
#ifndef LONG_IS_64BIT
  /* make sure that default accuracy is the same on 32/64bit */
  if (odd(prec)) prec += EXTRAPREC64;
#endif
  return prec;
}

/* thetainit using {an: n <= L}; if (m = 0 && easytheta), an2 is an * n^al */
static GEN
lfunthetainit0(GEN ldata, GEN tdom, GEN an2, long m,
    long bitprec, long extrabit)
{
  long prec = nbits2prec(bitprec);
  GEN tech, N = ldata_get_conductor(ldata);
  GEN K = gammamellininvinit(ldata, m, bitprec + extrabit);
  GEN R = lfunrtoR(ldata, prec);
  if (!tdom) tdom = gen_1;
  if (typ(tdom) != t_VEC)
  {
    double r, a;
    get_cone_fuzz(tdom, &r, &a);
    tdom = mkvec2(dbltor(r), a? dbltor(a): gen_0);
  }
  prec += maxss(EXTRAPREC64, nbits2extraprec(extrabit));
  tech = mkvecn(7, an2,K,R, stoi(bitprec), stoi(m), tdom,
                   gsqrt(ginv(N), prec_fix(prec)));
  return mkvec3(mkvecsmall(t_LDESC_THETA), ldata, tech);
}

/* tdom: 1) positive real number r, t real, t >= r; or
 *       2) [r,a], describing the cone |t| >= r, |arg(t)| <= a */
static GEN
lfunthetainit_i(GEN data, GEN tdom, long m, long bit)
{
  GEN ldata = lfunmisc_to_ldata_shallow(data);
  long b = 32, L = lfunthetacost(ldata, tdom, m, bit), prec = nbits2prec(bit);
  GEN ldatan = ldata_newprec(ldata, prec);
  GEN an = ldata_vecan(ldata_get_an(ldatan), L, prec);
  GEN Vga = ldata_get_gammavec(ldatan);
  if (m == 0 && Vgaeasytheta(Vga)) an = antwist(an, Vga, prec);
  if (typ(an) != t_VECSMALL) b = maxss(b, gexpo(an));
  return lfunthetainit0(ldatan, tdom, an, m, bit, b);
}

GEN
lfunthetainit(GEN ldata, GEN tdom, long m, long bitprec)
{
  pari_sp av = avma;
  GEN S = lfunthetainit_i(ldata, tdom? tdom: gen_1, m, bitprec);
  return gerepilecopy(av, S);
}

GEN
lfunan(GEN ldata, long L, long prec)
{
  pari_sp av = avma;
  GEN an ;
  ldata = ldata_newprec(lfunmisc_to_ldata_shallow(ldata), prec);
  an = gerepilecopy(av, ldata_vecan(ldata_get_an(ldata), L, prec));
  if (typ(an) != t_VEC) an = vecsmall_to_vec_inplace(an);
  return an;
}

static GEN
mulrealvec(GEN x, GEN y)
{
  if (is_vec_t(typ(x)) && is_vec_t(typ(y)))
    pari_APPLY_same(mulreal(gel(x,i),gel(y,i)))
  else
    return mulreal(x,y);
}
static GEN
gmulvec(GEN x, GEN y)
{
  if (is_vec_t(typ(x)) && is_vec_t(typ(y)))
    pari_APPLY_same(gmul(gel(x,i),gel(y,i)))
  else
    return gmul(x,y);
}
static GEN
gdivvec(GEN x, GEN y)
{
  if (is_vec_t(typ(x)) && is_vec_t(typ(y)))
    pari_APPLY_same(gdiv(gel(x,i),gel(y,i)))
  else
    return gdiv(x,y);
}

static GEN
gsubvec(GEN x, GEN y)
{
  if (is_vec_t(typ(x)) && !is_vec_t(typ(y)))
    pari_APPLY_same(gsub(gel(x,i),y))
  else
    return gsub(x,y);
}

/* return [1^(2/d), 2^(2/d),...,lim^(2/d)] */
static GEN
mkvroots(long d, long lim, long prec)
{
  if (d <= 4)
  {
    GEN v = cgetg(lim+1,t_VEC);
    long n;
    switch(d)
    {
      case 1:
        for (n=1; n <= lim; n++) gel(v,n) = sqru(n);
        return v;
      case 2:
        for (n=1; n <= lim; n++) gel(v,n) = utoipos(n);
        return v;
      case 4:
        for (n=1; n <= lim; n++) gel(v,n) = sqrtr(utor(n, prec));
        return v;
    }
  }
  return vecpowug(lim, gdivgu(gen_2,d), prec);
}

GEN
lfunthetacheckinit(GEN data, GEN t, long m, long bitprec)
{
  if (is_linit(data) && linit_get_type(data)==t_LDESC_THETA)
  {
    GEN tdom, thetainit = linit_get_tech(data);
    long bitprecnew = theta_get_bitprec(thetainit);
    long m0 = theta_get_m(thetainit);
    double r, al, rt, alt;
    if (m0 != m)
      pari_err_DOMAIN("lfuntheta","derivative order","!=", stoi(m),stoi(m0));
    if (bitprec > bitprecnew) goto INIT;
    get_cone(t, &rt, &alt);
    tdom = theta_get_tdom(thetainit);
    r = gtodouble(gel(tdom,1));
    al= gtodouble(gel(tdom,2)); if (rt >= r && alt <= al) return data;
  }
INIT:
  return lfunthetainit_i(data, t, m, bitprec);
}

static GEN
get_an(GEN an, long n)
{
  if (typ(an) == t_VECSMALL) { long a = an[n]; if (a) return stoi(a); }
  else { GEN a = gel(an,n); if (a && !gequal0(a)) return a; }
  return NULL;
}
/* x * an[n] */
static GEN
mul_an(GEN an, long n, GEN x)
{
  if (typ(an) == t_VECSMALL) { long a = an[n]; if (a) return gmulsg(a,x); }
  else { GEN a = gel(an,n); if (a && !gequal0(a)) return gmul(a,x); }
  return NULL;
}
/* 2*t^a * x **/
static GEN
mulT(GEN t, GEN a, GEN x, long prec)
{
  if (gequal0(a)) return gmul2n(x,1);
  return gmul(x, gmul2n(gequal1(a)? t: gpow(t,a,prec), 1));
}

static GEN
vecan_cmul(void *E, GEN P, long a, GEN x)
{
  (void)E;
  if (typ(P) == t_VECSMALL)
    return (a==0 || !P[a])? NULL: gmulsg(P[a], x);
  else
    return (a==0 || !gel(P,a))? NULL: gmul(gel(P,a), x);
}
/* d=2, 2 sum_{n <= N} a(n) (n t)^al q^n, q = exp(-2pi t),
 * an2[n] = a(n) * n^al */
static GEN
theta2_i(GEN an2, long N, GEN t, GEN al, long prec)
{
  GEN S, q, pi2 = Pi2n(1,prec);
  const struct bb_algebra *alg = get_Rg_algebra();
  setsigne(pi2,-1); q = gexp(gmul(pi2, t), prec);
  /* Brent-Kung in case the a_n are small integers */
  S = gen_bkeval(an2, N, q, 1, NULL, alg, vecan_cmul);
  return mulT(t, al, S, prec);
}
static GEN
theta2(GEN an2, long N, GEN t, GEN al, long prec)
{
  pari_sp av = avma;
  return gerepileupto(av, theta2_i(an2, N, t, al, prec));
}

/* d=1, 2 sum_{n <= N} a_n (n t)^al q^(n^2), q = exp(-pi t^2),
 * an2[n] is a_n n^al */
static GEN
theta1(GEN an2, long N, GEN t, GEN al, long prec)
{
  GEN q = gexp(gmul(negr(mppi(prec)), gsqr(t)), prec);
  GEN vexp = gsqrpowers(q, N), S = gen_0;
  pari_sp av = avma;
  long n;
  for (n = 1; n <= N; n++)
  {
    GEN c = mul_an(an2, n, gel(vexp,n));
    if (c)
    {
      S = gadd(S, c);
      if (gc_needed(av, 3)) S = gerepileupto(av, S);
    }
  }
  return mulT(t, al, S, prec);
}

/* If m > 0, compute m-th derivative of theta(t) = theta0(t/sqrt(N))
 * with absolute error 2^-bitprec; theta(t)=\sum_{n\ge1}a(n)K(nt/N^(1/2)) */
GEN
lfuntheta(GEN data, GEN t, long m, long bitprec)
{
  pari_sp ltop = avma;
  long limt, d;
  GEN isqN, vecan, Vga, ldata, theta, thetainit, S;
  long n, prec;

  theta = lfunthetacheckinit(data, t, m, bitprec);
  ldata = linit_get_ldata(theta);
  thetainit = linit_get_tech(theta);
  vecan = theta_get_an(thetainit);
  isqN = theta_get_isqrtN(thetainit);
  prec = maxss(realprec(isqN), nbits2prec(bitprec));
  t = gprec_w(t, prec);
  limt = lg(vecan)-1;
  if (theta == data)
    limt = minss(limt, lfunthetacost(ldata, t, m, bitprec));
  if (!limt)
  {
    set_avma(ltop); S = real_0_bit(-bitprec);
    if (!is_real_t(typ(t)) || !ldata_isreal(ldata))
      S = gerepilecopy(ltop, mkcomplex(S,S));
    return S;
  }
  t = gmul(t, isqN);
  Vga = ldata_get_gammavec(ldata);
  d = lg(Vga)-1;
  if (m == 0 && Vgaeasytheta(Vga))
  {
    if (theta_get_m(thetainit) > 0) vecan = antwist(vecan, Vga, prec);
    if (d == 1) S = theta1(vecan, limt, t, gel(Vga,1), prec);
    else        S = theta2_i(vecan, limt, t, vecmin(Vga), prec);
  }
  else
  {
    GEN K = theta_get_K(thetainit);
    GEN vroots = mkvroots(d, limt, prec);
    pari_sp av;
    t = gpow(t, gdivgu(gen_2,d), prec);
    S = gen_0; av = avma;
    for (n = 1; n <= limt; ++n)
    {
      GEN nt, an = get_an(vecan, n);
      if (!an) continue;
      nt = gmul(gel(vroots,n), t);
      if (m) an = gmul(an, powuu(n, m));
      S = gadd(S, gmul(an, gammamellininvrt(K, nt, bitprec)));
      if ((n & 0x1ff) == 0) S = gerepileupto(av, S);
    }
    if (m) S = gmul(S, gpowgs(isqN, m));
  }
  return gerepileupto(ltop, S);
}

/*******************************************************************/
/* Second part: Computation of L-Functions.                        */
/*******************************************************************/

struct lfunp {
  long precmax, Dmax, D, M, m0, nmax, d, vgaell;
  double k1, dc, dw, dh, MAXs, sub;
  GEN L, an, bn;
};

static void
lfunp_set(GEN ldata, long der, long bitprec, struct lfunp *S)
{
  const long derprec = (der > 1)? dbllog2(mpfact(der)): 0; /* log2(der!) */
  GEN Vga, N, L, k;
  long k1, d, m, M, flag, nmax;
  double a, A, E, hd, Ep, d2, suma, maxs, mins, sub, B0,B1;
  double logN2, logC, Lestimate, Mestimate;

  Vga = ldata_get_gammavec(ldata);
  S->d = d = lg(Vga)-1; d2 = d/2.;

  suma = gtodouble(sumVga(Vga));
  k = ldata_get_k(ldata);
  N = ldata_get_conductor(ldata);
  logN2 = log(gtodouble(N)) / 2;
  maxs = S->dc + S->dw;
  mins = S->dc - S->dw;
  S->MAXs = maxdd(maxs, gtodouble(k)-mins);

  /* we compute Lambda^(der)(s) / der!; need to compensate for L^(der)(s)
   * ln |gamma(s)| ~ -(pi/4) \sum_i |Im(s + a_i)|; max with 1: fudge factor */
  a = (M_PI/(4*M_LN2))*(d*S->dh + sumVgaimpos(Vga));
  S->D = (long)ceil(bitprec + derprec + maxdd(a, 1));
  E = M_LN2*S->D; /* D:= required absolute bitprec */

  Ep = E + maxdd(M_PI * S->dh * d2, (d*S->MAXs + suma - 1) * log(E));
  hd = d2*M_PI*M_PI / Ep;
  S->m0 = (long)ceil(M_LN2/hd);
  hd = M_LN2/S->m0;

  logC = d2*M_LN2 - log(d2)/2;
  k1 = maxdd(ldata_get_k1_dbl(ldata), 0.);
  S->k1 = k1; /* assume |a_n| << n^k1 with small implied constant */
  A = gammavec_expo(d, suma);

  sub = 0.;
  if (mins > 1)
  {
    GEN sig = dbltor(mins);
    sub += logN2*mins;
    if (gammaordinary(Vga, sig))
    {
      long ext;
      GEN gas = gammafactproduct(gammafactor(Vga), sig, &ext, LOWDEFAULTPREC);
      if (typ(gas) != t_SER)
      {
        double dg = dbllog2(gas);
        if (dg > 0) sub += dg * M_LN2;
      }
    }
  }
  S->sub = sub;
  M = 1000;
  L = cgetg(M+2, t_VECSMALL);
  a = S->k1 + A;

  B0 = 5 + E - S->sub + logC + S->k1*logN2; /* 5 extra bits */
  B1 = hd * (S->MAXs - S->k1);
  Lestimate = dblcoro526(a + S->MAXs - 2./d, d/2.,
    E - S->sub + logC - log(2*M_PI*hd) + S->MAXs*logN2);
  Mestimate = ((Lestimate > 0? log(Lestimate): 0) + logN2) / hd;
  nmax = 0;
  flag = 0;
  for (m = 0;; m++)
  {
    double x, H = logN2 - m*hd, B = B0 + m*B1;
    long n;
    x = dblcoro526(a, d/2., B);
    n = floor(x*exp(H));
    if (n > nmax) nmax = n;
    if (m > M) { M *= 2; L = vecsmall_lengthen(L,M+2); }
    L[m+1] = n;
    if (n == 0) { if (++flag > 2 && m > Mestimate) break; } else flag = 0;
  }
  m -= 2; while (m > 0 && !L[m]) m--;
  if (m == 0) { nmax = 1; L[1] = 1; m = 1; } /* can happen for tiny bitprec */
  setlg(L, m+1); S->M = m-1;
  S->L = L;
  S->nmax = nmax;

  S->Dmax = S->D + (long)ceil((S->M * hd * S->MAXs - S->sub) / M_LN2);
  if (S->Dmax < S->D) S->Dmax = S->D;
  S->precmax = nbits2prec(S->Dmax);
  if (DEBUGLEVEL > 1)
    err_printf("Dmax=%ld, D=%ld, M = %ld, nmax = %ld, m0 = %ld\n",
               S->Dmax,S->D,S->M,S->nmax, S->m0);
}

static GEN
lfuninit_pol(GEN v, GEN poqk, long prec)
{
  long m, M = lg(v) - 2;
  GEN pol = cgetg(M+3, t_POL);
  pol[1] = evalsigne(1) | evalvarn(0);
  gel(pol, 2) = gprec_w(gmul2n(gel(v,1), -1), prec);
  if (poqk)
    for (m = 2; m <= M+1; m++)
      gel(pol, m+1) = gprec_w(gmul(gel(poqk,m), gel(v,m)), prec);
  else
    for (m = 2; m <= M+1; m++)
      gel(pol, m+1) = gprec_w(gel(v,m), prec);
  return RgX_renormalize_lg(pol, M+3);
}

static void
worker_init(long q, GEN *an, GEN *bn, GEN *AB, GEN *A, GEN *B)
{
  if (typ(*bn) == t_INT) *bn = NULL;
  if (*bn)
  {
    *AB = cgetg(3, t_VEC);
    gel(*AB,1) = *A = cgetg(q+1, t_VEC);
    gel(*AB,2) = *B = cgetg(q+1, t_VEC);
    if (typ(an) == t_VEC) *an = RgV_kill0(*an);
    if (typ(bn) == t_VEC) *bn = RgV_kill0(*bn);
  }
  else
  {
    *B = NULL;
    *AB = *A = cgetg(q+1, t_VEC);
    if (typ(*an) == t_VEC) *an = RgV_kill0(*an);
  }
}
GEN
lfuninit_theta2_worker(long r, GEN L, GEN qk, GEN a, GEN di, GEN an, GEN bn)
{
  long q, m, prec = di[1], M = di[2], m0 = di[3], L0 = lg(an)-1;
  GEN AB, A, B;
  worker_init((M - r) / m0 + 1, &an, &bn, &AB, &A, &B);
  for (q = 0, m = r; m <= M; m += m0, q++)
  {
    GEN t = gel(qk, m+1);
    long N = minss(L[m+1],L0);
    gel(A, q+1) = theta2(an, N, t, a, prec); /* theta(exp(mh)) */
    if (bn) gel(B, q+1) = theta2(bn, N, t, a, prec);
  }
  return AB;
}

/* theta(exp(mh)) ~ sum_{n <= N} a(n) k[m,n] */
static GEN
an_msum(GEN an, long N, GEN vKm)
{
  pari_sp av = avma;
  GEN s = gen_0;
  long n;
  for (n = 1; n <= N; n++)
    if (gel(vKm,n))
    {
      GEN c = mul_an(an, n, gel(vKm,n));
      if (c) s = gadd(s, c);
    }
  return gerepileupto(av, s);
}

GEN
lfuninit_worker(long r, GEN K, GEN L, GEN peh2d, GEN vroots, GEN dr, GEN di,
                GEN an, GEN bn)
{
  pari_sp av0 = avma;
  long m, n, q, L0 = lg(an)-1;
  double sig0 = rtodbl(gel(dr,1)), sub2 = rtodbl(gel(dr,2));
  double k1 = rtodbl(gel(dr,3)), MAXs = rtodbl(gel(dr,4));
  long D = di[1], M = di[2], m0 = di[3];
  double M0 = sig0? sub2 / sig0: 1./0.;
  GEN AB, A, B, vK = cgetg(M/m0 + 2, t_VEC);

  for (q = 0, m = r; m <= M; m += m0, q++)
    gel(vK, q+1) = const_vec(L[m+1], NULL);
  worker_init(q, &an, &bn, &AB, &A, &B);
  for (m -= m0, q--; m >= 0; m -= m0, q--)
  {
    double c1 = D + ((m > M0)? m * sig0 - sub2 : 0);
    GEN vKm = gel(vK,q+1); /* conceptually K(m,n) */
    for (n = 1; n <= L[m+1]; n++)
    {
      GEN t2d, kmn;
      long nn, mm, qq, p = 0;
      double c, c2;
      pari_sp av;

      if (gel(vKm, n)) continue; /* done already */
      c = c1 + k1 * log2(n);
      /* n *= 2; m -= m0 => c += c2, provided m >= M0. Else c += k1 */
      c2 = k1 - MAXs;
      /* p = largest (absolute) accuracy to which we need K(m,n) */
      for (mm=m,nn=n; mm >= M0;)
      {
        if (nn <= L[mm+1] && (gel(an, nn) || (bn && gel(bn, nn))))
          if (c > 0) p = maxuu(p, (ulong)c);
        nn <<= 1;
        mm -= m0; if (mm >= M0) c += c2; else { c += k1; break; }
      }
      /* mm < M0 || nn > L[mm+1] */
      for (         ; mm >= 0; nn<<=1,mm-=m0,c+=k1)
        if (nn <= L[mm+1] && (gel(an, nn) || (bn && gel(bn, nn))))
          if (c > 0) p = maxuu(p, (ulong)c);
      if (!p) continue; /* a_{n 2^v} = 0 for all v in range */
      av = avma;
      t2d = mpmul(gel(vroots,n), gel(peh2d,m+1));/*(n exp(mh)/sqrt(N))^(2/d)*/
      kmn = gerepileupto(av, gammamellininvrt(K, t2d, p));
      for (qq=q,mm=m,nn=n; mm >= 0; nn<<=1,mm-=m0,qq--)
        if (nn <= L[mm+1]) gmael(vK, qq+1, nn) = kmn;
    }
  }
  for (q = 0, m = r; m <= M; m += m0, q++)
  {
    long N = minss(L0, L[m+1]);
    gel(A, q+1) = an_msum(an, N, gel(vK,q+1));
    if (bn) gel(B, q+1) = an_msum(bn, N, gel(vK,q+1));
  }
  return gerepileupto(av0, AB);
}
/* return A = [\theta(exp(mh)), m=0..M], theta(t) = sum a(n) K(n/sqrt(N) t),
 * h = log(2)/m0. If bn != NULL, return the pair [A, B] */
static GEN
lfuninit_ab(GEN theta, GEN h, struct lfunp *S)
{
  const long M = S->M, prec = S->precmax;
  GEN tech = linit_get_tech(theta), isqN = theta_get_isqrtN(tech);
  GEN an = S->an, bn = S->bn, va, vb;
  struct pari_mt pt;
  GEN worker;
  long m0, r, pending;

  if (S->vgaell)
  { /* d=2 and Vga = [a,a+1] */
    GEN a = vecmin(ldata_get_gammavec(linit_get_ldata(theta)));
    GEN qk = gpowers0(mpexp(h), M, isqN);
    m0 = minss(M+1, mt_nbthreads());
    worker = snm_closure(is_entry("_lfuninit_theta2_worker"),
                         mkvecn(6, S->L, qk, a, mkvecsmall3(prec, M, m0),
                                an, bn? bn: gen_0));
  }
  else
  {
    GEN vroots, peh2d, d2;
    double sig0 = S->MAXs / S->m0, sub2 = S->sub / M_LN2;
    /* For all 0<= m <= M, and all n <= L[m+1] such that a_n!=0, we compute
     *   k[m,n] = K(n exp(mh)/sqrt(N))
     * with ln(absolute error) <= E + max(mh sigma - sub, 0) + k1 * log(n).
     * N.B. we use the 'rt' variant and pass (n exp(mh)/sqrt(N))^(2/d).
     * Speedup: if n' = 2n and m' = m - m0 >= 0; then k[m,n] = k[m',n']. */
    vroots = mkvroots(S->d, S->nmax, prec); /* vroots[n] = n^(2/d) */
    d2 = gdivgu(gen_2, S->d);
    peh2d = gpowers0(gexp(gmul(d2,h), prec), M, gpow(isqN, d2, prec));
    m0 = S->m0; /* peh2d[m+1] = (exp(mh)/sqrt(N))^(2/d) */
    worker = snm_closure(is_entry("_lfuninit_worker"),
                         mkvecn(8, theta_get_K(tech), S->L, peh2d, vroots,
                                mkvec4(dbltor(sig0), dbltor(sub2),
                                       dbltor(S->k1), dbltor(S->MAXs)),
                                mkvecsmall3(S->D, M, m0),
                                an, bn? bn: gen_0));
    /* For each 0 <= m <= M, we will sum for n<=L[m+1] a(n) K(m,n)
     * bit accuracy for K(m,n): D + k1*log2(n) + 1_{m > M0} (m*sig0 - sub2)
     * We restrict m to arithmetic progressions r mod m0 to save memory and
     * allow parallelization */
  }
  va = cgetg(M+2, t_VEC);
  vb = bn? cgetg(M+2, t_VEC): NULL;
  mt_queue_start_lim(&pt, worker, m0);
  pending = 0;
  for (r = 0; r < m0 || pending; r++)
  { /* m = q m0 + r */
    GEN done, A, B;
    long q, m, workid;
    mt_queue_submit(&pt, r, r < m0 ? mkvec(utoi(r)): NULL);
    done = mt_queue_get(&pt, &workid, &pending);
    if (!done) continue;
    if (bn) { A = gel(done,1); B = gel(done,2); } else { A = done; B = NULL; }
    for (q = 0, m = workid; m <= M; m += m0, q++)
    {
      gel(va, m+1) = gel(A, q+1);
      if (bn) gel(vb, m+1) = gel(B, q+1);
    }
  }
  mt_queue_end(&pt);
  return bn? mkvec2(va, vb): va;
}

static void
parse_dom(double k, GEN dom, struct lfunp *S)
{
  long l = lg(dom);
  if (typ(dom)!=t_VEC) pari_err_TYPE("lfuninit [domain]", dom);
  if (l == 2)
  {
    S->dc = k/2.;
    S->dw = 0.;
    S->dh = gtodouble(gel(dom,1));
  }
  else if (l == 3)
  {
    S->dc = k/2.;
    S->dw = gtodouble(gel(dom,1));
    S->dh = gtodouble(gel(dom,2));
  }
  else if (l == 4)
  {
    S->dc = gtodouble(gel(dom,1));
    S->dw = gtodouble(gel(dom,2));
    S->dh = gtodouble(gel(dom,3));
  }
  else
  {
    pari_err_TYPE("lfuninit [domain]", dom);
    S->dc = S->dw = S->dh = 0; /*-Wall*/
  }
  if (S->dw < 0 || S->dh < 0) pari_err_TYPE("lfuninit [domain]", dom);
}

/* do we have dom \subset dom0 ? dom = [center, width, height] */
int
sdomain_isincl(double k, GEN dom, GEN dom0)
{
  struct lfunp S0, S;
  parse_dom(k, dom, &S);
  parse_dom(k, dom0, &S0);
  return S0.dc - S0.dw <= S.dc - S.dw
      && S0.dc + S0.dw >= S.dc + S.dw && S0.dh >= S.dh;
}

static int
checklfuninit(GEN linit, GEN dom, long der, long bitprec)
{
  GEN ldata = linit_get_ldata(linit);
  GEN domain = lfun_get_domain(linit_get_tech(linit));
  return domain_get_der(domain) >= der
    && domain_get_bitprec(domain) >= bitprec
    && sdomain_isincl(gtodouble(ldata_get_k(ldata)), dom, domain_get_dom(domain));
}

static GEN
ginvsqrtvec(GEN x, long prec)
{
  if (is_vec_t(typ(x)))
    pari_APPLY_same(ginv(gsqrt(gel(x,i), prec)))
  else return ginv(gsqrt(x, prec));
}

GEN
lfuninit_make(long t, GEN ldata, GEN tech, GEN domain)
{
  GEN Vga = ldata_get_gammavec(ldata);
  long d = lg(Vga)-1;
  GEN w2 = gen_1, k2 = gmul2n(ldata_get_k(ldata), -1);
  GEN expot = gdivgu(gadd(gmulsg(d, gsubgs(k2, 1)), sumVga(Vga)), 4);
  if (typ(ldata_get_dual(ldata))==t_INT)
  {
    GEN eno = ldata_get_rootno(ldata);
    long prec = nbits2prec( domain_get_bitprec(domain) );
    if (!isint1(eno)) w2 = ginvsqrtvec(eno, prec);
  }
  tech = mkvec3(domain, tech, mkvec4(k2, w2, expot, gammafactor(Vga)));
  return mkvec3(mkvecsmall(t), ldata, tech);
}

static void
lfunparams2(struct lfunp *S)
{
  GEN L = S->L, an = S->an, bn = S->bn;
  double pmax;
  long m, nan, nmax, neval, M = S->M;

  S->vgaell = 0;
  /* try to reduce parameters now we know the a_n (some may be 0) */
  if (typ(an) == t_VEC) an = RgV_kill0(an);
  nan = S->nmax; /* lg(an)-1 may be large than this */
  nmax = neval = 0;
  if (!bn)
    for (m = 0; m <= M; m++)
    {
      long n = minss(nan, L[m+1]);
      while (n > 0 && !gel(an,n)) n--;
      if (n > nmax) nmax = n;
      neval += n;
      L[m+1] = n; /* reduce S->L[m+1] */
    }
  else
  {
    if (typ(bn) == t_VEC) bn = RgV_kill0(bn);
    for (m = 0; m <= M; m++)
    {
      long n = minss(nan, L[m+1]);
      while (n > 0 && !gel(an,n) && !gel(bn,n)) n--;
      if (n > nmax) nmax = n;
      neval += n;
      L[m+1] = n; /* reduce S->L[m+1] */
    }
  }
  if (DEBUGLEVEL >= 1) err_printf("expected evaluations: %ld\n", neval);
  for (; M > 0; M--)
    if (L[M+1]) break;
  setlg(L, M+2);
  S->M = M;
  S->nmax = nmax;

  /* need K(n*exp(mh)/sqrt(N)) to absolute accuracy
   *   D + k1*log(n) + max(m * sig0 - sub2, 0) */
  pmax = S->D + S->k1 * log2(L[1]);
  if (S->MAXs)
  {
    double sig0 = S->MAXs/S->m0, sub2 = S->sub / M_LN2;
    for (m = ceil(sub2 / sig0); m <= S->M; m++)
    {
      double c = S->D + m*sig0 - sub2;
      if (S->k1 > 0) c += S->k1 * log2(L[m+1]);
      pmax = maxdd(pmax, c);
    }
  }
  S->Dmax = pmax;
  S->precmax = nbits2prec(pmax);
}

static GEN
lfun_init_theta(GEN ldata, GEN eno, struct lfunp *S)
{
  GEN an2, dual, tdom = NULL, Vga = ldata_get_gammavec(ldata);
  long L, prec = S->precmax;
  if (eno)
    L = S->nmax;
  else
  {
    tdom = dbltor(sqrt(0.5));
    L = maxss(S->nmax, lfunthetacost(ldata, tdom, 0, S->D));
  }
  dual = ldata_get_dual(ldata);
  S->an = ldata_vecan(ldata_get_an(ldata), L, prec);
  S->bn = typ(dual)==t_INT? NULL: ldata_vecan(dual, S->nmax, prec);
  if (!vgaell(Vga)) lfunparams2(S);
  else
  {
    S->an = antwist(S->an, Vga, prec);
    if (S->bn) S->bn = antwist(S->bn, Vga, prec);
    S->vgaell = 1;
  }
  an2 = lg(Vga)-1 == 1? antwist(S->an, Vga, prec): S->an;
  return lfunthetainit0(ldata, tdom, an2, 0, S->Dmax, 0);
}

GEN
lfuncost(GEN L, GEN dom, long der, long bit)
{
  pari_sp av = avma;
  GEN ldata = lfunmisc_to_ldata_shallow(L);
  GEN w, k = ldata_get_k(ldata);
  struct lfunp S;

  parse_dom(gtodouble(k), dom, &S);
  lfunp_set(ldata, der, bit, &S);
  w = ldata_get_rootno(ldata);
  if (isintzero(w)) /* for lfunrootres */
    S.nmax = maxss(S.nmax, lfunthetacost(ldata, dbltor(sqrt(0.5)), 0, bit+1));
  set_avma(av); return mkvecsmall2(S.nmax, S.Dmax);
}
GEN
lfuncost0(GEN L, GEN dom, long der, long bitprec)
{
  pari_sp av = avma;
  GEN C;

  if (is_linit(L))
  {
    GEN tech = linit_get_tech(L);
    GEN domain = lfun_get_domain(tech);
    dom = domain_get_dom(domain);
    der = domain_get_der(domain);
    bitprec = domain_get_bitprec(domain);
    if (linit_get_type(L) == t_LDESC_PRODUCT)
    {
      GEN v = lfunprod_get_fact(linit_get_tech(L)), F = gel(v,1);
      long i, l = lg(F);
      C = cgetg(l, t_VEC);
      for (i = 1; i < l; ++i)
        gel(C, i) = zv_to_ZV( lfuncost(gel(F,i), dom, der, bitprec) );
      return gerepileupto(av, C);
    }
  }
  if (!dom) pari_err_TYPE("lfuncost [missing s domain]", L);
  C = lfuncost(L,dom,der,bitprec);
  return gerepileupto(av, zv_to_ZV(C));
}

GEN
lfuninit(GEN lmisc, GEN dom, long der, long bitprec)
{
  pari_sp av = avma;
  GEN poqk, AB, R, h, theta, ldata, eno, r, domain, tech, k;
  struct lfunp S;

  if (is_linit(lmisc))
  {
    long t = linit_get_type(lmisc);
    if (t==t_LDESC_INIT || t==t_LDESC_PRODUCT)
    {
      if (checklfuninit(lmisc, dom, der, bitprec)) return lmisc;
      pari_warn(warner,"lfuninit: insufficient initialization");
    }
  }
  ldata = lfunmisc_to_ldata_shallow(lmisc);

  switch (ldata_get_type(ldata))
  {
  case t_LFUN_NF:
    {
      GEN T = gel(ldata_get_an(ldata), 2);
      return gerepilecopy(av, lfunzetakinit(T, dom, der, bitprec));
    }
  case t_LFUN_ABELREL:
    {
      GEN T = gel(ldata_get_an(ldata), 2);
      return gerepilecopy(av, lfunabelianrelinit(gel(T,1), gel(T,2), dom, der, bitprec));
    }
  }
  k = ldata_get_k(ldata);
  parse_dom(gtodouble(k), dom, &S);
  lfunp_set(ldata, der, bitprec, &S);
  ldata = ldata_newprec(ldata, nbits2prec(S.Dmax));
  r = ldata_get_residue(ldata);
  /* Note: all guesses should already have been performed (thetainit more
   * expensive than needed: should be either tdom = 1 or bitprec = S.D).
   * BUT if the root number / polar part do not have an algebraic
   * expression, there is no way to do this until we know the
   * precision, i.e. now. So we can't remove guessing code from here and
   * lfun_init_theta */
  if (r && isintzero(r)) eno = NULL;
  else
  {
    eno = ldata_get_rootno(ldata);
    if (isintzero(eno)) eno = NULL;
  }
  theta = lfun_init_theta(ldata, eno, &S);
  if (eno && !r)
    R = gen_0;
  else
  {
    GEN v = lfunrootres(theta, S.D);
    ldata = shallowcopy(ldata);
    gel(ldata, 6) = gel(v,3);
    r = gel(v,1);
    R = gel(v,2);
    if (isintzero(r)) setlg(ldata,7); else gel(ldata, 7) = r;
  }
  h = divru(mplog2(S.precmax), S.m0);
  /* exp(kh/2 . [0..M]) */
  poqk = gequal0(k) ? NULL
       : gpowers(gprec_w(mpexp(gmul2n(gmul(k,h), -1)), S.precmax), S.M);
  AB = lfuninit_ab(theta, h, &S);
  if (S.bn)
  {
    GEN A = gel(AB,1), B = gel(AB,2);
    A = lfuninit_pol(A, poqk, S.precmax);
    B = lfuninit_pol(B, poqk, S.precmax);
    AB = mkvec2(A, B);
  }
  else
    AB = lfuninit_pol(AB, poqk, S.precmax);
  tech = mkvec3(h, AB, R);
  domain = mkvec2(dom, mkvecsmall2(der, bitprec));
  return gerepilecopy(av, lfuninit_make(t_LDESC_INIT, ldata, tech, domain));
}

GEN
lfuninit0(GEN lmisc, GEN dom, long der, long bitprec)
{
  GEN z = lfuninit(lmisc, dom, der, bitprec);
  return z == lmisc? gcopy(z): z;
}

/* If s is a pole of Lambda, return polar part at s; else return NULL */
static GEN
lfunpoleresidue(GEN R, GEN s)
{
  long j;
  for (j = 1; j < lg(R); j++)
  {
    GEN Rj = gel(R, j), be = gel(Rj, 1);
    if (gequal(s, be)) return gel(Rj, 2);
  }
  return NULL;
}

/* Compute contribution of polar part at s when not a pole. */
static GEN
veccothderivn(GEN a, long n)
{
  long i;
  pari_sp av = avma;
  GEN c = pol_x(0), cp = mkpoln(3, gen_m1, gen_0, gen_1);
  GEN v = cgetg(n+2, t_VEC);
  gel(v, 1) = poleval(c, a);
  for(i = 2; i <= n+1; i++)
  {
    c = ZX_mul(ZX_deriv(c), cp);
    gel(v, i) = gdiv(poleval(c, a), mpfact(i-1));
  }
  return gerepilecopy(av, v);
}

static GEN
polepart(long n, GEN h, GEN C)
{
  GEN h2n = gpowgs(gdiv(h, gen_2), n-1);
  GEN res = gmul(h2n, gel(C,n));
  return odd(n)? res : gneg(res);
}

static GEN
lfunsumcoth(GEN R, GEN s, GEN h, long prec)
{
  long i,j;
  GEN S = gen_0;
  for (j = 1; j < lg(R); ++j)
  {
    GEN r = gel(R,j), be = gel(r,1), Rj = gel(r, 2);
    long e = valser(Rj);
    GEN z1 = gexpm1(gmul(h, gsub(s,be)), prec); /* exp(h(s-beta))-1 */
    GEN c1 = gaddgs(gdivsg(2, z1), 1); /* coth((h/2)(s-beta)) */
    GEN C1 = veccothderivn(c1, 1-e);
    for (i = e; i < 0; i++)
    {
      GEN Rbe = mysercoeff(Rj, i);
      GEN p1 = polepart(-i, h, C1);
      S = gadd(S, gmul(Rbe, p1));
    }
  }
  return gmul2n(S, -1);
}

static GEN lfunlambda_OK(GEN linit, GEN s, GEN sdom, long bitprec);
/* L is a t_LDESC_PRODUCT Linit */
static GEN
lfunlambda_product(GEN L, GEN s, GEN sdom, long bitprec)
{
  GEN ldata = linit_get_ldata(L), v = lfunprod_get_fact(linit_get_tech(L));
  GEN r = gen_1, F = gel(v,1), E = gel(v,2), C = gel(v,3), cs = conj_i(s);
  long i, l = lg(F), isreal = gequal(imag_i(s), imag_i(cs));
  for (i = 1; i < l; ++i)
  {
    GEN f = lfunlambda_OK(gel(F, i), s, sdom, bitprec);
    if( DEBUGLEVEL>=2) err_printf("lfunlambda(%ld): %Ps\n",i,f);
    if (typ(f)==t_VEC) f = RgV_prod(f);
    if (E[i]) r = gmul(r, gpowgs(f, E[i]));
    if (C[i])
    {
      GEN fc = isreal? f: conj_i(lfunlambda_OK(gel(F, i), cs, sdom, bitprec));
      r = gmul(r, gpowgs(fc, C[i]));
    }
  }
  return (ldata_isreal(ldata) && gequal0(imag_i(s)))? real_i(r): r;
}

/* s a t_SER */
static long
der_level(GEN s)
{ return signe(s)? lg(s)-3: valser(s)-1; }

/* s a t_SER; return coeff(s, X^0) */
static GEN
ser_coeff0(GEN s) { return simplify_shallow(polcoef_i(s, 0, -1)); }

static GEN
get_domain(GEN s, GEN *dom, long *der)
{
  GEN sa = s;
  *der = 0;
  switch(typ(s))
  {
    case t_POL:
    case t_RFRAC: s = toser_i(s);
    case t_SER:
      *der = der_level(s);
      sa = ser_coeff0(s);
  }
  *dom = mkvec3(real_i(sa), gen_0, gabs(imag_i(sa),DEFAULTPREC));
  return s;
}
/* assume s went through get_domain and s/bitprec belong to domain */
static GEN
lfunlambda_OK(GEN linit, GEN s, GEN sdom, long bitprec)
{
  GEN eno, ldata, tech, h, pol;
  GEN S, S0 = NULL, k2, cost;
  long prec, prec0;
  struct lfunp D, D0;

  if (linit_get_type(linit) == t_LDESC_PRODUCT)
    return lfunlambda_product(linit, s, sdom, bitprec);
  ldata = linit_get_ldata(linit);
  eno = ldata_get_rootno(ldata);
  tech = linit_get_tech(linit);
  h = lfun_get_step(tech); prec = realprec(h);
  /* try to reduce accuracy */
  parse_dom(0, sdom, &D0);
  parse_dom(0, domain_get_dom(lfun_get_domain(tech)), &D);
  if (0.8 * D.dh > D0.dh)
  {
    cost = lfuncost(linit, sdom, typ(s)==t_SER? der_level(s): 0, bitprec);
    prec0 = nbits2prec(cost[2]);
    if (prec0 < prec) { prec = prec0; h = gprec_w(h, prec); }
  }
  pol = lfun_get_pol(tech);
  s = gprec_w(s, prec);
  if (ldata_get_residue(ldata))
  {
    GEN R = lfun_get_Residue(tech);
    GEN Ra = lfunpoleresidue(R, s);
    if (Ra) return gprec_w(Ra, nbits2prec(bitprec));
    S0 = lfunsumcoth(R, s, h, prec);
  }
  k2 = lfun_get_k2(tech);
  if (typ(pol)==t_POL && typ(s) != t_SER && gequal(real_i(s), k2))
  { /* on critical line: shortcut */
    GEN polz, b = imag_i(s);
    polz = gequal0(b)? poleval(pol,gen_1): poleval(pol, expIr(gmul(h,b)));
    S = gadd(polz, gmulvec(eno, conj_i(polz)));
  }
  else
  {
    GEN z = gexp(gmul(h, gsub(s, k2)), prec);
    GEN zi = ginv(z), zc = conj_i(zi);
    if (typ(pol)==t_POL)
      S = gadd(poleval(pol, z), gmulvec(eno, conj_i(poleval(pol, zc))));
    else
      S = gadd(poleval(gel(pol,1), z), gmulvec(eno, poleval(gel(pol,2), zi)));
  }
  if (S0) S = gadd(S,S0);
  return gprec_w(gmul(S,h), nbits2prec(bitprec));
}
GEN
lfunlambda(GEN lmisc, GEN s, long bitprec)
{
  pari_sp av = avma;
  GEN linit, dom, z;
  long der;
  s = get_domain(s, &dom, &der);
  linit = lfuninit(lmisc, dom, der, bitprec);
  z = lfunlambda_OK(linit,s, dom, bitprec);
  return gerepilecopy(av, z);
}

static long
is_ser(GEN x)
{
  long t = typ(x);
  if (t == t_SER) return 1;
  if (!is_vec_t(t) || lg(x)==1) return 0;
  if (typ(gel(x,1))==t_SER) return 1;
  return 0;
}

static GEN
lfunser(GEN L)
{
  long v = valser(L);
  if (v > 0) return gen_0;
  if (v == 0) L = gel(L, 2);
  else
    setlg(L, minss(lg(L), 2-v));
  return L;
}

static GEN
lfunservec(GEN x)
{
  if (typ(x)==t_SER) return lfunser(x);
  pari_APPLY_same(lfunser(gel(x,i)))
}
static GEN
lfununext(GEN L)
{
  setlg(L, maxss(lg(L)-1, valser(L)? 2: 3));
  return normalizeser(L);
}
static GEN
lfununextvec(GEN x)
{
  if (typ(x)==t_SER) return lfununext(x);
  pari_APPLY_same(lfununext(gel(x,i)));
}

/* assume lmisc is an linit, s went through get_domain and s/bitprec belong
 * to domain */
static GEN
lfun_OK(GEN linit, GEN s, GEN sdom, long bitprec)
{
  GEN N, gas, S, FVga, res, ss = s;
  long prec = nbits2prec(bitprec), ext;

  FVga = lfun_get_factgammavec(linit_get_tech(linit));
  S = lfunlambda_OK(linit, s, sdom, bitprec);
  if (is_ser(S))
  {
    GEN r = typ(S)==t_SER ? S : gel(S,1);
    long d = lg(r) - 2 + fracgammadegree(FVga);
    if (typ(s) == t_SER)
      ss = sertoser(s, d);
    else
      ss = deg1ser_shallow(gen_1, s, varn(r), d);
  }
  gas = gammafactproduct(FVga, ss, &ext, prec);
  N = ldata_get_conductor(linit_get_ldata(linit));
  res = gdiv(S, gmul(gpow(N, gdivgu(ss, 2), prec), gas));
  if (typ(s) != t_SER && is_ser(res)) res = lfunservec(res);
  else if (ext) res = lfununextvec(res);
  return gprec_w(res, prec);
}

GEN
lfun(GEN lmisc, GEN s, long bitprec)
{
  pari_sp av = avma;
  GEN linit, dom, z;
  long der;
  s = get_domain(s, &dom, &der);
  if (!der && typ(s) == t_INT && !is_bigint(s))
  { /* special value ? */
    GEN ldata;
    long t, ss = itos(s);
    if (is_linit(lmisc))
      ldata = linit_get_ldata(lmisc);
    else
      lmisc = ldata = lfunmisc_to_ldata_shallow(lmisc);
    t = ldata_get_type(ldata);
    if (t == t_LFUN_KRONECKER || t == t_LFUN_ZETA)
    {
      long D = itos_or_0(gel(ldata_get_an(ldata), 2));
      if (D)
      {
        if (ss <= 0) return lfunquadneg(D, ss);
        /* ss > 0 */
        if ((!odd(ss) && D > 0) || (odd(ss) && D < 0))
        {
          long prec = nbits2prec(bitprec), q = labs(D);
          ss = 1 - ss; /* <= 0 */
          z = powrs(divrs(mppi(prec + EXTRAPREC64), q), 1-ss);
          z = mulrr(shiftr(z, -ss), sqrtr_abs(utor(q, prec)));
          z = gdiv(z, mpfactr(-ss, prec));
          if (smodss(ss, 4) > 1) togglesign(z);
          return gmul(z, lfunquadneg(D, ss));
        }
      }
    }
  }
  linit = lfuninit(lmisc, dom, der, bitprec);
  z = lfun_OK(linit, s, dom, bitprec);
  return gerepilecopy(av, z);
}

/* given a t_SER a+x*s(x), return x*s(x), shallow */
static GEN
sersplit1(GEN s, GEN *head)
{
  long i, l = lg(s);
  GEN y;
  *head = simplify_shallow(mysercoeff(s, 0));
  if (valser(s) > 0) return s;
  y = cgetg(l-1, t_SER); y[1] = s[1];
  setvalser(y, 1);
  for (i=3; i < l; i++) gel(y,i-1) = gel(s,i);
  return normalizeser(y);
}

/* order of pole of Lambda at s (0 if regular point) */
static long
lfunlambdaord(GEN linit, GEN s)
{
  GEN tech = linit_get_tech(linit);
  if (linit_get_type(linit)==t_LDESC_PRODUCT)
  {
    GEN v = lfunprod_get_fact(linit_get_tech(linit));
    GEN F = gel(v, 1), E = gel(v, 2), C = gel(v, 3);
    long i, ex = 0, l = lg(F);
    for (i = 1; i < l; i++)
      ex += lfunlambdaord(gel(F,i), s) * (E[i]+C[i]);
    return ex;
  }
  if (ldata_get_residue(linit_get_ldata(linit)))
  {
    GEN r = lfunpoleresidue(lfun_get_Residue(tech), s);
    if (r) return lg(r)-2;
  }
  return 0;
}

static GEN
derser(GEN res, long m)
{
  long v = valser(res);
  if (v > m) return gen_0;
  if (v >= 0)
    return gmul(mysercoeff(res, m), mpfact(m));
  else
    return derivn(res, m, -1);
}

static GEN
derservec(GEN x, long m) { pari_APPLY_same(derser(gel(x,i),m)) }

/* derivative of order m > 0 of L (flag = 0) or Lambda (flag = 1) */
static GEN
lfunderiv(GEN lmisc, long m, GEN s, long flag, long bitprec)
{
  pari_sp ltop = avma;
  GEN res, S = NULL, linit, dom;
  long der, prec = nbits2prec(bitprec);
  if (m <= 0) pari_err_DOMAIN("lfun", "D", "<=", gen_0, stoi(m));
  s = get_domain(s, &dom, &der);
  linit = lfuninit(lmisc, dom, der + m, bitprec);
  if (typ(s) == t_SER)
  {
    long v, l = lg(s)-1;
    GEN sh;
    if (valser(s) < 0) pari_err_DOMAIN("lfun","valuation", "<", gen_0, s);
    S = sersplit1(s, &sh);
    v = valser(S);
    s = deg1ser_shallow(gen_1, sh, varn(S), m + (l+v-1)/v);
  }
  else
  {
    long ex = lfunlambdaord(linit, s);
    /* HACK: pretend lfuninit was done to right accuracy */
    if (gequal0(s)) s = gen_0;
    s = deg1ser_shallow(gen_1, s, 0, m+1+ex);
  }
  res = flag ? lfunlambda_OK(linit, s, dom, bitprec):
               lfun_OK(linit, s, dom, bitprec);
  if (S)
    res = gsubst(derivn(res, m, -1), varn(S), S);
  else if (typ(res)==t_SER)
  {
    long v = valser(res);
    if (v > m) { set_avma(ltop); return gen_0; }
    if (v >= 0)
      res = gmul(mysercoeff(res, m), mpfact(m));
    else
      res = derivn(res, m, -1);
  }
  else if (is_ser(res))
    res = derservec(res, m);
  return gerepilecopy(ltop, gprec_w(res, prec));
}

GEN
lfunlambda0(GEN lmisc, GEN s, long der, long bitprec)
{
  return der? lfunderiv(lmisc, der, s, 1, bitprec)
            : lfunlambda(lmisc, s, bitprec);
}

GEN
lfun0(GEN lmisc, GEN s, long der, long bitprec)
{
  return der? lfunderiv(lmisc, der, s, 0, bitprec)
            : lfun(lmisc, s, bitprec);
}

GEN
lfunhardy(GEN lmisc, GEN t, long bitprec)
{
  pari_sp ltop = avma;
  long prec = nbits2prec(bitprec), d;
  GEN argz, z, linit, ldata, tech, dom, w2, k2, E, h, a, k;

  switch(typ(t))
  {
    case t_INT: case t_FRAC: case t_REAL: break;
    default: pari_err_TYPE("lfunhardy",t);
  }

  ldata = lfunmisc_to_ldata_shallow(lmisc);
  if (!is_linit(lmisc)) lmisc = ldata;
  k = ldata_get_k(ldata);
  d = ldata_get_degree(ldata);
  dom = mkvec3(gmul2n(k, -1), gen_0, gabs(t,LOWDEFAULTPREC));
  linit = lfuninit(lmisc, dom, 0, bitprec);
  tech = linit_get_tech(linit);
  w2 = lfun_get_w2(tech);
  k2 = lfun_get_k2(tech);
  E = lfun_get_expot(tech); /* 4E = d(k2 - 1) + real(vecsum(Vga)) */
  z = mkcomplex(k2, t);
  /* more accurate than garg: k/2 in Q */
  argz = gequal0(k2)? Pi2n(-1, prec): gatan(gdiv(t, k2), prec);
  prec = precision(argz);
  /* prec may have increased: don't lose accuracy if |z|^2 is exact */
  a = gsub(gmulsg(d, gmul(t, gmul2n(argz,-1))),
           gmul(E, glog(gnorm(z),prec)));
  h = lfunlambda_OK(linit, z, dom, bitprec);
  if (!isint1(w2) && typ(ldata_get_dual(ldata))==t_INT)
    h = mulrealvec(h, w2);
  if (typ(h) == t_COMPLEX && gexpo(imag_i(h)) < -(bitprec >> 1))
    h = real_i(h);
  return gerepileupto(ltop, gmul(h, gexp(a, prec)));
}

/* L = log(t); return  \sum_{i = 0}^{v-1}  R[-i-1] L^i/i! */
static GEN
theta_pole_contrib(GEN R, long v, GEN L)
{
  GEN s = mysercoeff(R,-v);
  long i;
  for (i = v-1; i >= 1; i--)
    s = gadd(mysercoeff(R,-i), gdivgu(gmul(s,L), i));
  return s;
}
/* subtract successively rather than adding everything then subtracting.
 * The polar part is "large" and suffers from cancellation: a little stabler
 * this way */
static GEN
theta_add_polar_part(GEN S, GEN R, GEN t, long prec)
{
  GEN logt = NULL;
  long j, l = lg(R);
  for (j = 1; j < l; j++)
  {
    GEN Rj = gel(R,j), b = gel(Rj,1), Rb = gel(Rj,2);
    long v = -valser(Rb);
    if (v > 1 && !logt) logt = glog(t, prec);
    S = gsub(S, gmul(theta_pole_contrib(Rb,v,logt), gpow(t,b,prec)));
  }
  return S;
}

static long
lfuncheckfeq_i(GEN theta, GEN thetad, GEN t0, GEN t0i, long bitprec)
{
  GEN ldata = linit_get_ldata(theta);
  GEN S0, S0i, w, eno;
  long prec = nbits2prec(bitprec);
  if (thetad)
    S0 = lfuntheta(thetad, t0, 0, bitprec);
  else
    S0 = conj_i(lfuntheta(theta, conj_i(t0), 0, bitprec));
  S0i = lfuntheta(theta, t0i, 0, bitprec);

  eno = ldata_get_rootno(ldata);
  if (ldata_get_residue(ldata))
  {
    GEN R = theta_get_R(linit_get_tech(theta));
    if (gequal0(R))
    {
      GEN v, r;
      long t = ldata_get_type(ldata);
      if (t == t_LFUN_NF || t == t_LFUN_ABELREL)
      { /* inefficient since theta not needed; no need to optimize for this
           (artificial) query [e.g. lfuncheckfeq(t_POL)] */
        GEN L = lfuninit(ldata,zerovec(3),0,bitprec);
        return lfuncheckfeq(L,t0,bitprec);
      }
      v = lfunrootres(theta, bitprec);
      r = gel(v,1);
      if (gequal0(eno)) eno = gel(v,3);
      R = lfunrtoR_i(ldata, r, eno, nbits2prec(bitprec));
    }
    S0i = theta_add_polar_part(S0i, R, t0, prec);
  }
  if (gequal0(S0i) || gequal0(S0)) pari_err_PREC("lfuncheckfeq");

  w = gdivvec(S0i, gmul(S0, gpow(t0, ldata_get_k(ldata), prec)));
  /* missing rootno: guess it */
  if (gequal0(eno)) eno = lfunrootno(theta, bitprec);
  w = gsubvec(w, eno);
  if (thetad) w = gdivvec(w, eno); /* |eno| may be large in non-dual case */
  return gexpo(w);
}

/* Check whether the coefficients, conductor, weight, polar part and root
 * number are compatible with the functional equation at t0 and 1/t0.
 * Different from lfunrootres. */
long
lfuncheckfeq(GEN lmisc, GEN t0, long bitprec)
{
  GEN ldata, theta, thetad, t0i;
  pari_sp av;

  if (is_linit(lmisc) && linit_get_type(lmisc)==t_LDESC_PRODUCT)
  {
    GEN v = lfunprod_get_fact(linit_get_tech(lmisc)), F = gel(v,1);
    long i, b = -bitprec, l = lg(F);
    for (i = 1; i < l; i++) b = maxss(b, lfuncheckfeq(gel(F,i), t0, bitprec));
    return b;
  }
  av = avma;
  if (!t0)
  { /* ~Pi/3 + I/7, some random complex number */
    t0 = mkcomplex(uutoQ(355,339), uutoQ(1,7));
    t0i = ginv(t0);
  }
  else if (gcmpgs(gnorm(t0), 1) < 0) { t0i = t0; t0 = ginv(t0); }
  else t0i = ginv(t0);
  /* |t0| >= 1 */
  theta = lfunthetacheckinit(lmisc, t0i, 0, bitprec);
  ldata = linit_get_ldata(theta);
  thetad = theta_dual(theta, ldata_get_dual(ldata));
  return gc_long(av, lfuncheckfeq_i(theta, thetad, t0, t0i, bitprec));
}

/*******************************************************************/
/*       Compute root number and residues                          */
/*******************************************************************/
/* round root number to \pm 1 if close to integer. */
static GEN
ropm1(GEN w, long prec)
{
  long e;
  GEN r;
  if (typ(w) == t_INT) return w;
  r = grndtoi(w, &e);
  return (e < -prec2nbits(prec)/2)? r: w;
}

/* theta for t=1/sqrt(2) and t2==2t simultaneously, saving 25% of the work.
 * Assume correct initialization (no thetacheck) */
static void
lfunthetaspec(GEN linit, long bitprec, GEN *pv, GEN *pv2)
{
  pari_sp av = avma, av2;
  GEN t, Vga, an, K, ldata, thetainit, v, v2, vroots;
  long L, prec, n, d;

  ldata = linit_get_ldata(linit);
  thetainit = linit_get_tech(linit);
  prec = nbits2prec(bitprec);
  Vga = ldata_get_gammavec(ldata); d = lg(Vga)-1;
  if (Vgaeasytheta(Vga))
  {
    GEN v2 = sqrtr(real2n(1, nbits2prec(bitprec)));
    GEN v = shiftr(v2,-1);
    *pv = lfuntheta(linit, v,  0, bitprec);
    *pv2= lfuntheta(linit, v2, 0, bitprec);
    return;
  }
  an = RgV_kill0( theta_get_an(thetainit) );
  L = lg(an)-1;
  /* to compute theta(1/sqrt(2)) */
  t = ginv(gsqrt(gmul2n(ldata_get_conductor(ldata), 1), prec));
  /* t = 1/sqrt(2N) */

  /* From then on, the code is generic and could be used to compute
   * theta(t) / theta(2t) without assuming t = 1/sqrt(2) */
  K = theta_get_K(thetainit);
  vroots = mkvroots(d, L, prec);
  t = gpow(t, gdivgu(gen_2, d), prec); /* rt variant: t->t^(2/d) */
  /* v = \sum_{n <= L, n odd} a_n K(nt) */
  for (v = gen_0, n = 1; n <= L; n+=2)
  {
    GEN tn, Kn, a = gel(an, n);

    if (!a) continue;
    av2 = avma;
    tn = gmul(t, gel(vroots,n));
    Kn = gammamellininvrt(K, tn, bitprec);
    v = gerepileupto(av2, gadd(v, gmul(a,Kn)));
  }
  /* v += \sum_{n <= L, n even} a_n K(nt), v2 = \sum_{n <= L/2} a_n K(2n t) */
  for (v2 = gen_0, n = 1; n <= L/2; n++)
  {
    GEN t2n, K2n, a = gel(an, n), a2 = gel(an,2*n);

    if (!a && !a2) continue;
    av2 = avma;
    t2n = gmul(t, gel(vroots,2*n));
    K2n = gerepileupto(av2, gammamellininvrt(K, t2n, bitprec));
    if (a) v2 = gadd(v2, gmul(a, K2n));
    if (a2) v = gadd(v,  gmul(a2,K2n));
  }
  *pv = v;
  *pv2 = v2;
  gerepileall(av, 2, pv,pv2);
}

static GEN
Rtor(GEN a, GEN R, GEN ldata, long prec)
{
  GEN FVga = gammafactor(ldata_get_gammavec(ldata));
  GEN Na = gpow(ldata_get_conductor(ldata), gdivgu(a,2), prec);
  long ext;
  return gdiv(R, gmul(Na, gammafactproduct(FVga, a, &ext, prec)));
}

/* v = theta~(t), vi = theta(1/t) */
static GEN
get_eno(GEN R, GEN k, GEN t, GEN v, GEN vi, long vx, long bitprec, long force)
{
  long prec = nbits2prec(bitprec);
  GEN a0, a1, S = deg1pol(gmul(gpow(t,k,prec), gneg(v)), vi, vx);

  S = theta_add_polar_part(S, R, t, prec);
  if (typ(S) != t_POL || degpol(S) != 1) return NULL;
  a1 = gel(S,3); if (!force && gexpo(a1) < -bitprec/4) return NULL;
  a0 = gel(S,2);
  return gdivvec(a0, gneg(a1));

}
/* Return w using theta(1/t) - w t^k \bar{theta}(t) = polar_part(t,w).
 * The full Taylor development of L must be known */
GEN
lfunrootno(GEN linit, long bitprec)
{
  GEN ldata, t, eno, v, vi, R, thetad;
  long c = 0, prec = nbits2prec(bitprec), vx = fetch_var();
  GEN k;
  pari_sp av;

  /* initialize for t > 1/sqrt(2) */
  linit = lfunthetacheckinit(linit, dbltor(sqrt(0.5)), 0, bitprec);
  ldata = linit_get_ldata(linit);
  k = ldata_get_k(ldata);
  R = ldata_get_residue(ldata)? lfunrtoR_eno(ldata, pol_x(vx), prec)
                              : cgetg(1, t_VEC);
  t = gen_1;
  v = lfuntheta(linit, t, 0, bitprec);
  thetad = theta_dual(linit, ldata_get_dual(ldata));
  vi = !thetad ? conj_i(v): lfuntheta(thetad, t, 0, bitprec);
  eno = get_eno(R,k,t,vi,v, vx, bitprec, 0);
  if (!eno && !thetad)
  { /* t = sqrt(2), vi = theta(1/t), v = theta(t) */
    lfunthetaspec(linit, bitprec, &vi, &v);
    t = sqrtr(utor(2, prec));
    eno = get_eno(R,k,t,conj_i(v),vi, vx, bitprec, 0);
  }
  av = avma;
  while (!eno)
  {
    t = addsr(1, shiftr(utor(pari_rand(), prec), -2-BITS_IN_LONG));
    /* t in [1,1.25[ */
    v = thetad? lfuntheta(thetad, t, 0, bitprec)
              : conj_i(lfuntheta(linit, t, 0, bitprec));
    vi = lfuntheta(linit, ginv(t), 0, bitprec);
    eno = get_eno(R,k,t,v,vi, vx, bitprec, c++ == 5);
    set_avma(av);
  }
  delete_var(); return ropm1(eno,prec);
}

/* Find root number and/or residues when L-function coefficients and
   conductor are known. For the moment at most a single residue allowed. */
GEN
lfunrootres(GEN data, long bitprec)
{
  pari_sp ltop = avma;
  GEN k, w, r, R, a, b, e, v, v2, be, ldata, linit;
  long prec;

  ldata = lfunmisc_to_ldata_shallow(data);
  r = ldata_get_residue(ldata);
  k = ldata_get_k(ldata);
  w = ldata_get_rootno(ldata);
  if (r) r = normalize_simple_pole(r, k);
  if (!r || residues_known(r))
  {
    if (isintzero(w)) w = lfunrootno(data, bitprec);
    if (!r)
      r = R = gen_0;
    else
      R = lfunrtoR_eno(ldata, w, nbits2prec(bitprec));
    return gerepilecopy(ltop, mkvec3(r, R, w));
  }
  linit = lfunthetacheckinit(data, dbltor(sqrt(0.5)), 0, bitprec);
  prec = nbits2prec(bitprec);
  if (lg(r) > 2) pari_err_IMPL("multiple poles in lfunrootres");
  /* Now residue unknown, and r = [[be,0]]. */
  be = gmael(r, 1, 1);
  if (ldata_isreal(ldata) && gequalm1(w))
    R = lfuntheta(linit, gen_1, 0, bitprec);
  else
  {
    GEN p2k = gpow(gen_2,k,prec);
    lfunthetaspec(linit, bitprec, &v2, &v);
    if (gequal(gmulsg(2, be), k)) pari_err_IMPL("pole at k/2 in lfunrootres");
    if (gequal(be, k))
    {
      a = conj_i(gsub(gmul(p2k, v), v2));
      b = subiu(p2k, 1);
      e = gmul(gsqrt(p2k, prec), gsub(v2, v));
    }
    else
    {
      GEN tk2 = gsqrt(p2k, prec);
      GEN tbe = gpow(gen_2, be, prec);
      GEN tkbe = gpow(gen_2, gdivgu(gsub(k, be), 2), prec);
      a = conj_i(gsub(gmul(tbe, v), v2));
      b = gsub(gdiv(tbe, tkbe), tkbe);
      e = gsub(gmul(gdiv(tbe, tk2), v2), gmul(tk2, v));
    }
    if (isintzero(w))
    { /* Now residue unknown, r = [[be,0]], and w unknown. */
      GEN t0  = mkfrac(utoi(11),utoi(10));
      GEN th1 = lfuntheta(linit, t0,  0, bitprec);
      GEN th2 = lfuntheta(linit, ginv(t0), 0, bitprec);
      GEN tbe = gpow(t0, gmulsg(2, be), prec);
      GEN tkbe = gpow(t0, gsub(k, be), prec);
      GEN tk2 = gpow(t0, k, prec);
      GEN c = conj_i(gsub(gmul(tbe, th1), th2));
      GEN d = gsub(gdiv(tbe, tkbe), tkbe);
      GEN f = gsub(gmul(gdiv(tbe, tk2), th2), gmul(tk2, th1));
      GEN D = gsub(gmul(a, d), gmul(b, c));
      w = gdiv(gsub(gmul(d, e), gmul(b, f)), D);
    }
    w = ropm1(w, prec);
    R = gdiv(gsub(e, gmul(a, w)), b);
  }
  r = normalize_simple_pole(Rtor(be, R, ldata, prec), be);
  R = lfunrtoR_i(ldata, r, w, prec);
  return gerepilecopy(ltop, mkvec3(r, R, w));
}

/*******************************************************************/
/*                           Zeros                                 */
/*******************************************************************/
struct lhardyz_t {
  long bitprec, prec;
  GEN linit;
};

static GEN
lfunhardyzeros(void *E, GEN t)
{
  struct lhardyz_t *S = (struct lhardyz_t*)E;
  GEN z = gprec_wensure(lfunhardy(S->linit, t, S->bitprec), S->prec);
  return typ(z) == t_VEC ? RgV_prod(z): z;
}

/* initialize for computation on critical line up to height h, zero
 * of order <= m */
static GEN
lfuncenterinit(GEN lmisc, double h, long m, long bitprec)
{
  if (m < 0)
  { /* choose a sensible default */
    if (!is_linit(lmisc) || linit_get_type(lmisc) != t_LDESC_INIT) m = 4;
    else
    {
      GEN domain = lfun_get_domain(linit_get_tech(lmisc));
      m = domain_get_der(domain);
    }
  }
  return lfuninit(lmisc, mkvec(dbltor(h)), m, bitprec);
}

long
lfunorderzero(GEN lmisc, long m, long bitprec)
{
  pari_sp ltop = avma;
  GEN eno, ldata, linit, k2;
  long G, c0, c, st;

  if (is_linit(lmisc) && linit_get_type(lmisc) == t_LDESC_PRODUCT)
  {
    GEN M = gmael(linit_get_tech(lmisc), 2,1);
    long i, l = lg(M);
    for (c=0, i=1; i < l; i++) c += lfunorderzero(gel(M,i), m, bitprec);
    return c;
  }
  linit = lfuncenterinit(lmisc, 0, m, bitprec);
  ldata = linit_get_ldata(linit);
  eno = ldata_get_rootno(ldata);
  k2 = gmul2n(ldata_get_k(ldata), -1);
  G = -bitprec/2;
  c0 = 0; st = 1;
  if (typ(eno) == t_VEC)
  {
    long i, l = lg(eno), cnt = l-1, s = 0;
    GEN v = zero_zv(l-1);
    if (ldata_isreal(ldata)) st = 2;
    for (c = c0; cnt; c += st)
    {
      GEN L = lfun0(linit, k2, c, bitprec);
      for (i = 1; i < l; i++)
      {
        if (v[i]==0 && gexpo(gel(L,i)) > G)
        {
          v[i] = c; cnt--; s += c;
        }
      }
    }
    return gc_long(ltop,s);
  }
  else
  {
    if (ldata_isreal(ldata)) { st = 2; if (!gequal1(eno)) c0 = 1; }
    for (c = c0;; c += st)
      if (gexpo(lfun0(linit, k2, c, bitprec)) > G) return gc_long(ltop, c);
  }
}

/* assume T1 * T2 > 0, T1 <= T2 */
static void
lfunzeros_i(struct lhardyz_t *S, GEN *pw, long *ct, GEN T1, GEN T2, long d,
            GEN cN, GEN pi2, GEN pi2div, long precinit, long prec)
{
  GEN T = T1, w = *pw;
  long W = lg(w)-1, s = gsigne(lfunhardyzeros(S, T1));
  for(;;)
  {
    pari_sp av = avma;
    GEN D, T0, z;
    D = gcmp(T, pi2) < 0? cN
                        : gadd(cN, gmulsg(d, glog(gdiv(T, pi2), prec)));
    D = gdiv(pi2div, D);
    for(;;)
    {
      long s0;
      T0 = T; T = gadd(T, D);
      if (gcmp(T, T2) >= 0) T = T2;
      s0 = gsigne(lfunhardyzeros(S, T));
      if (s0 != s) { s = s0; break; }
      if (T == T2) { setlg(w, *ct); *pw = w; return; }
    }
    z = zbrent(S, lfunhardyzeros, T0, T, prec); /* T <= T2 */
    gerepileall(av, 2, &T, &z);
    if (*ct > W) { W *= 2; w = vec_lengthen(w, W); }
    if (typ(z) == t_REAL) z  = rtor(z, precinit);
    gel(w, (*ct)++) = z;
  }
  setlg(w, *ct); *pw = w;
}
GEN
lfunzeros(GEN ldata, GEN lim, long divz, long bitprec)
{
  pari_sp ltop = avma;
  GEN linit, pi2, pi2div, cN, w, T, h1, h2;
  long i, d, NEWD, c, ct, s1, s2, prec, prec0 = nbits2prec(bitprec);
  double maxt;
  struct lhardyz_t S;

  if (is_linit(ldata) && linit_get_type(ldata) == t_LDESC_PRODUCT)
  {
    GEN M = gmael(linit_get_tech(ldata), 2,1);
    long l = lg(M);
    w = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(w,i) = lfunzeros(gel(M,i), lim, divz, bitprec);
    return gerepileupto(ltop, vecsort0(shallowconcat1(w), NULL, 0));
  }
  if (typ(lim) == t_VEC)
  {
    if (lg(lim) != 3 || gcmp(gel(lim,1),gel(lim,2)) >= 0)
      pari_err_TYPE("lfunzeros",lim);
    h1 = gel(lim,1);
    h2 = gel(lim,2);
    maxt = maxdd(fabs(gtodouble(h1)), fabs(gtodouble(h2)));
  }
  else
  {
    if (gcmp(lim,gen_0) <= 0) pari_err_TYPE("lfunzeros",lim);
    h1 = gen_0;
    h2 = lim;
    maxt = gtodouble(h2);
  }
  S.linit = linit = lfuncenterinit(ldata, maxt, -1, bitprec);
  S.bitprec = bitprec;
  S.prec = prec0;
  ldata = linit_get_ldata(linit);
  d = ldata_get_degree(ldata);

  NEWD = minss((long) ceil(bitprec + (M_PI/(4*M_LN2)) * d * maxt),
               lfun_get_bitprec(linit_get_tech(linit)));
  prec = nbits2prec(NEWD);
  cN = gdiv(ldata_get_conductor(ldata), gpowgs(Pi2n(-1, prec), d));
  cN = gexpo(cN) >= 0? gaddsg(d, gmulsg(2, glog(cN, prec))): utoi(d);
  pi2 = Pi2n(1, prec);
  pi2div = gdivgu(pi2, labs(divz));
  s1 = gsigne(h1);
  s2 = gsigne(h2);
  w = cgetg(100+1, t_VEC); c = 1; ct = 0; T = NULL;
  if (s1 <= 0 && s2 >= 0)
  {
    GEN r = ldata_get_residue(ldata);
    if (!r || gequal0(r))
    {
      ct = lfunorderzero(linit, -1, bitprec);
      if (ct) T = real2n(-prec2nbits(prec) / (2*ct), prec);
    }
  }
  if (s1 <= 0)
  {
    if (s1 < 0)
      lfunzeros_i(&S, &w, &c, h1, T? negr(T): h2,
                  d, cN, pi2, pi2div, prec0, prec);
    if (ct)
    {
      long n = lg(w)-1;
      if (c + ct >= n) w = vec_lengthen(w, n + ct);
      for (i = 1; i <= ct; i++) gel(w,c++) = gen_0;
    }
  }
  if (s2 > 0 && (T || s1 >= 0))
    lfunzeros_i(&S, &w, &c, T? T: h1, h2, d, cN, pi2, pi2div, prec0, prec);
  return gerepilecopy(ltop, w);
}

/*******************************************************************/
/*       Guess conductor                                           */
/*******************************************************************/
struct huntcond_t {
  GEN k;
  GEN theta, thetad;
  GEN *pM, *psqrtM, *pMd, *psqrtMd;
};

static void
condset(struct huntcond_t *S, GEN M, long prec)
{
  *(S->pM) = M;
  *(S->psqrtM) = gsqrt(ginv(M), prec);
  if (S->thetad != S->theta)
  {
    *(S->pMd) = *(S->pM);
    *(S->psqrtMd) = *(S->psqrtM);
  }
}

/* M should eventually converge to N, the conductor. L has no pole. */
static GEN
wrap1(void *E, GEN M)
{
  struct huntcond_t *S = (struct huntcond_t*)E;
  GEN thetainit, tk, p1, p1inv;
  GEN t = mkfrac(stoi(11), stoi(10));
  long prec, bitprec;

  thetainit = linit_get_tech(S->theta);
  bitprec = theta_get_bitprec(thetainit);
  prec = nbits2prec(bitprec);
  condset(S, M, prec);
  tk = gpow(t, S->k, prec);
  p1 = lfuntheta(S->thetad, t, 0, bitprec);
  p1inv = lfuntheta(S->theta, ginv(t), 0, bitprec);
  return glog(gabs(gmul(tk, gdiv(p1, p1inv)), prec), prec);
}

/* M should eventually converge to N, the conductor. L has a pole. */
static GEN
wrap2(void *E, GEN M)
{
  struct huntcond_t *S = (struct huntcond_t*)E;
  GEN t1k, t2k, p1, p1inv, p2, p2inv, thetainit, R;
  GEN t1 = mkfrac(stoi(11), stoi(10)), t2 = mkfrac(stoi(13), stoi(11));
  GEN t1be, t2be, t1bemk, t2bemk, t1kmbe, t2kmbe;
  GEN F11, F12, F21, F22, P1, P2, res;
  long prec, bitprec;
  GEN k = S->k;

  thetainit = linit_get_tech(S->theta);
  bitprec = theta_get_bitprec(thetainit);
  prec = nbits2prec(bitprec);
  condset(S, M, prec);

  p1 = lfuntheta(S->thetad, t1, 0, bitprec);
  p2 = lfuntheta(S->thetad, t2, 0, bitprec);
  p1inv = lfuntheta(S->theta, ginv(t1), 0, bitprec);
  p2inv = lfuntheta(S->theta, ginv(t2), 0, bitprec);
  t1k = gpow(t1, k, prec);
  t2k = gpow(t2, k, prec);
  R = theta_get_R(thetainit);
  if (typ(R) == t_VEC)
  {
    GEN be = gmael(R, 1, 1);
    t1be = gpow(t1, be, prec); t1bemk = gdiv(gsqr(t1be), t1k);
    t2be = gpow(t2, be, prec); t2bemk = gdiv(gsqr(t2be), t2k);
    t1kmbe = gdiv(t1k, t1be);
    t2kmbe = gdiv(t2k, t2be);
  }
  else
  { /* be = k */
    t1be = t1k; t1bemk = t1k; t1kmbe = gen_1;
    t2be = t2k; t2bemk = t2k; t2kmbe = gen_1;
  }
  F11 = conj_i(gsub(gmul(gsqr(t1be), p1), p1inv));
  F12 = conj_i(gsub(gmul(gsqr(t2be), p2), p2inv));
  F21 = gsub(gmul(t1k, p1), gmul(t1bemk, p1inv));
  F22 = gsub(gmul(t2k, p2), gmul(t2bemk, p2inv));
  P1 = gsub(gmul(t1bemk, t1be), t1kmbe);
  P2 = gsub(gmul(t2bemk, t2be), t2kmbe);
  res = gdiv(gsub(gmul(P2,F21), gmul(P1,F22)),
             gsub(gmul(P2,F11), gmul(P1,F12)));
  return glog(gabs(res, prec), prec);
}

/* If flag = 0 (default) return all conductors found as integers. If
flag = 1, return the approximations, not the integers. If flag = 2,
return all, even nonintegers. */

static GEN
checkconductor(GEN v, long bit, long flag)
{
  GEN w;
  long e, j, k, l = lg(v);
  if (flag == 2) return v;
  w = cgetg(l, t_VEC);
  for (j = k = 1; j < l; j++)
  {
    GEN N = grndtoi(gel(v,j), &e);
    if (e < -bit) gel(w,k++) = flag ? gel(v,j): N;
  }
  if (k == 2) return gel(w,1);
  setlg(w,k); return w;
}

static GEN
parse_maxcond(GEN maxN)
{
  GEN M;
  if (!maxN)
    M = utoipos(10000);
  else if (typ(maxN) == t_VEC)
  {
    if (!RgV_is_ZV(maxN)) pari_err_TYPE("lfunconductor",maxN);
    return ZV_sort_shallow(maxN);
  }
  else
    M = maxN;
  return (typ(M) == t_INT)? addiu(M, 1): gceil(M);
}

GEN
lfunconductor(GEN data, GEN maxcond, long flag, long bitprec)
{
  struct huntcond_t S;
  pari_sp av = avma;
  GEN ldata = lfunmisc_to_ldata_shallow(data);
  GEN ld, r, v, theta, thetad, M, tdom, t0 = NULL, t0i = NULL;
  GEN (*eval)(void *, GEN);
  long prec;
  M = parse_maxcond(maxcond);
  r = ldata_get_residue(ldata);
  if (typ(M) == t_VEC) /* select in list */
  {
    if (lg(M) == 1) { set_avma(av); return cgetg(1,t_VEC); }
    eval = NULL; tdom = dbltor(0.7);
  }
  else if (!r) { eval = wrap1; tdom = uutoQ(10,11); }
  else
  {
    if (typ(r) == t_VEC && lg(r) > 2)
      pari_err_IMPL("multiple poles in lfunconductor");
    eval = wrap2; tdom = uutoQ(11,13);
  }
  if (eval) bitprec += bitprec/2;
  prec = nbits2prec(bitprec);
  ld = shallowcopy(ldata);
  gel(ld, 5) = eval? M: veclast(M);
  theta = lfunthetainit_i(ld, tdom, 0, bitprec);
  thetad = theta_dual(theta, ldata_get_dual(ldata));
  gel(theta,3) = shallowcopy(linit_get_tech(theta));
  S.k = ldata_get_k(ldata);
  S.theta = theta;
  S.thetad = thetad? thetad: theta;
  S.pM = &gel(linit_get_ldata(theta),5);
  S.psqrtM = &gel(linit_get_tech(theta),7);
  if (thetad)
  {
    S.pMd = &gel(linit_get_ldata(thetad),5);
    S.psqrtMd = &gel(linit_get_tech(thetad),7);
  }
  if (!eval)
  {
    long i, besti = 0, beste = -10, l = lg(M);
    t0 = uutoQ(11,10); t0i = uutoQ(10,11);
    for (i = 1; i < l; i++)
    {
      pari_sp av2 = avma;
      long e;
      condset(&S, gel(M,i), prec);
      e = lfuncheckfeq_i(theta, thetad, t0, t0i, bitprec);
      set_avma(av2);
      if (e < beste) { beste = e; besti = i; }
      else if (e == beste) beste = besti = 0; /* tie: forget */
    }
    if (!besti) { set_avma(av); return cgetg(1,t_VEC); }
    return gerepilecopy(av, mkvec2(gel(M,besti), stoi(beste)));
  }
  v = solvestep((void*)&S, eval, ghalf, M, gen_2, 14, prec);
  return gerepilecopy(av, checkconductor(v, bitprec/2, flag));
}

/* assume chi primitive */
static GEN
znchargauss_i(GEN G, GEN chi, long bitprec)
{
  GEN z, q, F = znstar_get_N(G);
  long prec;

  if (equali1(F)) return gen_1;
  prec = nbits2prec(bitprec);
  q = sqrtr_abs(itor(F, prec));
  z = lfuntheta(mkvec2(G,chi), gen_1, 0, bitprec);
  if (gexpo(z) < 10 - bitprec)
  {
    if (equaliu(F,300))
    {
      GEN z = rootsof1u_cx(25, prec);
      GEN n = znconreyexp(G, chi);
      if (equaliu(n, 131)) return gmul(q, gpowgs(z,14));
      if (equaliu(n, 71)) return gmul(q, gpowgs(z,11));
    }
    if (equaliu(F,600))
    {
      GEN z = rootsof1u_cx(25, prec);
      GEN n = znconreyexp(G, chi);
      if (equaliu(n, 491)) return gmul(q, gpowgs(z,7));
      if (equaliu(n, 11)) return gmul(q, gpowgs(z,18));
    }
    pari_err_BUG("znchargauss [ Theta(chi,1) = 0 ]");
  }
  z = gmul(gdiv(z, conj_i(z)), q);
  if (zncharisodd(G,chi)) z = mulcxI(z);
  return z;
}
static GEN
Z_radical(GEN N, long *om)
{
  GEN P = gel(Z_factor(N), 1);
  *om = lg(P)-1; return ZV_prod(P);
}
GEN
znchargauss(GEN G, GEN chi, GEN a, long bitprec)
{
  GEN v, T, N, F, b0, b1, b2, bF, a1, aF, A, r, GF, tau, B, faB, u, S;
  long omb0, prec = nbits2prec(bitprec);
  pari_sp av = avma;

  if (typ(chi) != t_COL) chi = znconreylog(G,chi);
  T = znchartoprimitive(G, chi);
  GF  = gel(T,1);
  chi = gel(T,2); /* now primitive */
  N = znstar_get_N(G);
  F = znstar_get_N(GF);
  if (equalii(N,F)) b1 = bF = gen_1;
  else
  {
    v = Z_ppio(diviiexact(N,F), F);
    bF = gel(v,2); /* (N/F, F^oo) */
    b1 = gel(v,3); /* cofactor */
  }
  if (!a) a = a1 = aF = gen_1;
  else
  {
    if (typ(a) != t_INT) pari_err_TYPE("znchargauss",a);
    a = modii(a, N);
    if (!signe(a)) { set_avma(av); return is_pm1(F)? eulerphi(N): gen_0; }
    v = Z_ppio(a, F);
    aF = gel(v,2);
    a1 = gel(v,3);
  }
  if (!equalii(aF, bF)) { set_avma(av); return gen_0; }
  b0 = Z_radical(b1, &omb0);
  b2 = diviiexact(b1, b0);
  A = dvmdii(a1, b2, &r);
  if (r != gen_0) { set_avma(av); return gen_0; }
  B = gcdii(A,b0); faB = Z_factor(B); /* squarefree */
  S = eulerphi(mkvec2(B,faB));
  if (odd(omb0 + lg(gel(faB,1))-1)) S = negi(S); /* moebius(b0/B) * phi(B) */
  S = mulii(S, mulii(aF,b2));
  tau = znchargauss_i(GF, chi, bitprec);
  u = Fp_div(b0, A, F);
  if (!equali1(u))
  {
    GEN ord = zncharorder(GF, chi), z = rootsof1_cx(ord, prec);
    tau = gmul(tau, znchareval(GF, chi, u, mkvec2(z,ord)));
  }
  return gerepileupto(av, gmul(tau, S));
}
