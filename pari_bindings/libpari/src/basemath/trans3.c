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

/********************************************************************/
/**                                                                **/
/**                   TRANSCENDENTAL FONCTIONS                     **/
/**                          (part 3)                              **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_trans

#define HALF_E 1.3591409 /* exp(1) / 2 */

/***********************************************************************/
/**                                                                   **/
/**                       BESSEL FUNCTIONS                            **/
/**                                                                   **/
/***********************************************************************/

static GEN
_abs(GEN x)
{ return gabs(gtofp(x,LOWDEFAULTPREC), LOWDEFAULTPREC); }
/* can we use asymptotic expansion ? */
static int
bessel_asymp(GEN n, GEN z, long bit)
{
  GEN Z, N;
  long t = typ(n);
  if (!is_real_t(t) && t != t_COMPLEX) return 0;
  Z = _abs(z); N = gaddgs(_abs(n), 1);
  return gcmpgs(gdiv(Z, gsqr(N)), (bit+10)/2) >= 0; }

/* Region I: 0 < Arg z <= Pi, II: -Pi < Arg z <= 0 */
static int
regI(GEN z)
{
  long s = gsigne(imag_i(z));
  return (s > 0 || (s == 0 && gsigne(real_i(z)) < 0)) ? 1 : 2;
}
/* Region 1: Re(z) >= 0, 2: Re(z) < 0, Im(z) >= 0, 3: Re(z) < 0, Im(z) < 0 */
static int
regJ(GEN z)
{
  if (gsigne(real_i(z)) >= 0) return 1;
  return gsigne(imag_i(z)) >= 0 ? 2 : 3;
}

/* Hankel's expansions:
 * a_k(n) = \prod_{0 <= j < k} (4n^2 - (2j+1)^2)
 * C(k)[n,z] = a_k(n) / (k! (8 z)^k)
 * A(z)  = exp(-z) sum_{k >= 0} C(k)
 * A(-z) = exp(z) sum_{k >= 0} (-1)^k C(k)
 * J_n(z) ~ [1] (A(z/i) / r + A(-z/i) r) / sqrt(2Pi z)
 *          [2] (A(z/i) r^3 + A(-z/i) r) / sqrt(2Pi z)
 *          [3] (A(z/i) / r + A(-z/i) / r^3) / sqrt(2Pi z)
 * Y_n(z) ~ [1] i(A(z/i) / r + A(-z/i) r) / sqrt(2Pi z)
 *          [2] i(A(z/i) (r^3-2/r) + A(-z/i) r) / sqrt(2Pi z)
 *          [3] i(-A(z/i)/r + A(-z/i)(2r-1/r^3)) / sqrt(2Pi z)
 * K_n(z) ~ A(z) Pi / sqrt(2 Pi z)
 * I_n(z) ~ [I] (A(-z) + r^2 A(z)) / sqrt(2 Pi z)
 *          [II](A(-z) + r^(-2) A(z)) / sqrt(2 Pi z) */

/* set [A(z), A(-z), exp((2*nu+1)*I*Pi/4)] */
static void
hankel_ABr(GEN *pA, GEN *pB, GEN *pr, GEN n, GEN z, long bit)
{
  GEN E, P, C, Q = gen_0, zi = ginv(gmul2n(z, 3));
  GEN K = gaddgs(_abs(n), 1), n2 = gmul2n(gsqr(n),2);
  long prec = nbits2prec(bit), B = bit + 4, m;

  P = C = real_1_bit(bit);
  for (m = 1;; m += 2)
  {
    C = gmul(C, gdivgu(gmul(gsub(n2, sqru(2*m - 1)), zi), m));
    Q = gadd(Q, C);
    C = gmul(C, gdivgu(gmul(gsub(n2, sqru(2*m + 1)), zi), m + 1));
    P = gadd(P, C);
    if (gexpo(C) < -B && gcmpgs(K, m) <= 0) break;
  }
  E = gexp(z, prec);
  *pA = gdiv(gadd(P, Q), E);
  *pB = gmul(gsub(P, Q), E);
  *pr = gexp(mulcxI(gmul(gaddgs(gmul2n(n,1), 1), Pi2n(-2, prec))), prec);
}

/* sqrt(2*Pi*z) */
static GEN
sqz(GEN z, long bit)
{
  long prec = nbits2prec(bit);
  return gsqrt(gmul(Pi2n(1, prec), z), prec);
}

static GEN
besskasymp(GEN nu, GEN z, long bit)
{
  GEN A, B, r;
  long prec = nbits2prec(bit);
  hankel_ABr(&A,&B,&r, nu, z, bit);
  return gdiv(gmul(A, mppi(prec)), sqz(z, bit));
}

static GEN
bessiasymp(GEN nu, GEN z, long bit)
{
  GEN A, B, r, R, r2;
  hankel_ABr(&A,&B,&r, nu, z, bit);
  r2 = gsqr(r);
  R = regI(z) == 1 ? gmul(A, r2) : gdiv(A, r2);
  return gdiv(gadd(B, R), sqz(z, bit));
}

static GEN
bessjasymp(GEN nu, GEN z, long bit)
{
  GEN A, B, r, R;
  long reg = regJ(z);
  hankel_ABr(&A,&B,&r, nu, mulcxmI(z), bit);
  if (reg == 1) R = gadd(gdiv(A, r), gmul(B, r));
  else if (reg == 2) R = gadd(gmul(A, gpowgs(r, 3)), gmul(B, r));
  else R = gadd(gdiv(A, r), gdiv(B, gpowgs(r, 3)));
  return gdiv(R, sqz(z, bit));
}

static GEN
bessyasymp(GEN nu, GEN z, long bit)
{
  GEN A, B, r, R;
  long reg = regJ(z);
  hankel_ABr(&A,&B,&r, nu, mulcxmI(z), bit);
  if (reg == 1) R = gsub(gmul(B, r), gdiv(A, r));
  else if (reg == 2)
    R = gadd(gmul(A, gsub(gpowgs(r, 3), gmul2n(ginv(r), 1))), gmul(B, r));
  else
    R = gsub(gmul(B, gsub(gmul2n(r, 1), ginv(gpowgs(r, 3)))), gdiv(A, r));
  return gdiv(mulcxI(R), sqz(z, bit));
}

/* n! sum_{0 <= k <= m} x^k / (k!*(k+n)!) */
static GEN
_jbessel(GEN n, GEN x, long m)
{
  pari_sp av = avma;
  GEN s = gen_1;
  long k;

  for (k = m; k >= 1; k--)
  {
    s = gaddsg(1, gdiv(gmul(x,s), gmulgu(gaddgs(n, k), k)));
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"besselj");
      s = gerepileupto(av, s);
    }
  }
  return s;
}

/* max(2, L * approximate solution to x log x = B) */
static long
bessel_get_lim(double B, double L)
{ return maxss(2, L * exp(dbllambertW0(B))); }

static GEN vjbesselh(void* E, GEN z, long prec){return jbesselh((GEN)E,z,prec);}
static GEN vjbessel(void* E, GEN z, long prec) {return jbessel((GEN)E,z,prec);}
static GEN vibessel(void* E, GEN z, long prec) {return ibessel((GEN)E,z,prec);}
static GEN vnbessel(void* E, GEN z, long prec) {return nbessel((GEN)E,z,prec);}
static GEN vkbessel(void* E, GEN z, long prec) {return kbessel((GEN)E,z,prec);}

/* if J != 0 BesselJ, else BesselI. */
static GEN
jbesselintern(GEN n, GEN z, long J, long prec)
{
  const char *f = J? "besselj": "besseli";
  long i, ki;
  pari_sp av = avma;
  GEN y;

  switch(typ(z))
  {
    case t_INT: case t_FRAC: case t_REAL: case t_COMPLEX:
    {
      int flz0 = gequal0(z);
      long lim, k, precnew, bit;
      GEN p1, p2;
      double az, L;

      i = precision(z); if (i) prec = i;
      if (flz0 && gequal0(n)) return real_1(prec);
      bit = prec2nbits(prec);
      if (bessel_asymp(n, z, bit))
      {
        GEN R = J? bessjasymp(n, z, bit): bessiasymp(n, z, bit);
        if (typ(R) == t_COMPLEX && isexactzero(imag_i(n))
                                && gsigne(real_i(z)) > 0
                                && isexactzero(imag_i(z))) R = gcopy(gel(R,1));
        return gerepileupto(av, R);
      }
      p2 = gpow(gmul2n(z,-1),n,prec);
      p2 = gdiv(p2, ggamma(gaddgs(n,1),prec));
      if (flz0) return gerepileupto(av, p2);
      az = dblmodulus(z); L = HALF_E * az;
      precnew = prec;
      if (az >= 1.0) precnew += 1 + nbits2extraprec((long)(az/M_LN2));
      if (issmall(n,&ki)) {
        k = labs(ki);
        n = utoi(k);
      } else {
        i = precision(n);
        if (i && i < precnew) n = gtofp(n,precnew);
      }
      z = gtofp(z,precnew);
      lim = bessel_get_lim(prec2nbits_mul(prec,M_LN2/2) / L, L);
      z = gmul2n(gsqr(z),-2); if (J) z = gneg(z);
      p1 = gprec_wtrunc(_jbessel(n,z,lim), prec);
      return gerepileupto(av, gmul(p2,p1));
    }

    case t_PADIC: pari_err_IMPL(stack_strcat("p-adic ",f));
    default:
    {
      long v, k, m;
      if (!(y = toser_i(z))) break;
      if (issmall(n,&ki)) n = utoi(labs(ki));
      y = gmul2n(gsqr(y),-2); if (J) y = gneg(y);
      v = valser(y);
      if (v < 0) pari_err_DOMAIN(f, "valuation", "<", gen_0, z);
      if (v == 0) pari_err_IMPL(stack_strcat(f, " around a!=0"));
      m = lg(y) - 2;
      k = m - (v >> 1);
      if (k <= 0) { set_avma(av); return scalarser(gen_1, varn(z), v); }
      setlg(y, k+2); return gerepileupto(av, _jbessel(n, y, m));
    }
  }
  return trans_evalgen(f, (void*)n, J? vjbessel: vibessel, z, prec);
}
GEN
jbessel(GEN n, GEN z, long prec) { return jbesselintern(n,z,1,prec); }
GEN
ibessel(GEN n, GEN z, long prec) { return jbesselintern(n,z,0,prec); }

/* k > 0 */
static GEN
_jbesselh(long k, GEN z, long prec)
{
  GEN s, c, p0, p1, zinv = ginv(z);
  long i;

  gsincos(z,&s,&c,prec);
  p1 = gmul(zinv,s);
  p0 = p1; p1 = gmul(zinv,gsub(p0,c));
  for (i = 2; i <= k; i++)
  {
    GEN p2 = gsub(gmul(gmulsg(2*i-1,zinv), p1), p0);
    p0 = p1; p1 = p2;
  }
  return p1;
}

/* J_{n+1/2}(z) */
GEN
jbesselh(GEN n, GEN z, long prec)
{
  long k, i;
  pari_sp av;
  GEN y;

  if (typ(n)!=t_INT) pari_err_TYPE("jbesselh",n);
  k = itos(n);
  if (k < 0) return jbessel(gadd(ghalf,n), z, prec);

  switch(typ(z))
  {
    case t_INT: case t_FRAC: case t_REAL: case t_COMPLEX:
    {
      long pr;
      GEN p1;
      if (gequal0(z))
      {
        av = avma;
        p1 = gmul(gsqrt(gdiv(z,mppi(prec)),prec),gpowgs(z,k));
        p1 = gdiv(p1, mulu_interval(k+1, 2*k+1)); /* x k! / (2k+1)! */
        return gerepileupto(av, gmul2n(p1,2*k));
      }
      if ( (pr = precision(z)) ) prec = pr;
      if (bessel_asymp(n, z, prec2nbits(prec)))
        return jbessel(gadd(ghalf,n), z, prec);
      y = cgetc(prec); av = avma;
      p1 = gsqrt(gdiv(z, Pi2n(-1,prec)), prec);
      if (!k)
        p1 = gmul(p1, gsinc(z, prec));
      else
      {
        long bits = BITS_IN_LONG + 2*k * (log2(k) -  dbllog2(z));
        if (bits > 0)
        {
          prec += nbits2extraprec(bits);
          if (pr) z = gtofp(z, prec);
        }
        p1 = gmul(p1, _jbesselh(k,z,prec));
      }
      set_avma(av); return affc_fixlg(p1, y);
    }
    case t_PADIC: pari_err_IMPL("p-adic jbesselh function");
    default:
    {
      long t, v;
      av = avma; if (!(y = toser_i(z))) break;
      if (gequal0(y)) return gerepileupto(av, gpowgs(y,k));
      v = valser(y);
      if (v < 0) pari_err_DOMAIN("besseljh","valuation", "<", gen_0, z);
      t = lg(y)-2;
      if (v) y = sertoser(y, t + (2*k+1)*v);
      if (!k)
        y = gsinc(y,prec);
      else
      {
        GEN T, a = _jbesselh(k, y, prec);
        if (v) y = sertoser(y, t + k*v); /* lower precision */
        y = gdiv(a, gpowgs(y, k));
        T = cgetg(k+1, t_VECSMALL);
        for (i = 1; i <= k; i++) T[i] = 2*i+1;
        y = gmul(y, zv_prod_Z(T));
      }
      return gerepileupto(av, y);
    }
  }
  return trans_evalgen("besseljh",(void*)n, vjbesselh, z, prec);
}

static GEN
kbessel2(GEN nu, GEN x, long prec)
{
  pari_sp av = avma;
  GEN p1, a, x2 = gshift(x,1);

  a = gtofp(gaddgs(gshift(nu,1), 1), prec);
  p1 = hyperu(gshift(a,-1), a, x2, prec);
  p1 = gmul(gmul(p1, gpow(x2,nu,prec)), sqrtr(mppi(prec)));
  return gerepileupto(av, gmul(p1, gexp(gneg(x),prec)));
}

/* special case of hyperu */
static GEN
kbessel1(GEN nu, GEN gx, long prec)
{
  GEN x, y, zf, r, u, pi, nu2;
  long bit, k, k2, n2, n, l = (typ(gx)==t_REAL)? realprec(gx): prec;
  pari_sp av;

  if (typ(nu)==t_COMPLEX) return kbessel2(nu, gx, l);
  y = cgetr(l); av = avma;
  x = gtofp(gx, l);
  nu = gtofp(nu,l); nu2 = sqrr(nu);
  shiftr_inplace(nu2,2); togglesign(nu2); /* nu2 = -4nu^2 */
  n = (long) (prec2nbits_mul(l,M_LN2) + M_PI*fabs(rtodbl(nu))) / 2;
  bit = prec2nbits(l) - 1;
  l += EXTRAPREC64;
  pi = mppi(l); n2 = n<<1; r = gmul2n(x,1);
  if (cmprs(x, n) < 0)
  {
    pari_sp av2 = avma;
    GEN q, v, c, s = real_1(l), t = real_0(l);
    for (k = n2, k2 = 2*n2-1; k > 0; k--, k2 -= 2)
    {
      GEN ak = divri(addri(nu2, sqru(k2)), mulss(n2<<2, -k));
      s = addsr(1, mulrr(ak,s));
      t = addsr(k2,mulrr(ak,t));
      if (gc_needed(av2,3)) gerepileall(av2, 2, &s,&t);
    }
    shiftr_inplace(t, -1);
    q = utor(n2, l);
    zf = sqrtr(divru(pi,n2));
    u = gprec_wensure(mulrr(zf, s), l);
    v = gprec_wensure(divrs(addrr(mulrr(t,zf),mulrr(u,nu)),-n2), l);
    for(;;)
    {
      GEN p1, e, f, d = real_1(l);
      pari_sp av3;
      c = divur(5,q); if (expo(c) >= -1) c = real2n(-1,l);
      p1 = subsr(1, divrr(r,q)); if (cmprr(c,p1)>0) c = p1;
      togglesign(c); av3 = avma;
      e = u;
      f = v;
      for (k = 1;; k++)
      {
        GEN w = addrr(gmul2n(mulur(2*k-1,u), -1), mulrr(subrs(q,k),v));
        w = addrr(w, mulrr(nu, subrr(u,gmul2n(v,1))));
        u = divru(mulrr(q,v), k);
        v = divru(w,k);
        d = mulrr(d,c);
        e = addrr(e, mulrr(d,u));
        f = addrr(f, p1 = mulrr(d,v));
        if (expo(p1) - expo(f) <= 1-prec2nbits(realprec(p1))) break;
        if (gc_needed(av3,3)) gerepileall(av3,5,&u,&v,&d,&e,&f);
      }
      u = e;
      v = f;
      q = mulrr(q, addrs(c,1));
      if (expo(r) - expo(subrr(q,r)) >= bit) break;
      gerepileall(av2, 3, &u,&v,&q);
    }
    u = mulrr(u, gpow(divru(x,n),nu,prec));
  }
  else
  {
    GEN s, zz = ginv(gmul2n(r,2));
    pari_sp av2 = avma;
    s = real_1(l);
    for (k = n2, k2 = 2*n2-1; k > 0; k--, k2 -= 2)
    {
      GEN ak = divru(mulrr(addri(nu2, sqru(k2)), zz), k);
      s = subsr(1, mulrr(ak,s));
      if (gc_needed(av2,3)) s = gerepileuptoleaf(av2, s);
    }
    zf = sqrtr(divrr(pi,r));
    u = mulrr(s, zf);
  }
  affrr(mulrr(u, mpexp(negr(x))), y);
  set_avma(av); return y;
}

/*   sum_{k=0}^m x^k (H(k)+H(k+n)) / (k! (k+n)!)
 * + sum_{k=0}^{n-1} (-x)^(k-n) (n-k-1)!/k! */
static GEN
_kbessel(long n, GEN x, long m, long prec)
{
  GEN p1, p2, s, H;
  long k, M = m + n, exact = (M <= bit_accuracy(prec));
  pari_sp av;

  H = cgetg(M+2,t_VEC); gel(H,1) = gen_0;
  if (exact)
  {
    gel(H,2) = s = gen_1;
    for (k=2; k<=M; k++) gel(H,k+1) = s = gdivgu(gaddsg(1,gmulsg(k,s)),k);
  }
  else
  {
    gel(H,2) = s = real_1(prec);
    for (k=2; k<=M; k++) gel(H,k+1) = s = divru(addsr(1,mulur(k,s)),k);
  }
  s = gadd(gel(H,m+1), gel(H,M+1)); av = avma;
  for (k = m; k > 0; k--)
  {
    s = gadd(gadd(gel(H,k),gel(H,k+n)), gdiv(gmul(x,s),mulss(k,k+n)));
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"_kbessel");
      s = gerepileupto(av, s);
    }
  }
  p1 = exact? mpfact(n): mpfactr(n,prec);
  s = gdiv(s,p1);
  if (n)
  {
    x = gneg(ginv(x));
    p2 = gmulsg(n, gdiv(x,p1));
    s = gadd(s,p2);
    for (k=n-1; k>0; k--)
    {
      p2 = gmul(p2, gmul(mulss(k,n-k),x));
      s = gadd(s,p2);
    }
  }
  return s;
}

/* N = 1: Bessel N, else Bessel K */
static GEN
kbesselintern(GEN n, GEN z, long N, long prec)
{
  const char *f = N? "besseln": "besselk";
  long i, k, ki, lim, precnew, fl2, ex, bit;
  pari_sp av = avma;
  GEN p1, p2, y, p3, pp, pm, s, c;
  double az;

  switch(typ(z))
  {
    case t_INT: case t_FRAC: case t_REAL: case t_COMPLEX:
      if (gequal0(z)) pari_err_DOMAIN(f, "argument", "=", gen_0, z);
      i = precision(z); if (i) prec = i;
      i = precision(n); if (i && prec > i) prec = i;
      bit = prec2nbits(prec);
      if (bessel_asymp(n, z, bit))
      {
        GEN R = N? bessyasymp(n, z, bit): besskasymp(n, z, bit);
        if (typ(R) == t_COMPLEX && isexactzero(imag_i(n))
                                && gsigne(real_i(z)) > 0
                                && isexactzero(imag_i(z))) R = gcopy(gel(R,1));
        return gerepileupto(av, R);
      }
      /* heuristic threshold */
      if (!N && !gequal0(n) && gexpo(z) > bit/16 + gexpo(n))
        return kbessel1(n,z,prec);
      az = dblmodulus(z); precnew = prec;
      if (az >= 1) precnew += 1 + nbits2extraprec((long)((N?az:2*az)/M_LN2));
      z = gtofp(z, precnew);
      if (issmall(n,&ki))
      {
        GEN z2 = gmul2n(z, -1), Z;
        double B, L = HALF_E * az;
        k = labs(ki);
        B = prec2nbits_mul(prec,M_LN2/2) / L;
        if (!N) B += 0.367879; /* exp(-1) */
        lim = bessel_get_lim(B, L);
        Z = gsqr(z2); if (N) Z = gneg(Z);
        p1 = gmul(gpowgs(z2,k), _kbessel(k, Z, lim, precnew));
        p2 = gadd(mpeuler(precnew), glog(z2,precnew));
        p3 = jbesselintern(stoi(k),z,N,precnew);
        p2 = gsub(gmul2n(p1,-1),gmul(p2,p3));
        p2 = gprec_wtrunc(p2, prec);
        if (N)
        {
          p2 = gdiv(p2, Pi2n(-1,prec));
          if (ki >= 0 || !odd(k)) p2 = gneg(p2);
        } else
          if (odd(k)) p2 = gneg(p2);
        return gerepilecopy(av, p2);
      }

      n = gtofp(n, precnew);
      gsincos(gmul(n,mppi(precnew)), &s,&c,precnew);
      ex = gexpo(s);
      if (ex < 0) precnew += nbits2extraprec(N? -ex: -2*ex);
      if (i && i < precnew) {
        n = gtofp(n,precnew);
        z = gtofp(z,precnew);
        gsincos(gmul(n,mppi(precnew)), &s,&c,precnew);
      }

      pp = jbesselintern(n,      z,N,precnew);
      pm = jbesselintern(gneg(n),z,N,precnew);
      if (N)
        p1 = gsub(gmul(c,pp),pm);
      else
        p1 = gmul(gsub(pm,pp), Pi2n(-1,precnew));
      p1 = gdiv(p1, s);
      return gerepilecopy(av, gprec_wtrunc(p1,prec));

    case t_PADIC: pari_err_IMPL(stack_strcat("p-adic ",f));
    default:
      if (!(y = toser_i(z))) break;
      if (issmall(n,&ki))
      {
        long v, mv, k = labs(ki), m = lg(y)-2;
        y = gmul2n(gsqr(y),-2); if (N) y = gneg(y);
        v = valser(y);
        if (v < 0) pari_err_DOMAIN(f, "valuation", "<", gen_0, z);
        if (v == 0) pari_err_IMPL(stack_strcat(f, " around a!=0"));
        mv = m - (v >> 1);
        if (mv <= 0) { set_avma(av); return scalarser(gen_1, varn(z), v); }
        setlg(y, mv+2); return gerepilecopy(av, _kbessel(k, y, m, prec));
      }
      if (!issmall(gmul2n(n,1),&ki))
        pari_err_DOMAIN(f, "2n mod Z", "!=", gen_0, n);
      k = labs(ki); n = gmul2n(stoi(k),-1);
      fl2 = (k&3)==1;
      pm = jbesselintern(gneg(n), y, N, prec);
      if (N) p1 = pm;
      else
      {
        pp = jbesselintern(n, y, N, prec);
        p2 = gpowgs(y,-k); if (fl2 == 0) p2 = gneg(p2);
        p3 = gmul2n(diviiexact(mpfact(k + 1),mpfact((k + 1) >> 1)),-(k + 1));
        p3 = gdivgu(gmul2n(gsqr(p3),1),k);
        p2 = gmul(p2,p3);
        p1 = gsub(pp,gmul(p2,pm));
      }
      return gerepileupto(av, fl2? gneg(p1): gcopy(p1));
  }
  return trans_evalgen(f, (void*)n, N? vnbessel: vkbessel, z, prec);
}

GEN
kbessel(GEN n, GEN z, long prec) { return kbesselintern(n,z,0,prec); }
GEN
ybessel(GEN n, GEN z, long prec) { return kbesselintern(n,z,1,prec); }
/* J + iN */
GEN
hbessel1(GEN n, GEN z, long prec)
{
  pari_sp av = avma;
  GEN J = jbessel(n,z,prec);
  GEN Y = ybessel(n,z,prec);
  return gerepileupto(av, gadd(J, mulcxI(Y)));
}
/* J - iN */
GEN
hbessel2(GEN n, GEN z, long prec)
{
  pari_sp av = avma;
  GEN J = jbessel(n,z,prec);
  GEN Y = ybessel(n,z,prec);
  return gerepileupto(av, gadd(J, mulcxmI(Y)));
}

static GEN
besselrefine(GEN z, GEN nu, GEN (*B)(GEN,GEN,long), long bit)
{
  GEN z0 = gprec_w(z, DEFAULTPREC), nu1 = gaddgs(nu, 1), t;
  long e, n, c, j, prec = DEFAULTPREC;

  t = gdiv(B(nu1, z0, prec), B(nu, z0, prec));
  t = gadd(z0, gdiv(gsub(gsqr(z0), gsqr(nu)), gsub(gdiv(nu, z0), t)));
  e = gexpo(t) - 2 * gexpo(z0) - 1; if (e < 0) e = 0;
  n = expu(bit + 32 - e);
  c = 1 + e + ((bit - e) >> n);
  for (j = 1; j <= n; j++)
  {
    c = 2 * c - e;
    prec = nbits2prec(c); z = gprec_w(z, prec);
    t = gdiv(B(nu1, z, prec), B(nu, z, prec));
    z = gsub(z, ginv(gsub(gdiv(nu, z), t)));
  }
  return gprec_w(z, nbits2prec(bit));
}

/* solve tan(fi) - fi = y, y >= 0; Temme's method */
static double
fi(double y)
{
  double p, pp, r;
  if (y == 0) return 0;
  if (y > 100000) return M_PI/2;
  if (y < 1)
  {
    p = pow(3*y, 1.0/3); pp = p * p;
    p = p * (1 + pp * (-210 * pp + (27 - 2*pp)) / 1575);
  }
  else
  {
    p = 1 / (y + M_PI/2); pp = p * p;
    p = M_PI/2 - p*(1 + pp*(2310 + pp*(3003 + pp*(4818 + pp*(8591 + pp*16328)))) / 3465);
  }
  pp = (y + p) * (y + p); r = (p - atan(p + y)) / pp;
  return p - (1 + pp) * r * (1 + r / (p + y));
}

static GEN
besselzero(GEN nu, long n, GEN (*B)(GEN,GEN,long), long bit)
{
  pari_sp av = avma;
  long prec = nbits2prec(bit);
  int J = B == jbessel;
  GEN z;
  if (n <= 0) pari_err_DOMAIN("besselzero", "n", "<=", gen_0, stoi(n));
  if (n > LONG_MAX / 4) pari_err_OVERFLOW("besselzero");
  if (is_real_t(typ(nu)) && gsigne(nu) >= 0)
  { /* Temme */
    double x, c, b, a = gtodouble(nu), t = J? 0.25: 0.75;
    if (n >= 3*a - 8)
    {
      double aa = a*a, mu = 4*aa, mu2 = mu*mu, p, p0, p1, q1;
      p = 7 * mu - 31; p0 = mu-1;
      if (1 + p == p) /* p large */
        p1 = q1 = 0;
      else
      {
        p1 = 4 * (253 * mu2 - 3722 * mu + 17869) / (15 * p);
        q1 = 1.6 * (83 * mu2 - 982 * mu + 3779) / p;
      }
      b = (n + a/2 - t) * M_PI;
      c = 1 / (64 * b * b);
      x = b - p0 * (1 - p1 * c) / (8 * b * (1 - q1 * c));
    }
    else
    {
      double u, v, w, xx, bb = a >= 3? pow(a, -2./3): 1;
      if (n == 1)
        x = J? -2.33811: -1.17371;
      else
      {
        double pp1 = 5./48, qq1 = -5./36, y = 3./8 * M_PI;
        x = 4 * y * (n - t); v = 1 / (x*x);
        x = - pow(x, 2.0/3) * (1 + v * (pp1 + qq1 * v));
      }
      u = x * bb; v = fi(2.0/3 * pow(-u, 1.5));
      w = 1 / cos(v); xx = 1 - w*w; c = sqrt(u/xx);
      x = w * (a + c / (48*a*u) * (-5/u-c * (-10/xx + 6)));
    }
    z = dbltor(x);
  }
  else
  { /* generic, hope for the best */
    long a = 4 * n - (J? 1: 3);
    GEN b, m;
    b = gmul(mppi(prec), gmul2n(gaddgs(gmul2n(nu,1), a), -2));
    m = gmul2n(gsqr(nu),2);
    z = gsub(b, gdiv(gsubgs(m, 1), gmul2n(b, 3)));
  }
  return gerepilecopy(av, besselrefine(z, nu, B, bit));
}
GEN
besseljzero(GEN nu, long k, long b) { return besselzero(nu, k, jbessel, b); }
GEN
besselyzero(GEN nu, long k, long b) { return besselzero(nu, k, ybessel, b); }

/***********************************************************************/
/**                    INCOMPLETE GAMMA FUNCTION                      **/
/***********************************************************************/
/* mx ~ |x|, b = bit accuracy */
static int
gamma_use_asymp(GEN x, long b)
{
  long e;
  if (is_real_t(typ(x)))
  {
    pari_sp av = avma;
    return gc_int(av, gcmpgs(R_abs_shallow(x), 3*b / 4) >= 0);
  }
  e = gexpo(x); return e >= b || dblmodulus(x) >= 3*b / 4;
}
/* x a t_REAL */
static GEN
eint1r_asymp(GEN x, GEN expx, long prec)
{
  pari_sp av = avma, av2;
  GEN S, q, z, ix;
  long oldeq = LONG_MAX, esx = -prec2nbits(prec), j;

  if (realprec(x) < prec + EXTRAPREC64) x = rtor(x, prec+EXTRAPREC64);
  ix = invr(x); q = z = negr(ix);
  av2 = avma; S = addrs(q, 1);
  for (j = 2;; j++)
  {
    long eq = expo(q); if (eq < esx) break;
    if ((j & 3) == 0)
    { /* guard against divergence */
      if (eq > oldeq) return gc_NULL(av); /* regressing, abort */
      oldeq = eq;
    }
    q = mulrr(q, mulru(z, j)); S = addrr(S, q);
    if (gc_needed(av2, 1)) gerepileall(av2, 2, &S, &q);
  }
  if (DEBUGLEVEL > 2) err_printf("eint1: using asymp\n");
  S = expx? divrr(S, expx): mulrr(S, mpexp(negr(x)));
  return gerepileuptoleaf(av, mulrr(S, ix));
}
/* cf incgam_asymp(0, x); z = -1/x
 *   exp(-x)/x * (1 + z + 2! z^2 + ...) */
static GEN
eint1_asymp(GEN x, GEN expx, long prec)
{
  pari_sp av = avma, av2;
  GEN S, q, z, ix;
  long oldeq = LONG_MAX, esx = -prec2nbits(prec), j;

  if (typ(x) != t_REAL) x = gtofp(x, prec+EXTRAPREC64);
  if (typ(x) == t_REAL) return eint1r_asymp(x, expx, prec);
  ix = ginv(x); q = z = gneg_i(ix);
  av2 = avma; S = gaddgs(q, 1);
  for (j = 2;; j++)
  {
    long eq = gexpo(q); if (eq < esx) break;
    if ((j & 3) == 0)
    { /* guard against divergence */
      if (eq > oldeq) return gc_NULL(av); /* regressing, abort */
      oldeq = eq;
    }
    q = gmul(q, gmulgu(z, j)); S = gadd(S, q);
    if (gc_needed(av2, 1)) gerepileall(av2, 2, &S, &q);
  }
  if (DEBUGLEVEL > 2) err_printf("eint1: using asymp\n");
  S = expx? gdiv(S, expx): gmul(S, gexp(gneg_i(x), prec));
  return gerepileupto(av, gmul(S, ix));
}

/* eint1(x) = incgam(0, x); typ(x) = t_REAL, x > 0 */
static GEN
eint1p(GEN x, GEN expx)
{
  pari_sp av;
  long l = realprec(x), bit = prec2nbits(l), prec, i;
  double mx;
  GEN z, S, t, H, run;

  if (gamma_use_asymp(x, bit)
      && (z = eint1r_asymp(x, expx, l))) return z;
  mx = rtodbl(x);
  prec = l + nbits2extraprec((mx+log(mx))/M_LN2 + 10);
  bit = prec2nbits(prec);
  run = real_1(prec); x = rtor(x, prec);
  av = avma; S = z = t = H = run;
  for (i = 2; expo(S) - expo(t) <= bit; i++)
  {
    H = addrr(H, divru(run,i)); /* H = sum_{k<=i} 1/k */
    z = divru(mulrr(x,z), i);   /* z = x^(i-1)/i! */
    t = mulrr(z, H); S = addrr(S, t);
    if ((i & 0x1ff) == 0) gerepileall(av, 4, &z,&t,&S,&H);
  }
  return subrr(mulrr(x, divrr(S,expx? expx: mpexp(x))),
               addrr(mplog(x), mpeuler(prec)));
}
/* eint1(x) = incgam(0, x); typ(x) = t_REAL, x < 0
 * rewritten from code contributed by Manfred Radimersky */
static GEN
eint1m(GEN x, GEN expx)
{
  GEN p1, q, S, y, z = cgetg(3, t_COMPLEX);
  long l  = realprec(x), n  = prec2nbits(l), j;
  pari_sp av = avma;

  y  = rtor(x, l + EXTRAPREC64); setsigne(y,1); /* |x| */
  if (gamma_use_asymp(y, n))
  { /* ~eint1_asymp: asymptotic expansion */
    p1 = q = invr(y); S = addrs(q, 1);
    for (j = 2; expo(q) >= -n; j++) {
      q = mulrr(q, mulru(p1, j));
      S = addrr(S, q);
    }
    y  = mulrr(p1, expx? divrr(S, expx): mulrr(S, mpexp(y)));
  }
  else
  {
    p1 = q = S = y;
    for (j = 2; expo(q) - expo(S) >= -n; j++) {
      p1 = mulrr(y, divru(p1, j)); /* (-x)^j/j! */
      q = divru(p1, j);
      S = addrr(S, q);
    }
    y  = addrr(S, addrr(logr_abs(x), mpeuler(l)));
  }
  y = gerepileuptoleaf(av, y); togglesign(y);
  gel(z, 1) = y;
  y = mppi(l); setsigne(y, -1);
  gel(z, 2) = y; return z;
}

/* real(z*log(z)-z), z = x+iy */
static double
mygamma(double x, double y)
{
  if (x == 0.) return -(M_PI/2)*fabs(y);
  return (x/2)*log(x*x+y*y)-x-y*atan(y/x);
}

/* x^s exp(-x) */
static GEN
expmx_xs(GEN s, GEN x, GEN logx, long prec)
{
  GEN z;
  long ts = typ(s);
  if (ts == t_INT || (ts == t_FRAC && absequaliu(gel(s,2), 2)))
    z = gmul(gexp(gneg(x), prec), gpow(x, s, prec));
  else
    z = gexp(gsub(gmul(s, logx? logx: glog(x,prec+EXTRAPREC64)), x), prec);
  return z;
}

/* Not yet: doesn't work at low accuracy
#define INCGAM_CF
*/

#ifdef INCGAM_CF
/* Is s very close to a nonpositive integer ? */
static int
isgammapole(GEN s, long bitprec)
{
  pari_sp av = avma;
  GEN t = imag_i(s);
  long e, b = bitprec - 10;

  if (gexpo(t) > - b) return 0;
  s = real_i(s);
  if (gsigne(s) > 0 && gexpo(s) > -b) return 0;
  (void)grndtoi(s, &e); return gc_bool(av, e < -b);
}

/* incgam using the continued fraction. x a t_REAL or t_COMPLEX, mx ~ |x|.
 * Assume precision(s), precision(x) >= prec */
static GEN
incgam_cf(GEN s, GEN x, double mx, long prec)
{
  GEN ms, y, S;
  long n, i, j, LS, bitprec = prec2nbits(prec);
  double rs, is, m;

  if (typ(s) == t_COMPLEX)
  {
    rs = gtodouble(gel(s,1));
    is = gtodouble(gel(s,2));
  }
  else
  {
    rs = gtodouble(s);
    is = 0.;
  }
  if (isgammapole(s, bitprec)) LS = 0;
  else
  {
    double bit,  LGS = mygamma(rs,is);
    LS = LGS <= 0 ? 0: ceil(LGS);
    bit = (LGS - (rs-1)*log(mx) + mx)/M_LN2;
    if (bit > 0)
    {
      prec += nbits2extraprec((long)bit);
      x = gtofp(x, prec);
      if (isinexactreal(s)) s = gtofp(s, prec);
    }
  }
  /* |ln(2*gamma(s)*sin(s*Pi))| <= ln(2) + |lngamma(s)| + |Im(s)*Pi|*/
  m = bitprec*M_LN2 + LS + M_LN2 + fabs(is)*M_PI + mx;
  if (rs < 1) m += (1 - rs)*log(mx);
  m /= 4;
  n = (long)(1 + m*m/mx);
  y = expmx_xs(gsubgs(s,1), x, NULL, prec);
  if (rs >= 0 && bitprec >= 512)
  {
    GEN A = cgetg(n+1, t_VEC), B = cgetg(n+1, t_VEC);
    ms = gsubsg(1, s);
    for (j = 1; j <= n; ++j)
    {
      gel(A,j) = ms;
      gel(B,j) = gmulsg(j, gsubgs(s,j));
      ms = gaddgs(ms, 2);
    }
    S = contfraceval_inv(mkvec2(A,B), x, -1);
  }
  else
  {
    GEN x_s = gsub(x, s);
    pari_sp av2 = avma;
    S = gdiv(gsubgs(s,n), gaddgs(x_s,n<<1));
    for (i=n-1; i >= 1; i--)
    {
      S = gdiv(gsubgs(s,i), gadd(gaddgs(x_s,i<<1),gmulsg(i,S)));
      if (gc_needed(av2,3))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"incgam_cf");
        S = gerepileupto(av2, S);
      }
    }
    S = gaddgs(S,1);
  }
  return gmul(y, S);
}
#endif

static double
findextraincgam(GEN s, GEN x)
{
  double sig = gtodouble(real_i(s)), t = gtodouble(imag_i(s));
  double xr = gtodouble(real_i(x)), xi = gtodouble(imag_i(x));
  double exd = 0., Nx = xr*xr + xi*xi, D = Nx - t*t;
  long n;

  if (xr < 0)
  {
    long ex = gexpo(x);
    if (ex > 0 && ex > gexpo(s)) exd = sqrt(Nx)*log(Nx)/2; /* |x| log |x| */
  }
  if (D <= 0.) return exd;
  n = (long)(sqrt(D)-sig);
  if (n <= 0) return exd;
  return maxdd(exd, (n*log(Nx)/2 - mygamma(sig+n, t) + mygamma(sig, t)) / M_LN2);
}

/* use exp(-x) * (x^s/s) * sum_{k >= 0} x^k / prod(i=1, k, s+i) */
static GEN
incgamc_i(GEN s, GEN x, long *ptexd, long prec)
{
  GEN S, t, y;
  long l, n, i, exd;
  pari_sp av = avma, av2;

  if (gequal0(x))
  {
    if (ptexd) *ptexd = 0.;
    return gtofp(x, prec);
  }
  l = precision(x);
  if (!l) l = prec;
  n = -prec2nbits(l)-1;
  exd = (long)findextraincgam(s, x);
  if (ptexd) *ptexd = exd;
  if (exd > 0)
  {
    long p = l + nbits2extraprec(exd);
    x = gtofp(x, p);
    if (isinexactreal(s)) s = gtofp(s, p);
  }
  else x = gtofp(x, l+EXTRAPREC64);
  av2 = avma;
  S = gdiv(x, gaddsg(1,s));
  t = gaddsg(1, S);
  for (i=2; gexpo(S) >= n; i++)
  {
    S = gdiv(gmul(x,S), gaddsg(i,s)); /* x^i / ((s+1)...(s+i)) */
    t = gadd(S,t);
    if (gc_needed(av2,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"incgamc");
      gerepileall(av2, 2, &S, &t);
    }
  }
  y = expmx_xs(s, x, NULL, prec);
  return gerepileupto(av, gmul(gdiv(y,s), t));
}

GEN
incgamc(GEN s, GEN x, long prec)
{ return incgamc_i(s, x, NULL, prec); }

/* incgamma using asymptotic expansion:
 *   exp(-x)x^(s-1)(1 + (s-1)/x + (s-1)(s-2)/x^2 + ...) */
static GEN
incgam_asymp(GEN s, GEN x, long prec)
{
  pari_sp av = avma, av2;
  GEN S, q, cox, invx;
  long oldeq = LONG_MAX, eq, esx, j;
  int flint = (typ(s) == t_INT && signe(s) > 0);

  x = gtofp(x,prec+EXTRAPREC64);
  invx = ginv(x);
  esx = -prec2nbits(prec);
  av2 = avma;
  q = gmul(gsubgs(s,1), invx);
  S = gaddgs(q, 1);
  for (j = 2;; j++)
  {
    eq = gexpo(q); if (eq < esx) break;
    if (!flint && (j & 3) == 0)
    { /* guard against divergence */
      if (eq > oldeq) return gc_NULL(av); /* regressing, abort */
      oldeq = eq;
    }
    q = gmul(q, gmul(gsubgs(s,j), invx));
    S = gadd(S, q);
    if (gc_needed(av2, 1)) gerepileall(av2, 2, &S, &q);
  }
  if (DEBUGLEVEL > 2) err_printf("incgam: using asymp\n");
  cox = expmx_xs(gsubgs(s,1), x, NULL, prec);
  return gerepileupto(av, gmul(cox, S));
}

/* gasx = incgam(s-n,x). Compute incgam(s,x)
 * = (s-1)(s-2)...(s-n)gasx + exp(-x)x^(s-1) *
 *   (1 + (s-1)/x + ... + (s-1)(s-2)...(s-n+1)/x^(n-1)) */
static GEN
incgam_asymp_partial(GEN s, GEN x, GEN gasx, long n, long prec)
{
  pari_sp av;
  GEN S, q, cox, invx, s1 = gsubgs(s, 1), sprod;
  long j;
  cox = expmx_xs(s1, x, NULL, prec);
  if (n == 1) return gadd(cox, gmul(s1, gasx));
  invx = ginv(x);
  av = avma;
  q = gmul(s1, invx);
  S = gaddgs(q, 1);
  for (j = 2; j < n; j++)
  {
    q = gmul(q, gmul(gsubgs(s, j), invx));
    S = gadd(S, q);
    if (gc_needed(av, 2)) gerepileall(av, 2, &S, &q);
  }
  sprod = gmul(gmul(q, gpowgs(x, n-1)), gsubgs(s, n));
  return gadd(gmul(cox, S), gmul(sprod, gasx));
}

/* Assume s != 0; called when Re(s) <= 1/2 */
static GEN
incgamspec(GEN s, GEN x, GEN g, long prec)
{
  GEN q, S, cox = gen_0, P, sk, S1, S2, S3, F3, logx, mx;
  long n, esk, E, k = itos(ground(gneg(real_i(s)))); /* >= 0 */

  if (k && gexpo(x) > 0)
  {
    GEN xk = gdivgu(x, k);
    long bitprec = prec2nbits(prec);
    double d = (gexpo(xk) > bitprec)? bitprec*M_LN2: log(dblmodulus(xk));
    d = k * (d + 1) / M_LN2;
    if (d > 0) prec += nbits2extraprec((long)d);
    if (isinexactreal(s)) s = gtofp(s, prec);
  }
  x = gtofp(x, maxss(precision(x), prec) + EXTRAPREC64);
  sk = gaddgs(s, k); /* |Re(sk)| <= 1/2 */
  logx = glog(x, prec);
  mx = gneg(x);
  if (k == 0) { S = gen_0; P = gen_1; }
  else
  {
    long j;
    q = ginv(s); S = q; P = s;
    for (j = 1; j < k; j++)
    {
      GEN sj = gaddgs(s, j);
      q = gmul(q, gdiv(x, sj));
      S = gadd(S, q);
      P = gmul(P, sj);
    }
    cox = expmx_xs(s, x, logx, prec); /* x^s exp(-x) */
    S = gmul(S, gneg(cox));
  }
  if (k && gequal0(sk))
    return gadd(S, gdiv(eint1(x, prec), P));
  esk = gexpo(sk);
  if (esk > -7)
  {
    GEN a, b, PG = gmul(sk, P);
    if (g) g = gmul(g, PG);
    a = incgam0(gaddgs(sk,1), x, g, prec);
    if (k == 0) cox = expmx_xs(s, x, logx, prec);
    b = gmul(gpowgs(x, k), cox);
    return gadd(S, gdiv(gsub(a, b), PG));
  }
  E = prec2nbits(prec) + 1;
  if (gexpo(x) > 0)
  {
    long X = (long)(dblmodulus(x)/M_LN2);
    prec += 2*nbits2extraprec(X);
    x = gtofp(x, prec); mx = gneg(x);
    logx = glog(x, prec); sk = gtofp(sk, prec);
    E += X;
  }
  if (isinexactreal(sk)) sk = gtofp(sk, prec+EXTRAPREC64);
  /* |sk| < 2^-7 is small, guard against cancellation */
  F3 = gexpm1(gmul(sk, logx), prec);
  /* ( gamma(1+sk) - exp(sk log(x))) ) / sk */
  S1 = gdiv(gsub(ggamma1m1(sk, prec+EXTRAPREC64), F3), sk);
  q = x; S3 = gdiv(x, gaddsg(1,sk));
  for (n = 2; gexpo(q) - gexpo(S3) > -E; n++)
  {
    q = gmul(q, gdivgu(mx, n));
    S3 = gadd(S3, gdiv(q, gaddsg(n, sk)));
  }
  S2 = gadd(gadd(S1, S3), gmul(F3, S3));
  return gadd(S, gdiv(S2, P));
}

/* return |x| */
double
dblmodulus(GEN x)
{
  if (typ(x) == t_COMPLEX)
  {
    double a = gtodouble(gel(x,1));
    double b = gtodouble(gel(x,2));
    return sqrt(a*a + b*b);
  }
  else
    return fabs(gtodouble(x));
}

/* Driver routine. If g != NULL, assume that g=gamma(s,prec). */
GEN
incgam0(GEN s, GEN x, GEN g, long prec)
{
  pari_sp av;
  long E, l;
  GEN z, rs, is;

  if (gequal0(x)) return g? gcopy(g): ggamma(s,prec);
  if (gequal0(s)) return eint1(x, prec);
  l = precision(s); if (!l) l = prec;
  E = prec2nbits(l);
  if (gamma_use_asymp(x, E) ||
      (typ(s) == t_INT && signe(s) > 0 && gexpo(x) >= expi(s)))
    if ((z = incgam_asymp(s, x, l))) return z;
  av = avma; E++;
  rs = real_i(s);
  is = imag_i(s);
#ifdef INCGAM_CF
  /* Can one use continued fraction ? */
  if (gequal0(is) && gequal0(imag_i(x)) && gsigne(x) > 0)
  {
    double sd = gtodouble(rs), LB, UB;
    double xd = gtodouble(real_i(x));
    if (sd > 0) {
      LB = 15 + 0.1205*E;
      UB = 5 + 0.752*E;
    } else {
      LB = -6 + 0.1205*E;
      UB = 5 + 0.752*E + fabs(sd)/54.;
    }
    if (xd >= LB && xd <= UB)
    {
      if (DEBUGLEVEL > 2) err_printf("incgam: using continued fraction\n");
      return gerepileupto(av, incgam_cf(s, x, xd, prec));
    }
  }
#endif
  if (gsigne(rs) > 0 && gexpo(rs) >= -1)
  { /* use complementary incomplete gamma */
    long n, egs, exd, precg, es = gexpo(s);
    if (es < 0) {
      l += nbits2extraprec(-es) + 1;
      x = gtofp(x, l);
      if (isinexactreal(s)) s = gtofp(s, l);
    }
    n = itos(gceil(rs));
    if (n > 100)
    {
      GEN gasx;
      n -= 100;
      if (es > 0)
      {
        es = mygamma(gtodouble(rs) - n, gtodouble(is)) / M_LN2;
        if (es > 0)
        {
          l += nbits2extraprec(es);
          x = gtofp(x, l);
          if (isinexactreal(s)) s = gtofp(s, l);
        }
      }
      gasx = incgam0(gsubgs(s, n), x, NULL, prec);
      return gerepileupto(av, incgam_asymp_partial(s, x, gasx, n, prec));
    }
    if (DEBUGLEVEL > 2) err_printf("incgam: using power series 1\n");
    /* egs ~ expo(gamma(s)) */
    precg = g? precision(g): 0;
    egs = g? gexpo(g): (long)(mygamma(gtodouble(rs), gtodouble(is)) / M_LN2);
    if (egs > 0) {
      l += nbits2extraprec(egs) + 1;
      x = gtofp(x, l);
      if (isinexactreal(s)) s = gtofp(s, l);
      if (precg < l) g = NULL;
    }
    z = incgamc_i(s, x, &exd, l);
    if (exd > 0)
    {
      l += nbits2extraprec(exd);
      if (isinexactreal(s)) s = gtofp(s, l);
      if (precg < l) g = NULL;
    }
    else
    { /* gamma(s) negligible ? Compute to lower accuracy */
      long e = gexpo(z) - egs;
      if (e > 3)
      {
        E -= e;
        if (E <= 0) g = gen_0; else if (!g) g = ggamma(s, nbits2prec(E));
      }
    }
    /* worry about possible cancellation */
    if (!g) g = ggamma(s, maxss(l,precision(z)));
    return gerepileupto(av, gsub(g,z));
  }
  if (DEBUGLEVEL > 2) err_printf("incgam: using power series 2\n");
  return gerepileupto(av, incgamspec(s, x, g, l));
}

GEN
incgam(GEN s, GEN x, long prec) { return incgam0(s, x, NULL, prec); }

/* x a t_REAL */
GEN
mpeint1(GEN x, GEN expx)
{
  long s = signe(x);
  pari_sp av;
  GEN z;
  if (!s) pari_err_DOMAIN("eint1", "x","=",gen_0, x);
  if (s < 0) return eint1m(x, expx);
  z = cgetr(realprec(x));
  av = avma; affrr(eint1p(x, expx), z);
  set_avma(av); return z;
}

static GEN
cxeint1(GEN x, long prec)
{
  pari_sp av = avma, av2;
  GEN q, S, run, z, H;
  long n, E = prec2nbits(prec);

  if (gamma_use_asymp(x, E) && (z = eint1_asymp(x, NULL, prec))) return z;
  E++;
  if (gexpo(x) > 0)
  { /* take cancellation into account, log2(\sum |x|^n / n!) = |x| / log(2) */
    double dbx = dblmodulus(x);
    long X = (long)((dbx + log(dbx))/M_LN2 + 10);
    prec += nbits2extraprec(X);
    x = gtofp(x, prec); E += X;
  }
  if (DEBUGLEVEL > 2) err_printf("eint1: using power series\n");
  run = real_1(prec);
  av2 = avma;
  S = z = q = H = run;
  for (n = 2; gexpo(q) - gexpo(S) >= -E; n++)
  {
    H = addrr(H, divru(run, n)); /* H = sum_{k<=n} 1/k */
    z = gdivgu(gmul(x,z), n);   /* z = x^(n-1)/n! */
    q = gmul(z, H); S = gadd(S, q);
    if ((n & 0x1ff) == 0) gerepileall(av2, 4, &z, &q, &S, &H);
  }
  S = gmul(gmul(x, S), gexp(gneg_i(x), prec));
  return gerepileupto(av, gsub(S, gadd(glog(x, prec), mpeuler(prec))));
}

GEN
eint1(GEN x, long prec)
{
  switch(typ(x))
  {
    case t_COMPLEX: return cxeint1(x, prec);
    case t_REAL: break;
    default: x = gtofp(x, prec);
  }
  return mpeint1(x,NULL);
}

GEN
veceint1(GEN C, GEN nmax, long prec)
{
  if (!nmax) return eint1(C,prec);
  if (typ(nmax) != t_INT) pari_err_TYPE("veceint1",nmax);
  if (typ(C) != t_REAL) {
    C = gtofp(C, prec);
    if (typ(C) != t_REAL) pari_err_TYPE("veceint1",C);
  }
  if (signe(C) <= 0) pari_err_DOMAIN("veceint1", "argument", "<=", gen_0,C);
  return mpveceint1(C, NULL, itos(nmax));
}

/* j > 0, a t_REAL. Return sum_{m >= 0} a^m / j(j+1)...(j+m)).
 * Stop when expo(summand) < E; note that s(j-1) = (a s(j) + 1) / (j-1). */
static GEN
mp_sum_j(GEN a, long j, long E, long prec)
{
  pari_sp av = avma;
  GEN q = divru(real_1(prec), j), s = q;
  long m;
  for (m = 0;; m++)
  {
    if (expo(q) < E) break;
    q = mulrr(q, divru(a, m+j));
    s = addrr(s, q);
  }
  return gerepileuptoleaf(av, s);
}
/* Return the s_a(j), j <= J */
static GEN
sum_jall(GEN a, long J, long prec)
{
  GEN s = cgetg(J+1, t_VEC);
  long j, E = -prec2nbits(prec) - 5;
  gel(s, J) = mp_sum_j(a, J, E, prec);
  for (j = J-1; j; j--)
    gel(s,j) = divru(addrs(mulrr(a, gel(s,j+1)), 1), j);
  return s;
}

/* T a dense t_POL with t_REAL coeffs. Return T(n) [faster than poleval] */
static GEN
rX_s_eval(GEN T, long n)
{
  long i = lg(T)-1;
  GEN c = gel(T,i);
  for (i--; i>=2; i--) c = gadd(mulrs(c,n),gel(T,i));
  return c;
}

/* C>0 t_REAL, eC = exp(C). Return eint1(n*C) for 1<=n<=N. Absolute accuracy */
GEN
mpveceint1(GEN C, GEN eC, long N)
{
  const long prec = realprec(C);
  long Nmin = 15; /* >= 1. E.g. between 10 and 30, but little effect */
  GEN en, v, w = cgetg(N+1, t_VEC);
  pari_sp av0;
  double DL;
  long n, j, jmax, jmin;
  if (!N) return w;
  for (n = 1; n <= N; n++) gel(w,n) = cgetr(prec);
  av0 = avma;
  if (N < Nmin) Nmin = N;
  if (!eC) eC = mpexp(C);
  en = eC; affrr(eint1p(C, en), gel(w,1));
  for (n = 2; n <= Nmin; n++)
  {
    pari_sp av2;
    en = mulrr(en,eC); /* exp(n C) */
    av2 = avma;
    affrr(eint1p(mulru(C,n), en), gel(w,n));
    set_avma(av2);
  }
  if (Nmin == N) { set_avma(av0); return w; }

  DL = prec2nbits_mul(prec, M_LN2) + 5;
  jmin = ceil(DL/log((double)N)) + 1;
  jmax = ceil(DL/log((double)Nmin)) + 1;
  v = sum_jall(C, jmax, prec);
  en = powrs(eC, -N); /* exp(-N C) */
  affrr(eint1p(mulru(C,N), invr(en)), gel(w,N));
  for (j = jmin, n = N-1; j <= jmax; j++)
  {
    long limN = maxss((long)ceil(exp(DL/j)), Nmin);
    GEN polsh;
    setlg(v, j+1);
    polsh = RgV_to_RgX_reverse(v, 0);
    for (; n >= limN; n--)
    {
      pari_sp av2 = avma;
      GEN S = divri(mulrr(en, rX_s_eval(polsh, -n)), powuu(n,j));
      /* w[n+1] - exp(-n C) * polsh(-n) / (-n)^j */
      GEN c = odd(j)? addrr(gel(w,n+1), S) : subrr(gel(w,n+1), S);
      affrr(c, gel(w,n)); set_avma(av2);
      en = mulrr(en,eC); /* exp(-n C) */
    }
  }
  set_avma(av0); return w;
}

/* erfc via numerical integration : assume real(x)>=1 */
static GEN
cxerfc_r1(GEN x, long prec)
{
  GEN h, h2, eh2, denom, res, lambda;
  long u, v;
  const double D = prec2nbits_mul(prec, M_LN2);
  const long npoints = (long)ceil(D/M_PI)+1;
  pari_sp av = avma;
  {
    double t = exp(-2*M_PI*M_PI/D); /* ~exp(-2*h^2) */
    v = 30; /* bits that fit in both long and double mantissa */
    u = (long)floor(t*(1L<<v));
    /* define exp(-2*h^2) to be u*2^(-v) */
  }
  incrprec(prec);
  x = gtofp(x,prec);
  eh2 = sqrtr_abs(rtor(shiftr(dbltor(u),-v),prec));
  h2 = negr(logr_abs(eh2));
  h = sqrtr_abs(h2);
  lambda = gdiv(x,h);
  denom = gsqr(lambda);
  { /* res = h/x + 2*x*h*sum(k=1,npoints,exp(-(k*h)^2)/(lambda^2+k^2)); */
    GEN Uk; /* = exp(-(kh)^2) */
    GEN Vk = eh2;/* = exp(-(2k+1)h^2) */
    pari_sp av2 = avma;
    long k;
    /* k = 0 moved out for efficiency */
    denom = gaddsg(1,denom);
    Uk = Vk;
    Vk = mulur(u,Vk); shiftr_inplace(Vk, -v);
    res = gdiv(Uk, denom);
    for (k = 1; k < npoints; k++)
    {
      if ((k & 255) == 0) gerepileall(av2,4,&denom,&Uk,&Vk,&res);
      denom = gaddsg(2*k+1,denom);
      Uk = mpmul(Uk,Vk);
      Vk = mulur(u,Vk); shiftr_inplace(Vk, -v);
      res = gadd(res, gdiv(Uk, denom));
    }
  }
  res = gmul(res, gshift(lambda,1));
  /* 0 term : */
  res = gadd(res, ginv(lambda));
  res = gmul(res, gdiv(gexp(gneg(gsqr(x)), prec), mppi(prec)));
  if (rtodbl(real_i(x)) < sqrt(D))
  {
    GEN t = gmul(divrr(Pi2n(1,prec),h), x);
    res = gsub(res, gdivsg(2, cxexpm1(t, prec)));
  }
  return gerepileupto(av,res);
}

static GEN
sererfc(GEN x, long prec)
{
  GEN u, z = invr(sqrtr_abs(Pi2n(-2,prec)));
  setsigne(z, -1); /* -2/sqrt(Pi) */
  z = gmul(z, integser(gmul(derivser(x), gexp(gneg(gsqr(x)), prec))));
  u = polcoef_i(x, 0, varn(x));
  if (!gcmp0(u)) z = gadd(z, gerfc(u, prec));
  return z;
}

GEN
gerfc(GEN x, long prec)
{
  GEN z, xr, xi, res;
  long s;
  pari_sp av;

  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC: case t_COMPLEX:
      break;
    default:
      av = avma;
      if ((z = toser_i(x))) return gerepileupto(av, sererfc(z,prec));
      return trans_eval("erfc",gerfc,x,prec);
  }
  /* x a complex scalar */
  x = trans_fix_arg(&prec,&x,&xr,&xi, &av,&res);
  s = signe(xr);
  if (s > 0 || (s == 0 && signe(xi) >= 0)) {
    if (cmprs(xr, 1) > 0) /* use numerical integration */
      z = cxerfc_r1(x, prec);
    else
    { /* erfc(x) = incgam(1/2,x^2)/sqrt(Pi) */
      GEN sqrtpi = sqrtr(mppi(prec));
      z = incgam0(ghalf, gsqr(x), sqrtpi, prec);
      z = gdiv(z, sqrtpi);
    }
  }
  else
  { /* erfc(-x)=2-erfc(x) */
    /* FIXME could decrease prec
    long size = nbits2extraprec((imag(x)^2-real(x)^2)/log(2));
    prec = size > 0 ? prec : prec + size;
    */
    /* NOT gsubsg(2, ...) : would create a result of
     * huge accuracy if re(x)>>1, rounded to 2 by subsequent affc_fixlg... */
    z = gsub(real2n(1,prec+EXTRAPREC64), gerfc(gneg(x), prec));
  }
  set_avma(av); return affc_fixlg(z, res);
}

/***********************************************************************/
/**                                                                   **/
/**                      RIEMANN ZETA FUNCTION                        **/
/**                                                                   **/
/***********************************************************************/
static const double log2PI = 1.83787706641;

static double
get_xinf(double beta)
{
  const double maxbeta = 0.06415003; /* 3^(-2.5) */
  double x0, y0, x1;

  if (beta < maxbeta) return beta + pow(3*beta, 1.0/3.0);
  x0 = beta + M_PI/2.0;
  for(;;)
  {
    y0 = x0*x0;
    x1 = (beta+atan(x0)) * (1+y0) / y0 - 1/x0;
    if (0.99*x0 < x1) return x1;
    x0 = x1;
  }
}
/* optimize for zeta( s + it, prec ), assume |s-1| > 0.1
 * (if gexpo(u = s-1) < -5, we use the functional equation s->1-s) */
static void
optim_zeta(GEN S, long prec, long *pp, long *pn)
{
  double s, t, alpha, beta, n, B;
  long p;
  if (typ(S) == t_REAL) {
    s = rtodbl(S);
    t = 0.;
  } else {
    s = rtodbl(gel(S,1));
    t = fabs( rtodbl(gel(S,2)) );
  }

  B = prec2nbits_mul(prec, M_LN2);
  if (s > 0 && !t) /* positive real input */
  {
    beta = B + 0.61 + s*(log2PI - log(s));
    if (beta > 0)
    {
      p = (long)ceil(beta / 2.0);
      n = fabs(s + 2*p-1)/(2*M_PI);
    }
    else
    {
      p = 0;
      n = exp((B - M_LN2) / s);
    }
  }
  else if (s <= 0 || t < 0.01) /* s < 0 may occur if s ~ 0 */
  { /* TODO: the crude bounds below are generally valid. Optimize ? */
    double l,l2, la = 1.; /* heuristic */
    double rlog, ilog; dcxlog(s-1,t, &rlog,&ilog);
    l2 = (s - 0.5)*rlog - t*ilog; /* = Re( (S - 1/2) log (S-1) ) */
    l = (B - l2 + s*log2PI) / (2. * (1.+ log((double)la)));
    l2 = dabs(s, t)/2;
    if (l < l2) l = l2;
    p = (long) ceil(l); if (p < 2) p = 2;
    n = 1 + dabs(p+s/2.-.25, t/2) * la / M_PI;
  }
  else
  {
    double sn = dabs(s, t), L = log(sn/s);
    alpha = B - 0.39 + L + s*(log2PI - log(sn));
    beta = (alpha+s)/t - atan(s/t);
    p = 0;
    if (beta > 0)
    {
      beta = 1.0 - s + t * get_xinf(beta);
      if (beta > 0) p = (long)ceil(beta / 2.0);
    }
    else
      if (s < 1.0) p = 1;
    n = p? dabs(s + 2*p-1, t) / (2*M_PI) : exp((B-M_LN2+L) / s);
  }
  *pp = p;
  *pn = (long)ceil(n);
  if (*pp < 0 || *pn < 0) pari_err_OVERFLOW("zeta");
}

/* zeta(a*j+b), j=0..N-1, b>1, using sumalt. Johansonn's thesis, Algo 4.7.1 */
static GEN
veczetas(long a, long b, long N, long prec)
{
  const long n = ceil(2 + prec2nbits_mul(prec, M_LN2/1.7627));
  pari_sp av = avma;
  GEN c, d, z = zerovec(N);
  long j, k;
  c = d = int2n(2*n-1);
  for (k = n; k > 1; k--)
  {
    GEN u, t = divii(d, powuu(k,b));
    if (!odd(k)) t = negi(t);
    gel(z,1) = addii(gel(z,1), t);
    u = powuu(k,a);
    for (j = 1; j < N; j++)
    {
      t = divii(t,u); if (!signe(t)) break;
      gel(z,j+1) = addii(gel(z,j+1), t);
    }
    c = muluui(k,2*k-1,c);
    c = diviuuexact(c, 2*(n-k+1),n+k-1);
    d = addii(d,c);
    if (gc_needed(av,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"zetaBorwein, k = %ld", k);
      gerepileall(av, 3, &c,&d,&z);
    }
  }
  /* k = 1 */
  for (j = 1; j <= N; j++) gel(z,j) = addii(gel(z,j), d);
  d = addiu(d, 1);
  for (j = 0, k = b - 1; j < N; j++, k += a)
    gel(z,j+1) = rdivii(shifti(gel(z,j+1), k), subii(shifti(d,k), d), prec);
  return z;
}
/* zeta(a*j+b), j=0..N-1, b>1, using sumalt */
GEN
veczeta(GEN a, GEN b, long N, long prec)
{
  pari_sp av = avma;
  long n, j, k;
  GEN L, c, d, z = zerovec(N);
  if (typ(a) == t_INT && typ(b) == t_INT)
    return gerepilecopy(av, veczetas(itos(a),  itos(b), N, prec));
  n = ceil(2 + prec2nbits_mul(prec, M_LN2/1.7627));
  c = d = int2n(2*n-1);
  for (k = n; k; k--)
  {
    GEN u, t;
    L = logr_abs(utor(k, prec)); /* log(k) */
    t = gdiv(d, gexp(gmul(b, L), prec)); /* d / k^b */
    if (!odd(k)) t = gneg(t);
    gel(z,1) = gadd(gel(z,1), t);
    u = gexp(gmul(a, L), prec);
    for (j = 1; j < N; j++)
    {
      t = gdiv(t,u); if (gexpo(t) < 0) break;
      gel(z,j+1) = gadd(gel(z,j+1), t);
    }
    c = muluui(k,2*k-1,c);
    c = diviuuexact(c, 2*(n-k+1),n+k-1);
    d = addii(d,c);
    if (gc_needed(av,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"veczeta, k = %ld", k);
      gerepileall(av, 3, &c,&d,&z);
    }
  }
  L = mplog2(prec);
  for (j = 0; j < N; j++)
  {
    GEN u = gsubgs(gadd(b, gmulgu(a,j)), 1);
    GEN w = gexp(gmul(u, L), prec); /* 2^u */
    gel(z,j+1) = gdiv(gmul(gel(z,j+1), w), gmul(d,gsubgs(w,1)));
  }
  return gerepilecopy(av, z);
}

GEN
constzeta(long n, long prec)
{
  GEN o = zetazone, z;
  long l = o? lg(o): 0;
  pari_sp av;
  if (l > n)
  {
    long p = realprec(gel(o,1));
    if (p >= prec) return o;
  }
  n = maxss(n, l + 15);
  av = avma; z = veczetas(1, 2, n-1, prec);
  zetazone = gclone(vec_prepend(z, mpeuler(prec)));
  set_avma(av); guncloneNULL(o); return zetazone;
}

/* zeta(s) using sumalt, case h=0,N=1. Assume s > 1 */
static GEN
zetaBorwein(long s, long prec)
{
  pari_sp av = avma;
  const long n = ceil(2 + prec2nbits_mul(prec, M_LN2/1.7627));
  long k;
  GEN c, d, z = gen_0;
  c = d = int2n(2*n-1);
  for (k = n; k; k--)
  {
    GEN t = divii(d, powuu(k,s));
    z = odd(k)? addii(z,t): subii(z,t);
    c = muluui(k,2*k-1,c);
    c = diviuuexact(c, 2*(n-k+1),n+k-1);
    d = addii(d,c);
    if (gc_needed(av,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"zetaBorwein, k = %ld", k);
      gerepileall(av, 3, &c,&d,&z);
    }
  }
  return rdivii(shifti(z, s-1), subii(shifti(d,s-1), d), prec);
}

/* assume k != 1 */
GEN
szeta(long k, long prec)
{
  pari_sp av = avma;
  GEN z;

  if (!k) { z = real2n(-1, prec); setsigne(z,-1); return z; }
  if (k < 0)
  {
    if (!odd(k)) return gen_0;
    /* the one value such that k < 0 and 1 - k < 0, due to overflow */
    if ((ulong)k == (HIGHBIT | 1))
      pari_err_OVERFLOW("zeta [large negative argument]");
    k = 1-k;
    z = bernreal(k, prec); togglesign(z);
    return gerepileuptoleaf(av, divru(z, k));
  }
  /* k > 1 */
  if (k > prec2nbits(prec)+1) return real_1(prec);
  if (zetazone && realprec(gel(zetazone,1)) >= prec && lg(zetazone) > k)
    return rtor(gel(zetazone, k), prec);
  if (!odd(k))
  {
    GEN B;
    if (!bernzone) constbern(0);
    if (k < lg(bernzone))
      B = gel(bernzone, k>>1);
    else
    {
      if (bernbitprec(k) > prec2nbits(prec))
        return gerepileupto(av, invr(inv_szeta_euler(k, prec)));
      B = bernfrac(k);
    }
    /* B = B_k */
    z = gmul(powru(Pi2n(1, prec + EXTRAPREC64), k), B);
    z = divrr(z, mpfactr(k,prec));
    setsigne(z, 1); shiftr_inplace(z, -1);
  }
  else
  {
    double p = prec2nbits_mul(prec,0.393); /* bit / log_2(3+sqrt(8)) */
    p = log2(p * log(p));
    z = (p * k > prec2nbits(prec))? invr(inv_szeta_euler(k, prec))
                                  : zetaBorwein(k, prec);
  }
  return gerepileuptoleaf(av, z);
}

/* Ensure |1-s| >= 1/32 and (|s| <= 1/32 or real(s) >= 1/2) */
static int
zeta_funeq(GEN *ps)
{
  GEN s = *ps, u;
  if (typ(s) == t_REAL)
  {
    u = subsr(1, s);
    if (expo(u) >= -5
        && ((signe(s) > 0 && expo(s) >= -1) || expo(s) <= -5)) return 0;
  }
  else
  {
    GEN sig = gel(s,1);
    u = gsubsg(1, s);
    if (gexpo(u) >= -5
        && ((signe(sig) > 0 && expo(sig) >= -1) || gexpo(s) <= -5)) return 0;
  }
  *ps = u; return 1;
}
/* s0 a t_INT, t_REAL or t_COMPLEX.
 * If a t_INT, assume it's not a trivial case (i.e we have s0 > 1, odd) */
static GEN
czeta(GEN s0, long prec)
{
  GEN ms, s, u, y, res, tes, sig, tau, invn2, ns, Ns, funeq_factor = NULL;
  long i, nn, lim, lim2;
  pari_sp av0 = avma, av, av2;
  pari_timer T;

  if (DEBUGLEVEL>2) timer_start(&T);
  s = trans_fix_arg(&prec,&s0,&sig,&tau,&av,&res);
  if (typ(s0) == t_INT) return gerepileupto(av0, gzeta(s0, prec));
  if (zeta_funeq(&s)) /* s -> 1-s */
  { /* Gamma(s) (2Pi)^-s 2 cos(Pi s/2) [new s] */
    GEN t = gmul(ggamma(s,prec), pow2Pis(gsubgs(s0,1), prec));
    sig = real_i(s);
    funeq_factor = gmul2n(gmul(t, gsin(gmul(Pi2n(-1,prec),s0), prec)), 1);
  }
  if (gcmpgs(sig, prec2nbits(prec) + 1) > 0) { /* zeta(s) = 1 */
    if (!funeq_factor) { set_avma(av0); return real_1(prec); }
    return gerepileupto(av0, funeq_factor);
  }
  optim_zeta(s, prec, &lim, &nn);
  if (DEBUGLEVEL>2) err_printf("lim, nn: [%ld, %ld]\n", lim, nn);

  ms = gneg(s);
  if (umuluu_le(nn, prec, 10000000))
  {
    incrprec(prec); /* one extra word of precision */
    Ns = vecpowug(nn, ms, prec);
    ns = gel(Ns,nn); setlg(Ns, nn);
    y = gadd(gmul2n(ns, -1), RgV_sum(Ns));
  }
  else
  {
    Ns = dirpowerssum(nn, ms, 0, prec);
    incrprec(prec); /* one extra word of precision */
    ns = gpow(utor(nn, prec), ms, prec);
    y = gsub(Ns, gmul2n(ns, -1));
  }
  if (DEBUGLEVEL>2) timer_printf(&T,"sum from 1 to N");
  constbern(lim);
  if (DEBUGLEVEL>2) timer_start(&T);
  invn2 = divri(real_1(prec), sqru(nn)); lim2 = lim<<1;
  tes = bernfrac(lim2);
  {
    GEN s1, s2, s3, s4, s5;
    s2 = gmul(s, gsubgs(s,1));
    s3 = gmul2n(invn2,3);
    av2 = avma;
    s1 = gsubgs(gmul2n(s,1), 1);
    s4 = gmul(invn2, gmul2n(gaddsg(4*lim-2,s1),1));
    s5 = gmul(invn2, gadd(s2, gmulsg(lim2, gaddgs(s1, lim2))));
    for (i = lim2-2; i>=2; i -= 2)
    {
      s5 = gsub(s5, s4);
      s4 = gsub(s4, s3);
      tes = gadd(bernfrac(i), gdivgunextu(gmul(s5,tes), i+1));
      if (gc_needed(av2,3))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"czeta i = %ld", i);
        gerepileall(av2,3, &tes,&s5,&s4);
      }
    }
    u = gmul(gmul(tes,invn2), gmul2n(s2, -1));
    tes = gmulsg(nn, gaddsg(1, u));
  }
  if (DEBUGLEVEL>2) timer_printf(&T,"Bernoulli sum");
  /* y += tes n^(-s) / (s-1) */
  y = gadd(y, gmul(tes, gdiv(ns, gsubgs(s,1))));
  if (funeq_factor) y = gmul(y, funeq_factor);
  set_avma(av); return affc_fixlg(y,res);
}
/* v a t_VEC/t_COL */
int
RgV_is_arithprog(GEN v, GEN *a, GEN *b)
{
  pari_sp av = avma, av2;
  long i, n = lg(v)-1;
  if (n == 0) { *a = *b = gen_0; return 1; }
  *a = gel(v,1);
  if (n == 1) { * b = gen_0; return 1; }
  *b = gsub(gel(v,2), *a); av2 = avma;
  for (i = 2; i < n; i++)
    if (!gequal(*b, gsub(gel(v,i+1), gel(v,i)))) return gc_int(av,0);
  return gc_int(av2,1);
}

GEN
gzeta(GEN x, long prec)
{
  pari_sp av = avma;
  GEN y;
  if (gequal1(x)) pari_err_DOMAIN("zeta", "argument", "=", gen_1, x);
  switch(typ(x))
  {
    case t_INT:
      if (is_bigint(x))
      {
        if (signe(x) > 0) return real_1(prec);
        if (mod2(x) == 0) return real_0(prec);
        pari_err_OVERFLOW("zeta [large negative argument]");
      }
      return szeta(itos(x),prec);
    case t_REAL: case t_COMPLEX: return czeta(x,prec);
    case t_PADIC: return Qp_zeta(x);
    case t_VEC: case t_COL:
    {
      GEN a, b;
      long n = lg(x) - 1;
      if (n > 1 && RgV_is_arithprog(x, &a, &b))
      {
        if (!is_real_t(typ(a)) || !is_real_t(typ(b)) || gcmpgs(a, 1) <= 0)
        { set_avma(av); break; }
        a = veczeta(b, a, n, prec);
        settyp(a, typ(x)); return a;
      }
    }
    default:
      if (!(y = toser_i(x))) break;
      if (gequal1(y))
        pari_err_DOMAIN("zeta", "argument", "=", gen_1, y);
      return gerepileupto(av, lfun(gen_1,y,prec2nbits(prec)));
  }
  return trans_eval("zeta",gzeta,x,prec);
}

/***********************************************************************/
/**                                                                   **/
/**                    FONCTIONS POLYLOGARITHME                       **/
/**                                                                   **/
/***********************************************************************/

/* smallish k such that bernbitprec(K) > bit + Kdz, K = 2k+4 */
static long
get_k(double dz, long bit)
{
  long a, b;
  for (b = 128;; b <<= 1)
    if (bernbitprec(b) > bit + b*dz) break;
  if (b == 128) return 128;
  a = b >> 1;
  while (b - a > 64)
  {
    long c = (a+b) >> 1;
    if (bernbitprec(c) > bit + c*dz) b = c; else a = c;
  }
  return b >> 1;
}

/* m >= 2. Validity domain |log x| < 2*Pi, contains log |x| < 5.44,
 * Li_m(x = e^z) = sum_{n >= 0} zeta(m-n) z^n / n!
 *    with zeta(1) := H_{m-1} - log(-z) */
static GEN
cxpolylog(long m, GEN x, long prec)
{
  long li, n, k, ksmall, real;
  GEN vz, z, Z, h, q, s, S;
  pari_sp av;
  double dz;
  pari_timer T;

  if (gequal1(x)) return szeta(m,prec);
  /* x real <= 1 ==> Li_m(x) real */
  real = (typ(x) == t_REAL && (expo(x) < 0 || signe(x) <= 0));

  vz = constzeta(m, prec);
  z = glog(x,prec);
  /* n = 0 */
  q = gen_1; s = gel(vz, m);
  for (n=1; n < m-1; n++)
  {
    q = gdivgu(gmul(q,z),n);
    s = gadd(s, gmul(gel(vz,m-n), real? real_i(q): q));
  }
  /* n = m-1 */
    q = gdivgu(gmul(q,z),n); /* multiply by "zeta(1)" */
    h = gmul(q, gsub(harmonic(m-1), glog(gneg_i(z),prec)));
    s = gadd(s, real? real_i(h): h);
  /* n = m */
    q = gdivgu(gmul(q,z),m);
    s = gadd(s, gdivgs(real? real_i(q): q, -2)); /* zeta(0) = -1/2 */
  /* n = m+1 */
    q = gdivgu(gmul(q,z),m+1); /* = z^(m+1) / (m+1)! */
    s = gadd(s, gdivgs(real? real_i(q): q, -12)); /* zeta(-1) = -1/12 */

  li = -(prec2nbits(prec)+1);
  if (DEBUGLEVEL) timer_start(&T);
  dz = dbllog2(z) - log2PI; /*  ~ log2(|z|/2Pi) */
  /* sum_{k >= 1} zeta(-1-2k) * z^(2k+m+1) / (2k+m+1)!
   * = 2 z^(m-1) sum_{k >= 1} zeta(2k+2) * Z^(k+1) / (2k+2)..(2k+1+m)), where
   * Z = -(z/2Pi)^2. Stop at 2k = (li - (m-1)*Lz - m) / dz, Lz = log2 |z| */
  /* We cut the sum in two: small values of k first */
  Z = gsqr(z); av = avma;
  ksmall = get_k(dz, prec2nbits(prec));
  constbern(ksmall);
  for(k = 1; k < ksmall; k++)
  {
    GEN t = q = gdivgunextu(gmul(q,Z), 2*k+m); /* z^(2k+m+1)/(2k+m+1)! */
    if (real) t = real_i(t);
    t = gmul(t, gdivgu(bernfrac(2*k+2), 2*k+2)); /* - t * zeta(1-(2k+2)) */
    s = gsub(s, t);
    if (gexpo(t)  < li) return s;
    /* large values ? */
    if ((k & 0x1ff) == 0) gerepileall(av, 2, &s, &q);
  }
  if (DEBUGLEVEL>2) timer_printf(&T, "polylog: small k <= %ld", k);
  Z = gneg(gsqr(gdiv(z, Pi2n(1,prec))));
  q = gmul(gpowgs(z, m-1), gpowgs(Z, k+1)); /* Z^(k+1) * z^(m-1) */
  S = gen_0; av = avma;
  for(;; k++)
  {
    GEN t = q;
    long b;
    if (real) t = real_i(t);
    b = prec + gexpo(t) / BITS_IN_LONG; /* decrease accuracy */
    if (b == 2) break;
    /* t * zeta(2k+2) / (2k+2)..(2k+1+m) */
    t = gdiv(t, mulri(inv_szeta_euler(2*k+2, b),
                      mulu_interval(2*k+2, 2*k+1+m)));
    S = gadd(S, t); if (gexpo(t)  < li) break;
    q = gmul(q, Z);
    if ((k & 0x1ff) == 0) gerepileall(av, 2, &S, &q);
  }
  if (DEBUGLEVEL>2) timer_printf(&T, "polylog: large k <= %ld", k);
  return gadd(s, gmul2n(S,1));
}

static GEN
Li1(GEN x, long prec) { return gneg(glog(gsubsg(1, x), prec)); }

static GEN
polylog(long m, GEN x, long prec)
{
  long l, e, i, G, sx;
  pari_sp av, av1;
  GEN X, Xn, z, p1, p2, y, res;

  if (m < 0) pari_err_DOMAIN("polylog", "index", "<", gen_0, stoi(m));
  if (!m) return mkfrac(gen_m1,gen_2);
  if (gequal0(x)) return gcopy(x);
  if (m==1) { av = avma; return gerepileupto(av, Li1(x, prec)); }

  l = precision(x);
  if (!l) l = prec; else prec = l;
  res = cgetc(l); av = avma;
  x = gtofp(x, l+EXTRAPREC64);
  e = gexpo(gnorm(x));
  if (!e || e == -1) {
    y = cxpolylog(m,x,prec);
    set_avma(av); return affc_fixlg(y, res);
  }
  X = (e > 0)? ginv(x): x;
  G = -prec2nbits(l);
  av1 = avma;
  y = Xn = X;
  for (i=2; ; i++)
  {
    Xn = gmul(X,Xn); p2 = gdiv(Xn,powuu(i,m));
    y = gadd(y,p2);
    if (gexpo(p2) <= G) break;

    if (gc_needed(av1,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"polylog");
      gerepileall(av1,2, &y, &Xn);
    }
  }
  if (e < 0) { set_avma(av); return affc_fixlg(y, res); }

  sx = gsigne(imag_i(x));
  if (!sx)
  {
    if (m&1) sx = gsigne(gsub(gen_1, real_i(x)));
    else     sx = - gsigne(real_i(x));
  }
  z = divri(mppi(l), mpfact(m-1)); setsigne(z, sx);
  z = mkcomplex(gen_0, z);

  if (m == 2)
  { /* same formula as below, written more efficiently */
    y = gneg_i(y);
    if (typ(x) == t_REAL && signe(x) < 0)
      p1 = logr_abs(x);
    else
      p1 = gsub(glog(x,l), z);
    p1 = gmul2n(gsqr(p1), -1); /* = (log(-x))^2 / 2 */

    p1 = gadd(p1, divru(sqrr(mppi(l)), 6));
    p1 = gneg_i(p1);
  }
  else
  {
    GEN logx = glog(x,l), logx2 = gsqr(logx), vz = constzeta(m, l);
    p1 = mkfrac(gen_m1,gen_2);
    for (i = m-2; i >= 0; i -= 2)
      p1 = gadd(gel(vz, m-i), gmul(p1, gdivgunextu(logx2, i+1)));
    if (m&1) p1 = gmul(logx,p1); else y = gneg_i(y);
    p1 = gadd(gmul2n(p1,1), gmul(z,gpowgs(logx,m-1)));
    if (typ(x) == t_REAL && signe(x) < 0) p1 = real_i(p1);
  }
  y = gadd(y,p1);
  set_avma(av); return affc_fixlg(y, res);
}
static GEN
RIpolylog(long m, GEN x, long real, long prec)
{
  GEN y = polylog(m, x, prec);
  return real? real_i(y): imag_i(y);
}
GEN
dilog(GEN x, long prec) { return gpolylog(2, x, prec); }

/* x a floating point number, t_REAL or t_COMPLEX of t_REAL */
static GEN
logabs(GEN x)
{
  GEN y;
  if (typ(x) == t_COMPLEX)
  {
    y = logr_abs( cxnorm(x) );
    shiftr_inplace(y, -1);
  } else
    y = logr_abs(x);
  return y;
}

static GEN
polylogD(long m, GEN x, long flag, long prec)
{
  long fl = 0, k, l, m2;
  pari_sp av;
  GEN p1, p2, y;

  if (gequal0(x)) return gcopy(x);
  m2 = m&1;
  if (gequal1(x) && m>=2) return m2? szeta(m,prec): gen_0;
  av = avma; l = precision(x);
  if (!l) { l = prec; x = gtofp(x,l); }
  p1 = logabs(x);
  if (signe(p1) > 0) { x = ginv(x); fl = !m2; } else setabssign(p1);
  /* |x| <= 1, p1 = - log|x| >= 0 */
  p2 = gen_1;
  y = RIpolylog(m, x, m2, l);
  for (k = 1; k < m; k++)
  {
    GEN t = RIpolylog(m-k, x, m2, l);
    p2 = gdivgu(gmul(p2,p1), k); /* (-log|x|)^k / k! */
    y = gadd(y, gmul(p2, t));
  }
  if (m2)
  {
    p1 = flag? gdivgs(p1, -2*m): gdivgs(logabs(gsubsg(1,x)), m);
    y = gadd(y, gmul(p2, p1));
  }
  if (fl) y = gneg(y);
  return gerepileupto(av, y);
}

static GEN
polylogP(long m, GEN x, long prec)
{
  long fl = 0, k, l, m2;
  pari_sp av;
  GEN p1,y;

  if (gequal0(x)) return gcopy(x);
  m2 = m&1;
  if (gequal1(x) && m>=2) return m2? szeta(m,prec): gen_0;
  av = avma; l = precision(x);
  if (!l) { l = prec; x = gtofp(x,l); }
  p1 = logabs(x);
  if (signe(p1) > 0) { x = ginv(x); fl = !m2; setsigne(p1, -1); }
  /* |x| <= 1 */
  y = RIpolylog(m, x, m2, l);
  if (m==1)
  {
    shiftr_inplace(p1, -1); /* log |x| / 2 */
    y = gadd(y, p1);
  }
  else
  { /* m >= 2, \sum_{0 <= k <= m} 2^k B_k/k! (log |x|)^k Li_{m-k}(x),
       with Li_0(x) := -1/2 */
    GEN u, t = RIpolylog(m-1, x, m2, l);
    u = gneg_i(p1); /* u = 2 B1 log |x| */
    y = gadd(y, gmul(u, t));
    if (m > 2)
    {
      GEN p2;
      shiftr_inplace(p1, 1); /* 2log|x| <= 0 */
      constbern(m>>1);
      p1 = sqrr(p1);
      p2 = shiftr(p1,-1);
      for (k = 2; k < m; k += 2)
      {
        if (k > 2) p2 = gdivgunextu(gmul(p2,p1),k-1); /* 2^k/k! log^k |x|*/
        t = RIpolylog(m-k, x, m2, l);
        u = gmul(p2, bernfrac(k));
        y = gadd(y, gmul(u, t));
      }
    }
  }
  if (fl) y = gneg(y);
  return gerepileupto(av, y);
}

static GEN
gpolylog_i(void *E, GEN x, long prec)
{
  pari_sp av = avma;
  long i, n, v, m = (long)E;
  GEN a, y;

  if (m <= 0)
  {
    a = gmul(x, poleval(eulerianpol(-m, 0), x));
    return gerepileupto(av, gdiv(a, gpowgs(gsubsg(1, x), 1-m)));
  }
  switch(typ(x))
  {
    case t_REAL: case t_COMPLEX: return polylog(m,x,prec);
    case t_INTMOD: case t_PADIC: pari_err_IMPL( "padic polylogarithm");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (!m) { set_avma(av); return mkfrac(gen_m1,gen_2); }
      if (m==1) return gerepileupto(av, Li1(y, prec));
      if (gequal0(y)) return gerepilecopy(av, y);
      v = valser(y);
      if (v < 0) pari_err_DOMAIN("polylog","valuation", "<", gen_0, x);
      if (v > 0) {
        n = (lg(y)-3 + v) / v;
        a = zeroser(varn(y), lg(y)-2);
        for (i=n; i>=1; i--)
          a = gmul(y, gadd(a, powis(utoipos(i),-m)));
      } else { /* v == 0 */
        long vy = varn(y);
        GEN a0 = polcoef_i(y, 0, -1), t = gdiv(derivser(y), y);
        a = Li1(y, prec);
        for (i=2; i<=m; i++)
          a = gadd(gpolylog(i, a0, prec), integ(gmul(t, a), vy));
      }
      return gerepileupto(av, a);
  }
  return trans_evalgen("polylog", E, gpolylog_i, x, prec);
}
GEN
gpolylog(long m, GEN x, long prec) { return gpolylog_i((void*)m, x, prec); }

GEN
polylog0(long m, GEN x, long flag, long prec)
{
  switch(flag)
  {
    case 0: return gpolylog(m,x,prec);
    case 1: return polylogD(m,x,0,prec);
    case 2: return polylogD(m,x,1,prec);
    case 3: return polylogP(m,x,prec);
    default: pari_err_FLAG("polylog");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

GEN
upper_to_cx(GEN x, long *prec)
{
  long tx = typ(x), l;
  if (tx == t_QUAD) { x = quadtofp(x, *prec); tx = typ(x); }
  switch(tx)
  {
    case t_COMPLEX:
      if (gsigne(gel(x,2)) > 0) break; /*fall through*/
    case t_REAL: case t_INT: case t_FRAC:
      pari_err_DOMAIN("modular function", "Im(argument)", "<=", gen_0, x);
    default:
      pari_err_TYPE("modular function", x);
  }
  l = precision(x); if (l) *prec = l;
  return x;
}

/* sqrt(3)/2 */
static GEN
sqrt32(long prec) { GEN z = sqrtr_abs(utor(3,prec)); setexpo(z, -1); return z; }
/* exp(i x), x = k pi/12 */
static GEN
e12(ulong k, long prec)
{
  int s, sPi, sPiov2;
  GEN z, t;
  k %= 24;
  if (!k) return gen_1;
  if (k == 12) return gen_m1;
  if (k >12) { s = 1; k = 24 - k; } else s = 0; /* x -> 2pi - x */
  if (k > 6) { sPi = 1; k = 12 - k; } else sPi = 0; /* x -> pi  - x */
  if (k > 3) { sPiov2 = 1; k = 6 - k; } else sPiov2 = 0; /* x -> pi/2 - x */
  z = cgetg(3, t_COMPLEX);
  switch(k)
  {
    case 0: gel(z,1) = icopy(gen_1); gel(z,2) = gen_0; break;
    case 1: t = gmul2n(addrs(sqrt32(prec), 1), -1);
      gel(z,1) = sqrtr(t);
      gel(z,2) = gmul2n(invr(gel(z,1)), -2); break;

    case 2: gel(z,1) = sqrt32(prec);
            gel(z,2) = real2n(-1, prec); break;

    case 3: gel(z,1) = sqrtr_abs(real2n(-1,prec));
            gel(z,2) = rcopy(gel(z,1)); break;
  }
  if (sPiov2) swap(gel(z,1), gel(z,2));
  if (sPi) togglesign(gel(z,1));
  if (s)   togglesign(gel(z,2));
  return z;
}
/* z a t_FRAC */
static GEN
expIPifrac(GEN z, long prec)
{
  GEN n = gel(z,1), d = gel(z,2);
  ulong r, q = uabsdivui_rem(12, d, &r);
  if (!r) return e12(q * umodiu(n, 24), prec); /* 12 | d */
  n = centermodii(n, shifti(d,1), d);
  return expIr(divri(mulri(mppi(prec), n), d));
}
/* exp(i Pi z), z a t_INT or t_FRAC */
static GEN
expIPiQ(GEN z, long prec)
{
  if (typ(z) == t_INT) return mpodd(z)? gen_m1: gen_1;
  return expIPifrac(z, prec);
}

/* convert power of 2 t_REAL to rational */
static GEN
real2nQ(GEN x)
{
  long e = expo(x);
  GEN z;
  if (e < 0)
    z = mkfrac(signe(x) < 0? gen_m1: gen_1, int2n(-e));
  else
  {
    z = int2n(e);
    if (signe(x) < 0) togglesign_safe(&z);
  }
  return z;
}
/* x a real number */
GEN
expIPiR(GEN x, long prec)
{
  if (typ(x) == t_REAL && absrnz_equal2n(x)) x = real2nQ(x);
  switch(typ(x))
  {
    case t_INT:  return mpodd(x)? gen_m1: gen_1;
    case t_FRAC: return expIPifrac(x, prec);
  }
  return expIr(mulrr(mppi(prec), x));
}
/* z a t_COMPLEX */
GEN
expIPiC(GEN z, long prec)
{
  GEN pi, r, x, y;
  if (typ(z) != t_COMPLEX) return expIPiR(z, prec);
  x = gel(z,1);
  y = gel(z,2); if (gequal0(y)) return expIPiR(x, prec);
  pi = mppi(prec);
  r = gmul(pi, y); togglesign(r); r = mpexp(r); /* exp(-pi y) */
  if (typ(x) == t_REAL && absrnz_equal2n(x)) x = real2nQ(x);
  switch(typ(x))
  {
    case t_INT: if (mpodd(x)) togglesign(r);
                return r;
    case t_FRAC: return gmul(r, expIPifrac(x, prec));
  }
  return gmul(r, expIr(mulrr(pi, x)));
}
/* exp(I x y), more efficient for x in R, y pure imaginary */
GEN
expIxy(GEN x, GEN y, long prec) { return gexp(gmul(x, mulcxI(y)), prec); }

static GEN
qq(GEN x, long prec)
{
  long tx = typ(x);
  GEN y;

  if (is_scalar_t(tx))
  {
    if (tx == t_PADIC) return x;
    x = upper_to_cx(x, &prec);
    return cxtoreal(expIPiC(gmul2n(x,1), prec)); /* e(x) */
  }
  if (! ( y = toser_i(x)) ) pari_err_TYPE("modular function", x);
  return y;
}

/* return (y * X^d) + x. Assume d > 0, x != 0, valser(x) = 0 */
static GEN
ser_addmulXn(GEN y, GEN x, long d)
{
  long i, lx, ly, l = valser(y) + d; /* > 0 */
  GEN z;

  lx = lg(x);
  ly = lg(y) + l; if (lx < ly) ly = lx;
  if (l > lx-2) return gcopy(x);
  z = cgetg(ly,t_SER);
  for (i=2; i<=l+1; i++) gel(z,i) = gel(x,i);
  for (   ; i < ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i-l));
  z[1] = x[1]; return z;
}

/* q a t_POL s.t. q(0) != 0, v > 0, Q = x^v*q; return \prod_i (1-Q^i) */
static GEN
RgXn_eta(GEN q, long v, long lim)
{
  pari_sp av = avma;
  GEN qn, ps, y;
  ulong vps, vqn, n;

  if (!degpol(q) && isint1(gel(q,2))) return eta_ZXn(v, lim+v);
  y = qn = ps = pol_1(0);
  vps = vqn = 0;
  for(n = 0;; n++)
  { /* qn = q^n,  ps = (-1)^n q^(n(3n+1)/2),
     * vps, vqn valuation of ps, qn HERE */
    pari_sp av2 = avma;
    ulong vt = vps + 2*vqn + v; /* valuation of t at END of loop body */
    long k1, k2;
    GEN t;
    vqn += v; vps = vt + vqn; /* valuation of qn, ps at END of body */
    k1 = lim + v - vt + 1;
    k2 = k1 - vqn; /* = lim + v - vps + 1 */
    if (k1 <= 0) break;
    t = RgX_mul(q, RgX_sqr(qn));
    t = RgXn_red_shallow(t, k1);
    t = RgX_mul(ps,t);
    t = RgXn_red_shallow(t, k1);
    t = RgX_neg(t); /* t = (-1)^(n+1) q^(n(3n+1)/2 + 2n+1) */
    t = gerepileupto(av2, t);
    y = RgX_addmulXn_shallow(t, y, vt);
    if (k2 <= 0) break;

    qn = RgX_mul(qn,q);
    ps = RgX_mul(t,qn);
    ps = RgXn_red_shallow(ps, k2);
    y = RgX_addmulXn_shallow(ps, y, vps);

    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"eta, n = %ld", n);
      gerepileall(av, 3, &y, &qn, &ps);
    }
  }
  return y;
}

static GEN
inteta(GEN q)
{
  long tx = typ(q);
  GEN ps, qn, y;

  y = gen_1; qn = gen_1; ps = gen_1;
  if (tx==t_PADIC)
  {
    if (valp(q) <= 0) pari_err_DOMAIN("eta", "v_p(q)", "<=",gen_0,q);
    for(;;)
    {
      GEN t = gneg_i(gmul(ps,gmul(q,gsqr(qn))));
      y = gadd(y,t); qn = gmul(qn,q); ps = gmul(t,qn);
      t = y;
      y = gadd(y,ps); if (gequal(t,y)) return y;
    }
  }

  if (tx == t_SER)
  {
    ulong vps, vqn;
    long l = lg(q), v, n;
    pari_sp av;

    v = valser(q); /* handle valuation separately to avoid overflow */
    if (v <= 0) pari_err_DOMAIN("eta", "v_p(q)", "<=",gen_0,q);
    y = ser2pol_i(q, l); /* t_SER inefficient when input has low degree */
    n = degpol(y);
    if (n <= (l>>2))
    {
      GEN z = RgXn_eta(y, v, l-2);
      setvarn(z, varn(y)); return RgX_to_ser(z, l+v);
    }
    q = leafcopy(q); av = avma;
    setvalser(q, 0);
    y = scalarser(gen_1, varn(q), l+v);
    vps = vqn = 0;
    for(n = 0;; n++)
    { /* qn = q^n,  ps = (-1)^n q^(n(3n+1)/2) */
      ulong vt = vps + 2*vqn + v;
      long k;
      GEN t;
      t = gneg_i(gmul(ps,gmul(q,gsqr(qn))));
      /* t = (-1)^(n+1) q^(n(3n+1)/2 + 2n+1) */
      y = ser_addmulXn(t, y, vt);
      vqn += v; vps = vt + vqn;
      k = l+v - vps; if (k <= 2) return y;

      qn = gmul(qn,q); ps = gmul(t,qn);
      y = ser_addmulXn(ps, y, vps);
      setlg(q, k);
      setlg(qn, k);
      setlg(ps, k);
      if (gc_needed(av,3))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"eta");
        gerepileall(av, 3, &y, &qn, &ps);
      }
    }
  }
  {
    long l = -prec2nbits(precision(q));
    pari_sp av = avma;

    for(;;)
    {
      GEN t = gneg_i(gmul(ps,gmul(q,gsqr(qn))));
      /* qn = q^n
       * ps = (-1)^n q^(n(3n+1)/2)
       * t = (-1)^(n+1) q^(n(3n+1)/2 + 2n+1) */
      y = gadd(y,t); qn = gmul(qn,q); ps = gmul(t,qn);
      y = gadd(y,ps);
      if (gexpo(ps)-gexpo(y) < l) return y;
      if (gc_needed(av,3))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"eta");
        gerepileall(av, 3, &y, &qn, &ps);
      }
    }
  }
}

GEN
eta(GEN x, long prec)
{
  pari_sp av = avma;
  GEN z = inteta( qq(x,prec) );
  if (typ(z) == t_SER) return gerepilecopy(av, z);
  return gerepileupto(av, z);
}

/* s(h,k) = sum(n = 1, k-1, (n/k)*(frac(h*n/k) - 1/2))
 * Knuth's algorithm. h integer, k integer > 0, (h,k) = 1 */
GEN
sumdedekind_coprime(GEN h, GEN k)
{
  pari_sp av = avma;
  GEN s2, s1, p, pp;
  long s;
  if (lgefint(k) == 3 && uel(k,2) <= (2*(ulong)LONG_MAX) / 3)
  {
    ulong kk = k[2], hh = umodiu(h, kk);
    long s1, s2;
    GEN v;
    if (signe(k) < 0) { k = negi(k); hh = Fl_neg(hh, kk); }
    v = u_sumdedekind_coprime(hh, kk);
    s1 = v[1]; s2 = v[2];
    return gerepileupto(av, gdiv(addis(mulis(k,s1), s2), muluu(12, kk)));
  }
  s = 1;
  s1 = gen_0; p = gen_1; pp = gen_0;
  s2 = h = modii(h, k);
  while (signe(h)) {
    GEN r, nexth, a = dvmdii(k, h, &nexth);
    if (is_pm1(h)) s2 = s == 1? addii(s2, p): subii(s2, p);
    s1 = s == 1? addii(s1, a): subii(s1, a);
    s = -s;
    k = h; h = nexth;
    r = addii(mulii(a,p), pp); pp = p; p = r;
  }
  /* at this point p = original k */
  if (s == -1) s1 = subiu(s1, 3);
  return gerepileupto(av, gdiv(addii(mulii(p,s1), s2), muliu(p,12)));
}
/* as above, for ulong arguments.
 * k integer > 0, 0 <= h < k, (h,k) = 1. Returns [s1,s2] such that
 * s(h,k) = (s2 + k s1) / (12k). Requires max(h + k/2, k) < LONG_MAX
 * to avoid overflow, in particular k <= LONG_MAX * 2/3 is fine */
GEN
u_sumdedekind_coprime(long h, long k)
{
  long s = 1, s1 = 0, s2 = h, p = 1, pp = 0;
  while (h) {
    long r, nexth = k % h, a = k / h; /* a >= 1, a >= 2 if h = 1 */
    if (h == 1) s2 += p * s; /* occurs exactly once, last step */
    s1 += a * s;
    s = -s;
    k = h; h = nexth;
    r = a*p + pp; pp = p; p = r; /* p >= pp >= 0 */
  }
  /* in the above computation, p increases from 1 to original k,
   * -k/2 <= s2 <= h + k/2, and |s1| <= k */
  if (s < 0) s1 -= 3; /* |s1| <= k+3 ? */
  /* But in fact, |s2 + p s1| <= k^2 + 1/2 - 3k; if (s < 0), we have
   * |s2| <= k/2 and it follows that |s1| < k here as well */
  /* p = k; s(h,k) = (s2 + p s1)/12p. */
  return mkvecsmall2(s1, s2);
}
GEN
sumdedekind(GEN h, GEN k)
{
  pari_sp av = avma;
  GEN d;
  if (typ(h) != t_INT) pari_err_TYPE("sumdedekind",h);
  if (typ(k) != t_INT) pari_err_TYPE("sumdedekind",k);
  d = gcdii(h,k);
  if (!is_pm1(d))
  {
    h = diviiexact(h, d);
    k = diviiexact(k, d);
  }
  return gerepileupto(av, sumdedekind_coprime(h,k));
}

/* eta(x); assume Im x >> 0 (e.g. x in SL2's standard fundamental domain) */
static GEN
eta_reduced(GEN x, long prec)
{
  GEN z = expIPiC(gdivgu(x, 12), prec); /* e(x/24) */
  if (24 * gexpo(z) >= -prec2nbits(prec))
    z = gmul(z, inteta( gpowgs(z,24) ));
  return z;
}

/* x = U.z (flag = 1), or x = U^(-1).z (flag = 0)
 * Return [s,t] such that eta(z) = eta(x) * sqrt(s) * exp(I Pi t) */
static GEN
eta_correction(GEN x, GEN U, long flag)
{
  GEN a,b,c,d, s,t;
  long sc;
  a = gcoeff(U,1,1);
  b = gcoeff(U,1,2);
  c = gcoeff(U,2,1);
  d = gcoeff(U,2,2);
  /* replace U by U^(-1) */
  if (flag) {
    swap(a,d);
    togglesign_safe(&b);
    togglesign_safe(&c);
  }
  sc = signe(c);
  if (!sc) {
    if (signe(d) < 0) togglesign_safe(&b);
    s = gen_1;
    t = uutoQ(umodiu(b, 24), 12);
  } else {
    if (sc < 0) {
      togglesign_safe(&a);
      togglesign_safe(&b);
      togglesign_safe(&c);
      togglesign_safe(&d);
    } /* now c > 0 */
    s = mulcxmI(gadd(gmul(c,x), d));
    t = gadd(gdiv(addii(a,d),muliu(c,12)), sumdedekind_coprime(negi(d),c));
    /* correction : exp(I Pi (((a+d)/12c) + s(-d,c)) ) sqrt(-i(cx+d))  */
  }
  return mkvec2(s, t);
}

/* returns the true value of eta(x) for Im(x) > 0, using reduction to
 * standard fundamental domain */
GEN
trueeta(GEN x, long prec)
{
  pari_sp av = avma;
  GEN U, st, s, t;

  if (!is_scalar_t(typ(x))) pari_err_TYPE("trueeta",x);
  x = upper_to_cx(x, &prec);
  x = cxredsl2(x, &U);
  st = eta_correction(x, U, 1);
  x = eta_reduced(x, prec);
  s = gel(st, 1);
  t = gel(st, 2);
  x = gmul(x, expIPiQ(t, prec));
  if (s != gen_1) x = gmul(x, gsqrt(s, prec));
  return gerepileupto(av, x);
}

GEN
eta0(GEN x, long flag,long prec)
{ return flag? trueeta(x,prec): eta(x,prec); }

/* eta(q) = 1 + \sum_{n>0} (-1)^n * (q^(n(3n-1)/2) + q^(n(3n+1)/2)) */
static GEN
ser_eta(long prec)
{
  GEN e = cgetg(prec+2, t_SER), ed = e+2;
  long n, j;
  e[1] = evalsigne(1)|_evalvalser(0)|evalvarn(0);
  gel(ed,0) = gen_1;
  for (n = 1; n < prec; n++) gel(ed,n) = gen_0;
  for (n = 1, j = 0; n < prec; n++)
  {
    GEN s;
    j += 3*n-2; /* = n*(3*n-1) / 2 */;
    if (j >= prec) break;
    s = odd(n)? gen_m1: gen_1;
    gel(ed, j) = s;
    if (j+n >= prec) break;
    gel(ed, j+n) = s;
  }
  return e;
}

static GEN
coeffEu(GEN fa)
{
  pari_sp av = avma;
  return gerepileuptoint(av, mului(65520, usumdivk_fact(fa,11)));
}
/* E12 = 1 + q*E/691 */
static GEN
ser_E(long prec)
{
  GEN e = cgetg(prec+2, t_SER), ed = e+2;
  GEN F = vecfactoru_i(2, prec); /* F[n] = factoru(n+1) */
  long n;
  e[1] = evalsigne(1)|_evalvalser(0)|evalvarn(0);
  gel(ed,0) = utoipos(65520);
  for (n = 1; n < prec; n++) gel(ed,n) = coeffEu(gel(F,n));
  return e;
}
/* j = E12/Delta + 432000/691, E12 = 1 + q*E/691 */
static GEN
ser_j2(long prec, long v)
{
  pari_sp av = avma;
  GEN iD = gpowgs(ginv(ser_eta(prec)), 24); /* q/Delta */
  GEN J = gmul(ser_E(prec), iD);
  setvalser(iD,-1); /* now 1/Delta */
  J = gadd(gdivgu(J, 691), iD);
  J = gerepileupto(av, J);
  if (prec > 1) gel(J,3) = utoipos(744);
  setvarn(J,v); return J;
}

/* j(q) = \sum_{n >= -1} c(n)q^n,
 * \sum_{n = -1}^{N-1} c(n) (-10n \sigma_3(N-n) + 21 \sigma_5(N-n))
 * = c(N) (N+1)/24 */
static GEN
ser_j(long prec, long v)
{
  GEN j, J, S3, S5, F;
  long i, n;
  if (prec > 64) return ser_j2(prec, v);
  S3 = cgetg(prec+1, t_VEC);
  S5 = cgetg(prec+1,t_VEC);
  F = vecfactoru_i(1, prec);
  for (n = 1; n <= prec; n++)
  {
    GEN fa = gel(F,n);
    gel(S3,n) = mului(10, usumdivk_fact(fa,3));
    gel(S5,n) = mului(21, usumdivk_fact(fa,5));
  }
  J = cgetg(prec+2, t_SER),
  J[1] = evalvarn(v)|evalsigne(1)|evalvalser(-1);
  j = J+3;
  gel(j,-1) = gen_1;
  gel(j,0) = utoipos(744);
  gel(j,1) = utoipos(196884);
  for (n = 2; n < prec; n++)
  {
    pari_sp av = avma;
    GEN c, s3 = gel(S3,n+1), s5 = gel(S5,n+1);
    c = addii(s3, s5);
    for (i = 0; i < n; i++)
    {
      s3 = gel(S3,n-i); s5 = gel(S5,n-i);
      c = addii(c, mulii(gel(j,i), subii(s5, mului(i,s3))));
    }
    gel(j,n) = gerepileuptoint(av, diviuexact(muliu(c,24), n+1));
  }
  return J;
}

GEN
jell(GEN x, long prec)
{
  long tx = typ(x);
  pari_sp av = avma;
  GEN q, h, U;

  if (!is_scalar_t(tx))
  {
    long v;
    if (gequalX(x)) return ser_j(precdl, varn(x));
    q = toser_i(x); if (!q) pari_err_TYPE("ellj",x);
    v = fetch_var_higher();
    h = ser_j(lg(q)-2, v);
    h = gsubst(h, v, q);
    delete_var(); return gerepileupto(av, h);
  }
  if (tx == t_PADIC)
  {
    GEN p2, p1 = gdiv(inteta(gsqr(x)), inteta(x));
    p1 = gmul2n(gsqr(p1),1);
    p1 = gmul(x,gpowgs(p1,12));
    p2 = gaddsg(768,gadd(gsqr(p1),gdivsg(4096,p1)));
    p1 = gmulsg(48,p1);
    return gerepileupto(av, gadd(p2,p1));
  }
  /* Let h = Delta(2x) / Delta(x), then j(x) = (1 + 256h)^3 / h */
  x = upper_to_cx(x, &prec);
  x = cxredsl2(x, &U); /* forget about Ua : j has weight 0 */
  { /* cf eta_reduced, raised to power 24
     * Compute
     *   t = (inteta(q(2x)) / inteta(q(x))) ^ 24;
     * then
     *   h = t * (q(2x) / q(x) = t * q(x);
     * but inteta(q) costly and useless if expo(q) << 1  => inteta(q) = 1.
     * log_2 ( exp(-2Pi Im tau) ) < -prec2nbits(prec)
     * <=> Im tau > prec2nbits(prec) * log(2) / 2Pi */
    long C = (long)prec2nbits_mul(prec, M_LN2/(2*M_PI));
    q = expIPiC(gmul2n(x,1), prec); /* e(x) */
    if (gcmpgs(gel(x,2), C) > 0) /* eta(q(x)) = 1 : no need to compute q(2x) */
      h = q;
    else
    {
      GEN t = gdiv(inteta(gsqr(q)), inteta(q));
      h = gmul(q, gpowgs(t, 24));
    }
  }
  /* real_1 important ! gaddgs(, 1) could increase the accuracy ! */
  return gerepileupto(av, gdiv(gpowgs(gadd(gmul2n(h,8), real_1(prec)), 3), h));
}

static GEN
to_form(GEN a, GEN w, GEN C)
{ return mkvec3(a, w, diviiexact(C, a)); }
static GEN
form_to_quad(GEN f, GEN sqrtD)
{
  long a = itos(gel(f,1)), a2 = a << 1;
  GEN b = gel(f,2);
  return mkcomplex(gdivgs(b, -a2), gdivgs(sqrtD, a2));
}
static GEN
eta_form(GEN f, GEN sqrtD, GEN *s_t, long prec)
{
  GEN U, t = form_to_quad(redimagsl2(f, &U), sqrtD);
  *s_t = eta_correction(t, U, 0);
  return eta_reduced(t, prec);
}

/* eta(t/p)eta(t/q) / (eta(t)eta(t/pq)), t = (-w + sqrt(D)) / 2a */
GEN
double_eta_quotient(GEN a, GEN w, GEN D, long p, long q, GEN pq, GEN sqrtD)
{
  GEN C = shifti(subii(sqri(w), D), -2);
  GEN d, t, z, zp, zq, zpq, s_t, s_tp, s_tpq, s, sp, spq;
  long prec = realprec(sqrtD);

  z = eta_form(to_form(a, w, C), sqrtD, &s_t, prec);
  s = gel(s_t, 1);
  zp = eta_form(to_form(mului(p, a), w, C), sqrtD, &s_tp, prec);
  sp = gel(s_tp, 1);
  zpq = eta_form(to_form(mulii(pq, a), w, C), sqrtD, &s_tpq, prec);
  spq = gel(s_tpq, 1);
  if (p == q) {
    z = gdiv(gsqr(zp), gmul(z, zpq));
    t = gsub(gmul2n(gel(s_tp,2), 1),
             gadd(gel(s_t,2), gel(s_tpq,2)));
    if (sp != gen_1) z = gmul(z, sp);
  } else {
    GEN s_tq, sq;
    zq = eta_form(to_form(mului(q, a), w, C), sqrtD, &s_tq, prec);
    sq = gel(s_tq, 1);
    z = gdiv(gmul(zp, zq), gmul(z, zpq));
    t = gsub(gadd(gel(s_tp,2), gel(s_tq,2)),
             gadd(gel(s_t,2), gel(s_tpq,2)));
    if (sp != gen_1) z = gmul(z, gsqrt(sp, prec));
    if (sq != gen_1) z = gmul(z, gsqrt(sq, prec));
  }
  d = NULL;
  if (s != gen_1) d = gsqrt(s, prec);
  if (spq != gen_1) {
    GEN x = gsqrt(spq, prec);
    d = d? gmul(d, x): x;
  }
  if (d) z = gdiv(z, d);
  return gmul(z, expIPiQ(t, prec));
}

typedef struct { GEN u; long v, t; } cxanalyze_t;

/* Check whether a t_COMPLEX, t_REAL or t_INT z != 0 can be written as
 * z = u * 2^(v/2) * exp(I Pi/4 t), u > 0, v = 0,1 and -3 <= t <= 4.
 * Allow z t_INT/t_REAL to simplify handling of eta_correction() output */
static int
cxanalyze(cxanalyze_t *T, GEN z)
{
  GEN a, b;
  long ta, tb;

  T->v = 0;
  if (is_intreal_t(typ(z)))
  {
    T->u = mpabs_shallow(z);
    T->t = signe(z) < 0? 4: 0;
    return 1;
  }
  a = gel(z,1); ta = typ(a);
  b = gel(z,2); tb = typ(b);

  T->t = 0;
  if (ta == t_INT && !signe(a))
  {
    T->u = R_abs_shallow(b);
    T->t = gsigne(b) < 0? -2: 2;
    return 1;
  }
  if (tb == t_INT && !signe(b))
  {
    T->u = R_abs_shallow(a);
    T->t = gsigne(a) < 0? 4: 0;
    return 1;
  }
  if (ta != tb || ta == t_REAL) { T->u = z; return 0; }
  /* a,b both non zero, both t_INT or t_FRAC */
  if (ta == t_INT)
  {
    if (!absequalii(a, b)) return 0;
    T->u = absi_shallow(a);
    T->v = 1;
    if (signe(a) == signe(b))
    { T->t = signe(a) < 0? -3: 1; }
    else
    { T->t = signe(a) < 0? 3: -1; }
  }
  else
  {
    if (!absequalii(gel(a,2), gel(b,2)) || !absequalii(gel(a,1),gel(b,1)))
      return 0;
    T->u = absfrac_shallow(a);
    T->v = 1;
    a = gel(a,1);
    b = gel(b,1);
    if (signe(a) == signe(b))
    { T->t = signe(a) < 0? -3: 1; }
    else
    { T->t = signe(a) < 0? 3: -1; }
  }
  return 1;
}

/* z * sqrt(st_b) / sqrt(st_a) exp(I Pi (t + t0)). Assume that
 * sqrt2 = gsqrt(gen_2, prec) or NULL */
static GEN
apply_eta_correction(GEN z, GEN st_a, GEN st_b, GEN t0, GEN sqrt2, long prec)
{
  GEN t, s_a = gel(st_a, 1), s_b = gel(st_b, 1);
  cxanalyze_t Ta, Tb;
  int ca, cb;

  t = gsub(gel(st_b,2), gel(st_a,2));
  if (t0 != gen_0) t = gadd(t, t0);
  ca = cxanalyze(&Ta, s_a);
  cb = cxanalyze(&Tb, s_b);
  if (ca || cb)
  { /* compute sqrt(s_b) / sqrt(s_a) in a more efficient way:
     * sb = ub sqrt(2)^vb exp(i Pi/4 tb) */
    GEN u = gdiv(Tb.u,Ta.u);
    switch(Tb.v - Ta.v)
    {
      case -1: u = gmul2n(u,-1); /* fall through: write 1/sqrt2 = sqrt2/2 */
      case 1: u = gmul(u, sqrt2? sqrt2: sqrtr_abs(real2n(1, prec)));
    }
    if (!isint1(u)) z = gmul(z, gsqrt(u, prec));
    t = gadd(t, gmul2n(stoi(Tb.t - Ta.t), -3));
  }
  else
  {
    z = gmul(z, gsqrt(s_b, prec));
    z = gdiv(z, gsqrt(s_a, prec));
  }
  return gmul(z, expIPiQ(t, prec));
}

/* sqrt(2) eta(2x) / eta(x) */
GEN
weberf2(GEN x, long prec)
{
  pari_sp av = avma;
  GEN z, sqrt2, a,b, Ua,Ub, st_a,st_b;

  x = upper_to_cx(x, &prec);
  a = cxredsl2(x, &Ua);
  b = cxredsl2(gmul2n(x,1), &Ub);
  if (gequal(a,b)) /* not infrequent */
    z = gen_1;
  else
    z = gdiv(eta_reduced(b,prec), eta_reduced(a,prec));
  st_a = eta_correction(a, Ua, 1);
  st_b = eta_correction(b, Ub, 1);
  sqrt2 = sqrtr_abs(real2n(1, prec));
  z = apply_eta_correction(z, st_a, st_b, gen_0, sqrt2, prec);
  return gerepileupto(av, gmul(z, sqrt2));
}

/* eta(x/2) / eta(x) */
GEN
weberf1(GEN x, long prec)
{
  pari_sp av = avma;
  GEN z, a,b, Ua,Ub, st_a,st_b;

  x = upper_to_cx(x, &prec);
  a = cxredsl2(x, &Ua);
  b = cxredsl2(gmul2n(x,-1), &Ub);
  if (gequal(a,b)) /* not infrequent */
    z = gen_1;
  else
    z = gdiv(eta_reduced(b,prec), eta_reduced(a,prec));
  st_a = eta_correction(a, Ua, 1);
  st_b = eta_correction(b, Ub, 1);
  z = apply_eta_correction(z, st_a, st_b, gen_0, NULL, prec);
  return gerepileupto(av, z);
}
/* exp(-I*Pi/24) * eta((x+1)/2) / eta(x) */
GEN
weberf(GEN x, long prec)
{
  pari_sp av = avma;
  GEN z, t0, a,b, Ua,Ub, st_a,st_b;
  x = upper_to_cx(x, &prec);
  a = cxredsl2(x, &Ua);
  b = cxredsl2(gmul2n(gaddgs(x,1),-1), &Ub);
  if (gequal(a,b)) /* not infrequent */
    z = gen_1;
  else
    z = gdiv(eta_reduced(b,prec), eta_reduced(a,prec));
  st_a = eta_correction(a, Ua, 1);
  st_b = eta_correction(b, Ub, 1);
  t0 = mkfrac(gen_m1, utoipos(24));
  z = apply_eta_correction(z, st_a, st_b, t0, NULL, prec);
  if (typ(z) == t_COMPLEX && isexactzero(real_i(x)))
    z = gerepilecopy(av, gel(z,1));
  else
    z = gerepileupto(av, z);
  return z;
}
GEN
weber0(GEN x, long flag,long prec)
{
  switch(flag)
  {
    case 0: return weberf(x,prec);
    case 1: return weberf1(x,prec);
    case 2: return weberf2(x,prec);
    default: pari_err_FLAG("weber");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/* check |q| < 1 */
static GEN
check_unit_disc(const char *fun, GEN q, long prec)
{
  GEN Q = gtofp(q, prec), Qlow;
  Qlow = (prec > LOWDEFAULTPREC)? gtofp(Q,LOWDEFAULTPREC): Q;
  if (gcmp(gnorm(Qlow), gen_1) >= 0)
    pari_err_DOMAIN(fun, "abs(q)", ">=", gen_1, q);
  return Q;
}

GEN
theta(GEN q, GEN z, long prec)
{
  long l, n;
  pari_sp av = avma, av2;
  GEN s, c, snz, cnz, s2z, c2z, ps, qn, y, zy, ps2, k, zold;

  l = precision(q);
  n = precision(z); if (n && n < l) l = n;
  if (l) prec = l;
  z = gtofp(z, prec);
  q = check_unit_disc("theta", q, prec);
  zold = NULL; /* gcc -Wall */
  zy = imag_i(z);
  if (gequal0(zy)) k = gen_0;
  else
  {
    GEN lq = glog(q,prec); k = roundr(divrr(zy, real_i(lq)));
    if (signe(k)) { zold = z; z = gadd(z, mulcxmI(gmul(lq,k))); }
  }
  qn = gen_1;
  ps2 = gsqr(q);
  ps = gneg_i(ps2);
  gsincos(z, &s, &c, prec);
  s2z = gmul2n(gmul(s,c),1); /* sin 2z */
  c2z = gsubgs(gmul2n(gsqr(c),1), 1); /* cos 2z */
  snz = s;
  cnz = c; y = s;
  av2 = avma;
  for (n = 3;; n += 2)
  {
    long e;
    s = gadd(gmul(snz, c2z), gmul(cnz,s2z));
    qn = gmul(qn,ps);
    y = gadd(y, gmul(qn, s));
    e = gexpo(s); if (e < 0) e = 0;
    if (gexpo(qn) + e < -prec2nbits(prec)) break;

    ps = gmul(ps,ps2);
    c = gsub(gmul(cnz, c2z), gmul(snz,s2z));
    snz = s; /* sin nz */
    cnz = c; /* cos nz */
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"theta (n = %ld)", n);
      gerepileall(av2, 5, &snz, &cnz, &ps, &qn, &y);
    }
  }
  if (signe(k))
  {
    y = gmul(y, gmul(powgi(q,sqri(k)),
                     gexp(gmul(mulcxI(zold),shifti(k,1)), prec)));
    if (mod2(k)) y = gneg_i(y);
  }
  return gerepileupto(av, gmul(y, gmul2n(gsqrt(gsqrt(q,prec),prec),1)));
}

GEN
thetanullk(GEN q, long k, long prec)
{
  long l, n;
  pari_sp av = avma;
  GEN p1, ps, qn, y, ps2;

  if (k < 0)
    pari_err_DOMAIN("thetanullk", "k", "<", gen_0, stoi(k));
  l = precision(q);
  if (l) prec = l;
  q = check_unit_disc("thetanullk", q, prec);

  if (!odd(k)) { set_avma(av); return gen_0; }
  qn = gen_1;
  ps2 = gsqr(q);
  ps = gneg_i(ps2);
  y = gen_1;
  for (n = 3;; n += 2)
  {
    GEN t;
    qn = gmul(qn,ps);
    ps = gmul(ps,ps2);
    t = gmul(qn, powuu(n, k)); y = gadd(y, t);
    if (gexpo(t) < -prec2nbits(prec)) break;
  }
  p1 = gmul2n(gsqrt(gsqrt(q,prec),prec),1);
  if (k&2) y = gneg_i(y);
  return gerepileupto(av, gmul(p1, y));
}

/* q2 = q^2 */
static GEN
vecthetanullk_loop(GEN q2, long k, long prec)
{
  GEN ps, qn = gen_1, y = const_vec(k, gen_1);
  pari_sp av = avma;
  const long bit = prec2nbits(prec);
  long i, n;

  if (gexpo(q2) < -2*bit) return y;
  ps = gneg_i(q2);
  for (n = 3;; n += 2)
  {
    GEN t = NULL/*-Wall*/, P = utoipos(n), N2 = sqru(n);
    qn = gmul(qn,ps);
    ps = gmul(ps,q2);
    for (i = 1; i <= k; i++)
    {
      t = gmul(qn, P); gel(y,i) = gadd(gel(y,i), t);
      P = mulii(P, N2);
    }
    if (gexpo(t) < -bit) return y;
    if (gc_needed(av,2))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"vecthetanullk_loop, n = %ld",n);
      gerepileall(av, 3, &qn, &ps, &y);
    }
  }
}
/* [d^i theta/dz^i(q, 0), i = 1, 3, .., 2*k - 1] */
GEN
vecthetanullk(GEN q, long k, long prec)
{
  long i, l = precision(q);
  pari_sp av = avma;
  GEN p1, y;

  if (l) prec = l;
  q = check_unit_disc("vecthetanullk", q, prec);
  y = vecthetanullk_loop(gsqr(q), k, prec);
  p1 = gmul2n(gsqrt(gsqrt(q,prec),prec),1);
  for (i = 2; i <= k; i += 2) gel(y,i) = gneg_i(gel(y,i));
  return gerepileupto(av, gmul(p1, y));
}

/* [d^i theta/dz^i(q, 0), i = 1, 3, .., 2*k - 1], q = exp(2iPi tau) */
GEN
vecthetanullk_tau(GEN tau, long k, long prec)
{
  long i, l = precision(tau);
  pari_sp av = avma;
  GEN q4, y;

  if (l) prec = l;
  if (typ(tau) != t_COMPLEX || gsigne(gel(tau,2)) <= 0)
    pari_err_DOMAIN("vecthetanullk_tau", "imag(tau)", "<=", gen_0, tau);
  q4 = expIPiC(gmul2n(tau,-1), prec); /* q^(1/4) */
  y = vecthetanullk_loop(gpowgs(q4,8), k, prec);
  for (i = 2; i <= k; i += 2) gel(y,i) = gneg_i(gel(y,i));
  return gerepileupto(av, gmul(gmul2n(q4,1), y));
}

/* Return E_k(tau). Slow if tau is not in standard fundamental domain */
GEN
cxEk(GEN tau, long k, long prec)
{
  pari_sp av;
  GEN q, y, qn;
  long n, b, l = precision(tau);

  if (l) prec = l;
  b = bit_accuracy(prec);
  /* sum n^(k-1) x^n <= x(1 + (k!-1)x) / (1-x)^k (cf Eulerian polynomials)
   * S = \sum_{n > 0} n^(k-1) |q^n/(1-q^n)| <= x(1+(k!-1)x) / (1-x)^(k+1),
   * where x = |q| = exp(-2Pi Im(tau)) < 1. Neglegt 2/zeta(1-k) * S if
   * (2Pi)^k/(k-1)! x < 2^(-b-1) and k! x < 1. Use log2((2Pi)^k/(k-1)!) < 10 */
  if (gcmpgs(imag_i(tau), (M_LN2 / (2*M_PI)) * (b+1+10)) > 0)
    return real_1(prec);
  if (k == 2)
  { /* -theta^(3)(tau/2) / theta^(1)(tau/2). Assume that Im tau > 0 */
    y = vecthetanullk_loop(qq(tau,prec), 2, prec);
    return gdiv(gel(y,2), gel(y,1));
  }
  q = cxtoreal(expIPiC(gneg(gmul2n(tau,1)), prec));
  av = avma; y = gen_0; qn = q;
  for(n = 1;; n++)
  { /* compute y := sum_{n>0} n^(k-1) / (q^n-1) */
    GEN t = gdiv(powuu(n,k-1), gsubgs(qn,1));
    if (gequal0(t) || gexpo(t) <= - prec2nbits(prec) - 5) break;
    y = gadd(y, t);
    qn = gmul(q,qn);
    if (gc_needed(av,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"elleisnum");
      gerepileall(av, 2, &y,&qn);
    }
  }
  return gadd(gen_1, gmul(y, gdiv(gen_2, szeta(1-k, prec))));
}
