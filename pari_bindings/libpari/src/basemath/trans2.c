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
/**                   TRANSCENDENTAL FUNCTIONS                     **/
/**                          (part 2)                              **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_trans

GEN
trans_fix_arg(long *prec, GEN *s0, GEN *sig, GEN *tau, pari_sp *av, GEN *res)
{
  GEN p1, s = *s0 = cxtoreal(*s0);
  long l;
  l = precision(s); if (!l) l = *prec;
  if (l < LOWDEFAULTPREC) l = LOWDEFAULTPREC;
  *res = cgetc(l); *av = avma;
  if (typ(s) == t_COMPLEX)
  { /* s = sig + i t */
    s = cxtofp(s, l+EXTRAPREC64);
    *sig = gel(s,1);
    *tau = gel(s,2);
  }
  else /* real number */
  {
    *sig = s = gtofp(s, l+EXTRAPREC64);
    *tau = gen_0;
    p1 = trunc2nr(s, 0);
    if (!signe(subri(s,p1))) *s0 = p1;
  }
  *prec = l; return s;
}

/********************************************************************/
/**                                                                **/
/**                          ARCTANGENT                            **/
/**                                                                **/
/********************************************************************/
/* atan(b/a), real a and b, suitable for gerepileupto */
static GEN
atan2_agm(GEN a, GEN b, long prec)
{ return gel(logagmcx(mkcomplex(a, b), prec), 2); }
static GEN
mpatan(GEN x)
{
  long l, l1, l2, n, m, i, lp, e, s, sx = signe(x);
  pari_sp av0, av;
  double alpha, beta, delta;
  GEN y, p1, p2, p3, p4, p5, unr;
  int inv;

  if (!sx) return real_0_bit(expo(x));
  l = lp = realprec(x);
  if (absrnz_equal1(x)) { /* |x| = 1 */
    y = Pi2n(-2, l+EXTRAPREC64); if (sx < 0) setsigne(y,-1);
    return y;
  }
  if (l > AGM_ATAN_LIMIT)
  { av = avma; return gerepileuptoleaf(av, atan2_agm(gen_1, x, l)); }

  e = expo(x); inv = (e >= 0); /* = (|x| > 1 ) */
  if (e > 0) lp += nbits2extraprec(e);

  y = cgetr(lp); av0 = avma;
  p1 = rtor(x, l+EXTRAPREC64); setabssign(p1); /* p1 = |x| */
  if (inv) p1 = invr(p1);
  e = expo(p1);
  if (e < -100)
    alpha = 1.65149612947 - e; /* log_2(Pi) - e */
  else
    alpha = log2(M_PI / atan(rtodbl(p1)));
  beta = (double)(prec2nbits(l)>>1);
  delta = 1 + beta - alpha/2;
  if (delta <= 0) { n = 1; m = 0; }
  else
  {
    double fi = alpha-2;
    if (delta >= fi*fi)
    {
      double t = 1 + sqrt(delta);
      n = (long)t;
      m = (long)(t - fi);
    }
    else
    {
      n = (long)(1+beta/fi);
      m = 0;
    }
  }
  l2 = l + nbits2extraprec(m);
  p2 = rtor(p1, l2); av = avma;
  for (i=1; i<=m; i++)
  {
    p5 = addsr(1, sqrr(p2)); setprec(p5,l2);
    p5 = addsr(1, sqrtr_abs(p5)); setprec(p5,l2);
    affrr(divrr(p2,p5), p2); set_avma(av);
  }
  p3 = sqrr(p2); l1 = minss(LOWDEFAULTPREC+EXTRAPREC64, l2); /* l1 increases to l2 */;
  unr = real_1(l2); setprec(unr,l1);
  p4 = cgetr(l2); setprec(p4,l1);
  affrr(divru(unr,2*n+1), p4);
  s = 0; e = expo(p3); av = avma;
  for (i = n; i > 1; i--) /* n >= 1. i = 1 done outside for efficiency */
  {
    setprec(p3,l1); p5 = mulrr(p4,p3);
    l1 += dvmdsBIL(s - e, &s); if (l1 > l2) l1 = l2;
    setprec(unr,l1); p5 = subrr(divru(unr,2*i-1), p5);
    setprec(p4,l1); affrr(p5,p4); set_avma(av);
  }
  setprec(p3, l2); p5 = mulrr(p4,p3); /* i = 1 */
  setprec(unr,l2); p4 = subrr(unr, p5);

  p4 = mulrr(p2,p4); shiftr_inplace(p4, m);
  if (inv) p4 = subrr(Pi2n(-1, lp), p4);
  if (sx < 0) togglesign(p4);
  affrr_fixlg(p4,y); set_avma(av0); return y;
}

GEN
gatan(GEN x, long prec)
{
  pari_sp av;
  GEN a, y;

  switch(typ(x))
  {
    case t_REAL: return mpatan(x);
    case t_COMPLEX: /* atan(x) = -i atanh(ix) */
      if (ismpzero(gel(x,2))) return gatan(gel(x,1), prec);
      av = avma; return gerepilecopy(av, mulcxmI(gatanh(mulcxI(x),prec)));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valser(y) < 0) pari_err_DOMAIN("atan","valuation", "<", gen_0, x);
      if (lg(y)==2) return gerepilecopy(av, y);
      /* lg(y) > 2 */
      a = integser(gdiv(derivser(y), gaddsg(1,gsqr(y))));
      if (!valser(y)) a = gadd(a, gatan(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return trans_eval("atan",gatan,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                             ARCSINE                            **/
/**                                                                **/
/********************************************************************/
/* |x| < 1, x != 0 */
static GEN
mpasin(GEN x) {
  pari_sp av = avma;
  GEN z, a = sqrtr(subsr(1, sqrr(x)));
  if (realprec(x) > AGM_ATAN_LIMIT)
    z = atan2_agm(a, x, realprec(x));
  else
    z = mpatan(divrr(x, a));
  return gerepileuptoleaf(av, z);
}

static GEN mpacosh(GEN x);
GEN
gasin(GEN x, long prec)
{
  long sx;
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_REAL: sx = signe(x);
      if (!sx) return real_0_bit(expo(x));
      if (absrnz_equal1(x)) { /* |x| = 1 */
        if (sx > 0) return Pi2n(-1, realprec(x)); /* 1 */
        y = Pi2n(-1, realprec(x)); setsigne(y, -1); return y; /* -1 */
      }
      if (expo(x) < 0) return mpasin(x);
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = Pi2n(-1, realprec(x));
      gel(y,2) = mpacosh(x);
      if (sx < 0) togglesign(gel(y,1)); else togglesign(gel(y,2));
      return y;

    case t_COMPLEX: /* asin(z) = -i asinh(iz) */
      if (ismpzero(gel(x,2))) return gasin(gel(x,1), prec);
      av = avma;
      return gerepilecopy(av, mulcxmI(gasinh(mulcxI(x), prec)));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gequal0(y)) return gerepilecopy(av, y);
      /* lg(y) > 2*/
      if (valser(y) < 0) pari_err_DOMAIN("asin","valuation", "<", gen_0, x);
      p1 = gsubsg(1,gsqr(y));
      if (gequal0(p1))
      {
        GEN t = Pi2n(-1,prec);
        if (gsigne(gel(y,2)) < 0) setsigne(t, -1);
        return gerepileupto(av, scalarser(t, varn(y), valser(p1)>>1));
      }
      p1 = gdiv(derivser(y), gsqrt(p1,prec));
      a = integser(p1);
      if (!valser(y)) a = gadd(a, gasin(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return trans_eval("asin",gasin,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                             ARCCOSINE                          **/
/**                                                                **/
/********************************************************************/
static GEN
acos0(long e) { return Pi2n(-1, nbits2prec(e<0? -e: 1)); }

/* |x| < 1, x != 0 */
static GEN
mpacos(GEN x)
{
  pari_sp av = avma;
  GEN z, a = sqrtr(subsr(1, sqrr(x)));
  if (realprec(x) > AGM_ATAN_LIMIT)
    z = atan2_agm(x, a, realprec(x));
  else
  {
    z = mpatan(divrr(a, x));
    if (signe(x) < 0) z = addrr(mppi(realprec(z)), z);
  }
  return gerepileuptoleaf(av, z);
}

GEN
gacos(GEN x, long prec)
{
  long sx;
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_REAL: sx = signe(x);
      if (!sx) return acos0(expo(x));
      if (absrnz_equal1(x)) /* |x| = 1 */
        return sx > 0? real_0_bit( -(bit_prec(x)>>1) ) : mppi(realprec(x));
      if (expo(x) < 0) return mpacos(x);

      y = cgetg(3,t_COMPLEX); p1 = mpacosh(x);
      if (sx < 0) { gel(y,1) = mppi(realprec(x)); togglesign(p1); }
      else gel(y,1) = gen_0;
      gel(y,2) = p1; return y;

    case t_COMPLEX:
      if (ismpzero(gel(x,2))) return gacos(gel(x,1), prec);
      av = avma;
      p1 = gadd(x, mulcxI(gsqrt(gsubsg(1,gsqr(x)), prec)));
      y = glog(p1,prec); /* log(x + I*sqrt(1-x^2)) */
      return gerepilecopy(av, mulcxmI(y));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valser(y) < 0) pari_err_DOMAIN("acos","valuation", "<", gen_0, x);
      if (lg(y) > 2)
      {
        p1 = gsubsg(1,gsqr(y));
        if (gequal0(p1)) { set_avma(av); return zeroser(varn(y), valser(p1)>>1); }
        p1 = integser(gdiv(gneg(derivser(y)), gsqrt(p1,prec)));
        /*y(t) = 1+O(t)*/
        if (gequal1(gel(y,2)) && !valser(y)) return gerepileupto(av, p1);
      }
      else p1 = y;
      a = (lg(y)==2 || valser(y))? Pi2n(-1, prec): gacos(gel(y,2),prec);
      return gerepileupto(av, gadd(a,p1));
  }
  return trans_eval("acos",gacos,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                            ARGUMENT                            **/
/**                                                                **/
/********************************************************************/

/* we know that x and y are not both 0 */
static GEN
mparg(GEN x, GEN y)
{
  long prec, sx = signe(x), sy = signe(y);
  GEN z;

  if (!sy)
  {
    if (sx > 0) return real_0_bit(expo(y) - expo(x));
    return mppi(realprec(x));
  }
  prec = realprec(y); if (prec < realprec(x)) prec = realprec(x);
  if (!sx)
  {
    z = Pi2n(-1, prec); if (sy < 0) setsigne(z,-1);
    return z;
  }

  if (expo(x)-expo(y) > -2)
  {
    z = mpatan(divrr(y,x)); if (sx > 0) return z;
    return addrr_sign(z, signe(z), mppi(prec), sy);
  }
  z = mpatan(divrr(x,y));
  return addrr_sign(z, -signe(z), Pi2n(-1, prec), sy);
}

static GEN
rfix(GEN x,long prec)
{
  switch(typ(x))
  {
    case t_INT: return itor(x, prec);
    case t_FRAC: return fractor(x, prec);
    case t_REAL: break;
    default: pari_err_TYPE("rfix (conversion to t_REAL)",x);
  }
  return x;
}

static GEN
cxarg(GEN x, GEN y, long prec)
{
  pari_sp av = avma;
  x = rfix(x,prec);
  y = rfix(y,prec); return gerepileuptoleaf(av, mparg(x,y));
}

GEN
garg(GEN x, long prec)
{
  long l;
  if (gequal0(x)) pari_err_DOMAIN("arg", "argument", "=", gen_0, x);
  switch(typ(x))
  {
    case t_REAL: prec = realprec(x); /* fall through */
    case t_INT: case t_FRAC: return (gsigne(x)>0)? real_0(prec): mppi(prec);
    case t_COMPLEX:
      l = precision(x); if (l) prec = l;
      return cxarg(gel(x,1),gel(x,2),prec);
  }
  return trans_eval("arg",garg,x,prec);
}

/********************************************************************/
/**                                                                **/
/**                      HYPERBOLIC COSINE                         **/
/**                                                                **/
/********************************************************************/
/* 1 + x */
static GEN
mpcosh0(long e) { return e >= 0? real_0_bit(e): real_1_bit(-e); }
static GEN
mpcosh(GEN x)
{
  pari_sp av;
  GEN z;

  if (!signe(x)) return mpcosh0(expo(x));
  av = avma;
  z = mpexp(x); z = addrr(z, invr(z)); shiftr_inplace(z, -1);
  return gerepileuptoleaf(av, z);
}

GEN
gcosh(GEN x, long prec)
{
  pari_sp av;
  GEN y, p1;
  long v;

  switch(typ(x))
  {
    case t_REAL: return mpcosh(x);
    case t_COMPLEX:
      if (isintzero(gel(x,1))) return gcos(gel(x,2),prec);
      /* fall through */
    case t_PADIC:
      av = avma; p1 = gexp(x,prec); p1 = gadd(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gequal0(y) && valser(y) == 0) return gerepilecopy(av, y);
      v = valser(y);
      if (v > 0) y = sertoser(y, lg(y) - 2 + v);
      p1 = gexp(y,prec); p1 = gadd(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
  }
  return trans_eval("cosh",gcosh,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                       HYPERBOLIC SINE                          **/
/**                                                                **/
/********************************************************************/
static GEN
mpsinh0(long e) { return real_0_bit(e); }
static GEN
mpsinh(GEN x)
{
  pari_sp av;
  long ex = expo(x), lx;
  GEN z, res;

  if (!signe(x)) return mpsinh0(ex);
  lx = realprec(x); res = cgetr(lx); av = avma;
  if (ex < 1 - BITS_IN_LONG)
  { /* y = e^x-1; e^x - e^(-x) = y(1 + 1/(y+1)) */
    GEN y = mpexpm1(x);
    z = addrs(y,1); if (realprec(z) > lx+1) z = rtor(z,lx+1); /* e^x */
    z = mulrr(y, addsr(1,invr(z)));
  }
  else
  {
    z = mpexp(x);
    z = subrr(z, invr(z));
  }
  shiftr_inplace(z, -1);
  affrr(z, res); set_avma(av); return res;
}

void
mpsinhcosh(GEN x, GEN *s, GEN *c)
{
  pari_sp av;
  long lx, ex = expo(x);
  GEN z, zi, S, C;
  if (!signe(x))
  {
    *c = mpcosh0(ex);
    *s = mpsinh0(ex); return;
  }
  lx = realprec(x);
  *c = cgetr(lx);
  *s = cgetr(lx); av = avma;
  if (ex < 1 - BITS_IN_LONG)
  { /* y = e^x-1; e^x - e^(-x) = y(1 + 1/(y+1)) */
    GEN y = mpexpm1(x);
    z = addrs(y,1); if (realprec(z) > lx+1) z = rtor(z,lx+1); /* e^x */
    zi = invr(z); /* z = exp(x), zi = exp(-x) */
    S = mulrr(y, addsr(1,zi));
  }
  else
  {
    z = mpexp(x);
    zi = invr(z);
    S = subrr(z, zi);
  }
  C = addrr(z, zi);
  shiftr_inplace(S, -1); affrr(S, *s);
  shiftr_inplace(C, -1); affrr(C, *c); set_avma(av);
}

GEN
gsinh(GEN x, long prec)
{
  pari_sp av;
  GEN y, p1;

  switch(typ(x))
  {
    case t_REAL: return mpsinh(x);
    case t_COMPLEX:
      if (isintzero(gel(x,1))) retmkcomplex(gen_0, gsin(gel(x,2),prec));
      /* fall through */
    case t_PADIC:
      av = avma; p1 = gexp(x,prec); p1 = gsub(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gequal0(y) && valser(y) == 0) return gerepilecopy(av, y);
      p1 = gexp(y, prec); p1 = gsub(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
  }
  return trans_eval("sinh",gsinh,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                      HYPERBOLIC TANGENT                        **/
/**                                                                **/
/********************************************************************/

static GEN
mptanh(GEN x)
{
  long lx, s = signe(x);
  GEN y;

  if (!s) return real_0_bit(expo(x));
  lx = realprec(x);
  if (abscmprr(x, utor(prec2nbits(lx), LOWDEFAULTPREC)) >= 0) {
    y = real_1(lx);
  } else {
    pari_sp av = avma;
    long ex = expo(x);
    GEN t;
    if (ex < 1 - BITS_IN_LONG) x = rtor(x, lx + nbits2extraprec(-ex)-1);
    t = exp1r_abs(gmul2n(x,1)); /* exp(|2x|) - 1 */
    y = gerepileuptoleaf(av, divrr(t, addsr(2,t)));
  }
  if (s < 0) togglesign(y); /* tanh is odd */
  return y;
}

GEN
gtanh(GEN x, long prec)
{
  pari_sp av;
  GEN y, t;

  switch(typ(x))
  {
    case t_REAL: return mptanh(x);
    case t_COMPLEX:
      if (isintzero(gel(x,1))) retmkcomplex(gen_0, gtan(gel(x,2),prec));
      /* fall through */
    case t_PADIC:
      av = avma;
      t = gexp(gmul2n(x,1),prec);
      t = gdivsg(-2, gaddgs(t,1));
      return gerepileupto(av, gaddsg(1,t));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gequal0(y)) return gerepilecopy(av, y);
      t = gexp(gmul2n(y, 1),prec);
      t = gdivsg(-2, gaddgs(t,1));
      return gerepileupto(av, gaddsg(1,t));
  }
  return trans_eval("tanh",gtanh,x,prec);
}

static GEN
mpcotanh(GEN x)
{
  long lx, s = signe(x);
  GEN y;

  if (!s) pari_err_DOMAIN("cotan", "argument", "=", gen_0, x);

  lx = realprec(x);
  if (abscmprr(x, utor(prec2nbits(lx), LOWDEFAULTPREC)) >= 0) {
    y = real_1(lx);
  } else {
    pari_sp av = avma;
    long ex = expo(x);
    GEN t;
    if (ex < 1 - BITS_IN_LONG) x = rtor(x, lx + nbits2extraprec(-ex)-1);
    t = exp1r_abs(gmul2n(x,1)); /* exp(|2x|) - 1 */
    y = gerepileuptoleaf(av, divrr(addsr(2,t), t));
  }
  if (s < 0) togglesign(y); /* cotanh is odd */
  return y;
}

GEN
gcotanh(GEN x, long prec)
{
  pari_sp av;
  GEN y, t;

  switch(typ(x))
  {
    case t_REAL: return mpcotanh(x);
    case t_COMPLEX:
      if (isintzero(gel(x,1))) retmkcomplex(gen_0, gcotan(gel(x,2),prec));
      /* fall through */
    case t_PADIC:
      av = avma;
      t = gexpm1(gmul2n(x,1),prec);
      return gerepileupto(av, gaddsg(1, gdivsg(2,t)));
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gequal0(y)) return gerepilecopy(av, y);
      t = gexpm1(gmul2n(y,1),prec);
      return gerepileupto(av, gaddsg(1, gdivsg(2,t)));
  }
  return trans_eval("cotanh",gcotanh,x,prec);
}

/********************************************************************/
/**                                                                **/
/**                     AREA HYPERBOLIC SINE                       **/
/**                                                                **/
/********************************************************************/

/* x != 0 */
static GEN
mpasinh(GEN x)
{
  long lx = realprec(x), ex = expo(x);
  GEN z, res = cgetr(lx);
  pari_sp av = avma;
  if (ex < 1 - BITS_IN_LONG) x = rtor(x, lx + nbits2extraprec(-ex)-1);
  z = logr_abs( addrr_sign(x,1, sqrtr_abs( addrs(sqrr(x), 1) ), 1) );
  if (signe(x) < 0) togglesign(z);
  affrr(z, res); return gc_const(av, res);
}

GEN
gasinh(GEN x, long prec)
{
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_REAL:
      if (!signe(x)) return rcopy(x);
      return mpasinh(x);

    case t_COMPLEX: {
      GEN a, b, d;
      if (ismpzero(gel(x,2))) return gasinh(gel(x,1), prec);
      av = avma;
      if (ismpzero(gel(x,1))) /* avoid cancellation */
        return gerepilecopy(av, mulcxI(gasin(gel(x,2), prec)));
      d = gsqrt(gaddsg(1,gsqr(x)), prec); /* Re(d) >= 0 */
      a = gadd(d, x);
      b = gsub(d, x);
      /* avoid cancellation as much as possible */
      if (gprecision(a) < gprecision(b))
        y = gneg(glog(b,prec));
      else
        y = glog(a,prec);
      return gerepileupto(av, y); /* log (x + sqrt(1+x^2)) */
    }
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gequal0(y)) return gerepilecopy(av, y);
      if (valser(y) < 0) pari_err_DOMAIN("asinh","valuation", "<", gen_0, x);
      p1 = gaddsg(1,gsqr(y));
      if (gequal0(p1))
      {
        GEN t = PiI2n(-1,prec);
        if ( gsigne(imag_i(gel(y,2))) < 0 ) setsigne(gel(t,2), -1);
        return gerepileupto(av, scalarser(t, varn(y), valser(p1)>>1));
      }
      p1 = gdiv(derivser(y), gsqrt(p1,prec));
      a = integser(p1);
      if (!valser(y)) a = gadd(a, gasinh(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return trans_eval("asinh",gasinh,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                     AREA HYPERBOLIC COSINE                     **/
/**                                                                **/
/********************************************************************/

/* |x| >= 1, return ach(|x|) */
static GEN
mpacosh(GEN x)
{
  long lx = realprec(x), e;
  GEN z, res = cgetr(lx);
  pari_sp av = avma;
  e = expo(signe(x) > 0? subrs(x,1): addrs(x,1));
  if (e == -(long)HIGHEXPOBIT)
    return gc_const((pari_sp)(res + lx), real_0_bit(- bit_prec(x) >> 1));
  if (e < -5) x = rtor(x, realprec(x) + nbits2extraprec(-e));
  z = logr_abs( addrr_sign(x, 1, sqrtr( subrs(sqrr(x), 1) ), 1) );
  affrr(z, res); return gc_const(av, res);
}

GEN
gacosh(GEN x, long prec)
{
  pari_sp av;
  GEN y;

  switch(typ(x))
  {
    case t_REAL: {
      long s = signe(x), e = expo(x);
      GEN a, b;
      if (s > 0 && e >= 0) return mpacosh(x);
      /* x < 1 */
      y = cgetg(3,t_COMPLEX); a = gen_0;
      if (s == 0) b = acos0(e);
      else if (e < 0) b = mpacos(x); /* -1 < x < 1 */
      else {
        if (!absrnz_equal1(x)) a = mpacosh(x);
        b = mppi(realprec(x));
      }
      gel(y,1) = a;
      gel(y,2) = b; return y;
    }
    case t_COMPLEX: {
      GEN a, b, d;
      if (ismpzero(gel(x,2))) return gacosh(gel(x,1), prec);
      av = avma;
      d = gsqrt(gaddsg(-1,gsqr(x)), prec); /* Re(d) >= 0 */
      a = gadd(x, d);
      b = gsub(x, d);
      /* avoid cancellation as much as possible */
      if (gprecision(a) < gprecision(b))
        y = glog(b,prec);
      else
        y = glog(a,prec);
      /* y = \pm log(x + sqrt(x^2-1)) */
      if (gsigne(real_i(y)) < 0) y = gneg(y);
      return gerepileupto(av, y);
    }
    default: {
      GEN a, d;
      long v;
      av = avma; if (!(y = toser_i(x))) break;
      v = valser(y);
      if (v < 0) pari_err_DOMAIN("acosh","valuation", "<", gen_0, x);
      if (gequal0(y))
      {
        if (!v) return gerepilecopy(av, y);
        return gerepileupto(av, gadd(y, PiI2n(-1, prec)));
      }
      d = gsubgs(gsqr(y),1);
      if (gequal0(d)) { set_avma(av); return zeroser(varn(y), valser(d)>>1); }
      d = gdiv(derivser(y), gsqrt(d,prec));
      a = integser(d);
      if (v)
        d = PiI2n(-1, prec); /* I Pi/2 */
      else
      {
        d = gel(y,2); if (gequal1(d)) return gerepileupto(av,a);
        d = gacosh(d, prec);
      }
      return gerepileupto(av, gadd(d,a));
    }
  }
  return trans_eval("acosh",gacosh,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                     AREA HYPERBOLIC TANGENT                    **/
/**                                                                **/
/********************************************************************/

/* |x| < 1, x != 0 */
static GEN
mpatanh(GEN x)
{
  pari_sp av = avma;
  long e, s = signe(x);
  GEN z;
  z = s > 0? subsr(1,x): addsr(1,x); e = expo(z);
  if (e < -5)
  {
    x = rtor(x, realprec(x) + nbits2extraprec(-e)-1);
    z = s > 0? subsr(1,x): addsr(1,x); e = expo(z);
  }
  z = invr(z); shiftr_inplace(z, 1); /* 2/(1-|x|) */
  z = logr_abs( addrs(z,-1) ); if (s < 0) togglesign(z);
  shiftr_inplace(z, -1); return gerepileuptoleaf(av, z);
}

/* atanh(u/v) using binary splitting, 0 < u < v */
GEN
atanhuu(ulong u, ulong v, long prec)
{
  long i, nmax;
  GEN u2 = sqru(u), v2 = sqru(v);
  double d = ((double)v) / u;
  struct abpq_res R;
  struct abpq A;
  /* satisfies (2n+1) (v/u)^2n > 2^bitprec */
  nmax = (long)ceil(prec2nbits(prec) / (2*log2(d)));
  abpq_init(&A, nmax);
  A.a[0] = A.b[0] = gen_1;
  A.p[0] = utoipos(u);
  A.q[0] = utoipos(v);
  for (i = 1; i <= nmax; i++)
  {
    A.a[i] = gen_1;
    A.b[i] = utoipos((i<<1)+1);
    A.p[i] = u2;
    A.q[i] = v2;
  }
  abpq_sum(&R, 0, nmax, &A);
  return rdivii(R.T, mulii(R.B,R.Q),prec);
}
/* atanh(u/v) using binary splitting, 0 < u < v */
GEN
atanhui(ulong u, GEN v, long prec)
{
  long i, nmax;
  GEN u2 = sqru(u), v2 = sqri(v);
  double d = gtodouble(v) / u;
  struct abpq_res R;
  struct abpq A;
  /* satisfies (2n+1) (v/u)^2n > 2^bitprec */
  nmax = (long)ceil(prec2nbits(prec) / (2*log2(d)));
  abpq_init(&A, nmax);
  A.a[0] = A.b[0] = gen_1;
  A.p[0] = utoipos(u);
  A.q[0] = v;
  for (i = 1; i <= nmax; i++)
  {
    A.a[i] = gen_1;
    A.b[i] = utoipos((i<<1)+1);
    A.p[i] = u2;
    A.q[i] = v2;
  }
  abpq_sum(&R, 0, nmax, &A);
  return rdivii(R.T, mulii(R.B,R.Q),prec);
}

static void
err_atanh(GEN x, GEN bad) { pari_err_DOMAIN("atanh", "x", "=", bad, x); }

GEN
gatanh(GEN x, long prec)
{
  long sx;
  pari_sp av;
  GEN a, y, z;

  switch(typ(x))
  {
    case t_INT:
      sx = signe(x);
      if (!sx) return real_0(prec);
      z = cgetg(3, t_COMPLEX); av = avma;
      if (lgefint(x) == 3)
      {
        if (x[2] == 1) err_atanh(x, sx == 1? gen_1: gen_m1);
        a = atanhuu(1, x[2], prec);
      }
      else
        a = atanhui(1, x, prec);
      gel(z,1) = gerepileuptoleaf(av, a);
      gel(z,2) = Pi2n(-1, prec);
      togglesign(sx > 0? gel(z,2): gel(z,1));
      return z;
    case t_FRAC:
    {
      long ly, lz;

      y = gel(x,1); ly = lgefint(y);
      z = gel(x,2); lz = lgefint(z); if (ly > 3 && lz > 3) break;
      if (abscmpii(y, z) > 0) /* |y| > z; lz = 3 */
      {
        ulong u = z[2];
        z = cgetg(3, t_COMPLEX); av = avma;
        a = ly == 3? atanhuu(u, y[2], prec): atanhui(u, y, prec);
        gel(z,1) = gerepileuptoleaf(av, a);
        gel(z,2) = Pi2n(-1, prec);
        togglesign(signe(y) > 0? gel(z,2): gel(z,1));
      }
      else
      { /* |y| < z; ly = 3 */
        av = avma;
        a = lz == 3? atanhuu(y[2], z[2], prec): atanhui(y[2], z, prec);
        z = gerepileuptoleaf(av, a);
        if (signe(y) < 0) togglesign(z);
      }
      return z;
    }
    case t_REAL:
      sx = signe(x);
      if (!sx) return real_0_bit(expo(x));
      if (expo(x) < 0) return mpatanh(x);

      y = cgetg(3,t_COMPLEX);
      av = avma;
      z = subrs(x,1);
      if (!signe(z)) err_atanh(x, gen_1);
      z = invr(z); shiftr_inplace(z, 1); /* 2/(x-1)*/
      z = addrs(z,1);
      if (!signe(z)) err_atanh(x, gen_m1);
      z = logr_abs(z);
      shiftr_inplace(z, -1); /* (1/2)log((1+x)/(x-1)) */
      gel(y,1) = gerepileuptoleaf(av, z);
      gel(y,2) = Pi2n(-1, realprec(x));
      if (sx > 0) togglesign(gel(y,2));
      return y;

    case t_COMPLEX: /* 2/(1-z) - 1 = (1+z) / (1-z) */
      if (ismpzero(gel(x,2))) return gatanh(gel(x,1), prec);
      av = avma; z = glog( gaddgs(gdivsg(2,gsubsg(1,x)),-1), prec );
      return gerepileupto(av, gmul2n(z,-1));

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valser(y) < 0) pari_err_DOMAIN("atanh","valuation", "<", gen_0, x);
      z = gdiv(derivser(y), gsubsg(1,gsqr(y)));
      a = integser(z);
      if (!valser(y)) a = gadd(a, gatanh(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return trans_eval("atanh",gatanh,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                         EULER'S GAMMA                          **/
/**                                                                **/
/********************************************************************/
/* 0 < a < b */
static GEN
mulu_interval_step_i(ulong a, ulong b, ulong step)
{
  ulong k, l, N, n;
  long lx;
  GEN x;

  n = 1 + (b-a) / step;
  b -= (b-a) % step;
  /* step | b-a */
  lx = 1; x = cgetg(2 + n/2, t_VEC);
  N = b + a;
  for (k = a;; k += step)
  {
    l = N - k; if (l <= k) break;
    gel(x,lx++) = muluu(k,l);
  }
  if (l == k) gel(x,lx++) = utoipos(k);
  setlg(x, lx); return x;
}
static GEN
_mul(void *data, GEN x, GEN y)
{
  long prec = (long)data;
  /* switch to t_REAL ? */
  if (typ(x) == t_INT && lgefint(x) > prec) x = itor(x, prec);
  if (typ(y) == t_INT && lgefint(y) > prec) y = itor(y, prec);
  return mpmul(x, y);
}
static GEN
mulu_interval_step_prec(long l, long m, long s, long prec)
{
  GEN v = mulu_interval_step_i(l, m, s);
  return gen_product(v, (void*)prec, &_mul);
}

/* x * (i*(i+1)) */
static GEN
muliunextu(GEN x, ulong i)
{
  if (i & HIGHMASK) /* i(i+1) >= 2^BITS_IN_LONG*/
    return mulii(x, muluu(i, i+1));
  else
    return muliu(x, i*(i+1));
}
/* arg(s+it) */
double
darg(double s, double t)
{
  double x;
  if (!t) return (s>0)? 0.: M_PI;
  if (!s) return (t>0)? M_PI/2: -M_PI/2;
  x = atan(t/s);
  return (s>0)? x
              : ((t>0)? x+M_PI : x-M_PI);
}

void
dcxlog(double s, double t, double *a, double *b)
{
  *a = log(s*s + t*t) / 2; /* log |s| = Re(log(s)) */
  *b = darg(s,t);          /* Im(log(s)) */
}

double
dabs(double s, double t) { return sqrt( s*s + t*t ); }
double
dnorm(double s, double t) { return s*s + t*t; }

#if 0
/* x, z t_REAL. Compute unique x in ]-z,z] congruent to x mod 2z */
static GEN
red_mod_2z(GEN x, GEN z)
{
  GEN Z = gmul2n(z, 1), d = subrr(z, x);
  /* require little accuracy */
  if (!signe(d)) return x;
  setprec(d, nbits2prec(expo(d) - expo(Z)));
  return addrr(mulir(floorr(divrr(d, Z)), Z), x);
}
#endif

static GEN
negeuler(long prec) { GEN g = mpeuler(prec); setsigne(g, -1); return g; }
/* lngamma(1+z) = -Euler*z + sum_{i > 1} zeta(i)/i (-z)^i
 * at relative precision prec, |z| <= 1/2 is small */
static GEN
lngamma1(GEN z, long prec)
{ /* sum_{i > l} |z|^(i-1) = |z|^l / (1-|z|) < 2^-B
   * for l > (B+1) / |log2(|z|)| */
  long i, l = ceil((bit_accuracy(prec) + 1) / - dbllog2(z));
  GEN s, vz;

  if (l <= 1) return gmul(negeuler(prec), z);
  vz = constzeta(l, prec);
  for (i = l, s = gen_0; i > 0; i--)
  {
    GEN c = divru(gel(vz,i), i);
    if (odd(i)) setsigne(c, -1);
    s = gadd(gmul(s,z), c);
  }
  return gmul(z, s);
}
/* B_i / (i(i-1)), i even. Sometimes NOT reduced (but gadd/gmul won't care)!*/
static GEN
bern_unextu(long i)
{ GEN B = bernfrac(i); return mkfrac(gel(B,1), muliunextu(gel(B,2), i-1)); }
/* B_i / i, i even. Sometimes NOT reduced (but gadd/gmul won't care)!*/
static GEN
bern_u(long i)
{ GEN B = bernfrac(i); return mkfrac(gel(B,1), muliu(gel(B,2), i)); }
/* sum_{i > 0} B_{2i}/(2i(2i-1)) * a^(i-1) */
static GEN
lngamma_sum(GEN a, long N)
{
  pari_sp av = avma;
  GEN S = bern_unextu(2*N);
  long i;
  for (i = 2*N-2; i > 0; i -= 2)
  {
    S = gadd(bern_unextu(i), gmul(a,S));
    if (gc_needed(av,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"gamma: i = %ld", i);
      S = gerepileupto(av, S);
    }
  }
  return S;
}
/* sum_{i > 0} B_{2i}/(2i) * a^i */
static GEN
psi_sum(GEN a, long N)
{
  pari_sp av = avma;
  GEN S = bern_u(2*N);
  long i;
  for (i = 2*N-2; i > 0; i -= 2)
  {
    S = gadd(bern_u(i), gmul(a,S));
    if (gc_needed(av,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"psi: i = %ld", i);
      S = gerepileupto(av, S);
    }
  }
  return gmul(a,S);
}
static void
gamma_optim(double ssig, double st, long prec, long *plim, long *pN)
{
  double la, l,l2,u,v, rlogs, ilogs;
  long N = 1, lim;
  dcxlog(ssig,st, &rlogs,&ilogs);
  /* Re (s - 1/2) log(s) */
  u = (ssig - 0.5)*rlogs - st * ilogs;
  /* Im (s - 1/2) log(s) */
  v = (ssig - 0.5)*ilogs + st * rlogs;
  /* l2 = | (s - 1/2) log(s) - s + log(2Pi)/2 |^2 ~ |lngamma(s))|^2 */
  u = u - ssig + log(2.*M_PI)/2;
  v = v - st;
  l2 = u*u + v*v;
  if (l2 < 0.000001) l2 = 0.000001;
  l = (prec2nbits_mul(prec, M_LN2) - log(l2)/2) / 2.;
  if (l < 0) l = 0.;

  if (st > 1 && l > 0)
  {
    double t = st * M_PI / l;
    la = t * log(t);
    if (la < 4.) la = 4.;
    if (la > 150) la = t;
  }
  else
    la = 4.; /* heuristic */
  lim = (long)ceil(l / (1.+ log(la)));
  if (lim == 0) lim = 1;

  u = (lim-0.5) * la / M_PI;
  l2 = u*u - st*st;
  if (l2 > 0)
  {
    double t = ceil(sqrt(l2) - ssig);
    if (t > 1) N = (long)t;
  }
  *plim = lim; *pN = N;
}
/* do we use lngamma1 instead of Euler-Maclaurin ? */
static int
gamma_use_1(double s, double t, long prec, long *plim, long *pN)
{
  double a = s-1, d = fabs(a) + fabs(t);
  long k;
  if (d < 1e-16) return 1;
  gamma_optim(s, t, prec, plim, pN);
  if (d >= 0.5) return 0;
  k = bit_accuracy(prec) / -log2(dnorm(a, t)); /* 2k = lngamma1 bound */
  return (t ? k: 1.5*k) < *plim + *pN;
}
static GEN
cxgamma(GEN s0, int dolog, long prec)
{
  GEN s, a, y, res, sig, tau, B, nnx, pi, pi2;
  long i, esig, et, lim, N = 1;
  pari_sp av, av2;
  int funeq = 0;
  pari_timer T;

  if (DEBUGLEVEL>5) timer_start(&T);
  s = trans_fix_arg(&prec,&s0,&sig,&tau,&av,&res);

  esig = expo(sig);
  et = signe(tau)? expo(tau): 0;
  if ((signe(sig) <= 0 || esig < -1) && et <= 16)
  { /* s <--> 1-s */
    funeq = 1; s = gsubsg(1, s); sig = real_i(s);
  }

  /* find "optimal" parameters [lim, N] */
  if (esig > 300 || et > 300)
  { /* |s| is HUGE ! Play safe and avoid inf / NaN */
    GEN S, iS, l2, la, u;
    double logla, l;

    S = gprec_w(s,LOWDEFAULTPREC);
    /* l2 ~ |lngamma(s))|^2 */
    l2 = gnorm(gmul(S, glog(S, LOWDEFAULTPREC)));
    l = (prec2nbits_mul(prec, M_LN2) - rtodbl(glog(l2,LOWDEFAULTPREC))/2) / 2.;
    if (l < 0) l = 0.;

    iS = imag_i(S);
    if (et > 0 && l > 0)
    {
      GEN t = gmul(iS, dbltor(M_PI / l)), logt = glog(t,LOWDEFAULTPREC);
      la = gmul(t, logt);
      if      (gcmpgs(la, 3) < 0)   { logla = log(3.); la = stoi(3); }
      else if (gcmpgs(la, 150) > 0) { logla = rtodbl(logt); la = t; }
      else logla = rtodbl(mplog(la));
    }
    else
    {
      logla = log(3.); la = stoi(3);
    }
    lim = (long)ceil(l / (1.+ logla));
    if (lim == 0) lim = 1;

    u = gmul(la, dbltor((lim-0.5)/M_PI));
    l2 = gsub(gsqr(u), gsqr(iS));
    if (signe(l2) > 0)
    {
      l2 = gsub(gsqrt(l2,3), sig);
      if (signe(l2) > 0) N = itos( gceil(l2) );
    }
  }
  else
  { /* |s| is moderate. Use floats  */
    double ssig = rtodbl(sig);
    double st = typ(s) == t_REAL? 0.0: rtodbl(imag_i(s));

    if (gamma_use_1(ssig, st, prec, &lim, &N))
    { /* s ~ 1: loggamma(1+u) ~ - Euler * u, cancellation */
      if (funeq) /* s0 ~ 0: use lngamma(s0)+log(s0) = lngamma(s0+1) */
        y = dolog? gsub(lngamma1(s0,prec), glog(s0,prec))
                 : gdiv(gexp(lngamma1(s0,prec), prec), s0);
      else
      {
        if (isint1(s0))
        {
          set_avma(av);
          return dolog? real_0(prec): real_1(prec);
        }
        y = lngamma1(gsubgs(s0,1),prec);
        if (!dolog) y = gexp(y,prec);
      }
      set_avma(av); return affc_fixlg(y, res);
    }
  }
  if (DEBUGLEVEL>5) err_printf("lim, N: [%ld, %ld]\n",lim,N);
  incrprec(prec);

  av2 = avma;
  y = s;
  if (typ(s0) == t_INT)
  {
    ulong ss = itou_or_0(s0);
    if (signe(s0) <= 0)
      pari_err_DOMAIN("gamma","argument", "=",
                       strtoGENstr("nonpositive integer"), s0);
    if (!ss || ss + (ulong)N < ss) {
      for (i=1; i < N; i++)
      {
        y = mulri(y, addiu(s0, i));
        if (gc_needed(av2,3))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
          y = gerepileuptoleaf(av2, y);
        }
      }
    } else {
      for (i=1; i < N; i++)
      {
        y = mulru(y, ss + i);
        if (gc_needed(av2,3))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
          y = gerepileuptoleaf(av2, y);
        }
      }
    }
  }
  else
  { /* Compute lngamma mod 2 I Pi */
    GEN sq = gsqr(s);
    pari_sp av3 = avma;
    for (i = 1; i < N - 1; i += 2)
    {
      y = gmul(y, gaddsg(i*(i + 1), gadd(gmulsg(2*i + 1, s), sq)));
      if (gc_needed(av2,3))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
        y = gerepileupto(av3, y);
      }
    }
    if (!odd(N)) y = gmul(y, gaddsg(N - 1, s));
  }
  if (DEBUGLEVEL>5) timer_printf(&T,"product from 0 to N-1");
  constbern(lim);
  nnx = gaddgs(s, N); a = ginv(nnx);
  B = gadd(gsub(gmul(gsub(nnx, ghalf), glog(nnx,prec)), nnx),
           gmul(a, lngamma_sum(gsqr(a), lim)));
  if (DEBUGLEVEL>5) timer_printf(&T,"Bernoulli sum");

  pi = mppi(prec); pi2 = shiftr(pi, 1);
  if (dolog)
  {
    if (typ(s) == t_REAL)
    {
      if (!funeq) y = logr_abs(divrr(sqrtr(pi2), y));
      else
      {
        GEN T = shiftr(sqrtr(pi2),-1); /* sqrt(Pi/2) */
        /* s0 < 0, step (*) simplifies: imag(lngamma(s0)) = - Pi * floor(s0) */
        y = logr_abs(divrr(mulrr(y, T), mpsin(gmul(pi,s0))));
        y = mkcomplex(y, mulri(pi, gfloor(s0)));
        B = gneg(B);
      }
    }
    else
    { /* log(y), fixing imaginary part */
      long prec2 = LOWDEFAULTPREC;
      GEN k, s2 = gprec_w(s, prec2), y2 = garg(s2, prec2); /* ~ Im log(s) */
      for (i=1; i < N; i++) y2 = gadd(y2, garg(gaddgs(s2,i), prec2));
      y = glog(y, prec);
      k = ground( gdiv(gsub(y2, imag_i(y)), Pi2n(1,prec2)) );
      if (signe(k)) y = gadd(y, mulcxI(mulir(k, Pi2n(1, prec))));
      if (!funeq) y = gsub(shiftr(logr_abs(pi2),-1), y); /* y -> sqrt(2Pi)/y */
      else
      { /* recall that s = 1 - s0 */
        GEN T = shiftr(sqrtr(pi2),-1); /* sqrt(Pi/2) */
        /* (*) Compute log(sin(Pi s0)) so that it has branch cuts along
        * (-oo, 0] and [1, oo). To do this in a numerically stable way
        * we must compute the log first then mangle its imaginary part.
        * The rounding operation below is stable because we're rounding
        * a number which is already within 1/4 of an integer. */

        /* z = log(sin(Pi s0) / sqrt(Pi/2)) */
        GEN z = glog(gdiv(gsin(gmul(pi,s0),prec), T), prec);
        GEN b = shiftr(subrs(shiftr(sig, 1), 1), -2); /* (2 Re(s)-1) / 4 */
        y = gsub(y, z);
        if (gsigne(imag_i(s)) > 0) togglesign(b);
        z = roundr(gsub(gdiv(imag_i(z), pi2), b)); /* round( Im(z)/2Pi - b ) */
        if (signe(z)) { /* y += I*z, z a t_REAL */
          z = mulir(z, pi2);
          if (typ(y) == t_COMPLEX) gel(y,2) = gadd(gel(y,2), z);
          else y = mkcomplex(y, z);
        }
        B = gneg(B);
      }
    }
    y = gadd(B, y);
  }
  else
  {
    GEN sqrtpi2 = sqrtr(pi2);
    if (funeq)
    { /* y --> y Pi/(sin(Pi s) * sqrt(2Pi)) = y sqrt(Pi/2)/sin(Pi s) */
      y = gdiv(gmul(shiftr(sqrtpi2,-1),y), gsin(gmul(pi,s0), prec));
      /* don't use s above: sin(pi s0) = sin(pi s) and the former is
       * more accurate, esp. if s0 ~ 0 */
      B = gneg(B);
    }
    else /* y --> sqrt(2Pi) / y */
      y = gdiv(sqrtpi2, y);
    y = gmul(gexp(B, prec), y);
  }
  set_avma(av); return affc_fixlg(y, res);
}

/* Theory says n > C * b^1.5 / log(b). Timings:
 * b = 64*[1, 2, 3, 4, 5, 6, 7, 10, 20, 30, 50, 100, 200, 500];
 * n = [1450, 1930, 2750, 3400, 4070, 5000, 6000, 8800, 26000, 50000, 130000,
 *      380000, 1300000, 6000000]; */
static long
gamma2_n(long prec)
{
  long b = bit_accuracy(prec);
  if (b <=  64) return 1450;
  if (b <= 128) return 1930;
  if (b <= 192) return 2750;
  if (b <= 256) return 3400;
  if (b <= 320) return 4070;
  if (b <= 384) return 5000;
  if (b <= 448) return 6000;
  return 10.0 * b * sqrt(b) / log(b);
}

/* m even, Gamma((m+1) / 2) */
static GEN
gammahs(long m, long prec)
{
  GEN y = cgetr(prec), z;
  pari_sp av = avma;
  long ma = labs(m);

  if (ma > gamma2_n(prec))
  {
    z = stor(m + 1, prec); shiftr_inplace(z, -1);
    affrr(cxgamma(z,0,prec), y);
    set_avma(av); return y;
  }
  z = sqrtr( mppi(prec) );
  if (m)
  {
    GEN t = mulu_interval_step_prec(1, ma-1, 2, prec + EXTRAPREC64);
    if (m >= 0) z = mpmul(z,t);
    else
    {
      z = mpdiv(z,t);
      if ((m&3) == 2) setsigne(z,-1);
    }
    shiftr_inplace(z, -m/2);
  }
  affrr(z, y); set_avma(av); return y;
}
GEN
ggammah(GEN x, long prec)
{
  switch(typ(x))
  {
    case t_INT:
    {
      long k = itos_or_0(x);
      if (!k && signe(x)) pari_err_OVERFLOW("gamma");
      return gammahs(k * 2, prec);
    }
    case t_REAL: case t_COMPLEX: case t_PADIC: case t_SER: {
      pari_sp av = avma;
      return gerepileupto(av, ggamma(gadd(x,ghalf), prec));
    }
  }
  return trans_eval("gammah",ggammah,x,prec);
}

/* find n such that n+v_p(n!)>=k p^2/(p-1)^2 */
static long
nboft(long k, long p)
{
  pari_sp av = avma;
  long s, n;

  if (k <= 0) return 0;
  k = itou( gceil(gdiv(mului(k, sqru(p)), sqru(p-1))) );
  set_avma(av);
  for (s=0, n=0; n+s < k; n++, s += u_lval(n, p));
  return n;
}

/* Using Dwork's expansion, compute \Gamma(px+1)=-\Gamma(px) with x a unit.
 * See p-Adic Gamma Functions and Dwork Cohomology, Maurizio Boyarsky
 * Transactions of the AMS, Vol. 257, No. 2. (Feb., 1980), pp. 359-369.
 * Inspired by a GP script by Fernando Rodriguez-Villegas */
static GEN
gadw(GEN x, long p)
{
  pari_sp ltop = avma;
  GEN s, t, u = cgetg(p+1, t_VEC);
  long j, k, kp, n = nboft(precp(x)+valp(x)+1, p);

  t = s = gaddsg(1, zeropadic(gel(x,2), n));
  gel(u, 1) = s;
  gel(u, 2) = s;
  for (j = 2; j < p; ++j)
    gel(u, j+1) = gdivgu(gel(u, j), j);
  for (k = 1, kp = p; k < n; ++k, kp += p) /* kp = k*p */
  {
    GEN c;
    gel(u, 1) = gdivgu(gadd(gel(u, 1), gel(u, p)), kp);
    for (j = 1; j < p; ++j)
      gel(u, j+1) = gdivgu(gadd(gel(u, j), gel(u, j+1)), kp + j);

    t = gmul(t, gaddgs(x, k-1));
    c = leafcopy(gel(u,1)); setvalp(c, valp(c) + k); /* c = u[1] * p^k */
    s = gadd(s, gmul(c, t));
    if ((k&0xFL)==0) gerepileall(ltop, 3, &u,&s,&t);
  }
  return gneg(s);
}

/*Use Dwork expansion*/
/*This is a O(p*e*log(pe)) algorithm, should be used when p small
 * If p==2 this is a O(pe) algorithm. */
static GEN
Qp_gamma_Dwork(GEN x, long p)
{
  pari_sp ltop = avma;
  long k = padic_to_Fl(x, p);
  GEN p1;
  long j;
  long px = precp(x);
  if (p==2 && px)
  {
    x = shallowcopy(x);
    setprecp(x, px+1);
    gel(x,3) = shifti(gel(x,3),1);
  }
  if (k)
  {
    GEN x_k = gsubgs(x,k);
    x = gdivgu(x_k, p);
    p1 = gadw(x, p); if (!odd(k)) p1 = gneg(p1);
    for (j = 1; j < k; ++j) p1 = gmul(p1, gaddgs(x_k, j));
  }
  else
    p1 = gneg(gadw(gdivgu(x, p), p));
  return gerepileupto(ltop, p1);
}

/* Compute Qp_gamma using the definition. This is a O(x*M(log(pe))) algorithm.
 * This should be used if x is very small. */
static GEN
Qp_gamma_Morita(long n, GEN p, long e)
{
  pari_sp ltop=avma;
  GEN p2 = gaddsg((n&1)?-1:1, zeropadic(p, e));
  long i;
  long pp=is_bigint(p)? 0: itos(p);
  for (i = 2; i < n; i++)
    if (!pp || i%pp)
    {
      p2 = gmulgu(p2, i);
      if ((i&0xFL) == 0xFL)
        p2 = gerepileupto(ltop, p2);
    }
  return gerepileupto(ltop, p2);
}

/* x\in\N: Gamma(-x)=(-1)^(1+x+x\p)*Gamma(1+x) */
static GEN
Qp_gamma_neg_Morita(long n, GEN p, long e)
{
  GEN g = ginv(Qp_gamma_Morita(n+1, p, e));
  return ((n^sdivsi(n,p)) & 1)? g: gneg(g);
}

/* p-adic Gamma function for x a p-adic integer */
/* If n < p*e : use Morita's definition.
 * Else : use Dwork's expansion.
 * If both n and p are big : itos(p) will fail.
 * TODO: handle p=2 better (Qp_gamma_Dwork is slow for p=2). */
GEN
Qp_gamma(GEN x)
{
  GEN n, m, N, p = gel(x,2);
  long s, e = precp(x);
  if (absequaliu(p, 2) && e == 2) e = 1;
  if (valp(x) < 0) pari_err_DOMAIN("gamma","v_p(x)", "<", gen_0, x);
  n = gtrunc(x);
  m = gtrunc(gneg(x));
  N = cmpii(n,m)<=0?n:m;
  s = itos_or_0(N);
  if (s && cmpsi(s, muliu(p,e)) < 0) /* s < p*e */
    return (N == n) ? Qp_gamma_Morita(s,p,e): Qp_gamma_neg_Morita(s,p,e);
  return Qp_gamma_Dwork(x, itos(p));
}

static GEN
Qp_lngamma(GEN x)
{
  GEN s, y, Y;
  long v = valp(x), e, k, K;
  if (v >= 0) return Qp_log(Qp_gamma(x));
  e = precp(x) + v; K = (2 + (e + 4) / (-v)) >> 1;
  s = gen_0; Y = y = ginv(x); y = gsqr(y); constbern(K);
  for (k = 1; k <= K; k++)
  {
    s = gadd(s, gmul(gdivgunextu(bernfrac(2*k), 2*k-1), Y));
    if (k < K) Y = gmul(Y, y); /* x^(1-2k) */
  }
  return gadd(s, gsub(gmul(gsub(x, ghalf), Qp_log(x)), x));
}

/* gamma(1+x) - 1, |x| < 1 is "small" */
GEN
ggamma1m1(GEN x, long prec) { return gexpm1(lngamma1(x, prec), prec); }

/* lngamma(y) with 0 constant term, using (lngamma y)' = y' psi(y) */
static GEN
serlngamma0(GEN y, long prec)
{
  GEN t;
  if (valser(y)) pari_err_DOMAIN("lngamma","valuation", "!=", gen_0, y);
  t = derivser(y);
  /* don't compute psi if y'=0 */
  if (signe(t)) t = gmul(t, gpsi(y,prec));
  return integser(t);
}

static GEN
sergamma(GEN y, long prec)
{
  GEN z, y0, Y;
  if (lg(y) == 2) pari_err_DOMAIN("gamma", "argument", "=", gen_0,y);
  /* exp(lngamma) */
  if (valser(y) > 0) return gdiv(gexp(glngamma(gaddgs(y,1),prec),prec),y);
  y0 = simplify_shallow(gel(y,2));
  z = NULL; Y = y;
  if (isint(y0, &y0))
  { /* fun eq. avoids log singularity of lngamma at negative ints */
    long s = signe(y0);
    /* possible if y[2] is an inexact 0 */
    if (!s) return gdiv(gexp(glngamma(gaddgs(y,1),prec),prec),y);
    if (signe(y0) < 0) { Y = gsubsg(1, y); y0 = subsi(1, y0); }
    if (abscmpiu(y0, 50) < 0) z = mpfact(itos(y0)-1); /* more precise */
  }
  if (!z) z = ggamma(y0,prec);
  z = gmul(z, gexp(serlngamma0(Y,prec),prec));
  if (Y != y)
  {
    GEN pi = mppi(prec);
    z = gdiv(mpodd(y0)? pi: negr(pi),
             gmul(z, gsin(gmul(pi,serchop0(y)), prec)));
  }
  return z;
}

static GEN
sqrtu(ulong a, long prec) { return sqrtr_abs(utor(a, prec)); }
static GEN
cbrtu(ulong a, long prec) { return sqrtnr_abs(utor(a, prec), 3); }
/* N | 6 */
static GEN
ellkprime(long N, GEN s2, GEN s3)
{
  GEN z;
  switch(N)
  {
    case 1: return shiftr(s2, -1);
    case 2: return sqrtr_abs(shiftr(subrs(s2,1), 1));
    case 3: return shiftr(mulrr(s2, addrs(s3, 1)), -2);
    default: /* 6 */
      z = mulrr(subrr(s3,s2), subsr(2,s3));
      return mulrr(addsr(2,s2), sqrtr_abs(z));
  }
}

static GEN
ellKk(long N, GEN s2, GEN s3, long prec)
{ return gdiv(Pi2n(-1,prec), agm(ellkprime(N,s2,s3), gen_1, prec)); }

/* Gamma(1/3) */
static GEN
G3(GEN s2, GEN s3, long prec)
{
  GEN A = ellKk(3, s2,s3, prec), pi = mppi(prec);
  A = shiftr(divrs(powrs(mulrr(pi, A), 12), 27), 28);
  return sqrtnr_abs(A, 36);
}
/* Gamma(1/4) */
static GEN
G4(GEN s2, long prec)
{
  GEN A = ellKk(1, s2,NULL, prec), pi = mppi(prec);
  return shiftr(sqrtr_abs(mulrr(sqrtr_abs(pi), A)), 1);
}

/* Gamma(n / 24), n = 1,5,7,11 */
static GEN
Gn24(long n, GEN s2, GEN s3, long prec)
{
  GEN A, B, C, t, t1, t2, t3, t4, pi = mppi(prec);
  A = ellKk(1, s2,s3, prec);
  B = ellKk(3, s2,s3, prec);
  C = ellKk(6, s2,s3, prec);
  t1 = sqrtr_abs(mulur(3, addsr(2, s3)));
  t2 = sqrtr_abs(divrr(s3, pi));
  t2 = mulrr(t2, shiftr(mulrr(addrr(s2,s3), A), 2));
  t3 = mulrr(divur(3,pi), sqrr(B));
  t3 = mulrr(addsr(2,s2), sqrtnr_abs(shiftr(powrs(t3, 3), 8), 9));
  t4 = mulrr(mulrr(addsr(1, s2), subrr(s3, s2)), subsr(2, s3));
  t4 = mulrr(mulrr(mulur(384, t4), pi), sqrr(C));
  switch (n)
  {
    case 1: t = mulrr(mulrr(t1, t2), mulrr(t3, t4)); break;
    case 5: t = divrr(mulrr(t2, t4), mulrr(t1, t3)); break;
    case 7: t = divrr(mulrr(t3, t4), mulrr(t1, t2)); break;
    default:t = divrr(mulrr(t1, t4), mulrr(t2, t3)); break;
  }
  return sqrtnr_abs(t, 4);
}
/* sin(x/2) = sqrt((1-c) / 2) > 0 given c = cos(x) */
static GEN
sinx2(GEN c)
{ c = subsr(1, c); shiftr_inplace(c,-1); return sqrtr_abs(c); }
/* sin(Pi/12), given sqrt(3) */
static GEN
sin12(GEN s3)
{ GEN t = subsr(2, s3); shiftr_inplace(t, -2); return sqrtr_abs(t); }
/* cos(Pi/12) = sin(5Pi/12), given sqrt(3) */
static GEN
cos12(GEN s3)
{ GEN t = addsr(2, s3); shiftr_inplace(t, -2); return sqrtr_abs(t); }
/* 0 < n < d, (n, d) = 1, 2 < d | 24 */
static GEN
gammafrac24_s(long n, long d, long prec)
{
  GEN A, B, s2, s3, pi = mppi(prec);
  s2 = sqrtu(2, prec);
  s3 = d % 3? NULL: sqrtu(3, prec);
  switch(d)
  {
    case 3:
      A = G3(s2,s3,prec);
      if (n == 1) return A;
      return divrr(Pi2n(1, prec), mulrr(s3, A));
    case 4:
      A = G4(s2,prec);
      if (n == 1) return A;
      return divrr(mulrr(pi, s2), A);
    case 6:
      A = sqrr(G3(s2,s3,prec));
      A = mulrr(A, sqrtr_abs(divsr(3, pi)));
      A = divrr(A, cbrtu(2, prec));
      if (n == 1) return A;
      return divrr(Pi2n(1, prec), A);
    case 8:
      A = ellKk(1, s2,s3, prec);
      B = ellKk(2, s2,s3, prec);
      A = shiftr(sqrtr_abs(divrr(mulrr(addsr(1, s2), A), sqrtr_abs(pi))), 1);
      B = shiftr(mulrr(sqrtr_abs(gmul(subrs(s2, 1), mulrr(s2, pi))), B), 3);
      switch (n)
      {
        GEN t;
        case 1: return sqrtr_abs(mulrr(A, B));
        case 3: return sqrtr_abs(divrr(B, A));
        case 5: A = sqrtr_abs(divrr(B, A));
          t = sqrtr_abs(shiftr(addsr(1, shiftr(s2, -1)), -1)); /*sin(3Pi/8)*/
          return divrr(pi, mulrr(t, A));
        default: A = sqrtr_abs(mulrr(A, B));
          t = sqrtr_abs(shiftr(subsr(1, shiftr(s2, -1)), -1)); /*sin(Pi/8)*/
          return divrr(pi, mulrr(t, A));
      }
    case 12:
      A = G3(s2,s3,prec);
      B = G4(s2,prec);
      switch (n)
      {
        GEN t2;
        case 1: case 11:
          t2 = shiftr(mulur(27, powrs(divrr(addsr(1,s3), pi), 4)), -2);
          t2 = mulrr(sqrtnr_abs(t2, 8), mulrr(A, B));
          if (n == 1) return t2;
          return divrr(pi, mulrr(sin12(s3), t2));
        case 5: case 7:
          t2 = shiftr(divrs(powrs(mulrr(subrs(s3,1), pi), 4), 3), 2);
          t2 = mulrr(sqrtnr_abs(t2, 8), divrr(B, A));
          if (n == 5) return t2;
          return divrr(pi, mulrr(cos12(s3), t2));
      }
    default: /* n = 24 */
      if (n > 12)
      {
        GEN t;
        n = 24 - n;
        A = Gn24(n, s2,s3, prec);
        switch(n)
        { /* t = sin(n*Pi/24) */
          case 1: t = cos12(s3); t = sinx2(t); break;
          case 5: t = sin12(s3); t = sinx2(t); break;
          case 7: t = sin12(s3); togglesign(t); t = sinx2(t); break;
          default:t = cos12(s3); togglesign(t); t = sinx2(t); break; /* n=11 */
        }
        return divrr(pi, mulrr(A, t));
      }
      return Gn24(n, s2,s3, prec);
  }
}

/* (a,b) = 1. If 0 < x < b, m >= 0
gamma(x/b + m) = gamma(x/b) * mulu_interval_step(x, x+(m-1)*b, b) / b^m
gamma(x/b - m) = gamma(x/b) / mulu_interval_step(b-x, b*m-x, b) * (-b)^m */
static GEN
gammafrac24(GEN a, GEN b, long prec)
{
  pari_sp av;
  long A, B, m, x, bit;
  GEN z0, z, t;
  if (!(A = itos_or_0(a)) || !(B = itos_or_0(b)) || B > 24) return NULL;
  switch(B)
  {
    case 2: return gammahs(A-1, prec);
    case 3: case 4: case 6: case 8: case 12: case 24:
      m = A / B;
      x = A % B; /* = A - m*B */
      if (x < 0) { x += B; m--; } /* now 0 < x < B, A/B = x/B + m */
      bit = bit_accuracy(prec);
      /* Depending on B and prec, we must experimentally replace the 0.5
       * by 0.4 to 2.0 for optimal value. Play safe. */
      if (labs(m) > 0.5 * bit * sqrt(bit) / log(bit)) return NULL;
      z0 = cgetr(prec); av = avma;
      prec += EXTRAPREC64;
      z = gammafrac24_s(x, B, prec);
      if (m)
      {
        if (m > 0)
          t = mpdiv(mulu_interval_step(x, (m-1)*B + x, B), rpowuu(B,m,prec));
        else
        {
          m = -m;
          t = mpdiv(rpowuu(B,m,prec), mulu_interval_step(B-x, m*B - x, B));
          if (odd(m)) togglesign(t);
        }
        z = mpmul(z,t);
      }
      affrr(z, z0); set_avma(av); return z0;
  }
  return NULL;
}
GEN
ggamma(GEN x, long prec)
{
  pari_sp av;
  GEN y;

  switch(typ(x))
  {
    case t_INT:
      if (signe(x) <= 0)
        pari_err_DOMAIN("gamma","argument", "=",
                         strtoGENstr("nonpositive integer"), x);
      return mpfactr(itos(x) - 1, prec);

    case t_REAL: case t_COMPLEX:
      return cxgamma(x, 0, prec);

    case t_FRAC:
    {
      GEN a = gel(x,1), b = gel(x,2), c = gammafrac24(a, b, prec);
      if (c) return c;
      av = avma; c = subii(a,b);
      if (signe(a) < 0)
      { /* gamma will use functional equation x -> z = 1-x = -c/b >= 1/2.
         * Gamma(x) = Pi / (sin(Pi z) * Gamma(z)) */
        GEN z = mkfrac(negi(c), b), q = ground(z), r = gsub(z,q);
        GEN pi = mppi(prec); /* |r| <= 1/2 */
        z = fractor(z, prec+EXTRAPREC64);
        y = divrr(pi, mulrr(mpsin(gmul(pi, r)), cxgamma(z, 0, prec)));
        if (mpodd(q)) togglesign(y);
        return gerepileupto(av, y);
      }
      if (cmpii(shifti(a,1), b) < 0)
      { /* 0 < x < 1/2 gamma would use funeq: adding 1 is cheaper. */
        if (expi(a) - expi(b) < -3) /* close to 0 */
        {
          if (lgefint(b) >= prec) x = fractor(x,prec);
          y = mpexp(lngamma1(x, prec));
        }
        else
          y = cxgamma(fractor(mkfrac(addii(a,b), b), prec), 0, prec);
        return gerepileupto(av, gdiv(y, x));
      }
      if (expi(c) - expi(b) < -3)
      { /* x = 1 + c/b is close to 1 */
        x = mkfrac(c,b);
        if (lgefint(b) >= prec) x = fractor(x,prec);
        y = mpexp(lngamma1(x, prec));
      }
      else
        y = cxgamma(fractor(x, prec), 0, prec);
      return gerepileupto(av, y);
    }

    case t_PADIC: return Qp_gamma(x);
    default:
      av = avma; if (!(y = toser_i(x))) break;
      return gerepileupto(av, sergamma(y, prec));
  }
  return trans_eval("gamma",ggamma,x,prec);
}

static GEN
mpfactr_basecase(long n, long prec)
{
  GEN v = cgetg(expu(n) + 2, t_VEC);
  long k, prec2 = prec + EXTRAPREC64;
  GEN a;
  for (k = 1;; k++)
  {
    long m = n >> (k-1), l;
    if (m <= 2) break;
    l = (1 + (n >> k)) | 1;
    /* product of odd numbers in ]n / 2^k, 2 / 2^(k-1)] */
    a = mulu_interval_step_prec(l, m, 2, prec2);
    gel(v,k) = k == 1? a: gpowgs(a, k);
  }
  a = gel(v,--k); while (--k) a = mpmul(a, gel(v,k));
  if (typ(a) == t_INT) a = itor(a, prec); else a = gprec_wtrunc(a, prec);
  shiftr_inplace(a, factorial_lval(n, 2));
  return a;
}
/* Theory says n > C * b^1.5 / log(b). Timings:
 * b = [64, 128, 192, 256, 512, 1024, 2048, 4096, 8192, 16384]
 * n = [1930, 2650, 3300, 4270, 9000, 23000, 75000, 210000, 750000, 2400000] */
static long
mpfactr_n(long prec)
{
  long b = bit_accuracy(prec);
  if (b <=  64) return 1930;
  if (b <= 128) return 2650;
  if (b <= 192) return 3300;
  return b * sqrt(b);
}
static GEN
mpfactr_small(long n, long prec)
{
  GEN f = cgetr(prec);
  pari_sp av = avma;
  if (n < 410)
    affir(mpfact(n), f);
  else
    affrr(mpfactr_basecase(n, prec), f);
  set_avma(av); return f;
}
GEN
mpfactr(long n, long prec)
{
  GEN f = cgetr(prec);
  pari_sp av = avma;

  if (n < 410)
    affir(mpfact(n), f);
  else
  {
    long N = mpfactr_n(prec);
    GEN z = n <= N? mpfactr_basecase(n, prec)
                  : cxgamma(utor(n+1, prec), 0, prec);
    affrr(z, f);
  }
  set_avma(av); return f;
}

/* First a little worse than mpfactr_n because of the extra logarithm.
 * Asymptotically same. */
static ulong
lngamma_n(long prec)
{
  long b = bit_accuracy(prec);
  double N;
  if (b <=  64) return 1450UL;
  if (b <= 128) return 2010UL;
  if (b <= 192) return 2870UL;
  N = b * sqrt(b);
  if (b <= 256) return N/1.25;
  if (b <= 512) return N/1.2;
  if (b <= 2048) return N/1.1;
  return N;
}

GEN
glngamma(GEN x, long prec)
{
  pari_sp av = avma;
  GEN y, y0, t;

  switch(typ(x))
  {
    case t_INT:
    {
      ulong n;
      if (signe(x) <= 0)
        pari_err_DOMAIN("lngamma","argument", "=",
                         strtoGENstr("nonpositive integer"), x);
      n = itou_or_0(x);
      if (!n || n > lngamma_n(prec)) return cxgamma(x, 1, prec);
      return gerepileuptoleaf(av, logr_abs( mpfactr_small(n-1, prec) ));
    }
    case t_FRAC:
    {
      GEN a = gel(x,1), b = gel(x,2), c = gammafrac24(a, b, prec);
      long e;
      if (c) return glog(c, prec);
      c = subii(a,b); e = expi(b) - expi(c);
      if (signe(a) < 0)
      { /* gamma will use functional equation x -> z = 1-x = -c/b >= 1/2.
         * lngamma(x) = log |Pi / (sin(Pi z) * Gamma(z))| + I*Pi * floor(x) */
        GEN z = mkfrac(negi(c), b), q = ground(z), r = gsub(z,q);
        GEN pi = mppi(prec); /* |r| <= 1/2 */
        z = fractor(z, prec+EXTRAPREC64);
        y = subrr(logr_abs(divrr(pi, mpsin(gmul(pi,r)))), cxgamma(z, 1, prec));
        y = gadd(y, mkcomplex(gen_0, mulri(pi, gfloor(x))));
        return gerepileupto(av, y);
      }
      if (cmpii(shifti(a,1), b) < 0)
      { /* 0 < x < 1/2 gamma would use funeq: adding 1 is cheaper. */
        if (expi(a) - expi(b) < -3) /* close to 0 */
        {
          if (lgefint(b) >= prec) x = fractor(x,prec);
          y = lngamma1(x, prec);
        }
        else
          y = cxgamma(fractor(mkfrac(addii(a,b), b), prec), 1, prec);
        return gerepileupto(av, gsub(y, glog(x, prec)));
      }
      if (e > 3)
      {
        x = mkfrac(c,b);
        if (lgefint(b) >= prec) x = fractor(x,prec + nbits2nlong(e));
        y = lngamma1(x, prec);
      }
      else
      {
        x = fractor(x, e > 1? prec+EXTRAPREC64: prec);
        y = cxgamma(x, 1, prec);
      }
      return gerepileupto(av, y);
    }

    case t_REAL: case t_COMPLEX:
      return cxgamma(x, 1, prec);

    default:
      if (!(y = toser_i(x))) break;
      if (lg(y) == 2) pari_err_DOMAIN("lngamma", "argument", "=", gen_0,y);
      t = serlngamma0(y,prec);
      y0 = simplify_shallow(gel(y,2));
      /* no constant term if y0 = 1 or 2 */
      if (!isint(y0,&y0) || signe(y0) <= 0 || abscmpiu(y0,2) > 2)
        t = gadd(t, glngamma(y0,prec));
      return gerepileupto(av, t);

    case t_PADIC: return gerepileupto(av, Qp_lngamma(x));
  }
  return trans_eval("lngamma",glngamma,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                  PSI(x) = GAMMA'(x)/GAMMA(x)                   **/
/**                                                                **/
/********************************************************************/
static void
err_psi(GEN s)
{
  pari_err_DOMAIN("psi","argument", "=",
                  strtoGENstr("nonpositive integer"), s);
}
/* L ~ |log s|^2 */
static long
psi_lim(double L, double la, long prec)
{
  double d = (prec2nbits_mul(prec, 2*M_LN2) - log(L)) / (4*(1+log(la)));
  return (d < 2)? 2: 2 + (long)ceil(d);
}
/* max(|log (s + it - Euler)|, 1e-6) */
static double
dlogE(double s, double t)
{
  double rlog, ilog;
  dcxlog(s - 0.57721566, t, &rlog,&ilog);
  return maxdd(dnorm(rlog,ilog), 1e-6);
}
static GEN
cxpsi(GEN s0, long prec)
{
  pari_sp av, av2;
  GEN sum, z, a, res, sig, tau, s, unr, s2, sq;
  long lim, nn, k;
  const long la = 3;
  int funeq = 0;
  pari_timer T;

  if (DEBUGLEVEL>2) timer_start(&T);
  s = trans_fix_arg(&prec,&s0,&sig,&tau,&av,&res);
  if (signe(sig) <= 0) { funeq = 1; s = gsub(gen_1, s); sig = real_i(s); }
  if (typ(s0) == t_INT && signe(s0) <= 0) err_psi(s0);

  if (expo(sig) > 300 || (typ(s) == t_COMPLEX && gexpo(gel(s,2)) > 300))
  { /* |s| is HUGE. Play safe */
    GEN L, S = gprec_w(s,LOWDEFAULTPREC), rS = real_i(S), iS = imag_i(S);
    double l;
    lim = psi_lim(rtodbl(gnorm(glog(S,LOWDEFAULTPREC))), la, prec);
    l = (2*lim-1)*la / (2.*M_PI);
    L = gsub(dbltor(l*l), gsqr(iS));
    if (signe(L) < 0) L = gen_0;
    L = gsub(gsqrt(L, LOWDEFAULTPREC), rS);
    if (signe(L) > 0) nn = (long)ceil(rtodbl(L)); else nn = 1;
  }
  else
  {
    double l, rS = rtodbl(sig), iS = typ(s) == t_REAL? 0.0: rtodbl(imag_i(s));
    lim = psi_lim(dlogE(rS, iS), la, prec);
    l = (2*lim-1)*la / (2.*M_PI);
    l = l*l - iS*iS;
    if (l < 0.) l = 0.;
    nn = (long)ceil( sqrt(l) - rS );
    if (nn < 1) nn = 1;
  }
  if (DEBUGLEVEL>2) err_printf("lim, nn: [%ld, %ld]\n",lim,nn);
  incrprec(prec); unr = real_1(prec); /* one extra word of precision */
  s2 = gmul2n(s, 1); sq = gsqr(s);
  a = gdiv(unr, gaddgs(s, nn)); /* 1 / (s+n) */
  av2 = avma; sum = gmul2n(a, -1);
  for (k = 0; k < nn - 1; k += 2)
  {
    GEN tmp = gaddsg(k*(k + 1), gadd(gmulsg(2*k + 1, s), sq));
    sum = gadd(sum, gdiv(gaddsg(2*k + 1, s2), tmp));
    if ((k & 1023) == 0) sum = gerepileupto(av2, sum);
  }
  if (odd(nn)) sum = gadd(sum, gdiv(unr, gaddsg(nn - 1, s)));
  z = gsub(glog(gaddgs(s, nn), prec), sum);
  if (DEBUGLEVEL>2) timer_printf(&T,"sum from 0 to N - 1");
  constbern(lim);
  z = gsub(z, psi_sum(gsqr(a), lim));
  if (DEBUGLEVEL>2) timer_printf(&T,"Bernoulli sum");
  if (funeq)
  {
    GEN pi = mppi(prec);
    z = gadd(z, gmul(pi, gcotan(gmul(pi,s), prec)));
  }
  set_avma(av); return affc_fixlg(z, res);
}

/* n >= 0; return psi(1+x) + O(x^n), x = pol_x(v) */
GEN
psi1series(long n, long v, long prec)
{
  long i, l = n+3;
  GEN s = cgetg(l, t_SER), z = constzeta(n + 1, prec);

  s[1] = evalsigne(1)|evalvalser(0)|evalvarn(v);
  for (i = 1; i <= n+1; i++)
  {
    GEN c = gel(z,i); /* zeta(i) */
    gel(s,i+1) = odd(i)? negr(c): c;
  }
  return s;
}
/* T an RgX, return T(X + z0) + O(X^L) */
static GEN
tr(GEN T, GEN z0, long L)
{
  GEN s = RgX_to_ser(RgX_translate(T, z0), L+3);
  setvarn(s, 0); return s;
}
/* z0 a complex number with Re(z0) > 1/2; return psi(z0+x) + O(x^L)
 * using Luke's rational approximation for psi(x) */
static GEN
serpsiz0(GEN z0, long L, long v, long prec)
{
  pari_sp av;
  GEN A,A1,A2, B,B1,B2, Q;
  long n;
  n = gprecision(z0); if (n) prec = n;
  z0 = gtofp(z0, prec + EXTRAPREC64);
  /* Start from n = 3; in Luke's notation, A2 := A_{n-2}, A1 := A_{n-1},
   * A := A_n. Same for B */
  av = avma;
  A2= gdivgu(mkpoln(2, gen_1, utoipos(6)), 2);
  B2 = scalarpol_shallow(utoipos(4), 0);
  A1= gdivgu(mkpoln(3, gen_1, utoipos(82), utoipos(96)), 6);
  B1 = mkpoln(2, utoipos(8), utoipos(28));
  A = gdivgu(mkpoln(4, gen_1, utoipos(387), utoipos(2906), utoipos(1920)), 12);
  B = mkpoln(3, utoipos(14), utoipos(204), utoipos(310));
  A2= tr(A2,z0, L);
  B2= tr(B2,z0, L);
  A1= tr(A1,z0, L);
  B1= tr(B1,z0, L);
  A = tr(A, z0, L);
  B = tr(B, z0, L); Q = gdiv(A, B);
  /* work with z0+x as a variable */
  for (n = 4;; n++)
  {
    GEN Q0 = Q, a, b, r, c3,c2,c1,c0 = muluu(2*n-3, n+1);
    GEN u = subiu(muluu(n, 7*n-9), 6);
    GEN t = addiu(muluu(n, 7*n-19), 4);
    /* c1=(2*n-1)*(3*(n-1)*z+7*n^2-9*n-6);
     * c2=(2*n-3)*(z-n-1)*(-3*(n-1)*z+7*n^2-19*n+4);
     * c3=(2*n-1)*(n-3)*(z-n)*(z-(n+1))*(z+(n-4)); */
    c1 = deg1pol_shallow(muluu(3*(n-1),2*n-1), muliu(u,2*n-1), 0);
    c2 = ZX_mul(deg1pol_shallow(utoipos(2*n-3), negi(muluu(2*n-3,n+1)), 0),
                deg1pol_shallow(utoineg(3*(n-1)), t, 0));
    r = mkvec3(utoipos(n), utoipos(n+1), stoi(4-n));
    c3 = ZX_Z_mul(roots_to_pol(r,0), muluu(2*n-1,n-3));
    c1 = tr(c1, z0, L+3);
    c2 = tr(c2, z0, L+3);
    c3 = tr(c3, z0, L+3);

    /* A_{n+1}, B_{n+1} */
    a = gdiv(gadd(gadd(gmul(c1,A),gmul(c2,A1)),gmul(c3,A2)), c0);
    b = gdiv(gadd(gadd(gmul(c1,B),gmul(c2,B1)),gmul(c3,B2)), c0);
    Q = gdiv(a,b);
    if (gexpo(gsub(Q,Q0)) < -prec2nbits(prec)) break;
    A2 = A1; A1 = A; A = a;
    B2 = B1; B1 = B; B = b;
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"serpsiz0, n = %ld", n);
      gerepileall(av, 7, &A,&A1,&A2, &B,&B1,&B2, &Q);
    }
  }
  Q = gmul(Q, gmul2n(gsubsg(1, ginv(tr(pol_x(v),z0, L))), 1));
  setvarn(Q, v);
  return gadd(negeuler(prec), Q);
}
/* sum (-1)^k*H(m,k)x^k + O(x^L); L > 0;
 * H(m,k) = (-1)^{k * \delta_{m > 0}} sum_{1<=i<m} 1/i^(k+1) */
static GEN
Hseries(long m, long L, long v, long prec)
{
  long i, k, bit, l = L+3, M = m < 0? 1-m: m;
  pari_sp av = avma;
  GEN H = cgetg(l, t_SER);
  H[1] = evalsigne(1)|evalvarn(v)|evalvalser(0);
  prec += EXTRAPREC64;
  bit = -prec2nbits(prec);
  for(k = 2; k < l; k++) gel(H,k) = gen_1; /* i=1 */
  for (i = 2; i < M; i++)
  {
    GEN ik = invr(utor(i, prec));
    for (k = 2; k < l; k++)
    {
      if (k > 2) { ik = divru(ik, i); if (expo(ik) < bit) break; }
      gel(H,k) = gadd(gel(H,k), ik);
    }
    if (gc_needed(av,3))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Hseries, i = %ld/%ld", i,M);
      H = gerepilecopy(av, H);
    }
  }
  if (m > 0)
    for (k = 3; k < l; k+=2) togglesign_safe(&gel(H,k));
  return H;
}

static GEN
serpsi(GEN y, long prec)
{
  GEN Q = NULL, z0, Y = y, Y2;
  long L = lg(y)-2, v  = varn(y), vy = valser(y);

  if (!L) pari_err_DOMAIN("psi", "argument", "=", gen_0,y);
  if (vy < 0) pari_err_DOMAIN("psi", "series valuation", "<", gen_0,y);
  if (vy)
    z0 = gen_0;
  else
  {
    z0 = simplify_shallow(gel(y,2));
    (void)isint(z0, &z0);
  }
  if (typ(z0) == t_INT && !is_bigint(z0))
  {
    long m = itos(z0);
    if (abscmpiu(muluu(prec2nbits(prec),L), labs(m)) > 0)
    { /* psi(m+x) = psi(1+x) + sum_{1 <= i < m} 1/(i+x) for m > 0
                    psi(1+x) - sum_{0 <= i < -m} 1/(i+x) for m <= 0 */
      GEN H = NULL;
      if (m <= 0) L--; /* lose series accuracy due to 1/x term */
      if (L)
      {
        Q = psi1series(L, v, prec);
        if (m && m != 1) { H = Hseries(m, L, v, prec); Q = gadd(Q, H); }
        if (m <= 0) Q = gsub(Q, ginv(pol_x(v)));
      }
      else
      {
        Q = scalarser(gen_m1, v, 1);
        setvalser(Q,-1);
      }
    }
  }
  if (!Q)
  { /* use psi(1-y)=psi(y)+Pi*cotan(Pi*y) ? */
    if (gcmp(real_i(z0),ghalf) < 0) { z0 = gsubsg(1,z0); Y = gsubsg(1,y); }
    Q = serpsiz0(z0, L, v, prec);
  }
  Y2 = serchop0(Y); if (signe(Y2)) Q = gsubst(Q, v, Y2);
  /* psi(z0 + Y2) = psi(Y) */
  if (Y != y)
  { /* psi(y) = psi(Y) + Pi cotan(Pi Y) */
    GEN pi = mppi(prec);
    if (typ(z0) == t_INT) Y = Y2; /* in this case cotan(Pi*Y2) = cotan(Pi*Y) */
    Q = gadd(Q, gmul(pi, gcotan(gmul(pi,Y), prec)));
  }
  return Q;
}

static ulong
psi_n(ulong b)
{
  if (b <= 64) return 50;
  if (b <= 128) return 85;
  if (b <= 192) return 122;
  if (b <= 256) return 150;
  if (b <= 512) return 320;
  if (b <= 1024) return 715;
  return 0.010709 * pow((double)b, 1.631); /* 1.631 ~ log_3(6) */
}
GEN
gpsi(GEN x, long prec)
{
  pari_sp av;
  ulong n;
  GEN y;
  switch(typ(x))
  {
    case t_INT:
      if (signe(x) <= 0) err_psi(x);
      if (lgefint(x) > 3 || (n = itou(x)) > psi_n(prec2nbits(prec))) break;
      av = avma; y = mpeuler(prec);
      return gerepileuptoleaf(av, n == 1? negr(y): gsub(harmonic(n-1), y));
    case t_REAL: case t_COMPLEX: return cxpsi(x,prec);
    default:
      av = avma; if (!(y = toser_i(x))) break;
      return gerepileupto(av, serpsi(y,prec));
  }
  return trans_eval("psi",gpsi,x,prec);
}
