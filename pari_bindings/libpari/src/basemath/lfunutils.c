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
/**                 L-functions: Applications                      **/
/**                                                                **/
/********************************************************************/

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_lfun

static GEN
tag(GEN x, long t) { return mkvec2(mkvecsmall(t), x); }

/* v a t_VEC of length > 1 */
static int
is_tagged(GEN v)
{
  GEN T = gel(v,1);
  return (typ(T)==t_VEC && lg(T)==3 && typ(gel(T,1))==t_VECSMALL);
}
/* rough check */
static long
is_ldata(GEN L)
{
  long l = lg(L);
  return typ(L) == t_VEC && (l == 7 || l == 8);
}
/* thorough check */
static void
checkldata(GEN ldata)
{
  GEN vga, w, N;
#if 0 /* assumed already checked and true */
  if (!is_ldata(ldata) || !is_tagged(ldata)) pari_err_TYPE("checkldata", ldata);
#endif
  vga = ldata_get_gammavec(ldata);
  if (typ(vga) != t_VEC) pari_err_TYPE("checkldata [gammavec]",vga);
  w = gel(ldata, 4); /* FIXME */
  switch(typ(w))
  {
    case t_INT: case t_FRAC: break;
    case t_VEC: if (lg(w) == 3 && is_rational_t(typ(gel(w,1)))) break;
    default: pari_err_TYPE("checkldata [weight]",w);
  }
  N = ldata_get_conductor(ldata);
  if (typ(N) != t_INT) pari_err_TYPE("checkldata [conductor]",N);
}

/* tag as t_LFUN_GENERIC */
static void
lfuncreate_tag(GEN L)
{
  if (is_tagged(L)) return;
  gel(L,1) = tag(gel(L,1), t_LFUN_GENERIC);
  if (typ(gel(L,2)) != t_INT) gel(L,2) = tag(gel(L,2), t_LFUN_GENERIC);
}

/* shallow */
static GEN
closure2ldata(GEN C, long prec)
{
  GEN L = closure_callgen0prec(C, prec);
  if (is_ldata(L)) { checkldata(L); lfuncreate_tag(L); }
  else L = lfunmisc_to_ldata_shallow(L);
  return L;
}

/* data may be either an object (polynomial, elliptic curve, etc...)
 * or a description vector [an,sd,Vga,k,conductor,rootno,{poles}]. */
GEN
lfuncreate(GEN data)
{
  if (is_ldata(data))
  {
    GEN L = gcopy(data);
    lfuncreate_tag(L); checkldata(L); return L;
  }
  if (typ(data) == t_CLOSURE && closure_arity(data)==0)
  {
    pari_sp av = avma;
    GEN L = closure2ldata(data, DEFAULTPREC);
    gel(L,1) = tag(data, t_LFUN_CLOSURE0); return gerepilecopy(av, L);
  }
  return lfunmisc_to_ldata(data);
}

GEN
lfunparams(GEN L, long prec)
{
  pari_sp av = avma;
  GEN k, N, v;
  long p;

  if (!is_ldata(L) || !is_tagged(L)) L = lfunmisc_to_ldata_shallow(L);
  N = ldata_get_conductor(L);
  k = ldata_get_k(L);
  v = ldata_get_gammavec(L);
  p = gprecision(v);
  if (p > prec) v = gprec_wtrunc(v, prec);
  else if (p < prec)
  {
    GEN van = ldata_get_an(L), an = gel(van,2);
    long t = mael(van,1,1);
    if (t == t_LFUN_CLOSURE0) L = closure2ldata(an, prec);
  }
  return gerepilecopy(av, mkvec3(N, k, v));
}

/********************************************************************/
/**                     Simple constructors                        **/
/********************************************************************/
static GEN ldata_eulerf(GEN van, GEN p, long prec);

static GEN
vecan_conj(GEN an, long n, long prec)
{
  GEN p1 = ldata_vecan(gel(an,1), n, prec);
  return typ(p1) == t_VEC? conj_i(p1): p1;
}

static GEN
eulerf_conj(GEN an, GEN p, long prec)
{
  GEN p1 = ldata_eulerf(gel(an,1), p, prec);
  return conj_i(p1);
}

static GEN
vecan_mul(GEN an, long n, long prec)
{
  GEN p1 = ldata_vecan(gel(an,1), n, prec);
  GEN p2 = ldata_vecan(gel(an,2), n, prec);
  if (typ(p1) == t_VECSMALL) p1 = vecsmall_to_vec(p1);
  if (typ(p2) == t_VECSMALL) p2 = vecsmall_to_vec(p2);
  return dirmul(p1, p2);
}

static GEN
eulerf_mul(GEN an, GEN p, long prec)
{
  GEN p1 = ldata_eulerf(gel(an,1), p, prec);
  GEN p2 = ldata_eulerf(gel(an,2), p, prec);
  return gmul(p1, p2);
}

static GEN
lfunconvol(GEN a1, GEN a2)
{ return tag(mkvec2(a1, a2), t_LFUN_MUL); }

static GEN
vecan_div(GEN an, long n, long prec)
{
  GEN p1 = ldata_vecan(gel(an,1), n, prec);
  GEN p2 = ldata_vecan(gel(an,2), n, prec);
  if (typ(p1) == t_VECSMALL) p1 = vecsmall_to_vec(p1);
  if (typ(p2) == t_VECSMALL) p2 = vecsmall_to_vec(p2);
  return dirdiv(p1, p2);
}

static GEN
eulerf_div(GEN an, GEN p, long prec)
{
  GEN p1 = ldata_eulerf(gel(an,1), p, prec);
  GEN p2 = ldata_eulerf(gel(an,2), p, prec);
  return gdiv(p1, p2);
}

static GEN
lfunconvolinv(GEN a1, GEN a2)
{ return tag(mkvec2(a1,a2), t_LFUN_DIV); }

static GEN
lfunconj(GEN a1)
{ return tag(mkvec(a1), t_LFUN_CONJ); }

static GEN
lfuncombdual(GEN (*fun)(GEN, GEN), GEN ldata1, GEN ldata2)
{
  GEN a1 = ldata_get_an(ldata1), a2 = ldata_get_an(ldata2);
  GEN b1 = ldata_get_dual(ldata1), b2 = ldata_get_dual(ldata2);
  if (typ(b1)==t_INT && typ(b2)==t_INT)
    return utoi(signe(b1) || signe(b2));
  else
  {
    if (typ(b1)==t_INT) b1 = signe(b1) ? lfunconj(a1): a1;
    if (typ(b2)==t_INT) b2 = signe(b2) ? lfunconj(a2): a2;
    return fun(b1, b2);
  }
}

static GEN
vecan_twist(GEN an, long n, long prec)
{
  GEN p1 = ldata_vecan(gel(an,1), n, prec);
  GEN p2 = ldata_vecan(gel(an,2), n, prec);
  long i;
  GEN V;
  if (typ(p1) == t_VECSMALL) p1 = vecsmall_to_vec(p1);
  if (typ(p2) == t_VECSMALL) p2 = vecsmall_to_vec(p2);
  V = cgetg(n+1, t_VEC);
  for(i = 1; i <= n ; i++)
    gel(V, i) = gmul(gel(p1, i), gel(p2, i));
  return V;
}

static GEN
eulerf_twist(GEN an, GEN p, long prec)
{
  GEN p1 = ldata_eulerf(gel(an,1), p, prec);
  GEN p2 = ginv(ldata_eulerf(gel(an,2), p, prec));
  if (typ(p2)!=t_POL || degpol(p2)==0)
    return poleval(p1,pol_0(0));
  if (degpol(p2)!=1) pari_err_IMPL("lfuneuler");
  return poleval(p1,monomial(gneg(gel(p2,3)),1,0));
}

static GEN
vecan_shift(GEN an, long n, long prec)
{
  GEN p1 = ldata_vecan(gel(an,1), n, prec);
  GEN s = gel(an,2);
  long i;
  GEN V;
  if (typ(p1) == t_VECSMALL) p1 = vecsmall_to_vec(p1);
  V = cgetg(n+1, t_VEC);
  if (typ(s)==t_INT)
  {
    if (equali1(s))
      for(i = 1; i <= n ; i++)
      {
        GEN gi = gel(p1, i);
        gel(V, i) = gequal0(gi)? gi: gmulgu(gi, i);
      }
    else
      for(i = 1; i <= n ; i++)
      {
        GEN gi = gel(p1, i);
        gel(V, i) = gequal0(gi)? gi: gmul(gi, powgi(utoi(i), s));
      }
  }
  else
  {
    GEN D = dirpowers(n, s, prec);
    for(i = 1; i <= n ; i++)
      gel(V, i) = gmul(gel(p1,i), gel(D,i));
  }
  return V;
}

static GEN
eulerf_shift(GEN an, GEN p, long prec)
{
  GEN p1 = ldata_eulerf(gel(an,1), p, prec);
  GEN s = gel(an,2);
  return gsubst(p1, 0, monomial(gpow(p, s, prec), 1, 0));
}

static GEN
eulerf_hgm(GEN an, GEN p)
{
  GEN H = gel(an,1), t = gel(an,2);
  if (typ(t)==t_VEC && lg(t)==3)
  {
    GEN L = gel(t,2);
    long i, l = lg(L);
    t = gel(t,1);
    for (i = 1; i < l; i++) /* wild primes */
      if (equalii(p, gmael(L, i, 1))) break;
    if (i<l) return gmael(L,i,2);
  }
  return ginv(hgmeulerfactor(H, t, itos(p), NULL));
}

static GEN
deg1ser_shallow(GEN a1, GEN a0, long e)
{ return RgX_to_ser(deg1pol_shallow(a1, a0, 0), e+2); }
/* lfunrtopoles without sort */
static GEN
rtopoles(GEN r)
{
  long j, l = lg(r);
  GEN v = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
  {
    GEN rj = gel(r,j), a = gel(rj,1);
    gel(v,j) = a;
  }
  return v;
}
/* re = polar part; overestimate when re = gen_0 (unknown) */
static long
orderpole(GEN re) { return typ(re) == t_SER? -valser(re): 1; }
static GEN
lfunmulpoles(GEN ldata1, GEN ldata2, long bitprec)
{
  GEN r, k = ldata_get_k(ldata1), b1 = NULL, b2 = NULL;
  GEN r1 = ldata_get_residue(ldata1);
  GEN r2 = ldata_get_residue(ldata2);
  long i, j, l, L = 0;

  if (!r1 && !r2) return NULL;
  if (r1 && !is_vec_t(typ(r1))) r1 = mkvec(mkvec2(k, r1));
  if (r2 && !is_vec_t(typ(r2))) r2 = mkvec(mkvec2(k, r2));
  if (r1) { b1 = rtopoles(r1); L += lg(b1); }
  if (r2) { b2 = rtopoles(r2); L += lg(b2); }
  r = cgetg(L, t_VEC); j = 1;
  if (b1)
  {
    l = lg(b1);
    for (i = 1; i < l; i++)
    {
      GEN z, z1, z2, be = gmael(r1,i,1);
      long n, v = orderpole(gmael(r1,i,2));
      if (b2 && (n = RgV_isin(b2, be))) v += orderpole(gmael(r2,n,2));
      z = deg1ser_shallow(gen_1, be, 2 + v);
      z1 = lfun(ldata1,z,bitprec);
      z2 = lfun(ldata2,z,bitprec);
      gel(r,j++) = mkvec2(be, gmul(z1, z2));
    }
  }
  if (b2)
  {
    long l = lg(b2);
    for (i = 1; i < l; i++)
    {
      GEN z, z1, z2, be = gmael(r2,i,1);
      long n, v = orderpole(gmael(r2,i,2));
      if (b1 && (n = RgV_isin(b1, be))) continue; /* done already */
      z = deg1ser_shallow(gen_1, be, 2 + v);
      z1 = lfun(ldata1,z,bitprec);
      z2 = lfun(ldata2,z,bitprec);
      gel(r,j++) = mkvec2(be, gmul(z1, z2));
    }
  }
  setlg(r, j); return r;
}

static GEN
lfunmul_k(GEN ldata1, GEN ldata2, GEN k, long bitprec)
{
  GEN r, N, Vga, eno, a1a2, b1b2;
  r = lfunmulpoles(ldata1, ldata2, bitprec);
  N = gmul(ldata_get_conductor(ldata1), ldata_get_conductor(ldata2));
  Vga = shallowconcat(ldata_get_gammavec(ldata1), ldata_get_gammavec(ldata2));
  Vga = sort(Vga);
  eno = gmul(ldata_get_rootno(ldata1), ldata_get_rootno(ldata2));
  a1a2 = lfunconvol(ldata_get_an(ldata1), ldata_get_an(ldata2));
  b1b2 = lfuncombdual(lfunconvol, ldata1, ldata2);
  return mkvecn(r? 7: 6, a1a2, b1b2, Vga, k, N, eno, r);
}

GEN
lfunmul(GEN ldata1, GEN ldata2, long bitprec)
{
  pari_sp ltop = avma;
  GEN k;
  long prec = nbits2prec(bitprec);
  ldata1 = ldata_newprec(lfunmisc_to_ldata_shallow(ldata1), prec);
  ldata2 = ldata_newprec(lfunmisc_to_ldata_shallow(ldata2), prec);
  k = ldata_get_k(ldata1);
  if (!gequal(ldata_get_k(ldata2),k))
    pari_err_OP("lfunmul [weight]",ldata1, ldata2);
  return gerepilecopy(ltop, lfunmul_k(ldata1, ldata2, k, bitprec));
}

static GEN
lfundivpoles(GEN ldata1, GEN ldata2, long bitprec)
{
  long i, j, l;
  GEN be2, k  = ldata_get_k(ldata1);
  GEN r1 = ldata_get_residue(ldata1);
  GEN r2 = ldata_get_residue(ldata2), r;

  if (r1 && !is_vec_t(typ(r1))) r1 = mkvec(mkvec2(k, r1));
  if (r2 && !is_vec_t(typ(r2))) r2 = mkvec(mkvec2(k, r2));
  if (!r1) return NULL;
  l = lg(r1); r = cgetg(l, t_VEC);
  be2 = r2? rtopoles(r2): NULL;
  for (i = j = 1; j < l; j++)
  {
    GEN z, v = gel(r1,j), be = gel(v,1), s1 = gel(v,2);
    long n;
    if (be2 && (n = RgV_isin(be2, be)))
    {
      GEN s2 = gmael(r2,n,2); /* s1,s2: polar parts */
      if (orderpole(s1) == orderpole(s2)) continue;
    }
    z = gdiv(lfun(ldata1,be,bitprec), lfun(ldata2,be,bitprec));
    if (valser(z) < 0) gel(r,i++) = mkvec2(be, z);
  }
  if (i == 1) return NULL;
  setlg(r, i); return r;
}

static GEN
lfunvgasub(GEN v01, GEN v2)
{
  GEN v1 = shallowcopy(v01), v;
  long l1 = lg(v1), l2 = lg(v2), j1, j2, j;
  for (j2 = 1; j2 < l2; j2++)
  {
    for (j1 = 1; j1 < l1; j1++)
      if (gel(v1,j1) && gequal(gel(v1,j1), gel(v2,j2)))
      {
        gel(v1,j1) = NULL; break;
      }
    if (j1 == l1) pari_err_OP("lfunvgasub", v1, v2);
  }
  v = cgetg(l1-l2+1, t_VEC);
  for (j1 = j = 1; j1 < l1; j1++)
    if (gel(v1, j1)) gel(v,j++) = gel(v1,j1);
  return v;
}

GEN
lfundiv(GEN ldata1, GEN ldata2, long bitprec)
{
  pari_sp ltop = avma;
  GEN k, r, N, v, eno, a1a2, b1b2, eno2;
  long prec = nbits2prec(bitprec);
  ldata1 = ldata_newprec(lfunmisc_to_ldata_shallow(ldata1), prec);
  ldata2 = ldata_newprec(lfunmisc_to_ldata_shallow(ldata2), prec);
  k = ldata_get_k(ldata1);
  if (!gequal(ldata_get_k(ldata2),k))
    pari_err_OP("lfundiv [weight]",ldata1, ldata2);
  N = gdiv(ldata_get_conductor(ldata1), ldata_get_conductor(ldata2));
  if (typ(N) != t_INT) pari_err_OP("lfundiv [conductor]",ldata1, ldata2);
  r = lfundivpoles(ldata1, ldata2, bitprec);
  a1a2 = lfunconvolinv(ldata_get_an(ldata1), ldata_get_an(ldata2));
  b1b2 = lfuncombdual(lfunconvolinv, ldata1, ldata2);
  eno2 = ldata_get_rootno(ldata2);
  eno = isintzero(eno2)? gen_0: gdiv(ldata_get_rootno(ldata1), eno2);
  v = lfunvgasub(ldata_get_gammavec(ldata1), ldata_get_gammavec(ldata2));
  return gerepilecopy(ltop,  mkvecn(r? 7: 6, a1a2, b1b2, v, k, N, eno, r));
}

static GEN
gamma_imagchi(GEN gam, GEN w)
{
  long i, j, k=1, l;
  GEN g = cgetg_copy(gam, &l);
  gam = shallowcopy(gam);
  for (i = l-1; i>=1; i--)
  {
    GEN al = gel(gam, i);
    if (al)
    {
      GEN N = gadd(w,gmul2n(real_i(al),1));
      if (gcmpgs(N,2) > 0)
      {
        GEN bl = gsubgs(al, 1);
        for (j=1; j < i; j++)
          if (gel(gam,j) && gequal(gel(gam,j), bl))
          { gel(gam,j) = NULL; break; }
        if (j==i) return NULL;
        gel(g, k++) = al;
        gel(g, k++) = bl;
      } else if (gequal0(N))
        gel(g, k++) = gaddgs(al, 1);
      else if (gequal1(N))
        gel(g, k++) = gsubgs(al, 1);
      else return NULL;
    }
  }
  return sort(g);
}

GEN
lfuntwist(GEN ldata1, GEN chi, long bitprec)
{
  pari_sp ltop = avma;
  GEN k, L, N, N1, N2, a, a1, a2, b, b1, b2, gam, gam1, gam2;
  GEN ldata2;
  long d1, t;
  long prec = nbits2prec(bitprec);
  ldata1 = ldata_newprec(lfunmisc_to_ldata_shallow(ldata1), prec);
  ldata2 = lfunmisc_to_ldata_shallow(chi);
  t = ldata_get_type(ldata2);
  a1 = ldata_get_an(ldata1);
  a2 = ldata_get_an(ldata2);
  if (t == t_LFUN_ZETA)
    return gerepilecopy(ltop, ldata1);
  if (t != t_LFUN_CHIZ && t != t_LFUN_KRONECKER &&
    ( t != t_LFUN_CHIGEN || nf_get_degree(bnr_get_nf(gmael(a2,2,1))) != 1))
    pari_err_TYPE("lfuntwist", chi);
  N1 = ldata_get_conductor(ldata1);
  N2 = ldata_get_conductor(ldata2);
  if (!gequal1(gcdii(N1, N2)))
    pari_err_IMPL("lfuntwist (conductors not coprime)");
  k = ldata_get_k(ldata1);
  d1 = ldata_get_degree(ldata1);
  N = gmul(N1, gpowgs(N2, d1));
  gam1 = ldata_get_gammavec(ldata1);
  gam2 = ldata_get_gammavec(ldata2);
  if (gequal0(gel(gam2, 1)))
    gam = gam1;
  else
    gam = gamma_imagchi(ldata_get_gammavec(ldata1), gaddgs(k,-1));
  if (!gam) pari_err_IMPL("lfuntwist (gammafactors)");
  b1 = ldata_get_dual(ldata1);
  b2 = ldata_get_dual(ldata2);
  a = tag(mkvec2(a1, a2), t_LFUN_TWIST);
  if (typ(b1)==t_INT)
    b = signe(b1) && signe(b2) ? gen_0: gen_1;
  else
    b = tag(mkvec2(b1,lfunconj(a2)), t_LFUN_TWIST);
  L = mkvecn(6, a, b, gam, k, N, gen_0);
  return gerepilecopy(ltop, L);
}

static GEN
lfundualpoles(GEN ldata, GEN reno)
{
  long l, j;
  GEN k = ldata_get_k(ldata);
  GEN r = gel(reno,2), eno = gel(reno,3), R;
  R = cgetg_copy(r, &l);
  for (j = 1; j < l; j++)
  {
    GEN b = gmael(r,j,1), e = gmael(r,j,2);
    long v = varn(e);
    GEN E = gsubst(gdiv(e, eno), v, gneg(pol_x(v)));
    gel(R,l-j) = mkvec2(gsub(k,b), E);
  }
  return R;
}

static GEN
ginvvec(GEN x)
{
  if (is_vec_t(typ(x)))
    pari_APPLY_same(ginv(gel(x,i)))
  else
    return ginv(x);
}

GEN
lfundual(GEN L, long bitprec)
{
  pari_sp av = avma;
  long prec = nbits2prec(bitprec);
  GEN ldata = ldata_newprec(lfunmisc_to_ldata_shallow(L), prec);
  GEN a = ldata_get_an(ldata), b = ldata_get_dual(ldata);
  GEN e = ldata_get_rootno(ldata);
  GEN ldual, ad, bd, ed, Rd = NULL;
  if (typ(b) == t_INT)
  {
    ad = equali1(b) ? lfunconj(a): a;
    bd = b;
  }
  else { ad = b; bd = a; }
  if (lg(ldata)==8)
  {
    GEN reno = lfunrootres(ldata, bitprec);
    e = gel(reno,3);
    Rd = lfundualpoles(ldata, reno);
  }
  ed = isintzero(e) ? e: ginvvec(e);
  ldual = mkvecn(Rd ? 7:6, ad, bd, gel(ldata,3), gel(ldata,4), gel(ldata,5), ed, Rd);
  return gerepilecopy(av, ldual);
}

static GEN
RgV_Rg_translate(GEN x, GEN s)
{ pari_APPLY_same(gadd(gel(x,i),s)) }

static GEN
pole_translate(GEN x, GEN s, GEN Ns)
{
  x = shallowcopy(x);
  gel(x,1) = gadd(gel(x,1), s);
  if (Ns)
    gel(x,2) = gmul(gel(x,2), Ns);
  return x;
}

static GEN
poles_translate(GEN x, GEN s, GEN Ns)
{ pari_APPLY_same(pole_translate(gel(x,i), s, Ns)) }

/* r / x + O(1) */
static GEN
simple_pole(GEN r)
{
  GEN S;
  if (isintzero(r)) return gen_0;
  S = deg1ser_shallow(gen_0, r, 1);
  setvalser(S, -1); return S;
}

GEN
lfunshift(GEN ldata, GEN s, long flag, long bitprec)
{
  pari_sp ltop = avma;
  GEN k, k1, L, N, a, b, gam, eps, res;
  long prec = nbits2prec(bitprec);
  if (!is_rational_t(typ(s))) pari_err_TYPE("lfunshift",s);
  ldata = ldata_newprec(lfunmisc_to_ldata_shallow(ldata), prec);
  a = ldata_get_an(ldata);
  b = ldata_get_dual(ldata);
  gam = RgV_Rg_translate(ldata_get_gammavec(ldata), gneg(s));
  k = gadd(ldata_get_k(ldata), gmul2n(s, 1));
  k1 = gadd(ldata_get_k1(ldata), s);
  N = ldata_get_conductor(ldata);
  eps = ldata_get_rootno(ldata);
  res = ldata_get_residue(ldata);
  a = tag(mkvec2(a, s), t_LFUN_SHIFT);
  if (typ(b) != t_INT)
    b = tag(mkvec2(b, s), t_LFUN_SHIFT);
  if (res)
    switch(typ(res))
    {
    case t_VEC:
      res = poles_translate(res, s, NULL);
      break;
    case t_COL:
      res = poles_translate(res, s, gpow(N, gmul2n(s, -1), prec));
      break;
    default:
      res = mkvec(mkvec2(gsub(k, s), simple_pole(res)));
    }
  L = mkvecn(res ? 7: 6, a, b, gam, mkvec2(k, k1), N, eps, res);
  if (flag) L = lfunmul_k(ldata, L, gsub(k, s), bitprec);
  return gerepilecopy(ltop, L);
}

/*****************************************************************/
/*  L-series from closure                                        */
/*****************************************************************/
static GEN
localfactor(void *E, GEN p, long n)
{
  GEN s = closure_callgen2((GEN)E, p, utoi(n));
  return direuler_factor(s, n);
}
static GEN
vecan_closure(GEN a, long L, long prec)
{
  long ta = typ(a);
  GEN gL, Sbad = NULL;

  if (!L) return cgetg(1,t_VEC);
  if (ta == t_VEC)
  {
    long l = lg(a);
    if (l == 1) pari_err_TYPE("vecan_closure", a);
    ta = typ(gel(a,1));
    /* regular vector, return it */
    if (ta != t_CLOSURE) return vecslice(a, 1, minss(L,l-1));
    if (l != 3) pari_err_TYPE("vecan_closure", a);
    Sbad = gel(a,2);
    if (typ(Sbad) != t_VEC) pari_err_TYPE("vecan_closure", a);
    a = gel(a,1);
  }
  else if (ta != t_CLOSURE) pari_err_TYPE("vecan_closure", a);
  push_localprec(prec);
  gL = stoi(L);
  switch(closure_arity(a))
  {
    case 2:
      a = direuler_bad((void*)a, localfactor, gen_2, gL,gL, Sbad);
      break;
    case 1:
      a = closure_callgen1(a, gL);
      if (typ(a) != t_VEC) pari_err_TYPE("vecan_closure", a);
      break;
    default: pari_err_TYPE("vecan_closure [wrong arity]", a);
      a = NULL; /*LCOV_EXCL_LINE*/
  }
  pop_localprec(); return a;
}

static GEN
eulerf_closure(GEN a, GEN p, long prec)
{
  long ta = typ(a);
  GEN Sbad = NULL, f;

  if (ta == t_VEC)
  {
    long l = lg(a);
    if (l == 1) pari_err_TYPE("vecan_closure", a);
    ta = typ(gel(a,1));
    /* regular vector, return it */
    if (ta != t_CLOSURE) return NULL;
    if (l != 3) pari_err_TYPE("vecan_closure", a);
    Sbad = gel(a,2);
    if (typ(Sbad) != t_VEC) pari_err_TYPE("vecan_closure", a);
    a = gel(a,1);
  }
  else if (ta != t_CLOSURE) pari_err_TYPE("vecan_closure", a);
  push_localprec(prec);
  switch(closure_arity(a))
  {
    case 2:
      f = closure_callgen2(a, p, mkoo()); break;
    case 1:
      f = NULL; break;
    default:
      f = NULL; pari_err_TYPE("vecan_closure", a);
  }
  pop_localprec(); return f;
}

/*****************************************************************/
/*  L-series of Dirichlet characters.                            */
/*****************************************************************/

static GEN
lfunzeta(void)
{
  GEN zet = mkvecn(7, NULL, gen_0, NULL, gen_1, gen_1, gen_1, gen_1);
  gel(zet,1) = tag(gen_1, t_LFUN_ZETA);
  gel(zet,3) = mkvec(gen_0);
  return zet;
}

static GEN
vecan_Kronecker(GEN D, long n)
{
  GEN v = cgetg(n+1, t_VECSMALL);
  ulong Du = itou_or_0(D);
  long i, id, d = Du ? minuu(Du, n): n;
  for (i = 1; i <= d; i++) v[i] = krois(D,i);
  for (id = i; i <= n; i++,id++) /* periodic mod d */
  {
    if (id > d) id = 1;
    gel(v, i) = gel(v, id);
  }
  return v;
}

static GEN
lfunchiquad(GEN D)
{
  GEN r;
  D = coredisc(D);
  if (equali1(D)) return lfunzeta();
  if (!isfundamental(D)) pari_err_TYPE("lfunchiquad [not primitive]", D);
  r = mkvecn(6, NULL, gen_0, NULL, gen_1, NULL, gen_1);
  gel(r,1) = tag(icopy(D), t_LFUN_KRONECKER);
  gel(r,3) = mkvec(signe(D) < 0? gen_1: gen_0);
  gel(r,5) = mpabs(D);
  return r;
}

/* Begin Hecke characters. Here a character is assumed to be given by a
   vector on the generators of the ray class group clgp of CL_m(K).
   If clgp = [h,[d1,...,dk],[g1,...,gk]] with dk|...|d2|d1, a character chi
   is given by [a1,a2,...,ak] such that chi(gi)=\zeta_di^ai. */

/* Value of CHI on x, coprime to bnr.mod */
static GEN
chigeneval_i(GEN logx, GEN d, GEN nchi, GEN z, long prec)
{
  pari_sp av = avma;
  GEN e = FpV_dotproduct(nchi, logx, d);
  if (!is_vec_t(typ(z)))
    return gerepileupto(av, gpow(z, e, prec));
  else
  {
    ulong i = itou(e);
    set_avma(av); return gel(z, i+1);
  }
}

static GEN
chigenevalvec(GEN logx, GEN nchi, GEN z, long prec, long multi)
{
  GEN d = gel(nchi,1), x = gel(nchi, 2);
  if (multi)
    pari_APPLY_same(chigeneval_i(logx, d, gel(x,i), z, prec))
  else
    return chigeneval_i(logx, d, x, z, prec);
}

/* return x + yz; y != 0; z = 0,1 "often"; x = 0 "often" */
static GEN
gaddmul(GEN x, GEN y, GEN z)
{
  pari_sp av;
  if (typ(z) == t_INT)
  {
    if (!signe(z)) return x;
    if (equali1(z)) return gadd(x,y);
  }
  if (isintzero(x)) return gmul(y,z);
  av = avma;
  return gerepileupto(av, gadd(x, gmul(y,z)));
}

static GEN
gaddmulvec(GEN x, GEN y, GEN z, long multi)
{
  if (multi)
    pari_APPLY_same(gaddmul(gel(x,i),gel(y,i),gel(z,i)))
  else
    return gaddmul(x,y,z);
}

static GEN
mkvchi(GEN chi, long n)
{
  GEN v;
  if (lg(chi) > 1 && is_vec_t(typ(gel(chi,1))))
  {
    long d = lg(chi)-1;
    v = const_vec(n, zerovec(d));
    gel(v,1) = const_vec(d, gen_1);
  }
  else
    v = vec_ei(n, 1);
  return v;
}

static GEN
vecan_chiZ(GEN an, long n, long prec)
{
  forprime_t iter;
  GEN G = gel(an,1);
  GEN nchi = gel(an,2), gord = gel(nchi,1), chi = gel(nchi,2), z;
  GEN gp = cgetipos(3), v = mkvchi(chi, n);
  GEN N = znstar_get_N(G);
  long ord = itos_or_0(gord);
  ulong Nu = itou_or_0(N);
  long i, id, d = Nu ? minuu(Nu, n): n;
  long multichi= (lg(chi) > 1 && is_vec_t(typ(gel(chi,1))));
  ulong p;
  if (!multichi && ord && n > (ord>>4))
  {
    GEN w = ncharvecexpo(G, nchi);
    z = grootsof1(ord, prec);
    for (i = 1; i <= d; i++)
      if (w[i] >= 0) gel(v, i) = gel(z, w[i]+1);
  }
  else
  {
    z = rootsof1_cx(gord, prec);
    u_forprime_init(&iter, 2, d);
    while ((p = u_forprime_next(&iter)))
    {
      GEN ch;
      ulong k;
      if (!umodiu(N,p)) continue;
      gp[2] = p;
      ch = chigenevalvec(znconreylog(G, gp), nchi, z, prec, multichi);
      gel(v, p)  = ch;
      for (k = 2*p; k <= (ulong)d; k += p)
        gel(v, k) = gaddmulvec(gel(v, k), ch, gel(v, k/p), multichi);
    }
  }
  for (id = i = d+1; i <= n; i++,id++) /* periodic mod d */
  {
    if (id > d) id = 1;
    gel(v, i) = gel(v, id);
  }
  return v;
}

static GEN
eulerf_chiZ(GEN an, GEN p, long prec)
{
  GEN G = gel(an,1);
  GEN nchi = gel(an,2), gord = gel(nchi,1), chi = gel(nchi,2);
  long multichi= (lg(chi) > 1 && is_vec_t(typ(gel(chi,1))));
  GEN z = rootsof1_cx(gord, prec);
  GEN N = znstar_get_N(G);
  GEN ch = dvdii(N,p) ? gen_0: chigenevalvec(znconreylog(G, p), nchi, z, prec, multichi);
  return mkrfrac(gen_1, deg1pol_shallow(gneg(ch), gen_1,0));
}

static GEN
vecan_chigen(GEN an, long n, long prec)
{
  forprime_t iter;
  GEN bnr = gel(an,1), nf = bnr_get_nf(bnr);
  GEN nchi = gel(an,2), gord = gel(nchi,1), chi = gel(nchi,2), z;
  GEN gp = cgetipos(3), v = mkvchi(chi, n);
  GEN N = gel(bnr_get_mod(bnr), 1), NZ = gcoeff(N,1,1);
  long ord = itos_or_0(gord);
  long multichi= (lg(chi) > 1 && is_vec_t(typ(gel(chi,1))));
  ulong p;

  if (ord && n > (ord>>4))
    z = grootsof1(ord, prec);
  else
    z = rootsof1_cx(gord, prec);

  if (nf_get_degree(nf) == 1)
  {
    ulong Nu = itou_or_0(NZ);
    long i, id, d = Nu ? minuu(Nu, n): n;
    u_forprime_init(&iter, 2, d);
    while ((p = u_forprime_next(&iter)))
    {
      GEN ch;
      ulong k;
      if (!umodiu(NZ,p)) continue;
      gp[2] = p;
      ch = chigenevalvec(isprincipalray(bnr,gp), nchi, z, prec, multichi);
      gel(v, p)  = ch;
      for (k = 2*p; k <= (ulong)d; k += p)
        gel(v, k) = gaddmulvec(gel(v, k), ch, gel(v, k/p), multichi);
    }
    for (id = i = d+1; i <= n; i++,id++) /* periodic mod d */
    {
      if (id > d) id = 1;
      gel(v, i) = gel(v, id);
    }
  }
  else
  {
    GEN BOUND = stoi(n);
    u_forprime_init(&iter, 2, n);
    while ((p = u_forprime_next(&iter)))
    {
      GEN L;
      long j;
      int check = !umodiu(NZ,p);
      gp[2] = p;
      L = idealprimedec_limit_norm(nf, gp, BOUND);
      for (j = 1; j < lg(L); j++)
      {
        GEN pr = gel(L, j), ch;
        ulong k, q;
        if (check && idealval(nf, N, pr)) continue;
        ch = chigenevalvec(isprincipalray(bnr,pr), nchi, z, prec, multichi);
        q = upr_norm(pr);
        gel(v, q) = gadd(gel(v, q), ch);
        for (k = 2*q; k <= (ulong)n; k += q)
          gel(v, k) = gaddmulvec(gel(v, k), ch, gel(v, k/q), multichi);
      }
    }
  }
  return v;
}

static GEN
eulerf_chigen(GEN an, GEN p, long prec)
{
  GEN bnr = gel(an,1), nf = bnr_get_nf(bnr);
  GEN nchi = gel(an,2), gord = gel(nchi,1), chi = gel(nchi,2), z;
  GEN N = gel(bnr_get_mod(bnr), 1), NZ = gcoeff(N,1,1), f;
  long multichi= (lg(chi) > 1 && is_vec_t(typ(gel(chi,1))));

  z = rootsof1_cx(gord, prec);
  if (nf_get_degree(nf) == 1)
  {
    GEN ch;
    if (dvdii(NZ,p)) ch = gen_0;
    else
    {
      ch = chigenevalvec(isprincipalray(bnr,p), nchi, z, prec, multichi);
      if (typ(ch)==t_VEC) return NULL;
    }
    f = deg1pol_shallow(gneg(ch), gen_1, 0);
  }
  else
  {
    int check = dvdii(NZ,p);
    GEN L = idealprimedec(nf, p);
    long j, lL = lg(L);
    f = pol_1(0);
    for (j = 1; j < lL; j++)
    {
      GEN pr = gel(L, j), ch;
      if (check && idealval(nf, N, pr)) ch = gen_0;
      else
      ch = chigenevalvec(isprincipalray(bnr,pr), nchi, z, prec, multichi);
      if (typ(ch)==t_VEC) return NULL;
      f = gmul(f, gsub(gen_1, monomial(ch, pr_get_f(pr), 0)));
    }
  }
  return mkrfrac(gen_1,f);
}

static GEN
vec01(long r1, long r2)
{
  long d = r1+r2, i;
  GEN v = cgetg(d+1,t_VEC);
  for (i = 1; i <= r1; i++) gel(v,i) = gen_0;
  for (     ; i <= d;  i++) gel(v,i) = gen_1;
  return v;
}

/* true nf or t_POL */
static GEN
lfunzetak_i(GEN T)
{
  GEN Vga, N;
  long r1, r2;
  if (typ(T) == t_POL)
  {
    T = nfinit0(T, nf_NOLLL, DEFAULTPREC);
    if (lg(T) == 3) T = gel(T,1); /* [nf,change of var] */
  }
  nf_get_sign(T,&r1,&r2); Vga = vec01(r1+r2,r2);
  N = absi_shallow(nf_get_disc(T));
  return mkvecn(7, tag(T,t_LFUN_NF), gen_0, Vga, gen_1, N, gen_1, gen_0);
}
/* truen nf or t_POL */
static GEN
lfunzetak(GEN T)
{ pari_sp av = avma; return gerepilecopy(av, lfunzetak_i(T)); }

/* v = vector of normalized characters of order dividing o; renormalize
 * so that all have same apparent order o */
static GEN
char_renormalize(GEN v, GEN o)
{
  long i, l;
  GEN w = cgetg_copy(v, &l);
  for (i = 1; i < l; i++)
  {
    GEN C = gel(v,i), oc = gel(C,1), c = gel(C,2);
    if (!equalii(o, oc)) c = gmul(c, diviiexact(o, oc));
    gel(w,i) = c;
  }
  return w;
}
/* G is a bid of nftyp typ_BIDZ */
static GEN
lfunchiZ(GEN G, GEN CHI)
{
  pari_sp av = avma;
  GEN sig = NULL, N = bid_get_ideal(G), nchi, r;
  int real;
  long s;

  if (typ(N) != t_INT) pari_err_TYPE("lfunchiZ", G);
  if (typ(CHI) == t_VEC && !RgV_is_ZV(CHI))
  {
    GEN C, G0 = G, o = gen_1;
    long i, l = lg(CHI);
    nchi = cgetg(l, t_VEC);
    N = znconreyconductor(G, gel(CHI,1), &C);
    if (typ(N) != t_INT) G = znstar0(N, 1);
    s = zncharisodd(G, C);
    for (i = 1; i < l; i++)
    {
      if (i > 1)
      {
        if (!gequal(N, znconreyconductor(G0, gel(CHI,i), &C))
            || zncharisodd(G, C) != s)
          pari_err_TYPE("lfuncreate [different conductors]", CHI);
      }
      C = znconreylog_normalize(G, C);
      o = lcmii(o, gel(C,1)); /* lcm with charorder */
      gel(nchi,i) = C;
    }
    nchi = mkvec2(o, char_renormalize(nchi, o));
    if (typ(N) != t_INT) N = gel(N,1);
  }
  else
  {
    N = znconreyconductor(G, CHI, &CHI);
    if (typ(N) != t_INT)
    {
      if (equali1(gel(N,1))) { set_avma(av); return lfunzeta(); }
      G = znstar0(N, 1);
      N = gel(N,1);
    }
    /* CHI now primitive on G */
    switch(itou_or_0(zncharorder(G, CHI)))
    {
      case 1: set_avma(av); return lfunzeta();
      case 2: if (zncharisodd(G,CHI)) N = negi(N);
              return gerepileupto(av, lfunchiquad(N));
    }
    nchi = znconreylog_normalize(G, CHI);
    s = zncharisodd(G, CHI);
  }
  sig = mkvec(s? gen_1: gen_0);
  real = abscmpiu(gel(nchi,1), 2) <= 0;
  r = mkvecn(6, tag(mkvec2(G,nchi), t_LFUN_CHIZ),
                real? gen_0: gen_1, sig, gen_1, N, gen_0);
  return gerepilecopy(av, r);
}

static GEN
lfunchigen(GEN bnr, GEN CHI)
{
  pari_sp av = avma;
  GEN N, sig, Ldchi, nf, nchi, NN;
  long r1, r2, n1;
  int real;

  if (typ(CHI) == t_VEC && !RgV_is_ZV(CHI))
  {
    long map, i, l = lg(CHI);
    GEN bnr0 = bnr, D, chi = gel(CHI,1), o = gen_1;
    nchi = cgetg(l, t_VEC);
    bnr_char_sanitize(&bnr, &chi);
    D = cyc_normalize(bnr_get_cyc(bnr));
    N = bnr_get_mod(bnr);
    map = (bnr != bnr0);
    for (i = 1; i < l; i++)
    {
      if (i > 1)
      {
        chi = gel(CHI,i);
        if (!map)
        {
          if (!bnrisconductor(bnr, chi))
            pari_err_TYPE("lfuncreate [different conductors]", CHI);
        }
        else
        {
          if (!gequal(bnrconductor_raw(bnr0, chi), N))
            pari_err_TYPE("lfuncreate [different conductors]", CHI);
          chi = bnrchar_primitive_raw(bnr0, bnr, chi);
        }
      }
      chi = char_normalize(chi, D);
      o = lcmii(o, gel(chi,1)); /* lcm with charorder */
      gel(nchi,i) = chi;
    }
    nchi = mkvec2(o, char_renormalize(nchi, o));
  }
  else
  {
    bnr_char_sanitize(&bnr, &CHI);
    nchi = NULL; /* now CHI is primitive wrt bnr */
  }

  N = bnr_get_mod(bnr);
  nf = bnr_get_nf(bnr);
  n1 = lg(vec01_to_indices(gel(N,2))) - 1; /* vecsum(N[2]) */
  N = gel(N,1);
  NN = mulii(idealnorm(nf, N), absi_shallow(nf_get_disc(nf)));
  if (!nchi)
  {
    if (equali1(NN)) { set_avma(av); return lfunzeta(); }
    if (ZV_equal0(CHI)) return gerepilecopy(av, lfunzetak_i(bnr_get_nf(bnr)));
    nchi = char_normalize(CHI, cyc_normalize(bnr_get_cyc(bnr)));
  }
  real = abscmpiu(gel(nchi,1), 2) <= 0;
  nf_get_sign(nf, &r1, &r2);
  sig = vec01(r1+r2-n1, r2+n1);
  Ldchi = mkvecn(6, tag(mkvec2(bnr, nchi), t_LFUN_CHIGEN),
                    real? gen_0: gen_1, sig, gen_1, NN, gen_0);
  return gerepilecopy(av, Ldchi);
}

/* Find all characters of clgp whose kernel contain group given by HNF H.
 * Set *pcnj[i] to the conductor */
static GEN
chigenkerfind(GEN bnr, GEN H, GEN *pcnj)
{
  GEN res, cnj, L = bnrchar(bnr, H, NULL);
  long i, k, l = lg(L);

  res = cgetg(l, t_VEC);
  *pcnj = cnj = cgetg(l, t_VEC);
  for (i = k = 1; i < l; i++)
  {
    GEN chi = gel(L,i);
    gel(res, k) = chi;
    gel(cnj, k) = ZV_equal0(chi)? gen_0: bnrconductorofchar(bnr, chi);
    k++;
  }
  setlg(cnj, k);
  setlg(res, k); return res;
}

static GEN
vec_classes(GEN A, GEN F)
{
  GEN w = vec_equiv(F);
  long i, l = lg(w);
  GEN V = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(V,i) = vecpermute(A,gel(w,i));
  return V;
}

static GEN
abelrel_pfactor(GEN bnr, GEN pr, GEN U, GEN D, GEN h)
{
  GEN v = bnrisprincipalmod(bnr, pr, h, 0);
  GEN E = ZV_ZV_mod(ZM_ZC_mul(U, v), D);
  ulong o = itou(charorder(D, E)), f = pr_get_f(pr);
  return gpowgs(gsub(gen_1, monomial(gen_1, f * o, 0)), itou(h) / o);
}

static GEN
abelrel_factor(GEN bnr, GEN C, GEN p, GEN mod, GEN U, GEN D, GEN h)
{
  GEN nf = bnr_get_nf(bnr), F = pol_1(0), prid = idealprimedec(nf,p);
  GEN mod2 = shallowcopy(mod);
  long i, l = lg(prid);
  for (i = 1; i < l; i++)
  {
    GEN pr = gel(prid, i), Fpr;
    long v = idealval(nf,mod,pr);
    if (v > 0)
    {
      GEN bnr2, C2, U2, D2, h2;
      gel(mod2, 1) = idealdivpowprime(nf, gel(mod, 1), pr, utoi(v));
      bnr2 = bnrinitmod(bnr, mod2, 0, h);
      C2 = bnrmap(bnrmap(bnr, bnr2), C);
      D2 = ZM_snfall_i(C2, &U2, NULL, 1);
      h2 = ZV_prod(D2);
      Fpr = abelrel_pfactor(bnr2, pr, U2, D2, h2);
    }
    else
      Fpr = abelrel_pfactor(bnr, pr, U, D, h);
    F = ZX_mul(F, Fpr);
  }
  return gcopy(mkrfrac(gen_1, F));
}

static GEN
eulerf_abelrel(GEN an, GEN p)
{
  GEN bnr = gel(an,1), C = gel(an,2), mod = gel(an,3);
  GEN U, D = ZM_snfall_i(C, &U, NULL, 1), h = ZV_prod(D);
  return abelrel_factor(bnr, C, p, mod, U, D, h);
}

struct direuler_abelrel
{
  GEN bnr, C, mod, U, D, h;
};

static GEN
_direuler_abelrel(void *E, GEN p)
{
  struct direuler_abelrel *s = (struct direuler_abelrel*) E;
  return abelrel_factor(s->bnr, s->C, p, s->mod, s->U, s->D, s->h);
}

static GEN
vecan_abelrel(GEN an, long N)
{
  struct direuler_abelrel s;
  s.bnr = gel(an,1);
  s.C   = gel(an,2);
  s.mod = gel(an,3);
  s.D = ZM_snfall_i(s.C, &s.U, NULL, 1);
  s.h = ZV_prod(s.D);
  return direuler((void*)&s, _direuler_abelrel, gen_1, stoi(N), NULL);
}

static GEN
lfunabelrel_i(GEN bnr, GEN H, GEN mod)
{
  GEN NrD = bnrdisc(bnr, H, 0), N = absi_shallow(gel(NrD,3));
  long n = itos(gel(NrD,1)), r1 = itos(gel(NrD,2)), r2 = (n-r1)>>1;
  if (!mod) mod = bnrconductor(bnr, H, 0);
  return mkvecn(7, tag(mkvec3(bnr,H,mod),t_LFUN_ABELREL),
                   gen_0, vec01(r1+r2, r2), gen_1, N, gen_1, gen_0);
}
static GEN
lfunabelrel(GEN bnr, GEN H, GEN mod)
{ pari_sp av = avma; return gerepilecopy(av, lfunabelrel_i(bnr, H, mod)); }

GEN
lfunabelianrelinit(GEN bnr, GEN H, GEN dom, long der, long bitprec)
{
  GEN X, cnj, M, D,C ;
  long l, i;
  C = chigenkerfind(bnr, H, &cnj);
  C = vec_classes (C, cnj);
  X = cgetg_copy(C,&l);
  for (i = 1; i < l; ++i)
  {
    GEN chi = gel(C,i);
    GEN L = lfunchigen(bnr, lg(chi)==2 ? gel(chi,1): chi);
    gel(X,i) = lfuninit(L, dom, der, bitprec);
  }
  M = mkvec3(X, const_vecsmall(l-1, 1), const_vecsmall(l-1, 0));
  D = mkvec2(dom, mkvecsmall2(der, bitprec));
  return lfuninit_make(t_LDESC_PRODUCT, lfunabelrel_i(bnr, H, NULL), M, D);
}

/*****************************************************************/
/*                 Dedekind zeta functions                       */
/*****************************************************************/
/* true nf */
static GEN
dirzetak0(GEN nf, ulong N)
{
  GEN vect, c, c2, T = nf_get_pol(nf), index = nf_get_index(nf);
  pari_sp av = avma, av2;
  const ulong SQRTN = usqrt(N);
  ulong i, p, lx;
  long court[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3),0};
  forprime_t S;

  c  = cgetalloc(N+1, t_VECSMALL);
  c2 = cgetalloc(N+1, t_VECSMALL);
  c2[1] = c[1] = 1; for (i=2; i<=N; i++) c[i] = 0;
  u_forprime_init(&S, 2, N); av2 = avma;
  while ( (p = u_forprime_next(&S)) )
  {
    set_avma(av2);
    if (umodiu(index, p)) /* p does not divide index */
      vect = gel(Flx_degfact(ZX_to_Flx(T,p), p),1);
    else
    {
      court[2] = p;
      vect = idealprimedec_degrees(nf,court);
    }
    lx = lg(vect);
    if (p <= SQRTN)
      for (i=1; i<lx; i++)
      {
        ulong qn, q = upowuu(p, vect[i]); /* Norm P[i] */
        if (!q || q > N) break;
        memcpy(c2 + 2, c + 2, (N-1)*sizeof(long));
        /* c2[i] <- c[i] + sum_{k = 1}^{v_q(i)} c[i/q^k] for all i <= N */
        for (qn = q; qn <= N; qn *= q)
        {
          ulong k0 = N/qn, k, k2; /* k2 = k*qn */
          for (k = k0, k2 = k*qn; k > 0; k--, k2 -=qn) c2[k2] += c[k];
          if (q > k0) break; /* <=> q*qn > N */
        }
        swap(c, c2);
      }
    else /* p > sqrt(N): simpler */
      for (i=1; i<lx; i++)
      {
        ulong k, k2; /* k2 = k*p */
        if (vect[i] > 1) break;
        /* c2[i] <- c[i] + sum_{k = 1}^{v_q(i)} c[i/q^k] for all i <= N */
        for (k = N/p, k2 = k*p; k > 0; k--, k2 -= p) c[k2] += c[k];
      }
  }
  set_avma(av);
  pari_free(c2); return c;
}

static GEN
eulerf_zetak(GEN nf, GEN p)
{
  GEN v, f = pol_1(0);
  long i, l;
  if (dvdii(nf_get_index(nf), p)) /* p does not divide index */
    v = idealprimedec_degrees(nf,p);
  else
    v = gel(FpX_degfact(nf_get_pol(nf), p), 1);
  l = lg(v);
  for (i = 1; i < l; i++) f = ZX_sub(f, RgX_shift_shallow(f, v[i]));
  retmkrfrac(gen_1, ZX_copy(f));
}

GEN
dirzetak(GEN nf, GEN b)
{
  GEN z, c;
  long n;

  if (typ(b) != t_INT) pari_err_TYPE("dirzetak",b);
  if (signe(b) <= 0) return cgetg(1,t_VEC);
  nf = checknf(nf);
  n = itou_or_0(b); if (!n) pari_err_OVERFLOW("dirzetak");
  c = dirzetak0(nf, n);
  z = vecsmall_to_vec(c); pari_free(c); return z;
}

static GEN
linit_get_mat(GEN linit)
{
  if (linit_get_type(linit)==t_LDESC_PRODUCT)
    return lfunprod_get_fact(linit_get_tech(linit));
  else
    return mkvec3(mkvec(linit), mkvecsmall(1), mkvecsmall(0));
}

static GEN
lfunproduct(GEN ldata, GEN linit1, GEN linit2, GEN domain)
{
  GEN M1 = linit_get_mat(linit1);
  GEN M2 = linit_get_mat(linit2);
  GEN M3 = mkvec3(shallowconcat(gel(M1, 1), gel(M2, 1)),
                  vecsmall_concat(gel(M1, 2), gel(M2, 2)),
                  vecsmall_concat(gel(M1, 3), gel(M2, 3)));
  return lfuninit_make(t_LDESC_PRODUCT, ldata, M3, domain);
}
static GEN lfunzetakinit_artin(GEN nf, GEN gal, GEN dom, long der, long bit);
/* true nf */
static GEN
lfunzetakinit_quotient(GEN nf, GEN polk, GEN dom, long der, long bitprec)
{
  GEN ak, an, nfk, Vga, ldata, N, Lk, LKk, domain;
  long r1k, r2k, r1, r2;

  nf_get_sign(nf,&r1,&r2);
  nfk = nfinit(polk, nbits2prec(bitprec));
  Lk = lfunzetakinit(nfk, dom, der, bitprec); /* zeta_k */
  nf_get_sign(nfk,&r1k,&r2k);
  Vga = vec01((r1+r2) - (r1k+r2k), r2-r2k);
  N = absi_shallow(diviiexact(nf_get_disc(nf), nf_get_disc(nfk)));
  ak = nf_get_degree(nf)==1 ? tag(gen_1, t_LFUN_ZETA): tag(nfk, t_LFUN_NF);
  an = tag(mkvec2(tag(nf,t_LFUN_NF), ak), t_LFUN_DIV);
  ldata = mkvecn(6, an, gen_0, Vga, gen_1, N, gen_1);
  LKk = lfuninit(ldata, dom, der, bitprec); /* zeta_K/zeta_k */
  domain = mkvec2(dom, mkvecsmall2(der, bitprec));
  return lfunproduct(lfunzetak_i(nf), Lk, LKk, domain);
}
/* true nf */
GEN
lfunzetakinit(GEN nf, GEN dom, long der, long bitprec)
{
  long n, d = nf_get_degree(nf);
  GEN L, Q, R, G, T = nf_get_pol(nf);
  if (d == 1) return lfuninit(lfunzeta(), dom, der, bitprec);
  G = galoisinit(nf, NULL);
  if (isintzero(G))
  {
    GEN S = nfsubfields(nf, 0); n = lg(S)-1;
    return lfunzetakinit_quotient(nf, gmael(S,n-1,1), dom, der, bitprec);
  }
  if (!group_isabelian(galois_group(G)))
    return lfunzetakinit_artin(nf, G, dom, der, bitprec);
  Q = Buchall(pol_x(1), 0, nbits2prec(bitprec));
  T = shallowcopy(T); setvarn(T,0);
  R = rnfconductor0(Q, T, 1);
  L = lfunabelianrelinit(gel(R,2), gel(R,3), dom, der, bitprec);
  delete_var(); return L;
}

/***************************************************************/
/*             Elliptic Curves and Modular Forms               */
/***************************************************************/

static GEN
lfunellnf(GEN e)
{
  pari_sp av = avma;
  GEN ldata = cgetg(7, t_VEC), nf = ellnf_get_nf(e);
  GEN N = gel(ellglobalred(e), 1);
  long n = nf_get_degree(nf);
  gel(ldata, 1) = tag(e, t_LFUN_ELL);
  gel(ldata, 2) = gen_0;
  gel(ldata, 3) = vec01(n, n);
  gel(ldata, 4) = gen_2;
  gel(ldata, 5) = mulii(idealnorm(nf,N), sqri(nf_get_disc(nf)));
  gel(ldata, 6) = stoi(ellrootno_global(e));
  return gerepilecopy(av, ldata);
}

static GEN
lfunellQ(GEN e)
{
  pari_sp av = avma;
  GEN ldata = cgetg(7, t_VEC);
  gel(ldata, 1) = tag(ellanal_globalred(e, NULL), t_LFUN_ELL);
  gel(ldata, 2) = gen_0;
  gel(ldata, 3) = mkvec2(gen_0, gen_1);
  gel(ldata, 4) = gen_2;
  gel(ldata, 5) = ellQ_get_N(e);
  gel(ldata, 6) = stoi(ellrootno_global(e));
  return gerepilecopy(av, ldata); /* ellanal_globalred not gerepile-safe */
}

static GEN
lfunell(GEN e)
{
  long t = ell_get_type(e);
  switch(t)
  {
    case t_ELL_Q: return lfunellQ(e);
    case t_ELL_NF:return lfunellnf(e);
  }
  pari_err_TYPE("lfun",e);
  return NULL; /*LCOV_EXCL_LINE*/
}

static GEN
ellsympow_gamma(long m)
{
  GEN V = cgetg(m+2, t_VEC);
  long i = 1, j;
  if (!odd(m)) gel(V, i++) = stoi(-2*(m>>2));
  for (j = (m+1)>>1; j > 0; i+=2, j--)
  {
    gel(V,i)   = stoi(1-j);
    gel(V,i+1) = stoi(1-j+1);
  }
  return V;
}

static GEN
ellsympow_trace(GEN p, GEN t, long m)
{
  long k, n = m >> 1;
  GEN tp = gpowers0(sqri(t), n, odd(m)? t: NULL);
  GEN pp = gen_1, b = gen_1, r = gel(tp,n+1);
  for(k=1; k<=n; k++)
  {
    GEN s;
    pp = mulii(pp, p);
    b  = diviuexact(muliu(b, (m-(2*k-1))*(m-(2*k-2))), k*(m-(k-1)));
    s = mulii(mulii(b, gel(tp,1+n-k)), pp);
    r = odd(k) ? subii(r, s): addii(r, s);
  }
  return r;
}

static GEN
ellsympow_abelian(GEN p, GEN ap, long m, long o)
{
  pari_sp av = avma;
  long i, M, n = (m+1)>>1;
  GEN pk, tv, pn, pm, F, v;
  if (!odd(o))
  {
    if (odd(m)) return pol_1(0);
    M = m >> 1; o >>= 1;
  }
  else
    M = m * ((o+1) >> 1);
  pk = gpowers(p,n); pn = gel(pk,n+1);
  tv = cgetg(m+2,t_VEC);
  gel(tv, 1) = gen_2;
  gel(tv, 2) = ap;
  for (i = 3; i <= m+1; i++)
    gel(tv,i) = subii(mulii(ap,gel(tv,i-1)), mulii(p,gel(tv,i-2)));
  pm = odd(m)? mulii(gel(pk,n), pn): sqri(pn); /* cheap p^m */
  F = deg2pol_shallow(pm, gen_0, gen_1, 0);
  v = odd(m) ? pol_1(0): deg1pol_shallow(negi(pn), gen_1, 0);
  for (i = M % o; i < n; i += o) /* o | m-2*i */
  {
    gel(F,3) = negi(mulii(gel(tv,m-2*i+1), gel(pk,i+1)));
    v = ZX_mul(v, F);
  }
  return gerepilecopy(av, v);
}

static GEN
ellsympow(GEN E, ulong m, GEN p, long n)
{
  pari_sp av = avma;
  GEN ap = ellap(E, p);
  if (n <= 2)
  {
    GEN t = ellsympow_trace(p, ap, m);
    return deg1pol_shallow(t, gen_1, 0);
  }
  else
    return gerepileupto(av, RgXn_inv_i(ellsympow_abelian(p, ap, m, 1), n));
}

GEN
direllsympow_worker(GEN P, ulong X, GEN E, ulong m)
{
  pari_sp av = avma;
  long i, l = lg(P);
  GEN W = cgetg(l, t_VEC);
  for(i = 1; i < l; i++)
  {
    ulong p = uel(P,i);
    long d = ulogint(X, p) + 1; /* minimal d such that p^d > X */
    gel(W,i) = ellsympow(E, m, utoi(uel(P,i)), d);
  }
  return gerepilecopy(av, mkvec2(P,W));
}

static GEN
eulerf_bad(GEN bad, GEN p)
{
  long i, l = lg(bad);
  for (i = 1; i < l; i++)
    if (equalii(gmael(bad,i,1), p))
      return gmael(bad,i,2);
  return NULL;
}

static GEN
vecan_ellsympow(GEN an, long n)
{
  GEN nn = utoi(n), crvm = gel(an,1), bad = gel(an,2);
  GEN worker = snm_closure(is_entry("_direllsympow_worker"), crvm);
  return pardireuler(worker, gen_2, nn, nn, bad);
}

static GEN
eulerf_ellsympow(GEN an, GEN p)
{
  GEN crvm = gel(an,1), bad = gel(an,2), E = gel(crvm,1);
  GEN f = eulerf_bad(bad, p);
  if (f) return f;
  retmkrfrac(gen_1,ellsympow_abelian(p, ellap(E, p), itos(gel(crvm,2)), 1));
}

static long
ellsympow_betam(long o, long m)
{
  const long c3[]={3, -1, 1};
  const long c12[]={6, -2, 2, 0, 4, -4};
  const long c24[]={12, -2, -4, 6, 4, -10};
  if (!odd(o) && odd(m)) return 0;
  switch(o)
  {
    case 1:  return m+1;
    case 2:  return m+1;
    case 3:  case 6: return (m+c3[m%3])/3;
    case 4:  return m%4 == 0 ? (m+2)/2: m/2;
    case 8:  return m%4 == 0 ? (m+4)/4: (m-2)/4;
    case 12: return (m+c12[(m%12)/2])/6;
    case 24: return (m+c24[(m%12)/2])/12;
  }
  return 0;
}

static long
ellsympow_epsm(long o, long m) { return m + 1 - ellsympow_betam(o, m); }

static GEN
ellsympow_multred(GEN E, GEN p, long m, long vN, long *cnd, long *w)
{
  if (vN == 1 || !odd(m))
  {
    GEN s = (odd(m) && signe(ellap(E,p)) < 0)? gen_1: gen_m1;
    *cnd = m;
    *w = odd(m)? ellrootno(E, p): 1;
    return deg1pol_shallow(s, gen_1, 0);
  }
  else
  {
    *cnd = equaliu(p,2)? ((m+1)>>1) * vN: m+1;
    *w = (m & 3) == 1? ellrootno(E, p): 1;
    return pol_1(0);
  }
}

static GEN
ellsympow_nonabelian(GEN p, long m, long bet)
{
 GEN q = powiu(p, m >> 1), q2 = sqri(q), F;
 if (odd(m))
 {
   q2 = mulii(q2, p); /* p^m */
   return gpowgs(deg2pol_shallow(q2, gen_0, gen_1, 0), bet>>1);
 }
 togglesign_safe(&q2);
 F = gpowgs(deg2pol_shallow(q2, gen_0, gen_1, 0), bet>>1);
 if (!odd(bet)) return F;
 if (m%4 != 2) togglesign_safe(&q);
 return gmul(F, deg1pol_shallow(q, gen_1, 0));
}

static long
safe_Z_pvalrem(GEN n, GEN p, GEN *pr)
{ return signe(n)==0? -1: Z_pvalrem(n, p, pr); }

static GEN
c4c6_ap(GEN c4, GEN c6, GEN p)
{
  GEN N = Fp_ellcard(Fp_muls(c4, -27, p), Fp_muls(c6, -54, p), p);
  return subii(addiu(p, 1), N);
}

static GEN
ellsympow_abelian_twist(GEN E, GEN p, long m, long o)
{
  GEN ap, c4t, c6t, c4 = ell_get_c4(E), c6 = ell_get_c6(E);
  long v4 = safe_Z_pvalrem(c4, p, &c4t);
  long v6 = safe_Z_pvalrem(c6, p, &c6t);
  if (v6>=0 && (v4==-1 || 3*v4>=2*v6)) c6 = c6t;
  if (v4>=0 && (v6==-1 || 3*v4<=2*v6)) c4 = c4t;
  ap = c4c6_ap(c4, c6, p);
  return ellsympow_abelian(p, ap, m, o);
}

static GEN
ellsympow_goodred(GEN E, GEN p, long m, long *cnd, long *w)
{
  long o = 12/cgcd(12, Z_pval(ell_get_disc(E), p));
  long bet = ellsympow_betam(o, m);
  long eps = m + 1 - bet;
  *w = odd(m) && odd(eps>>1) ? ellrootno(E,p): 1;
  *cnd = eps;
  if (umodiu(p, o) == 1)
    return ellsympow_abelian_twist(E, p, m, o);
  else
    return ellsympow_nonabelian(p, m, bet);
}

static long
ellsympow_inertia3(GEN E, long vN)
{
  long vD = Z_lval(ell_get_disc(E), 3);
  if (vN==2) return vD%2==0 ? 2: 4;
  if (vN==4) return vD%4==0 ? 3: 6;
  if (vN==3 || vN==5) return 12;
  return 0;
}

static long
ellsympow_deltam3(long o, long m, long vN)
{
  if (o==3 || o==6) return ellsympow_epsm(3, m);
  if (o==12 && vN ==3) return (ellsympow_epsm(3, m))/2;
  if (o==12 && vN ==5) return (ellsympow_epsm(3, m))*3/2;
  return 0;
}

static long
ellsympow_isabelian3(GEN E)
{
  ulong c4 = umodiu(ell_get_c4(E),81), c6 = umodiu(ell_get_c6(E), 243);
  return (c4 == 27 || (c4%27==9 && (c6==108 || c6==135)));
}

static long
ellsympow_rootno3(GEN E, GEN p, long o, long m)
{
  const long  w6p[]={1,-1,-1,-1,1,1};
  const long  w6n[]={-1,1,-1,1,-1,1};
  const long w12p[]={1,1,-1,1,1,1};
  const long w12n[]={-1,-1,-1,-1,-1,1};
  long w = ellrootno(E, p), mm = (m%12)>>1;
  switch(o)
  {
    case 2: return m%4== 1 ? -1: 1;
    case 6:  return w == 1 ? w6p[mm]: w6n[mm];
    case 12: return w == 1 ? w12p[mm]: w12n[mm];
    default: return 1;
  }
}

static GEN
ellsympow_goodred3(GEN E, GEN F, GEN p, long m, long vN, long *cnd, long *w)
{
  long o = ellsympow_inertia3(E, vN);
  long bet = ellsympow_betam(o, m);
  *cnd = m + 1 - bet + ellsympow_deltam3(o, m, vN);
  *w = odd(m)? ellsympow_rootno3(E, p, o, m): 1;
  if (o==1 || o==2)
    return ellsympow_abelian(p, ellap(F, p), m, o);
  if ((o==3 || o==6) && ellsympow_isabelian3(F))
    return ellsympow_abelian(p, p, m, o);
  else
    return ellsympow_nonabelian(p, m, bet);
}

static long
ellsympow_inertia2(GEN F, long vN)
{
  long vM = itos(gel(elllocalred(F, gen_2),1));
  GEN c6 = ell_get_c6(F);
  long v6 = signe(c6) ? vali(c6): 24;
  if (vM==0) return vN==0 ? 1: 2;
  if (vM==2) return vN==2 ? 3: 6;
  if (vM==5) return 8;
  if (vM==8) return v6>=9? 8: 4;
  if (vM==3 || vN==7) return 24;
  return 0;
}

static long
ellsympow_deltam2(long o, long m, long vN)
{
  if ((o==2 || o==6) && vN==4) return ellsympow_epsm(2, m);
  if ((o==2 || o==6) && vN==6) return 2*ellsympow_epsm(2, m);
  if (o==4) return 2*ellsympow_epsm(4, m)+ellsympow_epsm(2, m);
  if (o==8 && vN==5) return ellsympow_epsm(8, m)+ellsympow_epsm(2, m)/2;
  if (o==8 && vN==6) return ellsympow_epsm(8, m)+ellsympow_epsm(2, m);
  if (o==8 && vN==8) return ellsympow_epsm(8, m)+ellsympow_epsm(4, m)+ellsympow_epsm(2, m);
  if (o==24 && vN==3) return (2*ellsympow_epsm(8, m)+ellsympow_epsm(2, m))/6;
  if (o==24 && vN==4) return (ellsympow_epsm(8, m)+ellsympow_epsm(2, m)*2)/3;
  if (o==24 && vN==6) return (ellsympow_epsm(8, m)+ellsympow_epsm(2, m)*5)/3;
  if (o==24 && vN==7) return (ellsympow_epsm(8, m)*10+ellsympow_epsm(2, m)*5)/6;
  return 0;
}

static long
ellsympow_isabelian2(GEN F)
{ return umodi2n(ell_get_c4(F),7) == 96; }

static long
ellsympow_rootno2(GEN E, long vN, long m, long bet)
{
  long eps2 = (m + 1 - bet)>>1;
  long eta = odd(vN) && m%8==3 ? -1 : 1;
  long w2 = odd(eps2) ? ellrootno(E, gen_2): 1;
  return eta == w2 ? 1 : -1;
}

static GEN
ellsympow_goodred2(GEN E, GEN F, GEN p, long m, long vN, long *cnd, long *w)
{
  long o = ellsympow_inertia2(F, vN);
  long bet = ellsympow_betam(o, m);
  *cnd = m + 1 - bet + ellsympow_deltam2(o, m, vN);
  *w = odd(m) ? ellsympow_rootno2(E, vN, m, bet): 1;
  if (o==1 || o==2)
    return ellsympow_abelian(p, ellap(F, p), m, o);
  if (o==4 && ellsympow_isabelian2(F))
    return ellsympow_abelian(p, p, m, o);
  else
    return ellsympow_nonabelian(p, m, bet);
}

static GEN
ellminimaldotwist(GEN E, GEN *pD)
{
  GEN D = ellminimaltwistcond(E), Et = elltwist(E, D), Etmin;
  if (pD) *pD = D;
  Etmin = ellminimalmodel(Et, NULL);
  obj_free(Et); return Etmin;
}

/* Based on
Symmetric powers of elliptic curve L-functions,
Phil Martin and Mark Watkins, ANTS VII
<http://magma.maths.usyd.edu.au/users/watkins/papers/antsVII.pdf>
with thanks to Mark Watkins. BA20180402
*/
static GEN
lfunellsympow(GEN e, ulong m)
{
  pari_sp av = avma;
  GEN B, N, Nfa, pr, ex, ld, bad, ejd, et, pole;
  long i, l, mero, w = (m&7)==1 || (m&7)==3 ? -1: 1;
  checkell_Q(e);
  e = ellminimalmodel(e, NULL);
  ejd = Q_denom(ell_get_j(e));
  mero = m==0 || (m%4==0 && ellQ_get_CM(e)<0);
  ellQ_get_Nfa(e, &N, &Nfa);
  pr = gel(Nfa,1);
  ex = gel(Nfa,2); l = lg(pr);
  if (ugcd(umodiu(N,6), 6) == 1)
    et = NULL;
  else
    et = ellminimaldotwist(e, NULL);
  B = gen_1;
  bad = cgetg(l, t_VEC);
  for (i=1; i<l; i++)
  {
    long vN = itos(gel(ex,i));
    GEN p = gel(pr,i), eul;
    long cnd, wp;
    if (dvdii(ejd, p))
      eul = ellsympow_multred(e, p, m, vN, &cnd, &wp);
    else if (equaliu(p, 2))
      eul = ellsympow_goodred2(e, et, p, m, vN, &cnd, &wp);
    else if (equaliu(p, 3))
      eul = ellsympow_goodred3(e, et, p, m, vN, &cnd, &wp);
    else
      eul = ellsympow_goodred(e, p, m, &cnd, &wp);
    gel(bad, i) = mkvec2(p, ginv(eul));
    B = mulii(B, powiu(p,cnd));
    w *= wp;
  }
  pole = mero ? mkvec(mkvec2(stoi(1+(m>>1)),gen_0)): NULL;
  ld = mkvecn(mero? 7: 6, tag(mkvec2(mkvec2(e,utoi(m)),bad), t_LFUN_SYMPOW_ELL),
        gen_0, ellsympow_gamma(m), stoi(m+1), B, stoi(w), pole);
  if (et) obj_free(et);
  return gerepilecopy(av, ld);
}

GEN
lfunsympow(GEN ldata, ulong m)
{
  ldata = lfunmisc_to_ldata_shallow(ldata);
  if (ldata_get_type(ldata) != t_LFUN_ELL)
    pari_err_IMPL("lfunsympow");
  return lfunellsympow(gel(ldata_get_an(ldata), 2), m);
}

static GEN
lfunmfspec_i(GEN lmisc, long bit)
{
  GEN linit, ldataf, v, ve, vo, om, op, B, dom;
  long k, k2, j;

  ldataf = lfunmisc_to_ldata_shallow(lmisc);
  if (!gequal(ldata_get_gammavec(ldataf), mkvec2(gen_0,gen_1)))
    pari_err_TYPE("lfunmfspec", lmisc);
  k = gtos(ldata_get_k(ldataf));
  if (k == 1) return mkvec2(cgetg(1, t_VEC), gen_1);
  dom = mkvec3(dbltor(k/2.), dbltor((k-2)/2.), gen_0);
  if (is_linit(lmisc) && linit_get_type(lmisc) == t_LDESC_INIT
      && sdomain_isincl((double)k, dom, lfun_get_dom(linit_get_tech(lmisc))))
    linit = lmisc;
  else
    linit = lfuninit(ldataf, dom, 0, bit);
  B = int2n(bit/4);
  v = cgetg(k, t_VEC);
  for (j = 1; j < k; j++) gel(v,j) = lfunlambda(linit, utoi(j), bit);
  om = gel(v,1);
  if (odd(k)) return mkvec2(bestappr(gdiv(v, om), B), om);

  k2 = k/2;
  ve = cgetg(k2, t_VEC);
  vo = cgetg(k2+1, t_VEC);
  gel(vo,1) = om;
  for (j = 1; j < k2; j++)
  {
    gel(ve,j) = gel(v,2*j);
    gel(vo,j+1) = gel(v,2*j+1);
  }
  if (k2 == 1) { om = gen_1;    op = gel(v,1); }
  else         { om = gel(v,2); op = gel(v,3); }
  if (maxss(gexpo(imag_i(om)), gexpo(imag_i(op))) > -bit/2)
    pari_err_TYPE("lfunmfspec", lmisc);
  ve = gdiv(ve, om);
  vo = gdiv(vo, op);
  return mkvec4(bestappr(ve,B), bestappr(vo,B), om, op);
}
GEN
lfunmfspec(GEN lmisc, long bit)
{
  pari_sp av = avma;
  return gerepilecopy(av, lfunmfspec_i(lmisc, bit));
}

static long
ellsymsq_bad2(GEN c4, GEN c6, long e)
{
  switch (e)
  {
    case 2: return 1;
    case 3: return 0;
    case 5: return 0;
    case 7: return 0;
    case 8:
      if (!umodi2n(c6,9)) return 0;
      return umodi2n(c4,7)==32 ? 1 : -1;
    default: return 0;
  }
}
static long
ellsymsq_bad3(GEN c4, GEN c6, long e)
{
  long c6_243, c4_81;
  switch (e)
  {
    case 2: return 1;
    case 3: return 0;
    case 5: return 0;
    case 4:
      c4_81 = umodiu(c4,81);
      if (c4_81 == 27) return -1;
      if (c4_81%27 != 9) return 1;
      c6_243 = umodiu(c6,243);
      return (c6_243==108 || c6_243==135)? -1: 1;
    default: return 0;
  }
}
static int
c4c6_testp(GEN c4, GEN c6, GEN p)
{ GEN p2 = sqri(p); return (dvdii(c6,p2) && !dvdii(c4,p2)); }
/* assume e = v_p(N) >= 2 */
static long
ellsymsq_badp(GEN c4, GEN c6, GEN p, long e)
{
  if (absequaliu(p, 2)) return ellsymsq_bad2(c4, c6, e);
  if (absequaliu(p, 3)) return ellsymsq_bad3(c4, c6, e);
  switch(umodiu(p, 12UL))
  {
    case 1: return -1;
    case 5: return c4c6_testp(c4,c6,p)? -1: 1;
    case 7: return c4c6_testp(c4,c6,p)?  1:-1;
    default:return 1; /* p%12 = 11 */
  }
}
static GEN
lfunellsymsqmintwist(GEN e)
{
  pari_sp av = avma;
  GEN N, Nfa, P, E, V, c4, c6, ld;
  long i, l, k;
  checkell_Q(e);
  e = ellminimalmodel(e, NULL);
  ellQ_get_Nfa(e, &N, &Nfa);
  c4 = ell_get_c4(e);
  c6 = ell_get_c6(e);
  P = gel(Nfa,1); l = lg(P);
  E = gel(Nfa,2);
  V = cgetg(l, t_VEC);
  for (i=k=1; i<l; i++)
  {
    GEN p = gel(P,i);
    long a, e = itos(gel(E,i));
    if (e == 1) continue;
    a = ellsymsq_badp(c4, c6, p, e);
    gel(V,k++) = mkvec2(p, stoi(a));
  }
  setlg(V, k);
  ld = lfunellsympow(e, 2);
  return gerepilecopy(av, mkvec2(ld, V));
}

static GEN
mfpeters(GEN ldata2, GEN fudge, GEN N, long k, long bitprec)
{
  GEN t, L = real_i(lfun(ldata2, stoi(k), bitprec));
  long prec = nbits2prec(bitprec);
  t = powrs(mppi(prec), k+1); shiftr_inplace(t, 2*k-1); /* Pi/2 * (4Pi)^k */
  return gmul(gdiv(gmul(mulii(N,mpfact(k-1)), fudge), t), L);
}

/* Assume E to be twist-minimal */
static GEN
lfunellmfpetersmintwist(GEN E, long bitprec)
{
  pari_sp av = avma;
  GEN symsq, veceuler, N = ellQ_get_N(E), fudge = gen_1;
  long j, k = 2;
  symsq = lfunellsymsqmintwist(E);
  veceuler = gel(symsq,2);
  for (j = 1; j < lg(veceuler); j++)
  {
    GEN v = gel(veceuler,j), p = gel(v,1), q = powis(p,1-k);
    long s = signe(gel(v,2));
    if (s) fudge = gmul(fudge, s==1 ? gaddsg(1, q): gsubsg(1, q));
  }
  return gerepileupto(av, mfpeters(gel(symsq,1),fudge,N,k,bitprec));
}

/* From Christophe Delaunay, http://delaunay.perso.math.cnrs.fr/these.pdf */
static GEN
elldiscfix(GEN E, GEN Et, GEN D)
{
  GEN N = ellQ_get_N(E), Nt = ellQ_get_N(Et);
  GEN P = gel(absZ_factor(D), 1);
  GEN f = gen_1;
  long i, l = lg(P);
  for (i=1; i < l; i++)
  {
    GEN r, p = gel(P,i);
    long v = Z_pval(N, p), vt = Z_pval(Nt, p);
    if (v <= vt) continue;
    /* v > vt */
    if (absequaliu(p, 2))
    {
      if (vt == 0 && v >= 4)
        r = shifti(subsi(9, sqri(ellap(Et, p))), v-3);  /* 9=(2+1)^2 */
      else if (vt == 1)
        r = gmul2n(utoipos(3), v-3);  /* not in Z if v=2 */
      else if (vt >= 2)
        r = int2n(v-vt);
      else
        r = gen_1; /* vt = 0, 1 <= v <= 3 */
    }
    else if (vt >= 1)
      r = gdiv(subiu(sqri(p), 1), p);
    else
      r = gdiv(mulii(subiu(p, 1), subii(sqri(addiu(p, 1)), sqri(ellap(Et, p)))), p);
    f = gmul(f, r);
  }
  return f;
}

GEN
lfunellmfpeters(GEN E, long bitprec)
{
  pari_sp ltop = avma;
  GEN D, Et = ellminimaldotwist(E, &D);
  GEN nor = lfunellmfpetersmintwist(Et, bitprec);
  GEN nor2 = gmul(nor, elldiscfix(E, Et, D));
  obj_free(Et); return gerepileupto(ltop, nor2);
}

/*************************************************************/
/*               Genus 2 curves                              */
/*************************************************************/

static void
Flv_diffnext(GEN d, ulong p)
{
  long j, n = lg(d)-1;
  for(j = n; j>=2; j--)
    uel(d,j) = Fl_add(uel(d,j), uel(d,j-1), p);
}

static GEN
Flx_difftable(GEN P, ulong p)
{
  long i, n = degpol(P);
  GEN V = cgetg(n+2, t_VECSMALL);
  uel(V, n+1) = Flx_constant(P);
  for(i = n; i >= 1; i--)
  {
    P = Flx_diff1(P, p);
    uel(V, i) = Flx_constant(P);
  }
  return V;
}

static long
Flx_genus2trace_naive(GEN H, ulong p)
{
  pari_sp av = avma;
  ulong i, j;
  long a, n = degpol(H);
  GEN k = const_vecsmall(p, -1), d;
  k[1] = 0;
  for (i=1, j=1; i < p; i += 2, j = Fl_add(j, i, p))
    k[j+1] = 1;
  a = n == 5 ? 0: k[1+Flx_lead(H)];
  d = Flx_difftable(H, p);
  for (i=0; i < p; i++)
  {
    a += k[1+uel(d,n+1)];
    Flv_diffnext(d, p);
  }
  set_avma(av);
  return a;
}

static GEN
dirgenus2(GEN Q, GEN p, long n)
{
  pari_sp av = avma;
  GEN f;
  if (n > 2)
    f = RgX_recip(hyperellcharpoly(gmul(Q,gmodulo(gen_1, p))));
  else
  {
    ulong pp = itou(p);
    GEN Qp = ZX_to_Flx(Q, pp);
    long t = Flx_genus2trace_naive(Qp, pp);
    f = deg1pol_shallow(stoi(t), gen_1, 0);
  }
  return gerepileupto(av, RgXn_inv_i(f, n));
}

GEN
dirgenus2_worker(GEN P, ulong X, GEN Q)
{
  pari_sp av = avma;
  long i, l = lg(P);
  GEN V = cgetg(l, t_VEC);
  for(i = 1; i < l; i++)
  {
    ulong p = uel(P,i);
    long d = ulogint(X, p) + 1; /* minimal d such that p^d > X */
    gel(V,i) = dirgenus2(Q, utoi(uel(P,i)), d);
  }
  return gerepilecopy(av, mkvec2(P,V));
}

static GEN
vecan_genus2(GEN an, long L)
{
  GEN Q = gel(an,1), bad = gel(an, 2);
  GEN worker = snm_closure(is_entry("_dirgenus2_worker"), mkvec(Q));
  return pardireuler(worker, gen_2, stoi(L), NULL, bad);
}

static GEN
eulerf_genus2(GEN an, GEN p)
{
  GEN Q = gel(an,1), bad = gel(an, 2);
  GEN f = eulerf_bad(bad, p);
  if (f) return f;
  f = RgX_recip(hyperellcharpoly(gmul(Q,gmodulo(gen_1, p))));
  return mkrfrac(gen_1,f);
}

static GEN
genus2_redmodel(GEN P, GEN p)
{
  GEN M = FpX_factor(P, p);
  GEN F = gel(M,1), E = gel(M,2);
  long i, k, r = lg(F);
  GEN U = scalarpol(leading_coeff(P), varn(P));
  GEN G = cgetg(r, t_COL);
  for (i=1, k=0; i<r; i++)
  {
    if (E[i]>1)
      gel(G,++k) = gel(F,i);
    if (odd(E[i]))
      U = FpX_mul(U, gel(F,i), p);
  }
  setlg(G,++k);
  return mkvec2(G,U);
}

static GEN
oneminusxd(long d)
{
  return gsub(gen_1, pol_xn(d, 0));
}

static GEN
ellfromeqncharpoly(GEN P, GEN Q, GEN p)
{
  long v;
  GEN E, F, t, y;
  v = fetch_var();
  y = pol_x(v);
  F = gsub(gadd(ZX_sqr(y), gmul(y, Q)), P);
  E = ellinit(ellfromeqn(F), p, DEFAULTPREC);
  delete_var();
  t = ellap(E, p);
  obj_free(E);
  return mkpoln(3, gen_1, negi(t), p);
}

static GEN
genus2_eulerfact(GEN P, GEN p)
{
  GEN Pp = FpX_red(P, p);
  GEN GU = genus2_redmodel(Pp, p);
  long d = 6-degpol(Pp), v = d/2, w = odd(d);
  GEN abe, tor;
  GEN ki, kp = pol_1(0), kq = pol_1(0);
  GEN F = gel(GU,1), Q = gel(GU,2);
  long dQ = degpol(Q), lF = lg(F)-1;

  abe = dQ >= 5 ? RgX_recip(hyperellcharpoly(gmul(Q,gmodulo(gen_1,p))))
      : dQ >= 3 ? RgX_recip(ellfromeqncharpoly(Q,gen_0,p))
                : pol_1(0);
  ki = dQ != 0 ? oneminusxd(1)
              : Fp_issquare(gel(Q,2),p) ? ZX_sqr(oneminusxd(1))
                                        : oneminusxd(2);
  if (lF)
  {
    long i;
    for(i=1; i <= lF; i++)
    {
      GEN Fi = gel(F, i);
      long d = degpol(Fi);
      GEN e = FpX_rem(Q, Fi, p);
      GEN kqf = lgpol(e)==0 ? oneminusxd(d):
                FpXQ_issquare(e, Fi, p) ? ZX_sqr(oneminusxd(d))
                                        : oneminusxd(2*d);
      kp = gmul(kp, oneminusxd(d));
      kq = gmul(kq, kqf);
    }
  }
  if (v)
  {
    GEN kqoo = w==1 ? oneminusxd(1):
               Fp_issquare(leading_coeff(Q), p)? ZX_sqr(oneminusxd(1))
                                              : oneminusxd(2);
    kp = gmul(kp, oneminusxd(1));
    kq = gmul(kq, kqoo);
  }
  tor = RgX_div(ZX_mul(oneminusxd(1), kq), ZX_mul(ki, kp));
  return ginv( ZX_mul(abe, tor) );
}

static GEN
F2x_genus2_find_trans(GEN P, GEN Q, GEN F)
{
  pari_sp av = avma;
  long i, d = F2x_degree(F), v = P[1];
  GEN M, C, V;
  M = cgetg(d+1, t_MAT);
  for (i=1; i<=d; i++)
  {
    GEN Mi = F2x_rem(F2x_add(F2x_shift(Q,i-1), monomial_F2x(2*i-2,v)), F);
    gel(M,i) = F2x_to_F2v(Mi, d);
  }
  C = F2x_to_F2v(F2x_rem(P, F), d);
  V = F2m_F2c_invimage(M, C);
  return gerepileuptoleaf(av, F2v_to_F2x(V, v));
}

static GEN
F2x_genus2_trans(GEN P, GEN Q, GEN H)
{
  return F2x_add(P,F2x_add(F2x_mul(H,Q), F2x_sqr(H)));
}

static GEN
F2x_genus_redoo(GEN P, GEN Q, long k)
{
  if (F2x_degree(P)==2*k)
  {
    long c = F2x_coeff(P,2*k-1), dQ = F2x_degree(Q);
    if ((dQ==k-1 && c==1) || (dQ<k-1 && c==0))
     return F2x_genus2_trans(P, Q, monomial_F2x(k, P[1]));
  }
  return P;
}

static GEN
F2x_pseudodisc(GEN P, GEN Q)
{
  GEN dP = F2x_deriv(P), dQ = F2x_deriv(Q);
  return F2x_gcd(Q, F2x_add(F2x_mul(P, F2x_sqr(dQ)), F2x_sqr(dP)));
}

static GEN
F2x_genus_red(GEN P, GEN Q)
{
  long dP, dQ;
  GEN F, FF;
  P = F2x_genus_redoo(P, Q, 3);
  P = F2x_genus_redoo(P, Q, 2);
  P = F2x_genus_redoo(P, Q, 1);
  dP = F2x_degree(P);
  dQ = F2x_degree(Q);
  FF = F = F2x_pseudodisc(P,Q);
  while(F2x_degree(F)>0)
  {
    GEN M = gel(F2x_factor(F),1);
    long i, l = lg(M);
    for(i=1; i<l; i++)
    {
      GEN R = F2x_sqr(gel(M,i));
      GEN H = F2x_genus2_find_trans(P, Q, R);
      P = F2x_div(F2x_genus2_trans(P, Q, H), R);
      Q = F2x_div(Q, gel(M,i));
    }
    F = F2x_pseudodisc(P, Q);
  }
  return mkvec4(P,Q,FF,mkvecsmall2(dP,dQ));
}

/* Number of solutions of x^2+b*x+c */
static long
F2xqX_quad_nbroots(GEN b, GEN c, GEN T)
{
  if (lgpol(b) > 0)
  {
    GEN d = F2xq_div(c, F2xq_sqr(b, T), T);
    return F2xq_trace(d, T)? 0: 2;
  }
  else
    return 1;
}

static GEN
genus2_eulerfact2(GEN PQ)
{
  GEN V = F2x_genus_red(ZX_to_F2x(gel(PQ, 1)), ZX_to_F2x(gel(PQ, 2)));
  GEN P = gel(V, 1), Q = gel(V, 2);
  GEN F = gel(V, 3), v = gel(V, 4);
  GEN abe, tor;
  GEN ki, kp = pol_1(0), kq = pol_1(0);
  long dP = F2x_degree(P), dQ = F2x_degree(Q), d = maxss(dP, 2*dQ);
  if (!lgpol(F)) return pol_1(0);
  ki = dQ!=0 || dP>0 ? oneminusxd(1):
      dP==-1 ? ZX_sqr(oneminusxd(1)): oneminusxd(2);
  abe = d>=5? RgX_recip(hyperellcharpoly(gmul(PQ,gmodulss(1,2)))):
        d>=3? RgX_recip(ellfromeqncharpoly(F2x_to_ZX(P), F2x_to_ZX(Q), gen_2)):
        pol_1(0);
  if (lgpol(F))
  {
    GEN M = gel(F2x_factor(F), 1);
    long i, lF = lg(M)-1;
    for(i=1; i <= lF; i++)
    {
      GEN Fi = gel(M, i);
      long d = F2x_degree(Fi);
      long nb  = F2xqX_quad_nbroots(F2x_rem(Q, Fi), F2x_rem(P, Fi), Fi);
      GEN kqf = nb==1 ? oneminusxd(d):
                nb==2 ? ZX_sqr(oneminusxd(d))
                      : oneminusxd(2*d);
      kp = gmul(kp, oneminusxd(d));
      kq = gmul(kq, kqf);
    }
  }
  if (maxss(v[1],2*v[2])<5)
  {
    GEN kqoo = v[1]>2*v[2] ? oneminusxd(1):
               v[1]<2*v[2] ? ZX_sqr(oneminusxd(1))
                           : oneminusxd(2);
    kp = gmul(kp, oneminusxd(1));
    kq = gmul(kq, kqoo);
  }
  tor = RgX_div(ZX_mul(oneminusxd(1),kq), ZX_mul(ki, kp));
  return ginv( ZX_mul(abe, tor) );
}

GEN
lfungenus2(GEN G)
{
  pari_sp ltop = avma;
  GEN Ldata;
  GEN gr = genus2red(G, NULL);
  GEN N  = gel(gr, 1), M = gel(gr, 2), PQ = gel(gr, 3), L = gel(gr, 4);
  GEN e, F = gadd(gsqr(gel(PQ, 2)), gmul2n(gel(PQ, 1), 2));
  long i, lL = lg(L), ram2;
  ram2 = absequaliu(gmael(M,1,1),2);
  if (ram2 && equalis(gmael(M,2,1),-1))
    pari_warn(warner,"unknown valuation of conductor at 2");
  e = cgetg(lL+(ram2?0:1), t_VEC);
  gel(e,1) = mkvec2(gen_2, ram2 ? genus2_eulerfact2(PQ)
           : ginv( RgX_recip(hyperellcharpoly(gmul(PQ,gmodulss(1,2))))) );
  for(i = ram2? 2: 1; i < lL; i++)
  {
    GEN Li = gel(L, i);
    GEN p = gel(Li, 1);
    gel(e, ram2 ? i: i+1) = mkvec2(p, genus2_eulerfact(F,p));
  }
  Ldata = mkvecn(6, tag(mkvec2(F,e), t_LFUN_GENUS2),
      gen_0, mkvec4(gen_0, gen_0, gen_1, gen_1), gen_2, N, gen_0);
  return gerepilecopy(ltop, Ldata);
}

/*************************************************************/
/*                        ETA QUOTIENTS                      */
/* An eta quotient is a matrix with 2 columns [m, r_m] with  */
/* m >= 1 representing f(\tau)=\prod_m\eta(m\tau)^{r_m}.     */
/*************************************************************/

/* eta(x^v) + O(x^L) */
GEN
eta_ZXn(long v, long L)
{
  long n, k = 0, v2 = 2*v, bn = v, cn = 0;
  GEN P;
  if (!L) return zeropol(0);
  P = cgetg(L+2,t_POL); P[1] = evalsigne(1);
  for(n = 0; n < L; n++) gel(P,n+2) = gen_0;
  for(n = 0;; n++, bn += v2, cn += v)
  { /* k = v * (3*n-1) * n / 2; bn = v * (2*n+1); cn = v * n */
    long k2;
    gel(P, k+2) = odd(n)? gen_m1: gen_1;
    k2 = k+cn; if (k2 >= L) break;
    k = k2;
    /* k = v * (3*n+1) * n / 2 */;
    gel(P, k+2) = odd(n)? gen_m1: gen_1;
    k2 = k+bn; if (k2 >= L) break;
    k = k2;
  }
  setlg(P, k+3); return P;
}
GEN
eta_product_ZXn(GEN eta, long L)
{
  pari_sp av = avma;
  GEN P = NULL, D = gel(eta,1), R = gel(eta,2);
  long i, l = lg(D);
  for (i = 1; i < l; ++i)
  {
    GEN Q = eta_ZXn(D[i], L);
    long r = R[i];
    if (r < 0) { Q = RgXn_inv_i(Q, L); r = -r; }
    if (r != 1) Q = RgXn_powu_i(Q, r, L);
    P = P? ZXn_mul(P, Q, L): Q;
    if (gc_needed(av,1) && i > 1)
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"eta_product_ZXn");
      P = gerepilecopy(av, P);
    }
  }
  return P;
}
static GEN
vecan_eta(GEN an, long L)
{
  long v = itos(gel(an, 3));
  GEN t;
  if (v > L) return zerovec(L);
  t = eta_product_ZXn(an, L - v);
  if (v) t = RgX_shift_shallow(t, v);
  return RgX_to_RgV(t, L);
}
/* return 1 if cuspidal, 0 if holomorphic, -1 otherwise */
static int
etacuspidal(GEN N, GEN k, GEN B, GEN R, GEN NB)
{
  long i, j, lD, l, cusp = 1;
  pari_sp av = avma;
  GEN D;
  if (gsigne(k) < 0) return -1;
  D = divisors(N); lD = lg(D); l = lg(B);
  for (i = 1; i < lD; i++)
  {
    GEN t = gen_0, d = gel(D,i);
    long s;
    for (j = 1; j < l; j++)
      t = addii(t, mulii(gel(NB,j), mulii(gel(R,j), sqri(gcdii(d, gel(B,j))))));
    s = signe(t);
    if (s < 0) return -1;
    if (!s) cusp = 0;
  }
  return gc_bool(av, cusp);
}
/* u | 24, level N = u*N0, N0 = lcm(B), NB[i] = N0/B[i] */
static int
etaselfdual(GEN B, GEN R, GEN NB, ulong u)
{
  pari_sp av = avma;
  long i, l = lg(B);
  for (i = 1; i < l; i++)
  {
    long j = ZV_search(B, muliu(gel(NB,i), u)); /* search for N / B[i] */
    set_avma(av); if (!j || !equalii(gel(R,i),gel(R,j))) return 0;
  }
  return 1;
}
/* return Nebentypus of eta quotient, k2 = 2*k integral */
static GEN
etachar(GEN B, GEN R, GEN k2)
{
  long i, l = lg(B);
  GEN P = gen_1;
  for (i = 1; i < l; ++i) if (mpodd(gel(R,i))) P = mulii(P, gel(B,i));
  switch(Mod4(k2))
  {
    case 0: break;
    case 2:  P = negi(P); break;
    default: P = shifti(P, 1); break;
  }
  return coredisc(P);
}
/* Return 0 if not on gamma_0(N). Sets conductor, modular weight, character,
 * canonical matrix, v_q(eta), sd = 1 iff self-dual, cusp = 1 iff cuspidal
 * [0 if holomorphic at all cusps, else -1] */
long
etaquotype(GEN *peta, GEN *pN, GEN *pk, GEN *CHI, long *pv, long *sd,
           long *cusp)
{
  GEN B, R, S, T, N, NB, eta = *peta;
  long l, i, u, S24;

  if (lg(eta) != 3) pari_err_TYPE("lfunetaquo", eta);
  switch(typ(eta))
  {
    case t_VEC: eta = mkmat2(mkcol(gel(eta,1)), mkcol(gel(eta,2))); break;
    case t_MAT: break;
    default: pari_err_TYPE("lfunetaquo", eta);
  }
  if (!RgV_is_ZVpos(gel(eta,1)) || !RgV_is_ZV(gel(eta,2)))
    pari_err_TYPE("lfunetaquo", eta);
  *peta = eta = famat_reduce(eta);
  B = gel(eta,1); l = lg(B); /* sorted in increasing order */
  R = gel(eta,2);
  N = ZV_lcm(B); NB = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(NB,i) = diviiexact(N, gel(B,i));
  S = gen_0; T = gen_0; u = 0;
  for (i = 1; i < l; ++i)
  {
    GEN b = gel(B,i), r = gel(R,i);
    S = addii(S, mulii(b, r));
    T = addii(T, r);
    u += umodiu(r,24) * umodiu(gel(NB,i), 24);
  }
  S = divis_rem(S, 24, &S24);
  if (S24) return 0; /* nonintegral valuation at oo */
  u = 24 / ugcd(24, u % 24);
  *pN = muliu(N, u); /* level */
  *pk = gmul2n(T,-1); /* weight */
  *pv = itos(S); /* valuation */
  if (cusp) *cusp = etacuspidal(*pN, *pk, B, R, NB);
  if (sd) *sd = etaselfdual(B, R, NB, u);
  if (CHI) *CHI = etachar(B, R, T);
  return 1;
}

GEN
lfunetaquo(GEN eta0)
{
  pari_sp ltop = avma;
  GEN Ldata, N, BR, k, eta = eta0;
  long v, sd, cusp;
  if (!etaquotype(&eta, &N, &k, NULL, &v, &sd, &cusp))
    pari_err_TYPE("lfunetaquo", eta0);
  if (!cusp) pari_err_IMPL("noncuspidal eta quotient");
  if (!sd) pari_err_IMPL("non self-dual eta quotient");
  if (typ(k) != t_INT) pari_err_TYPE("lfunetaquo [nonintegral weight]", eta0);
  BR = mkvec3(ZV_to_zv(gel(eta,1)), ZV_to_zv(gel(eta,2)), stoi(v - 1));
  Ldata = mkvecn(6, tag(BR,t_LFUN_ETA), gen_0, mkvec2(gen_0,gen_1), k,N, gen_1);
  return gerepilecopy(ltop, Ldata);
}

static GEN
vecan_qf(GEN Q, long L)
{
  GEN v, w = qfrep0(Q, utoi(L), 1);
  long i;
  v = cgetg(L+1, t_VEC);
  for (i = 1; i <= L; i++) gel(v,i) = utoi(2 * w[i]);
  return v;
}

long
qfiseven(GEN M)
{
  long i, l = lg(M);
  for (i=1; i<l; i++)
    if (mpodd(gcoeff(M,i,i))) return 0;
  return 1;
}

GEN
lfunqf(GEN M, long prec)
{
  pari_sp ltop = avma;
  long n;
  GEN k, D, d, Mi, Ldata, poles, eno, dual;

  if (typ(M) != t_MAT) pari_err_TYPE("lfunqf", M);
  if (!RgM_is_ZM(M))   pari_err_TYPE("lfunqf [not integral]", M);
  n = lg(M)-1;
  k = uutoQ(n,2);
  M = Q_primpart(M);
  Mi = ZM_inv(M, &d); /* d M^(-1) */
  if (!qfiseven(M)) { M = gmul2n(M, 1); d = shifti(d,1); }
  if (!qfiseven(Mi)){ Mi= gmul2n(Mi,1); d = shifti(d,1); }
  /* det(Mi) = d^n/det(M), D^2 = det(Mi)/det(M) */
  D = gdiv(gpow(d,k,prec), ZM_det(M));
  if (!issquareall(D, &eno)) eno = gsqrt(D, prec);
  dual = gequal1(D) ? gen_0: tag(Mi, t_LFUN_QF);
  poles = mkcol2(mkvec2(k, simple_pole(gmul2n(eno,1))),
                 mkvec2(gen_0, simple_pole(gen_m2)));
  Ldata = mkvecn(7, tag(M, t_LFUN_QF), dual,
       mkvec2(gen_0, gen_1), k, d, eno, poles);
  return gerepilecopy(ltop, Ldata);
}

/********************************************************************/
/**  Artin L function, based on a GP script by Charlotte Euvrard   **/
/********************************************************************/

static GEN
artin_charfromgens(GEN G, GEN M)
{
  GEN R, V, ord = gal_get_orders(G), grp = gal_get_group(G);
  long i, j, k, n = lg(ord)-1, m = lg(grp)-1;

  if (lg(M)-1 != n) pari_err_DIM("lfunartin");
  R = cgetg(m+1, t_VEC);
  gel(R, 1) = matid(lg(gel(M, 1))-1);
  for (i = 1, k = 1; i <= n; ++i)
  {
    long c = k*(ord[i] - 1);
    gel(R, ++k) = gel(M, i);
    for (j = 2; j <= c; ++j) gel(R, ++k) = gmul(gel(R,j), gel(M,i));
  }
  V = cgetg(m+1, t_VEC);
  for (i = 1; i <= m; i++) gel(V, gel(grp,i)[1]) = gtrace(gel(R,i));
  return V;
}

/* TODO move somewhere else? */
GEN
galois_get_conj(GEN G)
{
  GEN grp = gal_get_group(G);
  long i, k, r = lg(grp)-1;
  GEN b = zero_F2v(r);
  for (k = 2; k <= r; ++k)
  {
    GEN g = gel(grp,k);
    if (!F2v_coeff(b,g[1]) && g[g[1]]==1)
    {
      pari_sp av = avma;
      GEN F = galoisfixedfield(G, g, 1, -1);
      if (ZX_sturmpart(F, NULL) > 0) { set_avma(av); return g; }
      for (i = 1; i<=r; i++)
      {
        GEN h = gel(grp, i);
        long t = h[1];
        while (h[t]!=1) t = h[t];
        F2v_set(b, h[g[t]]);
      }
      set_avma(av);
    }
  }
  pari_err_BUG("galois_get_conj");
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN  cyclotoi(GEN v) { return simplify_shallow(lift_shallow(v)); }
static long cyclotos(GEN v) { return gtos(cyclotoi(v)); }
static long char_dim(GEN ch) { return cyclotos(gel(ch,1)); }

static GEN
artin_gamma(GEN N, GEN G, GEN ch)
{
  long a, t, d = char_dim(ch);
  if (nf_get_r2(N) == 0) return vec01(d, 0);
  a = galois_get_conj(G)[1];
  t = cyclotos(gel(ch,a));
  return vec01((d+t) / 2, (d-t) / 2);
}

static long
artin_dim(GEN ind, GEN ch)
{
  long n = lg(ch)-1;
  GEN elts = group_elts(ind, n);
  long i, d = lg(elts)-1;
  GEN s = gen_0;
  for(i=1; i<=d; i++)
    s = gadd(s, gel(ch, gel(elts,i)[1]));
  return gtos(gdivgu(cyclotoi(s), d));
}

static GEN
artin_ind(GEN elts, GEN ch, GEN p)
{
  long i, d = lg(elts)-1;
  GEN s = gen_0;
  for(i=1; i<=d; i++)
    s = gadd(s, gel(ch, gmul(gel(elts,i),p)[1]));
  return gdivgu(s, d);
}

static GEN
artin_ram(GEN nf, GEN gal, GEN aut, GEN pr, GEN ramg, GEN ch, long d)
{
  pari_sp av = avma;
  long i, v, n;
  GEN p, q, V, elts;
  if (d==0) return pol_1(0);
  n = degpol(gal_get_pol(gal));
  q = p = idealramfrobenius_aut(nf, gal, pr, ramg, aut);
  elts = group_elts(gel(ramg,2), n);
  v = fetch_var_higher();
  V = cgetg(d+2, t_POL);
  V[1] = evalsigne(1)|evalvarn(v);
  for(i=1; i<=d; i++)
  {
    gel(V,i+1) = artin_ind(elts, ch, q);
    q = gmul(q, p);
  }
  delete_var();
  V = RgXn_expint(RgX_neg(V),d+1);
  setvarn(V,0); return gerepileupto(av, ginv(V));
}

/* N true nf; [Artin conductor, vec of [p, Lp]] */
static GEN
artin_badprimes(GEN N, GEN G, GEN aut, GEN ch)
{
  pari_sp av = avma;
  long i, d = char_dim(ch);
  GEN P = gel(absZ_factor(nf_get_disc(N)), 1);
  long lP = lg(P);
  GEN B = cgetg(lP, t_VEC), C = cgetg(lP, t_VEC);

  for (i = 1; i < lP; ++i)
  {
    GEN p = gel(P, i), pr = idealprimedec_galois(N, p);
    GEN J = idealramgroups_aut(N, G, pr, aut);
    GEN G0 = gel(J,2); /* inertia group */
    long lJ = lg(J);
    long sdec = artin_dim(G0, ch);
    long ndec = group_order(G0);
    long j, v = ndec * (d - sdec);
    for (j = 3; j < lJ; ++j)
    {
      GEN Jj = gel(J, j);
      long s = artin_dim(Jj, ch);
      v += group_order(Jj) * (d - s);
    }
    gel(C, i) = powiu(p, v/ndec);
    gel(B, i) = mkvec2(p, artin_ram(N, G, aut, pr, J, ch, sdec));
  }
  return gerepilecopy(av, mkvec2(ZV_prod(C), B));
}

/* p does not divide nf.index */
static GEN
idealfrobenius_easy(GEN nf, GEN gal, GEN aut, GEN T, GEN p)
{
  long i, l = lg(aut), f = degpol(T);
  GEN D, Dzk, DzkT, DXp, grp = gal_get_group(gal);
  pari_sp av = avma;
  if (f==1) return gel(grp,1);
  Dzk = nf_get_zkprimpart(nf);
  D = modii(nf_get_zkden(nf), p);
  DzkT = RgV_to_RgM(FqV_red(Dzk, T, p), f);
  DXp = RgX_to_RgC(FpX_Frobenius(T, p), f);
  if (!equali1(D)) DXp = FpC_Fp_mul(DXp, D, p);
  for(i=1; i < l; i++)
  {
    GEN g = gel(grp,i);
    if (perm_orderu(g) == (ulong)f)
    {
      GEN A = FpM_FpC_mul(DzkT, gel(aut,g[1]), p);
      if (ZV_equal(A, DXp)) {set_avma(av); return g; }
    }
  }
  return NULL; /* LCOV_EXCL_LINE */
}
/* true nf; p divides nf.index, pr/p unramified */
static GEN
idealfrobenius_hard(GEN nf, GEN gal, GEN aut, GEN pr)
{
  long i, l = lg(aut), f = pr_get_f(pr);
  GEN modpr, p, T, X, Xp, pi, grp = gal_get_group(gal);
  pari_sp av = avma;
  if (f==1) return gel(grp,1);
  pi = pr_get_gen(pr);
  modpr = zkmodprinit(nf, pr);
  p = modpr_get_p(modpr);
  T = modpr_get_T(modpr);
  X = modpr_genFq(modpr);
  Xp = FpX_Frobenius(T, p);
  for (i = 1; i < l; i++)
  {
    GEN g = gel(grp,i);
    if (perm_orderu(g) == (ulong)f)
    {
      GEN S = gel(aut,g[1]);
      GEN A = nf_to_Fq(nf, zk_galoisapplymod(nf,X,S,p), modpr);
      /* sigma(X) = X^p (mod pr) and sigma(pi) in pr */
      if (ZX_equal(A, Xp) && (f == nf_get_degree(nf) ||
          ZC_prdvd(zk_galoisapplymod(nf,pi,S,p),pr))) { set_avma(av); return g; }
    }
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/* true nf */
static GEN
dirartin(GEN nf, GEN G, GEN V, GEN aut, GEN p, long n)
{
  pari_sp av = avma;
  GEN pr, frob;
  /* pick one maximal ideal in the conjugacy class above p */
  GEN T = nf_get_pol(nf);
  if (!dvdii(nf_get_index(nf), p))
  { /* simple case */
    GEN F = FpX_factor(T, p), P = gmael(F,1,1);
    frob = idealfrobenius_easy(nf, G, aut, P, p);
  }
  else
  {
    pr = idealprimedec_galois(nf,p);
    frob = idealfrobenius_hard(nf, G, aut, pr);
  }
  set_avma(av); return n ? RgXn_inv(gel(V, frob[1]), n): gel(V, frob[1]);
}

GEN
dirartin_worker(GEN P, ulong X, GEN nf, GEN G, GEN V, GEN aut)
{
  pari_sp av = avma;
  long i, l = lg(P);
  GEN W = cgetg(l, t_VEC);
  for(i = 1; i < l; i++)
  {
    ulong p = uel(P,i);
    long d = ulogint(X, p) + 1; /* minimal d such that p^d > X */
    gel(W,i) = dirartin(nf, G, V, aut, utoi(uel(P,i)), d);
  }
  return gerepilecopy(av, mkvec2(P,W));
}

static GEN
vecan_artin(GEN an, long L, long prec)
{
  GEN A, Sbad = gel(an,5);
  long n = itos(gel(an,6)), isreal = lg(an)<8 ? 0: !itos(gel(an,7));
  GEN worker = snm_closure(is_entry("_dirartin_worker"), vecslice(an,1,4));
  A = lift_shallow(pardireuler(worker, gen_2, stoi(L), NULL, Sbad));
  A = RgXV_RgV_eval(A, grootsof1(n, prec));
  if (isreal) A = real_i(A);
  return A;
}

static GEN
eulerf_artin(GEN an, GEN p, long prec)
{
  GEN nf = gel(an,1), G = gel(an,2), V = gel(an,3), aut = gel(an,4);
  GEN Sbad = gel(an,5);
  long n = itos(gel(an,6)), isreal = lg(an)<8 ? 0: !itos(gel(an,7));
  GEN f = eulerf_bad(Sbad, p);
  if (!f) f = mkrfrac(gen_1,dirartin(nf, G, V, aut, p, 0));
  f = gsubst(liftpol(f),1, rootsof1u_cx(n, prec));
  if (isreal) f = real_i(f);
  return f;
}

static GEN
char_expand(GEN conj, GEN ch)
{
  long i, l = lg(conj);
  GEN V = cgetg(l, t_COL);
  for (i=1; i<l; i++) gel(V,i) = gel(ch, conj[i]);
  return V;
}

static GEN
handle_zeta(long n, GEN ch, long *m)
{
  GEN c;
  long t, i, l = lg(ch);
  GEN dim = cyclotoi(vecsum(ch));
  if (typ(dim) != t_INT)
    pari_err_DOMAIN("lfunartin","chi","is not a", strtoGENstr("character"), ch);
  t = itos(dim);
  if (t < 0 || t % n)
    pari_err_DOMAIN("lfunartin","chi","is not a", strtoGENstr("character"), ch);
  if (t == 0) { *m = 0; return ch; }
  *m = t / n;
  c = cgetg(l, t_COL);
  for (i=1; i<l; i++)
    gel(c,i) = gsubgs(gel(ch,i), *m);
  return c;
}

static int
cyclo_is_real(GEN v, GEN ix)
{
  pari_sp av = avma;
  GEN w = poleval(lift_shallow(v), ix);
  return gc_bool(av, gequal(w, v));
}

static int
char_is_real(GEN ch, GEN mod)
{
  long i, l = lg(ch);
  GEN ix = QXQ_inv(pol_x(varn(mod)), mod);
  for (i=1; i<l; i++)
    if (!cyclo_is_real(gel(ch,i), ix)) return 0;
  return 1;
}

GEN
lfunartin(GEN nf, GEN gal, GEN ch, long o, long bitprec)
{
  pari_sp av = avma;
  GEN bc, V, aut, mod, Ldata = NULL, chx, cc, conj, repr;
  long tmult, var;
  nf = checknf(nf);
  checkgal(gal);
  var = gvar(ch);
  if (var == 0) pari_err_PRIORITY("lfunartin",ch,"=",0);
  if (var < 0) var = 1;
  if (!is_vec_t(typ(ch))) pari_err_TYPE("lfunartin", ch);
  cc = group_to_cc(gal);
  conj = gel(cc,2);
  repr = gel(cc,3);
  mod = mkpolmod(gen_1, polcyclo(o, var));
  if (lg(ch)>1 && typ(gel(ch,1))==t_MAT)
    chx = artin_charfromgens(gal, gmul(ch,mod));
  else
  {
    if (lg(repr) != lg(ch)) pari_err_DIM("lfunartin");
    chx = char_expand(conj, gmul(ch,mod));
  }
  chx = handle_zeta(nf_get_degree(nf), chx, &tmult);
  ch = shallowextract(chx, repr);
  if (!gequal0(chx))
  {
    GEN real = char_is_real(chx, gel(mod,1))? gen_0: gen_1;
    aut = nfgaloispermtobasis(nf, gal);
    V = gmul(char_expand(conj, galoischarpoly(gal, ch, o)), mod);
    bc = artin_badprimes(nf, gal, aut, chx);
    Ldata = mkvecn(6,
      tag(mkcoln(7, nf, gal, V, aut, gel(bc, 2), stoi(o), real), t_LFUN_ARTIN),
      real, artin_gamma(nf, gal, chx), gen_1, gel(bc,1), gen_0);
  }
  if (tmult==0 && Ldata==NULL) pari_err_TYPE("lfunartin",ch);
  if (tmult)
  {
    long i;
    if (Ldata==NULL) { Ldata = lfunzeta(); tmult--; }
    for(i=1; i<=tmult; i++)
      Ldata = lfunmul(Ldata, gen_1, bitprec);
  }
  return gerepilecopy(av, Ldata);
}

/* true nf */
static GEN
lfunzetakinit_artin(GEN nf, GEN gal, GEN dom, long der, long bitprec)
{
  GEN F, E, M, domain, To = galoischartable(gal), T = gel(To, 1);
  long i, o = itos(gel(To, 2)), l = lg(T);
  F = cgetg(l, t_VEC);
  E = cgetg(l, t_VECSMALL);
  for (i = 1; i < l; ++i)
  {
    GEN L = lfunartin(nf, gal, gel(T,i), o, bitprec);
    gel(F, i) = lfuninit(L, dom, der, bitprec);
    E[i] = char_dim(gel(T,i));
  }
  domain = mkvec2(dom, mkvecsmall2(der, bitprec));
  M = mkvec3(F, E, zero_zv(l-1));
  return lfuninit_make(t_LDESC_PRODUCT, lfunzetak_i(nf), M, domain);
}

/********************************************************************/
/**                    High-level Constructors                     **/
/********************************************************************/
enum { t_LFUNMISC_POL, t_LFUNMISC_CHIQUAD, t_LFUNMISC_CHICONREY,
       t_LFUNMISC_CHIGEN, t_LFUNMISC_ELLINIT, t_LFUNMISC_ETAQUO,
       t_LFUNMISC_GCHAR, t_LFUNMISC_ABELREL };
static long
lfundatatype(GEN data)
{
  switch(typ(data))
  {
    case t_INT: return t_LFUNMISC_CHIQUAD;
    case t_INTMOD: return t_LFUNMISC_CHICONREY;
    case t_POL: return t_LFUNMISC_POL;
    case t_VEC:
      switch(lg(data))
      {
        case 17: return t_LFUNMISC_ELLINIT;
        case 10: return t_LFUNMISC_POL;
        case 3:
          if (typ(gel(data,1)) != t_VEC) break;
          return is_gchar_group(gel(data,1))  ? t_LFUNMISC_GCHAR
                    : typ(gel(data,2))==t_MAT ? t_LFUNMISC_ABELREL
                                              : t_LFUNMISC_CHIGEN;
      }
      break;
  }
  return -1;
}
static GEN
lfunmisc_to_ldata_i(GEN ldata, long shallow)
{
  GEN x;
  if (is_linit(ldata)) ldata = linit_get_ldata(ldata);
  if (is_ldata(ldata) && is_tagged(ldata))
  {
    if (!shallow) ldata = gcopy(ldata);
    checkldata(ldata); return ldata;
  }
  x = checknf_i(ldata); if (x) return lfunzetak(x);
  switch (lfundatatype(ldata))
  {
    case t_LFUNMISC_POL: return lfunzetak(ldata);
    case t_LFUNMISC_CHIQUAD: return lfunchiquad(ldata);
    case t_LFUNMISC_CHICONREY:
    {
      GEN G = znstar0(gel(ldata,1), 1);
      return lfunchiZ(G, gel(ldata,2));
    }
    case t_LFUNMISC_CHIGEN:
    {
      GEN G = gel(ldata,1), chi = gel(ldata,2);
      switch(nftyp(G))
      {
        case typ_BIDZ: return lfunchiZ(G, chi);
        case typ_BNR: return lfunchigen(G, chi);
      }
    }
    break;
    case t_LFUNMISC_GCHAR: return lfungchar(gel(ldata,1), gel(ldata,2));
    case t_LFUNMISC_ABELREL:
      return lfunabelrel(gel(ldata,1), gel(ldata,2),
                         lg(ldata)==3? NULL: gel(ldata,3));
    case t_LFUNMISC_ELLINIT: return lfunell(ldata);
  }
  if (shallow != 2) pari_err_TYPE("lfunmisc_to_ldata",ldata);
  return NULL;
}

GEN
lfunmisc_to_ldata(GEN ldata)
{ return lfunmisc_to_ldata_i(ldata, 0); }

GEN
lfunmisc_to_ldata_shallow(GEN ldata)
{ return lfunmisc_to_ldata_i(ldata, 1); }

GEN
lfunmisc_to_ldata_shallow_i(GEN ldata)
{ return lfunmisc_to_ldata_i(ldata, 2); }

/********************************************************************/
/**                    High-level an expansion                     **/
/********************************************************************/
/* van is the output of ldata_get_an: return a_1,...a_L at precision prec */
GEN
ldata_vecan(GEN van, long L, long prec)
{
  GEN an = gel(van, 2);
  long t = mael(van,1,1);
  pari_timer ti;
  if (DEBUGLEVEL >= 1)
    err_printf("Lfun: computing %ld coeffs, prec %ld, type %ld\n", L, prec, t);
  if (DEBUGLEVEL >= 2) timer_start(&ti);
  switch (t)
  {
    long n;
    case t_LFUN_GENERIC:
      an = vecan_closure(an, L, prec);
      n = lg(an)-1;
      if (n < L)
      {
        pari_warn(warner, "#an = %ld < %ld, results may be imprecise", n, L);
        an = shallowconcat(an, zerovec(L-n));
      }
      break;
    case t_LFUN_CLOSURE0:
      pari_err_BUG("ldata_vecan: please call ldata_newprec");/*LCOV_EXCL_LINE*/
    case t_LFUN_ZETA: an = const_vecsmall(L, 1); break;
    case t_LFUN_NF:  an = dirzetak(an, stoi(L)); break;
    case t_LFUN_ELL:
      an = (ell_get_type(an) == t_ELL_Q) ? ellanQ_zv(an, L): ellan(an, L);
      break;
    case t_LFUN_KRONECKER: an = vecan_Kronecker(an, L); break;
    case t_LFUN_ABELREL: an = vecan_abelrel(an, L); break;
    case t_LFUN_CHIZ: an = vecan_chiZ(an, L, prec); break;
    case t_LFUN_CHIGEN: an = vecan_chigen(an, L, prec); break;
    case t_LFUN_HECKE: an = vecan_gchar(an, L, prec); break;
    case t_LFUN_ARTIN: an = vecan_artin(an, L, prec); break;
    case t_LFUN_ETA: an = vecan_eta(an, L); break;
    case t_LFUN_QF: an = vecan_qf(an, L); break;
    case t_LFUN_DIV: an = vecan_div(an, L, prec); break;
    case t_LFUN_MUL: an = vecan_mul(an, L, prec); break;
    case t_LFUN_CONJ: an = vecan_conj(an, L, prec); break;
    case t_LFUN_SYMPOW_ELL: an = vecan_ellsympow(an, L); break;
    case t_LFUN_GENUS2: an = vecan_genus2(an, L); break;
    case t_LFUN_HGM:
      an = hgmcoefs(gel(an,1), gel(an,2), L); break;
    case t_LFUN_MFCLOS:
    {
      GEN F = gel(an,1), E = gel(an,2), c = gel(an,3);
      an = mfcoefs(F,L,1) + 1; /* skip a_0 */
      an[0] = evaltyp(t_VEC)|evallg(L+1);
      an = mfvecembed(E, an);
      if (!isint1(c)) an = RgV_Rg_mul(an,c);
      break;
    }
    case t_LFUN_TWIST: an = vecan_twist(an, L, prec); break;
    case t_LFUN_SHIFT: an = vecan_shift(an, L, prec); break;
    default: pari_err_TYPE("ldata_vecan", van);
  }
  if (DEBUGLEVEL >= 2) timer_printf(&ti, "ldata_vecan");
  return an;
}

/* shallow */
GEN
ldata_newprec(GEN ldata, long prec)
{
  GEN van = ldata_get_an(ldata), an = gel(van, 2);
  long t = mael(van,1,1);
  switch (t)
  {
    case t_LFUN_CLOSURE0: return closure2ldata(an, prec);
    case t_LFUN_HECKE:
    {
      GEN gc = gel(an, 1), chiw = gel(an, 2);
      gc = gcharnewprec(gc, prec);
      return gchari_lfun(gc, chiw, gen_0); /* chi in internal format */
    }
    case t_LFUN_QF:
    {
      GEN eno = ldata_get_rootno(ldata);
      if (typ(eno)==t_REAL && realprec(eno) < prec) return lfunqf(an, prec);
      break;
    }
  }
  return ldata;
}

GEN
ldata_eulerf(GEN van, GEN p, long prec)
{
  GEN an = gel(van, 2), f = gen_0;
  long t = mael(van,1,1);
  switch (t)
  {
    case t_LFUN_GENERIC:
      f = eulerf_closure(an, p, prec); break;
    case t_LFUN_CLOSURE0:
      pari_err_BUG("ldata_vecan: please call ldata_newprec");/*LCOV_EXCL_LINE*/
    case t_LFUN_ZETA: f = mkrfrac(gen_1,deg1pol(gen_m1, gen_1,0)); break;
    case t_LFUN_NF:  f = eulerf_zetak(an, p); break;
    case t_LFUN_ELL: f = elleulerf(an, p); break;
    case t_LFUN_KRONECKER:
      f = mkrfrac(gen_1, deg1pol_shallow(stoi(-kronecker(an, p)), gen_1, 0)); break;
    case t_LFUN_ABELREL: f = eulerf_abelrel(an, p); break;
    case t_LFUN_CHIZ: f = eulerf_chiZ(an, p, prec); break;
    case t_LFUN_CHIGEN: f = eulerf_chigen(an, p, prec); break;
    case t_LFUN_HECKE: f = eulerf_gchar(an, p, prec); break;
    case t_LFUN_ARTIN: f = eulerf_artin(an, p, prec); break;
    case t_LFUN_DIV: f = eulerf_div(an, p, prec); break;
    case t_LFUN_MUL: f = eulerf_mul(an, p, prec); break;
    case t_LFUN_CONJ: f = eulerf_conj(an, p, prec); break;
    case t_LFUN_SYMPOW_ELL: f = eulerf_ellsympow(an, p); break;
    case t_LFUN_GENUS2: f = eulerf_genus2(an, p); break;
    case t_LFUN_TWIST: f = eulerf_twist(an, p, prec); break;
    case t_LFUN_SHIFT: f = eulerf_shift(an, p, prec); break;
    case t_LFUN_HGM: f = eulerf_hgm(an, p); break;
    default: f = NULL; break;
  }
  if (!f) pari_err_DOMAIN("lfuneuler", "L", "Euler product", strtoGENstr("unknown"), an);
  return f;
}

GEN
lfuneuler(GEN ldata, GEN p, long prec)
{
  pari_sp av = avma;
  if (typ(p)!=t_INT || signe(p)<=0) pari_err_TYPE("lfuneuler", p);
  ldata = ldata_newprec(lfunmisc_to_ldata_shallow(ldata), prec);
  return gerepilecopy(av, ldata_eulerf(ldata_get_an(ldata), p, prec));
}
