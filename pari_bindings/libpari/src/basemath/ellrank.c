/* Copyright (C) 2020  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_ellrank

static long LIM1 = 10000;
static long LIMTRIV = 10000;

/*******************************************************************/
/*               NFHILBERT and LOCAL SOLUBILITY                    */
/*       adapted from Denis Simon's original C implementation      */
/*******************************************************************/
/* p > 2, T ZX, p prime, x t_INT */
static long
lemma6(GEN T, GEN p, long nu, GEN x)
{
  long la, mu;
  GEN y = ZX_Z_eval(T, x);

  if (Zp_issquare(y, p)) return 1;
  la = Z_pval(y, p); y = ZX_Z_eval(ZX_deriv(T), x);
  if (!signe(y)) return la >= (nu<<1)? 0: -1;
  mu = Z_pval(y,p); if (la > (mu<<1)) return 1;
  return (la >= (nu<<1) && mu >= nu)? 0: -1;
}
static long
lemma7_aux(long nu, long la, long r)
{
  long nu2 = nu << 1;
  return (la >= nu2 || (la == nu2 - 2 && r == 1))? 0: -1;
}
/* p = 2, T ZX, x t_INT: return 1 = yes, -1 = no, 0 = inconclusive */
static long
lemma7(GEN T, long nu, GEN x)
{
  long r, la, mu;
  GEN y = ZX_Z_eval(T, x);

  if (!signe(y)) return 1;
  la = Z_lvalrem(y, 2, &y);
  r = Mod8(y); if (!odd(la) && r == 1) return 1;
  r &= 3; /* T(x) / 2^oo mod 4 */
  y = ZX_Z_eval(ZX_deriv(T), x);
  if (!signe(y)) return lemma7_aux(nu, la, r);
  mu = vali(y); if (la > mu<<1) return 1;
  if (nu <= mu) return lemma7_aux(nu, la, r);
  /* la <= 2mu, mu < nu */
  if (!odd(la) && mu + nu - la <= (r == 1? 2: 1)) return 1;
  return -1;
}

/* T a ZX, p a prime, pnu = p^nu, x0 t_INT */
static long
zpsol(GEN T, GEN p, long nu, GEN pnu, GEN x0)
{
  long i, res;
  pari_sp av = avma, btop;
  GEN x, pnup;

  res = absequaliu(p,2)? lemma7(T,nu,x0): lemma6(T,p,nu,x0);
  set_avma(av);
  if (res== 1) return 1;
  if (res==-1) return 0;
  x = x0; pnup = mulii(pnu,p);
  btop = avma;
  for (i=0; i < itos(p); i++)
  {
    x = addii(x,pnu);
    if (zpsol(T,p,nu+1,pnup,x)) return gc_long(av,1);
    if (gc_needed(btop, 2))
    {
      x = gerepileupto(btop, x);
      if (DEBUGMEM > 1)
        pari_warn(warnmem, "hyperell_locally_soluble: %ld/%Ps",i,p);
    }
  }
  return gc_long(av,0);
}

/* return 1 if equation y^2=T(x) has a rational p-adic solution (possibly
 * infinite), 0 otherwise. */
long
hyperell_locally_soluble(GEN T,GEN p)
{
  pari_sp av = avma;
  long res;
  if (typ(T)!=t_POL) pari_err_TYPE("hyperell_locally_soluble",T);
  if (typ(p)!=t_INT) pari_err_TYPE("hyperell_locally_soluble",p);
  RgX_check_ZX(T, "hyperell_locally_soluble");
  res = zpsol(T,p,0,gen_1,gen_0) || zpsol(RgX_recip_i(T), p, 1, p, gen_0);
  return gc_long(av, res);
}

/* is t a square in (O_K/pr) ? Assume v_pr(t) = 0 */
static long
quad_char(GEN nf, GEN t, GEN pr)
{
  GEN T, p, modpr = zk_to_Fq_init(nf, &pr,&T,&p);
  return Fq_issquare(nf_to_Fq(nf,t,modpr), T, p)? 1: -1;
}
/* quad_char(x), x in Z, nonzero mod p */
static long
Z_quad_char(GEN x, GEN pr)
{
  long f = pr_get_f(pr);
  if (!odd(f)) return 1;
  return kronecker(x, pr_get_p(pr));
}

/* (pr,2) = 1. return 1 if x in Z_K is a square in Z_{K_pr}, 0 otherwise.
 * modpr = zkmodprinit(nf,pr) */
static long
psquarenf(GEN nf,GEN x,GEN pr,GEN modpr)
{
  pari_sp av = avma;
  GEN p = pr_get_p(pr);
  long v;

  x = nf_to_scalar_or_basis(nf, x);
  if (typ(x) == t_INT) {
    if (!signe(x)) return 1;
    v = Z_pvalrem(x, p, &x) * pr_get_e(pr);
    if (v&1) return 0;
    v = (Z_quad_char(x, pr) == 1);
  } else {
    v = ZC_nfvalrem(x, pr, &x);
    if (v&1) return 0;
    v = (quad_char(nf, x, modpr) == 1);
  }
  return gc_long(av,v);
}

static long
ZV_iseven(GEN zlog)
{
  long i, l = lg(zlog);
  for (i = 1; i < l; i++)
    if (mpodd(gel(zlog,i))) return 0;
  return 1;
}

/* pr | 2, project to principal units (trivializes later discrete log) */
static GEN
to_principal_unit(GEN nf, GEN x, GEN pr, GEN sprk)
{
  if (pr_get_f(pr) != 1)
  {
    GEN prk = gel(sprk,3);
    x = nfpowmodideal(nf, x, gmael(sprk,5,1), prk);
  }
  return x;
}
/* pr | 2. Return 1 if x in Z_K is square in Z_{K_pr}, 0 otherwise */
static int
psquare2nf(GEN nf, GEN x, GEN pr, GEN sprk)
{
  long v = nfvalrem(nf, x, pr, &x);
  if (v == LONG_MAX) return 1; /* x = 0 */
  /* (x,pr) = 1 */
  if (odd(v)) return 0;
  x = to_principal_unit(nf, x, pr, sprk); /* = 1 mod pr */
  return ZV_iseven(sprk_log_prk1(nf, x, sprk));
}

/* pr above an odd prime */
static long
lemma6nf(GEN nf, GEN T, GEN pr, long nu, GEN x, GEN modpr)
{
  long la, mu;
  GEN y = nfpoleval(nf, T, x);

  if (psquarenf(nf,y,pr,modpr)) return 1;
  la = nfval(nf, y, pr); y = nfpoleval(nf, RgX_deriv(T), x);
  if (gequal0(y)) return la >= (nu<<1)? 0: -1;
  mu = nfval(nf, y, pr); if (la > (mu<<1)) return 1;
  return (la >= (nu<<1) && mu >= nu)? 0: -1;
}
/* pr above 2 */
static long
lemma7nf(GEN nf, GEN T, GEN pr, long nu, GEN x, GEN sprk)
{
  long i, res, la, mu, q, e, v;
  GEN M, y, gpx, loggx = NULL, gx = nfpoleval(nf, T, x);

  la = nfvalrem(nf, gx, pr, &gx); /* gx /= pi^la, pi a pr-uniformizer */
  if (la == LONG_MAX) return 1;
  if (!odd(la))
  {
    gx = to_principal_unit(nf, gx, pr, sprk); /* now 1 mod pr */
    loggx = sprk_log_prk1(nf, gx, sprk); /* cheap */
    if (ZV_iseven(loggx)) return 1;
  }
  gpx = nfpoleval(nf, RgX_deriv(T), x);
  mu = gequal0(gpx)? la+nu+1 /* oo */: nfval(nf,gpx,pr);

  if (la > (mu << 1)) return 1;
  if (nu > mu)
  {
    if (odd(la)) return -1;
    q = mu+nu-la; res = 1;
  }
  else
  {
    q = (nu << 1) - la;
    if (q <= 0) return 0;
    if (odd(la)) return -1;
    res = 0;
  }
  /* la even */
  e = pr_get_e(pr);
  if (q > e << 1)  return -1;
  if (q == 1) return res;

  /* gx = 1 mod pr; square mod pi^q ? */
  v = nfvalrem(nf, nfadd(nf, gx, gen_m1), pr, &y);
  if (v >= q) return res;
  /* 1 + pi^v y = (1 + pi^vz z)^2 mod pr^q ? v < q <= 2e => vz < e => vz = vy/2
   * => y = x^2 mod pr^(min(q-v, e+v/2)), (y,pr) = 1 */
  if (odd(v)) return -1;
  /* e > 1 */
  M = cgetg(2*e+1 - q + 1, t_VEC);
  for (i = q+1; i <= 2*e+1; i++) gel(M, i-q) = sprk_log_gen_pr(nf, sprk, i);
  M = ZM_hnfmodid(shallowconcat1(M), gen_2);
  return hnf_solve(M,loggx)? res: -1;
}
/* zinit either a sprk (pr | 2) or a modpr structure (pr | p odd).
   pnu = pi^nu, pi a uniformizer */
static long
zpsolnf(GEN nf,GEN T,GEN pr,long nu,GEN pnu,GEN x0,GEN repr,GEN zinit)
{
  long i, res;
  pari_sp av = avma;
  GEN pnup;

  res = typ(zinit) == t_VEC? lemma7nf(nf,T,pr,nu,x0,zinit)
                           : lemma6nf(nf,T,pr,nu,x0,zinit);
  set_avma(av);
  if (res== 1) return 1;
  if (res==-1) return 0;
  pnup = nfmul(nf, pnu, pr_get_gen(pr));
  nu++;
  for (i=1; i<lg(repr); i++)
  {
    GEN x = nfadd(nf, x0, nfmul(nf,pnu,gel(repr,i)));
    if (zpsolnf(nf,T,pr,nu,pnup,x,repr,zinit)) return gc_long(av,1);
  }
  return gc_long(av,0);
}

/* Let y = copy(x); y[k] := j; return y */
static GEN
ZC_add_coeff(GEN x, long k, long j)
{ GEN y = shallowcopy(x); gel(y, k) = utoipos(j); return y; }

/* system of representatives for Zk/pr */
static GEN
repres(GEN nf, GEN pr)
{
  long f = pr_get_f(pr), N = nf_get_degree(nf), p = itos(pr_get_p(pr));
  long i, j, k, pi, pf = upowuu(p, f);
  GEN perm = pr_basis_perm(nf, pr), rep = cgetg(pf+1,t_VEC);

  gel(rep,1) = zerocol(N);
  for (pi=i=1; i<=f; i++,pi*=p)
  {
    long t = perm[i];
    for (j=1; j<p; j++)
      for (k=1; k<=pi; k++) gel(rep, j*pi+k) = ZC_add_coeff(gel(rep,k), t, j);
  }
  return rep;
}

/* = 1 if equation y^2 = z^deg(T) * T(x/z) has a pr-adic rational solution
 * (possibly (1,y,0) = oo), 0 otherwise.
 * coeffs of T are algebraic integers in nf */
static long
locally_soluble(GEN nf,GEN T,GEN pr)
{
  GEN repr, zinit;

  if (typ(T)!=t_POL) pari_err_TYPE("nf_hyperell_locally_soluble",T);
  if (gequal0(T)) return 1;
  checkprid(pr);
  if (absequaliu(pr_get_p(pr), 2))
  { /* tough case */
    zinit = log_prk_init(nf, pr, 1+2*pr_get_e(pr), NULL);
    if (psquare2nf(nf,constant_coeff(T),pr,zinit)) return 1;
    if (psquare2nf(nf, leading_coeff(T),pr,zinit)) return 1;
  }
  else
  {
    zinit = zkmodprinit(nf, pr);
    if (psquarenf(nf,constant_coeff(T),pr,zinit)) return 1;
    if (psquarenf(nf, leading_coeff(T),pr,zinit)) return 1;
  }
  repr = repres(nf,pr); /* FIXME: inefficient if Npr is large */
  return zpsolnf(nf, T, pr, 0, gen_1, gen_0, repr, zinit) ||
         zpsolnf(nf, RgX_recip_i(T), pr, 1, pr_get_gen(pr),
                 gen_0, repr, zinit);
}
long
nf_hyperell_locally_soluble(GEN nf,GEN T,GEN pr)
{
  pari_sp av = avma;
  return gc_long(av, locally_soluble(nf, T, pr));
}

/* return a * denom(a)^2, as an 'liftalg' */
static GEN
den_remove(GEN nf, GEN a)
{
  GEN da;
  a = nf_to_scalar_or_basis(nf, a);
  switch(typ(a))
  {
    case t_INT: return a;
    case t_FRAC: return mulii(gel(a,1), gel(a,2));
    case t_COL:
      a = Q_remove_denom(a, &da);
      if (da) a = ZC_Z_mul(a, da);
      return nf_to_scalar_or_alg(nf, a);
    default: pari_err_TYPE("nfhilbert",a);
      return NULL;/*LCOV_EXCL_LINE*/
  }
}

static long
hilb2nf(GEN nf,GEN a,GEN b,GEN p)
{
  pari_sp av = avma;
  GEN pol;
  a = den_remove(nf, a);
  b = den_remove(nf, b);
  pol = mkpoln(3, a, gen_0, b);
  /* varn(nf.pol) = 0, pol is not a valid GEN  [as in Pol([x,x], x)].
   * But it is only used as a placeholder, hence it is not a problem */
  return gc_long(av, nf_hyperell_locally_soluble(nf,pol,p)? 1: -1);
}

/* local quadratic Hilbert symbol (a,b)_pr, for a,b (nonzero) in nf */
static long
nfhilbertp(GEN nf, GEN a, GEN b, GEN pr)
{
  GEN t;
  long va, vb;
  pari_sp av = avma;

  if (absequaliu(pr_get_p(pr), 2)) return hilb2nf(nf,a,b,pr);

  /* pr not above 2, compute t = tame symbol */
  va = nfval(nf,a,pr);
  vb = nfval(nf,b,pr);
  if (!odd(va) && !odd(vb)) return gc_long(av,1);
  /* Trick: pretend the exponent is 2, result is OK up to squares ! */
  t = famat_makecoprime(nf, mkvec2(a,b), mkvec2s(vb, -va),
                        pr, pr_hnf(nf, pr), gen_2);
  /* quad. symbol is image of t = (-1)^(v(a)v(b)) a^v(b) b^(-v(a))
   * by the quadratic character  */
  switch(typ(t))
  {
    default: /* t_COL */
      if (!ZV_isscalar(t)) break;
      t = gel(t,1); /* fall through */
    case t_INT:
      if (odd(va) && odd(vb)) t = negi(t);
      return gc_long(av,  Z_quad_char(t, pr));
  }
  if (odd(va) && odd(vb)) t = ZC_neg(t);
  return gc_long(av, quad_char(nf, t, pr));
}

/* Global quadratic Hilbert symbol (a,b):
 *  =  1 if X^2 - aY^2 - bZ^2 has a point in projective plane
 *  = -1 otherwise
 * a, b should be nonzero */
long
nfhilbert(GEN nf, GEN a, GEN b)
{
  pari_sp av = avma;
  long i, l;
  GEN S, S2, Sa, Sb, sa, sb;

  a = nf_to_scalar_or_basis(nf, a);
  b = nf_to_scalar_or_basis(nf, b);
  /* local solutions in real completions ? [ error in nfsign if arg is 0 ]*/
  sa = nfsign(nf, a);
  sb = nfsign(nf, b); l = lg(sa);
  for (i=1; i<l; i++)
    if (sa[i] && sb[i])
    {
      if (DEBUGLEVEL>3)
        err_printf("nfhilbert not soluble at real place %ld\n",i);
      return gc_long(av,-1);
    }

  /* local solutions in finite completions ? (pr | 2ab)
   * primes above 2 are toughest. Try the others first */
  Sa = idealfactor(nf, a);
  Sb = idealfactor(nf, b);
  S2 = idealfactor(nf, gen_2);
  S = merge_factor(Sa, Sb, (void*)&cmp_prime_ideal, &cmp_nodata);
  S = merge_factor(S,  S2, (void*)&cmp_prime_ideal, &cmp_nodata);
  S = gel(S,1);
  /* product of all hilbertp is 1 ==> remove one prime (above 2!) */
  for (i=lg(S)-1; i>1; i--)
    if (nfhilbertp(nf,a,b,gel(S,i)) < 0)
    {
      if (DEBUGLEVEL>3)
        err_printf("nfhilbert not soluble at finite place %Ps\n",S[i]);
      return gc_long(av,-1);
    }
  return gc_long(av,1);
}

long
nfhilbert0(GEN nf,GEN a,GEN b,GEN p)
{
  nf = checknf(nf);
  if (p) {
    checkprid(p);
    if (gequal0(a)) pari_err_DOMAIN("nfhilbert", "a", "=", gen_0, a);
    if (gequal0(b)) pari_err_DOMAIN("nfhilbert", "b", "=", gen_0, b);
    return nfhilbertp(nf,a,b,p);
  }
  return nfhilbert(nf,a,b);
}

/*******************************************************************/
/*                      HYPERELL_LOCAL_SOLVE                       */
/*******************************************************************/

/* Based on
T.A. Fisher and G.F. Sills
Local solubility and height bounds for coverings of elliptic curves
https://www.dpmms.cam.ac.uk/~taf1000/papers/htbounds.pdf
*/

static int
FpX_issquare(GEN q, GEN p)
{
  GEN F = FpX_factor_squarefree(q,p);
  long i, l = lg(F);
  for (i = 1; i < l; i+=2)
    if (degpol(gel(F,i)) > 0) return 0;
  return 1;
}

static GEN
hyperell_red(GEN q, GEN p)
{
  GEN Q;
  long v = ZX_pvalrem(q, p, &Q);
  if (v == 1) return q;
  return odd(v)? ZX_Z_mul(Q, p): Q;
}

static GEN
hyperell_reg_point(GEN q, GEN p)
{
  GEN Q, F;
  long i, l, v = ZX_pvalrem(q, p, &Q);
  if (v != 1) q = odd(v)? ZX_Z_mul(Q, p): Q;
  if (!odd(v))
  {
    GEN qr = FpX_red(q, p);
    if (!FpX_issquare(qr,p) || Fp_issquare(leading_coeff(qr), p))
      return mkvec2s(0,1);
  }
  F = FpX_roots(Q, p); l = lg(F);
  for (i = 1; i < l; i++)
  {
    GEN r = gel(F,i), s = hyperell_reg_point(ZX_affine(q, p, r), p);
    if (s)
      retmkvec2(addii(r,mulii(gel(s,1),p)), mulii(gel(s,2),p));
  }
  return NULL;
}

static GEN
hyperell_lift(GEN q, GEN x, GEN p)
{
  long e;
  GEN qy2 = ZX_Z_sub(q, sqri(p));
  for (e = 2; ; e<<=1)
  {
    pari_sp av = avma;
    GEN z = ZpX_liftroot(qy2, x, p, e);
    if (signe(z) == 0) z = powiu(p, e);
    if (Zp_issquare(poleval(q, z), p)) return z;
    set_avma(av);
  }
}

static GEN
affine_apply(GEN r, GEN x)
{
  return addii(mulii(gel(r,2),x), gel(r,1));
}

static GEN
Qp_hyperell_solve_odd(GEN q, GEN p)
{
  GEN qi = RgX_recip_shallow(q);
  GEN r = hyperell_reg_point(q,  p), qr = NULL, qrp = NULL;
  GEN s = hyperell_reg_point(qi, p), qs = NULL, qsp = NULL;
  if (!r && !s) return NULL;
  if (r)
  {
    qr = hyperell_red(ZX_affine(q,  gel(r,2), gel(r,1)), p);
    qrp = FpX_deriv(qr, p);
  }
  if (s)
  {
    qs = hyperell_red(ZX_affine(qi, gel(s,2), gel(s,1)), p);
    qsp = FpX_deriv(qs, p);
  }
  while(1)
  {
    GEN x = randomi(p);
    if (r)
    {
      GEN y2 = FpX_eval(qr, x, p);
      if (Fp_issquare(y2,p))
      {
         if (signe(y2))
           return affine_apply(r,x);
         if (signe(FpX_eval(qrp, x, p)))
         {
           x = hyperell_lift(qr, x, p);
           return affine_apply(r,x);
         }
      }
    }
    if (s)
    {
      GEN y2 = FpX_eval(qs, x, p);
      if (Fp_issquare(y2,p))
      {
         if (signe(x)==0) x = p;
         if (signe(y2))
           return ginv(affine_apply(s,x));
         if (signe(FpX_eval(qsp, x, p)))
         {
           GEN xl = hyperell_lift(qs, x, p);
           return ginv(affine_apply(s,xl));
         }
      }
    }
  }
}

static GEN
Q2_hyperell_lift(GEN p, GEN q, long x, long y)
{
  GEN T, h;
  long e;
  if (y==0) y = 2;
  T = ZX_sub(p, ZX_Z_add(ZX_mulu(q, y), sqru(y)));
  h = ZX_add(ZX_sqr(q), ZX_shifti(p, 2));
  for (e = 3; ; e++)
  {
    pari_sp av = avma;
    GEN r = ZpX_liftroot(T, utoi(x), gen_2, e);
    if (signe(r) == 0) r = int2n(e);
    if (Zp_issquare(poleval(h, r), gen_2)) return r;
    set_avma(av);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN
Q2_hyperell_regpoint(GEN P, GEN Q)
{
  long x;
  GEN p = ZX_to_F2x(P), dp = F2x_deriv(p);
  GEN q = ZX_to_F2x(Q), dq = F2x_deriv(q);

  for (x = 0; x <= 1; x++)
  {
    long px = F2x_eval(p,x), qx = F2x_eval(q,x);
    long dpx, dqx;
    if (qx == 1)
    {
      if(px == 0) return x==0 ? gen_2: gen_1;
      continue;
    }
    dpx = F2x_eval(dp,x);
    dqx = F2x_eval(dq,x);
    if (odd(dqx * px + dpx))
      return Q2_hyperell_lift(P, Q, x, px);
  }
  return NULL;
}

static GEN
Q2_hyperell_solve_affine(GEN p, GEN q)
{
  pari_sp av = avma;
  GEN R, p4, q4;
  long x;
  while(1)
  {
    GEN pp, p0, p1;
    long vp = ZX_lval(p, 2);
    long vq = ZX_lval(q, 2);
    long w = minss(vp>>1, vq);
    p = ZX_shifti(p, -2*w);
    q = ZX_shifti(q, -w);
    if (ZX_lval(q,2)==0) break;
    pp = RgX_splitting(p,2); p0 = gel(pp,1); p1 = gel(pp,2);
    if (ZX_lval(p1,2)==0 || ZX_lval(p0,2)>=1) break;
    p = ZX_sub(p, ZX_mul(p0, ZX_add(q, p0)));
    q = ZX_add(q, ZX_shifti(p0, 1));
  }
  R = Q2_hyperell_regpoint(p, q);
  if (R) return gerepileuptoint(av, R);
  p4 = ZX_to_Flx(p,4);
  q4 = ZX_to_Flx(q,4);
  for (x = 0; x <= 1; x++)
  {
    ulong px = Flx_eval(p4, x, 4);
    ulong qx = Flx_eval(q4, x, 4);
    if (px == 0 || (1+qx+3*px)%4==0)
    {
      GEN p2 = ZX_affine(p, gen_2, utoi(x));
      GEN q2 = ZX_affine(q, gen_2, utoi(x));
      GEN S = Q2_hyperell_solve_affine(p2, q2);
      if (S) return gerepileuptoint(av, addiu(shifti(S,1),x));
    }
  }
  return gc_NULL(av);
}

static GEN
Q2_hyperell_solve(GEN P)
{
  long v = varn(P);
  GEN S = Q2_hyperell_solve_affine(P, pol_0(v));
  if (!S) S = ginv(Q2_hyperell_solve_affine(RgX_recip(P), pol_0(v)));
  return S;
}

static GEN
hyperell_local_solve(GEN q, GEN p)
{
  if (equaliu(p,2))
    return Q2_hyperell_solve(q);
  return Qp_hyperell_solve_odd(q, p);
}

/*******************************************************************/
/*                         BINARY QUARTIC                          */
/*******************************************************************/
static int
Qp_issquare(GEN a, GEN p)
{
  GEN b = typ(a) == t_INT? a: mulii(gel(a,1), gel(a,2));
  return Zp_issquare(b, p);
}

static GEN
quartic_IJ(GEN g)
{
  GEN a = gel(g, 6), b = gel(g, 5), c = gel(g, 4), d = gel(g, 3), e = gel(g, 2);
  GEN ae = gmul(a,e), bd = gmul(b,d), c2 = gsqr(c);
  /* 12ae - 3bd + c^2 */
  GEN I = gadd(gsub(gmulsg(12, ae), gmulsg(3, bd)), c2);
  /* c(72ae + 9bd - 2c^2) - 27ad^2 - 27eb^2*/
  GEN J = gsub(gmul(c, gsub(gadd(gmulsg(72,ae), gmulsg(9,bd)), gmul2n(c2,1))),
               gmulsg(27, gadd(gmul(a, gsqr(d)), gmul(gsqr(b), e))));
  return mkvec2(I, J);
}

static GEN
quartic_hessiandd(GEN g)
{
  GEN a = gel(g, 6), b = gel(g, 5), c = gel(g, 4), d = gel(g, 3), e = gel(g, 2);
  GEN a8 = gmul2n(a, 3), p4 = gsub(gmulsg(3, gsqr(b)), gmul(a8, c));
  GEN p3 = gsub(gmul(b, c), gmul(gmulsg(6, a), d));
  GEN p2 = gsub(gmulsg(8, gsqr(c)), gmulsg(12, gadd(gmul(b, d), gmul(a8, e))));
  return deg2pol_shallow(gmulgu(p4,12), gmulgu(p3,24), p2, 1);
}

static GEN
quartic_cubic(GEN g, long v)
{
  GEN a = gel(g, 6), b = gel(g, 5), c = gel(g, 4);
  GEN a3 = gdivgu(a,3);
  return deg1pol(gmul2n(a3,2), gsub(gsqr(b),gmul2n(gmul(a3,c),3)), v);
}

static GEN
quarticinv_pol(GEN IJ)
{
  GEN I = gel(IJ,1), J = gel(IJ,2);
  return mkpoln(4, gen_1, gen_0, gmulgs(I,-3), J);
}
static GEN
quartic_H(GEN g, GEN *pT)
{
  GEN a = gel(g, 6), b = gel(g, 5), c = gel(g, 4);
  GEN IJ = quartic_IJ(g), I = gel(IJ, 1);
  GEN ddh = quartic_hessiandd(g);
  GEN ddg = deg2pol_shallow(gmulgu(a,12), gmulgu(b,6), gmulgu(c,2), 1);
  *pT = quarticinv_pol(IJ);
  return deg2pol_shallow(stoi(-8), gmul2n(ddg,2), gadd(ddh,gmul2n(I,3)),0);
}

static GEN
quartic_disc(GEN q)
{
  GEN IJ = quartic_IJ(q), I = gel(IJ,1), J = gel(IJ,2);
  return gsub(gmul2n(gpowgs(I, 3), 2), gsqr(J));
}

static GEN
quartic_minim_all(GEN F, GEN discF)
{
  GEN IJ = quartic_IJ(F), I = gel(IJ,1), J = gel(IJ,2);
  GEN g = Z_ppo(gcdii(I,J), gel(discF,1));
  GEN plist = ZV_sort_uniq_shallow(shallowconcat(gel(absZ_factor(g),1),gel(discF,2)));
  GEN W, C, PQ = hyperellminimalmodel(F, &C, plist);
  GEN P = gel(PQ,1), Q = gel(PQ,2);
  if (signe(Q)==0)
    W = mkvec2(P, C);
  else
    W = mkvec2(ZX_add(ZX_shifti(P,2),ZX_sqr(Q)),mkvec2(shifti(gel(C,1),-1),gel(C,2)));
  return W;
}

/*******************************************************************/
/*                         Cassels' pairing                        */
/*******************************************************************/

static GEN
nfsqrt(GEN nf, GEN z)
{
  long v = fetch_var_higher();
  GEN R = nfroots(nf, deg2pol_shallow(gen_m1, gen_0, z, v));
  delete_var();
  return lg(R)==1 ? NULL: gel(R,1);
}

static GEN
nfsqrt_safe(GEN nf, GEN z)
{
  GEN r = nfsqrt(nf, z);
  if (!r) pari_err_BUG("ellrank");
  return r;
}

static GEN
vecnfsqrtmod(GEN x, GEN P)
{ pari_APPLY_same(gmodulo(nfsqrt_safe(gel(x,i), P), gel(x,i))) }

static GEN
enfsqrt(GEN T, GEN P)
{
  GEN F = gel(ZX_factor(T),1);
  return liftpol(chinese1(vecnfsqrtmod(F,P)));
}

/* Quartic q, at most quadratic g s.t. lc(g) > 0. There exist a real r s.t.
 * q(r) > 0 and g(r) != 0. Return sign(g(r)) */
static int
cassels_oo_solve_i(GEN q, GEN g)
{
  long dg = degpol(g);
  GEN D, a, b, c;

  if (dg == 0 || signe(leading_coeff(q)) > 0) return 1;
  if (signe(gel(q,2)) > 0) return signe(gel(g,2));
  c = gel(g,2); b = gel(g,3);
  /* g = bx+c, b>0, is negative on I=]-oo,-c/b[: if q has a root there,
   * then g(r) < 0. Else it has the sign of q(oo) < 0 on I */
  if (dg == 1) return ZX_sturmpart(q, mkvec2(mkmoo(), gdiv(negi(c), b)))? -1: 1;
  a = gel(g,4); D = subii(sqri(b), shifti(mulii(a,c), 2)); /* g = ax^2+bx+c */
  if (signe(D) <= 0) return 1; /* sign(g) = 1 is constant */
  /* Rescale q and g: x->(x - b)/2a; roots of new g are \pm sqrt(D) */
  q = ZX_translate(ZX_rescale(q, shifti(a,1)), negi(b));
  /* Now g has sign -1 in I=[-sqrt(D),sqrt(D)] and 1 elsewhere.
   * Check if q vanishes in I <=> Graeffe(q) vanishes on [0,D].
   * If so or if q(0) > 0 we take r in there; else r is outside of I */
  return (signe(gel(q,2)) > 0 ||
          ZX_sturmpart(ZX_graeffe(q), mkvec2(gen_0, D)))? -1: 1;
}
static int
cassels_oo_solve(GEN q, GEN g)
{ pari_sp av = avma; return gc_int(av, cassels_oo_solve_i(q, g)); }

static GEN
cassels_Qp_solve(GEN q, GEN gam, GEN p)
{
  pari_sp av = avma;
  GEN a = hyperell_local_solve(q, p);
  GEN c = poleval(gam,a);
  long e;
  if (!gequal0(c)) return c;
  for (e = 2; ; e++)
  {
    GEN b = gadd(a, powiu(p,e));
    if (Qp_issquare(poleval(q, b), p))
    {
      c = poleval(gam, b);
      if (!gequal0(c)) return gerepileupto(av,c);
    }
  }
}

static GEN
to_ZX(GEN a, long v) { return typ(a)==t_INT? scalarpol_shallow(a,v): a; }

static GEN
quartic_findunit(GEN D, GEN q)
{
  GEN T = quarticinv_pol(quartic_IJ(q));
  while(1)
  {
    pari_sp av = avma;
    GEN z = quartic_cubic(q,0);
    if (signe(QXQ_norm(z,T)))
      return absequalii(quartic_disc(q), D)? q: ZX_shifti(q, 2);
    set_avma(av);
    q = ZX_translate(RgX_recip(q), gen_1);
  }
}

/* Crude implementation of an algorithm by Tom Fisher
 * On binary quartics and the Cassels-Tate pairing
 * https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf */

/* FD = primes | 2*3*5*7*D, q1,q2,q3 have discriminant D */
static long
casselspairing(GEN FD, GEN q1, GEN q2, GEN q3)
{
  pari_sp av = avma;
  GEN T, H = quartic_H(q1, &T);
  GEN z1 = quartic_cubic(q1, 0);
  GEN z2 = quartic_cubic(q2, 0);
  GEN z3 = quartic_cubic(q3, 0);
  GEN m = to_ZX(enfsqrt(T, QXQ_mul(QX_mul(z1,z2),z3,T)), 0);
  GEN Hm = RgXQ_mul(QXQ_div(m, z1, T), H, T); /* deg(Hm) >= 2 */
  GEN gam = to_ZX(Q_primpart(gel(Hm,4)),1);
  GEN a = leading_coeff(q2), Fa = gel(absZ_factor(a),1);
  GEN F = ZV_sort_uniq_shallow(shallowconcat1(mkvec2(Fa, FD)));
  long i, e = 0, lF = lg(F);
  if (signe(a) <= 0)
  {
    if (signe(leading_coeff(gam)) < 0) gam = ZX_neg(gam);
    if (cassels_oo_solve(q1, gam) < 0) e = 1;
  }
  for (i = 1; i < lF; i++)
  {
    GEN p = gel(F, i);
    GEN c = cassels_Qp_solve(q1, gam, p);
    if (hilbert(c, a, p) < 0) e = !e;
  }
  return gc_long(av,e);
}

static GEN
matcassels(GEN F, GEN M)
{
  long i, j, n = lg(M)-1;
  GEN C = zero_F2m_copy(n,n);
  pari_sp av = avma;
  for (i = 1; i <= n; i++)
  {
    GEN Mii = gcoeff(M,i,i);
    if (isintzero(Mii)) continue;
    for (j = 1; j < i; j++)
    {
      GEN Mjj = gcoeff(M,j,j);
      if (!isintzero(Mjj) && casselspairing(F, Mii, Mjj, gcoeff(M,i,j)))
      { F2m_set(C,i,j); F2m_set(C,j,i); }
    }
  }
  return gc_const(av, C);
}

/*******************************************************************/
/*                         ELLRANK                                 */
/*******************************************************************/
/* This section is a port by Bill Allombert of ellQ.gp by Denis Simon
 * Copyright (C) 2019 Denis Simon
 * Distributed under the terms of the GNU General Public License (GPL) */

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
/*    FUNCTIONS FOR POLYNOMIALS                \\ */
/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */

static GEN
ell2pol(GEN ell)
{ return mkpoln(4, gen_1, ell_get_a2(ell), ell_get_a4(ell), ell_get_a6(ell)); }

/* find point (x:y:z) on y^2 = pol, return [x,z]~ and set *py = y */
static GEN
projratpointxz(GEN pol, long lim, GEN *py)
{
  pari_timer ti;
  GEN P;
  if (issquareall(leading_coeff(pol), py)) return mkcol2(gen_1, gen_0);
  if (DEBUGLEVEL) timer_start(&ti);
  P = hyperellratpoints(pol, stoi(lim), 1);
  if (DEBUGLEVEL) timer_printf(&ti,"hyperellratpoints(%ld)",lim);
  if (lg(P) == 1) return NULL;
  P = gel(P,1); *py = gel(P,2); return mkcol2(gel(P,1), gen_1);
}

/* P a list of integers (actually primes) one of which divides x; return
 * the first one */
static GEN
first_divisor(GEN x, GEN P)
{
  long i, n = lg(P);
  for (i = 1; i < n; i++)
    if (dvdii(x, gel(P,i))) return gel(P,i);
  return gel(P,i);
}

/* find point (x:y:z) on y^2 = pol, return [x,z]~ and set *py = y */
static GEN
projratpointxz2(GEN pol, long lim, GEN *py)
{
  pari_sp av = avma;
  GEN list = mkvec(mkvec4(pol, matid(2), gen_1, gen_1));
  long i, j, c;

  for (i = 1, c = 1; i < lg(list); i++,c++)
  {
    GEN K, k, ff, co, p, M, C, r, pol, L = gel(list, i);
    long lr;

    list = vecsplice(list, i); i--;
    pol = Q_primitive_part(gel(L,1), &K);
    M = gel(L,2);
    K = K? mulii(gel(L,3), K): gel(L,3);
    C = gel(L,4);
    if (Z_issquareall(K, &k))
    {
      GEN xz, y, aux, U;
      if (c==1) continue;
      pol = ZX_hyperellred(pol, &U);
      if (DEBUGLEVEL>1) err_printf("  reduced quartic(%ld): Y^2 = %Ps\n", i, pol);
      xz = projratpointxz(pol, lim, &y); if (!xz) continue;
      *py = gmul(y, mulii(C, k));
      aux = RgM_RgC_mul(ZM2_mul(M, U), xz);
      if (gequal0(gel(aux, 2))) return mkcol2(gel(aux,1), gen_0);
      *py = gdiv(*py, gpowgs(gel(aux, 2), degpol(pol)>>1));
      return mkcol2(gdiv(gel(aux, 1), gel(aux, 2)), gen_1);
    }
    ff = Z_factor(K); co = core2(mkvec2(K, ff)); K = gel(co,1); /* > 1 */
    p = first_divisor(K, gel(ff,1));
    K = diviiexact(K, p);
    C = mulii(mulii(C, gel(co,2)), p);
    /* root at infinity */
    if (dvdii(leading_coeff(pol), p))
    {
      GEN U = mkmat2(gel(M,1), ZC_Z_mul(gel(M,2), p));
      if (equali1(content(U)))
      {
        GEN t = ZX_rescale(pol, p);
        list = vec_append(list, mkvec4(ZX_Z_divexact(t,p), U, K, C));
      }
    }
    r = FpC_center(FpX_roots(pol, p), p, shifti(p,-1)); lr = lg(r);
    for (j = 1; j < lr; j++)
    {
      GEN U = ZM2_mul(M, mkmat22(p, gel(r, j), gen_0, gen_1));
      if (equali1(content(U)))
      {
        GEN t = ZX_unscale_div(ZX_translate(pol, gel(r,j)), p);
        list = vec_append(list, mkvec4(t, U, K, C));
      }
    }
    if (gc_needed(av, 1)) gerepileall(av, 2, &pol, &list);
  }
  return NULL;
}

static GEN
polrootsmodpn(GEN pol, GEN p)
{
  pari_sp av = avma, av2;
  long j, l, i = 1, vd = Z_pval(ZX_disc(pol), p);
  GEN v, r, P;

  if (!vd) { set_avma(av); retmkvec(zerovec(2)); }
  pol = Q_primpart(pol);
  P = gpowers0(p, vd-1, p); av2 = avma;
  v = FpX_roots(pol, p); l = lg(v);
  for (j = 1; j < l; j++) gel(v,j) = mkvec2(gel(v,j), gen_1);
  while (i < lg(v))
  {
    GEN pol2, pe, roe = gel(v, i), ro = gel(roe,1);
    long e = itou(gel(roe,2));

    if (e >= vd) { i++; continue; }
    pe = gel(P, e);
    (void)ZX_pvalrem(ZX_affine(pol, pe, ro), p, &pol2);
    r = FpX_roots(pol2, p); l = lg(r);
    if (l == 1) { i++; continue; }
    for (j = 1; j < l; j++)
      gel(r, j) = mkvec2(addii(ro, mulii(pe, gel(r, j))), utoi(e + 1));
    /* roots with higher precision = ro + r*p^(e+1) */
    if (l > 2) v = shallowconcat(v, vecslice(r, 2, l-1));
    gel(v, i) = gel(r, 1);
    if (gc_needed(av2, 1)) gerepileall(av2, 1, &v);
  }
  if (lg(v) == 1) { set_avma(av); retmkvec(zerovec(2)); }
  return gerepilecopy(av, v);
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
/*    FUNCTIONS FOR LOCAL COMPUTATIONS         \\ */
/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */

/* z is integral; sprk true (pr | 2) [t_VEC] or modpr structure (pr | p odd)
 * [t_COL] */
static GEN
kpmodsquares1(GEN nf, GEN z, GEN sprk)
{
  GEN modpr = (typ(sprk) == t_VEC)? gmael(sprk, 4, 1): sprk;
  GEN pr = modpr_get_pr(modpr), p = pr_get_p(pr);
  long v = nfvalrem(nf, z, pr, &z);
  if (equaliu(p, 2))
  {
    GEN c;
    z = to_principal_unit(nf, z, pr, sprk); /* = 1 mod pr */
    c = ZV_to_Flv(sprk_log_prk1(nf, z, sprk), 2);
    return vecsmall_prepend(c, odd(v));
  }
  else
  {
    GEN T = modpr_get_T(modpr);
    long c = !Fq_issquare(nf_to_Fq(nf, z, modpr), T, p);
    return mkvecsmall2(odd(v), c);
  }
}

static GEN
kpmodsquares(GEN vnf, GEN z, GEN PP)
{
  pari_sp av = avma;
  long i, j, l = lg(vnf);
  GEN dz, vchar = cgetg(l, t_VEC);
  z = Q_remove_denom(z, &dz); if (dz) z = ZX_Z_mul(z, dz);
  for (i = 1; i < l; i++)
  {
    GEN nf = gel(vnf, i), pp = gel(PP, i);
    GEN kp, delta = ZX_rem(z, nf_get_pol(nf));
    long lp = lg(pp);
    kp = cgetg(lp, t_VEC);
    for (j = 1; j < lp; j++) gel(kp, j) = kpmodsquares1(nf, delta, gel(pp,j));
    gel(vchar, i) = shallowconcat1(kp);
  }
  return gerepileuptoleaf(av, shallowconcat1(vchar));
}

static GEN
veckpmodsquares(GEN x, GEN vnf, GEN PP)
{ pari_APPLY_type(t_MAT, kpmodsquares(vnf, gel(x, i), PP)) }

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
/*    GENERIC FUNCTIONS FOR ELLIPTIC CURVES    \\ */
/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */

static GEN
ellabs(GEN P)
{ return ell_is_inf(P) ? P: mkvec2(gel(P,1), Q_abs_shallow(gel(P,2))); }
static GEN
vecellabs(GEN x) { pari_APPLY_same(ellabs(gel(x,i))) }

/* y^2 = x^3 + K a2 x + K^2 a4 x + K^3 a6 */
static GEN
elltwistequation(GEN ell, GEN K)
{
  GEN K2, a2, a4, a6;
  if (!K || equali1(K)) return ell;
  K2 = sqri(K);
  a2 = mulii(ell_get_a2(ell), K);
  a4 = mulii(ell_get_a4(ell), K2);
  a6 = mulii(ell_get_a6(ell), mulii(K, K2));
  return ellinit(mkvec5(gen_0, a2, gen_0, a4, a6), NULL, DEFAULTPREC);
}

/* P=[x,y] a point on Ky^2 =  pol(x); [Kx, K^2y] point on Y^2 = K^3pol(X/K) */
static GEN
elltwistpoint(GEN P, GEN K, GEN K2)
{
  if (ell_is_inf(P)) return ellinf();
  return mkvec2(gmul(gel(P,1), K), gmul(gel(P,2), K2));
}

static GEN
elltwistpoints(GEN x, GEN K)
{
  GEN K2;
  if (!K || gequal1(K)) return x;
  K2 = gsqr(K);
  pari_APPLY_same(elltwistpoint(gel(x,i), K, K2))
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
/*    FUNCTIONS FOR NUMBER FIELDS              \\ */
/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */

/* return a set S2 of prime ideals disjoint from S such that
 * Cl_{S \cup S2}(K) has no p-torsion */
static GEN
bestS(GEN bnf,GEN S, ulong p)
{
  GEN v, S2, h = bnf_get_no(bnf), cyc = bnf_get_cyc(bnf);
  long i, lS2;
  ulong l, vD;
  forprime_t P;

  if (!dvdiu(h, p)) return cgetg(1,t_VEC);
  if (!S)
  {
    v = diagonal_shallow(cyc);
    vD = Z_lval(h, p);
  }
  else
  {
    long lS = lg(S);
    v = cgetg(lS,t_MAT);
    for (i = 1; i < lS; i++) gel(v,i) = isprincipal(bnf, gel(S,i));
    v = ZM_hnfmodid(v, cyc);
    vD = Z_lval(ZM_det(v), p); if (!vD) return cgetg(1, t_VEC);
  }
  S2 = cgetg(vD+2, t_VEC); lS2 = 1;
  u_forprime_init(&P,2,ULONG_MAX);
  while ((l = u_forprime_next(&P)))
  {
    pari_sp av = avma;
    GEN w, Sl = idealprimedec(bnf, utoi(l));
    long nSl = lg(Sl)-1;
    ulong vDl;
    for (i = 1; i < nSl; i++) /* remove one prime ideal */
    {
      w = ZM_hnf(shallowconcat(v, isprincipal(bnf, gel(Sl,i))));
      vDl = Z_lval(ZM_det(w), p);
      if (vDl < vD)
      {
        gel(S2,lS2++) = gel(Sl,i);
        vD = vDl; v = w; av = avma;
        if (!vD) { setlg(S2, lS2); return S2; }
      }
    }
    set_avma(av);
  }
  return NULL;/*LCOV_EXCL_LINE*/
}

static GEN
nfC_prV_val(GEN nf, GEN G, GEN P)
{
  long i, j, lG = lg(G), lP = lg(P);
  GEN M = cgetg(lG, t_MAT);
  for (i = 1; i < lG; i++)
  {
    GEN V = cgetg(lP, t_COL);
    for (j = 1; j < lP; j++)
      gel(V,j) = gpnfvalrem(nf, gel(G,i), gel(P,j), NULL);
    gel(M,i) = V;
  }
  return M;
}

static GEN
_factorbackmod(GEN nf, GEN g, GEN e, ulong p)
{
  GEN y = nffactorback(nf, g, e), den;
  GEN z = nfmul(nf, y, nfsqr(nf, idealredmodpower(nf, y, p, 0)));
  z = Q_remove_denom(z, &den);
  if (den)
  {
    if (p != 2) den = powiu(den, p-1);
    z = gmul(z, den);
  }
  return z;
}
static GEN
nfV_factorbackmod(GEN nf, GEN x, ulong p)
{
  long i, l = lg(x);
  GEN v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN y = gel(x,i), g = gel(y,1), e = gel(y,2);
    gel(v,i) = _factorbackmod(nf, g, ZV_to_Flv(e,p), p);
  }
  return v;
}
static GEN
nfV_zm_factorback(GEN nf, GEN G, GEN x, long p)
{ pari_APPLY_type(t_VEC, _factorbackmod(nf, G, gel(x,i), p)) }

static GEN
prV_ZM_factorback(GEN nf, GEN S, GEN x)
{ pari_APPLY_type(t_VEC,idealfactorback(nf, S, gel(x,i), 0)) }

/* shortcut for bnf = Q and p = 2 */
static GEN
bnfselmerQ(GEN S)
{
  GEN g = vec_prepend(prV_primes(S), gen_m1), e;
  long n = lg(S)-1;
  e = n? shallowconcat(zerocol(n), matid(n)): zeromat(0, 1);
  return mkvec3(g, e, const_vec(n + 1, gen_1));
}

static GEN
bnfselmer(GEN bnf, GEN S)
{
  const long p = 2;
  pari_sp av = avma;
  GEN nf = bnf_get_nf(bnf), S2, S3, e, f, e2, kerval, LS2gen, LS2fu, LS2all;
  long n = lg(S)-1, n3, n2all, r;

  S2 = bestS(bnf, S, p);
  S3 = shallowconcat(S, S2);
  LS2all = nfV_factorbackmod(nf, gel(bnfunits(bnf, S3), 1), p);
  n3 = lg(S3)-1; n2all = lg(LS2all)-1;
  LS2gen = vecslice(LS2all,1,n3);
  LS2fu  = vecslice(LS2all,n3+1, n2all-1);
  e2 = nfC_prV_val(nf, LS2gen, S2);
  kerval = Flm_ker(ZM_to_Flm(e2, p), p);
  LS2gen = nfV_zm_factorback(nf, LS2gen, kerval, p);
  e = nfC_prV_val(nf, LS2gen, S);
  e2 = ZM_divexactu(ZM_zm_mul(e2, kerval), p);
  f = prV_ZM_factorback(nf, S2, e2);
  LS2gen = shallowconcat(LS2fu, LS2gen);
  LS2gen = nfV_to_scalar_or_alg(nf, vec_prepend(LS2gen, bnf_get_tuU(bnf)));
  r = n2all - n3;
  e = shallowconcat(zeromat(n, r), e);
  f = shallowconcat(const_vec(r, gen_1), f);
  return gerepilecopy(av, mkvec3(LS2gen,e,f));
}

static GEN
get_kerval(GEN nf, GEN S2, GEN LS2gen)
{
  long i, j, lS2 = lg(S2), l = lg(LS2gen);
  GEN e = cgetg(l, t_MAT);
  for (i = 1; i < l; i++)
  {
    GEN v = cgetg(lS2, t_VECSMALL);
    for (j=1; j < lS2; j++) v[j] = odd(idealval(nf, gel(LS2gen, i), gel(S2,j)));
    gel(e, i) = Flv_to_F2v(v);
  }
  return F2m_ker(e);
}
static GEN
nf2selmer_quad(GEN nf, GEN S)
{
  pari_sp ltop = avma;
  GEN D = nf_get_disc(nf), factD = nf_get_ramified_primes(nf);
  GEN SlistQ = prV_primes(S), QS2gen, gen, Hlist, H, KerH, norms, LS2gen;
  GEN chpol, Q, kerval, S2, G, e, f, b, c, bad;
  long lS = lg(S), l, lHlist, i, j, k;

  QS2gen = vec_prepend(SlistQ, gen_m1);
  bad = ZV_sort_uniq_shallow(shallowconcat(factD, SlistQ));
  Hlist = ZV_search(bad, gen_2)? bad: vec_prepend(bad, gen_2);
  l = lg(QS2gen);
  lHlist = lg(Hlist);
  H = cgetg(l, t_MAT);
  for (j = 1; j < l; j++)
  {
    GEN v = cgetg(lHlist, t_VECSMALL);
    for (i = 1; i < lHlist; i++)
      v[i] = hilbertii(D, gel(QS2gen, j), gel(Hlist, i)) < 0;
    gel(H, j) = Flv_to_F2v(v);
  }
  KerH = F2m_ker(H); l = lg(KerH);
  norms = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(norms, i) = factorback2(QS2gen, F2c_to_ZC(gel(KerH, i)));
  LS2gen = cgetg(l, t_VEC);
  chpol = QXQ_charpoly(gel(nf_get_zk(nf), 2), nf_get_pol(nf), 0);
  b = gdivgu(negi(gel(chpol, 3)), 2);
  c = gel(chpol, 2);
  Q = mkmat3(mkcol3(gen_1, b, gen_0),
             mkcol3(b, c, gen_0),
             mkcol3(gen_0, gen_0, gen_0));
  for (k = 1; k < l; k++)
  {
    GEN sol;
    gcoeff(Q, 3, 3) = negi(gel(norms, k));
    sol = qfsolve(Q); /* must be solvable */
    sol = Q_primpart(mkcol2(gel(sol,1), gel(sol,2)));
    gel(LS2gen, k) = basistoalg(nf, sol);
  }
  if (equalis(D, -4)) G = bad;
  else
  {
    G = vecsplice(bad, ZV_search(bad, veclast(factD)));
    G = vec_prepend(G, gen_m1);
  }
  LS2gen = shallowconcat(G, LS2gen);
  l = lg(SlistQ); S2 = cgetg(l, t_VEC);
  if (l > 1)
  {
    for (i = 1; i < l; i++) gel(S2, i) = idealprimedec(nf, gel(SlistQ, i));
    S2 = setminus(shallowconcat1(S2), S);
  }
  kerval = get_kerval(nf, S2, LS2gen); l = lg(kerval);
  gen = cgetg(l, t_VEC);
  e = cgetg(l, t_MAT);
  f = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN id, ei, x = nffactorback(nf, LS2gen, F2v_to_Flv(gel(kerval, i)));
    gel(e,i) = ei = cgetg(lS, t_COL);
    for (j = 1; j < lS; j++) gel(ei,j) = stoi(idealval(nf, x, gel(S,j)));
    id = idealdiv(nf, x, idealfactorback(nf, S, gel(e,i), 0));
    if (!idealispower(nf, id, 2, &gel(f,i))) pari_err_BUG("nf2selmer_quad");
    gel(gen, i) = nf_to_scalar_or_alg(nf, x);
  }
  return gerepilecopy(ltop, mkvec3(gen, e, f));
}

static GEN
makevbnf(GEN ell, long prec)
{
  GEN vbnf, P = gel(ZX_factor(ell2pol(ell)), 1);
  long k, l = lg(P);
  vbnf = cgetg(l, t_VEC);
  for (k = 1; k < l; k++)
  {
    GEN t = gel(P,k);
    gel(vbnf,k) = degpol(t) <= 2? nfinit(t, prec): Buchall(t, nf_FORCE, prec);
  }
  return vbnf;
}

static GEN
kernorm(GEN G, GEN S, GEN pol, GEN signs)
{
  long i, j, lS = lg(S), lG = lg(G), lv = signs? lS+2: lS+1;
  GEN M = cgetg(lG, t_MAT);
  for (j = 1; j < lG; j++)
  {
    GEN v, N = QXQ_norm(gel(G,j), pol);
    gel(M, j) = v = cgetg(lv, t_VECSMALL);
    v[1] = gsigne(N) < 0;
    for (i = 1; i < lS; i++) v[i+1] = odd(Q_pvalrem(N, gel(S,i), &N));
    if (signs) v[i+1] = signs[j];
  }
  return Flm_ker(M, 2);
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
/*    FUNCTIONS FOR 2-DESCENT                  \\ */
/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
/* vector of t_VEC; return total number of entries */
static long
RgVV_nb(GEN v)
{
  long i, l = lg(v), n = 0;
  for (i = 1; i < l; i++) n += lg(gel(v,i)) - 1;
  return n;
}

/* return an Fp basis */
static GEN
elllocalimage(GEN pol, GEN K, GEN vnf, GEN p, GEN pp, GEN pts)
{
  long n, p2 = RgVV_nb(pp), prank = equaliu(p, 2)? p2: p2 - 1;
  GEN R = polrootsmodpn(pol, p), bound = addiu(p, 6);

  for(n = 1;; n++)
  {
    pari_sp btop;
    GEN x, y2, d;
    pts = Flm_image(pts, 2); if (lg(pts)-1 == prank) break;
    if ((n & 0xf) == 0) bound = mulii(bound, p);
    btop = avma;
    do
    {
      GEN r = gel(R, random_Fl(lg(R)-1)+1);
      long pprec = random_Fl(itou(gel(r, 2)) + 3) - 2; /* >= -2 */
      set_avma(btop);
      x = gadd(gel(r, 1), gmul(powis(p, pprec), randomi(bound)));
      y2 = gmul(K, poleval(pol, x));
    } while (gequal0(y2) || !Qp_issquare(y2, p));
    d = deg1pol_shallow(negi(K), gmul(K, x), 0);
    pts = vec_append(pts, kpmodsquares(vnf, d, pp));
  }
  return pts;
}

static GEN
ellLS2image(GEN pol, GEN vP, GEN K, GEN vpol, GEN vcrt)
{
  long i, l = lg(vP);
  GEN v;

  if (l == 1) return cgetg(1, t_VEC);
  v = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN P = gel(vP, i), x, t;
    if (ell_is_inf(P)) { gel(v, i) = gen_1; continue; }
    x = gel(P,1);
    t = deg1pol_shallow(negi(K), gmul(K, x), 0);
    if (gequal0(gel(P,2)))
    { /* 2-torsion, x now integer and a root of exactly one linear vpol[k]=T */
      long k, lp = lg(vpol);
      GEN a;
      for (k = 1; k < lp; k++)
      {
        GEN T = gel(vpol,k), z = gel(T,2);
        if (absequalii(x, z) && signe(x) == -signe(z)) break; /* T = X-x */
      }
      a = ZX_Z_eval(ZX_deriv(pol), x);
      t = gadd(a, gmul(gel(vcrt,k), gsub(t, a))); /* a mod T, t mod pol/T*/
    }
    gel(v, i) = t;
  }
  return v;
}

/* find small points on ell; 2-torsion points must be returned first */
static GEN
ellsearchtrivialpoints(GEN ell, GEN lim, GEN help)
{
  pari_sp av = avma;
  GEN tors2 = gel(elltors_psylow(ell,2),3);
  GEN triv = lim ? ellratpoints(ell, lim, 0): cgetg(1,t_VEC);
  if (help) triv = shallowconcat(triv, help);
  return gerepilecopy(av, shallowconcat(tors2, triv));
}

GEN
ellrankinit(GEN ell, long prec)
{
  pari_sp av = avma;
  GEN urst;
  checkell_Q(ell); ell = ellminimalbmodel(ell, &urst);
  return gerepilecopy(av, mkvec3(ell, urst, makevbnf(ell, prec)));
}

INLINE GEN
ZV_isneg(GEN x) { pari_APPLY_long(signe(gel(x,i)) < 0) }

static GEN
selmersign(GEN x, GEN vpol, GEN L)
{ pari_APPLY_same(ZV_isneg(nfeltsign(gel(x, i), RgX_rem(L, gel(vpol, i)), NULL))) }

static GEN
matselmersign(GEN vnf, GEN vpol, GEN x)
{ pari_APPLY_type(t_MAT, shallowconcat1(selmersign(vnf, vpol, gel(x, i)))) }

static GEN
_trace(GEN z, GEN T)
{
  long n = degpol(T)-1;
  if (degpol(z) < n) return gen_0;
  return gdiv(gel(z, 2+n), gel(T, 3+n));
}
static GEN
tracematrix(GEN zc, GEN b, GEN T)
{
  long i, j;
  GEN q = cgetg(4, t_MAT);
  for (j = 1; j <= 3; j++) gel(q,j) = cgetg(4, t_COL);
  for (j = 1; j <= 3; j++)
  {
    for (i = 1; i < j; i++) gcoeff(q,i,j) = gcoeff(q,j,i) =
      _trace(QXQ_mul(zc, QXQ_mul(gel(b,i), gel(b,j), T), T), T);
    gcoeff(q,i,i) = _trace(QXQ_mul(zc, QXQ_sqr(gel(b,i), T), T), T);
  }
  return q;
}

static GEN
RgXV_cxeval(GEN x, GEN r, GEN ri)
{ pari_APPLY_same(RgX_cxeval(gel(x,i), r, ri)) }

static GEN
redquadric(GEN base, GEN q2, GEN pol, GEN zc)
{
  long i, l, prec = nbits2prec(2*gexpo(q2)) + 1;
  GEN s = NULL, R = roots(pol, prec);
  l = lg(R);
  for (i = 1; i < l; ++i)
  {
    GEN r = gel(R,i), ri = gexpo(r) > 1? ginv(r): NULL;
    GEN b = RgXV_cxeval(base, r, ri), z = RgX_cxeval(zc, r, ri);
    GEN M = RgC_RgV_mulrealsym(RgV_Rg_mul(b, gabs(z, prec)), gconj(b));
    s = s? RgM_add(s, M): M;
  }
  return lllgram(s);
}

static GEN
RgX_homogenous_evaldeg(GEN P, GEN A, GEN B)
{
  long i, d = degpol(P), e = lg(B)-1;
  GEN s = gmul(gel(P, d+2), gel(B,e-d));
  for (i = d-1; i >= 0; i--)
    s = gadd(gmul(s, A), gmul(gel(B,e-i), gel(P,i+2)));
  return s;
}

static GEN
RgXV_homogenous_evaldeg(GEN x, GEN a, GEN b)
{ pari_APPLY_same(RgX_homogenous_evaldeg(gel(x,i), a, b)) }

static void
check_oncurve(GEN ell, GEN v)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++)
  {
    GEN P = gel(v, i);
    if (lg(P) != 3 || !oncurve(ell,P)) pari_err_TYPE("ellrank",P);
  }
}

static GEN
basis(GEN nf, GEN x, GEN crt, GEN pol)
{
  long i, l = lg(x);
  GEN b = cgetg(l, t_COL);
  for (i = 1; i < l; i++)
  {
    GEN z = nf_to_scalar_or_alg(nf, gel(x, i));
    gel(b, i) = grem(gsub(z, gmul(crt, z)), pol); /* z mod T, 0 mod (pol/T) */
  }
  return b;
}

static GEN
vecsmallbasis(GEN x, GEN vcrt, GEN pol)
{ pari_APPLY_same(basis(gel(x,i), nf_get_zk(gel(x,i)), gel(vcrt,i), pol)) }

static GEN
ZC_shifti(GEN x, long n) { pari_APPLY_type(t_COL, shifti(gel(x,i), n)) }

/* true nf */
static GEN
selmerbasis(GEN nf, GEN ek, GEN sqrtLS2, GEN factLS2,
            GEN badprimes, GEN crt, GEN pol)
{
  GEN sqrtzc = idealfactorback(nf, sqrtLS2, zv_neg(ek), 0);
  GEN E = ZC_shifti(ZM_zc_mul(factLS2, ek), -1);

  if (ZV_equal0(E))
    sqrtzc = idealhnf_shallow(nf, sqrtzc);
  else
    sqrtzc = idealmul(nf, sqrtzc, idealfactorback(nf, badprimes, ZC_neg(E), 0));
  return basis(nf, sqrtzc, crt, pol);
}

static long randu(void) { return random_Fl(127) - 63; }
static GEN
randS(GEN b)
{
  return gadd(gmulgs(gel(b,1), randu()),
              gadd(gmulgs(gel(b,2), randu()), gmulgs(gel(b,3), randu())));
}

static GEN
liftselmerinit(GEN expo, GEN vnf, GEN sqrtLS2, GEN factLS2,
               GEN badprimes, GEN vcrt, GEN pol)
{
  long n = lg(vnf)-1, k, t;
  GEN b = cgetg(n+1, t_VEC);
  for (k = t = 1; k <= n; k++)
  {
    GEN fak = gel(factLS2,k), ek;
    long m = lg(fak)-1;
    ek = vecslice(expo, t, t + m-1); t += m;
    gel(b,k) = selmerbasis(gel(vnf,k), ek, gel(sqrtLS2,k),
                           fak, gel(badprimes,k), gel(vcrt,k), pol);
  }
  return shallowconcat1(b);
}

static GEN
liftselmer_cover(GEN b, GEN expo, GEN LS2, GEN pol, GEN discF, GEN K)
{
  GEN P, Q, Q4, R, den, q0, q1, q2, xz, x, y, y2m, U, param, newb;
  GEN ttheta, tttheta, zc, polprime;
  GEN QM, zden;
  zc = RgXQV_factorback(LS2, expo, pol);
  if (typ(zc)==t_INT) zc = scalarpol(zc, varn(pol));
  ttheta = RgX_shift_shallow(pol,-2);
  tttheta = RgX_shift_shallow(pol, -1);
  polprime = ZX_deriv(pol);
  q2 = Q_primpart(tracematrix(zc, b, pol));
  U = redquadric(b, q2, pol, QXQ_div(zc, polprime, pol));
  q2 = qf_apply_RgM(q2, U);
  newb = RgV_RgM_mul(b, U);
  param = Q_primpart(qfparam(q2, qfsolve(q2), 1));
  param = RgM_to_RgXV_reverse(shallowtrans(param), 0);
  q1 = RgM_neg(tracematrix(QXQ_mul(zc, ttheta, pol), newb, pol));
  q1 = Q_remove_denom(qfeval(q1, param), &den);
  if (den) q1 = ZX_Z_mul(q1, den);
  if (!equali1(K)) q1 = RgX_Rg_mul(q1, K);
  QM = quartic_minim_all(q1, discF);
  q1 = gel(QM,1);
  zden = gmael(QM,2,1);
  Q = ZX_hyperellred(q1, &R);
  R = gmul(gmael(QM,2,2), R);
  if (DEBUGLEVEL>1) err_printf("  reduced quartic: Y^2 = %Ps\n", Q);
  xz = mkcol2(pol_x(0),gen_1);
  P = RgM_RgC_mul(R, xz); x = gel(P,1); y = gel(P,2);
  param = RgXV_homogenous_evaldeg(param, x, gpowers(y, 2));
  param = gmul(param, gdiv(den? mulii(K, den): K, zden));
  q0 = tracematrix(QXQ_mul(zc, tttheta, pol), newb, pol);
  x = gdiv(qfeval(q0, param), K);
  Q4 = gpowers(Q,4);
  y2m = gmul(RgX_homogenous_evaldeg(pol, x, Q4), Q);
  if (!issquareall(gdiv(y2m, K), &y))
    pari_err_BUG("liftselmer_cover");
  y = gdiv(y, gel(Q4,2));
  P = mkvec2(gdiv(gmul(K,x),pol_xn(2,1)),gdiv(gmul(gsqr(K),y),pol_xn(3,1)));
  return mkvec2(Q,P);
}

static GEN
liftselmer(GEN b, GEN expo, GEN sbase, GEN LS2, GEN pol, GEN discF, GEN K, long ntry, GEN *pt_Q)
{
  pari_sp av = avma, av2;
  long t, lim = ntry * LIM1;
  GEN ttheta, tttheta, z, polprime;
  hashtable h;
  hash_init_GEN(&h, ntry, ZX_equal, 1);
  z = RgXQV_factorback(LS2, expo, pol);
  ttheta = RgX_shift_shallow(pol,-2);
  tttheta = RgX_shift_shallow(pol, -1);
  polprime = ZX_deriv(pol); av2 = avma;
  for (t = 1; t <= ntry; t++, set_avma(av2))
  {
    GEN P, Q, Qk, R, den, q0, q1, q2, xz, x, y, zz, zc, U, param, newb, zden, QM;
    long idx;
    if (t == 1) zc = z;
    else
    {
      GEN r;
      do r = randS(sbase); while (degpol(QX_gcd(r, pol)));
      zc = QXQ_mul(z, QXQ_sqr(r, pol), pol);
    }
    q2 = Q_primpart(tracematrix(zc, b, pol));
    U = redquadric(b, q2, pol, QXQ_div(zc, polprime, pol));
    if (lg(U) < 4) continue;
    q2 = qf_apply_RgM(q2, U);
    newb = RgV_RgM_mul(b, U);

    param = Q_primpart(qfparam(q2, qfsolve(q2), 1));
    param = RgM_to_RgXV_reverse(shallowtrans(param), 0);
    q1 = RgM_neg(tracematrix(QXQ_mul(zc, ttheta, pol), newb, pol));
    q1 = Q_remove_denom(qfeval(q1, param), &den);
    if (den) q1 = ZX_Z_mul(q1, den);
    if (!equali1(K)) q1 = RgX_Rg_mul(q1, K);
    QM = quartic_minim_all(q1, discF);
    q1 = gel(QM,1);
    zden = gmael(QM,2,1);
    Q = ZX_hyperellred(q1, &R);
    R = gmul(gmael(QM,2,2), R);
    if (pt_Q) *pt_Q = Q;
    Qk = shallowcopy(Q);
    (void) ZX_canon_neg(Qk);
    if (hash_haskey_long(&h, (void*)Qk, &idx)) continue;
    hash_insert_long(&h, (void*)Qk, 1); av2 = avma;
    if (DEBUGLEVEL>1) err_printf("  reduced quartic: Y^2 = %Ps\n", Q);

    xz = projratpointxz(Q, lim, &zz);
    if (!xz)
    {
      xz = projratpointxz2(Q, lim, &zz);
      if (!xz)
      {
        if (pt_Q) return NULL; else continue;
      }
    }
    P = RgM_RgC_mul(R, xz); x = gel(P,1); y = gel(P,2);
    param = RgXV_homogenous_evaldeg(param, x, gpowers(y, 2));
    param = gmul(param, gdiv(den? mulii(K, den): K, gmul(zz, zden)));
    q0 = tracematrix(QXQ_mul(zc, tttheta, pol), newb, pol);
    x = gdiv(qfeval(q0, param), K);
    if (!issquareall(gdiv(poleval(pol, x), K), &y)) /* K y^2 = pol(x) */
      pari_err_BUG("ellrank");
    P = mkvec2(x, y);
    if (DEBUGLEVEL) err_printf("Found point: %Ps\n", P);
    if (pt_Q) *pt_Q = gen_0;
    return gerepilecopy(av, P);
  }
  return NULL;
}

static void
gtoset_inplace(GEN x)
{ gen_sort_inplace(x, (void*)&cmp_universal, cmp_nodata, NULL); }

/* FIXME: export */
static void
setlgall(GEN x, long L)
{
  long i, l = lg(x);
  for(i = 1; i < l; i++) setlg(gel(x,i), L);
}

static long
dim_selmer(GEN p, GEN pol, GEN K, GEN vnf, GEN LS2, GEN helpLS2,
           GEN *selmer, GEN *LS2chars, GEN *helpchars)
{
  pari_sp av;
  long dim, k, lvnf = lg(vnf);
  GEN X, L, LS2image, helpimage, pp = cgetg(lvnf, t_VEC);
  int pis2 = equaliu(p, 2);

  for (k = 1; k < lvnf; k++)
  {
    GEN v, nf = gel(vnf,k), PR = idealprimedec(nf, p);
    long j, l = lg(PR);
    gel(pp, k) = v = cgetg(l, t_VEC);
    for (j = 1; j < l; j++)
    {
      GEN pr = gel(PR,j);
      gel(v,j) = pis2? log_prk_init(nf, pr, 1 + 2 * pr_get_e(pr), NULL)
                     : zkmodprinit(nf, pr);
    }
  }
  LS2image = veckpmodsquares(LS2, vnf, pp);
  *LS2chars = vconcat(*LS2chars, LS2image);
  helpimage = veckpmodsquares(helpLS2, vnf, pp);
  *helpchars = vconcat(*helpchars, helpimage);
  av = avma;
  L = elllocalimage(pol, K, vnf, p, pp, helpimage);
  X = Flm_ker(shallowconcat(LS2image, L), 2); setlgall(X, lg(LS2image));
  /* intersect(LS2image, locim) = LS2image.X */
  *selmer = Flm_intersect_i(*selmer, shallowconcat(Flm_ker(LS2image,2), X), 2);
  *selmer = gerepileupto(av, Flm_image(*selmer, 2));
  dim = lg(*selmer)-1; return (dim == Flm_rank(helpimage,2))? dim: -1;
}

/* Assume there are 3 real roots, if K>0 return the smallest, otherwise the largest */
static long
get_row(GEN vnf, GEN K)
{
  long k, sK = signe(K), n = lg(vnf)-1;
  GEN R;
  if (n == 1) return sK > 0? 1: 3;
  if (n == 2)
  {
    GEN P = nf_get_pol(gel(vnf,2));
    GEN z = negi(constant_coeff(nf_get_pol(gel(vnf,1))));
    GEN y = poleval(P,z);
    GEN b = gel(P,3), a = gel(P,4);
    if (signe(y) != signe(a))
      /* 1 is between 2 and 3 */
      return sK > 0? 2: 3;
    else if (cmpii(mulii(z,mulis(a,-2)), b) == signe(a))
      return sK > 0? 1: 3;
    else
      return sK > 0? 2: 1;
  }
  R = cgetg(4, t_VEC);
  for (k = 1; k <= 3; k++) gel(R, k) = gel(nf_get_roots(gel(vnf,k)), 1);
  return sK > 0? vecindexmin(R): vecindexmax(R);
}

static GEN
ell2selmer(GEN ell, GEN ell_K, GEN help, GEN K, GEN vbnf,
           long effort, long flag, long prec)
{
  GEN KP, pol, vnf, vpol, vcrt, sbase, LS2, factLS2, sqrtLS2, signs;
  GEN selmer, helpLS2, LS2chars, helpchars, newselmer, factdisc, badprimes;
  GEN helplist, listpoints, etors2, p, covers, disc, discF;
  long i, k, n, tors2, mwrank, dim, nbpoints, lfactdisc, t, u, sha2 = 0;
  forprime_t T;

  pol = ell2pol(ell);
  help = ellsearchtrivialpoints(ell_K, flag ? NULL:muluu(LIMTRIV,effort+1), help);
  help = elltwistpoints(help, ginv(K)); /* points on K Y^2 = pol(X) */
  n = lg(vbnf) - 1; tors2 = n - 1;
  etors2 = vecslice(help,1, tors2);
  gtoset_inplace(etors2);
  KP = gel(absZ_factor(K), 1);
  disc = ZX_disc(pol);
  factdisc = mkvec3(KP, mkcol(gen_2), gel(absZ_factor(disc), 1));
  factdisc = ZV_sort_uniq_shallow(shallowconcat1(factdisc));
  discF = mkvec2(gmul(K,disc), factdisc);
  badprimes = cgetg(n+1, t_VEC);
  vnf = cgetg(n+1, t_VEC);
  vpol = cgetg(n+1, t_VEC);
  vcrt = cgetg(n+1, t_VEC);
  LS2 = cgetg(n+1, t_VEC);
  factLS2 = cgetg(n+1, t_VEC);
  sqrtLS2 = cgetg(n+1, t_VEC);
  for (k = 1; k <= n; k++)
  {
    GEN T, Tinv, id, f, sel, L, S, nf, NF = gel(vbnf,k);
    long i, lk;
    nf = (n == 1)? bnf_get_nf(NF): NF;
    gel(vnf, k) = nf;
    gel(vpol, k) = T = nf_get_pol(nf);
    Tinv = RgX_div(pol, gel(vpol, k));
    gel(vcrt, k) = QX_mul(QXQ_inv(T, Tinv), T); /* 0 mod T, 1 mod pol/T */

    id = idealadd(nf, nf_get_index(nf), ZX_deriv(T));
    f = nf_pV_to_prV(nf, KP); settyp(f, t_COL);
    f = mkvec3(gel(idealfactor(nf, Tinv), 1),
               gel(idealfactor(nf, id), 1), f);
    gel(badprimes, k) = S = gtoset(shallowconcat1(f));
    if (n == 1)
    {
      sel = bnfselmer(NF, S);
      obj_free(NF); /* units */
    }
    else if (degpol(T) == 1)
      sel = bnfselmerQ(S);
    else /* degree 2 */
      sel = nf2selmer_quad(NF, S);
    gel(LS2, k) = L = gel(sel, 1); lk = lg(L);
    gel(factLS2, k) = gel(sel, 2);
    gel(sqrtLS2, k) = gel(sel, 3);
    for (i = 1; i < lk; i++)
    {
      GEN z = gel(L,i); /* z mod T, 1 mod (pol/T) */
      gel(L,i) = grem(gadd(z, gmul(gsubsg(1,z), gel(vcrt,k))), pol);
    }
  }
  sbase = shallowconcat1(vecsmallbasis(vnf, vcrt, pol));
  if (DEBUGLEVEL>2) err_printf("   local badprimes = %Ps\n", badprimes);
  LS2 = shallowconcat1(LS2);
  helpLS2 = ellLS2image(pol, help, K, vpol, vcrt);
  LS2chars = matselmersign(vnf, vpol, LS2);
  helpchars = matselmersign(vnf, vpol, helpLS2);
  signs = NULL;
  if (signe(ell_get_disc(ell)) > 0) signs = Flm_row(LS2chars, get_row(vnf,K));
  selmer = kernorm(LS2, factdisc, pol, signs);
  forprime_init(&T, gen_2, NULL); lfactdisc = lg(factdisc); dim = -1;
  for (i = 1; dim < 0 && i < lfactdisc; i++)
    dim = dim_selmer(gel(factdisc,i), pol, K, vnf, LS2, helpLS2,
                     &selmer,&LS2chars,&helpchars);
  while (dim < 0 && Flm_rank(Flm_mul(LS2chars, selmer, 2), 2) < lg(selmer)-1)
  {
    while ((p = forprime_next(&T)) && ZV_search(factdisc, p));
    dim = dim_selmer(p, pol, K, vnf, LS2, helpLS2,
                     &selmer,&LS2chars,&helpchars);
  }
  helplist = gel(Flm_indexrank(helpchars,2), 2);
  help = shallowextract(help, helplist);
  helpchars = shallowextract(helpchars, helplist);
  helpLS2 = shallowextract(helpLS2, helplist);
  dim = lg(selmer)-1;
  if (DEBUGLEVEL) err_printf("Selmer rank: %ld\n", dim);
  newselmer = Flm_invimage(Flm_mul(LS2chars, selmer, 2), helpchars, 2);
  nbpoints = lg(help) - 1;
  if (flag==1)
  {
    GEN u = nbpoints? Flm_mul(selmer,Flm_suppl(newselmer,2), 2): selmer;
    long l = lg(u);
    GEN z = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(z,i) = RgXQV_factorback(LS2, gel(u,i), pol);
    return mkvec2(mkvec3(vnf,sbase,pol), z);
  }
  else if (flag==2)
  {
    GEN u = nbpoints ? Flm_mul(selmer,Flm_suppl(newselmer,2), 2): selmer;
    long l = lg(u);
    GEN P = cgetg(l, t_VEC), b;
    for (i = 1; i < l; i++)
    {
      b = liftselmerinit(gel(u,i), vnf, sqrtLS2, factLS2, badprimes, vcrt, pol);
      gel(P,i) = liftselmer_cover(b, gel(u,i), LS2, pol, discF, K);
    }
    return P;
  }
  listpoints = vec_lengthen(help, dim); /* points on ell */
  covers = zerovec(dim);
  for (i=1; i <= dim; i++)
  {
    GEN b, P, expo, vec = vecsmall_ei(dim, i);
    if (Flm_Flc_invimage(newselmer, vec, 2)) continue;
    expo = Flm_Flc_mul(selmer, vec, 2);
    b = liftselmerinit(expo, vnf, sqrtLS2, factLS2, badprimes, vcrt, pol);
    P = liftselmer(b, expo, sbase, LS2, pol, discF, K, 1, &gel(covers,i));
    if (P)
    {
      gel(listpoints, ++nbpoints) = P; /* new point on ell */
      gel(newselmer, nbpoints) = vec;
      setlg(newselmer, nbpoints+1);
    }
  }
  if (nbpoints < dim)
  {
    long i, j;
    GEN M = cgetg(dim+1, t_MAT), selker;
    GEN D = mulii(muliu(absi(disc), 27*4096), powiu(K,6));
    GEN FD = ZV_sort_uniq_shallow(shallowconcat1(mkvec2(mkcol3s(3,5,7), factdisc)));

    for (i = 1; i <= dim; i++) gel(M,i) = cgetg(dim+1, t_COL);
    for (i = 1; i <= dim; i++)
      for (j = 1; j <= i; j++)
      {
        GEN Q;
        if (isintzero(gel(covers,i)))
          Q = gen_0;
        else if (i==j)
          Q = quartic_findunit(D, gel(covers,i));
        else
        {
          GEN e = Flv_add(gel(selmer,i), gel(selmer,j), 2);
          GEN b = liftselmerinit(e, vnf, sqrtLS2, factLS2, badprimes, vcrt, pol);
          Q = quartic_findunit(D, gel(liftselmer_cover(b, e, LS2, pol, discF, K),1));
        }
        gmael(M,j,i) = gmael(M,i,j) = Q;
      }
    selker = F2m_to_Flm(F2m_ker(matcassels(FD, M)));
    sha2 = dim - (lg(selker)-1);
    dim = lg(selker)-1;
    for (t=1, u=1; nbpoints < dim && effort > 0; t++)
    {
      pari_sp btop = avma;
      GEN expo, b, P, vec;
      do vec = Flm_Flc_mul(selker,random_Flv(dim, 2), 2);
      while (zv_equal0(vec) || Flm_Flc_invimage(newselmer, vec, 2));
      expo = Flm_Flc_mul(selmer, vec, 2);
      b = liftselmerinit(expo, vnf, sqrtLS2, factLS2, badprimes, vcrt, pol);
      P = liftselmer(b, expo, sbase, LS2, pol, discF, K, u, NULL);
      if (P)
      {
        gel(listpoints, ++nbpoints) = P;
        gel(newselmer, nbpoints) = vec;
        setlg(newselmer, nbpoints+1);
      } else set_avma(btop);
      if (t == dim) { t = 0; u++; effort--; }
    }
  }
  setlg(listpoints, nbpoints+1);
  mwrank = nbpoints - tors2;
  if (odd(dim - nbpoints)) mwrank++;
  gtoset_inplace(listpoints);
  listpoints = setminus(listpoints, etors2);
  listpoints = elltwistpoints(listpoints, K);
  listpoints = vecellabs(ellQ_genreduce(ell_K, listpoints, NULL, prec));
  return mkvec4(utoi(mwrank), utoi(dim-tors2), utoi(sha2), listpoints);
}

GEN
ell2selmer_basis(GEN ell, GEN *cb, long prec)
{
  GEN E = ellminimalbmodel(ell, cb);
  GEN S = ell2selmer(E, E, NULL, gen_1, makevbnf(E, prec), 0, 1, prec);
  obj_free(E); return S;
}

static void
ellrank_get_nudur(GEN E, GEN F, GEN *nu, GEN *du, GEN *r)
{
  GEN ea2 = ell_get_a2(E), ea2t = ell_get_a2(F);
  GEN ec4 = ell_get_c4(E), ec4t = ell_get_c4(F);
  GEN ec6 = ell_get_c6(E), ec6t = ell_get_c6(F);
  GEN N, D, d;
  if (signe(ec4)==0)
  {
    if (!Z_ispowerall(mulii(ec6,sqri(ec6t)),3,&N))
      pari_err_TYPE("ellrank",F);
    D = ec6t;
  }
  else if (signe(ec6)==0)
  {
    if (!Z_issquareall(mulii(ec4,ec4t),&N))
      pari_err_TYPE("ellrank",F);
    D = ec4t;
  }
  else
  {
    GEN d46 = mulii(gcdii(ec4,ec4t),gcdii(ec6,ec6t));
    N = diviiexact(mulii(ec6,ec4t),d46);
    D = diviiexact(mulii(ec6t,ec4),d46);
  }
  d = gcdii(N, D);
  *nu = diviiexact(N, d); *du = diviiexact(D, d);
  *r  = diviuexact(subii(mulii(*nu,ea2t),mulii(*du,ea2)),3);
}

static GEN
ellrank_flag(GEN e, long effort, GEN help, long flag, long prec)
{
  pari_sp ltop = avma;
  GEN urst, v, vbnf, eK;
  GEN et = NULL, K = gen_1, nu = NULL, du = NULL, r = NULL, urstK = NULL;
  long newell = 0;

  if (lg(e)==3 && typ(e)==t_VEC) { et = gel(e,2); e = gel(e,1); }
  if (lg(e)==4 && typ(e)==t_VEC)
  {
    vbnf = gel(e,3); urst = gel(e,2);
    e = gel(e,1); checkell_Q(e);
    if (!ell_is_integral(e)) pari_err_TYPE("ellrank [nonintegral model]",e);
    if (signe(ell_get_a1(e))) pari_err_TYPE("ellrank [a1 != 0]", e);
    if (signe(ell_get_a3(e))) pari_err_TYPE("ell2rank [a3 != 0]", e);
  }
  else
  {
    checkell_Q(e);
    e = ellminimalbmodel(e, &urst);
    newell = 1;
    vbnf = makevbnf(e, prec);
  }
  if (et)
  {
    checkell_Q(et);
    if (!gequal(ell_get_j(et),ell_get_j(e))) pari_err_TYPE("ellrank",et);
    et = ellminimalbmodel(et, &urst);
    ellrank_get_nudur(e, et, &nu, &du, &r);
    K = mulii(nu, du);
    urstK = mkvec4(nu, mulii(nu,r), gen_0, gen_0);
  }
  if (help)
  {
    if (urst) help = ellchangepoint(help, urst);
    if (et) help = ellchangepointinv(help, urstK);
  }
  eK = elltwistequation(e, K);
  /* help is a vector of rational points [x,y] satisfying K y^2 = pol(x) */
  /* [Kx, K^2y] is on eK: Y^2 = K^3 pol(X/K)  */
  if (help) check_oncurve(eK, help);
  v = ell2selmer(e, eK, help, K, vbnf, effort, flag, prec);
  if (flag==0)
  {
    if (et)   gel(v,4) = ellchangepoint(gel(v,4), urstK);
    if (urst) gel(v,4) = ellchangepointinv(gel(v,4), urst);
  }
  else
  {
    long i, l = lg(v);
    for (i = 1; i < l; i++)
    {
      if (et)   gmael(v,i,2) = ellchangepoint(gmael(v,i,2), urstK);
      if (urst) gmael(v,i,2) = ellchangepointinv(gmael(v,i,2), urst);
    }
  }
  if (newell) obj_free(e);
  if (eK != e) obj_free(eK);
  return gerepilecopy(ltop, v);
}

GEN
ellrank(GEN e, long effort, GEN help, long prec)
{
  return ellrank_flag(e, effort, help, 0, prec);
}

GEN
ell2cover(GEN ell, long prec)
{
  return ellrank_flag(ell, 0, NULL, 2, prec);
}
