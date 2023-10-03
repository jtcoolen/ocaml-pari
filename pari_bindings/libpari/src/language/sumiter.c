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

GEN
iferrpari(GEN a, GEN b, GEN c)
{
  GEN res;
  struct pari_evalstate state;
  evalstate_save(&state);
  pari_CATCH(CATCH_ALL)
  {
    GEN E;
    if (!b&&!c) return gnil;
    E = evalstate_restore_err(&state);
    if (c)
    {
      push_lex(E,c);
      res = closure_evalnobrk(c);
      pop_lex(1);
      if (gequal0(res))
        pari_err(0, E);
    }
    if (!b) return gnil;
    push_lex(E,b);
    res = closure_evalgen(b);
    pop_lex(1);
    return res;
  } pari_TRY {
    res = closure_evalgen(a);
  } pari_ENDCATCH;
  return res;
}

/********************************************************************/
/**                                                                **/
/**                        ITERATIONS                              **/
/**                                                                **/
/********************************************************************/

static void
forparii(GEN a, GEN b, GEN code)
{
  pari_sp av, av0 = avma;
  GEN aa;
  if (gcmp(b,a) < 0) return;
  if (typ(b) != t_INFINITY) b = gfloor(b);
  aa = a = setloop(a);
  av=avma;
  push_lex(a,code);
  while (gcmp(a,b) <= 0)
  {
    closure_evalvoid(code); if (loop_break()) break;
    a = get_lex(-1);
    if (a == aa)
    {
      a = incloop(a);
      if (a != aa) { set_lex(-1,a); aa = a; }
    }
    else
    { /* 'code' modified a ! Be careful (and slow) from now on */
      a = gaddgs(a,1);
      if (gc_needed(av,1))
      {
        if (DEBUGMEM>1) pari_warn(warnmem,"forparii");
        a = gerepileupto(av,a);
      }
      set_lex(-1,a);
    }
  }
  pop_lex(1);  set_avma(av0);
}

void
forpari(GEN a, GEN b, GEN code)
{
  pari_sp ltop=avma, av;
  if (typ(a) == t_INT) { forparii(a,b,code); return; }
  b = gcopy(b); /* Kludge to work-around the a+(a=2) bug */
  av=avma;
  push_lex(a,code);
  while (gcmp(a,b) <= 0)
  {
    closure_evalvoid(code); if (loop_break()) break;
    a = get_lex(-1); a = gaddgs(a,1);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"forpari");
      a = gerepileupto(av,a);
    }
    set_lex(-1, a);
  }
  pop_lex(1); set_avma(ltop);
}

void
foreachpari(GEN x, GEN code)
{
  long i, l;
  switch(typ(x))
  {
    case t_LIST:
      x = list_data(x); /* FALL THROUGH */
      if (!x) return;
    case t_MAT: case t_VEC: case t_COL:
      break;
    default:
      pari_err_TYPE("foreach",x);
      return; /*LCOV_EXCL_LINE*/
  }
  clone_lock(x); l = lg(x);
  push_lex(gen_0,code);
  for (i = 1; i < l; i++)
  {
    set_lex(-1, gel(x,i));
    closure_evalvoid(code); if (loop_break()) break;
  }
  pop_lex(1); clone_unlock_deep(x);
}

/* is it better to sieve [a,b] or to factor individually ? */
static int
no_sieve(ulong a, ulong b)
{ return b - a < usqrt(b) / tridiv_boundu(b); }

/* 0 < a <= b. Using small consecutive chunks to 1) limit memory use, 2) allow
 * cheap early abort */
static int
forfactoredpos(ulong a, ulong b, GEN code)
{
  ulong x1, step = maxuu(2 * usqrt(b), 1024);
  pari_sp av = avma;
  if (no_sieve(a, b))
  {
    ulong n;
    for (n = a; n <= b; n++, set_avma(av))
    {
      GEN m = factoru(n);
      set_lex(-1, mkvec2(utoipos(n), Flm_to_ZM(m)));
      closure_evalvoid(code); if (loop_break()) return 1;
    }
    return 0;
  }
  for(x1 = a;; x1 += step, set_avma(av))
  { /* beware overflow, fuse last two bins (avoid a tiny remainder) */
    ulong j, lv, x2 = (b >= 2*step && b - 2*step >= x1)? x1-1 + step: b;
    GEN v = vecfactoru_i(x1, x2);
    lv = lg(v);
    for (j = 1; j < lv; j++)
    {
      ulong n = x1-1 + j;
      set_lex(-1, mkvec2(utoipos(n), Flm_to_ZM(gel(v,j))));
      closure_evalvoid(code);
      if (loop_break()) return 1;
    }
    if (x2 == b) break;
    set_lex(-1, gen_0);
  }
  return 0;
}

/* vector of primes to squarefree factorization */
static GEN
zv_to_ZM(GEN v)
{ return mkmat2(zc_to_ZC(v), const_col(lg(v)-1,gen_1)); }
/* vector of primes to negative squarefree factorization */
static GEN
zv_to_mZM(GEN v)
{
  long i, l = lg(v);
  GEN w = cgetg(l+1, t_COL);
  gel(w,1) = gen_m1; for (i = 1; i < l; i++) gel(w,i+1) = utoipos(v[i]);
  return mkmat2(w, const_col(l,gen_1));
}
/* 0 <= a <= b. Using small consecutive chunks to 1) limit memory use, 2) allow
 * cheap early abort */
static void
forsquarefreepos(ulong a, ulong b, GEN code)
{
  const ulong step = maxuu(1024, 2 * usqrt(b));
  pari_sp av = avma;
  ulong x1;
  if (no_sieve(a, b))
  {
    ulong n;
    for (n = a; n <= b; n++, set_avma(av))
    {
      GEN m = factoru(n);
      if (!uissquarefree_fact(m)) continue;
      set_lex(-1, mkvec2(utoipos(n), Flm_to_ZM(m)));
      closure_evalvoid(code); if (loop_break()) return;
    }
    return;
  }
  for(x1 = a;; x1 += step, set_avma(av))
  { /* beware overflow, fuse last two bins (avoid a tiny remainder) */
    ulong j, lv, x2 = (b >= 2*step && b - 2*step >= x1)? x1-1 + step: b;
    GEN v = vecfactorsquarefreeu(x1, x2);
    lv = lg(v);
    for (j = 1; j < lv; j++) if (gel(v,j))
    {
      ulong n = x1-1 + j;
      set_lex(-1, mkvec2(utoipos(n), zv_to_ZM(gel(v,j))));
      closure_evalvoid(code); if (loop_break()) return;
    }
    if (x2 == b) break;
    set_lex(-1, gen_0);
  }
}
/* 0 <= a <= b. Loop from -b, ... -a through squarefree integers */
static void
forsquarefreeneg(ulong a, ulong b, GEN code)
{
  const ulong step = maxuu(1024, 2 * usqrt(b));
  pari_sp av = avma;
  ulong x2;
  if (no_sieve(a, b))
  {
    ulong n;
    for (n = b; n >= a; n--, set_avma(av))
    {
      GEN m = factoru(n);
      if (!uissquarefree_fact(m)) continue;
      set_lex(-1, mkvec2(utoineg(n), zv_to_mZM(gel(m,1))));
      closure_evalvoid(code); if (loop_break()) return;
    }
    return;
  }
  for(x2 = b;; x2 -= step, set_avma(av))
  { /* beware overflow, fuse last two bins (avoid a tiny remainder) */
    ulong j, x1 = (x2 >= 2*step && x2-2*step >= a)? x2+1 - step: a;
    GEN v = vecfactorsquarefreeu(x1, x2);
    for (j = lg(v)-1; j > 0; j--) if (gel(v,j))
    {
      ulong n = x1-1 + j;
      set_lex(-1, mkvec2(utoineg(n), zv_to_mZM(gel(v,j))));
      closure_evalvoid(code); if (loop_break()) return;
    }
    if (x1 == a) break;
    set_lex(-1, gen_0);
  }
}
void
forsquarefree(GEN a, GEN b, GEN code)
{
  pari_sp av = avma;
  long s;
  if (typ(a) != t_INT) pari_err_TYPE("forsquarefree", a);
  if (typ(b) != t_INT) pari_err_TYPE("forsquarefree", b);
  if (cmpii(a,b) > 0) return;
  s = signe(a); push_lex(NULL,code);
  if (s < 0)
  {
    if (signe(b) <= 0)
      forsquarefreeneg(itou(b), itou(a), code);
    else
    {
      forsquarefreeneg(1, itou(a), code);
      forsquarefreepos(1, itou(b), code);
    }
  }
  else
    forsquarefreepos(itou(a), itou(b), code);
  pop_lex(1); set_avma(av);
}

/* convert factoru(n) to factor(-n); M pre-allocated factorization matrix
 * with (-1)^1 already set */
static void
Flm2negfact(GEN v, GEN M)
{
  GEN p = gel(v,1), e = gel(v,2), P = gel(M,1), E = gel(M,2);
  long i, l = lg(p);
  for (i = 1; i < l; i++)
  {
    gel(P,i+1) = utoipos(p[i]);
    gel(E,i+1) = utoipos(e[i]);
  }
  setlg(P,l+1);
  setlg(E,l+1);
}
/* 0 < a <= b, from -b to -a */
static int
forfactoredneg(ulong a, ulong b, GEN code)
{
  ulong x2, step = maxuu(2 * usqrt(b), 1024);
  GEN P, E, M;
  pari_sp av;

  P = cgetg(18, t_COL); gel(P,1) = gen_m1;
  E = cgetg(18, t_COL); gel(E,1) = gen_1;
  M = mkmat2(P,E);
  av = avma;
  if (no_sieve(a, b))
  {
    ulong n;
    for (n = b; n >= a; n--, set_avma(av))
    {
      GEN m = factoru(n);
      Flm2negfact(m, M);
      set_lex(-1, mkvec2(utoineg(n), M));
      closure_evalvoid(code); if (loop_break()) return 1;
    }
    return 0;
  }
  for (x2 = b;; x2 -= step, set_avma(av))
  { /* beware overflow, fuse last two bins (avoid a tiny remainder) */
    ulong j, x1 = (x2 >= 2*step && x2-2*step >= a)? x2+1 - step: a;
    GEN v = vecfactoru_i(x1, x2);
    for (j = lg(v)-1; j; j--)
    { /* run backward: from factor(x1..x2) to factor(-x2..-x1) */
      ulong n = x1-1 + j;
      Flm2negfact(gel(v,j), M);
      set_lex(-1, mkvec2(utoineg(n), M));
      closure_evalvoid(code); if (loop_break()) return 1;
    }
    if (x1 == a) break;
    set_lex(-1, gen_0);
  }
  return 0;
}
static int
eval0(GEN code)
{
  pari_sp av = avma;
  set_lex(-1, mkvec2(gen_0, mkmat2(mkcol(gen_0),mkcol(gen_1))));
  closure_evalvoid(code); set_avma(av);
  return loop_break();
}
void
forfactored(GEN a, GEN b, GEN code)
{
  pari_sp av = avma;
  long sa, sb, stop = 0;
  if (typ(a) != t_INT) pari_err_TYPE("forfactored", a);
  if (typ(b) != t_INT) pari_err_TYPE("forfactored", b);
  if (cmpii(a,b) > 0) return;
  push_lex(NULL,code);
  sa = signe(a);
  sb = signe(b);
  if (sa < 0)
  {
    stop = forfactoredneg((sb < 0)? uel(b,2): 1UL, itou(a), code);
    if (!stop && sb >= 0) stop = eval0(code);
    if (!stop && sb > 0) forfactoredpos(1UL, b[2], code);
  }
  else
  {
    if (!sa) stop = eval0(code);
    if (!stop && sb) forfactoredpos(sa? uel(a,2): 1UL, itou(b), code);
  }
  pop_lex(1); set_avma(av);
}
void
whilepari(GEN a, GEN b)
{
  pari_sp av = avma;
  for(;;)
  {
    GEN res = closure_evalnobrk(a);
    if (gequal0(res)) break;
    set_avma(av);
    closure_evalvoid(b); if (loop_break()) break;
  }
  set_avma(av);
}

void
untilpari(GEN a, GEN b)
{
  pari_sp av = avma;
  for(;;)
  {
    GEN res;
    closure_evalvoid(b); if (loop_break()) break;
    res = closure_evalnobrk(a);
    if (!gequal0(res)) break;
    set_avma(av);
  }
  set_avma(av);
}

static int negcmp(GEN x, GEN y) { return gcmp(y,x); }

void
forstep(GEN a, GEN b, GEN s, GEN code)
{
  long ss, i;
  pari_sp av, av0 = avma;
  GEN v = NULL;
  int (*cmp)(GEN,GEN);

  b = gcopy(b);
  s = gcopy(s); av = avma;
  switch(typ(s))
  {
    case t_VEC: case t_COL: ss = gsigne(vecsum(s)); v = s; break;
    case t_INTMOD:
      if (typ(a) != t_INT) a = gceil(a);
      a = addii(a, modii(subii(gel(s,2),a), gel(s,1)));
      s = gel(s,1);
    default: ss = gsigne(s);
  }
  if (!ss) pari_err_DOMAIN("forstep","step","=",gen_0,s);
  cmp = (ss > 0)? &gcmp: &negcmp;
  i = 0;
  push_lex(a,code);
  while (cmp(a,b) <= 0)
  {
    closure_evalvoid(code); if (loop_break()) break;
    if (v)
    {
      if (++i >= lg(v)) i = 1;
      s = gel(v,i);
    }
    a = get_lex(-1); a = gadd(a,s);

    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"forstep");
      a = gerepileupto(av,a);
    }
    set_lex(-1,a);
  }
  pop_lex(1); set_avma(av0);
}

static void
_fordiv(GEN a, GEN code, GEN (*D)(GEN))
{
  pari_sp av = avma;
  long i, l;
  GEN t = D(a);
  push_lex(gen_0,code); l = lg(t);
  for (i=1; i<l; i++)
  {
    set_lex(-1,gel(t,i));
    closure_evalvoid(code); if (loop_break()) break;
  }
  pop_lex(1); set_avma(av);
}
void
fordiv(GEN a, GEN code) { return _fordiv(a, code, &divisors); }
void
fordivfactored(GEN a, GEN code) { return _fordiv(a, code, &divisors_factored); }

/* Embedded for loops:
 *   fl = 0: execute ch (a), where a = (ai) runs through all n-uplets in
 *     [m1,M1] x ... x [mn,Mn]
 *   fl = 1: impose a1 <= ... <= an
 *   fl = 2:        a1 <  ... <  an
 */
/* increment and return d->a [over integers]*/
static GEN
_next_i(forvec_t *d)
{
  long i = d->n;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (;;) {
    if (cmpii(d->a[i], d->M[i]) < 0) {
      d->a[i] = incloop(d->a[i]);
      return (GEN)d->a;
    }
    d->a[i] = resetloop(d->a[i], d->m[i]);
    if (--i <= 0) return NULL;
  }
}
/* increment and return d->a [generic]*/
static GEN
_next(forvec_t *d)
{
  long i = d->n;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (;;) {
    d->a[i] = gaddgs(d->a[i], 1);
    if (gcmp(d->a[i], d->M[i]) <= 0) return (GEN)d->a;
    d->a[i] = d->m[i];
    if (--i <= 0) return NULL;
  }
}

/* nondecreasing order [over integers] */
static GEN
_next_le_i(forvec_t *d)
{
  long i = d->n;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (;;) {
    if (cmpii(d->a[i], d->M[i]) < 0)
    {
      d->a[i] = incloop(d->a[i]);
      /* m[i] < a[i] <= M[i] <= M[i+1] */
      while (i < d->n)
      {
        GEN t;
        i++;
        if (cmpii(d->a[i-1], d->a[i]) <= 0) continue;
        /* a[i] < a[i-1] <= M[i-1] <= M[i] */
        t = d->a[i-1]; if (cmpii(t, d->m[i]) < 0) t = d->m[i];
        d->a[i] = resetloop(d->a[i], t);/*a[i]:=max(a[i-1],m[i])*/
      }
      return (GEN)d->a;
    }
    d->a[i] = resetloop(d->a[i], d->m[i]);
    if (--i <= 0) return NULL;
  }
}
/* nondecreasing order [generic] */
static GEN
_next_le(forvec_t *d)
{
  long i = d->n;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (;;) {
    d->a[i] = gaddgs(d->a[i], 1);
    if (gcmp(d->a[i], d->M[i]) <= 0)
    {
      while (i < d->n)
      {
        GEN c;
        i++;
        if (gcmp(d->a[i-1], d->a[i]) <= 0) continue;
        /* M[i] >= M[i-1] >= a[i-1] > a[i] */
        c = gceil(gsub(d->a[i-1], d->a[i]));
        d->a[i] = gadd(d->a[i], c);
        /* a[i-1] <= a[i] < M[i-1] + 1 => a[i] < M[i]+1 => a[i] <= M[i] */
      }
      return (GEN)d->a;
    }
    d->a[i] = d->m[i];
    if (--i <= 0) return NULL;
  }
}
/* strictly increasing order [over integers] */
static GEN
_next_lt_i(forvec_t *d)
{
  long i = d->n;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (;;) {
    if (cmpii(d->a[i], d->M[i]) < 0)
    {
      d->a[i] = incloop(d->a[i]);
      /* m[i] < a[i] <= M[i] < M[i+1] */
      while (i < d->n)
      {
        pari_sp av;
        GEN t;
        i++;
        if (cmpii(d->a[i-1], d->a[i]) < 0) continue;
        av = avma;
        /* M[i] > M[i-1] >= a[i-1] */
        t = addiu(d->a[i-1],1); if (cmpii(t, d->m[i]) < 0) t = d->m[i];
        d->a[i] = resetloop(d->a[i], t);/*a[i]:=max(a[i-1]+1,m[i]) <= M[i]*/
        set_avma(av);
      }
      return (GEN)d->a;
    }
    d->a[i] = resetloop(d->a[i], d->m[i]);
    if (--i <= 0) return NULL;
  }
}
/* strictly increasing order [generic] */
static GEN
_next_lt(forvec_t *d)
{
  long i = d->n;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (;;) {
    d->a[i] = gaddgs(d->a[i], 1);
    if (gcmp(d->a[i], d->M[i]) <= 0)
    {
      while (i < d->n)
      {
        GEN c;
        i++;
        if (gcmp(d->a[i-1], d->a[i]) < 0) continue;
        /* M[i] > M[i-1] >= a[i-1] >= a[i] */
        c = addiu(gfloor(gsub(d->a[i-1], d->a[i])), 1); /* > a[i-1] - a[i] */
        d->a[i] = gadd(d->a[i], c);
        /* a[i-1] < a[i] <= M[i-1] + 1 => a[i] < M[i]+1 => a[i] <= M[i] */
      }
      return (GEN)d->a;
    }
    d->a[i] = d->m[i];
    if (--i <= 0) return NULL;
  }
}

/* on Z^n /(cyc Z^n) [over integers]
 * torsion (cyc>0) and free (cyc=0) components may be interleaved */
static GEN
_next_mod_cyc(forvec_t *d)
{ /* keep free components indices t1 < t2 last nonzero < t3 */
  long t, t1 = 0, t2 = 0, t3 = 0;
  if (d->first) { d->first = 0; return (GEN)d->a; }
  for (t = d->n; t > 0; t--)
  {
    if (signe(d->M[t]) > 0)
    { /* torsion component */
      d->a[t] = incloop(d->a[t]);
      if (cmpii(d->a[t], d->M[t]) < 0) return (GEN)d->a;
      d->a[t] = resetloop(d->a[t], gen_0);
    }
    else
    { /* set or update t1,t2,t3 */
      if (t2 && !t1) t1 = t;
      if (!t2 && signe(d->a[t])) t2 = t;
      if (!t2) t3 = t;
    }
  }
  if (!t3 && !t2) return NULL; /* no free component, stop */
  if (!t2) d->a[t3] = resetloop(d->a[t3], gen_m1);
  else if (!t3 && signe(d->a[t2]) < 0) togglesign(d->a[t2]);
  else if (signe(d->a[t2]) < 0)
  {
    d->a[t2] = incloop(d->a[t2]);
    d->a[t3] = resetloop(d->a[t3], gen_m1);
  }
  else if (!t1) { d->a[t2] = incloop(d->a[t2]); togglesign(d->a[t2]); }
  else
  {
    if (signe(d->a[t1]) < 0)
    { d->a[t2] = incloop(d->a[t2]); togglesign(d->a[t2]); }
    else
    { togglesign(d->a[t2]); d->a[t2] = incloop(d->a[t2]); }
    d->a[t1] = incloop(d->a[t1]);
  }
  return (GEN)d->a;
}
/* for forvec(v=[],) */
static GEN
_next_void(forvec_t *d)
{
  if (d->first) { d->first = 0; return (GEN)d->a; }
  return NULL;
}
static int
RgV_is_ZV_nonneg(GEN x)
{
  long i;
  for (i = lg(x)-1; i > 0; i--)
    if (typ(gel(x,i)) != t_INT || signe(gel(x, i)) < 0) return 0;
  return 1;
}
/* x assumed to be cyc vector, l>1 */
static int
forvec_mod_cyc_init(forvec_t *d, GEN x)
{
  long i, tx = typ(x), l = lg(x);
  d->a = (GEN*)cgetg(l,tx); /* current */
  d->M = (GEN*)cgetg(l,tx); /* cyc */
  for (i = 1; i < l; i++)
  {
    d->a[i] = setloop(gen_0);
    d->M[i] = setloop(gel(x, i));
  }
  d->first = 1;
  d->n = l-1;
  d->m = NULL;
  d->next = &_next_mod_cyc;
  return 1;
}

/* Initialize minima (m) and maxima (M); guarantee M[i] - m[i] integer and
 *   if flag = 1: m[i-1] <= m[i] <= M[i] <= M[i+1]
 *   if flag = 2: m[i-1] <  m[i] <= M[i] <  M[i+1],
 * for all i */
int
forvec_init(forvec_t *d, GEN x, long flag)
{
  long i, tx = typ(x), l = lg(x), t = t_INT;
  if (!is_vec_t(tx)) pari_err_TYPE("forvec [not a vector]", x);
  if (l > 1 && RgV_is_ZV_nonneg(x))
      return forvec_mod_cyc_init(d, x);
  d->first = 1;
  d->n = l - 1;
  d->a = (GEN*)cgetg(l,tx);
  d->m = (GEN*)cgetg(l,tx);
  d->M = (GEN*)cgetg(l,tx);
  if (l == 1) { d->next = &_next_void; return 1; }
  for (i = 1; i < l; i++)
  {
    GEN a, e = gel(x,i), m = gel(e,1), M = gel(e,2);
    tx = typ(e);
    if (! is_vec_t(tx) || lg(e)!=3)
      pari_err_TYPE("forvec [expected vector not of type [min,MAX]]",e);
    if (typ(m) != t_INT) t = t_REAL;
    if (i > 1) switch(flag)
    {
      case 1: /* a >= m[i-1] - m */
        a = gceil(gsub(d->m[i-1], m));
        if (typ(a) != t_INT) pari_err_TYPE("forvec",a);
        if (signe(a) > 0) m = gadd(m, a); else m = gcopy(m);
        break;
      case 2: /* a > m[i-1] - m */
        a = gfloor(gsub(d->m[i-1], m));
        if (typ(a) != t_INT) pari_err_TYPE("forvec",a);
        a = addiu(a, 1);
        if (signe(a) > 0) m = gadd(m, a); else m = gcopy(m);
        break;
      default: m = gcopy(m);
        break;
    }
    M = gadd(m, gfloor(gsub(M,m))); /* ensure M-m is an integer */
    if (gcmp(m,M) > 0) { d->a = NULL; d->next = &_next; return 0; }
    d->m[i] = m;
    d->M[i] = M;
  }
  if (flag == 1) for (i = l-2; i >= 1; i--)
  {
    GEN M = d->M[i], a = gfloor(gsub(d->M[i+1], M));
    if (typ(a) != t_INT) pari_err_TYPE("forvec",a);
    /* M[i]+a <= M[i+1] */
    if (signe(a) < 0) d->M[i] = gadd(M, a);
  }
  else if (flag == 2) for (i = l-2; i >= 1; i--)
  {
    GEN M = d->M[i], a = gceil(gsub(d->M[i+1], M));
    if (typ(a) != t_INT) pari_err_TYPE("forvec",a);
    a = subiu(a, 1);
    /* M[i]+a < M[i+1] */
    if (signe(a) < 0) d->M[i] = gadd(M, a);
  }
  if (t == t_INT) {
    for (i = 1; i < l; i++) {
      d->a[i] = setloop(d->m[i]);
      if (typ(d->M[i]) != t_INT) d->M[i] = gfloor(d->M[i]);
    }
  } else {
    for (i = 1; i < l; i++) d->a[i] = d->m[i];
  }
  switch(flag)
  {
    case 0: d->next = t==t_INT? &_next_i:    &_next; break;
    case 1: d->next = t==t_INT? &_next_le_i: &_next_le; break;
    case 2: d->next = t==t_INT? &_next_lt_i: &_next_lt; break;
    default: pari_err_FLAG("forvec");
  }
  return 1;
}
GEN
forvec_next(forvec_t *d) { return d->next(d); }

void
forvec(GEN x, GEN code, long flag)
{
  pari_sp av = avma;
  forvec_t T;
  GEN v;
  if (!forvec_init(&T, x, flag)) { set_avma(av); return; }
  push_lex((GEN)T.a, code);
  while ((v = forvec_next(&T)))
  {
    closure_evalvoid(code);
    if (loop_break()) break;
  }
  pop_lex(1); set_avma(av);
}

/********************************************************************/
/**                                                                **/
/**                              SUMS                              **/
/**                                                                **/
/********************************************************************/

GEN
somme(GEN a, GEN b, GEN code, GEN x)
{
  pari_sp av, av0 = avma;
  GEN p1;

  if (typ(a) != t_INT) pari_err_TYPE("sum",a);
  if (!x) x = gen_0;
  if (gcmp(b,a) < 0) return gcopy(x);

  b = gfloor(b);
  a = setloop(a);
  av=avma;
  push_lex(a,code);
  for(;;)
  {
    p1 = closure_evalnobrk(code);
    x=gadd(x,p1); if (cmpii(a,b) >= 0) break;
    a = incloop(a);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"sum");
      x = gerepileupto(av,x);
    }
    set_lex(-1,a);
  }
  pop_lex(1); return gerepileupto(av0,x);
}

static GEN
sum_init(GEN x0, GEN t)
{
  long tp = typ(t);
  GEN x;
  if (is_vec_t(tp))
  {
    x = const_vec(lg(t)-1, x0);
    settyp(x, tp);
  }
  else
    x = x0;
  return x;
}

GEN
suminf_bitprec(void *E, GEN (*eval)(void *, GEN), GEN a, long bit)
{
  long fl = 0, G = bit + 1;
  pari_sp av0 = avma, av;
  GEN x = NULL, _1;

  if (typ(a) != t_INT) pari_err_TYPE("suminf",a);
  a = setloop(a); av = avma;
  for(;;)
  {
    GEN t = eval(E, a);
    if (!x) _1 = x = sum_init(real_1_bit(bit), t);

    x = gadd(x,t);
    if (!gequal0(t) && gexpo(t) > gexpo(x)-G)
      fl = 0;
    else if (++fl == 3)
      break;
    a = incloop(a);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"suminf");
      gerepileall(av,2, &x, &_1);
    }
  }
  return gerepileupto(av0, gsub(x, _1));
}
GEN
suminf(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{ return suminf_bitprec(E, eval, a, prec2nbits(prec)); }
GEN
suminf0_bitprec(GEN a, GEN code, long bit)
{ EXPR_WRAP(code, suminf_bitprec(EXPR_ARG, a, bit)); }

GEN
sumdivexpr(GEN num, GEN code)
{
  pari_sp av = avma;
  GEN y = gen_0, t = divisors(num);
  long i, l = lg(t);

  push_lex(gen_0, code);
  for (i=1; i<l; i++)
  {
    set_lex(-1,gel(t,i));
    y = gadd(y, closure_evalnobrk(code));
  }
  pop_lex(1); return gerepileupto(av,y);
}

GEN
sumdivmultexpr(void *D, GEN (*fun)(void*, GEN), GEN num)
{
  pari_sp av = avma;
  GEN y = gen_1, P,E;
  int isint = divisors_init(num, &P,&E);
  long i, l = lg(P);
  GEN (*mul)(GEN,GEN);

  if (l == 1) return gc_const(av, gen_1);
  mul = isint? mulii: gmul;
  for (i=1; i<l; i++)
  {
    GEN p = gel(P,i), q = p, z = gen_1;
    long j, e = E[i];
    for (j = 1; j <= e; j++, q = mul(q, p))
    {
      z = gadd(z, fun(D, q));
      if (j == e) break;
    }
    y = gmul(y, z);
  }
  return gerepileupto(av,y);
}

GEN
sumdivmultexpr0(GEN num, GEN code)
{ EXPR_WRAP(code, sumdivmultexpr(EXPR_ARG, num)) }

/********************************************************************/
/**                                                                **/
/**                           PRODUCTS                             **/
/**                                                                **/
/********************************************************************/

GEN
produit(GEN a, GEN b, GEN code, GEN x)
{
  pari_sp av, av0 = avma;
  GEN p1;

  if (typ(a) != t_INT) pari_err_TYPE("prod",a);
  if (!x) x = gen_1;
  if (gcmp(b,a) < 0) return gcopy(x);

  b = gfloor(b);
  a = setloop(a);
  av=avma;
  push_lex(a,code);
  for(;;)
  {
    p1 = closure_evalnobrk(code);
    x = gmul(x,p1); if (cmpii(a,b) >= 0) break;
    a = incloop(a);
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"prod");
      x = gerepileupto(av,x);
    }
    set_lex(-1,a);
  }
  pop_lex(1); return gerepileupto(av0,x);
}

GEN
prodinf(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{
  pari_sp av0 = avma, av;
  long fl,G;
  GEN p1,x = real_1(prec);

  if (typ(a) != t_INT) pari_err_TYPE("prodinf",a);
  a = setloop(a);
  av = avma;
  fl=0; G = -prec2nbits(prec)-5;
  for(;;)
  {
    p1 = eval(E, a); if (gequal0(p1)) { x = p1; break; }
    x = gmul(x,p1); a = incloop(a);
    p1 = gsubgs(p1, 1);
    if (gequal0(p1) || gexpo(p1) <= G) { if (++fl==3) break; } else fl=0;
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"prodinf");
      x = gerepileupto(av,x);
    }
  }
  return gerepilecopy(av0,x);
}
GEN
prodinf1(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{
  pari_sp av0 = avma, av;
  long fl,G;
  GEN p1,p2,x = real_1(prec);

  if (typ(a) != t_INT) pari_err_TYPE("prodinf1",a);
  a = setloop(a);
  av = avma;
  fl=0; G = -prec2nbits(prec)-5;
  for(;;)
  {
    p2 = eval(E, a); p1 = gaddgs(p2,1);
    if (gequal0(p1)) { x = p1; break; }
    x = gmul(x,p1); a = incloop(a);
    if (gequal0(p2) || gexpo(p2) <= G) { if (++fl==3) break; } else fl=0;
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"prodinf1");
      x = gerepileupto(av,x);
    }
  }
  return gerepilecopy(av0,x);
}
GEN
prodinf0(GEN a, GEN code, long flag, long prec)
{
  switch(flag)
  {
    case 0: EXPR_WRAP(code, prodinf (EXPR_ARG, a, prec));
    case 1: EXPR_WRAP(code, prodinf1(EXPR_ARG, a, prec));
  }
  pari_err_FLAG("prodinf");
  return NULL; /* LCOV_EXCL_LINE */
}

GEN
prodeuler(void *E, GEN (*eval)(void *, GEN), GEN a, GEN b, long prec)
{
  pari_sp av, av0 = avma;
  GEN x = real_1(prec), prime;
  forprime_t T;

  av = avma;
  if (!forprime_init(&T, a,b)) return gc_const(av, x);

  av = avma;
  while ( (prime = forprime_next(&T)) )
  {
    x = gmul(x, eval(E, prime));
    if (gc_needed(av,1))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"prodeuler");
      x = gerepilecopy(av, x);
    }
  }
  return gerepilecopy(av0,x);
}
GEN
prodeuler0(GEN a, GEN b, GEN code, long prec)
{ EXPR_WRAP(code, prodeuler(EXPR_ARG, a, b, prec)); }
GEN
direuler0(GEN a, GEN b, GEN code, GEN c)
{ EXPR_WRAP(code, direuler(EXPR_ARG, a, b, c)); }

/********************************************************************/
/**                                                                **/
/**                       VECTORS & MATRICES                       **/
/**                                                                **/
/********************************************************************/

INLINE GEN
copyupto(GEN z, GEN t)
{
  if (is_universal_constant(z) || (z>(GEN)pari_mainstack->bot && z<=t))
    return z;
  else
    return gcopy(z);
}

GEN
vecexpr0(GEN vec, GEN code, GEN pred)
{
  switch(typ(vec))
  {
    case t_LIST:
    {
      if (list_typ(vec)==t_LIST_MAP)
        vec = mapdomain_shallow(vec);
      else
        vec = list_data(vec);
      if (!vec) return cgetg(1, t_VEC);
      break;
    }
    case t_VECSMALL:
      vec = vecsmall_to_vec(vec);
      break;
    case t_VEC: case t_COL: case t_MAT: break;
    default: pari_err_TYPE("[_|_<-_,_]",vec);
  }
  if (pred && code)
    EXPR_WRAP(code,vecselapply((void*)pred,&gp_evalbool,EXPR_ARGUPTO,vec))
  else if (code)
    EXPR_WRAP(code,vecapply(EXPR_ARGUPTO,vec))
  else
    EXPR_WRAP(pred,vecselect(EXPR_ARGBOOL,vec))
}

GEN
vecexpr1(GEN vec, GEN code, GEN pred)
{
  GEN v = vecexpr0(vec, code, pred);
  return lg(v) == 1? v: shallowconcat1(v);
}

GEN
vecteur(GEN nmax, GEN code)
{
  GEN y, c;
  long i, m = gtos(nmax);

  if (m < 0)  pari_err_DOMAIN("vector", "dimension", "<", gen_0, stoi(m));
  if (!code) return zerovec(m);
  c = cgetipos(3); /* left on stack */
  y = cgetg(m+1,t_VEC); push_lex(c, code);
  for (i=1; i<=m; i++)
  {
    c[2] = i;
    gel(y,i) = copyupto(closure_evalnobrk(code), y);
    set_lex(-1,c);
  }
  pop_lex(1); return y;
}

GEN
vecteursmall(GEN nmax, GEN code)
{
  pari_sp av;
  GEN y, c;
  long i, m = gtos(nmax);

  if (m < 0)  pari_err_DOMAIN("vectorsmall", "dimension", "<", gen_0, stoi(m));
  if (!code) return zero_zv(m);
  c = cgetipos(3); /* left on stack */
  y = cgetg(m+1,t_VECSMALL); push_lex(c,code);
  av = avma;
  for (i = 1; i <= m; i++)
  {
    c[2] = i;
    y[i] = gtos(closure_evalnobrk(code));
    set_avma(av);
    set_lex(-1,c);
  }
  pop_lex(1); return y;
}

GEN
vvecteur(GEN nmax, GEN n)
{
  GEN y = vecteur(nmax,n);
  settyp(y,t_COL); return y;
}

GEN
matrice(GEN nlig, GEN ncol, GEN code)
{
  GEN c1, c2, y;
  long i, m, n;

  n = gtos(nlig);
  m = ncol? gtos(ncol): n;
  if (m < 0)  pari_err_DOMAIN("matrix", "nbcols", "<", gen_0, stoi(m));
  if (n < 0)  pari_err_DOMAIN("matrix", "nbrows", "<", gen_0, stoi(n));
  if (!m) return cgetg(1,t_MAT);
  if (!code || !n) return zeromatcopy(n, m);
  c1 = cgetipos(3); push_lex(c1,code);
  c2 = cgetipos(3); push_lex(c2,NULL); /* c1,c2 left on stack */
  y = cgetg(m+1,t_MAT);
  for (i = 1; i <= m; i++)
  {
    GEN z = cgetg(n+1,t_COL);
    long j;
    c2[2] = i; gel(y,i) = z;
    for (j = 1; j <= n; j++)
    {
      c1[2] = j;
      gel(z,j) = copyupto(closure_evalnobrk(code), y);
      set_lex(-2,c1);
      set_lex(-1,c2);
    }
  }
  pop_lex(2); return y;
}

/********************************************************************/
/**                                                                **/
/**                         SUMMING SERIES                         **/
/**                                                                **/
/********************************************************************/
/* h = (2+2x)g'- g; g has t_INT coeffs */
static GEN
delt(GEN g, long n)
{
  GEN h = cgetg(n+3,t_POL);
  long k;
  h[1] = g[1];
  gel(h,2) = gel(g,2);
  for (k=1; k<n; k++)
    gel(h,k+2) = addii(mului(k+k+1,gel(g,k+2)), mului(k<<1,gel(g,k+1)));
  gel(h,n+2) = mului(n<<1, gel(g,n+1)); return h;
}

#ifdef _MSC_VER /* Bill Daly: work around a MSVC bug */
#pragma optimize("g",off)
#endif
/* P = polzagier(n,m)(-X), unnormalized (P(0) != 1) */
static GEN
polzag1(long n, long m)
{
  long d = n - m, i, k, d2, r, D;
  pari_sp av = avma;
  GEN g, T;

  if (d <= 0 || m < 0) return pol_0(0);
  d2 = d << 1; r = (m+1) >> 1, D = (d+1) >> 1;
  g = cgetg(d+2, t_POL);
  g[1] = evalsigne(1)|evalvarn(0);
  T = cgetg(d+1,t_VEC);
  /* T[k+1] = binomial(2d,2k+1), 0 <= k < d */
  gel(T,1) = utoipos(d2);
  for (k = 1; k < D; k++)
  {
    long k2 = k<<1;
    gel(T,k+1) = diviiexact(mulii(gel(T,k), muluu(d2-k2+1, d2-k2)),
                            muluu(k2,k2+1));
  }
  for (; k < d; k++) gel(T,k+1) = gel(T,d-k);
  gel(g,2) = gel(T,d); /* binomial(2d, 2(d-1)+1) */
  for (i = 1; i < d; i++)
  {
    pari_sp av2 = avma;
    GEN s, t = gel(T,d-i); /* binomial(2d, 2(d-1-i)+1) */
    s = t;
    for (k = d-i; k < d; k++)
    {
      long k2 = k<<1;
      t = diviiexact(mulii(t, muluu(d2-k2+1, d-k)), muluu(k2+1,k-(d-i)+1));
      s = addii(s, t);
    }
    /* g_i = sum_{d-1-i <= k < d}, binomial(2*d, 2*k+1)*binomial(k,d-1-i) */
    gel(g,i+2) = gerepileuptoint(av2, s);
  }
  /* sum_{0 <= i < d} g_i x^i * (x+x^2)^r */
  g = RgX_mulXn(gmul(g, gpowgs(deg1pol(gen_1,gen_1,0),r)), r);
  if (!odd(m)) g = delt(g, n);
  for (i = 1; i <= r; i++)
  {
    g = delt(ZX_deriv(g), n);
    if (gc_needed(av,4))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"polzag, i = %ld/%ld", i,r);
      g = gerepilecopy(av, g);
    }
  }
  return g;
}
GEN
polzag(long n, long m)
{
  pari_sp av = avma;
  GEN g = polzag1(n,m);
  if (lg(g) == 2) return g;
  g = ZX_z_unscale(polzag1(n,m), -1);
  return gerepileupto(av, RgX_Rg_div(g,gel(g,2)));
}

/*0.39322 > 1/log_2(3+sqrt(8))*/
static ulong
sumalt_N(long prec)
{ return (ulong)(0.39322*(prec2nbits(prec) + 7)); }

GEN
sumalt(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{
  ulong k, N;
  pari_sp av = avma, av2;
  GEN s, az, c, d;

  if (typ(a) != t_INT) pari_err_TYPE("sumalt",a);
  N = sumalt_N(prec);
  d = powru(addsr(3, sqrtr(utor(8,prec))), N);
  d = shiftr(addrr(d, invr(d)),-1);
  a = setloop(a);
  az = gen_m1; c = d;
  s = gen_0;
  av2 = avma;
  for (k=0; ; k++) /* k < N */
  {
    c = addir(az,c); s = gadd(s, gmul(c, eval(E, a)));
    if (k==N-1) break;
    az = diviuuexact(muluui((N-k)<<1,N+k,az), k+1, (k<<1)+1);
    a = incloop(a); /* in place! */
    if (gc_needed(av,4))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"sumalt, k = %ld/%ld", k,N-1);
      gerepileall(av2, 3, &az,&c,&s);
    }
  }
  return gerepileupto(av, gdiv(s,d));
}

GEN
sumalt2(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{
  long k, N;
  pari_sp av = avma, av2;
  GEN s, dn, pol;

  if (typ(a) != t_INT) pari_err_TYPE("sumalt",a);
  N = (long)(0.307073*(prec2nbits(prec) + 5)); /*0.307073 > 1/log_2(\beta_B)*/
  pol = ZX_div_by_X_1(polzag1(N,N>>1), &dn);
  a = setloop(a);
  N = degpol(pol);
  s = gen_0;
  av2 = avma;
  for (k=0; k<=N; k++)
  {
    GEN t = itor(gel(pol,k+2), prec+EXTRAPREC64);
    s = gadd(s, gmul(t, eval(E, a)));
    if (k == N) break;
    a = incloop(a); /* in place! */
    if (gc_needed(av,4))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"sumalt2, k = %ld/%ld", k,N-1);
      s = gerepileupto(av2, s);
    }
  }
  return gerepileupto(av, gdiv(s,dn));
}

GEN
sumalt0(GEN a, GEN code, long flag, long prec)
{
  switch(flag)
  {
    case 0: EXPR_WRAP(code, sumalt (EXPR_ARG,a,prec));
    case 1: EXPR_WRAP(code, sumalt2(EXPR_ARG,a,prec));
    default: pari_err_FLAG("sumalt");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/* For k > 0, set S[k*2^i] <- g(k*2^i), k*2^i <= N = #S.
 * Only needed with k odd (but also works for g even). */
static void
binsum(GEN S, ulong k, void *E, GEN (*f)(void *, GEN), GEN a,
        long G, long prec)
{
  long e, i, N = lg(S)-1, l = expu(N / k); /* k 2^l <= N < k 2^(l+1) */
  pari_sp av = avma;
  GEN t = real_0(prec); /* unused unless f(a + k <<l) = 0 */

  G -= l;
  if (!signe(a)) a = NULL;
  for (e = 0;; e++)
  { /* compute g(k 2^l) with absolute error ~ 2^(G-l) */
    GEN u, r = shifti(utoipos(k), l+e);
    if (a) r = addii(r, a);
    u = gtofp(f(E, r), prec);
    if (typ(u) != t_REAL) pari_err_TYPE("sumpos",u);
    if (!signe(u)) break;
    if (!e)
      t = u;
    else {
      shiftr_inplace(u, e);
      t = addrr(t,u); if (expo(u) < G) break;
      if ((e & 0x1ff) == 0) t = gerepileuptoleaf(av, t);
    }
  }
  gel(S, k << l) = t = gerepileuptoleaf(av, t);
  /* g(j) = 2g(2j) + f(a+j) for all j > 0 */
  for(i = l-1; i >= 0; i--)
  { /* t ~ g(2 * k*2^i) with error ~ 2^(G-i-1) */
    GEN u;
    av = avma; u = gtofp(f(E, a? addiu(a, k << i): utoipos(k << i)), prec);
    if (typ(u) != t_REAL) pari_err_TYPE("sumpos",u);
    t = addrr(gtofp(u,prec), mpshift(t,1)); /* ~ g(k*2^i) */
    gel(S, k << i) = t = gerepileuptoleaf(av, t);
  }
}
/* For k > 0, let g(k) := \sum_{e >= 0} 2^e f(a + k*2^e).
 * Return [g(k), 1 <= k <= N] */
static GEN
sumpos_init(void *E, GEN (*f)(void *, GEN), GEN a, long N, long prec)
{
  GEN S = cgetg(N+1,t_VEC);
  long k, G = -prec2nbits(prec) - 5;
  for (k=1; k<=N; k+=2) binsum(S,k, E,f, a,G,prec);
  return S;
}

GEN
sumpos(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{
  ulong k, N;
  pari_sp av = avma;
  GEN s, az, c, d, S;

  if (typ(a) != t_INT) pari_err_TYPE("sumpos",a);
  a = subiu(a,1);
  N = sumalt_N(prec);
  if (odd(N)) N++; /* extra precision for free */
  d = powru(addsr(3, sqrtr(utor(8,prec))), N);
  d = shiftr(addrr(d, invr(d)),-1);
  az = gen_m1; c = d;

  S = sumpos_init(E, eval, a, N, prec);
  s = gen_0;
  for (k=0; k<N; k++)
  {
    GEN t;
    c = addir(az,c);
    t = mulrr(gel(S,k+1), c);
    s = odd(k)? mpsub(s, t): mpadd(s, t);
    if (k == N-1) break;
    az = diviuuexact(muluui((N-k)<<1,N+k,az), k+1, (k<<1)+1);
  }
  return gerepileupto(av, gdiv(s,d));
}

GEN
sumpos2(void *E, GEN (*eval)(void *, GEN), GEN a, long prec)
{
  ulong k, N;
  pari_sp av = avma;
  GEN s, pol, dn, S;

  if (typ(a) != t_INT) pari_err_TYPE("sumpos2",a);
  a = subiu(a,1);
  N = (ulong)(0.31*(prec2nbits(prec) + 5));

  if (odd(N)) N++; /* extra precision for free */
  S = sumpos_init(E, eval, a, N, prec);
  pol = ZX_div_by_X_1(polzag1(N,N>>1), &dn);
  s = gen_0;
  for (k=0; k<N; k++)
  {
    GEN t = mulri(gel(S,k+1), gel(pol,k+2));
    s = odd(k)? mpsub(s,t): mpadd(s,t);
  }
  return gerepileupto(av, gdiv(s,dn));
}

GEN
sumpos0(GEN a, GEN code, long flag, long prec)
{
  switch(flag)
  {
    case 0: EXPR_WRAP(code, sumpos (EXPR_ARG,a,prec));
    case 1: EXPR_WRAP(code, sumpos2(EXPR_ARG,a,prec));
    default: pari_err_FLAG("sumpos");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/********************************************************************/
/**                                                                **/
/**            SEARCH FOR REAL ZEROS of an expression              **/
/**                                                                **/
/********************************************************************/
/* Brent's method, [a,b] bracketing interval */
GEN
zbrent(void *E, GEN (*eval)(void *, GEN), GEN a, GEN b, long prec)
{
  long sig, iter, itmax, bit, bit0;
  pari_sp av = avma;
  GEN c, d, e, fa, fb, fc;

  if (typ(a) == t_INFINITY && typ(b) != t_INFINITY) swap(a,b);
  if (typ(a) == t_INFINITY && typ(b) == t_INFINITY)
  {
    long s = gsigne(eval(E, real_0(prec))), r = 0;
    if (gidentical(gel(a,1), gel(b,1)))
      pari_err_DOMAIN("solve", "a and b", "=", a, mkvec2(a, b));
    a = real_m1(prec); /* domain = R */
    b = real_1(prec);
    for(;;)
    {
      fa = eval(E, a);
      fb = eval(E, b);
      if (gsigne(fa) != s)
      {
        if (r) b[1] = evalsigne(-1) | _evalexpo(r-1); else b = real_0(prec);
        break;
      }
      if (gsigne(fb) != s)
      {
        if (r) a[1] = evalsigne(1) | _evalexpo(r-1); else a = real_0(prec);
        break;
      }
      r++; setexpo(a, r); setexpo(b, r);
    }
    c = b;
    goto SOLVE;
  }
  if (typ(b) == t_INFINITY)
  { /* a real, b == [+-]oo */
    long s, r, minf = inf_get_sign(b) < 0;
    GEN inc;
    if (typ(a) != t_REAL || realprec(a) < prec) a = gtofp(a, prec);
    fa = eval(E, a);
    s = gsigne(fa);
    inc = minf ? real_m1(prec) : real_1(prec);
    r = gsigne(a) ? expo(a) : 0;
    for(;;)
    {
      setexpo(inc, r);
      b = addrr(a, inc); fb = eval(E, b);
      if (gsigne(fb) != s) break;
      a = b; fa = fb; r++;
    }
    if (minf) { c = a; swap(a, b); swap(fa, fb);} else c = b;
    goto SOLVE;
  }
  if (typ(a) != t_REAL || realprec(a) < prec) a = gtofp(a, prec);
  if (typ(b) != t_REAL || realprec(b) < prec) b = gtofp(b, prec);
  sig = cmprr(b, a);
  if (!sig) return gerepileupto(av, a);
  if (sig < 0) { c = a; swap(a, b); } else c = b;
  fa = eval(E, a);
  fb = eval(E, b);
  if (gsigne(fa)*gsigne(fb) > 0)
    pari_err_DOMAIN("solve", "f(a)f(b)", ">", gen_0, mkvec2(fa, fb));
SOLVE:
  bit0 = -prec2nbits(prec); bit = 3+bit0; itmax = 1 - 2*bit0;
  fc = fb; e = d = NULL;
  for (iter = 1; iter <= itmax; ++iter)
  { /* b = current best guess, a and c auxiliary points */
    long bit2, exb;
    GEN m;
    if (gsigne(fb)*gsigne(fc) > 0) { c = a; fc = fa; e = d = subrr(b, a); }
    if (gcmp(gabs(fc, 0), gabs(fb, 0)) < 0)
    { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    m = subrr(c, b); shiftr_inplace(m, -1);
    exb = expo(b);
    if (bit < exb)
    {
      bit2 = bit + exb - 1;
      if (expo(m) <= exb + bit0) break; /*SUCCESS*/
    }
    else
    { /* b ~ 0 */
      bit2 = 2*bit - 1;
      if (expo(m) <= bit2) break; /*SUCCESS*/
    }
    if (gequal0(fb)) break; /*SUCCESS*/

    if (expo(e) > bit2 && gexpo(fa) > gexpo(fb))
    { /* interpolation, m != 0, |f(c)| >= |f(b)|, f(b)f(c) < 0 */
      GEN min1, min2, p, q, s = gdiv(fb, fa);
      if (a == c || equalrr(a,c))
      {
        p = gmul2n(gmul(m, s), 1);
        q = gsubsg(1, s);
      }
      else
      {
        GEN r = gdiv(fb, fc);
        q = gdiv(fa, fc);
        p = gmul2n(gmul(gsub(q, r), gmul(m, q)), 1);
        p = gmul(s, gsub(p, gmul(subrr(b, a), gsubgs(r, 1))));
        q = gmul(gmul(gsubgs(q, 1), gsubgs(r, 1)), gsubgs(s, 1));
      }
      if (gsigne(p) > 0) q = gneg_i(q); else p = gneg_i(p);
      min1 = gsub(gmulsg(3, gmul(m,q)), gmul2n(gabs(q,0), bit2));
      min2 = gabs(gmul(e, q), 0);
      if (gcmp(gmul2n(p, 1), gmin_shallow(min1, min2)) < 0)
        { e = d; d = gdiv(p, q); } /* interpolation OK */
      else
        e = d = m; /* failed, use bisection */
    }
    else e = d = m; /* bound decreasing too slowly, use bisection */
    a = b; fa = fb;
    if (d == m) { b = addrr(c, b); shiftr_inplace(b,-1); }
    else if (gexpo(d) > bit2) b = gadd(b, d);
    else if (gsigne(m) > 0) b = addrr(b, real2n(bit2, LOWDEFAULTPREC));
    else                    b = subrr(b, real2n(bit2, LOWDEFAULTPREC));
    if (equalrr(a, b)) break;
    if (realprec(b) < prec) b = rtor(b, prec);
    fb = eval(E, b);
  }
  if (iter > itmax) pari_err_IMPL("solve recovery [too many iterations]");
  return gerepileuptoleaf(av, rcopy(b));
}

GEN
zbrent0(GEN a, GEN b, GEN code, long prec)
{ EXPR_WRAP(code, zbrent(EXPR_ARG, a, b, prec)); }

/* Find zeros of a function in the real interval [a,b] by interval splitting */
GEN
solvestep(void *E, GEN (*f)(void *,GEN), GEN a, GEN b, GEN step, long flag, long prec)
{
  const long ITMAX = 10;
  pari_sp av = avma;
  GEN fa, a0, b0;
  long sa0, it, bit = prec2nbits(prec) / 2, ct = 0, s = gcmp(a,b);

  if (!s) return gequal0(f(E, a)) ? gcopy(mkvec(a)): cgetg(1,t_VEC);
  if (s > 0) swap(a, b);
  if (flag&4)
  {
    if (gcmpgs(step,1)<=0) pari_err_DOMAIN("solvestep","step","<=",gen_1,step);
    if (gsigne(a) <= 0) pari_err_DOMAIN("solvestep","a","<=",gen_0,a);
  }
  else if (gsigne(step) <= 0)
    pari_err_DOMAIN("solvestep","step","<=",gen_0,step);
  a0 = a = gtofp(a, prec); fa = f(E, a);
  b0 = b = gtofp(b, prec); step = gtofp(step, prec);
  sa0 = gsigne(fa);
  if (gexpo(fa) < -bit) sa0 = 0;
  for (it = 0; it < ITMAX; it++)
  {
    pari_sp av2 = avma;
    GEN v = cgetg(1, t_VEC);
    long sa = sa0;
    a = a0; b = b0;
    while (gcmp(a,b) < 0)
    {
      GEN fc, c = (flag&4)? gmul(a, step): gadd(a, step);
      long sc;
      if (gcmp(c,b) > 0) c = b;
      fc = f(E, c); sc = gsigne(fc);
      if (gexpo(fc) < -bit) sc = 0;
      if (!sc || sa*sc < 0)
      {
        GEN z = sc? zbrent(E, f, a, c, prec): c;
        long e;
        (void)grndtoi(z, &e);
        if (e <= -bit) ct = 1;
        if ((flag&1) && ((!(flag&8)) || ct)) return gerepileupto(av, z);
        v = shallowconcat(v, z);
      }
      a = c; fa = fc; sa = sc;
      if (gc_needed(av2,1))
      {
        if (DEBUGMEM>1) pari_warn(warnmem,"solvestep");
        gerepileall(av2, 4, &a ,&fa, &v, &step);
      }
    }
    if ((!(flag&2) || lg(v) > 1) && (!(flag&8) || ct))
      return gerepilecopy(av, v);
    step = (flag&4)? sqrtnr(step,4): gmul2n(step, -2);
    gerepileall(av2, 2, &fa, &step);
  }
  pari_err_IMPL("solvestep recovery [too many iterations]");
  return NULL;/*LCOV_EXCL_LINE*/
}

GEN
solvestep0(GEN a, GEN b, GEN step, GEN code, long flag, long prec)
{ EXPR_WRAP(code, solvestep(EXPR_ARG, a,b, step, flag, prec)); }

/********************************************************************/
/**                     Numerical derivation                       **/
/********************************************************************/

struct deriv_data
{
  GEN code;
  GEN args;
  GEN def;
};

static GEN deriv_eval(void *E, GEN x, long prec)
{
 struct deriv_data *data=(struct deriv_data *)E;
 gel(data->args,1)=x;
 uel(data->def,1)=1;
 return closure_callgenvecdefprec(data->code, data->args, data->def, prec);
}

/* Rationale: (f(2^-e) - f(-2^-e) + O(2^-b)) / (2 * 2^-e) = f'(0) + O(2^-2e)
 * since 2nd derivatives cancel.
 *   prec(LHS) = b - e
 *   prec(RHS) = 2e, equal when  b = 3e = 3/2 b0 (b0 = required final bitprec)
 *
 * For f'(x), x far from 0: prec(LHS) = b - e - expo(x)
 * --> pr = 3/2 b0 + expo(x) */
GEN
derivnum(void *E, GEN (*eval)(void *, GEN, long), GEN x, long prec)
{
  long newprec, e, ex = gexpo(x), p = precision(x);
  long b0 = prec2nbits(p? p: prec), b = (long)ceil(b0 * 1.5 + maxss(0,ex));
  GEN eps, u, v, y;
  pari_sp av = avma;
  newprec = nbits2prec(b + BITS_IN_LONG);
  switch(typ(x))
  {
    case t_REAL:
    case t_COMPLEX:
      x = gprec_w(x, newprec);
  }
  e = b0/2; /* 1/2 required prec (in sig. bits) */
  b -= e; /* >= b0 */
  eps = real2n(-e, ex < -e? newprec: nbits2prec(b));
  u = eval(E, gsub(x, eps), newprec);
  v = eval(E, gadd(x, eps), newprec);
  y = gmul2n(gsub(v,u), e-1);
  return gerepilecopy(av, gprec_wtrunc(y, nbits2prec(b0)));
}

/* Fornberg interpolation algorithm for finite differences coefficients
* using 2N+1 equidistant grid points around 0 [ assume 2N even >= M ].
* Compute \delta[m]_{N,i} for all derivation orders m = 0..M such that
*   h^m * f^{(m)}(0) = \sum_{i = 0}^n delta[m]_{N,i}  f(a_i) + O(h^{N-m+1}),
* for step size h.
* Return a = [0,-1,1...,-N,N] and vector of vectors d: d[m+1][i+1]
* = w'(a_i) delta[m]_{2N,i}, i = 0..2N */
static void
FD(long M, long N2, GEN *pd, GEN *pa)
{
  GEN d, a, b, W, F;
  long N = N2>>1, m, i;

  F = cgetg(N2+2, t_VEC);
  a = cgetg(N2+2, t_VEC);
  b = cgetg(N+1, t_VEC);
  gel(a,1) = gen_0;
  for (i = 1; i <= N; i++)
  {
    gel(a,2*i)   = utoineg(i);
    gel(a,2*i+1) = utoipos(i);
    gel(b,i) = sqru(i);
  }
  /* w = \prod (X - a[i]) = x W(x^2) */
  W = roots_to_pol(b, 0);
  gel(F,1) = RgX_inflate(W,2);
  for (i = 1; i <= N; i++)
  {
    pari_sp av = avma;
    GEN r, U, S;
    U = RgX_inflate(RgX_div_by_X_x(W, gel(b,i), &r), 2);
    U = RgXn_red_shallow(U, M); /* higher terms not needed */
    U = RgX_shift_shallow(U,1); /* w(X) / (X^2-a[i]^2) mod X^(M+1) */
    S = ZX_sub(RgX_shift_shallow(U,1),
               ZX_Z_mul(U, gel(a,2*i+1)));
    S = gerepileupto(av, S);
    gel(F,2*i)   = S;
    gel(F,2*i+1) = ZX_z_unscale(S, -1);
  }
  /* F[i] = w(X) / (X-a[i]) + O(X^(M+1)) in Z[X] */
  d = cgetg(M+2, t_VEC);
  for (m = 0; m <= M; m++)
  {
    GEN v = cgetg(N2+2, t_VEC); /* coeff(F[i],X^m) */
    for (i = 0; i <= N2; i++) gel(v, i+1) = gmael(F, i+1, m+2);
    gel(d,m+1) = v;
  }
  *pd = d;
  *pa = a;
}

static void
chk_ord(long m)
{
  if (m < 0)
    pari_err_DOMAIN("derivnumk", "derivation order", "<", gen_0, stoi(m));
}
/* m! / N! for m in ind; vecmax(ind) <= N. Result not a GEN if ind contains 0. */
static GEN
vfact(GEN ind, long N, long prec)
{
  GEN v, iN;
  long i, l;
  ind = vecsmall_uniq(ind); chk_ord(ind[1]); l = lg(ind);
  iN = invr(itor(mulu_interval(ind[1] + 1, N), prec));
  v = const_vec(ind[l-1], NULL); gel(v, ind[1]) = iN;
  for (i = 2; i < l; i++)
    gel(v, ind[i]) = iN = mulri(iN, mulu_interval(ind[i-1] + 1, ind[i]));
  return v;
}

static GEN
chk_ind(GEN ind, long *M)
{
  *M = 0;
  switch(typ(ind))
  {
    case t_INT: ind = mkvecsmall(itos(ind)); break;
    case t_VECSMALL:
      if (lg(ind) == 1) return NULL;
      break;
    case t_VEC: case t_COL:
      if (lg(ind) == 1) return NULL;
      if (RgV_is_ZV(ind)) { ind = ZV_to_zv(ind); break; }
      /* fall through */
    default:
      pari_err_TYPE("derivnum", ind);
      return NULL; /*LCOV_EXCL_LINE*/
  }
  *M = vecsmall_max(ind); chk_ord(*M); return ind;
}
GEN
derivnumk(void *E, GEN (*eval)(void *, GEN, long), GEN x, GEN ind0, long prec)
{
  GEN A, C, D, DM, T, X, F, v, ind, t;
  long M, N, N2, fpr, p, i, pr, l, lA, e, ex, emin, emax, newprec;
  pari_sp av = avma;
  int allodd = 1;

  ind = chk_ind(ind0, &M); if (!ind) return cgetg(1, t_VEC);
  l = lg(ind); F = cgetg(l, t_VEC);
  if (!M) /* silly degenerate case */
  {
    X = eval(E, x, prec);
    for (i = 1; i < l; i++) { chk_ord(ind[i]); gel(F,i) = X; }
    if (typ(ind0) == t_INT) F = gel(F,1);
    return gerepilecopy(av, F);
  }
  N2 = 3*M - 1; if (odd(N2)) N2++;
  N = N2 >> 1;
  FD(M, N2, &D,&A); /* optimal if 'eval' uses quadratic time */
  C = vecbinomial(N2); DM = gel(D,M);
  T = cgetg(N2+2, t_VEC);
  /* (2N)! / w'(i) = (2N)! / w'(-i) = (-1)^(N-i) binom(2*N, N-i) */
  t = gel(C, N+1);
  gel(T,1) = odd(N)? negi(t): t;
  for (i = 1; i <= N; i++)
  {
    t = gel(C, N-i+1);
    gel(T,2*i) = gel(T,2*i+1) = odd(N-i)? negi(t): t;
  }
  N = N2 >> 1; emin = LONG_MAX; emax = 0;
  for (i = 1; i <= N; i++)
  {
    e = expi(gel(DM,i)) + expi(gel(T,i));
    if (e < 0) continue; /* 0 */
    if (e < emin) emin = e;
    else if (e > emax) emax = e;
  }

  p = precision(x);
  fpr = p ? prec2nbits(p): prec2nbits(prec);
  e = (fpr + 3*M*log2((double)M)) / (2*M);
  ex = gexpo(x);
  if (ex < 0) ex = 0; /* near 0 */
  pr = (long)ceil(fpr + e * M); /* ~ 3fpr/2 */
  newprec = nbits2prec(pr + (emax - emin) + ex + BITS_IN_LONG);
  switch(typ(x))
  {
    case t_REAL:
    case t_COMPLEX:
      x = gprec_w(x, newprec);
  }
  lA = lg(A); X = cgetg(lA, t_VEC);
  for (i = 1; i < l; i++)
    if (!odd(ind[i])) { allodd = 0; break; }
  /* if only odd derivation orders, the value at 0 (A[1]) is not needed */
  gel(X, 1) = gen_0;
  for (i = allodd? 2: 1; i < lA; i++)
  {
    GEN t = eval(E, gadd(x, gmul2n(gel(A,i), -e)), newprec);
    t = gmul(t, gel(T,i));
    if (!gprecision(t))
      t = is_scalar_t(typ(t))? gtofp(t, newprec): gmul(t, real_1(newprec));
    gel(X,i) = t;
  }

  v = vfact(ind, N2, nbits2prec(fpr + 32));
  for (i = 1; i < l; i++)
  {
    long m = ind[i];
    GEN t = RgV_dotproduct(gel(D,m+1), X);
    gel(F,i) = gmul(t, gmul2n(gel(v, m), e*m));
  }
  if (typ(ind0) == t_INT) F = gel(F,1);
  return gerepilecopy(av, F);
}
/* v(t') */
static long
rfrac_val_deriv(GEN t)
{
  long v = varn(gel(t,2));
  return gvaluation(deriv(t, v), pol_x(v));
}

GEN
derivfunk(void *E, GEN (*eval)(void *, GEN, long), GEN x, GEN ind0, long prec)
{
  pari_sp av;
  GEN ind, xp, ixp, F, G;
  long i, l, vx, M;
  if (!ind0) return derivfun(E, eval, x, prec);
  switch(typ(x))
  {
  case t_REAL: case t_INT: case t_FRAC: case t_COMPLEX:
    return derivnumk(E,eval, x, ind0, prec);
  case t_POL:
    ind = chk_ind(ind0,&M); if (!ind) return cgetg(1,t_VEC);
    xp = RgX_deriv(x);
    x = RgX_to_ser(x, precdl+2 + M * (1+RgX_val(xp)));
    break;
  case t_RFRAC:
    ind = chk_ind(ind0,&M); if (!ind) return cgetg(1,t_VEC);
    x = rfrac_to_ser_i(x, precdl+2 + M * (1+rfrac_val_deriv(x)));
    xp = derivser(x);
    break;
  case t_SER:
    ind = chk_ind(ind0,&M); if (!ind) return cgetg(1,t_VEC);
    xp = derivser(x);
    break;
  default: pari_err_TYPE("numerical derivation",x);
    return NULL; /*LCOV_EXCL_LINE*/
  }
  av = avma; vx = varn(x);
  ixp = M? ginv(xp): NULL;
  F = cgetg(M+2, t_VEC);
  gel(F,1) = eval(E, x, prec);
  for (i = 1; i <= M; i++) gel(F,i+1) = gmul(deriv(gel(F,i),vx), ixp);
  l = lg(ind); G = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    long m = ind[i]; chk_ord(m);
    gel(G,i) = gel(F,m+1);
  }
  if (typ(ind0) == t_INT) G = gel(G,1);
  return gerepilecopy(av, G);
}

GEN
derivfun(void *E, GEN (*eval)(void *, GEN, long), GEN x, long prec)
{
  pari_sp av = avma;
  GEN xp;
  long vx;
  switch(typ(x))
  {
  case t_REAL: case t_INT: case t_FRAC: case t_COMPLEX:
    return derivnum(E,eval, x, prec);
  case t_POL:
    xp = RgX_deriv(x);
    x = RgX_to_ser(x, precdl+2+ (1 + RgX_val(xp)));
    break;
  case t_RFRAC:
    x = rfrac_to_ser_i(x, precdl+2+ (1 + rfrac_val_deriv(x)));
    /* fall through */
  case t_SER:
    xp = derivser(x);
    break;
  default: pari_err_TYPE("formal derivation",x);
    return NULL; /*LCOV_EXCL_LINE*/
  }
  vx = varn(x);
  return gerepileupto(av, gdiv(deriv(eval(E, x, prec),vx), xp));
}

GEN
laurentseries(void *E, GEN (*f)(void*,GEN x, long), long M, long v, long prec)
{
  pari_sp av = avma;
  long d;

  if (v < 0) v = 0;
  d = maxss(M+1,1);
  for (;;)
  {
    long i, dr, vr;
    GEN s;
    s = cgetg(d+2, t_SER); s[1] = evalsigne(1) | evalvalser(1) | evalvarn(v);
    gel(s, 2) = gen_1; for (i = 3; i <= d+1; i++) gel(s, i) = gen_0;
    s = f(E, s, prec);
    if (typ(s) != t_SER || varn(s) != v) pari_err_TYPE("laurentseries", s);
    vr = valser(s);
    if (M < vr) { set_avma(av); return zeroser(v, M); }
    dr = lg(s) + vr - 3 - M;
    if (dr >= 0) return gerepileupto(av, s);
    set_avma(av); d -= dr;
  }
}
static GEN
_evalclosprec(void *E, GEN x, long prec)
{
  GEN s;
  push_localprec(prec); s = closure_callgen1((GEN)E, x);
  pop_localprec(); return s;
}
#define CLOS_ARGPREC __E, &_evalclosprec
GEN
laurentseries0(GEN f, long M, long v, long prec)
{
  if (typ(f) != t_CLOSURE || closure_arity(f) != 1 || closure_is_variadic(f))
    pari_err_TYPE("laurentseries",f);
  EXPR_WRAP(f, laurentseries(CLOS_ARGPREC,M,v,prec));
}

GEN
derivnum0(GEN a, GEN code, GEN ind, long prec)
{ EXPR_WRAP(code, derivfunk(EXPR_ARGPREC,a,ind,prec)); }

GEN
derivfun0(GEN args, GEN def, GEN code, long k, long prec)
{
  pari_sp av = avma;
  struct deriv_data E;
  GEN z;
  E.code=code; E.args=args; E.def=def;
  z = gel(derivfunk((void*)&E, deriv_eval, gel(args,1), mkvecs(k), prec),1);
  return gerepilecopy(av, z);
}

/********************************************************************/
/**                   Numerical extrapolation                      **/
/********************************************************************/
/* [u(n), u <= N] */
static GEN
get_u(void *E, GEN (*f)(void *, GEN, long), long N, long prec)
{
  long n;
  GEN u;
  if (f)
  {
    GEN v = f(E, utoipos(N), prec);
    u = cgetg(N+1, t_VEC);
    if (typ(v) != t_VEC || lg(v) != N+1) { gel(u,N) = v; v = NULL; }
    else
    {
      GEN w = f(E, gen_1, LOWDEFAULTPREC);
      if (typ(w) != t_VEC || lg(w) != 2) { gel(u,N) = v; v = NULL; }
    }
    if (v) u = v;
    else
      for (n = 1; n < N; n++) gel(u,n) = f(E, utoipos(n), prec);
  }
  else
  {
    GEN v = (GEN)E;
    long t = lg(v)-1;
    if (t < N) pari_err_COMPONENT("limitnum","<",stoi(N), stoi(t));
    u = vecslice(v, 1, N);
  }
  for (n = 1; n <= N; n++)
  {
    GEN un = gel(u,n);
    if (is_rational_t(typ(un))) gel(u,n) = gtofp(un, prec);
  }
  return u;
}

struct limit
{
  long prec; /* working accuracy */
  long N; /* number of terms */
  GEN na; /* [n^alpha, n <= N] */
  GEN coef; /* or NULL (alpha != 1) */
};

static GEN
_gi(void *E, GEN x)
{
  GEN A = (GEN)E, y = gsubgs(x, 1);
  if (gequal0(y)) return A;
  return gdiv(gsubgs(gpow(x, A, LOWDEFAULTPREC), 1), y);
}
static GEN
_g(void *E, GEN x)
{
  GEN D = (GEN)E, A = gel(D,1), T = gel(D,2);
  const long prec = LOWDEFAULTPREC;
  return gadd(glog(x,prec), intnum((void*)A, _gi, gen_0, gaddgs(x,1), T, prec));
}

/* solve log(b) + int_0^{b+1} (x^(1/a)-1) / (x-1) dx = 0, b in [0,1]
 * return -log_2(b), rounded up */
static double
get_accu(GEN a)
{
  pari_sp av = avma;
  const long prec = LOWDEFAULTPREC;
  const double We2 = 1.844434455794; /* (W(1/e) + 1) / log(2) */
  GEN b, T;
  if (!a) return We2;
  if (typ(a) == t_INT) switch(itos_or_0(a))
  {
    case 1: return We2;
    case 2: return 1.186955309668;
    case 3: return 0.883182331990;
  }
  else if (typ(a) == t_FRAC && equali1(gel(a,1))) switch(itos_or_0(gel(a,2)))
  {
    case 2: return 2.644090500290;
    case 3: return 3.157759214459;
    case 4: return 3.536383237500;
  }
  T = intnuminit(gen_0, gen_1, 0, prec);
  b = zbrent((void*)mkvec2(ginv(a), T), &_g, dbltor(1E-5), gen_1, prec);
  return gc_double(av, -dbllog2r(b));
}

static double
get_c(GEN a)
{
  double A = a? gtodouble(a): 1.0;
  if (A <= 0) pari_err_DOMAIN("limitnum","alpha","<=",gen_0, a);
  if (A >= 2) return 0.2270;
  if (A >= 1) return 0.3318;
  if (A >= 0.5) return 0.6212;
  if (A >= 0.3333) return 1.2;
  return 3; /* only tested for A >= 0.25 */
}
static void
limit_Nprec(struct limit *L, GEN alpha, long prec)
{
  long bit = prec2nbits(prec);
  L->N = ceil(get_c(alpha) * bit);
  L->prec = nbits2prec(bit + (long)ceil(get_accu(alpha) * L->N));
}
/* solve x - a log(x) = b; a, b >= 3 */
static double
solvedivlog(double a, double b) { return dbllemma526(a,1,1,b); }

/* #u > 1, prod_{j != i} u[i] - u[j] */
static GEN
proddiff(GEN u, long i)
{
  pari_sp av = avma;
  long l = lg(u), j;
  GEN p = NULL;
  if (i == 1)
  {
    p = gsub(gel(u,1), gel(u,2));
    for (j = 3; j < l; j++)
      p = gmul(p, gsub(gel(u,i), gel(u,j)));
  }
  else
  {
    p = gsub(gel(u,i), gel(u,1));
    for (j = 2; j < l; j++)
      if (j != i) p = gmul(p, gsub(gel(u,i), gel(u,j)));
  }
  return gerepileupto(av, p);
}
static GEN
vecpows(GEN u, long N)
{
  long i, l;
  GEN v = cgetg_copy(u, &l);
  for (i = 1; i < l; i++) gel(v,i) = gpowgs(gel(u,i), N);
  return v;
}

static void
limit_init(struct limit *L, GEN alpha, int asymp)
{
  long n, N = L->N, a = 0;
  GEN c, v, T = NULL;

  if (!alpha) a = 1;
  else if (typ(alpha) == t_INT)
  {
    a = itos_or_0(alpha);
    if (a > 2) a = 0;
  }
  else if (typ(alpha) == t_FRAC)
  {
    long na = itos_or_0(gel(alpha,1)), da = itos_or_0(gel(alpha,2));
    if (da && na && da <= 4 && na <= 4)
    { /* don't bother with other cases */
      long e = (N-1) % da, k = (N-1) / da;
      if (e) { N += da - e; k++; } /* N = 1 (mod d) => simpler ^ (n/d)(N-1) */
      L->N = N;
      T = vecpowuu(N, na * k);
    }
  }
  L->coef = v = cgetg(N+1, t_VEC);
  if (!a)
  {
    long prec2 = gprecision(alpha);
    GEN u;
    if (prec2 && prec2 < L->prec) alpha = gprec_w(alpha, L->prec);
    L->na = u = vecpowug(N, alpha, L->prec);
    if (!T) T = vecpows(u, N-1);
    for (n = 1; n <= N; n++) gel(v,n) = gdiv(gel(T,n), proddiff(u,n));
    return;
  }
  L->na = asymp? vecpowuu(N, a): NULL;
  c = mpfactr(N-1, L->prec);
  if (a == 1)
  {
    c = invr(c);
    gel(v,1) = c; if (!odd(N)) togglesign(c);
    for (n = 2; n <= N; n++) gel(v,n) = divru(mulrs(gel(v,n-1), n-1-N), n);
  }
  else
  { /* a = 2 */
    c = invr(mulru(sqrr(c), (N*(N+1)) >> 1));
    gel(v,1) = c; if (!odd(N)) togglesign(c);
    for (n = 2; n <= N; n++) gel(v,n) = divru(mulrs(gel(v,n-1), n-1-N), N+n);
  }
  T = vecpowuu(N, a*N);
  for (n = 2; n <= N; n++) gel(v,n) = mulri(gel(v,n), gel(T,n));
}

/* Zagier/Lagrange extrapolation */
static GEN
limitnum_i(struct limit *L, GEN u, long prec)
{ return gprec_w(RgV_dotproduct(u,L->coef), prec); }
GEN
limitnum(void *E, GEN (*f)(void *, GEN, long), GEN alpha, long prec)
{
  pari_sp av = avma;
  struct limit L;
  GEN u;
  limit_Nprec(&L, alpha, prec);
  limit_init(&L, alpha, 0);
  u = get_u(E, f, L.N, L.prec);
  return gerepilecopy(av, limitnum_i(&L, u, prec));
}
typedef GEN (*LIMIT_FUN)(void*,GEN,long);
static LIMIT_FUN get_fun(GEN u, const char *s)
{
  switch(typ(u))
  {
    case t_COL: case t_VEC: break;
    case t_CLOSURE: return gp_callprec;
    default: pari_err_TYPE(s, u);
  }
  return NULL;
}
GEN
limitnum0(GEN u, GEN alpha, long prec)
{ return limitnum((void*)u,get_fun(u, "limitnum"), alpha, prec); }

GEN
asympnum(void *E, GEN (*f)(void *, GEN, long), GEN alpha, long prec)
{
  const long MAX = 100;
  pari_sp av = avma;
  GEN u, A = cgetg(MAX+1, t_VEC);
  long i, B = prec2nbits(prec);
  double LB = 0.9*expu(B); /* 0.9 and 0.95 below are heuristic */
  struct limit L;
  limit_Nprec(&L, alpha, prec);
  if (alpha) LB *= gtodouble(alpha);
  limit_init(&L, alpha, 1);
  u = get_u(E, f, L.N, L.prec);
  for(i = 1; i <= MAX; i++)
  {
    GEN a, v, q, s = limitnum_i(&L, u, prec);
    long n;
    /* NOT bestappr: lindep properly ignores the lower bits */
    v = lindep_bit(mkvec2(gen_1, s), maxss((long)(0.95*floor(B - i*LB)), 32));
    if (lg(v) == 1) break;
    q = gel(v,2); if (!signe(q)) break;
    a = gdiv(negi(gel(v,1)), q);
    s = gsub(s, a);
    /* |s|q^2 > eps */
    if (!gequal0(s) && gexpo(s) + 2*expi(q) > -17) break;
    gel(A,i) = a;
    for (n = 1; n <= L.N; n++) gel(u,n) = gmul(gsub(gel(u,n), a), gel(L.na,n));
  }
  setlg(A,i); return gerepilecopy(av, A);
}
GEN
asympnumraw(void *E, GEN (*f)(void *,GEN,long), long LIM, GEN alpha, long prec)
{
  pari_sp av = avma;
  double c, d, al;
  long i, B;
  GEN u, A;
  struct limit L;

  if (LIM < 0) return cgetg(1, t_VEC);
  c = get_c(alpha);
  d = get_accu(alpha);
  al = alpha? gtodouble(alpha): 1.0;
  B = prec2nbits(prec);
  L.N = ceil(solvedivlog(c * al * LIM / M_LN2, c * B));
  L.prec = nbits2prec(ceil(B + L.N / c + d * L.N));
  limit_init(&L, alpha, 1);
  u = get_u(E, f, L.N, L.prec);
  A = cgetg(LIM+2, t_VEC);
  for(i = 0; i <= LIM; i++)
  {
    GEN a = RgV_dotproduct(u,L.coef);
    long n;
    for (n = 1; n <= L.N; n++)
      gel(u,n) = gprec_wensure(gmul(gsub(gel(u,n), a), gel(L.na,n)), L.prec);
    gel(A,i+1) = gprec_wtrunc(a, prec);
  }
  return gerepilecopy(av, A);
}
GEN
asympnum0(GEN u, GEN alpha, long prec)
{ return asympnum((void*)u,get_fun(u, "asympnum"), alpha, prec); }
GEN
asympnumraw0(GEN u, long LIM, GEN alpha, long prec)
{ return asympnumraw((void*)u,get_fun(u, "asympnumraw"), LIM, alpha, prec); }
