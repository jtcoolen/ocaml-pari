/* Copyright (C) 2021  The PARI group.

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

static GEN
strtoint(char *s)
{
  long neg = (*s == '-');
  GEN n;
  if (!neg) return strtoi(s);
  n = strtoi(s+1); togglesign(n); return n;
}
/* export ? */
static GEN
mkZVn(long n, ...)
{
  va_list ap;
  GEN x;
  long i;
  va_start(ap,n);
  x = cgetg(n+1, t_VEC);
  for (i=1; i <= n; i++) gel(x,i) = strtoint(va_arg(ap, char*));
  va_end(ap); return x;
}

static GEN /* -t */
pol_mx(long v) { return deg1pol_shallow(gen_m1, gen_0, v); }
static GEN QT4(long k, long v)
{ switch(k) {
  case 1: return mkpoln(5, gen_1, pol_x(v), stoi(-6), pol_mx(v), gen_1);
  case 2: return mkpoln(5, gen_1, gen_0, pol_x(v), gen_0, gen_1);
  case 3: return mkpoln(5, gen_1, pol_x(v), gen_0, pol_x(v), gen_1);
  default: return NULL; }}
static GEN QT5(long k, long v)
{ GEN a3, a2, a1, a0;
  switch(k) {
  case 1:
  a3 = mkpoln(4, stoi(-2), stoi(-6), stoi(-10), stoi(-10)); setvarn(a3,v);
  a2 = mkpoln(5, gen_1, stoi(5), stoi(11), stoi(15), stoi(5)); setvarn(a2,v);
  a1 = mkpoln(4, gen_1, stoi(4), stoi(10), stoi(10)); setvarn(a1,v);
  return mkpoln(6, gen_1, pol_xn(2,v), a3, a2, a1, gen_1);
  case 2:
  a3 = deg1pol_shallow(gen_m1, stoi(-50), v);
  a1 = deg1pol_shallow(stoi(5), stoi(625), v);
  a0 = deg1pol_shallow(stoi(-3), gen_0, v);
  return mkpoln(6, gen_1, gen_0, a3, pol_mx(v), a1, a0);
  case 3:
  a2 = deg1pol_shallow(stoi(5), gen_0, v);
  a0 = deg2pol_shallow(gen_1, gen_m1, stoi(16), v);
  return mkpoln(6, gen_1, gen_0, stoi(10), a2, stoi(-15), a0);
  default: return NULL; }}
static GEN QT6(long k, long v)
{ GEN a5, a4, a3, a2, a1, a0;
  switch(k) {
  case 1:
  a5 = deg1pol_shallow(gen_2, gen_0, v);
  a4 = deg1pol_shallow(stoi(-5), stoi(-15), v);
  a2 = deg1pol_shallow(stoi(5), gen_0, v);
  a1 = deg1pol_shallow(gen_m2, stoi(-6), v);
  return mkpoln(7, gen_1, a5, a4, stoi(20), a2, a1, gen_1);
  case 2:
  a0 = deg2pol_shallow(stoi(3), gen_0, stoi(4), v);
  return mkpoln(7, gen_1, gen_0, stoi(6), gen_0, stoi(9), gen_0, a0);
  case 3:
  return mkpoln(7, gen_1, gen_0, stoi(6), gen_0, stoi(9), gen_0, pol_mx(v));
  case 4:
  a2 = deg1pol_shallow(gen_1, stoi(-3), v);
  return mkpoln(7, gen_1, gen_0, pol_x(v), gen_0, a2, gen_0, gen_m1);
  case 5:
  a4 = deg1pol_shallow(gen_1, stoi(-6), v);
  a3 = deg1pol_shallow(gen_2, gen_m2, v);
  a2 = deg1pol_shallow(gen_1, stoi(9), v);
  return mkpoln(7, gen_1, gen_0, a4, a3, a2, stoi(6), gen_1);
  case 6:
  a2 = mkpoln(5,stoi(-12),gen_0,stoi(-36),gen_0,gen_0); setvarn(a2,v);
  a0 = mkpoln(7,stoi(16),gen_0,stoi(48),gen_0,gen_0,gen_0,gen_0); setvarn(a0,v);
  return mkpoln(7, gen_1, gen_0, gen_0, gen_0, a2, gen_0, a0);
  case 7: return mkpoln(7, gen_1,gen_0,gen_0,gen_0,pol_x(v),gen_0,gen_m1);
  case 8:
  a0 = deg2pol_shallow(stoi(3), gen_0, stoi(4), v);
  return mkpoln(7, gen_1,gen_0,stoi(-3),gen_0,gen_0,gen_0,a0);
  case 9:
  a4 = deg2pol_shallow(stoi(3), gen_0, stoi(-6), v);
  a3 = deg2pol_shallow(gen_m2, gen_0, stoi(4), v);
  return mkpoln(7, gen_1,gen_0,a4,a3,stoi(9),stoi(-12),stoi(4));
  case 10:
  a0 = deg2pol_shallow(gen_m1, gen_0, stoi(-1024), v);
  return mkpoln(7, gen_1,stoi(-12),stoi(36),gen_0,gen_0,gen_0,a0);
  case 11:
  a0 = deg1pol_shallow(gen_1, stoi(4), v);
  return mkpoln(7, gen_1,gen_0,stoi(-3),gen_0,gen_0,gen_0,a0);
  case 12:
  a5 = deg2pol_shallow(stoi(10), gen_0, stoi(-50), v);
  a4 = gtopoly(mkvecsmall5(55, 0,-550, 0, 1375), v);
  a3 = gtopoly(mkvecsmalln(7, 140L, 0L,-2100L, 0L, 10500L, 0L,-17500L), v);
  a2 = gtopoly(mkvecsmalln(9, 175L, 0L,-3500L, 0L, 26250L, 0L,-87500L, 0L, 109375L), v);
  a1 = gtopoly(mkvecsmalln(11, 106L, 0L,-1370L, 0L, 900L, 0L, 59500L, 0L,-308750L, 0L, 468750L), v);
  a0 = gtopoly(mkvecsmalln(13, 25L, 0L,-750L, 0L, 9375L, 0L,-62500L, 0L, 234375L, 0L,-468750L, 0L, 390625L), v);
  return mkpoln(7, gen_1,a5,a4,a3,a2,a1,a0);
  case 13:
  return mkpoln(7, gen_1,gen_m2,gen_1,gen_0,gen_0,gen_0,pol_mx(v));
  case 14:
  return mkpoln(7, gen_1,stoi(4),stoi(20),gen_0,gen_0,pol_mx(v),pol_x(v));
  default: return NULL; }}
static GEN QT7(long k, long v)
{ GEN a6, a5, a4, a3, a2, a1, a0;
  switch(k) {
  case 1:
  a6 = gtopoly(mkvecsmall4(1, 2,-1, 13), v);
  a5 = gtopoly(mkvecsmalln(6, 3L,-3L, 9L, 24L,-21L, 54L), v);
  a4 = gtopoly(mkvecsmalln(8, 3L,-9L, 27L,-22L, 6L, 84L,-121L, 75L), v);
  a3 = gtopoly(mkvecsmalln(10, 1L,-6L, 22L,-57L, 82L,-70L,-87L, 140L,-225L,-2L), v);
  a2 = gtopoly(mkvecsmalln(11, -1L, 5L,-25L, 61L,-126L, 117L,-58L,-155L, 168L,-80L,-44L), v);
  a1 = gtopoly(mkvecsmalln(11, -1L, 8L,-30L, 75L,-102L, 89L, 34L,-56L, 113L, 42L,-17L), v);
  a0 = gtopoly(mkvecsmalln(10, 1L,-7L, 23L,-42L, 28L, 19L,-60L,-2L, 16L,-1L), v);
  return mkpoln(8, gen_1,a6,a5,a4,a3,a2,a1,a0);
  case 2:
  a5 = gtopoly(mkvecsmall4(-147,-735,-441,-21), v);
  a4 = gtopoly(mkvecsmall5(-686,-3920,-4508,-1568,-70), v);
  a3 = gtopoly(mkvecsmalln(7, 7203L, 67914L, 183505L, 107996L, 8085L,-1862L,-105L), v);
  a2 = gtopoly(mkvecsmalln(8, 67228L, 547428L, 1373372L, 1227940L, 416500L, 38220L,-588L,-84L), v);
  a1 = gtopoly(mkvecsmalln(10, -117649L,-1563051L,-6809236L,-10708460L,-4050830L, 788214L, 402780L, 37828L, 343L,-35L), v);
  a0 = gtopoly(mkvecsmalln(11, -1647086L,-16893436L,-56197806L,-69977488L,-44893212L,-13304872L,-624652L, 103152L, 11466L, 196L,-6L), v);
  return mkpoln(8, gen_1,gen_0,a5,a4,a3,a2,a1,a0);
  case 3:
  a5 = gtopoly(mkvecsmalln(7, -21L,0L,-1176L,147L,-20580L,3969L,-107163L), v);
  a4 = gtopoly(mkvecsmalln(10, -21L,49L,-1715L,4214L,-51107L,129850L,-653905L,1648458L,-3000564L,6751269L), v);
  a3 = gtopoly(mkvecsmalln(13, 91L,98L,9849L,8673L,427133L,291354L,9385460L,4618152L,108334149L,35173278L,608864445L,114771573L,1275989841L), v);
  a2 = gtopoly(mkZVn(16, "112","-49","14651","-10682","800513","-821730","23571744","-30983190","401636536","-628991685","3929562693","-6832117530","20190045015","-35916751080","40831674912","-68903451414]"), v);
  a1 = gtopoly(mkZVn(19, "-84","-98","-14896","-16709","-1127098","-1228626","-47347279","-51034970","-1201635330","-1316073164","-18735012261","-21705143929","-173551408569","-224605199322","-861876002232","-1329675932088","-1728966234555","-3376269119286","0"), v);
  a0 = gtopoly(mkZVn(22, "-97","-14","-19803","-903","-1765232","84609","-89982172","11950757","-2882068329","588528171","-59885187418","15296374002","-801314604769","222442927665","-6560078164731","1705024373220","-28577589326937","5543939564730","-38647180304208","4961048501808","74415727527120","25115308040403"), v);
  return mkpoln(8, gen_1,gen_0,a5,a4,a3,a2,a1,a0);
  case 4:
  a4 = deg2pol_shallow(stoi(-7), gen_0, stoi(98), v);
  a3 = deg1pol_shallow(stoi(28), stoi(441), v);
  a2 = gtopoly(mkvecsmall4(-35,-112,-196, 343), v);
  a1 = deg2pol_shallow(stoi(7), stoi(196), stoi(1372), v);
  a0 = gtopoly(mkvecsmalln(6, -1L,-30L,-259L,-588L,-1372L, 0L), v);
  return mkpoln(8, gen_1,stoi(7),stoi(42),a4,a3,a2,a1,a0);
  case 5:
  a3 = deg1pol_shallow(stoi(12), stoi(7203), v);
  a2 = deg1pol_shallow(stoi(-30), gen_0, v);
  a1 = deg1pol_shallow(stoi(28), stoi(-117649), v);
  a0 = deg1pol_shallow(stoi(-9), gen_0, v);
  return mkpoln(8, gen_1,gen_0,stoi(-147),pol_mx(v),a3,a2,a1,a0);
  default: return NULL; }}

static GEN
nflistQTfile(long n, long t)
{
  pariFILE *F;
  GEN z;
  char *f = stack_sprintf("%s/nflistdata/%ld/%ld/QT", pari_datadir, n, t);
  F = pari_fopengz(f); if (!F) return NULL;
  z = gp_read_stream(F->file); pari_fclose(F); return z;
}

static GEN
nfmakeQT(long deg, long k, long v)
{
  long i, l;
  GEN P;
  switch(deg) {
  case 4: P = QT4(k, v); break;
  case 5: P = QT5(k, v); break;
  case 6: P = QT6(k, v); break;
  case 7: P = QT7(k, v); break;
  default: P = nflistQTfile(deg, k);
  }
  if (!P)
    pari_err_IMPL(stack_sprintf("group %ldT%ld in nflist / Q(T)", deg,k));
  if (deg <= 7) return P;
  l = lg(P);
  for (i = 1; i < l; i++)
  {
    GEN p = gel(P,i);
    if (typ(p) != t_INT) gel(P,i) = RgV_to_RgX(p, v);
  }
  return RgV_to_RgX(P, 0);
}

static GEN
nfmakeAnQT(long n, long v)
{
  GEN A, P = vec_ei(n + 1, 1);
  long s;
  if (odd(n))
  {
    s = (n & 3L) == 1? 1: -1;
    A = sqru(n-2); setsigne(A, s);
    gel(P,2) = monomial(sqru(n), 1, v);
    gel(P,n) = s > 0? gen_1: gen_m1;
    gel(P,n+1) = monomial(A, 1, v);
  }
  else
  {
    s = (n & 3L)? -1 : 1;
    gel(P,2) = utoineg(n);
    gel(P,n+1) = deg2pol_shallow(stoi(s), gen_0, powuu(n-1,n-1), v);
  }
  return RgV_to_RgX_reverse(P, 0);
}

static GEN
nfmakeSnQT(long n, long v)
{
  GEN P = vec_ei(n + 1, 1);
  gel(P,n) = pol_x(v);
  gel(P,n+1) = gen_1; return RgV_to_RgX_reverse(P, 0);
}

GEN
nflistQT(long n, long t, long v)
{
  if (varncmp(0,v) >= 0)
    pari_err(e_MISC, "incorrect variable in nflist / Q(T)");
  if (n == 1) return pol_x(0);
  if (n == 2) return deg2pol_shallow(gen_1, pol_mx(v), gen_1, 0);
  if (t == -1) return nfmakeSnQT(n, v);
  if (t == -2) return nfmakeAnQT(n, v);
  return nfmakeQT(n, t, v);
}
