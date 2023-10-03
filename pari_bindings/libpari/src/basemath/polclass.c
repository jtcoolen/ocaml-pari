/* Copyright (C) 2014  The PARI group.

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

#define DEBUGLEVEL DEBUGLEVEL_polclass

#define dbg_printf(lvl) if (DEBUGLEVEL >= (lvl) + 3) err_printf

static GEN
pcp_get_L(GEN G) { return gmael(G,1,1)+1; }
static GEN
pcp_get_n(GEN G) { return gmael(G,1,2)+1; }
static GEN
pcp_get_m(GEN G) { return gmael(G,1,4)+1; }
static GEN
pcp_get_o(GEN G) { return gmael(G,1,3)+1; }
static GEN
pcp_get_r(GEN G) { return gmael(G,1,5)+1; }
static GEN
pcp_get_orient_p(GEN G) { return gmael(G,3,1)+1; }
static GEN
pcp_get_orient_q(GEN G) { return gmael(G,3,2)+1; }
static GEN
pcp_get_orient_reps(GEN G) { return gmael(G,3,3)+1; }
static long
pcp_get_L0(GEN G) { return mael(G,2,1); }
static long
pcp_get_k(GEN G) { return mael(G,2,2); }
static long
pcp_get_h(GEN G) { return mael(G,4,1); }
static long
pcp_get_inv(GEN G) { return mael(G,4,2); }
static long
pcp_get_D(GEN G) { return mael(G,4,3); }
static long
pcp_get_D0(GEN G) { return mael(G,4,4); }
static long
pcp_get_u(GEN G) { return mael(G,4,5); }
static GEN
pcp_get_fau(GEN G) { return gel(G,5); }

/* SECTION: Functions dedicated to finding a j-invariant with a given
 * trace. */

static void
hasse_bounds(long *low, long *high, long p)
{ long u = usqrt(4*p); *low = p + 1 - u; *high = p + 1 + u; }

/* a / b : a and b are from factoru and b must divide a exactly */
INLINE GEN
famatsmall_divexact(GEN a, GEN b)
{
  GEN a1 = gel(a,1), a2 = gel(a,2), c1, c2;
  GEN b1 = gel(b,1), b2 = gel(b,2);
  long i, j, k, la = lg(a1);
  c1 = cgetg(la, t_VECSMALL);
  c2 = cgetg(la, t_VECSMALL);
  for (i = j = k = 1; j < la; j++)
  {
    c1[k] = a1[j];
    c2[k] = a2[j];
    if (a1[j] == b1[i]) { c2[k] -= b2[i++]; if (!c2[k]) continue; }
    k++;
  }
  setlg(c1, k);
  setlg(c2, k); return mkvec2(c1,c2);
}

/* This is Sutherland, 2009, TestCurveOrder.
 *
 * [a4, a6] and p specify an elliptic curve over FF_p.  N0,N1 are the two
 * possible curve orders (N0+N1 = 2p+2), and n0,n1 their factoru */
static long
test_curve_order(norm_eqn_t ne, ulong a4, ulong a6,
  long N0, long N1, GEN n0, GEN n1, long hasse_low, long hasse_high)
{
  pari_sp ltop = avma, av;
  ulong a4t, a6t, p = ne->p, pi = ne->pi, T = ne->T;
  long m0, m1;
  int swapped = 0, first = 1;

  if (p <= 11) {
    long card = (long)p + 1 - Fl_elltrace(a4, a6, p);
    return card == N0 || card == N1;
  }
  m0 = m1 = 1;
  for (av = avma;;)
  {
    long a1, x, n_s;
    GEN Q = random_Flj_pre(a4, a6, p, pi);
    if (m0 == 1)
      n_s = Flj_order_ufact(Q, N0, n0, a4, p, pi);
    else if (N0 % m0) n_s = 0;
    else
    { /* m0 | N0 */
      GEN fa0 = famatsmall_divexact(n0, factoru(m0));
      Q = Flj_mulu_pre(Q, m0, a4, p, pi);
      n_s = Flj_order_ufact(Q, N0 / m0, fa0, a4, p, pi);
    }
    if (n_s == 0)
    { /* If m0 divides N1 and m1 divides N0 and N0 < N1, then swap */
      if (swapped || N1 % m0 || N0 % m1) return gc_long(ltop,0);
      swapspec(n0, n1, N0, N1); swapped = 1; continue;
    }

    m0 *= n_s; a1 = (2 * p + 2) % m1;
    for (x = ceildivuu(hasse_low, m0) * m0; x <= hasse_high; x += m0)
      if ((x % m1) == a1 && x != N0 && x != N1) break;
    /* every x was either N0 or N1, so we return true */
    if (x > hasse_high) return gc_long(ltop,1);

    /* [a4, a6] is the given curve and [a4t, a6t] is its quadratic twist */
    if (first) { Fl_elltwist_disc(a4, a6, T, p, &a4t, &a6t); first = 0; }
    lswap(a4, a4t);
    lswap(a6, a6t);
    lswap(m0, m1); set_avma(av);
  }
}

static GEN
random_FleV(GEN x, GEN a6, ulong p, ulong pi)
{ pari_APPLY_type(t_VEC, random_Fle_pre(uel(x,i), uel(a6,i), p, pi)) }

/* START Code from AVSs "torcosts.h" */
struct torctab_rec {
  int m;
  int fix2, fix3;
  int N;
  int s2_flag;
  int t3_flag;
  double rating;
};

/* These costs assume p=2 mod 3, 3 mod 4 and not 1 mod N */
static struct torctab_rec torctab1[] = {
{ 11, 1, 1, 11, 1, 1, 0.047250 },
{ 33, 1, 0, 11, 1, 2, 0.047250 },
{ 22, 1, 1, 11, 3, 1, 0.055125 },
{ 66, 1, 0, 11, 3, 2, 0.055125 },
{ 11, 1, 0, 11, 1, 0, 0.058000 },
{ 13, 1, 1, 13, 1, 1, 0.058542 },
{ 39, 1, 0, 13, 1, 2, 0.058542 },
{ 22, 0, 1, 11, 2, 1, 0.061333 },
{ 66, 0, 0, 11, 2, 2, 0.061333 },
{ 22, 1, 0, 11, 3, 0, 0.061750 },
{ 14, 1, 1, 14, 3, 1, 0.062500 },
{ 42, 1, 0, 14, 3, 2, 0.062500 },
{ 26, 1, 1, 13, 3, 1, 0.064583 },
{ 78, 1, 0, 13, 3, 2, 0.064583 },
{ 28, 0, 1, 14, 4, 1, 0.065625 },
{ 84, 0, 0, 14, 4, 2, 0.065625 },
{ 7, 1, 1, 7, 1, 1, 0.068750 },
{ 13, 1, 0, 13, 1, 0, 0.068750 },
{ 21, 1, 0, 7, 1, 2, 0.068750 },
{ 26, 1, 0, 13, 3, 0, 0.069583 },
{ 17, 1, 1, 17, 1, 1, 0.069687 },
{ 51, 1, 0, 17, 1, 2, 0.069687 },
{ 11, 0, 1, 11, 0, 1, 0.072500 },
{ 33, 0, 0, 11, 0, 2, 0.072500 },
{ 44, 1, 0, 11, 130, 0, 0.072667 },
{ 52, 0, 1, 13, 4, 1, 0.073958 },
{ 156, 0, 0, 13, 4, 2, 0.073958 },
{ 34, 1, 1, 17, 3, 1, 0.075313 },
{ 102, 1, 0, 17, 3, 2, 0.075313 },
{ 15, 1, 0, 15, 1, 0, 0.075625 },
{ 13, 0, 1, 13, 0, 1, 0.076667 },
{ 39, 0, 0, 13, 0, 2, 0.076667 },
{ 44, 0, 0, 11, 4, 0, 0.076667 },
{ 30, 1, 0, 15, 3, 0, 0.077188 },
{ 22, 0, 0, 11, 2, 0, 0.077333 },
{ 34, 1, 0, 17, 3, 0, 0.077969 },
{ 17, 1, 0, 17, 1, 0, 0.078750 },
{ 14, 0, 1, 14, 0, 1, 0.080556 },
{ 28, 0, 0, 14, 4, 0, 0.080556 },
{ 42, 0, 0, 14, 0, 2, 0.080556 },
{ 7, 1, 0, 7, 1, 0, 0.080833 },
{ 9, 1, 0, 9, 1, 0, 0.080833 },
{ 68, 0, 1, 17, 4, 1, 0.081380 },
{ 204, 0, 0, 17, 4, 2, 0.081380 },
{ 52, 0, 0, 13, 4, 0, 0.082292 },
{ 10, 1, 1, 10, 3, 1, 0.084687 },
{ 17, 0, 1, 17, 0, 1, 0.084687 },
{ 51, 0, 0, 17, 0, 2, 0.084687 },
{ 20, 0, 1, 10, 4, 1, 0.085938 },
{ 60, 0, 0, 10, 4, 2, 0.085938 },
{ 19, 1, 1, 19, 1, 1, 0.086111 },
{ 57, 1, 0, 19, 1, 2, 0.086111 },
{ 68, 0, 0, 17, 4, 0, 0.088281 },
{ 38, 1, 1, 19, 3, 1, 0.089514 },
{ 114, 1, 0, 19, 3, 2, 0.089514 },
{ 20, 0, 0, 10, 4, 0, 0.090625 },
{ 36, 0, 0, 18, 4, 0, 0.090972 },
{ 26, 0, 0, 13, 2, 0, 0.091667 },
{ 11, 0, 0, 11, 0, 0, 0.092000 },
{ 19, 1, 0, 19, 1, 0, 0.092778 },
{ 38, 1, 0, 19, 3, 0, 0.092778 },
{ 14, 1, 0, 7, 3, 0, 0.092917 },
{ 18, 1, 0, 9, 3, 0, 0.092917 },
{ 76, 0, 1, 19, 4, 1, 0.095255 },
{ 228, 0, 0, 19, 4, 2, 0.095255 },
{ 10, 0, 1, 10, 0, 1, 0.096667 },
{ 13, 0, 0, 13, 0, 0, 0.096667 },
{ 30, 0, 0, 10, 0, 2, 0.096667 },
{ 19, 0, 1, 19, 0, 1, 0.098333 },
{ 57, 0, 0, 19, 0, 2, 0.098333 },
{ 17, 0, 0, 17, 0, 0, 0.100000 },
{ 23, 1, 1, 23, 1, 1, 0.100227 },
{ 69, 1, 0, 23, 1, 2, 0.100227 },
{ 7, 0, 1, 7, 0, 1, 0.100833 },
{ 21, 0, 0, 7, 0, 2, 0.100833 },
{ 76, 0, 0, 19, 4, 0, 0.102083 },
{ 14, 0, 0, 14, 0, 0, 0.102222 },
{ 18, 0, 0, 9, 2, 0, 0.102222 },
{ 5, 1, 1, 5, 1, 1, 0.103125 },
{ 46, 1, 1, 23, 3, 1, 0.104318 },
{ 138, 1, 0, 23, 3, 2, 0.104318 },
{ 23, 1, 0, 23, 1, 0, 0.105682 },
{ 46, 1, 0, 23, 3, 0, 0.106705 },
{ 92, 0, 1, 23, 4, 1, 0.109091 },
{ 276, 0, 0, 23, 4, 2, 0.109091 },
{ 19, 0, 0, 19, 0, 0, 0.110000 },
{ 23, 0, 1, 23, 0, 1, 0.112273 },
{ 69, 0, 0, 23, 0, 2, 0.112273 },
{ 7, 0, 0, 7, 0, 0, 0.113333 },
{ 9, 0, 0, 9, 0, 0, 0.113333 },
{ 92, 0, 0, 23, 4, 0, 0.113826 },
{ 16, 0, 1, 16, 0, 1, 0.118125 },
{ 48, 0, 0, 16, 0, 2, 0.118125 },
{ 5, 1, 0, 5, 1, 0, 0.121250 },
{ 15, 0, 0, 15, 0, 0, 0.121250 },
{ 10, 0, 0, 10, 0, 0, 0.121667 },
{ 23, 0, 0, 23, 0, 0, 0.123182 },
{ 12, 0, 0, 12, 0, 0, 0.141667 },
{ 5, 0, 1, 5, 0, 1, 0.145000 },
{ 16, 0, 0, 16, 0, 0, 0.145000 },
{ 8, 0, 1, 8, 0, 1, 0.151250 },
{ 29, 1, 1, 29, 1, 1, 0.153036 },
{ 87, 1, 0, 29, 1, 2, 0.153036 },
{ 25, 0, 0, 25, 0, 0, 0.155000 },
{ 58, 1, 1, 29, 3, 1, 0.156116 },
{ 174, 1, 0, 29, 3, 2, 0.156116 },
{ 29, 1, 0, 29, 1, 0, 0.157500 },
{ 58, 1, 0, 29, 3, 0, 0.157500 },
{ 116, 0, 1, 29, 4, 1, 0.161086 },
{ 29, 0, 1, 29, 0, 1, 0.163393 },
{ 87, 0, 0, 29, 0, 2, 0.163393 },
{ 116, 0, 0, 29, 4, 0, 0.163690 },
{ 5, 0, 0, 5, 0, 0, 0.170000 },
{ 8, 0, 0, 8, 0, 0, 0.170000 },
{ 29, 0, 0, 29, 0, 0, 0.171071 },
{ 31, 1, 1, 31, 1, 1, 0.186583 },
{ 93, 1, 0, 31, 1, 2, 0.186583 },
{ 62, 1, 1, 31, 3, 1, 0.189750 },
{ 186, 1, 0, 31, 3, 2, 0.189750 },
{ 31, 1, 0, 31, 1, 0, 0.191333 },
{ 62, 1, 0, 31, 3, 0, 0.192167 },
{ 124, 0, 1, 31, 4, 1, 0.193056 },
{ 31, 0, 1, 31, 0, 1, 0.195333 },
{ 93, 0, 0, 31, 0, 2, 0.195333 },
{ 124, 0, 0, 31, 4, 0, 0.197917 },
{ 2, 1, 1, 2, 3, 1, 0.200000 },
{ 6, 1, 0, 2, 3, 2, 0.200000 },
{ 31, 0, 0, 31, 0, 0, 0.206667 },
{ 4, 1, 1, 4, 130, 1, 0.214167 },
{ 6, 0, 0, 6, 0, 0, 0.226667 },
{ 3, 1, 0, 3, 1, 0, 0.230000 },
{ 4, 0, 1, 4, 0, 1, 0.241667 },
{ 4, 1, 0, 2, 130, 0, 0.266667 },
{ 4, 0, 0, 4, 0, 0, 0.283333 },
{ 3, 0, 0, 3, 0, 0, 0.340000 },
{ 1, 1, 1, 1, 1, 1, 0.362500 },
{ 2, 0, 1, 2, 0, 1, 0.386667 },
{ 1, 1, 0, 1, 1, 0, 0.410000 },
{ 2, 0, 0, 2, 0, 0, 0.453333 },
};

static struct torctab_rec torctab2[] = {
{ 11, 1, 1, 11, 1, 1, 0.047250 },
{ 33, 1, 0, 11, 1, 2, 0.047250 },
{ 22, 1, 1, 11, 3, 1, 0.055125 },
{ 66, 1, 0, 11, 3, 2, 0.055125 },
{ 13, 1, 1, 13, 1, 1, 0.057500 },
{ 39, 1, 0, 13, 1, 2, 0.057500 },
{ 11, 1, 0, 11, 1, 0, 0.058000 },
{ 22, 0, 1, 11, 2, 1, 0.061333 },
{ 66, 0, 0, 11, 2, 2, 0.061333 },
{ 14, 1, 1, 14, 3, 1, 0.061458 },
{ 42, 1, 0, 14, 3, 2, 0.061458 },
{ 22, 1, 0, 11, 3, 0, 0.061750 },
{ 26, 1, 1, 13, 3, 1, 0.064062 },
{ 78, 1, 0, 13, 3, 2, 0.064062 },
{ 28, 0, 1, 14, 4, 1, 0.065625 },
{ 84, 0, 0, 14, 4, 2, 0.065625 },
{ 13, 1, 0, 13, 1, 0, 0.066667 },
{ 26, 1, 0, 13, 3, 0, 0.069583 },
{ 17, 1, 1, 17, 1, 1, 0.069687 },
{ 51, 1, 0, 17, 1, 2, 0.069687 },
{ 11, 0, 1, 11, 0, 1, 0.070000 },
{ 33, 0, 0, 11, 0, 2, 0.070000 },
{ 7, 1, 1, 7, 1, 1, 0.070417 },
{ 21, 1, 0, 7, 1, 2, 0.070417 },
{ 15, 1, 0, 15, 1, 0, 0.072500 },
{ 52, 0, 1, 13, 4, 1, 0.073090 },
{ 156, 0, 0, 13, 4, 2, 0.073090 },
{ 34, 1, 1, 17, 3, 1, 0.074219 },
{ 102, 1, 0, 17, 3, 2, 0.074219 },
{ 7, 1, 0, 7, 1, 0, 0.076667 },
{ 13, 0, 1, 13, 0, 1, 0.076667 },
{ 39, 0, 0, 13, 0, 2, 0.076667 },
{ 44, 0, 0, 11, 4, 0, 0.076667 },
{ 17, 1, 0, 17, 1, 0, 0.077188 },
{ 22, 0, 0, 11, 2, 0, 0.077333 },
{ 34, 1, 0, 17, 3, 0, 0.077969 },
{ 30, 1, 0, 15, 3, 0, 0.080312 },
{ 14, 0, 1, 14, 0, 1, 0.080556 },
{ 28, 0, 0, 14, 4, 0, 0.080556 },
{ 42, 0, 0, 14, 0, 2, 0.080556 },
{ 9, 1, 0, 9, 1, 0, 0.080833 },
{ 68, 0, 1, 17, 4, 1, 0.081380 },
{ 204, 0, 0, 17, 4, 2, 0.081380 },
{ 52, 0, 0, 13, 4, 0, 0.082292 },
{ 10, 1, 1, 10, 3, 1, 0.083125 },
{ 20, 0, 1, 10, 4, 1, 0.083333 },
{ 60, 0, 0, 10, 4, 2, 0.083333 },
{ 17, 0, 1, 17, 0, 1, 0.084687 },
{ 51, 0, 0, 17, 0, 2, 0.084687 },
{ 19, 1, 1, 19, 1, 1, 0.084722 },
{ 57, 1, 0, 19, 1, 2, 0.084722 },
{ 11, 0, 0, 11, 0, 0, 0.087000 },
{ 68, 0, 0, 17, 4, 0, 0.088281 },
{ 38, 1, 1, 19, 3, 1, 0.090139 },
{ 114, 1, 0, 19, 3, 2, 0.090139 },
{ 36, 0, 0, 18, 4, 0, 0.090972 },
{ 19, 1, 0, 19, 1, 0, 0.091389 },
{ 26, 0, 0, 13, 2, 0, 0.091667 },
{ 13, 0, 0, 13, 0, 0, 0.092500 },
{ 38, 1, 0, 19, 3, 0, 0.092778 },
{ 14, 1, 0, 7, 3, 0, 0.092917 },
{ 18, 1, 0, 9, 3, 0, 0.092917 },
{ 20, 0, 0, 10, 4, 0, 0.095833 },
{ 76, 0, 1, 19, 4, 1, 0.096412 },
{ 228, 0, 0, 19, 4, 2, 0.096412 },
{ 17, 0, 0, 17, 0, 0, 0.096875 },
{ 19, 0, 1, 19, 0, 1, 0.098056 },
{ 57, 0, 0, 19, 0, 2, 0.098056 },
{ 23, 1, 1, 23, 1, 1, 0.100682 },
{ 69, 1, 0, 23, 1, 2, 0.100682 },
{ 7, 0, 1, 7, 0, 1, 0.100833 },
{ 21, 0, 0, 7, 0, 2, 0.100833 },
{ 30, 0, 0, 15, 2, 0, 0.100833 },
{ 76, 0, 0, 19, 4, 0, 0.102083 },
{ 14, 0, 0, 14, 0, 0, 0.102222 },
{ 5, 1, 1, 5, 1, 1, 0.103125 },
{ 46, 1, 1, 23, 3, 1, 0.104034 },
{ 138, 1, 0, 23, 3, 2, 0.104034 },
{ 23, 1, 0, 23, 1, 0, 0.104545 },
{ 7, 0, 0, 7, 0, 0, 0.105000 },
{ 10, 0, 1, 10, 0, 1, 0.105000 },
{ 16, 0, 1, 16, 0, 1, 0.105417 },
{ 48, 0, 0, 16, 0, 2, 0.105417 },
{ 46, 1, 0, 23, 3, 0, 0.106705 },
{ 18, 0, 0, 9, 2, 0, 0.107778 },
{ 92, 0, 1, 23, 4, 1, 0.108239 },
{ 276, 0, 0, 23, 4, 2, 0.108239 },
{ 19, 0, 0, 19, 0, 0, 0.110000 },
{ 23, 0, 1, 23, 0, 1, 0.111136 },
{ 69, 0, 0, 23, 0, 2, 0.111136 },
{ 9, 0, 0, 9, 0, 0, 0.113333 },
{ 10, 0, 0, 10, 0, 0, 0.113333 },
{ 92, 0, 0, 23, 4, 0, 0.113826 },
{ 5, 1, 0, 5, 1, 0, 0.115000 },
{ 15, 0, 0, 15, 0, 0, 0.115000 },
{ 23, 0, 0, 23, 0, 0, 0.120909 },
{ 8, 0, 1, 8, 0, 1, 0.126042 },
{ 24, 0, 0, 8, 0, 2, 0.126042 },
{ 16, 0, 0, 16, 0, 0, 0.127188 },
{ 8, 0, 0, 8, 0, 0, 0.141667 },
{ 25, 0, 1, 25, 0, 1, 0.144000 },
{ 5, 0, 1, 5, 0, 1, 0.151250 },
{ 12, 0, 0, 12, 0, 0, 0.152083 },
{ 29, 1, 1, 29, 1, 1, 0.153929 },
{ 87, 1, 0, 29, 1, 2, 0.153929 },
{ 25, 0, 0, 25, 0, 0, 0.155000 },
{ 58, 1, 1, 29, 3, 1, 0.155045 },
{ 174, 1, 0, 29, 3, 2, 0.155045 },
{ 29, 1, 0, 29, 1, 0, 0.156429 },
{ 58, 1, 0, 29, 3, 0, 0.157857 },
{ 116, 0, 1, 29, 4, 1, 0.158631 },
{ 116, 0, 0, 29, 4, 0, 0.163542 },
{ 29, 0, 1, 29, 0, 1, 0.164286 },
{ 87, 0, 0, 29, 0, 2, 0.164286 },
{ 29, 0, 0, 29, 0, 0, 0.169286 },
{ 5, 0, 0, 5, 0, 0, 0.170000 },
{ 31, 1, 1, 31, 1, 1, 0.187000 },
{ 93, 1, 0, 31, 1, 2, 0.187000 },
{ 62, 1, 1, 31, 3, 1, 0.188500 },
{ 186, 1, 0, 31, 3, 2, 0.188500 },
{ 31, 1, 0, 31, 1, 0, 0.191333 },
{ 62, 1, 0, 31, 3, 0, 0.192083 },
{ 124, 0, 1, 31, 4, 1, 0.193472 },
{ 31, 0, 1, 31, 0, 1, 0.196167 },
{ 93, 0, 0, 31, 0, 2, 0.196167 },
{ 124, 0, 0, 31, 4, 0, 0.197083 },
{ 2, 1, 1, 2, 3, 1, 0.200000 },
{ 6, 1, 0, 2, 3, 2, 0.200000 },
{ 31, 0, 0, 31, 0, 0, 0.205000 },
{ 6, 0, 0, 6, 0, 0, 0.226667 },
{ 3, 1, 0, 3, 1, 0, 0.230000 },
{ 4, 0, 1, 4, 0, 1, 0.241667 },
{ 4, 0, 0, 4, 0, 0, 0.283333 },
{ 3, 0, 0, 3, 0, 0, 0.340000 },
{ 1, 1, 1, 1, 1, 1, 0.362500 },
{ 2, 0, 1, 2, 0, 1, 0.370000 },
{ 1, 1, 0, 1, 1, 0, 0.385000 },
{ 2, 0, 0, 2, 0, 0, 0.453333 },
};

static struct torctab_rec torctab3[] = {
{ 66, 1, 0, 11, 3, 2, 0.040406 },
{ 33, 1, 0, 11, 1, 2, 0.043688 },
{ 78, 1, 0, 13, 3, 2, 0.045391 },
{ 132, 1, 0, 11, 130, 2, 0.046938 },
{ 39, 1, 0, 13, 1, 2, 0.047656 },
{ 102, 1, 0, 17, 3, 2, 0.049922 },
{ 42, 1, 0, 14, 3, 2, 0.050000 },
{ 51, 1, 0, 17, 1, 2, 0.051680 },
{ 132, 0, 0, 11, 4, 2, 0.052188 },
{ 156, 1, 0, 13, 130, 2, 0.053958 },
{ 156, 0, 0, 13, 4, 2, 0.054818 },
{ 84, 1, 0, 14, 130, 2, 0.055000 },
{ 15, 1, 0, 15, 1, 0, 0.056719 },
{ 204, 0, 0, 17, 4, 2, 0.057227 },
{ 114, 1, 0, 19, 3, 2, 0.057500 },
{ 11, 1, 0, 11, 1, 0, 0.058000 },
{ 66, 0, 0, 11, 2, 2, 0.058000 },
{ 57, 1, 0, 19, 1, 2, 0.059062 },
{ 30, 1, 0, 15, 3, 0, 0.059063 },
{ 84, 0, 0, 14, 4, 2, 0.060677 },
{ 22, 1, 0, 11, 3, 0, 0.061750 },
{ 78, 0, 0, 13, 2, 2, 0.063542 },
{ 228, 0, 0, 19, 4, 2, 0.063889 },
{ 21, 1, 0, 7, 1, 2, 0.065000 },
{ 138, 1, 0, 23, 3, 2, 0.065028 },
{ 69, 1, 0, 23, 1, 2, 0.066903 },
{ 13, 1, 0, 13, 1, 0, 0.068750 },
{ 102, 0, 0, 17, 2, 2, 0.068906 },
{ 26, 1, 0, 13, 3, 0, 0.069583 },
{ 51, 0, 0, 17, 0, 2, 0.070312 },
{ 60, 1, 0, 15, 130, 0, 0.071094 },
{ 276, 0, 0, 23, 4, 2, 0.071236 },
{ 39, 0, 0, 13, 0, 2, 0.071250 },
{ 33, 0, 0, 11, 0, 2, 0.072750 },
{ 44, 1, 0, 11, 130, 0, 0.073500 },
{ 60, 0, 0, 15, 4, 0, 0.073828 },
{ 9, 1, 0, 9, 1, 0, 0.074097 },
{ 30, 0, 0, 15, 2, 0, 0.075625 },
{ 57, 0, 0, 19, 0, 2, 0.075625 },
{ 7, 1, 0, 7, 1, 0, 0.076667 },
{ 44, 0, 0, 11, 4, 0, 0.076667 },
{ 22, 0, 0, 11, 2, 0, 0.077333 },
{ 17, 1, 0, 17, 1, 0, 0.078750 },
{ 34, 1, 0, 17, 3, 0, 0.078750 },
{ 69, 0, 0, 23, 0, 2, 0.079943 },
{ 28, 0, 0, 14, 4, 0, 0.080556 },
{ 42, 0, 0, 14, 0, 2, 0.080833 },
{ 52, 0, 0, 13, 4, 0, 0.082292 },
{ 14, 1, 1, 14, 3, 1, 0.083333 },
{ 36, 0, 0, 18, 4, 0, 0.083391 },
{ 18, 1, 0, 9, 3, 0, 0.085174 },
{ 68, 0, 0, 17, 4, 0, 0.089583 },
{ 15, 0, 0, 15, 0, 0, 0.090938 },
{ 19, 1, 0, 19, 1, 0, 0.091389 },
{ 26, 0, 0, 13, 2, 0, 0.091667 },
{ 11, 0, 0, 11, 0, 0, 0.092000 },
{ 13, 0, 0, 13, 0, 0, 0.092500 },
{ 38, 1, 0, 19, 3, 0, 0.092778 },
{ 14, 1, 0, 7, 3, 0, 0.092917 },
{ 18, 0, 0, 9, 2, 0, 0.093704 },
{ 174, 1, 0, 29, 3, 2, 0.095826 },
{ 20, 0, 0, 10, 4, 0, 0.095833 },
{ 96, 1, 0, 16, 133, 2, 0.096562 },
{ 21, 0, 0, 21, 0, 0, 0.096875 },
{ 87, 1, 0, 29, 1, 2, 0.096964 },
{ 17, 0, 0, 17, 0, 0, 0.100000 },
{ 348, 0, 0, 29, 4, 2, 0.100558 },
{ 76, 0, 0, 19, 4, 0, 0.100926 },
{ 14, 0, 0, 14, 0, 0, 0.102222 },
{ 9, 0, 0, 9, 0, 0, 0.103889 },
{ 46, 1, 0, 23, 3, 0, 0.105114 },
{ 23, 1, 0, 23, 1, 0, 0.105682 },
{ 48, 0, 0, 16, 0, 2, 0.106406 },
{ 87, 0, 0, 29, 0, 2, 0.107545 },
{ 19, 0, 0, 19, 0, 0, 0.107778 },
{ 7, 0, 0, 7, 0, 0, 0.113333 },
{ 10, 0, 0, 10, 0, 0, 0.113333 },
{ 92, 0, 0, 23, 4, 0, 0.113636 },
{ 12, 0, 0, 12, 0, 0, 0.114062 },
{ 5, 1, 0, 5, 1, 0, 0.115000 },
{ 186, 1, 0, 31, 3, 2, 0.115344 },
{ 93, 1, 0, 31, 1, 2, 0.118125 },
{ 23, 0, 0, 23, 0, 0, 0.120909 },
{ 93, 0, 0, 31, 0, 2, 0.128250 },
{ 16, 0, 0, 16, 0, 0, 0.138750 },
{ 25, 0, 0, 25, 0, 0, 0.155000 },
{ 58, 1, 0, 29, 3, 0, 0.155714 },
{ 29, 1, 0, 29, 1, 0, 0.158214 },
{ 3, 1, 0, 3, 1, 0, 0.163125 },
{ 116, 0, 0, 29, 4, 0, 0.163690 },
{ 5, 0, 0, 5, 0, 0, 0.170000 },
{ 6, 0, 0, 6, 0, 0, 0.170000 },
{ 8, 0, 0, 8, 0, 0, 0.170000 },
{ 29, 0, 0, 29, 0, 0, 0.172857 },
{ 31, 1, 0, 31, 1, 0, 0.191333 },
{ 62, 1, 0, 31, 3, 0, 0.191750 },
{ 124, 0, 0, 31, 4, 0, 0.197917 },
{ 31, 0, 0, 31, 0, 0, 0.201667 },
{ 3, 0, 0, 3, 0, 0, 0.236250 },
{ 4, 0, 0, 4, 0, 0, 0.262500 },
{ 2, 1, 1, 2, 3, 1, 0.317187 },
{ 1, 1, 0, 1, 1, 0, 0.410000 },
{ 2, 0, 0, 2, 0, 0, 0.453333 },
};

static struct torctab_rec torctab4[] = {
{ 66, 1, 0, 11, 3, 2, 0.041344 },
{ 33, 1, 0, 11, 1, 2, 0.042750 },
{ 78, 1, 0, 13, 3, 2, 0.045781 },
{ 39, 1, 0, 13, 1, 2, 0.046875 },
{ 264, 1, 0, 11, 131, 2, 0.049043 },
{ 42, 1, 0, 14, 3, 2, 0.050000 },
{ 102, 1, 0, 17, 3, 2, 0.050508 },
{ 51, 1, 0, 17, 1, 2, 0.051094 },
{ 528, 1, 0, 11, 132, 2, 0.052891 },
{ 132, 0, 0, 11, 4, 2, 0.052969 },
{ 168, 1, 0, 14, 131, 2, 0.053965 },
{ 156, 0, 0, 13, 4, 2, 0.054948 },
{ 336, 1, 0, 14, 132, 2, 0.056120 },
{ 15, 1, 0, 15, 1, 0, 0.056719 },
{ 66, 0, 0, 11, 2, 2, 0.057000 },
{ 114, 1, 0, 19, 3, 2, 0.057812 },
{ 11, 1, 0, 11, 1, 0, 0.058000 },
{ 204, 0, 0, 17, 4, 2, 0.058203 },
{ 57, 1, 0, 19, 1, 2, 0.058542 },
{ 84, 0, 0, 14, 4, 2, 0.059375 },
{ 30, 1, 0, 15, 3, 0, 0.061406 },
{ 22, 1, 0, 11, 3, 0, 0.063000 },
{ 78, 0, 0, 13, 2, 2, 0.063542 },
{ 138, 1, 0, 23, 3, 2, 0.064815 },
{ 21, 1, 0, 7, 1, 2, 0.065000 },
{ 228, 0, 0, 19, 4, 2, 0.065104 },
{ 69, 1, 0, 23, 1, 2, 0.066477 },
{ 13, 1, 0, 13, 1, 0, 0.068750 },
{ 102, 0, 0, 17, 2, 2, 0.068906 },
{ 51, 0, 0, 17, 0, 2, 0.069141 },
{ 26, 1, 0, 13, 3, 0, 0.070625 },
{ 276, 0, 0, 23, 4, 2, 0.071236 },
{ 39, 0, 0, 13, 0, 2, 0.071250 },
{ 33, 0, 0, 11, 0, 2, 0.072750 },
{ 60, 0, 0, 15, 4, 0, 0.073828 },
{ 9, 1, 0, 9, 1, 0, 0.074097 },
{ 57, 0, 0, 19, 0, 2, 0.074583 },
{ 30, 0, 0, 15, 2, 0, 0.075625 },
{ 44, 0, 0, 11, 4, 0, 0.076667 },
{ 17, 1, 0, 17, 1, 0, 0.077188 },
{ 22, 0, 0, 11, 2, 0, 0.077333 },
{ 69, 0, 0, 23, 0, 2, 0.080114 },
{ 36, 0, 0, 18, 4, 0, 0.080208 },
{ 34, 1, 0, 17, 3, 0, 0.080312 },
{ 28, 0, 0, 14, 4, 0, 0.080556 },
{ 7, 1, 0, 7, 1, 0, 0.080833 },
{ 52, 0, 0, 13, 4, 0, 0.082292 },
{ 42, 0, 0, 14, 0, 2, 0.082500 },
{ 14, 1, 1, 14, 3, 1, 0.083333 },
{ 15, 0, 0, 15, 0, 0, 0.086250 },
{ 18, 1, 0, 9, 3, 0, 0.087083 },
{ 26, 0, 0, 13, 2, 0, 0.088889 },
{ 68, 0, 0, 17, 4, 0, 0.089583 },
{ 48, 1, 0, 16, 132, 2, 0.089844 },
{ 19, 1, 0, 19, 1, 0, 0.091389 },
{ 11, 0, 0, 11, 0, 0, 0.092000 },
{ 38, 1, 0, 19, 3, 0, 0.092917 },
{ 18, 0, 0, 9, 2, 0, 0.093704 },
{ 14, 1, 0, 7, 3, 0, 0.095000 },
{ 96, 1, 0, 16, 133, 2, 0.095391 },
{ 20, 0, 0, 10, 4, 0, 0.095833 },
{ 174, 1, 0, 29, 3, 2, 0.095893 },
{ 13, 0, 0, 13, 0, 0, 0.096667 },
{ 17, 0, 0, 17, 0, 0, 0.096875 },
{ 21, 0, 0, 21, 0, 0, 0.096875 },
{ 87, 1, 0, 29, 1, 2, 0.097366 },
{ 48, 0, 0, 16, 0, 2, 0.097969 },
{ 24, 1, 0, 12, 131, 0, 0.098789 },
{ 76, 0, 0, 19, 4, 0, 0.100926 },
{ 348, 0, 0, 29, 4, 2, 0.101116 },
{ 14, 0, 0, 14, 0, 0, 0.102222 },
{ 9, 0, 0, 9, 0, 0, 0.103889 },
{ 23, 1, 0, 23, 1, 0, 0.104545 },
{ 46, 1, 0, 23, 3, 0, 0.105682 },
{ 12, 0, 0, 12, 0, 0, 0.106250 },
{ 87, 0, 0, 29, 0, 2, 0.108348 },
{ 19, 0, 0, 19, 0, 0, 0.110000 },
{ 7, 0, 0, 7, 0, 0, 0.113333 },
{ 10, 0, 0, 10, 0, 0, 0.113333 },
{ 92, 0, 0, 23, 4, 0, 0.113826 },
{ 186, 1, 0, 31, 3, 2, 0.116094 },
{ 93, 1, 0, 31, 1, 2, 0.116813 },
{ 23, 0, 0, 23, 0, 0, 0.120909 },
{ 5, 1, 0, 5, 1, 0, 0.121250 },
{ 93, 0, 0, 31, 0, 2, 0.127625 },
{ 16, 0, 0, 16, 0, 0, 0.132917 },
{ 8, 0, 0, 8, 0, 0, 0.141667 },
{ 25, 0, 0, 25, 0, 0, 0.152500 },
{ 58, 1, 0, 29, 3, 0, 0.157946 },
{ 29, 1, 0, 29, 1, 0, 0.158393 },
{ 116, 0, 0, 29, 4, 0, 0.162946 },
{ 3, 1, 0, 3, 1, 0, 0.163125 },
{ 29, 0, 0, 29, 0, 0, 0.169286 },
{ 5, 0, 0, 5, 0, 0, 0.170000 },
{ 6, 0, 0, 6, 0, 0, 0.170000 },
{ 31, 1, 0, 31, 1, 0, 0.191333 },
{ 62, 1, 0, 31, 3, 0, 0.192083 },
{ 124, 0, 0, 31, 4, 0, 0.196389 },
{ 31, 0, 0, 31, 0, 0, 0.205000 },
{ 3, 0, 0, 3, 0, 0, 0.255000 },
{ 4, 0, 0, 4, 0, 0, 0.262500 },
{ 2, 1, 1, 2, 3, 1, 0.325000 },
{ 1, 1, 0, 1, 1, 0, 0.385000 },
{ 2, 0, 0, 2, 0, 0, 0.420000 },
};

#define TWIST_DOUBLE_RATIO              (9.0/16.0)

static long
torsion_constraint(struct torctab_rec *torctab, long ltorc, double tormod[], long n, long m)
{
  long i, b = -1;
  double rb = -1.;
  for (i = 0 ; i < ltorc ; i++)
  {
    struct torctab_rec *ti = torctab + i;
    if ( ! (n%ti->m) && ( !ti->fix2 || (n%(2*ti->m)) ) && ( ! ti->fix3 || (n%(3*ti->m)) ) )
      if ( n == m || ( ! (m%ti->m) && ( !ti->fix2 || (m%(2*ti->m)) ) && ( ! ti->fix3 || (m%(3*ti->m)) ) ) )
      {
        double ri = ti->rating*tormod[ti->N];
        if ( b < 0 || ri < rb ) {  b = i; rb = ri; }
      }
  }
  if (b < 0) pari_err_BUG("find_rating");
  return b;
}

/* p > 3 prime */
static void
best_torsion_constraint(ulong p, long t, int *ptwist, ulong *ptor)
{
  struct torctab_rec *torctab;
  double tormod[32];
  long ltorc, n1, n2, b, b1, b2, b12, i;

  switch(p % 12)
  {
    case 11:torctab = torctab1; ltorc = numberof(torctab1); break;
    case 5: torctab = torctab2; ltorc = numberof(torctab2); break;
    case 7: torctab = torctab3; ltorc = numberof(torctab3); break;
    default/*1*/: torctab = torctab4; ltorc = numberof(torctab4); break;
  }
  for ( i = 0 ; i < 32 ; i++ ) tormod[i] = 1.0;
  if (p%5  == 1) tormod[5] = tormod[10] = tormod[15] = 6.0/5.0;
  if (p%7  == 1) tormod[7] = tormod[14] = 8.0/7.0;
  if (p%11 == 1) tormod[11] = 12.0/11.0;
  if (p%13 == 1) tormod[13] = 14.0/13.0;
  if (p%17 == 1) tormod[17] = 18.0/17.0;
  if (p%19 == 1) tormod[19] = 20.0/19.0;
  if (p%23 == 1) tormod[23] = 24.0/23.0;
  if (p%29 == 1) tormod[29] = 30.0/29.0;
  if (p%31 == 1) tormod[31] = 32.0/31.0;

  n1 = p+1-t;
  n2 = p+1+t;
  b1  = torsion_constraint(torctab, ltorc, tormod, n1, n1);
  b2  = torsion_constraint(torctab, ltorc, tormod, n2, n2);
  b12 = torsion_constraint(torctab, ltorc, tormod, n1, n2);
  if (b1 > b2) {
    if (torctab[b2].rating / TWIST_DOUBLE_RATIO > torctab[b12].rating)
      *ptwist = 3;
    else
      *ptwist = 2;
  } else
    if (torctab[b1].rating / TWIST_DOUBLE_RATIO > torctab[b12].rating)
      *ptwist = 3;
    else
      *ptwist = 1;
  b = *ptwist ==1 ? b1: *ptwist ==2 ? b2: b12;
  *ptor = torctab[b].N;
}

/* This is Sutherland 2009 Algorithm 1.1 */
static long
find_j_inv_with_given_trace(
  ulong *j_t, norm_eqn_t ne, long rho_inv, long max_curves)
{
  pari_sp ltop = avma, av;
  long tested = 0, t = ne->t, batch_size, N0, N1, hasse_low, hasse_high, i;
  GEN n0, n1, A4, A6, tx, ty;
  ulong p = ne->p, pi = ne->pi, p1 = p + 1, a4, a6, m, N;
  int twist;

  if (p == 2 || p == 3) {
    if (t == 0) pari_err_BUG("find_j_inv_with_given_trace");
    *j_t = t; return 1;
  }

  N0 = (long)p1 - t; n0 = factoru(N0);
  N1 = (long)p1 + t; n1 = factoru(N1);
  best_torsion_constraint(p, t, &twist, &m);
  switch(twist)
  {
    case 1: N = N0; break;
    case 2: N = N1; break;
    default: N = p1;
  }

  /* Select batch size so that we have roughly a 50% chance of finding
   * a good curve in a batch. */
  batch_size = 1.0 + rho_inv / (2.0 * m);
  A4 = cgetg(batch_size + 1, t_VECSMALL);
  A6 = cgetg(batch_size + 1, t_VECSMALL);
  tx = cgetg(batch_size + 1, t_VECSMALL);
  ty = cgetg(batch_size + 1, t_VECSMALL);

  dbg_printf(2)("  Selected torsion constraint m = %lu and batch "
                "size = %ld\n", m, batch_size);
  hasse_bounds(&hasse_low, &hasse_high, p);
  av = avma;
  while (max_curves <= 0 || tested < max_curves)
  {
    GEN Pp1, Pt;
    random_curves_with_m_torsion((ulong *)(A4 + 1), (ulong *)(A6 + 1),
                                 (ulong *)(tx + 1), (ulong *)(ty + 1),
                                 batch_size, m, p, pi);
    Pp1 = random_FleV(A4, A6, p, pi);
    Pt = gcopy(Pp1);
    FleV_mulu_pre_inplace(Pp1, N, A4, p, pi);
    if (twist >= 3) FleV_mulu_pre_inplace(Pt, t, A4,  p, pi);
    for (i = 1; i <= batch_size; ++i) {
      ++tested;
      a4 = A4[i];
      a6 = A6[i]; if (a4 == 0 || a6 == 0) continue;

      if (( (twist >= 3 && mael(Pp1,i,1) == mael(Pt,i,1))
         || (twist < 3 && umael(Pp1,i,1) == p))
          && test_curve_order(ne, a4, a6, N0,N1, n0,n1, hasse_low,hasse_high)) {
        *j_t = Fl_ellj_pre(a4, a6, p, pi);
        return gc_long(ltop, tested);
      }
    }
    set_avma(av);
  }
  return gc_long(ltop, tested);
}

/* SECTION: Functions for dealing with polycyclic presentations. */

static GEN
next_generator(GEN DD, long D, ulong u, long filter, GEN *genred, long *P)
{
  pari_sp av = avma;
  ulong p = (ulong)*P;
  while (1)
  {
    p = unextprime(p + 1);
    if (p > LONG_MAX) pari_err_BUG("next_generator");
    if (kross(D, (long)p) != -1 && u % p != 0 && filter % p != 0)
    {
      GEN gen = primeform_u(DD, p);
      /* If gen is in the principal class, skip it */
      *genred = qfbred_i(gen);
      if (!equali1(gel(*genred,1))) { *P = (long)p; return gen; }
      set_avma(av);
    }
  }
}

INLINE long *
evec_ri_mutate(long r[], long i)
{ return r + (i * (i - 1) >> 1); }

INLINE const long *
evec_ri(const long r[], long i)
{ return r + (i * (i - 1) >> 1); }

/* Reduces evec e so that e[i] < n[i] (assume e[i] >= 0) using pcp(n,r,k).
 * No check for overflow, this could be an issue for large groups */
INLINE void
evec_reduce(long e[], const long n[], const long r[], long k)
{
  long i, j, q;
  const long *ri;
  if (!k) return;
  for (i = k - 1; i > 0; i--) {
    if (e[i] >= n[i]) {
      q = e[i] / n[i];
      ri = evec_ri(r, i);
      for (j = 0; j < i; j++) e[j] += q * ri[j];
      e[i] -= q * n[i];
    }
  }
  e[0] %= n[0];
}

/* Computes e3 = log(a^e1*a^e2) in terms of the given polycyclic
 * presentation (here a denotes the implicit vector of generators) */
INLINE void
evec_compose(long e3[],
  const long e1[], const long e2[], const long n[],const long r[], long k)
{
    long i;
    for (i = 0; i < k; i++) e3[i] = e1[i] + e2[i];
    evec_reduce(e3, n, r, k);
}

/* Converts an evec to an integer index corresponding to the
 * multi-radix representation of the evec with moduli corresponding to
 * the subgroup orders m[i] */
INLINE long
evec_to_index(const long e[], const long m[], long k)
{
  long i, index = e[0];
  for (i = 1; i < k; i++) index += e[i] * m[i - 1];
  return index;
}

INLINE void
evec_copy(long f[], const long e[], long k)
{
  long i;
  for (i = 0; i < k; ++i) f[i] = e[i];
}

INLINE void
evec_clear(long e[], long k)
{
  long i;
  for (i = 0; i < k; ++i) e[i] = 0;
}

/* e1 and e2 may overlap */
/* Note that this function is not very efficient because it does not know the
 * orders of the elements in the presentation, only the relative orders */
INLINE void
evec_inverse(long e2[], const long e1[], const long n[], const long r[], long k)
{
  pari_sp av = avma;
  long i, *e3, *e4;

  e3 = new_chunk(k);
  e4 = new_chunk(k);
  evec_clear(e4, k);
  evec_copy(e3, e1, k);
  /* We have e1 + e4 = e3 which we maintain throughout while making e1
   * the zero vector */
  for (i = k - 1; i >= 0; i--) if (e3[i])
  {
    e4[i] += n[i] - e3[i];
    evec_reduce(e4, n, r, k);
    e3[i] = n[i];
    evec_reduce(e3, n, r, k);
  }
  evec_copy(e2, e4, k);
  set_avma(av);
}

/* e1 and e2 may overlap */
/* This is a faster way to compute inverses, if the presentation
 * element orders are known (these are specified in the array o, the
 * array n holds the relative orders) */
INLINE void
evec_inverse_o(
  long e2[],
  const long e1[], const long n[], const long o[], const long r[], long k)
{
  long j;
  for (j = 0; j < k; j++) e2[j] = (e1[j] ? o[j] - e1[j] : 0);
  evec_reduce(e2, n, r, k);
}

/* Computes the order of the group element a^e using the pcp (n,r,k) */
INLINE long
evec_order(const long e[], const long n[], const long r[], long k)
{
  pari_sp av = avma;
  long *f = new_chunk(k);
  long i, j, o, m;

  evec_copy(f, e, k);
  for (o = 1, i = k - 1; i >= 0; i--) if (f[i])
  {
    m = n[i] / ugcd(f[i], n[i]);
    for (j = 0; j < k; j++) f[j] *= m;
    evec_reduce(f, n, r, k);
    o *= m;
  }
  return gc_long(av,o);
}

/* Computes orders o[] for each generator using relative orders n[]
 * and power relations r[] */
INLINE void
evec_orders(long o[], const long n[], const long r[], long k)
{
  pari_sp av = avma;
  long i, *e = new_chunk(k);

  evec_clear(e, k);
  for (i = 0; i < k; i++) {
    e[i] = 1;
    if (i) e[i - 1] = 0;
    o[i] = evec_order(e, n, r, k);
  }
  set_avma(av);
}

INLINE int
evec_equal(const long e1[], const long e2[], long k)
{
  long j;
  for (j = 0; j < k; ++j)
    if (e1[j] != e2[j]) break;
  return j == k;
}

INLINE void
index_to_evec(long e[], long index, const long m[], long k)
{
  long i;
  for (i = k - 1; i > 0; --i) {
    e[i] = index / m[i - 1];
    index -= e[i] * m[i - 1];
  }
  e[0] = index;
}

INLINE void
evec_n_to_m(long m[], const long n[], long k)
{
  long i;
  m[0] = n[0];
  for (i = 1; i < k; ++i) m[i] = m[i - 1] * n[i];
}

/* Based on logfac() in Sutherland's classpoly package.
 * Ramanujan approximation to log(n!), accurate to O(1/n^3) */
INLINE double
logfac(long n)
{
  const double HALFLOGPI = 0.57236494292470008707171367567653;
  return n * log((double) n) - (double) n +
    log((double) n * (1.0 + 4.0 * n * (1.0 + 2.0 * n))) / 6.0 + HALFLOGPI;
}

/* This is based on Sutherland 2009, Lemma 8 (p31). */
static double
upper_bound_on_classpoly_coeffs(long D, long h, GEN qfinorms)
{
  double B, logbinom, lnMk, lnMh = 0, C = 2114.567, t = M_PI * sqrt((double)-D);
  ulong maxak = 0;
  long k, m;

  for (k = 1, B = 0.0; k <= h; ++k)
  {
    ulong ak = uel(qfinorms, k);
    double tk = t / ak;
    lnMk = tk + log(1.0 + C * exp(-tk));
    B += lnMk;
    if (ak > maxak) { maxak = ak; lnMh = lnMk; }
  }
  m = floor((h + 1)/(exp(lnMh) + 1.0));
  /* log(binom(h, m)); 0 unless D <= -1579751 */
  logbinom = (m > 0 && m < h)? logfac(h) - logfac(m) - logfac(h - m): 0;
  return (B + logbinom - m * lnMh) * (1 / M_LN2) + 2.0;
}

INLINE long
distinct_inverses(const long f[], const long ef[], const long ei[],
  const long n[], const long o[], const long r[], long k, long L0, long i)
{
  pari_sp av = avma;
  long j, *e2, *e3;

  if ( ! ef[i] || (L0 && ef[0])) return 0;
  for (j = i + 1; j < k; ++j)
    if (ef[j]) break;
  if (j < k) return 0;

  e2 = new_chunk(k);
  evec_copy(e2, ef, i);
  e2[i] = o[i] - ef[i];
  for (j = i + 1; j < k; ++j) e2[j] = 0;
  evec_reduce(e2, n, r, k);

  if (evec_equal(ef, e2, k)) return gc_long(av,0);

  e3 = new_chunk(k);
  evec_inverse_o(e3, ef, n, o, r, k);
  if (evec_equal(e2, e3, k)) return gc_long(av,0);

  if (f) {
    evec_compose(e3, f, ei, n, r, k);
    if (evec_equal(e2, e3, k)) return gc_long(av,0);

    evec_inverse_o(e3, e3, n, o, r, k);
    if (evec_equal(e2, e3, k)) return gc_long(av,0);
  }
  return gc_long(av,1);
}

INLINE long
next_prime_evec(long *qq, long f[], const long m[], long k,
  hashtable *tbl, long D, GEN DD, long u, long lvl, long ubound)
{
  pari_sp av = avma;
  hashentry *he;
  GEN P;
  long idx, q = *qq;

  do q = unextprime(q + 1);
  while (!(u % q) || kross(D, q) == -1 || !(lvl % q) || !(D % (q * q)));
  if (q > ubound) return 0;
  *qq = q;

  /* Get evec f corresponding to q */
  P = qfbred_i(primeform_u(DD, q));
  he = hash_search(tbl, P);
  if (!he) pari_err_BUG("next_prime_evec");
  idx = (long)he->val;
  index_to_evec(f, idx, m, k);
  return gc_long(av,1);
}

/* Return 1 on success, 0 on failure. */
static int
orient_pcp(GEN G, long *ni, long D, long u, hashtable *tbl)
{
  pari_sp av = avma;
  /* 199 seems to suffice, but can be increased if necessary */
  enum { MAX_ORIENT_P = 199 };
  const long *L = pcp_get_L(G), *n = pcp_get_n(G), *r = pcp_get_r(G), *m = pcp_get_m(G), *o = pcp_get_o(G);
  long i, *ps = pcp_get_orient_p(G), *qs = pcp_get_orient_q(G), *reps = pcp_get_orient_reps(G);
  long *ef, *e, *ei, *f, k = pcp_get_k(G), lvl = modinv_level(pcp_get_inv(G));
  GEN DD = stoi(D);
  long L0 = pcp_get_L0(G);

  memset(ps, 0, k * sizeof(long));
  memset(qs, 0, k * sizeof(long));
  memset(reps, 0, k * k * sizeof(long));

  for (i = 0; i < k; ++i) { ps[i] = -1; if (o[i] > 2) break; }
  for (++i; i < k; ++i) ps[i] = (o[i] > 2) ? 0 : -1; /* ps[i] = -!(o[i] > 2); */

  e = new_chunk(k);
  ei = new_chunk(k);
  f = new_chunk(k);

  for (i = 0; i < k; ++i) {
    long p;
    if (ps[i]) continue;
    p = L[i];
    ef = &reps[i * k];
    while (!ps[i]) {
      if (!next_prime_evec(&p, ef, m, k, tbl, D, DD, u, lvl, MAX_ORIENT_P))
        break;
      evec_inverse_o(ei, ef, n, o, r, k);
      if (!distinct_inverses(NULL, ef, ei, n, o, r, k, L0, i)) continue;
      ps[i] = p;
      qs[i] = 1;
    }
    if (ps[i]) continue;

    p = unextprime(L[i] + 1);
    while (!ps[i]) {
      long q;

      if (!next_prime_evec(&p, e, m, k, tbl, D, DD, u, lvl, MAX_ORIENT_P))
        break;
      evec_inverse_o(ei, e, n, o, r, k);

      q = L[i];
      while (!qs[i]) {
        if (!next_prime_evec(&q, f, m, k, tbl, D, DD, u, lvl, p - 1)) break;
        evec_compose(ef, e, f, n, r, k);
        if (!distinct_inverses(f, ef, ei, n, o, r, k, L0, i)) continue;
        ps[i] = p;
        qs[i] = q;
      }
    }
    if (!ps[i]) return 0;
  }
  if (ni) {
    GEN N = qfb_nform(D, *ni);
    hashentry *he = hash_search(tbl, N);
    if (!he) pari_err_BUG("orient_pcp");
    *ni = (long)he->val;
  }
  return gc_bool(av,1);
}

/* We must avoid situations where L_i^{+/-2} = L_j^2 (or = L_0*L_j^2
 * if ell0 flag is set), with |L_i| = |L_j| = 4 (or have 4th powers in
 * <L0> but not 2nd powers in <L0>) and j < i */
/* These cases cause problems when enumerating roots via gcds */
/* returns the index of the first bad generator, or -1 if no bad
 * generators are found */
static long
classgp_pcp_check_generators(const long *n, long *r, long k, long L0)
{
  pari_sp av = avma;
  long *e1, i, i0, j, s;
  const long *ei;

  s = !!L0;
  e1 = new_chunk(k);

  for (i = s + 1; i < k; i++) {
    if (n[i] != 2) continue;
    ei = evec_ri(r, i);
    for (j = s; j < i; j++)
      if (ei[j]) break;
    if (j == i) continue;
    for (i0 = s; i0 < i; i0++) {
      if ((4 % n[i0])) continue;
      evec_clear(e1, k);
      e1[i0] = 4;
      evec_reduce(e1, n, r, k);
      for (j = s; j < i; j++)
        if (e1[j]) break;
      if (j < i) continue; /* L_i0^4 is not trivial or in <L_0> */
      evec_clear(e1, k);
      e1[i0] = 2;
      evec_reduce(e1, n, r, k); /* compute L_i0^2 */
      for (j = s; j < i; j++)
        if (e1[j] != ei[j]) break;
      if (j == i) return i;
      evec_inverse(e1, e1, n, r, k); /* compute L_i0^{-2} */
      for (j = s; j < i; j++)
        if (e1[j] != ei[j]) break;
      if (j == i) return i;
    }
  }
  return gc_long(av,-1);
}

/* This is Sutherland 2009, Algorithm 2.2 (p16). */
static GEN
classgp_make_pcp(double *height, long *ni, long h, long D, long D0, ulong u,
  GEN Pu, GEN Eu, long inv, long Lfilter, long orient)
{
  const long MAX_GENS = 16, lvl = modinv_level(inv);
  pari_sp av, av2;
  long curr_p, nelts, L0, enum_cnt, GLfilter, i, k, L1, L2;
  GEN G, DD, ident, T, v, L, m, n, o, r, orient_p, orient_q, orient_reps;
  hashtable *tbl;

  L = zero_zv(MAX_GENS);
  m = zero_zv(MAX_GENS);
  n = zero_zv(MAX_GENS);
  o = zero_zv(MAX_GENS);
  r = zero_zv(MAX_GENS * (MAX_GENS-1) / 2);
  orient_p = zero_zv(MAX_GENS);
  orient_q = zero_zv(MAX_GENS);
  orient_reps = zero_zv(MAX_GENS*MAX_GENS);
  G = mkvec5(mkvec5(L, n, o, m, r), mkvecsmall3(0, 0, 0),
             mkvec3(orient_p, orient_q, orient_reps),
             mkvecsmall5(h, inv, D, D0, u), mkmat2(Pu, Eu));
  av = avma;
  if (h == 1)
  {
    *height = upper_bound_on_classpoly_coeffs(D, h, mkvecsmall(1));
    return gc_const(av, G); /* no need to set *ni when h = 1 */
  }
  if (!modinv_is_double_eta(inv) || !modinv_ramified(D, inv, &L0)) L0 = 0;
  enum_cnt = h / (1 + !!L0);
  GLfilter = ulcm(Lfilter, lvl);
  DD = stoi(D); av2 = avma;
  while (1) {
    k = 0;
    /* Hash table has an imaginary QFB as a key and the index of that
     * form in T as its value */
    tbl = hash_create(h, (ulong(*)(void*)) hash_GEN,
                         (int(*)(void*,void*))&gidentical, 1);
    ident = primeform(DD, gen_1);
    hash_insert(tbl, ident, (void*)0);

    T = vectrunc_init(h + 1);
    vectrunc_append(T, ident);
    nelts = 1;
    curr_p = 1;

    while (nelts < h) {
      GEN gamma_i = NULL, beta;
      hashentry *e;
      long N = lg(T)-1, ri = 1;

      if (k == MAX_GENS) pari_err_IMPL("classgp_pcp");

      if (nelts == 1 && L0) {
        curr_p = L0;
        gamma_i = qfb_nform(D, curr_p);
        beta = qfbred_i(gamma_i);
        if (equali1(gel(beta, 1))) { curr_p = 1; gamma_i  = NULL; }
      }
      if (!gamma_i)
        gamma_i = next_generator(DD, D, u, GLfilter, &beta, &curr_p);
      while ((e = hash_search(tbl, beta)) == NULL) {
        long j;
        for (j = 1; j <= N; ++j) {
          GEN t = qfbcomp_i(beta, gel(T, j));
          vectrunc_append(T, t);
          hash_insert(tbl, t, (void*)tbl->nb);
        }
        beta = qfbcomp_i(beta, gamma_i);
        ++ri;
      }
      if (ri > 1) {
        long j, si;
        L[k+1] = curr_p;
        n[k+1] = ri;
        nelts *= ri;

        /* This is to reset the curr_p counter when we have L0 != 0
         * in the first position of L. */
        if (curr_p == L0) curr_p = 1;

        N = 1;
        si = (long)e->val;
        for (j = 1; j <= k; ++j) {
          evec_ri_mutate(r, k)[j] = (si / N) % n[j];
          N *= n[j];
        }
        ++k;
      }
    }

    if ((i = classgp_pcp_check_generators(n+1, r+1, k, L0)) < 0) {
      mael(G,2,1) = L0;
      mael(G,2,2) = k;
      mael(G,2,3) = enum_cnt;
      evec_orders(o+1, n+1, r+1, k);
      evec_n_to_m(m+1, n+1, k);
      if (!orient || orient_pcp(G, ni, D, u, tbl)) break;
      GLfilter *= L[1];
    } else {
      GLfilter = umuluu_or_0(GLfilter, L[i+1]);
      if (!GLfilter) pari_err_IMPL("classgp_pcp");
    }
    set_avma(av2);
  }
  v = cgetg(h + 1, t_VECSMALL);
  v[1] = 1;
  for (i = 2; i <= h; ++i) uel(v,i) = itou(gmael(T,i,1));
  *height = upper_bound_on_classpoly_coeffs(D, enum_cnt, v);

  /* The norms of the last one or two generators. */
  L1 = L[k];
  L2 = k > 1 ? L[k - 1] : 1;
  /* 4 * L1^2 * L2^2 must fit in a ulong */
  if (2 * (1 + log2(L1) + log2(L2)) >= BITS_IN_LONG)
    pari_err_IMPL("classgp_pcp");
  if (L0 && (L[1] != L0 || o[1] != 2)) pari_err_BUG("classgp_pcp");
  return gc_const(av, G);
}

/* SECTION: Functions for calculating class polynomials. */
static const long SMALL_PRIMES[11] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };

/* Encodes prime divisors of smooth integers <= 1200
  P = primes(11); V = vector(31, i, vecsearch(P, i));
  vfactor(v) =
  { if (v == 1, return (0));
    my(q = factor(v)[,1]);
    if (vecmax(q) > 31, return (-1)); \\ not smooth
    sum(i = 1, #P, 1 << (V[q[i]]-1)); }
  vector(1200, v, vfactor(v)) */
static const long
SMOOTH_INTS[] = { -1, /* 0 */
0,1,2,1,4,3,8,1,2,5,16,3,32,9,6,1,64,3,128,5,10,17,256,3,4,33,2,9,512,7,1024,1,
18,65,12,3,-1,129,34,5,-1,11,-1,17,6,257,-1,3,8,5,66,33,-1,3,20,9,130,513,-1,7,
-1,1025,10,1,36,19,-1,65,258,13,-1,3,-1,-1,6,129,24,35,-1,5,2,-1,-1,11,68,-1,
514,17,-1,7,40,257,1026,-1,132,3,-1,9,18,5,-1,67,-1,33,14,-1,-1,3,-1,21,-1,9,-1,
131,260,513,34,-1,72,7,16,-1,-1,1025,4,11,-1,1,-1,37,-1,19,136,-1,6,65,-1,259,
-1,13,-1,-1,48,3,516,-1,10,-1,-1,7,-1,129,66,25,1028,35,-1,-1,-1,5,264,3,-1,-1,
22,-1,-1,11,32,69,130,-1,-1,515,12,17,-1,-1,-1,7,-1,41,-1,257,-1,1027,80,-1,10,
133,-1,3,-1,-1,38,9,-1,19,-1,5,-1,-1,520,67,-1,-1,258,33,144,15,-1,-1,-1,-1,-1,
3,1032,-1,-1,21,96,-1,-1,9,6,-1,-1,131,-1,261,26,513,-1,35,-1,-1,-1,73,-1,7,-1,
17,2,-1,12,-1,160,1025,-1,5,-1,11,272,-1,70,1,-1,-1,-1,37,514,-1,-1,19,-1,137,
-1,-1,-1,7,-1,65,42,-1,20,259,-1,-1,1026,13,-1,-1,-1,-1,134,49,-1,3,64,517,-1,
-1,-1,11,-1,-1,18,-1,288,7,-1,-1,-1,129,-1,67,-1,25,-1,1029,-1,35,-1,-1,14,-1,
-1,-1,528,5,-1,265,192,3,36,-1,-1,-1,-1,23,-1,-1,-1,-1,-1,11,-1,33,-1,69,1040,
131,8,-1,262,-1,-1,515,-1,13,34,17,-1,-1,-1,-1,74,-1,-1,7,128,-1,18,41,-1,-1,-1,
257,-1,-1,-1,1027,-1,81,6,-1,544,11,-1,133,-1,-1,-1,3,28,-1,-1,-1,-1,39,320,9,
-1,-1,-1,19,-1,-1,138,5,-1,-1,1056,-1,6,521,-1,67,-1,-1,-1,-1,-1,259,-1,33,-1,
145,-1,15,-1,-1,-1,-1,68,-1,-1,-1,50,-1,-1,3,-1,1033,518,-1,384,-1,-1,21,10,97,
-1,-1,-1,-1,-1,9,-1,7,-1,-1,-1,-1,44,131,-1,-1,66,261,-1,27,-1,513,1030,-1,-1,
35,-1,-1,-1,-1,-1,-1,132,73,-1,-1,-1,7,-1,-1,266,17,-1,3,-1,-1,-1,13,-1,-1,576,
161,22,1025,-1,-1,-1,5,-1,-1,-1,11,-1,273,34,-1,-1,71,-1,1,130,-1,-1,-1,-1,-1,
-1,37,-1,515,-1,-1,14,-1,1088,19,256,-1,-1,137,-1,-1,-1,-1,-1,-1,24,7,-1,-1,-1,
65,-1,43,-1,-1,-1,21,640,259,-1,-1,-1,-1,-1,1027,-1,13,82,-1,-1,-1,-1,-1,10,-1,
-1,135,-1,49,-1,-1,260,3,-1,65,-1,517,-1,-1,-1,-1,38,-1,-1,11,1152,-1,-1,-1,-1,
19,76,-1,-1,289,-1,7,-1,-1,-1,-1,20,-1,-1,129,522,-1,-1,67,-1,-1,-1,25,-1,-1,-1,
1029,258,-1,-1,35,4,-1,146,-1,-1,15,-1,-1,-1,-1,-1,-1,40,529,-1,5,-1,-1,-1,265,
-1,193,-1,3,-1,37,1034,-1,-1,-1,-1,-1,-1,-1,-1,23,-1,-1,98,-1,140,-1,768,-1,-1,
-1,-1,11,-1,-1,6,33,-1,-1,-1,69,-1,1041,-1,131,-1,9,-1,-1,-1,263,-1,-1,26,-1,-1,
515,-1,-1,-1,13,-1,35,-1,17,-1,-1,-1,-1,-1,-1,-1,-1,1280,75,52,-1,-1,-1,-1,7,-1,
129,-1,-1,516,19,-1,41,2,-1,-1,-1,-1,-1,14,257,-1,-1,-1,-1,162,-1,-1,1027,-1,-1,
-1,81,-1,7,-1,-1,-1,545,-1,11,-1,-1,274,133,-1,-1,-1,-1,70,-1,-1,3,-1,29,-1,-1,
-1,-1,1028,-1,-1,-1,-1,39,-1,321,514,9,-1,-1,-1,-1,-1,-1,-1,19,-1,-1,-1,-1,-1,
139,-1,5,-1,-1,-1,-1,268,1057,-1,-1,-1,7,-1,521,-1,-1,-1,67,-1,-1,42,-1,-1,-1,
-1,-1,22,-1,-1,259,-1,-1,-1,33,72,-1,-1,145,1026,-1,-1,15,512,-1,-1,-1,36,-1,24,
-1,-1,69,-1,-1,-1,-1,134,-1,-1,51,-1,-1,-1,-1,-1,3,-1,-1,66,1033,-1,519,-1,-1,
-1,385,12,-1,-1,-1,-1,21,-1,11,-1,97,-1,-1,-1,-1,-1,-1,18,-1,-1,-1,-1,9,290,-1,
1536,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,45,-1,131,-1,-1,-1,-1,-1,67,-1,261,-1,-1,-1,
27,-1,-1,-1,513,-1,1031,136,-1,-1,-1,84,35,-1,-1,-1,-1,-1,-1,-1,-1,14,-1,-1,-1,
-1,133,-1,73,-1,-1,-1,-1,530,-1,-1,7,1024,-1,-1,-1,-1,267,-1,17,194,-1,-1,3,-1,
-1,38,-1,-1,-1,-1,13,-1,-1,-1,-1,-1,577,-1,161,-1,23,-1,1025,-1,-1,-1,-1,-1,-1,
-1,5,56,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,273,-1,35,524,-1,-1,-1,-1,71,-1,-1,1042,1,
-1,131,-1,-1,10,-1,-1,-1,-1,-1,262,-1,-1,-1,-1,37,-1,-1,-1,515,148,-1,-1,-1,-1,
15,-1,-1,34,1089,-1,19,-1,257,-1,-1,-1,-1,-1,137,-1,-1,-1,-1,-1,-1,74,-1,-1,-1,
-1,-1,-1,25,-1,7,-1,-1,130,-1,1036,-1,-1,65,18,-1,-1,43,-1,-1,-1,-1,-1,-1,-1,21,
-1,641,-1,259,100,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1027,-1,-1,-1,13,-1,83,-1,-1,6,
-1,264,-1,-1,-1,546,-1,-1,11,-1,-1,-1,-1,-1,135,-1,-1,-1,49,-1,-1,-1,-1,-1,261,
-1,3,-1,-1,30,65,-1,-1,-1,517,-1,-1,-1,-1,-1,-1,-1,-1,-1,39,-1,-1,322,-1,-1,11,
-1,1153,-1,-1,-1,-1,40,-1,-1,-1,-1,19,-1,77,-1,-1,-1,-1,-1,289,138,-1,-1,7
};

/* Upper bound for H(v^2 d) / H(d) = \prod_{p | v} (p + 1) / (p - 1)
 * We actually store ceil(128 * bound)
P = primes(11); V = vector(31, i, vecsearch(P, i));
hbound(v) =
{ my(q = factor(v)[,1]); if (q && vecmax(q) > 31, return (0));
   ceil(prod(i = 1, #q, (q[i] + 1)/(q[i] - 1), 1.) * 128); }
vector(1200, v, hbound(v))
*/
static const long HURWITZ_RATIO[] = {
128,384,256,384,192,768,171,384,256,576,154,768,150,512,384,
384,144,768,143,576,342,461,140,768,192,448,256,512,138,1152,
137,384,308,432,256,768,0,427,299,576,0,1024,0,461,384,419,0,
768,171,576,288,448,0,768,231,512,285,412,0,1152,0,410,342,
384,224,922,0,432,280,768,0,768,0,0,384,427,205,896,0,576,
256,0,0,1024,216,0,275,461,0,1152,200,419,274,0,214,768,0,
512,308,576,0,864,0,448,512,0,0,768,0,692,0,512,0,854,210,
412,299,0,192,1152,154,0,0,410,192,1024,0,384,0,672,0,922,
190,0,384,432,0,838,0,768,0,0,180,768,206,0,342,0,0,1152,0,
427,288,615,205,896,0,0,0,576,187,768,0,0,461,0,0,1024,150,
648,285,0,0,823,256,461,0,0,0,1152,0,598,0,419,0,820,173,0,
342,640,0,768,0,0,448,512,0,922,0,576,0,0,183,864,0,0,280,
448,171,1536,0,0,0,0,0,768,183,0,0,692,168,0,0,512,384,0,
0,854,0,629,410,412,0,896,0,0,0,576,0,1152,0,461,256,0,256,
0,166,410,0,576,0,1024,168,0,432,384,0,0,0,672,275,0,0,922,
0,569,0,0,0,1152,0,432,399,0,231,838,0,0,274,768,0,0,0,0,
427,538,0,768,144,618,0,0,0,1024,0,0,308,0,163,1152,0,0,0,
427,0,864,0,615,0,615,0,896,0,0,512,0,0,0,165,576,0,559,
160,768,224,0,0,0,0,1383,0,0,0,0,0,1024,0,448,0,648,164,
854,171,0,419,0,0,823,0,768,299,461,0,0,0,0,384,0,0,1152,
143,0,308,598,0,0,0,419,0,0,0,820,0,519,384,0,160,1024,0,
640,0,0,0,768,308,0,0,0,0,1344,158,512,0,0,0,922,0,0,380,
576,0,0,160,0,384,549,0,864,0,0,0,0,0,838,0,448,0,512,0,
1536,0,0,0,0,216,0,0,0,359,0,0,768,0,547,412,0,156,0,0,
692,342,504,0,0,0,0,0,512,0,1152,0,0,0,0,299,854,0,0,288,
629,0,1229,0,412,410,0,0,896,0,0,0,0,0,0,214,576,0,0,0,
1152,0,0,373,461,0,768,0,0,0,768,0,0,155,498,461,410,0,0,
0,576,0,0,0,1024,0,503,299,0,0,1296,0,384,285,0,0,0,0,0,
0,672,0,823,0,0,512,0,154,922,140,0,0,569,0,0,0,0,0,0,
205,1152,0,0,0,432,0,1195,0,0,0,692,153,838,0,0,0,0,0,820,
0,768,346,0,0,0,0,0,342,0,0,1280,0,538,0,0,210,768,0,432,
0,618,0,0,0,0,448,0,0,1024,152,0,0,0,0,922,288,0,0,489,0,
1152,0,0,0,0,231,0,0,427,366,0,0,864,0,0,0,615,0,0,0,615,
280,0,0,896,192,0,342,0,0,1536,0,0,0,0,0,0,200,494,0,576,
0,0,0,559,0,480,0,768,0,672,365,0,0,0,0,0,0,0,0,1383,0,
0,336,0,285,0,150,0,0,0,0,1024,0,0,384,448,0,0,0,648,0,
492,0,854,0,512,0,0,0,1257,0,0,410,0,0,823,0,0,0,768,0,
896,0,461,0,0,0,0,0,0,0,0,149,1152,269,0,0,0,0,1152,0,
427,0,0,206,922,0,598,256,0,0,0,0,0,512,419,0,0,0,0,332,
0,0,820,0,0,0,519,0,1152,0,0,0,480,0,1024,0,0,336,640,0,
0,0,0,432,0,0,768,0,922,0,0,0,0,205,0,0,0,0,1344,0,472,
275,512,0,0,0,0,0,0,0,922,0,0,0,0,0,1138,0,576,0,0,0,0,
280,478,0,0,0,1152,0,549,0,0,0,864,0,0,399,0,0,0,0,0,461,
0,0,838,0,0,0,448,192,0,0,512,274,0,0,1536,138,0,0,0,224,
0,205,0,0,648,0,0,0,0,427,0,0,1076,0,0,0,0,0,768,0,0,
288,547,0,1235,0,0,0,466,256,0,0,0,0,692,0,1024,0,504,0,0,
0,0,0,0,308,0,0,0,0,512,326,0,147,1152,0,0,0,0,0,0,0,0,
0,896,0,854,0,0,0,0,0,864,0,629,0,0,0,1229,0,0,0,412,0,
1229,190,0,0,0,260,896,0,0,0,0,0,0,0,0,512,0,0,0,0,640,
0,576,0,0,0,0,330,0,0,1152,137,0,0,0,0,1118,0,461,320,0,
0,768,0,0,448,0,0,0,0,768,0,0,0,0,0,463,0,498,0,1383,0,
410,0,0,0,0,0,0,0,576,239,0,0,0,0,0,0,1024,0,0,0,503,0,
896,275,0,0,0,0,1296,0,0,328,384,0,854,0,0,342,0,0,0,0,0,
419,0,0,0,0,672,0,0,0,823,256,0,0,0,0,1536,0,0,299,461,0,
922,0,419,0,0,0,0,0,569,0,0,0,0,0,0,384,0,0,0,0,0,0,
615,0,1152,0,0,285,0,274,0,0,432,308,0,0,1195,0,0,0,0,0,
0,0,692,0,458,0,838,252,0,0,0,0,0,0,0,0,0,0,820,0,0,0,
768,0,1037,0,0,384,0,187,0,0,0,320,0,0,1024,0,0,0,0,0,
1280,0,0,0,538,0,0,0,0,0,629,0,768,0,0,615,432,0,0,0,618,
0,0,0,0,0,0,0,0,0,1344,0,0,315,0,0,1024,0,456,0,0,0,0,
200,0,0,0,0,922,0,864,0,0,0,0,0,489,380,0,0,1152
};

/* Hurwitz class number of Df^2, D < 0; h = classno(D), Faf = factoru(f) */
static double
hclassno_wrapper(long h, long D, GEN Faf)
{
  pari_sp av;
  if (lg(gel(Faf,1)) == 1) switch(D)
  {
    case -3: return 1. / 3;
    case -4: return 1. / 2;
    default: return (double)h;
  }
  av = avma;
  return (double)gc_long(av,  h * uhclassnoF_fact(Faf, D));
}

/* return factor(u*v) */
static GEN
factor_uv(GEN fau, ulong v, ulong vfactors)
{
  GEN P, E;
  long i;
  if (!vfactors) return fau;
  P = gel(fau,1);
  E = gel(fau,2);
  for (i = 0; vfactors; i++, vfactors >>= 1)
    if (vfactors & 1UL)
    {
      long p = SMALL_PRIMES[i];
      P = vecsmall_append(P, p);
      E = vecsmall_append(E, u_lvalrem(v, p, &v));
      if (v == 1) break;
    }
  return famatsmall_reduce(mkmat2(P, E));
}

/* This is Sutherland 2009, Algorithm 2.1 (p8); delta > 0 */
static GEN
select_classpoly_prime_pool(double min_bits, double delta, GEN G)
{ /* Sutherland defines V_MAX to be 1200 without saying why */
  const long V_MAX = 1200;
  double bits = 0.0, hurwitz, z;
  ulong t_size_lim, d = (ulong)-pcp_get_D(G);
  long ires, inv = pcp_get_inv(G);
  GEN fau = pcp_get_fau(G);
  GEN res, t_min; /* t_min[v] = lower bound for the t we look at for that v */
#ifdef LONG_IS_64BIT
  long L = pcp_get_L(G)[!!pcp_get_L0(G)];
#endif

  hurwitz = hclassno_wrapper(pcp_get_h(G), pcp_get_D0(G), fau);

  res = cgetg(128+1, t_VEC);
  ires = 1;
  /* Initialise t_min to be all 2's.  This avoids trace 0 and trace 1 curves */
  t_min = const_vecsmall(V_MAX-1, 2);

  /* maximum possible trace = sqrt(2^BIL + D) */
  t_size_lim = 2.0 * sqrt((double)((1UL << (BITS_IN_LONG - 2)) - (d >> 2)));

  for (z = d / (2.0 * hurwitz); ; z *= delta + 1.0)
  {
    double v_bound_aux = 4.0 * z * hurwitz; /* = 4 z H(d) */
    ulong v;
    dbg_printf(1)("z = %.2f\n", z);
    for (v = 1; v < V_MAX; v++)
#ifdef LONG_IS_64BIT
      if (L<=2 || v%L)
#endif
      {
        ulong p, t, t_max, vfactors, v2d, vd;
        double hurwitz_ratio_bound, max_p, H;
        long ires0;
        GEN faw;

        if ((long)(vfactors = SMOOTH_INTS[v]) < 0) continue;
        hurwitz_ratio_bound = HURWITZ_RATIO[v] / 128.0;
        vd = v * d;
        if (vd >= v_bound_aux * hurwitz_ratio_bound) break;
        v2d = v * vd;
        faw = factor_uv(fau, v, vfactors);
        H = hclassno_wrapper(pcp_get_h(G), pcp_get_D0(G), faw);
        /* t <= 2 sqrt(p) and p <= z H(v^2 d) and
         *   H(v^2 d) < v H(d) \prod_{p | v} (p+1)/(p-1)
         * This last term is v * hurwitz * hurwitz_ratio_bound. */
        max_p = z * v * hurwitz * hurwitz_ratio_bound;
        t_max = 2.0 * sqrt(mindd((1UL<<(BITS_IN_LONG-2)) - (v2d>>2), max_p));
        t = t_min[v]; if ((t & 1) != (v2d & 1)) t++;
        p = (t * t + v2d) >> 2;
        ires0 = ires;
        for (; t <= t_max; p += t+1, t += 2) /* 4p = t^2 + v^2*d */
          if (modinv_good_prime(inv,p) && uisprime(p))
          {
            if (ires == lg(res)) res = vec_lengthen(res, lg(res) << 1);
            gel(res, ires++) = mkvec2(mkvecsmall5(p,t,v,(long)(p/H),vfactors),
                faw);
            bits += log2(p);
          }
        t_min[v] = t;

        if (ires - ires0) {
          dbg_printf(2)("  Found %lu primes for v = %lu.\n", ires - ires0, v);
        }
        if (bits > min_bits) {
          dbg_printf(1)("Found %ld primes; total size %.2f bits.\n", ires-1,bits);
          setlg(res, ires); return res;
        }
      }
    if (uel(t_min,1) >= t_size_lim) {
      /* exhausted all solutions that fit in ulong */
      char *err = stack_sprintf("class polynomial of discriminant %ld", pcp_get_D(G));
      pari_err(e_ARCH, err);
    }
  }
}

static int
primecmp(void *data, GEN v1, GEN v2)
{ (void)data; return cmpss(gel(v1,1)[4], gel(v2,1)[4]); }

static long
height_margin(long inv, long D)
{
  (void)D;
  /* NB: avs just uses a height margin of 256 for everyone and everything. */
  if (inv == INV_F) return 64;  /* checked for discriminants up to -350000 */
  if (inv == INV_G2) return 5;
  if (inv != INV_J) return 256; /* TODO: This should be made more accurate */
  return 0;
}

static GEN
select_classpoly_primes(ulong *vfactors, ulong *biggest_v, double delta,
  GEN G, double height)
{
  const long k = 2;
  pari_sp av = avma;
  long i, s, D = pcp_get_D(G), inv = pcp_get_inv(G);
  ulong biggest_p;
  double prime_bits, min_prime_bits, b;
  GEN prime_pool;

  s = modinv_height_factor(inv);
  b = height / s + height_margin(inv, D);
  dbg_printf(1)("adjusted height = %.2f\n", b);
  min_prime_bits = k * b;

  prime_pool = select_classpoly_prime_pool(min_prime_bits, delta, G);

  /* FIXME: Apply torsion constraints */
  /* FIXME: Rank elts of res according to cost/benefit ratio */
  gen_sort_inplace(prime_pool, NULL, primecmp, NULL);
  prime_bits = 0.0; biggest_p = *biggest_v = *vfactors = 0;
  for (i = 1; i < lg(prime_pool); i++)
  {
    GEN q = gmael(prime_pool, i, 1);
    ulong p = q[1], v = q[3];
    *vfactors |= q[5];
    prime_bits += log2(p);
    if (p > biggest_p) biggest_p = p;
    if (v > *biggest_v) *biggest_v = v;
    if (prime_bits > b) break;
  }
  dbg_printf(1)("Selected %ld primes; largest is %lu ~ 2^%.2f\n",
             i, biggest_p, log2(biggest_p));
  return gerepilecopy(av, vecslice(prime_pool, 1, i));
}

/* This is Sutherland 2009 Algorithm 1.2. */
static long
oneroot_of_classpoly(
  ulong *j_endo, int *endo_cert, ulong j, norm_eqn_t ne, GEN jdb)
{
  pari_sp av = avma;
  long nfactors, L_bound, i;
  ulong p = ne->p, pi = ne->pi;
  GEN factors, vdepths;

  if (j == 0 || j == 1728 % p) pari_err_BUG("oneroot_of_classpoly");

  *endo_cert = 1;
  factors = gel(ne->faw, 1); nfactors = lg(factors) - 1;
  if (!nfactors) { *j_endo = j; return 1; }
  vdepths = gel(ne->faw, 2);

  /* FIXME: This should be bigger */
  L_bound = maxdd(log((double) -ne->D), (double)ne->v);

  /* Iterate over the primes L dividing w */
  for (i = 1; i <= nfactors; ++i) {
    pari_sp av2 = avma;
    GEN phi;
    long jlvl, lvl_diff, depth = vdepths[i], L = factors[i];
    if (L > L_bound) { *endo_cert = 0; break; }

    phi = polmodular_db_getp(jdb, L, p);

    /* TODO: See if I can reuse paths created in j_level_in_volcano()
     * later in {ascend,descend}_volcano(), perhaps by combining the
     * functions into one "adjust_level" function. */
    jlvl = j_level_in_volcano(phi, j, p, pi, L, depth);
    lvl_diff = z_lval(ne->u, L) - jlvl;
    if (lvl_diff < 0)
      /* j's level is less than v(u) so we must ascend */
      j = ascend_volcano(phi, j, p, pi, jlvl, L, depth, -lvl_diff);
    else if (lvl_diff > 0)
      /* otherwise j's level is greater than v(u) so we descend */
      j = descend_volcano(phi, j, p, pi, jlvl, L, depth, lvl_diff);
    set_avma(av2);
  }
  /* Prob(j has the wrong endomorphism ring) ~ \sum_{p|u_compl} 1/p
   * (and u_compl > L_bound), so just return it and rely on detection code in
   * enum_j_with_endo_ring().  Detection is that we hit a previously found
   * j-invariant earlier than expected.  OR, we evaluate class polynomials of
   * the suborders at j and if any are zero then j must be chosen again. */
  set_avma(av); *j_endo = j; return j != 0 && j != 1728 % p;
}

INLINE long
vecsmall_isin_skip(GEN v, long x, long k)
{
  long i, l = lg(v);
  for (i = k; i < l; ++i)
    if (v[i] == x) return i;
  return 0;
}

void
norm_eqn_set(norm_eqn_t ne, long D, long t, long u, long v, GEN faw, ulong p)
{
  ne->D = D;
  ne->u = u;
  ne->t = t;
  ne->v = v;
  ne->faw = faw;
  ne->p = p;
  ne->pi = get_Fl_red(ne->p);
  ne->s2 = Fl_2gener_pre(ne->p, ne->pi);
  /* select twisting parameter */
  do ne->T = random_Fl(p); while (krouu(ne->T, p) != -1);
}

INLINE ulong
Flv_powsum_pre(GEN v, ulong n, ulong p, ulong pi)
{
  long i, l = lg(v);
  ulong psum = 0;
  for (i = 1; i < l; ++i)
    psum = Fl_add(psum, Fl_powu_pre(uel(v,i), n, p, pi), p);
  return psum;
}

INLINE int
modinv_has_sign_ambiguity(long inv)
{
  switch (inv) {
  case INV_F:
  case INV_F3:
  case INV_W2W3E2:
  case INV_W2W7E2:
  case INV_W2W3:
  case INV_W2W5:
  case INV_W2W7:
  case INV_W3W3:
  case INV_W2W13:
  case INV_W3W7: return 1;
  }
  return 0;
}

INLINE int
modinv_units(int inv)
{ return modinv_is_double_eta(inv) || modinv_is_Weber(inv); }

INLINE int
adjust_signs(GEN js, ulong p, ulong pi, long inv, GEN T, long e)
{
  long negate = 0;
  long h = lg(js) - 1;
  if ((h & 1) && modinv_units(inv)) {
    ulong prod = Flv_prod_pre(js, p, pi);
    if (prod != p - 1) {
      if (prod != 1) pari_err_BUG("adjust_signs: constant term is not +/-1");
      negate = 1;
    }
  } else {
    ulong tp, t;
    tp = umodiu(T, p);
    t = Flv_powsum_pre(js, e, p, pi);
    if (t == 0) return 0;
    if (t != tp) {
      if (Fl_neg(t, p) != tp) pari_err_BUG("adjust_signs: incorrect trace");
      negate = 1;
    }
  }
  if (negate) Flv_neg_inplace(js, p);
  return 1;
}

static ulong
find_jinv(
  long *trace_tries, long *endo_tries, int *cert,
  norm_eqn_t ne, long inv, long rho_inv, GEN jdb)
{
  long found, ok = 1;
  ulong j, r;
  do {
    do {
      long tries;
      ulong j_t = 0;
      /* TODO: Set batch size according to expected number of tries and
       * experimental cost/benefit analysis. */
      tries = find_j_inv_with_given_trace(&j_t, ne, rho_inv, 0);
      if (j_t == 0)
        pari_err_BUG("polclass0: Couldn't find j-invariant with given trace.");
      dbg_printf(2)("  j-invariant %ld has trace +/-%ld (%ld tries, 1/rho = %ld)\n",
          j_t, ne->t, tries, rho_inv);
      *trace_tries += tries;

      found = oneroot_of_classpoly(&j, cert, j_t, ne, jdb);
      ++*endo_tries;
    } while (!found);

    if (modinv_is_double_eta(inv))
      ok = modfn_unambiguous_root(&r, inv, j, ne, jdb);
    else
      r = modfn_root(j, ne, inv);
  } while (!ok);
  return r;
}

static GEN
polclass_roots_modp(long *n_trace_curves,
  norm_eqn_t ne, long rho_inv, GEN G, GEN db)
{
  pari_sp av;
  ulong j = 0;
  long inv = pcp_get_inv(G), endo_tries = 0;
  int endo_cert;
  GEN res, jdb, fdb, vshape = factoru(ne->v);

  jdb = polmodular_db_for_inv(db, INV_J);
  fdb = polmodular_db_for_inv(db, inv);
  dbg_printf(2)("p = %ld, t = %ld, v = %ld\n", ne->p, ne->t, ne->v);
  av = avma;
  do {
    j = find_jinv(n_trace_curves,&endo_tries,&endo_cert, ne, inv, rho_inv, jdb);
    res = enum_roots(j, ne, fdb, G, vshape);
    if (!res && endo_cert) pari_err_BUG("polclass_roots_modp");
    if (res && !endo_cert && vecsmall_isin_skip(res, res[1], 2))
    {
      set_avma(av);
      res = NULL;
    }
  } while (!res);

  dbg_printf(2)("  j-invariant %ld has correct endomorphism ring "
             "(%ld tries)\n", j, endo_tries);
  dbg_printf(4)("  all such j-invariants: %Ps\n", res);
  return res;
}

INLINE int
modinv_inverted_involution(long inv)
{ return modinv_is_double_eta(inv); }

INLINE int
modinv_negated_involution(long inv)
{ /* determined by trial and error */
  return inv == INV_F || inv == INV_W3W5 || inv == INV_W3W7
    || inv == INV_W3W3 || inv == INV_W5W7;
}

/* Return true iff Phi_L(j0, j1) = 0. */
INLINE long
verify_edge(ulong j0, ulong j1, ulong p, ulong pi, long L, GEN fdb)
{
  pari_sp av = avma;
  GEN phi = polmodular_db_getp(fdb, L, p);
  GEN f = Flm_Fl_polmodular_evalx(phi, L, j1, p, pi);
  return gc_long(av, Flx_eval_pre(f, j0, p, pi) == 0);
}

INLINE long
verify_2path(
  ulong j1, ulong j2, ulong p, ulong pi, long L1, long L2, GEN fdb)
{
  pari_sp av = avma;
  GEN phi1 = polmodular_db_getp(fdb, L1, p);
  GEN phi2 = polmodular_db_getp(fdb, L2, p);
  GEN f = Flm_Fl_polmodular_evalx(phi1, L1, j1, p, pi);
  GEN g = Flm_Fl_polmodular_evalx(phi2, L2, j2, p, pi);
  GEN d = Flx_gcd(f, g, p);
  long n = degpol(d);
  if (n >= 2) n = Flx_nbroots(d, p);
  return gc_long(av, n);
}

static long
oriented_n_action(
  const long *ni, GEN G, GEN v, ulong p, ulong pi, GEN fdb)
{
  pari_sp av = avma;
  long i, j, k = pcp_get_k(G);
  long nr = k * (k - 1) / 2;
  const long *n = pcp_get_n(G), *m = pcp_get_m(G), *o = pcp_get_o(G), *r = pcp_get_r(G),
    *ps = pcp_get_orient_p(G), *qs = pcp_get_orient_q(G), *reps = pcp_get_orient_reps(G);
  long *signs = new_chunk(k);
  long *e = new_chunk(k);
  long *rels = new_chunk(nr);

  evec_copy(rels, r, nr);

  for (i = 0; i < k; ++i) {
    /* If generator doesn't require orientation, continue; power rels already
     * copied to *rels in initialisation */
    if (ps[i] <= 0) { signs[i] = 1; continue; }
    /* Get rep of orientation element and express it in terms of the
     * (partially) oriented presentation */
    for (j = 0; j < i; ++j) {
      long t = reps[i * k + j];
      e[j] = (signs[j] < 0 ? o[j] - t : t);
    }
    e[j] = reps[i * k + j];
    for (++j; j < k; ++j) e[j] = 0;
    evec_reduce(e, n, rels, k);
    j = evec_to_index(e, m, k);

    /* FIXME: These calls to verify_edge recalculate powers of v[0]
     * and v[j] over and over again, they also reduce Phi_{ps[i]} modulo p over
     * and over again.  Need to cache these things! */
    if (qs[i] > 1)
      signs[i] =
        (verify_2path(uel(v,1), uel(v,j+1), p, pi, ps[i], qs[i], fdb) ? 1 : -1);
    else
      /* Verify ps[i]-edge to orient ith generator */
      signs[i] =
        (verify_edge(uel(v,1), uel(v,j+1), p, pi, ps[i], fdb) ? 1 : -1);
    /* Update power relation */
    for (j = 0; j < i; ++j) {
      long t = evec_ri(r, i)[j];
      e[j] = (signs[i] * signs[j] < 0 ? o[j] - t : t);
    }
    while (j < k) e[j++] = 0;
    evec_reduce(e, n, rels, k);
    for (j = 0; j < i; ++j) evec_ri_mutate(rels, i)[j] = e[j];
    /* TODO: This is a sanity check, can be removed if everything is working */
    for (j = 0; j <= i; ++j) {
      long t = reps[i * k + j];
      e[j] = (signs[j] < 0 ? o[j] - t : t);
    }
    while (j < k) e[j++] = 0;
    evec_reduce(e, n, rels, k);
    j = evec_to_index(e, m, k);
    if (qs[i] > 1) {
      if (!verify_2path(uel(v,1), uel(v, j+1), p, pi, ps[i], qs[i], fdb))
        pari_err_BUG("oriented_n_action");
    } else {
      if (!verify_edge(uel(v,1), uel(v, j+1), p, pi, ps[i], fdb))
        pari_err_BUG("oriented_n_action");
    }
  }

  /* Orient representation of [N] relative to the torsor <signs, rels> */
  for (i = 0; i < k; ++i) e[i] = (signs[i] < 0 ? o[i] - ni[i] : ni[i]);
  evec_reduce(e, n, rels, k);
  return gc_long(av, evec_to_index(e,m,k));
}

/* F = double_eta_raw(inv) */
INLINE void
adjust_orientation(GEN F, long inv, GEN v, long e, ulong p, ulong pi)
{
  ulong j0 = uel(v, 1), je = uel(v, e);

  if (!modinv_j_from_2double_eta(F, inv, j0, je, p, pi)) {
    if (modinv_inverted_involution(inv)) Flv_inv_pre_inplace(v, p, pi);
    if (modinv_negated_involution(inv)) Flv_neg_inplace(v, p);
  }
}

static void
polclass_psum(
  GEN *psum, long *d, GEN roots, GEN primes, GEN pilist, ulong h, long inv)
{
  /* Number of consecutive CRT stabilisations before we assume we have
   * the correct answer. */
  enum { MIN_STAB_CNT = 3 };
  pari_sp av = avma, btop;
  GEN ps, psum_sqr, P;
  long i, e, stabcnt, nprimes = lg(primes) - 1;

  if ((h & 1) && modinv_units(inv)) { *psum = gen_1; *d = 0; return; }
  e = -1;
  ps = cgetg(nprimes+1, t_VECSMALL);
  do {
    e += 2;
    for (i = 1; i <= nprimes; ++i)
    {
      GEN roots_modp = gel(roots, i);
      ulong p = uel(primes, i), pi = uel(pilist, i);
      uel(ps, i) = Flv_powsum_pre(roots_modp, e, p, pi);
    }
    btop = avma;
    psum_sqr = Z_init_CRT(0, 1);
    P = gen_1;
    for (i = 1, stabcnt = 0; stabcnt < MIN_STAB_CNT && i <= nprimes; ++i)
    {
      ulong p = uel(primes, i), pi = uel(pilist, i);
      ulong ps2 = Fl_sqr_pre(uel(ps, i), p, pi);
      ulong stab = Z_incremental_CRT(&psum_sqr, ps2, &P, p);
      /* stabcnt = stab * (stabcnt + 1) */
      if (stab) ++stabcnt; else stabcnt = 0;
      if (gc_needed(av, 2)) gerepileall(btop, 2, &psum_sqr, &P);
    }
    if (stabcnt == 0 && nprimes >= MIN_STAB_CNT)
      pari_err_BUG("polclass_psum");
  } while (!signe(psum_sqr));

  if ( ! Z_issquareall(psum_sqr, psum)) pari_err_BUG("polclass_psum");

  dbg_printf(1)("Classpoly power sum (e = %ld) is %Ps; found with %.2f%% of the primes\n",
      e, *psum, 100 * (i - 1) / (double) nprimes);
  *psum = gerepileupto(av, *psum);
  *d = e;
}

static GEN
polclass_small_disc(long D, long inv, long vx)
{
  if (D == -3) return pol_x(vx);
  if (D == -4) {
    switch (inv) {
    case INV_J: return deg1pol(gen_1, stoi(-1728), vx);
    case INV_G2:return deg1pol(gen_1, stoi(-12), vx);
    default: /* no other invariants for which we can calculate H_{-4}(X) */
      pari_err_BUG("polclass_small_disc");
    }
  }
  return NULL;
}

static ulong
quadnegclassnou(long D, long *pD0, GEN *pP, GEN *pE)
{
  ulong d = (ulong)-D, d0 = coredisc2u_fact(factoru(d), -1, pP, pE);
  ulong h = uquadclassnoF_fact(d0, -1, *pP, *pE) * quadclassnos(-(long)d0);
  *pD0 = -(long)d0; return h;
}

GEN
polclass_worker(GEN q, GEN G, GEN db)
{
  GEN T = cgetg(3, t_VEC), z, P = gel(q,1), faw = gel(q,2);
  long n_curves_tested = 0, t = P[2], v = P[3], rho_inv = P[4];
  ulong p = (ulong)P[1];
  pari_sp av = avma;
  norm_eqn_t ne; norm_eqn_set(ne, pcp_get_D(G), t, pcp_get_u(G), v, faw, p);
  z = polclass_roots_modp(&n_curves_tested, ne, rho_inv, G, db);
  gel(T,1) = gerepileuptoleaf(av, z);
  gel(T,2) = mkvecsmall3(ne->p, ne->pi, n_curves_tested); return T;
}

GEN
polclass0(long D, long inv, long vx, GEN *db)
{
  pari_sp av = avma;
  GEN primes, H, plist, pilist, Pu, Eu;
  long n_curves_tested = 0, filter = 1, pending = 0, cnt = 0;
  long D0, nprimes, s, i, j, del, ni, orient, h, p1, p2, k;
  ulong u, vfactors, biggest_v;
  GEN G, worker, vec;
  double height;
  static const double delta = 0.5;
  struct pari_mt pt;

  if (D >= -4) return polclass_small_disc(D, inv, vx);

  h = quadnegclassnou(D, &D0, &Pu, &Eu);
  u = D == D0? 1: (ulong)sqrt(D / D0); /* u = \prod Pu[i]^Eu[i] */
  dbg_printf(1)("D = %ld, conductor = %ld, inv = %ld\n", D, u, inv);

  ni = modinv_degree(&p1, &p2, inv);
  orient = modinv_is_double_eta(inv) && kross(D, p1) && kross(D, p2);

  G = classgp_make_pcp(&height, &ni, h, D, D0, u, Pu, Eu, inv, filter, orient);
  primes = select_classpoly_primes(&vfactors, &biggest_v, delta, G, height);

  /* Prepopulate *db with all the modpolys we might need */
  if (u > 1)
  { /* L_bound in oneroot_of_classpoly() */
    long l = lg(Pu), maxL = maxdd(log((double) -D), (double)biggest_v);
    for (i = 1; i < l; i++)
    {
      long L = Pu[i];
      if (L > maxL) break;
      polmodular_db_add_level(db, L, INV_J);
    }
  }
  for (i = 0; vfactors; ++i) {
    if (vfactors & 1UL)
      polmodular_db_add_level(db, SMALL_PRIMES[i], INV_J);
    vfactors >>= 1;
  }
  if (p1 > 1) polmodular_db_add_level(db, p1, INV_J);
  if (p2 > 1) polmodular_db_add_level(db, p2, INV_J);
  s = !!pcp_get_L0(G); k = pcp_get_k(G);
  polmodular_db_add_levels(db, pcp_get_L(G) + s, k - s, inv);
  if (orient)
  {
    GEN orient_p = pcp_get_orient_p(G);
    GEN orient_q = pcp_get_orient_q(G);
    for (i = 0; i < k; ++i)
    {
      if (orient_p[i] > 1) polmodular_db_add_level(db, orient_p[i], inv);
      if (orient_q[i] > 1) polmodular_db_add_level(db, orient_q[i], inv);
    }
  }
  nprimes = lg(primes) - 1;
  worker = snm_closure(is_entry("_polclass_worker"), mkvec2(G, *db));
  H = cgetg(nprimes + 1, t_VEC);
  plist = cgetg(nprimes + 1, t_VECSMALL);
  pilist = cgetg(nprimes + 1, t_VECSMALL);
  vec = cgetg(2, t_VEC);
  dbg_printf(0)("Calculating class polynomial of disc %ld and inv %ld:", D, inv);
  mt_queue_start_lim(&pt, worker, nprimes);
  for (i = 1; i <= nprimes || pending; i++)
  {
    long workid;
    GEN done;
    if (i <= nprimes) gel(vec,1) = gel(primes,i);
    mt_queue_submit(&pt, i, i <= nprimes? vec: NULL);
    done = mt_queue_get(&pt, &workid, &pending);
    if (done)
    {
      gel(H, workid) = gel(done,1);
      uel(plist, workid) = mael(done,2,1);
      uel(pilist, workid) = mael(done,2,2);
      n_curves_tested += mael(done,2,3);
      cnt++;
      if (DEBUGLEVEL>2 && (cnt & 3L)==0)
        err_printf(" %ld%%", cnt*100/nprimes);
    }
  }
  mt_queue_end(&pt);
  dbg_printf(0)(" done\n");

  if (orient) {
    GEN nvec = new_chunk(k);
    GEN fdb = polmodular_db_for_inv(*db, inv);
    GEN F = double_eta_raw(inv);
    index_to_evec((long *)nvec, ni, pcp_get_m(G), k);
    for (i = 1; i <= nprimes; ++i) {
      GEN v = gel(H, i);
      ulong p = uel(plist, i), pi = uel(pilist, i);
      long oni = oriented_n_action(nvec, G, v, p, pi, fdb);
      adjust_orientation(F, inv, v, oni + 1, p, pi);
    }
  }

  if (modinv_has_sign_ambiguity(inv)) {
    GEN psum;
    long e;
    polclass_psum(&psum, &e, H, plist, pilist, h, inv);
    for (i = 1; i <= nprimes; ++i) {
      GEN v = gel(H, i);
      ulong p = uel(plist, i), pi = uel(pilist, i);
      if (!adjust_signs(v, p, pi, inv, psum, e))
        uel(plist, i) = 0;
    }
  }

  for (i = 1, j = 1, del = 0; i <= nprimes; ++i)
  {
    pari_sp av2 = avma;
    ulong p = uel(plist, i);
    long l;
    GEN v;
    if (!p) { del++; continue; }
    v = Flv_roots_to_pol(gel(H,i), p, vx); l = lg(v);
    *++v = evaltyp(t_VECSMALL) | evallg(l-1); /* Flx_to_Flv inplace */
    uel(plist, j) = p;
    gel(H, j++) = gerepileuptoleaf(av2, v);
  }
  setlg(H,nprimes+1-del);
  setlg(plist,nprimes+1-del);

  dbg_printf(1)("Total number of curves tested: %ld\n", n_curves_tested);
  H = ncV_chinese_center(H, plist, NULL);
  dbg_printf(1)("Result height: %.2f\n", dbllog2(gsupnorm(H, DEFAULTPREC)));
  return gerepilecopy(av, RgV_to_RgX(H, vx));
}

void
check_modinv(long inv)
{
  switch (inv) {
  case INV_J:
  case INV_F:
  case INV_F2:
  case INV_F3:
  case INV_F4:
  case INV_G2:
  case INV_W2W3:
  case INV_F8:
  case INV_W3W3:
  case INV_W2W5:
  case INV_W2W7:
  case INV_W3W5:
  case INV_W3W7:
  case INV_W2W3E2:
  case INV_W2W5E2:
  case INV_W2W13:
  case INV_W2W7E2:
  case INV_W3W3E2:
  case INV_W5W7:
  case INV_W3W13:
    break;
  default:
    pari_err_DOMAIN("polmodular", "inv", "invalid invariant", stoi(inv), gen_0);
  }
}

GEN
polclass(GEN DD, long inv, long vx)
{
  GEN db, H;
  long D;

  if (vx < 0) vx = 0;
  check_quaddisc_imag(DD, NULL, "polclass");
  check_modinv(inv);

  D = itos(DD);
  if (!modinv_good_disc(inv, D))
    pari_err_DOMAIN("polclass", "D", "incompatible with given invariant", stoi(inv), DD);

  db = polmodular_db_init(inv);
  H = polclass0(D, inv, vx, &db);
  gunclone_deep(db); return H;
}
