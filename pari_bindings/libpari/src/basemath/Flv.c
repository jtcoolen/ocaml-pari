/* Copyright (C) 2000-2019  The PARI group.

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
Flv_to_ZV(GEN x)
{ pari_APPLY_type(t_VEC, utoi(x[i])) }

GEN
Flc_to_ZC(GEN x)
{ pari_APPLY_type(t_COL, utoi(x[i])) }

GEN
Flm_to_ZM(GEN x)
{ pari_APPLY_type(t_MAT, Flc_to_ZC(gel(x,i))) }

GEN
Flc_to_ZC_inplace(GEN z)
{
  long i, l = lg(z);
  for (i=1; i<l; i++) gel(z,i) = utoi(z[i]);
  settyp(z, t_COL);
  return z;
}

GEN
Flm_to_ZM_inplace(GEN z)
{
  long i, l = lg(z);
  for (i=1; i<l; i++) Flc_to_ZC_inplace(gel(z, i));
  return z;
}

static GEN
Flm_solve_upper_1(GEN U, GEN B, ulong p, ulong pi)
{ return Flm_Fl_mul_pre(B, Fl_inv(ucoeff(U, 1, 1), p), p, pi); }

static GEN
Flm_rsolve_upper_2(GEN U, GEN B, ulong p, ulong pi)
{
  ulong a = ucoeff(U, 1, 1), b = ucoeff(U, 1, 2), d = ucoeff(U, 2, 2);
  ulong D = Fl_mul_pre(a, d, p, pi), Dinv = Fl_inv(D, p);
  ulong ainv = Fl_mul_pre(d, Dinv, p, pi);
  ulong dinv = Fl_mul_pre(a, Dinv, p, pi);
  GEN B1 = rowslice(B, 1, 1);
  GEN B2 = rowslice(B, 2, 2);
  GEN X2 = Flm_Fl_mul_pre(B2, dinv, p, pi);
  GEN X1 = Flm_Fl_mul_pre(Flm_sub(B1, Flm_Fl_mul_pre(X2, b, p, pi), p),
                      ainv, p, pi);
  return vconcat(X1, X2);
}

/* solve U*X = B,  U upper triangular and invertible */
static GEN
Flm_rsolve_upper_pre(GEN U, GEN B, ulong p, ulong pi)
{
  long n = lg(U) - 1, n1;
  GEN U2, U11, U12, U22, B1, B2, X1, X2, X;
  pari_sp av = avma;

  if (n == 0) return B;
  if (n == 1) return Flm_solve_upper_1(U, B, p, pi);
  if (n == 2) return Flm_rsolve_upper_2(U, B, p, pi);
  n1 = (n + 1)/2;
  U2 = vecslice(U, n1 + 1, n);
  U22 = rowslice(U2, n1 + 1, n);
  B2 = rowslice(B, n1 + 1, n);
  X2 = Flm_rsolve_upper_pre(U22, B2, p, pi);
  U12 = rowslice(U2, 1, n1);
  B1 = rowslice(B, 1, n1);
  B1 = Flm_sub(B1, Flm_mul_pre(U12, X2, p, pi), p);
  if (gc_needed(av, 1)) gerepileall(av, 2, &B1, &X2);
  U11 = matslice(U, 1,n1, 1,n1);
  X1 = Flm_rsolve_upper_pre(U11, B1, p, pi);
  X = vconcat(X1, X2);
  if (gc_needed(av, 1)) X = gerepilecopy(av, X);
  return X;
}

static GEN
Flm_lsolve_upper_2(GEN U, GEN B, ulong p, ulong pi)
{
  ulong a = ucoeff(U, 1, 1), b = ucoeff(U, 1, 2), d = ucoeff(U, 2, 2);
  ulong D = Fl_mul_pre(a, d, p, pi), Dinv = Fl_inv(D, p);
  ulong ainv = Fl_mul_pre(d, Dinv, p, pi);
  ulong dinv = Fl_mul_pre(a, Dinv, p, pi);
  GEN B1 = vecslice(B, 1, 1);
  GEN B2 = vecslice(B, 2, 2);
  GEN X1 = Flm_Fl_mul_pre(B1, ainv, p, pi);
  GEN X2 = Flm_Fl_mul_pre(Flm_sub(B2, Flm_Fl_mul_pre(X1, b, p, pi), p),
                      dinv, p, pi);
  return shallowconcat(X1, X2);
}

/* solve X*U = B,  U upper triangular and invertible */
static GEN
Flm_lsolve_upper_pre(GEN U, GEN B, ulong p, ulong pi)
{
  long n = lg(U) - 1, n1;
  GEN U2, U11, U12, U22, B1, B2, X1, X2, X;
  pari_sp av = avma;

  if (n == 0) return B;
  if (n == 1) return Flm_solve_upper_1(U, B, p, pi);
  if (n == 2) return Flm_lsolve_upper_2(U, B, p, pi);
  n1 = (n + 1)/2;
  U2 = vecslice(U, n1 + 1, n);
  U11 = matslice(U, 1,n1, 1,n1);
  U12 = rowslice(U2, 1, n1);
  U22 = rowslice(U2, n1 + 1, n);
  B1 = vecslice(B, 1, n1);
  B2 = vecslice(B, n1 + 1, n);
  X1 = Flm_lsolve_upper_pre(U11, B1, p, pi);
  B2 = Flm_sub(B2, Flm_mul_pre(X1, U12, p, pi), p);
  if (gc_needed(av, 1)) gerepileall(av, 3, &B2, &U22, &X1);
  X2 = Flm_lsolve_upper_pre(U22, B2, p, pi);
  X = shallowconcat(X1, X2);
  if (gc_needed(av, 1)) X = gerepilecopy(av, X);
  return X;
}

static GEN
Flm_rsolve_lower_unit_2(GEN L, GEN A, ulong p, ulong pi)
{
  GEN X1 = rowslice(A, 1, 1);
  GEN X2 = Flm_sub(rowslice(A, 2, 2),
                   Flm_Fl_mul_pre(X1, ucoeff(L, 2, 1), p, pi), p);
  return vconcat(X1, X2);
}

/* solve L*X = A,  L lower triangular with ones on the diagonal
* (at least as many rows as columns) */
static GEN
Flm_rsolve_lower_unit_pre(GEN L, GEN A, ulong p, ulong pi)
{
  long m = lg(L) - 1, m1, n;
  GEN L1, L11, L21, L22, A1, A2, X1, X2, X;
  pari_sp av = avma;

  if (m == 0) return zero_Flm(0, lg(A) - 1);
  if (m == 1) return rowslice(A, 1, 1);
  if (m == 2) return Flm_rsolve_lower_unit_2(L, A, p, pi);
  m1 = (m + 1)/2;
  n = nbrows(L);
  L1 = vecslice(L, 1, m1);
  L11 = rowslice(L1, 1, m1);
  L21 = rowslice(L1, m1 + 1, n);
  A1 = rowslice(A, 1, m1);
  X1 = Flm_rsolve_lower_unit_pre(L11, A1, p, pi);
  A2 = rowslice(A, m1 + 1, n);
  A2 = Flm_sub(A2, Flm_mul_pre(L21, X1, p, pi), p);
  if (gc_needed(av, 1)) gerepileall(av, 2, &A2, &X1);
  L22 = matslice(L, m1+1,n, m1+1,m);
  X2 = Flm_rsolve_lower_unit_pre(L22, A2, p, pi);
  X = vconcat(X1, X2);
  if (gc_needed(av, 1)) X = gerepilecopy(av, X);
  return X;
}

static GEN
Flm_lsolve_lower_unit_2(GEN L, GEN A, ulong p, ulong pi)
{
  GEN X2 = vecslice(A, 2, 2);
  GEN X1 = Flm_sub(vecslice(A, 1, 1),
                   Flm_Fl_mul_pre(X2, ucoeff(L, 2, 1), p, pi), p);
  return shallowconcat(X1, X2);
}

/* solve L*X = A,  L square lower triangular with ones on the diagonal */
static GEN
Flm_lsolve_lower_unit_pre(GEN L, GEN A, ulong p, ulong pi)
{
  long m = lg(L) - 1, m1;
  GEN L1, L2, L11, L21, L22, A1, A2, X1, X2, X;
  pari_sp av = avma;

  if (m <= 1) return A;
  if (m == 2) return Flm_lsolve_lower_unit_2(L, A, p, pi);
  m1 = (m + 1)/2;
  L2 = vecslice(L, m1 + 1, m);
  L22 = rowslice(L2, m1 + 1, m);
  A2 = vecslice(A, m1 + 1, m);
  X2 = Flm_lsolve_lower_unit_pre(L22, A2, p, pi);
  if (gc_needed(av, 1)) X2 = gerepilecopy(av, X2);
  L1 = vecslice(L, 1, m1);
  L21 = rowslice(L1, m1 + 1, m);
  A1 = vecslice(A, 1, m1);
  A1 = Flm_sub(A1, Flm_mul_pre(X2, L21, p, pi), p);
  L11 = rowslice(L1, 1, m1);
  if (gc_needed(av, 1)) gerepileall(av, 3, &A1, &L11, &X2);
  X1 = Flm_lsolve_lower_unit_pre(L11, A1, p, pi);
  X = shallowconcat(X1, X2);
  if (gc_needed(av, 1)) X = gerepilecopy(av, X);
  return X;
}

/* destroy A */
static long
Flm_CUP_basecase(GEN A, GEN *R, GEN *C, GEN *U, GEN *P, ulong p, ulong pi)
{
  long i, j, k, m = nbrows(A), n = lg(A) - 1, pr, pc;
  ulong u, v;

  if (P) *P = identity_perm(n);
  *R = cgetg(m + 1, t_VECSMALL);
  for (j = 1, pr = 0; j <= n; j++)
  {
    for (pr++, pc = 0; pr <= m; pr++)
    {
      for (k = j; k <= n; k++) { v = ucoeff(A, pr, k); if (!pc && v) pc = k; }
      if (pc) break;
    }
    if (!pc) break;
    (*R)[j] = pr;
    if (pc != j)
    {
      swap(gel(A, j), gel(A, pc));
      if (P) lswap((*P)[j], (*P)[pc]);
    }
    u = Fl_inv(ucoeff(A, pr, j), p);
    for (i = pr + 1; i <= m; i++)
    {
      v = Fl_mul_pre(ucoeff(A, i, j), u, p, pi);
      ucoeff(A, i, j) = v;
      v = Fl_neg(v, p);
      for (k = j + 1; k <= n; k++)
        ucoeff(A, i, k) = Fl_addmul_pre(ucoeff(A, i, k),
                                        ucoeff(A, pr, k), v, p, pi);
    }
  }
  setlg(*R, j);
  *C = vecslice(A, 1, j - 1);
  if (U) *U = rowpermute(A, *R);
  return j - 1;
}

static const long Flm_CUP_LIMIT = 8;

static long
Flm_CUP_pre(GEN A, GEN *R, GEN *C, GEN *U, GEN *P, ulong p, ulong pi)
{
  long m = nbrows(A), m1, n = lg(A) - 1, i, r1, r2, r;
  GEN R1, C1, U1, P1, R2, C2, U2, P2;
  GEN A1, A2, B2, C21, U11, U12, T21, T22;
  pari_sp av = avma;

  if (m < Flm_CUP_LIMIT || n < Flm_CUP_LIMIT)
    /* destroy A; not called at the outermost recursion level */
    return Flm_CUP_basecase(A, R, C, U, P, p, pi);
  m1 = (minss(m, n) + 1)/2;
  A1 = rowslice(A, 1, m1);
  A2 = rowslice(A, m1 + 1, m);
  r1 = Flm_CUP_pre(A1, &R1, &C1, &U1, &P1, p, pi);
  if (r1 == 0)
  {
    r2 = Flm_CUP_pre(A2, &R2, &C2, &U2, &P2, p, pi);
    *R = cgetg(r2 + 1, t_VECSMALL);
    for (i = 1; i <= r2; i++) (*R)[i] = R2[i] + m1;
    *C = vconcat(zero_Flm(m1, r2), C2);
    *U = U2;
    *P = P2;
    r = r2;
  }
  else
  {
    U11 = vecslice(U1, 1, r1);
    U12 = vecslice(U1, r1 + 1, n);
    T21 = vecslicepermute(A2, P1, 1, r1);
    T22 = vecslicepermute(A2, P1, r1 + 1, n);
    C21 = Flm_lsolve_upper_pre(U11, T21, p, pi);
    if (gc_needed(av, 1))
      gerepileall(av, 7, &R1, &C1, &P1, &U11, &U12, &T22, &C21);
    B2 = Flm_sub(T22, Flm_mul_pre(C21, U12, p, pi), p);
    r2 = Flm_CUP_pre(B2, &R2, &C2, &U2, &P2, p, pi);
    r = r1 + r2;
    *R = cgetg(r + 1, t_VECSMALL);
    for (i = 1; i <= r1; i++) (*R)[i] = R1[i];
    for (;      i <= r; i++)  (*R)[i] = R2[i - r1] + m1;
    *C = shallowconcat(vconcat(C1, C21),
                       vconcat(zero_Flm(m1, r2), C2));
    *U = shallowconcat(vconcat(U11, zero_Flm(r2, r1)),
                       vconcat(vecpermute(U12, P2), U2));
    *P = cgetg(n + 1, t_VECSMALL);
    for (i = 1; i <= r1; i++) (*P)[i] = P1[i];
    for (; i <= n; i++)       (*P)[i] = P1[P2[i - r1] + r1];
  }
  if (gc_needed(av, 1)) gerepileall(av, 4, R, C, U, P);
  return r;
}

static long
Flm_echelon_gauss(GEN A, GEN *R, GEN *C, ulong p, ulong pi)
{ return Flm_CUP_basecase(A, R, C, NULL, NULL, p, pi); }

/* complement of a strictly increasing subsequence of (1, 2, ..., n) */
static GEN
indexcompl(GEN v, long n)
{
  long i, j, k, m = lg(v) - 1;
  GEN w = cgetg(n - m + 1, t_VECSMALL);
  for (i = j = k = 1; i <= n; i++)
    if (j <= m && v[j] == i) j++; else w[k++] = i;
  return w;
}

/* column echelon form */
static long
Flm_echelon_pre(GEN A, GEN *R, GEN *C, ulong p, ulong pi)
{
  long j, j1, j2, m = nbrows(A), n = lg(A) - 1, n1, r, r1, r2;
  GEN A1, A2, R1, R1c, C1, R2, C2;
  GEN A12, A22, B2, C11, C21, M12;
  pari_sp av = avma;

  if (m < Flm_CUP_LIMIT || n < Flm_CUP_LIMIT)
    return Flm_echelon_gauss(Flm_copy(A), R, C, p, pi);

  n1 = (n + 1)/2;
  A1 = vecslice(A, 1, n1);
  A2 = vecslice(A, n1 + 1, n);
  r1 = Flm_echelon_pre(A1, &R1, &C1, p, pi);
  if (!r1) return Flm_echelon_pre(A2, R, C, p, pi);
  if (r1 == m) { *R = R1; *C = C1; return r1; }

  R1c = indexcompl(R1, m);
  C11 = rowpermute(C1, R1);
  C21 = rowpermute(C1, R1c);
  A12 = rowpermute(A2, R1);
  A22 = rowpermute(A2, R1c);
  M12 = Flm_rsolve_lower_unit_pre(C11, A12, p, pi);
  B2 = Flm_sub(A22, Flm_mul_pre(C21, M12, p, pi), p);
  r2 = Flm_echelon_pre(B2, &R2, &C2, p, pi);
  if (!r2) { *R = R1; *C = C1; r = r1; }
  else
  {
    R2 = perm_mul(R1c, R2);
    C2 = rowpermute(vconcat(zero_Flm(r1, r2), C2),
                    perm_inv(vecsmall_concat(R1, R1c)));
    r = r1 + r2;
    *R = cgetg(r + 1, t_VECSMALL);
    *C = cgetg(r + 1, t_MAT);
    for (j = j1 = j2 = 1; j <= r; j++)
      if (j2 > r2 || (j1 <= r1 && R1[j1] < R2[j2]))
      {
        gel(*C, j) = gel(C1, j1);
        (*R)[j] = R1[j1++];
      }
      else
      {
        gel(*C, j) = gel(C2, j2);
        (*R)[j] = R2[j2++];
      }
  }
  if (gc_needed(av, 1)) gerepileall(av, 2, R, C);
  return r;
}

static void /* assume m < p */
_Fl_addmul(GEN b, long k, long i, ulong m, ulong p, ulong pi)
{
  uel(b,k) = Fl_addmul_pre(uel(b, k), m, uel(b, i), p, pi);
}
static void /* same m = 1 */
_Fl_add(GEN b, long k, long i, ulong p)
{
  uel(b,k) = Fl_add(uel(b,k), uel(b,i), p);
}
static void /* assume m < p && SMALL_ULONG(p) && (! (b[i] & b[k] & HIGHMASK)) */
_Fl_addmul_OK(GEN b, long k, long i, ulong m, ulong p)
{
  uel(b,k) += m * uel(b,i);
  if (uel(b,k) & HIGHMASK) uel(b,k) %= p;
}
static void /* assume SMALL_ULONG(p) && (! (b[i] & b[k] & HIGHMASK)) */
_Fl_add_OK(GEN b, long k, long i, ulong p)
{
  uel(b,k) += uel(b,i);
  if (uel(b,k) & HIGHMASK) uel(b,k) %= p;
}

/* assume 0 <= a[i,j] < p */
static GEN
Fl_get_col_OK(GEN a, GEN b, long li, ulong p)
{
  GEN u = cgetg(li+1,t_VECSMALL);
  ulong m = uel(b,li) % p;
  long i,j;

  uel(u,li) = (m * ucoeff(a,li,li)) % p;
  for (i = li-1; i > 0; i--)
  {
    m = p - uel(b,i)%p;
    for (j = i+1; j <= li; j++) {
      if (m & HIGHBIT) m %= p;
      m += ucoeff(a,i,j) * uel(u,j); /* 0 <= u[j] < p */
    }
    m %= p;
    if (m) m = ((p-m) * ucoeff(a,i,i)) % p;
    uel(u,i) = m;
  }
  return u;
}
static GEN
Fl_get_col(GEN a, GEN b, long li, ulong p)
{
  GEN u = cgetg(li+1,t_VECSMALL);
  ulong m = uel(b,li) % p;
  long i,j;

  uel(u,li) = Fl_mul(m, ucoeff(a,li,li), p);
  for (i=li-1; i>0; i--)
  {
    m = b[i]%p;
    for (j = i+1; j <= li; j++)
      m = Fl_sub(m, Fl_mul(ucoeff(a,i,j), uel(u,j), p), p);
    if (m) m = Fl_mul(m, ucoeff(a,i,i), p);
    uel(u,i) = m;
  }
  return u;
}

static GEN
Flm_ker_gauss_OK(GEN x, ulong p, long deplin)
{
  GEN y, c, d;
  long i, j, k, r, t, m, n;
  ulong a;

  n = lg(x)-1;
  m=nbrows(x); r=0;

  c = zero_zv(m);
  d = cgetg(n+1, t_VECSMALL);
  a = 0; /* for gcc -Wall */
  for (k=1; k<=n; k++)
  {
    for (j=1; j<=m; j++)
      if (!c[j])
      {
        a = ucoeff(x,j,k) % p;
        if (a) break;
      }
    if (j > m)
    {
      if (deplin==1) {
        c = cgetg(n+1, t_VECSMALL);
        for (i=1; i<k; i++) c[i] = ucoeff(x,d[i],k) % p;
        c[k] = 1; for (i=k+1; i<=n; i++) c[i] = 0;
        return c;
      }
      r++; d[k] = 0;
    }
    else
    {
      ulong piv = p - Fl_inv(a, p); /* -1/a */
      c[j] = k; d[k] = j;
      ucoeff(x,j,k) = p-1;
      if (piv != 1)
        for (i=k+1; i<=n; i++) ucoeff(x,j,i) = (piv * ucoeff(x,j,i)) % p;
      for (t=1; t<=m; t++)
      {
        if (t == j) continue;

        piv = ( ucoeff(x,t,k) %= p );
        if (!piv) continue;
        if (piv == 1)
          for (i=k+1; i<=n; i++) _Fl_add_OK(gel(x,i),t,j, p);
        else
          for (i=k+1; i<=n; i++) _Fl_addmul_OK(gel(x,i),t,j,piv, p);
      }
    }
  }
  if (deplin==1) return NULL;

  y = cgetg(r+1, t_MAT);
  for (j=k=1; j<=r; j++,k++)
  {
    GEN C = cgetg(n+1, t_VECSMALL);

    gel(y,j) = C; while (d[k]) k++;
    for (i=1; i<k; i++)
      if (d[i])
        uel(C,i) = ucoeff(x,d[i],k) % p;
      else
        uel(C,i) = 0UL;
    uel(C,k) = 1UL; for (i=k+1; i<=n; i++) uel(C,i) = 0UL;
  }
  if (deplin == 2) {
    GEN pc = cgetg(n - r + 1, t_VECSMALL);  /* indices of pivot columns */
    for (i = j = 1; j <= n; j++)
      if (d[j]) pc[i++] = j;
    return mkvec2(y, pc);
  }
  return y;
}

/* in place, destroy x */
static GEN
Flm_ker_gauss(GEN x, ulong p, long deplin)
{
  GEN y, c, d;
  long i, j, k, r, t, m, n;
  ulong a, pi;
  n = lg(x)-1;
  if (!n) return cgetg(1,t_MAT);
  if (SMALL_ULONG(p)) return Flm_ker_gauss_OK(x, p, deplin);
  pi = get_Fl_red(p);

  m=nbrows(x); r=0;

  c = zero_zv(m);
  d = cgetg(n+1, t_VECSMALL);
  a = 0; /* for gcc -Wall */
  for (k=1; k<=n; k++)
  {
    for (j=1; j<=m; j++)
      if (!c[j])
      {
        a = ucoeff(x,j,k);
        if (a) break;
      }
    if (j > m)
    {
      if (deplin==1) {
        c = cgetg(n+1, t_VECSMALL);
        for (i=1; i<k; i++) c[i] = ucoeff(x,d[i],k);
        c[k] = 1; for (i=k+1; i<=n; i++) c[i] = 0;
        return c;
      }
      r++; d[k] = 0;
    }
    else
    {
      ulong piv = p - Fl_inv(a, p); /* -1/a */
      c[j] = k; d[k] = j;
      ucoeff(x,j,k) = p-1;
      if (piv != 1)
        for (i=k+1; i<=n; i++)
          ucoeff(x,j,i) = Fl_mul_pre(piv, ucoeff(x,j,i), p, pi);
      for (t=1; t<=m; t++)
      {
        if (t == j) continue;

        piv = ucoeff(x,t,k);
        if (!piv) continue;
        if (piv == 1)
          for (i=k+1; i<=n; i++) _Fl_add(gel(x,i),t,j,p);
        else
          for (i=k+1; i<=n; i++) _Fl_addmul(gel(x,i),t,j,piv,p, pi);
      }
    }
  }
  if (deplin==1) return NULL;

  y = cgetg(r+1, t_MAT);
  for (j=k=1; j<=r; j++,k++)
  {
    GEN C = cgetg(n+1, t_VECSMALL);

    gel(y,j) = C; while (d[k]) k++;
    for (i=1; i<k; i++)
      if (d[i])
        uel(C,i) = ucoeff(x,d[i],k);
      else
        uel(C,i) = 0UL;
    uel(C,k) = 1UL; for (i=k+1; i<=n; i++) uel(C,i) = 0UL;
  }
  if (deplin == 2) {
    GEN pc = cgetg(n - r + 1, t_VECSMALL);  /* indices of pivot columns */
    for (i = j = 1; j <= n; j++)
      if (d[j]) pc[i++] = j;
    return mkvec2(y, pc);
  }
  return y;
}

GEN
Flm_intersect_i(GEN x, GEN y, ulong p)
{
  long j, lx = lg(x);
  GEN z;

  if (lx==1 || lg(y)==1) return cgetg(1,t_MAT);
  z = Flm_ker(shallowconcat(x,y), p);
  for (j=lg(z)-1; j; j--) setlg(gel(z,j),lx);
  return Flm_mul(x,z,p);
}
GEN
Flm_intersect(GEN x, GEN y, ulong p)
{
  pari_sp av = avma;
  return gerepileupto(av, Flm_image(Flm_intersect_i(x, y, p), p));
}

static GEN
Flm_ker_echelon(GEN x, ulong p, long pivots) {
  pari_sp av = avma;
  GEN R, Rc, C, C1, C2, S, K;
  long n = lg(x) - 1, r;
  ulong pi = get_Fl_red(p);
  r = Flm_echelon_pre(Flm_transpose(x), &R, &C, p, pi);
  Rc = indexcompl(R, n);
  C1 = rowpermute(C, R);
  C2 = rowpermute(C, Rc);
  S = Flm_lsolve_lower_unit_pre(C1, C2, p, pi);
  K = vecpermute(shallowconcat(Flm_neg(S, p), matid_Flm(n - r)),
                 perm_inv(vecsmall_concat(R, Rc)));
  K = Flm_transpose(K);
  if (pivots)
    return gerepilecopy(av, mkvec2(K, R));
  return gerepileupto(av, K);
}

static GEN
Flm_deplin_echelon(GEN x, ulong p) {
  pari_sp av = avma;
  GEN R, Rc, C, C1, C2, s, v;
  long i, n = lg(x) - 1, r;
  ulong pi = get_Fl_red(p);
  r = Flm_echelon_pre(Flm_transpose(x), &R, &C, p, pi);
  if (r == n) return gc_NULL(av);
  Rc = indexcompl(R, n);
  i = Rc[1];
  C1 = rowpermute(C, R);
  C2 = rowslice(C, i, i);
  s = Flm_row(Flm_lsolve_lower_unit_pre(C1, C2, p, pi), 1);
  v = vecsmallpermute(vecsmall_concat(Flv_neg(s, p), vecsmall_ei(n - r, 1)),
                 perm_inv(vecsmall_concat(R, Rc)));
  return gerepileuptoleaf(av, v);
}

static GEN
Flm_ker_i(GEN x, ulong p, long deplin, long inplace) {
  if (lg(x) - 1 >= Flm_CUP_LIMIT && nbrows(x) >= Flm_CUP_LIMIT)
    switch(deplin) {
    case 0: return Flm_ker_echelon(x, p, 0);
    case 1: return Flm_deplin_echelon(x, p);
    case 2: return Flm_ker_echelon(x, p, 1);
    }
  return Flm_ker_gauss(inplace? x: Flm_copy(x), p, deplin);
}

GEN
Flm_ker_sp(GEN x, ulong p, long deplin) { return Flm_ker_i(x, p, deplin, 1); }
GEN
Flm_ker(GEN x, ulong p) { return Flm_ker_i(x, p, 0, 0); }
GEN
Flm_deplin(GEN x, ulong p) { return Flm_ker_i(x, p, 1, 0); }

/* in place, destroy a, SMALL_ULONG(p) is TRUE */
static ulong
Flm_det_gauss_OK(GEN a, long nbco, ulong p)
{
  long i,j,k, s = 1;
  ulong q, x = 1;

  for (i=1; i<nbco; i++)
  {
    for(k=i; k<=nbco; k++)
    {
      ulong c = ucoeff(a,k,i) % p;
      ucoeff(a,k,i) = c;
      if (c) break;
    }
    for(j=k+1; j<=nbco; j++) ucoeff(a,j,i) %= p;
    if (k > nbco) return ucoeff(a,i,i);
    if (k != i)
    { /* exchange the lines s.t. k = i */
      for (j=i; j<=nbco; j++) lswap(ucoeff(a,i,j), ucoeff(a,k,j));
      s = -s;
    }
    q = ucoeff(a,i,i);

    if (x & HIGHMASK) x %= p;
    x *= q;
    q = Fl_inv(q,p);
    for (k=i+1; k<=nbco; k++)
    {
      ulong m = ucoeff(a,i,k) % p;
      if (!m) continue;

      m = p - ((m*q)%p);
      for (j=i+1; j<=nbco; j++)
      {
        ulong c = ucoeff(a,j,k);
        if (c & HIGHMASK) c %= p;
        ucoeff(a,j,k) = c  + m*ucoeff(a,j,i);
      }
    }
  }
  if (x & HIGHMASK) x %= p;
  q = ucoeff(a,nbco,nbco);
  if (q & HIGHMASK) q %= p;
  x = (x*q) % p;
  if (s < 0 && x) x = p - x;
  return x;
}

/* in place, destroy a */
static ulong
Flm_det_gauss(GEN a, ulong p)
{
  long i,j,k, s = 1, nbco = lg(a)-1;
  ulong pi, q, x = 1;

  if (SMALL_ULONG(p)) return Flm_det_gauss_OK(a, nbco, p);
  pi = get_Fl_red(p);
  for (i=1; i<nbco; i++)
  {
    for(k=i; k<=nbco; k++)
      if (ucoeff(a,k,i)) break;
    if (k > nbco) return ucoeff(a,i,i);
    if (k != i)
    { /* exchange the lines s.t. k = i */
      for (j=i; j<=nbco; j++) lswap(ucoeff(a,i,j), ucoeff(a,k,j));
      s = -s;
    }
    q = ucoeff(a,i,i);

    x = Fl_mul_pre(x, q, p, pi);
    q = Fl_inv(q,p);
    for (k=i+1; k<=nbco; k++)
    {
      ulong m = ucoeff(a,i,k);
      if (!m) continue;

      m = Fl_mul_pre(m, q, p, pi);
      for (j=i+1; j<=nbco; j++)
        ucoeff(a,j,k) = Fl_sub(ucoeff(a,j,k), Fl_mul_pre(m,ucoeff(a,j,i), p, pi), p);
    }
  }
  if (s < 0) x = Fl_neg(x, p);
  return Fl_mul(x, ucoeff(a,nbco,nbco), p);
}

static ulong
Flm_det_CUP(GEN a, ulong p) {
  GEN R, C, U, P;
  long i, n = lg(a) - 1, r;
  ulong d;
  ulong pi = get_Fl_red(p);
  r = Flm_CUP_pre(a, &R, &C, &U, &P, p, pi);
  if (r < n)
    d = 0;
  else {
    d = perm_sign(P) == 1? 1: p-1;
    for (i = 1; i <= n; i++)
      d = Fl_mul_pre(d, ucoeff(U, i, i), p, pi);
  }
  return d;
}

static ulong
Flm_det_i(GEN x, ulong p, long inplace) {
  pari_sp av = avma;
  ulong d;
  if (lg(x) - 1 >= Flm_CUP_LIMIT)
    d = Flm_det_CUP(x, p);
  else
    d = Flm_det_gauss(inplace? x: Flm_copy(x), p);
  return gc_ulong(av, d);
}

ulong
Flm_det_sp(GEN x, ulong p) { return Flm_det_i(x, p, 1); }
ulong
Flm_det(GEN x, ulong p) { return Flm_det_i(x, p, 0); }

/* Destroy x */
static GEN
Flm_gauss_pivot(GEN x, ulong p, long *rr)
{
  GEN c,d;
  long i,j,k,r,t,n,m;

  n=lg(x)-1; if (!n) { *rr=0; return NULL; }

  m=nbrows(x); r=0;
  d=cgetg(n+1,t_VECSMALL);
  c = zero_zv(m);
  for (k=1; k<=n; k++)
  {
    for (j=1; j<=m; j++)
      if (!c[j])
      {
        ucoeff(x,j,k) %= p;
        if (ucoeff(x,j,k)) break;
      }
    if (j>m) { r++; d[k]=0; }
    else
    {
      ulong piv = p - Fl_inv(ucoeff(x,j,k), p);
      c[j]=k; d[k]=j;
      for (i=k+1; i<=n; i++)
        ucoeff(x,j,i) = Fl_mul(piv, ucoeff(x,j,i), p);
      for (t=1; t<=m; t++)
        if (!c[t]) /* no pivot on that line yet */
        {
          piv = ucoeff(x,t,k);
          if (piv)
          {
            ucoeff(x,t,k) = 0;
            for (i=k+1; i<=n; i++)
              ucoeff(x,t,i) = Fl_add(ucoeff(x,t,i),
                                     Fl_mul(piv,ucoeff(x,j,i),p),p);
          }
        }
      for (i=k; i<=n; i++) ucoeff(x,j,i) = 0; /* dummy */
    }
  }
  *rr = r; return gc_const((pari_sp)d, d);
}

static GEN
Flm_pivots_CUP(GEN x, ulong p, long *rr)
{
  long i, n = lg(x) - 1, r;
  GEN R, C, U, P, d = zero_zv(n);
  ulong pi = get_Fl_red(p);
  r = Flm_CUP_pre(x, &R, &C, &U, &P, p, pi);
  for(i = 1; i <= r; i++)
    d[P[i]] = R[i];
  *rr = n - r; return gc_const((pari_sp)d, d);
}

GEN
Flm_pivots(GEN x, ulong p, long *rr, long inplace)
{
  if (lg(x) - 1 >= Flm_CUP_LIMIT && nbrows(x) >= Flm_CUP_LIMIT)
    return Flm_pivots_CUP(x, p, rr);
  return Flm_gauss_pivot(inplace? x: Flm_copy(x), p, rr);
}

long
Flm_rank(GEN x, ulong p)
{
  pari_sp av = avma;
  long r;
  if (lg(x) - 1 >= Flm_CUP_LIMIT && nbrows(x) >= Flm_CUP_LIMIT) {
    GEN R, C;
    ulong pi = get_Fl_red(p);
    return gc_long(av, Flm_echelon_pre(x, &R, &C, p, pi));
  }
  (void) Flm_pivots(x, p, &r, 0);
  return gc_long(av, lg(x)-1 - r);
}

/* assume dim A >= 1, A invertible + upper triangular, 1s on diagonal,
 * reduced mod p */
static GEN
Flm_inv_upper_1_ind(GEN A, long index, ulong p)
{
  long n = lg(A)-1, i = index, j;
  GEN u = const_vecsmall(n, 0);
  u[i] = 1;
  if (SMALL_ULONG(p))
    for (i--; i>0; i--)
    {
      ulong m = ucoeff(A,i,i+1) * uel(u,i+1); /* j = i+1 */
      for (j=i+2; j<=n; j++)
      {
        if (m & HIGHMASK) m %= p;
        m += ucoeff(A,i,j) * uel(u,j);
      }
      u[i] = Fl_neg(m % p, p);
    }
  else
    for (i--; i>0; i--)
    {
      ulong m = Fl_mul(ucoeff(A,i,i+1),uel(u,i+1), p); /* j = i+1 */
      for (j=i+2; j<=n; j++) m = Fl_add(m, Fl_mul(ucoeff(A,i,j),uel(u,j),p), p);
      u[i] = Fl_neg(m, p);
    }
  return u;
}
static GEN
Flm_inv_upper_1(GEN A, ulong p)
{
  long i, l;
  GEN B = cgetg_copy(A, &l);
  for (i = 1; i < l; i++) gel(B,i) = Flm_inv_upper_1_ind(A, i, p);
  return B;
}
/* destroy a, b */
static GEN
Flm_gauss_sp_OK(GEN a, GEN b, ulong *detp, ulong p)
{
  long i, j, k, li, bco, aco = lg(a)-1, s = 1;
  ulong det = 1;
  GEN u;

  li = nbrows(a);
  bco = lg(b)-1;
  for (i=1; i<=aco; i++)
  {
    ulong invpiv;
    /* Fl_get_col wants 0 <= a[i,j] < p for all i,j */
    for (k = 1; k < i; k++) ucoeff(a,k,i) %= p;
    for (k = i; k <= li; k++)
    {
      ulong piv = ( ucoeff(a,k,i) %= p );
      if (piv)
      {
        ucoeff(a,k,i) = Fl_inv(piv, p);
        if (detp)
        {
          if (det & HIGHMASK) det %= p;
          det *= piv;
        }
        break;
      }
    }
    /* found a pivot on line k */
    if (k > li) return NULL;
    if (k != i)
    { /* swap lines so that k = i */
      s = -s;
      for (j=i; j<=aco; j++) swap(gcoeff(a,i,j), gcoeff(a,k,j));
      for (j=1; j<=bco; j++) swap(gcoeff(b,i,j), gcoeff(b,k,j));
    }
    if (i == aco) break;

    invpiv = p - ucoeff(a,i,i); /* -1/piv mod p */
    for (k=i+1; k<=li; k++)
    {
      ulong m = ( ucoeff(a,k,i) %= p );
      if (!m) continue;

      m = Fl_mul(m, invpiv, p);
      if (m == 1) {
        for (j=i+1; j<=aco; j++) _Fl_add_OK(gel(a,j),k,i, p);
        for (j=1;   j<=bco; j++) _Fl_add_OK(gel(b,j),k,i, p);
      } else {
        for (j=i+1; j<=aco; j++) _Fl_addmul_OK(gel(a,j),k,i,m, p);
        for (j=1;   j<=bco; j++) _Fl_addmul_OK(gel(b,j),k,i,m, p);
      }
    }
  }
  if (detp)
  {
    det %=  p;
    if (s < 0 && det) det = p - det;
    *detp = det;
  }
  u = cgetg(bco+1,t_MAT);
  for (j=1; j<=bco; j++) gel(u,j) = Fl_get_col_OK(a,gel(b,j), aco,p);
  return u;
}

/* destroy a, b */
static GEN
Flm_gauss_sp_i(GEN a, GEN b, ulong *detp, ulong p)
{
  long i, j, k, li, bco, aco = lg(a)-1, s = 1;
  ulong det = 1;
  GEN u;
  ulong pi;
  if (!aco) { if (detp) *detp = 1; return cgetg(1,t_MAT); }
  if (SMALL_ULONG(p)) return Flm_gauss_sp_OK(a, b, detp, p);
  pi = get_Fl_red(p);
  li = nbrows(a);
  bco = lg(b)-1;
  for (i=1; i<=aco; i++)
  {
    ulong invpiv;
    /* Fl_get_col wants 0 <= a[i,j] < p for all i,j */
    for (k = i; k <= li; k++)
    {
      ulong piv = ucoeff(a,k,i);
      if (piv)
      {
        ucoeff(a,k,i) = Fl_inv(piv, p);
        if (detp) det = Fl_mul_pre(det, piv, p, pi);
        break;
      }
    }
    /* found a pivot on line k */
    if (k > li) return NULL;
    if (k != i)
    { /* swap lines so that k = i */
      s = -s;
      for (j=i; j<=aco; j++) swap(gcoeff(a,i,j), gcoeff(a,k,j));
      for (j=1; j<=bco; j++) swap(gcoeff(b,i,j), gcoeff(b,k,j));
    }
    if (i == aco) break;

    invpiv = p - ucoeff(a,i,i); /* -1/piv mod p */
    for (k=i+1; k<=li; k++)
    {
      ulong m = ucoeff(a,k,i);
      if (!m) continue;

      m = Fl_mul_pre(m, invpiv, p, pi);
      if (m == 1) {
        for (j=i+1; j<=aco; j++) _Fl_add(gel(a,j),k,i, p);
        for (j=1;   j<=bco; j++) _Fl_add(gel(b,j),k,i, p);
      } else {
        for (j=i+1; j<=aco; j++) _Fl_addmul(gel(a,j),k,i,m, p, pi);
        for (j=1;   j<=bco; j++) _Fl_addmul(gel(b,j),k,i,m, p, pi);
      }
    }
  }
  if (detp)
  {
    if (s < 0 && det) det = p - det;
    *detp = det;
  }
  u = cgetg(bco+1,t_MAT);
  for (j=1; j<=bco; j++) gel(u,j) = Fl_get_col(a,gel(b,j), aco,p);
  return u;
}

static GEN
Flm_gauss_from_CUP(GEN b, GEN R, GEN C, GEN U, GEN P, ulong p, ulong pi, ulong *detp)
{
  GEN Y = Flm_rsolve_lower_unit_pre(rowpermute(C, R), rowpermute(b, R), p, pi);
  GEN X = rowpermute(Flm_rsolve_upper_pre(U, Y, p, pi), perm_inv(P));
  if (detp)
  {
    ulong d = perm_sign(P) == 1? 1: p-1;
    long i, r = lg(R);
    for (i = 1; i < r; i++)
      d = Fl_mul_pre(d, ucoeff(U, i, i), p, pi);
    *detp = d;
  }
  return X;
}

static GEN
Flm_gauss_CUP(GEN a, GEN b, ulong *detp, ulong p) {
  GEN R, C, U, P;
  long n = lg(a) - 1, r;
  ulong pi = get_Fl_red(p);
  if (nbrows(a) < n || (r = Flm_CUP_pre(a, &R, &C, &U, &P, p, pi)) < n)
    return NULL;
  return Flm_gauss_from_CUP(b, R, C, U, P, p, pi, detp);
}

GEN
Flm_gauss_sp(GEN a, GEN b, ulong *detp, ulong p) {
  pari_sp av = avma;
  GEN x;
  if (lg(a) - 1 >= Flm_CUP_LIMIT)
    x = Flm_gauss_CUP(a, b, detp, p);
  else
    x = Flm_gauss_sp_i(a, b, detp, p);
  if (!x) return gc_NULL(av);
  return gerepileupto(av, x);
}

GEN
Flm_gauss(GEN a, GEN b, ulong p) {
  pari_sp av = avma;
  GEN x;
  if (lg(a) - 1 >= Flm_CUP_LIMIT)
    x = Flm_gauss_CUP(a, b, NULL, p);
  else {
    a = RgM_shallowcopy(a);
    b = RgM_shallowcopy(b);
    x = Flm_gauss_sp(a, b, NULL, p);
  }
  if (!x) return gc_NULL(av);
  return gerepileupto(av, x);
}

static GEN
Flm_inv_i(GEN a, ulong *detp, ulong p, long inplace) {
  pari_sp av = avma;
  long n = lg(a) - 1;
  GEN b, x;
  if (n == 0) return cgetg(1, t_MAT);
  b = matid_Flm(nbrows(a));
  if (n >= Flm_CUP_LIMIT)
    x = Flm_gauss_CUP(a, b, detp, p);
  else {
    if (!inplace)
      a = RgM_shallowcopy(a);
    x = Flm_gauss_sp(a, b, detp, p);
  }
  if (!x) return gc_NULL(av);
  return gerepileupto(av, x);
}

GEN
Flm_inv_sp(GEN a, ulong *detp, ulong p) {
  return Flm_inv_i(a, detp, p, 1);
}

GEN
Flm_inv(GEN a, ulong p) {
  return Flm_inv_i(a, NULL, p, 0);
}

GEN
Flm_Flc_gauss(GEN a, GEN b, ulong p) {
  pari_sp av = avma;
  GEN z = Flm_gauss(a, mkmat(b), p);
  if (!z) return gc_NULL(av);
  if (lg(z) == 1) { set_avma(av); return cgetg(1,t_VECSMALL); }
  return gerepileuptoleaf(av, gel(z,1));
}

GEN
Flm_adjoint(GEN A, ulong p)
{
  pari_sp av = avma;
  GEN R, C, U, P, C1, U1, v, c, d;
  long r, i, q, n = lg(A)-1, m;
  ulong D;
  ulong pi = get_Fl_red(p);
  if (n == 0) return cgetg(1, t_MAT);
  r = Flm_CUP_pre(A, &R, &C, &U, &P, p, pi);
  m = nbrows(A);
  if (r == n)
  {
    GEN X = Flm_gauss_from_CUP(matid_Flm(m), R, C, U, P, p, pi, &D);
    return gerepileupto(av, Flm_Fl_mul_pre(X, D, p, pi));
  }
  if (r < n-1) return zero_Flm(n, m);
  for (q = n, i = 1; i < n; i++)
    if (R[i] != i) { q = i; break; }
  C1 = matslice(C, 1, q-1, 1, q-1);
  c = vecslice(Flm_row(C, q), 1, q-1);
  c = Flm_lsolve_lower_unit_pre(C1, Flm_transpose(mkmat(c)), p, pi);
  d = cgetg(m+1, t_VECSMALL);
  for (i=1; i<q; i++)    uel(d,i) = ucoeff(c,1,i);
  uel(d,q) = p-1;
  for (i=q+1; i<=m; i++) uel(d,i) = 0;
  U1 = vecslice(U, 1, n-1);
  v = gel(Flm_rsolve_upper_pre(U1, mkmat(gel(U,n)), p, pi),1);
  v = vecsmall_append(v, p-1);
  D = perm_sign(P) != (odd(q+n)?-1:1) ? p-1 : 1;
  for (i = 1; i <= n-1; i++)
    D = Fl_mul_pre(D, ucoeff(U1, i, i), p, pi);
  d = Flv_Fl_mul(d, D, p);
  return rowpermute(Flc_Flv_mul(v, d, p),perm_inv(P));
}

static GEN
Flm_invimage_CUP(GEN A, GEN B, ulong p) {
  pari_sp av = avma;
  GEN R, Rc, C, U, P, B1, B2, C1, C2, X, Y, Z;
  long r;
  ulong pi = get_Fl_red(p);
  r = Flm_CUP_pre(A, &R, &C, &U, &P, p, pi);
  Rc = indexcompl(R, nbrows(B));
  C1 = rowpermute(C, R);
  C2 = rowpermute(C, Rc);
  B1 = rowpermute(B, R);
  B2 = rowpermute(B, Rc);
  Z = Flm_rsolve_lower_unit_pre(C1, B1, p, pi);
  if (!gequal(Flm_mul_pre(C2, Z, p, pi), B2))
    return NULL;
  Y = vconcat(Flm_rsolve_upper_pre(vecslice(U, 1, r), Z, p, pi),
              zero_Flm(lg(A) - 1 - r, lg(B) - 1));
  X = rowpermute(Y, perm_inv(P));
  return gerepileupto(av, X);
}

GEN
Flm_invimage_i(GEN A, GEN B, ulong p)
{
  GEN d, x, X, Y;
  long i, j, nY, nA = lg(A)-1, nB = lg(B)-1;

  if (!nB) return cgetg(1, t_MAT);
  if (maxss(nA, nB) >= Flm_CUP_LIMIT && nbrows(B) >= Flm_CUP_LIMIT)
    return Flm_invimage_CUP(A, B, p);

  x = Flm_ker(shallowconcat(Flm_neg(A,p), B), p);
  /* AX = BY, Y in strict upper echelon form with pivots = 1.
   * We must find T such that Y T = Id_nB then X T = Z. This exists iff
   * Y has at least nB columns and full rank */
  nY = lg(x)-1;
  if (nY < nB) return NULL;
  Y = rowslice(x, nA+1, nA+nB); /* nB rows */
  d = cgetg(nB+1, t_VECSMALL);
  for (i = nB, j = nY; i >= 1; i--, j--)
  {
    for (; j>=1; j--)
      if (coeff(Y,i,j)) { d[i] = j; break; }
    if (!j) return NULL;
  }
  /* reduce to the case Y square, upper triangular with 1s on diagonal */
  Y = vecpermute(Y, d);
  x = vecpermute(x, d);
  X = rowslice(x, 1, nA);
  return Flm_mul(X, Flm_inv_upper_1(Y,p), p);
}
GEN
Flm_invimage(GEN A, GEN B, ulong p)
{
  pari_sp av = avma;
  GEN X = Flm_invimage_i(A,B,p);
  if (!X) return gc_NULL(av);
  return gerepileupto(av, X);
}

GEN
Flm_Flc_invimage(GEN A, GEN y, ulong p)
{
  pari_sp av = avma;
  long i, l = lg(A);
  GEN M, x;
  ulong t;

  if (l==1) return NULL;
  if (lg(y) != lgcols(A)) pari_err_DIM("Flm_Flc_invimage");
  M = cgetg(l+1,t_MAT);
  for (i=1; i<l; i++) gel(M,i) = gel(A,i);
  gel(M,l) = y; M = Flm_ker(M,p);
  i = lg(M)-1; if (!i) return gc_NULL(av);

  x = gel(M,i); t = x[l];
  if (!t) return gc_NULL(av);

  setlg(x,l); t = Fl_inv(Fl_neg(t,p),p);
  if (t!=1) x = Flv_Fl_mul(x, t, p);
  return gerepileuptoleaf(av, x);
}
