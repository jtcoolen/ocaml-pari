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

/* Self-Initializing Multi-Polynomial Quadratic Sieve, based on code developed
 * as part of the LiDIA project.
 *
 * Original version: Thomas Papanikolaou and Xavier Roblot
 * Extensively modified by The PARI group. */
/* Notation commonly used in this file, and sketch of algorithm:
 *
 * Given an odd integer N > 1 to be factored, we throw in a small odd squarefree
 * multiplier k so as to make kN = 1 mod 4 and to have many small primes over
 * which X^2 - kN splits.  We compute a factor base FB of such primes then
 * look for values x0 such that Q0(x0) = x0^2 - kN can be decomposed over FB,
 * up to a possible factor dividing k and a possible "large prime". Relations
 * involving the latter can be combined into full relations which don't; full
 * relations, by Gaussian elimination over F2 for the exponent vectors lead us
 * to an expression X^2 - Y^2 divisible by N and hopefully to a nontrivial
 * splitting when we compute gcd(X + Y, N).  Note that this can never
 * split prime powers.
 *
 * Candidates x0 are found by sieving along arithmetic progressions modulo the
 * small primes in FB and evaluation of candidates picks out those x0 where
 * many of these progressions coincide, resulting in a highly divisible Q0(x0).
 *
 * The Multi-Polynomial version improves this by choosing a modest subset of
 * FB primes (let A be their product) and forcing these to divide Q0(x).
 * Write Q(x) = Q0(2Ax + B) = (2Ax + B)^2 - kN = 4A(Ax^2 + Bx + C), where B is
 * suitably chosen.  For each A, there are 2^omega_A possible values for B
 * but we'll use only half of these, since the other half is easily covered by
 * exploiting the symmetry x -> -x of the original Q0. The "Self-Initializating"
 * bit refers to the fact that switching from one B to the next is fast, whereas
 * switching to the next A involves some recomputation (C is never needed).
 * Thus we quickly run through many polynomials sharing the same A.
 *
 * The sieve ranges over values x0 such that |x0| < M  (we use x = x0 + M
 * as array subscript).  The coefficients A are chosen so that A*M ~ sqrt(kN).
 * Then |B| is bounded by ~ (j+4)*A, and |C| = -C ~ (M/4)*sqrt(kN), so
 * Q(x0)/(4A) takes values roughly between -|C| and 3|C|.
 *
 * Refinements. We do not use the smallest FB primes for sieving, incorporating
 * them only after selecting candidates).  The substitution of 2Ax+B into
 * X^2 - kN, with odd B, forces 2 to occur; when kN is 1 mod 8, it occurs at
 * least to the 3rd power; when kN = 5 mod 8, it occurs exactly to the 2nd
 * power.  We never sieve on 2 and always pull out the power of 2 directly. The
 * prime factors of k show up whenever 2Ax + B has a factor in common with k;
 * we don't sieve on these either but easily recognize them in a candidate. */

#include "paricfg.h"
#ifdef HAS_SSE2
#include <emmintrin.h>
#endif

#include "pari.h"
#include "paripriv.h"

#define DEBUGLEVEL DEBUGLEVEL_mpqs

/** DEBUG **/
/* #define MPQS_DEBUG_VERBOSE 1 */
#include "mpqs.h"

#define REL_OFFSET 20
#define REL_MASK ((1UL<<REL_OFFSET)-1)
#define MAX_PE_PAIR 60

#ifdef HAS_SSE2
#define EXT0(a) ((ulong)__builtin_ia32_vec_ext_v2di((__v2di)(a), 0))
#define EXT1(a) ((ulong)__builtin_ia32_vec_ext_v2di((__v2di)(a), 1))
#define TEST(a) (EXT0(a) || EXT1(a))
typedef __v2di mpqs_bit_array;
const mpqs_bit_array mpqs_mask = { (long) 0x8080808080808080L, (long) 0x8080808080808080UL };
#else
/* Use ulong for the bit arrays */
typedef ulong mpqs_bit_array;
#define TEST(a) (a)

#ifdef LONG_IS_64BIT
const mpqs_bit_array mpqs_mask = 0x8080808080808080UL;
#else
const mpqs_bit_array mpqs_mask = 0x80808080UL;
#endif
#endif

static GEN rel_q(GEN c) { return gel(c,3); }
static GEN rel_Y(GEN c) { return gel(c,1); }
static GEN rel_p(GEN c) { return gel(c,2); }

static void
frel_add(hashtable *frel, GEN R)
{
  ulong h = hash_GEN(R);
  if (!hash_search2(frel, (void*)R, h))
    hash_insert2(frel, (void*)R, (void*)1, h);
}

/*********************************************************************/
/**                         INITIAL SIZING                          **/
/*********************************************************************/
/* # of decimal digits of argument */
static long
decimal_len(GEN N)
{ pari_sp av = avma; return gc_long(av, 1+logint(N, utoipos(10))); }

/* To be called after choosing k and putting kN into the handle:
 * Pick up the parameters for given size of kN in decimal digits and fill in
 * the handle. Return 0 when kN is too large, 1 when we're ok. */
static int
mpqs_set_parameters(mpqs_handle_t *h)
{
  long s, D;
  const mpqs_parameterset_t *P;

  h->digit_size_kN = D = decimal_len(h->kN);
  if (D > MPQS_MAX_DIGIT_SIZE_KN) return 0;
  P = &(mpqs_parameters[maxss(0, D - 9)]);
  h->tolerance   = P->tolerance;
  h->lp_scale    = P->lp_scale;
  /* make room for prime factors of k if any: */
  h->size_of_FB  = s = P->size_of_FB + h->_k->omega_k;
  /* for the purpose of Gauss elimination etc., prime factors of k behave
   * like real FB primes, so take them into account when setting the goal: */
  h->target_rels = (s >= 200 ? s + 10 : (mpqs_int32_t)(s * 1.05));
  h->M           = P->M;
  h->omega_A     = P->omega_A;
  h->no_B        = 1UL << (P->omega_A - 1);
  h->pmin_index1 = P->pmin_index1;
  /* certain subscripts into h->FB should also be offset by omega_k: */
  h->index0_FB   = 3 + h->_k->omega_k;
  if (DEBUGLEVEL >= 5)
  {
    err_printf("MPQS: kN = %Ps\n", h->kN);
    err_printf("MPQS: kN has %ld decimal digits\n", D);
    err_printf("\t(estimated memory needed: %4.1fMBy)\n",
               (s + 1)/8388608. * h->target_rels);
  }
  return 1;
}

/*********************************************************************/
/**                       OBJECT HOUSEKEEPING                       **/
/*********************************************************************/

/* factor base constructor. Really a home-grown memalign(3c) underneath.
 * We don't want FB entries to straddle L1 cache line boundaries, and
 * malloc(3c) only guarantees alignment adequate for all primitive data
 * types of the platform ABI - typically to 8 or 16 byte boundaries.
 * Also allocate the inv_A_H array.
 * The FB array pointer is returned for convenience */
static mpqs_FB_entry_t *
mpqs_FB_ctor(mpqs_handle_t *h)
{
  /* leave room for slots 0, 1, and sentinel slot at the end of the array */
  long size_FB_chunk = (h->size_of_FB + 3) * sizeof(mpqs_FB_entry_t);
  /* like FB, except this one does not have a sentinel slot at the end */
  long size_IAH_chunk = (h->size_of_FB + 2) * sizeof(mpqs_inv_A_H_t);
  char *fbp = (char*)stack_malloc(size_FB_chunk + 64);
  char *iahp = (char*)stack_malloc(size_IAH_chunk + 64);
  long fbl, iahl;

  h->FB_chunk = (void *)fbp;
  h->invAH_chunk = (void *)iahp;
  /* round up to next higher 64-bytes-aligned address */
  fbl = (((long)fbp) + 64) & ~0x3FL;
  /* and put the actual array there */
  h->FB = (mpqs_FB_entry_t *)fbl;

  iahl = (((long)iahp) + 64) & ~0x3FL;
  h->inv_A_H = (mpqs_inv_A_H_t *)iahl;
  return (mpqs_FB_entry_t *)fbl;
}

/* sieve array constructor;  also allocates the candidates array
 * and temporary storage for relations under construction */
static void
mpqs_sieve_array_ctor(mpqs_handle_t *h)
{
  long size = (h->M << 1) + 1;
  mpqs_int32_t size_of_FB = h->size_of_FB;

  h->sieve_array = (unsigned char *) stack_calloc_align(size, sizeof(mpqs_mask));
  h->sieve_array_end = h->sieve_array + size - 2;
  h->sieve_array_end[1] = 255; /* sentinel */
  h->candidates = (long *)stack_malloc(MPQS_CANDIDATE_ARRAY_SIZE * sizeof(long));
  /* whereas mpqs_self_init() uses size_of_FB+1, we just use the size as
   * it is, not counting FB[1], to start off the following estimate */
  if (size_of_FB > MAX_PE_PAIR) size_of_FB = MAX_PE_PAIR;
  /* and for tracking which primes occur in the current relation: */
  h->relaprimes = (long *) stack_malloc((size_of_FB << 1) * sizeof(long));
}

/* allocate GENs for current polynomial and self-initialization scratch data */
static void
mpqs_poly_ctor(mpqs_handle_t *h)
{
  mpqs_int32_t i, w = h->omega_A;
  h->per_A_pr = (mpqs_per_A_prime_t *)
                stack_calloc(w * sizeof(mpqs_per_A_prime_t));
  /* A is the product of w primes, each below word size.
   * |B| <= (w + 4) * A, so can have at most one word more
   * H holds residues modulo A: the same size as used for A is sufficient. */
  h->A = cgeti(w + 2);
  h->B = cgeti(w + 3);
  for (i = 0; i < w; i++) h->per_A_pr[i]._H = cgeti(w + 2);
}

/*********************************************************************/
/**                        FACTOR BASE SETUP                        **/
/*********************************************************************/
/* fill in the best-guess multiplier k for N. We force kN = 1 mod 4.
 * Caller should proceed to fill in kN
 * See Knuth-Schroeppel function in
 * Robert D. Silverman
 * The multiple polynomial quadratic sieve
 * Math. Comp. 48 (1987), 329-339
 * https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866119-8/
 */
static ulong
mpqs_find_k(mpqs_handle_t *h)
{
  const pari_sp av = avma;
  const long N_mod_8 = mod8(h->N), N_mod_4 = N_mod_8 & 3;
  long dl = decimal_len(h->N);
  long D = maxss(0, minss(dl,MPQS_MAX_DIGIT_SIZE_KN)-9);
  long MPQS_MULTIPLIER_SEARCH_DEPTH = mpqs_parameters[D].size_of_FB;
  forprime_t S;
  struct {
    const mpqs_multiplier_t *_k;
    long np; /* number of primes in factorbase so far for this k */
    double value; /* the larger, the better */
  } cache[MPQS_POSSIBLE_MULTIPLIERS];
  ulong MPQS_NB_MULTIPLIERS = dl < 40 ? 5 : MPQS_POSSIBLE_MULTIPLIERS;
  ulong p, i, nbk;

  for (i = nbk = 0; i < numberof(cand_multipliers); i++)
  {
    const mpqs_multiplier_t *cand_k = &cand_multipliers[i];
    long k = cand_k->k;
    double v;
    if ((k & 3) != N_mod_4) continue; /* want kN = 1 (mod 4) */
    v = -log((double)k)/2;
    if ((k & 7) == N_mod_8) v += M_LN2; /* kN = 1 (mod 8) */
    cache[nbk].np = 0;
    cache[nbk]._k = cand_k;
    cache[nbk].value = v;
    if (++nbk == MPQS_NB_MULTIPLIERS) break; /* enough */
  }
  /* next test is an impossible situation: kills spurious gcc-5.1 warnings
   * "array subscript is above array bounds" */
  if (nbk > MPQS_POSSIBLE_MULTIPLIERS) nbk = MPQS_POSSIBLE_MULTIPLIERS;
  u_forprime_init(&S, 2, ULONG_MAX);
  while ( (p = u_forprime_next(&S)) )
  {
    long kroNp = kroiu(h->N, p), seen = 0;
    if (!kroNp) return p;
    for (i = 0; i < nbk; i++)
    {
      long krokp;
      if (cache[i].np > MPQS_MULTIPLIER_SEARCH_DEPTH) continue;
      seen++;
      krokp = krouu(cache[i]._k->k % p, p);
      if (krokp == kroNp) /* kronecker(k*N, p)=1 */
      {
        cache[i].value += 2*log((double) p)/p;
        cache[i].np++;
      } else if (krokp == 0)
      {
        cache[i].value += log((double) p)/p;
        cache[i].np++;
      }
    }
    if (!seen) break; /* we're gone through SEARCH_DEPTH primes for all k */
  }
  if (!p) pari_err_OVERFLOW("mpqs_find_k [ran out of primes]");
  {
    long best_i = 0;
    double v = cache[0].value;
    for (i = 1; i < nbk; i++)
      if (cache[i].value > v) { best_i = i; v = cache[i].value; }
    h->_k = cache[best_i]._k; return gc_ulong(av,0);
  }
}

/* Create a factor base of 'size' primes p_i such that legendre(k*N, p_i) != -1
 * We could have shifted subscripts down from their historical arrangement,
 * but this seems too risky for the tiny potential gain in memory economy.
 * The real constraint is that the subscripts of anything which later shows
 * up at the Gauss stage must be nonnegative, because the exponent vectors
 * there use the same subscripts to refer to the same FB entries.  Thus in
 * particular, the entry representing -1 could be put into FB[0], but could
 * not be moved to FB[-1] (although mpqs_FB_ctor() could be easily adapted
 * to support negative subscripts).-- The historically grown layout is:
 * FB[0] is unused.
 * FB[1] is not explicitly used but stands for -1.
 * FB[2] contains 2 (always).
 * Before we are called, the size_of_FB field in the handle will already have
 * been adjusted by _k->omega_k, so there's room for the primes dividing k,
 * which when present will occupy FB[3] and following.
 * The "real" odd FB primes begin at FB[h->index0_FB].
 * FB[size_of_FB+1] is the last prime p_i.
 * FB[size_of_FB+2] is a sentinel to simplify some of our loops.
 * Thus we allocate size_of_FB+3 slots for FB.
 *
 * If a prime factor of N is found during the construction, it is returned
 * in f, otherwise f = 0. */

/* returns the FB array pointer for convenience */
static mpqs_FB_entry_t *
mpqs_create_FB(mpqs_handle_t *h, ulong *f)
{
  mpqs_FB_entry_t *FB = mpqs_FB_ctor(h);
  const pari_sp av = avma;
  mpqs_int32_t size = h->size_of_FB;
  long i;
  mpqs_uint32_t k = h->_k->k;
  forprime_t S;

  h->largest_FB_p = 0; /* -Wall */
  FB[2].fbe_p = 2;
  /* the fbe_logval and the fbe_sqrt_kN for 2 are never used */
  FB[2].fbe_flags = MPQS_FBE_CLEAR;
  for (i = 3; i < h->index0_FB; i++)
  { /* this loop executes h->_k->omega_k = 0, 1, or 2 times */
    mpqs_uint32_t kp = (ulong)h->_k->kp[i-3];
    if (MPQS_DEBUGLEVEL >= 7) err_printf(",<%lu>", (ulong)kp);
    FB[i].fbe_p = kp;
    /* we could flag divisors of k here, but no need so far */
    FB[i].fbe_flags = MPQS_FBE_CLEAR;
    FB[i].fbe_flogp = (float)log2((double) kp);
    FB[i].fbe_sqrt_kN = 0;
  }
  (void)u_forprime_init(&S, 3, ULONG_MAX);
  while (i < size + 2)
  {
    ulong p = u_forprime_next(&S);
    if (p > k || k % p)
    {
      ulong kNp = umodiu(h->kN, p);
      long kr = krouu(kNp, p);
      if (kr != -1)
      {
        if (kr == 0) { *f = p; return FB; }
        FB[i].fbe_p = (mpqs_uint32_t) p;
        FB[i].fbe_flags = MPQS_FBE_CLEAR;
        /* dyadic logarithm of p; single precision suffices */
        FB[i].fbe_flogp = (float)log2((double)p);
        /* cannot yet fill in fbe_logval because the scaling multiplier
         * depends on the largest prime in FB, as yet unknown */

        /* x such that x^2 = kN (mod p_i) */
        FB[i++].fbe_sqrt_kN = (mpqs_uint32_t)Fl_sqrt(kNp, p);
      }
    }
  }
  set_avma(av);
  if (MPQS_DEBUGLEVEL >= 7)
  {
    err_printf("MPQS: FB [-1,2");
    for (i = 3; i < h->index0_FB; i++) err_printf(",<%lu>", FB[i].fbe_p);
    for (; i < size + 2; i++) err_printf(",%lu", FB[i].fbe_p);
    err_printf("]\n");
  }

  FB[i].fbe_p = 0;              /* sentinel */
  h->largest_FB_p = FB[i-1].fbe_p; /* at subscript size_of_FB + 1 */

  /* locate the smallest prime that will be used for sieving */
  for (i = h->index0_FB; FB[i].fbe_p != 0; i++)
    if (FB[i].fbe_p >= h->pmin_index1) break;
  h->index1_FB = i;
  /* with our parameters this will never fall off the end of the FB */
  *f = 0; return FB;
}

/*********************************************************************/
/**                      MISC HELPER FUNCTIONS                      **/
/*********************************************************************/

/* Effect of the following:  multiplying the base-2 logarithm of some
 * quantity by log_multiplier will rescale something of size
 *    log2 ( sqrt(kN) * M / (largest_FB_prime)^tolerance )
 * to 232.  Note that sqrt(kN) * M is just A*M^2, the value our polynomials
 * take at the outer edges of the sieve interval.  The scale here leaves
 * a little wiggle room for accumulated rounding errors from the approximate
 * byte-sized scaled logarithms for the factor base primes which we add up
 * in the sieving phase.-- The threshold is then chosen so that a point in
 * the sieve has to reach a result which, under the same scaling, represents
 *    log2 ( sqrt(kN) * M / (largest_FB_prime)^tolerance )
 * in order to be accepted as a candidate. */
/* The old formula was...
 *   log_multiplier =
 *      127.0 / (0.5 * log2 (handle->dkN) + log2((double)M)
 *               - tolerance * log2((double)handle->largest_FB_p));
 * and we used to use this with a constant threshold of 128. */

/* NOTE: We used to divide log_multiplier by an extra factor 2, and in
 * compensation we were multiplying by 2 when the fbe_logp fields were being
 * filled in, making all those bytes even.  Tradeoff: the extra bit of
 * precision is helpful, but interferes with a possible sieving optimization
 * (artificially shift right the logp's of primes in A, and just run over both
 * arithmetical progressions  (which coincide in this case)  instead of
 * skipping the second one, to avoid the conditional branch in the
 * mpqs_sieve() loops).  We could still do this, but might lose a little bit
 * accuracy for those primes.  Probably no big deal. */
static void
mpqs_set_sieve_threshold(mpqs_handle_t *h)
{
  mpqs_FB_entry_t *FB = h->FB;
  double log_maxval, log_multiplier;
  long i;

  h->l2sqrtkN = 0.5 * log2(h->dkN);
  h->l2M = log2((double)h->M);
  log_maxval = h->l2sqrtkN + h->l2M - MPQS_A_FUDGE;
  log_multiplier = 232.0 / log_maxval;
  h->sieve_threshold = (unsigned char) (log_multiplier *
    (log_maxval - h->tolerance * log2((double)h->largest_FB_p))) + 1;
  /* That "+ 1" really helps - we may want to tune towards somewhat smaller
   * tolerances  (or introduce self-tuning one day)... */

  /* If this turns out to be <128, scream loudly.
   * That means that the FB or the tolerance or both are way too
   * large for the size of kN.  (Normally, the threshold should end
   * up in the 150...170 range.) */
  if (h->sieve_threshold < 128) {
    h->sieve_threshold = 128;
    pari_warn(warner,
        "MPQS: sizing out of tune, FB size or tolerance\n\ttoo large");
  }
  if (DEBUGLEVEL >= 5)
    err_printf("MPQS: sieve threshold: %ld\n",h->sieve_threshold);
  /* Now fill in the byte-sized approximate scaled logarithms of p_i */
  if (DEBUGLEVEL >= 5)
    err_printf("MPQS: computing logarithm approximations for p_i in FB\n");
  for (i = h->index0_FB; i < h->size_of_FB + 2; i++)
    FB[i].fbe_logval = (unsigned char) (log_multiplier * FB[i].fbe_flogp);
}

/* Given the partially populated handle, find the optimum place in the FB
 * to pick prime factors for A from.  The lowest admissible subscript is
 * index0_FB, but unless kN is very small, we stay away a bit from that.
 * The highest admissible is size_of_FB + 1, where the largest FB prime
 * resides.  The ideal corner is about (sqrt(kN)/M) ^ (1/omega_A),
 * so that A will end up of size comparable to sqrt(kN)/M;  experimentally
 * it seems desirable to stay slightly below this.  Moreover, the selection
 * of the individual primes happens to err on the large side, for which we
 * compensate a bit, using the (small positive) quantity MPQS_A_FUDGE.
 * We rely on a few auxiliary fields in the handle to be already set by
 * mqps_set_sieve_threshold() before we are called.
 * Return 1 on success, and 0 otherwise. */
static int
mpqs_locate_A_range(mpqs_handle_t *h)
{
  /* i will be counted up to the desirable index2_FB + 1, and omega_A is never
   * less than 3, and we want
   *   index2_FB - (omega_A - 1) + 1 >= index0_FB + omega_A - 3,
   * so: */
  long i = h->index0_FB + 2*(h->omega_A) - 4;
  double l2_target_pA;
  mpqs_FB_entry_t *FB = h->FB;

  h->l2_target_A = (h->l2sqrtkN - h->l2M - MPQS_A_FUDGE);
  l2_target_pA = h->l2_target_A / h->omega_A;

  /* find the sweet spot, normally shouldn't take long */
  while (FB[i].fbe_p && FB[i].fbe_flogp <= l2_target_pA) i++;

  /* check whether this hasn't walked off the top end... */
  /* The following should actually NEVER happen. */
  if (i > h->size_of_FB - 3)
  { /* this isn't going to work at all. */
    pari_warn(warner,
        "MPQS: sizing out of tune, FB too small or\n\tway too few primes in A");
    return 0;
  }
  h->index2_FB = i - 1; return 1;
  /* assert: index0_FB + (omega_A - 3) [the lowest FB subscript used in primes
   * for A]  + (omega_A - 2) <= index2_FB  [the subscript from which the choice
   * of primes for A starts, putting omega_A - 1 of them at or below index2_FB,
   * and the last and largest one above, cf. mpqs_si_choose_primes]. Moreover,
   * index2_FB indicates the last prime below the ideal size, unless (when kN
   * is tiny) the ideal size was too small to use. */
}

/*********************************************************************/
/**                       SELF-INITIALIZATION                       **/
/*********************************************************************/

#ifdef MPQS_DEBUG
/* Debug-only helper routine: check correctness of the root z mod p_i
 * by evaluating A * z^2 + B * z + C mod p_i  (which should be 0). */
static void
check_root(mpqs_handle_t *h, GEN mC, long p, long start)
{
  pari_sp av = avma;
  long z = start - ((long)(h->M) % p);
  if (umodiu(subii(mulsi(z, addii(h->B, mulsi(z, h->A))), mC), p))
  {
    err_printf("MPQS: p = %ld\n", p);
    err_printf("MPQS: A = %Ps\n", h->A);
    err_printf("MPQS: B = %Ps\n", h->B);
    err_printf("MPQS: C = %Ps\n", negi(mC));
    err_printf("MPQS: z = %ld\n", z);
    pari_err_BUG("MPQS: self_init: found wrong polynomial");
  }
  set_avma(av);
}
#endif

/* Increment *x > 0 to a larger value which has the same number of 1s in its
 * binary representation.  Wraparound can be detected by the caller as long as
 * we keep total_no_of_primes_for_A strictly less than BITS_IN_LONG.
 *
 * Changed switch to increment *x in all cases to the next larger number
 * which (a) has the same count of 1 bits and (b) does not arise from the
 * old value by moving a single 1 bit one position to the left  (which was
 * undesirable for the sieve). --GN based on discussion with TP */
INLINE void
mpqs_increment(mpqs_uint32_t *x)
{
  mpqs_uint32_t r1_mask, r01_mask, slider=1UL;

  switch (*x & 0x1F)
  { /* 32-way computed jump handles 22 out of 32 cases */
  case 29:
    (*x)++; break; /* shifts a single bit, but we postprocess this case */
  case 26:
    (*x) += 2; break; /* again */
  case 1: case 3: case 6: case 9: case 11:
  case 17: case 19: case 22: case 25: case 27:
    (*x) += 3; return;
  case 20:
    (*x) += 4; break; /* again */
  case 5: case 12: case 14: case 21:
    (*x) += 5; return;
  case 2: case 7: case 13: case 18: case 23:
    (*x) += 6; return;
  case 10:
    (*x) += 7; return;
  case 8:
    (*x) += 8; break; /* and again */
  case 4: case 15:
    (*x) += 12; return;
  default: /* 0, 16, 24, 28, 30, 31 */
    /* isolate rightmost 1 */
    r1_mask = ((*x ^ (*x - 1)) + 1) >> 1;
    /* isolate rightmost 1 which has a 0 to its left */
    r01_mask = ((*x ^ (*x + r1_mask)) + r1_mask) >> 2;
    /* simple cases.  Both of these shift a single bit one position to the
       left, and will need postprocessing */
    if (r1_mask == r01_mask) { *x += r1_mask; break; }
    if (r1_mask == 1) { *x += r01_mask; break; }
    /* General case: add r01_mask, kill off as many 1 bits as possible to its
     * right while at the same time filling in 1 bits from the LSB. */
    if (r1_mask == 2) { *x += (r01_mask>>1) + 1; return; }
    while (r01_mask > r1_mask && slider < r1_mask)
    {
      r01_mask >>= 1; slider <<= 1;
    }
    *x += r01_mask + slider - 1;
    return;
  }
  /* post-process cases which couldn't be finalized above */
  r1_mask = ((*x ^ (*x - 1)) + 1) >> 1;
  r01_mask = ((*x ^ (*x + r1_mask)) + r1_mask) >> 2;
  if (r1_mask == r01_mask) { *x += r1_mask; return; }
  if (r1_mask == 1) { *x += r01_mask; return; }
  if (r1_mask == 2) { *x += (r01_mask>>1) + 1; return; }
  while (r01_mask > r1_mask && slider < r1_mask)
  {
    r01_mask >>= 1; slider <<= 1;
  }
  *x += r01_mask + slider - 1;
}

/* self-init (1): advancing the bit pattern, and choice of primes for A.
 * On first call, h->bin_index = 0. On later occasions, we need to begin
 * by clearing the MPQS_FBE_DIVIDES_A bit in the fbe_flags of the former
 * prime factors of A (use per_A_pr to find them). Upon successful return, that
 * array will have been filled in, and the flag bits will have been turned on
 * again in the right places.
 * Return 1 when all is fine and 0 when we found we'd be using more bits to
 * the left in bin_index than we have matching primes in the FB. In the latter
 * case, bin_index will be zeroed out, index2_FB will be incremented by 2,
 * index2_moved will be turned on; the caller, after checking that index2_FB
 * has not become too large, should just call us again, which then succeeds:
 * we'll start again with a right-justified sequence of 1 bits in bin_index,
 * now interpreted as selecting primes relative to the new index2_FB. */
INLINE int
mpqs_si_choose_primes(mpqs_handle_t *h)
{
  mpqs_FB_entry_t *FB = h->FB;
  mpqs_per_A_prime_t *per_A_pr = h->per_A_pr;
  double l2_last_p = h->l2_target_A;
  mpqs_int32_t omega_A = h->omega_A;
  int i, j, v2, prev_last_p_idx;
  int room = h->index2_FB - h->index0_FB - omega_A + 4;
  /* The notion of room here (cf mpqs_locate_A_range() above) is the number
   * of primes at or below index2_FB which are eligible for A. We need
   * >= omega_A - 1 of them, and it is guaranteed by mpqs_locate_A_range() that
   * this many are available: the lowest FB slot used for A is never less than
   * index0_FB + omega_A - 3. When omega_A = 3 (very small kN), we allow
   * ourselves to reach all the way down to index0_FB; otherwise, we keep away
   * from it by at least one position.  For omega_A >= 4 this avoids situations
   * where the selection of the smaller primes here has advanced to a lot of
   * very small ones, and the single last larger one has soared away to bump
   * into the top end of the FB. */
  mpqs_uint32_t room_mask;
  mpqs_int32_t p;
  ulong bits;

  /* XXX also clear the index_j field here? */
  if (h->bin_index == 0)
  { /* first time here, or after increasing index2_FB, initialize to a pattern
     * of omega_A - 1 consecutive 1 bits. Caller has ensured that there are
     * enough primes for this in the FB below index2_FB. */
    h->bin_index = (1UL << (omega_A - 1)) - 1;
    prev_last_p_idx = 0;
  }
  else
  { /* clear out old flags */
    for (i = 0; i < omega_A; i++) MPQS_FLG(i) = MPQS_FBE_CLEAR;
    prev_last_p_idx = MPQS_I(omega_A-1);

    if (room > 30) room = 30;
    room_mask = ~((1UL << room) - 1);

    /* bump bin_index to next acceptable value. If index2_moved is off, call
     * mpqs_increment() once; otherwise, repeat until there's something in the
     * least significant 2 bits - to ensure that we never re-use an A which
     * we'd used before increasing index2_FB - but also stop if something shows
     * up in the forbidden bits on the left where we'd run out of bits or walk
     * beyond index0_FB + omega_A - 3. */
    mpqs_increment(&h->bin_index);
    if (h->index2_moved)
    {
      while ((h->bin_index & (room_mask | 0x3)) == 0)
        mpqs_increment(&h->bin_index);
    }
    /* did we fall off the edge on the left? */
    if ((h->bin_index & room_mask) != 0)
    { /* Yes. Turn on the index2_moved flag in the handle */
      h->index2_FB += 2; /* caller to check this isn't too large!!! */
      h->index2_moved = 1;
      h->bin_index = 0;
      if (MPQS_DEBUGLEVEL >= 5)
        err_printf("MPQS: wrapping, more primes for A now chosen near FB[%ld] = %ld\n",
                   (long)h->index2_FB,
                   (long)FB[h->index2_FB].fbe_p);
      return 0; /* back off - caller should retry */
    }
  }
  /* assert: we aren't occupying any of the room_mask bits now, and if
   * index2_moved had already been on, at least one of the two LSBs is on */
  bits = h->bin_index;
  if (MPQS_DEBUGLEVEL >= 6)
    err_printf("MPQS: new bit pattern for primes for A: 0x%lX\n", bits);

  /* map bits to FB subscripts, counting downward with bit 0 corresponding
   * to index2_FB, and accumulate logarithms against l2_last_p */
  j = h->index2_FB;
  v2 = vals((long)bits);
  if (v2) { j -= v2; bits >>= v2; }
  for (i = omega_A - 2; i >= 0; i--)
  {
    MPQS_I(i) = j;
    l2_last_p -= MPQS_LP(i);
    MPQS_FLG(i) |= MPQS_FBE_DIVIDES_A;
    bits &= ~1UL;
    if (!bits) break; /* i = 0 */
    v2 = vals((long)bits); /* > 0 */
    bits >>= v2; j -= v2;
  }
  /* Choose the larger prime.  Note we keep index2_FB <= size_of_FB - 3 */
  for (j = h->index2_FB + 1; (p = FB[j].fbe_p); j++)
    if (FB[j].fbe_flogp > l2_last_p) break;
  /* The following trick avoids generating a large proportion of duplicate
   * relations when the last prime falls into an area where there are large
   * gaps from one FB prime to the next, and would otherwise often be repeated
   * (so that successive A's would wind up too similar to each other). While
   * this trick isn't perfect, it gets rid of a major part of the potential
   * duplication. */
  if (p && j == prev_last_p_idx) { j++; p = FB[j].fbe_p; }
  MPQS_I(omega_A - 1) = p? j: h->size_of_FB + 1;
  MPQS_FLG(omega_A - 1) |= MPQS_FBE_DIVIDES_A;

  if (MPQS_DEBUGLEVEL >= 6)
  {
    err_printf("MPQS: chose primes for A");
    for (i = 0; i < omega_A; i++)
      err_printf(" FB[%ld]=%ld%s", (long)MPQS_I(i), (long)MPQS_AP(i),
                 i < omega_A - 1 ? "," : "\n");
  }
  return 1;
}

/* There are 4 parts to self-initialization, exercised at different times:
 * - choosing a new sqfree coef. A (selecting its prime factors, FB bookkeeping)
 * - doing the actual computations attached to a new A
 * - choosing a new B keeping the same A (much simpler)
 * - a small common bit that needs to happen in both cases.
 * As to the first item, the scheme works as follows: pick omega_A - 1 prime
 * factors for A below the index2_FB point which marks their ideal size, and
 * one prime above this point, choosing the latter so log2(A) ~ l2_target_A.
 * Lower prime factors are chosen using bit patterns of constant weight,
 * gradually moving away from index2_FB towards smaller FB subscripts.
 * If this bumps into index0_FB (for very small input), back up by increasing
 * index2_FB by two, and from then on choosing only bit patterns with either or
 * both of their bottom bits set, so at least one of the omega_A - 1 smaller
 * prime factor will be beyond the original index2_FB point. In this way we
 * avoid re-using the same A. (The choice of the upper "flyer" prime is
 * constrained by the size of the FB, which normally should never a problem.
 * For tiny kN, we might have to live with a nonoptimal choice.)
 *
 * Mathematically, we solve a quadratic (over F_p for each prime p in the FB
 * which doesn't divide A), a linear equation for each prime p | A, and
 * precompute differences between roots mod p so we can adjust the roots
 * quickly when we change B. See Thomas Sosnowski's Diplomarbeit. */
/* compute coefficients of sieving polynomial for self initializing variant.
 * Coefficients A and B are set (preallocated GENs) and several tables are
 * updated. */
static int
mpqs_self_init(mpqs_handle_t *h)
{
  const ulong size_of_FB = h->size_of_FB + 1;
  mpqs_FB_entry_t *FB = h->FB;
  mpqs_inv_A_H_t *inv_A_H = h->inv_A_H;
  const pari_sp av = avma;
  GEN p1, A = h->A, B = h->B;
  mpqs_per_A_prime_t *per_A_pr = h->per_A_pr;
  long i, j;

#ifdef MPQS_DEBUG
  err_printf("MPQS DEBUG: enter self init, avma = 0x%lX\n", (ulong)avma);
#endif
  if (++h->index_j == (mpqs_uint32_t)h->no_B)
  { /* all the B's have been used, choose new A; this is indicated by setting
     * index_j to 0 */
    h->index_j = 0;
    h->index_i++; /* count finished A's */
  }

  if (h->index_j == 0)
  { /* compute first polynomial with new A */
    GEN a, b, A2;
    if (!mpqs_si_choose_primes(h))
    { /* Ran out of room towards small primes, and index2_FB was raised. */
      if (size_of_FB - h->index2_FB < 4) return 0; /* Fail */
      (void)mpqs_si_choose_primes(h); /* now guaranteed to succeed */
    }
    /* bin_index and per_A_pr now populated with consistent values */

    /* compute A = product of omega_A primes given by bin_index */
    a = b = NULL;
    for (i = 0; i < h->omega_A; i++)
    {
      ulong p = MPQS_AP(i);
      a = a? muliu(a, p): utoipos(p);
    }
    affii(a, A);
    /* Compute H[i], 0 <= i < omega_A.  Also compute the initial
     * B = sum(v_i*H[i]), by taking all v_i = +1
     * TODO: following needs to be changed later for segmented FB and sieve
     * interval, where we'll want to precompute several B's. */
    for (i = 0; i < h->omega_A; i++)
    {
      ulong p = MPQS_AP(i);
      GEN t = divis(A, (long)p);
      t = remii(mulii(t, muluu(Fl_inv(umodiu(t, p), p), MPQS_SQRT(i))), A);
      affii(t, MPQS_H(i));
      b = b? addii(b, t): t;
    }
    affii(b, B); set_avma(av);

    /* ensure B = 1 mod 4 */
    if (mod2(B) == 0)
      affii(addii(B, mului(mod4(A), A)), B); /* B += (A % 4) * A; */

    A2 = shifti(A, 1);
    /* compute the roots z1, z2, of the polynomial Q(x) mod p_j and
     * initialize start1[i] with the first value p_i | Q(z1 + i p_j)
     * initialize start2[i] with the first value p_i | Q(z2 + i p_j)
     * The following loop does The Right Thing for primes dividing k (where
     * sqrt_kN is 0 mod p). Primes dividing A are skipped here, and are handled
     * further down in the common part of SI. */
    for (j = 3; (ulong)j <= size_of_FB; j++)
    {
      ulong s, mb, t, m, p, iA2, iA;
      if (FB[j].fbe_flags & MPQS_FBE_DIVIDES_A) continue;
      p = (ulong)FB[j].fbe_p;
      m = h->M % p;
      iA2 = Fl_inv(umodiu(A2, p), p); /* = 1/(2*A) mod p_j */
      iA = iA2 << 1; if (iA > p) iA -= p;
      mb = umodiu(B, p); if (mb) mb = p - mb; /* mb = -B mod p */
      s = FB[j].fbe_sqrt_kN;
      t = Fl_add(m, Fl_mul(Fl_sub(mb, s, p), iA2, p), p);
      FB[j].fbe_start1 = (mpqs_int32_t)t;
      FB[j].fbe_start2 = (mpqs_int32_t)Fl_add(t, Fl_mul(s, iA, p), p);
      for (i = 0; i < h->omega_A - 1; i++)
      {
        ulong h = umodiu(MPQS_H(i), p);
        MPQS_INV_A_H(i,j) = Fl_mul(h, iA, p); /* 1/A * H[i] mod p_j */
      }
    }
  }
  else
  { /* no "real" computation -- use recursive formula */
    /* The following exploits that B is the sum of omega_A terms +-H[i]. Each
     * time we switch to a new B, we choose a new pattern of signs; the
     * precomputation of the inv_A_H array allows us to change the two
     * arithmetic progressions equally fast. The choice of sign patterns does
     * not follow the bit pattern of the ordinal number of B in the current
     * cohort; rather, we use a Gray code, changing only one sign each time.
     * When the i-th rightmost bit of the new ordinal number index_j of B is 1,
     * the sign of H[i] is changed; the next bit to the left tells us whether
     * we should be adding or subtracting the difference term. We never need to
     * change the sign of H[omega_A-1] (the topmost one), because that would
     * just give us the same sieve items Q(x) again with the opposite sign
     * of x.  This is why we only precomputed inv_A_H up to i = omega_A - 2. */
    ulong p, v2 = vals(h->index_j); /* new starting positions for sieving */
    j = h->index_j >> v2;
    p1 = shifti(MPQS_H(v2), 1);
    if (j & 2)
    { /* j = 3 mod 4 */
      for (j = 3; (ulong)j <= size_of_FB; j++)
      {
        if (FB[j].fbe_flags & MPQS_FBE_DIVIDES_A) continue;
        p = (ulong)FB[j].fbe_p;
        FB[j].fbe_start1 = Fl_sub(FB[j].fbe_start1, MPQS_INV_A_H(v2,j), p);
        FB[j].fbe_start2 = Fl_sub(FB[j].fbe_start2, MPQS_INV_A_H(v2,j), p);
      }
      p1 = addii(B, p1);
    }
    else
    { /* j = 1 mod 4 */
      for (j = 3; (ulong)j <= size_of_FB; j++)
      {
        if (FB[j].fbe_flags & MPQS_FBE_DIVIDES_A) continue;
        p = (ulong)FB[j].fbe_p;
        FB[j].fbe_start1 = Fl_add(FB[j].fbe_start1, MPQS_INV_A_H(v2,j), p);
        FB[j].fbe_start2 = Fl_add(FB[j].fbe_start2, MPQS_INV_A_H(v2,j), p);
      }
      p1 = subii(B, p1);
    }
    affii(p1, B);
  }

  /* p=2 is a special case.  start1[2], start2[2] are never looked at,
   * so don't bother setting them. */

  /* compute zeros of polynomials that have only one zero mod p since p | A */
  p1 = diviiexact(subii(h->kN, sqri(B)), shifti(A, 2)); /* coefficient -C */
  for (i = 0; i < h->omega_A; i++)
  {
    ulong p = MPQS_AP(i), s = h->M + Fl_div(umodiu(p1, p), umodiu(B, p), p);
    FB[MPQS_I(i)].fbe_start1 = FB[MPQS_I(i)].fbe_start2 = (mpqs_int32_t)(s % p);
  }
#ifdef MPQS_DEBUG
  for (j = 3; j <= size_of_FB; j++)
  {
    check_root(h, p1, FB[j].fbe_p, FB[j].fbe_start1);
    check_root(h, p1, FB[j].fbe_p, FB[j].fbe_start2);
  }
#endif
  if (MPQS_DEBUGLEVEL >= 6)
    err_printf("MPQS: chose Q_%ld(x) = %Ps x^2 %c %Ps x + C\n",
               (long) h->index_j, h->A,
               signe(h->B) < 0? '-': '+', absi_shallow(h->B));
  set_avma(av);
#ifdef MPQS_DEBUG
  err_printf("MPQS DEBUG: leave self init, avma = 0x%lX\n", (ulong)avma);
#endif
  return 1;
}

/*********************************************************************/
/**                           THE SIEVE                             **/
/*********************************************************************/
/* p4 = 4*p, logp ~ log(p), B/E point to the beginning/end of a sieve array */
INLINE void
mpqs_sieve_p(unsigned char *B, unsigned char *E, long p4, long p,
             unsigned char logp)
{
  unsigned char *e = E - p4;
  /* Unrolled loop. It might be better to let the compiler worry about this
   * kind of optimization, based on its knowledge of whatever useful tricks the
   * machine instruction set architecture is offering */
  while (e - B >= 0) /* signed comparison */
  {
    (*B) += logp, B += p;
    (*B) += logp, B += p;
    (*B) += logp, B += p;
    (*B) += logp, B += p;
  }
  while (E - B >= 0) (*B) += logp, B += p;
}

INLINE void
mpqs_sieve_p1(unsigned char *B, unsigned char *E, long s1, long s2,
             unsigned char logp)
{
  while (E - B >= 0)
  {
    (*B) += logp, B += s1;
    if (E - B < 0) break;
    (*B) += logp, B += s2;
  }
}

INLINE void
mpqs_sieve_p2(unsigned char *B, unsigned char *E, long p4, long s1, long s2,
             unsigned char logp)
{
  unsigned char *e = E - p4;
  /* Unrolled loop. It might be better to let the compiler worry about this
   * kind of optimization, based on its knowledge of whatever useful tricks the
   * machine instruction set architecture is offering */
  while (e - B >= 0) /* signed comparison */
  {
    (*B) += logp, B += s1;
    (*B) += logp, B += s2;
    (*B) += logp, B += s1;
    (*B) += logp, B += s2;
    (*B) += logp, B += s1;
    (*B) += logp, B += s2;
    (*B) += logp, B += s1;
    (*B) += logp, B += s2;
  }
  while (E - B >= 0) {(*B) += logp, B += s1; if (E - B < 0) break; (*B) += logp, B += s2;}
}
static void
mpqs_sieve(mpqs_handle_t *h)
{
  long p, l = h->index1_FB;
  mpqs_FB_entry_t *FB = &(h->FB[l]);
  unsigned char *S = h->sieve_array, *Send = h->sieve_array_end;
  long size = h->M << 1, size4 = size >> 3;
  memset((void*)S, 0, size * sizeof(unsigned char));
  for (  ; (p = FB->fbe_p) && p <= size4; FB++) /* l++ */
  {
    unsigned char logp = FB->fbe_logval;
    long s1 = FB->fbe_start1, s2 = FB->fbe_start2;
    /* sieve with FB[l] from start1[l], and from start2[l] if s1 != s2 */
    if (s1 == s2) mpqs_sieve_p(S + s1, Send, p << 2, p, logp);
    else
    {
      if (s1>s2) lswap(s1,s2)
      mpqs_sieve_p2(S + s1, Send, p << 2, s2-s1,p+s1-s2, logp);
    }
  }
  for (   ; (p = FB->fbe_p) && p <= size; FB++) /* l++ */
  {
    unsigned char logp = FB->fbe_logval;
    long s1 = FB->fbe_start1, s2 = FB->fbe_start2;
    /* sieve with FB[l] from start1[l], and from start2[l] if s1 != s2 */
    if (s1 == s2) mpqs_sieve_p(S + s1, Send, p << 2, p, logp);
    else
    {
      if (s1>s2) lswap(s1,s2)
      mpqs_sieve_p1(S + s1, Send, s2-s1, p+s1-s2, logp);
    }
  }
  for (    ; (p = FB->fbe_p); FB++)
  {
    unsigned char logp = FB->fbe_logval;
    long s1 = FB->fbe_start1, s2 = FB->fbe_start2;
    if (s1 < size) S[s1] += logp;
    if (s2!=s1 && s2 < size) S[s2] += logp;
  }
}

/* Could use the fact that 4 | M, but let the compiler worry about unrolling. */
static long
mpqs_eval_sieve(mpqs_handle_t *h)
{
  long x = 0, count = 0, M2 = h->M << 1;
  unsigned char t = h->sieve_threshold;
  unsigned char *S = h->sieve_array;
  mpqs_bit_array * U = (mpqs_bit_array *) S;
  long *cand = h->candidates;
  const long sizemask = sizeof(mpqs_mask);

  /* Exploiting the sentinel, we don't need to check for x < M2 in the inner
   * while loop; more than makes up for the lack of explicit unrolling. */
  while (count < MPQS_CANDIDATE_ARRAY_SIZE - 1)
  {
    long j, y;
    while (!TEST(U[x]&mpqs_mask)) x++;
    y = x*sizemask;
    for (j=0; j<sizemask; j++, y++)
    {
      if (y >= M2)
        { cand[count] = 0; return count; }
      if (S[y]>=t)
        cand[count++] = y;
    }
    x++;
  }
  cand[count] = 0; return count;
}

/*********************************************************************/
/**                     CONSTRUCTING RELATIONS                      **/
/*********************************************************************/

/* only used for debugging */
static void
split_relp(GEN rel, GEN *prelp, GEN *prelc)
{
  long j, l = lg(rel);
  GEN relp, relc;
  *prelp = relp = cgetg(l, t_VECSMALL);
  *prelc = relc = cgetg(l, t_VECSMALL);
  for (j=1; j<l; j++)
  {
    relc[j] = rel[j] >> REL_OFFSET;
    relp[j] = rel[j] & REL_MASK;
  }
}

#ifdef MPQS_DEBUG
static GEN
mpqs_factorback(mpqs_handle_t *h, GEN relp)
{
  GEN N = h->N, Q = gen_1;
  long j, l = lg(relp);
  for (j = 1; j < l; j++)
  {
    long e = relp[j] >> REL_OFFSET, i = relp[j] & REL_MASK;
    if (i == 1) Q = Fp_neg(Q,N); /* special case -1 */
    else Q = Fp_mul(Q, Fp_powu(utoipos(h->FB[i].fbe_p), e, N), N);
  }
  return Q;
}
static void
mpqs_check_rel(mpqs_handle_t *h, GEN c)
{
  pari_sp av = avma;
  int LP = (lg(c) == 4);
  GEN rhs = mpqs_factorback(h, rel_p(c));
  GEN Y = rel_Y(c), Qx_2 = remii(sqri(Y), h->N);
  if (LP) rhs = modii(mulii(rhs, rel_q(c)), h->N);
  if (!equalii(Qx_2, rhs))
  {
    GEN relpp, relpc;
    split_relp(rel_p(c), &relpp, &relpc);
    err_printf("MPQS: %Ps : %Ps %Ps\n", Y, relpp,relpc);
    err_printf("\tQx_2 = %Ps\n", Qx_2);
    err_printf("\t rhs = %Ps\n", rhs);
    pari_err_BUG(LP? "MPQS: wrong large prime relation found"
                   : "MPQS: wrong full relation found");
  }
  PRINT_IF_VERBOSE(LP? "\b(;)": "\b(:)");
  set_avma(av);
}
#endif

static void
rel_to_ei(GEN ei, GEN relp)
{
  long j, l = lg(relp);
  for (j=1; j<l; j++)
  {
    long e = relp[j] >> REL_OFFSET, i = relp[j] & REL_MASK;
    ei[i] += e;
  }
}
static void
mpqs_add_factor(GEN relp, long *i, ulong ei, ulong pi)
{ relp[++*i] = pi | (ei << REL_OFFSET); }

static int
zv_is_even(GEN V)
{
  long i, l = lg(V);
  for (i=1; i<l; i++)
    if (odd(uel(V,i))) return 0;
  return 1;
}

static GEN
combine_large_primes(mpqs_handle_t *h, GEN rel1, GEN rel2)
{
  GEN new_Y, new_Y1, Y1 = rel_Y(rel1), Y2 = rel_Y(rel2);
  long l, lei = h->size_of_FB + 1, nb = 0;
  GEN ei, relp, iq, q = rel_q(rel1);

  if (!invmod(q, h->N, &iq)) return equalii(iq, h->N)? NULL: iq; /* rare */
  ei = zero_zv(lei);
  rel_to_ei(ei, rel_p(rel1));
  rel_to_ei(ei, rel_p(rel2));
  if (zv_is_even(ei)) return NULL;
  new_Y = modii(mulii(mulii(Y1, Y2), iq), h->N);
  new_Y1 = subii(h->N, new_Y);
  if (abscmpii(new_Y1, new_Y) < 0) new_Y = new_Y1;
  relp = cgetg(MAX_PE_PAIR+1,t_VECSMALL);
  if (odd(ei[1])) mpqs_add_factor(relp, &nb, 1, 1);
  for (l = 2; l <= lei; l++)
    if (ei[l]) mpqs_add_factor(relp, &nb, ei[l],l);
  setlg(relp, nb+1);
  if (DEBUGLEVEL >= 6)
  {
    GEN relpp, relpc, rel1p, rel1c, rel2p, rel2c;
    split_relp(relp,&relpp,&relpc);
    split_relp(rel1,&rel1p,&rel1c);
    split_relp(rel2,&rel2p,&rel2c);
    err_printf("MPQS: combining\n");
    err_printf("    {%Ps @ %Ps : %Ps}\n", q, Y1, rel1p, rel1c);
    err_printf("  * {%Ps @ %Ps : %Ps}\n", q, Y2, rel2p, rel2c);
    err_printf(" == {%Ps, %Ps}\n", relpp, relpc);
  }
#ifdef MPQS_DEBUG
  {
    pari_sp av1 = avma;
    if (!equalii(modii(sqri(new_Y), h->N), mpqs_factorback(h, relp)))
      pari_err_BUG("MPQS: combined large prime relation is false");
    set_avma(av1);
  }
#endif
  return mkvec2(new_Y, relp);
}

/* nc candidates */
static GEN
mpqs_eval_cand(mpqs_handle_t *h, long nc, hashtable *frel, hashtable *lprel)
{
  mpqs_FB_entry_t *FB = h->FB;
  GEN A = h->A, B = h->B;
  long *relaprimes = h->relaprimes, *candidates = h->candidates;
  long pi, i;
  int pii;
  mpqs_per_A_prime_t *per_A_pr = h->per_A_pr;

  for (i = 0; i < nc; i++)
  {
    pari_sp btop = avma;
    GEN Qx, Qx_part, Y, relp = cgetg(MAX_PE_PAIR+1,t_VECSMALL);
    long powers_of_2, p, x = candidates[i], nb = 0;
    int relaprpos = 0;
    long k;
    unsigned char thr = h->sieve_array[x];
    /* Y = 2*A*x + B, Qx = Y^2/(4*A) = Q(x) */
    Y = addii(mulis(A, 2 * (x - h->M)), B);
    Qx = subii(sqri(Y), h->kN); /* != 0 since N not a square and (N,k) = 1 */
    if (signe(Qx) < 0)
    {
      setabssign(Qx);
      mpqs_add_factor(relp, &nb, 1, 1); /* i = 1, ei = 1, pi */
    }
    /* Qx > 0, divide by powers of 2; we're really dealing with 4*A*Q(x), so we
     * always have at least 2^2 here, and at least 2^3 when kN = 1 mod 4 */
    powers_of_2 = vali(Qx);
    Qx = shifti(Qx, -powers_of_2);
    mpqs_add_factor(relp, &nb, powers_of_2, 2); /* i = 1, ei = 1, pi */
    /* When N is small, it may happen that N | Qx outright. In any case, when
     * no extensive prior trial division / Rho / ECM was attempted, gcd(Qx,N)
     * may turn out to be a nontrivial factor of N (not in FB or we'd have
     * found it already, but possibly smaller than the large prime bound). This
     * is too rare to check for here in the inner loop, but it will be caught
     * if such an LP relation is ever combined with another. */

    /* Pass 1 over odd primes in FB: pick up all possible divisors of Qx
     * including those sitting in k or in A, and remember them in relaprimes.
     * Do not yet worry about possible repeated factors, these will be found in
     * the Pass 2. Pass 1 recognizes divisors of A by their corresponding flags
     * bit in the FB entry. (Divisors of k are ignored at this stage.)
     * We construct a preliminary table of FB subscripts and "exponents" of FB
     * primes which divide Qx. (We store subscripts, not the primes themselves.)
     * We distinguish three cases:
     * 0) prime in A which does not divide Qx/A,
     * 1) prime not in A which divides Qx/A,
     * 2) prime in A which divides Qx/A.
     * Cases 1 and 2 need checking for repeated factors, kind 0 doesn't.
     * Cases 0 and 1 contribute 1 to the exponent in the relation, case 2
     * contributes 2.
     * Factors in common with k are simpler: if they occur, they occur
     * exactly to the first power, and this makes no difference in Pass 1,
     * so they behave just like every normal odd FB prime. */
    for (Qx_part = A, pi = 3; pi< h->index1_FB; pi++)
    {
      ulong p = FB[pi].fbe_p;
      long xp = x % p;
      /* Here we used that MPQS_FBE_DIVIDES_A = 1. */

      if (xp == FB[pi].fbe_start1 || xp == FB[pi].fbe_start2)
      { /* p divides Q(x)/A and possibly A, case 2 or 3 */
        ulong ei = FB[pi].fbe_flags & MPQS_FBE_DIVIDES_A;
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 1 + ei;
        Qx_part = muliu(Qx_part, p);
      }
    }
    for (  ; thr && (p = FB[pi].fbe_p); pi++)
    {
      long xp = x % p;
      /* Here we used that MPQS_FBE_DIVIDES_A = 1. */

      if (xp == FB[pi].fbe_start1 || xp == FB[pi].fbe_start2)
      { /* p divides Q(x)/A and possibly A, case 2 or 3 */
        ulong ei = FB[pi].fbe_flags & MPQS_FBE_DIVIDES_A;
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 1 + ei;
        Qx_part = muliu(Qx_part, p);
        thr -= FB[pi].fbe_logval;
      }
    }
    for (k = 0;  k< h->omega_A; k++)
    {
      long pi = MPQS_I(k);
      ulong p = FB[pi].fbe_p;
      long xp = x % p;
      if (!(xp == FB[pi].fbe_start1 || xp == FB[pi].fbe_start2))
      { /* p divides A but does not divide Q(x)/A, case 1 */
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 0;
      }
    }
    /* We have accumulated the known factors of Qx except for possible repeated
     * factors and for possible large primes.  Divide off what we have.
     * This is faster than dividing off A and each prime separately. */
    Qx = diviiexact(Qx, Qx_part);

#ifdef MPQS_DEBUG
    err_printf("MPQS DEBUG: eval loop 3, avma = 0x%lX\n", (ulong)avma);
#endif
    /* Pass 2: deal with repeated factors and store tentative relation. At this
     * point, the only primes which can occur again in the adjusted Qx are
     * those in relaprimes which are followed by 1 or 2. We must pick up those
     * followed by a 0, too. */
    PRINT_IF_VERBOSE("a");
    for (pii = 0; pii < relaprpos; pii += 2)
    {
      ulong r, ei = relaprimes[pii+1];
      GEN q;

      pi = relaprimes[pii];
      /* p | k (identified by its index before index0_FB)* or p | A (ei = 0) */
      if ((mpqs_int32_t)pi < h->index0_FB || ei == 0)
      {
        mpqs_add_factor(relp, &nb, 1, pi);
        continue;
      }
      p = FB[pi].fbe_p;
      /* p might still divide the current adjusted Qx. Try it. */
      switch(cmpiu(Qx, p))
      {
        case 0: ei++; Qx = gen_1; break;
        case 1:
          q = absdiviu_rem(Qx, p, &r);
          while (r == 0) { ei++; Qx = q; q = absdiviu_rem(Qx, p, &r); }
          break;
      }
      mpqs_add_factor(relp, &nb, ei, pi);
    }

#ifdef MPQS_DEBUG
    err_printf("MPQS DEBUG: eval loop 4, avma = 0x%lX\n", (ulong)avma);
#endif
    PRINT_IF_VERBOSE("\bb");
    setlg(relp, nb+1);
    if (is_pm1(Qx))
    {
      GEN rel = gerepilecopy(btop, mkvec2(absi_shallow(Y), relp));
#ifdef MPQS_DEBUG
      mpqs_check_rel(h, rel);
#endif
      frel_add(frel, rel);
    }
    else if (cmpiu(Qx, h->lp_bound) <= 0)
    {
      ulong q = itou(Qx);
      GEN rel = mkvec3(absi_shallow(Y),relp,Qx);
      GEN col = hash_haskey_GEN(lprel, (void*)q);
#ifdef MPQS_DEBUG
      mpqs_check_rel(h, rel);
#endif
      if (!col) /* relation up to large prime */
        hash_insert(lprel, (void*)q, (void*)gerepilecopy(btop,rel));
      else if ((rel = combine_large_primes(h, rel, col)))
      {
        if (typ(rel) == t_INT) return rel; /* very unlikely */
#ifdef MPQS_DEBUG
        mpqs_check_rel(h, rel);
#endif
        frel_add(frel, gerepilecopy(btop,rel));
      }
      else
        set_avma(btop);
    }
    else
    { /* TODO: check for double large prime */
      PRINT_IF_VERBOSE("\b.");
      set_avma(btop);
    }
  }
  PRINT_IF_VERBOSE("\n");
  return NULL;
}

/*********************************************************************/
/**                    FROM RELATIONS TO DIVISORS                   **/
/*********************************************************************/

/* create an F2m from a relations list */
static GEN
rels_to_F2Ms(GEN rel)
{
  long i, cols = lg(rel)-1;
  GEN m = cgetg(cols+1, t_VEC);
  for (i = 1; i <= cols; i++)
  {
    GEN relp = gmael(rel,i,2), rel2;
    long j, l = lg(relp), o = 0, k;
    for (j = 1; j < l; j++)
      if (odd(relp[j] >> REL_OFFSET)) o++;
    rel2 = cgetg(o+1, t_VECSMALL);
    for (j = 1, k = 1; j < l; j++)
      if (odd(relp[j] >> REL_OFFSET))
        rel2[k++] = relp[j] & REL_MASK;
    gel(m, i) = rel2;
  }
  return m;
}

static int
split(GEN *D, long *e)
{
  ulong mask;
  long flag;
  if (MR_Jaeschke(*D)) { *e = 1; return 1; } /* probable prime */
  if (Z_issquareall(*D, D))
  { /* squares could cost us a lot of time */
    if (DEBUGLEVEL >= 5) err_printf("MPQS: decomposed a square\n");
    *e = 2; return 1;
  }
  mask = 7;
  /* 5th/7th powers aren't worth the trouble. OTOH once we have the hooks for
   * dealing with cubes, higher powers can be handled essentially for free) */
  if ((flag = is_357_power(*D, D, &mask)))
  {
    if (DEBUGLEVEL >= 5)
      err_printf("MPQS: decomposed a %s power\n", uordinal(flag));
    *e = flag; return 1;
  }
  *e = 0; return 0; /* known composite */
}

/* return a GEN structure containing NULL but safe for gerepileupto */
static GEN
mpqs_solve_linear_system(mpqs_handle_t *h, hashtable *frel)
{
  mpqs_FB_entry_t *FB = h->FB;
  pari_sp av = avma;
  GEN rels = hash_keys(frel), N = h->N, r, c, res, ei, M, Ker;
  long i, j, nrows, rlast, rnext, rmax, rank;

  M = rels_to_F2Ms(rels);
  Ker = F2Ms_ker(M, h->size_of_FB+1); rank = lg(Ker)-1;
  if (DEBUGLEVEL >= 4)
  {
    if (DEBUGLEVEL >= 7)
      err_printf("\\\\ MPQS RELATION MATRIX\nFREL=%Ps\nKERNEL=%Ps\n",M, Ker);
    err_printf("MPQS: Gauss done: kernel has rank %ld, taking gcds...\n", rank);
  }
  if (!rank)
  { /* trivial kernel; main loop may look for more relations */
    if (DEBUGLEVEL >= 3)
      pari_warn(warner, "MPQS: no solutions found from linear system solver");
    return gc_NULL(av); /* no factors found */
  }

  /* Expect up to 2^rank pairwise coprime factors, but a kernel basis vector
   * may not contribute to the decomposition; r stores the factors and c what
   * we know about them (0: composite, 1: probably prime, >= 2: proper power) */
  ei = cgetg(h->size_of_FB + 2, t_VECSMALL);
  rmax = logint(N, utoi(3));
  if (rank <= BITS_IN_LONG-2)
    rmax = minss(rmax, 1L<<rank); /* max # of factors we can hope for */
  r = cgetg(rmax+1, t_VEC);
  c = cgetg(rmax+1, t_VECSMALL);
  rnext = rlast = 1;
  nrows = lg(M)-1;
  for (i = 1; i <= rank; i++)
  { /* loop over kernel basis */
    GEN X = gen_1, Y_prod = gen_1, X_plus_Y, D;
    pari_sp av2 = avma, av3;
    long done = 0; /* # probably-prime factors or powers whose bases we won't
                    * handle any further */
    memset((void *)(ei+1), 0, (h->size_of_FB + 1) * sizeof(long));
    for (j = 1; j <= nrows; j++)
      if (F2m_coeff(Ker, j, i))
      {
        GEN R = gel(rels,j);
        Y_prod = gerepileuptoint(av2, remii(mulii(Y_prod, gel(R,1)), N));
        rel_to_ei(ei, gel(R,2));
      }
    av3 = avma;
    for (j = 2; j <= h->size_of_FB + 1; j++)
      if (ei[j])
      {
        GEN q = utoipos(FB[j].fbe_p);
        if (ei[j] & 1) pari_err_BUG("MPQS (relation is a nonsquare)");
        X = remii(mulii(X, Fp_powu(q, (ulong)ei[j]>>1, N)), N);
        X = gerepileuptoint(av3, X);
      }
    if (MPQS_DEBUGLEVEL >= 1 && !dvdii(subii(sqri(X), sqri(Y_prod)), N))
    {
      err_printf("MPQS: X^2 - Y^2 != 0 mod N\n");
      err_printf("\tindex i = %ld\n", i);
      pari_warn(warner, "MPQS: wrong relation found after Gauss");
    }
    /* At this point, gcd(X-Y, N) * gcd(X+Y, N) = N:
     * 1) N | X^2 - Y^2, so it divides the LHS;
     * 2) let P be any prime factor of N. If P | X-Y and P | X+Y, then P | 2X
     * But X is a product of powers of FB primes => coprime to N.
     * Hence we work with gcd(X+Y, N) alone. */
    X_plus_Y = addii(X, Y_prod);
    if (rnext == 1)
    { /* we still haven't decomposed, and want both a gcd and its cofactor. */
      D = gcdii(X_plus_Y, N);
      if (is_pm1(D) || equalii(D,N)) { set_avma(av2); continue; }
      /* got something that works */
      if (DEBUGLEVEL >= 5)
        err_printf("MPQS: splitting N after %ld kernel vector%s\n",
                   i+1, (i? "s" : ""));
      gel(r,1) = diviiexact(N, D);
      gel(r,2) = D;
      rlast = rnext = 3;
      if (split(&gel(r,1), &c[1])) done++;
      if (split(&gel(r,2), &c[2])) done++;
      if (done == 2 || rmax == 2) break;
      if (DEBUGLEVEL >= 5)
        err_printf("MPQS: got two factors, looking for more...\n");
    }
    else
    { /* we already have factors */
      for (j = 1; j < rnext; j++)
      { /* loop over known-composite factors */
        /* skip probable primes and also roots of pure powers: they are a lot
         * smaller than N and should be easy to deal with later */
        if (c[j]) { done++; continue; }
        av3 = avma; D = gcdii(X_plus_Y, gel(r,j));
        if (is_pm1(D) || equalii(D, gel(r,j))) { set_avma(av3); continue; }
        /* got one which splits this factor */
        if (DEBUGLEVEL >= 5)
          err_printf("MPQS: resplitting a factor after %ld kernel vectors\n",
                     i+1);
        gel(r,j) = diviiexact(gel(r,j), D);
        gel(r,rnext) = D;
        if (split(&gel(r,j), &c[j])) done++;
        /* Don't increment done: happens later when we revisit c[rnext] during
         * the present inner loop. */
        (void)split(&gel(r,rnext), &c[rnext]);
        if (++rnext > rmax) break; /* all possible factors seen */
      } /* loop over known composite factors */

      if (rnext > rlast)
      {
        if (DEBUGLEVEL >= 5)
          err_printf("MPQS: got %ld factors%s\n", rlast - 1,
                     (done < rlast ? ", looking for more..." : ""));
        rlast = rnext;
      }
      /* break out if we have rmax factors or all current factors are probable
       * primes or tiny roots from pure powers */
      if (rnext > rmax || done == rnext - 1) break;
    }
  }
  if (rnext == 1) return gc_NULL(av); /* no factors found */

  /* normal case: convert to ifac format as described in ifactor1.c (value,
   * exponent, class [unknown, known composite, known prime]) */
  rlast = rnext - 1; /* # of distinct factors found */
  res = cgetg(3*rlast + 1, t_VEC);
  if (DEBUGLEVEL >= 6) err_printf("MPQS: wrapping up %ld factors\n", rlast);
  for (i = j = 1; i <= rlast; i++, j += 3)
  {
    long C  = c[i];
    icopyifstack(gel(r,i), gel(res,j)); /* factor */
    gel(res,j+1) = C <= 1? gen_1: utoipos(C); /* exponent */
    gel(res,j+2) = C ? NULL: gen_0; /* unknown or known composite */
    if (DEBUGLEVEL >= 6)
      err_printf("\tpackaging %ld: %Ps ^%ld (%s)\n", i, gel(r,i),
                 itos(gel(res,j+1)), (C? "unknown": "composite"));
  }
  return res;
}

/*********************************************************************/
/**               MAIN ENTRY POINT AND DRIVER ROUTINE               **/
/*********************************************************************/
static void
toolarge()
{ pari_warn(warner, "MPQS: number too big to be factored with MPQS,\n\tgiving up"); }

/* Factors N using the self-initializing multipolynomial quadratic sieve
 * (SIMPQS).  Returns one of the two factors, or (usually) a vector of factors
 * and exponents and information about which ones are still composite, or NULL
 * when we can't seem to make any headway. */
GEN
mpqs(GEN N)
{
  const long size_N = decimal_len(N);
  mpqs_handle_t H;
  GEN fact; /* will in the end hold our factor(s) */
  mpqs_FB_entry_t *FB; /* factor base */
  double dbg_target, DEFEAT;
  ulong p;
  pari_timer T;
  hashtable lprel, frel;
  pari_sp av = avma;

  if (DEBUGLEVEL >= 4) err_printf("MPQS: number to factor N = %Ps\n", N);
  if (size_N > MPQS_MAX_DIGIT_SIZE_KN) { toolarge(); return NULL; }
  if (DEBUGLEVEL >= 4)
  {
    timer_start(&T);
    err_printf("MPQS: factoring number of %ld decimal digits\n", size_N);
  }
  H.N = N;
  H.bin_index = 0;
  H.index_i = 0;
  H.index_j = 0;
  H.index2_moved = 0;
  p = mpqs_find_k(&H);
  if (p) return gc_utoipos(av,p);
  if (DEBUGLEVEL >= 5)
    err_printf("MPQS: found multiplier %ld for N\n", H._k->k);
  H.kN = muliu(N, H._k->k);
  if (!mpqs_set_parameters(&H)) { toolarge(); return NULL; }

  if (DEBUGLEVEL >= 5)
    err_printf("MPQS: creating factor base and allocating arrays...\n");
  FB = mpqs_create_FB(&H, &p);
  if (p) return gc_utoipos(av, p);
  mpqs_sieve_array_ctor(&H);
  mpqs_poly_ctor(&H);

  H.lp_bound = minss(H.largest_FB_p, MPQS_LP_BOUND);
  /* don't allow large primes to have room for two factors both bigger than
   * what the FB contains (...yet!) */
  H.lp_bound *= minss(H.lp_scale, H.largest_FB_p - 1);
  H.dkN = gtodouble(H.kN);
  /* compute the threshold and fill in the byte-sized scaled logarithms */
  mpqs_set_sieve_threshold(&H);
  if (!mpqs_locate_A_range(&H)) return NULL;
  if (DEBUGLEVEL >= 4)
  {
    err_printf("MPQS: sieving interval = [%ld, %ld]\n", -(long)H.M, (long)H.M);
    /* that was a little white lie, we stop one position short at the top */
    err_printf("MPQS: size of factor base = %ld\n", (long)H.size_of_FB);
    err_printf("MPQS: striving for %ld relations\n", (long)H.target_rels);
    err_printf("MPQS: coefficients A will be built from %ld primes each\n",
               (long)H.omega_A);
    err_printf("MPQS: primes for A to be chosen near FB[%ld] = %ld\n",
               (long)H.index2_FB, (long)FB[H.index2_FB].fbe_p);
    err_printf("MPQS: smallest prime used for sieving FB[%ld] = %ld\n",
               (long)H.index1_FB, (long)FB[H.index1_FB].fbe_p);
    err_printf("MPQS: largest prime in FB = %ld\n", (long)H.largest_FB_p);
    err_printf("MPQS: bound for `large primes' = %ld\n", (long)H.lp_bound);
    if (DEBUGLEVEL >= 5)
      err_printf("MPQS: sieve threshold = %u\n", (unsigned int)H.sieve_threshold);
    err_printf("MPQS: computing relations:");
  }

  /* main loop which
   * - computes polynomials and their zeros (SI)
   * - does the sieving
   * - tests candidates of the sieve array */

  /* Let (A, B_i) the current pair of coeffs. If i == 0 a new A is generated */
  H.index_j = (mpqs_uint32_t)-1;  /* increment below will have it start at 0 */

  dbg_target = H.target_rels / 100.;
  DEFEAT = H.target_rels * 1.5;
  hash_init_GEN(&frel, H.target_rels, gequal, 1);
  hash_init_ulong(&lprel,H.target_rels, 1);
  for(;;)
  {
    long tc;
    /* self initialization: compute polynomial and its zeros */
    if (!mpqs_self_init(&H))
    { /* have run out of primes for A; give up */
      if (DEBUGLEVEL >= 2)
        err_printf("MPQS: Ran out of primes for A, giving up.\n");
      return gc_NULL(av);
    }
    mpqs_sieve(&H);
    tc = mpqs_eval_sieve(&H);
    if (DEBUGLEVEL >= 6)
      err_printf("MPQS: found %lu candidate%s\n", tc, (tc==1? "" : "s"));
    if (tc)
    {
      fact = mpqs_eval_cand(&H, tc, &frel, &lprel);
      if (fact)
      { /* factor found during combining */
        if (DEBUGLEVEL >= 4)
        {
          err_printf("\nMPQS: split N whilst combining, time = %ld ms\n",
                     timer_delay(&T));
          err_printf("MPQS: found factor = %Ps\n", fact);
        }
        return gerepileupto(av, fact);
      }
    }
    if (DEBUGLEVEL >= 4 && frel.nb > dbg_target)
    {
      err_printf(" %ld%%", 100*frel.nb/ H.target_rels);
      if (DEBUGLEVEL >= 5) err_printf(" (%ld ms)", timer_delay(&T));
      dbg_target += H.target_rels / 100.;
    }
    if (frel.nb < (ulong)H.target_rels) continue; /* main loop */

    if (DEBUGLEVEL >= 4)
    {
      timer_start(&T);
      err_printf("\nMPQS: starting Gauss over F_2 on %ld distinct relations\n",
                 frel.nb);
    }
    fact = mpqs_solve_linear_system(&H, &frel);
    if (fact)
    { /* solution found */
      if (DEBUGLEVEL >= 4)
      {
        err_printf("\nMPQS: time in Gauss and gcds = %ld ms\n",timer_delay(&T));
        if (typ(fact) == t_INT) err_printf("MPQS: found factor = %Ps\n", fact);
        else
        {
          long j, nf = (lg(fact)-1)/3;
          err_printf("MPQS: found %ld factors =\n", nf);
          for (j = 1; j <= nf; j++)
            err_printf("\t%Ps%s\n", gel(fact,3*j-2), (j < nf)? ",": "");
        }
      }
      return gerepileupto(av, fact);
    }
    if (DEBUGLEVEL >= 4)
    {
      err_printf("\nMPQS: time in Gauss and gcds = %ld ms\n",timer_delay(&T));
      err_printf("MPQS: no factors found.\n");
      if (frel.nb < DEFEAT)
        err_printf("\nMPQS: restarting sieving ...\n");
      else
        err_printf("\nMPQS: giving up.\n");
    }
    if (frel.nb >= DEFEAT) return gc_NULL(av);
    H.target_rels += 10;
  }
}
