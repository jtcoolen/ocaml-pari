/* - debug support */

#ifdef MPQS_DEBUG_VERBOSE
#  ifndef MPQS_DEBUG
#  define MPQS_DEBUG
#  endif
#  define PRINT_IF_VERBOSE(x) err_printf(x)
#else
#  define PRINT_IF_VERBOSE(x)
#endif

#ifdef MPQS_DEBUG
#  define MPQS_DEBUGLEVEL 1000  /* infinity */
#else
#  define MPQS_DEBUGLEVEL DEBUGLEVEL
#endif

/* - non-configurable sizing parameters */

/* 'large primes' must be smaller than min(MPQS_LP_BOUND, largest_FB_p) */
#define MPQS_LP_BOUND              12500000 /* works for 32 and 64bit */

/* see mpqs_locate_A_range() for an explanation of the following.  I had
 * some good results with about -log2(0.85) but in the range I was testing,
 * this shifts the primes for A only by one position in the FB.  Don't go
 * over the top with this one... */
#define MPQS_A_FUDGE               0.15 /* ~ -log2(0.9) */

#define MPQS_CANDIDATE_ARRAY_SIZE  2000 /* max. this many cand's per poly */

/* - structures, types, and constants */

/* -- reasonably-sized integers */
#ifdef LONG_IS_64BIT
typedef int  mpqs_int32_t;
typedef unsigned int  mpqs_uint32_t;
#else
typedef long mpqs_int32_t;
typedef ulong mpqs_uint32_t;
#endif

/* -- factor base entries should occupy 32 bytes  (and we'll keep them
 * aligned, for good L1 cache hit rates).  Some of the entries will be
 * abused for e.g. -1 and (factors of) k instead for real factor base
 * primes, and for a sentinel at the end.  This is why __p is a signed
 * field.-- The two start fields depend on the current polynomial and
 * keep changing during sieving, the flags will also change depending on
 * the current A. */
/* Let (z1, z2) be the roots of Q(x) = A x^2 + Bx + C mod p_i; then
 * Q(z1 + p_i Z) == 0 mod p_i and Q(z2 + p_i Z) == 0 mod p_i;
 * start_1, start_2 are the positions where p_i divides Q(x) for the
 * first time, already adjusted for the fact that the sieving array,
 * nominally [-M, M], is represented by a 0-based C array of length
 * 2M + 1.  For the prime factors of A and those of k, the two roots
 * are equal mod p_i. */

#define MPQS_FB_ENTRY_PAD 32

typedef union mpqs_FB_entry {
  char __pad[MPQS_FB_ENTRY_PAD];
  struct {
    /* the prime p, the two arith. prog. mod p, sqrt(kN) mod p */
    mpqs_int32_t __p, __start1, __start2, __sqrt_kN;
    float __flogp; /* log(p) as a 4-byte float */
    unsigned char __val; /* 8-bit approx. scaled log for sieving */
    unsigned char __flags;
  } __entry;
} mpqs_FB_entry_t;

/* --- convenience accessor macros for the preceding: */
#define fbe_p           __entry.__p
#define fbe_flogp       __entry.__flogp
#define fbe_start1      __entry.__start1
#define fbe_start2      __entry.__start2
#define fbe_sqrt_kN     __entry.__sqrt_kN
#define fbe_logval      __entry.__val
#define fbe_flags       __entry.__flags

/* --- flag bits for fbe_flags: */
#define MPQS_FBE_CLEAR       0x0 /* no flags */

/* following used for odd FB primes, and applies to the divisors of A but not
 * those of k.  Must occupy the rightmost bit because we also use it as a
 * shift count after extracting it from the byte. */
#define MPQS_FBE_DIVIDES_A   0x1ul /* and Q(x) mod p only has degree 1 */

/* TODO (tentative): one bit to mark normal FB primes,
 * one to mark the factors of k,
 * one to mark primes used in sieving,
 * later maybe one to mark primes of which we'll be tracking the square,
 * one to mark primes currently in use for A;
 * once we segment the FB, one bit marking the members of the first segment */

/* -- multiplier k and attached quantities */
typedef struct mpqs_multiplier {
  mpqs_uint32_t k;       /* the multiplier (odd, squarefree) */
  mpqs_uint32_t omega_k; /* number (0, 1 or 2) of primes dividing k */
  mpqs_uint32_t kp[2]; /* prime factors of k, if any */
} mpqs_multiplier_t;

#define MPQS_POSSIBLE_MULTIPLIERS  15 /* how many values for k we'll try */
/* following must be in range of the cand_multipliers table below */

static const mpqs_multiplier_t cand_multipliers[] = {
  {  1, 0, { 0,  0}},
  {  3, 1, { 3,  0}},
  {  5, 1, { 5,  0}},
  {  7, 1, { 7,  0}},
  { 11, 1, {11,  0}},
  { 13, 1, {13,  0}},
  { 15, 2, { 3,  5}},
  { 17, 1, {17,  0}},
  { 19, 1, {19,  0}},
  { 21, 2, { 3,  7}},
  { 23, 1, {23,  0}},
  { 29, 1, {29,  0}},
  { 31, 1, {31,  0}},
  { 33, 2, { 3, 11}},
  { 35, 2, { 5,  7}},
  { 37, 1, {37,  0}},
  { 39, 2, { 3, 13}},
  { 41, 1, {41,  0}},
  { 43, 1, {43,  0}},
  { 47, 1, {47,  0}},
  { 51, 2, { 3, 17}},
  { 53, 1, {53,  0}},
  { 55, 2, { 5, 11}},
  { 57, 2, { 3, 19}},
  { 59, 1, {59,  0}},
  { 61, 1, {61,  0}},
  { 65, 2, { 5, 13}},
  { 67, 1, {67,  0}},
  { 69, 2, { 3, 23}},
  { 71, 1, {71,  0}},
  { 73, 1, {73,  0}},
  { 77, 2, { 7, 11}},
  { 79, 1, {79,  0}},
  { 83, 1, {83,  0}},
  { 85, 2, { 5, 17}},
  { 87, 2, { 3, 29}},
  { 89, 1, {89,  0}},
  { 91, 2, { 7, 13}},
  { 93, 2, { 3, 31}},
  { 95, 2, { 5, 19}},
  { 97, 1, {97,  0}}
};

/* -- the array of (Chinese remainder) idempotents which add/subtract up to
 * the middle coefficient B, and for convenience, the FB subscripts of the
 * primes in current use for A.  We keep these together since both arrays
 * are of the same size and are used at the same times. */
typedef struct mqps_per_A_prime {
  GEN _H;          /* summand for B */
  mpqs_int32_t _i; /* subscript into FB */
} mpqs_per_A_prime_t;

/* following cooperate with names of local variables in the self_init fcns.
 * per_A_pr must exist and be an alias for the eponymous handle pointer for
 * all of these, and FB must exist and correspond to the handle FB pointer
 * for all but the first two of them. */
#define MPQS_H(i) (per_A_pr[i]._H)
#define MPQS_I(i) (per_A_pr[i]._i)
#define MPQS_AP(i) (FB[MPQS_I(i)].fbe_p)
#define MPQS_LP(i) (FB[MPQS_I(i)].fbe_flogp)
#define MPQS_SQRT(i) (FB[MPQS_I(i)].fbe_sqrt_kN)
#define MPQS_FLG(i) (FB[MPQS_I(i)].fbe_flags)

/* -- the array of addends / subtrahends for changing polynomials during
 * self-initialization: (1/A) H[i] mod p_j, with i subscripting the inner
 * array in each entry, and j choosing the entry in an outer array.
 * Entries will occupy 64 bytes each no matter what (which imposes at most 17
 * prime factors for A; thus i will range from 0 to at most 15.) This wastes a
 * little memory for smaller N but makes it easier for compilers to generate
 * efficient code. */

/* NOTE: At present, memory locality vis-a-vis accesses to this array is good
 * in the slow (new A) branch of mpqs_self_init(), but poor in the fast
 * (same A, new B) branch, which now loops over the outer array index,
 * reading just one field of each inner array each time through the FB
 * loop.  This doesn't really harm, but will improve one day when we do
 * segmented sieve arrays with the attached segmented FB-range accesses. */
#define MPQS_MAX_OMEGA_A 17
typedef struct mpqs_inv_A_H {
  mpqs_uint32_t _i[MPQS_MAX_OMEGA_A - 1];
} mpqs_inv_A_H_t;

#define MPQS_INV_A_H(i,j) (inv_A_H[j]._i[i])

/* -- global handle to keep track of everything used through one factorization
 * attempt. The order of the fields is determined by keeping most frequently
 * used stuff near the beginning. */
typedef struct mpqs_handle {
  /* pointers */
  unsigned char *sieve_array;/* 0-based, representing [-M,M-1] */
  unsigned char *sieve_array_end; /* points at sieve_array[M-1] */
  mpqs_FB_entry_t *FB;       /* (aligned) FB array itself */
  long *candidates;          /* collects promising sieve subscripts */
  long *relaprimes;          /* prime/exponent pairs in a relation */
  mpqs_inv_A_H_t *inv_A_H;   /* self-init: (aligned) stepping array, and */
  mpqs_per_A_prime_t *per_A_pr; /* FB subscripts of primes in A etc. */

  /* other stuff that's being used all the time */
  mpqs_int32_t M;            /* sieving over |x| <= M */
  mpqs_int32_t size_of_FB;   /* # primes in FB (or dividing k) */
  /* the following three are in non-descending order, and the first two must
   * be adjusted for omega_k at the beginning */
  mpqs_int32_t index0_FB;    /* lowest subscript into FB of a "real" prime
                              * (i.e. other than -1, 2, factors of k) */
  mpqs_int32_t index1_FB;    /* lowest subscript into FB for sieving */
  mpqs_int32_t index2_FB;    /* primes for A are chosen relative to this */
  unsigned char index2_moved;/* true when we're starved for small A's */
  unsigned char sieve_threshold; /* distinguishes candidates in sieve */
  GEN N, kN;                 /* number to be factored, with multiplier */
  GEN A, B;                  /* leading, middle coefficient */
  mpqs_int32_t omega_A;      /* number of primes going into each A */
  mpqs_int32_t no_B;         /* number of B's for each A: 2^(omega_A-1) */
  double l2_target_A;        /* ~log2 of desired typical A */
  /* counters and bit pattern determining and numbering current polynomial: */
  mpqs_uint32_t bin_index;   /* bit pattern for selecting primes for A */
  mpqs_uint32_t index_i;     /* running count of A's */
  mpqs_uint32_t index_j;     /* B's ordinal number in A's cohort */

  /* further sizing parameters: */
  mpqs_int32_t target_rels;  /* target number of full relations */
  mpqs_int32_t largest_FB_p; /* largest prime in the FB */
  mpqs_int32_t pmin_index1;  /* lower bound for primes used for sieving */
  mpqs_int32_t lp_scale;     /* factor by which LPs may exceed FB primes */

  /* subscripts determining where to pick primes for A */
  /* FIXME: lp_bound might have to be mpqs_int64_t ? */
  long lp_bound;             /* cutoff for Large Primes */
  long digit_size_kN;
  const mpqs_multiplier_t *_k;  /* multiplier k and attached quantities */
  double tolerance;          /* controls the tightness of the sieve */
  double dkN;                /* - double prec. approximation of kN */
  double l2sqrtkN;           /* ~log2(sqrt(kN)) */
  double l2M;                /* ~log2(M) (cf. below) */
  /* TODO: need an index2_FB here to remember where to start picking primes */
  /* bookkeeping pointers to containers of aligned memory chunks: */
  void *FB_chunk;            /* (unaligned) chunk containing the FB */
  void *invAH_chunk;         /* (unaligned) chunk for self-init array */
} mpqs_handle_t;

/* -- sizing table entries */

/* For "tolerance", see mpqs_set_sieve_threshold(). The LP scale, for very
 * large kN, prevents us from accumulating vast amounts of LP relations with
 * little chance of hitting any particular large prime a second time and being
 * able to combine a full relation from two LP ones; however, the sieve
 * threshold (determined by the tolerance) already works against very large LPs
 * being produced. The present relations "database" can detect duplicate full
 * relations only during the sort/combine phases, so we must do some sort
 * points even for tiny kN where we do not admit large primes at all.
 * Some constraints imposed by the present implementation:
 * + omega_A should be at least 3, and no more than MPQS_MAX_OMEGA_A
 * + The size of the FB must be large enough compared to omega_A
 *   (about 2*omega_A + 3, but this is always true below) */
/* XXX Changes needed for segmented mode:
 * XXX When using it (kN large enough),
 * XXX - M must become a multiple of the (cache block) segment size
 * XXX   (or to keep things simple: a multiple of 32K)
 * XXX - we need index3_FB to separate (smaller) primes used for normal
 * XXX   sieving from larger ones used with transaction buffers
 * XXX   (and the locate_A_range and attached logic must be changed to
 * XXX   cap index2_FB below index3_FB instead of below size_of_FB)
 */
typedef struct mpqs_parameterset {
  float tolerance;          /* "mesh width" of the sieve */
  mpqs_int32_t lp_scale;    /* factor by which LPs may exceed FB primes */
  mpqs_int32_t M;           /* size of half the sieving interval */
  mpqs_int32_t size_of_FB;  /* #primes to use for FB (including 2) */
  mpqs_int32_t omega_A;     /* #primes to go into each A */
  /* Following is auto-adjusted to account for prime factors of k inserted
   * near the start of the FB. NB never ever sieve on the prime 2,which would
   * just contribute a constant at each sieve point. */
  mpqs_int32_t pmin_index1; /* lower bound for primes used for sieving */
} mpqs_parameterset_t;

/* - the table of sizing parameters itself */

/* indexed by size of kN in decimal digits, subscript 0 corresponding to
 * 9 (or fewer) digits */
static const mpqs_parameterset_t mpqs_parameters[] =
{ /*       tol lp_scl     M   szFB  oA pmx1 */
  {  /*9*/ 0.8,   1,    350,    19,  3,   5},
  { /*10*/ 0.8,   1,    300,    23,  3,   5},
  { /*11*/ 0.8,   1,   1000,    27,  3,   5},
  { /*12*/ 0.8,   1,   1100,    27,  3,   5},
  { /*13*/ 0.8,   1,   1400,    31,  3,   5},
  { /*14*/ 0.8,   1,   2200,    33,  3,   5},
  { /*15*/ 0.8,   1,   2300,    39,  3,   5},
  { /*16*/ 0.8,   1,   2900,    43,  3,   5},
  { /*17*/ 0.8,   1,   3200,    51,  3,   5},
  { /*18*/ 0.8,   1,   2800,    55,  3,   5},
  { /*19*/ 0.8,   1,   3400,    65,  3,   5},
  { /*20*/ 0.8,   1,   3400,    71,  3,   5},
  { /*21*/ 0.8,   1,   5400,    90,  3,   5},
  { /*22*/ 0.8,   1,   5700,    95,  3,   5},
  { /*23*/ 0.8,   1,   5700,   110,  3,   5},
  { /*24*/ 0.8,   1,   6000,   130,  4,   7},
  { /*25*/ 0.8,   1,   6500,   140,  4,   7},
  { /*26*/ 0.9,   1,   9000,   160,  4,   7},
  { /*27*/ 1.12,  1,  10000,   160,  4,   7},
  { /*28*/ 1.17,  1,  13000,   180,  4,  11},
  { /*29*/ 1.22,  1,  14000,   220,  4,  11},
  { /*30*/ 1.30,  1,  13000,   240,  4,  11},
  { /*31*/ 1.33,  1,  11000,   240,  4,  13},
  { /*32*/ 1.36,  1,  14000,   300,  5,  13},
  { /*33*/ 1.40,  1,  15000,   340,  5,  13},
  { /*34*/ 1.43,  1,  15000,   380,  5,  17},
  { /*35*/ 1.48, 30,  15000,   380,  5,  17},
  { /*36*/ 1.53, 45,  16000,   440,  5,  17},
  { /*37*/ 1.60, 60,  15000,   420,  6,  19},
  { /*38*/ 1.66, 70,  15000,   520,  6,  19},
  { /*39*/ 1.69, 80,  16000,   540,  6,  23},
  /* around here, the largest prime in FB becomes comparable to M in size */
  { /*40*/ 1.69, 80,  16000,   600,  6,  23},
  { /*41*/ 1.69, 80,  16000,   700,  6,  23},
  { /*42*/ 1.69, 80,  24000,   900,  6,  29},
  { /*43*/ 1.69, 80,  26000,  1000,  6,  29},
  { /*44*/ 1.69, 80,  18000,  1100,  7,  31},
  { /*45*/ 1.69, 80,  20000,  1200,  7,  31},
  { /*46*/ 1.69, 80,  22000,  1300,  7,  37},
  { /*47*/ 1.69, 80,  24000,  1400,  7,  37},
  { /*48*/ 1.69, 80,  24000,  1600,  7,  37},
  { /*49*/ 1.72, 80,  28000,  1900,  7,  41},
  { /*50*/ 1.75, 80,  36000,  2100,  7,  41},
  { /*51*/ 1.80, 80,  32000,  2100,  7,  43},
  { /*52*/ 1.85, 80,  44000,  2300,  7,  43},
  { /*53*/ 1.90, 80,  44000,  2600,  7,  47},
  { /*54*/ 1.95, 80,  40000,  2700,  7,  47},
  { /*55*/ 1.95, 80,  48000,  3200,  7,  53},
  { /*56*/ 1.95, 80,  56000,  3400,  7,  53},
  { /*57*/ 2.00, 80,  40000,  3000,  8,  53},
  { /*58*/ 2.05, 80,  64000,  3400,  8,  59},
  { /*59*/ 2.10, 80,  64000,  3800,  8,  59},
  { /*60*/ 2.15, 80,  80000,  4300,  8,  61},
  { /*61*/ 2.20, 80,  80000,  4800,  8,  61},
  { /*62*/ 2.25, 80,  80000,  4600,  8,  67},
  { /*63*/ 2.39, 80,  80000,  4800,  8,  67},
  { /*64*/ 2.30, 80,  88000,  5400,  8,  67},
  { /*65*/ 2.31, 80, 120000,  6600,  8,  71},
  { /*66*/ 2.32, 80, 120000,  6800,  8,  71},
  { /*67*/ 2.33, 80, 144000,  7600,  8,  73},
  { /*68*/ 2.34, 80, 144000,  9000,  8,  73},
  { /*69*/ 2.35, 80, 160000,  9500,  8,  79},
  { /*70*/ 2.36, 80, 176000, 10500,  8,  79},
  { /*71*/ 2.37, 80, 240000, 11000,  9,  79},
  { /*72*/ 2.38, 80, 240000, 12500,  9,  83},
  { /*73*/ 2.41, 80, 240000, 13000,  9,  83},
  { /*74*/ 2.46, 80, 256000, 13250,  9,  83},
  { /*75*/ 2.51, 80, 256000, 14500,  9,  89},
  { /*76*/ 2.56, 80, 256000, 15250,  9,  89},
  { /*77*/ 2.58, 80, 320000, 17000,  9,  89},
  { /*78*/ 2.60, 80, 320000, 18000,  9,  89},
  { /*79*/ 2.63, 80, 320000, 19500,  9,  97},
  { /*80*/ 2.65, 80, 448000, 21000,  9,  97},
  { /*81*/ 2.72, 80, 448000, 22000,  9,  97},
  { /*82*/ 2.77, 80, 448000, 24000,  9, 101},
  { /*83*/ 2.82, 80, 480000, 23000, 10, 101},
  { /*84*/ 2.84, 80, 480000, 24000, 10, 103},
  { /*85*/ 2.86, 80, 512000, 28000, 10, 103},
  { /*86*/ 2.88, 80, 448000, 29000, 10, 107},
  /* architectures with 1MBy L2 cache will become noticeably slower here
   * as 2*M exceeds that mark - to be addressed in a future version by
   * segmenting the sieve interval */
  { /*87*/ 2.90, 80, 512000, 32000, 10, 107},
  { /*88*/ 2.91, 80, 512000, 35000, 10, 109},
  { /*89*/ 2.92, 80, 512000, 38000, 10, 109},
  { /*90*/ 2.93, 80, 512000, 40000, 10, 113},
  { /*91*/ 2.94, 80, 770000, 32200, 10, 113},
  /* entries below due to Thomas Denny, never tested */
  { /*92*/ 3.6, 90, 2000000, 35000,  9, 113},
  { /*93*/ 3.7, 90, 2000000, 37000,  9, 113},
  { /*94*/ 3.7, 90, 2000000, 39500,  9, 127},
  { /*95*/ 3.7, 90, 2500000, 41500,  9, 127},
  { /*96*/ 3.8, 90, 2500000, 45000, 10, 127},
  { /*97*/ 3.8, 90, 2500000, 47500, 10, 131},
  { /*98*/ 3.7, 90, 3000000, 51000, 10, 131},
  { /*99*/ 3.8, 90, 3000000, 53000, 10, 133},
  {/*100*/ 3.8, 90,  875000, 50000, 10, 133},
  {/*101*/ 3.8, 90, 3500000, 54000, 10, 139},
  {/*102*/ 3.8, 90, 3500000, 57000, 10, 139},
  {/*103*/ 3.9, 90, 4000000, 61000, 10, 139},
  {/*104*/ 3.9, 90, 4000000, 66000, 10, 149},
  {/*105*/ 3.9, 90, 4000000, 70000, 10, 149},
  {/*106*/ 3.9, 90, 4000000, 75000, 10, 151},
  {/*107*/ 3.9, 90, 4000000, 80000, 10, 151},
};

#define MPQS_MAX_DIGIT_SIZE_KN 107
