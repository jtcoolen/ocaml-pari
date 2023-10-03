#include <immintrin.h>
 /* Use AVX 256 bit registers for the bit arrays */
typedef unsigned long bit_array __attribute__ ((vector_size (32)));
#define EXT0(a) ((unsigned long)a[0])
#define EXT(a,i) ((unsigned long)a[i])
#ifdef __AVX2__
#define TEST(a) ( _mm256_movemask_epi8(_mm256_cmpeq_epi8((__m256i)(a), (__m256i)RBA(0))) != 0xffffffffU )
#elif defined(__AVX__)
#define TEST(a) ( !_mm256_testz_si256((__m256i)(a), (__m256i)(a)) )
#else
#define TEST(a) (EXT(a,0) || EXT(a,1) || EXT(a,2) || EXT(a,3))
#endif
#define RBA(a) ((bit_array){((unsigned long) a), ((unsigned long) a), ((unsigned long) a), ((unsigned long) a)})
int main(void)
{
  bit_array x = RBA(1L), y = RBA(3L);
  unsigned long t = TEST(x&y);
  (void) t;
  return 0;
}
