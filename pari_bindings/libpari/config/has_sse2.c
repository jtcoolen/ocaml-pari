#include <emmintrin.h>
typedef __v2di bit_array;
#define AND(a,b) ((a)&(b))
#define EXT0(a) ((unsigned long)__builtin_ia32_vec_ext_v2di((__v2di)(a), 0))
#define EXT1(a) ((unsigned long)__builtin_ia32_vec_ext_v2di((__v2di)(a), 1))
#define TEST(a) (EXT0(a) || EXT1(a))
#define RBA(a) ((bit_array){((long) a), ((long) a)})

int main(void)
{
  bit_array x = RBA(1L), y = RBA(3L);
  unsigned long t = TEST(AND(x,y));
  (void) t;
  return 0;
}
