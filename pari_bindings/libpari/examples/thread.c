#include <pari/pari.h> /* Include PARI headers */

#include <pthread.h>   /* Include POSIX threads headers */

void *
mydet(void *arg)
{
  GEN F, M;
  /* Set up thread stack and get thread parameter */
  M = pari_thread_start((struct pari_thread*) arg);
  F = QM_det(M);
  /* Free memory used by the thread */
  pari_thread_close();
  return (void*)F;
}

void *
myfactor(void *arg)  /* same principle */
{
  GEN F, N;
  N = pari_thread_start((struct pari_thread*) arg);
  F = factor(N);
  pari_thread_close();
  return (void*)F;
}

int
main(void)
{
  long prec = DEFAULTPREC;
  GEN M1,M2, N1,N2, F1,F2, D1,D2;
  pthread_t th1, th2, th3, th4; /* POSIX-thread variables */
  struct pari_thread pth1, pth2, pth3, pth4; /* pari thread variables */

  /* Initialise the main PARI stack and global objects (gen_0, etc.) */
  pari_init(32000000,500000);
  /* Compute in the main PARI stack */
  N1 = addis(int2n(256), 1); /* 2^256 + 1 */
  N2 = subis(int2n(193), 1); /* 2^193 - 1 */
  M1 = mathilbert(149);
  M2 = mathilbert(150);
  /* Allocate pari thread structures */
  pari_thread_alloc(&pth1,8000000,N1);
  pari_thread_alloc(&pth2,8000000,N2);
  pari_thread_alloc(&pth3,32000000,M1);
  pari_thread_alloc(&pth4,32000000,M2);
  /* pthread_create() and pthread_join() are standard POSIX-thread
   * functions to start and get the result of threads. */
  pthread_create(&th1,NULL, &myfactor, (void*)&pth1);
  pthread_create(&th2,NULL, &myfactor, (void*)&pth2);
  pthread_create(&th3,NULL, &mydet,    (void*)&pth3);
  pthread_create(&th4,NULL, &mydet,    (void*)&pth4); /* Start 4 threads */
  pthread_join(th1,(void*)&F1);
  pthread_join(th2,(void*)&F2);
  pthread_join(th3,(void*)&D1);
  pthread_join(th4,(void*)&D2); /* Wait for termination, get the results */
  pari_printf("F1=%Ps\nF2=%Ps\nlog(D1)=%Ps\nlog(D2)=%Ps\n",
              F1,F2, glog(D1,prec),glog(D2,prec));
  pari_thread_free(&pth1);
  pari_thread_free(&pth2);
  pari_thread_free(&pth3);
  pari_thread_free(&pth4); /* clean up */
  return 0;
}
