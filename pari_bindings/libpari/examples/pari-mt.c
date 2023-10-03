#include <pari/pari.h>

GEN
Cworker(GEN d, long kind) { return kind? det(d): Z_factor(d); }

int
main(void)
{
  long i, taskid, pending;
  GEN M,N1,N2, in,out, done;
  struct pari_mt pt;
  entree ep = {"_worker",0,(void*)Cworker,20,"GL",""};
  /* initialize PARI, postponing parallelism initialization */
  pari_init_opts(8000000,500000, INIT_JMPm|INIT_SIGm|INIT_DFTm|INIT_noIMTm);
  pari_add_function(&ep); /* add Cworker function to gp */
  pari_mt_init(); /* ... THEN initialize parallelism */
  /* Create inputs and room for output in main PARI stack */
  N1 = addis(int2n(256), 1); /* 2^256 + 1 */
  N2 = subis(int2n(193), 1); /* 2^193 - 1 */
  M = mathilbert(80);
  in  = mkvec3(mkvec2(N1,gen_0), mkvec2(N2,gen_0), mkvec2(M,gen_1));
  out = cgetg(4,t_VEC);
  /* Initialize parallel evaluation of Cworker */
  mt_queue_start(&pt, strtofunction("_worker"));
  for (i = 1; i <= 3 || pending; i++)
  { /* submit job (in) and get result (out) */
    mt_queue_submit(&pt, i, i<=3? gel(in,i): NULL);
    done = mt_queue_get(&pt, &taskid, &pending);
    if (done) gel(out,taskid) = done;
  }
  mt_queue_end(&pt); /* end parallelism */
  output(out); pari_close(); return 0;
}
