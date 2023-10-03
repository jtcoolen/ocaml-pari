/* Copyright (C) 2013  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */
#include <mpi.h>
#include "pari.h"
#include "paripriv.h"
#include "../src/language/anal.h"
#include "mt.h"

#define DEBUGLEVEL DEBUGLEVEL_mt

static THREAD int pari_MPI_size, pari_MPI_rank;
static THREAD long nbreq = 0;

enum PMPI_cmd { PMPI_close, PMPI_worker, PMPI_work,
                PMPI_eval, PMPI_exportadd, PMPI_exportdel};

struct mt_mstate
{
  long n;
  int source;
  long nbint;
  long *workid;
};

static struct mt_mstate pari_mt_data;
static struct mt_mstate *pari_mt;

static void
send_long(long a, long dest)
{
  BLOCK_SIGINT_START
  MPI_Send(&a, 1, MPI_LONG, dest, 0, MPI_COMM_WORLD);
  BLOCK_SIGINT_END
}

static void
send_bcast_long(long a, MPI_Comm comm)
{
  BLOCK_SIGINT_START
  MPI_Bcast(&a, 1, MPI_LONG, 0, comm);
  BLOCK_SIGINT_END
}

static long
recv_bcast_long(MPI_Comm comm)
{
  long a;
  BLOCK_SIGINT_START
  MPI_Bcast(&a, 1, MPI_LONG, 0, comm);
  BLOCK_SIGINT_END
  return a;
}

static void
send_request(enum PMPI_cmd ecmd, long dest)
{
  send_long((long)ecmd, dest);
}

static void
send_request_all(enum PMPI_cmd ecmd, long n)
{
  long i;
  for (i=1; i<=n; i++)
    send_long((long)ecmd, i);
}

static void
send_GEN(GEN elt, int dest)
{
  pari_sp av = avma;
  GEN reloc = copybin_unlink(elt);
  GENbin *buf = copy_bin_canon(mkvec2(elt,reloc));
  long size = sizeof(GENbin) + buf->len*sizeof(ulong);
  {
    BLOCK_SIGINT_START
    MPI_Send(buf, size/sizeof(long), MPI_LONG, dest, 0, MPI_COMM_WORLD);
    BLOCK_SIGINT_END
  }
  pari_free(buf); set_avma(av);
}

static void
send_bcast_GEN(GEN elt, MPI_Comm comm)
{
  pari_sp av = avma;
  GEN reloc = copybin_unlink(elt);
  GENbin *buf = copy_bin_canon(mkvec2(elt,reloc));
  long size = sizeof(GENbin) + buf->len*sizeof(ulong);
  {
    send_bcast_long(size, comm);
    BLOCK_SIGINT_START
    MPI_Bcast(buf, size/sizeof(long), MPI_LONG, 0, comm);
    BLOCK_SIGINT_END
  }
  pari_free(buf); set_avma(av);
}

static void
send_bcast_vlong(long *a, long n, MPI_Comm comm)
{
  BLOCK_SIGINT_START
  MPI_Bcast(a, n, MPI_LONG, 0, comm);
  BLOCK_SIGINT_END
}

static void
send_request_GEN(enum PMPI_cmd ecmd, GEN elt, int dest)
{
  send_request(ecmd, dest);
  send_GEN(elt, dest);
}

static long
recvfrom_long(int src)
{
  long a;
  BLOCK_SIGINT_START
  MPI_Recv(&a, 1, MPI_LONG, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  BLOCK_SIGINT_END
  return a;
}

static enum PMPI_cmd
recvfrom_request(int src)
{
  return (enum PMPI_cmd) recvfrom_long(src);
}

static GENbin *
recvstatus_buf(int source, MPI_Status *status)
{
  int size;
  GENbin *buf;
  BLOCK_SIGINT_START

  MPI_Get_count(status, MPI_LONG, &size);
  buf = (GENbin *)pari_malloc(size*sizeof(long));
  MPI_Recv(buf, size, MPI_LONG, source, 0/* tag */,
          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  BLOCK_SIGINT_END
  return buf;
}

static GEN
recvstatus_GEN(int source, MPI_Status *status)
{
  GEN res;
  GENbin *buf = recvstatus_buf(source, status);
  buf->rebase = &shiftaddress_canon;
  res = bin_copy(buf);
  bincopy_relink(gel(res,1),gel(res,2));
  return gel(res,1);
}

static void
recvstatus_void(int source, MPI_Status *status)
{
  GENbin *buf = recvstatus_buf(source, status);
  free(buf);
}

static GEN
recvfrom_GEN(int src)
{
  MPI_Status status;
  BLOCK_SIGINT_START
  MPI_Probe(src, 0, MPI_COMM_WORLD, &status);
  BLOCK_SIGINT_END
  return recvstatus_GEN(src, &status);
}

static GEN
recv_bcast_GEN(MPI_Comm comm)
{
  GEN res;
  GENbin *buf;
  long size;

  size = recv_bcast_long(comm);
  buf = (GENbin *)pari_malloc(size);

  BLOCK_SIGINT_START
  MPI_Bcast(buf, size/sizeof(long), MPI_LONG, 0, comm);
  BLOCK_SIGINT_END

  buf->rebase = &shiftaddress_canon;
  res = bin_copy(buf);
  bincopy_relink(gel(res,1),gel(res,2));
  return gel(res,1);
}

static void
recv_bcast_vlong(long *a, long n, MPI_Comm comm)
{
  BLOCK_SIGINT_START
  MPI_Bcast(a, n, MPI_LONG, 0, comm);
  BLOCK_SIGINT_END
}

static GEN
recvany_GEN(int *source)
{
  MPI_Status status;
  BLOCK_SIGINT_START
  MPI_Probe(MPI_ANY_SOURCE, 0 /* tag */, MPI_COMM_WORLD, &status);
  *source = status.MPI_SOURCE;
  BLOCK_SIGINT_END
  return recvstatus_GEN(*source, &status);
}

static void
recvany_void(int *source)
{
  MPI_Status status;
  BLOCK_SIGINT_START
  MPI_Probe(MPI_ANY_SOURCE, 0 /* tag */, MPI_COMM_WORLD, &status);
  *source = status.MPI_SOURCE;
  BLOCK_SIGINT_END
  recvstatus_void(*source, &status);
}

static jmp_buf child_env;

static void
pari_MPI_child(void)
{
  pari_sp av = avma;
  ulong rsize = 0, vsize = 0;
  GEN worker = NULL, work, done;
  struct gp_context rec;
  pari_mt_nbthreads = 1;
  gp_context_save(&rec);
  if (setjmp(child_env))
  {
    send_GEN(pari_err_last(), 0);
    gp_context_restore(&rec);
  }
  while (1)
    switch (recvfrom_request(0))
    {
    case PMPI_worker:
      {
        MPI_Comm comm;
        long n = recv_bcast_long(MPI_COMM_WORLD), status = pari_MPI_rank <= n;
        MPI_Comm_split(MPI_COMM_WORLD, status, pari_MPI_rank, &comm);
        if (status==0)
        {
          MPI_Comm_free(&comm);
          break;
        }
        rsize = recv_bcast_long(comm);
        vsize = recv_bcast_long(comm);
        precreal = recv_bcast_long(comm);
        recv_bcast_vlong(varpriority-1,MAXVARN+2,comm);
        {
          pari_sp ltop = avma;
          GEN tab = recv_bcast_GEN(comm);
          if (!gequal(tab, primetab))
          {
            long i, l = lg(tab);
            GEN old = primetab, t = cgetg_block(l, t_VEC);
            for (i = 1; i < l; i++) gel(t,i) = gclone(gel(tab,i));
            primetab = t;
            gunclone_deep(old);
          }
          set_avma(ltop);
        }
        paristack_setsize(rsize, vsize);
        worker = recv_bcast_GEN(comm);
        MPI_Comm_free(&comm);
        gp_context_save(&rec);
        av = avma;
        break;
      }
    case PMPI_work:
      work = recvfrom_GEN(0);
      done = closure_callgenvec(worker, work);
      send_GEN(done, 0);
      set_avma(av);
      break;
    case PMPI_eval:
      (void) closure_evalgen(recv_bcast_GEN(MPI_COMM_WORLD));
      set_avma(av);
      break;
    case PMPI_exportadd:
    {
      GEN str = recv_bcast_GEN(MPI_COMM_WORLD);
      GEN val = recv_bcast_GEN(MPI_COMM_WORLD);
      entree *ep = fetch_entry(GSTR(str));
      export_add(ep->name, val);
      set_avma(av);
      break;
    }
    case PMPI_exportdel:
    {
      GEN str = recv_bcast_GEN(MPI_COMM_WORLD);
      entree *ep = fetch_entry(GSTR(str));
      export_del(ep->name);
      set_avma(av);
      break;
    }
    case PMPI_close:
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(0);
      break;
    }
}

void
mt_err_recover(long er)
{
  if (pari_MPI_rank) longjmp(child_env,er);
  else mtsingle_err_recover(er);
}

void
mt_break_recover(void)
{
  if (!pari_MPI_rank) mtsingle_err_recover(0);
}

void mt_sigint_block(void) { }
void mt_sigint_unblock(void) { }
void mt_sigint(void) {}
void mt_thread_init(void) { }

int
mt_is_parallel(void)
{
  return !!pari_mt;
}

int
mt_is_thread(void)
{
  return pari_MPI_rank ? 1 : mtsingle_is_thread();
}

long
mt_nbthreads(void)
{
  return pari_mt || pari_MPI_rank || pari_MPI_size <= 2 ? 1: pari_mt_nbthreads;
}

void
mt_export_add(const char *str, GEN val)
{
  pari_sp av = avma;
  long n = pari_MPI_size-1;
  GEN s;
  if (pari_mt || pari_MPI_rank)
    pari_err(e_MISC,"export not allowed during parallel sections");
  export_add(str, val);
  s = strtoGENstr(str);
  send_request_all(PMPI_exportadd, n);
  send_bcast_GEN(s, MPI_COMM_WORLD);
  send_bcast_GEN(val, MPI_COMM_WORLD);
  set_avma(av);
}

void
mt_export_del(const char *str)
{
  pari_sp av = avma;
  long n = pari_MPI_size-1;
  if (pari_MPI_rank)
    pari_err(e_MISC,"unexport not allowed during parallel sections");
  export_del(str);
  send_request_all(PMPI_exportdel, n);
  send_bcast_GEN(strtoGENstr(str), MPI_COMM_WORLD);
  set_avma(av);
}

void
mt_broadcast(GEN code)
{
  if (!pari_MPI_rank && !pari_mt)
  {
    send_request_all(PMPI_eval, pari_MPI_size-1);
    send_bcast_GEN(code, MPI_COMM_WORLD);
  }
}

void
pari_mt_init(void)
{
  int res = MPI_Init(0, NULL);
  if (res == MPI_SUCCESS)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &pari_MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &pari_MPI_rank);
    if (pari_MPI_rank) pari_MPI_child();
#ifdef _IOFBF
  /* HACK: most MPI implementation does not handle stdin well.
  stdinsize is sufficient for the largest test file to fit */
  setvbuf(stdin,pari_malloc(128*1024),_IOFBF,128*1024);
#endif
    if (!pari_mt_nbthreads)
      pari_mt_nbthreads = maxss(1, pari_MPI_size-1);
  }
  else
  {
    pari_MPI_size = 0;
    pari_MPI_rank = 0;
    pari_mt_nbthreads = 1;
  }
  pari_mt = NULL;
}

void
pari_mt_close(void)
{
  long i;
  if (!pari_MPI_rank)
    for (i = 1; i < pari_MPI_size; i++)
      send_request(PMPI_close, i);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}

static GEN
mtmpi_queue_get(struct mt_state *junk, long *workid, long *pending)
{
  struct mt_mstate *mt = pari_mt;
  GEN done;
  (void) junk;
  if (mt->nbint<=mt->n) { mt->source=mt->nbint; *pending = nbreq; return NULL; }
  done = recvany_GEN(&mt->source);
  nbreq--; *pending = nbreq;
  if (workid) *workid = mt->workid[mt->source];
  if (typ(done) == t_ERROR)
  {
    if (err_get_num(done)==e_STACK)
      pari_err(e_STACKTHREAD);
    else
      pari_err(0,done);
  }
  return done;
}

static void
mtmpi_queue_submit(struct mt_state *junk, long workid, GEN work)
{
  struct mt_mstate *mt = pari_mt;
  (void) junk;
  if (!work) { mt->nbint=mt->n+1; return; }
  if (mt->nbint<=mt->n) mt->nbint++;
  nbreq++;
  mt->workid[mt->source] = workid;
  send_request_GEN(PMPI_work, work, mt->source);
}

void
mt_queue_reset(void)
{
  struct mt_mstate *mt = pari_mt;
  if (DEBUGLEVEL>0 && nbreq)
    pari_warn(warner,"%ld discarded threads (MPI)",nbreq);
  for(  ;nbreq>0;  nbreq--) recvany_void(&mt->source);
  pari_free(mt->workid);
  pari_mt = NULL;
}

void
mt_queue_start_lim(struct pari_mt *pt, GEN worker, long lim)
{
  if (lim==0) lim = pari_mt_nbthreads;
  else        lim = minss(pari_mt_nbthreads, lim);
  if (pari_MPI_rank)
    mtsequential_queue_start(pt, worker);
  else if (pari_mt || pari_MPI_size <= 2 || lim <= 1)
    mtsingle_queue_start(pt, worker);
  else
  {
    MPI_Comm comm;
    struct mt_mstate *mt = &pari_mt_data;
    long n = minss(lim, pari_MPI_size-1);
    long mtparisize = GP_DATA->threadsize? GP_DATA->threadsize: pari_mainstack->rsize;
    long mtparisizemax = GP_DATA->threadsizemax;
    pari_mt = mt;
    mt->workid = (long*) pari_malloc(sizeof(long)*(n+1));
    send_request_all(PMPI_worker, pari_MPI_size-1);
    send_bcast_long(n, MPI_COMM_WORLD);
    MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &comm);
    send_bcast_long(mtparisize, comm);
    send_bcast_long(mtparisizemax, comm);
    send_bcast_long(get_localbitprec(), comm);
    send_bcast_vlong(varpriority-1,MAXVARN+2, comm);
    send_bcast_GEN(primetab, comm);
    send_bcast_GEN(worker, comm);
    MPI_Comm_free(&comm);
    mt->n = n;
    mt->nbint = 1;
    mt->source = 1;
    pt->get=&mtmpi_queue_get;
    pt->submit=&mtmpi_queue_submit;
    pt->end=&mt_queue_reset;
  }
}
