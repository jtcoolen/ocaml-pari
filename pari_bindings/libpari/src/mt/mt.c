/* $Id$

Copyright (C) 2013  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
#include "pari.h"
#include "paripriv.h"
#include "mt.h"

static GEN
mtsingle_queue_get(struct mt_state *mt, long *workid, long *pending)
{
  GEN done = mt->pending;
  if (workid) *workid = mt->workid;
  mt->pending = NULL; *pending = 0;
  return done;
}

static int single_is_thread = 0;
static long single_trace_level = 0;

static void
mtsingle_queue_submit(struct mt_state *mt, long workid, GEN work)
{
  int is_thread = single_is_thread;
  single_is_thread = 1;
  mt->pending = work? closure_callgenvec(mt->worker, work): NULL;
  single_is_thread = is_thread;
  mt->workid = workid;
}

static void
mtsingle_queue_end(void) {  }

static GEN
mtsequential_queue_get(struct mt_state *mt, long *workid, long *pending)
{
  GEN done = mt->pending;
  if (workid) *workid = mt->workid;
  mt->pending = NULL; *pending = 0;
  return done;
}

static void
mtsequential_queue_submit(struct mt_state *mt, long workid, GEN work)
{
  mt->pending = work? closure_callgenvec(mt->worker, work): NULL;
  mt->workid = workid;
}

static void
mtsequential_queue_end(void) {  }

int
mtsingle_is_thread(void) { return single_is_thread; }

void
mtsingle_err_recover(long er)
{
  (void) er;
  if (single_is_thread)
  {
    evalstate_set_trace(single_trace_level);
    single_is_thread = 0;
  }
}

void
mtsingle_queue_start(struct pari_mt *pt, GEN worker)
{
  pt->get = mtsingle_queue_get;
  pt->submit = mtsingle_queue_submit;
  pt->end = mtsingle_queue_end;
  pt->mt.worker = worker;
  pt->mt.pending = NULL;
  single_trace_level = evalstate_get_trace();
}

void
mtsequential_queue_start(struct pari_mt *pt, GEN worker)
{
  pt->get = mtsequential_queue_get;
  pt->submit = mtsequential_queue_submit;
  pt->end = mtsequential_queue_end;
  pt->mt.worker = worker;
  pt->mt.pending = NULL;
}

void
mt_queue_end(struct pari_mt *pt) { pt->end(); }

void
mt_queue_submit(struct pari_mt *pt, long workid, GEN work)
{ pt->submit(&pt->mt, workid, work); }

GEN
mt_queue_get(struct pari_mt *pt, long *workid, long *pending)
{ return pt->get(&pt->mt, workid, pending); }

void
mt_queue_start(struct pari_mt *pt, GEN worker)
{ return mt_queue_start_lim(pt, worker, 0); }

void
mtstate_save(struct pari_mtstate *mt)
{
  if (mt_is_parallel())
  {
    mt->is_thread = 0;
    mt->trace_level = 0;
    mt->pending_threads = 1;
  } else
  {
    mt->is_thread = single_is_thread;
    mt->trace_level = single_trace_level;
    mt->pending_threads = 0;
  }
}

void
mtstate_restore(struct pari_mtstate *mt)
{
  if (!mt_is_parallel())
  {
    single_is_thread = mt->is_thread;
    single_trace_level = mt->trace_level;
  }
  if (!mt->pending_threads && mt_is_parallel())
    mt_queue_reset();
}

void
mtstate_reset(void)
{
  single_is_thread = 0;
  single_trace_level = 0;
  if (mt_is_parallel())
    mt_queue_reset();
}
