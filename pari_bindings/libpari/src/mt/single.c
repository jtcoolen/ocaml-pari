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
#include "pari.h"
#include "paripriv.h"
#include "mt.h"

void mt_sigint_block(void) { }
void mt_sigint_unblock(void) { }
void mt_err_recover(long er) { mtsingle_err_recover(er); }
void mt_break_recover(void) { mtsingle_err_recover(0); }
void pari_mt_close(void) { }
void mt_queue_reset(void) { }
void mt_broadcast(GEN code) {(void) code;}
void mt_thread_init(void) { }

void
mt_sigint(void) {}

int
mt_is_parallel(void)
{
  return 0;
}

int
mt_is_thread(void)
{
  return mtsingle_is_thread();
}

long
mt_nbthreads(void)
{
  return 1;
}

void
mt_export_add(const char *str, GEN val)
{
  if (mtsingle_is_thread())
    pari_err(e_MISC,"export not allowed during parallel sections");
  export_add(str, val);
}

void
mt_export_del(const char *str)
{
  if (mtsingle_is_thread())
    pari_err(e_MISC,"unexport not allowed during parallel sections");
  export_del(str);
}

void
pari_mt_init(void)
{
  pari_mt_nbthreads = 1;
}

void
mt_queue_start_lim(struct pari_mt *pt, GEN worker, long lim)
{
  (void) lim;
  mtsingle_queue_start(pt, worker);
}
