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

BEGINEXTERN

void mtsingle_queue_start(struct pari_mt *pt, GEN worker);
void mtsequential_queue_start(struct pari_mt *pt, GEN worker);
int  mtsingle_is_thread(void);
void mtsingle_err_recover(long er);
void mt_queue_reset(void);
int  mt_is_parallel(void);

ENDEXTERN
