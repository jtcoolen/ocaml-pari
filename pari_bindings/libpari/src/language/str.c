/* Copyright (C) 2018  The PARI group.

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

/********************************************************************/
/**                                                                **/
/**                       CHARACTER STRINGS                        **/
/**                                                                **/
/********************************************************************/

/* Utillity functions */
char *
stack_strdup(const char *s)
{
  long n = strlen(s)+1;
  char *t = stack_malloc(n);
  memcpy(t,s,n); return t;
}
char *
stack_strcat(const char *s, const char *t)
{
  long ls = strlen(s), lt = strlen(t);
  long n = ls + lt + 1;
  char *u = stack_malloc(n);
  memcpy(u,     s, ls);
  memcpy(u + ls,t, lt+1); return u;
}

char *
pari_strdup(const char *s)
{
  long n = strlen(s)+1;
  char *t = (char*)pari_malloc(n);
  memcpy(t,s,n); return t;
}
char *
pari_strndup(const char *s, long n)
{
  char *t = (char*)pari_malloc(n+1);
  memcpy(t,s,n); t[n] = 0; return t;
}

/* return the first n0 chars of s as a GEN [s may not be 0-terminated] */
GEN
strntoGENstr(const char *s, long n0)
{
  long n = nchar2nlong(n0+1); /* +1 for trailing 0 */
  GEN x = cgetg(n+1, t_STR);
  char *t = GSTR(x);
  x[n] = 0; /* avoid uninitialized memory */
  strncpy(t, s, n0); t[n0] = 0; return x;
}

/* strntoGENstr would trigger gcc-8 stringop-truncation warning */
GEN
strtoGENstr(const char *s)
{
  long n0 = strlen(s) + 1, n = nchar2nlong(n0);
  GEN x = cgetg(n+1, t_STR);
  char *t = GSTR(x);
  x[n] = 0; strncpy(t, s, n0); return x;
}

GEN
chartoGENstr(char c)
{
  GEN x = cgetg(2, t_STR);
  char *t = GSTR(x);
  t[0] = c; t[1] = 0; return x;
}

const char *
type_name(long t)
{
  const char *s;
  switch(t)
  {
    case t_INT    : s="t_INT";     break;
    case t_REAL   : s="t_REAL";    break;
    case t_INTMOD : s="t_INTMOD";  break;
    case t_FRAC   : s="t_FRAC";    break;
    case t_FFELT  : s="t_FFELT";   break;
    case t_COMPLEX: s="t_COMPLEX"; break;
    case t_PADIC  : s="t_PADIC";   break;
    case t_QUAD   : s="t_QUAD";    break;
    case t_POLMOD : s="t_POLMOD";  break;
    case t_POL    : s="t_POL";     break;
    case t_SER    : s="t_SER";     break;
    case t_RFRAC  : s="t_RFRAC";   break;
    case t_QFB    : s="t_QFB";     break;
    case t_VEC    : s="t_VEC";     break;
    case t_COL    : s="t_COL";     break;
    case t_MAT    : s="t_MAT";     break;
    case t_LIST   : s="t_LIST";    break;
    case t_STR    : s="t_STR";     break;
    case t_VECSMALL:s="t_VECSMALL";break;
    case t_CLOSURE: s="t_CLOSURE"; break;
    case t_ERROR:   s="t_ERROR";   break;
    case t_INFINITY:s="t_INFINITY";break;
    default: pari_err_BUG("type"); s = NULL; /* LCOV_EXCL_LINE */
  }
  return s;
}

GEN
type0(GEN x)
{
  const char *s = type_name(typ(x));
  return strtoGENstr(s);
}

static char
ltoc(long n) {
  if (n <= 0 || n > 255)
    pari_err(e_MISC, "out of range in integer -> character conversion (%ld)", n);
  return (char)n;
}
static char
itoc(GEN x) { return ltoc(gtos(x)); }

GEN
pari_strchr(GEN g)
{
  long i, l, len, t = typ(g);
  char *s;
  GEN x;
  if (is_vec_t(t)) {
    l = lg(g); len = nchar2nlong(l);
    x = cgetg(len+1, t_STR); s = GSTR(x);
    for (i=1; i<l; i++) *s++ = itoc(gel(g,i));
  }
  else if (t == t_VECSMALL)
  {
    l = lg(g); len = nchar2nlong(l);
    x = cgetg(len+1, t_STR); s = GSTR(x);
    for (i=1; i<l; i++) *s++ = ltoc(g[i]);
  }
  else
    return chartoGENstr(itoc(g));
  *s = 0; return x;
}

GEN
strsplit(GEN x, GEN p)
{
  long i0, i, iv, ls, lt;
  char *s, *t;
  GEN v;
  if (typ(x) != t_STR) pari_err_TYPE("strsplit",x);
  s = GSTR(x); ls = strlen(s);
  if (!p) lt = 0;
  else
  {
    if (typ(p) != t_STR) pari_err_TYPE("strsplit",p);
    t = GSTR(p); lt = strlen(t);
  }
  if (!lt) /* empty separator: split by char */
  {
    v = cgetg(ls+1, t_VEC);
    for (i = 1; i <= ls; i++) gel(v,i) = chartoGENstr(s[i-1]);
    return v;
  }
  v = cgetg(ls + 2, t_VEC); iv = 1;
  for (i = i0 = 0; i < ls; i++)
    while (!strncmp(s + i, t, lt))
    {
      gel(v, iv++) = strntoGENstr(s + i0, i - i0);
      i += lt; i0 = i;
    }
  gel(v, iv++) = strntoGENstr(s + i0, i - i0);
  stackdummy((pari_sp)(v + iv), (pari_sp)(v + ls + 2));
  setlg(v, iv); return v;
}

GEN
strjoin(GEN v, GEN p)
{
  pari_sp av = avma;
  long i, l;
  GEN w;
  if (!is_vec_t(typ(v))) pari_err_TYPE("strjoin",v);
  if (p && typ(p) != t_STR) pari_err_TYPE("strjoin",p);
  l = lg(v);
  if (l == 1) return strtoGENstr("");
  if (l == 2)
  {
    char *s = GENtostr_unquoted(gel(v,1));
    return gerepileuptoleaf(av, strtoGENstr(s));
  }
  if (!p) p = strtoGENstr("");
  w = cgetg(2*l - 2, t_VEC);
  gel(w, 1) = gel(v, 1);
  for (i = 2; i < l; i++)
  {
    gel(w, 2*i-2) = p;
    gel(w, 2*i-1) = gel(v, i);
  }
  return gerepileuptoleaf(av, shallowconcat1(w));
}
