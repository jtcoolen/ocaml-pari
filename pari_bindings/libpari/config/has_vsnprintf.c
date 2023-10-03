#include <stdio.h>
#include <stdarg.h>

int main() { return 0; }
int f(int i,...) { char s[2]; va_list ap; va_start(ap,i);
  vsnprintf(s,2," ",ap); return 0;
}


