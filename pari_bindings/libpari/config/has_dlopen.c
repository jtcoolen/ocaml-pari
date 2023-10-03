#include <stdio.h>
#include <dlfcn.h>

int main()
{
  void *handle = dlopen(NULL,RTLD_LAZY);
  if (handle)
    printf("%p\n",dlsym(handle, "fun"));
  return 0;
}
