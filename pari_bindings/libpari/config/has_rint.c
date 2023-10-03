#include <math.h>
double (*f)(double) = rint;
int main(){ return f != rint; }
