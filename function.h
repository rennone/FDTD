#ifndef _FUNCTION_H
#define _FUNCTION_H
#include <complex.h>

static inline double max(double a, double b){
  return a > b ? a : b;
}

static inline double min(double a, double b){
  return a > b ? b : a;
}

//norm of complex
static inline double cnorm(double complex c){
  double re = creal(c);
  double im = cimag(c);
  return re*re + im*im;
}

static inline double complex cbilinear(double complex *p, double x, double y, int width, int height)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  int index = i*height + j;
  return p[index]*(1.0-dx)*(1.0-dy)
       + p[index+height]*dx*(1.0-dy)
       + p[index+1]*(1.0-dx)*dy
       + p[index+height+1]*dx*dy;
}

static inline double dbilinear(double *p, double x, double y, int width, int height)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  int index = i*height + j;
  return p[index]*(1.0-dx)*(1.0-dy)
       + p[index+height]*dx*(1.0-dy)
       + p[index+1]*(1.0-dx)*dy
       + p[index+height+1]*dx*dy;
}

#endif
