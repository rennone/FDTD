#ifndef _FDTD_TE_H
#define _FDTD_TE_H
#include <complex.h>

extern void(*fdtdTE_getUpdate(void))(void);
extern void(*fdtdTE_getFinish(void))(void);
extern void(*fdtdTE_getInit(void))(void);

extern double complex* fdtdTE_getEx(void);
extern double complex* fdtdTE_getEy(void);
extern double complex* fdtdTE_getHzx(void);
extern double complex* fdtdTE_getHzy(void);
extern double complex* fdtdTE_getHz(void);

extern double* fdtdTE_getEps();
#endif
