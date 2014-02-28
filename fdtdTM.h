#ifndef _FDTD_TM_H
#define _FDTD_TM_H
#include <complex.h>

extern void(*fdtdTM_getUpdate(void))(void);
extern void(*fdtdTM_getFinish(void))(void);
extern void(* fdtdTM_getInit(void))(void);

extern double complex* fdtdTM_getEzx(void);
extern double complex* fdtdTM_getEzy(void);
extern double complex* fdtdTM_getEz(void);
extern double complex* fdtdTM_getHy(void);
extern double complex* fdtdTM_getHx(void);

extern double* fdtdTM_getEps();
#endif
