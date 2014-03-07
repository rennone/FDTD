#ifndef _SIMULATOR_H
#define _SIMULATOR_H
#include<complex.h>

#include "models.h"

enum SOLVER{
  TE_2D,
  TM_2D,
  TM_UPML_2D,
  TE_UPML_2D,
};

extern void simulator_init(int width, int height , double h_u, int pml, double lambda, int waveAngle, int step, enum MODEL model, enum SOLVER solver);
extern void simulator_calc(void);
extern double complex* simulator_getDrawingData();
extern double * simulator_getEps();
#endif
