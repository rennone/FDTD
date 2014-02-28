#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
//#include "shelf.h"
//#include "nonshelf.h"


static double (*epsMethod)(double, double, int, int);

static void noModel(void)
{
  //no material
  epsMethod = noModel_EPS();
}

static void circleModel(void)
{
  //cylinder material whitch radius = lambda, origin = center of field
  epsMethod = circleModel_EPS(N_PX*0.5, N_PY*0.5, field_getLambda());
}

void setModel(enum MODEL model)
{
  switch(model){
  case NO_MODEL:
    noModel();
    break;
  case MIE_CYLINDER:
    circleModel();
    break;
  case SHELF :
  case NONSHELF:
  case LAYER:
    break;
  }
}

double models_eps(double x, double y, enum MODE mode){
  double epsilon = EPSILON_0_S;
  switch(mode){
  case D_X :
    epsilon = (*epsMethod)(x, y, 1, 0);
  case D_Y :
    epsilon = (*epsMethod)(x, y, 0, 1);
  case D_XY :
    epsilon = (*epsMethod)(x, y, 1, 1);
  }
  return epsilon;
}
