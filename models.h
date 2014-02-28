#ifndef _MODELS_H
#define _MODELS_H

enum MODEL{
  NO_MODEL,
  MIE_CYLINDER,
  SHELF,
  NONSHELF,
  LAYER
};

enum MODE{
  D_X, //x方向に線積分
  D_Y, //y方向に線積分
  D_XY //面積分 
};

extern void setModel(enum MODEL model);
extern double models_eps(double x, double y, enum MODE mode);

#endif
