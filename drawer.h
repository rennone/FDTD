#ifndef _DRAWER_H
#define _DRAWER_H
#include <complex.h>

//#define DEBUG

#ifdef DEBUG

enum COLOR_MODE{
  CREAL,
  CABS
};
extern void (*drawer_getDraw(void))(void);
extern void drawer_paintImage(int l, int b, int r, int t,int wid, int hei, double complex*, ...);
extern void drawer_paintModel(int l, int b, int r, int t,int wid, int hei, double *, ...);
extern void drawer_init(enum COLOR_MODE);
extern void drawer_finish(void);
extern void drawer_draw(void);
#endif //DEBUG

#endif //_DRAWER_H
