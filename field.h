#ifndef _FIELD_H
#define _FIELD_H
#include <stdio.h>
#include <complex.h>
#include "bool.h"

//入射波のモード
enum WAVE_MODE{
  POINT_LIGHT_IN_CENTER,  //中心に点光源
  SCATTER //散乱波
};

typedef struct ntffInfo{
  int top, bottom, left, right;
  int cx,cy;
  double RFperC;
  int step;
} ntffInfo;

//シミュレーション上の物理定数
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define C_S 0.7
static const double LIGHT_SPEED_S = 0.7;
static const double EPSILON_0_S = 1.0;
static const double MU_0_S = 1.0/0.7/0.7;
static const double Z_0_S = 1.42857142857; //√(1.0/0.7/0.7/1.0) = √(μ/ε);

extern int N_X;
extern int N_Y;
extern int N_CELL;
extern int H;
extern int N_PML;
extern int N_PX;
extern int N_PY;

//インデックスを取ってくる 
extern inline int ind(const int, const int);

//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
extern void setField(const int wid, const int hei, const int h, const int pml, const double lambda, const double step);

//pml用のσを取ってくる
extern double field_sigmaX(double x, double y);
extern double field_sigmaY(double x, double y);
extern double field_pmlCoef(double x, double y);
extern double field_pmlCoef_LXY(double x, double y);
extern inline double field_toCellUnit(const double);
extern inline double field_toPhisycalUnit(const double);
//extern inline double complex field_cbilinear(double complex *p, double x, double y);
//extern inline double field_dbilinear(double *p, double x, double y);
//---------------入射波---------------
extern double complex field_pointLight(void);
extern void field_defaultIncidence(double complex *p, double *eps);

//:NTFF
extern bool field_nextStep(void);


//:getter
extern double field_getLambda(void);
extern double field_getWaveAngle(void);
extern double field_getTime(void);
extern double field_getMaxTime(void);
extern ntffInfo field_getNTFFInfo(void);
extern double field_getOmega(void);
extern double field_getK(void);
extern double field_getRayCoef(void);

//:setter
extern void field_setDefaultIncidence(enum WAVE_MODE wm);

//output method
extern void field_outputElliptic(const char *fileName,double complex* data, double ox, double oy, double r); //
extern void field_outputAllData(const char *fileName,double complex* data); //

#endif
