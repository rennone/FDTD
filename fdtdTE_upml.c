#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fdtdTE_upml.h"
#include "field.h"
#include "models.h"

//Ex(i+0.5,j) -> Ex[i,j]
//Ey(i,j+0.5) -> Ey[i,j]
//Hz(i+0.5,j+0.5) -> Hz[i,j]
static double complex *Ex = NULL;
static double complex *Jx = NULL;
static double complex *Dx = NULL;

static double complex *Ey = NULL;
static double complex *Jy = NULL;
static double complex *Dy = NULL;

static double complex *Hz = NULL;
static double complex *Mz = NULL;
static double complex *Bz = NULL;

static double *C_JX = NULL, *C_JY = NULL, *C_MZ= NULL;
static double *C_JXHZ = NULL, *C_JYHZ = NULL, *C_MZEXEY= NULL;
static double *C_DX=NULL, *C_DY=NULL, *C_BZ=NULL;
static double *C_DXJX0=NULL, *C_DXJX1=NULL;
static double *C_DYJY0=NULL, *C_DYJY1=NULL;
static double *C_BZMZ0=NULL, *C_BZMZ1=NULL;

static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_HZ=NULL;

static const int timeStep = 500;
static double complex *Uz[360], *Wx[360], *Wy[360]; 
//------prototype--------//
static void update(void);
static void finish(void);
static void init(void);
static inline void calcE(void);
static inline void calcJD(void);
static inline void calcH(void);
static inline void calcMB(void);
//:public-------------------------------//

void (* fdtdTE_upml_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTE_upml_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTE_upml_getInit(void))(void)
{
  return init;
}

double complex* fdtdTE_upml_getEx(void){
  return Ex;
}

double complex* fdtdTE_upml_getEy(void){
  return Ey;
}

double complex* fdtdTE_upml_getHz(void)
{
  return Hz;
}

double* fdtdTE_upml_getEps()
{
  return EPS_EY;
}

//:private ----------------------------//
static inline double complex EX(const int i, const int j)
{
  return Ex[ind(i,j)];
}

static inline double complex EY(const int i, const int j)
{
  return Ey[ind(i,j)];
}

static inline double complex HZ(const int i, const int j)
{
  return Hz[ind(i,j)];
}

static inline double complex JX(const int i, const int j)
{
  return Jx[ind(i,j)];
}

static inline double complex JY(const int i, const int j)
{
  return Jy[ind(i,j)];
}

static inline double complex MZ(const int i, const int j)
{
  return Mz[ind(i,j)];
}

static inline double complex DX(const int i, const int j)
{
  return Dx[ind(i,j)];
}

static inline double complex DY(const int i, const int j)
{
  return Dy[ind(i,j)];
}

static inline double complex BZ(const int i, const int j)
{
  return Bz[ind(i,j)];
}

//Coefficient of Jx in calcJ
static inline double CJX(const int i, const int j)
{
  return C_JX[ind(i,j)];
}

//Coefficient of Jy in calcJ
static inline double CJY(const int i, const int j)
{
  return C_JY[ind(i,j)];
}

//Coefficient of Hx in calcJx
static inline double CJXHZ(const int i, const int j)
{
  return C_JXHZ[ind(i,j)];
}

//Coefficient of Hz in calcJy
static inline double CJYHZ(const int i, const int j)
{
  return C_JYHZ[ind(i,j)];
}

//Coefficient of Mz in calcM
static inline double CMZ(const int i, const int j)
{
  return C_MZ[ind(i,j)];
}

//Coefficient of Ex and Ey in calcM
static inline double CMZEXEY(const int i, const int j)
{
  return C_MZEXEY[ind(i,j)];
}

//Coefficient of Dx in calcD
static inline double CDX(const int i, const int j)
{
  return C_DX[ind(i,j)];
}

//Coefficient of Jx(n) in calcD
static inline double CDXJX0(const int i, const int j)
{
  return C_DXJX0[ind(i,j)];
}

//Coefficient of Jx(n+1) in calcD
static inline double CDXJX1(const int i, const int j)
{
  return C_DXJX1[ind(i,j)];
}

//Coefficient of Dy in calcD
static inline double CDY(const int i, const int j)
{
  return C_DY[ind(i,j)];
}

//Coefficient of Jy(n) in calcD
static inline double CDYJY0(const int i, const int j)
{
  return C_DYJY0[ind(i,j)];
}
//Coefficient of Jy(n+1) in calcD
static inline double CDYJY1(const int i, const int j)
{
  return C_DYJY1[ind(i,j)];
}

//Coefficient of Bz in calcB
static inline double CBZ(const int i, const int j)
{
  return C_BZ[ind(i,j)];
}

//Coefficient of Mz(n) in calcB
static inline double CBZMZ0(const int i, const int j)
{
  return C_BZMZ0[ind(i,j)];
}

//Coefficient of Mz(n+1) in calcB
static inline double CBZMZ1(const int i, const int j)
{
  return C_BZMZ1[ind(i,j)];
}

static inline double EPSHZ(const int i, const int j)
{
  return EPS_HZ[ind(i,j)];
}

static inline double EPSEX(const int i, const int j)
{
  return EPS_EX[ind(i,j)];
}

static inline double EPSEY(const int i, const int j)
{
  return EPS_EY[ind(i,j)];
}
//---------------------------------------------------//


//-----------------memory allocate-------------//
static void init(){  
  Ex = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Ey = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Hz = (double complex*)malloc(sizeof(double complex)*N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  Jx = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Jy = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Mz = (double complex*)malloc(sizeof(double complex)*N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  Dx = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Dy = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Bz = (double complex*)malloc(sizeof(double complex)*N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  C_JX = (double *)malloc(sizeof(double)*N_CELL);
  C_JY = (double *)malloc(sizeof(double)*N_CELL);
  C_MZ = (double *)malloc(sizeof(double)*N_CELL);  
  C_JXHZ = (double *)malloc(sizeof(double)*N_CELL);
  C_JYHZ = (double *)malloc(sizeof(double)*N_CELL);
  C_MZEXEY = (double *)malloc(sizeof(double)*N_CELL);

  C_DX = (double *)malloc(sizeof(double)*N_CELL);
  C_DY = (double *)malloc(sizeof(double)*N_CELL);
  C_BZ = (double *)malloc(sizeof(double)*N_CELL);  
  C_DXJX0 = (double *)malloc(sizeof(double)*N_CELL);
  C_DXJX1 = (double *)malloc(sizeof(double)*N_CELL);
  C_DYJY0 = (double *)malloc(sizeof(double)*N_CELL);
  C_DYJY1 = (double *)malloc(sizeof(double)*N_CELL);
  C_BZMZ0 = (double *)malloc(sizeof(double)*N_CELL);
  C_BZMZ1 = (double *)malloc(sizeof(double)*N_CELL);
  
  EPS_EX = (double *)malloc(sizeof(double)*N_CELL);
  EPS_EY = (double *)malloc(sizeof(double)*N_CELL);
  EPS_HZ = (double *)malloc(sizeof(double)*N_CELL);

  memset(Ex, 0, sizeof(double complex)*N_CELL);
  memset(Ey, 0, sizeof(double complex)*N_CELL);
  memset(Hz,0, sizeof(double complex)*N_CELL);
  memset(Jx, 0, sizeof(double complex)*N_CELL);
  memset(Jy, 0, sizeof(double complex)*N_CELL);
  memset(Mz,0, sizeof(double complex)*N_CELL);
  memset(Dx, 0, sizeof(double complex)*N_CELL);
  memset(Dy, 0, sizeof(double complex)*N_CELL);
  memset(Bz,0, sizeof(double complex)*N_CELL);


  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_y;
  double sig_ey_x, sig_ey_y;
  double sig_hz_x, sig_hz_y;
  double R = 1.0e-8;
  double M = 2.0;
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);
  int i,j;
  for(i=0; i<N_PX; i++){
    for(j=0; j<N_PY; j++){      
      int k = ind(i,j);
      EPS_EX[k] = models_eps(i+0.5,j, D_Y);
      EPS_EY[k] = models_eps(i,j+0.5, D_X);
      EPS_HZ[k] = 0.5*(models_eps(i+0.5,j+0.5, D_X) + models_eps(i+0.5,j+0.5, D_Y));

      sig_ex_x = sig_max*field_sigmaX(i+0.5,j);
      sig_ex_y = sig_max*field_sigmaY(i+0.5,j);
      sig_ey_x = sig_max*field_sigmaX(i,j+0.5);
      sig_ey_y = sig_max*field_sigmaY(i,j+0.5);
      sig_hz_x = sig_max*field_sigmaX(i+0.5,j+0.5);
      sig_hz_y = sig_max*field_sigmaY(i+0.5,j+0.5);

      double eps = EPSILON_0_S;
      double sig_z = 0;

      C_JX[k] = (2*eps - sig_ex_y)/(2*eps + sig_ex_y);
      C_JXHZ[k] = (2*eps)/(2*eps + sig_ex_y);
      C_DX[k] = (2*eps - sig_z) / (2*eps + sig_z);
      C_DXJX1[k] = (2*eps + sig_ex_x) / (2*eps + sig_z);
      C_DXJX0[k] = (2*eps - sig_ex_x) / (2*eps + sig_z);

      C_JY[k] = (2*eps - sig_z) / (2*eps + sig_z);
      C_JYHZ[k] = (2*eps)/(2*eps + sig_z);
      C_DY[k] = (2*eps - sig_ey_x) / (2*eps + sig_ey_x);
      C_DYJY1[k] = (2*eps + sig_ey_y) / (2*eps + sig_ey_x);
      C_DYJY0[k] = (2*eps - sig_ey_y) / (2*eps + sig_ey_x);

      C_MZ[k] = (2*eps - sig_hz_x) / (2*eps + sig_hz_x );
      C_MZEXEY[k] = (2*eps) / (2*eps + sig_hz_x);
      C_BZ[k] = (2*eps - sig_hz_y) / (2*eps + sig_hz_y);
      C_BZMZ1[k] = (2*eps + sig_z) / (2*eps + sig_hz_y);
      C_BZMZ0[k] = (2*eps - sig_z) / (2*eps + sig_hz_y);
    }
  }  
}

//---------------------メモリの解放--------------------//
static void finish(){
  if(Ex != NULL){    free(Ex); Ex = NULL;}  
  if(Ey != NULL){    free(Ey); Ey = NULL;}  
  if(Hz != NULL){    free(Hz); Hz = NULL;}

  if(C_JX!= NULL){    free(C_JX);  C_JX = NULL;}
  if(C_JXHZ!= NULL){   free(C_JXHZ); C_JXHZ = NULL;}
  if(C_DX!= NULL){   free(C_DX); C_DX = NULL;}
  if(C_DXJX0 != NULL){   free(C_DXJX0); C_DXJX0 = NULL;}
  if(C_DXJX1 != NULL){   free(C_DXJX1); C_DXJX1 = NULL;}
  
  if(C_JY!= NULL){    free(C_JY);  C_JY = NULL;}
  if(C_JYHZ!= NULL){   free(C_JYHZ); C_JYHZ = NULL;}
  if(C_DY!= NULL){   free(C_DY); C_DY = NULL;}
  if(C_DYJY0 != NULL){   free(C_DYJY0); C_DYJY0 = NULL;}
  if(C_DYJY1 != NULL){   free(C_DYJY1); C_DYJY1 = NULL;}
  
  if(C_MZ!= NULL){    free(C_MZ);  C_MZ = NULL;}
  if(C_MZEXEY!= NULL){   free(C_MZEXEY); C_MZEXEY = NULL;}
  if(C_BZ != NULL){   free(C_BZ); C_BZ = NULL;}
  if(C_BZMZ0 != NULL){   free(C_BZMZ0); C_BZMZ0 = NULL;}
  if(C_BZMZ1 != NULL){   free(C_BZMZ1); C_BZMZ1 = NULL;}
  
  if(EPS_EX != NULL)   free(EPS_EX);
  if(EPS_EY != NULL)   free(EPS_EY);
  if(EPS_HZ != NULL)   free(EPS_HZ);
}

//Standard Scattered Wave
static void scatteredWave(double complex *p, double *eps){
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s = field_getK();  
  double rad = field_getWaveAngle()*M_PI/180;	//ラジアン変換
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく
  
  int i,j;
  for(i=N_PML; i<N_X+N_PML; i++){
    for(j=N_PML; j<N_Y+N_PML; j++){
      double ikx = i*ks_cos + j*ks_sin; //k_s*(i*cos + j*sin)
      p[ind(i,j)] += ray_coef*(EPSILON_0_S/eps[ind(i,j)] - 1)*( cos(ikx-w_s*time) + I*sin(ikx-w_s*time) );
    }
  }
}

static inline void update(void)
{
  calcJD();

  calcE();

  scatteredWave(Ey, EPS_EY);
  
  calcMB();

  calcH();
}

static inline void calcJD(void)
{
  int i,j;
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      double complex nowJx = JX(i,j);
      Jx[ind(i,j)] = CJX(i,j)*JX(i,j) + CJXHZ(i,j)*(HZ(i,j) - HZ(i,j-1));
      Dx[ind(i,j)] = CDX(i,j)*DX(i,j) + CDXJX1(i,j)*JX(i,j) - CDXJX0(i,j)*nowJx;
    }
  }
  
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      double complex nowJy = JY(i,j);
      Jy[ind(i,j)] = CJY(i,j)*JY(i,j) + CJYHZ(i,j)*(-HZ(i,j) + HZ(i-1,j));
      Dy[ind(i,j)] = CDY(i,j)*DY(i,j) + CDYJY1(i,j)*JY(i,j) - CDYJY0(i,j)*nowJy;
    }
  }
}

static inline void calcE(void)
{
  int i,j;
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ex[ind(i,j)] = DX(i,j)/EPSEX(i,j);

  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ey[ind(i,j)] = DY(i,j)/EPSEY(i,j);
}

static inline void calcMB(void)
{
  int i,j;
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      double complex nowMz = MZ(i,j);
      Mz[ind(i,j)] = CMZ(i,j)*MZ(i,j) - CMZEXEY(i,j)*(EY(i+1,j) - EY(i,j) - EX(i,j+1) + EX(i,j));
      Bz[ind(i,j)] = CBZ(i,j)*BZ(i,j) + CBZMZ1(i,j)*MZ(i,j) - CBZMZ0(i,j)*nowMz; 
    }
  }
}

static inline void calcH(void)
{
  int i,j;
   for(i=1; i<N_PX-1; i++)     
    for(j=1; j<N_PY-1; j++)
      Hz[ind(i,j)] = BZ(i,j)/MU_0_S;   
}

//---------------------- ntff--------------------//
static inline void ntffCoef(double time, double timeShift, int *m, double *a, double *b, double *ab)
{  
  double t = time + timeShift;
  *m = floor(t + 0.5);
  *a = (0.5 + t - *m);
  *b = 1.0-*a;
  *ab = *a-*b;
}

static void ntff(void)
{
  const double C = LIGHT_SPEED_S;  
  const int cx = N_PX/2; //a center of field is origin
  const int cy = N_PY/2;  
  const int offset = 5;		// closed line offset
  const int tp = N_PY - N_PML - offset;	//上から-5
  const int bm = N_PML + offset;	//下から5
  const int rt = N_PX - N_PML - offset;	//右から-5
  const int lt = N_PML + offset;	//左から5
  const int len = tp-bm+1; //length
  const double RFperC = len*2;  //C = 1/√2,Rf=len*√2 => Rf/C = len*2  
  const double coef = 1.0/(4*M_PI*C);
  
  int m_e, m_h;
  double timeE = field_getTime() - 1;  //t - Δt  todo
  double timeH = field_getTime() - 0.5;  //t - Δt/2  todo
  double a_e, b_e, ab_e, a_h, b_h, ab_h;

  int ang;
  int i, j;
  for(ang=0; ang<360; ang++){
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);    
    //bottom side
    //normal vector n is (0,-1)
    //Js = n × H = (-Hz, 0, 0)  Ms = E × n = (0, 0, -Ex)

    for(i=lt; i<rt; i++){
      double r2x = i-cx+0.5, r2y = bm-cy;            //distance between origin to cell
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_e, &a_e, &b_e, &ab_e);
      
      double complex ex = -EX(i,bm);
      double complex hz = -0.5*( HZ(i,bm) + HZ(i,bm-1) );
          
      Uz[ang][m_e-1] += ex*b_e*coef;
      Wx[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ex*ab_e*coef;
      Wx[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ex*a_e*coef;
      Wx[ang][m_h+1] -= hz*a_h*coef;      
    }

    //top side
    //normal vector n is (0,1)
    //Js = n × H = (Hz, 0, 0)  Ms = E × n = (0, 0, Ex)
    for(i=lt; i<rt; i++){
      double r2x = i-cx+0.5, r2y = tp-cy;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;      
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      double complex ex = EX(i,bm);
      double complex hz = 0.5*( HZ(i,tp) + HZ(i,tp-1) );

      Uz[ang][m_e-1] += ex*b_e*coef;
      Wx[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ex*ab_e*coef;
      Wx[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ex*a_e*coef;
      Wx[ang][m_h+1] -= hz*a_h*coef;            
    }
    
    //left side
    //normal vector n is (-1,0)
    //Js = n × H = (0,Hz,0)    Ms = E × n = (0,0,Ey)
    for(j=bm; j<tp; j++){
      double r2x = lt-cx, r2y = j-cy+0.5;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);
      
      double complex ey = EY(lt,j);
      double complex hz = 0.5*( HZ(lt,j) + HZ(lt-1,j) );
      
      Uz[ang][m_e-1] += ey*b_e*coef;
      Wy[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ey*ab_e*coef;
      Wy[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ey*a_e*coef;
      Wy[ang][m_h+1] -= hz*a_h*coef;
    }

    //right side
    //normal vector n is (1,0)
    //Js = n × H = (0,-Hz,0)    Ms = E × n = (0,0,-Ey)    
    for(j=bm; j<tp; j++){
      double r2x = rt-cx, r2y = j-cy+0.5;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      double complex ey = -EY(lt,j);
      double complex hz = -0.5*( HZ(lt,j) + HZ(lt-1,j) );
      Uz[ang][m_e-1] += ey*b_e*coef;
      Wy[ang][m_h-1] += hz*b_h*coef;
      Uz[ang][m_e]   += ey*ab_e*coef;
      Wy[ang][m_h]   += hz*ab_h*coef;
      Uz[ang][m_e+1] -= ey*a_e*coef;
      Wy[ang][m_h+1] -= hz*a_h*coef;      
    }

  }
}
