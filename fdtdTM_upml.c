#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fdtdTM_upml.h"
#include "field.h"
#include "models.h"

//Ez(i    , j    ) -> Ez[i,j];
//Hx(i    , j+0.5) -> Hx[i,j];
//Hy(i+0.5, j    ) -> Hy[i,j];
#define M_PI 3.1415926535897932384626433832795

static double complex *Ez = NULL;
static double complex *Jz = NULL;
static double complex *Dz = NULL;

static double complex *Hx = NULL;
static double complex *Mx = NULL;
static double complex *Bx = NULL;

static double complex *Hy = NULL;
static double complex *My = NULL;
static double complex *By = NULL;

//--------------------for NTFF --------------------//
static double complex *Ux,*Uy,*Wz;
//------------------------------------------------//

static double *C_JZ = NULL, *C_MX = NULL, *C_MY = NULL;
static double *C_JZHXHY = NULL, *C_MXEZ = NULL, *C_MYEZ = NULL;
static double *C_DZ = NULL, *C_BX = NULL, *C_BY = NULL;
static double *C_DZJZ0 = NULL, *C_DZJZ1 = NULL;
static double *C_BXMX0 = NULL, *C_BXMX1 = NULL;
static double *C_BYMY0 = NULL, *C_BYMY1 = NULL;

static double *EPS_EZ = NULL, *EPS_HX = NULL, *EPS_HY = NULL;

static void update(void);
static void finish(void);
static void init(void);

static inline void calcJD(void);
static inline void calcE(void);

static inline void calcMB(void);
static inline void calcH(void);
static void ntff(void);
static void outputNTFF();
  
//:public
void (* fdtdTM_upml_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTM_upml_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTM_upml_getInit(void))(void)
{
  return init;
}

double complex* fdtdTM_upml_getHx(void){
  return Hx;
}

double complex* fdtdTM_upml_getHy(void){
  return Hy;
}

double complex* fdtdTM_upml_getEz(void){
  return Ez;
}

double* fdtdTM_upml_getEps()
{
  return EPS_EZ;
}
//:private
//--------------------getter--------------------//
static inline double complex EZ(const int i, const int j)
{
  return Ez[ind(i,j)];
}

static inline double complex HX(const int i, const int j)
{
  return Hx[ind(i,j)];
}

static inline double complex HY(const int i, const int j)
{
  return Hy[ind(i,j)];
}

static inline double complex DZ(const int i, const int j)
{
  return Dz[ind(i,j)];
}

static inline double complex BX(const int i, const int j)
{
  return Bx[ind(i,j)];
}

static inline double complex BY(const int i, const int j)
{
  return By[ind(i,j)];
}

static inline double complex JZ(const int i, const int j)
{
  return Jz[ind(i,j)];
}

static inline double complex MX(const int i, const int j)
{
  return Mx[ind(i,j)];
}

static inline double complex MY(const int i, const int j)
{
  return My[ind(i,j)];
}

//Coefficient of Jz in calcJ
static inline double CJZ(const int i, const int j)
{
  return C_JZ[ind(i,j)];
}

//Coefficient of Hx and Hy in calcJ
static inline double CJZHXHY(const int i, const int j)
{
  return C_JZHXHY[ind(i,j)];
}

//Coefficient of Dz in calcD
static inline double CDZ(const int i, const int j)
{
  return C_DZ[ind(i,j)];
}

//Coefficient of Jz(n) in calcD
static inline double CDZJZ0(const int i, const int j)
{
  return C_DZJZ0[ind(i,j)];
}

//Coefficient of Jz(n+1) in calcD
static inline double CDZJZ1(const int i, const int j)
{
  return C_DZJZ1[ind(i,j)];
}

//Coefficient of Mx in calcM
static inline double CMX (const int i, const int j)
{
  return  C_MX[ind(i,j)];
}

//Coefficient of Ez in calcM
static inline double CMXEZ(const int i, const int j)
{
  return C_MXEZ [ind(i,j)];
}

//Coefficient of Bx in calcB
static inline double CBX(const int i, const int j)
{
  return C_BX [ind(i,j)];
}

//Coefficient of Mx(n+1) in calcB
static inline double CBXMX1(const int i, const int j)
{
  return C_BXMX1 [ind(i,j)];
}

//Coefficient of Mx(n) in calcB
static inline double CBXMX0(const int i, const int j)
{
  return C_BXMX0[ind(i,j)];
}

//Coefficient of My in calcM
static inline double CMY (const int i, const int j)
{
  return  C_MY[ind(i,j)];
}

//Coefficient of Ez in calcM
static inline double CMYEZ(const int i, const int j)
{
  return C_MYEZ [ind(i,j)];
}

//Coefficient of By in calcB
static inline double CBY(const int i, const int j)
{
  return C_BY[ind(i,j)];
}

//Coefficient of My(n+1) in calcB
static inline double CBYMY1(const int i, const int j)
{
  return C_BYMY1[ind(i,j)];
}

//Coefficient of My(n) in calcB
static inline double CBYMY0(const int i, const int j)
{
  return C_BYMY0[ind(i,j)];
}

//Epsilon of Ez
static inline double EPSEZ(const int i, const int j)
{
  return EPS_EZ[ind(i,j)];
}

//Epsilon of Hx
static inline double EPSHX(const int i, const int j)
{
  return EPS_HX[ind(i,j)];
}

//Epsilon of Hy
static inline double EPSHY(const int i, const int j)
{
  return EPS_HY[ind(i,j)];
}

static inline void swap(double *a, double *b)
{
  double tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

//Initialize
static void init(void)
{
  Ez = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Dz = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Jz = (double complex*)malloc(sizeof(double complex)*N_CELL);
  
  Hx = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Mx = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Bx = (double complex*)malloc(sizeof(double complex)*N_CELL);
  
  Hy = (double complex*)malloc(sizeof(double complex)*N_CELL);
  My = (double complex*)malloc(sizeof(double complex)*N_CELL);
  By = (double complex*)malloc(sizeof(double complex)*N_CELL);

  
  int step = (field_getNTFFInfo()).step;
  Ux = (double complex*)malloc(sizeof(double complex)*360*step);
  Uy = (double complex*)malloc(sizeof(double complex)*360*step);
  Wz = (double complex*)malloc(sizeof(double complex)*360*step);
  
  C_JZ = (double *)malloc(sizeof(double)*N_CELL);
  C_MX = (double *)malloc(sizeof(double)*N_CELL);
  C_MY = (double *)malloc(sizeof(double)*N_CELL);

  C_DZ = (double *)malloc(sizeof(double)*N_CELL);
  C_BX = (double *)malloc(sizeof(double)*N_CELL);
  C_BY = (double *)malloc(sizeof(double)*N_CELL);

  C_JZHXHY = (double *)malloc(sizeof(double)*N_CELL);
  C_MXEZ = (double *)malloc(sizeof(double)*N_CELL);
  C_MYEZ = (double *)malloc(sizeof(double)*N_CELL);

  C_DZJZ0 = (double *)malloc(sizeof(double)*N_CELL);
  C_DZJZ1 = (double *)malloc(sizeof(double)*N_CELL);

  C_BXMX1 = (double *)malloc(sizeof(double)*N_CELL);
  C_BXMX0 = (double *)malloc(sizeof(double)*N_CELL);

  C_BYMY1 = (double *)malloc(sizeof(double)*N_CELL);
  C_BYMY0 = (double *)malloc(sizeof(double)*N_CELL);

  EPS_HY = (double *)malloc(sizeof(double)*N_CELL);
  EPS_HX = (double *)malloc(sizeof(double)*N_CELL);
  EPS_EZ = (double *)malloc(sizeof(double)*N_CELL);

  memset(Hx, 0, sizeof(double complex)*N_CELL);
  memset(Hy, 0, sizeof(double complex)*N_CELL);
  memset(Ez, 0, sizeof(double complex)*N_CELL);

  memset(Mx, 0, sizeof(double complex)*N_CELL);
  memset(My, 0, sizeof(double complex)*N_CELL);
  memset(Jz, 0, sizeof(double complex)*N_CELL);

  memset(Bx, 0, sizeof(double complex)*N_CELL);
  memset(By, 0, sizeof(double complex)*N_CELL);
  memset(Dz, 0, sizeof(double complex)*N_CELL);

  //Ez,, Hx, Hyそれぞれでσx,σyが違う(場所が違うから)
  double sig_ez_x, sig_ez_y;
  double sig_hx_x, sig_hx_y;
  double sig_hy_x, sig_hy_y;
  
  const double R = 1.0e-8;
  const double M = 2.0;
  //const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML/cos(M_PI/3)*log(R);
  int i,j;
  for(i=0; i<N_PX; i++){
    for(j=0; j<N_PY; j++){
      EPS_EZ[ind(i,j)] = models_eps(i,j, D_XY);     //todo D_X, D_Yにしなくていいのか?
      EPS_HX[ind(i,j)] = models_eps(i,j+0.5, D_XY);
      EPS_HY[ind(i,j)] = models_eps(i+0.5,j, D_XY);
      
      sig_ez_x = sig_max*field_sigmaX(i,j);
      sig_ez_y = sig_max*field_sigmaY(i,j);

      sig_hx_x = sig_max*field_sigmaX(i,j+0.5);
      sig_hx_y = sig_max*field_sigmaY(i,j+0.5);
      
      sig_hy_x = sig_max*field_sigmaX(i+0.5,j);
      sig_hy_y = sig_max*field_sigmaY(i+0.5,j);

      double sig_z = 0; // σz is zero in
      
      //Δt = 1  Κ_i = 1
      double eps = EPSILON_0_S;
      
      C_JZ[ind(i,j)]     = ( 2*eps - sig_ez_x) / (2*eps + sig_ez_x);
      C_JZHXHY[ind(i,j)] = ( 2*eps ) / (2*eps + sig_ez_x);
      C_DZ[ind(i,j)]     = ( 2*eps - sig_ez_y)  / (2*eps + sig_ez_y);      
      C_DZJZ1[ind(i,j)]  = ( 2*eps + sig_z)/(2*eps + sig_ez_y);
      C_DZJZ0[ind(i,j)]  = ( 2*eps - sig_z)/(2*eps + sig_ez_y);

      C_MX[ind(i,j)]    = (2*eps - sig_hx_y) / (2*eps + sig_hx_y);
      C_MXEZ[ind(i,j)]  = (2*eps) / (2*eps + sig_hx_y);
      C_BX[ind(i,j)]    = (2*eps - sig_z) / (2*eps + sig_z);
      C_BXMX1[ind(i,j)] = (2*eps + sig_hx_x) / (2*eps + sig_z);
      C_BXMX0[ind(i,j)] = (2*eps - sig_hx_x) / (2*eps + sig_z);
      
      C_MY[ind(i,j)]    = (2*eps - sig_z) / (2*eps + sig_z);
      C_MYEZ[ind(i,j)]  = (2*eps) / (2*eps + sig_z);      
      C_BY[ind(i,j)]    = (2*eps - sig_hy_x) / (2*eps + sig_hy_x);
      C_BYMY1[ind(i,j)] = (2*eps + sig_hy_y) / (2*eps + sig_hy_x);
      C_BYMY0[ind(i,j)] = (2*eps - sig_hy_y) / (2*eps + sig_hy_x);      
    }
  }
}

//Finish
static void finish(void)
{
  outputNTFF();
  
  if(Ez != NULL){   free(Ez); Ez = NULL;  }  
  if(Hx != NULL){   free(Hx); Hx = NULL;  }
  if(Hy != NULL){   free(Hy); Hy = NULL;  }
  
  if(Ux != NULL){   free(Ux); Ux = NULL;  }  
  if(Uy != NULL){   free(Uy); Uy = NULL;  }
  if(Wz != NULL){   free(Wz); Wz = NULL;  }
  
  if(C_JZ != NULL){ free(C_JZ); C_JZ = NULL;  }
  if(C_JZHXHY != NULL){ free(C_JZHXHY); C_JZHXHY = NULL;  }
  
  if(C_DZ != NULL){ free(C_DZ); C_DZ = NULL;  }
  if(C_DZJZ1 != NULL){ free(C_DZJZ1); C_DZJZ1 = NULL;  }
  if(C_DZJZ0 != NULL){ free(C_DZJZ0); C_DZJZ0 = NULL;  }

  if(C_MX != NULL){ free(C_MX); C_MX = NULL;  }
  if(C_MXEZ != NULL){ free(C_MXEZ); C_MXEZ = NULL;  }
  
  if(C_BX != NULL){ free(C_BX); C_BX = NULL;  }
  if(C_BXMX1 != NULL){ free(C_BXMX1); C_BXMX1 = NULL;  }
  if(C_BXMX0 != NULL){ free(C_BXMX0); C_BXMX0 = NULL;  }
  
  if(C_MY != NULL){ free(C_MY); C_MY = NULL;  }
  if(C_MYEZ != NULL){ free(C_MYEZ); C_MYEZ = NULL;  }

  if(C_BY != NULL){ free(C_BY); C_BY = NULL;  }
  if(C_BYMY1 != NULL){ free(C_BYMY1); C_BYMY1 = NULL;  }
  if(C_BYMY0 != NULL){ free(C_BYMY0); C_BYMY0 = NULL;  }  
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

//Update
static void update(void)
{
  calcJD();
  calcE();
  scatteredWave(Ez, EPS_EZ);
  calcMB();  
  calcH();

  ntff();
}

//calculate J and D
static inline void calcJD()
{
  int i,j;
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      /*
	double complex nowJz = JZ(i,j);
      Jz[ind(i,j)] = CJZ(i,j)*JZ(i,j) + CJZHXHY(i,j)*(+HY(i,j) - HY(i-1,j) - HX(i,j) + HX(i,j-1) );
      Dz[ind(i,j)] = CDZ(i,j)*DZ(i,j) + CDZJZ1(i,j)*JZ(i,j) - CDZJZ0(i,j)*nowJz;
      */
      const int k = ind(i,j);
      double complex nowJz = Jz[k];
      Jz[k] = C_JZ[k]*Jz[k] + C_JZHXHY[k]*(+Hy[k] - HY(i-1,j) - Hx[k] + HX(i,j-1) );
      Dz[k] = C_DZ[k]*Dz[k] + C_DZJZ1[k]*Jz[k] - C_DZJZ0[k]*nowJz;
    }
  }
}

//calculate E 
static inline void calcE()
{
  double epsilon = EPSILON_0_S;
  int i,j;
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      const int k = ind(i,j);
      Ez[k] = Dz[k]/EPS_EZ[k];
    }
  }
}

//calculate M and B
static inline void calcMB()
{
  int i,j;
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      /*
      double complex nowMx = MX(i,j);
      Mx[ind(i,j)] = CMX(i,j)*MX(i,j) - CMXEZ(i,j)*(EZ(i,j+1) - EZ(i,j));
      Bx[ind(i,j)] = CBX(i,j)*BX(i,j) + CBXMX1(i,j)*MX(i,j) - CBXMX0(i,j)*nowMx;
      */
      const int k = ind(i,j);
      double complex nowMx = Mx[k];
      Mx[k] = C_MX[k]*Mx[k] - C_MXEZ[k]*(EZ(i,j+1) - Ez[k]);
      Bx[k] = C_BX[k]*Bx[k] + C_BXMX1[k]*Mx[k] - C_BXMX0[k]*nowMx;
    }
  }
  
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      /*
      double complex nowMy = MY(i,j);
      My[ind(i,j)] = CMY(i,j)*MY(i,j) - CMYEZ(i,j)*(-EZ(i+1,j) + EZ(i,j));
      By[ind(i,j)] = CBY(i,j)*BY(i,j) + CBYMY1(i,j)*MY(i,j) - CBYMY0(i,j)*nowMy;
      */
      const int k = ind(i,j);
      double complex nowMy = My[k];
      My[k] = C_MY[k]*My[k] - C_MYEZ[k]*(-EZ(i+1,j) + Ez[k]);
      By[k] = C_BY[k]*By[k] + C_BYMY1[k]*My[k] - C_BYMY0[k]*nowMy;
    }
  }
}

//calculate H
static inline void calcH()
{
  int i,j;
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      const int k = ind(i,j);
      Hx[k] = Bx[k]/MU_0_S;
    }
  }
  
  for(i=1; i<N_PX-1; i++){
    for(j=1; j<N_PY-1; j++){
      const int k = ind(i,j);
      Hy[k] = By[k]/MU_0_S;
    }
  }
}

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

  ntffInfo nInfo = field_getNTFFInfo();
  const int tp = nInfo.top;
  const int bm = nInfo.bottom;
  const int rt = nInfo.right;
  const int lt = nInfo.left;
  const double RFperC = nInfo.RFperC;
  
  const double coef = 1.0/(4*M_PI*C);
  int m_e, m_h;
  double a_e, b_e, ab_e, a_h, b_h, ab_h;
  
  double timeE = field_getTime() - 1;  //t - Δt
  double timeH = field_getTime() - 0.5;  //t - Δt/2
  int ang;
  for(ang=0; ang<360; ang++){
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);
    int i,j;
    //bottom side
    //normal vector n is (0,-1)
    //Js = n × H = ( 0, 0, Hx)  Ms = E × n = (Ez, 0,  0)
    //Js -> W                   Ms -> U
    for(i=lt; i<rt; i++){
      double r2x = i-cx+0.5, r2y = bm-cy;            //distance between origin to cell
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;      
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);
      
      double complex ez = 0.5*(EZ(i,bm) + EZ(i+1, bm));
      double complex hx = 0.25*(HX(i,bm) + HX(i+1,bm) + HX(i,bm-1) + HX(i+1,bm-1) );

      int k = ang*nInfo.step;

      Ux[k+m_e-1] += ez*b_e*coef;
      Ux[k+m_e]   += ez*ab_e*coef;
      Ux[k+m_e+1] -= ez*a_e*coef;
      Wz[k+m_h-1] += hx*b_h*coef;
      Wz[k+m_h]   += hx*ab_h*coef;
      Wz[k+m_h+1] -= hx*a_h*coef;
      
    }

    //top side
    //normal vector n is (0,1)
    //Js = n × H = (0, 0,-Hx)  Ms = E × n = (-Ez, 0,  0)
    //Js -> W                   Ms -> U
    for(i=lt; i<rt; i++){
      double r2x = i-cx+0.5, r2y = tp-cy;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;      
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);
        
      double complex ez = -0.5*(EZ(i,tp) + EZ(i+1, tp));
      double complex hx = -0.25*(HX(i,tp) + HX(i+1,tp) + HX(i,tp-1) + HX(i+1,tp-1) );


      int k = ang*nInfo.step;
      Ux[k+m_e-1] += ez*b_e*coef;
      Ux[k+m_e]   += ez*ab_e*coef;
      Ux[k+m_e+1] -= ez*a_e*coef;
      Wz[k+m_h-1] += hx*b_h*coef;
      Wz[k+m_h]   += hx*ab_h*coef;
      Wz[k+m_h+1] -= hx*a_h*coef;
    }
    
    //left side
    //normal vector n is (-1,0)
    //Js = n × H = (0,0,-Hy)    Ms = E × n = (0,-Ez,0)
    //Js -> W                   Ms -> U
    for(j=bm; j<tp; j++){
      double r2x = lt-cx, r2y = j-cy+0.5;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);
      
      double complex ez = -0.5*(EZ(lt,j) + EZ(lt, j+1));
      double complex hy = -0.25*( HY(lt,j) + HY(lt-1,j) + HY(lt,j+1) + HY(lt-1,j+1) );

      int k = ang*nInfo.step;
      Uy[k+m_e-1] += ez*b_e*coef;
      Uy[k+m_e]   += ez*ab_e*coef;
      Uy[k+m_e+1] -= ez*a_e*coef;
      
      Wz[k+m_h-1] += hy*b_h*coef;      
      Wz[k+m_h]   += hy*ab_h*coef;      
      Wz[k+m_h+1] -= hy*a_h*coef;
    }

    //right side
    //normal vector n is (1,0)
    //Js = n × H = (0, 0,Hy)  Ms = E × n = ( 0,Ez,0)
    //Js -> W                   Ms -> U
    for(j=bm; j<tp; j++){
      double r2x = rt-cx, r2y = j-cy+0.5;
      double timeShift = -(r1x*r2x + r1y*r2y)/C + RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);      
      double complex ez = 0.5*(EZ(rt,j) + EZ(rt, j+1));
      double complex hy = 0.25*( HY(rt,j) + HY(rt-1,j) + HY(rt,j+1) + HY(rt-1,j+1) );
      int k = ang*nInfo.step;
      Uy[k+m_e-1] += ez*b_e*coef;
      Uy[k+m_e]   += ez*ab_e*coef;
      Uy[k+m_e+1] -= ez*a_e*coef;
      Wz[k+m_h-1] += hy*b_h*coef;
      Wz[k+m_h]   += hy*ab_h*coef;      
      Wz[k+m_h+1] -= hy*a_h*coef;
    }

  }
}



static void outputNTFF()
{
  ntffInfo nInfo = field_getNTFFInfo();
  double complex *EPH = (double complex*)malloc(sizeof(double complex)*nInfo.step);

  FILE *fp_r, *fp_i;
  char fileName[256] = {'\0'};
  char folderName[256] = {'\0'};
  sprintf(folderName, "TM");
  /*
  mode_t mode = S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR | S_IWGRP;
  if( mkdir(folderName, mode) !=0){
    //error process
  }
  */
  sprintf(fileName, "%s/Ezf_r.txt",folderName);  
  if( (fp_r=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }
  
  sprintf(fileName, "%s/Ezf_i.txt",folderName);  
  if( (fp_i=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  double theta = 0;
  int ang;
  for(ang=0; ang<360; ang++){
    double phi = ang*M_PI/180.0;
    double _cos = cos(phi), _sin = sin(phi);
    int k= ang*nInfo.step;
    int i=0;
    for(i=0; i < nInfo.step; i++){
      //todo θとφがたぶん逆
      double complex UTH = Ux[k+i]*_cos;
      double complex UPH = Uy[k+i];
      double complex WTH = -Wz[k+i]*_sin;      
      double complex WPH = 0;
      
      EPH[i] = -Z_0_S*WPH+UTH;
      fprintf(fp_r,"%lf \n", creal(EPH[i]));
      fprintf(fp_i,"%lf \n", cimag(EPH[i]));      
    }
  }
  
  fclose(fp_r);
  fclose(fp_i);   
  free(EPH);
}
