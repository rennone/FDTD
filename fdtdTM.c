#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "field.h"
#include "fdtdTM.h"
#include "models.h"
#include "function.h"

//   メンバ変数    //
static double complex *Ezx = NULL;
static double complex *Ezy = NULL;
static double complex *Ez = NULL;
static double complex *Hx = NULL;
static double complex *Hy = NULL;
static double *C_EZ = NULL, *C_EZLH = NULL, *C_HXLY = NULL, *C_HYLX = NULL;
static double *C_EZX = NULL, *C_EZY = NULL, *C_EZXLX = NULL, *C_EZYLY = NULL;
static double *C_HX = NULL, *C_HY = NULL;
static double *EPS_EZ = NULL, *EPS_HX = NULL, *EPS_HY = NULL;

//  メンバ関数      //
static void update(void);
static void finish(void);
static void init(void);
static inline void calcE(void);
static inline void calcH(void);
static void ntffFrequency(void);

//---------------public method---------------//
//---------------getter------------------//
void (* fdtdTM_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTM_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTM_getInit(void))(void)
{
  return init;
}

double complex* fdtdTM_getHx(void){
  return Hx;
}

double complex* fdtdTM_getHy(void){
  return Hy;
}

double complex* fdtdTM_getEzx(void){
  return Ezx;
}

double complex* fdtdTM_getEzy(void){
  return Ezy;
}

double complex* fdtdTM_getEz(void){
  return Ez;
}

double* fdtdTM_getEps()
{
  return EPS_EZ;
}
//---------------getter------------------//
//---------------public method---------------//


//---------------private getter----------------//
static inline double complex EZ(const int i, const int j){
  return Ez[ind(i,j)];
}

static inline double complex HX(const int i, const int j){
  return Hx[ind(i,j)];
}

static inline double complex HY(const int i, const int j){
  return Hy[ind(i,j)];
}

static inline double complex EZX(const int i, const int j){
  return Ezx[ind(i,j)];
}

static inline double complex EZY(const int i, const int j){
  return Ezy[ind(i,j)];
}

static inline double CEZX(const int i, const int j){
  return C_EZX[ind(i,j)];
}

static inline double CEZY(const int i, const int j){
  return C_EZY[ind(i,j)];
}

static inline double CEZXLX(const int i, const int j){
  return C_EZXLX[ind(i,j)];
}

static inline double CEZYLY(const int i, const int j){
  return C_EZYLY[ind(i,j)];
}

static inline double CEZ(const int i, const int j){
  return C_EZ[ind(i,j)];
}

static inline double CEZLH(const int i, const int j){
  return C_EZLH[ind(i,j)];
}

static inline double CHXLY(const int i, const int j){
  return C_HXLY[ind(i,j)];
}

static inline double CHYLX(const int i, const int j){
  return C_HYLX[ind(i,j)];
}

static inline double CHX(const int i, const int j){
  return C_HX[ind(i,j)];
}

static inline double CHY(const int i, const int j){
  return C_HY[ind(i,j)];
}

static inline double EPSEZ(const int i, const int j){
  return EPS_EZ[ind(i,j)];
}

static inline double EPSHX(const int i, const int j){
  return EPS_HX[ind(i,j)];
}

static inline double EPSHY(const int i, const int j){
  return EPS_HY[ind(i,j)];
}

//----------private getter--------------//

static void init(){
  Hx = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Hy = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Ezy = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Ezx = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Ez = (double complex*)malloc(sizeof(double complex)*N_CELL);

  C_HX = (double *)malloc(sizeof(double)*N_CELL);
  C_HY = (double *)malloc(sizeof(double)*N_CELL);
  C_EZ = (double *)malloc(sizeof(double)*N_CELL);
  C_EZX = (double *)malloc(sizeof(double)*N_CELL);
  C_EZY = (double *)malloc(sizeof(double)*N_CELL);

  C_HXLY = (double *)malloc(sizeof(double)*N_CELL);
  C_HYLX = (double *)malloc(sizeof(double)*N_CELL);
  C_EZLH = (double *)malloc(sizeof(double)*N_CELL);

  C_EZXLX = (double *)malloc(sizeof(double)*N_CELL);
  C_EZYLY = (double *)malloc(sizeof(double)*N_CELL);

  EPS_HY = (double *)malloc(sizeof(double)*N_CELL);
  EPS_HX = (double *)malloc(sizeof(double)*N_CELL);
  EPS_EZ = (double *)malloc(sizeof(double)*N_CELL);
  
  memset(Hx, 0, sizeof(double complex)*N_CELL);
  memset(Hy, 0, sizeof(double complex)*N_CELL);
  memset(Ezx,0, sizeof(double complex)*N_CELL);
  memset(Ezy,0, sizeof(double complex)*N_CELL);
  memset(Ez,0, sizeof(double complex)*N_CELL);
  
  //Ez,, Hx, Hyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ez_x, sig_ez_y, sig_ez_xx, sig_ez_yy;
  double sig_hx_x, sig_hx_y, sig_hx_xx, sig_hx_yy;
  double sig_hy_x, sig_hy_y, sig_hy_xx, sig_hy_yy;

  double R = 1.0e-8;
  double M = 2.0;
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);
  int i,j;
  for(i=0; i<N_PX; i++){
    for(j=0; j<N_PY; j++){
      EPS_EZ[ind(i,j)] = models_eps(i,j, D_XY);     //todo D_X, D_Yにしなくていいのか?
      EPS_HX[ind(i,j)] = models_eps(i,j+0.5, D_Y);
      EPS_HY[ind(i,j)] = models_eps(i+0.5,j, D_X);
      
      sig_ez_x = sig_max*field_sigmaX(i,j);
      sig_ez_xx = MU_0_S/EPSILON_0_S * sig_ez_x;      
      sig_ez_y = sig_max*field_sigmaY(i,j);
      sig_ez_yy = MU_0_S/EPSILON_0_S * sig_ez_y;

      sig_hx_x = sig_max*field_sigmaX(i,j+0.5);
      sig_hx_xx = MU_0_S/EPSILON_0_S * sig_hx_x;
      sig_hx_y = sig_max*field_sigmaY(i,j+0.5);
      sig_hx_yy = MU_0_S/EPSILON_0_S * sig_hx_y;

      sig_hy_x = sig_max*field_sigmaX(i+0.5,j);
      sig_hy_xx = MU_0_S/EPSILON_0_S * sig_hy_x;      
      sig_hy_y = sig_max*field_sigmaY(i+0.5,j);
      sig_hy_yy = MU_0_S/EPSILON_0_S * sig_hy_y;

      //Δt = 1, μ(i,j) = μ0
      C_EZX[ind(i,j)]   = field_pmlCoef(EPSEZ(i,j), sig_ez_x);
      C_EZXLX[ind(i,j)] = field_pmlCoef_LXY(EPSEZ(i,j), sig_ez_x);

      C_EZY[ind(i,j)]   = field_pmlCoef(EPSEZ(i,j), sig_ez_y);
      C_EZYLY[ind(i,j)] = field_pmlCoef_LXY(EPSEZ(i,j), sig_ez_y);

      C_HX[ind(i,j)]    = field_pmlCoef(MU_0_S, sig_hx_yy);
      C_HXLY[ind(i,j)]  = field_pmlCoef_LXY(MU_0_S, sig_hx_yy);

      C_HY[ind(i,j)]    = field_pmlCoef(MU_0_S, sig_hy_xx);
      C_HYLX[ind(i,j)]  = field_pmlCoef_LXY(MU_0_S, sig_hy_xx);
    }
  }  
}

static void finish(){
  ntffFrequency();
  
  if(Hx != NULL){    free(Hx); Hx = NULL;}
  if(Hy != NULL){    free(Hy); Hy = NULL;}
  if(Ezx != NULL){    free(Ezx); Ezx = NULL;}
  if(Ezy != NULL){    free(Ezy); Ezy = NULL;}
  if(Ez != NULL){    free(Ez); Ez = NULL;}

  if(C_HY!= NULL){    free(C_HY); C_HY = NULL;}
  if(C_HX!= NULL){    free(C_HX); C_HX = NULL;}
  if(C_HXLY!= NULL){    free(C_HXLY); C_HXLY = NULL;}
  if(C_HYLX!= NULL){    free(C_HYLX); C_HYLX = NULL;}
  if(C_EZLH!= NULL){    free(C_EZLH); C_EZLH = NULL;}
  if(C_EZX != NULL){    free(C_EZX); C_EZX = NULL;}
  if(C_EZY != NULL){    free(C_EZY); C_EZY = NULL;}
  if(C_EZXLX != NULL){    free(C_EZXLX); C_EZXLX = NULL;}
  if(C_EZYLY != NULL){    free(C_EZYLY); C_EZYLY = NULL;}
  
  if(EPS_HX != NULL){    free(EPS_HX); EPS_HX = NULL;}
  if(EPS_HY != NULL){    free(EPS_HY); EPS_HY = NULL;}
  if(EPS_EZ != NULL){    free(EPS_EZ); EPS_EZ = NULL;}
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
      p[ind(i,j)] += ray_coef*(EPSILON_0_S/eps[ind(i,j)] - 1)*(
        cos(ikx-w_s*(time+0.5)) + I*sin(ikx-w_s*(time+0.5))
        -cos(ikx-w_s*(time-0.5)) - I*sin(ikx-w_s*(time-0.5))
        );
    }
  }
}

static void update(){
  calcE();

  scatteredWave(Ezx, EPS_EZ);
//  scatteredWave(Ezy, EPS_EZ);

  calcH();
}

static inline void calcE()
{
  int i,j;
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ezx[ind(i,j)] = CEZX(i,j)*EZX(i,j) + CEZXLX(i,j)*(HY(i,j) - HY(i-1,j));

  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ezy[ind(i,j)] = CEZY(i,j)*EZY(i,j) - CEZYLY(i,j)*(HX(i,j) - HX(i,j-1));
  
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ez[ind(i,j)] = EZX(i,j) + EZY(i,j);
}

static inline void calcH()
{
  int i,j;
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Hx[ind(i,j)] = CHX(i,j)*HX(i,j) - CHXLY(i,j)*( EZX(i,j+1)-EZX(i,j) + EZY(i,j+1)-EZY(i,j));

  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Hy[ind(i,j)] = CHY(i,j)*HY(i,j) + CHYLX(i,j)*( EZX(i+1,j)-EZX(i,j) + EZY(i+1,j)-EZY(i,j) );  
}

static inline void ntffFrequency()
{
  double cx = N_PX/2;
  double cy = N_PY/2;

  double k_s = field_getK();

  double r0 = 1.0e6;
  double complex coef = csqrt( I*k_s/(8*M_PI*r0) ) * cexp(I*k_s*r0);
  
  ntffInfo nInfo = field_getNTFFInfo();  
  int tp = nInfo.top;     //上面
  int bm = nInfo.bottom; //下面
  int rt = nInfo.right;  //右
  int lt = nInfo.left;	 //左

  const int max_angle = 360;
  double complex ntffEz[max_angle];
  /*
  double complex debugNz[4][max_angle];
  double complex debugLx[2][max_angle];
  double complex debugLy[2][max_angle];
  */
  int ang;
  for(ang=0; ang<max_angle; ang++)
  {
    double rad = ang*M_PI/180.0;

    double rx  = cos(rad), ry = sin(rad);
    double r2x, r2y;

    double complex Nz = 0;
    double complex Lx = 0;
    double complex Ly = 0;
    double complex C_EZ, C_HX, C_HY;
    int i,j;	
    // (left,bottom) -> (right,bottom)
    // 法線ベクトルはn=(0, -1)
    for ( i=lt; i<rt; i++ )
    {
      r2x  =  i-cx;
      r2y  = bm-cy;
      C_EZ = EZ(i,bm);
      C_HX = 0.5 * ( HX(i,bm) + HX(i,bm-1) );
      
      double innerProd = rx*r2x + ry*r2y;  //内積
      Nz  += C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0   ,0,C_HX)
      Lx  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (C_EZ,0,   0)
    }
    /*
    debugNz[0][ang] = Nz;
    debugLx[0][ang] = Lx;
    */

    // (right,bottom) -> (right,top) n=(1,0)
    for ( j=bm; j<tp; j++ )
    {
      r2x  = rt-cx;
      r2y  =  j-cy;
      C_EZ = EZ(rt,j);
      C_HY = 0.5 * ( HY(rt,j) + HY(rt-1,j) );

      double innerProd = rx*r2x + ry*r2y;  //内積
      Nz  += C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,   0, C_HY)
      Ly  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0,C_EZ,    0)
    }
    /*
    debugNz[1][ang] = Nz - debugNz[0][ang];
    debugLy[0][ang] = Ly;
    */

    // (right,top) -> (left,top)  n=(0,1)
    for ( i=lt; i<rt; i++ )
    {
      r2x  =  i-cx;
      r2y  = tp-cy;
      C_EZ = EZ(i,tp);
      C_HX = 0.5 * ( HX(i,tp) + HX(i,tp-1) );

      double innerProd = rx*r2x  + ry*r2y;  //内積
      Nz   -= C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0,    0, -C_HX)
      Lx   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (-C_EZ,0,     0)
    }
    /*
    debugNz[2][ang] = Nz - debugNz[0][ang] - debugNz[1][ang];
    debugLx[1][ang] = Lx - debugLx[0][ang];
    */
    // (left,top) -> (left,bottom)
    for ( j=bm; j<tp; j++ )
    {
      r2x =   lt-cx; r2y = j-cy;
      C_EZ  = EZ(lt,j);
      C_HY  = 0.5 * ( HY(lt,j) + HY(lt-1,j));

      double innerProd = rx*r2x  + ry*r2y;  //内積      
      Nz   -= C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,     0, -C_HY)		
      Ly   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0, -C_EZ,     0)
    }
    /*
    debugNz[3][ang] = Nz - debugNz[0][ang] - debugNz[1][ang] - debugNz[2][ang];
    debugLy[1][ang] = Ly - debugLy[0][ang];
    */
    double complex Lphi = -Lx*sin(rad) + Ly*cos(rad); //極座標変換
    ntffEz[ang] = coef * ( Z_0_S*Nz + Lphi );
  }

  FILE *fpR = fopen("TMntffRe.txt", "w");
  FILE *fpI = fopen("TMntffIm.txt", "w");
  for(ang = 0; ang<max_angle; ang++)
  {
    fprintf(fpR, "%.18lf\n", creal(ntffEz[ang]));
    fprintf(fpI, "%.18lf\n", cimag(ntffEz[ang]));
  }
/*
  FILE *debugNzFp[4];
  debugNzFp[0] = fopen("debugNz0.txt", "w");
  debugNzFp[1] = fopen("debugNz1.txt", "w");
  debugNzFp[2] = fopen("debugNz2.txt", "w");
  debugNzFp[3] = fopen("debugNz3.txt", "w");

  FILE *debugLxFp[2];
  debugLxFp[0] = fopen("debugLx0.txt", "w");
  debugLxFp[1] = fopen("debugLx1.txt", "w");
  
  FILE *debugLyFp[2];
  debugLyFp[0] = fopen("debugLy0.txt", "w");
  debugLyFp[1] = fopen("debugLy1.txt", "w");
  
  for(ang=0; ang<max_angle; ang++)
  {
    for(int j=0; j<4; j++)
      fprintf(debugNzFp[j], "%lf %lf\n", creal(debugNz[j][ang]), cimag(debugNz[j][ang]) );

    for(int j=0; j<2; j++)
    {
      fprintf(debugLxFp[j], "%lf %lf\n", creal(debugLx[j][ang]), cimag(debugLx[j][ang]) );
      fprintf(debugLyFp[j], "%lf %lf\n", creal(debugLy[j][ang]), cimag(debugLy[j][ang]) );
    }
  }
*/
}
