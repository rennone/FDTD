#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulator.h"
#include "field.h"
#include "fdtdTE.h"
#include "fdtdTE_upml.h"
#include "fdtdTM.h"
#include "fdtdTM_upml.h"
#include "drawer.h"
#include "models.h"

//#define _OS_UNIX_

#ifdef _OS_UNIX_
#include <sys/time.h>
static struct timeval timer1, timer2;
#endif

static void (*update)() = NULL;        //update method のポインタ
static void (* finishMethod)() = NULL;
static void (* initMethod)() = NULL;

//データ取得用のポインタ
static double complex* (*getDataX)() = NULL;
static double complex* (*getDataY)() = NULL;
static double complex* (*getDataZ)() = NULL;

//上記の3つのどれかが入る
static double complex* (* getDrawDataMethod)() = NULL;

//εデータ取得用のポインタ
static double* (* getEpsMethod)() = NULL;

static char folderName[256];


static void finish(void);

void setTE(){
  update = fdtdTE_getUpdate();
  initMethod = fdtdTE_getInit();
  finishMethod = fdtdTE_getFinish();
  getDataX = fdtdTE_getEx;
  getDataY = fdtdTE_getEy;
  getDataZ = fdtdTE_getHz;

  getEpsMethod = fdtdTE_getEps;
  getDrawDataMethod = getDataY;
  strcpy(folderName, "TE/");
}

void setTM(){
  update = fdtdTM_getUpdate();
  initMethod = fdtdTM_getInit();
  finishMethod = fdtdTM_getFinish();
  getDataX = fdtdTM_getHx;
  getDataY = fdtdTM_getHy;
  getDataZ = fdtdTM_getEz;

  getEpsMethod = fdtdTM_getEps;
  getDrawDataMethod = getDataZ;
  strcpy(folderName, "TM/");
}

void setTMupml(){
  update = fdtdTM_upml_getUpdate();
  initMethod = fdtdTM_upml_getInit();
  finishMethod = fdtdTM_upml_getFinish();
  getDataX = fdtdTM_upml_getHx;
  getDataY = fdtdTM_upml_getHy;
  getDataZ = fdtdTM_upml_getEz;

  getEpsMethod = fdtdTM_upml_getEps;
  getDrawDataMethod = getDataZ;
  strcpy(folderName, "TMupml/");  
}

void setTEupml(){
  update = fdtdTE_upml_getUpdate();
  initMethod = fdtdTE_upml_getInit();
  finishMethod = fdtdTE_upml_getFinish();
  getDataX = fdtdTE_upml_getEx;
  getDataY = fdtdTE_upml_getEy;
  getDataZ = fdtdTE_upml_getHz;

  getEpsMethod = fdtdTE_upml_getEps;
  getDrawDataMethod = getDataY;
  strcpy(folderName, "TEupml/");  
}

void setSolver(enum SOLVER solver)
{
  switch(solver){
  case TE_2D :
    setTE();
    break;
  case TM_2D :
    setTM();
    break;
  case TM_UPML_2D:
    setTMupml();
    break;
  case TE_UPML_2D:
    setTEupml();
    break;
  default :
    break;
  }
  
  (*initMethod)(); //Solverの初期化, EPS, Coeffの設定  
}

void simulator_calc()
{
  (*update)();
  //時間を一つ進める
  
  if(field_nextStep())
  {
    finish();
  }
}

void simulator_init(int width, int height , double h_u, int pml, double lambda, int step,  enum MODEL modelType, enum SOLVER solverType)
{
  //横幅(nm), 縦幅(nm), 1セルのサイズ(nm), pmlレイヤの数, 波長(nm), 計算ステップ
  setField(width, height, h_u, pml, lambda, step); //必ず最初にこれ

  /*NO_MODEL. MIE_CYLINDER, SHELF, NONSHELF*/
  setModel(modelType);  //次にこれ,  モデル(散乱体)を定義

  setSolver(solverType);//Solverの設定と初期化

  #ifdef _OS_UNIX_
  gettimeofday(&timer1, NULL); //開始時間の取得
  #endif
}

static void finish(){
  printf("finish\n");

#ifdef _OS_UNIX_
  gettimeofday(&timer2,NULL);
  printf("time = %lf \n", timer2.tv_sec-timer1.tv_sec+(timer2.tv_usec-timer1.tv_usec)*1e-6);
  #endif
  
  //char fileName[256];
  //strcpy(fileName, folderName);
  //strcat(fileName, "mieTE.txt");
  //field_outputElliptic(fileName, (*getDataY)(), N_PX/2.0, N_PY/2.0, 1.2*field_getLambda());
  //field_outputAllData("./data.txt", (*getDataZ)());
  
  (*finishMethod)(); //メモリの解放等
  exit(0);
}

double complex* simulator_getDrawingData(void){
  return (* getDrawDataMethod)();
}

double * simulator_getEps(void)
{
  return (* getEpsMethod)();
}


