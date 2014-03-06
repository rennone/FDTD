#include "drawer.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef _USE_OPENGL
#include <GL/glew.h>
#include <GLUT/glut.h>
#endif

#include "simulator.h"

#include "field.h"

void display(void);
void idle(void);

void display(void)
{
  #ifdef _USE_OPENGL
  //todo gl method should be remove here, go to drawer.c
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  drawer_paintImage( N_PML, N_PML, N_X+N_PML, N_Y+N_PML, N_PX, N_PY, simulator_getDrawingData());   //テクスチャの更新
  drawer_paintModel( N_PML, N_PML, N_X+N_PML, N_Y+N_PML, N_PX, N_PY, simulator_getEps());
  drawer_draw();
    
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
  
  #endif
}

void idle(void)
{
  simulator_calc();
  
  #ifdef _USE_OPENGL
  glutPostRedisplay();  //再描画
  #endif
}

const int windowX = 100;
const int windowY = 100;
const int windowWidth = 500;
const int windowHeight = 500;

int main( int argc, char *argv[] )
{
  int    width  = 256; //横幅(nm)
  int    height = 256; //縦幅(nm)
  double   h_u  = 1;   //1セルの大きさ(nm)
  int       pml = 10;  //pmlレイヤの数
  double lambda = 60;  //波長(nm)
  int      step = 2000; //計算ステップ 
  enum MODEL   modelType = MIE_CYLINDER; // モデルの種類
  enum SOLVER solverType = TE_2D;        // 計算方法
  simulator_init(width, height, h_u, pml, lambda, step, modelType, solverType);    //simulator
    
#ifdef _USE_OPENGL  
    glutInit(&argc, argv);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("FDTD Simulator");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glewInit();


    drawer_init(CABS);  //drawerの初期化      /* CREAL : 振幅表示, CABS : 強度表示 */
    glutMainLoop();
#endif

#ifndef _USE_OPENGL
    //only calculate mode  
    while(1)
    {
      simulator_calc(); 
    }    
#endif
  
    return 1;
}