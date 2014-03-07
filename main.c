#include <stdio.h>
#include <stdlib.h>
#include "simulator.h"
#include "field.h"

#ifdef USE_OPENGL
#include "drawer.h"
#include <GL/glew.h>
#include <GLUT/glut.h>
void display()
{
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  drawer_paintImage( N_PML, N_PML, N_X+N_PML, N_Y+N_PML, N_PX, N_PY, simulator_getDrawingData());   //テクスチャの更新
  drawer_paintModel( N_PML, N_PML, N_X+N_PML, N_Y+N_PML, N_PX, N_PY, simulator_getEps());
  drawer_draw();
    
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
}

void idle()
{
  simulator_calc();  
  glutPostRedisplay();  //再描画
}
#endif

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
  int waveAngle = 0;
  enum MODEL   modelType = MIE_CYLINDER; // モデルの種類
  enum SOLVER solverType = TE_2D;        // 計算方法
  simulator_init(width, height, h_u, pml, lambda, waveAngle, step, modelType, solverType);    //simulator
    
#ifdef USE_OPENGL
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

#ifndef USE_OPENGL
    while(1)
    {
      simulator_calc(); //only calculate mode
    }    
#endif
  
    return 1;
}
