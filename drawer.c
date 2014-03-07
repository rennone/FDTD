#ifdef USE_OPENGL

#include "drawer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/glew.h>
#include <GLUT/glut.h>
#include "function.h"

typedef struct {
  GLfloat r,g,b;
}colorf;

static const int vertexNum = 4; //頂点数
#define TEX_NX 256
#define TEX_NY 256
static colorf texColor[TEX_NX][TEX_NY];
static GLuint ver_buf, tex_buf;
static GLuint texId;

static double (*colorMode)(double complex);
static void colorTransform(double p, colorf *c);

static GLfloat vertices[] =
  {-1.0f, -1.0f, 0.0f,
   +1.0f, -1.0f, 0.0f, 
   +1.0f, +1.0f, 0.0f, 
   -1.0f, +1.0f, 0.0f};

static GLfloat texCoords[] =
  { 0.0f, 0.0f,
    0.0f, 1.0f,
    1.0f, 1.0f,
    1.0f, 0.0f };

//--------------------prototype--------------------//
void drawer_paintImage(int l, int b,int r, int t, int wid,int hei, double complex*,...);
void drawer_draw();
//--------------------------------------//


//--------------public Method-----------------//
void (*drawer_getDraw(void))(void)
{
  return drawer_draw;
}
//--------------------------------------//

void drawer_init(enum COLOR_MODE cm)
{
  int i,j;
  for(i=0; i<TEX_NX; i++)
    for(j=0; j<TEX_NY; j++)
    {
      texColor[i][j].r = 0; texColor[i][j].g = 0;  texColor[i][j].b = 0;
    }
  
  if(cm == CREAL)
    colorMode = creal;
  else
    colorMode = cnorm;//cabs;

  glGenBuffers(1, &ver_buf);

  glBindBuffer(GL_ARRAY_BUFFER, ver_buf);
  glBufferData(GL_ARRAY_BUFFER, 3*vertexNum*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

  glGenBuffers(1, &tex_buf);
  glBindBuffer(GL_ARRAY_BUFFER, tex_buf);
  glBufferData(GL_ARRAY_BUFFER, 2*vertexNum*sizeof(GLfloat), texCoords, GL_STATIC_DRAW);

  //
  glEnable( GL_TEXTURE_2D );
  glGenTextures( 1, &texId );

  glActiveTexture( GL_TEXTURE0 );

  glBindTexture( GL_TEXTURE_2D, texId );
  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, TEX_NX, TEX_NY, 0, GL_RGB, GL_FLOAT, texColor);

    //min, maxフィルタ
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );    //min, maxフィルター
}

void drawer_draw()
{  
  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  //glEnable(GL_ALPHA_TEST);
  
  glBindBuffer( GL_ARRAY_BUFFER, ver_buf);
  glVertexPointer( 3, GL_FLOAT, 0, 0);

  glBindBuffer( GL_ARRAY_BUFFER, tex_buf);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);
  
  glBindTexture( GL_TEXTURE_2D, texId );
  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, TEX_NX, TEX_NY, 0, GL_RGB, GL_FLOAT, texColor);  
  glDrawArrays( GL_POLYGON, 0, vertexNum);
  
}

//todo 可変長引数を利用して, 複数のデータの平均で色を出すようにするべき?
void drawer_paintImage(int left, int bottom, int right, int top, int width, int height, double complex *phis, ...)
{
  colorf c;
  double complex cphi;
  double ux = 1.0*(right-left)/TEX_NX;
  double uy = 1.0*(top-bottom)/TEX_NY;
  double u = max(ux,uy);
  int i,j;
  double x,y;
  
  for(i=0,x=left; i<TEX_NX && x<right; i++, x+=u){
    for(j=0,y=bottom; j<TEX_NY && y<top; j++, y+=u){
      cphi = cbilinear(phis,x,y,width,height);
      colorTransform(colorMode(cphi), &c);
      texColor[i][j] = c;
    }
  }
}

//todo 可変長引数を利用して, 複数のデータの平均で色を出すようにするべき?
void drawer_paintModel(int left, int bottom, int right, int top, int width, int height, double *phis, ...)
{
  colorf c;
  double dphi;
  double ux = 1.0*(right-left)/TEX_NX;
  double uy = 1.0*(top-bottom)/TEX_NY;
  double u = max(ux,uy);
  int i,j;
  double x,y;
  
  for(i=0,x=left; i<TEX_NX && x<right; i++, x+=u){
    for(j=0,y=bottom; j<TEX_NY && y<top; j++, y+=u){
      dphi = dbilinear(phis, x,y,width,height);
      double n = 1-1.0/dphi;
      texColor[i][j].r -= n;
      texColor[i][j].g -= n;
      texColor[i][j].b -= n;
    }
  }
}

void drawer_finish()
{
  printf("drawer_finish not implemented\n");
}

//--------------------public Method--------------------//




//--------------------Color Trancform---------------------//
static void colorTransform(double phi, colorf *col)
{
  double range = 2.0; //波の振幅  
  double ab_phi = phi < 0 ? -phi : phi;
  
  double a = ab_phi < range ? (ab_phi <  range/3.0 ? 3.0/range*ab_phi : (-3.0/4.0/range*ab_phi+1.25) ) : 0.5;
  
  col->r = phi > 0 ? a:0;
  col->b = phi < 0 ? a:0;
  col->g = min(1.0, max(0.0, -3*ab_phi+2));
}

#endif //USE_OPENGL
