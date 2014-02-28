#include <math.h>
#include "circleModel.h"
#include "field.h"
static double radius;
static double epsilon;
static double posx;
static double posy;

static double eps(double, double, int, int);

double (*circleModel_EPS(double x, double y, double r))(double, double, int , int)
{
  radius = r;
  posx = x;
  posy = y;
  epsilon = 1.6*1.6*EPSILON_0_S;
  return eps;
}

//col : D_Xモード, row : D_Yモード
static double eps(double x, double y, int col, int row)
{
  if(x < N_PML || y < N_PML || x > N_X+N_PML || y > N_Y + N_PML)
    return EPSILON_0_S;

  double dx = x-posx;
  double dy = y-posy;
  //2乗距離
  double len = dx*dx+dy*dy;

  //中心との距離がr+1セル以上なら,そのセルは完全に媒質の外 
  if(len >= (radius+1)*(radius+1))
    return EPSILON_0_S;

  //中心との距離がr-1セル以下なら,そのセルは完全に媒質の外 
  if(len <= (radius-1)*(radius-1))
    return epsilon;

  //さらに32*32分割し媒質内と媒質外の数を求めepsilonを決定する
  double sum=0;
  double i,j;
  for(i=-16+0.5; i<16; i+=1){
    for(j=-16+0.5; j<16; j+=1){
      if(pow(dx+col*i/32.0, 2.0) + pow(dy+row*j/32.0, 2.0) <= radius*radius)
	sum+=1;
    }
  }
  sum /= 32.0*32.0;
  return epsilon*sum + EPSILON_0_S*(1-sum);
}



