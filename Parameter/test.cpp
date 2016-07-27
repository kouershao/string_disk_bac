#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PI 3.141592653589793

const int I = 200;
const int J = 200;


int main () {
  FILE *fp;
  char path[100];
  int i,j,k,n,l,m,s;
  int read_ok;
  double x,y;

  double radius[I+1];
  double coe_r[I];

  double lambda[J];
  double coe_theta[J];

  fp = fopen("roots_I_200.txt","w");
  for (i=0;i<=I;i++) {
    radius[i] = 1.0*i/I;
    fprintf(fp,"%.16e\n",radius[i]);
  }
  fclose(fp);
  
  fp = fopen("roots_J_201.txt","w");
  for (i=0;i<=J;i++) {
    lambda[i] = cos(PI*i/J);
    fprintf(fp,"%.16e\n",lambda[i]);
  }
  fclose(fp);

  fp = fopen("coe_r_200.txt","w");
  for (i=0;i<I;i++) {
    fprintf(fp,"%.16e\n",1.0/I);
  }
  fclose(fp);
  
  fp = fopen("coe_theta_201.txt","w");
  for (i=0;i<J;i++) {
    fprintf(fp,"%.16e\n",PI/J);
  }
  fclose(fp);


  return 0;
}
