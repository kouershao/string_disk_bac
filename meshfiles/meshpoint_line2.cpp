#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "disk.h"
#include "zernike.h"
#include "initial.h"
#include "rotation.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "point.h"

/*********************************************************************/

int main(int argc,char *argv[]) {
  int num_point;
  int num_cell;
  int i,j,k,n,i1,i2,i3,ix;
  double st_Rad,st_t,st_eta,ed_Rad,ed_t,ed_eta;
  double x[5];
  double q1,q2,q3,q4,q5;
  double tensor[5];
  double eg[3];
  double vec[3][3];
  int gi;
  double r,t,p;
  double beta;
  char fname[200];
  FILE *fp1;
  FILE *fp2;
  FILE *fp0;
  FILE *fp;
  char line[200];
  char line2[200];
  int grid_number;
  int max_grid_number = 1000000;
  double xx, yy, zz;
  double* sx;
  double* sy;
  double* sz;
  
  zernike_init(64,64,1000,1000);
  double *Anm = new double[5*Basis]();
  double *Qik = new double[5*Point]();
  
  eta = 1e2;
  Rad = 2;
  landau_t = -0.1;
  
  iput(Anm,landau_t,Rad,eta,N,M,argv[1],1,0);
  for (i=0;i<5;i++)
    calc_fik(Anm + i*Basis,Qik + i*Point);
  //for (i=0;i<5*Point;i++)
    //Qik[i] = Qik[i]/Qscale;

  //grid_number = ReadMeshPoints(sx, sy, sz);

  sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_I_%d_K_%d_eta_%.1e_%s_line.csv",DIR,func_type,boundary,landau_t,Rad,I,K,eta,argv[1]);
  fp1 = fopen(fname,"w");

  num_point = 1000;

  sx = (double*)malloc(sizeof(double)*max_grid_number);
  if (sx == NULL) exit(1);
  sy = (double*)malloc(sizeof(double)*max_grid_number);
  if (sy == NULL) exit(1);
  sz = (double*)malloc(sizeof(double)*max_grid_number);
  if (sz == NULL) exit(1);

  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = (double(gi) + 0.5)/double(num_point);
    yy = 0.0;
    zz = 0.0;

    r = sqrt(xx*xx + yy*yy + zz*zz);
    if (r <= 1e-14) {
      p = 0;
    }
    else {
      if (fabs(xx) < 1e-14 && yy >= 0)
  	p = PI/2;
      else if (fabs(xx) < 1e-14 && yy < 0)
  	p = PI*3/2;
      else if (xx > 0 && yy >= 0)
  	p = atan(yy/xx);
      else if (xx > 0 && yy < 0)
  	p = 2*PI + atan(yy/xx);
      else if (xx < 0 && yy >= 0)
  	p = PI + atan(yy/xx);
      else if (xx < 0 && yy < 0)
  	p = PI + atan(yy/xx);
    }
    q1 = interpolation(Qik,I,K,r,p,1000);
    q2 = interpolation(Qik + Point,I,K,r,p,1000);
    q3 = interpolation(Qik + 2*Point,I,K,r,p,1000);
    q4 = interpolation(Qik + 3*Point,I,K,r,p,1000);
    q5 = interpolation(Qik + 4*Point,I,K,r,p,1000);

    tensor[0] = q1;
    tensor[1] = q2;
    tensor[2] = q3;
    tensor[3] = q4;
    tensor[4] = q5;

    //    printf("%f %f %f\n", xx, yy, zz);
    //printf("%f %f\n", r, p);
    //getchar();

    
    QRforEig(tensor,eg,vec);
    sort(eg,i1,i2,i3);
    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
    fprintf(fp1,"%16.15lf, %16.15lf, %16.15lf\n", xx, 0.5*(eg[i3]+eg[i2]), 0.5*(eg[i3]-eg[i2]));
    printf("%16.15lf, %16.15lf\n", xx, 0.5*(eg[i1]+eg[i2]));
  }


  fclose(fp1);
  
  zer_destroy();
  delete [] Anm;
  delete [] Qik;
  free(sx);
  free(sy);
  free(sz);

  return 0;
}
