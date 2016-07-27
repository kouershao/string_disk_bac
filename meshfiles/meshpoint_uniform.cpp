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
  Rad = 5;
  landau_t = -1;
  
  iput(Anm,landau_t,Rad,eta,N,M,argv[1],1,0);
  for (i=0;i<5;i++)
    calc_fik(Anm + i*Basis,Qik + i*Point);
  for (i=0;i<5*Point;i++)
    Qik[i] = Qik[i]/Qscale;

  if(!(fp0 = fopen("visu/disk_mesh_4.vtk","r"))){
    printf("%s\n", "input file open error.");
    return -1;
  }
  //grid_number = ReadMeshPoints(sx, sy, sz);

  sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_%s_beta.vtk",DIR,func_type,boundary,landau_t,Rad,argv[1]);
  fp1 = fopen(fname,"w");

  sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_I_%d_K_%d_eta_%.1e_%s_point.csv",DIR,func_type,boundary,landau_t,Rad,I,K,eta,argv[1]);
  fp2 = fopen(fname,"w");

  fgets(line, 200, fp0);
  fputs(line, fp1);
  fgets(line, 200, fp0);
  fputs(line, fp1);
  fgets(line, 200, fp0);
  fputs(line, fp1);
  fgets(line, 200, fp0);
  fputs(line, fp1);

  fscanf(fp0, "%s %d %s\n", line, &num_point, line2);
  fprintf(fp1, "%s %d %s\n", line, num_point, line2);

  sx = (double*)malloc(sizeof(double)*max_grid_number);
  if (sx == NULL) exit(1);
  sy = (double*)malloc(sizeof(double)*max_grid_number);
  if (sy == NULL) exit(1);
  sz = (double*)malloc(sizeof(double)*max_grid_number);
  if (sz == NULL) exit(1);

  for (i=0; i<num_point; i++){
    fscanf(fp0, "%lf %lf %lf\n", &xx, &yy, &zz);
    fprintf(fp1, "%lf %lf %lf\n", xx, yy, zz);
    sx[i] = xx;
    sy[i] = yy;
    sz[i] = zz;
  }

  while (fgets(line, 200, fp0) != NULL)
    {
      fputs(line, fp1);
    }

  fprintf(fp1,"\n%s%d\n\n", "POINT_DATA ", num_point);

  fprintf(fp1,"%s\n", "SCALARS beta float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", beta);
    fprintf(fp2,"%16.15lf, %16.15lf, %16.15lf\n", xx, yy, beta);
    //printf("%16.15lf\n", beta);
  }

  fprintf(fp1,"%s\n", "SCALARS cl float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", eg[i1] - eg[i2]);
  }

  fprintf(fp1,"%s\n", "SCALARS lambda1 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", eg[i1]);
  }

  fprintf(fp1,"%s\n", "SCALARS lambda2 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", eg[i2]);
  }

  fprintf(fp1,"%s\n", "SCALARS lambda3 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", eg[i3]);
  }

  fprintf(fp1,"%s\n", "SCALARS q1 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", q1);
  }

  fprintf(fp1,"%s\n", "SCALARS trQ2 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", 2*(eg[i1]*eg[i1]+eg[i2]*eg[i2]+eg[i1]*eg[i2]));
  }

  fprintf(fp1,"%s\n", "SCALARS q3 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", q3);
  }

  fprintf(fp1,"%s\n", "SCALARS q5 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf\n", q5);
  }

  fprintf(fp1,"%s\n", "Vectors v1 float");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    xx = sx[gi];
    yy = sy[gi];
    zz = sz[gi];

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
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

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
    fprintf(fp1,"%16.15lf %16.15lf %16.15lf\n", vec[0][i3],vec[1][i3],vec[2][i3]);
  }

  fclose(fp1);
  fclose(fp0);
  fclose(fp2);
  
  zer_destroy();
  delete [] Anm;
  delete [] Qik;
  free(sx);
  free(sy);
  free(sz);

  return 0;
}
