#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "sphere.h"
#include "zernike.h"
#include "initial.h"

//#include "calcbeta.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "calc_vector.h"
#include "point.h"
#include "rotation.h"
#include "check.h"

int main(int argc,char *argv[]) {
  int num_point;
  int num_cell;
  char fname[200];
  FILE *fp1;
  FILE *fp0;
  int gi;
  int i;
  double x, y, z;
  double r,t,p;
  double q1,q2,q3,q4,q5;
  double tensor[5];
  double eg[3];
  double vec[3][3];
  double beta;
  char line[200];
  char line2[200];
  int n,j,k,i1,i2,i3;
  double* sx;
  double* sy;
  double* sz;
  int max_grid_number = 1000000;

  sphere_param_init(128,32,32,4,200,201,200);
  zernike_init();
  malloc_variables_init();
  iput(Anlm,landau_t,Rad,eta,32,32,4,argv[1]);
  for (i=0;i<5;i++)
    calc_fijk(Anlm + i*Basis,Qijk + i*Point);

  if(!(fp0 = fopen("visu/ball_mesh_4.vtk","r"))){
    printf("%s\n", "input file open error.");
    return -1;
  }
  //grid_number = ReadMeshPoints(sx, sy, sz);

  sprintf(fname,"visu/output/beta.vtk");
  fp1 = fopen(fname,"w");  

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
    fscanf(fp0, "%lf %lf %lf\n", &x, &y, &z);
    fprintf(fp1, "%lf %lf %lf\n", x, y, z);
    sx[i] = x;
    sy[i] = y;
    sz[i] = z;
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

    x = sx[gi];
    y = sy[gi];
    z = sz[gi];
    r = sqrt(x*x + y*y + z*z);
    if (r <= 1e-14) {
      t = 0;
      p = 0;
    }
    else {
      t = acos(z/r);
      if (fabs(x) < 1e-14 && y >= 0)
  	p = PI/2;
      else if (fabs(x) < 1e-14 && y < 0)
  	p = PI*3/2;
      else if (x > 0 && y >= 0)
  	p = atan(y/x);
      else if (x > 0 && y < 0)
  	p = 2*PI + atan(y/x);
      else if (x < 0 && y >= 0)
  	p = PI + atan(y/x);
      else if (x < 0 && y < 0)
  	p = PI + atan(y/x);
    }
    q1 = interpolation(Qijk,I,J,K,r,t,p,200);
    q2 = interpolation(Qijk + Point,I,J,K,r,t,p,200);
    q3 = interpolation(Qijk + 2*Point,I,J,K,r,t,p,200);
    q4 = interpolation(Qijk + 3*Point,I,J,K,r,t,p,200);
    q5 = interpolation(Qijk + 4*Point,I,J,K,r,t,p,200);
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
    
    QRforEig(tensor,eg,vec);
    sort(eg,i1,i2,i3);
    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
    fprintf(fp1,"%16.15lf\n", beta);
  }

  fprintf(fp1,"%s\n", "SCALARS cl float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    x = sx[gi];
    y = sy[gi];
    z = sz[gi];
    r = sqrt(x*x + y*y + z*z);
    if (r <= 1e-14) {
      t = 0;
      p = 0;
    }
    else {
      t = acos(z/r);
      if (fabs(x) < 1e-14 && y >= 0)
  	p = PI/2;
      else if (fabs(x) < 1e-14 && y < 0)
  	p = PI*3/2;
      else if (x > 0 && y >= 0)
  	p = atan(y/x);
      else if (x > 0 && y < 0)
  	p = 2*PI + atan(y/x);
      else if (x < 0 && y >= 0)
  	p = PI + atan(y/x);
      else if (x < 0 && y < 0)
  	p = PI + atan(y/x);
    }
    q1 = interpolation(Qijk,I,J,K,r,t,p,200);
    q2 = interpolation(Qijk + Point,I,J,K,r,t,p,200);
    q3 = interpolation(Qijk + 2*Point,I,J,K,r,t,p,200);
    q4 = interpolation(Qijk + 3*Point,I,J,K,r,t,p,200);
    q5 = interpolation(Qijk + 4*Point,I,J,K,r,t,p,200);
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
    
    QRforEig(tensor,eg,vec);
    sort(eg,i1,i2,i3);
    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
    fprintf(fp1,"%16.15lf\n", eg[i3]-eg[i2]);
  }


  fprintf(fp1,"%s\n", "SCALARS lambda1 float 1");
  fprintf(fp1,"%s\n", "LOOKUP_TABLE default");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    x = sx[gi];
    y = sy[gi];
    z = sz[gi];
    r = sqrt(x*x + y*y + z*z);
    if (r <= 1e-14) {
      t = 0;
      p = 0;
    }
    else {
      t = acos(z/r);
      if (fabs(x) < 1e-14 && y >= 0)
  	p = PI/2;
      else if (fabs(x) < 1e-14 && y < 0)
  	p = PI*3/2;
      else if (x > 0 && y >= 0)
  	p = atan(y/x);
      else if (x > 0 && y < 0)
  	p = 2*PI + atan(y/x);
      else if (x < 0 && y >= 0)
  	p = PI + atan(y/x);
      else if (x < 0 && y < 0)
  	p = PI + atan(y/x);
    }
    q1 = interpolation(Qijk,I,J,K,r,t,p,200);
    q2 = interpolation(Qijk + Point,I,J,K,r,t,p,200);
    q3 = interpolation(Qijk + 2*Point,I,J,K,r,t,p,200);
    q4 = interpolation(Qijk + 3*Point,I,J,K,r,t,p,200);
    q5 = interpolation(Qijk + 4*Point,I,J,K,r,t,p,200);
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
    
    QRforEig(tensor,eg,vec);
    sort(eg,i1,i2,i3);
    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
    fprintf(fp1,"%16.15lf\n", eg[i3]);
  }

  fprintf(fp1,"%s\n", "VECTORS v1 float");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    x = sx[gi];
    y = sy[gi];
    z = sz[gi];
    r = sqrt(x*x + y*y + z*z);
    if (r <= 1e-14) {
      t = 0;
      p = 0;
    }
    else {
      t = acos(z/r);
      if (fabs(x) < 1e-14 && y >= 0)
  	p = PI/2;
      else if (fabs(x) < 1e-14 && y < 0)
  	p = PI*3/2;
      else if (x > 0 && y >= 0)
  	p = atan(y/x);
      else if (x > 0 && y < 0)
  	p = 2*PI + atan(y/x);
      else if (x < 0 && y >= 0)
  	p = PI + atan(y/x);
      else if (x < 0 && y < 0)
  	p = PI + atan(y/x);
    }
    q1 = interpolation(Qijk,I,J,K,r,t,p,200);
    q2 = interpolation(Qijk + Point,I,J,K,r,t,p,200);
    q3 = interpolation(Qijk + 2*Point,I,J,K,r,t,p,200);
    q4 = interpolation(Qijk + 3*Point,I,J,K,r,t,p,200);
    q5 = interpolation(Qijk + 4*Point,I,J,K,r,t,p,200);
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
    
    QRforEig(tensor,eg,vec);
    sort(eg,i1,i2,i3);
    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
    fprintf(fp1,"%16.15lf %16.15lf %16.15lf\n", vec[0][i3],vec[1][i3],vec[2][i3]);
  }


  fprintf(fp1,"%s\n", "TENSORS Q float");
  for ( gi = 0; gi<num_point; gi++ ){
    // Replace the q11-q33 by elements of Q at points (sx[gi], sy[gi], sz[gi])
    // fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q11, q12, q13, q21, q22, q23, q31, q32, q33);

    x = sx[gi];
    y = sy[gi];
    z = sz[gi];
    r = sqrt(x*x + y*y + z*z);
    if (r <= 1e-14) {
      t = 0;
      p = 0;
    }
    else {
      t = acos(z/r);
      if (fabs(x) < 1e-14 && y >= 0)
  	p = PI/2;
      else if (fabs(x) < 1e-14 && y < 0)
  	p = PI*3/2;
      else if (x > 0 && y >= 0)
  	p = atan(y/x);
      else if (x > 0 && y < 0)
  	p = 2*PI + atan(y/x);
      else if (x < 0 && y >= 0)
  	p = PI + atan(y/x);
      else if (x < 0 && y < 0)
  	p = PI + atan(y/x);
    }
    q1 = interpolation(Qijk,I,J,K,r,t,p,200);
    q2 = interpolation(Qijk + Point,I,J,K,r,t,p,200);
    q3 = interpolation(Qijk + 2*Point,I,J,K,r,t,p,200);
    q4 = interpolation(Qijk + 3*Point,I,J,K,r,t,p,200);
    q5 = interpolation(Qijk + 4*Point,I,J,K,r,t,p,200);
    q1 = q1/Qscale;
    q2 = q2/Qscale;
    q3 = q3/Qscale;
    q4 = q4/Qscale;
    q5 = q5/Qscale;

    fprintf(fp1,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", q1+1.0/3.0,q2,q3,q2,q4+1.0/3.0,q5,q3,q5,-q1-q4+1.0/3.0);
  }
  
  fclose(fp0);
  fclose(fp1);

  zer_destroy();
  var_destroy();
  free(sx);
  free(sy);
  free(sz);

  return 0;
}
