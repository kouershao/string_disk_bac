#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "disk.h"
#include "zernike.h"
#include "initial.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "point.h"

/*********************************************************************/

int main(int argc,char *argv[]) {
  int i,j,k,n,i1,i2,i3,ix;
  double st_Rad,st_t,st_eta,ed_Rad,ed_t,ed_eta;
  double x[5];
  double eg[3];
  double vec[3][3];
  double beta;
  char fname[200];
  FILE *fp;
  int grid_number;
  int ii0 = 10, kk0 = 15;
  int ii1 = 2, kk1 = 2;
  int ii=100, kk=100;
  int counter_i;

  zernike_init(64,64,1000,1000);
  double *Anm = new double[5*Basis]();
  double *Qik = new double[5*Point]();
  
  eta = 1e5;
  Rad = 2;
  landau_t = -0.1;
  
  iput(Anm,landau_t,Rad,eta,N,M,argv[1],1,0);
  for (i=0;i<5;i++)
    calc_fik(Anm + i*Basis,Qik + i*Point);
  for (i=0;i<5*Point;i++)
    Qik[i] = Qik[i]/Qscale;

  sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_I_%d_K_%d_eta_%.1e_%s_vector.vtk",DIR,func_type,boundary,landau_t,Rad,I,K,eta,argv[1]);
  fp = fopen(fname,"w");

  fprintf(fp,"%s\n", "# vtk DataFile Version 3.0");
  fprintf(fp,"%s\n", "Tensor");
  fprintf(fp,"%s\n", "ASCII");
  fprintf(fp,"%s\n", "DATASET UNSTRUCTURED_GRID");

  grid_number = 0;
  //ii = ii0;
  //kk = kk0;
  //counter_i = ii0;
  //for (i=0;i<=I;i++) {
  //  if(counter_i == 0) {
  //    if (ii > ii1) ii--;
  //    counter_i = ii;
  //  }
  //  counter_i--;
    //if(i%ii == 0 && kk > kk1) kk--;
  //  printf("%d - %d\n", i, ii);

  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      if ( i!=0 && i%ii == 0 && k%kk == 0 )
	//if ( counter_i == 0 && k%kk == 0 )
	//printf("%d %d $d \n", i, k);
	//getchar();
	grid_number = grid_number + 1;
    }
  }
  //getchar();

  fprintf(fp,"%s%d%s\n", "POINTS ", grid_number, " float");
  ix = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      //if ( i!=0 && i%ii == 0 && k%kk == 0 ){
      if ( i!=0 && i%ii == 0 && k%kk == 0 ){
	fprintf(fp,"%16.15f %16.15f %16.15f\n", radius[i]*cos(phi[k]), radius[i]*sin(phi[k]), 0.0); 
	printf("%d %f %f\n", grid_number,  radius[i]*cos(phi[k]), radius[i]*sin(phi[k]));}
      ++ix;
    }
    //if(i!=0 && i%ii == 0 && ii > ii1) ii--;
    //if(i!=0 && i%ii == 0 && kk > kk1) kk--;
  }


  fprintf(fp,"\n%s%d\n\n", "POINT_DATA ", grid_number);

  fprintf(fp,"%s\n", "SCALARS beta float 1");
  fprintf(fp,"%s\n", "LOOKUP_TABLE default");
  ix = 0;
  //ii = ii0;
  //kk = kk0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      //if ( counter_i == 0 && k%kk == 0 ){
      if ( i!=0 && i%ii == 0 && k%kk == 0 ){
	QRforEig(x,eg,vec);
	sort(eg,i1,i2,i3);
	beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
	fprintf(fp,"%.6e\n",beta);
      }
      ++ix;
    }
    //if(i!=0 && i%ii == 0 && ii > ii1) ii--;
    //if(i!=0 && i%ii == 0 && kk > kk1) kk--;
  }

  fprintf(fp,"%s\n", "TENSORS Q float");
  ix = 0;
  //ii = ii0;
  //kk = kk0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      //if ( counter_i == 0 && k%kk == 0 ){
	      if ( i!=0 && i%ii == 0 && k%kk == 0 ){
	fprintf(fp,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", x[0]+1.0/3.0,x[1],x[2],x[1],x[3]+1.0/3.0,x[4],x[2],x[4],-x[0]-x[3]+1.0/3.0);
      }
      ++ix;
    }
    //if(i!=0 && i%ii == 0 && ii > ii1) ii--;
    //if(i!=0 && i%ii == 0 && kk > kk1) kk--;
  }

  fprintf(fp,"%s\n", "Vectors v float");
  ix = 0;
  //ii = ii0;
  //kk = kk0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      if ( i!=0 && i%ii == 0 && k%kk == 0 ){
	QRforEig(x,eg,vec);
	sort(eg,i1,i2,i3);
	fprintf(fp,"%16.15lf %16.15lf %16.15lf\n", vec[0][i3],vec[1][i3],vec[2][i3]);
      }
      ++ix;
    }
    //if(i!=0 && i%ii == 0 && ii > ii1) ii--;
    //if(i!=0 && i%ii == 0 && kk > kk1) kk--;
  }


  fclose(fp);
  
  zer_destroy();
  delete [] Anm;
  delete [] Qik;
  
  return 0;
}
