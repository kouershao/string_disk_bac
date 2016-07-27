#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <thread>
//#include <future>
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
  int ii=4, kk=4;
  
  zernike_init(64,64,128,200);
  double *Anm = new double[5*Basis]();
  double *Qik = new double[5*Point]();
  
   st_eta = 1e5;
  ed_eta = 1e5;
  st_Rad = 11;
  ed_Rad = 11;
  st_t = 1;
  ed_t = 1 ; 


  for (eta = st_eta;eta <= ed_eta;eta = eta*1.1) {

    for (Rad = st_Rad;Rad <= ed_Rad;Rad += 2) {
   
      for (landau_t = st_t;landau_t >= ed_t;landau_t -= 1) {
 //eta = 1e5;
  //Rad = 2;
  //landau_t = -1;
  
  iput(Anm,landau_t,Rad,eta,N,M,argv[1],1,0);
  for (i=0;i<5;i++)
    calc_fik(Anm + i*Basis,Qik + i*Point);
  for (i=0;i<5*Point;i++)
    Qik[i] = Qik[i]/Qscale;

  sprintf(fname,"%s/Drawpoint/%s%s_R_%.2f_t_%.2f_%s_point.vtk",DIR,func_type,boundary,Rad,landau_t,argv[1]);
  fp = fopen(fname,"w");
  printf("%s\n", fname);

  fprintf(fp,"%s\n", "# vtk DataFile Version 3.0");
  fprintf(fp,"%s\n", "Tensor");
  fprintf(fp,"%s\n", "ASCII");
  fprintf(fp,"%s\n", "DATASET UNSTRUCTURED_GRID");

  grid_number = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      if ( i%ii == 0 && k%kk == 0 )
	//printf("%d %d $d \n", i, k);
	//getchar();
	grid_number = grid_number + 1;
    }
  }

  fprintf(fp,"%s%d%s\n", "POINTS ", grid_number, " float");
  ix = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      if ( i%ii == 0 && k%kk == 0 ){
      fprintf(fp,"%16.15f %16.15f %16.15f\n", radius[i]*cos(phi[k]), radius[i]*sin(phi[k]), 0.0); 
      printf("%d %f %f\n", grid_number,  radius[i]*cos(phi[k]), radius[i]*sin(phi[k]));}
      ++ix;
    }
  }


  fprintf(fp,"\n%s%d\n\n", "POINT_DATA ", grid_number);

  fprintf(fp,"%s\n", "SCALARS beta float 1");
  fprintf(fp,"%s\n", "LOOKUP_TABLE default");
  ix = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      if ( i%ii == 0 && k%kk == 0 ){
	    QRforEig(x,eg,vec);
	    sort(eg,i1,i2,i3);
	    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
	    fprintf(fp,"%.6e\n",beta);
    }
      ++ix;
    }
  }

  fprintf(fp,"%s\n", "TENSORS Q float");
  ix = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      if ( i%ii == 0 && k%kk == 0 ){
	fprintf(fp,"%16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f \n %16.15f %16.15f %16.15f\n\n", x[0]+1.0/3.0,x[1],x[2],x[1],x[3]+1.0/3.0,x[4],x[2],x[4],-x[0]-x[3]+1.0/3.0);
    }
      ++ix;
    }
  }

  fprintf(fp,"%s\n", "Vectors v float");
  ix = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	x[n] = Qik[n*Point + ix];
      if ( i%ii == 0 && k%kk == 0 ){
	QRforEig(x,eg,vec);
	sort(eg,i1,i2,i3);
	fprintf(fp,"%16.15lf %16.15lf %16.15lf\n", vec[0][i3],vec[1][i3],vec[2][i3]);
    }
      ++ix;
    }
  }


  fclose(fp);

  
      }

    }
 
  }


  zer_destroy();
  delete [] Anm;
  delete [] Qik;
  
  return 0;
}


//int main (int argc, char *argv[])
//{

//	if(argc != 2) {
//		fprintf(stderr, "Usage: ./point four_<num>");
//		exit(1);
//	}
//	mm(argv[1]);
	

//	return 0;
//}

