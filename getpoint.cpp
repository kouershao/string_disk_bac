#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "disk.h"
#include "zernike.h"
#include "initial.h"
#include "calculate_Fb.h"

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
  
  zernike_init(64,64,128,100);
  Anm = new double[5*Basis]();
  Qik = new double[5*Point]();
  fbulk = new double[Point]();
  
  st_eta = 1e2;
  ed_eta = 1e2;
  st_Rad = 1;
  ed_Rad = 1;
  st_t = 1;
  ed_t = 1;
  
  for (eta = st_eta;eta <= ed_eta;eta = eta*1.1) {
    for (Rad = st_Rad;Rad <= ed_Rad;Rad += 1) {
      for (landau_t = st_t;landau_t >= ed_t;landau_t -= 1) {
	
	iput(Anm,landau_t,Rad,eta,N,M,argv[1],1,0);
	for (i=0;i<5;i++)
	  calc_fik(Anm + i*Basis,Qik + i*Point);
	
	calc_bulkenergy_landau_de(Qik,fbulk);
	
	// for (i=0;i<5*Point;i++)
	//   Qik[i] = Qik[i]/Qscale;
	
	sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_I_%d_K_%d_eta_%.1e_%s_point.txt",DIR,func_type,boundary,landau_t,Rad,I,K,eta,argv[1]);
	fp = fopen(fname,"w");
	ix = 0;
	for (i=0;i<=I;i++) {
	  for (k=0;k<K;k++) {
	    for (n=0;n<5;n++)
	      x[n] = Qik[n*Point + ix];
	    QRforEig(x,eg,vec);
	    sort(eg,i1,i2,i3);
	    beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
	    fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",radius[i],phi[k],eg[i1],eg[i2],eg[i3],vec[0][i3],vec[1][i3],vec[2][i3],beta,fbulk[ix]);
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
  delete [] fbulk;
  
  return 0;
}
