#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.141592653589793

double S0;

int main () {
  double landau_t,Rad,Energy,ep;
  double l1,l2,l3,f0;

  
  FILE *fp,*fp2;
  fp  = fopen("energy_half.txt","r");
  fp2 = fopen("energy_half_ep.txt","w");
  
  for (int i = 0;i < 46;i++) {
    fscanf(fp,"%lf %lf %lf\n",&landau_t,&Rad,&Energy);
    ep = 1/2.0/Rad/Rad;

    S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;
    l1 = +2.0/3*S0;
    l2 = -1.0/3*S0;
    l3 = -1.0/3*S0;
    f0 = landau_t/2*(l1*l1 + l2*l2 + l3*l3) - sqrt(6)*(l1*l1*l1 + l2*l2*l2 + l3*l3*l3) + 0.5*(l1*l1 + l2*l2 + l3*l3)*(l1*l1 + l2*l2 + l3*l3);
    
    fprintf(fp2,"%.2f %.6f %.8e\n",landau_t,ep,(Energy - f0*PI)/ep);
  }
  
  fclose(fp);
  fclose(fp2);
  
}
