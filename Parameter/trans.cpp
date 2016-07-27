#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DIR "/home/quyang/work/modprogram/"

const int I = 128;
const int J = 128;

int main () {
  FILE *fp;
  char path[100];
  int i,j,k,n,l,m,s;
  int read_ok;
  long double x,y;

  long double radius[I+1];
  long double coe_r[I];

  long double lambda[J];
  long double coe_theta[J];
  
  sprintf(path,"%s/Parameter/roots_%d.txt",DIR,I);
  if ((fp = fopen(path,"r")) == NULL) printf("Open roots_I Error!\n");
  read_ok = fscanf(fp,"%d\n",&i);
  for (i=0;i<I/2;i++) {
    read_ok = fscanf(fp,"%Lf %Lf\n",&x,&y);
    radius[I/2-1-i] = 0.5*(1-x);
    radius[I/2+i] = 0.5*(1+x);
    coe_r[I/2-1-i] = 0.5*y;
    coe_r[I/2+i] = 0.5*y;
  }
  fclose(fp);
  radius[I] = 1;

  sprintf(path,"%s/Parameter/roots_I_%d.txt",DIR,I);
  if ((fp = fopen(path,"w")) == NULL) printf("Open roots_I Error!\n");
  for (i=0;i<=I;i++) {
    fprintf(fp,"%21.20Lf\n",radius[i]);
  }
  fclose(fp);
  
  sprintf(path,"%s/Parameter/coe_r_%d.txt",DIR,I);
  if ((fp = fopen(path,"w")) == NULL) printf("Open roots_I Error!\n");
  for (i=0;i<I;i++) {
    fprintf(fp,"%21.20Lf\n",coe_r[i]);
  }
  fclose(fp);
  

  sprintf(path,"%s/Parameter/roots_%d.txt",DIR,J);
  if ((fp = fopen(path,"r")) == NULL) printf("Open Error!\n");
  read_ok = fscanf(fp,"%d\n",&j);
  for (j=0;j<J/2;j++) {
    read_ok = fscanf(fp,"%Lf %Lf\n",&x,&y);
    lambda[J/2-1-j] = -x;
    lambda[J/2+j] = x;
    coe_theta[J/2-1-j] = y;
    coe_theta[J/2+j] = y;
  }
  fclose(fp);

  sprintf(path,"%s/Parameter/roots_J_%d.txt",DIR,J);
  if ((fp = fopen(path,"w")) == NULL) printf("Open roots_J Error!\n");
  for (i=0;i<J;i++) {
    fprintf(fp,"%21.20Lf\n",lambda[i]);
  }
  fclose(fp);

  sprintf(path,"%s/Parameter/coe_theta_%d.txt",DIR,J);
  if ((fp = fopen(path,"w")) == NULL) printf("Open roots_J Error!\n");
  for (i=0;i<J;i++) {
    fprintf(fp,"%21.20Lf\n",coe_theta[i]);
  }
  fclose(fp);

  return 0;
}
