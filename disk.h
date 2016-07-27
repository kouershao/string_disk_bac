
const int RANDOM_MODE = 0;
const int RADIAL_MODE = 1;
const int SPLIT_MODE = 2;


double Qscale = 3.0*sqrt(6)/2*func_C/func_B;
double beig = Qscale/3.0;


/*********************************************************************/

double *radius,*phi;
double *xb,*yb;
double *coe_r;
double *Rnmr,*Xm;
double *coe_int_r,*coe_int_r2;
double *hmr;
double *Kim;
double *Knm;

double *FDI_p,*FDO_p,*FDI_c,*FDO_c;
fftw_plan FFTp_p,FFTq_p,FFTp_c,FFTq_c;

double *Anm;
double *Vnm;
double *grad_Fb;
double *grad_Fbeta;
double *grad_Fe;
double *grad_Fp;
double *grad_Energy;
double *Qik;
double *fbulk;
double *fbeta;
double *grad_fbulk;
double *grad_fbeta;

double *dRnm,*dXm;
double *coe_int_dr;
double *dr_qik,*dp_qik;
double *Q1x,*Q2x,*Q2y,*Q3x,*Q4y,*Q5y;
double *QL2,*grad_QL2_dr,*grad_QL2_dp;
double *grad_FeL2,*grad_FeL2_r,*grad_FeL2_p;

double *drdx,*drdy;
double *dpdx,*dpdy;

/*********************************************************************/

inline int inx(int i,int n,int m) {return i*(2*M-1)*N + (m+M-1)*N + n;}
