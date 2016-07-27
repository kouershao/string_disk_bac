
void MyNode::malloc_variables_init() {

  Anm = new double[5*Basis]();
  Vnm = new double[5*(2*M-1)*N]();
  grad_Fb = new double[5*Basis]();
  grad_Fbeta = new double[5*Basis]();
  grad_Fe = new double[5*Basis]();
  grad_Fp = new double[5*Basis]();
  grad_Energy = new double[5*Basis]();
  Qik = new double[5*Point]();
  fbulk = new double[innerPoint]();
  fbeta = new double[innerPoint]();
  grad_fbulk = new double[5*Point]();
  grad_fbeta = new double[5*Point]();

  drdx = new double[innerPoint]();
  drdy = new double[innerPoint]();
  dpdx = new double[innerPoint]();
  dpdy = new double[innerPoint]();
  
  /******************************************************************/
  
  for (int i=0;i<I;i++) {
    for (int k=0;k<K;k++) {
      drdx[i*K + k] = cos(phi[k]);
      drdy[i*K + k] = sin(phi[k]);
      dpdx[i*K + k] = -sin(phi[k])/radius[i];
      dpdy[i*K + k] = cos(phi[k])/radius[i];
    }
  }

  /******************************************************************/  
  
  dr_qik = new double[5*Point]();
  dp_qik = new double[5*Point]();

  Q1x = new double[innerPoint]();
  Q2x = new double[innerPoint]();
  Q2y = new double[innerPoint]();
  Q3x = new double[innerPoint]();
  Q4y = new double[innerPoint]();
  Q5y = new double[innerPoint]();

  QL2 = new double[innerPoint]();

  grad_QL2_dr = new double[5*Point]();
  grad_QL2_dp = new double[5*Point]();

  grad_FeL2 = new double[5*Basis]();
  grad_FeL2_r = new double[5*Basis]();
  grad_FeL2_p = new double[5*Basis]();
}

void MyNode::var_destroy() {
  delete [] Anm;
  delete [] Vnm;
  delete [] grad_Fb;
  delete [] grad_Fbeta;
  delete [] grad_Fe;
  delete [] grad_Fp;
  delete [] grad_Energy;
  delete [] Qik;
  delete [] fbulk;
  delete [] fbeta;
  delete [] grad_fbulk;
  delete [] grad_fbeta;

  delete [] drdx;
  delete [] drdy;
  delete [] dpdx;
  delete [] dpdy;
  
  delete [] dr_qik;
  delete [] dp_qik;
  
  delete [] Q1x;
  delete [] Q2x;
  delete [] Q2y;
  delete [] Q3x;
  delete [] Q4y;
  delete [] Q5y;

  delete [] QL2;
  delete [] grad_QL2_dr;
  delete [] grad_QL2_dp;
  
  delete [] grad_FeL2;
  delete [] grad_FeL2_r;
  delete [] grad_FeL2_p;
}

inline int MyNode::inx(int i,int n,int m) 
{
	return i*(2*M-1)*N + (m+M-1)*N + n;
}

void MyNode::tranA2V(double *A,double *V) 
{
  int i,n,m;
  for (i=0;i<5;i++)
    for (m=1-M;m<=M-1;m++)
      for (n=abs(m);n<N;n+=2)
	V[inx(i,n,m)] = (*A++);
}

void MyNode::tranQ2b(double *Q,double *Qb) 
{
  for (int k=0;k<K;k++)
    for (int n=0;n<5;++n)
      Qb[k*5+n] = Q[n*Point + I*K + k];
}

void MyNode::initial()  
{
  const int RANDOM_MODE = 0;
  const int RADIAL_MODE = 1;
  const int SPLIT_MODE = 2;
  int i,j,k,n;
  double *s = new double[I+1]();
  double *Q = new double[5*Point]();
  
  switch (mode) {
  case RANDOM_MODE:
    srand((unsigned)time(NULL));
    for (i=0;i<5*Basis;i++)
      Anm[i] = 0.2*(2.0*rand()/RAND_MAX - 1);
    break;
  }
  delete [] s;
  delete [] Q;
}
