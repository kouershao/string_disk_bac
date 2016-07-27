double interpolation(double *Q,int I,int K,double r,double p,int n) {
  int i,k;
  double dr,dp;
  double x,z;
  
  dr = 1.0/n;
  dp = 2.0*PI/n;
  i = r/dr;
  k = p/dp;
  if (i == I)
    i = I-1;
  x = r/dr - i;
  z = p/dp - k;

  if (k < K-1)
    return x*z*Q[(i+1)*K + (k+1)] + x*(1-z)*Q[(i+1)*K + k] + (1-x)*z*Q[i*K + (k+1)] + (1-x)*(1-z)*Q[i*K + k];
  else
    return x*z*Q[(i+1)*K + 0] + x*(1-z)*Q[(i+1)*K + k] + (1-x)*z*Q[i*K + 0] + (1-x)*(1-z)*Q[i*K + k];
}
