
double MyNode::calc_energy_Fe(double *A,double *V) {
  int i,n,s,l,m;
  double Felas = 0;
  tranA2V(A,V);
  for (m=1-M;m<=M-1;m++) {
    for (n=abs(m);n<N;n+=2) {
      for (s=abs(m);s<N;s+=2) {
	for (i=0;i<5;i++)
	  Felas += Knm[n*N*N + s*N + abs(m)]*V[inx(i,n,m)]*V[inx(i,s,m)];
	Felas += Knm[n*N*N + s*N + abs(m)]*V[inx(0,n,m)]*V[inx(3,s,m)];
      }
    }
  }
  return Felas;
}

void MyNode::calc_grad_energy_Fe(double *A,double *V,double *grad_Fe) {
  int i,n,s,l,m,ix;
  double vk;
  tranA2V(A,V);
  ix = 0;
  for (m=1-M;m<=M-1;m++) {
    for (n=abs(m);n<N;n+=2) {
      for (i=0;i<5;i++)
	grad_Fe[i*Basis + ix] = 0;
      for (s=abs(m);s<N;s+=2) {
	grad_Fe[0*Basis + ix] += Knm[n*N*N + s*N + abs(m)]*(2*V[inx(0,s,m)] + V[inx(3,s,m)]);
	grad_Fe[1*Basis + ix] += 2*Knm[n*N*N + s*N + abs(m)]*V[inx(1,s,m)];
	grad_Fe[2*Basis + ix] += 2*Knm[n*N*N + s*N + abs(m)]*V[inx(2,s,m)];
	grad_Fe[3*Basis + ix] += Knm[n*N*N + s*N + abs(m)]*(V[inx(0,s,m)] + 2*V[inx(3,s,m)]);
	grad_Fe[4*Basis + ix] += 2*Knm[n*N*N + s*N + abs(m)]*V[inx(4,s,m)];
      }
      ++ix;
    }
  }
}

double MyNode::calc_FeL2_and_gradFeL2(double *qnm) {
  int i,k;
  double intG;
  
  for (i=0;i<5;i++) {
    calc_dr_fik(qnm + i*Basis,dr_qik + i*Point);
    calc_dp_fik(qnm + i*Basis,dp_qik + i*Point);
  }
  
  for (int ix = 0;ix < innerPoint;ix++) {
    Q1x[ix] = dr_qik[ix]*drdx[ix] + dp_qik[ix]*dpdx[ix];
    Q2x[ix] = dr_qik[ix + Point]*drdx[ix] + dp_qik[ix + Point]*dpdx[ix];
    Q2y[ix] = dr_qik[ix + Point]*drdy[ix] + dp_qik[ix + Point]*dpdy[ix];
    Q3x[ix] = dr_qik[ix + 2*Point]*drdx[ix] + dp_qik[ix + 2*Point]*dpdx[ix];
    Q4y[ix] = dr_qik[ix + 3*Point]*drdy[ix] + dp_qik[ix + 3*Point]*dpdy[ix];
    Q5y[ix] = dr_qik[ix + 4*Point]*drdy[ix] + dp_qik[ix + 4*Point]*dpdy[ix];
    
    QL2[ix] = (Q1x[ix] + Q2y[ix])*(Q1x[ix] + Q2y[ix]) + (Q2x[ix] + Q4y[ix])*(Q2x[ix] + Q4y[ix]) + (Q3x[ix] + Q5y[ix])*(Q3x[ix] + Q5y[ix]);
  }
  
  intG = intdisk(QL2);

  for (int ix = 0;ix < innerPoint;ix++) {
    grad_QL2_dr[ix] = 2*(Q1x[ix] + Q2y[ix])*drdx[ix];
    grad_QL2_dp[ix] = 2*(Q1x[ix] + Q2y[ix])*dpdx[ix];
    
    grad_QL2_dr[Point + ix] = 2*(Q1x[ix] + Q2y[ix])*drdy[ix] + 2*(Q2x[ix] + Q4y[ix])*drdx[ix];
    grad_QL2_dp[Point + ix] = 2*(Q1x[ix] + Q2y[ix])*dpdy[ix] + 2*(Q2x[ix] + Q4y[ix])*dpdx[ix];
    
    grad_QL2_dr[2*Point + ix] = 2*(Q3x[ix] + Q5y[ix])*drdx[ix];
    grad_QL2_dp[2*Point + ix] = 2*(Q3x[ix] + Q5y[ix])*dpdx[ix];

    grad_QL2_dr[3*Point + ix] = 2*(Q2x[ix] + Q4y[ix])*drdy[ix];
    grad_QL2_dp[3*Point + ix] = 2*(Q2x[ix] + Q4y[ix])*dpdy[ix];

    grad_QL2_dr[4*Point + ix] = 2*(Q3x[ix] + Q5y[ix])*drdy[ix];
    grad_QL2_dp[4*Point + ix] = 2*(Q3x[ix] + Q5y[ix])*dpdy[ix];
  }

  for (int i=0;i<5;i++) {
    calc_dr_fnm(grad_QL2_dr + i*Point,grad_FeL2_r + i*Basis);
    calc_dp_fnm(grad_QL2_dp + i*Point,grad_FeL2_p + i*Basis);
  }
  
  for (int i=0;i<5*Basis;i++) {
    grad_FeL2[i] = grad_FeL2_r[i] + grad_FeL2_p[i];
  }
  
  return intG;
}
