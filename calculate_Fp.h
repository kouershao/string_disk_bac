
void MyNode::calc_fp(double *Qb, double *fp) 
{
  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;
  double nx[K],ny[K],nz[K];
  if (strcmp(boundary,"fixs") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(phi[k]);
      ny[k] = sin(phi[k]);
      nz[k] = 0;
    }
  }

  if (strcmp(boundary,"minus2") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(-phi[k]);
      ny[k] = sin(-phi[k]);
      nz[k] = 0;
    }
  }
  
  if (strcmp(boundary,"half") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(phi[k]/2);
      ny[k] = sin(phi[k]/2);
      nz[k] = 0;
    }
  }

  if (strcmp(boundary,"minus1") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(-phi[k]/2);
      ny[k] = sin(-phi[k]/2);
      nz[k] = 0;
    }
  }


  if (strcmp(boundary,"three") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(3*phi[k]/2);
      ny[k] = sin(3*phi[k]/2);
      nz[k] = 0;
    }
  }

  if (strcmp(boundary,"minus3") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(-3*phi[k]/2);
      ny[k] = sin(-3*phi[k]/2);
      nz[k] = 0;
    }
  }
  
  
  if (strcmp(boundary,"mobius") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(phi[k]/2)*cos(phi[k]);
      ny[k] = cos(phi[k]/2)*sin(phi[k]);
      nz[k] = sin(phi[k]/2);
    }
  }

  if (strcmp(boundary,"symy") == 0) {
    for (int k = 0;k <= K/4;++k) {
      nx[k] = cos(phi[k])*cos(phi[k]);
      ny[k] = cos(phi[k])*sin(phi[k]);
      nz[k] = sin(phi[k]);
    }
    for (int k = K/4;k <= K/2;++k) {
      nx[k] = -nx[K/2 - k];
      ny[k] = ny[K/2 - k];
      nz[k] = nz[K/2 - k];
    }
    for (int k = K/2;k < K;++k) {
      nx[k] = nx[K - k];
      ny[k] = -ny[K - k];
      nz[k] = nz[K - k];
    }
  }
  
  if (strcmp(boundary,"double") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(2*phi[k]);
      ny[k] = sin(2*phi[k]);
      nz[k] = 0;
    }
  }

  if (strcmp(boundary,"minus4") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = cos(-2*phi[k]);
      ny[k] = sin(-2*phi[k]);
      nz[k] = 0;
    }
  }
  
  if (strcmp(boundary,"tang") == 0) {
    for (int k = 0;k < K;++k) {
      nx[k] = -sin(phi[k]);
      ny[k] = cos(phi[k]);
      nz[k] = 0;
    }
  }
  
  for (int ix = 0;ix < K;++ix) {
    fp[ix*5+0] = Qb[ix*5+0] - S0*(nx[ix]*nx[ix] - 1.0/3);
    fp[ix*5+1] = Qb[ix*5+1] - S0*nx[ix]*ny[ix];
    fp[ix*5+2] = Qb[ix*5+2] - S0*nx[ix]*nz[ix];
    fp[ix*5+3] = Qb[ix*5+3] - S0*(ny[ix]*ny[ix] - 1.0/3);
    fp[ix*5+4] = Qb[ix*5+4] - S0*ny[ix]*nz[ix];
  }
}

double MyNode::calc_energy_Fp(double *Q) 
{
  int ix,j,k,n,dim;
  double Fpena = 0;
  double *fp;
  double *Qb = new double[K*5]();
  tranQ2b(Q,Qb);
  
  dim = 5;
  
  fp = new double[K*dim]();
  calc_fp(Qb,fp);
  ix = 0;
  for (k=0;k<K;k++)
    for (n=0;n<dim;++n)
      Fpena += 2*PI/K*fp[ix]*fp[ix++];
  
  delete [] Qb;
  delete [] fp;
  return Fpena;
}

void MyNode::calc_grad_energy_Fp(double *Q,double *grad_Fp) 
{
  int i,k,n,m,s;
  int ix;
  double *value;
  double *fp;
  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;
  
  double *Qb = new double[K*5]();
  tranQ2b(Q,Qb);
  
  double *Pm = new double[(2*M-1)]();
  double *der_fp = new double[5*K]();
  

  if (1) {
    fp = new double[K*5]();
    calc_fp(Qb,fp);
    for (n=0;n<5;++n)
      for (ix = 0;ix < K;++ix)
	der_fp[n*K + ix] = 2*fp[ix*5 + n];
  }
  
  for (s=0;s<5;++s) {
    ix = 0;
    for (m=1-M;m<=M-1;m++) {
      Pm[m+M-1] = 0;
      for (k=0;k<K;k++)
	Pm[m+M-1] += der_fp[s*K + k]*Xm[k*(2*M-1)+m+M-1];
      Pm[m+M-1] *= 2*PI/K;
      value = Rnmr + I*M/2*(N-M/2+1) + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4;
      for (n=abs(m);n<N;n+=2)
	grad_Fp[s*Basis + (ix++)] = Pm[m+M-1]*(*value++);
    }
  }
  
  delete [] Pm;
  delete [] der_fp;
  
  delete [] Qb;
  delete [] fp;
}
