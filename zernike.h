
/* fnm坐标顺序为(m,n)，一共有 3*N/2 + (2*N-M)*(M/2-1) 个 */
/* Rnmr坐标顺序为(i,m,n),储存的位置为：Rnmr[i*M/2*(N-M/2+1)+(2*N*m+2*n-m*m)/4] */
/* Xm坐标顺序为(k,m),储存的位置为：Xm[k*(2*M-1)+(m+M-1)] */
/* coe_int_r2(m,n,i)储存的位置为：coe_int_r2[(2*N*m+2*n-m*m)/4*I+i] */
/* hmr坐标顺序为(i,m),元素个数为(I+1)*(2*M-1) */
/* Kim坐标顺序为(m,i) */
/* 坐标的测试程序为test_idx.m */

/*********************************************************************/

void MyNode::sphere_param_init() 
{
  Basis = 3*N/2 + (2*N-M)*(M/2-1);
  Point = (I+1)*K;
  innerPoint = I*K;
}

void MyNode::zernike_init() 
{
  int i,k,n,m,s;
  int read_ok;
  int ix;
  double x,y;
  FILE *fp;
  char path[200];

  radius = new double[I+1]();
  phi = new double[K]();
  xb = new double[K]();
  yb = new double[K]();
  coe_r = new double[I]();
 
  Rnmr = new double[(I+1)*M/2*(N-M/2+1)]();
  Xm = new double[K*(2*M-1)]();

  dRnm = new double[(I+1)*M/2*(N-M/2+1)]();
  dXm = new double[K*(2*M-1)]();
  
  coe_int_r = new double[I]();
  coe_int_r2 = new double[M/2*(N-M/2+1)*I]();
  
  coe_int_dr = new double[M/2*(N-M/2+1)*I]();

  hmr = new double[(I+1)*(2*M-1)]();
  Kim = new double[I*(2*M-1)]();
 
  Knm = new double[N*N*N]();

  FDI_p = new double[K/2+1]();
  FDO_p = new double[K/2+1]();
  FFTp_p = fftw_plan_r2r_1d(K/2+1,FDI_p,FDO_p,FFTW_REDFT00,FFTW_MEASURE);
  FFTq_p = fftw_plan_r2r_1d(K/2-1,FDI_p,FDO_p,FFTW_RODFT00,FFTW_MEASURE);
  
  FDI_c = new double[K+1]();
  FDO_c = new double[K+1]();
  FFTp_c = fftw_plan_r2r_1d(K+1,FDI_c,FDO_c,FFTW_REDFT00,FFTW_MEASURE);
  FFTq_c = fftw_plan_r2r_1d(K-1,FDI_c,FDO_c,FFTW_RODFT00,FFTW_MEASURE);
  
  /*********************************************************************/
  
  sprintf(path,"%s/Parameter/roots_I_%d.txt",DIR,I);
  if ((fp = fopen(path,"r")) == NULL) printf("Open roots_I Error!\n");
  for (i=0;i<=I;i++)
    read_ok = fscanf(fp,"%lf\n",radius+i);
  fclose(fp);

  sprintf(path,"%s/Parameter/coe_r_%d.txt",DIR,I);
  if ((fp = fopen(path,"r")) == NULL) printf("Open coe_r Error!\n");
  for (i=0;i<I;i++)
    read_ok = fscanf(fp,"%lf\n",coe_r+i);
  fclose(fp);
  
  for (k=0;k<K;k++)
    phi[k] = 2.0*PI*k/K;
  
  /*********************************************************************/
  
  for (k=0;k<K;k++) {
    xb[k] = cos(phi[k]);
    yb[k] = sin(phi[k]);
  }
  
  /*********************************************************************/
  
  double *DataRnmr = new double[(I+1)*MaxN*MaxN]();
  sprintf(path,"%s/Parameter/Rnmr_N_%d_I_%d.txt",DIR,MaxN,I);
  if ((fp = fopen(path,"r")) == NULL) printf("Open Rnmr Error!\n");
  for (i=0;i<=I;i++)
    for (m=0;m<MaxN;m++)
      for (n=m;n<MaxN;n+=2)
	read_ok = fscanf(fp,"%lf\n",&DataRnmr[i*MaxN*MaxN+m*MaxN+n]);
  fclose(fp);

  ix = 0;
  for (i=0;i<=I;i++)
    for (m=0;m<M;m++)
      for (n=m;n<N;n+=2)
	Rnmr[ix++] = DataRnmr[i*MaxN*MaxN+m*MaxN+n];
  
  for (i=0;i<I;i++)
    coe_int_r[i] = coe_r[i]*radius[i];
  
  ix = 0;
  for (m=0;m<M;m++)
    for (n=m;n<N;n+=2)
      for (i=0;i<I;i++)
	coe_int_r2[ix++] = coe_int_r[i]*DataRnmr[i*MaxN*MaxN+m*MaxN+n];
  
  delete [] DataRnmr;
  
  double *DatadRnm = new double[(I+1)*MaxN*MaxN]();
  if (I < 200 && I >= 8) {
    sprintf(path,"%s/Parameter/dRnm_N_%d_I_%d.txt",DIR,MaxN,I);
    if ((fp = fopen(path,"r")) == NULL) printf("Open dRnm Error!\n");
    for (i=0;i<=I;i++)
      for (m=0;m<MaxN;m++)
	for (n=m;n<MaxN;n+=2)
	  read_ok = fscanf(fp,"%lf\n",&DatadRnm[i*MaxN*MaxN+m*MaxN+n]);
    fclose(fp);
    
    ix = 0;
    for (i=0;i<=I;i++)
      for (m=0;m<M;m++)
	for (n=m;n<N;n+=2)
	  dRnm[ix++] = DatadRnm[i*MaxN*MaxN+m*MaxN+n];
    
    ix = 0;
    for (m=0;m<M;m++)
      for (n=m;n<N;n+=2)
	for (i=0;i<I;i++)
	  coe_int_dr[ix++] = coe_int_r[i]*DatadRnm[i*MaxN*MaxN+m*MaxN+n];
  }
  delete [] DatadRnm;
  
  /******************************************************************/
  
  ix = 0;
  for (k=0;k<K;k++) {
    for (m=1-M;m<0;m++)
      Xm[ix++] = 1.0/sqrt(PI)*sin(abs(m)*phi[k]);
    
    Xm[ix++] = 1.0/sqrt(2.0*PI);
    
    for (m=1;m<=M-1;m++)
      Xm[ix++] = 1.0/sqrt(PI)*cos(m*phi[k]);
  }

  ix = 0;
  for (k=0;k<K;k++) {
    for (m=1-M;m<0;m++)
      dXm[ix++] = 1.0/sqrt(PI)*abs(m)*cos(abs(m)*phi[k]);
    
    dXm[ix++] = 0;
    
    for (m=1;m<=M-1;m++)
      dXm[ix++] = -1.0/sqrt(PI)*m*sin(m*phi[k]);
  }
  
  /******************************************************************/
  
  sprintf(path,"%s/Parameter/Knm_N_%d.txt",DIR,MaxN);
  if ((fp = fopen(path,"r")) == NULL) printf("Open Knm Error!\n");
  
  for (n=0;n<N;n++)
    for (s=0;s<N;s++)
      for (m=0;m<N;m++)
	Knm[n*N*N + s*N + m] = 0;
  
  for (n=0;n<N;n++)
    for (s=n;s>=0;s-=2)
      for (m=s;m>=0;m-=2)
  	read_ok = fscanf(fp,"%lf\n",&Knm[n*N*N + s*N + m]);
  fclose(fp);
  
  for (n=0;n<N;n++)
    for (s=n;s<N;s+=2)
      for (m=n;m>=0;m-=2)
  	Knm[n*N*N + s*N + m] = Knm[s*N*N + n*N + m];
}

/*********************************************************************/

void MyNode::zer_destroy() 
{
  delete [] radius;
  delete [] phi;
  delete [] xb;
  delete [] yb;
  delete [] coe_r;

  delete [] Rnmr;
  delete [] Xm;
  
  delete [] coe_int_r;
  delete [] coe_int_r2;
 
  delete [] dRnm;
  delete [] dXm;
  delete [] coe_int_dr;

  delete [] hmr;
  delete [] Kim;
  delete [] Knm;
  
  fftw_destroy_plan(FFTp_p);
  fftw_destroy_plan(FFTq_p);
  delete [] FDI_p;
  delete [] FDO_p;
  
  fftw_destroy_plan(FFTp_c);
  fftw_destroy_plan(FFTq_c);
  delete [] FDI_c;
  delete [] FDO_c;
}

/*********************************************************************/

void MyNode::calc_fik(double *fnm,double *fik) 
{

  int i,k,n,m;
  double *value;
  double *ix,*ix2;

  value = hmr;
  for (i=0;i<=I;i++) {
    ix = fnm;
    for (m=1-M;m<=M-1;m++) {
      ix2 = Rnmr + i*M/2*(N-M/2+1) + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4;
      *value = 0;
      for (n=abs(m);n<N;n+=2)
	*value += (*ix++)*(*ix2++);
      ++value;
    }
  }
  
  if (K >= 2*M) {
    for (i=0;i<=I;i++) {
      for (m=1-M;m<0;m++)
	hmr[i*(2*M-1)+m+M-1] *= 1.0/sqrt(PI);
      hmr[i*(2*M-1)+M-1] *= 1.0/sqrt(2.0*PI);
      for (m=1;m<=M-1;m++)
	hmr[i*(2*M-1)+m+M-1] *= 1.0/sqrt(PI);
    }

    for (i=0;i<=I;i++) {
      for (m=0;m<M;m++)
  	FDI_p[m] = hmr[i*(2*M-1)+m+M-1];
      for (m=M;m<K/2+1;m++)
  	FDI_p[m] = 0;
      fftw_execute(FFTp_p);
      for (k=0;k<K/2;k++)
  	fik[i*K+k] = 0.5*(FDO_p[k] + FDI_p[0]);
      
      for (m=1;m<M;m+=2)
  	FDI_p[m] = -FDI_p[m];
      fftw_execute(FFTp_p);
      for (k=0;k<K/2;k++)
  	fik[i*K+k+K/2] = 0.5*(FDO_p[k] + FDI_p[0]);
      
      for (m=0;m<M-1;m++)
  	FDI_p[m] = hmr[i*(2*M-1)+M-2-m];
      for (m=M-1;m<K/2+1;m++)
  	FDI_p[m] = 0;
      fftw_execute(FFTq_p);
      for (k=1;k<K/2;k++)
  	fik[i*K+k] += 0.5*FDO_p[k-1];
      
      for (m=0;m<M-1;m=m+2)
  	FDI_p[m] = -FDI_p[m];
      fftw_execute(FFTq_p);
      for (k=1;k<K/2;k++)
  	fik[i*K+k+K/2] += 0.5*FDO_p[k-1];
    }
  }
  else {
    for (i=0;i<=I;i++) {
      for (k=0;k<K;k++) {
    	fik[i*K+k] = 0;
    	for (m=0;m<2*M-1;m++) {
    	  fik[i*K+k] += hmr[i*(2*M-1)+m]*Xm[k*(2*M-1)+m];
    	}
      }
    }
  }

  //std::cout << "~~~~~~~" <<std::endl;
}

/*********************************************************************/

void MyNode::calc_dr_fik(double *fnm,double *dr_fik) 
{
  int i,k,n,m;
  double *value;
  double *ix,*ix2;
  
  value = hmr;
  for (i=0;i<I;i++) {
    ix = fnm;
    for (m=1-M;m<=M-1;m++) {
      ix2 = dRnm + i*M/2*(N-M/2+1) + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4;
      *value = 0;
      for (n=abs(m);n<N;n+=2)
	*value += (*ix++)*(*ix2++);
      ++value;
    }
  }

  if (K >= 2*M) {
    for (i=0;i<I;i++) {
      for (m=1-M;m<0;m++)
	hmr[i*(2*M-1)+m+M-1] *= 1.0/sqrt(PI);
      hmr[i*(2*M-1)+M-1] *= 1.0/sqrt(2.0*PI);
      for (m=1;m<=M-1;m++)
	hmr[i*(2*M-1)+m+M-1] *= 1.0/sqrt(PI);
    }
    
    for (i=0;i<I;i++) {
      for (m=0;m<M;m++)
	FDI_p[m] = hmr[i*(2*M-1)+m+M-1];
      for (m=M;m<K/2+1;m++)
	FDI_p[m] = 0;
      fftw_execute(FFTp_p);
      for (k=0;k<K/2;k++)
	dr_fik[i*K+k] = 0.5*(FDO_p[k] + FDI_p[0]);
	
      for (m=1;m<M;m+=2)
	FDI_p[m] = -FDI_p[m];
      fftw_execute(FFTp_p);
      for (k=0;k<K/2;k++)
	dr_fik[i*K+k+K/2] = 0.5*(FDO_p[k] + FDI_p[0]);
      
      for (m=0;m<M-1;m++)
	FDI_p[m] = hmr[i*(2*M-1)+M-2-m];
      for (m=M-1;m<K/2+1;m++)
	FDI_p[m] = 0;
      fftw_execute(FFTq_p);
      for (k=1;k<K/2;k++)
	dr_fik[i*K+k] += 0.5*FDO_p[k-1];
      
      for (m=0;m<M-1;m=m+2)
	FDI_p[m] = -FDI_p[m];
      fftw_execute(FFTq_p);
      for (k=1;k<K/2;k++)
	dr_fik[i*K+k+K/2] += 0.5*FDO_p[k-1];
    }
  }
  else {
    for (i=0;i<I;i++) {
      for (k=0;k<K;k++) {
	dr_fik[i*K+k] = 0;
	for (m=0;m<2*M-1;m++) {
	  dr_fik[i*K+k] += hmr[i*(2*M-1)+m]*Xm[k*(2*M-1)+m];
	}
      }
    }
  }
}

/*********************************************************************/

void MyNode::calc_dp_fik(double *fnm,double *dp_fik) 
{
  int i,k,n,m;
  double *value;
  double *ix,*ix2;
  
  value = hmr;
  for (i=0;i<I;i++) {
    ix = fnm;
    for (m=1-M;m<=M-1;m++) {
      ix2 = Rnmr + i*M/2*(N-M/2+1) + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4;
      *value = 0;
      for (n=abs(m);n<N;n+=2)
	*value += (*ix++)*(*ix2++);
      ++value;
    }
  }
  
  for (i=0;i<I;i++) {
    for (k=0;k<K;k++) {
      dp_fik[i*K+k] = 0;
      for (m=0;m<2*M-1;m++) {
	dp_fik[i*K+k] += hmr[i*(2*M-1)+m]*dXm[k*(2*M-1)+m];
      }
    }
  }
}

/*********************************************************************/

void MyNode::calc_fnm(double *fik,double *fnm) 
{
  int i,k,n,m;
  double *value;
  double *ix,*ix2;
  int idx;
  
  for (i=0;i<I;i++) {
    for (k=0;k<K;k++)
      FDI_c[k] = fik[i*K+k];
    fftw_execute(FFTp_c);
    for (m=0;m<M;m++)
      Kim[(m+M-1)*I+i] = 0.5*(FDO_c[2*m] + FDI_c[0]);
      
    for (k=0;k<K-1;k++)
      FDI_c[k] = fik[i*K+k+1];
    fftw_execute(FFTq_c);
    for (m=1;m<M;m++)
      Kim[(M-1-m)*I+i] = 0.5*FDO_c[2*m-1];
  }

  for (i=0;i<I;i++) {
    for (m=1-M;m<0;m++)
      Kim[(m+M-1)*I+i] *= 1.0/sqrt(PI);
    Kim[(M-1)*I+i] *= 1.0/sqrt(2.0*PI);
    for (m=1;m<=M-1;m++)
      Kim[(m+M-1)*I+i] *= 1.0/sqrt(PI);
  }
   
  value = fnm;
  idx = 0;
  for (m=1-M;m<=M-1;m++) {
    ix = coe_int_r2 + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4*I;
    for (n=abs(m);n<N;n+=2) {
      *value = 0;
      for (i=0;i<I;i++)
	*value += (*ix++)*Kim[idx*I+i];
      *value = 2*PI/K*(*value);
      ++value;
    }
    ++idx;
  }
}

/*********************************************************************/

void MyNode::calc_dr_fnm(double *fik,double *dr_fnm) 
{
  int i,k,n,m;
  double *value;
  double *ix,*ix2;
  int idx;
  
  for (i=0;i<I;i++) {
    for (k=0;k<K;k++)
      FDI_c[k] = fik[i*K+k];
    fftw_execute(FFTp_c);
    for (m=0;m<M;m++)
      Kim[(m+M-1)*I+i] = 0.5*(FDO_c[2*m] + FDI_c[0]);
    
    for (k=0;k<K-1;k++)
      FDI_c[k] = fik[i*K+k+1];
    fftw_execute(FFTq_c);
    for (m=1;m<M;m++)
      Kim[(M-1-m)*I+i] = 0.5*FDO_c[2*m-1];
  }

  for (i=0;i<I;i++) {
    for (m=1-M;m<0;m++)
      Kim[(m+M-1)*I+i] *= 1.0/sqrt(PI);
    Kim[(M-1)*I+i] *= 1.0/sqrt(2.0*PI);
    for (m=1;m<=M-1;m++)
      Kim[(m+M-1)*I+i] *= 1.0/sqrt(PI);
  }
  
  value = dr_fnm;
  idx = 0;
  for (m=1-M;m<=M-1;m++) {
    ix = coe_int_dr + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4*I;
    for (n=abs(m);n<N;n+=2) {
      *value = 0;
      for (i=0;i<I;i++)
	*value += (*ix++)*Kim[idx*I+i];
      *value = 2*PI/K*(*value);
      ++value;
    }
    ++idx;
  }
}

/*********************************************************************/

void MyNode::calc_dp_fnm(double *fik,double *dp_fnm) 
{
  int i,k,n,m;
  double *value;
  double *ix,*ix2;
  int idx;
  
  for (m=0;m<2*M-1;m++) {
    for (i=0;i<I;i++) {
      Kim[m*I+i] = 0;
      for (k=0;k<K;k++)
	Kim[m*I+i] += fik[i*K+k]*dXm[k*(2*M-1)+m];
    }
  }
  
  value = dp_fnm;
  idx = 0;
  for (m=1-M;m<=M-1;m++) {
    ix = coe_int_r2 + (2*N*abs(m)+2*abs(m)-abs(m)*abs(m))/4*I;
    for (n=abs(m);n<N;n+=2) {
      *value = 0;
      for (i=0;i<I;i++)
	*value += (*ix++)*Kim[idx*I+i];
      *value = 2*PI/K*(*value);
      ++value;
    }
    ++idx;
  }
}

/*********************************************************************/
double MyNode::intdisk(double *f) 
{
  int i,k;
  double intf;
  double *G1 = new double[K]();
  
  intf = 0;
  for (k=0;k<K;k++) {
    G1[k] = 0;
    for (i=0;i<I;i++)
      G1[k] += coe_int_r[i]*f[i*K+k];
  }
  for (k=0;k<K;k++) {
    intf += G1[k];
  }
  intf = intf*2*PI/K;
  
  delete [] G1;
  
  return intf;
}
