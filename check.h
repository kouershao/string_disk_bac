
void Node::check_integrate() {
  int i;
  double err = 0;

  double *A = new double[Basis]();
  double *B = new double[Basis]();
  double *Q = new double[Point]();

  srand((unsigned)time(NULL));
  for (i=0;i<Basis;i++)
    A[i] = 0.2*(2.0*rand()/RAND_MAX - 1);
  
  calc_fik(A,Q);
  calc_fnm(Q,B);

  for (i=0;i<Basis;i++)
    err += fabs(A[i] - B[i]);

  printf("Error = %.4e\n",err);

  delete [] A;
  delete [] B;
  delete [] Q;
}

void Node::check_interpolation(int dis,int judge,char number[],double ep) {
  int i,j,k,n;
  zernike_init(32,4,32,12);
  double *A = new double[5*Basis]();
  double *Q1 = new double[5*Point]();
  
  if (judge == 0)
    initial(A,0);
  else
    iput(A,landau_t,Rad,eta,N,M,number,1,0);

  for (i=0;i<5;i++)
    calc_fik(A + i*Basis,Q1 + i*Point);
  zer_destroy();
  
  zernike_init(N,M,dis,dis);
  double *Q2 = new double[5*Point]();
  for (i=0;i<5;i++)
    calc_fik(A + i*Basis,Q2 + i*Point);
  zer_destroy();
  
  zernike_init(32,4,32,12);
  double *Q3 = new double[5*Point]();
  int ix = 0;
  for (i=0;i<=I;i++) {
    for (k=0;k<K;k++) {
      for (n=0;n<5;n++)
	//	Q3[ix + n*Point] = interpolation(Q2 + n*(dis+1)*(dis+1)*dis,dis,dis+1,dis,radius[i],theta[j],phi[k],dis);
	++ix;
    }
  }
  zer_destroy();
  
  for (i=0;i<5*Point;i++)
    if (fabs(Q1[i]-Q3[i]) > ep)
      printf("%.4e %.4e %.4e\n",Q1[i],Q3[i],Q1[i]-Q3[i]);
  
  delete [] Q1;
  delete [] Q2;
  delete [] Q3;
  delete [] A;
}

void Node::perturbation() {
  double Qscale = 3.0*sqrt(6)/2*func_C/func_B;
  double beig = Qscale/3.0;
  int i;
  srand((unsigned)time(NULL));
  for (i=0;i<5;i++)
    calc_fik(Anm + i*Basis,Qik + i*Point);
  for (i=0;i<5*Point;i++)
    Qik[i] += 0.1*(2.0*rand()/RAND_MAX - 1)*Qscale;
  for (i=0;i<5;i++)
    calc_fnm(Qik + i*Point,Anm + i*Basis);
}
