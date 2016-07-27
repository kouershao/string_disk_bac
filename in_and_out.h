

void MyNode::getfname(char fname[], int N, int M) {
  char strR[20];
  char strN[20];
  char strM[20];
  
  if (Rad < 10)
    sprintf(strR,"0%.2f",Rad);
  else
    sprintf(strR,"%.2f",Rad);

  if (N < 10)
    sprintf(strN,"0%d",N);
  else
    sprintf(strN,"%d",N);

  if (M < 10)
    sprintf(strM,"0%d",M);
  else
    sprintf(strM,"%d",M);
  
  //sprintf(fname,"%s/Result/%s%s_t_%+.2f_R_%s_N_%s_M_%s_eta_%.1e_L21_%+.1f_%s.txt",DIR,func_type,boundary,t,strR,strN,strM,eta,L21,num);  
  sprintf(fname,"%s/Result/%s%s_t_%+.2f_R_%s_N_%s_M_%s_eta_%.1e_%s.txt",DIR,func_type,boundary,landau_t,strR,strN,strM,eta,suffix.str().c_str());
  printf("%s\n", fname);
}

/***************************************************************************/

void MyNode::iput(int FN,int FM) {
  int i,n,m,ix,err;
  int MaxN;
  FILE *fp;
  char path[200];
  //char *path;
  
  if (readfile == 0)
    initial();
  else {
    if (N >= FN)
      MaxN = N;
    else
      MaxN = FN;

    double *DataAnm = new double[5*MaxN*(2*MaxN-1)]();
    getfname(path, FN, FM);
    if ((fp = fopen(path,"r")) == NULL) printf("Open Anm Error!\n");
    for (i=0;i<5;i++)
      for (m=1-FM;m<=FM-1;m++)
	for (n=abs(m);n<FN;n+=2)
	{
		err = fscanf(fp,"%lf\n",&DataAnm[i*MaxN*(2*MaxN-1)+n*(2*MaxN-1)+m+MaxN-1]);

	}

	fclose(fp);

	ix = 0;

	for (i=0;i<5;i++)
	{
	  for (m=1-M;m<=M-1;m++)
	  {  
		  for (n=abs(m);n<N;n+=2)
		{
		Anm[ix++] = DataAnm[i*MaxN*(2*MaxN-1)+n*(2*MaxN-1)+m+MaxN-1];
		}
	  }
	}
    delete [] DataAnm;
  }
}

/***************************************************************************/

void MyNode::oput(int FN,int FM) {
  int i,n,m,ix,err;
  FILE* fp;
  char path[200];
  
  getfname(path, FN, FM);
  fp = fopen(path,"w");
  for (int i=0;i<5*Basis;i++)
    fprintf(fp,"%16.15e\n",Anm[i]);
  fclose(fp);
}

/***************************************************************************/
