
void MyNode::calc_bulkenergy_landau_de(double *qik,double *fbulk) 
{
  int i;
  double a,b,c;
  double trQ2,trQ3,trQ4;
  double q1,q2,q3,q4,q5;

  a = 0.5*landau_t;
  b = -sqrt(6);
  c = 0.5;

  for (i=0;i<innerPoint;i++) {
    q1 = qik[i];
    q2 = qik[Point + i];
    q3 = qik[2*Point + i];
    q4 = qik[3*Point + i];
    q5 = qik[4*Point + i];

    trQ2 = 2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    trQ3 = 3*(2*q2*q3*q5 - (q1*q1-q2*q2+q3*q3)*q4 + q1*(q2*q2-q4*q4-q5*q5));
    trQ4 = trQ2*trQ2;

    fbulk[i] = a*trQ2 + b*trQ3 + c*trQ4;
  }


  /* for (i=0;i<innerPoint;i++) { */
  /*   q1 = qik[i]; */
  /*   q2 = qik[Point + i]; */
    
  /*   trQ2 = (q1*q1 + q2*q2 + (q1+q2)*(q1+q2)); */
  /*   trQ3 = (q1*q1*q1 + q2*q2*q2 - (q1+q2)*(q1+q2)*(q1+q2)); */
  /*   trQ4 = trQ2*trQ2; */

  /*   fbulk[i] = a*trQ2 + b*trQ3 + c*trQ4; */
  /* } */


}

void MyNode::calc_grad_bulkenergy_landau_de(double *qik,double *grad_fbulk) 
{
  int i,n;
  double a,b,c;
  double trQ2,trQ3,trQ4;
  double q1,q2,q3,q4,q5;

  a = 0.5*landau_t;
  b = -sqrt(6);
  c = 0.5;
  
  for (i=0;i<innerPoint;i++) {
    q1 = qik[i];
    q2 = qik[Point + i];
    q3 = qik[2*Point + i];
    q4 = qik[3*Point + i];
    q5 = qik[4*Point + i];
    
    trQ2 = 2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    
    grad_fbulk[i] = a*(4*q1+2*q4) + b*3*(q2*q2-q4*q4-q5*q5-2*q1*q4) + c*2*trQ2*(4*q1+2*q4);
    grad_fbulk[Point + i] = a*4*q2 + b*6*(q3*q5 + q2*q4 + q1*q2) + c*2*trQ2*4*q2;
    grad_fbulk[2*Point + i] = a*4*q3 + b*6*(q2*q5 - q3*q4) + c*2*trQ2*4*q3;
    grad_fbulk[3*Point + i] = a*(4*q4+2*q1) + b*3*(q2*q2-q1*q1-q3*q3-2*q1*q4) + c*2*trQ2*(4*q4+2*q1);
    grad_fbulk[4*Point + i] = a*4*q5 + b*6*(q2*q3 - q1*q5) + c*2*trQ2*4*q5;
  }



  /* for (i=0;i<innerPoint;i++) { */
  /*   q1 = qik[i]; */
  /*   q2 = qik[Point + i]; */
    
  /*   trQ2 = (q1*q1 + q2*q2 + (q1+q2)*(q1+q2)); */
  /*   trQ3 = (q1*q1*q1 + q2*q2*q2 - (q1+q2)*(q1+q2)*(q1+q2)); */
  /*   trQ4 = trQ2*trQ2; */
    
  /*   grad_fbulk[i] = a*2*(2*q1 + q2) + b*3*(q1*q1 - (q1+q2)*(q1+q2)) + c*2*trQ2*2*(2*q1 + q2); */
  /*   grad_fbulk[Point + i] = a*2*(2*q2 + q1) + b*3*(q2*q2 - (q1+q2)*(q1+q2)) + c*2*trQ2*2*(2*q2 + q1); */
  /* } */

}

double MyNode::calc_energy_Fb(double *qnm,double *qik,double *fbulk) 
{
  int i,k;
  double intG;

  for (i=0;i<5;i++)
    calc_fik(qnm + i*Basis,qik + i*Point);
  
  calc_bulkenergy_landau_de(qik,fbulk);
  
  intG = intdisk(fbulk);
  
  return intG;
}
