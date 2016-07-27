
void MyNode::calc_beta(double *Q) 
{
  int i;
  double trQ2,trQ3,trQ4;
  double q1,q2,q3,q4,q5;

  for (i=0;i<innerPoint;i++) {
    q1 = Q[i];
    q2 = Q[Point + i];
    q3 = Q[2*Point + i];
    q4 = Q[3*Point + i];
    q5 = Q[4*Point + i];
    trQ2 = 2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    trQ3 = 3*(2*q2*q3*q5 - (q1*q1-q2*q2+q3*q3)*q4 + q1*(q2*q2-q4*q4-q5*q5));
    
    fbeta[i] = trQ2*trQ2*trQ2 - 6*trQ3*trQ3;
    grad_fbeta[i] = 6*trQ2*trQ2*(2*q1+q4) - 36*trQ3*(q2*q2-q4*q4-q5*q5-2*q1*q4);
    grad_fbeta[Point + i] = 12*trQ2*trQ2*q2 - 72*trQ3*(q3*q5 + q2*q4 + q1*q2);
    grad_fbeta[2*Point + i] = 12*trQ2*trQ2*q3 - 72*trQ3*(q2*q5 - q3*q4);
    grad_fbeta[3*Point + i] = 6*trQ2*trQ2*(2*q4+q1) - 36*trQ3*(q2*q2-q1*q1-q3*q3-2*q1*q4);
    grad_fbeta[4*Point + i] = 12*trQ2*trQ2*q5 - 72*trQ3*(q2*q3 - q1*q5);
  }
}

double MyNode::calc_energy_beta(double *Q) 
{
  int i,k;
  double intG;
  
  calc_beta(Q);
  intG = intdisk(fbeta);

  return intG;
}
