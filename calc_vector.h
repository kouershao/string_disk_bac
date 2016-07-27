
double MyNode::Norm(double *fnlm,int DIMM) 
{
  double norm;
  norm = 0;
  for (int i=0;i<DIMM;i++)
    norm = norm + fnlm[i]*fnlm[i];
  norm = sqrt(norm);
  return norm;
}
