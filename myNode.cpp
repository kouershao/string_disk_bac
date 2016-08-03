#include "myNode.h"
#include "zernike.h"
#include "initial.h"
#include "calculate_Fb.h"
#include "calculate_Fe.h"
#include "calculate_Fp.h"
#include "calcbeta.h"
#include "in_and_out.h"
#include "calc_vector.h"
//#include "check.h"


/*********************************************************************/
void MyNode::node_initialization()
{
	sphere_param_init(); 
	zernike_init();
	malloc_variables_init();
	lbfgs_parameter_init(&param);
	param.m = 10;
	param.epsilon = 1e-5;
	param.max_iterations = 100;
	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
//    check_integrate();
}

/*********************************************************************/

double MyNode::cal_F(double *A, double *V, double *Q) 
{
 
  double Fbulk,Felas,Fpena,Fbeta;
  Fbulk = calc_energy_Fb(A,Q,fbulk);

  if (anisotropy == 1)
    Felas = calc_energy_Fe(A,V) + 0.5*L21*calc_FeL2_and_gradFeL2(A);
  else
    Felas = calc_energy_Fe(A,V);
  
  Fpena = calc_energy_Fp(Q);

  Fbeta = 0;
  if (uniaxial == 1)
    Fbeta = calc_energy_beta(Q);
  
  return Fbulk + 1.0/Rad/Rad*Felas + eta*Fpena + eta2*Fbeta;
}

/*********************************************************************/

void MyNode::cal_dF(double *A, double *V, double *Q, double *grad_Energy) 
{
  calc_grad_bulkenergy_landau_de(Q,grad_fbulk);
  for (int i=0;i<5;i++)
    calc_fnm(grad_fbulk + i*Point,grad_Fb + i*Basis);
  calc_grad_energy_Fe(A,V,grad_Fe);
  calc_grad_energy_Fp(Q,grad_Fp);

  for (int i=0;i<5*Basis;i++)
    grad_Fbeta[i] = 0;
  
  if (uniaxial == 1) {
    for (int i=0;i<5;i++)
      calc_fnm(grad_fbeta + i*Point,grad_Fbeta + i*Basis);
  }
  
  if (anisotropy == 1) {
    for (int i=0;i<5*Basis;i++)
      grad_Energy[i] = grad_Fb[i] + 1.0/Rad/Rad*(grad_Fe[i] + 0.5*L21*grad_FeL2[i]) + eta*grad_Fp[i] + eta2*grad_Fbeta[i];
  }
  else {
    for (int i=0;i<5*Basis;i++)
      grad_Energy[i] = grad_Fb[i] + 1.0/Rad/Rad*grad_Fe[i] + eta*grad_Fp[i] + eta2*grad_Fbeta[i];
  }
}

/*********************************************************************/

//lbfgsfloatval_t MyNode::evaluate(void *instance,
		//const double *Bnm,
		//double *grad_Energy,
		//const int n,
		//const lbfgsfloatval_t step) {
	////double Fbulk,Felas,Fpena,Fbeta,Energy;
	//Energy = cal_F(Anm,Vnm,Qik);
	//cal_dF(Anm,Vnm,Qik,grad_Energy);
	//return Energy;
//}

//int MyNode::progress(void *instance,
		//const double *Anm,
		//const double *grad_Energy,
		//const double Energy,
		//const lbfgsfloatval_t xnorm,
		//const lbfgsfloatval_t gnorm,
		//const lbfgsfloatval_t step,
		//int n,
		//int k,
		//int ls) {         
	//if (k%100 == 0) 
	//{
		//printf("Iteration %d:  ",k);
		//printf("Energy = %16.15f  ",Energy);
		//printf("normdF = %16.15f  step = %16.15f\n",gnorm,step);
	//}
	//return 0;
//}

/*********************************************************************/
