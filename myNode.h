#ifndef _NODE_H_
#define _NODE_H_

#include <iostream>
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <fftw3.h>
#include <lbfgs.h>
#include <eigen3/Eigen/Dense>

#define DIR "."
#define PI 3.141592653589793
#define func_type "poly"


class MyNode
{
	public:
		MyNode(int arg1 = 128,int arg2 = 64,int arg3 = 64,int arg4 = 128,int arg5 = 128, double arg6 = 1e5, double arg7 = 1000, double arg8 = 1000, double arg9 = 0, double arg10 = 1, double arg11 = 1.2, double arg12 = 0, double arg13 = 1, const char* arg14 = "fixs"):
			MaxN(arg1), N(arg2), M(arg3), I(arg4), K(arg5), eta(arg6), eta2(arg7), eta3(arg8), L21(arg9), func_B(arg10), func_C(arg11), uniaxial(arg12), anisotropy(arg13), boundary(arg14)
			{
				//		initialization();
			}

		// static member
		const int MaxN, N, M, I, K;
		int Basis, Point, innerPoint; 
		double landau_t, Rad;

		// non-static
		int readfile, mode;
		std::stringstream suffix;
		double Energy, *grad_Energy, *Anm, *Vnm, *Qik;
		lbfgs_parameter_t param;

		void node_initialization();
		void cal_dF(double*, double*, double*, double*); 
		double cal_F(double*, double*, double*);
		void iput(int, int);
		void oput(int, int);
		double Norm(double*, int); 
		void var_destroy();
		void zer_destroy(); 

		int run(int i, int n, Eigen::MatrixXd &gradient)
		{
			Energy = cal_F(Anm, Vnm, Qik);
			cal_dF(Anm, Vnm, Qik, grad_Energy);	
			for( int j = 0; j < 5*Basis; j += 1)
			{
				gradient(i, j) = grad_Energy[j];
			}
			std::cout << i << " " << Energy*2*PI << " " << Norm(grad_Energy, n) << "|"; 
			//	std::cout << "node = " << i << " "; 
			return lbfgs(n, Anm, &Energy, _evaluate, _progress, this, &param);
		}


	protected:
		static lbfgsfloatval_t _evaluate(
				void *instance,
				const double *Bnm,
				double *grad_Energy,
				const int n,
				const lbfgsfloatval_t step
				)
		{
			return reinterpret_cast<MyNode*>(instance)->evaluate(Bnm, grad_Energy, n, step);
		}

		lbfgsfloatval_t evaluate(
				const double *Bnm,
				double *grad_Energy,
				const int n,
				const lbfgsfloatval_t step
				)
		{
			Energy = cal_F(Anm,Vnm,Qik);
			cal_dF(Anm,Vnm,Qik,grad_Energy);
			return Energy;
		}

		static int _progress(
				void *instance,
				const double *Anm,
				const double *grad_Energy,
				const double Energy,
				const lbfgsfloatval_t xnorm,
				const lbfgsfloatval_t gnorm,
				const lbfgsfloatval_t step,
				int n,
				int k,
				int ls
				)
		{
			return reinterpret_cast<MyNode*>(instance)->progress(Anm, grad_Energy, Energy, xnorm, gnorm, step, n, k, ls);
		}

		int progress(
				const double *Anm,
				const double *grad_Energy,
				const double Energy,
				const lbfgsfloatval_t xnorm,
				const lbfgsfloatval_t gnorm,
				const lbfgsfloatval_t step,
				int n,
				int k,
				int ls) {         
			//if (k % 100 == 0) 
			//{
				//printf("Iteration %d:  ",k);
				//printf("Energy = %16.15f  ",Energy*2*PI);
				//printf("normdF = %16.15f  step = %16.15f\n",gnorm,step);
			//}
			return 0;
		}


	private:
		const double eta, eta2, eta3, L21, func_B, func_C, uniaxial, anisotropy;
		const char* boundary;
		//double *V;
		//double *Q;
		//double *A;
		//void check_integrate();
		//void check_interpolation(int,int,char,double);
		//void perturbation();

		double *radius, *phi;
		double *xb,*yb;
		double *coe_r;
		double *Rnmr,*Xm;
		double *coe_int_r,*coe_int_r2;
		double *hmr;
		double *Kim;
		double *Knm;
		double *dRnm,*dXm;
		double *FDI_p,*FDO_p,*FDI_c,*FDO_c;
		fftw_plan FFTp_p,FFTq_p,FFTp_c,FFTq_c;

		//	double *Vnm;
		double *grad_Fb;
		double *grad_Fbeta;
		double *grad_Fe;
		double *grad_Fp;
		//	double *Qik;
		double *fbulk;
		double *fbeta;
		double *grad_fbulk;
		double *grad_fbeta;
		double *drdx,*drdy;
		double *dpdx,*dpdy;
		double *coe_int_dr;
		double *dr_qik,*dp_qik;
		double *Q1x,*Q2x,*Q2y,*Q3x,*Q4y,*Q5y;
		double *QL2,*grad_QL2_dr,*grad_QL2_dp;
		double *grad_FeL2,*grad_FeL2_r,*grad_FeL2_p;

		void getfname(char[], int, int);

		void sphere_param_init(); 
		void zernike_init(); 
		//		void zer_destroy(); 
		void calc_fik(double *,double *); 
		void calc_dr_fik(double *,double *); 
		void calc_dp_fik(double *,double *); 
		void calc_fnm(double *,double *); 
		void calc_dr_fnm(double *,double *); 
		void calc_dp_fnm(double *,double *); 
		double intdisk(double *); 

		void malloc_variables_init();
		//		void var_destroy();
		inline int inx(int,int,int); 
		void tranA2V(double *,double *);
		void tranQ2b(double *,double *); 
		void initial();  

		void calc_beta(double *);
		double calc_energy_beta(double *); 

		void calc_bulkenergy_landau_de(double *,double *); 
		void calc_grad_bulkenergy_landau_de(double *,double *); 
		double calc_energy_Fb(double *,double *,double *); 

		double calc_energy_Fe(double *,double *);
		void calc_grad_energy_Fe(double *,double *,double *);
		double calc_FeL2_and_gradFeL2(double *);

		double calc_energy_Fp(double *);
		void calc_fp(double *,double *);
		void calc_grad_energy_Fp(double *,double *);


};

#endif
