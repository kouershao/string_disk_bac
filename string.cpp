#include "string.h"
#include <gsl/gsl_errno.h>                                              
#include <gsl/gsl_spline.h>

void MyString::interp1()
{ 
	
	double xx, yy, x[m], y[m];
	int i, j;
	gsl_interp_accel *acc 
			= gsl_interp_accel_alloc();
	gsl_spline *spline 
			 = gsl_spline_alloc(gsl_interp_cspline, m);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			y[j] = u_new(j, i);
			x[j] = dist(j);
		}

		gsl_spline_init(spline, x, y, m);
		//for(  j = 0; j < m; j += 1)
		//{
			//u(j, i) = gsl_spline_eval (spline, g1(j), acc);
		//}
		
		for (j = 0, xx = x[0]; xx <= x[m-1]; xx += 1.0 / (m - 1), j++)
		{
			u(j, i) = gsl_spline_eval (spline, xx, acc);
		}

	}
	//for( int i = 0; i < m; i += 1)
	//{
		//for( int j = 0; j < n*n; j += 1)
		//{
			//if (u_new(i, j)>0){ 		
				//std::cout << i << " " << j << " "<< u(i, j) << std::endl;
			//}
		//}
		
	//}

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
}

void MyString::distance()
{
	Eigen::MatrixXd du(m-1, n);
	dist.setZero(m,1);
	du = u_new.bottomRows(m-1) - u_new.topRows(m-1);
	dist(0) = 0;
	for (int i = 1; i< m; i++)
	{
		dist(i) = ((du.array().square().rowwise().sum()).array().sqrt()).head(i).sum();
	}
	dist = dist.array()/dist(m-1);
//	std::cout << dist << std::endl;
}


void MyString::newstring()
{
	int ret;
	u_old = u;
	for( int i = 0; i < m; i += 1)
	{
		for( int j = 0; j < n; j += 1)
		{
			node.Anm[j] = u(i,j);
		}
		node.Energy = node.cal_F(node.Anm, node.Vnm, node.Qik);
		node.cal_dF(node.Anm, node.Vnm, node.Qik, node.grad_Energy);	
		
		std::cout << i << " " << node.Energy*2*PI << " " << node.Norm(node.grad_Energy, n) << "|"; 
	//	std::cout << "node=" << i << " "; 
		ret = node.run(n);
		for( int j = 0; j < n; j += 1)
		{
			u_new(i, j) = node.Anm[j];
		}
	}
	std::cout << std::endl;
}

void MyString::initialization()
{
	node.node_initialization();
	n = 5*node.Basis;
	u.setZero(m, n);
	u_new.setZero(m, n);
	u_old.setZero(m, n);
	node.readfile =1;node.mode = 0;

	for( int i = 0; i < m; i += 1)
	{
		//node.suffixe[0] = 'a';
	//std::cout << "~~~~" << std::endl;
		//node.suffixe[1] = char(i+1);
		node.suffixe.str("");
		node.suffixe << i+1;
		node.iput(64, 64);
		for( int j = 0; j < n; j += 1)
		{
			u(i, j) = node.Anm[j];
		}
	}
	//for(int  i = 0; i < n; i++)
	//{
		//u.col(i) = Eigen::VectorXd::LinSpaced(m*3, A1[i], A2[i]).segment(m,m);
	//}

//	std::cout << u.row(m-1) << std::endl;
}

double MyString::error()
{	
//	u_old = Eigen::node.MatrixXd::Zero(m, n*n);
	//for( int i = 0; i < m; i += 1)
	//{
		//for( int j = 0; j < n; j += 1)
		//{
			//if (u(i,j)-u_old(i,j)!=0){
			//std::cout << u(i, j)-u_old(i, j) << std::endl;
			//}	
		//}
		
	//}
	//std::cout << (u-u_old) << std::endl;
	double err = 0;
	//for( int i = 0; i < m; i += 1)
	//{
		//err = err + (u.row(i)-u_old.row(i)).norm();
	//}

	return (u-u_old).array().abs().maxCoeff();
	//Eigen::VectorXd err;
	//err.setZero(m);
	//for( int i = 0; i < m; i += 1)
	//{
		//err(i) = (u.row(i)-u_old.row(i)).norm();
	//}
	//return (err.maxCoeff() / h);
}

void MyString::result()
{
	//double Fbulk,Felas,Fpena,Fbeta,Energy,normdF;
	Eigen::VectorXd Ve(m);
	for( int i = 0; i < m; i += 1)
	{
		//node.suffixe[0] = 'a';
		//node.suffixe[1] = char(i+1);
		node.suffixe.str("");
		node.suffixe  << i+1; 
		for( int j = 0; j < n; j += 1)
		{
			node.Anm[j] = u_old(i, j);	
		}
		node.Energy = node.cal_F(node.Anm, node.Vnm, node.Qik);
		node.cal_dF(node.Anm, node.Vnm, node.Qik, node.grad_Energy);	
		fprintf(fp,"R = %.2f t = %.2f Energy = %16.15f normdF = %16.15f\n", node.Rad, node.landau_t,node.Energy*2*PI, node.Norm(node.grad_Energy,n));
		Ve(i) = node.Energy*2*PI;
		node.oput(node.N,node.M);
	}
		
	int num;
	Ve.maxCoeff(&num);
	
	std::cout << "Ve " << std::endl << Ve.transpose() << std::endl;
	std::cout << "saddle num " << std::endl << num+1 << std::endl;
}

void MyString::end()
{
  node.zer_destroy();
  node.var_destroy();
}
