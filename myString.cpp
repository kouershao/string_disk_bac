#include "myString.h"
#include <cmath>
#include <thread>
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
//	u = u_new;
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

	std::vector<std::thread> ths;
	for(int i = 0; i < m; i++)
		ths.push_back(std::thread(&MyString::inneriter, this, i));
	for(int i = 0; i < m; i++)
		ths[i].join();

	//thread.join
	Eigen::VectorXd a;
	a.setZero(n);
	for (int i = 1; i < m-1; i++)
	{
		a = u_new.row(i) - u_old.row(i);
		a = perp(a, i);
		u_new.row(i) = u_old.row(i) + a.transpose();
	}
	std::cout << std::endl;
}

int MyString::inneriter(int i)
{
	int ret = 0;
	for(int j = 0; j < n; j += 1)
		nodes[i].Anm[j] = u(i, j);
	if (i >= 1 && i < m-1){
	 ret = nodes[i].run(i, n, u, u_new, gradient, Ve);}
	for(int j = 0; j < n; j += 1)
		u_new(i, j) = nodes[i].Anm[j];
	return ret;
}

void MyString::initialization(double Rad, double t)
{
	nodes = std::vector<MyNode>(m);
	std::for_each(nodes.begin(), nodes.end(), [Rad, t](MyNode &node) {
		node.node_initialization();
		node.readfile = 1;
		node.mode = 0;
		node.Rad = Rad;
		node.landau_t = t;
	});

	n = 5 * nodes[0].Basis;
	u.setZero(m, n);
	u_new.setZero(m, n);
	u_old.setZero(m, n);
	gradient.setZero(m, n);
	Ve.setZero(m);
	for(int i = 0; i < m; i += 1)
	{
		nodes[i].suffix.str("");
		nodes[i].suffix << "newini" << i+1;
		nodes[i].iput(64, 64);
		for(int j = 0; j < n; j += 1)
		{
			u(i, j) = nodes[i].Anm[j];
		}
	}
	//u_old = u;
	//for(int  i = 0; i < n; i++)
	//{
	//    u.col(i) = Eigen::VectorXd::LinSpaced(m*3, A1[i], A2[i]).segment(m,m);
	//}

	//std::cout << u.row(m-1) << std::endl;
}
Eigen::VectorXd MyString::perp(Eigen::VectorXd vector, int i)
{
	Eigen::VectorXd tao;
	tao.setZero(n);
	if (Ve(i+1)>Ve(i)&&Ve(i)>Ve(i-1))
	{tao = u_old.row(i+1)-u_old.row(i);}
	else if (Ve(i+1)<Ve(i)&&Ve(i)<Ve(i-1))
	{tao = u_old.row(i)-u_old.row(i-1);}
	else
	{tao = u_old.row(i+1)-u_old.row(i-1);}
	tao = tao/(tao.norm());
	std::cout  << i << " angle " << 180/PI*acos(tao.dot(vector)/vector.norm()) << " | ";
	tao = (vector-tao.dot(vector)*tao);
	//std::cout << " tao.norm() " << tao.norm() << std::endl;
	return(tao);
}
int MyString::error(FILE *fp)
{	
	Eigen::VectorXd perpgrad;
	double sgrad = 0, sperp = 0;
	int stopnow = 0;
	std::cout << "node_error" << std::endl;
	for( int i = 1; i < m-1; i += 1)
	{
		perpgrad.setZero(n);
		perpgrad = perp(gradient.row(i), i);
		//std::cout << i << " " << (gradient.row(i).norm()) << " " << perpgrad.norm() << " | ";
		sgrad = sgrad + gradient.row(i).norm();
		sperp = sperp + perpgrad.norm();
	}
	std::cout << std::endl;
	double err = 0;
	for( int i = 0; i < m; i += 1)
	{
		err = err + (u.row(i)-u_old.row(i)).norm();
	}
	fprintf(fp,"normdF = %16.15f  normdF_normal = %16.15f  err = %16.15f\n", sgrad, sperp, err);
	fflush(fp);
	if (sgrad/10>sperp || err<1e-6) 
	{
		stopnow = 1;
	}
	return stopnow;
}
void MyString::result(FILE* fp)
{
	//double Fbulk,Felas,Fpena,Fbeta,Energy,normdF;
	Eigen::VectorXd Ve(m);
	for(int i = 0; i < m; i++)
	{
		nodes[i].suffix.str("");
		nodes[i].suffix  << "newinifinal" << i + 1; 
		for( int j = 0; j < n; j += 1)
		{
			nodes[i].Anm[j] = u_old(i, j);	
		}
		nodes[i].Energy = nodes[i].cal_F(nodes[i].Anm, nodes[i].Vnm, nodes[i].Qik);
		nodes[i].cal_dF(nodes[i].Anm, nodes[i].Vnm, nodes[i].Qik, nodes[i].grad_Energy);	
		fprintf(fp,"node = %d R = %.2f t = %.2f Energy = %16.15f normdF = %16.15f\n", i, nodes[i].Rad, nodes[i].landau_t,nodes[i].Energy*2*PI, nodes[i].Norm(nodes[i].grad_Energy,n));

		Ve(i) = nodes[i].Energy*2*PI;
		nodes[i].oput(nodes[i].N,nodes[i].M);
	}

	int num;
	Ve.maxCoeff(&num);

	std::cout << "Ve " << std::endl << Ve.transpose() << std::endl;
	std::cout << "saddle num " << std::endl << num+1 << std::endl;
}

void MyString::end()
{
	std::for_each(nodes.begin(), nodes.end(), [](MyNode &node) {
		node.zer_destroy();
		node.var_destroy();
	});
}
