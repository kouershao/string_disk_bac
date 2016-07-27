#ifndef _STRING_H_
#define _STRING_H_

#include <eigen3/Eigen/Dense>
#include "node.h"

class MyString
{
	public:
		MyString(int arg1 = 20, double arg2 = 0.01, double arg3 = 1e-8):
			m(arg1), diff_tol1(arg2), h(arg3) { }
		double diff_tol1;
		int n;
		MyNode node;
		
		void initialization();
		void newstring();
		void distance();
		void interp1();
		void result();
		void end();
		double error();
	
	private:
		FILE *fp;
		const int m;
		double h, saddle_error;
		Eigen::MatrixXd u, u_new, u_old;
		Eigen::VectorXd dist;
};

#endif
