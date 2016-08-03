#ifndef _MY_STRING_H_
#define _MY_STRING_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include "myNode.h"

class MyString
{
	public:
		MyString(int arg1 = 20, double arg2 = 1e-7, double arg3 = 1e-8):
			m(arg1), diff_tol1(arg2), h(arg3) { }
		double diff_tol1;
		int n;
		std::vector<MyNode> nodes;

		void initialization(double, double);
		void newstring();
		void distance();
		void interp1();
		void result(FILE*);
		void end();
		double error();

	private:
		FILE *fp;
		const int m;
		double h, saddle_error;
		Eigen::MatrixXd u, u_new, u_old, gradient;
		Eigen::VectorXd dist;
		int inneriter(int i);
};

#endif
