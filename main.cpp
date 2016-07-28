#include "myString.h"

/*********************************************************************/
int main(int argc, const char *argv[]) 
{
	double Rad, t, eta, st_Rad, st_t, st_eta, ed_Rad, ed_t, ed_eta;

	std::stringstream filename("");
	if(argc <= 1) { std::cout << "no parameter 1" << std::endl; return 1; }
	filename << "./Energy/energy_" << argv[1] << ".txt";
	FILE *fp = fopen(filename.str().c_str(),"w");

	//double x1[5],x2[5];
	//double eg1[3],eg2[3];
	//double vec1[3][3],vec2[3][3];
	//int i11,i21,i31,i12,i22,i32;

	//st_eta = 1e5;
	//ed_eta = 1e5;
	st_Rad = 11;
	ed_Rad = 11;
	st_t = 1;
	ed_t = 1;
	MyString S;
	// for (eta = st_eta;eta <= ed_eta;eta = eta*1.1) {
	for (Rad = st_Rad; Rad <= ed_Rad; Rad += 2) 
	{
		for (t = st_t; t >= ed_t; t -= 1) 
		{
			S.initialization(Rad, t);
			double diff1 = S.diff_tol1 + 1;
			int nstep = 0;
			for (int k = 1; k < 200; k++)
				//while (diff1 >= S.diff_tol1)
			{
				S.newstring();
				S.distance();
				S.interp1();
				diff1 = S.error();
				nstep = nstep + 1;
				//	if (nstep%10 == 0)
				{
					std::cout << "nstep " << nstep << "\t diff1 " << diff1 << std::endl; 
				}
			}
			S.result(fp);
			fprintf(fp,"nstep = %d diff1 = %16.15f", nstep, diff1);
			fflush(fp);
		}
	}
	//}
	S.end(); 
	fclose(fp);
	return 0;
}


