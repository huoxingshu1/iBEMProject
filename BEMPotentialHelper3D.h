#pragma once
#include "GlobalConstant.h"
#include "configPotential3D_elliptical.h"

namespace  BEMPotential3D {

	void Triangle_Coord(int ntr, int l_pi, int** cor_tri);

	void Trans_Tri(int ntr, int l_pi, double xsib, double etab, double& xsi, double& eta, double& Jacb);

	void Serendip_func(double xsi, double eta, double& N1, double& N2, double& N3, double& N4, double* N, double** Ni);

	void Normal_Jac(double xsi, double eta, double** Elcor, double& Jac, double* Vnorm);

	void Cartesian(double** Elcor, double& N1, double& N2, double& N3, double& N4, double* GCcor);

	double R(double* x);

	double R_01(double* x, int i);

	double R_02(double* x, int i, int j);

	double R_03(double* x, int i, int j, int k);

	double fR(double* x, int sign);

	double fR_01(double* x, int i, int sign);

	double fR_02(double* x, int i, int j, int sign);

	double fR_03(double* x, int i, int j, int k, int sign);

	double Green_01(double* x, int i);

	double Green_02(double* x, int i, int j);

	double Green_03(double* x, int i, int j, int k);

	void UK(double* x, double r, double& UP);

	void TK(double* x, double r, double& TP, double* Vnorm);

}