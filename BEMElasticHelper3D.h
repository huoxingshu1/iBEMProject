#pragma once
# include "GlobalConstant.h"
void Triangle_Coord(int ntr, int l_pi, int** cor_tri);
void Trans_Tri(int ntr, int l_pi, double xsib, double etab, double& xsi, double& eta, double& Jacb);
void Serendip_func(double xsi, double eta, double& N1, double& N2, double& N3, double& N4, double* N, double** Ni);
void Normal_Jac(double xsi, double eta, double** Elcor, double& Jac, double* Vnorm);
void Cartesian(double** Elcor, double& N1, double& N2, double& N3, double& N4, double* GCcor);
void UK(double* dxr, double r, double E, double nu, double** UP);
void TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm);
