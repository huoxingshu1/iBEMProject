/*
 this header file defines the basic computation of elements
 
 Triangle_Coord, Trans_Tri: transform of coordinates
 Serendip_func: Shape functions
 Normal_Jac: Jacobian of the elements
 Cartesian: Transform back to Cartesian
 UK, TK: Green's function for surface tranctions and displacements
 
 */
#pragma once
# include "GlobalConstant.h"
inline void  Cartesian(double** Elcor, double& N1, double& N2, double& N3, double* GCcor);
inline void Serendip_func(double xsi, double& N1, double& N2, double& N3, double** Ni, double* N);
inline void Normal_Jac(double xsi, double** Elcor, double& Jac, double* Vnorm);
inline double Compute_length(double xsi, double weight, double** Elcor);

inline void UK(double* dxr, double r, double E, double nu, double** UP);

inline void TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm);
