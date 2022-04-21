/* ----------------------------------------------------------------------
   iBEM - inclusion-based Boundary Element Method
   This software is distributed under the GNU General Public License
   version 3 or any later version.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
    This file is part of iBEM.

    iBEM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    iBEM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with iBEM.  If not, see <https://www.gnu.org/licenses/>.

------------------------------------------------------------------------- */


/*
 This file contains the integration of Green's function of
 Elliptical particles
 
 */

#pragma once
#include "GlobalConstant.h"
#include "Eigen/Dense"

extern double mu_0;
extern double nu0;
extern double mu_1;
extern double nu1;

using namespace Eigen;

//////////////////////////////////////////////////////////

// derivative of rho and direction normal
inline double rho_der(double *x,double a, double r, int j, int flag);

inline double rho_der2(double* x, double a, double r, int j, int m, int flag);

inline double rho_der3(double* x, double a, double r, int i, int j, int k);

inline double rho_der4(double* x, double a, double r, int i, int j, int k, int l);

inline double rho_der5(double* x, double a, double r, int i, int j, int k, int l, int s);

inline double dir_der(double* x, double a, double r, int p, int j);

inline double dir_der2(double* x, double a, double r, int i, int j, int m);

// integration for displacements

inline double Phi_1_1(double *x,double a, int p, int i);
inline double Phi_2_1(double* x, double a, int m, int n, int i);

inline double Psi_1_3(double* x, double a, int p, int i, int k, int l);

inline double Psi_2_3(double* x, double a, int m, int n, int i, int j, int k);

// integration for strain
inline double Phi_1_2(double* x, double a, int p, int i, int j);

inline double Phi_2_2(double* x, double a, int p,int q, int i, int j);

inline double Psi_1_4(double* x, double a, int p, int i, int j, int k, int l);

inline double Psi_2_4(double* x, double a, int p, int q, int i, int j, int k, int l);

// integration for strain der
inline double Phi_der2(double* x, double a, double r, int i, int j, int m);

inline double Phi_1_3(double* x, double a, int p, int i, int j, int k);

inline double Phi_2_3(double* x, double a, int p, int q, int i, int j, int k);

inline double Psi_der2(double* x, double a, double r, int i, int k, int l, int j, int m);

inline double Psi_1_5(double* x, double a, int p, int i, int j, int k, int l, int s);

inline double Psi_2_5(double* x, double a, int p, int q, int i, int j, int k, int l, int s);

// integration for strain der2
inline double Phi_der3(double* x, double a, double r, int i, int j, int k, int l);

inline double Phi_1_4(double* x, double a, int p, int i, int j, int k, int l);

inline double Phi_2_4(double* x, double a, int p, int q, int i, int j, int k, int l);

inline double Psi_der3(double* x, double a, double r, int i, int j, int k, int l, int s, int t);

inline double Psi_1_6(double* x, double a, int p, int i, int j, int k, int l, int s, int q);

inline double Psi_2_6(double* x, double a, int p, int q, int i, int j, int k, int l, int s, int t);

// define functions by V potential in Mura's book

inline double vxx_0(double* x, double a);
inline double vxx_1(double* x, double a, int i);
inline double vxx_2(double* x, double a, int i, int j);
inline double vxx_3(double* x, double a, int i, int j, int k);
inline double vxx_4(double* x, double a, int i, int j, int k, int l);
inline double vxxa_0(double* x, double a);


inline double vxxa_1(double* x, double a, int i);
inline double vxxa_2(double* x, double a, int i, int j);
inline double vxxa_3(double* x, double a, int i, int j, int k);
inline double vxxa_4(double* x, double a, int i, int j, int k, int l);
inline double vxxa_5(double* x, double a, int i, int j, int k, int l, int m);


inline double vx2axxx_0(double* x, double a);
inline double vx2axxx_1(double* x, double a, int i);
inline double vx2axxx_2(double* x, double a, int i, int j);
inline double vx2axxx_3(double* x, double a, int i, int j, int k);
inline double vx2axxx_4(double* x, double a, int i, int j, int k, int l);

inline double vxa_0(double* x, double a);
inline double vxa_1(double* x, double a, int i);
inline double vxa_2(double* x, double a, int i, int j);
inline double vxa_3(double* x, double a, int i, int j, int k);
inline double vxa_4(double* x, double a, int i, int j, int k, int l);
inline double vxa_5(double* x, double a, int i, int j, int k, int l, int s);

inline double vx_4(double* x, double a, int i, int j, int k, int l);
inline double vx_5(double* x, double a, int i, int j, int k, int l, int s);

// displacement
double Eshelby_displacement(double* x, double a, int i, int k, int l);
double Eshelby_displacement_linear(double* x, double a, int i, int k, int l, int p);
double Eshelby_displacement_quadratic(double* x, double a, int i, int k, int l, int p, int q);

// Eshelby tensor integration results by uniform strain
double Eshelby_strain_integration(int i, int j, int k, int l, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// Eshelby tensor integration results by linear strain
double Eshelby_strain_integration_linear(int i, int j, int k, int l, int p, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// Eshelby tensor integration results by quadratic strain
double Eshelby_strain_integration_quadratic(int i, int j, int k, int l, int p, int q, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// Eshelby tensor for der equivalent(uniform)
double Eshelby_strain_integration_der(int i, int j, int k, int l, int m, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// Eshelby tensor for der equivalent (linear)
double Eshelby_strain_integration_linear_der(int i, int j, int k, int l, int p,int m, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// Eshelby tensor for der equivalent (quadratic)
double Eshelby_strain_integration_quadratic_der(int i, int j, int k, int l, int p, int q, int m, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// Eshelby tensor for der2 equivalent(uniform)
double Eshelby_strain_integration_uniform_der2(int i, int j, int k, int l, int m, int s, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

double Eshelby_strain_integration_linear_der2(int i, int j, int k, int l, int p, int m,int s, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

double Eshelby_strain_integration_quadratic_der2(int i, int j, int k, int l, int p,int q, int m, int s, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s);

// constutive law
double cons_law(double nu, double mu, int i, int j, int k, int l);