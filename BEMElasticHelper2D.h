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
 this header file defines the basic computation of elements
 
 Triangle_Coord, Trans_Tri: transform of coordinates
 Serendip_func: Shape functions
 Normal_Jac: Jacobian of the elements
 Cartesian: Transform back to Cartesian
 UK, TK: Green's function for surface tranctions and displacements
 
 */
#pragma once
# include "GlobalConstant.h"
class BEM_2D {
public:
    void Cartesian(double** Elcor, double& N1, double& N2, double& N3, double* GCcor);
    void Serendip_func(double xsi, double& N1, double& N2, double& N3, double** Ni, double* N);
    void Normal_Jac(double xsi, double** Elcor, double& Jac, double* Vnorm);
    double Compute_length(double xsi, double weight, double** Elcor);

    void UK(double* dxr, double r, double E, double nu, double** UP);

    void TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm);

    // define the Green's function partial derivatives
    double Green_01(double* x, int i, int j, int k);
    double Green_02(double* x, int i, int j, int k, int l);
    double Green_03(double* x, int i, int j, int k, int l, int s);
    double Green_04(double* x, int i, int j, int k, int l, int s, int m);

};
