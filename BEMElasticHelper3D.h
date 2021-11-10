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

#pragma once
# include "GlobalConstant.h"
void Triangle_Coord(int ntr, int l_pi, int** cor_tri);
void Trans_Tri(int ntr, int l_pi, double xsib, double etab, double& xsi, double& eta, double& Jacb);
void Serendip_func(double xsi, double eta, double& N1, double& N2, double& N3, double& N4, double* N, double** Ni);
void Normal_Jac(double xsi, double eta, double** Elcor, double& Jac, double* Vnorm);
void Cartesian(double** Elcor, double& N1, double& N2, double& N3, double& N4, double* GCcor);
void UK(double* dxr, double r, double E, double nu, double** UP);
void TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm);
