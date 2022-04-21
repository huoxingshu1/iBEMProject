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

#include <cmath>
#include "GlobalConstant.h"
#include "configElastic3D_arbitrary.h"
#include "integratorhelper.h"
#include "integratorElastic3D_arbitrary.h"
#include "integratorElastic3DHelper.h"


/*
    Private members
*/
void integratorElastic3D_arbitrary::Polyhedral_Eshelby_displacement_all_10(int i, int j, int k, double** coord, double* x, double* result)
{
	double a, b, lp, lm, delta, para;
	double* Shape_0 = new double[4]; double* Shape_1 = new double[4];
	double* Shape_2 = new double[4]; double* Shape_3 = new double[4];
	double* Vnorm, * Vdir, * Svnorm, * vm, * vp; int om, op;
	Vnorm = new double[3]; Vdir = new double[3]; Svnorm = new double[3];
	vm = new double[3]; vp = new double[3];
	//double Vnorm[3], Vdir[3], Svnorm[3]; double vm[3], vp[3]; 
	para = 1.0 / (8.0 * Pi * (1.0 - nu0)); delta = 1E-7;
	// define the initial parameters
	// the sequence of the edge is not necessary to be assigned
	//int surface[4][3] = { {0,2,1},{1,2,3},{0,3,2},{0,1,3} };

	// provide the shape function for node "n"
	I3A.tet_shape(coord, 0, Shape_0); I3A.tet_shape(coord, 1, Shape_1); I3A.tet_shape(coord, 2, Shape_2); I3A.tet_shape(coord, 3, Shape_3);

	// when n = 0, 1, 2,3, shape function = Li * (2 * Li - 1)
	double* A_ijk = new double[10]; double* B_rest = new double[10];
	// define the shape functions in ten orders
	double** coeff = new double* [10];
	for (int h = 0; h < 10; h++)
		coeff[h] = new double[10];
	I3A.tet_shape_10_coeff_A(Shape_0, Shape_1, Shape_2, Shape_3, coeff);

	// define the summation parameters
	double Psi_ijk, Psi_0_ijk, Psi_1_ijk, Psi_2_ijk, Psi_01_ijk, Psi_02_ijk, Psi_12_ijk, Psi_00_ijk, Psi_11_ijk, Psi_22_ijk;
	Psi_ijk = Psi_0_ijk = Psi_1_ijk = Psi_2_ijk = Psi_01_ijk = Psi_02_ijk = Psi_12_ijk = Psi_00_ijk = Psi_11_ijk = Psi_22_ijk = 0.0;

	double fi_i, fi_0_i, fi_1_i, fi_2_i, fi_01_i, fi_02_i, fi_12_i, fi_00_i, fi_11_i, fi_22_i;
	fi_i = fi_0_i = fi_1_i = fi_2_i = fi_01_i = fi_02_i = fi_12_i = fi_00_i = fi_11_i = fi_22_i = 0.0;
	double fi_j, fi_0_j, fi_1_j, fi_2_j, fi_01_j, fi_02_j, fi_12_j, fi_00_j, fi_11_j, fi_22_j;
	fi_j = fi_0_j = fi_1_j = fi_2_j = fi_01_j = fi_02_j = fi_12_j = fi_00_j = fi_11_j = fi_22_j = 0.0;
	double fi_k, fi_0_k, fi_1_k, fi_2_k, fi_01_k, fi_02_k, fi_12_k, fi_00_k, fi_11_k, fi_22_k;
	fi_k = fi_0_k = fi_1_k = fi_2_k = fi_01_k = fi_02_k = fi_12_k = fi_00_k = fi_11_k = fi_22_k = 0.0;

	// given the two linear and quadratic potentials
	double* result_linear = new double[3]; double* result_quadratic = new double[6];
	for (int ppp = 0; ppp < 3; ppp++) {
		result_linear[ppp] = 0.0;
	}
	for (int ppp = 0; ppp < 6; ppp++) {
		result_quadratic[ppp] = 0.0;
	}


	for (int m = 0; m < 4; m++) {

		// calculate the surface normal and corresponding length a
		I3A.surface_normal(Svnorm, coord[surface[m][0]], coord[surface[m][1]], coord[surface[m][2]]);
		a = Svnorm[0] * (coord[surface[m][0]][0] - x[0]) + Svnorm[1] * (coord[surface[m][0]][1] - x[1]) + Svnorm[2] * (coord[surface[m][0]][2] - x[2]);


		if (abs(a) <= delta) a = 0.0;


		// start looping at each surface, 3 nodes at each surface
		for (int e = 0; e < 3; e++) {

			switch (e) {
			case 0:
				om = 1; op = 2;
				break;
			case 1:
				om = 2; op = 0;
				break;
			case 2:
				om = 0; op = 1;
				break;
			default:
				cout << "check the order of elements node numbering" << endl;
			}

			// assign the detailed coordinate to v+ and v- for the edge
			for (int g = 0; g < 3; g++) {
				vm[g] = coord[surface[m][om]][g];	vp[g] = coord[surface[m][op]][g];
			}

			// initialization of the specific TC geometry for b, l+ and l-
			b = lp = lm = 0.0;

			// assign the directional vector / normal vector of the e^{th} edge on m^{th} surface
			I3A.edge_direction_normal(Svnorm, Vnorm, Vdir, vm, vp);

			// assign the specific TC geometry
			for (int g = 0; g < 3; g++) {
				b = b + (vp[g] - x[g]) * Vnorm[g];
				lm = lm + (vm[g] - x[g]) * Vdir[g];
				lp = lp + (vp[g] - x[g]) * Vdir[g];
			}

			if (abs(b) <= delta) b = 0.0;
			if (abs(lm) <= delta) lm = 0.0;
			if (abs(lp) <= delta) lp = 0.0;

			// define the components in potentials
			double ij, ik, jk, kj, ijk, jk_5, jk_6;

			// uniform potentials
			Psi_ijk += U_P.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);
			fi_i += U_P.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir); fi_j += U_P.PHI_i(j, a, b, lm, lp, Svnorm, Vnorm, Vdir); fi_k += U_P.PHI_i(k, a, b, lm, lp, Svnorm, Vnorm, Vdir);

			// linear potentials in the combination form, the biharmonic
			L_P.PSI_p_ijk_combine(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x, result_linear, ij, ik, jk, kj, ijk, jk_5, jk_6);
			Psi_0_ijk += result_linear[0]; Psi_1_ijk += result_linear[1]; Psi_2_ijk += result_linear[2];

			// linear potentials, the harmonic
			//Psi_0_ijk += PSI_p_ijk(0, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_1_ijk += PSI_p_ijk(1, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_2_ijk += PSI_p_ijk(2, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_0_i += L_P.PHI_p_i(0, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_i += L_P.PHI_p_i(1, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_i += L_P.PHI_p_i(2, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_0_j += L_P.PHI_p_i(0, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_j += L_P.PHI_p_i(1, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_j += L_P.PHI_p_i(2, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_0_k += L_P.PHI_p_i(0, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_k += L_P.PHI_p_i(1, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_k += L_P.PHI_p_i(2, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

			// quadratic potentials, the biharmonic
			Q_P.PSI_pq_ijk_combine(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x, result_linear, result_quadratic, ij, ik, jk, kj, ijk, jk_5, jk_6);
			Psi_01_ijk += result_quadratic[0]; Psi_02_ijk += result_quadratic[1]; Psi_12_ijk += result_quadratic[2];
			Psi_00_ijk += result_quadratic[3]; Psi_11_ijk += result_quadratic[4]; Psi_22_ijk += result_quadratic[5];
			//Psi_01_ijk += PSI_pq_ijk(0, 1, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_02_ijk += PSI_pq_ijk(0, 2, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_12_ijk += PSI_pq_ijk(1, 2, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			//Psi_00_ijk += PSI_pq_ijk(0, 0, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_11_ijk += PSI_pq_ijk(1, 1, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_22_ijk += PSI_pq_ijk(2, 2, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

			//cout << Psi_01_ijk << endl;

			// quadratic potentials, the harmonic
			fi_01_i += Q_P.PHI_pq_i(0, 1, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_i += Q_P.PHI_pq_i(0, 2, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_i += Q_P.PHI_pq_i(1, 2, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_00_i += Q_P.PHI_pq_i(0, 0, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_i += Q_P.PHI_pq_i(1, 1, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_i += Q_P.PHI_pq_i(2, 2, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

			fi_01_j += Q_P.PHI_pq_i(0, 1, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_j += Q_P.PHI_pq_i(0, 2, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_j += Q_P.PHI_pq_i(1, 2, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_00_j += Q_P.PHI_pq_i(0, 0, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_j += Q_P.PHI_pq_i(1, 1, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_j += Q_P.PHI_pq_i(2, 2, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

			fi_01_k += Q_P.PHI_pq_i(0, 1, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_k += Q_P.PHI_pq_i(0, 2, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_k += Q_P.PHI_pq_i(1, 2, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_00_k += Q_P.PHI_pq_i(0, 0, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_k += Q_P.PHI_pq_i(1, 1, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_k += Q_P.PHI_pq_i(2, 2, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);


		} // end looping at each surface


	}
	// rewrite in forms of new double
	double* Psi_co = new double[10]; double* Phi_i_co = new double[10]; double* Phi_j_co = new double[10]; double* Phi_k_co = new double[10];
	Psi_co[0] = Psi_ijk; Psi_co[1] = Psi_0_ijk; Psi_co[2] = Psi_1_ijk; Psi_co[3] = Psi_2_ijk; Psi_co[4] = Psi_01_ijk;
	Psi_co[5] = Psi_02_ijk; Psi_co[6] = Psi_12_ijk; Psi_co[7] = Psi_00_ijk; Psi_co[8] = Psi_11_ijk; Psi_co[9] = Psi_22_ijk;

	Phi_i_co[0] = fi_i; Phi_i_co[1] = fi_0_i; Phi_i_co[2] = fi_1_i; Phi_i_co[3] = fi_2_i; Phi_i_co[4] = fi_01_i;
	Phi_i_co[5] = fi_02_i; Phi_i_co[6] = fi_12_i; Phi_i_co[7] = fi_00_i; Phi_i_co[8] = fi_11_i; Phi_i_co[9] = fi_22_i;

	Phi_j_co[0] = fi_j; Phi_j_co[1] = fi_0_j; Phi_j_co[2] = fi_1_j; Phi_j_co[3] = fi_2_j; Phi_j_co[4] = fi_01_j;
	Phi_j_co[5] = fi_02_j; Phi_j_co[6] = fi_12_j; Phi_j_co[7] = fi_00_j; Phi_j_co[8] = fi_11_j; Phi_j_co[9] = fi_22_j;

	Phi_k_co[0] = fi_k; Phi_k_co[1] = fi_0_k; Phi_k_co[2] = fi_1_k; Phi_k_co[3] = fi_2_k; Phi_k_co[4] = fi_01_k;
	Phi_k_co[5] = fi_02_k; Phi_k_co[6] = fi_12_k; Phi_k_co[7] = fi_00_k; Phi_k_co[8] = fi_11_k; Phi_k_co[9] = fi_22_k;

	//double Psi_co[10] = { Psi_ijk, Psi_0_ijk, Psi_1_ijk, Psi_2_ijk, Psi_01_ijk, Psi_02_ijk, Psi_12_ijk, Psi_00_ijk, Psi_11_ijk, Psi_22_ijk };
	//double Phi_i_co[10] = { fi_i, fi_0_i, fi_1_i, fi_2_i, fi_01_i, fi_02_i, fi_12_i, fi_00_i, fi_11_i, fi_22_i };
	//double Phi_j_co[10] = { fi_j, fi_0_j, fi_1_j, fi_2_j, fi_01_j, fi_02_j, fi_12_j, fi_00_j, fi_11_j, fi_22_j };
	//double Phi_k_co[10] = { fi_k, fi_0_k, fi_1_k, fi_2_k, fi_01_k, fi_02_k, fi_12_k, fi_00_k, fi_11_k, fi_22_k };

	for (int n = 0; n < 10; n++) {
		A_ijk[n] = 0.0; B_rest[n] = 0.0;
		for (int h = 0; h < 10; h++) {
			A_ijk[n] = A_ijk[n] + coeff[n][h] * Psi_co[h];

			B_rest[n] = B_rest[n] - 2.0 * nu0 * d[j][k] * Phi_i_co[h] * coeff[n][h]
				- 2.0 * (1.0 - nu0) * coeff[n][h] * (Phi_k_co[h] * d[j][i] + Phi_j_co[h] * d[i][k]);
		}

		result[n] = para * (A_ijk[n] + B_rest[n]);
	}

	delete[] Shape_0; delete[] Shape_1; delete[] Shape_2; delete[] Shape_3;
	delete[] A_ijk; delete[] B_rest; delete[] result_linear; delete[] result_quadratic;
	delete[] Psi_co; delete[] Phi_i_co; delete[] Phi_j_co; delete[] Phi_k_co;
	for (int s = 0; s < 10; s++) {
		delete coeff[s];
	}
	delete[] Svnorm; delete[] Vnorm; delete[] Vdir; delete[] vm; delete[] vp; delete[] coeff;
}

void integratorElastic3D_arbitrary::displacement_collect_normal(double** coord, double* x, int** index, double*** ret_u3)
{
	// initialization of the ret_u3 arrays for collecting the displacement field contribution of eigenstrain
	for (int i = 0; i < 10; i++) {
		for (int h = 0; h < 3; h++) {
			for (int m = 0; m < 6; m++) {
				ret_u3[i][h][m] = 0.0;
			}
		}
	}

	double* result = new double[10];

	// loop at the direction, then the m and n for summation
	for (int i = 0; i < 3; i++) {

		for (int m = 0; m < 3; m++) {

			for (int n = m; n < 3; n++) {

				// function of displacement 
				Polyhedral_Eshelby_displacement_all_10(i, m, n, coord, x, result);
				for (int s = 0; s < 10; s++) {
					ret_u3[s][i][index[m][n]] += map[m][n] * result[s];
				}

			}

		}


	}


	delete[] result;
}

void integratorElastic3D_arbitrary::Polyhedral_Eshelby_strain_all_10(int i, int j, int k, int l, double** coord, double* x, double* result)
{
	double a, b, lp, lm, delta, para;	double* Shape_0, * Shape_1, * Shape_2, * Shape_3;
	Shape_0 = new double[4]; Shape_1 = new double[4]; Shape_2 = new double[4]; Shape_3 = new double[4];

	double* Vnorm, * Vdir, * Svnorm, * vm, * vp;
	Vnorm = new double[3]; Vdir = new double[3]; Svnorm = new double[3]; vm = new double[3]; vp = new double[3];

	//double Vnorm[3], Vdir[3], Svnorm[3];
	//double vm[3], vp[3];
	int om, op;
	para = 1.0 / (8.0 * Pi * (1.0 - nu0)); delta = 1E-7;
	// define the initial parameters
	// the sequence of the edge is not necessary to be assigned
	//int surface[4][3] = { {0,2,1},{1,2,3},{0,3,2},{0,1,3} };

	// provide the shape function for node "n"
	I3A.tet_shape(coord, 0, Shape_0); I3A.tet_shape(coord, 1, Shape_1); I3A.tet_shape(coord, 2, Shape_2); I3A.tet_shape(coord, 3, Shape_3);

	// when n = 0, 1, 2,3, shape function = Li * (2 * Li - 1)
	double* A_ijkl = new double[10]; double* B_rest = new double[10];

	//double A_ijkl[10] = {}; double B_rest[10] = {};
	// define the shape functions in ten orders
	double** coeff = new double* [10];
	for (int h = 0; h < 10; h++) {
		coeff[h] = new double[10];
	}
	I3A.tet_shape_10_coeff_A(Shape_0, Shape_1, Shape_2, Shape_3, coeff);

	// components to be stored
	double Psi_klij, Psi_0_klij, Psi_1_klij, Psi_2_klij, Psi_01_klij, Psi_02_klij, Psi_12_klij, Psi_00_klij, Psi_11_klij, Psi_22_klij;
	Psi_klij = Psi_0_klij = Psi_1_klij = Psi_2_klij = Psi_01_klij = Psi_02_klij = Psi_12_klij = Psi_00_klij = Psi_11_klij = Psi_22_klij = 0.0;

	double fi_ij, fi_ik, fi_il, fi_jl, fi_jk;
	fi_ij = fi_ik = fi_il = fi_jk = fi_jl = 0.0;

	double fi_0_ij, fi_0_ik, fi_0_il, fi_0_jl, fi_0_jk;
	fi_0_ij = fi_0_ik = fi_0_il = fi_0_jl = fi_0_jk = 0.0;
	double fi_1_ij, fi_1_ik, fi_1_il, fi_1_jl, fi_1_jk;
	fi_1_ij = fi_1_ik = fi_1_il = fi_1_jl = fi_1_jk = 0.0;
	double fi_2_ij, fi_2_ik, fi_2_il, fi_2_jl, fi_2_jk;
	fi_2_ij = fi_2_ik = fi_2_il = fi_2_jl = fi_2_jk = 0.0;

	double fi_01_ij, fi_01_ik, fi_01_il, fi_01_jl, fi_01_jk;
	fi_01_ij = fi_01_ik = fi_01_il = fi_01_jl = fi_01_jk = 0.0;
	double fi_02_ij, fi_02_ik, fi_02_il, fi_02_jl, fi_02_jk;
	fi_02_ij = fi_02_ik = fi_02_il = fi_02_jl = fi_02_jk = 0.0;
	double fi_12_ij, fi_12_ik, fi_12_il, fi_12_jl, fi_12_jk;
	fi_12_ij = fi_12_ik = fi_12_il = fi_12_jl = fi_12_jk = 0.0;

	double fi_00_ij, fi_00_ik, fi_00_il, fi_00_jl, fi_00_jk;
	fi_00_ij = fi_00_ik = fi_00_il = fi_00_jl = fi_00_jk = 0.0;
	double fi_11_ij, fi_11_ik, fi_11_il, fi_11_jl, fi_11_jk;
	fi_11_ij = fi_11_ik = fi_11_il = fi_11_jl = fi_11_jk = 0.0;
	double fi_22_ij, fi_22_ik, fi_22_il, fi_22_jl, fi_22_jk;
	fi_22_ij = fi_22_ik = fi_22_il = fi_22_jl = fi_22_jk = 0.0;

	//double result_linear[3] = {}; double result_quadratic[6] = {};
	double* result_linear = new double[3]; double* result_quadratic = new double[6];
	// m -> 4, 4 surfaces
	for (int m = 0; m < 4; m++) {

		// calculate the surface normal and corresponding length a
		I3A.surface_normal(Svnorm, coord[surface[m][0]], coord[surface[m][1]], coord[surface[m][2]]);
		a = Svnorm[0] * (coord[surface[m][0]][0] - x[0]) + Svnorm[1] * (coord[surface[m][0]][1] - x[1]) + Svnorm[2] * (coord[surface[m][0]][2] - x[2]);


		if (abs(a) <= delta) a = 0.0;


		// start looping at each surface, 3 nodes at each surface
		for (int e = 0; e < 3; e++) {

			switch (e) {
			case 0:
				om = 1; op = 2;
				break;
			case 1:
				om = 2; op = 0;
				break;
			case 2:
				om = 0; op = 1;
				break;
			default:
				cout << "check the order of elements node numbering" << endl;
			}

			// assign the detailed coordinate to v+ and v- for the edge
			for (int g = 0; g < 3; g++) {
				vm[g] = coord[surface[m][om]][g];	vp[g] = coord[surface[m][op]][g];
			}

			// initialization of the specific TC geometry for b, l+ and l-
			b = lp = lm = 0.0;

			// assign the directional vector / normal vector of the e^{th} edge on m^{th} surface
			I3A.edge_direction_normal(Svnorm, Vnorm, Vdir, vm, vp);

			// assign the specific TC geometry
			for (int g = 0; g < 3; g++) {
				b = b + (vp[g] - x[g]) * Vnorm[g];
				lm = lm + (vm[g] - x[g]) * Vdir[g];
				lp = lp + (vp[g] - x[g]) * Vdir[g];
			}


			if (abs(b) <= delta) b = 0.0;
			if (abs(lm) <= delta) lm = 0.0;
			if (abs(lp) <= delta) lp = 0.0;

			int cal_flag = 1;

			// The calculation of Psi components, Psi_pq_klij
			if (a == 0 && b == 0) {
				if (lm == 0 && lp == 0) {

					cal_flag = 0;
				}
			}

			// write in a more expanding way

			//double Phi_0_a = Phi_u_II_a(a, b, lp) - Phi_u_II_a(a, b, lm);double Phi_0_b = Phi_u_II_b(a, b, lp) - Phi_u_II_b(a, b, lm);double Phi_0_le = Phi_u_II_le(a, b, lp) - Phi_u_II_le(a, b, lm);
			double ijk, ikl, ijl, jkl, kjl, ljk, ijkl, jkl_5, jkl_6;
			if (cal_flag != 0) {
				Psi_klij += U_P.PSI_ijkl(k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir);
				fi_ij += U_P.PHI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir); fi_ik += U_P.PHI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); fi_il += U_P.PHI_ij(i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); fi_jk += U_P.PHI_ij(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); fi_jl += U_P.PHI_ij(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
			}

			//Psi_0_klij += PSI_p_ijkl(0, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_1_klij += PSI_p_ijkl(1, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_2_klij += PSI_p_ijkl(2, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			//Psi_01_klij += PSI_pq_ijkl(0, 1, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_02_klij += PSI_pq_ijkl(0, 2, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_12_klij += PSI_pq_ijkl(1, 2, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			//Psi_00_klij += PSI_pq_ijkl(0, 0, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_11_klij += PSI_pq_ijkl(1, 1, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); Psi_22_klij += PSI_pq_ijkl(2, 2, k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			// linear potentials:
			L_P.PSI_p_ijkl_combine(k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x, result_linear, ijk, ikl, ijl, jkl, kjl, ljk, ijkl, jkl_5, jkl_6);
			Psi_0_klij += result_linear[0]; Psi_1_klij += result_linear[1]; Psi_2_klij += result_linear[2];
			// quadratic potentials:
			Q_P.PSI_pq_ijkl_combine(k, l, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x, result_linear, result_quadratic, ijk, ikl, ijl, jkl, kjl, ljk, ijkl, jkl_5, jkl_6);
			Psi_01_klij += result_quadratic[0]; Psi_02_klij += result_quadratic[1]; Psi_12_klij += result_quadratic[2];
			Psi_00_klij += result_quadratic[3]; Psi_11_klij += result_quadratic[4]; Psi_22_klij += result_quadratic[5];

			// the calculation of Phi components,

			fi_0_ij += L_P.PHI_p_ij(0, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_0_ik += L_P.PHI_p_ij(0, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_0_il += L_P.PHI_p_ij(0, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_0_jk += L_P.PHI_p_ij(0, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_0_jl += L_P.PHI_p_ij(0, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_1_ij += L_P.PHI_p_ij(1, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_ik += L_P.PHI_p_ij(1, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_il += L_P.PHI_p_ij(1, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_jk += L_P.PHI_p_ij(1, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_1_jl += L_P.PHI_p_ij(1, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_2_ij += L_P.PHI_p_ij(2, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_ik += L_P.PHI_p_ij(2, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_il += L_P.PHI_p_ij(2, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_jk += L_P.PHI_p_ij(2, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_2_jl += L_P.PHI_p_ij(2, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_01_ij += Q_P.PHI_pq_ij(0, 1, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_01_ik += Q_P.PHI_pq_ij(0, 1, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_01_il += Q_P.PHI_pq_ij(0, 1, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_01_jk += Q_P.PHI_pq_ij(0, 1, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_01_jl += Q_P.PHI_pq_ij(0, 1, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_02_ij += Q_P.PHI_pq_ij(0, 2, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_ik += Q_P.PHI_pq_ij(0, 2, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_il += Q_P.PHI_pq_ij(0, 2, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_jk += Q_P.PHI_pq_ij(0, 2, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_02_jl += Q_P.PHI_pq_ij(0, 2, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_12_ij += Q_P.PHI_pq_ij(1, 2, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_ik += Q_P.PHI_pq_ij(1, 2, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_il += Q_P.PHI_pq_ij(1, 2, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_jk += Q_P.PHI_pq_ij(1, 2, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_12_jl += Q_P.PHI_pq_ij(1, 2, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_00_ij += Q_P.PHI_pq_ij(0, 0, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_00_ik += Q_P.PHI_pq_ij(0, 0, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_00_il += Q_P.PHI_pq_ij(0, 0, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_00_jk += Q_P.PHI_pq_ij(0, 0, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_00_jl += Q_P.PHI_pq_ij(0, 0, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_11_ij += Q_P.PHI_pq_ij(1, 1, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_ik += Q_P.PHI_pq_ij(1, 1, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_il += Q_P.PHI_pq_ij(1, 1, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_jk += Q_P.PHI_pq_ij(1, 1, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_11_jl += Q_P.PHI_pq_ij(1, 1, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
			fi_22_ij += Q_P.PHI_pq_ij(2, 2, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_ik += Q_P.PHI_pq_ij(2, 2, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_il += Q_P.PHI_pq_ij(2, 2, i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_jk += Q_P.PHI_pq_ij(2, 2, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); fi_22_jl += Q_P.PHI_pq_ij(2, 2, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);


		} // end looping at each surface


	}

	// rewrite in forms of new double
	double* Psi_co = new double[10]; double* Phi_co_ij = new double[10]; double* Phi_co_ik = new double[10]; double* Phi_co_il = new double[10]; double* Phi_co_jk = new double[10]; double* Phi_co_jl = new double[10];

	Psi_co[0] = Psi_klij; Psi_co[1] = Psi_0_klij; Psi_co[2] = Psi_1_klij; Psi_co[3] = Psi_2_klij; Psi_co[4] = Psi_01_klij;
	Psi_co[5] = Psi_02_klij; Psi_co[6] = Psi_12_klij; Psi_co[7] = Psi_00_klij; Psi_co[8] = Psi_11_klij; Psi_co[9] = Psi_22_klij;

	Phi_co_ij[0] = fi_ij; Phi_co_ij[1] = fi_0_ij; Phi_co_ij[2] = fi_1_ij; Phi_co_ij[3] = fi_2_ij; Phi_co_ij[4] = fi_01_ij;
	Phi_co_ij[5] = fi_02_ij; Phi_co_ij[6] = fi_12_ij; Phi_co_ij[7] = fi_00_ij; Phi_co_ij[8] = fi_11_ij; Phi_co_ij[9] = fi_22_ij;

	Phi_co_ik[0] = fi_ik; Phi_co_ik[1] = fi_0_ik; Phi_co_ik[2] = fi_1_ik; Phi_co_ik[3] = fi_2_ik; Phi_co_ik[4] = fi_01_ik;
	Phi_co_ik[5] = fi_02_ik; Phi_co_ik[6] = fi_12_ik; Phi_co_ik[7] = fi_00_ik; Phi_co_ik[8] = fi_11_ik; Phi_co_ik[9] = fi_22_ik;

	Phi_co_il[0] = fi_il; Phi_co_il[1] = fi_0_il; Phi_co_il[2] = fi_1_il; Phi_co_il[3] = fi_2_il; Phi_co_il[4] = fi_01_il;
	Phi_co_il[5] = fi_02_il; Phi_co_il[6] = fi_12_il; Phi_co_il[7] = fi_00_il; Phi_co_il[8] = fi_11_il; Phi_co_il[9] = fi_22_il;

	Phi_co_jk[0] = fi_jk; Phi_co_jk[1] = fi_0_jk; Phi_co_jk[2] = fi_1_jk; Phi_co_jk[3] = fi_2_jk; Phi_co_jk[4] = fi_01_jk;
	Phi_co_jk[5] = fi_02_jk; Phi_co_jk[6] = fi_12_jk; Phi_co_jk[7] = fi_00_jk; Phi_co_jk[8] = fi_11_jk; Phi_co_jk[9] = fi_22_jk;

	Phi_co_jl[0] = fi_jl; Phi_co_jl[1] = fi_0_jl; Phi_co_jl[2] = fi_1_jl; Phi_co_jl[3] = fi_2_jl; Phi_co_jl[4] = fi_01_jl;
	Phi_co_jl[5] = fi_02_jl; Phi_co_jl[6] = fi_12_jl; Phi_co_jl[7] = fi_00_jl; Phi_co_jl[8] = fi_11_jl; Phi_co_jl[9] = fi_22_jl;

	for (int n = 0; n < 10; n++) {
		A_ijkl[n] = 0.0; B_rest[n] = 0.0;
		for (int h = 0; h < 10; h++) {
			A_ijkl[n] = A_ijkl[n] + coeff[n][h] * Psi_co[h];

			B_rest[n] = B_rest[n] - 2.0 * nu0 * d[k][l] * Phi_co_ij[h] * coeff[n][h]
				- (1.0 - nu0) * coeff[n][h] * (Phi_co_ik[h] * d[j][l] + Phi_co_il[h] * d[j][k] + Phi_co_jk[h] * d[i][l] + Phi_co_jl[h] * d[i][k]);
		}

		result[n] = para * (A_ijkl[n] + B_rest[n]);
	}

	delete[] Shape_0; delete[] Shape_1; delete[] Shape_2; delete[] Shape_3;
	delete[] Svnorm; delete[] Vnorm; delete[] Vdir; delete[] vm; delete[] vp;
	delete[] A_ijkl; delete[] B_rest; delete[] result_linear; delete[] result_quadratic;
	for (int s = 0; s < 10; s++) {
		delete coeff[s];
	}
	delete[] coeff;
	delete[] Psi_co; delete[] Phi_co_ij; delete[] Phi_co_ik; delete[] Phi_co_il; delete[] Phi_co_jk; delete[] Phi_co_jl;
}

void integratorElastic3D_arbitrary::Eshelby_collect_normal_quadratic(int i, int j, double** coord, double* x, int** index_E_ij, double*** ret_u4)
{
	int m, n, k, l, h;
	double* D_ijmn = new double[10]; double* result = new double[10];
	//double D_ijmn[10] = {};
	for (h = 0; h < 10; h++) {
		for (m = 0; m < 6; m++) {
			for (n = 0; n < 6; n++) {
				ret_u4[h][m][n] = 0.0;
			}
		}
	}


	for (m = 0; m < 3; m++) {

		for (n = m; n < 3; n++) {
			//double result[10] = {};
			for (int sqr = 0; sqr < 10; sqr++) {
				D_ijmn[sqr] = 0.0;
			}

			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					Polyhedral_Eshelby_strain_all_10(k, l, m, n, coord, x, result);
					for (int s = 0; s < 10; s++)
						D_ijmn[s] = D_ijmn[s] - (Constitutive(i, j, k, l) - Constitutive1(i, j, k, l)) * result[s];//Polyhedral_Eshelby_strain_all(k, l, m, n, h, coord, x);

				}


			}			// end of looping k
			for (h = 0; h < 10; h++)
				ret_u4[h][index_E_ij[i][j]][index_E_ij[m][n]] = ret_u4[h][index_E_ij[i][j]][index_E_ij[m][n]] + map[m][n] * D_ijmn[h];

		}


	}
	delete[] D_ijmn; delete[] result;

}

void integratorElastic3D_arbitrary::Eshleby_strain_polyhedra(int num_NN, int num_NE, double** A, Ref<MatrixXd> eigen_point, Ref<MatrixXi> eigen_connect, \
	int** index_E_ij, double*** changed_coord, double** number_weight, int* count_number)
{
	double**** record_weight = new double*** [num_NN];
	double*** returned_coord = new double** [num_NN];

	for (int t = 0; t < num_NN; t++) {
		record_weight[t] = new double** [number_in_array];
		returned_coord[t] = new double* [number_in_array];
	}

	for (int t = 0; t < num_NN; t++) {
		for (int m = 0; m < number_in_array; m++) {
			record_weight[t][m] = new double* [10];
			returned_coord[t][m] = new double[3];
		}
	}

	for (int t = 0; t < num_NN; t++)
		for (int m = 0; m < number_in_array; m++)
			for (int e = 0; e < 10; e++)
				record_weight[t][m][e] = new double[3];

	int* record_numbers = new int[num_NN];

# pragma omp parallel shared(A, record_weight, returned_coord, changed_coord, number_weight, count_number)

	{
		int t = 0;
		int m, n, p, r, h, k, i, j, s, y;
		int l, pp, ppp;
		double le, dist;
		// create the coordinate
		double** coord = new double* [10];
		double* x = new double[3];

		for (pp = 0; pp < 10; pp++) {
			coord[pp] = new double[3];
		}

		double** vertex_local = new double* [10];
		double*** ret_u4 = new double** [10];

		for (int sqr = 0; sqr < 10; sqr++) {
			ret_u4[sqr] = new double* [6];
			vertex_local[sqr] = new double[3];
		}

		for (int sqr = 0; sqr < 10; sqr++) {
			for (int sqq = 0; sqq < 6; sqq++) {
				ret_u4[sqr][sqq] = new double[6];
			}
		}

		////Assembly Matrix Tensor///////

#pragma omp for schedule(dynamic)

		for (s = 0; s < num_NN; s++) {

			x[0] = eigen_point(s, 0); x[1] = eigen_point(s, 1); x[2] = eigen_point(s, 2);
			int num_array = 0;

			I3A.gauss_10_point_shared(s, num_NE, eigen_point, eigen_connect, num_array, x, returned_coord, record_weight, changed_coord, number_weight, count_number);
			record_numbers[s] = num_array;


			for (int u = 0; u < num_array; u++) {

				x[0] = returned_coord[s][u][0]; x[1] = returned_coord[s][u][1]; x[2] = returned_coord[s][u][2];
				for (i = 0; i < 3; i++) {
					for (j = i; j < 3; j++) {

						for (pp = 0; pp < num_NE; pp++) {
							// assgin the coordinate of the nodes in the element
							for (int sqr = 0; sqr < 10; sqr++) {
								for (int sqq = 0; sqq < 3; sqq++)
									coord[sqr][sqq] = eigen_point(eigen_connect(pp, sqr) - 1, sqq);
							}
							// adding the cutoff distance, this keeps the same since the 4 points are the corner ones
							le = pow(6.0 * I3A.Volume_tet(coord), 1.0 / 3.0);
							dist = sqrt((x[0] - 0.25 * (coord[0][0] + coord[1][0] + coord[2][0] + coord[3][0])) * (x[0] - 0.25 * (coord[0][0] + coord[1][0] + coord[2][0] + coord[3][0]))
								+ (x[1] - 0.25 * (coord[0][1] + coord[1][1] + coord[2][1] + coord[3][1])) * (x[1] - 0.25 * (coord[0][1] + coord[1][1] + coord[2][1] + coord[3][1]))
								+ (x[2] - 0.25 * (coord[0][2] + coord[1][2] + coord[2][2] + coord[3][2])) * (x[2] - 0.25 * (coord[0][2] + coord[1][2] + coord[2][2] + coord[3][2])));
							if (dist >= 10.0 * le)	continue;


							Eshelby_collect_normal_quadratic(i, j, coord, x, index_E_ij, ret_u4);

#pragma omp critical
							{
								for (h = 0; h < 10; h++) {
									for (m = 0; m < 3; m++) {

										for (n = m; n < 3; n++) {

											A[index_E_ij[3 * s + i][j]][index_E_ij[3 * (eigen_connect(pp, h) - 1) + m][n]] += ret_u4[h][index_E_ij[i][j]][index_E_ij[m][n]] * record_weight[s][u][0][2];

										}


									}
								}



							}





						}


					}
				}
			}

			cout << "num = " << " " << s;
		}

		// delete and reset the coord
		for (int nnnn = 0; nnnn < 10; nnnn++) {
			delete coord[nnnn];
			delete vertex_local[nnnn];
			for (int ssss = 0; ssss < 6; ssss++) {
				delete ret_u4[nnnn][ssss];
			}
			delete ret_u4[nnnn];
		}

		delete[] vertex_local;
		delete[]coord;
		delete[] ret_u4;
		delete[] x;

	}

	int s, i, j, m, n;
	for (s = 0; s < num_NN; s++) {

		for (i = 0; i < 3; i++) {
			for (j = i; j < 3; j++) {
				for (m = 0; m < 3; m++) {
					for (n = m; n < 3; n++) {

						for (int y = 0; y < 10; y++) {
							for (int u = 0; u < record_numbers[s]; u++) {
								A[index_E_ij[3 * s + i][j]][index_E_ij[3 * int(record_weight[s][u][y][0]) + m][n]] += map[m][n] * Constitutive(i, j, m, n) * record_weight[s][u][y][1] * record_weight[s][u][y][2];
								
							}
						}

					}
				}
			}
		}
	}


	for (int nnnn = 0; nnnn < num_NN; nnnn++) {
		for (int ssss = 0; ssss < number_in_array; ssss++) {
			for (int qqqq = 0; qqqq < 10; qqqq++)
				delete record_weight[nnnn][ssss][qqqq];
			delete record_weight[nnnn][ssss];
			delete returned_coord[nnnn][ssss];
		}
		delete record_weight[nnnn];
		delete returned_coord[nnnn];
	}
	delete[] record_weight;
	delete[] returned_coord;



}

/*
    Public members
*/
void integratorElastic3D_arbitrary::addFieldToBEM(Config& config_)
{
    configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

    int& num_NN = config.num_NN; int& num_NE = config.num_NE;
    int& NN = config.NN;
    Ref<MatrixXd> HMAT = config.HMAT;

    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> eigen_point = config.eigen_point;
    Ref<MatrixXi> eigen_connect = config.eigen_connect;

# pragma omp parallel shared (HMAT) 
	{
		double* x = new double[3];
		double** coord = new double* [10];
		for (int ii = 0; ii < 10; ii++) {
			coord[ii] = new double[3];
		}

		int** index; int id = 0; int pp;

		index = new int* [3];
		for (int ii = 0; ii < 3; ii++) {
			index[ii] = new int[3];
		}


		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				index[ii][jj] = 0;
			}
		}

		for (int ii = 0; ii < 10; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				coord[ii][jj] = 0.0;
			}
		}

		// numbering the indices for id

		for (int ii = 0; ii < 3; ii++) {
			for (int jj = ii; jj < 3; jj++) {
				index[ii][jj] = id;
				id = id + 1;
			}
		}

		// create the array to store the computation of displacement 
		double*** ret_u3 = new double** [10];
		for (int ii = 0; ii < 10; ii++) {
			ret_u3[ii] = new double* [3];
		}
		for (int ii = 0; ii < 10; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				ret_u3[ii][jj] = new double[6];
			}
		}

		////////////////////////////////////////////////////////////////////      finish the initialization

#pragma omp for schedule(dynamic)

		for (int MM = 0; MM < num_NE; MM++) {

			// assign the coordinates of each element to 4 nodes, the MM^{th} element
			for (int ii = 0; ii < 10; ii++) {	// 10 nodes of tet
				for (int jj = 0; jj < 3; jj++) {	// 3D problem
					coord[ii][jj] = eigen_point(eigen_connect(MM, ii) - 1, jj);
				}
			}

			for (int jj = 0; jj < NN; jj++) {
				x[0] = NODES(jj, 0);
				x[1] = NODES(jj, 1);
				x[2] = NODES(jj, 2);

				// get_addu
				//get_addu_tensor_polyhedral(x, coord, DD_u3, index, kk);

				pp = 3 * jj;
				displacement_collect_normal(coord, x, index, ret_u3);

#pragma omp critical
				{
					for (int h = 0; h < 10; h++) {
						for (int qq = pp; qq < pp + 3; qq++) {

							for (int mm = 0; mm < 6; mm++) {
								HMAT(qq, 3 * int(NN) + mm + 6 * int(eigen_connect(MM, h) - 1)) += -ret_u3[h][qq - pp][mm];

							}

						}
					}
				}




			}




		}	// finishe looping


		for (int ii = 0; ii < 10; ii++) {
			for (int jj = 0; jj < 3; jj++) {
				delete ret_u3[ii][jj];
			}
			delete ret_u3[ii];
		}
		delete[] ret_u3;

		for (int ii = 0; ii < 3; ii++)
		{
			delete[] index[ii];

		}

		delete[] x;
		delete[] index;


		for (int i = 0; i < 10; i++) {
			delete[] coord[i];
		}
		delete[]coord;

	}


}

void integratorElastic3D_arbitrary::addFluxEquivalentInclusion(Config& config_)
{
	configElastic3D_arbitrary& config = dynamic_cast<configElastic3D_arbitrary&> (config_);

	int& num_NN = config.num_NN; int& num_NE = config.num_NE;

	int*& count_number = config.count_number;
	double**& number_weight = config.number_weight;
	double***& changed_coord = config.changed_coord;

	int& NN = config.NN;
	Ref<MatrixXd> HMAT = config.HMAT;
	int**& index_E_ij = config.index_E_ij;
	Ref<MatrixXd> NODES = config.NODES;
	Ref<MatrixXd> eigen_point = config.eigen_point;
	Ref<MatrixXi> eigen_connect = config.eigen_connect;

	double** A;
	const int row = 6 * num_NN;

	A = new double* [row];

	for (int i = 0; i < row; i++) {
		A[i] = new double[row];
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < row; j++) {
			A[i][j] = 0.0;
		}
	}

	// finish the initilization
	Eshleby_strain_polyhedra(num_NN, num_NE, A, eigen_point, eigen_connect, index_E_ij, changed_coord, number_weight, count_number);

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < row; j++)
		{
			HMAT(3 * NN + i, 3 * NN + j) = A[i][j];
		}
	}


	for (int i = 0; i < row; i++) {
		delete A[i];
	}
	delete[] A;


}

void integratorElastic3D_arbitrary::addFluxEquivalentFirstOrderInclusion(Config& config_) {}

void integratorElastic3D_arbitrary::addFluxEquivalentSecondOrderInclusion(Config& config_) {}
