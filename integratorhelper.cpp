
#include "integratorhelper.h"
#include "GlobalConstant.h"
#include <cmath>

using namespace Eigen;

double integratorhelper3D_arbitrary::Volume_tet(double** coord)
{
	MatrixXd AA = MatrixXd::Zero(4, 4);

	for (int i = 0; i < 4; i++) {
		AA(i, 0) = 1.0; AA(i, 1) = coord[i][0]; AA(i, 2) = coord[i][1]; AA(i, 3) = coord[i][2];
	}

	double result = AA.determinant();

	return (1.0 / 6.0) * abs(result);

	AA.resize(0, 0);
}

void integratorhelper3D_arbitrary::tet_shape(double** coord, int n, double* SH)
{
	// n is used to show which one is the I node
	MatrixXd AA = MatrixXd::Zero(4, 4);

	for (int i = 0; i < 4; i++) {
		AA(i, 0) = 1.0; AA(i, 1) = coord[i][0]; AA(i, 2) = coord[i][1]; AA(i, 3) = coord[i][2];
	}

	// inverse of the matrix creating the shape function
	MatrixXd inv_AA = AA.inverse();

	SH[0] = inv_AA(0, n);			// alpha
	SH[1] = inv_AA(1, n);			// beta
	SH[2] = inv_AA(2, n);			// gamma	
	SH[3] = inv_AA(3, n);			// chi

	AA.resize(0, 0); inv_AA.resize(0, 0);

}

double integratorhelper3D_arbitrary::tet_shape_10(double** coord, int n, double* x)
{
	// n is used to show which one is the I node
	double SH[4] = {}; double L;
	switch (n) {
	case 0:
		tet_shape(coord, n, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		L = (2.0 * L - 1.0) * L;
		break;
	case 1:
		tet_shape(coord, n, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		L = (2.0 * L - 1.0) * L;
		break;
	case 2:
		tet_shape(coord, n, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		L = (2.0 * L - 1.0) * L;
		break;

	case 3:
		tet_shape(coord, n, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		L = (2.0 * L - 1.0) * L;
		break;

	case 4:
		tet_shape(coord, 0, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		tet_shape(coord, 1, SH);
		L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2]);
		break;

	case 5:
		tet_shape(coord, 1, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		tet_shape(coord, 2, SH);
		L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2]);
		break;
	case 6:
		tet_shape(coord, 2, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		tet_shape(coord, 0, SH);
		L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2]);
		break;
	case 7:
		tet_shape(coord, 0, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		tet_shape(coord, 3, SH);
		L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2]);
		break;
	case 8:
		tet_shape(coord, 1, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		tet_shape(coord, 3, SH);
		L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2]);
		break;
	case 9:
		tet_shape(coord, 2, SH);
		L = SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2];
		tet_shape(coord, 3, SH);
		L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1] + SH[3] * x[2]);
		break;
	}

	return L;
}

/*---------------------------------------------------------------
	Geometry specifiction of normal and directional vector
-----------------------------------------------------------------*/

void integratorhelper3D_arbitrary::surface_normal(double* Svnorm, double* x1, double* x2, double* x3)
{
	// define the vector of the edge in counterclockwise direction
	double u1 = x2[0] - x1[0];
	double u2 = x2[1] - x1[1];
	double u3 = x2[2] - x1[2];

	double v1 = x3[0] - x2[0];
	double v2 = x3[1] - x2[1];
	double v3 = x3[2] - x2[2];

	// compute the cross product 
	double a1 = u2 * v3 - u3 * v2;
	double a2 = v1 * u3 - u1 * v3;
	double a3 = u1 * v2 - u2 * v1;

	double mod = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
	Svnorm[0] = a1 / mod; Svnorm[1] = a2 / mod; Svnorm[2] = a3 / mod;
}

void integratorhelper3D_arbitrary::edge_direction_normal(double* Svnorm, double* Vnorm, double* Vdir, double* x1, double* x2)
{
	double len = sqrt((x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[2] - x2[2]) * (x1[2] - x2[2]));

	// x2 = v+, x1 = v-
	Vdir[0] = (x2[0] - x1[0]) / len;
	Vdir[1] = (x2[1] - x1[1]) / len;
	Vdir[2] = (x2[2] - x1[2]) / len;

	// the normal vector = vdir \times Svnorm, as we want the outward direction
	Vnorm[0] = (Vdir[1] * Svnorm[2] - Vdir[2] * Svnorm[1]);
	Vnorm[1] = (Svnorm[0] * Vdir[2] - Vdir[0] * Svnorm[2]);
	Vnorm[2] = (Vdir[0] * Svnorm[1] - Vdir[1] * Svnorm[0]);

}

void integratorhelper3D_arbitrary::tet_shape_10_coeff_A(double* SH_0, double* SH_1, double* SH_2, double* SH_3, double** coeff)
{
	for (int i = 0; i < 10; i++)
		for (int h = 0; h < 10; h++)
			coeff[i][h] = 0.0;

	// define the first 4, vertex, 0
	coeff[0][0] = 2.0 * SH_0[0] * SH_0[0] - SH_0[0]; coeff[0][1] = 4.0 * (SH_0[0] * SH_0[1]) - SH_0[1]; coeff[0][2] = 4.0 * (SH_0[0] * SH_0[2]) - SH_0[2]; coeff[0][3] = 4.0 * (SH_0[0] * SH_0[3]) - SH_0[3];
	coeff[0][4] = 4.0 * SH_0[1] * SH_0[2]; coeff[0][5] = 4.0 * SH_0[1] * SH_0[3]; coeff[0][6] = 4.0 * SH_0[2] * SH_0[3]; coeff[0][7] = 2.0 * SH_0[1] * SH_0[1];
	coeff[0][8] = 2.0 * SH_0[2] * SH_0[2]; coeff[0][9] = 2.0 * SH_0[3] * SH_0[3];
	// vertex, 1
	coeff[1][0] = 2.0 * SH_1[0] * SH_1[0] - SH_1[0]; coeff[1][1] = 4.0 * (SH_1[0] * SH_1[1]) - SH_1[1]; coeff[1][2] = 4.0 * (SH_1[0] * SH_1[2]) - SH_1[2]; coeff[1][3] = 4.0 * (SH_1[0] * SH_1[3]) - SH_1[3];
	coeff[1][4] = 4.0 * SH_1[1] * SH_1[2]; coeff[1][5] = 4.0 * SH_1[1] * SH_1[3]; coeff[1][6] = 4.0 * SH_1[2] * SH_1[3]; coeff[1][7] = 2.0 * SH_1[1] * SH_1[1];
	coeff[1][8] = 2.0 * SH_1[2] * SH_1[2]; coeff[1][9] = 2.0 * SH_1[3] * SH_1[3];
	// vertex, 2
	coeff[2][0] = 2.0 * SH_2[0] * SH_2[0] - SH_2[0]; coeff[2][1] = 4.0 * (SH_2[0] * SH_2[1]) - SH_2[1]; coeff[2][2] = 4.0 * (SH_2[0] * SH_2[2]) - SH_2[2]; coeff[2][3] = 4.0 * (SH_2[0] * SH_2[3]) - SH_2[3];
	coeff[2][4] = 4.0 * SH_2[1] * SH_2[2]; coeff[2][5] = 4.0 * SH_2[1] * SH_2[3]; coeff[2][6] = 4.0 * SH_2[2] * SH_2[3]; coeff[2][7] = 2.0 * SH_2[1] * SH_2[1];
	coeff[2][8] = 2.0 * SH_2[2] * SH_2[2]; coeff[2][9] = 2.0 * SH_2[3] * SH_2[3];
	// vertex, 3
	coeff[3][0] = 2.0 * SH_3[0] * SH_3[0] - SH_3[0]; coeff[3][1] = 4.0 * (SH_3[0] * SH_3[1]) - SH_3[1]; coeff[3][2] = 4.0 * (SH_3[0] * SH_3[2]) - SH_3[2]; coeff[3][3] = 4.0 * (SH_3[0] * SH_3[3]) - SH_3[3];
	coeff[3][4] = 4.0 * SH_3[1] * SH_3[2]; coeff[3][5] = 4.0 * SH_3[1] * SH_3[3]; coeff[3][6] = 4.0 * SH_3[2] * SH_3[3]; coeff[3][7] = 2.0 * SH_3[1] * SH_3[1];
	coeff[3][8] = 2.0 * SH_3[2] * SH_3[2]; coeff[3][9] = 2.0 * SH_3[3] * SH_3[3];

	// the mid-points are assigned as 4 * Li * Lj
	// vertex 5, between point 2 and 3, 4 * L0 * L1
	coeff[4][0] = 4.0 * SH_0[0] * SH_1[0]; coeff[4][1] = 4.0 * (SH_0[0] * SH_1[1] + SH_1[0] * SH_0[1]); coeff[4][2] = 4.0 * (SH_0[0] * SH_1[2] + SH_1[0] * SH_0[2]); coeff[4][3] = 4.0 * (SH_0[0] * SH_1[3] + SH_1[0] * SH_0[3]);
	coeff[4][4] = 4.0 * (SH_0[1] * SH_1[2] + SH_1[1] * SH_0[2]); coeff[4][5] = 4.0 * (SH_0[1] * SH_1[3] + SH_1[1] * SH_0[3]); coeff[4][6] = 4.0 * (SH_0[2] * SH_1[3] + SH_1[2] * SH_0[3]); coeff[4][7] = 4.0 * SH_0[1] * SH_1[1];
	coeff[4][8] = 4.0 * SH_0[2] * SH_1[2]; coeff[4][9] = 4.0 * SH_0[3] * SH_1[3];

	// vertex 6, between point 1 and 3, 4 * L1 * L2
	coeff[5][0] = 4.0 * SH_1[0] * SH_2[0]; coeff[5][1] = 4.0 * (SH_1[0] * SH_2[1] + SH_2[0] * SH_1[1]); coeff[5][2] = 4.0 * (SH_1[0] * SH_2[2] + SH_2[0] * SH_1[2]); coeff[5][3] = 4.0 * (SH_1[0] * SH_2[3] + SH_2[0] * SH_1[3]);
	coeff[5][4] = 4.0 * (SH_1[1] * SH_2[2] + SH_2[1] * SH_1[2]); coeff[5][5] = 4.0 * (SH_1[1] * SH_2[3] + SH_2[1] * SH_1[3]); coeff[5][6] = 4.0 * (SH_1[2] * SH_2[3] + SH_2[2] * SH_1[3]); coeff[5][7] = 4.0 * SH_1[1] * SH_2[1];
	coeff[5][8] = 4.0 * SH_1[2] * SH_2[2]; coeff[5][9] = 4.0 * SH_1[3] * SH_2[3];

	// vertex 7, between point 1 and 2, 4 * L0 * L2
	coeff[6][0] = 4.0 * SH_0[0] * SH_2[0]; coeff[6][1] = 4.0 * (SH_0[0] * SH_2[1] + SH_2[0] * SH_0[1]); coeff[6][2] = 4.0 * (SH_0[0] * SH_2[2] + SH_2[0] * SH_0[2]); coeff[6][3] = 4.0 * (SH_0[0] * SH_2[3] + SH_2[0] * SH_0[3]);
	coeff[6][4] = 4.0 * (SH_0[1] * SH_2[2] + SH_2[1] * SH_0[2]); coeff[6][5] = 4.0 * (SH_0[1] * SH_2[3] + SH_2[1] * SH_0[3]); coeff[6][6] = 4.0 * (SH_0[2] * SH_2[3] + SH_2[2] * SH_0[3]); coeff[6][7] = 4.0 * SH_0[1] * SH_2[1];
	coeff[6][8] = 4.0 * SH_0[2] * SH_2[2]; coeff[6][9] = 4.0 * SH_0[3] * SH_2[3];

	// vertex 8, between point 1 and 4, 4 * L0 * L3
	coeff[7][0] = 4.0 * SH_0[0] * SH_3[0]; coeff[7][1] = 4.0 * (SH_0[0] * SH_3[1] + SH_3[0] * SH_0[1]); coeff[7][2] = 4.0 * (SH_0[0] * SH_3[2] + SH_3[0] * SH_0[2]); coeff[7][3] = 4.0 * (SH_0[0] * SH_3[3] + SH_3[0] * SH_0[3]);
	coeff[7][4] = 4.0 * (SH_0[1] * SH_3[2] + SH_3[1] * SH_0[2]); coeff[7][5] = 4.0 * (SH_0[1] * SH_3[3] + SH_3[1] * SH_0[3]); coeff[7][6] = 4.0 * (SH_0[2] * SH_3[3] + SH_3[2] * SH_0[3]); coeff[7][7] = 4.0 * SH_0[1] * SH_3[1];
	coeff[7][8] = 4.0 * SH_0[2] * SH_3[2]; coeff[7][9] = 4.0 * SH_0[3] * SH_3[3];

	// vertex 9, between point 2 and 4, 4 * L1 * L3
	coeff[8][0] = 4.0 * SH_1[0] * SH_3[0]; coeff[8][1] = 4.0 * (SH_1[0] * SH_3[1] + SH_3[0] * SH_1[1]); coeff[8][2] = 4.0 * (SH_1[0] * SH_3[2] + SH_3[0] * SH_1[2]); coeff[8][3] = 4.0 * (SH_1[0] * SH_3[3] + SH_3[0] * SH_1[3]);
	coeff[8][4] = 4.0 * (SH_1[1] * SH_3[2] + SH_3[1] * SH_1[2]); coeff[8][5] = 4.0 * (SH_1[1] * SH_3[3] + SH_3[1] * SH_1[3]); coeff[8][6] = 4.0 * (SH_1[2] * SH_3[3] + SH_3[2] * SH_1[3]); coeff[8][7] = 4.0 * SH_1[1] * SH_3[1];
	coeff[8][8] = 4.0 * SH_1[2] * SH_3[2]; coeff[8][9] = 4.0 * SH_1[3] * SH_3[3];

	// vertex 10, between point 3 and 4, 4 * L2 * L3
	coeff[9][0] = 4.0 * SH_2[0] * SH_3[0]; coeff[9][1] = 4.0 * (SH_2[0] * SH_3[1] + SH_3[0] * SH_2[1]); coeff[9][2] = 4.0 * (SH_2[0] * SH_3[2] + SH_3[0] * SH_2[2]); coeff[9][3] = 4.0 * (SH_2[0] * SH_3[3] + SH_3[0] * SH_2[3]);
	coeff[9][4] = 4.0 * (SH_2[1] * SH_3[2] + SH_3[1] * SH_2[2]); coeff[9][5] = 4.0 * (SH_2[1] * SH_3[3] + SH_3[1] * SH_2[3]); coeff[9][6] = 4.0 * (SH_2[2] * SH_3[3] + SH_3[2] * SH_2[3]); coeff[9][7] = 4.0 * SH_2[1] * SH_3[1];
	coeff[9][8] = 4.0 * SH_2[2] * SH_3[2]; coeff[9][9] = 4.0 * SH_2[3] * SH_3[3];

}

void integratorhelper3D_arbitrary::gauss_10_point_shared(int s, int num_NE, Ref<MatrixXd> eigen_point, Ref<MatrixXi> eigen_connect, \
	int& num_array, double* x, double*** returned_coord, double**** record_weight, double*** changed_coord, double** number_weight, int* count_number)
{
	// this is to judge how many elements share the same nodes judge_flag!
	num_array = 0; int* record = new int[number_in_array];
	for (int i = 0; i < num_NE; i++) {
		for (int j = 0; j < 10; j++) {
			if (s == eigen_connect(i, j) - 1) {
				record[num_array] = i;
				num_array++;
			}
		}
	}

	count_number[s] = num_array;

	double* SH = new double[10]; double* GS = new double[3]; double** vertex = new double* [10]; double dist, min;

	for (int i = 0; i < 10; i++) {
		vertex[i] = new double[3];
	}

	// given the 10 points shape function, construct the most near equations at the given nodes
	for (int n = 0; n < num_array; n++) {

		// specify the coordinates of n^th elements
		for (int q = 0; q < 10; q++) {
			for (int t = 0; t < 3; t++) {
				vertex[q][t] = eigen_point(eigen_connect(record[n], q) - 1, t);
			}
		}

		min = 500;

		for (int i = 0; i < 10; i++) {

			// initilize thr GS coordinate
			GS[0] = GS[1] = GS[2] = 0.0;

			// corner points
			SH[0] = (2.0 * Shape_L1_3D[i] - 1.0) * Shape_L1_3D[i]; SH[1] = (2.0 * Shape_L2_3D[i] - 1.0) * Shape_L2_3D[i];
			SH[2] = (2.0 * Shape_L3_3D[i] - 1.0) * Shape_L3_3D[i]; SH[3] = (2.0 * Shape_L4_3D[i] - 1.0) * Shape_L4_3D[i];
			// mid_points
			SH[4] = 4.0 * Shape_L1_3D[i] * Shape_L2_3D[i]; SH[5] = 4.0 * Shape_L2_3D[i] * Shape_L3_3D[i]; SH[6] = 4.0 * Shape_L3_3D[i] * Shape_L1_3D[i];
			SH[7] = 4.0 * Shape_L1_3D[i] * Shape_L4_3D[i]; SH[8] = 4.0 * Shape_L2_3D[i] * Shape_L4_3D[i]; SH[9] = 4.0 * Shape_L3_3D[i] * Shape_L4_3D[i];

			// given the specific coordinates of Gauss point
			for (int q = 0; q < 10; q++) {
				for (int j = 0; j < 3; j++) {
					GS[j] = GS[j] + SH[q] * vertex[q][j];
				}
			}

			// examine the distance to locate the most close Gauss point
			dist = sqrt((x[0] - GS[0]) * (x[0] - GS[0]) + (x[1] - GS[1]) * (x[1] - GS[1]) + (x[2] - GS[2]) * (x[2] - GS[2]));
			if (dist < min) {
				// renew the coordinate, later in the returned_coord;
				min = dist;
				returned_coord[s][n][0] = GS[0]; returned_coord[s][n][1] = GS[1]; returned_coord[s][n][2] = GS[2];
				changed_coord[s][n][0] = GS[0]; changed_coord[s][n][1] = GS[1]; changed_coord[s][n][2] = GS[2];
				// renew the weight of Shpae function and Weight of the total equations set up
				for (int q = 0; q < 10; q++) {
					record_weight[s][n][q][0] = eigen_connect(record[n], q) - 1;
					record_weight[s][n][q][1] = SH[q];
					record_weight[s][n][q][2] = Weight_3D[i];
				}
				number_weight[s][n] = Weight_3D[i];
			}


		} // end of looping poinsts of 10


	}


	delete[] SH; delete[]GS;
	for (int i = 0; i < 10; i++) {
		delete vertex[i];
	}
	delete[] vertex;


}

double integratorhelper3D_arbitrary::func_surface(double* ko, double* Sv, double* x)
{
	// later, check the results of four surface in some signs or not
	double result = Sv[0] * (x[0] - ko[0]) + Sv[1] * (x[1] - ko[1]) + Sv[2] * (x[2] - ko[2]);
	return result;
}

int integratorhelper3D_arbitrary::judge_tet_point(double* x, double** coord)
{
	int flag;

	// create the surface nodes combination
	int surface[4][3] = { {0,2,1},{1,2,3},{0,3,2},{0,1,3} };
	// this assumes the surface vertices already existing in the counter-clockwise order
	// thus, the unit normal vector is towards outside!

	// calculate the surface normal vector
	double Sv1[3], Sv2[3], Sv3[3], Sv4[3];

	surface_normal(Sv1, coord[surface[0][0]], coord[surface[0][1]], coord[surface[0][2]]);
	surface_normal(Sv2, coord[surface[1][0]], coord[surface[1][1]], coord[surface[1][2]]);
	surface_normal(Sv3, coord[surface[2][0]], coord[surface[2][1]], coord[surface[2][2]]);
	surface_normal(Sv4, coord[surface[3][0]], coord[surface[3][1]], coord[surface[3][2]]);

	// function of the surfaces
	double p1, p2, p3, p4;
	p1 = func_surface(coord[surface[0][1]], Sv1, x);
	p2 = func_surface(coord[surface[1][1]], Sv2, x);
	p3 = func_surface(coord[surface[2][1]], Sv3, x);
	p4 = func_surface(coord[surface[3][1]], Sv4, x);

	// they have to be in the some sign excepting zero!
	// due to the low number, out is impossible!
	if (p1 <= 1E-6 && p2 <= 1E-6 && p3 <= 1E-6 && p4 <= 1E-6) {
		flag = 0;
	}
	else
		flag = 1;  // outside

	// even sometimes, p1 - p4 could be zero. they cannot all be zero!


	return flag;
}