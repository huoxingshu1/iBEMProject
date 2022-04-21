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
#include "configElastic2D_arbitrary.h"
#include "Uniform_polygon_integral.h"
#include "Linear_Polygon_integral.h"
#include "Quadratic_Polygon_integral.h"
#include "integratorElastic2D_arbitrary.h"
#include "integratorElastic2DHelper.h"
/*
    Private members
*/



void integratorElastic2D_arbitrary::Polyhedral_Eshelby_displacement_all_6(int i, int j, int k, double** coord, double* x, double* result)
{
    double bI, lp, lm, delta, para;
    double Shape_0[3] = {}; double Shape_1[3] = {}; double Shape_2[3] = {};
    double Vnorm[2], Vdir[2]; double vm[2], vp[2]; int om, op;
    para = 1.0 / (8.0 * Pi * (1.0 - nu0)); delta = 1E-7;
    double* A_ijk = new double[6]; double* B_rest = new double[6];

    B2D_A.Triangle_shape(coord, 0, Shape_0);  B2D_A.Triangle_shape(coord, 1, Shape_1);  B2D_A.Triangle_shape(coord, 2, Shape_2);
    double** coeff = new double* [6];
    // Uniform_polygon Uni; Linear_polygon Lin; Quadratic_polygon Qua;
    for (int h = 0; h < 6; h++) {
        coeff[h] = new double[6];
    }

    B2D_A.Trangle_shape_6_coeff_6(Shape_0, Shape_1, Shape_2, coeff);

    // components to be stored
    double Psi_ijk, Psi_0_ijk, Psi_1_ijk, Psi_01_ijk, Psi_00_ijk, Psi_11_ijk;
    Psi_ijk = Psi_0_ijk = Psi_1_ijk = Psi_01_ijk = Psi_00_ijk = Psi_11_ijk = 0.0;

    double fi_i, fi_0_i, fi_1_i, fi_01_i, fi_00_i, fi_11_i;
    fi_i = fi_0_i = fi_1_i = fi_01_i = fi_00_i = fi_11_i = 0.0;
    double fj_j, fj_0_j, fj_1_j, fj_01_j, fj_00_j, fj_11_j;
    fj_j = fj_0_j = fj_1_j = fj_01_j = fj_00_j = fj_11_j = 0.0;
    double fk_k, fk_0_k, fk_1_k, fk_01_k, fk_00_k, fk_11_k;
    fk_k = fk_0_k = fk_1_k = fk_01_k = fk_00_k = fk_11_k = 0.0;

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

        for (int g = 0; g < 2; g++) {
            vm[g] = coord[om][g]; vp[g] = coord[op][g];
        }

        bI = lm = lp = 0.0;

        B2D_A.direction_normal(Vnorm, Vdir, vm, vp);

        for (int pp = 0; pp < 2; pp++) {
            bI += (vp[pp] - x[pp]) * Vnorm[pp];
            lm += (vm[pp] - x[pp]) * Vdir[pp];
            lp += (vp[pp] - x[pp]) * Vdir[pp];
        }                                       // compute the integral contour parameters

        if (abs(bI) < delta) bI = 0;
        if (abs(lm) < delta) lm = 0;
        if (abs(lp) < delta) lp = 0;

        // uniform potentials:
        Psi_ijk += U_P.PSI_3(i, j, k, bI, lm, lp, Vnorm, Vdir);
        fi_i += U_P.PHI_1(i, bI, lm, lp, Vnorm, Vdir); fj_j += U_P.PHI_1(j, bI, lm, lp, Vnorm, Vdir); fk_k += U_P.PHI_1(k, bI, lm, lp, Vnorm, Vdir);

        // linear potentials:
        Psi_0_ijk += L_P.PSI_p_3(0, i, j, k, bI, lm, lp, Vnorm, Vdir, x); Psi_1_ijk += L_P.PSI_p_3(1, i, j, k, bI, lm, lp, Vnorm, Vdir, x);
        fi_0_i += L_P.PHI_p_1(0, i, bI, lm, lp, Vnorm, Vdir, x); fj_0_j += L_P.PHI_p_1(0, j, bI, lm, lp, Vnorm, Vdir, x); fk_0_k += L_P.PHI_p_1(0, k, bI, lm, lp, Vnorm, Vdir, x);
        fi_1_i += L_P.PHI_p_1(1, i, bI, lm, lp, Vnorm, Vdir, x); fj_1_j += L_P.PHI_p_1(1, j, bI, lm, lp, Vnorm, Vdir, x); fk_1_k += L_P.PHI_p_1(1, k, bI, lm, lp, Vnorm, Vdir, x);


        // Quadratic potentials:
        Psi_01_ijk += Q_P.PSI_pq_3(0, 1, i, j, k, bI, lm, lp, Vnorm, Vdir, x); Psi_00_ijk += Q_P.PSI_pq_3(0, 0, i, j, k, bI, lm, lp, Vnorm, Vdir, x); Psi_11_ijk += Q_P.PSI_pq_3(1, 1, i, j, k, bI, lm, lp, Vnorm, Vdir, x);
        fi_01_i += Q_P.PHI_pq_1(0, 1, i, bI, lm, lp, Vnorm, Vdir, x); fi_00_i += Q_P.PHI_pq_1(0, 0, i, bI, lm, lp, Vnorm, Vdir, x); fi_11_i += Q_P.PHI_pq_1(1, 1, i, bI, lm, lp, Vnorm, Vdir, x);
        fj_01_j += Q_P.PHI_pq_1(0, 1, j, bI, lm, lp, Vnorm, Vdir, x); fj_00_j += Q_P.PHI_pq_1(0, 0, j, bI, lm, lp, Vnorm, Vdir, x); fj_11_j += Q_P.PHI_pq_1(1, 1, j, bI, lm, lp, Vnorm, Vdir, x);
        fk_01_k += Q_P.PHI_pq_1(0, 1, k, bI, lm, lp, Vnorm, Vdir, x); fk_00_k += Q_P.PHI_pq_1(0, 0, k, bI, lm, lp, Vnorm, Vdir, x); fk_11_k += Q_P.PHI_pq_1(1, 1, k, bI, lm, lp, Vnorm, Vdir, x);

    }

    double Psi_co[6] = { Psi_ijk, Psi_0_ijk, Psi_1_ijk, Psi_01_ijk, Psi_00_ijk, Psi_11_ijk };
    double Phi_i[6] = { fi_i, fi_0_i, fi_1_i, fi_01_i, fi_00_i, fi_11_i };
    double Phi_j[6] = { fj_j, fj_0_j, fj_1_j, fj_01_j, fj_00_j, fj_11_j };
    double Phi_k[6] = { fk_k, fk_0_k, fk_1_k, fk_01_k, fk_00_k, fk_11_k };


    for (int n = 0; n < 6; n++) {
        A_ijk[n] = 0.0; B_rest[n] = 0.0;
        for (int h = 0; h < 6; h++) {
            A_ijk[n] = A_ijk[n] + coeff[n][h] * Psi_co[h];

            B_rest[n] = B_rest[n] - (2.0 * nu0 * d[j][k] * Phi_i[h] + 2.0 * (1.0 - nu0) * (d[i][k] * Phi_j[h] + d[i][j] * Phi_k[h])) * coeff[n][h];
        }

        result[n] = para * (A_ijk[n] + B_rest[n]);
    }


    delete[] A_ijk; delete[] B_rest;

}

void integratorElastic2D_arbitrary::displacement_collect_normal(double** coord, double* x, int** index, double*** ret_u3)
{
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 2; j++) {
            for (int m = 0; m < 3; m++) {
                ret_u3[i][j][m] = 0.0;
            }
        }
    }

    // in the 6 elements
    for (int i = 0; i < 2; i++) {
        double result[6] = { };

        // third order tensor, imn

        for (int m = 0; m < 2; m++) {
            for (int n = m; n < 2; n++) {

                Polyhedral_Eshelby_displacement_all_6(i, m, n, coord, x, result);

                for (int s = 0; s < 6; s++) {
                    ret_u3[s][i][index[m][n]] = ret_u3[s][i][index[m][n]] + map[m][n] * result[s];
                }
            }
        }

    }

}

void integratorElastic2D_arbitrary::Eshelby_strain_polygons(int num_NN, int num_NE, double** A, Ref<MatrixXd> eigen_point, \
    Ref<MatrixXi> eigen_connect, int** index_E_ij, double*** changed_coord, double** number_weight, int* count_number)
{
    double**** record_weight = new double*** [num_NN];
    double*** returned_coord = new double** [num_NN];

    for (int t = 0; t < num_NN; t++) {
        record_weight[t] = new double** [number_in_array];
        returned_coord[t] = new double* [number_in_array];
    }

    for (int t = 0; t < num_NN; t++) {
        for (int m = 0; m < number_in_array; m++) {
            record_weight[t][m] = new double* [6];
            returned_coord[t][m] = new double[2];
        }
    }

    for (int t = 0; t < num_NN; t++) {
        for (int m = 0; m < number_in_array; m++) {
            for (int e = 0; e < 6; e++) {
                record_weight[t][m][e] = new double[3];
            }
        }
    }

    int* record_numbers = new int[num_NN];

# pragma omp parallel shared(A, record_numbers, record_weight, returned_coord, count_number, number_weight, changed_coord)
    {
        int t = 0; int m, n, p, r, h, k, i, j, s, l, y;
        double** coord = new double* [6];
        double* x = new double[2];
        for (int pp = 0; pp < 6; pp++) {
            coord[pp] = new double[2];
        }

        double*** ret_u4 = new double** [6];
        double*** ret_u3 = new double** [6];

        for (int sqr = 0; sqr < 6; sqr++) {
            ret_u4[sqr] = new double* [3];
            ret_u3[sqr] = new double* [3];
        }
        for (int sqq = 0; sqq < 6; sqq++)
            for (int sqr = 0; sqr < 3; sqr++) {
                ret_u4[sqq][sqr] = new double[3];
                ret_u3[sqq][sqr] = new double[2];
            }

        /// Assemble the matrix ///
        ////////////////////////////////////////////////////////////////////
# pragma omp for schedule(dynamic)
        for (s = 0; s < num_NN; s++) {

            x[0] = eigen_point(s, 0); x[1] = eigen_point(s, 1);

            int num_array = 0;

            gauss_6_point_shared(s, num_NE, eigen_point, eigen_connect, num_array, x, returned_coord, record_weight, changed_coord, number_weight, count_number);
            record_numbers[s] = num_array;
            for (int u = 0; u < num_array; u++) {
                x[0] = returned_coord[s][u][0]; x[1] = returned_coord[s][u][1];
                for (i = 0; i < 2; i++) {

                    for (j = i; j < 2; j++) {
                        for (int pp = 0; pp < num_NE; pp++) {

                            // assign the coordinate of the triangles
                            for (int sqr = 0; sqr < 6; sqr++) {
                                for (int sqq = 0; sqq < 2; sqq++)
                                    coord[sqr][sqq] = eigen_point(eigen_connect(pp, sqr) - 1, sqq);
                            }

                            // setting up the cutoff distance as 15 \times characteristic distance
                            double dist = sqrt(2.0 * B2D_A.Triangle_area(coord));
                            double center[2] = { (coord[0][0] + coord[1][0] + coord[2][0]) / 3.0, (coord[0][1] + coord[1][1] + coord[2][1]) / 3.0 };
                            double ddst = sqrt((center[0] - x[0]) * (center[0] - x[0]) + (center[1] - x[1]) * (center[1] - x[1]));
                            //if (ddst >= 10.0 * dist) {
                            //    continue;
                            //}


                            // continue;
                            Eshelby_collect_normal(i, j, coord, x, index_E_ij, ret_u4);



#pragma omp critical 

                            {
                                for (h = 0; h < 6; h++) {
                                    for (m = 0; m < 2; m++) {


                                        for (n = m; n < 2; n++) {

                                            A[index_E_ij[2 * s + i][j]][index_E_ij[2 * (eigen_connect(pp, h) - 1) + m][n]] += ret_u4[h][index_E_ij[i][j]][index_E_ij[m][n]] * record_weight[s][u][0][2];

                                        }

                                      


                                    }
                                }


                            }




                        }       // end of the pp

                    }


                }
            }
            cout << "num = " << " " << s;
        }
        for (int pp = 0; pp < 6; pp++) {
            delete coord[pp];
        }
        delete[] coord;
        for (int nnnn = 0; nnnn < 6; nnnn++) {
            for (int qqqq = 0; qqqq < 3; qqqq++) {
                delete ret_u4[nnnn][qqqq];
                delete ret_u3[nnnn][qqqq];
            }
            delete ret_u4[nnnn];
            delete ret_u3[nnnn];
        }

        delete[] ret_u4;
        delete[] ret_u3;

        delete[] x;
    }


    // avoid overlapping summation of the sym index.
    for (int s = 0; s < num_NN; s++) {
        for (int i = 0; i < 2; i++) {
            for (int j = i; j < 2; j++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = m; n < 2; n++) {
                        for (int y = 0; y < 6; y++) {
                            for (int u = 0; u < record_numbers[s]; u++) {
                                //A[index_E_ij[2 * s + i][j]][index_E_ij[2 * s + m][n]] += map[m][n] * cons_law(nu0, mu0, i, j, m, n);
                                A[index_E_ij[2 * s + i][j]][index_E_ij[2 * int(record_weight[s][u][y][0]) + m][n]] += map[m][n] * cons_law(nu0, mu_0, i, j, m, n) * record_weight[s][u][y][1] * record_weight[s][u][y][2];
                              

                            }
                        }
                    }
                }
            }
        }
    }

    /*
        Uncomment the following block to obtain output information of changed coordinates
    */
    /*
    ofstream myfile; myfile.open("changed_coord.txt"); ofstream myfile_2; myfile_2.open("weight_coord.txt"); ofstream myfile_3; myfile_3.open("surface_coord.txt");
    for (int nnnn = 0; nnnn < num_NN; nnnn++) {
        for (int u = 0; u < record_numbers[nnnn]; u++) {
            int judd = -1;
            myfile << "original coord" << " " << eigen_point(nnnn, 0) << " " << eigen_point(nnnn, 1) << endl;

            myfile << "changed coord" << " " << returned_coord[nnnn][u][0] << " " << returned_coord[nnnn][u][1] << endl;

            myfile_2 << "node number" << " ";
            for (int qqqq = 0; qqqq < 6; qqqq++) {
                myfile_2 << record_weight[nnnn][u][qqqq][0] << " ";
            }
            myfile_2 << endl; myfile_2 << "weights" << " ";

            for (int qqqq = 0; qqqq < 6; qqqq++) {
                myfile_2 << record_weight[nnnn][u][qqqq][1] << " ";
            }
            myfile_2 << endl;

            for (int qqqq = 0; qqqq < 6; qqqq++) {
                if (nnnn == record_weight[nnnn][u][qqqq][0])
                    judd = 1;
            }

            if (judd == -1)
                myfile_2 << "Not included!" << endl;
        }
    }

    for (int nnnn = 0; nnnn < num_NN; nnnn++) {
        myfile_3 << "number of shared elements" << record_numbers[nnnn] << endl;
    }
    */

    for (int nnnn = 0; nnnn < num_NN; nnnn++) {
        for (int ssss = 0; ssss < number_in_array; ssss++) {
            for (int qqqq = 0; qqqq < 6; qqqq++)
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

void integratorElastic2D_arbitrary::gauss_6_point_shared(int s, int num_NE, Ref<MatrixXd> eigen_point, Ref<MatrixXi> eigen_connect, int& num_array, \
    double* x, double*** returned_coord, double**** record_weight, double*** changed_coord, double** number_weight, int* count_number)
{

    num_array = 0; int* record = new int[number_in_array];
    for (int i = 0; i < num_NE; i++) {

        for (int j = 0; j < 6; j++) {
            if (s == eigen_connect(i, j) - 1) {
                record[num_array] = i; num_array++;
            }
        }
    }

    count_number[s] = num_array;

    // define the shape functions for quadratic;
    double* SH = new double[6]; double* GS = new double[2]; double** vertex = new double* [6]; double dist, min, min_x;
    for (int i = 0; i < 6; i++) {
        vertex[i] = new double[2];
    }

    // given the 6 points shape functions, construct the most near equations at the given point,

    for (int n = 0; n < num_array; n++) {

        // specify the coordinates of n^{th} elements
        for (int q = 0; q < 6; q++) {
            for (int t = 0; t < 2; t++) {
                vertex[q][t] = eigen_point(eigen_connect(record[n], q) - 1, t);
            }
        }

        min = 500.0; min_x = 500;

        for (int i = 0; i < 6; i++) {
            // initilize thr GS coordinate
            GS[0] = GS[1] = 0.0;
            // corner points
            SH[0] = (2.0 * Shape_L1_2D[i] - 1.0) * Shape_L1_2D[i]; SH[1] = (2.0 * Shape_L2_2D[i] - 1.0) * Shape_L2_2D[i]; SH[2] = (2.0 * Shape_L3_2D[i] - 1.0) * Shape_L3_2D[i];
            // mid-points
            SH[3] = 4.0 * Shape_L1_2D[i] * Shape_L2_2D[i]; SH[4] = 4.0 * Shape_L2_2D[i] * Shape_L3_2D[i]; SH[5] = 4.0 * Shape_L3_2D[i] * Shape_L1_2D[i];

            // given the specific coordinates of Gauss point
            for (int q = 0; q < 6; q++) {
                for (int j = 0; j < 2; j++) {
                    GS[j] = GS[j] + SH[q] * vertex[q][j];
                }
            }

            // examine the distance to locate the most close Gauss points
            dist = sqrt((x[0] - GS[0]) * (x[0] - GS[0]) + (x[1] - GS[1]) * (x[1] - GS[1]));
            if (dist <= min) {
                // renew the coordinate, later in the returned_coord;
                min = dist;
                returned_coord[s][n][0] = GS[0]; returned_coord[s][n][1] = GS[1];
                changed_coord[s][n][0] = GS[0]; changed_coord[s][n][1] = GS[1];
                // renew the weight of Shpae function and Weight of the total equations set up
                for (int q = 0; q < 6; q++) {
                    record_weight[s][n][q][0] = eigen_connect(record[n], q) - 1;
                    record_weight[s][n][q][1] = SH[q];
                    record_weight[s][n][q][2] = Weight_2D[i];
                }
                number_weight[s][n] = Weight_2D[i];

            }


        }

    }


    delete[] record; delete[] SH; delete[] GS;
    for (int i = 0; i < 6; i++) {
        delete vertex[i];
    }
    delete[] vertex;

}

void integratorElastic2D_arbitrary::Polyhedral_Eshelby_strain_all_6(int i, int j, int k, int l, double** coord, double* x, double* result)
{
    double bI, lp, lm, delta, para;
    double Shape_0[3] = {}; double Shape_1[3] = {}; double Shape_2[3] = {};
    double Vnorm[2], Vdir[2]; double vm[2], vp[2]; int om, op;
    para = 1.0 / (8.0 * Pi * (1.0 - nu0)); delta = 1E-7;
    double* A_ijkl = new double[6]; double* B_rest = new double[6];

    B2D_A.Triangle_shape(coord, 0, Shape_0);  B2D_A.Triangle_shape(coord, 1, Shape_1);  B2D_A.Triangle_shape(coord, 2, Shape_2);
    double** coeff = new double* [6];
    // Uniform_polygon Uni; Linear_polygon Lin; Quadratic_polygon Qua;
    for (int h = 0; h < 6; h++) {
        coeff[h] = new double[6];
    }

    B2D_A.Trangle_shape_6_coeff_6(Shape_0, Shape_1, Shape_2, coeff);

    // components to be stored
    double Psi_klij, Psi_0_klij, Psi_1_klij, Psi_01_klij, Psi_00_klij, Psi_11_klij;
    Psi_klij = Psi_0_klij = Psi_1_klij = Psi_01_klij = Psi_00_klij = Psi_11_klij = 0.0;

    double fi_ij, fi_ik, fi_il, fi_jl, fi_jk;
    fi_ij = fi_ik = fi_il = fi_jk = fi_jl = 0.0;
    double fi_0_ij, fi_0_ik, fi_0_il, fi_0_jl, fi_0_jk;
    fi_0_ij = fi_0_ik = fi_0_il = fi_0_jl = fi_0_jk = 0.0;
    double fi_1_ij, fi_1_ik, fi_1_il, fi_1_jl, fi_1_jk;
    fi_1_ij = fi_1_ik = fi_1_il = fi_1_jl = fi_1_jk = 0.0;
    double fi_01_ij, fi_01_ik, fi_01_il, fi_01_jl, fi_01_jk;
    fi_01_ij = fi_01_ik = fi_01_il = fi_01_jl = fi_01_jk = 0.0;
    double fi_00_ij, fi_00_ik, fi_00_il, fi_00_jl, fi_00_jk;
    fi_00_ij = fi_00_ik = fi_00_il = fi_00_jl = fi_00_jk = 0.0;
    double fi_11_ij, fi_11_ik, fi_11_il, fi_11_jl, fi_11_jk;
    fi_11_ij = fi_11_ik = fi_11_il = fi_11_jl = fi_11_jk = 0.0;

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

        for (int g = 0; g < 2; g++) {
            vm[g] = coord[om][g]; vp[g] = coord[op][g];
        }

        bI = lm = lp = 0.0;

        B2D_A.direction_normal(Vnorm, Vdir, vm, vp);

        for (int pp = 0; pp < 2; pp++) {
            bI += (vp[pp] - x[pp]) * Vnorm[pp];
            lm += (vm[pp] - x[pp]) * Vdir[pp];
            lp += (vp[pp] - x[pp]) * Vdir[pp];
        }                                       // compute the integral contour parameters

        if (abs(bI) < delta) bI = 0;
        if (abs(lm) < delta) lm = 0;
        if (abs(lp) < delta) lp = 0;


        // uniform potentials:
        Psi_klij += U_P.PSI_4(k, l, i, j, bI, lm, lp, Vnorm, Vdir);
        fi_ij += U_P.PHI_2(i, j, bI, lm, lp, Vnorm, Vdir);  fi_ik += U_P.PHI_2(i, k, bI, lm, lp, Vnorm, Vdir);
        fi_il += U_P.PHI_2(i, l, bI, lm, lp, Vnorm, Vdir); fi_jl += U_P.PHI_2(j, l, bI, lm, lp, Vnorm, Vdir); fi_jk += U_P.PHI_2(j, k, bI, lm, lp, Vnorm, Vdir);

        // linear potentials:
        Psi_0_klij += L_P.PSI_p_4(0, k, l, i, j, bI, lm, lp, Vnorm, Vdir, x);  Psi_1_klij += L_P.PSI_p_4(1, k, l, i, j, bI, lm, lp, Vnorm, Vdir, x);
        fi_0_ij += L_P.PHI_p_2(0, i, j, bI, lm, lp, Vnorm, Vdir, x); fi_1_ij += L_P.PHI_p_2(1, i, j, bI, lm, lp, Vnorm, Vdir, x);
        fi_0_ik += L_P.PHI_p_2(0, i, k, bI, lm, lp, Vnorm, Vdir, x); fi_1_ik += L_P.PHI_p_2(1, i, k, bI, lm, lp, Vnorm, Vdir, x);
        fi_0_il += L_P.PHI_p_2(0, i, l, bI, lm, lp, Vnorm, Vdir, x); fi_1_il += L_P.PHI_p_2(1, i, l, bI, lm, lp, Vnorm, Vdir, x);
        fi_0_jl += L_P.PHI_p_2(0, j, l, bI, lm, lp, Vnorm, Vdir, x); fi_1_jl += L_P.PHI_p_2(1, j, l, bI, lm, lp, Vnorm, Vdir, x);
        fi_0_jk += L_P.PHI_p_2(0, j, k, bI, lm, lp, Vnorm, Vdir, x); fi_1_jk += L_P.PHI_p_2(1, j, k, bI, lm, lp, Vnorm, Vdir, x);

        // Quadratic potentials:
        Psi_01_klij += Q_P.PSI_pq_4(0, 1, k, l, i, j, bI, lm, lp, Vnorm, Vdir, x); Psi_00_klij += Q_P.PSI_pq_4(0, 0, k, l, i, j, bI, lm, lp, Vnorm, Vdir, x); Psi_11_klij += Q_P.PSI_pq_4(1, 1, k, l, i, j, bI, lm, lp, Vnorm, Vdir, x);
        fi_01_ij += Q_P.PHI_pq_2(0, 1, i, j, bI, lm, lp, Vnorm, Vdir, x); fi_00_ij += Q_P.PHI_pq_2(0, 0, i, j, bI, lm, lp, Vnorm, Vdir, x); fi_11_ij += Q_P.PHI_pq_2(1, 1, i, j, bI, lm, lp, Vnorm, Vdir, x);
        fi_01_ik += Q_P.PHI_pq_2(0, 1, i, k, bI, lm, lp, Vnorm, Vdir, x); fi_00_ik += Q_P.PHI_pq_2(0, 0, i, k, bI, lm, lp, Vnorm, Vdir, x); fi_11_ik += Q_P.PHI_pq_2(1, 1, i, k, bI, lm, lp, Vnorm, Vdir, x);
        fi_01_il += Q_P.PHI_pq_2(0, 1, i, l, bI, lm, lp, Vnorm, Vdir, x); fi_00_il += Q_P.PHI_pq_2(0, 0, i, l, bI, lm, lp, Vnorm, Vdir, x); fi_11_il += Q_P.PHI_pq_2(1, 1, i, l, bI, lm, lp, Vnorm, Vdir, x);
        fi_01_jl += Q_P.PHI_pq_2(0, 1, j, l, bI, lm, lp, Vnorm, Vdir, x); fi_00_jl += Q_P.PHI_pq_2(0, 0, j, l, bI, lm, lp, Vnorm, Vdir, x); fi_11_jl += Q_P.PHI_pq_2(1, 1, j, l, bI, lm, lp, Vnorm, Vdir, x);
        fi_01_jk += Q_P.PHI_pq_2(0, 1, j, k, bI, lm, lp, Vnorm, Vdir, x); fi_00_jk += Q_P.PHI_pq_2(0, 0, j, k, bI, lm, lp, Vnorm, Vdir, x); fi_11_jk += Q_P.PHI_pq_2(1, 1, j, k, bI, lm, lp, Vnorm, Vdir, x);

    }

    double Psi_co[6] = { Psi_klij, Psi_0_klij, Psi_1_klij, Psi_01_klij, Psi_00_klij, Psi_11_klij };
    double Phi_co_ij[6] = { fi_ij, fi_0_ij, fi_1_ij, fi_01_ij, fi_00_ij, fi_11_ij };
    double Phi_co_ik[6] = { fi_ik, fi_0_ik, fi_1_ik, fi_01_ik, fi_00_ik, fi_11_ik };
    double Phi_co_il[6] = { fi_il, fi_0_il, fi_1_il, fi_01_il, fi_00_il, fi_11_il };
    double Phi_co_jk[6] = { fi_jk, fi_0_jk, fi_1_jk, fi_01_jk, fi_00_jk, fi_11_jk };
    double Phi_co_jl[6] = { fi_jl, fi_0_jl, fi_1_jl, fi_01_jl, fi_00_jl, fi_11_jl };

    for (int n = 0; n < 6; n++) {
        A_ijkl[n] = 0.0; B_rest[n] = 0.0;
        for (int h = 0; h < 6; h++) {
            A_ijkl[n] = A_ijkl[n] + coeff[n][h] * Psi_co[h];

            B_rest[n] = B_rest[n] - 2.0 * nu0 * d[k][l] * Phi_co_ij[h] * coeff[n][h]
                - (1.0 - nu0) * coeff[n][h] * (Phi_co_ik[h] * d[j][l] + Phi_co_il[h] * d[j][k] + Phi_co_jk[h] * d[i][l] + Phi_co_jl[h] * d[i][k]);
        }

        result[n] = para * (A_ijkl[n] + B_rest[n]);
    }


    delete[] A_ijkl; delete[] B_rest;
}

void integratorElastic2D_arbitrary::Eshelby_collect_normal(int i, int j, double** coord, double* x, int** index_E_ij, double*** ret_u4)
{

    int m, n, k, l, h;
    double D_ijmn[6] = {};
    for (h = 0; h < 6; h++) {
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                ret_u4[h][m][n] = 0.0;
            }
        }
    }



    for (m = 0; m < 2; m++) {

        for (n = m; n < 2; n++) {
            double result[6] = {};
            for (int sqr = 0; sqr < 6; sqr++)
                D_ijmn[sqr] = 0.0;

            for (k = 0; k < 2; k++) {
                for (l = 0; l < 2; l++) {
                    Polyhedral_Eshelby_strain_all_6(k, l, m, n, coord, x, result);
                    for (int s = 0; s < 6; s++) {
                        D_ijmn[s] = D_ijmn[s] - (cons_law(nu0, mu_0, i, j, k, l) - cons_law(nu1, mu_1, i, j, k, l)) * result[s];
                    }

                }


            }			// end of looping k
            for (h = 0; h < 6; h++)
                ret_u4[h][index_E_ij[i][j]][index_E_ij[m][n]] = ret_u4[h][index_E_ij[i][j]][index_E_ij[m][n]] + map[m][n] * D_ijmn[h];

        }


    }


}

/*
    Public members:
*/

// two empty functions
void integratorElastic2D_arbitrary::addFluxEquivalentFirstOrderInclusion(Config& config_) {};
void integratorElastic2D_arbitrary::addFluxEquivalentSecondOrderInclusion(Config& config_) {};

void integratorElastic2D_arbitrary::addFieldToBEM(Config& config_)
{
    configElastic2D_arbitrary& config = dynamic_cast<configElastic2D_arbitrary&> (config_);

    int& num_NN = config.num_NN; int& num_NE = config.num_NE;
    int& NN = config.NN;
    Ref<MatrixXd> HMAT = config.HMAT;

    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> eigen_point = config.eigen_point;
    Ref<MatrixXi> eigen_connect = config.eigen_connect;
 
#pragma omp parallel shared(HMAT)
    {
        double* x = new double[2];
        double** coord = new double* [6];
        for (int ii = 0; ii < 6; ii++) coord[ii] = new double[2];
        double*** ret_u3 = new double** [6];

        int** index;
        int id = 0;
        int pp;

        index = new int* [2];
        for (int ii = 0; ii < 2; ii++) index[ii] = new int[2];                          // 2D problem, index i,j = 1 or 2

        for (int ii = 0; ii < 2; ii++) {
            for (int jj = ii; jj < 2; jj++) {
                index[ii][jj] = 0;
            }
        }

        for (int ii = 0; ii < 6; ii++) {
            for (int jj = 0; jj < 2; jj++) {
                coord[ii][jj] = 0.0;
            }
        }

        // numbering the index for id
        for (int ii = 0; ii < 2; ii++) {
            for (int jj = ii; jj < 2; jj++) {
                index[ii][jj] = id; id++;                                              //index[0][0] = 0, ..., index[1][1] = 2
            }
        }

        for (int ii = 0; ii < 6; ii++) {
            ret_u3[ii] = new double* [2];
        }

        for (int ii = 0; ii < 6; ii++) {
            for (int jj = 0; jj < 2; jj++) {
                ret_u3[ii][jj] = new double[3];
            }
        }

        //////////////////////////////////////////////////////////////////////       finish the initialization

#pragma omp for schedule(dynamic)

        for (int MM = 0; MM < num_NE; MM++) {

            // assign the coordinate of the element
            for (int ii = 0; ii < 6; ii++) {
                for (int jj = 0; jj < 2; jj++) {
                    coord[ii][jj] = eigen_point(eigen_connect(MM, ii) - 1, jj);
                }
            }
            // obtain the coordinates of the nodes in the MM^{th} element
            for (int jj = 0; jj < NN; jj++) {
                x[0] = NODES(jj, 0);// -eigen_point(kk, 0);
                x[1] = NODES(jj, 1);// -eigen_point(kk, 1);                             // let x as the observing point, which is the nodes on the boundary


                // get_addu_tensor(x, radius[kk], DD_u3, index);                          // here, DD_u3 is Gikl, u = -Gikl * epsilon_{kl}       
                // get_addu_tensor_polygonal(x, coord, DD_u3, index, kk);

                pp = 2 * jj;                                                        // influence to which displacement field, x, y                    
                displacement_collect_normal(coord, x, index, ret_u3);


#pragma omp critical
                {
                    for (int h = 0; h < 6; h++) {
                        for (int qq = pp; qq < pp + 2; qq++) {

                            for (int mm = 0; mm < 3; mm++) {
                                HMAT(qq, 2 * int(NN) + 3 * int(eigen_connect(MM, h) - 1) + mm) += -ret_u3[h][qq - pp][mm];

                            }
                        }
                    }

                }
            }


        }// finish loop

        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 2; j++)
                delete ret_u3[i][j];
            delete ret_u3[i];
        }
        delete[] ret_u3;

        for (int ii = 0; ii < 2; ii++) {
            delete[] index[ii];
        }
        delete[] index;

        for (int i = 0; i < 6; i++) {
            delete[] coord[i];
        }
        delete[]coord;
    }


}

void integratorElastic2D_arbitrary::addFluxEquivalentInclusion(Config& config_)
{
    configElastic2D_arbitrary& config = dynamic_cast<configElastic2D_arbitrary&> (config_);

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

    double** A;                                             // create for place the eshelby tensors
    const int row = 3 * num_NN;
    A = new double* [row];

    for (int i = 0; i < row; i++) {
        A[i] = new double[row];
    }


    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            A[i][j] = 0.0; 
        }
    }

    // initilization

    Eshelby_strain_polygons(num_NN, num_NE, A, eigen_point, eigen_connect, index_E_ij, changed_coord, number_weight, count_number);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            HMAT(2 * NN + i, 2 * NN + j) = A[i][j];
        }
    }


    for (int i = 0; i < row; i++) {
        delete A[i];
    }

    delete[] A;

}

