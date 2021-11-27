
#include "configElastic2D.h"
#include "integratorElastic2D.h"
#include "integratorElastic2DHelper.h"

/*-----------------------------------------------------------------------
    functions: Eshelby's tensor for displacement, get_addu_tensor, Gikl
    u = -Gikl * epsilon_{kl}
-------------------------------------------------------------------------*/
void get_addu_tensor(double* x, double a, double** DD_u3, int** index)
{
    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            DD_u3[ii][jj] = 0.0;                                    // note i,j = 1,2, only 3 variables, ii, xy
        }
    }

    for (int i = 0; i < 2; i++) {

        for (int m = 0; m < 2; m++) {                               // m = k', n = l', i = i
            for (int n = m; n < 2; n++) {
                DD_u3[i][index[m][n]] = DD_u3[i][index[m][n]] + map[m][n] * Eshelby_displacement(x, a, i, m, n);
            }
        }
    }

}

/*--------------------------------------------------------------------------
    functions: Eshelby's tensor for linear eigen strain caused
    displacement field, get_addu1_tensor
    \int_{\omega} G_{ikl}x_p' dx'
    u_i = -M_{iklp} e_{klp}
---------------------------------------------------------------------------*/
void get_addu1_tensor(double* x, double a, double** DD_u4, int*** indexp)
{
    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 6; jj++) {
            DD_u4[ii][jj] = 0.0;
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int m = 0; m < 2; m++) {
            for (int n = m; n < 2; n++) {
                for (int p = 0; p < 2; p++) {
                     
                    DD_u4[i][indexp[m][n][p]] = DD_u4[i][indexp[m][n][p]] + map[m][n] * Eshelby_displacement_linear(x, a, i, m, n, p);
                        
                    
                }
            }
        }
    }
}

/*-----------------------------------------------------------------------
    functions: Eshelby's tensor for quadratic eigen strain caused
    displacement field, get_addu1_tensor
    \int_{\omega} G_{ikl}x_p' dx'
    u_i = -M_{iklpm} e_{klpm}
-------------------------------------------------------------------------*/
void get_addu2_tensor(double* x, double a, double** DD_u5, int**** indexpp)
{
    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 12; jj++) {
            DD_u5[ii][jj] = 0.0;
        }
    }

    for (int i = 0; i < 2; i++) {

        for (int m = 0; m < 2; m++) {
            for (int n = m; n < 2; n++) {
                for (int p = 0; p < 2; p++) {
                    for (int q = 0; q < 2; q++) {

                        DD_u5[i][indexpp[m][n][p][q]] = DD_u5[i][indexpp[m][n][p][q]] + map[m][n] * Eshelby_displacement_quadratic(x, a, i, m, n, p, q);
                    
                    }
                }
            }
        }

    }

}


void integratorElastic2D::addFieldToBEM(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> eigen_point = config.eigen_point;
    Ref<VectorXd> radius = config.radius;
    Ref<MatrixXd> NODES = config.NODES;
    int num = config.num;
    int NN = config.NN;
    
# pragma omp parallel shared (HMAT) {
        double* x = new double[2];
        double** DD_u3, ** DD_u4, ** DD_u5;

        int** index;
        int*** indexp;
        int**** indexpp;
        int id = 0;
        int pp;

        index = new int* [2];
        for (int ii = 0; ii < 2; ii++) index[ii] = new int[2];                          // 2D problem, index i,j = 1 or 2



        indexp = new int** [2];
        for (int ii = 0; ii < 2; ii++) {
            indexp[ii] = new int* [2];
            for (int jj = 0; jj < 2; jj++) {
                indexp[ii][jj] = new int[2];
            }
        }

        indexpp = new int*** [2];

        for (int ii = 0; ii < 2; ii++) {
            indexpp[ii] = new int** [2];
            for (int jj = 0; jj < 2; jj++) {
                indexpp[ii][jj] = new int* [2];
                for (int kk = 0; kk < 2; kk++) {
                    indexpp[ii][jj][kk] = new int[2];
                }
            }
        }

        for (int ii = 0; ii < 2; ii++) {
            for (int jj = 0; jj < 2; jj++) {
                index[ii][jj] = 0;
                for (int kk = 0; kk < 2; kk++) {
                    indexp[ii][jj][kk] = 0;
                    for (int ll = 0; ll < 2; ll++) {
                        indexpp[ii][jj][kk][ll] = 0;
                    }
                }
            }
        }


        // numbering the index for id
        for (int ii = 0; ii < 2; ii++) {
            for (int jj = ii; jj < 2; jj++) {
                index[ii][jj] = id; id++;                                              //index[0][0] = 0, ..., index[1][1] = 2
            }
        }

        id = 0;
        for (int ii = 0; ii < 2; ii++) {
            for (int jj = ii; jj < 2; jj++) {
                for (int kk = 0; kk < 2; kk++) {
                    indexp[ii][jj][kk] = id;
                    id++;
                }
            }
        }                                                                              // linear term has (3)*(2)

        id = 0;
        for (int ii = 0; ii < 2; ii++) {
            for (int jj = ii; jj < 2; jj++) {
                for (int kk = 0; kk < 2; kk++) {
                    for (int ll = 0; ll < 2; ll++) {
                        indexpp[ii][jj][kk][ll] = id;
                        id = id + 1;
                    }
                }
            }
        }                                                                              // quadratic term


        DD_u3 = new double* [2];
        DD_u4 = new double* [2];
        DD_u5 = new double* [2];

        for (int ii = 0; ii < 2; ii++) {
            DD_u3[ii] = new double[3];                                                // i,j = 1,2, uniform eigen strain has 3 to be determined
            DD_u4[ii] = new double[6];
            DD_u5[ii] = new double[12];
        }

        for (int ii = 0; ii < 2; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                DD_u3[ii][jj] = 0.0;
            }
            for (int kk = 0; kk < 6; kk++) {
                DD_u4[ii][kk] = 0.0;
            }
            for (int ll = 0; ll < 12; ll++) {
                DD_u5[ii][ll] = 0.0;
            }
        }

        //////////////////////////////////////////////////////////////////////       finish the initialization
#pragma omp for
    
        for (int kk = 0; kk < num; kk++) {
            for (int jj = 0; jj < NN; jj++) {
                x[0] = NODES(jj, 0) - eigen_point(kk, 0);
                x[1] = NODES(jj, 1) - eigen_point(kk, 1);

                get_addu_tensor(x, radius[kk], DD_u3, index);                          // here, DD_u3 is Gikl, u = -Gikl * epsilon_{kl}
                get_addu1_tensor(x,radius[kk],DD_u4,indexp);
                get_addu2_tensor(x, radius[kk], DD_u5, indexpp);

                pp = 2 * jj;                                                        // influence to which displacement field, x, y
                int tt = (3+6+12) * kk;                                                    // 3 * kk, 3 variables * number of inclusions

                for (int qq = pp; qq < pp + 2; qq++) {

                    for (int mm = 0; mm < 3; mm++) {
                        HMAT(qq, 2 * NN + mm + tt) = -DD_u3[qq - pp][mm];

                    }

                    for (int mm = 0; mm < 6; mm++) {
                        HMAT(qq, 2 * NN + mm + tt + 3) = -DD_u4[qq - pp][mm];
                    }

                    for (int mm = 0; mm < 12; mm++) {
                        HMAT(qq, 2 * NN + mm + tt + 9) = -DD_u5[qq - pp][mm];
                    }
                }
            }

        }                                                                            // finish loop

        for (int i = 0; i < 2; i++) {
            delete[] DD_u3[i];
            delete[] DD_u4[i];
            delete[] DD_u5[i];
        }
        delete[] DD_u3;
        delete[] DD_u4;
        delete[] DD_u5;

        for (int ii = 0; ii < 2; ii++) {
            delete[] index[ii];
            for (int jj = 0; jj < 2; jj++) {
                for (int kk = 0; kk < 2; kk++) {
                    delete[] indexpp[ii][jj][kk];
                }
                delete[] indexp[ii][jj];
            }
        }
        delete[] index;
        delete[] indexp;
        delete[] indexpp;
    }
}


void integratorElastic2D::addFluxEquivalentInclusion(Config& config_)
{
    
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> x_o = config.eigen_point;
    Ref<VectorXd> radius = config.radius;
    Ref<MatrixXd> NODES = config.NODES;
    int**& index_E_ij = config.index_E_ij;
    int***& index_E_ijk = config.index_E_ijk;
    int****& index_E_ijkl = config.index_E_ijkl;
    double** A;
    double &mu0 = config.mu_0; double &mu1 = config.mu_1;
    double &nu0 = config.nu0; double &nu1 = config.nu1;

    int num = config.num;
    int NN = config.NN;
    A = new double* [21 * num];

    for (int i = 0; i < 21 * num; i++)
        A[i] = new double[21 * num];

    for (int i = 0; i < 21 * num; i++){
        for (int j = 0; j < 21 * num; j++){
            A[i][j] = 0.0;
        }
    }
    
# pragma omp parallel shared(A,num,radius,x_o,index_E_ij,index_E_ijk,index_E_ijkl)
    {
        int t = 0; int m, n, p, r, h, k, i, j, s, l;
        double sym = 0.0;                                       // used for the diagnoal element
        int q, z;
        
#pragma omp for
        for (s = 0; s < num; s++) {
            for (i = 0; i < 2; i++) {
                for (j = i; j < 2; j++) {
                    for (h = 0; h < num; h++) {
                        for (m = 0; m < 2; m++) {
                            for (n = m; n < 2; n++) {
                                double D_ijmn = 0.0;
                                for (k = 0; k < 2; k++) {
                                    for (l = 0; l < 2; l++) {
                                        D_ijmn = D_ijmn - (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_v(k, l, m, n, s, h, radius, x_o, x_o);
                                    }
                                }
                                if (h == s) sym = 1.0;
                                else sym = 0.0;
                                A[index_E_ij[2 * s + i][j]][index_E_ij[2 * h + m][n]] = map[m][n] * (D_ijmn + cons_law(nu0, mu0, i, j, m, n) * sym);
                            
                                for (p = 0; p < 2; p++) {
                                    double D_ijmnp = 0.0;
                                    for (k = 0; k < 2; k++) {
                                        for (l = 0; l < 2; l++) {
                                            D_ijmnp = D_ijmnp - (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_linear(k,l,m,n,p,s,h,radius,x_o,x_o);

                                        }
                                    }
                                    A[index_E_ij[2 * s + i][j]][index_E_ijk[2 * h + m][n][p]] = map[m][n] * D_ijmnp;


                                }

                                for (p = 0; p < 2; p++) {
                                    for (q = 0; q < 2; q++) {
                                        double D_ijmnpq = 0.0;
                                        for (int k = 0; k < 2; k++) {
                                            for (l = 0; l < 2; l++) {

                                                D_ijmnpq = D_ijmnpq - (cons_law(nu0,mu0,i, j, k, l) - cons_law(nu1,mu1,i, j, k, l)) * Eshelby_strain_integration_quadratic(k, l, m, n, p, q, s, h, radius, x_o, x_o);

                                            }
                                        }
                                        A[index_E_ij[2 * s + i][j]][index_E_ijkl[2 * h + m][n][p][q]] = map[m][n] * D_ijmnpq;

                                    }
                                }
                            
                            
                            
                            
                            }
                        }

                    }

                }
            }
        }
        
    }
    
    for (int i = 0; i < 21 * num; i++)
    {
        for (int j = 0; j < 21 * num; j++)
        {
            HMAT(2 * NN + i, 2 * NN + j) += A[i][j];
        }
    }

    for (int i = 0; i < 21 * num; i++) {
        delete A[i];
    }
    delete[] A;

    
    
}

void integratorElastic2D::addFluxEquivalentFirstOrderInclusion(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> x_o = config.eigen_point;
    Ref<VectorXd> radius = config.radius;
    Ref<MatrixXd> NODES = config.NODES;
    int**& index_E_ij = config.index_E_ij;
    int***& index_E_ijk = config.index_E_ijk;
    int****& index_E_ijkl = config.index_E_ijkl;
    double** A;
    double &mu0 = config.mu_0; double &mu1 = config.mu_1;
    double &nu0 = config.nu0; double &nu1 = config.nu1;

    int num = config.num;
    int NN = config.NN;
    A = new double* [21 * num];

    for (int i = 0; i < 21 * num; i++)
        A[i] = new double[21 * num];

    for (int i = 0; i < 21 * num; i++){
        for (int j = 0; j < 21 * num; j++){
            A[i][j] = 0.0;
        }
    }
    
# pragma omp parallel shared(A,num,radius,x_o,index_E_ij,index_E_ijk,index_E_ijkl)
    {
        int t = 0; int m, n, p, r, h, k, i, j, s, l;
        double sym = 0.0;                                       // used for the diagnoal element
        int q, z;
        
#pragma omp for
        for (int s = 0; s < num; s++) {
            
            for (i = 0; i < 2; i++) {
                for (j = i; j < 2; j++) {
                    for (r = 0; r < 2; r++) {

                        for (h = 0; h < num; h++) {


                            for (m = 0; m < 2; m++) {
                                for (n = m; n < 2; n++) {
                                    double D_ijmn_r = 0.0;
                                    for (k = 0; k < 2; k++) {
                                        for (l = 0; l < 2; l++) {

                                            D_ijmn_r = D_ijmn_r - (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_der(k, l, m, n, r, s, h, radius, x_o, x_o);

                                        }
                                    }

                                    A[index_E_ijk[2 * s + i][j][r]][index_E_ij[2 * h + m][n]] = map[m][n] * D_ijmn_r;

                                    if (h == s)
                                        sym = 1.0;
                                    else
                                        sym = 0.0;

                                    for (p = 0; p < 2; p++) {
                                        double D_ijmnp_r = 0.0;
                                        for (k = 0; k < 2; k++) {
                                            for (l = 0; l < 2; l++) {

                                                D_ijmnp_r = D_ijmnp_r - (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_linear_der(k, l, m, n, p, r, s, h, radius, x_o, x_o);

                                            }
                                        }
                                        A[index_E_ijk[2 * s + i][j][r]][index_E_ijk[2 * h + m][n][p]] = map[m][n] * (D_ijmnp_r + sym * cons_law(nu0, mu0, i, j, m, n) * d[p][r]);
                                    }


                                    for (p = 0; p < 2; p++) {
                                        for (q = 0; q < 2; q++) {
                                            double D_ijmnpq_r = 0.0;
                                            for (int k = 0; k < 2; k++) {
                                                for (l = 0; l < 2; l++) {

                                                    D_ijmnpq_r = D_ijmnpq_r - (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_quadratic_der(k, l, m, n, p, q, r, s, h, radius, x_o, x_o);

                                                }
                                            }
                                            A[index_E_ijk[2 * s + i][j][r]][index_E_ijkl[2 * h + m][n][p][q]] = map[m][n] * D_ijmnpq_r;
                                        }
                                    }


                                }
                            }
                        }

                        //    t = t + 1;
                    }
                }
            }
            
        }
    }
    
    
    for (int i = 0; i < 21 * num; i++)
    {
        for (int j = 0; j < 21 * num; j++)
        {
            HMAT(2 * NN + i, 2 * NN + j) += A[i][j];
        }
    }

    for (int i = 0; i < 21 * num; i++) {
        delete A[i];
    }
    delete[] A;
}


void integratorElastic3D::addFluxEquivalentSecondOrderInclusion(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> x_o = config.eigen_point;
    Ref<VectorXd> radius = config.radius;
    Ref<MatrixXd> NODES = config.NODES;
    int**& index_E_ij = config.index_E_ij;
    int***& index_E_ijk = config.index_E_ijk;
    int****& index_E_ijkl = config.index_E_ijkl;
    double** A;
    double &mu0 = config.mu_0; double &mu1 = config.mu_1;
    double &nu0 = config.nu0; double &nu1 = config.nu1;

    int num = config.num;
    int NN = config.NN;
    A = new double* [21 * num];

    for (int i = 0; i < 21 * num; i++)
        A[i] = new double[21 * num];

    for (int i = 0; i < 21 * num; i++){
        for (int j = 0; j < 21 * num; j++){
            A[i][j] = 0.0;
        }
    }
    
# pragma omp parallel shared(A,num,radius,x_o,index_E_ij,index_E_ijk,index_E_ijkl)
    {
        int t = 0; int m, n, p, r, h, k, i, j, s, l;
        double sym = 0.0;                                       // used for the diagnoal element
        int q, z;
        
#pragma omp for
        for (int s = 0; s < num; s++) {
            for (int i = 0; i < 2; i++) {
                for (int j = i; j < 2; j++) {
                    for (r = 0; r < 2; r++) {
                        for (z = 0; z < 2; z++) {
                            for (h = 0; h < num; h++) {
                                for (m = 0; m < 2; m++) {
                                    for (n = m; n < 2; n++) {
                                        double D_ijmn_rz = 0.0;
                                        for (int k = 0; k < 2; k++) {
                                            for (l = 0; l < 2; l++) {

                                                D_ijmn_rz = D_ijmn_rz - 0.5 * (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_uniform_der2(k, l, m, n, r, z, s, h, radius, x_o, x_o);

                                            }
                                        }
                                        A[index_E_ijkl[2 * s + i][j][r][z]][index_E_ij[2 * h + m][n]] = map[m][n] * D_ijmn_rz;

                                        for (p = 0; p < 2; p++) {
                                            double D_ijmnp_rz = 0.0;
                                            for (int k = 0; k < 2; k++) {
                                                for (l = 0; l < 2; l++) {

                                                    D_ijmnp_rz = D_ijmnp_rz - 0.5 * (cons_law(nu0, mu0, i, j, k, l) - cons_law(nu1, mu1, i, j, k, l)) * Eshelby_strain_integration_linear_der2(k, l, m, n, p, r, z, s, h, radius, x_o, x_o);

                                                }
                                            }
                                            A[index_E_ijkl[2 * s + i][j][r][z]][index_E_ijk[2 * h + m][n][p]] = map[m][n] * D_ijmnp_rz;
                                        }

                                        if (h == s)
                                            sym = 1.0;
                                        else
                                            sym = 0.0;

                                        for (p = 0; p < 2; p++) {
                                            for (q = 0; q < 2; q++) {
                                                double D_ijmnpq_rz = 0.0;
                                                for (int k = 0; k < 2; k++) {
                                                    for (l = 0; l < 2; l++) {

                                                        D_ijmnpq_rz = D_ijmnpq_rz - 0.5 * (cons_law(nu0,mu0,i, j, k, l) - cons_law(nu1,mu1,i, j, k, l)) * Eshelby_strain_integration_quadratic_der2(k, l, m, n, p, q, r, z, s, h, radius, x_o, x_o);

                                                    }
                                                }

                                                A[index_E_ijkl[2 * s + i][j][r][z]][index_E_ijkl[2 * h + m][n][p][q]] = map[m][n] * (D_ijmnpq_rz + sym * cons_law(nu0,mu0,i, j, m, n) * d[p][r] * d[q][z]);

                                            }
                                        }

                                    }
                                }
                            }

                            //    t = t + 1;
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < 21 * num; i++)
    {
        for (int j = 0; j < 21 * num; j++)
        {
            HMAT(2 * NN + i, 2 * NN + j) += A[i][j];
        }
    }

    for (int i = 0; i < 21 * num; i++) {
        delete A[i];
    }
    delete[] A;
}
