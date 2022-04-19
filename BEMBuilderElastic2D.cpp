
#include "configElastic2D.h"
#include "integratorElastic2DHelper.h"
#include "BEMBuilderElastic2D.h"
#include "BEMElasticHelper2D.h"

/*
 add the essential functions combined with Green's function and their derivatives
 */
/*------------------------------------------------------------------------------------------
    Obtain the strain field and multiply with the constitutive law of inclusion
-------------------------------------------------------------------------------------------*/
void BEMBuilderElastic2D::sr_ijk_inclusion(double *Vnorm, double *point, double *GCcor, double ** Dst, double **Dsu)
{
    int t = 0;

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 2; j++){
            Dst[i][j] = Dsu[i][j] = 0.0;
        }
    }

    double x[2] = { GCcor[0] - point[0], GCcor[1] - point[1] };

    double TG1[2];
    double TG[4];
    double sgm[4]; // 4 = 2*2, i,j = 1,2

    for (int i = 0; i < 2; i++){ TG1[i] = 0.0; }
    for (int i = 0; i < 4; i++){ TG[i] = 0.0; }
    for (int i = 0; i < 4; i++){ sgm[i] = 0.0; }

    for (int n = 0; n < 2; n++){
        for (int i = 0; i < 2; i++){
            TG1[n] = TG1[n] + B2D.Green_01(x, n, i, i);
        }
    }

    int tt = 0;
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int l = 0; l < 2; l++){
                TG[tt] = TG[tt] + B2D.Green_02(x, l, i, l, j);
            }
            tt++;
        }
    }

    tt = 0;
    double lambda0 = (2.0 * mu_0 * nu0) / (1.0 - 2.0 * nu0);
    double lambda1 = (2.0 * mu_1 * nu1) / (1.0 - 2.0 * nu1);
    // to obtain the strain field by the original solution
    for (int m = 0; m < 2; m++){
        for (int n = 0; n < 2; n++){
            for (int i = 0; i < 2; i++){
                sgm[tt] = sgm[tt] + lambda0 * d[m][n] * TG[2 * i + i] + 2.0*mu_0* B2D.Green_02(x, m, i, n, i);
            }
            tt++;
        }
    }


    double ***DSU, ***DST;
    DSU = new double **[2];
    DST = new double **[2];

    for (int i = 0; i < 2; i++){
        DSU[i] = new double*[2];
        DST[i] = new double *[2];
        for (int j = 0; j < 2; j++){
            DSU[i][j] = new double[2];
            DST[i][j] = new double[2];
        }
    }

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                DSU[i][j][k] = DST[i][j][k] = 0.0;
            }
        }
    }

    for (int k = 0; k < 2; k++){
        for (int l = 0; l < 2; l++){
            for (int n = 0; n < 2; n++){
                for (int m = 0; m < 2; m++){
                    DST[k][l][n] += ((lambda0 - lambda1)*d[k][l] * sgm[2 * m + n] + (mu_0 - mu_1)*(2.0*mu_0*nu0 / (1 - 2 * nu0)*d[m][n] * (TG[2 * k + l] + TG[2 * l + k]) + mu_0* B2D.Green_02(x, m, k, n, l) + mu_0* B2D.Green_02(x, m, l, n, k) + mu_0* B2D.Green_02(x, n, k, m, l) + mu_0* B2D.Green_02(x, n, l, m, k)))* Vnorm[m];

                }

                DSU[k][l][n] += -((lambda0 - lambda1)*d[k][l] * TG1[n] + (mu_0 - mu_1)* B2D.Green_01(x, n, k, l) + (mu_0 - mu_1)* B2D.Green_01(x, n, l, k));
            }
        }
    }

    t = 0;
    for (int i = 0; i < 2; i++){
        for (int j = i; j < 2; j++){
            for (int m = 0; m < 2; m++){
                Dsu[t][m] = DSU[i][j][m];
                Dst[t][m] = DST[i][j][m];
            }
            t++;
        }
    }

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++)
        {
            delete[] DSU[i][j];
            delete[] DST[i][j];

        }
        delete[] DSU[i];
        delete[] DST[i];

    }

    delete[] DSU;
    delete[] DST;

}

/*----------------------------------------------------------------------------------------
        Obtain the strain field and multiply with the constitutive law of inclusion
------------------------------------------------------------------------------------------*/
void BEMBuilderElastic2D::sr_ijk_inclusion_linear(double* Vnorm, double* point, double* GCcor, double** Dstw, double** Dsuw)
{
    int t = 0;

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 2; j++) {
            Dstw[i][j] = Dsuw[i][j] = 0.0;
        }
    }

    double x[2];
    x[0] = GCcor[0] - point[0];
    x[1] = GCcor[1] - point[1];

    double TG1[4];
    double TG[8];
    double sgm[8];

    for (int i = 0; i < 2 * 2; i++) { TG1[i] = 0.0; }
    for (int i = 0; i < 2 * 2 * 2; i++) { TG[i] = 0.0; }
    for (int i = 0; i < 2 * 2 * 2; i++) { sgm[i] = 0.0; }

    int tt = 0;

    for (int n = 0; n < 2; n++) {
        for (int w = 0; w < 2; w++) {
            for (int i = 0; i < 2; i++) {
                TG1[tt] = TG1[tt] + B2D.Green_02(x, n, i, i, w);
            }
            tt = tt + 1;
        }
    }

    tt = 0;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int w = 0; w < 2; w++) {

                for (int l = 0; l < 2; l++) {
                    TG[tt] = TG[tt] + B2D.Green_03(x, l, i, l, j, w);
                }

                tt = tt + 1;
            }
        }
    }

    tt = 0;
    for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) {
            for (int w = 0; w < 2; w++) {

                for (int i = 0; i < 2; i++) {
                    sgm[tt] = sgm[tt] + 2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * TG[4 * i + 2 * i + w] + mu_0 * B2D.Green_03(x, m, i, n, i, w) + mu_0 * B2D.Green_03(x, n, i, m, i, w);
                }

                tt = tt + 1;
            }
        }
    }

    for (int k = 0; k < 2; k++) {
        for (int l = k; l < 2; l++) {
            for (int w = 0; w < 2; w++) {

                for (int n = 0; n < 2; n++) {

                    for (int m = 0; m < 2; m++) {

                        Dsuw[t][n] = Dsuw[t][n] - ((2.0 * mu_0 * nu0 / (1 - 2 * nu0) - 2.0 * mu_1 * nu1 / (1 - 2 * nu1)) * d[k][l] * sgm[4 * m + 2 * n + w] + (mu_0 - mu_1) * (2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * (TG[4 * k + 2 * l + w] + TG[4 * l + 2 * k + w]) + mu_0 * B2D.Green_03(x, m, k, n, l, w) + mu_0 * B2D.Green_03(x, m, l, n, k, w) + mu_0 * B2D.Green_03(x, n, k, m, l, w) + mu_0 * B2D.Green_03(x, n, l, m, k, w))) * Vnorm[m];
                    }


                    Dstw[t][n] = Dstw[t][n] + ((2.0 * mu_0 * nu0 / (1 - 2 * nu0) - 2.0 * mu_1 * nu1 / (1 - 2 * nu1)) * d[k][l] * TG1[2 * n + w] + (mu_0 - mu_1) * B2D.Green_02(x, n, k, l, w) + (mu_0 - mu_1) * B2D.Green_02(x, n, l, k, w));

                }

                t = t + 1;
            }
        }
    }


}

/*---------------------------------------------------------------------------------------
        Obtain the strain field and multiply with the constitutive law of inclusion
----------------------------------------------------------------------------------------*/
void BEMBuilderElastic2D::sr_ijk_inclusion_quadratic(double* Vnorm, double* point, double* GCcor, double** Dstwz, double** Dsuwz)
{
    int t = 0;

    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 2; j++) {
            Dstwz[i][j] = Dsuwz[i][j] = 0.0;
        }
    }

    double x[2];
    x[0] = GCcor[0] - point[0];
    x[1] = GCcor[1] - point[1];

    double TG1[2 * 2 * 2];
    double TG[2 * 2 * 2 * 2];
    double sgm[2 * 2 * 2 * 2];

    for (int i = 0; i < 2 * 2 * 2; i++) { TG1[i] = 0.0; }
    for (int i = 0; i < 2 * 2 * 2 * 2; i++) { TG[i] = 0.0; }
    for (int i = 0; i < 2 * 2 * 2 * 2; i++) { sgm[i] = 0.0; }

    int tt = 0;

    for (int n = 0; n < 2; n++) {
        for (int w = 0; w < 2; w++) {
            for (int z = 0; z < 2; z++) {

                for (int i = 0; i < 2; i++) {
                    TG1[tt] = TG1[tt] + B2D.Green_03(x, n, i, i, w, z);
                }
                tt++;
            }
        }
    }

    tt = 0;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int w = 0; w < 2; w++) {
                for (int z = 0; z < 2; z++) {

                    for (int l = 0; l < 2; l++) {
                        TG[tt] = TG[tt] + B2D.Green_04(x, l, i, l, j, w, z);
                    }

                    tt = tt + 1;
                }
            }
        }
    }

    tt = 0;
    for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 2; n++) {
            for (int w = 0; w < 2; w++) {
                for (int z = 0; z < 2; z++) {

                    for (int i = 0; i < 2; i++) {
                        sgm[tt] = sgm[tt] + 2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * TG[8 * i + 4 * i + 2 * w + z]
                                + mu_0 * B2D.Green_04(x, m, i, n, i, w, z)
                                + mu_0 * B2D.Green_04(x, n, i, m, i, w, z);
                    }

                    tt = tt + 1;
                }
            }
        }
    }

    for (int k = 0; k < 2; k++) {
        for (int l = k; l < 2; l++) {
            for (int w = 0; w < 2; w++) {
                for (int z = 0; z < 2; z++) {

                    for (int n = 0; n < 2; n++) {

                        for (int m = 0; m < 2; m++) {

                            Dsuwz[t][n] = Dsuwz[t][n] + ((2.0 * mu_0 * nu0 / (1 - 2 * nu0) - 2.0 * mu_1 * nu1 / (1 - 2 * nu1)) * d[k][l] * sgm[8 * m + 4 * n + 2 * w + z] + (mu_0 - mu_1) * (2.0 * mu_0 * nu0 / (1 - 2 * nu0) * d[m][n] * (TG[8 * k + 4 * l + 2 * w + z] + TG[8 * l + 4 * k + 2 * w + z]) + mu_0 * B2D.Green_04(x, m, k, n, l, w, z) + mu_0 * B2D.Green_04(x, m, l, n, k, w, z) + mu_0 * B2D.Green_04(x, n, k, m, l, w, z) + mu_0 * B2D.Green_04(x, n, l, m, k, w, z))) * Vnorm[m];
                        }


                        Dstwz[t][n] = Dstwz[t][n] - ((2.0 * mu_0 * nu0 / (1 - 2 * nu0) - 2.0 * mu_1 * nu1 / (1 - 2 * nu1)) * d[k][l] * TG1[4 * n + 2 * w + z] + (mu_0 - mu_1) * B2D.Green_03(x, n, k, l, w, z) + (mu_0 - mu_1) * B2D.Green_03(x, n, l, k, w, z));

                    }

                    t = t + 1;
                }
            }
        }
    }


}




/*---------------------------------------------------------------------------------------
        Integ2e_stress created for BEM in equilibrium equation
-----------------------------------------------------------------------------------------*/
void BEMBuilderElastic2D::Integ2e_stress(double** Elcor, double** DDsu, double** DDst, double* point)
{
    int i, m, k, ngp;
    double xsi, W, Jac, N1, N2, N3;

    double** Dsu, ** Dst, ** Ni, ** TDDsu, ** TDDst;

    int Ndof = 2;

    double* N = new double[3];
    double* Vnorm = new double[2];
    double* GCcor = new double[2];


    Dsu = new double* [3];
    Dst = new double* [3];
    Ni = new double* [2];

    for (i = 0; i < 3; i++) { Dsu[i] = new double[2]; Dst[i] = new double[2]; }

    for (i = 0; i < 2; i++) { Ni[i] = new double[6]; }


    ngp = 8;
    double w[8] = { 0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459, 0.3137066459, 0.3626837834, 0.3626837834 };
    double gp[8] = { 0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425 };

    for (m = 0; m < ngp; m++) {
        xsi = gp[m]; W = w[m];

        B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);
        B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);
        B2D.Cartesian(Elcor, N1, N2, N3, GCcor);

        //get_for_stress_BEM(Vnorm, point, GCcor, Dst, Dsu);
        sr_ijk_inclusion(Vnorm, point, GCcor, Dst, Dsu);


        //for (int mm = 0; mm < 3; mm++) {
          //  for (int ll = 0; ll < 6; ll++) {
            //    for (int nn = 0; nn < 2; nn++) {
              ///      DDsu[mm][ll] = DDsu[mm][ll] + Dsu[mm][nn] * Ni[nn][ll] * W * Jac;
                 //   DDst[mm][ll] = DDst[mm][ll] + Dst[mm][nn] * Ni[nn][ll] * W * Jac;
                //}

 //           }
        
        //}


        for (int mm = 0; mm < 3; mm++) {
            DDsu[mm][0] = DDsu[mm][0] + Dsu[mm][0] * N[0] * W * Jac;
            DDsu[mm][1] = DDsu[mm][1] + Dsu[mm][1] * N[0] * W * Jac;
            DDsu[mm][2] = DDsu[mm][2] + Dsu[mm][0] * N[1] * W * Jac;
            DDsu[mm][3] = DDsu[mm][3] + Dsu[mm][1] * N[1] * W * Jac;
            DDsu[mm][4] = DDsu[mm][4] + Dsu[mm][0] * N[2] * W * Jac;
            DDsu[mm][5] = DDsu[mm][5] + Dsu[mm][1] * N[2] * W * Jac;


            DDst[mm][0] = DDst[mm][0] + Dst[mm][0] * N[0] * W * Jac;
            DDst[mm][1] = DDst[mm][1] + Dst[mm][1] * N[0] * W * Jac;
            DDst[mm][2] = DDst[mm][2] + Dst[mm][0] * N[1] * W * Jac;
            DDst[mm][3] = DDst[mm][3] + Dst[mm][1] * N[1] * W * Jac;
            DDst[mm][4] = DDst[mm][4] + Dst[mm][0] * N[2] * W * Jac;
            DDst[mm][5] = DDst[mm][5] + Dst[mm][1] * N[2] * W * Jac;

        }
    }       // end loop with ngp

    delete[] N;
    delete[] Vnorm;
    delete[] GCcor;

    for (int i = 0; i < 3; i++) {
        delete Dsu[i];
        delete Dst[i];
    }

    for (int i = 0; i < 2; i++) {
        delete Ni[i];
    }
    delete[] Dsu;
    delete[] Dst;
    delete[] Ni;
}

/*---------------------------------------------------------------------------------------
        Compute stress for linear strain terms for EIM
-----------------------------------------------------------------------------------------*/
void BEMBuilderElastic2D::Integ2e_stress_linear(double** Elcor, double** DDsuw, double** DDstw, double* point)
{
    int i, m, k, ngp;

    double xsi, eta, W, Jac, N1, N2, N3, N4;
    double** Dsuw, ** Dstw, ** Ni;

    int Ndof = 2;

    double* N = new double[3];
    double* Vnorm = new double[2];
    double* GCcor = new double[2];


    Dsuw = new double* [6];
    Dstw = new double* [6];
    Ni = new double* [2];

    for (i = 0; i < 6; i++) { Dsuw[i] = new double[2]; Dstw[i] = new double[2]; }

    for (i = 0; i < 2; i++) { Ni[i] = new double[6]; }
    
    ngp = 8;
    double w[8] = { 0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459, 0.3137066459, 0.3626837834, 0.3626837834 };
    double gp[8] = { 0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425 };

    for (m = 0; m < ngp; m++) {
        xsi = gp[m]; W = w[m];
        B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);
        B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);
        B2D.Cartesian(Elcor, N1, N2, N3, GCcor);

        sr_ijk_inclusion_linear(Vnorm, point, GCcor, Dstw, Dsuw);

        for (int mm = 0; mm < 6; mm++) {
            for (int ll = 0; ll < 6; ll++) {
                for (int nn = 0; nn < 2; nn++) {
                    DDsuw[mm][ll] = DDsuw[mm][ll] + Dsuw[mm][nn] * Ni[nn][ll] * W * Jac;
                    DDstw[mm][ll] = DDstw[mm][ll] + Dstw[mm][nn] * Ni[nn][ll] * W * Jac;
                }
            }
        }

    }

    delete[] N;
    delete[] Vnorm;
    delete[] GCcor;

    for (int i = 0; i < 6; i++) {
        delete[] Dstw[i];
        delete[] Dsuw[i];
    }
    for (int i = 0; i < 2; i++) {
        delete[] Ni[i];
    }

    delete[] Ni;
    delete[] Dsuw;
    delete[] Dstw;


}

/*---------------------------------------------------------------------------------------
        Compute stress for quadratic strain terms for EIM
----------------------------------------------------------------------------------------*/
void BEMBuilderElastic2D::Integ2e_stress_quadratic(double** Elcor, double** DDsuwz, double** DDstwz, double* point)
{
    int i, m, k, ngp;

    double xsi, eta, W, Jac, N1, N2, N3, N4;
    double** Dsuwz, ** Dstwz, ** Ni;

    int Ndof = 2;

    double* N = new double[3];
    double* Vnorm = new double[2];
    double* GCcor = new double[2];


    Dsuwz = new double* [12];
    Dstwz = new double* [12];
    Ni = new double* [2];

    for (i = 0; i < 12; i++) { Dsuwz[i] = new double[2]; Dstwz[i] = new double[2]; }

    for (i = 0; i < 2; i++) { Ni[i] = new double[6]; }

    ngp = 8;
    double w[8] = { 0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459, 0.3137066459, 0.3626837834, 0.3626837834 };
    double gp[8] = { 0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425 };

    for (m = 0; m < ngp; m++) {
        xsi = gp[m]; W = w[m];
        B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);
        B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);
        B2D.Cartesian(Elcor, N1, N2, N3, GCcor);

        sr_ijk_inclusion_quadratic(Vnorm, point, GCcor, Dstwz, Dsuwz);

        for (int mm = 0; mm < 12; mm++) {
            for (int ll = 0; ll < 6; ll++) {
                for (int nn = 0; nn < 2; nn++) {
                    DDsuwz[mm][ll] = DDsuwz[mm][ll] + Dsuwz[mm][nn] * Ni[nn][ll] * W * Jac;
                    DDstwz[mm][ll] = DDstwz[mm][ll] + Dstwz[mm][nn] * Ni[nn][ll] * W * Jac;
                }
            }
        }

    }

    delete[] N;
    delete[] Vnorm;
    delete[] GCcor;

    for (int i = 0; i < 12; i++) {
        delete[] Dstwz[i];
        delete[] Dsuwz[i];
    }
    for (int i = 0; i < 2; i++) {
        delete[] Ni[i];
    }

    delete[] Ni;
    delete[] Dsuwz;
    delete[] Dstwz;
}


/*--------------------------------------------------------------------
        Integ3 for the 2D elastic boundary integrals
----------------------------------------------------------------------*/
void BEMBuilderElastic2D::Integ3(double E, double nu, int* sctrb, double** Elcor, double** dUe, double** dTe, Ref<MatrixXd> xP, int Ncol)
{
    /*----------------------------------------
     ngp: number of Gauss integration point
     judge: whether conincide with same node
     W: weight funtion of Gauss
    ------------------------------------------*/
    int i, n, m, k, ii, jj, judge, l_pi, iD, nD, ntr, Ntri, ngp;
    double xsi, W, r, Jac, xsib, N1, N2, N3, Jacb, rjudge;

    double** UP, ** TP, ** Ni;
    double epsilon = 0.0001;           // to judge whether is too close and we could neglect the TK
    int Ndof = 2;

    double* dxr = new double[2];       // 2 dof, hence 2 derivatives
    double* N = new double[3];         // 3 nodes
    double* Vnorm = new double[2];     // x,y 2D coordinate
    double* GCcor = new double[2];

    UP = new double* [2]; TP = new double* [2]; Ni = new double* [2];

    for (i = 0; i < 2; i++) {
        UP[i] = new double[2]; TP[i] = new double[2]; Ni[i] = new double[6];
    }
    double cons = (1.0 + nu) / (4.0 * pi * E * (1.0 - nu));
    double consp = cons * (3.0 - 4.0 * nu);

    ngp = 8;
    double w[8] = { 0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459, 0.3137066459, 0.3626837834, 0.3626837834 };
    double gp[8] = { 0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425 };

    int ngl = 4;                                   // the Gauss-laguerre integral created for log integral on the conincide points
    double gl[4] = { 0.0414484801, 0.245274914, 0.556165453, 0.848982394 };
    double wl[4] = { 0.383464068, 0.386875317, 0.190435126, 0.0392254871 };

    Ntri = 2;                                      // division of the integration on the element, if judge == 1, Ntri = 2, otherwise Ntri = 1

    judge = 0;                                     // judge whether we have the coincide nodes
    for (i = 0; i < Ncol; i++) {
        //judge = 0;
        //for (n = 0; n < 3; n++) {
        //    if ((sctrb[n] - 1) == i) {               // because the real nodes start from 1 not 0
            //    judge = 1; l_pi = n;                 // record which node in the element coincide
            //}
        //}


        //if (judge == 0) {
        for (m = 0; m < ngp; m++) {
            xsi = gp[m];                            // only 1 GP coordinate
            W = w[m];
            B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);  // function for shape functions
            B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);     // Jacobian and normal unit vector
            B2D.Cartesian(Elcor, N1, N2, N3, GCcor);     // change coordinates to the intrinsic way
            r = sqrt((GCcor[0] - xP(i, 0)) * (GCcor[0] - xP(i, 0)) + (GCcor[1] - xP(i, 1)) * (GCcor[1] - xP(i, 1)));

            dxr[0] = (GCcor[0] - xP(i, 0)) / r;
            dxr[1] = (GCcor[1] - xP(i, 1)) / r;

            B2D.UK(dxr, r, E, nu, UP);                  // For H matrix
            B2D.TK(dxr, r, E, nu, TP, Vnorm);           // For G matrix

            for (ii = 0; ii < Ndof; ii++) {
                iD = Ndof * i + ii;                 // insert position of the integration in the matrix

                for (jj = 0; jj < Ndof; jj++) {
                    for (n = 0; n < 3; n++) {
                        nD = Ndof * n + jj;
                        double r_check = sqrt((Elcor[n][0] - xP(i, 0)) * (Elcor[n][0] - xP(i, 0)) + (Elcor[n][1] - xP(i, 1)) * (Elcor[n][1] - xP(i, 1)));
                        if (r_check >= 0.0001) {
                            dUe[iD][nD] = dUe[iD][nD] + N[n] * UP[ii][jj] * Jac * W;
                            dTe[iD][nD] = dTe[iD][nD] + N[n] * TP[ii][jj] * Jac * W;
                        }
                        else {
                            dUe[iD][nD] = dUe[iD][nD] + N[n] * cons * dxr[jj] * dxr[ii] * Jac * W;
                        }
                    }


                }
            }


        }
    } // end of the node is not the same

       // }                                                 // end of the judge == 0 case

        //else {                                            // judge == 1, the coincide happens
    for (i = 0; i < Ncol; i++) {
        for (n = 0; n < 3; n++) {
            // check the need by r_check
            double r_check = sqrt((Elcor[n][0] - xP(i, 0)) * (Elcor[n][0] - xP(i, 0)) + (Elcor[n][1] - xP(i, 1)) * (Elcor[n][1] - xP(i, 1)));
            if (r_check >= 0.0001) continue;                // means no coincidence


            double Elength = 0.0; double xxi, wwi; double dxdxb; double conspp = 0.0;
            for (int kkk = 0; kkk < ngp; kkk++) {         // compute the element length
                xxi = gp[kkk]; wwi = w[kkk];
                Elength += B2D.Compute_length(xxi, wwi, Elcor);
            }

            if (n != 2) Ntri = 1;                        // if node == 3, 2 regions
            else Ntri = 2;
            for (ntr = 0; ntr < Ntri; ntr++) {
                for (m = 0; m < ngl; m++) {
                    if (n == 0) {
                        xsi = 2.0 * gl[m] - 1.0;            // intrinsic parameters
                        dxdxb = 2.0;
                    }
                    else if (n == 1) {
                        xsi = 1.0 - 2.0 * gl[m];
                        dxdxb = 2.0;
                    }
                    else {
                        dxdxb = 1.0;
                        if (ntr == 0) xsi = -gl[m];
                        else xsi = gl[m];
                    }

                    B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);
                    B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);

                    for (ii = 0; ii < Ndof; ii++) {
                        iD = Ndof * i + ii;
                        // for (jj = 0; jj < Ndof; jj++) {

                        // for (n = 0; n < 3; n++) {
                        nD = Ndof * n + ii;

                        dUe[iD][nD] = dUe[iD][nD] + N[n] * consp * Jac * dxdxb * wl[m];
                        // }
                        // }
                    }

                }

            }                                                                // complete the log integration
            // normal integration
            for (m = 0; m < ngp; m++) {

                if (n == 0) {
                    conspp = -log(Elength) * consp;
                }
                else if (n == 1) {
                    conspp = -log(Elength) * consp;
                }
                else if (n == 2) {
                    conspp = (log(2.0) - log(Elength)) * consp;
                }                                                                // the constant part based on the which n


                xsi = gp[m]; W = w[m];

                B2D.Serendip_func(xsi, N1, N2, N3, Ni, N);
                B2D.Normal_Jac(xsi, Elcor, Jac, Vnorm);

                for (ii = 0; ii < Ndof; ii++) {
                    iD = Ndof * i + ii;
                    //for (jj = 0; jj < Ndof; jj++) {

                    // for (n = 0; n < 3; n++) {
                    nD = Ndof * n + ii;

                    dUe[iD][nD] = dUe[iD][nD] + N[n] * conspp * Jac * W;
                    // }
                    //}
                }

            }
        }
        //}
    }


    // delete all the memories to avoid errors
    delete[] N;
    delete[] Vnorm;
    delete[] GCcor;
    for (int i = 0; i < 2; i++) {
        delete UP[i];
        delete TP[i];
        delete Ni[i];
    }
    delete[] UP;
    delete[] TP;
    delete[] Ni;
}

// construct conventional BEM part
void BEMBuilderElastic2D::addBEM(Config& config_)
{
    
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);
    
    double& EE = config.EE;
    double& XNU = config.XNU;
    int& NE = config.NE;
    int& NN = config.NN;
    Ref<MatrixXi> NConnect = config.NConnect;
    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> GMAT = config.GMAT;
    int Dimension = config.Dimension;
    
    // the diag matrix to cancel strong singularities
    MatrixXd Diag = MatrixXd::Zero(Dimension * NN, Dimension);
    
#pragma omp parallel shared(HMAT,GMAT)
    {
        double** dUe; double** dTe;
        dUe = new double* [Dimension * NN]; dTe = new double* [Dimension * NN];
        for (int j = 0; j < Dimension * NN; j++) {
            dUe[j] = new double[6]; dTe[j] = new double[6];              // 6 = 2 (dof) * 3 (node number for one element)
        }

        double** Elcor;
        Elcor = new double* [3];
        for (int j = 0; j < 3; j++) {
            Elcor[j] = new double[Dimension];                                    // 3 nodes, x,y = 2-D problem, record the coordinate information of each node
        }

        int* sctrb = new int[3];                                         // 3 nodes, record the global number of the nodes

        ////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for
        for (int i = 0; i < NE; i++) {
            sctrb[0] = NConnect(i, 0);
            sctrb[1] = NConnect(i, 1);
            sctrb[2] = NConnect(i, 2);                                    // record the global number of one element nodes

            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < Dimension; k++) {
                    Elcor[j][k] = NODES(sctrb[j] - 1, k);                 // record the coodinates of each node
                }

            }

            for (int ii = 0; ii < Dimension * NN; ii++) {
                for (int jj = 0; jj < 6; jj++) {
                    dUe[ii][jj] = 0.0; dTe[ii][jj] = 0.0;
                }
            }

            Integ3(EE, XNU, sctrb, Elcor, dUe, dTe, NODES, NN);

            for (int l = 0; l < Dimension * NN; l++) {
                for (int k = 0; k < 6; k++) {
                    GMAT(l, 6 * i + k) = dUe[l][k];
                }
            }

            // end for loop on NE

////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp critical
            {
                for (int jj = 0; jj < 3; jj++) {
                    int pp = 0;
                    pp = Dimension * (sctrb[jj] - 1);
                    for (int kk = 0; kk < Dimension; kk++) {
                        for (int ll = 0; ll < Dimension * NN; ll++) {

                            HMAT(ll, pp + kk) = HMAT(ll, pp + kk) + dTe[ll][Dimension * jj + kk];
                        }
                    }

                }

            }

        }               // end loop for NE


        for (int i = 0; i < Dimension * NN; i++) {
            delete dTe[i];
            delete dUe[i];
        }

        delete[] dTe;
        delete[] dUe;
    }                                                                   // end of pragma 1

    ////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < NN; i++) {
        for (int j = 0; j < Dimension; j++) {
            int n = i * Dimension + j;

            for (int k = 0; k < Dimension * NN; k++) {
                Diag(k, j) = Diag(k, j) + HMAT(k, n);
            }

        }
    }

    for (int m = 0; m < NN; m++) {
        int posi = Dimension * m;
        for (int n = posi; n < posi + Dimension; n++) {
            for (int l = posi; l < posi + Dimension; l++) {
                HMAT(n, l) = -Diag(n, l - posi);
            }
        }
    }
    Diag.resize(0, 0);
}


// the normal stress, 0th order stress equivalent conditions
void BEMBuilderElastic2D::addFluxEquivalentBEM(Config& config_)
{
    
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);

    double& EE = config.EE;
    double& XNU = config.XNU;
    int& NE = config.NE;
    int& NN = config.NN;
    int& num = config.num;
    Ref<MatrixXi> NConnect = config.NConnect;
    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> GMAT = config.GMAT;
    Ref<MatrixXd> eigen_point = config.eigen_point;
    int Dimension = config.Dimension;
    
#pragma omp parallel shared(HMAT,GMAT)
    {
        double** DDsu, ** DDst;                                             // one from the Hmat, and one from Gmat
        double* point = new double[Dimension];
        DDsu = new double* [3];
        DDst = new double* [3];                                             // 3 because the variable to be determined is 3 for 2D

        double** Elcor;
        Elcor = new double* [3];                                            // 3 node element

        for (int j = 0; j < 3; j++) Elcor[j] = new double[Dimension];               // 2D problem

        for (int i = 0; i < 3; i++) {
            DDsu[i] = new double[6]; DDst[i] = new double[6];
        }

        int* sctrb = new int[3];                                           // 3 node element

        ////////////////////////////////////////////////////////////////

#pragma omp for

        for (int ee = 0; ee < NE; ee++) {

            for (int kk = 0; kk < num; kk++) {

                for (int mm = 0; mm < 3; mm++) {
                    for (int ll = 0; ll < 6; ll++) {
                        DDsu[mm][ll] = 0.0; DDst[mm][ll] = 0.0;
                    }
                }

                point[0] = eigen_point(kk, 0); point[1] = eigen_point(kk, 1);
                sctrb[0] = NConnect(ee, 0); sctrb[1] = NConnect(ee, 1); sctrb[2] = NConnect(ee, 2);

                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < Dimension; k++) {
                        Elcor[j][k] = NODES(sctrb[j] - 1, k);
                    }
                }

                Integ2e_stress(Elcor, DDsu, DDst, point);
                //////////////////////////////////////////////////////////////////////////
#pragma omp critical
                {
                    for (int j = 0; j < 3; j++) {

                        int p = 0;
                        p = Dimension * (sctrb[j] - 1);

                        for (int k = 0; k < Dimension; k++) {
                            for (int l = 0; l < 3; l++) {

                                HMAT(Dimension * NN + 21 * kk + l, p + k) = HMAT(Dimension * NN + 21 * kk + l, p + k) - DDst[l][Dimension * j + k];
                            }
                        }
                    }

                }   // end critical

                for (int l = 0; l < 3; l++) {
                    for (int k = 0; k < 6; k++) {
                        GMAT(Dimension * NN + 21 * kk + l, 6 * ee + k) = DDsu[l][k];
                    }
                }

            }                                                // end loop kk



        }

        for (int i = 0; i < 3; i++) {
            delete DDsu[i];
            delete DDst[i];
        }

        for (int i = 0; i < 3; i++) {
            delete Elcor[i];
        }

        delete[] DDsu;
        delete[] DDst;
        delete[] sctrb;
        delete[] point;
        delete[] Elcor;

    }
    
}

//the normal stress, 1st order stress equivalent conditions
void BEMBuilderElastic2D::addFluxEquivalentFirstOrderBEM(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);

    double& EE = config.EE;
    double& XNU = config.XNU;
    int& NE = config.NE;
    int& NN = config.NN;
    int& num = config.num;
    Ref<MatrixXi> NConnect = config.NConnect;
    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> GMAT = config.GMAT;
    Ref<MatrixXd> eigen_point = config.eigen_point;
    int Dimension = config.Dimension;
    
#pragma omp parallel shared (HMAT,GMAT)
    {

        double** DDsuw, ** DDstw;
        double* point = new double[Dimension];

        DDsuw = new double* [6];
        DDstw = new double* [6];

        double** Elcor = new double* [3];
        for (int j = 0; j < 3; j++) {
            Elcor[j] = new double[Dimension];
        }


        for (int i = 0; i < 6; i++) {
            DDsuw[i] = new double[6]; DDstw[i] = new double[6];
        }

        int* sctrb = new int[3];

#pragma omp for
        for (int ee = 0; ee < NE; ee++) {

            for (int KK = 0; KK < num; KK++) {

                for (int mm = 0; mm < 6; mm++) {
                    for (int ll = 0; ll < 6; ll++) {
                        DDsuw[mm][ll] = DDstw[mm][ll] = 0.0;
                    }
                }           // initilization of the DDsutw tensors

                point[0] = eigen_point(KK, 0); point[1] = eigen_point(KK, 1);
                sctrb[0] = NConnect(ee, 0); sctrb[1] = NConnect(ee, 1); sctrb[2] = NConnect(ee, 2);

                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < Dimension; k++) {
                        Elcor[j][k] = NODES(sctrb[j] - 1, k);
                    }
                }
                

                Integ2e_stress_linear(Elcor, DDsuw, DDstw, point);


#pragma omp critical
                {

                    for (int j = 0; j < 3; j++)
                    {
                        int p = 0;
                        p = Dimension * (sctrb[j] - 1);

                        for (int k = 0; k < Dimension; k++) {

                            for (int l = 0; l < 6; l++) {

                                HMAT(Dimension * NN + 21 * KK + l + 3, p + k) = HMAT(Dimension * NN + 21 * KK + l + 3, p + k) - DDsuw[l][2 * j + k];
                            }
                        }
                    }


                }   // end of critical

                for (int l = 0; l < 6; l++) {
                    for (int k = 0; k < 6; k++) {
                        GMAT(Dimension * NN + 21 * KK + l + 3, 6 * ee + k) = DDstw[l][k];
                    }
                }

            }

        }       // end of looping elememt

        for (int i = 0; i < 6; i++) {
            delete DDsuw[i];
            delete DDstw[i];
        }

        for (int i = 0; i < 3; i++) {
            delete Elcor[i];
        }
        delete[] DDsuw;
        delete[] DDstw;
        delete[] sctrb;
        delete[] Elcor;
        delete[] point;

    }                       // end of pragma first

    
}

//the normal stress, 2nd order stress equivalent conditions
void BEMBuilderElastic2D::addFluxEquivalentSecondOrderBEM(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);

    double& EE = config.EE;
    double& XNU = config.XNU;
    int& NE = config.NE;
    int& NN = config.NN;
    int& num = config.num;
    Ref<MatrixXi> NConnect = config.NConnect;
    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> HMAT = config.HMAT;
    Ref<MatrixXd> GMAT = config.GMAT;
    Ref<MatrixXd> eigen_point = config.eigen_point;
    int Dimension = config.Dimension;
    
#pragma omp parallel shared (HMAT,GMAT)
    {

        double** DDsuwz, ** DDstwz;
        double* point = new double[Dimension];

        DDsuwz = new double* [12];
        DDstwz = new double* [12];

        double** Elcor = new double* [3];
        for (int j = 0; j < 3; j++) {
            Elcor[j] = new double[Dimension];
        }


        for (int i = 0; i < 12; i++) {
            DDsuwz[i] = new double[6]; DDstwz[i] = new double[6];
        }

        int* sctrb = new int[3];

#pragma omp for
        for (int ee = 0; ee < NE; ee++) {

            for (int KK = 0; KK < num; KK++) {

                for (int mm = 0; mm < 12; mm++) {
                    for (int ll = 0; ll < 6; ll++) {
                        DDsuwz[mm][ll] = DDstwz[mm][ll] = 0.0;
                    }
                }           // initilization of the DDsutw tensors

                point[0] = eigen_point(KK, 0); point[1] = eigen_point(KK, 1);
                sctrb[0] = NConnect(ee, 0); sctrb[1] = NConnect(ee, 1); sctrb[2] = NConnect(ee, 2);

                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < Dimension; k++) {
                        Elcor[j][k] = NODES(sctrb[j] - 1, k);
                    }
                }


                Integ2e_stress_quadratic(Elcor, DDsuwz, DDstwz, point);


#pragma omp critical
                {

                    for (int j = 0; j < 3; j++)
                    {
                        int p = 0;
                        p = Dimension * (sctrb[j] - 1);

                        for (int k = 0; k < Dimension; k++) {

                            for (int l = 0; l < 12; l++) {

                                HMAT(Dimension * NN + 21 * KK + l + 3 + 6, p + k) = HMAT(Dimension * NN + 21 * KK + l + 3 + 6, p + k) - 0.5 * DDsuwz[l][2 * j + k];
                            }
                        }
                    }


                }   // end of critical

                for (int l = 0; l < 12; l++) {
                    for (int k = 0; k < 6; k++) {
                        GMAT(Dimension * NN + 21 * KK + l + 3 + 6, 6 * ee + k) = 0.5 * DDstwz[l][k];
                    }
                }

            }

        }       // end of looping elememt

        for (int i = 0; i < 12; i++) {
            delete DDsuwz[i];
            delete DDstwz[i];
        }

        for (int i = 0; i < 3; i++) {
            delete Elcor[i];
        }
        delete[] DDsuwz;
        delete[] DDstwz;
        delete[] sctrb;
        delete[] Elcor;
        delete[] point;

    }                       // end of pragma first


}


// rearrange the matrix based on BCs

// judge the connections in the element
int  BEMBuilderElastic2D::ffloord(int i, int j)
{

    int d;

    if (i == 0 || i == 1) {
        d = 0;
    }

    if (i == 2 || i == 3) {
        d = 1;
    }

    if (i == 4 || i == 5) {
        d = 2;
    }


    return d;

}

/*----------------------------------------------------------------
        Construction for the Find function
        Find the node number with element connections
        Therefore, the parameters for the displacement field could
        be placed together to solve the BVP
------------------------------------------------------------------*/
int BEMBuilderElastic2D::Find(int h, int e, int i, double** record, Ref<MatrixXi> NConnect, int NE)
{
    int t, p;
    p = 2;
    t = ffloord(i, p);

    for (int j = 0; j < 6 * NE; j++) {
        if ((NConnect(e, t) + 0.1 * (i % 2 + 1)) == record[j][1]) {
            h = j;
            goto pao;
        }
    }

pao:
    ;
    return h;
}

void BEMBuilderElastic2D::ApplyBC(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);

    double& EE = config.EE;
    double& XNU = config.XNU;
    int& NE = config.NE;
    int& NN = config.NN;
    int& num = config.num;
    Ref<MatrixXi> NConnect = config.NConnect;
    Ref<MatrixXd> NODES = config.NODES;
    Ref<MatrixXd> BC   = config.BC;
    Ref<MatrixXd> t_bc = config.t_bc;
    MatrixXd& HMAT = config.HMAT;
    MatrixXd& GMAT = config.GMAT;
    MatrixXd& AA = config.A;
    VectorXd& BB = config.b;

    ////////////////////////////////////////////////////////////////
    
    record = new double*[NE * 6];               // recording 2 traction directions * 3 node per element
    for (int i = 0; i < 6 * NE; i++) {
        record[i] = new double[2];
    }

    for (int i = 0; i < 6 * NE; i++) {
        for (int j = 0; j < 2; j++) {
            record[i][j] = 0.0;
        }
    }
    
    /////////////////////////////////////////////////////////////////////
    int t, h; t = 0;

    for (int j = 0; j < NN * 2; j++) {
        if (BC(j, 0) == 1) {

            for (int i = 0; i < NN * 2 + 21 * num; i++) {
                AA(i, t) = HMAT(i, j);
            }
            t++;
        }
    }                                                                               // BC defines the boundary types for displacement
    // BC(xx,0) = 1, it is unknown, t records how many unknown for the nodal displacement
    for (int j = 0; j < NN * 2; j++) {
        if (BC(j, 0) == 0) {
            for (int i = 0; i < NN * 2 + 21 * num; i++) {
                BB(i) = BB(i) - HMAT(i, j) * BC(j, 1);                              // knowns to right
            }

        }
    }

    int mm = 0; int p = 2;
    //p = 3?

    for (int e = 0; e < NE; e++) {

        for (int i = 0; i < 6; i++) {

            if (t_bc(e, i) == 1) {
                for (int row = 0; row < NN * 2 + 21 * num; row++) {
                    BB(row) = BB(row) + t_bc(e, 6 + i) * GMAT(row, 6 * e + i);
                }
            }                                                                           // traction BC, 1 is known

            else if (t_bc(e, i) == 0) {
                h = -1;
                h = Find(h, e, i, record, NConnect, NE);

                if (h == -1) {
                    record[mm][0] = t; int tt = ffloord(i, p);
                    record[mm][1] = NConnect(e, tt) + 0.1 * (i % 2 + 1);
                    mm++;

                    for (int ii = 0; ii < NN * 2 + 21 * num; ii++) {
                        AA(ii, t) = -GMAT(ii, 6 * e + i);
                    }
                    t++;
                }

                else {

                    int cc = record[h][0];
                    for (int ii = 0; ii < NN * 2 + 21 * num; ii++) {
                        AA(ii, cc) = AA(ii, cc) - GMAT(ii, 6 * e + i);
                    }

                }
            }
        }


    }

    for (int i = 0; i < 2 * NN + 21 * num; i++) {
        for (int j = 2 * NN; j < 2 * NN + 21 * num; j++) {
            AA(i, j) = HMAT(i, j);
        }
    }



    config.HMAT.resize(0, 0);
    config.GMAT.resize(0, 0);
}

// obtain the UT vectors
void BEMBuilderElastic2D::getUTVector(Config& config_)
{
    configElastic2D& config = dynamic_cast<configElastic2D&> (config_);

    VectorXd& U = config.U;
    VectorXd& T = config.T;
    MatrixXd& BC = config.BC;
    MatrixXd& t_bc = config.t_bc;
    VectorXd& XX = config.X;
    MatrixXi& NConnect = config.NConnect;

    int& NN = config.NN;
    int& num = config.num;
    int& NE = config.NE;

    
    int t,h;
    /*--------------------------------
    Insert back the data
    -----------------------------------*/
    t = h = 0; int tt = 0;

    for (int i = 0; i < NN * 2; i++) {
        if (BC(i, 0) == 0) {
            U(i) = BC(i, 1);
        }
        else {
            U(i) = XX(t);
            t++;                                                    // t is the record number of solution
        }
    }

    for (int i = 2 * NN; i < NN * 2 + 21 * num; i++) {
        U(i) = XX(i);
    }

    int m = 0;

    for (int e = 0; e < NE; e++) {
        for (int i = 0; i < 6; i++) {

            if (t_bc(e, i) == 1) {
                T(m) = t_bc(e, i + 6);
                m++;
            }
            else {
                h = Find(h, e, i, record, NConnect, NE);
                tt = record[h][0];
                T(m) = XX(tt);
                m++;
            }

        }

    }
    
}
