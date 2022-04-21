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


# include "BEMElasticHelper2D_arbitrary.h"
# include "configElastic2D_arbitrary.h"
#include "GlobalConstant.h"

/*----------------------------------------------------------
        Apply shape function in local and global coordinate
------------------------------------------------------------*/
void BEM_2D_arbitrary::Cartesian(double** Elcor, double& N1, double& N2, double& N3, double* GCcor)
{
    for (int i = 0; i < 2; i++) {
        GCcor[i] = N1 * Elcor[0][i] + N2 * Elcor[1][i] + N3 * Elcor[2][i];
    }
}

/*-----------------------------------------------------------
        Serendip_func: calculate the shape function of
        the quadratic 2D element
------------------------------------------------------------*/
void BEM_2D_arbitrary::Serendip_func(double xsi, double& N1, double& N2, double& N3, double** Ni, double* N)
{
    // Ni is used for int_stress function
    // j is the product of node freedom 2 times nodes number 3
    // the shape functions will be used for internal area
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
            Ni[i][j] = 0.0;
        }
    }

    N1 = 0.5 * (1.0 - xsi); N2 = 0.5 * (1.0 + xsi);
    // since it is 3 node problem
    N3 = 1.0 - (xsi * xsi);
    N1 = N1 - 0.5 * N3; N2 = N2 - 0.5 * N3;
    N[0] = N1; N[1] = N2; N[2] = N3;

    // the Ni functions are created for easy use in strain field
    Ni[0][0] = Ni[1][1] = N1;
    Ni[0][2] = Ni[1][3] = N2;
    Ni[0][4] = Ni[1][5] = N3;

}

/*---------------------------------------------------------------
           Normal_jacobian function
           Calculate the derivative of the shape functions
           and jacobian, the length of the element
-----------------------------------------------------------------*/
void BEM_2D_arbitrary::Normal_Jac(double xsi, double** Elcor, double& Jac, double* Vnorm)
{
    double dN1_dxsi, dN2_dxsi, dN3_dxsi;
    double V_xsi[2]; // the vnorm

    dN3_dxsi = -2.0 * xsi;
    dN1_dxsi = -0.5 - 0.5 * dN3_dxsi;
    dN2_dxsi = 0.5 - 0.5 * dN3_dxsi;

    // 3 shape functions, x,y for positions
    for (int i = 0; i < 2; i++) {
        V_xsi[i] = dN1_dxsi * Elcor[0][i] + dN2_dxsi * Elcor[1][i] + dN3_dxsi * Elcor[2][i];
    }
    /*--------------------------------
    math: V[0] = dy/dxi, V[1] = -dx/dxi
    ----------------------------------*/
    // the normal vector (not unit)
    Vnorm[0] = V_xsi[1];
    Vnorm[1] = -V_xsi[0];
    // the third direction does not matter  Vnorm[2] = 0.0
    // Jacobian calculation
    Jac = sqrt(V_xsi[0] * V_xsi[0] + V_xsi[1] * V_xsi[1]);
    // unit outward normal vector
    Vnorm[0] = Vnorm[0] / Jac; Vnorm[1] = Vnorm[1] / Jac;

}

/*------------------------------------------------------------------
        Compute the element length using Gauss integral method
--------------------------------------------------------------------*/
double BEM_2D_arbitrary::Compute_length(double xsi, double weight, double** Elcor)
{
    // loop n times in the external loops
    double V_xsi[2]; double N1, N2, N3, dElength;
    double dN3_dxsi, dN1_dxsi, dN2_dxsi, Jac;
    N1 = 0.5 * (1.0 - xsi); N2 = 0.5 * (1.0 + xsi);
    // since it is 3 node problem
    N3 = 1.0 - (xsi * xsi);
    N1 = N1 - 0.5 * N3; N2 = N2 - 0.5 * N3;
    dN3_dxsi = -2.0 * xsi;
    dN1_dxsi = -0.5 - 0.5 * dN3_dxsi;
    dN2_dxsi = 0.5 - 0.5 * dN3_dxsi;
    for (int i = 0; i < 2; i++) {
        V_xsi[i] = dN1_dxsi * Elcor[0][i] + dN2_dxsi * Elcor[1][i] + dN3_dxsi * Elcor[2][i];
    }
    Jac = sqrt(V_xsi[0] * V_xsi[0] + V_xsi[1] * V_xsi[1]);
    dElength = Jac * weight;

    return dElength;
}


/*------------------------------------------------------------------
          The boundary integral origin function
          Tk (H),  Uk(G), plane strain problem
---------------------------------------------------------------------*/
void BEM_2D_arbitrary::UK(double* dxr, double r, double E, double nu, double** UP)
{
    double mu, c1, conr, c, clog;
    // lame constants
    mu = E / (2.0 * (1.0 + nu));
    c1 = 3.0 - 4.0 * nu;
    c = 1.0 / (8.0 * Pi * mu * (1.0 - nu));
    clog = -c1 * log(r);

    // the green's functions
    UP[0][0] = c * (clog + dxr[0] * dxr[0]);
    UP[0][1] = c * dxr[0] * dxr[1];
    UP[1][0] = UP[0][1];
    UP[1][1] = c * (clog + dxr[1] * dxr[1]);
}

/*---------------------------------------------------------------------
          The boundary integral origin function
          Tk (H),  Uk(G), plane strain problem
-----------------------------------------------------------------------*/
void BEM_2D_arbitrary::TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm)
{
    double c3, conr, Costh, c2;
    c3 = 1.0 - 2.0 * nu;
    c2 = 1.0 / (4.0 * Pi * (1.0 - nu));
    conr = c2 / r;
    Costh = Vnorm[0] * dxr[0] + Vnorm[1] * dxr[1];

    // dot product of the TK
    TP[0][0] = -(conr * (c3 + 2.0 * dxr[0] * dxr[0]) * Costh);
    TP[1][1] = -(conr * (c3 + 2.0 * dxr[1] * dxr[1]) * Costh);
    int i, j;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            if (i != j) {
                TP[i][j] = -(conr * (2.0 * dxr[i] * dxr[j] * Costh - c3 * (Vnorm[j] * dxr[i] - Vnorm[i] * dxr[j])));
            }
        }
    }

}

/*
    Define the partial derivatives of Green's function
*/
/*----------------------------------------------------------------------------------------------*/
double BEM_2D_arbitrary::Green_01(double* x, int i, int j, int k)
{
    // first derivative of Green's function
    double result, para;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    para = (1.0) / (8.0 * Pi * (1.0 - nu0) * mu_0);
    double cons = -(3.0 - 4.0 * nu0);
    double term1 = cons * ((d[i][j] * x[k]) / (r * r));
    double term2 = (d[i][k] * x[j] + d[j][k] * x[i]) / (r * r);
    double term3 = (-2.0 / (r * r * r * r)) * (x[i] * x[j] * x[k]);

    result = para * (term1 + term2 + term3);
    return result;

}

/*-----------------------------------------------------------------------------------------------*/
double BEM_2D_arbitrary::Green_02(double* x, int i, int j, int k, int l)
{
    // second derivative of the Green's function
    double result, para, cons;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    para = (1.0) / (8.0 * Pi * (1.0 - nu0) * mu_0);
    cons = -(3.0 - 4.0 * nu0);

    double term1 = cons * d[i][j] * ((d[k][l] / (r * r)) - (2.0 / pow(r, 4)) * x[k] * x[l]);
    double term2 = (1.0 / (r * r)) * (d[i][k] * d[j][l] + d[j][k] * d[i][l]);
    double term3 = (-2.0 / pow(r, 4)) * x[l] * (d[i][k] * x[j] + d[j][k] * x[i]);
    double term4 = (-2.0 / pow(r, 4)) * (d[i][l] * x[j] * x[k] + d[j][l] * x[i] * x[k] + d[k][l] * x[i] * x[j]);
    double term5 = (8.0 / pow(r, 6)) * x[i] * x[k] * x[l] * x[j];

    result = para * (term1 + term2 + term3 + term4 + term5);
    return result;
}

/*-----------------------------------------------------------------------------------------------*/
double BEM_2D_arbitrary::Green_03(double* x, int i, int j, int k, int l, int s)
{
    // third derivatives of the Green's function
    double result, para, cons;
    double term1, term2, term3, term4, term5, term6, term7, term8;
    result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    para = (1.0) / (8.0 * Pi * (1.0 - nu0) * mu_0);
    cons = -(3.0 - 4.0 * nu0);

    term1 = cons * d[i][j] * (-d[k][l] * 2.0 * x[s] * pow(r, -4.0) - 2.0 * (pow(r, -4.0) * (d[s][k] * x[l] + d[l][s] * x[k]) - 4.0 * pow(r, -6.0) * x[s] * x[k] * x[l]));
    term2 = -2.0 * pow(r, -4.0) * x[s] * (d[i][k] * d[j][l] + d[i][l] * d[j][k]);
    term8 = 8.0 * pow(r, -6.0) * (d[i][k] * x[j] * x[l] + d[j][k] * x[i] * x[l]) * x[s];
    term3 = -2.0 * pow(r, -4.0) * (d[i][k] * (d[l][s] * x[j] + d[j][s] * x[l]) + d[j][k] * (d[i][s] * x[l] + d[l][s] * x[i]));
    term4 = -2.0 * pow(r, -4.0) * (d[i][l] * (d[j][s] * x[k] + d[k][s] * x[j]) + d[j][l] * (d[s][i] * x[k] + d[s][k] * x[i]) + d[k][l] * (d[i][s] * x[j] + d[j][s] * x[i]));
    term5 = 8.0 * pow(r, -6.0) * x[s] * (d[i][l] * x[j] * x[k] + d[j][l] * x[i] * x[k] + d[k][l] * x[i] * x[j]);
    term6 = 8.0 * pow(r, -6.0) * (d[i][s] * x[j] * x[k] * x[l] + d[j][s] * x[i] * x[k] * x[l] + d[k][s] * x[i] * x[j] * x[l] + d[l][s] * x[i] * x[j] * x[k]);
    term7 = -48.0 * pow(r, -8.0) * x[i] * x[j] * x[k] * x[l] * x[s];

    result = para * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8);

    return result;
}

/*------------------------------------------------------------------------------------------------*/
double BEM_2D_arbitrary::Green_04(double* x, int i, int j, int k, int l, int s, int m)
{
    // fourth derivatives of the Green's function
    double result;
    double para = (1.0) / (8.0 * Pi * (1.0 - nu0) * mu_0);
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double term1 = -(3.0 - 4.0 * nu0) * d[i][j] * (
        -2.0 * d[k][l] * (pow(r, -4.0) * d[s][m] - 4.0 * x[s] * x[m] * pow(r, -6.0))
        - 2.0 * ((d[k][s] * d[l][m] + d[l][s] * d[k][m]) * pow(r, -4.0) - 4.0 * x[m] * pow(r, -6.0) * (d[k][s] * x[l] + d[l][s] * x[k]))
        + 8.0 * ((d[s][m] * x[k] * x[l] + d[k][m] * x[s] * x[l] + d[l][m] * x[s] * x[k]) * pow(r, -6.0) - 6.0 * (x[s] * x[k] * x[l] * x[m]) * pow(r, -8.0))
        );

    double term2 = -2.0 * (d[i][k] * d[j][l] + d[i][l] * d[j][k]) * (d[s][m] * pow(r, -4.0) - 4.0 * x[s] * x[m] * pow(r, -6.0));

    double term3 = 8.0 * pow(r, -6.0) * (d[s][m] * (d[i][k] * x[j] * x[l] + d[j][k] * x[l] * x[i]) + x[s] * (d[i][k] * (d[j][m] * x[l] + d[l][m] * x[j]) + d[j][k] * (d[i][m] * x[l] + d[l][m] * x[i])))
        - 48.0 * x[s] * x[m] * pow(r, -8.0) * (d[i][k] * x[j] * x[l] + d[j][k] * x[i] * x[l]);

    double term4 = -2.0 * ((d[i][k] * (d[l][s] * d[j][m] + d[j][s] * d[l][m]) + d[j][k] * (d[i][s] * d[l][m] + d[l][s] * d[i][m])) * pow(r, -4.0)
        - 4.0 * x[m] * pow(r, -6.0) * (d[i][k] * (d[l][s] * x[j] + d[j][s] * x[l]) + d[j][k] * (d[i][s] * x[l] + d[l][s] * x[i])));

    double term5 = -2.0 * ((d[i][l] * (d[j][s] * d[k][m] + d[k][s] * d[j][m]) + d[j][l] * (d[s][i] * d[k][m] + d[s][k] * d[i][m]) + d[k][l] * (d[i][s] * d[j][m] + d[j][s] * d[i][m])) * pow(r, -4.0)
        - 4.0 * (d[i][l] * (d[j][s] * x[k] + d[k][s] * x[j]) + d[j][l] * (d[s][i] * x[k] + d[s][k] * x[i]) + d[k][l] * (d[i][s] * x[j] + d[j][s] * x[i])) * x[m] * pow(r, -6.0));

    double term6 = 8.0 * ((d[s][m] * (d[i][l] * x[j] * x[k] + d[j][l] * x[i] * x[k] + d[k][l] * x[i] * x[j]) + x[s] * (d[i][l] * (d[j][m] * x[k] + d[k][m] * x[j]) + d[j][l] * (d[i][m] * x[k] + d[k][m] * x[i]) + d[k][l] * (d[i][m] * x[j] + d[j][m] * x[i]))) * pow(r, -6.0)
        - 6.0 * pow(r, -8.0) * x[s] * x[m] * (d[i][l] * x[j] * x[k] + d[j][l] * x[i] * x[k] + d[k][l] * x[i] * x[j]));

    double term7 = 8.0 * (
        (d[i][s] * (d[j][m] * x[k] * x[l] + d[k][m] * x[j] * x[l] + d[l][m] * x[k] * x[j])
            + d[j][s] * (d[i][m] * x[k] * x[l] + d[k][m] * x[i] * x[l] + d[l][m] * x[i] * x[k])
            + d[k][s] * (d[i][m] * x[j] * x[l] + d[j][m] * x[i] * x[l] + d[l][m] * x[i] * x[j])
            + d[l][s] * (d[i][m] * x[j] * x[k] + d[j][m] * x[i] * x[k] + d[k][m] * x[i] * x[j])) * pow(r, -6.0)
        - 6.0 * pow(r, -8.0) * x[m] * (d[i][s] * x[j] * x[k] * x[l] + d[j][s] * x[i] * x[k] * x[l] + d[k][s] * x[i] * x[j] * x[l] + d[l][s] * x[i] * x[j] * x[k])
        );

    double term8 = -48.0 * ((d[i][m] * x[j] * x[k] * x[l] * x[s]
        + d[j][m] * x[i] * x[k] * x[l] * x[s]
        + d[k][m] * x[i] * x[j] * x[l] * x[s]
        + d[l][m] * x[i] * x[j] * x[k] * x[s]
        + d[s][m] * x[i] * x[j] * x[k] * x[l]) * pow(r, -8.0)
        - 8.0 * pow(r, -10.0) * x[i] * x[j] * x[k] * x[l] * x[s] * x[m]);

    result = para * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8);
    return result;
}

double BEM_2D_arbitrary::Triangle_area(double** coor)
{
    double result = 0.5 * (coor[1][0] * coor[2][1] + coor[2][0] * coor[0][1] + coor[0][0] * coor[1][1]
        - coor[0][1] * coor[1][0] - coor[1][1] * coor[2][0] - coor[2][1] * coor[0][0]);
    return abs(result);
}

void BEM_2D_arbitrary::Triangle_shape(double** coor, int n, double* SH)
{
    // n is used to show which one is the I node
    int j, k;
    if (n == 0) {
        j = 1; k = 2;
    }
    else if (n == 1) {
        j = 2; k = 0;
    }
    else if (n == 2) {
        j = 0; k = 1;
    }
    // note that n,j,k is in counter clockwise sequence
    double area = Triangle_area(coor);
    double alpha, beta, gamma;
    alpha = coor[j][0] * coor[k][1] - coor[k][0] * coor[j][1];
    beta = coor[j][1] - coor[k][1];
    gamma = coor[k][0] - coor[j][0];

    SH[0] = alpha / (2.0 * area);
    SH[1] = beta / (2.0 * area);
    SH[2] = gamma / (2.0 * area);

}

double BEM_2D_arbitrary::Triangle_shape_6(double** coord, int n, double* x)
{
    double SH[3] = {}; double L;

    switch (n) {
    case 0:
        Triangle_shape(coord, n, SH);
        L = SH[0] + SH[1] * x[0] + SH[2] * x[1];
        L = (2.0 * L - 1.0) * L;
        break;
    case 1:
        Triangle_shape(coord, n, SH);
        L = SH[0] + SH[1] * x[0] + SH[2] * x[1];
        L = (2.0 * L - 1.0) * L;
        break;
    case 2:
        Triangle_shape(coord, n, SH);
        L = SH[0] + SH[1] * x[0] + SH[2] * x[1];
        L = (2.0 * L - 1.0) * L;
        break;
    case 3:
        Triangle_shape(coord, 0, SH);
        L = SH[0] + SH[1] * x[0] + SH[2] * x[1];
        Triangle_shape(coord, 1, SH);
        L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1]);
        break;
    case 4:
        Triangle_shape(coord, 1, SH);
        L = SH[0] + SH[1] * x[0] + SH[2] * x[1];
        Triangle_shape(coord, 2, SH);
        L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1]);
        break;
    case 5:
        Triangle_shape(coord, 2, SH);
        L = SH[0] + SH[1] * x[0] + SH[2] * x[1];
        Triangle_shape(coord, 0, SH);
        L = 4.0 * L * (SH[0] + SH[1] * x[0] + SH[2] * x[1]);
        break;

    }

    return L;
}

void BEM_2D_arbitrary::Trangle_shape_6_coeff_6(double* SH_0, double* SH_1, double* SH_2, double** coeff)
{
    for (int i = 0; i < 6; i++) {
        for (int h = 0; h < 6; h++) {
            coeff[i][h] = 0.0;
        }
    }

    // define the first 3 corner nodes, vertex 0
    coeff[0][0] = 2.0 * SH_0[0] * SH_0[0] - SH_0[0]; coeff[0][1] = 4.0 * (SH_0[0] * SH_0[1]) - SH_0[1]; coeff[0][2] = 4.0 * (SH_0[0] * SH_0[2]) - SH_0[2];
    coeff[0][3] = 4.0 * SH_0[1] * SH_0[2]; coeff[0][4] = 2.0 * SH_0[1] * SH_0[1]; coeff[0][5] = 2.0 * SH_0[2] * SH_0[2];

    // vertex 1,
    coeff[1][0] = 2.0 * SH_1[0] * SH_1[0] - SH_1[0]; coeff[1][1] = 4.0 * (SH_1[0] * SH_1[1]) - SH_1[1]; coeff[1][2] = 4.0 * (SH_1[0] * SH_1[2]) - SH_1[2];
    coeff[1][3] = 4.0 * SH_1[1] * SH_1[2]; coeff[1][4] = 2.0 * SH_1[1] * SH_1[1]; coeff[1][5] = 2.0 * SH_1[2] * SH_1[2];

    // vertex 2,
    coeff[2][0] = 2.0 * SH_2[0] * SH_2[0] - SH_2[0]; coeff[2][1] = 4.0 * (SH_2[0] * SH_2[1]) - SH_2[1]; coeff[2][2] = 4.0 * (SH_2[0] * SH_2[2]) - SH_2[2];
    coeff[2][3] = 4.0 * SH_2[1] * SH_2[2]; coeff[2][4] = 2.0 * SH_2[1] * SH_2[1]; coeff[2][5] = 2.0 * SH_2[2] * SH_2[2];

    // mid-points, vertex 3, (0 and 1)
    coeff[3][0] = 4.0 * SH_0[0] * SH_1[0]; coeff[3][1] = 4.0 * (SH_0[0] * SH_1[1] + SH_0[1] * SH_1[0]); coeff[3][2] = 4.0 * (SH_0[0] * SH_1[2] + SH_0[2] * SH_1[0]);
    coeff[3][3] = 4.0 * (SH_0[1] * SH_1[2] + SH_0[2] * SH_1[1]); coeff[3][4] = 4.0 * SH_0[1] * SH_1[1]; coeff[3][5] = 4.0 * SH_0[2] * SH_1[2];

    // mid-points, vertex 4, (1 and 2)
    coeff[4][0] = 4.0 * SH_2[0] * SH_1[0]; coeff[4][1] = 4.0 * (SH_2[0] * SH_1[1] + SH_2[1] * SH_1[0]); coeff[4][2] = 4.0 * (SH_2[0] * SH_1[2] + SH_2[2] * SH_1[0]);
    coeff[4][3] = 4.0 * (SH_2[1] * SH_1[2] + SH_2[2] * SH_1[1]); coeff[4][4] = 4.0 * SH_2[1] * SH_1[1]; coeff[4][5] = 4.0 * SH_2[2] * SH_1[2];

    // mid-points, vertex 5, (2 and 0)
    coeff[5][0] = 4.0 * SH_2[0] * SH_0[0]; coeff[5][1] = 4.0 * (SH_2[0] * SH_0[1] + SH_2[1] * SH_0[0]); coeff[5][2] = 4.0 * (SH_2[0] * SH_0[2] + SH_2[2] * SH_0[0]);
    coeff[5][3] = 4.0 * (SH_2[1] * SH_0[2] + SH_2[2] * SH_0[1]); coeff[5][4] = 4.0 * SH_2[1] * SH_0[1]; coeff[5][5] = 4.0 * SH_2[2] * SH_0[2];


}


void BEM_2D_arbitrary::direction_normal(double* Vnorm, double* Vdir, double* x1, double* x2)
{
    double len = sqrt((x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[0] - x2[0]) * (x1[0] - x2[0]));

    // x2 = v+, x1 = v-
    Vdir[0] = (x2[0] - x1[0]) / len;
    Vdir[1] = (x2[1] - x1[1]) / len;

    Vnorm[0] = Vdir[1];
    Vnorm[1] = -Vdir[0];

}

