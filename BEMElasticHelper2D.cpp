/*
 this file contains the basic element functions for BEM part
 */

# include "BEMElasticHelper2D.h"
# include "Eigen/Dense"

/*----------------------------------------------------------
        Apply shape function in local and global coordinate
------------------------------------------------------------*/
inline void  Cartesian(double** Elcor, double& N1, double& N2, double& N3, double* GCcor)
{
    for (int i = 0; i < 2; i++) {
        GCcor[i] = N1 * Elcor[0][i] + N2 * Elcor[1][i] + N3 * Elcor[2][i];
    }
}

/*-----------------------------------------------------------
        Serendip_func: calculate the shape function of
        the quadratic 2D element
------------------------------------------------------------*/
inline void Serendip_func(double xsi, double& N1, double& N2, double& N3, double** Ni, double* N)
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
inline void Normal_Jac(double xsi, double** Elcor, double& Jac, double* Vnorm)
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
inline double Compute_length(double xsi, double weight, double** Elcor)
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
inline void UK(double* dxr, double r, double E, double nu, double** UP)
{
    double mu, c1, conr, c, clog;
    // lame constants
    mu = E / (2.0 * (1.0 + nu));
    c1 = 3.0 - 4.0 * nu;
    c = 1.0 / (8.0 * pi * mu * (1.0 - nu));
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
inline void TK(double* dxr, double r, double E, double nu, double** TP, double* Vnorm)
{
    double c3, conr, Costh, c2;
    c3 = 1.0 - 2.0 * nu;
    c2 = 1.0 / (4.0 * pi * (1.0 - nu));
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
