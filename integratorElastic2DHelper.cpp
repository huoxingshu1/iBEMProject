
#include "configElastic2D.h"
#include "integratorElastic2D.h"
#include "integratorElastic2DHelper.h"

// derivative of rhos
/*---------------------------------------------------------------------------------
        Derivative of the ratio = a / r or its squared one
------------------------------------------------------------------------------------*/
inline double rho_der(double *x,double a, double r, int j, int flag)
{
    double result = 0.0;
    if (flag == 1)
        result = -a * x[j] * pow(r, -3);
    else if(flag == 2)
        result = -2.0 * a * a * x[j] * pow(r, -4);
    return result;
}


inline double rho_der2(double* x, double a, double r, int j, int m, int flag)
{
    double result = 0.0;
    if (flag == 1)
        result = a * (3.0 * x[m] * x[j] * pow(r, -5) - d[j][m] * pow(r, -3));
    else if (flag == 2)
        result = 2.0 * a * a * (4.0 * x[m] * x[j] * pow(r, -6) - d[j][m] * pow(r, -4));

    return result;
}

inline double rho_der3(double* x, double a, double r, int i, int j, int k)
{
    double result = 0.0;
    double term1 = (d[i][k] * x[j] + d[j][k] * x[i] + d[i][j] * x[k]) * pow(r, -6.0);
    double term2 = -6.0 * x[i] * x[j] * x[k] * pow(r, -8.0);
    result = 8.0 * a * a * (term1 + term2);
   
    return result;
}

inline double rho_der4(double* x, double a, double r, int i, int j, int k, int l)
{
    double result = 0.0;
    double term1 = pow(r, -6.0) * (d[i][k] * d[j][l] + d[j][k] * d[i][l] + d[i][j] * d[k][l]);
    double term2 = -6.0 * pow(r, -8.0) * (d[i][k] * x[j] + d[j][k] * x[i] + d[i][j] * x[k]) * x[l];
    double term3 = -6.0 * pow(r, -8.0) * (d[i][l] * x[j] * x[k] + d[j][l] * x[i] * x[k] + d[k][l] * x[i] * x[j]);
    double term4 = 48.0 * pow(r, -10.0) * x[i] * x[j] * x[k] * x[l];
    result = 8.0 * a * a * (term1 + term2 + term3 + term4);
    //result = 8.0 * pow(a, -4.0) * ((d[i][k] * d[j][l] + d[j][k] * d[i][l]) * pow(rho, 6.0)
      //  + (d[i][k] * x[j] + d[j][k] * x[i]) * 3.0 * pow(rho, 4.0) * rho_der(x, a, r, l, 2)
      //  + (d[i][l] * x[j] + d[j][l] * x[i]) * 3.0 * pow(rho, 4.0) * rho_der(x, a, r, k, 2)
      //  + x[i] * x[j] * 3.0 * (2.0 * pow(rho, 2.0) * rho_der(x, a, r, l, 2) * rho_der(x, a, r, k, 2) + 3.0 * pow(rho, 4.0) * rho_der2(x, a, r, k, l, 2))
       // ) - 2.0 * pow(a, -2.0) * d[i][j] * 2.0 * (
        //    rho_der(x, a, r, l, 2) * rho_der(x, a, r, k, 2) + pow(rho, 2.0) * rho_der2(x, a, r, k, l, 2)
        //    );
    return result;
}

inline double rho_der5(double* x, double a, double r, int i, int j, int k, int l, int s)
{
    double result = 0.0;
    double term1 = -6.0 * pow(r, -8.0) * (d[i][k] * d[j][l] + d[j][k] * d[i][l] + d[i][j] * d[k][l]) * x[s];
    double term2 = -6.0 * pow(r, -8.0) * ((d[i][k] * d[j][s] + d[j][k] * d[i][s] + d[i][j] * d[k][s]) * x[l] + d[l][s] * (d[i][k] * x[j] + d[j][k] * x[i] + d[i][j] * x[k]));
    term2 += -6.0 * (-8.0) * pow(r, -10.0) * x[l] * x[s] * (d[i][k] * x[j] + d[j][k] * x[i] + d[i][j] * x[k]);
    double term3 = -6.0 * pow(r, -8.0) * (d[i][l] * (d[j][s] * x[k] + d[k][s] * x[j]) + d[j][l] * (d[i][s] * x[k] + d[k][s] * x[i]) + d[k][l] * (d[i][s] * x[j] + d[j][s] * x[i]));
    term3 += -6.0 * (-8.0) * pow(r, -10.0) * x[s] * (d[i][l] * x[j] * x[k] + d[j][l] * x[i] * x[k] + d[k][l] * x[i] * x[j]);
    double term4 = 48.0 * pow(r, -10.0) * (d[i][s] * x[j] * x[k] * x[l] + d[j][s] * x[i] * x[k] * x[l] + d[k][s] * x[i] * x[j] * x[l] + d[l][s] * x[i] * x[j] * x[k]);
    term4 += 48.0 * (-10.0) * x[i] * x[j] * x[k] * x[l] * x[s] * pow(r, -12.0);
    result = 8.0 * a * a * (term1 + term2 + term3 + term4);
    // result = 8.0 * pow(a, -4.0) * (
        // (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * 3.0 * pow(rho, 4.0) * rho_der(x, a, r, s, 2)
        // + (d[i][k] * d[j][s] + d[j][k] * d[i][s]) * 3.0 * pow(rho, 4.0) * rho_der(x, a, r, l, 2)
        // + (d[i][k] * x[j] + d[j][k] * x[i]) * 3.0 * (2.0 * pow(rho, 2.0) * rho_der(x, a, r, s, 2) * rho_der(x, a, r, l, 2) + pow(rho, 4.0) * rho_der2(x, a, r, l, s, 2))
         //+ (d[i][l] * d[j][s] + d[j][l] * d[i][s]) * 3.0 * pow(rho, 4.0) * rho_der(x, a, r, k, 2)
         //+ (d[i][l] * x[s] + d[j][l] * x[i]) * 3.0 * (2.0 * pow(rho, 2.0) * rho_der(x, a, r, s, 2) * rho_der(x, a, r, k, 2) + pow(rho, 4.0) * rho_der2(x, a, r, k, s, 2))
        // + (d[i][s] * x[j] + d[j][s] * x[i]) * 3.0 * (2.0 * pow(rho, 2.0) * rho_der(x, a, r, l, 2) * rho_der(x, a, r, k, 2) + 3.0 * pow(rho, 4.0) * rho_der2(x, a, r, k, l, 2))
        // + 3.0 * x[i] * x[j] * (2.0 * (rho_der(x, a, r, s, 2) * rho_der(x, a, r, k, 2) * rho_der(x, a, r, l, 2) + pow(rho, 2.0) * (rho_der2(x, a, r, k, s, 2) * rho_der(x, a, r, l, 2) + rho_der2(x, a, r, l, s, 2) * rho_der(x, a, r, k, 2)))
         //    + 3.0 * (2.0 * pow(rho, 2.0) * rho_der(x, a, r, s, 2) * rho_der2(x, a, r, k, l, 2) + pow(rho, 4.0) * rho_der3(x, a, r, k, l, s)))
        // ) - 2.0 * pow(a, -2.0) * d[i][j] * 2.0 * (
         //    rho_der2(x, a, r, l, s, 2) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, k, s, 2) * rho_der(x, a, r, l, 2) + pow(rho, 2.0) * rho_der3(x, a, r, k, l, s) + rho_der(x, a, r, s, 2) * rho_der2(x, a, r, k, l, 2)
         //    );
    return result;
}

// derivatives of directional normal
/*------------------------------------------------------------------------------------
        Derivative of the directional np = x_p / |x|, 2 is second
-------------------------------------------------------------------------------------*/
inline double dir_der(double* x, double a, double r, int p, int j)
{
    // math: \partial x_j  n_p
    double result = 0.0;
    // derivatives of the n_p w.r.t x_j
    result = d[p][j] / r - x[p] * x[j] / (r * r * r);
    return result;
}

inline double dir_der2(double* x, double a, double r, int i, int j, int m)
{
    // math: n_{i,jm}
    double result = 0.0;
    result = -d[i][j] * x[m] * pow(r, -3) - (d[i][m] * x[j] + d[j][m] * x[i]) * pow(r, -3) + 3.0 * pow(r, -5) * x[i] * x[j] * x[m];
    return result;
}

// uniform Eigenstrain effects on the boundary nodes
double Eshelby_displacement(double* x, double a, int i, int k, int l)
{
    // x is the |x|, a is the radius of the inclusion, i is direciton, k,l are the dummy index
    double conp, conpp; double para, parap, parapp;
    double rho, nu, result; result = 0.0;
    nu = nu0;
    double rp = sqrt(x[0] * x[0] + x[1] * x[1]);
    rho = a / rp;

    conp = (1.0 - 2.0 * nu) / (4.0 * Pi * (1.0 - nu));
    conpp = -Pi * rho * rho;
    para = -(2.0) / (4.0 * Pi * (1.0 - nu));
    parap = -(Pi * a * rho) / (4.0);
    parapp = 4.0 * (1.0 - (rho * rho));

    double delta_kl, delta_il, delta_ik, ni, nk, nl;
    delta_kl = delta_il = delta_ik = 0.0;
    if (i == k) delta_ik = 1.0;
    if (i == l) delta_il = 1.0;
    if (k == l) delta_kl = 1.0;

    ni = x[i] / rp; nl = x[l] / rp; nk = x[k] / rp;

    result = conp * conpp * (delta_kl * x[i] - delta_il * x[k] - delta_ik * x[l]);
    result += para * (parap * (rho * rho * (delta_ik * nl + delta_il * nk + delta_kl * ni) + parapp * ni * nk * nl));

    return result;
}

// linear eigenstrain to boundary nodes
double Eshelby_displacement_linear(double* x, double a, int i, int k, int l, int p)
{
    // p is the xp' linear term
    double result = 0.0;
    // tensor in jupyter notebook is F_{iklp}
    double para = (1.0) / (8.0 * Pi * (1.0 - nu0));

    result = para * (Psi_1_3(x, a, p, i, k, l) - 2.0 * nu0 * d[k][l] * Phi_1_1(x, a, p, i) - 2.0 * (1.0 - nu0) * d[i][k] * Phi_1_1(x, a, p, l) - 2.0 * (1.0 - nu0) * d[i][l] * Phi_1_1(x, a, p, k));
    // the reason for a negative sign is due to the Eshelby tensor on the left of the matrix
    // this negative sign will be assgined at the outer functions

    return result;
}

// quadratic eigenstrain to boundary nodes
double Eshelby_displacement_quadratic(double* x, double a, int i, int k, int l, int p, int q)
{
    // p is the xp' quadratic term, q is the xq' quadratic term
    double result = 0.0;
    double para = (1.0) / (8.0 * Pi * (1.0 - nu0));

    result = para * (Psi_2_3(x, a, p, q, i, k, l)
            - 2.0 * nu0 * Phi_2_1(x, a, p, q, i) * d[k][l]
            - 2.0 * (1.0 - nu0) * Phi_2_1(x, a, p, q, k) * d[i][l]
            - 2.0 * (1.0 - nu0) * Phi_2_1(x, a, p, q, l) * d[i][k]);

    return result;
}


/*
 For uniform eigenstrains, because the functions are simple, we have directly integrated into the Eshelby_displacement function itself
 
 In the following, you can see linear and quadratic potentials for displacements
 
 */

inline double Phi_1_1(double *x,double a, int p, int i)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    if (a > r) {   // interior case
        result = 0.0;
        cout << " Check the solution on the linear displacement function Phi_1_1 due to the wrong ratio " << endl;
    }
    else {      // exterior case
        result = 0.0;
        double rho = a / r;
        double n[2] = { x[0] / r, x[1] / r };
        result = Pi * 0.5 * a * a * rho * rho * d[p][i] - Pi * a * a * rho * rho * n[i] * n[p];
    }
    return result;
}

inline double Phi_2_1(double* x, double a, int m, int n, int i)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (a >= r) {       // interior case
        result = 0.0;
        cout << "please check the overlapPing of the boundary of the inclusion" << endl;
    }

    else {          // exterior case
        double rho = a / r;
        result = pow(a, 2.0) * (pow(a, 2.0) * (d[m][i] * x[n] + d[n][i] * x[m]) * vxx_0(x, a)
            + pow(a, 2.0) * x[m] * x[n] * vxx_1(x, a, i)
            - d[m][n] * x[i] * vxa_0(x, a)

            );
    }

    return result;
}

inline double Psi_1_3(double* x, double a, int p, int i, int k, int l)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a > r) {    // interior case
        result = 0.0;
        cout << " Check the solution on the linear displacement function Psi_1_3 due to the wrong ratio " << endl;
    }
    else {      // exterior case
        double rho = a / r;
        double n[2] = { x[0] / r, x[1] / r };

        double term1 = a * a * (d[p][i] * d[k][l] + d[p][l] * d[i][k] + d[p][k] * d[i][l]) * Pi * (0.5 * pow(rho, 2.0) - pow(rho, 4.0) / 6.0);
        double term2 = (d[p][i] * x[k] + d[p][k] * x[i] + d[i][k] * x[p]) * Pi * x[l] * (-pow(rho, 4.0) + (2.0 / 3.0) * pow(rho, 6.0));
        double term3 = x[p] * x[i] * Pi * (-d[k][l] * pow(rho, 4.0) + 4.0 * n[k] * n[l] * pow(rho, 4.0));
        double term4 = (d[p][l] * x[i] + d[i][l] * x[p]) * x[k] * Pi * (-pow(rho, 4.0) + (2.0 / 3.0) * pow(rho, 6.0));
        double term5 = x[p] * x[i] * (2.0 / 3.0) * Pi * (d[k][l] * pow(rho, 6.0) - 6.0 * n[k] * n[l] * pow(rho, 6.0));
        
        result = term1 + term2 + term3 + term4 + term5;
    }
    return result;
}

inline double Psi_2_3(double* x, double a, int m, int n, int i, int j, int k)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {       // interior case
        result = 0.0;
        cout << "please check the overlapPing of the boundary of the inclusion" << endl;
    }

    else {          // exterior case
        result = (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * x[k] * vxxa_0(x, a);
        result += ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * vxxa_0(x, a);
        result += (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * vxxa_1(x, a, k);
        result += (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * vxxa_0(x, a);
        result += (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * vxxa_1(x, a, k);
        result += (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * vxxa_1(x, a, j);
        result += x[m] * x[i] * x[n] * vxxa_2(x, a, j, k);
        result *= pow(a, 4.0);

        result += pow(a, 2.0) * d[m][n] * (
                -d[i][j] * x[k] * vx2axxx_0(x, a)
                - (d[i][k] * x[j] + d[j][k] * x[i]) * vx2axxx_0(x, a)
                - x[i] * x[j] * vx2axxx_1(x, a, k));

    }
    
    return result;
}

/*------------------------------------------------------------------------------------
        Integral function of 3D cylinderical integral
        Phi_{p,ij}
--------------------------------------------------------------------------------------*/
inline double Phi_1_2(double* x, double a, int p, int i, int j)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a >= r) {   // interior case
        result = -Pi * (d[i][p] * x[j] + d[j][p] * x[i] + d[i][j] * x[p]);

    }
    else {      // exterior case
        double term1, term2, term3;
        term1 = term2 = term3 = 0.0;
        term1 = d[p][i] * rho_der(x, a, r, j, 2);
        term2 = d[p][j] * rho_der(x, a, r, i, 2);
        term3 = x[p] * rho_der2(x, a, r, i, j, 2);
        result = 0.5 * Pi * a * a * (term1 + term2 + term3);

    }
    return result;
}

inline double Phi_2_2(double* x, double a, int p,int q, int i, int j)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    int m = p; int n = q;
    result = pow(a, 2.0) * (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * vxx_0(x, a)
        + pow(a, 2.0) * (d[m][i] * x[n] + d[n][i] * x[m]) * vxx_1(x, a, j)
        + pow(a, 2.0) * (d[m][j] * x[n] + d[n][j] * x[m]) * vxx_1(x, a, i)
        + pow(a, 2.0) * x[m] * x[n] * vxx_2(x, a, i, j)
        - d[m][n] * (d[i][j] * vxa_0(x, a) + x[i] * vxa_1(x, a, j));
    
    result *= pow(a, 2.0);
    return result;
}

/*----------------------------------------------------------------------------------
        Integral function of 3D cylinderical integral
        Psi_{p,iklj}
------------------------------------------------------------------------------------*/
inline double Psi_1_4(double* x, double a, int p, int i, int j, int k, int l)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    double term1 = d[i][p] * d[j][k] * vxa_1(x, a, l) + d[i][p] * d[j][l] * vxa_1(x, a, k);
    double term2 = d[i][p] * x[j] * vxa_2(x, a, k, l) + (d[p][j] * d[i][k] + d[i][j] * d[p][k]) * vxa_1(x, a, l);
    double term3 = (d[p][j] * d[i][l] + d[i][j] * d[p][l]) * vxa_1(x, a, k) + (d[p][j] * x[i] + d[i][j] * x[p]) * vxa_2(x, a, k, l);
    double term4 = (d[p][k] * d[i][l] + d[i][k] * d[p][l]) * vxa_1(x, a, j);
    double term5 = (d[p][k] * x[i] + d[i][k] * x[p]) * vxa_2(x, a, j, l);
    double term6 = (d[p][l] * x[i] + d[i][l] * x[p]) * vxa_2(x, a, j, k) + x[p] * x[i] * vxa_3(x, a, j, k, l);
    result = a * a * (term1 + term2 + term3 + term4 + term5 + term6);

    return result;
}

inline double Psi_2_4(double* x, double a, int p,int q, int i, int j, int k, int l)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    int m = p; int n = q;
    result = pow(a, 4.0) * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * d[k][l] * vxxa_0(x, a)
        + (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * x[k] * vxxa_1(x, a, l)
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * vxxa_0(x, a)
        + ((d[m][i] * d[n][k] + d[m][k] * d[n][i]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * vxxa_1(x, a, l)
        + ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * vxxa_1(x, a, k)
        + (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * vxxa_2(x, a, k, l)
        + (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * vxxa_0(x, a)
        + (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * vxxa_1(x, a, l)
        + (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * vxxa_1(x, a, k)
        + (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * vxxa_2(x, a, k, l)
        + (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * vxxa_1(x, a, j)
        + (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * vxxa_2(x, a, j, l)
        + (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[m] * x[n]) * vxxa_2(x, a, j, k)
        + x[m] * x[n] * x[i] * vxxa_3(x, a, j, k, l)
        ) + pow(a, 2.0) * d[m][n] * (
            -d[i][j] * d[k][l] * vx2axxx_0(x, a)
            - d[i][j] * x[k] * vx2axxx_1(x, a, l)
            - (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * vx2axxx_0(x, a)
            - (d[i][k] * x[j] + d[j][k] * x[i]) * vx2axxx_1(x, a, l)
            - (d[i][l] * x[j] + d[j][l] * x[i]) * vx2axxx_1(x, a, k)
            - x[i] * x[j] * vx2axxx_2(x, a, k, l)

            );

    return result;
}

// integration for strain der
inline double Phi_der2(double* x, double a, double r, int i, int j, int m)
{
    double result = 0.0;
   // double cons = -Pi * a * a;
   // result = d[i][j] * (-2.0) * x[m] * pow(r, -4) - 2.0 * ((d[i][m] * x[j] + d[j][m] * x[i]) * pow(r, -4) - 4.0 * x[i] * x[j] * x[m] * pow(r, -6));
    double cons = -Pi;
    result = d[i][j] * rho_der(x, a, r, m, 2) + d[i][m] * rho_der(x, a, r, j, 2) + x[i] * rho_der2(x, a, r, j, m, 2);
    result *= cons;
    return result;
}

inline double Psi_der2(double* x, double a, double r, int i, int k, int l, int j, int m)
{
    double result = 0.0;    double rho = a / r;
    double n[2] = {x[0] / r, x[1] / r};
    double para1 = -Pi * a * rho_der2(x,a,r,j,m,1) / (8.0);
    double term1 = 2.0 * rho * rho * ((d[i][k] * n[l] + d[i][l] * n[k] + d[k][l] * n[i]) + 8.0 * (1.0 - pow(rho,2)) * n[i] * n[k] * n[l]);

    double para2 = -Pi * a * rho_der(x, a, r, j, 1) / (8.0);
    double term2 = 2.0 * rho_der(x, a, r, m, 2) * (d[i][k] * n[l] + d[i][l] * n[k] + d[k][l] * n[i]);
    term2 += 2.0 * rho * rho * (d[i][k] * dir_der(x, a, r, l, m) + d[i][l] * dir_der(x, a, r, k, m) + d[k][l] * dir_der(x, a, r, i, m));
    term2 -= 8.0 * rho_der(x, a, r, m, 2) * n[i] * n[k] * n[l];
    term2 += 8.0 * (1.0 - pow(rho, 2)) * (dir_der(x, a, r, i, m) * n[k] * n[l] + dir_der(x, a, r, k, m) * n[i] * n[l] + dir_der(x, a, r, l, m) * n[i] * n[k]);


    double para3 = -Pi * a * rho_der(x, a, r, m, 1) / (8.0);
    double term3 = 2.0 * rho_der(x, a, r, j, 2) * (d[i][k] * n[l] + d[i][l] * n[k] + d[k][l] * n[i]);
    term3 += 2.0 * rho * rho * (d[i][k] * dir_der(x, a, r, l, j) + d[i][l] * dir_der(x, a, r, k, j) + d[k][l] * dir_der(x, a, r, i, j));
    term3 -= 8.0 * rho_der(x, a, r, j, 2) * n[i] * n[k] * n[l];
    term3 += 8.0 * (1.0 - pow(rho, 2)) * (dir_der(x, a, r, i, j) * n[k] * n[l] + dir_der(x, a, r, k, j) * n[i] * n[l] + dir_der(x, a, r, l, j) * n[i] * n[k]);

    double para4 = -Pi * a * rho / (8.0);
    double term4 = 2.0 * rho_der2(x, a, r, j, m, 2) * (d[i][k] * n[l] + d[i][l] * n[k] + d[k][l] * n[i]);
    term4 += 2.0 * rho_der(x, a, r, j, 2) * (d[i][k] * dir_der(x, a, r, l, m) + d[i][l] * dir_der(x, a, r, k, m) + d[k][l] * dir_der(x, a, r, i, m));
    term4 += 2.0 * rho_der(x, a, r, m, 2) * (d[i][k] * dir_der(x, a, r, l, j) + d[i][l] * dir_der(x, a, r, k, j) + d[k][l] * dir_der(x, a, r, i, j));
    term4 += 2.0 * rho * rho * (d[i][k] * dir_der2(x, a, r, l, j, m) + d[i][l] * dir_der2(x, a, r, k, j, m) + d[k][l] * dir_der2(x, a, r, i, j, m));
    term4 -= 8.0 * rho_der2(x, a, r, j, m, 2) * n[i] * n[k] * n[l];
    term4 -= 8.0 * rho_der(x, a, r, j, 2) * (dir_der(x, a, r, i, m) * n[k] * n[l] + dir_der(x, a, r, k, m) * n[i] * n[l] + dir_der(x, a, r, l, m) * n[i] * n[k]);
    term4 -= 8.0 * rho_der(x, a, r, m, 2) * (dir_der(x, a, r, i, j) * n[k] * n[l] + dir_der(x, a, r, k, j) * n[i] * n[l] + dir_der(x, a, r, l, j) * n[i] * n[k]);
    term4 += 8.0 * (1.0 - pow(rho, 2)) * (dir_der2(x, a, r, i, j, m) * n[k] * n[l] + dir_der(x, a, r, i, j) * dir_der(x, a, r, k, m) * n[l] + dir_der(x, a, r, i, j) * n[k] * dir_der(x, a, r, l, m) + dir_der2(x, a, r, k, j, m) * n[i] * n[l] + dir_der(x, a, r, k, j) * dir_der(x, a, r, i, m) * n[l] + dir_der(x, a, r, k, j) * n[i] * dir_der(x, a, r, l, m) + dir_der2(x, a, r, l, j, m) * n[i] * n[k] + dir_der(x, a, r, l, j) * dir_der(x, a, r, i, m) * n[k] + dir_der(x, a, r, l, j) * n[i] * dir_der(x, a, r, k, m));
        
    result = para1 * term1 + para2 * term2 + para3 * term3 + para4 * term4;
    return result;
}


inline double Phi_1_3(double* x, double a, int p, int i, int j, int k)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {
        result = -Pi * (d[i][p] * d[j][k] + d[j][p] * d[i][k] + d[p][k] * d[i][j]);
    }
    else {
        result = 0.5 * Pi * pow(a, 2.0) * (d[i][p] * rho_der2(x, a, r, j, k, 2)
                + d[j][p] * rho_der2(x, a, r, i, k, 2)
                + d[p][k] * rho_der2(x, a, r, i, j, 2)
                + x[p] * rho_der3(x, a, r, i, j, k));
    }

    return result;
}

inline double Phi_2_3(double* x, double a, int p, int q, int i, int j, int k)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    int m = p; int n = q;
    result = pow(a, 4.0) * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * vxx_1(x, a, k) + (d[m][i] * d[n][k] + d[n][i] * d[m][k]) * vxx_1(x, a, j)
        + (d[m][j] * d[n][k] + d[n][j] * d[m][k]) * vxx_1(x, a, i) + (d[m][i] * x[n] + d[n][i] * x[m]) * vxx_2(x, a, j, k)
        + (d[m][j] * x[n] + d[n][j] * x[m]) * vxx_2(x, a, i, k) + (d[m][k] * x[n] + d[n][k] * x[m]) * vxx_2(x, a, i, j) + x[m] * x[n] * vxx_3(x, a, i, j, k))
        - pow(a, 2.0) * d[m][n] * (d[i][j] * vxa_1(x, a, k) + d[i][k] * vxa_1(x, a, j) + x[i] * vxa_2(x, a, j, k));

    return result;
}


inline double Psi_1_5(double* x, double a, int p, int i, int j, int k, int l, int s)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    double term1 = d[i][p] * d[j][k] * vxa_2(x, a, l, s);
    double term2 = d[i][p] * d[j][l] * vxa_2(x, a, k, s) + d[i][p] * d[j][s] * vxa_2(x, a, k, l) + d[i][p] * x[j] * vxa_3(x, a, k, l, s);
    double term3 = (d[p][j] * d[i][k] + d[i][j] * d[p][k]) * vxa_2(x, a, l, s) + (d[p][j] * d[i][l] + d[i][j] * d[p][l]) * vxa_2(x, a, k, s);
    double term4 = (d[p][j] * d[i][s] + d[i][j] * d[p][s]) * vxa_2(x, a, k, l) + (d[p][j] * x[i] + d[i][j] * x[p]) * vxa_3(x, a, k, l, s);
    double term5 = (d[p][k] * d[i][l] + d[i][k] * d[p][l]) * vxa_2(x, a, j, s) + (d[p][k] * d[i][s] + d[i][k] * d[p][s]) * vxa_2(x, a, j, l);
    double term6 = (d[p][k] * x[i] + d[i][k] * x[p]) * vxa_3(x, a, j, l, s) + (d[p][l] * d[i][s] + d[i][l] * d[p][s]) * vxa_2(x, a, j, k);
    double term7 = (d[p][l] * x[i] + d[i][l] * x[p]) * vxa_3(x, a, j, k, s) + (d[p][s] * x[i] + d[i][s] * x[p]) * vxa_3(x, a, j, k, l);
    double term8 = x[p] * x[i] * vxa_4(x, a, j, k, l, s);

    result = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8);
    result *= pow(a, 2.0);
    return result;
}

inline double Psi_2_5(double* x, double a, int p, int q, int i, int j, int k, int l, int s)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    int m = p; int n = q;
    result = pow(a, 4.0) * (
        (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * (d[k][l] * vxxa_1(x, a, s) + d[k][s] * vxxa_1(x, a, l) + x[k] * vxxa_2(x, a, l, s))
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * vxxa_1(x, a, s)
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][k]) * vxxa_1(x, a, l)
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * vxxa_2(x, a, l, s)
        + ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][l]) * vxxa_1(x, a, k)
        + ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * vxxa_2(x, a, k, s)
        + ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][s]) * vxxa_2(x, a, k, l)
        + (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * vxxa_3(x, a, k, l, s)
        + (d[m][s] * x[n] * x[i] + d[n][s] * x[m] * x[i] + d[i][s] * x[m] * x[n]) * vxxa_3(x, a, j, k, l) + x[m] * x[n] * x[i] * vxxa_4(x, a, j, k, l, s)
        + (d[m][l] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][l] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][l] * (d[n][s] * x[m] + d[m][s] * x[n])) * vxxa_2(x, a, j, k)
        + (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[n] * x[m]) * vxxa_3(x, a, j, k, s)
        + (d[m][k] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][k] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][k] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * vxxa_1(x, a, j)
        + (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * vxxa_2(x, a, j, s)
        + (d[m][k] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][k] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][k] * (d[m][s] * x[n] + d[n][s] * x[m])) * vxxa_2(x, a, j, l)
        + (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * vxxa_3(x, a, j, l, s)
        + (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * vxxa_1(x, a, s)
        + (d[m][j] * (d[n][k] * d[i][s] + d[i][k] * d[n][s]) + d[n][j] * (d[m][k] * d[i][s] + d[i][k] * d[m][s]) + d[i][j] * (d[m][k] * d[n][s] + d[n][k] * d[m][s])) * vxxa_1(x, a, l)
        + (d[m][j] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][j] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][j] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * vxxa_1(x, a, k)
        + (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * vxxa_2(x, a, l, s)
        + (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * vxxa_2(x, a, k, s)
        + (d[m][j] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][j] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][j] * (d[m][s] * x[n] + d[n][s] * x[m])) * vxxa_2(x, a, k, l)
        + (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * vxxa_3(x, a, k, l, s)

        ) + pow(a, 2.0) * d[m][n] * (
            -d[i][j] * d[k][l] * vx2axxx_1(x, a, s)
            - d[i][j] * d[k][s] * vx2axxx_1(x, a, l)
            - d[i][j] * x[k] * vx2axxx_2(x, a, l, s)
            - (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * vx2axxx_1(x, a, s)
            - (d[i][k] * d[j][s] + d[j][k] * d[i][s]) * vx2axxx_1(x, a, l)
            - (d[i][k] * x[j] + d[j][k] * x[i]) * vx2axxx_2(x, a, l, s)
            - (d[i][l] * d[j][s] + d[j][l] * d[i][s]) * vx2axxx_1(x, a, k)
            - (d[i][l] * x[j] + d[j][l] * x[i]) * vx2axxx_2(x, a, k, s)
            - (d[i][s] * x[j] + d[j][s] * x[i]) * vx2axxx_2(x, a, k, l)
            - x[i] * x[j] * vx2axxx_3(x, a, k, l, s)

            );

    return result;
}

inline double Phi_der3(double* x, double a, double r, int i, int j, int k, int l)
{
    double result = 0.0;
    if (a >= r) {
        result = 0.0;
    }
    else {
        result = -2.0 * Pi * (d[i][j] * rho_der2(x, a, r, k, l, 2) + rho_der3(x, a, r, j, k, l) * x[i] + rho_der2(x, a, r, j, k, 2) * d[i][l] + d[i][k] * rho_der2(x, a, r, j, l, 2));
    }
    return result;
}

inline double Phi_1_4(double* x, double a, int p, int i, int j, int k, int l)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {
        result = 0.0;
    }
    else {
        result = 0.5 * Pi * pow(a, 2.0) * (d[i][p] * rho_der3(x, a, r, j, k, l)
            + d[j][p] * rho_der3(x, a, r, i, k, l)
            + d[p][k] * rho_der3(x, a, r, i, j, l)
            + d[p][l] * rho_der3(x, a, r, i, j, k)
            + x[p] * rho_der4(x, a, r, i, j, k, l));
    }

    return result;
}

inline double Phi_2_4(double* x, double a, int p, int q, int i, int j, int k, int l)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    int m = p; int n = q;
    result = pow(a, 4.0) * ((d[m][i] * d[n][j] + d[n][i] * d[m][j]) * vxx_2(x, a, k, l) + (d[m][i] * d[n][k] + d[n][i] * d[m][k]) * vxx_2(x, a, j, l)
        + (d[m][j] * d[n][k] + d[n][j] * d[m][k]) * vxx_2(x, a, i, l) + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * vxx_2(x, a, j, k)
        + (d[m][j] * d[n][l] + d[n][j] * d[m][l]) * vxx_2(x, a, i, k) + (d[m][k] * d[n][l] + d[n][k] * d[m][l]) * vxx_2(x, a, i, j)
        + (d[m][i] * x[n] + d[n][i] * x[m]) * vxx_3(x, a, j, k, l) + (d[m][j] * x[n] + d[n][j] * x[m]) * vxx_3(x, a, i, k, l)
        + (d[m][k] * x[n] + d[n][k] * x[m]) * vxx_3(x, a, i, j, l) + (d[m][l] * x[n] + d[n][l] * x[m]) * vxx_3(x, a, i, j, k) + x[m] * x[n] * vxx_4(x, a, i, j, k, l))
        - pow(a, 2.0) * d[m][n] * (d[i][j] * vxa_2(x, a, k, l) + d[i][k] * vxa_2(x, a, j, l) + x[i] * vxa_3(x, a, j, k, l) + d[i][l] * vxa_2(x, a, j, k));

    return result;
}

inline double Psi_der3(double* x, double a, double r, int i, int j, int k, int l, int s, int t)
{
    double result = 0.0;
    result = d[i][j] * vx_4(x, a, k, l, s, t) + d[i][k] * vx_4(x, a, j, l, s, t) + d[i][l] * vx_4(x, a, j, k, s, t) + d[i][s] * vx_4(x, a, j, k, l, t) + d[i][t] * vx_4(x, a, j, k, l, s) + x[i] * vx_5(x, a, j, k, l, s, t);

    return result;
}

inline double Psi_1_6(double* x, double a, int p, int i, int j, int k, int l, int s, int q)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    double term1 = d[i][p] * d[j][k] * vxa_3(x, a, l, s, q) + d[i][p] * d[j][l] * vxa_3(x, a, k, s, q);
    double term2 = d[i][p] * d[j][s] * vxa_3(x, a, k, l, q) + d[i][p] * d[j][q] * vxa_3(x, a, k, l, s) + d[i][p] * x[j] * vxa_4(x, a, k, l, s, q);
    double term3 = (d[p][j] * d[i][k] + d[i][j] * d[p][k]) * vxa_3(x, a, l, s, q) + (d[p][j] * d[i][l] + d[i][j] * d[p][l]) * vxa_3(x, a, k, s, q);
    double term4 = (d[p][j] * d[i][s] + d[i][j] * d[p][s]) * vxa_3(x, a, k, l, q) + (d[p][j] * d[i][q] + d[i][j] * d[p][q]) * vxa_3(x, a, k, l, s);
    double term5 = (d[p][j] * x[i] + d[i][j] * x[p]) * vxa_4(x, a, k, l, s, q) + (d[p][k] * d[i][l] + d[i][k] * d[p][l]) * vxa_3(x, a, j, s, q);
    double term6 = (d[p][k] * d[s][i] + d[i][k] * d[p][s]) * vxa_3(x, a, j, l, q) + (d[p][k] * d[q][i] + d[i][k] * d[p][q]) * vxa_3(x, a, j, l, s);
    double term7 = (d[p][k] * x[i] + d[i][k] * x[p]) * vxa_4(x, a, j, l, s, q) + (d[p][l] * d[i][s] + d[i][l] * d[p][s]) * vxa_3(x, a, j, k, q);
    double term8 = (d[p][l] * d[i][q] + d[i][l] * d[p][q]) * vxa_3(x, a, j, k, s) + (d[p][l] * x[i] + d[i][l] * x[p]) * vxa_4(x, a, j, k, s, q);
    double term9 = (d[p][s] * d[q][i] + d[i][s] * d[p][q]) * vxa_3(x, a, j, k, l) + (d[p][s] * x[i] + d[i][s] * x[p]) * vxa_4(x, a, j, k, l, q);
    double term10 = (d[p][q] * x[i] + d[i][q] * x[p]) * vxa_4(x, a, j, k, l, s) + x[p] * x[i] * vxa_5(x, a, j, k, l, s, q);
    result = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10;
    result *= pow(a, 2.0);
    return result;
}

inline double Psi_2_6(double* x, double a, int p, int q, int i, int j, int k, int l, int s, int t)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    int m = p; int n = q;
    result = pow(a, 4.0) * (
        (d[m][i] * d[n][j] + d[n][i] * d[m][j]) * (d[k][l] * vxxa_2(x, a, s, t) + d[k][s] * vxxa_2(x, a, l, t) + d[k][t] * vxxa_2(x, a, l, s) + x[k] * vxxa_3(x,a,l,s,t))
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][l] + (d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][k]) * vxxa_2(x, a, s, t)
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][k]) * vxxa_2(x, a, l, t)
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][k]) * vxxa_2(x, a, l, s)
        + ((d[m][i] * d[n][k] + d[n][i] * d[m][k]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][k]) * vxxa_3(x, a, l, s, t)
        + ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][s] + (d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][l]) * vxxa_2(x, a, k, t)
        + ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][l]) * vxxa_2(x, a, k, s)
        + ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * d[j][t] + (d[m][i] * d[n][t] + d[n][i] * d[m][t]) * d[j][s]) * vxxa_2(x, a, k, l)
        + ((d[m][i] * d[n][l] + d[n][i] * d[m][l]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][l]) * vxxa_3(x,a,k,s,t)
        + ((d[m][i] * d[n][s] + d[n][i] * d[m][s]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][s]) * vxxa_3(x,a,k,l,t)
        + ((d[m][i] * d[n][t] + d[n][i] * d[m][t]) * x[j] + (d[m][i] * x[n] + d[n][i] * x[m]) * d[j][t]) * vxxa_3(x,a,k,l,s)
        + (d[m][i] * x[n] + d[n][i] * x[m]) * x[j] * vxxa_4(x, a, k, l, s, t)
        + (d[m][s] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][s] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][s] * (d[m][t] * x[n] + d[n][t] * x[m])) * vxxa_3(x,a,j,k,l)
        + (d[m][s] * x[n] * x[i] + d[n][s] * x[m] * x[i] + d[i][s] * x[m] * x[n]) * vxxa_4(x, a, j, k, l, t) + x[m] * x[n] * x[i] * vxxa_5(x, a, j, k, l, s, t)
        + (d[m][t] * x[n] * x[i] + d[n][t] * x[m] * x[i] + d[i][t] * x[m] * x[n]) * vxxa_4(x,a,j,k,l,s)
        + (d[m][l] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][l] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][l] * (d[n][s] * d[m][t] + d[m][s] * d[n][t])) * vxxa_2(x,a,j,k)
        + (d[m][l] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][l] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][l] * (d[n][s] * x[m] + d[m][s] * x[n])) * vxxa_3(x, a, j, k, t)
        + (d[m][l] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][l] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][l] * (d[n][t] * x[m] + d[m][t] * x[n])) * vxxa_3(x,a,j,k,s)
        + (d[m][l] * x[n] * x[i] + d[n][l] * x[m] * x[i] + d[i][l] * x[n] * x[m]) * vxxa_4(x, a, j, k, s, t)
        + (d[m][k] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][k] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][k] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * vxxa_2(x, a, j, t)
        + (d[m][k] * (d[n][l] * d[i][t] + d[i][l] * d[n][t]) + d[n][k] * (d[m][l] * d[i][t] + d[i][l] * d[m][t]) + d[i][k] * (d[m][l] * d[n][t] + d[n][l] * d[m][t])) * vxxa_2(x, a, j, s)
        + (d[m][k] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][k] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][k] * (d[m][l] * x[n] + d[n][l] * x[m])) * vxxa_3(x,a,j,s,t)
        + (d[m][k] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][k] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][k] * (d[m][s] * d[n][t] + d[n][s] * d[m][t])) * vxxa_2(x, a, j, l)
        + (d[m][k] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][k] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][k] * (d[m][s] * x[n] + d[n][s] * x[m])) * vxxa_3(x,a,j,l,t)
        + (d[m][k] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][k] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][k] * (d[m][t] * x[n] + d[n][t] * x[m])) * vxxa_3(x,a,j,l,s)
        + (d[m][k] * x[n] * x[i] + d[n][k] * x[m] * x[i] + d[i][k] * x[m] * x[n]) * vxxa_4(x, a, j, l, s, t)
        + (d[m][j] * (d[n][k] * d[i][l] + d[i][k] * d[n][l]) + d[n][j] * (d[m][k] * d[i][l] + d[i][k] * d[m][l]) + d[i][j] * (d[m][k] * d[n][l] + d[n][k] * d[m][l])) * vxxa_2(x, a, s, t)
        + (d[m][j] * (d[n][k] * d[i][s] + d[i][k] * d[n][s]) + d[n][j] * (d[m][k] * d[i][s] + d[i][k] * d[m][s]) + d[i][j] * (d[m][k] * d[n][s] + d[n][k] * d[m][s])) * vxxa_2(x, a, l, t)
        + (d[m][j] * (d[n][l] * d[i][s] + d[i][l] * d[n][s]) + d[n][j] * (d[m][l] * d[i][s] + d[i][l] * d[m][s]) + d[i][j] * (d[m][l] * d[n][s] + d[n][l] * d[m][s])) * vxxa_2(x, a, k, t)
        + (d[m][j] * (d[n][k] * d[i][t] + d[i][k] * d[n][t]) + d[n][j] * (d[m][k] * d[i][t] + d[i][k] * d[m][t]) + d[i][j] * (d[m][k] * d[n][t] + d[n][k] * d[m][t])) * vxxa_2(x, a, l, s)
        + (d[m][j] * (d[n][k] * x[i] + d[i][k] * x[n]) + d[n][j] * (d[m][k] * x[i] + d[i][k] * x[m]) + d[i][j] * (d[m][k] * x[n] + d[n][k] * x[m])) * vxxa_3(x, a, l, s, t)
        + (d[m][j] * (d[n][l] * x[i] + d[i][l] * x[n]) + d[n][j] * (d[m][l] * x[i] + d[i][l] * x[m]) + d[i][j] * (d[m][l] * x[n] + d[n][l] * x[m])) * vxxa_3(x, a, k, s, t)
        + (d[m][j] * (d[n][l] * d[i][t] + d[i][l] * d[n][t]) + d[n][j] * (d[m][l] * d[i][t] + d[i][l] * d[m][t]) + d[i][j] * (d[m][l] * d[n][t] + d[n][l] * d[m][t])) * vxxa_2(x, a, k, s)
        + (d[m][j] * (d[n][s] * d[i][t] + d[i][s] * d[n][t]) + d[n][j] * (d[m][s] * d[i][t] + d[i][s] * d[m][t]) + d[i][j] * (d[m][s] * d[n][t] + d[n][s] * d[m][t])) * vxxa_2(x, a, k, l)
        + (d[m][j] * (d[n][s] * x[i] + d[i][s] * x[n]) + d[n][j] * (d[m][s] * x[i] + d[i][s] * x[m]) + d[i][j] * (d[m][s] * x[n] + d[n][s] * x[m])) * vxxa_3(x, a, k, l, t)
        + (d[m][j] * (d[n][t] * x[i] + d[i][t] * x[n]) + d[n][j] * (d[m][t] * x[i] + d[i][t] * x[m]) + d[i][j] * (d[m][t] * x[n] + d[n][t] * x[m])) * vxxa_3(x,a,k,l,s)
        + (d[m][j] * x[n] * x[i] + d[n][j] * x[m] * x[i] + d[i][j] * x[m] * x[n]) * vxxa_4(x, a, k, l, s, t)

        ) + pow(a, 2.0) * d[m][n] * (
            -d[i][j] * d[k][l] * vx2axxx_2(x, a, s, t)
            - d[i][j] * d[k][s] * vx2axxx_2(x, a, l, t)
            - d[i][j] * d[k][t] * vx2axxx_2(x,a,l,s)
            - d[i][j] * x[k] * vx2axxx_3(x, a, l, s, t)
            - (d[i][k] * d[j][l] + d[j][k] * d[i][l]) * vx2axxx_2(x, a, s, t)
            - (d[i][k] * d[j][s] + d[j][k] * d[i][s]) * vx2axxx_2(x, a, l, t)
            - (d[i][k] * d[j][t] + d[j][k] * d[i][t]) * vx2axxx_2(x, a, l, s)
            - (d[i][k] * x[j] + d[j][k] * x[i]) * vx2axxx_3(x,a,l,s,t)
            - (d[i][l] * d[j][s] + d[j][l] * d[i][s]) * vx2axxx_2(x, a, k, t)
            - (d[i][l] * d[j][t] + d[j][l] * d[i][t]) * vx2axxx_2(x, a, k, s)
            - (d[i][l] * x[j] + d[j][l] * x[i]) * vx2axxx_3(x,a,k,s,t)
            - (d[i][s] * d[j][t] + d[j][s] * d[i][t]) * vx2axxx_2(x,a,k,l)
            - (d[i][s] * x[j] + d[j][s] * x[i]) * vx2axxx_3(x, a, k, l, t)
            - (d[i][t] * x[j] + d[j][t] * x[i]) * vx2axxx_3(x, a, k, l, s)
            - x[i] * x[j] * vx2axxx_4(x, a, k, l, s, t)

            );

    return result;
}

// define functions by V potentials in Mura's book
/*-----------------------------------------------------------------------
     Essential functions for V_xx, orginal
------------------------------------------------------------------------*/
inline double vxx_0(double* x, double a)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (a >= r) {
        result = 0.5 * Pi * (pow(a, -2.0) - (2.0 / 3.0) * pow(a, -4.0) * pow(r, 2.0));
    }
    else {
        result = 0.5 * Pi * (1.0 / 3.0) * pow(a, 2.0) * pow(r, -4.0);
    }

    return result;
}

/*-----------------------------------------------------------------------
     Essential functions for V_xx, 1st derivative
------------------------------------------------------------------------*/
inline double vxx_1(double* x, double a, int i)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (a >= r) {
        result = -Pi * (2.0 / 3.0) * pow(a, -4.0) * x[i];
    }
    else {
        double rho = a / r;
        result = Pi * (1.0 / 3.0) * pow(a, -2.0) * pow(rho, 2.0) * rho_der(x, a, r, i, 2);
    }

    return result;
}

/*-----------------------------------------------------------------------
     Essential functions for V_xx, 2nd derivative
------------------------------------------------------------------------*/
inline double vxx_2(double* x, double a, int i, int j)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (a >= r) {
        result = -Pi * (2.0 / 3.0) * pow(a, -4.0) * d[i][j];
    }
    else {
        double rho = a / r;
        result = Pi * (1.0 / 3.0) * pow(a, -2.0) * (rho_der(x,a,r,j,2) * rho_der(x, a, r, i, 2) + rho_der2(x,a,r,i,j,2) * rho * rho);
    }

    return result;
}

/*-----------------------------------------------------------------------
     Essential functions for V_xx, 3rd derivative
------------------------------------------------------------------------*/
inline double vxx_3(double* x, double a, int i, int j, int k)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (a >= r) {
        result = 0.0;
    }
    else {
        double rho = a / r;
        result = Pi * (1.0 / 3.0) * pow(a, -2.0) * (rho_der2(x, a, r, j, k, 2) * rho_der(x, a, r, i, 2)
            + rho_der(x, a, r, j, 2) * rho_der2(x, a, r, i, k, 2)
            + rho_der3(x, a, r, i, j, k) * rho * rho
            + rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2));
    }

    return result;
}

/*-----------------------------------------------------------------------
     Essential functions for V_xx, 4th derivative
------------------------------------------------------------------------*/
inline double vxx_4(double* x, double a, int i, int j, int k, int l)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);

    if (a >= r) {
        result = 0.0;
    }
    else {
        double rho = a / r;
        result = Pi * (1.0 / 3.0) * pow(a, -2.0) * (rho_der3(x, a, r, j, k, l) * rho_der(x, a, r, i, 2)
            + rho_der2(x, a, r, j, k, 2) * rho_der2(x, a, r, i, l, 2)
            + rho_der2(x, a, r, j, l, 2) * rho_der2(x, a, r, i, k, 2)
            + rho_der(x, a, r, j, 2) * rho_der3(x, a, r, i, k, l)
            + rho_der4(x, a, r, i, j, k, l) * rho * rho
            + rho_der3(x, a, r, i, j, k) * rho_der(x, a, r, l, 2)
            + rho_der3(x, a, r, i, j, l) * rho_der(x, a, r, k, 2)
            + rho_der2(x, a, r, i, j, 2) * rho_der2(x, a, r, k, l, 2));
    }

    return result;
}

/*------------------------------------------------------------------------
        Essential functions for the V_xx - a^2 Vxxx, orginal
---------------------------------------------------------------------------*/
inline double vxxa_0(double* x, double a)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = Pi * (1.0 / 6.0) * (pow(a, -2.0) - 0.5 * pow(a, -4.0) * pow(r, 2.0));
    }
    else {      // exterior case
        result = Pi * (1.0 / 6.0) * (pow(a, 2.0) * pow(r, -4.0) - 0.5 * pow(a, 4.0) * pow(r, -6.0));
    }
    return result;
}

/*-----------------------------------------------------------------------------
        Essential functions for the V_xx - a^2 Vxxx, 1^st derivative
---------------------------------------------------------------------------------*/
inline double vxxa_1(double* x, double a, int i)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = -Pi * (1.0 / 6.0) * x[i] * pow(a, -4.0);
    }
    else {      // exterior case
        double rho = a / r;
       
        result = (1.0 / 12.0) * Pi * pow(a, -2.0) * rho_der(x, a, r, i, 2) * (4.0 * pow(rho, 2.0) - 3.0 * pow(rho, 4.0));
    }
    return result;
}

/*---------------------------------------------------------------------------
       Essential functions for the V_xx - a^2 Vxxx, 2^nd derivative
-----------------------------------------------------------------------------*/
inline double vxxa_2(double* x, double a, int i, int j)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = -Pi * (1.0 / 6.0) * d[i][j] * pow(a, -4.0);
    }
    else {      // exterior case

        double rho = a / r;
        result = (1.0 / 12.0) * Pi * pow(a, -2.0) * (
            rho_der2(x, a, r, i, j, 2) * (4.0 * pow(rho, 2.0) - 3.0 * pow(rho, 4.0))
            + rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2) * (4.0 - 6.0 * pow(rho, 2.0))
            );
    }
    return result;
}

/*---------------------------------------------------------------------------
       Essential functions for the V_xx - a^2 Vxxx, 3^rd derivative
-----------------------------------------------------------------------------*/
inline double vxxa_3(double* x, double a, int i, int j, int k)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = 0.0;
    }
    else {      // exterior case
       
        double rho = a / r;
        result = (1.0 / 12.0) * Pi * pow(a, -2.0) * (
            rho_der3(x, a, r, i, j, k) * (4.0 * pow(rho, 2.0) - 3.0 * pow(rho, 4.0))
            + rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2) * (4.0 - 6.0 * pow(rho, 2.0))
            + (rho_der2(x, a, r, i, k, 2) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, j, k, 2) * rho_der(x, a, r, i, 2)) * (4.0 - 6.0 * pow(rho, 2.0))
            - 6.0 * rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2) * rho_der(x, a, r, k, 2)
            );
    }
    return result;
}

/*---------------------------------------------------------------------------
       Essential functions for the V_xx - a^2 Vxxx, 4^th derivative
-----------------------------------------------------------------------------*/
inline double vxxa_4(double* x, double a, int i, int j, int k, int l)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = 0.0;
    }
    else {      // exterior case
      
        double rho = a / r;
        result = (1.0 / 12.0) * Pi * pow(a, -2.0) * (
            rho_der4(x, a, r, i, j, k, l) * (4.0 * pow(rho, 2.0) - 3.0 * pow(rho, 4.0))
            + rho_der3(x, a, r, i, j, k) * rho_der(x, a, r, l, 2) * (4.0 - 6.0 * pow(rho, 2.0))
            + (rho_der3(x, a, r, i, j, l) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, i, j, 2) * rho_der2(x, a, r, k, l, 2)) * (4.0 - 6.0 * pow(rho, 2.0))
            - 6.0 * rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2) * rho_der(x, a, r, l, 2)
            + (rho_der3(x, a, r, i, k, l) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, i, k, 2) * rho_der2(x, a, r, j, l, 2) + rho_der2(x, a, r, i, l, 2) * rho_der2(x, a, r, j, k, 2) + rho_der3(x, a, r, j, k, l) * rho_der(x, a, r, i, 2)) * (4.0 - 6.0 * pow(rho, 2.0))
            - 6.0 * (rho_der2(x, a, r, i, k, 2) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, j, k, 2) * rho_der(x, a, r, i, 2)) * rho_der(x, a, r, l, 2)
            - 6.0 * (rho_der2(x, a, r, i, l, 2) * rho_der(x, a, r, j, 2) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, j, l, 2) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, k, l, 2) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2))
            );
    }
    return result;
}

/*---------------------------------------------------------------------------
       Essential functions for the V_xx - a^2 Vxxx, 5^th derivative
-----------------------------------------------------------------------------*/
inline double vxxa_5(double* x, double a, int i, int j, int k, int l, int m)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = 0.0;
    }
    else {      // exterior case
        
        double rho = a / r;
        result = (1.0 / 12.0) * Pi * pow(a, -2.0) * (
            rho_der5(x, a, r, i, j, k, l, m) * (4.0 * pow(rho, 2.0) - 3.0 * pow(rho, 4.0))
            + rho_der4(x, a, r, i, j, k, l) * rho_der(x, a, r, m, 2) * (4.0 - 6.0 * pow(rho, 2.0))
            + (rho_der4(x, a, r, i, j, k, m) * rho_der(x, a, r, l, 2) + rho_der3(x, a, r, i, j, k) * rho_der2(x, a, r, l, m, 2)) * (4.0 - 6.0 * pow(rho, 2.0))
            - 6.0 * rho_der3(x, a, r, i, j, k) * rho_der(x, a, r, l, 2) * rho_der(x, a, r, m, 2)
            + (rho_der4(x, a, r, i, j, l, m) * rho_der(x, a, r, k, 2) + rho_der3(x, a, r, i, j, l) * rho_der2(x, a, r, k, m, 2) + rho_der3(x, a, r, i, j, m) * rho_der2(x, a, r, k, l, 2) + rho_der2(x, a, r, i, j, 2) * rho_der3(x, a, r, k, l, m)) * (4.0 - 6.0 * pow(rho, 2.0))
            - 6.0 * (rho_der3(x, a, r, i, j, l) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, i, j, 2) * rho_der2(x, a, r, k, l, 2)) * rho_der(x, a, r, m, 2)
            - 6.0 * (rho_der3(x, a, r, i, j, m) * rho_der(x, a, r, k, 2) * rho_der(x, a, r, l, 2) + rho_der2(x, a, r, i, j, 2) * rho_der2(x, a, r, k, m, 2) * rho_der(x, a, r, l, 2) + rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2) * rho_der2(x, a, r, l, m, 2))
            + (rho_der4(x, a, r, i, k, l, m) * rho_der(x, a, r, j, 2) + rho_der3(x, a, r, i, k, l) * rho_der2(x, a, r, j, m, 2) + rho_der3(x, a, r, i, k, m) * rho_der2(x, a, r, j, l, 2) + rho_der2(x, a, r, i, k, 2) * rho_der3(x, a, r, j, l, m) + rho_der3(x, a, r, i, l, m) * rho_der2(x, a, r, j, k, 2) + rho_der2(x, a, r, i, l, 2) * rho_der3(x, a, r, j, k, m) + rho_der4(x, a, r, j, k, l, m) * rho_der(x, a, r, i, 2) + rho_der3(x, a, r, j, k, l) * rho_der2(x, a, r, i, m, 2)) * (4.0 - 6.0 * pow(rho, 2.0))
            - 6.0 * (rho_der3(x, a, r, i, k, l) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, i, k, 2) * rho_der2(x, a, r, j, l, 2) + rho_der2(x, a, r, i, l, 2) * rho_der2(x, a, r, j, k, 2) + rho_der3(x, a, r, j, k, l) * rho_der(x, a, r, i, 2)) * rho_der(x, a, r, m, 2)
            - 6.0 * (rho_der3(x, a, r, i, k, m) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, i, k, 2) * rho_der2(x, a, r, j, m, 2) + rho_der3(x, a, r, j, k, m) * rho_der(x, a, r, i, 2) + rho_der2(x, a, r, j, k, 2) * rho_der2(x, a, r, i, m, 2)) * rho_der(x, a, r, l, 2)
            - 6.0 * (rho_der2(x, a, r, i, k, 2) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, j, k, 2) * rho_der(x, a, r, i, 2)) * rho_der2(x, a, r, l, m, 2)
            - 6.0 * (rho_der3(x, a, r, i, l, m) * rho_der(x, a, r, j, 2) * rho_der(x, a, r, k, 2)
                + rho_der2(x, a, r, i, l, 2) * (rho_der2(x, a, r, j, m, 2) * rho_der(x, a, r, k, 2) + rho_der(x, a, r, j, 2) * rho_der2(x, a, r, k, m, 2))
                + rho_der3(x, a, r, j, l, m) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, k, 2)
                + rho_der2(x, a, r, j, l, 2) * (rho_der2(x, a, r, i, m, 2) * rho_der(x, a, r, k, 2) + rho_der(x, a, r, i, 2) * rho_der2(x, a, r, k, m, 2))
                + rho_der3(x, a, r, k, l, m) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2)
                + rho_der2(x, a, r, k, l, 2) * (rho_der2(x, a, r, i, m, 2) * rho_der(x, a, r, j, 2) + rho_der(x, a, r, i, 2) * rho_der2(x, a, r, j, m, 2)))
            );
    }
    return result;
}

/*---------------------------------------------------------------------------
    Essential functions for the V_x - 2a^2 V_xx + a^4 V_xxx, orginal
----------------------------------------------------------------------------*/
inline double vx2axxx_0(double* x, double a)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = Pi * ((1.0 / 3.0) - (1.0 / 12.0) * pow(a, -2.0) * pow(r, 2.0));
    }
    else {      // exterior case
        double rho = a / r;
        result = Pi * (0.5 * pow(rho, 2.0)
                 - (1.0 / 3.0) * pow(rho, 4.0)
                + (1.0 / 12.0) * pow(rho, 6.0));
    }
    return result;
}

/*--------------------------------------------------------------------------
    Essential functions for the V_x - 2a^2 V_xx + a^4 V_xxx, 1st derivative
----------------------------------------------------------------------------*/
inline double vx2axxx_1(double* x, double a, int i)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = -Pi * (1.0 / 6.0) * x[i] * pow(a, -2.0);
    }
    else {      // exterior case
        double rho = a / r;
        result = Pi * (0.5 * rho_der(x, a, r, i, 2)
            - (2.0 / 3.0) * rho * rho * rho_der(x, a, r, i, 2)
            + (1.0 / 4.0) * pow(rho, 4.0) * rho_der(x, a, r, i, 2));
    }
    return result;
}

/*--------------------------------------------------------------------------
    Essential functions for the V_x - 2a^2 V_xx + a^4 V_xxx, 2nd derivative
----------------------------------------------------------------------------*/
inline double vx2axxx_2(double* x, double a, int i, int j)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = -Pi * (1.0 / 6.0) * d[i][j] * pow(a, -2.0);
    }
    else {      // exterior case
        double rho = a / r;
        result = Pi * (
            rho_der2(x, a, r, i, j, 2) * (0.5 - (2.0 / 3.0) * rho * rho + 0.25 * pow(rho, 4.0))
            + rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2) * ((-2.0 / 3.0) + 0.5 * rho * rho)
            );
    }
    return result;
}

/*--------------------------------------------------------------------------
    Essential functions for the V_x - 2a^2 V_xx + a^4 V_xxx, 3rd derivative
----------------------------------------------------------------------------*/
inline double vx2axxx_3(double* x, double a, int i, int j, int k)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = 0.0;
    }
    else {      // exterior case
        double rho = a / r;
        result = Pi * (
            rho_der3(x, a, r, i, j, k) * (0.5 - (2.0 / 3.0) * rho * rho + 0.25 * pow(rho, 4.0))
            + rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2) * ((-2.0 / 3.0) + 0.5 * rho * rho)
            + (rho_der2(x, a, r, i, k, 2) * rho_der(x, a, r, j, 2) + rho_der(x, a, r, i, 2) * rho_der2(x, a, r, j, k, 2)) * ((-2.0 / 3.0) + 0.5 * rho * rho)
            + 0.5 * rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2) * rho_der(x, a, r, k, 2)
            );
    }
    return result;
}

/*--------------------------------------------------------------------------
    Essential functions for the V_x - 2a^2 V_xx + a^4 V_xxx, 4^st derivative
----------------------------------------------------------------------------*/
inline double vx2axxx_4(double* x, double a, int i, int j, int k, int l)
{
    double result = 0.0;
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {   // interior case
        result = -Pi * (1.0 / 6.0) * x[i] * pow(a, -2.0);
    }
    else {      // exterior case
        double rho = a / r;

        result = Pi * (
            rho_der4(x, a, r, i, j, k, l) * (0.5 - (2.0 / 3.0) * rho * rho + 0.25 * pow(rho, 4.0))
            + rho_der3(x, a, r, i, j, k) * ((-2.0 / 3.0) + 0.5 * rho * rho) * rho_der(x, a, r, l, 2)
            + (rho_der3(x, a, r, i, j, l) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, i, j, 2) * rho_der2(x, a, r, k, l, 2)) * ((-2.0 / 3.0) + 0.5 * rho * rho)
            + 0.5 * rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2) * rho_der(x, a, r, l, 2)
            + (rho_der3(x, a, r, i, k, l) * rho_der(x, a, r, j, 2) + rho_der2(x, a, r, i, k, 2) * rho_der2(x, a, r, j, l, 2) + rho_der2(x, a, r, i, l, 2) * rho_der2(x, a, r, j, k, 2) + rho_der(x, a, r, i, 2) * rho_der3(x, a, r, j, k, l)) * ((-2.0 / 3.0) + 0.5 * rho * rho)
            + 0.5 * (rho_der2(x, a, r, i, k, 2) * rho_der(x, a, r, j, 2) + rho_der(x, a, r, i, 2) * rho_der2(x, a, r, j, k, 2)) * rho_der(x, a, r, l, 2)
            + 0.5 * (rho_der2(x, a, r, i, l, 2) * rho_der(x, a, r, j, 2) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, j, l, 2) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, k, 2) + rho_der2(x, a, r, k, l, 2) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2))
            );
    }
    return result;
}

// definetion of Vxa potentials for harmonic ones
inline double vxa_0(double* x, double a)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a >= r) {
        result = Pi * (0.5 - (1.0 / 6.0) * pow(a, -2.0) * pow(r, 2.0));
    }
    else {
        double rho = a / r;
        result = Pi * (0.5 * pow(rho, 2.0) - (1.0 / 6.0) * pow(rho, 4.0));
    }
    return result;

}

/*--------------------------------------------------------------------------------
        Essential derivative of V_x - a^2 V_xx: first order
----------------------------------------------------------------------------------*/
inline double vxa_1(double* x, double a, int i)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;
    
    if (a >= r) {
        result = (-Pi / 3.0) * pow(a, -2.0) * x[i];
    }
    else {
        double rho = a / r;
        result = Pi * rho_der(x, a, r, i, 2) * (0.5 - (1.0 / 3.0) * pow(rho, 2.0));
    }
    return result;
}

/*--------------------------------------------------------------------------------
        Essential derivative of V_x - a^2 V_xx: second order
----------------------------------------------------------------------------------*/
inline double vxa_2(double* x, double a, int i, int j)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a >= r) {
        result = (-Pi / 3.0) * pow(a, -2.0) * d[i][j];
    }
    else {
        double rho = a / r;
        result = Pi * (rho_der2(x, a, r, i, j, 2) * (0.5 - (1.0 / 3.0) * pow(rho, 2.0))
                - (1.0 / 3.0) * rho_der(x, a, r, i, 2) * rho_der(x, a, r, j, 2));

    }
    return result;
}

/*--------------------------------------------------------------------------------
        Essential derivative of V_x - a^2 V_xx: Third order
----------------------------------------------------------------------------------*/
inline double vxa_3(double* x, double a, int i, int j, int k)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a >= r) {
        result = 0.0;
    }
    else {
        double rho = a / r;
        result = Pi * rho_der3(x, a, r, i, j, k) * (0.5 - (1.0 / 3.0) * pow(rho, 2.0));
        result += Pi * (-1.0 / 3.0) * (rho_der2(x, a, r, i, j, 2) * rho_der(x, a, r, k, 2)
                + rho_der2(x, a, r, i, k, 2) * rho_der(x, a, r, j, 2)
                + rho_der2(x, a, r, j, k, 2) * rho_der(x, a, r, i, 2));
    }
    return result;
}

/*--------------------------------------------------------------------------------
        Essential derivative of V_x - a^2 V_xx: Fourth order
----------------------------------------------------------------------------------*/
inline double vxa_4(double* x, double a, int i, int j, int k, int l)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a >= r) {
        result = 0.0;
    }
    else {
        double rho = a / r;
        result = Pi * rho_der4(x, a, r, i, j, k,l) * (0.5 - (1.0 / 3.0) * pow(rho, 2.0));
        result += Pi * (-1.0 / 3.0) * (
            rho_der3(x, a, r, i, j, k) * rho_der(x, a, r, l, 2)
            + rho_der3(x, a, r, i, j, l) * rho_der(x, a, r, k, 2)
            + rho_der2(x, a, r, i, j, 2) * rho_der2(x, a, r, k, l, 2)
            + rho_der3(x, a, r, i, k, l) * rho_der(x, a, r, j, 2)
            + rho_der2(x, a, r, i, k, 2) * rho_der2(x, a, r, j, l, 2)
            + rho_der3(x, a, r, j, k, l) * rho_der(x, a, r, i, 2)
            + rho_der2(x, a, r, j, k, 2) * rho_der2(x, a, r, i, l, 2));
    }
    return result;
}

/*--------------------------------------------------------------------------------
        Essential derivative of V_x - a^2 V_xx: Fourth order
----------------------------------------------------------------------------------*/
inline double vxa_5(double* x, double a, int i, int j, int k, int l, int s)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double result = 0.0;

    if (a >= r) {
        result = 0.0;
    }
    else {
        double rho = a / r;
        result = Pi * rho_der5(x, a, r, i, j, k, l, s) * (0.5 - (1.0 / 3.0) * pow(rho, 2.0));
        result += Pi * (-1.0 / 3.0) * (
            rho_der4(x, a, r, i, j, k, s) * rho_der(x, a, r, l, 2)
            + rho_der3(x, a, r, i, j, k) * rho_der2(x, a, r, l, s, 2)
            + rho_der4(x, a, r, i, j, k, l) * rho_der(x, a, r, s, 2)
            + rho_der4(x, a, r, i, j, l, s) * rho_der(x, a, r, k, 2)
            + rho_der3(x, a, r, i, j, l) * rho_der2(x, a, r, k, s, 2)
            + rho_der3(x, a, r, i, j, s) * rho_der2(x, a, r, k, l, 2)
            + rho_der2(x, a, r, i, j, 2) * rho_der3(x, a, r, k, l, s)
            + rho_der4(x, a, r, i, k, l, s) * rho_der(x, a, r, j, 2)
            + rho_der3(x, a, r, i, k, l) * rho_der2(x, a, r, j, s, 2)
            + rho_der3(x, a, r, i, k, s) * rho_der2(x, a, r, j, l, 2)
            + rho_der2(x, a, r, i, k, 2) * rho_der3(x, a, r, j, l, s)
            + rho_der4(x, a, r, j, k, l, s) * rho_der(x, a, r, i, 2)
            + rho_der3(x, a, r, j, k, l) * rho_der2(x, a, r, i, s, 2)
            + rho_der3(x, a, r, j, k, s) * rho_der2(x, a, r, i, l, 2))
            + rho_der2(x, a, r, j, k, 2) * rho_der3(x, a, r, i, l, s);
    }
    return result;
}

/*-----------------------------------------------------------------------------------
        Essential derivative of the V - a^2 V_x for uniform term, fourth derivative
-------------------------------------------------------------------------------------*/
inline double vx_4(double* x, double a, int i, int j, int k, int l)
{
    double result, r;
    r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {
        result = 0.0;
    }
    else {
        result = -2.0 * Pi * (d[i][j] * rho_der2(x, a, r, k, l, 2) + rho_der2(x, a, r, j, l, 2) * d[i][k] + rho_der3(x, a, r, j, k, l) * x[i] + d[i][l] * rho_der2(x, a, r, j, k,2));
        result += -0.5 * Pi * a * a * rho_der4(x, a, r, i, j, k, l);
    }
    return result;
}

/*-----------------------------------------------------------------------------------
        Essential derivative of the V - a^2 V_x for uniform term, fifth derivative
-------------------------------------------------------------------------------------*/
inline double vx_5(double* x, double a, int i, int j, int k, int l, int s)
{
    double result, r;
    r = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (a >= r) {
        result = 0.0;
    }
    else {
        result = -2.0 * Pi * (d[i][j] * rho_der3(x, a, r, k, l, s) + rho_der3(x, a, r, j, l, s) * d[i][k] + rho_der4(x, a, r, j, k, l, s) * x[i] + d[i][s] * rho_der3(x, a, r, j, k, l) + d[i][l] * rho_der3(x,a,r,j,k,s));
        result += -0.5 * Pi * a * a * rho_der5(x, a, r, i, j, k, l, s);
    }
    return result;
}

/*
 // Eshelby tensor integration results by uniform strain
 */
double Eshelby_strain_integration(int i, int j, int k, int l, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);

    double a = radius[sec];
    double rp = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    double result = 0.0;

    double delta_ij, delta_ik, delta_il, delta_jk, delta_jl, delta_kl;
    delta_ij = delta_ik = delta_il = delta_jk = delta_jl = delta_kl = 0.0;

    if (i == j) delta_ij = 1.0;
    if (i == k) delta_ik = 1.0;
    if (i == l) delta_il = 1.0;
    if (j == k) delta_jk = 1.0;
    if (j == l) delta_jl = 1.0;
    if (k == l) delta_kl = 1.0;

    // judge the integration range, the Eshelby's tensor is always about the metrix material
    double nu = nu0; double mu = mu_0;
    if (a >= rp) {                  // interior integration
        double consp = (4.0 * nu - 1.0) / (8.0 * (1.0 - nu));
        double conspp = (3.0 - 4.0 * nu) / (8.0 * (1.0 - nu));
        result = consp * delta_ij * delta_kl + conspp * (delta_ik * delta_jl + delta_il * delta_jk);
    }
    else {                          // exterior integration
        double dxr[2]; dxr[0] = x1[0] / rp; dxr[1] = x1[1] / rp;

        double consp = (a * a) / (8.0 * (1.0 - nu) * pow(rp, 4));
        double term1 = (4.0 * nu * rp * rp - 2.0 * rp * rp + a * a) * delta_ij * delta_kl;
        double term2 = (-4.0 * nu * rp * rp + 2.0 * rp * rp + a * a) * (delta_ik * delta_jl + delta_il * delta_jk);
        double term3 = -4.0 * (2.0 * nu * rp * rp - rp * rp + a * a) * delta_kl * dxr[i] * dxr[j];
        double term4 = 4.0 * (rp * rp - a * a) * delta_ij * dxr[k] * dxr[l];
        double term5 = 4.0 * (nu * rp * rp - a * a) * (delta_ik * dxr[j] * dxr[l] + delta_il * dxr[j] * dxr[k] + delta_jk * dxr[i] * dxr[l] + delta_jl * dxr[i] * dxr[k]);
        double term6 = 8.0 * (3.0 * a * a - 2.0 * rp * rp) * dxr[i] * dxr[j] * dxr[k] * dxr[l];

        result = consp * (term1 + term2 + term3 + term4 + term5 + term6);
    }
    return result;
}

/*------------------------------------------------------------------------------------
        Function: Eshelby tensor for strain field integration results (linear)
-------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_linear(int i, int j, int k, int l, int p, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;
    double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);

    double para = (1.0) / (8.0 * Pi * (1.0 - nu0));

    // the 0.5 is due to make the tensor major symmetry
    result = para * (Psi_1_4(x1, a, p, k, l, i, j)
            - 2.0 * nu0 * d[k][l] * Phi_1_2(x1, a, p, i, j)
            - (1.0 - nu0) * d[i][k] * Phi_1_2(x1, a, p, l, j)
            - (1.0 - nu0) * d[j][k] * Phi_1_2(x1, a, p, l, i)
            - (1.0 - nu0) * d[i][l] * Phi_1_2(x1, a, p, k, j)
            - (1.0 - nu0) * d[j][l] * Phi_1_2(x1, a, p, k, i));
    
    return result;
}

/*---------------------------------------------------------------------------------
        Function: Eshelby tensor for strain field integration results (quadratic)
------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_quadratic(int i, int j, int k, int l, int p, int q, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;
    double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);

    double para = (1.0) / (8.0 * Pi * (1.0 - nu0));

    result = para * (Psi_2_4(x1, a, p, q, k, l, i, j)
        - 2.0 * nu0 * d[k][l] * Phi_2_2(x1, a, p, q, i, j)
        - (1.0 - nu0) * d[i][k] * Phi_2_2(x1, a, p, q, l, j)
        - (1.0 - nu0) * d[j][k] * Phi_2_2(x1, a, p, q, l, i)
        - (1.0 - nu0) * d[i][l] * Phi_2_2(x1, a, p, q, k, j)
        - (1.0 - nu0) * d[j][l] * Phi_2_2(x1, a, p, q, k, i));

    return result;
}

// Eshelby tensor for der equivalent (uniform)
double Eshelby_strain_integration_der(int i, int j, int k, int l, int m, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;  double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double rp = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    if (a <= rp) {
        double para1 = (1.0 - 2.0 * nu0) / (8.0 * Pi * (1.0 - nu0));
        double para2 = (1.0) / (4.0 * Pi * (1.0 - nu0));
        double term1, term2;        // taking partial derivative w,r,t x_m
        
        term1 = d[i][l] * Phi_der2(x1, a, rp, k, j, m) + d[i][k] * Phi_der2(x1, a, rp, l, j, m) + d[j][l] * Phi_der2(x1, a, rp, k, i, m) + d[j][k] * Phi_der2(x1, a, rp, l, i, m) - d[k][l] * Phi_der2(x1, a, rp, i, j, m) - d[k][l] * Phi_der2(x1, a, rp, j, i, m);
        term2 = Psi_der2(x1, a, rp, i, k, l, j, m) + Psi_der2(x1, a, rp, j, k, l, i, m);
        result =-(para1 * term1 + para2 * term2);
    }
    else // interior has no derivatives since it is x not related
        result = 0.0;

    return result;
}

/*------------------------------------------------------------------------------------
        Function: first order derivative of Eshelby's tensor linear
--------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_linear_der(int i, int j, int k, int l, int p,int m, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;  double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    double para = 1.0 / (8.0 * Pi * (1.0 - nu0));


    result = para * (Psi_1_5(x1, a, p, k, l, i, j, m)
        - 2.0 * nu0 * d[k][l] * Phi_1_3(x1, a, p, i, j, m)
        - (1.0 - nu0) * d[i][k] * Phi_1_3(x1, a, p, l, j, m)
        - (1.0 - nu0) * d[j][k] * Phi_1_3(x1, a, p, l, i, m)
        - (1.0 - nu0) * d[i][l] * Phi_1_3(x1, a, p, k, j, m)
        - (1.0 - nu0) * d[j][l] * Phi_1_3(x1, a, p, k, i, m));

    return result;
}

/*------------------------------------------------------------------------------------
        Function: first order derivative of Eshelby's tensor linear
--------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_quadratic_der(int i, int j, int k, int l, int p, int q, int m, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;  double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    double para = 1.0 / (8.0 * Pi * (1.0 - nu0));


    result = para * (Psi_2_5(x1, a, p, q, k, l, i, j, m)
        - 2.0 * nu0 * d[k][l] * Phi_2_3(x1, a, p, q, i, j, m)
        - (1.0 - nu0) * d[i][k] * Phi_2_3(x1, a, p, q, l, j, m)
        - (1.0 - nu0) * d[j][k] * Phi_2_3(x1, a, p, q, l, i, m)
        - (1.0 - nu0) * d[i][l] * Phi_2_3(x1, a, p, q, k, j, m)
        - (1.0 - nu0) * d[j][l] * Phi_2_3(x1, a, p, q, k, i, m));

    return result;
}

/*------------------------------------------------------------------------------------
        Function: second order derivative of Eshelby's tensor uniform
--------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_uniform_der2(int i, int j, int k, int l, int m, int s, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;  double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    double para = 1.0 / (8.0 * Pi * (1.0 - nu0));


    result = para * (Psi_der3(x1, a, r, k, l, i, j, m, s)
        - 2.0 * nu0 * d[k][l] * Phi_der3(x1, a, r, i, j, m, s)
        - (1.0 - nu0) * d[i][k] * Phi_der3(x1, a, r, l, j, m, s)
        - (1.0 - nu0) * d[j][k] * Phi_der3(x1, a, r, l, i, m, s)
        - (1.0 - nu0) * d[i][l] * Phi_der3(x1, a, r, k, j, m, s)
        - (1.0 - nu0) * d[j][l] * Phi_der3(x1, a, r, k, i, m, s));

    return result;
}

/*------------------------------------------------------------------------------------
        Function: first order derivative of Eshelby's tensor linear
--------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_linear_der2(int i, int j, int k, int l, int p, int m,int s, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;  double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    double para = 1.0 / (8.0 * Pi * (1.0 - nu0));


    result = para * (Psi_1_6(x1, a, p, k, l, i, j, m, s)
        - 2.0 * nu0 * d[k][l] * Phi_1_4(x1, a, p, i, j, m, s)
        - (1.0 - nu0) * d[i][k] * Phi_1_4(x1, a, p, l, j, m, s)
        - (1.0 - nu0) * d[j][k] * Phi_1_4(x1, a, p, l, i, m, s)
        - (1.0 - nu0) * d[i][l] * Phi_1_4(x1, a, p, k, j, m, s)
        - (1.0 - nu0) * d[j][l] * Phi_1_4(x1, a, p, k, i, m, s));

    return result;
}

/*------------------------------------------------------------------------------------
        Function: first order derivative of Eshelby's tensor linear
--------------------------------------------------------------------------------------*/
double Eshelby_strain_integration_quadratic_der2(int i, int j, int k, int l, int p,int q, int m, int s, int fir, int sec, Ref<VectorXd> radius, Ref<MatrixXd> x_f, Ref<MatrixXd> x_s)
{
    double result = 0.0;  double x1[2];
    x1[0] = x_f(fir, 0) - x_s(sec, 0); x1[1] = x_f(fir, 1) - x_s(sec, 1);
    double a = radius[sec];
    double r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
    double para = 1.0 / (8.0 * Pi * (1.0 - nu0));


    result = para * (Psi_2_6(x1, a, p, q, k, l, i, j, m, s)
        - 2.0 * nu0 * d[k][l] * Phi_2_4(x1, a, p, q, i, j, m, s)
        - (1.0 - nu0) * d[i][k] * Phi_2_4(x1, a, p, q, l, j, m, s)
        - (1.0 - nu0) * d[j][k] * Phi_2_4(x1, a, p, q, l, i, m, s)
        - (1.0 - nu0) * d[i][l] * Phi_2_4(x1, a, p, q, k, j, m, s)
        - (1.0 - nu0) * d[j][l] * Phi_2_4(x1, a, p, q, k, i, m, s));

    return result;
}

/*------------------------------------------------------------------------------------
        Constitutive law with given parameters, nu and mu
-------------------------------------------------------------------------------------*/
double cons_law(double nu, double mu, int i, int j, int k, int l)
{
    double delta_ij, delta_ik, delta_il, delta_jk, delta_jl, delta_kl;
    delta_ij = delta_ik = delta_il = delta_jk = delta_jl = delta_kl = 0.0;

    double lambda = (2.0 * mu * nu) / (1.0 - 2.0 * nu);
    if (i == j) delta_ij = 1.0;
    if (i == k) delta_ik = 1.0;
    if (i == l) delta_il = 1.0;
    if (j == k) delta_jk = 1.0;
    if (j == l) delta_jl = 1.0;
    if (k == l) delta_kl = 1.0;

    double result = 0.0;
    result = lambda * delta_ij * delta_kl + mu * (delta_il * delta_jk + delta_ik * delta_jl);
    return result;
}