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
#include "Linear_Polygon_integral.h"
#include "Uniform_Polygon_integral.h"

double Linear_Polygon::Psi_p_b(double b, double lm, double lp)
{
    double result = 0.5 * b * (lm * lm - lp * lp - (b * b + lm * lm) * log(b * b + lm * lm) + (b * b + lp * lp) * log(b * b + lp * lp));
    return result;
}

double Linear_Polygon::Psi_p_lp(double b, double lp)
{
    double result = 0.5 * lp * (b * b + lp * lp) * (-1.0 + log(b * b + lp * lp));
    return result;
}

double Linear_Polygon::Psi_p_lm(double b, double lm)
{
    double result = -0.5 * lm * (b * b + lm * lm) * (-1.0 + log(b * b + lm * lm));
    return result;
}

double Linear_Polygon::PSI_p_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double reuslt = 0.0;
    return -Vnorm[j] * Psi_p_b(b, lm, lp) - Vdir[j] * (Psi_p_lp(b, lp) + Psi_p_lm(b, lm));
}

double Linear_Polygon::Psi_p_2b(double b, double lm, double lp)
{
    return 0.5 * (lm * lm - lp * lp - (3.0 * b * b + lm * lm) * log(b * b + lm * lm)
        + (3.0 * b * b + lp * lp) * log(b * b + lp * lp));
}

double Linear_Polygon::Psi_p_lp_b(double b, double lp)
{
    double result = b * lp * log(b * b + lp * lp);
    return result;
}

double Linear_Polygon::Psi_p_lm_b(double b, double lm)
{
    double result = -b * lm * log(b * b + lm * lm);
    return result;
}

double Linear_Polygon::Psi_p_2lp(double b, double lp)
{
    double result = 0.5 * (-b * b - lp * lp + (b * b + 3.0 * lp * lp) * log(b * b + lp * lp));
    return result;
}

double Linear_Polygon::Psi_p_2lm(double b, double lm)
{
    double result = 0.5 * (b * b + lm * lm - (b * b + 3.0 * lm * lm) * log(b * b + lm * lm));
    return result;
}

double Linear_Polygon::PSI_p_jk(int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = Vnorm[j] * Vnorm[k] * Psi_p_2b(b, lm, lp)
        + (Psi_p_lm_b(b, lm) + Psi_p_lp_b(b, lp)) * (Vnorm[j] * Vdir[k] + Vnorm[k] * Vdir[j])
        + (Psi_p_2lm(b, lm) + Psi_p_2lp(b, lp)) * Vdir[j] * Vdir[k];

    return result;
}
/*-------------------------------------------------------------------*/
// PSI_{p,jkl}

double Linear_Polygon::Psi_p_3b(double b, double lm, double lp)
{
    double result = 0.0;

    result = 2.0 * pow(b, 3.0) * (lm * lm - lp * lp) / ((b * b + lm * lm) * (b * b + lp * lp));
    result += -3.0 * b * log(b * b + lm * lm);
    result += 3.0 * b * log(b * b + lp * lp);
    if (b == 0) result = 0.0;
    return result;
}

double Linear_Polygon::Psi_p_lp_2b(double b, double lp)
{
    double result = 0.0;

    result = 2.0 * pow(b, 2.0) / (pow(b, 2.0) + pow(lp, 2.0)) + log(pow(b, 2.0) + pow(lp, 2.0));
    result *= lp;
    return result;
}

double Linear_Polygon::Psi_p_lm_2b(double b, double lm)
{
    double result = 0.0;

    result = 2.0 * pow(b, 2.0) / (pow(b, 2.0) + pow(lm, 2.0)) + log(pow(b, 2.0) + pow(lm, 2.0));
    result *= -lm;
    return result;
}

double Linear_Polygon::Psi_p_2lp_b(double b, double lp)
{
    double result = 0.0;

    result = 2.0 * pow(lp, 2.0) / (pow(b, 2.0) + pow(lp, 2.0)) + log(pow(b, 2.0) + pow(lp, 2.0));
    result *= b;
    return result;
}

double Linear_Polygon::Psi_p_2lm_b(double b, double lm)
{
    double result = 0.0;

    result = 2.0 * pow(lm, 2.0) / (pow(b, 2.0) + pow(lm, 2.0)) + log(pow(b, 2.0) + pow(lm, 2.0));
    result *= -b;
    return result;
}

double Linear_Polygon::Psi_p_3lp(double b, double lp)
{
    double result = 0.0;

    result = 2.0 * pow(lp, 3.0) / (pow(b, 2.0) + pow(lp, 2.0)) + 3.0 * lp * log(pow(b, 2.0) + pow(lp, 2.0));

    return result;
}

double Linear_Polygon::Psi_p_3lm(double b, double lm)
{
    double result = 0.0;

    result = -2.0 * pow(lm, 3.0) / (pow(b, 2.0) + pow(lm, 2.0)) - 3.0 * lm * log(pow(b, 2.0) + pow(lm, 2.0));

    return result;
}

double Linear_Polygon::PSI_p_jkl(int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = -Psi_p_3b(b, lm, lp) * Vnorm[j] * Vnorm[k] * Vnorm[l] - (Psi_p_lm_2b(b, lm) + Psi_p_lp_2b(b, lp)) * (Vnorm[j] * Vnorm[k] * Vdir[l] + Vnorm[j] * Vdir[k] * Vnorm[l] + Vdir[j] * Vnorm[k] * Vnorm[l]);
    result += -(Psi_p_2lp_b(b, lp) + Psi_p_2lm_b(b, lm)) * (Vnorm[j] * Vdir[k] * Vdir[l] + Vdir[j] * Vnorm[k] * Vdir[l] + Vdir[j] * Vdir[k] * Vnorm[l]);
    result += -(Psi_p_3lp(b, lp) + Psi_p_3lm(b, lm)) * Vdir[j] * Vdir[k] * Vdir[l];
    return result;
}


// harmonic potentials:
double Linear_Polygon::PHI_p(double b, double lm, double lp)
{
    double result = 0.5 * (lm * lm - lp * lp - (b * b + lm * lm) * log(b * b + lm * lm) + (b * b + lp * lp) * log(b * b + lp * lp));
    return result;
}

double Linear_Polygon::Phi_p_b(double b, double lm, double lp)
{
    double result = 0.0;
    result = b * (-log(pow(lm, 2.0) + pow(b, 2.0)) + log(pow(b, 2.0) + pow(lp, 2.0)));
    if (b == 0) result = 0.0;
    return result;
}

double Linear_Polygon::Phi_p_lp(double b, double lp)
{
    double result = 0.0;
    result = lp * log(pow(b, 2.0) + pow(lp, 2.0));
    return result;
}

double Linear_Polygon::Phi_p_lm(double b, double lm)
{
    double result = 0.0;
    result = -lm * log(pow(b, 2.0) + pow(lm, 2.0));
    return result;
}

double Linear_Polygon::PHI_p_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = -Vnorm[j] * Phi_p_b(b, lm, lp) - Vdir[j] * (Phi_p_lp(b, lp) + Phi_p_lm(b, lm));

    return result;
}

double Linear_Polygon::PSI_p_2(int p, int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
    double result = 0.0;
    double d[2][2] = { 1,0,0,1 };
    result = Vnorm[i] * Vdir[p] * PSI_p_j(j, b, lm, lp, Vnorm, Vdir);
    result += (d[i][p]) * Vnorm[j] * U_P.Psi(b, lm, lp);
    result += (d[j][p] - Vnorm[j] * Vnorm[p]) * Vnorm[i] * U_P.Psi(b, lm, lp);
    result += (x[p] + b * Vnorm[p]) * Vnorm[i] * U_P.PSI_j(j, b, lm, lp, Vnorm, Vdir);
    return result;

}

double Linear_Polygon::PSI_p_3(int p, int i, int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
    //Uniform_polygon Uni;

    double result = 0.0;
    double d[2][2] = { 1,0,0,1 };
    result = Vnorm[i] * Vdir[p] * PSI_p_jk(j, k, b, lm, lp, Vnorm, Vdir);
    result += (d[i][p]) * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir) * Vnorm[j];
    result += (d[j][p] - Vnorm[j] * Vnorm[p]) * U_P.PSI_j(k, b, lm, lp, Vnorm, Vdir) * Vnorm[i];
    result += (d[k][p] - Vnorm[k] * Vnorm[p]) * U_P.PSI_j(j, b, lm, lp, Vnorm, Vdir) * Vnorm[i];
    result += (x[p] + b * Vnorm[p]) * U_P.PSI_jk(j, k, b, lm, lp, Vnorm, Vdir) * Vnorm[i];
    return result;
}


double Linear_Polygon::PSI_p_4(int p, int i, int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
    //Uniform_polygon Uni;

    double result = 0.0;
    double d[2][2] = { 1,0,0,1 };

    result = Vnorm[i] * Vdir[p] * PSI_p_jkl(j, k, l, b, lm, lp, Vnorm, Vdir);
    result += (d[i][p]) * Vnorm[j] * U_P.PSI_jk(k, l, b, lm, lp, Vnorm, Vdir);
    result += (d[p][j] - Vnorm[j] * Vnorm[p]) * Vnorm[i] * U_P.PSI_jk(k, l, b, lm, lp, Vnorm, Vdir);
    result += (d[p][k] - Vnorm[k] * Vnorm[p]) * Vnorm[i] * U_P.PSI_jk(j, l, b, lm, lp, Vnorm, Vdir);
    result += (d[p][l] - Vnorm[l] * Vnorm[p]) * Vnorm[i] * U_P.PSI_jk(j, k, b, lm, lp, Vnorm, Vdir);
    result += (x[p] + b * Vnorm[p]) * Vnorm[i] * U_P.PSI_jkl(j, k, l, b, lm, lp, Vnorm, Vdir);

    return result;


}

double Linear_Polygon::PHI_p_o(int p, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
    double result = x[p] * U_P.PHI_o(b, lm, lp);

    result += b * b / 9.0 * (
        -8 * lm + 8.0 * lp + 6.0 * b * atan(lm / b) - 6.0 * b * atan(lp / b)
        + 3.0 * lm * log(b * b + lm * lm) - 3.0 * lp * log(b * b + lp * lp)
        ) * Vnorm[p];

    result += 1.0 / 18.0 * b * (
        -5.0 * lm * lm + 5.0 * lp * lp + 3.0 * (b * b + lm * lm) * log(b * b + lm * lm)
        - 3.0 * (b * b + lp * lp) * log(b * b + lp * lp)
        ) * Vdir[p];

    return result;

}

double Linear_Polygon::PHI_p_1(int p, int i, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
    //Uniform_polygon Uni;
    double result = 0.0;
    double d[2][2] = { 1,0,0,1 };
    result = Vdir[p] * Vnorm[i] * PHI_p(b, lm, lp);
    result += (d[i][p]) * U_P.PHI_o(b, lm, lp);
    result += (x[p] + b * Vnorm[p]) * Vnorm[i] * U_P.PHI(b, lm, lp);

    return result;
}

double Linear_Polygon::PHI_p_2(int p, int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir, double* x)
{
    //Uniform_polygon Uni;
    double result = 0.0;
    double d[2][2] = { 1,0,0,1 };
    result = Vdir[p] * Vnorm[i] * PHI_p_j(j, b, lm, lp, Vnorm, Vdir);
    result += (d[i][p]) * Vnorm[j] * U_P.PHI(b, lm, lp);
    result += (d[p][j] - Vnorm[p] * Vnorm[j]) * Vnorm[i] * U_P.PHI(b, lm, lp);
    result += (x[p] + b * Vnorm[p]) * Vnorm[i] * U_P.PHI_j(j, b, lm, lp, Vnorm, Vdir);
    return result;


}