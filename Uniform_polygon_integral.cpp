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

#include "Uniform_Polygon_integral.h"
#include <cmath>

// define the functions
double Uniform_Polygon::PHI_o(double b, double lm, double lp)
{
    double result = 0.0;
    result = (-1.5 * lm + 1.5 * lp + b * atan(lm / b) - b * atan(lp / b) + 0.5 * lm * log(b * b + lm * lm) - 0.5 * lp * log(b * b + lp * lp));
    result *= b;
    if (b == 0) result = 0;
    return result;
}


double Uniform_Polygon::PHI(double b, double lm, double lp)
{
    double result = 0.0;
    if (b != 0)
        result = 2.0 * b * (atan(lp / b) - atan(lm / b)) + lp * log(b * b + lp * lp) - lm * log(b * b + lm * lm) - 2.0 * (lp - lm);
    else if (b == 0)
        result = lp * log(b * b + lp * lp) - lm * log(b * b + lm * lm) - 2.0 * (lp - lm);
    return result;
}

double Uniform_Polygon::Phi_b(double b, double lm, double lp)
{
    double result = 0.0;
    if (b != 0) {
        result = 2.0 * (atan(lp / b) - atan(lm / b));
    }
    else
        result = 0.0;

    return result;
}

double Uniform_Polygon::Phi_lp(double b, double lp)
{
    return log(b * b + lp * lp);
}

double Uniform_Polygon::Phi_lm(double b, double lm)
{
    return -log(b * b + lm * lm);
}

double Uniform_Polygon::PHI_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    return -Vnorm[j] * Phi_b(b, lm, lp) - Vdir[j] * (Phi_lp(b, lp) + Phi_lm(b, lm));
}

// Biharmonic potentials
// no derivative form, in the Green;s theorem
double Uniform_Polygon::Psi(double b, double lm, double lp)
{
    double result = (1.0 / 18.0) * (
        21.0 * b * b * (lm - lp) + 5.0 * (lm * lm * lm - lp * lp * lp)
        + 12.0 * b * b * b * (-atan(lm / b) + atan(lp / b))
        - 3.0 * lm * (3.0 * b * b + lm * lm) * log(b * b + lm * lm)
        + 3.0 * lp * (3.0 * b * b + lp * lp) * log(b * b + lp * lp)
        );
    if (b == 0) {
        result = (1.0 / 18.0) * (5.0 * lm * lm * lm - 5.0 * lp * lp * lp
            - 3.0 * lm * lm * lm * log(lm * lm) + 3.0 * lp * lp * lp * log(lp * lp));
    }
    return result;
}


// first derivatives:
double Uniform_Polygon::Psi_b(double b, double lm, double lp)
{
    double result = (b / 3.0) * (-7.0 * lp + (lp * b * b) / (b * b + lp * lp) + pow(lp, 3.0) / (b * b + lp * lp)
        + 7.0 * lm - (lm * b * b) / (b * b + lm * lm) - pow(lm, 3.0) / (b * b + lm * lm) +
        6.0 * b * atan(lp / b) - 6.0 * b * atan(lm / b) + 3.0 * lp * log(b * b + lp * lp) - 3.0 * lm * log(b * b + lm * lm)
        );
    if (b == 0) result = 0.0;
    return result;
}

double Uniform_Polygon::Psi_lp(double b, double lp)
{
    double p1 = 0.0;
    p1 = (-0.5 * pow(b, 4.0) - pow(b, 2.0) * pow(lp, 2.0) - 0.5 * pow(lp, 4.0) + (0.5 * pow(b, 4.0) + pow(b, 2.0) * pow(lp, 2.0) + 0.5 * pow(lp, 4.0)) * log(b * b + lp * lp)) / (b * b + lp * lp);
    return p1;
}

double Uniform_Polygon::Psi_lm(double b, double lp)
{
    double result = -(-0.5 * pow(b, 4.0) - pow(b, 2.0) * pow(lp, 2.0) - 0.5 * pow(lp, 4.0) + (0.5 * pow(b, 4.0) + pow(b, 2.0) * pow(lp, 2.0) + 0.5 * pow(lp, 4.0)) * log(b * b + lp * lp)) / (b * b + lp * lp);
    return result;
}


// second derivatives
double Uniform_Polygon::Psi_2b(double b, double lm, double lp)
{
    double result = 0.0;

    result = -7.0 * lp - (6.0 * pow(b, 4.0) * lp / ((b * b + lp * lp) * (b * b + lp * lp)));
    result += -6.0 * pow(b, 2.0) * pow(lp, 3.0) / ((b * b + lp * lp) * (b * b + lp * lp));
    result += 7.0 * pow(b, 2.0) * lp / (b * b + lp * lp);
    result += pow(lp, 3.0) / (b * b + lp * lp);
    result += 7.0 * lm + (6.0 * pow(b, 4.0) * lm / ((b * b + lm * lm) * (b * b + lm * lm)));
    result += 6.0 * pow(b, 2.0) * pow(lm, 3.0) / ((b * b + lm * lm) * (b * b + lm * lm));
    result += -7.0 * pow(b, 2.0) * lm / (b * b + lm * lm);
    result += -pow(lm, 3.0) / (b * b + lm * lm);
    result += 12.0 * b * atan(lp / b) - 12.0 * b * atan(lm / b);
    result += 3.0 * lp * log(b * b + lp * lp) - 3.0 * lm * log(b * b + lm * lm);
    result *= (1.0 / 3.0);

    if (b == 0) result = 2.0 * lm - 2.0 * lp - lm * log(lm * lm) + lp * log(lp * lp);

    return result;
}

double Uniform_Polygon::Psi_lp_b(double b, double lp)
{
    double result = 0.0;
    result = (pow(b, 5.0) + 2.0 * pow(b, 3.0) * pow(lp, 2.0) + b * pow(lp, 4.0)) / ((b * b + lp * lp) * (b * b + lp * lp));
    result *= log(b * b + lp * lp);
    if (b == 0) result = 0.0;
    return result;
}

double Uniform_Polygon::Psi_lm_b(double b, double lm)
{
    double result = 0.0;
    result = (pow(b, 5.0) + 2.0 * pow(b, 3.0) * pow(lm, 2.0) + b * pow(lm, 4.0)) / ((b * b + lm * lm) * (b * b + lm * lm));
    result *= -log(b * b + lm * lm);
    if (b == 0) result = 0.0;
    return result;

}

double Uniform_Polygon::Psi_2lp(double b, double lp)
{
    double result = lp * log(b * b + lp * lp);
    return result;
}

double Uniform_Polygon::Psi_2lm(double b, double lm)
{
    double result = -lm * log(b * b + lm * lm);
    return result;
}


double Uniform_Polygon::PSI_jk(int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = Vnorm[j] * Vnorm[k] * Psi_2b(b, lm, lp)
        + (Psi_lm_b(b, lm) + Psi_lp_b(b, lp)) * (Vnorm[j] * Vdir[k] + Vnorm[k] * Vdir[j])
        + (Psi_2lm(b, lm) + Psi_2lp(b, lp)) * Vdir[j] * Vdir[k];

    return result;
}

// third derivatives
double Uniform_Polygon::Psi_3b(double b, double lm, double lp)
{
    double result = 0.0;
    result = 8.0 * pow(b, 5.0) * lp / pow((pow(b, 2.0) + pow(lp, 2.0)), 3.0) +
        (8.0 / 3.0) * pow(b, 3.0) * pow(lp, 3.0) / pow((pow(b, 2.0) + pow(lp, 2.0)), 3.0) -
        (16.0 / 3.0) * b * pow(lp, 5.0) / pow((pow(b, 2.0) + pow(lp, 2.0)), 3.0) -
        18.0 * pow(b, 3.0) * lp / pow((pow(b, 2.0) + pow(lp, 2.0)), 2.0) -
        (14.0 / 3.0) * b * pow(lp, 3.0) / pow((pow(b, 2.0) + pow(lp, 2.0)), 2.0) +
        8.0 * b * lp / (pow(b, 2.0) + pow(lp, 2.0));

    result += -8.0 * pow(b, 5.0) * lm / pow((pow(b, 2.0) + pow(lm, 2.0)), 3.0) -
        (8.0 / 3.0) * pow(b, 3.0) * pow(lm, 3.0) / pow((pow(b, 2.0) + pow(lm, 2.0)), 3.0) +
        (16.0 / 3.0) * b * pow(lm, 5.0) / pow((pow(b, 2.0) + pow(lm, 2.0)), 3.0) +
        18.0 * pow(b, 3.0) * lm / pow((pow(b, 2.0) + pow(lm, 2.0)), 2.0) +
        (14.0 / 3.0) * b * pow(lm, 3.0) / pow((pow(b, 2.0) + pow(lm, 2.0)), 2.0) -
        8.0 * b * lm / (pow(b, 2.0) + pow(lm, 2.0));

    result += 4.0 * (atan(lp / b) - atan(lm / b));
    if (b == 0) {
        if (lp * lm > 0 || lp * lm == 0) {
            result = 0.0;
        }
        else if (lp * lm < 0) {
            result = 0;// 4.0 * pi;
           // cout << "b = 0" << endl;
        }

    }

    return result;
}

double Uniform_Polygon::Psi_3lp(double b, double lp)
{
    double result = 0.0;
    result = 2.0 * pow(b, 4.0) * pow(lp, 2.0) + 4.0 * pow(b, 2.0) * pow(lp, 4.0) +
        2.0 * pow(lp, 6.0) +
        (pow(b, 6.0) + 3.0 * pow(b, 4.0) * pow(lp, 2.0) + 3.0 * pow(b, 2.0) * pow(lp, 4.0) + pow(lp, 6.0)) * log(pow(b, 2.0) + pow(lp, 2.0));

    result *= (1.0 / pow((pow(b, 2.0) + pow(lp, 2.0)), 3.0));
    if (b == 0) result = 2.0 + log(lp * lp);
    return result;
}

double Uniform_Polygon::Psi_3lm(double b, double lm)
{
    double result = 0.0;
    result = 2.0 * pow(b, 4.0) * pow(lm, 2.0) + 4.0 * pow(b, 2.0) * pow(lm, 4.0) +
        2.0 * pow(lm, 6.0) +
        (pow(b, 6.0) + 3.0 * pow(b, 4.0) * pow(lm, 2.0) + 3.0 * pow(b, 2.0) * pow(lm, 4.0) + pow(lm, 6.0)) * log(pow(b, 2.0) + pow(lm, 2.0));

    result *= -(1.0 / pow((pow(b, 2.0) + pow(lm, 2.0)), 3.0));

    if (b == 0) result = -2.0 - log(lm * lm);

    return result;
}

double Uniform_Polygon::Psi_2lp_b(double b, double lp)
{
    double result = 0.0;

    result = 2.0 * pow(b, 5.0) * lp + 4.0 * pow(b, 3.0) * pow(lp, 3.0) + 2.0 * b * pow(lp, 5.0);
    result *= (1.0 / pow((pow(b, 2.0) + pow(lp, 2.0)), 3.0));
    if (b == 0) result = 0;
    return result;
}

double Uniform_Polygon::Psi_2lm_b(double b, double lm)
{
    double result = 0.0;

    result = 2.0 * pow(b, 5.0) * lm + 4.0 * pow(b, 3.0) * pow(lm, 3.0) + 2.0 * b * pow(lm, 5.0);
    result *= -(1.0 / pow((pow(b, 2.0) + pow(lm, 2.0)), 3.0));
    if (b == 0) result = 0;
    return result;
}

double Uniform_Polygon::Psi_lp_2b(double b, double lp)
{
    double result = 0.0;

    result = 2.0 * pow(b, 6.0) + 4.0 * pow(b, 4.0) * pow(lp, 2.0) + 2.0 * pow(b, 2.0) * pow(lp, 4.0) +
        (pow(b, 6.0) + 3.0 * pow(b, 4.0) * pow(lp, 2.0) + 3.0 * pow(b, 2.0) * pow(lp, 4.0) + pow(lp, 6.0)) * log(pow(b, 2.0) + pow(lp, 2.0));

    result *= (1.0 / pow((pow(b, 2.0) + pow(lp, 2.0)), 3.0));
    //if (b == 0) result = 0;
    return result;
}

double Uniform_Polygon::Psi_lm_2b(double b, double lm)
{
    double result = 0.0;

    result = 2.0 * pow(b, 6.0) + 4.0 * pow(b, 4.0) * pow(lm, 2.0) + 2.0 * pow(b, 2.0) * pow(lm, 4.0) +
        (pow(b, 6.0) + 3.0 * pow(b, 4.0) * pow(lm, 2.0) + 3.0 * pow(b, 2.0) * pow(lm, 4.0) + pow(lm, 6.0)) * log(pow(b, 2.0) + pow(lm, 2.0));

    result *= -(1.0 / pow((pow(b, 2.0) + pow(lm, 2.0)), 3.0));
    //if (b == 0) result = 0;
    return result;
}


double Uniform_Polygon::PSI_j(int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = -Vnorm[j] * Psi_b(b, lm, lp) - Vdir[j] * (Psi_lp(b, lp) + Psi_lm(b, lm));
    return result;
}


double Uniform_Polygon::PSI_jkl(int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = -Psi_3b(b, lm, lp) * Vnorm[j] * Vnorm[k] * Vnorm[l] - (Psi_lm_2b(b, lm) + Psi_lp_2b(b, lp)) * (Vnorm[j] * Vnorm[k] * Vdir[l] + Vnorm[j] * Vdir[k] * Vnorm[l] + Vdir[j] * Vnorm[k] * Vnorm[l]);
    result += -(Psi_2lp_b(b, lp) + Psi_2lm_b(b, lm)) * (Vnorm[j] * Vdir[k] * Vdir[l] + Vdir[j] * Vnorm[k] * Vdir[l] + Vdir[j] * Vdir[k] * Vnorm[l]);
    result += -(Psi_3lp(b, lp) + Psi_3lm(b, lm)) * Vdir[j] * Vdir[k] * Vdir[l];

    return result;
}


/*-------------------------------------
   Combined potentials
-------------------------------------*/
double Uniform_Polygon::PHI_1(int i, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    return Vnorm[i] * PHI(b, lm, lp);
}


double Uniform_Polygon::PHI_2(int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    return Vnorm[i] * (-Vnorm[j] * Phi_b(b, lm, lp) - Vdir[j] * (Phi_lp(b, lp) + Phi_lm(b, lm)));
}

double Uniform_Polygon::PSI_1(int i, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    return Vnorm[i] * Psi(b, lm, lp);
}

double Uniform_Polygon::PSI_2(int i, int j, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    return Vnorm[i] * PSI_j(j, b, lm, lp, Vnorm, Vdir);
}

double Uniform_Polygon::PSI_3(int i, int j, int k, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = Vnorm[i] * (Vnorm[j] * Vnorm[k] * Psi_2b(b, lm, lp)
        + (Psi_lm_b(b, lm) + Psi_lp_b(b, lp)) * (Vnorm[j] * Vdir[k] + Vnorm[k] * Vdir[j])
        + (Psi_2lm(b, lm) + Psi_2lp(b, lp)) * Vdir[j] * Vdir[k]);

    return result;
}

double Uniform_Polygon::PSI_4(int i, int j, int k, int l, double b, double lm, double lp, double* Vnorm, double* Vdir)
{
    double result = 0.0;

    result = -Psi_3b(b, lm, lp) * Vnorm[j] * Vnorm[k] * Vnorm[l] - (Psi_lm_2b(b, lm) + Psi_lp_2b(b, lp)) * (Vnorm[j] * Vnorm[k] * Vdir[l] + Vnorm[j] * Vdir[k] * Vnorm[l] + Vdir[j] * Vnorm[k] * Vnorm[l]);
    result += -(Psi_2lp_b(b, lp) + Psi_2lm_b(b, lm)) * (Vnorm[j] * Vdir[k] * Vdir[l] + Vdir[j] * Vnorm[k] * Vdir[l] + Vdir[j] * Vdir[k] * Vnorm[l]);
    result += -(Psi_3lp(b, lp) + Psi_3lm(b, lm)) * Vdir[j] * Vdir[k] * Vdir[l];

    return Vnorm[i] * result;
}
