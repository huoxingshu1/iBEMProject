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
#include "Uniform_Polyhedral_integral.h"

double Uniform_Polyhedral::Psi_I1(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / 6.0) * (b * le * dis - 2.0 * pow(abs(a), 3.0) * atan(le / b)
			+ b * (3 * a * a + b * b) * atanh(le / dis)
			+ 2.0 * pow(a, 3.0) * atan(a * le / (b * dis))

			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = (1.0 / 6.0) * (b * le * dis + b * b * b * atanh(le / dis));
		}
		else {
			result = 0.0;
		}
	}
	else {
		result = 0.0;
	}

	return result;
}

// derivatives of Psi_I1, used for Eshelby's tensor
//double Uniform_Polyhedral::Psi_I1_
double Uniform_Polyhedral::Psi_I1_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = a * (
			-abs(a) * atan(le / b) + a * atan(a * le / (b * dis)) + b * atanh(le / dis)
			);
	}
	else {
		result = 0.0;
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (6.0 * (b * b + le * le))) * (
			le * (-3.0 * a * a + dis * dis) * dis + 2.0 * le * pow(abs(a), 3.0) + 3.0 * (a * a + b * b) * (b * b + le * le) * atanh(le / dis)
			);

	}
	else if (a == 0) {
		if (b != 0) {
			result = (1.0 / 6.0) * (le * dis + 3.0 * b * b * atanh(le / dis));
		}
		else {
			result = le * abs(le) / 6.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 1.0 / (6.0 * le) * (
				(-2.0 * a * a + le * le) * dis + 2.0 * pow(abs(a), 3.0) + 3.0 * a * a * le * atanh(le / dis)
				);
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * (
			dis * dis * dis - abs(a) * abs(a) * abs(a)
			) / (3.0 * (b * b + le * le));
	}
	else if (a == 0) {
		result = b * dis / 3.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (dis * dis * dis - abs(a) * abs(a) * abs(a)) / (3.0 * b);
		}
		else {
			result = 0;
		}
	}


	return result;
}

// second partial derivative of Psi_I1
 /*double Uniform_Polyhedral::Psi_I1_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);

	if (a < 0 && b != 0 && le != 0) {
		result = 2.0 * a * (atan(le / b) + atan(a * le / (b * dis))) + b * atanh(le / dis);
	}

	else if (a > 0 && b != 0 && le != 0) {
		result = 2.0 * a * (-atan(le / b) + atan(a * le / (b * dis))) + b * atanh(le / dis);
	}

	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * atanh(le / dis);
		}
		else if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
		else {
			result = 0.0;
		}
	}

	else if (a != 0) {
		result = 0;
	}

	return result;
}*/

/*double Uniform_Polyhedral::Psi_I1_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a < 0 && b != 0 && le != 0) {
		result = a * (le / (a - dis) + atanh(le / dis));
	}

	else if (a > 0 && b != 0 && le != 0) {
		result = a * (le * (a - dis) / (b * b + le * le) + atanh(le / dis));
	}

	else if (a == 0) {
		result = 0;
	}

	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			if (a < 0) {
				result = a * (le / (a - dis) + atanh(le / dis));
			}
			else if (a > 0) {
				result = a * ((a - dis) / le + atanh(le / dis));
			}
		}
	}


	return result;
}*/

double Uniform_Polyhedral::Psi_I1_a_a(double a, double b, double le)
{
	double result; double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = -2.0 * abs(a) * atan(le / b) + 2.0 * a * atan(a * le / (b * dis)) + b * atanh(le / dis);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * atanh(le / dis);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		result = 0.0;
	}

	return result;
}


double Uniform_Polyhedral::Psi_I1_a_b(double a, double b, double le)
{
	double result; double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = a * (le * (abs(a) - dis) / (be * be) + atanh(le / dis));
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = a / le * (-dis + abs(a) + le * atanh(le / dis));
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * (dis - abs(a)) / (b * b + le * le);
	}
	else if (a == 0) {
		result = 0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (dis - abs(a)) / b;
		}
		else {  // as b ==0, results are zero{
			result = 0;
		}
	}
	return result;
}

double Uniform_Polyhedral::Psi_I1_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = (b) / (3.0 * (b * b + le * le) * (b * b + le * le)) * (
			-le * (-3.0 * a * a + dis * dis) * dis - 2.0 * le * pow(abs(a), 3.0) + 3.0 * pow(b * b + le * le, 2.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = -b * le / (3.0 * dis) + b * atanh(le / dis);
		}
		else {
			result = 0;
		}
	}
	else if (a != 0) {
		result = 0;
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_b_le(double a, double b, double le)
{
	double result = 0.0;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (3.0 * pow(b * b + le * le, 2.0))) * (
			dis * (2.0 * pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + pow(le, 4.0) + a * a * (-b * b + le * le))
			+ (b * b - le * le) * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = (2.0 * b * b + le * le) / (3.0 * dis);
		}
		else if (b == 0 && le != 0) {
			result = abs(le) / 3.0;
		}
		else if (b != 0 && le == 0) {
			result = 2.0 * abs(b) / 3.0;
		}
		else {
			result = 0.0;
		}
	}

	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (-(a * a - 2.0 * b * b) * dis + pow(abs(a), 3.0)) / (3.0 * b * b);
		}
		else if (b == 0 && le != 0) {
			result = (dis * dis * dis - pow(abs(a), 3.0)) / (3.0 * le * le);
		}
		else if (b == 0 && le == 0) {
			result = 0.5 * abs(a);
		}

	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_le_le(double a, double b, double le)
{
	double result = 0.0;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b * le) / (3.0 * dis * pow(b * b + le * le, 2.0)) * (
			-2.0 * pow(a, 4.0) - a * a * (b * b + le * le) + pow(b * b + le * le, 2.0) + 2.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * le / (3.0 * dis);
		}
		else {
			result = 0;
		}
	}
	else if (a != 0) {
		result = 0;
	}

	return result;
}


// third derivatives 
double Uniform_Polyhedral::Psi_I1_a_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a < 0 && b != 0 && le != 0) {
		result = a * b * le / ((a * a + b * b) * dis) + 2.0 * atan(le / b) + 2.0 * atan(a * le / (b * dis));
	}
	else if (a > 0 && b != 0 && le != 0) {
		result = a * b * le / ((a * a + b * b) * dis) - 2.0 * atan(le / b) + 2.0 * atan(a * le / (b * dis));
	}
	else {
		result = 0;
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_a_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a < 0 && b != 0 && le != 0) {
		result = -le * (2.0 * a * a + b * b + a * b * b / dis) / ((a * a + b * b) * (dis - a)) + atanh(le / dis);
	}
	else if (a > 0 && b != 0 && le != 0) {
		result = 2.0 * a * le / (b * b + le * le) + (le / dis) * (-2.0 + b * b / (a * a + b * b) - 2.0 * a * a / (b * b + le * le)) + atanh(le / dis);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -le / dis + atanh(le / dis);
		}
		else {   // include b = 0, le != 0, this is infinite
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = 2.0 * (abs(a) - dis) / le + atanh(le / dis);
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}
	return result;
}

double Uniform_Polyhedral::Psi_I1_a_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b / (b * b + le * le)) * (-2.0 * abs(a) + dis + (a * a) / dis);
	}
	else if (a == 0) {
		if (b != 0) {
			result = b / dis;
		}
		else {
			result = 0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = ((2.0 * a * a + b * b) / dis - 2.0 * abs(a)) / b;
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}
	return result;
}

double Uniform_Polyhedral::Psi_I1_a_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (a * b * le / ((a * a + b * b) * pow(b * b + le * le, 2.0) * dis)) * (
			2.0 * pow(a, 4.0) - le * le * (b * b + le * le) + a * a * (3.0 * b * b + le * le) - 2.0 * (a * a + b * b) * dis * abs(a)
			);
	}
	else {
		result = 0;
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_a_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (a / (pow(b * b + le * le, 2.0) * dis)) * (
			-pow(a * b, 2.0) + (a * a + b * b) * le * le + pow(le, 4.0) + (b * b - le * le) * dis * abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (-a * a / dis + abs(a)) / (b * b);
		}
		else if (b == 0 && le != 0) {
			result = a * (dis - abs(a)) / (le * le);
		}
		else if (b == 0 && le == 0) {
			result = 0.5 * a / abs(a);
		}
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_a_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = -a * b * le * (
			a * a + dis * dis - 2.0 * dis * abs(a)
			) / (pow(b * b + le * le, 2.0) * dis);
	}
	else {
		result = 0;
	}
	return result;

}

double Uniform_Polyhedral::Psi_I1_b_b_b(double a, double b, double le)
{
	double result; double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (3.0 * (a * a + b * b) * pow(b * b + le * le, 3.0) * dis)) * (
			2.0 * pow(a, 6.0) * le * (-3.0 * b * b + le * le) - b * b * le * pow(b * b + le * le, 2.0) * (3.0 * b * b + 4.0 * le * le)
			- a * a * le * (b * b + le * le) * (3.0 * pow(b, 4.0) + pow(le, 4.0)) + pow(a, 4.0) * (-9.0 * pow(b, 4.0) * le + pow(le, 5.0))
			)
			- 2.0 * le * (-3.0 * b * b + le * le) * pow(abs(a), 3.0) / (3.0 * pow(b * b + le * le, 3.0))
			+ atanh(le / dis);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -le * (3.0 * b * b + 4.0 * le * le) / (3.0 * dis * dis * dis) + atanh(le / dis);
		}
		else { // infinite
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = (1.0 / (3.0 * pow(le, 3.0) * dis)) * (2.0 * pow(a, 4.0) + a * a * le * le - pow(le, 4.0)) - 2.0 * pow(abs(a), 3.0) / (3.0 * le * le * le) + atanh(le / dis);
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_b_b_le(double a, double b, double le)
{
	double result; double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b / (3.0 * pow(b * b + le * le, 3.0) * dis)) * (
			2.0 * pow(a, 4.0) * (b * b - 3.0 * le * le) + a * a * (b * b - 3.0 * le * le) * (b * b + le * le)
			+ pow(b * b + le * le, 2.0) * (2.0 * b * b + 3.0 * le * le) - 2.0 * (b * b - 3.0 * le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = b * (2.0 * b * b + 3.0 * le * le) / (3.0 * dis * dis * dis);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (2.0 * pow(a, 4.0) + pow(a * b, 2.0) + 2.0 * pow(b, 4.0) - 2.0 * dis * pow(abs(a), 3.0)) / (3.0 * pow(b, 3.0) * dis);
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}

	return result;
}

double Uniform_Polyhedral::Psi_I1_b_le_le(double a, double b, double le)
{
	double result; double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (le / (3.0 * pow(b * b + le * le, 3.0) * dis)) * (
			pow(a, 4.0) * (6.0 * b * b - 2.0 * le * le) + a * a * (3.0 * b * b - le * le) * (b * b + le * le)
			+ le * le * pow(b * b + le * le, 2.0) + 2.0 * (-3.0 * b * b + le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = pow(le, 3.0) / (3.0 * dis * dis * dis);
		}
		else if (le != 0) {
			result = 0.0; // le / (3.0 * abs(le));
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = (-2.0 * pow(a, 4.0) - pow(a * le, 2.0) + pow(le, 4.0) + 2.0 * dis * pow(abs(a), 3.0)) / (3.0 * dis * pow(le, 3.0));
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}

	}
	return result;
}

double Uniform_Polyhedral::Psi_I1_le_le_le(double a, double b, double le)
{
	double result; double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b / (3.0 * pow(b * b + le * le, 3.0) * dis)) * (
			-2.0 * pow(a, 4.0) * (b * b - 3.0 * le * le) - a * a * (b * b - 3.0 * le * le) * (b * b + le * le)
			+ b * b * pow(b * b + le * le, 2.0) + 2.0 * (b * b - 3.0 * le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = pow(b, 3.0) / (3.0 * dis * dis * dis);
		}
		else {
			result = 0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (-2.0 * pow(a, 4.0) - pow(a * b, 2.0) + pow(b, 4.0) + 2.0 * dis * pow(abs(a), 3.0)) / (3.0 * dis * pow(b, 3.0));
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}

	return result;
}

////////////////////////////////////////////////////
// the original Phi potential without the Green;s theorem
double Uniform_Polyhedral::Phi_u(double a, double b, double le)
{
	double result = 0.0; //double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * a * (
			a * abs(b) * asin(a * le / sqrt((a * a + b * b) * (b * b + le * le))) / b
			- abs(a) * atan(le / b) + b * atanh(le / sqrt(a * a + b * b + le * le))
			);
	}
	else if (a == 0) {		// this include a == b == le == 0
		result = 0;
	}
	else if (a != 0) {
		result = 0.0;
	}

	return result;
}


double Uniform_Polyhedral::Phi_u_II(double a, double b, double le)
{
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = (a * b / abs(b)) * asin(a * le / sqrt((a * a + b * b) * (b * b + le * le))) - abs(a) * atan(le / b) + b * atanh(le / sqrt(a * a + b * b + le * le));
	}
	else if (a == 0) {
		if (b != 0) {	// le may equal zero, all includes
			result = b * atanh(le / sqrt(b * b + le * le));
		}
		else {  // b ==0, le != 0
			result = 0.0;
		}
	}

	else if (a != 0) {
		result = 0.0;
	}

	else { // a, b, le equals zero
		result = 0.0;
	}

	return result;
}

double Uniform_Polyhedral::Phi_u_II_a(double a, double b, double le)
{
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = -a / abs(a) * atan(le / b) + asin(a * le * abs(b) / (b * sqrt((a * a + b * b) * (b * b + le * le))));
	}
	else if (a == 0) {
		result = 0.0;
	}

	else if (a != 0) {

		result = 0.0;
	}

	else { // a, b, le equals zero
		result = 0.0;
	}
	return result;
}

double Uniform_Polyhedral::Phi_u_II_b(double a, double b, double le)
{
	double result; double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (b * b + le * le)) * (
			-le * dis + le * abs(a) + (b * b + le * le) * atanh(le / dis));
	}
	else if (a == 0) {
		if (b != 0) {	// le may equal zero, all includes
			result = -le / dis + atanh(le / dis);
		}
		else { // b ==0, le != 0
			result = 0.0;
		}
	}

	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.0;
		}
		else if (b == 0 && le == 0) {
			result = 0.0;
		}
		else if (b == 0 && le != 0) {
			result = (-dis + abs(a)) / le + atanh(le / dis);
		}
	}

	else { // a, b, le equals zero
		result = 0.0;
	}
	return result;
}

double Uniform_Polyhedral::Phi_u_II_le(double a, double b, double le)
{
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = b * (a * a + b * b + le * le - sqrt(a * a + b * b + le * le) * abs(a)) / ((b * b + le * le) * sqrt(a * a + b * b + le * le));
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b / sqrt(b * b + le * le);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (sqrt(a * a + b * b) - abs(a)) / b;
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}

	return result;
}

// define the combined potentials

double Uniform_Polyhedral::I1(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Phi_u_II(a, b, lp) - Phi_u_II(a, b, lm));
}

double Uniform_Polyhedral::I1_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j] * (Phi_u_II_a(a, b, lp) - Phi_u_II_a(a, b, lm))
		- Vnorm[j] * (Phi_u_II_b(a, b, lp) - Phi_u_II_b(a, b, lm))
		- Vdir[j] * (Phi_u_II_le(a, b, lp) - Phi_u_II_le(a, b, lm))
		);
	return result;
}

double Uniform_Polyhedral::I4(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Psi_I1(a, b, lp) - Psi_I1(a, b, lm));
}

double Uniform_Polyhedral::I4_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = ((-Svnorm[j]) * (Psi_I1_a(a, b, lp) - Psi_I1_a(a, b, lm))
		+ (-Vnorm[j]) * (Psi_I1_b(a, b, lp) - Psi_I1_b(a, b, lm))
		+ (-Vdir[j]) * (Psi_I1_le(a, b, lp) - Psi_I1_le(a, b, lm))
		);
	return result;
}

double Uniform_Polyhedral::I4_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I1_a_a(a, b, lp) - Psi_I1_a_a(a, b, lm)) - Vnorm[k] * (Psi_I1_a_b(a, b, lp) - Psi_I1_a_b(a, b, lm)) - Vdir[k] * (Psi_I1_a_le(a, b, lp) - Psi_I1_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I1_a_b(a, b, lp) - Psi_I1_a_b(a, b, lm)) - Vnorm[k] * (Psi_I1_b_b(a, b, lp) - Psi_I1_b_b(a, b, lm)) - Vdir[k] * (Psi_I1_b_le(a, b, lp) - Psi_I1_b_le(a, b, lm)))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I1_a_le(a, b, lp) - Psi_I1_a_le(a, b, lm)) - Vnorm[k] * (Psi_I1_b_le(a, b, lp) - Psi_I1_b_le(a, b, lm)) - Vdir[k] * (Psi_I1_le_le(a, b, lp) - Psi_I1_le_le(a, b, lm)))
		);

	return result;
}

double Uniform_Polyhedral::I4_jkl(int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result;
	result = (-Svnorm[j]) * (-Svnorm[k] * (-Svnorm[l] * (Psi_I1_a_a_a(a, b, lp) - Psi_I1_a_a_a(a, b, lm)) - Vnorm[l] * (Psi_I1_a_a_b(a, b, lp) - Psi_I1_a_a_b(a, b, lm)) - Vdir[l] * (Psi_I1_a_a_le(a, b, lp) - Psi_I1_a_a_le(a, b, lm)))
		- Vnorm[k] * (-Svnorm[l] * (Psi_I1_a_a_b(a, b, lp) - Psi_I1_a_a_b(a, b, lm)) - Vnorm[l] * (Psi_I1_a_b_b(a, b, lp) - Psi_I1_a_b_b(a, b, lm)) - Vdir[l] * (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)))
		- Vdir[k] * (-Svnorm[l] * (Psi_I1_a_a_le(a, b, lp) - Psi_I1_a_a_le(a, b, lm)) - Vnorm[l] * (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)) - Vdir[l] * (Psi_I1_a_le_le(a, b, lp) - Psi_I1_a_le_le(a, b, lm)))) // end of the first

		+ (-Vnorm[j]) * (-Svnorm[k] * (-Svnorm[l] * (Psi_I1_a_a_b(a, b, lp) - Psi_I1_a_a_b(a, b, lm)) - Vnorm[l] * (Psi_I1_a_b_b(a, b, lp) - Psi_I1_a_b_b(a, b, lm)) - Vdir[l] * (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)))
			- Vnorm[k] * (-Svnorm[l] * (Psi_I1_a_b_b(a, b, lp) - Psi_I1_a_b_b(a, b, lm)) - Vnorm[l] * (Psi_I1_b_b_b(a, b, lp) - Psi_I1_b_b_b(a, b, lm)) - Vdir[l] * (Psi_I1_b_b_le(a, b, lp) - Psi_I1_b_b_le(a, b, lm)))
			- Vdir[k] * (-Svnorm[l] * (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)) - Vnorm[l] * (Psi_I1_b_b_le(a, b, lp) - Psi_I1_b_b_le(a, b, lm)) - Vdir[l] * (Psi_I1_b_le_le(a, b, lp) - Psi_I1_b_le_le(a, b, lm))))  // end of the second

		+ (-Vdir[j]) * (-Svnorm[k] * (-Svnorm[l] * (Psi_I1_a_a_le(a, b, lp) - Psi_I1_a_a_le(a, b, lm)) - Vnorm[l] * (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)) - Vdir[l] * (Psi_I1_a_le_le(a, b, lp) - Psi_I1_a_le_le(a, b, lm)))
			- Vnorm[k] * (-Svnorm[l] * (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)) - Vnorm[l] * (Psi_I1_b_b_le(a, b, lp) - Psi_I1_b_b_le(a, b, lm)) - Vdir[l] * (Psi_I1_b_le_le(a, b, lp) - Psi_I1_b_le_le(a, b, lm)))
			- Vdir[k] * (-Svnorm[l] * (Psi_I1_a_le_le(a, b, lp) - Psi_I1_a_le_le(a, b, lm)) - Vnorm[l] * (Psi_I1_b_le_le(a, b, lp) - Psi_I1_b_le_le(a, b, lm)) - Vdir[l] * (Psi_I1_le_le_le(a, b, lp) - Psi_I1_le_le_le(a, b, lm))));

	return result;

}

// the complete set of potentials

double Uniform_Polyhedral::PHI(double a, double b, double lm, double lp)
{
	return Phi_u(a, b, lp) - Phi_u(a, b, lm);
}

/* double Uniform_Polyhedral::PHI(double a, double b, double lm, double lp)
 {
	 double result = 0.0;
	 if (b != 0)
		 result = 0.5 * a * (-abs(a) * (atan(lp / b) - atan(lm / b))
			 + b * (log(lp + sqrt(a * a + b * b + lp * lp)) - log(lm + sqrt(a * a + b * b + lm * lm)))
			 + 1.0 * a * (atan((a * lp) / (b * sqrt(a * a + b * b + lp * lp))) - atan((a * lm) / (b * sqrt(a * a + b * b + lm * lm))))
			 );
	 else {
		 result = 0.0;
	 }
	 return result;
 }*/


double Uniform_Polyhedral::PHI_i(int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = -(Svnorm[i]) * (Phi_u_II(a, b, lp) - Phi_u_II(a, b, lm));
	return result;

}

double Uniform_Polyhedral::PHI_ij(int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = -(Svnorm[i]) * (-(Phi_u_II_a(a, b, lp) - Phi_u_II_a(a, b, lm)) * Svnorm[j]
		- (Phi_u_II_b(a, b, lp) - Phi_u_II_b(a, b, lm)) * Vnorm[j]
		- (Phi_u_II_le(a, b, lp)) * Vdir[j] + (Phi_u_II_le(a, b, lm)) * Vdir[j]
		);

	return result;

}

double Uniform_Polyhedral::PSI_i(int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = -(Svnorm[i]) * (Psi_I1(a, b, lp) - Psi_I1(a, b, lm));
	return result;
}

double Uniform_Polyhedral::PSI_ij(int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = -(Svnorm[i]) * ((-Svnorm[j]) * (Psi_I1_a(a, b, lp) - Psi_I1_a(a, b, lm))
		+ (-Vnorm[j]) * (Psi_I1_b(a, b, lp) - Psi_I1_b(a, b, lm))
		+ (-Vdir[j]) * (Psi_I1_le(a, b, lp) - Psi_I1_le(a, b, lm))
		);
	return result;
}

double Uniform_Polyhedral::PSI_ijk(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = -(Svnorm[i]) * (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I1_a_a(a, b, lp) - Psi_I1_a_a(a, b, lm)) - Vnorm[k] * (Psi_I1_a_b(a, b, lp) - Psi_I1_a_b(a, b, lm)) - Vdir[k] * (Psi_I1_a_le(a, b, lp) - Psi_I1_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I1_a_b(a, b, lp) - Psi_I1_a_b(a, b, lm)) - Vnorm[k] * (Psi_I1_b_b(a, b, lp) - Psi_I1_b_b(a, b, lm)) - Vdir[k] * ((Psi_I1_b_le(a, b, lp) - Psi_I1_b_le(a, b, lm))))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I1_a_le(a, b, lp) - Psi_I1_a_le(a, b, lm)) - Vnorm[k] * ((Psi_I1_b_le(a, b, lp) - Psi_I1_b_le(a, b, lm))) - Vdir[k] * ((Psi_I1_le_le(a, b, lp) - Psi_I1_le_le(a, b, lm))))
		);
	return result;
}

double Uniform_Polyhedral::PSI_ijkl(int l, int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	// revise it in a concise form
	double result = (-(Psi_I1_a_a_a(a, b, lp) - Psi_I1_a_a_a(a, b, lm)) * Svnorm[i] * Svnorm[j] * Svnorm[k] - (Psi_I1_b_b_b(a, b, lp) - Psi_I1_b_b_b(a, b, lm)) * Vnorm[i] * Vnorm[j] * Vnorm[k]
		- (Psi_I1_le_le_le(a, b, lp) - Psi_I1_le_le_le(a, b, lm)) * Vdir[i] * Vdir[j] * Vdir[k] - (Psi_I1_a_a_b(a, b, lp) - Psi_I1_a_a_b(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vnorm[k] + Svnorm[i] * Vnorm[j] * Svnorm[k] + Vnorm[i] * Svnorm[j] * Svnorm[k])
		- (Psi_I1_a_a_le(a, b, lp) - Psi_I1_a_a_le(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vdir[k] + Svnorm[i] * Vdir[j] * Svnorm[k] + Vdir[i] * Svnorm[j] * Svnorm[k]) - (Psi_I1_a_b_b(a, b, lp) - Psi_I1_a_b_b(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vnorm[k] + Svnorm[j] * Vnorm[i] * Vnorm[k] + Svnorm[k] * Vnorm[j] * Vnorm[i])
		- (Psi_I1_b_b_le(a, b, lp) - Psi_I1_b_b_le(a, b, lm)) * (Vdir[i] * Vnorm[j] * Vnorm[k] + Vdir[j] * Vnorm[i] * Vnorm[k] + Vdir[k] * Vnorm[i] * Vnorm[j]) - (Psi_I1_a_le_le(a, b, lp) - Psi_I1_a_le_le(a, b, lm)) * (Svnorm[i] * Vdir[j] * Vdir[k] + Svnorm[j] * Vdir[i] * Vdir[k] + Svnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I1_b_le_le(a, b, lp) - Psi_I1_b_le_le(a, b, lm)) * (Vnorm[i] * Vdir[j] * Vdir[k] + Vnorm[j] * Vdir[i] * Vdir[k] + Vnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I1_a_b_le(a, b, lp) - Psi_I1_a_b_le(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vdir[k] + Svnorm[i] * Vnorm[k] * Vdir[j] + Svnorm[j] * Vnorm[i] * Vdir[k] + Svnorm[k] * Vnorm[i] * Vdir[j] + Svnorm[j] * Vnorm[k] * Vdir[i] + Svnorm[k] * Vnorm[j] * Vdir[i])) * (-Svnorm[l]);


	return result;
}