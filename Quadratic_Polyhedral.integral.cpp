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

#include "Quadratic_Polyhedral_integral.h"
#include "GlobalConstant.h"
#include <cmath>

// second derivatives of Psi_I4
double Quadratic_Polyhedral::Psi_I4_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (3.0 * be * be)) * (
			-4.0 * a * a * b * le * dis + 4.0 * pow(abs(a), 3.0) * (b * le + be * be * atan(le / b))
			+ be * be * (-4.0 * pow(a, 3.0) * atan(a * le / (b * dis)) + pow(b, 3.0) * atanh(le / dis))
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * b * b * atanh(le / be) / 3.0;
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

double Quadratic_Polyhedral::Psi_I4_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (a * b * b) / (3.0 * pow(be, 4.0)) * (
			-le * (-3.0 * a * a + dis * dis) * dis - 2.0 * le * pow(abs(a), 3.0) + 3.0 * pow(be, 4.0) * atanh(le / dis)
			);
	}
	else {
		result = 0.0;
	}


	return result;
}

double Quadratic_Polyhedral::Psi_I4_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * pow(b, 3.0) / (3.0 * pow(be, 4.0) * dis) * (
			-2.0 * pow(a, 4.0) - a * a * be * be + pow(be, 4.0) + 2.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (-2.0 * pow(a, 4.0) - a * a * b * b + pow(b, 4.0) + 2.0 * dis * pow(abs(a), 3.0)) / (3.0 * b * dis);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I4_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / (15.0 * pow(be, 6.0)) * (
			le * dis * (4.0 * pow(a, 4.0) * (-b * b + le * le) + 3.0 * pow(be, 4.0) * (2.0 * b * b + 3.0 * le * le) + 2.0 * a * a * (pow(b, 4.0) - pow(le, 4.0)))
			+ 4.0 * pow(a, 4.0) * le * (b * b - le * le) * abs(a) + 15.0 * (a * a + 2.0 * b * b) * pow(be, 6.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = 0.2 / (be * be) * (
				b * le * be * (2.0 * b * b + 3.0 * le * le) + 10.0 * b * b * b * be * be * atanh(le / be)
				);
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

double Quadratic_Polyhedral::Psi_I4_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b / (15.0 * pow(be, 6.0) * dis) * (
			2.0 * pow(a, 6.0) * (b * b - 3.0 * le * le) + 3.0 * pow(be, 6.0) * (4.0 * b * b + 3.0 * le * le)
			+ a * a * pow(be, 4.0) * (11.0 * b * b + 12.0 * le * le) + pow(a, 4.0) * (b * b - 3.0 * le * le) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * b * (3.0 * be * be + b * b) / (5.0 * be);
		}
		else if (b != 0 && le == 0) {
			result = 0.8 * b * b * abs(b);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 2.0 * pow(a, 6.0) + pow(a * a * b, 2.0) + 11.0 * pow(a * b * b, 2.0)
				+ 12.0 * pow(b, 6.0) - 2.0 * pow(a, 4.0) * dis * abs(a);
			result = result / (15.0 * b * b * dis);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I4_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b * b * le / (15.0 * pow(be, 6.0) * dis) * (
			8.0 * pow(a, 6.0) - a * a * pow(be, 4.0) + 3.0 * pow(be, 6.0) + 4.0 * pow(a, 4.0) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * b * b * le / (5.0 * be);
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



// third derivatives of Psi_I4
double Quadratic_Polyhedral::Psi_I4_a_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * (
			-b * le / ((a * a + b * b) * be * be * dis) * (
				4.0 * pow(a, 4.0) + 3.0 * pow(b * be, 2.0) + 4.0 * a * a * (b * b + be * be)
				)
			- 4.0 * a * atan(a * le / (b * dis)) + 4.0 * a * (b * le / (be * be) + atan(le / b)) * a / abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		result = 0.0;
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I4_a_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b / (3.0 * pow(be, 4.0)) * (
			le / ((a * a + b * b) * dis) * (
				8.0 * pow(a, 6.0) - 4.0 * pow(a * le * be, 2.0) - b * b * pow(be, 4.0) + 4.0 * pow(a, 4.0) * (3.0 * b * b + le * le)
				)
			- 8.0 * le * pow(abs(a), 3.0) + 3.0 * pow(be, 4.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -b * b * le / (3.0 * be) + b * b * atanh(le / be);
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

double Quadratic_Polyhedral::Psi_I4_a_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = pow(b, 3.0) / (3.0 * pow(be, 4.0) * dis) * (
			-8.0 * pow(a, 4.0) - 4.0 * a * a * be * be + pow(be, 4.0) + 8.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * b / (3.0 * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (-8.0 * pow(a, 4.0) - 4.0 * a * a * b * b + pow(b, 4.0) + 8.0 * dis * pow(abs(a), 3.0)) / (3.0 * b * dis);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I4_a_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / 3.0) * a * b * (
			le / ((a * a + b * b) * pow(be, 6.0) * dis) * (
				4.0 * pow(a, 6.0) * (-b * b + le * le) - b * b * pow(be, 4.0) * (4.0 * b * b + 5.0 * le * le)
				+ 2.0 * pow(a, 4.0) * (-3.0 * pow(b, 4.0) + 2.0 * pow(b * le, 2.0) + pow(le, 4.0))
				- a * a * be * be * (3.0 * pow(b, 4.0) + pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
				)
			+ 4.0 * le * (b * b - le * le) * pow(abs(a), 3.0) / pow(be, 6.0)
			+ 6.0 * atanh(le / dis)
			);
	}
	else {
		result = 0.0;
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I4_a_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * b / (3.0 * pow(be, 6.0) * dis) * (
			2.0 * pow(a, 4.0) * (b * b - 3.0 * le * le) + a * a * be * be * (b * b - 3.0 * le * le)
			+ pow(be, 4.0) * (2.0 * be * be + le * le) - 2.0 * (b * b - 3.0 * le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a / (3.0 * b * b * dis) * (2.0 * pow(a, 4.0) + a * a * b * b + 2.0 * pow(b, 4.0) - 2.0 * dis * pow(abs(a), 3.0));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I4_a_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = -a * pow(b, 3.0) * le / (3.0 * pow(be, 6.0) * dis) * (
			-8.0 * pow(a, 4.0) - 4.0 * a * a * be * be + pow(be, 4.0) + 8.0 * dis * pow(abs(a), 3.0)
			);
	}
	else {
		result = 0.0;
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I4_b_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = -4.0 * pow(a, 4.0) * abs(a) * le / (15.0 * pow(be, 8.0)) * (3.0 * pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + pow(le, 4.0))
			+ (a * a + b * b * 6.0) * atanh(le / dis)
			+ le / (15.0 * (a * a + b * b) * pow(be, 8.0) * dis) * (
				-3.0 * b * b * pow(be, 6.0) * (6.0 * pow(b, 4.0) + 4.0 * pow(b * le, 2.0) - 3.0 * pow(le, 4.0))
				+ 4.0 * pow(a, 8.0) * (3.0 * pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + pow(le, 4.0))
				+ 2.0 * pow(a, 6.0) * (9.0 * pow(b, 6.0) - 21.0 * pow(b * b * le, 2.0) - 5.0 * pow(b * le * le, 2.0) + pow(le, 6.0))
				+ pow(a, 4.0) * (12.0 * pow(b, 8.0) + 21.0 * pow(b * b * b * le, 2.0) + 37.0 * pow(b * le, 4.0) + 35.0 * pow(b * le * le * le, 2.0) + 7.0 * pow(le, 8.0))
				+ a * a * pow(be, 4.0) * (3.0 * pow(b, 6.0) + 19.0 * pow(b * b * le, 2.0) + 19.0 * pow(b * le * le, 2.0) + 9.0 * pow(le, 6.0))
				);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -le / (5.0 * be * be * be) * (6.0 * pow(b, 4.0) + 4.0 * b * b * le * le - 3.0 * pow(le, 4.0)) + 6.0 * b * b * atanh(le / be);
		}
		else {
			result = 0.6 * le * abs(le);
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.0;
		}
		else if (b == 0 && le != 0) {
			result = a * a * atanh(le / dis) - 4.0 / (15.0 * le * le * le) * pow(a, 4.0) * abs(a) + (4.0 * pow(a, 6.0) + 2.0 * a * a * a * a * le * le + 7.0 * a * a * le * le * le * le + 9.0 * pow(le, 6.0)) / (15.0 * le * le * le * dis);
		}
		else if (b == 0 && le == 0) {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I4_b_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / (15.0 * pow(be, 8.0) * dis) * (
			-4.0 * pow(a, 6.0) * (pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0))
			+ 3.0 * pow(be, 6.0) * (12.0 * pow(b, 4.0) + 19.0 * pow(b * le, 2.0) + 6.0 * pow(le, 4.0))
			+ a * a * pow(be, 4.0) * (23.0 * pow(b, 4.0) + 41 * pow(b * le, 2.0) + 24.0 * pow(le, 4.0))
			- 2.0 * pow(a, 4.0) * (pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0)) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * (12.0 * b * b * b * b + 19.0 * b * b * le * le + 6.0 * le * le * le * le) / (5.0 * be * be * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (-4.0 * pow(a, 6.0) - 2.0 * a * a * a * a * b * b + 23.0 * a * a * b * b * b * b + 36.0 * pow(b, 6.0) + 4.0 * a * a * a * a * dis * abs(a)) / (15.0 * b * b * b * dis);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I4_b_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b * le / (5.0 * pow(be, 8.0) * dis) * (
			8.0 * pow(a, 6.0) * (-b * b + le * le) + a * a * (b * b - le * le) * pow(be, 4.0)
			+ pow(be, 6.0) * (2.0 * be * be + le * le) - 4.0 * pow(a, 4.0) * (b * b - le * le) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * le * (2.0 * b * b + 3.0 * le * le) / (5.0 * be * be * be);
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

double Quadratic_Polyhedral::Psi_I4_le_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = pow(b, 3.0) / (15.0 * pow(be, 8.0) * dis) * (
			8.0 * pow(a, 6.0) * (b * b - 5.0 * le * le) - a * a * (b * b - 5.0 * le * le) * pow(be, 4.0)
			+ 3.0 * b * b * pow(be, 6.0) + 4.0 * pow(a, 4.0) * (b * b - 5.0 * le * le) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = pow(b, 5.0) / (5.0 * be * be * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (8.0 * pow(a, 6.0) + 4.0 * a * a * a * a * b * b - a * a * b * b * b * b + 3.0 * pow(b, 6.0) - 8.0 * a * a * a * a * dis * abs(a)) / (15.0 * b * b * b * dis);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

// second derivatives of Psi_I5
double Quadratic_Polyhedral::Psi_I5_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 1.0 / (30.0 * be * be) * (
			10.0 * b * b * dis * (3.0 * a * a + dis * dis) - (b * b - le * le) * 20.0 * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = b * b * be / 3.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 2.0 * pow(abs(a), 3.0) / 3.0;
		}
		else if (b != 0 && le == 0) {
			result = (dis * (4.0 * a * a + b * b) - 2.0 * pow(abs(a), 3.0)) / 3.0;
		}
		else if (b == 0 && le == 0) {
			result = 2.0 * pow(abs(a), 3.0) / 3.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I5_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b / (3.0 * pow(be, 4.0)) * (
			dis * (3.0 * pow(b, 4.0) + 5.0 * pow(b * le, 2.0) + 2.0 * le * le * (a * a + le * le)) - 2.0 * le * le * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * b * dis;
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I5_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * b * le / (3.0 * pow(be, 4.0)) * (
			(-3.0 * a * a + dis * dis) * dis + 2.0 * pow(abs(a), 3.0)
			);
	}
	else {
		result = 0.0;
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I5_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 1.0 / (15.0 * pow(be, 6.0)) * (
			dis * (2.0 * pow(a, 4.0) * le * le * (-3.0 * b * b + le * le) + pow(be, 4.0) * (20.0 * pow(b, 4.0) + 19.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
				+ a * a * be * be * (5.0 * pow(b, 4.0) + 13.0 * pow(b * le, 2.0) + 4.0 * pow(le, 4.0)))
			- 2.0 * pow(a, 4.0) * le * le * (-3.0 * b * b + le * le) * abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = (20.0 * b * b * b * b + 19.0 * b * b * le * le + 2.0 * le * le * le * le) / (15.0 * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = dis * (a * a + 4.0 * b * b) / 3.0;
		}
		else if (b == 0 && le != 0) {
			result = 2.0 / (15.0 * le * le) * (
				pow(dis, 5.0) - pow(abs(a), 5.0)
				);
		}
		else if (b == 0 && le == 0) {
			result = abs(a) * a * a / 3.0;
		}

	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * le / (15.0 * pow(be, 6.0)) * (
			dis * (4.0 * pow(a, 4.0) * (b * b - le * le) + 3.0 * pow(be, 4.0) * (3.0 * b * b + 2.0 * le * le) + 2.0 * a * a * (-pow(b, 4.0) + pow(le, 4.0)))
			+ 4.0 * pow(a, 4.0) * (-b * b + le * le) * abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * le * (2.0 * be * be + b * b) * 0.2 / be;
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

double Quadratic_Polyhedral::Psi_I5_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b / (15.0 * pow(be, 6.0)) * (
			dis * (-2.0 * pow(a, 4.0) * (b * b - 3.0 * le * le) + a * a * (b * b - 3.0 * le * le) * be * be + 3.0 * pow(be, 4.0) * (be * be + le * le))
			+ 2.0 * pow(a, 4.0) * (b * b - 3.0 * le * le) * abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * (be * be + le * le) * 0.2 / be;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 1.0 / (15.0 * b * b) * (dis * (-2.0 * pow(a, 4.0) + a * a * b * b + 3.0 * pow(b, 4.0)) + 2.0 * pow(abs(a), 5.0));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}


// third derivatives of Psi_I5
double Quadratic_Polyhedral::Psi_I5_a_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a / (be * be) * (
			b * b * (4.0 * a * a + 3.0 * be * be) / dis + 2.0 * a * (-b * b + le * le) * a / abs(a)
			);
	}
	else if (a == 0) {
		result = 0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (4.0 * a * a + 3.0 * b * b - 2.0 * a * dis * a / abs(a)) / dis;
		}
		else if (b == 0 && le != 0) {
			result = 2.0 * a * a * a / abs(a);
		}
		else if (b == 0 && le == 0) {
			result = 2.0 * a * a * a / abs(a);
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_a_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / (3.0 * pow(be, 4.0) * dis) * (
			3.0 * pow(b, 4.0) * (2.0 * a * a + b * b) + 8.0 * pow(a * a + b * b, 2.0) * le * le + (10.0 * a * a + 7.0 * b * b) * pow(le, 4.0)
			+ 2.0 * pow(le, 6.0) - 8.0 * le * le * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * (3.0 * b * b + 2.0 * le * le) / (3.0 * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = b * (2.0 * a * a + b * b) / dis;
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_a_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b * le / (3.0 * pow(be, 4.0) * dis) * (
			-8.0 * pow(a, 4.0) - 4.0 * a * a * be * be
			+ pow(be, 4.0) + 8.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * le / (3.0 * be);
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

double Quadratic_Polyhedral::Psi_I5_a_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a / (3.0 * pow(be, 6.0) * dis) * (
			2.0 * pow(a * a * le, 2.0) * (-3.0 * b * b + le * le) + pow(be, 4.0) * (6.0 * pow(b, 4.0) + 9.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
			+ a * a * be * be * (3.0 * pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + 4.0 * pow(le, 4.0))
			- 2.0 * le * le * (-3.0 * b * b + le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (a * a + 2.0 * b * b) / dis;
		}
		else if (b == 0 && le != 0) {
			result = 2.0 * a * dis * dis * dis / (3.0 * le * le) - 2.0 * a * pow(abs(a), 3.0) / (3.0 * le * le);
		}
		else if (b == 0 && le == 0) {
			result = a * abs(a);
		}

	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_a_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * le / (3.0 * pow(be, 6.0) * dis) * (
			4.0 * pow(a, 4.0) * (b * b - le * le) + pow(be, 4.0) * (be * be + le * le)
			+ 2.0 * a * a * (pow(b, 4.0) - pow(le, 4.0)) + 4.0 * (le * le - b * b) * dis * pow(abs(a), 3.0)
			);
	}
	else {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_a_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * b / (3.0 * pow(be, 6.0) * dis) * (
			-2.0 * pow(a, 4.0) * (b * b - 3.0 * le * le) - a * a * (b * b - 3.0 * le * le) * (b * b + le * le)
			+ b * b * pow(be, 4.0) + 2.0 * (b * b - 3.0 * le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (-2.0 * pow(a, 4.0) - a * a * b * b + pow(b, 4.0) + 2.0 * dis * pow(abs(a), 3.0)) / (3.0 * b * b * dis);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_b_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / (5.0 * pow(be, 8.0) * dis) * (
			8.0 * pow(a, 6.0) * (b * b - le * le) * le * le + pow(be, 6.0) * (20.0 * pow(b, 4.0) + 33.0 * pow(b * le, 2.0) + 12.0 * pow(le, 4.0))
			+ a * a * pow(be, 4.0) * (15.0 * pow(b, 4.0) + 29.0 * pow(b * le, 2.0) + 16.0 * pow(le, 4.0))
			+ 4.0 * pow(a, 4.0) * le * le * (b * b - le * le) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * (20.0 * b * b * b * b + 33.0 * b * b * le * le + 12.0 * le * le * le * le) / (5.0 * be * be * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = b * (3.0 * a * a + 4.0 * b * b) / dis;
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I5_b_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = le / (15.0 * pow(be, 8.0) * dis) * (
			-4.0 * pow(a, 6.0) * (3.0 * pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + pow(le, 4.0))
			+ 3.0 * pow(be, 6.0) * (6.0 * pow(b, 4.0) + 9.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
			+ a * a * pow(be, 4.0) * (9.0 * pow(b, 4.0) + 11.0 * pow(b * le, 2.0) + 8.0 * pow(le, 4.0))
			- 2.0 * pow(a, 4.0) * (3.0 * pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + pow(le, 4.0)) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = le / (5.0 * be * be * be) * (6.0 * b * b * b * b + 9.0 * b * b * le * le + 2.0 * le * le * le * le);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 4.0 / (15.0 * le * le * le) * pow(a, 4.0) * abs(a) - 2.0 * (2.0 * a * a - 3.0 * le * le) * dis * dis * dis / (15.0 * le * le * le);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I5_b_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / (15.0 * pow(be, 8.0) * dis) * (
			4.0 * pow(a, 6.0) * (pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0))
			+ 3.0 * pow(be, 6.0) * (3.0 * pow(b, 4.0) + 6.0 * pow(b * le, 2.0) + 4.0 * pow(le, 4.0))
			+ a * a * pow(be, 4.0) * (7.0 * pow(b, 4.0) + 19.0 * pow(b * le, 2.0) + 6.0 * pow(le, 4.0))
			+ 2.0 * pow(a, 4.0) * (pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0)) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b / (5.0 * be * be * be) * (3.0 * pow(b, 4.0) + 6.0 * b * b * le * le + 4.0 * pow(le, 4.0));
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (4.0 * pow(a, 6.0) + 2.0 * a * a * a * a * b * b + 7.0 * a * a * b * b * b * b + 9.0 * pow(b, 6.0) - 4.0 * a * a * a * a * dis * abs(a)) / (15.0 * b * b * b * dis);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I5_le_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * b * le / (5.0 * pow(be, 8.0) * dis) * (
			8.0 * pow(a, 6.0) * (b * b - le * le) - a * a * (b * b - le * le) * pow(be, 4.0)
			+ pow(be, 6.0) * (2.0 * be * be + b * b) + 4.0 * pow(a, 4.0) * (b * b - le * le) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * le * (3.0 * b * b + 2.0 * le * le) / (5.0 * be * be * be);
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

// second derivatives of Psi_pq_I6
double Quadratic_Polyhedral::Psi_I6_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 1.0 / (6.0 * be * be) * (
			b * le * dis * (7.0 * a * a + dis * dis) + 8.0 * pow(abs(a), 3.0) * (-b * le + be * be * atan(le / b))
			- be * be * (8.0 * pow(a, 3.0) * atan(a * le / (b * dis)) + b * (9.0 * a * a + b * b) * atanh(le / dis))
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = (b * le * be - b * b * b * atanh(le / be)) / 6.0;
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

double Quadratic_Polyhedral::Psi_I6_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a / (6.0 * pow(be, 4.0)) * (
			le * dis * (3.0 * pow(b, 4.0) + 4.0 * (a * a + b * b) * le * le + pow(le, 4.0))
			- 4.0 * pow(abs(a) * le, 3.0)
			- 3.0 * (a * a + b * b) * pow(be, 4.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = a / (6.0 * pow(le, 4.0)) * (
				le * dis * (4.0 * a * a * le * le + pow(le, 4.0)) - 4.0 * le * le * le * pow(abs(a), 3.0)
				- 3.0 * a * a * pow(le, 4.0) * atanh(le / dis)
				);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Quadratic_Polyhedral::Psi_I6_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * le * le / (3.0 * pow(be, 4.0) * dis) * (
			-2.0 * pow(a, 4.0) - a * a * be * be
			+ pow(be, 4.0) + 2.0 * dis * pow(abs(a), 3.0)
			);
	}
	else {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / 30.0 * (
			le * dis / (pow(be, 6.0)) * (
				-16.0 * pow(a * a * le, 2.0) + be * be * (15.0 * pow(b, 4.0) + 8.0 * (a * a + 3.0 * b * b) * le * le + 9.0 * pow(le, 4.0))
				)
			+ 16.0 * pow(a * le, 3.0) * a * abs(a) / pow(be, 6.0)
			- 15.0 * (a * a + b * b) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = (5.0 * b * b * b * le + 3.0 * b * le * le * le - 5.0 * b * b * b * be * atanh(le / be)) / (10.0 * be);
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

double Quadratic_Polyhedral::Psi_I6_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = le * le / (15.0 * pow(be, 6.0)) * (
			dis * (
				pow(a, 4.0) * (6.0 * b * b - 2.0 * le * le) + 3.0 * pow(be, 4.0) * (2.0 * b * b + le * le)
				+ a * a * (-3.0 * pow(b, 4.0) - 2.0 * pow(b * le, 2.0) + pow(le, 4.0))
				)
			+ 2.0 * pow(a, 4.0) * (-3.0 * b * b + le * le) * abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = le * le * (2.0 * b * b + le * le) * 0.2 / be;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = dis * (-2.0 * a * a * a * a * le * le + a * a * le * le * le * le + 3.0 * pow(le, 6.0)) + 2.0 * pow(a, 4.0) * le * le * abs(a);
			result = result / (15.0 * le * le * le * le);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}


double Quadratic_Polyhedral::Psi_I6_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * le / (15.0 * pow(be, 6.0)) * (
			dis * (
				4.0 * pow(a, 4.0) * (-b * b + le * le) + 3.0 * pow(be, 4.0) * (2.0 * b * b + 3.0 * le * le)
				+ 2.0 * a * a * (pow(b, 4.0) - pow(le, 4.0)))
			+ 4.0 * pow(a, 4.0) * (b * b - le * le) * abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * le * (2.0 * be * be + le * le) * 0.2 / be;
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

// third derivatives of Psi_pq_I6
double Quadratic_Polyhedral::Psi_I6_a_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a / (be * be) * (
			b * le * (4.0 * a * a + 3.0 * be * be) / dis + 4.0 * abs(a) * (-b * le + be * be * atan(le / b))
			- be * be * (4.0 * a * atan(a * le / (b * dis)) + 3.0 * b * atanh(le / dis))
			);
	}
	else {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_a_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 1.0 / (6.0 * pow(be, 4.0) * dis) * (
			16.0 * pow(a, 4.0) * pow(le, 3.0) + le * pow(be, 4.0) * (3.0 * b * b + le * le)
			+ a * a * le * be * be * (9.0 * be * be + 8.0 * le * le)
			- dis * (
				16.0 * pow(abs(a) * le, 3.0) + 3.0 * (3.0 * a * a + b * b) * pow(be, 4.0) * atanh(le / dis)
				)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = (3.0 * b * b * le + le * le * le - 3.0 * b * b * dis * atanh(le / dis)) / (6.0 * dis);
		}
		else {
			result = le * abs(le) / 6.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = (dis * (16.0 * a * a + le * le) - 16.0 * pow(abs(a), 3.0) - 9.0 * a * a * le * atanh(le / dis)) / (6.0 * le);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_a_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * le * le / (3.0 * pow(be, 4.0) * dis) * (
			-8.0 * pow(a, 4.0) - 4.0 * a * a * be * be + pow(be, 4.0)
			+ 8.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * le * le / (3.0 * be);
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

double Quadratic_Polyhedral::Psi_I6_a_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = -a * b / (3.0 * pow(be, 6.0) * dis) * (
			-3.0 * pow(b, 6.0) * le + 2.0 * (4.0 * pow(a, 4.0) + 2.0 * pow(a * b, 2.0) - 5.0 * pow(b, 4.0)) * pow(le, 3.0)
			+ (4.0 * a * a - 11.0 * b * b) * pow(le, 5.0) - 4.0 * pow(le, 7.0)
			+ dis * (-8.0 * pow(le * abs(a), 3.0) + 3.0 * pow(be, 6.0) * atanh(le / dis))
			);
	}
	else {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_a_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * le * le / (3.0 * pow(be, 6.0) * dis) * (
			pow(a, 4.0) * (6.0 * b * b - 2.0 * le * le) + a * a * (3.0 * b * b - le * le) * (b * b + le * le)
			+ le * le * pow(be, 4.0) + 2.0 * (-3.0 * b * b + le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = a / (3.0 * le * le * dis) * (-2.0 * pow(a, 4.0) - a * a * le * le + pow(le, 4.0) + 2.0 * dis * pow(abs(a), 3.0));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_a_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * le / (3.0 * pow(be, 6.0) * dis) * (
			4.0 * pow(a, 4.0) * (-b * b + le * le) + pow(be, 4.0) * (be * be + b * b)
			+ 2.0 * a * a * (-pow(b, 4.0) + pow(le, 4.0)) + 4.0 * (b * b - le * le) * dis * pow(abs(a), 3.0)
			);
	}
	else {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_b_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / 30.0) * (
			16.0 * pow(a, 4.0) * pow(le, 3.0) * (-5.0 * b * b + le * le) * abs(a) / pow(be, 8.0) - 15.0 * (a * a + 3.0 * b * b) * atanh(le / dis)
			+ le / (dis * pow(be, 8.0)) * (
				16.0 * pow(a, 6.0) * (5.0 * pow(b * le, 2.0) - pow(le, 4.0))
				+ 3.0 * pow(be, 6.0) * (15.0 * pow(b, 4.0) + 20.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0))
				+ a * a * pow(be, 4.0) * (15.0 * pow(b, 4.0) + 20.0 * pow(b * le, 2.0) + 17.0 * pow(le, 4.0))
				+ 8.0 * pow(a, 4.0) * (5.0 * pow(b * b * le, 2.0) + 4.0 * pow(b * le * le, 2.0) - pow(le, 6.0))
				)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = le / (10.0 * be * be * be) * (15.0 * b * b * b * b + 20.0 * b * b * le * le + 3.0 * le * le * le * le) - 1.5 * b * b * atanh(le / dis);
		}
		else {
			result = 0.3 * le * abs(le);
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 1.0 / (30.0 * le * le * le) * (dis * (-16.0 * pow(a, 4.0) + 8.0 * a * a * le * le + 9.0 * pow(le, 4.0))
				+ 16.0 * pow(a, 4.0) * abs(a) - 15.0 * a * a * le * le * le * atanh(le / dis));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_b_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b * le * le / (5.0 * pow(be, 8.0) * dis) * (
			8.0 * pow(a, 6.0) * (le * le - b * b) + a * a * (b * b - le * le) * pow(be, 4.0)
			+ pow(be, 6.0) * (2.0 * b * b + 3.0 * le * le) - 4.0 * pow(a, 4.0) * (b * b - le * le) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * le * le * (2.0 * b * b + 3.0 * le * le) / (5.0 * be * be * be);
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

double Quadratic_Polyhedral::Psi_I6_b_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = le / (15.0 * pow(be, 8.0) * dis) * (
			4.0 * pow(a, 6.0) * (3.0 * pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + pow(le, 4.0))
			+ 3.0 * pow(be, 6.0) * (4.0 * pow(b, 4.0) + 6.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0))
			+ a * a * pow(be, 4.0) * (6.0 * pow(b, 4.0) + 19.0 * pow(b * le, 2.0) + 7.0 * pow(le, 4.0))
			+ 2.0 * pow(a, 4.0) * (3.0 * pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + pow(le, 4.0)) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = le * (4.0 * b * b * b * b + 6.0 * b * b * le * le + 3.0 * pow(le, 4.0)) / (5.0 * be * be * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = (4.0 * pow(a, 6.0) + 2.0 * a * a * a * a * le * le + 7.0 * a * a * le * le * le * le + 9.0 * pow(le, 6.0) - 4.0 * a * a * a * a * dis * abs(a)) / (15.0 * le * le * le * dis);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Psi_I6_le_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = b / (15.0 * pow(be, 8.0) * dis) * (
			-4.0 * pow(a, 6.0) * (pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0))
			+ 3.0 * pow(be, 6.0) * (2.0 * pow(b, 4.0) + 9.0 * pow(b * le, 2.0) + 6.0 * pow(le, 4.0))
			+ a * a * pow(be, 4.0) * (8.0 * pow(b, 4.0) + 11.0 * pow(b * le, 2.0) + 9.0 * pow(le, 4.0))
			- 2.0 * pow(a, 4.0) * (pow(b, 4.0) - 8.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0)) * (be * be - 2.0 * dis * abs(a))
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * (2.0 * pow(b, 4.0) + 9.0 * b * b * le * le + 6.0 * le * le * le * le) / (5.0 * be * be * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = -2.0 * (2.0 * a * a - 3.0 * b * b) * dis * dis * dis / (15.0 * b * b * b) + 4.0 * pow(a, 4.0) * abs(a) / (15.0 * b * b * b);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I4_III(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (6.0 * be * be)) * (
			-2.0 * a * a * b * le * dis + 2.0 * a * a * abs(a) * (b * le + be * be * (atan(le / b) - atan(le * abs(a) / (b * dis))))
			+ 2.0 * pow(b, 3.0) * be * be * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * b * b * atanh(le / be) / 3.0;
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

double Quadratic_Polyhedral::Phi_I4_III_a(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * (
			-b * le * dis / pow(be, 2.0) + abs(a) * (b * le / (be * be) + atan(le / b)) - a * atan(a * le / (b * dis))
			);
	}
	else {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I4_III_b(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b * b / (3.0 * pow(be, 4.0))) * (
			-le * (-3.0 * a * a + dis * dis) * dis - 2.0 * le * pow(abs(a), 3.0) + 3.0 * pow(be, 4.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -b * b * le / (3.0 * be) + b * b * atanh(le / dis);
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

double Quadratic_Polyhedral::Phi_I4_III_le(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (pow(b, 3.0) / (3.0 * pow(be, 4.0) * dis)) * (
			-2.0 * pow(a, 4.0) - a * a * be * be + pow(be, 4.0) + 2.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * b / (be * 3.0);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = ((-2.0 * a * a + b * b) * dis + 2.0 * pow(abs(a), 3.0)) / (3.0 * b);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I5_III(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b * b / (3.0 * be * be)) * (pow(dis, 3.0) - pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		result = b * b * be / 3.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (1.0 / 3.0) * (dis * dis * dis - abs(a) * abs(a) * abs(a));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}


double Quadratic_Polyhedral::Phi_I5_III_a(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * b * (dis - abs(a)) / (be * be);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * (dis - abs(a));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I5_III_b(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b / (3.0 * pow(be, 4.0))) * (
			dis * (3.0 * pow(b, 4.0) + 5.0 * pow(b * le, 2.0) + 2.0 * le * le * (a * a + le * le)) - 2.0 * le * le * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * (3.0 * b * b + 2.0 * le * le) / (3.0 * be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = b * dis;
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I5_III_le(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b * b * le / (3.0 * pow(be, 4.0))) * (
			(-3.0 * a * a + dis * dis) * dis + 2.0 * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * b * le / (3.0 * be);
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

double Quadratic_Polyhedral::Phi_I6_III(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (6.0 * be * be)) * (
			b * le * dis * (dis * dis + a * a) + pow(abs(a), 3.0) * (-2.0 * b * le + 2.0 * be * be * (atan(le / b) - atan(abs(a) * le / (b * dis))))
			- b * (3.0 * a * a + b * b) * be * be * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = (b * le * be - b * b * b * atanh(le / be)) / 6.0;
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

double Quadratic_Polyhedral::Phi_I6_III_a(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = -(a / (be * be)) * (
			-b * le * dis + abs(a) * (b * le - be * be * atan(le / b))
			+ be * be * (a * atan(a * le / (b * dis)) + b * atanh(le / dis))
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		result = 0.0;
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I6_III_b(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / (6.0 * pow(be, 4.0))) * (
			le * dis * (3.0 * pow(b, 4.0) + 4.0 * (a * a + b * b) * le * le + pow(le, 4.0)) - 4.0 * pow(le * abs(a), 3.0) - 3.0 * (a * a + b * b) * pow(be, 4.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = (1.0 / 6.0) * (le * (3.0 * b * b + le * le) / be - 3.0 * b * b * atanh(le / be));
		}
		else {
			result = le * abs(le) / 6.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = (dis * (4.0 * a * a + le * le) - 4.0 * pow(abs(a), 3.0) - 3.0 * a * a * le * atanh(le / dis)) / (6.0 * le);
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Quadratic_Polyhedral::Phi_I6_III_le(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (b * le * le / (3.0 * pow(be, 4.0) * dis)) * (
			-2.0 * pow(a, 4.0) - a * a * be * be
			+ pow(be, 4.0) + 2.0 * dis * pow(abs(a), 3.0)
			);
	}
	else if (a == 0) {
		if (b != 0 || le != 0) {
			result = b * le * le / (3.0 * be);
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

double Quadratic_Polyhedral::I7(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Phi_I4_III(a, b, lp) - Phi_I4_III(a, b, lm));
}

double Quadratic_Polyhedral::I8(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Phi_I5_III(a, b, lp) - Phi_I5_III(a, b, lm));
}

double Quadratic_Polyhedral::I9(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Phi_I6_III(a, b, lp) - Phi_I6_III(a, b, lm));
}

double Quadratic_Polyhedral::I7_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j] * (Phi_I4_III_a(a, b, lp) - Phi_I4_III_a(a, b, lm))
		- Vnorm[j] * (Phi_I4_III_b(a, b, lp) - Phi_I4_III_b(a, b, lm))
		- Vdir[j] * (Phi_I4_III_le(a, b, lp) - Phi_I4_III_le(a, b, lm))
		);


	return result;
}

double Quadratic_Polyhedral::I8_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j] * (Phi_I5_III_a(a, b, lp) - Phi_I5_III_a(a, b, lm))
		- Vnorm[j] * (Phi_I5_III_b(a, b, lp) - Phi_I5_III_b(a, b, lm))
		- Vdir[j] * (Phi_I5_III_le(a, b, lp) - Phi_I5_III_le(a, b, lm))
		);


	return result;
}


double Quadratic_Polyhedral::I9_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j] * (Phi_I6_III_a(a, b, lp) - Phi_I6_III_a(a, b, lm))
		- Vnorm[j] * (Phi_I6_III_b(a, b, lp) - Phi_I6_III_b(a, b, lm))
		- Vdir[j] * (Phi_I6_III_le(a, b, lp) - Phi_I6_III_le(a, b, lm))
		);


	return result;
}


double Quadratic_Polyhedral::I10_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I4_a_a(a, b, lp) - Psi_I4_a_a(a, b, lm)) - Vnorm[k] * (Psi_I4_a_b(a, b, lp) - Psi_I4_a_b(a, b, lm)) - Vdir[k] * (Psi_I4_a_le(a, b, lp) - Psi_I4_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I4_a_b(a, b, lp) - Psi_I4_a_b(a, b, lm)) - Vnorm[k] * (Psi_I4_b_b(a, b, lp) - Psi_I4_b_b(a, b, lm)) - Vdir[k] * ((Psi_I4_b_le(a, b, lp) - Psi_I4_b_le(a, b, lm))))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I4_a_le(a, b, lp) - Psi_I4_a_le(a, b, lm)) - Vnorm[k] * ((Psi_I4_b_le(a, b, lp) - Psi_I4_b_le(a, b, lm))) - Vdir[k] * ((Psi_I4_le_le(a, b, lp) - Psi_I4_le_le(a, b, lm))))
		);

	return result;
}

double Quadratic_Polyhedral::I11_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I5_a_a(a, b, lp) - Psi_I5_a_a(a, b, lm)) - Vnorm[k] * (Psi_I5_a_b(a, b, lp) - Psi_I5_a_b(a, b, lm)) - Vdir[k] * (Psi_I5_a_le(a, b, lp) - Psi_I5_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I5_a_b(a, b, lp) - Psi_I5_a_b(a, b, lm)) - Vnorm[k] * (Psi_I5_b_b(a, b, lp) - Psi_I5_b_b(a, b, lm)) - Vdir[k] * ((Psi_I5_b_le(a, b, lp) - Psi_I5_b_le(a, b, lm))))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I5_a_le(a, b, lp) - Psi_I5_a_le(a, b, lm)) - Vnorm[k] * ((Psi_I5_b_le(a, b, lp) - Psi_I5_b_le(a, b, lm))) - Vdir[k] * ((Psi_I5_le_le(a, b, lp) - Psi_I5_le_le(a, b, lm))))
		);

	return result;
}

double Quadratic_Polyhedral::I12_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I6_a_a(a, b, lp) - Psi_I6_a_a(a, b, lm)) - Vnorm[k] * (Psi_I6_a_b(a, b, lp) - Psi_I6_a_b(a, b, lm)) - Vdir[k] * (Psi_I6_a_le(a, b, lp) - Psi_I6_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I6_a_b(a, b, lp) - Psi_I6_a_b(a, b, lm)) - Vnorm[k] * (Psi_I6_b_b(a, b, lp) - Psi_I6_b_b(a, b, lm)) - Vdir[k] * ((Psi_I6_b_le(a, b, lp) - Psi_I6_b_le(a, b, lm))))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I6_a_le(a, b, lp) - Psi_I6_a_le(a, b, lm)) - Vnorm[k] * ((Psi_I6_b_le(a, b, lp) - Psi_I6_b_le(a, b, lm))) - Vdir[k] * ((Psi_I6_le_le(a, b, lp) - Psi_I6_le_le(a, b, lm))))
		);

	return result;
}

double Quadratic_Polyhedral::I10_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{

	double result = -(Psi_I4_a_a_a(a, b, lp) - Psi_I4_a_a_a(a, b, lm)) * Svnorm[i] * Svnorm[j] * Svnorm[k] - (Psi_I4_b_b_b(a, b, lp) - Psi_I4_b_b_b(a, b, lm)) * Vnorm[i] * Vnorm[j] * Vnorm[k]
		- (Psi_I4_le_le_le(a, b, lp) - Psi_I4_le_le_le(a, b, lm)) * Vdir[i] * Vdir[j] * Vdir[k] - (Psi_I4_a_a_b(a, b, lp) - Psi_I4_a_a_b(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vnorm[k] + Svnorm[i] * Vnorm[j] * Svnorm[k] + Vnorm[i] * Svnorm[j] * Svnorm[k])
		- (Psi_I4_a_a_le(a, b, lp) - Psi_I4_a_a_le(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vdir[k] + Svnorm[i] * Vdir[j] * Svnorm[k] + Vdir[i] * Svnorm[j] * Svnorm[k]) - (Psi_I4_a_b_b(a, b, lp) - Psi_I4_a_b_b(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vnorm[k] + Svnorm[j] * Vnorm[i] * Vnorm[k] + Svnorm[k] * Vnorm[j] * Vnorm[i])
		- (Psi_I4_b_b_le(a, b, lp) - Psi_I4_b_b_le(a, b, lm)) * (Vdir[i] * Vnorm[j] * Vnorm[k] + Vdir[j] * Vnorm[i] * Vnorm[k] + Vdir[k] * Vnorm[i] * Vnorm[j]) - (Psi_I4_a_le_le(a, b, lp) - Psi_I4_a_le_le(a, b, lm)) * (Svnorm[i] * Vdir[j] * Vdir[k] + Svnorm[j] * Vdir[i] * Vdir[k] + Svnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I4_b_le_le(a, b, lp) - Psi_I4_b_le_le(a, b, lm)) * (Vnorm[i] * Vdir[j] * Vdir[k] + Vnorm[j] * Vdir[i] * Vdir[k] + Vnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I4_a_b_le(a, b, lp) - Psi_I4_a_b_le(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vdir[k] + Svnorm[i] * Vnorm[k] * Vdir[j] + Svnorm[j] * Vnorm[i] * Vdir[k] + Svnorm[k] * Vnorm[i] * Vdir[j] + Svnorm[j] * Vnorm[k] * Vdir[i] + Svnorm[k] * Vnorm[j] * Vdir[i]);

	return result;
}


double Quadratic_Polyhedral::I11_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = -(Psi_I5_a_a_a(a, b, lp) - Psi_I5_a_a_a(a, b, lm)) * Svnorm[i] * Svnorm[j] * Svnorm[k] - (Psi_I5_b_b_b(a, b, lp) - Psi_I5_b_b_b(a, b, lm)) * Vnorm[i] * Vnorm[j] * Vnorm[k]
		- (Psi_I5_le_le_le(a, b, lp) - Psi_I5_le_le_le(a, b, lm)) * Vdir[i] * Vdir[j] * Vdir[k] - (Psi_I5_a_a_b(a, b, lp) - Psi_I5_a_a_b(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vnorm[k] + Svnorm[i] * Vnorm[j] * Svnorm[k] + Vnorm[i] * Svnorm[j] * Svnorm[k])
		- (Psi_I5_a_a_le(a, b, lp) - Psi_I5_a_a_le(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vdir[k] + Svnorm[i] * Vdir[j] * Svnorm[k] + Vdir[i] * Svnorm[j] * Svnorm[k]) - (Psi_I5_a_b_b(a, b, lp) - Psi_I5_a_b_b(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vnorm[k] + Svnorm[j] * Vnorm[i] * Vnorm[k] + Svnorm[k] * Vnorm[j] * Vnorm[i])
		- (Psi_I5_b_b_le(a, b, lp) - Psi_I5_b_b_le(a, b, lm)) * (Vdir[i] * Vnorm[j] * Vnorm[k] + Vdir[j] * Vnorm[i] * Vnorm[k] + Vdir[k] * Vnorm[i] * Vnorm[j]) - (Psi_I5_a_le_le(a, b, lp) - Psi_I5_a_le_le(a, b, lm)) * (Svnorm[i] * Vdir[j] * Vdir[k] + Svnorm[j] * Vdir[i] * Vdir[k] + Svnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I5_b_le_le(a, b, lp) - Psi_I5_b_le_le(a, b, lm)) * (Vnorm[i] * Vdir[j] * Vdir[k] + Vnorm[j] * Vdir[i] * Vdir[k] + Vnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I5_a_b_le(a, b, lp) - Psi_I5_a_b_le(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vdir[k] + Svnorm[i] * Vnorm[k] * Vdir[j] + Svnorm[j] * Vnorm[i] * Vdir[k] + Svnorm[k] * Vnorm[i] * Vdir[j] + Svnorm[j] * Vnorm[k] * Vdir[i] + Svnorm[k] * Vnorm[j] * Vdir[i]);


	return result;
}

double Quadratic_Polyhedral::I12_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = -(Psi_I6_a_a_a(a, b, lp) - Psi_I6_a_a_a(a, b, lm)) * Svnorm[i] * Svnorm[j] * Svnorm[k] - (Psi_I6_b_b_b(a, b, lp) - Psi_I6_b_b_b(a, b, lm)) * Vnorm[i] * Vnorm[j] * Vnorm[k]
		- (Psi_I6_le_le_le(a, b, lp) - Psi_I6_le_le_le(a, b, lm)) * Vdir[i] * Vdir[j] * Vdir[k] - (Psi_I6_a_a_b(a, b, lp) - Psi_I6_a_a_b(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vnorm[k] + Svnorm[i] * Vnorm[j] * Svnorm[k] + Vnorm[i] * Svnorm[j] * Svnorm[k])
		- (Psi_I6_a_a_le(a, b, lp) - Psi_I6_a_a_le(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vdir[k] + Svnorm[i] * Vdir[j] * Svnorm[k] + Vdir[i] * Svnorm[j] * Svnorm[k]) - (Psi_I6_a_b_b(a, b, lp) - Psi_I6_a_b_b(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vnorm[k] + Svnorm[j] * Vnorm[i] * Vnorm[k] + Svnorm[k] * Vnorm[j] * Vnorm[i])
		- (Psi_I6_b_b_le(a, b, lp) - Psi_I6_b_b_le(a, b, lm)) * (Vdir[i] * Vnorm[j] * Vnorm[k] + Vdir[j] * Vnorm[i] * Vnorm[k] + Vdir[k] * Vnorm[i] * Vnorm[j]) - (Psi_I6_a_le_le(a, b, lp) - Psi_I6_a_le_le(a, b, lm)) * (Svnorm[i] * Vdir[j] * Vdir[k] + Svnorm[j] * Vdir[i] * Vdir[k] + Svnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I6_b_le_le(a, b, lp) - Psi_I6_b_le_le(a, b, lm)) * (Vnorm[i] * Vdir[j] * Vdir[k] + Vnorm[j] * Vdir[i] * Vdir[k] + Vnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I6_a_b_le(a, b, lp) - Psi_I6_a_b_le(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vdir[k] + Svnorm[i] * Vnorm[k] * Vdir[j] + Svnorm[j] * Vnorm[i] * Vdir[k] + Svnorm[k] * Vnorm[i] * Vdir[j] + Svnorm[j] * Vnorm[k] * Vdir[i] + Svnorm[k] * Vnorm[j] * Vdir[i]);

	return result;
}


// complete set of potentials

double Quadratic_Polyhedral::PHI_pq_i(int p, int q, int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;

	result = Svnorm[p] * Svnorm[q] * a * a * U_PL.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (-Svnorm[i]) * a * (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * L_PL.I2(a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (-Svnorm[i]) * a * (Vdir[p] * Svnorm[q] + Vdir[q] * Svnorm[p]) * L_PL.I3(a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (-Svnorm[i]) * Vnorm[p] * Vnorm[q] * I7(a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (-Svnorm[i]) * (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * I8(a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (-Svnorm[i]) * Vdir[p] * Vdir[q] * I9(a, b, lm, lp, Svnorm, Vnorm, Vdir);


	result = result - (d[p][i] * x[q] + d[q][i] * x[p]) * U_PL.PHI(a, b, lm, lp)
		- x[p] * x[q] * U_PL.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[p][i] * L_PL.PHI_p(q, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ x[p] * L_PL.PHI_p_i(q, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[q][i] * L_PL.PHI_p(p, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ x[q] * L_PL.PHI_p_i(p, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

	return result;
}

double Quadratic_Polyhedral::PHI_pq_ij(int p, int q, int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	result = Svnorm[p] * Svnorm[q] * (
		2.0 * a * (-Svnorm[j]) * U_PL.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * a * U_PL.PHI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir))
		+ (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * (
			Svnorm[i] * Svnorm[j] * L_PL.I2(a, b, lm, lp, Svnorm, Vnorm, Vdir)
			- a * Svnorm[i] * L_PL.I2_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
			)
		+ (Vdir[p] * Svnorm[q] + Vdir[q] * Svnorm[p]) * (
			Svnorm[i] * Svnorm[j] * L_PL.I3(a, b, lm, lp, Svnorm, Vnorm, Vdir)
			- a * Svnorm[i] * L_PL.I3_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
			);
	result += Vnorm[p] * Vnorm[q] * (-Svnorm[i]) * I7_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * I8_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vdir[p] * Vdir[q] * I9_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	result += -x[p] * x[q] * U_PL.PHI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * x[q] + d[q][j] * x[p]) * U_PL.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[i][p] * x[q] + d[i][q] * x[p]) * U_PL.PHI_i(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * d[q][i] + d[q][j] * d[p][i]) * U_PL.PHI(a, b, lm, lp)
		+ d[i][p] * L_PL.PHI_p_i(q, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[p][j] * L_PL.PHI_p_i(q, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ x[p] * L_PL.PHI_p_ij(q, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[i][q] * L_PL.PHI_p_i(p, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[q][j] * L_PL.PHI_p_i(p, i, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ x[q] * L_PL.PHI_p_ij(p, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
	return result;
}

double Quadratic_Polyhedral::PSI_pq_ijk(int p, int q, int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;

	result = Svnorm[p] * Svnorm[q] * (
		2.0 * Svnorm[j] * Svnorm[k] * U_PL.PSI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ 2.0 * (-Svnorm[j]) * a * U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ 2.0 * (-Svnorm[k]) * a * U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * a * U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result += -(Svnorm[i]) * (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * (
		-Svnorm[j] * L_PL.I5_j(k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[k] * L_PL.I5_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * L_PL.I5_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result += -(Svnorm[i]) * (Vdir[p] * Svnorm[q] + Svnorm[p] * Vdir[q]) * (
		-Svnorm[j] * L_PL.I6_j(k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[k] * L_PL.I6_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * L_PL.I6_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result += -(Svnorm[i]) * (Vnorm[p] * Vnorm[q] * I10_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * I11_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ Vdir[p] * Vdir[q] * I12_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result = result - x[p] * x[q] * U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir) - (d[p][k] * x[q] + d[q][k] * x[p]) * U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * U_PL.PSI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir) - (d[p][j] * x[q] + d[q][j] * x[p]) * U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[i][p] * d[q][k] + d[i][q] * d[p][k]) * U_PL.PSI_i(j, a, b, lm, lp, Svnorm, Vnorm, Vdir) - (d[i][p] * x[q] + d[i][q] * x[p]) * U_PL.PSI_ij(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * d[q][i] + d[q][j] * d[p][i]) * U_PL.PSI_i(k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[i][p] * L_PL.PSI_p_ij(q, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + d[p][j] * L_PL.PSI_p_ij(q, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[p][k] * L_PL.PSI_p_ij(q, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + x[p] * L_PL.PSI_p_ijk(q, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[i][q] * L_PL.PSI_p_ij(p, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + d[q][j] * L_PL.PSI_p_ij(p, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[q][k] * L_PL.PSI_p_ij(p, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + x[q] * L_PL.PSI_p_ijk(p, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

	return result;
}

double Quadratic_Polyhedral::PSI_pq_ijkl(int p, int q, int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;

	result = (Svnorm[p] * Svnorm[q]) * 2.0 * (
		Svnorm[j] * Svnorm[k] * U_PL.PSI_ij(i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ Svnorm[j] * Svnorm[l] * U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ Svnorm[k] * Svnorm[l] * U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * (-Svnorm[j] * U_PL.PSI_ijk(i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
			- Svnorm[k] * U_PL.PSI_ijk(i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
			- Svnorm[l] * U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir))
		+ a * a * 0.5 * U_PL.PSI_ijkl(i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result += -Svnorm[i] * (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * (
		-Svnorm[j] * L_PL.I5_jk(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[k] * L_PL.I5_jk(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[l] * L_PL.I5_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * L_PL.I5_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result += -Svnorm[i] * (Vdir[p] * Svnorm[q] + Svnorm[p] * Vdir[q]) * (
		-Svnorm[j] * L_PL.I6_jk(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[k] * L_PL.I6_jk(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[l] * L_PL.I6_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ a * L_PL.I6_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);

	result += -(Svnorm[i]) * (Vnorm[p] * Vnorm[q] * I10_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * I11_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ Vdir[p] * Vdir[q] * I12_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		);


	result += -x[p] * x[q] * U_PL.PSI_ijkl(i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][l] * x[q] + d[q][l] * x[p]) * U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][k] * x[q] + d[q][k] * x[p]) * U_PL.PSI_ijk(i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * x[q] + d[q][j] * x[p]) * U_PL.PSI_ijk(i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][i] * x[q] + d[q][i] * x[p]) * U_PL.PSI_ijk(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][k] * d[q][l] + d[q][k] * d[p][l]) * U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * U_PL.PSI_ij(i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * d[q][l] + d[q][j] * d[p][l]) * U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[i][p] * d[q][k] + d[i][q] * d[p][k]) * U_PL.PSI_ij(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[i][p] * d[q][l] + d[q][i] * d[p][l]) * U_PL.PSI_ij(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (d[p][j] * d[q][i] + d[q][j] * d[p][i]) * U_PL.PSI_ij(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[i][p] * L_PL.PSI_p_ijk(q, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + d[p][j] * L_PL.PSI_p_ijk(q, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[p][k] * L_PL.PSI_p_ijk(q, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + d[p][l] * L_PL.PSI_p_ijk(q, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ x[p] * L_PL.PSI_p_ijkl(q, i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[i][q] * L_PL.PSI_p_ijk(p, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + d[q][j] * L_PL.PSI_p_ijk(p, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ d[q][k] * L_PL.PSI_p_ijk(p, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) + d[q][l] * L_PL.PSI_p_ijk(p, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x)
		+ x[q] * L_PL.PSI_p_ijkl(p, i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

	return result;
}

void Quadratic_Polyhedral::PSI_pq_ijk_combine(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, \
	double* x, double* linear_result, double* result, double ij, double ik, double jk, double kj, double ijk, double jk_5, double jk_6)
{
	int p, q;

	// uniform potentials
	double ii = U_PL.PSI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir); double jj = U_PL.PSI_i(j, a, b, lm, lp, Svnorm, Vnorm, Vdir); double kk = U_PL.PSI_i(k, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	// linear potentials:
	double j_5 = L_PL.I5_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir); double j_6 = L_PL.I6_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	double k_5 = L_PL.I5_j(k, a, b, lm, lp, Svnorm, Vnorm, Vdir); double k_6 = L_PL.I6_j(k, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	// rewrite in new double:
	double* linear_ij, * linear_ik, * linear_jk;
	linear_ij = new double[3]; linear_ik = new double[3]; linear_jk = new double[3];

	linear_ij[0] = L_PL.PSI_p_ij(0, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_ij[1] = L_PL.PSI_p_ij(1, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_ij[2] = L_PL.PSI_p_ij(2, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
	linear_ik[0] = L_PL.PSI_p_ij(0, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_ik[1] = L_PL.PSI_p_ij(1, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_ik[2] = L_PL.PSI_p_ij(2, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
	linear_jk[0] = L_PL.PSI_p_ij(0, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_jk[1] = L_PL.PSI_p_ij(1, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_jk[2] = L_PL.PSI_p_ij(2, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);


	//double linear_ij[3] = { PSI_p_ij(0, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ij(1, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x),PSI_p_ij(2, i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };
	//double linear_ik[3] = { PSI_p_ij(0, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ij(1, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x),PSI_p_ij(2, i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };
	//double linear_jk[3] = { PSI_p_ij(0, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ij(1, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x),PSI_p_ij(2, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };

	// quadratic potentials:
	double jk_10 = I10_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); double jk_11 = I11_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); double jk_12 = I12_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	for (int s = 0; s < 6; s++) {
		switch (s) {
		case 0:
			p = 0; q = 1;
			break;
		case 1:
			p = 0; q = 2;
			break;
		case 2:
			p = 1; q = 2;
			break;
		case 3:
			p = 0; q = 0;
			break;
		case 4:
			p = 1; q = 1;
			break;
		case 5:
			p = 2; q = 2;
			break;
		}

		result[s] = Svnorm[p] * Svnorm[q] * (2.0 * Svnorm[j] * Svnorm[k] * ii + 2.0 * (-Svnorm[j]) * a * ik + 2.0 * (-Svnorm[k]) * a * ij + a * a * ijk)
			- (Svnorm[i]) * (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * (-Svnorm[j] * k_5 - Svnorm[k] * j_5 + a * jk_5)
			- (Svnorm[i]) * (Vdir[p] * Svnorm[q] + Svnorm[p] * Vdir[q]) * (-Svnorm[j] * k_6 - Svnorm[k] * j_6 + a * jk_6)
			- Svnorm[i] * (Vnorm[p] * Vnorm[q] * jk_10 + (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * jk_11 + Vdir[p] * Vdir[q] * jk_12)
			- x[p] * x[q] * ijk - (d[p][k] * x[q] + d[q][k] * x[p]) * ij - (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * ii
			- (d[p][j] * x[q] + d[q][j] * x[p]) * ik - (d[i][p] * d[q][k] + d[i][q] * d[p][k]) * jj - (d[i][p] * x[q] + d[i][q] * x[p]) * jk
			- (d[p][j] * d[q][i] + d[q][j] * d[p][i]) * kk + d[i][p] * linear_jk[q] + d[p][j] * linear_ik[q] + d[p][k] * linear_ij[q] + x[p] * linear_result[q]
			+ d[i][q] * linear_jk[p] + d[q][j] * linear_ik[p] + d[q][k] * linear_ij[p] + x[q] * linear_result[p];
	}

	delete[] linear_ij; delete[] linear_jk; delete[] linear_ik;

}


void Quadratic_Polyhedral::PSI_pq_ijkl_combine(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* linear_result, double* result, double ijk, double ikl, double ijl, \
	double jkl, double kjl, double ljk, double ijkl, double jkl_5, double jkl_6)
{
	double il = U_PL.PSI_ij(i, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); double ik = U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); double ij = U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	double jl = U_PL.PSI_ij(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); double jk = U_PL.PSI_ij(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); double kl = U_PL.PSI_ij(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	double jkl_10 = I10_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); double jkl_11 = I11_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); double jkl_12 = I12_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	//rewrite it in the form of new double;
	double* linear_result_jkl = new double[3]; double* linear_result_ikl = new double[3]; double* linear_result_ijl = new double[3]; double* linear_result_ijk = new double[3];

	linear_result_jkl[0] = L_PL.PSI_p_ijk(0, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_jkl[1] = L_PL.PSI_p_ijk(1, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_jkl[2] = L_PL.PSI_p_ijk(2, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
	linear_result_ikl[0] = L_PL.PSI_p_ijk(0, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_ikl[1] = L_PL.PSI_p_ijk(1, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_ikl[2] = L_PL.PSI_p_ijk(2, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
	linear_result_ijl[0] = L_PL.PSI_p_ijk(0, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_ijl[1] = L_PL.PSI_p_ijk(1, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_ijl[2] = L_PL.PSI_p_ijk(2, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);
	linear_result_ijk[0] = L_PL.PSI_p_ijk(0, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_ijk[1] = L_PL.PSI_p_ijk(1, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x); linear_result_ijk[2] = L_PL.PSI_p_ijk(2, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x);

	//double linear_result_jkl[3] = { PSI_p_ijk(0, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) , PSI_p_ijk(1, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ijk(2, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };
	//double linear_result_ikl[3] = { PSI_p_ijk(0, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) , PSI_p_ijk(1, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ijk(2, i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };
	//double linear_result_ijl[3] = { PSI_p_ijk(0, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) , PSI_p_ijk(1, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ijk(2, i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };
	//double linear_result_ijk[3] = { PSI_p_ijk(0, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) , PSI_p_ijk(1, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x), PSI_p_ijk(2, i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir, x) };

	double jk_5, jl_5, kl_5, jk_6, jl_6, kl_6;
	int p, q;
	jk_5 = L_PL.I5_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); jl_5 = L_PL.I5_jk(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); kl_5 = L_PL.I5_jk(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	jk_6 = L_PL.I6_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); jl_6 = L_PL.I6_jk(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); kl_6 = L_PL.I6_jk(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	for (int s = 0; s < 6; s++) {
		switch (s) {
		case 0:
			p = 0; q = 1;
			break;
		case 1:
			p = 0; q = 2;
			break;
		case 2:
			p = 1; q = 2;
			break;
		case 3:
			p = 0; q = 0;
			break;

		case 4:
			p = 1; q = 1;
			break;
		case 5:
			p = 2; q = 2;
			break;
			//default:
			//	cout << " check the input of p and q" << endl;
		}

		result[s] = (Svnorm[p] * Svnorm[q]) * 2.0 * (Svnorm[j] * Svnorm[k] * il + Svnorm[j] * Svnorm[l] * ik + Svnorm[k] * Svnorm[l] * ij + a * (-Svnorm[j] * ikl - Svnorm[k] * ijl - Svnorm[l] * ijk) + a * a * 0.5 * ijkl)
			- Svnorm[i] * (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * (-Svnorm[j] * kl_5 - Svnorm[k] * jl_5 - Svnorm[l] * jk_5 + a * jkl_5)
			- Svnorm[i] * (Vdir[p] * Svnorm[q] + Svnorm[p] * Vdir[q]) * (-Svnorm[j] * kl_6 - Svnorm[k] * jl_6 - Svnorm[l] * jk_6 + a * jkl_6)
			- (Svnorm[i]) * (Vnorm[p] * Vnorm[q] * jkl_10 + (Vnorm[p] * Vdir[q] + Vnorm[q] * Vdir[p]) * jkl_11 + Vdir[p] * Vdir[q] * jkl_12)
			- x[p] * x[q] * ijkl - (d[p][l] * x[q] + d[q][l] * x[p]) * ijk - (d[p][k] * x[q] + d[q][k] * x[p]) * ijl - (d[p][j] * x[q] + d[q][j] * x[p]) * ikl - (d[p][i] * x[q] + d[q][i] * x[p]) * jkl
			- (d[p][k] * d[q][l] + d[q][k] * d[p][l]) * ij - (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * il - (d[p][j] * d[q][l] + d[q][j] * d[p][l]) * ik - (d[i][p] * d[q][k] + d[i][q] * d[p][k]) * jl - (d[i][p] * d[q][l] + d[q][i] * d[p][l]) * jk - (d[p][j] * d[q][i] + d[q][j] * d[p][i]) * kl
			+ d[i][p] * linear_result_jkl[q] + d[p][j] * linear_result_ikl[q] + d[p][k] * linear_result_ijl[q] + d[p][l] * linear_result_ijk[q] + d[i][q] * linear_result_jkl[p] + d[q][j] * linear_result_ikl[p] + d[q][k] * linear_result_ijl[p] + d[q][l] * linear_result_ijk[p]
			+ x[p] * linear_result[q] + x[q] * linear_result[p];
	}

	delete[] linear_result_ijk; delete[] linear_result_ijl; delete[] linear_result_jkl; delete[] linear_result_ikl;
}

void Quadratic_Polyhedral::PSI_pq_ijkl_combine_II(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* linear_result, double* result)
{
	double jk_5, jl_5, kl_5, jk_6, jl_6, kl_6;
	int p, q;
	jk_5 = L_PL.I5_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); jl_5 = L_PL.I5_jk(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); kl_5 = L_PL.I5_jk(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	jk_6 = L_PL.I6_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir); jl_6 = L_PL.I6_jk(j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir); kl_6 = L_PL.I6_jk(k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	for (int s = 0; s < 6; s++) {
		switch (s) {
		case 0:
			p = 0; q = 1;
			break;
		case 1:
			p = 0; q = 2;
			break;
		case 2:
			p = 1; q = 2;
			break;
		case 3:
			p = 0; q = 0;
			break;

		case 4:
			p = 1; q = 1;
			break;
		case 5:
			p = 2; q = 2;
			break;
			//default:
			//	cout << " check the input of p and q" << endl;
		}

		result[s] = -Svnorm[i] * (Vnorm[p] * Svnorm[q] + Svnorm[p] * Vnorm[q]) * (-Svnorm[j] * kl_5 - Svnorm[k] * jl_5 - Svnorm[l] * jk_5)
			- Svnorm[i] * (Vdir[p] * Svnorm[q] + Svnorm[p] * Vdir[q]) * (-Svnorm[j] * kl_6 - Svnorm[k] * jl_6 - Svnorm[l] * jk_6);
	}

}