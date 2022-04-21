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
#include "Linear_Polyhedral_integral.h"
#include "GlobalConstant.h"

double Linear_Polyhedral::Psi_I2(double a, double b, double le)
{
	double result = 0.0;
	double sq = sqrt(1 + pow(le / b, 2.0));
	double dis = sqrt(a * a + b * b + le * le);
	result = 0.125 * le * pow(b, 2.0) * dis
		- 0.125 * le * pow(a, 5.0) * asinh(b * sq / a) / (b * sq * abs(a))
		+ 0.125 * le * pow(a * a + b * b, 2.0) * pow(abs(le), -1.0) * (
			-log(sqrt(a * a + b * b) / abs(le)) + log(1.0 + dis / abs(le))
			);

	return result;
}

// first partial derivative of the Psi_p
double Linear_Polyhedral::Psi_I2_a(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * a * (
			-a * le * abs(a) * asinh(be / a) / be
			+ (a * a + b * b) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 0.5 * a * a * (-le * abs(a) * asinh(be / a) / be
				+ a * atanh(le / dis));
		}
		else {
			result = 0.0;
		}
	}
	return result;
}

double Linear_Polyhedral::Psi_I2_b(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * pow(be, -3.0) * b * (
			le * be * dis * (-a * a + 2.0 * be * be) + pow(a, 3.0) * le * abs(a) * asinh(be / a)
			+ 4.0 * (a * a + b * b) * pow(be, 3.0) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = 0.25 * (b * le * be + 2.0 * b * b * b * atanh(le / be));
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

double Linear_Polyhedral::Psi_I2_le(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * pow(be, -3.0) * b * b * (
			be * dis * (a * a + 2.0 * be * be) - pow(a, 3.0) * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		result = 0.25 * b * b * be;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.125 * (dis * (a * a + 2.0 * b * b) - a * a * a * abs(a) * asinh(b / a) / b);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}



// Second partial derivative of the Psi_p
double Linear_Polyhedral::Psi_I2_a_a(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * (
			-3.0 * a * le * abs(a) * asinh(be / a) / be
			+ (3.0 * a * a + b * b) * atanh(le / dis)
			//- a * a * le / dis
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.5 * b * b * atanh(le / be);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 1.5 * a * (-abs(a) * asinh(le / a) + a * atanh(le / dis));
		}
		else {
			result = 0;
		}
	}
	return result;
}

double Linear_Polyhedral::Psi_I2_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * pow(be, -3.0) * a * b * (
			-le * be * dis + a * le * abs(a) * asinh(be / a) + 2.0 * pow(be, 3.0) * atanh(le / dis)
			);
	}
	else {
		result = 0.0;
	}


	return result;
}

double Linear_Polyhedral::Psi_I2_a_le(double a, double b, double le)
{
	double result = 0.0;
	double dis = sqrt(a * a + b * b + le * le); double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b * b / (2.0 * pow(be, 3.0)) * (
			be * dis - abs(a) * a * asinh(be / a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.5 * a / abs(b) * (be * dis - a * abs(a) * asinh(be / a));
		}
		else {
			result = 0.0;
		}
	}

	return result;
}


double Linear_Polyhedral::Psi_I2_b_b(double a, double b, double le)
{
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * (
			le * dis * (a * a * (2.0 * b * b - le * le) + 2.0 * le * le * be * be) / pow(be, 4.0)
			+ pow(a, 3.0) * le * (-2.0 * b * b + le * le) * abs(a) * asinh(be / a) / pow(be, 5.0)
			+ 4.0 * (a * a + 3.0 * b * b) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * (pow(le, 3.0) / be + 6.0 * b * b * atanh(le / be));
		}
		else {
			result = 0.25 * le * abs(le);
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = 0.125 * (
				dis * (-a * a + 2.0 * le * le) / le + a * a * a * be * abs(a) * asinh(be / a) / (le * le * le)
				+ 4.0 * a * a * atanh(le / dis)
				);
		}
		else if (b == 0 && le == 0) {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I2_b_le(double a, double b, double le)
{
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * b * pow(be, -5.0) * (
			be * dis * (-a * a * b * b + 6.0 * pow(b, 4.0) + 2.0 * (a * a + 5.0 * b * b) * le * le + 4.0 * pow(le, 4.0))
			+ pow(a, 3.0) * (b * b - 2.0 * le * le) * abs(a) * asinh(be / a)

			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * b * (3.0 * b * b + 2.0 * le * le) / be;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.125 / (b * be) * (be * dis * (-a * a + 6.0 * b * b) + a * a * a * abs(a) * asinh(be / a));
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I2_le_le(double a, double b, double le)
{

	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * pow(be, -5.0) * b * b * le * (
			dis * be * (-3.0 * a * a + 2.0 * be * be) + 3.0 * pow(a, 3.0) * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * b * b * le / be;
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


// third derivatives of Psi_I2

double Linear_Polyhedral::Psi_I2_a_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {

		result = a * b * b * le / (dis * (a * a + b * b)) - 3.0 * le * abs(a) * asinh(be / a) / be
			+ 3.0 * a * atanh(le / dis);
	}
	else if (a == 0) {
		result = 0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = -3.0 * le * abs(a) * asinh(be / a) / be + 3.0 * a * atanh(le / dis);
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I2_a_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * (
			-le / ((a * a + b * b) * be * be * dis) * (3.0 * pow(a, 4.0) + 6.0 * pow(a * b, 2.0) + pow(b, 4.0) + (3.0 * a * a + b * b) * le * le)
			+ 3.0 * a * le * abs(a) * asinh(be / a) * pow(be, -3.0)
			+ 2.0 * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -0.5 * b * le / be + b * atanh(le / be);
		}
		else {
			result = 0;
		}
	}
	else if (a != 0) {
		result = 0;
	}
	/*else if (a != 0) {
		result = 0.0;
	}*/
	return result;
}

double Linear_Polyhedral::Psi_I2_a_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * b * pow(be, -3.0) * pow(dis, -2.0) * (
			be * dis * (2.0 * a * a + dis * dis) - 3.0 * a * dis * dis * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = b * b * 0.5 / be;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.5 * ((3.0 * a * a + b * b) / dis - 3.0 * a * abs(a) * asinh(b / a) / b);
		}
		else {
			result = 0;
		}
	}
	return result;
}

double Linear_Polyhedral::Psi_I2_a_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * a * (
			8.0 * atanh(le / dis)
			- (4.0 * le / (dis * pow(be, 4.0) * (a * a + b * b))) * (
				2.0 * pow(b, 6.0) + 5.0 * pow(b * b * le, 2.0) + 3.0 * pow(b * le * le, 2.0)
				+ pow(a, 4.0) * (-2.0 * b * b + le * le) + a * a * (-2.0 * pow(b, 4.0) + 2.0 * pow(b * le, 2.0) + pow(le, 4.0))

				)
			+ 4.0 * a * le * (-2.0 * b * b + le * le) * abs(a) * asinh(be / a) / pow(be, 5.0)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = 0.125 * a * (
				8.0 * atanh(le / dis) - 4.0 * dis / le + 4.0 * a * be * abs(a) * asinh(be / a) / pow(le, 3.0)
				);
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I2_a_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * a * b * pow(be, -6.0) * pow(dis, -1.0) * (
			be * be * (pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0)
				- a * a * (b * b - 2.0 * le * le)
				)
			+ a * (b * b - 2.0 * le * le) * dis * be * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.5 * a / (b * b * dis) * (-a * a * b + b * b * b + a * dis * abs(a) * asinh(b / a));
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I2_a_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = -0.5 * a * b * b * le * pow(be, -5.0) * (
			(be / dis) * (dis * dis + 2.0 * a * a) - 3.0 * a * abs(a) * asinh(be / a)
			);
	}
	else {
		result = 0;
	}

	return result;

}

double Linear_Polyhedral::Psi_I2_b_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 3.0 * b * atanh(le / dis) + 0.375 * b * pow(a, 3.0) * le * (2.0 * b * b - 3.0 * le * le) * abs(a) * asinh(be / a) / pow(be, 7.0);
		result += -0.125 * b / ((a * a + b * b) * pow(be, 6.0) * dis) * le * (
			pow(a, 6.0) * (6.0 * b * b - 9.0 * le * le) + pow(a, 4.0) * (2.0 * b * b - 3.0 * le * le) * (4.0 * b * b + le * le)
			+ 2.0 * b * b * pow(be, 4.0) * (6.0 * b * b + 7.0 * le * le)
			+ a * a * be * be * (6.0 * pow(b, 4.0) + 7.0 * b * b * le * le + 6.0 * pow(le, 4.0))
			);
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = b * 0.25 / (be * be * be) * (-6.0 * b * b * le - 7.0 * le * le * le + 12.0 * be * be * be * atanh(le / be));
		}
		else {
			result = 0;
		}
	}

	else if (a != 0)
		result = 0;

	return result;

}

double Linear_Polyhedral::Psi_I2_b_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * pow(be, -7.0) * (
			(be / dis) * (pow(a, 4.0) * (2.0 * pow(b, 4.0) - 11.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
				+ 2.0 * pow(be, 4.0) * (6.0 * pow(b, 4.0) + 9.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
				+ a * a * be * be * (6.0 * pow(b, 4.0) + 7.0 * pow(b * le, 2.0) + 6.0 * pow(le, 4.0))
				)
			+ pow(a, 3.0) * (-2.0 * pow(b, 4.0) + 11.0 * pow(b * le, 2.0) - 2.0 * pow(le, 4.0)) * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * pow(be, -3.0) * (
				6.0 * pow(b, 4.0) + 9.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0)
				);
		}
		else {
			result = 0.5 * abs(le);
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.25 * pow(abs(b), -3.0) * (
				abs(b) / dis * (pow(a, 4.0) + 3.0 * a * a * b * b + 6.0 * pow(b, 4.0)) - pow(a, 3.0) * abs(a) * asinh(abs(b) / a)
				);
		}
		else if (b == 0 && le != 0) {
			result = 0.25 * pow(abs(le), -3.0) * (
				abs(le) / dis * (pow(a, 4.0) + 3.0 * a * a * le * le + 2.0 * pow(le, 4.0)) - pow(a, 3.0) * abs(a) * asinh(abs(le) / a)
				);
		}
		else if (b == 0 && le == 0) {
			result = 2.0 / 3.0 * abs(a);
		}
	}

	return result;

}

double Linear_Polyhedral::Psi_I2_b_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * b * le * (
			pow(be, -6.0) * pow(dis, -1.0) * (
				pow(a, 4.0) * (9.0 * b * b - 6.0 * le * le) + 2.0 * pow(be, 4.0) * (be * be + le * le)
				+ a * a * (3.0 * pow(b, 4.0) + pow(b * le, 2.0) - 2.0 * pow(le, 4.0))
				)
			+ pow(be, -7.0) * 3.0 * pow(a, 3.0) * (-3.0 * b * b + 2.0 * le * le) * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * pow(be, -3.0) * b * le * (be * be + le * le);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		result = 0;
	}


	return result;

}

double Linear_Polyhedral::Psi_I2_le_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * b * b * (
			pow(be, -6.0) * pow(dis, -1.0) * (
				-3.0 * pow(a, 4.0) * (b * b - 4.0 * le * le) + 2.0 * b * b * pow(be, 4.0)
				+ a * a * (-pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + 4.0 * pow(le, 4.0))
				)
			+ pow(be, -7.0) * 3.0 * pow(a, 3.0) * (b * b - 4.0 * le * le) * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * pow(be, -3.0) * pow(b, 4.0);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.125 / (b * b * b) * (b * dis * (-3.0 * a * a + 2.0 * b * b) + 3.0 * a * a * a * abs(a) * asinh(b / a));
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
	}

	return result;

}

// third kind area integral by P * Sin[t]
double Linear_Polyhedral::Psi_I3(double a, double b, double le)
{
	double sq = sqrt(1 + pow(le / b, 2.0));
	double dis = sqrt(a * a + b * b + le * le);

	double result = b * dis * dis * (5.0 * a * a + 2.0 * (b * b + le * le));
	result += 3.0 * pow(a, 5.0) * pow(abs(a), -1.0) * dis * asinh(b * sq / a) / sq;

	return result / (24.0 * dis);

	return result;
}

double Linear_Polyhedral::Psi_I3_a(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = a * b / (2.0 * dis) * (
			dis * dis + a * dis / be * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * b / (2.0 * dis) * (
				dis * dis + a * dis / be * abs(a) * asinh(be / a)
				);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_b(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (1.0 / 24.0) * pow(dis, -3.0) * (
			b * b * dis * dis * (10 * a * a + 4.0 * be * be + 4.0 * dis * dis + 3.0 * pow(a, 4.0) / (be * be) - 3.0 * pow(a, 5.0) * asinh(be / a) * abs(a) / (dis * pow(be, 3.0)))
			+ (a * a + le * le) * (dis * dis * (5.0 * a * a + 2.0 * be * be) + 3.0 * pow(a, 5.0) * dis / be * asinh(be / a) / abs(a))
			);
	}
	else if (a == 0) {
		result = (1.0 / 12.0) * be * (4.0 * b * b + le * le);
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = (5.0 * pow(a, 4.0) + 7.0 * a * a * le * le + 2.0 * le * le * le * le + 3.0 * pow(a, 5.0) * sqrt(1.0 / (a * a) + 1.0 / (le * le)) * asinh(abs(le) / a)) / (24.0 * dis);
		}
		else if (b != 0 && le == 0) {
			result = dis * dis * dis / 3.0;
		}
		else {
			result = abs(a) * abs(a) * abs(a) / 3.0;
		}
	}
	return result;
}

double Linear_Polyhedral::Psi_I3_le(double a, double b, double le)
{
	double result;
	double be = sqrt(b * b + le * le); double dis = sqrt(a * a + b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * b * le * pow(be, -2.0) * pow(dis, -3.0) * pow(abs(a), -1.0) * (
			pow(dis, 4.0) * (a * a + 2.0 * be * be) * abs(a)
			- pow(a, 5.0) * (be * dis + a * a * dis / be) * asinh(be / a)
			);
	}
	else if (a == 0) {
		result = 0.25 * b * le * be;
	}
	else if (a != 0) {
		result = 0.0;
	}
	return result;
}

double Linear_Polyhedral::Psi_I3_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * pow(dis, -1.0) * b * (
			dis * dis + 3.0 * a * (dis / be) * abs(a) * asinh(be / a)
			);
	}
	else if (a == 0) {
		result = 0.5 * b * be;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (0.5 * b / dis) * (dis * dis + 3.0 * a * abs(a) * sqrt(1.0 + (a * a / (b * b))) * asinh(be / a));
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * pow(be, -4.0) * pow(dis, -1.0) * (
			a * be * be * dis * dis * (b * b + be * be)
			+ pow(a, 4.0) * le * le * be * dis * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * dis;
		}
		else if (b == 0 && le != 0) {
			result = 0.5 / (le * le * dis) * (a * le * le * dis * dis + pow(a, 4.0) * sqrt(le * le + (pow(le, 4.0) / (a * a))) * asinh(be / a));
		}
		else if (b == 0 && le == 0) {
			result = a * abs(a);
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (0.5 * a * b * le / (pow(be, 4.0) * dis * abs(a))) * (
			be * be * dis * dis * abs(a)
			- pow(a, 3.0) * be * dis * asinh(be / a)

			);
	}
	else {
		result = 0;
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (0.125 / (pow(be, 6.0) * dis)) * (
			b * dis * dis * be * be * (8.0 * pow(b, 4.0) + (3.0 * a * a + 14.0 * b * b) * le * le + 6.0 * pow(le, 4.0))
			- 3.0 * pow(a, 5.0) * b * le * le * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * b * (3.0 * be * be + b * b) / be;
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

double Linear_Polyhedral::Psi_I3_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (0.125 / (pow(be, 6.0) * dis)) * (
			le * dis * dis * be * be * (a * a * (-2.0 * b * b + le * le) + 2.0 * be * be * (be * be + b * b))
			- pow(a, 5.0) * le * (-2.0 * b * b + le * le) * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * le * (be * be + b * b) / be;
		}
		else {
			result = 0.25 * le * abs(le);
		}

	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 0.125 / (le * le * le * dis) * (
				pow(a * a * le, 2.0) + 3.0 * pow(a * le * le, 2.0) + 2.0 * pow(le, 6.0) - pow(a, 5.0) * sqrt(le * le + (pow(le, 4.0) / (a * a))) * asinh(be / a)
				);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (0.125 / (pow(be, 6.0) * dis)) * (
			b * be * be * dis * dis * (a * a * (b * b - 2.0 * le * le) + 2.0 * be * be * (b * b + 2.0 * le * le))
			- pow(a, 5.0) * b * (b * b - 2.0 * le * le) * dis * be * asinh(be / a) / abs(a)
			);
	}

	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * b * (2.0 * le * le + b * b) / be;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (0.125 / (b * b * b * dis)) * (
				pow(a * a * b, 2.0) + 3.0 * pow(a * b * b, 2.0) + 2.0 * pow(b, 6.0)
				- pow(a, 5.0) * sqrt(b * b + (pow(b, 4.0) / (a * a))) * asinh(be / a)
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

double Linear_Polyhedral::Psi_I3_a_a_a(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = -a * b / (be * be * dis) * (
			be * be - 3.0 * be * dis * asinh(be / a) * a / abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = -a / (b * dis) * (b * b - 3.0 * be * dis * asinh(be / a) * a / abs(a));
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_a_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 / (dis * pow(be, 4.0)) * (
			be * be * (2.0 * pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + pow(le, 4.0) + a * a * (4.0 * b * b + le * le))
			+ 3.0 * pow(a, 3.0) * pow(le, 2.0) * be * dis * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.5 * (be * be + b * b) / be;
		}
		else {
			result = 0.5 * abs(le);
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (2.0 * a * a + b * b) / dis;
		}
		else if (b == 0 && le != 0) {
			result = (0.5 / (le * le * dis)) * (a * a * le * le + pow(le, 4.0) + 3.0 * a * a * a * sqrt(le * le + pow(le, 4.0) / (a * a)) * asinh(abs(le) / a));
		}
		else if (b == 0 && le == 0) {
			result = 2.0 * abs(a);
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_a_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * le / (dis * pow(be, 4.0)) * (
			be * be * (3.0 * a * a + b * b + le * le) - 3.0 * pow(a, 3.0) * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.5 * b * le / dis;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		result = 0;
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * a * b / (dis * pow(be, 6.0)) * (
			be * be * (2.0 * pow(b, 4.0) + 5.0 * pow(b * le, 2.0) + 3.0 * le * le * (a * a + le * le))
			- 3.0 * pow(a, 3.0) * le * le * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * b / dis;
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 / (dis * pow(be, 6.0)) * (
			a * le * be * be * (-2.0 * pow(a * b, 2.0) + (a * a + b * b) * le * le + pow(le, 4.0))
			- pow(a, 4.0) * le * (-2.0 * b * b + le * le) * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = (0.5 / (le * le * le * dis)) * (a * le * le * dis * dis - pow(a, 4.0) * sqrt(le * le + pow(le, 4.0) / (a * a)) * asinh(abs(le) / a));
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_a_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * a / (dis * pow(be, 6.0)) * (
			be * be * (a * a * (b * b - 2.0 * le * le) + b * b * (b * b + le * le))
			- pow(a, 3.0) * (b * b - 2.0 * le * le) * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (0.5 / (b * b * b * dis)) * (a * b * b * dis * dis - pow(a, 4.0) * sqrt(b * b + pow(b, 4.0) / (a * a)) * asinh(abs(b) / a));
		}
		else {
			result = 0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_b_b_b(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 / (dis * pow(be, 8.0)) * (
			be * be * (3.0 * pow(a * a * le, 2.0) * (-4.0 * b * b + le * le)
				+ 2.0 * pow(be, 4.0) * (8.0 * pow(b, 4.0) + 12.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0))
				+ pow(a * be, 2.0) * (8.0 * pow(b, 4.0) + 12.0 * pow(b * le, 2.0) + 9.0 * pow(le, 4.0))
				)
			- 3.0 * pow(a, 5.0) * le * le * (-4.0 * b * b + le * le) * dis * be * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * (8.0 * pow(b, 4.0) + 12.0 * pow(b * le, 2.0) + 3.0 * pow(le, 4.0)) / pow(be, 3.0);
		}
		else {
			result = 0.75 * abs(le);
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (a * a + 2.0 * b * b) / dis;
		}
		else if (b == 0 && le != 0) {
			result = (0.125 / (pow(le, 4.0) * dis)) * (3.0 * pow(a * a * le, 2.0) + 9.0 * pow(a * le * le, 2.0) + 6.0 * pow(le, 6.0)
				- 3.0 * pow(a, 5.0) * sqrt(le * le + pow(le, 4.0) / (a * a)) * asinh(abs(le) / a));
		}
		else if (b == 0 && le == 0) {
			result = abs(a);
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_b_b_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * b * le / (dis * pow(be, 8.0)) * (
			be * be * (
				pow(a, 4.0) * (6.0 * b * b - 9.0 * le * le) + a * a * (2.0 * b * b - 3.0 * le * le) * be * be
				+ 2.0 * pow(be, 4.0) * (2.0 * b * b + 3.0 * le * le)
				)
			+ 3.0 * pow(a, 5.0) * dis * be * (-2.0 * b * b + 3.0 * le * le) * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * b * le * (2.0 * b * b + 3.0 * le * le) / pow(be, 3.0);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		result = 0;
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_b_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 / (dis * pow(be, 8.0)) * (
			be * be * (
				pow(a, 4.0) * (-2.0 * pow(b, 4.0) + 11.0 * pow(b * le, 2.0) - 2.0 * pow(le, 4.0))
				+ 2.0 * pow(be, 4.0) * (2.0 * pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
				+ a * a * pow(be, 2.0) * (2.0 * pow(b, 4.0) + 9.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0))
				)
			+ pow(a, 5.0) * dis * be * (2.0 * pow(b, 4.0) - 11.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0)) * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * pow(be, -3.0) * (2.0 * pow(b, 4.0) + 3.0 * pow(b * le, 2.0) + 2.0 * pow(le, 4.0));
		}
		else {
			result = 0.5 * abs(le);
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = (0.25 / (pow(b, 4.0) * dis)) * (-a * a * a * a * b * b + a * a * b * b * b * b + 2.0 * pow(b, 6.0) + pow(a, 5.0) * sqrt(b * b + pow(b, 4.0) / (a * a)) * asinh(abs(b) / a));
		}
		else if (b == 0 && le != 0) {
			result = (0.25 / (pow(le, 4.0) * dis)) * (-a * a * a * a * le * le + a * a * le * le * le * le + 2.0 * pow(le, 6.0) + pow(a, 5.0) * sqrt(le * le + pow(le, 4.0) / (a * a)) * asinh(abs(le) / a));
		}
		else if (b == 0 && le == 0) {
			result = abs(a) / 3.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Psi_I3_le_le_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.125 * b * le / (dis * pow(be, 8.0)) * (
			be * be * (
				2.0 * pow(be, 4.0) * (3.0 * b * b + 2.0 * le * le) + pow(a, 4.0) * (-9.0 * b * b + 6.0 * le * le)
				- a * a * (3.0 * pow(b, 4.0) + pow(b * le, 2.0) - 2.0 * pow(le, 4.0))
				)
			+ 3.0 * pow(a, 5.0) * dis * be * (3.0 * b * b - 2.0 * le * le) * asinh(be / a) / abs(a)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.25 * pow(be, -3.0) * b * le * (3.0 * b * b + 2.0 * le * le);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		result = 0;
	}

	return result;
}

/*--------------------------------------------------
	Components of the Phi_p linear one
--------------------------------------------------*/
// use Phi_u2_II_a instead

///////////////////////////////////////////////////////

 // Non- derivative form
// the Original form of Linear Phi potentials
double Linear_Polyhedral::Phi_Lin_H(double a, double b, double le)
{
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = a * a / (3.0 * b) * (
			a * abs(b) * asin(a * le / sqrt((a * a + b * b) * (b * b + le * le)))
			+ b * (-abs(a) * atan(le / b) + b * atanh(le / sqrt(a * a + b * b + le * le)))
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

double Linear_Polyhedral::Phi_Lin_pcos(double a, double b, double le)
{
	double result; double be = sqrt(b * b + le * le);
	double dis = sqrt(a * a + b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = a / 6.0 * (
			-a * a * le * atanh(be / dis) / be
			+ 0.5 * (a * a + b * b) * (log(le + dis) - log(dis - le))
			);
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 || le != 0) {
			result = a / 6.0 * (
				-a * a * le * atanh(be / dis) / be
				+ 0.5 * (a * a + b * b) * (log(le + dis) - log(dis - le))
				);
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Phi_Lin_psin(double a, double b, double le)
{
	double result; double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = a * b / 6.0 * (dis + a * a / be * atanh(be / dis));
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b == 0 && le != 0) {
			result = 0.0;
		}
		else if (b != 0 && le == 0) {
			result = a * b / 6.0 * (dis + a * a / be * atanh(be / dis));
		}
		else if (b == 0 && le == 0) {
			result = 0.0;
		}
	}
	return result;

}

double Linear_Polyhedral::Phi_I2_III(double a, double b, double le)
{
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	double result;
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * (
			-a * a * le * atanh(be / dis) / be + (a * a + b * b) * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = 0.5 * b * b * atanh(le / be);
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


double Linear_Polyhedral::Phi_I2_III_a(double a, double b, double le)
{
	double result = 0;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = a * (-le * atanh(be / dis) / be + atanh(le / dis));
	}
	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		result = 0;
	}



	return result;
}

double Linear_Polyhedral::Phi_I2_III_b(double a, double b, double le)
{
	double result = 0;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * (
			-le * dis / (be * be) + a * a * le * atanh(be / dis) / pow(be, 3.0) + 2.0 * atanh(le / dis)
			);
	}
	else if (a == 0) {
		if (b != 0) {
			result = -b * le / (2.0 * be) + b * atanh(le / be);
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

double Linear_Polyhedral::Phi_I2_III_le(double a, double b, double le)
{
	double result;
	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * b * dis * (
			be * be - a * a * be / dis * atanh(be / dis)
			) / pow(be, 4.0);
	}
	else if (a == 0) {
		// we known b, le cannot be zero at the same time
		if (b != 0 || le != 0) {
			result = b * b / be * 0.5;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.5 * dis * (
				b * b - a * a * be / dis * atanh(be / dis)
				) / (b * b);
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

double Linear_Polyhedral::Phi_I3_III(double a, double b, double le)
{
	double result = 0;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * (dis + a * a * atanh(be / dis) / be);
	}
	else if (a == 0) {
		result = 0.5 * b * be;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = 0.5 * (b * dis + a * a * atanh(b / dis));
		}
		else {
			result = 0.0;
		}
	}

	return result;
}

double Linear_Polyhedral::Phi_I3_III_a(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = a * b * atanh(be / dis) / be;
	}

	else if (a == 0) {
		result = 0.0;
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = a * b * atanh(be / dis) / be;
		}
		else if (b == 0 && le != 0) {
			result = 0;
		}
		else if (b == 0 && le == 0) {
			result = 0;
		}
	}
	else {
		result = 0.0;
	}
	return result;
}

double Linear_Polyhedral::Phi_I3_III_b(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);
	if (a != 0 && b != 0 && le != 0) {
		result = (be * dis * (2.0 * b * b + le * le) + a * a * le * le * atanh(be / dis)) / (2.0 * pow(be, 3.0));
	}
	else if (a == 0) {
		if (b != 0 && le != 0) {
			result = 0.5 * (2.0 * b * b + le * le) / be;
		}
		else if (b != 0 && le == 0) {
			result = b * b / abs(b);
		}
		else if (b == 0 && le != 0) {
			result = le * le * 0.5 / abs(le);
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		if (b != 0 && le == 0) {
			result = dis;
		}
		else if (b == 0 && le != 0) {
			result = 0.5 * dis + 0.5 * a * a * atanh(be / dis) / abs(le);
		}
		else if (b == 0 && le == 0) {
			result = abs(a);
		}
	}
	else {
		result = 0;
	}
	return result;
}

double Linear_Polyhedral::Phi_I3_III_le(double a, double b, double le)
{
	double result;	double dis = sqrt(a * a + b * b + le * le);
	double be = sqrt(b * b + le * le);

	if (a != 0 && b != 0 && le != 0) {
		result = 0.5 * b * le * (
			be * dis - a * a * atanh(be / dis)
			) / pow(be, 3.0);
	}
	else if (a == 0) { // knowing that b and le cannot be zero at the same time
		if (b != 0 || le != 0) {
			result = 0.5 * b * le / be;
		}
		else {
			result = 0.0;
		}
	}
	else if (a != 0) {
		result = 0;
	}
	else {
		result = 0;
	}
	return result;
}


// define the combined potentials

double Linear_Polyhedral::I2(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Phi_I2_III(a, b, lp) - Phi_I2_III(a, b, lm));
}

double Linear_Polyhedral::I3(double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	return (Phi_I3_III(a, b, lp) - Phi_I3_III(a, b, lm));
}



double Linear_Polyhedral::I2_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j] * (Phi_I2_III_a(a, b, lp) - Phi_I2_III_a(a, b, lm))
		- Vnorm[j] * (Phi_I2_III_b(a, b, lp) - Phi_I2_III_b(a, b, lm))
		- Vdir[j] * (Phi_I2_III_le(a, b, lp) - Phi_I2_III_le(a, b, lm))
		);

	return result;
}

double Linear_Polyhedral::I3_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j] * (Phi_I3_III_a(a, b, lp) - Phi_I3_III_a(a, b, lm))
		- Vnorm[j] * (Phi_I3_III_b(a, b, lp) - Phi_I3_III_b(a, b, lm))
		- Vdir[j] * (Phi_I3_III_le(a, b, lp) - Phi_I3_III_le(a, b, lm))
		);
	return result;
}


double Linear_Polyhedral::I5_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result;
	// using the numerical differentiation this time
	result = (-Svnorm[j]) * (Psi_I2_a(a, b, lp) - Psi_I2_a(a, b, lm));
	result += (-Vnorm[j]) * (Psi_I2_b(a, b, lp) - Psi_I2_b(a, b, lm));
	result += (-Vdir[j]) * (Psi_I2_le(a, b, lp) - Psi_I2_le(a, b, lm));

	return result;

}

double Linear_Polyhedral::I6_j(int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (-Svnorm[j]) * (Psi_I3_a(a, b, lp) - Psi_I3_a(a, b, lm));
	result += (-Vnorm[j]) * (Psi_I3_b(a, b, lp) - Psi_I3_b(a, b, lm));
	result += (-Vdir[j]) * (Psi_I3_le(a, b, lp) - Psi_I3_le(a, b, lm));

	return result;
}

double Linear_Polyhedral::I5_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	//double result = 0.0;
	double result = (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I2_a_a(a, b, lp) - Psi_I2_a_a(a, b, lm)) - Vnorm[k] * (Psi_I2_a_b(a, b, lp) - Psi_I2_a_b(a, b, lm)) - Vdir[k] * (Psi_I2_a_le(a, b, lp) - Psi_I2_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I2_a_b(a, b, lp) - Psi_I2_a_b(a, b, lm)) - Vnorm[k] * (Psi_I2_b_b(a, b, lp) - Psi_I2_b_b(a, b, lm)) - Vdir[k] * ((Psi_I2_b_le(a, b, lp) - Psi_I2_b_le(a, b, lm))))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I2_a_le(a, b, lp) - Psi_I2_a_le(a, b, lm)) - Vnorm[k] * ((Psi_I2_b_le(a, b, lp) - Psi_I2_b_le(a, b, lm))) - Vdir[k] * ((Psi_I2_le_le(a, b, lp) - Psi_I2_le_le(a, b, lm))))
		);

	return result;
}

double Linear_Polyhedral::I6_jk(int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = 0.0;
	result = (
		(-Svnorm[j]) * (-Svnorm[k] * (Psi_I3_a_a(a, b, lp) - Psi_I3_a_a(a, b, lm)) - Vnorm[k] * (Psi_I3_a_b(a, b, lp) - Psi_I3_a_b(a, b, lm)) - Vdir[k] * (Psi_I3_a_le(a, b, lp) - Psi_I3_a_le(a, b, lm)))
		+ (-Vnorm[j]) * (-Svnorm[k] * (Psi_I3_a_b(a, b, lp) - Psi_I3_a_b(a, b, lm)) - Vnorm[k] * (Psi_I3_b_b(a, b, lp) - Psi_I3_b_b(a, b, lm)) - Vdir[k] * ((Psi_I3_b_le(a, b, lp) - Psi_I3_b_le(a, b, lm))))
		+ (-Vdir[j]) * (-Svnorm[k] * (Psi_I3_a_le(a, b, lp) - Psi_I3_a_le(a, b, lm)) - Vnorm[k] * ((Psi_I3_b_le(a, b, lp) - Psi_I3_b_le(a, b, lm))) - Vdir[k] * ((Psi_I3_le_le(a, b, lp) - Psi_I3_le_le(a, b, lm))))
		);

	return result;
}



double Linear_Polyhedral::I5_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{
	double result = -(Psi_I2_a_a_a(a, b, lp) - Psi_I2_a_a_a(a, b, lm)) * Svnorm[i] * Svnorm[j] * Svnorm[k] - (Psi_I2_b_b_b(a, b, lp) - Psi_I2_b_b_b(a, b, lm)) * Vnorm[i] * Vnorm[j] * Vnorm[k]
		- (Psi_I2_le_le_le(a, b, lp) - Psi_I2_le_le_le(a, b, lm)) * Vdir[i] * Vdir[j] * Vdir[k] - (Psi_I2_a_a_b(a, b, lp) - Psi_I2_a_a_b(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vnorm[k] + Svnorm[i] * Vnorm[j] * Svnorm[k] + Vnorm[i] * Svnorm[j] * Svnorm[k])
		- (Psi_I2_a_a_le(a, b, lp) - Psi_I2_a_a_le(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vdir[k] + Svnorm[i] * Vdir[j] * Svnorm[k] + Vdir[i] * Svnorm[j] * Svnorm[k]) - (Psi_I2_a_b_b(a, b, lp) - Psi_I2_a_b_b(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vnorm[k] + Svnorm[j] * Vnorm[i] * Vnorm[k] + Svnorm[k] * Vnorm[j] * Vnorm[i])
		- (Psi_I2_b_b_le(a, b, lp) - Psi_I2_b_b_le(a, b, lm)) * (Vdir[i] * Vnorm[j] * Vnorm[k] + Vdir[j] * Vnorm[i] * Vnorm[k] + Vdir[k] * Vnorm[i] * Vnorm[j]) - (Psi_I2_a_le_le(a, b, lp) - Psi_I2_a_le_le(a, b, lm)) * (Svnorm[i] * Vdir[j] * Vdir[k] + Svnorm[j] * Vdir[i] * Vdir[k] + Svnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I2_b_le_le(a, b, lp) - Psi_I2_b_le_le(a, b, lm)) * (Vnorm[i] * Vdir[j] * Vdir[k] + Vnorm[j] * Vdir[i] * Vdir[k] + Vnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I2_a_b_le(a, b, lp) - Psi_I2_a_b_le(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vdir[k] + Svnorm[i] * Vnorm[k] * Vdir[j] + Svnorm[j] * Vnorm[i] * Vdir[k] + Svnorm[k] * Vnorm[i] * Vdir[j] + Svnorm[j] * Vnorm[k] * Vdir[i] + Svnorm[k] * Vnorm[j] * Vdir[i]);

	return result;
}


double Linear_Polyhedral::I6_jkl(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir)
{

	double result = -(Psi_I3_a_a_a(a, b, lp) - Psi_I3_a_a_a(a, b, lm)) * Svnorm[i] * Svnorm[j] * Svnorm[k] - (Psi_I3_b_b_b(a, b, lp) - Psi_I3_b_b_b(a, b, lm)) * Vnorm[i] * Vnorm[j] * Vnorm[k]
		- (Psi_I3_le_le_le(a, b, lp) - Psi_I3_le_le_le(a, b, lm)) * Vdir[i] * Vdir[j] * Vdir[k] - (Psi_I3_a_a_b(a, b, lp) - Psi_I3_a_a_b(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vnorm[k] + Svnorm[i] * Vnorm[j] * Svnorm[k] + Vnorm[i] * Svnorm[j] * Svnorm[k])
		- (Psi_I3_a_a_le(a, b, lp) - Psi_I3_a_a_le(a, b, lm)) * (Svnorm[i] * Svnorm[j] * Vdir[k] + Svnorm[i] * Vdir[j] * Svnorm[k] + Vdir[i] * Svnorm[j] * Svnorm[k]) - (Psi_I3_a_b_b(a, b, lp) - Psi_I3_a_b_b(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vnorm[k] + Svnorm[j] * Vnorm[i] * Vnorm[k] + Svnorm[k] * Vnorm[j] * Vnorm[i])
		- (Psi_I3_b_b_le(a, b, lp) - Psi_I3_b_b_le(a, b, lm)) * (Vdir[i] * Vnorm[j] * Vnorm[k] + Vdir[j] * Vnorm[i] * Vnorm[k] + Vdir[k] * Vnorm[i] * Vnorm[j]) - (Psi_I3_a_le_le(a, b, lp) - Psi_I3_a_le_le(a, b, lm)) * (Svnorm[i] * Vdir[j] * Vdir[k] + Svnorm[j] * Vdir[i] * Vdir[k] + Svnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I3_b_le_le(a, b, lp) - Psi_I3_b_le_le(a, b, lm)) * (Vnorm[i] * Vdir[j] * Vdir[k] + Vnorm[j] * Vdir[i] * Vdir[k] + Vnorm[k] * Vdir[i] * Vdir[j])
		- (Psi_I3_a_b_le(a, b, lp) - Psi_I3_a_b_le(a, b, lm)) * (Svnorm[i] * Vnorm[j] * Vdir[k] + Svnorm[i] * Vnorm[k] * Vdir[j] + Svnorm[j] * Vnorm[i] * Vdir[k] + Svnorm[k] * Vnorm[i] * Vdir[j] + Svnorm[j] * Vnorm[k] * Vdir[i] + Svnorm[k] * Vnorm[j] * Vdir[i]);

	return result;
}


// complete set of potentials

double Linear_Polyhedral::PHI_p(int p, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = x[p] * U_PL.PHI(a, b, lm, lp);

	result += Svnorm[p] * (Phi_Lin_H(a, b, lp) - Phi_Lin_H(a, b, lm))
		+ Vnorm[p] * (Phi_Lin_pcos(a, b, lp) - Phi_Lin_pcos(a, b, lm))
		+ Vdir[p] * (Phi_Lin_psin(a, b, lp) - Phi_Lin_psin(a, b, lm));

	return result;
}


double Linear_Polyhedral::PHI_p_i(int p, int i, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	result = (a * Svnorm[p] + x[p]) * U_PL.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vnorm[p] * I2(a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vdir[p] * I3(a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[i][p] * U_PL.PHI(a, b, lm, lp);

	return result;
}

double Linear_Polyhedral::PHI_p_ij(int p, int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	result = (-Svnorm[j] * Svnorm[p] + d[j][p]) * U_PL.PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (a * Svnorm[p] + x[p]) * U_PL.PHI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vnorm[p] * I2_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vdir[p] * I3_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[i][p] * U_PL.PHI_i(j, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	//+ d[j][p] * PHI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	//+ x[p] * PHI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	return result;
}

double Linear_Polyhedral::PSI_p_ij(int p, int i, int j, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	result = (-Svnorm[j] * Svnorm[p] + d[j][p]) * U_PL.PSI_i(i, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (a * Svnorm[p] + x[p]) * U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (Svnorm[i]) * Vnorm[p] * I5_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- (Svnorm[i]) * Vdir[p] * I6_j(j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[i][p] * U_PL.PSI_i(j, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	return result;
}

double Linear_Polyhedral::PSI_p_ijk(int p, int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0.0;
	result = (-Svnorm[i] * Svnorm[p] + d[i][p]) * U_PL.PSI_ij(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Svnorm[p] * U_PL.PSI_ij(k, j, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (Svnorm[p] * a + x[p]) * U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vnorm[p] * I5_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vdir[p] * I6_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[p][j] * U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[p][k] * U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	return result;
}

double Linear_Polyhedral::PSI_p_ijkl(int p, int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x)
{
	double result = 0;
	result = (-Svnorm[i] * Svnorm[p] + d[i][p]) * U_PL.PSI_ijk(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Svnorm[p] * U_PL.PSI_ijk(k, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Svnorm[p] * U_PL.PSI_ijk(l, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ (Svnorm[p] * a + x[p]) * U_PL.PSI_ijkl(i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vnorm[p] * I5_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		- Svnorm[i] * Vdir[p] * I6_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		//+ d[i][p] * PSI_ijk(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[p][j] * U_PL.PSI_ijk(i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[p][k] * U_PL.PSI_ijk(i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir)
		+ d[p][l] * U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	//+ x[p] * PSI_ijkl(i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	return result;
}

void Linear_Polyhedral::PSI_p_ijk_combine(int i, int j, int k, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* result, double& ij, \
	double& ik, double& jk, double& kj, double& ijk, double& jk_5, double& jk_6)
{
	ij = U_PL.PSI_ij(i, j, a, b, lm, lp, Svnorm, Vnorm, Vdir);  ik = U_PL.PSI_ij(i, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);  jk = U_PL.PSI_ij(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	kj = U_PL.PSI_ij(k, j, a, b, lm, lp, Svnorm, Vnorm, Vdir);  ijk = U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	jk_5 = I5_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);  jk_6 = I6_jk(j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	for (int p = 0; p < 3; p++) {
		result[p] = (-Svnorm[i] * Svnorm[p] + d[i][p]) * jk - Svnorm[i] * Svnorm[p] * kj + (Svnorm[p] * a + x[p]) * ijk
			- Svnorm[i] * Vnorm[p] * jk_5 - Svnorm[i] * Vdir[p] * jk_6 + d[p][j] * ik + d[p][k] * ij;
	}

}


void Linear_Polyhedral::PSI_p_ijkl_combine(int i, int j, int k, int l, double a, double b, double lm, double lp, double* Svnorm, double* Vnorm, double* Vdir, double* x, double* result, double& ijk, double& ikl, double& ijl, \
	double& jkl, double& kjl, double& ljk, double& ijkl, double& jkl_5, double& jkl_6)
{
	ikl = U_PL.PSI_ijk(i, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);  ijl = U_PL.PSI_ijk(i, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);  ijk = U_PL.PSI_ijk(i, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	jkl = U_PL.PSI_ijk(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);  kjl = U_PL.PSI_ijk(k, j, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	ljk = U_PL.PSI_ijk(l, j, k, a, b, lm, lp, Svnorm, Vnorm, Vdir);  ijkl = U_PL.PSI_ijkl(i, j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);
	jkl_5 = I5_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);  jkl_6 = I6_jkl(j, k, l, a, b, lm, lp, Svnorm, Vnorm, Vdir);

	for (int p = 0; p < 3; p++) {
		result[p] = (-Svnorm[i] * Svnorm[p] + d[i][p]) * jkl - Svnorm[i] * Svnorm[p] * kjl - Svnorm[i] * Svnorm[p] * ljk
			+ (Svnorm[p] * a + x[p]) * ijkl - Svnorm[i] * Vnorm[p] * jkl_5 - Svnorm[i] * Vdir[p] * jkl_6 + d[p][j] * ikl + d[p][k] * ijl + d[p][l] * ijk;
	}

}