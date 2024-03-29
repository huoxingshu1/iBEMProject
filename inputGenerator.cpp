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

#include "inputGenerator.h"

// this source file specifies the types of PDEs to solve
shared_ptr<Config> inputGenerator(int ID) {
	
	switch (ID){
		case 1:
			return make_shared<configElastic3D>();
			break;
		case 2:
			return make_shared<configElastic2D>();
			break;
		
		case 3:
			return make_shared<configElastic3D_arbitrary>();

		case 4: 
			return make_shared<configElastic2D_arbitrary>();
			break;
		case 5:
			return make_shared<configPotential3D_elliptical>();
			break;
		default:
			cout << "You have not specified the types of PDEs, check CALCULATIONTYPE.txt" << endl;
	}


	return NULL;
}
