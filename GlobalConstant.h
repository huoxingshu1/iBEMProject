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

#pragma once
const double Pi = 3.141592653589793;
const int d[3][3] = { 1,0,0,0,1,0,0,0,1 };
const int map[3][3] = { 1,2,2,0,1,2,0,0,1 };
const int surface[4][3] = { {0,2,1},{1,2,3},{0,3,2},{0,1,3} };
/*
    Add Gaussian integral points to use subdomain collcation method for 2D arbitrary cases
*/
const double Shape_L1_2D[6] = { 0.8168475729, 0.0915762135, 0.0915762135, \
    0.1081030181, 0.4459484909,  0.4459484909 };
const double Shape_L2_2D[6] = { 0.0915762135, 0.8168475729, 0.0915762135, \
    0.4459484909, 0.1081030181, 0.4459484909 };
const double Shape_L3_2D[6] = { 0.0915762135, 0.0915762135, 0.8168475729, \
    0.4459484909, 0.4459484909, 0.1081030181 };
const double Weight_2D[6] = { 0.1099517436, 0.1099517436, 0.1099517436,\
    0.2233815896, 0.2233815896, 0.2233815896 };

/*
    Add Gaussian integral points to use subdomain collocation method for 3D arbitrary cases
*/
const double Shape_L1_3D[10] = { 0.7784952948213300,  0.0738349017262234,  0.0738349017262234,  0.0738349017262234, \
0.4062443438840510, 0.4062443438840510, 0.406244343884051, 0.0937556561159491, 0.0937556561159491, 0.0937556561159491 };

const double Shape_L2_3D[10] = { 0.0738349017262234,  0.7784952948213300, 0.0738349017262234, 0.0738349017262234, \
0.4062443438840510, 0.0937556561159491, 0.0937556561159491, 0.4062443438840510, 0.4062443438840510, 0.0937556561159491 };

const double Shape_L3_3D[10] = { 0.0738349017262234, 0.0738349017262234, 0.7784952948213300, 0.0738349017262234, 0.0937556561159491, \
0.4062443438840510, 0.0937556561159491, 0.4062443438840510, 0.0937556561159491, 0.4062443438840510 };

const double Shape_L4_3D[10] = { 0.0738349017262234, 0.0738349017262234, 0.0738349017262234, 0.7784952948213300, 0.0937556561159491, \
0.0937556561159491, 0.4062443438840510, 0.0937556561159491, 0.4062443438840510, 0.4062443438840510 };

const double Weight_3D[10] = { 0.0476331348432089,  0.0476331348432089,  0.0476331348432089,  0.0476331348432089 ,\
0.1349112434378610, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610 };