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
#include "integrator.h"

class integratorElastic3D : public integrator {
public:
	void addFieldToBEM(Config& config_);
	void addFluxEquivalentInclusion(Config& config_);
	void addFluxEquivalentFirstOrderInclusion(Config& config_);
	void addFluxEquivalentSecondOrderInclusion(Config& config_);
};

//void integratorElastic3D::addFieldToBEM(Config& config_) {};
//void integratorElastic3D::addFluxEquivalentInclusion(Config& config_) {};
//void integratorElastic3D::addFluxEquivalentFirstOrderInclusion(Config& config_) {};
//void integratorElastic3D::addFluxEquivalentSecondOrderInclusion(Config& config_) {};