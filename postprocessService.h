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
#include "Config.h"
#include "postProcessor.h"

class postprocessService {
public:
	Config& config;
	postProcessor& PostProcessor;
	postprocessService(Config& config_, postProcessor& postProcessor_) : config{ config_ }, PostProcessor{ postProcessor_ } {};
	void run();
	void getField();  // compute the field value for interested point;
	void getFlux();   // compute the flux value for interested point;

};

void postprocessService::run() {
	this->getField();
	this->getFlux();
}

void postprocessService::getField() {
	PostProcessor.getField(config);
}


void postprocessService::getFlux() {
	PostProcessor.getFlux(config);
}