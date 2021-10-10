#pragma once
#include "Config.h"
class postProcessor {
public:
	virtual void getFlux(Config& config_) {}; // compute the field value for interested point;
	virtual void getField(Config& config_) {}; // compute the flux value for interested point;
};
