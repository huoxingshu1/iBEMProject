#pragma once
#include "postProcessor.h"
#include "configElastic3D.h"
#include "Config.h"

class postProcessorElastic3D : public postProcessor {
	void getFlux( Config& config_);
	void getField( Config& config_);

};




