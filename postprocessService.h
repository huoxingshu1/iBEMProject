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