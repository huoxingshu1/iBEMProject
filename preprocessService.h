#pragma once
#include "Config.h"
class preprocessService {
public:
	Config& config;
	preprocessService(Config& config_) : config{ config_ } {};
	void run();

};

void preprocessService::run() {
	config.setUp();                        // set up configuration
}