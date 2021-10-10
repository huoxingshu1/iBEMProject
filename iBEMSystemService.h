#pragma once
#include "Config.h"
#include "BEMBuilder.h"
#include "integrator.h"

class iBEMSystemService {
public:
	Config& config;
	BEMBuilder& bembuilder;
	integrator& Integrator;
	iBEMSystemService(Config& config_, BEMBuilder& bembuilder_, integrator& integrator_) : config{ config_ }, bembuilder{ bembuilder_ }, Integrator{ integrator_ }{};
	void run();
	void tradition_BEM();             // build up traditioanl boundary integral equations
	void addFieldToBEM();        // build up effects from inclusion to boundary integral equations
	void addFluxEquivalentBEM(); // build up constant term of flux equivalent equations from boundary integral
	void addFluxEquivalentFirstOrderBEM();  // build up first order term of flux equivalent equations from boundary integral
	void addFluxEquivalentSecondOrderBEM();  // build up second order term of flux equivalent equations from boundary integral
	void addFluxEquivalentInclusion();  // build up  stress quivalent equations for includsions (equivalent includsion method)


};

void iBEMSystemService::run() {
	this->tradition_BEM();
	this->addFieldToBEM();
	this->addFluxEquivalentBEM();
	this->addFluxEquivalentFirstOrderBEM();
	this->addFluxEquivalentSecondOrderBEM();
	this->addFluxEquivalentInclusion();

}


void iBEMSystemService::tradition_BEM() {
	bembuilder.addBEM(config);
}

void iBEMSystemService::addFieldToBEM() {
	Integrator.addFieldToBEM(config);

}

void iBEMSystemService::addFluxEquivalentBEM() {
	bembuilder.addFluxEquivalentBEM(config);

}

void iBEMSystemService::addFluxEquivalentFirstOrderBEM() {
	bembuilder.addFluxEquivalentFirstOrderBEM(config);

}

void iBEMSystemService::addFluxEquivalentSecondOrderBEM() {
	bembuilder.addFluxEquivalentSecondOrderBEM(config);

}

void iBEMSystemService::addFluxEquivalentInclusion() {
	Integrator.addFluxEquivalentInclusion(config);
	Integrator.addFluxEquivalentFirstOrderInclusion(config);
	Integrator.addFluxEquivalentSecondOrderInclusion(config);

}