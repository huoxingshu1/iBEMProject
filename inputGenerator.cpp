#include "inputGenerator.h"

shared_ptr<Config> inputGenerator(int ID) {

	if (ID == 1) {
		return make_shared<configElastic3D>();
	}

	return NULL;
}