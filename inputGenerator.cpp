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
		default:
			cout << "You have not specified the types of PDEs, check CALCULATIONTYPE.txt" << endl;
	}


	return NULL;
}
