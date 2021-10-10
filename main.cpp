#include<string>
#include "ibemRunner.h"
using namespace std;



int main() {

	//define the target service!!!!! {"preprocessService", "iBEMSystemService", "systemSolverService", "postprocessService" }

	string mainService = "iBEMSystemService";
	ibemRunner runner{};
	runner.run(mainService);
	
	return 0;


}


