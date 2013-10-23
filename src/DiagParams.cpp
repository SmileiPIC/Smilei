#include "DiagParams.h"
#include "Tools.h"
#include <cmath>

using namespace std;

DiagParams::DiagParams() {
}

void DiagParams::parseInputData(InputData &ifile, PicParams& params) {
	ifile.extract("every",scalar_every,"diagnostic scalar");
	ifile.extract("every",map_every,"diagnostic map");
	
}

	/*******************************************************************************************************************
	 caclulate useful parameters
	 ******************************************************************************************************************/
void DiagParams::compute() {
}

void DiagParams::print() {

}

