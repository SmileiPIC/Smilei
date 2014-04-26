#include "DiagnosticPhase2DxPx.h"

using namespace std;

DiagnosticPhase2DxPx::DiagnosticPhase2DxPx(phaseStructure phaseStruct) : DiagnosticPhase2D(phaseStruct) {

	if (phaseStruct.pos_num.size() >0 && phaseStruct.mom_num.size() >0) {
		my_data.allocateDims(phaseStruct.pos_num[0],phaseStruct.mom_num[0]);
	} else {
		ERROR("must define pos_ and mom_ stuff");
	}
	//!\todo add more checks here (TV MC MG)
	firstmin = phaseStruct.pos_min[0];
	firstmax = phaseStruct.pos_max[0];
	firstnum = phaseStruct.pos_num[0];
	secondmin = phaseStruct.mom_min[0];
	secondmax = phaseStruct.mom_max[0];
	secondnum = phaseStruct.mom_num[0];
		
}

void DiagnosticPhase2DxPx::doSomething(partStruct& my_part) {
	if (my_part.pos[0] > firstmin && my_part.pos[0] < firstmax && my_part.mom[0] > secondmin && my_part.mom[0] < secondmax) {
		//!\todo check if useful to have projector here
		int i = firstnum*(my_part.pos[0]-firstmin)/(firstmax-firstmin);
		int j = secondnum*(my_part.mom[0]-secondmin)/(secondmax-secondmin);
		my_data(i,j)+=my_part.weight;		
	}
}

