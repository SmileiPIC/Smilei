#include "DiagnosticPhasePosMom.h"

using namespace std;

DiagnosticPhasePosMom::DiagnosticPhasePosMom(phaseStructure phaseStruct, const unsigned int directionPosition, const unsigned int directionMomentum) : 
DiagnosticPhase(phaseStruct), 
my_dirPos(directionPosition), 
my_dirMom(directionMomentum) {
    
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

void DiagnosticPhasePosMom::run(partStruct& my_part) {
	if (my_part.pos[my_dirPos] > (2*M_PI)*firstmin && my_part.pos[my_dirPos] < (2*M_PI)*firstmax && my_part.mom[my_dirMom] > secondmin && my_part.mom[my_dirMom] < secondmax) {
		//!\todo check if useful to have projector here
		int i = firstnum*(my_part.pos[my_dirPos]-(2*M_PI)*firstmin)/((2*M_PI)*(firstmax-firstmin));
		int j = secondnum*(my_part.mom[my_dirMom]-secondmin)/(secondmax-secondmin);
		my_data(i,j)+=my_part.weight;		
	}
}
