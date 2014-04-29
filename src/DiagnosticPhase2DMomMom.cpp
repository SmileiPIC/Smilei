#include "DiagnosticPhase2DMomMom.h"

using namespace std;

DiagnosticPhase2DMomMom::DiagnosticPhase2DMomMom(phaseStructure phaseStruct, const unsigned int directionMomentum1, const unsigned int directionMomentum2) : 
DiagnosticPhase2D(phaseStruct), 
my_dirMom1(directionMomentum1), 
my_dirMom2(directionMomentum2) {
    
	if (phaseStruct.mom_num.size() > 1) {
		my_data.allocateDims(phaseStruct.mom_num[0],phaseStruct.mom_num[1]);
	} else {
		ERROR("must define mom_ stuff (2 vals each)");
	}
	//!\todo add more checks here (TV MC MG)
	firstmin = phaseStruct.mom_min[0];
	firstmax = phaseStruct.mom_max[0];
	firstnum = phaseStruct.mom_num[0];
	secondmin = phaseStruct.mom_min[1];
	secondmax = phaseStruct.mom_max[1];
	secondnum = phaseStruct.mom_num[1];    
}

void DiagnosticPhase2DMomMom::doSomething(partStruct& my_part) {
	if (my_part.mom[my_dirMom1] > firstmin && my_part.mom[my_dirMom1] < firstmax && my_part.mom[my_dirMom2] > secondmin && my_part.mom[my_dirMom2] < secondmax) {
		//!\todo check if useful to have projector here
		int i = firstnum*(my_part.mom[my_dirMom1]-firstmin)/(firstmax-firstmin);
		int j = secondnum*(my_part.mom[my_dirMom2]-secondmin)/(secondmax-secondmin);
		my_data(i,j)+=my_part.weight;		
	}
}
