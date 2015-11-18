#include "DiagnosticPhaseMomMom.h"

using namespace std;

DiagnosticPhaseMomMom::DiagnosticPhaseMomMom(Params &params, unsigned int n_phase, const unsigned int directionMomentum1, const unsigned int directionMomentum2) :
DiagnosticPhase(params, n_phase),
my_dirMom1(directionMomentum1), 
my_dirMom2(directionMomentum2) {
}

void DiagnosticPhaseMomMom::run(partStruct& my_part) {
	if (my_part.mom[my_dirMom1] > firstmin && my_part.mom[my_dirMom1] < firstmax && my_part.mom[my_dirMom2] > secondmin && my_part.mom[my_dirMom2] < secondmax) {
		//!\todo check if useful to have projector here
		int i = firstnum*(my_part.mom[my_dirMom1]-firstmin)/(firstmax-firstmin);
		int j = secondnum*(my_part.mom[my_dirMom2]-secondmin)/(secondmax-secondmin);
		my_data(i,j)+=my_part.weight;		
	}
}
