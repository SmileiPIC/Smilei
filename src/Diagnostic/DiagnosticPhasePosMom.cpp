#include "DiagnosticPhasePosMom.h"

using namespace std;

DiagnosticPhasePosMom::DiagnosticPhasePosMom(Params &params, unsigned int n_phase, const unsigned int directionPosition, const unsigned int directionMomentum) :
DiagnosticPhase(params, n_phase),
my_dirPos(directionPosition), 
my_dirMom(directionMomentum) {
}

void DiagnosticPhasePosMom::run(partStruct& my_part) {
	if (my_part.pos[my_dirPos] > firstmin && my_part.pos[my_dirPos] < firstmax && my_part.mom[my_dirMom] > secondmin && my_part.mom[my_dirMom] < secondmax) {
		//!\todo check if useful to have projector here
		int i = firstnum*(my_part.pos[my_dirPos]-firstmin)/(firstmax-firstmin);
		int j = secondnum*(my_part.mom[my_dirMom]-secondmin)/(secondmax-secondmin);
		my_data(i,j)+=my_part.weight;		
	}
}
