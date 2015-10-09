#include "DiagnosticPhasePosLor.h"

using namespace std;

DiagnosticPhasePosLor::DiagnosticPhasePosLor(Params &params, unsigned int n_phase, const unsigned int directionPosition) :
DiagnosticPhase(params, n_phase),
my_dirPos(directionPosition) {
}

void DiagnosticPhasePosLor::run(partStruct& my_part) {
 	if (my_part.pos[my_dirPos] > firstmin && my_part.pos[my_dirPos] < firstmax) {
        double lor_fact=my_part.lor_fact();
        if(lor_fact > secondmin && lor_fact < secondmax) {
            
            //!\todo check if useful to have projector here
            int i = firstnum * (my_part.pos[my_dirPos] - firstmin) / (firstmax - firstmin);
            int j = secondnum*(lor_fact-secondmin)/(secondmax-secondmin);
            my_data(i,j)+=my_part.weight;		
        }
    }
}
