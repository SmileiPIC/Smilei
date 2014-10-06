#include "DiagnosticPhasePosLor.h"

using namespace std;

DiagnosticPhasePosLor::DiagnosticPhasePosLor(phaseStructure phaseStruct, const unsigned int directionPosition) : 
DiagnosticPhase(phaseStruct), 
my_dirPos(directionPosition) {
    
	if (phaseStruct.pos_num.size() >0 && phaseStruct.lor_num.size() >0) {
		my_data.allocateDims(phaseStruct.pos_num[0],phaseStruct.lor_num[0]);
	} else {
		ERROR("must define pos_ and lor_ stuff");
	}
	//!\todo add more checks here (TV MC MG)
	firstmin = phaseStruct.pos_min[0];
	firstmax = phaseStruct.pos_max[0];
	firstnum = phaseStruct.pos_num[0];
	secondmin = phaseStruct.lor_min[0];
	secondmax = phaseStruct.lor_max[0];
	secondnum = phaseStruct.lor_num[0];    
}

void DiagnosticPhasePosLor::run(partStruct& my_part) {
 	if (my_part.pos[my_dirPos] > (2*M_PI)*firstmin && my_part.pos[my_dirPos] < (2*M_PI)*firstmax) {
        double lor_fact=my_part.lor_fact();
        if(lor_fact > secondmin && lor_fact < secondmax) {
            
            //!\todo check if useful to have projector here
            int i = firstnum*(my_part.pos[my_dirPos]-(2*M_PI)*firstmin)/((2*M_PI)*(firstmax-firstmin));
            int j = secondnum*(lor_fact-secondmin)/(secondmax-secondmin);
            my_data(i,j)+=my_part.weight;		
        }
    }
}
