#ifndef DiagnosticPhase2DMomMom_H
#define DiagnosticPhase2DMomMom_H

#include "DiagnosticPhase2D.h"

class SmileiMPI;

//! class projection on the plane momentum-momentum
class DiagnosticPhase2DMomMom : public DiagnosticPhase2D {
    
public:
    //! the component of the 1st momentum (0,1,2) and the component of the 2nd momentum (0,1,2) are chosen at construction
    DiagnosticPhase2DMomMom(phaseStructure phaseStruct, const unsigned int directionMomentum1, const unsigned int directionMomentum2);	

	void run(partStruct&);
private:
    //! component of the momentum for the first axe
    const unsigned int my_dirMom1;
    //! component of the momentum for the second axe
    const unsigned int my_dirMom2;
};

#endif
