#ifndef DiagnosticPhase2DPosMom_H
#define DiagnosticPhase2DPosMom_H

#include "DiagnosticPhase2D.h"

class SmileiMPI;

//! class projection on the plane position-momentum
class DiagnosticPhase2DPosMom : public DiagnosticPhase2D {
    
public:
    //! the component of the position (0,1,2) and the component of the momentum (0,1,2) are chosen at construction
    DiagnosticPhase2DPosMom(phaseStructure phaseStruct, const unsigned int directionPosition, const unsigned int directionMomentum);	

	void doSomething(partStruct&);
private:
    //! component of the position for the first axe
    const unsigned int my_dirPos;
    //! component of the momentum for the second axe
    const unsigned int my_dirMom;
};

#endif
