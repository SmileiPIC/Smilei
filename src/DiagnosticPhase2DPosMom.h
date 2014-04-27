#ifndef DiagnosticPhase2DPosMom_H
#define DiagnosticPhase2DPosMom_H

#include "DiagnosticPhase2D.h"

class SmileiMPI;

//! class projection on the plane position-momentum
//! the component of the position (0,1,2) and the component of the momentum (0,1,2) are chosen at construction
class DiagnosticPhase2DPosMom : public DiagnosticPhase2D {
    
public:
    //! the component of the position (0,1,2) and the component of the momentum (0,1,2) are passed by directionPosition and directionMomentum
    DiagnosticPhase2DPosMom(phaseStructure phaseStruct, const unsigned int directionPosition, const unsigned int directionMomentum);	

	void doSomething(partStruct&);
private:
    //! component of the position for the first axe
    const unsigned int my_dirPos;
    //! component of the momentum for the first axe
    const unsigned int my_dirMom;
};

#endif
