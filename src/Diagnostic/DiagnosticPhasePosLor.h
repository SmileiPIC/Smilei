#ifndef DiagnosticPhasePosLor_H
#define DiagnosticPhasePosLor_H

#include "DiagnosticPhase.h"

//! class projection on the plane position-Lorentz_factor
//! the component of the position (0,1,2) is chosen at construction
class DiagnosticPhasePosLor : public DiagnosticPhase {
    
public:
    //! the component of the position (0,1,2) is passed by directionPosition
    DiagnosticPhasePosLor(phaseStructure phaseStruct, const unsigned int directionPosition);	

	void run(partStruct&);
private:
    //! component of the position for the first axe
    const unsigned int my_dirPos;
};

#endif
