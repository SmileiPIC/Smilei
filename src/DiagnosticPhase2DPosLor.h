#ifndef DiagnosticPhase2DPosLor_H
#define DiagnosticPhase2DPosLor_H

#include "DiagnosticPhase2D.h"

class SmileiMPI;

//! class projection on the plane position-Lorentz_factor
//! the component of the position (0,1,2) is chosen at construction
class DiagnosticPhase2DPosLor : public DiagnosticPhase2D {
    
public:
    //! the component of the position (0,1,2) is passed by directionPosition
    DiagnosticPhase2DPosLor(phaseStructure phaseStruct, const unsigned int directionPosition);	

	void run(partStruct&);
private:
    //! component of the position for the first axe
    const unsigned int my_dirPos;
};

#endif
