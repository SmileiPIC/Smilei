#ifndef DiagnosticPhase2DxPx_H
#define DiagnosticPhase2DxPx_H

#include "DiagnosticPhase2D.h"

class SmileiMPI;

//! class projection on the plane x-Px
class DiagnosticPhase2DxP : public DiagnosticPhase2D {
    
public:
    DiagnosticPhase2DxP(phaseStructure phaseStruct, const unsigned int direction);	
	void doSomething(partStruct&);
private:
    const unsigned int my_dir;
};

#endif
