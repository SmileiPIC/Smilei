#ifndef DiagnosticPhase2DxPx_H
#define DiagnosticPhase2DxPx_H

#include "DiagnosticPhase2D.h"

class SmileiMPI;

//! class projection on the plane x-Px
class DiagnosticPhase2DxPx : public DiagnosticPhase2D {

public:
    DiagnosticPhase2DxPx(phaseStructure phaseStruct);	
	void doSomething(partStruct&);
};
#endif
