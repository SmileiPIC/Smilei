#ifndef DiagnosticPhase2DxPx_H
#define DiagnosticPhase2DxPx_H

#include "DiagnosticPhase2D.h"

class DiagnosticPhase2DxPx : public DiagnosticPhase2D {

public:

    DiagnosticPhase2DxPx(phaseStructure phaseStruct, hid_t gid);
    ~DiagnosticPhase2DxPx();
	
	void doSomething(partStruct);

};
#endif
