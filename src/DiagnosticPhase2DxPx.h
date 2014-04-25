#ifndef DiagnosticPhase2DxPx_H
#define DiagnosticPhase2DxPx_H

#include "DiagnosticPhase2D.h"

class DiagnosticPhase2DxPx : public DiagnosticPhase2D {

public:

    DiagnosticPhase2DxPx(phaseStructure phaseStruct, hid_t gid);
    ~DiagnosticPhase2DxPx();
	
	void doSomething(short charge, double weight, double mom_x, double mom_y, double mom_z, double pos_x, double pos_y = 0, double pos_z=0);

};
#endif
