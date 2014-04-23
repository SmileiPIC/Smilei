#ifndef DiagnosticPhase1DxPx_H
#define DiagnosticPhase1DxPx_H

#include "DiagnosticPhase1D.h"


class DiagnosticPhase1DxPx : public DiagnosticPhase1D {

public:

    DiagnosticPhase1DxPx(phaseStructure phaseStruct, hid_t gid);
    ~DiagnosticPhase1DxPx();

};
#endif
