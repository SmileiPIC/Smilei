#ifndef DiagnosticPhase1D_H
#define DiagnosticPhase1D_H

#include "DiagnosticPhase.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

class DiagnosticPhase1D : public DiagnosticPhase {

public:

    DiagnosticPhase1D(phaseStructure phaseStruct, hid_t gid): DiagnosticPhase(phaseStruct, gid){};
    ~DiagnosticPhase1D(){};

};
#endif
