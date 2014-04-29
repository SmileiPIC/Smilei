#ifndef DiagnosticPhase1D_H
#define DiagnosticPhase1D_H

#include "DiagnosticPhase.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

//! this class holds all the phase projections that can be represented as 1d matrix (aka vector)
class DiagnosticPhase1D : public DiagnosticPhase {
//! \todo check is is useful and ceck if we would need a DiagnosticPhase3D
public:
    DiagnosticPhase1D(phaseStructure phaseStruct): DiagnosticPhase(phaseStruct){};

};
#endif
