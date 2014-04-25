#ifndef DiagnosticPhase2D_H
#define DiagnosticPhase2D_H

#include "DiagnosticPhase.h"
#include "Field2D.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

class DiagnosticPhase2D : public DiagnosticPhase {

public:

    DiagnosticPhase2D(phaseStructure phaseStruct, hid_t gid);
    ~DiagnosticPhase2D();

	Field2D my_data;
	
	
};
#endif
