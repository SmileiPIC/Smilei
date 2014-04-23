#include "DiagnosticPhase1DxPx.h"

using namespace std;

DiagnosticPhase1DxPx::~DiagnosticPhase1DxPx() {
}

DiagnosticPhase1DxPx::DiagnosticPhase1DxPx(phaseStructure phaseStruct, hid_t gid) : DiagnosticPhase1D(phaseStruct, gid) {
	DEBUG("here>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
}
