#include "DiagnosticPhase2DxPx.h"

using namespace std;

DiagnosticPhase2DxPx::~DiagnosticPhase2DxPx() {
}

DiagnosticPhase2DxPx::DiagnosticPhase2DxPx(phaseStructure phaseStruct, hid_t gid) : DiagnosticPhase2D(phaseStruct, gid) {
	DEBUG("here>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
	
	if (phaseStruct.pos_num.size() >0 && phaseStruct.mom_num.size() >0) {
		my_data.allocateDims(phaseStruct.pos_num[0],phaseStruct.mom_num[0]);
	} else {
		ERROR("must define pos_ and mom_ stuff");
	}

}

void DiagnosticPhase2DxPx::doSomething(partStruct) {
	
	DEBUG("here");
}
