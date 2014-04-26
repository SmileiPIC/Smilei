#include "DiagnosticPhase2DxPx.h"

using namespace std;

DiagnosticPhase2DxPx::~DiagnosticPhase2DxPx() {
}

DiagnosticPhase2DxPx::DiagnosticPhase2DxPx(phaseStructure phaseStruct) : DiagnosticPhase2D(phaseStruct) {

	if (phaseStruct.pos_num.size() >0 && phaseStruct.mom_num.size() >0) {
		my_data.allocateDims(phaseStruct.pos_num[0],phaseStruct.mom_num[0]);
	} else {
		ERROR("must define pos_ and mom_ stuff");
	}

	xmin = phaseStruct.pos_min[0];
	xmax = phaseStruct.pos_max[0];
	xnum = phaseStruct.pos_num[0];
	pxmin = phaseStruct.mom_min[0];
	pxmax = phaseStruct.mom_max[0];
	pxnum = phaseStruct.mom_num[0];
		
}

void DiagnosticPhase2DxPx::doSomething(partStruct& my_part) {
	if (my_part.pos[0] > xmin && my_part.pos[0] < xmax && my_part.mom[0] > pxmin && my_part.mom[0] < pxmax) {
		int i = xnum*(my_part.pos[0]-xmin)/(xmax-xmin);
		int j = pxnum*(my_part.mom[0]-pxmin)/(pxmax-pxmin);
		my_data(i,j)+=my_part.weight;		
//	} else {
//		DEBUG(my_part.pos[0] << " " << my_part.mom[0]);
	}
}

