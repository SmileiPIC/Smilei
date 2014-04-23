#include "DiagnosticPhase.h"

using namespace std;

DiagnosticPhase::~DiagnosticPhase() {
}

DiagnosticPhase::DiagnosticPhase(phaseStructure phaseStruct, hid_t gid) :
groupID(gid),
my_species(phaseStruct.species)
{
	every=phaseStruct.every;
	if (every==0) ERROR("every cannot be zero");
}

void DiagnosticPhase::close() {
	H5Gclose(groupID);
}

