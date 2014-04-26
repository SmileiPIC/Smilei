#include "DiagnosticPhase.h"

using namespace std;

#include "SmileiMPI.h"

DiagnosticPhase::~DiagnosticPhase() {
}

DiagnosticPhase::DiagnosticPhase(phaseStructure phaseStruct, hid_t gid) :
groupID(gid),
my_species(phaseStruct.species)
{
	every=phaseStruct.every;
	if (every==0) ERROR("every cannot be zero");
}

void DiagnosticPhase::close(SmileiMPI* smpi) {
	if (smpi->isMaster()) {
		DEBUG("here");
		H5Gclose(groupID);
	}
}

