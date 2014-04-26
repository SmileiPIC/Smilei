#include "DiagnosticPhase.h"

using namespace std;

#include "SmileiMPI.h"

DiagnosticPhase::~DiagnosticPhase() {
}

DiagnosticPhase::DiagnosticPhase(phaseStructure phaseStruct) :
my_species(phaseStruct.species)
{
	every=phaseStruct.every;
	if (every==0) ERROR("every cannot be zero");
}


