#include "Diagnostic.h"
#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "DiagnosticScalar.h"

#include <string>

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi ) : DiagScal(params, smpi) {
	everyScalar = diagparams->scalar_every;
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies){	
	if (everyScalar && timestep % everyScalar == 0) {
		DiagScal.run(timestep, EMfields, vecSpecies);
	}	
}
