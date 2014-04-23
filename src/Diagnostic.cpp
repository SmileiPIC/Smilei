#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "Interpolator.h"
#include "DiagnosticScalar.h"

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi, Interpolator *interp) :
diagScal(params, smpi),
probe0D(params, diagparams, smpi),
interp_(interp),
diagPhase(params, diagparams, smpi)
{
    everyScalar = diagparams->scalar_every;
    everyMap = diagparams->map_every;
    everyProbe0D = diagparams->probe0d_every;
}

void Diagnostic::closeAll () {
    probe0D.close();
	diagPhase.close();
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
    if (everyScalar && timestep % everyScalar == 0) {
        diagScal.run(timestep, EMfields, vecSpecies);
    }
    if (everyProbe0D && timestep % everyProbe0D == 0) {
        probe0D.run(timestep, EMfields, interp_);
    }
	diagPhase.run(timestep, vecSpecies);
}

