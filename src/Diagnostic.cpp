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
everyScalar (diagparams->scalar_every),
probe0D(params, diagparams, smpi),
probe1D(params, diagparams, smpi),
probe2D(params, diagparams, smpi),
interp_(interp),
diagPhaseSpace(params, diagparams, smpi)
{
}

Diagnostic::~Diagnostic () {
    probe0D.close();
    probe1D.close();
    probe2D.close();
	diagPhaseSpace.close();
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
    if (everyScalar && timestep % everyScalar == 0) {
        diagScal.run(timestep, EMfields, vecSpecies);
    }

    probe0D.runAll(timestep, EMfields, interp_);
    probe1D.runAll(timestep, EMfields, interp_);
    probe2D.runAll(timestep, EMfields, interp_);
        
	diagPhaseSpace.run(timestep, vecSpecies);
}

