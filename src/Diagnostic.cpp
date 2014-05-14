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
interp_(interp),
diagPhaseSpace(params, diagparams, smpi)
{
}

Diagnostic::~Diagnostic () {
    probe0D.close();
    probe1D.close();
	diagPhaseSpace.close();
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies) {
    if (everyScalar && timestep % everyScalar == 0) {
        diagScal.run(timestep, EMfields, vecSpecies);
    }

    for (unsigned int i=0; i<probe0D.every.size(); i++) {
        if (probe0D.every[i] && timestep % probe0D.every[i] == 0) {
            probe0D.run(i, EMfields, interp_);
        }
    }
    
    for (unsigned int i=0; i<probe1D.every.size(); i++) {
        if (probe1D.every[i] && timestep % probe1D.every[i] == 0) {
            probe1D.run(i, EMfields, interp_);
        }
    }
	diagPhaseSpace.run(timestep, vecSpecies);
}

