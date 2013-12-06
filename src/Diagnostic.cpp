#include "Diagnostic.h"
#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "DiagnosticScalar.h"

#include <string>

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi) : 
diagScal(params, smpi),
probe0d(params, smpi, diagparams->ps_coor)
{
	everyScalar = diagparams->scalar_every;
	everyMap = diagparams->map_every;
    everyProbe0D = diagparams->probe0d_every;

}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator* interp){	
	if (everyScalar && timestep % everyScalar == 0) {
		diagScal.run(timestep, EMfields, vecSpecies);
	}	
//	if (everyProbe0D && timestep % everyProbe0D == 0) {
//		probe0d.run(timestep, EMfields, interp);
//	}	
}

