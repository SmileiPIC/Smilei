#include "Diagnostic.h"
#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "Interpolator.h"
#include "DiagnosticScalar.h"

#include <string>

using namespace std;

<<<<<<< HEAD
Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi, Interpolator* interp ) : 
diagScal(params, smpi),
probe0d(),
interp_(interp)
=======
Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi) : 
diagScal(params, smpi),
probe0d(params, smpi, diagparams->ps_coor)
>>>>>>> d9430c98fb22fc65390629334be5bfeebc514228
{
	everyScalar = diagparams->scalar_every;
	everyMap = diagparams->map_every;
    everyProbe0D = diagparams->probe0d_every;
<<<<<<< HEAD
=======

>>>>>>> d9430c98fb22fc65390629334be5bfeebc514228
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator* interp){	
	if (everyScalar && timestep % everyScalar == 0) {
		diagScal.run(timestep, EMfields, vecSpecies);
<<<<<<< HEAD
	}	
	if (everyProbe0D && timestep % everyProbe0D == 0) {
		probe0d.run(timestep, EMfields, interp_);
=======
>>>>>>> d9430c98fb22fc65390629334be5bfeebc514228
	}	
//	if (everyProbe0D && timestep % everyProbe0D == 0) {
//		probe0d.run(timestep, EMfields, interp);
//	}	
}

