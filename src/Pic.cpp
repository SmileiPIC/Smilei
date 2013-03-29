#include "Pic.h"
#include "Tools.h"
#include "PicParams.h"

#include "Species.h"
#include "Species_norm.h"
#include "Species_rrll.h"
#include "ElectroMagn.h"
#include "ElectroMagn1D.h"

#include "Pusher.h"
#include "PusherBoris.h"
#include "PusherSklv.h"
#include "Interpolator.h"
#include "Interpolator1D2Order.h"
#include "Interpolator1D3Order.h"
#include "Projector.h"
#include "Projector1D2Order.h"

#include "Field3D.h"
#include "Field1D.h"

#include <ctime>
#include <cstdlib>

#include <iostream>

using namespace std;

int main (int argc, char* argv[]) {
	MESSAGE("------------------------------------------");
	MESSAGE(" Version : " << __VERSION DEBUGEXEC(<< " DEBUG") << " Compiled : " << __DATE__ << " " << __TIME__);
	MESSAGE("------------------------------------------");
	if (argc<2) ERROR("No namelists given!");
	for (int arg=1; arg<argc;arg++) {
		
		string namelist=argv[arg];
		/*******************************************************************************************************************
		 Simulation init
		 ******************************************************************************************************************/
		MESSAGE(" Namelist  : " << namelist);
		MESSAGE("------------------------------------------");

		/*******************************************************************************************************************
		 Simulation parameters
		 ******************************************************************************************************************/
		PicParams params(namelist);// this variable will hold the imput parameters from file
		
		//! this will randomizie each simulation done just in the release mode
        //! \todo{check if one wants to run the same again: save seed}
		RELEASEEXEC(srand (time(NULL)));
		
		/*******************************************************************************************************************
		 Variable declaration
		 ******************************************************************************************************************/
		vector<Species*> vecSpecies; // vector of species
		ElectroMagn* champs = NULL; // fields
		
		// operators
		Interpolator* Interp = NULL;
		Projector* Proj = NULL;
		
		// initialise Species
		ofstream ofile("dump", ios::out);
		
		vecSpecies.resize(params.n_species);
		for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
			MESSAGE(0,"Initializing Species "<<ispec);
			Species* sp = NULL;
			if (params.species_param[ispec].dynamics_type=="norm") {
				sp = new Species_norm(&params, ispec);
			} else if (params.species_param[ispec].dynamics_type=="rrll") {
				sp = new Species_rrll(&params, ispec);
			}
			sp->dump(ofile);
			ofile << endl;
			vecSpecies[ispec] = sp;
			MESSAGE(0,sp->getNbrOfParticles() << " Particles of species " << ispec);
		}// END for ispec
		
		// allocate
		if ( params.geometry == "1d3v" ) {
			champs = new ElectroMagn1D(&params);
			if ( params.interpolation_order == 2 )
			{
				Interp = new Interpolator1D2Order(&params);
				Proj   = new Projector1D2Order(&params);
			}
		}
		else {
			ERROR( "Unknwon geometry : " << params.geometry );
		}
		
		champs->initRho(vecSpecies, Proj);
		champs->initMaxwell();
		
		// ------------------------------------------------------------------
		// ------------------------------------------------------------------
		// ------------------------------------------------------------------
		
		//! \todo{clarify this}
		double time_dual = 0.; //-params.timestep/2.;
		MESSAGE(0,"Start time loop");
		
		for (unsigned int itps=1 ; itps <= params.n_time ; itps++) {
			
			//calculate new time
			time_dual += params.timestep ;
			if (itps % 100 == 0) {
				MESSAGE(1,"Time: " << time_dual << " " << itps);
				champs->dump();
			}
			
			//put density and currents to 0
			champs->initRhoJ();
			
			//plasma
			for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
				DEBUG(2, "Dynamic Species "<<ispec );
				vecSpecies[ispec]->dynamic(time_dual, champs, Interp, Proj);
			}
			champs->chargeConserving();
			champs->solveMaxwell(time_dual, params.timestep);
		}
		
		MESSAGE(0,"End time loop");
		
		for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
			vecSpecies[ispec]->dump(ofile);
			ofile << endl;
		}
		
		
		champs->dump();
		
		
		
		/*******************************************************************************************************************
		 cleanup
		 ******************************************************************************************************************/
		
		delete Proj;
		delete Interp;
		delete champs;
		for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
		vecSpecies.clear();

		MESSAGE("------------------------------------------");
		MESSAGE("END " << namelist);
		MESSAGE("------------------------------------------");
	}
	return 0;
}


