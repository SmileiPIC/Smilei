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

#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"

#include <ctime>
#include <cstdlib>

#include <iostream>

using namespace std;

int main (int argc, char* argv[]) {

	//! \todo{convert SmileiMPI in SmileiMPI_Cart1D}
 	// SmileiMPI smpi( &argc, &argv )
 	SmileiMPI_Cart1D smpi( &argc, &argv );

	if ( smpi.isMaster() ) {
		MESSAGE("------------------------------------------");
		MESSAGE(" Version : " << __VERSION DEBUGEXEC(<< " DEBUG") << " Compiled : " << __DATE__ << " " << __TIME__);
		MESSAGE("------------------------------------------");
	}

	if (argc<2) ERROR("No namelists given!");

	for (int arg=1; arg<argc;arg++) {
		
		string namelist=argv[arg];
		/*******************************************************************************************************************
		 Simulation init
		 ******************************************************************************************************************/
		if ( smpi.isMaster() ) {
			MESSAGE(" Namelist  : " << namelist);
			MESSAGE("------------------------------------------");
		}

		/*******************************************************************************************************************
		 Simulation parameters
		 ******************************************************************************************************************/
		PicParams params;

		// Process 0 read namelist, then broadcast
		if ( smpi.isMaster() ) 
			params.parseFile(namelist); // this variable will hold the imput parameters from file
		if ( smpi.isMaster() ) params.print();
		smpi.bcast( params );

		// Creation of a cartesian topology
		smpi.createTopology();
		
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
		
		//! this will randomizie each simulation done just in the release mode
	        //! \todo{check if one wants to run the same again: save seed}
		RELEASEEXEC(srand (time(NULL)));
		
		vecSpecies.resize(params.n_species);
		for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
			PMESSAGE( 0, smpi.getRank(), "Initializing Species "<<ispec);
			Species* sp = NULL;
			if (params.species_param[ispec].dynamics_type=="norm") {
				sp = new Species_norm(&params, ispec, &smpi);
			} else if (params.species_param[ispec].dynamics_type=="rrll") {
				sp = new Species_rrll(&params, ispec, &smpi);
			}
			vecSpecies[ispec] = sp;
			sp->dump(ofile);
			ofile << endl;

			smpi.exchangeParticles(vecSpecies[ispec], &params);
			//MESSAGE( 0, "\tProcess " << smpi.getRank() << " : " << sp->getNbrOfParticles() << " Particles of species " << ispec );
			PMESSAGE( 0, smpi.getRank(), sp->getNbrOfParticles() << " Particles of species " << ispec );
		}// END for ispec


		// allocate
		if ( params.geometry == "1d3v" ) {
		  champs = new ElectroMagn1D(&params, &smpi);
			if ( params.interpolation_order == 2 )
			{
				Interp = new Interpolator1D2Order(&params, &smpi);
				Proj   = new Projector1D2Order(&params, &smpi);
			}
		}
		else {
			ERROR( "Unknwon geometry : " << params.geometry );
		}

		// Init rho by pro all particles of subdomain -> local stuff
		champs->initRho(vecSpecies, Proj);
		smpi.sumRho( champs );
		//! \todo{FalseNot //, current algorithm is instrinsically sequential}
		smpi.initMaxwellPara( champs );		//champs->initMaxwell();
		
		// ------------------------------------------------------------------
		// ------------------------------------------------------------------
		// ------------------------------------------------------------------
		
		//! \todo{clarify this}
		double time_dual = 0.; //-params.timestep/2.;
		if ( smpi.isMaster() ) MESSAGE( 0, "Start time loop" );

		// t1-t0  = elapsed time in simulation time loop
		double t0, t1;
		t0 = MPI_Wtime();

		for (unsigned int itps=1 ; itps <= params.n_time ; itps++) {
			//calculate new time
			time_dual += params.timestep ;
			if (itps % 100 == 0) {
				if ( smpi.isMaster() ) MESSAGE( 1, "Time: " << time_dual << " " << itps );
				//if ( smpi.isMaster() ) champs->dump();
			}
			
			//put density and currents to 0
			champs->initRhoJ();
			
			//plasma
			for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
				if ( smpi.isMaster() ) DEBUG( 2, "Dynamic Species " << ispec );
				vecSpecies[ispec]->dynamic(time_dual, champs, Interp, Proj, &smpi);
				smpi.exchangeParticles(vecSpecies[ispec], &params);
				DEBUG( 2, "\tProcess " << smpi.getRank() << " : " << vecSpecies[ispec]->getNbrOfParticles() << " Particles of species " << ispec << " in loop" );
			}
			smpi.sumDensities( champs );

			//! \todo{Not //, current algorithm is instrinsically sequential}
			smpi.chargeConservingPara( champs );	//champs->chargeConserving();
			champs->solveMaxwell(time_dual, params.timestep, &smpi);

		}

		smpi.barrier();
		t1 = MPI_Wtime();
		if ( smpi.isMaster() ) MESSAGE(0, "Time in time loop : " << t1-t0 );
		if ( smpi.isMaster() ) MESSAGE(0, "End time loop");
		
		if ( smpi.isMaster() ) {
			for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
				vecSpecies[ispec]->dump(ofile);
				ofile << endl;
			}
		}

		//! \todo{Not //, processes write sequentially to validate. OK in 1D}
		smpi.writePlasma( vecSpecies, "dump_new" );  
		
		if ( smpi.isMaster() ) champs->dump();
		
		//! \todo{Not //, processes write sequentially to validate. OK in 1D}
		smpi.writeField( champs->Ex_, "fex_new" );
		smpi.writeField( champs->Ey_, "fey_new" );
		smpi.writeField( champs->Ez_, "fez_new" );
		smpi.writeField( champs->Bx_, "fbx_new" );
		smpi.writeField( champs->By_, "fby_new" );
		smpi.writeField( champs->Bz_, "fbz_new" );
		smpi.writeField( champs->Jx_, "fjx_new" );
		smpi.writeField( champs->Jy_, "fjy_new" );
		smpi.writeField( champs->Jz_, "fjz_new" );		
		smpi.writeField( champs->rho_, "rho_new" );

		/*******************************************************************************************************************
		 cleanup
		 ******************************************************************************************************************/
		delete Proj;
		delete Interp;
		delete champs;
		for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
		vecSpecies.clear();

		if ( smpi.isMaster() ) {
			MESSAGE("------------------------------------------");
			MESSAGE("END " << namelist);
			MESSAGE("------------------------------------------");
		}
	}

	return 0;
}


