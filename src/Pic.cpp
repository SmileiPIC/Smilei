////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                                                ////
////                                                                                                                ////
////                                   PARTICLE-IN-CELL CODE SMILEI                                                 ////
////                    Simulation of Matter Irradiated by Laser at Extreme Intensity                               ////
////                                                                                                                ////
////                          Cooperative OpenSource Object-Oriented Project                                        ////
////                                      from the Plateau de Saclay                                                ////
////                                          started January 2013                                                  ////
////                                                                                                                ////
////                                                                                                                ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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


// ------------------------------------------------------------------------------------------------------------------ //
//                                                   MAIN CODE
// ------------------------------------------------------------------------------------------------------------------ //
int main (int argc, char* argv[])
{
  std::cout.setf( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
	//! \todo{convert SmileiMPI in SmileiMPI_Cart1D}
 	// SmileiMPI smpi( &argc, &argv )
 	SmileiMPI_Cart1D smpi( &argc, &argv );
   
	// -------------------------
	// Simulation Initialisation
	// -------------------------
     
	// Check for namelist (input file)
	string namelist;
	if (argc<2) ERROR("No namelists given!");
	namelist=argv[1];

	// Send information on current simulation
	if ( smpi.isMaster() ) {
		MESSAGE("------------------------------------------");
		MESSAGE(" Version : " << __VERSION DEBUGEXEC(<< " DEBUG") << " Compiled : " << __DATE__ << " " << __TIME__);
		MESSAGE("------------------------------------------");
		MESSAGE(" Namelist  : " << namelist);
		MESSAGE("------------------------------------------");
	}

	// Read simulation parameters
	PicParams params;

	// Process 0 read namelist, then broadcast
	if ( smpi.isMaster() ) 
		params.parseFile(namelist); // this variable will hold the imput parameters from file
	if ( smpi.isMaster() ) params.print();
	smpi.bcast( params );

	// Creation of a cartesian topology
	smpi.createTopology();

	// Randomize the seed for simulations running in release mode
	//! \todo{Save the seed in case one wants to re-run the exact same simulation (MG)}
	RELEASEEXEC(srand (time(NULL)));

	// -------------------------------------------
	// Declaration of the main objects & operators
	// -------------------------------------------
    
	// object containing the electromagnetic fields (virtual)
	ElectroMagn* EMfields = NULL;
    
	// interpolation operator (virtual)
	Interpolator* Interp = NULL;
    
	// projection operator (virtual)
	Projector* Proj = NULL;
    
	// vector of Species (virtual)
	vector<Species*> vecSpecies;
    
	// species "dump" file
	//! \todo{Check if we keep going like that (MG)}
	ofstream ofile("dump", ios::out);

	// ------------------------------------------------------------------------------------
	// Initialize the vecSpecies object containing all information of the different Species
	// ------------------------------------------------------------------------------------

	vecSpecies.resize(params.n_species);
	for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
		PMESSAGE( 0, smpi.getRank(), "Initializing Species "<<ispec);
		Species* sp = NULL;
		if (params.species_param[ispec].dynamics_type=="norm") {
			// Species with Boris dynamics
			sp = new Species_norm(&params, ispec, &smpi);
		} else if (params.species_param[ispec].dynamics_type=="rrll") {
			// Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
			sp = new Species_rrll(&params, ispec, &smpi);
		}//endif

		//save temporary species sp in vecSpecies
		vecSpecies[ispec] = sp;

		//dump species at time 0
		sp->dump(ofile); ofile << endl;

		//PMESSAGE( 0, smpi.getRank(), sp->getNbrOfParticles() << " Particles of species " << ispec );
		smpi.exchangeParticles(vecSpecies[ispec], &params);
		PMESSAGE( 0, smpi.getRank(), sp->getNbrOfParticles() << " Particles of species " << ispec );
	}// END for ispec

	// ----------------------------------------------------------------------------
	// Initialize the electromagnetic fields and interpolation-projection operators
	// according to the simulation geometry
	// ----------------------------------------------------------------------------
	if ( params.geometry == "1d3v" ) {
		// ---------------
		// 1d3v Simulation
		// ---------------
		EMfields = new ElectroMagn1D(&params, &smpi);
		if ( params.interpolation_order == 2 )
			{
				Interp = new Interpolator1D2Order(&params, &smpi);
				Proj   = new Projector1D2Order(&params, &smpi);
			}
	}
	else {
		ERROR( "Unknwon geometry : " << params.geometry );
	}//endif params.geometry


	// -----------------------------------
	// Inialize the electromagnetic fields
	// -----------------------------------   
	//!\todo{Check & describe what is done here (MG)}
	// Init rho by pro all particles of subdomain -> local stuff
	EMfields->initRho(vecSpecies, Proj);
	smpi.sumRho( EMfields );

	//! \todo{FalseNot //, current algorithm is instrinsically sequential}
	smpi.solvePoissonPara( EMfields );		//champs->initMaxwell();


	// ------------------------------------------------------------------------
	// Initialise the simulation times time_prim at n=0 and time_dual at n=-1/2
	// ------------------------------------------------------------------------
	// time at integer time-steps (primal grid)
	double time_prim = 0.;
	// time at half-integer time-steps (dual grid)
	double time_dual = -0.5 * params.timestep;

	// ------------------------------------------------------------------
	//                     HERE STARTS THE PIC LOOP
	// ------------------------------------------------------------------
	if ( smpi.isMaster() ) MESSAGE(0,"Time-Loop is started: number of time-steps n_time =" << params.n_time);
	// t1-t0  = elapsed time in simulation time loop
	double t0, t1;
	t0 = MPI_Wtime();
	for (unsigned int itime=1 ; itime <= params.n_time ; itime++) {
        
		// calculate new times
		// -------------------
		time_prim += params.timestep;
		time_dual += params.timestep; 

		// send message at given time-steps
		// --------------------------------
		//!\todo{Introduce a control parameter in PicParams (MG)}
		if ( (itime % 100 == 0) &&  ( smpi.isMaster() ) )
			MESSAGE(1,"Time: " << time_dual << " " << itime);

		// put density and currents to 0
		// -----------------------------
		EMfields->initRhoJ();

		// apply the PIC method
		// --------------------
		// for all particles of all species (see dunamic in Species.cpp)
		// (1) interpolate the fields at the particle position
		// (2) move the particle
		// (3) calculate the currents (charge conserving method)
		for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
			if ( smpi.isMaster() ) DEBUG(2, "Dynamic Species "<<ispec );
			vecSpecies[ispec]->dynamic(time_dual, EMfields, Interp, Proj, &smpi);
			smpi.exchangeParticles(vecSpecies[ispec], &params);
			//DEBUG( 2, "\tProcess " << smpi.getRank() << " : " << vecSpecies[ispec]->getNbrOfParticles() << " Particles of species " << ispec << " in loop" );
		}
		//EMfields->dump(&params);
		smpi.sumDensities( EMfields );

		// calculate the longitudinal current using the charge conservation equation
		//! \todo{Not //, current algorithm is instrinsically sequential}
		//smpi.chargeConservingPara( EMfields);	//EMfields->chargeConserving();
		
		// solve Maxwell's equations
		EMfields->solveMaxwell(time_dual, params.timestep, &smpi);


		/*if (itime == 10000) {
		  smpi.writeField( EMfields->Ex_, "fex_new" );
		  smpi.writeField( EMfields->Ey_, "fey_new" );
		  smpi.writeField( EMfields->Ez_, "fez_new" );
		  smpi.writeField( EMfields->Bx_, "fbx_new" );
		  smpi.writeField( EMfields->By_, "fby_new" );
		  smpi.writeField( EMfields->Bz_, "fbz_new" );
		  smpi.writeField( EMfields->Jx_, "fjx_new" );
		  smpi.writeField( EMfields->Jy_, "fjy_new" );
		  smpi.writeField( EMfields->Jz_, "fjz_new" );		
		  smpi.writeField( EMfields->rho_, "rho_new" );

		  return 0;
		}*/

	        // call the various diagnostics
		// ----------------------------
		if (itime % 5000 == 0) {
			if ( smpi.isMaster() ) MESSAGE(1,"diags at " << time_dual << " " << itime);
			//EMfields->dump(&params);
			for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
			  //vecSpecies[ispec]->dump(ofile);
				ofile << endl;
			}
		}
		
	}//END of the time loop	

	smpi.barrier();
	t1 = MPI_Wtime();
	if ( smpi.isMaster() ) MESSAGE(0, "Time in time loop : " << t1-t0 );
	if ( smpi.isMaster() ) MESSAGE(0, "End time loop, time dual = " << time_dual);
	// ------------------------------------------------------------------
	//                      HERE ENDS THE PIC LOOP
	// ------------------------------------------------------------------
	
		
	// ------------------------------------------------------------------
	//                      Temporary validation diagnostics
	// ------------------------------------------------------------------
	if ( smpi.isMaster() ) {
		for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
			vecSpecies[ispec]->dump(ofile);
			ofile << endl;
		}
	}
	//! \todo{Not //, processes write sequentially to validate. OK in 1D}
	smpi.writePlasma( vecSpecies, "dump_new" );  
		
	//if ( smpi.isMaster() ) 
	  EMfields->dump(&params);	
	//! \todo{Not //, processes write sequentially to validate. OK in 1D}
	smpi.writeField( EMfields->Ex_, "fex_new" );
	smpi.writeFieldPrim( EMfields->Ey_, "fey_new" );
	smpi.writeFieldPrim( EMfields->Ez_, "fez_new" );
	smpi.writeFieldPrim( EMfields->Bx_, "fbx_new" );
	smpi.writeField( EMfields->By_, "fby_new" );
	smpi.writeField( EMfields->Bz_, "fbz_new" );
	smpi.writeField( EMfields->Jx_, "fjx_new" );
	smpi.writeFieldPrim( EMfields->Jy_, "fjy_new" );
	smpi.writeFieldPrim( EMfields->Jz_, "fjz_new" );		
	smpi.writeFieldPrim( EMfields->rho_, "rho_new" );

	// ------------------------------
	//  Cleanup & End the simulation
	// ------------------------------
	delete Proj;
	delete Interp;
	delete EMfields;	for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
	vecSpecies.clear();
    
	if ( smpi.isMaster() ) {
		MESSAGE("------------------------------------------");
		MESSAGE("END " << namelist);
		MESSAGE("------------------------------------------");
	}

	return 0;
    
}//END MAIN 
