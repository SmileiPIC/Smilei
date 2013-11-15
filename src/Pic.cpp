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
#include "PicParams.h"
#include "InputData.h"
#include "DiagParams.h"

#include "SmileiMPIFactory.h"
#include "SmileiIOFactory.h"

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include <ctime>
#include <cstdlib>
#include <iostream>

#include "Diagnostic.h"

using namespace std;


// ------------------------------------------------------------------------------------------------------------------ //
//                                                   MAIN CODE
// ------------------------------------------------------------------------------------------------------------------ //
int main (int argc, char* argv[])
{
	std::cout.setf( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed

	// Define 2 MPI environment :
	//  - smpiData : to broadcast input data, unknown geometry
	//  - smpi (defined later) : to compute/exchange data, specific to a geometry
 	SmileiMPI* smpiData = new SmileiMPI( &argc, &argv );
   
	// -------------------------
	// Simulation Initialization
	// -------------------------
     
	// Check for namelist (input file)
	string namelist;
	if (argc<2) ERROR("No namelists given!");
	namelist=argv[1];

	// Send information on current simulation
	if ( smpiData->isMaster() ) startingMessage(namelist);

	// Parse the namelist file (no check!)
	InputData input_data(namelist);
	DEBUGEXEC(input_data.write(namelist+".debug","parsed namelist"));

	// Read simulation parameters
	PicParams params;
	DiagParams diag_params;

	// Process 0 read namelist, then broadcast to all process
	if ( smpiData->isMaster() ) {
		params.parseInputData(input_data); // this variable will hold the input parameters from file
		diag_params.parseInputData(input_data, params); // this variable will hold the diagnostics parameters from file
	}
	
	// Brodcast importa parameters to all nodes
	smpiData->bcast( params );
	if ( smpiData->isMaster() ) params.print();

	smpiData->bcast( diag_params );


	// Geometry known, MPI environment specified
	SmileiMPI* smpi = SmileiMPIFactory::create(params, smpiData);
	SmileiIO*  sio  = SmileiIOFactory::create(params, smpi);


	// Create diagnostic
	Diagnostic diags(&params,&diag_params, smpi);
	
	
	// Randomize the seed for simulations running in release mode
	//! \todo{Save the seed in case one wants to re-run the exact same simulation (MG)}
	RELEASEEXEC(srand (time(NULL)));

	// -------------------------------------------
	// Declaration of the main objects & operators
	// -------------------------------------------

	// ------------------------------------------------------------------------------------
	// Initialize the vecSpecies object containing all information of the different Species
	// ------------------------------------------------------------------------------------
	// vector of Species (virtual)
	vector<Species*> vecSpecies = SpeciesFactory::createVector(params, smpi);
	// dump species at time 0
	sio->writePlasma( vecSpecies, 0., smpi );

	// ----------------------------------------------------------------------------
	// Initialize the electromagnetic fields and interpolation-projection operators
	// according to the simulation geometry
	// ----------------------------------------------------------------------------
	// object containing the electromagnetic fields (virtual)
	ElectroMagn* EMfields = ElectroMagnFactory::create(params, smpi);

	// interpolation operator (virtual)
	Interpolator* Interp = InterpolatorFactory::create(params, smpi);
    
	// projection operator (virtual)
	Projector* Proj = ProjectorFactory::create(params, smpi);
    
	// -----------------------------------
	// Initialize the electromagnetic fields
	// -----------------------------------   
	//!\todo{Check & describe what is done here (MG)}
	// Init rho by pro all particles of subdomain -> local stuff
	EMfields->initRho(vecSpecies, Proj);
    
	//smpi->sumRho( EMfields );
    smpi->sumDensities( EMfields );
    
	//! \todo{FalseNot //, current algorithm is instrinsicaly sequential}
	smpi->solvePoissonPara( EMfields );		//champs->initMaxwell();

    smpi->barrier();
    
    
	// ------------------------------------------------------------------------
	// Initialize the simulation times time_prim at n=0 and time_dual at n=-1/2
	// ------------------------------------------------------------------------
	// time at integer time-steps (primal grid)
	double time_prim = 0.;
	// time at half-integer time-steps (dual grid)
	double time_dual = -0.5 * params.timestep;

	// ------------------------------------------------------------------
	//                     HERE STARTS THE PIC LOOP
	// ------------------------------------------------------------------
	if ( smpi->isMaster() ) MESSAGE(0,"Time-Loop is started: number of time-steps n_time =" << params.n_time);
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
		if ( (itime % 100 == 0) &&  ( smpi->isMaster() ) )
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
// 			if ( smpi->isMaster() ) DEBUG(2, "Dynamic Species " << ispec );
			vecSpecies[ispec]->dynamic(time_dual, EMfields, Interp, Proj, smpi);
			smpi->exchangeParticles(vecSpecies[ispec], ispec, &params);
		}
		smpi->sumDensities( EMfields );

		// solve Maxwell's equations
		EMfields->solveMaxwell(time_dual, smpi);

        // call the various diagnostics
		// ----------------------------
		
		diags.runAllDiags(itime, EMfields, vecSpecies);
	}//END of the time loop	

	smpi->barrier();
	t1 = MPI_Wtime();
	if ( smpi->isMaster() ) MESSAGE(0, "End time loop, time dual = " << time_dual);
	if ( smpi->isMaster() ) MESSAGE(0, "Time in time loop : " << t1-t0 );
	// ------------------------------------------------------------------
	//                      HERE ENDS THE PIC LOOP
	// ------------------------------------------------------------------
	
		
	// ------------------------------------------------------------------
	//                      Temporary validation diagnostics
	// ------------------------------------------------------------------

	// 1 HDF5 file per process
	sio->writePlasma( vecSpecies, time_dual, smpi );
		
	//EMfields->initRho(vecSpecies, Proj);
	//smpi->sumRho( EMfields );

	//EMfields->dump(&params);  	// Sequential results, 1 file per process
	if (params.nDim_field == 1) { // If 1D
			//! \todo{Not //, processes write sequentially to validate. OK in 1D}
		smpi->writeFields( EMfields );
	}
	else { // If 2D
		sio->writeFields( EMfields );
		sio->writeFieldsPP( EMfields, time_dual, smpi->getRank() );
	}

 
	// ------------------------------
	//  Cleanup & End the simulation
	// ------------------------------
	delete Proj;
	delete Interp;
	delete EMfields;
	for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
	vecSpecies.clear();
  
	delete sio;
	if ( smpi->isMaster() ) {
		MESSAGE("------------------------------------------");
		MESSAGE("END " << namelist);
		MESSAGE("------------------------------------------");
	}
	delete smpi;
	delete smpiData;

	return 0;
    
}//END MAIN 
