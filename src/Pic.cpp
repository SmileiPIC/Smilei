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

#include <ctime>
#include <cstdlib>
#include <unistd.h>

#include <iostream>
#include <iomanip>

#include "PicParams.h"
#include "InputData.h"
#include "DiagParams.h"

#include "SmileiMPIFactory.h"
#include "SmileiIOFactory.h"

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "Diagnostic.h"

#include "SimWindow.h"

#include "Timer.h"
#include <omp.h>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    std::cout.setf( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
    
    // Define 2 MPI environments :
    //  - smpiData : to broadcast input data, unknown geometry
    //  - smpi (defined later) : to compute/exchange data, specific to a geometry
    SmileiMPI *smpiData= new SmileiMPI(&argc, &argv );
    
    // -------------------------
    // Simulation Initialization
    // ------------------------- 
    
    argc -= optind;
    argv += optind;
    
    // Check for namelist (input file)
    if (argc<1) ERROR("No namelists given!");
    string namelist=argv[0];
    
    // Send information on current simulation
    if ( smpiData->isMaster() ) startingMessage(namelist);
    
    // Parse the namelist file (no check!)
    InputData input_data;
    if ( smpiData->isMaster() ) input_data.readFile(namelist);    

    // broadcast file and parse it and randomize
    smpiData->bcast(input_data);    
    
    // Read simulation & diagnostics parameters
    PicParams params(input_data);
    smpiData->init(params);
    smpiData->barrier();
    
    // Print out the data parameters
    MESSAGE("----------------------------------------------");
    MESSAGE("Input data info");
    MESSAGE("----------------------------------------------");
    if ( smpiData->isMaster() ) params.print();
    
    
    // Geometry known, MPI environment specified
    MESSAGE("----------------------------------------------");
    MESSAGE("Creating MPI & IO environments");
    MESSAGE("----------------------------------------------");
    
    SmileiMPI* smpi = SmileiMPIFactory::create(params, smpiData);
    SmileiIO*  sio  = SmileiIOFactory::create(params, smpi);
    
    
    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------
    MESSAGE("----------------------------------------------");
    MESSAGE("Creating EMfields/Interp/Proj/Diags");
    MESSAGE("----------------------------------------------");
    
    // Initialize the electromagnetic fields and interpolation-projection operators
    // according to the simulation geometry
    // ----------------------------------------------------------------------------
    // object containing the electromagnetic fields (virtual)
    ElectroMagn* EMfields = ElectroMagnFactory::create(params, smpi);
    
    // interpolation operator (virtual)
    Interpolator* Interp = InterpolatorFactory::create(params, smpi);
    
    // projection operator (virtual)
    Projector* Proj = ProjectorFactory::create(params, smpi);
    
    // Create diagnostics
    Diagnostic diags(params,input_data, smpi);
    
    
    // ---------------------------
    // Initialize Species & Fields
    // ---------------------------
    MESSAGE("----------------------------------------------");
    MESSAGE("Initializing particles, fields & moving-window");
    MESSAGE("----------------------------------------------");
    
    // Initialize the vecSpecies object containing all information of the different Species
    // ------------------------------------------------------------------------------------
    
    // vector of Species (virtual)
    vector<Species*> vecSpecies = SpeciesFactory::createVector(params, smpi);

    // ----------------------------------------------------------------------------
    // Define Moving Window & restart
    // ----------------------------------------------------------------------------
    
    SimWindow* simWindow = NULL;
    if (params.res_space_win_x)
        simWindow = new SimWindow(params);
    smpi->barrier();
    
    unsigned int stepStart=0, stepStop=params.n_time;
    
    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        DEBUG(vecSpecies.size());
        sio->restartAll( EMfields,  stepStart, vecSpecies, smpi, simWindow, params, input_data);

	if (simWindow) {
	    simWindow->setOperators(vecSpecies, Interp, Proj, smpi);
	    simWindow->operate(vecSpecies, EMfields, Interp, Proj, smpi );
	}
	    
    } else {
        // Initialize the electromagnetic fields
        // -----------------------------------
        // Init rho and J by projecting all particles of subdomain
        EMfields->initRhoJ(vecSpecies, Proj);
        
        // Sum rho and J on ghost domains
        smpi->sumRhoJ( EMfields );
        
        // Init electric field (Ex/1D, + Ey/2D)
        if (smpiData->isMaster()) {
            MESSAGE("----------------------------------------------");
            MESSAGE("Solving Poisson at time t = 0");
            MESSAGE("----------------------------------------------");
        }
        EMfields->solvePoisson(smpi);
        
        
        MESSAGE("----------------------------------------------");
        MESSAGE("Running diags at time t = 0");
        MESSAGE("----------------------------------------------");
        // run diagnostics at time-step 0
        diags.runAllDiags(0, EMfields, vecSpecies, Interp, smpi);
        // temporary EM fields dump in Fields.h5
        sio->writeAllFieldsSingleFileTime( EMfields, 0 );
        // temporary EM fields dump in Fields_avg.h5
        sio->writeAvgFieldsSingleFileTime( EMfields, 0 );
        // temporary particle dump at time 0
        sio->writePlasma( vecSpecies, 0., smpi );
    }
    
    

    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // ------------------------------------------------------------------------
	
    // time at integer time-steps (primal grid)
    double time_prim = stepStart * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (stepStart +0.5) * params.timestep;
	
    // Count timer
    int ntimer(5);
    Timer timer[ntimer];
    timer[0].init(smpi, "global");
    timer[1].init(smpi, "particles");
    timer[2].init(smpi, "maxwell");
    timer[3].init(smpi, "diagnostics");
    timer[4].init(smpi, "densities");
    
    
	// ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    MESSAGE("Time-Loop is started: number of time-steps n_time = " << params.n_time);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
	
    for (unsigned int itime=stepStart+1 ; itime <= stepStop ; itime++) {
        
        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        // send message at given time-steps
        // --------------------------------
        timer[0].update();
        
        //double timElapsed=smpiData->time_seconds();
		if ( (itime % diags.params.print_every == 0) &&  ( smpi->isMaster() ) )
            MESSAGE(1, "t= "         << setw(11)                     << time_dual/(2*M_PI)
                    << " it= "       << setw(log10(params.n_time)+1) << itime  << "/" << params.n_time
                    << " sec: "      << setw(9)                      << timer[0].getTime()
                    << " E= "        << setw(9)                      << diags.getScalar("Etot")
                    //<< " Epart=  "   << setw(9)                      << diags.getScalar("Eparticles")
                    //<< " Efield= "   << setw(9)                      << diags.getScalar("EFields")
                    << " E_bal(%)= " << 100.0*diags.getScalar("Ebal_norm") );

        
        
        // put density and currents to 0 + save former density
        // ---------------------------------------------------
        EMfields->restartRhoJ();
        
        
        // apply the PIC method
        // --------------------
        // for all particles of all species (see dunamic in Species.cpp)
        // (1) interpolate the fields at the particle position
        // (2) move the particle
        // (3) calculate the currents (charge conserving method)
        timer[1].restart();
#pragma omp parallel shared (EMfields,time_dual,vecSpecies,smpi)
        {
            int tid(0);
#ifdef _OMP
            tid = omp_get_thread_num();
#endif
            for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
		vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, smpi, &params);
            }
            for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
#pragma omp barrier
                //#pragma omp master
                {
                    // Loop on dims to manage exchange in corners
                    for ( int iDim = 0 ; iDim<params.nDim_particle ; iDim++ )
                        smpi->exchangeParticles(vecSpecies[ispec], ispec, &params, tid);
                }
#pragma omp barrier
                    vecSpecies[ispec]->sort_part(params.cell_length[0]);
            }
        }
        timer[1].update();
        
		//!\todo To simplify : sum global and per species densities
        timer[4].restart();
        smpi->sumRhoJ( EMfields );
        EMfields->computeTotalRhoJ();
        timer[4].update();
        
        // solve Maxwell's equations
        timer[2].restart();
        EMfields->solveMaxwell(itime, time_dual, smpi, params, simWindow);
        timer[2].update();
        
        // incrementing averaged electromagnetic fields
        EMfields->incrementAvgFields(itime, diags.params.ntime_step_avg);
        
        // call the various diagnostics
        // ----------------------------
		
        // run all diagnostics
        timer[3].restart();
        diags.runAllDiags(itime, EMfields, vecSpecies, Interp, smpi);
        
        // temporary EM fields dump in Fields.h5
        if  ((diags.params.fieldDump_every != 0) && (itime % diags.params.fieldDump_every == 0))
            sio->writeAllFieldsSingleFileTime( EMfields, itime );
        
        // temporary EM fields dump in Fields.h5
        if  ((diags.params.avgfieldDump_every != 0) && (itime % diags.params.avgfieldDump_every == 0))
            sio->writeAvgFieldsSingleFileTime( EMfields, itime );
        
        // temporary particles dump (1 HDF5 file per process)
        if  ((diags.params.particleDump_every != 0) && (itime % diags.params.particleDump_every == 0))
            sio->writePlasma( vecSpecies, time_dual, smpi );
        
        if (sio->dump(EMfields, itime,  vecSpecies, smpi, simWindow, params, input_data)) break;
        
        timer[3].update();
		
        if ( simWindow && simWindow->isMoving(time_dual) ) {
            simWindow->operate(vecSpecies, EMfields, Interp, Proj, smpi );
        }
        
    }//END of the time loop
    
    smpi->barrier();
    
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    if ( smpi->isMaster() ) {
        MESSAGE("End time loop, time dual = " << time_dual);
        MESSAGE("-----------------------------------------------------------------------------------------------------");
    }
    
    //double timElapsed=smpiData->time_seconds();
    //if ( smpi->isMaster() ) MESSAGE(0, "Time in time loop : " << timElapsed );
    timer[0].update();
    MESSAGE(0, "Time in time loop : " << timer[0].getTime() );
    if ( smpi->isMaster() )
        for (int i=1 ; i<ntimer ; i++) timer[i].print();
    
    double coverage(0.);
    for (int i=1 ; i<ntimer ; i++) coverage += timer[i].getTime();
    if ( smpi->isMaster() ) MESSAGE(0, "\t" << setw(12) << "Coverage\t" << coverage/timer[0].getTime()*100. << " %" );
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    // temporary EM fields dump in Fields.h5
    if  ( (diags.params.fieldDump_every != 0) && (params.n_time % diags.params.fieldDump_every != 0) )
        sio->writeAllFieldsSingleFileTime( EMfields, params.n_time );
    
    // temporary time-averaged EM fields dump in Fields_avg.h5
    if  ( (diags.params.avgfieldDump_every != 0) && (params.n_time % diags.params.avgfieldDump_every != 0) )
        sio->writeAvgFieldsSingleFileTime( EMfields, params.n_time );
    
    // temporary particles dump (1 HDF5 file per process)
    if  ( (diags.params.particleDump_every != 0) && (params.n_time % diags.params.particleDump_every != 0) )
        sio->writePlasma( vecSpecies, time_dual, smpi );
    
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
        MESSAGE("-----------------------------------------------------------------------------------------------------");
        MESSAGE("END " << namelist);
        MESSAGE("-----------------------------------------------------------------------------------------------------");
    }
    delete smpi;
    delete smpiData;
    return 0;
    
}//END MAIN


// Printing starting message
void startingMessage(std::string inputfile) {
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    MESSAGE(" Version  : " << __VERSION DEBUGEXEC(<< " DEBUG") << " Compiled : " << __DATE__ << " " << __TIME__);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    MESSAGE(" Namelist : " << inputfile);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
}


