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

#include "Smilei.h"

#include <ctime>
#include <cstdlib>
#include <unistd.h>

#include <iostream>
#include <iomanip>

#include "InputData.h"
#include "PicParams.h"
#include "LaserParams.h"

#include "SmileiMPIFactory.h"
#include "SmileiIOFactory.h"

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "DiagParams.h"
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

    // Check for namelist (input file)
    if (argc<2) ERROR("No namelists given!");
    string namelist=argv[1];
    
    // Send information on current simulation
    
    MESSAGE("                   _            __     ");
    MESSAGE(" ___           _  | |        _  \\ \\    ");
    MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |   Version  :  " << __VERSION);
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |   Compiled :  " << __DATE__ << " " << __TIME__);
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |   Namelist :  " << namelist);
    MESSAGE("                                /_/    ");
    
    // Parse the namelist file (no check!)
    InputData input_data;
    if ( smpiData->isMaster() ) input_data.readFile(namelist);    

    // broadcast file and parse it and randomize
    smpiData->bcast(input_data);    
    
    MESSAGE("----------------------------------------------");
    MESSAGE("Input data info");
    MESSAGE("----------------------------------------------");
    // Read simulation & diagnostics parameters
    PicParams params(input_data);
    smpiData->init(params);
    smpiData->barrier();
    if ( smpiData->isMaster() ) params.print();
    smpiData->barrier();
    LaserParams laser_params(params, input_data);
    smpiData->barrier();
    DiagParams diag_params(params, input_data);
    
    
    // Geometry known, MPI environment specified
    MESSAGE("----------------------------------------------");
    MESSAGE("Creating MPI & IO environments");
    MESSAGE("----------------------------------------------");
    SmileiMPI* smpi = SmileiMPIFactory::create(params, smpiData);
    SmileiIO*  sio  = SmileiIOFactory::create(params, diag_params, smpi);
#ifdef _OMP
    int nthds(0);
#pragma omp parallel shared(nthds)
    {
        nthds = omp_get_num_threads();
    }
    if (smpi->isMaster())
        MESSAGE("\tOpenMP : Number of thread per MPI process : " << nthds );
#else
    if (smpi->isMaster()) MESSAGE("\tOpenMP : Disabled");
#endif
    
    
    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------

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
    int start_moving(0);
    if (params.nspace_win_x)
        simWindow = new SimWindow(params);

    MESSAGE("----------------------------------------------");
    MESSAGE("Creating EMfields/Interp/Proj/Diags");
    MESSAGE("----------------------------------------------");
    
    // Initialize the electromagnetic fields and interpolation-projection operators
    // according to the simulation geometry
    // ----------------------------------------------------------------------------

    // object containing the electromagnetic fields (virtual)
    ElectroMagn* EMfields = ElectroMagnFactory::create(params, laser_params, smpi);
    
    // interpolation operator (virtual)
    Interpolator* Interp = InterpolatorFactory::create(params, smpi);
    
    // projection operator (virtual)
    Projector* Proj = ProjectorFactory::create(params, smpi);
    
    // Create diagnostics
    Diagnostic *Diags =new Diagnostic(params,diag_params, smpi);    
    
    smpi->barrier();
    
    unsigned int stepStart=0, stepStop=params.n_time;
    
    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        DEBUG(vecSpecies.size());
        sio->restartAll( EMfields,  stepStart, vecSpecies, smpi, simWindow, params, input_data);

        double restart_time_dual = (stepStart +0.5) * params.timestep;
	if ( simWindow && ( simWindow->isMoving(restart_time_dual) ) ) {
	    simWindow->setOperators(vecSpecies, Interp, Proj, smpi);
	    simWindow->operate(vecSpecies, EMfields, Interp, Proj, smpi , params);
	}
	    
    } else {
        // Initialize the electromagnetic fields
        // -----------------------------------
        // Init rho and J by projecting all particles of subdomain
        EMfields->initRhoJ(vecSpecies, Proj);
        
        // Sum rho and J on ghost domains
        smpi->sumRhoJ( EMfields );
        for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
            smpi->sumRhoJs(EMfields, ispec, true);
        }
        
        // Init electric field (Ex/1D, + Ey/2D)
        MESSAGE("----------------------------------------------");
        MESSAGE("Solving Poisson at time t = 0");
        MESSAGE("----------------------------------------------");
	if (!EMfields->isRhoNull(smpi)) 
	    EMfields->solvePoisson(smpi);
        
        
        MESSAGE("----------------------------------------------");
        MESSAGE("Running diags at time t = 0");
        MESSAGE("----------------------------------------------");
        // run diagnostics at time-step 0
        Diags->runAllDiags(0, EMfields, vecSpecies, Interp, smpi);
        // temporary EM fields dump in Fields.h5
        sio->writeAllFieldsSingleFileTime( EMfields, 0 );
        // temporary EM fields dump in Fields_avg.h5
        if (diag_params.ntime_step_avg!=0)
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
    int ntimer(6);
    Timer timer[ntimer];
    timer[0].init(smpi, "global");
    timer[1].init(smpi, "particles");
    timer[2].init(smpi, "maxwell");
    timer[3].init(smpi, "diagnostics");
    timer[4].init(smpi, "densities");
    timer[5].init(smpi, "Mov window");
    
    
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
	if ( (itime % diag_params.print_every == 0) &&  ( smpi->isMaster() ) )
            MESSAGE(1,"t = "          << setw(7) << setprecision(2)   << time_dual/params.conv_fac
                    << "   it = "       << setw(log10(params.n_time)+1) << itime  << "/" << params.n_time
                    << "   sec = "      << setw(7) << setprecision(2)   << timer[0].getTime()
                    << "   E = "        << std::scientific << setprecision(4)<< Diags->getScalar("Etot")
                    << "   Epart = "        << std::scientific << setprecision(4)<< Diags->getScalar("Eparticles")
                    << "   Elost = "        << std::scientific << setprecision(4)<< Diags->getScalar("Elost")
                    << "   E_bal(%) = " << setw(6) << std::fixed << setprecision(2)   << 100.0*Diags->getScalar("Ebal_norm") );

        // put density and currents to 0 + save former density
        // ---------------------------------------------------
        EMfields->restartRhoJ();
        
        
        // apply the PIC method
        // --------------------
        // for all particles of all species (see dynamic in Species.cpp)
        // (1) interpolate the fields at the particle position
        // (2) move the particle
        // (3) calculate the currents (charge conserving method)
        timer[1].restart();
#pragma omp parallel shared (EMfields,time_dual,vecSpecies,smpi,params)
        {
            int tid(0);
#ifdef _OMP
            tid = omp_get_thread_num();
#endif
            for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
                if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ){
                    EMfields->restartRhoJs(ispec, time_dual > params.species_param[ispec].time_frozen);
                    vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, smpi, params, simWindow);
                }
            }
            for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
#pragma omp barrier
                if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ){
		    // Loop on dims to manage exchange in corners
		    for ( int iDim = 0 ; iDim<params.nDim_particle ; iDim++ )
			smpi->exchangeParticles(vecSpecies[ispec], ispec, params, tid, iDim);
#pragma omp barrier
                        vecSpecies[ispec]->sort_part();
                }
            }
        }
        timer[1].update();
        
	//!\todo To simplify : sum global and per species densities
        timer[4].restart();
        smpi->sumRhoJ( EMfields );
        for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
            if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ) smpi->sumRhoJs(EMfields, ispec, time_dual > params.species_param[ispec].time_frozen);
        }
        EMfields->computeTotalRhoJ();
        timer[4].update();
        
        // solve Maxwell's equations
        timer[2].restart();
        EMfields->solveMaxwell(itime, time_dual, smpi, params, simWindow);
        timer[2].update();
        
        // incrementing averaged electromagnetic fields
        if (diag_params.ntime_step_avg) EMfields->incrementAvgFields(itime, diag_params.ntime_step_avg);
        
        // call the various diagnostics
        // ----------------------------
		
        // run all diagnostics
        timer[3].restart();
        Diags->runAllDiags(itime, EMfields, vecSpecies, Interp, smpi);
        
        // temporary EM fields dump in Fields.h5
        if  ((diag_params.fieldDump_every != 0) && (itime % diag_params.fieldDump_every == 0))
            sio->writeAllFieldsSingleFileTime( EMfields, itime );
        
        // temporary EM fields dump in Fields.h5
        if  (diag_params.ntime_step_avg!=0)
            if ((diag_params.avgfieldDump_every != 0) && (itime % diag_params.avgfieldDump_every == 0))
                sio->writeAvgFieldsSingleFileTime( EMfields, itime );
        
#ifdef _IO_PARTICLE
        // temporary particles dump (1 HDF5 file per process)
        if  ((diag_params.particleDump_every != 0) && (itime % diag_params.particleDump_every == 0))
            sio->writePlasma( vecSpecies, time_dual, smpi );
#endif
        
        if (sio->dump(EMfields, itime,  vecSpecies, smpi, simWindow, params, input_data)) break;
        
        timer[3].update();
		
        timer[5].restart();
        if ( simWindow && simWindow->isMoving(time_dual) ) {
            start_moving++;
            if ((start_moving==1) && (smpi->isMaster()) ) {
		MESSAGE(">>> Window starts moving");
            }
            simWindow->operate(vecSpecies, EMfields, Interp, Proj, smpi, params);
        }
        timer[5].update();
        
    }//END of the time loop
    
    smpi->barrier();
    
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    MESSAGE("End time loop, time dual = " << time_dual/params.conv_fac);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    
    //double timElapsed=smpiData->time_seconds();
    //if ( smpi->isMaster() ) MESSAGE(0, "Time in time loop : " << timElapsed );
    timer[0].update();
    MESSAGE(0, "Time in time loop : " << timer[0].getTime() );
    if ( smpi->isMaster() )
        for (int i=1 ; i<ntimer ; i++) timer[i].print(timer[0].getTime());
    
    double coverage(0.);
    for (int i=1 ; i<ntimer ; i++) coverage += timer[i].getTime();
    MESSAGE(0, "\t" << setw(12) << "Coverage\t" << coverage/timer[0].getTime()*100. << " %" );
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    // temporary EM fields dump in Fields.h5
    if  ( (diag_params.fieldDump_every != 0) && (params.n_time % diag_params.fieldDump_every != 0) )
        sio->writeAllFieldsSingleFileTime( EMfields, params.n_time );
    // temporary time-averaged EM fields dump in Fields_avg.h5
    if  (diag_params.ntime_step_avg!=0)
        if  ( (diag_params.avgfieldDump_every != 0) && (params.n_time % diag_params.avgfieldDump_every != 0) )
            sio->writeAvgFieldsSingleFileTime( EMfields, params.n_time );
#ifdef _IO_PARTICLE
    // temporary particles dump (1 HDF5 file per process)
    if  ( (diag_params.particleDump_every != 0) && (params.n_time % diag_params.particleDump_every != 0) )
        sio->writePlasma( vecSpecies, time_dual, smpi );
#endif    

    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    delete Proj;
    delete Interp;
    delete EMfields;
    delete Diags;
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    MESSAGE("END " << namelist);
    MESSAGE("-----------------------------------------------------------------------------------------------------");

    delete sio;
    delete smpi;
    delete smpiData;
    if (params.nspace_win_x)
        delete simWindow;
    
    return 0;
    
}//END MAIN



