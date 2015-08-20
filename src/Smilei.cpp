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

#include "Params.h"

#include "SmileiMPIFactory.h"
#include "SmileiIOFactory.h"

#include "SpeciesFactory.h"
#include "Collisions.h"
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
    cout.setf( ios::fixed,  ios::floatfield ); // floatfield set to fixed
    
    // Define 2 MPI environments :
    //  - smpiData : to broadcast input data, unknown geometry
    //  - smpi (defined later) : to compute/exchange data, specific to a geometry
    SmileiMPI *smpiData= new SmileiMPI(&argc, &argv );
    
    // -------------------------
    // Simulation Initialization
    // ------------------------- 
    
    // Check for namelists (input files)
    vector<string> namelists(argv + 1, argv + argc);
    if (namelists.size()==0) ERROR("No namelists given!");
    
    // Send information on current simulation
    
    MESSAGE("                   _            _");
    MESSAGE(" ___           _  | |        _  \\ \\    ");
    MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |");
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |   Version  : " << __VERSION);
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |   Date     : " << __COMMITDATE);
    MESSAGE("                                /_/    ");
    
    TITLE("Input data info");
    
    // Read simulation & diagnostics parameters
    Params params(smpiData,namelists);
    smpiData->init(params);
    smpiData->barrier();
    if ( smpiData->isMaster() ) params.print();
    smpiData->barrier();
    
    
    // Geometry known, MPI environment specified
    TITLE("General MPI environement");
    SmileiMPI* smpi = SmileiMPIFactory::create(params, smpiData);
    
    // Create diagnostics
    TITLE("Creating Diagnostics");
    Diagnostic Diags(params, smpi);
    
    //Create mpi i/o environment
    TITLE("MPI input output environment");
    SmileiIO*  sio  = SmileiIOFactory::create(params, Diags, smpi);
    
    // setup OpenMP
    TITLE("OpenMP");
#ifdef _OMP
    int nthds(0);
#pragma omp parallel shared(nthds)
    {
        nthds = omp_get_num_threads();
    }
    if (smpi->isMaster())
        MESSAGE(1,"Number of thread per MPI process : " << nthds );
#else
    if (smpi->isMaster()) MESSAGE(1,"Disabled");
#endif
    
    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------
    
    // ---------------------------
    // Initialize Species & Fields
    // ---------------------------
    TITLE("Initializing particles, fields & moving-window");
    
    // Initialize the vecSpecies object containing all information of the different Species
    // ------------------------------------------------------------------------------------
    
    // vector of Species (virtual)
    vector<Species*> vecSpecies = SpeciesFactory::createVector(params, smpi);
    
    // Initialize the collisions (vector of collisions)
    // ------------------------------------------------------------------------------------
    vector<Collisions*> vecCollisions = Collisions::create(params, vecSpecies, smpi);
    
    // ----------------------------------------------------------------------------
    // Define Moving Window & restart
    // ----------------------------------------------------------------------------
    
    SimWindow* simWindow = NULL;
    int start_moving(0);
    if (params.nspace_win_x)
        simWindow = new SimWindow(params);
    
    TITLE("Creating EMfields/Interp/Proj");
    
    // Initialize the electromagnetic fields and interpolation-projection operators
    // according to the simulation geometry
    // ----------------------------------------------------------------------------
    
    // object containing the electromagnetic fields (virtual)
    ElectroMagn* EMfields = ElectroMagnFactory::create(params, smpi);
    
    // interpolation operator (virtual)
    Interpolator* Interp = InterpolatorFactory::create(params, smpi);
    
    // projection operator (virtual)
    Projector* Proj = ProjectorFactory::create(params, smpi);
    
    smpi->barrier();
    
    unsigned int stepStart=0, stepStop=params.n_time;
    
    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        DEBUG(vecSpecies.size());
        sio->restartAll( EMfields,  stepStart, vecSpecies, smpi, simWindow, params);
        
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
        for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++) {
            smpi->sumRhoJs(EMfields, ispec, true);  // only if !isTestParticles
        }
        
        if (!EMfields->isRhoNull(smpi))  {
            // Init electric field (Ex/1D, + Ey/2D)
            TITLE("Solving Poisson at time t = 0");
            Timer ptimer;
            ptimer.init(smpi, "global");
            ptimer.restart();
            EMfields->solvePoisson(smpi);
            ptimer.update();
            MESSAGE("Time in Poisson : " << ptimer.getTime() );
        }
        
        TITLE("Applying external fields at time t = 0");
        EMfields->applyExternalFields(smpi);
        
        TITLE("Running diags at time t = 0");
        Diags.runAllDiags(0, EMfields, vecSpecies, Interp, smpi);        
        // temporary EM fields dump in Fields.h5
        sio->writeAllFieldsSingleFileTime( &(EMfields->allFields), 0, 0 );
        // temporary EM fields dump in Fields_avg.h5
        if (Diags.ntime_step_avg!=0)
            sio->writeAllFieldsSingleFileTime( &(EMfields->allFields_avg), 0, 1 );
        // temporary particle dump at time 0
        sio->writePlasma( vecSpecies, 0., smpi );
        
        for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
            if ( (vecSpecies[ispec]->particles.isTestParticles) ) {
                sio->initWriteTestParticles(vecSpecies[ispec], ispec, 0, params, smpi);
                sio->writeTestParticles(vecSpecies[ispec], ispec, 0, params, smpi);
                //MPI_Finalize();
                //return 0;
            }
        }
        
    }
    
    
    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    params.cleanup();
    
    
    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // ------------------------------------------------------------------------
    
    // time at integer time-steps (primal grid)
    double time_prim = stepStart * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (stepStart +0.5) * params.timestep;
    
    // Count timer
    vector<Timer> timer(9);
    
    timer[0].init(smpi, "Global");
    timer[1].init(smpi, "Particles");
    timer[2].init(smpi, "Maxwell");
    timer[3].init(smpi, "Diagnostics");
    timer[4].init(smpi, "Densities");
    timer[5].init(smpi, "Mov window");
    timer[6].init(smpi, "Fields");
    timer[7].init(smpi, "AvgFields");
    timer[8].init(smpi, "Collisions");
    
    
    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE("Time-Loop is started: number of time-steps n_time = " << params.n_time);
    
    for (unsigned int itime=stepStart+1 ; itime <= stepStop ; itime++) {
        
        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        // send message at given time-steps
        // --------------------------------
        timer[0].update();
        
        if ( (itime % Diags.print_every == 0) &&  ( smpi->isMaster() ) ) {
            
            MESSAGE(1,"t = "          << setw(7) << setprecision(2)   << time_dual
                    << "   it = "     << setw(log10(params.n_time)+1) << itime  << "/" << params.n_time
                    << "   sec = "    << setw(7) << setprecision(2)   << timer[0].getTime()
                    << "   Utot = "   << scientific << setprecision(4)<< Diags.getScalar("Utot")
                    << "   Uelm = "   << scientific << setprecision(4)<< Diags.getScalar("Uelm")
                    << "   Ukin = "   << scientific << setprecision(4)<< Diags.getScalar("Ukin")
                    << "   Ubal(%) = "<< setw(6) << fixed << setprecision(2) << 100.0*Diags.getScalar("Ubal_norm")
                    );
            
            //!\todo (MG to JD/AD) We should clear this. Either not print it OR add it to the former stream
            if (simWindow) {
                double Uinj_mvw = Diags.getScalar("Uelm_inj_mvw") + Diags.getScalar("Ukin_inj_mvw");
                double Uout_mvw = Diags.getScalar("Uelm_out_mvw") + Diags.getScalar("Ukin_out_mvw");
                MESSAGE(1, "\t\t Uinj_mvw = " << scientific << setprecision(4) << Uinj_mvw
                        << "     Uout_mvw = " << scientific << setprecision(4) << Uout_mvw
                        );
            }//simWindow
            
        }//itime
        
        
        // put density and currents to 0 + save former density
        // ---------------------------------------------------
        EMfields->restartRhoJ();
        
        
        timer[8].restart();
        // apply collisions if requested
        // -----------------------------
        if (Collisions::debye_length_required)
            Collisions::calculate_debye_length(params,vecSpecies);
        for (unsigned int icoll=0 ; icoll<vecCollisions.size(); icoll++)
            vecCollisions[icoll]->collide(params,vecSpecies,itime);
        timer[8].update();
        
        
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
            for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++) {
                if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ){
                    EMfields->restartRhoJs(ispec, time_dual > params.species_param[ispec].time_frozen); // if (!isTestParticles)
                    vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, smpi, params, simWindow);
                }
            }
            for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++) {
#pragma omp barrier
                if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ){
                    // Loop on dims to manage exchange in corners
                    for ( int iDim = 0 ; iDim<(int)params.nDim_particle ; iDim++ )
                        smpi->exchangeParticles(vecSpecies[ispec], ispec, params, tid, iDim);
#pragma omp barrier
                        vecSpecies[ispec]->sort_part(); // Should we sort test particles ?? (JD)
                }
            }
        }
        timer[1].update();
        
        
        //!\todo To simplify : sum global and per species densities
        timer[4].restart();
        smpi->sumRhoJ( EMfields );
        for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++) {
            if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ) smpi->sumRhoJs(EMfields, ispec, time_dual > params.species_param[ispec].time_frozen);
        }
        EMfields->computeTotalRhoJ();
        timer[4].update();
        
        // solve Maxwell's equations
        if( time_dual > params.time_fields_frozen ) {
            timer[2].restart();
            EMfields->solveMaxwell(itime, time_dual, smpi, params, simWindow);
            timer[2].update();
        }
        
        // incrementing averaged electromagnetic fields
        if (Diags.ntime_step_avg) EMfields->incrementAvgFields(itime, Diags.ntime_step_avg);
        
        
        // call the various diagnostics
        // ----------------------------
        
        // run all diagnostics
        timer[3].restart();
        for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
            if ( (vecSpecies[ispec]->particles.isTestParticles)  )
                sio->writeTestParticles(vecSpecies[ispec], ispec, itime, params, smpi);
        }
        
        
        Diags.runAllDiags(itime, EMfields, vecSpecies, Interp, smpi);
        timer[3].update();
        
        timer[6].restart();
        // temporary EM fields dump in Fields.h5
        if  ((Diags.fieldDump_every != 0) && (itime % Diags.fieldDump_every == 0))
            sio->writeAllFieldsSingleFileTime( &(EMfields->allFields), itime, 0 );
        timer[6].update();
        
        timer[7].restart();
        // temporary EM fields dump in Fields_avg.h5
        if  ((Diags.ntime_step_avg!=0) &&
             (Diags.avgfieldDump_every != 0) && 
             (itime % Diags.avgfieldDump_every == 0)) {
            sio->writeAllFieldsSingleFileTime( &(EMfields->allFields_avg), itime, 1 );
        }
        timer[7].update();
        
#ifdef _IO_PARTICLE
        // temporary particles dump (1 HDF5 file per process)
        if  ((Diags.particleDump_every != 0) && (itime % Diags.particleDump_every == 0))
            sio->writePlasma( vecSpecies, time_dual, smpi );
#endif
        
        if (sio->dump(EMfields, itime, vecSpecies, smpi, simWindow, params)) break;
        
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
    MESSAGE("End time loop, time dual = " << time_dual);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    
    //double timElapsed=smpiData->time_seconds();
    //if ( smpi->isMaster() ) MESSAGE("Time in time loop : " << timElapsed );
    timer[0].update();
    TITLE("Time profiling :");
    
    double coverage(0.);
    for (unsigned int i=1 ; i<timer.size() ; i++) coverage += timer[i].getTime();
    MESSAGE("Time in time loop :\t" << timer[0].getTime() << "\t("<<coverage/timer[0].getTime()*100.<< "% coverage)" );
    if ( smpi->isMaster() )
        for (unsigned int i=1 ; i<timer.size() ; i++) timer[i].print(timer[0].getTime());
    
    Diags.printTimers(smpi,timer[3].getTime());
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    // temporary EM fields dump in Fields.h5
    if  ( (Diags.fieldDump_every != 0) && (params.n_time % Diags.fieldDump_every != 0) )
        sio->writeAllFieldsSingleFileTime( &(EMfields->allFields), params.n_time, 0 );
    // temporary time-averaged EM fields dump in Fields_avg.h5
    if  (Diags.ntime_step_avg!=0)
        if  ( (Diags.avgfieldDump_every != 0) && (params.n_time % Diags.avgfieldDump_every != 0) )
            sio->writeAllFieldsSingleFileTime( &(EMfields->allFields_avg), params.n_time, 1 );
#ifdef _IO_PARTICLE
    // temporary particles dump (1 HDF5 file per process)
    if  ( (Diags.particleDump_every != 0) && (params.n_time % Diags.particleDump_every != 0) )
        sio->writePlasma( vecSpecies, time_dual, smpi );
#endif    
    
    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    delete Proj;
    delete Interp;
    delete EMfields;
    for(unsigned int i=0; i<vecCollisions.size(); i++) delete vecCollisions[i];
    vecCollisions.clear();
    
    Diags.closeAll(smpi);
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    
    if (params.nspace_win_x)
        delete simWindow;
    
    PyTools::closePython();

    TITLE("END");
    delete sio;
    delete smpi;
    delete smpiData;
    
    return 0;
    
}//END MAIN



