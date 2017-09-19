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

#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <omp.h>

#include "Smilei.h"
#include "SmileiMPI_test.h"
#include "Params.h"
#include "PatchesFactory.h"
#include "Checkpoint.h"
#include "Solver.h"
#include "SimWindow.h"
#include "Diagnostic.h"
#include "Timers.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    cout.setf( ios::fixed,  ios::floatfield ); // floatfield set to fixed
    
    // -------------------------
    // Simulation Initialization
    // ------------------------- 
    
    // Create MPI environment :
#ifdef SMILEI_TESTMODE
    SmileiMPI_test * smpi = new SmileiMPI_test( &argc, &argv );
#else
    SmileiMPI * smpi = new SmileiMPI( &argc, &argv );
#endif
    
    MESSAGE("                   _            _");
    MESSAGE(" ___           _  | |        _  \\ \\   Version : " << __VERSION);
    MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |   ");
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |");
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |  ");
    MESSAGE("                                /_/    ");
    MESSAGE("");
    
    // Read and print simulation parameters
    TITLE("Reading the simulation parameters");
    Params params(smpi,vector<string>(argv + 1, argv + argc));
    OpenPMDparams openPMD(params);
    
    // Initialize MPI environment with simulation parameters
    TITLE("Initializing MPI");
    smpi->init(params);
    
    // Create timers
    Timers timers(smpi);
    
    // Print in stdout MPI, OpenMP, patchs parameters
    params.print_parallelism_params(smpi);
    
    TITLE("Initializing the restart environment");
    Checkpoint checkpoint(params, smpi);
    
    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // Update in "if restart" if necessary
    // ------------------------------------------------------------------------
    
    // time at integer time-steps (primal grid)
    double time_prim = 0;
    // time at half-integer time-steps (dual grid)
    double time_dual = 0.5 * params.timestep;
    
    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------
    // --------------------
    // Define Moving Window
    // --------------------
    TITLE("Initializing moving window");
    SimWindow* simWindow = new SimWindow(params);
    
    // ---------------------------------------------------
    // Initialize patches (including particles and fields)
    // ---------------------------------------------------
    TITLE("Initializing particles & fields");
    VectorPatch vecPatches;
    
    // reading from dumped file the restart values
    if (params.restart) {
        
        // smpi->patch_count recomputed in restartAll
        // vecPatches allocated in restartAll according to patch_count saved
        checkpoint.restartAll( vecPatches, smpi, simWindow, params, openPMD);
        
        // time at integer time-steps (primal grid)
        time_prim = checkpoint.this_run_start_step * params.timestep;
        // time at half-integer time-steps (dual grid)
        time_dual = (checkpoint.this_run_start_step +0.5) * params.timestep;
        
        TITLE("Initializing diagnostics");
        vecPatches.initAllDiags( params, smpi );
        
    } else {
       
        vecPatches = PatchesFactory::createVector(params, smpi, openPMD, 0);
        
        if( ! smpi->test_mode ) {
            
            // Initialize the electromagnetic fields
            // -------------------------------------
            vecPatches.computeCharge();
                vecPatches.sumDensities(params, time_dual, timers, 0, simWindow);
            
            // Apply antennas
            // --------------
            vecPatches.applyAntennas(0.5 * params.timestep);
            
            // Init electric field (Ex/1D, + Ey/2D)
            if (!vecPatches.isRhoNull(smpi) && params.solve_poisson == true) {
                TITLE("Solving Poisson at time t = 0");
                Timer ptimer("global");
                ptimer.init(smpi);
                ptimer.restart();
                
                vecPatches.solvePoisson( params, smpi );
                ptimer.update();
                MESSAGE("Time in Poisson : " << ptimer.getTime() );
            }
            
            vecPatches.dynamics(params, smpi, simWindow, time_dual, timers, 0);
            timers.particles.reboot();
            timers.syncPart .reboot();
            
            vecPatches.sumDensities(params, time_dual, timers, 0, simWindow );
            timers.densities.reboot();
            timers.syncDens .reboot();
            
            TITLE("Applying external fields at time t = 0");
            vecPatches.applyExternalFields();
            
            vecPatches.finalize_and_sort_parts(params, smpi, simWindow, time_dual, timers, 0);
            timers.syncPart .reboot();
            
            TITLE("Initializing diagnostics");
            vecPatches.initAllDiags( params, smpi );
            TITLE("Running diags at time t = 0");
            vecPatches.runAllDiags(params, smpi, 0, timers, simWindow);
            timers.diags.reboot();
        }
    }
    
    TITLE("Species creation summary");
    vecPatches.printNumberOfParticles( smpi );
    
    timers.global.reboot();
    
    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    params.cleanup(smpi);
    
    // ------------------------------------------------------------------------
    // Check memory consumption
    // ------------------------------------------------------------------------
    TITLE("Memory consumption");
    vecPatches.check_memory_consumption( smpi );
    
    
    // If test mode enable, code stops here
    if( smpi->test_mode ) {
        delete simWindow;
        PyTools::closePython();
        TITLE("END TEST MODE");
        delete smpi;
        return 0;
    }
    
/*tommaso
    // save latestTimeStep (used to test if we are at the latest timestep when running diagnostics at run's end)
    unsigned int latestTimeStep=checkpoint.this_run_start_step;
*/
    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    
    TITLE("Time-Loop started: number of time-steps n_time = " << params.n_time);
    if ( smpi->isMaster() ) params.print_timestep_headers();
    
    #pragma omp parallel shared (time_dual,smpi,params, vecPatches, simWindow, checkpoint)
    {
        
        unsigned int itime=checkpoint.this_run_start_step+1;
        while ( (itime <= params.n_time) && (!checkpoint.exit_asap) ) {
            
            // calculate new times
            // -------------------
            #pragma omp single
            {
                time_prim += params.timestep;
                time_dual += params.timestep;
            }
            
            // apply collisions if requested
            vecPatches.applyCollisions(params, itime, timers);
            
            // (1) interpolate the fields at the particle position
            // (2) move the particle
            // (3) calculate the currents (charge conserving method)
            vecPatches.dynamics(params, smpi, simWindow, time_dual, timers, itime);
            
            // Sum densities
            vecPatches.sumDensities(params, time_dual, timers, itime, simWindow );
            
            // apply currents from antennas
            vecPatches.applyAntennas(time_dual);
            
            // solve Maxwell's equations
            if( time_dual > params.time_fields_frozen )
                vecPatches.solveMaxwell( params, simWindow, itime, time_dual, timers );

            vecPatches.finalize_and_sort_parts(params, smpi, simWindow, time_dual, timers, itime);

            // call the various diagnostics
            vecPatches.runAllDiags(params, smpi, itime, timers, simWindow);
            
            timers.movWindow.restart();
            simWindow->operate(vecPatches, smpi, params, itime, time_dual);
            timers.movWindow.update();
            
            // ----------------------------------------------------------------------
            // Validate restart  : to do
            // Restart patched moving window : to do
            #pragma omp master
            checkpoint.dump(vecPatches, itime, smpi, simWindow, params);
            #pragma omp barrier
            // ----------------------------------------------------------------------        
            
        
            if ((params.balancing_every > 0) && (smpi->getSize()!=1) ) {
                if (( itime%params.balancing_every == 0 )) {
                    timers.loadBal.restart();
                    #pragma omp single
                    vecPatches.load_balance( params, time_dual, smpi, simWindow, itime );
                    timers.loadBal.update( params.printNow( itime ) );
                }
            }
        
            // print message at given time-steps
            // --------------------------------
            if ( smpi->isMaster() &&  params.printNow( itime ) )
                params.print_timestep(itime, time_dual, timers.global); //contain a timer.update !!!

            itime++;
            
        }//END of the time loop

    } //End omp parallel region

    
    smpi->barrier();
    
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE("End time loop, time dual = " << time_dual);
    timers.global.update();
    
    TITLE("Time profiling : (print time > 0.001%)");
    timers.profile(smpi);
    
/*tommaso
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    if (latestTimeStep==params.n_time)
        vecPatches.runAllDiags(params, smpi, &diag_flag, params.n_time, timer, simWindow);
*/
    
    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    vecPatches.close( smpi );
    smpi->barrier(); // Don't know why but sync needed by HDF5 Phasespace managment
    delete simWindow;
    PyTools::closePython();
    TITLE("END");
    delete smpi;
    
    return 0;
    
}//END MAIN

// ---------------------------------------------------------------------------------------------------------------------
//                                               END MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
