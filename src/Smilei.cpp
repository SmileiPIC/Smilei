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
#include "SyncVectorPatch.h"
#include "Checkpoint.h"
#include "Solver.h"
#include "SimWindow.h"
#include "Diagnostic.h"
#include "Region.h"
#include "DoubleGrids.h"
#include "DoubleGridsAM.h"
#include "Timers.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
    cout.setf( ios::fixed,  ios::floatfield ); // floatfield set to fixed

    // -------------------------
    // Simulation Initialization
    // -------------------------

    // Create MPI environment :

#ifdef SMILEI_TESTMODE
    SmileiMPI_test smpi( &argc, &argv );
#else
    SmileiMPI smpi( &argc, &argv );
#endif

    MESSAGE( "                   _            _" );
    MESSAGE( " ___           _  | |        _  \\ \\   Version : " << __VERSION );
    MESSAGE( "/ __|  _ __   (_) | |  ___  (_)  | |   " );
    MESSAGE( "\\__ \\ | '  \\   _  | | / -_)  _   | |" );
    MESSAGE( "|___/ |_|_|_| |_| |_| \\___| |_|  | |  " );
    MESSAGE( "                                /_/    " );
    MESSAGE( "" );

    // Read and print simulation parameters
    TITLE( "Reading the simulation parameters" );
    Params params( &smpi, vector<string>( argv + 1, argv + argc ) );
    OpenPMDparams openPMD( params );
    PyTools::setIteration( 0 );

    // Need to move it here because of domain decomposition need in smpi->init(_patch_count)
    //     abstraction of Hilbert curve
    VectorPatch vecPatches( params );
    Region region( params );

    // Initialize MPI environment with simulation parameters
    TITLE( "Initializing MPI" );
    smpi.init( params, vecPatches.domain_decomposition_ );

    // Create timers
    Timers timers( &smpi );

    // Print in stdout MPI, OpenMP, patchs parameters
    params.print_parallelism_params( &smpi );

    TITLE( "Initializing the restart environment" );
    Checkpoint checkpoint( params, &smpi );

    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // Update in "if restart" if necessary
    // ------------------------------------------------------------------------

    // time at integer time-steps (primal grid)
    double time_prim = 0;
    // time at half-integer time-steps (dual grid)
    double time_dual = 0.5 * params.timestep;

    // --------------------
    // Define Moving Window
    // --------------------
    SimWindow *simWindow = new SimWindow( params );

    // ------------------------------------------------------------------------
    // Init nonlinear inverse Compton scattering
    // ------------------------------------------------------------------------
    RadiationTables radiation_tables_;

    // ------------------------------------------------------------------------
    // Create MultiphotonBreitWheelerTables object for multiphoton
    // Breit-Wheeler pair creation
    // ------------------------------------------------------------------------
    MultiphotonBreitWheelerTables MultiphotonBreitWheelerTables;

    // ---------------------------------------------------
    // Special test mode
    // ---------------------------------------------------
    if( smpi.test_mode ) {
        executeTestMode( vecPatches, region, &smpi, simWindow, params, checkpoint, openPMD, &radiation_tables_ );
        return 0;
    }

    // ---------------------------------------------------------------------
    // Init and compute tables for radiation effects
    // (nonlinear inverse Compton scattering)
    // ---------------------------------------------------------------------
    radiation_tables_.initialization( params, &smpi);

    // ---------------------------------------------------------------------
    // Init and compute tables for multiphoton Breit-Wheeler pair creation
    // ---------------------------------------------------------------------
    MultiphotonBreitWheelerTables.initialization( params, &smpi );

    // reading from dumped file the restart values
    if( params.restart ) {
        // smpi.patch_count recomputed in readPatchDistribution
        checkpoint.readPatchDistribution( &smpi, simWindow );
        // allocate patches according to smpi.patch_count
        PatchesFactory::createVector( vecPatches, params, &smpi, openPMD, &radiation_tables_, checkpoint.this_run_start_step+1, simWindow->getNmoved() );
        
        // allocate region according to dump
        if( params.multiple_decomposition ) {
            TITLE( "Create SDMD grids" );
            // read region hindex
            checkpoint.readRegionDistribution( region );

            // Build params.map_rank contains MPI ranks assuming that regions are distributed linearly
            int target_map[smpi.getSize()];
            MPI_Allgather(&(region.vecPatch_.refHindex_), 1, MPI_INT,
                          target_map, 1, MPI_INT,
                          MPI_COMM_WORLD);
            region.define_regions_map(target_map, &smpi, params);

            // params.map_rank used to defined regions neighborood
            region.build( params, &smpi, vecPatches, openPMD, false );
            region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
            region.identify_missing_patches( &smpi, vecPatches, params );
        }
        
        // vecPatches data read in restartAll according to smpi.patch_count

        // if (params.multiple_decomposition) {
        //     region.vecPatch_.refHindex_ = smpi.getRank();
        //     region.build( params, &smpi, vecPatches, openPMD, false );
        //     region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
        //     region.identify_missing_patches( &smpi, vecPatches, params );

        //     region.reset_fitting( &smpi, params );

        //     region.clean();
        //     region.reset_mapping();

        //     region.build( params, &smpi, vecPatches, openPMD, false );
        //     region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
        //     region.identify_missing_patches( &smpi, vecPatches, params );
        // }

        checkpoint.restartAll( vecPatches, region, &smpi, simWindow, params, openPMD );
        vecPatches.sortAllParticles( params );
        
        TITLE( "Minimum memory consumption (does not include all temporary buffers)" );
        vecPatches.checkMemoryConsumption( &smpi, &region.vecPatch_ );
        
        // Patch reconfiguration for the adaptive vectorization
        if( params.has_adaptive_vectorization ) {
            vecPatches.configuration( params, timers, 0 );
        }
        
        // time at integer time-steps (primal grid)
        time_prim = checkpoint.this_run_start_step * params.timestep;
        // time at half-integer time-steps (dual grid)
        time_dual = ( checkpoint.this_run_start_step +0.5 ) * params.timestep;
        
        TITLE( "Open files & initialize diagnostics" );
        vecPatches.initAllDiags( params, &smpi );
        
    } else {
        
        PatchesFactory::createVector( vecPatches, params, &smpi, openPMD, &radiation_tables_, 0 );
        vecPatches.sortAllParticles( params );

        // Create SDMD grids 
        if( params.multiple_decomposition ) {
            TITLE( "Create SDMD grids" );
            region.vecPatch_.refHindex_ = smpi.getRank();
            region.build( params, &smpi, vecPatches, openPMD, false );
            region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
            region.identify_missing_patches( &smpi, vecPatches, params );
            //cout << smpi.getRank() << "\t - local : " << region.local_patches_.size()
            //     <<  "\t - missing : " << region.missing_patches_.size()
            //     <<  "\t - additional : " << region.additional_patches_.size() << endl;
            
            region.reset_fitting( &smpi, params );
            region.clean();
            region.reset_mapping();
            
            region.build( params, &smpi, vecPatches, openPMD, false );
            region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
            region.identify_missing_patches( &smpi, vecPatches, params );
            //cout << smpi.getRank() << "\t - local : " << region.local_patches_.size()
            //     <<  "\t - missing : " << region.missing_patches_.size()
            //     <<  "\t - additional : " << region.additional_patches_.size() << endl;
        }
        
        TITLE( "Minimum memory consumption (does not include all temporary buffers)" );
        vecPatches.checkMemoryConsumption( &smpi, &region.vecPatch_ );
        
        TITLE( "Initial fields setup" );
        
        // Solve "Relativistic Poisson" problem (including proper centering of fields)
        // Note: the mean gamma for initialization will be computed for all the species
        // whose fields are initialized at this iteration
        if( params.solve_relativistic_poisson == true ) {
            MESSAGE( 1, "Solving relativistic Poisson at time t = 0" );
            vecPatches.runRelativisticModule( time_prim, params, &smpi,  timers );
        }
        
        vecPatches.computeCharge();
        vecPatches.sumDensities( params, time_dual, timers, 0, simWindow, &smpi );
        
        // Init electric field (Ex/1D, + Ey/2D)
        if( params.solve_poisson == true && !vecPatches.isRhoNull( &smpi ) ) {
            MESSAGE( 1, "Solving Poisson at time t = 0" );
            vecPatches.runNonRelativisticPoissonModule( params, &smpi,  timers );
        }
        
        MESSAGE( 1, "Applying external fields at time t = 0" );
        vecPatches.applyExternalFields();
        vecPatches.saveExternalFields( params );
        
        MESSAGE( 1, "Applying prescribed fields at time t = 0" );
        vecPatches.applyPrescribedFields( time_prim );
        
        MESSAGE( 1, "Applying antennas at time t = 0" );
        vecPatches.applyAntennas( 0.5 * params.timestep );
        
        // Patch reconfiguration
        if( params.has_adaptive_vectorization ) {
            vecPatches.configuration( params, timers, 0 );
        }
        
        // if Laser Envelope is used, execute particles and envelope sections of ponderomotive loop
        if( params.Laser_Envelope_model ) {
            MESSAGE( 1, "Initialize envelope" );
            vecPatches.initNewEnvelope( params );
        }
        
        // Project charge and current densities (and susceptibility if envelope is used) only for diags at t=0
        vecPatches.projectionForDiags( params, &smpi, simWindow, time_dual, timers, 0 );
        
        // If Laser Envelope is used, comm and synch susceptibility at t=0
        if( params.Laser_Envelope_model ) {
            vecPatches.sumSusceptibility( params, time_dual, timers, 0, simWindow, &smpi );
        }
        
        // Comm and synch charge and current densities
        vecPatches.sumDensities( params, time_dual, timers, 0, simWindow, &smpi );
        
        // rotational cleaning on a single global region
        if( params.initial_rotational_cleaning ) {
            TITLE( "Rotational cleaning" );
            Region region_global( params );
            region_global.build( params, &smpi, vecPatches, openPMD, true );
            region_global.identify_additional_patches( &smpi, vecPatches, params, simWindow );
            region_global.identify_missing_patches( &smpi, vecPatches, params );
            for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  ) {
                DoubleGridsAM::syncFieldsOnRegion( vecPatches, region_global, params, &smpi, imode );
            }
            if( params.is_pxr && smpi.isMaster()) {
                region_global.coupling( params, true );
            }
            for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  ) {
                DoubleGridsAM::syncFieldsOnPatches( region_global, vecPatches, params, &smpi, timers, 0, imode );
            }
            vecPatches.setMagneticFieldsForDiagnostic( params );
            region_global.clean();
            
            if( params.multiple_decomposition ) {
                // Need to upload corrected data on Region
                for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  ) {
                    DoubleGridsAM::syncFieldsOnRegion( vecPatches, region, params, &smpi, imode );
                    // Need to fill all ghost zones, not covered by patches ghost zones
                    SyncVectorPatch::exchangeE( params, region.vecPatch_, imode, &smpi );
                    SyncVectorPatch::exchangeB( params, region.vecPatch_, imode, &smpi );
                }
            }
        }
        
        TITLE( "Open files & initialize diagnostics" );
        vecPatches.initAllDiags( params, &smpi );
        TITLE( "Running diags at time t = 0" );

#ifdef _OMPTASKS
            vecPatches.runAllDiagsTasks( params, &smpi, 0, timers, simWindow );
#else
            vecPatches.runAllDiags( params, &smpi, 0, timers, simWindow );
#endif
    }
    
    TITLE( "Species creation summary" );
    vecPatches.printGlobalNumberOfParticlesPerSpecies( &smpi );
    
    if( params.is_pxr ){
        if( params.multiple_decomposition ) {
            region.coupling( params, false );
        } else {
            vecPatches( 0 )->EMfields->MaxwellAmpereSolver_->coupling( params, vecPatches( 0 )->EMfields );
        }
    }
    
    if( params.is_spectral && params.geometry != "AMcylindrical") {
        vecPatches.saveOldRho( params );
    }
    
    timers.reboot();
    timers.global.reboot();
    
    // ------------------------------------------------------------------------
    // Check expected disk usage
    // ------------------------------------------------------------------------
    TITLE( "Expected disk usage (approximate)" );
    vecPatches.checkExpectedDiskUsage( &smpi, params, checkpoint );
    
    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE( "Keeping or closing the python runtime environment" );
    params.cleanup( &smpi );

    /*tommaso
        // save latestTimeStep (used to test if we are at the latest timestep when running diagnostics at run's end)
        unsigned int latestTimeStep=checkpoint.this_run_start_step;
    */
    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------

    TITLE( "Time-Loop started: number of time-steps n_time = " << params.n_time );
    if( smpi.isMaster() ) {
        params.print_timestep_headers( &smpi );
    }
    
    int count_dlb = 0;
    
    unsigned int itime=checkpoint.this_run_start_step+1;
    while( ( itime <= params.n_time ) && ( !checkpoint.exit_asap ) ) {

        #pragma omp parallel shared (time_dual,smpi,params, vecPatches, region, simWindow, checkpoint, itime)
        {

            // calculate new times
            // -------------------
            #pragma omp single
            {
                time_prim += params.timestep;
                time_dual += params.timestep;
                if( params.keep_python_running_ ) {
                    PyTools::setIteration( itime ); // sets python variable "Main.iteration" for users
                }
            }

            // Patch reconfiguration
            if( params.has_adaptive_vectorization && params.adaptive_vecto_time_selection->theTimeIsNow( itime ) ) {
                vecPatches.reconfiguration( params, timers, itime );
            }

            // apply collisions if requested
            vecPatches.applyCollisions( params, itime, timers );

            // Solve "Relativistic Poisson" problem (including proper centering of fields)
            // for species who stop to be frozen
            // Note: the mean gamma for initialization will be computed for all the species
            // whose fields are initialized at this iteration
            if( params.solve_relativistic_poisson == true ) {
                vecPatches.runRelativisticModule( time_prim, params, &smpi,  timers );
            }

            // Reset global charge and currents densities to zero and computes rho old before moving particles
            if ( params.geometry == "AMcylindrical" && params.is_spectral )
                vecPatches.computeCharge(true);

            // (1) interpolate the fields at the particle position
            // (2) move the particle
            // (3) calculate the currents (charge conserving method)
            vecPatches.dynamics( params, &smpi, simWindow, radiation_tables_,
                                 MultiphotonBreitWheelerTables,
                                 time_dual, timers, itime );

            // if Laser Envelope is used, execute particles and envelope sections of ponderomotive loop
            if( params.Laser_Envelope_model ) {
                vecPatches.runEnvelopeModule( params, &smpi, simWindow, time_dual, timers, itime );
            } // end condition if Laser Envelope Model is used

            // Sum densities
            vecPatches.sumDensities( params, time_dual, timers, itime, simWindow, &smpi );

            // apply currents from antennas
            vecPatches.applyAntennas( time_dual );

        } //End omp parallel region

        // solve Maxwell's equations
        if (!params.multiple_decomposition) {
            if( time_dual > params.time_fields_frozen ) {
                #pragma omp parallel shared (time_dual,smpi,params, vecPatches, region, simWindow, checkpoint, itime)
                {
                    // de-apply prescribed fields if requested
                    if ( vecPatches(0)->EMfields->prescribedFields.size() ) {
                        vecPatches.resetPrescribedFields();
                    }
                    vecPatches.solveMaxwell( params, simWindow, itime, time_dual, timers, &smpi );
                }
                
            }
        }
        else { //if ( params.multiple_decomposition ) {
            if( time_dual > params.time_fields_frozen ) {
                if ( params.geometry != "AMcylindrical" )
                    DoubleGrids::syncCurrentsOnRegion( vecPatches, region, params, &smpi, timers, itime );
                else {
                    for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  )
                        DoubleGridsAM::syncCurrentsOnRegion( vecPatches, region, params, &smpi, timers, itime, imode );
                }
                region.vecPatch_.diag_flag = false;

                //here filter + divergence cleaning
                if ( params.is_spectral && params.geometry == "AMcylindrical") {
                    timers.densitiesCorrection.restart();
                    region.vecPatch_( 0 )->EMfields->MaxwellAmpereSolver_->densities_correction( region.vecPatch_( 0 )->EMfields );
                    timers.densitiesCorrection.update();
                }


                timers.syncDens.restart();
                if( params.geometry != "AMcylindrical" )
                    SyncVectorPatch::sumRhoJ( params, region.vecPatch_, &smpi, timers, itime ); // MPI
                else
                    for( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
                        SyncVectorPatch::sumRhoJ( params, region.vecPatch_, imode, &smpi, timers, itime );
                    }
                timers.syncDens.update( params.printNow( itime ) );


                // apply external time fields if requested
                if( region.vecPatch_(0)->EMfields->prescribedFields.size() ) {
                    region.vecPatch_.applyPrescribedFields( time_prim );
                }

                region.solveMaxwell( params, simWindow, itime, time_dual, timers, &smpi );
                if ( params.geometry != "AMcylindrical" )
                    DoubleGrids::syncFieldsOnPatches( region, vecPatches, params, &smpi, timers, itime );
                else {
                    for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  )
                        DoubleGridsAM::syncFieldsOnPatches( region, vecPatches, params, &smpi, timers, itime, imode );
                }
            }
            if( vecPatches.diag_flag ) {

                if (!params.is_spectral) {
                    if ( params.geometry != "AMcylindrical" )
                        DoubleGrids::syncBOnPatches( region, vecPatches, params, &smpi, timers, itime );
                    else {
                        for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  )
                            DoubleGridsAM::syncBOnPatches( region, vecPatches, params, &smpi, timers, itime, imode );
                    }

                    // Currents and densities not corrected on regions
                    #pragma omp parallel shared (time_dual,smpi,params, vecPatches, region, simWindow, checkpoint, itime)
                    {
                        if( params.geometry != "AMcylindrical" ) {
                            SyncVectorPatch::sumRhoJ( params, vecPatches, &smpi, timers, itime ); // MPI
                        }
                        else {
                            for( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
                                SyncVectorPatch::sumRhoJ( params, vecPatches, imode, &smpi, timers, itime );
                            }
                        }
                    }
                }
                else {
                    // Just need to cp Bm in B for all patches
                    vecPatches.setMagneticFieldsForDiagnostic( params );

                    // Currents and densities could have been corrected on regions
                    if ( params.geometry != "AMcylindrical" ) {
                        DoubleGrids::syncCurrentsOnPatches( region, vecPatches, params, &smpi, timers, itime );
                    }
                    else {
                        for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  )
                            DoubleGridsAM::syncCurrentsOnPatches( region, vecPatches, params, &smpi, timers, itime, imode );
                    }
                }
            }
            bool old = (params.geometry == "AMcylindrical" && params.is_spectral);
            region.vecPatch_.resetRhoJ(old);
        }

        #pragma omp parallel shared (time_dual,smpi,params, vecPatches, region, simWindow, checkpoint, itime)
        {
            // finalize particle exchanges and sort particles
            vecPatches.finalizeAndSortParticles( params, &smpi, simWindow,
                                                 time_dual, timers, itime );

            // Particle merging
            vecPatches.mergeParticles(params, &smpi, time_dual,timers, itime );

            // Particle injection from the boundaries
            vecPatches.injectParticlesFromBoundaries(params, timers, itime );

            // Clean buffers and resize arrays
            vecPatches.cleanParticlesOverhead(params, timers, itime );

            // Finalize field synchronization and exchanges
            vecPatches.finalizeSyncAndBCFields( params, &smpi, simWindow, time_dual, timers, itime );
            
            if( !params.multiple_decomposition ) {
                if( time_dual > params.time_fields_frozen ) {
                    // Standard fields operations (maxwell + comms + boundary conditions) are completed
                    // apply prescribed fields can be considered if requested
                    if( vecPatches(0)->EMfields->prescribedFields.size() ) {
                        #pragma omp single
                        vecPatches.applyPrescribedFields( time_prim );
                        #pragma omp barrier
                    }
                }
            }
            
            // call the various diagnostics
#ifdef _OMPTASKS
            vecPatches.runAllDiagsTasks( params, &smpi, itime, timers, simWindow );
#else
            vecPatches.runAllDiags( params, &smpi, itime, timers, simWindow );
#endif

            timers.movWindow.restart();
            simWindow->shift( vecPatches, &smpi, params, itime, time_dual, region );

            if (itime == simWindow->getAdditionalShiftsIteration() ) {
                int adjust = simWindow->isMoving(time_dual)?0:1;
                for (unsigned int n=0;n < simWindow->getNumberOfAdditionalShifts()-adjust; n++)
                    simWindow->shift( vecPatches, &smpi, params, itime, time_dual, region );
            }
            timers.movWindow.update();
            // ----------------------------------------------------------------------
            // Validate restart  : to do
            // Restart patched moving window : to do
            #pragma omp master
            checkpoint.dump( vecPatches, region, itime, &smpi, simWindow, params );
            #pragma omp barrier
            // ----------------------------------------------------------------------
            
        } //End omp parallel region
        
        if( params.has_load_balancing && params.load_balancing_time_selection->theTimeIsNow( itime ) ) {
            count_dlb++;
            if (params.multiple_decomposition && count_dlb%5 ==0 ) {
                if ( params.geometry != "AMcylindrical" ) {
                    DoubleGrids::syncBOnPatches( region, vecPatches, params, &smpi, timers, itime );
                } else {
                    for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  ) {
                        DoubleGridsAM::syncBOnPatches( region, vecPatches, params, &smpi, timers, itime, imode );
                    }
                }
            }
            
            timers.loadBal.restart();
            #pragma omp single
            vecPatches.loadBalance( params, time_dual, &smpi, simWindow, itime );
            timers.loadBal.update( params.printNow( itime ) );
            
            if( params.multiple_decomposition ) {
                
                if( count_dlb%5 == 0 ) {
                    region.reset_fitting( &smpi, params );
                    region.clean();
                    region.reset_mapping();
                    region.build( params, &smpi, vecPatches, openPMD, false );
                    if( params.is_pxr ) {
                        region.coupling( params, false );
                    }
                    region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
                    region.identify_missing_patches( &smpi, vecPatches, params );
                    
                    if ( params.geometry != "AMcylindrical" ) {
                        DoubleGrids::syncFieldsOnRegion( vecPatches, region, params, &smpi );
                    } else {
                        for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  ) {
                            DoubleGridsAM::syncFieldsOnRegion( vecPatches, region, params, &smpi, imode );
                        }
                    }
                    
                } else {
                    region.reset_mapping();
                    region.identify_additional_patches( &smpi, vecPatches, params, simWindow );
                    region.identify_missing_patches( &smpi, vecPatches, params );
                }
            }
        }
        
        // print message at given time-steps
        // --------------------------------
        if( params.printNow( itime ) ) {
            double npart = vecPatches.getGlobalNumberOfParticles( &smpi );
            params.print_timestep( &smpi, itime, time_dual, timers.global, npart ); //contains a timer.update !!!
            
            #pragma omp master
            timers.consolidate( &smpi );
            #pragma omp barrier
        }
        
        itime++;
        
    }//END of the time loop
    
    smpi.barrier();

    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE( "End time loop, time dual = " << time_dual );
    timers.global.update();

    TITLE( "Time profiling : (print time > 0.001%)" );
    timers.profile( &smpi );

    smpi.barrier();

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
    if (params.multiple_decomposition) {
        region.clean();
    }
    vecPatches.close( &smpi );
    smpi.barrier(); // Don't know why but sync needed by HDF5 Phasespace managment
    delete simWindow;
    PyTools::closePython();
    TITLE( "END" );

    return 0;

}//END MAIN

// ---------------------------------------------------------------------------------------------------------------------
//                                               END MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------


int executeTestMode( VectorPatch &vecPatches,
                     Region &region,
                     SmileiMPI *smpi,
                     SimWindow *simWindow,
                     Params &params,
                     Checkpoint &checkpoint,
                     OpenPMDparams &openPMD,
                     RadiationTables * radiation_tables_ )
{
    int itime = 0;
    int moving_window_movement = 0;

    if( params.restart ) {
        checkpoint.readPatchDistribution( smpi, simWindow );
        itime = checkpoint.this_run_start_step+1;
        moving_window_movement = simWindow->getNmoved();
    }

    PatchesFactory::createVector( vecPatches, params, smpi, openPMD, radiation_tables_, itime, moving_window_movement );

    if( params.restart ) {
        if (params.multiple_decomposition) {
            checkpoint.readRegionDistribution( region );
            region.build( params, smpi, vecPatches, openPMD, false );
        }
        checkpoint.restartAll( vecPatches, region, smpi, simWindow, params, openPMD );
    }

    if( params.print_expected_disk_usage ) {
        TITLE( "Expected disk usage (approximate)" );
        vecPatches.checkExpectedDiskUsage( smpi, params, checkpoint );
    }

    // If test mode enable, code stops here
    TITLE( "Keeping or closing the python runtime environment" );
    params.cleanup( smpi );
    delete simWindow;
    PyTools::closePython();
    TITLE( "END TEST MODE" );

    return 0;
}
