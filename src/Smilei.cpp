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

#include "PatchesFactory.h"
#include "Checkpoint.h"

#include "Solver.h"

#include "SimWindow.h"

#include "Diagnostic.h"

#include "Timer.h"
#include <omp.h>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    cout.setf( ios::fixed,  ios::floatfield ); // floatfield set to fixed
    
    // Define MPI environment :
    SmileiMPI *smpiData= new SmileiMPI(&argc, &argv );
    
    // -------------------------
    // Simulation Initialization
    // ------------------------- 
    
    // Send information on current simulation
    MESSAGE("                   _            _");
    MESSAGE(" ___           _  | |        _  \\ \\    ");
    MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |   Version : " << __VERSION);
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |   Date    : " );//<< __COMMITDATE);
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |   " );//<< (string(__CONFIG).size()? "Config  : ":"") << __CONFIG);
    MESSAGE("                                /_/    ");
    
    TITLE("Input data info");
    
    // Read simulation parameters
    Params params(smpiData,vector<string>(argv + 1, argv + argc));
    smpiData->init(params);
    smpiData->barrier();
    if ( smpiData->isMaster() ) params.print();
    smpiData->barrier();
    
    // Initialize timers
    vector<Timer> timer;
    initialize_timers(timer, smpiData);
    
    // Print in stdout MPI, OpenMP, patchs parameters
    print_parallelism_params(params, smpiData);
    
    TITLE("Restart environments");
    Checkpoint checkpoint(params, smpiData);
    
    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // Update in "if restart" if necessary
    // ------------------------------------------------------------------------
    
    unsigned int stepStart=0, stepStop=params.n_time;
    
    // time at integer time-steps (primal grid)
    double time_prim = stepStart * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (stepStart +0.5) * params.timestep;
    // Do we initially do diags or not ?
    int diag_flag = 1;
    
    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------
    // --------------------
    // Define Moving Window
    // --------------------
    TITLE("Initializing moving window");
    int start_moving(0);
    SimWindow* simWindow = new SimWindow(params);
    params.hasWindow = simWindow->isActive();
    
    // ---------------------------------------------------
    // Initialize patches (including particles and fields)
    // ---------------------------------------------------
    TITLE("Initializing particles & fields");
    VectorPatch vecPatches = PatchesFactory::createVector(params, smpiData);
    
    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        checkpoint.restartAll( vecPatches, stepStart, smpiData, simWindow, params);
        
        // time at integer time-steps (primal grid)
        time_prim = checkpoint.this_run_start_step * params.timestep;
        // time at half-integer time-steps (dual grid)
        time_dual = (checkpoint.this_run_start_step +0.5) * params.timestep;
        
        double restart_time_dual = (checkpoint.this_run_start_step +0.5) * params.timestep;
        time_dual = restart_time_dual;
        // A revoir !
        if ( simWindow->isMoving(restart_time_dual) ) {
            simWindow->operate(vecPatches, smpiData, params);
        }
        //smpiData->recompute_patch_count( params, vecPatches, restart_time_dual );
        
    } else {
        
        // Initialize the electromagnetic fields
        // -----------------------------------
        vecPatches.dynamics(params, smpiData, simWindow, &diag_flag, time_dual, timer);
        timer[1].reboot();
        timer[8].reboot();
        
        vecPatches.sumDensities( &diag_flag, timer );
        timer[4].reboot();
        timer[9].reboot();
        
        #pragma omp single
        {
        if( vecPatches.hasAntennas )
            TITLE("Applying antennas at time t = " << 0.5 * params.timestep);
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
                vecPatches(ipatch)->EMfields->applyAntennas(smpiData, 0.5 * params.timestep); // smpi useless
        }
        
        // Init electric field (Ex/1D, + Ey/2D)
        if (!vecPatches.isRhoNull(smpiData)) {
            TITLE("Solving Poisson at time t = 0");
            Timer ptimer;
            ptimer.init(smpiData, "global");
            ptimer.restart();
            vecPatches.solvePoisson( params, smpiData );
            ptimer.update();
            MESSAGE("Time in Poisson : " << ptimer.getTime() );
        }
        
        TITLE("Applying external fields at time t = 0");
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
            vecPatches(ipatch)->EMfields->applyExternalFields( vecPatches(ipatch) ); // Must be patch
        
        TITLE("Running diags at time t = 0");
        vecPatches.runAllDiags(params, smpiData, &diag_flag, 0, timer);
        timer[3].reboot();
        timer[6].reboot();
    
    }
    timer[0].reboot();
    
    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    params.cleanup(smpiData);
    
    // ------------------------------------------------------------------------
    // Check memory consumption
    // ------------------------------------------------------------------------
    check_memory_consumption( vecPatches, smpiData );
    
    double old_print_time(0.), this_print_time;
     
    // save latestTimeStep (used to test if we are at the latest timestep when running diagnostics at run's end)
    unsigned int latestTimeStep=checkpoint.this_run_start_step;
    bool exit(false);
    
    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    
    TITLE("Time-Loop started: number of time-steps n_time = " << params.n_time);
    for (unsigned int itime=checkpoint.this_run_start_step+1 ; itime <= stepStop ; itime++) {
        
        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        if ( vecPatches.fieldTimeIsNow(itime) ) diag_flag = 1;
        
        // send message at given time-steps
        // --------------------------------
        timer[0].update();
        
        if ( vecPatches.printScalars( itime ) &&  ( smpiData->isMaster() ) ) {
            old_print_time = this_print_time;
            this_print_time=timer[0].getTime();
            ostringstream my_msg;
            my_msg << setw(log10(params.n_time)+1) << itime <<
            "/"     << setw(log10(params.n_time)+1) << params.n_time <<
            " t="          << scientific << setprecision(3)   << time_dual <<
            " sec "    << scientific << setprecision(1)   << this_print_time <<
            " ("    << scientific << setprecision(4)   << this_print_time - old_print_time << ")" <<
            "  Utot= "   << scientific << setprecision(4)<< vecPatches.getScalar("Utot") <<
            "  Uelm= "   << scientific << setprecision(4)<< vecPatches.getScalar("Uelm") <<
            "  Ukin= "   << scientific << setprecision(4)<< vecPatches.getScalar("Ukin") <<
            "  Ubal(%)= "<< scientific << fixed << setprecision(2) << 100.0*vecPatches.getScalar("Ubal_norm");
            
            if ( simWindow->isActive() ) {
                double Uinj_mvw = vecPatches.getScalar("Uelm_inj_mvw") + vecPatches.getScalar("Ukin_inj_mvw");
                double Uout_mvw = vecPatches.getScalar("Uelm_out_mvw") + vecPatches.getScalar("Ukin_out_mvw");
                my_msg << "  Uinj_mvw = " << scientific << setprecision(4) << Uinj_mvw <<
                "  Uout_mvw = " << scientific << setprecision(4) << Uout_mvw;

            }//simWindow

            MESSAGE(my_msg.str());
        }//itime
        
        
        // put density and currents to 0 + save former density
        // ---------------------------------------------------
        
        timer[10].restart();
        // apply collisions if requested
        // -----------------------------
        if (Collisions::debye_length_required)
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
                Collisions::calculate_debye_length(params,vecPatches(ipatch)->vecSpecies);
        for (unsigned int icoll=0 ; icoll<vecPatches(0)->vecCollisions.size(); icoll++)
            vecPatches(0)->vecCollisions[icoll]->createTimestep(itime);
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            for (unsigned int icoll=0 ; icoll<vecPatches(ipatch)->vecCollisions.size(); icoll++)
                vecPatches(ipatch)->vecCollisions[icoll]->collide(params,vecPatches(ipatch),itime);
        timer[10].update();
        
        /*******************************************/
        /********** Move particles *****************/
        /*******************************************/
        #pragma omp parallel shared (time_dual,smpiData,params, vecPatches, simWindow)
        {
            // apply the PIC method
            // --------------------
            // for all particles of all species (see dynamic in Species.cpp)
            // (1) interpolate the fields at the particle position
            // (2) move the particle
            // (3) calculate the currents (charge conserving method)
            vecPatches.dynamics(params, smpiData, simWindow, &diag_flag, time_dual, timer);
            
            /*******************************************/
            /*********** Sum densities *****************/
            /*******************************************/
            vecPatches.sumDensities( &diag_flag, timer );
            
            // apply currents from antennas
            #pragma omp single
            {
            if( vecPatches.hasAntennas )
                for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
                    vecPatches(ipatch)->EMfields->applyAntennas(smpiData, time_dual);
            }
            
            /*******************************************/
            /*********** Maxwell solver ****************/
            /*******************************************/
            
            // solve Maxwell's equations
            if( time_dual > params.time_fields_frozen )
                vecPatches.solveMaxwell( params, simWindow, itime, time_dual, timer );
            
            // call the various diagnostics
            // ----------------------------
            #pragma omp master
            vecPatches.runAllDiags(params, smpiData, &diag_flag, itime, timer);
            #pragma omp barrier
            
            
            // ----------------------------------------------------------------------
            // Validate restart  : to do
            // Restart patched moving window : to do
            // Break in an OpenMP region
            #pragma omp master
            exit = checkpoint.dump(vecPatches, itime, smpiData, simWindow, params);
            #pragma omp barrier
            // ----------------------------------------------------------------------        
            
        } //End omp parallel region

        if (exit) break;
        
        timer[5].restart();
        if ( simWindow->isMoving(time_dual) ) {
            start_moving++;
            if ((start_moving==1) && (smpiData->isMaster()) ) {
                MESSAGE(">>> Window starts moving");
            }
            simWindow->operate(vecPatches, smpiData, params);
        }
        timer[5].update();
        
        
        
        if ((itime%params.balancing_every == 0)&&(smpiData->getSize()!=1)) {
            timer[7].restart();
            vecPatches.load_balance( params, time_dual, smpiData, simWindow );
            timer[7].update();
        }
        
        latestTimeStep = itime;
        
    }//END of the time loop
    
    smpiData->barrier();
    
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE("End time loop, time dual = " << time_dual);
    
    //double timElapsed=smpiData->time_seconds();
    //if ( smpiData->isMaster() ) MESSAGE(0, "Time in time loop : " << timElapsed );
    timer[0].update();
    TITLE("Time profiling :");
    double coverage(0.);
    for (unsigned int i=1 ; i<timer.size() ; i++) coverage += timer[i].getTime();
    MESSAGE("Time in time loop :\t" << timer[0].getTime() << "\t"<<coverage/timer[0].getTime()*100.<< "% coverage" );
    if ( smpiData->isMaster() )
        for (unsigned int i=1 ; i<timer.size() ; i++) timer[i].print(timer[0].getTime());
    
    //WARNING( "Diabled vecPatches.Diagnostics->printTimers(vecPatches(0), timer[3].getTime());" );
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    if (latestTimeStep==params.n_time)
        vecPatches.runAllDiags(params, smpiData, &diag_flag, params.n_time, timer);
    
    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    vecPatches.close( smpiData );
    
    MPI_Barrier(MPI_COMM_WORLD); // Don't know why but sync needed by HDF5 Phasespace managment
    
    delete simWindow;
    
    PyTools::closePython();
    
    TITLE("END");
    delete smpiData;
    
    return 0;
    
}//END MAIN

// ---------------------------------------------------------------------------------------------------------------------
//                                               END MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------


void print_parallelism_params(Params& params, SmileiMPI* smpi)
{
    TITLE("MPI");
    MESSAGE(1,"Number of MPI process : " << smpi->getSize() );
    MESSAGE(1,"Number of patches : " );
    for (unsigned int iDim=0 ; iDim<params.nDim_field ; iDim++) 
        MESSAGE(2, "dimension " << iDim << " - number_of_patches : " << params.number_of_patches[iDim] );

    MESSAGE(1, "Patch size :");
    for (unsigned int iDim=0 ; iDim<params.nDim_field ; iDim++) 
        MESSAGE(2, "dimension " << iDim << " - n_space : " << params.n_space[iDim] << " cells.");        

    MESSAGE(1, "Dynamic load balancing frequency: every " << params.balancing_every << " iterations." );

    // setup OpenMP
    TITLE("OpenMP");
#ifdef _OPENMP
//    int nthds(0);
//#pragma omp parallel shared(nthds)
//    {
//        nthds = omp_get_num_threads();
//    }
    if (smpi->isMaster())
        MESSAGE(1,"Number of thread per MPI process : " << omp_get_max_threads() );
#else
    if (smpi->isMaster()) MESSAGE("Disabled");
#endif

} // End print_parallelism_params


void check_memory_consumption(VectorPatch& vecPatches, SmileiMPI* smpi)
{
    TITLE("Memory consumption");
    
    int particlesMem(0);
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
        for (unsigned int ispec=0 ; ispec<vecPatches(ipatch)->vecSpecies.size(); ispec++)
            particlesMem += vecPatches(ipatch)->vecSpecies[ispec]->getMemFootPrint();
    MESSAGE( 1, "(Master) Species part = " << (int)( (double)particlesMem / 1024./1024.) << " Mo" );

    double dParticlesMem = (double)particlesMem / 1024./1024./1024.;
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dParticlesMem, &dParticlesMem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( 1, setprecision(3) << "Global Species part = " << dParticlesMem << " Go" );

    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&particlesMem, &particlesMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( 1, "Max Species part = " << (int)( (double)particlesMem / 1024./1024.) << " Mb" );
    
    // fieldsMem contains field per species
    int fieldsMem(0);
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
        fieldsMem = vecPatches(ipatch)->EMfields->getMemFootPrint();
    MESSAGE( 1, "(Master) Fields part = " << (int)( (double)fieldsMem / 1024./1024.) << " Mo" );

    double dFieldsMem = (double)fieldsMem / 1024./1024./1024.;
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dFieldsMem, &dFieldsMem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( 1, setprecision(3) << "Global Fields part = " << dFieldsMem << " Go" );
    
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&fieldsMem, &fieldsMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( 1, "Max Fields part = " << (int)( (double)fieldsMem / 1024./1024.) << " Mb" );

    // Read value in /proc/pid/status
    //Tools::printMemFootPrint( "End Initialization" );

} // End check_memory_consumption


void initialize_timers(vector<Timer>& timer, SmileiMPI* smpi)
{
    // GC IDRIS : "Timer timer[ntimer];" to "Timer timer[8];"
    int ntimer(13);
    timer.resize(ntimer);
    // The entire time loop
    timer[0].init(smpi, "Global");
    // Call dynamics + restartRhoJ(s)
    timer[1].init(smpi, "Particles");
    // Maxwell
    timer[2].init(smpi, "Maxwell");
    // Diags.runAllDiags + MPI & Patch sync
    timer[3].init(smpi, "Diagnostics");
    // Local sum of rho, Jxyz
    timer[4].init(smpi, "Densities");
    // Moving Window
    timer[5].init(smpi, "Mov window");
    // Dump fields (including average)
    timer[6].init(smpi, "Diag fields");
    // Load balancing
    timer[7].init(smpi, "Load balacing");
    // Call exchangeParticles (MPI & Patch sync)
    timer[8].init(smpi, "Sync Particles");
    // Call sumRhoJ(s), exchangeB (MPI & Patch sync)
    timer[9].init(smpi, "Sync Fields");
    // Call to Collisions methods
    timer[10].init(smpi, "Collisions");
    // If necessary the following timers can be reintroduced
    //timer[11].init(smpi, "Fields");
    //timer[12].init(smpi, "AvgFields");



} // End initialize_timers
