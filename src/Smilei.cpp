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
#include "PartWall.h"
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
    
    // Send information on current simulation
    MESSAGE("                   _            _");
    MESSAGE(" ___           _  | |        _  \\ \\    ");
    MESSAGE("/ __|  _ __   (_) | |  ___  (_)  | |   Version : " << __VERSION);
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |   Date    : " << __COMMITDATE);
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |   " << (string(__CONFIG).size()? "Config  : ":"") << __CONFIG);
    MESSAGE("                                /_/    ");
    
    TITLE("Input data info");
    
    // Check for namelists (input files)
    vector<string> namelists(argv + 1, argv + argc);
    if (namelists.size()==0) ERROR("No namelists given!");
    
    // Read simulation & diagnostics parameters
    Params params(smpiData,namelists);
    smpiData->init(params);
    smpiData->barrier();
    if ( smpiData->isMaster() ) params.print();
    smpiData->barrier();
    
    
    // Geometry known, MPI environment specified
    TITLE("General MPI environement");
    SmileiMPI* smpi = SmileiMPIFactory::create(params, smpiData);
    
    // ---------------------------
    // Initialize Species
    // ---------------------------
    TITLE("Initializing species");
    
    // Initialize the vecSpecies vector containing all information of the different Species
    // ------------------------------------------------------------------------------------
    vector<Species*> vecSpecies = SpeciesFactory::createVector(params, smpi);

    // Create diagnostics
    TITLE("Creating Diagnostics");
    Diagnostic diags(params, vecSpecies, smpi);
    
    //Create mpi i/o environment
    TITLE("MPI input output environment");
    SmileiIO*  sio  = SmileiIOFactory::create(params, diags, smpi);
    
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
    // Initialize Fields
    // ---------------------------
    TITLE("Initializing fields & moving-window");

    // Initialize the collisions (vector of collisions)
    // ------------------------------------------------------------------------------------
    vector<Collisions*> vecCollisions = Collisions::create(params, vecSpecies, smpi);
    
    // Initialize the particle walls
    vector<PartWall*> vecPartWall = PartWall::create(params, smpi);
    
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
    ElectroMagn* EMfields = ElectroMagnFactory::create(params, vecSpecies, smpi);
    
    // interpolation operator (virtual)
    Interpolator* Interp = InterpolatorFactory::create(params, smpi);
    
    // projection operator (virtual)
    Projector* Proj = ProjectorFactory::create(params, smpi);
    smpi->barrier();
    
    unsigned int stepStop=params.n_time;
    
    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        DEBUG(vecSpecies.size());
        sio->restartAll( EMfields, vecSpecies, smpi, simWindow, params, diags);
                
        double restart_time_dual = (sio->this_run_start_step +0.5) * params.timestep;
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
        for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
            smpi->sumRhoJs(EMfields, ispec, true);  // only if !isTest
        }
        
        TITLE("Applying antennas at time t = " << 0.5 * params.timestep);
        EMfields->applyAntennas(smpi, 0.5 * params.timestep);
        
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
        diags.runAllDiags(0, EMfields, vecSpecies, Interp, smpi);
        // temporary EM fields dump in Fields.h5
        sio->writeAllFieldsSingleFileTime(EMfields->allFields, 0, 0 );
        // temporary EM fields dump in Fields_avg.h5
        if (diags.ntime_step_avg!=0)
            sio->writeAllFieldsSingleFileTime(EMfields->allFields_avg, 0, 1 );
        
    }

    // ------------------------------------------------------------------------
    // Check memory consumption
    // ------------------------------------------------------------------------
    TITLE("Memory consumption");
    
    int particlesMem(0);
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++)
        particlesMem += vecSpecies[ispec]->getMemFootPrint();
    MESSAGE( "(Master) Species part = " << (int)( (double)particlesMem / 1024./1024.) << " Mb" );
    
    double dParticlesMem = (double)particlesMem / 1024./1024./1024.;
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dParticlesMem, &dParticlesMem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( setprecision(3) << "Global Species part = " << dParticlesMem << " Gb" );
    
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&particlesMem, &particlesMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( "Max Species part = " << (int)( (double)particlesMem / 1024./1024.) << " Mb" );
    
    // fieldsMem contains field per species
    int fieldsMem = EMfields->getMemFootPrint();
    MESSAGE( "(Master) Fields part = " << (int)( (double)fieldsMem / 1024./1024.) << " Mb" );
    
    double dFieldsMem = (double)fieldsMem / 1024./1024./1024.;
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dFieldsMem, &dFieldsMem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( setprecision(3) << "Global Fields part = " << dFieldsMem << " Gb" );
    
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&fieldsMem, &fieldsMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( "Max Fields part = " << (int)( (double)fieldsMem / 1024./1024.) << " Mb" );

    // Read value in /proc/pid/status
    //Tools::printMemFootPrint( "End Initialization" );
    
    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    params.cleanup(smpi);
    
    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // ------------------------------------------------------------------------
    
    // time at integer time-steps (primal grid)
    double time_prim = sio->this_run_start_step * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (sio->this_run_start_step +0.5) * params.timestep;
    
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
    
    // save latestTimeStep (used to test if we are at the latest timestep when running diagnostics at run's end)
    unsigned int latestTimeStep=sio->this_run_start_step;
    
    TITLE("Time-Loop started: number of time-steps n_time = " << params.n_time);
    double old_print_time=timer[0].getTime();
    for (unsigned int itime=sio->this_run_start_step+1 ; itime <= stepStop ; itime++) {
        
        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        // send message at given time-steps
        // --------------------------------
        timer[0].update();
        
        if ( (itime % diags.print_every == 0) &&  ( smpi->isMaster() ) ) {
            
            double this_print_time=timer[0].getTime();

            ostringstream my_msg;
            my_msg << setw(log10(params.n_time)+1) << itime <<
            "/"     << setw(log10(params.n_time)+1) << params.n_time <<
            "  t "          << scientific << setprecision(3)   << time_dual <<
            "  sec "    << scientific << setprecision(1)   << this_print_time <<
            "("    << scientific << setprecision(1)   << this_print_time - old_print_time << ")" <<
            "  Utot "   << scientific << setprecision(4)<< diags.getScalar("Utot") <<
            "  Uelm "   << scientific << setprecision(4)<< diags.getScalar("Uelm") <<
            "  Ukin "   << scientific << setprecision(4)<< diags.getScalar("Ukin") <<
            "  Ubal "<< scientific << fixed << setprecision(2) << 100.0*diags.getScalar("Ubal_norm") << "%";
            
            old_print_time=this_print_time;
            
            if (simWindow) {
                double Uinj_mvw = diags.getScalar("Uelm_inj_mvw") + diags.getScalar("Ukin_inj_mvw");
                double Uout_mvw = diags.getScalar("Uelm_out_mvw") + diags.getScalar("Ukin_out_mvw");
                my_msg << " Uinj_mvw " << scientific << setprecision(4) << Uinj_mvw <<
                " Uout_mvw " << scientific << setprecision(4) << Uout_mvw;

            }//simWindow

            MESSAGE(my_msg.str());
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
#pragma omp parallel shared (EMfields,time_dual,smpi,params)
        {
            int tid(0);
#ifdef _OMP
            tid = omp_get_thread_num();
#endif
            for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
                if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ){
                    EMfields->restartRhoJs(ispec, time_dual > vecSpecies[ispec]->time_frozen); // if (!isTest)
                    vecSpecies[ispec]->dynamics(time_dual, EMfields, Interp, Proj, smpi, params, simWindow, vecPartWall);
                }
            }
            for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
#pragma omp barrier
                if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ){
                    // Loop on dims to manage exchange in corners
                    for ( int iDim = 0 ; iDim<(int)params.nDim_particle ; iDim++ )
                        smpi->exchangeParticles(vecSpecies[ispec], params, tid, iDim);
#pragma omp barrier
                        vecSpecies[ispec]->sort_part(); // Should we sort test particles ?? (JD)
                }
            }
        }
        timer[1].update();
        
        //!\todo To simplify : sum global and per species densities
        timer[4].restart();
        smpi->sumRhoJ( EMfields );
        for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
            if ( vecSpecies[ispec]->isProj(time_dual, simWindow) ) smpi->sumRhoJs(EMfields, ispec, time_dual > vecSpecies[ispec]->time_frozen);
        }
        EMfields->computeTotalRhoJ();
        
        // apply currents from antennas
        EMfields->applyAntennas(smpi, time_dual);
        
        timer[4].update();
        
        // solve Maxwell's equations
        if( time_dual > params.time_fields_frozen ) {
            timer[2].restart();
            EMfields->solveMaxwell(itime, time_dual, smpi, params, simWindow);
            timer[2].update();
        }
        
        // incrementing averaged electromagnetic fields
        if (diags.ntime_step_avg) EMfields->incrementAvgFields(itime, diags.ntime_step_avg);
        
        // call the various diagnostics
        // ----------------------------
        
        // run all diagnostics
        timer[3].restart();
        diags.runAllDiags(itime, EMfields, vecSpecies, Interp, smpi);
        timer[3].update();
        
        timer[6].restart();
        // temporary EM fields dump in Fields.h5
        if  ((diags.fieldDump_every != 0) && (itime % diags.fieldDump_every == 0)) {
            sio->writeAllFieldsSingleFileTime( EMfields->allFields, itime, 0 );
        }
        timer[6].update();
        
        timer[7].restart();
        // temporary EM fields dump in Fields_avg.h5
        if  ((diags.ntime_step_avg!=0) &&
             (diags.avgfieldDump_every != 0) && 
             (itime % diags.avgfieldDump_every == 0)) {
            sio->writeAllFieldsSingleFileTime(EMfields->allFields_avg, itime, 1 );
        }
        timer[7].update();
                
        if (sio->dump(EMfields, itime, vecSpecies, smpi, simWindow, params, diags)) break;
        
        timer[5].restart();
        if ( simWindow && simWindow->isMoving(time_dual) ) {
            start_moving++;
            if ((start_moving==1) && (smpi->isMaster()) ) {
                MESSAGE(">>> Window starts moving");
            }
            simWindow->operate(vecSpecies, EMfields, Interp, Proj, smpi, params);
        }
        timer[5].update();
        
        latestTimeStep = itime;
        
    }//END of the time loop
    
    smpi->barrier();
    
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE("End time loop, time dual = " << time_dual);
    
    //double timElapsed=smpiData->time_seconds();
    //if ( smpi->isMaster() ) MESSAGE("Time in time loop : " << timElapsed );
    timer[0].update();
    TITLE("Time profiling :");
    
    double coverage(0.);
    for (unsigned int i=1 ; i<timer.size() ; i++) coverage += timer[i].getTime();
    MESSAGE("Time in time loop :\t" << timer[0].getTime() << "\t"<<coverage/timer[0].getTime()*100.<< "% coverage" );
    if ( smpi->isMaster() )
        for (unsigned int i=1 ; i<timer.size() ; i++) timer[i].print(timer[0].getTime());
    
    diags.printTimers(smpi,timer[3].getTime());
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    if (latestTimeStep==params.n_time) {
    // temporary EM fields dump in Fields.h5
    if  ( (diags.fieldDump_every != 0) && (params.n_time % diags.fieldDump_every != 0) )
        sio->writeAllFieldsSingleFileTime(EMfields->allFields, params.n_time, 0 );
    // temporary time-averaged EM fields dump in Fields_avg.h5
    if  (diags.ntime_step_avg!=0)
        if  ( (diags.avgfieldDump_every != 0) && (params.n_time % diags.avgfieldDump_every != 0) )
            sio->writeAllFieldsSingleFileTime(EMfields->allFields_avg, params.n_time, 1 );
}
    
    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    delete Proj;
    delete Interp;
    delete EMfields;
    for(unsigned int i=0; i<vecCollisions.size(); i++) delete vecCollisions[i];
    vecCollisions.clear();
    
    diags.closeAll(smpi);
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    
    if (params.nspace_win_x)
        delete simWindow;
    
    PyTools::closePython();

    TITLE("END");
    delete sio;
    delete smpi;
    delete smpiData;
    
    return 0;
    
}//END MAIN



