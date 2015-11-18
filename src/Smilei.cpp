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
    MESSAGE("\\__ \\ | '  \\   _  | | / -_)  _   | |   Date    : " << __COMMITDATE);
    MESSAGE("|___/ |_|_|_| |_| |_| \\___| |_|  | |   " << (string(__CONFIG).size()? "Config  : ":"") << __CONFIG);
    MESSAGE("                                /_/    ");
    
    TITLE("Input data info");
    
    // Read simulation & diagnostics parameters
    Params params(smpiData,vector<string>(argv + 1, argv + argc));
    smpiData->init(params);
    smpiData->barrier();
    if ( smpiData->isMaster() ) params.print();
    smpiData->barrier();
    
    
    // setup OpenMP
    TITLE("OpenMP");
#ifdef _OPENMP
    int nthds(0);
#pragma omp parallel shared(nthds)
    {
        nthds = omp_get_num_threads();
    }
    if (smpiData->isMaster())
	MESSAGE(1,"Number of thread per MPI process : " << omp_get_max_threads() );
#else
    if (smpiData->isMaster()) MESSAGE("Disabled");
#endif

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
    // ----------------------------------------------------------------------------
    // Define Moving Window & restart
    // ----------------------------------------------------------------------------
    SimWindow* simWindow = NULL;
    int start_moving(0);
    if (params.nspace_win_x)
        simWindow = new SimWindow(params);
    
    // ---------------------------
    // Initialize Species & Fields
    // ---------------------------
    TITLE("Initializing particles, fields & moving-window");
    
    VectorPatch vecPatches = PatchesFactory::createVector(params, smpiData);
    vecPatches.initProbesDiags(params, 0);
    vecPatches.initDumpFields(params, 0);


    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        DEBUG(vecSpecies.size());
        checkpoint.restartAll( vecPatches, stepStart, smpiData, simWindow, params);

	// time at integer time-steps (primal grid)
	time_prim = checkpoint.this_run_start_step * params.timestep;
	// time at half-integer time-steps (dual grid)
	time_dual = (checkpoint.this_run_start_step +0.5) * params.timestep;
    
        double restart_time_dual = (checkpoint.this_run_start_step +0.5) * params.timestep;
	time_dual = restart_time_dual;
        // A revoir !
	if ( simWindow ) {
	    simWindow->setOperators(vecPatches);
	    if ( simWindow->isMoving(restart_time_dual) ) {
	        simWindow->operate(vecPatches, smpiData, params);
	    }
	}
        //smpiData->recompute_patch_count( params, vecPatches, restart_time_dual );
	

    } else {
	
        // Initialize the electromagnetic fields
        // -----------------------------------
        // Init rho and J by projecting all particles of subdomain
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
	    vecPatches(ipatch)->EMfields->restartRhoJs();
	    vecPatches(ipatch)->dynamics(time_dual, params, simWindow, diag_flag); //include test
	}
	for (unsigned int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size(); ispec++) {
	    if ( vecPatches(0)->vecSpecies[ispec]->isProj(time_dual, simWindow) )
		vecPatches.exchangeParticles(ispec, params, smpiData ); // Included sort_part
	}
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
	    vecPatches(ipatch)->EMfields->computeTotalRhoJ(); // Per species in global, Attention if output -> Sync / per species fields
	vecPatches.sumRhoJ( diag_flag ); // MPI
	for (unsigned int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size(); ispec++)
	    vecPatches.sumRhoJs( ispec ); // MPI
        diag_flag = 0;

	TITLE("Applying antennas at time t = " << 0.5 * params.timestep);
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
	    vecPatches(ipatch)->EMfields->applyAntennas(smpiData, 0.5 * params.timestep); // smpi useless

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
	
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
	    vecPatches(ipatch)->Diags->runAllDiags(0, vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, vecPatches(ipatch)->Interp);
	vecPatches.computeGlobalDiags(0);
	smpiData->computeGlobalDiags( vecPatches(0)->Diags, 0);
 
	// temporary EM fields dump in Fields.h5
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)  {
	    if (ipatch==0) vecPatches(ipatch)->sio->createTimeStepInSingleFileTime( 0, vecPatches.Diags );
	    vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields->allFields, 0, 0 );
	}
        // temporary EM fields dump in Fields_avg.h5
        if (vecPatches.Diags->ntime_step_avg!=0)
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		//if (ipatch==0) vecPatches(ipatch)->sio->createTimeStepInSingleFileTime( vecPatches.Diags );
		vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields->allFields_avg, 0, 1 );
	    }
	//for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	//vecPatches(ipatch)->EMfields->restartRhoJs();
	diag_flag = 0 ;

#ifdef _TRACKPARTICLES
	// Test particles need initialization now (after vecSpecies has been created)
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    for (unsigned int i=0 ; i<vecPatches(ipatch)->Diags->vecDiagnosticTrackParticles.size(); i++)
		vecPatches(ipatch)->Diags->vecDiagnosticTrackParticles[i]->init(vecPatches(ipatch)->vecSpecies, smpiData);
#endif


    }
    
    // ------------------------------------------------------------------------
    // Check memory consumption
    // ------------------------------------------------------------------------
    TITLE("Memory consumption");
    
    int particlesMem(0);
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	for (unsigned int ispec=0 ; ispec<vecPatches(ipatch)->vecSpecies.size(); ispec++)
	    particlesMem += vecPatches(ipatch)->vecSpecies[ispec]->getMemFootPrint();
    MESSAGE( "(Master) Species part = " << (int)( (double)particlesMem / 1024./1024.) << " Mo" );

    double dParticlesMem = (double)particlesMem / 1024./1024./1024.;
    MPI_Reduce( smpiData->isMaster()?MPI_IN_PLACE:&dParticlesMem, &dParticlesMem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( setprecision(3) << "Global Species part = " << dParticlesMem << " Go" );

    MPI_Reduce( smpiData->isMaster()?MPI_IN_PLACE:&particlesMem, &particlesMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( "Max Species part = " << (int)( (double)particlesMem / 1024./1024.) << " Mb" );
    
    // fieldsMem contains field per species
    int fieldsMem(0);
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	fieldsMem = vecPatches(ipatch)->EMfields->getMemFootPrint();
    MESSAGE( "(Master) Fields part = " << (int)( (double)fieldsMem / 1024./1024.) << " Mo" );

    double dFieldsMem = (double)fieldsMem / 1024./1024./1024.;
    MPI_Reduce( smpiData->isMaster()?MPI_IN_PLACE:&dFieldsMem, &dFieldsMem, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( setprecision(3) << "Global Fields part = " << dFieldsMem << " Go" );
    
    MPI_Reduce( smpiData->isMaster()?MPI_IN_PLACE:&fieldsMem, &fieldsMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( "Max Fields part = " << (int)( (double)fieldsMem / 1024./1024.) << " Mb" );

    // Read value in /proc/pid/status
    //Tools::printMemFootPrint( "End Initialization" );
    
    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    params.cleanup(smpiData);


int partperMPI;
int balancing_freq = 150;
int npatchmoy=0, npartmoy=0;
	
    params.cleanup(smpiData);
    
    // Count timer
    int ntimer(13);
    vector<Timer> timer(ntimer);
    timer[0].init(smpiData, "Global");
    timer[1].init(smpiData, "Particles");
    timer[2].init(smpiData, "Maxwell");
    timer[3].init(smpiData, "Diagnostics");
    timer[4].init(smpiData, "Densities");
    timer[5].init(smpiData, "Mov window");
    timer[6].init(smpiData, "Diag fields");
    timer[7].init(smpiData, "Load balacing");
    timer[8].init(smpiData, "Sync Particles");
    timer[9].init(smpiData, "Sync Fields");
    timer[10].init(smpiData, "Fields");
    timer[11].init(smpiData, "AvgFields");
    timer[12].init(smpiData, "Collisions");



    // Action to send to other MPI procs when an action is required
    int mpisize,itime2dump(-1),todump(0); 
    double starttime = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    MPI_Request action_srequests[mpisize];
    MPI_Request action_rrequests;
    MPI_Status action_status[2];

    
    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    
    // save latestTimeStep (used to test if we are at the latest timestep when running diagnostics at run's end)
    unsigned int latestTimeStep=checkpoint.this_run_start_step;
    
    TITLE("Time-Loop started: number of time-steps n_time = " << params.n_time);
    double old_print_time=timer[0].getTime();
    for (unsigned int itime=checkpoint.this_run_start_step+1 ; itime <= stepStop ; itime++) {
        
        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        if  ((vecPatches.Diags->fieldDump_every != 0) && (itime % vecPatches.Diags->fieldDump_every == 0)) diag_flag = 1;

        // send message at given time-steps
        // --------------------------------
        timer[0].update();
        
        if ( (itime % vecPatches.Diags->print_every == 0) &&  ( smpiData->isMaster() ) ) {
            double this_print_time=timer[0].getTime();
            ostringstream my_msg;
            my_msg << setw(log10(params.n_time)+1) << itime <<
            "/"     << setw(log10(params.n_time)+1) << params.n_time <<
	    "t = "          << scientific << setprecision(3)   << time_dual <<
            "  sec "    << scientific << setprecision(1)   << this_print_time <<
            "("    << scientific << setprecision(1)   << this_print_time - old_print_time << ")" <<
            "   Utot = "   << scientific << setprecision(4)<< vecPatches.Diags->getScalar("Utot") <<
            "   Uelm = "   << scientific << setprecision(4)<< vecPatches.Diags->getScalar("Uelm") <<
            "   Ukin = "   << scientific << setprecision(4)<< vecPatches.Diags->getScalar("Ukin") <<
            "   Ubal(%) = "<< scientific << fixed << setprecision(2) << 100.0*vecPatches.Diags->getScalar("Ubal_norm");
            
            if (simWindow) {
                double Uinj_mvw = vecPatches.Diags->getScalar("Uelm_inj_mvw") + vecPatches.Diags->getScalar("Ukin_inj_mvw");
                double Uout_mvw = vecPatches.Diags->getScalar("Uelm_out_mvw") + vecPatches.Diags->getScalar("Ukin_out_mvw");
                my_msg << "   Uinj_mvw = " << scientific << setprecision(4) << Uinj_mvw <<
                "   Uout_mvw = " << scientific << setprecision(4) << Uout_mvw;

            }//simWindow

            MESSAGE(my_msg.str());
        }//itime
        
        
        // put density and currents to 0 + save former density
        // ---------------------------------------------------
        
        timer[8].restart();
        // apply collisions if requested
        // -----------------------------
        if (Collisions::debye_length_required)
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
		Collisions::calculate_debye_length(params,vecPatches(ipatch)->vecSpecies);
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    for (unsigned int icoll=0 ; icoll<vecPatches(ipatch)->vecCollisions.size(); icoll++)
		vecPatches(ipatch)->vecCollisions[icoll]->collide(params,vecPatches(ipatch)->vecSpecies,itime);
        timer[8].update();
        
        // apply the PIC method
        // --------------------
        // for all particles of all species (see dynamic in Species.cpp)
        // (1) interpolate the fields at the particle position
        // (2) move the particle
        // (3) calculate the currents (charge conserving method)

	/*******************************************/
	/********** Move particles *****************/
	/*******************************************/
#pragma omp parallel shared (time_dual,smpiData,params, vecPatches, simWindow)
        {
	    timer[1].restart();
            if (diag_flag){
                #pragma omp for schedule(static)
                for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
		    vecPatches(ipatch)->EMfields->restartRhoJs();
            }
            #pragma omp for schedule(static)
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
		vecPatches(ipatch)->EMfields->restartRhoJ();

            #pragma omp for schedule(runtime)
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
	      vecPatches(ipatch)->dynamics(time_dual, params, simWindow, diag_flag); // include test -> Add , vecPatches(ipatch)->vecPartWall in call to Species::dynamics
	    }
	    timer[1].update();

	    timer[8].restart();
	    for (unsigned int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size(); ispec++) {
		if ( vecPatches(0)->vecSpecies[ispec]->isProj(time_dual, simWindow) ){
		    vecPatches.exchangeParticles(ispec, params, smpiData ); // Included sort_part
		}
	    }
	    timer[8].update();


	    /*******************************************/
	    /*********** Sum densities *****************/
	    /*******************************************/
	    timer[4].restart();
	    if  (diag_flag){
                #pragma omp for
		for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		    vecPatches(ipatch)->EMfields->computeTotalRhoJ(); // Per species in global, Attention if output -> Sync / per species fields
		}
	    }
	    timer[4].update();

	    timer[9].restart();
	    vecPatches.sumRhoJ( diag_flag ); // MPI

            if(diag_flag){
	        for (unsigned int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size(); ispec++) {
	            vecPatches.sumRhoJs( ispec ); // MPI
	        }
            }
	    //cout << "End sumrho" << endl;
	    timer[9].update();

	    // apply currents from antennas
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++){
		vecPatches(ipatch)->EMfields->applyAntennas(smpiData, time_dual);
        
		/*******************************************/
		/*********** Maxwell solver ****************/
		/*******************************************/
        
		// solve Maxwell's equations
		if( time_dual > params.time_fields_frozen ) {
		    timer[2].restart();
		    // saving magnetic fields (to compute centered fields used in the particle pusher)
                    #pragma omp for schedule(static)
		    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++){
			//Stores B at time n in B_m.
			vecPatches(ipatch)->EMfields->saveMagneticFields();
			// Computes Ex_, Ey_, Ez_ on all points. E is already synchronized because J has been synchronized before.
			vecPatches(ipatch)->EMfields->solveMaxwellAmpere();
		    }
		    //vecPatches.exchangeE();
                    #pragma omp for schedule(static)
		    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++){
			// Computes Bx_, By_, Bz_ at time n+1 on interior points.
			//vecPatches(ipatch)->EMfields->solveMaxwellFaraday();
			(*vecPatches(ipatch)->EMfields->MaxwellFaradaySolver_)(vecPatches(ipatch)->EMfields);
			// Applies boundary conditions on B
			vecPatches(ipatch)->EMfields->boundaryConditions(itime, time_dual, vecPatches(ipatch), params, simWindow);
		    }
		    //Synchronize B fields between patches.
		    timer[2].update();
		    timer[9].restart();
		    vecPatches.exchangeB();
		    timer[9].update();
		    timer[2].restart();
		    // Computes B at time n+1/2 using B and B_m.
                    #pragma omp for schedule(static)
		    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
			vecPatches(ipatch)->EMfields->centerMagneticFields();

		    timer[2].update();

		}
	    }
        
        // incrementing averaged electromagnetic fields
        if (vecPatches.Diags->ntime_step_avg)
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		vecPatches(ipatch)->EMfields->incrementAvgFields(itime, vecPatches.Diags->ntime_step_avg);
	    }

        // call the various diagnostics
        // ----------------------------

        #pragma omp master
        {		
        // temporary EM fields dump in Fields.h5
        timer[6].restart();
        if  (diag_flag){
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		if (ipatch==0) vecPatches(ipatch)->sio->createTimeStepInSingleFileTime( itime, vecPatches.Diags );
                vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields->allFields, itime, 0 );
		// temporary EM fields dump in Fields_avg.h5
		if  ((vecPatches.Diags->ntime_step_avg!=0) &&
		     (vecPatches.Diags->avgfieldDump_every != 0) && 
		     (itime % vecPatches.Diags->avgfieldDump_every == 0)) {
		    vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields->allFields_avg, itime, 1 );
		}
                vecPatches(ipatch)->EMfields->restartRhoJs();
            }
            diag_flag = 0 ;
	}
	timer[6].update();
 

        // run all diagnostics
        timer[3].restart();
        //for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
        //    if ( vecSpecies[ispec]->particles.isTestParticles ) {
        //        sio->writeTestParticles0(vecSpecies[ispec], params, smpi);
        //    }
        //}
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    vecPatches(ipatch)->Diags->runAllDiags(itime, vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, vecPatches(ipatch)->Interp);
	vecPatches.computeGlobalDiags(itime); // Only scalars reduction for now 
	smpiData->computeGlobalDiags( vecPatches(0)->Diags, itime); // Only scalars reduction for now 
	timer[3].update();
	} // end of pragma omp master
        

	// ----------------------------------------------------------------------
	// Validate restart  : to do
	// Restart patched moving window : to do
	checkpoint.dump(vecPatches, itime, smpiData, simWindow, params, vecPatches.Diags);
	// Break in an OpenMP region
	//if (checkpoint.dump(vecPatches, itime, smpiData, simWindow, params, vecPatches.Diags)) break; 
	 /*if  (smpiData->isMaster()){
	    if (!todump && checkpoint.dump( itime, MPI_Wtime() - starttime, params ) ){
                // Send the action to perform at next iteration
                itime2dump = itime + 1; 
                for (unsigned int islave=0; islave < mpisize; islave++) 
                    MPI_Isend(&itime2dump,1,MPI_INT,islave,0,MPI_COMM_WORLD,&action_srequests[islave]);
                todump = 1;
            }
        } else {
            MPI_Iprobe(0,0,MPI_COMM_WORLD,&todump,&action_status[0]); // waiting for a control message from master (rank=0)
            //Receive action
            if( todump ){
                MPI_Recv(&itime2dump,1,MPI_INT,0,0,MPI_COMM_WORLD,&action_status[1]);
                todump = 0;
            }
        }

        if(itime==itime2dump){
            checkpoint.dumpAll( vecPatches, itime, smpiData, simWindow, params, input_data);
            todump = 0;
	    // Warning: you can not use a break to exit an openMP structure. We have to find another way to implement the following.
            //if (params.exit_after_dump ) break;
	    }*/
	// ----------------------------------------------------------------------        

		
        timer[5].restart();
        if ( simWindow && simWindow->isMoving(time_dual) ) {
            #pragma omp single
            {
            start_moving++;
            if ((start_moving==1) && (smpiData->isMaster()) ) {
		MESSAGE(">>> Window starts moving");
            }
            }
            simWindow->operate(vecPatches, smpiData, params);
        }
        timer[5].update();


        } //End omp parallel region

	if ((itime%balancing_freq == 0)&&(smpiData->smilei_sz!=1)) {
            timer[7].restart();
            partperMPI = 0;
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++){
                for (unsigned int ispec=0 ; ispec < vecPatches(0)->vecSpecies.size() ; ispec++)
                    partperMPI += vecPatches(ipatch)->vecSpecies[ispec]->getNbrOfParticles();
            }
            partperMPI = 0;

	    smpiData->recompute_patch_count( params, vecPatches, 0. );


	    vecPatches.createPatches(params, smpiData, simWindow);
	    vecPatches.exchangePatches_new(smpiData, params);

	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++){
                for (unsigned int ispec=0 ; ispec < vecPatches(0)->vecSpecies.size() ; ispec++)
                    partperMPI += vecPatches(ipatch)->vecSpecies[ispec]->getNbrOfParticles();
            }
            npatchmoy += vecPatches.size();
            npartmoy += partperMPI;
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
    cout << "npart moy = " << npartmoy << " npatch moy = " << npatchmoy << endl;
    double coverage(0.);
    for (unsigned int i=1 ; i<timer.size() ; i++) coverage += timer[i].getTime();
    MESSAGE("Time in time loop :\t" << timer[0].getTime() << "\t"<<coverage/timer[0].getTime()*100.<< "% coverage" );
    if ( smpiData->isMaster() )
        for (unsigned int i=1 ; i<timer.size() ; i++) timer[i].print(timer[0].getTime());
    
    vecPatches.Diags->printTimers(vecPatches(0), timer[3].getTime());
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    if (latestTimeStep==params.n_time) {
	if  ( (vecPatches.Diags->fieldDump_every != 0) && (params.n_time % vecPatches.Diags->fieldDump_every != 0) )
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		if (ipatch==0) vecPatches(ipatch)->sio->createTimeStepInSingleFileTime( params.n_time, vecPatches.Diags );
		vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields->allFields, params.n_time, 0 );
	    }
	// temporary time-averaged EM fields dump in Fields_avg.h5
	if  (vecPatches.Diags->ntime_step_avg!=0) {
	    if  ( (vecPatches.Diags->avgfieldDump_every != 0) && (params.n_time % vecPatches.Diags->avgfieldDump_every != 0) )
		for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		    vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields->allFields_avg, params.n_time, 1 );
		}
	}

    }

    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    vecPatches.finalizeProbesDiags(params, stepStop);
    vecPatches.finalizeDumpFields(params, stepStop);

    for (unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++) delete vecPatches(ipatch);
    vecPatches.clear();

    if (params.nspace_win_x)
        delete simWindow;
    
    PyTools::closePython();

    TITLE("END");
    delete smpiData;
    
    return 0;
    
}//END MAIN



