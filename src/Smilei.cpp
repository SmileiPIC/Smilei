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

#include "PatchesFactory.h"
#include "Checkpoint.h"

#include "DiagParams.h"

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
    SmileiMPI* smpi = NULL;
    Checkpoint checkpoint(params, diag_params);

#ifdef _OPENMP
    if (smpiData->isMaster())
        MESSAGE("\tOpenMP : Number of thread per MPI process : " << omp_get_max_threads() );
#else
    if (smpiData->isMaster()) MESSAGE("\tOpenMP : Disabled");
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
    vector<Species*> vecSpecies(0);
    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // ------------------------------------------------------------------------

    unsigned int stepStart=0, stepStop=params.n_time;
	
    // time at integer time-steps (primal grid)
    double time_prim = stepStart * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (stepStart +0.5) * params.timestep;
    // Do we initially do diags or not ?
    int diag_flag = 1;

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
    ElectroMagn* EMfields = NULL;
    
    // interpolation operator (virtual)
    Interpolator* Interp = NULL;
    
    // projection operator (virtual)
    Projector* Proj = NULL;
    
    // Create diagnostics
    Diagnostic *Diags = NULL;

    
    VectorPatch vecPatches = PatchesFactory::createVector(params, diag_params, laser_params, smpiData);
    vecPatches.initProbesDiags(params, diag_params, 0);
    vecPatches.initDumpFields(params, diag_params, 0);
    
    // reading from dumped file the restart values
    if (params.restart) {
        MESSAGE(1, "READING fields and particles for restart");
        DEBUG(vecSpecies.size());
        checkpoint.restartAll( vecPatches, stepStart, smpiData, simWindow, params, input_data);

        double restart_time_dual = (stepStart +0.5) * params.timestep;
        // A revoir !
	//if ( simWindow && ( simWindow->isMoving(restart_time_dual) ) ) {
	//    simWindow->setOperators(vecSpecies, Interp, Proj, smpiData);
	//    simWindow->operate(vecSpecies, EMfields, Interp, Proj, smpiData , params);
	//}
        //smpiData->recompute_patch_count( params, vecPatches, restart_time_dual );
	
    } else {
        // Initialize the electromagnetic fields
        // -----------------------------------
        // Init rho and J by projecting all particles of subdomain
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
	    vecPatches(ipatch)->EMfields->restartRhoJs();
	    vecPatches(ipatch)->dynamics(time_dual, params, simWindow, diag_flag); //include test
	}
	for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
	    if ( vecPatches(0)->vecSpecies[ispec]->isProj(time_dual, simWindow) )
		vecPatches.exchangeParticles(ispec, params, smpiData ); // Included sort_part
	}
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) 
	    vecPatches(ipatch)->EMfields->computeTotalRhoJ(); // Per species in global, Attention if output -> Sync / per species fields
	vecPatches.sumRhoJ( diag_flag ); // MPI
	for (unsigned int ispec=0 ; ispec<params.n_species; ispec++)
	    vecPatches.sumRhoJs( ispec ); // MPI
        diag_flag = 0;

 
#ifdef _PATCHINGROGRESS
        // Init electric field (Ex/1D, + Ey/2D)
	/*if (!EMfields->isRhoNull(smpi)) {
	    MESSAGE("----------------------------------------------");
	    MESSAGE("Solving Poisson at time t = 0");
	    MESSAGE("----------------------------------------------");    
	    EMfields->solvePoisson(smpi);
	}*/
	vecPatches.solvePoisson( params, smpiData );
#endif
  
	
        MESSAGE("----------------------------------------------");
        MESSAGE("Running diags at time t = 0");
        MESSAGE("----------------------------------------------");
        // run diagnostics at time-step 0	
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    vecPatches(ipatch)->Diags->runAllDiags(0, vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, vecPatches(ipatch)->Interp, smpi);
	vecPatches.computeGlobalDiags(0);
	smpiData->computeGlobalDiags( vecPatches(0)->Diags, 0);
	
	
        for (unsigned int ispec=0 ; ispec<params.n_species; ispec++)
	  MESSAGE(1,"Species " << ispec << " (" << params.species_param[ispec].species_type << ") created with " << (int)vecPatches(0)->Diags->getScalar("N_"+params.species_param[ispec].species_type) << " particles" );
#ifdef _DEBUGPATCH
	for (int ipatch = 0 ; ipatch<vecPatches.size() ; ipatch++)
	    cout << (int)vecPatches(ipatch)->vecSpecies[0]->getNbrOfParticles() << " particles on " << vecPatches(ipatch)->Hindex() << endl;
#endif
	
        // temporary EM fields dump in Fields.h5
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
	    if (ipatch==0) vecPatches(ipatch)->sio->createTimeStepInSingleFileTime( 0, diag_params );
	    vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields, 0 );
	    // temporary EM fields dump in Fields_avg.h5
	    if (diag_params.ntime_step_avg!=0)
		vecPatches(ipatch)->sio->writeAvgFieldsSingleFileTime( vecPatches(ipatch)->EMfields, 0 );
	    vecPatches(ipatch)->EMfields->restartRhoJs();
	}
	diag_flag = 0 ;
    }

	
    // Count timer
    int ntimer(7);
    Timer timer[ntimer];
    timer[0].init(smpiData, "global");
    timer[1].init(smpiData, "particles");
    timer[2].init(smpiData, "maxwell");
    timer[3].init(smpiData, "diagnostics");
    timer[4].init(smpiData, "densities");
    timer[5].init(smpiData, "Mov window");
    timer[6].init(smpiData, "diag fields");

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
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    MESSAGE("Time-Loop is started: number of time-steps n_time = " << params.n_time);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
	
    for (unsigned int itime=stepStart+1 ; itime <= stepStop ; itime++) {

        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        if  ((diag_params.fieldDump_every != 0) && (itime % diag_params.fieldDump_every == 0)) diag_flag = 1;

        // send message at given time-steps
        // --------------------------------
        timer[0].update();
        
        //double timElapsed=smpiData->time_seconds();
	if ( (itime % diag_params.print_every == 0) &&  ( smpiData->isMaster() ) ) {
            MESSAGE(1,"t = "          << setw(7) << setprecision(2)   << time_dual/params.conv_fac
                    << "   it = "       << setw(log10(params.n_time)+1) << itime  << "/" << params.n_time
                    << "   sec = "      << setw(7) << setprecision(2)   << timer[0].getTime() 
                    << "   E = "        << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("Etot") 
		    << "   Epart = "        << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("Eparticles")
		    << "   EFields = "        << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("EFields")
//		    << "   Elost = "        << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("Elost") 
                  << "   E_bal(%) = " << setw(6) << std::fixed << setprecision(2) 
		    << 100.0*vecPatches(0)->Diags->getScalar("Ebal_norm")
		    );
	    if (simWindow) 
		MESSAGE(1, "\t\t MW Elost = " << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("Emw_lost")
			<< "     MW Eadd  = " << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("Emw_part")
			<< "     MW Elost (fields) = " << std::scientific << setprecision(4)<< vecPatches(0)->Diags->getScalar("Emw_lost_fields")
			<< setw(6) << std::fixed << setprecision(2) );
	}

        // put density and currents to 0 + save former density
        // ---------------------------------------------------
        
        
        // apply the PIC method
        // --------------------
        // for all particles of all species (see dynamic in Species.cpp)
        // (1) interpolate the fields at the particle position
        // (2) move the particle
        // (3) calculate the currents (charge conserving method)

	/*******************************************/
	/********** Move particles *****************/
	/*******************************************/
#pragma omp parallel shared (EMfields,time_dual,smpi,params, vecPatches, simWindow)
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

            //cout << "Starting dynamics" << endl;
            #pragma omp for schedule(runtime)
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		vecPatches(ipatch)->dynamics(time_dual, params, simWindow, diag_flag); // include test
	    }
	    // Inter Patch exchange
                //for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
	        //    if ( vecPatches(0)->vecSpecies[ispec]->isProj(time_dual, simWindow) ){
	        //        vecPatches.exchangeParticles(ispec, params); 
                //    }
                //}
                for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
	            if ( vecPatches(0)->vecSpecies[ispec]->isProj(time_dual, simWindow) ){
	                vecPatches.exchangeParticles(ispec, params, smpiData ); // Included sort_part
                            //if (itime%200 == 0) {
	            	    //    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	            	    //        vecPatches(ipatch)->vecSpecies[ispec]->count_sort_part(params);
                            //}
	            }
	        }

	    timer[1].update();


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
        //cout << "End computeTrho" << endl;
        //Synchronize J and posisbly Rho between patches.
	vecPatches.sumRhoJ( diag_flag ); // MPI
	#pragma omp master
        {

            if(diag_flag){
	        for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
	            vecPatches.sumRhoJs( ispec ); // MPI
	        }
            }
        }
        //cout << "End sumrho" << endl;
        timer[4].update();


	/*******************************************/
	/*********** Maxwell solver ****************/
	/*******************************************/
        
        // solve Maxwell's equations
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
	    vecPatches(ipatch)->EMfields->solveMaxwellFaraday();
            // Applies boundary conditions on B
	    vecPatches(ipatch)->EMfields->boundaryConditions(itime, time_dual, vecPatches(ipatch), params, simWindow);
        }
        //Synchronize B fields between patches.
	vecPatches.exchangeB();
        // Computes B at time n+1/2 using B and B_m.
        #pragma omp for schedule(static)
	for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    vecPatches(ipatch)->EMfields->centerMagneticFields();

        timer[2].update();

        
#ifdef _TOBEPATCHED
        // incrementing averaged electromagnetic fields
        if (diag_params.ntime_step_avg) EMfields->incrementAvgFields(itime, diag_params.ntime_step_avg);
#endif        
        // call the various diagnostics
        // ----------------------------

        #pragma omp master
        {		
        // temporary EM fields dump in Fields.h5
        timer[6].restart();
        if  (diag_flag){
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		if (ipatch==0) vecPatches(ipatch)->sio->createTimeStepInSingleFileTime( itime, diag_params );
                vecPatches(ipatch)->sio->writeAllFieldsSingleFileTime( vecPatches(ipatch)->EMfields, itime );
		// temporary EM fields dump in Fields_avg.h5
		if (diag_params.ntime_step_avg!=0) //  if ((diag_params.avgfieldDump_every != 0) && (itime % diag_params.avgfieldDump_every == 0))
		    vecPatches(ipatch)->sio->writeAvgFieldsSingleFileTime( vecPatches(ipatch)->EMfields, 0 );
                vecPatches(ipatch)->EMfields->restartRhoJs();
            }
            diag_flag = 0 ;
            }
	timer[6].update();

        // run all diagnostics
        timer[3].restart();
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	    vecPatches(ipatch)->Diags->runAllDiags(itime, vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, vecPatches(ipatch)->Interp, smpiData);
	vecPatches.computeGlobalDiags(itime); // Only scalars reduction for now 
	smpiData->computeGlobalDiags( vecPatches(0)->Diags, itime); // Only scalars reduction for now 
	timer[3].update();
        }

	// ----------------------------------------------------------------------
	// Validate restart  : to do
	// Restart patched moving window : to do
        //if  (smpiData->isMaster()){
	//    if (!todump && checkpoint.dump( itime, MPI_Wtime() - starttime, params ) ){
        //        // Send the action to perform at next iteration
        //        itime2dump = itime + 1; 
        //        for (unsigned int islave=0; islave < mpisize; islave++) 
        //            MPI_Isend(&itime2dump,1,MPI_INT,islave,0,MPI_COMM_WORLD,&action_srequests[islave]);
        //        todump = 1;
        //    }
        //} else {
        //    MPI_Iprobe(0,0,MPI_COMM_WORLD,&todump,&action_status[0]); // waiting for a control message from master (rank=0)
        //    //Receive action
        //    if( todump ){
        //        MPI_Recv(&itime2dump,1,MPI_INT,0,0,MPI_COMM_WORLD,&action_status[1]);
        //        todump = 0;
        //    }
        //}

        //if(itime==itime2dump){
        //    checkpoint.dumpAll( vecPatches, itime, smpiData, simWindow, params, input_data);
        //    todump = 0;
        //    if (params.exit_after_dump ) break;
        //}
	// ----------------------------------------------------------------------        

		
//#ifdef _INPROGRESS
        timer[5].restart();
        if ( simWindow && simWindow->isMoving(time_dual) ) {
            start_moving++;
            if ((start_moving==1) && (smpiData->isMaster()) ) {
		MESSAGE(">>> Window starts moving");
            }
            simWindow->operate(vecPatches, smpiData, params, diag_params, laser_params);
        }
        timer[5].update();
//#endif


        } //End omp parallel region

#ifdef _TESTPATCHEXCH
	if (itime==61) {
	    /*if (smpiData->getRank()==0)
		vecPatches.send_patch_id_.push_back(3);
	    else if (smpiData->getRank()==1) 
	    vecPatches.recv_patch_id_.push_back(3);*/

	    //smpiData->patch_count[0] = 5;
	    //smpiData->patch_count[1] = 3;
	    smpiData->recompute_patch_count( params, vecPatches, 0. );
	    smpiData->patch_count[0] = 5;
	    smpiData->patch_count[1] = 3;


	    //cout << "patch_count modified" << endl;
	    vecPatches.createPatches(params, diag_params, laser_params, smpiData, simWindow);
	    //cout << "patch created" << endl;
	    vecPatches.setNbrParticlesToExch(smpiData);
	    //cout << "nbr particles exchanged" << endl;
	    vecPatches.exchangePatches(smpiData);
	    //cout << "data exchanged" << endl;

	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
		cout << smpiData->getRank() << " run patchId = " << vecPatches(ipatch)->Hindex() << endl;
	    
	}
#endif
    }//END of the time loop
    
    smpiData->barrier();
    
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    MESSAGE("End time loop, time dual = " << time_dual/params.conv_fac);
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    
    //double timElapsed=smpiData->time_seconds();
    //if ( smpiData->isMaster() ) MESSAGE(0, "Time in time loop : " << timElapsed );
    timer[0].update();
    MESSAGE(0, "Time in time loop : " << timer[0].getTime() );
    if ( smpiData->isMaster() )
        for (int i=1 ; i<ntimer ; i++) timer[i].print(timer[0].getTime());
    
    double coverage(0.);
    for (int i=1 ; i<ntimer ; i++) coverage += timer[i].getTime();
    MESSAGE(0, "\t" << setw(12) << "Coverage\t" << coverage/timer[0].getTime()*100. << " %" );
    
    
    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------
    
    // temporary EM fields dump in Fields.h5
//    if  ( (diag_params.fieldDump_every != 0) && (params.n_time % diag_params.fieldDump_every != 0) )
//        sio->writeAllFieldsSingleFileTime( EMfields, params.n_time );
//    // temporary time-averaged EM fields dump in Fields_avg.h5
//    if  (diag_params.ntime_step_avg!=0)
//        if  ( (diag_params.avgfieldDump_every != 0) && (params.n_time % diag_params.avgfieldDump_every != 0) )
//            sio->writeAvgFieldsSingleFileTime( EMfields, params.n_time );

    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    vecPatches.finalizeProbesDiags(params, diag_params, stepStop);
    vecPatches.finalizeDumpFields(params, diag_params, stepStop);

    for (unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++) delete vecPatches(ipatch);
    vecPatches.clear();


    if (Proj) delete Proj;
    if (Interp) delete Interp;
    if (EMfields) delete EMfields;
    if (Diags) delete Diags;

    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    
    MESSAGE("-----------------------------------------------------------------------------------------------------");
    MESSAGE("END " << namelist);
    MESSAGE("-----------------------------------------------------------------------------------------------------");

    if (smpi) delete smpi;
    delete smpiData;
    if (params.nspace_win_x)
        delete simWindow;
    
    return 0;
    
}//END MAIN



