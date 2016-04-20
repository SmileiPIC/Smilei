
#include "SimWindow.h"
#include "Params.h"
#include "Species.h"
#include "ElectroMagn.h"
#include "Interpolator.h"
#include "Projector.h"
#include "SmileiMPI.h"
#include "VectorPatch.h"
#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "DiagsVectorPatch.h"
#include <iostream>
#include <omp.h>
#include <fstream>

using namespace std;

SimWindow::SimWindow(Params& params)
{
    nspace_win_x_ = params.nspace_win_x;
    cell_length_x_   = params.cell_length[0];
    x_moved = 0.;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    n_moved = 0;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    vx_win_ = params.vx_win; 
    int nthds(1);
    #ifdef _OPENMP
        nthds = omp_get_max_threads();
    #endif
    patch_to_be_created.resize(nthds);

    t_move_win_ = params.t_move_win;
    MESSAGE(1,"Moving window is active:");
    MESSAGE(2,"nspace_win_x_ : " << nspace_win_x_);
    MESSAGE(2,"cell_length_x_ : " << cell_length_x_);
    MESSAGE(2,"vx_win_ : " << vx_win_);
    MESSAGE(2,"t_move_win_ : " << t_move_win_);
    
}

SimWindow::~SimWindow()
{
}


void SimWindow::operate(VectorPatch& vecPatches, SmileiMPI* smpi, Params& params)
{

        x_moved += cell_length_x_*params.n_space[0];
        n_moved += params.n_space[0];

    // Store current number of patch on current MPI process
    // Don't move during this process
    int nPatches( vecPatches.size() );
    int nSpecies  ( vecPatches(0)->vecSpecies.size() );
    int nmessage = 10+2*nSpecies;
    vector<int> nbrOfPartsSend(nSpecies,0);
    vector<int> nbrOfPartsRecv(nSpecies,0);

    // Delete western patch
    double energy_field_lost(0.);
    vector<double> energy_part_lost( vecPatches(0)->vecSpecies.size(), 0. );


    vecPatches.closeAllDiags(smpi);


    hid_t globalFile    = vecPatches(0)->sio->global_file_id_;
    hid_t globalFileAvg = vecPatches(0)->sio->global_file_id_avg;


    // Shift the patches, new patches will be created directly with their good patchid
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
	vecPatches(ipatch)->neighbor_[0][1] = vecPatches(ipatch)->hindex;
        vecPatches(ipatch)->hindex = vecPatches(ipatch)->neighbor_[0][0];
    }
    // Init new patches (really new and received)
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {

        if ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][1] ) {
            int patchid = vecPatches(ipatch)->neighbor_[0][1];
            Patch* newPatch = PatchesFactory::create(params, smpi, patchid, n_moved );
            vecPatches.patches_.push_back( newPatch );
        }
    }

    for ( int ipatch = nPatches-1 ; ipatch >= 0 ; ipatch--) {

        // Patch à supprimer
        if ( vecPatches(ipatch)->isWestern() ) {

	    // Compute energy lost 
	    energy_field_lost += vecPatches(ipatch)->EMfields->computeNRJ();
	    for ( int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size() ; ispec++ )
		energy_part_lost[ispec] += vecPatches(ipatch)->vecSpecies[ispec]->computeNRJ();

            vecPatches(ipatch)->sio->setFiles(0,0);
            delete  vecPatches.patches_[ipatch];
            vecPatches.patches_[ipatch] = NULL;
	    vecPatches.patches_.erase( vecPatches.patches_.begin() + ipatch );

        }
    }

    // Sync / Patches done for these diags -> Store in patch master 
    vecPatches(0)->EMfields->storeNRJlost( energy_field_lost );
    for ( int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size() ; ispec++ )
	vecPatches(0)->vecSpecies[ispec]->storeNRJlost( energy_part_lost[ispec] );

    nPatches = vecPatches.size();

    // Sort patch by hindex (to avoid deadlock)
    //bool stop;
    int jpatch(nPatches-1);
    do {
        for ( int ipatch = 0 ; ipatch<jpatch ; ipatch++  ) {
            if ( vecPatches(ipatch)->hindex > vecPatches(jpatch)->hindex ) {
                Patch* tmp = vecPatches(ipatch);
                vecPatches.patches_[ipatch] = vecPatches.patches_[jpatch];
        	vecPatches.patches_[jpatch] = tmp;
            }
        }
        jpatch--;
    } while(jpatch>=0);

            // Patch à envoyer
            for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
                //if my MPI left neighbor is not me AND I'm not a newly created patch, send me !
                if ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][0] && vecPatches(ipatch)->hindex == vecPatches(ipatch)->neighbor_[0][0] ) {
                    smpi->isend( vecPatches(ipatch), vecPatches(ipatch)->MPI_neighbor_[0][0], vecPatches(ipatch)->hindex*nmessage );
		    //cout << vecPatches(ipatch)->MPI_me_ << " send : " << vecPatches(ipatch)->vecSpecies[0]->getNbrOfParticles() << " & " << vecPatches(ipatch)->vecSpecies[1]->getNbrOfParticles() << endl;
                }
            }
            // Patch à recevoir
            for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
                //if my MPI right neighbor is not me AND my MPI right neighbor exists AND I am a newly created patch, I receive !
                if ( ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][1] ) && ( vecPatches(ipatch)->MPI_neighbor_[0][1] != MPI_PROC_NULL )  && (vecPatches(ipatch)->neighbor_[0][0] != vecPatches(ipatch)->hindex) ){
                    smpi->recv( vecPatches(ipatch), vecPatches(ipatch)->MPI_neighbor_[0][1], vecPatches(ipatch)->hindex*nmessage, params );
		    //cout << vecPatches(ipatch)->MPI_me_ << " recv : " << vecPatches(ipatch)->vecSpecies[0]->getNbrOfParticles() << " & " << vecPatches(ipatch)->vecSpecies[1]->getNbrOfParticles() << endl;
                }
            }

    //Wait for all send to be completed by the receivers too.
    //MPI_Barrier(MPI_COMM_WORLD);
    smpi->barrier();

    // Suppress after exchange to not distrub patch position during exchange
    for ( int ipatch = nPatches-1 ; ipatch >= 0 ; ipatch--) {
        if ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][0] && vecPatches(ipatch)->hindex == vecPatches(ipatch)->neighbor_[0][0] ) {

            vecPatches(ipatch)->sio->setFiles(0,0);
            delete vecPatches.patches_[ipatch];
            vecPatches.patches_[ipatch] = NULL;
	    vecPatches.patches_.erase( vecPatches.patches_.begin() + ipatch );

        }

    }
    nPatches = vecPatches.size();


    // Finish shifting the patches, new patches will be created directly with their good patches
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
	if (vecPatches(ipatch)->neighbor_[0][0] != vecPatches(ipatch)->hindex) continue;
	    
	//Compute missing part of the new neighborhood tables.
	vecPatches(ipatch)->Pcoordinates[0]--;

	int xcall = vecPatches(ipatch)->Pcoordinates[0]-1;
	int ycall = vecPatches(ipatch)->Pcoordinates[1]-1;
	if (params.bc_em_type_x[0]=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
	if (params.bc_em_type_y[0]=="periodic" && ycall <0) ycall += (1<<params.mi[1]);
	vecPatches(ipatch)->neighbor_[1][0] = generalhilbertindex(params.mi[0] , params.mi[1], vecPatches(ipatch)->Pcoordinates[0], ycall)       ;
	vecPatches(ipatch)->neighbor_[0][0] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, vecPatches(ipatch)->Pcoordinates[1]);
	ycall = vecPatches(ipatch)->Pcoordinates[1]+1;
	if (params.bc_em_type_y[0]=="periodic" && ycall >= 1<<params.mi[1]) ycall -= (1<<params.mi[1]);
	vecPatches(ipatch)->neighbor_[1][1] = generalhilbertindex(params.mi[0] , params.mi[1], vecPatches(ipatch)->Pcoordinates[0], ycall) ;
	
    }

    for (int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
        vecPatches(ipatch)->updateMPIenv(smpi);
    }

    for (int ipatch=0 ; ipatch<nPatches ; ipatch++ )
	vecPatches(ipatch)->EMfields->laserDisabled();


    // 
    vecPatches.openAllDiags(params,smpi);

    //vecPatches.definePatchDiagsMaster();
    DiagsVectorPatch::definePatchDiagsMaster( vecPatches, globalFile, globalFileAvg ); // sio
    DiagsVectorPatch::updatePatchFieldDump( vecPatches, params );                      // sio

    vecPatches.set_refHindex() ;
    vecPatches.update_field_list() ;

    return;

}



bool SimWindow::isMoving(double time_dual)
{
    return ( (nspace_win_x_) && ((time_dual - t_move_win_)*vx_win_ > x_moved) );
}

void SimWindow::setOperators(VectorPatch& vecPatches)
{

    for (unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++) {

	vecPatches(ipatch)->updateMvWinLimits( x_moved, n_moved );

	for (unsigned int ispec=0 ; ispec<vecPatches(ipatch)->vecSpecies.size(); ispec++) {
	    vecPatches(ipatch)->vecSpecies[ispec]->updateMvWinLimits(x_moved);
	}

	vecPatches(ipatch)->Interp->setMvWinLimits( vecPatches(ipatch)->getCellStartingGlobalIndex(0) );
	vecPatches(ipatch)->Proj->setMvWinLimits  ( vecPatches(ipatch)->getCellStartingGlobalIndex(0) );

    }


}

void SimWindow::operate_arnaud(VectorPatch& vecPatches, SmileiMPI* smpi, Params& params)
{
  int xcall, ycall, h0;
    Patch* mypatch;
    int tid(0), nthds(1), tag, Rneighbor, Lneighbor;
    #ifdef _OPENMP
        tid = omp_get_thread_num();
        nthds = omp_get_num_threads();
    #endif

    //Initialization for inter-process communications

    h0 = vecPatches(0)->hindex;
    int nSpecies( vecPatches(0)->vecSpecies.size() );
    int nDim_Parts( vecPatches(0)->vecSpecies[0]->particles->dimension() );
    vector<int> nbrOfPartsSend(nSpecies,0);
    vector<int> nbrOfPartsRecv(nSpecies,0);
    vector < vector<int> > store_npart_recv;
    std::vector<Patch*> send_patches_;
    MPI_Request srequest[8+3*nSpecies];//Number of calls made to MPI_Isend for each patch exchanged.
    //vector <MPI_Request*> srequests;

    #pragma omp for schedule(static)
    for (unsigned int i=0; i< nthds; i++){
        patch_to_be_created[i].clear();
    }
    vecPatches_old.resize(vecPatches.size());

    #pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++) {
        vecPatches_old[ipatch] = vecPatches(ipatch);
    } //Barrier at the end of this omp for is important to prevent an update of x_moved before resolution of isMoving in the main loop.
    #pragma omp single
    {
        x_moved += cell_length_x_*params.n_space[0];
        n_moved += params.n_space[0];
    }

    #pragma omp for schedule(runtime)
    for (unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++) {
         mypatch = vecPatches_old[ipatch];

        //If my right neighbor does not belong to me ...
        if (mypatch->MPI_neighbor_[0][1] != mypatch->MPI_me_)
            // Store it as a patch to be created later.
            patch_to_be_created[tid].push_back(ipatch); //(shared omp vector of int)

        //If my left neighbor does not belong to me ...
        if (mypatch->MPI_neighbor_[0][0] != mypatch->MPI_me_) {
            //... I might have to MPI send myself to the left...
            if (mypatch->MPI_neighbor_[0][0] != MPI_PROC_NULL){
                send_patches_.push_back(mypatch); // Stores pointers to patches to be sent in send_patches_ (private omp vector of pointers to patches)
                //Tag is the index of the patch after reception (left neighbour index because it is sent to the left).
                tag = mypatch->neighbor_[0][0];
                for (int ispec=0 ; ispec<nSpecies ; ispec++) {
	            nbrOfPartsSend[ispec] = mypatch->vecSpecies[ispec]->getNbrOfParticles();
	        }
                store_npart_recv.push_back(nbrOfPartsSend);
                //Left neighbour to which the patch should be sent to.
                Lneighbor = mypatch->MPI_neighbor_[0][0];
                //Sends the number of particles first
       //         srequests.push_back(srequest);
//	        smpi->send( store_npart_recv.back(), Lneighbor, tag, &srequests.back()[0] );
                //Then sends the patch
//	        smpi->send( mypatch, Lneighbor, tag, &srequests.back()[1] );
            } else {
            //Erase my data right away. Data of sent patches will be erased later (after the receive has been completed)
                //mypatch->Diags->probes.setFile(0);
                //mypatch->sio->setFiles(0,0);
                //delete (mypatch);
                //for (unsigned int ispec=0 ; ispec<mypatch->vecSpecies.size(); ispec++) delete (mypatch->vecSpecies[ispec]);
	        //mypatch->vecSpecies.clear();
                //delete (mypatch->EMfields);
                //delete (mypatch->Interp);
                //delete (mypatch->Proj);

                // Compute energy lost 
	        energy_field_lost += mypatch->EMfields->computeNRJ();
	        for ( int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size() ; ispec++ )
	            energy_part_lost[ispec] += mypatch->vecSpecies[ispec]->computeNRJ();

                mypatch->sio->setFiles(0,0);
                delete  vecPatches_old.patches_[ipatch];
                vecPatches.patches_[ipatch] = NULL; //Not necessary but for safety measures.

            }
        } else { //In case my left neighbor does belong to me:
            // I become my left neighbor.
            //Nothing to do on global indexes or min_locals. The Patch structure remains the same and unmoved.
            //Update hindex and coordinates.
            mypatch->neighbor_[0][1] =  mypatch->hindex;
            mypatch->hindex = mypatch->neighbor_[0][0];
            mypatch->Pcoordinates[0] -= 1;

            //Compute missing part of the new neighborhood tables.
            xcall = mypatch->Pcoordinates[0]-1;
            ycall = mypatch->Pcoordinates[1]-1;
            if (params.bc_em_type_x[0]=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
            if (params.bc_em_type_y[0]=="periodic" && ycall <0) ycall += (1<<params.mi[1]);
	    mypatch->neighbor_[0][0] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, mypatch->Pcoordinates[1]);
	    mypatch->neighbor_[1][0] = generalhilbertindex(params.mi[0] , params.mi[1], mypatch->Pcoordinates[0], ycall );
            ycall = mypatch->Pcoordinates[1]+1;
            if (params.bc_em_type_y[0]=="periodic" && ycall >= 1<<params.mi[1]) ycall -= (1<<params.mi[1]);
	    mypatch->neighbor_[1][1] = generalhilbertindex(params.mi[0] , params.mi[1], mypatch->Pcoordinates[0], ycall );
            
            //Reflect the changes on MPI neighbors
            mypatch->MPI_neighbor_[0][0] = smpi->hrank(mypatch->neighbor_[0][0]);
            mypatch->MPI_neighbor_[1][0] = smpi->hrank(mypatch->neighbor_[1][0]);
            mypatch->MPI_neighbor_[1][1] = smpi->hrank(mypatch->neighbor_[1][1]);

            //And finally put the patch at the correct rank in vecPatches.
            vecPatches.patches_[mypatch->hindex - h0 ] = mypatch ; 
             
       }

    }//End loop on Patches. This barrier matters.

    //Creation of new Patches if necessary
    //The "new" operator must be included in a single area otherwise conflicts arise for unknown reasons.
    #pragma omp single
    {
         for (unsigned int i=0; i<nthds; i++){
             for (unsigned int j=0; j< patch_to_be_created[i].size(); j++){
                 vecPatches.patches_[patch_to_be_created[i][j]] = new Patch(params, laser_params, smpi, h0 + patch_to_be_created[i][j], n_moved);
                 //stores all indices of patch_to_be_created in a single vector.
                 if (i>0) patch_to_be_created[0].push_back(patch_to_be_created[i][j]);
             }
         }
    } // This barrier is important.

    //Initialization of new Patches if necessary.
    #pragma omp for schedule(runtime)
    for (unsigned int ipatch = 0 ; ipatch < patch_to_be_created[0].size() ; ipatch++) {
         mypatch = vecPatches(patch_to_be_created[0][ipatch]);
         Rneighbor = mypatch->MPI_neighbor_[0][1];
         //If I receive something from my right neighbour:
         if (Rneighbor != MPI_PROC_NULL){
             smpi->recv( &nbrOfPartsRecv, Rneighbor, mypatch->Hindex() );
	     for (int ispec=0 ; ispec<nSpecies ; ispec++)
	         mypatch->vecSpecies[ispec]->particles->initialize( nbrOfPartsRecv[ispec], nDim_Parts );
             smpi->recv( mypatch, Rneighbor, mypatch->Hindex() );
         }
         // And else, nothing to do.
    } //This barrier matters. 

    //Each thread erases data of sent patches
    for (int j= send_patches_.size()-1; j>=0; j--){
        mypatch = send_patches_[j];
        for (unsigned int ispec=0 ; ispec<mypatch->vecSpecies.size(); ispec++) delete (mypatch->vecSpecies[ispec]);
	mypatch->vecSpecies.clear();
        delete (mypatch->EMfields);
        delete (mypatch->Interp);
        delete (mypatch->Proj);
        
        //send_patches_[j]->Diags->probes.setFile(0);
        //send_patches_[j]->sio->setFiles(0,0);
        //delete send_patches_[j];
    }
}

