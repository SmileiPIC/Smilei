
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
#include <iostream>
#include <omp.h>
#include <fstream>
#include <limits>

using namespace std;

SimWindow::SimWindow(Params& params)
{
    // ------------------------
    // Moving window parameters
    // ------------------------
    active = false;
    time_start = numeric_limits<double>::max();
    velocity_x = 1.;
    
    if( PyTools::nComponents("MovingWindow") ) {
        active = true;
        
        PyTools::extract("time_start",time_start, "MovingWindow");
        
        PyTools::extract("velocity_x",velocity_x, "MovingWindow");
    }
    
    cell_length_x_   = params.cell_length[0];
    x_moved = 0.;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    n_moved = 0 ;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    
    if( active ) {
        MESSAGE(1,"Moving window is active:");
        MESSAGE(2,"velocity_x : " << velocity_x);
        MESSAGE(2,"time_start : " << time_start);
        params.hasWindow = true;
    } else {
        params.hasWindow = false;
    }
    
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
    //int nmessage = 14+2*nSpecies;
    int nmax_laser = 4;
    int nmessage = 2*nSpecies+(2+params.nDim_particle)*vecPatches(0)->probes.size()+
        9+vecPatches(0)->EMfields->antennas.size()+4*nmax_laser;
    vector<int> nbrOfPartsSend(nSpecies,0);
    vector<int> nbrOfPartsRecv(nSpecies,0);
    
    double energy_field_lost(0.);
    vector<double> energy_part_lost( vecPatches(0)->vecSpecies.size(), 0. );

    
    // Shift the patches, new patches will be created directly with their good patchid and BC
    for (int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
        if ( vecPatches(ipatch)->isEastern() )
            for (int ispec=0 ; ispec<nSpecies ; ispec++)
                vecPatches(ipatch)->vecSpecies[ispec]->disableEast();
        vecPatches(ipatch)->neighbor_[0][1] = vecPatches(ipatch)->hindex;
        vecPatches(ipatch)->hindex = vecPatches(ipatch)->neighbor_[0][0];
    }
    // Init new patches (really new and received)
    for (int ipatch = 0 ; ipatch < nPatches ; ipatch++) {

        if ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][1] ) {
            int patchid = vecPatches(ipatch)->neighbor_[0][1]; //Because we just set neighbor_[0][1]= hindex in previous loop.
            Patch* newPatch = PatchesFactory::clone(vecPatches(0),params, smpi, patchid, n_moved );
            vecPatches.patches_.push_back( newPatch );
        }
    }

    for ( int ipatch = nPatches-1 ; ipatch >= 0 ; ipatch--) {

        // Patch à supprimer
        //if I'm western  AND I'm not a newly created patch (because we start at nPatches-1), delete me !
        if ( vecPatches(ipatch)->isWestern() ) {

            // Compute energy lost 
            energy_field_lost += vecPatches(ipatch)->EMfields->computeNRJ();
            for ( unsigned int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size() ; ispec++ )
                energy_part_lost[ispec] += vecPatches(ipatch)->vecSpecies[ispec]->computeNRJ();

            delete  vecPatches.patches_[ipatch];
            vecPatches.patches_[ipatch] = NULL;
            vecPatches.patches_.erase( vecPatches.patches_.begin() + ipatch );

        }
    }

    // Sync / Patches done for these diags -> Store in patch master 
    vecPatches(0)->EMfields->storeNRJlost( energy_field_lost );
    for ( unsigned int ispec=0 ; ispec<vecPatches(0)->vecSpecies.size() ; ispec++ )
        vecPatches(0)->vecSpecies[ispec]->storeNRJlost( energy_part_lost[ispec] );
    
    //! \todo Removed the following block because Probes are not transferred together with the patches
    //
    //// Store offset in file for current MPI process
    ////   Sould be store in the diagnostic itself
    //vector<int> offset(vecPatches(0)->probes.size());
    //for (unsigned int iprobe=0;iprobe<vecPatches(0)->probes.size();iprobe++)
    //    offset[iprobe] = vecPatches(0)->probes[iprobe]->offset_in_file;

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

    //Cut off laser before exchanging any patches to avoid deadlock
    for (int ipatch=0 ; ipatch<nPatches ; ipatch++)
        vecPatches(ipatch)->EMfields->laserDisabled();

    // Patch à envoyer
    for (int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
        //if my MPI left neighbor is not me AND I'm not a newly created patch, send me !
        if ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][0] && (int)vecPatches(ipatch)->hindex == vecPatches(ipatch)->neighbor_[0][0] ) {
            //cout << vecPatches(ipatch)->MPI_me_ << " send : " << vecPatches(ipatch)->hindex << " to " << vecPatches(ipatch)->MPI_neighbor_[0][0]<<" with tag " << vecPatches(ipatch)->hindex*nmessage << endl;
            smpi->isend( vecPatches(ipatch), vecPatches(ipatch)->MPI_neighbor_[0][0], vecPatches(ipatch)->hindex*nmessage, params );
        }
    }
    // Patch à recevoir
    for (int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
        //if my MPI right neighbor is not me AND my MPI right neighbor exists AND I am a newly created patch, I receive !
        if ( ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][1] ) && ( vecPatches(ipatch)->MPI_neighbor_[0][1] != MPI_PROC_NULL )  && (vecPatches(ipatch)->neighbor_[0][0] != (int)vecPatches(ipatch)->hindex) ){
            //cout << vecPatches(ipatch)->MPI_me_ << " recv : " << vecPatches(ipatch)->hindex << " from " << vecPatches(ipatch)->MPI_neighbor_[0][1] <<" with tag " << vecPatches(ipatch)->hindex*nmessage << endl;
            smpi->recv( vecPatches(ipatch), vecPatches(ipatch)->MPI_neighbor_[0][1], vecPatches(ipatch)->hindex*nmessage, params );
        }
    }

    //Wait for all send to be completed by the receivers too.
    //MPI_Barrier(MPI_COMM_WORLD);
    smpi->barrier();

    // Suppress after exchange to not distrub patch position during exchange
    for ( int ipatch = nPatches-1 ; ipatch >= 0 ; ipatch--) {
        if ( vecPatches(ipatch)->MPI_me_ != vecPatches(ipatch)->MPI_neighbor_[0][0] && (int)vecPatches(ipatch)->hindex == vecPatches(ipatch)->neighbor_[0][0] ) {

            delete vecPatches.patches_[ipatch];
            vecPatches.patches_[ipatch] = NULL;
            vecPatches.patches_.erase( vecPatches.patches_.begin() + ipatch );

        }

    }
    nPatches = vecPatches.size();


    // Finish shifting the patches, new patches will be created directly with their good patches
    for (int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
        if (vecPatches(ipatch)->neighbor_[0][0] != (int)vecPatches(ipatch)->hindex) continue;
            
        //For now also need to update neighbor_, corner_neighbor and their MPI counterparts even if these will be obsolete eventually.
        //vecPatches(ipatch)->corner_neighbor_[1][0]= vecPatches(ipatch)->neighbor_[1][0]; //useless
        vecPatches(ipatch)->neighbor_[1][0]=        vecPatches(ipatch)->corner_neighbor_[0][0];

        //vecPatches(ipatch)->corner_neighbor_[1][1]= vecPatches(ipatch)->neighbor_[1][1]; //useless
        vecPatches(ipatch)->neighbor_[1][1]=        vecPatches(ipatch)->corner_neighbor_[0][1];

        if (params.nDim_field == 3) {
            vecPatches(ipatch)->neighbor_[2][0]=    vecPatches(ipatch)->corner_neighbor_[1][0];
            vecPatches(ipatch)->neighbor_[2][1]=    vecPatches(ipatch)->corner_neighbor_[1][1];
        }

        //Compute missing part of the new neighborhood tables.
        vecPatches(ipatch)->Pcoordinates[0]--;

        
        if (params.nDim_field == 2) {

            int xcall = vecPatches(ipatch)->Pcoordinates[0]-1;
            int ycall = vecPatches(ipatch)->Pcoordinates[1]-1;
            if (params.bc_em_type_x[0]=="periodic" && xcall <0) xcall += (1<<params.mi[0]);
            if (params.bc_em_type_y[0]=="periodic" && ycall <0) ycall += (1<<params.mi[1]);
            vecPatches(ipatch)->corner_neighbor_[0][0] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, ycall);
            ycall = vecPatches(ipatch)->Pcoordinates[1];
            vecPatches(ipatch)->neighbor_[0][0] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, vecPatches(ipatch)->Pcoordinates[1]);
            ycall = vecPatches(ipatch)->Pcoordinates[1]+1;
            if (params.bc_em_type_y[0]=="periodic" && ycall >= 1<<params.mi[1]) ycall -= (1<<params.mi[1]);
            vecPatches(ipatch)->corner_neighbor_[0][1] = generalhilbertindex(params.mi[0] , params.mi[1], xcall, ycall);

        } else if (params.nDim_field == 3) {

            int xcall = vecPatches(ipatch)->Pcoordinates[0]-1;
            int ycall = vecPatches(ipatch)->Pcoordinates[1];
            int zcall = vecPatches(ipatch)->Pcoordinates[2];
            vecPatches(ipatch)->corner_neighbor_[0][0] = generalhilbertindex(params.mi[0] , params.mi[1], params.mi[2], xcall, ycall-1, zcall);
            vecPatches(ipatch)->neighbor_[0][0] =        generalhilbertindex(params.mi[0] , params.mi[1], params.mi[2], xcall, ycall, zcall);
            vecPatches(ipatch)->corner_neighbor_[0][1] = generalhilbertindex(params.mi[0] , params.mi[1], params.mi[2], xcall, ycall+1, zcall);
            vecPatches(ipatch)->corner_neighbor_[1][0] = generalhilbertindex(params.mi[0] , params.mi[1], params.mi[2], xcall, ycall, zcall-1);
            vecPatches(ipatch)->corner_neighbor_[1][1] = generalhilbertindex(params.mi[0] , params.mi[1], params.mi[2], xcall, ycall, zcall+1);

        }
        
    }
    
    for (int ipatch=0 ; ipatch<nPatches ; ipatch++){
        vecPatches(ipatch)->updateMPIenv(smpi);
        if ( vecPatches(ipatch)->isWestern() )
            for (int ispec=0 ; ispec<nSpecies ; ispec++)
                vecPatches(ipatch)->vecSpecies[ispec]->setWestBoundaryCondition(); 
    }
    
    vecPatches.set_refHindex() ;
    vecPatches.update_field_list() ;
    
    // \todo Temporary change: instead of moving probes, we re-create them.
    // This should allow load balancing and moving window to work for probes.
    //
    //vecPatches.move_probes(params, x_moved);
    //for (unsigned int iprobe=0;iprobe<vecPatches(0)->probes.size();iprobe++)
    //    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
    //        vecPatches(ipatch)->probes[iprobe]->offset_in_file = offset[iprobe];
    //
    for (unsigned int idiag = 0 ; idiag < vecPatches.localDiags.size() ; idiag++) {
        DiagnosticProbes* diagProbes = dynamic_cast<DiagnosticProbes*>(vecPatches.localDiags[idiag]);
        if ( diagProbes ) {
            diagProbes->patchesHaveMoved = true;
            diagProbes->x_moved = x_moved;
        }
    }
    
    return;
    
}

bool SimWindow::isMoving(double time_dual)
{
    return ((time_dual - time_start)*velocity_x > x_moved);
}

void SimWindow::operate_arnaud(VectorPatch& vecPatches, SmileiMPI* smpi, Params& params)
{
    int xcall, ycall, h0;
    double energy_field_lost(0.);
    vector<double> energy_part_lost( vecPatches(0)->vecSpecies.size(), 0. );
    Patch* mypatch;
    int tid(0), nthds(1);// Rneighbor, Lneighbor;
    #ifdef _OPENMP
        tid = omp_get_thread_num();
        nthds = omp_get_num_threads();
    #endif

    //Initialization for inter-process communications
    h0 = vecPatches(0)->hindex;
    int nPatches = vecPatches.size();
    int nSpecies( vecPatches(0)->vecSpecies.size() );
    int nmessage = 2*nSpecies+14;
    std::vector<Patch*> delete_patches_, update_patches_;

    // Inits in a single to minize sync
    //#pragma omp single
    {
        for (unsigned int i=0; i< nthds; i++)
            patch_to_be_created[i].clear();
        vecPatches_old.resize(nPatches);
        x_moved += cell_length_x_*params.n_space[0];
        n_moved += params.n_space[0];
        //Handle diags
        vecPatches.closeAllDiags(smpi);
    }

    //#pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++)
        vecPatches_old[ipatch] = vecPatches(ipatch);

    //#pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
         mypatch = vecPatches_old[ipatch];

        //If my right neighbor does not belong to me ...
        if (mypatch->MPI_neighbor_[0][1] != mypatch->MPI_me_){
            // Store it as a patch to be created later.
            patch_to_be_created[tid].push_back(ipatch); //(shared omp vector of int)
        }

        //If my left neighbor does not belong to me ...
        if (mypatch->MPI_neighbor_[0][0] != mypatch->MPI_me_) {
            delete_patches_.push_back(mypatch); // Stores pointers to patches to be deleted later (private omp vector of pointers to patches)
            //... I might have to MPI send myself to the left...
            if (mypatch->MPI_neighbor_[0][0] != MPI_PROC_NULL){
                //Left neighbour to which the patch should be sent to.
                //Lneighbor = mypatch->MPI_neighbor_[0][0];
                smpi->isend( mypatch, mypatch->MPI_neighbor_[0][0] , (mypatch->neighbor_[0][0]) * nmessage );
            }
         } else { //In case my left neighbor does belong to me:
            // I become my left neighbor.
            //Nothing to do on global indexes or min_locals. The Patch structure remains the same and unmoved.
            //Update hindex and coordinates.

            mypatch->Pcoordinates[0] -= 1;
            mypatch->neighbor_[0][1] =  mypatch->hindex;
            mypatch->hindex = mypatch->neighbor_[0][0];
	    mypatch->tmp_neighbor_[0][0] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[0][0];
	    mypatch->tmp_neighbor_[1][0] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[1][0];
	    mypatch->tmp_neighbor_[1][1] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[1][1];
            mypatch->MPI_neighbor_[0][1] = mypatch->MPI_me_ ;
	    mypatch->tmp_MPI_neighbor_[0][0] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[0][0];
	    mypatch->tmp_MPI_neighbor_[1][0] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[1][0];
	    mypatch->tmp_MPI_neighbor_[1][1] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[1][1];
            update_patches_.push_back(mypatch); // Stores pointers to patches to be deleted later (private omp vector of pointers to patches)

            //And finally put the patch at the correct rank in vecPatches.
            vecPatches.patches_[mypatch->hindex - h0 ] = mypatch ; 
             
       }

    }//End loop on Patches. This barrier matters.

    //Creation of new Patches if necessary
    //The "new" operator must be included in a single area otherwise conflicts arise for unknown reasons.
    //#pragma omp single
    {
         //#pragma omp for schedule(static)
         for (unsigned int i=0; i<nthds; i++){
             for (unsigned int j=0; j< patch_to_be_created[i].size(); j++){
                 //vecPatches.patches_[patch_to_be_created[i][j]] = new Patch(params, laser_params, smpi, h0 + patch_to_be_created[i][j], n_moved);
                 vecPatches.patches_[patch_to_be_created[i][j]] = PatchesFactory::create(params, smpi, h0 + patch_to_be_created[i][j], n_moved );
                 //stores all indices of patch_to_be_created in a single vector.
                 if (i>0) patch_to_be_created[0].push_back(patch_to_be_created[i][j]);
             }
         }
    } // This barrier is important.


    //Initialization of new Patches if necessary.
    //#pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < patch_to_be_created[0].size() ; ipatch++) {
         mypatch = vecPatches(patch_to_be_created[0][ipatch]);
         //Rneighbor = mypatch->MPI_neighbor_[0][1];
         //If I receive something from my right neighbour:
         if (mypatch->MPI_neighbor_[0][1] != MPI_PROC_NULL)
             smpi->recv( mypatch, mypatch->MPI_neighbor_[0][1], (mypatch->Hindex())*nmessage, params );
         // And else, nothing to do.
    } //This barrier matters. 

    

    //#pragma omp for schedule(static)
    for (int ipatch=0 ; ipatch<nPatches ; ipatch++ )
	vecPatches(ipatch)->EMfields->laserDisabled();


    //#pragma omp single
    { 
        vecPatches.openAllDiags(params,smpi);

        vecPatches.set_refHindex() ;
        vecPatches.update_field_list() ;
    }
    //#pragma omp for schedule(static)
    for (int j=0; j < update_patches_.size(); j++){
        mypatch = update_patches_[j];
	mypatch->MPI_neighbor_[0][0] = mypatch->tmp_MPI_neighbor_[0][0];
	mypatch->MPI_neighbor_[1][0] = mypatch->tmp_MPI_neighbor_[1][0];
	mypatch->MPI_neighbor_[1][1] = mypatch->tmp_MPI_neighbor_[1][1];
	mypatch->neighbor_[0][0] = mypatch->tmp_neighbor_[0][0];
	mypatch->neighbor_[1][0] = mypatch->tmp_neighbor_[1][0];
	mypatch->neighbor_[1][1] = mypatch->tmp_neighbor_[1][1];
    }
    smpi->barrier();
    //Each thread erases data of sent patches
    for (int j=0; j < delete_patches_.size(); j++){
        mypatch = delete_patches_[j];

       	energy_field_lost += mypatch->EMfields->computeNRJ();
        for ( int ispec=0 ; ispec<nSpecies ; ispec++ )
            energy_part_lost[ispec] += mypatch->vecSpecies[ispec]->computeNRJ();

        delete  mypatch;
    }

}
