
#include "SimWindow.h"
#include "Params.h"
#include "Species.h"
#include "ElectroMagn.h"
#include "Interpolator.h"
#include "Projector.h"
#include "SmileiMPI.h"
#include "VectorPatch.h"
#include "DiagnosticProbes.h"
#include "DiagnosticTrack.h"
#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include <iostream>
#include <omp.h>
#include <fstream>
#include <limits>
#include "ElectroMagnBC_Factory.h"

using namespace std;

SimWindow::SimWindow(Params& params)
{

    // ------------------------
    // Moving window parameters
    // ------------------------
    active = false;
    time_start = numeric_limits<double>::max();
    velocity_x = 1.;

#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#else
    max_threads = 1;
#endif
    patch_to_be_created.resize(max_threads);

    if( PyTools::nComponents("MovingWindow") ) {
        active = true;

        PyTools::extract("time_start",time_start, "MovingWindow");

        PyTools::extract("velocity_x",velocity_x, "MovingWindow");
    }

    cell_length_x_   = params.cell_length[0];
    x_moved = 0.;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    n_moved = 0 ;      //The window has not moved at t=0. Warning: not true anymore for restarts.


    if( active ) {
        if (velocity_x != 0. && params.EM_BCs[0][0] == "periodic")
            ERROR("Periodic topology in the moving window direction is neither encouraged nor supported");

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

bool SimWindow::isMoving(double time_dual)
{
    return active && ((time_dual - time_start)*velocity_x > x_moved);
}

void SimWindow::operate(VectorPatch& vecPatches, SmileiMPI* smpi, Params& params, unsigned int itime, double time_dual)
{
    if ( ! isMoving(time_dual) ) return;

    unsigned int h0;
    double energy_field_lost(0.);
    std::vector<double> energy_part_lost( vecPatches(0)->vecSpecies.size(), 0. );
    Patch* mypatch;

    //Initialization for inter-process communications
    h0 = vecPatches(0)->hindex;
    unsigned int nPatches = vecPatches.size();
    unsigned int nSpecies( vecPatches(0)->vecSpecies.size() );
    int nmessage( vecPatches.nrequests );

    std::vector<Patch*> delete_patches_, update_patches_, send_patches_;

#ifdef _OPENMP
    int my_thread = omp_get_thread_num();
#else
    int my_thread = 0;
#endif
    (patch_to_be_created[my_thread]).clear();


    #pragma omp single
    {
        if (n_moved == 0) MESSAGE(">>> Window starts moving");

        vecPatches_old.resize(nPatches);
        n_moved += params.n_space[0];
    }
    //Cut off laser before exchanging any patches to avoid deadlock and store pointers in vecpatches_old.
    #pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++){
        vecPatches_old[ipatch] = vecPatches(ipatch);
        vecPatches(ipatch)->EMfields->laserDisabled();
    }


    #pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
         mypatch = vecPatches_old[ipatch];
        //If my right neighbor does not belong to me store it as a patch to
        // be created later.
        if (mypatch->MPI_neighbor_[0][1] != mypatch->MPI_me_){
            (patch_to_be_created[my_thread]).push_back(ipatch);
        }

        // Do not sent Xmax conditions
        if ( mypatch->isXmax() && mypatch->EMfields->emBoundCond[1] )
            mypatch->EMfields->emBoundCond[1]->disableExternalFields();

        //If my left neighbor does not belong to me ...
        if (mypatch->MPI_neighbor_[0][0] != mypatch->MPI_me_) {
            delete_patches_.push_back(mypatch); // Stores pointers to patches to be deleted later
            //... I might have to MPI send myself to the left...
            if (mypatch->MPI_neighbor_[0][0] != MPI_PROC_NULL){
                send_patches_.push_back(mypatch); // Stores pointers to patches to be sent later
                smpi->isend( vecPatches_old[ipatch], vecPatches_old[ipatch]->MPI_neighbor_[0][0] , (vecPatches_old[ipatch]->neighbor_[0][0]) * nmessage, params );
            }
        } else { //In case my left neighbor belongs to me:
            // I become my left neighbor.
            //Update hindex and coordinates.

            if ( mypatch->isXmax() )
                for (unsigned int ispec=0 ; ispec<nSpecies ; ispec++)
                    mypatch->vecSpecies[ispec]->disableXmax();
            mypatch->Pcoordinates[0] -= 1;
            mypatch->neighbor_[0][1] =  mypatch->hindex;
            mypatch->hindex = mypatch->neighbor_[0][0];
            mypatch->MPI_neighbor_[0][1] = mypatch->MPI_me_ ;
            //stores indices in tmp buffers so that original values can be read by other patches.
            mypatch->tmp_neighbor_[0][0] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[0][0];
            mypatch->tmp_MPI_neighbor_[0][0] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[0][0];
            for (unsigned int idim = 1; idim < params.nDim_field ; idim++){
                mypatch->tmp_neighbor_[idim][0] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[idim][0];
                mypatch->tmp_neighbor_[idim][1] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[idim][1];
                mypatch->tmp_MPI_neighbor_[idim][0] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[idim][0];
                mypatch->tmp_MPI_neighbor_[idim][1] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[idim][1];
            }
            update_patches_.push_back(mypatch); // Stores pointers to patches that will need to update some neighbors from tmp_neighbors.

            //And finally put the patch at the correct rank in vecPatches.
            vecPatches.patches_[mypatch->hindex - h0 ] = mypatch ;

        }
    }//End loop on Patches. This barrier matters.
    // At this point, all isends have been done and the list of patches to delete at the end is complete.
    // The lists of patches to create and patches to update is also complete.

    //Creation of new Patches
    for (unsigned int j = 0; j < patch_to_be_created[my_thread].size();  j++){
        //create patch without particle.
        #pragma omp critical
        mypatch = PatchesFactory::clone(vecPatches(0),params, smpi, vecPatches.domain_decomposition_, h0 + patch_to_be_created[my_thread][j], n_moved, false );

        // Do not receive Xmin condition
        if ( mypatch->isXmin() && mypatch->EMfields->emBoundCond[0] )
            mypatch->EMfields->emBoundCond[0]->disableExternalFields();

        mypatch->finalizeMPIenvironment(params);
        //Position new patch
        vecPatches.patches_[patch_to_be_created[my_thread][j]] = mypatch ;
        //Receive Patch if necessary
        if (mypatch->MPI_neighbor_[0][1] != MPI_PROC_NULL){
            smpi->recv( mypatch, mypatch->MPI_neighbor_[0][1], (mypatch->hindex)*nmessage, params );
            patch_to_be_created[my_thread][j] = nPatches ; //Mark no needs of particles
        }

        // Create Xmin condition which could not be received
        if ( mypatch->isXmin() ){
            for (auto& embc:mypatch->EMfields->emBoundCond) {
                if (embc) delete embc;
            }
            mypatch->EMfields->emBoundCond = ElectroMagnBC_Factory::create(params, mypatch);
            mypatch->EMfields->laserDisabled();
        }

        mypatch->EMfields->laserDisabled();
        mypatch->EMfields->updateGridSize(params, mypatch);
    }

    //Wait for sends to be completed

    #pragma omp for schedule(static)
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++){
        if (vecPatches_old[ipatch]->MPI_neighbor_[0][0] !=  vecPatches_old[ipatch]->MPI_me_ && vecPatches_old[ipatch]->MPI_neighbor_[0][0] != MPI_PROC_NULL){
            smpi->waitall( vecPatches_old[ipatch] );
        }
    }

    //Update the correct neighbor values
    for (unsigned int j=0; j < update_patches_.size(); j++){
        mypatch = update_patches_[j];
        mypatch->MPI_neighbor_[0][0] = mypatch->tmp_MPI_neighbor_[0][0];
        mypatch->neighbor_[0][0] = mypatch->tmp_neighbor_[0][0];
        for (unsigned int idim = 1; idim < params.nDim_field ; idim++){
            mypatch->MPI_neighbor_[idim][0] = mypatch->tmp_MPI_neighbor_[idim][0];
            mypatch->MPI_neighbor_[idim][1] = mypatch->tmp_MPI_neighbor_[idim][1];
            mypatch->neighbor_[idim][0] = mypatch->tmp_neighbor_[idim][0];
            mypatch->neighbor_[idim][1] = mypatch->tmp_neighbor_[idim][1];
        }

        mypatch->updateTagenv(smpi);
        if ( mypatch->isXmin() ){
            for (unsigned int ispec=0 ; ispec<nSpecies ; ispec++)
                mypatch->vecSpecies[ispec]->setXminBoundaryCondition();
        }
        if (mypatch->has_an_MPI_neighbor())
            mypatch->createType(params);
        else

            mypatch->cleanType();

        if ( mypatch->isXmin() ){
            for (auto& embc:mypatch->EMfields->emBoundCond) {
                if (embc) delete embc;
            }
            mypatch->EMfields->emBoundCond = ElectroMagnBC_Factory::create(params, mypatch);
            mypatch->EMfields->laserDisabled();
        }
        if ( mypatch->wasXmax( params ) ){
            for (auto& embc:mypatch->EMfields->emBoundCond) {
                if (embc) delete embc;
            }
            mypatch->EMfields->emBoundCond = ElectroMagnBC_Factory::create(params, mypatch);
            mypatch->EMfields->laserDisabled();
            mypatch->EMfields->updateGridSize(params, mypatch);

        }
    }

    //Wait for sends to be completed
    for (unsigned int j=0; j < send_patches_.size(); j++){
            smpi->waitall( send_patches_[j] );
    }

    #pragma omp barrier

    //Fill necessary patches with particles
    #pragma omp master
    {
        for (int ithread=0; ithread < max_threads ; ithread++){
            for (unsigned int j=0; j< (patch_to_be_created[ithread]).size(); j++){
                //If new particles are required
                if (patch_to_be_created[ithread][j] != nPatches){
                    mypatch = vecPatches.patches_[patch_to_be_created[ithread][j]];
                    for (unsigned int ispec=0 ; ispec<nSpecies ; ispec++)
                        mypatch->vecSpecies[ispec]->createParticles(params.n_space, params, mypatch, 0 );
                    // We define the IDs of the new particles
                    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ )
                        if( DiagnosticTrack* track = dynamic_cast<DiagnosticTrack*>(vecPatches.localDiags[idiag]) )
                            track->setIDs( mypatch );
                }
            }
        }
    }

    #pragma omp single nowait
    {
        x_moved += cell_length_x_*params.n_space[0];
        vecPatches.update_field_list() ;
        //update list fields for species diag too ??

        // Tell that the patches moved this iteration (needed for probes)
        vecPatches.lastIterationPatchesMoved = itime;
    }

    std::vector<double> poynting[2];
    poynting[0].resize(params.nDim_field,0.0);
    poynting[1].resize(params.nDim_field,0.0);

    //Delete useless patches
    for (unsigned int j=0; j < delete_patches_.size(); j++){
        mypatch = delete_patches_[j];

        //if (mypatch->isXmin()) {
        //    energy_field_lost += mypatch->EMfields->computeNRJ();
        //    for ( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ )
        //        energy_part_lost[ispec] += mypatch->vecSpecies[ispec]->computeNRJ();
        //}

        for (unsigned int jp=0; jp<2;jp++) //directions (xmin/xmax, ymin/ymax, zmin/zmax)
            for (unsigned int i=0 ; i<params.nDim_field ; i++){ //axis 0=x, 1=y, 2=z
                poynting[jp][i] += mypatch->EMfields->poynting[jp][i];
            }


        delete  mypatch;
    }

    // SUM energy_field_lost, energy_part_lost and poynting / All threads
    #pragma omp critical
    {
        vecPatches(0)->EMfields->storeNRJlost( energy_field_lost );
        for ( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ )
            vecPatches(0)->vecSpecies[ispec]->storeNRJlost( energy_part_lost[ispec] );

        for (unsigned int j=0; j<2;j++) //directions (xmin/xmax, ymin/ymax, zmin/zmax)
            for (unsigned int i=0 ; i< params.nDim_field ; i++) //axis 0=x, 1=y, 2=z
                vecPatches(0)->EMfields->poynting[j][i] += poynting[j][i];
    }



}
