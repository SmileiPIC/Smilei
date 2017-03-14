
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

bool SimWindow::isMoving(double time_dual)
{
    return ((time_dual - time_start)*velocity_x > x_moved);
}

void SimWindow::operate(VectorPatch& vecPatches, SmileiMPI* smpi, Params& params)
{
    unsigned int h0;
    double energy_field_lost(0.);
    std::vector<double> energy_part_lost( vecPatches(0)->vecSpecies.size(), 0. );
    std::vector<unsigned int> patch_to_be_created;
    Patch* mypatch;
    
    //Initialization for inter-process communications
    h0 = vecPatches(0)->hindex;
    unsigned int nPatches = vecPatches.size();
    unsigned int nSpecies( vecPatches(0)->vecSpecies.size() );
    std::vector<Patch*> delete_patches_, update_patches_, send_patches_;
    int nmessage( vecPatches.nrequests );
    
    vecPatches_old.resize(nPatches);
    x_moved += cell_length_x_*params.n_space[0];
    n_moved += params.n_space[0];
    
    //Cut off laser before exchanging any patches to avoid deadlock and store pointers in vecpatches_old.
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++){
        vecPatches_old[ipatch] = vecPatches(ipatch);
        vecPatches(ipatch)->EMfields->laserDisabled();
    }
    
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++) {
         mypatch = vecPatches_old[ipatch];
    
        //If my right neighbor does not belong to me store it as a patch to be created later.
        if (mypatch->MPI_neighbor_[0][1] != mypatch->MPI_me_)
            patch_to_be_created.push_back(ipatch); 
    
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
            for (unsigned int idim = 1; idim < params.nDim_particle ; idim++){
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
    
    //Creation of new Patches if necessary
    //Use clone instead of create ??
    //These patches are created with correct parameters.
    for (unsigned int j=0; j< patch_to_be_created.size(); j++){
    //for (int j=1; j >= 0 ; j--){
        mypatch = PatchesFactory::clone(vecPatches(0),params, smpi, h0 + patch_to_be_created[j], n_moved );
        if (mypatch->MPI_neighbor_[0][1] != MPI_PROC_NULL){
            smpi->recv( mypatch, mypatch->MPI_neighbor_[0][1], (mypatch->hindex)*nmessage, params );
        }
        else { // Must force particles creation, see in SpeciesFactory :
            // if (params.restart)
            //     thisSpecies->particles->initialize( 0, params.nDim_particle );
            if (params.restart)
                for (unsigned int ispec=0 ; ispec<nSpecies ; ispec++)
                    mypatch->vecSpecies[ispec]->createParticles(params.n_space, params, mypatch, 0 );
            // We define the IDs of the new particles
            for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ )
                if( DiagnosticTrack* track = dynamic_cast<DiagnosticTrack*>(vecPatches.localDiags[idiag]) )
                    track->setIDs( mypatch );
        }
        mypatch->EMfields->laserDisabled();
        vecPatches.patches_[patch_to_be_created[j]] = mypatch ;
    }
    
    //Update the correct neighbor values
    for (unsigned int j=0; j < update_patches_.size(); j++){
        mypatch = update_patches_[j];
        mypatch->MPI_neighbor_[0][0] = mypatch->tmp_MPI_neighbor_[0][0];
        mypatch->neighbor_[0][0] = mypatch->tmp_neighbor_[0][0];
        for (unsigned int idim = 1; idim < params.nDim_particle ; idim++){
            mypatch->MPI_neighbor_[idim][0] = mypatch->tmp_MPI_neighbor_[idim][0];
            mypatch->MPI_neighbor_[idim][1] = mypatch->tmp_MPI_neighbor_[idim][1];
            mypatch->neighbor_[idim][0] = mypatch->tmp_neighbor_[idim][0];
            mypatch->neighbor_[idim][1] = mypatch->tmp_neighbor_[idim][1];
        }
    }
    
    //Wait for sends to be completed
    for (unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++){ 
        if (vecPatches_old[ipatch]->MPI_neighbor_[0][0] !=  vecPatches_old[ipatch]->MPI_me_ && vecPatches_old[ipatch]->MPI_neighbor_[0][0] != MPI_PROC_NULL){
            smpi->waitall( vecPatches_old[ipatch] );
        }
    }
    
    for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++){
        vecPatches(ipatch)->updateTagenv(smpi);
        if ( vecPatches(ipatch)->isXmin() ){
            for (unsigned int ispec=0 ; ispec<nSpecies ; ispec++)
                vecPatches(ipatch)->vecSpecies[ispec]->setXminBoundaryCondition(); 
        }
        if (vecPatches(ipatch)->has_an_MPI_neighbor())
            vecPatches(ipatch)->createType(params);
        else
            vecPatches(ipatch)->cleanType();
    }
    
    //Should be useless
    vecPatches.update_field_list() ;
    //update list fields for species diag too ??
    
    for (unsigned int idiag = 0 ; idiag < vecPatches.localDiags.size() ; idiag++) {
        DiagnosticProbes* diagProbes = dynamic_cast<DiagnosticProbes*>(vecPatches.localDiags[idiag]);
        if ( diagProbes ) {
            diagProbes->patchesHaveMoved = true;
            diagProbes->x_moved = x_moved;
        }
    }
    
    //Delete useless patches
    for (unsigned int j=0; j < delete_patches_.size(); j++){
        mypatch = delete_patches_[j];
        
        energy_field_lost += mypatch->EMfields->computeNRJ();
        for ( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ )
            energy_part_lost[ispec] += mypatch->vecSpecies[ispec]->computeNRJ();
        
        delete  mypatch;
    }

}
