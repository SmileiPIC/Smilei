
#include "SimWindow.h"
#include "Params.h"
#include "Species.h"
#include "SpeciesVAdaptiveMixedSort.h"
#include "SpeciesVAdaptive.h"
#include "SpeciesV.h"
#include "ElectroMagn.h"
#include "Interpolator.h"
#include "Projector.h"
#include "SmileiMPI.h"
#include "VectorPatch.h"
#include "Region.h"
#include "DiagnosticProbes.h"
#include "DiagnosticTrack.h"
#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include <iostream>
#include <omp.h>
#include <fstream>
#include <limits>
#include "ElectroMagnBC_Factory.h"
#include "EnvelopeBC_Factory.h"
#include "DoubleGrids.h"
#include "SyncVectorPatch.h"

using namespace std;

SimWindow::SimWindow( Params &params )
{

    // ------------------------
    // Moving window parameters
    // ------------------------
    active = false;
    time_start = numeric_limits<double>::max();
    velocity_x = 1.;
    number_of_additional_shifts = 0;
    additional_shifts_time = 0.;
    
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#else
    max_threads = 1;
#endif
    patch_to_be_created.resize( max_threads );
    patch_particle_created.resize( max_threads );
    
    if( PyTools::nComponents( "MovingWindow" ) ) {
        active = true;
        
        TITLE( "Initializing moving window" );
        
        PyTools::extract( "time_start", time_start, "MovingWindow"  );
        PyTools::extract( "velocity_x", velocity_x, "MovingWindow"  );
        PyTools::extract( "number_of_additional_shifts", number_of_additional_shifts, "MovingWindow"  );
        PyTools::extract( "additional_shifts_time", additional_shifts_time, "MovingWindow"  );
    }
    
    cell_length_x_   = params.cell_length[0];
    n_space_x_       = params.patch_size_[0];
    additional_shifts_iteration = floor(additional_shifts_time / params.timestep + 0.5);
    x_moved = 0.;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    n_moved = 0 ;      //The window has not moved at t=0. Warning: not true anymore for restarts.
    
    
    if( active ) {
        //if( velocity_x != 0. && params.EM_BCs[0][0] == "periodic" ) {
        //    ERROR( "Periodic topology in the moving window direction is neither encouraged nor supported" );
        //}
        
        MESSAGE( 1, "Moving window is active:" );
        MESSAGE( 2, "velocity_x : " << velocity_x );
        MESSAGE( 2, "time_start : " << time_start );
        if (number_of_additional_shifts > 0){
            MESSAGE( 2, "number_of_additional_shifts : " << number_of_additional_shifts );
            MESSAGE( 2, "additional_shifts_time : " << additional_shifts_time );
        }
        params.hasWindow = true;
    } else {
        params.hasWindow = false;
    }
    
}

SimWindow::~SimWindow()
{
}

bool SimWindow::isMoving( double time_dual )
{
    return active && ( ( time_dual - time_start )*velocity_x > x_moved - number_of_additional_shifts*cell_length_x_*n_space_x_*(time_dual>additional_shifts_time) );
}

void SimWindow::shift( VectorPatch &vecPatches, SmileiMPI *smpi, Params &params, unsigned int itime, double time_dual, Region& region )
{
    if( ! isMoving( time_dual ) && itime != additional_shifts_iteration ) {
        return;
    }
 
    unsigned int h0;
    Patch *mypatch;
    
    //Initialization for inter-process communications
    h0 = vecPatches( 0 )->hindex;
    unsigned int nPatches = vecPatches.size();
    unsigned int nSpecies( vecPatches( 0 )->vecSpecies.size() );
    int nmessage( vecPatches.nrequests );
    
    std::vector<Patch *> delete_patches_, update_patches_, send_patches_;
    
#ifdef _OPENMP
    int my_thread = omp_get_thread_num();
#else
    int my_thread = 0;
#endif
    
#ifdef _NO_MPI_TM
    #pragma omp master
    {
#endif
    
        ( patch_to_be_created[my_thread] ).clear();
        ( patch_particle_created[my_thread] ).clear();
        
#ifndef _NO_MPI_TM
        #pragma omp single
#endif
        {
            if( n_moved == 0 ) {
                MESSAGE( ">>> Window starts moving" );
            }

            vecPatches_old.resize( nPatches );
            n_moved += params.patch_size_[0];
        }
        //Cut off laser before exchanging any patches to avoid deadlock and store pointers in vecpatches_old.
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#endif
        for( unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++ ) {
            vecPatches_old[ipatch] = vecPatches( ipatch );
            vecPatches( ipatch )->EMfields->laserDisabled();
        }
        
        
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static) private(mypatch)
#endif
        for( unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++ ) {
            mypatch = vecPatches_old[ipatch];
            //If my right neighbor does not belong to me store it as a patch to
            // be created later.
            if( mypatch->MPI_neighbor_[0][1] != mypatch->MPI_me_ ) {
                ( patch_to_be_created[my_thread] ).push_back( ipatch );
                ( patch_particle_created[my_thread] ).push_back( true );
            }
            
            // Do not sent Xmax conditions
            if( mypatch->isXmax() && mypatch->EMfields->emBoundCond[1] ) {
                mypatch->EMfields->emBoundCond[1]->disableExternalFields();
            }
            
            //If my left neighbor does not belong to me ...
            if( mypatch->MPI_neighbor_[0][0] != mypatch->MPI_me_ ) {
                delete_patches_.push_back( mypatch ); // Stores pointers to patches to be deleted later
                //... I might have to MPI send myself to the left...
                if( mypatch->MPI_neighbor_[0][0] != MPI_PROC_NULL ) {
                    if ( vecPatches_old[ipatch]->Pcoordinates[0]!=0 ) {
                        send_patches_.push_back( mypatch ); // Stores pointers to patches to be sent later
                        int Href_receiver = 0;
                        for (int irk = 0; irk < mypatch->MPI_neighbor_[0][0]; irk++) Href_receiver += smpi->patch_count[irk];
                        // The tag is the patch number in the receiver vector of patches in order to avoid too large tags not supported by some MPI versions.
                        smpi->isend( vecPatches_old[ipatch], vecPatches_old[ipatch]->MPI_neighbor_[0][0], ( vecPatches_old[ipatch]->neighbor_[0][0] - Href_receiver ) * nmessage, params, false );
                    }
                }
            } else { //In case my left neighbor belongs to me:
                // I become my left neighbor.
                //Update hindex and coordinates.
                
                if( mypatch->isXmax() )
                    for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                        mypatch->vecSpecies[ispec]->disableXmax();
                    }
                mypatch->Pcoordinates[0] -= 1;
                mypatch->neighbor_[0][1] =  mypatch->hindex;
                mypatch->hindex = mypatch->neighbor_[0][0];
                mypatch->MPI_neighbor_[0][1] = mypatch->MPI_me_ ;
                //stores indices in tmp buffers so that original values can be read by other patches.
                mypatch->tmp_neighbor_[0][0] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[0][0];
                mypatch->tmp_MPI_neighbor_[0][0] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[0][0];
                for( unsigned int idim = 1; idim < params.nDim_field ; idim++ ) {
                    mypatch->tmp_neighbor_[idim][0] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[idim][0];
                    mypatch->tmp_neighbor_[idim][1] = vecPatches_old[mypatch->hindex - h0 ]->neighbor_[idim][1];
                    mypatch->tmp_MPI_neighbor_[idim][0] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[idim][0];
                    mypatch->tmp_MPI_neighbor_[idim][1] = vecPatches_old[mypatch->hindex - h0 ]->MPI_neighbor_[idim][1];
                }
                update_patches_.push_back( mypatch ); // Stores pointers to patches that will need to update some neighbors from tmp_neighbors.
                
                //And finally put the patch at the correct rank in vecPatches.
                vecPatches.patches_[mypatch->hindex - h0 ] = mypatch ;
                
            }
        }//End loop on Patches. This barrier matters.
        // At this point, all isends have been done and the list of patches to delete at the end is complete.
        // The lists of patches to create and patches to update is also complete.
        
        //Creation of new Patches
        for( unsigned int j = 0; j < patch_to_be_created[my_thread].size();  j++ ) {
            //create patch without particle.
#ifndef _NO_MPI_TM
            #pragma omp critical
#endif
                mypatch = PatchesFactory::clone( vecPatches( 0 ), params, smpi, vecPatches.domain_decomposition_, h0 + patch_to_be_created[my_thread][j], n_moved, false );
            
            // Do not receive Xmin condition
            if( mypatch->isXmin() && mypatch->EMfields->emBoundCond[0] ) {
                mypatch->EMfields->emBoundCond[0]->disableExternalFields();
            }
            
            mypatch->finalizeMPIenvironment( params );
            //Position new patch
            vecPatches.patches_[patch_to_be_created[my_thread][j]] = mypatch ;
            //Receive Patch if necessary
            if( mypatch->MPI_neighbor_[0][1] != MPI_PROC_NULL ) {
                if ( mypatch->Pcoordinates[0]!=params.number_of_patches[0]-1 ) {
                    // The tag is the patch number in the receiver vector of patches in order to avoid too large tags not supported by some MPI versions.
                    smpi->recv( mypatch, mypatch->MPI_neighbor_[0][1], ( mypatch->hindex - vecPatches.refHindex_ )*nmessage, params, false );
                    patch_particle_created[my_thread][j] = false ; //Mark no needs of particles
                }
            }
            
            // Create Xmin condition which could not be received
            if( mypatch->isXmin() ) {
                for( auto &embc:mypatch->EMfields->emBoundCond ) {
                    if( embc ) {
                        delete embc;
                    }
                }
#ifndef _NO_MPI_TM
                #pragma omp critical
#endif
                    mypatch->EMfields->emBoundCond = ElectroMagnBC_Factory::create( params, mypatch );
                if (mypatch->EMfields->envelope){
                    for( auto &embc:mypatch->EMfields->envelope->EnvBoundCond ) {
                        if( embc ) {
                            delete embc;
                        }
                    }
                    mypatch->EMfields->envelope->EnvBoundCond = EnvelopeBC_Factory::create( params, mypatch );
                }

                mypatch->EMfields->laserDisabled();
                if (!params.multiple_decomposition){
                    mypatch->EMfields->emBoundCond[0]->apply(mypatch->EMfields, time_dual, mypatch);
                    if (mypatch->EMfields->envelope) mypatch->EMfields->envelope->EnvBoundCond[0]->apply(mypatch->EMfields->envelope, mypatch->EMfields, mypatch);
                }
            }
            
            mypatch->EMfields->laserDisabled();
            mypatch->EMfields->updateGridSize( params, mypatch );
        }
        
        //Wait for sends to be completed
        
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#endif
        for( unsigned int ipatch = 0 ; ipatch < nPatches ; ipatch++ ) {
            if( vecPatches_old[ipatch]->MPI_neighbor_[0][0] !=  vecPatches_old[ipatch]->MPI_me_ && vecPatches_old[ipatch]->MPI_neighbor_[0][0] != MPI_PROC_NULL ) {
                smpi->waitall( vecPatches_old[ipatch] );
            }
        }
        
        //Update the correct neighbor values
        for( unsigned int j=0; j < update_patches_.size(); j++ ) {
            mypatch = update_patches_[j];
            mypatch->MPI_neighbor_[0][0] = mypatch->tmp_MPI_neighbor_[0][0];
            mypatch->neighbor_[0][0] = mypatch->tmp_neighbor_[0][0];
            for( unsigned int idim = 1; idim < params.nDim_field ; idim++ ) {
                mypatch->MPI_neighbor_[idim][0] = mypatch->tmp_MPI_neighbor_[idim][0];
                mypatch->MPI_neighbor_[idim][1] = mypatch->tmp_MPI_neighbor_[idim][1];
                mypatch->neighbor_[idim][0] = mypatch->tmp_neighbor_[idim][0];
                mypatch->neighbor_[idim][1] = mypatch->tmp_neighbor_[idim][1];
            }
            
            if( mypatch->isXmin() ) {
                for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                    mypatch->vecSpecies[ispec]->setXminBoundaryCondition();
                }
            }
            
            if( mypatch->isXmin() ) {
                for( auto &embc:mypatch->EMfields->emBoundCond ) {
                    if( embc ) {
                        delete embc;
                    }
                }
#ifndef _NO_MPI_TM
                #pragma omp critical
#endif
                    mypatch->EMfields->emBoundCond = ElectroMagnBC_Factory::create( params, mypatch );
                if (mypatch->EMfields->envelope){
                    for( auto &embc:mypatch->EMfields->envelope->EnvBoundCond ) {
                        if( embc ) {
                            delete embc;
                        }
                    }
                    mypatch->EMfields->envelope->EnvBoundCond = EnvelopeBC_Factory::create( params, mypatch );
                }

                mypatch->EMfields->laserDisabled();
                if (!params.multiple_decomposition){
                    mypatch->EMfields->emBoundCond[0]->apply(mypatch->EMfields, time_dual, mypatch);
                    if (mypatch->EMfields->envelope) mypatch->EMfields->envelope->EnvBoundCond[0]->apply(mypatch->EMfields->envelope, mypatch->EMfields, mypatch);
                }
            }
            if( mypatch->wasXmax( params ) ) {
                for( auto &embc:mypatch->EMfields->emBoundCond ) {
                    if( embc ) {
                        delete embc;
                    }
                }
#ifndef _NO_MPI_TM
                #pragma omp critical
#endif
                    mypatch->EMfields->emBoundCond = ElectroMagnBC_Factory::create( params, mypatch );

                if (mypatch->EMfields->envelope){
                    for( auto &embc:mypatch->EMfields->envelope->EnvBoundCond ) {
                        if( embc ) {
                            delete embc;
                        }
                    }
                    mypatch->EMfields->envelope->EnvBoundCond = EnvelopeBC_Factory::create( params, mypatch );
                }

                mypatch->EMfields->laserDisabled();
                mypatch->EMfields->updateGridSize( params, mypatch );
                
            }
        }
        
        //Wait for sends to be completed
        for( unsigned int j=0; j < send_patches_.size(); j++ ) {
            smpi->waitall( send_patches_[j] );
        }
        
#ifndef _NO_MPI_TM
        #pragma omp barrier
#endif
        
        //Fill necessary patches with particles
#ifndef _NO_MPI_TM
        #pragma omp master
#endif
        {
            for( int ithread=0; ithread < max_threads ; ithread++ ) {
                for( unsigned int j=0; j< ( patch_to_be_created[ithread] ).size(); j++ ) {
                
                    // Current newly created patch
                    mypatch = vecPatches.patches_[patch_to_be_created[ithread][j]];
                    
                    // If new particles are required
                    if( patch_particle_created[ithread][j] ) {
                        for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                            ParticleCreator particle_creator;
                            particle_creator.associate(mypatch->vecSpecies[ispec]);
                            
                            // Aera for particle creation
                            struct SubSpace init_space;
                            init_space.cell_index_[0] = 0;
                            init_space.cell_index_[1] = 0;
                            init_space.cell_index_[2] = 0;
                            init_space.box_size_[0]   = params.patch_size_[0];
                            init_space.box_size_[1]   = params.patch_size_[1];
                            init_space.box_size_[2]   = params.patch_size_[2];
                            
                            particle_creator.create( init_space, params, mypatch, 0 );
                            
                        }
                        
                        mypatch->EMfields->applyExternalFields( mypatch );
                        if( params.save_magnectic_fields_for_SM ) {
                            mypatch->EMfields->saveExternalFields( mypatch );
                        }
                        
                    } // end test patch_particle_created[ithread][j]
                } // end j loop
            } // End ithread loop
        } // End omp master region
#ifndef _NO_MPI_TM
        #pragma omp barrier
#endif
        
        //Fill necessary patches with particles
        if( ( params.vectorization_mode == "on" ) ) {

            //#pragma omp master
            //{
#ifndef _NO_MPI_TM
            #pragma omp for schedule(static) private(mypatch)
#endif
            for( int ithread=0; ithread < max_threads ; ithread++ ) {
                for( unsigned int j=0; j< ( patch_to_be_created[ithread] ).size(); j++ ) {
                
                    // Current newly created patch
                    mypatch = vecPatches.patches_[patch_to_be_created[ithread][j]];
                    
                    // If new particles are required
                        for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                            mypatch->vecSpecies[ispec]->computeParticleCellKeys( params );
                            mypatch->vecSpecies[ispec]->sortParticles( params );
                        }
                } // end j loop
            } // End ithread loop
            //}
        }
        
        // Adaptive vectorization mode -- mixed sort
        else if( params.vectorization_mode == "adaptive_mixed_sort" ) {
#ifndef _NO_MPI_TM
            #pragma omp for schedule(static) private(mypatch)
#endif
            for( int ithread=0; ithread < max_threads ; ithread++ ) {
                for( unsigned int j=0; j< ( patch_to_be_created[ithread] ).size(); j++ ) {
                
                    // Current newly created patch
                    mypatch = vecPatches.patches_[patch_to_be_created[ithread][j]];
                    
                    // If new particles are required
                    if( patch_particle_created[ithread][j] ) {
                        for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                            mypatch->vecSpecies[ispec]->configuration( params, mypatch );
                        }
                    }
                    // Patches that have received particles from another patch
                    // without the creation of new particles
                    else { // (patch_particle_created[ithread][j] == false)
                        for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                            // For the adaptive vectorization, we partially reconfigure the patch
                            // We do not have to sort, but operators may have to be reconfigured
                            mypatch->vecSpecies[ispec]->computeParticleCellKeys( params );
                            dynamic_cast<SpeciesVAdaptiveMixedSort *>( mypatch->vecSpecies[ispec] )->reconfigure_operators( params, mypatch );
                            // sorting will be necessary when clrw compatible
                            //mypatch->vecSpecies[ispec]->sortParticles(params);
                        }
                    } // end test patch_particle_created[ithread][j]
                } // end j loop
            } // End ithread loop
        }
        // Adaptive vectorization mode -- always sort
        else if( params.vectorization_mode == "adaptive" ) {
            //#pragma omp master
            //{
#ifndef _NO_MPI_TM
            #pragma omp for schedule(static) private(mypatch)
#endif
            for( int ithread=0; ithread < max_threads ; ithread++ ) {
                for( unsigned int j=0; j< ( patch_to_be_created[ithread] ).size(); j++ ) {
                
                    // Current newly created patch
                    mypatch = vecPatches.patches_[patch_to_be_created[ithread][j]];
                    
                    // If new particles are required
                    if( patch_particle_created[ithread][j] ) {
                        for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                            mypatch->vecSpecies[ispec]->computeParticleCellKeys( params );
                            mypatch->vecSpecies[ispec]->configuration( params, mypatch );
                            mypatch->vecSpecies[ispec]->sortParticles( params );
                            
                        }
                    }
                    // Patches that have received particles from another patch
                    // without the creation of new particles
                    else { // (patch_particle_created[ithread][j] == false)
                        for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                            // For the adaptive vectorization, we partially reconfigure the patch
                            // We do not have to sort, but operators may have to be reconfigured
                            mypatch->vecSpecies[ispec]->computeParticleCellKeys( params );
                            dynamic_cast<SpeciesVAdaptive *>( mypatch->vecSpecies[ispec] )->reconfigure_operators( params, mypatch );
                        }
                    } // end test patch_particle_created[ithread][j]
                } // end j loop
            } // End ithread loop
        }
        
        // Diagnostic Track Particles
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static) private(mypatch)
#endif
        for( int ithread=0; ithread < max_threads ; ithread++ ) {
            for( unsigned int j=0; j< ( patch_to_be_created[ithread] ).size(); j++ ) {
            
                // Current newly created patch
                mypatch = vecPatches.patches_[patch_to_be_created[ithread][j]];
                
                // If new particles are required
                if( patch_particle_created[ithread][j] ) {
                    // We define the IDs of the new particles
                    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ )
                        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
                            track->setIDs( mypatch );
                        }
                } // end test patch_particle_created[ithread][j]
            } // end j loop
        } // End ithread loop
#ifndef _NO_MPI_TM
        #pragma omp single nowait
#endif
        {
            x_moved += cell_length_x_*params.patch_size_[0];
            vecPatches.updateFieldList( smpi ) ;
            //update list fields for species diag too ??
            
            // Tell that the patches moved this iteration (needed for probes)
            vecPatches.lastIterationPatchesMoved = itime;
        }
        
        std::vector<double> poynting[2];
        poynting[0].resize( params.nDim_field, 0.0 );
        poynting[1].resize( params.nDim_field, 0.0 );
        
        double energy_field_out( 0. );
        double energy_field_inj( 0. );
        std::vector<double> energy_part_out( nSpecies, 0. );
        std::vector<double> energy_part_inj( nSpecies, 0. );
        std::vector<double> ukin_new( nSpecies, 0. );
        std::vector<double> ukin_bc ( nSpecies, 0. );
        std::vector<double> urad    ( nSpecies, 0. );
        
        // Delete useless patches while saving some scalar quantities
        for( unsigned int j=0; j < delete_patches_.size(); j++ ) {
            mypatch = delete_patches_[j];
            
            if( mypatch->isXmin() ) {
                for( unsigned int jp=0; jp<2; jp++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
                    for( unsigned int i=0 ; i<params.nDim_field ; i++ ) { //axis 0=x, 1=y, 2=z
                        poynting[jp][i] += mypatch->EMfields->poynting[jp][i];
                    }
                }
                
                energy_field_out += mypatch->EMfields->nrj_mw_out + mypatch->EMfields->computeEnergy();
                energy_field_inj += mypatch->EMfields->nrj_mw_inj;
                for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                    energy_part_out[ispec] += mypatch->vecSpecies[ispec]->nrj_mw_out + mypatch->vecSpecies[ispec]->computeEnergy();
                    energy_part_inj[ispec] += mypatch->vecSpecies[ispec]->nrj_mw_inj;
                    ukin_new[ispec] += mypatch->vecSpecies[ispec]->nrj_new_part_;
                    ukin_bc [ispec] += mypatch->vecSpecies[ispec]->nrj_bc_lost;
                    urad    [ispec] += mypatch->vecSpecies[ispec]->nrj_radiated_;
                }
            }
            
            delete  mypatch;
        }
        
        // Also account for new patches in scalars
        for( unsigned int j=0; j< patch_to_be_created[my_thread].size(); j++ ) {
            mypatch = vecPatches.patches_[patch_to_be_created[my_thread][j]];
            
            if( mypatch->isXmax() ) {
                energy_field_inj += mypatch->EMfields->computeEnergy();
                for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                    energy_part_inj[ispec] += mypatch->vecSpecies[ispec]->computeEnergy();
                }
            }
        }
        
        // Add the scalars to the patch #0
#ifndef _NO_MPI_TM
        #pragma omp critical
#endif
        {
            for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
                for( unsigned int i=0 ; i< params.nDim_field ; i++ ) { //axis 0=x, 1=y, 2=z
                    vecPatches( 0 )->EMfields->poynting[j][i] += poynting[j][i];
                }
            }
            
            vecPatches( 0 )->EMfields->nrj_mw_out += energy_field_out;
            vecPatches( 0 )->EMfields->nrj_mw_inj += energy_field_inj;
            for( unsigned int ispec=0 ; ispec<nSpecies ; ispec++ ) {
                vecPatches( 0 )->vecSpecies[ispec]->nrj_mw_out    += energy_part_out[ispec];
                vecPatches( 0 )->vecSpecies[ispec]->nrj_mw_inj    += energy_part_inj[ispec];
                vecPatches( 0 )->vecSpecies[ispec]->nrj_new_part_ += ukin_new[ispec];
                vecPatches( 0 )->vecSpecies[ispec]->nrj_bc_lost   += ukin_bc [ispec];
                vecPatches( 0 )->vecSpecies[ispec]->nrj_radiated_ += urad    [ispec];
            }
        }

#ifdef _NO_MPI_TM
    } // end omp master
#endif

    #pragma omp barrier
    #pragma omp master
    {
        // for( unsigned int i=0; i<vecPatches.size(); i++ ){
        //     MESSAGE(vecPatches(i)->vecSpecies[0]->getNrjOutMW());
        // }
        if (params.multiple_decomposition) {
            if ( params.geometry != "AMcylindrical" )
                operate(region, vecPatches, smpi, params, time_dual);
            else {
                operate(region, vecPatches, smpi, params, time_dual, params.nmodes);
            }
        }
    }
    #pragma omp barrier


    if (params.multiple_decomposition) {
        if ( params.geometry != "AMcylindrical" ) {
            // warkaround for !params.full_B_exchange (in 3D, with SM some border elements are not computed)
            SyncVectorPatch::exchangeE( params, region.vecPatch_, smpi );
            SyncVectorPatch::finalizeexchangeE( params, region.vecPatch_ );
            SyncVectorPatch::exchangeB( params, region.vecPatch_, smpi );
            SyncVectorPatch::finalizeexchangeB( params, region.vecPatch_ );
        }
        else {
            for (unsigned int imode = 0 ; imode < params.nmodes ; imode++  ) {
                SyncVectorPatch::exchangeE( params, region.vecPatch_, imode, smpi );
                SyncVectorPatch::exchangeB( params, region.vecPatch_, imode, smpi );
            }
        }
    }

}

void SimWindow::operate(Region& region,  VectorPatch&, SmileiMPI*, Params& params, double time_dual)
{
    region.patch_->exchangeField_movewin( region.patch_->EMfields->Ex_, params.patch_size_[0] );
    region.patch_->exchangeField_movewin( region.patch_->EMfields->Ey_, params.patch_size_[0] );
    region.patch_->exchangeField_movewin( region.patch_->EMfields->Ez_, params.patch_size_[0] );
    
    if (region.patch_->EMfields->Bx_->data_!= region.patch_->EMfields->Bx_m->data_) {
        region.patch_->exchangeField_movewin( region.patch_->EMfields->Bx_, params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region.patch_->EMfields->By_, params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region.patch_->EMfields->Bz_, params.patch_size_[0] );
    }
    
    region.patch_->exchangeField_movewin( region.patch_->EMfields->Bx_m, params.patch_size_[0] );
    region.patch_->exchangeField_movewin( region.patch_->EMfields->By_m, params.patch_size_[0] );
    region.patch_->exchangeField_movewin( region.patch_->EMfields->Bz_m, params.patch_size_[0] );

    if (params.is_spectral) {
        region.patch_->exchangeField_movewin( region.patch_->EMfields->rho_, params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region.patch_->EMfields->rhoold_, params.patch_size_[0] );
    }

    for( unsigned int bcId=2; bcId<2*params.nDim_field; bcId++ ){
        if( (dynamic_cast<ElectroMagnBC2D_PML *>( region.patch_->EMfields->emBoundCond[bcId] ) || dynamic_cast<ElectroMagnBC3D_PML *>( region.patch_->EMfields->emBoundCond[bcId] )) ){
            if( dynamic_cast<ElectroMagnBC2D_PML *>( region.patch_->EMfields->emBoundCond[bcId] )){
                ElectroMagnBC2D_PML *embc = static_cast<ElectroMagnBC2D_PML *>( region.patch_->EMfields->emBoundCond[bcId] );
                exchangePML_movewin( region, embc, params.patch_size_[0] );
            } else {
                ElectroMagnBC3D_PML *embc = static_cast<ElectroMagnBC3D_PML *>( region.patch_->EMfields->emBoundCond[bcId] );
                exchangePML_movewin( region, embc, params.patch_size_[0] );
            }
        }
    }

    //DoubleGrids::syncFieldsOnRegion( vecPatches, region, params, smpi );

    region.patch_->EMfields->laserDisabled();
    region.patch_->EMfields->emBoundCond[0]->apply(region.patch_->EMfields, time_dual, region.patch_);
    region.patch_->EMfields->emBoundCond[1]->apply(region.patch_->EMfields, time_dual, region.patch_);
    // External fields

    //mypatch->EMfields->emBoundCond[1]->disableExternalFields();

    // Deadlock if moving window & load balancing enabled
    //     Recompute patch distribution does not change
    //if (params.multiple_decomposition) {
    //    region.reset_mapping();
    //    region.identify_additional_patches( smpi, vecPatches, params );
    //    region.identify_missing_patches( smpi, vecPatches, params );
    //}


}


void SimWindow::operate(Region& region,  VectorPatch&, SmileiMPI*, Params& params, double time_dual, unsigned int nmodes)
{
    ElectroMagnAM * region_fields = static_cast<ElectroMagnAM *>( region.patch_->EMfields );
   
    for (unsigned int imode = 0; imode < nmodes; imode++){
        region.patch_->exchangeField_movewin( region_fields->El_[imode], params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region_fields->Er_[imode], params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region_fields->Et_[imode], params.patch_size_[0] );
        
        if (region_fields->Bl_[imode]->cdata_!= region_fields->Bl_m[imode]->cdata_) {
            region.patch_->exchangeField_movewin( region_fields->Bl_[imode], params.patch_size_[0] );
            region.patch_->exchangeField_movewin( region_fields->Br_[imode], params.patch_size_[0] );
            region.patch_->exchangeField_movewin( region_fields->Bt_[imode], params.patch_size_[0] );
        }

        region.patch_->exchangeField_movewin( region_fields->Bl_m[imode], params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region_fields->Br_m[imode], params.patch_size_[0] );
        region.patch_->exchangeField_movewin( region_fields->Bt_m[imode], params.patch_size_[0] );

        if (params.is_spectral) {
            region.patch_->exchangeField_movewin( region_fields->rho_AM_[imode], params.patch_size_[0] );
            region.patch_->exchangeField_movewin( region_fields->rho_old_AM_[imode], params.patch_size_[0] );
        }

        if( dynamic_cast<ElectroMagnBCAM_PML *>( region.patch_->EMfields->emBoundCond[3] )){
            ElectroMagnBCAM_PML *embc = static_cast<ElectroMagnBCAM_PML *>( region_fields->emBoundCond[3] );
            if (embc->Hl_[imode]) {
                region.patch_->exchangeField_movewin( embc->Hl_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Hr_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Ht_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Bl_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Br_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Bt_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->El_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Er_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Et_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Dl_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Dr_[imode], params.patch_size_[0] );
                region.patch_->exchangeField_movewin( embc->Dt_[imode], params.patch_size_[0] );
            }

        }
    }

    //DoubleGrids::syncFieldsOnRegion( vecPatches, region, params, smpi );

    region_fields->laserDisabled();
    region.patch_->EMfields->emBoundCond[0]->apply(region.patch_->EMfields, time_dual, region.patch_);
    region.patch_->EMfields->emBoundCond[1]->apply(region.patch_->EMfields, time_dual, region.patch_);
    // External fields

    //mypatch->EMfields->emBoundCond[1]->disableExternalFields();

    // Deadlock if moving window & load balancing enabled
    //     Recompute patch distribution does not change
    //if (params.multiple_decomposition) {
    //    region.reset_mapping();
    //    region.identify_additional_patches( smpi, vecPatches, params );
    //    region.identify_missing_patches( smpi, vecPatches, params );
    //}
}

template <typename Tpml>
void  SimWindow::exchangePML_movewin( Region& region, Tpml embc, int clrw ) {
                if (embc->Hx_) {
                    region.patch_->exchangeField_movewin( embc->Hx_, clrw );
                    region.patch_->exchangeField_movewin( embc->Hy_, clrw );
                    region.patch_->exchangeField_movewin( embc->Hz_, clrw );
                    region.patch_->exchangeField_movewin( embc->Bx_, clrw );
                    region.patch_->exchangeField_movewin( embc->By_, clrw );
                    region.patch_->exchangeField_movewin( embc->Bz_, clrw );
                    region.patch_->exchangeField_movewin( embc->Ex_, clrw );
                    region.patch_->exchangeField_movewin( embc->Ey_, clrw );
                    region.patch_->exchangeField_movewin( embc->Ez_, clrw );
                    region.patch_->exchangeField_movewin( embc->Dx_, clrw );
                    region.patch_->exchangeField_movewin( embc->Dy_, clrw );
                    region.patch_->exchangeField_movewin( embc->Dz_, clrw );
                }
}

