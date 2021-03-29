
#include "Species_taskomp.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

#include <omp.h>

// IDRIS
#include <cstring>
// IDRIS
#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "RadiationFactory.h"
#include "MultiphotonBreitWheelerFactory.h"
#include "MergingFactory.h"
#include "PartBoundCond.h"
#include "PartWall.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "Profile.h"
#include "ElectroMagnAM.h"
#include "Projector.h"
#include "ProjectorFactory.h"
#include "ParticleCreator.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

#include "DiagnosticTrack.h"

// necessary for the static_cast
#include "ProjectorAM2Order.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species_taskomp
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species_taskomp::Species_taskomp( Params &params, Patch *patch )
    : Species( params, patch )
{
    Nbins = particles->first_index.size();
    // // Init tags for the task dependencies of the particle operations
    // bin_has_interpolated                   = new int[Nbins+1]; // the last element is used to manage the Multiphoton Breit Wheeler dependency
    // bin_has_ionized                        = new int[Nbins];
    // bin_has_radiated                       = new int[Nbins];
    // bin_has_done_Multiphoton_Breit_Wheeler = new int[Nbins];
    // bin_has_pushed                         = new int[Nbins];
    // bin_has_done_particles_BC              = new int[Nbins];
    // bin_has_projected                      = new int[Nbins];

    nrj_lost_per_bin                       = new double[Nbins];
    nrj_radiation_per_bin                  = new double[Nbins];

    // //! buffers for currents and charge
    // b_Jx.resize(Nbins);
    // b_Jy.resize(Nbins);
    // b_Jz.resize(Nbins);
    // b_rho.resize(Nbins);
    // 
    // for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
    //     // allocate current-buffers, then put to zero their content
    //     b_Jx[ibin]  = new double[size_proj_buffer_Jx ];
    //     b_Jy[ibin]  = new double[size_proj_buffer_Jy ];
    //     b_Jz[ibin]  = new double[size_proj_buffer_Jz ];
    //     b_rho[ibin] = new double[size_proj_buffer_rho];
    // }

}//END Species creator



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species_taskomp::~Species_taskomp()
{
    // delete Push;
    // delete Interp;
    // delete Proj;
    //
    //
    // if( partBoundCond ) {
    //     delete partBoundCond;
    // }
    // if( particles_per_cell_profile_ ) {
    //     delete particles_per_cell_profile_;
    // }
    // if( charge_profile_ ) {
    //     delete charge_profile_;
    // }
    // if( density_profile_ ) {
    //     delete density_profile_;
    // }
    // for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
    //     delete velocity_profile_[i];
    // }
    // for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
    //     delete temperature_profile_[i];
    // }
    // if( ionization_rate_!=Py_None ) {
    //     Py_DECREF( ionization_rate_ );
    // }

    // if (bin_has_interpolated != NULL){
    //     delete bin_has_interpolated;
    // }
    // if (bin_has_ionized != NULL){
    //     delete bin_has_ionized;
    // }
    // if (bin_has_radiated != NULL){
    //     delete bin_has_radiated;
    // }
    // if (bin_has_done_Multiphoton_Breit_Wheeler != NULL){
    //     delete bin_has_done_Multiphoton_Breit_Wheeler;
    // }
    // if (bin_has_pushed != NULL){
    //     delete bin_has_pushed;
    // }
    // if (bin_has_done_particles_BC != NULL){
    //     delete bin_has_done_particles_BC;
    // }
    // if (bin_has_projected != NULL){
    //     delete bin_has_projected;
    // }
    if (nrj_lost_per_bin != NULL){
        delete nrj_lost_per_bin;
    }
    if (nrj_radiation_per_bin != NULL){
        delete nrj_radiation_per_bin;
    }

    // for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
    //     // delete buffers
    //     delete[] b_Jx[ibin];
    //     delete[] b_Jy[ibin];
    //     delete[] b_Jz[ibin];
    //     delete[] b_rho[ibin];
    // }

}

void Species_taskomp::dynamics( double time_dual, unsigned int ispec,
                        ElectroMagn *EMfields,
                        Params &params, bool diag_flag,
                        PartWalls *partWalls,
                        Patch *patch, SmileiMPI *smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                        vector<Diagnostic *> &localDiags )
{
}

void Species_taskomp::dynamicsWithTasks( double time_dual, unsigned int ispec,
                        ElectroMagn *EMfields,
                        Params &params, bool diag_flag,
                        PartWalls *partWalls,
                        Patch *patch, SmileiMPI *smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                        vector<Diagnostic *> &localDiags, int buffer_id )
{
    
#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif
    int bin_size0 = b_dim[0];
    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
        nrj_lost_per_bin[ibin] = 0.;
        nrj_radiation_per_bin[ibin] = 0.;
    }
    // Init tags for the task dependencies of the particle operations
    int *bin_has_interpolated                   = new int[Nbins+1]; // the last element is used to manage the Multiphoton Breit Wheeler dependency
    int *bin_has_ionized                        = new int[Nbins];
    int *bin_has_radiated                       = new int[Nbins];
    int *bin_has_done_Multiphoton_Breit_Wheeler = new int[Nbins];
    int *bin_has_pushed                         = new int[Nbins];
    int *bin_has_done_particles_BC              = new int[Nbins];
    int *bin_has_projected                      = new int[Nbins];

    int *bin_can_radiate = new int[Nbins];;
    int *bin_can_push = new int[Nbins];;

    if (Radiate){ // if Radiation True ... 
        if (!Ionize) { 
            // ... radiate only after ionization if present ...
            bin_can_radiate = bin_has_interpolated;
        } else { 
            // ... radiate directly after interpolation if ionization is not present ...
            bin_can_radiate = bin_has_ionized;
        }
        // ... and push only after radiation 
        bin_can_push = bin_has_radiated;
    } else { // if Radiation False ...
        if (Ionize){ 
            // ... push after ionization if present
            bin_can_push = bin_has_ionized;
        } else { 
            // ... push directly after interpolation if ionization is not present
            // A Species with mass = 0 cannot Ionize or Radiate, thus this this is the used dependency array.
            // Remember that the element ibin = Nbins of bin_has_interpolated 
            // is used to manage the pusher dependency on the photon cleaning
            bin_can_push = bin_has_interpolated;       
        }
    }

    #pragma omp taskgroup
    {
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_  || Ionize ) { // if moving particle or it can be ionized
        // resize the dynamics buffers to treat all the particles in this Patch ipatch and Species ispec
        smpi->dynamics_resize( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
            {

            if ( params.geometry != "AMcylindrical" ){
                // Reset densities sub-buffers - each of these buffers stores a grid density on the ibin physical space
                // This must be done before Projection and before Ionization (because of the ionization currents)
                for (unsigned int i = 0; i < size_proj_buffer_Jx; i++)  b_Jx[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jy; i++)  b_Jy[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jz; i++)  b_Jz[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
            } else {
                // Reset densities sub-buffers - each of these buffers stores a grid density on the ibin physical space
                // This must be done before Projection and before Ionization (because of the ionization currents)
                for (unsigned int i = 0; i < size_proj_buffer_Jl; i++)  b_Jl[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jr; i++)  b_Jr[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jt; i++)  b_Jt[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
            }
                
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), buffer_id );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[0*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } //end task Interpolator
        } // end ibin loop for Interpolator

        // Ionization
        if( Ionize ) {        
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_ionized[ibin]) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_ionized[ibin])
#endif
                {
    
#ifdef  __DETAILED_TIMERS
                ithread = omp_get_thread_num();
                timer = MPI_Wtime();
#endif
                vector<double> *Epart = &( smpi->dynamics_Epart[buffer_id] );
                Ionize->ionizationTunnelWithTasks( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, patch, Proj, ibin, ibin*clrw, b_Jx [ibin], b_Jy [ibin], b_Jz [ibin] );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            
                } // end task Ionize bin
            } // end ibin loop for Ionize
        } // end Ionize

            if( time_dual>time_frozen_ ){ // if moving particle push

                // Radiation losses
                if( Radiate ) {
                    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                        #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_can_radiate[ibin]) depend(out:bin_has_radiated[ibin]) private(ithread,timer)
#else
                        #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_can_radiate[ibin]) depend(out:bin_has_radiated[ibin])
#endif
                        {
                    

#ifdef  __DETAILED_TIMERS
                        ithread = omp_get_thread_num();
                        timer = MPI_Wtime();
#endif

                        // Radiation process
                        ( *Radiate )( *particles, photon_species_, smpi,
                                      RadiationTables,
                                      nrj_radiation_per_bin[ibin],
                                      particles->first_index[ibin],
                                      particles->last_index[ibin], buffer_id, ibin );

                        // Update scalar variable for diagnostics
                        // nrj_radiation += Radiate->getRadiatedEnergy();

                        // Update the quantum parameter chi
                        // Radiate->computeParticlesChi( *particles,
                        //                               smpi,
                        //                               first_index[ibin],
                        //                               last_index[ibin],
                        //                               ithread );

#ifdef  __DETAILED_TIMERS
                        patch->patch_timers_[5*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif

                    
                        } // end task Radiate bin
                    } // end ibin loop for Radiate
                } // end if Radiate

                if( Multiphoton_Breit_Wheeler_process ) {
                    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                        #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_done_Multiphoton_Breit_Wheeler[ibin]) private(ithread,timer)
#else
                        #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_done_Multiphoton_Breit_Wheeler[ibin])
#endif
                        {

#ifdef  __DETAILED_TIMERS
                        ithread = omp_get_thread_num();
                        timer = MPI_Wtime();
#endif
                        // Pair generation process
                        ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                        smpi,
                                                        MultiphotonBreitWheelerTables,
                                                        particles->first_index[ibin], particles->last_index[ibin], 
                                                        buffer_id, ibin );

                        // // Update scalar variable for diagnostics
                        // // We reuse nrj_radiation for the pairs
                        nrj_radiation_per_bin[ibin] += Multiphoton_Breit_Wheeler_process->getPairEnergyOfBin(ibin);

                        // Update the photon quantum parameter chi of all photons
                        Multiphoton_Breit_Wheeler_process->compute_thread_chiph( *particles,
                                smpi,
                                particles->first_index[ibin],
                                particles->last_index[ibin],
                                buffer_id );    

#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[6*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                        } // end Multiphoton Breit Wheeler on ibin
                    } // end ibin task for Multiphoton Breit Wheeler
                    #pragma omp taskwait
#ifdef  __DETAILED_TIMERS
                    #pragma omp task default(shared) depend(in:bin_has_done_Multiphoton_Breit_Wheeler[0:(Nbins-1)]) private(ithread,timer) depend(out:bin_has_interpolated[Nbins]) 
#else
                    #pragma omp task default(shared) depend(in:bin_has_done_Multiphoton_Breit_Wheeler[0:(Nbins-1)]) depend(out:bin_has_interpolated[Nbins]) 
#endif
                    {
#ifdef  __DETAILED_TIMERS
                    ithread = omp_get_thread_num();
                    timer = MPI_Wtime();
#endif
                    // clean decayed photons from arrays 
                    // this loop must not be parallelized unless race conditions are prevented
                    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                        Multiphoton_Breit_Wheeler_process->decayed_photon_cleaning(
                            *particles, smpi, ibin, particles->first_index.size(), &particles->first_index[0], &particles->last_index[0], buffer_id );              
                    } // end ibin loop to clean decayed photons
#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[6*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                    } // end task for photon cleaning for all bins
                } else {
                  // empty task for the pusher dependency
                    #pragma omp task default(shared) depend(out:bin_has_interpolated[Nbins])
                    {
                    // Remember that bin_has_interpolated[Nbins] 
                    // is used to manage the pusher dependency on the photon cleaning 
                    } 
                }// end if Multiphoton_Breit_Wheeler_process

                for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                    #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_can_push[ibin],bin_can_push[Nbins]) depend(out:bin_has_pushed[ibin]) private(ithread,timer)
#else
                    #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_can_push[ibin],bin_can_push[Nbins]) depend(out:bin_has_pushed[ibin])
#endif
                    {
#ifdef  __DETAILED_TIMERS
                    ithread = omp_get_thread_num();
                    timer = MPI_Wtime();
#endif

                    // Push the particles and the photons
                    ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], buffer_id );
                    //particles->testMove( particles->first_index[ibin], particles->last_index[ibin], params );

#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[1*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                    } // end task for Push on ibin
                } // end ibin loop for Push
        } // end if moving particle, radiate and push
    } // end if moving particle or it can be ionized 

    if( time_dual>time_frozen_){ // do not apply particles BC nor projection on frozen particles     
        
        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
#endif

            {
            double ener_iPart( 0. );

#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif

            // Apply wall and boundary conditions
            if( mass_>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    (*partWalls)[iwall]->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, buffer_id, ener_iPart );
                    nrj_lost_per_bin[ibin] += mass_ * ener_iPart;
                }
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, buffer_id, ener_iPart );
                nrj_lost_per_bin[ibin] += mass_ * ener_iPart;
                
            } else if( mass_==0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    (*partWalls)[iwall]->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, buffer_id, ener_iPart );
                    nrj_lost_per_bin[ibin] += ener_iPart;
                }
                
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, buffer_id, ener_iPart );
                nrj_lost_per_bin[ibin] += ener_iPart;
                
            }

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task for particles BC on ibin
        } // end ibin loop for particles BC

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin])
#else
            #pragma omp task default(shared) firstprivate(ibin,bin_size0) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin])
#endif
            {
                
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif
                
            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                if (params.geometry != "AMcylindrical"){
                    Proj->currentsAndDensityWrapperOnBuffers( b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin], 
                                                              ibin*clrw, *particles, smpi, 
                                                              particles->first_index[ibin], particles->last_index[ibin], 
                                                              buffer_id, diag_flag, params.is_spectral, ispec );
                } else {
                    ProjectorAM2Order *ProjAM = static_cast<ProjectorAM2Order *>(Proj);
                    ProjAM->currentsAndDensityWrapperOnAMBuffers( EMfields, b_Jl[ibin], b_Jr[ibin], b_Jt[ibin], b_rhoAM[ibin], 
                                                                ibin*clrw, bin_size0, *particles, smpi, 
                                                                particles->first_index[ibin], particles->last_index[ibin], 
                                                                buffer_id, diag_flag);
                }
            } // end condition on test and mass

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[2*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            }//end task for Proj of ibin
         }// end ibin loop for Proj

        // reduction of the lost energy in each ibin 
        // the taskgroup before ensures that it is done after the particles BC
#ifdef  __DETAILED_TIMERS
        #pragma omp task default(shared) private(ithread,timer) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
#else
        #pragma omp task default(shared) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
#endif
        {
        // reduce the energy lost with BC per bin
        for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
           nrj_bc_lost += nrj_lost_per_bin[ibin];
        }

        // sum the radiated energy
        // The taskgroup above ensures that this is done after the radiation method
        if( Radiate ) {
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
            ithread = omp_get_thread_num();
#endif

            for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
               nrj_radiation += nrj_radiation_per_bin[ibin];
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[5*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
        } // end if Radiate
        } // end task for lost/radiated energy reduction          

     } // end if moving particle

     }// end taskgroup for all the Interp, Push, Particles BC and Projector tasks
  
//      // reduction of the lost energy in each ibin 
//      // the taskgroup before ensures that it is done after the particles BC
// #ifdef  __DETAILED_TIMERS
//      #pragma omp task default(shared) private(ithread,timer) 
// #else
//      #pragma omp task default(shared)
// #endif
//      {
//      // reduce the energy lost with BC per bin
//      for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
//         nrj_bc_lost += nrj_lost_per_bin[ibin];
//      }
// 
//      // sum the radiated energy
//      // The taskgroup above ensures that this is done after the radiation method
//      if( Radiate ) {
// #ifdef  __DETAILED_TIMERS
//          timer = MPI_Wtime();
//          ithread = omp_get_thread_num();
// #endif
// 
//          for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
//             nrj_radiation += nrj_radiation_per_bin[ibin];
//          }
// #ifdef  __DETAILED_TIMERS
//          patch->patch_timers_[5*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
// #endif
//      } // end if Radiate
//      } // end task for lost/radiated energy reduction
// 
//      // smpi->reduce_dynamics_buffer_size( buffer_id, params.geometry=="AMcylindrical" );

} // end dynamicsWithTasks
