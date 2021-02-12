
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
#include "ElectroMagn2D.h"
#include "Projector2D2Order.h"
#include "ElectroMagn3D.h"
#include "Projector3D2Order.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species_taskomp
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species_taskomp::Species_taskomp( Params &params, Patch *patch )
    : Species( params, patch )
{

    // Init tags for the task dependencies of the particle operations
    bin_has_pushed            = new int[particles->first_index.size()];
    bin_has_done_particles_BC = new int[particles->first_index.size()];
    bin_has_projected         = new int[particles->first_index.size()];

    //! buffers for currents and charge
    b_Jx.resize(particles->first_index.size());
    b_Jy.resize(particles->first_index.size());
    b_Jz.resize(particles->first_index.size());
    b_rho.resize(particles->first_index.size());


    for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
        // allocate current-buffers, then put to zero their content
        b_Jx[ibin]  = new double[size_proj_buffer_Jx ];
        b_Jy[ibin]  = new double[size_proj_buffer_Jy ];
        b_Jz[ibin]  = new double[size_proj_buffer_Jz ];
        b_rho[ibin] = new double[size_proj_buffer_rho];
    }

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

    if (bin_has_pushed != NULL){
        delete bin_has_pushed;
        delete bin_has_done_particles_BC;
        delete bin_has_projected;
    }

    for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
        // delete buffers
        delete[] b_Jx[ibin];
        delete[] b_Jy[ibin];
        delete[] b_Jz[ibin];
        delete[] b_rho[ibin];
    }

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
    int tid( 0 );

    

#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif

    unsigned int iPart;

    int Nbins = particles->first_index.size();
    std::vector<double> nrj_lost_per_bin( Nbins, 0. );
    #pragma omp taskgroup
    {
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle
        // resize the dynamics buffers to treat all the particles in this Patch ipatch and Species ispec
        smpi->dynamics_resize( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_pushed[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_pushed[ibin])
#endif
            {

            // Reset densities sub-buffers - each of these buffers stores a grid density on the ibin physical space
            // This must be done before Projection and before Ionization (because of the ionization currents)
            for (int i = 0; i < size_proj_buffer_Jx; i++)  b_Jx[ibin][i]  = 0.0;
            for (int i = 0; i < size_proj_buffer_Jy; i++)  b_Jy[ibin][i]  = 0.0;
            for (int i = 0; i < size_proj_buffer_Jz; i++)  b_Jz[ibin][i]  = 0.0;
            for (int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i] = 0.0;
                
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), buffer_id );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[0*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif


            // Ionization
            if( Ionize ) {          
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
                vector<double> *Epart = &( smpi->dynamics_Epart[buffer_id] );
                Ionize->ionizationTunnelWithTasks( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, patch, Proj, ibin, ibin*clrw, b_Jx [ibin], b_Jy [ibin], b_Jz [ibin] );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[4*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end Ionize

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Push the particles and the photons
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], buffer_id );
            //particles->testMove( particles->first_index[ibin], particles->last_index[ibin], params );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[1*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task for Interp+Push on ibin
        } // end ibin loop for Interp+Push
      

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
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

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin])
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin])
#endif
            {
                
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif

            // // Reset densities sub-buffers - each of these buffers stores a grid density on the ibin physical space
            // for (int i = 0; i < size_proj_buffer_Jx; i++)  b_Jx [ibin][i]  = 0.0;
            // for (int i = 0; i < size_proj_buffer_Jy; i++)  b_Jy [ibin][i]  = 0.0;
            // for (int i = 0; i < size_proj_buffer_Jz; i++)  b_Jz [ibin][i]  = 0.0;
            // for (int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i] = 0.0;
                
            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                if (params.geometry == "2Dcartesian"){
                    Projector2D2Order *Proj2D = static_cast<Projector2D2Order *>(Proj);
                    Proj2D->currentsAndDensityWrapperOnBuffers( b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin], ibin*clrw, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], buffer_id, diag_flag, params.is_spectral, ispec );
                } else if (params.geometry == "3Dcartesian"){
                    Projector3D2Order *Proj3D = static_cast<Projector3D2Order *>(Proj);
                    Proj3D->currentsAndDensityWrapperOnBuffers( b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin], ibin*clrw, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], buffer_id, diag_flag, params.is_spectral, ispec );
                } else {ERROR("Task strategy not yet implemented in 1Dcartesian or AMcylindrical geometries");}
            } // end condition on test and mass

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[2*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            }//end task for Proj of ibin
         }// end ibin loop for Proj

     } // end if moving particle
     }// end taskgroup for all the Interp, Push, Particles BC and Projector tasks
  
     // reduction of the lost energy in each ibin 
     // the taskgroup before ensures that it is done after the particles BC
#ifdef  __DETAILED_TIMERS
     #pragma omp task default(shared) private(ithread,timer) 
#else
     #pragma omp task default(shared)
#endif
     {
     
     // reduce the energy lost with BC per bin
     for( unsigned int ibin=0 ; ibin<nrj_lost_per_bin.size() ; ibin++ ) {
        nrj_bc_lost += nrj_lost_per_bin[ibin];
     }

     // unite the electrons created by ionization in a unique list for the patch 
     // The taskgroup above ensures that this is done after the ionization method
     if( Ionize ) {          
#ifdef  __DETAILED_TIMERS
         timer = MPI_Wtime();
         ithread = omp_get_thread_num();
#endif
         Ionize->joinNewElectrons(particles->first_index.size());

#ifdef  __DETAILED_TIMERS
         patch->patch_timers_[4*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
     } // end Ionize     

     
     } // end task for lost energy reduction

     // smpi->reduce_dynamics_buffer_size( buffer_id, params.geometry=="AMcylindrical" );

} // end dynamicsWithTasks
