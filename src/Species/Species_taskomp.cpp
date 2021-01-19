
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
    int ithread, tid( 0 );
// #ifdef _OPENMP
//     ithread = omp_get_thread_num();
// #else
//     ithread = 0;
// #endif

    //ithread = 0;
    ithread = buffer_id;

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    //std::vector<double> nrj_lost_per_thd( 1, 0. );
    int Nbins = particles->first_index.size();
    std::vector<double> nrj_lost_per_bin( Nbins, 0. );

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle
    
        //smpi->dynamics_resize( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );
        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_pushed[ibin])
            {
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[0] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Push the particles and the photons
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
            //particles->testMove( particles->first_index[ibin], particles->last_index[ibin], params );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[1] += MPI_Wtime() - timer;
#endif
            } // end task
        } //ibin
      
        } // end first if moving particle 



        if( time_dual>time_frozen_){ // do not apply particles BC nor project frozen particles

            #pragma omp taskgroup
            {
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
                {
                double ener_iPart( 0. );

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        (*partWalls)[iwall]->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                        nrj_lost_per_bin[ibin] += mass_ * ener_iPart;
                    }
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                    nrj_lost_per_bin[ibin] += mass_ * ener_iPart;
                
                } else if( mass_==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        (*partWalls)[iwall]->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                        nrj_lost_per_bin[ibin] += ener_iPart;
                    }
                
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                    nrj_lost_per_bin[ibin] += ener_iPart;
                
                }

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[3] += MPI_Wtime() - timer;
#endif
                } // end task
            } // ibin

            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
                //START EXCHANGE PARTICLES OF THE CURRENT BIN ?
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin])
                {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Reset densities sub-buffers
                for (int i = 0; i < size_proj_buffer_Jx; i++)  b_Jx [ibin][i]  = 0.0;
                for (int i = 0; i < size_proj_buffer_Jy; i++)  b_Jy [ibin][i]  = 0.0;
                for (int i = 0; i < size_proj_buffer_Jz; i++)  b_Jz [ibin][i]  = 0.0;
                for (int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i] = 0.0;
                
                // Project currents if not a Test species and charges as well if a diag is needed.
                // Do not project if a photon
                if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                    if (params.geometry == "2Dcartesian"){
                        Projector2D2Order *Proj2D = static_cast<Projector2D2Order *>(Proj);
                        Proj2D->currentsAndDensityWrapperOnBuffers( b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin], ibin*clrw, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
                    } else if (params.geometry == "3Dcartesian"){
                        Projector3D2Order *Proj3D = static_cast<Projector3D2Order *>(Proj);
                        Proj3D->currentsAndDensityWrapperOnBuffers( b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin], ibin*clrw, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
                    } else {ERROR("Task strategy not yet implemented in 1Dcartesian or AMcylindrical geometries");}
                
                } // end condition on test and mass

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[2] += MPI_Wtime() - timer;
#endif
                }//end task
            }// ibin

        }// end taskgroup  

     } // end second if moving particle


     // #pragma omp task default(shared)
     // {
         for( unsigned int ibin=0 ; ibin<nrj_lost_per_bin.size() ; ibin++ ) {
             nrj_bc_lost += nrj_lost_per_bin[ibin];
         }
     // } // end task

     //smpi->reduce_dynamics_buffer_size( ithread, params.geometry=="AMcylindrical" );


//////// Projection for frozen particles

    // if(time_dual <= time_frozen_ && diag_flag &&( !particles->is_test ) ) { //immobile particle (at the moment only project density)
    //     if( params.geometry != "AMcylindrical" ) {
    //         double *b_rho=nullptr;
    //         for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
    //             b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
    //             for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
    //                 Proj->basic( b_rho, ( *particles ), iPart, 0 );
    //             }
    //         }
    //     } else {
    //         int n_species = patch->vecSpecies.size();
    //         complex<double> *b_rho=nullptr;
    //         ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
    //         for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
    //             int ifield = imode*n_species+ispec;
    //             b_rho = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
    //             for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
    //                 for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
    //                     Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
    //                 }
    //             }
    //         }
    //     }
    // } // End projection for frozen particles

} // end dynamicsWithTasks

