
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
                        vector<Diagnostic *> &localDiags )
{
    int ithread, tid( 0 );
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize) { // moving particle
    
        smpi->dynamics_resize( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );
        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[0] += MPI_Wtime() - timer;
#endif

            // Ionization
            if( Ionize ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                ( *Ionize )( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, patch, Proj );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
            }

            if( time_dual<=time_frozen_ ) continue; // Do not push frozen particles

            // Radiation losses
            if( Radiate ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Radiation process
                ( *Radiate )( *particles, photon_species_, smpi,
                              RadiationTables,
                              nrj_radiation,
                              particles->first_index[ibin],
                              particles->last_index[ibin], ithread );

                // Update scalar variable for diagnostics
                // nrj_radiation += Radiate->getRadiatedEnergy();

                // Update the quantum parameter chi
                // Radiate->computeParticlesChi( *particles,
                //                               smpi,
                //                               first_index[ibin],
                //                               last_index[ibin],
                //                               ithread );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[5] += MPI_Wtime() - timer;
#endif

            }


            // Multiphoton Breit-Wheeler
            if( Multiphoton_Breit_Wheeler_process ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Pair generation process
                ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                        smpi,
                                                        MultiphotonBreitWheelerTables,
                                                        particles->first_index[ibin], particles->last_index[ibin], ithread );

                // Update scalar variable for diagnostics
                // We reuse nrj_radiation for the pairs
                nrj_radiation += Multiphoton_Breit_Wheeler_process->getPairEnergy();

                // Update the photon quantum parameter chi of all photons
                Multiphoton_Breit_Wheeler_process->compute_thread_chiph( *particles,
                        smpi,
                        particles->first_index[ibin],
                        particles->last_index[ibin],
                        ithread );
                
                // Suppression of the decayed photons into pairs
                Multiphoton_Breit_Wheeler_process->decayed_photon_cleaning(
                    *particles, smpi, ibin, particles->first_index.size(), &particles->first_index[0], &particles->last_index[0], ithread );
                    
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[6] += MPI_Wtime() - timer;
#endif

            }

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Push the particles and the photons
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
            //particles->testMove( particles->first_index[ibin], particles->last_index[ibin], params );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[1] += MPI_Wtime() - timer;
#endif

        } //ibin
        
        if( time_dual>time_frozen_){ // do not apply particles BC nor project frozen particles
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
                double ener_iPart( 0. );

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        (*partWalls)[iwall]->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                        nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                    }
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                    nrj_lost_per_thd[tid] += mass_ * ener_iPart;

                } else if( mass_==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        (*partWalls)[iwall]->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                        nrj_lost_per_thd[tid] += ener_iPart;
                    }

                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                    nrj_lost_per_thd[tid] += ener_iPart;

                }

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[3] += MPI_Wtime() - timer;
#endif

                //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Project currents if not a Test species and charges as well if a diag is needed.
                // Do not project if a photon
                if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                    Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
                }

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[2] += MPI_Wtime() - timer;
#endif

            }// ibin
        } // end if moving particle

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

//        // Add the ionized electrons to the electron species
//        if (Ionize)
//            electron_species->importParticles( params, patch, Ionize->new_electrons, localDiags );
//
//        // Radiation losses
//        if (Radiate)
//        {
//            // If creation of macro-photon, we add them to photon_species
//            if (photon_species)
//            {
//                photon_species->importParticles(params,
//                                                patch,
//                                                Radiate->new_photons_,
//                                                localDiags);
//            }
//        }
//
//        // Multiphoton Breit-Wheeler
//        if (Multiphoton_Breit_Wheeler_process)
//        {
//
//            // Addition of the electron-positron particles
//            for (int k=0; k<2; k++) {
//                mBW_pair_species[k]->importParticles(params,
//                                             patch,
//                                             Multiphoton_Breit_Wheeler_process->new_pair[k],
//                                             localDiags);
//            }
//        }

    } //End if moving or ionized particles

    if(time_dual <= time_frozen_ && diag_flag &&( !particles->is_test ) ) { //immobile particle (at the moment only project density)
        if( params.geometry != "AMcylindrical" ) {
            double *b_rho=nullptr;
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            }
        } else {
            int n_species = patch->vecSpecies.size();
            complex<double> *b_rho=nullptr;
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                int ifield = imode*n_species+ispec;
                b_rho = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                    for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }
        }
    } // End projection for frozen particles
}

