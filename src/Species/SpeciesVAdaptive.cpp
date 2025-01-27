#include "SpeciesVAdaptive.h"

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
#include "PartBoundCond.h"
#include "PartWall.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "Profile.h"

#include "Projector.h"
#include "ProjectorFactory.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

#include "DiagnosticTrack.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
SpeciesVAdaptive::SpeciesVAdaptive( Params &params, Patch *patch ) :
    SpeciesV( params, patch )
{
    initCluster( params, patch );
    npack_ = 0 ;
    packsize_ = 0;
}//END SpeciesVAdaptive creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
SpeciesVAdaptive::~SpeciesVAdaptive()
{
}

//! Method calculating the Particle dynamics (interpolation, pusher, projection)
//! without vectorized operators but with the cell sorting algorithm
void SpeciesVAdaptive::scalarDynamics( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields, Params &params, bool diag_flag,
        PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi,
        RadiationTables &RadiationTables,
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables )
{

    const int ithread = Tools::getOMPThreadNum();

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    unsigned int iPart;

    int tid( 0 );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize ) {
        // moving particle

        smpi->resizeBuffers( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );

        //Prepare for sorting
        for( unsigned int i=0; i<count.size(); i++ ) {
            count[i] = 0;
        }

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif


        smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 0);
        // Interpolate the fields at the particle position
        Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread, particles->first_index[0] );
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 0);

#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[0] += MPI_Wtime() - timer;
#endif

        // Ionization
        if( Ionize ) {
            // Interpolate the fields at the particle position
            //for (unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++)
            //    (*Interp)(EMfields, *particles, smpi, &(particles->first_index[scell]), &(particles->last_index[scell]), ithread );

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 5);
            for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {


#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                ( *Ionize )( particles, particles->first_index[scell], particles->last_index[scell], Epart, patch, Proj );
#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4] += MPI_Wtime() - timer;
#endif

            }
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 5);

        }
            // if( time_dual<=time_frozen_ ) continue; // Do not push nor project frozen particles

    if( time_dual>time_frozen_ ) { // do not radiate, compute multiphoton MBW, push, nor apply particles BC, nor project frozen particles

            // Radiation losses
            if( Radiate ) {

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 6);
                for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {

#ifdef  __DETAILED_TIMERS
                    timer = MPI_Wtime();
#endif
                // Radiation process
                ( *Radiate )( *particles,
                              radiated_photons_,
                              smpi,
                              RadiationTables, nrj_radiated_,
                              particles->first_index[scell], particles->last_index[scell], ithread );

                // // Update scalar variable for diagnostics
                // nrj_radiated_ += Radiate->getRadiatedEnergy();
                //
                // // Update the quantum parameter chi
                // Radiate->computeParticlesChi( *particles,
                //                               smpi,
                //                               first_index[scell],
                //                               last_index[scell],
                //                               ithread );
#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[5] += MPI_Wtime() - timer;
#endif

            }
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 6);
        } // end if Radiate

        // Multiphoton Breit-Wheeler
        if( Multiphoton_Breit_Wheeler_process ) {

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 7);
            for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {


#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                // Pair generation process
                // We reuse nrj_radiated_ for the pairs
                ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                        smpi,
                                                        mBW_pair_particles_,
                                                        mBW_pair_species_,
                                                        MultiphotonBreitWheelerTables,
                                                        nrj_radiated_,
                                                        particles->first_index[scell], particles->last_index[scell], ithread );

                // Update the photon quantum parameter chi of all photons
                Multiphoton_Breit_Wheeler_process->computeThreadPhotonChi( *particles,
                        smpi,
                        particles->first_index[scell],
                        particles->last_index[scell],
                        ithread );

                // Suppression of the decayed photons into pairs
                Multiphoton_Breit_Wheeler_process->removeDecayedPhotons(
                    *particles, smpi, scell, particles->first_index.size(), &particles->first_index[0], &particles->last_index[0], ithread );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[6] += MPI_Wtime() - timer;
#endif
            }
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 7);

        } // end if Multiphoton Breit Wheeler

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            size_t start = 0, stop = particles->last_index.back(), n = stop - start;
            vector<vector<double>> pold;
            particles->prepareInterpolatedFields( pold, start, n );
            
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 1);
            // Push the particles and the photons
            ( *Push )( *particles, smpi, 0, particles->last_index.back(), ithread, 0. );
            
            // Copy interpolated fields to persistent buffers if requested
            if( particles->interpolated_fields_ ) {
                particles->copyInterpolatedFields( &( smpi->dynamics_Epart[ithread][start] ), &( smpi->dynamics_Bpart[ithread][start] ), pold, start, n, n, mass_ );
            }
            
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 1);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[1] += MPI_Wtime() - timer;
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 2);
            for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {
                double energy_lost( 0. );
                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        (*partWalls )[iwall]->apply( this, particles->first_index[scell], particles->last_index[scell], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                        nrj_lost_per_thd[tid] += mass_ * energy_lost;
                    }

                    partBoundCond->apply( this, particles->first_index[scell], particles->last_index[scell], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                    nrj_lost_per_thd[tid] += mass_ * energy_lost;

                } else if( mass_==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        ( *partWalls )[iwall]->apply( this, particles->first_index[scell], particles->last_index[scell], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                        nrj_lost_per_thd[tid] += energy_lost;
                    }

                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    partBoundCond->apply( this, particles->first_index[scell], particles->last_index[scell], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                    nrj_lost_per_thd[tid] += energy_lost;

                } // end if mass_ > 0
            } // end loop on cells
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 2);


            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 0, 11);
            // for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {
            //
            //     // Apply wall and boundary conditions
            //     if( mass_>0 ) {
            //
            //         for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
            //             if ( particles->cell_keys[iPart] >= 0 ) {
            //                 //Compute cell_keys of remaining particles
            //                 for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
            //                     particles->cell_keys[iPart] *= this->length_[i];
            //                     particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
            //                 }
            //                 //First reduction of the count sort algorithm. Lost particles are not included.
            //                 count[particles->cell_keys[iPart]] ++;
            //             }
            //         }
            //
            //     } else if( mass_==0 ) {
            //
            //         for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
            //             if ( particles->cell_keys[iPart] >= 0 ) {
            //                  //Compute cell_keys of remaining particles
            //                 for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
            //                     particles->cell_keys[iPart] *= this->length_[i];
            //                     particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
            //                 }
            //                 //First reduction of the count sort algorithm. Lost particles are not included.
            //                 count[particles->cell_keys[iPart]] ++;
            //             }
            //         }
            //     } // end if mass_ > 0
            // } // end loop on cells

            // Cell keys
            computeParticleCellKeys( params );

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread, 1, 11);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3] += MPI_Wtime() - timer;
#endif

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,3);
                Proj->currentsAndDensityWrapper(
                    EMfields, *particles, smpi, particles->first_index[0],
                    particles->last_index.back(),
                    ithread,
                    diag_flag,
                    params.is_spectral,
                    ispec
            );
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,3);

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[2] += MPI_Wtime() - timer;
#endif

            }
        } // end if moving particle

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

    }  // end if moving particle or ionize

    if (time_dual <= time_frozen_ && diag_flag &&( !particles->is_test ) ) { //immobile particle (at the moment only project density)


        if( params.geometry != "AMcylindrical" ) {

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,3);
            double *b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
            for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell ++ ) { //Loop for projection on buffer_proj
                for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            }
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,3);

        } else { // AM case
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            int n_species = patch->vecSpecies.size();

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,3);
            for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                int ifield = imode*n_species+ispec;
                complex<double> *b_rho = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell ++ ) { //Loop for projection on buffer_proj
                    for( int iPart=particles->first_index[scell] ; iPart<particles->last_index[scell]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,3);
        }
    }

}//END scalarDynamics

#ifdef _OMPTASKS
void SpeciesVAdaptive::scalarDynamicsTasks( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields, Params &params, bool diag_flag,
        PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi,
        RadiationTables &RadiationTables,
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables, int buffer_id )
{

#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif

    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
        nrj_lost_per_bin[ibin] = 0.;
        radiated_energy_per_bin[ibin] = 0.;
    }

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    int bin_size0 = b_dim[0]; // used for AM
    // Init tags for the task dependencies of the particle operations
    int *bin_has_interpolated                   = new int[Nbins];
    int *bin_has_pushed                         = new int[Nbins];
    int *bin_has_done_particles_BC              = new int[Nbins];


    if( npack_==0 ) {
        npack_    = 1;
        packsize_ = ( f_dim1-2*oversize[1] );

        //if ( ( (long int)particles->last_index.back() < (long int)60000 ) || (Radiate) || (Ionize) || (Multiphoton_Breit_Wheeler_process) )
        packsize_ *= ( f_dim0-2*oversize[0] );
        //else
        //    npack_ *= (f_dim0-2*oversize[0]);

        if( nDim_field == 3 ) {
            packsize_ *= ( f_dim2-2*oversize[2] );
        }
    }

    //unsigned int iPart;

    std::vector<double> nrj_lost_per_thd( 1, 0. );
    int ipack = 0;
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------

    // why was this used if the array is resized later?
    smpi->resizeBuffers( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

    //Prepare for sorting
    for( unsigned int i=0; i<count.size(); i++ ) {
        count[i] = 0;
    }

    int nparts_in_pack = particles->last_index[( ipack+1 ) * packsize_-1 ];
    smpi->resizeBuffers( buffer_id, nDim_field, nparts_in_pack );

    #pragma omp taskgroup
    {
    if( time_dual>time_frozen_ || Ionize ) { // moving particle



        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ){
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
            {

            if ( params.geometry != "AMcylindrical" ){
                // Reset densities sub-buffers - each of these buffers store a grid density on the ibin physical space
                // This must be done before Projection and before Ionization (because of the ionization currents)
                for (unsigned int i = 0; i < size_proj_buffer_Jx; i++)  b_Jx[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jy; i++)  b_Jy[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jz; i++)  b_Jz[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
            } else {
                // Reset densities sub-buffers - each of these buffers store a grid density on the ibin physical space
                // This must be done before Projection and before Ionization (because of the ionization currents)
                for (unsigned int i = 0; i < size_proj_buffer_Jl; i++)  b_Jl[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jr; i++)  b_Jr[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jt; i++)  b_Jt[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
            }

#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif


            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,0);
            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[first_cell_of_bin[ibin]] ),
                                   &( particles->last_index[last_cell_of_bin[ibin]] ),
                                   buffer_id, particles->first_index[0] );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[0*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } //end task Interpolator

        } // end ibin loop for Interpolator

        // Ionization
        if( Ionize ) {
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
                {

#ifdef  __DETAILED_TIMERS
                ithread = Tools::getOMPThreadNum();
                timer = MPI_Wtime();
#endif

                double *bJx, *bJy, *bJz;
                if (params.geometry != "AMcylindrical"){
                    bJx         = b_Jx[ibin];
                    bJy         = b_Jy[ibin];
                    bJz         = b_Jz[ibin];
                } else {
                    bJx         = NULL;
                    bJy         = NULL;
                    bJz         = NULL;
                }

                vector<double> *Epart = &( smpi->dynamics_Epart[buffer_id] );

                // Loop over scell is not performed since ionization operator is not vectorized
                // Instead, it is applied to all particles in the cells pertaining to ibin
                // for( unsigned int scell = first_cell_of_bin[ibin] ; scell <= last_cell_of_bin[ibin] ; scell++ ){
                //     Ionize->ionizationTunnelWithTasks( particles, particles->first_index[scell], particles->last_index[scell], Epart, patch, Proj, ibin, ibin*cluster_width_, bJx, bJy, bJz );
                // } // end cell loop for Interpolator

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,5);
                Ionize->ionizationTunnelWithTasks( particles, particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]], Epart, patch, Proj, ibin, ibin*cluster_width_, bJx, bJy, bJz );
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,5);


#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif

                } // end task Ionize bin
            } // end ibin loop for Ionize
        } // end Ionize

     } // end if moving particle or Ionization

     if( time_dual>time_frozen_ ) {

         // Radiation losses
         if( Radiate ) {
             for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                 #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
                 #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
                 {


#ifdef  __DETAILED_TIMERS
                 ithread = Tools::getOMPThreadNum();
                 timer = MPI_Wtime();
#endif

                 // Loop over scell is not performed since radiation operator is not completely vectorized
                 // Instead, it is applied to all particles in the cells pertaining to ibin
                 // for( unsigned int scell = first_cell_of_bin[ibin] ; scell <= last_cell_of_bin[ibin] ; scell++ ){
                 //     // Radiation process
                 //     ( *Radiate )( *particles, photon_species_, smpi,
                 //                   RadiationTables,
                 //                   radiated_energy_per_bin[ibin],
                 //                   particles->first_index[scell],
                 //                   particles->last_index[scell], buffer_id, ibin );
                 // }


                 smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,6);
                 // Radiation process
                 ( *Radiate )( *particles,
                               radiated_photons_,
                               smpi,
                               RadiationTables,
                               radiated_energy_per_bin[ibin],
                               particles->first_index[first_cell_of_bin[ibin]],
                               particles->last_index[last_cell_of_bin[ibin]], buffer_id, ibin );
                 smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,6);

                 // Update scalar variable for diagnostics
                 // radiated_energy += Radiate->getRadiatedEnergy();

                 // Update the quantum parameter chi
                 // Radiate->computeParticlesChi( *particles,
                 //                               smpi,
                 //                               first_index[ibin],
                 //                               last_index[ibin],
                 //                               ithread );

#ifdef  __DETAILED_TIMERS
                 patch->patch_timers_[5*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif


                 } // end task Radiate bin
             } // end ibin loop for Radiate
        } // end if Radiate

        if( Multiphoton_Breit_Wheeler_process ) {
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
                {

#ifdef  __DETAILED_TIMERS
                ithread = Tools::getOMPThreadNum();
                timer = MPI_Wtime();
#endif

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,7);
                // for( unsigned int scell = first_cell_of_bin[ibin] ; scell <= last_cell_of_bin[ibin] ; scell++ ){
                // Pair generation process
                double radiated_energy_bin = 0;
                ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                        smpi,
                                                        mBW_pair_particles_,
                                                        mBW_pair_species_,
                                                        MultiphotonBreitWheelerTables,
                                                        radiated_energy_bin,
                                                        particles->first_index[first_cell_of_bin[ibin]],
                                                        particles->last_index[last_cell_of_bin[ibin]],
                                                        buffer_id, ibin );

                radiated_energy_per_bin[ibin] = radiated_energy_bin;
                // Update the photon quantum parameter chi of all photons
                Multiphoton_Breit_Wheeler_process->computeThreadPhotonChi( *particles,
                                                        smpi,
                                                        particles->first_index[first_cell_of_bin[ibin]],
                                                        particles->last_index[last_cell_of_bin[ibin]],
                                                        buffer_id );


                // } // end scell
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,7);

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[6*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
                } // end Multiphoton Breit Wheeler on ibin
            } // end ibin task for Multiphoton Breit Wheeler
            #pragma omp taskwait
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) depend(out:bin_has_interpolated[0:(Nbins-1)]) private(ithread,timer)
#else
            #pragma omp task default(shared) depend(out:bin_has_interpolated[0:(Nbins-1)])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,7);
            // clean decayed photons from arrays
            // this loop must not be parallelized unless race conditions are prevented
            for( int scell = first_cell_of_bin[0] ; scell <= last_cell_of_bin[Nbins-1] ; scell++ ){
                // Suppression of the decayed photons into pairs
                Multiphoton_Breit_Wheeler_process->removeDecayedPhotons(
                    *particles, smpi, scell, particles->first_index.size(), &particles->first_index[0], &particles->last_index[0], buffer_id  );
            } // end scell

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,7);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[6*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task for photon cleaning for all bins
        } // end if Multiphoton Breit Wheeler

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_pushed[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_pushed[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            size_t start = particles->first_index[first_cell_of_bin[ibin]], stop = particles->last_index[last_cell_of_bin[ibin]], n = stop - start;
            vector<vector<double>> pold;
            particles->prepareInterpolatedFields( pold, start, n );
            
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);
            // Push the particles and the photons
            ( *Push )( *particles, smpi, start, stop, buffer_id, particles->first_index[0] );
            
            // Copy interpolated fields to persistent buffers if requested
            if( particles->interpolated_fields_ ) {
                particles->copyInterpolatedFields( &( smpi->dynamics_Epart[buffer_id][start] ), &( smpi->dynamics_Bpart[buffer_id][start] ), pold, start, n, smpi->getBufferSize(buffer_id), mass_ );
            }
            
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[1*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task for Push on ibin
        } // end ibin loop for Push

        // Particles BC and keys
        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
#endif
            {

#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif
            // Variables to compute cell_keys for the sorting
            unsigned int length[3];
            length[0]=0;
            length[1]=params.patch_size_[1]+1;
            length[2]=params.patch_size_[2]+1;


            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,2);
            for( int scell = first_cell_of_bin[ibin] ; scell <= last_cell_of_bin[ibin] ; scell++ ) {

                double energy_lost( 0. );
                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        ( *partWalls )[iwall]->apply( this, particles->first_index[ipack*packsize_+scell], particles->last_index[ipack*packsize_+scell], smpi->dynamics_invgf[buffer_id], patch->rand_, energy_lost );
                    }

                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore

                    partBoundCond->apply( this, particles->first_index[ipack*packsize_+scell], particles->last_index[ipack*packsize_+scell], smpi->dynamics_invgf[buffer_id], patch->rand_, energy_lost );
                    nrj_lost_per_bin[ibin] += mass_ * energy_lost;

                } else if( mass_==0 ) {

                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                      ( *partWalls )[iwall]->apply( this, particles->first_index[ipack*packsize_+scell], particles->last_index[ipack*packsize_+scell], smpi->dynamics_invgf[buffer_id], patch->rand_, energy_lost );
                        nrj_lost_per_bin[ibin] += energy_lost;
                    }

                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore

                    partBoundCond->apply( this, particles->first_index[ipack*packsize_+scell], particles->last_index[ipack*packsize_+scell], smpi->dynamics_invgf[buffer_id], patch->rand_, energy_lost );
                    nrj_lost_per_bin[ibin] += energy_lost;

                } // end if condition on mass

            } // end scell loop
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,2);


            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,11);
            for( int scell = first_cell_of_bin[ibin] ; scell <= last_cell_of_bin[ibin] ; scell++ ) {
                // Apply wall and boundary conditions
                if( mass_>0 ) {

                    for( int iPart=particles->first_index[ipack*packsize_+scell] ; ( int )iPart<particles->last_index[ipack*packsize_+scell]; iPart++ ) {
                        if ( particles->cell_keys[iPart] >= 0 ) {
                            //Compute cell_keys of remaining particles
                            for( unsigned int i = 0 ; i<nDim_field; i++ ) {
                                particles->cell_keys[iPart] *= this->length_[i];
                                particles->cell_keys[iPart] += round( ((this)->*(distance[i]))(particles, i, iPart) * dx_inv_[i] );
                            }
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            // count[particles->cell_keys[iPart]] ++;
                        }
                    }

                } else if( mass_==0 ) {

                    for( int iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                        if ( particles->cell_keys[iPart] >= 0 ) {
                            //Compute cell_keys of remaining particles
                            for( unsigned int i = 0 ; i<nDim_field; i++ ) {
                                particles->cell_keys[iPart] *= length[i];
                                particles->cell_keys[iPart] += round( ((this)->*(distance[i]))(particles, i, iPart) * dx_inv_[i] );
                            }
                            // count[particles->cell_keys[iPart]] ++;
                        }
                    }
                } // end if condition on mass

            } // end scell loop
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,11);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task for particles BC and cell_keys on ibin
        } // end ibin loop for particles BC


        // Project currents if not a Test species and charges as well if a diag is needed.
        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer) depend(in:bin_has_done_particles_BC[ibin])
#else
            #pragma omp task default(shared) firstprivate(ibin,bin_size0) depend(in:bin_has_done_particles_BC[ibin])
#endif
            {

#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif


            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                if (params.geometry != "AMcylindrical"){
                    Proj->currentsAndDensityWrapperOnBuffers(
                          b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin],
                          ibin*cluster_width_, *particles, smpi,
                          particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]],
                          buffer_id, diag_flag, params.is_spectral, ispec, 0 );
                } else {
                    // for( unsigned int scell = first_cell_of_bin[ibin] ; scell < last_cell_of_bin[ibin] ; scell++ ) {
                    // Proj->currentsAndDensityWrapperOnAMBuffers( EMfields, b_Jl[ibin], b_Jr[ibin], b_Jt[ibin], b_rhoAM[ibin],
                    //                                             ibin*cluster_width_, bin_size0, *particles, smpi,
                    //                                             particles->first_index[scell], particles->last_index[scell],
                    //                                             buffer_id, diag_flag);
                    // } // end scell loop
                } // end if AM
            } // end condition on test and mass
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[2*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            }//end task for Proj of ibin
        }// end ibin loop for Proj

//         // reduction of the lost energy in each ibin
//         // the dependency ensures that it is done after the particles BC
// #ifdef  __DETAILED_TIMERS
//         #pragma omp task default(shared) private(ithread,timer) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
// #else
//         #pragma omp task default(shared) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
// #endif
//         {
//         // reduce the energy lost with BC per bin
//         for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
//            nrj_bc_lost += nrj_lost_per_bin[ibin];
//         }
//
//         // sum the radiated energy / energy converted in pairs
//         // The dependencies above ensure that this is done after the Radiation and MultiPhoton Breit Wheeler methods
//         if( Radiate || Multiphoton_Breit_Wheeler_process) {
// #ifdef  __DETAILED_TIMERS
//             timer = MPI_Wtime();
//             ithread = Tools::getOMPThreadNum();
// #endif
//
//             for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
//                radiated_energy += radiated_energy_per_bin[ibin];
//             }
// #ifdef  __DETAILED_TIMERS
//             patch->patch_timers_[5*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
// #endif
//         } // end if Radiate or Multiphoton_Breit_Wheeler_process
//         } // end task for lost/radiated energy reduction



    } //End if moving particles

    if(time_dual <= time_frozen_ && diag_flag &&( !particles->is_test ) ) { //immobile particle (at the moment only project density)

        if (Ionize){ // project only after ionizing
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer) depend(in:bin_has_interpolated[ibin])
#else
                #pragma omp task default(shared) firstprivate(ibin,bin_size0) depend(in:bin_has_interpolated[ibin])
#endif
                {

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
                for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;

#ifdef  __DETAILED_TIMERS
                ithread = Tools::getOMPThreadNum();
                timer = MPI_Wtime();
#endif
                for( int iPart=particles->first_index[first_cell_of_bin[ibin]] ; iPart<particles->last_index[last_cell_of_bin[ibin]]; iPart++ ) {
                    Proj->basic( b_rho[ibin], ( *particles ), iPart, 0, ibin*cluster_width_ );
                }
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[2*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
                } // end task
            } // end ibin
        } else { // without ionization, no dependency is needed for the projection of frozen particles

            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin,bin_size0)
#endif
                {

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
                for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
#ifdef  __DETAILED_TIMERS
                ithread = Tools::getOMPThreadNum();
                timer = MPI_Wtime();
#endif
                for( int iPart=particles->first_index[first_cell_of_bin[ibin]] ; iPart<particles->last_index[last_cell_of_bin[ibin]]; iPart++ ) {
                    Proj->basic( b_rho[ibin], ( *particles ), iPart, 0, ibin*cluster_width_ );
                }
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[2*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
                } // end task
            } // end ibin
        } // end if Ionize
    } // end frozen and not test particle
    } // end taskgroup

    if (time_dual>time_frozen_){
        // reduction of the lost energy in each ibin
        // the dependency ensures that it is done after the particles BC
// #ifdef  __DETAILED_TIMERS
//         #pragma omp task default(shared) private(ithread,timer) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
// #else
//         #pragma omp task default(shared) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
// #endif
#ifdef  __DETAILED_TIMERS
        #pragma omp task default(shared) private(ithread,timer)
#else
        #pragma omp task default(shared)
#endif
        {
        // reduce the energy lost with BC per bin
        for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
           nrj_bc_lost += nrj_lost_per_bin[ibin];
        }

        // sum the radiated energy / energy converted in pairs
        // The dependencies above ensure that this is done after the Radiation and MultiPhoton Breit Wheeler methods
        if( Radiate || Multiphoton_Breit_Wheeler_process) {
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
            ithread = Tools::getOMPThreadNum();
#endif

            for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
               nrj_radiated_ += radiated_energy_per_bin[ibin];
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[5*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
        } // end if Radiate or Multiphoton_Breit_Wheeler_process
        } // end task for lost/radiated energy reduction
    } // end if moving particle

} // end scalarDynamicsTasks
#endif //endif tasks are used

// -----------------------------------------------------------------------------
//! Compute part_cell_keys at patch creation.
//! This operation is normally done in the pusher to avoid additional particles pass.
// -----------------------------------------------------------------------------
/*void SpeciesVAdaptive::computeParticleCellKeys(Params &params)
{

    unsigned int ip, nparts;
    int IX;
    double X;
    unsigned int length[3];

    //Number of particles before exchange
    nparts = particles->size();

    // Cell_keys is resized at the current number of particles
    particles->cell_keys.resize(nparts);

    length[0]=0;
    length[1]=params.patch_size_[1]+1;
    length[2]=params.patch_size_[2]+1;

    #pragma omp simd
    for (ip=0; ip < nparts ; ip++){
    // Counts the # of particles in each cell (or sub_cell) and store it in sparticles->last_index.
        for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
            X = particles->position(ipos,ip)-min_loc_vec[ipos];
            IX = round(X * dx_inv_[ipos] );
            particles->cell_keys[ip] = particles->cell_keys[ip] * length[ipos] + IX;
        }
    }

    // Reduction of the number of particles per cell in count
    for (ip=0; ip < nparts ; ip++)
        count[particles->cell_keys[ip]] ++ ;
}*/


// -----------------------------------------------------------------------------
//! This function reconfigures the type of species according
//! to the vectorization mode
// -----------------------------------------------------------------------------
void SpeciesVAdaptive::reconfiguration( Params &params, Patch *patch )
{

    //unsigned int ncell;
    bool reasign_operators = false;
    //float ratio_number_of_vecto_cells;
    float vecto_time = 0.;
    float scalar_time = 0.;

    //split cell into smaller sub_cells for refined sorting
    // cell = (params.patch_size_[0]+1);
    //for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.patch_size_[i]+1);

    // --------------------------------------------------------------------
    // Metrics 1 - based on the ratio of vectorized cells
    // Compute the number of cells that contain more than 8 particles
    //ratio_number_of_vecto_cells = SpeciesMetrics::get_ratio_number_of_vecto_cells(count,8);

    // Test metrics, if necessary we reasign operators
    //if ( (ratio_number_of_vecto_cells > 0.5 && this->vectorized_operators == false)
    //  || (ratio_number_of_vecto_cells < 0.5 && this->vectorized_operators == true))
    //{
    //    reasign_operators = true;
    //}
    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // Metrics 2 - based on the evaluation of the computational time
    (*part_comp_time_)( count,
                    vecto_time,
                    scalar_time );

    if( ( vecto_time <= scalar_time && this->vectorized_operators == false )
            || ( vecto_time > scalar_time && this->vectorized_operators == true ) ) {
        reasign_operators = true;
    }
    // --------------------------------------------------------------------

    // Operator reasignment if required by the metrics
    if( reasign_operators ) {

        // The type of operator is changed
        this->vectorized_operators = !this->vectorized_operators;

        /*MESSAGE(1,"> Species " << this->name_ << " reconfiguration (" << this->vectorized_operators
                  << ") in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ")"
                  << " of MPI process "<< patch->MPI_me_);*/

        this->reconfigure_operators( params, patch );

    }

    /*std::cout << " bin number: " << particles->first_index.size()
              << " nb particles: " << particles->last_index[particles->last_index.size()-1]
              << '\n';*/

}

// -----------------------------------------------------------------------------
//! This function configures the type of species according to the default mode
//! regardless the number of particles per cell
// -----------------------------------------------------------------------------
void SpeciesVAdaptive::defaultConfigure( Params &params, Patch *patch )
{

    // Setup the species state regardless the number of particles per cell
    this->vectorized_operators = ( params.adaptive_default_mode == "on" );

    // Configure the species regardless the number of particles per cell
    this->reconfigure_operators( params, patch );

}

// -----------------------------------------------------------------------------
//! This function configures the type of species according
//! to the vectorization mode
// -----------------------------------------------------------------------------
void SpeciesVAdaptive::configuration( Params &params, Patch *patch )
{
    //float ratio_number_of_vecto_cells;
    float vecto_time = 0.;
    float scalar_time = 0.;

    // Species with particles
    if( particles->size() > 0 ) {

        // --------------------------------------------------------------------
        // Metrics 2 - based on the evaluation of the computational time
        (*part_comp_time_)( count,
                        vecto_time,
                        scalar_time );

        if( vecto_time <= scalar_time ) {
            this->vectorized_operators = true;
        } else if( vecto_time > scalar_time ) {
            this->vectorized_operators = false;
        }
    }

    // Default mode where there is no particles
    else {
        this->vectorized_operators = ( params.adaptive_default_mode == "on" );
    }

    // --------------------------------------------------------------------

    this->reconfigure_operators( params, patch );
}

// -----------------------------------------------------------------------------
//! This function reconfigures the species operators after evaluating
//! the best mode from the particle distribution
// -----------------------------------------------------------------------------
void SpeciesVAdaptive::reconfigure_operators( Params &params, Patch *patch )
{
    // Destroy current operators
    delete Interp;
    //delete Push;
    delete Proj;

    // Reassign the correct Interpolator
    Interp = InterpolatorFactory::create( params, patch, this->vectorized_operators );
    // Reassign the correct Pusher to Push
    //Push = PusherFactory::create(params, this);
    // Reassign the correct Projector
    Proj = ProjectorFactory::create( params, patch, this->vectorized_operators );
}


void SpeciesVAdaptive::scalarPonderomotiveUpdateSusceptibilityAndMomentum( double time_dual, 
        ElectroMagn *EMfields,
        Params &params, 
        Patch *patch, SmileiMPI *smpi )
{

    int ithread;
#ifdef _OPENMP
    ithread = Tools::getOMPThreadNum();
#else
    ithread = 0;
#endif

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->resizeBuffers( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif

        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);
        Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread );
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);

#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[7] += MPI_Wtime() - timer;
#endif

        // Ionization
        if (Ionize){
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
            vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
            vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );
            vector<double> *Phipart = &( smpi->dynamics_PHIpart[ithread] );

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,5);
            Interp->envelopeFieldForIonization( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread );
            Ionize->envelopeIonization( particles, ( particles->first_index[0] ), ( particles->last_index[particles->last_index.size()-1] ), Epart, EnvEabs_part, EnvExabs_part, Phipart, patch, Proj );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,5);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[4] += MPI_Wtime() - timer;
#endif
        }
        // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[0], particles->last_index.back(), ithread );
#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[8] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif

        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);
        // Push only the particle momenta
        ( *Push )( *particles, smpi, 0, particles->last_index.back(), ithread );
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[9] += MPI_Wtime() - timer;
#endif

    } else { // immobile particle
    } //END if time vs. time_frozen_
} // ponderomotiveUpdateSusceptibilityAndMomentum

#ifdef _OMPTASKS
void SpeciesVAdaptive::scalarPonderomotiveUpdateSusceptibilityAndMomentumTasks( double time_dual, 
        ElectroMagn *EMfields,
        Params &params, 
        Patch *patch, SmileiMPI *smpi, int buffer_id )
{

#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif
    int bin_size0 = b_dim[0];

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    // for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
    //     nrj_lost_per_bin[ibin] = 0.;
    //     radiated_energy_per_bin[ibin] = 0.;
    // }
    // Init tags for the task dependencies of the particle operations
    int *bin_has_interpolated                   = new int[Nbins];
    int *bin_has_ionized                        = new int[Nbins];
    int *bin_has_projected_chi                  = new int[Nbins];


    #pragma omp taskgroup
    {

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize) { // moving particle

        smpi->resizeBuffers( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );
        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) { // loop on ibin
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
            {
            if ( params.geometry != "AMcylindrical" ){
                // Reset susceptibility sub-buffer - this buffer stores a grid susceptibility on the ibin physical space
                // This must be done before Projection of susceptibility
                for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_Chi[ibin][i]   = 0.0;
            } else {
                // Reset susceptibility sub-buffer - this buffer stores a grid susceptibility on the ibin physical space
                // This must be done before Projection of susceptibility
                for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_ChiAM[ibin][i] = 0.0;
            }
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,0);
            // Interpolate the fields and envelope at the particle position
            Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( particles->first_index[first_cell_of_bin[ibin]] ), &( particles->last_index[last_cell_of_bin[ibin]] ), buffer_id );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[7*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task interp
        } // end ibin

        // Ionization
        if( Ionize ) {
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) { // loop on ibin
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
                {

#ifdef  __DETAILED_TIMERS
                ithread = Tools::getOMPThreadNum();
                timer = MPI_Wtime();
#endif
                vector<double> *Epart = &( smpi->dynamics_Epart[buffer_id] );
                vector<double> *EnvEabs_part = &( smpi->dynamics_EnvEabs_part[buffer_id] );
                vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[buffer_id] );
                vector<double> *Phipart = &( smpi->dynamics_PHIpart[buffer_id] );

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,5);
                Interp->envelopeFieldForIonization( EMfields, *particles, smpi, &( particles->first_index[first_cell_of_bin[ibin]] ), &( particles->last_index[last_cell_of_bin[ibin]] ), buffer_id );
                Ionize->envelopeIonization( particles, particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]], Epart, EnvEabs_part, EnvExabs_part, Phipart, patch, Proj, 0 );
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,5);

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
                } // end task ionize
            } // end ibin
        }  // end Ionize

        if( time_dual>time_frozen_) {
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) { // loop on ibin
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_projected_chi[ibin]) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_projected_chi[ibin])
#endif
                {
                // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
                ithread = Tools::getOMPThreadNum();
                timer = MPI_Wtime();
#endif

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
                if (params.geometry != "AMcylindrical"){
                    Proj->susceptibilityOnBuffer( EMfields, b_Chi[ibin],
                                                  ibin*cluster_width_, bin_size0,
                                                  *particles, mass_, smpi,
                                                  particles->first_index[ibin], particles->last_index[ibin],
                                                  buffer_id );
                } else { // vectorized envelope not implemente in AM
                    // Proj->susceptibilityOnBuffer( EMfields, b_ChiAM[ibin],
                    //                               ibin*cluster_width_, bin_size0,
                    //                               *particles, mass_, smpi,
                    //                               particles->first_index[ibin], particles->last_index[ibin],
                    //                               buffer_id );
                } // end if geometry
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);
#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[8*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif

                } // end task susceptibility
            } // end ibin
        }

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) { // loop on ibin
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_projected_chi[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_projected_chi[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);
            // Push only the particle momenta
            ( *Push )( *particles, smpi, particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]], buffer_id);
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[9*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task susceptibility
        } // end ibin

    } else { // immobile particle
    } //END if time vs. time_frozen_

    } // end taskgroup

} // ponderomotiveUpdateSusceptibilityAndMomentumTasks
#endif // endif tasks are used

void SpeciesVAdaptive::scalarPonderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag, PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi )
{

    int ithread;
#ifdef _OPENMP
    ithread = Tools::getOMPThreadNum();
#else
    ithread = 0;
#endif

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    int tid( 0 );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle updated position
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->resizeBuffers( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        // Interpolate the ponderomotive potential and its gradient at the particle position, present and previous timestep
#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif

        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,0);
        Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread );
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);

#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[10] += MPI_Wtime() - timer;
#endif

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif

        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);
        // Push only the particle position
        ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[0], particles->last_index.back(), ithread );
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[11] += MPI_Wtime() - timer;
#endif

        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,2);
        for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {
            double energy_lost( 0. );
            // Apply wall and boundary conditions
            if( mass_>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    ( *partWalls )[iwall]->apply( this, particles->first_index[scell], particles->last_index[scell], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                    nrj_lost_per_thd[tid] += mass_ * energy_lost;
                }

                // Boundary Condition may be physical or due to domain decomposition
                partBoundCond->apply( this, particles->first_index[scell], particles->last_index[scell], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                nrj_lost_per_thd[tid] += mass_ * energy_lost;

            } else if( mass_==0 ) {
                ERROR( "Particles with zero mass cannot interact with envelope" );
            } // end mass_ = 0? condition
        }
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,2);


        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,11);
        if( mass_>0 ) {
            computeParticleCellKeys( params );
        } else if( mass_==0 ) {
            ERROR_NAMELIST( "Particles with zero mass cannot interact with envelope",
            LINK_NAMELIST + std::string("#laser-envelope-model"));
        } // end mass_ = 0? condition
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,11);

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif

        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
        if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
            Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, particles->first_index[0], particles->last_index.back(), ithread, diag_flag, params.is_spectral, ispec );
        }
        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
        patch->patch_timers_[12] += MPI_Wtime() - timer;
#endif

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

    } // end case of moving particle
    else { // immobile particle


        smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
        if( diag_flag &&( !particles->is_test ) ) {
            double *b_rho=nullptr;
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                // only 3D is implemented actually
                if( nDim_field==2 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                }
                if( nDim_field==3 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                } else if( nDim_field==1 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                }
                for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                } //End loop on particles
            }//End loop on bins
          smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);
        } // end condition on diag and not particle test

    }//END if time vs. time_frozen_
} // End ponderomotive_position_update

#ifdef _OMPTASKS
void SpeciesVAdaptive::scalarPonderomotiveUpdatePositionAndCurrentsTasks( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag, PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi, int buffer_id )
{


#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif
    int bin_size0 = b_dim[0];

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
        nrj_lost_per_bin[ibin] = 0.;
        // radiated_energy_per_bin[ibin] = 0.;
    }
    // Init tags for the task dependencies of the particle operations
    int *bin_has_interpolated               = new int[Nbins];
    int *bin_has_pushed                     = new int[Nbins];
    int *bin_has_done_particles_BC          = new int[Nbins];

    #pragma omp taskgroup
    {

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ ) {

        smpi->resizeBuffers( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        //Prepare for sorting
        for( unsigned int i=0; i<count.size(); i++ ) {
            count[i] = 0;
        }

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(out:bin_has_interpolated[ibin])
#endif
            {
            if ( params.geometry != "AMcylindrical" ){
                // Reset densities sub-buffers - each of these buffers store a grid density on the ibin physical space
                // This must be done before Projection and before Ionization (because of the ionization currents)
                for (unsigned int i = 0; i < size_proj_buffer_Jx; i++)  b_Jx[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jy; i++)  b_Jy[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jz; i++)  b_Jz[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
            } else {
                // Reset densities sub-buffers - each of these buffers store a grid density on the ibin physical space
                // This must be done before Projection and before Ionization (because of the ionization currents)
                for (unsigned int i = 0; i < size_proj_buffer_Jl; i++)  b_Jl[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jr; i++)  b_Jr[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_Jt; i++)  b_Jt[ibin][i]    = 0.0;
                for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
            }
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,0);
            Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( particles->first_index[first_cell_of_bin[ibin]] ), &( particles->last_index[last_cell_of_bin[ibin]] ), buffer_id );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[10*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task interpolate
        } // end ibin

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_pushed[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_pushed[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]], buffer_id );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[11*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task
        } // end ibin

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif
            // Apply wall and boundary conditions
            if( mass_>0 ) {

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,2);
                double energy_lost=0;
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    ( *partWalls )[iwall]->apply( this, particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]], smpi->dynamics_invgf[buffer_id], patch->rand_, energy_lost );
                    nrj_lost_per_bin[ibin] += mass_ * energy_lost;
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                partBoundCond->apply( this, particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]], smpi->dynamics_invgf[last_cell_of_bin[ibin]], patch->rand_, energy_lost );
                nrj_lost_per_bin[ibin] += mass_ * energy_lost;
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,2);


                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,11);
                for( int iPart=particles->first_index[first_cell_of_bin[ibin]] ; iPart<particles->last_index[last_cell_of_bin[ibin]]; iPart++ ) {
                    if ( particles->cell_keys[iPart] >= 0 ) {
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        for( int i = 0 ; i<( int )nDim_field; i++ ) {
                            particles->cell_keys[iPart] *= length_[i];
                            particles->cell_keys[iPart] += round( ((this)->*(distance[i]))(particles, i, iPart) * dx_inv_[i] );
                        }
                        // count[particles->cell_keys[iPart]] ++; //First reduction of the count sort algorithm. Lost particles are not included.
                    }
                } // end iPart loop
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,11);

            } else if( mass_==0 ) {
                ERROR( "Particles with zero mass cannot interact with envelope" );

            } // end mass_ = 0? condition
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task
        } // end ibin

        for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_done_particles_BC[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_done_particles_BC[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = Tools::getOMPThreadNum();
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                if (params.geometry != "AMcylindrical"){
                    Proj->currentsAndDensityWrapperOnBuffers( b_Jx[ibin], b_Jy[ibin], b_Jz[ibin], b_rho[ibin],
                                                              ibin*cluster_width_, *particles, smpi,
                                                              particles->first_index[first_cell_of_bin[ibin]], particles->last_index[last_cell_of_bin[ibin]],
                                                              buffer_id, diag_flag, params.is_spectral, ispec );
                  } else { // vectorized envelope not implemented in AM
                    // Proj->currentsAndDensityWrapperOnAMBuffers( EMfields, b_Jl[ibin], b_Jr[ibin], b_Jt[ibin], b_rhoAM[ibin],
                    //                                             ibin*cluster_width_, bin_size0, *particles, smpi,
                    //                                             particles->first_index[ibin], particles->last_index[ibin],
                    //                                             buffer_id, diag_flag);
                  } // end if AM
            } // end condition on test and mass
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task
        } // end ibin

        // // reduction of the lost energy in each ibin
        // // the dependency ensures that it is done after the particles BC
        // //#pragma omp task default(shared) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
        // // using depend out on particles BC a segfault is caused - check why this happens
        // #pragma omp task default(shared) depend(in:bin_has_projected[0:(Nbins-1)])
        // {
        // // reduce the energy lost with BC per bin
        // for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
        //     nrj_bc_lost += nrj_lost_per_bin[ibin];
        // } // end ibin
        // } // end task for lost/radiated energy reduction

    } else { // immobile particle
        if( diag_flag &&( !particles->is_test ) ) {
            if( params.geometry != "AMcylindrical" ) {

                for( unsigned int ibin = 0 ; ibin < Nbins ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                    #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer)
#else
                    #pragma omp task default(shared) firstprivate(ibin)
#endif
                    {
#ifdef  __DETAILED_TIMERS
                    ithread = Tools::getOMPThreadNum();
                    timer = MPI_Wtime();
#endif


                    smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
                    for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
                    for( int iPart=particles->first_index[first_cell_of_bin[ibin]] ; iPart<particles->last_index[last_cell_of_bin[ibin]]; iPart++ ) {
                        Proj->basic( b_rho[ibin], ( *particles ), iPart, 0, ibin*cluster_width_ );
                    }
                    smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[3*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
#endif
                    } // end task projection for frozen or test
                } // end ibin
            } else { // vectorized envelope not implemented in AM geometry

//                 for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
// #ifdef  __DETAILED_TIMERS
//                     #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer)
// #else
//                     #pragma omp task default(shared) firstprivate(ibin,bin_size0)
// #endif
//                     {
// #ifdef  __DETAILED_TIMERS
//                     ithread = Tools::getOMPThreadNum();
//                     timer = MPI_Wtime();
// #endif
//                     for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
//                     int imode = 0; // only mode 0 is used with envelope
//                     for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
//                         Proj->basicForComplexOnBuffer( b_rhoAM[ibin], ( *particles ), iPart, 0, imode, bin_size0, ibin*cluster_width_ );
//                     } // end loop on particles
// #ifdef  __DETAILED_TIMERS
//                     patch->patch_timers_[3*patch->number_of_threads_ + ithread] += MPI_Wtime() - timer;
// #endif
//                     } // end task projection for frozen or test
//                 } //end ibin
            } // end if on geometry
        } // end condition on diag and not particle test

    } // end if moving particle

    } // end taskgroup

    if (time_dual>time_frozen_){
        // reduction of the lost energy in each ibin
        // the dependency ensures that it is done after the particles BC
        //#pragma omp task default(shared) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
        // using depend out on particles BC a segfault is caused - check why this happens
        // #pragma omp task default(shared) depend(in:bin_has_projected[0:(Nbins-1)])
        #pragma omp task default(shared)
        {
        // reduce the energy lost with BC per bin
        for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
            nrj_bc_lost += nrj_lost_per_bin[ibin];
        } // end ibin
        } // end task for lost/radiated energy reduction
    } // end if moving particle


} // End ponderomotive_position_updateTasks
#endif // endif tasks are used
