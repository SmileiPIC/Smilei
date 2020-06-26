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

#include "SpeciesMetrics.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
SpeciesVAdaptive::SpeciesVAdaptive( Params &params, Patch *patch ) :
    SpeciesV( params, patch )
{
    initCluster( params );
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
        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
        vector<Diagnostic *> &localDiags )
{
    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    // Reset list of particles to exchange
    clearExchList();

    int tid( 0 );
    double ener_iPart( 0. );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize ) {
        // moving particle

        smpi->dynamics_resize( ithread, nDim_particle, particles->last_index.back() );

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

        // Interpolate the fields at the particle position
        Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread, particles->first_index[0] );

#ifdef  __DETAILED_TIMERS
        patch->patch_timers[0] += MPI_Wtime() - timer;
#endif

        // Interpolate the fields at the particle position
        //for (unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++)
        //    (*Interp)(EMfields, *particles, smpi, &(particles->first_index[scell]), &(particles->last_index[scell]), ithread );
        for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {

            // Ionization
            if( Ionize ) {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                ( *Ionize )( particles, particles->first_index[scell], particles->last_index[scell], Epart, patch, Proj );
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
            }

            if( time_dual<=time_frozen_ ) continue; // Do not push nor project frozen particles

            // Radiation losses
            if( Radiate ) {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                // Radiation process
                ( *Radiate )( *particles, this->photon_species_, smpi,
                              RadiationTables, nrj_radiation,
                              particles->first_index[scell], particles->last_index[scell], ithread );

                // // Update scalar variable for diagnostics
                // nrj_radiation += Radiate->getRadiatedEnergy();
                //
                // // Update the quantum parameter chi
                // Radiate->computeParticlesChi( *particles,
                //                               smpi,
                //                               first_index[scell],
                //                               last_index[scell],
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
                                                        particles->first_index[scell], particles->last_index[scell], ithread );

                // Update scalar variable for diagnostics
                // We reuse nrj_radiation for the pairs
                nrj_radiation += Multiphoton_Breit_Wheeler_process->getPairEnergy();

                // Update the photon quantum parameter chi of all photons
                Multiphoton_Breit_Wheeler_process->compute_thread_chiph( *particles,
                        smpi,
                        particles->first_index[scell],
                        particles->last_index[scell],
                        ithread );

                // Suppression of the decayed photons into pairs
                Multiphoton_Breit_Wheeler_process->decayed_photon_cleaning(
                    *particles, smpi, scell, particles->first_index.size(), &particles->first_index[0], &particles->last_index[0], ithread );
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[6] += MPI_Wtime() - timer;
#endif
            }
        }

    if( time_dual>time_frozen_ ) { // do not push, nor apply particles BC, nor project frozen particles

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push the particles and the photons
            ( *Push )( *particles, smpi, 0, particles->last_index.back(), ithread, 0. );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[1] += MPI_Wtime() - timer;
            timer = MPI_Wtime();
#endif

            // Computation of the particle cell keys for all particles
            // this->compute_bin_cell_keys(params,0, particles->last_index.back());

            for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {
                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                            }
                        }
                    }

                    for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                            particles->cell_keys[iPart] = -1;
                        } else {
                            //Compute cell_keys of remaining particles
                            for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
                                particles->cell_keys[iPart] *= this->length_[i];
                                particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
                            }
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            count[particles->cell_keys[iPart]] ++;
                        }
                    }

                } else if( mass_==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += ener_iPart;
                            }
                        }
                    }

                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += ener_iPart;
                            particles->cell_keys[iPart] = -1;
                        } else {
                             //Compute cell_keys of remaining particles
                            for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
                                particles->cell_keys[iPart] *= this->length_[i];
                                particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
                            }
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            count[particles->cell_keys[iPart]] ++;
                        }
                    }
                } // end if mass_ > 0
            } // end loop on cells

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[3] += MPI_Wtime() - timer;
#endif

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                Proj->currentsAndDensityWrapper(
                    EMfields, *particles, smpi, particles->first_index[0],
                    particles->last_index.back(),
                    ithread,
                    diag_flag,
                    params.is_spectral,
                    ispec
            );
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[2] += MPI_Wtime() - timer;
#endif

            }
        } // end if moving particle

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

    }  // end if moving particle or ionize

    if(time_dual <= time_frozen_ && diag_flag &&( !particles->is_test ) ) { //immobile particle (at the moment only project density)

        double *b_rho=nullptr;
        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj

            b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
            for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                Proj->basic( b_rho, ( *particles ), iPart, 0 );
            }
        }
    }

}//END scalarDynamics


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
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;

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
    // cell = (params.n_space[0]+1);
    //for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.n_space[i]+1);

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
    SpeciesMetrics::get_computation_time( count,
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
        SpeciesMetrics::get_computation_time( this->count,
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


void SpeciesVAdaptive::scalarPonderomotiveUpdateSusceptibilityAndMomentum( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{

    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->dynamics_resize( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread );
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[7] += MPI_Wtime() - timer;
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
            Interp->envelopeFieldForIonization( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread );
            Ionize->envelopeIonization( particles, ( particles->first_index[0] ), ( particles->last_index[particles->last_index.size()-1] ), Epart, EnvEabs_part, EnvExabs_part, Phipart, patch, Proj );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
        }
        // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[0], particles->last_index.back(), ithread );
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[8] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        // Push only the particle momenta
        ( *Push )( *particles, smpi, 0, particles->last_index.back(), ithread );
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[9] += MPI_Wtime() - timer;
#endif

    } else { // immobile particle
    } //END if time vs. time_frozen_
} // ponderomotiveUpdateSusceptibilityAndMomentum


void SpeciesVAdaptive::scalarPonderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag, PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{

    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    // Reset list of particles to exchange - WARNING Should it be reset?
    clearExchList();

    int tid( 0 );
    double ener_iPart( 0. );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle updated position
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->dynamics_resize( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        //Prepare for sorting
        for( unsigned int i=0; i<count.size(); i++ ) {
            count[i] = 0;
        }



        // Interpolate the ponderomotive potential and its gradient at the particle position, present and previous timestep
#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread );
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[10] += MPI_Wtime() - timer;
#endif

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        // Push only the particle position
        ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[0], particles->last_index.back(), ithread );
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[11] += MPI_Wtime() - timer;
#endif

        for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {
            // Apply wall and boundary conditions
            if( mass_>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                            nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for( iPart=particles->first_index[scell] ; ( int )iPart<particles->last_index[scell]; iPart++ ) {
                    if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                        particles->cell_keys[iPart] = -1;
                    } else {
                        //Compute cell_keys of remaining particles
                        for( unsigned int i = 0 ; i<nDim_particle; i++ ) {
                            particles->cell_keys[iPart] *= this->length_[i];
                            particles->cell_keys[iPart] += round( ( particles->position( i, iPart )-min_loc_vec[i] ) * dx_inv_[i] );
                        }
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        count[particles->cell_keys[iPart]] ++;
                    }

                }

            } else if( mass_==0 ) {
                ERROR( "Particles with zero mass cannot interact with envelope" );
            } // end mass_ = 0? condition
        }

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
            Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, particles->first_index[0], particles->last_index.back(), ithread, diag_flag, params.is_spectral, ispec );
        }
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[12] += MPI_Wtime() - timer;
#endif

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

    } // end case of moving particle
    else { // immobile particle

        if( Ionize ) {
            smpi->dynamics_resize( ithread, nDim_particle, particles->last_index.back() );

            //Point to local thread dedicated buffers
            //Still needed for ionization
            vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[0] ), &( particles->last_index[particles->last_index.size()-1] ), ithread, particles->first_index[0] );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[0] += MPI_Wtime() - timer;
#endif

            // Interpolate the fields at the particle position
            //for (unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++)
            //    (*Interp)(EMfields, *particles, smpi, &(particles->first_index[scell]), &(particles->last_index[scell]), ithread );
            for( unsigned int scell = 0 ; scell < particles->first_index.size() ; scell++ ) {

                // Ionization
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                ( *Ionize )( particles, particles->first_index[scell], particles->last_index[scell], Epart, patch, Proj );
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
            }// end loop on scells
        }// end if ionize

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
        } // end condition on diag and not particle test

    }//END if time vs. time_frozen_
} // End ponderomotive_position_update
