
#include "Species.h"

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

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species::Species( Params &params, Patch *patch ) :
    c_part_max_( 1 ),
    ionization_rate_( Py_None ),
    pusher_name_( "boris" ),
    radiation_model_( "none" ),
    time_frozen_( 0 ),
    radiating_( false ),
    relativistic_field_initialization_( false ),
    time_relativistic_initialization_( 0 ),
    multiphoton_Breit_Wheeler_( 2, "" ),
    ionization_model( "none" ),
    density_profile_type_( "none" ),
    charge_profile_( NULL ),
    density_profile_( NULL ),
    velocity_profile_( 3, NULL ),
    temperature_profile_( 3, NULL ),
    particles_per_cell_profile_( NULL ),
    max_charge_( 0. ),
    particles( &particles_sorted[0] ),
    position_initialization_array_( NULL ),
    momentum_initialization_array_( NULL ),
    n_numpy_particles_( 0 ),
    position_initialization_on_species_( false ),
    position_initialization_on_species_index( -1 ),
    electron_species( NULL ),
    electron_species_index( -1 ),
    photon_species_( NULL ),
    //photon_species_index(-1),
    radiation_photon_species( "" ),
    mBW_pair_creation_sampling_( 2, 1 ),
    clrw( params.clrw ),
    oversize( params.oversize ),
    cell_length( params.cell_length ),
    min_loc_vec( patch->getDomainLocalMin() ),
    tracking_diagnostic( 10000 ),
    nDim_particle( params.nDim_particle ),
    nDim_field(    params.nDim_field  ),
    merging_time_selection_( 0 )
{
    regular_number_array_.clear();
    partBoundCond = NULL;
    min_loc = patch->getDomainLocalMin( 0 );
    merging_method_ = "none";

    PI2 = 2.0 * M_PI;
    PI_ov_2 = 0.5*M_PI;

    dx_inv_[0] = 1./cell_length[0];
    dx_inv_[1] = 1./cell_length[1];
    dx_inv_[2] = 1./cell_length[2];

    initCluster( params );
    inv_nDim_particles = 1./( ( double )nDim_particle );

    length_[0]=0;
    length_[1]=params.n_space[1]+1;
    length_[2]=params.n_space[2]+1;

    merge_momentum_cell_size_.resize(3);

    merge_min_momentum_cell_length_.resize(3);

}//END Species creator

void Species::initCluster( Params &params )
{
    // Arrays of the min and max indices of the particle bins
    particles->first_index.resize( params.n_space[0]/clrw );
    particles->last_index.resize( params.n_space[0]/clrw );

    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)

    //Primal dimension of fields.
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;

    b_dim.resize( params.nDim_field, 1 );
    if( nDim_particle == 1 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0];
        f_dim1 = 1;
        f_dim2 = 1;
    }
    if( nDim_particle == 2 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] =  f_dim1;
        f_dim2 = 1;
    }
    if( nDim_particle == 3 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] = f_dim1;
        b_dim[2] = f_dim2;
    }

    //Initialize specMPI
    MPI_buffer_.allocate( nDim_field );

    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    new_particles_energy_ = 0.;
    nrj_radiation = 0.;

}//END initCluster

// -----------------------------------------------------------------------------
//! This function enables to resize the number of bins
// -----------------------------------------------------------------------------
void Species::resizeCluster( Params &params )
{

    // We keep the current number of particles
    int npart = particles->size();
    int size = params.n_space[0]/clrw;

    // Arrays of the min and max indices of the particle bins
    particles->first_index.resize( size );
    particles->last_index.resize( size );

    // We redistribute the particles between the bins
    int quotient = npart / size; // Fixed part for all bin
    int remainder = npart - quotient*size; // To be distributed among the first bin

    for( int ibin=0 ; ibin < size ; ibin++ ) {
        if( ibin < remainder ) {
            particles->first_index[ibin] = ibin*quotient + ibin;
            particles->last_index[ibin] = particles->first_index[ibin] + quotient + 1;
        } else {
            particles->first_index[ibin] = ibin*quotient + remainder;
            particles->last_index[ibin] = particles->first_index[ibin] + quotient;
        }
    }

    //std::cout << "size: " << size << " " << npart << " " << particles->first_index[0] << " " << particles->last_index[0] << '\n';

    // Recommended: A sorting process may be needed for best porfermance after this step

}// end resizeCluster

// Initialize the operators (Push, Ionize, PartBoundCond)
// This must be separate from the parameters because the Species cloning copies
// the parameters but not the operators.
void Species::initOperators( Params &params, Patch *patch )
{

    // interpolation operator (virtual)
    Interp = InterpolatorFactory::create( params, patch, this->vectorized_operators && !params.cell_sorting ); // + patchId -> idx_domain_begin (now = ref smpi)
    
    // assign the correct Pusher to Push
    Push = PusherFactory::create( params, this );
    if( this->ponderomotive_dynamics ) {
        Push_ponderomotive_position = PusherFactory::create_ponderomotive_position_updater( params, this );
    }

    // projection operator (virtual)
    Proj = ProjectorFactory::create( params, patch, this->vectorized_operators && !params.cell_sorting );  // + patchId -> idx_domain_begin (now = ref smpi)
    
    // Assign the Ionization model (if needed) to Ionize
    //  Needs to be placed after ParticleCreator() because requires the knowledge of max_charge_
    // \todo pay attention to restart
    Ionize = IonizationFactory::create( params, this );

    // Create the radiation model
    Radiate = RadiationFactory::create( params, this, patch->rand_ );

    // Create the multiphoton Breit-Wheeler model
    Multiphoton_Breit_Wheeler_process = MultiphotonBreitWheelerFactory::create( params, this, patch->rand_  );

    // assign the correct Merging method to Merge
    Merge = MergingFactory::create( params, this, patch->rand_ );

    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, this, patch );
    for( unsigned int iDim=0 ; iDim < nDim_field ; iDim++ ) {
        for( unsigned int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++ ) {
            MPI_buffer_.partRecv[iDim][iNeighbor].initialize( 0, ( *particles ) );
            MPI_buffer_.partSend[iDim][iNeighbor].initialize( 0, ( *particles ) );
            MPI_buffer_.part_index_send[iDim][iNeighbor].resize( 0 );
            MPI_buffer_.part_index_recv_sz[iDim][iNeighbor] = 0;
            MPI_buffer_.part_index_send_sz[iDim][iNeighbor] = 0;
        }
    }
    typePartSend.resize( nDim_field*2, MPI_DATATYPE_NULL );
    typePartRecv.resize( nDim_field*2, MPI_DATATYPE_NULL );
    exchangePatch = MPI_DATATYPE_NULL;

}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete Push;
    delete Interp;
    delete Proj;

    if( Merge ) {
        delete Merge;
    }

    if( Ionize ) {
        delete Ionize;
    }
    if( Radiate ) {
        delete Radiate;
    }
    if( Multiphoton_Breit_Wheeler_process ) {
        delete Multiphoton_Breit_Wheeler_process;
    }
    if( partBoundCond ) {
        delete partBoundCond;
    }
    if( particles_per_cell_profile_ ) {
        delete particles_per_cell_profile_;
    }
    if( charge_profile_ ) {
        delete charge_profile_;
    }
    if( density_profile_ ) {
        delete density_profile_;
    }
    for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
        delete velocity_profile_[i];
    }
    for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
        delete temperature_profile_[i];
    }
    if( ionization_rate_!=Py_None ) {
        Py_DECREF( ionization_rate_ );
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - perform ionization
//   - perform the radiation reaction
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics( double time_dual, unsigned int ispec,
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

    // Reset list of particles to exchange
    clearExchList();

    double ener_iPart( 0. );
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

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                            }
                        }
                    }
                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                            //}
                            //else if ( partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            //std::cout<<"removed particle position"<< particles->position(0,iPart)<<" , "<<particles->position(1,iPart)<<" ,"<<particles->position(2,iPart)<<std::endl;
                        }
                    }

                } else if( mass_==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                            double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                            if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                                nrj_lost_per_thd[tid] += ener_iPart;
                            }
                        }
                    }

                    // Boundary Condition may be physical or due to domain decomposition
                    // apply returns 0 if iPart is not in the local domain anymore
                    //        if omp, create a list per thread
                    for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                        if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                            addPartInExchList( iPart );
                            nrj_lost_per_thd[tid] += ener_iPart;
                        }
                    }

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
} //END dynamics


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - perform ionization
//   - perform the radiation reaction
//   - perform the multiphoton Breit-Wheeler
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::scalarDynamics( double time_dual, unsigned int ispec,
                               ElectroMagn *EMfields,
                               Params &params, bool diag_flag,
                               PartWalls *partWalls,
                               Patch *patch, SmileiMPI *smpi,
                               RadiationTables &RadiationTables,
                               MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                               vector<Diagnostic *> &localDiags )
{

}

void Species::projectionForDiags( double time_dual, unsigned int ispec,
                                    ElectroMagn *EMfields,
                                    Params &params, bool diag_flag,
                                    Patch *patch, SmileiMPI *smpi )
{
    if( diag_flag &&( !particles->is_test ) ) {

        if( params.geometry != "AMcylindrical" ) {
            double *buf[4];

            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj

                buf[0] = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                buf[1] = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
                buf[2] = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
                buf[3] = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;

                for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                    for( unsigned int quantity=0; quantity < 4; quantity++ ) {
                        Proj->basic( buf[quantity], ( *particles ), iPart, quantity );
                    }
                } //End loop on particles
            }//End loop on bins
        } else { // AM case
            complex<double> *buf[4];
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            int n_species = patch->vecSpecies.size();
            for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                int ifield = imode*n_species+ispec;

                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj

                    buf[0] = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                    buf[1] = emAM->Jl_s [ifield] ? &( *emAM->Jl_s [ifield] )( 0 ) : &( *emAM->Jl_[imode] )( 0 ) ;
                    buf[2] = emAM->Jr_s [ifield] ? &( *emAM->Jr_s [ifield] )( 0 ) : &( *emAM->Jr_[imode] )( 0 ) ;
                    buf[3] = emAM->Jt_s [ifield] ? &( *emAM->Jt_s [ifield] )( 0 ) : &( *emAM->Jt_[imode] )( 0 ) ;

                    for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                        for( unsigned int quantity=0; quantity < 4; quantity++ ) {
                            Proj->basicForComplex( buf[quantity], ( *particles ), iPart, quantity, imode );
                        }
                    } //End loop on particles
                }//End loop on bins
            } //End loop on modes
        }
        
    }
}

// -----------------------------------------------------------------------------
//! For all particles of the species, import the new particles generated
//! from these different physical processes:
//! - ionization
//! - radiation reaction
//! - multiphoton Breit-Wheeler
// -----------------------------------------------------------------------------
void Species::dynamicsImportParticles( double time_dual, unsigned int ispec,
        Params &params,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{
    // Add the ionized electrons to the electron species (possible even if ion is frozen)
    if( Ionize ) {
        electron_species->importParticles( params, patch, Ionize->new_electrons, localDiags );
                 }

    // if moving particle
    if( time_dual>time_frozen_ ) { // moving particle

        // Radiation losses
        if( Radiate ) {
            // If creation of macro-photon, we add them to photon_species
            if( photon_species_ ) {
                photon_species_->importParticles( params,
                                                 patch,
                                                 Radiate->new_photons_,
                                                 localDiags );
            }
        }

        // Multiphoton Breit-Wheeler
        if( Multiphoton_Breit_Wheeler_process ) {
            // Addition of the electron-positron particles
            for( int k=0; k<2; k++ ) {
                mBW_pair_species[k]->importParticles( params,
                                                      patch,
                                                      Multiphoton_Breit_Wheeler_process->new_pair[k],
                                                      localDiags );
            }
        }
    }//END if time vs. time_frozen_
}




// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - increment the charge (projection)
//   - used at initialisation for Poisson (and diags if required, not for now dynamics )
// ---------------------------------------------------------------------------------------------------------------------
void Species::computeCharge( unsigned int ispec, ElectroMagn *EMfields )
{
    // -------------------------------
    // calculate the particle charge
    // -------------------------------
    if( ( !particles->is_test ) ) {
        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
            // Not for now, else rho is incremented twice. Here and dynamics. Must add restartRhoJs and manage independantly diags output
            //b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(bin_start) : &(*EMfields->rho_)(bin_start);
            if( !dynamic_cast<ElectroMagnAM *>( EMfields ) ) {
                double *b_rho = &( *EMfields->rho_ )( 0 );

                for( unsigned int iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            } else {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
                unsigned int Nmode = emAM->rho_AM_.size();
                for( unsigned int imode=0; imode<Nmode; imode++ ) {
                    complex<double> *b_rho = &( *emAM->rho_AM_[imode] )( 0 );
                    for( unsigned int iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }
        }

    }
}//END computeCharge


// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::sortParticles( Params &params, Patch * patch )
{

    int ndim = params.nDim_field;
    int idim;
    //cleanupSentParticles(ispec, indexes_of_particles_to_exchange);

    //We have stored in indexes_of_particles_to_exchange the list of all particles that needs to be removed.
    /********************************************************************************/
    // Delete Particles included in the index of particles to exchange. Assumes indexes are sorted.
    /********************************************************************************/
    int ii, iPart;

    // Push lost particles at the end of bins
    for( unsigned int ibin = 0 ; ibin < particles->last_index.size() ; ibin++ ) {
        ii = indexes_of_particles_to_exchange.size()-1;
        if( ii >= 0 ) { // Push lost particles to the end of the bin
            iPart = indexes_of_particles_to_exchange[ii];
            while( iPart >= particles->last_index[ibin] && ii > 0 ) {
                ii--;
                iPart = indexes_of_particles_to_exchange[ii];
            }
            while( iPart == particles->last_index[ibin]-1 && iPart >= particles->first_index[ibin] && ii > 0 ) {
                particles->last_index[ibin]--;
                ii--;
                iPart = indexes_of_particles_to_exchange[ii];
            }
            while( iPart >= particles->first_index[ibin] && ii > 0 ) {
                particles->overwriteParticle( particles->last_index[ibin]-1, iPart );
                particles->last_index[ibin]--;
                ii--;
                iPart = indexes_of_particles_to_exchange[ii];
            }
            //On traite la dernière particule (qui peut aussi etre la premiere)
            if( iPart >= particles->first_index[ibin] && iPart < particles->last_index[ibin] ) {
                particles->overwriteParticle( particles->last_index[ibin]-1, iPart );
                particles->last_index[ibin]--;
            }
        }
    }


    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for( int unsigned ibin = 1 ; ibin < particles->last_index.size() ; ibin++ ) { //First bin don't need to be shifted
        ii = particles->first_index[ibin]-particles->last_index[ibin-1]; // Shift the bin in memory by ii slots.
        iPart = min( ii, particles->last_index[ibin]-particles->first_index[ibin] ); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
        if( iPart > 0 ) {
            particles->overwriteParticle( particles->last_index[ibin]-iPart, particles->last_index[ibin-1], iPart );
        }
        particles->last_index[ibin] -= ii;
        particles->first_index[ibin] = particles->last_index[ibin-1];
    }



    int nmove, lmove; // local, OK
    int shift[particles->last_index.size()+1];//how much we need to shift each bin in order to leave room for the new particle
    double dbin;

    dbin = params.cell_length[0]*params.clrw; //width of a bin.
    for( unsigned int j=0; j<particles->last_index.size()+1 ; j++ ) {
        shift[j]=0;
    }


    int nbNeighbors_ = 2;
    int n_part_recv;

    indexes_of_particles_to_exchange.clear();
    particles->eraseParticleTrail( particles->last_index.back() );

    //Evaluation of the necessary shift of all bins.2
    //idim=0
    shift[1] += MPI_buffer_.part_index_recv_sz[0][0];//Particles coming from xmin all go to bin 0 and shift all the other bins.
    shift[particles->last_index.size()] += MPI_buffer_.part_index_recv_sz[0][1];//Used only to count the total number of particles arrived.
    //idim>0
    for( idim = 1; idim < ndim; idim++ ) {
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            n_part_recv = MPI_buffer_.part_index_recv_sz[idim][iNeighbor];
            for( unsigned int j=0; j<( unsigned int )n_part_recv ; j++ ) {
                //We first evaluate how many particles arrive in each bin.
                ii = int( ( MPI_buffer_.partRecv[idim][iNeighbor].position( 0, j )-min_loc )/dbin ); //bin in which the particle goes.
                shift[ii+1]++; // It makes the next bins shift.
            }
        }
    }


    //Must be done sequentially
    for( unsigned int j=1; j<particles->last_index.size()+1; j++ ) { //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
        shift[j]+=shift[j-1];
    }
    //Make room for new particles
    if( shift[particles->last_index.size()] ) {
        //! vecor::resize of Charge crashed ! Temporay solution : push_back / Particle
        //particles->initialize( particles->size()+shift[particles->last_index.size()], particles->Position.size() );
        for( int inewpart=0 ; inewpart<shift[particles->last_index.size()] ; inewpart++ ) {
            particles->createParticle();
        }
    }

    //Shift bins, must be done sequentially
    for( unsigned int j=particles->last_index.size()-1; j>=1; j-- ) {
        int n_particles = particles->last_index[j]-particles->first_index[j]; //Nbr of particle in this bin
        nmove = min( n_particles, shift[j] ); //Nbr of particles to move
        lmove = max( n_particles, shift[j] ); //How far particles must be shifted
        if( nmove>0 ) {
            particles->overwriteParticle( particles->first_index[j], particles->first_index[j]+lmove, nmove );
        }
        particles->first_index[j] += shift[j];
        particles->last_index[j] += shift[j];
    }

    //Space has been made now to write the arriving particles into the correct bins
    //idim == 0  is the easy case, when particles arrive either in first or last bin.
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        n_part_recv = MPI_buffer_.part_index_recv_sz[0][iNeighbor];
        //if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
        if( ( n_part_recv!=0 ) ) {
            ii = iNeighbor*( particles->last_index.size()-1 ); //0 if iNeighbor=0(particles coming from Xmin) and particles->last_index.size()-1 otherwise.
            MPI_buffer_.partRecv[0][iNeighbor].overwriteParticle( 0, *particles, particles->last_index[ii], n_part_recv );
            particles->last_index[ii] += n_part_recv ;
        }
    }
    //idim > 0; this is the difficult case, when particles can arrive in any bin.
    for( idim = 1; idim < ndim; idim++ ) {
        //if (idim!=iDim) continue;
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            n_part_recv = MPI_buffer_.part_index_recv_sz[idim][iNeighbor];
            //if ( (neighbor_[idim][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
            if( ( n_part_recv!=0 ) ) {
                for( unsigned int j=0; j<( unsigned int )n_part_recv; j++ ) {
                    ii = int( ( MPI_buffer_.partRecv[idim][iNeighbor].position( 0, j )-min_loc )/dbin ); //bin in which the particle goes.
                    MPI_buffer_.partRecv[idim][iNeighbor].overwriteParticle( j, *particles, particles->last_index[ii] );
                    particles->last_index[ii] ++ ;
                }
            }
        }
    }


    //The width of one bin is cell_length[0] * clrw.

    int p1, p2, first_index_init;
    unsigned int bin;
    double limit;


    //Backward pass
    for( bin=0; bin<particles->first_index.size()-1; bin++ ) { //Loop on the bins.
        limit = min_loc + ( bin+1 )*cell_length[0]*clrw;
        p1 = particles->last_index[bin]-1;
        //If first particles change bin, they do not need to be swapped.
        while( p1 == particles->last_index[bin]-1 && p1 >= particles->first_index[bin] ) {
            if( particles->position( 0, p1 ) >= limit ) {
                particles->last_index[bin]--;
            }
            p1--;
        }
        //         Now particles have to be swapped
        for( p2 = p1 ; p2 >= particles->first_index[bin] ; p2-- ) { //Loop on the bin's particles.
            if( particles->position( 0, p2 ) >= limit ) {
                //This particle goes up one bin.
                particles->swapParticle( p2, particles->last_index[bin]-1 );
                particles->last_index[bin]--;
            }
        }
    }
    //Forward pass + Rebracketting
    for( bin=1; bin<particles->first_index.size(); bin++ ) { //Loop on the bins.
        limit = min_loc + bin*cell_length[0]*clrw;
        first_index_init = particles->first_index[bin];
        p1 = particles->first_index[bin];
        while( p1 == particles->first_index[bin] && p1 < particles->last_index[bin] ) {
            if( particles->position( 0, p1 ) < limit ) {
                particles->first_index[bin]++;
            }
            p1++;
        }
        for( p2 = p1 ; p2 < particles->last_index[bin] ; p2++ ) { //Loop on the bin's particles.
            if( particles->position( 0, p2 ) < limit ) {
                //This particle goes down one bin.
                particles->swapParticle( p2, particles->first_index[bin] );
                particles->first_index[bin]++;
            }
        }

        //Rebracketting
        //Number of particles from bin going down is: particles->first_index[bin]-first_index_init.
        //Number of particles from bin-1 going up is: first_index_init-particles->last_index[bin-1].
        //Total number of particles we need to swap is the min of both.
        p2 = min( particles->first_index[bin]-first_index_init, first_index_init-particles->last_index[bin-1] );
        if( p2 >0 ) {
            particles->swapParticle( particles->last_index[bin-1], particles->first_index[bin]-p2, p2 );
        }
        particles->last_index[bin-1] += particles->first_index[bin] - first_index_init;
        particles->first_index[bin] = particles->last_index[bin-1];
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function configures the type of species according to the default mode
//! regardless the number of particles per cell
// ---------------------------------------------------------------------------------------------------------------------
void Species::defaultConfigure( Params &param, Patch *patch )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function configures the species according to the vectorization mode
// ---------------------------------------------------------------------------------------------------------------------
void Species::configuration( Params &param, Patch *patch )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function reconfigures the species operators after evaluating
//! the best mode from the particle distribution
// ---------------------------------------------------------------------------------------------------------------------
void Species::reconfiguration( Params &param, Patch *patch )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Sort particles using the count sort method
// ---------------------------------------------------------------------------------------------------------------------
void Species::countSortParticles( Params &params )
{
    unsigned int ip, npart, ixy, tot, oc, nxy, token;
    int ix, iy;
    double x, y;

    nxy = params.n_space[0]*params.n_space[1];
    token = ( particles == &particles_sorted[0] );

    int indices[nxy];
    npart = particles->size();
    //particles_sorted = particles ;
    particles_sorted[token].initialize( npart, *particles );

    for( unsigned int i=0; i < nxy ; i++ ) {
        indices[i] = 0 ;
    }

    // first loop counts the # of particles in each cell
    for( ip=0; ip < npart; ip++ ) {
        x = particles->position( 0, ip )-min_loc;
        y = particles->position( 1, ip )-min_loc_vec[1];

        ix = floor( x * dx_inv_[0] ) ;
        iy = floor( y * dx_inv_[1] ) ;

        ixy = iy + ix*params.n_space[1];


        indices[ixy] ++;
    }

    // second loop convert the count array in cumulative sum
    tot=0;
    for( ixy=0; ixy < nxy; ixy++ ) {
        oc = indices[ixy];
        indices[ixy] = tot;
        tot += oc;
    }
    
    // last loop puts the particles and update the count array
    for( ip=0; ip < npart; ip++ ) {
        x = particles->position( 0, ip )-min_loc;
        y = particles->position( 1, ip )-min_loc_vec[1];

        ix = floor( x * dx_inv_[1] ) ;
        iy = floor( y * dx_inv_[2] ) ;

        ixy = iy + ix*params.n_space[1];
        particles->overwriteParticle( ip, particles_sorted[token], indices[ixy] );
        indices[ixy]++;
    }

    particles = &particles_sorted[token] ;

}

// Move all particles from another species to this one
void Species::importParticles( Params &params, Patch *patch, Particles &source_particles, vector<Diagnostic *> &localDiags )
{
    unsigned int npart = source_particles.size(), nbin=particles->first_index.size();
    double inv_cell_length = 1./ params.cell_length[0];

    // If this species is tracked, set the particle IDs
    if( particles->tracked ) {
        dynamic_cast<DiagnosticTrack *>( localDiags[tracking_diagnostic] )->setIDs( source_particles );
    }

    // Move particles
    vector<int> src_bin_keys( npart, 0 );
    for( unsigned int i=0; i<npart; i++ ) {
        // Copy particle to the correct bin
        src_bin_keys[i] = source_particles.position( 0, i )*inv_cell_length - ( patch->getCellStartingGlobalIndex( 0 ) + params.oversize[0] );
        src_bin_keys[i] /= params.clrw;
    }

    vector<int> bin_count( nbin, 0 );
    for( unsigned int ip=0; ip < npart ; ip++ )
        bin_count[src_bin_keys[ip]] ++;

    // sort new parts par bins
    int istart = 0;
    int istop  = bin_count[0];

    for ( int ibin = 0 ; ibin < (int)nbin ; ibin++ ) {
        if (bin_count[ibin]!=0) {
            for( int ip=istart; ip < istop ; ip++ ) {
                if ( src_bin_keys[ip] == ibin )
                    continue;
                else { // rearrange particles
                    int ip_swap = istop;
                    while (( src_bin_keys[ip_swap] != ibin ) && (ip_swap<(int)npart))
                        ip_swap++;
                    source_particles.swapParticle(ip, ip_swap);
                    int tmp = src_bin_keys[ip];
                    src_bin_keys[ip] = src_bin_keys[ip_swap];
                    src_bin_keys[ip_swap] = tmp;
                } // rearrange particles
            } // end loop on particles of a cell

            // inject in main data structure per cell
            source_particles.copyParticles( istart, bin_count[ibin],
                                        *particles,
                                        particles->first_index[ibin] );
            particles->last_index[ibin] += bin_count[ibin];
            for ( unsigned int idx=ibin+1 ; idx<particles->last_index.size() ; idx++ ) {
                particles->first_index[idx] += bin_count[ibin];
                particles->last_index[idx]  += bin_count[ibin];
            }

        }
        // update istart/istop fot the next cell
        istart += bin_count[ibin];
        if ( ibin != (int)nbin-1  )
            istop  += bin_count[ibin+1];
        else
            istop = npart;

    } // End cell loop

    source_particles.clear();
}


// ------------------------------------------------
// Set position when using restart & moving window
// patch are initialized with t0 position
// ------------------------------------------------
//void Species::updateMvWinLimits(double x_moved)
//{
//    partBoundCond->updateMvWinLimits(x_moved);
//    min_loc += x_moved;
//
//} // End updateMvWinLimits


//Do we have to project this species ?
bool Species::isProj( double time_dual, SimWindow *simWindow )
{

    return time_dual > time_frozen_  || ( simWindow->isMoving( time_dual ) || Ionize ) ;

    //Recompute frozen particles density if
    //moving window is activated, actually moving at this time step, and we are not in a density slope.
    /*    bool isproj =(time_dual > species_param.time_frozen_  ||
                 (simWindow && simWindow->isMoving(time_dual) &&
                     (species_param.species_geometry == "gaussian" ||
                         (species_param.species_geometry == "trapezoidal" &&
                            //Before end of density ramp up.
                            (simWindow->getXmoved() < species_param.vacuum_length[0] + species_param.dens_length_x[1] + oversize[0]*cell_length[0] ||
                            //After begining of density ramp down.
                            simWindow->getXmoved() +  simWindow->getNspace_win_x()*cell_length[0] > species_param.vacuum_length[0] + species_param.dens_length_x[1]+ species_param.dens_length_x[0]
                            )
                        )
                    )
                )
            );
            return isproj;*/
    //return time_dual > species_param.time_frozen_  || (simWindow && simWindow->isMoving(time_dual)) ;
}

void Species::disableXmax()
{
    partBoundCond->bc_xmax   = NULL;
}

void Species::setXminBoundaryCondition()
{
    partBoundCond->bc_xmin   = &remove_particle;
}

// ---------------------------------------------------------------------------------------------------------------------
// Particle merging cell by cell
// ---------------------------------------------------------------------------------------------------------------------
void Species::mergeParticles( double time_dual, unsigned int ispec,
                              Params &params,
                              Patch *patch, SmileiMPI *smpi,
                              std::vector<Diagnostic *> &localDiags ) {}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
//   - calculate the new momentum
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotiveUpdateSusceptibilityAndMomentum( double time_dual, unsigned int ispec,
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
    if( time_dual>time_frozen_ || Ionize) { // moving particle

        smpi->dynamics_resize( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) { // loop on ibin

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[7] += MPI_Wtime() - timer;
#endif

            // Ionization
            if( Ionize ) {
            
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
                vector<double> *EnvEabs_part = &( smpi->dynamics_EnvEabs_part[ithread] );
                vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );
                vector<double> *Phipart = &( smpi->dynamics_PHIpart[ithread] );
                Interp->envelopeFieldForIonization( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
                Ionize->envelopeIonization( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, EnvEabs_part, EnvExabs_part, Phipart, patch, Proj );
                
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif            
            }
            
            if( time_dual<=time_frozen_ ) continue; // Do not push nor project frozen particles

            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[8] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle momenta
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[9] += MPI_Wtime() - timer;
#endif

        } // end loop on ibin
    } else { // immobile particle
    } //END if time vs. time_frozen_
} // ponderomotiveUpdateSusceptibilityAndMomentum

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotiveProjectSusceptibility( double time_dual, unsigned int ispec,
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

        smpi->dynamics_resize( ithread, nDim_particle, particles->last_index.back(), false );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) { // loop on ibin

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[7] += MPI_Wtime() - timer;
#endif

            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[8] += MPI_Wtime() - timer;
#endif


        } // end loop on ibin
    } else { // immobile particle
    } //END if time vs. time_frozen_
} // ponderomotiveProjectSusceptibility


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the ponderomotive potential and its gradient at the particle position, for present and previous timestep
//   - calculate the new particle position
//   - particles BC
//   - project charge and current density
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
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

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {

            // Interpolate the ponderomotive potential and its gradient at the particle position, present and previous timestep
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[10] += MPI_Wtime() - timer;
#endif

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[11] += MPI_Wtime() - timer;
#endif

            // Apply wall and boundary conditions
            if( mass_>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                            nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                    if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                    }
                }

            } else if( mass_==0 ) {
                ERROR( "Particles with zero mass cannot interact with envelope" );
                // for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                //     for (iPart=particles->first_index[ibin] ; (int)iPart<particles->last_index[ibin]; iPart++ ) {
                //         double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                //         if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, dtgf, ener_iPart)) {
                //                 nrj_lost_per_thd[tid] += ener_iPart;
                //         }
                //     }
                // }
                //
                // // Boundary Condition may be physical or due to domain decomposition
                // // apply returns 0 if iPart is not in the local domain anymore
                // //        if omp, create a list per thread
                // for (iPart=particles->first_index[ibin] ; (int)iPart<particles->last_index[ibin]; iPart++ ) {
                //     if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                //         addPartInExchList( iPart );
                //         nrj_lost_per_thd[tid] += ener_iPart;
                //     }
                //  }

            } // end mass_ = 0? condition

            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[12] += MPI_Wtime() - timer;
#endif

        } // end ibin loop

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

    } // end case of moving particle
    else { // immobile particle

        if( diag_flag &&( !particles->is_test ) ) {
            double *b_rho=nullptr;
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj

                if( params.geometry != "AMcylindrical" ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                    for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basic( b_rho, ( *particles ), iPart, 0 );
                    } 
                } else {
                    int n_species = patch->vecSpecies.size();
                    complex<double> *b_rho=nullptr;
                    ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );

                    for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                        int ifield = imode*n_species+ispec;
                        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                            b_rho = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                            for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                                Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                            } // end loop on particles
                        } //end loop for projection on buffer_proj

                    } // end loop on modes
                } // end if on geometry

            }//End loop on bins
        } // end condition on diag and not particle test
    }//END if time vs. time_frozen_
} // End ponderomotive_position_update

void Species::check( Patch *patch, std::string title )
{
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;
    double sum_px = 0;
    double sum_py = 0;
    double sum_pz = 0;
    double sum_w = 0;
    unsigned int sum_ck = 0;
    for( unsigned int ip=0; ip < particles->size() ; ip++ ) {
        sum_x += particles->position( 0, ip );
        sum_y += particles->position( 1, ip );
        sum_z += particles->position( 2, ip );
        sum_px += particles->momentum( 0, ip );
        sum_py += particles->momentum( 1, ip );
        sum_pz += particles->momentum( 1, ip );
        sum_w += particles->weight( ip );
        sum_ck += particles->cell_keys[ip];
    }
    std::cerr << "Check sum at " << title
              << " for "<< this->name_
              << " in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ") "
              << " mpi process " << patch->MPI_me_
              << " - mode: " << this->vectorized_operators
              << " - nb bin: " << particles->first_index.size()
              << " - nbp: " << particles->size()
              << " - w: " << sum_w
              << " - x: " << sum_x
              << " - y: " << sum_y
              << " - z: " << sum_z
              << " - px: " << sum_px
              << " - py: " << sum_py
              << " - pz: " << sum_pz
              << " - ck: " << sum_ck
              << '\n';
}

void Species::eraseWeightlessParticles()
{
    unsigned int nbins = particles->first_index.size();
    unsigned int i = 0, available_i = 0;
    
    // Loop all particles, bin per bin
    // Overwrite over earlier particles to erase them
    for( unsigned int ibin = 0; ibin < nbins; ibin++ ) {
        particles->first_index[ibin] = available_i;
        while( i < (unsigned int) particles->last_index[ibin] ) {
            if( particles->weight(i) > 0. ) {
                if( i > available_i ) {
                    particles->overwriteParticle( i, available_i );
                }
                available_i ++;
            }
            i++;
        }
        particles->last_index[ibin] = available_i;
    }
    
    // Remove trailing particles
    particles->eraseParticleTrail( available_i );
}
