
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

// necessary for the static_cast
#include "ProjectorAM2Order.h"

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
    file_position_npart_( 0 ),
    file_momentum_npart_( 0 ),
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

    particles_to_move = new Particles();

    // this variable is used for tasks
    geometry = params.geometry;

}//END Species creator

void Species::initCluster( Params &params )
{
    // Arrays of the min and max indices of the particle bins
    particles->first_index.resize( params.n_space[0]/clrw );
    particles->last_index.resize( params.n_space[0]/clrw );
    Nbins = params.n_space[0]/clrw ;

    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)

    //Primal dimension of fields.
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;

    //Dual dimension of fields.
    f_dim0_d =  params.n_space[0] + 2 * oversize[0] +2;
    f_dim1_d =  params.n_space[1] + 2 * oversize[1] +2;
    f_dim2_d =  params.n_space[2] + 2 * oversize[2] +2;

    b_dim.resize( params.nDim_field, 1 );
    if( nDim_particle == 1 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0];
        b_dim[1] =  1;
        f_dim1 = 1;
        f_dim2 = 1;
        f_dim1_d = 1;
        f_dim2_d = 1;
    }
    if( nDim_particle == 2 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] =  f_dim1;
        f_dim2 = 1;
        f_dim2_d = 1;
    }
    if( nDim_particle == 3 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] = f_dim1;
        b_dim[2] = f_dim2;
    }

#ifdef _OMPTASKS
        nrj_lost_per_bin                       = new double[Nbins];
        nrj_radiation_per_bin                  = new double[Nbins];      
        if (params.geometry != "AMcylindrical" ){
            //! buffers for currents and charge
            b_Jx.resize(Nbins);
            b_Jy.resize(Nbins);
            b_Jz.resize(Nbins);
            b_rho.resize(Nbins);
            if (params.Laser_Envelope_model){ 
                b_Chi.resize(Nbins);
            }

            size_proj_buffer_rho = b_dim[0]*b_dim[1]*f_dim2;
            size_proj_buffer_Jx  = b_dim[0]*b_dim[1]*f_dim2;
            size_proj_buffer_Jy  = b_dim[0]*f_dim1_d*f_dim2;
            size_proj_buffer_Jz  = b_dim[0]*b_dim[1]*f_dim2_d;

            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                // allocate current-buffers, then put to zero their content
                b_Jx[ibin]  = new double[size_proj_buffer_Jx ];
                b_Jy[ibin]  = new double[size_proj_buffer_Jy ];
                b_Jz[ibin]  = new double[size_proj_buffer_Jz ];
                b_rho[ibin] = new double[size_proj_buffer_rho];
                if (params.Laser_Envelope_model){ // Chi has the same size of rho
                    b_Chi[ibin] = new double[size_proj_buffer_rho];
                }
            }
        } else { // AM geometry
            //! buffers for currents and charge
            size_proj_buffer_rhoAM = b_dim[0]*b_dim[1]    * params.nmodes; // used for rhoAM
            size_proj_buffer_Jl    = b_dim[0]*b_dim[1]    * params.nmodes; // used for Jl
            size_proj_buffer_Jr    = b_dim[0]*f_dim1_d    * params.nmodes; // used for Jr
            size_proj_buffer_Jt    = b_dim[0]*b_dim[1]    * params.nmodes; // used for Jt

            //! buffers for currents and charge
            b_Jl.resize(Nbins);
            b_Jr.resize(Nbins);
            b_Jt.resize(Nbins);
            b_rhoAM.resize(Nbins);
            if (params.Laser_Envelope_model){
                b_ChiAM.resize(Nbins);
            }

            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                // allocate current-buffers, then put to zero their content
                b_Jl[ibin]    = new std::complex<double>[size_proj_buffer_Jl   ];
                b_Jr[ibin]    = new std::complex<double>[size_proj_buffer_Jr   ];
                b_Jt[ibin]    = new std::complex<double>[size_proj_buffer_Jt   ];
                b_rhoAM[ibin] = new std::complex<double>[size_proj_buffer_rhoAM];
                if (params.Laser_Envelope_model){ // Chi has the same size of rho
                    b_ChiAM[ibin] = new double[size_proj_buffer_rhoAM];
                }
            }
        } // end condition on geometry
#endif

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


// Create the particles once the namelist is read
void Species::initParticles( Params &params, Patch *patch, bool with_particles, Particles * like_particles )
{
    if( params.restart || !with_particles ) {
        
        if( like_particles ) {
            particles->initialize( 0, *like_particles );
        } else {
            particles->initialize( 0, params.nDim_particle, params.keep_position_old );
        }
        
    } else {
        
        // Area for particle creation
        struct SubSpace init_space;
        init_space.cell_index_[0] = 0;
        init_space.cell_index_[1] = 0;
        init_space.cell_index_[2] = 0;
        init_space.box_size_[0]   = params.n_space[0];
        init_space.box_size_[1]   = params.n_space[1];
        init_space.box_size_[2]   = params.n_space[2];
        
        // Creation of the particle creator and association to the new species
        ParticleCreator particle_creator;
        particle_creator.associate(this);
        particle_creator.create( init_space, params, patch, 0 );
        
    }
    
}

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

    particles_to_move->initialize( 0, *particles );

}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete particles_to_move;

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

#ifdef _OMPTASKS
    if (nrj_lost_per_bin != NULL){
        delete nrj_lost_per_bin;
    }
    if (nrj_radiation_per_bin != NULL){
        delete nrj_radiation_per_bin;
    }     
    if (geometry != "AMcylindrical"){
        if (b_Jx[0]){
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                // delete buffers
                delete[] b_Jx[ibin];
                delete[] b_Jy[ibin];
                delete[] b_Jz[ibin];
                delete[] b_rho[ibin];
            }
        }
        if (Push_ponderomotive_position){
            if (b_Chi[0]){
                for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                    // delete buffers
                    delete[] b_Chi[ibin];
                }
            }
        }
    } else {
        if (b_Jl[0]){
            for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                // delete buffer
                delete[] b_Jl[ibin];
                delete[] b_Jr[ibin];
                delete[] b_Jt[ibin];
                delete[] b_rhoAM[ibin];
            }
        }
        if (Push_ponderomotive_position){
            if (b_ChiAM[0]){
                for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
                    // delete buffer
                    delete[] b_ChiAM[ibin];
                }
            }
        }
    }
#endif
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
            patch->patch_timers_[0] += MPI_Wtime() - timer;
#endif

            // Ionization
            if( Ionize ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                ( *Ionize )( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, patch, Proj );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4] += MPI_Wtime() - timer;
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
                patch->patch_timers_[5] += MPI_Wtime() - timer;
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
                patch->patch_timers_[6] += MPI_Wtime() - timer;
#endif

            }

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Push the particles and the photons
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
            //particles->testMove( particles->first_index[ibin], particles->last_index[ibin], params );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[1] += MPI_Wtime() - timer;
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
                    if(!params.is_spectral){
                        partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                        nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                    }

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
                patch->patch_timers_[3] += MPI_Wtime() - timer;
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
                patch->patch_timers_[2] += MPI_Wtime() - timer;
#endif
                if(params.is_spectral && mass_>0){
                    partBoundCond->apply( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], this, ithread, ener_iPart );
                    nrj_lost_per_thd[tid] += mass_ * ener_iPart;
                }

            }// ibin
        } // end if moving particle

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

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

void Species::dynamicsTasks( double time_dual, unsigned int ispec,
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

    int *bin_can_radiate = new int[Nbins];
    int *bin_can_push = new int[Nbins];

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
                Ionize->ionizationTunnelWithTasks( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, patch, Proj, ibin, ibin*clrw, bJx, bJy, bJz );

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
                } // end if AM
            } // end condition on test and mass

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[2*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            }//end task for Proj of ibin
         }// end ibin loop for Proj

        // reduction of the lost energy in each ibin 
        // the dependency ensures that it is done after the particles BC
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

        // sum the radiated energy / energy converted in pairs
        // The dependencies above ensure that this is done after the Radiation and MultiPhoton Breit Wheeler methods
        if( Radiate || Multiphoton_Breit_Wheeler_process) {
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
        } // end if Radiate or Multiphoton_Breit_Wheeler_process
        } // end task for lost/radiated energy reduction          

     } else { // immobile particle
         if( diag_flag &&( !particles->is_test ) ) {
             if (Ionize) { // a dependency is needed to project after ionization
                 if( params.geometry != "AMcylindrical" ) {

                    for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                        #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer) depend(in:bin_has_ionized[ibin])
#else
                        #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_ionized[ibin])
#endif
                        {
#ifdef  __DETAILED_TIMERS
                        ithread = omp_get_thread_num();
                        timer = MPI_Wtime();
#endif
                        for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
                        for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                            Proj->basic( b_rho[ibin], ( *particles ), iPart, 0, ibin*clrw );
                        }
#ifdef  __DETAILED_TIMERS
                        patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                        } // end task projection for frozen or test
                    } // end ibin
                } else { // AM geometry

                    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                        #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer) depend(in:bin_has_ionized[ibin])
#else
                        #pragma omp task default(shared) firstprivate(ibin,bin_size0) depend(in:bin_has_ionized[ibin]) 
#endif
                        {
#ifdef  __DETAILED_TIMERS
                        ithread = omp_get_thread_num();
                        timer = MPI_Wtime();
#endif
                    
                        for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
                        for( int imode = 0; imode<params.nmodes; imode++ ) {
                            for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                                Proj->basicForComplexOnBuffer( b_rhoAM[ibin], ( *particles ), iPart, 0, imode, bin_size0, ibin*clrw );
                            } // end loop on particles
                        } // end imode loop
#ifdef  __DETAILED_TIMERS
                        patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                        } // end task projection for frozen or test
                    } //end ibin
                } // end if on geometry

            } else { // no ionization, no dependency for tasks needed

                 if( params.geometry != "AMcylindrical" ) {

                    for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                        #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer) 
#else
                        #pragma omp task default(shared) firstprivate(ibin)
#endif
                        {
#ifdef  __DETAILED_TIMERS
                        ithread = omp_get_thread_num();
                        timer = MPI_Wtime();
#endif
                        for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
                        for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                            Proj->basic( b_rho[ibin], ( *particles ), iPart, 0, ibin*clrw );
                        }
#ifdef  __DETAILED_TIMERS
                        patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                        } // end task projection for frozen or test
                    } // end ibin
                } else { // AM geometry

                    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                        #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer) 
#else
                        #pragma omp task default(shared) firstprivate(ibin,bin_size0)
#endif
                        {
#ifdef  __DETAILED_TIMERS
                        ithread = omp_get_thread_num();
                        timer = MPI_Wtime();
#endif
                    
                        for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
                        for( int imode = 0; imode<params.nmodes; imode++ ) {
                            for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                                Proj->basicForComplexOnBuffer( b_rhoAM[ibin], ( *particles ), iPart, 0, imode, bin_size0, ibin*clrw );
                            } // end loop on particles
                        } // end imode loop
#ifdef  __DETAILED_TIMERS
                        patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                        } // end task projection for frozen or test
                    } //end ibin
                } // end if on geometry

            } // end condition on ionization
        } // end condition on diag and not particle test

     } // end moving particle
     

     }// end taskgroup for all the Interp, Push, Particles BC and Projector tasks
  

//      // smpi->reduce_dynamics_buffer_size( buffer_id, params.geometry=="AMcylindrical" );

} // end dynamicsTasks



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
void Species::computeCharge( unsigned int ispec, ElectroMagn *EMfields, bool old /*=false*/ )
{
    // -------------------------------
    // calculate the particle charge
    // -------------------------------
    if( ( !particles->is_test ) ) {
        if( !dynamic_cast<ElectroMagnAM *>( EMfields ) ) {
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                double *b_rho = &( *EMfields->rho_ )( 0 );

                for( unsigned int iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            }
        } else {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            unsigned int Nmode = emAM->rho_AM_.size();
            for( unsigned int imode=0; imode<Nmode; imode++ ) {
                unsigned int ifield = imode*(*EMfields).n_species+ispec;
                complex<double> *b_rho = old ? &( *emAM->rho_old_AM_[imode] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 );
                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                    for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }
        }
    }
}//END computeCharge


void Species::extractParticles()
{
    //particles->last_index[0] -= dev_particles->extractParticles( particles_to_move );

    // Count cell_keys = -1
    //count( particles->cell_keys.begin(), particles->cell_keys.end(), -1 );
    //nparts_to_move_ = thrust::count(thrust::device, nvidia_cell_keys.begin(), nvidia_cell_keys.begin()+nparts, -1);
    // resize particles_to_move
    // ....resize( nparts_to_move )
    // copy in particles_to_move if cell_keys = -1
    //thrust::copy_if(thrust::device, iter, iter+nparts, nvidia_cell_keys.begin(), iter_copy, count_if_out());

    particles_to_move->clear();
    for ( int ipart=0 ; ipart<(int)(getNbrOfParticles()) ; ipart++ ) {
        if ( particles->cell_keys[ipart] == -1 ) {
            particles->copyParticle( ipart, *particles_to_move );
        }
    }

}

void Species::injectParticles( Params &params )
{
    // through particles_to_move ... not in scalar but :

    //thrust::remove_if( ... remove_if_out() ); cell_keys < 0
    // resize particles (include in sortParticles)
    //thrust::copy_n(thrust::device, iter_copy, nparts_add, iter+nparts); (include in sortParticles)

}


// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::sortParticles( Params &params, Patch * patch )
{
    injectParticles( params );

    int ndim = params.nDim_field;
    int idim;

    int total_number_part_recv = 0;
    //Merge all MPI_buffer_.partRecv in particles_to_move
    for( int idim = 0; idim < ndim; idim++ ) {
        for( int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++ ) {
            int n_part_recv = MPI_buffer_.part_index_recv_sz[idim][iNeighbor];
            if( ( n_part_recv!=0 ) ) {
                 // insert n_part_recv in particles_to_move from 0
                //MPI_buffer_.partRecv[idim][iNeighbor].copyParticles( 0, n_part_recv, *particles_to_move, 0 );
                total_number_part_recv += n_part_recv;
                //particles->last_index[particles->last_index.size()-1] += n_part_recv;
                //particles->cell_keys.resize(particles->cell_keys.size()+n_part_recv);
            }
        }
    }
    //cout << "\t Species id : " << species_number_ << " - nparticles recv : " << blabla << endl;


    // Sort to adapt do cell_keys usage
    std::vector<int> indexes_of_particles_to_exchange;
    for ( int ipart=0 ; ipart< (int)(getNbrOfParticles()) ; ipart++ ) {
        if ( particles->cell_keys[ipart] == -1 ) {
            indexes_of_particles_to_exchange.push_back( ipart );
        }
    }
    //cout << "\t Species id : " << species_number_ << " - nparticles send : " << indexes_of_particles_to_exchange.size() << endl;

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

    //particles->cell_keys.resize( particles->size() );
    particles->resizeCellKeys(particles->size());
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
    //particles->cell_keys.resize( particles->size() );
    particles->resizeCellKeys( particles->size() );

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
    partBoundCond->bc_xmax   = &internal_sup;
}

void Species::setXminBoundaryCondition()
{
    partBoundCond->bc_xmin   = &remove_particle_inf;
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
            patch->patch_timers_[7] += MPI_Wtime() - timer;
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
                patch->patch_timers_[4] += MPI_Wtime() - timer;
#endif
            }
            
            if( time_dual<=time_frozen_ ) continue; // Do not push nor project frozen particles

            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[8] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle momenta
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[9] += MPI_Wtime() - timer;
#endif

        } // end loop on ibin
    } else { // immobile particle
    } //END if time vs. time_frozen_
} // ponderomotiveUpdateSusceptibilityAndMomentum

void Species::ponderomotiveUpdateSusceptibilityAndMomentumTasks( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags, int buffer_id )
{
#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif
    int bin_size0 = b_dim[0];
    // for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
    //     nrj_lost_per_bin[ibin] = 0.;
    //     nrj_radiation_per_bin[ibin] = 0.;
    // }
    // Init tags for the task dependencies of the particle operations
    int *bin_has_interpolated                   = new int[Nbins]; // the last element is used to manage the Multiphoton Breit Wheeler dependency
    int *bin_has_ionized                        = new int[Nbins];
    int *bin_has_projected_chi                  = new int[Nbins];

    int *bin_can_project_chi = new int[Nbins];


    if (Ionize){ 
        // ... project chi after ionization if present
        bin_can_project_chi = bin_has_ionized;
    } else { 
        // ... project chi directly after interpolation if ionization is not present
        bin_can_project_chi = bin_has_interpolated;       
    }


    #pragma omp taskgroup
    {

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize) { // moving particle

        smpi->dynamics_resize( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );
        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) { // loop on ibin
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
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif

            // Interpolate the fields and envelope at the particle position
            Interp->fieldsAndEnvelopeForTasks( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), buffer_id );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[7*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task interp
        } // end ibin

        // Ionization
        if( Ionize ) {
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) { // loop on ibin
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
                vector<double> *EnvEabs_part = &( smpi->dynamics_EnvEabs_part[buffer_id] );
                vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[buffer_id] );
                vector<double> *Phipart = &( smpi->dynamics_PHIpart[buffer_id] );
                Interp->envelopeFieldForIonizationTasks( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), buffer_id );
                Ionize->envelopeIonization( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, EnvEabs_part, EnvExabs_part, Phipart, patch, Proj, ibin, 0 );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                } // end task ionize
            } // end ibin
        }  // end Ionize      
      
        if( time_dual>time_frozen_) {
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) { // loop on ibin
#ifdef  __DETAILED_TIMERS
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_can_project_chi[ibin]) depend(out:bin_has_projected_chi[ibin]) private(ithread,timer)
#else
                #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_can_project_chi[ibin]) depend(out:bin_has_projected_chi[ibin])
#endif
                {
                // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
                ithread = omp_get_thread_num();
                timer = MPI_Wtime();
#endif
                
                if (params.geometry != "AMcylindrical"){
                    Proj->susceptibilityOnBuffer( EMfields, b_Chi[ibin], 
                                                      ibin*clrw, *particles, mass_, smpi, 
                                                      particles->first_index[ibin], particles->last_index[ibin], 
                                                      buffer_id );
                } else {
                    ProjectorAM2Order *ProjAM = static_cast<ProjectorAM2Order *>(Proj);
                    ProjAM->susceptibilityOnAMBuffer( EMfields, b_ChiAM[ibin], 
                                                      ibin*clrw, bin_size0,
                                                      *particles, mass_, smpi, 
                                                      particles->first_index[ibin], particles->last_index[ibin], 
                                                      buffer_id );
                } // end if geometry
#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[8*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                } // end task susceptibility
            } // end ibin
        }

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) { // loop on ibin
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_projected_chi[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_projected_chi[ibin]) 
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif
            // Push only the particle momenta
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], buffer_id);
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[9*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task susceptibility
        } // end ibin
    
    } else { // immobile particle
    } //END if time vs. time_frozen_

    } // end taskgroup
 

} // ponderomotiveUpdateSusceptibilityAndMomentumTasks

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
            patch->patch_timers_[7] += MPI_Wtime() - timer;
#endif

            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[8] += MPI_Wtime() - timer;
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

    int tid( 0 );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle updated position
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->dynamics_resize( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
            double ener_iPart( 0. );

            // Interpolate the ponderomotive potential and its gradient at the particle position, present and previous timestep
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[10] += MPI_Wtime() - timer;
#endif

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[11] += MPI_Wtime() - timer;
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
                ERROR( "Particles with zero mass cannot interact with envelope" );

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
            patch->patch_timers_[12] += MPI_Wtime() - timer;
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

void Species::ponderomotiveUpdatePositionAndCurrentsTasks( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag, PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags, int buffer_id )
{
#ifdef  __DETAILED_TIMERS
    double timer;
    int ithread;
#endif
    int bin_size0 = b_dim[0];

    for( unsigned int ibin = 0 ; ibin < Nbins ; ibin++ ) {
        nrj_lost_per_bin[ibin] = 0.;
        // nrj_radiation_per_bin[ibin] = 0.;
    }
    // Init tags for the task dependencies of the particle operations
    int *bin_has_interpolated               = new int[Nbins]; // the last element is used to manage the Multiphoton Breit Wheeler dependency
    int *bin_has_pushed                     = new int[Nbins];
    int *bin_has_done_particles_BC          = new int[Nbins];
    int *bin_has_projected                  = new int[Nbins];

    #pragma omp taskgroup
    {

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ ) {

        smpi->dynamics_resize( buffer_id, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
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
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif
            Interp->timeCenteredEnvelopeForTasks( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), buffer_id );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[10*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task interpolate
        } // end ibin

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_pushed[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_interpolated[ibin]) depend(out:bin_has_pushed[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], buffer_id );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[11*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task
        } // end ibin

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_pushed[ibin]) depend(out:bin_has_done_particles_BC[ibin])
#endif
            {
#ifdef  __DETAILED_TIMERS
            ithread = omp_get_thread_num();
            timer = MPI_Wtime();
#endif
            // Apply wall and boundary conditions
            if( mass_>0 ) {
                double ener_iPart;
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
                ERROR( "Particles with zero mass cannot interact with envelope" );
            
            } // end mass_ = 0? condition
#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task
        } // end ibin

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
#ifdef  __DETAILED_TIMERS
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin]) private(ithread,timer)
#else
            #pragma omp task default(shared) firstprivate(ibin) depend(in:bin_has_done_particles_BC[ibin]) depend(out:bin_has_projected[ibin])
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
                  } // end if AM
            } // end condition on test and mass

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
            } // end task
        } // end ibin

        // reduction of the lost energy in each ibin 
        // the dependency ensures that it is done after the particles BC
        //#pragma omp task default(shared) depend(in:bin_has_done_particles_BC[0:(Nbins-1)])
        // using depend out on particles BC a segfault is caused - check why this happens
        #pragma omp task default(shared) depend(in:bin_has_projected[0:(Nbins-1)])
        {
        // reduce the energy lost with BC per bin
        for( unsigned int ibin=0 ; ibin < Nbins ; ibin++ ) {
            nrj_bc_lost += nrj_lost_per_bin[ibin];
        } // end ibin
        } // end task for lost/radiated energy reduction 

    } else { // immobile particle
        if( diag_flag &&( !particles->is_test ) ) {
            if( params.geometry != "AMcylindrical" ) {

                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                    #pragma omp task default(shared) firstprivate(ibin) private(ithread,timer)
#else
                    #pragma omp task default(shared) firstprivate(ibin)
#endif
                    {
#ifdef  __DETAILED_TIMERS
                    ithread = omp_get_thread_num();
                    timer = MPI_Wtime();
#endif
                    for (unsigned int i = 0; i < size_proj_buffer_rho; i++) b_rho[ibin][i]   = 0.0;
                    for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basic( b_rho[ibin], ( *particles ), iPart, 0, ibin*clrw );
                    }
#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                    } // end task projection for frozen or test
                } // end ibin
            } else { // AM geometry
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );

                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
#ifdef  __DETAILED_TIMERS
                    #pragma omp task default(shared) firstprivate(ibin,bin_size0) private(ithread,timer)
#else
                    #pragma omp task default(shared) firstprivate(ibin,bin_size0) 
#endif
                    {
#ifdef  __DETAILED_TIMERS
                    ithread = omp_get_thread_num();
                    timer = MPI_Wtime();
#endif
                    for (unsigned int i = 0; i < size_proj_buffer_rhoAM; i++) b_rhoAM[ibin][i] = 0.0;
                    int imode = 0; // only mode 0 is used with envelope
                    for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basicForComplexOnBuffer( b_rhoAM[ibin], ( *particles ), iPart, 0, imode, bin_size0, ibin*clrw );
                    } // end loop on particles
#ifdef  __DETAILED_TIMERS
                    patch->patch_timers_[3*patch->thread_number_ + ithread] += MPI_Wtime() - timer;
#endif
                    } // end task projection for frozen or test
                } //end ibin
            } // end if on geometry
        } // end condition on diag and not particle test

    } // end if moving particle



    } // end taskgroup

} // End ponderomotiveUpdatePositionAndCurrentsTasks

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
