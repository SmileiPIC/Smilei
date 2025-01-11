#include "Species.h"

#include <omp.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>

#include "BoundaryConditionType.h"
#include "DiagnosticTrack.h"
#include "ElectroMagn.h"
#include "ElectroMagnAM.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "IonizationFactory.h"
#include "MergingFactory.h"
#include "MultiphotonBreitWheelerFactory.h"
#include "PartBoundCond.h"
#include "PartCompTimeFactory.h"
#include "PartWall.h"
#include "ParticleCreator.h"
#include "ParticlesFactory.h"
#include "Patch.h"
#include "Profile.h"
#include "Projector.h"
#include "ProjectorFactory.h"
#include "PusherFactory.h"
#include "RadiationFactory.h"
#include "SimWindow.h"
#include "Tools.h"
#include "gpu.h"


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
    iter_relativistic_initialization_( 0 ),
    ionization_model_( "none" ),
    density_profile_type_( "none" ),
    charge_profile_( NULL ),
    density_profile_( NULL ),
    velocity_profile_( 3, NULL ),
    radial_velocity_profile_( false ),
    temperature_profile_( 3, NULL ),
    particles_per_cell_profile_( NULL ),
    max_charge_( 0. ),
    // particles( &particles_sorted[0] ),
    file_position_npart_( 0 ),
    file_momentum_npart_( 0 ),
    position_initialization_array_( NULL ),
    momentum_initialization_array_( NULL ),
    n_numpy_particles_( 0 ),
    position_initialization_on_species_( false ),
    position_initialization_on_species_index_( -1 ),
    electron_species( NULL ),
    photon_species_( nullptr ),
    //photon_species_index(-1),
    radiation_photon_species( "" ),
    radiated_photons_( nullptr ),
    //mBW_pair_creation_sampling_( {1,1} ),
    mBW_pair_species_names_( 2, "" ),
    cluster_width_( params.cluster_width_ ),
    oversize( params.oversize ),
    cell_length( params.cell_length ),
    min_loc_vec( patch->getDomainLocalMin() ),
    tracking_diagnostic( 10000 ),
    nDim_particle( params.nDim_particle ),
    nDim_field(    params.nDim_field  ),
    merging_time_selection_( 0 )
{
    // &particles_sorted[0]
    particles         = ParticlesFactory::create( params, *patch );

    regular_number_array_.clear();
    partBoundCond = NULL;
    min_loc = patch->getDomainLocalMin( 0 );
    merging_method_ = "none";

    PI2 = 2.0 * M_PI;
    PI_ov_2 = 0.5*M_PI;

    dx_inv_[0] = 1./cell_length[0];
    dx_inv_[1] = 1./cell_length[1];
    dx_inv_[2] = 1./cell_length[2];

    initCluster( params, patch );
    inv_nDim_particles = 1./( ( double )nDim_particle );

    length_[0]=0;
    length_[1]=params.patch_size_[1]+1;
    length_[2]=params.patch_size_[2]+1;

    merge_momentum_cell_size_.resize(3);

    merge_min_momentum_cell_length_.resize(3);

    mBW_pair_species_index_[0] = -1;
    mBW_pair_species_index_[1] = -1;

    mBW_pair_creation_sampling_[0] = 1;
    mBW_pair_creation_sampling_[1] = 1;

}//END Species creator

void Species::initCluster( Params &params, Patch *patch )
{
    // NOTE: On GPU we dont use first_index, it would contain redundant data but
    // we are forced to initialize it due to ParticleCreator::create() and the
    // way the "structural" code of Smilei depends on it (maybe it could have
    // been delegated to the operators).

    // Arrays of the min and max indices of the particle bins
    particles->first_index.resize( params.patch_size_[0]/cluster_width_ );
    particles->last_index.resize( params.patch_size_[0]/cluster_width_ );
    Nbins = params.patch_size_[0]/cluster_width_ ;

    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)

    //Primal dimension of fields.
    f_dim0 =  params.patch_size_[0] + 2 * oversize[0] +1;
    f_dim1 =  params.patch_size_[1] + 2 * oversize[1] +1;
    f_dim2 =  params.patch_size_[2] + 2 * oversize[2] +1;

    //Dual dimension of fields.
    f_dim0_d =  params.patch_size_[0] + 2 * oversize[0] +2;
    f_dim1_d =  params.patch_size_[1] + 2 * oversize[1] +2;
    f_dim2_d =  params.patch_size_[2] + 2 * oversize[2] +2;

    b_dim.resize( 3, 1 );
    if( nDim_particle == 1 ) {
        b_dim[0] = ( 1 + cluster_width_ ) + 2 * oversize[0];
        b_dim[1] =  1;
        f_dim1 = 1;
        f_dim2 = 1;
        f_dim1_d = 1;
        f_dim2_d = 1;
    }
    if( nDim_particle == 2 ) {
        b_dim[0] = ( 1 + cluster_width_ ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] =  f_dim1;
        f_dim2 = 1;
        f_dim2_d = 1;
    }
    if( nDim_particle == 3 ) {
        b_dim[0] = ( 1 + cluster_width_ ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] = f_dim1;
        b_dim[2] = f_dim2;
    }

    //Initialize specMPI
    MPI_buffer_.allocate( params, patch );

    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_out = 0.;
    nrj_mw_inj = 0.;
    nrj_new_part_ = 0.;
    nrj_radiated_ = 0.;

}//END initCluster

// -----------------------------------------------------------------------------
//! This function enables to resize the number of bins
// -----------------------------------------------------------------------------
void Species::resizeCluster( Params &params )
{

    // We keep the current number of particles (from the vector size)
    int npart = particles->hostVectorSize();
    int size = params.patch_size_[0]/cluster_width_;

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

    // Area for particle creation
    struct SubSpace init_space;
    init_space.cell_index_[0] = 0;
    init_space.cell_index_[1] = 0;
    init_space.cell_index_[2] = 0;
    init_space.box_size_[0]   = params.patch_size_[0];
    init_space.box_size_[1]   = params.patch_size_[1];
    init_space.box_size_[2]   = params.patch_size_[2];

    // Creation of the particle creator
    ParticleCreator particle_creator;
    // Associate the ceator to the current species (this)
    particle_creator.associate(this);

    // If restart from a checkpoint or without particle creation
    if( params.restart || !with_particles ) {

        if( like_particles ) {
            particles->initialize( 0, *like_particles );
        } else {
            particles->initialize( 0, params.nDim_particle, params.keep_position_old );
        }

        // Compute only `max_charge_`
        particle_creator.createChargeProfile( init_space, patch);

    } else {

        // Create profiles and particles
        particle_creator.create( init_space, params, patch, 0 );

    }

}

// Initialize the operators (Push, Ionize, PartBoundCond)
// This must be separate from the parameters because the Species cloning copies
// the parameters but not the operators.
void Species::initOperators( Params &params, Patch *patch )
{

    // interpolation operator (virtual)
    Interp = InterpolatorFactory::create( params, patch, this->vectorized_operators ); // + patchId -> idx_domain_begin (now = ref smpi)

    // assign the correct Pusher to Push
    Push = PusherFactory::create( params, this );
    if( params.Laser_Envelope_model ) {
        Push_ponderomotive_position = PusherFactory::create_ponderomotive_position_updater( params, this );
    }

    // projection operator (virtual)
    Proj = ProjectorFactory::create( params, patch, this->vectorized_operators );  // + patchId -> idx_domain_begin (now = ref smpi)

    // Assign the Ionization model (if needed) to Ionize
    //  Needs to be placed after ParticleCreator() because requires the knowledge of max_charge_
    // \todo pay attention to restart
    Ionize = IonizationFactory::create( params, this );

    // Create the radiation model
    Radiate = RadiationFactory::create( params, this, patch->rand_ );

    // Create the multiphoton Breit-Wheeler model
    Multiphoton_Breit_Wheeler_process = MultiphotonBreitWheelerFactory::create( params, this, patch->rand_  );

    // assign the correct Merging method to Merge
    Merge = MergingFactory::create( this, patch->rand_ );

    // Evaluation of the particle computation time
    if (params.has_adaptive_vectorization ) {
        part_comp_time_ = PartCompTimeFactory::create( params );
    }

    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, this, patch );
    for( unsigned int iDim=0 ; iDim < nDim_field ; iDim++ ) {
        for( unsigned int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++ ) {
            MPI_buffer_.partRecv[iDim][iNeighbor]->initialize( 0, ( *particles ) );
            MPI_buffer_.partSend[iDim][iNeighbor]->initialize( 0, ( *particles ) );
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
    delete particles;
    
    delete Push;
    delete Interp;
    delete Proj;
    delete Merge;
    delete Ionize;
    delete Radiate;
    delete part_comp_time_;
    delete Multiphoton_Breit_Wheeler_process;
    delete partBoundCond;
    delete particles_per_cell_profile_;
    delete charge_profile_;
    delete density_profile_;
    for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
        delete velocity_profile_[i];
    }
    for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
        delete temperature_profile_[i];
    }
    if( ionization_rate_!=Py_None ) {
        Py_DECREF( ionization_rate_ );
    }
    delete radiated_photons_;
    for (int k=0 ; k<2 ; k++) {
        delete mBW_pair_particles_[k];
    }
    
    delete birth_records_;
    
}

#if defined( SMILEI_ACCELERATOR_GPU )
//! Prepare the species Current and Rho grids on Device
void
Species::prepareSpeciesCurrentAndChargeOnDevice( 
    unsigned int ispec,
    ElectroMagn * EMfields)
{

    unsigned int Jx_size;
    unsigned int Jy_size;
    unsigned int Jz_size;
    unsigned int rho_size;

    double * __restrict__ Jx_s = nullptr;
    double * __restrict__ Jy_s = nullptr;
    double * __restrict__ Jz_s = nullptr;
    double * __restrict__ rho_s = nullptr;

    if (EMfields->Jx_s[ispec]) {
        Jx_size             = EMfields->Jx_s[ispec]->size();
        Jx_s  = EMfields->Jx_s[ispec]->data() ;
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( Jx_s, Jx_size );
    }
    if (EMfields->Jy_s[ispec]) {
        Jy_size             = EMfields->Jy_s[ispec]->size();
        Jy_s  = EMfields->Jy_s[ispec]->data() ;
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( Jy_s, Jy_size );
    }
    if (EMfields->Jz_s[ispec]) {
        Jz_size             = EMfields->Jz_s[ispec]->size();
        Jz_s  = EMfields->Jz_s[ispec]->data() ;
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( Jz_s, Jz_size );
    }
    if (EMfields->rho_s[ispec]) {
        rho_size             = EMfields->rho_s[ispec]->size();
        rho_s  = EMfields->rho_s[ispec]->data() ;
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( rho_s, rho_size );
    }


#if defined( SMILEI_ACCELERATOR_GPU_OACC )
        #pragma acc parallel present( Jx_s[0:Jx_size],     \
                                        Jy_s[0:Jy_size], \
                                        Jz_s[0:Jz_size],   \
                                        rho_s[0:rho_size] )  
        {
#endif
        if (Jx_s) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<Jx_size; i++ ) {
                Jx_s[i] = 0;
            }
        }
        if (Jy_s) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<Jy_size; i++ ) {
                Jy_s[i] = 0;
            }
        }
        if (Jz_s) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<Jz_size; i++ ) {
                Jz_s[i] = 0;
            }
        }
        if (rho_s) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<rho_size; i++ ) {
                rho_s[i] = 0;
            }
        }
#if defined( SMILEI_ACCELERATOR_GPU_OACC )  
        } // end parallel region
#endif
}

//! Deallocate species Current (J) and Charge (Rho) arrays on Device
void
Species::deleteSpeciesCurrentAndChargeOnDevice(
    unsigned int ispec,
    ElectroMagn * EMfields)
{
    if (EMfields->Jx_s[ispec]) {
        // double *const __restrict__ pointer  = EMfields->Jx_s[ispec]->data() ;
        // const int size                      = EMfields->Jx_s[ispec]->size();
        // smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( pointer, size );
        EMfields->Jx_s[ispec]->deleteOnDevice();
    }
    if (EMfields->Jy_s[ispec]) {
        // double *const __restrict__ pointer  = EMfields->Jy_s[ispec]->data() ;
        // const int size                      = EMfields->Jy_s[ispec]->size();
        // smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( pointer, size );
        EMfields->Jy_s[ispec]->deleteOnDevice();
    }
    if (EMfields->Jz_s[ispec]) {
        // double *const __restrict__ pointer  = EMfields->Jz_s[ispec]->data() ;
        // const int size                      = EMfields->Jz_s[ispec]->size();
        // smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( pointer, size );
        EMfields->Jz_s[ispec]->deleteOnDevice();
    }
    if (EMfields->rho_s[ispec]) {
        // double *const __restrict__ pointer  = EMfields->rho_s[ispec]->data();
        // const int size                      = EMfields->rho_s[ispec]->size();
        // smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( pointer, size );
        EMfields->rho_s[ispec]->deleteOnDevice();
    }
}


void Species::allocateParticlesOnDevice()
{
    particles->initializeDataOnDevice();
    
    // The first send/recv buffers are also on device
    MPI_buffer_.partSend[0][0]->initializeDataOnDevice();
    MPI_buffer_.partSend[0][1]->initializeDataOnDevice();
    MPI_buffer_.partRecv[0][0]->initializeDataOnDevice();
    MPI_buffer_.partRecv[0][1]->initializeDataOnDevice();

    // Create photon species on the device
    if( radiation_model_ == "mc" && photon_species_ ) {
        radiated_photons_->initializeDataOnDevice();
    }

    // Create pair species on the device
    if( mBW_pair_species_[0] && mBW_pair_species_[1] ) {
        mBW_pair_particles_[0]->initializeDataOnDevice();
        mBW_pair_particles_[1]->initializeDataOnDevice();
    }
}


//! Copy particles from host to device
void
Species::copyParticlesFromHostToDevice()
{
    particles->copyFromHostToDevice();
}

#endif // end if SMILEI_ACCELERATOR_GPU

// ---------------------------------------------------------------------------------------------------------------------
//! Method calculating the Particle dynamics (interpolation, pusher, projection and more)
//! For all particles of the species
//!   - interpolate the fields at the particle position
//!   - perform ionization
//!   - perform the radiation reaction
//!   - calculate the new velocity
//!   - calculate the new position
//!   - apply the boundary conditions
//!   - increment the currents (projection)
void Species::dynamics( double time_dual,
                        unsigned int ispec,
                        ElectroMagn *EMfields,
                        Params &params,
                        bool diag_flag,
                        PartWalls *partWalls,
                        Patch *patch, SmileiMPI *smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables )
{
    int tid( 0 );

    const int ithread = Tools::getOMPThreadNum();

    unsigned int iPart;

    std::vector<double> nrj_lost_per_thd( 1, 0. );
    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize) { // moving particle

        // Prepare temporary buffers for this iteration
#if defined( SMILEI_ACCELERATOR_GPU )
        smpi->resizeDeviceBuffers( ithread,
                                   nDim_field,
                                   particles->numberOfParticles() );
#else
        smpi->resizeBuffers( ithread, nDim_field, particles->numberOfParticles(), params.geometry == "AMcylindrical" );
#endif

        // Prepare particles buffers for multiphoton Breit-Wheeler
        if( Multiphoton_Breit_Wheeler_process ) {

            patch->startFineTimer(mBW_timer_id_);

#if defined( SMILEI_ACCELERATOR_GPU_OACC) 
            static_cast<nvidiaParticles*>(mBW_pair_particles_[0])->deviceResize( particles->deviceSize() * Multiphoton_Breit_Wheeler_process->getPairCreationSampling(0) );
            static_cast<nvidiaParticles*>(mBW_pair_particles_[0])->resetCellKeys();
            static_cast<nvidiaParticles*>(mBW_pair_particles_[1])->deviceResize( particles->deviceSize() * Multiphoton_Breit_Wheeler_process->getPairCreationSampling(1) );
            static_cast<nvidiaParticles*>(mBW_pair_particles_[1])->resetCellKeys();
#else
            mBW_pair_particles_[0]->reserve(particles->numberOfParticles() * Multiphoton_Breit_Wheeler_process->getPairCreationSampling(0));
            mBW_pair_particles_[1]->reserve(particles->numberOfParticles() * Multiphoton_Breit_Wheeler_process->getPairCreationSampling(1));
#endif

            patch->stopFineTimer(mBW_timer_id_);
        }

#if defined( SMILEI_ACCELERATOR_GPU )
        // Make sure some bin preconditions are respected
        SMILEI_ASSERT( particles->first_index.size() == 1 );
        SMILEI_ASSERT( particles->last_index.size() >= 1 );
        SMILEI_ASSERT( particles->last_index.back() == particles->last_index[0] );
#endif
        for( unsigned int ibin = 0 ; ibin < particles->numberOfBins() ; ibin++ ) {

            patch->startFineTimer(interpolation_timer_id_);
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,0);

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);
            patch->stopFineTimer(interpolation_timer_id_);

            // Ionization
            if( Ionize ) {

                patch->startFineTimer(4);
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,5);
                ( *Ionize )( particles, particles->first_index[ibin], particles->last_index[ibin], &smpi->dynamics_Epart[ithread], patch, Proj );

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,5);
                patch->stopFineTimer(4);
            }

            if( time_dual<=time_frozen_ ) continue; // Do not push frozen particles

            // Radiation losses
            if( Radiate ) {

                patch->startFineTimer(5);

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,6);
                // Radiation process
                ( *Radiate )( *particles,
                              radiated_photons_,
                              smpi,
                              RadiationTables,
                              nrj_radiated_,
                              particles->first_index[ibin],
                              particles->last_index[ibin], ithread );
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,6);
                // Update scalar variable for diagnostics
                // nrj_radiated_ += Radiate->getRadiatedEnergy();

                // Update the quantum parameter chi
                // Radiate->computeParticlesChi( *particles,
                //                               smpi,
                //                               first_index[ibin],
                //                               last_index[ibin],
                //                               ithread );

                patch->stopFineTimer(5);

            }


            // Multiphoton Breit-Wheeler
            if( Multiphoton_Breit_Wheeler_process ) {

                patch->startFineTimer(6);
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,7);

                // Pair generation process
                // We reuse nrj_radiated_ for the pairs
                ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                        smpi,
                                                        mBW_pair_particles_,
                                                        mBW_pair_species_,
                                                        MultiphotonBreitWheelerTables,
                                                        nrj_radiated_,
                                                        particles->first_index[ibin],
                                                        particles->last_index[ibin], ithread );

                // Update the photon quantum parameter chi of all photons
                // Multiphoton_Breit_Wheeler_process->computeThreadPhotonChi( *particles,
                //         smpi,
                //         particles->first_index[ibin],
                //         particles->last_index[ibin],
                //         ithread );

                // Suppression of the decayed photons into pairs

                // Multiphoton_Breit_Wheeler_process->removeDecayedPhotons(
                //     *particles, smpi, ibin, particles->first_index.size(), &particles->first_index[0], &particles->last_index[0], ithread );

                Multiphoton_Breit_Wheeler_process->removeDecayedPhotonsWithoutBinCompression(
                    *particles, smpi, ibin,
                    &particles->first_index[0],
                    &particles->last_index[0],
                    ithread );


                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,7);
                patch->stopFineTimer(6);
            }

        } //ibin
        
        if( time_dual>time_frozen_){ // do not apply particles BC nor project frozen particles
            
            // Compression of the bins if necessary
            if( Multiphoton_Breit_Wheeler_process ) {

#ifdef SMILEI_ACCELERATOR_GPU_OACC
                removeTaggedParticles(smpi,
                                    &particles->first_index[0],
                                    &particles->last_index[0],
                                    ithread,
                                    false);
#else
                // Remove Photons while keeping the first index of each bin
                // Concerns as well the smpi buffers
                removeTaggedParticlesPerBin(smpi, ithread, false);

                // Delete the gap between the bins
                // Concerns as well the smpi buffers
                compress(smpi, ithread, true);
#endif
                patch->startFineTimer(6);

                // Remove Particles while keeping the first index of each bin
                // Concerns as well the smpi buffers
                removeTaggedParticlesPerBin(smpi, ithread, false);

                // Delete the gap between the bins
                // Concerns as well the smpi buffers
                compress(smpi, ithread, true);

                patch->stopFineTimer(6);

            }
            
            
// #ifdef  __DETAILED_TIMERS
//             timer = MPI_Wtime();
// #endif
            patch->startFineTimer(1);

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);

            size_t start = 0, stop = particles->size(), n = stop - start;
            vector<vector<double>> pold;
            particles->prepareInterpolatedFields( pold, start, n );

            // Push the particles and the photons
            ( *Push )( *particles, smpi, 0, particles->last_index.back(), ithread );
            //particles->testMove( particles->first_index[ibin], particles->last_index[ibin], params );
            
            // Copy interpolated fields to persistent buffers if requested
            if( particles->interpolated_fields_ ) {
                particles->copyInterpolatedFields( &( smpi->dynamics_Epart[ithread][start] ), &( smpi->dynamics_Bpart[ithread][start] ), pold, start, n, smpi->getBufferSize(ithread), mass_ );
            }
            
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

            patch->stopFineTimer(1);

// #ifdef  __DETAILED_TIMERS
//                 patch->patch_timers_[1] += MPI_Wtime() - timer;
// #endif

            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
                double energy_lost( 0. );

                patch->startFineTimer(3);

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,2);
                // Apply wall and boundary conditions
                if( mass_>0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        ( *partWalls )[iwall]->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                        nrj_lost_per_thd[tid] += mass_ * energy_lost;
                    }
                    // Boundary Condition may be physical or due to domain decomposition
                    if(!params.is_spectral){
                        partBoundCond->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                        nrj_lost_per_thd[tid] += mass_ * energy_lost;
                    }

                } else if( mass_==0 ) {
                    for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                        ( *partWalls )[iwall]->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                        nrj_lost_per_thd[tid] += energy_lost;
                    }
                    // Boundary Condition may be physical or due to domain decomposition
                    partBoundCond->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                    nrj_lost_per_thd[tid] += energy_lost;
                }
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,2);

                patch->stopFineTimer(3);

                //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

                patch->startFineTimer(2);
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);

                // Project currents if not a Test species and charges as well if a diag is needed.
                // Do not project if a photon
                if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                    Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
                }

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);
                patch->stopFineTimer(2);

                if(params.is_spectral && mass_>0){
                    partBoundCond->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                    nrj_lost_per_thd[tid] += mass_ * energy_lost;
                }

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

            patch->startFineTimer(2);
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);

            double *b_rho=nullptr;
            for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                // TODO(Etienne M): DIAGS. The projector needs to work on valid data. Currently, in GPU mode, it'll read
                // outdated particles data because basic() is always done on CPU. We need to pull the GPU data to the
                // host.
                b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                for( iPart=particles->first_index[ibin] ; ( int )iPart<particles->last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            }

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);
            patch->stopFineTimer(2);

        } else {
            int n_species = patch->vecSpecies.size();
            complex<double> *b_rho=nullptr;
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );

            patch->startFineTimer(2);
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);

            for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                int ifield = imode*n_species+ispec;
                b_rho = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                    for( int iPart=particles->first_index[ibin] ; iPart<particles->last_index[ibin]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);
            patch->stopFineTimer(2);

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
void Species::scalarDynamics( double, unsigned int ,
                               ElectroMagn *,
                               Params &, bool ,
                               PartWalls *,
                               Patch *, SmileiMPI *,
                               RadiationTables &,
                               MultiphotonBreitWheelerTables & )
{

}

void Species::projectionForDiags( unsigned int ispec,
                                  ElectroMagn *EMfields,
                                  Params &params, bool diag_flag,
                                  Patch *patch )
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
void Species::dynamicsImportParticles( double time_dual, Params &params, Patch *patch, vector<Diagnostic *> &localDiags )
{
    // Add the ionized electrons to the electron species (possible even if ion is frozen)
    if( Ionize ) {
        electron_species->importParticles( params, patch, Ionize->new_electrons, localDiags, time_dual, Ionize );
    }

    // if moving particle
    if( time_dual>time_frozen_ ) {

        // Radiation losses
        if( Radiate && photon_species_ ) {
            // If creation of macro-photon, we add them to photon_species
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                // We first erase empty slots in the buffer of photons
                // radiation_photons_->cell_keys is used as a mask
            static_cast<nvidiaParticles*>(radiated_photons_)->eraseLeavingParticles();
#endif
            photon_species_->importParticles( params, patch, *radiated_photons_, localDiags, time_dual );

#ifdef SMILEI_ACCELERATOR_GPU_OACC
            // We explicitely clear the device Particles
            static_cast<nvidiaParticles*>(radiated_photons_)->deviceClear();
#endif
        }

        // Multiphoton Breit-Wheeler
        if( Multiphoton_Breit_Wheeler_process ) {

            // Addition of the electron-positron particles
            for( int k=0; k<2; k++ ) {

#ifdef SMILEI_ACCELERATOR_GPU_OACC
                // We first erase empty slots in the buffer of photons
                // radiation_photons_->cell_keys is used as a mask
                static_cast<nvidiaParticles*>(mBW_pair_particles_[k])->eraseLeavingParticles();
#endif

                mBW_pair_species_[k]->importParticles( params, patch, *mBW_pair_particles_[k], localDiags, time_dual );
                
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                // We explicitely clear the device Particles
                static_cast<nvidiaParticles*>(mBW_pair_particles_[k])->deviceClear();
#endif
            }

        }
    }
}




// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - increment the charge (projection)
//   - used at initialisation for Poisson (and diags if required, not for now dynamics )
// ---------------------------------------------------------------------------------------------------------------------
void Species::computeCharge( ElectroMagn *EMfields, bool old /*=false*/ )
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


// ---------------------------------------------------------------------------------------------------------------------
//! Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::sortParticles( Params &params )
{

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )

    // -----------------------------
    // GPU version
    
    // Merge all MPI_buffer_.partRecv in the first one
    Particles * first_buffer = MPI_buffer_.partRecv[0][0];
    for( auto &partRecvs: MPI_buffer_.partRecv ) {
        for( auto partRecv: partRecvs ) {
            if( partRecv != first_buffer && partRecv->size() > 0 ) {
                partRecv->copyParticles( 0, partRecv->size(), *first_buffer, first_buffer->size() );
                partRecv->clear();
            }
        }
    }
    
    first_buffer->copyFromHostToDevice();
    
    particles->importAndSortParticles( first_buffer );
    
#else

    // --------------------------
    // CPU version

    // injectParticles( params );

    int ndim = params.nDim_field;
    int idim;

    // Sort to adapt do cell_keys usage
    std::vector<int> indexes_of_particles_to_exchange;
    for ( int ipart=0 ; ipart< (int)(getNbrOfParticles()) ; ipart++ ) {
        if ( particles->cell_keys[ipart] < 0 ) {
            indexes_of_particles_to_exchange.push_back( ipart );
        }
    }

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
            particles->overwriteParticle( particles->last_index[ibin]-iPart,
                                          particles->last_index[ibin-1],
                                          iPart, false );
        }
        particles->last_index[ibin] -= ii;
        particles->first_index[ibin] = particles->last_index[ibin-1];
    }



    int nmove, lmove; // local, OK
    int shift[particles->last_index.size()+1];//how much we need to shift each bin in order to leave room for the new particle
    double dbin;

    dbin = params.cell_length[0]*params.cluster_width_; //width of a bin.
    for( unsigned int j=0; j<particles->last_index.size()+1 ; j++ ) {
        shift[j]=0;
    }


    int nbNeighbors_ = 2;
    int n_part_recv;

    particles->eraseParticleTrail( particles->numberOfParticles() );

    //Evaluation of the necessary shift of all bins.2
    //idim=0
    shift[1] += MPI_buffer_.partRecv[0][0]->size();//Particles coming from xmin all go to bin 0 and shift all the other bins.
    shift[particles->last_index.size()] += MPI_buffer_.partRecv[0][1]->size();//Used only to count the total number of particles arrived.
    //idim>0
    for( idim = 1; idim < ndim; idim++ ) {
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            n_part_recv = MPI_buffer_.partRecv[idim][iNeighbor]->size();
            for( unsigned int j=0; j<( unsigned int )n_part_recv ; j++ ) {
                //We first evaluate how many particles arrive in each bin.
                ii = int( ( MPI_buffer_.partRecv[idim][iNeighbor]->position( 0, j )-min_loc )/dbin ); //bin in which the particle goes.
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
            particles->overwriteParticle( particles->first_index[j], particles->first_index[j]+lmove, nmove, false );
        }
        particles->first_index[j] += shift[j];
        particles->last_index[j] += shift[j];
    }

    //Space has been made now to write the arriving particles into the correct bins
    //idim == 0  is the easy case, when particles arrive either in first or last bin.
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        n_part_recv = MPI_buffer_.partRecv[0][iNeighbor]->size();
        //if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
        if( ( n_part_recv!=0 ) ) {
            ii = iNeighbor*( particles->last_index.size()-1 ); //0 if iNeighbor=0(particles coming from Xmin) and particles->last_index.size()-1 otherwise.
            MPI_buffer_.partRecv[0][iNeighbor]->overwriteParticle( 0, *particles, particles->last_index[ii], n_part_recv );
            particles->last_index[ii] += n_part_recv ;
        }
    }
    //idim > 0; this is the difficult case, when particles can arrive in any bin.
    for( idim = 1; idim < ndim; idim++ ) {
        //if (idim!=iDim) continue;
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            n_part_recv = MPI_buffer_.partRecv[idim][iNeighbor]->size();
            //if ( (neighbor_[idim][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
            if( ( n_part_recv!=0 ) ) {
                for( unsigned int j=0; j<( unsigned int )n_part_recv; j++ ) {
                    ii = int( ( MPI_buffer_.partRecv[idim][iNeighbor]->position( 0, j )-min_loc )/dbin ); //bin in which the particle goes.
                    MPI_buffer_.partRecv[idim][iNeighbor]->overwriteParticle( j, *particles, particles->last_index[ii] );
                    particles->last_index[ii] ++ ;
                }
            }
        }
    }


    //The width of one bin is cell_length[0] * cluster_width_.

    int p1, p2, first_index_init;
    unsigned int bin;
    double limit;


    //Backward pass
    for( bin=0; bin<particles->first_index.size()-1; bin++ ) { //Loop on the bins.
        limit = min_loc + ( bin+1 )*cell_length[0]*cluster_width_;
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
        limit = min_loc + bin*cell_length[0]*cluster_width_;
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
#endif
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function configures the type of species according to the default mode
//! regardless the number of particles per cell
// ---------------------------------------------------------------------------------------------------------------------
void Species::defaultConfigure( Params &, Patch * )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function configures the species according to the vectorization mode
// ---------------------------------------------------------------------------------------------------------------------
void Species::configuration( Params &, Patch * )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function reconfigures the species operators after evaluating
//! the best mode from the particle distribution
// ---------------------------------------------------------------------------------------------------------------------
void Species::reconfiguration( Params &, Patch * )
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

    nxy = params.patch_size_[0]*params.patch_size_[1];
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

        ixy = iy + ix*params.patch_size_[1];


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

        ixy = iy + ix*params.patch_size_[1];
        particles->overwriteParticle( ip, particles_sorted[token], indices[ixy] );
        indices[ixy]++;
    }

    particles = &particles_sorted[token] ;

}

// Move all particles from another species to this one
void Species::importParticles( Params &params, Patch *patch, Particles &source_particles, vector<Diagnostic *> &localDiags, double time_dual, Ionization *I )
{
#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )
    // ---------------------------------------------------
    // GPU version
    // Warning: the GPU version does not handle bin and sorting
    // Warning: the current GPU version does not handle tracked particles

    // Inject particles from source_particles
    particles->last_index.back() += particles->addParticles( &source_particles );
    particles->last_index[0] = particles->last_index.back();
    source_particles.resize( 0 );
    
#else
    // ---------------------------------------------------
    // CPU version


    const unsigned int npart     = source_particles.size();
    const unsigned int nbin      = particles->numberOfBins();
    const double inv_cell_length = 1./ params.cell_length[0];

    // If this species is tracked, set the particle IDs
    if( particles->tracked ) {
        dynamic_cast<DiagnosticTrack *>( localDiags[tracking_diagnostic] )->setIDs( source_particles );
    }
    
    // If there is a diagnostic for recording particle birth, then copy new particles to the buffer
    if( birth_records_ ) {
        birth_records_->update( source_particles, npart, time_dual, I );
    }
    
    // Move particles
    vector<int> src_bin_keys( npart, 0 );
    for( unsigned int i=0; i<npart; i++ ) {
        // Copy particle to the correct bin
        src_bin_keys[i] = source_particles.position( 0, i )*inv_cell_length - ( patch->getCellStartingGlobalIndex( 0 ) + params.oversize[0] );
        src_bin_keys[i] /= params.cluster_width_;
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

    // Clear all particles
    source_particles.clear();

    // Put capacity to 0 to save memory
    source_particles.shrinkToFit(true);

#endif

}

// ---------------------------------------------------------------------------------------------------------------------
//! This method eliminates the space gap between the bins
//! (presence of empty particles between the bins)
// ---------------------------------------------------------------------------------------------------------------------
void Species::compress(SmileiMPI *smpi, int ithread, bool compute_cell_keys) {

    // std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    // std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    //std::vector<double> *gamma = &( smpi->dynamics_invgf[ithread] );
    //std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    // std::vector<double> *deltaold = &( smpi->dynamics_deltaold[ithread] );

    // std::vector<std::complex<double>> *thetaold = NULL;
    // if ( smpi->dynamics_eithetaold.size() )
    //     thetaold = &( smpi->dynamics_eithetaold[ithread] );

    const int nparts = smpi->dynamics_Epart[ithread].size()/3;

#ifdef SMILEI_ACCELERATOR_GPU_OACC

    double *const __restrict__ weight =  particles->getPtrWeight();

    double *const __restrict__ position_x = particles->getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_particle > 1 ? particles->getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_particle > 2 ? particles->getPtrPosition( 2 ) : nullptr;

    double *const __restrict__ momentum_x = particles->getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles->getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles->getPtrMomentum(2);

    short *const __restrict__ charge = particles->getPtrCharge();

    double *const __restrict__ chi = particles->getPtrChi();
    double *const __restrict__ tau = particles->getPtrTau();
#endif

    double *const __restrict__ Ex = &( ( smpi->dynamics_Epart[ithread] )[0*nparts] );
    double *const __restrict__ Ey = &( ( smpi->dynamics_Epart[ithread] )[1*nparts] );
    double *const __restrict__ Ez = &( ( smpi->dynamics_Epart[ithread] )[2*nparts] );

    double *const __restrict__ Bx = &( ( smpi->dynamics_Bpart[ithread] )[0*nparts] );
    double *const __restrict__ By = &( ( smpi->dynamics_Bpart[ithread] )[1*nparts] );
    double *const __restrict__ Bz = &( ( smpi->dynamics_Bpart[ithread] )[2*nparts] );

    double *const __restrict__ gamma    = &( smpi->dynamics_invgf[ithread][0] );
    double *const __restrict__ deltaold = &( smpi->dynamics_deltaold[ithread][0] );
    int    *const __restrict__ iold     = &( smpi->dynamics_iold[ithread][0] );

    std::complex<double> * __restrict__ thetaold = nullptr;
    if ( smpi->dynamics_eithetaold.size() ) {
        thetaold = &(smpi->dynamics_eithetaold[ithread][0]);
    }

    // std::cerr << nparts << " " << Epart->size() << std::endl;

    const int nbin = particles->numberOfBins();

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel \
    present(Ex[0:nparts],Ey[0:nparts],Ez[0:nparts], \
    Bx[0:nparts], By[0:nparts], Bz[0:nparts], \
    deltaold[0:3*nparts], \
    iold[0:3*nparts], \
    gamma[0:nparts]) \
    deviceptr( \
        position_x,position_y,position_z, \
        momentum_x,momentum_y,momentum_z, \
        charge,weight,tau,chi)
    {
    #pragma acc loop seq
#endif

    for (auto ibin = 0 ; ibin < nbin-1 ; ibin++) {

        // Removal of the photons
        const unsigned int bin_gap = particles->first_index[ibin+1] - particles->last_index[ibin];

        if( bin_gap > 0 ) {

            // Determine first index and number of particles to copy.
            // We copy from first index to the end to limit the number of copy (more efficient than copying the full bin to keep the same order)

            // Compute the number of particles
            unsigned int copy_particle_number = 0;

            // Index from where we move the particles
            unsigned int copy_first_index = particles->last_index[ibin+1] - bin_gap;
            // Total number of particles in the bin [ibin+1]
            unsigned int particle_number = particles->last_index[ibin+1] - particles->first_index[ibin+1];

            // if copy_first_index < particles->first_index[ibin+1], it means that the empty space is larger than the number of particles in ibin
            // then we move the full bin
            // Else we only move the particles from copy_first_index to last_index[ibin+1]
            if( (int) copy_first_index < particles->first_index[ibin+1] ) {
                copy_first_index = particles->first_index[ibin+1];
                copy_particle_number = particle_number;
            } else {
                copy_particle_number = bin_gap;
            }

            if (copy_particle_number>0) {

#ifndef SMILEI_ACCELERATOR_GPU_OACC
                particles->overwriteParticle(copy_first_index, particles->last_index[ibin], copy_particle_number, compute_cell_keys );
#else
                for (auto ipart = 0 ; ipart < copy_particle_number ; ipart ++) {
                    const auto ipart_l = copy_first_index + ipart;
                    const auto ipart_r = particles->last_index[ibin] + ipart;
                    weight[ipart_l] = weight[ipart_r];
                    position_x[ipart_l] = position_x[ipart_r];
                    if( nDim_particle > 1 ) {
                        position_y[ipart_l] = position_y[ipart_r];
                        if( nDim_particle > 2 ) {
                            position_z[ipart_l] = position_z[ipart_r];
                        }
                    }
                    momentum_x[ipart_l] = momentum_x[ipart_r];
                    momentum_y[ipart_l] = momentum_y[ipart_r];
                    momentum_z[ipart_l] = momentum_z[ipart_r];
                    charge[ipart_l] = charge[ipart_r];
                    if( particles->has_quantum_parameter ) {
                        chi[ipart_l] = chi[ipart_r];
                    }
                    if( particles->has_Monte_Carlo_process ) {
                        tau[ipart_l] = tau[ipart_r];
                    }
                }
#endif

                for( unsigned int ipart = 0 ; ipart < copy_particle_number ; ipart ++ ) {
                    Ex[copy_first_index + ipart] = Ex[particles->last_index[ibin] + ipart];
                    Ey[copy_first_index + ipart] = Ey[particles->last_index[ibin] + ipart];
                    Ez[copy_first_index + ipart] = Ez[particles->last_index[ibin] + ipart];
                }

                for( unsigned int ipart = 0 ; ipart < copy_particle_number ; ipart ++ ) {
                    Bx[copy_first_index + ipart] = Bx[particles->last_index[ibin] + ipart];
                    By[copy_first_index + ipart] = By[particles->last_index[ibin] + ipart];
                    Bz[copy_first_index + ipart] = Bz[particles->last_index[ibin] + ipart];
                }

                for( unsigned int ipart = 0 ; ipart < copy_particle_number ; ipart ++ ) {
                    gamma[copy_first_index + ipart] = gamma[particles->last_index[ibin] + ipart];
                }

                for( unsigned int ipart = 0 ; ipart < copy_particle_number ; ipart ++ ) {
                    for ( int iDim=particles->dimension()-1; iDim>=0 ; iDim-- ) {
                        iold[iDim*nparts + copy_first_index + ipart] = iold[iDim*nparts + particles->last_index[ibin] + ipart];
                    }
                }

                for( unsigned int ipart = 0 ; ipart < copy_particle_number ; ipart ++ ) {
                    for ( int iDim=particles->dimension()-1; iDim>=0 ; iDim-- ) {
                        deltaold[iDim*nparts + copy_first_index + ipart] = deltaold[iDim*nparts + particles->last_index[ibin] + ipart];
                    }
                }

#ifndef SMILEI_ACCELERATOR_GPU_OACC
                if (thetaold) {
                    for( unsigned int ipart = 0 ; ipart < copy_particle_number ; ipart ++ ) {
                        thetaold[copy_first_index + ipart] = thetaold[particles->last_index[ibin] + ipart];
                    }
                }
#endif
            }
            //particles->eraseParticle( particles->last_index[ibin], bin_gap, true );

            // const int nparts = Epart->size()/3;

            // Erase bufferised data
            // for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
            //     Epart->erase(Epart->begin()+iDim*nparts+particles->last_index[ibin],Epart->begin()+iDim*nparts+particles->last_index[ibin]+bin_gap);
            //     Bpart->erase(Bpart->begin()+iDim*nparts+particles->last_index[ibin],Bpart->begin()+iDim*nparts+particles->last_index[ibin]+bin_gap);
            // }
            // for ( int iDim=particles->dimension()-1; iDim>=0 ; iDim-- ) {
            //     iold->erase(iold->begin()+iDim*nparts+particles->last_index[ibin],iold->begin()+iDim*nparts+particles->last_index[ibin]+bin_gap);
            //     deltaold->erase(deltaold->begin()+iDim*nparts+particles->last_index[ibin],deltaold->begin()+iDim*nparts+particles->last_index[ibin]+bin_gap);
            // }
            // gamma->erase(gamma->begin()+0*nparts+particles->last_index[ibin],gamma->begin()+0*nparts+particles->last_index[ibin]+bin_gap);
            //
            // if (thetaold) {
            //     thetaold->erase(thetaold->begin()+0*nparts+particles->last_index[ibin],thetaold->begin()+0*nparts+particles->last_index[ibin]+bin_gap);
            // }

            // for( int ii=ibin+1; ii<nbin; ii++ ) {
            //     particles->first_index[ii] -= bin_gap;
            //     particles->last_index[ii] -= bin_gap;
            // }

            particles->first_index[ibin+1] = particles->last_index[ibin];
            particles->last_index[ibin+1] = particles->first_index[ibin+1] + particle_number;

        }
    }

#ifdef SMILEI_ACCELERATOR_GPU_OACC
} // end parallel region
#endif

    // Old particles (deleted particles) are now at the end of the vectors
    // Erase trailing particles
    particles->eraseParticleTrail( particles->last_index[nbin-1], true );
    // smpi->eraseBufferParticleTrail( particles->dimension(), particles->last_index[nbin-1], ithread );
}

//! This method removes particles with a negative weight
//! without changing the bin first index
//! Bins are therefore potentially seperated by empty particle slots
void Species::removeTaggedParticlesPerBin(
    SmileiMPI *smpi,
    int ithread,
    bool compute_cell_keys)
{
    // Buffers for particles
    double *const __restrict__ Epart     = smpi->dynamics_Epart[ithread].data();
    double *const __restrict__ Bpart     = smpi->dynamics_Bpart[ithread].data();
    double *const __restrict__ gamma     = smpi->dynamics_invgf[ithread].data();
    int *const __restrict__ iold         = smpi->dynamics_iold[ithread].data();
    double *const __restrict__ deltaold  = smpi->dynamics_deltaold[ithread].data();

    std::complex<double> * __restrict__ thetaold = NULL;
    if ( smpi->dynamics_eithetaold.size() )
        thetaold = smpi->dynamics_eithetaold[ithread].data();

    const int nparts = smpi->getBufferSize(ithread);

    // Weight shortcut
    double *const __restrict__ weight =  particles->getPtrWeight();

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    double *const __restrict__ position_x = particles->getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_particle > 1 ? particles->getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_particle > 2 ? particles->getPtrPosition( 2 ) : nullptr;

    double *const __restrict__ momentum_x = particles->getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles->getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles->getPtrMomentum(2);

    short *const __restrict__ charge = particles->getPtrCharge();

    double *const __restrict__ chi = particles->getPtrChi();
    double *const __restrict__ tau = particles->getPtrTau();
#endif

    // Total number of bins / cells
    const int nbin = particles->numberOfBins();

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc parallel  \
    present(Epart[0:nparts*3],\
    Bpart[0:nparts*3], \
    gamma[0:nparts], \
    iold[0:nparts*nDim_particle], \
    deltaold[0:nparts*nDim_particle]) \
    deviceptr( \
        position_x,position_y,position_z, \
        momentum_x,momentum_y,momentum_z, \
        charge,weight,tau,chi)
    {
    #pragma acc loop gang worker
#endif

    // loop over the bins
    for (auto ibin = 0 ; ibin < nbin; ibin++) {

        if( particles->last_index[ibin] > particles->first_index[ibin] ) {

            //int nb_deleted_photon;

            // Backward loop over the photons to find the first existing photon
            int last_photon_index = particles->last_index[ibin]-1; // Index of the last existing photon (weight > 0)
            int first_photon_index = particles->first_index[ibin]; // Index of the first photon
            while( ( last_photon_index >= particles->first_index[ibin] )
                    && ( weight[last_photon_index] <= 0 ) ) {
                last_photon_index--;
            }
            while( ( first_photon_index < particles->last_index[ibin] )
                    && ( weight[first_photon_index] > 0 ) ) {
                first_photon_index++;
            }
            // At this level, last_photon_index is the position of the last still-existing photon (weight > 0)
            // that will not be erased

            // Backward loop over the photons to fill holes in the photon particle array (at the bin level only)
            for( int ipart=last_photon_index-1 ; ipart>=particles->first_index[ibin]; ipart-- ) {
                if( weight[ipart] <= 0 ) {
                    if( ipart < last_photon_index ) {
                        // The last existing photon comes to the position of
                        // the deleted photon
#ifndef SMILEI_ACCELERATOR_GPU_OACC
                        particles->overwriteParticle( last_photon_index, ipart, compute_cell_keys );
#else
                        weight[ipart] = weight[last_photon_index];
                        position_x[ipart] = position_x[last_photon_index];
                        if( nDim_particle > 1 ) {
                            position_y[ipart] = position_y[last_photon_index];
                            if( nDim_particle > 2 ) {
                                position_z[ipart] = position_z[last_photon_index];
                            }
                        }
                        momentum_x[ipart] = momentum_x[last_photon_index];
                        momentum_y[ipart] = momentum_y[last_photon_index];
                        momentum_z[ipart] = momentum_z[last_photon_index];
                        charge[ipart] = charge[last_photon_index];
                        if( particles->has_quantum_parameter ) {
                            chi[ipart] = chi[last_photon_index];
                        }
                        if( particles->has_Monte_Carlo_process ) {
                            tau[ipart] = tau[last_photon_index];
                        }

#endif
                        // Overwrite bufferised data
                        for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
                            Epart[iDim*nparts+ipart] = Epart[iDim*nparts+last_photon_index];
                            Bpart[iDim*nparts+ipart] = Bpart[iDim*nparts+last_photon_index];
                        }
                        for ( int iDim=nDim_particle-1 ; iDim>=0 ; iDim-- ) {
                            iold[iDim*nparts+ipart] = iold[iDim*nparts+last_photon_index];
                            deltaold[iDim*nparts+ipart] = deltaold[iDim*nparts+last_photon_index];
                        }
                        gamma[ipart] = gamma[0*nparts+last_photon_index];

#ifndef SMILEI_ACCELERATOR_GPU_OACC
                        if (thetaold) {
                            thetaold[0*nparts+ipart] = thetaold[0*nparts+last_photon_index];
                        }
#endif
                        last_photon_index --;
                    }
                }
            } // end for ipart

            // Update of the bin boundaries
            // const unsigned int nb_deleted_photon = last_index[ibin]-last_photon_index-1;

            // We suppress the deleted photons
            if( last_photon_index + 1 < particles->last_index[ibin] ) {
                particles->last_index[ibin] = last_photon_index+1;

                // std::cerr
                //         << " ibin: " << ibin
                //         << " - first_index: " << first_index[ibin]
                //         << " - last_index: " << last_index[ibin]
                //         << " - nb_deleted_photon: " << nb_deleted_photon
                //         << std::endl;
            }
        } // if last_index[ibin] > first_index[ibin]
    } // end loop over the bins

#ifdef SMILEI_ACCELERATOR_GPU_OACC
    } // end parallel region
#endif
}

//! This method removes particles with a negative weight
//! when a single bin is used
#ifdef SMILEI_ACCELERATOR_GPU_OACC
void Species::removeTaggedParticles(
    SmileiMPI *smpi,
    int *const first_index,
    int *const last_index,
    int ithread,
    bool compute_cell_keys)
{

    unsigned int new_n_parts = 0;
    unsigned int nb_deleted  = 0;

    // Buffers for particles
    double *const __restrict__ Epart     = smpi->dynamics_Epart[ithread].data();
    double *const __restrict__ Bpart     = smpi->dynamics_Bpart[ithread].data();
    double *const __restrict__ gamma     = smpi->dynamics_invgf[ithread].data();
    int    *const __restrict__ iold      = smpi->dynamics_iold[ithread].data();
    double *const __restrict__ deltaold  = smpi->dynamics_deltaold[ithread].data();

    std::complex<double> * __restrict__ thetaold = NULL;
    if ( smpi->dynamics_eithetaold.size() )
        thetaold = smpi->dynamics_eithetaold[ithread].data();

    const int nparts = smpi->getBufferSize(ithread);
    const int nparts_thetaold = nparts * (thetaold ? 1: 0);

    // Weight shortcut
    double *const __restrict__ weight =  particles->getPtrWeight();

    double *const __restrict__ position_x = particles->getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_particle > 1 ? particles->getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_particle > 2 ? particles->getPtrPosition( 2 ) : nullptr;

    double *const __restrict__ momentum_x = particles->getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles->getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles->getPtrMomentum(2);

    short *const __restrict__ charge = particles->getPtrCharge();

    double *const __restrict__ chi = particles->has_quantum_parameter ? particles->getPtrChi() : nullptr;
    double *const __restrict__ tau = particles->has_Monte_Carlo_process ? particles->getPtrTau() : nullptr;

    // Only if there are particles
    if( nparts > 0 ) {

    int last_moving_index = *last_index-1; // Index of the last existing photon (weight > 0)

//    #pragma acc kernels
    #pragma acc serial  \
    present(Epart[0:nparts*3],\
    Bpart[0:nparts*3], \
    gamma[0:nparts], \
    iold[0:nparts*nDim_particle], \
    deltaold[0:nparts*nDim_particle], \
    thetaold[0:nparts_thetaold]) \
    deviceptr( \
        position_x,position_y,position_z, \
        momentum_x,momentum_y,momentum_z, \
        charge,weight,tau,chi)
    {
        //int nb_deleted_photon;

        // Backward loop over the tagged particles to find the first existing photon
        int first_moving_index = *first_index ; // Index of the first photon
        while( ( last_moving_index >= *first_index )
                && ( weight[last_moving_index] <= 0 ) ) {
            last_moving_index--;
        }
        while( ( first_moving_index < *last_index )
                && ( weight[first_moving_index] > 0 ) ) {
            first_moving_index++;
        }
        // At this level, last_photon_index is the position of the last still-existing photon (weight > 0)
        // that will not be erased

        // Backward loop over the tagged particles to fill holes in the photon particle array (at the bin level only)
//#ifdef SMILEI_ACCELERATOR_GPU_OACC
//        #pragma acc loop seq
//#endif
        for( int ipart=last_moving_index-1 ; ipart>=*first_index; ipart-- ) {
            if( weight[ipart] <= 0 ) {
                if( ipart < last_moving_index ) {
                    // The last existing photon comes to the position of
                    // the deleted photon

                    weight[ipart] = weight[last_moving_index];
                    position_x[ipart] = position_x[last_moving_index];
                    if( nDim_particle > 1 ) {
                        position_y[ipart] = position_y[last_moving_index];
                        if( nDim_particle > 2 ) {
                            position_z[ipart] = position_z[last_moving_index];
                        }
                    }
                    momentum_x[ipart] = momentum_x[last_moving_index];
                    momentum_y[ipart] = momentum_y[last_moving_index];
                    momentum_z[ipart] = momentum_z[last_moving_index];
                    charge[ipart] = charge[last_moving_index];
                    if( particles->has_quantum_parameter ) {
                        chi[ipart] = chi[last_moving_index];
                    }
                    if( particles->has_Monte_Carlo_process ) {
                        tau[ipart] = tau[last_moving_index];
                    }

                    // Overwrite bufferised data
                    for ( int iDim=2 ; iDim>=0 ; iDim-- ) {
                        Epart[iDim*nparts+ipart] = Epart[iDim*nparts+last_moving_index];
                        Bpart[iDim*nparts+ipart] = Bpart[iDim*nparts+last_moving_index];
                    }
                    for ( int iDim=nDim_particle-1 ; iDim>=0 ; iDim-- ) {
                        iold[iDim*nparts+ipart] = iold[iDim*nparts+last_moving_index];
                        deltaold[iDim*nparts+ipart] = deltaold[iDim*nparts+last_moving_index];
                    }
                    gamma[ipart] = gamma[0*nparts+last_moving_index];

                    if (thetaold) {
                        thetaold[0*nparts+ipart] = thetaold[0*nparts+last_moving_index];
                    }

                    last_moving_index --;
                }
            }
        } // end for ipart

        // We suppress the deleted photons

        // Removal of the photons
        unsigned int new_n_parts = last_moving_index + 1;
        unsigned int nb_deleted = *last_index-new_n_parts;

        // Update the bin size
        *last_index = new_n_parts;

        // Update the buffers (remove empty space between them)
        for (auto ip = 0 ; ip < nb_deleted ; ip++) {
            for ( int idim=1 ; idim<2 ; idim++ ) {
                Epart[idim*new_n_parts+ip] = Epart[idim*nparts+nparts-1-ip];
                Bpart[idim*new_n_parts+ip] = Bpart[idim*nparts+nparts-1-ip];
            }
            for ( int idim=1; idim < nDim_particle; idim++ ) {
                iold[idim*new_n_parts+ip] = iold[idim*nparts+nparts-1-ip];
                deltaold[idim*new_n_parts+ip] = deltaold[idim*nparts+nparts-1-ip];
            }
        }

    } // end openacc region

    if( nb_deleted > 0 ) {
        // Update the size of the particle vectors
        static_cast<nvidiaParticles*>(particles)->deviceResize(new_n_parts);
    }
    } // if nparts > 0

}
#endif


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
void Species::mergeParticles( double /*time_dual*/ )
{
}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
//   - calculate the new momentum
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotiveUpdateSusceptibilityAndMomentum( double time_dual,
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

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ || Ionize) { // moving particle

        smpi->resizeBuffers( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->numberOfBins() ; ibin++ ) { // loop on ibin

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,0);
            Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,0);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[7] += MPI_Wtime() - timer;
#endif

            // Ionization
            if( Ionize ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,5);
                vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
                vector<double> *EnvEabs_part = &( smpi->dynamics_EnvEabs_part[ithread] );
                vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );
                vector<double> *Phipart = &( smpi->dynamics_PHIpart[ithread] );
                Interp->envelopeFieldForIonization( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
                Ionize->envelopeIonization( particles, particles->first_index[ibin], particles->last_index[ibin], Epart, EnvEabs_part, EnvExabs_part, Phipart, patch, Proj );
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,5);

#ifdef  __DETAILED_TIMERS
                patch->patch_timers_[4] += MPI_Wtime() - timer;
#endif
            }

            if( time_dual<=time_frozen_ ) continue; // Do not push nor project frozen particles

            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,3);
            Proj->susceptibility( EMfields, *particles, mass_, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,3);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[8] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),0,1);
            // Push only the particle momenta
            ( *Push )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, Tools::getOMPThreadNum(),1,1);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[9] += MPI_Wtime() - timer;
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
void Species::ponderomotiveProjectSusceptibility( double time_dual, 
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

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->resizeBuffers( ithread, nDim_particle, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->numberOfBins() ; ibin++ ) { // loop on ibin

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
    } //END if time vs. time_frozen_
    
    SMILEI_UNUSED( patch );
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
        Patch *patch, SmileiMPI *smpi )
{

    const int ithread = Tools::getOMPThreadNum();

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    bool diag_PartEventTracing {false};

# ifdef _PARTEVENTTRACING
    diag_PartEventTracing = smpi->diagPartEventTracing( time_dual, params.timestep);
# endif

    //unsigned int iPart;

    int tid( 0 );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle updated position
    // -------------------------------
    if( time_dual>time_frozen_ ) { // moving particle

        smpi->resizeBuffers( ithread, nDim_field, particles->last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin++ ) {
            double energy_lost( 0. );

            // Interpolate the ponderomotive potential and its gradient at the particle position, present and previous timestep
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,0);
            Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( particles->first_index[ibin] ), &( particles->last_index[ibin] ), ithread );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,0);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[10] += MPI_Wtime() - timer;
#endif

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,1);
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread );
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,1);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers_[11] += MPI_Wtime() - timer;
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,2);
            // Apply wall and boundary conditions
            if( mass_>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    ( *partWalls )[iwall]->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                    nrj_lost_per_thd[tid] += mass_ * energy_lost;
                }

                // Boundary Condition may be physical or due to domain decomposition
                partBoundCond->apply( this, particles->first_index[ibin], particles->last_index[ibin], smpi->dynamics_invgf[ithread], patch->rand_, energy_lost );
                nrj_lost_per_thd[tid] += mass_ * energy_lost;

            } else if( mass_==0 ) {
                ERROR( "Particles with zero mass cannot interact with envelope" );

            } // end mass_ = 0? condition
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,2);

            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,3);
            if( ( !particles->is_test ) && ( mass_ > 0 ) ) {
                Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, particles->first_index[ibin], particles->last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
            }
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,3);

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

            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,3);
            double *b_rho=nullptr;
            if( params.geometry != "AMcylindrical" ) {
                b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                for( unsigned int ibin = 0 ; ibin < particles->first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                    for( unsigned int iPart = (unsigned int)(particles->first_index[ibin]) ; iPart < (unsigned int)(particles->last_index[ibin]); iPart++ ) {
                        Proj->basic( b_rho, ( *particles ), iPart, 0 );
                    }
                }//End loop on bins
            smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,3);

            } else {

                smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,0,3);
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
                    }//End loop on bins
                } // end loop on modes
                smpi->traceEventIfDiagTracing(diag_PartEventTracing, ithread,1,3);

            } // end if on geometry
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
              << " in patch (" << patch->Pcoordinates[0];
    for( unsigned int idim = 1; idim<patch->Pcoordinates.size(); idim++ ) {
        std::cerr << "," <<  patch->Pcoordinates[idim];
    }
    std::cerr << ") "
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

//! Erase all particles with zero weight
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
