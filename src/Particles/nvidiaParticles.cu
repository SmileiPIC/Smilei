// -----------------------------------------------------------------------------
//
//! \file nvidiaParticles.cu
//
//! \brief contains the nvidiaParticles class methods
//! Extension of the Class Particles for GPU
//
// -----------------------------------------------------------------------------
// #if defined _GPU

// #include <thrust/host_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <iostream>

#include "nvidiaParticles.h"

//! Structure with specific function count_if_out for thrust::tuple operator
//! Return True if the entry is -1 as in the cell keys vector for instance
struct count_if_out
{
    __host__ __device__
  bool operator()(const int x)
  {
    return (x  == -1);
  }
};

//! Structure with specific function thrust::tuple operator
struct remove_if_out
{
    typedef thrust::tuple<double,double,double,double,double,double,double,short,int> Tuple;
    __host__ __device__
    bool operator()(const Tuple&  t)
    {
        // thrust::get<8>(t) = cell_keys
        return (thrust::get<8>(t) < 0);

    }
};

nvidiaParticles::nvidiaParticles() : Particles()
{
    gpu_nparts_ = 0;
}

// -----------------------------------------------------------------------------
//! Set capacity of Particles vectors
void nvidiaParticles::deviceReserve( unsigned int reserved_particles, unsigned int nDim)
{
    for (int idim=0;idim<nvidia_position_.size();idim++)
        nvidia_position_[idim].reserve( reserved_particles );
    for (int idim=0;idim<3;idim++) {
        nvidia_momentum_[idim].reserve( reserved_particles );
    }
    nvidia_weight_.reserve( reserved_particles );
    nvidia_charge_.reserve( reserved_particles );
    if( isQuantumParameter ) {
        nvidia_chi_.reserve( reserved_particles );
    }
    if( isMonteCarlo ) {
        nvidia_tau_.reserve( reserved_particles );
    }
    nvidia_cell_keys_.reserve( reserved_particles );
}

// -----------------------------------------------------------------------------
//! Set capacity of Particles vectors based on already used dimension on CPU
void nvidiaParticles::deviceReserve( unsigned int reserved_particles )
{
    deviceReserve(reserved_particles, (unsigned int) (Position.size()));
}

// ---------------------------------------------------------------------------------------------------------------------
// Reset of Particles vectors
// Cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void nvidiaParticles::deviceClear()
{
    for( unsigned int iprop=0 ; iprop<nvidia_double_prop_.size() ; iprop++ ) {
        nvidia_double_prop_[iprop]->clear();
    }

    for( unsigned int iprop=0 ; iprop<nvidia_short_prop_.size() ; iprop++ ) {
        nvidia_short_prop_[iprop]->clear();
    }

    //cell_keys.clear();
    gpu_nparts_ = 0;
}

// -----------------------------------------------------------------------------
//! Initialize the particle properties on devide as a mirror of the host definition
// -----------------------------------------------------------------------------
void nvidiaParticles::initializeDataOnDevice()
{
    if (!Position.size()) {
        std::cout << "Set CPU first" << std::endl;
        return ;
    }

    int ndim = Position.size();

    double res_fac = 1.;
    nvidia_position_.resize( Position.size() );
    for (int idim=0;idim<ndim;idim++) {
        nvidia_position_[idim].reserve( res_fac*Position[idim].size() );
        nvidia_position_[idim].resize( Position[idim].size() );
    }
    nvidia_momentum_.resize( Momentum.size() );
    for (int idim=0;idim<Momentum.size();idim++) {
        nvidia_momentum_[idim].reserve( res_fac*Momentum[idim].size() );
        nvidia_momentum_[idim].resize( Momentum[idim].size() );
    }
    nvidia_weight_.reserve( res_fac*Weight.size() );
    nvidia_weight_.resize( Weight.size() );
    nvidia_charge_.reserve( res_fac*Charge.size() );
    nvidia_charge_.resize( Charge.size() );

    if( isQuantumParameter ) {
        nvidia_chi_.reserve( res_fac*Chi.size() );
        nvidia_chi_.resize( Chi.size() );
    }

    if( isMonteCarlo ) {
        nvidia_tau_.reserve( res_fac*Tau.size() );
        nvidia_tau_.resize( Tau.size() );
    }

    nvidia_cell_keys_.reserve( res_fac*Charge.size() );
    nvidia_cell_keys_.resize( Charge.size() );
    gpu_nparts_ = Charge.size();

    if (gpu_nparts_!=0)
        syncGPU();
    else {
        deviceReserve(100);

        /*for (int idim=0;idim<Position.size();idim++)
            nvidia_position_[idim].reserve( 100 );
        for (int idim=0;idim<Momentum.size();idim++)
            nvidia_momentum_[idim].reserve( 100 );
        nvidia_weight_.reserve( 100 );
        nvidia_charge_.reserve( 100 );
        if( isQuantumParameter ) {
            nvidia_chi_.reserve( 100 );
        }
        if( isMonteCarlo ) {
            nvidia_tau_.reserve( 100 );
        }
        nvidia_cell_keys_.reserve( 100 );*/
    }

    // Initialize the list of pointers

    if( nvidia_double_prop_.empty() ) {  // do this just once

        for( unsigned int i=0 ; i< ndim ; i++ ) {
            nvidia_double_prop_.push_back( &( nvidia_position_[i] ) );
        }

        for( unsigned int i=0 ; i< 3 ; i++ ) {
            nvidia_double_prop_.push_back( &( nvidia_momentum_[i] ) );
        }

        nvidia_double_prop_.push_back( &nvidia_weight_ );

        nvidia_short_prop_.push_back( &nvidia_charge_ );

        // Quantum parameter (for QED effects):
        // - if radiation reaction (continuous or discontinuous)
        // - if multiphoton-Breit-Wheeler if photons
        if( isQuantumParameter ) {
            nvidia_double_prop_.push_back( &nvidia_chi_ );
        }

        // Optical Depth for Monte-Carlo processes:
        // - if the discontinuous (Monte-Carlo) radiation reaction
        // is activated, tau is the incremental optical depth to emission
        if( isMonteCarlo ) {
            nvidia_double_prop_.push_back( &nvidia_tau_ );
        }
    }
}

// Copy device to host
void nvidiaParticles::syncCPU()
{
    for (int idim=0;idim<Position.size();idim++) {
        Position[idim].resize( gpu_nparts_ );
        thrust::copy((nvidia_position_[idim]).begin(), (nvidia_position_[idim]).begin()+gpu_nparts_, (Position[idim]).begin());
    }
    for (int idim=0;idim<Momentum.size();idim++) {
        Momentum[idim].resize( gpu_nparts_ );
        thrust::copy((nvidia_momentum_[idim]).begin(), (nvidia_momentum_[idim]).begin()+gpu_nparts_, (Momentum[idim]).begin());
    }
    Weight.resize( gpu_nparts_ );
    thrust::copy((nvidia_weight_).begin(), (nvidia_weight_).begin()+gpu_nparts_, (Weight).begin());
    Charge.resize( gpu_nparts_ );
    thrust::copy((nvidia_charge_).begin(), (nvidia_charge_).begin()+gpu_nparts_, (Charge).begin());
    if (isQuantumParameter) {
        Chi.resize( gpu_nparts_ );
        thrust::copy((nvidia_chi_).begin(), (nvidia_chi_).begin()+gpu_nparts_, (Chi).begin());
    }
    if (isMonteCarlo) {
        Tau.resize( gpu_nparts_ );
        thrust::copy((nvidia_tau_).begin(), (nvidia_tau_).begin()+gpu_nparts_, (Tau).begin());
    }

}

//! Send the particles from host to device
void nvidiaParticles::syncGPU()
{
    for (int idim=0;idim<Position.size();idim++) {
        nvidia_position_[idim].resize( Position[idim].size() );
        thrust::copy((Position[idim]).begin(), (Position[idim]).end(), (nvidia_position_[idim]).begin());
    }
    for (int idim=0;idim<Momentum.size();idim++) {
        nvidia_momentum_[idim].resize( Momentum[idim].size() );
        thrust::copy((Momentum[idim]).begin(), (Momentum[idim]).end(), (nvidia_momentum_[idim]).begin());
    }
    nvidia_weight_.resize( Weight.size() );
    thrust::copy((Weight).begin(), (Weight).end(), (nvidia_weight_).begin());
    nvidia_charge_.resize( Charge.size() );
    thrust::copy((Charge).begin(), (Charge).end(), (nvidia_charge_).begin());
    if (isQuantumParameter) {
        nvidia_chi_.resize( Chi.size() );
        thrust::copy((Chi).begin(), (Chi).end(), (nvidia_chi_).begin());
    }
    if (isMonteCarlo) {
        nvidia_tau_.resize( Tau.size() );
        thrust::copy((Tau).begin(), (Tau).end(), (nvidia_tau_).begin());
    }
    gpu_nparts_ = Charge.size();
}

// -----------------------------------------------------------------------------
//! Extract particles from the Particles object and put 
//! them in the Particles object `particles_to_move`
// -----------------------------------------------------------------------------
void nvidiaParticles::extractParticles( Particles* particles_to_move )
{
    // TODO(Etienne M): Async copies would be great here! Launch all the kernels
    // and wait at the end.

    const int nparts                   = gpu_nparts_;
    const int position_dimension_count = nvidia_position_.size();

    int nparts_to_move_ = thrust::count( thrust::device,
                                         std::cbegin( nvidia_cell_keys_ ),
                                         std::cbegin( nvidia_cell_keys_ ) + nparts,
                                         -1 /* represent a particle out of boundary */ );

    nvidiaParticles* cp_parts = static_cast<nvidiaParticles*>( particles_to_move );
    cp_parts.gpu_nparts_      = nparts_to_move_;

    // Resize it, if too small (copy_if do not resize)
    if (nparts_to_move > cp_parts->nvidia_weight_.size() ) {
        for( int i = 0; i < position_dimension_count; ++i ) {
            cp_parts.nvidia_position_[i].resize( nparts_to_move_ );
        }
        cp_parts->nvidia_momentum_[0].resize( nparts_to_move );
        cp_parts->nvidia_momentum_[1].resize( nparts_to_move );
        cp_parts->nvidia_momentum_[2].resize( nparts_to_move );
        cp_parts->nvidia_weight_.resize( nparts_to_move );
        cp_parts->nvidia_charge_.resize( nparts_to_move );
        if (isQuantumParameter) {
            cp_parts->nvidia_chi_.resize( nparts_to_move );
        }
        if (isMonteCarlo) {
            cp_parts->nvidia_tau_.resize( nparts_to_move );
        }

        cp_parts->nvidia_cell_keys_.resize( nparts_to_move );
    }

    // Iterator of the main data structure
    // Note: https://nvidia.github.io/thrust/api/classes/classthrust_1_1zip__iterator.html#class-thrustzip_iterator
    const auto source_iterator_first      = thrust::make_zip_iterator( thrust::make_tuple( std::begin( nvidia_position_[0] ),
                                                                                           std::begin( nvidia_momentum_[0] ),
                                                                                           std::begin( nvidia_momentum_[1] ),
                                                                                           std::begin( nvidia_momentum_[2] ),
                                                                                           std::begin( nvidia_weight_ ),
                                                                                           std::begin( nvidia_charge_ ) ) );
    const auto source_iterator_last       = source_iterator_first + nparts; // std::advance
    const auto destination_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( std::begin( cp_parts.nvidia_position_[0] ),
                                                                                           std::begin( cp_parts.nvidia_momentum_[0] ),
                                                                                           std::begin( cp_parts.nvidia_momentum_[1] ),
                                                                                           std::begin( cp_parts.nvidia_momentum_[2] ),
                                                                                           std::begin( cp_parts.nvidia_weight_ ),
                                                                                           std::begin( cp_parts.nvidia_charge_ ) ) );

    // Copy send particles in dedicated data structure if nvidia_cell_keys_=0 (currently = 1 if keeped, new PartBoundCond::apply(...))
    thrust::copy_if( thrust::device,
                     source_iterator_first,
                     source_iterator_last,
                     // Copy depending on count_if_out()(nvidia_cell_keys_[i])
                     std::cbegin( nvidia_cell_keys_ ),
                     destination_iterator_first,
                     count_if_out() );

    // Copy the other position values depending on the simulation's grid
    // dimensions
    for( int i = 1; i < position_dimension_count; ++i ) {
        thrust::copy_if( thrust::device,
                         std::cbegin( nvidia_position_[i] ),
                         std::cbegin( nvidia_position_[i] ) + nparts,
                         std::cbegin( nvidia_cell_keys_ ),
                         std::begin( cp_parts.nvidia_position_[i] ),
                         count_if_out() );
    }

    // Special treatment for chi if radiation emission
    if( isQuantumParameter ) {
        thrust::copy_if( thrust::device,
                         std::cbegin( nvidia_chi_ ),
                         std::cbegin( nvidia_chi_ ) + nparts,
                         std::cbegin( nvidia_cell_keys_ ),
                         std::begin( cp_parts.nvidia_chi_ ),
                         count_if_out() );
    }

    if( isMonteCarlo ) {
        thrust::copy_if( thrust::device,
                         std::cbegin( nvidia_tau_ ),
                         std::cbegin( nvidia_tau_ ) + nparts,
                         std::cbegin( nvidia_cell_keys_ ),
                         std::begin( cp_parts.nvidia_tau_ ),
                         count_if_out() );
    }
    particles_to_move->syncCPU();
}

// -----------------------------------------------------------------------------
//! Erase particles leaving the patch object on device
// -----------------------------------------------------------------------------
int nvidiaParticles::eraseLeavingParticles() {
    
    const int position_dimension_count = nvidia_position_.size();
    
    const int nparts = gpu_nparts_;

    const int nparts_to_remove = thrust::count(thrust::device, 
                                               nvidia_cell_keys_.begin(),                    nvidia_cell_keys_.begin()+nparts, 
                                               -1);

    if( nparts_to_remove > 0 ) {

        const auto first_particle = thrust::make_zip_iterator( thrust::make_tuple( std::begin( nvidia_position_[0] ),
                                                                                   std::begin( nvidia_momentum_[0] ),
                                                                                   std::begin( nvidia_momentum_[1] ),
                                                                                   std::begin( nvidia_momentum_[2] ),
                                                                                   std::begin( nvidia_weight_ ),
                                                                                   std::begin( nvidia_charge_ ) ) );

        const auto last_particle  = first_particle + nparts;

        // Remove particles which leaves current patch
        thrust::remove_if( thrust::device,
                           first_particle,
                           last_particle,
                           std::cbegin( nvidia_cell_keys_ ),
                           count_if_out() );

        // Remove the other position values depending on the simulation's grid
        // dimensions
        for( int i = 1; i < position_dimension_count; ++i ) {
            thrust::remove_if( thrust::device,
                               std::begin( nvidia_position_[i] ),
                               std::begin( nvidia_position_[i] ) + nparts,
                               std::cbegin( nvidia_cell_keys_ ),
                               count_if_out() );
        }

        if( isQuantumParameter ) {
            thrust::remove_if( thrust::device,
                               std::begin( nvidia_chi_ ),
                               std::begin( nvidia_chi_ ) + nparts,
                               std::cbegin( nvidia_cell_keys_ ),
                               count_if_out() );
        }
        if( isMonteCarlo ) {
            thrust::remove_if( thrust::device,
                               std::begin( nvidia_tau_ ),
                               std::begin( nvidia_tau_ ) + nparts,
                               std::cbegin( nvidia_cell_keys_ ),
                               count_if_out() );
        }


        //std::cerr << "Removed particles: " << nparts - nvidia_weight_.size() << " " << nparts_to_remove << std::endl;

        // Update current number of particles
        gpu_nparts_ -= nparts_to_remove;

        // Reisze data structures (remove if does not resize)
        for( int i = 0; i < position_dimension_count; ++i ) {
            nvidia_position_[i].resize( tot_parts );
        }
        nvidia_momentum_[0].resize( gpu_nparts_ );
        nvidia_momentum_[1].resize( gpu_nparts_ );
        nvidia_momentum_[2].resize( gpu_nparts_ );
        nvidia_weight_.resize( gpu_nparts_ );
        nvidia_charge_.resize( gpu_nparts_ );
        nvidia_cell_keys_.resize( gpu_nparts_ );
        if (isQuantumParameter) {
            nvidia_chi_.resize( gpu_nparts_ );
        }
        if (isMonteCarlo) {
            nvidia_tau_.resize( gpu_nparts_ );
        }
        
    }

    return gpu_nparts_;
    
}


// -----------------------------------------------------------------------------
//! Inject particles from particles_to_inject object and put 
//! them in the Particles object
// -----------------------------------------------------------------------------
int nvidiaParticles::injectParticles( Particles* particles_to_inject )
{

    int nparts = gpu_nparts_;

    // Just resize cell keys, no need to remove
    // nvidia_cell_keys_.resize(gpu_nparts_);

    // Manage the recv data structure
    nvidiaParticles* cp_parts = static_cast<nvidiaParticles*>( particles_to_inject );

    const int nparts_add = cp_parts.gpu_nparts_;
    const int tot_parts  = nparts + nparts_add;

    // Resize main data structure, if too small (copy_n do not resize)
    if( tot_parts > nvidia_weight_.size() ) {
        for( int i = 0; i < position_dimension_count; ++i ) {
            nvidia_position_[i].resize( tot_parts );
        }

        nvidia_momentum_[0].resize( tot_parts );
        nvidia_momentum_[1].resize( tot_parts );
        nvidia_momentum_[2].resize( tot_parts );
        nvidia_weight_.resize( tot_parts );
        nvidia_charge_.resize( tot_parts );
        nvidia_cell_keys_.resize( tot_parts );

        if( isQuantumParameter ) {
            nvidia_chi_.resize( tot_parts );
        }

        if( isMonteCarlo ) {
            nvidia_tau_.resize( tot_parts );
        }
    }

    const auto source_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( std::cbegin( cp_parts.nvidia_position_[0] ),
                                                                                      std::cbegin( cp_parts.nvidia_momentum_[0] ),
                                                                                      std::cbegin( cp_parts.nvidia_momentum_[1] ),
                                                                                      std::cbegin( cp_parts.nvidia_momentum_[2] ),
                                                                                      std::cbegin( cp_parts.nvidia_weight_ ),
                                                                                      std::cbegin( cp_parts.nvidia_charge_ ) ) );

    // Iterator of the main data structure (once it has been resized)
    const auto destination_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( std::begin( nvidia_position_[0] ),
                                                                                           std::begin( nvidia_momentum_[0] ),
                                                                                           std::begin( nvidia_momentum_[1] ),
                                                                                           std::begin( nvidia_momentum_[2] ),
                                                                                           std::begin( nvidia_weight_ ),
                                                                                           std::begin( nvidia_charge_ ) ) ) +
                                            // /!\ Note the advance by nparts
                                            nparts;

    // Copy recv particles in main data structure
    thrust::copy_n( thrust::device,
                    source_iterator_first,
                    nparts_add,
                    destination_iterator_first );

    // Remove the other position values depending on the simulation's grid
    // dimensions
    for( int i = 1; i < position_dimension_count; ++i ) {
        thrust::copy_n( thrust::device,
                        std::cbegin( cp_parts.nvidia_position_[i] ),
                        nparts_add,
                        std::begin( nvidia_position_[i] ) + nparts );
    }

    if( isQuantumParameter ) {
        thrust::copy_n( thrust::device,
                        std::cbegin( cp_parts.nvidia_chi_ ),
                        nparts_add,
                        std::begin( nvidia_chi_ ) + nparts );
    }

    if( isMonteCarlo ) {
        thrust::copy_n( thrust::device,
                        std::cbegin( cp_parts.nvidia_tau_ ),
                        nparts_add,
                        std::begin( nvidia_tau_ ) + nparts );
    }

    // Resize below useless : nvidia_cell_keys_ resized if necessary above, cell_keys not used on cpu
    // nvidia_cell_keys_.resize( gpu_nparts_ );
    // cell_keys.resize       ( gpu_nparts_ );

    // No more particles to move
    cp_parts.gpu_nparts_ = 0;

    // Update number of particles
    gpu_nparts_ += nparts_add;

    return nparts_add;

}

// ---------------------------------------------------------------------------------------------------------------------
//! Create n_additional_particles new particles at the end of vectors
//! Fill the new elements with 0
// ---------------------------------------------------------------------------------------------------------------------
void nvidiaParticles::createParticles( int n_additional_particles )
{
    int n_particles = gpu_nparts_;
    int new_size = n_particles + n_additional_particles;
    for( unsigned int iprop=0 ; iprop<nvidia_double_prop_.size() ; iprop++ ) {
        ( *nvidia_double_prop_[iprop] ).resize(new_size);
         thrust::fill(( *nvidia_double_prop_[iprop] ).begin() + n_particles, ( *nvidia_double_prop_[iprop] ).begin() + new_size, 0);
    }

    for( unsigned int iprop=0 ; iprop<nvidia_short_prop_.size() ; iprop++ ) {
        ( *nvidia_short_prop_[iprop] ).resize(new_size);
        thrust::fill(( *nvidia_short_prop_[iprop] ).begin() + n_particles, ( *nvidia_short_prop_[iprop] ).begin() + new_size, 0);
    }
    
    // 
    // for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
    //     ( *nvidia_uint64_prop[iprop] ).resize( n_particles+n_additional_particles );
    // }

    nvidia_cell_keys_.resize( new_size );
    thrust::fill(nvidia_cell_keys_.begin() + n_particles, nvidia_cell_keys_.begin() + new_size, -1);

    gpu_nparts_ = new_size;

}

extern "C" {
void* CreateGPUParticles() {
    return new nvidiaParticles();
}
}

// #endif
