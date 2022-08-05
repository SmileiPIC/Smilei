// -----------------------------------------------------------------------------
//
//! \file nvidiaParticles.cu
//
//! \brief contains the nvidiaParticles class methods
//! Extension of the Class Particles for GPU
//
// -----------------------------------------------------------------------------

#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include "nvidiaParticles.h"

//! Structure with specific function count_if_out for thrust::tuple operator
//! Return True if the entry is -1 as in the cell keys vector for instance
struct count_if_out
{
    __host__ __device__ bool
    operator()( const short x )
    {
        return x == -1;
    }
};

nvidiaParticles::nvidiaParticles( const Params& parameters )
    : Particles{}
    , gpu_nparts_{}
{
    // EMPTY
}

void nvidiaParticles::allocateDimensions( unsigned int nDim )
{
    nvidia_position_.resize( nDim );
    nvidia_momentum_.resize( 3 );
}

void nvidiaParticles::reserveParticles( unsigned int particle_count )
{
    for( unsigned int idim = 0; idim < nvidia_position_.size(); idim++ ) {
        nvidia_position_[idim].reserve( particle_count );
    }

    for( unsigned int idim = 0; idim < 3; idim++ ) {
        nvidia_momentum_[idim].reserve( particle_count );
    }

    nvidia_weight_.reserve( particle_count );
    nvidia_charge_.reserve( particle_count );

    if( isQuantumParameter ) {
        nvidia_chi_.reserve( particle_count );
    }

    if( isMonteCarlo ) {
        nvidia_tau_.reserve( particle_count );
    }

    nvidia_cell_keys_.reserve( particle_count );
}

void nvidiaParticles::allocateParticles( unsigned int particle_count )
{
    for( int idim = 0; idim < nvidia_position_.size(); idim++ ) {
        nvidia_position_[idim].resize( particle_count );
    }

    for( int idim = 0; idim < 3; idim++ ) {
        nvidia_momentum_[idim].resize( particle_count );
    }

    nvidia_weight_.resize( particle_count );
    nvidia_charge_.resize( particle_count );

    if( isQuantumParameter ) {
        nvidia_chi_.resize( particle_count );
    }

    if( isMonteCarlo ) {
        nvidia_tau_.resize( particle_count );
    }

    nvidia_cell_keys_.resize( particle_count );
}

// ---------------------------------------------------------------------------------------------------------------------
// Reset of Particles vectors
// Cell keys not affected
// ---------------------------------------------------------------------------------------------------------------------
void nvidiaParticles::deviceClear()
{
    for( unsigned int iprop = 0; iprop < nvidia_double_prop_.size(); iprop++ ) {
        nvidia_double_prop_[iprop]->clear();
    }

    for( unsigned int iprop = 0; iprop < nvidia_short_prop_.size(); iprop++ ) {
        nvidia_short_prop_[iprop]->clear();
    }

    // TODO(Etienne M): Clear cell keys too ?

    gpu_nparts_ = 0;
}

// -----------------------------------------------------------------------------
//! Initialize the particle properties on devide as a mirror of the host definition
// -----------------------------------------------------------------------------
void nvidiaParticles::initializeDataOnDevice()
{
    SMILEI_ASSERT( Position.size() > 0 );

    // "Over-reserve" to minimise the cost of future reallcoation, it might be 
    // interesting to set the value to 1.2F ~~
    const auto kGrowthFactor      = 1.0F;
    const auto kPositionDimension = Position.size();
    const auto kHostParticleCount = Position[0].size();

    allocateDimensions( kPositionDimension );
    reserveParticles( static_cast<unsigned int>( static_cast<float>( kHostParticleCount ) * kGrowthFactor ) );
    allocateParticles( kHostParticleCount );

    gpu_nparts_ = kHostParticleCount;

    if( gpu_nparts_ == 0 ) {
        reserveParticles( 100 );
    } else {
        syncGPU();
    }

    if( nvidia_double_prop_.empty() ) {
        // Initialize the list of pointers
        // Done only once

        for( unsigned int i = 0; i < kPositionDimension; i++ ) {
            nvidia_double_prop_.push_back( &nvidia_position_[i] );
        }

        for( unsigned int i = 0; i < 3; i++ ) {
            nvidia_double_prop_.push_back( &nvidia_momentum_[i] );
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

//! Copy the particles from host to device
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
    // TODO(Etienne M): We are doing extra work. We could use something like
    // partition to output the invalidated particles in particles_to_move and
    // keep the good ones.

    // Manage the send data structure
    nvidiaParticles* const cp_parts                 = static_cast<nvidiaParticles*>( particles_to_move );
    const int              nparts                   = gpu_nparts_;
    const int              position_dimension_count = nvidia_position_.size();

    const int nparts_to_move = thrust::count_if( thrust::device,
                                                 std::cbegin( nvidia_cell_keys_ ),
                                                 std::cbegin( nvidia_cell_keys_ ) + nparts,
                                                 count_if_out() );

    cp_parts->gpu_nparts_ = nparts_to_move;

    // Resize it, if too small (copy_if do not resize)
    cp_parts->allocateParticles( nparts_to_move );

    // Iterator of the main data structure
    // Note: https://nvidia.github.io/thrust/api/classes/classthrust_1_1zip__iterator.html#class-thrustzip_iterator
    const auto source_iterator_first      = thrust::make_zip_iterator( thrust::make_tuple( std::begin( nvidia_position_[0] ),
                                                                                           std::begin( nvidia_momentum_[0] ),
                                                                                           std::begin( nvidia_momentum_[1] ),
                                                                                           std::begin( nvidia_momentum_[2] ),
                                                                                           std::begin( nvidia_weight_ ),
                                                                                           std::begin( nvidia_charge_ ) ) );
    const auto source_iterator_last       = source_iterator_first + nparts; // std::advance
    const auto destination_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( std::begin( cp_parts->nvidia_position_[0] ),
                                                                                           std::begin( cp_parts->nvidia_momentum_[0] ),
                                                                                           std::begin( cp_parts->nvidia_momentum_[1] ),
                                                                                           std::begin( cp_parts->nvidia_momentum_[2] ),
                                                                                           std::begin( cp_parts->nvidia_weight_ ),
                                                                                           std::begin( cp_parts->nvidia_charge_ ) ) );

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
                         std::begin( cp_parts->nvidia_position_[i] ),
                         count_if_out() );
    }

    // Special treatment for chi if radiation emission
    if( isQuantumParameter ) {
        thrust::copy_if( thrust::device,
                         std::cbegin( nvidia_chi_ ),
                         std::cbegin( nvidia_chi_ ) + nparts,
                         std::cbegin( nvidia_cell_keys_ ),
                         std::begin( cp_parts->nvidia_chi_ ),
                         count_if_out() );
    }

    if( isMonteCarlo ) {
        thrust::copy_if( thrust::device,
                         std::cbegin( nvidia_tau_ ),
                         std::cbegin( nvidia_tau_ ) + nparts,
                         std::cbegin( nvidia_cell_keys_ ),
                         std::begin( cp_parts->nvidia_tau_ ),
                         count_if_out() );
    }

    particles_to_move->syncCPU();
}

int nvidiaParticles::eraseLeavingParticles()
{
    const int position_dimension_count = nvidia_position_.size();
    const int nparts                   = gpu_nparts_;
    const int nparts_to_remove         = thrust::count_if( thrust::device,
                                                           nvidia_cell_keys_.begin(),
                                                           nvidia_cell_keys_.begin() + nparts,
                                                           count_if_out() );

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

        // Update current number of particles
        gpu_nparts_ -= nparts_to_remove;

        // Resize data structures (remove if does not resize)
        allocateParticles( gpu_nparts_ );
    }

    return /* Note: the minus matters! */ -nparts_to_remove;
}

int nvidiaParticles::injectParticles( Particles* particles_to_inject )
{
    const int nparts = gpu_nparts_;

    // Manage the recv data structure
    nvidiaParticles* const cp_parts = static_cast<nvidiaParticles*>( particles_to_inject );

    const int nparts_add = cp_parts->gpu_nparts_;
    const int tot_parts  = nparts + nparts_add;

    const int position_dimension_count = nvidia_position_.size();

    // Resize main data structure, if too small (copy_n do not resize)
    allocateParticles( tot_parts );

    const auto source_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( std::cbegin( cp_parts->nvidia_position_[0] ),
                                                                                      std::cbegin( cp_parts->nvidia_momentum_[0] ),
                                                                                      std::cbegin( cp_parts->nvidia_momentum_[1] ),
                                                                                      std::cbegin( cp_parts->nvidia_momentum_[2] ),
                                                                                      std::cbegin( cp_parts->nvidia_weight_ ),
                                                                                      std::cbegin( cp_parts->nvidia_charge_ ) ) );

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
                        std::cbegin( cp_parts->nvidia_position_[i] ),
                        nparts_add,
                        std::begin( nvidia_position_[i] ) + nparts );
    }

    if( isQuantumParameter ) {
        thrust::copy_n( thrust::device,
                        std::cbegin( cp_parts->nvidia_chi_ ),
                        nparts_add,
                        std::begin( nvidia_chi_ ) + nparts );
    }

    if( isMonteCarlo ) {
        thrust::copy_n( thrust::device,
                        std::cbegin( cp_parts->nvidia_tau_ ),
                        nparts_add,
                        std::begin( nvidia_tau_ ) + nparts );
    }

    // No more particles to move
    cp_parts->gpu_nparts_ = 0;

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
    thrust::fill( nvidia_cell_keys_.begin() + n_particles, nvidia_cell_keys_.begin() + new_size, -1 );

    gpu_nparts_ = new_size;
}

extern "C" {
    void* CreateGPUParticles( const void* parameters )
    {
        return new nvidiaParticles{ *static_cast<const Params*>( parameters ) };
    }
}
