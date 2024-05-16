// -----------------------------------------------------------------------------
//
//! \file nvidiaParticles.cu
//
//! \brief contains the nvidiaParticles class methods
//! Extension of the Class Particles for GPU
//
// -----------------------------------------------------------------------------

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/count.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/gather.h>


#include "Patch.h"
#include "gpu.h"
#include "nvidiaParticles.h"

// TODO(Etienne M): The makefile does not recognise this file and doesn't compute
// it's dependencies. If you make a modification in one of the header this file
// includes, you must `touch` this file. IF you dont do that you'll have ABI/ODR
// issues (!).

// Language: "in cell" means the number of cells for that, conversely, in cluster means 
// the number of clusters as a unit of length, etc.

////////////////////////////////////////////////////////////////////////////////
// Cell key manipulation functor definition
////////////////////////////////////////////////////////////////////////////////

//! Predicate for cell_keys
//! Return True if the entry is equal to `code`
template<int code>
struct cellKeyEquals
{
    constexpr __host__ __device__ bool
    operator()( const int& x ) const
    {
        return x == code;
    }
};

template<int key>
struct cellKeyBelow
{
    constexpr __host__ __device__ bool
    operator()( const int& x ) const
    {
        return x < key;
    }
};

namespace detail {

    ////////////////////////////////////////////////////////////////////////////////
    // Cluster manipulation functor definition
    ////////////////////////////////////////////////////////////////////////////////

    //! Cluster manipulation functionalities common to all dimension.
    //! NOTE: This only focus on GPU data manipulation. The host data shall
    //! not be handled here !
    //!
    struct Cluster
    {
    public:
        //! Same type as what is used in nvidia_cell_keys_
        //!
        using IDType         = int;
        using SizeType       = unsigned int;
        using DifferenceType = int;

    public:
        //! Compute the cell key for all the particles (not only a subset).
        //!
        static inline void
        computeParticleClusterKey( nvidiaParticles& particle_container,
                                   const Params&    parameters,
                                   const Patch&     a_parent_patch );

        //! precondition:
        //!     - nvidia_cell_keys_ shall be sorted in non decreasing order
        //!     - last_index.data() is a pointer mapped to GPU via
        //!       HostDeviceMemoryManagement
        //!
        static inline void
        computeBinIndex( nvidiaParticles& particle_container );

        //! Sorting by cluster and binning
        //!
        //! precondition:
        //!     - particle_container is already sorted by cluster or
        //!       particle_container is not sorted anymore (after a push) but
        //!       still contains the old cluster key untouched.
        //!       PartBoundCond::apply will set the keys to zero !
        //!
        static inline void
        importAndSortParticles( nvidiaParticles& particle_container,
                                nvidiaParticles& particle_to_inject,
                                const Params&    parameters,
                                const Patch&     a_parent_patch );

    protected:
        template <typename InputIterator,
                  typename ClusterType>
        static void
        doComputeParticleClusterKey( InputIterator first,
                                     InputIterator last,
                                     ClusterType   cluster_type );

    };


    template <Cluster::DifferenceType kClusterWidth>
    struct Cluster2D : public Cluster
    {
    public:
    public:
        Cluster2D( double   inverse_x_cell_dimension,
                   double   inverse_y_cell_dimension,
                   SizeType local_x_dimension_in_cell,
                   SizeType local_y_dimension_in_cell,
                   int CellStartingGlobalIndex_for_x,
                   int CellStartingGlobalIndex_for_y);

        //! Compute the cell key of a_particle. a_particle shall be a tuple (from a
        //! zipiterator).
        //! The first value of a_particle is the cell key value, the other values are
        //! the positions x and y.
        //!
        template <typename Tuple>
        __host__ __device__ IDType
        Index( const Tuple& a_particle ) const;

        //! Compute the cell key of a particle range.
        //!
        static void
        computeParticleClusterKey( nvidiaParticles& particle_container,
                                   const Params&    parameters,
                                   const Patch&     a_parent_patch );

    public:
        double   inverse_of_x_cell_dimension_;
        double   inverse_of_y_cell_dimension_;
        SizeType local_y_dimension_in_cluster_;
        int CellStartingGlobalIndex_for_x_;
        int CellStartingGlobalIndex_for_y_;
    };

    template <Cluster::DifferenceType kClusterWidth>
    struct Cluster3D : public Cluster
    {
    public:
    public:
        Cluster3D( double   inverse_x_cell_dimension,
                   double   inverse_y_cell_dimension,
                   double   inverse_z_cell_dimension,
                   SizeType local_x_dimension_in_cell,
                   SizeType local_y_dimension_in_cell,
                   SizeType local_z_dimension_in_cell,
                   int CellStartingGlobalIndex_for_x,
                   int CellStartingGlobalIndex_for_y,
                   int CellStartingGlobalIndex_for_z);

        //! Compute the cell key of a_particle. a_particle shall be a tuple (from a
        //! zipiterator).
        //! The first value of a_particle is the cell key value, the other values are
        //! the positions x and y.
        //!
        template <typename Tuple>
        __host__ __device__ IDType
        Index( const Tuple& a_particle ) const;

        //! Compute the cell key of a particle range.
        //!
        static void
        computeParticleClusterKey( nvidiaParticles& particle_container,
                                   const Params&    parameters,
                                   const Patch&     a_parent_patch );

    public:
        double   inverse_of_x_cell_dimension_;
        double   inverse_of_y_cell_dimension_;
        double   inverse_of_z_cell_dimension_;
        SizeType local_y_dimension_in_cluster_;
        SizeType local_z_dimension_in_cluster_;
        int CellStartingGlobalIndex_for_x_;
        int CellStartingGlobalIndex_for_y_;
        int CellStartingGlobalIndex_for_z_;
    };


    //! This functor assign a cluster key to a_particle.
    //!
    template <typename ClusterType>
    class AssignClusterIndex
    {
    public:
    public:
        AssignClusterIndex( ClusterType cluster_type )
            : cluster_type_{ cluster_type }
        {
            // EMPTY
        }

        template <typename Tuple>
        __host__ __device__ void
        operator()( Tuple& a_particle ) const
        {
            thrust::get<0>( a_particle ) /* cluster key */ = cluster_type_.Index( a_particle );
        }

    protected:
        ClusterType cluster_type_;
    };


    ////////////////////////////////////////////////////////////////////////////////
    // Cluster manipulation functor method definitions
    ////////////////////////////////////////////////////////////////////////////////

    inline void
    Cluster::computeParticleClusterKey( nvidiaParticles& particle_container,
                                        const Params&    parameters,
                                        const Patch&     a_parent_patch )
    {
        // This is where we do a runtime dispatch depending on the simulation's
        // dimensions.

        switch( particle_container.dimension() ) {
            case 2: {
                Cluster2D<Params::getGPUClusterWidth( 2 )>::computeParticleClusterKey( particle_container,
                                                                                                parameters,
                                                                                                a_parent_patch );
                break;
            }
            case 3: {
                Cluster3D<Params::getGPUClusterWidth( 3 )>::computeParticleClusterKey( particle_container,
                                                                                                parameters,
                                                                                                a_parent_patch );
                break;
            }
            default:
                // Not implemented, only Cartesian 2D or 3D for the moment
                SMILEI_ASSERT( false );
                break;
        }
    }

    inline void
    Cluster::computeBinIndex( nvidiaParticles& particle_container )
    {
        SMILEI_GPU_ASSERT_MEMORY_IS_ON_DEVICE( particle_container.last_index.data() );

        Cluster::IDType* bin_upper_bound = smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( particle_container.last_index.data() );

        // SMILEI_ASSERT( thrust::is_sorted( thrust::device,
        //                                   static_cast<const IDType*>( particle_container.getPtrCellKeys() ),
        //                                   static_cast<const IDType*>( particle_container.getPtrCellKeys() ) + particle_container.deviceSize() ) );

        // NOTE: On some benchmark, I found this upper_bound usage faster than the counting_iterator (by a lot(!) ~x3, but
        // it's so fast anyway..)

        // thrust::upper_bound( thrust::device,
        //                      nvidia_cell_keys_.cbegin(), nvidia_cell_keys_.cend(),
        //                      key_bound_to_search.cbegin(), key_bound_to_search.cend(),
        //                      bin_upper_bound );

        // NOTE: A particle is in a bin if the index of the bin is the same integer value as the particle's cell key.
        // The particles are sorted by cell key. We can do a simple binary search to find the upper bound of a bin.
        //
        thrust::upper_bound( thrust::device,
                             static_cast<const IDType*>( particle_container.getPtrCellKeys() ),
                             static_cast<const IDType*>( particle_container.getPtrCellKeys() ) + particle_container.deviceSize(),
                             thrust::counting_iterator<Cluster::IDType>{ static_cast<Cluster::IDType>( 0 ) },
                             thrust::counting_iterator<Cluster::IDType>{ static_cast<Cluster::IDType>( particle_container.last_index.size() ) },
                             bin_upper_bound );

        // SMILEI_ASSERT( thrust::is_sorted( thrust::device,
        //                                   bin_upper_bound,
        //                                   bin_upper_bound + particle_container.last_index.size() ) );
    }

    inline void
    Cluster::importAndSortParticles( nvidiaParticles& particle_container,
                                     nvidiaParticles& particle_to_inject,
                                     const Params&    parameters,
                                     const Patch&     a_parent_patch )
    {
        // Remove out of bound particles
        const auto erased_count = particle_container.eraseParticlesByPredicate( cellKeyBelow<0>() );
        
        const auto initial_count = particle_container.deviceSize() - erased_count;
        const auto inject_count  = particle_to_inject.deviceSize();
        const auto new_count     = initial_count + inject_count;
        
        // Resize particles
        // NOTE: We really want a non-initializing vector here!
        // It's possible to give a custom allocator to thrust::device_vector.
        // Create one with construct(<>) as a noop and derive from
        // thrust::device_malloc_allocator. For now we do an explicit resize.
        particle_container.softReserve( new_count );
        particle_container.resize( new_count );
        
        // Combine imported particles to main particles
        particle_container.copyParticles( &particle_to_inject, initial_count );
        
        // Compute keys of particles
        computeParticleClusterKey( particle_container, parameters, a_parent_patch );
        
        // Use particle_to_inject as a buffer
        particle_to_inject.softReserve( new_count );
        particle_to_inject.resize( new_count );
        
        // Sort particles using thrust::gather, according to the sorting map
        particle_container.sortParticleByKey( particle_to_inject );
        
        // Recompute bins
        computeBinIndex( particle_container );
        
        // This free generates a lot of memory fragmentation. If we enable it we
        // reduce significantly the memory usage over time but a memory spike
        // will still be present. Unfortunately, this free generates soo much
        // fragmentation (like the one above) that at some point the GPU memory
        // allocator will fail!
        // particle_to_inject.free();
    }

    template <typename InputIterator,
              typename ClusterType>
    void
    Cluster::doComputeParticleClusterKey( InputIterator first,
                                          InputIterator last,
                                          ClusterType   cluster_type )
    {
        thrust::for_each( thrust::device,
                          first, last,
                          AssignClusterIndex<ClusterType>{ cluster_type } );
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Cluster2D method definitions
    ////////////////////////////////////////////////////////////////////////////////

    template <Cluster::DifferenceType kClusterWidth>
    Cluster2D<kClusterWidth>::Cluster2D( double   inverse_x_cell_dimension,
                                         double   inverse_y_cell_dimension,
                                         SizeType local_x_dimension_in_cell,
                                         SizeType local_y_dimension_in_cell,
                                         int CellStartingGlobalIndex_for_x, int CellStartingGlobalIndex_for_y )
        : inverse_of_x_cell_dimension_{ inverse_x_cell_dimension }
        , inverse_of_y_cell_dimension_{ inverse_y_cell_dimension }
        , local_y_dimension_in_cluster_{ local_y_dimension_in_cell / kClusterWidth }
        , CellStartingGlobalIndex_for_x_{CellStartingGlobalIndex_for_x}
        , CellStartingGlobalIndex_for_y_{CellStartingGlobalIndex_for_y}
    {
        // EMPTY
    }

    template <Cluster::DifferenceType kClusterWidth>
    Cluster3D<kClusterWidth>::Cluster3D( double   inverse_x_cell_dimension,
                                         double   inverse_y_cell_dimension,
                                         double   inverse_z_cell_dimension,
                                         SizeType local_x_dimension_in_cell,
                                         SizeType local_y_dimension_in_cell,
                                         SizeType local_z_dimension_in_cell,
                                         int CellStartingGlobalIndex_for_x,
                                         int CellStartingGlobalIndex_for_y, int CellStartingGlobalIndex_for_z )
        : inverse_of_x_cell_dimension_{ inverse_x_cell_dimension }
        , inverse_of_y_cell_dimension_{ inverse_y_cell_dimension }
        , inverse_of_z_cell_dimension_{ inverse_z_cell_dimension }
        , local_y_dimension_in_cluster_{ local_y_dimension_in_cell / kClusterWidth }
        , local_z_dimension_in_cluster_{ local_z_dimension_in_cell / kClusterWidth }
        , CellStartingGlobalIndex_for_x_{CellStartingGlobalIndex_for_x}
        , CellStartingGlobalIndex_for_y_{CellStartingGlobalIndex_for_y}
        , CellStartingGlobalIndex_for_z_{CellStartingGlobalIndex_for_z}
    {
        // EMPTY
    }

    template <Cluster::DifferenceType kClusterWidth>
    template <typename Tuple>
    __host__ __device__ typename Cluster2D<kClusterWidth>::IDType
    Cluster2D<kClusterWidth>::Index( const Tuple& a_particle ) const
    {
        const SizeType local_x_particle_coordinate_in_cell = static_cast<SizeType>( thrust::get<1>( a_particle ) *
                                                                                    inverse_of_x_cell_dimension_ ) -
                                                             CellStartingGlobalIndex_for_x_;
        const SizeType local_y_particle_coordinate_in_cell = static_cast<SizeType>( thrust::get<2>( a_particle ) *
                                                                                    inverse_of_y_cell_dimension_ ) -
                                                             CellStartingGlobalIndex_for_y_;

        // These divisions will be optimized.
        // The integer division rounding behavior is expected.

        // NOTE: Flat tiles have been studied but were not as efficient for the
        // projection. The square provides the minimal perimeter (and thus ghost
        // cell amount) for a given area.
        static constexpr SizeType x_cluster_dimension_in_cell = kClusterWidth;
        static constexpr SizeType y_cluster_dimension_in_cell = kClusterWidth;

        const SizeType local_x_particle_cluster_coordinate_in_cluster = local_x_particle_coordinate_in_cell / x_cluster_dimension_in_cell;
        const SizeType local_y_particle_cluster_coordinate_in_cluster = local_y_particle_coordinate_in_cell / y_cluster_dimension_in_cell;

        const SizeType y_stride = local_y_dimension_in_cluster_;

        // The indexing order is: x * ywidth * zwidth + y * zwidth + z
        const SizeType cluster_index = local_x_particle_cluster_coordinate_in_cluster * y_stride +
                                       local_y_particle_cluster_coordinate_in_cluster;

        return static_cast<IDType>( cluster_index );
    }
    
    template <Cluster::DifferenceType kClusterWidth>
    template <typename Tuple>
    __host__ __device__ typename Cluster3D<kClusterWidth>::IDType
    Cluster3D<kClusterWidth>::Index( const Tuple& a_particle ) const
    {
        const SizeType local_x_particle_coordinate_in_cell = static_cast<SizeType>( thrust::get<1>( a_particle ) *
                                                                                    inverse_of_x_cell_dimension_ ) -
                                                             CellStartingGlobalIndex_for_x_;
        const SizeType local_y_particle_coordinate_in_cell = static_cast<SizeType>( thrust::get<2>( a_particle ) *
                                                                                    inverse_of_y_cell_dimension_ ) -
                                                             CellStartingGlobalIndex_for_y_;
        const SizeType local_z_particle_coordinate_in_cell = static_cast<SizeType>( thrust::get<3>( a_particle ) *
                                                                                    inverse_of_z_cell_dimension_ ) -
                                                             CellStartingGlobalIndex_for_z_;

        // These divisions will be optimized.
        // The integer division rounding behavior is expected.

        // NOTE: Flat tiles have been studied but were not as efficient for the
        // projection. The square provides the minimal perimeter (and thus ghost
        // cell amount) for a given area.
        static constexpr SizeType x_cluster_dimension_in_cell = kClusterWidth;
        static constexpr SizeType y_cluster_dimension_in_cell = kClusterWidth;
        static constexpr SizeType z_cluster_dimension_in_cell = kClusterWidth;

        const SizeType local_x_particle_cluster_coordinate_in_cluster = local_x_particle_coordinate_in_cell / x_cluster_dimension_in_cell;
        const SizeType local_y_particle_cluster_coordinate_in_cluster = local_y_particle_coordinate_in_cell / y_cluster_dimension_in_cell;
        const SizeType local_z_particle_cluster_coordinate_in_cluster = local_z_particle_coordinate_in_cell / z_cluster_dimension_in_cell;

        const SizeType y_stride = local_y_dimension_in_cluster_;
        const SizeType z_stride = local_z_dimension_in_cluster_;

        // The indexing order is: x * ywidth * zwidth + y * zwidth + z
        const SizeType cluster_index = local_x_particle_cluster_coordinate_in_cluster * z_stride * y_stride +
                                       local_y_particle_cluster_coordinate_in_cluster * z_stride +
                                       local_z_particle_cluster_coordinate_in_cluster;

        return static_cast<IDType>( cluster_index );
    }

    template <Cluster::DifferenceType kClusterWidth>
    void
    Cluster2D<kClusterWidth>::computeParticleClusterKey( nvidiaParticles& particle_container,
                                                         const Params&    parameters,
                                                         const Patch&     a_parent_patch )
    {
        const auto first = thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrCellKeys(),
                                                                          static_cast<const double*>( particle_container.getPtrPosition( 0 ) ),
                                                                          static_cast<const double*>( particle_container.getPtrPosition( 1 ) ) ) );
        const auto last  = first + particle_container.deviceSize();
        int CellStartingGlobalIndex_for_x = a_parent_patch.getCellStartingGlobalIndex_noGC(0);
        int CellStartingGlobalIndex_for_y = a_parent_patch.getCellStartingGlobalIndex_noGC(1);
        doComputeParticleClusterKey( first, last,
                                     Cluster2D<Params::getGPUClusterWidth( 2 )>{ parameters.res_space[0],
                                                                                          parameters.res_space[1],
                                                                                          parameters.patch_size_[0],
                                                                                          parameters.patch_size_[1],
                                                                                          CellStartingGlobalIndex_for_x,
                                                                                          CellStartingGlobalIndex_for_y } );
    }

    template <Cluster::DifferenceType kClusterWidth>
    void
    Cluster3D<kClusterWidth>::computeParticleClusterKey( nvidiaParticles& particle_container,
                                                         const Params&    parameters,
                                                         const Patch&     a_parent_patch )
    {
        const auto first = thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrCellKeys(),
                                                                          static_cast<const double*>( particle_container.getPtrPosition( 0 ) ),
                                                                          static_cast<const double*>( particle_container.getPtrPosition( 1 ) ),
                                                                          static_cast<const double*>( particle_container.getPtrPosition( 2 ) ) ) );
        const auto last  = first + particle_container.deviceSize();
        int CellStartingGlobalIndex_for_x = a_parent_patch.getCellStartingGlobalIndex_noGC(0);
        int CellStartingGlobalIndex_for_y = a_parent_patch.getCellStartingGlobalIndex_noGC(1);
        int CellStartingGlobalIndex_for_z = a_parent_patch.getCellStartingGlobalIndex_noGC(2);
        doComputeParticleClusterKey( first, last,
                                     Cluster3D<Params::getGPUClusterWidth( 3 )>{ parameters.res_space[0],
                                                                                          parameters.res_space[1],
                                                                                          parameters.res_space[2],
                                                                                          parameters.patch_size_[0],
                                                                                          parameters.patch_size_[1],
                                                                                          parameters.patch_size_[2],
                                                                                          CellStartingGlobalIndex_for_x,
                                                                                          CellStartingGlobalIndex_for_y,
                                                                                          CellStartingGlobalIndex_for_z } );
    }

} // namespace detail


////////////////////////////////////////////////////////////////////////////////
// nvidiaParticles method definitions
////////////////////////////////////////////////////////////////////////////////

nvidiaParticles::nvidiaParticles( const Params& parameters,
                                  const Patch&  a_parent_patch )
    : Particles{}
    , parameters_{ &parameters }
    , parent_patch_{ &a_parent_patch }
    , gpu_nparts_{}
{
    // EMPTY
}

nvidiaParticles::~nvidiaParticles() {
    // Manage last_index if allocated on GPU
    if (smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( last_index.data() )) {
        smilei::tools::gpu::HostDeviceMemoryManagement::DeviceFree( last_index );
    }
}

void nvidiaParticles::resizeDimensions( unsigned int nDim )
{
    nvidia_position_.resize( nDim );
    nvidia_momentum_.resize( 3 );
}

void nvidiaParticles::softReserve( unsigned int particle_count, float growth_factor  )
{
    if( particle_count <= deviceCapacity() ) {
        // Dont reserve, for now we have enough capacity.
        return;
    }

    const unsigned int new_capacity = static_cast<unsigned int>( particle_count * growth_factor );

    for( unsigned int idim = 0; idim < nvidia_position_.size(); idim++ ) {
        nvidia_position_[idim].reserve( new_capacity );
    }

    for( unsigned int idim = 0; idim < 3; idim++ ) {
        nvidia_momentum_[idim].reserve( new_capacity );
    }

    nvidia_weight_.reserve( new_capacity );
    nvidia_charge_.reserve( new_capacity );

    if( has_quantum_parameter ) {
        nvidia_chi_.reserve( new_capacity );
    }

    if( has_Monte_Carlo_process ) {
        nvidia_tau_.reserve( new_capacity );
    }

    if( tracked ) {
        nvidia_id_.reserve( new_capacity );
    }

    nvidia_cell_keys_.reserve( new_capacity );
}

void nvidiaParticles::reserve( unsigned int particle_count )
{
    for( unsigned int idim = 0; idim < nvidia_position_.size(); idim++ ) {
        nvidia_position_[idim].reserve( particle_count );
    }

    for( unsigned int idim = 0; idim < 3; idim++ ) {
        nvidia_momentum_[idim].reserve( particle_count );
    }

    nvidia_weight_.reserve( particle_count );
    nvidia_charge_.reserve( particle_count );

    if( has_quantum_parameter ) {
        nvidia_chi_.reserve( particle_count );
    }

    if( has_Monte_Carlo_process ) {
        nvidia_tau_.reserve( particle_count );
    }

    if( tracked ) {
        nvidia_id_.reserve( particle_count );
    }

    nvidia_cell_keys_.reserve( particle_count );
}

void nvidiaParticles::resize( unsigned int particle_count )
{

    // TODO(Etienne M): Use non-initializing vector/allocator (dont pay the cost
    // of what you dont use) ?

    for( int idim = 0; idim < nvidia_position_.size(); idim++ ) {
        nvidia_position_[idim].resize( particle_count );
    }

    for( int idim = 0; idim < 3; idim++ ) {
        nvidia_momentum_[idim].resize( particle_count );
    }

    nvidia_weight_.resize( particle_count );
    nvidia_charge_.resize( particle_count );

    if( has_quantum_parameter ) {
        nvidia_chi_.resize( particle_count );
    }

    if( has_Monte_Carlo_process ) {
        nvidia_tau_.resize( particle_count );
    }

    if( tracked ) {
        nvidia_id_.resize( particle_count );
    }

    nvidia_cell_keys_.resize( particle_count );

    gpu_nparts_ = particle_count;
}

void nvidiaParticles::free()
{
    for( auto& a_vector : nvidia_position_ ) {
        thrust::device_vector<double> a_dummy_vector{};
        std::swap( a_vector, a_dummy_vector );
    }

    for( auto& a_vector : nvidia_momentum_ ) {
        thrust::device_vector<double> a_dummy_vector{};
        std::swap( a_vector, a_dummy_vector );
    }

    {
        thrust::device_vector<double> a_dummy_vector{};
        std::swap( nvidia_weight_, a_dummy_vector );
    }

    {
        thrust::device_vector<short> a_dummy_vector{};
        std::swap( nvidia_charge_, a_dummy_vector );
    }

    if( has_quantum_parameter ) {
        thrust::device_vector<double> a_dummy_vector{};
        std::swap( nvidia_chi_, a_dummy_vector );
    }

    if( has_Monte_Carlo_process ) {
        thrust::device_vector<double> a_dummy_vector{};
        std::swap( nvidia_tau_, a_dummy_vector );
    }

    if( tracked ) {
        thrust::device_vector<uint64_t> a_dummy_vector{};
        std::swap( nvidia_id_, a_dummy_vector );
    }

    {
        thrust::device_vector<int> a_dummy_vector{};
        std::swap( nvidia_cell_keys_, a_dummy_vector );
    }

    gpu_nparts_ = 0;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Resize particle vectors
// ---------------------------------------------------------------------------------------------------------------------
void nvidiaParticles::deviceResize( unsigned int new_size )
{
    for( unsigned int iprop=0 ; iprop<nvidia_double_prop_.size() ; iprop++ ) {
        ( *nvidia_double_prop_[iprop] ).resize(new_size);
    }

    for( unsigned int iprop=0 ; iprop<nvidia_short_prop_.size() ; iprop++ ) {
        ( *nvidia_short_prop_[iprop] ).resize(new_size);
    }

    //
    // for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
    //     ( *nvidia_uint64_prop[iprop] ).resize( n_particles+n_additional_particles );
    // }

    if (tracked) {
        nvidia_id_.resize( new_size );
    }

    nvidia_cell_keys_.resize( new_size );

    gpu_nparts_ = new_size;
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

    if (tracked) {
        nvidia_id_.clear();
    }

    gpu_nparts_ = 0;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Reset cell_keys to default value
// ---------------------------------------------------------------------------------------------------------------------
void nvidiaParticles::resetCellKeys(void)
{
    thrust::fill(nvidia_cell_keys_.begin(), nvidia_cell_keys_.begin() + gpu_nparts_, -1);
}

// -----------------------------------------------------------------------------
//! Initialize the particle properties on device as a mirror of the host definition
// -----------------------------------------------------------------------------
void nvidiaParticles::initializeDataOnDevice()
{
    SMILEI_ASSERT( Position.size() > 0 );
    // The world shall end if we call this function multiple times
    SMILEI_ASSERT( nvidia_double_prop_.empty() );

    const auto kPositionDimension = Position.size();

    // We sure that we have as many say, position dimension as the base class.
    resizeDimensions( kPositionDimension );

    // Initialize the list of pointers

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
    if( has_quantum_parameter ) {
        nvidia_double_prop_.push_back( &nvidia_chi_ );
    }

    // Optical Depth for Monte-Carlo processes:
    // - if the discontinuous (Monte-Carlo) radiation reaction
    // is activated, tau is the incremental optical depth to emission
    if( has_Monte_Carlo_process ) {
        nvidia_double_prop_.push_back( &nvidia_tau_ );
    }

    const auto kHostParticleCount = Position[0].size();

    if( kHostParticleCount == 0 ) {
        // Should we reserve some space ?
        // reserve( 100 );
    } else {
        copyFromHostToDevice();
    }

    if( prepareBinIndex() < 0 ) {
        // Either we deal with a simulation with unsupported space dimensions
        // (1D/AM) or we are not using OpenMP or we are dealing with particle
        // object without allocated bin (particle_to_move for instance).
        // We'll use the old, naive, unsorted particles injection
        // implementation.

        // Dont call setHostBinIndex. For particle that have binning this is a
        // redundant call. But for the particle that should not get binned
        // (ie: particle_to_move) , this is a bug (!) and will trigger an
        // assertion.

        // setHostBinIndex();
    } else {

        // At this point, a copy of the host particles and last_index is on the
        // device and we know we support the space dimension.

        detail::Cluster::computeParticleClusterKey( *this, *parameters_, *parent_patch_ );

        // The particles are not correctly sorted when created.
        sortParticleByKey();

        detail::Cluster::computeBinIndex( *this );
        setHostBinIndex();
    }
}

// -------------------------------------------------------------------------------------------------
//! Copy particle IDs from host to device
// -------------------------------------------------------------------------------------------------
void nvidiaParticles::initializeIDsOnDevice()
{
    nvidia_id_.resize( Id.size() );
    thrust::copy((Id).begin(), (Id).end(), (nvidia_id_).begin());
}

// -------------------------------------------------------------------------------------------------
//! Copy the particles from host to device
// -------------------------------------------------------------------------------------------------
void nvidiaParticles::copyFromHostToDevice()
{
    resize( Position[0].size() );

    for( int idim = 0; idim < Position.size(); idim++ ) {
        thrust::copy( Position[idim].begin(), Position[idim].end(), nvidia_position_[idim].begin() );
    }

    for( int idim = 0; idim < Momentum.size(); idim++ ) {
        thrust::copy( Momentum[idim].begin(), Momentum[idim].end(), nvidia_momentum_[idim].begin() );
    }

    thrust::copy( Weight.begin(), Weight.end(), nvidia_weight_.begin() );

    thrust::copy( Charge.begin(), Charge.end(), nvidia_charge_.begin() );

    if( has_quantum_parameter ) {
        thrust::copy( Chi.begin(), Chi.end(), nvidia_chi_.begin() );
    }

    if( has_Monte_Carlo_process ) {
        thrust::copy( Tau.begin(), Tau.end(), nvidia_tau_.begin() );
    }

    if( tracked ) {
        thrust::copy( Id.begin(), Id.end(), nvidia_id_.begin() );
    }
}

// -------------------------------------------------------------------------------------------------
//! Copy device to host
// -------------------------------------------------------------------------------------------------
void nvidiaParticles::copyFromDeviceToHost( bool copy_keys )
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
    if (has_quantum_parameter) {
        Chi.resize( gpu_nparts_ );
        thrust::copy((nvidia_chi_).begin(), (nvidia_chi_).begin()+gpu_nparts_, (Chi).begin());
    }
    if (has_Monte_Carlo_process) {
        Tau.resize( gpu_nparts_ );
        thrust::copy((nvidia_tau_).begin(), (nvidia_tau_).begin()+gpu_nparts_, (Tau).begin());
    }
    if (tracked) {
        Id.resize( gpu_nparts_ );
        thrust::copy((nvidia_id_).begin(), (nvidia_id_).begin()+gpu_nparts_, (Id).begin());
    }
    if (copy_keys) {
        cell_keys.resize( gpu_nparts_ );
        thrust::copy((nvidia_cell_keys_).begin(), (nvidia_cell_keys_).begin()+gpu_nparts_, (cell_keys).begin());
    }
}

unsigned int nvidiaParticles::deviceCapacity() const
{
    SMILEI_ASSERT( nvidia_momentum_.size() >= 1 );
    // Could be any particle component that we know will be used in any case.
    return nvidia_momentum_[0].capacity();
}

// -----------------------------------------------------------------------------
//! Move leaving particles to the buffer
// -----------------------------------------------------------------------------
void nvidiaParticles::copyLeavingParticlesToBuffer( Particles* buffer )
{
    copyParticlesByPredicate( buffer, cellKeyBelow<-1>() );
    buffer->copyFromDeviceToHost( true );
}


//! Copy particles which statisfy some predicate
template<typename Predicate>
void nvidiaParticles::copyParticlesByPredicate( Particles* buffer, Predicate pred )
{
    // TODO(Etienne M): We are doing extra work. We could use something like
    // std::partition to output the invalidated particles in buffer
    // and keep the good ones. This would help us avoid the std::remove_if in
    // the particle injection and sorting algorithm.
    
    // Count particles satisfying the predicate
    const auto keys = getPtrCellKeys();
    const int nparts_to_copy = thrust::count_if( thrust::device, keys, keys + gpu_nparts_, pred );
    
    // Resize destination buffer (copy_if does not resize)
    nvidiaParticles* const dest = static_cast<nvidiaParticles*>( buffer );
    dest->resize( nparts_to_copy );
    
    if( nparts_to_copy ) {
        // Copy the particles to the destination
        for( int ip = 0; ip < getNDoubleProp(); ip++ ) {
            const auto in = getPtrDoubleProp( ip );
            const auto out = dest->getPtrDoubleProp( ip );
            thrust::copy_if( thrust::cuda::par_nosync, in, in + gpu_nparts_, keys, out, pred );
        }
        for( int ip = 0; ip < getNShortProp(); ip++ ) {
            const auto in = getPtrShortProp( ip );
            const auto out = dest->getPtrShortProp( ip );
            thrust::copy_if( thrust::cuda::par_nosync, in, in + gpu_nparts_, keys, out, pred );
        }
        if( tracked ) {
            const auto in = getPtrId();
            const auto out = dest->getPtrId();
            thrust::copy_if( thrust::cuda::par_nosync, in, in + gpu_nparts_, keys, out, pred );
        }
        cudaDeviceSynchronize();
    }
}

void nvidiaParticles::copyParticles( Particles* particles_to_inject )
{
    const auto nparts = gpu_nparts_;
    nvidiaParticles* to_inject = static_cast<nvidiaParticles*>( particles_to_inject );
    resize( nparts + to_inject->gpu_nparts_ );
    copyParticles( to_inject, nparts );
}

void nvidiaParticles::copyParticles( nvidiaParticles* particles_to_inject, size_t offset )
{
    // Copy the particles to the destination
    for( int ip = 0; ip < getNDoubleProp(); ip++ ) {
        const auto in = particles_to_inject->getPtrDoubleProp( ip );
        const auto out = getPtrDoubleProp( ip );
        thrust::copy_n( thrust::cuda::par_nosync, in, particles_to_inject->gpu_nparts_, out + offset );
    }
    for( int ip = 0; ip < getNShortProp(); ip++ ) {
        const auto in = particles_to_inject->getPtrShortProp( ip );
        const auto out = getPtrShortProp( ip );
        thrust::copy_n( thrust::cuda::par_nosync, in, particles_to_inject->gpu_nparts_, out + offset );
    }
    if( tracked ) {
        const auto in = particles_to_inject->getPtrId();
        const auto out = getPtrId();
        thrust::copy_n( thrust::cuda::par_nosync, in, particles_to_inject->gpu_nparts_, out + offset );
    }
    cudaDeviceSynchronize();
}

// -----------------------------------------------------------------------------
//! Erase `npart` particles from `ipart`
// -----------------------------------------------------------------------------
//void nvidiaParticles::eraseParticleOnDevice(int ipart, int npart) {
//
//    const auto first_particle = thrust::make_zip_iterator( thrust::make_tuple( std::begin( nvidia_position_[0] ),
//                                                                               std::begin( nvidia_momentum_[0] ),
//                                                                               std::begin( nvidia_momentum_[1] ),
//                                                                               std::begin( nvidia_momentum_[2] ),
//                                                                               std::begin( nvidia_weight_ ),
//                                                                               std::begin( nvidia_charge_ ) ) );
//
//    // Remove the other position values depending on the simulation's grid
//    // dimensions
//    for( int i = 1; i < position_dimension_count; ++i ) {
//        thrust::remove_if( thrust::device,
//                           std::begin( nvidia_position_[i] ),
//                           std::begin( nvidia_position_[i] ) + nparts,
//                           std::cbegin( nvidia_cell_keys_ ),
//                           cellKeyEquals<-1>() );
//    }
//
//}

// -----------------------------------------------------------------------------
//! Erase particles leaving the patch on device
// -----------------------------------------------------------------------------
int nvidiaParticles::eraseLeavingParticles()
{
    const auto nremoved = eraseParticlesByPredicate( cellKeyBelow<0>() );
    resize( gpu_nparts_ - nremoved );
    return nremoved;
}

//! "Erase" particles but does not resize the arrays!
template<typename Predicate>
int nvidiaParticles::eraseParticlesByPredicate( Predicate pred )
{
    const auto keys = getPtrCellKeys();
    const int nparts_to_remove = thrust::count_if( thrust::device, keys, keys + gpu_nparts_, pred );
    
    // Copy the particles to the destination
    // Using more memory, we could use the faster remove_copy_if
    // NOTE: remove_if is stable.
    for( int ip = 0; ip < getNDoubleProp(); ip++ ) {
        const auto in = getPtrDoubleProp( ip );
        thrust::remove_if( thrust::cuda::par_nosync, in, in + gpu_nparts_, keys, pred );
    }
    for( int ip = 0; ip < getNShortProp(); ip++ ) {
        const auto in = getPtrShortProp( ip );
        thrust::remove_if( thrust::cuda::par_nosync, in, in + gpu_nparts_, keys, pred );
    }
    if( tracked ) {
        const auto in = getPtrId();
        thrust::remove_if( thrust::cuda::par_nosync, in, in + gpu_nparts_, keys, pred );
    }
    cudaDeviceSynchronize();
    
    return nparts_to_remove;
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

    // for( unsigned int iprop=0 ; iprop<uint64_prop.size() ; iprop++ ) {
    //     ( *nvidia_uint64_prop[iprop] ).resize( n_particles+n_additional_particles );
    // }

    if (tracked) {
        nvidia_id_.resize( new_size );
        thrust::fill( nvidia_id_.begin() + n_particles, nvidia_id_.begin() + new_size, 0 );
    }

    nvidia_cell_keys_.resize( new_size );
    thrust::fill( nvidia_cell_keys_.begin() + n_particles, nvidia_cell_keys_.begin() + new_size, -1 );

    gpu_nparts_ = new_size;
}

//! Import Particles and sort depending if Binning is available or not
void nvidiaParticles::importAndSortParticles( Particles* particles_to_inject )
{
    if( parameters_->isGPUParticleBinningAvailable() ) {
        detail::Cluster::importAndSortParticles( *static_cast<nvidiaParticles*>( this ),
                                                 *static_cast<nvidiaParticles*>( particles_to_inject ),
                                                 *parameters_,
                                                 *parent_patch_ );
    } else {
        // When GPU particle binning is not supported, fallback to a naive implementation
        naiveImportAndSortParticles( static_cast<nvidiaParticles*>( particles_to_inject ) );
    }

    setHostBinIndex();
}

//! Sort by cell_keys_
//! This version synchronizes for every vector, but uses less buffers
void nvidiaParticles::sortParticleByKey()
{
    // Make a sorting map using the cell keys (like numpy.argsort)
    thrust::device_vector<int> index( gpu_nparts_ );
    thrust::sequence( thrust::device, index.begin(), index.end() );
    thrust::sort_by_key( thrust::device, nvidia_cell_keys_.begin(), nvidia_cell_keys_.end(), index.begin() );
    
    // Sort particles using thrust::gather, according to the sorting map
    thrust::device_vector<double> buffer( gpu_nparts_ );
    for( int ip = 0; ip < getNDoubleProp(); ip++ ) {
        thrust::gather( thrust::device, index.begin(), index.end(), getPtrDoubleProp( ip ), buffer.begin() );
        swapDoubleProp( ip, buffer );
    }
    buffer.clear();
    thrust::device_vector<short> buffer_short( gpu_nparts_ );
    for( int ip = 0; ip < getNShortProp(); ip++ ) {
        thrust::gather( thrust::device, index.begin(), index.end(), getPtrShortProp( ip ), buffer_short.begin() );
        swapShortProp( ip, buffer_short );
    }
    buffer_short.clear();
    if( tracked ) {
        thrust::device_vector<uint64_t> buffer_uint64( gpu_nparts_ );
        thrust::gather( thrust::device, index.begin(), index.end(), getPtrId(), buffer_uint64.begin() );
        swapId( buffer_uint64 );
        buffer_uint64.clear();
    }
}

//! Sort by cell_keys_
//! This version is asynchronous, but requires a buffer of equal size to be provided
void nvidiaParticles::sortParticleByKey( nvidiaParticles& buffer )
{
    // Make a sorting map using the cell keys (like numpy.argsort)
    thrust::device_vector<int> index( gpu_nparts_ );
    thrust::sequence( thrust::device, index.begin(), index.end() );
    thrust::sort_by_key( thrust::device, nvidia_cell_keys_.begin(), nvidia_cell_keys_.end(), index.begin() );
    
    // Sort particles using thrust::gather, according to the sorting map
    for( int ip = 0; ip < getNDoubleProp(); ip++ ) {
        thrust::gather( thrust::cuda::par_nosync, index.begin(), index.end(), getPtrDoubleProp( ip ), buffer.getPtrDoubleProp( ip ) );
    }
    for( int ip = 0; ip < getNShortProp(); ip++ ) {
        thrust::gather( thrust::cuda::par_nosync, index.begin(), index.end(), getPtrShortProp( ip ), buffer.getPtrShortProp( ip ) );
    }
    if( tracked ) {
        thrust::gather( thrust::cuda::par_nosync, index.begin(), index.end(), getPtrId(), buffer.getPtrId() );
    }
    cudaDeviceSynchronize();
    
    swap( buffer );
}

int nvidiaParticles::prepareBinIndex()
{
    if( first_index.size() == 0 ) {
        // Some Particles object do not have allocated bins, we skip theses.
        return -1;
    }

    const int kGPUBinCount = parameters_->getGPUBinCount();

    if( kGPUBinCount < 0 ) {
        // Unsupported space dimension or the offloading technology is not
        // supported, dont do GPU binning.
        return -1;
    }

    // We completely ignore/discard/overwrite what's done in
    // ParticleCreator::create regarding binning.
    // NOTE: maybe ParticleCreator::create should not be doing the particle
    // binning and should only be responsible for particle initialization (pos,
    // momentum etc.).
    // We are forced to deal with first_index even though its completely
    // redundant as long as the bins are dense (no holes).

    const auto particle_count = last_index.back();

    first_index.resize( 1 );
    last_index.resize( kGPUBinCount );

    // By definition it should be zero, so this is a redundant assignment
    first_index.back() = 0;

    // Dont try to allocate 2 times, even if it's harmless, that would be a bug!
    SMILEI_ASSERT( !smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( last_index.data() ) );

    // We'll need last_index to be on the GPU.

    // TODO(Etienne M): FREE. If we have load balancing or other patch
    // creation/destruction available (which is not the case on GPU ATM),
    // we should be taking care of freeing this GPU memory.
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( last_index );

    return 0;
}

void nvidiaParticles::setHostBinIndex()
{
    // TODO(Etienne M): You may want to inject, create etc. into a non binned
    // nvidiaParticles object (without allocated first/last_index). For now, we
    // assert it does not happen. I think a fix only requires:
    //  if( last_index.empty() ) { return; }
    //
    SMILEI_ASSERT( !last_index.empty() );

    last_index.back() = deviceSize();
    last_index[0]     = last_index.back();
}

void nvidiaParticles::naiveImportAndSortParticles( nvidiaParticles* particles_to_inject )
{
    // Erase particles that leaves this patch
    eraseLeavingParticles();

    // Inject newly arrived particles in particles_to_inject
    const size_t current_size = gpu_nparts_;
    resize( current_size + particles_to_inject->size() );
    copyParticles( particles_to_inject, current_size );
    particles_to_inject->clear();
}

extern "C"
{
    void* CreateGPUParticles( const void* parameters, const void* a_parent_patch )
    {
        return new nvidiaParticles{ *static_cast<const Params*>( parameters ),
                                    *static_cast<const Patch*>( a_parent_patch ) };
    }
}
