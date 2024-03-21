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

//! Structure with specific function count_if_out for thrust::tuple operator
//! Return True if the entry is -1 as in the cell keys vector for instance
struct count_if_out
{
    constexpr __host__ __device__ bool
    operator()( const int& x ) const
    {
        return x == -1;
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

        //! Sort the particle on GPU by their cluster/cell key.
        //!
        static inline void
        sortParticleByKey( nvidiaParticles& particle_container,
                           const Params&    parameters );

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

        template <typename RandomAccessIterator0,
                  typename RandomAccessIterator1>
        static void
        doSortParticleByKey( RandomAccessIterator0 key_first,
                             RandomAccessIterator0 key_last,
                             RandomAccessIterator1 value_first );

        template <typename ClusterType,
                  typename ParticleIteratorProvider,
                  typename ParticleNoKeyIteratorProvider>
        static void
        doImportAndSortParticles( nvidiaParticles&              particle_container,
                                  nvidiaParticles&              particle_to_inject,
                                  ClusterType                   cluster_type,
                                  ParticleIteratorProvider      particle_iterator_provider,
                                  ParticleNoKeyIteratorProvider particle_no_key_iterator_provider );
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

        static void
        sortParticleByKey( nvidiaParticles& particle_container,
                           const Params&    parameters );

        static void
        importAndSortParticles( nvidiaParticles& particle_container,
                                nvidiaParticles& particle_to_inject,
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

        static void
        sortParticleByKey( nvidiaParticles& particle_container,
                           const Params&    parameters );

        static void
        importAndSortParticles( nvidiaParticles& particle_container,
                                nvidiaParticles& particle_to_inject,
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


    //! This functor assign a cluster key to a_particle.
    //!
    template <typename ClusterType>
    struct OutOfClusterPredicate
    {
    public:
    public:
        OutOfClusterPredicate( ClusterType cluster_type )
            : cluster_type_{ cluster_type }
        {
            // EMPTY
        }

        template <typename Tuple>
        __host__ __device__ bool
        operator()( const Tuple& a_particle ) const
        {
            // NOTE: its ub to set the cluster key to wrongly keyed particles
            // now..
            return thrust::get<0>( a_particle ) /* cluster key */ != cluster_type_.Index( a_particle );
        }

    protected:
        ClusterType cluster_type_;
    };


    //! If the particle's cell/cluster key is -1 it means that it needs to be
    //! evicted.
    //!
    struct OutOfBoundaryPredicate
    {
        template <typename Tuple>
        __host__ __device__ bool
        operator()( const Tuple& a_particle ) const
        {
            return thrust::get<0>( a_particle ) /* cluster key */ == -1;
        }
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
    Cluster::sortParticleByKey( nvidiaParticles& particle_container,
                                const Params&    parameters )
    {
        // This is where we do a runtime dispatch depending on the simulation's
        // dimensions.

        switch( particle_container.dimension() ) {
            case 2: {
                Cluster2D<Params::getGPUClusterWidth( 2 )>::sortParticleByKey( particle_container,
                                                                                        parameters );
                break;
            }
            case 3: {
                Cluster3D<Params::getGPUClusterWidth( 3 )>::sortParticleByKey( particle_container,
                                                                                        parameters );
                break;
            }
            default:
                // Not implemented, only Cartesian 2D or 3D for the moment
                SMILEI_ASSERT( false );
                break;
        }
    }

    // Functor to apply stride
    struct StrideFunctor {
        int stride_;
    
        __host__ __device__
        StrideFunctor(int stride) : stride_(stride) {}
    
        __host__ __device__
        int operator()(int x) const {
            return stride_ - 1 + x * stride_;
        }
    };
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
        
        //kClusterwidth should not be hard coded here
        int Ncells_per_cluster = 4;
        switch( particle_container.dimension() ) {
            case 2: {
                Ncells_per_cluster *= 4;
                break;
            }
            case 3: {
                Ncells_per_cluster *= 4*4;
                break;
            }
        }

        // Create counting iterator starting from 'zero'
        thrust::counting_iterator<int> basic_count(0);
        // Transform iterator applying the stride and changing the first element
        thrust::transform_iterator<StrideFunctor, thrust::counting_iterator<int>, int> cluster_count(basic_count, StrideFunctor(Ncells_per_cluster));
        thrust::upper_bound( thrust::device,
                             static_cast<const IDType*>( particle_container.getPtrCellKeys() ),
                             static_cast<const IDType*>( particle_container.getPtrCellKeys() ) + particle_container.deviceSize(),
                             //thrust::counting_iterator<Cluster::IDType>{ static_cast<Cluster::IDType>( 0 ) },
                             //thrust::counting_iterator<Cluster::IDType>{ static_cast<Cluster::IDType>( particle_container.last_index.size() ) },
                             cluster_count,
                             cluster_count + particle_container.last_index.size(),
                             bin_upper_bound
                             //,
                             //[](auto x, auto y)-> bool {
                             //// Divide x by 16 and compare with y
                             //return (x / 16) < y;
                             //}
 );

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
        // This is where we do a runtime dispatch depending on the simulation's
        // dimensions.

        switch( particle_container.dimension() ) {
            case 2: {
                Cluster2D<Params::getGPUClusterWidth( 2 /* 2D */ )>::importAndSortParticles( particle_container,
                                                                                             particle_to_inject,
                                                                                             parameters,
                                                                                             a_parent_patch );
                break;
            }
            case 3: {
                Cluster3D<Params::getGPUClusterWidth( 3 /* 2D */ )>::importAndSortParticles( particle_container,
                                                                                             particle_to_inject,
                                                                                             parameters,
                                                                                             a_parent_patch );
                break;
            }

            default:
                // Not implemented, only 2D for the moment
                SMILEI_ASSERT( false );
                break;
        }
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

    template <typename RandomAccessIterator0,
              typename RandomAccessIterator1>
    void
    Cluster::doSortParticleByKey( RandomAccessIterator0 key_first,
                                  RandomAccessIterator0 key_last,
                                  RandomAccessIterator1 value_first )
    {
        thrust::sort_by_key( thrust::device,
                             key_first, key_last,
                             value_first );
    }

    template <typename ClusterType,
              typename ParticleIteratorProvider,
              typename ParticleNoKeyIteratorProvider>
    void
    Cluster::doImportAndSortParticles( nvidiaParticles&              particle_container,
                                       nvidiaParticles&              particle_to_inject,
                                       ClusterType                   cluster_type,
                                       ParticleIteratorProvider      particle_iterator_provider,
                                       ParticleNoKeyIteratorProvider particle_no_key_iterator_provider )
    {
        const auto first_particle = particle_iterator_provider( particle_container );

        auto last_particle = first_particle +
                             particle_container.deviceSize(); // Obviously, we use half open ranges

        // Remove out of bound particles
        // Using more memory, we could use the faster remove_copy_if
        // NOTE: remove_if is stable.
        last_particle = thrust::remove_if( thrust::device,
                                           first_particle,
                                           last_particle,
                                           OutOfBoundaryPredicate{} );

        // Idea 1: - remove_copy_if instead of copy_if
        //         - sort(the_particles_to_inject)
        //         - merge
        //         - compute bins
        // NOTE: This method consumes a lot of memory ! O(N)

        const auto new_particle_to_inject_count  = particle_to_inject.deviceSize();
        const auto current_local_particles_count = std::distance( first_particle, last_particle );
        const auto new_particle_count            = new_particle_to_inject_count + current_local_particles_count;

        // NOTE: We really want a non-initializing vector here!
        // It's possible to give a custom allocator to thrust::device_vector.
        // Create one with construct(<>) as a noop and derive from
        // thrust::device_malloc_allocator. For now we do an explicit resize.
        particle_to_inject.softReserve( new_particle_count );
        particle_to_inject.resize( new_particle_count ); // We probably invalidated the iterators

        // Copy out of cluster/tile/chunk particles
        // partition_copy is way slower than copy_if/remove_copy_if on rocthrust
        // https://github.com/ROCmSoftwarePlatform/rocThrust/issues/247

        const auto first_particle_to_inject = particle_iterator_provider( particle_to_inject );

        // NOTE: copy_if/remove_copy_if are stable.
        const auto partitioned_particles_bounds_true  = thrust::copy_if( thrust::device,
                                                                         first_particle, last_particle,
                                                                         // Dont overwrite the particle_to_inject (at the start of the array)
                                                                         first_particle_to_inject + new_particle_to_inject_count,
                                                                         OutOfClusterPredicate<ClusterType>{ cluster_type } );
        const auto partitioned_particles_bounds_false = thrust::remove_copy_if( thrust::device,
                                                                                first_particle, last_particle,
                                                                                // Do the copy with a destination
                                                                                // starting from partitioned_particles_bounds_true
                                                                                partitioned_particles_bounds_true,
                                                                                OutOfClusterPredicate<ClusterType>{ cluster_type } );

        // Compute or recompute the cluster index of the particle_to_inject
        // NOTE:
        // - we can "save" some work here if cluster index is already computed
        // for the new particles to inject (not the one we got with copy_if).
        //
        doComputeParticleClusterKey( first_particle_to_inject,
                                     partitioned_particles_bounds_true,
                                     cluster_type );

        const auto first_particle_to_inject_no_key = particle_no_key_iterator_provider( particle_to_inject );
        const auto particle_to_rekey_count         = std::distance( first_particle_to_inject,
                                                                    partitioned_particles_bounds_true );

        doSortParticleByKey( particle_to_inject.getPtrCellKeys(),
                             particle_to_inject.getPtrCellKeys() + particle_to_rekey_count,
                             first_particle_to_inject_no_key );

        // This free generates a lot of memory fragmentation.
        // particle_container.free();
        // Same as for particle_to_inject, non-initializing vector is best.
        particle_container.softReserve( new_particle_count );
        particle_container.resize( new_particle_count );

        // Merge by key
        // NOTE: Dont merge in place on GPU. That means we need an other large buffer!
        //
        thrust::merge_by_key( thrust::device,
                              particle_to_inject.getPtrCellKeys(),                           // Input range 1, first key
                              particle_to_inject.getPtrCellKeys() + particle_to_rekey_count, // Input range 1, last key
                              particle_to_inject.getPtrCellKeys() + particle_to_rekey_count, // Input range 2, first key
                              particle_to_inject.getPtrCellKeys() + new_particle_count,      // Input range 2, last key
                              first_particle_to_inject_no_key,                               // Input range 1, first value
                              first_particle_to_inject_no_key + particle_to_rekey_count,     // Input range 2, first value
                              particle_container.getPtrCellKeys(),                           // Output range first key
                              particle_no_key_iterator_provider( particle_container ) );     // Output range first value

        // Recompute bins
        computeBinIndex( particle_container );

        // This free generates a lot of memory fragmentation. If we enable it we
        // reduce significantly the memory usage over time but a memory spike
        // will still be present. Unfortunately, this free generates soo much
        // fragmentation (like the one above) that at some point the GPU memory
        // allocator will fail!
        // particle_to_inject.free();
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
    // arrange cell index to have a proper sort for collisions
    {
        //Sort is done per cell for collisions but arranged in larger tiles for projector performance.

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
        // We assume same kClusterWidth in all dimensions
        //static constexpr SizeType x_cluster_dimension_in_cell = kClusterWidth;
        //static constexpr SizeType y_cluster_dimension_in_cell = kClusterWidth;
        static constexpr SizeType cluster_size_in_cell = kClusterWidth*kClusterWidth;

        //const SizeType local_x_particle_cluster_coordinate_in_cluster = local_x_particle_coordinate_in_cell / x_cluster_dimension_in_cell;
        const SizeType local_x_particle_cluster_coordinate_in_cluster = local_x_particle_coordinate_in_cell / kClusterWidth;
        //const SizeType local_y_particle_cluster_coordinate_in_cluster = local_y_particle_coordinate_in_cell / y_cluster_dimension_in_cell;
        const SizeType local_y_particle_cluster_coordinate_in_cluster = local_y_particle_coordinate_in_cell / kClusterWidth;

        //const SizeType y_stride = local_y_dimension_in_cluster_;

        // The indexing order is: x * ywidth * zwidth + y * zwidth + z
        //const SizeType cluster_index = local_x_particle_cluster_coordinate_in_cluster * y_stride +
        //                               local_y_particle_cluster_coordinate_in_cluster;
        const SizeType cluster_index = (local_x_particle_cluster_coordinate_in_cluster * local_y_dimension_in_cluster_
                                       + local_y_particle_cluster_coordinate_in_cluster ) * cluster_size_in_cell
                                       + kClusterWidth * (local_x_particle_coordinate_in_cell % kClusterWidth) 
                                       + local_y_particle_coordinate_in_cell % kClusterWidth;

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
        static constexpr SizeType cluster_size_in_cell = kClusterWidth*kClusterWidth*kClusterWidth;

        const SizeType local_x_particle_cluster_coordinate_in_cluster = local_x_particle_coordinate_in_cell / x_cluster_dimension_in_cell;
        const SizeType local_y_particle_cluster_coordinate_in_cluster = local_y_particle_coordinate_in_cell / y_cluster_dimension_in_cell;
        const SizeType local_z_particle_cluster_coordinate_in_cluster = local_z_particle_coordinate_in_cell / z_cluster_dimension_in_cell;

        const SizeType y_stride = local_y_dimension_in_cluster_;
        const SizeType z_stride = local_z_dimension_in_cluster_;

        // The indexing order is: x * ywidth * zwidth + y * zwidth + z
        //const SizeType cluster_index = local_x_particle_cluster_coordinate_in_cluster * z_stride * y_stride +
        //                               local_y_particle_cluster_coordinate_in_cluster * z_stride +
        //                               local_z_particle_cluster_coordinate_in_cluster;
        const SizeType cluster_index = ( local_x_particle_cluster_coordinate_in_cluster * local_y_dimension_in_cluster_ * local_z_dimension_in_cluster_
                                       + local_y_particle_cluster_coordinate_in_cluster * local_z_dimension_in_cluster_ 
                                       + local_z_particle_cluster_coordinate_in_cluster) * cluster_size_in_cell
                                       + kClusterWidth * kClusterWidth * (local_x_particle_coordinate_in_cell % kClusterWidth) 
                                       +                 kClusterWidth * (local_y_particle_coordinate_in_cell % kClusterWidth)
                                       +                                  local_z_particle_coordinate_in_cell % kClusterWidth;
        // It is not the cluster index anymore. The name of this variable should be changed
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

    template <Cluster::DifferenceType kClusterWidth>
    void
    Cluster2D<kClusterWidth>::sortParticleByKey( nvidiaParticles& particle_container,
                                                 const Params& )
    {
        // This is where we do a runtime dispatch depending on the simulation's
        // qed/radiation settings.

        // NOTE: For now we support dont support qed/radiations. Performance
        // comes from specialization.

        // TODO(Etienne M): Find a better way to dispatch at runtime. This is
        // complex to read and to maintain.

        if( particle_container.has_quantum_parameter ) {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                SMILEI_ASSERT( false );
            }
        } else {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                // The appropriate thrust::zip_iterator for the current
                // simulation's parameters

                const auto value_first = thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrPosition( 0 ),
                                                                                        particle_container.getPtrPosition( 1 ),
                                                                                        particle_container.getPtrMomentum( 0 ),
                                                                                        particle_container.getPtrMomentum( 1 ),
                                                                                        particle_container.getPtrMomentum( 2 ),
                                                                                        particle_container.getPtrWeight(),
                                                                                        particle_container.getPtrCharge() ) );

                doSortParticleByKey( particle_container.getPtrCellKeys(),
                                     particle_container.getPtrCellKeys() + particle_container.deviceSize(),
                                     value_first );
            }
        }
    }

    template <Cluster::DifferenceType kClusterWidth>
    void
    Cluster3D<kClusterWidth>::sortParticleByKey( nvidiaParticles& particle_container,
                                                 const Params& )
    {
        // This is where we do a runtime dispatch depending on the simulation's
        // qed/radiation settings.

        // NOTE: For now we support dont support qed/radiations. Performance
        // comes from specialization.

        // TODO(Etienne M): Find a better way to dispatch at runtime. This is
        // complex to read and to maintain.

        if( particle_container.has_quantum_parameter ) {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                SMILEI_ASSERT( false );
            }
        } else {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                // The appropriate thrust::zip_iterator for the current
                // simulation's parameters

                if (particle_container.tracked) {
                    const auto value_first = thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrPosition( 0 ),
                                                                                            particle_container.getPtrPosition( 1 ),
                                                                                            particle_container.getPtrPosition( 2 ),
                                                                                            particle_container.getPtrMomentum( 0 ),
                                                                                            particle_container.getPtrMomentum( 1 ),
                                                                                            particle_container.getPtrMomentum( 2 ),
                                                                                            particle_container.getPtrWeight(),
                                                                                            particle_container.getPtrCharge(),
                                                                                            particle_container.getPtrId() ) );
                    doSortParticleByKey( particle_container.getPtrCellKeys(),
                                         particle_container.getPtrCellKeys() + particle_container.deviceSize(),
                                         value_first );

                }
                else {
                    const auto value_first = thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrPosition( 0 ),
                                                                                            particle_container.getPtrPosition( 1 ),
                                                                                            particle_container.getPtrPosition( 2 ),
                                                                                            particle_container.getPtrMomentum( 0 ),
                                                                                            particle_container.getPtrMomentum( 1 ),
                                                                                            particle_container.getPtrMomentum( 2 ),
                                                                                            particle_container.getPtrWeight(),
                                                                                            particle_container.getPtrCharge() ) );
                    doSortParticleByKey( particle_container.getPtrCellKeys(),
                                         particle_container.getPtrCellKeys() + particle_container.deviceSize(),
                                         value_first );
                }
            }
        }
    }

    template <Cluster::DifferenceType kClusterWidth>
    void
    Cluster2D<kClusterWidth>::importAndSortParticles( nvidiaParticles& particle_container,
                                                      nvidiaParticles& particle_to_inject,
                                                      const Params&    parameters,
                                                      const Patch&     a_parent_patch )
    {
        // This is where we do a runtime dispatch depending on the simulation's
        // qed/radiation settings.

        // NOTE: For now we support dont support qed/radiations. Performance
        // comes from specialization.

        // TODO(Etienne M): Find a better way to dispatch at runtime. This is
        // complex to read and to maintain.
        int CellStartingGlobalIndex_for_x = a_parent_patch.getCellStartingGlobalIndex_noGC(0);
        int CellStartingGlobalIndex_for_y = a_parent_patch.getCellStartingGlobalIndex_noGC(1);

        const Cluster2D cluster_manipulator{ parameters.res_space[0],
                                             parameters.res_space[1],
                                             parameters.patch_size_[0],
                                             parameters.patch_size_[1],
                                             CellStartingGlobalIndex_for_x, CellStartingGlobalIndex_for_y};

        if( particle_container.has_quantum_parameter ) {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                SMILEI_ASSERT( false );
            }
        } else {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                // Returns the appropriate thrust::zip_iterator for the
                // current simulation's parameters
                const auto particle_iterator_provider = []( nvidiaParticles& particle_container ) {
                    return thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrCellKeys(),
                                                                          particle_container.getPtrPosition( 0 ),
                                                                          particle_container.getPtrPosition( 1 ),
                                                                          particle_container.getPtrMomentum( 0 ),
                                                                          particle_container.getPtrMomentum( 1 ),
                                                                          particle_container.getPtrMomentum( 2 ),
                                                                          particle_container.getPtrWeight(),
                                                                          particle_container.getPtrCharge() ) );
                };

                const auto particle_no_key_iterator_provider = []( nvidiaParticles& particle_container ) {
                    return thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrPosition( 0 ),
                                                                          particle_container.getPtrPosition( 1 ),
                                                                          particle_container.getPtrMomentum( 0 ),
                                                                          particle_container.getPtrMomentum( 1 ),
                                                                          particle_container.getPtrMomentum( 2 ),
                                                                          particle_container.getPtrWeight(),
                                                                          particle_container.getPtrCharge() ) );
                };

                doImportAndSortParticles( particle_container,
                                          particle_to_inject,
                                          cluster_manipulator,
                                          particle_iterator_provider,
                                          particle_no_key_iterator_provider );
            }
        }
    }

    template <Cluster::DifferenceType kClusterWidth>
    void
    Cluster3D<kClusterWidth>::importAndSortParticles( nvidiaParticles& particle_container,
                                                      nvidiaParticles& particle_to_inject,
                                                      const Params&    parameters,
                                                      const Patch&     a_parent_patch )
    {
        // This is where we do a runtime dispatch depending on the simulation's
        // qed/radiation settings.

        // NOTE: For now we support dont support qed/radiations. Performance
        // comes from specialization.

        // TODO(Etienne M): Find a better way to dispatch at runtime. This is
        // complex to read and to maintain.
        int CellStartingGlobalIndex_for_x = a_parent_patch.getCellStartingGlobalIndex_noGC(0);
        int CellStartingGlobalIndex_for_y = a_parent_patch.getCellStartingGlobalIndex_noGC(1);
        int CellStartingGlobalIndex_for_z = a_parent_patch.getCellStartingGlobalIndex_noGC(2);

        const Cluster3D cluster_manipulator{ parameters.res_space[0],
                                             parameters.res_space[1],
                                             parameters.res_space[2],
                                             parameters.patch_size_[0],
                                             parameters.patch_size_[1],
                                             parameters.patch_size_[2],
                                             CellStartingGlobalIndex_for_x, 
                                             CellStartingGlobalIndex_for_y, CellStartingGlobalIndex_for_z};

        if( particle_container.has_quantum_parameter ) {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                SMILEI_ASSERT( false );
            }
        } else {
            if( particle_container.has_Monte_Carlo_process ) {
                SMILEI_ASSERT( false );
            } else {
                // Returns the appropriate thrust::zip_iterator for the
                // current simulation's parameters
                if (particle_container.tracked) {
                    const auto particle_iterator_provider = []( nvidiaParticles& particle_container ) {
                        return thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrCellKeys(),
                                                                          particle_container.getPtrPosition( 0 ),
                                                                          particle_container.getPtrPosition( 1 ),
                                                                          particle_container.getPtrPosition( 2 ),
                                                                          particle_container.getPtrMomentum( 0 ),
                                                                          particle_container.getPtrMomentum( 1 ),
                                                                          particle_container.getPtrMomentum( 2 ),
                                                                          particle_container.getPtrWeight(),
                                                                          particle_container.getPtrCharge(),
                                                                          particle_container.getPtrId() ) );
                    };
                    const auto particle_no_key_iterator_provider = []( nvidiaParticles& particle_container ) {
                        return thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrPosition( 0 ),
                                                                          particle_container.getPtrPosition( 1 ),
                                                                          particle_container.getPtrPosition( 2 ),
                                                                          particle_container.getPtrMomentum( 0 ),
                                                                          particle_container.getPtrMomentum( 1 ),
                                                                          particle_container.getPtrMomentum( 2 ),
                                                                          particle_container.getPtrWeight(),
                                                                          particle_container.getPtrCharge(),
                                                                          particle_container.getPtrId() ) );
                    };
                    doImportAndSortParticles( particle_container,
                                          particle_to_inject,
                                          cluster_manipulator,
                                          particle_iterator_provider,
                                          particle_no_key_iterator_provider );
                }
                else {
                    const auto particle_iterator_provider = []( nvidiaParticles& particle_container ) {
                        return thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrCellKeys(),
                                                                          particle_container.getPtrPosition( 0 ),
                                                                          particle_container.getPtrPosition( 1 ),
                                                                          particle_container.getPtrPosition( 2 ),
                                                                          particle_container.getPtrMomentum( 0 ),
                                                                          particle_container.getPtrMomentum( 1 ),
                                                                          particle_container.getPtrMomentum( 2 ),
                                                                          particle_container.getPtrWeight(),
                                                                          particle_container.getPtrCharge() ) );
                    };

                    const auto particle_no_key_iterator_provider = []( nvidiaParticles& particle_container ) {
                        return thrust::make_zip_iterator( thrust::make_tuple( particle_container.getPtrPosition( 0 ),
                                                                          particle_container.getPtrPosition( 1 ),
                                                                          particle_container.getPtrPosition( 2 ),
                                                                          particle_container.getPtrMomentum( 0 ),
                                                                          particle_container.getPtrMomentum( 1 ),
                                                                          particle_container.getPtrMomentum( 2 ),
                                                                          particle_container.getPtrWeight(),
                                                                          particle_container.getPtrCharge() ) );
                    };

                    doImportAndSortParticles( particle_container,
                                          particle_to_inject,
                                          cluster_manipulator,
                                          particle_iterator_provider,
                                          particle_no_key_iterator_provider );
                }
            }
        }
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
        detail::Cluster::sortParticleByKey( *this, *parameters_ );

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
void nvidiaParticles::copyFromDeviceToHost()
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
}

unsigned int nvidiaParticles::deviceCapacity() const
{
    SMILEI_ASSERT( nvidia_momentum_.size() >= 1 );
    // Could be any particle component that we know will be used in any case.
    return nvidia_momentum_[0].capacity();
}

// -----------------------------------------------------------------------------
//! Extract particles from the Particles object and put
//! them in the Particles object `particles_to_move`
// -----------------------------------------------------------------------------
void nvidiaParticles::extractParticles( Particles* particles_to_move )
{
    // TODO(Etienne M): We are doing extra work. We could use something like
    // std::partition to output the invalidated particles in particles_to_move
    // and keep the good ones. This would help us avoid the std::remove_if in
    // the particle injection and sorting algorithm.

    // Manage the send data structure
    nvidiaParticles* const cp_parts                 = static_cast<nvidiaParticles*>( particles_to_move );
    const int              nparts                   = gpu_nparts_;
    const int              position_dimension_count = nvidia_position_.size();

    const int nparts_to_move = thrust::count_if( thrust::device,
                                                 nvidia_cell_keys_.cbegin(),
                                                 nvidia_cell_keys_.cbegin() + nparts,
                                                 count_if_out() );

    // Resize it, if too small (copy_if do not resize)
    cp_parts->resize( nparts_to_move );

    // Iterator of the main data structure
    // NOTE: https://nvidia.github.io/thrust/api/classes/classthrust_1_1zip__iterator.html#class-thrustzip_iterator
    const auto source_iterator_first      = thrust::make_zip_iterator( thrust::make_tuple( nvidia_position_[0].begin(),
                                                                                           nvidia_momentum_[0].begin(),
                                                                                           nvidia_momentum_[1].begin(),
                                                                                           nvidia_momentum_[2].begin(),
                                                                                           nvidia_weight_.begin(),
                                                                                           nvidia_charge_.begin() ) );
    const auto source_iterator_last       = source_iterator_first + nparts; // std::advance
    const auto destination_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( cp_parts->nvidia_position_[0].begin(),
                                                                                           cp_parts->nvidia_momentum_[0].begin(),
                                                                                           cp_parts->nvidia_momentum_[1].begin(),
                                                                                           cp_parts->nvidia_momentum_[2].begin(),
                                                                                           cp_parts->nvidia_weight_.begin(),
                                                                                           cp_parts->nvidia_charge_.begin() ) );

    // Copy send particles in dedicated data structure if nvidia_cell_keys_=0 (currently = 1 if keeped, new PartBoundCond::apply(...))
    thrust::copy_if( thrust::device,
                     source_iterator_first,
                     source_iterator_last,
                     // Copy depending on count_if_out()(nvidia_cell_keys_[i])
                     nvidia_cell_keys_.cbegin(),
                     destination_iterator_first,
                     count_if_out() );

    // Copy the other position values depending on the simulation's grid
    // dimensions
    for( int i = 1; i < position_dimension_count; ++i ) {
        thrust::copy_if( thrust::device,
                         nvidia_position_[i].cbegin(),
                         nvidia_position_[i].cbegin() + nparts,
                         nvidia_cell_keys_.cbegin(),
                         cp_parts->nvidia_position_[i].begin(),
                         count_if_out() );
    }

    // Special treatment for chi if radiation emission
    if( has_quantum_parameter ) {
        thrust::copy_if( thrust::device,
                         nvidia_chi_.cbegin(),
                         nvidia_chi_.cbegin() + nparts,
                         nvidia_cell_keys_.cbegin(),
                         cp_parts->nvidia_chi_.begin(),
                         count_if_out() );
    }

    if( has_Monte_Carlo_process ) {
        thrust::copy_if( thrust::device,
                         nvidia_tau_.cbegin(),
                         nvidia_tau_.cbegin() + nparts,
                         nvidia_cell_keys_.cbegin(),
                         cp_parts->nvidia_tau_.begin(),
                         count_if_out() );
    }

    if( tracked ) {
        thrust::copy_if( thrust::device,
                         nvidia_id_.cbegin(),
                         nvidia_id_.cbegin() + nparts,
                         nvidia_cell_keys_.cbegin(),
                         cp_parts->nvidia_id_.begin(),
                         count_if_out() );
    }

    particles_to_move->copyFromDeviceToHost();
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
//                           count_if_out() );
//    }
//
//}

// -----------------------------------------------------------------------------
//! Erase particles leaving the patch object on device
// -----------------------------------------------------------------------------
int nvidiaParticles::eraseLeavingParticles()
{
    const int position_dimension_count = nvidia_position_.size();
    const int nparts                   = gpu_nparts_;
    const int nparts_to_remove         = thrust::count_if( thrust::device,
                                                           nvidia_cell_keys_.begin(),
                                                           nvidia_cell_keys_.begin() + nparts,
                                                           count_if_out() );


    if( nparts_to_remove > 0 ) {
        const auto first_particle = thrust::make_zip_iterator( thrust::make_tuple( nvidia_position_[0].begin(),
                                                                                   nvidia_momentum_[0].begin(),
                                                                                   nvidia_momentum_[1].begin(),
                                                                                   nvidia_momentum_[2].begin(),
                                                                                   nvidia_weight_.begin(),
                                                                                   nvidia_charge_.begin() ) );

        const auto last_particle = first_particle + nparts;

        // Remove particles which leaves current patch
        thrust::remove_if( thrust::device,
                           first_particle,
                           last_particle,
                           nvidia_cell_keys_.cbegin(),
                           count_if_out() );

        // Remove the other position values depending on the simulation's grid
        // dimensions
        for( int i = 1; i < position_dimension_count; ++i ) {
            thrust::remove_if( thrust::device,
                               nvidia_position_[i].begin(),
                               nvidia_position_[i].begin() + nparts,
                               nvidia_cell_keys_.cbegin(),
                               count_if_out() );
        }

        if( has_quantum_parameter ) {
            thrust::remove_if( thrust::device,
                               nvidia_chi_.begin(),
                               nvidia_chi_.begin() + nparts,
                               nvidia_cell_keys_.cbegin(),
                               count_if_out() );
        }

        if( has_Monte_Carlo_process ) {
            thrust::remove_if( thrust::device,
                               nvidia_tau_.begin(),
                               nvidia_tau_.begin() + nparts,
                               nvidia_cell_keys_.cbegin(),
                               count_if_out() );
        }

        if( tracked ) {
            thrust::remove_if( thrust::device,
                               nvidia_id_.begin(),
                               nvidia_id_.begin() + nparts,
                               nvidia_cell_keys_.cbegin(),
                               count_if_out() );
        }

        // Update current number of particles
        gpu_nparts_ -= nparts_to_remove;

        // Resize data structures (remove_if does not resize)
        resize( gpu_nparts_ );
    }

    return nparts_to_remove;
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
    resize( tot_parts );

    const auto source_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( cp_parts->nvidia_position_[0].cbegin(),
                                                                                      cp_parts->nvidia_momentum_[0].cbegin(),
                                                                                      cp_parts->nvidia_momentum_[1].cbegin(),
                                                                                      cp_parts->nvidia_momentum_[2].cbegin(),
                                                                                      cp_parts->nvidia_weight_.cbegin(),
                                                                                      cp_parts->nvidia_charge_.cbegin() ) );

    // Iterator of the main data structure (once it has been resized)
    const auto destination_iterator_first = thrust::make_zip_iterator( thrust::make_tuple( nvidia_position_[0].begin(),
                                                                                           nvidia_momentum_[0].begin(),
                                                                                           nvidia_momentum_[1].begin(),
                                                                                           nvidia_momentum_[2].begin(),
                                                                                           nvidia_weight_.begin(),
                                                                                           nvidia_charge_.begin() ) ) +
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
                        cp_parts->nvidia_position_[i].cbegin(),
                        nparts_add,
                        nvidia_position_[i].begin() + nparts );
    }

    if( has_quantum_parameter ) {
        thrust::copy_n( thrust::device,
                        cp_parts->nvidia_chi_.cbegin(),
                        nparts_add,
                        nvidia_chi_.begin() + nparts );
    }

    if( has_Monte_Carlo_process ) {
        thrust::copy_n( thrust::device,
                        cp_parts->nvidia_tau_.cbegin(),
                        nparts_add,
                        nvidia_tau_.begin() + nparts );
    }

    if( tracked ) {
        thrust::copy_n( thrust::device,
                        cp_parts->nvidia_id_.cbegin(),
                        nparts_add,
                        nvidia_id_.begin() + nparts );
    }

    // No more particles to move
    cp_parts->resize( 0 );

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

int nvidiaParticles::prepareBinIndex()
{
    if( first_index.size() == 0 ) {
        // Some Particles object like particles_to_move do not have allocated
        // bins, we skip theses.
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
    injectParticles( particles_to_inject );
}

extern "C"
{
    void* CreateGPUParticles( const void* parameters, const void* a_parent_patch )
    {
        return new nvidiaParticles{ *static_cast<const Params*>( parameters ),
                                    *static_cast<const Patch*>( a_parent_patch ) };
    }
}
