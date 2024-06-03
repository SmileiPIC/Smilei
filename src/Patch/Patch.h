#ifndef PATCH_H
#define PATCH_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <limits.h>

#include "Random.h"
#include "Params.h"
#include "SmileiMPI.h"
#include "PartWall.h"
#include "ParticleInjector.h"
#include "Interpolator.h"
#include "Projector.h"

// Fine timer ids
#define interpolation_timer_id_     0
#define push_timer_id_              1
#define projection_timer_id_        2
#define cell_keys_timer_id_         3
#define ionization_timer_id_        4
#define radiation_timer_id_         5
#define mBW_timer_id_               6
#define interp_fields_env_timer_id_ 7

class DomainDecomposition;
class Diagnostic;
class SimWindow;
class BinaryProcesses;

//! Class Patch :
//!   - data container
//!   - sub MPI domain + MPI methods
//! Collection of patch = MPI domain
class Patch
{
    friend class SmileiMPI;
    friend class VectorPatch;
    friend class SimWindow;
    friend class SyncVectorPatch;
    friend class AsyncMPIbuffers;
public:
    //! Constructor for Patch
    Patch( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch );
    //! Cloning Constructor for Patch
    Patch( Patch *patch, Params &params, SmileiMPI *smpi, unsigned int ipatch );
    
    //! First initialization step for patches
    void initStep1( Params &params );
    //! Second initialization step for patches
    virtual void initStep2( Params &params, DomainDecomposition *domain_decomposition ) = 0;
    //! Third initialization step for patches
    void initStep3( Params &params, SmileiMPI *smpi, unsigned int n_moved );
    //! Last creation step
    void finishCreation( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition );
    //! Last cloning step
    void finishCloning( Patch *patch, Params &params, SmileiMPI *smpi, unsigned int n_moved, bool with_particles );
    
    //! Finalize MPI environment : especially requests array for non blocking communications
    void finalizeMPIenvironment( Params &params );
    
    void setLocationAndAllocateFields( Params &params, DomainDecomposition *domain_decomposition, VectorPatch &vecPatch );
   
    //! Destructor for Patch
    virtual ~Patch();
    
    // Main PIC objects : data & operators
    // -----------------------------------
    
    //! Species, Particles, of the current Patch
    std::vector<Species *> vecSpecies;
    //! Electromagnetic fields and densities (E, B, J, rho) of the current Patch
    ElectroMagn *EMfields;
    
    //! Optional internal boundary condifion on Particles
    PartWalls *partWalls;
    //! Optional binary processes operators
    std::vector<BinaryProcesses *> vecBPs;
    
    //! Injectors of the current patch
    std::vector<ParticleInjector *> particle_injector_vector_;
    
    //! "fake" particles for the probe diagnostics
    std::vector<ProbeParticles *> probes;
    
    //! Classical interpolator for the probe diagnostic only
    Interpolator *probesInterp;
    
    
    // Geometrical description
    // -----------------------
    
    //!Hilbert index of the patch. Number of the patch along the Hilbert curve.
    unsigned int hindex;
    
    //!Cartesian coordinates of the patch. X,Y,Z of the Patch according to its Hilbert index.
    std::vector<unsigned int> Pcoordinates;
    
    std::vector<unsigned int> size_;
    std::vector<unsigned int> oversize;
    
    // Detailed timers (at the patch level)
    // -----------------------

    // Initialize timers
    // 0 - Interpolation
    // 1 - Pusher
    // 2 - Projection
    // 3 - exchange init + cell_keys
    // 4 - ionization
    // 5 - radiation
    // 6 - Breit-Wheeler
    // 7 - Interp Fields_Env
    // 8 - Proj Susceptibility
    // 9 - Push Momentum
    // 10 - Interp Env_Old
    // 11 - Proj Currents
    // 12 - Push Pos
    // 13 - Sorting

#ifdef  __DETAILED_TIMERS

    // OpenMP properties
    // -----------------------
    
    int number_of_threads_;
    
    // Detailed timers
    // -----------------------
    
    //! Timers for the patch
    std::vector<double> patch_timers_;

    //! temporary timers
    std::vector<double> patch_tmp_timers_;

#endif

#ifdef __DETAILED_TIMERS
    inline void __attribute__((always_inline)) startFineTimer(unsigned int index) {
#ifdef _OMPTASKS
        const int ithread = Tools::getOMPThreadNum();
        patch_tmp_timers_[index * number_of_threads_ + ithread] = MPI_Wtime();
#else
        patch_tmp_timers_[index] = MPI_Wtime();
#endif
#else
    inline void __attribute__((always_inline)) startFineTimer(unsigned int) {
#endif
    }
    
#ifdef  __DETAILED_TIMERS
    inline void __attribute__((always_inline)) stopFineTimer(unsigned int index) {
#ifdef _OMPTASKS
        const int ithread = Tools::getOMPThreadNum();   
        patch_timers_[index * number_of_threads_ + ithread] += MPI_Wtime() - patch_tmp_timers_[index * number_of_threads_ + ithread];
#else
        patch_timers_[index] += MPI_Wtime() - patch_tmp_timers_[index];
#endif
#else
    inline void __attribute__((always_inline)) stopFineTimer(unsigned int) {
#endif
    }

    // Random number generator.
    Random * rand_;
    
    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------
    
    //! Clean the MPI buffers for communications
    void cleanMPIBuffers( int ispec, Params &params );
    //! manage Idx of particles per direction,
    void copyExchParticlesToBuffers( int ispec, Params &params );
    //! init comm  nbr of particles
    void exchNbrOfParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! finalize comm / nbr of particles, init exch / particles
    void endNbrOfParticles( int ispec, int iDim );
    //! extract particles from main data structure to buffers, init exch / particles
    void prepareParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! effective exchange of particles
    void exchParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! finalize exch / particles
    void waitExchParticles( int ispec, int iDim );
    //! Treat diagonalParticles
    void cornersParticles( int ispec, Params &params, int iDim );
    //! inject particles received in main data structure and particles sorting
    void importAndSortParticles( int ispec, Params &params );
    //! clean memory resizing particles structure
    void cleanParticlesOverhead( Params &params );
    //! delete Particles included in the index of particles to exchange. Assumes indexes are sorted.
    void cleanupSentParticles( int ispec, std::vector<int> *indexes_of_particles_to_exchange );

#ifdef SMILEI_ACCELERATOR_GPU
    //! Allocate and copy all the field grids on device
    void allocateAndCopyFieldsOnDevice();

    //! Allocate all field grids on device
    void allocateFieldsOnDevice();

    //! Copy All field grids from device to host
    void copyFieldsFromDeviceToHost();

    //! Copy All fields from host to device
    void copyFieldsFromHostToDevice();

    //! Deallocate field grids on device
    void deleteFieldsOnDevice();
#endif

    //! init comm / sum densities
    virtual void initSumField( Field *field, int iDim, SmileiMPI *smpi, bool devPtr = false );
    //! init comm / sum densities
    virtual void initSumFieldComplex( Field *, int, SmileiMPI * ) {};
    //! finalize comm / sum densities
    virtual void finalizeSumField( Field *field, int iDim );
    
    //! init comm / exchange fields in direction iDim only
    virtual void initExchange( Field *field, int iDim, SmileiMPI *smpi, bool devPtr = false );
    //! init comm / exchange complex fields in direction iDim only
    virtual void initExchangeComplex( Field *field, int iDim, SmileiMPI *smpi );
    //! finalize comm / exchange fields
    virtual void finalizeExchange( Field *field, int iDim );
    
    virtual void exchangeField_movewin ( Field* field, int clrw ) = 0;
    
    // Create MPI_Datatype to exchange fields
    virtual void createType2( Params &params ) = 0;
    virtual void cleanType() = 0;
    
    // Geometrical methods
    // --------------------
    
    //! Return the hibert index of current patch
    inline unsigned int Hindex()
    {
        return  hindex;
    }
    
    //! Method to identify the rank 0 MPI process
    inline bool isMaster()
    {
        return ( hindex==0 );
    }
    
    //! Should be pure virtual, see child classes
    inline bool isXmin()
    {
        return isBoundary( 0, 0 );
    }
    //! Should be pure virtual, see child classes
    inline bool isXmax()
    {
        return isBoundary( 0, 1 );
    }
    //! Should be pure virtual, see child classes
    inline bool isYmin()
    {
        return isBoundary( 1, 0 );
    }
    //! Should be pure virtual, see child classes
    inline bool isYmax()
    {
        return isBoundary( 1, 1 );
    }
    //! Should be pure virtual, see child classes
    inline bool isZmin()
    {
        return isBoundary( 2, 0 );
    }
    //! Should be pure virtual, see child classes
    inline bool isZmax()
    {
        return isBoundary( 2, 1 );
    }
    //! Determine wether the patch is at the domain boundary
    inline bool isAnyBoundary()
    {
        bool flag = false;
        for( unsigned int i = 0 ; i < (unsigned int) nDim_fields_ ; i++ ) {
            flag = flag || isBoundary( i, 0 ) || isBoundary( i, 1 ) ;
        }
        return flag;
    }
    //! Define old xmax patch for moiving window,(non periodic eature)
    inline bool wasXmax( Params &params )
    {
        return Pcoordinates[0]+1 ==  params.number_of_patches[0]-1;
    }
    
    //! Test neighbbor's patch Id to apply or not a boundary condition
    inline bool isBoundary( unsigned int axis, unsigned int min_max )
    {
        return neighbor_[axis][min_max] == MPI_PROC_NULL;
    }
    //! Test neighbbor's patch Id to apply or not a boundary condition (xmin,xmax,ymin,ymax,...)
    inline bool isBoundary( unsigned int iboundary )
    {
        unsigned int axis = iboundary / 2;
        unsigned int min_max = iboundary % 2;
        return isBoundary( axis, min_max );
    }
    
    //! Compute MPI rank of neigbors patch regarding neigbors patch Ids
    void updateMPIenv( SmileiMPI *smpi );
    
    // Test who is MPI neighbor of current patch
    inline bool is_a_MPI_neighbor( int iDim, int iNeighbor )
    {
        return( ( neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL ) && ( MPI_neighbor_[iDim][iNeighbor]!=MPI_me_ ) );
    }
    
    inline bool has_an_MPI_neighbor()
    {
        for( unsigned int iDim=0 ; iDim<MPI_neighbor_.size() ; iDim++ ) {
            if( ( MPI_neighbor_[iDim][0] != MPI_me_ ) && ( MPI_neighbor_[iDim][0]!= MPI_PROC_NULL ) ) {
                return true;
            }
            if( ( MPI_neighbor_[iDim][1] != MPI_me_ ) && ( MPI_neighbor_[iDim][1]!= MPI_PROC_NULL ) ) {
                return true;
            }
        }
        return false;
    }
    
    inline bool has_an_MPI_neighbor( int iDim )
    {
        {
            if( ( MPI_neighbor_[iDim][0] != MPI_me_ ) && ( MPI_neighbor_[iDim][0]!= MPI_PROC_NULL ) ) {
                return true;
            }
            if( ( MPI_neighbor_[iDim][1] != MPI_me_ ) && ( MPI_neighbor_[iDim][1]!= MPI_PROC_NULL ) ) {
                return true;
            }
        }
        return false;
    }
    
    inline bool has_an_local_neighbor( int iDim )
    {
        {
            if( ( MPI_neighbor_[iDim][0] == MPI_me_ ) && ( MPI_neighbor_[iDim][0]!= MPI_PROC_NULL ) ) {
                return true;
            }
            if( ( MPI_neighbor_[iDim][1] == MPI_me_ ) && ( MPI_neighbor_[iDim][1]!= MPI_PROC_NULL ) ) {
                return true;
            }
        }
        return false;
    }
    
    
    //! Return real (excluding oversize) min coordinates (ex : rank 0 returns 0.) for direction i
    //! @see min_local_
    inline double getDomainLocalMin( int i ) const
    {
        return min_local_[i];
    }
    //! Return real (excluding oversize) max coordinates for direction i
    //! @see min_local_
    inline double getDomainLocalMax( int i ) const
    {
        return max_local_[i];
    }
    //! Return global starting (including oversize, ex : rank 0 returns -oversize) index for direction i
    //! \param i direction
    //! @see cell_starting_global_index
    inline int    getCellStartingGlobalIndex( int i ) const
    {
        return cell_starting_global_index[i];
    }
    //! Set global starting index for direction i
    //! @see cell_starting_global_index
    inline int    &getCellStartingGlobalIndex( int i )
    {
        return cell_starting_global_index[i];
    }

    //! Return global starting ( NOT including oversize, ex : rank 0 returns 0) index for direction i
    //! \param i direction
    //! @see cell_starting_global_index_noGC
    inline int    getCellStartingGlobalIndex_noGC( int i ) const
    {
        return cell_starting_global_index_noGC[i];
    }
    
    //! Set global starting index for direction i
    //! @see cell_starting_global_index_noGC
    inline int    &getCellStartingGlobalIndex_noGC( int i )
    {
        return cell_starting_global_index_noGC[i];
    }

    //! Set real min coordinate for direction i
    //! @see min_local_
    inline double &getDomainLocalMin( int i )
    {
        return min_local_[i];
    }
    //! Set real max coordinate for direction i
    //! @see max_local_
    inline double &getDomainLocalMax( int i )
    {
        return max_local_[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local_
    inline std::vector<double> getDomainLocalMin() const
    {
        return min_local_;
    }
    
    //! Return the volume (or surface or length depending on simulation dimension)
    //! of one cell at the position of a given particle
    virtual double getPrimalCellVolume( Particles *p, unsigned int ipart, Params &params ) = 0;
    
    //! Given several arrays (x,y,z for instance), return indices of points in patch
    virtual std::vector<unsigned int> indicesInDomain( double **position, unsigned int n_particles ) = 0;
    
    //! Set geometry data in case of moving window restart
    //! \param x_moved difference on coordinates regarding t0 geometry
    //! \param idx_moved number of displacement of the window
    //inline void updateMvWinLimits(double x_moved, int idx_moved) {
    //    min_local_[0] += x_moved;
    //    max_local_[0] += x_moved;
    //    //cell_starting_global_index[0] = (idx_moved-oversize[0]);
    //    cell_starting_global_index[0] += (idx_moved);
    //}
    
    //! Update Poyting quantities depending on location of the patch
    virtual void computePoynting();
    
    //! MPI rank of current patch
    int MPI_me_;
    
    //! The debye length, computed for collisions
    std::vector<double> debye_length_squared;
    
    //! The patch geometrical center
    std::vector<double> center_;
    //! The patch geometrical maximal radius (from its center)
    double radius;
    
    std::vector<MPI_Request> requests_;
    
    bool is_small = true;

    void copySpeciesBinsInLocalDensities(int ispec, int clrw, Params &params, bool diag_flag);
    void copySpeciesBinsInLocalSusceptibility(int ispec, int clrw, Params &params, bool diag_flag);
        
protected:
    // Complementary members for the description of the geometry
    // ---------------------------------------------------------
    
    //! Store number of space dimensions for the fields
    int nDim_fields_;
    
    //! Number of MPI process per direction in the cartesian topology (2)
    int nbNeighbors_;
    
    //! Hilbert index of neighbors patch
    std::vector< std::vector<int> > neighbor_, tmp_neighbor_;
    
    
    //! MPI rank of neighbors patch
    std::vector< std::vector<int> > MPI_neighbor_, tmp_MPI_neighbor_;
    
    //! "Real" min limit of local sub-subdomain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local_;
    //! "Real" max limit of local sub-subdomain (ghost data not concerned)
    std::vector<double> max_local_;
    //! cell_starting_global_index : index of 1st cell of local sub-subdomain in the global domain.
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    std::vector<int> cell_starting_global_index_noGC;
    double cell_volume;
    
    //! Buffers for exchange
    std::vector<int> buffer_vecto;
    std::vector<double> buffer_scalars_particles;
    std::vector<double> buffer_scalars_fields;
        
};


//! Return a unique id to identify all MPI communications
//!  - 2 MPI process can have several communications in the same direction for the same operation
//!  - the communication is identientified using :
//!      - hilbert index of the sender + idir + ineighbor
inline int buildtag( int hindex, int send, int recv )
{
    std::stringstream stag( "" );
    stag << hindex << send  << recv;
    long long int tag( 0 );
    stag >> tag;
    return ( int )( tag );
}

inline int buildtag( int hindex, int send, int recv, int tagp )
{
    std::stringstream stag( "" );
    //stag << hindex << send  << recv << tagp;
    stag << hindex << send*2+ recv << tagp;
    long long int tag( 0 );
    stag >> tag;
    return ( int )( tag );
}


#endif
