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

class DomainDecomposition;
class Collisions;
class Diagnostic;
class SimWindow;

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
    Patch( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved );
    //! Cloning Constructor for Patch
    Patch( Patch *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles );
    
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
   
    //Copy positions of particles from source species to species which are initialized on top of another one.
    void copyPositions( std::vector<Species *> vecSpecies_to_update);

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
    //! Optional binary collisions operators
    std::vector<Collisions *> vecCollisions;
    
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
    
    // Detailed timers
    // -----------------------
    
#ifdef  __DETAILED_TIMERS
    //! Timers for the patch
    std::vector<double> patch_timers;
#endif
    
    // Random number generator.
    Random * rand_;
    
    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------
    
    //! Clean the MPI buffers for communications
    void cleanMPIBuffers( int ispec, Params &params );
    //! manage Idx of particles per direction,
    void initExchParticles( SmileiMPI *smpi, int ispec, Params &params );
    //! init comm  nbr of particles
    void exchNbrOfParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! finalize comm / nbr of particles, init exch / particles
    void endNbrOfParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! extract particles from main data structure to buffers, init exch / particles
    void prepareParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! effective exchange of particles
    void exchParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! finalize exch / particles
    void finalizeExchParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! Treat diagonalParticles
    void cornersParticles( SmileiMPI *smpi, int ispec, Params &params, int iDim, VectorPatch *vecPatch );
    //! inject particles received in main data structure and particles sorting
    void importAndSortParticles( SmileiMPI *smpi, int ispec, Params &params, VectorPatch *vecPatch );
    //! clean memory resizing particles structure
    void cleanParticlesOverhead( Params &params );
    //! delete Particles included in the index of particles to exchange. Assumes indexes are sorted.
    void cleanupSentParticles( int ispec, std::vector<int> *indexes_of_particles_to_exchange );
    
    //! init comm / sum densities
    virtual void initSumField( Field *field, int iDim, SmileiMPI *smpi ) = 0;
    //! finalize comm / sum densities
    virtual void finalizeSumField( Field *field, int iDim ) = 0;
    //! init comm / sum densities
    virtual void initSumFieldComplex( Field *field, int iDim, SmileiMPI *smpi ) = 0;
    //! finalize comm / sum densities
    virtual void finalizeSumFieldComplex( Field *field, int iDim ) = 0;
    
    //! init comm / exchange fields in direction iDim only
    virtual void initExchange( Field *field, int iDim, SmileiMPI *smpi ) = 0;
    //! init comm / exchange complex fields in direction iDim only
    virtual void initExchangeComplex( Field *field, int iDim, SmileiMPI *smpi ) = 0;
    //! finalize comm / exchange fields
    virtual void finalizeExchange( Field *field, int iDim ) = 0;
    //! finalize comm / exchange complex fields in direction iDim only
    virtual void finalizeExchangeComplex( Field *field, int iDim ) = 0;
    
    virtual void exchangeField_movewin ( Field* field, int clrw ) = 0;
    
    // Create MPI_Datatype to exchange fields
    virtual void createType( Params &params ) = 0;
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
        return locateOnBorders( 0, 0 );
    }
    //! Should be pure virtual, see child classes
    inline bool isXmax()
    {
        return locateOnBorders( 0, 1 );
    }
    //! Should be pure virtual, see child classes
    inline bool isYmin()
    {
        return locateOnBorders( 1, 0 );
    }
    //! Should be pure virtual, see child classes
    inline bool isYmax()
    {
        return locateOnBorders( 1, 1 );
    }
    //! Should be pure virtual, see child classes
    inline bool isZmin()
    {
        return locateOnBorders( 2, 0 );
    }
    //! Should be pure virtual, see child classes
    inline bool isZmax()
    {
        return locateOnBorders( 2, 1 );
    }
    //! Determine wether the patch is at the domain boundary
    inline bool isBoundary()
    {
        bool flag = false;
        for (int i = 0 ; i < nDim_fields_ ; i++) {
            flag = flag || locateOnBorders( i, 0 ) || locateOnBorders( i, 1 ) ;
        }
        return flag;
    }
    //! Define old xmax patch for moiving window,(non periodic eature)
    inline bool wasXmax( Params &params )
    {
        return Pcoordinates[0]+1 ==  params.number_of_patches[0]-1;
    }
    
    //! Test neighbbor's patch Id to apply or not a boundary condition
    inline bool locateOnBorders( int dir, int way )
    {
        if( neighbor_[dir][way] == MPI_PROC_NULL ) {
            return true;
        }
        return false;
    }
    
    //! Compute MPI rank of neigbors patch regarding neigbors patch Ids
    void updateMPIenv( SmileiMPI *smpi );
    void updateTagenv( SmileiMPI *smpi );
    
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
    //! @see min_local
    inline double getDomainLocalMin( int i ) const
    {
        return min_local[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 returns 0.) for direction i
    //! @see min_local
    inline double getDomainLocalMax( int i ) const
    {
        return max_local[i];
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
    //! Set real min coordinate for direction i
    //! @see min_local
    inline double &getDomainLocalMin( int i )
    {
        return min_local[i];
    }
    //! Set real max coordinate for direction i
    //! @see max_local
    inline double &getDomainLocalMax( int i )
    {
        return max_local[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline std::vector<double> getDomainLocalMin() const
    {
        return min_local;
    }
    
    //! Return the volume (or surface or length depending on simulation dimension)
    //! of one cell at the position of a given particle
    virtual double getPrimalCellVolume( Particles *p, unsigned int ipart, Params &params ) = 0;
    
    //! Set geometry data in case of moving window restart
    //! \param x_moved difference on coordinates regarding t0 geometry
    //! \param idx_moved number of displacement of the window
    //inline void updateMvWinLimits(double x_moved, int idx_moved) {
    //    min_local[0] += x_moved;
    //    max_local[0] += x_moved;
    //    //cell_starting_global_index[0] = (idx_moved-oversize[0]);
    //    cell_starting_global_index[0] += (idx_moved);
    //}
    
    //! MPI rank of current patch
    int MPI_me_;
    
    //! The debye length, computed for collisions
    std::vector<double> debye_length_squared;
    
    //! The patch geometrical center
    std::vector<double> center;
    //! The patch geometrical maximal radius (from its center)
    double radius;
    
    std::vector<MPI_Request> requests_;
    
    bool is_small = true;
    
        
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
    std::vector<double> min_local;
    //! "Real" max limit of local sub-subdomain (ghost data not concerned)
    std::vector<double> max_local;
    //! cell_starting_global_index : index of 1st cell of local sub-subdomain in the global domain.
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    
    std::vector<unsigned int> oversize;
    
    double cell_volume;
    
    //! Buffers for exchange
    std::vector<int> buffer_vecto;
    std::vector<double> buffer_scalars;
        
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
