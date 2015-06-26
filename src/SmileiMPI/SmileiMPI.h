#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <string>
#include <vector>

#include <mpi.h>

#include "PicParams.h"
#include "Tools.h"
#include "Particles.h"

class PicParams;
class DiagParams;
class Diagnostic;
class DiagnosticScalar;
class Species;
class VectorPatch;

class ElectroMagn;
class Field;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPI
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI {
public:
    friend class Patch;

    //! Create intial MPI environment
    SmileiMPI( int* argc, char*** argv );
    //! Create MPI environment for the data geometry from 
    //! \param smpi the initil MPI environment
    SmileiMPI(SmileiMPI *smpi);
    //! Default creator for SmileiMPI
    SmileiMPI() {};
    //! Destructor for SmileiMPI
    virtual ~SmileiMPI();

    //! Initialize geometry members dimension from
    //! \param params Parameters
    //! @see oversize
    //! @see cell_starting_global_index
    //! @see min_local
    //! @see max_local
    //! @see n_space_global
    void init( PicParams& params );
    // Initialize the patch_count vector. Patches are distributed in order to balance the load between MPI processes.
    void init_patch_count( PicParams& params );
    // Recompute the patch_count vector. Browse patches and redistribute them in order to balance the load between MPI processes.
    void recompute_patch_count( PicParams& params, VectorPatch& vecpatches, double time_dual );
    //! Broadcast to all process
    //! \param idata read data
    void bcast( InputData& idata );

    //! Create MPI communicator
    virtual void createTopology( PicParams& params ) {};
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    //! See child classes
    virtual void exchangeParticles(Species* species, int ispec, PicParams& params, int tnum, int iDim) {};

    MPI_Datatype createMPIparticles( Particles* particles, int nbrOfProp );


    //! Create MPI_Datatype to exchange/sum fields on ghost data
    //! See child classes
    virtual void createType( PicParams& params ) {};

    //! Exchange all electric fields on borders
    void exchangeE( ElectroMagn* EMfields );
    //! Exchange all magnectic fields on borders
    void exchangeB( ElectroMagn* EMfields );
    //! Exchange all centered magnectic fields on borders
    void exchangeBm( ElectroMagn* EMfields );

    void exchangeAvg( ElectroMagn* EMfields );

    //! Exchange clrw columns of electric fields towards the west
    void exchangeE( ElectroMagn* EMfields, unsigned int clrw );
    //! Exchange clrw columns of magnetic field towards the west
    void exchangeB( ElectroMagn* EMfields, unsigned int clrw);
    //! Exchange clrw columns of centered magnetic field towards the west
    void exchangeBm( ElectroMagn* EMfields, unsigned int clrw );

    //! Sum rho and densities on 2 x oversize[]
    void sumRho( ElectroMagn* EMfields );
    //! Sum rho and all J on the shared domain between processors
    //! 2 x oversize + 1 ( + 1 if direction is dual )
    void sumRhoJ( ElectroMagn* EMfields );
    //! Sum rho_s and all J_s on the shared domain between processors
    void sumRhoJs( ElectroMagn* EMfields, int ispec, bool currents );

    //! Basic method to exchange a field, defined in child class
    virtual void exchangeField ( Field* field ) {};
    //! Basic method to exchange a field towards the west, defined in child class
    virtual void exchangeField_movewin ( Field* field, unsigned int clrw ) {};
    //! Basic method to sum a field, defined in child class
    virtual void sumField      ( Field* field ) {};


    // --------------------------------------------------
    // ------------- PATCH EXCHANGE METHODS -------------
    // --------------------------------------------------
    void send(Patch* patch, int to);
    void recv(Patch* patch, int from);

    void send(int value, int to);
    void recv(int value, int from);

    void send(Species* species, int to);
    void recv(Species* species, int from);
    void send(Particles* particles, int to);
    void recv(Particles* partictles, int from);
    void send(std::vector<int> vec, int to);
    void recv(std::vector<int> vec, int from);

    void send(ElectroMagn* fields, int to);
    void recv(ElectroMagn* fields, int from);
    void send(Field* field, int to);
    void recv(Field* field, int from);
    // --------------------------------------------------
    // ------ END OF PATCH EXCHANGE METHODS -------------
    // --------------------------------------------------

    void computeGlobalDiags(Diagnostic* diags, int timestep);
    void computeGlobalDiags(DiagnosticScalar& scalars, int timestep);

    //! Method to identify the rank 0 MPI process
    inline bool isMaster() {
        return (smilei_rk==0);
    }
    //! Method to synchronize MPI process in the current MPI communicator
    inline void barrier() {
        MPI_Barrier( SMILEI_COMM_WORLD );
    }
    //! Return MPI_Comm_rank
    inline int getRank() {
        return smilei_rk;
    }
    //! Return MPI_Comm_size
    inline int getSize() {
        return smilei_sz;
    }
    //! Return global starting (including oversize, ex : rank 0 returns -oversize) index for direction i
    //! \param i direction
    //! @see cell_starting_global_index
    inline int    getCellStartingGlobalIndex(int i) const {
        return cell_starting_global_index[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline double getDomainLocalMin(int i) const {
        return min_local[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline std::vector<double> getDomainLocalMin() const {
        return min_local;
    }
    //! Return real (excluding oversize) max coordinates for direction i
    //! @see max_local
    inline double getDomainLocalMax(int i) const {
        return max_local[i];
    }

    //! Set geometry data in case of moving window restart
    //! \param x_moved difference on coordinates regarding t0 geometry
    //! \param idx_moved number of displacement of the window
    inline void updateMvWinLimits(double x_moved, int idx_moved) {
	min_local[0] += x_moved;
	max_local[0] += x_moved;
	cell_starting_global_index[0] = (idx_moved-oversize[0]);
    }

    //! Set global starting index for direction i
    //! @see cell_starting_global_index
    inline int&    getCellStartingGlobalIndex(int i)  {
        return cell_starting_global_index[i];
    }
    //! Set real min coordinate for direction i
    //! @see min_local
    inline double& getDomainLocalMin(int i)  {
        return min_local[i];
    }
    //! Set real max coordinate for direction i
    //! @see max_local
    inline double& getDomainLocalMax(int i)  {
        return max_local[i];
    }

    //! Temporary storage of particles Id to exchange, merge from per thread storage
    //! A single communication per direction managed by thread master in OpenMP
    //! @see Species::indexes_of_particles_to_exchange_per_thrd
    std::vector<int>                 indexes_of_particles_to_exchange;

    //! Should be pure virtual, see child classes 
    virtual bool isEastern(){WARNING("Problem");return false;}
    //! Should be pure virtual, see child classes 
    virtual bool isWestern(){WARNING("Problem");return false;}
    //! Should be pure virtual, see child classes 
    virtual bool isSouthern(){WARNING("Problem");return false;}
    //! Should be pure virtual, see child classes 
    virtual bool isNorthern(){WARNING("Problem");return false;}

    //! Real (exclunding oversize) global number of cells (res_space x sim_length)
    std::vector<int> n_space_global;
    //! Number of MPI process in the current communicator
    int smilei_sz;
    //! MPI process Id in the current communicator
    int smilei_rk;

    //! For patch decomposition
    Particles interParticles;
    std::vector<int>  patch_count, target_patch_count;  //Number of patches owned by each mpi process.
    int hrank(int h); // Returns the rank of the MPI process currently owning patch h.

protected:
    //! Global MPI Communicator
    MPI_Comm SMILEI_COMM_WORLD;

    int* periods_;

    //! Sort particles to exchange per side (2), contains indexes
    //! Reinitialized per direction
    std::vector<int> buff_index_send[2];
    //! buff_index_recv_sz : number of particles to recv per side (2)
    //! Reinitialized per direction
    int buff_index_recv_sz[2];

    //! Size of ghost data (= 2 x oversize + 1 + 1 if dual direction), depend on :
    //!    - projection/interpolation order
    //!    - rate of particles exchange (to implement)
    std::vector<unsigned int> oversize;
    //! cell_starting_global_index : index of 1st cell of local subdomain in the global domain
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    //! "Real" min limit of local domain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local;
    //! "Real" max limit of local domain (ghost data not concerned)
    std::vector<double> max_local;

private:
    // Broadcast a string in current communicator
    void bcast( std::string& val );

};

#endif

