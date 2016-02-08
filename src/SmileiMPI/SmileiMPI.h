#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <string>
#include <vector>

#include <mpi.h>

#include "Params.h"
#include "Tools.h"
#include "Particles.h"
#include "Field.h"

class Params;
class Species;
class VectorPatch;

class ElectroMagn;
class Field;
class Diagnostic;
class DiagnosticScalar;
class DiagnosticPhaseSpace;
class DiagnosticParticles;

#define SMILEI_COMM_DUMP_TIME 1312

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPI
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI {
    friend class SmileiIO;
    friend class Checkpoint;
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
    void init( Params& params );
    // Initialize the patch_count vector. Patches are distributed in order to balance the load between MPI processes.
    void init_patch_count( Params& params );
    // Recompute the patch_count vector. Browse patches and redistribute them in order to balance the load between MPI processes.
    void recompute_patch_count( Params& params, VectorPatch& vecpatches, double time_dual );

    MPI_Datatype createMPIparticles( Particles* particles );



    // --------------------------------------------------
    // ------------- PATCH EXCHANGE METHODS -------------
    // --------------------------------------------------
    void send(Patch* patch, int to  , int hindex);
    void isend(Patch* patch, int to  , int hindex);
    void recv(Patch* patch, int from, int hindex);
    void new_recv(Patch* patch, int from, int hindex, Params& params);

    void send(Species* species, int to  , int hindex);
    void recv(Species* species, int from, int hindex);
    void send(Particles* particles, int to   , int hindex);
    void isend(Particles* particles, int to   , int hindex, MPI_Datatype datatype);
    void recv(Particles* partictles, int from, int hindex);
    void new_recv(Particles* partictles, int from, int hindex, MPI_Datatype datatype);
    void send(std::vector<int> vec, int to  , int hindex);
    void isend(std::vector<int>* vec, int to  , int hindex);
    void recv(std::vector<int> *vec, int from, int hindex);

    void send(ElectroMagn* fields, int to  , int hindex);
    void isend(ElectroMagn* fields, int to  , int hindex);
    void recv(ElectroMagn* fields, int from, int hindex);
    void send(Field* field, int to  , int hindex);
    void isend(Field* field, int to  , int hindex);
    void recv(Field* field, int from, int hindex);
    void send( Diagnostic* diags, int to  , int hindex );
    void isend( Diagnostic* diags, int to  , int hindex );
    void recv( Diagnostic* diags, int from, int hindex );
    // --------------------------------------------------
    // ------ END OF PATCH EXCHANGE METHODS -------------
    // --------------------------------------------------

    void computeGlobalDiags(Diagnostic* diags, int timestep);
    void computeGlobalDiags(DiagnosticScalar& scalars, int timestep);
    void computeGlobalDiags(DiagnosticPhaseSpace& phases, int timestep);
    void computeGlobalDiags(DiagnosticParticles* diagParticles, int timestep);

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

    // Global buffers for vectorization of Species::dynamics
    // -----------------------------------------------------
    
    //! value of the Efield 
    std::vector<std::vector<LocalFields>> dynamics_Epart;
    //! value of the Bfield
    std::vector<std::vector<LocalFields>> dynamics_Bpart;
    //! gamma factor
    std::vector<std::vector<double>> dynamics_gf;
    //! iold_pos
    std::vector<std::vector<int>> dynamics_iold;
    //! delta_old_pos
    std::vector<std::vector<double>> dynamics_deltaold;

    inline void dynamics_resize(int ithread, int ndim_part, int npart ){
        dynamics_Epart[ithread].resize(npart);
        dynamics_Bpart[ithread].resize(npart);
        dynamics_gf[ithread].resize(npart);
        dynamics_iold[ithread].resize(ndim_part*npart);
        dynamics_deltaold[ithread].resize(ndim_part*npart);
    }


    //! Number of MPI process in the current communicator
    int smilei_sz;
    //! MPI process Id in the current communicator
    int smilei_rk;

    //! For patch decomposition
    std::vector<int>  patch_count, target_patch_count;  //Number of patches owned by each mpi process.
    int hrank(int h); // Returns the rank of the MPI process currently owning patch h.

    inline int globalNbrParticles(Species* species, int locNbrParticles) {
	int nParticles(0);
	MPI_Reduce( &locNbrParticles, &nParticles, 1, MPI_INT, MPI_SUM, 0, SMILEI_COMM_WORLD );
	return nParticles;
    }


    // Broadcast a string in current communicator
    void bcast( std::string& val );
    // Broadcast an int in current communicator
    void bcast( int& val );

protected:
    //! Global MPI Communicator
    MPI_Comm SMILEI_COMM_WORLD;

    int* periods_;

};

#endif

