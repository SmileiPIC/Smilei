#ifndef SMILEIMPI_TEST_H
#define SMILEIMPI_TEST_H

#include "SmileiMPI.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Child class SmileiMPI_test for test mode only
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI_test : public SmileiMPI {
    friend class Checkpoint;
    friend class PatchesFactory;
    friend class Patch;
    friend class VectorPatch;

public:
    
    //! Create empty MPI environment
    SmileiMPI_test( int nMPI, int nOMP );
    //! Destructor for SmileiMPI_test
    ~SmileiMPI_test();
    
    // Broadcast nothing
    void bcast( std::string& val ) override {};
    void bcast( int& val ) override {};
    
    //! Fake initialize  MPI (per process) environment
    void init( Params& params ) override;
    
    // Initialize the patch_count vector. Patches are distributed in order to balance the load between MPI processes.
    void init_patch_count( Params& params ) override;
    // Recompute the patch_count vector. Browse patches and redistribute them in order to balance the load between MPI processes.
    void recompute_patch_count( Params& params, VectorPatch& vecpatches, double time_dual ) override {};
    
    // Creates nothing
    MPI_Datatype createMPIparticles( Particles* particles ) override { return NULL; };
    
    // Fake communication functions
    void isend(Patch* patch, int to  , int hindex, Params& params) override {};
    void waitall(Patch* patch) override {};
    void recv (Patch* patch, int from, int hindex, Params& params) override {};
    void isend(Particles* particles, int to   , int hindex, MPI_Datatype datatype, MPI_Request& request) override {};
    void recv (Particles* partictles, int from, int hindex, MPI_Datatype datatype) override {};
    void isend(std::vector<int>* vec, int to  , int hindex, MPI_Request& request) override {};
    void recv (std::vector<int> *vec, int from, int hindex) override {};
    void isend(std::vector<double>* vec, int to  , int hindex, MPI_Request& request) override {};
    void recv (std::vector<double> *vec, int from, int hindex) override {};
    void isend(ElectroMagn* fields, int to  , int maxtag, std::vector<MPI_Request>& requests, int mpi_tag) override {};
    void recv (ElectroMagn* fields, int from, int hindex) override {};
    void isend(Field* field, int to  , int hindex, MPI_Request& request) override {};
    void recv (Field* field, int from, int hindex) override {};
    void isend( ProbeParticles* probe, int to  , int hindex, unsigned int ) override {};
    void recv ( ProbeParticles* probe, int from, int hindex, unsigned int ) override {};
    
    // Resize nothing
    inline void dynamics_resize(int ithread, int ndim_part, int npart ) override {};
    
    // Compute nothing
    inline int globalNbrParticles(Species* species, int locNbrParticles) override { return 0; };
    
    //! Method to synchronize nothing
    inline void barrier() override {};
};

#endif

