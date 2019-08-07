#ifndef SMILEIMPI_TEST_H
#define SMILEIMPI_TEST_H

#include "SmileiMPI.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Child class SmileiMPI_test for test mode only
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI_test : public SmileiMPI
{
    friend class Checkpoint;
    friend class PatchesFactory;
    friend class Patch;
    friend class VectorPatch;
    
public:

    //! Create empty MPI environment
    SmileiMPI_test( int *argc, char ***argv );
    //! Destructor for SmileiMPI_test
    ~SmileiMPI_test();
    
    //! Fake initialize  MPI (per process) environment
    void init( Params &params, DomainDecomposition *domain_decomposition ) override;
    
    // Initialize the patch_count vector. Patches are distributed in order to balance the load between MPI processes.
    void init_patch_count( Params &params, DomainDecomposition *domain_decomposition ) override;
    // Recompute the patch_count vector. Browse patches and redistribute them in order to balance the load between MPI processes.
    
};

#endif

