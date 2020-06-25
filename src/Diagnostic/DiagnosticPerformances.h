#ifndef DIAGNOSTICPERFORMANCES_H
#define DIAGNOSTICPERFORMANCES_H

#include "Diagnostic.h"
#include "VectorPatch.h"

class DiagnosticPerformances : public Diagnostic
{
    friend class SmileiMPI;
    
public :

    //! Default constructor
    DiagnosticPerformances( Params &params, SmileiMPI *smpi );
    //! Default destructor
    ~DiagnosticPerformances() override;
    
    void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    bool prepare( int itime ) override;
    
    void run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers ) override;
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        return 0;
    };
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
    //! Set the hdf5 spaces for 2D arrays with one column selected per proc
    void setHDF5spaces( hid_t &filespace, hid_t &memspace, unsigned int height, unsigned int width, unsigned int column );
    
private :

    //! Number of MPI ranks
    hsize_t mpi_size_;
    
    //! MPI rank
    hsize_t mpi_rank_;
    
    //! HDF5 link to the group corresponding to one iteration
    bool has_group;
    std::string group_name;
    
    //! HDF5 shapes of datasets
    H5Space filespace_double, filespace_uint;
    H5Space memspace_double, memspace_uint;
    
    //! Total number of patches
    unsigned int tot_number_of_patches;
    
    //! Dimension of the patches
    unsigned int ndim;
    
    //! Whether the adaptive vectorization is active
    bool has_adaptive_vectorization;
    
    //! Whether to output patch information
    bool patch_information;
    
    //! Number of cells per patch
    unsigned int ncells_per_patch;
    
    double timestep, cell_load, frozen_particle_load;
};

#endif
