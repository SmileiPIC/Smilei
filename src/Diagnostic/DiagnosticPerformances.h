#ifndef DIAGNOSTICPERFORMANCES_H
#define DIAGNOSTICPERFORMANCES_H

#include "Diagnostic.h"
#include "VectorPatch.h"

class DiagnosticPerformances : public Diagnostic {
    friend class SmileiMPI;

public :

    //! Default constructor
    DiagnosticPerformances( Params & params, SmileiMPI* smpi );
    //! Default destructor
    ~DiagnosticPerformances() override;

    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;

    void closeFile() override;

    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;

    bool prepare( int itime ) override;

    void run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime, SimWindow* simWindow, Timers & timers ) override;

    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        return 0;
    };

    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint(int istart, int istop, Patch* patch) override;

    //! Set the hdf5 spaces for 2D arrays with one column selected per proc
    void setHDF5spaces(hid_t &filespace, hid_t &memspace, unsigned int height, unsigned int width, unsigned int column);

    //! Set the hdf5 spaces for 1D arrays for patches
    void set_HDF5_patch_spaces(hid_t &filespace,
                               hid_t &memspace,
                               std::vector <unsigned int> &patches_per_mpi,
                               const unsigned int mpi_rank);

private :

    //! HDF5 link to the group corresponding to one iteration
    hid_t iteration_group_id;

    //! HDF5 shapes of datasets
    hid_t filespace_double, memspace_double;
    hid_t filespace_uint  , memspace_uint  ;
    //! HDF5 shapes of patch datasets
    hid_t filespace_patches, memspace_patches;

    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;

    //! Variable to store the status of a dataset (whether it exists or not)
    htri_t status;

    //! Number of MPI ranks
    hsize_t mpi_size;

    //! MPI rank
    unsigned int mpi_rank_;

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
