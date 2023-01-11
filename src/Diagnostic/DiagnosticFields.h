#ifndef DIAGNOSTICFIELDS_H
#define DIAGNOSTICFIELDS_H

#include "Diagnostic.h"

class DiagnosticFields  : public Diagnostic
{

public :

    DiagnosticFields( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticFields() override;
    
    virtual void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    virtual void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    virtual bool prepare( int itime ) override;
    
    virtual void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) = 0;
    
    virtual void run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers ) override;
    
    virtual H5Write writeField( H5Write*, std::string ) = 0;
    
    virtual bool needsRhoJs( int itime ) override;
    
    void findSubgridIntersection( unsigned int subgrid_start,
                                  unsigned int subgrid_stop,
                                  unsigned int subgrid_step,
                                  unsigned int zone_begin,
                                  unsigned int zone_end,
                                  unsigned int &istart_in_zone,
                                  unsigned int &istart_in_file,
                                  unsigned int &nsteps );
    void findSubgridIntersection1( hsize_t idim,
                                   hsize_t &zone_begin,
                                   hsize_t &zone_npoints,
                                   hsize_t &start_in_zone );
                                  
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        return 0;
    };
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
protected :

    //! Index of this diag
    unsigned int diag_n;
    
    //! Indexes of the fields to be dumped
    std::vector<unsigned int> fields_indexes;
    //! Names of the fields to be dumped
    std::vector<std::string> fields_names;
    
    //! Number of timesteps for time averaging
    int time_average;
    
    //! Inverse of the time average
    double time_average_inv;
    
    //! Subgrid requested
    std::vector<unsigned int> subgrid_start_, subgrid_stop_, subgrid_step_;
    
    //! Number of cells to skip in each direction
    std::vector<unsigned int> patch_offset_in_grid;
    //! Number of cells in each direction
    std::vector<unsigned int> patch_size_;
    //! Buffer for the output of a field
    std::vector<double> data;
    
    //! 1st patch index of vecPatches
    unsigned int refHindex;
    
    H5Write *data_group_, *iteration_group_;
    H5Space *filespace, *memspace;
    
    //! Total number of patches
    int tot_number_of_patches;
    
    //! Copy patch field to current "data" buffer
    virtual void getField( Patch *patch, unsigned int ) = 0;
    
    //! Variable to store the status of a dataset (whether it exists or not)
    bool status;
    
    //! Total size of the data in file (for estimating disk usage)
    unsigned int total_dataset_size;
    
    //! True if this diagnostic requires the pre-calculation of the particle J & Rho
    bool hasRhoJs;
    
    //! Save the field type (needed for OpenPMD units dimensionality)
    std::vector<unsigned int> field_type;
};

#endif

