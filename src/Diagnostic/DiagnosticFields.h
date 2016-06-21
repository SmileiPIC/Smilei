#ifndef DIAGNOSTICFIELDS_H
#define DIAGNOSTICFIELDS_H

#include "Diagnostic.h"


class DiagnosticFields  : public Diagnostic {

public :
    
    DiagnosticFields( Params &params, SmileiMPI* smpi, Patch* patch, int );
    ~DiagnosticFields() override;
    
    virtual void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    virtual void closeFile() override;
    
    virtual void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    virtual bool prepare( int timestep ) override;
    
    virtual void setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches ) = 0;
    
    virtual void run( Patch* patch, int timestep ) {};
    virtual void run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep ) override;
    
    virtual bool write(int timestep) {};
    
    virtual void finish(int, VectorPatch& ) {};
    
    virtual void writeField(hid_t, int) = 0;
    
protected :
    //! Indexes of the fields to be dumped
    std::vector<int> fields_indexes;
    //! Names of the fields to be dumped
    std::vector<std::string> fields_names;
    
    //! Number of timesteps for time averaging
    int time_average;
    
    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;
    
    //! Number of cells to skip in each direction
    std::vector<unsigned int> patch_offset_in_grid;
    //! Number of cells in each direction
    std::vector<unsigned int> patch_size;
    //! Number of cells in a patch
    unsigned int total_patch_size;
    //! Buffer for the output of a field
    std::vector<double> data;
    
    //! 1st patch index of vecPatches
    unsigned int refHindex;
    
    hid_t timestep_group_id, filespace, memspace;
    
    //! Total number of patches
    int tot_number_of_patches;
    
    //! Copy patch field to current "data" buffer
    virtual void getField( Patch* patch, int ) = 0;
    
    //! Temporary dataset that is used for folding the 2D hilbert curve
    hid_t tmp_dset_id;
};

#endif

