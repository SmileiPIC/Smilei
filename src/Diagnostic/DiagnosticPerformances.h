#ifndef DIAGNOSTICPERFORMANCES_H
#define DIAGNOSTICPERFORMANCES_H

#include "Diagnostic.h"
#include "VectorPatch.h"

class DiagnosticPerformances : public Diagnostic {
    friend class SmileiMPI;

public :
    
    //! Default constructor
    DiagnosticPerformances( SmileiMPI* smpi );
    //! Default destructor
    ~DiagnosticPerformances() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    bool prepare( int timestep ) override;
    
    void run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime, SimWindow* simWindow, Timers & timers ) override;
    
    void writeQuantity( double       quantity, const char* name, hid_t gid, hid_t create_plist );
    void writeQuantity( unsigned int quantity, const char* name, hid_t gid, hid_t create_plist );
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        return 0;
    };
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint(int istart, int istop, Patch* patch) override;
    
private :
    
    hid_t iteration_group_id, filespace, memspace;
    
    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;
    
    //! Variable to store the status of a dataset (whether it exists or not)
    htri_t status;
    
    hsize_t mpi_size;
    
};

#endif

