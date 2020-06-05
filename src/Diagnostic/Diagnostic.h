#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "H5.h"
#include "Patch.h"
#include "Timers.h"

class Params;
class OpenPMDparams;
class SmileiMPI;
class VectorPatch;
class TimeSelection;

class Diagnostic
{

public :

    Diagnostic( ) : openPMD_( NULL ) {};
    Diagnostic( OpenPMDparams *o, std::string diag_type, int idiag ) : openPMD_( o ) {
        PyTools::extract( "name", diag_name_, diag_type, idiag );
    };
    virtual ~Diagnostic() {};
    
    //! Opens the file. Only by MPI master for global diags. Only by patch master for local diags.
    virtual void openFile( Params &params, SmileiMPI *smpi, bool newfile ) = 0;
    //! Closes the file. Only by MPI master for global diags. Only by patch master for local diags.
    virtual void closeFile() = 0;
    
    //! Misc init.
    virtual void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) {};
    
    //! Prepares the diag and check whether it is time to run. Only by MPI master for global diags. Only by patch master for local diags.
    virtual bool prepare( int timestep ) = 0;
    
    //! Runs the diag for a given patch for global diags.
    virtual void run( Patch *patch, int timestep, SimWindow *simWindow ) {};
    
    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI *smpi, VectorPatch &vecPatches, int timestep, SimWindow *simWindow, Timers &timers ) {};
    
    //! Writes out a global diag diag.
    virtual void write( int timestep, SmileiMPI *smpi ) {};
    
    //! Tells whether this diagnostic requires the pre-calculation of the particle J & Rho
    virtual bool needsRhoJs( int timestep )
    {
        return false;
    };
    
    //! Time selection for writing the diagnostic
    TimeSelection *timeSelection;
    
    //! Time selection for flushing the file
    TimeSelection *flush_timeSelection;
    
    //! Get memory footprint of current diagnostic
    virtual int getMemFootPrint() = 0;
    
    //! Get disk footprint of current diagnostic
    virtual uint64_t getDiskFootPrint( int istart, int istop, Patch *patch )
    {
        return 0.;
    };
    
    //! this is the file name
    std::string filename;
    
    bool theTimeIsNow;
    
    std::string name() {
        return diag_name_;
    };
    
protected :

    //! Id of the file for one diagnostic
    hid_t fileId_;
    
    //! Pointer to all parameters needed for openPMD compatibility
    OpenPMDparams *openPMD_;
    
    //! Label of the diagnostic (for post-processing)
    std::string diag_name_;
};

#endif

