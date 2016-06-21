#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class Diagnostic {

public :

    Diagnostic() {};
    virtual ~Diagnostic() {};
    
    //! Opens the file. Only by MPI master for global diags. Only by patch master for local diags.
    virtual void openFile( Params& params, SmileiMPI* smpi, bool newfile ) = 0;
    //! Closes the file. Only by MPI master for global diags. Only by patch master for local diags.
    virtual void closeFile() = 0;
    
    //! Splits the file for local diags only. Only by patch master.
    virtual void setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches ) {};
    
    //! Misc init.
    virtual void init(SmileiMPI* smpi, VectorPatch& vecPatches) {};
    
    //! Prepares the diag and check whether it is time to run. Only by MPI master for global diags. Only by patch master for local diags.
    virtual bool prepare( int timestep ) = 0;
    
    //! Runs the diag for a given patch for global diags.
    virtual void run( Patch* patch, int timestep ) {};
    
    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep ) {};
    
    //! Writes out the diag.
    virtual bool write(int timestep) {};
    
    //! Does some more work after writing
    virtual void finish(int, VectorPatch& ) {};
    
    hid_t getFileId() {
        return fileId_;
    }
    void setFileId( hid_t fileId ) {
        fileId_ = fileId;
    }
    
    //! Time selection
    TimeSelection * timeSelection;
    
    //! this is the file name
    std::string filename;
    std::string type_;
protected :
    int diagId_;
    
    //! Id of the file for one diagnostic
    hid_t fileId_;
};

#endif

