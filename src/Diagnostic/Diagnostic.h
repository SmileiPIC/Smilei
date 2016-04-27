#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class Diagnostic {

public :

    Diagnostic() {};
    ~Diagnostic() {};
    
    //! Opens the file. Only by MPI master for global diags. Only by patch master for local diags.
    virtual void openFile( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches, bool newfile ) = 0;
    //! Closes the file. Only by MPI master for global diags. Only by patch master for local diags.
    virtual void closeFile() = 0;
    
    //! Splits the file for local diags only. Only by patch master.
    virtual void setFileSplitting( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches ) {
        ERROR("Should not happen");
    };
    
    //! Prepares the diag and check whether it is time to run. Only by MPI master for global diags. Only by patch master for local diags.
    virtual bool prepare( Patch* patch, int timestep ) = 0;
    
    //! Runs the diag. By all patches.
    virtual void run( Patch* patch, int timestep ) = 0;
    
    //! Writes out the diag. By all patches.
    virtual void write(int timestep) = 0;
    
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
    int probeId_;
    
    //! Id of the file for one diagnostic
    hid_t fileId_;
};

#endif

