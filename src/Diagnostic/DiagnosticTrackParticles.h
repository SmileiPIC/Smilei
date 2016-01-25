/*
-----------------------------------------------------------------------
WRITE PARTICLE DIAGNOSTICS
-----------------------------------------------------------------------
*/



#ifndef DiagnosticTrackParticles_H
#define DiagnosticTrackParticles_H

#include <cmath>

#include "Params.h"
#include "Species.h"
#include "Particles.h"
#include "SmileiMPI.h"
#include "H5.h"


// Class for the writable particles diagnostics
class DiagnosticTrackParticles {

public:

    DiagnosticTrackParticles(Params&, Patch* patch, Species*);
    void createFile(Params&, Patch* patch, Species*, Diagnostic*);
    
    ~DiagnosticTrackParticles(){};

    void createFile(int nParticles, Params& params);
    void setGlobalNbrParticles(int nParticles) { dims[1] = nParticles;}

    //! Runs the diag (writes to file) at each timestep
    void run( int );
    void open();
    void close();
    void setFile( hid_t, hid_t );
    //! HDF5 memory space (dimensions of the array in memory)
    //hid_t mem_space;
    hid_t fid_;
    hid_t file_access_;

private:
    
    //! HDF5 file transfer protocol
    hid_t transfer;
    //! HDF5 file space (dimensions of the array in file)
    hsize_t dims[2];
    
    //! Pointer to the species used
    Species* species;
    
    //! Number of spatial dimensions
    int nDim_particle;
    
    //! Adds one row in a HDF5 file, within a given dataspace
    template <class T> void append(hid_t, std::string, std::vector<T>, hid_t, int, hid_t, std::vector<hsize_t>&);

    int iter;
    
};




#endif
