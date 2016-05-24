#ifndef DIAGNOSTICTRACK_H
#define DIAGNOSTICTRACK_H

#include "Diagnostic.h"

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


class DiagnosticTrack : public Diagnostic {

public :
    //! Default constructor
    DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, int diagId, int );
    //! Cloning constructor
    DiagnosticTrack(DiagnosticTrack* track, Patch* patch);
    //! Default destructor
    ~DiagnosticTrack() ;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    void setFileSplitting( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches ) override;
    
    void closeFile() override;
    
    bool prepare( int timestep ) override;
    
    void run( Patch* patch, int timestep ) override;
    
    bool write(int timestep) override;
    
    void setGlobalNbrParticles(int totNbrParts) {
        nbrParticles_ = totNbrParts;
    }
    
private :
    //! Pointer to the species used
    Species* species;
    int speciesId_;
    
    //! Size of the diag (number of particles)
    int nbrParticles_;
     
    //! HDF5 file transfer protocol
    hid_t transfer;
    //! HDF5 file space (dimensions of the array in file)
    hsize_t dims[2];
     
    //! Number of spatial dimensions
    int nDim_particle;
    
    // iterator for dataset extension
    int iter;
    
    template <class T> void append( hid_t fid, std::string name, T & property,  hid_t  mem_space, int nParticles, hid_t type, std::vector<hsize_t> &locator);
};

#endif

