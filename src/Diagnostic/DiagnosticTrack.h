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
    ~DiagnosticTrack() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    bool prepare( int timestep ) override;
    
    void run( Patch* patch, int timestep ) {};
    void run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep ) override;
    
    bool write(int timestep) {};
    
    virtual void finish(int, VectorPatch& ) {};
    
    
private :
    //! Flag to test whether IDs have been set already
    bool IDs_done;
    
    //! Index of the species used
    int speciesId_;
    
    //! Size of the diag (total number of particles)
    int nbrParticles_;
     
    //! HDF5 file transfer protocol
    hid_t transfer;
    //! HDF5 file space (dimensions of the array in file)
    hsize_t dims[2];
    //! HDF5 memory space (dimensions of the current particle array in memory)
    hid_t mem_space;
     
    //! Number of spatial dimensions
    int nDim_particle;
    
    //! list of datasets to be added to the file
    std::vector<std::string> datasets;
    //! list of data types for each dataset
    std::vector<hid_t> datatypes;
    
    //! Current particle partition among the patches own by current MPI
    std::vector<int> patch_start;
    
    //! Buffer for the output of double array
    std::vector<double> data_double;
    //! Buffer for the output of short array
    std::vector<short> data_short;
    //! Buffer for the output of uint array
    std::vector<unsigned int> data_uint;
    
};

#endif

