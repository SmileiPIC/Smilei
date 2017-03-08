#ifndef DIAGNOSTICTRACK_H
#define DIAGNOSTICTRACK_H

#include "Diagnostic.h"

class Patch;
class Params;
class SmileiMPI;


class DiagnosticTrack : public Diagnostic {

public :
    //! Default constructor
    DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, int );
    //! Cloning constructor
    DiagnosticTrack(DiagnosticTrack* track, Patch* patch);
    //! Default destructor
    ~DiagnosticTrack() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    bool prepare( int itime ) override;
    
    void run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime ) override;
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        return 0;
    }
    
    //! Last ID assigned to a particle by this MPI domain
    uint64_t latest_Id;
    
    //! Index of the species used
    int speciesId_;
    

private :
    
    //! Flag to test whether IDs have been set already
    bool IDs_done;
    
    //! HDF5 file transfer protocol
    hid_t transfer;
     
    //! Number of spatial dimensions
    unsigned int nDim_particle;
    
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
    //! Buffer for the output of uint64 array
    std::vector<uint64_t> data_uint64;
    
};

#endif

