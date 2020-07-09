#ifndef DIAGNOSTICTRACK_H
#define DIAGNOSTICTRACK_H

#include "Diagnostic.h"

class Patch;
class Params;
class SmileiMPI;


class DiagnosticTrack : public Diagnostic
{

public :
    //! Default constructor
    DiagnosticTrack( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int, unsigned int, OpenPMDparams & );
    //! Default destructor
    ~DiagnosticTrack() override;
    
    void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    bool prepare( int itime ) override;
    
    void run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers ) override;
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        return 0;
    }
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
    //! Fills a buffer with the required particle property
    template<typename T> void fill_buffer( VectorPatch &vecPatches, unsigned int iprop, std::vector<T> &buffer );
    
    //! Write a scalar dataset with the given buffer
    template<typename T> void write_scalar( H5Write*, std::string, T &, hid_t, H5Space*, H5Space*, unsigned int );
    
    //! Write a vector component dataset with the given buffer
    template<typename T> void write_component( H5Write*, std::string, T &, hid_t, H5Space*, H5Space*, unsigned int );
    
    //! Set a given patch's particles with the required IDs
    void setIDs( Patch * );
    
    //! Set a given particles with the required IDs
    void setIDs( Particles & );
    
    //! Index of the species used
    unsigned int speciesId_;
    
    //! Last ID assigned to a particle by this MPI domain
    uint64_t latest_Id;
    
    //! Flag to test whether IDs have been set already
    bool IDs_done;
    
private :
    
    H5Write *data_group;
    
    //! Number of spatial dimensions
    unsigned int nDim_particle;
    
    //! Current particle partition among the patches own by current MPI
    std::vector<unsigned int> patch_start;
    
    //! Tells whether this diag includes a particle filter
    bool has_filter;
    
    //! Tells whether this diag includes a particle filter
    PyObject *filter;
    
    //! Selection of the filtered particles in each patch
    std::vector<std::vector<unsigned int> > patch_selection;
    
    //! Buffer for the output of double array
    std::vector<double> data_double;
    //! Buffer for the output of short array
    std::vector<short> data_short;
    //! Buffer for the output of uint64 array
    std::vector<uint64_t> data_uint64;
    
    //! Approximate total number of particles
    double npart_total;
    
    //! Number of particles shared among patches in this proc
    uint32_t nParticles_local;
    
    //! Booleans to determine which attributes to write out
    std::vector<bool> write_position;
    std::vector<bool> write_momentum;
    bool write_charge;
    bool write_weight;
    bool write_chi   ;
    std::vector<bool> write_E;
    std::vector<bool> write_B;
    bool interpolate;
    bool write_any_position;
    bool write_any_momentum;
    bool write_any_E;
    bool write_any_B;
    
};

#endif

