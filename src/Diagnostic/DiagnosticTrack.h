#ifndef DIAGNOSTICTRACK_H
#define DIAGNOSTICTRACK_H

#include "DiagnosticParticleList.h"

class Patch;
class Params;
class SmileiMPI;

class DiagnosticTrack : public DiagnosticParticleList
{

public :
    //! Default constructor
    DiagnosticTrack( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int, unsigned int, OpenPMDparams & );
    //! Default destructor
    ~DiagnosticTrack() override;
    
    void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
    //! Returns the Particles object of interest in a given patch
    Particles * getParticles( Patch * patch ) override;
    
    //! Prepare all HDF5 groups, datasets and spaces
    H5Space * prepareH5( SimWindow *simWindow, SmileiMPI *smpi, int itime, uint32_t nParticles_local, uint64_t nParticles_global, uint64_t offset ) override;
    
    //! Close HDF5 groups, datasets and spaces
    void deleteH5() override;
    
    //! Modify the filtered particles (apply new ID)
    void modifyFiltered( VectorPatch &, unsigned int ) override;
    
    //! Write a dataset
    void write_scalar_uint64( H5Write * location, std::string name, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_scalar_short ( H5Write * location, std::string name, short    &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_scalar_double( H5Write * location, std::string name, double   &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_component_uint64( H5Write * location, std::string name, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_component_short ( H5Write * location, std::string name, short    &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_component_double( H5Write * location, std::string name, double   &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    
    //! Set a given patch's particles with the required IDs (used at initialization & simWindow)
    void setIDs( Patch * );
    
    //! Set a given particles with the required IDs (used by importParticles)
    void setIDs( Particles & );
    
    //! Last ID assigned to a particle by this MPI domain
    uint64_t latest_Id;
    
    //! Flag to test whether IDs have been set already
    bool IDs_done;
    
private :
    
    H5Write * data_group_;
};

#endif

