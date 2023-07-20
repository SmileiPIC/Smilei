#ifndef DIAGNOSTICNEWPARTICLES_H
#define DIAGNOSTICNEWPARTICLES_H

#include "DiagnosticParticleList.h"

class Patch;
class Params;
class SmileiMPI;

class DiagnosticNewParticles : public DiagnosticParticleList
{

public :
    //! Default constructor
    DiagnosticNewParticles( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int, unsigned int, OpenPMDparams & );
    //! Default destructor
    ~DiagnosticNewParticles() override;
    
    void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
    //! Returns the Particles object of interest in a given patch
    Particles * getParticles( Patch * patch ) override;
    
    //! Prepare all HDF5 groups, datasets and spaces
    H5Space * prepareH5( SimWindow *simWindow, SmileiMPI *smpi, int itime, uint32_t nParticles_local, uint64_t nParticles_global, uint64_t offset ) override;
    
    //! Write extra particles properties
    void writeOther( VectorPatch &, size_t, H5Space *, H5Space * ) override;
    
    //! Write a dataset
    void write_scalar_uint64( H5Write * location, std::string name, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_scalar_short ( H5Write * location, std::string name, short    &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_scalar_double( H5Write * location, std::string name, double   &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_component_uint64( H5Write * location, std::string name, uint64_t &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_component_short ( H5Write * location, std::string name, short    &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    void write_component_double( H5Write * location, std::string name, double   &buffer, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) override;
    
    H5Write * newDataset( H5Write &group, std::string name, hid_t dtype, H5Space &file_space, unsigned int unit_type ) {
        H5Write * d = new H5Write( &group, name, dtype, &file_space );
        openPMD_->writeComponentAttributes( *d, unit_type );
        return d;
    };
    
    
private :
    
    //! Number of particles previously written in the file
    hsize_t nParticles_written = 0;
    
    //! Number of times the file was previously written
    hsize_t nTimes_written = 0;
    
    //! Storage for the number of particles at each output iteration
    H5Write * iteration_npart_ = nullptr;
};

#endif

