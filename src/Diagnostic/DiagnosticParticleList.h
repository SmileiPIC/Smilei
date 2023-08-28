#ifndef DIAGNOSTICPARTICLELIST_H
#define DIAGNOSTICPARTICLELIST_H

#include "VectorPatch.h"
#include "Diagnostic.h"

class Patch;
class Params;
class SmileiMPI;


class DiagnosticParticleList : public Diagnostic
{

public :
    //! Default constructor
    DiagnosticParticleList( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, std::string, std::string, unsigned int, OpenPMDparams & );
    //! Default destructor
    ~DiagnosticParticleList() override;
    
    virtual void openFile( Params &params, SmileiMPI *smpi ) = 0;
    
    virtual void closeFile() = 0;
    
    virtual void init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches ) = 0;
    
    bool prepare( int itime ) override;
    
    void run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers ) override;
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        return 0;
    }
    
    //! Get disk footprint of current diagnostic
    virtual uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) = 0;
    
    //! Returns the Particles object of interest in a given patch
    virtual Particles * getParticles( Patch * patch ) = 0;
    
    //! Prepare all HDF5 groups, datasets and spaces
    virtual H5Space * prepareH5( SimWindow *simWindow, SmileiMPI *smpi, int itime, uint32_t nParticles_local, uint64_t nParticles_global, uint64_t offset ) = 0;
    
    //! Close HDF5 groups, datasets and spaces
    virtual void deleteH5() {};
    
    //! Modify the filtered particles
    virtual void modifyFiltered( VectorPatch &, unsigned int ) {};
    
    //! Write extra particles properties
    virtual void writeOther( VectorPatch &, size_t, H5Space *, H5Space * ) {};
    
    //! Fills a buffer with the required particle property
    template<typename T> void fill_buffer( VectorPatch &vecPatches, size_t iprop, std::vector<T> &buffer )
    {
        const size_t nPatches = vecPatches.size();
        std::vector<T> *property = NULL;
        
        #pragma omp barrier
        if( has_filter ) {
            #pragma omp for schedule(runtime)
            for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
                const size_t patch_nParticles = patch_selection[ipatch].size();
                getParticles( vecPatches( ipatch ) )->getProperty( iprop, property );
                size_t i=0;
                size_t j=patch_start[ipatch];
                while( i < patch_nParticles ) {
                    buffer[j] = ( *property )[patch_selection[ipatch][i]];
                    i++;
                    j++;
                }
            }
        } else {
            #pragma omp for schedule(runtime)
            for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
                Particles * p = getParticles( vecPatches( ipatch ) );
                const size_t patch_nParticles = p->numberOfParticles();
                p->getProperty( iprop, property );
                std::copy( property->begin(), property->begin() + patch_nParticles, buffer.begin() + patch_start[ipatch] );
            }
        }
    };

    //! Write a dataset
    virtual void write_scalar_uint64( H5Write * location, std::string name, uint64_t &, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) = 0;
    virtual void write_scalar_short ( H5Write * location, std::string name, short    &, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) = 0;
    virtual void write_scalar_double( H5Write * location, std::string name, double   &, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) = 0;
    virtual void write_component_uint64( H5Write * location, std::string name, uint64_t &, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) = 0;
    virtual void write_component_short ( H5Write * location, std::string name, short    &, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) = 0;
    virtual void write_component_double( H5Write * location, std::string name, double   &, H5Space *file_space, H5Space *mem_space, unsigned int unit_type ) = 0;
    
    //! Index of the species used
    unsigned int species_index_;
    
    //! Name of the species used
    std::string species_name_;
    
protected:
    
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
    
    //! HDF5 locations where attributes must be written
    bool write_id_ = false; H5Write* loc_id_ = nullptr;
    bool write_charge_ = false; H5Write* loc_charge_ = nullptr;
    bool write_weight_ = false; H5Write* loc_weight_ = nullptr;
    bool write_chi_    = false; H5Write* loc_chi_    = nullptr;
    std::vector<bool> write_position_ = {false, false, false}; std::vector<H5Write*> loc_position_ = {nullptr, nullptr, nullptr};
    std::vector<bool> write_momentum_ = {false, false, false}; std::vector<H5Write*> loc_momentum_ = {nullptr, nullptr, nullptr};
    std::vector<bool> write_E_ = {false, false, false}; std::vector<H5Write*> loc_E_ = {nullptr, nullptr, nullptr};
    std::vector<bool> write_B_ = {false, false, false}; std::vector<H5Write*> loc_B_ = {nullptr, nullptr, nullptr};
    std::vector<bool> write_W_ = {false, false, false}; std::vector<H5Write*> loc_W_ = {nullptr, nullptr, nullptr};
    bool interpolate_ = false;
    bool write_any_position_;
    bool write_any_momentum_;
    bool write_any_E_;
    bool write_any_B_;
    bool write_any_W_;
    
    H5Write* loc_birth_time_ = nullptr;
};

#endif

