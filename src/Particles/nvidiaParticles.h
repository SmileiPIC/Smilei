// -----------------------------------------------------------------------------
//
//! \file nvidiaParticles.h
//
//! \brief contains the nvidiaParticles class description
//
//! The nvidiaParticles inherits from the Particles class to deal with GPUs.
//! It uses NVIDIA/AMD thrust::device_vector instead of std::vector
//
// -----------------------------------------------------------------------------

#ifndef NVIDIAPARTICLES_H
#define NVIDIAPARTICLES_H

#include <thrust/device_vector.h>

#include "Params.h"
#include "Particles.h"

/*! \class nvidiaParticles
    \brief Particle class for GPUs
*/
class nvidiaParticles : public Particles
{
public:
    //! Constructor for Particle
    nvidiaParticles(const Params& parameters);

    //! Destructor for nvidiaParticles
    virtual ~nvidiaParticles() {};

    //! Allocate the right amount of position and momentum dimensions
    void allocateDimensions( unsigned int nDim );

    //! Reserve space for particle_count particles. Must be called after 
    //! allocateDimensions()
    void reserveParticles( unsigned int particle_count );

    //! Allocate particle_count particles. Must be called after
    //! allocateDimensions()
    void allocateParticles( unsigned int particle_count );

    //! Reset Particles vectors
    void deviceClear();

    //! Initialize the particle properties on devide as a mirror of the host definition
    // 
    void initializeDataOnDevice() override;
    
    //! Send the particles from host to device
    void syncGPU() override;
    
    //! Update the particles from device to host
    void syncCPU() override;

    //! Position vector on device
    std::vector< thrust::device_vector<double> > nvidia_position_;
    
    //! Momentum vector on device 
    std::vector< thrust::device_vector<double> > nvidia_momentum_;

    //! Weight
    thrust::device_vector<double> nvidia_weight_;

    //! Charge on GPU
    thrust::device_vector<short>  nvidia_charge_;

    //! cell_keys of the particle
    thrust::device_vector<int> nvidia_cell_keys_;

    //! Quantum parameter
    thrust::device_vector<double> nvidia_chi_;

    //! Monte-Carlo parameter
    thrust::device_vector<double> nvidia_tau_;

    //! List of double* arrays
    std::vector< thrust::device_vector<double>* > nvidia_double_prop_;
    
    //! List of short* arrays
    std::vector< thrust::device_vector<short>* > nvidia_short_prop_;

    double* getPtrPosition( int idim ) override {
        return thrust::raw_pointer_cast( nvidia_position_[idim].data() );
    };
    double* getPtrMomentum( int idim ) override {
        return thrust::raw_pointer_cast( nvidia_momentum_[idim].data() );
    };
    double* getPtrWeight() override {
        return thrust::raw_pointer_cast( nvidia_weight_.data() );
    };
    short * getPtrCharge() override {
        return thrust::raw_pointer_cast( nvidia_charge_.data() );
    };
    double * getPtrChi() override {
        return thrust::raw_pointer_cast( nvidia_chi_.data() );
    };
    double * getPtrTau() override {
        return thrust::raw_pointer_cast( nvidia_tau_.data() );
    };
    int * getPtrCellKeys() override {
        return thrust::raw_pointer_cast( nvidia_cell_keys_.data() );
    };

    //! Get number of particles
    unsigned int gpu_size() const override
    {
        return gpu_nparts_;
    }

    // -----------------------------------------------------------------------------
    //! Extract particles from the Particles object and put 
    //! them in the Particles object `particles_to_move`
    // -----------------------------------------------------------------------------
    void extractParticles( Particles* particles_to_move ) override;
    
    // -----------------------------------------------------------------------------
    //! Erase particles leaving the patch object on device and returns the number of particle removed
    // -----------------------------------------------------------------------------
    int eraseLeavingParticles() override;
    
    // -----------------------------------------------------------------------------
    //! Inject particles from particles_to_move into *this and return he number of particle added
    // -----------------------------------------------------------------------------
    int injectParticles( Particles* particles_to_inject ) override;

    // ---------------------------------------------------------------------------------------------------------------------
    //! Create n_additional_particles new particles at the end of vectors
    //! Fill the new elements with 0
    // ---------------------------------------------------------------------------------------------------------------------
    void createParticles( int n_additional_particles ) override;

    //! See the Particles class for documentation.
    void importAndSortParticles( const Particles* particles_to_inject ) override;

    // Number of particles on device
    int gpu_nparts_;
};

#endif
