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


////////////////////////////////////////////////////////////////////////////////
// nvidiaParticles definition
////////////////////////////////////////////////////////////////////////////////

/*! \class nvidiaParticles
    \brief Particle class for GPUs
*/
class nvidiaParticles : public Particles
{
public:
    //! Constructor for Particle
    nvidiaParticles( const Params& parameters, const Patch& a_parent_patch );

    //! Destructor for nvidiaParticles
    ~nvidiaParticles();

    //! Allocate the right amount of position and momentum dimensions
    void resizeDimensions( unsigned int nDim );

    //! Reserve space for (particle_count * growth_factor) particles only if 
    //! particle_count >= deviceCapacity(). Must be called after
    //! allocateDimensions()
    void softReserve( unsigned int particle_count, float growth_factor = 1.3F );

    //! Reserve space for particle_count particles. Must be called after
    //! allocateDimensions()
    void reserve( unsigned int particle_count );

    //! Allocate particle_count particles. Must be called after
    //! allocateDimensions()
    //! Set the size (deviceSize) of nvidiaParticles to particle_count.
    //!
    void resize( unsigned int particle_count );

    //! Assures that the memory holden by the nvidia_[position|momentum|weight|
    //! charge|chi|tau|cell_keys]_ is freed. This is not something you can
    //! achieve via a naive resize.
    //! The pointers in nvidia_[double|short]_prop_ are not invalidated.
    //!
    void free();

    //! Resize Particle vectors on device
    void deviceResize(unsigned int new_size);

    //! Remove all particles
    void deviceClear();
    
    //! Reset cell_keys to default value
    void resetCellKeys();
    
    //! Initialize the particle properties on device as a mirror of the host definition
    void initializeDataOnDevice() override;

    //! Initialize particle IDs on device
    void initializeIDsOnDevice() override;

    //! Send the particles from host to device
    void copyFromHostToDevice() override;
    
    //! Update the particles from device to host
    void copyFromDeviceToHost( bool copy_keys = false ) override;

    unsigned int deviceCapacity() const override;

    inline unsigned int deviceSize() const override
    {
        return gpu_nparts_;
    }

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
    uint64_t * getPtrId() override {
        return thrust::raw_pointer_cast( nvidia_id_.data() );
    };

    void swapPosition( int idim, thrust::device_vector<double> &new_vector ) {
        nvidia_position_[idim].swap( new_vector );
    };
    void swapMomentum( int idim, thrust::device_vector<double> &new_vector ) {
        nvidia_momentum_[idim].swap( new_vector );
    };
    void swapWeight( thrust::device_vector<double> &new_vector ) {
        nvidia_weight_.swap( new_vector );
    };
    void swapChi( thrust::device_vector<double> &new_vector ) {
        nvidia_chi_.swap( new_vector );
    };
    void swapCharge( thrust::device_vector<short> &new_vector ) {
        nvidia_charge_.swap( new_vector );
    };
    void swapTau( thrust::device_vector<double> &new_vector ) {
        nvidia_tau_.swap( new_vector );
    };
    void swapCellKeys( thrust::device_vector<int> &new_vector ) {
        nvidia_cell_keys_.swap( new_vector );
    };
    void swapId( thrust::device_vector<uint64_t> &new_vector ) {
        nvidia_id_.swap( new_vector );
    };
    

    // -----------------------------------------------------------------------------
    //! Move leaving particles to the buffers
    // -----------------------------------------------------------------------------
    void copyLeavingParticlesToBuffer( Particles* buffer ) override;
    
    template<typename Predicate>
    void copyParticlesByPredicate( Particles* buffer, Predicate pred );
    
    // -----------------------------------------------------------------------------
    //! Erase particles leaving the patch object on device and returns the number of particle removed
    // -----------------------------------------------------------------------------
    int eraseLeavingParticles() override;
    
    template<typename Predicate>
    int eraseParticlesByPredicate( Predicate pred );
    
    // -----------------------------------------------------------------------------
    //! Inject particles from particles_to_inject into *this and return the number of particle added
    // -----------------------------------------------------------------------------
    int injectParticles( Particles* particles_to_inject ) override;

    // ---------------------------------------------------------------------------------------------------------------------
    //! Create n_additional_particles new particles at the end of vectors
    //! Fill the new elements with 0
    // ---------------------------------------------------------------------------------------------------------------------
    void createParticles( int n_additional_particles ) override;

    //! See the Particles class for documentation.
    void importAndSortParticles( Particles* particles_to_inject ) override;

protected:
    //! Redefine first_index and last_index according to the binning algorithm
    //! used on GPU.
    //! We simulate one bin on the host and use multiple bin on the GPU.
    //!
    //! NOTE: The bin partitioning initially implemented in Smilei works
    //! decently on CPU. The typical binning is done by sorting particles
    //! according to their x coordinate.
    //!
    //! On GPU we need large volume of data to process. We'd like a sorting in
    //! cluster/tiles such that it would be possible to load the fields of such
    //! cluster into GPU shared/LDS memory.
    //! To represent which particle is in which bin, we re-use last_index. This
    //! array contains "n_cell_per_patch/getGPUClusterCellVolume()" values, each
    //! representing the end of the bin. ie: last_index[0] is the last particle
    //! index (excluded), of bin 0 and the first particle (included) of bin 1.
    //! The last value of last_index is the number of particle of a given
    //! species in this patch.
    //!
    //! One should always use first_index.size() to get the number of "host
    //! visible bin" but note that it may be different than the "true"
    //! number of bin used on the device. It's done this way only to trick Smilei
    //! into passing large particle quantities to the operators because they know
    //! best how to process the data (in bin or not in bin). The Smilei
    //! "structural" code should not have to deal with the fine details (most of
    //! the time).
    //!
    //! In GPU mode, there is only one "host visible bin", it contains all the
    //! particles of a species on a patch.
    //!
    //! In GPU mode, the new interface for first_index and last_index is:
    //! first_index.size() : number of bin on the host (always 1)
    //! last_index.size()  : number of bin on the device (should be seldom used
    //!                      on the host!)
    //! last_index.back(), last_index[last_index.size()-1] and last_index[0]:
    //! number of particles (NOTE: last_index[0] is mandatory to simulate one bin
    //! on the host).
    //!
    //! Everything else is UNDEFINED (!), dont use it. This workaround with
    //! first_index/last_index is not pretty but requires few modifications to
    //! Smilei.
    //!
    //! The cluster width is hard coded in Params::getGPUClusterWidth.
    //!
    //! If the function succeed, last_index is allocated on GPU.
    //!
    //! returns -1 on error (binning is not supported for this Particles
    //! object).
    //!
    int prepareBinIndex();

    //! Set last_index.back() and last_index[0] to match the number of GPU
    //! particle (deviceSize()).
    //!
    void setHostBinIndex() override;

    //! Memcpy of the particle at the end. No sorting or binning.
    //!
    void naiveImportAndSortParticles( nvidiaParticles* particles_to_inject );

    //! Position vector on device
    std::vector<thrust::device_vector<double>> nvidia_position_;

    //! Momentum vector on device
    std::vector<thrust::device_vector<double>> nvidia_momentum_;

    //! Weight
    thrust::device_vector<double> nvidia_weight_;

    //! Charge on GPU
    thrust::device_vector<short> nvidia_charge_;

    //! cell_keys of the particle
    thrust::device_vector<int> nvidia_cell_keys_;

    //! Quantum parameter
    thrust::device_vector<double> nvidia_chi_;

    //! Monte-Carlo parameter
    thrust::device_vector<double> nvidia_tau_;

    //! Particle IDs
    thrust::device_vector<uint64_t> nvidia_id_;

    //! List of double* arrays
    std::vector<thrust::device_vector<double>*> nvidia_double_prop_;

    //! List of short* arrays
    std::vector<thrust::device_vector<short>*> nvidia_short_prop_;

    const Params* parameters_;
    //! We are interested in having the patch coordinates. This allows us to
    //! compute a bin index relative to the patch which in turns, makes the bin
    //! index computation much easier.
    const Patch* parent_patch_;

    //! Number of particles on device
    int gpu_nparts_;
};

#endif
