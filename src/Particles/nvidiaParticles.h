
#ifndef NVIDIAPARTICLES_H
#define NVIDIAPARTICLES_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include "Particles.h"

typedef thrust::device_vector<double>::iterator Diter;
typedef thrust::device_vector<short>::iterator  Siter;
// typedef a tuple of these iterators
typedef thrust::tuple<Diter, Diter, Diter, Diter, Diter, Diter, Diter, Siter> IteratorParticles;
// typedef the zip_iterator of this tuple
typedef thrust::zip_iterator<IteratorParticles> ZipIterParts;

class nvidiaParticles : public Particles
{
public:
    //! Constructor for Particle
    nvidiaParticles();

    //! Destructor for Particle
    virtual ~nvidiaParticles() {};

    void initializeDataOnDevice() override;
    void syncGPU() override;
    void syncCPU() override;

    std::vector< thrust::device_vector<double> > nvidia_position_;
    std::vector< thrust::device_vector<double> > nvidia_momentum_;

    //! Weight
    thrust::device_vector<double> nvidia_weight_;

    //! Charge on GPU
    thrust::device_vector<short>  nvidia_charge_;

    //! cell_keys of the particle
    thrust::device_vector<int> nvidia_cell_keys_;

    //! Quantum parameter
    thrust::device_vector<double> nvidia_chi_;

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
    int * getPtrCellKeys() override {
        return thrust::raw_pointer_cast( nvidia_cell_keys_.data() );
    };

    //! Get number of particules
    unsigned int gpu_size() const override
    {
        return gpu_nparts_;
    }

    void extractParticles( Particles* particles_to_move ) override;
    int injectParticles( Particles* particles_to_move ) override;

    int gpu_nparts_;
    int nparts_to_move_;

};

#endif
