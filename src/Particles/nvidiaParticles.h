
#ifndef NVIDIAPARTICLES_H
#define NVIDIAPARTICLES_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include "Particles.h"

typedef thrust::host_vector<double>::iterator Diter;
typedef thrust::host_vector<short>::iterator  Siter;
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

    void initGPU() override;
    void syncGPU() override;
    void syncCPU() override;

    std::vector< thrust::host_vector<double> > nvidiaPosition;
    std::vector< thrust::host_vector<double> > nvidiaMomentum;
    thrust::host_vector<double> nvidiaWeight;
    thrust::host_vector<short>  nvidiaCharge;

    //! cell_keys of the particle
    thrust::host_vector<int> nvidia_cell_keys;

    double* getPtrPosition( int idim ) override {
        return thrust::raw_pointer_cast( nvidiaPosition[idim].data() );
    };
    double* getPtrMomentum( int idim ) override {
        return thrust::raw_pointer_cast( nvidiaMomentum[idim].data() );
    };
    double* getPtrWeight() override {
        return thrust::raw_pointer_cast( nvidiaWeight.data() );
    };
    short * getPtrCharge() override {
        return thrust::raw_pointer_cast( nvidiaCharge.data() );
    };
    int * getPtrCellKeys() override {
        return thrust::raw_pointer_cast( nvidia_cell_keys.data() );
    };

    //! Get number of particules
    inline unsigned int gpu_size() const
    {
        return gpu_nparts_;
    }

    void extractParticles( Particles* particles_to_move ) override;
    int injectParticles( Particles* particles_to_move ) override;

    int gpu_nparts_;
    int nparts_to_move_;

};

#endif
