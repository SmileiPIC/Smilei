#ifndef ASYNCMPIBUFFERS_H
#define ASYNCMPIBUFFERS_H

#include <mpi.h>
#include <vector>
#include <array>

#include "Particles.h"


class AsyncMPIbuffers {
public:
    AsyncMPIbuffers();
    ~AsyncMPIbuffers();

    virtual void allocate(int nDim_field);
    
    //! ndim vectors of 2 sent requests (1 per direction) 
    std::vector< std::vector<MPI_Request> > srequest;
    //! ndim vectors of 2 received requests (1 per direction) 
    std::vector< std::vector<MPI_Request> > rrequest;
};

class SpeciesMPIbuffers : public AsyncMPIbuffers {
public:
    SpeciesMPIbuffers();
    ~SpeciesMPIbuffers();

    void allocate(int nDim_field) override;

    //! ndim vectors of 2 sent packets of particles (1 per direction) 
    std::vector< std::vector<Particles > > partRecv;
    //! ndim vectors of 2 received packets of particles (1 per direction) 
    std::vector< std::vector<Particles > > partSend;

    //! ndim vectors of 2 vectors of index particles to send (1 per direction) 
    //!   - not sent
    //    - used to sort Species::indexes_of_particles_to_exchange built in Species::dynamics
    std::vector< std::vector< std::vector<int> > > part_index_send;
    //! ndim vectors of 2 numbers of particles to send (1 per direction) 
    std::vector< std::vector< int > > part_index_send_sz;
    //! ndim vectors of 2 numbers of particles to receive (1 per direction) 
    std::vector< std::vector< int > > part_index_recv_sz;

};

#endif

