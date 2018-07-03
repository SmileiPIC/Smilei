#ifndef ASYNCMPIBUFFERS_H
#define ASYNCMPIBUFFERS_H

#include <mpi.h>
#include <vector>
#include <complex>

#include "Particles.h"

class Field;
class Patch;

class AsyncMPIbuffers {
public:
    AsyncMPIbuffers();
    ~AsyncMPIbuffers();

    virtual void allocate(unsigned int nDim_field);

    virtual void allocate(unsigned int nDim_field, Field* f, std::vector<unsigned int>& oversize);
    virtual void iallocate(unsigned int nDim_field, Field* f, std::vector<unsigned int>& oversize);
    void defineTags(Patch* patch, int tag ) ;
    
    //! ndim vectors of 2 sent requests (1 per direction) 
    std::vector< std::vector<MPI_Request> > srequest;
    //! ndim vectors of 2 received requests (1 per direction) 
    std::vector< std::vector<MPI_Request> > rrequest;
    std::vector< double >  buf[3][2];
    std::vector< std::complex<double> >  ibuf[3][2];

    std::vector< std::vector<int> > send_tags_, recv_tags_;

};

class SpeciesMPIbuffers : public AsyncMPIbuffers {
public:
    SpeciesMPIbuffers();
    ~SpeciesMPIbuffers();

    void allocate(unsigned int nDim_field) ;

    //! ndim vectors of 2 sent packets of particles (1 per direction) 
    std::vector< std::vector<Particles > > partRecv;
    //! ndim vectors of 2 received packets of particles (1 per direction) 
    std::vector< std::vector<Particles > > partSend;

    //! ndim vectors of 2 vectors of index particles to send (1 per direction) 
    //!   - not sent
    //    - used to sort Species::indexes_of_particles_to_exchange built in Species::dynamics
    std::vector< std::vector< std::vector<int> > > part_index_send;
    //! ndim vectors of 2 numbers of particles to send (1 per direction) 
    std::vector< std::vector< unsigned int > > part_index_send_sz;
    //! ndim vectors of 2 numbers of particles to receive (1 per direction) 
    std::vector< std::vector< unsigned int > > part_index_recv_sz;

};

#endif

