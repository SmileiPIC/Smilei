#ifndef ASYNCMPIBUFFERS_H
#define ASYNCMPIBUFFERS_H

#include <mpi.h>
#include <vector>
#include <complex>

#include "Particles.h"

class Field;
class Patch;
class SmileiMPI;

class AsyncMPIbuffers
{
public:
    AsyncMPIbuffers();
    ~AsyncMPIbuffers();
    
    void allocate( unsigned int nDim_field );
    
    void defineTags( Patch *patch, SmileiMPI *smpi, int tag ) ;
    
    //! ndim vectors of 2 sent requests (1 per direction)
    std::vector< std::vector<MPI_Request> > srequest;
    //! ndim vectors of 2 received requests (1 per direction)
    std::vector< std::vector<MPI_Request> > rrequest;
    std::vector< double >  buf[3][2];
    std::vector< std::complex<double> >  ibuf[3][2];
    
    std::vector< std::vector<int> > send_tags_, recv_tags_;
    
};

class SpeciesMPIbuffers : public AsyncMPIbuffers
{
public:
    SpeciesMPIbuffers();
    ~SpeciesMPIbuffers();
    
    void allocate( Params &params, Patch *patch ) ;
    
    //! ndim vectors of 2 sent packets of particles (1 per direction)
    std::vector< std::vector<Particles* > > partRecv;
    //! ndim vectors of 2 received packets of particles (1 per direction)
    std::vector< std::vector<Particles* > > partSend;
    
    //! ndim vectors of 2 numbers of particles to send (1 per direction)
    std::vector< std::vector< unsigned int > > partSendSize;
    //! ndim vectors of 2 numbers of particles to receive (1 per direction)
    std::vector< std::vector< unsigned int > > partRecvSize;
    
};

#endif

