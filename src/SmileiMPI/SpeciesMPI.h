#ifndef SPECIESMPI_H
#define SPECIESMPI_H

#include <mpi.h>

#include "Particles.h"

class SpeciesMPI {
 public:
    SpeciesMPI();
    ~SpeciesMPI();
    void init();

    //friend class SmileiMPI_Cart2D;


    Particles patchVectorRecv[2][2];
    Particles patchVectorSend[2][2];
    //Particles cornerVectorRecv[2][2];
    //Particles cornerVectorSend[2][2];

    std::vector<int> patch_buff_index_send[2][2];
    //std::vector<int> corner_buff_index_send[2][2];
    int patch_buff_index_send_sz[2][2];
    //int corner_buff_index_send_sz[2][2];
    int patch_buff_index_recv_sz[2][2];
    //int corner_buff_index_recv_sz[2][2];

    MPI_Request patch_srequest[2][2];
    MPI_Request patch_rrequest[2][2];
    //MPI_Request corner_srequest[2][2];
    //MPI_Request corner_rrequest[2][2];

};

#endif

