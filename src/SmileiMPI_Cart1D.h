#ifndef SMILEIMPI_CART1D_H
#define SMILEIMPI_CART1D_H

#include <vector>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"

class Species;

class SmileiMPI_Cart1D : public SmileiMPI {
public:
    friend class SmileiIO_Cart1D;


    SmileiMPI_Cart1D( int* argc, char*** argv );
    SmileiMPI_Cart1D(SmileiMPI *smpi);
    virtual ~SmileiMPI_Cart1D();

    virtual void whoami() {
        std::cout << "SmileiMPI_Cart1D" << std::endl;
    }

    //! Create MPI communicator
    virtual void createTopology(PicParams& params);
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    //! exchangeParticles implements particles sorting
    virtual void exchangeParticles(Species* species, int ispec, PicParams* params, int tnum);

    //! Create MPI_Datatype to exchange/sum fields on ghost data
    //! Useless if 1D
    void createType( PicParams& params ) {};

    //! Create MPI_Datatype to exchange all properties of particle in 1 communication
    MPI_Datatype createMPIparticles( Particles* particles, int nbrOfProp );


    virtual void exchangeField ( Field* field );
    virtual void sumField      ( Field* field );

    inline int getProcCoord(int i) {
        return coords_[i];
    }
    inline int getNbrOfProcs(int i) {
        return number_of_procs[i];
    }

    inline bool isWester() {
        return (coords_[0]==0);
    }
    inline bool isEaster() {
        return (coords_[0]==number_of_procs[0]-1);
    }

    int extrem_ranks[1][2];

protected:
    MPI_Comm SMILEI_COMM_1D;

    int ndims_;
    int* number_of_procs;

    // Cartesian ...
    int* coords_;
    int* periods_;
    int reorder_;

    int nbNeighbors_;     // Per direction, ie = 2
    int neighbor_[3][2];	// 

};

#endif

