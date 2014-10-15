#ifndef SMILEIMPI_CART2D_H
#define SMILEIMPI_CART2D_H

#include <vector>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"

class Species;

class SmileiMPI_Cart2D : public SmileiMPI {
public:
    friend class SmileiIO_Cart2D;

    SmileiMPI_Cart2D( int* argc, char*** argv );
    SmileiMPI_Cart2D(SmileiMPI *smpi);
    virtual ~SmileiMPI_Cart2D();

    virtual void whoami() {
        std::cout << "SmileiMPI_Cart2D" << std::endl;
    }

    //! Create MPI communicator
    virtual void createTopology(PicParams& params);
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    //! exchangeParticles = IexchangeParticles
    virtual void exchangeParticles(Species* species, int ispec, PicParams& params, int tnum);

    //! Create MPI_Datatype to exchange/sum fields on ghost data
    void createType( PicParams& params );

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

    inline bool isWestern() {
        return ((coords_[0]==0)&&(periods_[0]==0));
    }
    inline bool isEastern() {
        return ((coords_[0]==number_of_procs[0]-1)&&(periods_[0]==0));
    }
    inline bool isSouthern() {
        return ((coords_[1]==0)&&(periods_[1]==0));
    }
    inline bool isNorthern() {
        return ((coords_[1]==number_of_procs[1]-1)&&(periods_[1]==0));
    }

    int extrem_ranks[2][2];


protected:
    MPI_Comm SMILEI_COMM_2D;

    int ndims_;
    int* number_of_procs;

    // Cartesian ...
    int* coords_;
    int* periods_;
    int reorder_;

    int nbNeighbors_;     // Per direction, ie = 2
    int neighbor_[3][2];	// 

    // MPI_Datatype [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntype_   [2][2][2];
    MPI_Datatype ntypeSum_[2][2][2];


};

#endif
