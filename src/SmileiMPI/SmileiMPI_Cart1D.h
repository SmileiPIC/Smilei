#ifndef SMILEIMPI_CART1D_H
#define SMILEIMPI_CART1D_H

#include <vector>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"

class Species;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPI_Cart1D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI_Cart1D : public SmileiMPI {
public:
    friend class SmileiIO_Cart1D;


    //! Create intial MPI environment
    SmileiMPI_Cart1D( int* argc, char*** argv );
    //! Create MPI environment for the data geometry from 
    //! \param smpi the initil MPI environment
    SmileiMPI_Cart1D(SmileiMPI *smpi);
    //! Destructor for SmileiMPI
    virtual ~SmileiMPI_Cart1D();

    //! Create MPI communicator
    virtual void createTopology(Params& params);
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    virtual void exchangeParticles(Species* species, Params& params, int tnum, int iDim);

    //! Create MPI_Datatype to exchange/sum fields on ghost data
    //! Useless if 1D, data are contigous
    void createType( Params& params ) {};

    //! Create MPI_Datatype to exchange all properties of particle in 1 communication
    MPI_Datatype createMPIparticles( Particles* particles, int nbrOfProp );


    //! Basic method to exchange a field,
    virtual void exchangeField ( Field* field );
    //! Basic method to exchange several cells of field toward the X direction
    virtual void exchangeField_movewin ( Field* field, int clrw );
    //! Basic method to sum a field
    virtual void sumField      ( Field* field );

    //! Return coordinates in the cartesian MPI communicator
    //! \param i direction
    inline int getProcCoord(int i) {
        return coords_[i];
    }
    //! Return number of MPI process in the cartesian MPI communicator
    //! \param i direction
    inline int getNbrOfProcs(int i) {
        return number_of_procs[i];
    }

    //! Identify western MPI process, for boundary condition
    inline bool isWestern() {
        return ((coords_[0]==0)&&(periods_[0]==0));

    }
    //! Identify eastern MPI process, for boundary condition
    inline bool isEastern() {
        return ((coords_[0]==number_of_procs[0]-1)&&(periods_[0]==0));
    }
    //! Identify corner MPI ranks (1D, 2 sides) 
    int extrem_ranks[1][2];

protected:
    //! 1D Cartesian communicator
    MPI_Comm SMILEI_COMM_1D;
    //! Number of dimensions
    int ndims_;
    //! Number of MPI process per direction in the cartesian topology
    int* number_of_procs;

    //! Array of coordinates in the cartesian topology
    int* coords_;
    //! Periodicity of the geometry
    int* periods_;
    //! Reorder MPI rank (not)
    int reorder_;
    //! Number of neighbors per directions (=2)
    int nbNeighbors_;
    //! Id of neighbors, per direction (up to 3), per side (2)
    int neighbor_[3][2];

};

#endif

