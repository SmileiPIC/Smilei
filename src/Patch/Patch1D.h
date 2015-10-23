#ifndef PATCH1D_H
#define PATCH1D_H

#include "Patch.h"

class Diagnostic;
class DiagnosticScalar;
class SimWindow;

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch1D : public Patch
{
public:
    //! Constructor for Patch
    Patch1D(PicParams& params, DiagParams &diag_params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved);

    //! Destructor for Patch
    virtual ~Patch1D() {};

    virtual void initSumField( Field* field, int iDim );
    virtual void finalizeSumField( Field* field, int iDim );

    virtual void initExchange( Field* field );
    virtual void finalizeExchange( Field* field );
    virtual void initExchange( Field* field, int iDim );
    virtual void finalizeExchange( Field* field, int iDim );

    virtual void createType( PicParams& params );
    //! MPI_Datatype to exchange [ndims_][iDim=0 prim/dial]
    MPI_Datatype ntypeSum_[2][2];

    MPI_Datatype ntype_[2][2];
    
    // Use a buffer per direction to exchange data before summing
    Field1D buf[1][2];

    inline bool isWestern()  { return locateOnBorders(0, 0); }
    inline bool isEastern()  { return locateOnBorders(0, 1); }
    /*inline bool isSouthern() { return locateOnBorders(1, 0); }
      inline bool isNorthern() { return locateOnBorders(1, 1); }*/



    //! Return MPI rank of this->hrank +/- 1
    //! Should be replaced by an analytic formula
    inline int getMPIRank(int hrank_pm1) {
	if  (hrank_pm1 == neighbor_[0][0]) return MPI_neighbor_[0][0];
	else if  (hrank_pm1 == neighbor_[0][1]) return MPI_neighbor_[0][1];
	//else if  (hrank_pm1 == neighbor_[1][0]) return MPI_neighbor_[1][0];
	//else if  (hrank_pm1 == neighbor_[1][1]) return MPI_neighbor_[1][1];
	else
	    return MPI_PROC_NULL;
    }



};

#endif
