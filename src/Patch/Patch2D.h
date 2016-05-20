#ifndef PATCH2D_H
#define PATCH2D_H

#include "Patch.h"
#include "Field2D.h"

class SimWindow;

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch2D : public Patch
{
public:
    //! Constructor for Patch
    Patch2D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved);
    //! Cloning Constructor for Patch
    Patch2D(Patch2D* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved);

    void initStep2(Params& params) override;
    
    //! Destructor for Patch
    ~Patch2D() override {};


    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------

    //! init comm / sum densities
    void initSumField( Field* field, int iDim ) override;
    //! finalize comm / sum densities
    void finalizeSumField( Field* field, int iDim ) override;

    //! init comm / exchange fields
    void initExchange( Field* field ) override;
    //! finalize comm / exchange fields
    void finalizeExchange( Field* field ) override;
    //! init comm / exchange fields in direction iDim only
    void initExchange( Field* field, int iDim ) override;
    //! finalize comm / exchange fields in direction iDim only
    void finalizeExchange( Field* field, int iDim ) override;

    // Create MPI_Datatype to exchange fields
    void createType( Params& params ) override;

    //! MPI_Datatype to sum [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_[2][2][2];
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    //!   - +1 : an additional type to exchange clrw lines
    MPI_Datatype ntype_[3][2][2];
    // Use a buffer per direction to exchange data before summing
    Field2D buf[2][2];

    //! Should be pure virtual, see Patch class
    inline bool isWestern()  { return locateOnBorders(0, 0); }
    //! Should be pure virtual, see Patch class
    inline bool isEastern()  { return locateOnBorders(0, 1); }
    //! Should be pure virtual, see Patch class
    inline bool isSouthern() { return locateOnBorders(1, 0); }
    //! Should be pure virtual, see Patch class
    inline bool isNorthern() { return locateOnBorders(1, 1); }

    //! Return MPI rank of this->hrank +/- 1
    //! Should be replaced by an analytic formula
    inline int getMPIRank(int hrank_pm1) override {
	if  (hrank_pm1 == neighbor_[0][0]) return MPI_neighbor_[0][0];
	else if  (hrank_pm1 == neighbor_[0][1]) return MPI_neighbor_[0][1];
	else if  (hrank_pm1 == neighbor_[1][0]) return MPI_neighbor_[1][0];
	else if  (hrank_pm1 == neighbor_[1][1]) return MPI_neighbor_[1][1];
	else
	    return MPI_PROC_NULL;
    }



};

#endif
