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
    Patch2D(Patch2D* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved, bool with_particles);

    void initStep2(Params& params)  ;
    
    //! Destructor for Patch
    ~Patch2D()   ;


    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------

    //! init comm / sum densities
    void initSumField( Field* field, int iDim )  ;
    void reallyinitSumField( Field* field, int iDim )  ;
    //! finalize comm / sum densities
    void finalizeSumField( Field* field, int iDim )  ;
    void reallyfinalizeSumField( Field* field, int iDim )  ;

    //! init comm / exchange fields
    void initExchange( Field* field )  ;
    //! finalize comm / exchange fields
    void finalizeExchange( Field* field )  ;
    //! init comm / exchange fields in direction iDim only
    void initExchange( Field* field, int iDim )  ;
    //! finalize comm / exchange fields in direction iDim only
    void finalizeExchange( Field* field, int iDim )  ;

    // Create MPI_Datatype to exchange fields
    void createType( Params& params )  ;

    //! MPI_Datatype to sum [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_[2][2][2];
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    //!   - +1 : an additional type to exchange clrw lines
    MPI_Datatype ntype_[3][2][2];
    // Use a buffer per direction to exchange data before summing
    Field2D buf[2][2];



};

#endif
