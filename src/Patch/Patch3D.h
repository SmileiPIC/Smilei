#ifndef PATCH3D_H
#define PATCH3D_H

#include "Patch.h"
#include "Field3D.h"

class SimWindow;

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch3D : public Patch
{
public:
    //! Constructor for Patch
    Patch3D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved);
    //! Cloning Constructor for Patch
    Patch3D(Patch3D* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved, bool with_particles);

    void initStep2(Params& params) override final;
    
    //! Destructor for Patch
    ~Patch3D() override  final;


    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------

    //! init comm / sum densities
    void initSumField( Field* field, int iDim ) override final;
    void reallyinitSumField( Field* field, int iDim ) override final;
    //! finalize comm / sum densities
    void finalizeSumField( Field* field, int iDim ) override final;
    void reallyfinalizeSumField( Field* field, int iDim ) override final;

    //! init comm / exchange fields
    void initExchange( Field* field ) override final;
    //! finalize comm / exchange fields
    void finalizeExchange( Field* field ) override final;
    //! init comm / exchange fields in direction iDim only
    void initExchange( Field* field, int iDim ) override final;
    //! finalize comm / exchange fields in direction iDim only
    void finalizeExchange( Field* field, int iDim ) override final;


};

#endif
