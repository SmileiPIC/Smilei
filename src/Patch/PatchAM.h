#ifndef PATCHAM_H
#define PATCHAM_H

#include "Patch.h"
#include "cField2D.h"

class SimWindow;

//! Class Patch : sub MPI domain
//!     Collection of patch = MPI domain
class PatchAM : public Patch
{
public:
    //! Constructor for Patch
    PatchAM( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved );
    //! Cloning Constructor for Patch
    PatchAM( PatchAM *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles );
    
    void initStep2( Params &params, DomainDecomposition *domain_decomposition ) override final;
    
    //! Destructor for Patch
    ~PatchAM() override  final;
    
    //! Return the volume (or surface or length depending on simulation dimension)
    //! of one cell at the position of a given particle
    double getPrimalCellVolume( Particles *p, unsigned int ipart, Params &params ) override final
    {
        ERROR( "getPrimalCellVolume not implemented in geometry AM" );
        return cell_volume;
    };
    
    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------
    
    //! init comm / sum densities
    void initSumFieldComplex( Field *field, int iDim, SmileiMPI *smpi ) override final;
    //! finalize comm / sum densities
    void finalizeSumFieldComplex( Field *field, int iDim ) override final;
    //! init comm / sum densities
    void initSumField( Field *field, int iDim, SmileiMPI *smpi ) override final;
    //! finalize comm / sum densities
    void finalizeSumField( Field *field, int iDim ) override final;
    
    //! init comm / exchange fields in direction iDim only
    void initExchange( Field *field, int iDim, SmileiMPI *smpi ) override final;
    //! init comm / exchange complex fields in direction iDim only
    void initExchangeComplex( Field *field, int iDim, SmileiMPI *smpi ) override final;
    //! finalize comm / exchange fields in direction iDim only
    void finalizeExchange( Field *field, int iDim ) override final;
    //! finalize comm / exchange fields in direction iDim only
    void finalizeExchangeComplex( Field *field, int iDim ) override final;
    
    void exchangeField_movewin( Field* field, int clrw ) override final;
    
    // Create MPI_Datatype to exchange fields
    void createType( Params &params ) override final;
    void createType2( Params &params ) override final;
    void cleanType() override final;
    
    //! MPI_Datatype to sum [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_[2][2][2];
    //! MPI_Datatype to sum [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_complex_[2][2][2];
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntype_[3][2][2];
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntype_complex_[3][2][2];

    //! Inverse r coordinate
    std::vector<double> invR, invRd;
    
    
    
};

#endif
