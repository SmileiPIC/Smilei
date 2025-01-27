#ifndef PATCHAM_H
#define PATCHAM_H

#include "Patch.h"
#include "cField2D.h"

class SimWindow;

//! Class Patch : sub MPI domain
//!     Collection of patch = MPI domain
class PatchAM final : public Patch
{
public:
    //! Constructor for Patch
    PatchAM( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved );
    //! Cloning Constructor for Patch
    PatchAM( PatchAM *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles );
    
    void initStep2( Params &params, DomainDecomposition *domain_decomposition ) override final;
    void initInvR( Params &params );
    
    //! Destructor for Patch
    ~PatchAM() override  final;
    
    //! Return the volume (or surface or length depending on simulation dimension)
    //! of one cell at the position of a given particle
    double getPrimalCellVolume( Particles *p, unsigned int ipart, Params &params ) override final
    {
        double factor = 1.;
        
        double halfcell = 0.5 * params.cell_length[0];
        if( p->position(0,ipart) - getDomainLocalMin(0) < halfcell
         || getDomainLocalMax(0) - p->position(0,ipart) < halfcell ) {
             factor *= 0.5;
        }
        
        halfcell = 0.5 * params.cell_length[1];
        double radius = sqrt(p->position(1,ipart)*p->position(1,ipart) + p->position(2,ipart)*p->position(2,ipart));
        if( radius - getDomainLocalMin(1) < halfcell
          || getDomainLocalMax(1) - radius < halfcell ) {
             factor *= 0.5;
        }

        return factor * cell_volume * radius;
    };
    
    //! Given several arrays (x,y,z), return indices of points in patch
    std::vector<unsigned int> indicesInDomain( double **position, unsigned int n_particles ) override
    {
        std::vector<unsigned int> indices( 0 );
        for( unsigned int ip = 0; ip < n_particles; ip++ ) {
            double distance =  sqrt( position[1][ip]*position[1][ip]+position[2][ip]*position[2][ip] );
            if( position[0][ip] >= getDomainLocalMin( 0 ) && position[0][ip] < getDomainLocalMax( 0 )
             && distance >= getDomainLocalMin( 1 ) && distance < getDomainLocalMax( 1 ) ) {
                indices.push_back( ip );
            }
        }
        return indices;
    };
    
    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------
    
    //! init comm / sum densities
    void initSumFieldComplex( Field *field, int iDim, SmileiMPI *smpi ) override final;
    
    void exchangeField_movewin( Field* field, int clrw ) override final;
    
    // Create MPI_Datatype to exchange fields
    void createType2( Params &params ) override final;
    void cleanType() override final;
    
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntype_complex_[2][2];

    //! Inverse r coordinate
    std::vector<double> invR, invRd;
    
    //! Update Poyting quantities depending on location of the patch
    void computePoynting() override;
    
};

#endif
