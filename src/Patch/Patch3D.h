#ifndef PATCH3D_H
#define PATCH3D_H

#include "Patch.h"
#include "Field3D.h"
#include "cField3D.h"

class SimWindow;

//! Class Patch : sub MPI domain
//!     Collection of patch = MPI domain
class Patch3D final : public Patch
{
public:
    //! Constructor for Patch
    Patch3D( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved );
    //! Cloning Constructor for Patch
    Patch3D( Patch3D *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles );
    
    void initStep2( Params &params, DomainDecomposition *domain_decomposition ) override final;
    
    //! Destructor for Patch
    ~Patch3D() override  final;
    
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
        if( p->position(1,ipart) - getDomainLocalMin(1) < halfcell
         || getDomainLocalMax(1) - p->position(1,ipart) < halfcell ) {
             factor *= 0.5;
        }
        
        halfcell = 0.5 * params.cell_length[2];
        if( p->position(2,ipart) - getDomainLocalMin(2) < halfcell
         || getDomainLocalMax(2) - p->position(2,ipart) < halfcell ) {
             factor *= 0.5;
        }
        
        return factor * cell_volume;
    };
    
    //! Given several arrays (x,y,z), return indices of points in patch
    std::vector<unsigned int> indicesInDomain( double **position, unsigned int n_particles ) override
    {
        std::vector<unsigned int> indices( 0 );
        for( unsigned int ip = 0; ip < n_particles; ip++ ) {
            if( position[0][ip] >= getDomainLocalMin( 0 ) && position[0][ip] < getDomainLocalMax( 0 )
             && position[1][ip] >= getDomainLocalMin( 1 ) && position[1][ip] < getDomainLocalMax( 1 )
             && position[2][ip] >= getDomainLocalMin( 2 ) && position[2][ip] < getDomainLocalMax( 2 ) ) {
                indices.push_back( ip );
            }
        }
        return indices;
    };
    
    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------
    
    void exchangeField_movewin( Field* field, int clrw ) override final;
    
    // Create MPI_Datatype to exchange fields
    void createType2( Params &params ) override final;
    void cleanType() override final;
    
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntype_[2][2][2];
    
    
    
};

#endif
