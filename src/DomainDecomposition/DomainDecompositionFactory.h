#ifndef DOMAINDECOMPOSITIONFACTORY_H
#define DOMAINDECOMPOSITIONFACTORY_H

#include "DomainDecomposition.h"
// Patches decomposition along tht Hilbert curve
#include "HilbertDomainDecomposition.h"
// Patches decomposition along a linearized curve
#include "LinearizedDomainDecomposition.h"
// Domain decomposition (linearized)
#include "RegionDomainDecomposition.h"

class DomainDecompositionFactory
{
public:
    static DomainDecomposition *create( Params &params )
    {
        DomainDecomposition *domain_decomposition = NULL;
        
        TITLE( "Patch arrangement : " );
        
        
        if( params.patch_arrangement=="hilbertian" ) {
            if( ( params.geometry == "1Dcartesian" ) ) {
                domain_decomposition = new HilbertDomainDecomposition1D( params );
            } else if( ( params.geometry == "2Dcartesian" ) || ( params.geometry == "AMcylindrical" ) ) {
                domain_decomposition = new HilbertDomainDecomposition2D( params );
            } else if( ( params.geometry == "3Dcartesian" ) ) {
                domain_decomposition = new HilbertDomainDecomposition3D( params );
            } else {
                ERROR( "Unknown geometry" );
            }
        } else {
        
            bool enable_diagField( true );
            
            int nmpi = 1;
            MPI_Comm_size( MPI_COMM_WORLD, &nmpi );
            
            int npatches( 1 );
            for( unsigned int iDim=0 ; iDim<params.nDim_field ; iDim++ ) {
                npatches *= params.number_of_patches[iDim];
            }
            int npatches_per_rank = npatches / nmpi;
            
            if( ( npatches%nmpi!=0 )||( params.patch_arrangement == "linearized_YX" )||( params.patch_arrangement == "linearized_ZYX" ) ) {
                enable_diagField = false;
            } else {
                if( params.patch_arrangement == "linearized_XY" ) {
                
                    if( ( params.number_of_patches[1]%npatches_per_rank!=0 )
                            && ( npatches_per_rank%params.number_of_patches[1]!=0 ) ) {
                        enable_diagField = false;
                    }
                    
                } else if( params.patch_arrangement == "linearized_XYZ" ) {
                
                    if( ( params.number_of_patches[2]%npatches_per_rank!=0 )
                            && ( npatches_per_rank%params.number_of_patches[2]!=0 )
                            && ( npatches_per_rank%( params.number_of_patches[1]*params.number_of_patches[2]!=0 ) ) ) {
                        enable_diagField = false;
                    }
                } else {
                    ERROR( params.patch_arrangement << " not supported" );
                }
            }
            
            if( !enable_diagField ) {
                WARNING( "DiagFields not reliable because of the patch arrangement !!!" );
            }
            
            if( ( params.geometry == "1Dcartesian" ) ) {
                domain_decomposition = new LinearizedDomainDecomposition1D( params );
            } else if( ( params.geometry == "2Dcartesian" )  || ( params.geometry == "AMcylindrical" ) ) {
                if( params.patch_arrangement=="linearized_XY" ) {
                    domain_decomposition = new LinearizedDomainDecomposition2D( params );
                } else if( params.patch_arrangement=="linearized_YX" ) {
                    domain_decomposition = new LinearizedDomainDecomposition2D_YX( params );
                }
            } else if( ( params.geometry == "3Dcartesian" ) ) {
                if( params.patch_arrangement=="linearized_XYZ" ) {
                    domain_decomposition = new LinearizedDomainDecomposition3D( params );
                } else if( params.patch_arrangement=="linearized_ZYX" ) {
                    domain_decomposition = new LinearizedDomainDecomposition3D_ZYX( params );
                }
            } else {
                ERROR( "Unknown geometry" );
            }
        }
        
        return domain_decomposition;
    }
    
    static DomainDecomposition *createGlobal( Params &params )
    {
        DomainDecomposition *domain_decomposition = NULL;
        
        if( ( params.geometry == "1Dcartesian" ) ) {
            domain_decomposition = new RegionDomainDecomposition1D( params );
        } else if( ( params.geometry == "2Dcartesian" ) || ( params.geometry == "AMcylindrical" ) ) {
            domain_decomposition = new RegionDomainDecomposition2D( params );
        } else if( ( params.geometry == "3Dcartesian" ) ) {
            domain_decomposition = new RegionDomainDecomposition3D( params );
        } else {
            ERROR( "Unknown geometry" );
        }
        
        return domain_decomposition;
    }
    
    
};

#endif
