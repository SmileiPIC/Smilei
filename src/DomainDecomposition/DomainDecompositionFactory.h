#ifndef DOMAINDECOMPOSITIONFACTORY_H
#define DOMAINDECOMPOSITIONFACTORY_H

#include "DomainDecomposition.h"
#include "HilbertDomainDecomposition.h"
#include "GlobalDomainDecomposition.h"
#include "CartesianDomainDecomposition.h"

class DomainDecompositionFactory {
public:
    static DomainDecomposition* create(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if (params.patch_decomposition=="hilbert") {
            if ( ( params.geometry == "1Dcartesian" ) )
                domain_decomposition = new HilbertDomainDecomposition1D( params );
            else if ( ( params.geometry == "2Dcartesian" ) || ( params.geometry == "3drz" ) ) 
                domain_decomposition = new HilbertDomainDecomposition2D( params );
            else if ( ( params.geometry == "3Dcartesian" ) ) 
                domain_decomposition = new HilbertDomainDecomposition3D( params );
            else
                ERROR( "Unknown geometry" );
        }
        else if(params.patch_decomposition=="cartesian") {
            if ( ( params.geometry == "1Dcartesian" ) )
                domain_decomposition = new CartesianDomainDecomposition1D( params );
            else if ( ( params.geometry == "2Dcartesian" )  || ( params.geometry == "3drz" )  ) {
                if (params.patch_orientation!="YX")
                    domain_decomposition = new CartesianDomainDecomposition2D( params );
                else 
                    domain_decomposition = new CartesianDomainDecomposition2D_YX( params );
            }
            else if ( ( params.geometry == "3Dcartesian" ) ) {
                if (params.patch_orientation!="ZYX")
                    domain_decomposition = new CartesianDomainDecomposition3D( params );
                else
                    domain_decomposition = new CartesianDomainDecomposition3D_ZYX( params );
            }
            else
                ERROR( "Unknown geometry" );
        }
        else
            ERROR( "Unknown geometry" );

        return domain_decomposition;
    }

    static DomainDecomposition* createGlobal(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if ( ( params.geometry == "1Dcartesian" ) )
            domain_decomposition = new GlobalDomainDecomposition1D( params );
        else if ( ( params.geometry == "2Dcartesian" ) || ( params.geometry == "3drz" ) ) 
            domain_decomposition = new GlobalDomainDecomposition2D( params );
        else if ( ( params.geometry == "3Dcartesian" ) ) 
            domain_decomposition = new GlobalDomainDecomposition3D( params );
        else
            ERROR( "Unknown geometry" );

        return domain_decomposition;
    }


};

#endif
