#ifndef DOMAINDECOMPOSITIONFACTORY_H
#define DOMAINDECOMPOSITIONFACTORY_H

#include "DomainDecomposition.h"
#include "HilbertDomainDecomposition.h"
#include "CartesianDomainDecomposition.h"
#include "NewDomainDecomposition.h"

class DomainDecompositionFactory {
public:
    static DomainDecomposition* create(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if ( ( params.geometry == "1Dcartesian" ) )
            domain_decomposition = new HilbertDomainDecomposition1D( params );
        else if ( ( params.geometry == "2Dcartesian" ) ) 
            domain_decomposition = new HilbertDomainDecomposition2D( params );
        else if ( ( params.geometry == "3Dcartesian" ) ) 
            domain_decomposition = new HilbertDomainDecomposition3D( params );
        else
            ERROR( "Unknown geometry" );

        return domain_decomposition;
    }

    static DomainDecomposition* createNew(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if ( ( params.geometry == "1Dcartesian" ) )
            domain_decomposition = new NewDomainDecomposition1D( params );
        else if ( ( params.geometry == "2Dcartesian" ) ) 
            domain_decomposition = new NewDomainDecomposition2D( params );
        else if ( ( params.geometry == "3Dcartesian" ) ) 
            domain_decomposition = new NewDomainDecomposition3D( params );
        else
            ERROR( "Unknown geometry" );

        return domain_decomposition;
    }

    static DomainDecomposition* createGlobal(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if ( ( params.geometry == "1Dcartesian" ) )
            domain_decomposition = new CartesianDomainDecomposition1D( params );
        else if ( ( params.geometry == "2Dcartesian" ) ) 
            domain_decomposition = new CartesianDomainDecomposition2D( params );
        else if ( ( params.geometry == "3Dcartesian" ) ) 
            domain_decomposition = new CartesianDomainDecomposition3D( params );
        else
            ERROR( "Unknown geometry" );

        return domain_decomposition;
    }


};

#endif
