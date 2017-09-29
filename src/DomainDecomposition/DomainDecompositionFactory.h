#ifndef DOMAINDECOMPOSITIONFACTORY_H
#define DOMAINDECOMPOSITIONFACTORY_H

#include "DomainDecomposition.h"
#include "HilbertDomainDecomposition.h"
#include "CartesianDomainDecomposition.h"

class DomainDecompositionFactory {
public:
    static DomainDecomposition* create(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if ( ( params.geometry == "1d3v" ) )
            domain_decomposition = new HilbertDomainDecomposition1D( params );
        else if ( ( params.geometry == "2d3v" ) ) 
            domain_decomposition = new HilbertDomainDecomposition2D( params );
        else if ( ( params.geometry == "3d3v" ) ) 
            domain_decomposition = new HilbertDomainDecomposition3D( params );

        return domain_decomposition;
    }

    static DomainDecomposition* createGlobal(Params& params) {
        DomainDecomposition* domain_decomposition = NULL;

        if ( ( params.geometry == "1d3v" ) )
            domain_decomposition = new CartesianDomainDecomposition1D( params );
        else if ( ( params.geometry == "2d3v" ) ) 
            domain_decomposition = new CartesianDomainDecomposition2D( params );
        else if ( ( params.geometry == "3d3v" ) ) 
            domain_decomposition = new CartesianDomainDecomposition3D( params );

        return domain_decomposition;
    }


};

#endif
