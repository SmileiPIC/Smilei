#ifndef GEOMETRYFACTORY_H
#define GEOMETRYFACTORY_H

#include "DomainDecomposition.h"
#include "HilbertDomainDecomposition.h"
#include "CartesianDomainDecomposition.h"

class DomainDecompositionFactory {
public:
    static DomainDecomposition* create(Params& params) {
        DomainDecomposition* geometry = NULL;

        if ( ( params.geometry == "1d3v" ) )
            geometry = new HilbertDomainDecomposition1D( params );
        else if ( ( params.geometry == "2d3v" ) ) 
            geometry = new HilbertDomainDecomposition2D( params );
        else if ( ( params.geometry == "3d3v" ) ) 
            geometry = new HilbertDomainDecomposition3D( params );

        return geometry;
    }

    static DomainDecomposition* createGlobal(Params& params) {
        DomainDecomposition* geometry = NULL;

        if ( ( params.geometry == "1d3v" ) )
            geometry = new CartesianDomainDecomposition1D( params );
        else if ( ( params.geometry == "2d3v" ) ) 
            geometry = new CartesianDomainDecomposition2D( params );
        else if ( ( params.geometry == "3d3v" ) ) 
            geometry = new CartesianDomainDecomposition3D( params );

        return geometry;
    }


};

#endif
