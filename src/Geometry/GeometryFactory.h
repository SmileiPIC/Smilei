#ifndef GEOMETRYFACTORY_H
#define GEOMETRYFACTORY_H

#include "Geometry.h"
#include "HilbertGeometry.h"

class GeometryFactory {
public:
    static Geometry* create(Params& params) {
        Geometry* geometry = NULL;

        if ( ( params.geometry == "1d3v" ) )
            geometry = new HilbertGeometry1D( params );
        else if ( ( params.geometry == "2d3v" ) ) 
            geometry = new HilbertGeometry2D( params );
        else if ( ( params.geometry == "3d3v" ) ) 
            geometry = new HilbertGeometry3D( params );

        return geometry;
    }
};

#endif
