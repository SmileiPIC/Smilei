
#include "HilbertGeometry.h"
#include "Hilbert_functions.h"


HilbertGeometry::HilbertGeometry( Params& params )
    : Geometry( params )
{
    mi_.resize( params.nDim_field );

}


HilbertGeometry2D::HilbertGeometry2D( Params& params )
    : HilbertGeometry( params )
{
}


HilbertGeometry2D::~HilbertGeometry2D( )
{
}


// generalhilbertindex
unsigned int HilbertGeometry2D::getDomainId( std::vector<unsigned int> Coordinates )
{
    return generalhilbertindex( mi_[0], mi_[1], Coordinates[0], Coordinates[1] );

}


// generalhilbertindexinv
std::vector<unsigned int> HilbertGeometry2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    generalhilbertindexinv( mi_[0], mi_[1], &coords[0], &coords[1], Id );
    return coords;

}


HilbertGeometry3D::HilbertGeometry3D( Params& params )
    : HilbertGeometry( params )
{
}


HilbertGeometry3D::~HilbertGeometry3D( )
{
}


// generalhilbertindex
unsigned int HilbertGeometry3D::getDomainId( std::vector<unsigned int> Coordinates )
{
    return generalhilbertindex( mi_[0], mi_[1], mi_[2], Coordinates[0], Coordinates[1], Coordinates[2] );

}


// generalhilbertindexinv
std::vector<unsigned int> HilbertGeometry3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    generalhilbertindexinv( mi_[0], mi_[1], mi_[2], &coords[0], &coords[1], &coords[2], Id );
    return coords;

}
