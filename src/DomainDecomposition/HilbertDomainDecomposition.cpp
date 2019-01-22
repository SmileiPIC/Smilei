
#include "HilbertDomainDecomposition.h"


HilbertDomainDecomposition::HilbertDomainDecomposition( Params &params )
    : DomainDecomposition( params )
{
    ndomain_ = params.number_of_patches;
    mi_ = params.mi;
    
}


HilbertDomainDecomposition1D::HilbertDomainDecomposition1D( Params &params )
    : HilbertDomainDecomposition( params )
{
}


HilbertDomainDecomposition1D::~HilbertDomainDecomposition1D( )
{
}


// generalhilbertindex
unsigned int HilbertDomainDecomposition1D::getDomainId( std::vector<int> Coordinates )
{
    return generalhilbertindex( mi_[0], 0, Coordinates[0], 0 );
    
}


// generalhilbertindexinv
std::vector<unsigned int> HilbertDomainDecomposition1D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 1, 0 );
    coords[0] = Id;
    return coords;
    
}


HilbertDomainDecomposition2D::HilbertDomainDecomposition2D( Params &params )
    : HilbertDomainDecomposition( params )
{
}


HilbertDomainDecomposition2D::~HilbertDomainDecomposition2D( )
{
}


// generalhilbertindex
unsigned int HilbertDomainDecomposition2D::getDomainId( std::vector<int> Coordinates )
{
    return generalhilbertindex( mi_[0], mi_[1], Coordinates[0], Coordinates[1] );
    
}


// generalhilbertindexinv
std::vector<unsigned int> HilbertDomainDecomposition2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    generalhilbertindexinv( mi_[0], mi_[1], &coords[0], &coords[1], Id );
    return coords;
    
}


HilbertDomainDecomposition3D::HilbertDomainDecomposition3D( Params &params )
    : HilbertDomainDecomposition( params )
{
}


HilbertDomainDecomposition3D::~HilbertDomainDecomposition3D( )
{
}


// generalhilbertindex
unsigned int HilbertDomainDecomposition3D::getDomainId( std::vector<int> Coordinates )
{
    return generalhilbertindex( mi_[0], mi_[1], mi_[2], Coordinates[0], Coordinates[1], Coordinates[2] );
    
}


// generalhilbertindexinv
std::vector<unsigned int> HilbertDomainDecomposition3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    generalhilbertindexinv( mi_[0], mi_[1], mi_[2], &coords[0], &coords[1], &coords[2], Id );
    return coords;
    
}
