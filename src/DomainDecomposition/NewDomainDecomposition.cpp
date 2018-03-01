
#include "NewDomainDecomposition.h"


NewDomainDecomposition::NewDomainDecomposition( Params& params )
    : DomainDecomposition( params )
{
    ndomain_ = params.number_of_patches;

}


NewDomainDecomposition1D::NewDomainDecomposition1D( Params& params )
    : NewDomainDecomposition( params )
{
}


NewDomainDecomposition1D::~NewDomainDecomposition1D( )
{
}


// generalhilbertindex
unsigned int NewDomainDecomposition1D::getDomainId( std::vector<int> Coordinates )
{
    if ( Coordinates[0] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[0] >= (int)ndomain_[0] )
        return MPI_PROC_NULL;
    else
        return Coordinates[0];

}


// generalhilbertindexinv
std::vector<unsigned int> NewDomainDecomposition1D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 1, 0 );
    coords[0] = Id;
    return coords;

}


NewDomainDecomposition2D::NewDomainDecomposition2D( Params& params )
    : NewDomainDecomposition( params )
{
}


NewDomainDecomposition2D::~NewDomainDecomposition2D( )
{
}


// generalhilbertindex
unsigned int NewDomainDecomposition2D::getDomainId( std::vector<int> Coordinates )
{
    if ( Coordinates[0] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[0] >= (int)ndomain_[0] )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] >= (int)ndomain_[1] )
        return MPI_PROC_NULL;
    else
        return ( Coordinates[0]*ndomain_[1] + Coordinates[1] );

}


// generalhilbertindexinv
std::vector<unsigned int> NewDomainDecomposition2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    coords[0] = Id/ndomain_[1];
    coords[1] = Id%ndomain_[1];
    return coords;

}


NewDomainDecomposition3D::NewDomainDecomposition3D( Params& params )
    : NewDomainDecomposition( params )
{
}


NewDomainDecomposition3D::~NewDomainDecomposition3D( )
{
}


// generalhilbertindex
unsigned int NewDomainDecomposition3D::getDomainId( std::vector<int> Coordinates )
{
    if ( Coordinates[0] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[0] >= (int)ndomain_[0] )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] >= (int)ndomain_[1] )
        return MPI_PROC_NULL;
    else if ( Coordinates[2] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[2] >= (int)ndomain_[2] )
        return MPI_PROC_NULL;
    else
        return ( Coordinates[0]*ndomain_[1]*ndomain_[2] + Coordinates[1]*ndomain_[2] + Coordinates[0] );

}


// generalhilbertindexinv
std::vector<unsigned int> NewDomainDecomposition3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    coords[0] = Id/(ndomain_[1]*ndomain_[2]);
    coords[1] = (Id%(ndomain_[1]*ndomain_[2]))/ndomain_[2];
    coords[2] = (Id%(ndomain_[1]*ndomain_[2]))%ndomain_[2];
    return coords;

}
