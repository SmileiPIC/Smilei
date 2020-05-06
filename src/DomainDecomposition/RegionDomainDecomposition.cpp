#include "RegionDomainDecomposition.h"


RegionDomainDecomposition::RegionDomainDecomposition( Params &params )
    : DomainDecomposition( params )
{
    ndomain_ = params.number_of_region;
    block_size_.resize( ndomain_.size(), 2 );
}


RegionDomainDecomposition1D::RegionDomainDecomposition1D( Params &params )
    : RegionDomainDecomposition( params )
{
}


RegionDomainDecomposition1D::~RegionDomainDecomposition1D( )
{
}


// generalhilbertindex
unsigned int RegionDomainDecomposition1D::getDomainId( std::vector<int> Coordinates )
{
    if( Coordinates[0] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[0] >= ( int )ndomain_[0] )
    
    {
        return MPI_PROC_NULL;
    } else {
        return Coordinates[0];
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> RegionDomainDecomposition1D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 1, 0 );
    coords[0] = Id;
    return coords;
    
}


RegionDomainDecomposition2D::RegionDomainDecomposition2D( Params &params )
    : RegionDomainDecomposition( params )
{
}


RegionDomainDecomposition2D::~RegionDomainDecomposition2D( )
{
}


// generalhilbertindex
unsigned int RegionDomainDecomposition2D::getDomainId( std::vector<int> Coordinates )
{
    if( Coordinates[0] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[0] >= ( int )ndomain_[0] ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[1] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[1] >= ( int )ndomain_[1] ) {
        return MPI_PROC_NULL;
    } else {
        return Coordinates[0]*ndomain_[1]+Coordinates[1];
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> RegionDomainDecomposition2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    coords[0] = ( double )Id/( double )ndomain_[1];
    coords[1] = Id - coords[0]*ndomain_[1];
    return coords;
    
}


RegionDomainDecomposition3D::RegionDomainDecomposition3D( Params &params )
    : RegionDomainDecomposition( params )
{
}


RegionDomainDecomposition3D::~RegionDomainDecomposition3D( )
{
}


// generalhilbertindex
unsigned int RegionDomainDecomposition3D::getDomainId( std::vector<int> Coordinates )
{
    if( Coordinates[0] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[0] >= ( int )ndomain_[0] ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[1] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[1] >= ( int )ndomain_[1] ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[2] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[2] >= ( int )ndomain_[2] ) {
        return MPI_PROC_NULL;
    } else {
        return Coordinates[0]*ndomain_[1]*ndomain_[2]+Coordinates[1]*ndomain_[2]+Coordinates[2];
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> RegionDomainDecomposition3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    coords[0] = ( double )Id/( double )ndomain_[1]/( double )ndomain_[2];
    coords[1] = ( double )( Id - coords[0]*ndomain_[1]*ndomain_[2] ) / ( double )ndomain_[2];
    coords[2] = Id - coords[0]*ndomain_[1]*ndomain_[2] - coords[1]*ndomain_[2];
    return coords;
    
}
