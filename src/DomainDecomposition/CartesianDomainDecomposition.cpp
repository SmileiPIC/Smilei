#include "CartesianDomainDecomposition.h"


CartesianDomainDecomposition::CartesianDomainDecomposition( Params &params )
    : DomainDecomposition( params )
{
    ndomain_ = params.number_of_domain;
    block_size_.resize( ndomain_.size(), 2 );
}


CartesianDomainDecomposition1D::CartesianDomainDecomposition1D( Params &params )
    : CartesianDomainDecomposition( params )
{
}


CartesianDomainDecomposition1D::~CartesianDomainDecomposition1D( )
{
}


// generalhilbertindex
unsigned int CartesianDomainDecomposition1D::getDomainId( std::vector<int> Coordinates )
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
std::vector<unsigned int> CartesianDomainDecomposition1D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 1, 0 );
    coords[0] = Id;
    return coords;
    
}


CartesianDomainDecomposition2D::CartesianDomainDecomposition2D( Params &params )
    : CartesianDomainDecomposition( params )
{
}


CartesianDomainDecomposition2D::~CartesianDomainDecomposition2D( )
{
}


// generalhilbertindex
unsigned int CartesianDomainDecomposition2D::getDomainId( std::vector<int> Coordinates )
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
std::vector<unsigned int> CartesianDomainDecomposition2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    coords[0] = ( double )Id/( double )ndomain_[1];
    coords[1] = Id - coords[0]*ndomain_[1];
    return coords;
    
}


CartesianDomainDecomposition3D::CartesianDomainDecomposition3D( Params &params )
    : CartesianDomainDecomposition( params )
{
}


CartesianDomainDecomposition3D::~CartesianDomainDecomposition3D( )
{
}


// generalhilbertindex
unsigned int CartesianDomainDecomposition3D::getDomainId( std::vector<int> Coordinates )
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
std::vector<unsigned int> CartesianDomainDecomposition3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    coords[0] = ( double )Id/( double )ndomain_[1]/( double )ndomain_[2];
    coords[1] = ( double )( Id - coords[0]*ndomain_[1]*ndomain_[2] ) / ( double )ndomain_[2];
    coords[2] = Id - coords[0]*ndomain_[1]*ndomain_[2] - coords[1]*ndomain_[2];
    return coords;
    
}
