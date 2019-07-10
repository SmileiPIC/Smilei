
#include "LinearizedDomainDecomposition.h"


LinearizedDomainDecomposition::LinearizedDomainDecomposition( Params &params )
    : DomainDecomposition( params )
{
    ndomain_ = params.number_of_patches;
    
}


LinearizedDomainDecomposition1D::LinearizedDomainDecomposition1D( Params &params )
    : LinearizedDomainDecomposition( params )
{
}


LinearizedDomainDecomposition1D::~LinearizedDomainDecomposition1D( )
{
}


// generalhilbertindex
unsigned int LinearizedDomainDecomposition1D::getDomainId( std::vector<int> Coordinates )
{
    if( Coordinates[0] < 0 ) {
        return MPI_PROC_NULL;
    } else if( Coordinates[0] >= ( int )ndomain_[0] ) {
        return MPI_PROC_NULL;
    } else {
        return Coordinates[0];
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> LinearizedDomainDecomposition1D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 1, 0 );
    coords[0] = Id;
    return coords;
    
}


LinearizedDomainDecomposition2D::LinearizedDomainDecomposition2D( Params &params )
    : LinearizedDomainDecomposition( params )
{
}


LinearizedDomainDecomposition2D::~LinearizedDomainDecomposition2D( )
{
}


// generalhilbertindex
unsigned int LinearizedDomainDecomposition2D::getDomainId( std::vector<int> Coordinates )
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
        return ( Coordinates[0]*ndomain_[1] + Coordinates[1] );
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> LinearizedDomainDecomposition2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    coords[0] = Id/ndomain_[1];
    coords[1] = Id%ndomain_[1];
    return coords;
    
}


LinearizedDomainDecomposition2D_YX::LinearizedDomainDecomposition2D_YX( Params &params )
    : LinearizedDomainDecomposition( params )
{
}


LinearizedDomainDecomposition2D_YX::~LinearizedDomainDecomposition2D_YX( )
{
}


// generalhilbertindex
unsigned int LinearizedDomainDecomposition2D_YX::getDomainId( std::vector<int> Coordinates )
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
        return ( Coordinates[0] + Coordinates[1]*ndomain_[0] );
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> LinearizedDomainDecomposition2D_YX::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    coords[0] = Id%ndomain_[0];
    coords[1] = Id/ndomain_[0];
    return coords;
    
}



LinearizedDomainDecomposition3D::LinearizedDomainDecomposition3D( Params &params )
    : LinearizedDomainDecomposition( params )
{
}


LinearizedDomainDecomposition3D::~LinearizedDomainDecomposition3D( )
{
}


// generalhilbertindex
unsigned int LinearizedDomainDecomposition3D::getDomainId( std::vector<int> Coordinates )
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
        return ( Coordinates[0]*ndomain_[1]*ndomain_[2] + Coordinates[1]*ndomain_[2] + Coordinates[2] );
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> LinearizedDomainDecomposition3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    coords[0] = Id/( ndomain_[1]*ndomain_[2] );
    coords[1] = ( Id%( ndomain_[1]*ndomain_[2] ) )/ndomain_[2];
    coords[2] = ( Id%( ndomain_[1]*ndomain_[2] ) )%ndomain_[2];
    return coords;
    
}


LinearizedDomainDecomposition3D_ZYX::LinearizedDomainDecomposition3D_ZYX( Params &params )
    : LinearizedDomainDecomposition( params )
{
}


LinearizedDomainDecomposition3D_ZYX::~LinearizedDomainDecomposition3D_ZYX( )
{
}


// generalhilbertindex
unsigned int LinearizedDomainDecomposition3D_ZYX::getDomainId( std::vector<int> Coordinates )
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
        return ( Coordinates[0] + Coordinates[1]*ndomain_[0] + Coordinates[2]*ndomain_[0]*ndomain_[1] );
    }
    
}


// generalhilbertindexinv
std::vector<unsigned int> LinearizedDomainDecomposition3D_ZYX::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    coords[0] = ( Id%( ndomain_[0]*ndomain_[1] ) )%ndomain_[0];
    coords[1] = ( Id%( ndomain_[0]*ndomain_[1] ) )/ndomain_[0];
    coords[2] = Id/( ndomain_[0]*ndomain_[1] );
    return coords;
    
}
