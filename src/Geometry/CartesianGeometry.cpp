#include "CartesianGeometry.h"


CartesianGeometry::CartesianGeometry( Params& params )
    : Geometry( params )
{
    ndomain_ = params.number_of_patches;
    block_size_.resize( ndomain_.size(), 2 );
}


CartesianGeometry1D::CartesianGeometry1D( Params& params )
    : CartesianGeometry( params )
{
}


CartesianGeometry1D::~CartesianGeometry1D( )
{
}


// generalhilbertindex
unsigned int CartesianGeometry1D::getDomainId( std::vector<unsigned int> Coordinates )
{
    if ( Coordinates[0] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[0] >= ndomain_[0] )
        return MPI_PROC_NULL;
    else
        return Coordinates[0];

}


// generalhilbertindexinv
std::vector<unsigned int> CartesianGeometry1D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 1, 0 );
    coords[0] = Id;
    return coords;

}


CartesianGeometry2D::CartesianGeometry2D( Params& params )
    : CartesianGeometry( params )
{
}


CartesianGeometry2D::~CartesianGeometry2D( )
{
}


// generalhilbertindex
unsigned int CartesianGeometry2D::getDomainId( std::vector<unsigned int> Coordinates )
{
    if ( Coordinates[0] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[0] >= ndomain_[0] )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] >= ndomain_[1] )
        return MPI_PROC_NULL;
    else
        return Coordinates[0]*ndomain_[1]+Coordinates[1];

}


// generalhilbertindexinv
std::vector<unsigned int> CartesianGeometry2D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 2, 0 );
    coords[0] = (double)Id/(double)ndomain_[1];
    coords[1] = Id - coords[0]*ndomain_[1];
    return coords;

}


CartesianGeometry3D::CartesianGeometry3D( Params& params )
    : CartesianGeometry( params )
{
}


CartesianGeometry3D::~CartesianGeometry3D( )
{
}


// generalhilbertindex
unsigned int CartesianGeometry3D::getDomainId( std::vector<unsigned int> Coordinates )
{
    if ( Coordinates[0] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[0] >= ndomain_[0] )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[1] >= ndomain_[1] )
        return MPI_PROC_NULL;
    else if ( Coordinates[2] < 0 )
        return MPI_PROC_NULL;
    else if ( Coordinates[2] >= ndomain_[2] )
        return MPI_PROC_NULL;
    else
        return Coordinates[0]*ndomain_[1]*ndomain_[2]+Coordinates[1]*ndomain_[2]+Coordinates[2];

}


// generalhilbertindexinv
std::vector<unsigned int> CartesianGeometry3D::getDomainCoordinates( unsigned int Id )
{
    std::vector<unsigned int> coords( 3, 0 );
    coords[0] = (double)Id/(double)ndomain_[1]/(double)ndomain_[2];
    coords[1] = (double)(Id - coords[0]*ndomain_[1]*ndomain_[2]) / (double)ndomain_[2];
    coords[2] = Id - coords[0]*ndomain_[1]*ndomain_[2] - coords[1]*ndomain_[2];
    return coords;

}
