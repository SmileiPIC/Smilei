#ifndef CARTESIANGEOMETRY_H
#define CARTESIANGEOMETRY_H

#include "Geometry.h"

class CartesianGeometry : public Geometry
{
public:
    CartesianGeometry( Params& params );
    virtual ~CartesianGeometry( ) {};

    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;

protected:
    std::vector<unsigned int> block_size_;

};


class CartesianGeometry1D : public CartesianGeometry
{
public:
    CartesianGeometry1D( Params& params );
    ~CartesianGeometry1D( ) override final;

    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class CartesianGeometry2D : public CartesianGeometry
{
public:
    CartesianGeometry2D( Params& params );
    ~CartesianGeometry2D( ) override final;

    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class CartesianGeometry3D : public CartesianGeometry
{
public:
    CartesianGeometry3D( Params& params );
    ~CartesianGeometry3D( ) override final;

    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
