#ifndef HILBERTGEOMETRY_H
#define HILBERTGEOMETRY_H

#include "Geometry.h"

class HilbertGeometry : public Geometry
{
public:
    HilbertGeometry( Params& params );
    virtual ~HilbertGeometry( ) {};

    virtual unsigned int getDomainId( std::vector<unsigned int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;

protected:
    std::vector<unsigned int> mi_;

};


class HilbertGeometry2D : public HilbertGeometry
{
    HilbertGeometry2D( Params& params );
    ~HilbertGeometry2D( ) override final;

    unsigned int getDomainId( std::vector<unsigned int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class HilbertGeometry3D : public HilbertGeometry
{
    HilbertGeometry3D( Params& params );
    ~HilbertGeometry3D( ) override final;

    unsigned int getDomainId( std::vector<unsigned int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
