#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

#include "Params.h"

class Geometry
{
public:
    Geometry( Params& params ) {};
    virtual ~Geometry( ) {};

    std::vector<unsigned int> ndomain_;

    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;

};

#endif

