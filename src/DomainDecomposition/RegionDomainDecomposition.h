#ifndef REGIONDOMAINDECOMPOSITION_H
#define REGIONDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"

class RegionDomainDecomposition : public DomainDecomposition
{
public:
    RegionDomainDecomposition( Params &params );
    virtual ~RegionDomainDecomposition( ) {};
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
protected:
    std::vector<unsigned int> block_size_;
    
};


class RegionDomainDecomposition1D : public RegionDomainDecomposition
{
public:
    RegionDomainDecomposition1D( Params &params );
    ~RegionDomainDecomposition1D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class RegionDomainDecomposition2D : public RegionDomainDecomposition
{
public:
    RegionDomainDecomposition2D( Params &params );
    ~RegionDomainDecomposition2D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class RegionDomainDecomposition3D : public RegionDomainDecomposition
{
public:
    RegionDomainDecomposition3D( Params &params );
    ~RegionDomainDecomposition3D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
