#ifndef GLOBALDOMAINDECOMPOSITION_H
#define GLOBALDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"

class CartesianDomainDecomposition : public DomainDecomposition
{
public:
    CartesianDomainDecomposition( Params &params );
    virtual ~CartesianDomainDecomposition( ) {};
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
protected:
    std::vector<unsigned int> block_size_;
    
};


class CartesianDomainDecomposition1D : public CartesianDomainDecomposition
{
public:
    CartesianDomainDecomposition1D( Params &params );
    ~CartesianDomainDecomposition1D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class CartesianDomainDecomposition2D : public CartesianDomainDecomposition
{
public:
    CartesianDomainDecomposition2D( Params &params );
    ~CartesianDomainDecomposition2D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class CartesianDomainDecomposition3D : public CartesianDomainDecomposition
{
public:
    CartesianDomainDecomposition3D( Params &params );
    ~CartesianDomainDecomposition3D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
