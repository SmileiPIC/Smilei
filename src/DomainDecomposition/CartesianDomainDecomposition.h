#ifndef CARTESIANDOMAINDECOMPOSITION_H
#define CARTESIANDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"

class CartesianDomainDecomposition : public DomainDecomposition
{
public:
    CartesianDomainDecomposition( Params &params );
    virtual ~CartesianDomainDecomposition( ) {};
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
protected:
    std::vector<unsigned int> mi_;
    
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

class CartesianDomainDecomposition2D_YX : public CartesianDomainDecomposition
{
public:
    CartesianDomainDecomposition2D_YX( Params &params );
    ~CartesianDomainDecomposition2D_YX( ) override final;
    
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

class CartesianDomainDecomposition3D_ZYX : public CartesianDomainDecomposition
{
public:
    CartesianDomainDecomposition3D_ZYX( Params &params );
    ~CartesianDomainDecomposition3D_ZYX( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
