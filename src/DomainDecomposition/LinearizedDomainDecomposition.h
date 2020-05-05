#ifndef CARTESIANDOMAINDECOMPOSITION_H
#define CARTESIANDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"

class LinearizedDomainDecomposition : public DomainDecomposition
{
public:
    LinearizedDomainDecomposition( Params &params );
    virtual ~LinearizedDomainDecomposition( ) {};
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
protected:
    std::vector<unsigned int> mi_;
    
};


class LinearizedDomainDecomposition1D : public LinearizedDomainDecomposition
{
public:
    LinearizedDomainDecomposition1D( Params &params );
    ~LinearizedDomainDecomposition1D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class LinearizedDomainDecomposition2D : public LinearizedDomainDecomposition
{
public:
    LinearizedDomainDecomposition2D( Params &params );
    ~LinearizedDomainDecomposition2D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

class LinearizedDomainDecomposition2D_YX : public LinearizedDomainDecomposition
{
public:
    LinearizedDomainDecomposition2D_YX( Params &params );
    ~LinearizedDomainDecomposition2D_YX( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class LinearizedDomainDecomposition3D : public LinearizedDomainDecomposition
{
public:
    LinearizedDomainDecomposition3D( Params &params );
    ~LinearizedDomainDecomposition3D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

class LinearizedDomainDecomposition3D_ZYX : public LinearizedDomainDecomposition
{
public:
    LinearizedDomainDecomposition3D_ZYX( Params &params );
    ~LinearizedDomainDecomposition3D_ZYX( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
