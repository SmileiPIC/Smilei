#ifndef GLOBALDOMAINDECOMPOSITION_H
#define GLOBALDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"

class GlobalDomainDecomposition : public DomainDecomposition
{
public:
    GlobalDomainDecomposition( Params &params );
    virtual ~GlobalDomainDecomposition( ) {};
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
protected:
    std::vector<unsigned int> block_size_;
    
};


class GlobalDomainDecomposition1D : public GlobalDomainDecomposition
{
public:
    GlobalDomainDecomposition1D( Params &params );
    ~GlobalDomainDecomposition1D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class GlobalDomainDecomposition2D : public GlobalDomainDecomposition
{
public:
    GlobalDomainDecomposition2D( Params &params );
    ~GlobalDomainDecomposition2D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class GlobalDomainDecomposition3D : public GlobalDomainDecomposition
{
public:
    GlobalDomainDecomposition3D( Params &params );
    ~GlobalDomainDecomposition3D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
