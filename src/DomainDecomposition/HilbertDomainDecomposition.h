#ifndef HILBERTDOMAINDECOMPOSITION_H
#define HILBERTDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"
#include "Hilbert_functions.h"

class HilbertDomainDecomposition : public DomainDecomposition
{
public:
    HilbertDomainDecomposition( Params &params );
    virtual ~HilbertDomainDecomposition( ) {};
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
protected:
    std::vector<unsigned int> mi_;
    
};


class HilbertDomainDecomposition1D final : public HilbertDomainDecomposition
{
public:
    HilbertDomainDecomposition1D( Params &params );
    ~HilbertDomainDecomposition1D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class HilbertDomainDecomposition2D final : public HilbertDomainDecomposition
{
public:
    HilbertDomainDecomposition2D( Params &params );
    ~HilbertDomainDecomposition2D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class HilbertDomainDecomposition3D final : public HilbertDomainDecomposition
{
public:
    HilbertDomainDecomposition3D( Params &params );
    ~HilbertDomainDecomposition3D( ) override final;
    
    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
