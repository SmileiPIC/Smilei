#ifndef DOMAINDECOMPOSITION_H
#define DOMAINDECOMPOSITION_H

#include <vector>

#include "Params.h"

class DomainDecomposition
{
public:
    DomainDecomposition( Params & ) {};
    virtual ~DomainDecomposition( ) {};
    
    std::vector<unsigned int> ndomain_;
    
    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;
    
};

#endif

