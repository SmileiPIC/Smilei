#ifndef NEWDOMAINDECOMPOSITION_H
#define NEWDOMAINDECOMPOSITION_H

#include "DomainDecomposition.h"

class NewDomainDecomposition : public DomainDecomposition
{
public:
    NewDomainDecomposition( Params& params );
    virtual ~NewDomainDecomposition( ) {};

    virtual unsigned int getDomainId( std::vector<int> Coordinates ) = 0;
    virtual std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) = 0;

protected:
    std::vector<unsigned int> mi_;

};


class NewDomainDecomposition1D : public NewDomainDecomposition
{
public:
    NewDomainDecomposition1D( Params& params );
    ~NewDomainDecomposition1D( ) override final;

    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class NewDomainDecomposition2D : public NewDomainDecomposition
{
public:
    NewDomainDecomposition2D( Params& params );
    ~NewDomainDecomposition2D( ) override final;

    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};


class NewDomainDecomposition3D : public NewDomainDecomposition
{
public:
    NewDomainDecomposition3D( Params& params );
    ~NewDomainDecomposition3D( ) override final;

    unsigned int getDomainId( std::vector<int> Coordinates ) override final;
    std::vector<unsigned int> getDomainCoordinates( unsigned int Id ) override final;
};

#endif
