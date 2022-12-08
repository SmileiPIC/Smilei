#ifndef SOLVER_H
#define SOLVER_H

#include "Params.h"

class ElectroMagn;
class LaserEnvelope;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver
//  --------------------------------------------------------------------------------------------------------------------
class Solver
{

public:
    //! Creator for Solver
    Solver() {};
    virtual ~Solver() {};
    
    virtual void coupling( Params &, ElectroMagn *, bool = false ) {};
    virtual void uncoupling() {};
    virtual void rotational_cleaning( ElectroMagn * ) {};
    virtual void densities_correction( ElectroMagn * ) {};
    //! Overloading of () operator
    virtual void operator()( ElectroMagn * ) = 0;

    virtual void setDomainSizeAndCoefficients( int, int, std::vector<unsigned int>, int, int, int*, int*, Patch* ) {ERROR("Not using PML");};
    virtual void compute_E_from_D( ElectroMagn *, int, int, std::vector<unsigned int>, unsigned int, unsigned int ) {ERROR("Not using PML");};
    virtual void compute_H_from_B( ElectroMagn *, int, int, std::vector<unsigned int>, unsigned int, unsigned int ) {ERROR("Not using PML");};
    virtual void compute_A_from_G( LaserEnvelope *, int, int, std::vector<unsigned int>, unsigned int, unsigned int ) {ERROR("Not using PML");};

protected:

};//END class

class NullSolver : public Solver
{

public:
    NullSolver() : Solver() {};
    virtual ~NullSolver() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn * ) {};
    
protected:

};//END class


#endif

