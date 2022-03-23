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
    Solver( Params &params ) {};
    virtual ~Solver() {};
    
    virtual void coupling( Params &params, ElectroMagn *EMfields, bool full_domain = false ) {};
    virtual void uncoupling() {};
    virtual void rotational_cleaning( ElectroMagn *EMfields ) {};
    virtual void densities_correction( ElectroMagn *EMfields ) {};
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;

    virtual void setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch ) {ERROR("Not using PML");};
    virtual void compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, int solver_min, int solver_max ) {ERROR("Not using PML");};
    virtual void compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, int solver_min, int solver_max ) {ERROR("Not using PML");};
    virtual void compute_A_from_G( LaserEnvelope *envelope, int iDim, int min_or_max, int solver_min, int solver_max ) {ERROR("Not using PML");};

protected:

};//END class

class NullSolver : public Solver
{

public:
    NullSolver( Params &params ) : Solver( params ) {};
    virtual ~NullSolver() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) {};
    
protected:

};//END class


#endif

