#ifndef SOLVER_H
#define SOLVER_H

#include "Params.h"

class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver
//  --------------------------------------------------------------------------------------------------------------------
class Solver
{

public:
    //! Creator for Solver
    Solver(Params &params) {};
    virtual ~Solver() {};

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields) = 0;

protected:

};//END class

#endif

