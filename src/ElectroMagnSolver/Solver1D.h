#ifndef SOLVER1D_H
#define SOLVER1D_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver1D
//  --------------------------------------------------------------------------------------------------------------------
class Solver1D : public Solver
{

public:
    //! Creator for Solver
    Solver1D( Params &params ) : Solver()
    {
        dt = params.timestep;
        dt_ov_dx = params.timestep / params.cell_length[0];
    };
    virtual ~Solver1D() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    double dt;
    double dt_ov_dx;
    
};//END class

#endif

