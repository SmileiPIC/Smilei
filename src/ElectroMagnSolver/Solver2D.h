#ifndef SOLVER2D_H
#define SOLVER2D_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver2D
//  --------------------------------------------------------------------------------------------------------------------
class Solver2D : public Solver
{

public:
    //! Creator for Solver
    Solver2D( Params &params ) : Solver()
    {
        dt = params.timestep;
        dx = params.cell_length[0];
        dy = params.cell_length[1];
        dt_ov_dx = params.timestep / params.cell_length[0];
        dt_ov_dy = params.timestep / params.cell_length[1];
    };
    virtual ~Solver2D() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    double dt;
    double dx;
    double dy;
    double dt_ov_dy;
    double dt_ov_dx;
    
};//END class

#endif

