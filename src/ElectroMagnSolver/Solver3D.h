#ifndef SOLVER3D_H
#define SOLVER3D_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Solver3D
//  --------------------------------------------------------------------------------------------------------------------
class Solver3D : public Solver
{

public:
    //! Creator for Solver
    Solver3D( Params &params ) : Solver()
    {
        dt = params.timestep;
        dx = params.cell_length[0];
        dy = params.cell_length[1];
        dz = params.cell_length[2];
        dt_ov_dx = params.timestep / params.cell_length[0];
        dt_ov_dy = params.timestep / params.cell_length[1];
        dt_ov_dz = params.timestep / params.cell_length[2];
        
    };
    virtual ~Solver3D() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    double dt;
    double dx;
    double dy;
    double dz;
    double dt_ov_dx;
    double dt_ov_dy;
    double dt_ov_dz;
    
};//END class

#endif

