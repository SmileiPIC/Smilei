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
    Solver3D( Params &params ) : Solver( params )
    {
        nx_p = params.n_space[0] * params.global_factor[0]+1+2*params.oversize[0];
        nx_d = params.n_space[0] * params.global_factor[0]+2+2*params.oversize[0];
        ny_p = params.n_space[1] * params.global_factor[1]+1+2*params.oversize[1];
        ny_d = params.n_space[1] * params.global_factor[1]+2+2*params.oversize[1];
        nz_p = params.n_space[2] * params.global_factor[2]+1+2*params.oversize[2];
        nz_d = params.n_space[2] * params.global_factor[2]+2+2*params.oversize[2];
        
        dt = params.timestep;
        dt_ov_dx = params.timestep / params.cell_length[0];
        dt_ov_dy = params.timestep / params.cell_length[1];
        dt_ov_dz = params.timestep / params.cell_length[2];
        
    };
    virtual ~Solver3D() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    unsigned int nx_p;
    unsigned int nx_d;
    unsigned int ny_p;
    unsigned int ny_d;
    unsigned int nz_p;
    unsigned int nz_d;
    double dt;
    double dt_ov_dx;
    double dt_ov_dy;
    double dt_ov_dz;
    
};//END class

#endif

