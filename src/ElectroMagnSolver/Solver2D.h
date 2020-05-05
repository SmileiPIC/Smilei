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
    Solver2D( Params &params ) : Solver( params )
    {
        std::vector<unsigned int> n_space(params.n_space);
        if (params.uncoupled_grids)
            n_space = params.n_space_region;

        std::vector<unsigned int> oversize(params.oversize);
        if (params.uncoupled_grids)
            oversize = params.region_oversize;

        nx_p = n_space[0] +1+2*oversize[0];
        nx_d = n_space[0] +2+2*oversize[0];
        ny_p = n_space[1] +1+2*oversize[1];
        ny_d = n_space[1] +2+2*oversize[1];

        dt = params.timestep;
        dt_ov_dx = params.timestep / params.cell_length[0];
        dt_ov_dy = params.timestep / params.cell_length[1];
    };
    virtual ~Solver2D() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    unsigned int nx_p;
    unsigned int nx_d;
    unsigned int ny_p;
    unsigned int ny_d;
    double dt;
    double dt_ov_dy;
    double dt_ov_dx;
    
};//END class

#endif

