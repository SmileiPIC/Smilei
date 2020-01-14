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
    Solver1D( Params &params ) : Solver( params )
    {
        std::vector<unsigned int> n_space(params.n_space);
        if (params.uncoupled_grids)
            n_space = params.n_space_region;
        nx_p = n_space[0] +1+2*params.oversize[0];
        nx_d = n_space[0] +2+2*params.oversize[0];
        
        dt = params.timestep;
        dt_ov_dx = params.timestep / params.cell_length[0];
    };
    virtual ~Solver1D() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    unsigned int nx_p;
    unsigned int nx_d;
    double dt;
    double dt_ov_dx;
    
};//END class

#endif

