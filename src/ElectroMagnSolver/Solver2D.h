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
    Solver2D(PicParams &params) : Solver(params) {
	nx_p = params.n_space[0]+1+2*params.oversize[0];
	nx_d = params.n_space[0]+2+2*params.oversize[0];
	ny_p = params.n_space[1]+1+2*params.oversize[1];
	ny_d = params.n_space[1]+2+2*params.oversize[1];

	dt_ov_dx = params.timestep / params.cell_length[0];
	dt_ov_dy = params.timestep / params.cell_length[1];

    };
    virtual ~Solver2D() {};

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields) = 0;

protected:
    int nx_p;
    int nx_d;
    int ny_p;
    int ny_d;
    double dt_ov_dy;
    double dt_ov_dx;

};//END class

#endif

