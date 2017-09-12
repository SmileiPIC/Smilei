#ifndef SOLVERRZ_H
#define SOLVERRZ_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SolverRZ
//  --------------------------------------------------------------------------------------------------------------------
class SolverRZ : public Solver
{

public:
    //! Creator for Solver
    SolverRZ(Params &params) : Solver(params) {
	nx_p = params.n_space[0]+1+2*params.oversize[0];
	nx_d = params.n_space[0]+2+2*params.oversize[0];
	ny_p = params.n_space[1]+1+2*params.oversize[1];
	ny_d = params.n_space[1]+2+2*params.oversize[1];

    dt = params.timestep;
	dt_ov_dx = params.timestep / params.cell_length[0];
	dt_ov_dy = params.timestep / params.cell_length[1];

    };
    virtual ~SolverRZ() {};

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields) = 0;

protected:
    unsigned int nx_p;
    unsigned int nx_d;
    unsigned int ny_p;
    unsigned int ny_d;
    double dt;
    double dt_ov_dx;
    double dt_ov_dy;

};//END class

#endif

