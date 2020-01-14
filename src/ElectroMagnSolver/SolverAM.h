#ifndef SOLVERAM_H
#define SOLVERAM_H

#include "Solver.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SolverAM
//  --------------------------------------------------------------------------------------------------------------------
class SolverAM : public Solver
{

public:
    //! Creator for Solver
    SolverAM( Params &params ) : Solver( params )
    {
        std::vector<unsigned int> oversize(params.oversize);
        if (params.uncoupled_grids)
            oversize = params.region_oversize;
        nl_p = params.n_space[0]+1+2*oversize[0];
        nl_d = params.n_space[0]+2+2*oversize[0];
        nr_p = params.n_space[1]+1+2*oversize[1];
        nr_d = params.n_space[1]+2+2*oversize[1];
        
        std::vector<unsigned int> n_space(params.n_space);
        if (params.uncoupled_grids)
            n_space = params.n_space_region;

        nl_p = n_space[0] +1+2*oversize[0];
        nl_d = n_space[0] +2+2*oversize[0]-(params.is_pxr);
        nr_p = n_space[1] +1+2*oversize[1];
        nr_d = n_space[1] +2+2*oversize[1]-(params.is_pxr);


        Nmode= params.nmodes;
        dt = params.timestep;
        dr = params.cell_length[1];
        dt_ov_dl = params.timestep / params.cell_length[0];
        dt_ov_dr = params.timestep / params.cell_length[1];
        
    };
    virtual ~SolverAM() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    unsigned int nl_p;
    unsigned int nl_d;
    unsigned int nr_p;
    unsigned int nr_d;
    unsigned int Nmode;
    double dt;
    double dr;
    double dt_ov_dl;
    double dt_ov_dr;
    
};//END class

#endif

