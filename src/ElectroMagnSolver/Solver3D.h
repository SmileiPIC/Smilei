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
        std::vector<unsigned int> n_space(params.n_space);
        if (params.uncoupled_grids)
            n_space = params.n_space_region;

        std::vector<unsigned int> oversize(params.oversize);
        if (params.uncoupled_grids)
            oversize = params.region_oversize;
        
	nx_p = n_space[0] +1+2*oversize[0];
	nx_d = n_space[0] +2+2*oversize[0]-(params.is_pxr);
	ny_p = n_space[1] +1+2*oversize[1];
	ny_d = n_space[1] +2+2*oversize[1]-(params.is_pxr);
	nz_p = n_space[2] +1+2*oversize[2];
	nz_d = n_space[2] +2+2*oversize[2]-(params.is_pxr);
        
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

