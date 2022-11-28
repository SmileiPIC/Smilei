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
    SolverAM( Params &params ) : Solver()
    {
        Nmode= params.nmodes;
        dt = params.timestep;
        dl = params.cell_length[0];
        dr = params.cell_length[1];
        dt_ov_dl = params.timestep / params.cell_length[0];
        dt_ov_dr = params.timestep / params.cell_length[1];
        
    };
    virtual ~SolverAM() {};
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) = 0;
    
protected:
    unsigned int Nmode;
    double dt;
    double dl;
    double dr;
    double dt_ov_dl;
    double dt_ov_dr;
    
};//END class

#endif

