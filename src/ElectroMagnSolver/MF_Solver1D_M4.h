#ifndef MF_SOLVER1D_M4_H
#define MF_SOLVER1D_M4_H

#include "Solver1D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver1D_M4 : public Solver1D
{

public:
    //! Creator for MF_Solver1D_M4
    MF_Solver1D_M4( Params &params );
    virtual ~MF_Solver1D_M4();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // Parameters for the Maxwell-Faraday solver
    double dx;
    double delta_x;
    double alpha_x;

    double Ax  ;
    double Dx  ;
    
protected:

};//END class

#endif

