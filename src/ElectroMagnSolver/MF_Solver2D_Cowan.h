#ifndef MF_SOLVER2D_COWAN_H
#define MF_SOLVER2D_COWAN_H

#include "Solver2D.h" 
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver2D_Cowan : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Cowan
    MF_Solver2D_Cowan(Params &params);
    virtual ~MF_Solver2D_Cowan();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);
 
    // Parameters for the Maxwell-Faraday solver
    double dx;
    double dy;
    double beta_x;
    double beta_y;
    double Beta_x;
    double Beta_y;

protected:

};//END class

#endif
