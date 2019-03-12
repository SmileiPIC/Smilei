#ifndef MF_SOLVER2D_LEHE_H
#define MF_SOLVER2D_LEHE_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver2D_Lehe : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Lehe
    MF_Solver2D_Lehe( Params &params );
    virtual ~MF_Solver2D_Lehe();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // Parameters for the Maxwell-Faraday solver
    double dx;
    double dy;
    double delta_x;
    double beta_xy;
    double beta_yx;
    double alpha_x;
    double alpha_y;
    
protected:

};//END class

#endif
