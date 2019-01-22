#ifndef MF_SOLVER3D_LEHE_H
#define MF_SOLVER3D_LEHE_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver3D_Lehe : public Solver3D
{

public:
    //! Creator for MF_Solver3D_Lehe
    MF_Solver3D_Lehe( Params &params );
    virtual ~MF_Solver3D_Lehe();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // Parameters for the Maxwell-Faraday solver
    double dx;
    double dy;
    double dz;
    double delta_x;
    double beta_xy;
    double beta_xz;
    double beta_yx;
    double alpha_x;
    double alpha_y;
    
protected:

};//END class

#endif
