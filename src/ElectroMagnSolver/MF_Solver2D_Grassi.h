#ifndef MF_SOLVER2D_GRASSI_H
#define MF_SOLVER2D_GRASSI_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver2D_Grassi : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Yee
    MF_Solver2D_Grassi( Params &params );
    virtual ~MF_Solver2D_Grassi();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // Parameters for the Maxwell-Faraday solver
    double dt_ov_dx;
    double dt_ov_dy;
    double dx;
    double dy;
    double Ax;
    double Ay;
    double Dx;
    double Dy;
    
protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;
    
};//END class

#endif

