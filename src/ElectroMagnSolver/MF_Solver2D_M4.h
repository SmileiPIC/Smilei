#ifndef MF_SOLVER2D_M4_H
#define MF_SOLVER2D_M4_H

#include "Solver2D.h" 
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver2D_M4 : public Solver2D
{

public:
    //! Creator for MF_Solver2D_M4
    MF_Solver2D_M4(Params &params);
    virtual ~MF_Solver2D_M4();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);
 
    // Parameters for the Maxwell-Faraday solver
    double dx;
    double dy;
    double delta_x;
    double delta_y;
    double beta_xy;
    double beta_yx;
    double alpha_x;
    double alpha_y;
    double Ax;
    double Ay;
    double Bx;
    double By;
    double Dx;
    double Dy;

protected:
    // Check if time filter is applied or not
         bool isEFilterApplied;

};//END class

#endif
