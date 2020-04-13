#ifndef MF_SOLVER2D_BOUCHARD_H
#define MF_SOLVER2D_BOUCHARD_H

#include "Solver2D.h" 
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver2D_Bouchard : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Bouchard
    MF_Solver2D_Bouchard(Params &params);
    virtual ~MF_Solver2D_Bouchard();

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
