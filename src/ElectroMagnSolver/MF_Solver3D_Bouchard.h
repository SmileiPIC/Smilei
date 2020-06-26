#ifndef MF_SOLVER3D_BOUCHARD_H
#define MF_SOLVER3D_BOUCHARD_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver3D_Bouchard : public Solver3D
{

public:
    //! Creator for MF_Solver3D_Bouchard
    MF_Solver3D_Bouchard( Params &params );
    virtual ~MF_Solver3D_Bouchard();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // Parameters for the Maxwell-Faraday solver
    double dx;
    double dy;
    double dz;
    double delta_x;
    double delta_y;
    double delta_z;
    double beta_xy;
    double beta_yx;
    double beta_xz;
    double beta_zx;
    double beta_yz;
    double beta_zy;
    double alpha_x;
    double alpha_y;
    double alpha_z;

    double Ax  ;
    double Ay  ;
    double Az  ;
    double Bxy ;
    double Byx ;
    double Bxz ;
    double Bzx ;
    double Byz ;
    double Bzy ;
    double Dx  ;
    double Dy  ;
    double Dz  ;
    
protected:

};//END class

#endif
