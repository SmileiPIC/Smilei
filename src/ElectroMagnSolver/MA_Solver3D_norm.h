#ifndef MA_SOLVER3D_NORM_H
#define MA_SOLVER3D_NORM_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver3D_norm : public Solver3D
{

public:
    MA_Solver3D_norm( Params &params );
    virtual ~MA_Solver3D_norm();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:

};//END class

#endif

