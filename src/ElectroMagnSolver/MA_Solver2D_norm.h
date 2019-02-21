#ifndef MA_SOLVER2D_NORM_H
#define MA_SOLVER2D_NORM_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver2D_norm : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Yee
    MA_Solver2D_norm( Params &params );
    virtual ~MA_Solver2D_norm();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:

};//END class

#endif

