#ifndef MA_SOLVER1D_NORM_H
#define MA_SOLVER1D_NORM_H

#include "Solver1D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver1D_norm : public Solver1D
{

public:
    //! Creator for MF_Solver1D_Yee
    MA_Solver1D_norm( Params &params );
    virtual ~MA_Solver1D_norm();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:

};//END class

#endif

