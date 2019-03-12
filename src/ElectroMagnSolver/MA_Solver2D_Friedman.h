#ifndef MA_SOLVER2D_FRIEDMAN_H
#define MA_SOLVER2D_FRIEDMAN_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver2D_Friedman : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Yee
    MA_Solver2D_Friedman( Params &params );
    virtual ~MA_Solver2D_Friedman();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    //! parameter for time-filtering
    double ftheta;
    double alpha;
    double beta;
    double delta;
    
protected:

};//END class

#endif

