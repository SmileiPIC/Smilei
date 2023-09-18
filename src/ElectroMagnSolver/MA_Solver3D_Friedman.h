#ifndef MA_SOLVER3D_FRIEDMAN_H
#define MA_SOLVER3D_FRIEDMAN_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver3D_Friedman : public Solver3D
{

public:
    MA_Solver3D_Friedman( Params &params );
    virtual ~MA_Solver3D_Friedman();

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
