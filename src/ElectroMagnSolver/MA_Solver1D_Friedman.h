#ifndef MA_SOLVER1D_FRIEDMAN_H
#define MA_SOLVER1D_FRIEDMAN_H

#include "Solver1D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_Solver1D_Friedman : public Solver1D
{

public:
    //! Creator for MF_Solver1D_Yee
    MA_Solver1D_Friedman( Params &params );
    virtual ~MA_Solver1D_Friedman();

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
