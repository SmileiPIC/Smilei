#ifndef MA_SOLVERAM_FRIEDMAN_H
#define MA_SOLVERAM_FRIEDMAN_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_SolverAM_Friedman : public SolverAM
{

public:
    MA_SolverAM_Friedman( Params &params );
    virtual ~MA_SolverAM_Friedman();

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
