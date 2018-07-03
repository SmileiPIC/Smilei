#ifndef MA_SOLVERRZ_NORM_H
#define MA_SOLVERRZ_NORM_H

#include "SolverRZ.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_SolverRZ_norm : public SolverRZ
{

public:
    MA_SolverRZ_norm(Params &params);
    virtual ~MA_SolverRZ_norm();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);

protected:

};//END class

#endif

