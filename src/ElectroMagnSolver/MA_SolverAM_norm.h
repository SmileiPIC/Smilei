#ifndef MA_SOLVERRZ_NORM_H
#define MA_SOLVERRZ_NORM_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_SolverAM_norm : public SolverAM
{

public:
    MA_SolverAM_norm(Params &params);
    virtual ~MA_SolverAM_norm();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);

protected:

};//END class

#endif

