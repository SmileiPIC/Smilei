#ifndef MA_SOLVERAM_TERZANI_H
#define MA_SOLVERAM_TERZANI_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MA_SolverAM_Terzani : public SolverAM
{

public:
    MA_SolverAM_Terzani( Params &params );
    virtual ~MA_SolverAM_Terzani();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // coefficient necessary to reduce the dispersion
    double delta;
    
protected:

};//END class

#endif

