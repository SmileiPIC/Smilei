#ifndef MF_SOLVERAM_TERZANI_H
#define MF_SOLVERAM_TERZANI_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_SolverAM_Terzani : public SolverAM
{

public:
    //! Creator for MF_SolverAM_Terzani
    MF_SolverAM_Terzani( Params &params );
    virtual ~MF_SolverAM_Terzani();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // coefficient necessary to reduce the dispersion
    double delta;
    
protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;
    
};//END class

#endif

