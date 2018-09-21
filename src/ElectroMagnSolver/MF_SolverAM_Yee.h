#ifndef MF_SOLVERRZ_YEE_H
#define MF_SOLVERRZ_YEE_H

#include "SolverAM.h" 
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_SolverAM_Yee : public SolverAM
{

public:
    //! Creator for MF_SolverAM_Yee
    MF_SolverAM_Yee(Params &params);
    virtual ~MF_SolverAM_Yee();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);

protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;

};//END class

#endif

