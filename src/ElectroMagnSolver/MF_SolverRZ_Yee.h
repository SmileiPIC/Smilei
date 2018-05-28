#ifndef MF_SOLVERRZ_YEE_H
#define MF_SOLVERRZ_YEE_H

#include "SolverRZ.h" 
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_SolverRZ_Yee : public SolverRZ
{

public:
    //! Creator for MF_SolverRZ_Yee
    MF_SolverRZ_Yee(Params &params);
    virtual ~MF_SolverRZ_Yee();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);

protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;

};//END class

#endif

