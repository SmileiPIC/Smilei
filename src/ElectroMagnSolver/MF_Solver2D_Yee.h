#ifndef MF_SOLVER2D_YEE_H
#define MF_SOLVER2D_YEE_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver2D_Yee : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Yee
    MF_Solver2D_Yee( Params &params );
    virtual ~MF_Solver2D_Yee();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;
    
};//END class

#endif

