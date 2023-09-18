#ifndef MF_SOLVER3D_YEE_H
#define MF_SOLVER3D_YEE_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver3D_Yee : public Solver3D
{

public:
    //! Creator for MF_Solver3D_Yee
    MF_Solver3D_Yee( Params &params );
    virtual ~MF_Solver3D_Yee();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );

protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;

};//END class

#endif
