#ifndef MF_SOLVER1D_YEE_H
#define MF_SOLVER1D_YEE_H

#include "Solver1D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_Solver1D_Yee : public Solver1D
{

public:
    //! Creator for MF_Solver1D_Yee
    MF_Solver1D_Yee( Params &params );
    virtual ~MF_Solver1D_Yee();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
protected:

};//END class

#endif

