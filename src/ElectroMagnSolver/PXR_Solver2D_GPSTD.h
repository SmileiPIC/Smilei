#ifndef PXR_SOLVER2D_GPSTD_H
#define PXR_SOLVER2D_GPSTD_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PXR_Solver2D_GPSTD : public Solver2D
{

public:
    PXR_Solver2D_GPSTD( Params &params );
    virtual ~PXR_Solver2D_GPSTD();
    
    void coupling( Params &params, ElectroMagn *EMfields, bool full_domain = false ) override;
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) override;
    
protected:

};//END class

#endif

