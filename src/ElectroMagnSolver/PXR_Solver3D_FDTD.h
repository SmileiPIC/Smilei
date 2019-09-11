#ifndef PXR_SOLVER3D_FDTD_H
#define PXR_SOLVER3D_FDTD_H


#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PXR_Solver3D_FDTD : public Solver3D
{

public:
    PXR_Solver3D_FDTD( Params &params );
    virtual ~PXR_Solver3D_FDTD();
    
    void coupling( Params &params, ElectroMagn *EMfields, bool full_domain = false ) override;
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) override;
    
protected:

};//END class

#endif

