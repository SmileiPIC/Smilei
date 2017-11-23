#ifndef PXR_SOLVER3D_GPSTD_H
#define PXR_SOLVER3D_GPSTD_H

#include "Solver3D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PXR_Solver3D_GPSTD : public Solver3D
{

public:
    PXR_Solver3D_GPSTD(Params &params);
    virtual ~PXR_Solver3D_GPSTD();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn* fields);

protected:

};//END class

#endif

