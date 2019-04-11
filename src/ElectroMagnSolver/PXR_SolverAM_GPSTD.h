#ifndef PXR_SOLVERAM_GPSTD_H
#define PXR_SOLVERAM_GPSTD_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PXR_SolverAM_GPSTD : public SolverAM
{

public:
    PXR_SolverAM_GPSTD( Params &params );
    virtual ~PXR_SolverAM_GPSTD();
    
    void coupling( Params &params, ElectroMagn *EMfields ) override;
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) override;
    
protected:

};//END class

#endif

