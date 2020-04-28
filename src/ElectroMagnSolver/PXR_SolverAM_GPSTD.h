#ifndef PXR_SOLVERAM_GPSTD_H
#define PXR_SOLVERAM_GPSTD_H

#include "SolverAM.h"
class ElectroMagn;
class cField3D;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PXR_SolverAM_GPSTD : public SolverAM
{

public:
    PXR_SolverAM_GPSTD( Params &params );
    virtual ~PXR_SolverAM_GPSTD();
    
    void coupling( Params &params, ElectroMagn *EMfields, bool full_domain = false ) override;
    void uncoupling() override;
    void rotational_cleaning( ElectroMagn *fields ) override;
    void densities_correction( ElectroMagn *fields ) override;
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) override;
    void _2Dvectors_to_3D( ElectroMagn *fields );
    void _3D_to_2Dvectors( ElectroMagn *fields );

protected:
    cField3D* El_pxr;
    cField3D* Er_pxr;
    cField3D* Et_pxr;
    cField3D* Bl_pxr;
    cField3D* Br_pxr;
    cField3D* Bt_pxr;
    cField3D* Jl_pxr;
    cField3D* Jr_pxr;
    cField3D* Jt_pxr;
    cField3D* rho_pxr;
    cField3D* rhoold_pxr;

};//END class

#endif

