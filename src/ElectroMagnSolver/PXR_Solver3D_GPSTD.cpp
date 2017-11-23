
#include "PXR_Solver3D_GPSTD.h"

#include "ElectroMagn.h"
#include "Field3D.h"
#include "interface.h"

PXR_Solver3D_GPSTD::PXR_Solver3D_GPSTD(Params &params)
: Solver3D(params)
{
}

PXR_Solver3D_GPSTD::~PXR_Solver3D_GPSTD()
{
}

void PXR_Solver3D_GPSTD::operator() ( ElectroMagn* fields )
{
    duplicate_field_into_pxr( fields );
    push_psatd_ebfield_3d_();
    duplicate_field_into_smilei( fields );

}

