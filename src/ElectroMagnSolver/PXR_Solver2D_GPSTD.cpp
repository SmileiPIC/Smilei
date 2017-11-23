
#include "PXR_Solver2D_GPSTD.h"

#include "ElectroMagn.h"
#include "Field2D.h"
#include "interface.h"

PXR_Solver2D_GPSTD::PXR_Solver2D_GPSTD(Params &params)
: Solver2D(params)
{
}

PXR_Solver2D_GPSTD::~PXR_Solver2D_GPSTD()
{
}

void PXR_Solver2D_GPSTD::operator() ( ElectroMagn* fields )
{
    duplicate_field_into_pxr( fields );
    push_psatd_ebfield_3d_();
    duplicate_field_into_smilei( fields );

}

