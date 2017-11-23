
#include "PXR_Solver3D_FDTD.h"

#include "ElectroMagn.h"
#include "Field3D.h"
#include "interface.h"

PXR_Solver3D_FDTD::PXR_Solver3D_FDTD(Params &params)
: Solver3D(params)
{
}

PXR_Solver3D_FDTD::~PXR_Solver3D_FDTD()
{
}

void PXR_Solver3D_FDTD::operator() ( ElectroMagn* fields )
{
    duplicate_field_into_pxr( fields );
    solve_maxwell_fdtd_pxr();
    duplicate_field_into_smilei( fields );

}

