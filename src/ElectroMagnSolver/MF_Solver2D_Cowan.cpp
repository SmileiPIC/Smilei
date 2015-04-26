
#include "MF_Solver2D_Cowan.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MF_Solver2D_Cowan::MF_Solver2D_Cowan(PicParams &params)
    : Solver2D(params)
{
}

MF_Solver2D_Cowan::~MF_Solver2D_Cowan()
{
}

void MF_Solver2D_Cowan::operator() ( ElectroMagn* fields )
{
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(fields->Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(fields->Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(fields->Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(fields->Bx_);
    Field2D* By2D = static_cast<Field2D*>(fields->By_);
    Field2D* Bz2D = static_cast<Field2D*>(fields->Bz_);

    // ...
}

