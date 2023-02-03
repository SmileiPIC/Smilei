
#include "ElectroMagnBC2D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"


ElectroMagnBC2D::ElectroMagnBC2D( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC( params, patch, i_boundary )
{
}

ElectroMagnBC2D::~ElectroMagnBC2D()
{
}
