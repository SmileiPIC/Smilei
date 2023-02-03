
#include "ElectroMagnBC1D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field1D.h"


ElectroMagnBC1D::ElectroMagnBC1D( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC( params, patch, i_boundary )
{
}

ElectroMagnBC1D::~ElectroMagnBC1D()
{
}
