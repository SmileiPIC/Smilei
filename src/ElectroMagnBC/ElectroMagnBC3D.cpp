
#include "ElectroMagnBC3D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"


ElectroMagnBC3D::ElectroMagnBC3D( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC( params, patch, i_boundary )
{
}

ElectroMagnBC3D::~ElectroMagnBC3D()
{
}

