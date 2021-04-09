
#include "ElectroMagnBCAM.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "cField2D.h"


ElectroMagnBCAM::ElectroMagnBCAM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC( params, patch, i_boundary )
{
}

ElectroMagnBCAM::~ElectroMagnBCAM()
{
}


void ElectroMagnBCAM::applyBConEdges( ElectroMagn *EMfields, Patch *patch )
{
}
