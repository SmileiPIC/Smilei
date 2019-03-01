#include "Projector.h"

#include "Params.h"
#include "Patch.h"

Projector::Projector( Params &params, Patch *patch )
    : inv_cell_volume( 1. / params.cell_volume )
{
}

