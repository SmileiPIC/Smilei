#include "Interpolator2D.h"

#include <cmath>
#include <iostream>

#include "Patch.h"

using namespace std;

Interpolator2D::Interpolator2D( Params &params, Patch *patch )
    : Interpolator( params, patch )
{

    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    
}

