#include "Interpolator2D.h"

#include <cmath>
#include <iostream>

#include "Patch.h"

using namespace std;

Interpolator2D::Interpolator2D( Patch *patch )
    : Interpolator()
{

    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    
}

