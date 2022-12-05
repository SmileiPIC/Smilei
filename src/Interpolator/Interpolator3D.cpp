#include "Interpolator3D.h"

#include <cmath>
#include <iostream>

#include "Patch.h"

using namespace std;

Interpolator3D::Interpolator3D( Patch *patch )
    : Interpolator()
{

    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    k_domain_begin = patch->getCellStartingGlobalIndex( 2 );
    
}

