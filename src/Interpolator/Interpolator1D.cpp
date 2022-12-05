#include "Interpolator1D.h"

#include <cmath>
#include <iostream>

#include "Patch.h"

using namespace std;

Interpolator1D::Interpolator1D( Patch *patch )
    : Interpolator()
{

    index_domain_begin =  patch->getCellStartingGlobalIndex( 0 );
    
}

