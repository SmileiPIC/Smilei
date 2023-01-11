#include "InterpolatorAM.h"

#include <cmath>
#include <iostream>

#include "Patch.h"

using namespace std;

InterpolatorAM::InterpolatorAM( Patch *patch )
    : Interpolator()
{

    i_domain_begin_ = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin_ = patch->getCellStartingGlobalIndex( 1 );
    
}

