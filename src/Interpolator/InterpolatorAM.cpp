#include "InterpolatorAM.h"

#include <cmath>
#include <iostream>

#include "Patch.h"

using namespace std;

InterpolatorAM::InterpolatorAM(Params &params, Patch* patch)
  : Interpolator(params, patch) {

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);

}

