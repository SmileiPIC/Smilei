#include "Interpolator1D.h"

#include <cmath>
#include <iostream>

#include "SmileiMPI_Cart1D.h"
#include "Patch.h"

using namespace std;

Interpolator1D::Interpolator1D(PicParams &params, SmileiMPI *smpi, Patch* patch)
  : Interpolator(params, smpi, patch) {

    index_domain_begin = patch->cell_starting_global_index[0];

}

