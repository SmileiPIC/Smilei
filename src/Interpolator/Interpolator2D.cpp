#include "Interpolator2D.h"

#include <cmath>
#include <iostream>

#include "SmileiMPI_Cart2D.h"
#include "Patch.h"

using namespace std;

Interpolator2D::Interpolator2D(PicParams &params, SmileiMPI *smpi, Patch* patch)
  : Interpolator(params, smpi, patch) {

    i_domain_begin = patch->cell_starting_global_index[0];
    j_domain_begin = patch->cell_starting_global_index[1];

}

