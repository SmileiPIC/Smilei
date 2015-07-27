#ifndef INTERPOLATORFACTORY_H
#define INTERPOLATORFACTORY_H

#include "Interpolator.h"
#include "Interpolator1D2Order.h"
#include "Interpolator1D3Order.h"
#include "Interpolator1D4Order.h"
#include "Interpolator2D2Order.h"
#include "Interpolator2D4Order.h"

#include "PicParams.h"
#include "Patch.h"

#include "Tools.h"

class InterpolatorFactory {
public:
  static Interpolator* create(PicParams& params, Patch *patch) {
        Interpolator* Interp = NULL;
        // ---------------
        // 1d3v simulation
        // ---------------
        if ( ( params.geometry == "1d3v" ) && ( params.interpolation_order == 2 ) ) {
            Interp = new Interpolator1D2Order(params, patch);
        }
        else if ( ( params.geometry == "1d3v" ) && ( params.interpolation_order == 4 ) ) {
            Interp = new Interpolator1D4Order(params, patch);
        }
        // ---------------
        // 2d3v simulation
        // ---------------
        else if ( ( params.geometry == "2d3v" ) && ( params.interpolation_order == 2 ) ) {
            Interp = new Interpolator2D2Order(params, patch);
        }
        else {
            ERROR( "Unknwon parameters : " << params.geometry << ", Order : " << params.interpolation_order );
        }

        return Interp;
    }

};

#endif
