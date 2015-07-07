#ifndef ELECTROMAGNFACTORY_H
#define ELECTROMAGNFACTORY_H

#include "ElectroMagn.h"
#include "ElectroMagn1D.h"
#include "ElectroMagn2D.h"

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Patch.h"

#include "Tools.h"

class ElectroMagnFactory {
public:
    static ElectroMagn* create(PicParams& params,  LaserParams &laser_params, Patch* patch) {
        ElectroMagn* EMfields = NULL;
        if ( params.geometry == "1d3v" ) {
	  EMfields = new ElectroMagn1D(params, laser_params, patch);
        }
        else if ( params.geometry == "2d3v" ) {
            EMfields = new ElectroMagn2D(params, laser_params, patch);
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return EMfields;
    }

};

#endif

