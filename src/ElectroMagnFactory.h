
#ifndef ELECTROMAGNFACTORY_H
#define ELECTROMAGNFACTORY_H

#include "ElectroMagn.h"
#include "ElectroMagn1D.h"
#include "ElectroMagn2D.h"

#include "PicParams.h"
#include "SmileiMPI.h"

#include "Tools.h"

class ElectroMagnFactory {
public:
    static ElectroMagn* create(PicParams& params, SmileiMPI* smpi) {
        ElectroMagn* EMfields = NULL;
        if ( params.geometry == "1d3v" ) {
            EMfields = new ElectroMagn1D(&params, smpi);
        }
        else if ( params.geometry == "2d3v" ) {
            EMfields = new ElectroMagn2D(&params, smpi);
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return EMfields;
    }

};

#endif

