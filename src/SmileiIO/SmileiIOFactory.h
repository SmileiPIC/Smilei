#ifndef SMILEIIOFACTORY_H
#define SMILEIIOFACTORY_H

#include "SmileiIO.h"
#include "SmileiIO_Cart1D.h"
#include "SmileiIO_Cart2D.h"

#include "Params.h"
#include "Patch.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIOFactory
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIOFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate IO environment for the geometry 
    //! \param params : Parameters
    //! \param smpi : MPI environment
    //  --------------------------------------------------------------------------------------------------------------------
    static SmileiIO* create(Params& params, Patch* patch) {
        SmileiIO* sio = NULL;
        if ( params.geometry == "1d3v" ) {
            sio = new  SmileiIO_Cart1D(params, patch);
        }
        else if ( params.geometry == "2d3v" ) {
            sio = new  SmileiIO_Cart2D(params, patch);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }

        return sio;
    }

};

#endif

