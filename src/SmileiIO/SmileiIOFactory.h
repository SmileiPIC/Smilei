#ifndef SMILEIIOFACTORY_H
#define SMILEIIOFACTORY_H

#include "SmileiIO.h"
#include "SmileiIO_Cart1D.h"
#include "SmileiIO_Cart2D.h"

#include "PicParams.h"
#include "SmileiMPI.h"

#include "Tools.h"

class SmileiIOFactory {
public:
    static SmileiIO* create(PicParams& params, DiagParams& diagParams, SmileiMPI* smpi) {
        SmileiIO* sio = NULL;
        if ( params.geometry == "1d3v" ) {
            sio = new  SmileiIO_Cart1D(params, diagParams, smpi);
        }
        else if ( params.geometry == "2d3v" ) {
            sio = new  SmileiIO_Cart2D(params, diagParams, smpi);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }

//    	// Creation of a cartesian topology
//    	smpi->createTopology(params);
//
//
//    	if ( params.geometry == "2d3v" ) smpi->createType(params);

        return sio;
    }

};

#endif

