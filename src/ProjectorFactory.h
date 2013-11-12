
#ifndef PROJECTORFACTORY_H
#define PROJECTORFACTORY_H

#include "Projector.h"
#include "Projector1D2Order.h"
#include "Projector2D2Order.h"

#include "PicParams.h"
#include "SmileiMPI.h"

#include "Tools.h"

class ProjectorFactory {
public:
    static Projector* create(PicParams& params, SmileiMPI* smpi) {
    	Projector* Proj = NULL;
    	if ( ( params.geometry == "1d3v" ) && ( params.interpolation_order == (unsigned int)2 ) ) {
    		Proj = new Projector1D2Order(&params, smpi);
    	}
    	else if ( ( params.geometry == "2d3v" ) && ( params.interpolation_order == (unsigned int)2 ) ) {
    		Proj = new Projector2D2Order(&params, smpi);
   		}
    	else {
    	    ERROR( "Unknwon parameters : " << params.geometry << ", Order : " << params.interpolation_order );
    	}

    	return Proj;
    }

};

#endif
