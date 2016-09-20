#ifndef PROJECTORFACTORY_H
#define PROJECTORFACTORY_H

#include "Projector.h"
#include "Projector1D2Order.h"
#include "Projector1D4Order.h"
#include "Projector2D2Order.h"
#include "Projector2D4Order.h"
#include "Projector3D2Order.h"

#include "Params.h"
#include "Patch.h" 

#include "Tools.h"

class ProjectorFactory {
public:
  static Projector* create(Params& params, Patch* patch) {
        Projector* Proj = NULL;
        // ---------------
        // 1d3v simulation
        // ---------------
        if ( ( params.geometry == "1d3v" ) && ( params.interpolation_order == (unsigned int)2 ) ) {
            Proj = new Projector1D2Order(params, patch);
        }
        else if ( ( params.geometry == "1d3v" ) && ( params.interpolation_order == (unsigned int)4 ) ) {
            Proj = new Projector1D4Order(params, patch);
        }
        // ---------------
        // 2d3v simulation
        // ---------------
        else if ( ( params.geometry == "2d3v" ) && ( params.interpolation_order == (unsigned int)2 ) ) {
            Proj = new Projector2D2Order(params, patch);
        }
        // ---------------
        // 3d3v simulation
        // ---------------
        else if ( ( params.geometry == "3d3v" ) && ( params.interpolation_order == (unsigned int)2 ) ) {
            Proj = new Projector3D2Order(params, patch);
        }
        else {
            ERROR( "Unknwon parameters : " << params.geometry << ", Order : " << params.interpolation_order );
        }

        return Proj;
    }

};

#endif
