
#ifndef ELECTROMAGNBC_FACTORY_H
#define ELECTROMAGNBC_FACTORY_H

#include "ElectroMagnBC.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC2D_Long_SM.h"
#include "ElectroMagnBC2D_Trans_SM.h"
#include "ElectroMagnBC2D_Trans_Damping.h"

#include "PicParams.h"


class ElectroMagnBC_Factory {
public:
    static std::vector<ElectroMagnBC*> create(PicParams& params, LaserParams &laser_params) {
        std::vector<ElectroMagnBC*> emBoundCond;

        if ( params.geometry == "1d3v" ) {
            emBoundCond.resize(1, NULL);
            if ( params.bc_em_type_long == "silver-muller" )
		emBoundCond[0] = new ElectroMagnBC1D_SM(params, laser_params);
            else if ( params.bc_em_type_long != "periodic" ) {
                // If periodic : !applied -> NULL
                ERROR( "Unknwon boundary condition : " << params.bc_em_type_long );
            }
        }

        else if ( params.geometry == "2d3v" ) {
            emBoundCond.resize(2, NULL);

            if ( params.bc_em_type_long == "silver-muller" ) {
                //Boundary in the X direction is set to Silver-Muller
                emBoundCond[0] = new ElectroMagnBC2D_Long_SM(params, laser_params);
            }
            else if ( params.bc_em_type_long != "periodic" ) {
                // If periodic : !applied -> NULL
                ERROR( "Unknwon boundary condition : " << params.bc_em_type_long );
            }


            if ( params.bc_em_type_trans == "silver-muller" ) {
                //Boundary in the Y direction is set to Silver-Muller.
                emBoundCond[1] = new ElectroMagnBC2D_Trans_SM(params, laser_params);
            }
            else if ( params.bc_em_type_trans == "damping" ) {
                // Boundary in the Y direction is set to damping if they are not periodic. 
		emBoundCond[1] = new ElectroMagnBC2D_Trans_Damping(params, laser_params); 
            }
            else if ( params.bc_em_type_trans != "periodic" ) {
                // If periodic : !applied -> NULL
                ERROR( "Unknwon boundary condition : " << params.bc_em_type_trans );
	    }
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return emBoundCond;
    }

};

#endif

