
#ifndef FIELDSBC_FACTORY_H
#define FIELDSBC_FACTORY_H

#include "FieldsBC.h"
#include "FieldsBC1D_SM.h"
#include "FieldsBC2D_Long_SM.h"
#include "FieldsBC2D_Trans_SM.h"
#include "FieldsBC2D_Trans_Damping.h"

#include "PicParams.h"


class FieldsBC_Factory {
public:
    static std::vector<FieldsBC*> create(PicParams& params, LaserParams &laser_params) {
        std::vector<FieldsBC*> fieldsBoundCond;

        if ( params.geometry == "1d3v" ) {
            fieldsBoundCond.resize(1, NULL);
            if ( params.bc_em_type_long == "silver-muller" )
		fieldsBoundCond[0] = new FieldsBC1D_SM(params, laser_params);
            else if ( params.bc_em_type_long != "periodic" ) {
                // If periodic : !applied -> NULL
                ERROR( "Unknwon boundary condition : " << params.bc_em_type_long );
            }
        }

        else if ( params.geometry == "2d3v" ) {
            fieldsBoundCond.resize(2, NULL);

            if ( params.bc_em_type_long == "silver-muller" ) {
                //Boundary in the X direction is set to Silver-Muller
                fieldsBoundCond[0] = new FieldsBC2D_Long_SM(params, laser_params);
            }
            else if ( params.bc_em_type_long != "periodic" ) {
                // If periodic : !applied -> NULL
                ERROR( "Unknwon boundary condition : " << params.bc_em_type_long );
            }


            if ( params.bc_em_type_trans == "silver-muller" ) {
                //Boundary in the Y direction is set to Silver-Muller.
                fieldsBoundCond[1] = new FieldsBC2D_Trans_SM(params, laser_params);
            }
            else if ( params.bc_em_type_trans == "damping" ) {
                // Boundary in the Y direction is set to damping if they are not periodic. 
		fieldsBoundCond[1] = new FieldsBC2D_Trans_Damping(params, laser_params); 
            }
            else if ( params.bc_em_type_trans != "periodic" ) {
                // If periodic : !applied -> NULL
                ERROR( "Unknwon boundary condition : " << params.bc_em_type_trans );
	    }
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return fieldsBoundCond;
    }

};

#endif

