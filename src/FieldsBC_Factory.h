
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
    static std::vector<FieldsBC*> create(PicParams& params) {
        std::vector<FieldsBC*> fieldsBoundCond;
        if ( params.geometry == "1d3v" ) {
	    fieldsBoundCond.resize(1, NULL);
            fieldsBoundCond[0] = new FieldsBC1D_SM(params);
        }
        else if ( params.geometry == "2d3v" ) {
	    fieldsBoundCond.resize(1, NULL);

	    fieldsBoundCond[0] = new FieldsBC2D_Long_SM(params); //Boundary in the X direction is set to Silver-Muller.
	    if (!params.use_transverse_periodic) {
		fieldsBoundCond.resize(2, NULL);
                // Boundary in the Y direction is set to damping if they are not periodic. 
		//fieldsBoundCond[1] = new FieldsBC2D_Trans_Damping(&params); 
                // Boundary in the Y direction is set to SM if they are not periodic.
	        fieldsBoundCond[1] = new FieldsBC2D_Trans_SM(params); //Boundary in the Y direction is set to Silver-Muller.
	    }
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return fieldsBoundCond;
    }

};

#endif

