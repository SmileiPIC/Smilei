
#ifndef FIELDSBC_FACTORY_H
#define FIELDSBC_FACTORY_H

#include "FieldsBC.h"
#include "FieldsBC1D.h"
#include "FieldsBC2D.h"
#include "FieldsBC2D_Damping.h"

#include "PicParams.h"

class FieldsBC_Factory {
public:
    static std::vector<FieldsBC*> create(PicParams& params) {
        std::vector<FieldsBC*> fieldsBoundCond;
        if ( params.geometry == "1d3v" ) {
	    fieldsBoundCond.resize(1, NULL);
            fieldsBoundCond[0] = new FieldsBC1D(&params);
        }
        else if ( params.geometry == "2d3v" ) {
	    fieldsBoundCond.resize(1, NULL);

	    fieldsBoundCond[0] = new FieldsBC2D(&params);
	    if (!params.use_transverse_periodic) {
		fieldsBoundCond.resize(2, NULL);
		fieldsBoundCond[1] = new FieldsBC2D_Damping(&params);
	    }
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return fieldsBoundCond;
    }

};

#endif

