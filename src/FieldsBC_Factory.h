
#ifndef FIELDSBC_FACTORY_H
#define FIELDSBC_FACTORY_H

#include "FieldsBC.h"
#include "FieldsBC1D.h"
#include "FieldsBC2D.h"
#include "FieldsBC2D_Damping.h"

#include "PicParams.h"

class FieldsBC_Factory {
public:
    static FieldsBC* create(PicParams& params) {
        FieldsBC* fieldsBoundCond = NULL;
        if ( params.geometry == "1d3v" ) {
            fieldsBoundCond = new FieldsBC1D(&params);
        }
        else if ( params.geometry == "2d3v" ) {
            if (params.use_transverse_periodic)
		fieldsBoundCond = new FieldsBC2D(&params);
	    else // Outgoing
                fieldsBoundCond = new FieldsBC2D_Damping(&params);
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return fieldsBoundCond;
    }

};

#endif

