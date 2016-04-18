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
        
        PyTools::extract("fieldsToDump", sio->fieldsToDump);
        
        int ntime_step_avg(0);
        if( !PyTools::extract("ntime_step_avg", ntime_step_avg) )
            sio->dumpAvgFields_ = false;
        else
            sio->dumpAvgFields_ = true;
        
        sio->field_timeSelection    = new TimeSelection( PyTools::extract_py("fieldDump_every"),    "Fields" );
        sio->avgfield_timeSelection = new TimeSelection( PyTools::extract_py("avgfieldDump_every"), " Average fields" );
        
        return sio;
    }
    
    
    // Cloning 
    static SmileiIO* clone( SmileiIO* sio, Params& params, Patch* patch )
    {
        SmileiIO* newSio = NULL;
        if ( params.geometry == "1d3v" ) {
            newSio = new  SmileiIO_Cart1D(params, patch);
        }
        else if ( params.geometry == "2d3v" ) {
            newSio = new  SmileiIO_Cart2D(params, patch);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }
        
        newSio->fieldsToDump   = sio->fieldsToDump;
        newSio->dumpAvgFields_ = sio->dumpAvgFields_;
        newSio->field_timeSelection    = new TimeSelection(sio->field_timeSelection);
        newSio->avgfield_timeSelection = new TimeSelection(sio->field_timeSelection);
        return newSio;
    }

};

#endif

