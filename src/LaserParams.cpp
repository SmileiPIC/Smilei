#include "LaserParams.h"

#include <cmath>

#include "Tools.h"

using namespace std;

LaserParams::LaserParams(PicParams& params, InputData &ifile) {
	
    // -----------------
    // Lasers properties
    // -----------------
    n_laser=0;
    while (ifile.existGroup("laser",n_laser)) {
        
        
        ifile.extract("boxSide",tmpLaser.boxSide,"laser",0,n_laser);
        if ( (tmpLaser.boxSide!="west") && (tmpLaser.boxSide!="east") ) {
            ERROR("At the moment laser can enter only from West/East sides: boxSide \""
                  << tmpLaser.boxSide << "\" not defined");
        }
        LaserStructure tmpLaser;
        
        ifile.extract("a0",tmpLaser.a0,"laser",0,n_laser);
        ifile.extract("delta",tmpLaser.delta,"laser",0,n_laser);
        
        ifile.extract("focus",tmpLaser.focus,"laser",0,n_laser)
        ifile.extract("angle",tmpLaser.angle ,"laser",0,n_laser);
        
        
        ifile.extract("time_profile",tmpLaser.time_profile ,"laser",0,n_laser);
        ifile.extract("int_params",tmpLaser.int_params ,"laser",0,n_laser);
        ifile.extract("double_params",tmpLaser.double_params ,"laser",0,n_laser);
        ifile.extract("transv_profile",tmpLaser.transv_profile ,"laser",0,n_laser);
        ifile.extract("int_params_transv",tmpLaser.int_params_transv ,"laser",0,n_laser);
        ifile.extract("double_params_transv",tmpLaser.double_params_transv ,"laser",0,n_laser);
        
        for (unsigned int i=0; i<tmpLaser.double_params.size(); i++)
            tmpLaser.double_params[i] *= 2.0*M_PI;
        for (unsigned int i=0; i<tmpLaser.double_params_transv.size(); i++)
            tmpLaser.double_params_transv[i] *= 2.0*M_PI;
        
        laser_param.push_back(tmpLaser);
        n_laser++;
    }
    
    
    // Laser related parameters
    // ------------------------
    MESSAGE("Laser related parameters");
    MESSAGE(1,"n_laser        : " << n_laser);
    for ( unsigned int i=0 ; i<n_laser ; i++ ) {
        MESSAGE(2,"laser " << i << ": (boxSide, a0) : (" << laser_param[i].boxSide <<  ", " << laser_param[i].a0 <<  ")");
    }
    
	
}

