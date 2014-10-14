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
        LaserStructure tmpLaser;
        
        // side from which the laser enters the simulation box (only west/east at the moment)
        ifile.extract("boxSide",tmpLaser.boxSide,"laser",0,n_laser);
        if ( (tmpLaser.boxSide!="west") && (tmpLaser.boxSide!="east") ) {
            ERROR("At the moment laser can enter only from West/East sides: boxSide \""
                  << tmpLaser.boxSide << " not defined");
        }
        
        // laser intensity
        ifile.extract("a0",tmpLaser.a0,"laser",0,n_laser);
        
        // laser ellipticity/polarization parameter
        ifile.extract("delta",tmpLaser.delta,"laser",0,n_laser);
        
        // position of the laser focus
        tmpLaser.isFocused = ifile.extract("focus",tmpLaser.focus,"laser",0,n_laser);
        
        // incident angle
        ifile.extract("angle",tmpLaser.angle ,"laser",0,n_laser);
        
        // laser time-profile & associated parameters
        ifile.extract("time_profile",tmpLaser.time_profile ,"laser",0,n_laser);
        ifile.extract("int_params",tmpLaser.int_params ,"laser",0,n_laser);
        ifile.extract("double_params",tmpLaser.double_params ,"laser",0,n_laser);
        
        // laser transverse-profile & associated parameters
        ifile.extract("transv_profile",tmpLaser.transv_profile ,"laser",0,n_laser);
        ifile.extract("int_params_transv",tmpLaser.int_params_transv ,"laser",0,n_laser);
        ifile.extract("double_params_transv",tmpLaser.double_params_transv ,"laser",0,n_laser);
        bool delayExists = ifile.extract("delay",tmpLaser.delay ,"laser",0,n_laser);
        
        
        // -------------------------------------
        // Printing out laser related parameters
        // -------------------------------------
        MESSAGE("Laser related parameters");
        MESSAGE(1,"n_laser        : " << n_laser);
        for ( unsigned int i=0 ; i<n_laser ; i++ ) {
            MESSAGE(2,"laser " << i << ": (boxSide, a0) : (" << laser_param[i].boxSide <<  ", " << laser_param[i].a0 <<  ")");
        }
        
        
        // -----------------------------------------------------------------
        // normalization (from wavelength-related units to normalized units)
        // -----------------------------------------------------------------
        for (unsigned int i=0; i<tmpLaser.double_params.size(); i++)
            tmpLaser.double_params[i] *= params.conv_fac;
        for (unsigned int i=0; i<tmpLaser.double_params_transv.size(); i++)
            tmpLaser.double_params_transv[i] *= params.conv_fac;
        if ( (tmpLaser.angle!=0) || (tmpLaser.isFocused) ) {
            for (unsigned int i=0; i<tmpLaser.focus.size(); i++)
                tmpLaser.focus[i] *= params.conv_fac;
        }

        
        // -----------------------------------------------------------------
        // tests on the laser parameters (when arbitrary focus or incidence)
        // -----------------------------------------------------------------
        // at the moment, only Gaussian beam can be used when laser is focused or with a not-normal incidence
        // note that there is also no choice about polarization (for along Bz)
        if ( (tmpLaser.boxSide=="west") && ((tmpLaser.angle!=0) || (tmpLaser.isFocused)) ){
            
            if ( (tmpLaser.int_params_transv.size()==0) ) {
                WARNING("A default cut-off (3 sigma) is applied on laser " << n_laser << " transverse profile");
                tmpLaser.int_params_transv.resize(1);
                tmpLaser.int_params_transv[0] = 3;
            }
            if ( !delayExists ) {
                tmpLaser.delay = tmpLaser.int_params_transv[0]*tmpLaser.double_params_transv[0]
                * sin(tmpLaser.angle*M_PI/180.0);
                if (tmpLaser.delay!=0)
                    WARNING("Introduction of a time-delay of " << tmpLaser.delay << " on laser " << n_laser);
            }
            if ( ((tmpLaser.angle!=0) || (tmpLaser.isFocused)) && (tmpLaser.transv_profile!="focused") ) {
                WARNING("Laser "<<n_laser<<" transv_profile redefined as focused (Gaussian) and delta = "<<tmpLaser.delta<< " ignored");
                tmpLaser.transv_profile = "focused";
            }
        }//test on laser focus/angle
        
        laser_param.push_back(tmpLaser);
        n_laser++;
    }
    
	
}

