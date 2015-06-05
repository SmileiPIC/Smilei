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
        
        // laser angular frequency (default=1)
        tmpLaser.omega0=1.0;
        ifile.extract("omega0",tmpLaser.omega0,"laser",0,n_laser);
        
        // laser temporal chirp (default=0)
        tmpLaser.tchirp=0.0;
        ifile.extract("tchirp",tmpLaser.tchirp,"laser",0,n_laser);
        
        // laser ellipticity/polarization parameter
        ifile.extract("delta",tmpLaser.delta,"laser",0,n_laser);
        
        // position of the laser focus
        tmpLaser.isFocused = ifile.extract("focus",tmpLaser.focus,"laser",0,n_laser);
        
        // incident angle
        ifile.extract("angle",tmpLaser.angle ,"laser",0,n_laser);
        
        // laser time-profile & associated parameters
        ifile.extract("time_profile",tmpLaser.profile_time.profile ,"laser",0,n_laser);
        ifile.extract("int_params",tmpLaser.profile_time.int_params ,"laser",0,n_laser);
        ifile.extract("double_params",tmpLaser.profile_time.double_params ,"laser",0,n_laser);
        
        // laser transverse-profile & associated parameters
        ifile.extract("transv_profile",tmpLaser.profile_transv.profile ,"laser",0,n_laser);
        ifile.extract("int_params_transv",tmpLaser.profile_transv.int_params ,"laser",0,n_laser);
        ifile.extract("double_params_transv",tmpLaser.profile_transv.double_params ,"laser",0,n_laser);
        
        bool delayExists = ifile.extract("delay",tmpLaser.delay ,"laser",0,n_laser);
        
        
        // -------------------------------------
        // Printing out laser related parameters
        // -------------------------------------
        if (n_laser==0) // just print "Laser related parameters" once
            MESSAGE("Laser related parameters");
        MESSAGE(1, "n_laser: " << n_laser << " : (boxSide, a0) : (" << tmpLaser.boxSide
                << ", " << tmpLaser.a0 <<  ")");
        //for ( unsigned int i=0 ; i<n_laser ; i++ ) {
        //    MESSAGE(2,"laser " << i << ": (boxSide, a0) : (" << laser_param[i].boxSide <<  ", " << laser_param[i].a0 <<  ")");
        //}
        
        
        // -----------------------------------------------------------------
        // normalization (from wavelength-related units to normalized units)
        // -----------------------------------------------------------------
        for (unsigned int i=0; i<tmpLaser.profile_time.double_params.size(); i++)
            tmpLaser.profile_time.double_params[i] *= params.conv_fac;
        for (unsigned int i=0; i<tmpLaser.profile_transv.double_params.size(); i++)
            tmpLaser.profile_transv.double_params[i] *= params.conv_fac;
        
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
            
            if ( (tmpLaser.profile_transv.int_params.size()==0) ) {
                WARNING("A default cut-off (3 sigma) is applied on laser " << n_laser << " transverse profile");
                tmpLaser.profile_transv.int_params.resize(1);
                tmpLaser.profile_transv.int_params[0] = 3;
            }
            
            if ( !delayExists ) {
                double theta   = tmpLaser.angle * M_PI/180.0;
                double xfoc    = tmpLaser.focus[0];
                double yfoc    = tmpLaser.focus[1];
                double ylas    = yfoc - tan(theta)*xfoc;
                // estimating the effect of diffraction at the border
                double bwaist  = 0.5/sqrt(log(2.0))*tmpLaser.profile_transv.double_params[0];
                double zeta    = tmpLaser.focus[0]/cos(theta);       // here, this could be improved (only approximated)
                double z2ovLr2 = pow(zeta,2)/pow(bwaist,4);
                double waist   = bwaist * sqrt(1.0+z2ovLr2) * tmpLaser.profile_transv.int_params[0];
                // computing the entering point of the laser & delay
                if (theta<0) {
                    double ylas_max = ylas + waist*sqrt(1.0+pow(tan(theta),2));
                    if (ylas_max > params.sim_length[1]) WARNING("Possible problem (simulation box size) with laser " << n_laser);
                    tmpLaser.delay = -ylas_max * sin(theta);
                } else {
                    double ylas_min = ylas - waist*sqrt(1.0+pow(tan(theta),2));
                    if (ylas_min < 0.0) WARNING("Possible problem (simulation box size) with laser " << n_laser);
                    tmpLaser.delay = -ylas_min * sin(theta);
                }
                // send a warning if delay is introduced
                if (tmpLaser.delay!=0)
                    WARNING("Introduction of a time-delay: " << tmpLaser.delay/params.conv_fac << " (in input units) on laser " << n_laser);
            }
            
            if ( ((tmpLaser.angle!=0) || (tmpLaser.isFocused)) && (tmpLaser.profile_transv.profile!="focused") ) {
                WARNING("Laser "<<n_laser<<" transv_profile redefined as focused (Gaussian) and delta = "<<tmpLaser.delta<< " ignored");
                tmpLaser.profile_transv.profile = "focused";
            }
        }//test on laser focus/angle
        
        laser_param.push_back(tmpLaser);
        n_laser++;
    }
    
	
}

