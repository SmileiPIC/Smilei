#include "LaserParams.h"

#include <cmath>

#include "Tools.h"
#include "Params.h"

using namespace std;

LaserParams::LaserParams(Params& params) {
	
    // -----------------
    // Lasers properties
    // -----------------
    n_laser=PyTools::nComponents("Laser");
    for (unsigned int ilaser = 0; ilaser < n_laser; ilaser++) {
        LaserStructure tmpLaser;
        
        // side from which the laser enters the simulation box (only west/east at the moment)
        PyTools::extract("boxSide",tmpLaser.boxSide,"Laser",ilaser);
        if ( (tmpLaser.boxSide!="west") && (tmpLaser.boxSide!="east") ) {
            ERROR("At the moment laser can enter only from West/East sides: boxSide \""
                  << tmpLaser.boxSide << " not defined");
        }
        
        // laser intensity
        if( !PyTools::extract("a0",tmpLaser.a0,"Laser",ilaser)) {
            ERROR("Need parameter 'a0' for the laser");
        }
        
        // laser angular frequency (default=1)
        tmpLaser.omega0=1.0;
        PyTools::extract("omega0",tmpLaser.omega0,"Laser",ilaser);
        
        // laser temporal chirp (default=0)
        tmpLaser.tchirp=0.0;
        PyTools::extract("tchirp",tmpLaser.tchirp,"Laser",ilaser);
        
        // laser ellipticity/polarization parameter
        PyTools::extract("delta",tmpLaser.delta,"Laser",ilaser);
        
        // position of the laser focus
        tmpLaser.isFocused = PyTools::extract("focus",tmpLaser.focus,"Laser",ilaser);
        
        // incident angle
        tmpLaser.angle = 0.;
        PyTools::extract("angle",tmpLaser.angle ,"Laser",ilaser);
        
        // laser time-profile & associated parameters
        PyTools::extract("time_profile",tmpLaser.profile_time.profile ,"Laser",ilaser);
        if (tmpLaser.profile_time.profile.empty()) {
            PyObject *mypy = PyTools::extract_py("time_profile","Laser",ilaser);
            if (mypy && PyCallable_Check(mypy)) {
                tmpLaser.profile_time.py_profile=mypy;
                tmpLaser.profile_time.profile="python";
            } else {
                ERROR("Laser: time_profile parameter not understood");
            }
        } else {
            PyTools::extract("int_params",tmpLaser.profile_time.int_params ,"Laser",ilaser);
            PyTools::extract("double_params",tmpLaser.profile_time.double_params ,"Laser",ilaser);            
        }
        
        
        // laser transverse-profile & associated parameters
        if (params.geometry!="1d3v") {  //transv_profile is not define in 1d3v
            
            PyTools::extract("transv_profile",tmpLaser.profile_transv.profile ,"Laser",ilaser);
            if (tmpLaser.profile_transv.profile.empty()) {
                PyObject *mypy = PyTools::extract_py("transv_profile","Laser",ilaser);
                if (mypy && PyCallable_Check(mypy)) {
                    tmpLaser.profile_transv.py_profile=mypy;
                    tmpLaser.profile_transv.profile="python";
                } else {
                    ERROR("Laser: transv_profile not defined or not existing");
                }
            } else {
                PyTools::extract("int_params_transv",tmpLaser.profile_transv.int_params ,"Laser",ilaser);
                PyTools::extract("double_params_transv",tmpLaser.profile_transv.double_params ,"Laser",ilaser);
            }
        }//geometry
        
        //! \todo FP to MG : should we remove this now ? I commented it out because it is not used.
        //bool delayExists = ifile.extract("delay",tmpLaser.delay ,"Laser",ilaser);        
        PyTools::extract("delay",tmpLaser.delay ,"Laser",ilaser);
        
        // -----------------------------------------------------------------
        // tests on the laser parameters (when arbitrary focus or incidence)
        // -----------------------------------------------------------------
        // at the moment, only Gaussian beam can be used when laser is focused or with a not-normal incidence
        // note that there is also no choice about polarization (for along Bz)
        if ( (tmpLaser.boxSide=="west") && ((tmpLaser.angle!=0) || (tmpLaser.isFocused)) ){
            
            if ( (tmpLaser.profile_transv.int_params.size()==0) ) {
                WARNING("A default cut-off (3 sigma) is applied on laser " << ilaser << " transverse profile");
                tmpLaser.profile_transv.int_params.resize(1);
                tmpLaser.profile_transv.int_params[0] = 3;
            }
            
            if ( tmpLaser.delay==0. ) {
                tmpLaser.profile_transv.double_params.resize(1);
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
                    if (ylas_max > params.sim_length[1])
                        WARNING("Possible problem (simulation box size) with laser " << ilaser);
                    tmpLaser.delay = -ylas_max * sin(theta);
                } else {
                    double ylas_min = ylas - waist*sqrt(1.0+pow(tan(theta),2));
                    if (ylas_min < 0.0) WARNING("Possible problem (simulation box size) with laser " << ilaser);
                    tmpLaser.delay = -ylas_min * sin(theta);
                }
                // send a warning if delay is introduced
                if (tmpLaser.delay!=0)
                    WARNING("Introduction of a time-delay: " << tmpLaser.delay << " (in input units) on laser " << ilaser);
            }
            
            if ( ((tmpLaser.angle!=0) || (tmpLaser.isFocused)) && (tmpLaser.profile_transv.profile!="focused") ) {
                WARNING("Laser "<<ilaser<<" transv_profile redefined as focused (Gaussian) and delta = "<<tmpLaser.delta<< " ignored");
                tmpLaser.profile_transv.profile = "focused";
            }
        }//test on laser focus/angle
        
        laser_param.push_back(tmpLaser);
    }
    
    // -------------------------------------
    // Printing out laser related parameters
    // -------------------------------------
    MESSAGE("Laser related parameters");
    MESSAGE(1,"n_laser        : " << n_laser);
    for ( unsigned int i=0 ; i<n_laser ; i++ ) {
        MESSAGE(2,"laser " << i << ": (boxSide, a0) : (" << laser_param[i].boxSide <<  ", " << laser_param[i].a0 <<  ")");
    }
    
}

