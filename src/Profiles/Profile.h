#ifndef Profile_H
#define Profile_H


#include <vector>
#include <string>
#include <complex>
#include "SmileiMPI.h"
#include "Tools.h"
#include "PyTools.h"
#include "H5.h"
#include "Function.h"
#include "Field.h"
#include "Field3D.h"
#include "cField.h"
#include <cmath>

class Profile
{
public:
    //! Default constructor
    Profile( PyObject *, unsigned int nvariables, std::string name, Params &params, bool try_numpy=false, bool try_file=false, bool time_variable=false );
    //! Cloning constructor
    Profile( Profile * );
    //! Default destructor
    ~Profile();
    
    //! Get the value of the profile at some location (spatial)
    inline double valueAt( std::vector<double> coordinates )
    {
        return function_->valueAt( coordinates );
    };
    //! Get the value of the profile at some location (temporal)
    inline double valueAt( double time )
    {
        return function_->valueAt( time );
    };
    //! Get the value of the profile at some location (spatio-temporal)
    inline double valueAt( std::vector<double> coordinates, double time )
    {
        return function_->valueAt( coordinates, time );
    };
    //! Get the complex value of the profile at some location (spatio-temporal)
    inline std::complex<double> complexValueAt( std::vector<double> coordinates, double time )
    {
        return function_->complexValueAt( coordinates, time );
    };
    
    //! Get/add the value of the profile at several locations
    //! mode = 0 : set values
    //! mode = 1 : ADD values
    //! mode = 2 : set values at given time
    //! mode = 3 : ADD values at given time
    void valuesAt( std::vector<Field *> &coordinates, std::vector<double> global_origin, Field &ret, int mode = 0, double time = 0. );
    
    //! Get/add the complex value of the profile at several locations
    //! mode = 0 : set values
    //! mode = 1 : ADD values
    //! mode = 2 : set values at given time
    //! mode = 3 : ADD values at given time
    void complexValuesAt( std::vector<Field *> &coordinates, cField &ret, int mode = 0, double time = 0. );
    
    //! Get the complex value of the profile at several locations (spatial + times)
    void complexValuesAtTimes( std::vector<Field *> &coordinates, Field *time, cField &ret );
    
    //! Get info on the loaded profile, to be printed later
    std::string getInfo()
    {
        std::ostringstream info( "" );
        info << nvariables_ << "D";
        
        if( ! profileName_.empty() ) {
            info << " built-in profile `" << profileName_ << "`" ;
        } else if( uses_file_ ) {
            info << " from file `" << filename_ << "`";
        } else {
            info << " user-defined function";
            if( uses_numpy_ ) {
                info << " (uses numpy)";
            }
        }
        
        if( function_ ) {
            info << function_->getInfo();
        }
        
        return info.str();
    };

    //! Get profile name
    std::string getProfileName()
    {
        return profileName_;
    }

private:
    
    //! Name of the profile, in the case of a built-in profile
    std::string profileName_;
    
    //! Object that holds the information on the profile function
    Function *function_;
    
    //! Number of variables for the profile function
    int nvariables_;
    
    //! Whether the profile is using numpy
    bool uses_numpy_;
    
    //! Whether the profile is taken from a file
    bool uses_file_;
    std::string filename_;
    
};//END class Profile




#endif
