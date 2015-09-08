#include "ElectroMagn.h"
#include "Profile.h"

using namespace std;

struct ExtFieldStructure;

// Default constructor.
// Applies to profiles for species (density, velocity and temperature profiles)
Profile::Profile(ProfileStructure & pp, string geometry) :
profile_param(pp)
{    
    if      (geometry == "1d3v") dim = 1;
    else if (geometry == "2d3v") dim = 2;
    else {
        ERROR( "Unsupported geometry : " << geometry);
    }
    
}


// Special constructor.
// Applies to external field profiles
Profile::Profile(ExtFieldStructure & pp, int ndim):
dim(ndim),
profile_param(static_cast<ProfileStructure> (pp))
{
}

double Profile::valueAt (vector<double> x_cell) {
    
    if        ( dim == 1 ) {
        return PyTools::runPyFunction(profile_param.py_profile, x_cell[0]);
    } else if ( dim == 2 ) {
        return PyTools::runPyFunction(profile_param.py_profile, x_cell[0], x_cell[1]);
    } else if ( dim == 3 ) {
        return PyTools::runPyFunction(profile_param.py_profile, x_cell[0], x_cell[1], x_cell[2]);
    }
    
    return 0;
};
