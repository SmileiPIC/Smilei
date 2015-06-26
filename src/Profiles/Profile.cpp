#include "Profile.h"


using namespace std;

// Default constructor.
// Applies to profiles for species (density, velocity and temperature profiles)
/*MG150609
Profile::Profile(ProfileStructure & pp, string geometry, double convfac) :
conv_fac(convfac),
profile_param(pp),
factor(1.)
{
    // Launch the initialization
    init(profile_param, geometry);
    
}*/
Profile::Profile(ProfileStructure & pp, string geometry) :
profile_param(pp)
{
    
    // Launch the initialization
    init(profile_param, geometry);
    
}


// Special constructor.
// Applies to external field profiles
/*MG150609 Profile::Profile(ExtFieldStructure & pp, string geometry, double convfac) :
conv_fac(convfac)
{
    profile_param = static_cast<ProfileStructure> (pp);
    
    // define the factor
    factor = pp.magnitude;
    // define the vacuum_length as zeros
    profile_param.vacuum_length.resize(dim);
    for( int i=0; i<dim; i++) profile_param.vacuum_length[i]=0.;
    
    // Launch the initialization
    init(profile_param, geometry);
    
}*/
Profile::Profile(ExtFieldStructure & pp, string geometry)
{
    // Convert ExtFieldStructure in ProfileStructure
    profile_param = static_cast<ProfileStructure> (pp);
    
    // Launch the initialization
    init(profile_param, geometry);
    
}


void Profile::init(ProfileStructure & pp, string geometry)
{

    if      (geometry == "1d3v") dim = 1;
    else if (geometry == "2d3v") dim = 2;
    else {
        ERROR( "Unsupported geometry : " << geometry);
    }

}


Profile::~Profile()
{
}


double Profile::valueAt (vector<double> x_cell) {
    
    if        ( dim == 1 ) {
        //*MG150609 return PyTools::runPyFunction(profile_param.py_profile, x_cell[0]/conv_fac);
        return PyTools::runPyFunction(profile_param.py_profile, x_cell[0]);
    } else if ( dim == 2 ) {
        //*MG150609 return PyTools::runPyFunction(profile_param.py_profile, x_cell[0]/conv_fac, x_cell[1]/conv_fac);
        return PyTools::runPyFunction(profile_param.py_profile, x_cell[0], x_cell[1]);
    }
    
    return 0;
};
