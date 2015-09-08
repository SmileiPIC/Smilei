#ifndef Profile_H
#define Profile_H

#include <cmath>
#include <vector>
#include "PyTools.h"
#include "SmileiMPI.h"
#include "Tools.h"


// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each Profile
// ---------------------------------------------------------------------------------------------------------------------
struct ProfileStructure {
    
    //! Magnitude of the profile if constant profile
    double profile;
    
    //! in case profile is give in Python
    PyObject *py_profile;
    
};

struct ExtFieldStructure;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Profile
//  --------------------------------------------------------------------------------------------------------------------
class Profile
{
public:
    //! Default constructor (for species profiles)
    Profile(ProfileStructure& , int);
    
    //! Alternate constructor (for external fields profiles)
    Profile(ExtFieldStructure&, int);
    
    //! Default destructor
    ~Profile(){};
    
    //! Some initialization
    void init();
        
    //! Function to get the value of the profile at some location
    double valueAt(std::vector<double>);
    
private:
    int nvariables;
    
    double (*Evaluate)(PyObject *, std::vector<double>);
    
protected:
    ProfileStructure  profile_param;
        
};//END class

#endif
