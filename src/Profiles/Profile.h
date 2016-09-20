#ifndef Profile_H
#define Profile_H

#include <vector>
#include <string>
#include "SmileiMPI.h"
#include "Tools.h"
#include "PyTools.h"
#include "Function.h"


class Profile
{
public:
    //! Default constructor
    Profile(PyObject*, unsigned int, std::string);
    //! Cloning constructor
    Profile(Profile*);
    //! Default destructor
    ~Profile();
    
    //! Get the value of the profile at some location (spatial)
    inline double valueAt(std::vector<double> coordinates) {
        return function->valueAt(coordinates);
    };
    //! Get the value of the profile at some location (temporal)
    inline double valueAt(double time) {
        return function->valueAt(time);
    };
    //! Get the value of the profile at some location (spatio-temporal)
    inline double valueAt(std::vector<double> coordinates, double time) {
        return function->valueAt(coordinates, time);
    };
    
    //! Get info on the loaded profile, to be printed later
    inline std::string getInfo() { return info; };
    
    //! Name of the profile, in the case of a built-in profile
    std::string profileName;
    
private:
    //! Object that holds the information on the profile function
    Function * function;
    
    //! String containing some info on the profile
    std::string info;
    
    //! Number of variables for the profile function
    int nvariables;
    
};//END class Profile




#endif
