#ifndef Profile_H
#define Profile_H

#include "PyTools.h"
#include <cmath>
#include <vector>
#include "SmileiMPI.h"
#include "Tools.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class Profile
//  --------------------------------------------------------------------------------------------------------------------
class Profile
{
public:
    //! Default constructor
    Profile(PyObject* , unsigned int, std::string);

    //! Default destructor
    ~Profile(){};
    
    //! Function to get the value of the profile at some location
    double valueAt(std::vector<double>);
    
protected:
    PyObject *py_profile;
    
private:
    double (*Evaluate)(PyObject *, std::vector<double>);
    
};//END class

#endif
