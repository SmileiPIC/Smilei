#ifndef Profile_H
#define Profile_H

#include <cmath>
#include "PicParams.h"
#include "ExtFieldParams.h"
#include "PyTools.h"
#include "SmileiMPI.h"
#include "Tools.h"




//  --------------------------------------------------------------------------------------------------------------------
//! Class Profile
//  --------------------------------------------------------------------------------------------------------------------
class Profile
{
public:
    //! Default constructor (for species profiles)
    Profile(ProfileStructure& , std::string);
    
    //! Alternate constructor (for external fields profiles)
    Profile(ExtFieldStructure&, std::string);
    
    //! Default destructor
    ~Profile();
    
    //! Complementary function to all constructors
    void init(ProfileStructure& , std::string);
    
    //! Function to get the value of the profile at some location
    double valueAt(std::vector<double>);
    
private:
    int dim;
    
protected:
    ProfileStructure  profile_param;
        
};//END class

#endif
