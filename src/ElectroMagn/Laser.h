
#ifndef Laser_H
#define Laser_H

#include "PyTools.h"
#include "Profile.h"

#include <vector>
#include <string>

class Params;
class Patch;


// Class for choosing specific profiles
class LaserProfile {
public:
    LaserProfile() {};
    ~LaserProfile() {};
    virtual double getAmplitude(std::vector<double> pos, double t) {return 0.;};
    virtual std::string getInfo() { return "?"; };
};


//  --------------------------------------------------------------------------------------------------------------------
// Class for holding all information about one laser
//  --------------------------------------------------------------------------------------------------------------------
class Laser {
public:
    Laser(Params &params, int ilaser, Patch* patch);
    ~Laser();

    //! Gets the amplitude from both time and space profiles (By)
    inline double getAmplitude0(std::vector<double> pos, double t) {
        return profiles[0]->getAmplitude(pos, t);
    }
    //! Gets the amplitude from both time and space profiles (Bz)
    inline double getAmplitude1(std::vector<double> pos, double t) {
        return profiles[1]->getAmplitude(pos, t);
    }
    
    //! Side (west/east) from which the laser enters the box
    std::string boxSide;
    
    //! Disables the laser
    void disable();

private:
    //! Space and time profiles (Bx and By)
    std::vector<LaserProfile*> profiles;
    
};



// Laser profile for separable space and time
class LaserProfileSeparable : public LaserProfile {
public:
    LaserProfileSeparable(Profile * spaceProfile, Profile * timeProfile)
     : spaceProfile(spaceProfile), timeProfile(timeProfile) {};
    
    inline double getAmplitude(std::vector<double> pos, double t) {
        return timeProfile->valueAt(t) * spaceProfile->valueAt(pos);
    }
private:
    Profile * spaceProfile;
    Profile * timeProfile;
};

// Laser profile for non-separable space and time
class LaserProfileNonSeparable : public LaserProfile {
public:
    LaserProfileNonSeparable(Profile * spaceAndTimeProfile)
     : spaceAndTimeProfile(spaceAndTimeProfile) {};
    
    inline double getAmplitude(std::vector<double> pos, double t) {
        return spaceAndTimeProfile->valueAt(pos, t);
    }
private:
    Profile * spaceAndTimeProfile;
};

// Null laser profile
class LaserProfileNULL : public LaserProfile {
public:
    LaserProfileNULL() {};
    
    inline double getAmplitude(std::vector<double> pos, double t) {
        return 0.;
    }
};


#endif
