
#ifndef Laser_H
#define Laser_H

#include "PyTools.h"
#include "Profile.h"
#include "Field.h"
#include "Field1D.h"
#include "Field2D.h"

#include <vector>
#include <string>
#include <cmath>

class Params;
class Patch;


// Class for choosing specific profiles
class LaserProfile {
public:
    LaserProfile() {};
    ~LaserProfile() {};
    virtual double getAmplitude(std::vector<double> pos, double t, int j) {return 0.;};
    virtual std::string getInfo() { return "?"; };
};


//  --------------------------------------------------------------------------------------------------------------------
// Class for holding all information about one laser
//  --------------------------------------------------------------------------------------------------------------------
class Laser {
public:
    //! Normal laser constructor
    Laser(Params &params, int ilaser, Patch* patch);
    //! Cloning laser constructor
    Laser(Laser*);
    ~Laser();
    void clean();
    
    //! Gets the amplitude from both time and space profiles (By)
    inline double getAmplitude0(std::vector<double> pos, double t, int j) {
        return profiles[0]->getAmplitude(pos, t, j);
    }
    //! Gets the amplitude from both time and space profiles (Bz)
    inline double getAmplitude1(std::vector<double> pos, double t, int j) {
        return profiles[1]->getAmplitude(pos, t, j);
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
    LaserProfileSeparable(double, Profile*, Profile*, Profile*, Profile*, Params&, Patch*, bool);
    
    double getAmplitude(std::vector<double> pos, double t, int j);
private:
    double omega;
    Profile *timeProfile, *chirpProfile;
    Field *space_envelope, *phase;
};

// Laser profile for non-separable space and time
class LaserProfileNonSeparable : public LaserProfile {
public:
    LaserProfileNonSeparable(Profile * spaceAndTimeProfile)
     : spaceAndTimeProfile(spaceAndTimeProfile) {};
    
    inline double getAmplitude(std::vector<double> pos, double t, int j) {
        return spaceAndTimeProfile->valueAt(pos, t);
    }
private:
    Profile * spaceAndTimeProfile;
};

// Null laser profile
class LaserProfileNULL : public LaserProfile {
public:
    LaserProfileNULL() {};
    
    inline double getAmplitude(std::vector<double> pos, double t, int j) {
        return 0.;
    }
};


#endif
