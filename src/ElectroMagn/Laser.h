
#ifndef Laser_H
#define Laser_H

#include "PyTools.h"
#include "Profile.h"
#include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"

#include <vector>
#include <string>
#include <cmath>

class Params;
class Patch;


// Class for choosing specific profiles
class LaserProfile
{
    friend class SmileiMPI;
public:
    LaserProfile() {};
    virtual ~LaserProfile() {};
    virtual double getAmplitude( std::vector<double> pos, double t, int j, int k ) = 0;
    virtual std::complex<double> getAmplitudecomplex( std::vector<double>, double, int, int )
    {
        return 0.;
    };

    virtual std::string getInfo()
    {
        return "?";
    };
    virtual void createFields( Params &, Patch *, ElectroMagn * ) {};
    virtual void initFields( Params &, Patch *, ElectroMagn * ) {};
};


//  --------------------------------------------------------------------------------------------------------------------
// Class for holding all information about one laser
//  --------------------------------------------------------------------------------------------------------------------
class Laser
{
    friend class SmileiMPI;
    friend class Patch;
    friend class Checkpoint;
public:
    //! Normal laser constructor
    Laser( Params &params, int ilaser, Patch *patch, bool verbose = false );
    //! Cloning laser constructor
    Laser( Laser *, Params & );
    ~Laser();
    void clean();
    
    //! Gets the amplitude from both time and space profiles (By)
    inline double getAmplitude0( std::vector<double> pos, double t, int j, int k )
    {
        return profiles[0]->getAmplitude( pos, t, j, k );
    }
    //! Gets the amplitude from both time and space profiles (Bz)
    inline double getAmplitude1( std::vector<double> pos, double t, int j, int k )
    {
        return profiles[1]->getAmplitude( pos, t, j, k );
    }

    inline std::complex<double> getAmplitudecomplexN( std::vector<double> pos, double t, int j, int k, int imode )
    {
        return profiles[imode]->getAmplitudecomplex( pos, t, j, k );
    }
    
    void createFields( Params &params, Patch *patch, ElectroMagn *EMfields )
    {
        profiles[0]->createFields( params, patch, EMfields );
        profiles[1]->createFields( params, patch, EMfields );
    };
    void initFields( Params &params, Patch *patch, ElectroMagn *EMfields )
    {
        profiles[0]->initFields( params, patch, EMfields );
        profiles[1]->initFields( params, patch, EMfields );
    };
    
    //! Boundary from which the laser enters the box (xmin, xmax, ymin, ymax, zmin, zmax)
    unsigned int i_boundary_;
    
    //! Disables the laser
    void disable();

    //! True if spatio-temporal profile (Bx and By)
    std::vector<bool> spacetime;
    
protected:
    //! Space and time profiles (Bx and By)
    std::vector<LaserProfile *> profiles;
    
private:
    
    //! Non empty if laser profile read from a file
    std::string file;
    
};



// Laser profile for separable space and time
class LaserProfileSeparable : public LaserProfile
{
    friend class SmileiMPI;
    friend class Patch;
public:
    LaserProfileSeparable( double, Profile *, Profile *, Profile *, Profile *, double, bool, unsigned int );
    LaserProfileSeparable( LaserProfileSeparable * );
    ~LaserProfileSeparable();
    void createFields( Params &params, Patch *patch, ElectroMagn *EMfields ) override;
    void initFields( Params &params, Patch *patch, ElectroMagn *EMfields ) override;
    double getAmplitude( std::vector<double> pos, double t, int j, int k ) override;
protected:
    Field *space_envelope, *phase;
private:
    bool primal_;
    double omega_;
    Profile *timeProfile_, *chirpProfile_, *spaceProfile_, *phaseProfile_;
    double delay_phase_;
    unsigned int axis_;
};

// Laser profile for non-separable space and time
class LaserProfileNonSeparable : public LaserProfile
{
    friend class SmileiMPI;
public:
    LaserProfileNonSeparable( Profile *spaceAndTimeProfile )
        : spaceAndTimeProfile_( spaceAndTimeProfile ) {};
    LaserProfileNonSeparable( LaserProfileNonSeparable *lp )
        : spaceAndTimeProfile_( new Profile( lp->spaceAndTimeProfile_ ) ) {};
    ~LaserProfileNonSeparable();
    inline double getAmplitude( std::vector<double> pos, double t, int, int ) override
    {
        double amp;
        amp = spaceAndTimeProfile_->valueAt( pos, t );
        return amp;
    }

    inline std::complex<double> getAmplitudecomplex( std::vector<double> pos, double t, int, int ) override
    {
        std::complex<double> amp;
        amp = spaceAndTimeProfile_->complexValueAt( pos, t );
        return amp;
    }

private:
    Profile *spaceAndTimeProfile_;
};

// Laser profile from a file (see LaserOffset)
class LaserProfileFile : public LaserProfile
{
    friend class SmileiMPI;
    friend class Checkpoint;
public:
    LaserProfileFile( std::string file_, Profile *ep_, bool pr_, unsigned int axis )
        : magnitude( NULL ), phase( NULL ), file( file_ ), extraProfile( ep_ ), primal_( pr_ ), axis_( axis ) {};
    LaserProfileFile( LaserProfileFile *lp )
        : magnitude( NULL ), phase( NULL ), file( lp->file ), extraProfile( new Profile( lp->extraProfile ) ), primal_( lp->primal_ ), axis_( lp->axis_ ) {};
    ~LaserProfileFile();
    void createFields( Params &params, Patch *patch, ElectroMagn *EMfields ) override;
    void initFields( Params &params, Patch *patch, ElectroMagn *EMfields ) override;
    double getAmplitude( std::vector<double> pos, double t, int j, int k ) override;
protected:
    Field3D *magnitude, *phase;
    std::vector<double> omega;
private:
    std::string file;
    Profile *extraProfile;
    bool primal_;
    unsigned int axis_;
};

// Null laser profile
class LaserProfileNULL : public LaserProfile
{
public:
    LaserProfileNULL() {};
    ~LaserProfileNULL() {};
    
    inline double getAmplitude( std::vector<double>, double, int, int ) override
    {
        return 0.;
    }
};


#endif
