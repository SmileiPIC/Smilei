#include <cmath>

#include "Tools.h"
#include "Params.h"
#include "Laser.h"
#include "Patch.h"

using namespace std;


Laser::Laser(Params &params, int ilaser, Patch* patch)
{
    
    ostringstream name("");
    name << "Laser #" << ilaser;
    string errorPrefix = name.str();
    ostringstream info("");
    
    // side from which the laser enters the simulation box (only xmin/xmax at the moment)
    PyTools::extract("box_side",box_side,"Laser",ilaser);
    if ( box_side!="xmin" && box_side!="xmax" ) {
        ERROR(errorPrefix << ": box_side must be `xmin` or `xmax`");
    }
    
    // Profiles
    profiles.resize(0);
    PyObject *chirp_profile=nullptr, *time_profile=nullptr;
    vector<PyObject*>  space_profile, phase_profile, space_time_profile;
    bool time, space, omega, chirp, phase, space_time;
    double omega_value(0);
    Profile *p, *pchirp, *pchirp2, *ptime, *ptime2, *pspace1, *pspace2, *pphase1, *pphase2;
    pchirp2 = NULL;
    ptime2  = NULL;
    omega      = PyTools::extract("omega",omega_value,"Laser",ilaser);
    chirp      = PyTools::extract_pyProfile("chirp_profile"     , chirp_profile, "Laser", ilaser);
    time       = PyTools::extract_pyProfile("time_envelope"     , time_profile , "Laser", ilaser);
    space      = PyTools::extract2Profiles ("space_envelope"    , ilaser, space_profile     );
    phase      = PyTools::extract2Profiles ("phase"             , ilaser, phase_profile     );
    space_time = PyTools::extract2Profiles ("space_time_profile", ilaser, space_time_profile);
    
    spacetime.resize(2, false);
    if( space_time ) {
        
        spacetime[0] = (bool)(space_time_profile[0]);
        spacetime[1] = (bool)(space_time_profile[1]);
        
        if( time || space || omega || chirp || phase ) {
            name.str("");
            name << (time ?"time_envelope ":"")
                 << (space?"space_envelope ":"")
                 << (omega?"omega ":"")
                 << (chirp?"chirp_profile ":"")
                 << (phase?"phase ":"");
            WARNING(errorPrefix << ": space-time profile defined, dismissing " << name.str() );
        }
        
        info << "\t\t" << errorPrefix << ": space-time profile" << endl;
        
        // By
        name.str("");
        name << "Laser[" << ilaser <<"].space_time_profile[0]";
        if( spacetime[0] ) {
            p = new Profile(space_time_profile[0], params.nDim_field, name.str());
            profiles.push_back( new LaserProfileNonSeparable(p) );
            info << "\t\t\tfirst  axis : " << p->getInfo() << endl;
        } else {
            profiles.push_back( new LaserProfileNULL() );
            info << "\t\t\tfirst  axis : zero" << endl;
        }
        // Bz
        name.str("");
        name << "Laser[" << ilaser <<"].space_time_profile[1]";
        if( spacetime[1] ) {
            p = new Profile(space_time_profile[1], params.nDim_field, name.str());
            profiles.push_back( new LaserProfileNonSeparable(p) );
            info << "\t\t\tsecond axis : " << p->getInfo();
        } else {
            profiles.push_back( new LaserProfileNULL() );
            info << "\t\t\tsecond axis : zero";
        }
        
    } else {
        
        if( !time )
            ERROR(errorPrefix << ": missing `time_envelope`");
        if( !space )
            ERROR(errorPrefix << ": missing `space_envelope`");
        if( !omega )
            ERROR(errorPrefix << ": missing `omega`");
        if( !chirp )
            ERROR(errorPrefix << ": missing `chirp_profile`");
        if( !phase )
            ERROR(errorPrefix << ": missing `phase`");
        
        info << "\t\t" << errorPrefix << ": separable profile" << endl;
        
        unsigned int space_dims = (params.geometry=="3Dcartesian" ? 2 : 1);
        
        // omega
        info << "\t\t\tomega              : " << omega_value << endl;
        
        // chirp
        name.str("");
        name << "Laser[" << ilaser <<"].chirp_profile";
        pchirp = new Profile(chirp_profile, 1, name.str());
        pchirp2 = new Profile(chirp_profile, 1, name.str());
        info << "\t\t\tchirp_profile      : " << pchirp->getInfo();
        
        // time envelope
        name.str("");
        name << "Laser[" << ilaser <<"].time_envelope";
        ptime = new Profile(time_profile, 1, name.str());
        ptime2 = new Profile(time_profile, 1, name.str());
        info << endl << "\t\t\ttime envelope      : " << ptime->getInfo();
         
        // space envelope (By)
        name.str("");
        name << "Laser[" << ilaser <<"].space_envelope[0]";
        pspace1 = new Profile(space_profile[0], space_dims, name .str());
        info << endl << "\t\t\tspace envelope (y) : " << pspace1->getInfo();
        
        // space envelope (Bz)
        name.str("");
        name << "Laser[" << ilaser <<"].space_envelope[1]";
        pspace2 = new Profile(space_profile[1], space_dims, name .str());
        info << endl << "\t\t\tspace envelope (z) : " << pspace2->getInfo();
        
        // phase (By)
        name.str("");
        name << "Laser[" << ilaser <<"].phase[0]";
        pphase1 = new Profile(phase_profile[0], space_dims, name.str());
        info << endl << "\t\t\tphase          (y) : " << pphase1->getInfo();
        
        // phase (Bz)
        name.str("");
        name << "Laser[" << ilaser <<"].phase[1]";
        pphase2 = new Profile(phase_profile[1], space_dims, name.str());
        info << endl << "\t\t\tphase          (z) : " << pphase2->getInfo();
        
        // delay phase
        vector<double> delay_phase(2, 0.);
        PyTools::extract("delay_phase",delay_phase,"Laser",ilaser);
        info << endl << "\t\tdelay phase      (y) : " << delay_phase[0];
        info << endl << "\t\tdelay phase      (z) : " << delay_phase[1];
        
        // Create the LaserProfiles
        profiles.push_back( new LaserProfileSeparable(omega_value, pchirp , ptime , pspace1, pphase1, delay_phase[0], true ) );
        profiles.push_back( new LaserProfileSeparable(omega_value, pchirp2, ptime2, pspace2, pphase2, delay_phase[1], false) );
    
    }
    
    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
    }
}


// Cloning constructor
Laser::Laser(Laser* laser, Params& params)
{
    box_side   = laser->box_side;
    spacetime = laser->spacetime;
    profiles.resize(0);
    if( spacetime[0] || spacetime[1] ) {
        if( spacetime[0] ) {
            profiles.push_back( new LaserProfileNonSeparable(static_cast<LaserProfileNonSeparable*>(laser->profiles[0])) );
        } else {
            profiles.push_back( new LaserProfileNULL() );
        }
        if( spacetime[1] ) {
            profiles.push_back( new LaserProfileNonSeparable(static_cast<LaserProfileNonSeparable*>(laser->profiles[1])) );
        } else {
            profiles.push_back( new LaserProfileNULL() );
        }
    } else {
        profiles.push_back( new LaserProfileSeparable(static_cast<LaserProfileSeparable*>(laser->profiles[0])) );
        profiles.push_back( new LaserProfileSeparable(static_cast<LaserProfileSeparable*>(laser->profiles[1])) );
    }
}


Laser::~Laser()
{
    delete profiles[0];
    delete profiles[1];
}

void Laser::disable()
{
    
    profiles[0] = new LaserProfileNULL();
    profiles[1] = new LaserProfileNULL();
    
}


// Separable laser profile constructor
LaserProfileSeparable::LaserProfileSeparable(
    double omega, Profile* chirpProfile, Profile* timeProfile,
    Profile* spaceProfile, Profile* phaseProfile, double delay_phase, bool primal
):
    primal       ( primal       ),
    omega        ( omega        ),
    timeProfile  ( timeProfile  ),
    chirpProfile ( chirpProfile ),
    spaceProfile ( spaceProfile ),
    phaseProfile ( phaseProfile ),
    delay_phase  ( delay_phase  )
{
    space_envelope = NULL;
    phase = NULL;
}
// Cloning constructor
LaserProfileSeparable::LaserProfileSeparable(LaserProfileSeparable * lp) :
    primal       ( lp->primal ),
    omega        ( lp->omega  ),
    timeProfile  ( new Profile(lp->timeProfile ) ),
    chirpProfile ( new Profile(lp->chirpProfile) ),
    spaceProfile ( new Profile(lp->spaceProfile) ),
    phaseProfile ( new Profile(lp->phaseProfile) ),
    delay_phase  ( lp->delay_phase )
{
    space_envelope = NULL;
    phase = NULL;
}
//Destructor
LaserProfileSeparable::~LaserProfileSeparable()
{
    if(timeProfile   ) delete timeProfile;
    if(chirpProfile  ) delete chirpProfile;
        
    if(spaceProfile  ) delete spaceProfile;
    if(phaseProfile  ) delete phaseProfile;
    if(space_envelope) delete space_envelope;
    if(phase         ) delete phase;
}


void LaserProfileSeparable::createFields(Params& params, Patch* patch)
{
    vector<unsigned int> dim(2);
    dim[0] = 1;
    dim[1] = 1;
    
    if( params.geometry!="1Dcartesian" && params.geometry!="2Dcartesian" && params.geometry!="3Dcartesian" )
        ERROR("Unknown geometry in laser");
    
    if( params.geometry!="1Dcartesian" ) {
        unsigned int ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
        unsigned int ny_d = ny_p+1;
        dim[0] = primal ? ny_p : ny_d;
        
        if( params.geometry!="2Dcartesian" ) {
            unsigned int nz_p = params.n_space[2]*params.global_factor[2]+1+2*params.oversize[2];
            unsigned int nz_d = nz_p+1;
            dim[1] = primal ? nz_d : nz_p;
        }
    }
    
    space_envelope = new Field2D(dim);
    phase          = new Field2D(dim);
}

void LaserProfileSeparable::initFields(Params& params, Patch* patch)
{
    if( params.geometry=="1Dcartesian" ) {
        
        // Assign profile (only one point in 1D)
        vector<double> pos(1);
        pos[0] = 0.;
        (*space_envelope)(0,0) = spaceProfile->valueAt(pos);
        (*phase         )(0,0) = phaseProfile->valueAt(pos);
        
    } else if( params.geometry=="2Dcartesian" ) {
        
        unsigned int ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
        unsigned int ny_d = ny_p+1;
        double dy = params.cell_length[1];
        vector<unsigned int> dim(1);
        dim[0] = primal ? ny_p : ny_d;
        
        // Assign profile
        vector<double> pos(1);
        pos[0] = patch->getDomainLocalMin(1) - ((primal?0.:0.5) + params.oversize[1])*dy;
        for (unsigned int j=0 ; j<dim[0] ; j++) {
            (*space_envelope)(j,0) = spaceProfile->valueAt(pos);
            (*phase         )(j,0) = phaseProfile->valueAt(pos);
            pos[0] += dy;
        }
        
    } else if( params.geometry=="3Dcartesian" ) {
        
        unsigned int ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
        unsigned int ny_d = ny_p+1;
        unsigned int nz_p = params.n_space[2]*params.global_factor[2]+1+2*params.oversize[2];
        unsigned int nz_d = nz_p+1;
        double dy = params.cell_length[1];
        double dz = params.cell_length[2];
        vector<unsigned int> dim(2);
        dim[0] = primal ? ny_p : ny_d;
        dim[1] = primal ? nz_d : nz_p;
        
        // Assign profile
        vector<double> pos(2);
        pos[0] = patch->getDomainLocalMin(1) - ((primal?0.:0.5) + params.oversize[1])*dy;
        for (unsigned int j=0 ; j<dim[0] ; j++) {
            pos[1] = patch->getDomainLocalMin(2) - ((primal?0.5:0.) + params.oversize[2])*dz;
            for (unsigned int k=0 ; k<dim[1] ; k++) {
                (*space_envelope)(j,k) = spaceProfile->valueAt(pos);
                (*phase         )(j,k) = phaseProfile->valueAt(pos);
                pos[1] += dz;
            }
            pos[0] += dy;
        }
    }
}



// Amplitude of a separable laser profile
double LaserProfileSeparable::getAmplitude(std::vector<double> pos, double t, int j, int k)
{
    double amp;
    #pragma omp critical
    {
        double omega_ = omega * chirpProfile->valueAt(t);
        double phi = (*phase)(j, k);
        amp = timeProfile->valueAt(t-(phi+delay_phase)/omega_) * (*space_envelope)(j, k) * sin( omega_*t - phi );
    }
    return amp;
}

//Destructor
LaserProfileNonSeparable::~LaserProfileNonSeparable()
{
    if(spaceAndTimeProfile) delete spaceAndTimeProfile;
}


