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
    
    // side from which the laser enters the simulation box (only west/east at the moment)
    PyTools::extract("boxSide",boxSide,"Laser",ilaser);
    if ( boxSide!="west" && boxSide!="east" ) {
        ERROR(errorPrefix << ": boxSide must be `west` or `east`");
    }
    
    // Profiles
    profiles.resize(0);
    PyObject * space_profile;
    vector<PyObject*> time_profile, space_time_profile;
    bool time, space, space_time;
    Profile *p1, *p2;
    space      = PyTools::extract_pyProfile("space_profile"     , space_profile, "Laser", ilaser);
    time       = PyTools::extract2Profiles ("time_profile"      , ilaser, time_profile      );
    space_time = PyTools::extract2Profiles ("space_time_profile", ilaser, space_time_profile);
    
    if( space_time ) {
        
        if( time ) {
            ERROR(errorPrefix << ": cannot have both time and space-time profiles");
        }
        if( space ) {
            ERROR(errorPrefix << ": cannot have both space and space-time profiles");
        }
        
        info << "\t\t" << errorPrefix << ": space-time profile" << endl;
        
        // By
        name.str("");
        name << "Laser[" << ilaser <<"].space_time_profile[0]";
        if( space_time_profile[0] ) {
            p1 = new Profile(space_time_profile[0], params.nDim_field, name.str());
            profiles.push_back( new LaserProfileNonSeparable(p1) );
            info << "\t\t\tfirst  axis : " << p1->getInfo() << endl;
        } else {
            profiles.push_back( new LaserProfileNULL() );
            info << "\t\t\tfirst  axis : zero" << endl;
        }
        // Bz
        name.str("");
        name << "Laser[" << ilaser <<"].space_time_profile[1]";
        if( space_time_profile[1] ) {
            p1 = new Profile(space_time_profile[1], params.nDim_field, name.str());
            profiles.push_back( new LaserProfileNonSeparable(p1) );
            info << "\t\t\tsecond axis : " << p1->getInfo();
        } else {
            profiles.push_back( new LaserProfileNULL() );
            info << "\t\t\tsecond axis : zero";
        }
        
    } else {
        
        if( !time or !space ) {
            ERROR(errorPrefix << ": needs both time and space profiles");
        }
        
        info << "\t\t" << errorPrefix << ": separate space and time profiles" << endl;
        
        // space
        name.str("");
        name << "Laser[" << ilaser <<"].space_profile[0]";
        if( params.geometry=="3d3v" ) {
            p1 = new Profile(space_profile, 2, name .str());
        } else {
            p1 = new Profile(space_profile, 1, name .str());
        }
        info << "\t\t\tspace  : " << p1->getInfo() << endl;
        
        // By
        ostringstream name1("");
        name1 << "Laser[" << ilaser <<"].time_profile[0]";
        if( time_profile[0] ) {
            p2 = new Profile(time_profile[0], 1, name1.str());
            profiles.push_back( new LaserProfileSeparable( p1, p2 ) );
            info << "\t\t\ttime (1) : " << p2->getInfo() << endl;
        } else {
            profiles.push_back( new LaserProfileNULL() );
            info << "\t\t\ttime (1) : zero" << endl;
        }
        // Bz
        name1.str("");
        name1 << "Laser[" << ilaser <<"].time_profile[1]";
        if( time_profile[1] ) {
            p2 = new Profile(time_profile[1], 1, name1.str());
            profiles.push_back( new LaserProfileSeparable( p1, p2 ) );
            info << "\t\t\ttime (2) : " << p2->getInfo();
        } else {
            profiles.push_back( new LaserProfileNULL() );
            info << "\t\t\ttime (2) : zero";
        }
    
    }
    
    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
    }
}


Laser::~Laser()
{
    
    delete profiles[0];
    delete profiles[1];
    
}

void Laser::disable()
{
    
    delete profiles[0];
    delete profiles[1];
    profiles[0] = new LaserProfileNULL();
    profiles[1] = new LaserProfileNULL();
    
}
