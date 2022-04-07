#include "Tools.h"
#include "Params.h"
#include "Laser.h"
#include "Patch.h"
#include "H5.h"

#include <cmath>
#include <string>

using namespace std;


Laser::Laser( Params &params, int ilaser, Patch *patch, bool verbose )
{
    ostringstream name( "" );
    name << "Laser #" << ilaser;
    string errorPrefix = name.str();
    ostringstream info( "" );
    
    // side from which the laser enters the simulation box
    string box_side;
    PyTools::extract( "box_side", box_side, "Laser", ilaser );
    unsigned int normal_axis = 0;
    if( box_side == "xmin" || box_side == "xmax" ) {
        i_boundary_ = 0;
        normal_axis = 0;
    } else if( params.geometry == "1Dcartesian" ) {
        ERROR_NAMELIST( errorPrefix << ": in 1Dcartesian geometry, box_side must be `xmin` or `xmax`",
        LINK_NAMELIST + std::string("#laser"));
    } else if( params.geometry == "AMcylindrical" ) {
        ERROR_NAMELIST( errorPrefix << ": in AMcylindrical geometry, box_side must be `xmin` or `xmax`",
        LINK_NAMELIST + std::string("#laser"));
    } else if( box_side == "ymin" || box_side == "ymax" ) {
        i_boundary_ = 2;
        normal_axis = 1;
    } else if( params.geometry == "2Dcartesian" ) {
        ERROR_NAMELIST( errorPrefix << ": in 2Dcartesian geometry, box_side must be `xmin`, `xmax`, `ymin` or `ymax`",
        LINK_NAMELIST + std::string("#laser"));
    } else if( box_side == "zmin" || box_side == "zmax" ) {
        i_boundary_ = 4;
        normal_axis = 2;
    } else {
        ERROR_NAMELIST( errorPrefix << ": box_side must be `xmin`, `xmax`, `ymin`, `ymax`, `zmin` or `zmax`",
        LINK_NAMELIST + std::string("#laser"));
    }
    
    if( box_side.substr(1) == "max" ) {
        i_boundary_ ++;
    }
    
    // Profiles
    profiles.resize( 0 );
    PyObject *chirp_profile=nullptr, *time_profile=nullptr;
    vector<PyObject *>  space_profile, phase_profile, space_time_profile;
    double omega( 0 );
    file = "";
    PyTools::extract( "omega", omega, "Laser", ilaser );
    bool has_chirp         = PyTools::extract_pyProfile( "chirp_profile", chirp_profile, "Laser", ilaser );
    bool has_time          = PyTools::extract_pyProfile( "time_envelope", time_profile, "Laser", ilaser );
    bool has_space         = PyTools::extract2Profiles( "space_envelope", ilaser, space_profile );
    bool has_phase         = PyTools::extract2Profiles( "phase", ilaser, phase_profile );
    bool has_space_time    = PyTools::extract2Profiles( "space_time_profile", ilaser, space_time_profile );
    bool has_space_time_AM = PyTools::extract2NProfiles( "space_time_profile_AM", ilaser, space_time_profile );
    bool has_file          = PyTools::extractOrNone( "file", file, "Laser", ilaser );

    if( (has_space_time ||  has_space_time_AM) && has_file ) {
        ERROR_NAMELIST( errorPrefix << ": `space_time_profile` and `file` cannot both be set",
        LINK_NAMELIST + std::string("#laser"));
    }

    if(( has_space_time_AM) && params.geometry!="AMcylindrical"){
        ERROR_NAMELIST( errorPrefix << ": AM profiles can only be used in `AMcylindrical` geometry",
        LINK_NAMELIST + std::string("#laser"));
    }

    unsigned int space_dims     = ( params.geometry=="3Dcartesian" ? 2 : 1 );
    unsigned int spacetime_size = ( has_space_time_AM ? 2*params.nmodes+1 : 2 );//+1 to force spacetime_size to be always >2 in AM geometry.

    spacetime.resize( spacetime_size, false ); //Needs to be resized even if non separable profiles are not used.

    if( has_space_time || has_space_time_AM ) {

        info << "\t\t" << errorPrefix << ": space-time profile " << endl;
        if( has_time || has_space || has_chirp || has_phase ) {
            name.str( "" );
            name << ( has_time ?"time_envelope ":"" )
                 << ( has_space?"space_envelope ":"" )
                 << ( "omega " )
                 << ( has_chirp?"chirp_profile ":"" )
                 << ( has_phase?"phase ":"" );
            WARNING( errorPrefix << ": space-time profile defined, dismissing " << name.str() );
        }

        if(has_space_time_AM && space_time_profile.size() < 2*params.nmodes ) {
            WARNING( errorPrefix << ": not all modes are specified in the namelist. Unfilled higher order modes are considered null. " );
        }

        for (unsigned int i=0; i<space_time_profile.size(); i++){
            spacetime[i] = ( bool )( space_time_profile[i] );
        }

        for (unsigned int imode=0; imode<spacetime_size/2; imode++){ // only imode=0 if not AM
            // First axis (By or Br)
            name.str( "" );
            name << "Laser[" << ilaser <<"].space_time_profile["<< 2*imode << "]";
            if( spacetime[2*imode] ) {
                Profile *p = new Profile( space_time_profile[2*imode], params.nDim_field, name.str(), params );
                profiles.push_back( new LaserProfileNonSeparable( p ) );
                info << "\t\t\tfirst  component : " << p->getInfo();
                if (has_space_time_AM) info << " mode " << imode ;
                info << endl;
            } else {
                profiles.push_back( new LaserProfileNULL() );
                info << "\t\t\tfirst  component : zero" ;
                if (has_space_time_AM) info << " mode " << imode ;
                info << endl;
            }
            // Second axis (Bz or Bt)
            name.str( "" );
            name << "Laser[" << ilaser <<"].space_time_profile[" << 2*imode+1 << "]";
            if( spacetime[2*imode+1] ) {
                Profile *p = new Profile( space_time_profile[2*imode+1], params.nDim_field, name.str(), params );
                profiles.push_back( new LaserProfileNonSeparable( p ) );
                info << "\t\t\tsecond component : " << p->getInfo() ;
                if (has_space_time_AM) info << " mode " << imode ;
                info << endl;
            } else {
                profiles.push_back( new LaserProfileNULL() );
                info << "\t\t\tsecond component : zero" ;
                if (has_space_time_AM) info << " mode " << imode ;
                info << endl;
            }
        }


    } else if( has_file ) {

        info << "\t\t" << errorPrefix << endl;
        info << "\t\t\tData in file : " << file << endl;
        
        Profile *ptime1 = nullptr;
        Profile *ptime2 = nullptr;
        if( PyTools::extract_pyProfile( "_extra_envelope", time_profile, "Laser", ilaser ) ) {
            // extra envelope
            name.str( "" );
            name << "Laser[" << ilaser <<"].extra_envelope";
            ptime1 = new Profile( time_profile, space_dims+1, name.str(), params );
            ptime2 = new Profile( time_profile, space_dims+1, name.str(), params );
            info << "\t\t\tExtra envelope: " << ptime1->getInfo();
        } else {
            ERROR_NAMELIST( errorPrefix << ": `extra_envelope` missing or not understood",
            LINK_NAMELIST + std::string("#laser"));
        }

        profiles.push_back( new LaserProfileFile( file, ptime1, true , normal_axis ) );
        profiles.push_back( new LaserProfileFile( file, ptime2, false, normal_axis ) );

    } else {

        if( ! has_time ) {
            ERROR_NAMELIST( errorPrefix << ": missing `time_envelope`",
            LINK_NAMELIST + std::string("#laser"));
        }
        if( ! has_space ) {
            ERROR_NAMELIST( errorPrefix << ": missing `space_envelope`",
            LINK_NAMELIST + std::string("#laser"));
        }
        if( ! has_chirp ) {
            ERROR_NAMELIST( errorPrefix << ": missing `chirp_profile`",
            LINK_NAMELIST + std::string("#laser") );
        }
        if( ! has_phase ) {
            ERROR_NAMELIST( errorPrefix << ": missing `phase`",
            LINK_NAMELIST + std::string("#laser"));
        }

        info << "\t\t" << errorPrefix << ": separable profile" << endl;

        // omega
        info << "\t\t\tomega              : " << omega << endl;
        
        // chirp
        name.str( "" );
        name << "Laser[" << ilaser <<"].chirp_profile";
        Profile *pchirp1 = new Profile( chirp_profile, 1, name.str(), params );
        Profile *pchirp2 = new Profile( chirp_profile, 1, name.str(), params );
        info << "\t\t\tchirp_profile      : " << pchirp1->getInfo();

        // time envelope
        name.str( "" );
        name << "Laser[" << ilaser <<"].time_envelope";
        Profile *ptime1 = new Profile( time_profile, 1, name.str(), params );
        Profile *ptime2 = new Profile( time_profile, 1, name.str(), params );
        info << endl << "\t\t\ttime envelope      : " << ptime1->getInfo();

        // space envelope (By)
        name.str( "" );
        name << "Laser[" << ilaser <<"].space_envelope[0]";
        Profile *pspace1 = new Profile( space_profile[0], space_dims, name .str(), params );
        info << endl << "\t\t\tspace envelope (y) : " << pspace1->getInfo();

        // space envelope (Bz)
        name.str( "" );
        name << "Laser[" << ilaser <<"].space_envelope[1]";
        Profile *pspace2 = new Profile( space_profile[1], space_dims, name .str(), params );
        info << endl << "\t\t\tspace envelope (z) : " << pspace2->getInfo();

        // phase (By)
        name.str( "" );
        name << "Laser[" << ilaser <<"].phase[0]";
        Profile *pphase1 = new Profile( phase_profile[0], space_dims, name.str(), params );
        info << endl << "\t\t\tphase          (y) : " << pphase1->getInfo();

        // phase (Bz)
        name.str( "" );
        name << "Laser[" << ilaser <<"].phase[1]";
        Profile *pphase2 = new Profile( phase_profile[1], space_dims, name.str(), params );
        info << endl << "\t\t\tphase          (z) : " << pphase2->getInfo();

        // delay phase
        vector<double> delay_phase( 2, 0. );
        PyTools::extractV( "delay_phase", delay_phase, "Laser", ilaser );
        info << endl << "\t\tdelay phase      (y) : " << delay_phase[0];
        info << endl << "\t\tdelay phase      (z) : " << delay_phase[1];

        // Create the LaserProfiles
        profiles.push_back( new LaserProfileSeparable( omega, pchirp1, ptime1, pspace1, pphase1, delay_phase[0], true , normal_axis ) );
        profiles.push_back( new LaserProfileSeparable( omega, pchirp2, ptime2, pspace2, pphase2, delay_phase[1], false, normal_axis ) );

    }

    // Display info
    if( patch->isMaster() && verbose ) {
        MESSAGE( info.str() );
    }
}


// Cloning constructor
Laser::Laser( Laser *laser, Params &params )
{
    i_boundary_  = laser->i_boundary_;
    spacetime = laser->spacetime;
    file      = laser->file;
    bool spacetime_defined=false;

    for (unsigned int i=0;i<spacetime.size();i++){
        if(spacetime[i]) spacetime_defined=true;
    }

    profiles.resize( 0 );
    if( spacetime_defined ) {
        for (unsigned int i=0;i<spacetime.size();i++){
            if( spacetime[i] ) {
                profiles.push_back( new LaserProfileNonSeparable( static_cast<LaserProfileNonSeparable *>( laser->profiles[i] ) ) );
            } else {
                profiles.push_back( new LaserProfileNULL() );
            }
        }
    } else if( file != "" ) {
        profiles.push_back( new LaserProfileFile( static_cast<LaserProfileFile *>( laser->profiles[0] ) ) );
        profiles.push_back( new LaserProfileFile( static_cast<LaserProfileFile *>( laser->profiles[1] ) ) );
    } else {
        profiles.push_back( new LaserProfileSeparable( static_cast<LaserProfileSeparable *>( laser->profiles[0] ) ) );
        profiles.push_back( new LaserProfileSeparable( static_cast<LaserProfileSeparable *>( laser->profiles[1] ) ) );
    }
}


Laser::~Laser()
{
    for (unsigned int i=0;i<profiles.size();i++){
        delete profiles[i];
    }
}

void Laser::disable()
{
    for (unsigned int i=0;i<profiles.size();i++){
        profiles[i] = new LaserProfileNULL();
    }
}


// Separable laser profile constructor
LaserProfileSeparable::LaserProfileSeparable(
    double omega, Profile *chirpProfile, Profile *timeProfile,
    Profile *spaceProfile, Profile *phaseProfile, double delay_phase, bool primal, unsigned int axis
):
    primal_( primal ),
    omega_( omega ),
    timeProfile_( timeProfile ),
    chirpProfile_( chirpProfile ),
    spaceProfile_( spaceProfile ),
    phaseProfile_( phaseProfile ),
    delay_phase_( delay_phase ),
    axis_( axis )
{
    space_envelope = NULL;
    phase = NULL;
}
// Separable laser profile cloning constructor
LaserProfileSeparable::LaserProfileSeparable( LaserProfileSeparable *lp ) :
    primal_( lp->primal_ ),
    omega_( lp->omega_ ),
    timeProfile_( new Profile( lp->timeProfile_ ) ),
    chirpProfile_( new Profile( lp->chirpProfile_ ) ),
    spaceProfile_( new Profile( lp->spaceProfile_ ) ),
    phaseProfile_( new Profile( lp->phaseProfile_ ) ),
    delay_phase_( lp->delay_phase_ ),
    axis_( lp->axis_ )
{
    space_envelope = NULL;
    phase = NULL;
}
// Separable laser profile destructor
LaserProfileSeparable::~LaserProfileSeparable()
{
    if( timeProfile_ ) {
        delete timeProfile_;
    }
    if( chirpProfile_ ) {
        delete chirpProfile_;
    }

    if( spaceProfile_ ) {
        delete spaceProfile_;
    }
    if( phaseProfile_ ) {
        delete phaseProfile_;
    }
    if( space_envelope ) {
        delete space_envelope;
    }
    if( phase ) {
        delete phase;
    }
}


void LaserProfileSeparable::createFields( Params &params, Patch *patch )
{
    // Region size for SDMD
    std::vector<unsigned int> n_space( params.n_space );
    std::vector<unsigned int> oversize( params.oversize );
    if( params.multiple_decomposition && patch->vecSpecies.empty() ) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }
    
    vector<unsigned int> dim = { 1, 1 };
    
    if( params.geometry!="1Dcartesian" && params.geometry!="2Dcartesian" && params.geometry!="3Dcartesian" && params.geometry!="AMcylindrical" ) {
        ERROR_NAMELIST( "Unknown geometry in laser",
        LINK_NAMELIST + std::string("#laser"));
    }
    
    // Size in first direction
    if( params.geometry=="2Dcartesian" || params.geometry=="3Dcartesian" ) {
        unsigned int ax1 = ( axis_ == 0 ) ? 1 : 0;
        unsigned int n_p = n_space[ax1] + 1 + 2*oversize[ax1];
        unsigned int n_d = n_p + 1;
        dim[0] = primal_ ? n_p : n_d;
    } else if( params.geometry=="AMcylindrical" ) {
        unsigned int nr_p = n_space[1] + 1 + 2*oversize[1];
        unsigned int nr_d = nr_p + 1;
        dim[0] = nr_p + nr_d;
    }
    
    // Size in second direction
    if( params.geometry=="3Dcartesian" ) {
        unsigned int ax2 = ( axis_ == 2 ) ? 1 : 2;
        unsigned int n_p = n_space[ax2] + 1 + 2*oversize[ax2];
        unsigned int n_d = n_p + 1;
        dim[1] = primal_ ? n_d : n_p;
    }
    
    //Create laser fields
    space_envelope = new Field2D( dim );
    phase          = new Field2D( dim );
}

void LaserProfileSeparable::initFields( Params &params, Patch *patch )
{
    // Region size for SDMD
    std::vector<unsigned int> n_space(params.n_space);
    std::vector<unsigned int> oversize(params.oversize);
    if( params.multiple_decomposition && patch->vecSpecies.empty() ) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }
    
    if( params.geometry=="1Dcartesian" ) {
        
        // Assign profile (only one point in 1D)
        vector<double> pos( 1 );
        pos[0] = 0.;
        ( *space_envelope )( 0, 0 ) = spaceProfile_->valueAt( pos );
        ( *phase )( 0, 0 ) = phaseProfile_->valueAt( pos );
        
    } else if( params.geometry=="2Dcartesian" ) {
        
        unsigned int ax1 = ( axis_ == 0 ) ? 1 : 0;
        unsigned int n_p = n_space[ax1] + 1 + 2*oversize[ax1];
        unsigned int n_d = n_p + 1;
        double d = params.cell_length[ax1];
        unsigned int dim = primal_ ? n_p : n_d;
        
        // Assign profile
        vector<double> pos( 1 );
        pos[0] = patch->getDomainLocalMin( ax1 ) - ( ( primal_?0.:0.5 ) + oversize[ax1] )*d;
        for( unsigned int j=0 ; j<dim ; j++ ) {
            ( *space_envelope )( j, 0 ) = spaceProfile_->valueAt( pos );
            ( *phase )( j, 0 ) = phaseProfile_->valueAt( pos );
            pos[0] += d;
        }
        
    } else if( params.geometry=="AMcylindrical" ) {
        
        unsigned int nr_p = n_space[1]+1+2*oversize[1];
        unsigned int nr_d = nr_p+1;
        double dr = params.cell_length[1];
        unsigned int dim = nr_p + nr_d; // Need to account for both primal and dual positions
        
        // Assign profile
        vector<double> pos( 1 );
        for( unsigned int j=0 ; j<dim ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( j*0.5 - 0.5 - oversize[1] )*dr ; // Increment half cells
            ( *space_envelope )( j, 0 ) = spaceProfile_->valueAt( pos );
            ( *phase )( j, 0 ) = phaseProfile_->valueAt( pos );
        }
        
    } else if( params.geometry=="3Dcartesian" ) {
        
        unsigned int ax1 = ( axis_ == 0 ) ? 1 : 0;
        unsigned int ax2 = ( axis_ == 2 ) ? 1 : 2;
        unsigned int n1_p = n_space[ax1] + 1 + 2*oversize[ax1];
        unsigned int n1_d = n1_p + 1;
        unsigned int n2_p = n_space[ax2] + 1 + 2*oversize[ax2];
        unsigned int n2_d = n2_p + 1;
        double d1 = params.cell_length[ax1];
        double d2 = params.cell_length[ax2];
        unsigned int dim1 = primal_ ? n1_p : n1_d;
        unsigned int dim2 = primal_ ? n2_d : n2_p;
        
        // Assign profile
        vector<double> pos( 2 );
        pos[0] = patch->getDomainLocalMin( ax1 ) - ( ( primal_?0.:0.5 ) + oversize[ax1] )*d1;
        for( unsigned int j=0 ; j<dim1 ; j++ ) {
            pos[1] = patch->getDomainLocalMin( ax2 ) - ( ( primal_?0.5:0. ) + oversize[ax2] )*d2;
            for( unsigned int k=0 ; k<dim2 ; k++ ) {
                ( *space_envelope )( j, k ) = spaceProfile_->valueAt( pos );
                ( *phase )( j, k ) = phaseProfile_->valueAt( pos );
                pos[1] += d2;
            }
            pos[0] += d1;
        }
    }
}

// Amplitude of a separable laser profile
double LaserProfileSeparable::getAmplitude( std::vector<double> pos, double t, int j, int k )
{
    double amp;
    #pragma omp critical
    {
        double omega = omega_ * chirpProfile_->valueAt( t );
        double phi = ( *phase )( j, k );
        amp = timeProfile_->valueAt( t-( phi+delay_phase_ )/omega ) * ( *space_envelope )( j, k ) * sin( omega*t - phi );
    }
    return amp;
}

//Destructor
LaserProfileNonSeparable::~LaserProfileNonSeparable()
{
    if( spaceAndTimeProfile_ ) {
        delete spaceAndTimeProfile_;
    }
}


void LaserProfileFile::createFields( Params &params, Patch *patch )
{
    if( params.geometry!="2Dcartesian" && params.geometry!="3Dcartesian" ) {
        ERROR_NAMELIST( "Unknown geometry in LaserOffset (cartesian 2D or 3D only)",
        LINK_NAMELIST + std::string("#laser") );
    }
    
    magnitude = new Field3D();
    phase     = new Field3D();
}

void LaserProfileFile::initFields( Params &params, Patch *patch )
{
    // This is handled in checkpoints when restarting
    if( params.restart ) {
        return;
    }
    
    // Region size for SDMD
    std::vector<unsigned int> n_space(params.n_space);
    std::vector<unsigned int> oversize(params.oversize);
    if( params.multiple_decomposition && patch->vecSpecies.empty() ) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }
    
    unsigned int ndim = 2;
    if( params.geometry=="3Dcartesian" ) {
        ndim = 3;
    }
    
    // Define the part of the array to obtain
    unsigned int ax1 = ( axis_ == 0 ) ? 1 : 0;
    vector<hsize_t> dim( ndim ), offset( ndim );
    hsize_t n1_tot = params.n_space_global[ax1] + 2 + 2*params.oversize[ax1];
    hsize_t n1_p = n_space[ax1] + 1 + 2*oversize[ax1];
    hsize_t n1_d = n1_p + 1;
    dim[0] = primal_ ? n1_p : n1_d;
    dim[1] = 1;
    //offset[0] = patch->getCellStartingGlobalIndex( 1 ) + params.oversize[1];
    offset[1] = 0;
    
    if( ndim == 3 ) {
        unsigned int ax2 = ( axis_ == 2 ) ? 1 : 2;
        hsize_t n2_p = n_space[ax2] + 1 + 2*oversize[ax2];
        hsize_t n2_d = n2_p + 1;
        dim[1] = primal_ ? n2_d : n2_p;
        offset[1] = patch->getCellStartingGlobalIndex( ax2 ) + oversize[ax2];
    }
    
    // Open file
    H5Read f( file );
    // Obtain the omega dataset containing the different values of omega
    if( f.has( "omega" ) ) {
        f.vect( "omega", omega, true );
    } else {
        ERROR_NAMELIST( "File " << file << " does not contain the `omega` dataset",
        LINK_NAMELIST + std::string("#laser"));
    }
    
    // Allocate the fields
    magnitude->allocateDims( dim[0], dim[1], omega.size() );
    phase    ->allocateDims( dim[0], dim[1], omega.size() );
    
    // Read the datasets for the magnitude and phase of the field
    dim[ndim-1] = omega.size();
    string magnitude_name = primal_?"magnitude1":"magnitude2";
    if( f.has( magnitude_name ) ) {
        vector<hsize_t> shape = f.shape( magnitude_name );
        offset[0] = shape[0] - n1_tot + patch->getCellStartingGlobalIndex( ax1 ) + params.oversize[ax1];
        H5Space filespace( shape, offset, dim );
        H5Space memspace( dim );
        f.array( magnitude_name, magnitude->data_[0], &filespace, &memspace );
    } else {
        magnitude->put_to( 0. );
    }
    string phase_name = primal_?"phase1":"phase2";
    if( f.has( phase_name ) ) {
        vector<hsize_t> shape = f.shape( phase_name );
        offset[0] = shape[0] - n1_tot + patch->getCellStartingGlobalIndex( ax1 ) + params.oversize[ax1];
        H5Space filespace( shape, offset, dim );
        H5Space memspace( dim );
        f.array( phase_name, phase->data_[0], &filespace, &memspace );
    } else {
        phase->put_to( 0. );
    }
}

// Amplitude of a laser profile from a file (see LaserOffset)
double LaserProfileFile::getAmplitude( std::vector<double> pos, double t, int j, int k )
{
    double amp = 0;
    unsigned int n = omega.size();
    for( unsigned int i=0; i<n; i++ ) {
        amp += ( *magnitude )( j, k, i ) * cos( omega[i] * t + ( *phase )( j, k, i ) );
    }
    #pragma omp critical
    {
        amp *= extraProfile->valueAt( pos, t );
    }
    return amp;
}

//Destructor
LaserProfileFile::~LaserProfileFile()
{
    #pragma omp critical
    {
        if( magnitude )
        {
            delete magnitude   ;
            magnitude   =NULL;
        }
        if( phase )
        {
            delete phase       ;
            phase       =NULL;
        }
        if( extraProfile )
        {
            delete extraProfile;
            extraProfile=NULL;
        }
    }
}
