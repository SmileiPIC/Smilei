#include "Tools.h"
#include "Params.h"
#include "Laser.h"
#include "Patch.h"
#include "H5.h"

#include <cmath>
#include <string>

using namespace std;


Laser::Laser( Params &params, int ilaser, Patch *patch )
{
    ostringstream name( "" );
    name << "Laser #" << ilaser;
    string errorPrefix = name.str();
    ostringstream info( "" );

    // side from which the laser enters the simulation box (only xmin/xmax at the moment)
    PyTools::extract( "box_side", box_side, "Laser", ilaser );
    if( box_side!="xmin" && box_side!="xmax" ) {
        ERROR( errorPrefix << ": box_side must be `xmin` or `xmax`" );
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
        ERROR( errorPrefix << ": `space_time_profile` and `file` cannot both be set" );
    }

    if(( has_space_time_AM) && params.geometry!="AMcylindrical"){
        ERROR( errorPrefix << ": AM profiles can only be used in `AMcylindrical` geometry" );
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
                Profile *p = new Profile( space_time_profile[2*imode], params.nDim_field, name.str() );
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
                Profile *p = new Profile( space_time_profile[2*imode+1], params.nDim_field, name.str() );
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
        
        Profile *ptime, *ptime2;
        if( PyTools::extract_pyProfile( "_extra_envelope", time_profile, "Laser", ilaser ) ) {
            // extra envelope
            name.str( "" );
            name << "Laser[" << ilaser <<"].extra_envelope";
            ptime  = new Profile( time_profile, space_dims+1, name.str() );
            ptime2 = new Profile( time_profile, space_dims+1, name.str() );
            info << "\t\t\tExtra envelope: " << ptime->getInfo();
        } else {
            ERROR( errorPrefix << ": `extra_envelope` missing or not understood" );
        }

        profiles.push_back( new LaserProfileFile( file, ptime, true ) );
        profiles.push_back( new LaserProfileFile( file, ptime2, false ) );

    } else {

        if( ! has_time ) {
            ERROR( errorPrefix << ": missing `time_envelope`" );
        }
        if( ! has_space ) {
            ERROR( errorPrefix << ": missing `space_envelope`" );
        }
        if( ! has_chirp ) {
            ERROR( errorPrefix << ": missing `chirp_profile`" );
        }
        if( ! has_phase ) {
            ERROR( errorPrefix << ": missing `phase`" );
        }

        info << "\t\t" << errorPrefix << ": separable profile" << endl;

        // omega
        info << "\t\t\tomega              : " << omega << endl;
        
        // chirp
        name.str( "" );
        name << "Laser[" << ilaser <<"].chirp_profile";
        Profile *pchirp = new Profile( chirp_profile, 1, name.str() );
        Profile *pchirp2 = new Profile( chirp_profile, 1, name.str() );
        info << "\t\t\tchirp_profile      : " << pchirp->getInfo();

        // time envelope
        name.str( "" );
        name << "Laser[" << ilaser <<"].time_envelope";
        Profile *ptime = new Profile( time_profile, 1, name.str() );
        Profile *ptime2 = new Profile( time_profile, 1, name.str() );
        info << endl << "\t\t\ttime envelope      : " << ptime->getInfo();

        // space envelope (By)
        name.str( "" );
        name << "Laser[" << ilaser <<"].space_envelope[0]";
        Profile *pspace1 = new Profile( space_profile[0], space_dims, name .str() );
        info << endl << "\t\t\tspace envelope (y) : " << pspace1->getInfo();

        // space envelope (Bz)
        name.str( "" );
        name << "Laser[" << ilaser <<"].space_envelope[1]";
        Profile *pspace2 = new Profile( space_profile[1], space_dims, name .str() );
        info << endl << "\t\t\tspace envelope (z) : " << pspace2->getInfo();

        // phase (By)
        name.str( "" );
        name << "Laser[" << ilaser <<"].phase[0]";
        Profile *pphase1 = new Profile( phase_profile[0], space_dims, name.str() );
        info << endl << "\t\t\tphase          (y) : " << pphase1->getInfo();

        // phase (Bz)
        name.str( "" );
        name << "Laser[" << ilaser <<"].phase[1]";
        Profile *pphase2 = new Profile( phase_profile[1], space_dims, name.str() );
        info << endl << "\t\t\tphase          (z) : " << pphase2->getInfo();

        // delay phase
        vector<double> delay_phase( 2, 0. );
        PyTools::extractV( "delay_phase", delay_phase, "Laser", ilaser );
        info << endl << "\t\tdelay phase      (y) : " << delay_phase[0];
        info << endl << "\t\tdelay phase      (z) : " << delay_phase[1];

        // Create the LaserProfiles
        profiles.push_back( new LaserProfileSeparable( omega, pchirp, ptime, pspace1, pphase1, delay_phase[0], true ) );
        profiles.push_back( new LaserProfileSeparable( omega, pchirp2, ptime2, pspace2, pphase2, delay_phase[1], false ) );

    }

    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
    }
}


// Cloning constructor
Laser::Laser( Laser *laser, Params &params )
{
    box_side  = laser->box_side;
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
    Profile *spaceProfile, Profile *phaseProfile, double delay_phase, bool primal
):
    primal_( primal ),
    omega_( omega ),
    timeProfile_( timeProfile ),
    chirpProfile_( chirpProfile ),
    spaceProfile_( spaceProfile ),
    phaseProfile_( phaseProfile ),
    delay_phase_( delay_phase )
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
    delay_phase_( lp->delay_phase_ )
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
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids && (patch->vecSpecies.size() == 0) ) // If not species on the patch, cartesian decomposition
        n_space = params.n_space_region;

    std::vector<unsigned int> oversize(params.oversize);
    if (params.uncoupled_grids && (patch->vecSpecies.size() == 0) ) // If not species on the patch, cartesian decomposition
        oversize = params.region_oversize;

    vector<unsigned int> dim( 2 );
    dim[0] = 1;
    dim[1] = 1;

    if( params.geometry!="1Dcartesian" && params.geometry!="2Dcartesian" && params.geometry!="3Dcartesian" && params.geometry!="AMcylindrical" ) {
        ERROR( "Unknown geometry in laser" );
    }

    // dim[0] for 2D and 3D Cartesian
    if( params.geometry=="2Dcartesian" || params.geometry=="3Dcartesian" ) {
        unsigned int ny_p = n_space[1]+1+2*oversize[1];
        unsigned int ny_d = ny_p+1;
        dim[0] = primal_ ? ny_p : ny_d;
    }

    // dim[0] for LRT
    if( params.geometry=="AMcylindrical" ) {
        unsigned int nr_p = n_space[1]+1+2*oversize[1];
        unsigned int nr_d = nr_p+1;
        dim[0] = nr_p + nr_d;
    }

    // dim[1] for 3D Cartesian
    if( params.geometry=="3Dcartesian" ) {
        unsigned int nz_p = n_space[2]+1+2*oversize[2];
        unsigned int nz_d = nz_p+1;
        dim[1] = primal_ ? nz_d : nz_p;
    }

    //Create laser fields
    space_envelope = new Field2D( dim );
    phase          = new Field2D( dim );
}

void LaserProfileSeparable::initFields( Params &params, Patch *patch )
{
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids && (patch->vecSpecies.size() == 0) ) // If not species on the patch, cartesian decomposition
        n_space = params.n_space_region;

    std::vector<unsigned int> oversize(params.oversize);
    if (params.uncoupled_grids && (patch->vecSpecies.size() == 0) ) // If not species on the patch, cartesian decomposition
        oversize = params.region_oversize;

    if( params.geometry=="1Dcartesian" ) {

        // Assign profile (only one point in 1D)
        vector<double> pos( 1 );
        pos[0] = 0.;
        ( *space_envelope )( 0, 0 ) = spaceProfile_->valueAt( pos );
        ( *phase )( 0, 0 ) = phaseProfile_->valueAt( pos );

    } else if( params.geometry=="2Dcartesian" ) {
        
        unsigned int ny_p = n_space[1]+1+2*oversize[1];
        unsigned int ny_d = ny_p+1;
        double dy = params.cell_length[1];
        vector<unsigned int> dim( 1 );
        dim[0] = primal_ ? ny_p : ny_d;

        // Assign profile
        vector<double> pos( 1 );
        pos[0] = patch->getDomainLocalMin( 1 ) - ( ( primal_?0.:0.5 ) + oversize[1] )*dy;
        for( unsigned int j=0 ; j<dim[0] ; j++ ) {
            ( *space_envelope )( j, 0 ) = spaceProfile_->valueAt( pos );
            ( *phase )( j, 0 ) = phaseProfile_->valueAt( pos );
            pos[0] += dy;
        }

    } else if( params.geometry=="AMcylindrical" ) {
    
        unsigned int nr_p = n_space[1]+1+2*oversize[1];
        unsigned int nr_d = nr_p+1;
        double dr = params.cell_length[1];
        vector<unsigned int> dim( 1 );
        dim[0] = nr_p + nr_d; // Need to account for both primal and dual positions

        // Assign profile
        vector<double> pos( 1 );
        for( unsigned int j=0 ; j<dim[0] ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( j*0.5 - 0.5 - oversize[1] )*dr ; // Increment half cells
            ( *space_envelope )( j, 0 ) = spaceProfile_->valueAt( pos );
            ( *phase )( j, 0 ) = phaseProfile_->valueAt( pos );
        }

    } else if( params.geometry=="3Dcartesian" ) {
        unsigned int ny_p = n_space[1]+1+2*oversize[1];
        unsigned int ny_d = ny_p+1;
        unsigned int nz_p = n_space[2]+1+2*oversize[2];
        unsigned int nz_d = nz_p+1;
        double dy = params.cell_length[1];
        double dz = params.cell_length[2];
        vector<unsigned int> dim( 2 );
        dim[0] = primal_ ? ny_p : ny_d;
        dim[1] = primal_ ? nz_d : nz_p;

        // Assign profile
        vector<double> pos( 2 );
        pos[0] = patch->getDomainLocalMin( 1 ) - ( ( primal_?0.:0.5 ) + oversize[1] )*dy;
        for( unsigned int j=0 ; j<dim[0] ; j++ ) {
            pos[1] = patch->getDomainLocalMin( 2 ) - ( ( primal_?0.5:0. ) + oversize[2] )*dz;
            for( unsigned int k=0 ; k<dim[1] ; k++ ) {
                ( *space_envelope )( j, k ) = spaceProfile_->valueAt( pos );
                ( *phase )( j, k ) = phaseProfile_->valueAt( pos );
                pos[1] += dz;
            }
            pos[0] += dy;
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
        ERROR( "Unknown geometry in LaserOffset (cartesian 2D or 3D only)" );
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
    
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids && (patch->vecSpecies.size() == 0) ) // If not species on the patch, cartesian decomposition
        n_space = params.n_space_region;

    std::vector<unsigned int> oversize(params.oversize);
    if (params.uncoupled_grids && (patch->vecSpecies.size() == 0) ) // If not species on the patch, cartesian decomposition
        oversize = params.region_oversize;

    unsigned int ndim = 2;
    if( params.geometry=="3Dcartesian" ) {
        ndim = 3;
    }

    // Define the part of the array to obtain
    vector<hsize_t> dim( 3 ), offset( 3 ), fulldim( 3 );
    hsize_t ny_tot = params.n_space_global[1]+2+2*params.oversize[1];
    hsize_t ny_p = n_space[1]+1+2*oversize[1];
    hsize_t ny_d = ny_p+1;
    dim[0] = primal_ ? ny_p : ny_d;
    dim[1] = 1;
    //offset[0] = patch->getCellStartingGlobalIndex( 1 ) + params.oversize[1];
    offset[1] = 0;
    offset[2] = 0;

    if( ndim == 3 ) {
        hsize_t nz_p = n_space[2]+1+2*oversize[2];
        hsize_t nz_d = nz_p+1;
        dim[1] = primal_ ? nz_d : nz_p;
        offset[1] = patch->getCellStartingGlobalIndex( 2 ) + oversize[2];
    }

    // Open file
    ifstream f( file.c_str() );
    if( !f.good() ) {
        ERROR( "File " << file << " not found" );
    }
    if( H5Fis_hdf5( file.c_str() ) <= 0 ) {
        ERROR( "File " << file << " is not HDF5" );
    }
    hid_t fid = H5Fopen( file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

    // Obtain the omega dataset containing the different values of omega
    hid_t pid = H5Pcreate( H5P_DATASET_ACCESS );
    hssize_t npoints;
    if( H5Lexists( fid, "omega", H5P_DEFAULT ) >0 ) {
        hid_t did = H5Dopen( fid, "omega", pid );
        hid_t filespace = H5Dget_space( did );
        npoints = H5Sget_simple_extent_npoints( filespace );
        omega.resize( npoints );
        H5Dread( did, H5T_NATIVE_DOUBLE, filespace, filespace, H5P_DEFAULT, &omega[0] );
        H5Sclose( filespace );
        H5Dclose( did );
    } else {
        ERROR( "File " << file << " does not contain the `omega` dataset" );
    }

    // Allocate the fields
    magnitude->allocateDims( dim[0], dim[1], npoints );
    phase    ->allocateDims( dim[0], dim[1], npoints );

    // Obtain the datasets for the magnitude and phase of the field
    dim[ndim-1] = npoints;
    hid_t memspace = H5Screate_simple( ndim, &dim[0], NULL );
    string magnitude_name = primal_?"magnitude1":"magnitude2";
    if( H5Lexists( fid, magnitude_name.c_str(), H5P_DEFAULT ) >0 ) {
        hid_t did = H5Dopen( fid, magnitude_name.c_str(), pid );
        hid_t filespace = H5Dget_space( did );
        H5Sget_simple_extent_dims( filespace, &fulldim[0], NULL );
        offset[0] = fulldim[0] - ny_tot + patch->getCellStartingGlobalIndex( 1 ) + params.oversize[1];
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &offset[0], NULL, &dim[0], NULL );
        H5Dread( did, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &magnitude->data_[0] );
        H5Sclose( filespace );
        H5Dclose( did );
    } else {
        magnitude->put_to( 0. );
    }
    string phase_name = primal_?"phase1":"phase2";
    if( H5Lexists( fid, phase_name.c_str(), H5P_DEFAULT ) >0 ) {
        hid_t did = H5Dopen( fid, phase_name.c_str(), pid );
        hid_t filespace = H5Dget_space( did );
        H5Sget_simple_extent_dims( filespace, &fulldim[0], NULL );
        offset[0] = fulldim[0] - ny_tot + patch->getCellStartingGlobalIndex( 1 ) + params.oversize[1];
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &offset[0], NULL, &dim[0], NULL );
        H5Dread( did, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &phase->data_[0] );
        H5Sclose( filespace );
        H5Dclose( did );
    } else {
        phase->put_to( 0. );
    }

    H5Sclose( memspace );
    H5Pclose( pid );
    H5Fclose( fid );
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
