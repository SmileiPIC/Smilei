#include "Tools.h"
#include "Params.h"
#include "Laser.h"
#include "Patch.h"
#include "H5.h"

#include <cmath>
#include <string>

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
    bool has_time, has_space, has_omega, has_chirp, has_phase, has_space_time, has_file;
    double omega(0);
    Profile *p, *pchirp, *pchirp2, *ptime, *ptime2, *pspace1, *pspace2, *pphase1, *pphase2;
    pchirp2 = NULL;
    ptime2  = NULL;
    file = "";
    has_omega      = PyTools::extract("omega",omega,"Laser",ilaser);
    has_chirp      = PyTools::extract_pyProfile("chirp_profile"     , chirp_profile, "Laser", ilaser);
    has_time       = PyTools::extract_pyProfile("time_envelope"     , time_profile , "Laser", ilaser);
    has_space      = PyTools::extract2Profiles ("space_envelope"    , ilaser, space_profile     );
    has_phase      = PyTools::extract2Profiles ("phase"             , ilaser, phase_profile     );
    has_space_time = PyTools::extract2Profiles ("space_time_profile", ilaser, space_time_profile);
    has_file       = PyTools::extract("file",file,"Laser",ilaser);
    
    if( has_space_time && has_file )
        ERROR(errorPrefix << ": `space_time_profile` and `file` cannot both be set");
    
    
    spacetime.resize(2, false);
    if( has_space_time ) {
        
        spacetime[0] = (bool)(space_time_profile[0]);
        spacetime[1] = (bool)(space_time_profile[1]);
        
        if( has_time || has_space || has_omega || has_chirp || has_phase ) {
            name.str("");
            name << (has_time ?"time_envelope ":"")
                 << (has_space?"space_envelope ":"")
                 << (has_omega?"omega ":"")
                 << (has_chirp?"chirp_profile ":"")
                 << (has_phase?"phase ":"");
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
        
    } else if( has_file ) {
        
        profiles.push_back( new LaserProfileFile(file, true ) );
        profiles.push_back( new LaserProfileFile(file, false) );
        
    } else {
        
        if( ! has_time )
            ERROR(errorPrefix << ": missing `time_envelope`");
        if( ! has_space )
            ERROR(errorPrefix << ": missing `space_envelope`");
        if( ! has_omega )
            ERROR(errorPrefix << ": missing `omega`");
        if( ! has_chirp )
            ERROR(errorPrefix << ": missing `chirp_profile`");
        if( ! has_phase )
            ERROR(errorPrefix << ": missing `phase`");
        
        info << "\t\t" << errorPrefix << ": custom profile" << endl;
        
        unsigned int space_dims = (params.geometry=="3Dcartesian" ? 2 : 1);
 
        // omega
        info << "\t\t\tomega              : " << omega << endl;
        
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
        
        // Create the LaserProfiles
        profiles.push_back( new LaserProfileSeparable(omega, pchirp, ptime, pspace1, pphase1, true ) );
        profiles.push_back( new LaserProfileSeparable(omega, pchirp2, ptime2, pspace2, pphase2, false) );
    
    }
    
    // Display info
    if( patch->isMaster() ) {
        MESSAGE( info.str() );
    }
}


// Cloning constructor
Laser::Laser(Laser* laser, Params& params)
{
    box_side  = laser->box_side;
    spacetime = laser->spacetime;
    file      = laser->file;
    
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
    } else if( file != "" ) {
        profiles.push_back( new LaserProfileFile(file, true ) );
        profiles.push_back( new LaserProfileFile(file, false) );
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
    Profile* spaceProfile, Profile* phaseProfile, bool primal
):
    primal       ( primal       ),
    omega        ( omega        ),
    timeProfile  ( timeProfile  ),
    chirpProfile ( chirpProfile ),
    spaceProfile ( spaceProfile ),
    phaseProfile ( phaseProfile )
{
    space_envelope = NULL;
    phase = NULL;
}
// Separable laser profile cloning constructor
LaserProfileSeparable::LaserProfileSeparable(LaserProfileSeparable * lp) :
    primal       ( lp->primal ),
    omega        ( lp->omega  ),
    timeProfile  ( new Profile(lp->timeProfile ) ),
    chirpProfile ( new Profile(lp->chirpProfile) ),
    spaceProfile ( new Profile(lp->spaceProfile) ),
    phaseProfile ( new Profile(lp->phaseProfile) )
{
    space_envelope = NULL;
    phase = NULL;
}
// Separable laser profile destructor
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
        double t0 = (*phase)(j, k) / omega_;
        amp = timeProfile->valueAt(t-t0) * (*space_envelope)(j, k) * sin( omega_*(t - t0) );
    }
    return amp;
}

//Destructor
LaserProfileNonSeparable::~LaserProfileNonSeparable()
{
    if(spaceAndTimeProfile) delete spaceAndTimeProfile;
}


void LaserProfileFile::createFields(Params& params, Patch* patch)
{
    if( params.geometry!="2Dcartesian" && params.geometry!="3Dcartesian" )
        ERROR("Unknown geometry in LaserOffset (cartesian 2D or 3D only)");
    
    vector<unsigned int> dim(3);
    
    unsigned int ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
    unsigned int ny_d = ny_p+1;
    dim[0] = primal ? ny_p : ny_d;
    dim[1] = 1;
    
    if( params.geometry=="3Dcartesian" ) {
        unsigned int nz_p = params.n_space[2]*params.global_factor[2]+1+2*params.oversize[2];
        unsigned int nz_d = nz_p+1;
        dim[1] = primal ? nz_d : nz_p;
    }
    
    // Obtain the size of the omega dataset in the file
    if( H5Fis_hdf5( file.c_str() ) <= 0 )
        ERROR("File " << file << " not found");
    hid_t fid = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t pid = H5Pcreate( H5P_DATASET_ACCESS );
    hid_t did = H5Dopen( fid, "omega", pid);
    if( did<0 ) {
        ERROR("File " << file << " does not contain the `omega` dataset");
    } else {
        hid_t filespace = H5Dget_space( did );
        hssize_t npoints = H5Sget_simple_extent_npoints( filespace );
        dim[2] = npoints;
        H5Sclose(filespace);
        H5Dclose(did);
    }
    H5Pclose(pid);
    H5Fclose(fid);
    
    magnitude = new Field3D(dim);
    phase     = new Field3D(dim);
}

void LaserProfileFile::initFields(Params& params, Patch* patch)
{
    unsigned int ndim = 2;
    if( params.geometry=="3Dcartesian" ) ndim = 3;
    
    // Define the part of the array to obtain
    vector<hsize_t> dim(3), offset(3);
    hsize_t ny_p = params.n_space[1]*params.global_factor[1]+1+2*params.oversize[1];
    hsize_t ny_d = ny_p+1;
    dim[0] = primal ? ny_p : ny_d;
    offset[0] = patch->getCellStartingGlobalIndex(1) + 4;
    offset[1] = 0;
    offset[2] = 0;
    
    if( ndim == 3 ) {
        hsize_t nz_p = params.n_space[2]*params.global_factor[2]+1+2*params.oversize[2];
        hsize_t nz_d = nz_p+1;
        dim[1] = primal ? nz_d : nz_p;
        offset[1] = patch->getCellStartingGlobalIndex(2) + 4;
    }
    
    // Open file
    if( H5Fis_hdf5( file.c_str() ) <= 0 )
        ERROR("File " << file << " not found");
    hid_t fid = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    // Obtain the omega dataset containing the different values of omega
    hid_t pid = H5Pcreate( H5P_DATASET_ACCESS );
    hid_t did = H5Dopen( fid, "omega", pid );
    hssize_t npoints;
    if( did<0 ) {
        ERROR("File " << file << " does not contain the `omega` dataset");
    } else {
        hid_t filespace = H5Dget_space( did );
        npoints = H5Sget_simple_extent_npoints( filespace );
        omega.resize( npoints );
        H5Dread( did, H5T_NATIVE_DOUBLE, filespace, filespace, H5P_DEFAULT, &omega[0] );
        H5Sclose(filespace);
        H5Dclose(did);
    }
    
    dim[ndim-1] = npoints;
    hid_t memspace = H5Screate_simple( ndim, &dim[0], NULL );
    
    // Obtain the datasets for the magnitude and phase of the field
    did = H5Dopen( fid, (primal?"magnitude1":"magnitude2"), pid );
    if( did >=0 ) {
        hid_t filespace = H5Dget_space( did );
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dim[0], NULL );
        H5Dread( did, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &magnitude->data_[0] );
        H5Sclose(filespace);
        H5Dclose(did);
    } else {
        magnitude->put_to(0.);
    }
    did = H5Dopen( fid, (primal?"phase1":"phase2"), pid );
    if( did >=0 ) {
        hid_t filespace = H5Dget_space( did );
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dim[0], NULL );
        H5Dread( did, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &phase->data_[0] );
        H5Sclose(filespace);
        H5Dclose(did);
    } else {
        phase->put_to(0.);
    }
    
    H5Sclose(memspace);
    H5Pclose(pid);
    H5Fclose(fid);
}

// Amplitude of a laser profile from a file (see LaserOffset)
double LaserProfileFile::getAmplitude(std::vector<double> pos, double t, int j, int k)
{
    double amp = 0;
    //#pragma omp critical
    //{
    //    timeProfile->valueAt(t)
    //}
    unsigned int n = omega.size();
    for( unsigned int i=0; i<n; i++ ) {
        amp += (*magnitude)(j,k,i) * sin( omega[i] * t + (*phase)(j,k,i) );
    }
    return amp;
}

//Destructor
LaserProfileFile::~LaserProfileFile()
{
    if(magnitude) delete magnitude;
    if(phase    ) delete phase    ;
}
