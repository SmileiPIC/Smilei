#include "OpenPMDparams.h"

using namespace std;

OpenPMDparams::OpenPMDparams(Params& p):
    axisLabels(1),
    fieldBoundary(10),
    fieldBoundaryParameters(13),
    particleBoundary(0), 
    particleBoundaryParameters(0),
    params(&p)
{
    
    version = "1.0.0";
    
    extension = 0;
    //extension = 1; // ED-PIC extension. Not supported yet
    
    string xyz = "xyz";
    gridGlobalOffset.resize( params->nDim_field );
    gridOffset      .resize( params->nDim_field );
    position        .resize( params->nDim_field );
    gridSpacing = params->cell_length;
    for( unsigned int idim=0; idim<params->nDim_field; idim++ ) {
        axisLabels.addString( xyz.substr(idim, 1) );
        gridGlobalOffset[idim] = 0.;
        gridOffset      [idim] = 0.;
        position        [idim] = 0.;
    }
    unitDimension.resize( 4 );
    for( unsigned int field_type=0; field_type<4; field_type++ ) {
        unitDimension[field_type].resize(7, 0.);
        if        ( field_type == 0 ) { // Electric field
            unitDimension[field_type][0] = 1.;
            unitDimension[field_type][1] = 1.;
            unitDimension[field_type][2] = -3.;
            unitDimension[field_type][3] = -1.;
        } else if ( field_type == 1 ) { // Magnetic field
            unitDimension[field_type][1] = 1.;
            unitDimension[field_type][2] = -2.;
            unitDimension[field_type][3] = -1.;
        } else if ( field_type == 2 ) { // Current
            unitDimension[field_type][0] = -2.;
            unitDimension[field_type][3] = 1.;
        } else if ( field_type == 3 ) { // Density
            unitDimension[field_type][0] = -3.;
        }
    }
    fieldSolverParameters = "";
    if       ( params->maxwell_sol == "Yee" ) {
        fieldSolver = "Yee";
    } else if( params->maxwell_sol == "Lehe" ) {
        fieldSolver = "Lehe";
    } else {
        fieldSolver = "other";
        fieldSolverParameters = params->maxwell_sol;
    }
    vector<string> bc_em_type;
    bc_em_type.push_back(params->bc_em_type_x[0]);
    bc_em_type.push_back(params->bc_em_type_x[1]);
    if( params->nDim_field > 1 ) {
        bc_em_type.push_back(params->bc_em_type_y[0]);
        bc_em_type.push_back(params->bc_em_type_y[1]);
        if( params->nDim_field > 2 ) {
            bc_em_type.push_back(params->bc_em_type_z[0]);
            bc_em_type.push_back(params->bc_em_type_z[1]);
        }
    }
    for( unsigned int i=0; i<bc_em_type.size(); i++ ) {
        if( bc_em_type[i] == "periodic" ) {
            fieldBoundary          .addString( "periodic" );
            fieldBoundaryParameters.addString( "periodic" );
        } else if( bc_em_type[i] == "reflective" ) {
            fieldBoundary          .addString( "reflecting" );
            fieldBoundaryParameters.addString( "reflecting" );
        } else if( bc_em_type[i] == "silver-muller" ) {
            fieldBoundary          .addString( "open" );
            fieldBoundaryParameters.addString( "silver-muller");
        } else {
            ERROR(" impossible boundary condition ");
        }
        particleBoundary          .addString( "" );
        particleBoundaryParameters.addString( "" );
    }
    currentSmoothing = "none";
    currentSmoothingParameters = "";
    if( params->currentFilter_int > 0 ) {
        currentSmoothing = "Binomial";
        ostringstream t("");
        t << "numPasses="<<params->currentFilter_int;
        currentSmoothingParameters = t.str();
    }
}

void OpenPMDparams::writeRootAttributes( hid_t location, string meshesPath, string particlesPath )
{
    H5::attr( location, "openPMDextension", extension, H5T_NATIVE_UINT32);
    H5::attr( location, "openPMD", version);
    H5::attr( location, "basePath", "/data/%T/");
    H5::attr( location, "software", "Smilei");
    H5::attr( location, "softwareVersion", __VERSION);
    H5::attr( location, "date", getLocalTime());
    H5::attr( location, "iterationEncoding", "groupBased");
    H5::attr( location, "iterationFormat", "/data/%T/");
    H5::attr( location, "meshesPath", meshesPath);
    H5::attr( location, "particlesPath", particlesPath);
}

void OpenPMDparams::writeBasePathAttributes( hid_t location, unsigned int itime )
{
    H5::attr( location, "time", (double)(itime * params->timestep));
    H5::attr( location, "dt", (double)params->timestep);
    H5::attr( location, "timeUnitSI", 0.); // not relevant
}

void OpenPMDparams::writeParticlesPathAttributes( hid_t location )
{
}

void OpenPMDparams::writeMeshesPathAttributes( hid_t location )
{
    H5::attr( location, "fieldSolver", fieldSolver);
    H5::attr( location, "fieldSolverParameters", fieldSolverParameters);
    H5::attr( location, "fieldBoundary", fieldBoundary);
    H5::attr( location, "fieldBoundaryParameters", fieldBoundaryParameters);
    H5::attr( location, "particleBoundary", particleBoundary);
    H5::attr( location, "particleBoundaryParameters", particleBoundaryParameters);
    H5::attr( location, "currentSmoothing", currentSmoothing);
    H5::attr( location, "currentSmoothingParameters", currentSmoothingParameters);
    H5::attr( location, "chargeCorrection", "none");
    H5::attr( location, "chargeCorrectionParameters", "");
    H5::attr( location, "fieldSmoothing", "none");
    H5::attr( location, "fieldSmoothingParameters", "");
}

void OpenPMDparams::writeFieldAttributes( hid_t location, unsigned int field_type )
{
    H5::attr( location, "geometry", "cartesian");
    H5::attr( location, "dataOrder", "C");
    H5::attr( location, "axisLabels", axisLabels);
    H5::attr( location, "gridSpacing", gridSpacing);
    H5::attr( location, "gridGlobalOffset", gridGlobalOffset);
    H5::attr( location, "gridOffset", gridOffset);
    H5::attr( location, "gridUnitSI", 0.);      
    H5::attr( location, "unitSI", 0.);      
    H5::attr( location, "unitDimension", unitDimension[field_type]);
    H5::attr( location, "timeOffset", 0.);
}

void OpenPMDparams::writeFieldRecordAttributes( hid_t location )
{
    H5::attr( location, "position", position);
}

// WARNING: do not change the format. It is required for OpenPMD compatibility.
string OpenPMDparams::getLocalTime() {
    time_t t = time(0);
    struct tm * now = localtime( & t );
    char buffer[25];
    
    strftime(buffer, 25, "%Y-%m-%d %H:%M:%S %z", now);

    string s(buffer, 25);
    return s;
}

