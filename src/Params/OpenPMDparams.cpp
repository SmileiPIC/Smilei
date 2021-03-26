#include "OpenPMDparams.h"

using namespace std;

OpenPMDparams::OpenPMDparams( Params &p ):
    axisLabels( 1 ),
    fieldBoundary( 10 ),
    fieldBoundaryParameters( 13 ),
    particleBoundary( 0 ),
    particleBoundaryParameters( 0 ),
    params( &p )
{

    version = "1.0.0";
    
    extension = 0;
    //extension = 1; // ED-PIC extension. Not supported yet
    
    // Grid parameters
    string xyz = "xyz";
    gridGlobalOffset.resize( params->nDim_field );
    gridOffset      .resize( params->nDim_field );
    position        .resize( params->nDim_field );
    gridSpacing     .resize( params->nDim_field );
    for( unsigned int idim=0; idim<params->nDim_field; idim++ ) {
        axisLabels.addString( xyz.substr( idim, 1 ) );
        gridGlobalOffset[idim] = 0.;
        gridOffset      [idim] = 0.;
        position        [idim] = 0.;
        gridSpacing     [idim] = params->cell_length[idim];
    }
    fieldSolverParameters = "";
    if( params->maxwell_sol == "Yee" ) {
        fieldSolver = "Yee";
    } else if( params->maxwell_sol == "Lehe" ) {
        fieldSolver = "Lehe";
    } else {
        fieldSolver = "other";
        fieldSolverParameters = params->maxwell_sol;
    }
    patchSize = params->n_space;
    
    // Units
    unitDimension.resize( SMILEI_NUNITS );
    unitSI.resize( SMILEI_NUNITS );
    double Wr = params->reference_angular_frequency_SI;
    if( Wr==0. ) {
        Wr=1.;
    }
    for( unsigned int unit_type=0; unit_type<SMILEI_NUNITS; unit_type++ ) {
        unitDimension[unit_type].resize( 7, 0. );
        if( unit_type == SMILEI_UNIT_NONE ) {         // dimensionless
            unitSI[unit_type] = 1.;
        } else if( unit_type == SMILEI_UNIT_EFIELD ) {
            unitDimension[unit_type][0] = 1.;
            unitDimension[unit_type][1] = 1.;
            unitDimension[unit_type][2] = -3.;
            unitDimension[unit_type][3] = -1.;
            unitSI[unit_type] = 1.704508807123e-3 * Wr; // me * c * Wr / e
        } else if( unit_type == SMILEI_UNIT_BFIELD ) {
            unitDimension[unit_type][1] = 1.;
            unitDimension[unit_type][2] = -2.;
            unitDimension[unit_type][3] = -1.;
            unitSI[unit_type] = 5.685629380e-12 * Wr; // me * Wr / e
        } else if( unit_type == SMILEI_UNIT_CURRENT ) {
            unitDimension[unit_type][0] = -2.;
            unitDimension[unit_type][3] = 1.;
            unitSI[unit_type] = 1.5092041114e-14 * Wr*Wr; // e0 * me * c * Wr^2 / e
        } else if( unit_type == SMILEI_UNIT_DENSITY ) {
            unitDimension[unit_type][0] = -3.;
            unitSI[unit_type] = 3.14207756427e-4 * Wr*Wr; // e0 * me * Wr^2 / e^2
        } else if( unit_type == SMILEI_UNIT_POSITION ) {
            unitDimension[unit_type][0] = -3.;
            unitSI[unit_type] = 299792458. / Wr; // c / Wr
        } else if( unit_type == SMILEI_UNIT_MOMENTUM ) {
            unitDimension[unit_type][0] = -3.;
            unitSI[unit_type] = 2.7309240656e-22; // me * c
        } else if( unit_type == SMILEI_UNIT_CHARGE ) {
            unitDimension[unit_type][0] = -3.;
            unitSI[unit_type] = 1.602176565e-19; // e
        } else if( unit_type == SMILEI_UNIT_TIME ) {
            unitDimension[unit_type][2] = 1.;
            unitSI[unit_type] = 1. / Wr; // 1 / Wr
        }
    }
    
    // Boundary conditions
    for( unsigned int i=0; i<params->EM_BCs.size(); i++ ) {
        for( unsigned int j=0; j<2; j++ ) {
            if( params->EM_BCs[i][j] == "periodic" ) {
                fieldBoundary          .addString( "periodic" );
                fieldBoundaryParameters.addString( "periodic" );
            } else if( params->EM_BCs[i][j] == "reflective" ) {
                fieldBoundary          .addString( "reflecting" );
                fieldBoundaryParameters.addString( "reflecting" );
            } else if( params->EM_BCs[i][j] == "silver-muller" ) {
                fieldBoundary          .addString( "open" );
                fieldBoundaryParameters.addString( "silver-muller" );
            } else if( params->EM_BCs[i][j] == "buneman" ) {
                fieldBoundary          .addString( "open" );
                fieldBoundaryParameters.addString( "buneman" );
            } else if( params->EM_BCs[i][j].substr(0,4) == "ramp" ) {
                fieldBoundary          .addString( "open" );
                fieldBoundaryParameters.addString( params->EM_BCs[i][j] );
            } else {
                ERROR( " impossible boundary condition " );
            }
            particleBoundary          .addString( "" );
            particleBoundaryParameters.addString( "" );
        }
    }
    
    // Other parameters
    currentSmoothing = "none";
    currentSmoothingParameters = "";
    if (params->currentFilter_passes.size() > 0){
        if( *std::max_element(std::begin(params->currentFilter_passes), std::end(params->currentFilter_passes)) > 0 ) {
            currentSmoothing = "Binomial";
            ostringstream t( "" );
            t << "numPasses="<<*std::max_element(std::begin(params->currentFilter_passes), std::end(params->currentFilter_passes));
            currentSmoothingParameters = t.str();
        }
    }
}

void OpenPMDparams::writeRootAttributes( H5Write &location, string meshesPath, string particlesPath )
{
    location.attr( "openPMDextension", extension, H5T_NATIVE_UINT32 );
    location.attr( "openPMD", version );
    location.attr( "basePath", "/data/%T/" );
    location.attr( "software", "Smilei" );
    location.attr( "softwareVersion", __VERSION );
    location.attr( "date", getLocalTime() );
    location.attr( "iterationEncoding", "groupBased" );
    location.attr( "iterationFormat", "/data/%T/" );
    location.attr( "meshesPath", meshesPath );
    location.attr( "particlesPath", particlesPath );
}

void OpenPMDparams::writeBasePathAttributes( H5Write &location, unsigned int itime )
{
    location.attr( "time", ( double )( itime * params->timestep ) );
    location.attr( "dt", ( double )params->timestep );
    location.attr( "timeUnitSI", unitSI[SMILEI_UNIT_TIME] );
}

void OpenPMDparams::writeParticlesAttributes( H5Write &location )
{
}

void OpenPMDparams::writeMeshesAttributes( H5Write &location )
{
    location.attr( "patchSize", patchSize ); // this one is not openPMD
    location.attr( "fieldSolver", fieldSolver );
    location.attr( "fieldSolverParameters", fieldSolverParameters );
    location.attr( "fieldBoundary", fieldBoundary );
    location.attr( "fieldBoundaryParameters", fieldBoundaryParameters );
    location.attr( "particleBoundary", particleBoundary );
    location.attr( "particleBoundaryParameters", particleBoundaryParameters );
    location.attr( "currentSmoothing", currentSmoothing );
    location.attr( "currentSmoothingParameters", currentSmoothingParameters );
    location.attr( "chargeCorrection", "none" );
    location.attr( "chargeCorrectionParameters", "" );
    location.attr( "fieldSmoothing", "none" );
    location.attr( "fieldSmoothingParameters", "" );
}

void OpenPMDparams::writeFieldAttributes( H5Write &location, vector<unsigned int> subgrid_start, vector<unsigned int> subgrid_step )
{
    location.attr( "geometry", "cartesian" );
    location.attr( "dataOrder", "C" );
    location.attr( "axisLabels", axisLabels );
    if( subgrid_start.size() == 0 ) {
        location.attr( "gridSpacing", gridSpacing );
        location.attr( "gridGlobalOffset", gridGlobalOffset );
    } else {
        unsigned int ndim = subgrid_start.size();
        vector<double> subgridSpacing( ndim );
        vector<double> subgridOffset( ndim );
        for( unsigned int i=0; i<ndim; i++ ) {
            subgridSpacing[i] = gridSpacing [i] * subgrid_step [i];
            subgridOffset [i] = gridSpacing [i] * subgrid_start[i];
        }
        location.attr( "gridSpacing", subgridSpacing );
        location.attr( "gridGlobalOffset", subgridOffset );
    }
    //location.attr( "gridOffset", gridOffset);
    location.attr( "gridUnitSI", unitSI[SMILEI_UNIT_POSITION] );
}

void OpenPMDparams::writeSpeciesAttributes( H5Write &location )
{
}

void OpenPMDparams::writeRecordAttributes( H5Write &location, unsigned int unit_type )
{
    location.attr( "unitDimension", unitDimension[unit_type] );
    location.attr( "timeOffset", 0. );
}

void OpenPMDparams::writeFieldRecordAttributes( H5Write &location )
{
    location.attr( "position", position );
}

void OpenPMDparams::writeComponentAttributes( H5Write &location, unsigned int unit_type )
{
    location.attr( "unitSI", unitSI[unit_type] );
}



// WARNING: do not change the format. It is required for OpenPMD compatibility.
string OpenPMDparams::getLocalTime()
{
    time_t t = time( 0 );
    struct tm *now = localtime( & t );
    char buffer[25];
    
    strftime( buffer, 25, "%Y-%m-%d %H:%M:%S %z", now );
    
    string s( buffer, 25 );
    return s;
}

