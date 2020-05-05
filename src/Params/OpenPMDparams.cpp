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
            } else if( params->EM_BCs[i][j] == "zero" ) {
                fieldBoundary          .addString( "open" );
                fieldBoundaryParameters.addString( "zero" );
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

void OpenPMDparams::writeRootAttributes( hid_t location, string meshesPath, string particlesPath )
{
    H5::attr( location, "openPMDextension", extension, H5T_NATIVE_UINT32 );
    H5::attr( location, "openPMD", version );
    H5::attr( location, "basePath", "/data/%T/" );
    H5::attr( location, "software", "Smilei" );
    H5::attr( location, "softwareVersion", __VERSION );
    H5::attr( location, "date", getLocalTime() );
    H5::attr( location, "iterationEncoding", "groupBased" );
    H5::attr( location, "iterationFormat", "/data/%T/" );
    H5::attr( location, "meshesPath", meshesPath );
    H5::attr( location, "particlesPath", particlesPath );
}

void OpenPMDparams::writeBasePathAttributes( hid_t location, unsigned int itime )
{
    H5::attr( location, "time", ( double )( itime * params->timestep ) );
    H5::attr( location, "dt", ( double )params->timestep );
    H5::attr( location, "timeUnitSI", unitSI[SMILEI_UNIT_TIME] );
}

void OpenPMDparams::writeParticlesAttributes( hid_t location )
{
}

void OpenPMDparams::writeMeshesAttributes( hid_t location )
{
    H5::attr( location, "patchSize", patchSize ); // this one is not openPMD
    H5::attr( location, "fieldSolver", fieldSolver );
    H5::attr( location, "fieldSolverParameters", fieldSolverParameters );
    H5::attr( location, "fieldBoundary", fieldBoundary );
    H5::attr( location, "fieldBoundaryParameters", fieldBoundaryParameters );
    H5::attr( location, "particleBoundary", particleBoundary );
    H5::attr( location, "particleBoundaryParameters", particleBoundaryParameters );
    H5::attr( location, "currentSmoothing", currentSmoothing );
    H5::attr( location, "currentSmoothingParameters", currentSmoothingParameters );
    H5::attr( location, "chargeCorrection", "none" );
    H5::attr( location, "chargeCorrectionParameters", "" );
    H5::attr( location, "fieldSmoothing", "none" );
    H5::attr( location, "fieldSmoothingParameters", "" );
}

void OpenPMDparams::writeFieldAttributes( hid_t location, vector<unsigned int> subgrid_start, vector<unsigned int> subgrid_step )
{
    H5::attr( location, "geometry", "cartesian" );
    H5::attr( location, "dataOrder", "C" );
    H5::attr( location, "axisLabels", axisLabels );
    if( subgrid_start.size() == 0 ) {
        H5::attr( location, "gridSpacing", gridSpacing );
        H5::attr( location, "gridGlobalOffset", gridGlobalOffset );
    } else {
        unsigned int ndim = subgrid_start.size();
        vector<double> subgridSpacing( ndim );
        vector<double> subgridOffset( ndim );
        for( unsigned int i=0; i<ndim; i++ ) {
            subgridSpacing[i] = gridSpacing [i] * subgrid_step [i];
            subgridOffset [i] = gridSpacing [i] * subgrid_start[i];
        }
        H5::attr( location, "gridSpacing", subgridSpacing );
        H5::attr( location, "gridGlobalOffset", subgridOffset );
    }
    //H5::attr( location, "gridOffset", gridOffset);
    H5::attr( location, "gridUnitSI", unitSI[SMILEI_UNIT_POSITION] );
}

void OpenPMDparams::writeSpeciesAttributes( hid_t location )
{
}

void OpenPMDparams::writeRecordAttributes( hid_t location, unsigned int unit_type )
{
    H5::attr( location, "unitDimension", unitDimension[unit_type] );
    H5::attr( location, "timeOffset", 0. );
}

void OpenPMDparams::writeFieldRecordAttributes( hid_t location )
{
    H5::attr( location, "position", position );
}

void OpenPMDparams::writeComponentAttributes( hid_t location, unsigned int unit_type )
{
    H5::attr( location, "unitSI", unitSI[unit_type] );
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

