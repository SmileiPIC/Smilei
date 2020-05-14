#include "PyTools.h"
#include <iomanip>

#include "DiagnosticParticleBinningBase.h"
#include "HistogramFactory.h"


using namespace std;


// Constructor
DiagnosticParticleBinningBase::DiagnosticParticleBinningBase(
    Params &params,
    SmileiMPI *smpi,
    Patch *patch,
    int diagId,
    string diagName,
    bool time_accumulate_,
    PyObject *deposited_quantity,
    vector<string> excluded_axes
) : Diagnostic( nullptr, Tools::merge( "Diag", diagName ), diagId )
{
    fileId_ = 0;
    int idiag = diagId;
    time_accumulate = time_accumulate_;
    
    string pyDiag = Tools::merge( "Diag", diagName );
    string errorPrefix = Tools::merge( pyDiag, " #", to_string( idiag ) );
    
    // get parameter "deposited_quantity" that determines the quantity to sum in the output array
    if( deposited_quantity == nullptr ) {
        deposited_quantity = PyTools::extract_py( "deposited_quantity", pyDiag, idiag );
        PyTools::checkPyError();
    }
    
    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py( "every", pyDiag, idiag ),
        errorPrefix
    );
    
    // get parameter "flush_every" which describes a timestep selection for flushing the file
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py( "flush_every", pyDiag, idiag ),
        errorPrefix
    );
    
    // get parameter "time_average" that determines the number of timestep to average the outputs
    if( ! time_accumulate ) {
        time_average = 1;
        PyTools::extract( "time_average", time_average, pyDiag, idiag );
        if( time_average < 1 ) {
            time_average=1;
        }
        if( time_average > timeSelection->smallestInterval() ) {
            ERROR( errorPrefix << ": `time_average` is incompatible with `every`" );
        }
    }
    
    // get parameter "species" that determines the species to use (can be a list of species)
    vector<string> species_names;
    if( ! PyTools::extractV( "species", species_names, pyDiag, idiag ) ) {
        ERROR( errorPrefix << ": parameter `species` required" );
    }
    // verify that the species exist, remove duplicates and sort by number
    species = params.FindSpecies( patch->vecSpecies, species_names );
    
//    // Temporarily set the spatial min and max to the simulation box size
//    spatial_min.resize( params.nDim_particle, 0. );
//    spatial_max = params.grid_length;
    
    // get parameter "axes" that adds axes to the diagnostic
    // Each axis should contain several items:
    //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
    vector<PyObject *> pyAxes=PyTools::extract_pyVec( "axes", pyDiag, idiag );
    
    // Create the Histogram object based on the extracted parameters above
    histogram = HistogramFactory::create( params, deposited_quantity, pyAxes, species, patch, excluded_axes, errorPrefix );
    total_axes = histogram->axes.size();
    dims.resize( total_axes );
    for( int iaxis=0; iaxis<total_axes; iaxis++ ) {
        dims[iaxis] = histogram->axes[iaxis]->nbins;
    }
    
//    // Get info on the spatial extent
//    for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
//        if( histogram->axes[i]->type == "x" ) {
//            spatial_min[0] = histogram->axes[i]->min;
//            spatial_max[0] = histogram->axes[i]->max;
//        } else if( histogram->axes[i]->type == "y" ) {
//            spatial_min[1] = histogram->axes[i]->min;
//            spatial_max[1] = histogram->axes[i]->max;
//        } else if( histogram->axes[i]->type == "z" ) {
//            spatial_min[2] = histogram->axes[i]->min;
//            spatial_max[2] = histogram->axes[i]->max;
//        }
//    }
    
    // Calculate the size of the output array
    uint64_t total_size = 1;
    for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
        total_size *= histogram->axes[i]->nbins;
    }
    if( total_size > 2147483648 ) { // 2^31
        ERROR( errorPrefix << ": too many points (" << total_size << " > 2^31)" );
    }
    output_size = ( unsigned int ) total_size;
    
    // Output info on diagnostics
    if( smpi->isMaster() ) {
        ostringstream mystream( "" );
        mystream.str( "" );
        mystream << species_names[0];
        for( unsigned int i=1; i<species_names.size(); i++ ) {
            mystream << "," << species_names[i];
        }
        MESSAGE( 1, "Created "<<diagName<<" #" << idiag << ": species " << mystream.str() );
        for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
            MESSAGE( 2, histogram->axes[i]->info() );
        }
        
        // init HDF files (by master, only if it doesn't yet exist)
        mystream.str( "" ); // clear
        mystream << diagName << idiag << ".h5";
        filename = mystream.str();
    }
    
} // END DiagnosticParticleBinning::DiagnosticParticleBinning


DiagnosticParticleBinningBase::~DiagnosticParticleBinningBase()
{
    delete histogram;

    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticParticleBinning::~DiagnosticParticleBinning


// Called only by patch master of process master
void DiagnosticParticleBinningBase::openFile( Params &params, SmileiMPI *smpi, bool newfile )
{
    if( !smpi->isMaster() || fileId_>0 ) {
        return;
    }
    
    if( newfile ) {
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        // write all parameters as HDF5 attributes
        H5::attr( fileId_, "Version", string( __VERSION ) );
        H5::attr( fileId_, "name", diag_name_ );
        H5::attr( fileId_, "deposited_quantity", histogram->deposited_quantity );
        if( ! time_accumulate ) {
            H5::attr( fileId_, "time_average", time_average );
        }
        // write all species
        ostringstream mystream( "" );
        mystream.str( "" ); // clear
        for( unsigned int i=0 ; i < species.size() ; i++ ) {
            mystream << species[i] << " ";
        }
        H5::attr( fileId_, "species", mystream.str() );
        // write each axis
        for( unsigned int iaxis=0 ; iaxis < histogram->axes.size() ; iaxis++ ) {
            mystream.str( "" ); // clear
            mystream << "axis" << iaxis;
            string str1 = mystream.str();
            mystream.str( "" ); // clear
            mystream << histogram->axes[iaxis]->type << " " << histogram->axes[iaxis]->min << " " << histogram->axes[iaxis]->max << " "
                     << histogram->axes[iaxis]->nbins << " " << histogram->axes[iaxis]->logscale << " " << histogram->axes[iaxis]->edge_inclusive << " [";
            for( unsigned int idim=0; idim<histogram->axes[iaxis]->coefficients.size(); idim++ ) {
                mystream << histogram->axes[iaxis]->coefficients[idim];
                if( idim<histogram->axes[iaxis]->coefficients.size()-1 ) {
                    mystream << ",";
                }
            }
            mystream << "]";
            string str2 = mystream.str();
            H5::attr( fileId_, str1, str2 );
        }
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    } else {
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
    }
}


void DiagnosticParticleBinningBase::closeFile()
{
    if( fileId_!=0 ) {
        H5Fclose( fileId_ );
        fileId_ = 0;
    }
    
} // END closeFile


bool DiagnosticParticleBinningBase::prepare( int timestep )
{
    // Get the previous timestep of the time selection
    int previousTime = timeSelection->previousTime( timestep );
    
    // Leave if the timestep is not the good one
    if( timestep - previousTime >= time_average ) {
        return false;
    }
    
    // Allocate memory for the output array (already done if time-averaging)
    data_sum.resize( output_size );
    
    // if first time, erase output array
    if( timestep == previousTime ) {
        fill( data_sum.begin(), data_sum.end(), 0. );
    }
    
    return true;
    
} // END prepare


// run one particle binning diagnostic
void DiagnosticParticleBinningBase::run( Patch *patch, int timestep, SimWindow *simWindow )
{

    vector<int> int_buffer;
    vector<double> double_buffer;
    unsigned int npart;
    
//    // Update spatial_min and spatial_max if needed
//    if( simWindow ) {
//        bool did_move = false;
//        for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
//            if( histogram->axes[i]->type == "moving_x" ) {
//                spatial_min[0] = histogram->axes[i]->min + simWindow->getXmoved();
//                spatial_max[0] = histogram->axes[i]->max + simWindow->getXmoved();
//                did_move = true;
//            }
//        }
//        if( ! did_move ) {
//            spatial_max[0] += simWindow->getXmoved() - spatial_min[0];
//            spatial_min[0] = simWindow->getXmoved();
//        }
//    }
    
    // loop species
    for( unsigned int ispec=0 ; ispec < species.size() ; ispec++ ) {
    
        Species *s = patch->vecSpecies[species[ispec]];
        npart = s->particles->size();
        int_buffer   .resize( npart );
        double_buffer.resize( npart );
        
        fill( int_buffer.begin(), int_buffer.end(), 0 );
        
        histogram->digitize( s, double_buffer, int_buffer, simWindow );
        histogram->valuate( s, double_buffer, int_buffer );
        histogram->distribute( double_buffer, int_buffer, data_sum );
        
    }
    
} // END run

bool DiagnosticParticleBinningBase::writeNow( int timestep ) {
    return timestep - timeSelection->previousTime() == time_average-1;
}

// Now the data_sum has been filled
// if needed now, store result to hdf file
void DiagnosticParticleBinningBase::write( int timestep, SmileiMPI *smpi )
{
    if( !smpi->isMaster() || !writeNow( timestep ) ) {
        return;
    }
    
    // if time_average, then we need to divide by the number of timesteps
    if( !time_accumulate && time_average > 1 ) {
        double coeff = 1./( ( double )time_average );
        for( unsigned int i=0; i<output_size; i++ ) {
            data_sum[i] *= coeff;
        }
    }
    
    // make name of the array
    ostringstream mystream( "" );
    mystream.str( "" );
    mystream << "timestep" << setw( 8 ) << setfill( '0' ) << timestep;
    
    // write the array if it does not exist already
    if( ! H5Lexists( fileId_, mystream.str().c_str(), H5P_DEFAULT ) ) {
        // Create file space
        hid_t sid = H5Screate_simple( total_axes, &dims[0], NULL );
        hid_t pid = H5Pcreate( H5P_DATASET_CREATE );
        // create dataset
        hid_t did = H5Dcreate( fileId_, mystream.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid, H5P_DEFAULT );
        // write vector in dataset
        H5Dwrite( did, H5T_NATIVE_DOUBLE, sid, sid, H5P_DEFAULT, &data_sum[0] );
        // close all
        H5Dclose( did );
        H5Pclose( pid );
        H5Sclose( sid );
    }
    
    if( flush_timeSelection->theTimeIsNow( timestep ) ) {
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
    
    if( ! time_accumulate ) {
        // Clear the array
        clear();
        data_sum.resize( 0 );
    }
} // END write


//! Clear the array
void DiagnosticParticleBinningBase::clear()
{
    data_sum.resize( 0 );
    vector<double>().swap( data_sum );
}


// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticParticleBinningBase::getDiskFootPrint( int istart, int istop, Patch *patch )
{
    uint64_t footprint = 0;
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );
    
    // Add necessary global headers approximately
    footprint += 1500;
    
    // Add necessary timestep headers approximately
    footprint += ndumps * 640;
    
    // Add size of each dump
    footprint += ndumps * ( uint64_t )( output_size ) * 8;
    
    return footprint;
}

