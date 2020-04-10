#include "PyTools.h"
#include <iomanip>

#include "DiagnosticPerformances.h"


using namespace std;

const unsigned int n_quantities_double = 15;
const unsigned int n_quantities_uint   = 4;

// Constructor
DiagnosticPerformances::DiagnosticPerformances( Params &params, SmileiMPI *smpi )
{
    fileId_ = 0;
    timestep = params.timestep;
    cell_load = params.cell_load;
    frozen_particle_load = params.frozen_particle_load;
    tot_number_of_patches = params.tot_number_of_patches;
    
    ostringstream name( "" );
    name << "Diagnostic performances";
    string errorPrefix = name.str();
    
    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py( "every", "DiagPerformances" ),
        name.str()
    );
    
    // get parameter "flush_every" which describes a timestep selection for flushing the file
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py( "flush_every", "DiagPerformances" ),
        name.str()
    );
    
    // Get patch information flag
    PyTools::extract( "patch_information", patch_information, "DiagPerformances"  );
    
    // Output info on diagnostics
    if( smpi->isMaster() ) {
        MESSAGE( 1, "Created performances diagnostic" );
    }
    filename = "Performances.h5";
    
    mpi_size = smpi->getSize();
    mpi_rank_ = smpi->getRank();
    ndim     = params.nDim_field;
    has_adaptive_vectorization = params.has_adaptive_vectorization;
    
    // Define the HDF5 file and memory spaces
    setHDF5spaces( filespace_double, memspace_double, n_quantities_double, mpi_size, mpi_rank_ );
    setHDF5spaces( filespace_uint, memspace_uint, n_quantities_uint, mpi_size, mpi_rank_ );
    
    // Define HDF5 file access
    write_plist = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( write_plist, H5FD_MPIO_COLLECTIVE );
    
    // Calculate the number of cells per patch
    ncells_per_patch = 1;
    for( unsigned int idim = 0; idim < params.nDim_field; idim++ ) {
        ncells_per_patch *= params.n_space[idim]+2*params.oversize[idim];
    }
    
} // END DiagnosticPerformances::DiagnosticPerformances


DiagnosticPerformances::~DiagnosticPerformances()
{
    H5Pclose( write_plist );
    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticPerformances::~DiagnosticPerformances


// Called only by patch master of process master
void DiagnosticPerformances::openFile( Params &params, SmileiMPI *smpi, bool newfile )
{
    if( fileId_>0 ) {
        return;
    }
    
    if( newfile ) {
        hid_t pid = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( pid, MPI_COMM_WORLD, MPI_INFO_NULL );
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid );
        H5Pclose( pid );
        
        // write all parameters as HDF5 attributes
        H5::attr( fileId_, "MPI_SIZE", smpi->getSize() );
        H5::attr( fileId_, "patch_arrangement", params.patch_arrangement );
        
        vector<string> quantities_uint( n_quantities_uint );
        quantities_uint[0] = "hindex"                    ;
        quantities_uint[1] = "number_of_cells"           ;
        quantities_uint[2] = "number_of_particles"       ;
        quantities_uint[3] = "number_of_frozen_particles";
        H5::attr( fileId_, "quantities_uint", quantities_uint );
        
        vector<string> quantities_double( n_quantities_double );
        quantities_double[ 0] = "total_load"      ;
        quantities_double[ 1] = "timer_global"    ;
        quantities_double[ 2] = "timer_particles" ;
        quantities_double[ 3] = "timer_maxwell"   ;
        quantities_double[ 4] = "timer_densities" ;
        quantities_double[ 5] = "timer_collisions";
        quantities_double[ 6] = "timer_movWindow" ;
        quantities_double[ 7] = "timer_loadBal"   ;
        quantities_double[ 8] = "timer_syncPart"  ;
        quantities_double[ 9] = "timer_syncField" ;
        quantities_double[10] = "timer_syncDens"  ;
        quantities_double[11] = "timer_diags"     ;
        quantities_double[12] = "timer_grids"     ;
        quantities_double[13] = "timer_total"     ;
        quantities_double[14] = "memory_total"    ;
        H5::attr( fileId_, "quantities_double", quantities_double );
        
    } else {
        // Open the existing file
        hid_t pid = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( pid, MPI_COMM_WORLD, MPI_INFO_NULL );
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose( pid );
    }
}


void DiagnosticPerformances::closeFile()
{
    if( filespace_uint>0 ) {
        H5Sclose( filespace_uint );
    }
    if( memspace_uint >0 ) {
        H5Sclose( memspace_uint );
    }
    if( filespace_double>0 ) {
        H5Sclose( filespace_double );
    }
    if( memspace_double >0 ) {
        H5Sclose( memspace_double );
    }
    if( fileId_  >0 ) {
        H5Fclose( fileId_ );
    }
    fileId_ = 0;
} // END closeFile



void DiagnosticPerformances::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}


bool DiagnosticPerformances::prepare( int itime )
{
    if( timeSelection->theTimeIsNow( itime ) ) {
        return true;
    } else {
        return false;
    }
} // END prepare


void DiagnosticPerformances::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers )
{

    #pragma omp master
    {
        // Create group for this iteration
        ostringstream name_t;
        name_t.str( "" );
        name_t << setfill( '0' ) << setw( 10 ) << itime;
        status = H5Lexists( fileId_, name_t.str().c_str(), H5P_DEFAULT );
        if( status==0 )
        {
            iteration_group_id = H5::group( fileId_, name_t.str().c_str() );
        }
        // Warning if file unreachable
        if( status < 0 )
        {
            WARNING( "Performances diagnostic could not write" );
        }
    }
    #pragma omp barrier
    
    // Do not output diag if this iteration has already been written or if problem with file
    if( status != 0 ) {
        return;
    }
    
    #pragma omp master
    {
        // Calculate the loads
        unsigned int number_of_patches = vecPatches.size();
        unsigned int number_of_cells = ncells_per_patch * number_of_patches;
        unsigned int number_of_species = vecPatches( 0 )->vecSpecies.size();
        unsigned int number_of_particles=0, number_of_frozen_particles=0;
        double time = itime * timestep;
        for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
            for( unsigned int ispecies = 0; ispecies < number_of_species; ispecies++ ) {
                if( time < vecPatches( ipatch )->vecSpecies[ispecies]->time_frozen_ ) {
                    number_of_frozen_particles += vecPatches( ipatch )->vecSpecies[ispecies]->getNbrOfParticles();
                } else {
                    number_of_particles += vecPatches( ipatch )->vecSpecies[ispecies]->getNbrOfParticles();
                }
            }
        }
        double total_load =
            ( ( double )number_of_particles )
            + ( ( double )number_of_frozen_particles ) * frozen_particle_load
            + ( ( double )number_of_cells ) * cell_load;
            
        hid_t create_plist = H5Pcreate( H5P_DATASET_CREATE );
        
        // Fill the vector for uint quantities
        vector<unsigned int> quantities_uint( n_quantities_uint );
        quantities_uint[0] = vecPatches( 0 )->Hindex()   ;
        quantities_uint[1] = number_of_cells           ;
        quantities_uint[2] = number_of_particles       ;
        quantities_uint[3] = number_of_frozen_particles;
        
        // Write uints to file
        hid_t dset_uint  = H5Dcreate( iteration_group_id, "quantities_uint", H5T_NATIVE_UINT, filespace_uint, H5P_DEFAULT, create_plist, H5P_DEFAULT );
        H5Dwrite( dset_uint, H5T_NATIVE_UINT, memspace_uint, filespace_uint, write_plist, &quantities_uint[0] );
        H5Dclose( dset_uint );
        
        // Fill the vector for double quantities
        vector<double> quantities_double( n_quantities_double );
        quantities_double[ 0] = total_load                 ;
        quantities_double[ 1] = timers.global    .getTime();
        quantities_double[ 2] = timers.particles .getTime();
        quantities_double[ 3] = timers.maxwell   .getTime();
        quantities_double[ 4] = timers.densities .getTime();
        quantities_double[ 5] = timers.collisions.getTime();
        quantities_double[ 6] = timers.movWindow .getTime();
        quantities_double[ 7] = timers.loadBal   .getTime();
        quantities_double[ 8] = timers.syncPart  .getTime();
        quantities_double[ 9] = timers.syncField .getTime();
        quantities_double[10] = timers.syncDens  .getTime();
        // Specific for diags timer because we are within a diag
        double timer_diags = MPI_Wtime() - timers.diags.last_start_ + timers.diags.time_acc_;
        quantities_double[11] = timer_diags;
        quantities_double[12] = timers.grids     .getTime();
        // All timers are summed in timer_total
        double timer_total =
            quantities_double[ 2] + quantities_double[3] + quantities_double[ 4]
            + quantities_double[ 5] + quantities_double[6] + quantities_double[ 7]
            + quantities_double[ 8] + quantities_double[9] + quantities_double[10]
            + quantities_double[11] + quantities_double[12];
        quantities_double[13] = timer_total;
        
        quantities_double[14] = Tools::getMemFootPrint();
        
        // Write doubles to file
        hid_t dset_double  = H5Dcreate( iteration_group_id, "quantities_double", H5T_NATIVE_DOUBLE, filespace_double, H5P_DEFAULT, create_plist, H5P_DEFAULT );
        H5Dwrite( dset_double, H5T_NATIVE_DOUBLE, memspace_double, filespace_double, write_plist, &quantities_double[0] );
        H5Dclose( dset_double );
        
        // Patch information
        if( patch_information ) {
        
            // Creation of the group
            hid_t patch_group = H5::group( iteration_group_id, "patches" );
            
            // Prepare the hyperslab
            hsize_t matrix  = tot_number_of_patches;
            hsize_t start = vecPatches(0)->hindex;
            hsize_t count = 1;
            hsize_t block = vecPatches.size();
            hid_t filespace_patches = H5Screate_simple( 1, &matrix , NULL );
            hid_t memspace_patches  = H5Screate_simple( 1, &block, NULL );
            H5Sselect_hyperslab( filespace_patches, H5S_SELECT_SET, &start, NULL, &count, &block );
            
            // Gather x patch position in a buffer
            vector <unsigned int> buffer( number_of_patches );
            for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                buffer[ipatch] = vecPatches( ipatch )->Pcoordinates[0];
            }
            
            // Write x patch position to file
            hid_t dset_patches  = H5Dcreate( patch_group, "x", H5T_NATIVE_UINT, filespace_patches, H5P_DEFAULT, create_plist, H5P_DEFAULT );
            H5Dwrite( dset_patches, H5T_NATIVE_UINT, memspace_patches, filespace_patches, write_plist, &buffer[0] );
            H5Dclose( dset_patches );
            
            if( ndim > 1 ) {
                // Gather y patch position in a buffer
                for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                    buffer[ipatch] = vecPatches( ipatch )->Pcoordinates[1];
                }
                
                // Write y patch position to file
                dset_patches  = H5Dcreate( patch_group, "y", H5T_NATIVE_UINT, filespace_patches, H5P_DEFAULT, create_plist, H5P_DEFAULT );
                H5Dwrite( dset_patches, H5T_NATIVE_UINT, memspace_patches, filespace_patches, write_plist, &buffer[0] );
                H5Dclose( dset_patches );
            }
            
            if( ndim > 2 ) {
                // Gather z patch position in a buffer
                for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                    buffer[ipatch] = vecPatches( ipatch )->Pcoordinates[2];
                }
                
                // Write z patch position to file
                dset_patches  = H5Dcreate( patch_group, "z", H5T_NATIVE_UINT, filespace_patches, H5P_DEFAULT, create_plist, H5P_DEFAULT );
                H5Dwrite( dset_patches, H5T_NATIVE_UINT, memspace_patches, filespace_patches, write_plist, &buffer[0] );
                H5Dclose( dset_patches );
            }
            
            // Gather patch hindex in a buffer
            for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                buffer[ipatch] = vecPatches( ipatch )->hindex;
            }
            // Write patch index to file
            dset_patches  = H5Dcreate( patch_group, "index", H5T_NATIVE_UINT, filespace_patches, H5P_DEFAULT, create_plist, H5P_DEFAULT );
            H5Dwrite( dset_patches, H5T_NATIVE_UINT, memspace_patches, filespace_patches, write_plist, &buffer[0] );
            H5Dclose( dset_patches );
            
            // Creation and treatment of the species groups
            hid_t species_group;
            for( unsigned int ispecies = 0; ispecies < number_of_species; ispecies++ ) {
                species_group = H5::group( patch_group, vecPatches( 0 )->vecSpecies[ispecies]->name_ );
                
                // Vectorization properties
                if( has_adaptive_vectorization ) {
                    // Gather patch vectorization status in a buffer
                    for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                        buffer[ipatch] = ( unsigned int )( vecPatches( ipatch )->vecSpecies[ispecies]->vectorized_operators );
                    }
                    // Write patch vectorization status  to file
                    dset_patches  = H5Dcreate( species_group, "vecto", H5T_NATIVE_UINT, filespace_patches, H5P_DEFAULT, create_plist, H5P_DEFAULT );
                    H5Dwrite( dset_patches, H5T_NATIVE_UINT, memspace_patches, filespace_patches, write_plist, &buffer[0] );
                    H5Dclose( dset_patches );
                }
                
                // Close patch group
                H5Gclose( species_group );
            }
            
            // Write MPI process the owns the patch
            // Gather patch hindex in a buffer
            for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                buffer[ipatch] = mpi_rank_;
            }
            // Write mpi rank to file
            dset_patches  = H5Dcreate( patch_group, "mpi_rank", H5T_NATIVE_UINT, filespace_patches, H5P_DEFAULT, create_plist, H5P_DEFAULT );
            H5Dwrite( dset_patches, H5T_NATIVE_UINT, memspace_patches, filespace_patches, write_plist, &buffer[0] );
            H5Dclose( dset_patches );
            
            H5Sclose( filespace_patches );
            H5Sclose( memspace_patches );
            
            // Close patch group
            H5Gclose( patch_group );
        }
        
        // Close and flush
        H5Pclose( create_plist );
        H5Gclose( iteration_group_id );
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        }
    }
    
} // END run


// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticPerformances::getDiskFootPrint( int istart, int istop, Patch *patch )
{
    uint64_t footprint = 0;
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );
    
    // Add necessary global headers approximately
    footprint += 1000;
    
    // Add necessary timestep headers approximately
    footprint += ndumps * 800;
    
    // Add necessary dataset headers approximately
    footprint += ndumps * 2 * 600;
    
    // Add size of each dump
    footprint += ndumps * ( uint64_t )( mpi_size ) * ( uint64_t )( n_quantities_double * sizeof( double ) + n_quantities_uint * sizeof( unsigned int ) );
    
    return footprint;
}

void DiagnosticPerformances::setHDF5spaces( hid_t &filespace, hid_t &memspace, unsigned int height, unsigned int width, unsigned int column )
{
    hsize_t matrix [2] = {height, width};
    hsize_t portion[2] = {height, 1};
    filespace = H5Screate_simple( 2, &matrix [0], NULL );
    memspace  = H5Screate_simple( 2, &portion[0], NULL );
    hsize_t start[2] = {0, column};
    hsize_t count[2] = {1, 1};
    hsize_t block[2] = {height, 1};
    H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &start[0], NULL, &count[0], &block[0] );
}
