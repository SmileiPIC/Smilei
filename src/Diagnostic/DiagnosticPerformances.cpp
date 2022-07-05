#include "PyTools.h"
#include <iomanip>

#include "DiagnosticPerformances.h"


using namespace std;

const unsigned int n_quantities_double = 19;
const unsigned int n_quantities_uint   = 4;

// Constructor
DiagnosticPerformances::DiagnosticPerformances( Params &params, SmileiMPI *smpi )
: mpi_size_( smpi->getSize() ),
  mpi_rank_( smpi->getRank() ),
  filespace_double( {n_quantities_double, mpi_size_}, {0, mpi_rank_}, {n_quantities_double, 1} ),
  filespace_uint  ( {n_quantities_uint  , mpi_size_}, {0, mpi_rank_}, {n_quantities_uint  , 1} ),
  memspace_double( { n_quantities_double, 1 }, {}, {} ),
  memspace_uint  ( { n_quantities_uint  , 1 }, {}, {} )
{
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
    
    ndim     = params.nDim_field;
    has_adaptive_vectorization = params.has_adaptive_vectorization;
    
    // Calculate the number of cells per patch
    ncells_per_patch = 1;
    for( unsigned int idim = 0; idim < params.nDim_field; idim++ ) {
        ncells_per_patch *= params.n_space[idim]+2*params.oversize[idim];
    }
    
} // END DiagnosticPerformances::DiagnosticPerformances


DiagnosticPerformances::~DiagnosticPerformances()
{
    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticPerformances::~DiagnosticPerformances


// Called only by patch master of process master
void DiagnosticPerformances::openFile( Params &params, SmileiMPI *smpi )
{
    if( file_ ) {
        return;
    }
    
    file_ = new H5Write( filename, &smpi->world() );
    
    // write all parameters as HDF5 attributes
    file_->attr( "MPI_SIZE", smpi->getSize() );
    file_->attr( "patch_arrangement", params.patch_arrangement );
    
    vector<string> quantities_uint( n_quantities_uint );
    quantities_uint[0] = "hindex"                    ;
    quantities_uint[1] = "number_of_cells"           ;
    quantities_uint[2] = "number_of_particles"       ;
    quantities_uint[3] = "number_of_frozen_particles";
    file_->attr( "quantities_uint", quantities_uint );
    
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
    quantities_double[15] = "memory_peak"    ;
    quantities_double[16] = "timer_envelope"     ;
    quantities_double[17] = "timer_syncSusceptibility"     ;
    quantities_double[18] = "timer_partMerging"     ;
    file_->attr( "quantities_double", quantities_double );
    
    file_->flush();
}


void DiagnosticPerformances::closeFile()
{
    if( file_ ) {
        delete file_;
        file_ = NULL;
    }
} // END closeFile



void DiagnosticPerformances::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // create the file
    openFile( params, smpi );
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
        group_name = name_t.str();
        has_group = file_->has( group_name );
    }
    #pragma omp barrier
    
    // Do not output diag if this iteration has already been written or if problem with file
    if( has_group ) {
        return;
    }
    
    #pragma omp master
    {
        H5Write iteration_group = file_->group( group_name );
        
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
        
        // Fill the vector for uint quantities
        vector<unsigned int> quantities_uint( n_quantities_uint );
        quantities_uint[0] = vecPatches( 0 )->Hindex()   ;
        quantities_uint[1] = number_of_cells           ;
        quantities_uint[2] = number_of_particles       ;
        quantities_uint[3] = number_of_frozen_particles;
        
        // Write uints to file
        iteration_group.array( "quantities_uint", quantities_uint[0], &filespace_uint, &memspace_uint );
        
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
        
        quantities_double[14] = Tools::getMemFootPrint(0);
        quantities_double[15] = Tools::getMemFootPrint(1);
        quantities_double[16] = timers.envelope         .getTime();
        quantities_double[17] = timers.susceptibility   .getTime();
        quantities_double[18] = timers.particleMerging  .getTime();
        
        // Write doubles to file
        iteration_group.array( "quantities_double", quantities_double[0], &filespace_double, &memspace_double );
        
        // Patch information
        if( patch_information ) {
        
            // Creation of the group
            H5Write patch_group = iteration_group.group( "patches" );
            
            // Prepare the hyperslab
            hsize_t size = tot_number_of_patches;
            hsize_t offset = vecPatches(0)->hindex;
            hsize_t npoints = vecPatches.size();
            
            // Gather x patch position in a buffer
            vector <unsigned int> buffer( number_of_patches );
            for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                buffer[ipatch] = vecPatches( ipatch )->Pcoordinates[0];
            }
            // Write x patch position to file
            patch_group.vect( "x", buffer[0], size, H5T_NATIVE_UINT, offset, npoints );
            
            if( ndim > 1 ) {
                // Gather y patch position in a buffer
                for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                    buffer[ipatch] = vecPatches( ipatch )->Pcoordinates[1];
                }
                // Write y patch position to file
                patch_group.vect( "y", buffer[0], size, H5T_NATIVE_UINT, offset, npoints );
            }
            
            if( ndim > 2 ) {
                // Gather z patch position in a buffer
                for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                    buffer[ipatch] = vecPatches( ipatch )->Pcoordinates[2];
                }
                // Write z patch position to file
                patch_group.vect( "z", buffer[0], size, H5T_NATIVE_UINT, offset, npoints );
            }
            
            // Gather patch hindex in a buffer
            for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                buffer[ipatch] = vecPatches( ipatch )->hindex;
            }
            // Write patch index to file
            patch_group.vect( "index", buffer[0], size, H5T_NATIVE_UINT, offset, npoints );
            
            // Creation and treatment of the species groups
            for( unsigned int ispecies = 0; ispecies < number_of_species; ispecies++ ) {
                H5Write species_group = patch_group.group( vecPatches( 0 )->vecSpecies[ispecies]->name_ );
                
                // Vectorization properties
                if( has_adaptive_vectorization ) {
                    // Gather patch vectorization status in a buffer
                    for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                        buffer[ipatch] = ( unsigned int )( vecPatches( ipatch )->vecSpecies[ispecies]->vectorized_operators );
                    }
                    // Write patch vectorization status  to file
                    species_group.vect( "vecto", buffer[0], size, H5T_NATIVE_UINT, offset, npoints );
                }
            }
            
            // Write MPI process the owns the patch
            // Gather patch hindex in a buffer
            for( unsigned int ipatch=0; ipatch < number_of_patches; ipatch++ ) {
                buffer[ipatch] = mpi_rank_;
            }
            // Write mpi rank to file
            patch_group.vect( "mpi_rank", buffer[0], size, H5T_NATIVE_UINT, offset, npoints );
            
        }
        
        // Close and flush
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            file_->flush();
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
    footprint += ndumps * ( uint64_t )( mpi_size_ ) * ( uint64_t )( n_quantities_double * sizeof( double ) + n_quantities_uint * sizeof( unsigned int ) );
    
    return footprint;
}
