#include "PyTools.h"
#include <iomanip>

#include "DiagnosticPerformances.h"


using namespace std;


// Constructor
DiagnosticPerformances::DiagnosticPerformances( Params & params, SmileiMPI* smpi )
{
    fileId_ = 0;
    timestep = params.timestep;
    cell_load = params.cell_load;
    frozen_particle_load = params.frozen_particle_load;
    
    ostringstream name("");
    name << "Diagnostic performances";
    string errorPrefix = name.str();
    
    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py("every", "DiagPerformances"),
        name.str()
    );
    
    // get parameter "flush_every" which describes a timestep selection for flushing the file
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py("flush_every", "DiagPerformances"),
        name.str()
    );
    
    // Output info on diagnostics
    if ( smpi->isMaster() ) {
        MESSAGE(1,"Created performances diagnostic");
    }
    filename = "Performances.h5";
    
    // Initialize stuff for the handling of HDF5 dump
    mpi_size = smpi->getSize();
    hsize_t one = 1;
    filespace = H5Screate_simple(1, &mpi_size, NULL);
    memspace  = H5Screate_simple(1, &one, NULL );
    hsize_t start=smpi->getRank(), count=1, block=1;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, &block );
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
    
    // Calculate the number of cells per patch
    ncells_per_patch = 1;
    for (unsigned int idim = 0; idim < params.nDim_field; idim++)
        ncells_per_patch *= params.n_space[idim]+2*params.oversize[idim];
    
} // END DiagnosticPerformances::DiagnosticPerformances


DiagnosticPerformances::~DiagnosticPerformances()
{
    H5Pclose( write_plist );
    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticPerformances::~DiagnosticPerformances


// Called only by patch master of process master
void DiagnosticPerformances::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if( fileId_>0 ) return;
    
    if ( newfile ) {
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // write all parameters as HDF5 attributes
        H5::attr(fileId_, "MPI_SIZE", smpi->getSize());
    }
    else {
        // Open the existing file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose(pid);
    }
}


void DiagnosticPerformances::closeFile()
{
    if ( filespace>0 ) H5Sclose( filespace );
    if ( memspace >0 ) H5Sclose( memspace );
    if ( fileId_  >0 ) H5Fclose( fileId_ );
    fileId_ = 0;
} // END closeFile



void DiagnosticPerformances::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}


bool DiagnosticPerformances::prepare( int itime )
{
    if( timeSelection->theTimeIsNow(itime) ) {
        return true;
    } else {
        return false;
    }
} // END prepare


void DiagnosticPerformances::run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime, SimWindow* simWindow, Timers & timers )
{
    
    #pragma omp master
    {
        // Create group for this iteration
        ostringstream name_t;
        name_t.str("");
        name_t << setfill('0') << setw(10) << itime;
        status = H5Lexists(fileId_, name_t.str().c_str(), H5P_DEFAULT);
        if( status==0 )
           iteration_group_id = H5::group(fileId_, name_t.str().c_str());
        // Warning if file unreachable
        if( status < 0 ) WARNING("Performances diagnostic could not write");
    }
    #pragma omp barrier
    
    // Do not output diag if this iteration has already been written or if problem with file
    if( status != 0 ) return;
    
    #pragma omp master
    {
        // Calculate the loads
        unsigned int number_of_patches = vecPatches.size();
        unsigned int number_of_cells = ncells_per_patch * number_of_patches;
        unsigned int number_of_species = vecPatches(0)->vecSpecies.size();
        unsigned int number_of_particles=0, number_of_frozen_particles=0;
        double time = itime * timestep;
        for(unsigned int ipatch=0; ipatch < number_of_patches; ipatch++){
            for (unsigned int ispecies = 0; ispecies < number_of_species; ispecies++) {
                if( time < vecPatches(ipatch)->vecSpecies[ispecies]->time_frozen ) {
                    number_of_frozen_particles += vecPatches(ipatch)->vecSpecies[ispecies]->getNbrOfParticles();
                } else {
                    number_of_particles += vecPatches(ipatch)->vecSpecies[ispecies]->getNbrOfParticles();
                }
            }
        }
        double total_load =
            ((double)number_of_particles)
            + ((double)number_of_frozen_particles) * frozen_particle_load
            + ((double)number_of_cells) * cell_load;
        
        // Write all quantities to file
        hid_t plist_id = H5Pcreate( H5P_DATASET_CREATE );
        writeQuantity( vecPatches(0)->Hindex()    , "hindex"                    , iteration_group_id, plist_id);
        writeQuantity( number_of_cells            , "number_of_cells"           , iteration_group_id, plist_id);
        writeQuantity( number_of_particles        , "number_of_particles"       , iteration_group_id, plist_id);
        writeQuantity( number_of_frozen_particles , "number_of_frozen_particles", iteration_group_id, plist_id);
        writeQuantity( total_load                 , "total_load"                , iteration_group_id, plist_id);
        writeQuantity( timers.global    .getTime(), "timer_global"              , iteration_group_id, plist_id);
        writeQuantity( timers.particles .getTime(), "timer_particles"           , iteration_group_id, plist_id);
        writeQuantity( timers.maxwell   .getTime(), "timer_maxwell"             , iteration_group_id, plist_id);
        writeQuantity( timers.densities .getTime(), "timer_densities"           , iteration_group_id, plist_id);
        writeQuantity( timers.collisions.getTime(), "timer_collisions"          , iteration_group_id, plist_id);
        writeQuantity( timers.movWindow .getTime(), "timer_movWindow"           , iteration_group_id, plist_id);
        writeQuantity( timers.loadBal   .getTime(), "timer_loadBal"             , iteration_group_id, plist_id);
        writeQuantity( timers.syncPart  .getTime(), "timer_syncPart"            , iteration_group_id, plist_id);
        writeQuantity( timers.syncField .getTime(), "timer_syncField"           , iteration_group_id, plist_id);
        writeQuantity( timers.syncDens  .getTime(), "timer_syncDens"            , iteration_group_id, plist_id);
        double timer_total =
              timers.particles .getTime()
            + timers.maxwell   .getTime()
            + timers.densities .getTime()
            + timers.collisions.getTime()
            + timers.movWindow .getTime()
            + timers.loadBal   .getTime()
            + timers.syncPart  .getTime()
            + timers.syncField .getTime()
            + timers.syncDens  .getTime();
        // Specific for diags timer because we are within a diag
        double timer_diags = MPI_Wtime() - timers.diags.last_start_ + timers.diags.time_acc_;
        writeQuantity( timer_diags, "timer_diags", iteration_group_id, plist_id);
        // Now add the total timer
        timer_total += timer_diags;
        writeQuantity( timer_total, "timer_total", iteration_group_id, plist_id);
        
        H5Pclose(plist_id);
        H5Gclose(iteration_group_id);
        
        if( flush_timeSelection->theTimeIsNow(itime) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
    
} // END run


void DiagnosticPerformances::writeQuantity( double quantity, const char* name, hid_t gid, hid_t create_plist )
{
    hid_t dset_id  = H5Dcreate( gid, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, create_plist, H5P_DEFAULT);
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &quantity );
    H5Dclose( dset_id );
}
void DiagnosticPerformances::writeQuantity( unsigned int quantity, const char* name, hid_t gid, hid_t create_plist )
{
    hid_t dset_id  = H5Dcreate( gid, name, H5T_NATIVE_UINT, filespace, H5P_DEFAULT, create_plist, H5P_DEFAULT);
    H5Dwrite( dset_id, H5T_NATIVE_UINT, memspace, filespace, write_plist, &quantity );
    H5Dclose( dset_id );
}



// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticPerformances::getDiskFootPrint(int istart, int istop, Patch* patch)
{
    uint64_t footprint = 0;
    unsigned int n_double = 13;
    unsigned int n_uint = 4;
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore(istop) - timeSelection->howManyTimesBefore(istart);
    
    // Add necessary global headers approximately
    footprint += 1000;
    
    // Add necessary timestep headers approximately
    footprint += ndumps * 800;
    
    // Add necessary dataset headers approximately
    footprint += ndumps * (n_double + n_uint) * 350;
    
    // Add size of each dump
    footprint += ndumps * (uint64_t)(mpi_size) * (uint64_t)(n_double * sizeof(double) + n_uint * sizeof(unsigned int));
    
    return footprint;
}

