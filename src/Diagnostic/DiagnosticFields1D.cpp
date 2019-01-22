
#include "DiagnosticFields1D.h"

#include "Params.h"
#include "Patch.h"
#include "Field1D.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields1D::DiagnosticFields1D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &openPMD )
    : DiagnosticFields( params, smpi, vecPatches, ndiag, openPMD )
{
    // Calculate the offset in the local grid
    patch_offset_in_grid.resize( 1 );
    patch_offset_in_grid[0] = params.oversize[0]+1;
    
    // Calculate the patch size
    total_patch_size = params.n_space[0];
    
    // define space in file and in memory
    // All patch write n_space elements except, patch 0 which write n_space+1
    unsigned int global_size = tot_number_of_patches * total_patch_size + 1;
    // Take subgrid into account
    unsigned int istart, istart_in_file, nsteps;
    findSubgridIntersection(
        subgrid_start_[0], subgrid_stop_[0], subgrid_step_[0],
        0, global_size,
        istart, istart_in_file, nsteps
    );
    hsize_t file_size = nsteps;
    one_patch_buffer_size = nsteps;
    filespace = H5Screate_simple( 1, &file_size, NULL );
    memspace  = H5Screate_simple( 1, &file_size, NULL );
}

DiagnosticFields1D::~DiagnosticFields1D()
{
}

void DiagnosticFields1D::setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int total_vecPatches_size = total_patch_size * vecPatches.size();
    // One more cell on the left
    if( smpi->isMaster() ) {
        total_vecPatches_size++;
    }
    
    // Calculate the intersection between the local grid and the subgrid
    unsigned int MPI_begin = total_patch_size * refHindex;
    if( !smpi->isMaster() ) {
        MPI_begin++;
    }
    unsigned int MPI_end = MPI_begin + total_vecPatches_size;
    unsigned int istart_in_MPI, nsteps;
    findSubgridIntersection(
        subgrid_start_[0], subgrid_stop_[0], subgrid_step_[0],
        MPI_begin, MPI_end,
        istart_in_MPI, MPI_start_in_file, nsteps
    );
    
    if( nsteps > 0 ) {
        data.resize( nsteps );
        
        // Define offset and size for HDF5 file
        hsize_t offset[1], block[1], count[1];
        offset[0] = MPI_start_in_file;
        block[0] = nsteps;
        count[0] = 1;
        // Select portion of the file where this MPI will write to
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        
        // define space in memory
        offset[0] = 0;
        H5Sselect_hyperslab( memspace, H5S_SELECT_SET, offset, NULL, count, block );
    } else {
        data.resize( 0 );
        H5Sselect_none( filespace );
        H5Sselect_none( memspace );
    }
}


// Copy patch field to current "data" buffer
void DiagnosticFields1D::getField( Patch *patch, unsigned int ifield )
{
    // Get current field
    Field1D *field;
    if( time_average>1 ) {
        field = static_cast<Field1D *>( patch->EMfields->allFields_avg[diag_n][ifield] );
    } else {
        field = static_cast<Field1D *>( patch->EMfields->allFields[fields_indexes[ifield]] );
    }
    
    // Calculate the intersection between the patch grid and the subgrid
    unsigned int patch_begin = total_patch_size * patch->Hindex();
    unsigned int patch_end   = patch_begin + total_patch_size + 1;
    unsigned int ix, iout, nsteps;
    if( patch->Hindex() != 0 ) {
        patch_begin++;
    }
    findSubgridIntersection(
        subgrid_start_[0], subgrid_stop_[0], subgrid_step_[0],
        patch_begin, patch_end,
        ix, iout, nsteps
    );
    ix += patch_offset_in_grid[0];
    if( patch->Hindex() == 0 ) {
        ix--;
    }
    iout -= MPI_start_in_file;
    unsigned int ix_max = ix + nsteps * subgrid_step_[0];
    
    // Copy this patch field into buffer
    while( ix < ix_max ) {
        data[iout] = ( *field )( ix ) * time_average_inv;
        ix += subgrid_step_[0];
        iout++;
    }
    
    if( time_average>1 ) {
        field->put_to( 0.0 );
    }
}


// Write current buffer to file
void DiagnosticFields1D::writeField( hid_t dset_id, int itime )
{

    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &( data[0] ) );
    
}

