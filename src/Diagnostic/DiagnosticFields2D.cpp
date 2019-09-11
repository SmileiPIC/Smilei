
#include "DiagnosticFields2D.h"

#include <sstream>
#include <cmath>
#include <limits>

#include "Params.h"
#include "Patch.h"
#include "Field2D.h"
#include "VectorPatch.h"
#include "DomainDecomposition.h"
#include "Hilbert_functions.h"
#include "LinearizedDomainDecomposition.h"

using namespace std;

DiagnosticFields2D::DiagnosticFields2D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &openPMD )
    : DiagnosticFields( params, smpi, vecPatches, ndiag, openPMD )
{

    // Calculate the offset in the local grid
    patch_offset_in_grid.resize( 2 );
    patch_offset_in_grid[0] = params.oversize[0]+1;
    patch_offset_in_grid[1] = params.oversize[1]+1;
    
    // Calculate the patch size
    patch_size.resize( 2 );
    patch_size[0] = params.n_space[0];
    patch_size[1] = params.n_space[1];
    
    // define space in file and in memory for the first (1D) write
    // We assign, for each patch, the maximum buffer size necessary to fit the subgrid
    unsigned int istart_in_patch[2], istart_in_file[2], nsteps[2];
    for( unsigned int i=0; i<2; i++ ) {
        findSubgridIntersection(
            subgrid_start_[i], subgrid_stop_[i], subgrid_step_[i],
            subgrid_start_[i], subgrid_start_[i]+patch_size[i]+1,
            istart_in_patch[i], istart_in_file[i], nsteps[i]
        );
    }
    one_patch_buffer_size = nsteps[0] * nsteps[1];
    hsize_t file_size = ( hsize_t )one_patch_buffer_size * ( hsize_t )tot_number_of_patches;
    filespace_firstwrite = H5Screate_simple( 1, &file_size, NULL );
    memspace_firstwrite  = H5Screate_simple( 1, &file_size, NULL );
    
    if( smpi->test_mode ) {
        return;
    }
    
    // Define a second portion of the grid, which is unrelated to the current
    // composition of vecPatches. It is used for a second writing of the file
    // in order to fold the Hilbert curve. This new portion is necessarily
    // rectangular for efficient writing.
    int nproc = smpi->getSize(), iproc = smpi->getRank();
    int npatch = params.tot_number_of_patches;
    int npatch_local = 1<<int( log2( ( ( double )npatch )/nproc ) );
    if( dynamic_cast<LinearizedDomainDecomposition2D *>( vecPatches.domain_decomposition_ ) ) {
        npatch_local = vecPatches.size();
    }
    int first_proc_with_less_patches = ( npatch-npatch_local*nproc )/npatch_local;
    int first_patch_of_this_proc;
    if( iproc < first_proc_with_less_patches ) {
        npatch_local *= 2;
        first_patch_of_this_proc = npatch_local*iproc;
    } else {
        first_patch_of_this_proc = npatch_local*( first_proc_with_less_patches+iproc );
    }
    // Define space in file for re-reading
    filespace_reread = H5Screate_simple( 1, &file_size, NULL );
    hsize_t offset = ( hsize_t )one_patch_buffer_size * ( hsize_t )first_patch_of_this_proc;
    hsize_t block  = ( hsize_t )one_patch_buffer_size * ( hsize_t )npatch_local;
    hsize_t count  = 1;
    H5Sselect_hyperslab( filespace_reread, H5S_SELECT_SET, &offset, NULL, &count, &block );
    // Define space in memory for re-reading
    memspace_reread = H5Screate_simple( 1, &block, NULL );
    data_reread.resize( block );
    // Define the list of patches for re-writing
    rewrite_npatch = ( unsigned int )npatch_local;
    rewrite_patch.resize( rewrite_npatch );
    rewrite_xmin=numeric_limits<int>::max();
    rewrite_ymin=numeric_limits<int>::max();
    unsigned int rewrite_xmax=0, rewrite_ymax=0;
    for( unsigned int h=0; h<rewrite_npatch; h++ ) {
        std::vector<unsigned int> xcall( 2, 0 );
        xcall = vecPatches.domain_decomposition_->getDomainCoordinates( first_patch_of_this_proc+h );
        if( xcall[0]<rewrite_xmin ) {
            rewrite_xmin=xcall[0];
        }
        if( xcall[0]>rewrite_xmax ) {
            rewrite_xmax=xcall[0];
        }
        if( xcall[1]<rewrite_ymin ) {
            rewrite_ymin=xcall[1];
        }
        if( xcall[1]>rewrite_ymax ) {
            rewrite_ymax=xcall[1];
        }
        rewrite_patch[h] = xcall;
    }
    rewrite_npatchx = rewrite_xmax - rewrite_xmin + 1;
    rewrite_npatchy = rewrite_ymax - rewrite_ymin + 1;
    // Define space in file for re-writing
    hsize_t final_array_size[2], offset2[2], block2[2], count2[2];
    final_array_size[0] = params.number_of_patches[0] * params.n_space[0] + 1;
    final_array_size[1] = params.number_of_patches[1] * params.n_space[1] + 1;
    offset2[0] = rewrite_xmin * params.n_space[0] + ( ( rewrite_xmin==0 )?0:1 );
    offset2[1] = rewrite_ymin * params.n_space[1] + ( ( rewrite_ymin==0 )?0:1 );
    block2 [0] = rewrite_npatchx * params.n_space[0] + ( ( rewrite_xmin==0 )?1:0 );
    block2 [1] = rewrite_npatchy * params.n_space[1] + ( ( rewrite_ymin==0 )?1:0 );
    count2 [0] = 1;
    count2 [1] = 1;
    // Take subgrid into account
    unsigned int istart[2];
    total_dataset_size = 1;
    for( unsigned int i=0; i<2; i++ ) {
        findSubgridIntersection(
            subgrid_start_[i], subgrid_stop_[i], subgrid_step_[i],
            0, final_array_size[i],
            istart[i], rewrite_start_in_file[i], rewrite_size[i]
        );
        final_array_size[i] = rewrite_size[i];
        total_dataset_size *= final_array_size[i];
        findSubgridIntersection(
            subgrid_start_[i], subgrid_stop_[i], subgrid_step_[i],
            offset2[i], offset2[i] + block2[i],
            istart[i], rewrite_start_in_file[i], rewrite_size[i]
        );
        offset2[i] = rewrite_start_in_file[i];
        block2 [i] = rewrite_size[i];
    }
    filespace = H5Screate_simple( 2, final_array_size, NULL );
    if( rewrite_size[0]==0 || rewrite_size[1]==0 ) {
        H5Sselect_none( filespace );
    } else {
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset2, NULL, count2, block2 );
    }
    // Define space in memory for re-writing
    memspace = H5Screate_simple( 2, block2, NULL );
    data_rewrite.resize( rewrite_size[0]*rewrite_size[1] );
    
    // Define the chunk size (necessary above 2^28 points)
    const hsize_t max_size = 4294967295/2/sizeof( double );
    // For the first write
    dcreate_firstwrite = H5Pcreate( H5P_DATASET_CREATE );
    if( file_size > max_size ) {
        hsize_t n_chunks = 1 + ( file_size-1 ) / max_size;
        hsize_t chunk_size = file_size / n_chunks;
        if( n_chunks * chunk_size < file_size ) {
            chunk_size++;
        }
        H5Pset_layout( dcreate_firstwrite, H5D_CHUNKED );
        H5Pset_chunk( dcreate_firstwrite, 1, &chunk_size );
    }
    // For the second write
    hsize_t final_size = final_array_size[0]
                         *final_array_size[1];
    if( final_size > max_size ) {
        hsize_t n_chunks = 1 + ( final_size-1 ) / max_size;
        hsize_t chunk_size[2];
        chunk_size[0] = final_array_size[0] / n_chunks;
        chunk_size[1] = final_array_size[1];
        if( n_chunks * chunk_size[0] < final_array_size[0] ) {
            chunk_size[0]++;
        }
        H5Pset_layout( dcreate, H5D_CHUNKED );
        H5Pset_chunk( dcreate, 2, chunk_size );
    }
    
    tmp_dset_id=0;
}

DiagnosticFields2D::~DiagnosticFields2D()
{
    H5Pclose( dcreate_firstwrite );
}


void DiagnosticFields2D::setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int buffer_size = one_patch_buffer_size * vecPatches.size();
    
    // Resize the data
    data.resize( buffer_size );
    
    // Define offset and size for HDF5 file
    hsize_t offset = one_patch_buffer_size * refHindex;
    hsize_t block  = buffer_size;
    hsize_t count  = 1;
    // Select portion of the file where this MPI will write to
    H5Sselect_hyperslab( filespace_firstwrite, H5S_SELECT_SET, &offset, NULL, &count, &block );
    // define space in memory
    H5Sset_extent_simple( memspace_firstwrite, 1, &block, &block );
    
    // Create/Open temporary dataset
    status = H5Lexists( fileId_, "tmp", H5P_DEFAULT );
    if( status == 0 ) {
        tmp_dset_id  = H5Dcreate( fileId_, "tmp", H5T_NATIVE_DOUBLE, filespace_firstwrite, H5P_DEFAULT, dcreate_firstwrite, H5P_DEFAULT );
    } else {
        hid_t pid = H5Pcreate( H5P_DATASET_ACCESS );
        tmp_dset_id = H5Dopen( fileId_, "tmp", pid );
        H5Pclose( pid );
    }
}


// Copy patch field to current "data" buffer
void DiagnosticFields2D::getField( Patch *patch, unsigned int ifield )
{
    // Get current field
    Field2D *field;
    if( time_average>1 ) {
        field = static_cast<Field2D *>( patch->EMfields->allFields_avg[diag_n][ifield] );
    } else {
        field = static_cast<Field2D *>( patch->EMfields->allFields[fields_indexes[ifield]] );
    }
    
    // Find the intersection between this patch and the subgrid
    unsigned int istart_in_patch[2], istart_in_file[2], nsteps[2], patch_begin[2], patch_end[2];
    for( unsigned int i=0; i<2; i++ ) {
        patch_begin[i] = patch->Pcoordinates[i] * patch_size[i];
        patch_end  [i] = patch_begin[i] + patch_size[i] + 1;
        if( patch->Pcoordinates[i] != 0 ) {
            patch_begin[i]++;
        }
        findSubgridIntersection(
            subgrid_start_[i], subgrid_stop_[i], subgrid_step_[i],
            patch_begin[i], patch_end[i],
            istart_in_patch[i], istart_in_file[i], nsteps[i]
        );
        istart_in_patch[i] += patch_offset_in_grid[i];
        if( patch->Pcoordinates[i] == 0 ) {
            istart_in_patch[i]--;
        }
    }
    
    // Copy field to the "data" buffer
    unsigned int ix_max = istart_in_patch[0] + subgrid_step_[0]*nsteps[0];
    unsigned int iy_max = istart_in_patch[1] + subgrid_step_[1]*nsteps[1];
    unsigned int iout = one_patch_buffer_size * ( patch->Hindex()-refHindex );
    for( unsigned int ix = istart_in_patch[0]; ix < ix_max; ix += subgrid_step_[0] ) {
        for( unsigned int iy = istart_in_patch[1]; iy < iy_max; iy += subgrid_step_[1] ) {
            data[iout] = ( *field )( ix, iy ) * time_average_inv;
            iout++;
        }
    }
    
    if( time_average>1 ) {
        field->put_to( 0.0 );
    }
}


// Write current buffer to file
void DiagnosticFields2D::writeField( hid_t dset_id, int itime )
{

    // Write the buffer in a temporary location
    H5Dwrite( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_firstwrite, filespace_firstwrite, write_plist, &( data[0] ) );
    
    // Read the file with the previously defined partition
    H5Dread( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_reread, filespace_reread, write_plist, &( data_reread[0] ) );
    
    // Fold the data according to the Hilbert curve
    unsigned int read_position, write_position, write_skip;
    unsigned int istart_in_patch[2], istart_in_file[2], nsteps[2], patch_begin[2], patch_end[2];
    
    for( unsigned int h=0; h<rewrite_npatch; h++ ) {
    
        for( unsigned int i=0; i<2; i++ ) {
            patch_begin[i] = rewrite_patch[h][i] * patch_size[i];
            patch_end  [i] = patch_begin[i] + patch_size[i] + 1;
            if( patch_begin[i] != 0 ) {
                patch_begin[i]++;
            }
            findSubgridIntersection(
                subgrid_start_[i], subgrid_stop_[i], subgrid_step_[i],
                patch_begin[i], patch_end[i],
                istart_in_patch[i], istart_in_file[i], nsteps[i]
            );
            istart_in_file[i] -= rewrite_start_in_file[i];
        }
        
        read_position = one_patch_buffer_size * h;
        write_position = istart_in_file[1] + istart_in_file[0] * rewrite_size[1];
        write_skip = rewrite_size[1] - nsteps[1];
        for( unsigned int ix=0; ix<nsteps[0]; ix++ ) {
            for( unsigned int iy=0; iy<nsteps[1]; iy++ ) {
                data_rewrite[write_position] = data_reread[read_position];
                read_position ++;
                write_position++;
            }
            write_position += write_skip;
        }
        
    }
    
    // Rewrite the file with the previously defined partition
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &( data_rewrite[0] ) );
    
}

