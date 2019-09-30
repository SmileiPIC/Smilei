
#include "DiagnosticFieldsAM.h"

#include <sstream>
#include <cmath>
#include <limits>

#include "Params.h"
#include "Patch.h"
#include "cField2D.h"
#include "VectorPatch.h"
#include "Hilbert_functions.h"

using namespace std;

DiagnosticFieldsAM::DiagnosticFieldsAM( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &openPMD )
    : DiagnosticFields( params, smpi, vecPatches, ndiag, openPMD )
{
    factor_ = 2;
    // Be able to dump Env_Chi

    // Default : AM diag fields are complex
    vector<bool> isReal;
    // 1. Check if real or complex :
    Field* field = NULL;
    for( unsigned int ifield=0; ifield < fields_indexes.size(); ifield++ ) {
        if( time_average>1 ) {
            field = vecPatches(0)->EMfields->allFields_avg[diag_n][ifield];
        } else {
            field = vecPatches(0)->EMfields->allFields[fields_indexes[ifield]];
        }
        if ( !dynamic_cast<cField*>( field ) )
            isReal.push_back( true );
        else
            isReal.push_back( false );
    }
    // 2. convert to double :
    int nRealFields = std::count( isReal.begin(), isReal.end(), true );
    if ( nRealFields == fields_indexes.size() )
        factor_ /= 2;
    else if ( nRealFields != 0 ) {
        ERROR( "Mixing types in diag field : " << ndiag );
    }

    // Calculate the offset in the local grid
    patch_offset_in_grid.resize( 2 );
    patch_offset_in_grid[0] = params.oversize[0];
    patch_offset_in_grid[1] = params.oversize[1];
    
    // Calculate the patch size
    patch_size.resize( 2 );
    patch_size[0] = params.n_space[0] + 1;
    patch_size[1] = params.n_space[1] + 1;
    one_patch_buffer_size = patch_size[0] * patch_size[1];
    
    
    // define space in file
    hsize_t global_size[1];
    global_size[0] = tot_number_of_patches * one_patch_buffer_size;
    hsize_t iglobal_size[1];
    iglobal_size[0] = factor_ * global_size[0];
    filespace_firstwrite = H5Screate_simple( 1, iglobal_size, NULL );
    memspace_firstwrite = H5Screate_simple( 1, iglobal_size, NULL ); // redefined later
    
    // Define a second subset of the grid, which is unrelated to the current
    // composition of vecPatches. It is used for a second writing of the file
    // in order to fold the Hilbert curve. This new subset is necessarily
    // rectangular for efficient writing.
    hsize_t offset[1], block[1], count[1];
    hsize_t ioffset[1], iblock[1];
    int nproc = smpi->getSize(), iproc = smpi->getRank();
    int npatch = params.tot_number_of_patches;
    int npatch_local = 1<<int( log2( ( ( double )npatch )/nproc ) );
    int first_proc_with_less_patches = ( npatch-npatch_local*nproc )/npatch_local;
    int first_patch_of_this_proc;
    if( iproc < first_proc_with_less_patches ) {
        npatch_local *= 2;
        first_patch_of_this_proc = npatch_local*iproc;
    } else {
        first_patch_of_this_proc = npatch_local*( first_proc_with_less_patches+iproc );
    }
    // Define space in file for re-reading
    filespace_reread = H5Screate_simple( 1, iglobal_size, NULL );
    offset[0] = one_patch_buffer_size * first_patch_of_this_proc;
    block [0] = one_patch_buffer_size * npatch_local;
    ioffset[0] = factor_ * offset[0];
    iblock [0] = factor_ * block [0];
    count [0] = 1;
    H5Sselect_hyperslab( filespace_reread, H5S_SELECT_SET, ioffset, NULL, count, iblock );
    // Define space in memory for re-reading
    memspace_reread = H5Screate_simple( 1, iblock, NULL );
    if (factor_==2)
        idata_reread.resize( block[0] );
    else if (factor_==1)
        data_reread.resize( block[0] );

    // Define the list of patches for re-writing
    rewrite_npatch = ( unsigned int )npatch_local;
    rewrite_patches_x.resize( rewrite_npatch );
    rewrite_patches_y.resize( rewrite_npatch );
    rewrite_xmin=numeric_limits<int>::max();
    rewrite_ymin=numeric_limits<int>::max();
    unsigned int rewrite_xmax=0, rewrite_ymax=0, x, y;
    for( unsigned int h=0; h<rewrite_npatch; h++ ) {
        generalhilbertindexinv( params.mi[0],  params.mi[1], &x, &y, first_patch_of_this_proc+h );
        if( x<rewrite_xmin ) {
            rewrite_xmin=x;
        }
        if( x>rewrite_xmax ) {
            rewrite_xmax=x;
        }
        if( y<rewrite_ymin ) {
            rewrite_ymin=y;
        }
        if( y>rewrite_ymax ) {
            rewrite_ymax=y;
        }
        rewrite_patches_x[h] = x;
        rewrite_patches_y[h] = y;
    }
    rewrite_npatchx = rewrite_xmax - rewrite_xmin + 1;
    rewrite_npatchy = rewrite_ymax - rewrite_ymin + 1;
    // Define space in file for re-writing
    hsize_t final_array_size[2], offset2[2], block2[2], count2[2];
    hsize_t ifinal_array_size[2], ioffset2[2], iblock2[2];
    final_array_size[0] = params.number_of_patches[0] * params.n_space[0] + 1;
    final_array_size[1] = params.number_of_patches[1] * params.n_space[1] + 1;
    ifinal_array_size[0] =     final_array_size[0];
    ifinal_array_size[1] = factor_ * final_array_size[1];
    total_dataset_size = ifinal_array_size[0] * ifinal_array_size[1];
    filespace = H5Screate_simple( 2, ifinal_array_size, NULL );
    offset2[0] = rewrite_xmin * params.n_space[0] + ( ( rewrite_xmin==0 )?0:1 );
    offset2[1] = rewrite_ymin * params.n_space[1] + ( ( rewrite_ymin==0 )?0:1 );
    block2 [0] = rewrite_npatchx * params.n_space[0] + ( ( rewrite_xmin==0 )?1:0 );
    block2 [1] = rewrite_npatchy * params.n_space[1] + ( ( rewrite_ymin==0 )?1:0 );
    ioffset2[0] = offset2[0];
    ioffset2[1] = factor_ * offset2[1];
    iblock2 [0] = block2 [0];
    iblock2 [1] = factor_ * block2 [1];
    count2 [0] = 1;
    count2 [1] = 1;
    H5Sselect_hyperslab( filespace, H5S_SELECT_SET, ioffset2, NULL, count2, iblock2 );
    // Define space in memory for re-writing
    memspace = H5Screate_simple( 2, iblock2, NULL );
    if (factor_==2)
        idata_rewrite.resize( block2[0]*block2[1] );
    else if (factor_==1)
        data_rewrite.resize( block2[0]*block2[1] );
    
    tmp_dset_id=0;
}

DiagnosticFieldsAM::~DiagnosticFieldsAM()
{
}


void DiagnosticFieldsAM::setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int total_vecPatches_size = one_patch_buffer_size * vecPatches.size();
    
    // Resize the data
    if (factor_==2)
        idata.resize( total_vecPatches_size );
    else if (factor_==1)
        data.resize( total_vecPatches_size );
    
    // Define offset and size for HDF5 file
    hsize_t offset[1], block[1], count[1];
    offset[0] = one_patch_buffer_size * refHindex;
    block [0] = total_vecPatches_size;
    count [0] = 1;
    hsize_t ioffset[1], iblock[1];
    ioffset[0] = factor_ * offset[0];
    iblock [0] = factor_ * block [0];
    // Select portion of the file where this MPI will write to
    H5Sselect_hyperslab( filespace_firstwrite, H5S_SELECT_SET, ioffset, NULL, count, iblock );
    // define space in memory
    H5Sset_extent_simple( memspace_firstwrite, 1, iblock, iblock );
    
    // Create/Open temporary dataset
    status = H5Lexists( fileId_, "tmp", H5P_DEFAULT );
    if( status == 0 ) {
        hid_t pid = H5Pcreate( H5P_DATASET_CREATE );
        tmp_dset_id  = H5Dcreate( fileId_, "tmp", H5T_NATIVE_DOUBLE, filespace_firstwrite, H5P_DEFAULT, pid, H5P_DEFAULT );
        H5Pclose( pid );
    } else {
        hid_t pid = H5Pcreate( H5P_DATASET_ACCESS );
        tmp_dset_id = H5Dopen( fileId_, "tmp", pid );
        H5Pclose( pid );
    }
}


// Copy patch field to current "data" buffer
void DiagnosticFieldsAM::getField( Patch *patch, unsigned int ifield )
{
    // Copy field to the "data" buffer
    if (factor_==2) { // Complex
        // Get current field
        cField2D *field;
        if( time_average>1 ) {
            field = static_cast<cField2D *>( patch->EMfields->allFields_avg[diag_n][ifield] );
        } else {
            field = static_cast<cField2D *>( patch->EMfields->allFields[fields_indexes[ifield]] );
        }

        unsigned int ix = patch_offset_in_grid[0];
        unsigned int ix_max = ix + patch_size[0];
        unsigned int iy;
        unsigned int iy_max = patch_offset_in_grid[1] + patch_size[1];
        unsigned int iout = one_patch_buffer_size * ( patch->Hindex()-refHindex );
        while( ix < ix_max ) {
            iy = patch_offset_in_grid[1];
            while( iy < iy_max ) {
                idata[iout] = ( *field )( ix, iy ) * time_average_inv;
                iout++;
                iy++;
            }
            ix++;
        }
        if( time_average>1 ) {
            field->put_to( 0.0 );
        }
    }
    else {  // Double
        // Get current field
        Field2D *field;
        if( time_average>1 ) {
            field = static_cast<Field2D *>( patch->EMfields->allFields_avg[diag_n][ifield] );
        } else {
            field = static_cast<Field2D *>( patch->EMfields->allFields[fields_indexes[ifield]] );
        }
        unsigned int ix = patch_offset_in_grid[0];
        unsigned int ix_max = ix + patch_size[0];
        unsigned int iy;
        unsigned int iy_max = patch_offset_in_grid[1] + patch_size[1];
        unsigned int iout = one_patch_buffer_size * ( patch->Hindex()-refHindex );
        while( ix < ix_max ) {
            iy = patch_offset_in_grid[1];
            while( iy < iy_max ) {
                data[iout] = ( *field )( ix, iy ) * time_average_inv;
                iout++;
                iy++;
            }
            ix++;
        }
        if( time_average>1 ) {
            field->put_to( 0.0 );
        }
    }

}


// Write current buffer to file
void DiagnosticFieldsAM::writeField( hid_t dset_id, int itime )
{

    if (factor_==2) {
        // Write the buffer in a temporary location
        H5Dwrite( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_firstwrite, filespace_firstwrite, write_plist, &( idata[0] ) );

        // Read the file with the previously defined partition
        H5Dread( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_reread, filespace_reread, write_plist, &( idata_reread[0] ) );
    
        // Fold the data according to the Hilbert curve
        unsigned int read_position, write_position, write_skip_y, sx, sy;
    
        unsigned int write_sizey  = ( rewrite_npatchy*( patch_size[1]-1 ) + ( ( rewrite_ymin==0 )?1:0 ) );
    
        read_position = 0;
        for( unsigned int h=0; h<rewrite_npatch; h++ ) {
            int write_position0 = ( rewrite_patches_y[h]-rewrite_ymin )*( patch_size[1]-1 )
                + write_sizey *( ( rewrite_patches_x[h]-rewrite_xmin )*( patch_size[0]-1 ) );
                              
            write_skip_y = ( rewrite_npatchy - 1 )*( patch_size[1]-1 );
        
            sx = patch_size[0];
            sy = patch_size[1];
        
            if( rewrite_patches_y[h]!=0 ) {
                if( rewrite_ymin==0 ) {
                    write_position0++;
                    write_skip_y++;
                }
                sy--;
            }
            if( rewrite_patches_x[h]!=0 ) {
                read_position += patch_size[1];
                if( rewrite_xmin==0 ) {
                    write_position0 += write_sizey;
                }
                sx--;
            }

            write_position = write_position0;
            for( unsigned int ix=0; ix<sx; ix++ ) {
                if( rewrite_patches_y[h]!=0 ) {
                    read_position ++;
                }
                for( unsigned int iy=0; iy<sy; iy++ ) {
                    //data_rewrite[write_position] += h+1;
                    idata_rewrite[write_position] = idata_reread[read_position];
                    read_position ++;
                    write_position++;

                }
                write_position += write_skip_y;
            }
        }
        
        // Rewrite the file with the previously defined partition
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &( idata_rewrite[0] ) );
    }
    else if (factor_==1) {
        // Write the buffer in a temporary location
        H5Dwrite( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_firstwrite, filespace_firstwrite, write_plist, &( data[0] ) );

        // Read the file with the previously defined partition
        H5Dread( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_reread, filespace_reread, write_plist, &( data_reread[0] ) );

        // Fold the data according to the Hilbert curve
        unsigned int read_position, write_position, write_skip_y, sx, sy;

        unsigned int write_sizey  = ( rewrite_npatchy*( patch_size[1]-1 ) + ( ( rewrite_ymin==0 )?1:0 ) );

        read_position = 0;
        for( unsigned int h=0; h<rewrite_npatch; h++ ) {
            int write_position0 = ( rewrite_patches_y[h]-rewrite_ymin )*( patch_size[1]-1 )
                + write_sizey *( ( rewrite_patches_x[h]-rewrite_xmin )*( patch_size[0]-1 ) );

            write_skip_y = ( rewrite_npatchy - 1 )*( patch_size[1]-1 );

            sx = patch_size[0];
            sy = patch_size[1];

            if( rewrite_patches_y[h]!=0 ) {
                if( rewrite_ymin==0 ) {
                    write_position0++;
                    write_skip_y++;
                }
                sy--;
            }
            if( rewrite_patches_x[h]!=0 ) {
                read_position += patch_size[1];
                if( rewrite_xmin==0 ) {
                    write_position0 += write_sizey;
                }
                sx--;
            }

            write_position = write_position0;
            for( unsigned int ix=0; ix<sx; ix++ ) {
                if( rewrite_patches_y[h]!=0 ) {
                    read_position ++;
                }
                for( unsigned int iy=0; iy<sy; iy++ ) {
                    //data_rewrite[write_position] += h+1;
                    data_rewrite[write_position] = data_reread[read_position];
                    read_position ++;
                    write_position++;
                
                }
                write_position += write_skip_y;
            }
        }
    
        // Rewrite the file with the previously defined partition
        H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &( data_rewrite[0] ) );
    }
    
}

