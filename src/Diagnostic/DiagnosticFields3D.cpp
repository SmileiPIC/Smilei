
#include "DiagnosticFields3D.h"

#include <sstream>
#include <cmath> 
#include <limits>

#include "Params.h"
#include "Patch.h"
#include "Field3D.h"
#include "VectorPatch.h"
#include "Hilbert_functions.h"

using namespace std;

DiagnosticFields3D::DiagnosticFields3D( Params &params, SmileiMPI* smpi, Patch* patch, int ndiag )
    : DiagnosticFields( params, smpi, patch, ndiag )
{
    
    // Calculate the offset in the local grid
    patch_offset_in_grid.resize(3);
    patch_offset_in_grid[0] = params.oversize[0];
    patch_offset_in_grid[1] = params.oversize[1];
    patch_offset_in_grid[2] = params.oversize[2];
    
    // Calculate the patch size
    patch_size.resize(3);
    patch_size[0] = params.n_space[0] + 1;
    patch_size[1] = params.n_space[1] + 1;
    patch_size[2] = params.n_space[2] + 1;
    total_patch_size = patch_size[0] * patch_size[1] * patch_size[2];
    
    
    // define space in file
    hsize_t global_size[1];
    global_size[0] = tot_number_of_patches * total_patch_size;
    filespace_firstwrite = H5Screate_simple(1, global_size, NULL);
    memspace_firstwrite = H5Screate_simple(1, global_size, NULL ); // redefined later
    
    // Define a second subset of the grid, which is unrelated to the current 
    // composition of vecPatches. It is used for a second writing of the file
    // in order to fold the Hilbert curve. This new subset is necessarily
    // rectangular for efficient writing.
    hsize_t offset[1], block[1], count[1];
    int nproc = smpi->getSize(), iproc = smpi->getRank();
    int npatch = params.tot_number_of_patches;
    int npatch_local = 1<<int(log2( ((double)npatch)/nproc ));
    int first_proc_with_less_patches = (npatch-npatch_local*nproc)/npatch_local;
    int first_patch_of_this_proc;
    if( iproc < first_proc_with_less_patches ) {
        npatch_local *= 2;
        first_patch_of_this_proc = npatch_local*iproc;
    } else {
        first_patch_of_this_proc = npatch_local*(first_proc_with_less_patches+iproc);
    }
    // Define space in file for re-reading
    filespace_reread = H5Screate_simple(1, global_size, NULL);
    offset[0] = total_patch_size * first_patch_of_this_proc;
    block [0] = total_patch_size * npatch_local;
    count [0] = 1;
    H5Sselect_hyperslab(filespace_reread, H5S_SELECT_SET, offset, NULL, count, block);
    // Define space in memory for re-reading
    memspace_reread = H5Screate_simple(1, block, NULL);
    data_reread.resize( block[0] );
    cout << "--------------------------------------------" << endl;
    cout << " data_reread size : " << block[0] << endl;
    // Define the list of patches for re-writing
    rewrite_npatch = (unsigned int)npatch_local;
    rewrite_patches_x.resize( rewrite_npatch );
    rewrite_patches_y.resize( rewrite_npatch );
    rewrite_patches_z.resize( rewrite_npatch );
    rewrite_xmin=numeric_limits<int>::max();
    rewrite_ymin=numeric_limits<int>::max();
    rewrite_zmin=numeric_limits<int>::max();
    unsigned int rewrite_xmax=0, rewrite_ymax=0, rewrite_zmax=0, x, y, z;
    for( unsigned int h=0; h<rewrite_npatch; h++) {
        generalhilbertindexinv(params.mi[0],  params.mi[1],  params.mi[2], &x, &y, &z, first_patch_of_this_proc+h);
        if(x<rewrite_xmin) rewrite_xmin=x;
        if(x>rewrite_xmax) rewrite_xmax=x;
        if(y<rewrite_ymin) rewrite_ymin=y;
        if(y>rewrite_ymax) rewrite_ymax=y;
        if(z<rewrite_zmin) rewrite_zmin=z;
        if(z>rewrite_zmax) rewrite_zmax=z;
        rewrite_patches_x[h] = x;
        rewrite_patches_y[h] = y;
        rewrite_patches_z[h] = z;
    }
    rewrite_npatchx = rewrite_xmax - rewrite_xmin + 1;
    rewrite_npatchy = rewrite_ymax - rewrite_ymin + 1;
    rewrite_npatchz = rewrite_zmax - rewrite_zmin + 1;
    // Define space in file for re-writing
    hsize_t final_array_size[3], offset2[3], block2[3], count2[3];
    final_array_size[0] = params.number_of_patches[0] * params.n_space[0] + 1;
    final_array_size[1] = params.number_of_patches[1] * params.n_space[1] + 1;
    final_array_size[2] = params.number_of_patches[2] * params.n_space[2] + 1;
    filespace = H5Screate_simple(3, final_array_size, NULL);
    offset2[0] = rewrite_xmin * params.n_space[0] + ((rewrite_xmin==0)?0:1);
    offset2[1] = rewrite_ymin * params.n_space[1] + ((rewrite_ymin==0)?0:1);
    offset2[2] = rewrite_zmin * params.n_space[2] + ((rewrite_zmin==0)?0:1);
    block2 [0] = rewrite_npatchx * params.n_space[0] + ((rewrite_xmin==0)?1:0);
    block2 [1] = rewrite_npatchy * params.n_space[1] + ((rewrite_ymin==0)?1:0);
    block2 [2] = rewrite_npatchz * params.n_space[2] + ((rewrite_zmin==0)?1:0);
    count2 [0] = 1;
    count2 [1] = 1;
    count2 [2] = 1;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset2, NULL, count2, block2);
    // Define space in memory for re-writing
    memspace = H5Screate_simple(3, block2, NULL);
    cout << " data_rewrite size : " << block2[0]*block2[1]*block2[2] << endl;
    cout << "--------------------------------------------" << endl;
    data_rewrite.resize( block2[0]*block2[1]*block2[2] );
    
    tmp_dset_id=0;
}

DiagnosticFields3D::~DiagnosticFields3D()
{
    H5Sclose( filespace_firstwrite );
    H5Sclose( memspace_firstwrite );
    H5Sclose( filespace_reread );
    H5Sclose( memspace_reread );
}


void DiagnosticFields3D::setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int total_vecPatches_size = total_patch_size * vecPatches.size();
    
    // Resize the data
    data.resize(total_vecPatches_size);
    
    // Define offset and size for HDF5 file
    hsize_t offset[1], block[1], count[1];
    offset[0] = total_patch_size * refHindex;
    block [0] = total_vecPatches_size;
    count [0] = 1;
    // Select portion of the file where this MPI will write to
    H5Sselect_hyperslab(filespace_firstwrite, H5S_SELECT_SET, offset, NULL, count, block);
    // define space in memory
    H5Sset_extent_simple(memspace_firstwrite, 1, block, block );
    
    // Create/Open temporary dataset
    htri_t status = H5Lexists(fileId_, "tmp", H5P_DEFAULT);
    if( status == 0 ) {
        hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
        tmp_dset_id  = H5Dcreate( fileId_, "tmp", H5T_NATIVE_DOUBLE, filespace_firstwrite, H5P_DEFAULT, pid, H5P_DEFAULT);
        H5Pclose(pid);
    } else {
        hid_t pid = H5Pcreate(H5P_DATASET_ACCESS);
        tmp_dset_id = H5Dopen( fileId_, "tmp", pid );
        H5Pclose(pid);
    }
}


// Copy patch field to current "data" buffer
void DiagnosticFields3D::getField( Patch* patch, unsigned int field_index )
{
    // Get current field
    Field3D* field;
    if( time_average>1 ) {
        field = static_cast<Field3D*>(patch->EMfields->allFields_avg[field_index]);
    } else {
        field = static_cast<Field3D*>(patch->EMfields->allFields    [field_index]);
    }
    // Copy field to the "data" buffer
    unsigned int ix = patch_offset_in_grid[0];
    unsigned int ix_max = ix + patch_size[0];
    unsigned int iy;
    unsigned int iy_max = patch_offset_in_grid[1] + patch_size[1];
    unsigned int iz;
    unsigned int iz_max = patch_offset_in_grid[2] + patch_size[2];
    unsigned int iout = total_patch_size * (patch->Hindex()-refHindex);

    while( ix < ix_max ) {
        iy = patch_offset_in_grid[1];
        while( iy < iy_max ) {
            iz = patch_offset_in_grid[2];
            while( iz < iz_max ) {
                data[iout] = (*field)(ix, iy, iz);
                iout++;
                iz++;
            }
            iy++;
        }
        ix++;

    }
    
    if( time_average>1 ) field->put_to(0.0);
}


// Write current buffer to file
void DiagnosticFields3D::writeField( hid_t dset_id, int timestep ) {
    
    // Write the buffer in a temporary location
    H5Dwrite( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_firstwrite, filespace_firstwrite, write_plist, &(data[0]) );
    
    // Read the file with the previously defined partition
    H5Dread( tmp_dset_id, H5T_NATIVE_DOUBLE, memspace_reread, filespace_reread, write_plist, &(data_reread[0]) );
    
    // Fold the data according to the Hilbert curve
    unsigned int read_position, write_position, read_skipZ, read_skipYZ, write_skip_y, write_skip_z, sx, sy, sz;

    unsigned int write_sizez  =  ( rewrite_npatchz*(patch_size[2]-1) + ((rewrite_zmin==0)?1:0) );
    unsigned int write_sizeyz = ( rewrite_npatchy*(patch_size[1]-1) + ((rewrite_ymin==0)?1:0) ) * write_sizez;


    MESSAGE( "Before pre-processing loop, rewrite_npatch = " << rewrite_npatch );
    read_position = 0;
    for( unsigned int h=0; h<rewrite_npatch; h++ ) {

        write_position =    (rewrite_patches_z[h]-rewrite_zmin)*(patch_size[2]-1) 
            + write_sizez *((rewrite_patches_y[h]-rewrite_ymin)*(patch_size[1]-1))
            + write_sizeyz*((rewrite_patches_x[h]-rewrite_xmin)*(patch_size[0]-1));
        cout << "write_position 0 =" << write_position  << endl;
        cout << "read_position 0 =" << read_position  << endl;

        write_skip_z = (rewrite_npatchz - 1)*(patch_size[2]-1);
        write_skip_y = (rewrite_npatchy - 1)*(patch_size[1]-1);


        read_skipZ  = 0;
        read_skipYZ = 0;
        sx = patch_size[0];
        sy = patch_size[1];
        sz = patch_size[2];

        if( rewrite_patches_z[h]!=0 ) {
            read_skipZ++;
            if( rewrite_zmin==0 ) {
                write_position++;
                write_skip_z++;
            } 
            sy--;
        }
        cout << "write_position 1 =" << write_position  << endl;
        cout << "read_position 1 =" << read_position  << endl;
        if( rewrite_patches_y[h]!=0 ) {
            read_skipYZ += 1;//patch_size[2];
            if( rewrite_ymin==0 ) {
                write_position += write_sizez;
                write_skip_y++;
                write_skip_y *= write_skip_z; // needs up to date write_skip_z
            } 
            sy--;
        }
        cout << "write_position 2 =" << write_position  << endl;
        cout << "read_position 2 =" << read_position  << endl;
        if( rewrite_patches_x[h]!=0 ) {
            read_position += patch_size[1]*patch_size[2];
            if( rewrite_xmin==0 ) write_position += write_sizeyz;
            sx--;
        }
        cout << "write_position 3 =" << write_position  << endl;
        cout << "read_position 3 =" << read_position  << endl;
        MESSAGE( "Before data_rewrite" );

        unsigned int localmax(0);
        unsigned int localmin(10000000);
        unsigned int rlocalmax(0);
        unsigned int rlocalmin(10000000);
        for( unsigned int ix=0; ix<sx; ix++ ) {
            for( unsigned int iy=0; iy<sy; iy++ ) {
                for( unsigned int iz=0; iz<sz; iz++ ) {
                    //data_rewrite[write_position] = data_reread[read_position];
                    if (write_position > localmax) localmax = write_position;
                    if (write_position < localmin) localmin = write_position;
                    if (read_position > rlocalmax) rlocalmax = read_position;
                    if (read_position < rlocalmin) rlocalmin = read_position;
                    read_position ++;
                    write_position++;

                }
                read_position  += read_skipZ;
                write_position += write_skip_z;
            }
            //read_position  += read_skipYZ;
            write_position += write_skip_y;
        }
        MESSAGE( "Re-processing done, max write_position = " << localmax << " - localmin : " << localmin );
        MESSAGE( "Re-processing done, max read_position = " << rlocalmax << " - localmin : " << rlocalmin );
    }

    // Rewrite the file with the previously defined partition
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(data_rewrite[0]) );
    
}

