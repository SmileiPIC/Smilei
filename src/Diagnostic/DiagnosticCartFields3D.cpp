
#include "DiagnosticCartFields3D.h"

#include <sstream>
#include <cmath>
#include <limits>

#include "Params.h"
#include "Patch.h"
#include "Field3D.h"
#include "VectorPatch.h"
#include "DomainDecomposition.h"
#include "Hilbert_functions.h"

using namespace std;

DiagnosticCartFields3D::DiagnosticCartFields3D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &openPMD )
    : DiagnosticCartFields( params, smpi, vecPatches, ndiag, openPMD )
{

    // Calculate the offset in the local grid
    patch_offset_in_grid.resize( 3 );
    patch_offset_in_grid[0] = params.oversize[0];
    patch_offset_in_grid[1] = params.oversize[1];
    patch_offset_in_grid[2] = params.oversize[2];
    
    // Calculate the patch size
    patch_size.resize( 3 );
    patch_size[0] = params.n_space[0]*params.global_factor[0] + 1;
    patch_size[1] = params.n_space[1]*params.global_factor[1] + 1;
    patch_size[2] = params.n_space[2]*params.global_factor[2] + 1;
    total_patch_size = patch_size[0] * patch_size[1] * patch_size[2];
    
    // Define a second portion of the grid, which is unrelated to the current
    // composition of vecPatches. It is used for a second writing of the file
    // in order to fold the Hilbert curve. This new portion is necessarily
    // rectangular for efficient writing.
    int nproc = smpi->getSize(), iproc = smpi->getRank();
    int npatch = params.tot_number_of_patches;
    int npatch_local = 1<<int( log2( ( ( double )npatch )/nproc ) );
    int first_proc_with_less_patches = ( npatch-npatch_local*nproc )/npatch_local;
    if( iproc < first_proc_with_less_patches ) {
        npatch_local *= 2;
    }
    // Define the list of patches for re-writing
    rewrite_npatch = ( unsigned int )npatch_local;
    rewrite_patches_x.resize( rewrite_npatch );
    rewrite_patches_y.resize( rewrite_npatch );
    rewrite_patches_z.resize( rewrite_npatch );
    rewrite_xmin=numeric_limits<int>::max();
    rewrite_ymin=numeric_limits<int>::max();
    rewrite_zmin=numeric_limits<int>::max();
    unsigned int rewrite_xmax=0, rewrite_ymax=0, rewrite_zmax=0;
    for( unsigned int h=0; h<rewrite_npatch; h++ ) {
        std::vector<unsigned int> xcall( 3, 0 );
        
        xcall = vecPatches( 0 )->Pcoordinates;
        
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
        if( xcall[2]<rewrite_zmin ) {
            rewrite_zmin=xcall[2];
        }
        if( xcall[2]>rewrite_zmax ) {
            rewrite_zmax=xcall[2];
        }
        rewrite_patches_x[h] = xcall[0];
        rewrite_patches_y[h] = xcall[1];
        rewrite_patches_z[h] = xcall[2];
    }
    rewrite_npatchx = rewrite_xmax - rewrite_xmin + 1;
    rewrite_npatchy = rewrite_ymax - rewrite_ymin + 1;
    rewrite_npatchz = rewrite_zmax - rewrite_zmin + 1;
    // Define space in file for re-writing
    vector<hsize_t> final_array_size(3), offset2(3), block2(3);
    final_array_size[0] = params.number_of_patches[0] * params.n_space[0] + 1;
    final_array_size[1] = params.number_of_patches[1] * params.n_space[1] + 1;
    final_array_size[2] = params.number_of_patches[2] * params.n_space[2] + 1;
    offset2[0] = rewrite_xmin * params.n_space[0]*params.global_factor[0] + ( ( rewrite_xmin==0 )?0:1 );
    offset2[1] = rewrite_ymin * params.n_space[1]*params.global_factor[1] + ( ( rewrite_ymin==0 )?0:1 );
    offset2[2] = rewrite_zmin * params.n_space[2]*params.global_factor[2] + ( ( rewrite_zmin==0 )?0:1 );
    block2 [0] = rewrite_npatchx * params.n_space[0]*params.global_factor[0] + ( ( rewrite_xmin==0 )?1:0 );
    block2 [1] = rewrite_npatchy * params.n_space[1]*params.global_factor[1] + ( ( rewrite_ymin==0 )?1:0 );
    block2 [2] = rewrite_npatchz * params.n_space[2]*params.global_factor[2] + ( ( rewrite_zmin==0 )?1:0 );
    filespace = new H5Space( final_array_size, offset2, block2 );
    // Define space in memory for re-writing
    memspace = new H5Space( block2 );
    data_rewrite.resize( block2[0]*block2[1]*block2[2], 0. );
    
}

DiagnosticCartFields3D::~DiagnosticCartFields3D()
{
}


void DiagnosticCartFields3D::setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int total_vecPatches_size = total_patch_size * vecPatches.size();
    
    // Resize the data
    data.resize( total_vecPatches_size, 0. );
    
}


// Copy patch field to current "data" buffer
void DiagnosticCartFields3D::getField( Patch *patch, unsigned int ifield )
{
    // Get current field
    Field3D *field;
    if( time_average>1 ) {
        field = static_cast<Field3D *>( patch->EMfields->allFields_avg[diag_n][ifield] );
    } else {
        field = static_cast<Field3D *>( patch->EMfields->allFields[fields_indexes[ifield]] );
    }
    
    // Copy field to the "data" buffer
    unsigned int ix_max = patch_offset_in_grid[0] + patch_size[0];
    unsigned int iy_max = patch_offset_in_grid[1] + patch_size[1];
    unsigned int iz_max = patch_offset_in_grid[2] + patch_size[2];
    unsigned int iout = total_patch_size * ( patch->Hindex()-refHindex );
    for( unsigned int ix = patch_offset_in_grid[0]; ix < ix_max; ix++ ) {
        for( unsigned int iy = patch_offset_in_grid[1]; iy < iy_max; iy++ ) {
            for( unsigned int iz = patch_offset_in_grid[2]; iz < iz_max; iz++ ) {
                data[iout] = ( *field )( ix, iy, iz ) * time_average_inv;
                iout++;
            }
        }
    }
//    for (unsigned int ix = patch_offset_in_grid[0]; ix < ix_max; ix++){
//        for (unsigned int iy = patch_offset_in_grid[1]; iy < iy_max; iy++){
//            memcpy( data_pt, &((*field)(ix, iy, iz)), patch_size[2]*sizeof(double));
//            data_pt += patch_size[2];
//        }
//    }

    if( time_average>1 ) {
        field->put_to( 0.0 );
    }
}


// Write current buffer to file
H5Write DiagnosticCartFields3D::writeField( H5Write * loc, string name, int itime )
{

    
    // Fold the data according to the Hilbert curve
    unsigned int read_position, write_position, write_skip_y, write_skip_z, sx, sy, sz;
    unsigned int write_sizez  = ( rewrite_npatchz*( patch_size[2]-1 ) + ( ( rewrite_zmin==0 )?1:0 ) );
    unsigned int write_sizeyz = ( rewrite_npatchy*( patch_size[1]-1 ) + ( ( rewrite_ymin==0 )?1:0 ) ) * write_sizez;
    
    read_position = 0;
    for( unsigned int h=0; h<rewrite_npatch; h++ ) {
        int write_position0 = ( rewrite_patches_z[h]-rewrite_zmin )*( patch_size[2]-1 )
                              + write_sizez *( ( rewrite_patches_y[h]-rewrite_ymin )*( patch_size[1]-1 ) )
                              + write_sizeyz*( ( rewrite_patches_x[h]-rewrite_xmin )*( patch_size[0]-1 ) );
                              
        write_skip_z = ( rewrite_npatchz - 1 )*( patch_size[2]-1 );
        write_skip_y = write_sizeyz;
        
        sx = patch_size[0];
        sy = patch_size[1];
        sz = patch_size[2];
        
        if( rewrite_patches_z[h]!=0 ) {
            if( rewrite_zmin==0 ) {
                write_position0++;
                write_skip_z++;
            }
            sz--;
        }
        if( rewrite_patches_y[h]!=0 ) {
            if( rewrite_ymin==0 ) {
                write_position0 += write_sizez;
            }
            sy--;
        }
        if( rewrite_patches_x[h]!=0 ) {
            read_position += patch_size[1]*patch_size[2];
            if( rewrite_xmin==0 ) {
                write_position0 += write_sizeyz;
            }
            sx--;
        }
        
        write_position = write_position0;
        for( unsigned int ix=0; ix<sx; ix++ ) {
            if( rewrite_patches_y[h]!=0 ) {
                read_position+=patch_size[2];
            }
            for( unsigned int iy=0; iy<sy; iy++ ) {
                if( rewrite_patches_z[h]!=0 ) {
                    read_position ++;
                }
                for( unsigned int iz=0; iz<sz; iz++ ) {
                    //data_rewrite[write_position] += h+1;
                    data_rewrite[write_position] = data[read_position];
                    read_position ++;
                    write_position++;
                }
                write_position += write_skip_z;
            }
            write_position = write_position0+( ix+1 )*write_skip_y;
        }
        
    }
    
    // Rewrite the file with the previously defined partition
    return loc->array( name, data_rewrite[0], H5T_NATIVE_DOUBLE, filespace, memspace );
}

