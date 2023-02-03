
#include "DiagnosticFields3D.h"

#include <sstream>
#include <cmath>
#include <limits>

#include "Params.h"
#include "Patch.h"
#include "Field3D.h"
#include "VectorPatch.h"
#include "DomainDecomposition.h"
#include "Hilbert_functions.h"
#include "LinearizedDomainDecomposition.h"

using namespace std;

DiagnosticFields3D::DiagnosticFields3D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &openPMD )
    : DiagnosticFields( params, smpi, vecPatches, ndiag, openPMD )
{
    
    // Calculate the offset in the local grid
    patch_offset_in_grid = { params.oversize[0]+1, params.oversize[1]+1, params.oversize[2]+1 };
   
    // Calculate the patch size
    patch_size_ = { params.patch_size_[0], params.patch_size_[1], params.patch_size_[2] };
    
    // Get the full size of the array in file
    vector<hsize_t> final_array_size(3);
    // Take subgrid into account
    for( unsigned int i=0; i<3; i++ ) {
        hsize_t start = 0;
        final_array_size[i] = params.number_of_patches[i] * params.patch_size_[i] + 1;
        findSubgridIntersection1( i, start, final_array_size[i], start );
    }
    // Define the chunk size (necessary above 2^28 points)
    const hsize_t max_size = 4294967295/2/sizeof( double );
    hsize_t final_size = final_array_size[0] * final_array_size[1] * final_array_size[2];
    vector<hsize_t> chunk_size;
    if( final_size > max_size ) {
        hsize_t n_chunks = 1 + ( final_size-1 ) / max_size;
        chunk_size.resize( 3 );
        chunk_size[0] = final_array_size[0] / n_chunks;
        chunk_size[1] = final_array_size[1];
        chunk_size[2] = final_array_size[2];
        if( n_chunks * chunk_size[0] < final_array_size[0] ) {
            chunk_size[0]++;
        }
    }
    
    filespace = new H5Space( final_array_size, {}, {}, chunk_size );
    memspace = new H5Space( 1 );
    
    // info for log output
    total_dataset_size = final_size;
}

DiagnosticFields3D::~DiagnosticFields3D()
{
}


struct PatchIXYZ {
    unsigned int i;
    unsigned int x;
    unsigned int y;
    unsigned int z;
};
bool patch_sorting( PatchIXYZ patch1_ixyz, PatchIXYZ patch2_ixyz ) {
    if( patch1_ixyz.x == patch2_ixyz.x ) {
        if( patch1_ixyz.y == patch2_ixyz.y ) {
            return patch1_ixyz.z < patch2_ixyz.z;
        } else {
            return patch1_ixyz.y < patch2_ixyz.y;
        }
    } else {
        return patch1_ixyz.x < patch2_ixyz.x;
    }
}

void DiagnosticFields3D::setFileSplitting( SmileiMPI *, VectorPatch &vecPatches )
{
    H5Sselect_none( filespace->sid_ );
    
    // Get all patch coordinates and sort them along Z then Y then X
    vector<PatchIXYZ> patch_ixyz( vecPatches.size() );
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        patch_ixyz[ipatch].i = ipatch;
        patch_ixyz[ipatch].x = vecPatches( ipatch )->Pcoordinates[0];
        patch_ixyz[ipatch].y = vecPatches( ipatch )->Pcoordinates[1];
        patch_ixyz[ipatch].z = vecPatches( ipatch )->Pcoordinates[2];
    }
    sort( patch_ixyz.begin(), patch_ixyz.end(), patch_sorting );
    
    buffer_skip_x.resize( vecPatches.size() );
    buffer_skip_y.resize( vecPatches.size() );
    buffer_skip_z.resize( vecPatches.size() );
    unsigned int current_z_skip = 0;
    vector<hsize_t> offset( 3 ), npoints( 3 ), count( 3, 1 ), start_in_patch( 3 );
    
    // Loop coordinates: X first, then Y, then Z
    // This loop does several things:
    //   - Add patches to the filespace one by one (this makes a combination of hyperslabs)
    //   - Calculate the way the data should be stored in the buffer so that it matches the filespace
    unsigned int i = 0;
    while( i < patch_ixyz.size() ) {
        
        // For this slab of patches at a given X, find the number of points in X
        offset[0] = patch_ixyz[i].x * patch_size_[0] + ( ( patch_ixyz[i].x==0 )?0:1 );
        npoints[0] = patch_size_[0] + ( ( patch_ixyz[i].x==0 )?1:0 );
        findSubgridIntersection1( 0, offset[0], npoints[0], start_in_patch[0] );
        
        // Now iterate on the patches along Y & Z that share the same X
        unsigned int i0 = i, npoints_yz = 0;
        while( i < patch_ixyz.size() && patch_ixyz[i].x == patch_ixyz[i0].x ) {
            
            // For this line of patches at a given Y, find the number of points in Y
            offset[1] = patch_ixyz[i].y * patch_size_[1] + ( ( patch_ixyz[i].y==0 )?0:1 );
            npoints[1] = patch_size_[1] + ( ( patch_ixyz[i].y==0 )?1:0 );
            findSubgridIntersection1( 1, offset[1], npoints[1], start_in_patch[1] );
            
            // Now iterate on the patches along Z that share the same Y & X
            unsigned int i1 = i, npoints_z = 0;
            while( i < patch_ixyz.size() && patch_ixyz[i].y == patch_ixyz[i1].y && patch_ixyz[i].x == patch_ixyz[i0].x ) {
                
                unsigned int ipatch = patch_ixyz[i].i;
                // Find the number of points along Z for this patch
                offset[2] = patch_ixyz[i].z * patch_size_[2] + ( ( patch_ixyz[i].z==0 )?0:1 );
                npoints[2] = patch_size_[2] + ( ( patch_ixyz[i].z==0 )?1:0 );
                findSubgridIntersection1( 2, offset[2], npoints[2], start_in_patch[2] );
                
                // Add this patch to the filespace
                H5Sselect_hyperslab( filespace->sid_, H5S_SELECT_OR, &offset[0], NULL, &count[0], &npoints[0] );
                
                // Calculate the initial skip when writing this patch to the buffer
                buffer_skip_z[ipatch] = current_z_skip + npoints_yz + npoints_z;
                
                // Also store, temporarily, the number of Z points for this patch alone
                buffer_skip_y[ipatch] = npoints[2];
                
                // Sum Z points to get the total Y skip when writing to the buffer
                npoints_z += npoints[2];
                
                i++;
            }
            
            // Now calculate the Y skip for each patch of that line
            for( unsigned int j=i1; j<i; j++ ) {
                unsigned int ipatch = patch_ixyz[j].i;
                buffer_skip_y[ipatch] = npoints_z - buffer_skip_y[ipatch];
                
                // Also store, temporarily, the number of Z&Y points for this line of patches alone
                buffer_skip_x[ipatch] = npoints[1] * npoints_z;
            }
            
             // Sum Z&Y points to get the total X skip when writing to the buffer
            npoints_yz += npoints[1] * npoints_z;
            
        }
        
        // Now calculate the X skip for each patch of that slab
        for( unsigned int j=i0; j<i; j++ ) {
            unsigned int ipatch = patch_ixyz[j].i;
            buffer_skip_x[ipatch] = npoints_yz - buffer_skip_x[ipatch];
        }
        current_z_skip += npoints[0] * npoints_yz;
    }
    delete memspace;
    memspace = new H5Space( current_z_skip );
    data.resize( current_z_skip );
}


// Copy patch field to current "data" buffer
void DiagnosticFields3D::getField( Patch *patch, unsigned int ifield )
{
    // Get current field
    Field3D *field;
    if( time_average>1 ) {
        field = static_cast<Field3D *>( patch->EMfields->allFields_avg[diag_n][ifield] );
    } else {
        field = static_cast<Field3D *>( patch->EMfields->allFields[fields_indexes[ifield]] );
    }
    
    // Find the intersection between this patch and the subgrid
    hsize_t patch_begin[3], patch_npoints[3], start_in_patch[3];
    for( unsigned int i=0; i<3; i++ ) {
        patch_begin  [i] = patch->Pcoordinates[i] * patch_size_[i] + ( ( patch->Pcoordinates[i]==0 )?0:1 );
        patch_npoints[i] = patch_size_[i] + ( ( patch->Pcoordinates[i]==0 )?1:0 );
        findSubgridIntersection1( i, patch_begin[i], patch_npoints[i], start_in_patch[i] );
        start_in_patch[i] += patch_offset_in_grid[i] - ( ( patch->Pcoordinates[i]==0 )?1:0 );
    }
    
    // Copy field to the "data" buffer
    unsigned int ix_max = start_in_patch[0] + subgrid_step_[0]*patch_npoints[0];
    unsigned int iy_max = start_in_patch[1] + subgrid_step_[1]*patch_npoints[1];
    unsigned int iz_max = start_in_patch[2] + subgrid_step_[2]*patch_npoints[2];
    unsigned int iout = buffer_skip_z[patch->Hindex()-refHindex];
    unsigned int stepy_out = buffer_skip_y[patch->Hindex()-refHindex];
    unsigned int stepx_out = buffer_skip_x[patch->Hindex()-refHindex];
    for( unsigned int ix = start_in_patch[0]; ix < ix_max; ix += subgrid_step_[0] ) {
        for( unsigned int iy = start_in_patch[1]; iy < iy_max; iy += subgrid_step_[1] ) {
            for( unsigned int iz = start_in_patch[2]; iz < iz_max; iz += subgrid_step_[2] ) {
                data[iout] = ( *field )( ix, iy, iz ) * time_average_inv;
                iout++;
            }
            iout += stepy_out;
        }
        iout += stepx_out;
    }
    
    if( time_average>1 ) {
        field->put_to( 0.0 );
    }
}


// Write current buffer to file
H5Write DiagnosticFields3D::writeField( H5Write * loc, string name )
{
    return loc->array( name, data[0], filespace, memspace, false, file_datatype_ );
}

