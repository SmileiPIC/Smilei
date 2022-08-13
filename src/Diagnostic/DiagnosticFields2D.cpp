
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
    patch_offset_in_grid = { params.oversize[0]+1, params.oversize[1]+1 };
    
    // Calculate the patch size
    patch_size = { params.n_space[0], params.n_space[1] };
    
    // Get the full size of the array in file
    vector<hsize_t> final_array_size(2);
    final_array_size[0] = params.number_of_patches[0] * params.n_space[0] + 1;
    final_array_size[1] = params.number_of_patches[1] * params.n_space[1] + 1;
    // Take subgrid into account
    for( unsigned int i=0; i<2; i++ ) {
        hsize_t start = 0;
        findSubgridIntersection1( i, start, final_array_size[i], start );
    }
    // Define the chunk size (necessary above 2^28 points)
    const hsize_t max_size = 4294967295/2/sizeof( double );
    hsize_t final_size = final_array_size[0] * final_array_size[1];
    if( final_size > max_size ) {
        hsize_t n_chunks = 1 + ( final_size-1 ) / max_size;
        chunk_size.resize( 2 );
        chunk_size[0] = final_array_size[0] / n_chunks;
        chunk_size[1] = final_array_size[1];
        if( n_chunks * chunk_size[0] < final_array_size[0] ) {
            chunk_size[0]++;
        }
    } else {
        chunk_size.resize( 0 );
    }
    
    filespace = new H5Space( final_array_size, {}, {}, chunk_size );
}

DiagnosticFields2D::~DiagnosticFields2D()
{
}


struct PatchIXY {
    unsigned int i;
    unsigned int x;
    unsigned int y;
};
bool patch_sorting( PatchIXY patch1_ixy, PatchIXY patch2_ixy ) {
    if( patch1_ixy.x == patch2_ixy.x ) {
        return patch1_ixy.y < patch2_ixy.y;
    } else {
        return patch1_ixy.x < patch2_ixy.x;
    }
}

void DiagnosticFields2D::setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches )
{
    H5Sselect_none( filespace->sid_ );
    
    // Get all patch coordinates and sort them along Y then X
    vector<PatchIXY> patch_ixy( vecPatches.size() );
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        patch_ixy[ipatch].i = ipatch;
        patch_ixy[ipatch].x = vecPatches( ipatch )->Pcoordinates[0];
        patch_ixy[ipatch].y = vecPatches( ipatch )->Pcoordinates[1];
    }
    sort( patch_ixy.begin(), patch_ixy.end(), patch_sorting );
    
    // Group patches by the same X
    vector<vector<PatchIXY> > patch_ixy_grouped;
    int current_x = -1;
    for( unsigned int j=0 ; j<patch_ixy.size() ; j++ ) {
        if( (int)( patch_ixy[j].x ) > current_x ) {
            vector<PatchIXY> a( 1, patch_ixy[j] );
            patch_ixy_grouped.push_back( a );
            current_x = patch_ixy[j].x;
        } else {
            patch_ixy_grouped.back().push_back( patch_ixy[j] );
        }
    }
    
    buffer_skip_x.resize( vecPatches.size() );
    buffer_skip_y.resize( vecPatches.size() );
    unsigned int current_y_skip = 0;
    vector<hsize_t> offset( 2 ), npoints( 2 ), count( 2, 1 ), start_in_patch( 2 );
    
    // Loop coordinates: X first, then Y
    // This loop does several things:
    //   - Add patches to the filespace one by one (this makes a combination of hyperslabs)
    //   - Calculate the way the data should be stored in the buffer so that it matches the filespace
    for( unsigned int ix=0; ix<patch_ixy_grouped.size(); ix++ ) {
        
        // For this line of patches at a given X, find the number of points in X
        Patch * p0 = vecPatches(  patch_ixy_grouped[ix][0].i );
        offset[0] = p0->Pcoordinates[0] * patch_size[0] + ( ( p0->Pcoordinates[0]==0 )?0:1 );
        npoints[0] = patch_size[0] + ( ( p0->Pcoordinates[0]==0 )?1:0 );
        findSubgridIntersection1( 0, offset[0], npoints[0], start_in_patch[0] );
        
        // Now iterate on the patches along Y
        unsigned int npoints_y = 0, npoints_x = npoints[0];
        for( unsigned int iy=0; iy<patch_ixy_grouped[ix].size(); iy++ ) {
            
            unsigned int ipatch = patch_ixy_grouped[ix][iy].i;
            // Find the number of points along Y for this patch
            Patch * p = vecPatches( ipatch );
            offset[1] = p->Pcoordinates[1] * patch_size[1] + ( ( p->Pcoordinates[1]==0 )?0:1 );
            npoints[1] = patch_size[1] + ( ( p->Pcoordinates[1]==0 )?1:0 );
            findSubgridIntersection1( 1, offset[1], npoints[1], start_in_patch[1] );
            
            // Add this patch to the filespace
            H5Sselect_hyperslab( filespace->sid_, H5S_SELECT_OR, &offset[0], NULL, &count[0], &npoints[0] );
            
            // Calculate the initial skip when writing this patch to the buffer
            buffer_skip_y[ipatch] = current_y_skip + npoints_y;
            
            // Also store, temporarily, the number of Y points for this patch alone
            buffer_skip_x[ipatch] = npoints[1];
            
            // Sum Y points to get the total X skip when writing to the buffer
            npoints_y += npoints[1];
        }
        
        // Now calculate the X skip for each patch
        for( unsigned int iy=0; iy<patch_ixy_grouped[ix].size(); iy++ ) {
            unsigned int ipatch = patch_ixy_grouped[ix][iy].i;
            buffer_skip_x[ipatch] = npoints_y - buffer_skip_x[ipatch];
        }
        current_y_skip += npoints_x * npoints_y;
    }
    memspace = new H5Space( current_y_skip );
    data.resize( current_y_skip );
    
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
    hsize_t patch_begin[2], patch_npoints[2], start_in_patch[2];
    for( unsigned int i=0; i<2; i++ ) {
        patch_begin  [i] = patch->Pcoordinates[i] * patch_size[i] + ( ( patch->Pcoordinates[i]==0 )?0:1 );
        patch_npoints[i] = patch_size[i] + ( ( patch->Pcoordinates[i]==0 )?1:0 );
        findSubgridIntersection1( i, patch_begin[i], patch_npoints[i], start_in_patch[i] );
        start_in_patch[i] += patch_offset_in_grid[i] - ( ( patch->Pcoordinates[i]==0 )?1:0 );
    }
    
    // Copy field to the "data" buffer
    unsigned int ix_max = start_in_patch[0] + subgrid_step_[0]*patch_npoints[0];
    unsigned int iy_max = start_in_patch[1] + subgrid_step_[1]*patch_npoints[1];
    unsigned int iout = buffer_skip_y[patch->Hindex()-refHindex];
    unsigned int step_out = buffer_skip_x[patch->Hindex()-refHindex];
    for( unsigned int ix = start_in_patch[0]; ix < ix_max; ix += subgrid_step_[0] ) {
        for( unsigned int iy = start_in_patch[1]; iy < iy_max; iy += subgrid_step_[1] ) {
            data[iout] = ( *field )( ix, iy ) * time_average_inv;
            iout++;
        }
        iout += step_out;
    }
    
    if( time_average>1 ) {
        field->put_to( 0.0 );
    }
}


// Write current buffer to file
H5Write DiagnosticFields2D::writeField( H5Write * loc, std::string name, int itime )
{
    return loc->array( name, data[0], filespace, memspace );
}

