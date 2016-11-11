
#include "DiagnosticFields1D.h"

#include "Params.h"
#include "Patch.h"
#include "Field1D.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields1D::DiagnosticFields1D( Params &params, SmileiMPI* smpi, Patch* patch, int ndiag )
    : DiagnosticFields( params, smpi, patch, ndiag )
{
    // Calculate the offset in the local grid
    patch_offset_in_grid.resize(1);
    patch_offset_in_grid[0] = params.oversize[0];
    
    // Calculate the patch size
    patch_size.resize(1);
    patch_size[0] = params.n_space[0];
    total_patch_size = patch_size[0];
    
    // define space in file
    hsize_t global_size[1];
    global_size[0] = tot_number_of_patches * total_patch_size + 1;
    filespace = H5Screate_simple(1, global_size, NULL);
    memspace  = H5Screate_simple(1, global_size, NULL ); // redefined later

}

DiagnosticFields1D::~DiagnosticFields1D()
{
}

void DiagnosticFields1D::setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int total_vecPatches_size = total_patch_size * vecPatches.size();
    if( vecPatches(0)->hindex==0 ) total_vecPatches_size++; // One more cell on the left
    
    // Resize the data
    data.resize(total_vecPatches_size);
    
    // Define offset and size for HDF5 file
    hsize_t offset[1], block[1], count[1];
    offset[0] = total_patch_size * refHindex;
    if( vecPatches(0)->hindex!=0 ) offset[0]++;
    block[0] = total_vecPatches_size;
    count[0] = 1;
    // Select portion of the file where this MPI will write to
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
    // define space in memory
    H5Sset_extent_simple(memspace, 1, block, block );
}


// Copy patch field to current "data" buffer
void DiagnosticFields1D::getField( Patch* patch, unsigned int field_index )
{
    // Get current field
    Field1D* field;
    if( time_average>1 ) {
        field = static_cast<Field1D*>(patch->EMfields->allFields_avg[field_index]);
    } else {
        field = static_cast<Field1D*>(patch->EMfields->allFields    [field_index]);
    }
    // Copy field to the "data" buffer
    unsigned int ix = patch_offset_in_grid[0] ;
    unsigned int iout = patch_size[0] * (patch->Hindex()-refHindex);
    memcpy(&data[iout], &(*field)(ix), patch_size[0]*sizeof(double)); 

    if( time_average>1 ) field->put_to(0.0);
}


// Write current buffer to file
void DiagnosticFields1D::writeField(hid_t dset_id, int timestep)
{
    
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(data[0]) );
    
}

