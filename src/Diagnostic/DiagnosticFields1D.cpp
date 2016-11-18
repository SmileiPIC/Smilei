
#include "DiagnosticFields1D.h"

#include "Params.h"
#include "Patch.h"
#include "Field1D.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields1D::DiagnosticFields1D( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, int ndiag )
    : DiagnosticFields( params, smpi, vecPatches, ndiag )
{
    // Calculate the offset in the local grid
    patch_offset_in_grid.resize(1);
    patch_offset_in_grid[0] = params.oversize[0];
    
    // Calculate the patch size
    patch_size.resize(1);
    patch_size[0] = params.n_space[0]+1;
    total_patch_size = patch_size[0];
    
    // define space in file
    hsize_t global_size[1];

    // All patch write n_space elements except, patch 0 which write n_space+1
    global_size[0] = tot_number_of_patches * (total_patch_size-1) + 1;

    filespace = H5Screate_simple(1, global_size, NULL);
    memspace  = H5Screate_simple(1, global_size, NULL ); // redefined later

}

DiagnosticFields1D::~DiagnosticFields1D()
{
}

void DiagnosticFields1D::setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches )
{
    // Calculate the total size of the array in this proc
    unsigned int total_vecPatches_size = (total_patch_size-1) * vecPatches.size();
    // One more cell on the left
    if (smpi->isMaster()) total_vecPatches_size++;
    
    // Resize the data (must contain data[0] even if not used for all process, so +1)
    data.resize(total_vecPatches_size+1);
    
    // Define offset and size for HDF5 file
    hsize_t offset[1], block[1], count[1];
    offset[0] = (total_patch_size-1) * refHindex;
    // non master patch/process write 1 point less, so start 1 point later
    if (!smpi->isMaster()) offset[0]++;

    block[0] = total_vecPatches_size;
    count[0] = 1;
    // Select portion of the file where this MPI will write to
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
    
    // define space in memory
    offset[0] = 0;
    // if not patch/process 0, data[0] not used (see in getField, iout++)
    if (!smpi->isMaster()) offset[0]++;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, block);
}


// Copy patch field to current "data" buffer
void DiagnosticFields1D::getField( Patch* patch, unsigned int ifield )
{
    // Get current field
    Field1D* field;
    if( time_average>1 ) {
        field = static_cast<Field1D*>(patch->EMfields->allFields_avg[diag_n][ifield]);
    } else {
        field = static_cast<Field1D*>(patch->EMfields->allFields[fields_indexes[ifield]]);
    }
    // Copy field to the "data" buffer
    unsigned int ix = patch_offset_in_grid[0];
    unsigned int ix_max = ix + patch_size[0];

    // if not patch 0, then offset = offset+1
    if (patch->hindex!=0) ix++;

    unsigned int iout = (total_patch_size-1) * (patch->Hindex()-refHindex);
    // patch 0 really write total_patch_size
    if (patch->hindex!=0) iout++; 

    while( ix < ix_max ) {
        data[iout] = (*field)(ix) * time_average_inv;
        ix++;
        iout++;
    }
    
    if( time_average>1 ) field->put_to(0.0);
}


// Write current buffer to file
void DiagnosticFields1D::writeField(hid_t dset_id, int timestep)
{
    
    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(data[0]) );
    
}

