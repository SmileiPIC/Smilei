
#include "DiagnosticFields2D.h"

#include <sstream>

#include "Params.h"
#include "Patch.h"
#include "Field2D.h"

using namespace std;

DiagnosticFields2D::DiagnosticFields2D( Params &params, SmileiMPI* smpi, Patch* patch, int ndiag )
    : DiagnosticFields( params, smpi, patch, ndiag )
{
}

DiagnosticFields2D::~DiagnosticFields2D()
{
}


// Copy patch field to current "data" buffer
void DiagnosticFields2D::getField( Patch* patch, int field_index )
{
    // Get current field
    Field2D* field;
    if( time_average>1 ) {
        field = static_cast<Field2D*>(patch->EMfields->allFields_avg[field_index]);
    } else {
        field = static_cast<Field2D*>(patch->EMfields->allFields    [field_index]);
    }
    // Copy field to the "data" buffer
    unsigned int ix = patch_offset[0];
    unsigned int ix_max = ix + patch_size[0];
    unsigned int iy;
    unsigned int iy_max = patch_offset[1] + patch_size[1];
    unsigned int iout = total_patch_size * (patch->Hindex()-refHindex);
    while( ix < ix_max ) {
        iy = patch_offset[1];
        while( iy < iy_max ) {
            data[iout] = (*field)(ix, iy);
            iout++;
            iy++;
        }
        ix++;
    }
    
    if( time_average>1 ) field->put_to(0.0);
}

