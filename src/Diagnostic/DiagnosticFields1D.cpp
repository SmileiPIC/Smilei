
#include "DiagnosticFields1D.h"

#include "Params.h"
#include "Patch.h"
#include "Field1D.h"

using namespace std;

DiagnosticFields1D::DiagnosticFields1D( Params &params, SmileiMPI* smpi, Patch* patch, int ndiag )
    : DiagnosticFields( params, smpi, patch, ndiag )
{
}

DiagnosticFields1D::~DiagnosticFields1D()
{
}


// Copy patch field to current "data" buffer
void DiagnosticFields1D::getField( Patch* patch, int field_index )
{
    // Get current field
    Field1D* field;
    if( time_average>1 ) {
        field = static_cast<Field1D*>(patch->EMfields->allFields_avg[field_index]);
    } else {
        field = static_cast<Field1D*>(patch->EMfields->allFields    [field_index]);
    }
    // Copy field to the "data" buffer
    unsigned int ix = patch_offset[0];
    unsigned int ix_max = ix + patch_size[0];
    unsigned int iout = total_patch_size * (patch->Hindex()-refHindex);
    while( ix < ix_max ) {
        data[iout] = (*field)(ix);
        ix++;
        iout++;
    }
    
    if( time_average>1 ) field->put_to(0.0);
}
