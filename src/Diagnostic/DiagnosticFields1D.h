
#ifndef DIAGNOSTICFIELDS1D_H
#define DIAGNOSTICFIELDS1D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields1D : public DiagnosticFields
{
public:
    DiagnosticFields1D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticFields1D();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    
    void writeField( hid_t, int ) override;
private:
    unsigned int MPI_start_in_file, total_patch_size;
};

#endif
