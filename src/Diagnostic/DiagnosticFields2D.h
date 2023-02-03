
#ifndef DIAGNOSTICFIELDS2D_H
#define DIAGNOSTICFIELDS2D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields2D : public DiagnosticFields
{
public:
    DiagnosticFields2D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticFields2D();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    
    H5Write writeField( H5Write*, std::string ) override;
    
private:

    std::vector<unsigned int> buffer_skip_x, buffer_skip_y;
};

#endif
