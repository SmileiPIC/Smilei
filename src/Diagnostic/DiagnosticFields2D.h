
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
    
    H5Write writeField( H5Write*, std::string, int ) override;
    
private:

    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_npatchx, rewrite_npatchy;
    unsigned int rewrite_size[2], rewrite_start_in_file[2];
    std::vector<std::vector<unsigned int> > rewrite_patch;
};

#endif
