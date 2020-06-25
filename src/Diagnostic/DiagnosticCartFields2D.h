
#ifndef DIAGNOSTICCARTFIELDS2D_H
#define DIAGNOSTICCARTFIELDS2D_H

#include <string>
#include <vector>

#include "DiagnosticCartFields.h"
#include "Tools.h"

class DiagnosticCartFields2D : public DiagnosticCartFields
{
public:
    DiagnosticCartFields2D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticCartFields2D();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    
    H5Write writeField( H5Write*, std::string, int ) override;
    
private:

    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_npatchx, rewrite_npatchy;
    std::vector<unsigned int> rewrite_patches_x, rewrite_patches_y;
};

#endif
