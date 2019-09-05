
#ifndef DIAGNOSTICFIELDS3D_H
#define DIAGNOSTICFIELDS3D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields3D : public DiagnosticFields
{
public:
    DiagnosticFields3D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticFields3D();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    
    void writeField( hid_t, int ) override;
    
private:

    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_zmin, rewrite_npatchx, rewrite_npatchy, rewrite_npatchz;
    unsigned int rewrite_size[3], rewrite_start_in_file[3];
    std::vector<std::vector<unsigned int> > rewrite_patch;
};

#endif
