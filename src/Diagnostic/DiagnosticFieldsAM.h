
#ifndef DIAGNOSTICFIELDSAM_H
#define DIAGNOSTICFIELDSAM_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFieldsAM : public DiagnosticFields
{
public:
    DiagnosticFieldsAM( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticFieldsAM();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    
    void writeField( hid_t, int ) override;
    
private:

    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_npatchx, rewrite_npatchy;
    std::vector<unsigned int> rewrite_patches_x, rewrite_patches_y;
    
    std::vector<std::complex<double>> idata_reread, idata_rewrite, idata;
    
};

#endif
