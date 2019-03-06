
#ifndef DIAGNOSTICCARTFIELDS3D_H
#define DIAGNOSTICCARTFIELDS3D_H

#include <string>
#include <vector>

#include "DiagnosticCartFields.h"
#include "Tools.h"

class DiagnosticCartFields3D : public DiagnosticCartFields
{
public:
    DiagnosticCartFields3D( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticCartFields3D();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    
    void writeField( hid_t, int ) override;
    
private:

    // Tools for re-reading and re-writing the file in a folded pattern
    hid_t filespace_reread, filespace_firstwrite, memspace_reread, memspace_firstwrite;
    std::vector<double> data_reread, data_rewrite;
    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_zmin, rewrite_npatchx, rewrite_npatchy, rewrite_npatchz;
    std::vector<unsigned int> rewrite_patches_x, rewrite_patches_y, rewrite_patches_z;
};

#endif
