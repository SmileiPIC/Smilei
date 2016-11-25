
#ifndef DIAGNOSTICFIELDS3D_H
#define DIAGNOSTICFIELDS3D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields3D : public DiagnosticFields {
public:
    DiagnosticFields3D( Params &params, SmileiMPI* smpi, VectorPatch &vecPatches, int );
    ~DiagnosticFields3D();
    
    void setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches ) ;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch* patch, unsigned int ) ;
    
    void writeField(hid_t, int) ;

private:
    
    // Tools for re-reading and re-writing the file in a folded pattern
    hid_t filespace_reread, filespace_firstwrite, memspace_reread, memspace_firstwrite;
    std::vector<double> data_reread, data_rewrite;
    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_zmin, rewrite_npatchx, rewrite_npatchy, rewrite_npatchz;
    std::vector<unsigned int> rewrite_patches_x, rewrite_patches_y, rewrite_patches_z;
};

#endif
