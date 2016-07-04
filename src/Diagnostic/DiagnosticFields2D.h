
#ifndef DIAGNOSTICFIELDS2D_H
#define DIAGNOSTICFIELDS2D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields2D : public DiagnosticFields {
public:
    DiagnosticFields2D( Params &params, SmileiMPI* smpi, Patch* patch, int );
    ~DiagnosticFields2D();
    
    void setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch* patch, unsigned int ) override;
    
    void writeField(hid_t, int) override;

private:
    
    // Tools for re-reading and re-writing the file in a folded pattern
    hid_t filespace_reread, filespace_firstwrite, memspace_reread, memspace_firstwrite;
    std::vector<double> data_reread, data_rewrite;
    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_npatchx, rewrite_npatchy;
    std::vector<unsigned int> rewrite_patches_x, rewrite_patches_y;
};

#endif
