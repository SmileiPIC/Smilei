
#ifndef DIAGNOSTICFIELDS2D_H
#define DIAGNOSTICFIELDS2D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields2D : public DiagnosticFields {
public:
    DiagnosticFields2D( Params &params, SmileiMPI* smpi, VectorPatch &vecPatches, int );
    ~DiagnosticFields2D();
    
    void setFileSplitting( SmileiMPI* smpi, VectorPatch& vecPatches ) ;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch* patch, unsigned int ) ;
    
    void writeField(hid_t, int) ;

private:
    
    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_npatchx, rewrite_npatchy;
    std::vector<unsigned int> rewrite_patches_x, rewrite_patches_y;
};

#endif
