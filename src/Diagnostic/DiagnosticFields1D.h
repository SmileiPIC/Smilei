
#ifndef DIAGNOSTICFIELDS1D_H
#define DIAGNOSTICFIELDS1D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields1D : public DiagnosticFields {
public:
    DiagnosticFields1D( Params &params, SmileiMPI* smpi, Patch* patch, int );
    ~DiagnosticFields1D();
    
    //! Copy patch field to current "data" buffer
    void getField( Patch* patch, int ) override;

    //! Basic write field on its own file (debug)
    void write( Field* field );


};

#endif
