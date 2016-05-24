
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
    
    //! Copy patch field to current "data" buffer
    void getField( Patch* patch, int ) override;

    //! Basic write field on its own file (debug)
    void write( Field* field );


};

#endif
