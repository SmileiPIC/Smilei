#ifndef DiagnosticProbe0D_H
#define DiagnosticProbe0D_H

#include "DiagnosticProbe.h"

//! this class holds the point probe
class DiagnosticProbe0D : public DiagnosticProbe {
    
public:
    
    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe0D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
    ~DiagnosticProbe0D(){};

};
#endif
