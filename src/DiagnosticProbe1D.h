#ifndef DiagnosticProbe1D_H
#define DiagnosticProbe1D_H

#include "DiagnosticProbe.h"

//! this class holds the point probe
class DiagnosticProbe1D : public DiagnosticProbe {
    
public:
    
    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
    ~DiagnosticProbe1D(){};
    
};
#endif
