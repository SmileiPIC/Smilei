#ifndef DiagnosticProbe2D_H
#define DiagnosticProbe2D_H

#include "DiagnosticProbe.h"

//! this class holds the point probe
class DiagnosticProbe2D : public DiagnosticProbe {
    
public:
    
    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe2D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
    ~DiagnosticProbe2D(){};
    
};
#endif
