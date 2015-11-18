#ifndef DiagnosticPhaseMomMom_H
#define DiagnosticPhaseMomMom_H

#include "DiagnosticPhase.h"

//! class projection on the plane momentum-momentum
class DiagnosticPhaseMomMom : public DiagnosticPhase {
    
public:
    //! the component of the 1st momentum (0,1,2) and the component of the 2nd momentum (0,1,2) are chosen at construction
    DiagnosticPhaseMomMom(Params &params, unsigned int n_phase, const unsigned int directionMomentum1, const unsigned int directionMomentum2);	

	void run(partStruct&);
private:
    //! component of the momentum for the first axe
    const unsigned int my_dirMom1;
    //! component of the momentum for the second axe
    const unsigned int my_dirMom2;
};

#endif
