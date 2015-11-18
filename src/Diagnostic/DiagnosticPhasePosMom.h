#ifndef DiagnosticPhasePosMom_H
#define DiagnosticPhasePosMom_H

#include "DiagnosticPhase.h"

//! class projection on the plane position-momentum
class DiagnosticPhasePosMom : public DiagnosticPhase {
    
public:
    //! the component of the position (0,1,2) and the component of the momentum (0,1,2) are chosen at construction
    DiagnosticPhasePosMom(Params &params, unsigned int n_phase, const unsigned int directionPosition, const unsigned int directionMomentum);

	void run(partStruct&);
private:
    //! component of the position for the first axe
    const unsigned int my_dirPos;
    //! component of the momentum for the second axe
    const unsigned int my_dirMom;
};

#endif
