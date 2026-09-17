#ifndef MF_SOLVER3D_TERZANI_H
#define MF_SOLVER3D_TERZANI_H

#include "Solver3D.h"
class ElectroMagn;

class MF_Solver3D_Terzani : public Solver3D
{

public:
    //! Creator for MF_Solver3D_Terzani
    MF_Solver3D_Terzani( Params &params );
    virtual ~MF_Solver3D_Terzani();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );

    // coefficient necessary to reduce the dispersion
    double delta;

protected:

};//END class

#endif

