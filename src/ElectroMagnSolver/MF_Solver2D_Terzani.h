#ifndef MF_SOLVER2D_TERZANI_H
#define MF_SOLVER2D_TERZANI_H

#include "Solver2D.h"
class ElectroMagn;

class MF_Solver2D_Terzani : public Solver2D
{

public:
    //! Creator for MF_Solver2D_Yee
    MF_Solver2D_Terzani( Params &params );
    virtual ~MF_Solver2D_Terzani();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    // coefficient necessary to reduce the dispersion
    double delta;
    
protected:
    // Check if time filter is applied or not
    bool isEFilterApplied;
    
};//END class

#endif

