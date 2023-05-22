#ifndef MF_SOLVERAM_LEHE_H
#define MF_SOLVERAM_LEHE_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class MF_SolverAM_Lehe : public SolverAM
{

public:
    //! Creator for MF_SolverAM_Yee
    MF_SolverAM_Lehe( Params &params );
    virtual ~MF_SolverAM_Lehe();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );
    
    double delta_l;
    double beta_rl;
    double beta_tl;
    double alpha_l;
    double alpha_t;
    double alpha_r;
    
protected:
    
};//END class

#endif

