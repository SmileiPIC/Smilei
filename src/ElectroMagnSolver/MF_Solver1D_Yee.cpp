
#include "MF_Solver1D_Yee.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MF_Solver1D_Yee::MF_Solver1D_Yee(PicParams &params)
    : Solver1D(params)
{
}

MF_Solver1D_Yee::~MF_Solver1D_Yee()
{
}

void MF_Solver1D_Yee::operator() ( ElectroMagn* fields )
{
    Field1D* Ey1D   = static_cast<Field1D*>(fields->Ey_);
    Field1D* Ez1D   = static_cast<Field1D*>(fields->Ez_);
    Field1D* By1D   = static_cast<Field1D*>(fields->By_);
    Field1D* Bz1D   = static_cast<Field1D*>(fields->Bz_);
    
    // ---------------------
    // Solve Maxwell-Faraday
    // ---------------------
    // NB: bx is given in 1d and defined when initializing the fields (here put to 0)
    // Transverse fields  by & bz are defined on the dual grid
    //for (unsigned int ix=1 ; ix<nx_p ; ix++) {
    for (int ix=1 ; ix<nx_d-1 ; ix++) {
        (*By1D)(ix)= (*By1D)(ix) + dt_ov_dx * ( (*Ez1D)(ix) - (*Ez1D)(ix-1)) ;
        (*Bz1D)(ix)= (*Bz1D)(ix) - dt_ov_dx * ( (*Ey1D)(ix) - (*Ey1D)(ix-1)) ;
    }
}

