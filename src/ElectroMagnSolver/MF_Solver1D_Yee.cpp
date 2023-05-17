
#include "MF_Solver1D_Yee.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MF_Solver1D_Yee::MF_Solver1D_Yee( Params &params )
    : Solver1D( params )
{
    isEFilterApplied = params.Friedman_filter;
}

MF_Solver1D_Yee::~MF_Solver1D_Yee()
{
}

void MF_Solver1D_Yee::operator()( ElectroMagn *fields )
{
    // const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    
    // Static-cast of the fields
    Field1D* Ey1D;
    Field1D* Ez1D;
    if (isEFilterApplied) {
        Ey1D = static_cast<Field1D*>(fields->filter_->Ey_[0]);
        Ez1D = static_cast<Field1D*>(fields->filter_->Ez_[0]);
    } else {
        Ey1D = static_cast<Field1D*>(fields->Ey_);
        Ez1D = static_cast<Field1D*>(fields->Ez_);
    }
    Field1D *By1D   = static_cast<Field1D *>( fields->By_ );
    Field1D *Bz1D   = static_cast<Field1D *>( fields->Bz_ );
    // ---------------------
    // Solve Maxwell-Faraday
    // ---------------------
    // NB: bx is given in 1d and defined when initializing the fields (here put to 0)
    // Transverse fields  by & bz are defined on the dual grid
    //for (unsigned int ix=1 ; ix<nx_p ; ix++) {
    for( unsigned int ix=1 ; ix<nx_d-1 ; ix++ ) {
        ( *By1D )( ix )= ( *By1D )( ix ) + dt_ov_dx * ( ( *Ez1D )( ix ) - ( *Ez1D )( ix-1 ) ) ;
        ( *Bz1D )( ix )= ( *Bz1D )( ix ) - dt_ov_dx * ( ( *Ey1D )( ix ) - ( *Ey1D )( ix-1 ) ) ;
    }
}
