#include "ElectroMagnBCAM_zero.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
//#include "Field2D.h"
#include "cField2D.h"
#include "Tools.h"
#include "Laser.h"
#include <complex>
#include "dcomplex.h"

using namespace std;

ElectroMagnBCAM_zero::ElectroMagnBCAM_zero( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBCAM( params, patch, _min_max )
{
    //Number of modes
    Nmode= params.nmodes;
    
    if (params.is_pxr)
        pxr_offset = params.oversize[0];
    else
        pxr_offset = 0;
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_zero::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if (pxr_offset) {
        if (nl_p==nl_d) {
            nl_p--;
            nr_p--;
        }
    }

    // Loop on imode
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        // Static cast of the fields
        cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
        
        if( min_max == 0 && patch->isXmin() ) {
            for( unsigned int j=0 ; j<nr_p ; j++ ) {
                //x= Xmin
                unsigned int i=0;
                ( *Br )( i, j ) = 0.;
            }
            for( unsigned int j=0.; j<nr_d ; j++ ) {
                //x=Xmin
                unsigned int i=0;
                ( *Bt )( i, j ) = 0.;
            }

        } else if( min_max == 1 && patch->isXmax() ) {
            for( unsigned int j=0. ; j<nr_p ; j++ ) {
                unsigned int i= nl_p;
                ( *Br )( i, j ) = 0.;
            }
            for( unsigned int j=0 ; j<nr_d ; j++ ) {
                unsigned int i= nl_p;
                ( *Bt )( i, j ) = 0.;
            }
        }
    }
}
