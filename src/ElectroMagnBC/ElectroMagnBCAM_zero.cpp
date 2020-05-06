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
    Nmode = params.nmodes;
    number_of_damping_cells.push_back(params.number_of_damping_cells[0]);  
 
    if (params.uncoupled_grids)
        region_oversize_l = params.region_oversize[0];
    else
        region_oversize_l = 0;
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Zero Boundary Conditions
//
//
//   zero - damping area - unmodified   - region - unmodified - damping area - zero
//   <--------------------------------->           <------------------------------>
//          ghost cells                                       ghost cells                 
//          <-------------------------->           <------------------------> 
//                  damping cells                         damping cells       
//                        <------------>           <---------> 
//                        damping cells/2          damping cells/2       
//
//
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_zero::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    // Loop on imode
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        // Static cast of the fields
        cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
        
        if( min_max == 0 && patch->isXmin() ) {
            //x= Xmin
            int cell_start_dump =  region_oversize_l - number_of_damping_cells[0]/2;
            int cell_stop_dump  =  cell_start_dump - number_of_damping_cells[0]/2 ;
            double pi_ov_ndc = M_PI / number_of_damping_cells[0];
            double damping_phase = 0.;

            for (int i=cell_start_dump; i > cell_stop_dump; i--){
                //cos^2 damping over number_of_damping_cells/2 cells.
                damping_phase += pi_ov_ndc;
                double damp_coeff = cos(damping_phase);
                damp_coeff *= damp_coeff; 
            
                for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                    ( *El )( i, j ) *= damp_coeff;
                    ( *Er )( i, j ) *= damp_coeff;
                    ( *Et )( i, j ) *= damp_coeff;
                    ( *Bl )( i, j ) *= damp_coeff;
                    ( *Br )( i, j ) *= damp_coeff;
                    ( *Bt )( i, j ) *= damp_coeff;
                }
            }

            for (int i= 0; i<=cell_stop_dump; i++){
                for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                    ( *El )( i, j ) = 0.;
                    ( *Er )( i, j ) = 0.;
                    ( *Et )( i, j ) = 0.;
                    ( *Bl )( i, j ) = 0.;
                    ( *Br )( i, j ) = 0.;
                    ( *Bt )( i, j ) = 0.;
                }
            }

        } else if( min_max == 1 && patch->isXmax() ) {

            int cell_start_dump = nl_p - region_oversize_l + number_of_damping_cells[0]/2;
            int cell_stop_dump  =  cell_start_dump + number_of_damping_cells[0]/2 ;
            double pi_ov_ndc = M_PI / number_of_damping_cells[0];
            double damping_phase = 0.;
 
            for (int i=cell_start_dump; i < cell_stop_dump; i++){
                //cos^2 damping over number_of_damping_cells/2 cells.
                damping_phase += pi_ov_ndc;
                double damp_coeff  = cos(damping_phase);
                damp_coeff *= damp_coeff; 

                for( unsigned int j=0. ; j<nr_p ; j++ ) {
                    ( *El )( i, j ) *= damp_coeff;
                    ( *Er )( i, j ) *= damp_coeff;
                    ( *Et )( i, j ) *= damp_coeff;
                    ( *Bl )( i, j ) *= damp_coeff;
                    ( *Br )( i, j ) *= damp_coeff;
                    ( *Bt )( i, j ) *= damp_coeff;
                }
            }

            for (unsigned int i=cell_stop_dump; i < nl_p; i++){
                for( unsigned int j=0. ; j<nr_p ; j++ ) {
                    ( *El )( i, j ) = 0.;
                    ( *Er )( i, j ) = 0.;
                    ( *Et )( i, j ) = 0.;
                    ( *Bl )( i, j ) = 0.;
                    ( *Br )( i, j ) = 0.;
                    ( *Bt )( i, j ) = 0.;
                }
            }
        }
    }
}
