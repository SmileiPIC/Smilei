#ifndef PARTCOMPTIMEFACTORY_H
#define PARTCOMPTIMEFACTORY_H

#include "PartCompTime.h"
#include "PartCompTime3D2Order.h"
#include "PartCompTime3D4Order.h"

#include "Params.h"
#include "Tools.h"

class PartCompTimeFactory
{
public:
    static  PartCompTime*create( Params &params )
    {
        
        PartCompTime * part_comp_time = NULL;
        // ---------------
        // 1Dcartesian simulation
        // ---------------
        if( ( params.geometry == "1Dcartesian" ) && ( params.interpolation_order == 2 ) ) {
            
            // part_comp_time = new PartCompTime1D2Order( params, patch );
            
        } else if( ( params.geometry == "1Dcartesian" ) && ( params.interpolation_order == 4 ) ) {
            
            // part_comp_time = new PartCompTime1D4Order( params, patch );
            
        }
        // ---------------
        // 2Dcartesian simulation
        // ---------------
        else if( ( params.geometry == "2Dcartesian" ) && ( params.interpolation_order == 2 ) ) {
            
            // part_comp_time = new PartCompTime2D2Order( params, patch );

        } else if( ( params.geometry == "2Dcartesian" ) && ( params.interpolation_order == 4 ) ) {

            // part_comp_time = new PartCompTime2D4Order( params, patch );

        }
        // ---------------
        // 3Dcartesian simulation
        // ---------------
        else if( ( params.geometry == "3Dcartesian" ) && ( params.interpolation_order == 2 ) ) {
            
            part_comp_time = new PartCompTime3D2Order();
            
        } else if( ( params.geometry == "3Dcartesian" ) && ( params.interpolation_order == 4 ) ) {
            
            part_comp_time = new PartCompTime3D4Order();
            
        }
        // ---------------
        // AM simulation
        // ---------------
        else if( params.geometry == "AMcylindrical" ) {
            if ( !params.is_spectral){

                // part_comp_time = new PartCompTimeAM2Order( params );

            } else {
                // part_comp_time = new PartCompTimeAM1Order( params );
            }
        }
        else {
            ERROR_NAMELIST( "Unknwon parameters setting the particle computing time evaluation operator: "
             << params.geometry << ", Order : " << params.interpolation_order,
            LINK_NAMELIST + std::string("#main-variables"));
        }
        
        return part_comp_time;
    } // end PartCompTime::create
};

#endif
