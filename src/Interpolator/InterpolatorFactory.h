#ifndef INTERPOLATORFACTORY_H
#define INTERPOLATORFACTORY_H

#include "Interpolator.h"
#include "Interpolator1D2Order.h"
#include "Interpolator1D3Order.h"
#include "Interpolator1D4Order.h"
#include "Interpolator2D2Order.h"
#include "Interpolator2D4Order.h"
#include "Interpolator3D2Order.h"
#include "Interpolator3D4Order.h"
#include "InterpolatorAM1Order.h"
#include "InterpolatorAM2Order.h"

#ifdef _VECTO
#include "Interpolator2D2OrderV.h"
#include "Interpolator3D2OrderV.h"
#include "Interpolator3D4OrderV.h"
#endif

#include "Params.h"
#include "Patch.h"

#include "Tools.h"

class InterpolatorFactory
{
public:
    static Interpolator *create( Params &params, Patch *patch, bool vectorization )
    {
        Interpolator *Interp = NULL;
        // ---------------
        // 1Dcartesian simulation
        // ---------------
        if( ( params.geometry == "1Dcartesian" ) && ( params.interpolation_order == 2 ) ) {
            Interp = new Interpolator1D2Order( params, patch );
        } else if( ( params.geometry == "1Dcartesian" ) && ( params.interpolation_order == 4 ) ) {
            Interp = new Interpolator1D4Order( params, patch );
        }
        // ---------------
        // 2Dcartesian simulation
        // ---------------
        else if( ( params.geometry == "2Dcartesian" ) && ( params.interpolation_order == 2 ) ) {
            if( !vectorization ) {
                Interp = new Interpolator2D2Order( params, patch );
            }
#ifdef _VECTO
            else {
                Interp = new Interpolator2D2OrderV( params, patch );
            }
#endif
        } else if( ( params.geometry == "2Dcartesian" ) && ( params.interpolation_order == 4 ) ) {
            Interp = new Interpolator2D4Order( params, patch );
        }
        // ---------------
        // 3Dcartesian simulation
        // ---------------
        else if( ( params.geometry == "3Dcartesian" ) && ( params.interpolation_order == 2 ) ) {
            if( !vectorization ) {
                Interp = new Interpolator3D2Order( params, patch );
            }
#ifdef _VECTO
            else {
                Interp = new Interpolator3D2OrderV( params, patch );
            }
#endif
        } else if( ( params.geometry == "3Dcartesian" ) && ( params.interpolation_order == 4 ) ) {
            if( !vectorization ) {
                Interp = new Interpolator3D4Order( params, patch );
            }
#ifdef _VECTO
            else {
                Interp = new Interpolator3D4OrderV( params, patch );
            }
#endif
        }
        // ---------------
        // AM simulation
        // ---------------
        else if( params.geometry == "AMcylindrical" ) {
            if ( !params.is_spectral){
                Interp = new InterpolatorAM2Order( params, patch );
            } else {
                Interp = new InterpolatorAM1Order( params, patch );
            }
        } 
        else {
            ERROR( "Unknwon parameters : " << params.geometry << ", Order : " << params.interpolation_order );
        }
        
        return Interp;
    } // end InterpolatorFactory::create
    
    
};

#endif
