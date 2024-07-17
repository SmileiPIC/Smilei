#ifndef PROJECTORFACTORY_H
#define PROJECTORFACTORY_H

#include "Projector.h"
#include "Projector1D2Order.h"
#include "Projector1D2OrderGPU.h"
#include "Projector1D4Order.h"
#include "Projector2D2Order.h"
#include "Projector2D2OrderGPU.h"
#include "Projector2D4Order.h"
#include "Projector3D2Order.h"
#include "Projector3D2OrderGPU.h"
#include "Projector3D4Order.h"
#include "ProjectorAM1Order.h"
#include "ProjectorAM1OrderRuyten.h"
#include "ProjectorAM2Order.h"

#include "Projector2D2OrderV.h"
#include "Projector2D4OrderV.h"
#include "Projector3D2OrderV.h"
#include "Projector3D4OrderV.h"
#include "ProjectorAM2OrderV.h"
#include "ProjectorAM1OrderRuytenV.h"

#include "Params.h"
#include "Patch.h"
#include "Tools.h"

class ProjectorFactory
{
public:
    static Projector *create( Params &params, Patch *patch, bool vectorization )
    {
        Projector *Proj = NULL;
        // ---------------
        // 1Dcartesian simulation
        // ---------------
        if( ( params.geometry == "1Dcartesian" ) && ( params.interpolation_order == ( unsigned int )2 ) ) {
            #if defined( SMILEI_ACCELERATOR_GPU )
                Proj = new Projector1D2OrderGPU( params, patch );
            #else
                Proj = new Projector1D2Order( params, patch );
            #endif
        } else if( ( params.geometry == "1Dcartesian" ) && ( params.interpolation_order == ( unsigned int )4 ) ) {
            Proj = new Projector1D4Order( params, patch );
        }
        // ---------------
        // 2Dcartesian simulation
        // ---------------
        else if( ( params.geometry == "2Dcartesian" ) && ( params.interpolation_order == ( unsigned int )2 ) ) {
            if( !vectorization ) {
                #if defined( SMILEI_ACCELERATOR_GPU )
                    Proj = new Projector2D2OrderGPU( params, patch );
                #else
                    Proj = new Projector2D2Order( params, patch );
                #endif
            }
            else {
                Proj = new Projector2D2OrderV( params, patch );
            }
        } else if( ( params.geometry == "2Dcartesian" ) && ( params.interpolation_order == ( unsigned int )4 ) ) {
            if( !vectorization ) {
                Proj = new Projector2D4Order( params, patch );
            }
            else {
                Proj = new Projector2D4OrderV( params, patch );
            }
        }
        // ---------------
        // 3Dcartesian simulation
        // ---------------
        else if( ( params.geometry == "3Dcartesian" ) && ( params.interpolation_order == ( unsigned int )2 ) ) {
            if( !vectorization ) {
                #if defined( SMILEI_ACCELERATOR_GPU )
                    Proj = new Projector3D2OrderGPU( params, patch );
                #else
                    Proj = new Projector3D2Order( params, patch );
                #endif
            }
            else {
                Proj = new Projector3D2OrderV( params, patch );
            }
        } else if( ( params.geometry == "3Dcartesian" ) && ( params.interpolation_order == ( unsigned int )4 ) ) {
            if( !vectorization ) {
                Proj = new Projector3D4Order( params, patch );
            }
            else {
                Proj = new Projector3D4OrderV( params, patch );
            }

        // ---------------
        // AM simulation
        // ---------------
       } else if( params.geometry == "AMcylindrical" ) {
            if (params.is_spectral){
                Proj = new ProjectorAM1Order( params, patch );
            } else {
                if( !vectorization ) {
                    if ( params.interpolation_order == 1 ) {
                        Proj = new ProjectorAM1OrderRuyten( params, patch );
                    } else {
                        Proj = new ProjectorAM2Order( params, patch );
                    }
                }
                else {
                    if ( params.interpolation_order == 1 ) {
                        Proj = new ProjectorAM1OrderRuytenV( params, patch );
                    } else {
                        Proj = new ProjectorAM2OrderV( params, patch );
                    }
                }
            }
        } else {
            ERROR_NAMELIST( "Unknwon parameters : " << params.geometry << ", Order : " << params.interpolation_order,
                LINK_NAMELIST + std::string("#interpolation_order")
            );
        }

        return Proj;
    }

};
#endif
