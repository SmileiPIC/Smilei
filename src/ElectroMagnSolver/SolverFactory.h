#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include "MA_Solver1D_norm.h"
#include "MA_Solver1D_Friedman.h"
#include "MA_Solver2D_norm.h"
#include "MA_Solver2D_Friedman.h"
#include "MA_Solver3D_norm.h"
#include "MA_Solver3D_Friedman.h"
#include "MA_SolverAM_norm.h"
#include "MA_SolverAM_Friedman.h"
#include "MF_Solver1D_Yee.h"
#include "MF_Solver2D_Yee.h"
#include "MF_Solver3D_Yee.h"
#include "MF_SolverAM_Yee.h"
#include "MF_Solver2D_Grassi.h"
#include "MF_Solver2D_GrassiSpL.h"
#include "MF_Solver2D_Bouchard.h"
#include "MF_Solver3D_Bouchard.h"
#include "MF_Solver2D_Cowan.h"
#include "MF_Solver2D_Lehe.h"
#include "MF_Solver3D_Lehe.h"
#include "MF_SolverAM_Lehe.h"
#include "MF_SolverAM_Terzani.h"

#include "MF_Solver1D_M4.h"
#include "MF_Solver2D_M4.h"
#include "MF_Solver3D_M4.h"

#include "PXR_Solver2D_GPSTD.h"
#include "PXR_Solver3D_FDTD.h"
#include "PXR_Solver3D_GPSTD.h"
#include "PXR_SolverAM_GPSTD.h"

#include "PML_Solver2D_Bouchard.h"
#include "PML_Solver2D_Yee.h"
#include "PML_Solver3D_Bouchard.h"
#include "PML_Solver3D_Yee.h"
#include "PML_SolverAM.h"
#include "PML_Solver2D_Envelope.h"
#include "PML_Solver3D_Envelope.h"
#include "PML_SolverAM_Envelope.h"
#include "PML_SolverAM_EnvelopeReducedDispersion.h"

#include "Params.h"

#include "Tools.h"

class SolverFactory
{
public:

    // create Maxwell-Ampere solver
    // -----------------------------
    static Solver *createMA( Params &params )
    {
        Solver *solver = NULL;

        if( params.geometry == "1Dcartesian" ) {

            if( params.Friedman_filter ) {
                if (params.maxwell_sol != "Yee") ERROR( "Only Yee Maxwell solver is compatible with Friedman filter in 1Dcartesian geometry" );
                solver = new MA_Solver1D_Friedman( params );
            } else {
                solver = new MA_Solver1D_norm( params );
            }

        } else if( params.geometry == "2Dcartesian" ) {

            if( params.is_spectral ) {
                solver = new PXR_Solver2D_GPSTD( params );
            } else if( params.Friedman_filter ) {
                if ( (params.maxwell_sol != "Yee") && (params.maxwell_sol != "Bouchard") && (params.maxwell_sol != "Grassi") && (params.maxwell_sol != "GrassiSpL") ){
                    ERROR( "Only Yee, Bouchard, Grassi and GrassiSpL Maxwell solvers are compatible with Friedman filter in 2Dcartesian geometry" );
                }
                solver = new MA_Solver2D_Friedman( params );
            } else {
                solver = new MA_Solver2D_norm( params );
            }

        } else if( params.geometry == "3Dcartesian" ) {

            if( params.is_spectral ) {
                if( params.is_pxr ) {
                    solver = new PXR_Solver3D_GPSTD( params );
                } else {
                    ERROR( "Spectral solver not available without Picsar" );
                }
            } else {
                if( params.is_pxr ) {
                    solver = new PXR_Solver3D_FDTD( params );
                } else {
                    if( params.Friedman_filter ) {
                      if (params.maxwell_sol != "Yee") ERROR( "Only Yee Maxwell solver is compatible with Friedman filter in 3Dcartesian geometry" );
                      solver = new MA_Solver3D_Friedman( params );
                    } else {
                      solver = new MA_Solver3D_norm( params );
                    }

                }
            }

        } else if( params.geometry == "AMcylindrical" ) {

            if( params.is_pxr ) {
                solver = new PXR_SolverAM_GPSTD( params );
            } else {
                if( params.Friedman_filter ) {
                    if (params.maxwell_sol != "Yee") ERROR( "Only Yee Maxwell solver is compatible with Friedman filter in AMcylindrical geometry" );
                    solver = new MA_SolverAM_Friedman( params );
                } else {
                    solver = new MA_SolverAM_norm( params );
                }

            }

        }

        if( !solver ) {
            ERROR( "Unknown Maxwell-Ampere solver " );
        }

        return solver;
    };

    // Create Maxwell-Faraday solver
    // -----------------------------
    static Solver *createMF( Params &params )
    {
        Solver *solver = NULL;

        // Create the required solver for Faraday's Equation
        // -------------------------------------------------
        if( params.geometry == "1Dcartesian" ) {

            if( params.maxwell_sol == "Yee" ) {
                solver = new MF_Solver1D_Yee( params );
            } else if( params.maxwell_sol == "M4" ) {
                solver = new MF_Solver1D_M4( params );
            }

        } else if( params.geometry == "2Dcartesian" ) {

            if( params.maxwell_sol == "Yee" ) {
                solver = new MF_Solver2D_Yee( params );
            } else if( params.maxwell_sol == "Grassi" ) {
                solver = new MF_Solver2D_Grassi( params );
            } else if( params.maxwell_sol == "GrassiSpL" ) {
                solver = new MF_Solver2D_GrassiSpL( params );
            } else if( params.maxwell_sol == "Bouchard" ) {
                solver = new MF_Solver2D_Bouchard( params );
            } else if( params.maxwell_sol == "Cowan" ) {
                solver = new MF_Solver2D_Cowan( params );
            } else if( params.maxwell_sol == "Lehe" ) {
                solver = new MF_Solver2D_Lehe( params );
            } else if( params.maxwell_sol == "M4" ) {
                solver = new MF_Solver2D_M4( params );
            } else if( params.is_spectral ) {
                solver = new NullSolver();
            }

        } else if( params.geometry == "3Dcartesian" ) {

            if( params.maxwell_sol == "Yee" ) {
                solver = new MF_Solver3D_Yee( params );
            } else if( params.maxwell_sol == "Lehe" ) {
                solver = new MF_Solver3D_Lehe( params );
            } else if( params.maxwell_sol == "Bouchard" ) {
                solver = new MF_Solver3D_Bouchard( params );
            } else if( params.maxwell_sol == "M4" ) {
                solver = new MF_Solver3D_M4( params );
            } else if( params.is_pxr ) {
                solver = new NullSolver();
            }

        } else if( params.geometry == "AMcylindrical" ) {

            if( params.maxwell_sol == "Yee" ) {
                solver = new MF_SolverAM_Yee( params );
            } if( params.maxwell_sol == "Lehe" ) {
                solver = new MF_SolverAM_Lehe( params );
            } else if( params.maxwell_sol == "Terzani" ) {
                solver = new MF_SolverAM_Terzani( params );
            } else if( params.is_pxr ) {
                solver = new NullSolver();
            }

        }

        if( !solver ) {
            ERROR( "Unknwon solver '" << params.maxwell_sol << "' for geometry '" << params.geometry <<"'" );
        }

        return solver;
    };


    // Create PML solver
    // -----------------------------
    static Solver *createPML( Params &params )
    {
        Solver *solver = NULL;
        if( params.geometry == "2Dcartesian" ) {
            if (params.maxwell_sol=="Bouchard"){
                solver = new PML_Solver2D_Bouchard( params );
            }
            else {//(params.maxwell_sol=="Yee")
                solver = new PML_Solver2D_Yee( params );
            }
        }
        else if( params.geometry == "3Dcartesian" ) {
            if (params.maxwell_sol=="Bouchard"){
                solver = new PML_Solver3D_Bouchard( params );
            }
            else {
                solver = new PML_Solver3D_Yee( params );
            }
        }
        else if( params.geometry == "AMcylindrical" ) {
            solver = new PML_SolverAM( params );
        }
        else {
            ERROR( "PML configuration not implemented" );
        }

        return solver;
    }

    // create PML solver for envelope
    // -----------------------------
    static Solver *createPMLenvelope( Params &params )
    {
        Solver *solver = NULL;
        if( params.geometry == "2Dcartesian" ) {
            if (params.Laser_Envelope_model){
                if (params.envelope_solver == "explicit") {
                    solver = new PML_Solver2D_Envelope( params );
                }
                else if (params.envelope_solver == "explicit_reduced_dispersion") {
                    ERROR( "PML configuration not implemented yet" );
                }
                else {
                    ERROR( "PML configuration not implemented" );
                }
            }
        }
        else if( params.geometry == "3Dcartesian" ) {
            if (params.Laser_Envelope_model) {
                if (params.envelope_solver == "explicit") {
                    solver = new PML_Solver3D_Envelope( params );
                }
                else if (params.envelope_solver == "explicit_reduced_dispersion") {
                    ERROR( "PML configuration not implemented yet" );
                }
                else {
                    ERROR( "PML configuration not implemented" );
                }
            }
        }
        else if( params.geometry == "AMcylindrical" ) {
            if (params.Laser_Envelope_model){
                if (params.envelope_solver == "explicit") {
                    solver = new PML_SolverAM_Envelope( params );
                }
                else if (params.envelope_solver == "explicit_reduced_dispersion") {
                    solver = new PML_SolverAM_EnvelopeReducedDispersion( params );
                }
                else {
                    ERROR( "PML configuration not implemented" );
                }
            }
        }
        else {
            ERROR( "PML configuration not implemented" );
        }

        return solver;
    }

};

#endif
