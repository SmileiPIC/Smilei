#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include "MA_Solver1D_norm.h"
#include "MA_Solver2D_norm.h"
#include "MA_Solver2D_Friedman.h"
#include "MA_Solver3D_norm.h"
#include "MA_SolverAM_norm.h"
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

#include "PXR_Solver2D_GPSTD.h"
#include "PXR_Solver3D_FDTD.h"
#include "PXR_Solver3D_GPSTD.h"
#include "PXR_SolverAM_GPSTD.h"

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
            solver = new MA_Solver1D_norm( params );
        } else if( params.geometry == "2Dcartesian" ) {
            if( params.is_pxr == false ) {
                if( params.is_spectral ) {
                    WARNING( "PS solveur are not available without Picsar" );
                }
                if( params.Friedman_filter ) {
                    solver = new MA_Solver2D_Friedman( params );
                } else {
                    solver = new MA_Solver2D_norm( params );
                }
            } else if( params.is_spectral == true ) {
                solver = new PXR_Solver2D_GPSTD( params );
            }
            
        } else if( params.geometry == "3Dcartesian" ) {
            if( params.is_pxr == false ) {
                if( params.is_spectral ) {
                    WARNING( "PS solveur are not available without Picsar" );
                }
                solver = new MA_Solver3D_norm( params );
            } else if( ( params.is_pxr == true ) && ( params.is_spectral == false ) ) {
                solver = new PXR_Solver3D_FDTD( params );
            } else if( ( params.is_pxr == true ) && ( params.is_spectral == true ) ) {
                solver = new PXR_Solver3D_GPSTD( params );
            }
        } else if( params.geometry == "AMcylindrical" ) {
            if( params.is_pxr == false )
                solver = new MA_SolverAM_norm( params );
            else
                solver = new PXR_SolverAM_GPSTD( params );
        }
        
        if( !solver ) {
            ERROR( "Unknwon Maxwell-Ampere solver " );
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
            }
        } else if( params.geometry == "2Dcartesian" ) {
            if( params.is_pxr == false ) {
            
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
                }
            } else if( params.is_spectral == true ) {
                solver = new NullSolver( params );
            }
            
        } else if( params.geometry == "3Dcartesian" ) {
            if( params.is_pxr == false ) {
                if( params.maxwell_sol == "Yee" ) {
                    solver = new MF_Solver3D_Yee( params );
                } else if( params.maxwell_sol == "Lehe" ) {
                    solver = new MF_Solver3D_Lehe( params );
                } else if( params.maxwell_sol == "Bouchard" ) {
                    solver = new MF_Solver3D_Bouchard( params );
                }
            } else {
                solver = new NullSolver( params );
            }
        } else if( params.geometry == "AMcylindrical" ) {
            if( params.is_pxr == false ) {
                if( params.maxwell_sol == "Yee" ) {
                    solver = new MF_SolverAM_Yee( params );
                }
            }
            else
                solver = new NullSolver( params );
        }
        
        if( !solver ) {
            ERROR( "Unknwon solver '" << params.maxwell_sol << "' for geometry '" << params.geometry <<"'" );
        }
        
        return solver;
    };
    
};

#endif
