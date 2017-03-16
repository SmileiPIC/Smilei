#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include "MA_Solver1D_norm.h"
#include "MA_Solver2D_norm.h"
#include "MA_Solver2D_Friedman.h"
#include "MA_Solver3D_norm.h"
#include "MF_Solver1D_Yee.h"
#include "MF_Solver2D_Yee.h"
#include "MF_Solver3D_Yee.h"
#include "MF_Solver2D_Grassi.h"
#include "MF_Solver2D_GrassiSpL.h"
#include "MF_Solver2D_Cowan.h"
#include "MF_Solver2D_Lehe.h"

#include "Params.h"

#include "Tools.h"

class SolverFactory {
public:
    
    // create Maxwell-Ampere solver
    // -----------------------------
    static Solver* createMA(Params& params) {
        Solver* solver = NULL;
        DEBUG(params.maxwell_sol);
        
        if ( params.geometry == "1d3v" ) {
            solver = new MA_Solver1D_norm(params);
        } else if ( params.geometry == "2d3v" ) {
            if (params.Friedman_filter) {
                solver = new MA_Solver2D_Friedman(params);
            } else {
                solver = new MA_Solver2D_norm(params);
            }
        } else if ( params.geometry == "3d3v" ) {
            solver = new MA_Solver3D_norm(params);
        }
        
        if (!solver)
            ERROR( "Unknwon Maxwell-Ampere solver ");
        
        return solver;
    };
    
    // Create Maxwell-Faraday solver
    // -----------------------------
    static Solver* createMF(Params& params) {
        Solver* solver = NULL;
        DEBUG(params.maxwell_sol);
        
        // Create the required solver for Faraday's Equation
        // -------------------------------------------------
        if ( params.geometry == "1d3v" ) {
            if (params.maxwell_sol == "Yee") {
                solver = new MF_Solver1D_Yee(params);
            }
        } else if ( params.geometry == "2d3v" ) {
            if (params.maxwell_sol == "Yee") {
                solver = new MF_Solver2D_Yee(params);
            } else if (params.maxwell_sol == "Grassi") {
                solver = new MF_Solver2D_Grassi(params);
            } else if (params.maxwell_sol == "GrassiSpL") {
                solver = new MF_Solver2D_GrassiSpL(params);
            } else if (params.maxwell_sol == "Cowan") {
                solver = new MF_Solver2D_Cowan(params);
            } else if(params.maxwell_sol == "Lehe" ){
                solver = new MF_Solver2D_Lehe(params);
            }
            
        } else if ( params.geometry == "3d3v" ) {
            if (params.maxwell_sol == "Yee") {
                solver = new MF_Solver3D_Yee(params);
            }
        }
        
        if (!solver)
            ERROR( "Unknwon solver '" << params.maxwell_sol << "' for geometry '" << params.geometry <<"'" );
        
        return solver;
    };
    
};

#endif
