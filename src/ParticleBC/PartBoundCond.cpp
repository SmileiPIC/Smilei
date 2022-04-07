#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "Patch.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond( Params &params, Species *species, Patch *patch ) :
    isAM( params.geometry == "AMcylindrical" )
{
    // number of dimensions for the particle
    //!\todo (MG to JD) isn't it always 3?
    nDim_particle = params.nDim_particle;
    nDim_field = params.nDim_field;
    // Absolute global values
    double x_min_global = 0;
    double x_max_global = params.cell_length[0]*( params.n_space_global[0] );
    double y_min_global = 0;
    double y_max_global = params.cell_length[1]*( params.n_space_global[1] );
    double z_min_global = 0;
    double z_max_global = params.cell_length[2]*( params.n_space_global[2] );
    
    bc_xmin  = NULL;
    bc_xmax  = NULL;
    bc_ymin  = NULL;
    bc_ymax  = NULL;
    bc_zmin  = NULL;
    bc_zmax  = NULL;

    dt_ = params.timestep;
    
    // -----------------------------
    // Define limits of local domain
    if( params.EM_BCs[0][0]=="periodic" || params.hasWindow ) {
        x_min = patch->getDomainLocalMin( 0 );
        x_max = patch->getDomainLocalMax( 0 );
    } else {
        x_min = max( x_min_global, patch->getDomainLocalMin( 0 ) );
        x_max = min( x_max_global, patch->getDomainLocalMax( 0 ) );
    }
    
    if( nDim_particle > 1 ) {
        if( params.EM_BCs[1][0]=="periodic" ) {
            y_min = patch->getDomainLocalMin( 1 );
            y_max = patch->getDomainLocalMax( 1 );
        } else {
            y_min = max( y_min_global, patch->getDomainLocalMin( 1 ) );
            y_max = min( y_max_global, patch->getDomainLocalMax( 1 ) );
        }
        
        if( ( nDim_particle > 2 ) && ( !isAM ) ) {
            if( params.EM_BCs[2][0]=="periodic" ) {
                z_min = patch->getDomainLocalMin( 2 );
                z_max = patch->getDomainLocalMax( 2 );
            } else {
                z_min = max( z_min_global, patch->getDomainLocalMin( 2 ) );
                z_max = min( z_max_global, patch->getDomainLocalMax( 2 ) );
            }
        }
    }
    
    // Can be done after parsing
    // Check for inconsistencies between EM and particle BCs
    if( ! species->particles->tracked ) {
        for( unsigned int iDim=0; iDim<( unsigned int )nDim_field; iDim++ ) {
            if(     (  ( ( params.EM_BCs[iDim][0]=="periodic" )&&( species->boundary_conditions_[iDim][0]!="periodic" ) )
                    || ( ( params.EM_BCs[iDim][1]=="periodic" )&&( species->boundary_conditions_[iDim][1]!="periodic" ) ))
                 && (params.is_spectral==false)) {
                WARNING( "For species " << species->name_ << ", periodic EM "<<"xyz"[iDim]<<"-boundary conditions require particle BCs to be periodic." );
            }
        }
    }
    // ----------------------------------------------
    // Define the kind of applied boundary conditions
    // ----------------------------------------------
    
    void ( *remove_inf )( Species*, int, int, int, double, double, std::vector<double>&, Random *rand, double& );
    if( species->mass_ == 0 ) {
        remove_inf = &remove_photon_inf;
    } else {
        remove_inf = &remove_particle_inf;
    }
    void ( *remove_sup)( Species*, int, int, int, double, double, std::vector<double>&, Random *rand, double& );
    if( species->mass_ == 0 ) {
        remove_sup = &remove_photon_sup;
    } else {
        remove_sup = &remove_particle_sup;
    }
    // Xmin
    bc_xmin = &internal_inf;
    if( species->boundary_conditions_[0][0] == "reflective" ) {
        if( patch->isXmin() ) {
            bc_xmin = &reflect_particle_inf;
        }
    } else if( species->boundary_conditions_[0][0] == "remove" ) {
        if( patch->isXmin() ) {
            bc_xmin = remove_inf;
        }
    } else if( species->boundary_conditions_[0][0] == "stop" ) {
        if( patch->isXmin() ) {
            bc_xmin = &stop_particle_inf;
        }
    } else if( species->boundary_conditions_[0][0] == "thermalize" ) {
        if( patch->isXmin() ) {
            bc_xmin = &thermalize_particle_inf;
        }
    } else if( species->boundary_conditions_[0][0] == "periodic" ) {
        // Nothing to do
    } else {
        ERROR( "Xmin boundary condition `"<<species->boundary_conditions_[0][0]<<"` unknown" );
    }
    
    // Xmax
    bc_xmax = &internal_sup;
    if( species->boundary_conditions_[0][1] == "reflective" ) {
        if( patch->isXmax() ) {
            bc_xmax = &reflect_particle_sup;
        }
    } else if( species->boundary_conditions_[0][1] == "remove" ) {
        if( patch->isXmax() ) {
            bc_xmax = remove_sup;
        }
    } else if( species->boundary_conditions_[0][1] == "stop" ) {
        if( patch->isXmax() ) {
            bc_xmax = &stop_particle_sup;
        }
    } else if( species->boundary_conditions_[0][1] == "thermalize" ) {
        if( patch->isXmax() ) {
            bc_xmax = &thermalize_particle_sup;
        }
    } else if( species->boundary_conditions_[0][1] == "periodic" ) {
        // Nothing to do
    } else {
        ERROR( "Xmax boundary condition `"<<species->boundary_conditions_[0][1]<<"`  unknown" );
    }
    
    
    if( ( nDim_particle > 1 ) && ( !isAM ) ) {
        // Ymin
        bc_ymin = &internal_inf;
        if( species->boundary_conditions_[1][0] == "reflective" ) {
            if( patch->isYmin() ) {
                bc_ymin = &reflect_particle_inf;
            }
        } else if( species->boundary_conditions_[1][0] == "remove" ) {
            if( patch->isYmin() ) {
                bc_ymin = remove_inf;
            }
        } else if( species->boundary_conditions_[1][0] == "stop" ) {
            if( patch->isYmin() ) {
                bc_ymin = &stop_particle_inf;
            }
        } else if( species->boundary_conditions_[1][0] == "thermalize" ) {
            if( patch->isYmin() ) {
                bc_ymin = &thermalize_particle_inf;
            }
        } else if( species->boundary_conditions_[1][0] == "periodic" ) {
            // Nothing to do
        } else {
            ERROR( "Ymin boundary condition `"<< species->boundary_conditions_[1][0] << "` unknown" );
        }
        
        // Ymax
        bc_ymax = &internal_sup;
        if( species->boundary_conditions_[1][1] == "reflective" ) {
            if( patch->isYmax() ) {
                bc_ymax = &reflect_particle_sup;
            }
        } else if( species->boundary_conditions_[1][1] == "remove" ) {
            if( patch->isYmax() ) {
                bc_ymax = remove_sup;
            }
        } else if( species->boundary_conditions_[1][1] == "stop" ) {
            if( patch->isYmax() ) {
                bc_ymax = &stop_particle_sup;
            }
        } else if( species->boundary_conditions_[1][1] == "thermalize" ) {
            if( patch->isYmax() ) {
                bc_ymax = &thermalize_particle_sup;
            }
        } else if( species->boundary_conditions_[1][1] == "periodic" ) {
            // Nothing to do
        } else {
            ERROR( "Ymax boundary condition `"<< species->boundary_conditions_[1][1] <<"` undefined" );
        }
        
        
        if( nDim_particle > 2 ) {
            bc_zmin = &internal_inf;
            if( species->boundary_conditions_[2][0] == "reflective" ) {
                if( patch->isZmin() ) {
                    bc_zmin = &reflect_particle_inf;
                }
            } else if( species->boundary_conditions_[2][0] == "remove" ) {
                if( patch->isZmin() ) {
                    bc_zmin = remove_inf;
                }
            } else if( species->boundary_conditions_[2][0] == "stop" ) {
                if( patch->isZmin() ) {
                    bc_zmin = &stop_particle_inf;
                }
            } else if( species->boundary_conditions_[2][0] == "thermalize" ) {
                if( patch->isZmin() ) {
                    bc_zmin = &thermalize_particle_inf;
                }
            } else if( species->boundary_conditions_[2][0] == "periodic" ) {
                // Nothing to do
            } else {
                ERROR( "Zmin boundary condition `"<< species->boundary_conditions_[2][0] << "` unknown" );
            }
            
            bc_zmax = &internal_sup;
            if( species->boundary_conditions_[2][1] == "reflective" ) {
                if( patch->isZmax() ) {
                    bc_zmax = &reflect_particle_sup;
                }
            } else if( species->boundary_conditions_[2][1] == "remove" )  {
                if( patch->isZmax() ) {
                    bc_zmax = remove_sup;
                }
            } else if( species->boundary_conditions_[2][1] == "stop" ) {
                if( patch->isZmax() ) {
                    bc_zmax = &stop_particle_sup;
                }
            } else if( species->boundary_conditions_[2][1] == "thermalize" ) {
                if( patch->isZmax() ) {
                    bc_zmax = &thermalize_particle_sup;
                }
            } else if( species->boundary_conditions_[2][1] == "periodic" ) {
                // Nothing to do
            } else {
                ERROR( "Zmax boundary condition `"<< species->boundary_conditions_[2][1] << "` unknown" );
            }
            
        }//nDim_particle>2
        
    }//nDim_particle>1
    else if( isAM ) {
        
        // Ymax
        bc_ymin = &internal_inf_AM;
        bc_ymax = &internal_sup_AM;
        if( species->boundary_conditions_[1][1] == "remove" ) {
            if( patch->isYmax() ) {
                bc_ymax = &remove_particle_AM;
            }
        } else if( species->boundary_conditions_[1][1] == "reflective" ) {
            if( patch->isYmax() ) {
                bc_ymax = &refl_particle_AM;
            }
        }
        else {
            ERROR_NAMELIST( 
                "Only Remove and reflective boundary conditions can be applied to particles in AM geometry ",
                LINK_NAMELIST + std::string("#species") );
        }
    }
}


PartBoundCond::~PartBoundCond()
{
}
