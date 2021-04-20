#include "ElectroMagnBC3D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Tools.h"

using namespace std;


ElectroMagnBC3D_refl::ElectroMagnBC3D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC3D( params, patch, i_boundary )
{
    // oversize
    if (!params.multiple_decomposition) {
        oversize_x = params.oversize[0];
        oversize_y = params.oversize[1];
        oversize_z = params.oversize[2];
    }
    else {
        oversize_x = params.region_oversize[0];
        oversize_y = params.region_oversize[1];
        oversize_z = params.region_oversize[2];
    }
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_refl::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if( i_boundary_ == 0 && patch->isXmin() ) {
    
        // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
        
        // Static cast of the fields
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        // for By^(d,p,d)
        for( unsigned int i=oversize_x; i>0; i-- ) {
            for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
                for( unsigned int k=0 ; k<n_d[2] ; k++ ) {
                    ( *By3D )( i-1, j, k ) = ( *By3D )( i, j, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
        for( unsigned int i=oversize_x; i>0; i-- ) {
            for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
                for( unsigned int k=0 ; k<n_p[2] ; k++ ) {
                    ( *Bz3D )( i-1, j, k ) = ( *Bz3D )( i, j, k );
                }
            }
        }
        
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
    
        // Static cast of the fields
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        // for By^(d,p,d)
        for( unsigned int i=n_d[0]-oversize_x; i<n_d[0]; i++ ) {
            for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
                for( unsigned int k=0 ; k<n_d[2] ; k++ ) {
                    ( *By3D )( i, j, k ) = ( *By3D )( i-1, j, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
        for( unsigned int i=n_d[0]-oversize_x; i<n_d[0]; i++ ) {
            for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
                for( unsigned int k=0 ; k<n_p[2] ; k++ ) {
                    ( *Bz3D )( i, j, k ) = ( *Bz3D )( i-1, j, k );
                }
            }
        }
        
    } else if( i_boundary_ == 2 && patch->isYmin() ) {
    
        // Static cast of the fields
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for( unsigned int i=0; i<n_p[0]; i++ ) {
            for( unsigned int j=oversize_y ; j>0 ; j-- ) {
                for( unsigned int k=0; k<n_d[2] ; k++ ) {
                    ( *Bx3D )( i, j-1, k ) = ( *Bx3D )( i, j, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
        for( unsigned int i=0; i<n_d[0]; i++ ) {
            for( unsigned int j=oversize_y ; j>0 ; j-- ) {
                for( unsigned int k=0; k<n_p[2] ; k++ ) {
                    ( *Bz3D )( i, j-1, k ) = ( *Bz3D )( i, j, k );
                }
            }
        }
        
    } else if( i_boundary_ == 3 && patch->isYmax() ) {
    
        // Static cast of the fields
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for( unsigned int i=0; i<n_p[0]; i++ ) {
            for( unsigned int j=n_d[1]-oversize_y; j<n_d[1] ; j++ ) {
                for( unsigned int k=0; k<n_d[2] ; k++ ) {
                    ( *Bx3D )( i, j, k ) = ( *Bx3D )( i, j-1, k );
                }
            }
        }
        
        // for Bz^(d,d,p)
        for( unsigned int i=0; i<n_d[0]; i++ ) {
            for( unsigned int j=n_d[1]-oversize_y; j<n_d[1] ; j++ ) {
                for( unsigned int k=0; k<n_p[2] ; k++ ) {
                    ( *Bz3D )( i, j, k ) = ( *Bz3D )( i, j-1, k );
                }
            }
        }
        
    } else if( i_boundary_==4 && patch->isZmin() ) {
    
        // Static cast of the fields
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for( unsigned int i=0; i<n_p[0]; i++ ) {
            for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
                for( unsigned int k=oversize_z ; k>0 ; k-- ) {
                    ( *Bx3D )( i, j, k-1 ) = ( *Bx3D )( i, j, k );
                }
            }
        }
        
        // for By^(d,p,d)
        for( unsigned int i=0; i<n_d[0]; i++ ) {
            for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
                for( unsigned int k=oversize_z ; k>0 ; k-- ) {
                    ( *By3D )( i, j, k-1 ) = ( *By3D )( i, j, k );
                }
            }
        }
        
    } else if( i_boundary_==5 && patch->isZmax() ) {
    
        // Static cast of the fields
        Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
        Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
        
        // FORCE CONSTANT MAGNETIC FIELDS
        
        // for Bx^(p,d,d)
        for( unsigned int i=0; i<n_p[0]; i++ ) {
            for( unsigned int j=0 ; j<n_d[1] ; j++ ) {
                for( unsigned int k=n_d[2]-oversize_z; k<n_d[2] ; k++ ) {
                    ( *Bx3D )( i, j, k ) = ( *Bx3D )( i, j, k-1 );
                }
            }
        }
        
        // for By^(d,p,d)
        for( unsigned int i=0; i<n_d[0]; i++ ) {
            for( unsigned int j=0 ; j<n_p[1] ; j++ ) {
                for( unsigned int k=n_d[2]-oversize_z; k<n_d[2] ; k++ ) {
                    ( *By3D )( i, j, k ) = ( *By3D )( i, j, k-1 );
                }
            }
        }
        
    }
}

