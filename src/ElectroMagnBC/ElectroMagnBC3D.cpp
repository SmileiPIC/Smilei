
#include "ElectroMagnBC3D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"


ElectroMagnBC3D::ElectroMagnBC3D( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC( params, patch, _min_max )
{
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids)
        n_space = params.n_space_region;

    std::vector<unsigned int> oversize(params.oversize);
    if (params.uncoupled_grids)
        oversize = params.region_oversize;

    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = nx_p+1-params.is_pxr;
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = ny_p+1-params.is_pxr;
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = n_space[2]+1+2*oversize[2];
    nz_d = nz_p+1-params.is_pxr;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params.cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = params.cell_length[1];
    dt_ov_dy = dt/dy;
    dy_ov_dt = 1.0/dt_ov_dy;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the z-direction)
    dz       = params.cell_length[2];
    dt_ov_dz = dt/dz;
    dz_ov_dt = 1.0/dt_ov_dz;
    
}

ElectroMagnBC3D::~ElectroMagnBC3D()
{
}


void ElectroMagnBC3D::applyBConEdges( ElectroMagn *EMfields, Patch *patch )
{
    // Static cast of the fields
    //Field3D* Ex3D = static_cast<Field3D*>(EMfields->Ex_);
    //Field3D* Ey3D = static_cast<Field3D*>(EMfields->Ey_);
    //Field3D* Ez3D = static_cast<Field3D*>(EMfields->Ez_);
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
    
    double one_ov_dbeta ;
    
    if( patch->isXmin() ) {
    
        unsigned int i = 0;
        if( patch->isYmin() ) {
            // Xmin/Ymin
            // edge 0 : By[0,0,k] + beta(-x)Bx[0,0,k] = S(-x)
            // edge 8 : Bx[0,0,k] + beta(-y)By[0,0,k] = S(-y)
            if( EMfields->beta_edge[0]*EMfields->beta_edge[8] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[0]*EMfields->beta_edge[8] );
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    ( *By3D )( i, 0, k ) = ( EMfields->S_edge[0][k] - EMfields->beta_edge[0]* EMfields->S_edge[8][k] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, 0, k ) =   EMfields->S_edge[8][k] - EMfields->beta_edge[8]*( *By3D )( i, 0, k ) ;
                }
            } else {
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    //(*By3D)(i,0,k  ) = (EMfields->S_edge[0][k] + EMfields->S_edge[8][k])/2./(1+EMfields->beta_edge[8]);
                    //(*Bx3D)(i,0,k  ) = (EMfields->S_edge[0][k] + EMfields->S_edge[8][k])/2./(1+EMfields->beta_edge[0]);
                    ( *By3D )( i, 0, k ) = 0.;
                    ( *Bx3D )( i, 0, k ) = 0.;
                }
            }
        }//End Xmin Ymin edge
        
        if( patch->isYmax() ) {
            // Xmin/Ymax
            //edge 1 : By[0,ny_p-1,k] + beta(-x)Bx[0,ny_p,k] = S(-x)
            //edge 12 : Bx[0,ny_p,k] + beta(-y)By[0,ny_p-1,k] = S(-y)
            if( EMfields->beta_edge[1]*EMfields->beta_edge[12] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[1]*EMfields->beta_edge[12] );
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    ( *By3D )( i, ny_p-1, k ) = ( EMfields->S_edge[1][k] - EMfields->beta_edge[1]* EMfields->S_edge[12][k] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, ny_p, k ) =   EMfields->S_edge[12][k] - EMfields->beta_edge[12]*( *By3D )( i, ny_p-1, k ) ;
                }
            } else {
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    //(*By3D)(i,ny_p-1,k  ) = (EMfields->S_edge[1][k] + EMfields->S_edge[12][k])/2./(1+EMfields->beta_edge[12]);
                    //(*Bx3D)(i,ny_p,k    ) = (EMfields->S_edge[1][k] + EMfields->S_edge[12][k])/2./(1+EMfields->beta_edge[1]);
                    ( *By3D )( i, ny_p-1, k ) = 0.;
                    ( *Bx3D )( i, ny_p, k ) = 0.;
                }
            }
        }// End Xmin Ymax edge
        
        if( patch->isZmin() ) {
            // Xmin/Zmin
            // edge 2  : Bz[0,j,0] + beta(-x)Bx[0,j,0] = S(-x)
            // edge 16 : Bx[0,j,0] + beta(-z)Bz[0,j,0] = S(-z)
            if( EMfields->beta_edge[2]*EMfields->beta_edge[16] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[2]*EMfields->beta_edge[16] );
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    ( *Bz3D )( i, j, 0 ) = ( EMfields->S_edge[2][j] - EMfields->beta_edge[2]* EMfields->S_edge[16][j] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, j, 0 ) =   EMfields->S_edge[16][j] - EMfields->beta_edge[16]*( *Bz3D )( i, j, 0 ) ;
                }
            } else {
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    //(*Bz3D)(i,j,0  ) = (EMfields->S_edge[2][j] + EMfields->S_edge[16][j])/2./(1+EMfields->beta_edge[16]);
                    //(*Bx3D)(i,j,0  ) = (EMfields->S_edge[2][j] + EMfields->S_edge[16][j])/2./(1+EMfields->beta_edge[2]);
                    ( *Bz3D )( i, j, 0 ) = 0.;
                    ( *Bx3D )( i, j, 0 ) = 0.;
                }
            }
        } // End Xmin Zmin edge
        
        if( patch->isZmax() ) {
            // Xmin/Zmax
            // edge 3  : Bz[0,j,nz_p-1] + beta(-x)Bx[0,j,nz_p] = S(-x)
            // edge 20 : Bx[0,j,nz_p] + beta(+z)Bz[0,j,nz_p-1] = S(+z)
            if( EMfields->beta_edge[3]*EMfields->beta_edge[20] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[3]*EMfields->beta_edge[20] );
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    ( *Bz3D )( i, j, nz_p-1 ) = ( EMfields->S_edge[3][j] - EMfields->beta_edge[3]* EMfields->S_edge[20][j] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, j, nz_p ) =   EMfields->S_edge[20][j] - EMfields->beta_edge[20]*( *Bz3D )( i, j, nz_p-1 ) ;
                }
            } else {
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    //(*Bz3D)(i,j,nz_p-1  ) = (EMfields->S_edge[3][j] + EMfields->S_edge[20][j])/2./(1+EMfields->beta_edge[20]);
                    //(*Bx3D)(i,j,nz_p    ) = (EMfields->S_edge[3][j] + EMfields->S_edge[20][j])/2./(1+EMfields->beta_edge[3]);
                    ( *Bz3D )( i, j, nz_p-1 ) = 0.;
                    ( *Bx3D )( i, j, nz_p ) = 0.;
                }
            }
        }//End Xmin/Zmax edge
    } //End series of Xmin edges
    
    if( patch->isXmax() ) {
    
        unsigned int i = nx_p - 1;
        if( patch->isYmin() ) {
            // Xmax/Ymin
            // edge 4 : By[nx_p,0,k] + beta(+x)Bx[nx_p-1,0,k] = S(+x)
            // edge 9 : Bx[nx_p-1,0,k] + beta(-y)By[nx_p-1,0,k] = S(-y)
            if( EMfields->beta_edge[4]*EMfields->beta_edge[9] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[4]*EMfields->beta_edge[9] );
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    ( *By3D )( i+1, 0, k ) = ( EMfields->S_edge[4][k] - EMfields->beta_edge[4]* EMfields->S_edge[9][k] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, 0, k ) =   EMfields->S_edge[9][k] - EMfields->beta_edge[9]*( *By3D )( i+1, 0, k ) ;
                }
            } else {
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    //(*By3D)(i+1,0,k  ) = (EMfields->S_edge[4][k] + EMfields->S_edge[9][k])/2./(1+EMfields->beta_edge[9]);
                    //(*Bx3D)(i,0,k    ) = (EMfields->S_edge[4][k] + EMfields->S_edge[9][k])/2./(1+EMfields->beta_edge[4]);
                    ( *By3D )( i+1, 0, k ) = 0.;
                    ( *Bx3D )( i, 0, k ) = 0.;
                }
            }
        }//End Xmax Ymin edge
        
        if( patch->isYmax() ) {
            // Xmax/Ymax
            //edge 5 :  By[nx_p,ny_p-1,k] + beta(+x)Bx[nx_p-1,ny_p,k] = S(+x)
            //edge 13 : Bx[nx_p-1,ny_p,k] + beta(-y)By[nx_p,ny_p-1,k] = S(-y)
            if( EMfields->beta_edge[5]*EMfields->beta_edge[13] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[5]*EMfields->beta_edge[13] );
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    ( *By3D )( i+1, ny_p-1, k ) = ( EMfields->S_edge[5][k] - EMfields->beta_edge[5]* EMfields->S_edge[13][k] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, ny_p, k ) =   EMfields->S_edge[13][k] - EMfields->beta_edge[13]*( *By3D )( i+1, ny_p-1, k ) ;
                }
            } else {
                for( unsigned int k=patch->isZmin() ; k<nz_d-patch->isZmax() ; k++ ) {
                    //(*By3D)(i+1,ny_p-1,k  ) = (EMfields->S_edge[5][k] + EMfields->S_edge[13][k])/2./(1+EMfields->beta_edge[13]);
                    //(*Bx3D)(i,ny_p,k      ) = (EMfields->S_edge[5][k] + EMfields->S_edge[13][k])/2./(1+EMfields->beta_edge[5]);
                    ( *By3D )( i+1, ny_p-1, k ) = 0.;
                    ( *Bx3D )( i, ny_p, k ) = 0.;
                }
            }
        }// End Xmax Ymax edge
        
        if( patch->isZmin() ) {
            // Xmax/Zmin
            // edge 6  : Bz[nx_p,j,0] + beta(+x)Bx[nx_p-1,j,0] = S(+x)
            // edge 17 : Bx[nx_p-1,j,0] + beta(-z)Bz[nx_p,j,0] = S(-z)
            if( EMfields->beta_edge[6]*EMfields->beta_edge[17] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[6]*EMfields->beta_edge[17] );
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    ( *Bz3D )( i+1, j, 0 ) = ( EMfields->S_edge[6][j] - EMfields->beta_edge[6]* EMfields->S_edge[17][j] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, j, 0 ) =   EMfields->S_edge[17][j] - EMfields->beta_edge[17]*( *Bz3D )( i+1, j, 0 ) ;
                }
            } else {
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    //(*Bz3D)(i+1,j,0  ) = (EMfields->S_edge[6][j] + EMfields->S_edge[17][j])/2./(1+EMfields->beta_edge[17]);
                    //(*Bx3D)(i,j,0    ) = (EMfields->S_edge[6][j] + EMfields->S_edge[17][j])/2./(1+EMfields->beta_edge[6]);
                    ( *Bz3D )( i+1, j, 0 ) = 0.;
                    ( *Bx3D )( i, j, 0 ) = 0.;
                }
            }
        } // End Xmax Zmin edge
        
        if( patch->isZmax() ) {
            // Xmax/Zmax
            // edge 7  : Bz[nx_p,j,nz_p-1] + beta(+x)Bx[nx_p-1,j,nz_p] = S(+x)
            // edge 21 : Bx[nx_p-1,j,nz_p] + beta(+z)Bz[nx_p,j,nz_p-1] = S(+z)
            if( EMfields->beta_edge[7]*EMfields->beta_edge[21] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[7]*EMfields->beta_edge[21] );
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    ( *Bz3D )( i+1, j, nz_p-1 ) = ( EMfields->S_edge[7][j] - EMfields->beta_edge[7]* EMfields->S_edge[21][j] ) * one_ov_dbeta ;
                    ( *Bx3D )( i, j, nz_p ) =   EMfields->S_edge[21][j] - EMfields->beta_edge[21]*( *Bz3D )( i+1, j, nz_p-1 ) ;
                }
            } else {
                for( unsigned int j=patch->isYmin() ; j<ny_d-patch->isYmax() ; j++ ) {
                    //(*Bz3D)(i+1,j,nz_p-1  ) = (EMfields->S_edge[7][j] + EMfields->S_edge[21][j])/2./(1+EMfields->beta_edge[21]);
                    //(*Bx3D)(i,j,nz_p      ) = (EMfields->S_edge[7][j] + EMfields->S_edge[21][j])/2./(1+EMfields->beta_edge[7]);
                    ( *Bz3D )( i+1, j, nz_p-1 ) = 0.;
                    ( *Bx3D )( i, j, nz_p ) = 0.;
                }
            }
        }//End Xmax/Zmax edge
    } //End series of Xmax edges
    
    if( patch->isYmin() ) {
        unsigned int j = 0;
        if( patch->isZmin() ) { //Ymin/Zmin
            // edge 10 : Bz[i,0,0] + beta(-y)By[i,0,0] = S(-y)
            // edge 18 : By[i,0,0] + beta(-z)Bz[i,0,0] = S(-z)
            if( EMfields->beta_edge[10]*EMfields->beta_edge[18] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[10]*EMfields->beta_edge[18] );
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    ( *Bz3D )( i, j, 0 ) = ( EMfields->S_edge[10][i] - EMfields->beta_edge[10]* EMfields->S_edge[18][i] ) * one_ov_dbeta ;
                    ( *By3D )( i, j, 0 ) =   EMfields->S_edge[18][i] - EMfields->beta_edge[18]*( *Bz3D )( i, j, 0 ) ;
                }
            } else {
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    //(*Bz3D)(i,j,0   ) = (EMfields->S_edge[10][i] + EMfields->S_edge[18][i])/2./(1+EMfields->beta_edge[18]);
                    //(*By3D)(i,j,0   ) = (EMfields->S_edge[10][i] + EMfields->S_edge[18][i])/2./(1+EMfields->beta_edge[10]);
                    ( *Bz3D )( i, j, 0 ) = 0.;
                    ( *By3D )( i, j, 0 ) = 0.;
                }
            }
        }
        
        if( patch->isZmax() ) {
            // Ymin/Zmax
            //edge 11 : Bz[i,0,nz_p-1] + beta(-y)By[i,0,nz_p]   = S(-y)
            //edge 22 : By[i,0,nz_p]   + beta(+z)Bz[i,0,nz_p-1] = S(+z)
            if( EMfields->beta_edge[11]*EMfields->beta_edge[22] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[11]*EMfields->beta_edge[22] );
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    ( *Bz3D )( i, j, nz_p-1 ) = ( EMfields->S_edge[11][i] - EMfields->beta_edge[11]* EMfields->S_edge[22][i] ) * one_ov_dbeta ;
                    ( *By3D )( i, j, nz_p ) =   EMfields->S_edge[22][i] - EMfields->beta_edge[22]*( *Bz3D )( i, j, nz_p-1 ) ;
                    
                }
            } else {
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    //(*Bz3D)(i,j,nz_p-1   ) = (EMfields->S_edge[11][i] + EMfields->S_edge[22][i])/2./(1+EMfields->beta_edge[22]);
                    //(*By3D)(i,j,nz_p     ) = (EMfields->S_edge[11][i] + EMfields->S_edge[22][i])/2./(1+EMfields->beta_edge[11]);
                    ( *Bz3D )( i, j, nz_p-1 ) = 0.;
                    ( *By3D )( i, j, nz_p ) = 0.;
                }
            }
        } //End Ymin /Zmax edge
    } //End series of Ymin edges
    
    if( patch->isYmax() ) {
        unsigned int j = ny_p - 1;
        if( patch->isZmin() ) { //Ymax/Zmin
            // edge 14 : Bz[i,ny_p,0] + beta(+y)By[i,ny_p-1,0] = S(+y)
            // edge 19 : By[i,ny_p-1,0] + beta(-z)Bz[i,ny_p,0] = S(-z)
            if( EMfields->beta_edge[14]*EMfields->beta_edge[19] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[14]*EMfields->beta_edge[19] );
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    ( *Bz3D )( i, j+1, 0 ) = ( EMfields->S_edge[14][i] - EMfields->beta_edge[14]* EMfields->S_edge[19][i] ) * one_ov_dbeta ;
                    ( *By3D )( i, j, 0 )   =   EMfields->S_edge[19][i] - EMfields->beta_edge[19]*( *Bz3D )( i, j+1, 0 ) ;
                }
            } else {
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    //(*Bz3D)(i,j+1,0  ) = (EMfields->S_edge[14][i] + EMfields->S_edge[19][i])/2./(1+EMfields->beta_edge[19]);
                    //(*By3D)(i,j,0  )   = (EMfields->S_edge[14][i] + EMfields->S_edge[19][i])/2./(1+EMfields->beta_edge[14]);
                    ( *Bz3D )( i, j+1, 0 ) = 0.;
                    ( *By3D )( i, j, 0 )   = 0.;
                }
            }
        }//End Ymax /Zmin edge
        
        if( patch->isZmax() ) {
            // Ymax/Zmax
            //edge 15 : Bz[i,ny_p,nz_p-1] + beta(+y)By[i,ny_p-1,nz_p]   = S(+y)
            //edge 23 : By[i,ny_p-1,nz_p]   + beta(+z)Bz[i,ny_p,nz_p-1] = S(+z)
            if( EMfields->beta_edge[15]*EMfields->beta_edge[23] != 1.0 ) {
                one_ov_dbeta = 1./( 1. - EMfields->beta_edge[15]*EMfields->beta_edge[23] );
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    ( *Bz3D )( i, j+1, nz_p-1 ) = ( EMfields->S_edge[15][i] - EMfields->beta_edge[15]* EMfields->S_edge[23][i] ) * one_ov_dbeta ;
                    ( *By3D )( i, j, nz_p ) =   EMfields->S_edge[23][i] - EMfields->beta_edge[23]*( *Bz3D )( i, j+1, nz_p-1 ) ;
                    
                }
            } else {
                for( unsigned int i=patch->isXmin() ; i<nx_d-patch->isXmax() ; i++ ) {
                    //(*Bz3D)(i,j+1,nz_p-1  ) = (EMfields->S_edge[15][i] + EMfields->S_edge[23][i])/2./(1+EMfields->beta_edge[23]);
                    //(*By3D)(i,j,nz_p    )   = (EMfields->S_edge[15][i] + EMfields->S_edge[23][i])/2./(1+EMfields->beta_edge[15]);
                    ( *Bz3D )( i, j+1, nz_p-1 ) = 0.;
                    ( *By3D )( i, j, nz_p )   = 0.;
                }
            }
        } //End Ymax /Zmax edge
    } //End series of Ymax edges
}
