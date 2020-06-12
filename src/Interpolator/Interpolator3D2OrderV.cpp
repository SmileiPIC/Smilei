#include "Interpolator3D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "LaserEnvelope.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D2OrderV::Interpolator3D2OrderV( Params &params, Patch *patch ) : Interpolator3D( params, patch )
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    dz_inv_ = 1.0/params.cell_length[2];
    D_inv[0] = 1.0/params.cell_length[0];
    D_inv[1] = 1.0/params.cell_length[1];
    D_inv[2] = 1.0/params.cell_length[2];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd OrderV Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D2OrderV::fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc )
{
}

void Interpolator3D2OrderV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }
    
    int idxO[3];
    double idx[3];
    //Primal indices are the same for all particles
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = ( int )idx[0] - i_domain_begin ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = ( int )idx[1] - j_domain_begin ;
    idx[2]  = round( particles.position( 2, *istart ) * D_inv[2] );
    idxO[2] = ( int )idx[2] - k_domain_begin ;
    
    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_m );
    
    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );
    
    double *Epart[3], *Bpart[3];
    
    double coeff[3][2][3][32];
    int dual[3][32]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0, delta_primal < 0) or if it is +1 (dual=1, delta_primal>=0).
    
    int vecSize = 32;
    
    int cell_nparts( ( int )iend[0]-( int )istart[0] );
    int nbVec = ( iend[0]-istart[0]+( cell_nparts-1 )-( ( iend[0]-istart[0]-1 )&( cell_nparts-1 ) ) ) / vecSize;
    
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }
    
    for( int iivect=0 ; iivect<nbVec; iivect++ ) {
        int ivect = vecSize*iivect;
        
        int np_computed( 0 );
        if( cell_nparts > vecSize ) {
            np_computed = vecSize;
            cell_nparts -= vecSize;
        } else {
            np_computed = cell_nparts;
        }
        
        double *deltaO[3]; //Delta is the distance of the particle from its primal node in cell size. Delta is in [-0.5, +0.5[
        deltaO[0] = &( smpi->dynamics_deltaold[ithread][0        + ivect + istart[0] - ipart_ref] );
        deltaO[1] = &( smpi->dynamics_deltaold[ithread][nparts   + ivect + istart[0] - ipart_ref] );
        deltaO[2] = &( smpi->dynamics_deltaold[ithread][2*nparts + ivect + istart[0] - ipart_ref] );
        
        for( unsigned int k=0; k<3; k++ ) {
            Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts-ipart_ref+ivect+istart[0]] );
            Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts-ipart_ref+ivect+istart[0]] );
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double delta2, delta;
            
            for( int i=0; i<3; i++ ) { // for X/Y
                //delta primal = distance to primal node
                delta   = particles.position( i, ipart+ivect+istart[0] )*D_inv[i] - idx[i];
                delta2  = delta*delta;
                coeff[i][0][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
                coeff[i][0][1][ipart]    = ( 0.75 - delta2 );
                coeff[i][0][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
                //store delta primal in global array
                //deltaO[i][ipart-ipart_ref+ivect+istart[0]] = delta;
                deltaO[i][ipart] = delta;
                dual [i][ipart] = ( delta >= 0. );
                
                //delta dual = distance to dual node
                delta   = delta - dual[i][ipart] + 0.5 ;
                delta2  = delta*delta;
                
                coeff[i][1][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
                coeff[i][1][1][ipart]    = ( 0.75 - delta2 );
                coeff[i][1][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            }
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double *coeffyp = &( coeff[1][0][1][ipart] );
            double *coeffyd = &( coeff[1][1][1][ipart] );
            double *coeffxd = &( coeff[0][1][1][ipart] );
            double *coeffxp = &( coeff[0][0][1][ipart] );
            double *coeffzp = &( coeff[2][0][1][ipart] );
            double *coeffzd = &( coeff[2][1][1][ipart] );
            
            //Ex(dual, primal, primal)
            double interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[0][ipart] )*( *Ex3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0][ipart]*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
                    }
                }
            }
            Epart[0][ipart] = interp_res;
            
            //Ey(primal, dual, primal)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[1][ipart] )*( *Ey3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1][ipart]*( *Ey3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
                    }
                }
            }
            Epart[1][ipart] = interp_res;
            
            
            //Ez(primal, primal, dual)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] )*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2][ipart]*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
                    }
                }
            }
            Epart[2][ipart] = interp_res;
            
            
            //Bx(primal, dual , dual )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] ) * ( ( 1-dual[1][ipart] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1][ipart]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) )
                                        +    dual[2][ipart]  * ( ( 1-dual[1][ipart] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[1][ipart]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) ) );
                    }
                }
            }
            Bpart[0][ipart] = interp_res;
            
            //By(dual, primal, dual )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] ) * ( ( 1-dual[0][ipart] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0][ipart]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                        +    dual[2][ipart]  * ( ( 1-dual[0][ipart] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[0][ipart]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) ) );
                    }
                }
            }
            Bpart[1][ipart] = interp_res;
            
            //Bz(dual, dual, prim )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0][ipart]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                        +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) + dual[0][ipart]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) ) );
                    }
                }
            }
            Bpart[2][ipart] = interp_res;
        }
    }
} // END Interpolator3D2OrderV


void Interpolator3D2OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    // iend not used for now
    // probes are interpolated one by one for now
    
    int ipart = *istart;
    int nparts( particles.size() );
    
    
    double *Epart[3], *Bpart[3];
    for( unsigned int k=0; k<3; k++ ) {
        Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts] );
        Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts] );
    }
    
    int idxO[3];
    double idx[3];
    //Primal indices are the same for all particles
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = ( int )idx[0] - i_domain_begin ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = ( int )idx[1] - j_domain_begin ;
    idx[2]  = round( particles.position( 2, *istart ) * D_inv[2] );
    idxO[2] = ( int )idx[2] - k_domain_begin ;
    
    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_ );
    Field3D *Jx3D = static_cast<Field3D *>( EMfields->Jx_ );
    Field3D *Jy3D = static_cast<Field3D *>( EMfields->Jy_ );
    Field3D *Jz3D = static_cast<Field3D *>( EMfields->Jz_ );
    Field3D *rho3D = static_cast<Field3D *>( EMfields->rho_ );
    
    double coeff[3][2][3];
    int dual[3]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0, delta_primal < 0) or if it is +1 (dual=1, delta_primal>=0).
    
    double delta2, delta;
    
    for( int i=0; i<3; i++ ) { // for X/Y
        //delta primal = distance to primal node
        delta   = particles.position( i, ipart )*D_inv[i] - idx[i];
        delta2  = delta*delta;
        coeff[i][0][0]    =  0.5 * ( delta2-delta+0.25 );
        coeff[i][0][1]    = ( 0.75 - delta2 );
        coeff[i][0][2]    =  0.5 * ( delta2+delta+0.25 );
        dual [i] = ( delta >= 0. );
        
        //delta dual = distance to dual node
        delta   = delta - dual[i] + 0.5 ;
        delta2  = delta*delta;
        
        coeff[i][1][0]    =  0.5 * ( delta2-delta+0.25 );
        coeff[i][1][1]    = ( 0.75 - delta2 );
        coeff[i][1][2]    =  0.5 * ( delta2+delta+0.25 );
    }
    
    double *coeffyp = &( coeff[1][0][1] );
    double *coeffyd = &( coeff[1][1][1] );
    double *coeffxd = &( coeff[0][1][1] );
    double *coeffxp = &( coeff[0][0][1] );
    double *coeffzp = &( coeff[2][0][1] );
    double *coeffzd = &( coeff[2][1][1] );
    
    //Ex(dual, primal, primal)
    double interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[0] )*( *Ex3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
            }
        }
    }
    Epart[0][ipart] = interp_res;
    
    
    //Ey(primal, dual, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[1] )*( *Ey3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Ey3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
            }
        }
    }
    Epart[1][ipart] = interp_res;
    
    
    //Ez(primal, primal, dual)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] )*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2]*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
            }
        }
    }
    Epart[2][ipart] = interp_res;
    
    
    //Bx(primal, dual , dual )
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] ) * ( ( 1-dual[1] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) )
                                +    dual[2]  * ( ( 1-dual[1] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[1]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) ) );
            }
        }
    }
    Bpart[0][ipart] = interp_res;
    
    //By(dual, primal, dual )
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] ) * ( ( 1-dual[0] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                +    dual[2]  * ( ( 1-dual[0] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[0]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) ) );
            }
        }
    }
    Bpart[1][ipart] = interp_res;
    
    //Bz(dual, dual, prim )
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[1] ) * ( ( 1-dual[0] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                +    dual[1]  * ( ( 1-dual[0] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) + dual[0]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) ) );
            }
        }
    }
    Bpart[2][ipart] = interp_res;
    
    
    //Jx(dual, primal, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[0] )*( *Jx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Jx3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
            }
        }
    }
    JLoc->x = interp_res;
    
    
    //Jy(primal, dual, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[1] )*( *Jy3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Jy3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
            }
        }
    }
    JLoc->y = interp_res;
    
    
    //Jz(primal, primal, dual)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] )*( *Jz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2]*( *Jz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
            }
        }
    }
    JLoc->z = interp_res;
    
    //Rho(primal, primal, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) * ( *rho3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc );
            }
        }
    }
    ( *RhoLoc ) = interp_res;
    
}




void Interpolator3D2OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }
    
    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );
    
    double *Epart[3], *Bpart[3], *Phipart[1], *GradPhipart[3];
    
    double *deltaO[3];
    deltaO[0] = &( smpi->dynamics_deltaold[ithread][0] );
    deltaO[1] = &( smpi->dynamics_deltaold[ithread][nparts] );
    deltaO[2] = &( smpi->dynamics_deltaold[ithread][2*nparts] );
    
    for( unsigned int k=0; k<3; k++ ) {
    
        if( k==0 ) {   // scalar field, only one component
            Phipart[k]     = &( smpi->dynamics_PHIpart[ithread][k*nparts] );
        }
        
        Epart[k]       = &( smpi->dynamics_Epart[ithread][k*nparts] );
        Bpart[k]       = &( smpi->dynamics_Bpart[ithread][k*nparts] );
        GradPhipart[k] = &( smpi->dynamics_GradPHIpart[ithread][k*nparts] );
    }
    
    int idx[3], idxO[3];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;
    idx[2]  = round( particles.position( 2, *istart ) * D_inv[2] );
    idxO[2] = idx[2] - k_domain_begin -1 ;
    
    Field3D *Ex3D       = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D       = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D       = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D       = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D       = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D       = static_cast<Field3D *>( EMfields->Bz_m );
    Field3D *Phi3D      = static_cast<Field3D *>( EMfields->envelope->Phi_ );
    Field3D *GradPhix3D = static_cast<Field3D *>( EMfields->envelope->GradPhix_ );
    Field3D *GradPhiy3D = static_cast<Field3D *>( EMfields->envelope->GradPhiy_ );
    Field3D *GradPhiz3D = static_cast<Field3D *>( EMfields->envelope->GradPhiz_ );
    
    
    
    double coeff[3][2][3][32];
    int dual[3][32]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
    
    int vecSize = 32;
    
    int cell_nparts( ( int )iend[0]-( int )istart[0] );
    int nbVec = ( iend[0]-istart[0]+( cell_nparts-1 )-( ( iend[0]-istart[0]-1 )&( cell_nparts-1 ) ) ) / vecSize;
    
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }
    
    for( int iivect=0 ; iivect<nbVec; iivect++ ) {
        int ivect = vecSize*iivect;
        
        int np_computed( 0 );
        if( cell_nparts > vecSize ) {
            np_computed = vecSize;
            cell_nparts -= vecSize;
        } else {
            np_computed = cell_nparts;
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double delta0, delta;
            double delta2;
            
            
            for( int i=0; i<3; i++ ) { // for X/Y
                delta0 = particles.position( i, ipart+ivect+istart[0] )*D_inv[i];
                dual [i][ipart] = ( delta0 - ( double )idx[i] >=0. );
                
                for( int j=0; j<2; j++ ) { // for dual
                
                    delta   = delta0 - ( double )idx[i] + ( double )j*( 0.5-dual[i][ipart] );
                    delta2  = delta*delta;
                    
                    coeff[i][j][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
                    coeff[i][j][1][ipart]    = ( 0.75 - delta2 );
                    coeff[i][j][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
                    
                    if( j==0 ) {
                        deltaO[i][ipart-ipart_ref+ivect+istart[0]] = delta;
                    }
                }
            }
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double *coeffyp = &( coeff[1][0][1][ipart] );
            double *coeffyd = &( coeff[1][1][1][ipart] );
            double *coeffxd = &( coeff[0][1][1][ipart] );
            double *coeffxp = &( coeff[0][0][1][ipart] );
            double *coeffzp = &( coeff[2][0][1][ipart] );
            double *coeffzd = &( coeff[2][1][1][ipart] );
            
            //Ex(dual, primal, primal)
            double interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[0][ipart] )*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) + dual[0][ipart]*( *Ex3D )( idxO[0]+2+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) );
                    }
                }
            }
            Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            
            //Ey(primal, dual, primal)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[1][ipart] )*( *Ey3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) + dual[1][ipart]*( *Ey3D )( idxO[0]+1+iloc, idxO[1]+2+jloc, idxO[2]+1+kloc ) );
                    }
                }
            }
            Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            
            //Ez(primal, primal, dual)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        //interp_res += *(coeffxp+iloc*32) * *(coeffyd+jloc*32) * *(coeffzp+kloc*32) *
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] )*( *Ez3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) + dual[2][ipart]*( *Ez3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+2+kloc ) );
                    }
                }
            }
            Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            
            //Bx(primal, dual , dual )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        //interp_res += *(coeffxp+iloc*32) * *(coeffyd+jloc*32) * *(coeffzp+kloc*32) *
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] ) * ( ( 1-dual[1][ipart] )*( *Bx3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) + dual[1][ipart]*( *Bx3D )( idxO[0]+1+iloc, idxO[1]+2+jloc, idxO[2]+1+kloc ) )
                                        +    dual[2][ipart]  * ( ( 1-dual[1][ipart] )*( *Bx3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+2+kloc ) + dual[1][ipart]*( *Bx3D )( idxO[0]+1+iloc, idxO[1]+2+jloc, idxO[2]+2+kloc ) ) );
                    }
                }
            }
            Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            //By(dual, primal, dual )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        //interp_res += *(coeffxp+iloc*32) * *(coeffyd+jloc*32) * *(coeffzp+kloc*32) *
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] ) * ( ( 1-dual[0][ipart] )*( *By3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) + dual[0][ipart]*( *By3D )( idxO[0]+2+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) )
                                        +    dual[2][ipart]  * ( ( 1-dual[0][ipart] )*( *By3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+2+kloc ) + dual[0][ipart]*( *By3D )( idxO[0]+2+iloc, idxO[1]+1+jloc, idxO[2]+2+kloc ) ) );
                    }
                }
            }
            Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            //Bz(dual, dual, prim )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) + dual[0][ipart]*( *Bz3D )( idxO[0]+2+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) )
                                        +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+2+jloc, idxO[2]+1+kloc ) + dual[0][ipart]*( *Bz3D )( idxO[0]+2+iloc, idxO[1]+2+jloc, idxO[2]+1+kloc ) ) );
                    }
                }
            }
            Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of Phi^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *Phi3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            Phipart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of GradPhiX^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *GradPhix3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            GradPhipart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of GradPhiY^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *GradPhiy3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            GradPhipart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of GradPhiZ^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *GradPhiz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            GradPhipart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
        }
    }
    //}
    
    
}


void Interpolator3D2OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }
    
    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );
    
    double *Phi_mpart[1], *GradPhi_mpart[3];
    
    double *deltaO[3];
    deltaO[0] = &( smpi->dynamics_deltaold[ithread][0] );
    deltaO[1] = &( smpi->dynamics_deltaold[ithread][nparts] );
    deltaO[2] = &( smpi->dynamics_deltaold[ithread][2*nparts] );
    
    for( unsigned int k=0; k<3; k++ ) {
    
        if( k==0 ) {   // scalar field, only one component
            Phi_mpart[k]     = &( smpi->dynamics_PHI_mpart[ithread][k*nparts] );
        }
        
        GradPhi_mpart[k] = &( smpi->dynamics_GradPHI_mpart[ithread][k*nparts] );
    }
    
    int idx[3], idxO[3];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;
    idx[2]  = round( particles.position( 2, *istart ) * D_inv[2] );
    idxO[2] = idx[2] - k_domain_begin -1 ;
    
    Field3D *Phi_m3D      = static_cast<Field3D *>( EMfields->envelope->Phi_m );
    Field3D *GradPhi_mx3D = static_cast<Field3D *>( EMfields->envelope->GradPhix_m );
    Field3D *GradPhi_my3D = static_cast<Field3D *>( EMfields->envelope->GradPhiy_m );
    Field3D *GradPhi_mz3D = static_cast<Field3D *>( EMfields->envelope->GradPhiz_m );
    
    
    double coeff[3][2][3][32];
    int vecSize = 32;
    
    int cell_nparts( ( int )iend[0]-( int )istart[0] );
    int nbVec = ( iend[0]-istart[0]+( cell_nparts-1 )-( ( iend[0]-istart[0]-1 )&( cell_nparts-1 ) ) ) / vecSize;
    
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }
    
    for( int iivect=0 ; iivect<nbVec; iivect++ ) {
        int ivect = vecSize*iivect;
        
        int np_computed( 0 );
        if( cell_nparts > vecSize ) {
            np_computed = vecSize;
            cell_nparts -= vecSize;
        } else {
            np_computed = cell_nparts;
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double delta0, delta;
            double delta2;
            
            
            for( int i=0; i<3; i++ ) { // for X/Y
                delta0 = particles.position( i, ipart+ivect+istart[0] )*D_inv[i];
                delta   = delta0 - ( double )idx[i] ;
                delta2  = delta*delta;
                
                coeff[i][0][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
                coeff[i][0][1][ipart]    = ( 0.75 - delta2 );
                coeff[i][0][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
                
                deltaO[i][ipart-ipart_ref+ivect+istart[0]] = delta;
            }
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double *coeffyp = &( coeff[1][0][1][ipart] );
            double *coeffxp = &( coeff[0][0][1][ipart] );
            double *coeffzp = &( coeff[2][0][1][ipart] );
            
            // ----------------------------------
            // ----  PHI_m  and Grad Phi_m   ----
            // ----------------------------------
            // Interpolation of Phi_m^(p,p,p)
            
            double interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *Phi_m3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            Phi_mpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of GradPhi_mX^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *GradPhi_mx3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            GradPhi_mpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of GradPhi_mY^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *GradPhi_my3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            GradPhi_mpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            // Interpolation of GradPhi_mZ^(p,p,p)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) * ( *GradPhi_mz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
                    }
                }
            }
            GradPhi_mpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
        }
    }
    //}
    
    
}

// probes like diagnostic !
void Interpolator3D2OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    // probes are interpolated one by one for now
    
    int idx[3], idxO[3];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, ipart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, ipart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;
    idx[2]  = round( particles.position( 2, ipart ) * D_inv[2] );
    idxO[2] = idx[2] - k_domain_begin -1 ;
    
    Field3D *Env_A_abs_3D  = static_cast<Field3D *>( EMfields->Env_A_abs_ );
    Field3D *Env_Chi_3D    = static_cast<Field3D *>( EMfields->Env_Chi_ );
    Field3D *Env_E_abs_3D  = static_cast<Field3D *>( EMfields->Env_E_abs_ );
    Field3D *Env_Ex_abs_3D = static_cast<Field3D *>( EMfields->Env_Ex_abs_ );
    
    double coeff[3][2][3];
    
    double delta0, delta;
    double delta2;
    
    
    for( int i=0; i<3; i++ ) { // for X/Y
        delta0 = particles.position( i, ipart )*D_inv[i];
        
        delta   = delta0 - ( double )idx[i] ;
        delta2  = delta*delta;
        
        coeff[i][0][0]    =  0.5 * ( delta2-delta+0.25 );
        coeff[i][0][1]    = ( 0.75 - delta2 );
        coeff[i][0][2]    =  0.5 * ( delta2+delta+0.25 );
    }
    
    
    double *coeffyp = &( coeff[1][0][1] );
    double *coeffxp = &( coeff[0][0][1] );
    double *coeffzp = &( coeff[2][0][1] );
    
    // Interpolation of Env_A_abs^(p,p,p) (absolute value of envelope A)
    double interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) * ( *Env_A_abs_3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
            }
        }
    }
    *Env_A_abs_Loc= interp_res;
    
    // Interpolation of Env_Chi^(p,p,p)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) * ( *Env_Chi_3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
            }
        }
    }
    *Env_Chi_Loc= interp_res;
    
    // Interpolation of Env_E_abs^(p,p,p) (absolute value of envelope E transverse)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) * ( *Env_E_abs_3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
            }
        }
    }
    *Env_E_abs_Loc= interp_res;

    // Interpolation of Env_Ex_abs^(p,p,p) (absolute value of envelope Ex)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            for( int kloc=-1 ; kloc<2 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) * ( *Env_Ex_abs_3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc );
            }
        }
    }
    *Env_Ex_abs_Loc= interp_res;
    
}


// Interpolator on another field than the basic ones
void Interpolator3D2OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    ERROR( "Single field 3D2O interpolator not available in vectorized mode" );
}
