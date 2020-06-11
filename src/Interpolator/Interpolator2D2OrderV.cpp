#include "Interpolator2D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D2OrderV::Interpolator2D2OrderV( Params &params, Patch *patch ) : Interpolator2D( params, patch )
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    D_inv[0] = 1.0/params.cell_length[0];
    D_inv[1] = 1.0/params.cell_length[1];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd OrderV Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D2OrderV::fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc )
{
}

void Interpolator2D2OrderV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }
    
    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );
    
    double *Epart[3], *Bpart[3];
    
    double *deltaO[2];
    deltaO[0] = &( smpi->dynamics_deltaold[ithread][0] );
    deltaO[1] = &( smpi->dynamics_deltaold[ithread][nparts] );
    
    for( unsigned int k=0; k<3; k++ ) {
        Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts] );
        Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts] );
    }
    
    int idx[2], idxO[2];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;
    
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );
    
    double coeff[2][2][3][32];
    int dual[2][32]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
    
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
            
            for( int i=0; i<2; i++ ) { // for X/Y
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
            
            //Ex(dual, primal)
            double interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) *
                                  ( ( 1-dual[0][ipart] )*( *Ex2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0][ipart]*( *Ex2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) );
                }
            }
            Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            //Ey(primal, dual)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) *
                                  ( ( 1-dual[1][ipart] )*( *Ey2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[1][ipart]*( *Ey2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) );
                }
            }
            Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            
            //Ez(primal, primal)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * ( *Ez2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
                }
            }
            Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            //Bx(primal, dual)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) *
                                  ( ( ( 1-dual[1][ipart] )*( *Bx2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[1][ipart]*( *Bx2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) ) );
                }
            }
            Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            //By(dual, primal )
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) *
                                  ( ( ( 1-dual[0][ipart] )*( *By2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0][ipart]*( *By2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) ) );
                }
            }
            Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            //Bz(dual, dual)
            interp_res = 0.;
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += *( coeffxd+iloc*32 ) * *( coeffyd+jloc*32 ) *
                                  ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0][ipart]*( *Bz2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) )
                                    +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) + dual[0][ipart]*( *Bz2D )( idxO[0]+2+iloc, idxO[1]+2+jloc ) ) );
                }
            }
            Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
        }
    }
    
} // END Interpolator2D2OrderV

void Interpolator2D2OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    
    int idx[2], idxO[2];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;
    
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );
    
    double coeff[2][2][3];
    int dual[2]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
    
    double delta0, delta;
    double delta2;
    
    for( int i=0; i<2; i++ ) { // for X/Y
        delta0 = particles.position( i, ipart )*D_inv[i];
        dual [i] = ( delta0 - ( double )idx[i] >=0. );
        
        for( int j=0; j<2; j++ ) { // for dual
        
            delta   = delta0 - ( double )idx[i] + ( double )j*( 0.5-dual[i] );
            delta2  = delta*delta;
            
            coeff[i][j][0]    =  0.5 * ( delta2-delta+0.25 );
            coeff[i][j][1]    = ( 0.75 - delta2 );
            coeff[i][j][2]    =  0.5 * ( delta2+delta+0.25 );
            
        }
    }
    
    
    double *coeffyp = &( coeff[1][0][1] );
    double *coeffyd = &( coeff[1][1][1] );
    double *coeffxd = &( coeff[0][1][1] );
    double *coeffxp = &( coeff[0][0][1] );
    
    //Ex(dual, primal)
    double interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) *
                          ( ( 1-dual[0] )*( *Ex2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0]*( *Ex2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) );
        }
    }
    Epart[0][ipart] = interp_res;
    
    //Ey(primal, dual)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) *
                          ( ( 1-dual[1] )*( *Ey2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[1]*( *Ey2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) );
        }
    }
    Epart[1][ipart] = interp_res;
    
    
    //Ez(primal, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * ( *Ez2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
        }
    }
    Epart[2][ipart] = interp_res;
    
    //Bx(primal, dual)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) *
                          ( ( ( 1-dual[1] )*( *Bx2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[1]*( *Bx2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) ) );
        }
    }
    Bpart[0][ipart] = interp_res;
    
    //By(dual, primal )
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) *
                          ( ( ( 1-dual[0] )*( *By2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0]*( *By2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) ) );
        }
    }
    Bpart[1][ipart] = interp_res;
    
    //Bz(dual, dual)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxd+iloc*1 ) * *( coeffyd+jloc*1 ) *
                          ( ( 1-dual[1] ) * ( ( 1-dual[0] )*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0]*( *Bz2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) )
                            +    dual[1]  * ( ( 1-dual[0] )*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) + dual[0]*( *Bz2D )( idxO[0]+2+iloc, idxO[1]+2+jloc ) ) );
        }
    }
    Bpart[2][ipart] = interp_res;
    
    Field2D *Jx2D = static_cast<Field2D *>( EMfields->Jx_ );
    Field2D *Jy2D = static_cast<Field2D *>( EMfields->Jy_ );
    Field2D *Jz2D = static_cast<Field2D *>( EMfields->Jz_ );
    Field2D *rho2D = static_cast<Field2D *>( EMfields->rho_ );
    
    //Jx(dual, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) *
                          ( ( 1-dual[0] )*( *Jx2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[0]*( *Jx2D )( idxO[0]+2+iloc, idxO[1]+1+jloc ) );
        }
    }
    JLoc->x = interp_res;
    
    //Jy(primal, dual)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) *
                          ( ( 1-dual[1] )*( *Jy2D )( idxO[0]+1+iloc, idxO[1]+1+jloc ) + dual[1]*( *Jy2D )( idxO[0]+1+iloc, idxO[1]+2+jloc ) );
        }
    }
    JLoc->y = interp_res;
    
    
    //Jz(primal, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * ( *Jz2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
        }
    }
    JLoc->z = interp_res;
    
    //Rho(primal, primal)
    interp_res = 0.;
    for( int iloc=-1 ; iloc<2 ; iloc++ ) {
        for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * ( *rho2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
        }
    }
    ( *RhoLoc ) = interp_res;
    
}


// Interpolator on another field than the basic ones
void Interpolator2D2OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    ERROR( "Single field 2D2O interpolator not available in vectorized mode" );
}

void Interpolator2D2OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Vectorized interpolation for the envelope model is not implemented for 2D geometry" );
} // END Interpolator2D2Order


void Interpolator2D2OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Vectorized interpolation for the envelope model is not implemented for 2D geometry" );
} // END Interpolator2D2Order


void Interpolator2D2OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Vectorized interpolation for the envelope model is not implemented for 2D geometry" );
} // END Interpolator2D2Order


