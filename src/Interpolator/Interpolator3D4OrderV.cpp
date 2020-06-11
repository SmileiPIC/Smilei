#include "Interpolator3D4OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D4OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D4OrderV::Interpolator3D4OrderV( Params &params, Patch *patch ) : Interpolator3D( params, patch )
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    dz_inv_ = 1.0/params.cell_length[2];
    D_inv[0] = 1.0/params.cell_length[0];
    D_inv[1] = 1.0/params.cell_length[1];
    D_inv[2] = 1.0/params.cell_length[2];
    
    //double defined for use in coefficients
    dble_1_ov_384 = 1.0/384.0;
    dble_1_ov_48 = 1.0/48.0;
    dble_1_ov_16 = 1.0/16.0;
    dble_1_ov_12 = 1.0/12.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_1_ov_4 = 1.0/4.0;
    dble_1_ov_6 = 1.0/6.0;
    dble_115_ov_192 = 115.0/192.0;
    dble_5_ov_8 = 5.0/8.0;
    
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd OrderV Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D4OrderV::fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc )
{
}

void Interpolator3D4OrderV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }
    
    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );
    
    double *Epart[3], *Bpart[3];
    
    double *deltaO[3];
    deltaO[0] = &( smpi->dynamics_deltaold[ithread][0] );
    deltaO[1] = &( smpi->dynamics_deltaold[ithread][nparts] );
    deltaO[2] = &( smpi->dynamics_deltaold[ithread][2*nparts] );
    
    for( unsigned int k=0; k<3; k++ ) {
        Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts] );
        Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts] );
    }
    
    int idx[3], idxO[3];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin  ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin  ;
    idx[2]  = round( particles.position( 2, *istart ) * D_inv[2] );
    idxO[2] = idx[2] - k_domain_begin  ;
    
    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_m );
    
    double coeff[3][2][5][32];
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
            double delta2, delta3, delta4;
            
            
            for( int i=0; i<3; i++ ) { // for X/Y
                delta0 = particles.position( i, ipart+ivect+istart[0] )*D_inv[i];
                dual [i][ipart] = ( delta0 - ( double )idx[i] >=0. );
                
                for( int j=0; j<2; j++ ) { // for dual
                
                    delta   = delta0 - ( double )idx[i] + ( double )j*( 0.5-dual[i][ipart] );
                    delta2  = delta*delta;
                    delta3  = delta2*delta;
                    delta4  = delta3*delta;
                    
                    coeff[i][j][0][ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
                    coeff[i][j][1][ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
                    coeff[i][j][2][ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
                    coeff[i][j][3][ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
                    coeff[i][j][4][ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
                    
                    if( j==0 ) {
                        deltaO[i][ipart-ipart_ref+ivect+istart[0]] = delta;
                    }
                }
            }
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double *coeffyp = &( coeff[1][0][2][ipart] );
            double *coeffyd = &( coeff[1][1][2][ipart] );
            double *coeffxd = &( coeff[0][1][2][ipart] );
            double *coeffxp = &( coeff[0][0][2][ipart] );
            double *coeffzp = &( coeff[2][0][2][ipart] );
            double *coeffzd = &( coeff[2][1][2][ipart] );
            
            //Ex(dual, primal, primal)
            double interp_res = 0.;
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                        interp_res += *( coeffxd+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[0][ipart] )*( *Ex3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0][ipart]*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
                    }
                }
            }
            Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            
            //Ey(primal, dual, primal)
            interp_res = 0.;
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzp+kloc*32 ) *
                                      ( ( 1-dual[1][ipart] )*( *Ey3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1][ipart]*( *Ey3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
                    }
                }
            }
            Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            
            
            //Ez(primal, primal, dual)
            interp_res = 0.;
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                        interp_res += *( coeffxp+iloc*32 ) * *( coeffyp+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                      ( ( 1-dual[2][ipart] )*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2][ipart]*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
                    }
                }
            }
            Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;

        }

        interp_Bx( idxO, np_computed, &(coeff[0][0][2][0]), &(coeff[1][1][2][0]), &(coeff[2][1][2][0]), &(dual[1][0]), &(dual[2][0]), Bx3D, &(Bpart[0][ivect+istart[0]-ipart_ref]) );
        interp_By( idxO, np_computed, &(coeff[0][1][2][0]), &(coeff[1][0][2][0]), &(coeff[2][1][2][0]), &(dual[0][0]), &(dual[2][0]), By3D, &(Bpart[1][ivect+istart[0]-ipart_ref]) );
        interp_Bz( idxO, np_computed, &(coeff[0][1][2][0]), &(coeff[1][1][2][0]), &(coeff[2][0][2][0]), &(dual[0][0]), &(dual[1][0]), Bz3D, &(Bpart[2][ivect+istart[0]-ipart_ref]) );
        
    }
} // END Interpolator3D4OrderV


void Interpolator3D4OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    
    int idx[3], idxO[3];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin  ;
    idx[1]  = round( particles.position( 1, *istart ) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin  ;
    idx[2]  = round( particles.position( 2, *istart ) * D_inv[2] );
    idxO[2] = idx[2] - k_domain_begin  ;
    
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
    
    double coeff[3][2][5];
    int dual[3]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
    
    double delta0, delta;
    double delta2, delta3, delta4;
    
    for( int i=0; i<3; i++ ) { // for X/Y
        delta0 = particles.position( i, ipart )*D_inv[i];
        dual [i] = ( delta0 - ( double )idx[i] >=0. );
        
        for( int j=0; j<2; j++ ) { // for dual
        
            delta   = delta0 - ( double )idx[i] + ( double )j*( 0.5-dual[i] );
            delta2  = delta*delta;
            delta3  = delta2*delta;
            delta4  = delta3*delta;
            
            coeff[i][j][0] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeff[i][j][1] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[i][j][2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeff[i][j][3] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[i][j][4] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            
        }
    }
    
    double *coeffyp = &( coeff[1][0][2] );
    double *coeffyd = &( coeff[1][1][2] );
    double *coeffxd = &( coeff[0][1][2] );
    double *coeffxp = &( coeff[0][0][2] );
    double *coeffzp = &( coeff[2][0][2] );
    double *coeffzd = &( coeff[2][1][2] );
    
    //Ex(dual, primal, primal)
    double interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[0] )*( *Ex3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
            }
        }
    }
    Epart[0][ipart] = interp_res;
    
    
    //Ey(primal, dual, primal)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[1] )*( *Ey3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Ey3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
            }
        }
    }
    Epart[1][ipart] = interp_res;
    
    //Ez(primal, primal, dual)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] )*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2]*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
            }
        }
    }
    Epart[2][ipart] = interp_res;
    
    
    //Bx(primal, dual , dual )
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] ) * ( ( 1-dual[1] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) )
                                +    dual[2]  * ( ( 1-dual[1] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[1]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) ) );
            }
        }
    }
    Bpart[0][ipart] = interp_res;
    
    //By(dual, primal, dual )
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] ) * ( ( 1-dual[0] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                +    dual[2]  * ( ( 1-dual[0] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[0]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) ) );
            }
        }
    }
    Bpart[1][ipart] = interp_res;
    
    //Bz(dual, dual, prim )
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[1] ) * ( ( 1-dual[0] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                +    dual[1]  * ( ( 1-dual[0] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) + dual[0]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) ) );
            }
        }
    }
    Bpart[2][ipart] = interp_res;
    
    
    //Jx(dual, primal, primal)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[0] )*( *Jx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Jx3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
            }
        }
    }
    JLoc->x = interp_res;
    
    
    //Jy(primal, dual, primal)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzp+kloc*1 ) *
                              ( ( 1-dual[1] )*( *Jy3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Jy3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
            }
        }
    }
    JLoc->y = interp_res;
    
    
    //Jz(primal, primal, dual)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] )*( *Jz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2]*( *Jz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
            }
        }
    }
    JLoc->z = interp_res;
    
    //Rho(primal, primal, primal)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxp+iloc*1 ) * *( coeffyp+jloc*1 ) * *( coeffzp+kloc*1 ) * ( *rho3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc );
            }
        }
    }
    ( *RhoLoc ) = interp_res;
    
}



void Interpolator3D4OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator3D4OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator3D4OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}



// Interpolator on another field than the basic ones
void Interpolator3D4OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    ERROR( "Single field 3D4O interpolator not available in vectorized mode" );
}

void Interpolator3D4OrderV::interp_Bx( int* idxO, int np_computed, double *coeffxp, double *coeffyd, double *coeffzd, int *dualy, int* dualz, Field3D *Bx3D, double *Bpart ) {

    #pragma omp simd
    for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        double *coeff_yd = &( coeffyd[ipart] );
        double *coeff_xp = &( coeffxp[ipart] );
        double *coeff_zd = &( coeffzd[ipart] );
        //Bx(primal, dual , dual )
        double interp_res = 0.;
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeff_xp+iloc*32 ) * *( coeff_yd+jloc*32 ) * *( coeff_zd+kloc*32 ) *
                        ( ( 1-dualz[ipart] ) * ( ( 1-dualy[ipart] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dualy[ipart]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) )
                          +    dualz[ipart]  * ( ( 1-dualy[ipart] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dualy[ipart]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) ) );
                }
            }
        }
        Bpart[ipart] = interp_res;
    }
}

void Interpolator3D4OrderV::interp_By( int* idxO, int np_computed, double *coeffxd, double *coeffyp, double *coeffzd, int *dualx, int* dualz, Field3D *By3D, double *Bpart )
{
    #pragma omp simd
    for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        double *coeff_yp = &( coeffyp[ipart] );
        double *coeff_xd = &( coeffxd[ipart] );
        double *coeff_zd = &( coeffzd[ipart] );
        //By(dual, primal, dual )
        double interp_res = 0.;
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeff_xd+iloc*32 ) * *( coeff_yp+jloc*32 ) * *( coeff_zd+kloc*32 ) *
                        ( ( 1-dualz[ipart] ) * ( ( 1-dualx[ipart] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dualx[ipart]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                          +    dualz[ipart]  * ( ( 1-dualx[ipart] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dualx[ipart]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) ) );
                }
            }
        }
        Bpart[ipart] = interp_res;
    }
}

void Interpolator3D4OrderV::interp_Bz( int* idxO, int np_computed, double *coeffxd, double *coeffyd, double *coeffzp, int *dualx, int* dualy, Field3D *Bz3D, double *Bpart )
{
    #pragma omp simd
    for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        double *coeff_yd = &( coeffyd[ipart] );
        double *coeff_xd = &( coeffxd[ipart] );
        double *coeff_zp = &( coeffzp[ipart] );
        //Bz(dual, dual, prim )
        double interp_res = 0.;
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeff_xd+iloc*32 ) * *( coeff_yd+jloc*32 ) * *( coeff_zp+kloc*32 ) *
                        ( ( 1-dualy[ipart] ) * ( ( 1-dualx[ipart] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dualx[ipart]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                          +    dualy[ipart]  * ( ( 1-dualx[ipart] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) + dualx[ipart]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) ) );
                }
            }
        }
        Bpart[ipart] = interp_res;
    }
}
