#include "Interpolator3DWT4OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3DWT4OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3DWT4OrderV::Interpolator3DWT4OrderV( Params &params, Patch *patch ) : Interpolator3DWT4Order( params, patch )
{
    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    d_inv_[2] = 1.0/params.cell_length[2];
    dt_ov_D[0] = params.timestep/params.cell_length[0];
    dt_ov_D[1] = params.timestep/params.cell_length[1];
    dt_ov_D[2] = params.timestep/params.cell_length[2];
    dt2_ov_D2[0] = dt_ov_D[0] * dt_ov_D[0];
    dt2_ov_D2[1] = dt_ov_D[1] * dt_ov_D[1];
    dt2_ov_D2[2] = dt_ov_D[2] * dt_ov_D[2];
    D_ov_96dt[0] = 1.0/96.0/dt_ov_D[0];
    D_ov_96dt[1] = 1.0/96.0/dt_ov_D[1];
    D_ov_96dt[2] = 1.0/96.0/dt_ov_D[2];

    //double defined for use in coefficients
    dble_1_ov_6 = 1.0/6.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_115_ov_192 = 115.0/192.0;
}

// ---------------------------------------------------------------------------------------------------------------------
// 4th OrderV WT Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
// void Interpolator3DWT4OrderV::fields( ElectroMagn *EMfields, Particles &particles, int ipart, double *ELoc, double *BLoc )
// {
// }

void Interpolator3DWT4OrderV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }

    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );

    double * __restrict__ Epart[3];
    double * __restrict__ Bpart[3];

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
    idx[0]  = round( particles.position( 0, *istart ) * d_inv_[0] );
    idxO[0] = idx[0] - i_domain_begin  ;
    idx[1]  = round( particles.position( 1, *istart ) * d_inv_[1] );
    idxO[1] = idx[1] - j_domain_begin  ;
    idx[2]  = round( particles.position( 2, *istart ) * d_inv_[2] );
    idxO[2] = idx[2] - k_domain_begin  ;

    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_m );

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ position_z = particles.getPtrPosition(2);

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
            double var1, var2, var3, var4, var5;

            // i = 0

            delta0 = position_x[ipart+ivect+istart[0]] *d_inv_[0];
            dual [0][ipart] = ( delta0 - ( double )idx[0] >=0. );

            // j = 0

            delta   = delta0 - ( double )idx[0] ;
            var1 = dt_ov_D[0] - delta;
            var1 = var1 * var1;
            var1 = D_ov_96dt[0] * var1 * var1;
            var2 = dt_ov_D[0] + delta;
            var2 = var2 * var2;
            var2 = D_ov_96dt[0] * var2 * var2;
            var3 = copysign( var1, dt_ov_D[0]-delta ) + copysign( var2, dt_ov_D[0]+delta);
            var4 =  dble_1_ov_24 * ((( 3.0+delta ) * delta  - ( 3.0-dt2_ov_D2[0] )) * delta + (1.0+dt2_ov_D2[0] ));
            var5 =  dble_1_ov_24 * ((( 3.0-delta ) * delta  + ( 3.0-dt2_ov_D2[0] )) * delta + (1.0+dt2_ov_D2[0] ));

            coeff[0][0][0][ipart] = var3 + var1 - var2;
            coeff[0][0][1][ipart] = 4.0 * ( var4-var3 );
            coeff[0][0][2][ipart] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
            coeff[0][0][3][ipart] = 4.0 * ( var5-var3 );
            coeff[0][0][4][ipart] = var3 + var2 - var1;

            deltaO[0][ipart-ipart_ref+ivect+istart[0]] = delta;

            // j = 1

            delta   = delta0 - ( double )idx[0] + ( 0.5-dual[0][ipart] );
            var1 = dble_11_ov_24 * delta;
            var2 = delta * delta;
            var3 = dble_1_ov_6 * var2;
            var4 = 0.5 - delta;
            var4 = var4 * var4;
            var5 = 0.5 + delta;
            var5 = var5 * var5;

            coeff[0][1][0][ipart] = dble_1_ov_24 * var4 * var4;
            coeff[0][1][1][ipart] = dble_19_ov_96 - var1 + var3 * ( 1.5+delta-var2 );
            coeff[0][1][2][ipart] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
            coeff[0][1][3][ipart] = dble_19_ov_96 + var1 + var3 * ( 1.5-delta -var2 );
            coeff[0][1][4][ipart] = dble_1_ov_24 * var5 * var5;

            // i = 1

            delta0 = position_y[ipart+ivect+istart[0]] *d_inv_[1];
            dual [1][ipart] = ( delta0 - ( double )idx[1] >=0. );

            // j = 0

            delta   = delta0 - ( double )idx[1];
            var1 = dt_ov_D[1] - delta;
            var1 = var1 * var1;
            var1 = D_ov_96dt[1] * var1 * var1;
            var2 = dt_ov_D[1] + delta;
            var2 = var2 * var2;
            var2 = D_ov_96dt[1] * var2 * var2;
            var3 = copysign( var1, dt_ov_D[1]-delta ) + copysign( var2, dt_ov_D[1]+delta);
            var4 =  dble_1_ov_24 * ((( 3.0+delta ) * delta  - ( 3.0-dt2_ov_D2[1] )) * delta + (1.0+dt2_ov_D2[1] ));
            var5 =  dble_1_ov_24 * ((( 3.0-delta ) * delta  + ( 3.0-dt2_ov_D2[1] )) * delta + (1.0+dt2_ov_D2[1] ));

            coeff[1][0][0][ipart] = var3 + var1 - var2;
            coeff[1][0][1][ipart] = 4.0 * ( var4-var3 );
            coeff[1][0][2][ipart] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
            coeff[1][0][3][ipart] = 4.0 * ( var5-var3 );
            coeff[1][0][4][ipart] = var3 + var2 - var1;

            deltaO[1][ipart-ipart_ref+ivect+istart[0]] = delta;

            // j = 1

            delta   = delta0 - ( double )idx[1] + ( 0.5-dual[1][ipart] );
            var1 = dble_11_ov_24 * delta;
            var2 = delta * delta;
            var3 = dble_1_ov_6 * var2;
            var4 = 0.5 - delta;
            var4 = var4 * var4;
            var5 = 0.5 + delta;
            var5 = var5 * var5;

            coeff[1][1][0][ipart] = dble_1_ov_24 * var4 * var4;
            coeff[1][1][1][ipart] = dble_19_ov_96 - var1 + var3 * ( 1.5+delta-var2 );
            coeff[1][1][2][ipart] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
            coeff[1][1][3][ipart] = dble_19_ov_96 + var1 + var3 * ( 1.5-delta -var2 );
            coeff[1][1][4][ipart] = dble_1_ov_24 * var5 * var5;

            // i = 2

            delta0 = position_z[ipart+ivect+istart[0]] *d_inv_[2];
            dual [2][ipart] = ( delta0 - ( double )idx[2] >=0. );

            // j = 0

            delta   = delta0 - ( double )idx[2];
            var1 = dt_ov_D[2] - delta;
            var1 = var1 * var1;
            var1 = D_ov_96dt[2] * var1 * var1;
            var2 = dt_ov_D[2] + delta;
            var2 = var2 * var2;
            var2 = D_ov_96dt[2] * var2 * var2;
            var3 = copysign( var1, dt_ov_D[2]-delta ) + copysign( var2, dt_ov_D[2]+delta);
            var4 =  dble_1_ov_24 * ((( 3.0+delta ) * delta  - ( 3.0-dt2_ov_D2[2] )) * delta + (1.0+dt2_ov_D2[2] ));
            var5 =  dble_1_ov_24 * ((( 3.0-delta ) * delta  + ( 3.0-dt2_ov_D2[2] )) * delta + (1.0+dt2_ov_D2[2] ));

            coeff[2][0][0][ipart] = var3 + var1 - var2;
            coeff[2][0][1][ipart] = 4.0 * ( var4-var3 );
            coeff[2][0][2][ipart] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
            coeff[2][0][3][ipart] = 4.0 * ( var5-var3 );
            coeff[2][0][4][ipart] = var3 + var2 - var1;

            deltaO[2][ipart-ipart_ref+ivect+istart[0]] = delta;

            // j = 1

            delta   = delta0 - ( double )idx[2] + 1*( 0.5-dual[2][ipart] );
            var1 = dble_11_ov_24 * delta;
            var2 = delta * delta;
            var3 = dble_1_ov_6 * var2;
            var4 = 0.5 - delta;
            var4 = var4 * var4;
            var5 = 0.5 + delta;
            var5 = var5 * var5;

            coeff[2][1][0][ipart] = dble_1_ov_24 * var4 * var4;
            coeff[2][1][1][ipart] = dble_19_ov_96 - var1 + var3 * ( 1.5+delta-var2 );
            coeff[2][1][2][ipart] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
            coeff[2][1][3][ipart] = dble_19_ov_96 + var1 + var3 * ( 1.5-delta -var2 );
            coeff[2][1][4][ipart] = dble_1_ov_24 * var5 * var5;

            // --------------------------------------
            // The code above is the hand-unrollling of this loop
            //
            // for( int i=0; i<3; i++ ) { // for X/Y
            //     delta0 = particles.position( i, ipart+ivect+istart[0] )*d_inv_[i];
            //     dual [i][ipart] = ( delta0 - ( double )idx[i] >=0. );
            //
            //     for( int j=0; j<2; j++ ) { // for dual
            //
            //         delta   = delta0 - ( double )idx[i] + ( double )j*( 0.5-dual[i][ipart] );
            //         delta2  = delta*delta;
            //         delta3  = delta2*delta;
            //         delta4  = delta3*delta;
            //
            //         coeff[i][j][0][ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            //         coeff[i][j][1][ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            //         coeff[i][j][2][ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            //         coeff[i][j][3][ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            //         coeff[i][j][4][ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            //
            //         if( j==0 ) {
            //             deltaO[i][ipart-ipart_ref+ivect+istart[0]] = delta;
            //         }
            //     }
            // }
        }

        double interp_res = 0;

        #if defined __INTEL_COMPILER

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {

                double *coeffypt = &( coeff[1][0][2][ipart] );
                double *coeffyd  = &( coeff[1][1][2][ipart] );
                double *coeffxd  = &( coeff[0][1][2][ipart] );
                double *coeffxpt = &( coeff[0][0][2][ipart] );
                double *coeffzpt = &( coeff[2][0][2][ipart] );
                double *coeffzd  = &( coeff[2][1][2][ipart] );

                //Ex(dual, primal WT, primal WT)
                interp_res = 0.;
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += *( coeffxd+iloc*32 ) * *( coeffypt+jloc*32 ) * *( coeffzpt+kloc*32 ) *
                                          ( ( 1-dual[0][ipart] )*( *Ex3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0][ipart]*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
                        }
                    }
                }
                Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;


                //Ey(primal WT, dual, primal WT)
                interp_res = 0.;
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += *( coeffxpt+iloc*32 ) * *( coeffyd+jloc*32 ) * *( coeffzpt+kloc*32 ) *
                                          ( ( 1-dual[1][ipart] )*( *Ey3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1][ipart]*( *Ey3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
                        }
                    }
                }
                Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;


                //Ez(primal WT, primal WT, dual)
                interp_res = 0.;
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += *( coeffxpt+iloc*32 ) * *( coeffypt+jloc*32 ) * *( coeffzd+kloc*32 ) *
                                          ( ( 1-dual[2][ipart] )*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2][ipart]*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
                        }
                    }
                }
                Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

            interp_Bx( idxO, np_computed, &(coeff[0][0][2][0]), &(coeff[1][1][2][0]), &(coeff[2][1][2][0]), &(dual[1][0]), &(dual[2][0]), Bx3D, &(Bpart[0][ivect+istart[0]-ipart_ref]) );
            interp_By( idxO, np_computed, &(coeff[0][1][2][0]), &(coeff[1][0][2][0]), &(coeff[2][1][2][0]), &(dual[0][0]), &(dual[2][0]), By3D, &(Bpart[1][ivect+istart[0]-ipart_ref]) );
            interp_Bz( idxO, np_computed, &(coeff[0][1][2][0]), &(coeff[1][1][2][0]), &(coeff[2][0][2][0]), &(dual[0][0]), &(dual[1][0]), Bz3D, &(Bpart[2][ivect+istart[0]-ipart_ref]) );


        #else

            double * __restrict__ coeffxpt2 = &( coeff[0][0][2][0] );
            double * __restrict__ coeffxd2  = &( coeff[0][1][2][0] );

            double * __restrict__ coeffypt2 = &( coeff[1][0][2][0] );
            double * __restrict__ coeffyd2  = &( coeff[1][1][2][0] );

            double * __restrict__ coeffzpt2 = &( coeff[2][0][2][0] );
            double * __restrict__ coeffzd2  = &( coeff[2][1][2][0] );

            double field_buffer[6][6][6];

            // ---------------------------------------------------------------------
            //Ex(dual, primal WT, primal WT)

            // Field buffers for vectorization (required on A64FX)
            for( int iloc=-2 ; iloc<4 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                        field_buffer[iloc+2][jloc+2][kloc+2] = Ex3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
                    }
                }
            }

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {

                interp_res = 0.;

                UNROLL_S(5)
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    UNROLL_S(5)
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        UNROLL_S(5)
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += coeffxd2[ipart+iloc*32] * coeffypt2[ipart+jloc*32] * coeffzpt2[ipart+kloc*32] *
                                ( ( 1-dual[0][ipart] )* field_buffer[iloc+2][jloc+2][kloc+2]
                                + dual[0][ipart]*field_buffer[iloc+3][jloc+2][kloc+2] );
                        }
                    }
                }
                Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

            // ---------------------------------------------------------------------
            //Ey(primal WT, dual, primal WT)

            // Field buffers for vectorization (required on A64FX)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<4 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                        field_buffer[iloc+2][jloc+2][kloc+2] = Ey3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
                    }
                }
            }

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {

                interp_res = 0.;

                UNROLL_S(5)
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    UNROLL_S(5)
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        UNROLL_S(5)
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += coeffxpt2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] * coeffzpt2[ipart+kloc*32] *
                                ( ( 1-dual[1][ipart] )* field_buffer[iloc+2][jloc+2][kloc+2]
                                + dual[1][ipart]*field_buffer[2+iloc][jloc+3][kloc+2] );
                        }
                    }
                }
                Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

            // ---------------------------------------------------------------------
            //Ez(primal WT, primal WT, dual)

            // Field buffers for vectorization (required on A64FX)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<4 ; kloc++ ) {
                        field_buffer[iloc+2][jloc+2][kloc+2] = Ez3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
                    }
                }
            }

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {

                interp_res = 0.;

                UNROLL_S(5)
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    UNROLL_S(5)
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        UNROLL_S(5)
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += coeffxpt2[ipart+iloc*32] * coeffypt2[ipart+jloc*32] * coeffzd2[ipart+kloc*32] *
                                ( ( 1-dual[2][ipart] )* field_buffer[iloc+2][jloc+2][kloc+2] +
                                dual[2][ipart]*field_buffer[2+iloc][jloc+2][kloc+3] );
                        }
                    }
                }
                Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

            // ---------------------------------------------------------------------
            //Bx(primal WT, dual , dual )

            // interp_Bx( idxO, np_computed, &(coeff[0][0][2][0]), &(coeff[1][1][2][0]), &(coeff[2][1][2][0]),
            //             &(dual[1][0]), &(dual[2][0]), Bx3D, &(Bpart[0][ivect+istart[0]-ipart_ref]) );

            // Field buffers for vectorization (required on A64FX)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<4 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<4 ; kloc++ ) {
                        field_buffer[iloc+2][jloc+2][kloc+2] = Bx3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
                    }
                }
            }

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                interp_res = 0.;
                UNROLL_S(5)
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    UNROLL_S(5)
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        UNROLL_S(5)
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += coeffxpt2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] * coeffzd2[ipart+kloc*32] *
                                ( ( 1-dual[2][ipart] ) * ( ( 1-dual[1][ipart] )*field_buffer[2+iloc][2+jloc][2+kloc] + dual[1][ipart]*field_buffer[2+iloc][3+jloc][2+kloc] )
                                  +    dual[2][ipart]  * ( ( 1-dual[1][ipart] )*field_buffer[2+iloc][2+jloc][3+kloc] + dual[1][ipart]*field_buffer[2+iloc][3+jloc][3+kloc] ) );
                        }
                    }
                }
                Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

            // ---------------------------------------------------------------------
            //By(dual, primal WT, dual )

            // interp_By( idxO, np_computed, &(coeff[0][1][2][0]), &(coeff[1][0][2][0]), &(coeff[2][1][2][0]),
            //             &(dual[0][0]), &(dual[2][0]), By3D, &(Bpart[1][ivect+istart[0]-ipart_ref]) );

            // Field buffers for vectorization (required on A64FX)
            for( int iloc=-2 ; iloc<4 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<4 ; kloc++ ) {
                        field_buffer[iloc+2][jloc+2][kloc+2] = By3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
                    }
                }
            }

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {

                interp_res = 0.;
                UNROLL_S(5)
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    UNROLL_S(5)
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        UNROLL_S(5)
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += coeffxd2[ipart+iloc*32] * coeffypt2[ipart+jloc*32] * coeffzd2[ipart+kloc*32] *
                                ( ( 1-dual[2][ipart] ) * ( ( 1-dual[0][ipart] )*field_buffer[2+iloc][2+jloc][2+kloc] + dual[0][ipart]*field_buffer[3+iloc][2+jloc][2+kloc] )
                                  +    dual[2][ipart]  * ( ( 1-dual[0][ipart] )*field_buffer[2+iloc][2+jloc][3+kloc] + dual[0][ipart]*field_buffer[3+iloc][2+jloc][3+kloc] ) );
                        }
                    }
                }
                Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

            // ---------------------------------------------------------------------
            //Bz(dual, dual, primal WT )

            // interp_Bz( idxO, np_computed, &(coeff[0][1][2][0]), &(coeff[1][1][2][0]), &(coeff[2][0][2][0]),
            //             &(dual[0][0]), &(dual[1][0]), Bz3D, &(Bpart[2][ivect+istart[0]-ipart_ref]) );

            // Field buffers for vectorization (required on A64FX)
            for( int iloc=-2 ; iloc<4 ; iloc++ ) {
                for( int jloc=-2 ; jloc<4 ; jloc++ ) {
                    for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                        field_buffer[iloc+2][jloc+2][kloc+2] = Bz3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
                    }
                }
            }

            #pragma omp simd
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {

                interp_res = 0.;
                UNROLL_S(5)
                for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                    UNROLL_S(5)
                    for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        UNROLL_S(5)
                        for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                            interp_res += coeffxd2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] * coeffzpt2[ipart+kloc*32] *
                                ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*field_buffer[2+iloc][2+jloc][2+kloc] + dual[0][ipart]*field_buffer[3+iloc][2+jloc][2+kloc] )
                                  +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*field_buffer[2+iloc][3+jloc][2+kloc] + dual[0][ipart]*field_buffer[3+iloc][3+jloc][2+kloc] ) );
                        }
                    }
                }
                Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
            }

         #endif


        // #pragma omp simd
        // for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        //
        //     interp_res = 0.;
        //     UNROLL_S(5)
        //     for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        //         UNROLL_S(5)
        //         for( int jloc=-2 ; jloc<3 ; jloc++ ) {
        //             UNROLL_S(5)
        //             for( int kloc=-2 ; kloc<3 ; kloc++ ) {
        //                 interp_res += coeffxp2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] * coeffzd2[ipart+kloc*32] *
        //                     ( ( 1-dual[2][ipart] ) * ( ( 1-dual[1][ipart] )*Bx3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc]
        //                      + dual[1][ipart]*Bx3D->data_3D[idxO[0]+iloc][idxO[1]+jloc+1][idxO[2]+kloc] )
        //                       +    dual[2][ipart]  * ( ( 1-dual[1][ipart] )*Bx3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc+1] +
        //                        dual[1][ipart]*Bx3D->data_3D[idxO[0]+iloc][idxO[1]+jloc+1][idxO[2]+kloc+1] ) );
        //             }
        //         }
        //     }
        //     Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        //
        //     interp_res = 0.;
        //     UNROLL_S(5)
        //     for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        //         UNROLL_S(5)
        //         for( int jloc=-2 ; jloc<3 ; jloc++ ) {
        //             UNROLL_S(5)
        //             for( int kloc=-2 ; kloc<3 ; kloc++ ) {
        //                 interp_res += coeffxd2[ipart+iloc*32] * coeffyp2[ipart+jloc*32] * coeffzd2[ipart+kloc*32] *
        //                     ( ( 1-dual[2][ipart] ) * ( ( 1-dual[0][ipart] )*By3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc]
        //                     + dual[0][ipart]*By3D->data_3D[idxO[0]+iloc+1][idxO[1]+jloc][idxO[2]+kloc] )
        //                       +    dual[2][ipart]  * ( ( 1-dual[0][ipart] )*By3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc+1]
        //                       + dual[0][ipart]*By3D->data_3D[idxO[0]+iloc+1][idxO[1]+jloc][idxO[2]+kloc+1] ) );
        //             }
        //         }
        //     }
        //     Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        //
        //     interp_res = 0.;
        //     UNROLL_S(5)
        //     for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        //         UNROLL_S(5)
        //         for( int jloc=-2 ; jloc<3 ; jloc++ ) {
        //             UNROLL_S(5)
        //             for( int kloc=-2 ; kloc<3 ; kloc++ ) {
        //                 interp_res += coeffxd2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] * coeffzp2[ipart+kloc*32] *
        //                     ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*Bz3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc]
        //                     + dual[0][ipart]*Bz3D->data_3D[idxO[0]+iloc+1][idxO[1]+jloc][idxO[2]+kloc] )
        //                       +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*Bz3D->data_3D[idxO[0]+iloc][idxO[1]+jloc+1][idxO[2]+kloc]
        //                       + dual[0][ipart]*Bz3D->data_3D[idxO[0]+iloc+1][idxO[1]+jloc+1][idxO[2]+kloc] ) );
        //             }
        //         }
        //     }
        //     Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        // }
    }
} // END Interpolator3DWT4OrderV


void Interpolator3DWT4OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    idx[0]  = round( particles.position( 0, *istart ) * d_inv_[0] );
    idxO[0] = idx[0] - i_domain_begin  ;
    idx[1]  = round( particles.position( 1, *istart ) * d_inv_[1] );
    idxO[1] = idx[1] - j_domain_begin  ;
    idx[2]  = round( particles.position( 2, *istart ) * d_inv_[2] );
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

    double coeff[3][3][5];
    int dual[3]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).

    double delta0, delta;
    double var1, var2, var3, var4, var5;

    for( int i=0; i<3; i++ ) { // for X/Y
        delta0 = particles.position( i, ipart )*d_inv_[i];
        dual [i] = ( delta0 - ( double )idx[i] >=0. );

        // For primal
        delta   = delta0 - ( double )idx[i];
        var1 = dble_11_ov_24 * delta;
        var2 = delta * delta;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - delta;
        var4 = var4 * var4;
        var5 = 0.5 + delta;
        var5 = var5 * var5;
        coeff[i][0][0] = dble_1_ov_24 * var4 * var4;
        coeff[i][0][1] = dble_19_ov_96 - var1 + var3 * ( 1.5+delta-var2 );
        coeff[i][0][2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeff[i][0][3] = dble_19_ov_96 + var1 + var3 * ( 1.5-delta -var2 );
        coeff[i][0][4] = dble_1_ov_24 * var5 * var5;

        // For primal WT
        var1 = dt_ov_D[i] - delta;
        var1 = var1 * var1;
        var1 = D_ov_96dt[i] * var1 * var1;
        var2 = dt_ov_D[i] + delta;
        var2 = var2 * var2;
        var2 = D_ov_96dt[i] * var2 * var2;
        var3 = copysign( var1, dt_ov_D[i]-delta ) + copysign( var2, dt_ov_D[i]+delta);
        var4 =  dble_1_ov_24 * ((( 3.0+delta ) * delta  - ( 3.0-dt2_ov_D2[i] )) * delta + (1.0+dt2_ov_D2[i] ));
        var5 =  dble_1_ov_24 * ((( 3.0-delta ) * delta  + ( 3.0-dt2_ov_D2[i] )) * delta + (1.0+dt2_ov_D2[i] ));
        coeff[i][1][0] = var3 + var1 - var2;
        coeff[i][1][1] = 4.0 * ( var4-var3 );
        coeff[i][1][2] = 1.0 + 6.0 * var3 - 4.0 * ( var4+var5 );
        coeff[i][1][3] = 4.0 * ( var5-var3 );
        coeff[i][1][4] = var3 + var2 - var1;

        // For dual
        delta   = delta0 - ( double )idx[i] + ( 0.5-dual[i] );
        var1 = dble_11_ov_24 * delta;
        var2 = delta * delta;
        var3 = dble_1_ov_6 * var2;
        var4 = 0.5 - delta;
        var4 = var4 * var4;
        var5 = 0.5 + delta;
        var5 = var5 * var5;
        coeff[i][2][0] = dble_1_ov_24 * var4 * var4;
        coeff[i][2][1] = dble_19_ov_96 - var1 + var3 * ( 1.5+delta-var2 );
        coeff[i][2][2] = dble_115_ov_192 + 0.25 * var2 * ( -2.5+var2 );
        coeff[i][2][3] = dble_19_ov_96 + var1 + var3 * ( 1.5-delta -var2 );
        coeff[i][2][4] = dble_1_ov_24 * var5 * var5;
            
    }

    double *coeffyp  = &( coeff[1][0][2] );
    double *coeffypt = &( coeff[1][1][2] );
    double *coeffyd  = &( coeff[1][2][2] );
    double *coeffxd  = &( coeff[0][2][2] );
    double *coeffxpt = &( coeff[0][1][2] );
    double *coeffxp  = &( coeff[0][0][2] );
    double *coeffzp  = &( coeff[2][0][2] );
    double *coeffzpt = &( coeff[2][1][2] );
    double *coeffzd  = &( coeff[2][2][2] );
    
    //Ex(dual, primal WT, primal WT)
    double interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffypt+jloc*1 ) * *( coeffzpt+kloc*1 ) *
                              ( ( 1-dual[0] )*( *Ex3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *Ex3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) );
            }
        }
    }
    Epart[0][ipart] = interp_res;


    //Ey(primal WT, dual, primal WT)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxpt+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzpt+kloc*1 ) *
                              ( ( 1-dual[1] )*( *Ey3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Ey3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) );
            }
        }
    }
    Epart[1][ipart] = interp_res;

    //Ez(primal WT, primal WT, dual)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxpt+iloc*1 ) * *( coeffypt+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] )*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[2]*( *Ez3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) );
            }
        }
    }
    Epart[2][ipart] = interp_res;


    //Bx(primal WT, dual , dual )
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxpt+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] ) * ( ( 1-dual[1] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[1]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) )
                                +    dual[2]  * ( ( 1-dual[1] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[1]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) ) );
            }
        }
    }
    Bpart[0][ipart] = interp_res;

    //By(dual, primal WT, dual )
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffypt+jloc*1 ) * *( coeffzd+kloc*1 ) *
                              ( ( 1-dual[2] ) * ( ( 1-dual[0] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dual[0]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                                +    dual[2]  * ( ( 1-dual[0] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dual[0]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) ) );
            }
        }
    }
    Bpart[1][ipart] = interp_res;

    //Bz(dual, dual, primal WT)
    interp_res = 0.;
    for( int iloc=-2 ; iloc<3 ; iloc++ ) {
        for( int jloc=-2 ; jloc<3 ; jloc++ ) {
            for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                interp_res += *( coeffxd+iloc*1 ) * *( coeffyd+jloc*1 ) * *( coeffzpt+kloc*1 ) *
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



void Interpolator3DWT4OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator3DWT4OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator3DWT4OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}



// Interpolator on another field than the basic ones
void Interpolator3DWT4OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    ERROR( "Single field 3D4O interpolator not available in vectorized mode" );
}

void Interpolator3DWT4OrderV::interp_Bx( int* idxO, int np_computed, double *coeffxpt, double *coeffyd, double *coeffzd, int *dualy, int* dualz, Field3D *Bx3D, double *Bpart ) {

    #pragma omp simd
    for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        double *coeff_yd = &( coeffyd[ipart] );
        double *coeff_xpt = &( coeffxpt[ipart] );
        double *coeff_zd = &( coeffzd[ipart] );
        //Bx(primal, dual , dual )
        double interp_res = 0.;
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeff_xpt+iloc*32 ) * *( coeff_yd+jloc*32 ) * *( coeff_zd+kloc*32 ) *
                        ( ( 1-dualz[ipart] ) * ( ( 1-dualy[ipart] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dualy[ipart]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) )
                          +    dualz[ipart]  * ( ( 1-dualy[ipart] )*( *Bx3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dualy[ipart]*( *Bx3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+1+kloc ) ) );
                }
            }
        }
        Bpart[ipart] = interp_res;
    }
}

void Interpolator3DWT4OrderV::interp_By( int* idxO, int np_computed, double *coeffxd, double *coeffypt, double *coeffzd, int *dualx, int* dualz, Field3D *By3D, double *Bpart )
{
    #pragma omp simd
    for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        double *coeff_ypt = &( coeffypt[ipart] );
        double *coeff_xd = &( coeffxd[ipart] );
        double *coeff_zd = &( coeffzd[ipart] );
        //By(dual, primal, dual )
        double interp_res = 0.;
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeff_xd+iloc*32 ) * *( coeff_ypt+jloc*32 ) * *( coeff_zd+kloc*32 ) *
                        ( ( 1-dualz[ipart] ) * ( ( 1-dualx[ipart] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dualx[ipart]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                          +    dualz[ipart]  * ( ( 1-dualx[ipart] )*( *By3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) + dualx[ipart]*( *By3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+1+kloc ) ) );
                }
            }
        }
        Bpart[ipart] = interp_res;
    }
}

void Interpolator3DWT4OrderV::interp_Bz( int* idxO, int np_computed, double *coeffxd, double *coeffyd, double *coeffzpt, int *dualx, int* dualy, Field3D *Bz3D, double *Bpart )
{
    #pragma omp simd
    for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        double *coeff_yd = &( coeffyd[ipart] );
        double *coeff_xd = &( coeffxd[ipart] );
        double *coeff_zpt = &( coeffzpt[ipart] );
        //Bz(dual, dual, prim )
        double interp_res = 0.;
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                for( int kloc=-2 ; kloc<3 ; kloc++ ) {
                    interp_res += *( coeff_xd+iloc*32 ) * *( coeff_yd+jloc*32 ) * *( coeff_zpt+kloc*32 ) *
                        ( ( 1-dualy[ipart] ) * ( ( 1-dualx[ipart] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+jloc, idxO[2]+kloc ) + dualx[ipart]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+jloc, idxO[2]+kloc ) )
                          +    dualy[ipart]  * ( ( 1-dualx[ipart] )*( *Bz3D )( idxO[0]+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) + dualx[ipart]*( *Bz3D )( idxO[0]+1+iloc, idxO[1]+1+jloc, idxO[2]+kloc ) ) );
                }
            }
        }
        Bpart[ipart] = interp_res;
    }
}

// void Interpolator3D4OrderV::interp_Bx( int * __restrict__ idxO, int np_computed,
//                                         double * __restrict__ coeffxp,
//                                         double * __restrict__ coeffyd,
//                                         double * __restrict__ coeffzd,
//                                         int * __restrict__ dualy,
//                                         int * __restrict__ dualz,
//                                         Field3D *Bx3D,
//                                         double * __restrict__ Bpart ) {
//
//     double field_buffer[6][6][6];
//
//     // Field buffers for vectorization (required on A64FX)
//     for( int iloc=-2 ; iloc<3 ; iloc++ ) {
//         for( int jloc=-2 ; jloc<4 ; jloc++ ) {
//             for( int kloc=-2 ; kloc<4 ; kloc++ ) {
//                 field_buffer[iloc+2][jloc+2][kloc+2] = Bx3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
//             }
//         }
//     }
//
//     #pragma omp simd
//     for( int ipart=0 ; ipart<np_computed; ipart++ ) {
//         //Bx(primal, dual , dual )
//         double interp_res = 0.;
//         UNROLL_S(5)
//         for( int iloc=-2 ; iloc<3 ; iloc++ ) {
//             UNROLL_S(5)
//             for( int jloc=-2 ; jloc<3 ; jloc++ ) {
//                 UNROLL_S(5)
//                 for( int kloc=-2 ; kloc<3 ; kloc++ ) {
//                     interp_res += coeffxp[ipart+iloc*32] * coeffyd[ipart+jloc*32] * coeffzd[ipart+kloc*32] *
//                         ( ( 1-dualz[ipart] ) * ( ( 1-dualy[ipart] )*field_buffer[2+iloc][2+jloc][2+kloc] + dualy[ipart]*field_buffer[2+iloc][3+jloc][2+kloc] )
//                           +    dualz[ipart]  * ( ( 1-dualy[ipart] )*field_buffer[2+iloc][2+jloc][3+kloc] + dualy[ipart]*field_buffer[2+iloc][3+jloc][3+kloc] ) );
//                 }
//             }
//         }
//         Bpart[ipart] = interp_res;
//     }
// }
//
// void Interpolator3D4OrderV::interp_By( int * __restrict__  idxO,
//                                         int np_computed,
//                                         double * __restrict__  coeffxd,
//                                         double * __restrict__  coeffyp,
//                                         double * __restrict__  coeffzd,
//                                         int * __restrict__  dualx,
//                                         int * __restrict__  dualz,
//                                         Field3D * __restrict__ By3D,
//                                         double * __restrict__ Bpart )
// {
//
//     double field_buffer[6][6][6];
//
//     // Field buffers for vectorization (required on A64FX)
//     for( int iloc=-2 ; iloc<4 ; iloc++ ) {
//         for( int jloc=-2 ; jloc<3 ; jloc++ ) {
//             for( int kloc=-2 ; kloc<4 ; kloc++ ) {
//                 field_buffer[iloc+2][jloc+2][kloc+2] = By3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
//             }
//         }
//     }
//
//     #pragma omp simd
//     for( int ipart=0 ; ipart<np_computed; ipart++ ) {
//         //By(dual, primal, dual )
//         double interp_res = 0.;
//         UNROLL_S(5)
//         for( int iloc=-2 ; iloc<3 ; iloc++ ) {
//             UNROLL_S(5)
//             for( int jloc=-2 ; jloc<3 ; jloc++ ) {
//                 UNROLL_S(5)
//                 for( int kloc=-2 ; kloc<3 ; kloc++ ) {
//                     interp_res += coeffxd[ipart+iloc*32] * coeffyp[ipart+jloc*32] * coeffzd[ipart+kloc*32] *
//                         ( ( 1-dualz[ipart] ) * ( ( 1-dualx[ipart] )*field_buffer[2+iloc][2+jloc][2+kloc] + dualx[ipart]*field_buffer[3+iloc][2+jloc][2+kloc] )
//                           +    dualz[ipart]  * ( ( 1-dualx[ipart] )*field_buffer[2+iloc][2+jloc][3+kloc] + dualx[ipart]*field_buffer[3+iloc][2+jloc][3+kloc] ) );
//                 }
//             }
//         }
//         Bpart[ipart] = interp_res;
//     }
// }
//
// void Interpolator3D4OrderV::interp_Bz( int * __restrict__ idxO,
//                                         int np_computed,
//                                         double * __restrict__ coeffxd,
//                                         double * __restrict__ coeffyd,
//                                         double * __restrict__ coeffzp,
//                                         int * __restrict__ dualx,
//                                         int * __restrict__ dualy,
//                                         Field3D * __restrict__ Bz3D,
//                                         double * __restrict__ Bpart )
// {
//
//     double field_buffer[6][6][6];
//
//     // Field buffers for vectorization (required on A64FX)
//     for( int iloc=-2 ; iloc<4 ; iloc++ ) {
//         for( int jloc=-2 ; jloc<4 ; jloc++ ) {
//             for( int kloc=-2 ; kloc<3 ; kloc++ ) {
//                 field_buffer[iloc+2][jloc+2][kloc+2] = Bz3D->data_3D[idxO[0]+iloc][idxO[1]+jloc][idxO[2]+kloc];
//             }
//         }
//     }
//
//     #pragma omp simd
//     for( int ipart=0 ; ipart<np_computed; ipart++ ) {
//
//         //Bz(dual, dual, prim )
//         double interp_res = 0.;
//         UNROLL_S(5)
//         for( int iloc=-2 ; iloc<3 ; iloc++ ) {
//             UNROLL_S(5)
//             for( int jloc=-2 ; jloc<3 ; jloc++ ) {
//                 UNROLL_S(5)
//                 for( int kloc=-2 ; kloc<3 ; kloc++ ) {
//                     interp_res += coeffxd[ipart+iloc*32] * coeffyd[ipart+jloc*32] * coeffzp[ipart+kloc*32] *
//                         ( ( 1-dualy[ipart] ) * ( ( 1-dualx[ipart] )*field_buffer[2+iloc][2+jloc][2+kloc] + dualx[ipart]*field_buffer[3+iloc][2+jloc][2+kloc] )
//                           +    dualy[ipart]  * ( ( 1-dualx[ipart] )*field_buffer[2+iloc][3+jloc][2+kloc] + dualx[ipart]*field_buffer[3+iloc][3+jloc][2+kloc] ) );
//                 }
//             }
//         }
//         Bpart[ipart] = interp_res;
//     }
// }
