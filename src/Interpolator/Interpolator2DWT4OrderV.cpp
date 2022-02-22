#include "Interpolator2DWT4OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2DWT4OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2DWT4OrderV::Interpolator2DWT4OrderV( Params &params, Patch *patch ) : Interpolator2DWT4Order( params, patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    dt_ov_D[0] = params.timestep/params.cell_length[0];
    dt_ov_D[1] = params.timestep/params.cell_length[1];
    dt2_ov_D2[0] = dt_ov_D[0] * dt_ov_D[0];
    dt2_ov_D2[1] = dt_ov_D[1] * dt_ov_D[1];
    D_ov_96dt[0] = 1.0/96.0/dt_ov_D[0];
    D_ov_96dt[1] = 1.0/96.0/dt_ov_D[1];

    //double defined for use in coefficients
    dble_1_ov_6 = 1.0/6.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_115_ov_192 = 115.0/192.0;
}


void Interpolator2DWT4OrderV::fieldsWrapper(  ElectroMagn *EMfields,
                                              Particles &particles,
                                              SmileiMPI *smpi,
                                              int *istart, int *iend, int ithread, int ipart_ref )
{
    // std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    // std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    // std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    // std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    //
    // //Loop on bin particles
    // int nparts( particles.size() );
    // for( int ipart=*istart ; ipart<*iend; ipart++ ) {
    //     //Interpolation on current particle
    //     fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );
    //     //Buffering of iol and delta
    //     ( *iold )[ipart+0*nparts]  = ip_;
    //     ( *iold )[ipart+1*nparts]  = jp_;
    //     ( *delta )[ipart+0*nparts] = deltax;
    //     ( *delta )[ipart+1*nparts] = deltay;
    // }


    // ---------------------------------------------

    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }

    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );

    double * __restrict__ Epart[3];
    double * __restrict__ Bpart[3];

    double *deltaO[2];
    deltaO[0] = &( smpi->dynamics_deltaold[ithread][0] );
    deltaO[1] = &( smpi->dynamics_deltaold[ithread][nparts] );

    for( unsigned int k=0; k<3; k++ ) {
        Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts] );
        Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts] );
    }

    int idx[2], idxO[2];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * d_inv_[0] );
    idxO[0] = idx[0] - i_domain_begin  ;
    idx[1]  = round( particles.position( 1, *istart ) * d_inv_[1] );
    idxO[1] = idx[1] - j_domain_begin  ;

    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);

    double coeff[2][2][5][32];

    // Size ndim. Boolean indicating if the part has a dual indice
    // equal to the primal one (dual=0) or if it is +1 (dual=1).
    int dual[2][32];

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

            // Declaration and calculation of the coefficient for interpolation
            double delta0, delta, var1, var2, var3, var4, var5;

            // i = 0

            delta0 = position_x[ipart+ivect+istart[0]] *d_inv_[0];
            dual [0][ipart] = ( delta0 - ( double )idx[0] >=0. );

            // coefficient x primal grid (i=0, j=0)

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

            // std::cerr << "ipart: " << ipart +ivect+istart[0]
            //           << " x: " << particles.position( 0, ipart+ivect+istart[0] )
            //           << " y: " << particles.position( 1, ipart+ivect+istart[0] )
            //           << std::endl;
            //
            // std::cerr << " delta: " << delta
            //           << " coeffxp: " << coeff[0][0][0][ipart]
            //           << " " << coeff[0][0][1][ipart]
            //           << " " << coeff[0][0][2][ipart]
            //           << " " << coeff[0][0][3][ipart]
            //           << " " << coeff[0][0][4][ipart]
            //           << std::endl;

            // j = 1

            // deltax   = delta0 - ( double )idx[0] + 0.5;
            // delta2  = deltax*deltax;
            // delta3  = delta2*deltax;
            // delta4  = delta3*deltax;
            // coeffxd_[0] = dble_1_ov_384   - dble_1_ov_48  * deltax  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            // coeffxd_[1] = dble_19_ov_96   - dble_11_ov_24 * deltax  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            // coeffxd_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            // coeffxd_[3] = dble_19_ov_96   + dble_11_ov_24 * deltax  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            // coeffxd_[4] = dble_1_ov_384   + dble_1_ov_48  * deltax  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;


            // coefficient x dual grid (i=0, j=1)

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

            // std::cerr
            //           << " delta: " << delta
            //           << " coeffxd: " << coeff[0][1][0][ipart]
            //           << " " << coeff[0][1][1][ipart]
            //           << " " << coeff[0][1][2][ipart]
            //           << " " << coeff[0][1][3][ipart]
            //           << " " << coeff[0][1][4][ipart]
            //           << std::endl;


            // i = 1

            delta0 = position_y[ipart+ivect+istart[0]] *d_inv_[1];
            dual [1][ipart] = ( delta0 - ( double )idx[1] >=0. );

            // j = 0

            // deltay   = ypn - ( double )jp_;
            // delta2  = deltay*deltay;
            // delta3  = delta2*deltay;
            // delta4  = delta3*deltay;
            // coeffyp_[0] = dble_1_ov_384   - dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            // coeffyp_[1] = dble_19_ov_96   - dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            // coeffyp_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            // coeffyp_[3] = dble_19_ov_96   + dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            // coeffyp_[4] = dble_1_ov_384   + dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            delta   = delta0 - ( double )idx[1] ;
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

            // std::cerr
            //           << " delta: " << delta
            //           << " coeffyp: " << coeff[1][0][0][ipart]
            //           << " " << coeff[1][0][1][ipart]
            //           << " " << coeff[1][0][2][ipart]
            //           << " " << coeff[1][0][3][ipart]
            //           << " " << coeff[1][0][4][ipart]
            //           << std::endl;

            // j = 1

            // deltay   = ypn - ( double )jd_ + 0.5;
            // delta2  = deltay*deltay;
            // delta3  = delta2*deltay;
            // delta4  = delta3*deltay;
            // coeffyd_[0] = dble_1_ov_384   - dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            // coeffyd_[1] = dble_19_ov_96   - dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            // coeffyd_[2] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            // coeffyd_[3] = dble_19_ov_96   + dble_11_ov_24 * deltay  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            // coeffyd_[4] = dble_1_ov_384   + dble_1_ov_48  * deltay  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;


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

            // std::cerr
            //           << " delta: " << delta
            //           << " coeffyd: " << coeff[1][1][0][ipart]
            //           << " " << coeff[1][1][1][ipart]
            //           << " " << coeff[1][1][2][ipart]
            //           << " " << coeff[1][1][3][ipart]
            //           << " " << coeff[1][1][4][ipart]
            //           << std::endl;
            //
            // std::cerr << " xpn: " << particles.position( 0, *istart ) * d_inv_[0]
            //           << " ypn: " << particles.position( 1, *istart ) * d_inv_[1]
            //           << " idx0[0]: " << idxO[0]
            //           << " idx0[1]: " << idxO[1] << std::endl;

            //!\todo CHECK if this is correct for both primal & dual grids !!!
            // First index for summation
            // ip_ = ip_ - i_domain_begin;
            // id_ = id_ - i_domain_begin;
            // jp_ = jp_ - j_domain_begin;
            // jd_ = jd_ - j_domain_begin;

        }

        double * __restrict__ coeffxpt2 = &( coeff[0][0][2][0] );
        double * __restrict__ coeffxd2  = &( coeff[0][1][2][0] );

        double * __restrict__ coeffypt2 = &( coeff[1][0][2][0] );
        double * __restrict__ coeffyd2  = &( coeff[1][1][2][0] );

        // Local buffer to store the field components
        double field_buffer[6][6];

        // Interpolation result
        double interp_res = 0;

        // ---------------------------------------------------------------------
        //Ex(dual, primal WT)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-2 ; iloc<4 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                field_buffer[iloc+2][jloc+2] = ( *Ex2D )( idxO[0]+iloc, idxO[1]+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            //Ex(dual, primal WT)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffypt2[ipart+jloc*32] *
                                  ( ( 1-dual[0][ipart] )*field_buffer[iloc+2][jloc+2]
                                  + dual[0][ipart]*field_buffer[iloc+3][jloc+2] );
                }
            }
            Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        }

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-2 ; iloc<3 ; iloc++ ) {
            for( int jloc=-2 ; jloc<4 ; jloc++ ) {
                field_buffer[iloc+2][jloc+2] = ( *Ey2D )( idxO[0]+iloc, idxO[1]+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            //Ey(primal WT, dual)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        interp_res += coeffxpt2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                                      ( ( 1-dual[1][ipart] )*( *Ey2D )( idxO[0]+iloc, idxO[1]+jloc)
                                      + dual[1][ipart]*( *Ey2D )( idxO[0]+iloc, idxO[1]+1+jloc ) );
                }
            }
            Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Ez(primal WT, primal WT)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxpt2[ipart+iloc*32] * coeffypt2[ipart+jloc*32]  *
                                    ( *Ez2D )( idxO[0]+iloc, idxO[1]+jloc );
                }
            }
            Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Bx(primal WT, dual)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        interp_res += coeffxpt2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                            (  ( 1-dual[1][ipart] )*( *Bx2D )( idxO[0]+iloc, idxO[1]+jloc )
                            + dual[1][ipart]*( *Bx2D )( idxO[0]+iloc, idxO[1]+1+jloc ) );
                }
            }
            Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //By(dual, primal WT)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffypt2[ipart+jloc*32] *
                             ( ( 1-dual[0][ipart] )*( *By2D )( idxO[0]+iloc, idxO[1]+jloc )
                            + dual[0][ipart]*( *By2D )( idxO[0]+1+iloc, idxO[1]+jloc) );
                }
            }
            Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Bz(dual, dual )
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                            ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*( *Bz2D )( idxO[0]+iloc, idxO[1]+jloc) + dual[0][ipart]*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+jloc) )
                              +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*( *Bz2D )( idxO[0]+iloc, idxO[1]+1+jloc) + dual[0][ipart]*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+1+jloc) ) );
                }
            }
            Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;

        }

        // for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        //     std::cerr
        //               << " Ex: " << Epart[0][ipart-ipart_ref+ivect+istart[0]]
        //               << " Ey: " << Epart[1][ipart-ipart_ref+ivect+istart[0]]
        //               << " Ez: " << Epart[2][ipart-ipart_ref+ivect+istart[0]]
        //               << " Bx: " << Bpart[0][ipart-ipart_ref+ivect+istart[0]]
        //               << " By: " << Bpart[1][ipart-ipart_ref+ivect+istart[0]]
        //               << " Bz: " << Bpart[2][ipart-ipart_ref+ivect+istart[0]]
        //               << std::endl;
        // }

    }

}

// Interpolator on another field than the basic ones
void Interpolator2DWT4OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    ERROR( "Single field 2D2O interpolator not available in vectorized mode" );
}

void Interpolator2DWT4OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator2DWT4OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator2DWT4OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
