#include "Interpolator2D4OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D4OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D4OrderV::Interpolator2D4OrderV( Params &params, Patch *patch ) : Interpolator2D( params, patch )
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];


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
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D4OrderV::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
{
    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*dx_inv_;
    double ypn = particles.position( 1, ipart )*dy_inv_;
    // Calculate coeffs
    coeffs( xpn, ypn );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], Ex2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp_[2], &coeffyp_[2], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], By2D, id_, jp_ );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], Bz2D, id_, jd_ );
} // END Interpolator2D4OrderV

void Interpolator2D4OrderV::fieldsWrapper(  ElectroMagn *EMfields,
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

            // Normalized particle position
            double xpn = particles.position( 0, ipart )*dx_inv_;
            double ypn = particles.position( 1, ipart )*dy_inv_;

            // Declaration and calculation of the coefficient for interpolation
            double delta0, delta, delta2, delta3, delta4;

            // i = 0

            delta0 = position_x[ipart+ivect+istart[0]] *D_inv[0];
            dual [0][ipart] = ( delta0 - ( double )idx[0] >=0. );

            // j = 0

            delta   = delta0 - ( double )idx[0] ;
            delta2  = delta*delta;
            delta3  = delta2*delta;
            delta4  = delta3*delta;

            coeff[0][0][0][ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeff[0][0][1][ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[0][0][2][ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeff[0][0][3][ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[0][0][4][ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            deltaO[0][ipart-ipart_ref+ivect+istart[0]] = delta;

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


            delta   = delta0 - ( double )idx[0] + ( 0.5-dual[0][ipart] );
            delta2  = delta*delta;
            delta3  = delta2*delta;
            delta4  = delta3*delta;

            coeff[0][1][0][ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeff[0][1][1][ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[0][1][2][ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeff[0][1][3][ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[0][1][4][ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            // i = 1

            delta0 = position_y[ipart+ivect+istart[0]] *D_inv[1];
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

            delta   = delta0 - ( double )idx[1] + 0*( 0.5-dual[1][ipart] );
            delta2  = delta*delta;
            delta3  = delta2*delta;
            delta4  = delta3*delta;

            coeff[1][0][0][ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeff[1][0][1][ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[1][0][2][ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeff[1][0][3][ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[1][0][4][ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            deltaO[1][ipart-ipart_ref+ivect+istart[0]] = delta;

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


            delta   = delta0 - ( double )idx[1] + 1*( 0.5-dual[1][ipart] );
            delta2  = delta*delta;
            delta3  = delta2*delta;
            delta4  = delta3*delta;

            coeff[1][1][0][ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            coeff[1][1][1][ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[1][1][2][ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4 * delta4;
            coeff[1][1][3][ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4 * delta2  - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            coeff[1][1][4][ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;


            //!\todo CHECK if this is correct for both primal & dual grids !!!
            // First index for summation
            // ip_ = ip_ - i_domain_begin;
            // id_ = id_ - i_domain_begin;
            // jp_ = jp_ - j_domain_begin;
            // jd_ = jd_ - j_domain_begin;

        }

        double * __restrict__ coeffxp2 = &( coeff[0][0][2][0] );
        double * __restrict__ coeffxd2 = &( coeff[0][1][2][0] );

        double * __restrict__ coeffyp2 = &( coeff[1][0][2][0] );
        double * __restrict__ coeffyd2 = &( coeff[1][1][2][0] );

        double * __restrict__ coeffzp2 = &( coeff[2][0][2][0] );
        double * __restrict__ coeffzd2 = &( coeff[2][1][2][0] );

        // Local buffer to store the field components
        double field_buffer[6][6];

        // Interpolation result
        double interp_res = 0;

        // ---------------------------------------------------------------------
        //Ex(dual, primal)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-2 ; iloc<4 ; iloc++ ) {
            for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *Ex2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            //Ex(dual, primal)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyp2[ipart+jloc*32] *
                                  ( ( 1-dual[0][ipart] )*field_buffer[2+iloc][2+jloc]
                                  + dual[0][ipart]*field_buffer[3+iloc][2+jloc] );
                }
            }
            Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Ey(primal, dual)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        interp_res += coeffxp2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                                      ( ( 1-dual[1][ipart] )*( *Ey2D )( idxO[0]+iloc, idxO[1]+jloc)
                                      + dual[1][ipart]*( *Ey2D )( idxO[0]+iloc, idxO[1]+1+jloc ) );
                }
            }
            Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Ez(primal, primal)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxp2[ipart+iloc*32] * coeffyp2[ipart+jloc*32]  *
                                    ( *Ez2D )( idxO[0]+iloc, idxO[1]+jloc );
                }
            }
            Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Bx(primal, dual)
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                        interp_res += coeffxp2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                            (  ( 1-dual[1][ipart] )*( *Bx2D )( idxO[0]+iloc, idxO[1]+jloc )
                            + dual[1][ipart]*( *Bx2D )( idxO[0]+iloc, idxO[1]+1+jloc ) );
                }
            }
            Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //By(dual, primal )
            interp_res = 0.;
            UNROLL_S(5)
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                UNROLL_S(5)
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyp2[ipart+jloc*32] *
                             ( ( 1-dual[0][ipart] )*( *By2D )( idxO[0]+iloc, idxO[1]+jloc )
                            + dual[0][ipart]*( *By2D )( idxO[0]+1+iloc, idxO[1]+jloc) );
                }
            }
            Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;

            //Bz(dual, dual )
            interp_res = 0.;
            for( int iloc=-2 ; iloc<3 ; iloc++ ) {
                for( int jloc=-2 ; jloc<3 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                            ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*( *Bz2D )( idxO[0]+iloc, idxO[1]+jloc) + dual[0][ipart]*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+jloc) )
                              +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*( *Bz2D )( idxO[0]+iloc, idxO[1]+1+jloc) + dual[0][ipart]*( *Bz2D )( idxO[0]+1+iloc, idxO[1]+1+jloc) ) );
                }
            }
            Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;

        }

    }

}

void Interpolator2D4OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
{

    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    // Interpolate E, B
    // Compute coefficient for ipart position
    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );
    Field2D *Jx2D = static_cast<Field2D *>( EMfields->Jx_ );
    Field2D *Jy2D = static_cast<Field2D *>( EMfields->Jy_ );
    Field2D *Jz2D = static_cast<Field2D *>( EMfields->Jz_ );
    Field2D *Rho2D= static_cast<Field2D *>( EMfields->rho_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*dx_inv_;
    double ypn = particles.position( 1, ipart )*dy_inv_;
    // Calculate coeffs
    coeffs( xpn, ypn );

    int nparts( particles.size() );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) =  compute( &coeffxd_[2], &coeffyp_[2], Ex2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp_[2], &coeffyp_[2], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], By2D, id_, jp_ );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], Bz2D, id_, jd_ );
    // Interpolation of Jx^(d,p)
    JLoc->x = compute( &coeffxd_[2], &coeffyp_[2], Jx2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    JLoc->y = compute( &coeffxp_[2], &coeffyd_[2], Jy2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    JLoc->z = compute( &coeffxp_[2], &coeffyp_[2], Jz2D, ip_, jp_ );
    // Interpolation of Rho^(p,p)
    ( *RhoLoc ) = compute( &coeffxp_[2], &coeffyp_[2], Rho2D, ip_, jp_ );
}

// Interpolator on another field than the basic ones
void Interpolator2D4OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    ERROR( "Single field 2D2O interpolator not available in vectorized mode" );
}

void Interpolator2D4OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator2D4OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator2D4OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
