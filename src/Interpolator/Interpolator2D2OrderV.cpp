#include "Interpolator2D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
//! Creator for Interpolator2D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D2OrderV::Interpolator2D2OrderV( Params &params, Patch *patch ) : Interpolator2D2Order( params, patch )
{
    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];

}

// -----------------------------------------------------------------------------
//! Wrapper called by the particle dynamics section
// -----------------------------------------------------------------------------
void Interpolator2D2OrderV::fieldsWrapper(  ElectroMagn *EMfields,
                                            Particles &particles,
                                            SmileiMPI *smpi,
                                            int *istart,
                                            int *iend,
                                            int ithread,
                                            unsigned int,
                                            int ipart_ref )
{
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
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, *istart ) * d_inv_[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;

    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);

    double coeff[2][2][3][32];

     // Size ndim. Boolean indicating if the part
     // has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
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

            double delta0, delta;
            double delta2;

            // for( int i=0; i<2; i++ ) { // for X/Y
            //     delta0 = particles.position( i, ipart+ivect+istart[0] )*d_inv_[i];
            //     dual [i][ipart] = ( delta0 - ( double )idx[i] >=0. );
            //
            //     for( int j=0; j<2; j++ ) { // for dual
            //
            //         delta   = delta0 - ( double )idx[i] + ( double )j*( 0.5-dual[i][ipart] );
            //         delta2  = delta*delta;
            //
            //         coeff[i][j][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            //         coeff[i][j][1][ipart]    = ( 0.75 - delta2 );
            //         coeff[i][j][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            //
            //         if( j==0 ) {
            //             deltaO[i][ipart-ipart_ref+ivect+istart[0]] = delta;
            //         }
            //
            //     }
            // }

            // i = 0

            delta0 = position_x[ipart+ivect+istart[0]] *d_inv_[0];
            dual [0][ipart] = ( delta0 - ( double )idx[0] >=0. );

            // j = 0

            delta   = delta0 - ( double )idx[0] ;
            delta2  = delta*delta;
            coeff[0][0][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[0][0][1][ipart]    = ( 0.75 - delta2 );
            coeff[0][0][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            deltaO[0][ipart-ipart_ref+ivect+istart[0]] = delta;

            // j = 1

            delta   = delta0 - ( double )idx[0] + ( 0.5-dual[0][ipart] );
            delta2  = delta*delta;
            coeff[0][1][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[0][1][1][ipart]    = ( 0.75 - delta2 );
            coeff[0][1][2][ipart]    =  0.5 * ( delta2+delta+0.25 );

            // i = 1

            delta0 = position_y[ipart+ivect+istart[0]] *d_inv_[1];
            dual [1][ipart] = ( delta0 - ( double )idx[1] >=0. );

            // j = 0

            delta   = delta0 - ( double )idx[1] + ( double )0*( 0.5-dual[1][ipart] );
            delta2  = delta*delta;
            coeff[1][0][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[1][0][1][ipart]    = ( 0.75 - delta2 );
            coeff[1][0][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            deltaO[1][ipart-ipart_ref+ivect+istart[0]] = delta;

            // j = 1

            delta   = delta0 - ( double )idx[1] + ( double )1*( 0.5-dual[1][ipart] );
            delta2  = delta*delta;
            coeff[1][1][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[1][1][1][ipart]    = ( 0.75 - delta2 );
            coeff[1][1][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
        }

        // Coefficient pointer on primal and dual nodes
        double * __restrict__ coeffyp2 = &( coeff[1][0][1][0] );
        double * __restrict__ coeffyd2 = &( coeff[1][1][1][0] );
        double * __restrict__ coeffxd2 = &( coeff[0][1][1][0] );
        double * __restrict__ coeffxp2 = &( coeff[0][0][1][0] );

        // Local buffer to store the field components
        double field_buffer[4][4];

        double interp_res = 0.;

        //Ex(dual, primal)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-1 ; iloc<3 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *Ex2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            //Ex(dual, primal)
            interp_res = 0.;
            UNROLL_S(3)
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                UNROLL_S(3)
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyp2[ipart+jloc*32] *
                                  ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][1+jloc]
                                  + dual[0][ipart]*field_buffer[2+iloc][1+jloc] );
                }
            }
            Epart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        }

        //Ey(primal, dual)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<3 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *Ey2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            interp_res = 0.;
            UNROLL_S(3)
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                UNROLL_S(3)
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += coeffxp2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                                  ( ( 1-dual[1][ipart] )*field_buffer[1+iloc][1+jloc]
                                  + dual[1][ipart]*field_buffer[1+iloc][2+jloc] );
                }
            }
            Epart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        }

        //Ez(primal, primal)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *Ez2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            interp_res = 0.;
            UNROLL_S(3)
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                UNROLL_S(3)
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += coeffxp2[ipart+iloc*32] * coeffyp2[ipart+jloc*32] * field_buffer[1+iloc][1+jloc];
                }
            }
            Epart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        }

        //Bx(primal, dual)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<3 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *Bx2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            interp_res = 0.;
            UNROLL_S(3)
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                UNROLL_S(3)
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += coeffxp2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                                  ( ( 1-dual[1][ipart] )*field_buffer[1+iloc][1+jloc] +
                                  dual[1][ipart]*field_buffer[1+iloc][2+jloc] );
                }
            }
            Bpart[0][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        }

        //By(dual, primal )

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-1 ; iloc<3 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *By2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            interp_res = 0.;
            UNROLL_S(3)
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                UNROLL_S(3)
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyp2[ipart+jloc*32] *
                                  ( ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][1+jloc] +
                                  dual[0][ipart]*field_buffer[2+iloc][1+jloc] ) );
                }
            }
            Bpart[1][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        }

        //Bz(dual, dual)

        // Field buffers for vectorization (required on A64FX)
        for( int iloc=-1 ; iloc<3 ; iloc++ ) {
            for( int jloc=-1 ; jloc<3 ; jloc++ ) {
                field_buffer[iloc+1][jloc+1] = ( *Bz2D )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            }
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            interp_res = 0.;
            UNROLL_S(3)
            for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                UNROLL_S(3)
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    interp_res += coeffxd2[ipart+iloc*32] * coeffyd2[ipart+jloc*32] *
                                  ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][1+jloc] + dual[0][ipart]*field_buffer[2+iloc][1+jloc] )
                                    +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][2+jloc] + dual[0][ipart]*field_buffer[2+iloc][2+jloc] ) );
                }
            }
            Bpart[2][ipart-ipart_ref+ivect+istart[0]] = interp_res;
        } // end ipart
    }

} // END Interpolator2D2OrderV


// -----------------------------------------------------------------------------
//
//! Interpolation of all fields and currents for a single particles
//! located at istart.
//! This version is not vectorized.
//! The input parameter iend not used for now, probes are interpolated one by one for now.
//
// -----------------------------------------------------------------------------
void Interpolator2D2OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
{

    int ipart = *istart;
    int nparts( particles.numberOfParticles() );

    double *Epart[3], *Bpart[3];

    for( unsigned int k=0; k<3; k++ ) {
        Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts] );
        Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts] );
    }

    int idx[2], idxO[2];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position( 0, *istart ) * d_inv_[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position( 1, *istart ) * d_inv_[1] );
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
        delta0 = particles.position( i, ipart )*d_inv_[i];
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
void Interpolator2D2OrderV::oneField( Field **, Particles &, int *, int *, double *, double *, double *, double * )
{
    ERROR( "Single field 2D2O interpolator not available in vectorized mode" );
}

void Interpolator2D2OrderV::fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Vectorized interpolation for the envelope model is not implemented for 2D geometry" );
} // END Interpolator2D2Order


void Interpolator2D2OrderV::timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Vectorized interpolation for the envelope model is not implemented for 2D geometry" );
} // END Interpolator2D2Order


void Interpolator2D2OrderV::envelopeAndSusceptibility( ElectroMagn *, Particles &, int , double *, double *, double *, double * )
{
    ERROR( "Vectorized interpolation for the envelope model is not implemented for 2D geometry" );
} // END Interpolator2D2Order
