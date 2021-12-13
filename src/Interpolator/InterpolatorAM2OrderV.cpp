#include "InterpolatorAM2OrderV.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include "ElectroMagn.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "Particles.h"
#include "LaserEnvelope.h"
#include <complex>
#include "dcomplex.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for InterpolatorAM2OrderV
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorAM2OrderV::InterpolatorAM2OrderV( Params &params, Patch *patch ) : InterpolatorAM( params, patch )
{

    nmodes_ = params.nmodes;
    D_inv_[0] = 1.0/params.cell_length[0];
    D_inv_[1] = 1.0/params.cell_length[1];
}

//Function used in Probes
//void InterpolatorAM2OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
//{
//    int ipart = *istart;
//    
//    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
//    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );
//    
//    // Interpolate E, B
//    cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
//    cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
//    cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
//    cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
//    cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
//    cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];
//    cField2D *Jl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jl_[0];
//    cField2D *Jr = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jr_[0];
//    cField2D *Jt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jt_[0];
//    cField2D *Rho= ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_AM_[0];
//    
//    // Normalized particle position
//    double xpn = particles.position( 0, ipart ) * D_inv_[0];
//    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
//    double rpn = r * D_inv_[1];
//    complex<double> exp_mm_theta = 1. ;
//    
//    // Calculate coeffs
//    coeffs( xpn, rpn );
//    
//    int nparts( particles.size() );
//    
//    // Interpolation of El^(d,p)
//    *( ELoc+0*nparts ) = std::real( compute( &coeffld_[1], &coeffrp_[1], El, id_, jp_ ) );
//    // Interpolation of Er^(p,d)
//    *( ELoc+1*nparts ) = std::real( compute( &coefflp_[1], &coeffrd_[1], Er, ip_, jd_ ) );
//    // Interpolation of Et^(p,p)
//    *( ELoc+2*nparts ) = std::real( compute( &coefflp_[1], &coeffrp_[1], Et, ip_, jp_ ) );
//    // Interpolation of Bl^(p,d)
//    *( BLoc+0*nparts ) = std::real( compute( &coefflp_[1], &coeffrd_[1], Bl, ip_, jd_ ) );
//    // Interpolation of Br^(d,p)
//    *( BLoc+1*nparts ) = std::real( compute( &coeffld_[1], &coeffrp_[1], Br, id_, jp_ ) );
//    // Interpolation of Bt^(d,d)
//    *( BLoc+2*nparts ) = std::real( compute( &coeffld_[1], &coeffrd_[1], Bt, id_, jd_ ) );
//    // Interpolation of Jl^(d,p,p)
//    JLoc->x = std::real( compute( &coeffld_[1], &coeffrp_[1], Jl, id_, jp_ ) );
//    // Interpolation of Jr^(p,d,p)
//    JLoc->y = std::real( compute( &coefflp_[1], &coeffrd_[1], Jr, ip_, jd_ ) );
//    // Interpolation of Jt^(p,p,d)
//    JLoc->z = std::real( compute( &coefflp_[1], &coeffrp_[1], Jt, ip_, jp_ ) );
//    // Interpolation of Rho^(p,p,p)
//    ( *RhoLoc ) = std::real( compute( &coefflp_[1], &coeffrp_[1], Rho, ip_, jp_ ) );
//   
//    if (r > 0){ 
//        exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
//    } else {
//        exp_m_theta_ = 1. ;
//    }
//    for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
//        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
//        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
//        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
//        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
//        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
//        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
//        Jl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jl_[imode];
//        Jr = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jr_[imode];
//        Jt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jt_[imode];
//        Rho= ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_AM_[imode];
//        
//        exp_mm_theta *= exp_m_theta_ ;
//        
//        *( ELoc+0*nparts ) += std::real( compute( &coeffld_[1], &coeffrp_[1], El, id_, jp_ ) * exp_mm_theta ) ;
//        *( ELoc+1*nparts ) += std::real( compute( &coefflp_[1], &coeffrd_[1], Er, ip_, jd_ ) * exp_mm_theta ) ;
//        *( ELoc+2*nparts ) += std::real( compute( &coefflp_[1], &coeffrp_[1], Et, ip_, jp_ ) * exp_mm_theta ) ;
//        *( BLoc+0*nparts ) += std::real( compute( &coefflp_[1], &coeffrd_[1], Bl, ip_, jd_ ) * exp_mm_theta ) ;
//        *( BLoc+1*nparts ) += std::real( compute( &coeffld_[1], &coeffrp_[1], Br, id_, jp_ ) * exp_mm_theta ) ;
//        *( BLoc+2*nparts ) += std::real( compute( &coeffld_[1], &coeffrd_[1], Bt, id_, jd_ ) * exp_mm_theta ) ;
//        JLoc->x += std::real( compute( &coeffld_[1], &coeffrp_[1], Jl, id_, jp_ ) * exp_mm_theta ) ;
//        JLoc->y += std::real( compute( &coefflp_[1], &coeffrd_[1], Jr, ip_, jd_ ) * exp_mm_theta ) ;
//        JLoc->z += std::real( compute( &coefflp_[1], &coeffrp_[1], Jt, ip_, jp_ ) * exp_mm_theta ) ;
//        ( *RhoLoc ) += std::real( compute( &coefflp_[1], &coeffrp_[1], Rho, ip_, jp_ )* exp_mm_theta ) ;
//    }
//    double delta2 = std::real( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta_ ) * *( ELoc+2*nparts );
//    *( ELoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::real( exp_m_theta_ ) * *( ELoc+2*nparts );
//    *( ELoc+1*nparts ) = delta2 ;
//    delta2 = std::real( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta_ ) *  *( BLoc+2*nparts );
//    *( BLoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::real( exp_m_theta_ ) * *( BLoc+2*nparts );
//    *( BLoc+1*nparts ) = delta2 ;
//    delta2 = std::real( exp_m_theta_ ) * JLoc->y + std::imag( exp_m_theta_ ) * JLoc->z;
//    JLoc->z = -std::imag( exp_m_theta_ ) * JLoc->y + std::real( exp_m_theta_ ) * JLoc->z;
//    JLoc->y = delta2 ;
//    
//}
//
//// Interpolator on another field than the basic ones
//void InterpolatorAM2OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *Jxloc, double *Jyloc, double *Jzloc, double *Rholoc )
//{
//
//    // **field points to the first field of the species of interest in EM->allFields
//    // They are ordered as Jx0, Jy0, Jz0, Rho0, Jx1, Jy1, Jz1, Rho1, etc.
//
//    
//    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
//        double xpn = particles.position( 0, ipart )*D_inv_[0];
//        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
//        double rpn = r * D_inv_[1];
//        coeffs( xpn, rpn);
//        complex<double> exp_m_theta_ = 1., exp_mm_theta = 1. ;
//        if (r > 0) {
//            exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
//        }
//
//
//
//        double Jx_ = 0., Jy_ = 0., Jz_ = 0., Rho_ = 0.;
//        for( unsigned int imode = 0; imode < nmodes_ ; imode++ ) {
//            cField2D *Jl  = static_cast<cField2D *>( *(field+4*imode+0) );
//            cField2D *Jr  = static_cast<cField2D *>( *(field+4*imode+1) );
//            cField2D *Jt  = static_cast<cField2D *>( *(field+4*imode+2) );
//            cField2D *Rho = static_cast<cField2D *>( *(field+4*imode+3) );
//            Jx_  += std::real( compute( &coeffld_[1], &coeffrp_[1], Jl , id_, jp_ ) * exp_mm_theta );
//            Jy_  += std::real( compute( &coefflp_[1], &coeffrd_[1], Jr , ip_, jd_ ) * exp_mm_theta );
//            Jz_  += std::real( compute( &coefflp_[1], &coeffrp_[1], Jt , ip_, jp_ ) * exp_mm_theta );
//            Rho_ += std::real( compute( &coefflp_[1], &coeffrp_[1], Rho, ip_, jp_ ) * exp_mm_theta );
//
//            exp_mm_theta *= exp_m_theta_;
//        }
//        Jxloc [ipart] = Jx_;
//        Jyloc [ipart] = std::real( exp_m_theta_ ) * Jy_ + std::imag( exp_m_theta_ ) * Jz_;
//        Jzloc [ipart] = -std::imag( exp_m_theta_ ) * Jy_ + std::real( exp_m_theta_ ) * Jz_;
//        Rholoc[ipart] = Rho_;
//    }
//}

void InterpolatorAM2OrderV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    if( istart[0] == iend[0] ) {
        return;    //Don't treat empty cells.
    }

    int nparts( ( smpi->dynamics_invgf[ithread] ).size() );

    double * __restrict__ Epart[3];
    double * __restrict__ Bpart[3];

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ position_z = particles.getPtrPosition(2);


    double * __restrict__ deltaO[2]; //Delta is the distance of the particle from its primal node in cell size. Delta is in [-0.5, +0.5[
    std::complex<double> * __restrict__ eitheta_old; //eithetaold stores exp(i theta) of the particle before pusher. 

    int idx[2], idxO[2];
    //Primal indices are constant over the all cell
    idx[0]  = round( position_x[*istart] * D_inv_[0] );
    idxO[0] = idx[0] - i_domain_begin_ -1 ;
    idx[1] = round( sqrt( position_y[*istart]*position_y[*istart] + position_z[*istart]*position_z[*istart] ) * D_inv_[1] ) ;
    idxO[1] = idx[1] - j_domain_begin_ -1 ;

    double coeff[2][2][3][32];
    double dual[2][32]; // Size ndim. Boolean converted into double indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
    
    int vecSize = 32;
    double delta, delta2; 

   
    int cell_nparts( ( int )iend[0]-( int )istart[0] );

    std::vector<complex<double>> exp_m_theta_( vecSize), exp_mm_theta( vecSize,  1.) ;                                                          //exp(-i theta), exp(-i m theta)

    //Loop on groups of vecSize particles
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed( min( cell_nparts-ivect, vecSize ) );
        deltaO[0]   =  &(   smpi->dynamics_deltaold[ithread][0        + ivect + istart[0] - ipart_ref] );
        deltaO[1]   =  &(   smpi->dynamics_deltaold[ithread][nparts   + ivect + istart[0] - ipart_ref] );
        eitheta_old =  &( smpi->dynamics_eithetaold[ithread][           ivect + istart[0] - ipart_ref] );

        for( unsigned int k=0; k<3; k++ ) {
            Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts-ipart_ref+ivect+istart[0]] );
            Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts-ipart_ref+ivect+istart[0]] );
        }

        #pragma omp simd private(delta2, delta)
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            int ipart2 = ipart+ivect+istart[0];
            double r = sqrt( position_y[ipart2]*position_y[ipart2] + position_z[ipart2]*position_z[ipart2] );
            exp_m_theta_[ipart] = ( position_y[ipart2] - Icpx * position_z[ipart2] ) / r ;
            eitheta_old[ipart] =  2.*std::real(exp_m_theta_[ipart]) - exp_m_theta_[ipart] ;  //exp(i theta)

            // i= 0 ==> X
            //             j=0 primal
            delta = position_x[ipart2]*D_inv_[0] - (double)idx[0];
            delta2  = delta*delta;
            coeff[0][0][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[0][0][1][ipart]    = ( 0.75 - delta2 );
            coeff[0][0][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            deltaO[0][ipart] = delta;

            //              j=1 dual
            dual [0][ipart] = ( delta >= 0. );
            //delta dual = distance to dual node
            delta   = delta - dual[0][ipart] + 0.5 ;
            delta2  = delta*delta;
            coeff[0][1][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[0][1][1][ipart]    = ( 0.75 - delta2 );
            coeff[0][1][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            
            // i= 1 ==> Y
            //             j=0 primal
            delta = r * D_inv_[1] - (double)idx[1];
            delta2  = delta*delta;
            coeff[1][0][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[1][0][1][ipart]    = ( 0.75 - delta2 );
            coeff[1][0][2][ipart]    =  0.5 * ( delta2+delta+0.25 );
            deltaO[1][ipart] = delta;

            //              j=1 dual
            dual [1][ipart] = ( delta >= 0. );
            //delta dual = distance to dual node
            delta   = delta - dual[1][ipart] + 0.5 ;
            delta2  = delta*delta;
            coeff[1][1][0][ipart]    =  0.5 * ( delta2-delta+0.25 );
            coeff[1][1][1][ipart]    = ( 0.75 - delta2 );
            coeff[1][1][2][ipart]    =  0.5 * ( delta2+delta+0.25 );


        }

        double interp_res;
        double * __restrict__ coeffld = &( coeff[0][1][1][0] );
        double * __restrict__ coefflp = &( coeff[0][0][1][0] );
        double * __restrict__ coeffrp = &( coeff[1][0][1][0] );
        double * __restrict__ coeffrd = &( coeff[1][1][1][0] );


        // Local buffer to store the field components
        std::complex<double> field_buffer[4][4];

        for( unsigned int imode = 0; imode < nmodes_ ; imode++ ) {
            // Static cast of the electromagnetic fields
            cField2D * __restrict__ El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D * __restrict__ Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D * __restrict__ Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D * __restrict__ Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
            cField2D * __restrict__ Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
            cField2D * __restrict__ Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];



            // Field buffers for vectorization (required on A64FX)
            //for( int iloc=-1 ; iloc<3 ; iloc++ ) {
            //    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
            //        field_buffer[iloc+1][jloc+1] = ( *El )( idxO[0]+1+iloc, idxO[1]+1+jloc );
            //    }
            //}

            #pragma omp simd private(interp_res)
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            
                
                //El(dual, primal)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res += std::real( coeffld[ipart + iloc*32] * coeffrp[ipart + jloc*32] *
                                      ( ( 1.-dual[0][ipart] )*( *El )( idxO[0]+1+iloc, idxO[1]+1+jloc ) 
                                           + dual[0][ipart]  *( *El )( idxO[0]+2+iloc, idxO[1]+1+jloc ) ) * exp_mm_theta[ipart]) ;
                                      //( ( 1.-dual[0][ipart] )*field_buffer[1+iloc][1+jloc] 
                                      //     + dual[0][ipart]  *field_buffer[2+iloc][1+jloc] ) * exp_mm_theta[ipart]) ;
                    }
                }
                Epart[0][ipart] = interp_res;
           }

           //for( int iloc=-1 ; iloc<2 ; iloc++ ) {
           //    for( int jloc=-1 ; jloc<3 ; jloc++ ) {
           //        field_buffer[iloc+1][jloc+1] = ( *Er )( idxO[0]+1+iloc, idxO[1]+1+jloc );
           //    }
           //}
                
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
               //Er(primal, dual)
               interp_res = 0.;
               UNROLL_S(3) 
               for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                   UNROLL_S(3) 
                   for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                       interp_res +=  std::real (*( coefflp+ipart+iloc*32 ) * *( coeffrd+ipart+jloc*32 ) *
                                     ( ( 1-dual[1][ipart] )*( *Er )( idxO[0]+1+iloc, idxO[1]+1+jloc ) 
                                         + dual[1][ipart]  *( *Er )( idxO[0]+1+iloc, idxO[1]+2+jloc ) ) );
                                     //( ( 1-dual[1][ipart] )*field_buffer[1+iloc][1+jloc] 
                                     //    + dual[1][ipart]  *field_buffer[1+iloc][2+jloc] )  * exp_mm_theta[ipart]);
                   }
               }
               Epart[1][ipart] = interp_res;
           }
                
           //for( int iloc=-1 ; iloc<2 ; iloc++ ) {
           //    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
           //        field_buffer[iloc+1][jloc+1] = ( *Et )( idxO[0]+1+iloc, idxO[1]+1+jloc );
           //    }
           //}
                
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Et(primal, primal)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res +=  std::real(  *( coefflp+ipart+iloc*32 ) * *( coeffrp+ipart+jloc*32 ) * 
                                                    ( *Et )( idxO[0]+1+iloc, idxO[1]+1+jloc ) );
                                                    //field_buffer[1+iloc][1+jloc]  * exp_mm_theta[ipart]);
                    }
                }
                Epart[2][ipart] = interp_res;
           }
                
           //for( int iloc=-1 ; iloc<2 ; iloc++ ) {
           //    for( int jloc=-1 ; jloc<3 ; jloc++ ) {
           //        field_buffer[iloc+1][jloc+1] = ( *Bl )( idxO[0]+1+iloc, idxO[1]+1+jloc );
           //    }
           //}

           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Bl(primal, dual)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res +=  std::real(  *( coefflp+ipart+iloc*32 ) * *( coeffrd+ipart+jloc*32 ) *
                                      ( ( ( 1-dual[1][ipart] )*( *Bl )( idxO[0]+1+iloc, idxO[1]+1+jloc ) 
                                            + dual[1][ipart]  *( *Bl )( idxO[0]+1+iloc, idxO[1]+2+jloc ) )
                                      //( ( ( 1-dual[1][ipart] ) * field_buffer[1+iloc][1+jloc] 
                                      //      + dual[1][ipart]   * field_buffer[1+iloc][2+jloc] )
                                                )  * exp_mm_theta[ipart] );
                    }
                }
                Bpart[0][ipart] = interp_res;
           }
                
           //for( int iloc=-1 ; iloc<3 ; iloc++ ) {
           //    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
           //        field_buffer[iloc+1][jloc+1] = ( *Br )( idxO[0]+1+iloc, idxO[1]+1+jloc );
           //    }
           //}
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Br(dual, primal )
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res +=  std::real(  *( coeffld+ipart+iloc*32 ) * *( coeffrp+ipart+jloc*32 ) *
                                      ( ( ( 1-dual[0][ipart] )*( *Br )( idxO[0]+1+iloc, idxO[1]+1+jloc ) 
                                            + dual[0][ipart] * ( *Br )( idxO[0]+2+iloc, idxO[1]+1+jloc ) )
                                      //( ( ( 1-dual[0][ipart] )* field_buffer[ 1+iloc][1+jloc ] 
                                      //      + dual[0][ipart]  * field_buffer[ 2+iloc][1+jloc ] )
                                                )  * exp_mm_theta[ipart]);
                    }
                }
                Bpart[1][ipart] = interp_res;
           }
                
           //for( int iloc=-1 ; iloc<3 ; iloc++ ) {
           //    for( int jloc=-1 ; jloc<3 ; jloc++ ) {
           //        field_buffer[iloc+1][jloc+1] = ( *Bt )( idxO[0]+1+iloc, idxO[1]+1+jloc );
           //    }
           //}
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Bt(dual, dual)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res +=  std::real(  *( coeffld+ipart+iloc*32 ) * *( coeffrd+ipart+jloc*32 ) *
                                      ( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*( *Bt )( idxO[0]+1+iloc, idxO[1]+1+jloc ) 
                                                                     + dual[0][ipart]*( *Bt )( idxO[0]+2+iloc, idxO[1]+1+jloc ) )
                                        +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*( *Bt )( idxO[0]+1+iloc, idxO[1]+2+jloc ) 
                                                                     + dual[0][ipart]*( *Bt )( idxO[0]+2+iloc, idxO[1]+2+jloc ) )
                                      //( ( 1-dual[1][ipart] ) * ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][1+jloc] 
                                      //                               + dual[0][ipart]*field_buffer[2+iloc][1+jloc] )
                                      //  +    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][2+jloc] 
                                      //                               + dual[0][ipart]*field_buffer[2+iloc][2+jloc] )
                                                ) * exp_mm_theta[ipart] );
                    }
                }
                Bpart[2][ipart] = interp_res;
                exp_mm_theta[ipart] *= exp_m_theta_[ipart]; //prepare for next mode
           }
       } //end loop on modes

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            //Translate field into the cartesian y,z coordinates
            double delta2 = std::real( exp_m_theta_[ipart] ) * Epart[1][ipart] + std::imag( exp_m_theta_[ipart] ) * Epart[2][ipart];
            Epart[2][ipart] = -std::imag( exp_m_theta_[ipart] ) * Epart[1][ipart] + std::real( exp_m_theta_[ipart] ) * Epart[2][ipart];
            Epart[1][ipart] = delta2 ;
            delta2 = std::real( exp_m_theta_[ipart] ) * Bpart[1][ipart] + std::imag( exp_m_theta_[ipart] ) * Bpart[2][ipart];
            Bpart[2][ipart] = -std::imag( exp_m_theta_[ipart] ) * Bpart[1][ipart] + std::real( exp_m_theta_[ipart] ) * Bpart[2][ipart];
            Bpart[1][ipart] = delta2 ;
        }


    } //end loop on ivec
}

//void InterpolatorAM2OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
//{
//    
//    
//    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
//    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
//
//    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
//    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );
//
//    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
//    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
//    
//    // Interpolate E, B
//    cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
//    cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
//    cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
//    cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
//    cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
//    cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];
//
//    // Static cast of the envelope fields
//    Field2D *Phi = static_cast<Field2D *>( EMfields->envelope->Phi_ );
//    Field2D *GradPhil = static_cast<Field2D *>( EMfields->envelope->GradPhil_ );
//    Field2D *GradPhir = static_cast<Field2D *>( EMfields->envelope->GradPhir_ );
//    
//    // auxiliary quantities    
//    double delta2, xpn, r, rpn;
//    int nparts = particles.size() ;
//
//    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
//
//        // Normalized particle position
//        xpn = particles.position( 0, ipart ) * D_inv_[0];
//        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
//        rpn = r * D_inv_[1];
//    
//        // Calculate coeffs
//        coeffs( xpn, rpn );
//    
//        
//        // only mode 0 is used
//
//        // Interpolation of El^(d,p)
//        ( *Epart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffld_[1], &coeffrp_[1], El, id_, jp_ ) );
//        // Interpolation of Er^(p,d)
//        ( *Epart ) [ 1*nparts+ipart ]       = std::real( compute( &coefflp_[1], &coeffrd_[1], Er, ip_, jd_ ) );
//        // Interpolation of Et^(p,p)
//        ( *Epart ) [ 2*nparts+ipart ]       = std::real( compute( &coefflp_[1], &coeffrp_[1], Et, ip_, jp_ ) );
//        // Interpolation of Bl^(p,d)
//        ( *Bpart ) [ 0*nparts+ipart ]       = std::real( compute( &coefflp_[1], &coeffrd_[1], Bl, ip_, jd_ ) );
//        // Interpolation of Br^(d,p)
//        ( *Bpart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffld_[1], &coeffrp_[1], Br, id_, jp_ ) );
//        // Interpolation of Bt^(d,d)
//        ( *Bpart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffld_[1], &coeffrd_[1], Bt, id_, jd_ ) );
//        // Interpolation of Phi^(p,p)
//        ( *PHIpart ) [ 0*nparts+ipart ]     = compute( &coefflp_[1], &coeffrp_[1], Phi, ip_, jp_ ) ;
//        // Interpolation of GradPhil^(p,p)
//        ( *GradPHIpart ) [ 0*nparts+ipart ] = compute( &coefflp_[1], &coeffrp_[1], GradPhil, ip_, jp_ ) ;
//        // Interpolation of GradPhir^(p,p)
//        ( *GradPHIpart ) [ 1*nparts+ipart ] = compute( &coefflp_[1], &coeffrp_[1], GradPhir, ip_, jp_ ) ;
//        // GradPhit = 0 in cylindr_ical symmetry
//        ( *GradPHIpart ) [ 2*nparts+ipart ] = 0.;
//   
//        if (r > 0){ 
//            exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
//        } else {
//            exp_m_theta_ = 1. ;
//        }
//
//        // project on x,y,z, remember that GradPhit = 0 in cylindr_ical symmetry
//        delta2 = std::real( exp_m_theta_ ) * ( *Epart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta_ ) * ( *Epart ) [ 2*nparts+ipart ];
//        ( *Epart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_ ) * ( *Epart ) [ 1*nparts+ipart ] + std::real( exp_m_theta_ ) * ( *Epart ) [ 2*nparts+ipart ];
//        ( *Epart ) [ 1*nparts+ipart ] = delta2 ;
//        delta2 = std::real( exp_m_theta_ ) * ( *Bpart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta_ ) *  ( *Bpart ) [ 2*nparts+ipart ];
//        ( *Bpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_ ) * ( *Bpart ) [ 1*nparts+ipart ] + std::real( exp_m_theta_ ) * ( *Bpart ) [ 2*nparts+ipart ];
//        ( *Bpart ) [ 1*nparts+ipart ] = delta2 ;
//
//        delta2 = std::real( exp_m_theta_ ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ; 
//        ( *GradPHIpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_ ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ;
//        ( *GradPHIpart ) [ 1*nparts+ipart ] = delta2 ;
//
//        //Buffering of iold and delta
//        ( *iold )[ipart+0*nparts]  = ip_;
//        ( *iold )[ipart+1*nparts]  = jp_;
//        ( *delta )[ipart+0*nparts] = deltax;
//        ( *delta )[ipart+1*nparts] = deltar;
//        
//
//    }
//    
//} // END InterpolatorAM2OrderV
//
//
//void InterpolatorAM2OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
//{
//    // // Static cast of the envelope fields
//    // Static cast of the envelope fields
//    Field2D *Phi_m2Dcyl = static_cast<Field2D *>( EMfields->envelope->Phi_m );
//    Field2D *GradPhil_m2Dcyl = static_cast<Field2D *>( EMfields->envelope->GradPhil_m );
//    Field2D *GradPhir_m2Dcyl = static_cast<Field2D *>( EMfields->envelope->GradPhir_m );
//    
//    std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
//    std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );
//    
//    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
//    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
//    
//    double r, delta2, xpn, rpn;
//    //Loop on bin particles
//    int nparts =  particles.size() ;
//    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
//    
//        // Normalized particle position
//        xpn = particles.position( 0, ipart ) * D_inv_[0];
//        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
//        rpn = r * D_inv_[1];
//    
//        // Compute coefficients
//        coeffs( xpn, rpn );
//
//        // only mode 0 is used
//    
//        // -------------------------
//        // Interpolation of Phi_m^(p,p)
//        // -------------------------
//        ( *PHI_mpart )[ipart] = compute( &coefflp_[1], &coeffrp_[1], Phi_m2Dcyl, ip_, jp_ );
//    
//        // -------------------------
//        // Interpolation of GradPhi_m^(p,p), l component
//        // -------------------------
//        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coefflp_[1], &coeffrp_[1], GradPhil_m2Dcyl, ip_, jp_ );
//    
//        // -------------------------
//        // Interpolation of GradPhi_m^(p,p), r component
//        // -------------------------
//        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coefflp_[1], &coeffrp_[1], GradPhir_m2Dcyl, ip_, jp_ );
//    
//        // -------------------------
//        // Interpolation of GradPhi_m^(p,p), theta component
//        // -------------------------
//        ( *GradPHI_mpart )[ipart+2*nparts] = 0.; // zero with cylindr_ical symmetry
//
//
//        if (r > 0){ 
//            exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
//        } else {
//            exp_m_theta_ = 1. ;
//        }
//
//
//        // project on x,y,z, remember that GradPhit = 0 in cylindr_ical symmetry
//        delta2 = std::real( exp_m_theta_ ) * ( *GradPHI_mpart ) [ 1*nparts+ipart ] ; 
//        ( *GradPHI_mpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_ ) * ( *GradPHI_mpart ) [ 1*nparts+ipart ] ;
//        ( *GradPHI_mpart ) [ 1*nparts+ipart ] = delta2 ;
//
//        //Buffering of iold and delta
//        ( *iold )[ipart+0*nparts]  = ip_;
//        ( *iold )[ipart+1*nparts]  = jp_;
//      
//        ( *delta )[ipart+0*nparts] = deltax;
//        ( *delta )[ipart+1*nparts] = deltar;
//      
//    
//    }
//    
//} // END InterpolatorAM2OrderV
//
//
//void InterpolatorAM2OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
//{
//    // Static cast of the electromagnetic fields
//    Field2D *Env_A_abs_2Dcyl  = static_cast<Field2D *>( EMfields->Env_A_abs_ );
//    Field2D *Env_Chi_2Dcyl    = static_cast<Field2D *>( EMfields->Env_Chi_ );
//    Field2D *Env_E_abs_2Dcyl  = static_cast<Field2D *>( EMfields->Env_E_abs_ );
//    Field2D *Env_Ex_abs_2Dcyl = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );
//
//    // Normalized particle position
//    double xpn = particles.position( 0, ipart ) * D_inv_[0];
//    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
//    double rpn = r * D_inv_[1];
//
//    // Compute coefficients
//    coeffs( xpn, rpn );
//    
//    // -------------------------
//    // Interpolation of Env_A_abs_^(p,p)
//    // -------------------------
//    *( Env_A_abs_Loc ) = compute( &coefflp_[1], &coeffrp_[1], Env_A_abs_2Dcyl, ip_, jp_);
//  
//    // -------------------------  
//    // Interpolation of Env_Chi_^(p,p)
//    // -------------------------
//    *( Env_Chi_Loc ) = compute( &coefflp_[1], &coeffrp_[1], Env_Chi_2Dcyl, ip_, jp_);
//   
//    // -------------------------
//    // Interpolation of Env_E_abs_^(p,p)
//    // -------------------------
//    *( Env_E_abs_Loc ) = compute( &coefflp_[1], &coeffrp_[1], Env_E_abs_2Dcyl, ip_, jp_);
//
//    // -------------------------
//    // Interpolation of Env_Ex_abs_^(p,p)
//    // -------------------------
//    *( Env_Ex_abs_Loc ) = compute( &coefflp_[1], &coeffrp_[1], Env_Ex_abs_2Dcyl, ip_, jp_);
//  
//} // END InterpolatorAM2OrderV
//
//void InterpolatorAM2OrderV::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
//{
//    // Static cast of the envelope fields
//    Field2D *EnvEabs  = static_cast<Field2D*>( EMfields->Env_E_abs_ );
//    Field2D *EnvExabs = static_cast<Field2D*>( EMfields->Env_Ex_abs_ );
//    
//    std::vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
//    std::vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );
// 
//    double xpn,rpn,r;
//   
//    //Loop on bin particles
//    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
//
//        // Normalized particle position
//        xpn = particles.position( 0, ipart ) * D_inv_[0];
//        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
//        rpn = r * D_inv_[1];
//                                     
//        // Compute coefficients
//        coeffs( xpn, rpn );
// 
//        // only mode 0 is used
//    
//        // ---------------------------------
//        // Interpolation of Env_E_abs^(p,p)
//        // ---------------------------------
//        ( *EnvEabs_part )[ipart] = compute( &coefflp_[1], &coeffrp_[1], EnvEabs, ip_, jp_ );
//  
//        // ---------------------------------
//        // Interpolation of Env_Ex_abs^(p,p)
//        // ---------------------------------
//        ( *EnvExabs_part )[ipart] = compute( &coefflp_[1], &coeffrp_[1], EnvExabs, ip_, jp_ );
//    
//    }
//    
//    
//} // END InterpolatorAM2OrderV




