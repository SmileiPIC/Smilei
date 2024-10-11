#include "InterpolatorAM1OrderRuytenV.h"

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
// Creator for InterpolatorAM1OrderRuytenV
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorAM1OrderRuytenV::InterpolatorAM1OrderRuytenV( Params &params, Patch *patch ) : InterpolatorAM( patch )
{

    nmodes_ = params.nmodes;
    D_inv_[0] = 1.0/params.cell_length[0];
    D_inv_[1] = 1.0/params.cell_length[1];
    nscellr_ = params.patch_size_[1] + 1;
    oversize_[0] = params.oversize[0];
    oversize_[1] = params.oversize[1];
}

void InterpolatorAM1OrderRuytenV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell, int ipart_ref )
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
    idx[0] = scell/nscellr_+oversize_[0]+i_domain_begin_;
    idxO[0] = idx[0] - i_domain_begin_ -1 ;
    idx[1] = ( scell%nscellr_ )+oversize_[1]+j_domain_begin_;
    idxO[1] = idx[1] - j_domain_begin_ -1 ;

    double coeff[2][2][3][32];
    double dual[1][32]; // Size ndim. Boolean converted into double indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).
    
    int vecSize = 32;
    double delta, delta2; 

   
    int cell_nparts( ( int )iend[0]-( int )istart[0] );

    std::vector<complex<double>> exp_m_theta_( vecSize), exp_mm_theta( vecSize) ;                                                          //exp(-i theta), exp(-i m theta)

    //Loop on groups of vecSize particles
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed( min( cell_nparts-ivect, vecSize ) );
        deltaO[0]   =  &(   smpi->dynamics_deltaold[ithread][0        + ivect + istart[0] - ipart_ref] );
        deltaO[1]   =  &(   smpi->dynamics_deltaold[ithread][nparts   + ivect + istart[0] - ipart_ref] );
        eitheta_old =  &( smpi->dynamics_eithetaold[ithread][           ivect + istart[0] - ipart_ref] );


        #pragma omp simd private(delta2, delta)
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            int ipart2 = ipart+ivect+istart[0];
            double r = sqrt( position_y[ipart2]*position_y[ipart2] + position_z[ipart2]*position_z[ipart2] );
            exp_m_theta_[ipart] = ( position_y[ipart2] - Icpx * position_z[ipart2] ) / r ;
            exp_mm_theta[ipart] = 1. ;
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
                         

            double x_n, ypn;

            ypn = r * D_inv_[1];
            delta = ypn  - (double)idx[1];
            int rshift = (delta >= 0);
            x_n = (double)idx[1]-1. + (double)rshift;

            coeff[1][0][rshift][ipart] = (x_n+1-ypn)*(5*x_n + 2 - ypn)/(4.*x_n + 2.);
            coeff[1][0][rshift+1][ipart] = 1. - coeff[1][0][rshift][ipart];
            coeff[1][0][(rshift+2)%3][ipart] = 0.;
            deltaO[1][ipart] = delta;

            //              j=1 dual
            double x_np1 = ( double )idx[1] + 0.5 ; // This is known thanks to the cell sorting
            x_n = ( (( double )idx[1] - 0.5 > 0.) ? ( double )idx[1] - 0.5 : 0. ); 

            coeff[1][1][0][ipart] = (x_np1-ypn)*(5*x_n + 2 - ypn)/(4.*x_n + 2.);
            coeff[1][1][1][ipart] = 1. - coeff[1][0][0][ipart];
            //coeff[1][0][2][ipart] = 0.; // This coefficient is not supposed to be used
        }

        double interp_res;
        double * __restrict__ coeffld = &( coeff[0][1][1][0] );
        double * __restrict__ coefflp = &( coeff[0][0][1][0] );
        double * __restrict__ coeffrp = &( coeff[1][0][1][0] );
        double * __restrict__ coeffrd = &( coeff[1][1][1][0] );

        for( unsigned int k=0; k<3; k++ ) {
            Epart[k]= &( smpi->dynamics_Epart[ithread][k*nparts-ipart_ref+ivect+istart[0]] );
            Bpart[k]= &( smpi->dynamics_Bpart[ithread][k*nparts-ipart_ref+ivect+istart[0]] );
            #pragma omp simd 
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                Epart[k][ipart] = 0.;
                Bpart[k][ipart] = 0.;
            }
        }

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
            for( int iloc=-1 ; iloc<3 ; iloc++ ) {
                for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                    field_buffer[iloc+1][jloc+1] = ( *El )( idxO[0]+1+iloc, idxO[1]+1+jloc );
                }
            }

            #pragma omp simd private(interp_res)
            for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            
                
                //El(dual, primal)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res += std::real( coeffld[ipart + iloc*32] * coeffrp[ipart + jloc*32] *
                                      ( ( 1.-dual[0][ipart] )*field_buffer[1+iloc][1+jloc] 
                                           + dual[0][ipart]  *field_buffer[2+iloc][1+jloc] ) 
                                            * exp_mm_theta[ipart]) ; 
                    }
                }
                Epart[0][ipart] += interp_res;
           }

           for( int iloc=-1 ; iloc<2 ; iloc++ ) {
               for( int jloc=-1 ; jloc<3 ; jloc++ ) {
                   field_buffer[iloc+1][jloc+1] = ( *Er )( idxO[0]+1+iloc, idxO[1]+1+jloc );
               }
           }
                
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
               //Er(primal, dual)
               interp_res = 0.;
               UNROLL_S(3) 
               for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                   UNROLL_S(2) 
                   for( int jloc=-1 ; jloc<1 ; jloc++ ) {
                       interp_res +=  std::real (*( coefflp+ipart+iloc*32 ) * *( coeffrd+ipart+jloc*32 ) *
                                     ( 
                                       //( 1-dual[1][ipart] )*
                                       field_buffer[1+iloc][1+jloc] 
                                       //  + dual[1][ipart]  *field_buffer[1+iloc][2+jloc]
                                        )
                                           * exp_mm_theta[ipart]);
                   }
               }
               Epart[1][ipart] += interp_res;
           }
                
           for( int iloc=-1 ; iloc<2 ; iloc++ ) {
               for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                   field_buffer[iloc+1][jloc+1] = ( *Et )( idxO[0]+1+iloc, idxO[1]+1+jloc );
               }
           }
                
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Et(primal, primal)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res +=  std::real(  *( coefflp+ipart+iloc*32 ) * *( coeffrp+ipart+jloc*32 ) * 
                                                    field_buffer[1+iloc][1+jloc]  
                                                      * exp_mm_theta[ipart]);
                    }
                }
                Epart[2][ipart] += interp_res;
           }
                
           for( int iloc=-1 ; iloc<2 ; iloc++ ) {
               for( int jloc=-1 ; jloc<3 ; jloc++ ) {
                   field_buffer[iloc+1][jloc+1] = ( *Bl )( idxO[0]+1+iloc, idxO[1]+1+jloc );
               }
           }

           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Bl(primal, dual)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(2) 
                    for( int jloc=-1 ; jloc<1 ; jloc++ ) {
                        interp_res +=  std::real(  *( coefflp+ipart+iloc*32 ) * *( coeffrd+ipart+jloc*32 ) *
                                      ( ( 
                                          //( 1-dual[1][ipart] ) * 
                                          field_buffer[1+iloc][1+jloc] 
                                          //  + dual[1][ipart]   * field_buffer[1+iloc][2+jloc]
                                           )
                                                )  * exp_mm_theta[ipart] );
                    }
                }
                Bpart[0][ipart] += interp_res;
           }
                
           for( int iloc=-1 ; iloc<3 ; iloc++ ) {
               for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                   field_buffer[iloc+1][jloc+1] = ( *Br )( idxO[0]+1+iloc, idxO[1]+1+jloc );
               }
           }
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Br(dual, primal )
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(3) 
                    for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                        interp_res +=  std::real(  *( coeffld+ipart+iloc*32 ) * *( coeffrp+ipart+jloc*32 ) *
                                      ( ( ( 1-dual[0][ipart] )* field_buffer[ 1+iloc][1+jloc ] 
                                            + dual[0][ipart]  * field_buffer[ 2+iloc][1+jloc ] )
                                                )  * exp_mm_theta[ipart]);
                    }
                }
                Bpart[1][ipart] += interp_res;
           }
                
           for( int iloc=-1 ; iloc<3 ; iloc++ ) {
               for( int jloc=-1 ; jloc<3 ; jloc++ ) {
                   field_buffer[iloc+1][jloc+1] = ( *Bt )( idxO[0]+1+iloc, idxO[1]+1+jloc );
               }
           }
           #pragma omp simd private(interp_res)
           for( int ipart=0 ; ipart<np_computed; ipart++ ) {
                //Bt(dual, dual)
                interp_res = 0.;
                UNROLL_S(3) 
                for( int iloc=-1 ; iloc<2 ; iloc++ ) {
                    UNROLL_S(2) 
                    for( int jloc=-1 ; jloc<1 ; jloc++ ) {
                        interp_res +=  std::real(  *( coeffld+ipart+iloc*32 ) * *( coeffrd+ipart+jloc*32 ) *
                                      ( 
                                         //( 1-dual[1][ipart] ) * 
                                         ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][1+jloc] 
                                                                     + dual[0][ipart]*field_buffer[2+iloc][1+jloc] )
                                        //+    dual[1][ipart]  * ( ( 1-dual[0][ipart] )*field_buffer[1+iloc][2+jloc] 
                                        //                             + dual[0][ipart]*field_buffer[2+iloc][2+jloc] )
                                                ) * exp_mm_theta[ipart] );
                    }
                }
                Bpart[2][ipart] += interp_res;
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

