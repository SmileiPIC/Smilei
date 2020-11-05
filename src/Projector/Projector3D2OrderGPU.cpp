#include "Projector3D2OrderGPU.h"

#include <cmath>
#include <iostream>
#ifdef _GPU
#include <accelmath.h>
#include <openacc.h>
#endif
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D2OrderGPU
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2OrderGPU::Projector3D2OrderGPU( Params &params, Patch *patch ) : Projector3D( params, patch )
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    dz_inv_   = 1.0/params.cell_length[2];
    dz_ov_dt  = params.cell_length[2] / params.timestep;
    
    nprimz = params.n_space[2] + 2*params.oversize[2] + 1;
    nprimy = params.n_space[1] + 2*params.oversize[1] + 1;
    
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    k_domain_begin = patch->getCellStartingGlobalIndex( 2 );
    
    DEBUG( "cell_length "<< params.cell_length[0] );
    
    pxr = !params.is_pxr;
    
    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D2OrderGPU
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2OrderGPU::~Projector3D2OrderGPU()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderGPU::currents( ElectroMagn *EMfields, Particles &particles, int istart, int iend, double *invgf, int *iold, double *deltaold )
{
    double* Jx =  &( *EMfields->Jx_ )( 0 );
    double* Jy =  &( *EMfields->Jy_ )( 0 );
    double* Jz =  &( *EMfields->Jz_ )( 0 );

    double* position_x = particles.getPtrPosition(0);
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);
    short  *charge = particles.getPtrCharge();
    double *weight = particles.getPtrWeight();

    int nparts = particles.last_index.back();

    int sizeofEx = EMfields->Jx_->globalDims_;
    int sizeofEy = EMfields->Jy_->globalDims_;
    int sizeofEz = EMfields->Jz_->globalDims_;

    if (iend==istart)
        return;

    int packsize = nparts;
    int npack = (iend-istart)/packsize;
    if ( npack*packsize!=iend-istart )
        npack++;
#ifndef _GPU
    double* Sx0 = (double*) malloc( 5*packsize*sizeof(double) ); //acc_malloc
    double* Sy0 = (double*) malloc( 5*packsize*sizeof(double) ); // "
    double* Sz0 = (double*) malloc( 5*packsize*sizeof(double) );
    double* DSx = (double*) malloc( 5*packsize*sizeof(double) );
    double* DSy = (double*) malloc( 5*packsize*sizeof(double) );
    double* DSz = (double*) malloc( 5*packsize*sizeof(double) );
    double* sumX = (double*) malloc( 5*packsize*sizeof(double) ); 
#else
    double* Sx0 = (double*) acc_malloc( 5*packsize*sizeof(double) );
    double* Sy0 = (double*) acc_malloc( 5*packsize*sizeof(double) );
    double* Sz0 = (double*) acc_malloc( 5*packsize*sizeof(double) );
    double* DSx = (double*) acc_malloc( 5*packsize*sizeof(double) );
    double* DSy = (double*) acc_malloc( 5*packsize*sizeof(double) );
    double* DSz = (double*) acc_malloc( 5*packsize*sizeof(double) );
    double* sumX = (double*) acc_malloc( 5*packsize*sizeof(double) );
#endif
    for (int ipack=0 ; ipack<npack ; ipack++) {
        int istart_pack = istart+ipack*packsize;
        int iend_pack = (ipack+1)*packsize;
        iend_pack = istart + min( iend-istart, iend_pack );
#ifdef _GPU
        #pragma acc parallel present( iold[0:3*nparts], deltaold[0:3*nparts] ) deviceptr( position_x, position_y, position_z, Sx0, Sy0, Sz0, DSx, DSy, DSz )
        #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - (istart+ipack*packsize);

            // -------------------------------------
            // Variable declaration & initialization
            // -------------------------------------
    
            // (x,y,z) components of the current density for the macro-particle
            double charge_weight = inv_cell_volume * ( double )( charge[ ipart ] )*weight[ ipart ];
            double crx_p = charge_weight*dx_ov_dt;
            double cry_p = charge_weight*dy_ov_dt;
            double crz_p = charge_weight*dz_ov_dt;
    
            // variable declaration
            double xpn, ypn, zpn;
            double delta, delta2;
            // arrays used for the Esirkepov projection method
    
            for( unsigned int i=0; i<5; i++ ) {
                DSx[ipart_pack+i*packsize] = 0.;
                DSy[ipart_pack+i*packsize] = 0.;
                DSz[ipart_pack+i*packsize] = 0.;
            }
        
            // --------------------------------------------------------
            // Locate particles & Calculate Esirkepov coef. S, DS and W
            // --------------------------------------------------------
    
            // locate the particle on the primal grid at former time-step & calculate coeff. S0
            delta = deltaold[0*nparts+ipart];
            delta2 = delta*delta;
            Sx0[ipart_pack+0*packsize] = 0.;
            Sx0[ipart_pack+1*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sx0[ipart_pack+2*packsize] = 0.75-delta2;
            Sx0[ipart_pack+3*packsize] = 0.5 * ( delta2+delta+0.25 );
            Sx0[ipart_pack+4*packsize] = 0.;
    
            delta = deltaold[1*nparts+ipart];
            delta2 = delta*delta;
            Sy0[ipart_pack+0*packsize] = 0.;
            Sy0[ipart_pack+1*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sy0[ipart_pack+2*packsize] = 0.75-delta2;
            Sy0[ipart_pack+3*packsize] = 0.5 * ( delta2+delta+0.25 );
            Sy0[ipart_pack+4*packsize] = 0.;
    
            delta = deltaold[2*nparts+ipart];
            delta2 = delta*delta;
            Sz0[ipart_pack+0*packsize] = 0.;
            Sz0[ipart_pack+1*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sz0[ipart_pack+2*packsize] = 0.75-delta2;
            Sz0[ipart_pack+3*packsize] = 0.5 * ( delta2+delta+0.25 );
            Sz0[ipart_pack+4*packsize] = 0.;
    
            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            xpn = position_x[ ipart ] * dx_inv_;
            int ip = round( xpn );
            int ipo = iold[0*nparts+ipart];
            int ip_m_ipo = ip-ipo-i_domain_begin;
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            DSx[ipart_pack+(ip_m_ipo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            DSx[ipart_pack+(ip_m_ipo+2)*packsize] = 0.75-delta2;
            DSx[ipart_pack+(ip_m_ipo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );
    
            ypn = position_y[ ipart ] * dy_inv_;
            int jp = round( ypn );
            int jpo = iold[1*nparts+ipart];
            int jp_m_jpo = jp-jpo-j_domain_begin;
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            DSy[ipart_pack+(jp_m_jpo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            DSy[ipart_pack+(jp_m_jpo+2)*packsize] = 0.75-delta2;
            DSy[ipart_pack+(jp_m_jpo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );
    
            zpn = position_z[ ipart ] * dz_inv_;
            int kp = round( zpn );
            int kpo = iold[2*nparts+ipart];
            int kp_m_kpo = kp-kpo-k_domain_begin;
            delta  = zpn - ( double )kp;
            delta2 = delta*delta;
            DSz[ipart_pack+(kp_m_kpo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            DSz[ipart_pack+(kp_m_kpo+2)*packsize] = 0.75-delta2;
            DSz[ipart_pack+(kp_m_kpo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );
    
            // computes Esirkepov coefficients
            for( unsigned int i=0; i < 5; i++ ) {
                DSx[ipart_pack+i*packsize] -= Sx0[ipart_pack+i*packsize];
                DSy[ipart_pack+i*packsize] -= Sy0[ipart_pack+i*packsize];
                DSz[ipart_pack+i*packsize] -= Sz0[ipart_pack+i*packsize];
            }
    
            // ---------------------------
            // Calculate the total current
            // ---------------------------
    
            iold[ipart+0*nparts] -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
            // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
            iold[ipart+1*nparts] -= 2;
            iold[ipart+2*nparts] -= 2;
        }
    
        // Jx^(d,p,p)
#ifdef _GPU
        #pragma acc parallel deviceptr( DSx, sumX )
        #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - ipack*packsize;

            sumX[ipart_pack+0*packsize] = 0.;
            for( unsigned int k=1 ; k<5 ; k++ ) {
                sumX[ipart_pack+k*packsize] = sumX[ipart_pack+(k-1)*packsize]-DSx[ ipart_pack+(k-1)*packsize ];
            }
        }

#ifdef _GPU
        #pragma acc parallel present( iold[0:3*nparts], Jx[0:sizeofEx] ) deviceptr( charge, weight, Sy0, Sz0, DSy, DSz, sumX ) vector_length(8)
        #pragma acc loop gang worker
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - ipack*packsize;

            double charge_weight = inv_cell_volume * ( double )( charge[ ipart ] )*weight[ ipart ];
            double crx_p = charge_weight*dx_ov_dt;

            int  z_size0 = nprimz;
            int yz_size0 = nprimz*nprimy;
            int linindex0 = iold[ipart+0*nparts]*yz_size0+iold[ipart+1*nparts]*z_size0+iold[ipart+2*nparts];
#ifdef _GPU
            #pragma acc loop vector
#endif
            for( int k=0 ; k<5 ; k++ ) {
                for( int j=0 ; j<5 ; j++ ) {
                    double tmp = crx_p * ( Sy0[ipart_pack+j*packsize]*Sz0[ipart_pack+k*packsize] + 0.5*DSy[ipart_pack+j*packsize]*Sz0[ipart_pack+k*packsize] + 0.5*DSz[ipart_pack+k*packsize]*Sy0[ipart_pack+j*packsize] + one_third*DSy[ipart_pack+j*packsize]*DSz[ipart_pack+k*packsize] );
                    int idx = linindex0 + j*z_size0 + k;
                    for( int i=1 ; i<5 ; i++ ) {
                        double val = sumX[ipart_pack+(i)*packsize] * tmp;
                        int jdx = idx + i*yz_size0;
#ifdef _GPU
                        #pragma acc atomic
#endif
                        Jx [ jdx ] += val;
                    }
                }
            }//i
        }

        // Jy^(p,d,p)
#ifdef _GPU
        #pragma acc parallel deviceptr( DSy, sumX )
        #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - ipack*packsize;

            sumX[ipart_pack+0*packsize] = 0.;
            for( unsigned int k=1 ; k<5 ; k++ ) {
                sumX[ipart_pack+k*packsize] = sumX[ipart_pack+(k-1)*packsize]-DSy[ ipart_pack+(k-1)*packsize ];
            }
        }

#ifdef _GPU
        #pragma acc parallel present( iold[0:3*nparts], Jy[0:sizeofEy] ) deviceptr( charge, weight, Sx0, Sz0, DSx, DSz, sumX ) vector_length(8)
        #pragma acc loop gang worker
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - ipack*packsize;


            double charge_weight = inv_cell_volume * ( double )( charge[ ipart ] )*weight[ ipart ];
            double cry_p = charge_weight*dy_ov_dt;

            int  z_size1 = nprimz;
            int yz_size1 = nprimz*( nprimy+1 );
            int linindex1 = iold[ipart+0*nparts]*yz_size1+iold[ipart+1*nparts]*z_size1+iold[ipart+2*nparts];
#ifdef _GPU
            #pragma acc loop vector
#endif
            for( int k=0 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    double tmp = cry_p * ( Sz0[ipart_pack+k*packsize]*Sx0[ipart_pack+i*packsize] + 0.5*DSz[ipart_pack+k*packsize]*Sx0[ipart_pack+i*packsize] + 0.5*DSx[ipart_pack+i*packsize]*Sz0[ipart_pack+k*packsize] + one_third*DSz[ipart_pack+k*packsize]*DSx[ipart_pack+i*packsize] );
                    int idx = linindex1 + k + i*yz_size1;
                    for( int j=1 ; j<5 ; j++ ) {
                        double val = sumX[ipart_pack+(j)*packsize] * tmp;
                        int jdx = idx + j*z_size1;
#ifdef _GPU
                        #pragma acc atomic
#endif
                        Jy [ jdx ] += val;
                    }
                }
            }//i
        }

        // Jz^(p,p,d)
#ifdef _GPU
        #pragma acc parallel deviceptr( DSz, sumX )
        #pragma acc loop gang worker vector 
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - ipack*packsize;

            sumX[ipart_pack+0*packsize] = 0.;
            for( unsigned int k=1 ; k<5 ; k++ ) {
                sumX[ipart_pack+k*packsize] = sumX[ipart_pack+(k-1)*packsize]-DSz[ ipart_pack+(k-1)*packsize ];
            }
        }

#ifdef _GPU
        #pragma acc parallel present( iold[0:3*nparts], Jz[0:sizeofEz] ) deviceptr( charge, weight, Sx0, Sy0, DSx, DSy, sumX ) vector_length(4)
        #pragma acc loop gang worker
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            int ipart_pack = ipart - ipack*packsize;

            double charge_weight = inv_cell_volume * ( double )( charge[ ipart ] )*weight[ ipart ];
            double crz_p = charge_weight*dz_ov_dt;

            int z_size2 =  nprimz+1;
            int yz_size2 = ( nprimz+1 )*nprimy;
            int linindex2 = iold[ipart+0*nparts]*yz_size2+iold[ipart+1*nparts]*z_size2+iold[ipart+2*nparts];
#ifdef _GPU
            #pragma acc loop vector
#endif
            for( int k=1 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    for( int j=0 ; j<5 ; j++ ) {
                        double tmp = crz_p*( Sx0[ipart_pack+i*packsize]*Sy0[ipart_pack+j*packsize] + 0.5*DSx[ipart_pack+i*packsize]*Sy0[ipart_pack+j*packsize] + 0.5*DSy[ipart_pack+j*packsize]*Sx0[ipart_pack+i*packsize] + one_third*DSx[ipart_pack+i*packsize]*DSy[ipart_pack+j*packsize] );
                        int idx = linindex2 + j*z_size2 + i*yz_size2;
                        double val = sumX[ipart_pack+(k)*packsize] * tmp;
                        int jdx = idx + k;
#ifdef _GPU
                        #pragma acc atomic
#endif
                        Jz[ jdx ] += val;
                    }
                }
            }//i


        } // End for ipart

    } // End for ipack
#ifndef _GPU
    free( Sx0 ); // acc_free
    free( Sy0 );
    free( Sz0 );
    free( DSx );
    free( DSy );
    free( DSz );
    free( sumX );
#else
    acc_free( Sx0 ); // acc_free
    acc_free( Sy0 );
    acc_free( Sz0 );
    acc_free( DSx );
    acc_free( DSy );
    acc_free( DSz );
    acc_free( sumX );
#endif

} // END Project local current densities (Jx, Jy, Jz, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderGPU::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold )
{
    int nparts = particles.size();
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double crx_p = charge_weight*dx_ov_dt;
    double cry_p = charge_weight*dy_ov_dt;
    double crz_p = charge_weight*dz_ov_dt;
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double Sx0[5], Sx1[5], Sy0[5], Sy1[5], Sz0[5], Sz1[5], DSx[5], DSy[5], DSz[5];
    double tmpJx[5][5], tmpJy[5][5], tmpJz[5][5];
    
    for( unsigned int i=0; i<5; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    for( unsigned int j=0; j<5; j++ )
        for( unsigned int k=0; k<5; k++ ) {
            tmpJx[j][k] = 0.;
        }
    for( unsigned int i=0; i<5; i++ )
        for( unsigned int k=0; k<5; k++ ) {
            tmpJy[i][k] = 0.;
        }
    for( unsigned int i=0; i<5; i++ )
        for( unsigned int j=0; j<5; j++ ) {
            tmpJz[i][j] = 0.;
        }
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sx0[0] = 0.;
    Sx0[1] = 0.5 * ( delta2-delta+0.25 );
    Sx0[2] = 0.75-delta2;
    Sx0[3] = 0.5 * ( delta2+delta+0.25 );
    Sx0[4] = 0.;
    
    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sy0[0] = 0.;
    Sy0[1] = 0.5 * ( delta2-delta+0.25 );
    Sy0[2] = 0.75-delta2;
    Sy0[3] = 0.5 * ( delta2+delta+0.25 );
    Sy0[4] = 0.;
    
    delta = deltaold[2*nparts];
    delta2 = delta*delta;
    Sz0[0] = 0.;
    Sz0[1] = 0.5 * ( delta2-delta+0.25 );
    Sz0[2] = 0.75-delta2;
    Sz0[3] = 0.5 * ( delta2+delta+0.25 );
    Sz0[4] = 0.;
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn );
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin;
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sx1[ip_m_ipo+1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[ip_m_ipo+2] = 0.75-delta2;
    Sx1[ip_m_ipo+3] = 0.5 * ( delta2+delta+0.25 );
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn );
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin;
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sy1[jp_m_jpo+1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[jp_m_jpo+2] = 0.75-delta2;
    Sy1[jp_m_jpo+3] = 0.5 * ( delta2+delta+0.25 );
    
    zpn = particles.position( 2, ipart ) * dz_inv_;
    int kp = round( zpn );
    int kpo = iold[2*nparts];
    int kp_m_kpo = kp-kpo-k_domain_begin;
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    Sz1[kp_m_kpo+1] = 0.5 * ( delta2-delta+0.25 );
    Sz1[kp_m_kpo+2] = 0.75-delta2;
    Sz1[kp_m_kpo+3] = 0.5 * ( delta2+delta+0.25 );
    
    // computes Esirkepov coefficients
    for( unsigned int i=0; i < 5; i++ ) {
        DSx[i] = Sx1[i] - Sx0[i];
        DSy[i] = Sy1[i] - Sy0[i];
        DSz[i] = Sz1[i] - Sz0[i];
    }
    
    // ---------------------------
    // Calculate the total current
    // ---------------------------
    
    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;
    kpo -= 2;
    
    int iloc, jloc, kloc, linindex;
    
    // Jx^(d,p,p)
    for( unsigned int i=1 ; i<5 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                tmpJx[j][k] -= crx_p * DSx[i-1] * ( Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*DSz[k]*Sy0[j] + one_third*DSy[j]*DSz[k] );
                kloc = k+kpo;
                linindex = iloc*nprimz*nprimy+jloc*nprimz+kloc;
                Jx [linindex] += tmpJx[j][k];
            }
        }
    }//i
    
    // Jy^(p,d,p)
    for( unsigned int i=0 ; i<5 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=1 ; j<5 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                tmpJy[i][k] -= cry_p * DSy[j-1] * ( Sz0[k]*Sx0[i] + 0.5*DSz[k]*Sx0[i] + 0.5*DSx[i]*Sz0[k] + one_third*DSz[k]*DSx[i] );
                kloc = k+kpo;
                linindex = iloc*nprimz*( nprimy+1*pxr )+jloc*nprimz+kloc;
                Jy [linindex] += tmpJy[i][k]; //
            }
        }
    }//i
    
    // Jz^(p,p,d)
    for( unsigned int i=0 ; i<5 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=1 ; k<5 ; k++ ) {
                tmpJz[i][j] -= crz_p * DSz[k-1] * ( Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*DSy[j]*Sx0[i] + one_third*DSx[i]*DSy[j] );
                kloc = k+kpo;
                linindex = iloc*( nprimz+1*pxr )*nprimy+jloc*( nprimz+1*pxr )+kloc;
                Jz [linindex] += tmpJz[i][j]; //
            }
        }
    }//i
    
    // Rho^(p,p,p)
    for( unsigned int i=0 ; i<5 ; i++ ) {
        iloc = i+ipo;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            jloc = j+jpo;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                kloc = k+kpo;
                linindex = iloc*nprimz*nprimy+jloc*nprimz+kloc;
                rho[linindex] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i
    
} // END Project local densities (Jx, Jy, Jz, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderGPU::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{
    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc, jloc;
    int ny( nprimy ), nz( nprimz ), nyz;
    // (x,y,z) components of the current density for the macro-particle
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2;
    double Sx1[5], Sy1[5], Sz1[5]; // arrays used for the Esirkepov projection method
    
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    if( type > 0 ) {
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
            ny ++;
        } else {
            charge_weight *= particles.momentum( 2, ipart );
            nz ++;
        }
    }
    nyz = ny*nz;
    
// Initialize all current-related arrays to zero
    for( unsigned int i=0; i<5; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn + 0.5*( type==1 ) );
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sx1[1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[2] = 0.75-delta2;
    Sx1[3] = 0.5 * ( delta2+delta+0.25 );
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sy1[1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[2] = 0.75-delta2;
    Sy1[3] = 0.5 * ( delta2+delta+0.25 );
    
    zpn = particles.position( 2, ipart ) * dz_inv_;
    int kp = round( zpn + 0.5*( type==3 ) );
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    Sz1[1] = 0.5 * ( delta2-delta+0.25 );
    Sz1[2] = 0.75-delta2;
    Sz1[3] = 0.5 * ( delta2+delta+0.25 );
    
    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin + 2;
    jp -= j_domain_begin + 2;
    kp -= k_domain_begin + 2;
    
    for( unsigned int i=0 ; i<5 ; i++ ) {
        iloc = ( i+ip ) * nyz;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            jloc = ( jp+j ) * nz;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                rhoj[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i
    
} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderGPU::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field3D *Jx3D  = static_cast<Field3D *>( Jx );
    Field3D *Jy3D  = static_cast<Field3D *>( Jy );
    Field3D *Jz3D  = static_cast<Field3D *>( Jz );
    
    
    //Declaration of local variables
    int ip, id, jp, jd, kp, kd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double zpn, zpmzkp, zpmzkp2, zpmzkd, zpmzkd2;
    double Sxp[3], Sxd[3], Syp[3], Syd[3], Szp[3], Szd[3];
    
    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;
    
    //Locate particle on the grid
    xpn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    ypn    = particles.position( 1, ipart ) * dy_inv_; // normalized distance to the first node
    zpn    = particles.position( 2, ipart ) * dz_inv_; // normalized distance to the first node
    
    // x-primal index
    ip      = round( xpn );                  // x-index of the central node
    xpmxip  = xpn - ( double )ip;            // normalized distance to the nearest grid point
    xpmxip2 = xpmxip*xpmxip;                 // square of the normalized distance to the nearest grid point
    
    // x-dual index
    id      = round( xpn+0.5 );              // x-index of the central node
    xpmxid  = xpn - ( double )id + 0.5;      // normalized distance to the nearest grid point
    xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    
    // y-primal index
    jp      = round( ypn );                  // y-index of the central node
    ypmyjp  = ypn - ( double )jp;            // normalized distance to the nearest grid point
    ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    
    // y-dual index
    jd      = round( ypn+0.5 );              // y-index of the central node
    ypmyjd  = ypn - ( double )jd + 0.5;      // normalized distance to the nearest grid point
    ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    
    // z-primal index
    kp      = round( zpn );                  // z-index of the central node
    zpmzkp  = zpn - ( double )kp;            // normalized distance to the nearest grid point
    zpmzkp2 = zpmzkp*zpmzkp;                 // square of the normalized distance to the nearest grid point
    
    // z-dual index
    kd      = round( zpn+0.5 );              // z-index of the central node
    zpmzkd  = zpn - ( double )kd + 0.5;      // normalized distance to the nearest grid point
    zpmzkd2 = zpmzkd*zpmzkd;                 // square of the normalized distance to the nearest grid point
    
    Sxp[0] = 0.5 * ( xpmxip2-xpmxip+0.25 );
    Sxp[1] = ( 0.75-xpmxip2 );
    Sxp[2] = 0.5 * ( xpmxip2+xpmxip+0.25 );
    
    Sxd[0] = 0.5 * ( xpmxid2-xpmxid+0.25 );
    Sxd[1] = ( 0.75-xpmxid2 );
    Sxd[2] = 0.5 * ( xpmxid2+xpmxid+0.25 );
    
    Syp[0] = 0.5 * ( ypmyjp2-ypmyjp+0.25 );
    Syp[1] = ( 0.75-ypmyjp2 );
    Syp[2] = 0.5 * ( ypmyjp2+ypmyjp+0.25 );
    
    Syd[0] = 0.5 * ( ypmyjd2-ypmyjd+0.25 );
    Syd[1] = ( 0.75-ypmyjd2 );
    Syd[2] = 0.5 * ( ypmyjd2+ypmyjd+0.25 );
    
    Szp[0] = 0.5 * ( zpmzkp2-zpmzkp+0.25 );
    Szp[1] = ( 0.75-zpmzkp2 );
    Szp[2] = 0.5 * ( zpmzkp2+zpmzkp+0.25 );
    
    Szd[0] = 0.5 * ( zpmzkd2-zpmzkd+0.25 );
    Szd[1] = ( 0.75-zpmzkd2 );
    Szd[2] = 0.5 * ( zpmzkd2+zpmzkd+0.25 );
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    kp  -= k_domain_begin;
    kd  -= k_domain_begin;
    
    for( unsigned int i=0 ; i<3 ; i++ ) {
        int iploc=ip+i-1;
        int idloc=id+i-1;
        for( unsigned int j=0 ; j<3 ; j++ ) {
            int jploc=jp+j-1;
            int jdloc=jd+j-1;
            for( unsigned int k=0 ; k<3 ; k++ ) {
                int kploc=kp+k-1;
                int kdloc=kd+k-1;
                // Jx^(d,p,p)
                ( *Jx3D )( idloc, jploc, kploc ) += Jx_ion * Sxd[i]*Syp[j]*Szp[k];
                // Jy^(p,d,p)
                ( *Jy3D )( iploc, jdloc, kploc ) += Jy_ion * Sxp[i]*Syd[j]*Szp[k];
                // Jz^(p,p,d)
                ( *Jz3D )( iploc, jploc, kdloc ) += Jz_ion * Sxp[i]*Syp[j]*Szd[k];
            }//k
        }//j
    }//i
    
    
    
} // END Project global current densities (ionize)

//Wrapper for projection
void Projector3D2OrderGPU::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
{
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    Jx_  =  &( *EMfields->Jx_ )( 0 );
    Jy_  =  &( *EMfields->Jy_ )( 0 );
    Jz_  =  &( *EMfields->Jz_ )( 0 );
    rho_ =  &( *EMfields->rho_ )( 0 );
    
    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if( !diag_flag ) {
        if( !is_spectral ) {        
            currents( EMfields, particles, istart, iend, &( *invgf )[0], &( *iold )[0], &( *delta )[0] );
        } else {
            for( int ipart=istart ; ipart<iend; ipart++ ) {
                currentsAndDensity( Jx_, Jy_, Jz_, rho_, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
            }
        }
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        double *b_Jx  = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
        double *b_Jy  = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
        double *b_Jz  = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;
        double *b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
        currents( EMfields, particles, istart, iend, &( *invgf )[0], &( *iold )[0], &( *delta )[0] );
        //for( int ipart=istart ; ipart<iend; ipart++ ) {
        //    currentsAndDensity( b_Jx, b_Jy, b_Jz, b_rho, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart] );
        //}
    }
    
}
// Projector for susceptibility used as source term in envelope equation
void Projector3D2OrderGPU::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )

{
    double *Chi_envelope = &( *EMfields->Env_Chi_ )( 0 );
    
    std::vector<double> *Epart       = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Phipart     = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPhipart = &( smpi->dynamics_GradPHIpart[ithread] );
    std::vector<double> *inv_gamma_ponderomotive = &( smpi->dynamics_inv_gamma_ponderomotive[ithread] );
    
    
    int iloc, jloc;
    
    double momentum[3];
    
    double gamma_ponderomotive, gamma0, gamma0_sq;
    double charge_over_mass_dts2, charge_sq_over_mass_sq_dts4, charge_sq_over_mass_sq;
    double pxsm, pysm, pzsm;
    double one_over_mass=1./species_mass;
    
    int nparts = particles.size();
    double *Ex       = &( ( *Epart )[0*nparts] );
    double *Ey       = &( ( *Epart )[1*nparts] );
    double *Ez       = &( ( *Epart )[2*nparts] );
    double *Phi      = &( ( *Phipart )[0*nparts] );
    double *GradPhix = &( ( *GradPhipart )[0*nparts] );
    double *GradPhiy = &( ( *GradPhipart )[1*nparts] );
    double *GradPhiz = &( ( *GradPhipart )[2*nparts] );
    
    for( int ipart=istart ; ipart<iend; ipart++ ) {//Loop on bin particles
    
    
        charge_over_mass_dts2       = ( double )( particles.charge( ipart ) )*dts2*one_over_mass;
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_sq_dts4 = ( double )( particles.charge( ipart ) )*( double )( particles.charge( ipart ) )*dts4*one_over_mass*one_over_mass;
        // (charge over mass)^2
        charge_sq_over_mass_sq      = ( double )( particles.charge( ipart ) )*( double )( particles.charge( ipart ) )*one_over_mass*one_over_mass;
        
        for( int i = 0 ; i<3 ; i++ ) {
            momentum[i] = particles.momentum( i, ipart );
        }
        
        // compute initial ponderomotive gamma
        gamma0_sq = 1. + momentum[0]*momentum[0]+ momentum[1]*momentum[1] + momentum[2]*momentum[2] + *( Phi+ipart )*charge_sq_over_mass_sq ;
        gamma0    = sqrt( gamma0_sq ) ;
        
        // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
        pxsm = ( gamma0 * charge_over_mass_dts2*( *( Ex+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhix+ipart ) ) ) * momentum[0] / gamma0_sq;
        pysm = ( gamma0 * charge_over_mass_dts2*( *( Ey+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhiy+ipart ) ) ) * momentum[1] / gamma0_sq;
        pzsm = ( gamma0 * charge_over_mass_dts2*( *( Ez+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhiz+ipart ) ) ) * momentum[2] / gamma0_sq;
        
        // update of gamma ponderomotive
        gamma_ponderomotive = gamma0 + ( pxsm+pysm+pzsm )*0.5 ;
        // buffer inverse of ponderomotive gamma to use it in ponderomotive momentum pusher
        ( *inv_gamma_ponderomotive )[ipart] = 1./gamma_ponderomotive;
        
        // susceptibility for the macro-particle
        double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*( double )( particles.charge( ipart ) )*particles.weight( ipart )*one_over_mass/gamma_ponderomotive;
        
        // variable declaration
        double xpn, ypn, zpn;
        double delta, delta2;
        double Sx1[5], Sy1[5], Sz1[5]; // arrays used for the Esirkepov projection method
        
        // Initialize all current-related arrays to zero
        for( unsigned int i=0; i<5; i++ ) {
            Sx1[i] = 0.;
            Sy1[i] = 0.;
            Sz1[i] = 0.;
        }
        
        // --------------------------------------------------------
        // Locate particles & Calculate Esirkepov coef. S, DS and W
        // --------------------------------------------------------
        
        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        xpn = particles.position( 0, ipart ) * dx_inv_;
        int ip = round( xpn );
        delta  = xpn - ( double )ip;
        delta2 = delta*delta;
        Sx1[1] = 0.5 * ( delta2-delta+0.25 );
        Sx1[2] = 0.75-delta2;
        Sx1[3] = 0.5 * ( delta2+delta+0.25 );
        
        ypn = particles.position( 1, ipart ) * dy_inv_;
        int jp = round( ypn );
        delta  = ypn - ( double )jp;
        delta2 = delta*delta;
        Sy1[1] = 0.5 * ( delta2-delta+0.25 );
        Sy1[2] = 0.75-delta2;
        Sy1[3] = 0.5 * ( delta2+delta+0.25 );
        
        zpn = particles.position( 2, ipart ) * dz_inv_;
        int kp = round( zpn );
        delta  = zpn - ( double )kp;
        delta2 = delta*delta;
        Sz1[1] = 0.5 * ( delta2-delta+0.25 );
        Sz1[2] = 0.75-delta2;
        Sz1[3] = 0.5 * ( delta2+delta+0.25 );
        
        // ---------------------------
        // Calculate the total susceptibility
        // ---------------------------
        ip -= i_domain_begin + 2;
        jp -= j_domain_begin + 2;
        kp -= k_domain_begin + 2;
        
        for( unsigned int i=0 ; i<5 ; i++ ) { // i loop
            iloc = ( i+ip )*nprimz*nprimy;
            for( unsigned int j=0 ; j<5 ; j++ ) { // j loop
                jloc = ( jp+j )*nprimz;
                for( unsigned int k=0 ; k<5 ; k++ ) { // k loop
                    Chi_envelope[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
                } // end k loop
            } // end j loop
        } // end i loop
        
    }
    
}
