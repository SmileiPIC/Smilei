#include "Projector2D4OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D4OrderV
// ---------------------------------------------------------------------------------------------------------------------
Projector2D4OrderV::Projector2D4OrderV( Params &params, Patch *patch ) : Projector2D( params, patch )
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt_  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt_  = params.cell_length[1] / params.timestep;

    i_domain_begin_ = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin_ = patch->getCellStartingGlobalIndex( 1 );

    nscelly_ = params.patch_size_[1] + 1;

    oversize[0] = params.oversize[0];
    oversize[1] = params.oversize[1];

    nprimy = nscelly_ + 2*oversize[1];

    dq_inv[0] = dx_inv_;
    dq_inv[1] = dy_inv_;

    DEBUG( "cell_length "<< params.cell_length[0] );

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D4OrderV
// ---------------------------------------------------------------------------------------------------------------------
Projector2D4OrderV::~Projector2D4OrderV()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4OrderV::currentsAndDensity( double * __restrict__ Jx,
                                             double * __restrict__ Jy,
                                             double * __restrict__ Jz,
                                             double * __restrict__ rho,
                                             Particles &particles,
                                             unsigned int istart,
                                             unsigned int iend,
                                             double * __restrict__ invgf,
                                             int * __restrict__ iold,
                                             double * __restrict__ deltaold,
                                             unsigned int buffer_size,
                                             int ipart_ref, int bin_shift )
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int ipo = iold[0];
    int jpo = iold[1];

    int ipom2 = ipo-3;
    int jpom2 = jpo-3;

    int vecSize = 8;
    unsigned int bsize = 7*7*vecSize;

    double bJx[bsize] __attribute__( ( aligned( 64 ) ) );

    double DSx[56] __attribute__( ( aligned( 64 ) ) );
    double DSy[56] __attribute__( ( aligned( 64 ) ) );

    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);

    double * __restrict__ weight     = particles.getPtrWeight();
    short  * __restrict__ charge     = particles.getPtrCharge();

    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );
    int nbVec = ( iend-istart+( cell_nparts-1 )-( ( iend-istart-1 )&( cell_nparts-1 ) ) ) / vecSize;
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }

    // Jx, Jy, Jz
    currents( Jx, Jy, Jz, particles, istart, iend, invgf, iold, deltaold, buffer_size, ipart_ref, bin_shift );

    // rho^(p,p,d)
    cell_nparts = ( int )iend-( int )istart;
    #pragma omp simd
    for( unsigned int j=0; j<bsize; j++ ) {
        bJx[j] = 0.;
    }

    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed( min( cell_nparts-ivect, vecSize ) );

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            //                            X                                 //
            double pos = position_x[ivect+ipart+istart] * dx_inv_;
            int cell = round( pos );
            int cell_shift = cell-ipo-i_domain_begin_;
            double delta  = pos - ( double )cell;
            double delta2 = delta*delta;
            double delta3 = delta2*delta;
            double delta4 = delta3*delta;
            double S0 = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            double S1 = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            double S2 = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
            double S3 = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            double S4 = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            double m1 = ( cell_shift == -1 );
            double c0 = ( cell_shift ==  0 );
            double p1 = ( cell_shift ==  1 );
            DSx [          ipart] = m1 * S0                                        ;
            DSx [  vecSize+ipart] = c0 * S0 + m1 * S1                              ;
            DSx [2*vecSize+ipart] = p1 * S0 + c0 * S1 + m1* S2                     ;
            DSx [3*vecSize+ipart] =           p1 * S1 + c0* S2 + m1 * S3           ;
            DSx [4*vecSize+ipart] =                     p1* S2 + c0 * S3 + m1 * S4 ;
            DSx [5*vecSize+ipart] =                              p1 * S3 + c0 * S4 ;
            DSx [6*vecSize+ipart] =                                        p1 * S4 ;
            //                            Y                                 //
            pos = position_y[ivect+ipart+istart] * dy_inv_;
            cell = round( pos );
            cell_shift = cell-jpo-j_domain_begin_;
            delta  = pos - ( double )cell;
            delta2 = delta*delta;
            delta3 = delta2*delta;
            delta4 = delta3*delta;
            S0 = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            S1 = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            S2 = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
            S3 = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            S4 = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            m1 = ( cell_shift == -1 );
            c0 = ( cell_shift ==  0 );
            p1 = ( cell_shift ==  1 );
            DSy [          ipart] = m1 * S0                                        ;
            DSy [  vecSize+ipart] = c0 * S0 + m1 * S1                              ;
            DSy [2*vecSize+ipart] = p1 * S0 + c0 * S1 + m1* S2                     ;
            DSy [3*vecSize+ipart] =           p1 * S1 + c0* S2 + m1 * S3           ;
            DSy [4*vecSize+ipart] =                     p1* S2 + c0 * S3 + m1 * S4 ;
            DSy [5*vecSize+ipart] =                              p1 * S3 + c0 * S4 ;
            DSy [6*vecSize+ipart] =                                        p1 * S4 ;

            charge_weight[ipart] = inv_cell_volume * ( double )( charge[ivect+istart+ipart] )*
            weight[ivect+istart+ipart];
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            UNROLL_S(7)
            for( unsigned int i=0 ; i<7 ; i++ ) {
                UNROLL_S(7)
                for( unsigned int j=0 ; j<7 ; j++ ) {
                    int index( ( i*7 + j )*vecSize+ipart );
                    bJx [ index ] +=  charge_weight[ipart] * DSx[i*vecSize+ipart]*DSy[j*vecSize+ipart];
                }
            }//i
        } // END ipart (compute coeffs)

    }

    int iloc0 = ipom2*nprimy+jpom2;
    int iloc = iloc0;
    for( unsigned int i=0 ; i<7 ; i++ ) {
        #pragma omp simd
        for( unsigned int j=0 ; j<7 ; j++ ) {
            double tmpRho = 0.;
            int ilocal = ( ( i )*7+j )*vecSize;
            UNROLL(8)
            for( int ipart=0 ; ipart<8; ipart++ ) {
                tmpRho +=  bJx[ilocal+ipart];
            }
            rho [iloc + j] +=  tmpRho;
        }
        iloc += nprimy;
    }

} // END Project local current densities at dag timestep.


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4OrderV::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int /*bin_shift*/ )
{
    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3

    int iloc;
    int ny( nprimy );
    // (x,y,z) components of the current density for the macro-particle
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
        }
    }

    // variable declaration
    double xpn, ypn;
    double delta, delta2, delta3, delta4;
    // arrays used for the Esirkepov projection method
    double  Sx1[7], Sy1[7];

    for( unsigned int i=0; i<7; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
    }

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip        = round( xpn + 0.5 * ( type==1 ) );                       // index of the central node
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;

    Sx1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sx1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sx1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sx1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    delta3 = delta2*delta;
    delta4 = delta3*delta;

    Sy1[1] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
    Sy1[2] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[3] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
    Sy1[4] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
    Sy1[5] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ip -= i_domain_begin_ + 3;
    jp -= j_domain_begin_ + 3;

    for( unsigned int i=0 ; i<7 ; i++ ) {
        iloc = ( i+ip )*ny+jp;
        for( unsigned int j=0 ; j<7 ; j++ ) {
            rhoj[iloc+j] += charge_weight * Sx1[i]*Sy1[j];
        }
    }//i

} // END Project local current densities (Frozen species)


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization (WARNING: Not Vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4OrderV::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz,
                                             Particles &particles, int ipart, LocalFields Jion )
{
    Field2D *Jx2D  = static_cast<Field2D *>( Jx );
    Field2D *Jy2D  = static_cast<Field2D *>( Jy );
    Field2D *Jz2D  = static_cast<Field2D *>( Jz );


    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxip3, xpmxip4, xpmxid, xpmxid2, xpmxid3, xpmxid4;
    double ypn, ypmyjp, ypmyjp2, ypmyjp3, ypmyjp4, ypmyjd, ypmyjd2, ypmyjd3, ypmyjd4;
    double Sxp[5], Sxd[5], Syp[5], Syd[5];

    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;

    //Locate particle on the grid
    xpn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    ypn    = particles.position( 1, ipart ) * dy_inv_; // normalized distance to the first node

    // x-primal index
    ip      = round( xpn );                  // x-index of the central node
    xpmxip  = xpn - ( double )ip;            // normalized distance to the nearest grid point
    xpmxip2 = xpmxip*xpmxip;                 // square of the normalized distance to the nearest grid point
    xpmxip3 = xpmxip2*xpmxip;                // cube
    xpmxip4 = xpmxip2*xpmxip2;               // fourth-power

    // x-dual index
    id      = round( xpn+0.5 );              // x-index of the central node
    xpmxid  = xpn - ( double )id + 0.5;      // normalized distance to the nearest grid point
    xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    xpmxid3 = xpmxid2*xpmxid;                // cube
    xpmxid4 = xpmxid2*xpmxid2;               // fourth-power

    // y-primal index
    jp      = round( ypn );                  // y-index of the central node
    ypmyjp  = ypn - ( double )jp;            // normalized distance to the nearest grid point
    ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    ypmyjp3 = ypmyjp2*ypmyjp;                // cube
    ypmyjp4 = ypmyjp2*ypmyjp2;               // fourth-power

    // y-dual index
    jd      = round( ypn+0.5 );              // y-index of the central node
    ypmyjd  = ypn - ( double )jd + 0.5;      // normalized distance to the nearest grid point
    ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    ypmyjd3 = ypmyjd2*ypmyjd;                // cube
    ypmyjd4 = ypmyjd2*ypmyjd2;               // fourth-power

    Sxp[0] = dble_1_ov_384   - dble_1_ov_48  * xpmxip  + dble_1_ov_16 * xpmxip2 - dble_1_ov_12 * xpmxip3 + dble_1_ov_24 * xpmxip4;
    Sxp[1] = dble_19_ov_96   - dble_11_ov_24 * xpmxip  + dble_1_ov_4  * xpmxip2 + dble_1_ov_6  * xpmxip3 - dble_1_ov_6  * xpmxip4;
    Sxp[2] = dble_115_ov_192 - dble_5_ov_8   * xpmxip2 + dble_1_ov_4  * xpmxip4;
    Sxp[3] = dble_19_ov_96   + dble_11_ov_24 * xpmxip  + dble_1_ov_4  * xpmxip2 - dble_1_ov_6  * xpmxip3 - dble_1_ov_6  * xpmxip4;
    Sxp[4] = dble_1_ov_384   + dble_1_ov_48  * xpmxip  + dble_1_ov_16 * xpmxip2 + dble_1_ov_12 * xpmxip3 + dble_1_ov_24 * xpmxip4;

    Sxd[0] = dble_1_ov_384   - dble_1_ov_48  * xpmxid  + dble_1_ov_16 * xpmxid2 - dble_1_ov_12 * xpmxid3 + dble_1_ov_24 * xpmxid4;
    Sxd[1] = dble_19_ov_96   - dble_11_ov_24 * xpmxid  + dble_1_ov_4  * xpmxid2 + dble_1_ov_6  * xpmxid3 - dble_1_ov_6  * xpmxid4;
    Sxd[2] = dble_115_ov_192 - dble_5_ov_8   * xpmxid2 + dble_1_ov_4  * xpmxid4;
    Sxd[3] = dble_19_ov_96   + dble_11_ov_24 * xpmxid  + dble_1_ov_4  * xpmxid2 - dble_1_ov_6  * xpmxid3 - dble_1_ov_6  * xpmxid4;
    Sxd[4] = dble_1_ov_384   + dble_1_ov_48  * xpmxid  + dble_1_ov_16 * xpmxid2 + dble_1_ov_12 * xpmxid3 + dble_1_ov_24 * xpmxid4;

    Syp[0] = dble_1_ov_384   - dble_1_ov_48  * ypmyjp  + dble_1_ov_16 * ypmyjp2 - dble_1_ov_12 * ypmyjp3 + dble_1_ov_24 * ypmyjp4;
    Syp[1] = dble_19_ov_96   - dble_11_ov_24 * ypmyjp  + dble_1_ov_4  * ypmyjp2 + dble_1_ov_6  * ypmyjp3 - dble_1_ov_6  * ypmyjp4;
    Syp[2] = dble_115_ov_192 - dble_5_ov_8   * ypmyjp2 + dble_1_ov_4  * ypmyjp4;
    Syp[3] = dble_19_ov_96   + dble_11_ov_24 * ypmyjp  + dble_1_ov_4  * ypmyjp2 - dble_1_ov_6  * ypmyjp3 - dble_1_ov_6  * ypmyjp4;
    Syp[4] = dble_1_ov_384   + dble_1_ov_48  * ypmyjp  + dble_1_ov_16 * ypmyjp2 + dble_1_ov_12 * ypmyjp3 + dble_1_ov_24 * ypmyjp4;

    Syd[0] = dble_1_ov_384   - dble_1_ov_48  * ypmyjd  + dble_1_ov_16 * ypmyjd2 - dble_1_ov_12 * ypmyjd3 + dble_1_ov_24 * ypmyjd4;
    Syd[1] = dble_19_ov_96   - dble_11_ov_24 * ypmyjd  + dble_1_ov_4  * ypmyjd2 + dble_1_ov_6  * ypmyjd3 - dble_1_ov_6  * ypmyjd4;
    Syd[2] = dble_115_ov_192 - dble_5_ov_8   * ypmyjd2 + dble_1_ov_4  * ypmyjd4;
    Syd[3] = dble_19_ov_96   + dble_11_ov_24 * ypmyjd  + dble_1_ov_4  * ypmyjd2 - dble_1_ov_6  * ypmyjd3 - dble_1_ov_6  * ypmyjd4;
    Syd[4] = dble_1_ov_384   + dble_1_ov_48  * ypmyjd  + dble_1_ov_16 * ypmyjd2 + dble_1_ov_12 * ypmyjd3 + dble_1_ov_24 * ypmyjd4;

    ip  -= i_domain_begin_;
    id  -= i_domain_begin_;
    jp  -= j_domain_begin_;
    jd  -= j_domain_begin_;

    for (unsigned int i=0 ; i<5 ; i++) {
        int iploc=ip+i-2;
        int idloc=id+i-2;
        for (unsigned int j=0 ; j<5 ; j++) {
            int jploc=jp+j-2;
            int jdloc=jd+j-2;
            // Jx^(d,p)
            (*Jx2D)(idloc,jploc) += Jx_ion * Sxd[i]*Syp[j];
            // Jy^(p,d)
            (*Jy2D)(iploc,jdloc) += Jy_ion * Sxp[i]*Syd[j];
            // Jz^(p,p)
            (*Jz2D)(iploc,jploc) += Jz_ion * Sxp[i]*Syp[j];
        }
    }

} // END Project global current densities (ionize)


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector vectorized
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4OrderV::currents( double * __restrict__ Jx,
                                   double * __restrict__ Jy,
                                   double * __restrict__ Jz,
                                   Particles &particles,
                                   unsigned int istart, unsigned int iend,
                                   double * __restrict__ invgf,
                                   int * __restrict__ iold,
                                   double * __restrict__ deltaold,
                                   unsigned int buffer_size,
                                   int ipart_ref, int /*bin_shift*/ )
{

    // std::cerr << "Projection" << std::endl;

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int ipo = iold[0];
    int jpo = iold[1];

    int ipom2 = ipo-3;
    int jpom2 = jpo-3;

    int vecSize = 8;
    unsigned int bsize = 7*7*vecSize;

    double bJx[bsize] __attribute__( ( aligned( 64 ) ) );
    double bJy[bsize] __attribute__( ( aligned( 64 ) ) );
    double bJz[bsize] __attribute__( ( aligned( 64 ) ) );

    double Sx0_buff_vect[48] __attribute__( ( aligned( 64 ) ) );
    double Sy0_buff_vect[48] __attribute__( ( aligned( 64 ) ) );

    double Sx1_buff_vect[56] __attribute__( ( aligned( 64 ) ) );
    double Sy1_buff_vect[56] __attribute__( ( aligned( 64 ) ) );

    double DSx[56] __attribute__( ( aligned( 64 ) ) );
    double DSy[56] __attribute__( ( aligned( 64 ) ) );

    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    double crz_p[8]         __attribute__( ( aligned( 64 ) ) );

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);
    double * __restrict__ weight     = particles.getPtrWeight();
    short  * __restrict__ charge     = particles.getPtrCharge();

    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );
    int nbVec = ( iend-istart+( cell_nparts-1 )-( ( iend-istart-1 )&( cell_nparts-1 ) ) ) / vecSize;
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }

    // __________________________________________________________
    // For debugging

    // double Jxpart[cell_nparts][7][7];
    // double Jypart[cell_nparts][7][7];
    // double Jzpart[cell_nparts][7][7];

    // for (int ip = 0 ; ip < cell_nparts ; ip++) {
    //     for (int i = 0 ; i < 7 ; i++) {
    //         for (int j = 0 ; j < 7 ; j++) {
    //             Jxpart[ip][i][j] = 0;
    //             Jypart[ip][i][j] = 0;
    //             Jzpart[ip][i][j] = 0;
    //         }
    //     }
    // }
    // __________________________________________________________


    // Jx^(d,p,p)

    #pragma omp simd
    for( unsigned int j=0; j<bsize; j++ ) {
        bJx[j] = 0.;
        bJy[j] = 0.;
        bJz[j] = 0.;
    }

    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed( min( cell_nparts-ivect, vecSize ) );

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            // locate the particle on the primal grid at former time-step & calculate coeff. S0
            //                            X                                 //
            double delta = deltaold[ivect+ipart-ipart_ref+istart];
            double delta2 = delta*delta;
            double delta3 = delta2*delta;
            double delta4 = delta3*delta;

            Sx0_buff_vect[          ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            Sx0_buff_vect[  vecSize+ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            Sx0_buff_vect[2*vecSize+ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
            Sx0_buff_vect[3*vecSize+ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            Sx0_buff_vect[4*vecSize+ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            Sx0_buff_vect[5*vecSize+ipart] = 0.;


            //                            Y                                 //
            delta = deltaold[ivect+ipart-ipart_ref+istart+buffer_size];
            delta2 = delta*delta;
            delta3 = delta2*delta;
            delta4 = delta3*delta;

            Sy0_buff_vect[          ipart] = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            Sy0_buff_vect[  vecSize+ipart] = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            Sy0_buff_vect[2*vecSize+ipart] = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
            Sy0_buff_vect[3*vecSize+ipart] = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            Sy0_buff_vect[4*vecSize+ipart] = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            Sy0_buff_vect[5*vecSize+ipart] = 0.;



            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            //                            X                                 //
            double pos = position_x[ivect+ipart+istart] * dx_inv_;
            int cell = round( pos );
            int cell_shift = cell-ipo-i_domain_begin_;
            delta  = pos - ( double )cell;
            delta2 = delta*delta;
            delta3 = delta2*delta;
            delta4 = delta3*delta;
            double S0 = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            double S1 = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            double S2 = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
            double S3 = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            double S4 = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            double m1 = ( cell_shift == -1 );
            double c0 = ( cell_shift ==  0 );
            double p1 = ( cell_shift ==  1 );

            Sx1_buff_vect [          ipart] = m1 * S0                                        ;
            Sx1_buff_vect [  vecSize+ipart] = c0 * S0 + m1 * S1                              ;
            Sx1_buff_vect [2*vecSize+ipart] = p1 * S0 + c0 * S1 + m1* S2                     ;
            Sx1_buff_vect [3*vecSize+ipart] =           p1 * S1 + c0* S2 + m1 * S3           ;
            Sx1_buff_vect [4*vecSize+ipart] =                     p1* S2 + c0 * S3 + m1 * S4 ;
            Sx1_buff_vect [5*vecSize+ipart] =                              p1 * S3 + c0 * S4 ;
            Sx1_buff_vect [6*vecSize+ipart] =                                        p1 * S4 ;

            DSx [          ipart] = Sx1_buff_vect [          ipart]                                  ;
            DSx [  vecSize+ipart] = Sx1_buff_vect [  vecSize+ipart] - Sx0_buff_vect[0*vecSize+ipart] ;
            DSx [2*vecSize+ipart] = Sx1_buff_vect [2*vecSize+ipart] - Sx0_buff_vect[1*vecSize+ipart] ;
            DSx [3*vecSize+ipart] = Sx1_buff_vect [3*vecSize+ipart] - Sx0_buff_vect[2*vecSize+ipart] ;
            DSx [4*vecSize+ipart] = Sx1_buff_vect [4*vecSize+ipart] - Sx0_buff_vect[3*vecSize+ipart] ;
            DSx [5*vecSize+ipart] = Sx1_buff_vect [5*vecSize+ipart] - Sx0_buff_vect[4*vecSize+ipart] ;
            DSx [6*vecSize+ipart] = Sx1_buff_vect [6*vecSize+ipart]                                  ;

            // DSx [          ipart] = m1 * S0                                                                         ;
            // DSx [  vecSize+ipart] = c0 * S0 + m1 * S1                              - Sx0_buff_vect[          ipart] ;
            // DSx [2*vecSize+ipart] = p1 * S0 + c0 * S1 + m1* S2                     - Sx0_buff_vect[  vecSize+ipart] ;
            // DSx [3*vecSize+ipart] =           p1 * S1 + c0* S2 + m1 * S3           - Sx0_buff_vect[2*vecSize+ipart] ;
            // DSx [4*vecSize+ipart] =                     p1* S2 + c0 * S3 + m1 * S4 - Sx0_buff_vect[3*vecSize+ipart] ;
            // DSx [5*vecSize+ipart] =                              p1 * S3 + c0 * S4 - Sx0_buff_vect[4*vecSize+ipart] ;
            // DSx [6*vecSize+ipart] =                                        p1 * S4                                  ;

            //                            Y                                 //
            pos = position_y[ivect+ipart+istart] * dy_inv_;
            cell = round( pos );
            cell_shift = cell-jpo-j_domain_begin_;
            delta  = pos - ( double )cell;
            delta2 = delta*delta;
            delta3 = delta2*delta;
            delta4 = delta3*delta;

            S0 = dble_1_ov_384   - dble_1_ov_48  * delta  + dble_1_ov_16 * delta2
                                 - dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;
            S1 = dble_19_ov_96   - dble_11_ov_24 * delta  + dble_1_ov_4  * delta2
                                 + dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            S2 = dble_115_ov_192 - dble_5_ov_8   * delta2 + dble_1_ov_4  * delta4;
            S3 = dble_19_ov_96   + dble_11_ov_24 * delta  + dble_1_ov_4  * delta2
                                 - dble_1_ov_6  * delta3 - dble_1_ov_6  * delta4;
            S4 = dble_1_ov_384   + dble_1_ov_48  * delta  + dble_1_ov_16 * delta2
                                 + dble_1_ov_12 * delta3 + dble_1_ov_24 * delta4;

            m1 = ( cell_shift == -1 );
            c0 = ( cell_shift ==  0 );
            p1 = ( cell_shift ==  1 );

            Sy1_buff_vect [          ipart] = m1 * S0                                        ;
            Sy1_buff_vect [  vecSize+ipart] = c0 * S0 + m1 * S1                              ;
            Sy1_buff_vect [2*vecSize+ipart] = p1 * S0 + c0 * S1 + m1* S2                     ;
            Sy1_buff_vect [3*vecSize+ipart] =           p1 * S1 + c0* S2 + m1 * S3           ;
            Sy1_buff_vect [4*vecSize+ipart] =                     p1* S2 + c0 * S3 + m1 * S4 ;
            Sy1_buff_vect [5*vecSize+ipart] =                              p1 * S3 + c0 * S4 ;
            Sy1_buff_vect [6*vecSize+ipart] =                                        p1 * S4 ;

            DSy [          ipart] = Sy1_buff_vect [          ipart]                                  ;                                                                     ;
            DSy [  vecSize+ipart] = Sy1_buff_vect [  vecSize+ipart] - Sy0_buff_vect[0*vecSize+ipart] ;
            DSy [2*vecSize+ipart] = Sy1_buff_vect [2*vecSize+ipart] - Sy0_buff_vect[1*vecSize+ipart] ;
            DSy [3*vecSize+ipart] = Sy1_buff_vect [3*vecSize+ipart] - Sy0_buff_vect[2*vecSize+ipart] ;
            DSy [4*vecSize+ipart] = Sy1_buff_vect [4*vecSize+ipart] - Sy0_buff_vect[3*vecSize+ipart] ;
            DSy [5*vecSize+ipart] = Sy1_buff_vect [5*vecSize+ipart] - Sy0_buff_vect[4*vecSize+ipart] ;
            DSy [6*vecSize+ipart] = Sy1_buff_vect [6*vecSize+ipart]                                  ;

            // DSy [          ipart] = m1 * S0                                                                         ;
            // DSy [  vecSize+ipart] = c0 * S0 + m1 * S1                              - Sy0_buff_vect[          ipart] ;
            // DSy [2*vecSize+ipart] = p1 * S0 + c0 * S1 + m1* S2                     - Sy0_buff_vect[  vecSize+ipart] ;
            // DSy [3*vecSize+ipart] =           p1 * S1 + c0* S2 + m1 * S3           - Sy0_buff_vect[2*vecSize+ipart] ;
            // DSy [4*vecSize+ipart] =                     p1* S2 + c0 * S3 + m1 * S4 - Sy0_buff_vect[3*vecSize+ipart] ;
            // DSy [5*vecSize+ipart] =                              p1 * S3 + c0 * S4 - Sy0_buff_vect[4*vecSize+ipart] ;
            // DSy [6*vecSize+ipart] =                                        p1 * S4                                  ;

            charge_weight[ipart] = inv_cell_volume * ( double )( charge[ivect+istart+ipart] )*weight[ivect+istart+ipart];

            crz_p[ipart] = charge_weight[ipart]*one_third*momentum_z[ivect+istart+ipart]
                                               * invgf[ivect+istart+ipart];

        }


        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            //optrpt complains about the following loop but not unrolling it actually seems to give better result.
            double crx_p = charge_weight[ipart]*dx_ov_dt_;

            double sum[7];
            sum[0] = 0.;
            UNROLL_S(6)
            for( unsigned int k=1 ; k<7 ; k++ ) {
                sum[k] = sum[k-1]-DSx[( k-1 )*vecSize+ipart];
            }

            double tmp( crx_p * ( 0.5*DSy[ipart] ) ); //ok
            UNROLL_S(6)
            for( unsigned int i=1 ; i<7 ; i++ ) {
                bJx [( ( i )*7 )*vecSize+ipart] += sum[i]*tmp;
                //Jxpart[ivect + ipart][i][0] += sum[i]*tmp;
            }

            UNROLL_S(6)
            for( unsigned int j=1 ; j<7 ; j++ ) { // ok
                tmp = crx_p * ( Sy0_buff_vect[(j-1)*vecSize+ipart] + 0.5*DSy[j*vecSize+ipart] );
                UNROLL_S(6)
                for( unsigned int i=1 ; i<7 ; i++ ) {
                    bJx [ ( i*7+j )*vecSize+ipart ] += sum[i]*tmp;
                    //Jxpart[ivect + ipart][i][j] += sum[i]*tmp;
                }
            }//j

        } // END ipart (compute coeffs)

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            //optrpt complains about the following loop but not unrolling it actually seems to give better result.
            double cry_p = charge_weight[ipart]*dy_ov_dt_;

            double sum[7];
            sum[0] = 0.;
            UNROLL_S(6)
            for( unsigned int k=1 ; k<7 ; k++ ) {
                sum[k] = sum[k-1]-DSy[( k-1 )*vecSize+ipart];
            }

            double tmp( cry_p *0.5*DSx[ipart] );
            UNROLL_S(6)
            for( unsigned int j=1 ; j<7 ; j++ ) {
                bJy [j*vecSize+ipart] += sum[j]*tmp;
                //Jypart[ivect + ipart][0][j] += sum[j]*tmp;
            }


            UNROLL_S(6)
            for( unsigned int i=1 ; i<7 ; i++ ) {
                tmp = cry_p * ( Sx0_buff_vect[(i-1)*vecSize+ipart] + 0.5*DSx[i*vecSize+ipart] );

                UNROLL_S(6)
                for( unsigned int j=1 ; j<7 ; j++ ) {
                    bJy [ ( i*7+j )*vecSize+ipart ] += sum[j]*tmp;
                    //Jypart[ivect + ipart][i][j] += sum[j]*tmp;
                }
            }//i

        } // END ipart (compute coeffs)

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {

            bJz [ipart] += crz_p[ipart] * Sx1_buff_vect[ipart] * Sy1_buff_vect[ipart];
            double tmp( crz_p[ipart] * Sy1_buff_vect[ipart] );
            UNROLL(6)
            for( unsigned int i=1 ; i<7 ; i++ ) {
                bJz [( ( i )*7 )*vecSize+ipart] += tmp * ( 0.5*Sx0_buff_vect[(i-1)*vecSize+ipart]
                                                       + Sx1_buff_vect[i*vecSize+ipart] );
                //Jzpart[ivect+ipart][i][0] += tmp * ( 0.5*Sx0_buff_vect[(i-1)*vecSize+ipart]
                //                                       + Sx1_buff_vect[i*vecSize+ipart] );
            }

            tmp = crz_p[ipart] * Sx1_buff_vect[ipart];
            UNROLL(6)
            for( unsigned int j=1; j<7 ; j++ ) {
                bJz [j*vecSize+ipart] +=  tmp * ( 0.5*Sy0_buff_vect[(j-1)*vecSize+ipart]
                                              + Sy1_buff_vect[j*vecSize+ipart] );
                //Jzpart[ivect+ipart][0][j] +=  tmp * ( 0.5*Sy0_buff_vect[(j-1)*vecSize+ipart]
                //                              + Sy1_buff_vect[j*vecSize+ipart] );
            }

            UNROLL(6)
            for( unsigned int i=1 ; i<7 ; i++ ) {
                double tmp0( crz_p[ipart] * ( 0.5*Sx0_buff_vect[(i-1)*vecSize+ipart] + Sx1_buff_vect[i*vecSize+ipart] ) );
                double tmp1( crz_p[ipart] * ( 0.5*Sx1_buff_vect[i*vecSize+ipart] + Sx0_buff_vect[(i-1)*vecSize+ipart] ) );
                UNROLL(6)
                for( unsigned int j=1; j<7 ; j++ ) {
                    bJz [( ( i )*7+j )*vecSize+ipart] += ( Sy0_buff_vect[(j-1)*vecSize+ipart]* tmp1
                                                       + Sy1_buff_vect[j*vecSize+ipart]* tmp0 );
                    // Jzpart[ivect+ipart][i][j] += ( Sy0_buff_vect[(j-1)*vecSize+ipart]* tmp1
                    //                                    + Sy1_buff_vect[j*vecSize+ipart]* tmp0 );
                }
            }

        } // END ipart (compute coeffs)

    } // END ivect


    int iglobal0 = ipom2*nprimy+jpom2;
    int iglobal  = iglobal0;
    for( unsigned int i=1 ; i<7 ; i++ ) {
        iglobal += nprimy;
        #pragma omp simd
        for( unsigned int j=0 ; j<7 ; j++ ) {
            double tmpJx = 0.;
            int ilocal = ( ( i )*7+j )*vecSize;
            UNROLL(8)
            for( int ipart=0 ; ipart<8; ipart++ ) {
                tmpJx += bJx [ilocal+ipart];
            }
            Jx[iglobal+j] += tmpJx;
        }
    }

    iglobal = iglobal0+ipom2;
    for( unsigned int i=0 ; i<7 ; i++ ) {
        #pragma omp simd
        for( unsigned int j=1 ; j<7 ; j++ ) {
            double tmpJy = 0.;
            int ilocal = ( ( i )*7+j )*vecSize;
            UNROLL(8)
            for( int ipart=0 ; ipart<8; ipart++ ) {
                tmpJy += bJy [ilocal+ipart];
            }
            Jy[iglobal+j] += tmpJy;
        }
        iglobal += ( nprimy+1 );
    }

    // ________________________________________________________________________
    // For debugging

    // for (int ip=0 ; ip < cell_nparts ; ip++) {
    //     double summx = 0;
    //     double summy = 0;
    //     double summz = 0;
    //     std::cerr << "ipart: " << istart + ip
    //               << " x: " << position_x[ip+istart]
    //               << " y: " << position_y[ip+istart];
    //     for( unsigned int j=0 ; j<7 ; j++ ) {
    //         for( unsigned int i=0 ; i<7 ; i++ ) {
    //             summx += Jxpart[ip][i][j];
    //             summy += Jypart[ip][i][j];
    //             summz += Jzpart[ip][i][j];
    //         }
    //     }
    //     std::cerr << " sumx : " << summx
    //               << " sumy : " << summy
    //               << " sumz : " << summz;
    //     // for( unsigned int j=0 ; j<7 ; j++ ) {
    //     //     for( unsigned int i=0 ; i<7 ; i++ ) {
    //     //         std::cerr << " " << Jxpart[ip][i][j];
    //     //     }
    //     // }
    //     std::cerr << std::endl;
    // }

    // ________________________________________________________________________

    iglobal = iglobal0;
    for( unsigned int i=0 ; i<7 ; i++ ) {
        #pragma omp simd
        for( unsigned int j=0 ; j<7 ; j++ ) {
            double tmpJz( 0. );
            int ilocal = ( i*7+j )*vecSize;
            UNROLL(8)
            for( int ipart=0 ; ipart<8; ipart++ ) {
                tmpJz  +=  bJz [ilocal+ipart];
            }
            Jz[iglobal+j]  +=  tmpJz;
        }//i
        iglobal += nprimy;
    } // ipart

} // END Project vectorized


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4OrderV::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int scell, int ipart_ref )
{
    if( istart == iend ) {
        return;    //Don't treat empty cells.
    }

    //Independent of cell. Should not be here
    //{
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    //}
    int iold[2];

    iold[0] = scell/nscelly_+oversize[0];
    iold[1] = ( scell%nscelly_ )+oversize[1];

    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if( !diag_flag ) {
        if( !is_spectral ) {
            double *b_Jx =  &( *EMfields->Jx_ )( 0 );
            double *b_Jy =  &( *EMfields->Jy_ )( 0 );
            double *b_Jz =  &( *EMfields->Jz_ )( 0 );
            currents( b_Jx, b_Jy, b_Jz, particles,  istart, iend, invgf->data(), iold, &( *delta )[0], invgf->size(), ipart_ref );
        } else {
            ERROR( "TO DO with rho" );
        }

        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        double *b_Jx  = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
        double *b_Jy  = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
        double *b_Jz  = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;
        double *b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
        currentsAndDensity( b_Jx, b_Jy, b_Jz, b_rho, particles,  istart, iend,
                            invgf->data(), iold, &( *delta )[0], invgf->size(), ipart_ref );
    }
}


void Projector2D4OrderV::susceptibility( ElectroMagn *, Particles &, double , SmileiMPI *, int, int, int, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
