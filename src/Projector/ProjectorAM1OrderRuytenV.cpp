#include "ProjectorAM1OrderRuytenV.h"

#include <cmath>
#include <iostream>
#include <complex>
#include "dcomplex.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"
#include "PatchAM.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ProjectorAM1OrderRuytenV
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM1OrderRuytenV::ProjectorAM1OrderRuytenV( Params &params, Patch *patch ) : ProjectorAM( params, patch )
{
    dt = params.timestep;
    dr = params.cell_length[1];
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt_  = params.cell_length[0] / params.timestep;
    one_ov_dt  = 1.0 / params.timestep;
    dr_inv_   = 1.0/dr;
    dr_ov_dt_  = dr / dt;

    i_domain_begin_ = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin_ = patch->getCellStartingGlobalIndex( 1 );

    nscellr_ = params.patch_size_[1] + 1;
    oversize_[0] = params.oversize[0];
    oversize_[1] = params.oversize[1];
    nprimr_ = nscellr_ + 2*oversize_[1];
    npriml_ = params.patch_size_[0] + 1 + 2*oversize_[0];

    Nmode_=params.nmodes;
    dq_inv_[0] = dl_inv_;
    dq_inv_[1] = dr_inv_;

    invR_ = &((static_cast<PatchAM *>( patch )->invR)[0]);
    invRd_ = &((static_cast<PatchAM *>( patch )->invRd)[0]);

    DEBUG( "cell_length "<< params.cell_length[0] );

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorAM1OrderRuytenV
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM1OrderRuytenV::~ProjectorAM1OrderRuytenV()
{
}



// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1OrderRuytenV::currentsAndDensity( ElectroMagnAM *emAM,
                                   Particles &particles,
                                   unsigned int istart,
                                   unsigned int iend,
                                   double * __restrict__ invgf,
                                   int * __restrict__ iold,
                                   double * __restrict__ deltaold,
                                   std::complex<double> * __restrict__ array_eitheta_old,
                                   int npart_total,
                                   int ipart_ref,
                                   int ispec )
{

    currents( emAM, particles,  istart, iend, invgf, iold, deltaold, array_eitheta_old, npart_total, ipart_ref, ispec+1 ); //ispec+1 is passed as a marker of diag

    int ipo = iold[0];
    int jpo = iold[1];
    int ipom2 = ipo-2;
    int jpom2 = jpo-2;

    int vecSize = 8;
    int bsize = 5*5*vecSize*Nmode_;

    std::complex<double> brho[bsize] __attribute__( ( aligned( 64 ) ) );

    double Sl0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double Sr0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double DSl[40] __attribute__( ( aligned( 64 ) ) );
    double DSr[40] __attribute__( ( aligned( 64 ) ) );
    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    double r_bar[8] __attribute__( ( aligned( 64 ) ) );
    complex<double> * __restrict__ rho;

    double *invR_local = &(invR_[jpom2]);

    // Pointer for GPU and vectorization on ARM processors
    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ position_z = particles.getPtrPosition(2);
    double * __restrict__ weight     = particles.getPtrWeight();
    short  * __restrict__ charge     = particles.getPtrCharge();
    int * __restrict__ cell_keys  = particles.getPtrCellKeys();

    #pragma omp simd
    for( unsigned int j=0; j<200*Nmode_; j++ ) {
        brho[j] = 0.;
    }

    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );

    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed = min( cell_nparts-ivect, vecSize );
        int istart0 = ( int )istart + ivect;
        complex<double> e_bar[8], e_delta_m1[8];

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            compute_distances( position_x, position_y, position_z, cell_keys, npart_total, ipart, istart0, ipart_ref, deltaold, array_eitheta_old, iold, Sl0_buff_vect, Sr0_buff_vect, DSl, DSr, r_bar, e_bar, e_delta_m1 );
            charge_weight[ipart] = inv_cell_volume * ( double )( charge[istart0+ipart] )*weight[istart0+ipart];
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeRho( ipart, charge_weight, DSl, DSr, Sl0_buff_vect, Sr0_buff_vect, brho, invR_local, e_bar);
        }
    }

    int iloc0 = ipom2*nprimr_+jpom2;
    for( unsigned int imode=0; imode<( unsigned int )Nmode_; imode++ ) {
        rho =  &( *emAM->rho_AM_[imode] )( 0 );
        unsigned int n_species = emAM->rho_AM_s.size() / Nmode_;
        unsigned int ifield = imode*n_species+ispec;
        rho  = emAM->rho_AM_s    [ifield] ? &( * ( emAM->rho_AM_s    [ifield] ) )( 0 ) : &( *emAM->rho_AM_    [imode] )( 0 ) ;
        int iloc = iloc0;
        for( unsigned int i=0 ; i<5 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=0 ; j<5 ; j++ ) {
                complex<double> tmprho( 0. );
                int ilocal = ( i*5+j )*vecSize;
                UNROLL(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmprho += brho [200*imode + ilocal+ipart];
                }
                rho[iloc+j] += tmprho;
            }
            iloc += nprimr_;
        }
    }

} // END Project local current densities at dag timestep.

// ---------------------------------------------------------------------------------------------------------------------
//! Project for diags and frozen species -
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1OrderRuytenV::basicForComplex( complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode )
{
    //Warning : this function is not charge conserving.
    // This function also assumes that particles position is evaluated at the same time as currents which is usually not true (half time-step difference).
    // It will therefore fail to evaluate the current accurately at t=0 if a plasma is already in the box.



    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int iloc, nr( nprimr_ );
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );

    if( type > 0 ) { //if current density
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
        if( type == 1 ) { //if Jl
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) { //if Jr
            charge_weight *= ( particles.momentum( 1, ipart )*particles.position( 1, ipart ) + particles.momentum( 2, ipart )*particles.position( 2, ipart ) )/ r ;
            nr++;
        } else { //if Jt
            charge_weight *= ( -particles.momentum( 1, ipart )*particles.position( 2, ipart ) + particles.momentum( 2, ipart )*particles.position( 1, ipart ) ) / r ;
        }
    }

    complex<double> e_theta = ( particles.position( 1, ipart ) + Icpx*particles.position( 2, ipart ) )/r;
    complex<double> C_m = 1.;
    if( imode > 0 ) {
        C_m = 2.;
    }
    for( unsigned int i=0; i<( unsigned int )imode; i++ ) {
        C_m *= e_theta;
    }

    double xpn, ypn;
    double delta, delta2;
    double Sl1[5], Sr1[5];

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dl_inv_;
    int ip = round( xpn + 0.5 * ( type==1 ) );
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sl1[1] = 0.5 * ( delta2-delta+0.25 );
    Sl1[2] = 0.75-delta2;
    Sl1[3] = 0.5 * ( delta2+delta+0.25 );
    ypn = r * dr_inv_ ;
    int jp = round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sr1[1] = 0.5 * ( delta2-delta+0.25 );
    Sr1[2] = 0.75-delta2;
    Sr1[3] = 0.5 * ( delta2+delta+0.25 );

    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin_ + 2;
    jp -= j_domain_begin_ + 2;

    if( type != 2 ) {
        for( unsigned int i=1 ; i<4 ; i++ ) {
            iloc = ( i+ip )*nr+jp;
            for( unsigned int j=1 ; j<4 ; j++ ) {
                rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR_[j+jp];
            }
        }//i
    } else {
        for( unsigned int i=1 ; i<4 ; i++ ) {
            iloc = ( i+ip )*nr+jp;
            for( unsigned int j=1 ; j<4 ; j++ ) {
                rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invRd_[j+jp];
            }
        }//i
    }
} // END Project for diags local current densities

// Apply boundary conditions on axis for currents and densities
void ProjectorAM1OrderRuytenV::axisBC(ElectroMagnAM *emAM, bool diag_flag )
{

   for (unsigned int imode=0; imode < Nmode_; imode++){

       std::complex<double> *rhoj = &( *emAM->rho_AM_[imode] )( 0 );
       std::complex<double> *Jl = &( *emAM->Jl_[imode] )( 0 );
       std::complex<double> *Jr = &( *emAM->Jr_[imode] )( 0 );
       std::complex<double> *Jt = &( *emAM->Jt_[imode] )( 0 );

       apply_axisBC(rhoj, Jl, Jr, Jt, imode, diag_flag);
   }

   if (diag_flag){
       unsigned int n_species = emAM->Jl_s.size() / Nmode_;
       for( unsigned int imode = 0 ; imode < emAM->Jl_.size() ; imode++ ) {
           for( unsigned int ispec = 0 ; ispec < n_species ; ispec++ ) {
               unsigned int ifield = imode*n_species+ispec;
               complex<double> *Jl  = emAM->Jl_s    [ifield] ? &( * ( emAM->Jl_s    [ifield] ) )( 0 ) : NULL ;
               complex<double> *Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : NULL ;
               complex<double> *Jt  = emAM->Jt_s    [ifield] ? &( * ( emAM->Jt_s    [ifield] ) )( 0 ) : NULL ;
               complex<double> *rho = emAM->rho_AM_s[ifield] ? &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) : NULL ;
               apply_axisBC( rho , Jl, Jr, Jt, imode, diag_flag );
           }
       }
   }
}

void ProjectorAM1OrderRuytenV::apply_axisBC(std::complex<double> *rhoj,std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, unsigned int imode, bool diag_flag )
{
   // Mode 0 contribution "below axis" is added.
   // Non zero modes are substracted because a particle sitting exactly on axis has a non defined theta and can not contribute to a theta dependent mode. 
   double sign = (imode == 0) ? 1 : -1 ;

   if (diag_flag && rhoj) {
       for( unsigned int i=2 ; i<npriml_*nprimr_+2; i+=nprimr_ ) {
           //Fold rho
           for( unsigned int j=1 ; j<3; j++ ) {
               rhoj[i+j] += sign * rhoj[i-j];
           }
           //Apply BC
           if (imode > 0){
               rhoj[i] = 0.;
               rhoj[i-1]  = - rhoj[i+1]; // Zero Jl mode > 0 on axis.
           } else {
               rhoj[i] = rhoj[i+1]; //This smoothing is just for cosmetics on the picture, rho has no influence on the results.
               rhoj[i-1]  = rhoj[i+1]; // Non zero Jl mode > 0 on axis.
           }
       }
   }

   if (Jl) {
       for( unsigned int i=2 ; i<(npriml_+1)*nprimr_+2; i+=nprimr_ ) {
           //Fold Jl
           for( unsigned int j=1 ; j<3; j++ ) {
               Jl [i+j] +=  sign * Jl[i-j]; //Add even modes, substract odd modes since el(theta=0 = el(theta=pi) at all r.
            }
            if (imode > 0){
                Jl [i] = 0. ;
                Jl[i-1]   =  -Jl[i+1]; // Zero Jl mode > 0 on axis.
           } else {
		//Jl mode 0 on axis should be left as is. It looks over estimated but it might be necessary to conserve a correct divergence and a proper evaluation on the field on axis.
                Jl [i-1] =  Jl [i+1] ; // Non zero Jl mode 0 on axis.
           }
       }
   }

   if (Jt && Jr) {
       for( unsigned int i=0 ; i<npriml_; i++ ) {
           int iloc = i*nprimr_+2;
           int ilocr = i*(nprimr_+1)+3;
           //Fold Jt
           for( unsigned int j=1 ; j<3; j++ ) {
               Jt [iloc+j] += sign * Jt[iloc-j]; 
           }
           for( unsigned int j=0 ; j<3; j++ ) {
               Jr [ilocr+2-j] += sign * Jr [ilocr-3+j];
           }

           if (imode == 1){
               Jt [iloc]= -Icpx/8.*( 9.*Jr[ilocr]- Jr[ilocr+1]);// Jt mode 1 = -I Jr mode 1 on axis to keep div(J) = 0.
               Jr [ilocr-1] = Jr [ilocr]; // Jr mode 1 is non zero on axis.
           } else{
               Jt [iloc] = 0. ; // only mode 1 is non zero on axis
               Jt [iloc-1] = -Jt [iloc+1]; // only mode 1 is non zero on axis
               Jr [ilocr-1] = -Jr [ilocr]; // only mode 1 is non zero on axis
           }
       }
   }
}

void ProjectorAM1OrderRuytenV::axisBCEnvChi( double *EnvChi )
{
    double sign = 1.;
    int imode = 0;
    for (int i=0; i< imode; i++) sign *= -1;
    if (EnvChi) {
        for( unsigned int i=2 ; i<npriml_*nprimr_+2; i+=nprimr_ ) {
            //Fold EnvChi
            //for( unsigned int j=1 ; j<3; j++ ) {
            //    EnvChi[i+j] += sign * EnvChi[i-j];
            //    EnvChi[i-j]  = sign * EnvChi[i+j];
            //}
            //EnvChi[i] = (4.*EnvChi[i+1] - EnvChi[i+2])/3.;

            EnvChi[i]   = EnvChi[i+1];
            for( unsigned int j=1 ; j<3; j++ ) {
                EnvChi[i-j]  = sign * EnvChi[i+j];
            }

        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization (WARNING: Not Vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1OrderRuytenV::ionizationCurrents( Field */*Jx*/, Field */*Jy*/, Field */*Jz*/, Particles &/*particles*/, int /*ipart*/, LocalFields /*Jion*/ )
{
/*    Field2D *Jx2D  = static_cast<Field2D *>( Jx );
    Field2D *Jy2D  = static_cast<Field2D *>( Jy );
    Field2D *Jz2D  = static_cast<Field2D *>( Jz );


    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double Sxp[3], Sxd[3], Syp[3], Syd[3];

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

    ip  -= i_domain_begin_;
    id  -= i_domain_begin_;
    jp  -= j_domain_begin_;
    jd  -= j_domain_begin_;

    for( unsigned int i=0 ; i<3 ; i++ ) {
        int iploc=ip+i-1;
        int idloc=id+i-1;
        for( unsigned int j=0 ; j<3 ; j++ ) {
            int jploc=jp+j-1;
            int jdloc=jd+j-1;
            // Jx^(d,p)
            ( *Jx2D )( idloc, jploc ) += Jx_ion * Sxd[i]*Syp[j];
            // Jy^(p,d)
            ( *Jy2D )( iploc, jdloc ) += Jy_ion * Sxp[i]*Syd[j];
            // Jz^(p,p)
            ( *Jz2D )( iploc, jploc ) += Jz_ion * Sxp[i]*Syp[j];
        }
    }//i*/


} // END Project global current densities (ionize)


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector vectorized
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1OrderRuytenV::currents( ElectroMagnAM *emAM,
                                   Particles &particles,
                                   unsigned int istart,
                                   unsigned int iend,
                                   double * __restrict__ invgf,
                                   int * __restrict__ iold,
                                   double * __restrict__ deltaold,
                                   std::complex<double> * __restrict__ array_eitheta_old,
                                   int npart_total,
                                   int ipart_ref,
                                   int ispec )
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int ipo = iold[0];
    int jpo = iold[1];
    int ipom2 = ipo-2;
    int jpom2 = jpo-2;

    int vecSize = 8;
    int bsize = 5*5*vecSize*Nmode_;

    std::complex<double> bJl[bsize] __attribute__( ( aligned( 64 ) ) );
    std::complex<double> bJr[bsize] __attribute__( ( aligned( 64 ) ) );
    std::complex<double> bJt[bsize] __attribute__( ( aligned( 64 ) ) );

    double Sl0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double Sr0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double DSl[40] __attribute__( ( aligned( 64 ) ) );
    double DSr[40] __attribute__( ( aligned( 64 ) ) );
    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    double r_bar[8] __attribute__( ( aligned( 64 ) ) );
    complex<double> * __restrict__ Jl;
    complex<double> * __restrict__ Jr;
    complex<double> * __restrict__ Jt;

    double *invR_local = &(invR_[jpom2]);
    double *invRd_local = &(invRd_[jpom2]);

    // Pointer for GPU and vectorization on ARM processors
    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ position_z = particles.getPtrPosition(2);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);
    double * __restrict__ weight     = particles.getPtrWeight();
    short  * __restrict__ charge     = particles.getPtrCharge();
    int * __restrict__ cell_keys  = particles.getPtrCellKeys();

    #pragma omp simd
    for( unsigned int j=0; j<200*Nmode_; j++ ) {
        bJl[j] = 0.;
        bJr[j] = 0.;
        bJt[j] = 0.;
    }

    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );

    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed = min( cell_nparts-ivect, vecSize );
        int istart0 = ( int )istart + ivect;
        complex<double> e_bar[8], e_delta_m1[8];

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            compute_distances( position_x, position_y, position_z, cell_keys, npart_total, ipart, istart0, ipart_ref, deltaold, array_eitheta_old, iold, Sl0_buff_vect, Sr0_buff_vect, DSl, DSr, r_bar, e_bar, e_delta_m1 );
            charge_weight[ipart] = inv_cell_volume * ( double )( charge[istart0+ipart] )*weight[istart0+ipart];
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeJl( ipart, charge_weight, DSl, DSr, Sr0_buff_vect, bJl, dl_ov_dt_, invR_local, e_bar);
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeJr( ipart, charge_weight, DSl, DSr, Sl0_buff_vect, bJr, one_ov_dt, invRd_local, e_bar, jpom2);
        }

        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeJt( ipart, &momentum_y[istart0], &momentum_z[istart0], charge_weight, &invgf[istart0-ipart_ref], DSl, DSr, Sl0_buff_vect, Sr0_buff_vect, bJt, invR_local, r_bar, e_bar, e_delta_m1, one_ov_dt);
        }
    } //End ivect

    int iloc0 = ipom2*nprimr_+jpom2;

    for( unsigned int imode=0; imode<( unsigned int )Nmode_; imode++ ) {
        if (ispec == 0) { // When diags are not needed ispec is always 0.
            Jl =  &( *emAM->Jl_[imode] )( 0 );
        } else { // When diags are needed ispec+1 is passed.
            unsigned int n_species = emAM->Jl_s.size() / Nmode_;
            unsigned int ifield = imode*n_species+ispec-1;
            Jl  = emAM->Jl_s    [ifield] ? &( * ( emAM->Jl_s    [ifield] ) )( 0 ) : &( *emAM->Jl_    [imode] )( 0 ) ;
        }
        int iloc = iloc0;
        for( unsigned int i=1 ; i<5 ; i++ ) {
            iloc += nprimr_;
            #pragma omp simd
            for( unsigned int j=0 ; j<5 ; j++ ) {
                complex<double> tmpJl( 0. );
                int ilocal = ( i*5+j )*vecSize;
                UNROLL(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpJl += bJl [200*imode + ilocal+ipart];
                }
                Jl[iloc+j] += tmpJl;
            }
        }
    }


    for( unsigned int imode=0; imode<( unsigned int )Nmode_; imode++ ) {
        if (ispec == 0) { // When diags are not needed ispec is always 0.
            Jr =  &( *emAM->Jr_[imode] )( 0 );
        } else { // When diags are needed ispec+1 is passed.
            unsigned int n_species = emAM->Jr_s.size() / Nmode_;
            unsigned int ifield = imode*n_species+ispec-1;
            Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : &( *emAM->Jr_    [imode] )( 0 ) ;
        }

        int iloc = iloc0 + ipom2 + 1;
        for( unsigned int i=0 ; i<5 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=0 ; j<4 ; j++ ) {
                complex<double> tmpJr( 0. );
                int ilocal = ( i*5+j+1 )*vecSize;
                UNROLL(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpJr += bJr [200*imode + ilocal+ipart];
                }
                Jr[iloc+j] += tmpJr;
            }
            iloc += nprimr_+1;
        }
    }

    for( unsigned int imode=0; imode<( unsigned int )Nmode_; imode++ ) {
        if (ispec == 0) { // When diags are not needed ispec is always 0.
            Jt =  &( *emAM->Jt_[imode] )( 0 );
        } else { // When diags are needed ispec+1 is passed.
            unsigned int n_species = emAM->Jt_s.size() / Nmode_;
            unsigned int ifield = imode*n_species+ispec-1;
            Jt  = emAM->Jt_s    [ifield] ? &( * ( emAM->Jt_s    [ifield] ) )( 0 ) : &( *emAM->Jt_    [imode] )( 0 ) ;
        }
        int iloc = iloc0;
        for( unsigned int i=0 ; i<5 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=0 ; j<5 ; j++ ) {
                complex<double> tmpJt( 0. );
                int ilocal = ( i*5+j )*vecSize;
                UNROLL(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpJt += bJt [200*imode + ilocal+ipart];
                }
                Jt[iloc+j] += tmpJt;
            }
            iloc += nprimr_;
        }
    }
} // END Projection currents vectorized


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1OrderRuytenV::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread,  bool diag_flag, bool is_spectral, int ispec, int scell, int ipart_ref )
{
    if( istart == iend ) {
        return;    //Don't treat empty cells.
    }

    //Independent of cell. Should not be here
    //{
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    std::vector<std::complex<double>> *array_eitheta_old = &( smpi->dynamics_eithetaold[ithread] );
    ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );


    //}
    int iold[2];
    iold[0] = scell/nscellr_+oversize_[0];
    iold[1] = ( scell%nscellr_ )+oversize_[1];


    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if( !diag_flag ) {
        if( !is_spectral ) {
            currents( emAM, particles,  istart, iend, invgf->data(), iold, delta->data(), array_eitheta_old->data(), invgf->size(), ipart_ref );
        } else {
            ERROR( "Vectorized projection is not supported in spectral AM" );
        }

        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        //double *b_Jx  = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
        //double *b_Jy  = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
        //double *b_Jz  = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;
        //double *b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
        currentsAndDensity( emAM, particles, istart, iend, invgf->data(), iold, delta->data(), array_eitheta_old->data(), invgf->size(), ipart_ref, ispec );
    }
}

// Project susceptibility
void ProjectorAM1OrderRuytenV::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )
{
    double dts2 = dt/2.;
    double dts4 = dt/4.;
    double * __restrict__ Chi_envelope = &( *EMfields->Env_Chi_ )( 0 ) ;

    int iold[2];
    iold[0] = icell/nscellr_+oversize_[0];
    iold[1] = ( icell%nscellr_ )+oversize_[1];

    int ipom2 = iold[0]-2;
    int jpom2 = iold[1]-2;

    int vecSize = 8;
    int bsize = 5*5*vecSize; // Chi has only one mode //*Nmode_;

    std::vector<double> *Epart       = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Phipart     = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPhipart = &( smpi->dynamics_GradPHIpart[ithread] );
    double * __restrict__ inv_gamma_ponderomotive = &( smpi->dynamics_inv_gamma_ponderomotive[ithread][0] );


    int nparts = smpi->dynamics_invgf[ithread].size();
    double * __restrict__ Ex       = &( ( *Epart )[0*nparts] );
    double * __restrict__ Ey       = &( ( *Epart )[1*nparts] );
    double * __restrict__ Ez       = &( ( *Epart )[2*nparts] );
    double * __restrict__ Phi      = &( ( *Phipart )[0*nparts] );
    double * __restrict__ GradPhix = &( ( *GradPhipart )[0*nparts] );
    double * __restrict__ GradPhiy = &( ( *GradPhipart )[1*nparts] );
    double * __restrict__ GradPhiz = &( ( *GradPhipart )[2*nparts] );

    double bChi[bsize] __attribute__( ( aligned( 64 ) ) );

    double Sl1[32] __attribute__( ( aligned( 64 ) ) );
    double Sr1[32] __attribute__( ( aligned( 64 ) ) );
    // double Sl0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    // double Sr0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    // double DSl[40] __attribute__( ( aligned( 64 ) ) );
    // double DSr[40] __attribute__( ( aligned( 64 ) ) );
    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    // double r_bar[8] __attribute__( ( aligned( 64 ) ) );

    // Pointer for GPU and vectorization on ARM processors
    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = particles.getPtrPosition(1);
    double * __restrict__ position_z = particles.getPtrPosition(2);
    double * __restrict__ momentum_x = particles.getPtrMomentum(0);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);
    double * __restrict__ weight     = particles.getPtrWeight();
    short  * __restrict__ charge     = particles.getPtrCharge();
    //int * __restrict__ cell_keys  = particles.getPtrCellKeys();

    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );

    double one_over_mass=1./species_mass;

    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {

        int np_computed = min( cell_nparts-ivect, vecSize );
        int istart0 = ( int )istart + ivect;

        #pragma omp simd
        for( unsigned int j=0; j<200*Nmode_; j++ ) {
            bChi[j] = 0.;
        }

       #pragma omp simd
       for( int ipart=0 ; ipart<np_computed; ipart++ ) {

           double gamma_ponderomotive, gamma0, gamma0_sq;
           double charge_over_mass_dts2, charge_sq_over_mass_sq_dts4, charge_sq_over_mass_sq;
           double pxsm, pysm, pzsm;

           int ipart2 = istart-ipart_ref+ipart;

           double c = charge[istart0+ipart];

           charge_over_mass_dts2       = c *dts2*one_over_mass;
           // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
           charge_sq_over_mass_sq_dts4 = c*c*dts4*one_over_mass*one_over_mass;
           // (charge over mass)^2
           charge_sq_over_mass_sq      = c*c*one_over_mass*one_over_mass;

           // compute initial ponderomotive gamma
           gamma0_sq = (1. + momentum_x[istart0+ipart]*momentum_x[istart0+ipart]
                     + momentum_y[istart0+ipart]*momentum_y[istart0+ipart]
                     + momentum_z[istart0+ipart]*momentum_z[istart0+ipart]
                     + Phi[istart0-ipart_ref+ipart]*charge_sq_over_mass_sq);

           gamma0    = sqrt( gamma0_sq ) ;

           // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
           pxsm = ( gamma0 * charge_over_mass_dts2*( Ex[ipart2] ) - charge_sq_over_mass_sq_dts4*(  GradPhix[ipart2] ) ) * momentum_x[istart0+ipart] / gamma0_sq;
           pysm = ( gamma0 * charge_over_mass_dts2*( Ey[ipart2] ) - charge_sq_over_mass_sq_dts4*(  GradPhiy[ipart2] ) ) * momentum_y[istart0+ipart] / gamma0_sq;
           pzsm = ( gamma0 * charge_over_mass_dts2*( Ez[ipart2] ) - charge_sq_over_mass_sq_dts4*(  GradPhiz[ipart2] ) ) * momentum_z[istart0+ipart] / gamma0_sq;

           // update of gamma ponderomotive
           gamma_ponderomotive = gamma0 + ( pxsm+pysm+pzsm )*0.5 ;
           // buffer inverse of ponderomotive gamma to use it in ponderomotive momentum pusher
           inv_gamma_ponderomotive[istart0 + ipart - ipart_ref] = 1./gamma_ponderomotive;

           // susceptibility for the macro-particle
           charge_weight[ipart-ipart_ref] = c*c*inv_cell_volume * weight[istart0+ipart-ipart_ref]*one_over_mass*inv_gamma_ponderomotive[istart0 + ipart - ipart_ref] ;

           // variable declaration
           double xpn, rpn;
           double delta, delta2;

           // Initialize all current-related arrays to zero
           Sl1[ipart] = 0.;
           Sr1[ipart] = 0.;

           Sl1[vecSize+ipart] = 0.;
           Sr1[vecSize+ipart] = 0.;

           Sl1[2*vecSize+ipart] = 0.;
           Sr1[2*vecSize+ipart] = 0.;

           // --------------------------------------------------------
           // Locate particles & Calculate Esirkepov coef. S, DS and W
           // --------------------------------------------------------

           // locate the particle on the primal grid at current time-step & calculate coeff. S1
          xpn = position_x[istart0+ipart-ipart_ref] * dl_inv_;
          int ip = round( xpn );
          delta  = xpn - ( double )ip;
          delta2 = delta*delta;
          Sl1[0*vecSize+ipart] = 0.5 * ( delta2-delta+0.25 );
          Sl1[1*vecSize+ipart] = 0.75-delta2;
          Sl1[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );

          rpn = position_y[istart0+ipart-ipart_ref] * position_y[istart0+ipart-ipart_ref] ;
          rpn += position_z[istart0+ipart-ipart_ref] * position_z[istart0+ipart-ipart_ref] ;
          rpn = sqrt(rpn) * dr_inv_ ;
          int jp = round( rpn );
          delta  = rpn - ( double )jp;

          double x_n = floor(rpn);
          double coeff0 = (x_n+1-rpn)*(5*x_n + 2 - rpn)/(4.*x_n + 2.);

          jp -= j_domain_begin_ + 2;
          Sr1[0*vecSize+ipart] = coeff0 * (double)(delta < 0.);
          Sr1[2*vecSize+ipart] = (1.-coeff0) * (double)(delta >= 0.);
          Sr1[1*vecSize+ipart] = 1. - Sr1[0*vecSize+ipart] - Sr1[2*vecSize+ipart];

          Sr1[0*vecSize+ipart] *=  invR_[1+jp];
          Sr1[1*vecSize+ipart] *=  invR_[2+jp];
          Sr1[2*vecSize+ipart] *=  invR_[3+jp];

      } // end ipart loop

      #pragma omp simd
      for( int ipart=0 ; ipart<np_computed; ipart++ ) {
          UNROLL_S(5)
          for( unsigned int i=0 ; i<3 ; i++ ) {
              UNROLL_S(5)
              for ( unsigned int j=0; j<3 ; j++ ) {
                  int index( ( i*5 + j )*vecSize+ipart );//cout <<ipart<<" "<<i<<" "<<j<<endl;
                  bChi [index] += charge_weight[ipart]* Sl1[i*vecSize+ipart]*Sr1[j*vecSize+ipart] ;
              }
          }
      } // end ipart loop

      // ---------------------------
      // Compute the total charge
      // ---------------------------
      int iloc0 = ipom2*nprimr_+jpom2;
      int iloc = iloc0;
      for( unsigned int i=0 ; i<5 ; i++ ) {
          #pragma omp simd
          for( unsigned int j=0 ; j<5 ; j++ ) {
              double tmpChi( 0. );
              int ilocal = ( i*5+j )*vecSize;
              UNROLL(8)
              for( int ipart=0 ; ipart<8; ipart++ ) {
                  tmpChi += bChi [ilocal+ipart];
              }
              Chi_envelope[iloc+j] += tmpChi;
          }
          iloc += nprimr_;
      }

    } // end ivect

} // end ProjectorAM1OrderRuytenV::susceptibility
