#include "ProjectorAM2Order.h"

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
// Constructor for ProjectorAM2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM2Order::ProjectorAM2Order( Params &params, Patch *patch ) : ProjectorAM( params, patch )
{
    dt = params.timestep;
    dr = params.cell_length[1];
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt_  = params.cell_length[0] / params.timestep;
    dr_ov_dt_  = params.cell_length[1] / params.timestep;
    dr_inv_   = 1.0 / dr;
    one_ov_dt  = 1.0 / params.timestep;
    Nmode_=params.nmodes;
    i_domain_begin_ = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin_ = patch->getCellStartingGlobalIndex( 1 );

    nprimr_ = params.patch_size_[1] + 2*params.oversize[1] + 1;
    npriml_ = params.patch_size_[0] + 2*params.oversize[0] + 1;

    invR_ = &((static_cast<PatchAM *>( patch )->invR)[0]);
    invRd_ = &((static_cast<PatchAM *>( patch )->invRd)[0]);

    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorAM2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM2Order::~ProjectorAM2Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for all modes
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::currents(   ElectroMagnAM *emAM,
                                    Particles &particles,
                                    unsigned int ipart,
                                    double invgf,
                                    int *iold,
                                    double *deltaold,
                                    std::complex<double> *array_eitheta_old,
                                    bool diag_flag, int ispec)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------   int iloc,
    int nparts= particles.size();
    int iloc, jloc, linindex;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double crl_p = charge_weight*dl_ov_dt_;
    double crr_p = charge_weight*one_ov_dt;

    // variable declaration
    double xpn, ypn;
    double delta, delta2;
    // arrays used for the Esirkepov projection method
    double  Sl0[5], Sl1[5], Sr0[5], Sr1[5], DSl[5], DSr[5];
    complex<double>  Jl_p[5], Jr_p[5];
    complex<double> e_delta, e_delta_m1, e_bar, e_bar_m1, C_m = 1.; //, C_m_old;
    complex<double> *Jl, *Jr, *Jt, *rho;

    for( unsigned int i=0; i<5; i++ ) {
        Sl1[i] = 0.;
        Sr1[i] = 0.;
    }
    Sl0[0] = 0.;
    Sl0[4] = 0.;
    Sr0[0] = 0.;
    Sr0[4] = 0.;
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at former time-step & calculate coeff. S0
    delta = deltaold[0*nparts];
    delta2 = delta*delta;
    Sl0[1] = 0.5 * ( delta2-delta+0.25 );
    Sl0[2] = 0.75-delta2;
    Sl0[3] = 0.5 * ( delta2+delta+0.25 );

    delta = deltaold[1*nparts];
    delta2 = delta*delta;
    Sr0[1] = 0.5 * ( delta2-delta+0.25 );
    Sr0[2] = 0.75-delta2;
    Sr0[3] = 0.5 * ( delta2+delta+0.25 );
    //calculate exponential coefficients

    double rp = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    std::complex<double> eitheta_old = array_eitheta_old[0];
    std::complex<double> eitheta = ( particles.position( 1, ipart ) + Icpx * particles.position( 2, ipart ) ) / rp ; //exp(i theta)
    e_bar = 1.;
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dl_inv_;
    int ip = round( xpn );
    int ipo = iold[0*nparts];
    int ip_m_ipo = ip-ipo-i_domain_begin_;
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sl1[ip_m_ipo+1] = 0.5 * ( delta2-delta+0.25 );
    Sl1[ip_m_ipo+2] = 0.75-delta2;
    Sl1[ip_m_ipo+3] = 0.5 * ( delta2+delta+0.25 );

    ypn = rp *dr_inv_ ;
    int jp = round( ypn );
    int jpo = iold[1*nparts];
    int jp_m_jpo = jp-jpo-j_domain_begin_;
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sr1[jp_m_jpo+1] = 0.5 * ( delta2-delta+0.25 );
    Sr1[jp_m_jpo+2] = 0.75-delta2;
    Sr1[jp_m_jpo+3] = 0.5 * ( delta2+delta+0.25 );

    for( unsigned int i=0; i < 5; i++ ) {
        DSl[i] = Sl1[i] - Sl0[i];
        DSr[i] = Sr1[i] - Sr0[i];
    }

    double r_bar = ((jpo + j_domain_begin_ + deltaold[1*nparts])*dr + rp) * 0.5; // r at t = t0 - dt/2

    e_delta_m1 = std::sqrt(eitheta * std::conj(eitheta_old)); //ei^(theta-theta_old)/2. std::sqrt keeps the root with positive real part which is what we need here.
    e_bar_m1 = eitheta_old * e_delta_m1;                      //ei^(theta+theta_old)/2.

    ipo -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
    // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
    jpo -= 2;

    double *invR__local = &(invR_[jpo]);

    // ------------------------------------------------
    // Local current created by the particle
    // calculate using the charge conservation equation
    // ------------------------------------------------

    // ---------------------------
    // Calculate the total current
    // ---------------------------

    //initial value of crt_p for imode = 0.
    complex<double> crt_p= charge_weight*( particles.momentum( 2, ipart )* real(e_bar_m1) - particles.momentum( 1, ipart )*imag(e_bar_m1) ) * invgf;

    // Compute everything independent of theta
    double tmpJl[5];
    for( unsigned int j=0 ; j<5 ; j++ ) {
        tmpJl[j] = crl_p * ( Sr0[j] + 0.5*DSr[j] )* invR__local[j];
    }
    Jl_p[0]= 0.;
    for( unsigned int i=1 ; i<5 ; i++ ) {
        Jl_p[i]= Jl_p[i-1] - DSl[i-1];
    }


    double Vd[5];
    double tmpJr[5];
    for( int j=3 ; j>=0 ; j-- ) {
        jloc = j+jpo+1;
        Vd[j] = abs( jloc + j_domain_begin_ + 0.5 )* invRd_[jloc]*dr ;
        tmpJr[j] = crr_p * DSr[j+1] * invRd_[jpo+j+1]*dr;
    }
    Jr_p[4]= 0.;
    for( int j=3 ; j>=0 ; j-- ) {
        Jr_p[j] =  Jr_p[j+1] * Vd[j] + tmpJr[j];
    }

    e_delta = 0.5;

    //Compute division by R in advance for Jt and rho evaluation.
    for( unsigned int j=0 ; j<5 ; j++ ) {
        Sr0[j] *= invR__local[j];
        Sr1[j] *= invR__local[j];
    }

    for( unsigned int imode=0; imode<( unsigned int )Nmode_; imode++ ) {

        if (imode > 0){
            e_delta *= e_delta_m1;
            e_bar *= e_bar_m1;
            C_m = 2. * e_bar ; //multiply modes > 0 by 2 and C_m = 1 otherwise.
            crt_p = charge_weight*Icpx*e_bar / ( dt*( double )imode )*2.*r_bar;
        }

        // Add contribution J_p to global array
        if (!diag_flag){
            Jl =  &( *emAM->Jl_[imode] )( 0 );
            Jr =  &( *emAM->Jr_[imode] )( 0 );
            Jt =  &( *emAM->Jt_[imode] )( 0 );
        } else {
            unsigned int n_species = emAM->Jl_s.size() / Nmode_;
            unsigned int ifield = imode*n_species+ispec;
            Jl  = emAM->Jl_s    [ifield] ? &( * ( emAM->Jl_s    [ifield] ) )( 0 ) : &( *emAM->Jl_    [imode] )( 0 ) ;
            Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : &( *emAM->Jr_    [imode] )( 0 ) ;
            Jt  = emAM->Jt_s    [ifield] ? &( * ( emAM->Jt_s    [ifield] ) )( 0 ) : &( *emAM->Jt_    [imode] )( 0 ) ;
            rho = emAM->rho_AM_s[ifield] ? &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;

            for( unsigned int i=0 ; i<5 ; i++ ) {
                iloc = ( i+ipo )*nprimr_;
                for( unsigned int j=0 ; j<5 ; j++ ) {
                    jloc = j+jpo;
                    linindex = iloc+jloc;
                    rho [linindex] += C_m*charge_weight* Sl1[i]*Sr1[j];
                }
            }//i
        }

        // Jl^(d,p)
        for( unsigned int i=1 ; i<5 ; i++ ) {
            iloc = ( i+ipo )*nprimr_+jpo;
            for( unsigned int j=0 ; j<5 ; j++ ) {
                linindex = iloc+j;
                Jl [linindex] += C_m * Jl_p[i]*tmpJl[j] ;
            }
        }//i

        // Jr^(p,d)
        for( unsigned int i=0 ; i<5 ; i++ ) {
            iloc = ( i+ipo )*( nprimr_+1 )+jpo+1;
            for( unsigned int j=0 ; j<4 ; j++ ) {
                linindex = iloc+j;
                Jr [linindex] += C_m * ( Sl0[i] + 0.5*DSl[i] ) * Jr_p[j] ;
            }
        }//i

        // Jt^(p,p)
        for( unsigned int i=0 ; i<5 ; i++ ) {
            iloc = ( i+ipo )*nprimr_ + jpo;
            for( unsigned int j=0 ; j<5 ; j++ ) {
                linindex = iloc+j;
                Jt [linindex] -= crt_p*(Sr1[j]*Sl1[i]*( e_delta-1. ) - Sr0[j]*Sl0[i]*(std::conj(e_delta) - 1.));
            }
        }

        if (imode == 0) e_delta = 1. ; //Restore e_delta correct initial value.
    }// end loop on modes

} // END Project local current densities (Jl, Jr, Jt, sort)

// ---------------------------------------------------------------------------------------------------------------------
//! Project for diags and frozen species -
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::basicForComplex( complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode )
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
void ProjectorAM2Order::axisBC(ElectroMagnAM *emAM, bool diag_flag )
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

void ProjectorAM2Order::apply_axisBC(std::complex<double> *rhoj,std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, unsigned int imode, bool diag_flag )
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

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM2Order::ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion )
{

    return;

    cField2D *JlAM  = static_cast<cField2D *>( Jl );
    cField2D *JrAM  = static_cast<cField2D *>( Jr );
    cField2D *JtAM  = static_cast<cField2D *>( Jt );


    //Declaration of local variables
    int ip, id, jp, jd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double Slp[3], Sld[3], Srp[3], Srd[3];

    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jl_ion = Jion.x * weight;
    double Jr_ion = Jion.y * weight;
    double Jt_ion = Jion.z * weight;

    //Locate particle on the grid
    xpn    = particles.position( 0, ipart ) * dl_inv_; // normalized distance to the first node
    ypn = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) )*dr_inv_ ;
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

    Slp[0] = 0.5 * ( xpmxip2-xpmxip+0.25 );
    Slp[1] = ( 0.75-xpmxip2 );
    Slp[2] = 0.5 * ( xpmxip2+xpmxip+0.25 );

    Sld[0] = 0.5 * ( xpmxid2-xpmxid+0.25 );
    Sld[1] = ( 0.75-xpmxid2 );
    Sld[2] = 0.5 * ( xpmxid2+xpmxid+0.25 );

    Srp[0] = 0.5 * ( ypmyjp2-ypmyjp+0.25 );
    Srp[1] = ( 0.75-ypmyjp2 );
    Srp[2] = 0.5 * ( ypmyjp2+ypmyjp+0.25 );

    Srd[0] = 0.5 * ( ypmyjd2-ypmyjd+0.25 );
    Srd[1] = ( 0.75-ypmyjd2 );
    Srd[2] = 0.5 * ( ypmyjd2+ypmyjd+0.25 );

    ip  -= i_domain_begin_;
    id  -= i_domain_begin_;
    jp  -= j_domain_begin_;
    jd  -= j_domain_begin_;

    for( unsigned int i=0 ; i<3 ; i++ ) {
        //int iploc=ip+i-1;
        int idloc=id+i-1;
        for( unsigned int j=0 ; j<3 ; j++ ) {
            int jploc=jp+j-1;
            //int jdloc=jd+j-1;
            if( jploc+ j_domain_begin_ ==0 ) {
                // Jl^(d,p)
                ( *JlAM )( idloc, jploc ) += Jl_ion*8. /dr * Sld[i]*Srp[j];
                ( *JrAM )( idloc, jploc ) += Jr_ion*8. /dr * Slp[i]*Srd[j];
                ( *JtAM )( idloc, jploc ) += Jt_ion*8. /dr * Slp[i]*Srp[j]; //A corriger dualite et repliement
            } else {
                ( *JlAM )( idloc, jploc ) += Jl_ion /( ( jploc+ j_domain_begin_ )*dr ) * Sld[i]*Srp[j];
                ( *JrAM )( idloc, jploc ) += Jr_ion /( ( jploc+ j_domain_begin_ )*dr ) * Slp[i]*Srd[j];
                ( *JtAM )( idloc, jploc ) += Jt_ion /( ( jploc+ j_domain_begin_ )*dr ) * Slp[i]*Srp[j];
            }

        }
    }//i


} // END Project global current densities (ionize)

//------------------------------------//
//Wrapper for projection
void ProjectorAM2Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool /*is_spectral*/, int ispec, int /*icell*/, int /*ipart_ref*/ )
{

    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    std::vector<std::complex<double>> *array_eitheta_old = &( smpi->dynamics_eithetaold[ithread] );
    ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );

    for( int ipart=istart ; ipart<iend; ipart++ ) {
        currents( emAM, particles,  ipart, ( *invgf )[ipart], &( *iold )[ipart], &( *delta )[ipart], &( *array_eitheta_old )[ipart], diag_flag, ispec);
    }
}

// Projector for susceptibility used as source term in envelope equation
void ProjectorAM2Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int /*icell*/, int /*ipart_ref*/ )

{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    double *Chi_envelope = &( *EMfields->Env_Chi_ )( 0 );

    std::vector<double> *Epart       = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Phipart     = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPhipart = &( smpi->dynamics_GradPHIpart[ithread] );
    std::vector<double> *inv_gamma_ponderomotive = &( smpi->dynamics_inv_gamma_ponderomotive[ithread] );


    double gamma_ponderomotive, gamma0, gamma0_sq;
    double charge_over_mass_dts2, charge_sq_over_mass_sq_dts4, charge_sq_over_mass_sq;
    double pxsm, pysm, pzsm;
    double one_over_mass=1./species_mass;
    double momentum[3];

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

        int iloc, nr( nprimr_ );

        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );

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

        //complex<double> e_theta = ( particles.position( 1, ipart ) + Icpx*particles.position( 2, ipart ) )/r;
        double C_m = 1.; // only mode 0


        double xpn, ypn;
        double delta, delta2;
        double Sl1[5], Sr1[5];


        // locate the particle on the primal grid at current time-step & calculate coeff. S1
        xpn = particles.position( 0, ipart ) * dl_inv_;
        int ip = round( xpn );
        delta  = xpn - ( double )ip;
        delta2 = delta*delta;
        Sl1[1] = 0.5 * ( delta2-delta+0.25 );
        Sl1[2] = 0.75-delta2;
        Sl1[3] = 0.5 * ( delta2+delta+0.25 );
        ypn = r * dr_inv_ ;
        int jp = round( ypn );
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

        for( unsigned int i=1 ; i<4 ; i++ ) {
            iloc = ( i+ip )*nr+jp;
            for( unsigned int j=1 ; j<4 ; j++ ) {
                    Chi_envelope [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR_[j+jp];
            }
        }//i



    }

}

void ProjectorAM2Order::axisBCEnvChi( double *EnvChi )
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

return;
}
