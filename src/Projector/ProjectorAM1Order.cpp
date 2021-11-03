#include "ProjectorAM1Order.h"

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

//ProjectorAM1Order is written for exclusive use of spectral solvers with non staggered grids, primal size grids with half cell length shift along r.

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ProjectorAM1Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM1Order::ProjectorAM1Order( Params &params, Patch *patch ) : ProjectorAM( params, patch )
{
    dt = params.timestep;
    dr = params.cell_length[1];
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt_  = params.cell_length[0] / params.timestep;
    dr_ov_dt_  = params.cell_length[1] / params.timestep;
    dr_inv_   = 1.0 / dr;
    one_ov_dt  = 1.0 / params.timestep;
    Nmode=params.nmodes;
    i_domain_begin_ = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin_ = patch->getCellStartingGlobalIndex( 1 );

    oversize_[1] = params.oversize[1]; 
    nprimr_ = params.n_space[1] + 2*params.oversize[1] + 1;
    npriml_ = params.n_space[0] + 2*params.oversize[0] + 1;

    invR = &((static_cast<PatchAM *>( patch )->invR)[0]);
    invRd = &((static_cast<PatchAM *>( patch )->invRd)[0]);

    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorAM1Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorAM1Order::~ProjectorAM1Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Non charge conserving projector for diags at t=0, frozen species or compute charge 
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::basicForComplex( complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode )
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
   
    // This projection for currents is used only in cases where position=position_old.
    // Warning: will fail evaluating the current at t=0 if a plasma is already in the box.
    int iloc, nr( nprimr_ );
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    double m1powerimode = 1.;
    
    if( type > 0 ) { //if current density
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
        if( type == 1 ) { //if Jl
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) { //if Jr
            m1powerimode = -1.;
            charge_weight *= ( particles.momentum( 1, ipart )*particles.position( 1, ipart ) + particles.momentum( 2, ipart )*particles.position( 2, ipart ) ) / r ;
        } else { //if Jt
            m1powerimode = -1.;
            charge_weight *= ( -particles.momentum( 1, ipart )*particles.position( 2, ipart ) + particles.momentum( 2, ipart )*particles.position( 1, ipart ) ) / r ;
        }
    }
    
    complex<double> C_m = 1.;
    if( imode > 0 ) {
        double theta = atan2( particles.position( 2, ipart ) , particles.position( 1, ipart ) );// theta at t = t0
        complex<double> e_theta = std::polar( 1.0, theta );
        C_m = 2.;
        for( unsigned int i=0; i<( unsigned int )imode; i++ ) {
            C_m *= e_theta;
            m1powerimode *= -1.;
        }
    }
    
    double xpn, rpn;
    double delta;
    double Sl1[2], Sr1[2];
    
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = (particles.position( 0, ipart ) ) * dl_inv_ ;
    int ip = floor( xpn );
    delta  = xpn - ( double )ip;
    Sl1[0] = 1. - delta ;
    Sl1[1] = delta;

    rpn = r * dr_inv_ - 0.5 ; //-0.5 because of cells being shifted by dr/2
    int jp = floor( rpn );
    delta  = rpn - ( double )jp;
    Sr1[0] = 1. - delta;
    Sr1[1] = delta;
    
    if (rpn < 0.){ // If particle is between 0 and dr/2.
        jp = 0;
        Sr1[0] = Sr1[1] + m1powerimode*Sr1[0];
        Sr1[1] = 0.; 
    }
 
    ip -= i_domain_begin_ ;
    jp -= j_domain_begin_ ;
    // j_domain_begin_ should always be zero in spectral since no paralellization along r.
    
    for( unsigned int i=0 ; i<2 ; i++ ) {
        iloc = ( i+ip )*nr+jp;
        for( unsigned int j=0 ; j<2 ; j++ ) {
            rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR[j+jp];
        }
    }//i
} // END Project for diags local current densities

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents and charge densities for all modes, not charge conserving
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::currents( ElectroMagnAM *emAM, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, std::complex<double> *array_eitheta_old, bool diag_flag, int ispec)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    int nparts= particles.size();
    int ip[2], jp[2];
    int iloc[2], linindex;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double m1powerimode = 1.;
    
    // variable declaration
    double xpn[2], rp, rpn[2], delta;
    //double xpn_rho, rp_rho, rpn_rho; //Rho is not computed at the same particle position as J.

    double  Sl1[2][2], Sr1[2][3];
    complex<double> e_theta[2] = { 1., 1. };
    complex<double> C_m[2] = { 1., 1. };
    complex<double> *Jl, *Jr, *Jt, *rho;
    
    std::complex<double> theta_old = array_eitheta_old[0]; // theta at t = t0 - dt
    rp = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    std::complex<double> eitheta = ( particles.position( 1, ipart ) + Icpx * particles.position( 2, ipart ) ) / rp ; //exp(i theta)

    std::complex<double> e_delta_m1 = std::sqrt(eitheta * (2.*std::real(theta_old) - theta_old)); // std::sqrt keeps the root with positive real part which is what we need here.

    theta_old += e_delta_m1; // eitheta at t = t0 - dt/2
    e_theta[0] = theta_old;
    e_theta[1] = eitheta;

    double crl_p =  ( particles.momentum( 0, ipart )) *invgf;
    double crt_p =  ( particles.momentum( 2, ipart )*real(e_theta[0]) - particles.momentum( 1, ipart )*imag(e_theta[0]) ) * invgf;
    double crr_p =  ( particles.momentum( 1, ipart )*real(e_theta[0]) + particles.momentum( 2, ipart )*imag(e_theta[0]) ) * invgf;

    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn[0] = i_domain_begin_ + iold[0*nparts] + deltaold[0*nparts];  // position at t=t0
    xpn[1] = particles.position( 0, ipart ) * dl_inv_ ;             // position at t=t0+dt
    xpn[0] = 0.5*(xpn[0]+xpn[1]);                                   // position at t=t0+dt/2
    rpn[0] = j_domain_begin_ + iold[1*nparts] + deltaold[1*nparts];
    rpn[1] = rp * dr_inv_ - 0.5 ;
    rpn[0] = 0.5*(rpn[0]+rpn[1]);

    for (int irho=0; irho < 2; irho++){ //irho=0 for currents, irho=1 for charge density

        ip[irho] = floor( xpn[irho] );
        delta  = xpn[irho] - ( double )ip[irho];
        Sl1[irho][0] = 1. - delta;
        Sl1[irho][1] = delta;
        
        jp[irho] = floor( rpn[irho] );
        delta  = rpn[irho] - ( double )jp[irho];
        Sr1[irho][0] = 1. - delta;
        Sr1[irho][1] = delta;

        if (rpn[irho] < 0.){ // If particle is between 0 and dr/2.
            jp[irho] = 0;
            Sr1[irho][2] = Sr1[irho][0]; // This part deposited "below" axis must be brought back with a sign depending on the mode and on the quantity
            Sr1[irho][0] = Sr1[irho][1];
            Sr1[irho][1] = 0.; 
        } else {
            Sr1[irho][2] = 0.;
        }

        ip[irho]  -= i_domain_begin_ ;
        jp[irho]  -= j_domain_begin_ ;

    }

    for( unsigned int imode=0; imode<( unsigned int )Nmode; imode++ ) {
        if( imode == 1 ) {
            C_m[0] = 2.;
            C_m[1] = 2.;
        }
        if( imode > 0 ) {
            C_m[0] *= e_theta[0];
            C_m[1] *= e_theta[1];
            m1powerimode *= -1.;
        }
        
       unsigned int n_species = emAM->Jl_s.size() / Nmode;
       unsigned int ifield = imode*n_species+ispec;
        if (!diag_flag){
            Jl =  &( *emAM->Jl_[imode] )( 0 );
            Jr =  &( *emAM->Jr_[imode] )( 0 );
            Jt =  &( *emAM->Jt_[imode] )( 0 );
            rho = &( *emAM->rho_AM_[imode] )( 0 ) ; // In spectral, always project density
        } else {
            //Jl  = emAM->Jl_s    [ifield] ? &( * ( emAM->Jl_s    [ifield] ) )( 0 ) : &( *emAM->Jl_    [imode] )( 0 ) ;
            //Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : &( *emAM->Jr_    [imode] )( 0 ) ;
            //Jt  = emAM->Jt_s    [ifield] ? &( * ( emAM->Jt_s    [ifield] ) )( 0 ) : &( *emAM->Jt_    [imode] )( 0 ) ;
            //rho = emAM->rho_AM_s[ifield] ? &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
            Jl  =  &( * ( emAM->Jl_s    [ifield] ) )( 0 ) ;
            Jr  =  &( * ( emAM->Jr_s    [ifield] ) )( 0 ) ;
            Jt  =  &( * ( emAM->Jt_s    [ifield] ) )( 0 ) ;
            rho =  &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) ;

        }

        // Jr^(p,p) Jt^(p,p) Jl^(p,p) Rho^(p,p)
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc[0] = ( i+ip[0] )* nprimr_ + jp[0];
            iloc[1] = ( i+ip[1] )* nprimr_ + jp[1];
            for( unsigned int j=0 ; j<2 ; j++ ) {
                linindex = iloc[0]+j;
                complex<double> increment =  C_m[0]*charge_weight* Sl1[0][i]*Sr1[0][j]*invR[jp[0]+j];
                Jl [linindex] += crl_p * increment ;
                Jr [linindex] += crr_p * increment ;
                Jt [linindex] += crt_p * increment ;
                rho [iloc[1]+j] += C_m[1]*charge_weight* Sl1[1][i]*Sr1[1][j]*invR[jp[1]+j];
            }
        }//i

        //Correction below axis for currents
        if (Sr1[0][2] != 0.) {
            for( unsigned int i=0 ; i<2 ; i++ ) {
                iloc[0] = ( i+ip[0] )* nprimr_ + jp[0];
                linindex = iloc[0];
                complex<double> increment =  C_m[0]*charge_weight* Sl1[0][i]*m1powerimode*Sr1[0][2]*invR[jp[0]];
                Jl [linindex] += crl_p * increment ;
                Jr [linindex] -= crr_p * increment ;
                Jt [linindex] -= crt_p * increment ;
            }
        }
        //Correction below axis for density
        if (Sr1[1][2] != 0.) {
            for( unsigned int i=0 ; i<2 ; i++ ) {
                iloc[1] = ( i+ip[1] )* nprimr_ + jp[1];
                rho [iloc[1]] += C_m[1]*charge_weight* Sl1[1][i]*m1powerimode*Sr1[1][2]*invR[jp[1]];
            }
        }

    }// end loop on modes
    
} // END Project local current and charge densities (rho, Jl, Jr, Jt)

void ProjectorAM1Order::axisBC(ElectroMagnAM *emAM, bool diag_flag )
{
    
    for (unsigned int imode=0; imode < Nmode; imode++){ 

       std::complex<double> *rho     = &( *emAM->rho_AM_[imode] )( 0 );
       std::complex<double> *rho_old = &( *emAM->rho_old_AM_[imode] )( 0 );
       std::complex<double> *jl      = &( *emAM->Jl_[imode] )( 0 );
       std::complex<double> *jr      = &( *emAM->Jr_[imode] )( 0 );
       std::complex<double> *jt      = &( *emAM->Jt_[imode] )( 0 );
       apply_axisBC(rho, imode, 0);
       apply_axisBC(rho_old, imode, 0);
       apply_axisBC(jl, imode, 0);
       apply_axisBC(jr, imode, 1);
       apply_axisBC(jt, imode, 1);
    }

    if (diag_flag){
        unsigned int n_species = emAM->Jl_s.size() / Nmode;
        for( unsigned int imode = 0 ; imode < emAM->Jl_.size() ; imode++ ) {
            for( unsigned int ispec = 0 ; ispec < n_species ; ispec++ ) {
                unsigned int ifield = imode*n_species+ispec;
                complex<double> *rho = emAM->rho_AM_s[ifield] ? &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) : NULL ;
                complex<double> *jl = emAM->Jl_s[ifield] ? &( * ( emAM->Jl_s[ifield] ) )( 0 ) : NULL ;
                complex<double> *jr = emAM->Jr_s[ifield] ? &( * ( emAM->Jt_s[ifield] ) )( 0 ) : NULL ;
                complex<double> *jt = emAM->Jt_s[ifield] ? &( * ( emAM->Jr_s[ifield] ) )( 0 ) : NULL ;
                apply_axisBC( rho, imode, 0);
                apply_axisBC( jl, imode, 0);
                apply_axisBC( jr, imode, 1);
                apply_axisBC( jt, imode, 1);
            }
        }
    }

}

void ProjectorAM1Order::apply_axisBC(std::complex<double> *rho, unsigned int imode, unsigned int nonzeromode)
{
    //If pointer is NULL, nothing to do
    if(!rho) return;

    const double one_ov_9  = 1./9.; 
    const double one_ov_16 = 1./16.; 
    //non zero mode is not zero on axis. Mode 0 for rho and Jl, mode 1 for Jr and Jt
    if (imode == nonzeromode){
        for( unsigned int i=oversize_[1] ; i<npriml_*nprimr_+oversize_[1]; i+=nprimr_ ) {
            rho[i] = (25.*rho[i+1] - 9.*rho[i+2])*one_ov_16;
        }//i
    } else { //m !=  non zero mode
        // rho_m[r=0] = 0 and drho_m/dr[r=0] = 0 when m is even and !=0 when m is odd. 
        // quadratic interpolation for even m and linear interpolation when m is odd
        const double slope = (imode%2==0 ? one_ov_9 : 1./3.);
        for( unsigned int i=oversize_[1] ; i<npriml_*nprimr_+oversize_[1]; i+=nprimr_ ) {
            rho[i] = rho[i+1] * slope;
        }//i
    }

    return;
}

//------------------------------------//
//Wrapper for projection
void ProjectorAM1Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
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

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization NOT DONE YET
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion )
{
}


// Projector for susceptibility used as source term in envelope equation
void ProjectorAM1Order::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref )

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
    
    
        double xpn, rpn;
        double delta;
        double Sl1[2], Sr1[2];
    
    
        // locate the particle on the primal and dual grid at current time-step & calculate coeff. S1
        xpn = particles.position( 0, ipart ) * dl_inv_;
        int ip = int( xpn );
        delta  = xpn - ( double )ip;
        Sl1[0] = delta ;
        Sl1[1] = 1.-delta;
        rpn = r * dr_inv_ ;
        int jp = int( rpn );
        delta  = rpn - ( double )jp;
        Sr1[0] = delta;
        Sr1[1] = 1.-delta;
    
        ip -= i_domain_begin_ ;
        jp -= j_domain_begin_ ;
    
    
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )*nr+jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                    Chi_envelope [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR[j+jp];
            }
        }//i
    


    }
    
}



void ProjectorAM1Order::axisBCEnvChi( double *EnvChi )
{
    if(EnvChi == NULL)
        return;

    //const double one_ov_9  = 1./9.; 
    const double one_ov_16 = 1./16.; 
    
    for( unsigned int i=oversize_[1] ; i<npriml_*nprimr_+oversize_[1]; i+=nprimr_ ) {
          
        EnvChi[i] = (25.*EnvChi[i+1] - 9.*EnvChi[i+2])*one_ov_16;
    }//i
    
return;
}

