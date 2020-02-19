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
    dl_ov_dt  = params.cell_length[0] / params.timestep;
    dr_ov_dt  = params.cell_length[1] / params.timestep;
    dr_inv_   = 1.0 / dr;
    one_ov_dt  = 1.0 / params.timestep;
    one_nineth = 1./9.;
    Nmode=params.nmodes;
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );

    oversizeR = params.oversize[1]; 
    nprimr = params.n_space[1] + 2*params.oversize[1] + 1;
    npriml = params.n_space[0] + 2*params.oversize[0] + 1;

    invR = &((static_cast<PatchAM *>( patch )->invR)[0]);
    invRd = &((static_cast<PatchAM *>( patch )->invRd)[0]);

    if (!params.is_pxr)
        staggered = true;
    else
        staggered = false;
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
//! Non charge conserving projector for diags, frozen species or spectral solver - mode >= 0
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::basicForComplex( complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode )
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc, nr( nprimr );
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    
    if( type > 0 ) { //if current density
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
        if( type == 1 ) { //if Jl
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) { //if Jr
            charge_weight *= ( particles.momentum( 1, ipart )*particles.position( 1, ipart ) + particles.momentum( 2, ipart )*particles.position( 2, ipart ) ) / r ;
        } else { //if Jt
            charge_weight *= ( -particles.momentum( 1, ipart )*particles.position( 2, ipart ) + particles.momentum( 2, ipart )*particles.position( 1, ipart ) ) / r ;
        }
    }
    
    complex<double> C_m = 1.;
    if( imode > 0 ) {
        complex<double> e_theta = ( particles.position( 1, ipart ) + Icpx*particles.position( 2, ipart ) )/r;
        C_m = 2.;
        for( unsigned int i=0; i<( unsigned int )imode; i++ ) {
            C_m *= e_theta;
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
        Sr1[0] = Sr1[1];
        Sr1[1] = 0.; 
    }
 
    ip -= i_domain_begin ;
    jp -= j_domain_begin ;
    // j_domain_begin should always be zero in spectral since no paralellization along r.
    
    for( unsigned int i=0 ; i<2 ; i++ ) {
        iloc = ( i+ip )*nr+jp;
        for( unsigned int j=0 ; j<2 ; j++ ) {
            rhoj [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR[j+jp];
        }
    }//i
} // END Project for diags local current densities

// Apply boundary conditions on axis for Rho frozen particles
void ProjectorAM1Order::axisBCfrozen( complex<double> *rhoj,  int imode )
{

    return;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents for all modes, not charge conserving
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorAM1Order::currents( ElectroMagnAM *emAM, Particles &particles, unsigned int ipart, double invgf, bool diag_flag, int ispec)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    int iloc, jloc, linindex;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    
    // variable declaration
    double xpn, rpn;
    double delta;

    double  Sl1[2], Sr1[2],Sl1d[2], Sr1d[2];
    complex<double> e_theta, C_m = 1.; 
    complex<double> *Jl, *Jr, *Jt, *rho;
    
    double rp = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) );
    double crl_p =  ( particles.momentum( 0, ipart )) *invgf;
    double crt_p =  ( particles.momentum( 2, ipart )*particles.position( 1, ipart ) - particles.momentum( 1, ipart )*particles.position( 2, ipart ) )/rp * invgf;
    double crr_p =  ( particles.momentum( 1, ipart )*particles.position( 1, ipart ) + particles.momentum( 2, ipart )*particles.position( 2, ipart ) )/rp * invgf;
    e_theta = ( particles.position( 1, ipart ) + Icpx*particles.position( 2, ipart ) )/rp;
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dl_inv_ ;
    int ip = floor( xpn );
    delta  = xpn - ( double )ip;
    Sl1[0] = 1. - delta;
    Sl1[1] = delta;
    
    rpn = rp *dr_inv_ -0.5 ; //-0.5 because grid is shifted by dr/2
    int jp = floor( rpn );
    delta  = rpn - ( double )jp;
    Sr1[0] = 1. - delta;
    Sr1[1] = delta;

    if (rpn < 0.){ // If particle is between 0 and dr/2.
        jp = 0;
        Sr1[0] = Sr1[1];
        Sr1[1] = 0.; // Only account for deposition above axis. Symetry is handled in interpolation.
    }

    ip  -= i_domain_begin ;
    jp  -= j_domain_begin ;

    for( unsigned int imode=0; imode<( unsigned int )Nmode; imode++ ) {
        if( imode == 1 ) {
            C_m = 2.;
        }
        if( imode > 0 ) {
            C_m *= e_theta;
        }
        
        if (!diag_flag){
            Jl =  &( *emAM->Jl_[imode] )( 0 );
            Jr =  &( *emAM->Jr_[imode] )( 0 );
            Jt =  &( *emAM->Jt_[imode] )( 0 );
            rho = &( *emAM->rho_AM_[imode] )( 0 ) ; // In spectral, always project density
        } else {
            unsigned int n_species = emAM->Jl_s.size() / Nmode;
            unsigned int ifield = imode*n_species+ispec;
            Jl  = emAM->Jl_s    [ifield] ? &( * ( emAM->Jl_s    [ifield] ) )( 0 ) : &( *emAM->Jl_    [imode] )( 0 ) ;
            Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : &( *emAM->Jr_    [imode] )( 0 ) ;
            Jt  = emAM->Jt_s    [ifield] ? &( * ( emAM->Jt_s    [ifield] ) )( 0 ) : &( *emAM->Jt_    [imode] )( 0 ) ;
            rho = emAM->rho_AM_s[ifield] ? &( * ( emAM->rho_AM_s[ifield] ) )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;

        }

        // Jr^(p,p) Jt^(p,p) Jl^(p,p)
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )* nprimr + jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                linindex = iloc+j;
                complex<double> increment =  C_m*charge_weight* Sl1[i]*Sr1[j]*invR[jp+j];
                rho [linindex] += increment;
                Jl [linindex] += crl_p * increment ;
                Jr [linindex] += crr_p * increment ;
                Jt [linindex] += crt_p * increment ;
            }
        }//i
    }// end loop on modes
    
} // END Project local current densities (Jl, Jr, Jt)

//------------------------------------//
//Wrapper for projection
void ProjectorAM1Order::currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref )
{
        
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );

    for( int ipart=istart ; ipart<iend; ipart++ ) {
        currents( emAM, particles,  ipart, ( *invgf )[ipart], diag_flag, ispec);
    }

    //Boundary conditions for Jr mode 0 near axis
    if (emAM->isYmin ) {
        complex<double> *Jr; 
        if (!diag_flag){
            Jr =  &( *emAM->Jr_[0] )( 0 );
        } else {
            unsigned int n_species = emAM->Jl_s.size() / Nmode;
            unsigned int ifield = ispec;
            Jr  = emAM->Jr_s    [ifield] ? &( * ( emAM->Jr_s    [ifield] ) )( 0 ) : &( *emAM->Jr_    [0] )( 0 ) ;
        }

        // Jr_0[dr/2.] = Jr_0[3dr/2.] / 9.
        for( unsigned int i=0 ; i<npriml; i++ ) {
            int ilocr = i*nprimr+oversizeR;
            Jr [ilocr] = Jr [ilocr+1] * one_nineth;
        }//i
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

        int iloc, nr( nprimr );
    
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
    
        ip -= i_domain_begin ;
        jp -= j_domain_begin ;
    
    
        for( unsigned int i=0 ; i<2 ; i++ ) {
            iloc = ( i+ip )*nr+jp;
            for( unsigned int j=0 ; j<2 ; j++ ) {
                    Chi_envelope [iloc+j] += C_m*charge_weight* Sl1[i]*Sr1[j] * invR[j+jp];
            }
        }//i
    


    }




    
}

