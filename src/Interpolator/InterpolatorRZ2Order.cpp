#include "InterpolatorRZ2Order.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include "ElectroMagn.h"
#include "ElectroMagn3DRZ.h"
#include "cField2D.h"
#include "Particles.h"
#include <complex>
#include "dcomplex.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for InterpolatorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorRZ2Order::InterpolatorRZ2Order(Params &params, Patch *patch) : InterpolatorRZ(params, patch)
{

    dl_inv_ = 1.0/params.cell_length[0];
    dr_inv_ = 1.0/params.cell_length[1];
    nmodes = params.nmodes;

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc)
{

    //Treat mode 0 first

    // Static cast of the electromagnetic fields
    cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[0];
    cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[0];
    cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[0];
    cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_m[0];
    cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_m[0];
    cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_m[0];


    // Normalized particle position
    double xpn = particles.position(0, ipart) * dl_inv_;
    double r = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart)) ;
    double rpn = r * dr_inv_;
    //MESSAGE("rpn "<< rpn);
    exp_m_theta = ( particles.position(1, ipart) - Icpx * particles.position(2, ipart) ) / r ;   //exp(-i theta)
    
    complex<double> exp_mm_theta = 1. ;                                                                          //exp(-i m theta)
    //MESSAGE("expm "<< exp_m_theta);

    // Indexes of the central nodes
    ip_ = round(xpn);
    id_ = round(xpn+0.5);
    jp_ = round(rpn);
    jd_ = round(rpn+0.5);


    // Declaration and calculation of the coefficient for interpolation
    double delta2;

    deltax   = xpn - (double)id_ + 0.5;
    delta2  = deltax*deltax;
    coeffxd_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxd_[1] = 0.75 - delta2;
    coeffxd_[2] = 0.5 * (delta2+deltax+0.25);

    deltax   = xpn - (double)ip_;
    delta2  = deltax*deltax;
    coeffxp_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * (delta2+deltax+0.25);

    deltar   = rpn - (double)jd_ + 0.5;
    delta2  = deltar*deltar;
    coeffyd_[0] = 0.5 * (delta2-deltar+0.25);
    coeffyd_[1] = 0.75 - delta2;
    coeffyd_[2] = 0.5 * (delta2+deltar+0.25);

    deltar   = rpn - (double)jp_;
    delta2  = deltar*deltar;
    coeffyp_[0] = 0.5 * (delta2-deltar+0.25);
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * (delta2+deltar+0.25);

    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    id_ = id_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    jd_ = jd_ - j_domain_begin;

//Here we assume that mode 0 is real !!

    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    *(ELoc+0*nparts) = std::real (compute( &coeffxd_[1], &coeffyp_[1], ElRZ, id_, jp_));

    // -------------------------
    // Interpolation of Er^(p,d)
    // -------------------------
    *(ELoc+1*nparts) = std::real (compute( &coeffxp_[1], &coeffyd_[1], ErRZ, ip_, jd_));

    // -------------------------
    // Interpolation of Et^(p,p)
    // -------------------------
    *(ELoc+2*nparts) = std::real (compute( &coeffxp_[1], &coeffyp_[1], EtRZ, ip_, jp_));

    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    *(BLoc+0*nparts) = std::real (compute( &coeffxp_[1], &coeffyd_[1], BlRZ, ip_, jd_));

    // -------------------------
    // Interpolation of Br^(d,p)
    // -------------------------
    *(BLoc+1*nparts) = std::real (compute( &coeffxd_[1], &coeffyp_[1], BrRZ, id_, jp_));

    // -------------------------
    // Interpolation of Bt^(d,d)
    // -------------------------
    *(BLoc+2*nparts) = std::real (compute( &coeffxd_[1], &coeffyd_[1], BtRZ, id_, jd_));
    //MESSAGE("Elocx "<< *(ELoc+0*nparts));
    //MESSAGE("Elocy "<< *(ELoc+1*nparts));
    //MESSAGE("Elocz "<< *(ELoc+2*nparts));
    //MESSAGE("Blocx "<< *(BLoc+0*nparts));
    //MESSAGE("Blocy "<< *(BLoc+1*nparts));
    //MESSAGE("Blocz "<< *(BLoc+2*nparts));
 
    for (unsigned int imode = 1; imode < nmodes ; imode++){

        cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[imode];
        cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[imode];
        cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[imode];
        cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_m[imode];
        cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_m[imode];
        cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_m[imode];

        exp_mm_theta *= exp_m_theta ;
        
        *(ELoc+0*nparts) += std::real ( compute( &coeffxd_[1], &coeffyp_[1], ElRZ, id_, jp_) * exp_mm_theta ) ;
        *(ELoc+1*nparts) += std::real ( compute( &coeffxp_[1], &coeffyd_[1], ErRZ, ip_, jd_) * exp_mm_theta ) ;
        *(ELoc+2*nparts) += std::real ( compute( &coeffxp_[1], &coeffyp_[1], EtRZ, ip_, jp_) * exp_mm_theta ) ;
        *(BLoc+0*nparts) += std::real ( compute( &coeffxp_[1], &coeffyd_[1], BlRZ, ip_, jd_) * exp_mm_theta ) ;
        *(BLoc+1*nparts) += std::real ( compute( &coeffxd_[1], &coeffyp_[1], BrRZ, id_, jp_) * exp_mm_theta ) ;
        *(BLoc+2*nparts) += std::real ( compute( &coeffxd_[1], &coeffyd_[1], BtRZ, id_, jd_) * exp_mm_theta ) ;
        //MESSAGE("Elocx "<< *(ELoc+0*nparts));
        //MESSAGE("Elocy "<< *(ELoc+1*nparts));
        //MESSAGE("Elocz "<< *(ELoc+2*nparts));
        //MESSAGE("Blocx "<< *(BLoc+0*nparts));
        //MESSAGE("Blocy "<< *(BLoc+1*nparts));
        //MESSAGE("Blocz "<< *(BLoc+2*nparts));
        
    }

} // END InterpolatorRZ2Order

void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)
{

    // Interpolate E, B
    cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[0];
    cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[0];
    cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[0];
    cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_m[0];
    cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_m[0];
    cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_m[0];
    cField2D* JlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jl_[0]; 
    cField2D* JrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jr_[0]; 
    cField2D* JtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jt_[0]; 
    cField2D* RhoRZ= (static_cast<ElectroMagn3DRZ*>(EMfields))->rho_RZ_[0];
    
    
    // Normalized particle position
    double xpn = particles.position(0, ipart) * dl_inv_;
    double r = sqrt (particles.position(1, ipart)*particles.position(1, ipart)+particles.position(2, ipart)*particles.position(2, ipart)) ;
    double rpn = r * dr_inv_;
    exp_m_theta = ( particles.position(1, ipart) - Icpx * particles.position(2, ipart) ) / r ;
    complex<double> exp_mm_theta = 1. ;
    
    
    // Indexes of the central nodes
    ip_ = round(xpn);
    id_ = round(xpn+0.5);
    jp_ = round(rpn);
    jd_ = round(rpn+0.5);
    
    
    // Declaration and calculation of the coefficient for interpolation
    double delta2;
    
    deltax   = xpn - (double)id_ + 0.5;
    delta2  = deltax*deltax;
    coeffxd_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxd_[1] = 0.75 - delta2;
    coeffxd_[2] = 0.5 * (delta2+deltax+0.25);
    
    deltax   = xpn - (double)ip_;
    delta2  = deltax*deltax;
    coeffxp_[0] = 0.5 * (delta2-deltax+0.25);
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * (delta2+deltax+0.25);
    
    deltar   = rpn - (double)jd_ + 0.5;
    delta2  = deltar*deltar;
    coeffyd_[0] = 0.5 * (delta2-deltar+0.25);
    coeffyd_[1] = 0.75 - delta2;
    coeffyd_[2] = 0.5 * (delta2+deltar+0.25);
    
    deltar   = rpn - (double)jp_;
    delta2  = deltar*deltar;
    coeffyp_[0] = 0.5 * (delta2-deltar+0.25);
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * (delta2+deltar+0.25);
    
    
    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    id_ = id_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    jd_ = jd_ - j_domain_begin;
    
    
    // -------------------------
    // Interpolation of Ex^(d,p)
    // -------------------------
    (*ELoc).x = std::real(compute( &coeffxd_[1], &coeffyp_[1], ElRZ, id_, jp_));
    
    // -------------------------
    // Interpolation of Ey^(p,d)
    // -------------------------
    (*ELoc).y = std::real(compute( &coeffxp_[1], &coeffyd_[1], ErRZ, ip_, jd_));
    
    // -------------------------
    // Interpolation of Ez^(p,p)
    // -------------------------
    (*ELoc).z = std::real(compute( &coeffxp_[1], &coeffyp_[1], EtRZ, ip_, jp_));
    
    // -------------------------
    // Interpolation of Bx^(p,d)
    // -------------------------
    (*BLoc).x = std::real(compute( &coeffxp_[1], &coeffyd_[1], BlRZ, ip_, jd_));
    
    // -------------------------
    // Interpolation of By^(d,p)
    // -------------------------
    (*BLoc).y = std::real(compute( &coeffxd_[1], &coeffyp_[1], BrRZ, id_, jp_));
    
    // -------------------------
    // Interpolation of Bz^(d,d)
    // -------------------------
    (*BLoc).z = std::real(compute( &coeffxd_[1], &coeffyd_[1], BtRZ, id_, jd_));
    // -------------------------
    // Interpolation of Jx^(d,p,p)
    // -------------------------
    (*JLoc).x = std::real(compute( &coeffxd_[1], &coeffyp_[1], JlRZ, id_, jp_));
    
    // -------------------------
    // Interpolation of Jy^(p,d,p)
    // -------------------------
    (*JLoc).y = std::real(compute( &coeffxp_[1], &coeffyd_[1], JrRZ, ip_, jd_));
    
    // -------------------------
    // Interpolation of Jz^(p,p,d)
    // -------------------------
    (*JLoc).z = std::real(compute( &coeffxp_[1], &coeffyp_[1], JtRZ, ip_, jp_));
    
    // -------------------------
    // Interpolation of Rho^(p,p,p)
    // -------------------------
    (*RhoLoc) = std::real(compute( &coeffxp_[1], &coeffyp_[1], RhoRZ, ip_, jp_));

    for (unsigned int imode = 1; imode < nmodes ; imode++){

        cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[imode];
        cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[imode];
        cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[imode];
        cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_m[imode];
        cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_m[imode];
        cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_m[imode];
        cField2D* JlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jl_[imode]; 
        cField2D* JrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jr_[imode]; 
        cField2D* JtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Jt_[imode]; 
        cField2D* RhoRZ= (static_cast<ElectroMagn3DRZ*>(EMfields))->rho_RZ_[imode];

        exp_mm_theta *= exp_m_theta ;
        
        (*ELoc).x += std::real ( compute( &coeffxd_[1], &coeffyp_[1], ElRZ, id_, jp_) * exp_mm_theta ) ;
        (*ELoc).y += std::real ( compute( &coeffxp_[1], &coeffyd_[1], ErRZ, ip_, jd_) * exp_mm_theta ) ;
        (*ELoc).z += std::real ( compute( &coeffxp_[1], &coeffyp_[1], EtRZ, ip_, jp_) * exp_mm_theta ) ;
        (*BLoc).x += std::real ( compute( &coeffxp_[1], &coeffyd_[1], BlRZ, ip_, jd_) * exp_mm_theta ) ;
        (*BLoc).y += std::real ( compute( &coeffxd_[1], &coeffyp_[1], BrRZ, id_, jp_) * exp_mm_theta ) ;
        (*BLoc).z += std::real ( compute( &coeffxd_[1], &coeffyd_[1], BtRZ, id_, jd_) * exp_mm_theta ) ;
        (*JLoc).x += std::real ( compute( &coeffxd_[1], &coeffyp_[1], JlRZ, id_, jp_) * exp_mm_theta ) ;
        (*JLoc).y += std::real ( compute( &coeffxp_[1], &coeffyd_[1], JrRZ, ip_, jd_) * exp_mm_theta ) ;
        (*JLoc).z += std::real ( compute( &coeffxp_[1], &coeffyp_[1], JtRZ, ip_, jp_) * exp_mm_theta ) ;
        (*RhoLoc) += std::real ( compute( &coeffxp_[1], &coeffyp_[1], RhoRZ, ip_, jp_)* exp_mm_theta ) ;
        //MESSAGE("Elocx "<< (*ELoc).x);
        //MESSAGE("Elocy "<< (*ELoc).y);
        //MESSAGE("Elocz "<< (*ELoc).z);
        //MESSAGE("Blocx "<<(*BLoc).x);
        //MESSAGE("Blocy "<<(*BLoc).y);
        //MESSAGE("Blocz "<<(*BLoc).z);
    }

}

void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread)
{
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<std::complex<double>> *exp_m_theta_old = &(smpi->dynamics_thetaold[ithread]);

    //Loop on bin particles
    int nparts( particles.size() );
    for (int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        (*this)(EMfields, particles, ipart, nparts, &(*Epart)[ipart], &(*Bpart)[ipart]);
        //Buffering of iol and delta
        (*iold)[ipart+0*nparts]  = ip_;
        (*iold)[ipart+1*nparts]  = jp_;
        (*delta)[ipart+0*nparts] = deltax;
        (*delta)[ipart+1*nparts] = deltar;
        (*exp_m_theta_old)[ipart] = exp_m_theta;
        //MESSAGE("exp_m_theta "<< exp_m_theta);
    }

}


// Interpolator specific to tracked particles. A selection of particles may be provided
void InterpolatorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> * selection)
{
    ERROR("To Do");
}
