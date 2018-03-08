#include "ProjectorRZ2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn3DRZ.h"
#include "cField2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ProjectorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorRZ2Order::ProjectorRZ2Order (Params& params, Patch* patch) : ProjectorRZ(params, patch)
{
    dl_inv_   = 1.0/params.cell_length[0];
    dl_ov_dt  = params.cell_length[0] / params.timestep;
    dr_inv_   = 1.0/params.cell_length[1];
    dr_ov_dt  = params.cell_length[1] / params.timestep;
    
    one_third = 1.0/3.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ProjectorRZ2Order
// ---------------------------------------------------------------------------------------------------------------------
ProjectorRZ2Order::~ProjectorRZ2Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents (sort)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;

    
    
} // END Project local current densities (Jl, Jr, Jt, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (complex<double>* Jl, complex<double>* Jr, complex<double>* Jt, complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    // (x,y,z) components of the current density for the macro-particle
    
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);
    double crl_p = charge_weight*dl_ov_dt;
    double crr_p = charge_weight*dr_ov_dt;

   

} // END Project local densities (Jl, Jr, Jt, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{
    //Warning : this function is used for frozen species only. It is assumed that position = position_old !!!

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    
    int iloc,jloc;
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = (double)(particles.charge(ipart))*particles.weight(ipart);

   

} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void ProjectorRZ2Order::operator() (Field* Jl, Field* Jr, Field* Jt, Particles &particles, int ipart, LocalFields Jion)
{
    cField2D* Jl3D  = static_cast<cField2D*>(Jl);
    cField2D* Jr3D  = static_cast<cField2D*>(Jr);
    cField2D* Jt3D  = static_cast<cField2D*>(Jt);
    
    // weighted currents
    double Jx_ion = Jion.x * particles.weight(ipart);
    double Jy_ion = Jion.y * particles.weight(ipart);
    double Jz_ion = Jion.z * particles.weight(ipart);
    


} // END Project global current densities (ionize)

//Wrapper for projection
void ProjectorRZ2Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec)
{
    if (is_spectral)
        ERROR("Not implemented");

    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    
    int dim1 = EMfields->dimPrim[1];
    int dim2 = EMfields->dimPrim[2];

    ElectroMagn3DRZ* emRZ = static_cast<ElectroMagn3DRZ*>( EMfields );

    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if (!diag_flag){ 

        
        // Loop on modes ?
        int imode = 0;

        complex<double>* b_Jl =  &(*emRZ->Jl_[imode] )(ibin*clrw* dim1   * dim2   );
        complex<double>* b_Jr =  &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1)* dim2   );
        complex<double>* b_Jt =  &(*emRZ->Jt_[imode] )(ibin*clrw* dim1   *(dim2+1));
        for ( int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jl , b_Jr , b_Jt , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[3*ipart], &(*delta)[3*ipart]);
            
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        // Loop on modes ?
        int imode = 0;

        complex<double>* b_Jl  = emRZ->Jl_s [ispec] ? &(* static_cast<cField2D*>(emRZ->Jl_s [ispec]) )(ibin*clrw* dim1   *dim2) : &(*emRZ->Jl_[imode] )(ibin*clrw* dim1   *dim2) ;
        complex<double>* b_Jr  = emRZ->Jr_s [ispec] ? &(* static_cast<cField2D*>(emRZ->Jr_s [ispec]) )(ibin*clrw*(dim1+1)*dim2) : &(*emRZ->Jr_[imode] )(ibin*clrw*(dim1+1)*dim2) ;
        complex<double>* b_Jt  = emRZ->Jt_s [ispec] ? &(* static_cast<cField2D*>(emRZ->Jt_s [ispec]) )(ibin*clrw*dim1*(dim2+1)) : &(*emRZ->Jt_[imode] )(ibin*clrw*dim1*(dim2+1)) ;
        complex<double>* b_rho = emRZ->rho_s[ispec] ? &(* static_cast<cField2D*>(emRZ->rho_s[ispec]) )(ibin*clrw* dim1   *dim2) : &(*emRZ->rho_RZ_[imode])(ibin*clrw* dim1   *dim2) ;
        for ( int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jl , b_Jr , b_Jt ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[3*ipart], &(*delta)[3*ipart]);
    }

}
