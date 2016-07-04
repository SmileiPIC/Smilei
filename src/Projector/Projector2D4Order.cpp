#include "Projector2D4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D4Order::Projector2D4Order (Params& params, Patch* patch) : Projector2D(params, patch)
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;

    one_third = 1.0/3.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector2D4Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D4Order::~Projector2D4Order()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double gf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaol)
{
    ERROR("Not defined");
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities & charge : diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double gf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaol)
{
    ERROR("Not defined");
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{
    ERROR("Not defined");
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization
// ---------------------------------------------------------------------------------------------------------------------
void  Projector2D4Order::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    ERROR("Projection of ionization current not yet defined for 2D 4nd order");
}


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D4Order::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, int diag_flag, std::vector<unsigned int> &b_dim, int ispec)
{
    std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    std::vector<double> *gf = &(smpi->dynamics_gf[ithread]);
    
    int dim1 = EMfields->dimPrim[1];
    
    if (diag_flag == 0){ 
        double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw*dim1);
        double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1));
        double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw*dim1);
        for (int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jx , b_Jy , b_Jz , particles,  ipart, (*gf)[ipart], ibin*clrw, b_dim, &(*iold)[2*ipart], &(*delta)[2*ipart]);
    } else {
        double* b_Jx =  &(*EMfields->Jx_s[ispec] )(ibin*clrw*dim1);
        double* b_Jy =  &(*EMfields->Jy_s[ispec] )(ibin*clrw*(dim1+1));
        double* b_Jz =  &(*EMfields->Jz_s[ispec] )(ibin*clrw*dim1);
        double* b_rho = &(*EMfields->rho_s[ispec])(ibin*clrw*dim1);
        for (int ipart=istart ; ipart<iend; ipart++ )
            (*this)(b_Jx , b_Jy , b_Jz ,b_rho, particles,  ipart, (*gf)[ipart], ibin*clrw, b_dim, &(*iold)[2*ipart], &(*delta)[2*ipart]);
    }
}
