// ----------------------------------------------------------------------------
//! \file Nlics.cpp
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles.
//
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
//
// ----------------------------------------------------------------------------

#include "Nlics.h"

#include <cstring>
#include <iostream>
#include <fstream>

#include <cmath>

#include "userFunctions.h"

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Nlics
// ---------------------------------------------------------------------------------------------------------------------
Nlics::Nlics(Params& params, Species * species)
{
    // Dimension position
    nDim_ = params.nDim_particle;

    // Time step
    dt    = params.timestep;

    // Inverse of the species mass
    one_over_mass_ = 1./species->mass;
}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Nlics
// ---------------------------------------------------------------------------------------------------------------------
Nlics::~Nlics()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Overloading of the operator (): perform the Discontinuous radiation reaction
//! induced by the nonlinear inverse Compton scattering
// ---------------------------------------------------------------------------------------------------------------------
void Nlics::operator() (Particles &particles,
        SmileiMPI* smpi,
        int istart,
        int iend,
        int ithread)
{
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    double charge_over_mass2;
    double upx, upy, upz, us2;
    double alpha, s, T2 ;
    double Tx, Ty, Tz;
    double pxsm, pysm, pzsm;

    // Temporary quantum parameter
    double chipa;

    // Temporary Lorentz factor
    double gamma;

    // Time to emission
    double emission_time;

    // time spent in the iteration
    double local_it_time;

    // Number of Monte-Carlo iteration
    int mc_it_nb;

    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );
    double* position[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position[i] =  &( particles.position(i,0) );
#ifdef  __DEBUG
    double* position_old[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position_old[i] =  &( particles.position_old(i,0) );
#endif

    // Charge
    short* charge = &( particles.charge(0) );

    // Optical depth for the Monte-Carlo process
    double* tau = &( particles.tau(0));

    // Optical depth for the Monte-Carlo process
    double* chi = &( particles.chi(0));

    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*pow(one_over_mass_,2.);

        // Init local variables
        emission_time = 0;
        local_it_time = 0;
        mc_it_nb = 0;

        // Monte-Carlo Manager inside the time step
        while ((local_it_time < dt)
             &&(mc_it_nb < mc_it_nb_max))
        {

            // Gamma
            gamma = 1./(*invgf)[ipart];

            // Computation of the Lorentz invariant quantum parameter
            chipa = Nlics::compute_chipa(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*Epart)[ipart].x,(*Epart)[ipart].y,(*Epart)[ipart].z,
                     (*Bpart)[ipart].x,(*Bpart)[ipart].y,(*Bpart)[ipart].z);


        }

    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the Lorentz invariant quantum parameter for particles
//
//! \param charge_over_mass2 charge divided by the square of the mass
//! \param px particle x momentum
//! \param gamma particle Lorentz factor
//! \param Ex x component of the particle electric field
//! \param Bx x component of the particle magnetic field
// ---------------------------------------------------------------------------------------------------------------------
double Nlics::compute_chipa(double & charge_over_mass2,
                             double & px, double & py, double & pz,
                             double & gamma,
                             double & Ex, double & Ey, double & Ez,
                             double & Bx, double & By, double & Bz)
{
    return 0;
}
