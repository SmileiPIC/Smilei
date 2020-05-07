// ----------------------------------------------------------------------------
//! \file RadiationDiagRadiationSpectrum.cpp
//
//! \brief This class is for the classical continuous radiation loss with
//!        the Landau-Lifshitz model.
//!        This model does not include a quantum correction.
//
//! References:
//! L. D. Landau and E. M. Lifshitz, The classical theory of fields, 1947
//! F. Niel et al., 2017
//
// ----------------------------------------------------------------------------

#include "RadiationDiagRadiationSpectrum.h"

// -----------------------------------------------------------------------------
//! Constructor for RadiationNLandauLifshitz
//! Inherited from Radiation
// -----------------------------------------------------------------------------
RadiationDiagRadiationSpectrum::RadiationDiagRadiationSpectrum(Params& params,
                                                 Species * species, Random * rand )
      : Radiation(params, species, rand)
{
}

// -----------------------------------------------------------------------------
//! Destructor for RadiationDiagRadiationSpectrum
// -----------------------------------------------------------------------------
RadiationDiagRadiationSpectrum::~RadiationDiagRadiationSpectrum()
{
}

// -----------------------------------------------------------------------------
//! Overloading of the operator (): perform the Landau-Lifshitz classical
//! radiation reaction
//
//! \param particles   particle object containing the particle properties
//! \param smpi        MPI properties
//! \param RadiationTables Cross-section data tables and useful functions
//                     for nonlinear inverse Compton scattering
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
//! \param radiated_energy     overall energy radiated during the call to this method
// -----------------------------------------------------------------------------
void RadiationDiagRadiationSpectrum::operator() (
      Particles &particles,
      Species *photon_species,
      SmileiMPI *smpi,
      RadiationTables &RadiationTables,
      double          &radiated_energy,
      int istart,
      int iend,
      int ithread, int ipart_ref)

{

    // _______________________________________________________________
    // Parameters
    std::vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    int nparts = particles.size();
    double* Ex = &( (*Epart)[0*nparts] );
    double* Ey = &( (*Epart)[1*nparts] );
    double* Ez = &( (*Epart)[2*nparts] );
    double* Bx = &( (*Bpart)[0*nparts] );
    double* By = &( (*Bpart)[1*nparts] );
    double* Bz = &( (*Bpart)[2*nparts] );

    // Charge divided by the square of the mass
    double charge_over_mass2;

    // 1/mass^2
    const double one_over_mass_2 = std::pow(one_over_mass_,2.);

    // Temporary Lorentz factor
    double gamma;

    // Momentum shortcut
    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    // Charge shortcut
    short* charge = &( particles.charge(0) );

    // Optical depth for the Monte-Carlo process
    double* chi = &( particles.chi(0));

    // _______________________________________________________________
    // Computation

    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass2 = (double)(charge[ipart])*one_over_mass_2;

        // Gamma
        gamma = sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                             + momentum[1][ipart]*momentum[1][ipart]
                             + momentum[2][ipart]*momentum[2][ipart]);

        // Computation of the Lorentz invariant quantum parameter
        chi[ipart] = Radiation::computeParticleChi(charge_over_mass2,
                     momentum[0][ipart],momentum[1][ipart],momentum[2][ipart],
                     gamma,
                     (*(Ex+ipart)),(*(Ey+ipart)),(*(Ez+ipart)),
                     (*(Bx+ipart)),(*(By+ipart)),(*(Bz+ipart)) );

    }
}
