// ----------------------------------------------------------------------------
//! \file RadiationMonteCarlo.h
//
//! \brief This class performs the Nonlinear Inverse Compton Scattering
//! on particles.
//
//! \details This header contains the definition of the class RadiationMonteCarlo.
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#ifndef RADIATIONMONTECARLO_H
#define RADIATIONMONTECARLO_H

#include "RadiationTables.h"
#include "Radiation.h"
#include "userFunctions.h"

//----------------------------------------------------------------------------------------------------------------------
//! RadiationMonteCarlo class: holds parameters and functions to apply the
//! nonlinear inverse Compton scattering on Particles.
//----------------------------------------------------------------------------------------------------------------------
class RadiationMonteCarlo : public Radiation
{

public:

    //! Constructor for RadiationMonteCarlo
    RadiationMonteCarlo( Params &params, Species *species, Random * rand  );
    
    //! Destructor for RadiationMonteCarlo
    ~RadiationMonteCarlo();
    
    // ---------------------------------------------------------------------
    //! Overloading of () operator: perform the Discontinuous radiation
    //! reaction induced by the nonlinear inverse Compton scattering
    //! \param particles   particle object containing the particles
    //! \param photon_species species that will receive emitted photons
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param RadiationTables Cross-section data tables and useful functions
    //                     for nonlinear inverse Compton scattering
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    //! \param radiated_energy     overall energy radiated during the call to this method
    // ---------------------------------------------------------------------
    virtual void operator()(
        Particles &particles,
        Species *photon_species,
        SmileiMPI *smpi,
        RadiationTables &RadiationTables,
        double          &radiated_energy,
        int             istart,
        int             iend,
        int             ithread,
        int             ipart_ref = 0
       );
        
    // ---------------------------------------------------------------------
    //! Perform the phoon emission (creation of a super-photon
    //! and slow down of the emitting particle)
    //! \param ipart              particle index
    //! \param particle_chi          particle quantum parameter
    //! \param particle_gamma          particle gamma factor
    //! \param position           particle position
    //! \param momentum           particle momentum
    //! \param RadiationTables    Cross-section data tables and useful functions
    //                        for nonlinear inverse Compton scattering
    // ---------------------------------------------------------------------
    double photonEmission( int ipart,
                         double &particle_chi,
                         double &particle_gamma,
                         double *position[3],
                         double *momentum[3],
                         double *weight,
                         Species *photon_species,
                         RadiationTables &RadiationTables );
                         
protected:

    // ________________________________________
    // General parameters
    
    //! Number of photons emitted per event for statisctics purposes
    int radiation_photon_sampling_;
    
    //! Threshold on the photon Lorentz factor under which the macro-photon
    //! is not generated but directly added to the energy scalar diags
    //! This enable to limit emission of useless low-energy photons
    double radiation_photon_gamma_threshold_;
    
    //! Inverse number of photons emitted per event for statisctics purposes
    double inv_radiation_photon_sampling_;
    
    //! Max number of Monte-Carlo iteration
    const int max_monte_carlo_iterations_ = 100;
    
    //! Espilon to check when tau is near 0
    const double epsilon_tau_ = 1e-100;
    
private:

};

#endif
